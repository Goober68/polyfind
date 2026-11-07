"""The relaxed-ion dielectric tensor of a packed cell from its Gamma-point lattice dynamics, beside
the rigid-chain one the prediction sets report.

    python examples/gamma_dielectric.py --preset pvdf-dft-valence --systems beta
    python examples/gamma_dielectric.py --preset pvdf-dft-valence-belt2 --variant belt --geometry pbe --systems beta

``mechanics.dielectric_tensor(relax=True)`` lets the chains reorient under the field at fixed cell,
the ionic part a rigid-chain parametrisation can express; DFPT's relaxed-ion tensor lets every
atom move.  This script builds the E1 reference exactly as ``examples/e1_predictions.py`` does
(the fit path's geometry), gives the packer the valence stretch and bend terms of
``pvdf-dft-valence`` pinned at that geometry (zero force at the reference, the valence fit's
stiffness for the stretches and the pendant bends) with the backbone angles carrying the
preset's own restraint about the built angle instead (the term the refinement balanced the
lateral pressure against, so the reference stays stationary), relaxes every atom at fixed cell,
and forms

    eps_0 = eps_inf + (4 pi k_e / V) Z^T H^+ Z

from the all-atom Hessian ``H`` (``polyfind.phonon.cell_energy_and_grad``, central differences),
the Born tensors ``Z`` (``polyfind.infrared.born_charges_by_field``) and the clamped-ion tensor
``eps_inf`` (``polyfind.infrared.electronic_dielectric``), with the three translations projected
out.  Belt sites (massless, chargeless) are not independent coordinates: each follows its
fluorine rigidly, exact for the rigid-chain modes that carry the lateral response and an
fluorines (it lies in the CF2 plane), which carries it exactly through any rigid motion of the
group.  The per-mode decomposition says which modes carry each axis.  Output: ``docs/predictions/gamma_dielectric_<preset>.json``.
"""
import argparse
import copy
import json
import time
from dataclasses import replace

import numpy as np
from scipy.linalg import eigh, null_space
from scipy.optimize import minimize

from polyfind import infrared as IR
from polyfind import mechanics as M
from polyfind import pendants as PD
from polyfind import phonon as PH
from polyfind.ewald import EwaldSpec
from polyfind.forcefield import SimpleFF
from polyfind.pack import COULOMB, MASS, chain_valence
from polyfind.polarizability import Polarizable
from polyfind.polymers import PE, PVDF, belt_variant

SYSTEMS = {"beta": ("beta-PVDF (TTTT)", PVDF, [0, 0]), "pe": ("PE (control)", PE, [0]),
           "alpha": ("alpha-PVDF (TGTG')", PVDF, [0, 1, 0, 2])}
KW = {"pack_kw": {"table_cache_dir": None}, "lj_cutoff": "force", "coulomb": "ewald", "ewald": EwaldSpec(),
      "polarizable": Polarizable()}


def build(polymer, seq, label):
    refine_kw = {k: KW[k] for k in ("lj_cutoff", "coulomb", "ewald", "polarizable")}
    ref, rr = M.refined_reference(polymer, seq, label=label, refine_kw=refine_kw, **KW)
    return M.relax_reference(ref), rr


K_BOND_MIN = None  # --k-bond-min: a floor on the pinned stretch constants (5.22)


def pinned_valence_packer(ref, polymer, k_backbone: float):
    """The reference's packer with valence terms stationary at its geometry; now :func:`polyfind.pendants.pinned_valence_packer`."""
    return PD.pinned_valence_packer(ref, polymer, k_backbone, k_bond_min=K_BOND_MIN)


def site_constraint(pk, Pn):
    """Belt sites as affine combinations of their CF2 atoms; now :func:`polyfind.pendants.belt_site_constraint`."""
    return PD.belt_site_constraint(pk, Pn)


def gamma_dielectric(ref, polymer, k_backbone: float, h: float = 1e-3) -> dict:
    """The relaxed-ion tensor of ``ref`` by the Gamma route (module docstring); everything the report needs."""
    params = np.asarray(ref.params, dtype=float)
    rigid = M.dielectric_tensor(ref)
    pk, n_backbone, theta_ref = pinned_valence_packer(ref, polymer, k_backbone)
    e_cell, e_pack, rel_ka = PH.check_known_answer(pk, params)
    P0, latn = PH.placed_coordinates(pk, params)
    P0 = np.asarray(P0, dtype=float)
    latn = np.asarray(latn, dtype=float)
    V = abs(float(np.linalg.det(latn)))
    J, atoms = site_constraint(pk, P0)
    Na = len(atoms)
    T = np.zeros((3 * Na, 3))
    for d in range(3):
        T[d::3, d] = 1.0
    B = null_space(T.T)  # (3Na, 3Na-3): the translation-free displacements
    m = B.shape[1]

    def eg(v):
        P = P0 + (J @ (B @ v)).reshape(-1, 3)
        E, g = PH.cell_energy_and_grad(pk, params, P, latn)
        return float(E), B.T @ (J.T @ np.asarray(g, dtype=float).ravel())

    def hessian(v):
        H = np.zeros((m, m))
        for k in range(m):
            e = np.zeros(m)
            e[k] = h
            H[:, k] = (eg(v + e)[1] - eg(v - e)[1]) / (2 * h)
        asym = float(np.linalg.norm(H - H.T) / max(np.linalg.norm(H), 1e-300))
        return 0.5 * (H + H.T), asym

    E0, g0 = eg(np.zeros(m))
    res = minimize(eg, np.zeros(m), jac=True, method="L-BFGS-B", options={"maxiter": 3000, "gtol": 1e-8, "ftol": 1e-14})
    v = res.x
    E1, g1 = eg(v)
    H, asym = hessian(v)
    for _ in range(4):  # Newton polish on the finite-difference Hessian
        if np.abs(g1).max() < 1e-7:
            break
        w, U = np.linalg.eigh(H)
        keep = w > 1e-6 * w.max()
        step = -(U[:, keep] @ ((U[:, keep].T @ g1) / w[keep]))
        if np.abs(step).max() > 0.05:
            step *= 0.05 / np.abs(step).max()
        E2, g2 = eg(v + step)
        if np.abs(g2).max() < np.abs(g1).max():
            v, E1, g1 = v + step, E2, g2
            H, asym = hessian(v)
        else:
            break
    Pn = P0 + (J @ (B @ v)).reshape(-1, 3)
    disp = float(np.abs(Pn - P0).max())
    # Born tensors and the clamped-ion tensor at the relaxed geometry
    Z = IR.born_charges_by_field(pk, params, Pn, latn)  # (N, 3, 3): Z[i, a, b] = d mu_a / d u_ib
    Zflat = np.zeros((3 * pk.N, 3))
    for a in range(3):
        Zflat[:, a] = Z[:, a, :].reshape(-1)
    Zb = B.T @ (J.T @ Zflat)  # (m, 3)
    eps_inf = IR.electronic_dielectric(pk, params, Pn, latn)
    # eps_0 = eps_inf + (4 pi k_e / V) Z^T H^+ Z, decomposed over the mass-weighted modes
    masses = np.repeat(np.array([MASS[e] for e in list(pk.elements) * pk.n_chains])[atoms], 3)
    Mb = B.T @ (masses[:, None] * B)
    lam, W = eigh(H, Mb)  # W is Mb-orthonormal
    pref = 4.0 * np.pi * COULOMB / V
    contrib = np.array([pref * np.outer(Zb.T @ W[:, i], Zb.T @ W[:, i]) / lam[i] if lam[i] > 1e-9 * lam.max() else np.zeros((3, 3))
                        for i in range(m)])
    eps_ion = contrib.sum(axis=0)
    eps0 = eps_inf + eps_ion
    freq = np.sign(lam) * np.sqrt(np.abs(lam)) * PH.SQRT_KCAL_A2_AMU_TO_CM1
    order = np.argsort(-np.array([np.trace(c) for c in contrib]))
    modes = [{"freq_cm1": float(freq[i]), "d_eps_diag": np.diag(contrib[i]).round(4).tolist()} for i in order[:8]]
    # the Born tensors themselves, per site (the inert sites carry zeros), for the comparison with DFPT's
    el = list(pk.elements) * pk.n_chains
    born = {"frame": "x polar, y lateral, z chain; Z[i][a][b] = d mu_a / d u_ib in e, at the relaxed geometry",
            "elements": el, "positions_A": Pn.tolist(), "Z": Z.tolist(),
            "sum_rule_max_abs": float(np.abs(Z.sum(axis=0)).max()),
            "by_element_mean_diag": {e: np.mean([np.diag(Z[i]) for i in range(pk.N) if el[i] == e], axis=0).round(4).tolist()
                                     for e in sorted(set(el)) if e != "X"}}
    return {"cell": {"a": float(params[0]), "b": float(params[1]), "gamma_deg": float(params[2]), "c": float(ref.c), "volume_A3": V},
            "born": born,
            "k_backbone": float(k_backbone), "refined_backbone_angle_deg": theta_ref, "n_backbone_angles": n_backbone,
            "n_bonds": pk._valence.n_bonds, "n_angles": pk._valence.n_angles, "n_free": m, "n_atoms": Na, "n_sites": pk.N,
            "relaxation": {"max_force_before": float(np.abs(g0).max()), "max_force_after": float(np.abs(g1).max()),
                           "max_displacement_A": disp, "energy_before": E0, "energy_after": E1, "hessian_asymmetry": asym,
                           "known_answer_rel": rel_ka},
            "frequencies_cm1": np.sort(freq).tolist(), "n_negative_modes": int((freq < -1.0).sum()),
            "eps_inf": eps_inf.tolist(), "eps_0_rigid_chain": rigid.relaxed.tolist(), "eps_0_gamma": eps0.tolist(),
            "eps_ion_gamma": eps_ion.tolist(), "modes_by_contribution": modes}


def report(label: str, out: dict, seconds: float) -> str:
    c, r = out["cell"], out["relaxation"]
    lines = [f"\n=== {label}: a {c['a']:.4f} b {c['b']:.4f} gamma {c['gamma_deg']:.2f} c {c['c']:.4f} A, V {c['volume_A3']:.2f} A^3  ({seconds:.0f} s)",
             f"pinned valence terms: {out['n_bonds']} bonds, {out['n_angles']} angles, {out['n_backbone_angles']} backbone angles at "
             f"k = {out['k_backbone']:.1f} kcal/mol/rad^2 about the built angle (refined {out['refined_backbone_angle_deg']:.2f} deg); "
             f"known answer rel {r['known_answer_rel']:.1e}",
             f"all-atom relaxation: max |force| {r['max_force_before']:.2e} -> {r['max_force_after']:.2e} kcal/(mol A); energy "
             f"{r['energy_before']:.6f} -> {r['energy_after']:.6f}; max displacement {r['max_displacement_A']:.4f} A; {out['n_free']} free "
             f"coordinates ({out['n_atoms']} atoms of {out['n_sites']} sites); Hessian asymmetry {r['hessian_asymmetry']:.1e}",
             f"lowest frequencies (cm^-1): {' '.join(f'{f:.1f}' for f in out['frequencies_cm1'][:8])}; negative: {out['n_negative_modes']}",
             f"eps_inf                 diag {np.diag(out['eps_inf']).round(4).tolist()}",
             f"eps_0 rigid-chain route diag {np.diag(out['eps_0_rigid_chain']).round(4).tolist()}   (mechanics.dielectric_tensor)",
             f"eps_0 Gamma-mode route  diag {np.diag(out['eps_0_gamma']).round(4).tolist()}   ionic part {np.diag(out['eps_ion_gamma']).round(4).tolist()}",
             "modes carrying the ionic part (frequency, d eps diag):"]
    lines += [f"  {md['freq_cm1']:8.1f} cm^-1  {md['d_eps_diag']}" for md in out["modes_by_contribution"]]
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--preset", default="pvdf-dft-valence")
    ap.add_argument("--variant", default="plain", choices=("plain", "belt"))
    ap.add_argument("--geometry", default="textbook", choices=("textbook", "pbe"))
    ap.add_argument("--systems", default="beta")
    ap.add_argument("--h", type=float, default=1e-3, help="finite-difference step of the Hessian (A)")
    ap.add_argument("--k-bond-min", type=float, default=None, help="floor on the pinned stretch constants, kcal/mol/A^2 in the kernel's 1/2 k convention (the valence fit's C-C constant sits at its bound, 100; a real C-C bond is about 540): a probe of what that softness does to the relaxed-ion tensor (5.22)")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    global K_BOND_MIN
    K_BOND_MIN = args.k_bond_min
    from polyfind.fitting import PRESET_PARAMETERS
    from polyfind.refine import resolve_angle_stiffness
    params_ff = PRESET_PARAMETERS[args.preset]
    out_path = args.out or f"docs/predictions/gamma_dielectric_{args.preset}.json"
    results = {"preset": args.preset, "variant": args.variant, "geometry": args.geometry,
               "frame": "x polar, y lateral, z chain", "systems": {}}
    with params_ff.applied():
        k_backbone = resolve_angle_stiffness(None, None)
        for key in args.systems.split(","):
            label, polymer, seq = SYSTEMS[key]
            if args.variant == "belt":
                polymer = belt_variant(polymer, geometry=args.geometry)
            t0 = time.time()
            ref, _ = build(polymer, seq, label)
            out = gamma_dielectric(ref, polymer, k_backbone, h=args.h)
            out["label"] = label
            print(report(f"{label} [{args.preset}, {args.variant}{', ' + args.geometry if args.variant == 'belt' else ''}]", out, time.time() - t0), flush=True)
            results["systems"][key] = out
            with open(out_path, "w") as fh:
                json.dump(results, fh, indent=1)
    print(f"\nwritten {out_path}")


if __name__ == "__main__":
    main()
