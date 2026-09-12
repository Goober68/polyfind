"""Born effective charges of beta-PVDF: measure the model's, fit the charge flux to the reference's.

    python examples/fit_born_flux.py born                  # our Born tensor vs periodic DFPT, nothing fitted
    python examples/fit_born_flux.py fit                   # fit the flux to the transverse components
    python examples/fit_born_flux.py response [--preset P] # d33, d31, C33, cell, |P| with a preset
    python examples/fit_born_flux.py ordering --preset P   # E(alpha) - E(beta) with a preset
    python examples/fit_born_flux.py pe --preset P         # the polyethylene null with a preset
    python examples/fit_born_flux.py all

Reproduces ``docs/ELECTROMECHANICS.md`` section 5.9.  The reference is the sibling project's
periodic PBE-D3 DFPT Born tensor (``--data``, default the sarco checkout next to this one),
whose transverse components (``aa``, ``bb``, ``ab``, ``ba``; a = long axis, b = polar) are the
usable ones: the chain-axis (``cc``) components carry an imposed acoustic-sum correction of
0.052 e per atom on a raw residual of -0.62 e and are **reported as a check, never fitted**.

**What is fitted, and to what.**  At fixed geometry the cell dipole is linear in every flux
coefficient (the charges are, and the induced dipoles are linear in the charges), so the Born
tensor is too and the fit is an ordinary linear least squares with exact sensitivity columns.
The piezoelectric coefficients never enter: they stay out of sample.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np

from polyfind import mechanics as M
from polyfind.born import (COMPONENTS, TRANSVERSE, BornReference, born_charges, canonical,
                           comparison_table, load_reference)
from polyfind.ewald import EwaldSpec
from polyfind.fitting import FITTED_VALENCE, KJ_PER_KCAL
from polyfind.forcefield import SimpleFF
from polyfind.lattice_table import clear_pair_table_cache
from polyfind.pack import chain_flux
from polyfind.polarizability import Polarizable
from polyfind.polymers import PE, PVDF

T, GP, GM = 0, 1, 2
BETA = (PVDF, [T, T], "PVDF beta (TTTT)")
ALPHA = (PVDF, [T, GP, T, GM], "PVDF alpha (TGTG')")
MEASURED = {"d33": -32.0, "d31": 20.0, "C33_dft": 315.9, "P_dft": (0.176, 0.188)}
EW = dict(coulomb="ewald", ewald=EwaldSpec())
DEFAULT_DATA = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "..", "sarco", "materials",
                                             "gpu_bundle", "periodic_reference", "beta_pvdf", "born_results.json"))
TYPES = ("C(F2)", "F", "C(H2)", "H")
KEYS = (("C", "H", "angle"), ("C", "H", "bond"), ("C", "F", "angle"), ("C", "F", "bond"))


def flux_ff(coeffs: dict) -> SimpleFF:
    """``pvdf-dft-valence`` with ``{(e_i, e_j): (k_angle, k_bond)}`` flux coefficients."""
    base = SimpleFF.from_preset("pvdf-dft-valence")
    kw = {k: v for k, v in base.__dict__.items() if k != "_cache"}
    kw["charge_flux"] = tuple((a, b, float(ka), float(kb)) for (a, b), (ka, kb) in coeffs.items())
    return SimpleFF(**kw)


def preset_coeffs(name: str) -> dict:
    ff = SimpleFF.from_preset(name)
    return dict(ff.flux_table())


def beta_reference(flx: SimpleFF | None, pol: Polarizable | None):
    """beta-PVDF packed, refined and relaxed over cell and shape: ``(packer, params, shape, rel)``."""
    val = SimpleFF.from_preset("pvdf-dft-valence")
    clear_pair_table_cache()
    with FITTED_VALENCE.applied():
        polymer, seq, label = BETA
        ref, rr = M.refined_reference(polymer, seq, label=label, valence=val, refine_kw={"valence": val},
                                      charge_flux=flx, pack_kw={"table_cache_dir": None},
                                      polarizable=pol, **EW)
        shape = M.shape_of(polymer, ref, angles=rr.angles)
        rel = M.relax_deformable(ref, shape, ref.params, shape.x0, free=M._CELL_FREE, c_target=None, maxiter=300)
        ref.packer.update_chain(rel.chain)
    return ref.packer, rel.params, shape, rel


def base_charges(packer) -> np.ndarray:
    """The zero-flux (bond-charge-increment, off-site projected) charges at the packer's geometry (n,)."""
    from polyfind.pack import _repeat_images

    zero = flux_ff({("C", "H"): (0.0, 0.0), ("C", "F"): (0.0, 0.0)})
    tpl, where = _repeat_images(packer.chain)
    top = zero.flux_topology(tpl.elements, tpl.bonds, images=where)
    return top.charges(packer.chain.coords, packer.chain.c)[0]


def set_flux(packer, flx: SimpleFF | None) -> None:
    """Swap the packer's flux topology (and hence its charges) at fixed geometry.

    A relaxed chain's own ``charges`` attribute is whatever the line group attached and is
    never read by a fluxing packer, so the zero-flux charges are re-derived from the
    increments rather than taken from it.
    """
    packer.charge_flux = flx
    packer._flux = chain_flux(packer.chain, flx)
    if packer._flux is None:
        packer._set_charge_tables(base_charges(packer))
    packer.update_chain(packer.chain)


def static_charges(packer, base: bool = False) -> dict:
    """Mean charge per type (e): the packer's own (fluxed, at its geometry) or the zero-flux base."""
    from polyfind.born import atom_labels
    labels = atom_labels(packer.chain)
    q = base_charges(packer) if base else np.asarray(packer._q_cell, dtype=float)[:packer.n]
    return {lab: float(np.mean([q[i] for i, l in enumerate(labels) if l == lab])) for lab in dict.fromkeys(labels)}


def born_of(packer, params, h: float = 1e-4):
    b = born_charges(packer, params, h=h)
    return b, canonical(b, cz=float(packer.chain.c))


def observations(ref: BornReference, corrected: bool = True, transpose: bool = False) -> np.ndarray:
    """The transverse reference components in the fit's order (see :func:`rows`)."""
    out = []
    for lab in TYPES:
        for nm in TRANSVERSE:
            if nm in ("ab", "ba") and lab.startswith("C"):
                continue
            n2 = {"ab": "ba", "ba": "ab"}.get(nm, nm) if transpose else nm
            out.append(ref.component(lab, n2, corrected))
    return np.array(out)


def rows() -> list:
    return [(lab, nm) for lab in TYPES for nm in TRANSVERSE if not (nm in ("ab", "ba") and lab.startswith("C"))]


def model_vector(can) -> np.ndarray:
    return np.array([can.component(lab, nm) for lab, nm in rows()])


def sensitivity(packer, params, base_coeffs: dict, keys=KEYS, unit: float = 1.0) -> tuple:
    """``(z0, A, cc0, Acc)``: the zero-flux Born vector and the exact unit-coefficient columns.

    Linearity in the coefficients is asserted at a random point rather than assumed.
    """
    zero = {(a, b): (0.0, 0.0) for a, b, _ in keys}
    set_flux(packer, None)
    _, can0 = born_of(packer, params)
    z0 = model_vector(can0)
    cc0 = np.array([can0.component(lab, "cc") for lab in TYPES])
    cols, ccs = [], []
    for a, b, which in keys:
        co = {k: list(v) for k, v in zero.items()}
        co[(a, b)][0 if which == "angle" else 1] = unit
        set_flux(packer, flux_ff({k: tuple(v) for k, v in co.items()}))
        _, can = born_of(packer, params)
        cols.append((model_vector(can) - z0) / unit)
        ccs.append((np.array([can.component(lab, "cc") for lab in TYPES]) - cc0) / unit)
    A, Acc = np.array(cols).T, np.array(ccs).T
    # linearity check at a random point
    rng = np.random.default_rng(7)
    k = rng.uniform(-1.0, 1.0, len(keys))
    co = {kk: [0.0, 0.0] for kk in zero}
    for (a, b, which), v in zip(keys, k):
        co[(a, b)][0 if which == "angle" else 1] = v
    set_flux(packer, flux_ff({kk: tuple(v) for kk, v in co.items()}))
    _, can = born_of(packer, params)
    lin = float(np.abs(model_vector(can) - (z0 + A @ k)).max())
    set_flux(packer, flux_ff(base_coeffs))
    return z0, A, cc0, Acc, lin


def lsq(A, r):
    k, *_ = np.linalg.lstsq(A, r, rcond=None)
    res = r - A @ k
    return k, float(np.sqrt(np.mean(res ** 2))), res


def coeffs_from(k, keys=KEYS) -> dict:
    co = {(a, b): [0.0, 0.0] for a, b, _ in keys}
    for (a, b, which), v in zip(keys, k):
        co[(a, b)][0 if which == "angle" else 1] = float(v)
    return {kk: tuple(v) for kk, v in co.items()}


def stage_born(args, packer=None, params=None) -> dict:
    ref = load_reference(args.data)
    if packer is None:
        packer, params, shape, rel = beta_reference(SimpleFF.from_preset(args.preset), Polarizable())
        print(f"# reference structure: {args.preset} + polarizable, a={params[0]:.3f} b={params[1]:.3f} "
              f"c={packer.chain.c:.4f}, shape {rel.x[0]:+.3f} deg, grad {rel.shape_gradient:.1e}")
    print(f"# DFPT reference: {ref.protocol}; raw acoustic sum diag {np.diag(ref.asr_raw)}; "
          f"eps_inf diag {np.diag(ref.eps_inf)}")
    out = {"preset": args.preset}
    t0 = time.time()
    b, can = born_of(packer, params)
    b2, can2 = born_of(packer, params, h=2e-4)
    fd = float(np.abs(b.Z - b2.Z).max())
    print(f"# our Born tensor, flux + induced dipoles: {time.time() - t0:.1f} s, |asr| max {np.abs(b.asr).max():.1e} e, "
          f"step-halving change {fd:.1e} e")
    qb, qf = static_charges(packer, base=True), static_charges(packer)
    print("# static charges (e): base increments " + ", ".join(f"{l} {qb[l]:+.3f}" for l in TYPES)
          + "; with the flux at this geometry " + ", ".join(f"{l} {qf[l]:+.3f}" for l in TYPES))
    print(comparison_table(can, ref, TYPES, static=qf))
    out["full"] = {lab: can.mean[lab].tolist() for lab in TYPES}
    # decomposition: which piece carries what
    pol, flx = packer.polarizable, packer.charge_flux
    pieces = {}
    packer.polarizable = None
    set_flux(packer, None)
    pieces["static"] = born_of(packer, params)[1]
    set_flux(packer, flx)
    pieces["static+flux"] = born_of(packer, params)[1]
    packer.polarizable = pol
    set_flux(packer, None)
    pieces["static+induced"] = born_of(packer, params)[1]
    set_flux(packer, flx)
    print("\n# decomposition, diagonal components (aa, bb, cc) in e:")
    print(f"{'type':8s} {'static':>22s} {'static+flux':>22s} {'static+induced':>22s} {'all':>22s} {'DFPT asr':>22s}")
    for lab in TYPES:
        cells = [pieces["static"], pieces["static+flux"], pieces["static+induced"], can]
        fmt = lambda m: "(" + ", ".join(f"{m[i, i]:+6.3f}" for i in range(3)) + ")"  # noqa: E731
        print(f"{lab:8s} " + " ".join(f"{fmt(c.mean[lab]):>22s}" for c in cells)
              + f" {fmt(ref.corrected[lab]):>22s}")
    out["pieces"] = {k: {lab: v.mean[lab].tolist() for lab in TYPES} for k, v in pieces.items()}
    # the group sums that a rigid pendant group's translation sees
    print("\n# group Born sums (the response to a rigid translation of the group), diagonal, e:")
    for grp, members in (("CF2", (("C(F2)", 1), ("F", 2))), ("CH2", (("C(H2)", 1), ("H", 2)))):
        ours = sum(w * can.mean[l] for l, w in members)
        theirs = sum(w * ref.corrected[l] for l, w in members)
        print(f"  {grp}: ours ({', '.join(f'{ours[i, i]:+.3f}' for i in range(3))})  "
              f"DFPT ({', '.join(f'{theirs[i, i]:+.3f}' for i in range(3))})")
    # raw vs corrected, transverse only
    dt = max(abs(ref.corrected[l][a, b] - ref.raw[l][a, b]) for l in TYPES for nm, a, b in COMPONENTS if nm in TRANSVERSE)
    dc = max(abs(ref.corrected[l][2, 2] - ref.raw[l][2, 2]) for l in TYPES)
    print(f"\n# reference raw vs corrected: transverse components move by at most {dt:.2e} e, chain-axis by {dc:.3f} e")
    out["raw_vs_corrected"] = {"transverse_max": dt, "chain_max": dc}
    return out


def stage_fit(args) -> dict:
    ref = load_reference(args.data)
    base = preset_coeffs(args.preset)
    packer, params, shape, rel = beta_reference(SimpleFF.from_preset(args.preset), Polarizable())
    print(f"# structure: {args.preset} + polarizable, a={params[0]:.3f} b={params[1]:.3f} c={packer.chain.c:.4f}")
    z0, A, cc0, Acc, lin = sensitivity(packer, params, base)
    print(f"# sensitivity columns for {[f'{a}-{b} {w}' for a, b, w in KEYS]}; linearity residual {lin:.1e} e")
    sv = np.linalg.svd(A, compute_uv=False)
    print(f"# design matrix {A.shape}: singular values {np.array2string(sv, precision=3)}")
    y = observations(ref)
    r = y - z0
    labels = [f"{l}:{n}" for l, n in rows()]
    print(f"# {len(y)} transverse observations ({len(TYPES)} types x aa, bb; F, H x ab, ba), 2 tied by the "
          f"acoustic sum rule the model satisfies identically -> {len(y) - 2} independent")
    out = {"structure": {"a": float(params[0]), "b": float(params[1]), "c": float(packer.chain.c)},
           "observations": dict(zip(labels, y.tolist())), "z0": dict(zip(labels, z0.tolist())),
           "columns": {f"{a}-{b} {w}": dict(zip(labels, A[:, j].tolist())) for j, (a, b, w) in enumerate(KEYS)},
           "fits": {}}
    rms0 = float(np.sqrt(np.mean(r ** 2)))
    print(f"# zero flux (static + induced dipoles): rms {rms0:.4f} e")
    kcur = np.array([base.get((a, b), (0.0, 0.0))[0 if w == "angle" else 1] for a, b, w in KEYS])
    rcur = r - A @ kcur
    print(f"# current preset {args.preset}: rms {np.sqrt(np.mean(rcur ** 2)):.4f} e")

    def report(name, idx, yv=y, tag="asr"):
        k, rms, res = lsq(A[:, idx], yv - z0)
        full = np.zeros(len(KEYS))
        full[idx] = k
        print(f"  {name:34s} {len(idx)} params for {len(yv) - 2} independent obs ({(len(yv) - 2) / len(idx):.1f}:1)  "
              f"rms {rms:.4f} e  k = " + ", ".join(f"{KEYS[j][0]}-{KEYS[j][1]} {KEYS[j][2][0]}={full[j]:+.4f}" for j in idx))
        out["fits"][f"{name}|{tag}"] = {"k": dict(zip([f"{a}-{b} {w}" for a, b, w in KEYS], full.tolist())),
                                        "rms": rms, "residuals": dict(zip(labels, res.tolist()))}
        return full, res

    print("\n# fits to the ASR-corrected transverse components:")
    k4, res4 = report("angle + bond, both bonds", [0, 1, 2, 3])
    k2b, _ = report("bond only, both bonds", [1, 3])
    k2f, _ = report("angle + bond, C-F only", [2, 3])
    k3, _ = report("C-H angle, C-F angle + bond", [0, 2, 3])
    print("  per-component residuals of the 4-parameter fit (model - DFPT, e):")
    for lab, v in zip(labels, res4):
        print(f"    {lab:10s} {-v:+.4f}")
    cc4 = cc0 + Acc @ k4
    print("  chain-axis components, NOT fitted, as a check (model | DFPT asr | DFPT raw):")
    for lab, v in zip(TYPES, cc4):
        print(f"    {lab:8s} cc {v:+.4f} | {ref.component(lab, 'cc'):+.4f} | {ref.component(lab, 'cc', False):+.4f}")
    print("\n# sensitivity of the 4-parameter fit:")
    report("raw (uncorrected) reference", [0, 1, 2, 3], observations(ref, corrected=False), tag="raw")
    report("off-diagonals transposed", [0, 1, 2, 3], observations(ref, transpose=True), tag="transposed")
    print("  leave one atom type out (fit on the other three, predict the held-out one):")
    for held in TYPES:
        keep = [i for i, (l, _) in enumerate(rows()) if l != held]
        drop = [i for i, (l, _) in enumerate(rows()) if l == held]
        k, rms_in, _ = lsq(A[keep], r[keep])
        pred = z0[drop] + A[drop] @ k
        rms_out = float(np.sqrt(np.mean((pred - y[drop]) ** 2)))
        print(f"    without {held:6s}: k = ({', '.join(f'{v:+.3f}' for v in k)})  in-sample rms {rms_in:.4f}  "
              f"held-out rms {rms_out:.4f} e  [" + ", ".join(f"{labels[i]} {pred[j]:+.3f} vs {y[i]:+.3f}" for j, i in enumerate(drop)) + "]")
        out["fits"][f"loo {held}"] = {"k": k.tolist(), "rms_in": rms_in, "rms_out": rms_out}
    out["fitted"] = {f"{a}-{b}": [float(v) for v in coeffs_from(k4)[(a, b)]] for a, b, _ in KEYS[::2]}
    return out


def beta_response(flx: SimpleFF, pol: Polarizable | None, args) -> dict:
    """d33, d31 (both routes), C33, cell and |P| for a flux preset, then its Born tensor at that structure."""
    t0 = time.time()
    packer, params, shape, rel = beta_reference(flx, pol)
    ref = M.Reference(packer=packer, params=params, c=float(packer.chain.c), label=BETA[2])
    with FITTED_VALENCE.applied():
        resp = M.electromechanical_response(ref, polymer=BETA[0], shape=shape)
    P, el, pz = resp.polarization, resp.elastic, resp.piezo
    sgn = 1.0 if P[0] >= 0 else -1.0
    out = {"flux": dict((f"{a}-{b}", v) for (a, b), v in flx.flux_table().items()) if flx is not None else None,
           "polarizable": None if pol is None else pol.label(),
           "shape_x": float(rel.x[0]), "shape_grad": float(rel.shape_gradient),
           "on_cap": bool(abs(abs(rel.x[0]) - shape.max_angle_change) < 1e-6),
           "a": float(resp.reference.params[0]), "b": float(resp.reference.params[1]),
           "gamma": float(resp.reference.params[2]), "c": float(resp.reference.c),
           "absP": float(np.linalg.norm(P)), "P": [float(v) for v in P],
           "C11": float(el.C[0, 0]), "C22": float(el.C[1, 1]), "C33": float(el.C[2, 2]),
           "C33_energy": float(el.c33_from_energy), "C13": float(el.C[0, 2]),
           "e_x": [float(v) for v in pz.e[0]],
           "d33_film": float(sgn * pz.d_improper[0, 0]), "d32_film": float(sgn * pz.d_improper[0, 1]),
           "d31_film": float(sgn * pz.d_improper[0, 2]),
           "d33_proper": float(sgn * pz.d_from_e[0, 0]), "d31_proper": float(sgn * pz.d_from_e[0, 2]),
           "d33_direct": float(sgn * pz.d_direct[0, 0]), "d31_direct": float(sgn * pz.d_direct[0, 2]),
           "d33_film_direct": float(sgn * (pz.d_direct[0, 0] + pz.d_improper[0, 0] - pz.d_from_e[0, 0])),
           "d31_film_direct": float(sgn * (pz.d_direct[0, 2] + pz.d_improper[0, 2] - pz.d_from_e[0, 2])),
           "route_rel_diff": float(pz.relative_difference), "res_zz": float(el.residual_stress[2]),
           "E_per_monomer": float(resp.reference.packer.energy(resp.reference.params[None])[0])
           / (2 * resp.reference.packer.chain.n_monomers)}
    # the Born tensor at the response's own (relaxed) structure
    if getattr(args, "data", None) and os.path.exists(args.data):
        refb = load_reference(args.data)
        b, can = born_of(resp.reference.packer, resp.reference.params)
        out["born"] = {lab: can.mean[lab].tolist() for lab in TYPES}
        y = observations(refb)
        out["born_rms_transverse"] = float(np.sqrt(np.mean((model_vector(can) - y) ** 2)))
    out["seconds"] = time.time() - t0
    return out


def response_row(r: dict) -> str:
    flag = "  NOT A MEASUREMENT" if (r["on_cap"] or r["route_rel_diff"] > 0.05) else ""
    born = "" if "born_rms_transverse" not in r else f"  Born rms {r['born_rms_transverse']:.4f} e"
    return (f"a={r['a']:.3f} b={r['b']:.3f} c={r['c']:.4f} |P|={r['absP']:.4f} C11={r['C11']:6.2f} C33={r['C33']:7.2f} "
            f"d33={r['d33_film']:+7.3f} (direct {r['d33_film_direct']:+7.3f}) d31={r['d31_film']:+6.3f} "
            f"(direct {r['d31_film_direct']:+6.3f}) proper d33={r['d33_proper']:+7.3f} d31={r['d31_proper']:+6.3f} "
            f"routes={r['route_rel_diff'] * 100:.2f}% shape {r['shape_x']:+.2f} deg{born}{flag}")


def ordering(flx: SimpleFF | None, pol: Polarizable | None) -> dict:
    """E(alpha) - E(beta), kJ/mol per monomer, rigid refined chains, cells relaxed, Ewald."""
    t0 = time.time()
    clear_pair_table_cache()
    out = {}
    with FITTED_VALENCE.applied():
        for key, (polymer, seq, label) in (("alpha", ALPHA), ("beta", BETA)):
            ref, rr = M.refined_reference(polymer, seq, label=label, pack_kw={"table_cache_dir": None},
                                          charge_flux=flx, polarizable=pol, **EW)
            ref = M.relax_reference(ref)
            pk = ref.packer
            e = float(pk.energy(ref.params[None])[0]) / (pk.n_chains * pk.chain.n_monomers)
            out[key] = {"E_per_monomer": e, "a": float(ref.params[0]), "b": float(ref.params[1]),
                        "c": float(ref.c), "absP": float(np.linalg.norm(ref.polarization()))}
    out["alpha_minus_beta_kj"] = (out["alpha"]["E_per_monomer"] - out["beta"]["E_per_monomer"]) * KJ_PER_KCAL
    out["seconds"] = time.time() - t0
    return out


def pe_null(flx: SimpleFF | None, pol: Polarizable | None) -> dict:
    t0 = time.time()
    clear_pair_table_cache()
    val = SimpleFF.from_preset("pvdf-dft-valence")
    with FITTED_VALENCE.applied():
        ref, rr = M.refined_reference(PE, [T], label="PE", valence=val, refine_kw={"valence": val},
                                      charge_flux=flx, pack_kw={"table_cache_dir": None}, polarizable=pol, **EW)
        shape = M.shape_of(PE, ref, angles=rr.angles)
        resp = M.electromechanical_response(ref, polymer=PE, shape=shape)
    return {"max_abs_e": float(np.abs(resp.piezo.e).max()), "max_abs_d": float(np.abs(resp.piezo.d_from_e).max()),
            "max_abs_d_direct": float(np.abs(resp.piezo.d_direct).max()), "absP": float(np.linalg.norm(resp.polarization)),
            "C33": float(resp.elastic.C[2, 2]), "a": float(resp.reference.params[0]),
            "b": float(resp.reference.params[1]), "c": float(resp.reference.c), "seconds": time.time() - t0}


def stage_pendant(args) -> dict:
    """Would letting the pendant geometry relax under axial strain change ``e_x,zz``?  A diagnostic.

    The deformable path holds every pendant rigid on its carbon, so under ``eps_zz`` the CF2 and
    CH2 groups translate as blocks and the polar dipole can only respond through the *group*
    Born sums, which the DFPT tensor puts at 0.08-0.15 e.  Here the F-C-F and H-C-H angles (and,
    with ``--pendant all``, the C-F and C-H lengths) are relaxed against the fitted valence terms
    at each strain, by rebuilding the line group from a polymer whose pendant geometry is the
    minimiser, and ``e_x,zz`` is differenced along that relaxed path.  Nothing here is shipped:
    it measures how much of the missing response the rigid pendant is hiding.
    """
    import dataclasses

    from scipy.optimize import minimize

    flx = SimpleFF.from_preset(args.preset)
    packer, params, shape, rel = beta_reference(flx, Polarizable())
    ref = M.Reference(packer=packer, params=params, c=float(packer.chain.c), label=BETA[2])
    with FITTED_VALENCE.applied():
        ref, shape = M.relax_reference_deformable(ref, shape)
        el = M.elastic_constants(ref, step=2e-3, shape=shape)
    tors, angs = shape.conformation(shape.x0)
    polymer = BETA[0]
    names = []
    for k in range(len(polymer.backbone)):
        if args.pendant in ("angles", "all"):
            names.append(("sub_angle", k))
        if args.pendant in ("bonds", "all"):
            names.append(("sub_bond", k))
    x0 = np.array([float(getattr(polymer.backbone[k], f)) for f, k in names])
    label = [f"{f}[{polymer.backbone[k].substituent}]" for f, k in names]
    step = 2e-3
    calls = [0]

    def polymer_at(x):
        bb = list(polymer.backbone)
        for (f, k), v in zip(names, x):
            bb[k] = dataclasses.replace(bb[k], **{f: float(v)})
        return dataclasses.replace(polymer, backbone=tuple(bb), sequence=())  # the monomer is rebuilt

    def state(x, eps_zz: float):
        calls[0] += 1
        eps6 = np.zeros(6)
        eps6[2] = eps_zz
        sh = M.shape_of(polymer_at(x), ref, torsions=tors, angles=angs, frame=shape.frame)
        with FITTED_VALENCE.applied():
            return M.deformable_state(ref, sh, eps6)

    def relaxed(eps_zz: float, start):
        res = minimize(lambda x: state(x, eps_zz).energy, start, method="Nelder-Mead",
                       options={"xatol": 2e-4, "fatol": 1e-9, "maxiter": 400, "initial_simplex": None})
        return np.asarray(res.x, dtype=float), state(res.x, eps_zz)

    t0 = time.time()
    rigid = {s: state(x0, s * step) for s in (-1.0, 0.0, 1.0)}
    e_rigid = (rigid[1.0].m - rigid[-1.0].m) / (2.0 * step)
    x_ref, st_ref = relaxed(0.0, x0)
    xp, st_p = relaxed(step, x_ref)
    xm, st_m = relaxed(-step, x_ref)
    e_relaxed = (st_p.m - st_m.m) / (2.0 * step)
    S = el.S
    conv = M.C_PER_M2_PER_GPA_TO_PC_PER_N
    # the axial column's contribution to the film d33 (poling axis x = Voigt 1): e_x,zz S_zz,xx
    d33_axial_rigid = float(e_rigid[0] * S[2, 0] * conv)
    d33_axial_relaxed = float(e_relaxed[0] * S[2, 0] * conv)
    out = {"preset": args.preset, "pendant": args.pendant, "coordinates": label, "x_polymer": x0.tolist(),
           "x_relaxed_reference": x_ref.tolist(), "x_relaxed_plus": xp.tolist(), "x_relaxed_minus": xm.tolist(),
           "energy_gain_kcal_per_cell": float(rigid[0.0].energy - st_ref.energy),
           "m_reference_rigid": rigid[0.0].m.tolist(), "m_reference_relaxed": st_ref.m.tolist(),
           "e_x_zz_rigid": float(e_rigid[0]), "e_x_zz_pendant_relaxed": float(e_relaxed[0]),
           "e_zz_rigid": e_rigid.tolist(), "e_zz_pendant_relaxed": e_relaxed.tolist(),
           "d33_axial_column_rigid": d33_axial_rigid, "d33_axial_column_pendant_relaxed": d33_axial_relaxed,
           "S_zz_xx_per_GPa": float(S[2, 0]), "state_evaluations": calls[0], "seconds": time.time() - t0}
    print(f"# pendant relaxation under eps_zz = +/-{step}, {args.preset} + polarizable, cell fixed at the reference")
    print(f"  coordinates {label}: polymer {np.array2string(x0, precision=4)} -> relaxed at eps=0 "
          f"{np.array2string(x_ref, precision=4)} (+eps {np.array2string(xp, precision=4)}, -eps {np.array2string(xm, precision=4)}); "
          f"energy gain {out['energy_gain_kcal_per_cell']:.4f} kcal/mol per cell")
    print(f"  |P| at eps=0: rigid {np.linalg.norm(rigid[0.0].m):.4f}, pendant-relaxed {np.linalg.norm(st_ref.m):.4f} C/m^2")
    print(f"  e_x,zz (C/m^2): rigid pendants {e_rigid[0]:+.4f}, pendants relaxed {e_relaxed[0]:+.4f}; "
          f"axial column of film d33 (pC/N): {d33_axial_rigid:+.3f} -> {d33_axial_relaxed:+.3f}  "
          f"[{calls[0]} strain states, {out['seconds']:.0f} s]")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["born", "fit", "response", "ordering", "pe", "pendant", "all"])
    ap.add_argument("--pendant", default="angles", choices=["angles", "bonds", "all"],
                    help="pendant: which pendant coordinates to relax")
    ap.add_argument("--data", default=DEFAULT_DATA, help="path to the provider's born_results.json")
    ap.add_argument("--preset", default="pvdf-dft-valence-flux", help="flux preset to measure (born/fit: the structure)")
    ap.add_argument("--compare", default=None, help="response: a second preset to run beside --preset")
    ap.add_argument("--json", default=None)
    args = ap.parse_args()
    rows_out = {}
    if args.mode in ("born", "all"):
        print(f"# === our Born tensor at the {args.preset} structure, before any fit ===")
        rows_out["born"] = stage_born(args)
    if args.mode in ("fit", "all"):
        print(f"\n# === fit of the flux to the transverse DFPT Born components, at the {args.preset} structure ===")
        rows_out["fit"] = stage_fit(args)
    if args.mode in ("response", "all"):
        print("\n# === d33, d31 (film convention, poling along +P), C33, cell, |P|; out of sample ===")
        names = [args.preset] + ([args.compare] if args.compare else [])
        rows_out["response"] = {}
        for name in names:
            r = beta_response(SimpleFF.from_preset(name), Polarizable(), args)
            rows_out["response"][name] = r
            print(f"{name:32s} " + response_row(r) + f"  {r['seconds']:.0f} s", flush=True)
            print(f"    vs measured d33 {MEASURED['d33']:+.0f}, d31 {MEASURED['d31']:+.0f}; C33 {r['C33']:.1f} "
                  f"(energy route {r['C33_energy']:.1f}) vs DFT {MEASURED['C33_dft']}; |P| {r['absP']:.4f} vs DFT "
                  f"{MEASURED['P_dft'][0]}-{MEASURED['P_dft'][1]}; E/mon {r['E_per_monomer']:.4f}", flush=True)
    if args.mode in ("ordering", "all"):
        r = ordering(SimpleFF.from_preset(args.preset), Polarizable())
        rows_out["ordering"] = r
        print(f"\n# E(alpha) - E(beta) with {args.preset} + polarizable: {r['alpha_minus_beta_kj']:+.3f} kJ/mol per monomer "
              f"(accepted -6.5 .. -2.6); beta a={r['beta']['a']:.3f} b={r['beta']['b']:.3f} |P|={r['beta']['absP']:.4f}; "
              f"alpha a={r['alpha']['a']:.3f} b={r['alpha']['b']:.3f} |P|={r['alpha']['absP']:.4f}  {r['seconds']:.0f} s")
    if args.mode == "pendant":
        rows_out["pendant"] = stage_pendant(args)
    if args.mode in ("pe", "all"):
        r = pe_null(SimpleFF.from_preset(args.preset), Polarizable())
        rows_out["pe"] = r
        print(f"\n# polyethylene null with {args.preset} + polarizable: max|e| {r['max_abs_e']:.1e} C/m^2, "
              f"max|d| {r['max_abs_d']:.1e} (direct {r['max_abs_d_direct']:.1e}) pC/N, |P| {r['absP']:.1e}, "
              f"C33 {r['C33']:.1f}, cell {r['a']:.3f} x {r['b']:.3f} x {r['c']:.4f}  {r['seconds']:.0f} s")
    if args.json:
        with open(args.json, "w") as f:
            json.dump(rows_out, f, indent=1)


if __name__ == "__main__":
    sys.exit(main())
