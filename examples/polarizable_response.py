"""Induced dipoles against the dielectric constant, then against d_33 and d_31.

Reproduces ``docs/ELECTROMECHANICS.md`` section 5.8.  The electrostatic-scale scan of
section 5.7 showed the piezoelectric coefficients wanting a response the static charges
cannot supply without destroying the energetics; this adds the response that does not live
in the static charges -- atomic polarizability (:mod:`polyfind.polarizability`) -- and measures
it in the order the brief demanded: the dielectric tensor first, against periodic DFPT, with
literature polarizabilities and nothing fitted; then the piezoelectric coefficients, which
stay out of sample.

Usage (``PYTHONPATH=src``)::

    python examples/polarizable_response.py dielectric    # eps_inf and eps_0 of beta, both potentials
    python examples/polarizable_response.py beta          # d_33, d_31, C_33, cell, |P| with and without
    python examples/polarizable_response.py ordering      # E(alpha) - E(beta) with and without
    python examples/polarizable_response.py pe            # the polyethylene null
    python examples/polarizable_response.py all
    python examples/polarizable_response.py beta --scale 0.5 1.0 1.5   # sensitivity, not a fit

Every comparison is Ewald against Ewald: induced dipoles need the lattice sum, so the
fixed-charge baseline here is the Ewald one (``docs/BENCHMARK.md`` "What Ewald costs"), and
the truncated-kernel numbers of section 5.6 are quoted beside it for reference only.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from dataclasses import replace

import numpy as np

from polyfind import mechanics as M
from polyfind.ewald import EwaldSpec
from polyfind.fitting import FITTED_VALENCE, KJ_PER_KCAL
from polyfind.forcefield import SimpleFF
from polyfind.lattice_table import clear_pair_table_cache
from polyfind.polarizability import Polarizable
from polyfind.polymers import PE, PVDF

T, GP, GM = 0, 1, 2
BETA = (PVDF, [T, T], "PVDF beta (TTTT)")
ALPHA = (PVDF, [T, GP, T, GM], "PVDF alpha (TGTG')")
# Periodic PBE-D3 DFPT on beta-PVDF (clamped-ion electronic tensor), axes x = a (8.58 A),
# y = b (polar), z = c (chain).  The packer's frame has x = polar, y = the long axis, so the
# first two entries swap when compared.
DFPT_EPS_INF_THEIR_AXES = (2.2526, 2.2354, 2.6005)
DFPT_EPS_INF_PACKER = (2.2354, 2.2526, 2.6005)
MEASURED = {"d33": -32.0, "d31": 20.0, "eps0_low_frequency": 10.0, "C33_dft": 315.9}
EW = dict(coulomb="ewald", ewald=EwaldSpec())


def beta_point(pol: Polarizable | None, cap: float | None = None, dielectric: bool = True) -> dict:
    """beta-PVDF on ``pvdf-dft-valence`` + ``-flux`` under Ewald, with or without induced dipoles."""
    t0 = time.time()
    val = SimpleFF.from_preset("pvdf-dft-valence")
    flx = SimpleFF.from_preset("pvdf-dft-valence-flux")
    clear_pair_table_cache()
    with FITTED_VALENCE.applied():
        polymer, seq, label = BETA
        ref, rr = M.refined_reference(polymer, seq, label=label, valence=val, refine_kw={"valence": val},
                                      charge_flux=flx, pack_kw={"table_cache_dir": None},
                                      polarizable=pol, **EW)
        shape = M.shape_of(polymer, ref, angles=rr.angles)
        if cap is not None:
            shape = replace(shape, max_angle_change=cap)
        rel = M.relax_deformable(ref, shape, ref.params, shape.x0, free=M._CELL_FREE, c_target=None, maxiter=300)
        resp = M.electromechanical_response(ref, polymer=polymer, shape=shape)
        P, el, pz = resp.polarization, resp.elastic, resp.piezo
        sgn = 1.0 if P[0] >= 0 else -1.0
        out = {"polarizable": None if pol is None else pol.label(), "cap": shape.max_angle_change,
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
               "route_rel_diff": float(pz.relative_difference), "res_zz": float(el.residual_stress[2]),
               "E_per_monomer": float(resp.reference.packer.energy(resp.reference.params[None])[0])
               / (2 * resp.reference.packer.chain.n_monomers)}
        if dielectric:
            de = M.dielectric_tensor(resp.reference, relax=True, shape=shape)
            out["eps_inf"] = [float(de.clamped[i, i]) for i in range(3)]
            out["eps_inf_offdiag"] = float(np.abs(de.clamped - np.diag(np.diag(de.clamped))).max())
            out["eps_0_fixed_cell"] = [float(de.relaxed[i, i]) for i in range(3)]
        out["seconds"] = time.time() - t0
        return out


def beta_row(r: dict) -> str:
    flag = "  NOT A MEASUREMENT" if (r["on_cap"] or r["route_rel_diff"] > 0.05) else ""
    eps = "" if "eps_inf" not in r else "  eps_inf=(" + ", ".join(f"{v:.3f}" for v in r["eps_inf"]) + ")"
    return (f"{(r['polarizable'] or 'fixed charges'):44s} a={r['a']:.3f} b={r['b']:.3f} c={r['c']:.4f} "
            f"|P|={r['absP']:.4f} C11={r['C11']:6.2f} C33={r['C33']:7.2f} "
            f"d33={r['d33_film']:+7.3f} d31={r['d31_film']:+6.3f} routes={r['route_rel_diff'] * 100:.2f}%"
            f"{eps}{flag}")


def ordering(pol: Polarizable | None) -> dict:
    """E(alpha) - E(beta) in kJ/mol per monomer: the first acceptance test, rigid chains, Ewald."""
    t0 = time.time()
    clear_pair_table_cache()
    out = {"polarizable": None if pol is None else pol.label()}
    with FITTED_VALENCE.applied():
        for key, (polymer, seq, label) in (("alpha", ALPHA), ("beta", BETA)):
            ref, rr = M.refined_reference(polymer, seq, label=label, pack_kw={"table_cache_dir": None},
                                          polarizable=pol, **EW)
            ref = M.relax_reference(ref)
            pk = ref.packer
            e = float(pk.energy(ref.params[None])[0]) / (pk.n_chains * pk.chain.n_monomers)
            out[key] = {"E_per_monomer": e, "a": float(ref.params[0]), "b": float(ref.params[1]),
                        "c": float(ref.c), "absP": float(np.linalg.norm(ref.polarization()))}
    out["alpha_minus_beta_kj"] = (out["alpha"]["E_per_monomer"] - out["beta"]["E_per_monomer"]) * KJ_PER_KCAL
    out["seconds"] = time.time() - t0
    return out


def pe_null(pol: Polarizable | None) -> dict:
    t0 = time.time()
    clear_pair_table_cache()
    with FITTED_VALENCE.applied():
        ref, _ = M.refined_reference(PE, [T], label="PE", pack_kw={"table_cache_dir": None}, polarizable=pol, **EW)
        resp = M.electromechanical_response(ref)
        de = M.dielectric_tensor(resp.reference, relax=False)
    return {"polarizable": None if pol is None else pol.label(),
            "max_abs_e": float(np.abs(resp.piezo.e).max()), "max_abs_d": float(np.abs(resp.piezo.d_from_e).max()),
            "max_abs_d_direct": float(np.abs(resp.piezo.d_direct).max()),
            "eps_inf": [float(de.clamped[i, i]) for i in range(3)], "C11": float(resp.elastic.C[0, 0]),
            "a": float(resp.reference.params[0]), "b": float(resp.reference.params[1]), "seconds": time.time() - t0}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["dielectric", "beta", "ordering", "pe", "all"])
    ap.add_argument("--scale", type=float, nargs="*", default=[1.0],
                    help="multiply every polarizability (sensitivity only; the shipped value is 1)")
    ap.add_argument("--thole", default="exp", choices=["exp", "cubic"])
    ap.add_argument("--thole-a", type=float, default=None)
    ap.add_argument("--cap", type=float, default=None)
    ap.add_argument("--json", default=None)
    args = ap.parse_args()

    def pol(scale):
        kw = dict(scale=scale, thole_form=args.thole)
        if args.thole_a is not None:
            kw["thole_a"] = args.thole_a
        elif args.thole == "cubic":
            kw["thole_a"] = 0.39
        return Polarizable(**kw)

    rows = []
    if args.mode in ("dielectric", "beta", "all"):
        print("# beta-PVDF, pvdf-dft-valence(-flux), Ewald tinfoil; DFPT eps_inf (packer axes x=polar, y, z=chain) = "
              + ", ".join(f"{v:.4f}" for v in DFPT_EPS_INF_PACKER))
        for p in [None] + [pol(s) for s in args.scale]:
            r = beta_point(p, cap=args.cap)
            rows.append(r)
            print(beta_row(r), flush=True)
            if "eps_inf" in r:
                print(f"    eps_inf = diag({', '.join(f'{v:.4f}' for v in r['eps_inf'])}), offdiag {r['eps_inf_offdiag']:.1e};  "
                      f"eps_0 (relaxed-ion, fixed cell) = diag({', '.join(f'{v:.4f}' for v in r['eps_0_fixed_cell'])});  "
                      f"measured low-frequency ~{MEASURED['eps0_low_frequency']:g}", flush=True)
            print(f"    d33 proper {r['d33_proper']:+.3f} (direct {r['d33_direct']:+.3f}), d31 proper {r['d31_proper']:+.3f} "
                  f"(direct {r['d31_direct']:+.3f}); film d33 {r['d33_film']:+.3f}, d31 {r['d31_film']:+.3f} vs measured "
                  f"{MEASURED['d33']:+.0f}, {MEASURED['d31']:+.0f};  C33 {r['C33']:.1f} (energy route {r['C33_energy']:.1f}) "
                  f"vs DFT {MEASURED['C33_dft']};  shape {r['shape_x']:+.3f} deg, grad {r['shape_grad']:.1e}; "
                  f"E/mon {r['E_per_monomer']:.4f}; {r['seconds']:.0f} s", flush=True)
    if args.mode in ("ordering", "all"):
        print("# E(alpha) - E(beta), kJ/mol per monomer, rigid refined chains, cells relaxed, Ewald (accepted -6.5 .. -2.6)")
        for p in [None] + [pol(s) for s in args.scale]:
            r = ordering(p)
            rows.append(r)
            print(f"{(r['polarizable'] or 'fixed charges'):44s} gap={r['alpha_minus_beta_kj']:+7.3f}  "
                  f"beta a={r['beta']['a']:.3f} b={r['beta']['b']:.3f} |P|={r['beta']['absP']:.4f}  "
                  f"alpha a={r['alpha']['a']:.3f} b={r['alpha']['b']:.3f} |P|={r['alpha']['absP']:.4f}  {r['seconds']:.0f} s", flush=True)
    if args.mode in ("pe", "all"):
        print("# polyethylene null: every piezoelectric coefficient must stay at machine zero")
        for p in [None] + [pol(s) for s in args.scale]:
            r = pe_null(p)
            rows.append(r)
            print(f"{(r['polarizable'] or 'fixed charges'):44s} max|e|={r['max_abs_e']:.1e} C/m^2  max|d|={r['max_abs_d']:.1e} "
                  f"(direct {r['max_abs_d_direct']:.1e}) pC/N  eps_inf=({', '.join(f'{v:.3f}' for v in r['eps_inf'])})  "
                  f"a={r['a']:.3f} b={r['b']:.3f} C11={r['C11']:.2f}  {r['seconds']:.0f} s", flush=True)
    if args.json:
        with open(args.json, "w") as f:
            json.dump(rows, f, indent=1)


if __name__ == "__main__":
    sys.exit(main())
