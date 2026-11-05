"""polyfind's pre-registered predictions for E1, the intrinsic electrostriction references.

    python examples/e1_predictions.py [--out docs/predictions/e1_polyfind.json] [--systems beta,pe]

For beta-PVDF, PE, the alpha TGTG' cell and the 70/30 VDF/TrFE model cell, on the fitted
polarizable potential (``pvdf-dft-valence`` with induced dipoles and Ewald, rigid chain,
force-shifted Lennard-Jones, zero kelvin), in the definitions DFPT and the Berry phase give:
the relaxed-ion stiffness ``C_IJ`` (GPa; Voigt 1 = xx polar, 2 = yy lateral, 6 = xy; the chain
axis is rigid on this path, so there is no 3 column), the electronic and relaxed-ion
dielectric tensors at the stress-free cell, the proper piezoelectric ``e_iJ`` (C/m^2), the
true relaxed-ion dielectric tensor at +/-1% normal strain along the polar and lateral axes
with its central-difference derivative ``d eps_ii / d eta_J``, the true polarization and
stress at each strained cell (so the improper ``dP/d eta`` follows beside the proper ``e``),
and ``M`` and ``Q`` by both routes.  Frame: x polar, y lateral, z chain.  Dipole's periodic
reference uses x = a (lateral), y = b (polar), z = chain: swap x and y to compare.
"""
import argparse
import json
import time

import numpy as np

from polyfind import electrostriction as ES
from polyfind import mechanics as M
from polyfind.polymers import PE, PVDF, VDF_TRFE_7_3

T, GP, GM = 0, 1, 2
SYSTEMS = {
    "beta": ("beta-PVDF (TTTT)", PVDF, [T, T]),
    "pe": ("PE (control)", PE, [T]),
    "alpha": ("alpha-PVDF (TGTG')", PVDF, [T, GP, T, GM]),
    "vdftrfe": ("VDF/TrFE 70/30 model (TTTT)", VDF_TRFE_7_3, [T] * 20),
}


def build(polymer, seq, label):
    from polyfind.ewald import EwaldSpec
    from polyfind.polarizability import Polarizable
    kw = {"pack_kw": {"table_cache_dir": None}, "lj_cutoff": "force", "coulomb": "ewald", "ewald": EwaldSpec(),
          "polarizable": Polarizable()}
    ref, _ = M.refined_reference(polymer, seq, label=label, **kw)
    return M.relax_reference(ref)


true_dielectric = ES.true_dielectric


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="docs/predictions/e1_polyfind.json")
    ap.add_argument("--systems", default="beta,pe,alpha,vdftrfe")
    ap.add_argument("--strain", type=float, default=0.01)
    ap.add_argument("--skip-electrostriction", action="store_true", help="omit the M/Q routes (the costly part on a large cell)")
    ap.add_argument("--skip-converse", action="store_true", help="omit the converse piezoelectric route")
    ap.add_argument("--merge", action="store_true", help="update the systems in an existing output file instead of replacing it")
    args = ap.parse_args()
    from polyfind.fitting import FITTED_VALENCE
    results = {"potential": "pvdf-dft-valence, induced dipoles, Ewald, rigid chain, force-shifted LJ, 0 K",
               "frame": "x polar, y lateral, z chain; Voigt 1 xx, 2 yy, 6 xy; Dipole's frame has x and y swapped",
               "strain": args.strain, "systems": {}}
    if args.merge:
        try:
            with open(args.out) as fh:
                results = json.load(fh)
        except FileNotFoundError:
            pass
    with FITTED_VALENCE.applied():
        for key in args.systems.split(","):
            label, polymer, seq = SYSTEMS[key]
            t0 = time.time()
            ref = build(polymer, seq, label)
            el = M.elastic_constants(ref)
            di = M.dielectric_tensor(ref)
            pz = M.piezoelectric(ref, el, converse=not args.skip_converse)
            es = None if args.skip_electrostriction else ES.electrostriction(ref, el)
            zero = ES._state(ref, np.zeros(6), None, None, ES.TIGHT)
            out = {"label": label, "sequence": seq,
                   "cell": {"a": float(ref.params[0]), "b": float(ref.params[1]), "gamma_deg": float(ref.params[2]),
                            "c": float(ref.c), "volume_A3": float(ref.volume)},
                   "C_GPa": {"voigt": list(int(k) + 1 for k in el.reachable), "block": el.block.tolist()},
                   "eps_inf": di.clamped.tolist(), "eps_relaxed_ion": di.relaxed.tolist(),
                   "e_proper_C_m2": {"rows": "x y z", "cols_voigt": list(int(k) + 1 for k in el.reachable), "values": pz.e.tolist()},
                   "d_direct_pm_V": pz.d_direct.tolist(),
                   "P0_C_m2": zero.polarization.tolist(), "stress0_GPa": zero.stress.tolist(),
                   "strained": {}}
            if es is not None:
                out.update({"M_strain_m2_V2": es.M_strain.tolist(), "M_polarization_m2_V2": es.M_polarization.tolist(),
                            "Q_m4_C2": es.Q.tolist(), "chi_free": es.chi_free.tolist(), "chi_clamped": es.chi_clamped.tolist()})
            else:
                out["electrostriction"] = "not computed (--skip-electrostriction)"
            print(f"\n=== {label}: a {ref.params[0]:.3f} b {ref.params[1]:.3f} gamma {ref.params[2]:.2f} c {ref.c:.4f} A, "
                  f"V {ref.volume:.2f} A^3   (build and tensors {time.time() - t0:.0f} s)", flush=True)
            print(f"C (GPa, Voigt {out['C_GPa']['voigt']}):\n{np.array2string(el.block, precision=2, suppress_small=True)}")
            print(f"eps_inf diag {np.diag(di.clamped).round(4).tolist()}  eps_0 (relaxed-ion) diag {np.diag(di.relaxed).round(4).tolist()}")
            print(f"e proper (C/m^2):\n{np.array2string(pz.e, precision=4, suppress_small=True)}")
            print(f"P0 {np.round(zero.polarization, 4).tolist()} C/m^2")
            for J, name in ((0, "xx (polar)"), (1, "yy (lateral)")):
                if J not in el.reachable:
                    continue
                t1 = time.time()
                cells = {}
                for sgn in (1.0, -1.0):
                    e6 = np.zeros(6)
                    e6[J] = sgn * args.strain
                    eps = true_dielectric(ref, e6)
                    st = ES._state(ref, e6, None, None, ES.TIGHT)
                    cells[f"{sgn * args.strain:+.3f}"] = {"eps": eps.tolist(), "P_C_m2": st.polarization.tolist(),
                                                           "stress_GPa": st.stress.tolist(), "energy": float(st.energy)}
                plus, minus = cells[f"{args.strain:+.3f}"], cells[f"{-args.strain:+.3f}"]
                deps = (np.array(plus["eps"]) - np.array(minus["eps"])) / (2.0 * args.strain)
                dp = (np.array(plus["P_C_m2"]) - np.array(minus["P_C_m2"])) / (2.0 * args.strain)
                out["strained"][name] = {"cells": cells, "d_eps_d_eta": deps.tolist(), "dP_true_d_eta_C_m2": dp.tolist()}
                print(f"strain along {name}: d eps_ii/d eta = {np.diag(deps).round(3).tolist()}; "
                      f"dP_true/d eta = {dp.round(4).tolist()} C/m^2; eps at +1%: {np.diag(plus['eps']).round(4).tolist()}   ({time.time() - t1:.0f} s)", flush=True)
            if es is not None:
                print(es.table())
            results["systems"][key] = out
            with open(args.out, "w") as fh:
                json.dump(results, fh, indent=1)
            print(f"[{label} done, {time.time() - t0:.0f} s; written {args.out}]", flush=True)


if __name__ == "__main__":
    main()
