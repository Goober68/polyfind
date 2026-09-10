"""Elastic, piezoelectric and actuator response of the reference polymorphs.

Runs :mod:`polyfind.mechanics` on PVDF beta, alpha and gamma and on polyethylene as a null
control, under the illustrative potential and under a fitted preset, and prints the table
``docs/ELECTROMECHANICS.md`` records.

    python examples/electromechanics.py                        # both potentials, all four cases
    python examples/electromechanics.py --preset illustrative  # the built-in potential only
    python examples/electromechanics.py --axial                # add the axial-strain diagnosis

A few seconds per case: the packing and the refinement dominate, and the whole response
sweep costs about a hundred kernel rows on top of them.
"""
import argparse
import time

import numpy as np

from polyfind import mechanics as M
from polyfind.polymers import PE, PVDF

T, GP, GM = 0, 1, 2
CASES = [
    (PVDF, [T, T], "PVDF beta (TTTT)"),
    (PVDF, [T, GP, T, GM], "PVDF alpha (TGTG')"),
    (PVDF, [T, T, T, GP, T, T, T, GM], "PVDF gamma (T3GT3G')"),
    (PE, [T], "PE (all-trans, null control)"),
]


def run_case(polymer, seq, label, axial=False):
    t0 = time.time()
    ref, rr = M.refined_reference(polymer, seq, label=label)
    t_build = time.time() - t0
    resp = M.electromechanical_response(ref, polymer=polymer, axial=axial)
    return resp, rr, t_build


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--preset", default="both", choices=["both", "illustrative", "fitted"],
                    help="which potential(s) to run")
    ap.add_argument("--axial", action="store_true", help="also run the axial-strain diagnosis")
    args = ap.parse_args()

    runs = []
    if args.preset in ("both", "illustrative"):
        runs.append(("illustrative", None))
    if args.preset in ("both", "fitted"):
        from polyfind.fitting import FITTED_DFT

        runs.append(("pvdf-dft-fit", FITTED_DFT))

    for pot_name, params in runs:
        print("=" * 100)
        print(f"potential: {pot_name}")
        if params is not None:
            print("  " + params.describe())
        print("=" * 100)
        ctx = params.applied() if params is not None else None
        if ctx is not None:
            ctx.__enter__()
        try:
            for polymer, seq, label in CASES:
                resp, rr, t_build = run_case(polymer, seq, label, axial=args.axial)
                print()
                print(f"--- {label}   (pack+refine {t_build:.1f} s)")
                print(f"    refined cell a={resp.reference.params[0]:.3f} b={resp.reference.params[1]:.3f} "
                      f"gamma={resp.reference.params[2]:.2f} c={resp.reference.c:.4f}  "
                      f"residual sigma (GPa) = {np.round(resp.elastic.residual_stress, 4)}")
                print("    " + resp.summary().replace("\n", "\n    "))
                if resp.axial is not None:
                    ax = resp.axial
                    print(f"    axial: eps_zz reachable over {ax.strain_range[0] * 100:+.2f}% .. "
                          f"{ax.strain_range[1] * 100:+.2f}% via the backbone angles; "
                          f"sigma_zz(rigid chain) = {ax.sigma_zz_rigid:+.3f} GPa")
                    print("    axial: C_33 (GPa) against the bend constant k (kcal/mol/rad^2): "
                          + ", ".join(f"k={k:g}: {v:.1f}" for k, v in sorted(ax.c33.items())))
                    print("    axial: sigma_zz (GPa) along the same path:                     "
                          + ", ".join(f"k={k:g}: {v:+.2f}" for k, v in sorted(ax.sigma_zz.items())))
                    print(f"    axial verdict: computable={ax.computable}; the invented bend constant "
                          f"supplies {ax.bend_fraction * 100:.0f}% of C_33 at k=105")
        finally:
            if ctx is not None:
                ctx.__exit__(None, None, None)


if __name__ == "__main__":
    main()
