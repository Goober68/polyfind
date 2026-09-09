"""Reproduce the validation numbers quoted in DESIGN.md.

Part A packs the *experimental* chain conformations of PE and of alpha-, beta-
and gamma-PVDF at ideal torsion angles and compares the predicted cells with
the literature cells.  Part B runs the full funnel for PVDF with the built-in
illustrative potential.  Run with ``python examples/pvdf_polymorphs.py``
(a few minutes on 4 CPU cores; seconds with a GPU backend).
"""
import time

import numpy as np

from polyfind.pack import pack, periodic_chain
from polyfind.pipeline import EXPERIMENTAL_CELLS, PipelineConfig, run_pipeline
from polyfind.polymers import PE, PVDF, THREE_STATE
from polyfind.refine import refine_crystal

T, GP, GM = 0, 1, 2
CASES = [
    (PE, [T], "orthorhombic (all-trans)"),
    (PVDF, [T, T], "beta (TTTT)"),
    (PVDF, [T, GP, T, GM], "alpha/delta (TGTG')"),
    (PVDF, [T, T, T, GP, T, T, T, GM], "gamma/epsilon (T3GT3G')"),
]


def main():
    print("Part A: packing the known chain conformations (ideal angles)")
    print(f"{'chain':<26s} {'predicted a x b x c (A)':>26s} {'rho':>6s}   {'experiment':>22s} {'rho':>6s}")
    for poly, seq, label in CASES:
        chain = periodic_chain(poly, seq, THREE_STATE)
        t0 = time.time()
        res = pack(chain, n_random=3000, n_refine=4, rng=np.random.default_rng(0))[0]
        a, b = sorted([res.a, res.b])
        ea, eb, ec, erho = EXPERIMENTAL_CELLS[poly.name][label]
        ea, eb = sorted([ea, eb])
        print(f"{poly.name + ' ' + res.chain:<26s} {a:6.2f} x {b:6.2f} x {res.c:5.2f}{'':>4s} {res.density:6.3f}   {ea:6.2f} x {eb:6.2f} x {ec:5.2f} {erho:6.3f}   ({time.time() - t0:.0f} s)")

    print("\nPart B: full funnel for PVDF with the illustrative potential")
    cfg = PipelineConfig(polymer="pvdf", top_k_pack=2, n_random=2000, n_refine=3, refine_maxfev=1500, n_chains=2000)
    res = run_pipeline(cfg, verbose=True)
    print()
    print(res.report())


# Windows' multiprocessing spawn model re-imports this module in each worker process
# (polyfind.pipeline's ProcessPoolExecutor is used by run_pipeline); guarding the example's
# work under __main__ keeps that re-import side-effect-free.
if __name__ == "__main__":
    main()
