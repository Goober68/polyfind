"""The intrinsic electrostriction of beta-PVDF and the two-state switch of alpha to beta.

    python examples/electrostriction.py                 # illustrative potential, rigid chain, fixed charges
    python examples/electrostriction.py --polarizable   # fitted potential, Ewald, induced dipoles (slower)

Prints the field electrostriction tensor by the finite-field and Maxwell routes
(``docs/ELECTROSTRICTION_RESEARCH.md`` tier 1), the odd part against the piezoelectric
module, the polarization form ``Q``, and then the geometry and threshold of the
TGTG'-to-TTTT switch with the alpha-beta margin of the DFT fit as ``dE``.
"""
import argparse
import time

import numpy as np

from polyfind import electrostriction as ES
from polyfind import mechanics as M
from polyfind.polymers import PE, PVDF, THREE_STATE

T, GP, GM = 0, 1, 2


def build(polymer, seq, label, polarizable: bool):
    kw = {"pack_kw": {"table_cache_dir": None}, "lj_cutoff": "force"}  # see electrostriction() on the cutoff
    if polarizable:
        from polyfind.ewald import EwaldSpec
        from polyfind.polarizability import Polarizable
        kw.update(coulomb="ewald", ewald=EwaldSpec(), polarizable=Polarizable())
    ref, rr = M.refined_reference(polymer, seq, label=label, **kw)
    return M.relax_reference(ref), rr


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--polarizable", action="store_true")
    ap.add_argument("--field", type=float, default=0.005, help="V/A; 0.005 = 50 MV/m")
    args = ap.parse_args()
    if args.polarizable:
        from polyfind.fitting import FITTED_VALENCE
        ctx = FITTED_VALENCE.applied()
    else:
        import contextlib
        ctx = contextlib.nullcontext()
    with ctx:
        for polymer, seq, label in [(PVDF, [T, T], "beta-PVDF"), (PE, [T], "PE (null control)")]:
            t0 = time.time()
            ref, _ = build(polymer, seq, label, args.polarizable)
            el = M.elastic_constants(ref)
            es = ES.electrostriction(ref, el, field=args.field)
            pz = M.piezoelectric(ref, el, field=args.field)
            print(f"\n--- {label}: a={ref.params[0]:.3f} b={ref.params[1]:.3f} gamma={ref.params[2]:.2f} "
                  f"c={ref.c:.4f}  ({time.time() - t0:.1f} s)")
            print(es.table())
            idx = list(es.reachable)
            print("odd part / E (pm/V) against piezoelectric.d_direct, max |diff|: "
                  f"{np.nanmax(np.abs(es.linear[:, idx] - pz.d_direct)):.3g} of max |d| {np.nanmax(np.abs(pz.d_direct)):.3g}")
            if polymer is PVDF:
                from polyfind.helix import helix_parameters
                P = ref.polarization()
                hp = helix_parameters(PVDF, np.array(seq), THREE_STATE)
                n_chains = getattr(ref.packer, "n_chains", 2)
                n_mon = hp.monomers_per_period * hp.periods_per_repeat * n_chains
                dmu = ES.dipole_per_monomer(abs(P[0]), ref.volume, n_mon)
                ts = ES.transition_strain(PVDF, [T, GP, T, GM], [T, T], THREE_STATE)
                print(f"\nswitch TGTG' -> TTTT: rise per monomer {ts['rise_ground']:.3f} -> {ts['rise_polar']:.3f} A, "
                      f"axial strain {ts['axial'] * 100:+.1f}%")
                print(f"beta polarization {abs(P[0]):.3f} C/m^2 over {n_mon} monomers per cell -> {dmu:.3f} e.A per monomer")
                print("  dE (kcal/mol/monomer)  threshold (MV/m)   switched fraction at 100 MV/m, 300 K: n=1 / n=20")
                for de in (1.08, 0.27, 0.05):
                    one = ES.TwoState(delta_e=de, delta_mu=dmu, strain=ts["axial"], n=1)
                    coop = ES.TwoState(delta_e=de, delta_mu=dmu, strain=ts["axial"], n=20)
                    print(f"  {de:8.2f}             {one.threshold * 1e4:10.1f}          "
                          f"{one.fraction(0.01):.3g} / {coop.fraction(0.01):.3g}")


if __name__ == "__main__":
    main()
