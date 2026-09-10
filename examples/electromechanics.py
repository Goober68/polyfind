"""Elastic, piezoelectric and actuator response of the reference polymorphs.

Runs :mod:`polyfind.mechanics` on PVDF beta, alpha and gamma and on polyethylene as a null
control, and prints the tables ``docs/ELECTROMECHANICS.md`` records.

    python examples/electromechanics.py                        # every potential, all four cases
    python examples/electromechanics.py --preset illustrative  # the built-in potential only
    python examples/electromechanics.py --preset valence       # the one preset with valence terms
    python examples/electromechanics.py --axial                # add the rigid axial diagnosis
    python examples/electromechanics.py --stiffness-sweep      # is C_33 still an invented constant?
    python examples/electromechanics.py --cutoff               # the Lennard-Jones cutoff sweep

The three potentials are not interchangeable and the table says which is which:

* ``illustrative`` and ``pvdf-dft-fit`` are **rigid**: no bond or angle terms anywhere, so
  the chain cannot deform, ``C_33`` is ``nan`` and every diagonal column of ``e`` is zero.
* ``pvdf-dft-valence`` carries fitted stretch and bend terms.  They are forwarded into the
  lattice kernel (``CrystalPacker(valence=...)``, opt-in), which gives the chain a restoring
  force and makes ``eps_zz`` and the diagonal columns computable.  That path costs a
  constrained relaxation per strain state and is a few seconds per case rather than a
  fraction of one.

A few seconds per case for the rigid runs; beta 5 s, alpha 20 s and gamma 100 s for the
deformable one, the packing and refinement included.
"""
import argparse
import sys
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
VOIGT = {0: "1", 1: "2", 2: "3", 5: "6"}


def build(polymer, seq, label, valence=None, angle_stiffness=None):
    """Pack, refine and wrap one polymorph; returns ``(reference, refine result, shape, seconds)``."""
    t0 = time.time()
    refine_kw = None
    if valence is not None or angle_stiffness is not None:
        refine_kw = {"valence": valence}
        if angle_stiffness is not None:
            refine_kw["angle_stiffness"] = angle_stiffness
    ref, rr = M.refined_reference(polymer, seq, label=label, valence=valence, refine_kw=refine_kw)
    shape = None if valence is None else M.shape_of(polymer, ref, angles=rr.angles)
    return ref, rr, shape, time.time() - t0


def show(resp, rr, t_build):
    ref = resp.reference
    print()
    print(f"--- {resp.label}   (pack+refine {t_build:.1f} s)")
    print(f"    refined cell a={ref.params[0]:.3f} b={ref.params[1]:.3f} gamma={ref.params[2]:.2f} "
          f"c={ref.c:.4f}  residual sigma (GPa) = {np.round(resp.elastic.residual_stress, 6)}")
    if resp.shape is not None:
        print(f"    shape parameters ({resp.shape.n}): {resp.shape.labels()}")
    print("    " + resp.summary().replace("\n", "\n    "))
    if resp.piezo.d_improper.size:
        print("    d_improper (pC/N; NOT a piezoelectric constant -- the dimensional term):")
        for i, r in enumerate("xyz"):
            print(f"      {r}  " + "  ".join(f"{v:+9.3f}" for v in resp.piezo.d_improper[i]))
    if resp.axial is not None:
        ax = resp.axial
        print(f"    axial (rigid diagnosis): eps_zz reachable over {ax.strain_range[0] * 100:+.2f}% .. "
              f"{ax.strain_range[1] * 100:+.2f}% via the backbone angles; "
              f"sigma_zz(rigid chain) = {ax.sigma_zz_rigid:+.3f} GPa")
        print("    axial: C_33 (GPa) against the invented bend constant k (kcal/mol/rad^2): "
              + ", ".join(f"k={k:g}: {v:.1f}" for k, v in sorted(ax.c33.items())))
        print("    axial: sigma_zz (GPa) along the same path:                                "
              + ", ".join(f"k={k:g}: {v:+.2f}" for k, v in sorted(ax.sigma_zz.items())))
        print(f"    axial verdict: computable={ax.computable}; the invented bend constant supplies "
              f"{ax.bend_fraction * 100:.0f}% of C_33 at k=105")
    sys.stdout.flush()


def stiffness_sweep(valence, ks=(0.0, 52.5, 105.0, 210.0), step=2e-3):
    """C_33 against ``refine_crystal``'s invented bend constant, on the deformable path.

    The test that the number has stopped being a report of that constant.  ``k`` moves the
    structure the refinement hands over; the relaxation against the *fitted* valence terms
    then has to put it back, and what is printed is whether it does.
    """
    print("=" * 100)
    print("C_33 (GPa) against refine_crystal's invented angle_stiffness k, valence path")
    print("=" * 100)
    for polymer, seq, label in CASES:
        row = []
        for k in ks:
            ref, rr, shape, _ = build(polymer, seq, label, valence=valence, angle_stiffness=k)
            ref, shape = M.relax_reference_deformable(ref, shape)
            el = M.elastic_constants(ref, step=step, shape=shape)
            row.append((k, rr.angles.copy(), rr.result.c, ref.c, el.C[2, 2], el.c33_from_energy))
        print(f"  {label}")
        for k, angles, c_ref, c_rel, c33, c33e in row:
            print(f"    k={k:6.1f}  refined angles {np.round(angles, 3)}  c_refined={c_ref:.5f} -> "
                  f"c_relaxed={c_rel:.5f}   C_33={c33:9.4f}  (energy route {c33e:9.4f})")
        v = np.array([r[4] for r in row])
        print(f"    spread over the sweep: {v.max() - v.min():.3e} GPa "
              f"({(v.max() - v.min()) / v.mean() * 100:.2e}% of the mean)")
        sys.stdout.flush()


def cutoff_sweep(cases=("PVDF beta (TTTT)",), cutoffs=(8.0, 10.0, 12.0, 16.0, 20.0),
                 steps=(5e-4, 2e-3, 4e-3)):
    """The Lennard-Jones cutoff, both shift forms, against the elastic constants."""
    print("=" * 100)
    print("Lennard-Jones cutoff and shift form (illustrative potential, rigid chain)")
    print("=" * 100)
    print(f"{'case':22s} {'shift':7s} {'rc':>4s} {'a':>8s} {'b':>8s} {'E/mon':>9s} {'C11':>8s} "
          f"{'C22':>8s} {'C12':>8s} {'C66':>8s} {'C22 spread':>11s}")
    for polymer, seq, label in CASES:
        if cases and label not in cases:
            continue
        ref, rr, _, _ = build(polymer, seq, label)
        for shift in ("energy", "force"):
            for rc in cutoffs:
                r = M.relax_reference(M.reference_from_chain(ref.packer.chain, ref.params, n_chains=2,
                                                             cutoff=rc, lj_cutoff=shift))
                c22 = [M.elastic_constants(r, step=s).C[1, 1] for s in steps]
                el = M.elastic_constants(r, step=steps[len(steps) // 2])
                emon = float(r.packer.energy(r.params[None])[0]) / (2 * r.packer.chain.n_monomers)
                print(f"{label:22s} {shift:7s} {rc:4.0f} {r.params[0]:8.4f} {r.params[1]:8.4f} "
                      f"{emon:9.4f} {el.C[0, 0]:8.3f} {el.C[1, 1]:8.3f} {el.C[0, 1]:8.3f} "
                      f"{el.C[5, 5]:8.3f} {max(c22) - min(c22):11.3f}")
                sys.stdout.flush()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--preset", default="all",
                    choices=["all", "illustrative", "fitted", "valence"],
                    help="which potential(s) to run")
    ap.add_argument("--axial", action="store_true", help="also run the rigid axial diagnosis")
    ap.add_argument("--stiffness-sweep", action="store_true",
                    help="sweep refine_crystal's invented bend constant on the deformable path")
    ap.add_argument("--cutoff", action="store_true", help="run the Lennard-Jones cutoff sweep")
    args = ap.parse_args()

    from polyfind.forcefield import SimpleFF

    runs = []
    if args.preset in ("all", "illustrative"):
        runs.append(("illustrative", None, None))
    if args.preset in ("all", "fitted"):
        from polyfind.fitting import FITTED_DFT

        runs.append(("pvdf-dft-fit", FITTED_DFT, None))
    if args.preset in ("all", "valence"):
        from polyfind.fitting import FITTED_VALENCE

        runs.append(("pvdf-dft-valence", FITTED_VALENCE, SimpleFF.from_preset("pvdf-dft-valence")))

    if args.cutoff:
        cutoff_sweep(cases=tuple(label for _, _, label in CASES))
        return

    for pot_name, params, valence in runs:
        print("=" * 100)
        print(f"potential: {pot_name}"
              + ("   [valence terms forwarded into the lattice kernel: the chain can deform]"
                 if valence is not None else "   [rigid chain]"))
        if params is not None:
            print("  " + params.describe())
        print("=" * 100)
        ctx = params.applied() if params is not None else None
        if ctx is not None:
            ctx.__enter__()
        try:
            if args.stiffness_sweep and valence is not None:
                stiffness_sweep(valence)
                continue
            for polymer, seq, label in CASES:
                ref, rr, shape, t_build = build(polymer, seq, label, valence=valence)
                resp = M.electromechanical_response(ref, polymer=polymer, axial=args.axial, shape=shape)
                show(resp, rr, t_build)
        finally:
            if ctx is not None:
                ctx.__exit__(None, None, None)


if __name__ == "__main__":
    main()
