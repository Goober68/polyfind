"""Fit the built-in potential to first-principles torsion scans and conformers.

DESIGN.md 5.7 fitted ``SimpleFF`` to crystal data and it did not generalise: twelve
observations, several of them not independent, against five parameters, of which two were
structurally unidentifiable.  The diagnosis was about the *data* -- the quantities that
would discriminate, torsion profiles and conformer energies, were not in it.  This script
is the answer to that: 27 parameters against several hundred relative energies and a few
thousand torsional forces from PBE-D3 calculations on VDF-based oligomers, split by
chemistry rather than by frame.

    set POLYFIND_TRAINSET=...\\trainset_v1_566.xyz    (or export, on a POSIX shell)

    python examples/fit_dft.py              # fit, held-out table, acceptance tests
    python examples/fit_dft.py --preset     # score the recorded preset; no fitting
    python examples/fit_dft.py --baseline   # the acceptance tests for the unfitted potential
    python examples/fit_dft.py --sweep      # the force-weight sweep, held out
    python examples/fit_dft.py --no-accept  # skip the acceptance tests (they cost minutes)

Set ``POLYFIND_TABLE_CACHE`` to a directory to keep the packing screen tables between
runs; the acceptance tests pay 5-25 s per conformation to build them the first time.

The last block it prints is the one that matters.  The held-out energy error says whether
the fit generalises across chemistries; the three acceptance tests say whether it fixes
the failures the package already knows about, none of which the objective ever sees.  A
fit that improves the objective and fixes none of them has not earned its place.
"""
import argparse
import os
import time

import numpy as np

from polyfind import fitting as F


def show_parameters(x, reference=None):
    reference = F.REF_X0 if reference is None else reference
    print(f"\n{'parameter':<14s} {'default':>10s} {'fitted':>10s}   (bounds)")
    for name, a, b, (lo, hi) in zip(F.REF_VARIABLES, reference, x, F.REF_BOUNDS):
        at = " <- at bound" if min(abs(b - lo), abs(b - hi)) < 1e-6 else ""
        print(f"{name:<14s} {a:10.3f} {b:10.3f}   [{lo:5.2f}, {hi:5.2f}]{at}")
    p = F.reference_ff_parameters(x)
    print("\nas SimpleFF parameters for PVDF:")
    print("  " + p.describe())


def show_acceptance(label, params):
    t0 = time.time()
    result = F.acceptance_tests(params)
    print(f"\n{label}:")
    print(F.acceptance_table(result))
    n = sum(result[k] for k in ("alpha_beta_pass", "ris_pass", "alpha_antipolar_pass"))
    print(f"  {n} of 3 pass  ({time.time() - t0:.0f} s)")
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--preset", action="store_true", help="score the recorded preset, do not fit")
    ap.add_argument("--baseline", action="store_true", help="acceptance tests for the unfitted potential")
    ap.add_argument("--sweep", action="store_true", help="force-weight sweep and stop")
    ap.add_argument("--no-accept", action="store_true", help="skip the acceptance tests")
    ap.add_argument("--force-weight", type=float, default=F.REF_FORCE_WEIGHT)
    ap.add_argument("--starts", type=int, default=3)
    args = ap.parse_args()

    if args.baseline:
        show_acceptance("the illustrative potential (SimpleFF() as shipped)", F.ILLUSTRATIVE)
        return

    if not os.environ.get(F.REFERENCE_ENV):
        raise SystemExit(f"set ${F.REFERENCE_ENV} to the extended-XYZ reference set")

    frames = F.load_reference()
    train_sys, test_sys = F.split_systems(frames)
    print(f"{len(frames)} frames, {len(train_sys) + len(test_sys)} systems "
          f"(systems with an unsaturated or ring backbone are dropped; see load_reference)")
    print(f"train ({len(train_sys):2d}): " + " ".join(train_sys))
    print(f"held out ({len(test_sys):2d}): " + " ".join(test_sys))
    design = F.ReferenceDesign([f for f in frames if f.system in train_sys])
    held = F.ReferenceDesign([f for f in frames if f.system in test_sys])

    d = F.permittivity_degeneracy(design)
    print(f"\npermittivity / charge-scale degeneracy: scaling every charge by {d['scale']} and "
          f"eps_r by its square\nmoves no frame energy by more than {d['max_energy_change']:.2e} "
          f"kcal/mol and no torque by more than {d['max_torque_change']:.2e}.\n"
          "It is exact in the functional form, so relaxed off-ideal geometries do not lift it "
          "and\neps_r is held at 1 with the charges carrying the whole electrostatic scale.")

    if args.sweep:
        print(f"\nforce weight sweep (kcal^2/mol^2 of energy error per kcal^2/(mol rad)^2 of torque)")
        print(f"{'weight':>8s} {'train E':>9s} {'held-out E':>11s} {'train tau':>10s} {'held tau':>9s}")
        for w in (0.0, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0):
            fit = F.fit_reference(frames, force_weight=w, n_starts=args.starts, verbose=False)
            print(f"{w:8g} {fit.train['energy_rms']:9.3f} {fit.test['energy_rms']:11.3f} "
                  f"{fit.train['torque_rms']:10.2f} {fit.test['torque_rms']:9.2f}")
        return

    if args.preset:
        x = F.REF_FITTED_X
        print("\nscoring the recorded preset (no fitting)")
    else:
        print(f"\nfitting {len(F.REF_VARIABLES)} parameters, force weight {args.force_weight}")
        fit = F.fit_reference(frames, force_weight=args.force_weight, n_starts=args.starts)
        x = fit.x
        print(f"{fit.n_evaluations} objective evaluations, {fit.seconds:.0f} s")

    print(F.reference_table(x, design, held))
    tr, te = design.errors(x), held.errors(x)
    tr0, te0 = design.errors(F.REF_X0), held.errors(F.REF_X0)
    print(f"\nenergy RMS (kcal/mol):  train {tr0['energy_rms']:.3f} -> {tr['energy_rms']:.3f}, "
          f"HELD OUT {te0['energy_rms']:.3f} -> {te['energy_rms']:.3f}")
    print(f"torque RMS (kcal/mol/rad): train {tr0['torque_rms']:.2f} -> {tr['torque_rms']:.2f}, "
          f"HELD OUT {te0['torque_rms']:.2f} -> {te['torque_rms']:.2f}")
    print("Target: about 1 kcal/mol on conformer energies would be a good classical potential.")
    show_parameters(x)

    if not args.no_accept:
        show_acceptance("the illustrative potential (SimpleFF() as shipped)", F.ILLUSTRATIVE)
        show_acceptance("the fit to first-principles data", F.reference_ff_parameters(x))


if __name__ == "__main__":
    main()
