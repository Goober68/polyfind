"""Fit the potential with bond and angle terms, full Cartesian forces and an off-site charge.

``examples/fit_dft.py`` fitted the rigid form -- torsion, Lennard-Jones, atom-centred point
charges -- to first-principles conformers, and DESIGN.md 5.9 records where that left
things: the fit generalised, one acceptance test of three was fixed, and the binding
constraint had moved from the data to the *functional form*.  Two gaps were named.  It had
no bond or angle terms, so 93% of the reference force signal was unrepresentable and an
all-trans chain could not relieve a contact by opening an angle; and its electrostatics was
fixed atom-centred point charges, which the fit tried to escape by driving fluorine's charge
positive.  This script closes both and re-measures everything that was measured before.

    set POLYFIND_TRAINSET=...\\trainset_v1_566.xyz    (or export, on a POSIX shell)

    python examples/fit_valence.py              # fit, held-out table, acceptance tests
    python examples/fit_valence.py --preset     # score the recorded preset; no fitting
    python examples/fit_valence.py --sweep      # the force-weight sweep, held out
    python examples/fit_valence.py --strain     # the all-trans strain, angles free (no data needed)
    python examples/fit_valence.py --no-accept  # skip the acceptance tests (they cost a minute)

Set ``POLYFIND_TABLE_CACHE`` to a directory to keep the packing screen tables between runs.

What to read.  The held-out energy error says whether the fit generalises across chemistries
and the held-out *force* error says whether the forces are finally being used; the three
acceptance tests, none of which the objective sees, say whether any of it fixed a failure the
package already knew about.  docs/VALENCE_FIT.md has the numbers this script reproduces and
is blunt about which of them did not move.
"""
import argparse
import os
import time

import numpy as np

from polyfind import fitting as F
from polyfind.forcefield import SimpleFF, relax_backbone_angles
from polyfind.polymers import get_polymer

STRAIN_POLYMERS = ("pvdf", "pvdc", "cfe", "cdfe", "an", "vdcn", "fanome")


def show_parameters(x):
    print(f"\n{'parameter':<14s} {'start':>10s} {'fitted':>10s}   (bounds)")
    for name, a, b, (lo, hi) in zip(F.VAL_VARIABLES, F.VAL_X0, x, F.VAL_BOUNDS):
        at = " <- at bound" if min(abs(b - lo), abs(b - hi)) < 1e-6 * (hi - lo) else ""
        print(f"{name:<14s} {a:10.3f} {b:10.3f}   [{lo:7.2f}, {hi:7.2f}]{at}")


def show_acceptance(label, params):
    t0 = time.time()
    result = F.acceptance_tests(params)
    print(f"\n{label}:")
    print(F.acceptance_table(result))
    n = sum(result[k] for k in ("alpha_beta_pass", "ris_pass", "alpha_antipolar_pass"))
    print(f"  {n} of 3 pass  ({time.time() - t0:.0f} s)")
    return result


def show_strain(params):
    """All-trans strain with the backbone angles frozen and with them free.

    The convention is the one ``polymers.py`` and ``tests/test_pvdc.py`` already use: the
    Lennard-Jones energy of a ten-bond all-trans oligomer.  ``polyfind.polymers`` records
    9 kcal/mol for PVDF against 313 for PVDC and about 1e6 for FANOME, and concludes that
    for anything past PVDF the all-trans RIS reference is a state the polymer never
    occupies.  With an angle term the angles can open, which is what a real chain does.
    """
    ff = params.simple_ff()
    plain = SimpleFF()
    print(f"\n{'polymer':<9s} {'LJ rigid':>10s} {'LJ rigid':>10s} {'LJ relaxed':>11s} "
          f"{'E relaxed-rigid':>16s}   backbone angles (deg)")
    print(f"{'':<9s} {'(unfitted)':>10s} {'(fitted)':>10s}")
    for name in STRAIN_POLYMERS:
        p = get_polymer(name)
        dih = np.full(10, 180.0)
        r = relax_backbone_angles(p, ff, dih)
        from polyfind.chain import build_chain
        lj0 = plain.components(build_chain(p, dih))["lj"]
        angles = ", ".join(f"{a:.1f}" for a in r["angles"])
        print(f"{name:<9s} {lj0:10.1f} {r['lj0']:10.1f} {r['lj']:11.1f} "
              f"{r['energy'] - r['energy0']:16.1f}   {angles}")
    print("PVDC's measured backbone angles are 123 deg at CH2 and 114 at CCl2 "
          "(docs/REFERENCES.md).")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--preset", action="store_true", help="score the recorded preset, do not fit")
    ap.add_argument("--sweep", action="store_true", help="force-weight sweep and stop")
    ap.add_argument("--strain", action="store_true", help="the strain table only; needs no data")
    ap.add_argument("--no-accept", action="store_true", help="skip the acceptance tests")
    ap.add_argument("--force-weight", type=float, default=F.VAL_FORCE_WEIGHT)
    ap.add_argument("--starts", type=int, default=3)
    args = ap.parse_args()

    if args.strain:
        show_strain(F.FITTED_VALENCE)
        return

    if not os.environ.get(F.REFERENCE_ENV):
        raise SystemExit(f"set ${F.REFERENCE_ENV} to the extended-XYZ reference set")

    frames = F.load_reference()
    train_sys, test_sys = F.split_systems(frames)
    design = F.ValenceDesign([f for f in frames if f.system in train_sys])
    held = F.ValenceDesign([f for f in frames if f.system in test_sys])
    print(f"{len(frames)} frames over {len(train_sys) + len(test_sys)} systems; "
          f"train {design.n_frames} frames / {len(train_sys)} systems, "
          f"held out {held.n_frames} / {len(test_sys)}")
    print(f"{design.bond_r.size} bonds and {design.angle_theta.size} angles carry a valence term, "
          f"{3 * design.n_atoms} force components enter the objective")

    d = F.valence_degeneracies(design)
    print(f"\ndegeneracies, measured rather than asserted:\n"
          f"  charge scale against permittivity: {d['charge_permittivity_max_energy_change']:.1e} "
          f"kcal/mol -- still exact, so eps_r stays at 1\n"
          f"  fluorine charge against its offset at fixed dipole: "
          f"{d['charge_offset_max_energy_change']:.2f} kcal/mol -- not a degeneracy\n"
          f"  largest angle carrying a bend term: {d['max_bent_angle_deg']:.0f} deg "
          f"(a harmonic is never asked to work near 180)")

    if args.sweep:
        print("\nforce weight sweep (kcal^2/mol^2 of energy error per kcal^2/(mol A)^2 of force)")
        print(f"{'weight':>8s} {'train E':>9s} {'held E':>9s} {'train F':>9s} {'held F':>9s} {'r':>7s}")
        for w in (0.0, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0):
            fit = F.fit_valence(designs=(design, held), force_weight=w, n_starts=args.starts,
                                verbose=False)
            print(f"{w:8g} {fit.train['energy_rms']:9.3f} {fit.test['energy_rms']:9.3f} "
                  f"{fit.train['force_rms']:9.2f} {fit.test['force_rms']:9.2f} "
                  f"{fit.test['force_corr']:7.3f}")
        print(f"the reference forces themselves have RMS "
              f"{held.errors(F.VAL_X0)['force_ref_rms']:.2f} kcal/(mol A), which is what a model "
              f"predicting zero would score")
        return

    if args.preset:
        x = F.VAL_FITTED_X
        print("\nscoring the recorded preset (no fitting)")
    else:
        print(f"\nfitting {len(F.VAL_VARIABLES)} parameters, force weight {args.force_weight}")
        fit = F.fit_valence(designs=(design, held), force_weight=args.force_weight,
                            n_starts=args.starts)
        x = fit.x
        print(f"{fit.n_evaluations} objective evaluations, {fit.seconds:.0f} s")
        print(fit.table())

    for label, vec in (("illustrative, rigid", F.rigid_vector(F.REF_X0)),
                       ("pvdf-dft-fit, rigid", F.rigid_vector(F.REF_FITTED_X)),
                       ("this fit", x)):
        tr, te = design.errors(vec), held.errors(vec)
        print(f"{label:<22s} energy {tr['energy_rms']:.3f} / {te['energy_rms']:.3f}   "
              f"force {tr['force_rms']:6.2f} / {te['force_rms']:6.2f}   "
              f"torque {te['torque_rms']:6.2f}   (train / held out)")
    show_parameters(x)
    show_strain(F.valence_ff_parameters(x))

    if not args.no_accept:
        show_acceptance("pvdf-dft-fit (the previous fit)", F.FITTED_DFT)
        show_acceptance("this fit", F.valence_ff_parameters(x))


if __name__ == "__main__":
    main()
