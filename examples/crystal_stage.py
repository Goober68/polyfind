"""Run the crystal stage of the fit (polyfind.crystalfit) and write its result.

    python examples/crystal_stage.py [--maxfev 240] [--out docs/predictions/crystal_stage_v1.json] [--coulomb ewald]

The output file is also the run's checkpoint: the best evaluation so far is written there each time
it improves (``"complete": false``), and the final write marks it complete.  A run cut short resumes
from its checkpoint with ``--x0``.
"""
import argparse
import sys

from polyfind import crystalfit as C


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--maxfev", type=int, default=240)
    ap.add_argument("--out", default="docs/predictions/crystal_stage_v1.json")
    ap.add_argument("--coulomb", default="ewald", choices=("ewald", "dsf"))
    ap.add_argument("--stage", type=int, default=1, choices=(1, 2, 3, 4, 6, 7), help="2: induced dipoles, force-shifted LJ, H-F pair term, stiffness anchors; 3: plus charge penetration; 4: fluorine belt sites; 6: the belt form with Dipole's polar Berry slope as a further anchor; 7: stage 6 plus the backbone-angle restraint's stiffness as a variable, Dipole's backbone angles as anchors (run with --geometry pbe)")
    ap.add_argument("--x0", default=None, help="JSON result file to resume from (its best x)")
    ap.add_argument("--ridge", type=float, default=None, help="weight of the ridge toward the valence values (crystalfit.RIDGE, default 1.0): a probe of how far the crystal terms would leave the valence fit if the prior were weaker")
    ap.add_argument("--method", default="powell", choices=("powell", "trf"), help="powell: bounded Powell (the stages so far); trf: bounded least squares on the same objective with a forward-difference Jacobian (one evaluation per free variable per iteration)")
    ap.add_argument("--geometry", default="textbook", choices=("textbook", "pbe"), help="stages 4-7: the chain geometry of the belt-site polymer (pbe: Dipole's PBE-D3 bond lengths and angles)")
    ap.add_argument("--pendant-tilts", action="store_true", help="evaluate the alpha-beta gap with each polymorph's CH2/CF2 groups relaxed as rigid tilts about their carbons (polyfind.pendants.tilt_relaxation; docs/CRYSTAL_STAGE.md 5.22): the relaxation DFT's geometry has and the rigid-pendant refinement lacks")
    args = ap.parse_args()
    C.BELT_GEOMETRY = args.geometry
    x0 = None
    if args.x0:
        import json
        with open(args.x0) as fh:
            x0 = json.load(fh)["x"]
    if args.ridge is not None:
        C.RIDGE = float(args.ridge)
    extra = {"pendant_tilts": True} if args.pendant_tilts else {}
    best, history = C.fit(maxfev=args.maxfev, coulomb=args.coulomb, log=sys.stdout, stage=args.stage, x0=x0, method=args.method,
                          checkpoint=args.out, **extra)
    print(f"\nbest of {len(history)} evaluations:\n" + C.table(best))
    print("parameters:", C.parameters(best.x).describe())
    C.save(best, args.out, evaluations=len(history), complete=True)
    print(f"written {args.out}")


if __name__ == "__main__":
    main()
