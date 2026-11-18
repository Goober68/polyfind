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
    ap.add_argument("--stage", type=int, default=1, choices=(1, 2, 3, 4, 6, 7, 11, 12, 13), help="13: stage 12 with a charge site on every C-H bond and a charge drawn along the backbone, the chains built with the sites (crystalfit.STAGE13); 12: stage 11 with the buffered 14-7 wall in place of 12-6, its first buffer a variable (crystalfit.STAGE12); 11: stage 7 plus the cost of moving chains along each other against DFT (crystalfit.REGISTRY); 2: induced dipoles, force-shifted LJ, H-F pair term, stiffness anchors; 3: plus charge penetration; 4: fluorine belt sites; 6: the belt form with Dipole's polar Berry slope as a further anchor; 7: stage 6 plus the backbone-angle restraint's stiffness as a variable, Dipole's backbone angles as anchors (run with --geometry pbe)")
    ap.add_argument("--x0", default=None, help="JSON result file to resume from (its best x)")
    ap.add_argument("--ridge", type=float, default=None, help="weight of the ridge toward the valence values (crystalfit.RIDGE, default 1.0): a probe of how far the crystal terms would leave the valence fit if the prior were weaker")
    ap.add_argument("--method", default="powell", choices=("powell", "trf"), help="powell: bounded Powell (the stages so far); trf: bounded least squares on the same objective with a forward-difference Jacobian (one evaluation per free variable per iteration)")
    ap.add_argument("--geometry", default="textbook", choices=("textbook", "pbe"), help="stages 4-7: the chain geometry of the belt-site polymer (pbe: Dipole's PBE-D3 bond lengths and angles)")
    ap.add_argument("--pendant-tilts", action="store_true", help="evaluate the alpha-beta gap with each polymorph's CH2/CF2 groups relaxed as rigid tilts about their carbons (polyfind.pendants.tilt_relaxation; docs/CRYSTAL_STAGE.md 5.22): the relaxation DFT's geometry has and the rigid-pendant refinement lacks")
    ap.add_argument("--set", default="", help="start values over --x0's, name=value comma list (a stage's variable names)")
    ap.add_argument("--wide", action="store_true", help="stage 11: the diameters' bounds opened (crystalfit.BOUNDS11_WIDE)")
    ap.add_argument("--rows", action="store_true", help="stage 13: the rows from outside the module (examples/stage13_rows.py): the rigid cells, "
                    "the single-chain potentials and the conformers")
    ap.add_argument("--frozen", default=None, help="stage 13 with --rows: the reversed-unit candidate's relaxed cells (examples/chain_split.py polyfind "
                    "--freeze), on which its rows are evaluated; without it those rows are left out")
    ap.add_argument("--typed", action="store_true", help="stage 13: the CF2-CF2 and CH2-CH2 bonds' own torsion triples as six more variables "
                    "(crystalfit.VARIABLES13_TYPED), started at the CF2-CH2 bond's")
    args = ap.parse_args()
    C.BELT_GEOMETRY = args.geometry
    x0 = None
    if args.x0:
        import json
        with open(args.x0) as fh:
            x0 = json.load(fh)["x"]
    if args.set:
        import numpy as np
        names = {7: C.VARIABLES7, 11: C.VARIABLES7, 12: C.VARIABLES12, 13: C.VARIABLES13_TYPED}[args.stage]
        base = list(x0) if x0 is not None else list({7: C.X0_7, 11: C.X0_7, 12: C.X0_12, 13: C.X0_13}[args.stage])
        if args.stage in (12, 13) and len(base) == 13:
            base.append(float(C.X0_12[13]))
        if args.stage == 13 and len(base) == 14:
            base += [0.0, 0.0]
        if args.stage == 13 and len(base) == 16:
            base += [float(v) for v in C.X0_13[16:19]]
        if args.stage == 13 and len(base) == 19:
            base += [1.0, 1.0]
        if args.typed and len(base) == 21:
            base += base[16:19] * 2
        for item in args.set.split(","):
            k, v = item.split("=")
            base[names.index(k)] = float(v)
        x0 = np.array(base)
    elif args.typed and x0 is not None and len(x0) == 21:
        x0 = list(x0) + list(x0[16:19]) * 2
    if args.ridge is not None:
        C.RIDGE = float(args.ridge)
    extra = {"pendant_tilts": True} if args.pendant_tilts else {}
    if args.wide:
        extra["bounds"] = C.BOUNDS11_WIDE
    if args.rows:
        import importlib.util
        from pathlib import Path

        spec = importlib.util.spec_from_file_location("stage13_rows", Path(__file__).with_name("stage13_rows.py"))
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        extra["extra_rows"] = mod.Rows(log=print, frozen=args.frozen)
    best, history = C.fit(maxfev=args.maxfev, coulomb=args.coulomb, log=sys.stdout, stage=args.stage, x0=x0, method=args.method,
                          checkpoint=args.out, **extra)
    print(f"\nbest of {len(history)} evaluations:\n" + C.table(best))
    print("parameters:", C.parameters(best.x).describe())
    C.save(best, args.out, evaluations=len(history), complete=True)
    print(f"written {args.out}")


if __name__ == "__main__":
    main()
