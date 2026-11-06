"""Run the crystal stage of the fit (polyfind.crystalfit) and write its result.

    python examples/crystal_stage.py [--maxfev 240] [--out docs/predictions/crystal_stage_v1.json] [--coulomb ewald]
"""
import argparse
import sys

from polyfind import crystalfit as C


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--maxfev", type=int, default=240)
    ap.add_argument("--out", default="docs/predictions/crystal_stage_v1.json")
    ap.add_argument("--coulomb", default="ewald", choices=("ewald", "dsf"))
    ap.add_argument("--stage", type=int, default=1, choices=(1, 2, 3, 4), help="2: induced dipoles, force-shifted LJ, H-F pair term, stiffness anchors; 3: plus charge penetration; 4: fluorine belt sites")
    ap.add_argument("--x0", default=None, help="JSON result file to resume from (its best x)")
    args = ap.parse_args()
    x0 = None
    if args.x0:
        import json
        with open(args.x0) as fh:
            x0 = json.load(fh)["x"]
    best, history = C.fit(maxfev=args.maxfev, coulomb=args.coulomb, log=sys.stdout, stage=args.stage, x0=x0)
    print(f"\nbest of {len(history)} evaluations:\n" + C.table(best))
    print("parameters:", C.parameters(best.x).describe())
    C.save(best, args.out)
    print(f"written {args.out}")


if __name__ == "__main__":
    main()
