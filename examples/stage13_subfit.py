"""A few of stage 13's variables fitted to a few kinds of its outside rows, the rest held.

The joint fit (examples/crystal_stage.py --stage 13 --rows) moves every variable against every
row and costs a quarter of an hour an iteration.  Some variables touch few rows: the torsion
triples of the two bonds a reversed unit makes change nothing in a plain crystal, a rigid cell
or a chain's potential, only the lone chains that have those bonds, the reversed-unit
candidate's cells and the conformer set.  This fits such a subset against the rows it touches,
with the ridge of ``polyfind.crystalfit`` on the varied variables, and writes a vector the other
tools read like any crystal-stage result.

    POLYFIND_TRAINSET=... python examples/stage13_subfit.py --x0 FIT.json --vary V1_FF,V2_FF,V3_FF,V1_HH,V2_HH,V3_HH \\
        --rows lone,candidate,conformers --frozen CELLS.json --out OUT.json
"""
import argparse
import importlib.util
import json
import time
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

from polyfind import crystalfit as C

KINDS = {"rigid": "rigid", "potential": "potential", "lone": "lone", "reversed": "reversed_unit", "candidate": "candidate", "conformers": "conformers"}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--x0", required=True, help="a crystal-stage result; a 21-variable one is given the six typed torsion variables")
    ap.add_argument("--vary", required=True, help="variable names of crystalfit.VARIABLES13_TYPED, comma list")
    ap.add_argument("--rows", default="lone,candidate,conformers", help=f"kinds of rows, of {sorted(KINDS)}")
    ap.add_argument("--only", default="", help="keep only rows whose name contains one of these (comma list)")
    ap.add_argument("--frozen", default=None, help="the candidate's frozen cells (examples/chain_split.py polyfind --freeze)")
    ap.add_argument("--maxfev", type=int, default=12)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    C.BELT_GEOMETRY = "pbe"
    spec = importlib.util.spec_from_file_location("stage13_rows", Path(__file__).with_name("stage13_rows.py"))
    R = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(R)
    kinds = [k for k in args.rows.split(",") if k]
    rows = R.Rows(log=print, frozen=args.frozen, conformers="conformers" in kinds)
    x0 = np.array(json.loads(Path(args.x0).read_text())["x"], dtype=float)
    if len(x0) == 21:
        x0 = np.concatenate([x0, x0[16:19], x0[16:19]])
    names, _, bounds, _ = C._config(x0)
    vary = [v for v in args.vary.split(",") if v]
    idx = [names.index(v) for v in vary]
    only = [s for s in args.only.split(",") if s]
    history = []

    def evaluate(v):
        x = x0.copy()
        x[idx] = v
        prm = C.parameters(x)
        out = {}
        for k in kinds:
            out.update(getattr(rows, KINDS[k])(prm))
        if only:
            out = {k: r for k, r in out.items() if any(s in k for s in only)}
        ridge = np.sqrt(C.RIDGE) * C.ridge_parts(x)
        z = np.concatenate([[r[2] for r in out.values()], ridge])
        history.append((float(z @ z), x.copy(), out, float(ridge @ ridge)))
        return z

    def table(out):
        return "\n".join(f"   {k:52s} {r[0]:+9.4f}  DFT {r[1]:+9.4f}  z {r[2]:+.2f}" for k, r in out.items() if not k.startswith("rigid "))

    t0 = time.time()
    evaluate(x0[idx])
    print(f"start: objective {history[0][0]:.2f} (ridge {history[0][3]:.2f}); " + ", ".join(f"{n} {x0[i]:+.4f}" for n, i in zip(vary, idx)))
    print(table(history[0][2]), flush=True)
    steps = np.array([C.FD_STEPS[n] for n in vary])
    sol = least_squares(evaluate, x0[idx], bounds=(bounds[idx, 0], bounds[idx, 1]), diff_step=steps / np.maximum(np.abs(x0[idx]), 1e-12),
                        x_scale=np.maximum(0.3, np.abs(x0[idx])), max_nfev=args.maxfev) if args.maxfev > 0 else None
    best = min(history, key=lambda h: h[0])
    print(f"end: objective {best[0]:.2f} (ridge {best[3]:.2f}) after {len(history)} evaluations, {time.time() - t0:.0f} s; "
          + ", ".join(f"{n} {best[1][i]:+.4f}" for n, i in zip(vary, idx)))
    print(table(best[2]))
    prm = C.parameters(best[1])
    print(prm.describe())
    Path(args.out).write_text(json.dumps({
        "variables": list(names), "x": [float(v) for v in best[1]], "objective": best[0], "ridge": best[3], "varied": vary, "rows": kinds, "from": args.x0,
        "frozen": args.frozen, "evaluations": len(history), "stopped": None if sol is None else str(sol.message),
        "residuals": {k: {"polyfind": float(r[0]), "dft": float(r[1]), "z": float(r[2])} for k, r in best[2].items()}}, indent=1))
    print(f"written {args.out}")


if __name__ == "__main__":
    main()
