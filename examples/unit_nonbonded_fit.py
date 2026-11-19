"""A unit's own non-bonded terms, fitted to DFT's margins of the smallest cell that holds the unit, on top of a preset.

examples/unit_dft_check.py splits a unit's error into what is inside a chain and what is between
chains.  Inside the chain it is the unit's bonds (examples/unit_torsion_fit.py).  Between chains
it is the unit's own atoms: for the nitrile, nitrogen's wall, which every preset has at its UFF
values, and the charge the C-N bond moves, which was never a target.  This fits those to DFT's
all-trans-minus-helix energies of the cell, one for each helix placement DFT was asked for
(``docs/predictions/row_cells_dft_ROW.json``), polyfind's cells relaxed again at every step with
the shear free and the pendants tilted (examples/rows_under_preset.py, two minutes a step).

    python examples/unit_nonbonded_fit.py --base pvdf-dft-chain1-cfe --row vdf2-trfe1-an1-4 --increment C-N --element N \\
        --start C-N=0.75 --work DIR --out docs/predictions/unit_nonbonded_an_chain1.json

The output names its base and the terms to put on it: every tool that takes a preset takes it
(``polyfind.fitting.parameters_of``).  A fit to one cell is not a test on it: the test is a second
cell of the unit that the fit has not seen.
"""
import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

from polyfind import fitting as F
from polyfind.polymers import UFF_LJ

HERE = Path(__file__).resolve().parent
PRED = HERE.parent / "docs" / "predictions"
TOL = 0.05                          # kcal/mol per monomer: a row's tolerance
STEP = {"q": 0.1, "x": 0.15, "lnD": 0.3}   # one unit of each scaled variable: e, A, and a factor of 1.35 in the depth
RIDGE = {"q": 1.0, "x": 0.4, "lnD": float(np.log(2.0))}   # how far each may sit from its start before it costs one row's tolerance
BOUNDS = {"q": (-1.5, 1.5), "x": (2.6, 4.6), "lnD": (float(np.log(0.005)), float(np.log(0.5)))}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--base", required=True, help="a preset's name, or a file polyfind.fitting.parameters_of reads")
    ap.add_argument("--row", required=True)
    ap.add_argument("--increment", action="append", default=[], help="a bond whose charge increment is fitted, as C-N")
    ap.add_argument("--element", action="append", default=[], help="an element whose wall (distance and depth) is fitted")
    ap.add_argument("--start", action="append", default=[], help="NAME=VALUE: start a variable away from the base's value (C-N=0.75, x_N=3.5, D_N=0.05)")
    ap.add_argument("--dft", type=Path, default=None, help="default docs/predictions/row_cells_dft_ROW.json")
    ap.add_argument("--polished", type=Path, default=PRED / "polished_margins_belt4_dft_rows.json")
    ap.add_argument("--references", type=Path, default=PRED / "dilute_kink_belt4_dft_rows.json")
    ap.add_argument("--work", type=Path, required=True, help="a folder for each step's parameter file and cells")
    ap.add_argument("--maxfev", type=int, default=6, help="steps of the least-squares search; each costs one evaluation and one more a variable")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    args.work.mkdir(parents=True, exist_ok=True)
    base = F.parameters_of(args.base)
    dft = json.loads((args.dft or PRED / f"row_cells_dft_{args.row}.json").read_text())
    target = {dft["lines"][m["helix"]]["conformation"]: float(m["dft"]) for m in dft["margins"].values()}
    inc0 = {(a, b): q for a, b, q in base.charge_increments}
    lj0 = {e: (x, d) for e, x, d in base.lj}
    names, kinds, start = [], [], []
    for pair in args.increment:
        a, b = pair.split("-")
        names.append(pair), kinds.append("q"), start.append(inc0[(a, b)] if (a, b) in inc0 else -inc0.get((b, a), 0.0))
    for e in args.element:
        x, d = lj0.get(e, UFF_LJ[e])
        names.extend([f"x_{e}", f"D_{e}"]), kinds.extend(["x", "lnD"]), start.extend([float(x), float(np.log(d))])
    home = np.array(start, dtype=float)       # the base's own values: where the ridge holds
    first = home.copy()
    for item in args.start:
        name, value = item.split("=")
        k = names.index(name)
        first[k] = np.log(float(value)) if kinds[k] == "lnD" else float(value)
    scale = np.array([STEP[k] for k in kinds])
    history = []

    def spec_of(p):
        out = {"base": args.base, "charge_increments": {}, "lj": {}}
        for name, kind, v in zip(names, kinds, p):
            if kind == "q":
                out["charge_increments"][name] = float(v)
            elif kind == "x":
                out["lj"][name[2:]] = [float(v), None]
            else:
                out["lj"][name[2:]][1] = float(np.exp(v))
        return out

    def cells_at(p, tag):
        spec, out = args.work / f"{tag}.json", args.work / f"{tag}.margins.json"
        spec.write_text(json.dumps(spec_of(p)))
        if not out.exists():
            subprocess.run([sys.executable, "-W", "ignore", str(HERE / "rows_under_preset.py"), "--preset", str(spec), "--rows", args.row,
                            "--polished", str(args.polished), "--references", str(args.references), "--out", str(out)], stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT, check=True)
        entry = json.loads(out.read_text())["rows"][args.row]
        e = {c: v["sheared"] + v.get("tilt_sheared", 0.0) for c, v in entry["cells"].items()}
        return {c: e["TTTT"] - e[c] for c in e if c != "TTTT"}, entry["margin"]

    def evaluate(u):
        p = home + (np.asarray(u) - 1.0) * scale
        mine, margin = cells_at(p, f"step_{len(history):03d}")
        rows = np.array([(mine[c] - target[c]) / TOL for c in target])
        ridge = np.array([(v - h) / RIDGE[k] for v, h, k in zip(p, home, kinds)])
        z = np.concatenate([rows, ridge])
        history.append((float(z @ z), p.copy(), mine, margin, float(ridge @ ridge)))
        print(f"  {len(history):3d}: " + ", ".join(f"{n} {np.exp(v) if k == 'lnD' else v:+.4f}" for n, k, v in zip(names, kinds, p)) + f": objective {z @ z:8.2f}; "
              + ", ".join(f"{c} {mine[c]:+.3f} ({target[c]:+.3f})" for c in target), flush=True)
        return z

    def show(h):
        return ", ".join(f"{c} {h[2][c]:+.3f}" + (f" (DFT {target[c]:+.3f})" if c in target else "") for c in h[2]) + f"; margin {h[3]:+.3f}"

    t0 = time.time()
    print(f"{args.base}, {args.row}; DFT all-trans minus helix: " + ", ".join(f"{c} {v:+.3f}" for c, v in target.items()))
    evaluate(np.ones(len(names)))
    at_base = history[0]
    u0 = 1.0 + (first - home) / scale
    lo = 1.0 + (np.array([BOUNDS[k][0] for k in kinds]) - home) / scale
    hi = 1.0 + (np.array([BOUNDS[k][1] for k in kinds]) - home) / scale
    sol = least_squares(evaluate, u0, bounds=(lo, hi), diff_step=0.5, x_scale=1.0, max_nfev=args.maxfev) if args.maxfev > 0 else None
    best = min(history, key=lambda h: h[0])
    print(f"at the base: {show(at_base)}; objective {at_base[0]:.2f}")
    print(f"fitted:      {show(best)}; objective {best[0]:.2f} (ridge {best[4]:.2f}) after {len(history)} evaluations, {time.time() - t0:.0f} s")
    spec = spec_of(best[1])
    out = {**spec, "what": "non-bonded terms of a unit's own atoms to put on the base (charge increments in e; lj: the wall's distance in A and depth in kcal/mol), "
                           f"fitted to DFT's all-trans-minus-helix energies of {args.row} (examples/unit_nonbonded_fit.py)",
           "row": args.row, "variables": names, "at_the_base_values": {n: float(np.exp(v) if k == "lnD" else v) for n, k, v in zip(names, kinds, home)},
           "objective": best[0], "objective_at_the_base": at_base[0], "ridge": best[4], "evaluations": len(history), "stopped": None if sol is None else str(sol.message),
           "margin": best[3], "margin_at_the_base": at_base[3],
           "rows": {c: {"polyfind": float(best[2][c]), "dft": target.get(c), "at_the_base": float(at_base[2][c]), "fitted": c in target} for c in best[2]}}
    args.out.write_text(json.dumps(out, indent=1))
    print(f"written {args.out}")


if __name__ == "__main__":
    main()
