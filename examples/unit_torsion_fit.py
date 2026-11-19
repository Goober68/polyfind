"""A torsion triple for the bonds one unit makes, fitted to DFT's lone chains that have them, on top of a preset.

The lattice's torsion term is one triple for every backbone bond unless a bond type has its own
(``FFParameters.torsion_types``).  A unit whose DFT check puts its error inside the chain
(examples/unit_dft_check.py: a lone chain cut from the relaxed cell) gets a triple for its own
bonds here, as the reversed unit's CF2-CF2 and CH2-CH2 bonds did in the crystal stage
(docs/CRYSTAL_STAGE.md 5.29).  The rows are DFT's lone-chain conformers
(examples/chain_conformers.py collect), polyfind held at DFT's torsions with its angles, its
repeat and its pendants free; each conformer's energy is against its chain's all-trans one.  A
ridge holds the triple to the preset's own where the rows leave it open.

    python examples/unit_torsion_fit.py --base pvdf-dft-chain1 --type CH2-CFCl --dft docs/predictions/chain_conformers_cfe_pbe_d3bj.json \\
        --out docs/predictions/unit_torsion_cfe_chain1.json

The output names its base and the triples to put on it: every tool that takes a preset takes it
(``polyfind.fitting.parameters_of``).
"""
import argparse
import importlib.util
import json
import time
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

from polyfind import fitting as F

HERE = Path(__file__).resolve().parent
TOL = (0.05, 0.05)   # kcal/mol per monomer, and the part of a row's own DFT energy: a row's tolerance is their sum in quadrature
SIGMA = 2.0          # kcal/mol: the ridge on each coefficient's distance from the preset's own triple


def _load(name):
    spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--base", required=True, help="a preset's name, or a file polyfind.fitting.parameters_of reads")
    ap.add_argument("--type", required=True, action="append", help=f"a bond type of {F.LATTICE_TORSION_TYPES}; give it again for a second type fitted with it")
    ap.add_argument("--dft", type=Path, required=True, nargs="+", help="lone-chain conformers (examples/chain_conformers.py collect)")
    ap.add_argument("--only", default="", help="keep only conformers whose name contains one of these (comma list)")
    ap.add_argument("--leave", default="", help="conformers evaluated and shown but not fitted (comma list of name parts)")
    ap.add_argument("--maxfev", type=int, default=40)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    CC, SP = _load("chain_conformers"), _load("shear_polish")
    base = F.parameters_of(args.base)
    rows = {}
    for f in args.dft:
        rows.update(json.loads(f.read_text())["conformers"])
    only = [s for s in args.only.split(",") if s]
    if only:  # a chain's zero stays with it
        zeros = {r.get("zero", CC.REFERENCE) for k, r in rows.items() if any(s in k for s in only)}
        rows = {k: r for k, r in rows.items() if any(s in k for s in only) or k in zeros}
    leave = [s for s in args.leave.split(",") if s]
    fitted = [k for k, r in rows.items() if k != r.get("zero", CC.REFERENCE) and not any(s in k for s in leave)]
    plain = np.array(base.torsion, dtype=float)
    own = {t: np.array(v) for t, *v in base.torsion_types}
    x0 = np.concatenate([own.get(t, plain) for t in args.type])
    history = []

    def evaluate(x):
        prm = F.typed_on(base, {t: x[3 * k:3 * k + 3] for k, t in enumerate(args.type)})
        previous = SP.C.SITES_FOR
        SP.C.SITES_FOR = prm
        try:
            with prm.applied():
                mine = CC.polyfind_conformers(SP, SP.C, rows)
        finally:
            SP.C.SITES_FOR = previous
        table = {}
        for k, r in rows.items():
            zero = r.get("zero", CC.REFERENCE)
            if k != zero:
                e, dft = mine[k]["energy"] - mine[zero]["energy"], r["energy"]
                table[k] = (e, dft, (e - dft) / np.hypot(TOL[0], TOL[1] * dft))
        ridge = (np.asarray(x) - x0) / SIGMA
        z = np.concatenate([[table[k][2] for k in fitted], ridge])
        history.append((float(z @ z), np.array(x, dtype=float), table, float(ridge @ ridge)))
        return z

    def show(table):
        return "\n".join(f"   {k:16s} polyfind {v[0]:+8.3f}  DFT {v[1]:+8.3f}  off {v[0] - v[1]:+7.3f}{'' if k in fitted else '   (not fitted)'}" for k, v in table.items())

    t0 = time.time()
    evaluate(x0)
    print(f"{args.base}, {', '.join(args.type)} at ({', '.join(f'{v:+.4f}' for v in x0)}): objective {history[0][0]:.2f}")
    print(show(history[0][2]), flush=True)
    sol = least_squares(evaluate, x0, diff_step=0.02, x_scale=1.0, max_nfev=args.maxfev) if args.maxfev > 0 else None
    best = min(history, key=lambda h: h[0])
    print(f"fitted ({', '.join(f'{v:+.4f}' for v in best[1])}): objective {best[0]:.2f} (ridge {best[3]:.2f}) after {len(history)} evaluations, {time.time() - t0:.0f} s")
    print(show(best[2]))
    out = {"base": args.base, "torsion_types": {t: [float(v) for v in best[1][3 * k:3 * k + 3]] for k, t in enumerate(args.type)},
           "what": "torsion triples (V1, V2, V3; kcal/mol) to put on the base, fitted to DFT's lone-chain conformers (examples/unit_torsion_fit.py)",
           "dft": [str(f) for f in args.dft], "fitted": fitted, "objective": best[0], "objective_at_the_base": history[0][0], "ridge": best[3],
           "evaluations": len(history), "stopped": None if sol is None else str(sol.message),
           "rows": {k: {"polyfind": float(v[0]), "dft": float(v[1]), "at_the_base": float(history[0][2][k][0]), "fitted": k in fitted} for k, v in best[2].items()}}
    args.out.write_text(json.dumps(out, indent=1))
    print(f"written {args.out}")


if __name__ == "__main__":
    main()
