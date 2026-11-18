"""A candidate's cells under a crystal-stage fit: its margin with the cell held upright and with the shear free.

The acceptance test of a refitted potential on the reversed-unit candidate (docs/CRYSTAL_STAGE.md
5.25 to 5.27).  DFT's numbers for `vdf5-rvdf1-6` (docs/MLIP_SECOND_OPINION.md 3.8): all-trans
minus helix +0.070 kcal/mol per monomer with the all-trans cell held upright, -0.261 with
everything relaxed; the all-trans cell gains 0.35 by shearing to (0.30, -0.90) monomers and the
helix cell gains nothing.

Each polished cell of the row is taken as stored (examples/polish_margins.py, relaxed under
belt4), relaxed under the fit's parameters with the cell upright, then with the shear free from
several starts, and the pendant groups of each end cell are relaxed as rigid tilts.  The chains
are built with the fit's bond sites and each backbone bond takes its own type's torsion triple
(``--torsions``: the other three types, from the conformer refit).

    python examples/candidate_acceptance.py --fit docs/predictions/crystal_stage_v13.json
"""
import argparse
import importlib.util
import json
import time
from dataclasses import replace
from pathlib import Path

import numpy as np

from polyfind import crystalfit as C
from polyfind import fitting as F
from polyfind import mixedcell as MC

SPEC = importlib.util.spec_from_file_location("sp", Path(__file__).with_name("shear_polish.py"))
SP = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(SP)
DFT = {"vdf5-rvdf1-6": {"margin_upright": +0.070, "margin_relaxed": -0.261, "all_trans_shear_gain": -0.346, "all_trans_shear_monomers": (0.30, -0.90)}}


def parameters(fit: dict, torsions: dict | None) -> F.FFParameters:
    """The fit's parameter set, the three other torsion types from ``torsions`` (a valence vector)."""
    prm = C.parameters(np.array(fit["x"], dtype=float))
    if torsions is not None:
        block = np.array(torsions["x"], dtype=float)[:12].reshape(4, 3)
        prm = replace(prm, torsion_types=tuple((t, *(float(v) for v in block[k])) for k, t in enumerate(F.REF_TORSION_TYPES) if k))
    return prm


def relaxed(model, x, free):
    x1, info = model.relax(x, free=free, closed=[0, 1], maxiter=800, patience=60, patience_tol=1e-6)
    if info["stalled"] or info["gradient_max"] > 1e-3:
        x1, info = model.relax_on_surface(x1, free=free, closed=[0, 1], gtol=1e-4, maxiter=300, patience=60)
    return x1, info


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--fit", required=True, help="a crystal-stage result file")
    ap.add_argument("--torsions", default="", help="a valence vector file for the three other torsion types (docs/predictions/valence_torsions_*.json)")
    ap.add_argument("--row", default="vdf5-rvdf1-6")
    ap.add_argument("--confs", default="", help="comma list; default every polished cell of the row")
    ap.add_argument("--no-tilts", action="store_true")
    ap.add_argument("--out", default="")
    args = ap.parse_args()
    fit = json.loads(Path(args.fit).read_text())
    prm = parameters(fit, json.loads(Path(args.torsions).read_text()) if args.torsions else None)
    print(prm.describe())
    row = args.row
    pol, ref = json.load(open(SP.POLISHED)), json.load(open(SP.REFERENCES))
    confs = [c for c in args.confs.split(",") if c] or list(pol["rows"][row]["cells"])
    out = {"fit": args.fit, "torsions": args.torsions, "row": row, "cells": {}}
    SP.C.SITES_FOR = prm
    with prm.applied():
        letters, _ = SP.DK.sequence(row)
        polymer = SP.C.make_polymer(row, letters)
        for conf in confs:
            t0 = time.time()
            model, x0, ch = SP.sheared_model(polymer, SP.stored_cell(row, conf, pol, ref), check=False)
            per = 2 * ch.n_monomers
            cm = float(model.state(x0).lattice[3]) / ch.n_monomers
            everything, rigid = SP.free_sets(model)
            upright = np.setdiff1d(everything, model.index("shear"))
            xu, _ = relaxed(model, x0, upright)
            eu = model.physical_energy(xu)
            best = None
            line = DFT.get(row, {}).get("all_trans_shear_monomers", (0.3, -0.9))
            for start in (None, (line[0] * cm, line[1] * cm), (0.5 * line[0] * cm, 0.5 * line[1] * cm), (-line[0] * cm, -line[1] * cm), (0.0, 0.5 * cm), (0.5 * cm, 0.0)):
                x = xu.copy()
                if start is None:
                    x[model.index("shear")] += [0.05, 0.03]
                else:
                    x[model.index("shear")] = start
                    x, _ = model.relax(x, free=rigid, closed=[], maxiter=300)
                x1, _ = relaxed(model, x, everything)
                e1 = model.physical_energy(x1)
                print(f"   {conf} from {'the upright cell' if start is None else np.round(np.array(start) / cm, 3).tolist()}: {(e1 - eu) / per:+.4f} at shear "
                      f"{np.round(model.state(x1).shear / cm, 3).tolist()} monomers", flush=True)
                if best is None or e1 < best[0] - 1e-9:
                    best = (e1, x1)
            su, sb = model.state(xu), model.state(best[1])
            rec = {"upright": eu / per, "sheared": best[0] / per, "gain": (best[0] - eu) / per, "shear_monomers": (sb.shear / cm).tolist(),
                   "cell_upright": su.lattice.tolist(), "cell_sheared": sb.lattice.tolist()}
            if not args.no_tilts:
                rec["tilt_upright"] = float(MC.tilt_relaxation(model, xu)[0]) / per
                rec["tilt_sheared"] = float(MC.tilt_relaxation(model, best[1])[0]) / per
            out["cells"][conf] = rec
            print(f"{row} {conf}: upright {rec['upright']:.4f}, shear free {rec['sheared']:.4f} ({rec['gain']:+.4f} at {np.round(rec['shear_monomers'], 3).tolist()} monomers)"
                  + (f"; pendant tilts {rec['tilt_upright']:+.4f}, {rec['tilt_sheared']:+.4f}" if not args.no_tilts else "") + f"  ({time.time() - t0:.0f} s)", flush=True)
    cells = out["cells"]
    if "TTTT" in cells and len(cells) > 1:
        tt = cells["TTTT"]
        tilt = (lambda r, k: r.get(f"tilt_{k}", 0.0))
        tg_u = min((v for k, v in cells.items() if k != "TTTT"), key=lambda v: v["upright"] + tilt(v, "upright"))
        tg_s = min((v for k, v in cells.items() if k != "TTTT"), key=lambda v: v["sheared"] + tilt(v, "sheared"))
        out["margin_upright_rigid"] = tt["upright"] - min(v["upright"] for k, v in cells.items() if k != "TTTT")
        out["margin_sheared_rigid"] = tt["sheared"] - min(v["sheared"] for k, v in cells.items() if k != "TTTT")
        out["margin_upright"] = tt["upright"] + tilt(tt, "upright") - tg_u["upright"] - tilt(tg_u, "upright")
        out["margin_sheared"] = tt["sheared"] + tilt(tt, "sheared") - tg_s["sheared"] - tilt(tg_s, "sheared")
        out["all_trans_shear_gain"] = tt["sheared"] + tilt(tt, "sheared") - tt["upright"] - tilt(tt, "upright")
        d = DFT.get(row)
        print(f"all-trans minus helix: upright {out['margin_upright']:+.4f} (rigid pendants {out['margin_upright_rigid']:+.4f}), shear free {out['margin_sheared']:+.4f} "
              f"(rigid {out['margin_sheared_rigid']:+.4f}); the all-trans cell's gain {out['all_trans_shear_gain']:+.4f}"
              + (f"   DFT: {d['margin_upright']:+.3f}, {d['margin_relaxed']:+.3f}, {d['all_trans_shear_gain']:+.3f}" if d else ""))
        if d:
            out["dft"] = d
    if args.out:
        Path(args.out).write_text(json.dumps(out, indent=1))
        print(f"written {args.out}")


if __name__ == "__main__":
    main()
