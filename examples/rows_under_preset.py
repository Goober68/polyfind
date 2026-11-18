"""Every row's cells under a preset other than the one they were polished with: the margin, cell upright and with the shear free.

The polished cells of the 35 rows (examples/polish_margins.py) are belt4's.  Under another
preset each is relaxed again from there with the cell upright, then with the cell free to shear
along its chains from several starts, and the pendant groups of both are relaxed as rigid tilts
(as examples/candidate_acceptance.py does for one row).  The margin is all-trans minus the
lowest helix placement, kcal/mol per monomer, pendants tilted.

This is the first of four numbers a row is ranked on (CLAUDE.md): the cells' lengths along the
chain, which are the strain of the switch, are stored beside it; the loss, the rate and the
blocked stress are not here.

    python examples/rows_under_preset.py --preset pvdf-dft-chain1 --out docs/predictions/margins_chain1.json [--rows ROW,...] [--minutes 100]

The output file is also the checkpoint: a run goes on from the cells it holds.
"""
import argparse
import importlib.util
import json
import time
from pathlib import Path

import numpy as np

from polyfind import mixedcell as MC
from polyfind.fitting import PRESET_PARAMETERS

SPEC = importlib.util.spec_from_file_location("candidate_acceptance", Path(__file__).with_name("candidate_acceptance.py"))
CA = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CA)
SP = CA.SP
STARTS = ((0.3, -0.9), (0.15, -0.45), (0.0, 0.5))  # shears to start from besides a small kick off the upright cell, in monomers


def cell_under(model, x0, ch, starts=STARTS, tilts: bool = True) -> dict:
    """One stored cell relaxed under the parameter set in force: upright, then with the shear free
    from a kick and from each of ``starts`` (the lowest kept).  Energies in kcal/mol per monomer."""
    per = 2 * ch.n_monomers
    cm = float(model.state(x0).lattice[3]) / ch.n_monomers
    everything, rigid = SP.free_sets(model)
    upright = np.setdiff1d(everything, model.index("shear"))
    xu, iu = CA.relaxed(model, x0, upright)
    eu = model.physical_energy(xu)
    best = None
    for start in (None, *starts):
        x = xu.copy()
        if start is None:
            x[model.index("shear")] += [0.05, 0.03]
        else:
            x[model.index("shear")] = [start[0] * cm, start[1] * cm]
            x, _ = model.relax(x, free=rigid, closed=[], maxiter=300)
        x1, i1 = CA.relaxed(model, x, everything)
        e1 = model.physical_energy(x1)
        if best is None or e1 < best[0] - 1e-9:
            best = (e1, x1, i1)
    su, sb = model.state(xu), model.state(best[1])
    rec = {"upright": eu / per, "sheared": best[0] / per, "gain": (best[0] - eu) / per, "shear_monomers": [float(v) / cm for v in sb.shear],
           "cell_upright": su.lattice.tolist(), "cell_sheared": sb.lattice.tolist(), "c_per_monomer_upright": float(su.lattice[3]) / ch.n_monomers,
           "c_per_monomer_sheared": float(sb.lattice[3]) / ch.n_monomers, "monomers_per_chain": int(ch.n_monomers),
           "gradient_max": [float(iu["gradient_max"]), float(best[2]["gradient_max"])], "state_upright": [float(v) for v in xu],
           "state_sheared": [float(v) for v in best[1]], "flips": [bool(f) for f in model.flips]}
    if tilts:
        rec["tilt_upright"] = float(MC.tilt_relaxation(model, xu)[0]) / per
        rec["tilt_sheared"] = float(MC.tilt_relaxation(model, best[1])[0]) / per
    return rec


def margins(cells: dict) -> dict:
    """All-trans minus the lowest helix placement, upright and with the shear free, pendants tilted where they were."""
    e = lambda r, k: r[k] + r.get(f"tilt_{k}", 0.0)  # noqa: E731
    helix = [c for c in cells if c != "TTTT"]
    up, sh = min(helix, key=lambda c: e(cells[c], "upright")), min(helix, key=lambda c: e(cells[c], "sheared"))
    tt = cells["TTTT"]
    return {"helix_upright": up, "helix_sheared": sh, "margin_upright": e(tt, "upright") - e(cells[up], "upright"), "margin": e(tt, "sheared") - e(cells[sh], "sheared"),
            "all_trans_shear_gain": e(tt, "sheared") - e(tt, "upright"), "helix_shear_gain": e(cells[sh], "sheared") - e(cells[sh], "upright"),
            # the switch's strain along the chains: the all-trans cell's length a monomer against the helix cell's
            "strain_along_chain": tt["c_per_monomer_sheared"] / cells[sh]["c_per_monomer_sheared"] - 1.0}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--preset", required=True)
    ap.add_argument("--rows", default="", help="comma list; default every row of the polished table")
    ap.add_argument("--minutes", type=float, default=0.0, help="start no new cell after this long")
    ap.add_argument("--no-tilts", action="store_true")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    out_path = Path(args.out)
    pol, ref = json.load(open(SP.POLISHED)), json.load(open(SP.REFERENCES))
    rows = [r for r in args.rows.split(",") if r] or list(pol["rows"])
    disk = json.loads(out_path.read_text()) if out_path.exists() else {}
    disk.update(preset=args.preset, energy="kcal/mol per monomer; margin = all-trans minus the lowest helix placement, shear free, pendants tilted",
                what="each polished cell (polished_margins_belt4.json) relaxed again under the preset: upright, then free to shear along its chains from "
                     f"a kick and from {list(STARTS)} monomers; pendant tilts of both (examples/rows_under_preset.py)")
    disk.setdefault("rows", {})
    began = time.time()
    prm = PRESET_PARAMETERS[args.preset]
    SP.C.SITES_FOR = prm
    with prm.applied():
        for row in rows:
            letters, _ = SP.DK.sequence(row)
            polymer = SP.C.make_polymer(row, letters)
            entry = disk["rows"].setdefault(row, {"letters": letters, "cells": {}})
            for conf in pol["rows"][row]["cells"]:
                if conf in entry["cells"]:
                    continue
                if args.minutes and time.time() - began > 60 * args.minutes:
                    out_path.write_text(json.dumps(disk, indent=1))
                    print(f"the {args.minutes:g} minutes are used; stopped before {row} {conf}", flush=True)
                    return
                t0 = time.time()
                model, x0, ch = SP.sheared_model(polymer, SP.stored_cell(row, conf, pol, ref), check=False)
                rec = cell_under(model, x0, ch, tilts=not args.no_tilts)
                rec["seconds"] = time.time() - t0
                entry["cells"][conf] = rec
                out_path.write_text(json.dumps(disk, indent=1))
                print(f"{row:28s} {conf:6s}: upright {rec['upright'] + rec.get('tilt_upright', 0.0):.4f}, shear free {rec['sheared'] + rec.get('tilt_sheared', 0.0):.4f} "
                      f"at {np.round(rec['shear_monomers'], 3).tolist()} monomers; {rec['seconds']:.0f} s", flush=True)
            if "TTTT" in entry["cells"] and all(c in entry["cells"] for c in pol["rows"][row]["cells"]) and len(entry["cells"]) > 1:
                entry.update(margins(entry["cells"]))
                out_path.write_text(json.dumps(disk, indent=1))
                print(f"  -> {row}: margin {entry['margin']:+.4f} (upright {entry['margin_upright']:+.4f}); helix {entry['helix_sheared']}; "
                      f"strain along the chains {100 * entry['strain_along_chain']:+.1f}%", flush=True)
    print(f"written {out_path}")


if __name__ == "__main__":
    main()
