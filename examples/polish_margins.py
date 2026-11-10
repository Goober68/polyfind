"""The switching margin of a sequence with both of its cells relaxed alike: every torsion free, closure exact.

examples/candidate_screen.py takes each conformation's cell from the packing search, which
refines a recognised state pattern (all-trans, TGTG') under the line group of that pattern: one
torsion value per symmetry class.  For PVDF that is the molecule's own symmetry.  For a chiral
(isotactic copolymer) sequence it is not: the all-trans chain is held planar, where the relaxed
chain is deflected by up to ten degrees per bond and lies lower by up to 0.2 kcal/mol per monomer,
which is the size of the margins the screen ranks by.

This script polishes both cells of a row in the mixed cell (examples/dilute_kink.py ``polish``:
every torsion and backbone angle of both chains free, the cell free, each chain periodic exactly),
relaxes the pendant tilts of each as the screen does, and reports the margin TT - TG per monomer
beside the screen's.  The TG conformer is the one the screen found lowest for the row.
"""
import argparse
import importlib.util
import json
import time
from pathlib import Path

import numpy as np

from polyfind import mixedcell as MC
from polyfind.fitting import PRESET_PARAMETERS

DS = importlib.util.spec_from_file_location("dk", Path(__file__).with_name("dilute_kink.py"))
DK = importlib.util.module_from_spec(DS)
DS.loader.exec_module(DK)
C = DK.C


def letters_of(t) -> str:
    """G+ / G- / T per torsion, run-length compressed."""
    names = ["G+" if 0.0 < v < 120.0 else ("G-" if -120.0 < v < 0.0 else "T") for v in DK.wrap(t)]
    return " ".join(names)


def screen_row(row: str, window: float):
    """``(tt, {name: tg row})`` of the screen's cells for ``row``, from whichever screen file has
    it: the all-trans cell and every TG conformer the screen ran that lies within ``window``
    kcal/mol per monomer of its lowest (the polish can reorder conformers that close)."""
    for f in sorted(Path("docs/predictions").glob("candidate_screen_belt4*.json")):
        seqs = json.load(open(f)).get("sequences", {})
        if row in seqs and seqs[row].get("derived"):
            confs = {k: v for k, v in seqs[row]["conformations"].items() if v and "energy_per_monomer_kcal" in v}
            tgs = {k: v for k, v in confs.items() if k != "TTTT"}
            low = min(v["energy_per_monomer_kcal"] for v in tgs.values())
            return confs["TTTT"], {k: v for k, v in tgs.items() if v["energy_per_monomer_kcal"] <= low + window}
    raise KeyError(f"{row} is in no screen file with a derived margin")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--rows", default="vdf9-trfe4-cfe1-14,vdf10-trfe4-14")
    ap.add_argument("--tg", default="", help="comma list of TG conformer names (candidate_screen.CONFORMATION_UNITS); default: the screen's conformers within --window of its lowest")
    ap.add_argument("--window", type=float, default=0.3, help="kcal/mol per monomer above the screen's lowest TG conformer within which a conformer is polished too")
    ap.add_argument("--preset", default="pvdf-dft-valence-belt4")
    ap.add_argument("--no-tilts", action="store_true")
    ap.add_argument("--references", default="docs/predictions/dilute_kink_belt4.json", help="the polished all-trans cells (shared with dilute_kink.py)")
    ap.add_argument("--out", default="docs/predictions/polished_margins_belt4.json")
    args = ap.parse_args()
    out_path, ref_path = Path(args.out), Path(args.references)
    log = lambda s: print(s, flush=True)  # noqa: E731

    def locked_update(path: Path, fn):
        import fcntl
        with open(str(path) + ".lock", "w") as lock:
            fcntl.flock(lock, fcntl.LOCK_EX)
            disk = json.load(open(path)) if path.exists() else {}
            fn(disk)
            with open(path, "w") as fh:
                json.dump(disk, fh, indent=1)

    with PRESET_PARAMETERS[args.preset].applied():
        for row in [r for r in args.rows.split(",") if r]:
            letters, note = DK.sequence(row)
            polymer = C.make_polymer(row, letters)
            s_tt, s_tgs = screen_row(row, args.window)
            tg_names = [t for t in args.tg.split(",") if t] or sorted(s_tgs, key=lambda k: s_tgs[k]["energy_per_monomer_kcal"])
            P = polymer.bonds_per_repeat
            s_low = min(s_tgs.values(), key=lambda v: v["energy_per_monomer_kcal"])
            prior = (json.load(open(out_path)) if out_path.exists() else {}).get("rows", {}).get(row, {})
            entry = {"letters": letters, "note": note, "tg_conformers": tg_names,
                     "screen": {"margin": s_tt["energy_per_monomer_kcal"] - s_low["energy_per_monomer_kcal"],
                                "margin_rigid": s_tt["energy_rigid_per_monomer_kcal"] - s_low["energy_rigid_per_monomer_kcal"]},
                     "cells": dict(prior.get("cells", {}))}
            cells = entry["cells"]
            for conf in ["TTTT"] + tg_names:
                if conf in cells and "energy_per_monomer" in cells[conf]:
                    continue
                t0 = time.time()
                if conf == "TTTT":
                    stored = (json.load(open(ref_path)) if ref_path.exists() else {}).get("all_trans_references", {}).get(row)
                    d = DK.all_trans_reference(polymer, row, stored)
                    if stored is None or "state" not in stored:
                        locked_update(ref_path, lambda disk, d=d: disk.setdefault("all_trans_references", {}).__setitem__(row, DK.as_json(d)))
                else:
                    stored = (json.load(open(out_path)) if out_path.exists() else {}).get("references", {}).get(f"{row} {conf}")
                    d = DK.refined_cell(polymer, row, C.conformation(conf, P), f"{row} {conf}", stored)
                    if stored is None or "state" not in stored:
                        locked_update(out_path, lambda disk, d=d, conf=conf: disk.setdefault("references", {}).__setitem__(f"{row} {conf}", DK.as_json(d)))
                model, x, ch = DK.make_model(polymer, d, 1, 1)
                n_mon = ch.n_monomers
                st = model.state(x)
                ideal = np.array([[180.0, 60.0, -60.0][v] for v in ([0] * P if conf == "TTTT" else C.conformation(conf, P))])
                rec = {"search_rigid_per_monomer": float(d["search_energy"]) / (2 * n_mon), "rigid_per_monomer": float(d["energy"]) / (2 * n_mon),
                       "polish_gain_per_monomer": (float(d["search_energy"]) - float(d["energy"])) / (2 * n_mon),
                       "search_closure_residual_deg": float(d["search_closure_residual"][0]), "cell": st.lattice.tolist(), "c_per_monomer": float(st.lattice[3]) / n_mon,
                       "torsions": [letters_of(st.tors[s]) for s in range(2)],
                       "largest_torsion_off_ideal_deg": float(np.abs(DK.wrap(st.tors - ideal[None])).max())}
                if not args.no_tilts:
                    t1 = time.time()
                    dE, tilts, _ = MC.tilt_relaxation(model, x)
                    rec["tilt_per_monomer"] = float(dE) / (2 * n_mon)
                    rec["largest_tilt_deg"] = float(np.abs(tilts).max()) if np.size(tilts) else 0.0
                    rec["tilt_seconds"] = time.time() - t1
                else:
                    rec["tilt_per_monomer"] = 0.0
                rec["energy_per_monomer"] = rec["rigid_per_monomer"] + rec["tilt_per_monomer"]
                rec["seconds"] = time.time() - t0
                cells[conf] = rec
                log(f"{row:24s} {conf:6s}: rigid {rec['rigid_per_monomer']:.4f} ({rec['polish_gain_per_monomer']:.4f} below the search's {rec['search_rigid_per_monomer']:.4f}), "
                    f"tilts {rec['tilt_per_monomer']:+.4f}: {rec['energy_per_monomer']:.4f} kcal/mol per monomer; cell {np.round(st.lattice, 3).tolist()}; "
                    f"torsions up to {rec['largest_torsion_off_ideal_deg']:.1f} deg off the ideal pattern  ({time.time() - t0:.0f} s)")
                locked_update(out_path, lambda disk: disk.setdefault("rows", {}).__setitem__(row, entry))
            tt = cells["TTTT"]
            tg_name = min((k for k in cells if k != "TTTT"), key=lambda k: cells[k]["energy_per_monomer"])
            tg = cells[tg_name]
            entry["tg_conformer"] = tg_name
            entry["margin"] = tt["energy_per_monomer"] - tg["energy_per_monomer"]
            entry["margin_rigid"] = tt["rigid_per_monomer"] - tg["rigid_per_monomer"]
            entry["axial_transition_strain"] = tt["c_per_monomer"] / tg["c_per_monomer"] - 1.0
            log(f"  -> {row}: margin TT-TG {entry['margin']:+.3f} kcal/mol per monomer polished (rigid {entry['margin_rigid']:+.3f}; lowest TG {tg_name}); "
                f"the screen's {entry['screen']['margin']:+.3f} (rigid {entry['screen']['margin_rigid']:+.3f})")
            locked_update(out_path, lambda disk: (disk.update(preset=args.preset, energy="kcal/mol per monomer; margin = all-trans minus TG, positive when TG is lower"),
                                                  disk.setdefault("rows", {}).__setitem__(row, entry)))
    log(f"written {out_path}")


if __name__ == "__main__":
    main()
