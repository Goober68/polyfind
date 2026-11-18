"""The polished cells of a row relaxed again with the cell free to shear along its chains.

examples/polish_margins.py relaxes each two-chain cell with every torsion and backbone angle free
and the cell free in ``a``, ``b``, ``gamma`` and ``c``: an upright cell, in which a chain and its
image one lattice vector away are at the same height.  A variable-cell DFT relaxation has two more
freedoms, the components of ``a`` and ``b`` along the chains, and the reversed-unit candidate's
all-trans cell used them (docs/MLIP_SECOND_OPINION.md 3.8).  This script gives polyfind's cell the
same two (``polyfind.mixedcell.MixedModel(shear=True)``) and relaxes each stored polished cell
again from where it is, at no shear: the minimum polyfind's own energy has nearest the upright
cell.  It reports, per cell, the energy gained and the shear taken, and per row the rigid margin
all-trans minus TG before and after.

It is a local relaxation, started a twentieth of an angstrom off the upright cell in both shears
(``--kick``): many of the stored cells have a symmetry that makes no shear a stationary point, a
minimum or not, and a relaxation started on it would stay.  The energy is periodic in each shear
with the chain's repeat and has other minima; ``--starts`` adds starts at other shears (fractions
of a monomer along ``a`` and ``b``), each relaxed first as rigid chains and then with everything
free, and the lowest is kept.

    python examples/shear_polish.py --rows vdf5-rvdf1-6
    python examples/shear_polish.py --minutes 100          # every polished row, resumable
"""
import argparse
import importlib.util
import json
import time
from pathlib import Path

import numpy as np

from polyfind import mixedcell as MC
from polyfind import pack as pack_mod
from polyfind.fitting import PRESET_PARAMETERS

DS = importlib.util.spec_from_file_location("dk", Path(__file__).with_name("dilute_kink.py"))
DK = importlib.util.module_from_spec(DS)
DS.loader.exec_module(DK)
C = DK.C
POLISHED = Path("docs/predictions/polished_margins_belt4.json")
REFERENCES = Path("docs/predictions/dilute_kink_belt4.json")


def stored_cell(row: str, conf: str, polished: dict, references: dict) -> dict:
    """The polished upright cell of ``row`` in conformation ``conf``, arrays as arrays."""
    stored = references["all_trans_references"][row] if conf == "TTTT" else polished["references"][f"{row} {conf}"]
    return {k: (np.asarray(v, dtype=float) if isinstance(v, list) else v) for k, v in stored.items()}


def sheared_model(polymer, d: dict):
    """``(model, x, chain)``: the polished two-chain cell ``d`` in a model whose cell may shear, at no shear."""
    ch = DK.reference_chain(polymer, d)
    pk = pack_mod.CrystalPacker(ch, n_chains=2, **DK.KW)
    st = MC.MixedState.from_vector(d["state"], 2, len(d["tors"]))
    st.shear = np.zeros(2)
    model = MC.MixedModel(polymer, "ref", pk, np.asarray(d["flips"]) > 0.5, [ch, ch], angle_stiffness=DK.angle_stiffness(), shear=True)
    x = st.vector()
    if abs(model.physical_energy(x) - float(d["energy"])) > 1e-6:
        raise RuntimeError("the cell that may shear does not reproduce the polished upright cell at no shear")
    return model, x, ch


def free_sets(model):
    """``(everything, rigid chains)``: the variables of the polish (one chain's position held, which
    removes the cell's translation) with the two shears, and the same without the chains' shapes."""
    everything = np.setdiff1d(np.arange(model.n_var), np.concatenate([model.index("xy", 1), model.index("z", 1)]))
    rigid = np.concatenate([model.index("lattice"), model.index("shear"), model.index("phi"), model.index("xy", 0), model.index("z", 0)])
    return everything, rigid


def relax_sheared(model, x0, starts=(), maxiter: int = 800, kick: float = 0.05):
    """The cell relaxed with the shear free from ``x0`` moved off the upright cell by ``kick`` A
    (``(1, 0.6) * kick`` in the two shears: off every symmetry that could hold it), and from each
    extra start ``(sa, sb)`` in A (rigid chains first there, then everything); returns
    ``(x, info, tried)`` of the lowest."""
    everything, rigid = free_sets(model)
    tried = []
    best = None
    for k, start in enumerate([None] + list(starts)):
        x = np.array(x0, dtype=float)
        if start is None:
            x[model.index("shear")] += np.array([1.0, 0.6]) * kick
        else:
            x[model.index("shear")] = start
            x, _ = model.relax(x, free=rigid, closed=[], maxiter=300)
        # the optimiser's own test can go on failing at a minimum of a large cell (MixedModel.relax): stop when the energy has
        x1, info = model.relax(x, free=everything, closed=[0, 1], maxiter=maxiter, patience=60, patience_tol=1e-6)
        e1 = model.physical_energy(x1)
        tried.append({"start_shear_A": [float(v) for v in x[model.index("shear")]] if start is None else [float(v) for v in start], "energy": float(e1),
                      "shear_A": [float(v) for v in model.state(x1).shear]})
        if best is None or e1 < best[0] - 1e-9:
            best = (e1, x1, info)
    return best[1], best[2], tried


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--rows", default="", help="comma list; default: every row of the polished table")
    ap.add_argument("--cells", default="", help="comma list of conformations; default: every polished cell of the row")
    ap.add_argument("--starts", default="", help="extra starting shears as fractions of a monomer, 'fa:fb,fa:fb'")
    ap.add_argument("--kick", type=float, default=0.05, help="the first start's shear off the upright cell, A")
    ap.add_argument("--tilts", action="store_true", help="relax the pendant tilts of each sheared cell as polish_margins.py does")
    ap.add_argument("--preset", default="pvdf-dft-valence-belt4")
    ap.add_argument("--maxiter", type=int, default=800)
    ap.add_argument("--minutes", type=float, default=0.0, help="start no new cell after this long; the file holds what is done")
    ap.add_argument("--out", default="docs/predictions/sheared_margins_belt4.json")
    args = ap.parse_args()
    out_path = Path(args.out)
    polished, references = json.load(open(POLISHED)), json.load(open(REFERENCES))
    rows = [r for r in args.rows.split(",") if r] or list(polished["rows"])
    fractions = [tuple(float(v) for v in s.split(":")) for s in args.starts.split(",") if s]
    log = lambda s: print(s, flush=True)  # noqa: E731
    began = time.time()
    disk = json.load(open(out_path)) if out_path.exists() else {}
    disk.update(preset=args.preset, energy="kcal/mol per monomer, rigid pendants unless a cell has tilt_per_monomer; margin = all-trans minus TG",
                what="each polished upright cell (polished_margins_belt4.json) relaxed again with the cell free to shear along its chains, from no shear"
                     + (f" and from the shears {args.starts} (fractions of a monomer)" if fractions else ""))
    disk.setdefault("rows", {})

    def save():
        with open(out_path, "w") as fh:
            json.dump(disk, fh, indent=1)

    with PRESET_PARAMETERS[args.preset].applied():
        for row in rows:
            letters, _ = DK.sequence(row)
            polymer = C.make_polymer(row, letters)
            entry = disk["rows"].setdefault(row, {"letters": letters, "cells": {}})
            confs = [c for c in args.cells.split(",") if c] or list(polished["rows"][row]["cells"])
            for conf in confs:
                if conf in entry["cells"]:
                    continue
                if args.minutes and time.time() - began > 60 * args.minutes:
                    save()
                    log(f"the {args.minutes:g} minutes are used; stopped before {row} {conf}")
                    return
                t0 = time.time()
                model, x0, ch = sheared_model(polymer, stored_cell(row, conf, polished, references))
                per = 2 * ch.n_monomers
                st0 = model.state(x0)
                cm = float(st0.lattice[3]) / ch.n_monomers
                e0, g0 = model.value_and_grad(x0)
                e0 = model.physical_energy(x0)
                x1, info, tried = relax_sheared(model, x0, [(fa * cm, fb * cm) for fa, fb in fractions], args.maxiter, args.kick)
                st1 = model.state(x1)
                e1 = model.physical_energy(x1)
                rec = {"upright_per_monomer": e0 / per, "sheared_per_monomer": e1 / per, "gain_per_monomer": (e1 - e0) / per,
                       "shear_A": [float(v) for v in st1.shear], "shear_monomers": [float(v) / cm for v in st1.shear],
                       "upright_shear_gradient_kcal_per_A_per_cell": [float(v) for v in g0[model.index("shear")]],
                       "cell_upright": st0.lattice.tolist(), "cell_sheared": st1.lattice.tolist(), "c_per_monomer": float(st1.lattice[3]) / ch.n_monomers,
                       "monomers_per_chain": int(ch.n_monomers), "evaluations": int(info["evaluations"]), "gradient_max": float(info["gradient_max"]),
                       "message": info["message"], "stalled": bool(info["stalled"]), "state": [float(v) for v in x1], "flips": [bool(f) for f in model.flips]}
                if fractions:
                    rec["starts"] = [{**t, "energy_per_monomer": t["energy"] / per} for t in tried]
                if args.tilts:
                    dE, tilts, _ = MC.tilt_relaxation(model, x1)
                    rec["tilt_per_monomer"] = float(dE) / per
                    rec["largest_tilt_deg"] = float(np.abs(tilts).max()) if np.size(tilts) else 0.0
                rec["seconds"] = time.time() - t0
                entry["cells"][conf] = rec
                save()
                log(f"{row:28s} {conf:6s}: {rec['upright_per_monomer']:.4f} -> {rec['sheared_per_monomer']:.4f} ({rec['gain_per_monomer']:+.4f}) kcal/mol per monomer; "
                    f"shear {np.round(rec['shear_monomers'], 3).tolist()} monomers; {rec['evaluations']} ev, {rec['seconds']:.0f} s; {info['message']}")
            cells = entry["cells"]
            tgs = [c for c in cells if c != "TTTT"]
            if "TTTT" in cells and tgs and all(c in cells for c in polished["rows"][row]["cells"]):
                up = min(tgs, key=lambda c: cells[c]["upright_per_monomer"])
                sh = min(tgs, key=lambda c: cells[c]["sheared_per_monomer"])
                entry["tg_conformer_upright"], entry["tg_conformer_sheared"] = up, sh
                entry["margin_rigid_upright"] = cells["TTTT"]["upright_per_monomer"] - cells[up]["upright_per_monomer"]
                entry["margin_rigid_sheared"] = cells["TTTT"]["sheared_per_monomer"] - cells[sh]["sheared_per_monomer"]
                entry["margin_change"] = entry["margin_rigid_sheared"] - entry["margin_rigid_upright"]
                save()
                log(f"  -> {row}: rigid margin {entry['margin_rigid_upright']:+.4f} upright, {entry['margin_rigid_sheared']:+.4f} with the shear free "
                    f"({entry['margin_change']:+.4f}); lowest TG {up} -> {sh}")
    log(f"written {out_path}")


if __name__ == "__main__":
    main()
