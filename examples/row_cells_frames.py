"""The cells a preset gives a row, written as atoms for DFT.

examples/rows_under_preset.py stores each relaxed cell's state.  This writes those cells as
``NAME.start.xyz`` for ``examples/mlip_frames.py label``: the first step of testing a row's
margin against DFT (the cells as polyfind has them; then relaxed with the potential tilted to
DFT's forces, ``examples/mlip_train.py relax --dft-gradient``, and labelled again).  The pendant
groups are as the chain builder puts them: the rigid tilts polyfind relaxes them by are an
energy in the stored record, not a geometry, and DFT's relaxation moves them itself.

    python examples/row_cells_frames.py --margins docs/predictions/margins_chain1.json --row vdf3-trfe1-4 --kmesh "4 2 2" --out DIR
"""
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np

from polyfind.fitting import PRESET_PARAMETERS

SPEC = importlib.util.spec_from_file_location("shear_polish", Path(__file__).with_name("shear_polish.py"))
SP = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(SP)
ATOMS = ("C", "H", "F", "Cl", "N", "O")  # the sites that are atoms; the others are the potential's own (belt sites)


def safe(conf: str) -> str:
    """A conformation's name as part of a file name: TG'TG -> TGpTG."""
    return conf.replace("'", "p")


def atoms_of(model, x):
    """The cell of ``model`` at ``x`` as an ase.Atoms of its real atoms, chain after chain in the builder's order."""
    from ase import Atoms

    latn, _, _, Pn, _ = model.build(np.asarray(x, dtype=float))
    el = np.array([str(e) for e in list(model.pk.elements) * model.nc])
    keep = np.array([e in ATOMS for e in el])
    a = Atoms(symbols=list(el[keep]), positions=np.asarray(Pn, dtype=float)[keep], cell=np.asarray(latn, dtype=float), pbc=True)
    a.wrap()
    return a


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--margins", required=True, help="a file of examples/rows_under_preset.py")
    ap.add_argument("--row", required=True)
    ap.add_argument("--confs", default="", help="comma list; default every cell of the row")
    ap.add_argument("--which", default="sheared", choices=("sheared", "upright"), help="the cell free to shear along its chains, or held upright")
    ap.add_argument("--kmesh", required=True, help='the k-point mesh of the DFT, e.g. "4 2 2"')
    ap.add_argument("--tag", default="", help="a word for the frames' names; default the preset's last part")
    ap.add_argument("--out", required=True, help="working directory for the DFT, outside the repository")
    args = ap.parse_args()
    from ase.io import write

    disk = json.loads(Path(args.margins).read_text())
    entry = disk["rows"][args.row]
    preset = disk["preset"]
    tag = args.tag or preset.split("-")[-1]
    pol, ref = json.load(open(disk.get("polished", SP.POLISHED))), json.load(open(disk.get("references", SP.REFERENCES)))
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    prm = PRESET_PARAMETERS[preset]
    SP.C.SITES_FOR = prm
    with prm.applied():
        polymer = SP.C.make_polymer(args.row, entry["letters"])
        for conf in [c for c in args.confs.split(",") if c] or list(entry["cells"]):
            rec = entry["cells"][conf]
            model, _, ch = SP.sheared_model(polymer, SP.stored_cell(args.row, conf, pol, ref), check=False)
            x = np.asarray(rec[f"state_{args.which}"], dtype=float)
            a = atoms_of(model, x)
            energy = rec[args.which] + rec.get(f"tilt_{args.which}", 0.0)
            name = f"row_{args.row}_{safe(conf)}_{tag}"
            a.info.update(family="rowcell", frame=name, n_monomers=model.nc * ch.n_monomers, kmesh=args.kmesh, conformation=conf,
                          source=f"polyfind {preset}, {args.row} {conf}, cell {args.which}, pendants as built",
                          polyfind_kcal_per_monomer=float(energy))
            write(out / f"{name}.start.xyz", a, format="extxyz")
            print(f"{name}: {len(a)} atoms, cell {np.round(a.cell.cellpar(), 3).tolist()}, polyfind {energy:+.4f} kcal/mol per monomer")
    print(f"written {out}; label with examples/mlip_frames.py label --frames {out}")


if __name__ == "__main__":
    main()
