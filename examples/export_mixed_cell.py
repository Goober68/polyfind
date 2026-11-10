"""Write a polished two-chain cell (or its tiling) as an extended-XYZ file a DFT code can read.

The cells of examples/dilute_kink.py, kink_migration.py and polish_margins.py are stored as the
variable vector of a mixed cell.  This script rebuilds the atoms (real atoms only: the repulsion
belt sites of the fitted potential are not atoms) and writes them with the lattice vectors in the
comment line (``Lattice="ax ay az bx by bz cx cy cz"``), c along z.

    python examples/export_mixed_cell.py --row vdf7-trfe3 --cell TTTT --out docs/predictions/polished_vdf7-trfe3_TTTT.xyz
    python examples/export_mixed_cell.py --row vdf10 --cell kink@0 --out docs/predictions/polished_vdf10_kink0.xyz
"""
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np

from polyfind.fitting import PRESET_PARAMETERS

DS = importlib.util.spec_from_file_location("dk", Path(__file__).with_name("dilute_kink.py"))
DK = importlib.util.module_from_spec(DS)
DS.loader.exec_module(DK)
ATOMS = ("C", "H", "F", "Cl", "N", "O")


def stored_cell(row: str, cell: str) -> dict:
    """The stored polished cell: ``TTTT`` (dilute_kink_belt4.json), a TG conformer name
    (polished_margins_belt4.json) or ``kink@<offset>`` (kink_migration_belt4.json)."""
    base = Path("docs/predictions")
    if cell == "TTTT":
        return json.load(open(base / "dilute_kink_belt4.json"))["all_trans_references"][row]
    if cell.startswith("kink@"):
        return json.load(open(base / "kink_migration_belt4.json"))["kinked_references"][f"{row}@{cell[5:]}"]
    return json.load(open(base / "polished_margins_belt4.json"))["references"][f"{row} {cell}"]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--row", required=True)
    ap.add_argument("--cell", default="TTTT")
    ap.add_argument("--nx", type=int, default=1)
    ap.add_argument("--ny", type=int, default=1)
    ap.add_argument("--preset", default="pvdf-dft-valence-belt4")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    letters, note = DK.sequence(args.row)
    with PRESET_PARAMETERS[args.preset].applied():
        polymer = DK.C.make_polymer(args.row, letters)
        d = {k: (np.asarray(v, dtype=float) if isinstance(v, list) else v) for k, v in stored_cell(args.row, args.cell).items()}
        model, x, ch = DK.make_model(polymer, d, args.nx, args.ny)
        latn, _, _, Pn, _ = model.build(x)
        el = np.array([str(e) for e in list(model.pk.elements) * model.nc])
        keep = np.array([e in ATOMS for e in el])
        e_cell = model.physical_energy(x)
    P, el = Pn[keep], el[keep]
    frac = P @ np.linalg.inv(latn)
    P = (frac - np.floor(frac)) @ latn  # wrapped into the cell
    lat = " ".join(f"{v:.8f}" for v in latn.ravel())
    lines = [str(len(P)), f'Lattice="{lat}" Properties=species:S:1:pos:R:3 pbc="T T T" source="polyfind {args.preset} {args.row} {args.cell} '
                          f'{args.nx}x{args.ny}, {model.nc} chains of {ch.n_monomers} monomers; model energy {e_cell:.4f} kcal/mol per cell, rigid pendants"']
    lines += [f"{e:2s} {p[0]:14.8f} {p[1]:14.8f} {p[2]:14.8f}" for e, p in zip(el, P)]
    Path(args.out).write_text("\n".join(lines) + "\n")
    counts = {str(e): int((el == e).sum()) for e in sorted(set(el.tolist()))}
    print(f"{args.out}: {len(P)} atoms {counts}; cell {np.round(model.state(x).lattice, 4).tolist()}")


if __name__ == "__main__":
    main()
