"""Starts for a cell whose chain axis tilts as it relaxes (the shear along the chains arriving through the third lattice vector).

A kink sheet started upright in the potential may never come to rest: its chain axis tilts a little further every leg of the
relaxation, by forces of 0.01 to 0.02 eV/A, and the energy falls a third of a kcal/mol per kink pair a leg
(docs/MLIP_SECOND_OPINION.md 3.20). Two starts settle what such a slide means:

    python examples/cell_tilt.py like --flat ALL_TRANS.start.xyz --like SHEET.start.xyz --frame NAME --out NAME.start.xyz
        the all-trans cell given the sheet's tilt (the third vector's components along the first two, in their units): the
        control. If it rests tilted and lower than the upright cell, the sheet was falling with its reference.

    python examples/cell_tilt.py push --prev LEG_A.start.xyz --now LEG_B.start.xyz --times 4 --frame NAME --out NAME.start.xyz
        LEG_B with its chain axis moved TIMES times further than it moved from LEG_A to LEG_B: a start beyond the slide's
        end, which comes back to it, so the end is found from both sides instead of walked.

Fractional coordinates and the repeat length along the chain are kept; the relaxer (examples/mlip_train.py relax) does the rest.
"""
from __future__ import annotations

import argparse

import numpy as np

MONOMER_A = 2.56  # A of chain per monomer, for reading a shear in monomers


def tilt(cell: np.ndarray) -> np.ndarray:
    """The third lattice vector's components across the chains, in units of the first two vectors."""
    return np.linalg.solve(np.array([cell[0][:2], cell[1][:2]]).T, cell[2][:2])


def shear_along_chains(cell: np.ndarray) -> np.ndarray:
    """How far the neighbours along the first two vectors sit along the chain axis, in monomers."""
    axis = cell[2] / np.linalg.norm(cell[2])
    return np.array([cell[0] @ axis, cell[1] @ axis]) / MONOMER_A


def with_axis(atoms, axis: np.ndarray, frame: str, how: str):
    """The cell with its third vector turned to ``axis`` at the same length, fractional coordinates kept."""
    cell = atoms.cell.array.copy()
    cell[2] = axis * (np.linalg.norm(cell[2]) / np.linalg.norm(axis))
    out = atoms.copy()
    out.set_cell(cell, scale_atoms=True)
    out.info["frame"] = frame
    out.info["source"] = f"{atoms.info.get('source', '')}; {how}"
    for key in ("note", "potential_eV_per_monomer"):  # the relaxed energy belongs to the cell this was made from
        out.info.pop(key, None)
    return out


def main():
    from ase.io import read, write

    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("like", help="an all-trans cell given another cell's tilt of the chain axis")
    s.add_argument("--flat", required=True)
    s.add_argument("--like", required=True)
    s = sub.add_parser("push", help="a sliding cell pushed on along its slide")
    s.add_argument("--prev", required=True)
    s.add_argument("--now", required=True)
    s.add_argument("--times", type=float, required=True)
    for s in sub.choices.values():
        s.add_argument("--frame", required=True)
        s.add_argument("--out", required=True)
    args = ap.parse_args()
    if args.cmd == "like":
        a, b = read(args.flat), read(args.like)
        d = tilt(b.cell.array) - tilt(a.cell.array)
        cell = a.cell.array
        new = with_axis(a, cell[2] + d[0] * cell[0] + d[1] * cell[1], args.frame,
                        f"the chain axis tilted {d[0]:+.3f} a, {d[1]:+.3f} b to match {b.info['frame']}")
    else:
        p, a = read(args.prev), read(args.now)
        new = with_axis(a, a.cell.array[2] + args.times * (a.cell.array[2] - p.cell.array[2]), args.frame,
                        f"from {a.info['frame']}, the chain axis pushed {args.times:g} times its slide since {p.info['frame']} further")
    write(args.out, new, format="extxyz")
    c = new.cell.array
    print(f"{args.frame}: tilt {np.round(tilt(c), 3)} of a, b; shear along the chains {np.round(shear_along_chains(c), 2)} monomers; "
          f"cell {np.round(new.cell.cellpar(), 2)}")


if __name__ == "__main__":
    main()
