"""Starts for a cell whose chain axis tilts as it relaxes (the shear along the chains arriving through the third lattice vector).

A kink sheet started upright in the potential may never come to rest: its chain axis tilts a little further every leg of the
relaxation, by forces of 0.01 to 0.02 eV/A, and the energy falls a third of a kcal/mol per kink pair a leg
(docs/MLIP_SECOND_OPINION.md 3.20). Two starts settle what such a slide means:

    python examples/cell_tilt.py like --flat ALL_TRANS.start.xyz --like SHEET.start.xyz --frame NAME --out NAME.start.xyz
        the all-trans cell given the sheet's shear along the chains: the control. If it rests there and lower than the
        upright cell, the sheet was falling with its reference.

    python examples/cell_tilt.py push --prev LEG_A.start.xyz --now LEG_B.start.xyz --times 4 --frame NAME --out NAME.start.xyz
        LEG_B with its shear along the chains moved TIMES times further than it moved from LEG_A to LEG_B: a start beyond
        the slide's end, which comes back to it, so the end is found from both sides instead of walked.

The shear along the chains is how far the neighbours along the first two lattice vectors sit along the chain axis, in monomers.
It does not depend on how the cell lies in space, which the third vector's own components do (a relaxed cell is free to turn).
Only the chain axis is turned: the first two vectors, the repeat length and the fractional coordinates are kept, and the relaxer
(examples/mlip_train.py relax) does the rest.
"""
from __future__ import annotations

import argparse

import numpy as np

MONOMER_A = 2.56  # A of chain per monomer, for reading a shear in monomers


def shear_along_chains(cell: np.ndarray) -> np.ndarray:
    """How far the neighbours along the first two vectors sit along the chain axis (the third vector), in monomers."""
    axis = cell[2] / np.linalg.norm(cell[2])
    return np.array([cell[0] @ axis, cell[1] @ axis]) / MONOMER_A


def axis_for(cell: np.ndarray, shear: np.ndarray) -> np.ndarray:
    """The unit chain axis that gives the first two vectors this shear along the chains, on the third vector's side of their plane."""
    a, b = cell[0], cell[1]
    x, y = np.linalg.solve(np.array([[a @ a, a @ b], [a @ b, b @ b]]), np.asarray(shear, float) * MONOMER_A)
    in_plane = x * a + y * b
    if in_plane @ in_plane >= 1.0:
        raise ValueError(f"no chain axis gives a shear of {shear} monomers to this cell")
    normal = np.cross(a, b)
    normal *= np.sign(normal @ cell[2]) / np.linalg.norm(normal)
    return in_plane + np.sqrt(1.0 - in_plane @ in_plane) * normal


def with_shear(atoms, shear: np.ndarray, frame: str, how: str):
    """The cell with its chain axis turned to give this shear along the chains: repeat length and fractional coordinates kept."""
    cell = atoms.cell.array.copy()
    cell[2] = axis_for(cell, shear) * np.linalg.norm(cell[2])
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
    s = sub.add_parser("like", help="an all-trans cell given another cell's shear along the chains")
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
        shear = shear_along_chains(b.cell.array)
        new = with_shear(a, shear, args.frame, f"the chain axis turned to the shear along the chains of {b.info['frame']} ({shear[0]:+.2f}, {shear[1]:+.2f} monomers)")
    else:
        p, a = read(args.prev), read(args.now)
        now = shear_along_chains(a.cell.array)
        shear = now + args.times * (now - shear_along_chains(p.cell.array))
        new = with_shear(a, shear, args.frame, f"from {a.info['frame']}, the shear along the chains pushed {args.times:g} times its slide since "
                                               f"{p.info['frame']} further ({shear[0]:+.2f}, {shear[1]:+.2f} monomers)")
    write(args.out, new, format="extxyz")
    print(f"{args.frame}: shear along the chains {np.round(shear_along_chains(new.cell.array), 2)} monomers; cell {np.round(new.cell.cellpar(), 2)}")


if __name__ == "__main__":
    main()
