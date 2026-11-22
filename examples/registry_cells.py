"""Other registries and settings of a two-chain periodic cell (docs/MLIP_SECOND_OPINION.md 3.24).

A regular cell puts every chain's neighbours at one height.  3.23 found, by accident of a scan, that the film composition's
all-trans cell with its neighbours a monomer further along the chains is 0.84 kcal/mol per chain below the upright one in
DFT.  This makes such cells on purpose, from a labelled cell:

  sheared   the first two cell vectors given SA and SB monomers along the third (the chains' axis), the atoms left where
            they are: every chain keeps its place and its neighbours along those edges sit SA and SB monomers along
  turned    every chain turned about its own axis by an angle, the cell kept

    python -I examples/registry_cells.py build --labels FILES --cell LABEL:tt=FRAME LABEL:hx=FRAME ... \\
        --variant LABEL:tt:SA,SB,TURN ... [--like UPRIGHT LOW --variant LABEL:tt:like ...] --out STARTS.extxyz
    python -I examples/registry_cells.py table --starts STARTS.extxyz --relaxed DIR --labels FILES --out TABLE.json \\
        --frames-out RELAXED.extxyz
    python -I examples/registry_cells.py choose --table TABLE.json [--within 0.3] [--below LABEL:tt=KCAL ...] \
        [--unsettled-below KCAL] [--skip FRAME ...]
    python -I examples/registry_cells.py unsettled --table TABLE.json --below KCAL
    python -I examples/registry_cells.py margins --labels FILES --cell LABEL:tt=FRAME LABEL:hx=FRAME ... \
        [--also LABEL:tt=FRAME ...] --row LABEL=ROW ... --out MARGINS.json

`build` writes the starts in the order given (the order they are relaxed in).  `table` reads each start's relaxed cell
(examples/mlip_train.py relax: NAME_in_POTENTIAL.start.xyz) and gives its energy against the labelled cell it was made from,
kcal/mol per chain, with what it came to rest as.  `choose` prints two lines of names, the card's and the processors': the
lowest cell of each composition and conformation that came to rest, kept its conformation and is within --within of its
labelled cell (or below --below), the lowest first, dealt to the two lines in turn.  `margins` takes every labelled cell
of each composition (the two given, those of --also, and every labelled frame this tool made from them) and gives the
margin, all-trans less helix per monomer, between the lowest all-trans cell and the lowest helix cell in DFT, beside the
margin between the two cells given: the upright cells' margin.
"""

import argparse
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

KCAL = 23.0605
CHAINS = 2
BONDS = {("C", "C"): 1.75, ("C", "H"): 1.3, ("C", "F"): 1.6, ("C", "Cl"): 2.0, ("C", "N"): 1.4}


def _path_tool():
    spec = importlib.util.spec_from_file_location("mlip_path", Path(__file__).resolve().parent / "mlip_path.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def chains(atoms):
    """The chains as index arrays (the connected components of the bonds), and every atom's bonds with their cell offsets."""
    from ase.neighborlist import neighbor_list
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components

    n = len(atoms)
    i, j, S = neighbor_list("ijS", atoms, BONDS)
    count, label = connected_components(coo_matrix((np.ones(len(i)), (i, j)), shape=(n, n)), directed=False)
    if count != CHAINS:
        raise SystemExit(f"{count} chains found where {CHAINS} are expected")
    bonds = {k: [] for k in range(n)}
    for a, b, s in zip(i, j, S):
        bonds[a].append((b, s))
    return [np.flatnonzero(label == c) for c in range(count)], bonds


def whole(atoms):
    """Each chain's atoms where a walk along its bonds puts them (no chain cut by a cell face across the chains), the chains'
    members, and the unit vector along the chains."""
    cell = atoms.cell.array
    pos = atoms.positions.copy()
    members, bonds = chains(atoms)
    out = pos.copy()
    for atoms_of in members:
        un = {atoms_of[0]: pos[atoms_of[0]]}          # a spanning tree of the bonds: each atom placed once
        stack = [atoms_of[0]]
        while stack:
            k = stack.pop()
            for m, s in bonds[k]:
                if m not in un:
                    un[m] = un[k] + (pos[m] + s @ cell - pos[k])
                    stack.append(m)
        out[atoms_of] = np.array([un[k] for k in atoms_of])
    return out, members, cell[2] / np.linalg.norm(cell[2])


def turned(atoms, angle):
    """Every chain turned by `angle` degrees about the line through its centre along the third cell vector; the cell kept."""
    pos, members, c = whole(atoms)
    th = np.radians(angle)
    for atoms_of in members:
        r = pos[atoms_of]
        p0 = r.mean(axis=0)
        d = r - p0
        par = np.outer(d @ c, c)
        perp = d - par
        pos[atoms_of] = p0 + par + perp * np.cos(th) + np.cross(c, perp) * np.sin(th)
    b = atoms.copy()
    b.positions = pos
    b.wrap()
    return b


def sheared(atoms, sa, sb):
    """The first two cell vectors moved `sa` and `sb` monomers along the third, every chain left whole where it is: its
    neighbours along those two edges then sit `sa` and `sb` monomers along the chains."""
    pos, _, _ = whole(atoms)
    cell = atoms.cell.array.copy()
    monomer = cell[2] / (int(atoms.info["n_monomers"]) / CHAINS)
    cell[0] += sa * monomer
    cell[1] += sb * monomer
    b = atoms.copy()
    b.positions = pos
    b.set_cell(cell, scale_atoms=False)
    b.wrap()
    return b


def spread(atoms, factor):
    """The cell widened across the chains by `factor`, each chain carried whole with its centre: a start whose chains touch
    is given room, and the relaxation takes it back."""
    pos, members, c = whole(atoms)
    cell = atoms.cell.array
    wide = cell.copy()
    for row in (0, 1):
        along = (cell[row] @ c) * c
        wide[row] = along + factor * (cell[row] - along)
    for atoms_of in members:
        centre = pos[atoms_of].mean(axis=0)
        pos[atoms_of] += np.linalg.solve(cell.T, centre) @ wide - centre
    b = atoms.copy()
    b.positions = pos
    b.set_cell(wide, scale_atoms=False)
    b.wrap()
    return b


def lateral(atoms):
    """The first two cell vectors across the chains (their parts along the third taken out), and the unit vector along the chains."""
    cell = atoms.cell.array
    c = cell[2] / np.linalg.norm(cell[2])
    return np.array([cell[i] - (cell[i] @ c) * c for i in (0, 1)]), c


def like(atoms, upright, low):
    """`atoms` (a composition's upright cell) given what takes another composition's `upright` cell to its `low` cell: every
    chain turned by the same angle, the two cell vectors across the chains given the same combination of themselves (the
    same strain of the lattice), and the same monomers along the chains; each chain carried whole with its centre."""
    from ase.io import write  # noqa: F401  (ase is needed for the cells this takes)

    turn = dipole_from_a(low) - dipole_from_a(upright)
    b = turned(atoms, turn)
    up, _ = lateral(upright)
    lo, _ = lateral(low)
    combination = lo @ np.linalg.pinv(up)                        # low's vectors across the chains as combinations of upright's
    mine, c = lateral(b)
    tool = _path_tool()
    shear = np.array(tool._shear(atoms)) + np.array(tool._shear(low)) - np.array(tool._shear(upright))
    monomer = np.linalg.norm(b.cell.array[2]) / (int(b.info["n_monomers"]) / CHAINS)
    pos, members, _ = whole(b)
    cell = b.cell.array
    new = cell.copy()
    new[:2] = combination @ mine + np.outer(shear * monomer, c)
    for atoms_of in members:
        centre = pos[atoms_of].mean(axis=0)
        f = np.linalg.solve(cell.T, centre)
        along = (centre @ c) * c                                  # the chain keeps its height: only its place across the chains follows the lattice
        moved = f @ new
        pos[atoms_of] += (moved - (moved @ c) * c) - (centre - along)
    out = atoms.copy()
    out.positions = pos
    out.set_cell(new, scale_atoms=False)
    out.wrap()
    return out, float(turn)


def closest_between_chains(atoms, reach=2.6):
    """The shortest distance between atoms of different chains, A (None when none is within `reach`)."""
    from ase.neighborlist import neighbor_list

    members, _ = chains(atoms)
    label = np.zeros(len(atoms), dtype=int)
    label[members[1]] = 1
    i, j, d = neighbor_list("ijd", atoms, reach)
    apart = d[label[i] != label[j]]
    return float(apart.min()) if len(apart) else None


def dipole_from_a(atoms):
    """The direction, in degrees from the first cell vector across the chains, of the sum of every C-H bond less every C-F
    and C-Cl bond: where the chains' dipoles point.  (None for a cell with no net sum, like a helix cell's.)"""
    from ase.neighborlist import neighbor_list

    cell = atoms.cell.array
    c = cell[2] / np.linalg.norm(cell[2])
    i, j, D = neighbor_list("ijD", atoms, 1.75)
    v = np.zeros(3)
    for a, b, d in zip(i, j, D):
        if atoms.symbols[a] == "C" and atoms.symbols[b] in ("F", "H", "Cl"):
            v += d * (-1.0 if atoms.symbols[b] in ("F", "Cl") else 1.0)
    v -= c * (v @ c)
    if np.linalg.norm(v) / len(atoms) < 0.05:
        return None
    a = cell[0] - c * (cell[0] @ c)
    a /= np.linalg.norm(a)
    v /= np.linalg.norm(v)
    return float(np.degrees(np.arctan2(np.cross(a, v) @ c, a @ v)))


def gauche(atoms, tool=None):
    """How many backbone torsions are gauche (within 120 degrees of cis)."""
    tool = tool or _path_tool()
    return int((np.abs(tool.torsions_deg(atoms, tool.backbone_torsions(atoms))) < 120.0).sum())


def _ints(value):
    """A frame's "1 -1" as two integers, whether the frame file's reader left it a string or made it an array."""
    return [int(v) for v in (value.split() if isinstance(value, str) else np.atleast_1d(value).ravel())]


def _tag(v):
    return "z0" if v == 0 else f"{'p' if v > 0 else 'm'}{abs(v)}"


def name_of(label, conformation, sa, sb, turn):
    return f"reg_{label}_{conformation}_{_tag(sa)}{_tag(sb)}_t{_tag(turn)}"


def build(args):
    from ase.io import read, write

    labelled = {str(a.info["frame"]): a for f in args.labels for a in read(f, index=":")}
    cells = {}
    for item in args.cell:
        key, frame = item.split("=", 1)
        cells[key] = labelled[frame]
    out = []
    for item in args.variant:
        label, conformation, spec = item.split(":")
        start = cells[f"{label}:{conformation}"]
        if spec == "like":           # what takes --like's first cell to its second, given to this cell: turn, lattice and monomers along the chains
            b, turn = like(start, labelled[args.like[0]], labelled[args.like[1]])
            sa, sb = (int(round(v)) for v in np.array(_path_tool()._shear(b)) - np.array(_path_tool()._shear(start)))
            turn = int(round((turn + 180) % 360 - 180))
            name = f"reg_{label}_{conformation}_like_{args.like_name}"
        else:
            sa, sb, turn = (int(v) for v in spec.split(","))
            b = start.copy()
            if turn:
                b = turned(b, turn)
            if sa or sb:
                b = sheared(b, sa, sb)
            name = name_of(label, conformation, sa, sb, turn)
        room = 1.0
        while (closest_between_chains(b) or 9.0) < args.contact and room < 1.2:   # chains that touch after the move: widen the cell until they do not
            room += 0.02
            b = spread(b, 1.02)
        b.info = {k: start.info[k] for k in ("n_monomers", "kmesh") if k in start.info}
        b.info.update(frame=name, family="registry", conformation={"tt": "all-trans", "hx": "helix"}[conformation], registry_of=str(start.info["frame"]),
                      registry_label=label, registry_conformation=conformation, registry_shear=f"{sa} {sb}", registry_turn=turn, registry_widened=round(room, 2),
                      source=f"examples/registry_cells.py: {start.info['frame']}" + (f", given what takes {args.like[0]} to {args.like[1]}:" if spec == "like" else "")
                             + (f", every chain turned {turn:+d} degrees about its axis" if turn else "")
                             + (f", its neighbours along the two cell edges {sa:+d} and {sb:+d} monomers along the chains" if sa or sb else "")
                             + (f", the cell widened across the chains by {room:.2f}" if room > 1.0 else ""),
                      note="a start, not relaxed")
        for k in list(b.arrays):
            if k not in ("numbers", "positions"):
                del b.arrays[k]
        b.calc = None
        out.append(b)
        close = closest_between_chains(b)
        print(f"{name}: from {start.info['frame']}; the closest contact between chains {'beyond 2.6' if close is None else f'{close:.2f}'} A"
              + (f", the cell widened by {room:.2f}" if room > 1.0 else ""))
    names = [str(a.info["frame"]) for a in out]
    if len(set(names)) != len(names):
        raise SystemExit("a variant is given twice")
    write(args.out, out, format="extxyz")
    print(f"{len(out)} starts written to {args.out}")


def table(args):
    from ase.io import read, write

    tool = _path_tool()
    labelled = {str(a.info["frame"]): a for f in args.labels for a in read(f, index=":")}
    starts = read(args.starts, index=":")
    rows, frames = {}, []
    for start in starts:
        name = str(start.info["frame"])
        found = sorted(Path(args.relaxed).glob(f"{name}_in_*.start.xyz"))
        if not found:
            rows[name] = {"relaxed": False}
            continue
        a = read(found[0])
        ref = labelled[str(start.info["registry_of"])]
        n = int(a.info["n_monomers"]) / CHAINS
        note = str(a.info.get("note", ""))
        conformation = str(start.info["registry_conformation"])
        g, g_ref = gauche(a, tool), gauche(ref, tool)
        row = {"relaxed": True, "frame": str(a.info["frame"]), "label": str(start.info["registry_label"]), "conformation": conformation, "of": str(start.info["registry_of"]),
               "shear_asked": _ints(start.info["registry_shear"]), "turn_asked": int(start.info["registry_turn"]),
               "converged": "NOT converged" not in note, "steps": int(note.split("(")[-1].split(" steps")[0]) if " steps" in note else None,
               "energy": (float(a.info["potential_eV_per_monomer"]) - float(ref.info["potential_eV_per_monomer"])) * n * KCAL,
               "shear_along_chains_monomers": tool._shear(a), "shear_along_chains_monomers_labelled": tool._shear(ref),
               "dipole_from_a_deg": dipole_from_a(a), "dipole_from_a_deg_labelled": dipole_from_a(ref),
               "gauche_torsions": g, "gauche_torsions_labelled": g_ref, "conformation_kept": g == 0 if conformation == "tt" else abs(g - g_ref) <= 2,
               "cell_lengths": [float(v) for v in a.cell.cellpar()[:3]]}
        rows[name] = row
        frames.append(a)
        print(f"{name:34s} {row['energy']:+7.2f} {'at rest' if row['converged'] else 'NOT at rest'} ({row['steps']} steps), shear {row['shear_along_chains_monomers'][0]:+.2f} "
              f"{row['shear_along_chains_monomers'][1]:+.2f}, gauche {g} (labelled {g_ref})" + ("" if row["dipole_from_a_deg"] is None else f", dipole {row['dipole_from_a_deg']:+.0f} deg from a"))
    out = {"what": "cells made from labelled cells by examples/registry_cells.py (neighbours moved whole monomers along the chains, chains turned about their axes) and "
                   "relaxed in the potential: energy is kcal/mol per chain against the labelled cell each was made from, in the potential with D3",
           "starts": str(args.starts), "rows": rows}
    Path(args.out).write_text(json.dumps(out, indent=1) + "\n")
    if args.frames_out:
        write(args.frames_out, frames, format="extxyz")
    print(f"{sum(r['relaxed'] for r in rows.values())} of {len(rows)} relaxed; written {args.out}")


def chosen(rows, within, below, unsettled_below=None, skip=()):
    """The lowest cell of each composition and conformation that is worth a DFT label, the lowest first.  A cell that has
    not come to rest counts only when it is already below `unsettled_below` (its energy can only fall further)."""
    best = {}
    for name, r in rows.items():
        if not (r.get("relaxed") and r["conformation_kept"]) or r["frame"] in skip:
            continue
        if not r["converged"] and (unsettled_below is None or r["energy"] > unsettled_below):
            continue
        key = f"{r['label']}:{r['conformation']}"
        if r["energy"] > below.get(key, within):
            continue
        if key not in best or r["energy"] < best[key][1]:
            best[key] = (r["frame"], r["energy"])
    return sorted(best.values(), key=lambda v: v[1])


def choose(args):
    rows = json.loads(Path(args.table).read_text())["rows"]
    below = {item.split("=")[0]: float(item.split("=")[1]) for item in args.below}
    picked = chosen(rows, args.within, below, args.unsettled_below, set(args.skip))
    print(" ".join(name for name, _ in picked[0::2]))
    print(" ".join(name for name, _ in picked[1::2]))
    for name, e in picked:
        print(f"# {name}: {e:+.2f}", file=sys.stderr)


def unsettled(args):
    """The starts whose cells have not come to rest and are below --below, the lowest first: the ones worth relaxing on."""
    rows = json.loads(Path(args.table).read_text())["rows"]
    left = sorted((r["energy"], name) for name, r in rows.items() if r.get("relaxed") and not r["converged"] and r["conformation_kept"] and r["energy"] < args.below)
    print(" ".join(name for _, name in left))


def margins(args):
    from ase.io import read

    labelled = {}
    for f in args.labels:
        for a in read(f, index=":"):
            labelled[str(a.info["frame"])] = a
    groups = {}
    for item in list(args.cell) + list(args.also):
        key, frame = item.split("=", 1)
        groups.setdefault(key, []).append(frame)
    given = {item.split("=", 1)[0]: item.split("=", 1)[1] for item in args.cell}
    for name, a in labelled.items():      # the cells this tool made, labelled since: each knows what it was made from
        if "registry_label" in a.info:
            groups.setdefault(f"{a.info['registry_label']}:{a.info['registry_conformation']}", []).append(name)
    rows = dict(item.split("=", 1) for item in args.row)
    out = {"what": "margins of fourteen-monomer cells, all-trans less helix in kcal/mol per monomer, in DFT (PBE-D3(BJ) single points on cells the potential relaxed): "
                   "margin_upright between the two labelled cells the tables' search built, margin between the lowest labelled all-trans cell and the lowest labelled helix "
                   "cell of any registry. cells: kcal/mol per chain against the upright cell of the same conformation, DFT and the potential",
           "compositions": {}, "rows": {}}
    for label, row in rows.items():
        rec = {"row": row, "cells": {}}
        low = {}
        for conformation in ("tt", "hx"):
            ref = labelled[given[f"{label}:{conformation}"]]
            n = int(ref.info["n_monomers"]) / CHAINS
            cells = {}
            for frame in dict.fromkeys(groups[f"{label}:{conformation}"]):
                a = labelled[frame]
                cells[frame] = {"dft": (float(a.info["REF_energy"]) - float(ref.info["REF_energy"])) / CHAINS * KCAL,
                                "potential": (float(a.info["potential_eV_per_monomer"]) - float(ref.info["potential_eV_per_monomer"])) * n * KCAL,
                                "dft_eV_per_monomer": float(a.info["REF_energy"]) / int(a.info["n_monomers"]),
                                "dft_largest_force_eV_per_A": float(np.abs(a.arrays["REF_forces"]).max()), "upright": frame == given[f"{label}:{conformation}"]}
            rec["cells"][conformation] = cells
            low[conformation] = min(cells, key=lambda k: cells[k]["dft"])
            rec[f"lowest_{conformation}"] = low[conformation]
        e = lambda conformation, frame: rec["cells"][conformation][frame]["dft_eV_per_monomer"]
        rec["margin_upright"] = (e("tt", given[f"{label}:tt"]) - e("hx", given[f"{label}:hx"])) * KCAL
        rec["margin"] = (e("tt", low["tt"]) - e("hx", low["hx"])) * KCAL
        out["compositions"][label] = rec
        out["rows"][row] = {"margin": rec["margin"], "margin_upright": rec["margin_upright"]}
        print(f"{label:8s} margin {rec['margin']:+.3f} (upright cells {rec['margin_upright']:+.3f}); lowest all-trans {low['tt']} "
              f"({rec['cells']['tt'][low['tt']]['dft']:+.2f} per chain), lowest helix {low['hx']} ({rec['cells']['hx'][low['hx']]['dft']:+.2f})")
    Path(args.out).write_text(json.dumps(out, indent=1) + "\n")
    print(f"written {args.out}")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="step", required=True)
    s = sub.add_parser("build")
    s.add_argument("--labels", nargs="+", required=True, help="files of labelled frames")
    s.add_argument("--cell", nargs="+", required=True, metavar="LABEL:tt|hx=FRAME", help="a composition's labelled all-trans (tt) or helix (hx) cell")
    s.add_argument("--variant", nargs="+", required=True, metavar="LABEL:tt|hx:SA,SB,TURN", help="monomers along the chains for the two cell edges, and degrees every chain "
                   "is turned; or LABEL:tt:like, the cell given what takes --like's first cell to its second")
    s.add_argument("--like", nargs=2, default=None, metavar="FRAME", help="a composition's upright cell and its low cell (frames of --labels): the turn, the lattice's strain and "
                   "the monomers along the chains between them are what a `like` variant is given")
    s.add_argument("--like-name", default="film", help="the name a `like` variant carries")
    s.add_argument("--contact", type=float, default=1.9, help="A: a start whose chains come closer than this is widened across the chains until they do not")
    s.add_argument("--out", required=True)
    s = sub.add_parser("table")
    s.add_argument("--starts", required=True)
    s.add_argument("--relaxed", required=True, help="the directory of examples/mlip_train.py relax")
    s.add_argument("--labels", nargs="+", required=True)
    s.add_argument("--out", required=True)
    s.add_argument("--frames-out", default="", help="the relaxed cells as one file of frames, for a DFT single point each")
    s = sub.add_parser("choose")
    s.add_argument("--table", required=True)
    s.add_argument("--within", type=float, default=0.3, help="kcal/mol per chain above its labelled cell a cell may be and still be asked of DFT")
    s.add_argument("--below", nargs="*", default=[], metavar="LABEL:tt|hx=KCAL", help="for this composition and conformation, only a cell below this")
    s.add_argument("--unsettled-below", type=float, default=None, help="a cell that has not come to rest is taken too when it is already below this (default: never)")
    s.add_argument("--skip", nargs="*", default=[], help="relaxed frames not to choose (in a DFT line already)")
    s = sub.add_parser("unsettled")
    s.add_argument("--table", required=True)
    s.add_argument("--below", type=float, required=True)
    s = sub.add_parser("margins")
    s.add_argument("--labels", nargs="+", required=True, help="files of DFT-labelled frames")
    s.add_argument("--cell", nargs="+", required=True, metavar="LABEL:tt|hx=FRAME", help="each composition's upright all-trans and helix cell")
    s.add_argument("--also", nargs="*", default=[], metavar="LABEL:tt|hx=FRAME", help="other labelled cells of a composition and conformation that this tool did not make")
    s.add_argument("--row", nargs="+", required=True, metavar="LABEL=ROW", help="the tables' row a composition is")
    s.add_argument("--out", required=True)
    args = p.parse_args()
    {"build": build, "table": table, "choose": choose, "unsettled": unsettled, "margins": margins}[args.step](args)


if __name__ == "__main__":
    main()
