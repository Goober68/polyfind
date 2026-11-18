"""A relaxed cell's energy split into its chains alone and what is between them, DFT against polyfind.

DFT's cells of a row (each relaxed, one single point each) say what a conformer's cell costs and
not where the cost sits.  One more single point per cell says it: a chain of the cell taken out
with its atoms where the cell has them, alone in a 16 A box.  The cell minus its chains is the
energy between chains; the chain alone is everything inside it.  Both are at DFT's own minimum
of the cell, so neither carries the strain a single point on polyfind's geometry does
(docs/CRYSTAL_STAGE.md 5.28).

    python examples/chain_split.py frames --cells CELL.dft.xyz ... --out DIR
    python examples/mlip_frames.py label --frames DIR --pw PW --pseudo PSEUDO
    python examples/chain_split.py dft --frames DIR --pairs NAME=REFERENCE ... --out FILE.json
    python examples/chain_split.py polyfind --row vdf5-rvdf1-6 --fit name=FILE.json --out FILE.json
"""
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np

BOX = 16.0     # A, the lateral edge of a lone chain's box (examples/chain_potential.py)
EV = 23.060548
BOND = {"C": 0.77, "H": 0.37, "F": 0.71}


def chains_of(atoms):
    """``[(indices, unwrapped positions, translation)]``: the bonded chains of a cell, each unwrapped
    along its bonds, and the lattice vector that takes it onto itself."""
    from ase.neighborlist import neighbor_list

    i, j, S = neighbor_list("ijS", atoms, {(a, b): 1.2 * (BOND[a] + BOND[b]) for a in BOND for b in BOND})
    cell = atoms.cell.array
    nb = {k: [] for k in range(len(atoms))}
    for a, b, s in zip(i, j, S):
        nb[int(a)].append((int(b), s))
    seen, out = {}, []
    for start in range(len(atoms)):
        if start in seen:
            continue
        seen[start] = np.zeros(3, dtype=int)
        order, queue, period = [start], [start], None
        while queue:
            a = queue.pop()
            for b, s in nb[a]:
                image = seen[a] + s
                if b not in seen:
                    seen[b] = image
                    order.append(b)
                    queue.append(b)
                elif period is None and np.any(image != seen[b]):
                    period = image - seen[b]
        if period is None:
            raise ValueError("a bonded group of the cell does not close on its own image: not a chain")
        idx = np.array(sorted(order))
        T = period @ cell
        out.append((idx, atoms.positions[idx] + np.array([seen[k] for k in idx]) @ cell, T if T @ cell[2] > 0 else -T))
    return out


def lone_chain(atoms, idx, P, T):
    """The chain as a frame of its own: its translation along z, its axis through the middle of a 16 A box."""
    from ase import Atoms

    z = T / np.linalg.norm(T)
    x = np.cross([0.0, 1.0, 0.0], z)
    x /= np.linalg.norm(x)
    R = np.array([x, np.cross(z, x), z])
    Q = P @ R.T
    Q[:, :2] += BOX / 2.0 - Q[:, :2].mean(axis=0)
    c = float(np.linalg.norm(T))
    return Atoms(symbols=[atoms.symbols[k] for k in idx], positions=Q, cell=[BOX, BOX, c], pbc=True)


def signature(a):
    """The sorted distances between a chain's atoms, each pair at its nearest image along the chain: equal
    for equal chains, however each was cut out of its cell."""
    P = a.positions
    d = [np.linalg.norm(P[:, None] - P[None] - np.array([0.0, 0.0, k * a.cell[2, 2]]), axis=-1) for k in (-2, -1, 0, 1, 2)]
    return np.sort(np.min(d, axis=0).ravel())


def frames(args):
    from ase.io import read, write

    args.out.mkdir(parents=True, exist_ok=True)
    for f in args.cells:
        cell = read(f)
        name = cell.info["frame"]
        lone = [lone_chain(cell, *c) for c in chains_of(cell)]
        spread = max(np.abs(signature(a) - signature(lone[0])).max() if len(a) == len(lone[0]) else np.inf for a in lone)
        keep = lone if (args.all or spread > args.same) else lone[:1]
        for k, a in enumerate(keep):
            c = a.cell[2, 2]
            a.info.update(frame=f"{name}_chain{k}", family="chain", n_monomers=len(a) // 6, kmesh=f"1 1 {max(2, int(round(20.0 / c)))}",
                          source=f"chain {k} of {name}, alone", cell_frame=name, chains_in_cell=len(lone),
                          cell_eV=float(cell.info["REF_energy"]), cell_monomers=int(cell.info["n_monomers"]))
            write(args.out / f"{name}_chain{k}.start.xyz", a, format="extxyz")
        print(f"{name}: {len(lone)} chains of {len(lone[0])} atoms, repeat {lone[0].cell[2, 2]:.4f} A; their distances differ by at most {spread:.1e} A: "
              f"{'one stands for all' if len(keep) == 1 else 'each is written'}")


def split(frames_dir: Path, name: str) -> dict:
    """``{cell, chain, between}`` in kcal/mol per monomer, and the largest DFT force on the lone chain (eV/A)."""
    from ase.io import read

    lone = [read(f) for f in sorted(frames_dir.glob(f"{name}_chain*.dft.xyz"))]
    if not lone:
        raise FileNotFoundError(f"no labelled chain of {name} in {frames_dir}")
    info = lone[0].info
    if len(lone) not in (1, info["chains_in_cell"]):
        raise ValueError(f"{name}: {len(lone)} of its {info['chains_in_cell']} chains are labelled")
    chain = float(np.mean([a.info["REF_energy"] / a.info["n_monomers"] for a in lone])) * EV
    cell = info["cell_eV"] / info["cell_monomers"] * EV
    return {"cell": cell, "chain": chain, "between": cell - chain, "largest_force_on_the_lone_chain": float(max(np.abs(a.arrays["REF_forces"]).max() for a in lone)),
            "repeat": float(lone[0].cell[2, 2]), "monomers_in_a_chain": int(info["n_monomers"])}


def dft(args):
    out = {"what": "DFT cells split into their chains alone (in a 16 A box, atoms where the cell has them) and what is between the chains; "
                   "kcal/mol per monomer, each cell against its reference", "pairs": {}}
    for item in args.pairs:
        name, ref = item.split("=")
        a, b = split(args.frames, name), split(args.frames, ref)
        out["pairs"][f"{name} minus {ref}"] = {"cell": a["cell"] - b["cell"], "inside_a_chain": a["chain"] - b["chain"], "between_chains": a["between"] - b["between"],
                                              "absolute_between": [a["between"], b["between"]], "largest_force_on_the_lone_chains": [a["largest_force_on_the_lone_chain"], b["largest_force_on_the_lone_chain"]]}
        r = out["pairs"][f"{name} minus {ref}"]
        print(f"{name} minus {ref}: cell {r['cell']:+.3f} = inside a chain {r['inside_a_chain']:+.3f} + between chains {r['between_chains']:+.3f} "
              f"(between: {a['between']:+.3f}, {b['between']:+.3f}; largest force on the lone chains {a['largest_force_on_the_lone_chain']:.2f}, {b['largest_force_on_the_lone_chain']:.2f} eV/A)")
    if args.out:
        Path(args.out).write_text(json.dumps(out, indent=1))
        print(f"written {args.out}")


def _load(name):
    spec = importlib.util.spec_from_file_location(name, Path(__file__).with_name(f"{name}.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def polyfind_split(model, x) -> dict:
    """A cell's energy under the parameter set in force, split: the terms of one chain's own
    coordinates, the lattice sum with the chains on a 16 A lattice (inside a chain and its images
    as DFT's box has them), and the rest of the lattice sum (between chains)."""
    here, box = {}, {}
    model.value_and_grad(x, here)
    y = np.array(x, dtype=float)
    y[0], y[1], y[2] = BOX, 2.0 * BOX, 90.0
    model.value_and_grad(y, box)
    own = sum(v for k, v in here.items() if k != "lattice")
    return {"cell": here["lattice"] + own, "chain": box["lattice"] + own, "between": here["lattice"] - box["lattice"], "terms": {k: float(v) for k, v in here.items()}}


def polyfind(args):
    from polyfind import crystalfit as C
    from polyfind import fitting as F
    from polyfind import mixedcell as MC

    CA = _load("candidate_acceptance")
    SP = CA.SP
    pol, ref = json.load(open(SP.POLISHED)), json.load(open(SP.REFERENCES))
    out = json.loads(Path(args.out).read_text()) if args.out and Path(args.out).exists() else {"what": "polyfind's relaxed cells of a row split as the DFT cells are (examples/chain_split.py)", "sets": {}}
    sets = {}
    for item in args.fit:
        name, path = item.split("=", 1)
        sets[name] = F.PRESET_PARAMETERS[path] if path in F.PRESET_PARAMETERS else C.parameters(np.array(json.loads(Path(path).read_text())["x"], dtype=float))
    confs = [c for c in args.confs.split(",") if c]
    for name, prm in sets.items():
        rows = {}
        previous = SP.C.SITES_FOR
        SP.C.SITES_FOR = prm
        try:
            with prm.applied():
                letters, _ = SP.DK.sequence(args.row)
                polymer = SP.C.make_polymer(args.row, letters)
                for conf in confs:
                    model, x0, ch = SP.sheared_model(polymer, SP.stored_cell(args.row, conf, pol, ref), check=False)
                    per = 2 * ch.n_monomers
                    everything, _ = SP.free_sets(model)
                    upright = np.setdiff1d(everything, model.index("shear"))
                    xu, _ = CA.relaxed(model, x0, upright)
                    s = polyfind_split(model, xu)
                    tilt = float(MC.tilt_relaxation(model, xu)[0]) / per
                    rows[conf] = {"cell": s["cell"] / per, "chain": s["chain"] / per, "between": s["between"] / per, "pendant_tilts": tilt,
                                  "terms": {k: v / per for k, v in s["terms"].items()}}
                    print(f"{name} {conf}: cell {rows[conf]['cell']:.4f} = chain {rows[conf]['chain']:.4f} + between {rows[conf]['between']:+.4f}; pendant tilts {tilt:+.4f}", flush=True)
        finally:
            SP.C.SITES_FOR = previous
        base = rows[confs[0]]
        for conf in confs[1:]:
            r = rows[conf]
            print(f"{name}: {conf} minus {confs[0]}: cell {r['cell'] - base['cell']:+.3f} = inside a chain {r['chain'] - base['chain']:+.3f} + between chains {r['between'] - base['between']:+.3f}"
                  f"; with the pendants tilted {r['cell'] + r['pendant_tilts'] - base['cell'] - base['pendant_tilts']:+.3f}")
        out["sets"][name] = {"row": args.row, "cells": rows}
    if args.out:
        Path(args.out).write_text(json.dumps(out, indent=1))
        print(f"written {args.out}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="step", required=True)
    s = sub.add_parser("frames")
    s.add_argument("--cells", nargs="+", required=True, help="DFT-labelled cells (*.dft.xyz of examples/mlip_frames.py)")
    s.add_argument("--out", type=Path, required=True)
    s.add_argument("--all", action="store_true", help="write every chain of a cell even when they are the same")
    s.add_argument("--same", type=float, default=0.01, help="A: chains whose sorted distances agree to this are one")
    s = sub.add_parser("dft")
    s.add_argument("--frames", type=Path, required=True)
    s.add_argument("--pairs", nargs="+", required=True, help="CELL=REFERENCE, frame names of the cells")
    s.add_argument("--out", default="")
    s = sub.add_parser("polyfind")
    s.add_argument("--row", default="vdf5-rvdf1-6")
    s.add_argument("--confs", default="TG'TG,TTTT,GTG'T", help="the first is the reference")
    s.add_argument("--fit", action="append", default=[], help="name=crystal-stage result file, or name=preset")
    s.add_argument("--out", default="")
    args = ap.parse_args()
    {"frames": frames, "dft": dft, "polyfind": polyfind}[args.step](args)


if __name__ == "__main__":
    main()
