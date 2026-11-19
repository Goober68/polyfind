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
BOND = {"C": 0.77, "H": 0.37, "F": 0.71, "Cl": 0.99, "N": 0.75}


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


def backbone_torsions(atoms) -> dict:
    """``{dihedrals, bonds}`` of a lone chain's frame: one backbone dihedral a bond of the repeat, in
    degrees, walking up the chain, and the two groups each turns between (``CF2-CH2``)."""
    ST = _load("chain_stretch")
    order, _ = ST.backbone(atoms)
    n = sum(1 for e in atoms.symbols if e == "C")
    d = atoms.get_all_distances(mic=True)
    kind = []
    for c in order:
        near = [atoms.symbols[k] for k in np.flatnonzero(d[c] < 1.6) if atoms.symbols[k] != "C"]
        kind.append("C" + "".join(f"{e}{near.count(e) if near.count(e) > 1 else ''}" for e in ("H", "F") if e in near))
    return {"dihedrals": ST.dihedrals(atoms)[:n], "bonds": [f"{kind[k + 1]}-{kind[k + 2]}" for k in range(n)]}


def split(frames_dir: Path, name: str, relaxed_dir: Path | None = None) -> dict:
    """``{cell, chain, between}`` in kcal/mol per monomer and the largest DFT force on the lone chain
    (eV/A); with ``relaxed_dir``, also the lone chain relaxed (``lone``, its last labelled step
    there) and the cell against it (``binding``)."""
    from ase.io import read

    lone = [read(f) for f in sorted(frames_dir.glob(f"{name}_chain*.dft.xyz")) if f.name[len(name) + 6:-8].isdigit()]
    if not lone:
        raise FileNotFoundError(f"no labelled chain of {name} in {frames_dir}")
    info = lone[0].info
    if len(lone) not in (1, info["chains_in_cell"]):
        raise ValueError(f"{name}: {len(lone)} of its {info['chains_in_cell']} chains are labelled")
    chain = float(np.mean([a.info["REF_energy"] / a.info["n_monomers"] for a in lone])) * EV
    cell = info["cell_eV"] / info["cell_monomers"] * EV
    out = {"cell": cell, "chain": chain, "between": cell - chain, "largest_force_on_the_lone_chain": float(max(np.abs(a.arrays["REF_forces"]).max() for a in lone)),
           "repeat": float(lone[0].cell[2, 2]), "monomers_in_a_chain": int(info["n_monomers"]), **backbone_torsions(lone[0])}
    if relaxed_dir is not None:
        steps = sorted(relaxed_dir.glob(f"{name}_chain0_in_*.dft.xyz"), key=lambda f: len(f.name))
        if steps:
            r = read(steps[-1])
            e = r.info["REF_energy"] / r.info["n_monomers"] * EV
            out.update(lone=e, binding=cell - e, relaxation=e - chain, lone_steps=len(steps), lone_repeat=float(r.cell[2, 2]),
                       lone_predicted=r.info["dft_predicted_eV_per_monomer"] * EV - e,
                       lone_largest_force=float(np.abs(r.arrays["REF_forces"]).max()), lone_stress_zz_GPa=float(np.asarray(r.info["REF_stress"]).reshape(3, 3)[2, 2] * 160.2177))
    return out


def dft(args):
    out = {"what": "DFT cells split into their chains alone (in a 16 A box) and the rest; kcal/mol per monomer, each cell against its reference. "
                   "inside_a_chain / between_chains: the chain with its atoms where the cell has them. lone_chain / binding: the chain relaxed alone "
                   "(examples/mlip_train.py relax --dft-gradient, then a single point), and the cell against that.", "cells": {}, "pairs": {}}
    cells = {}
    for item in args.pairs:
        for name in item.split("="):
            if name not in cells:
                cells[name] = split(args.frames, name, args.relaxed)
    out["cells"] = cells
    for item in args.pairs:
        name, ref = item.split("=")
        a, b = cells[name], cells[ref]
        out["pairs"][f"{name} minus {ref}"] = r = {"cell": a["cell"] - b["cell"], "inside_a_chain": a["chain"] - b["chain"], "between_chains": a["between"] - b["between"]}
        text = f"{name} minus {ref}: cell {r['cell']:+.3f} = inside a chain {r['inside_a_chain']:+.3f} + between chains {r['between_chains']:+.3f}"
        if "lone" in a and "lone" in b:
            r.update(lone_chain=a["lone"] - b["lone"], binding=a["binding"] - b["binding"])
            text += f"; = the chains relaxed alone {r['lone_chain']:+.3f} + binding {r['binding']:+.3f}"
        print(text)
    for name, c in cells.items():
        print(f"  {name}: between {c['between']:+.3f}, largest force on the lone chain {c['largest_force_on_the_lone_chain']:.2f} eV/A"
              + (f"; relaxed alone {c['relaxation']:+.3f} (predicted {c['relaxation'] + c['lone_predicted']:+.3f}), binding {c['binding']:+.3f}, "
                 f"largest force left {c['lone_largest_force']:.3f} eV/A, repeat {c['repeat']:.4f} -> {c['lone_repeat']:.4f} A" if "lone" in c else ""))
    if args.out:
        Path(args.out).write_text(json.dumps(out, indent=1))
        print(f"written {args.out}")


def _load(name):
    spec = importlib.util.spec_from_file_location(name, Path(__file__).with_name(f"{name}.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def polyfind_split(SP, model, x, ch, tilts=True) -> dict:
    """A cell's energy under the parameter set in force, split as DFT's is: each chain alone on a
    16 A square lattice with its coordinates as the cell has them (one-chain packers, the same
    terms), and the cell minus those.  kcal/mol per cell.  With ``tilts`` true the pendant groups
    are relaxed as rigid tilts in the cell first (the relaxation DFT's cell has) and the lone
    chains carry the tilted groups; with an array (degrees, as ``tilts_deg`` came back) the groups
    are set to those tilts and nothing is relaxed: a cell frozen at an earlier parameter set's
    geometry, which is right to first order in how far the parameters have moved since."""
    from polyfind import mixedcell as MC
    from polyfind import pendants as PD

    polymer = model.polymer
    latn, frames, offs, P0, built = model.build(x)
    cz, n = float(latn[2, 2]), model.n
    cell = model.physical_energy(x)
    pk1 = SP.pack_mod.CrystalPacker(ch, n_chains=1, **SP.DK.KW)
    st = model.state(x)
    lone = 0.0
    for k in range(model.nc):
        one = MC.MixedModel(polymer, model.name, pk1, model.flips[k:k + 1], [ch], angle_stiffness=model.angle_stiffness)
        lone += one.physical_energy(np.concatenate([[BOX, BOX, 90.0, st.lattice[3]], [0.0, 0.0, st.place[k, 2], 0.0], st.tors[k], st.angs[k]]))
    out = {"cell": cell, "chain": lone, "between": cell - lone}
    if tilts is not False and tilts is not None:
        vals = [PD.pinned_chain_valence(b[0], polymer, None, cz=cz)[0] for b in built]
        J, atoms = PD.belt_site_constraint(model.pk, P0, cz=cz)
        if tilts is True:
            dE, deg, groups = MC.tilt_relaxation(model, x)
        else:
            deg = np.asarray(tilts, dtype=float)
            el = np.asarray(list(model.pk.elements) * model.pk.n_chains)
            groups = PD.pendant_groups(latn, P0[atoms], [str(e) for e in el[atoms]])
            if deg.shape != (len(groups), 3):
                raise ValueError(f"tilts for {deg.shape[0]} pendant groups, and the cell has {len(groups)}")
        Pt = P0 + (J @ PD.group_displacements(groups, np.radians(deg).ravel(), len(atoms))).reshape(-1, 3)
        if tilts is not True:
            dE = float(MC.placed_energy_and_grad(model.pk, frames, Pt, latn, valences=vals)[0]) - float(MC.placed_energy_and_grad(model.pk, frames, P0, latn, valences=vals)[0])
        box = np.diag([BOX, BOX, cz])
        moved = 0.0
        for k in range(model.nc):
            sl = slice(k * n, (k + 1) * n)
            e = [MC.placed_energy_and_grad(pk1, frames[k:k + 1], P[sl], box, valences=vals[k:k + 1])[0] for P in (P0, Pt)]
            moved += float(e[1]) - float(e[0])
        out.update(tilt_cell=float(dE), tilt_chain=moved, cell_tilted=cell + float(dE), chain_tilted=lone + moved, between_tilted=cell + float(dE) - lone - moved,
                   tilts_deg=np.asarray(deg, dtype=float))
    return out


def polyfind_lone(SP, model, x, ch, tilts: bool = True) -> dict:
    """Each chain of the cell relaxed alone on the 16 A lattice under the parameter set in force:
    its torsions, backbone angles and repeat free from where the cell has them, then (``tilts``)
    its pendant groups.  ``{lone, repeats, torsions}``, the energy summed over the chains."""
    from polyfind import mixedcell as MC

    pk1 = SP.pack_mod.CrystalPacker(ch, n_chains=1, **SP.DK.KW)
    st = model.state(x)
    total, repeats, torsions = 0.0, [], []
    for k in range(model.nc):
        one = MC.MixedModel(model.polymer, model.name, pk1, model.flips[k:k + 1], [ch], angle_stiffness=model.angle_stiffness)
        y = np.concatenate([[BOX, BOX, 90.0, st.lattice[3]], [0.0, 0.0, st.place[k, 2], 0.0], st.tors[k], st.angs[k]])
        free = np.concatenate([one.index("c"), one.index("internal")])
        y1, info = one.relax(y, free=free, closed=[0], maxiter=800, patience=60, patience_tol=1e-6)
        if info["stalled"] or info["gradient_max"] > 1e-3:
            y1, info = one.relax_on_surface(y1, free=free, closed=[0], gtol=1e-4, maxiter=300, patience=60)
        e = one.physical_energy(y1)
        if tilts:
            e += float(MC.tilt_relaxation(one, y1)[0])
        total += e
        s1 = one.state(y1)
        repeats.append(float(s1.lattice[3]))
        torsions.append([float(v) for v in s1.tors[0]])
    return {"lone": total, "repeats": repeats, "torsions": torsions}


PLAIN = "pvdf"   # --row: plain PVDF's two packings (examples/interchain_frames.py) in place of a row's polished cells


def polyfind(args):
    from polyfind import crystalfit as C
    from polyfind import fitting as F

    CA = _load("candidate_acceptance")
    SP = CA.SP
    plain = args.row == PLAIN
    if plain:
        IF = _load("interchain_frames")
        table = json.loads((Path(__file__).resolve().parents[1] / "docs" / "predictions" / "interchain_belt4.json").read_text())["packings"]
    else:
        pol, ref = json.load(open(args.polished or SP.POLISHED)), json.load(open(args.references or SP.REFERENCES))
    out = json.loads(Path(args.out).read_text()) if args.out and Path(args.out).exists() else {"what": "polyfind's relaxed upright cells split as the DFT cells are (examples/chain_split.py); kcal/mol per monomer", "sets": {}}
    sets = {}
    for item in args.fit:
        name, path = item.split("=", 1)
        sets[name] = F.PRESET_PARAMETERS[path] if path in F.PRESET_PARAMETERS else C.parameters(np.array(json.loads(Path(path).read_text())["x"], dtype=float))
    confs = [c for c in (args.confs or ("beta,alpha" if plain else "TG'TG,TTTT,GTG'T")).split(",") if c]
    for name, prm in sets.items():
        rows, frozen = {}, {}
        previous = SP.C.SITES_FOR
        SP.C.SITES_FOR = prm
        try:
            with prm.applied():
                if not plain:
                    letters, _ = SP.DK.sequence(args.row)
                    polymer = SP.C.make_polymer(args.row, letters)
                for conf in confs:
                    if plain:
                        model, x0, ch, per = IF.model_of(conf, table[conf].get("cell"), check=False)
                    else:
                        model, x0, ch = SP.sheared_model(polymer, SP.stored_cell(args.row, conf, pol, ref), check=False)
                        per = 2 * ch.n_monomers
                    everything, _ = SP.free_sets(model)
                    xu, _ = CA.relaxed(model, x0, np.setdiff1d(everything, model.index("shear")))
                    whole = polyfind_split(SP, model, xu, ch)
                    frozen[conf] = {"x": [float(v) for v in xu], "tilts_deg": whole.pop("tilts_deg").tolist()}
                    rows[conf] = r = {k: v / per for k, v in whole.items()}
                    lone = polyfind_lone(SP, model, xu, ch)
                    r.update(lone=lone["lone"] / per, binding=r["cell_tilted"] - lone["lone"] / per, lone_repeats=lone["repeats"], lone_torsions=lone["torsions"],
                             repeat=float(model.state(xu).lattice[3]), torsions=[[float(v) for v in t] for t in model.state(xu).tors])
                    print(f"{name} {conf}: cell {r['cell_tilted']:.4f} = chain {r['chain_tilted']:.4f} + between {r['between_tilted']:+.4f}"
                          f" (pendant tilts: cell {r['tilt_cell']:+.4f}, of which inside the chains {r['tilt_chain']:+.4f}); "
                          f"the chains relaxed alone {r['lone']:.4f} ({r['lone'] - r['chain_tilted']:+.4f}), binding {r['binding']:+.4f}, repeat {r['repeat']:.3f} -> {lone['repeats'][0]:.3f} A", flush=True)
        finally:
            SP.C.SITES_FOR = previous
        base = rows[confs[0]]
        pairs = {}
        for conf in confs[1:]:
            r = rows[conf]
            pairs[f"{conf} minus {confs[0]}"] = q = {"cell": r["cell_tilted"] - base["cell_tilted"], "inside_a_chain": r["chain_tilted"] - base["chain_tilted"],
                                                     "between_chains": r["between_tilted"] - base["between_tilted"], "lone_chain": r["lone"] - base["lone"],
                                                     "binding": r["binding"] - base["binding"]}
            print(f"{name}: {conf} minus {confs[0]}: cell {q['cell']:+.3f} = inside a chain {q['inside_a_chain']:+.3f} + between chains {q['between_chains']:+.3f}"
                  f"; = the chains relaxed alone {q['lone_chain']:+.3f} + binding {q['binding']:+.3f}")
        out["sets"].setdefault(name, {})[args.row] = {"cells": rows, "pairs": pairs}
        if args.freeze:  # the relaxed cells themselves, for rows evaluated on them as the parameters move (examples/stage13_rows.py)
            state = json.loads(Path(args.freeze).read_text()) if Path(args.freeze).exists() else {}
            state[args.row] = {"relaxed_under": name, "cells": frozen}
            Path(args.freeze).write_text(json.dumps(state))
            print(f"written {args.freeze}")
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
    s.add_argument("--relaxed", type=Path, default=None, help="directory of the lone chains relaxed and labelled")
    s.add_argument("--out", default="")
    s = sub.add_parser("polyfind")
    s.add_argument("--row", default="vdf5-rvdf1-6", help=f"a row with polished cells, or '{PLAIN}' for plain PVDF's beta and alpha")
    s.add_argument("--confs", default="", help="the first is the reference; default TG'TG,TTTT,GTG'T (for plain PVDF beta,alpha)")
    s.add_argument("--fit", action="append", default=[], help="name=crystal-stage result file, or name=preset")
    s.add_argument("--out", default="")
    s.add_argument("--freeze", default="", help="write the relaxed cells (variables and pendant tilts) of the last --fit here")
    s.add_argument("--polished", default="", help="another table of polished cells (examples/polish_margins.py --out)")
    s.add_argument("--references", default="", help="its all-trans cells (examples/polish_margins.py --references)")
    args = ap.parse_args()
    {"frames": frames, "dft": dft, "polyfind": polyfind}[args.step](args)


if __name__ == "__main__":
    main()
