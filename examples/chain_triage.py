"""A first read of a chain that is not PVDF's: what DFT can say of it from small cells.

The PI, 2026-11-23: "Nothing more to be done on pvdf until we have better numbers, but can we
take a realistic stab at other chains?"  polyfind's presets and the trained potential are of the
PVDF family and say nothing of an amide or a urea.  DFT does, and a chain alone in a box, or one
chain to a cell of a crystal, is small enough for this machine.  This script makes such cells for
any chain written as a repeat unit, relaxes them in DFT and reads from each:

    the length of a repeat along the chain (a state's stroke against another state's),
    the dipole of a repeat across the chain (the charge density's moment; the field's grip),
    the force that stretches the chain (stiffness without an assumed cross-section),
    and, for a cell of a crystal, the volume of a repeat (so grip per volume, C/m2).

What it does not give: a margin between two states to the few hundredths of a kcal/mol a field
of 55 V/um can tip, a rate in a packed cell, or a loss.  Those took a potential of the chain's own
and weeks on PVDF.

    python examples/chain_triage.py build --chains pvdf,nylon3 --out /opt/dft/triage
    python examples/chain_triage.py dft --frames /opt/dft/triage --pw /opt/mm/qe/bin/pw.x --pseudo /opt/pseudo
    python examples/chain_triage.py table --frames /opt/dft/triage --out docs/predictions/chain_triage.json

A chain is its repeat unit as a SMILES string whose first atom and atom ``tail`` are the two ends
of the backbone within the repeat (``m`` heavy atoms a repeat), and a state is the backbone's
torsion angles, a pattern repeated along the bonds, with the number of repeats in one period.
"""
import argparse
import json
import os
import subprocess
import time
from pathlib import Path

import numpy as np

E_ANGSTROM_TO_DEBYE = 4.803204712
E_PER_A2_TO_C_PER_M2 = 16.02176634
EV_PER_A_TO_NN = 1.602176634
RY_BOHR3_TO_EV_A3 = 13.605693123 / 0.529177210903**3
KCAL = 23.060548

# element: (file, mass, valence charge of the pseudopotential)
PSEUDO = {
    "H": ("H_ONCV_PBE-1.0.oncvpsp.upf", 1.00794, 1),
    "C": ("C.pbe-n-kjpaw_psl.1.0.0.UPF", 12.0107, 4),
    "F": ("F.oncvpsp.upf", 18.998403163, 7),
    "N": ("N.oncvpsp.upf", 14.0067, 5),
    "O": ("O.pbe-kjpaw.UPF", 15.9994, 6),
    "Cl": ("Cl.pbe-n-rrkjus_psl.1.0.0.UPF", 35.453, 7),
}

# The Berry phase of this program (pw.x 7.5, lberry) is right only with norm-conserving files: on
# the lone PVDF chain, where the charge density's own moment is 0.0066 C/m2 over the box, it gives
# 0.0066 with a norm-conserving carbon and is a third of a quantum off with the PAW carbon above and
# with an ultrasoft one.  So polarisation runs take these instead, on the geometry the files above
# relaxed.  The carbon is the DFT suite's own (C.pbe-mt_gipaw.UPF); no norm-conserving PBE oxygen
# comes with it, and examples/chain_triage_oxygen.in makes one with the suite's generator (ld1.x).
# A water molecule's dipole with it is 1.826 D against 1.831 D with the PAW oxygen.
PSEUDO_BERRY = {**PSEUDO, "C": ("C.pbe-mt_gipaw.UPF", 12.0107, 4), "O": ("O.pbe-tm-triage.UPF", 15.9994, 6)}
NOT_NORM_CONSERVING = ("Cl",)

T, G, GM = 180.0, 60.0, -60.0
CHAINS = {
    # the check: both shapes of PVDF's chain, whose rise and dipole are known
    "pvdf": {"repeat": "CC(F)(F)", "m": 4, "tail": 1, "what": "poly(vinylidene fluoride), the check",
             "states": {"tt": ([T], 1), "tgtg": ([T, G, T, GM], 2)}},
    # an amide every third backbone atom after the nitrogen: all its dipoles one way when extended
    "nylon3": {"repeat": "NCCC(=O)", "m": 5, "tail": 3, "what": "nylon-3, poly(beta-alanine)",
               "states": {"tt": ([T], 1)}},
    "nylon5": {"repeat": "NCCCCC(=O)", "m": 7, "tail": 5, "what": "nylon-5",
               "states": {"tt": ([T], 1)}},
    # a urea and one, or three, methylenes: the largest dipole to the least volume of the three kinds
    "urea1": {"repeat": "NC(=O)NC", "m": 5, "tail": 4, "what": "polyurea-1, poly(methylene urea)",
              "states": {"tt": ([T], 1)}},
    "urea3": {"repeat": "NC(=O)NCCC", "m": 7, "tail": 6, "what": "polyurea-3",
              "states": {"tt": ([T], 1)}},
}


def oligomer(chain: dict, n: int, torsions, seed: int = 7):
    """``(symbols, positions, backbone, bonds)`` of ``n`` repeats with the backbone torsions set."""
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdMolTransforms

    mol = Chem.AddHs(Chem.MolFromSmiles(chain["repeat"] * n))
    if AllChem.EmbedMolecule(mol, randomSeed=seed) != 0:
        AllChem.EmbedMolecule(mol, randomSeed=seed, useRandomCoords=True)
    path = list(Chem.GetShortestPath(mol, 0, (n - 1) * chain["m"] + chain["tail"]))
    conf = mol.GetConformer()
    # the embedded bond angles are rough (104 to 122 degrees along one chain).  A force field with
    # every backbone torsion held trans makes the repeats alike; the state's torsions are then set
    # rigidly, so that a period can be cut out (a force field does not move a crowded helix)
    quads = [(path[i - 1], path[i], path[i + 1], path[i + 2]) for i in range(1, len(path) - 2)]  # bond i joins path[i], path[i + 1]
    for q in quads:
        rdMolTransforms.SetDihedralDeg(conf, *q, 180.0)
    props = AllChem.MMFFGetMoleculeProperties(mol)
    field = AllChem.MMFFGetMoleculeForceField(mol, props) if props is not None else AllChem.UFFGetMoleculeForceField(mol)
    hold = field.MMFFAddTorsionConstraint if props is not None else field.UFFAddTorsionConstraint
    for q in quads:
        hold(*q, False, 180.0, 180.0, 1.0e4)
    for _ in range(20):
        if field.Minimize(maxIts=2000) == 0:
            break
    for i, q in enumerate(quads, start=1):
        rdMolTransforms.SetDihedralDeg(conf, *q, float(torsions[i % len(torsions)]))
    symbols = [a.GetSymbol() for a in mol.GetAtoms()]
    bonds = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds()]
    return symbols, np.array(conf.GetPositions()), path, bonds


def periodic(chain: dict, state: str, box: float = 14.0, n: int = 8):
    """One period of the chain along z in a box ``box`` wide, cut from the middle of an oligomer."""
    from ase import Atoms

    torsions, period = chain["states"][state]
    m = chain["m"]
    symbols, pos, path, bonds = oligomer(chain, n, torsions)
    first = (n - period) // 2
    heavy = list(range(first * m, (first + period) * m))
    # backbone atoms only: two fluorines on one carbon may be numbered the other way round a period on
    step = np.array([pos[j + period * m] - pos[j] for j in heavy if j in path])
    t = step.mean(axis=0)
    spread = float(np.abs(step - t).max())
    # a planar zigzag whose bond angles differ along the repeat is an arc, not a line (an amide's turn
    # it 8 to 13 degrees a repeat): a period cut from it is a little off, and the potential settles it
    if spread > 1.2:
        raise ValueError(f"{period} repeats of this state are not a period of the chain (the repeats differ by {spread:.2f} A)")
    keep = heavy + [b if symbols[b] == "H" else a for a, b in bonds
                    if (a in heavy and symbols[b] == "H") or (b in heavy and symbols[a] == "H")]
    z = t / np.linalg.norm(t)
    x = np.cross(z, [0.0, 0.0, 1.0] if abs(z[2]) < 0.9 else [1.0, 0.0, 0.0])
    x /= np.linalg.norm(x)
    frame = np.array([x, np.cross(z, x), z])
    p = pos[keep] @ frame.T
    p[:, :2] += box / 2 - p[:, :2].mean(axis=0)
    a = Atoms([symbols[k] for k in keep], positions=p, cell=[box, box, float(np.linalg.norm(t))], pbc=True)
    a.info.update(triage_chain="", triage_state=state, triage_repeats=period, triage_kind="chain")
    return a


def model(path, dtype: str = "float64"):
    """A MACE-architecture model file on the processor (the pretrained ones are for starts only).

    ``float32`` is the precision the pretrained files were made in and about twice as fast: enough
    for a band whose images DFT then reads, where one energy of 36 atoms takes five seconds."""
    import torch
    from mace.calculators import MACECalculator

    load = torch.jit.load
    torch.jit.load = lambda f, *args, **kwargs: load(f, *args, **{**kwargs, "map_location": "cpu"})
    try:
        return MACECalculator(model_paths=str(path), device="cpu", default_dtype=dtype)
    finally:
        torch.jit.load = load


def settle(atoms, calc, fmax: float = 0.03, steps: int = 400, free=(0, 0, 1, 0, 0, 0)):
    """Relax the atoms and the cell's free strains in a pretrained potential: a start for DFT."""
    from ase.filters import FrechetCellFilter
    from ase.optimize import LBFGS

    a = atoms.copy()
    a.calc = calc
    opt = LBFGS(FrechetCellFilter(a, mask=list(free)), logfile=None)
    opt.run(fmax=fmax, steps=steps)
    a.info.update(atoms.info)
    a.info["triage_start_steps"] = int(opt.nsteps)
    a.calc = None
    return a


def kmesh(atoms) -> str:
    """About 20 A of sampling along a periodic edge, one point across a box's vacuum."""
    lengths = atoms.cell.lengths()
    vacuum = atoms.info.get("triage_kind") == "chain"
    return " ".join(str(1 if (vacuum and k < 2) else max(2, round(20.0 / lengths[k]))) for k in range(3))


def build(args):
    from ase.io import write

    calc = model(args.model)
    args.out.mkdir(parents=True, exist_ok=True)
    for name in args.chains.split(","):
        chain = CHAINS[name]
        for state in chain["states"]:
            a = periodic(chain, state, box=args.box)
            made = a.cell.lengths()[2] / a.info["triage_repeats"]
            a = settle(a, calc)
            a.info["triage_chain"] = name
            a.info["kmesh"] = kmesh(a)
            write(args.out / f"{name}_{state}.start.xyz", a, format="extxyz")
            print(f"{name}_{state}: {len(a)} atoms, {a.info['triage_repeats']} repeat(s), rise a repeat {made:.3f} A as made, "
                  f"{a.cell.lengths()[2] / a.info['triage_repeats']:.3f} A settled ({a.info['triage_start_steps']} steps), k {a.info['kmesh']}", flush=True)


def lattice_start(chain, rng, sides=(3.9, 6.2), angles=(60.0, 120.0)):
    """One chain to a cell: a random pair of edges across the chain, each with a shift along it."""
    from ase import Atoms

    c = chain.cell.lengths()[2]
    a, b = rng.uniform(*sides, size=2)
    gamma = np.radians(rng.uniform(*angles))
    turn = rng.uniform(0.0, 2 * np.pi)
    rot = np.array([[np.cos(turn), -np.sin(turn), 0.0], [np.sin(turn), np.cos(turn), 0.0], [0.0, 0.0, 1.0]])
    pos = chain.positions.copy()
    pos[:, :2] -= pos[:, :2].mean(axis=0)
    cell = np.array([[a, 0.0, rng.uniform(-c / 2, c / 2)], [b * np.cos(gamma), b * np.sin(gamma), rng.uniform(-c / 2, c / 2)], [0.0, 0.0, c]])
    out = Atoms(chain.get_chemical_symbols(), positions=pos @ rot.T, cell=cell, pbc=True)
    out.info.update(chain.info)
    out.info["triage_kind"] = "crystal"
    return out


def nearest(atoms) -> float:
    """The shortest distance between any two atoms of the cell, images counted (a collapse reads below 0.85 A)."""
    from ase.neighborlist import neighbor_list

    d = neighbor_list("d", atoms, 1.5)
    return float(d.min()) if len(d) else 9.9


def closest(atoms) -> float:
    """The shortest distance between atoms of different chains (images across the chain)."""
    from ase.neighborlist import neighbor_list

    d, shift = neighbor_list("dS", atoms, 3.0)
    across = (shift[:, 0] != 0) | (shift[:, 1] != 0)
    return float(d[across].min()) if across.any() else 9.9


def pack(args):
    """Cells of one chain each, every chain's dipole the same way: found in a pretrained potential.

    Random lattices, each relaxed whole (atoms and cell); the lowest few go on to DFT.  A cell of
    one chain is the polar packing by construction: whether the chain would rather pack with its
    neighbours turned round is a margin this does not ask."""
    from ase.io import read, write

    calc = model(args.model)
    rng = np.random.default_rng(args.seed)
    src = Path(args.frames)
    for name in args.names.split(","):
        relaxed = src / f"{name}.relaxed.xyz"
        chain = read(relaxed if relaxed.exists() else src / f"{name}.start.xyz")
        found = []
        tried = 0
        while len(found) < args.starts and tried < 40 * args.starts:
            tried += 1
            start = lattice_start(chain, rng)
            if closest(start) < 1.7:
                continue
            try:
                a = settle(start, calc, fmax=0.05, steps=args.steps, free=(1, 1, 1, 1, 1, 1))
            except (RuntimeError, ValueError, FloatingPointError) as err:  # a start that falls apart in the potential
                print(f"{name}: a start failed ({type(err).__name__})", flush=True)
                continue
            a.calc = calc
            e = float(a.get_potential_energy())
            a.calc = None
            found.append((e, a))
            print(f"{name}: start {len(found)}: {e:.4f} eV, {a.get_volume() / a.info['triage_repeats']:.1f} A3 a repeat, "
                  f"edges {' '.join(f'{v:.2f}' for v in a.cell.lengths())}, closest {closest(a):.2f} A", flush=True)
        found.sort(key=lambda r: r[0])
        kept = []
        for e, a in found:  # one of each: cells within 2 meV and 1% in volume are the same packing
            if all(abs(e - e2) > 0.002 or abs(a.get_volume() / a2.get_volume() - 1) > 0.01 for e2, a2 in kept):
                kept.append((e, a))
        for rank, (e, a) in enumerate(kept[:args.keep], start=1):
            a.info["triage_packing_rank"] = rank
            a.info["triage_packing_ev_above"] = e - kept[0][0]
            a.info["kmesh"] = kmesh(a)
            write(src / f"{name}_x{rank}.start.xyz", a, format="extxyz")
            print(f"{name}_x{rank}: {e - kept[0][0]:+.4f} eV above the lowest, {a.get_volume() / a.info['triage_repeats']:.1f} A3 a repeat, k {a.info['kmesh']}", flush=True)


def bond_step(atoms):
    """The lattice step across the chains of the shortest hydrogen bond between chains, or None.

    ``(s1, s2)``: the acceptor (O or N) is in the image ``s1`` cells along the first edge and ``s2``
    along the second from the hydrogen that points at it."""
    from ase.neighborlist import neighbor_list

    symbols = atoms.get_chemical_symbols()
    i, j, d, shift = neighbor_list("ijdS", atoms, 2.5)
    found = [(dd, int(ss[0]), int(ss[1])) for ii, jj, dd, ss in zip(i, j, d, shift)
             if symbols[ii] == "H" and symbols[jj] in ("O", "N") and (ss[0] or ss[1])]
    if not found:
        return None
    _, s1, s2 = min(found)
    return (s1, s2) if (s1, s2) > (0, 0) or (s1 == 0 and s2 > 0) else (-s1, -s2)


def rebased(atoms):
    """The same crystal with its first edge along the hydrogen bonds between chains.

    Chains bonded in a row then follow each other along the first edge and rows are stacked along
    the second, whatever pair of edges the packing search happened to end on.  A cell without
    hydrogen bonds between its chains is returned as it is."""
    step = bond_step(atoms)
    if step is None or step == (1, 0):
        return atoms
    s1, s2 = step
    x, y, g = 1, 0, s1  # the extended Euclid: s1 x + s2 y = gcd
    x1, y1, r = 0, 1, s2
    while r:
        q = g // r
        x, x1, y, y1, g, r = x1, x - q * x1, y1, y - q * y1, r, g - q * r
    if abs(g) != 1:
        raise ValueError(f"the hydrogen bond's step {step} is not one row")
    x, y = x * g, y * g
    cell = atoms.cell.array
    first, second = s1 * cell[0] + s2 * cell[1], -y * cell[0] + x * cell[1]
    second = second - np.round(second @ first / (first @ first)) * first  # the least skewed second edge
    b = atoms.copy()
    b.info.update(atoms.info)
    # the atoms stay where they are, one whole chain: wrapped one by one into the new cell they would
    # be pieces of different images, and a step between chains could no longer be read off them
    b.set_cell([first, second, cell[2]], scale_atoms=False)
    b.info["kmesh"] = kmesh(b)
    return b


def turned(atoms, edge: int):
    """The cell doubled along one edge across the chain, its second chain turned half round its own axis.

    The chain keeps its sense and its period; what it carries across itself is reversed, so the two
    chains' polarisation across the chains cancels.  The cell is first put with its first edge along
    the hydrogen bonds (``rebased``): doubled along that edge (0) every second chain of a bonded row
    is turned and the row loses its bonds; doubled along the other (1) whole rows face opposite ways
    and keep them."""
    atoms = rebased(atoms)
    n = len(atoms)
    b = atoms.repeat([2 if edge == k else 1 for k in range(3)])
    along = atoms.cell.array[2] / np.linalg.norm(atoms.cell.array[2])
    d = b.positions[n:] - b.positions[n:].mean(axis=0)
    b.positions[n:] = b.positions[n:].mean(axis=0) + 2 * np.outer(d @ along, along) - d
    mesh = mesh_of(atoms)
    mesh[edge] = max(1, (mesh[edge] + 1) // 2)
    b.info.update(atoms.info)
    b.info.update(triage_repeats=2 * int(atoms.info["triage_repeats"]), triage_kind="crystal", kmesh=" ".join(str(v) for v in mesh))
    return b


def opposed(args):
    """Each relaxed one-chain cell against the same chains with every second one turned round.

    Writes ``NAME_o1`` and ``NAME_o2`` (doubled along the first and the second edge), settled in the
    pretrained potential; ``dft`` relaxes them, and ``table`` gives each one's energy a repeat above
    its polar cell: below zero, the chain would rather not be polar."""
    from ase.io import read, write

    calc = model(args.model)
    src = Path(args.frames)
    for name in args.names.split(","):
        a = read(src / f"{name}.relaxed.xyz")
        for edge in (0, 1):
            b = turned(a, edge)
            mesh = b.info["kmesh"]
            try:
                c = settle(b, calc, fmax=0.05, steps=args.steps, free=(1, 1, 1, 1, 1, 1))
            except (RuntimeError, ValueError, FloatingPointError) as err:
                print(f"{name}_o{edge + 1}: the start fell apart in the potential ({type(err).__name__})", flush=True)
                continue
            c.calc = calc
            energy = float(c.get_potential_energy()) / c.info["triage_repeats"]
            c.calc = None
            c.info.update(kmesh=mesh, triage_opposed_of=name, triage_start_ev_a_repeat=energy)
            write(src / f"{name}_o{edge + 1}.start.xyz", c, format="extxyz")
            print(f"{name}_o{edge + 1}: {len(c)} atoms, {energy:.4f} eV a repeat in the potential, {c.get_volume() / c.info['triage_repeats']:.1f} A3 a repeat "
                  f"(the polar cell {a.get_volume() / a.info['triage_repeats']:.1f}), closest {closest(c):.2f} A, k {mesh}", flush=True)


def turned_by(atoms, edge: int, angle: float):
    """``turned`` with the second chain turned by ``angle`` (radians) instead of half round."""
    atoms = rebased(atoms)
    n = len(atoms)
    b = turned(atoms, edge)
    start = atoms.repeat([2 if edge == k else 1 for k in range(3)]).positions[n:]
    along = atoms.cell.array[2] / np.linalg.norm(atoms.cell.array[2])
    d = start - start.mean(axis=0)
    b.positions[n:] = start.mean(axis=0) + d * np.cos(angle) + np.cross(along, d) * np.sin(angle) + np.outer(d @ along, along) * (1 - np.cos(angle))
    return b


def spun(row, per_chain: int, chain: int, angle: float, along):
    """A row cell with one of its chains turned by ``angle`` about its own axis, in place."""
    part = slice(chain * per_chain, (chain + 1) * per_chain)
    d = row.positions[part] - row.positions[part].mean(axis=0)
    row.positions[part] = row.positions[part].mean(axis=0) + d * np.cos(angle) + np.cross(along, d) * np.sin(angle) + np.outer(d @ along, along) * (1 - np.cos(angle))
    return row


def band(images, calc, fmax: float, steps: int, end_steps: int = 150):
    """Relax the two ends and settle the images between them as a band, the cell held: in a pretrained potential.

    A loose band on purpose (it is the start of a DFT reading of each image, not the reading), and
    it says where it has got to: on one processor a band of 36 atoms takes tens of minutes."""
    from ase.calculators.singlepoint import SinglePointCalculator
    from ase.mep import NEB
    from ase.optimize import FIRE, LBFGS

    t0 = time.time()
    for im in images:
        im.calc = calc
    for end in (images[0], images[-1]):
        LBFGS(end, logfile=None).run(fmax=0.05, steps=end_steps)
        # an end does not move again: keep its energy, or the band asks the potential for it at every step
        end.calc = SinglePointCalculator(end, energy=float(end.get_potential_energy()), forces=end.get_forces())
    print(f"  ends relaxed, {time.time() - t0:.0f} s", flush=True)
    neb = NEB(images, climb=True, method="improvedtangent", allow_shared_calculator=True)
    opt = FIRE(neb, logfile=None)
    for _ in opt.irun(fmax=fmax, steps=steps):
        if opt.nsteps and opt.nsteps % 20 == 0:
            print(f"  band step {opt.nsteps}, {time.time() - t0:.0f} s", flush=True)
    return [float(im.get_potential_energy()) for im in images], int(opt.nsteps)


def wall(args):
    """The smallest move of a row of bonded chains: a wall between its two senses stepping one chain on.

    Four chains of a row to a cell (the first edge, along the hydrogen bonds, taken four times),
    two facing one way and two turned: two walls.  The third chain is turned back, which moves one
    wall by a chain and leaves as many bonds made as before.  The top of that path is what one chain
    pays to change sides with a wall beside it; the two ends differ only by how near the walls are."""
    from ase.io import read, write

    calc = model(args.model, dtype="float32")
    src = Path(args.frames)
    for name in args.names.split(","):
        a = rebased(read(src / f"{name}.relaxed.xyz"))
        n = len(a)
        along = a.cell.array[2] / np.linalg.norm(a.cell.array[2])
        mesh = mesh_of(a)
        mesh[0] = max(1, round(mesh[0] / 4))
        images = []
        for t in np.linspace(np.pi, 0.0, args.images):
            row = a.repeat([4, 1, 1])
            spun(row, n, 3, np.pi, along)
            spun(row, n, 2, t, along)
            row.info.update(a.info)
            row.info.update(triage_repeats=4 * int(a.info["triage_repeats"]), triage_kind="crystal", kmesh=" ".join(str(v) for v in mesh))
            images.append(row)
        energies, nsteps = band(images, calc, args.fmax, args.steps)
        reps = int(a.info["triage_repeats"])
        out = []
        for k, im in enumerate(images):
            c = im.copy()
            c.info.update(images[0].info)
            c.info.update(triage_turn_of=name, triage_turn_edge=1, triage_turn_image=k, triage_turn_ev=energies[k] - energies[0], triage_turning_repeats=reps)
            out.append(c)
        write(src / f"{name}_wall.path.xyz", out, format="extxyz")
        print(f"{name}_wall: {len(images)} images of {len(images[0])} atoms, {nsteps} steps, top {(max(energies) - energies[0]) * KCAL / reps:.2f} kcal/mol a repeat of the chain that turns "
              f"in the potential, end {(energies[-1] - energies[0]) * KCAL / reps:+.2f}; profile {' '.join(f'{(e - energies[0]) * KCAL / reps:.2f}' for e in energies)}", flush=True)


def methylated(row, per_chain: int, chain: int):
    """A row cell with the first N-H of one chain made N-CH3: that chain gives one hydrogen bond fewer.

    The hydrogen is replaced where it stands by a carbon at 1.46 A, with three hydrogens round it;
    the new atoms go at the end of the cell, so the chains' atoms keep their places in the list."""
    from ase import Atoms

    symbols = row.get_chemical_symbols()
    part = range(chain * per_chain, (chain + 1) * per_chain)
    pos = row.positions
    for n_at in part:
        if symbols[n_at] != "N":
            continue
        near = [h for h in part if symbols[h] == "H" and np.linalg.norm(pos[h] - pos[n_at]) < 1.2]
        if near:
            break
    else:
        raise ValueError("no N-H on this chain")
    h_at = near[0]
    u = (pos[h_at] - pos[n_at]) / np.linalg.norm(pos[h_at] - pos[n_at])
    carbon = pos[n_at] + 1.46 * u
    w = np.cross(u, [0.0, 0.0, 1.0] if abs(u[2]) < 0.9 else [1.0, 0.0, 0.0])
    w /= np.linalg.norm(w)
    hydrogens = [carbon + 1.09 * (u / 3 + np.sqrt(8) / 3 * (np.cos(t) * w + np.sin(t) * np.cross(u, w))) for t in (0.0, 2 * np.pi / 3, 4 * np.pi / 3)]
    out = row.copy()
    out.info.update(row.info)
    out.symbols[h_at] = "C"
    out.positions[h_at] = carbon
    out += Atoms("HHH", positions=hydrogens)
    return out


def nucleus(args):
    """The first move of a row at rest: one chain of four turning round, alone, and what a defect does to it.

    All four chains face one way; the third is turned, which makes two walls where there were none.
    With ``--methyl`` the second chain, whose N-H holds the third, has that N-H made N-CH3 first and
    the cell is let settle with it.  Beside the PVDF numbers: this is a kink forming, ``wall`` is its
    hop.  The end of the path is how far from level the turned chain sits; the top is the hump."""
    from ase.io import read, write

    calc = model(args.model, dtype="float32")
    src = Path(args.frames)
    for name in args.names.split(","):
        a = rebased(read(src / f"{name}.relaxed.xyz"))
        n = len(a)
        along = a.cell.array[2] / np.linalg.norm(a.cell.array[2])
        mesh = mesh_of(a)
        mesh[0] = max(1, round(mesh[0] / 4))
        base = a.repeat([4, 1, 1])
        base.info.update(a.info)
        base.info.update(triage_repeats=4 * int(a.info["triage_repeats"]), triage_kind="crystal", kmesh=" ".join(str(v) for v in mesh))
        tag = "nuc"
        if args.methyl:
            tag = "nucm"
            # room for the methyl before anything is relaxed: it is put where the next chain's oxygen
            # is held (C to O 1.4 A), and from that start the potential does not part the chains, it
            # folds the cell up (17.5 x 3.3 x 3.8 A with atoms 0.1 A apart, the first time).  The two
            # chains beyond the methyl are moved along the row and the row's edge grows by as much.
            base = methylated(base, n, 1)
            cell = base.cell.array.copy()
            shift = args.gap * cell[0] / np.linalg.norm(cell[0])
            base.positions[2 * n:4 * n] += shift
            cell[0] += shift
            base.set_cell(cell, scale_atoms=False)
            base = settle(base, model(args.model), fmax=0.05, steps=args.settle, free=(1, 1, 1, 1, 1, 1))
            if nearest(base) < 0.85:
                raise ValueError(f"the cell with the methyl did not settle: two atoms {nearest(base):.2f} A apart")
            along = base.cell.array[2] / np.linalg.norm(base.cell.array[2])
        images = []
        for t in np.linspace(0.0, np.pi, args.images):
            row = base.copy()
            row.info.update(base.info)
            images.append(spun(row, n, 2, t, along))
        energies, nsteps = band(images, calc, args.fmax, args.steps)
        for which, im in (("start", images[0]), ("end", images[-1])):
            if nearest(im) < 0.85:
                raise ValueError(f"the {which} of the path has two atoms {nearest(im):.2f} A apart: not a structure")
        reps = int(a.info["triage_repeats"])
        out = []
        for k, im in enumerate(images):
            c = im.copy()
            c.info.update(base.info)
            c.info.update(triage_turn_of=name, triage_turn_edge=1, triage_turn_image=k, triage_turn_ev=energies[k] - energies[0], triage_turning_repeats=reps)
            out.append(c)
        write(src / f"{name}_{tag}.path.xyz", out, format="extxyz")
        print(f"{name}_{tag}: {len(images)} images of {len(images[0])} atoms, {nsteps} steps, top {(max(energies) - energies[0]) * KCAL / reps:.2f} kcal/mol a repeat of the chain that turns "
              f"in the potential, end {(energies[-1] - energies[0]) * KCAL / reps:+.2f}; profile {' '.join(f'{(e - energies[0]) * KCAL / reps:.2f}' for e in energies)}", flush=True)


def turn(args):
    """How locked a chain is: the way from the polar cell to the one with every second chain turned.

    The cell is the polar one doubled and is held; the second chain is turned about its own axis in
    equal steps, the two ends are relaxed and the path between them settled as a band, all in the
    pretrained potential.  ``points`` then asks DFT for the energy of every image.  The top of the
    path, for each repeat of the chain that turns, is what it costs to turn a whole row of chains at
    once: an upper mark for the smallest move, which in a real crystal is a twist that travels."""
    from ase.io import read, write

    calc = model(args.model, dtype="float32")
    src = Path(args.frames)
    for item in args.names.split(","):
        name, edge = item.rsplit(":", 1)
        edge = int(edge) - 1
        a = read(src / f"{name}.relaxed.xyz")
        images = [turned_by(a, edge, t) for t in np.linspace(0.0, np.pi, args.images)]
        energies, nsteps = band(images, calc, args.fmax, args.steps)
        reps = int(a.info["triage_repeats"])
        out = []
        for k, im in enumerate(images):
            c = im.copy()
            c.info.update(images[0].info)
            c.info.update(triage_turn_of=name, triage_turn_edge=edge + 1, triage_turn_image=k, triage_turn_ev=energies[k] - energies[0], triage_turning_repeats=reps)
            out.append(c)
        write(src / f"{name}_turn{edge + 1}.path.xyz", out, format="extxyz")
        top = max(energies) - energies[0]
        print(f"{name}_turn{edge + 1}: {len(images)} images, {nsteps} steps, top {top * KCAL / reps:.2f} kcal/mol a repeat of the turning chain in the potential, "
              f"end {(energies[-1] - energies[0]) * KCAL / reps:+.2f}; profile {' '.join(f'{(e - energies[0]) * KCAL / reps:.2f}' for e in energies)}", flush=True)


def points(args):
    """DFT's energy of every image of a path ``turn`` made (the cell and the atoms as they are)."""
    from ase.io import read

    src = Path(args.frames)
    for name in args.names.split(","):
        images = read(src / f"{name}.path.xyz", index=":")
        work = src / "work" / name
        reps = int(images[0].info.get("triage_turning_repeats", int(images[0].info["triage_repeats"]) // 2))  # of the chain that turns
        energies = []
        t0 = time.time()
        for k, im in enumerate(images):
            out = work / f"image_{k}.out"
            if not (out.exists() and "JOB DONE" in out.read_text(errors="replace")):
                run(args.pw, args.np, pw_input(im, args.pseudo, work / f"scratch_{k}", "scf"), work, f"image_{k}")
            text = out.read_text(errors="replace")
            ry = [float(line.split("=")[1].split()[0]) for line in text.splitlines() if line.startswith("!")]
            energies.append(ry[-1] * 13.605693123 if ry else float("nan"))
        rel = [(e - energies[0]) * KCAL / reps for e in energies]
        record = {"name": name, "kind": "path", "of": images[0].info["triage_turn_of"], "edge": int(images[0].info["triage_turn_edge"]),
                  "turning_repeats": reps, "energy_ev": energies, "kcal_per_turning_repeat": rel, "top_kcal_per_turning_repeat": float(np.nanmax(rel)),
                  "potential_kcal_per_turning_repeat": [float(im.info["triage_turn_ev"]) * KCAL / reps for im in images]}
        (src / f"{name}.path.json").write_text(json.dumps(record, indent=1))
        print(f"{name}: DFT on {len(images)} images, top {np.nanmax(rel):.2f} kcal/mol a repeat of the turning chain, end {rel[-1]:+.2f}; "
              f"profile {' '.join(f'{v:.2f}' for v in rel)}, {time.time() - t0:.0f} s", flush=True)


def pw_input(atoms, pseudo_dir, outdir, calculation: str, free: str = "z", conv_thr: float = 1e-9, ecut=(90, 360), files=None):
    """The anchor's settings (PBE, D3(BJ) two-body, 90/360 Ry) for a relaxation of this cell."""
    files = files or PSEUDO
    symbols = atoms.get_chemical_symbols()
    species = sorted(set(symbols), key=lambda s: -symbols.count(s))
    lines = [
        "&CONTROL", f"  calculation = '{calculation}'", "  prefix = 'frame'", f"  pseudo_dir = '{pseudo_dir}'",
        f"  outdir = '{outdir}'", "  tstress = .true.", "  tprnfor = .true.", "  disk_io = 'low'", "  nstep = 200",
        "  forc_conv_thr = 5.0d-4", "  etot_conv_thr = 1.0d-5", "/", "",
        "&SYSTEM", "  ibrav = 0", f"  nat = {len(atoms)}", f"  ntyp = {len(species)}", "  input_dft = 'PBE'",
        "  occupations = 'fixed'", f"  ecutwfc = {ecut[0]}", f"  ecutrho = {ecut[1]}", "  vdw_corr = 'grimme-d3'",
        "  dftd3_version = 4", "  dftd3_threebody = .false.", "/", "",
        "&ELECTRONS", f"  conv_thr = {conv_thr:.1e}".replace("e", "d"), "  electron_maxstep = 200",
        "  mixing_mode = 'plain'", "  mixing_beta = 0.3", "/", "",
    ]
    if calculation in ("relax", "vc-relax"):
        lines += ["&IONS", "  ion_dynamics = 'bfgs'", "/", ""]
    if calculation == "vc-relax":
        # a chain's box: 0.05 kbar over 14 A squared is 0.02 nN along the chain; a crystal: 0.2 kbar is 0.2% of strain
        lines += ["&CELL", "  cell_dynamics = 'bfgs'", f"  cell_dofree = '{free}'", f"  press_conv_thr = {0.05 if free == 'z' else 0.2}", "/", ""]
    lines += ["ATOMIC_SPECIES"] + [f"{s} {files[s][1]} {files[s][0]}" for s in species]
    lines += ["", "CELL_PARAMETERS angstrom"] + [" ".join(f"{v:.10f}" for v in row) for row in atoms.cell.array]
    lines += ["", "ATOMIC_POSITIONS angstrom"]
    lines += [f"{s} {p[0]:.10f} {p[1]:.10f} {p[2]:.10f}" for s, p in zip(symbols, atoms.positions)]
    mesh = atoms.info["kmesh"]  # a string when made, an array once read back from the file
    mesh = mesh.split() if isinstance(mesh, str) else np.atleast_1d(mesh)
    lines += ["", "K_POINTS automatic", " ".join(str(int(v)) for v in mesh) + " 0 0 0", ""]
    return "\n".join(lines)


def run(binary, np_, text, work, name, pools=None):
    work.mkdir(parents=True, exist_ok=True)
    (work / f"{name}.in").write_text(text)
    mpirun = Path(binary).with_name("mpirun")  # the launcher beside the program, as examples/mlip_frames.py takes it
    cmd = [str(mpirun if mpirun.exists() else "mpirun"), "--allow-run-as-root", "--oversubscribe", "-np", str(np_), str(binary)]
    cmd += (["-nk", str(pools)] if pools else []) + ["-in", f"{name}.in"]
    env = {**os.environ, "OMP_NUM_THREADS": "1", "PATH": f"{Path(binary).parent}:{os.environ.get('PATH', '/usr/bin:/bin')}"}
    with open(work / f"{name}.out", "w") as out:
        subprocess.run(cmd, cwd=work, stdout=out, stderr=subprocess.STDOUT, check=False, env=env)
    return work / f"{name}.out"


def finished(out: Path) -> bool:
    if not out.exists():
        return False
    text = out.read_text(errors="replace")
    return "JOB DONE" in text and ("bfgs converged" in text or "calculation = 'scf'" in text or "End of self-consistent" in text) \
        and "convergence NOT achieved" not in text.split("!")[-1]


def read_final(out: Path):
    """The last geometry of a finished run with its energy (eV) and stress (eV/A3, tension positive)."""
    from ase.io import read

    a = read(out, format="espresso-out", index=-1)
    text = out.read_text(errors="replace")
    ry = [float(l.split("=")[1].split()[0]) for l in text.splitlines() if l.startswith("!")]
    block = text.rsplit("total   stress", 1)[1].splitlines()[1:4]
    # the code prints the stress pressure-signed (a cell that wants to grow reads positive), Ry/bohr3 first
    pressure = np.array([[float(v) for v in row.split()[:3]] for row in block]) * RY_BOHR3_TO_EV_A3
    return a, ry[-1] * 13.605693123, -pressure


def moment(cube: Path, atoms):
    """The cell's dipole (e A) from the valence charge density and the ions' valence charges.

    Meaningful only along an edge the density does not cross: across the box of a chain."""
    from ase.io.cube import read_cube_data

    rho, _ = read_cube_data(str(cube))  # electrons per bohr3 on a grid over the cell
    n = np.array(rho.shape)
    cell = atoms.cell.array
    dv = abs(np.linalg.det(cell)) / n.prod() / 0.529177210903**3
    frac = [np.arange(n[k]) / n[k] for k in range(3)]
    electrons = float(rho.sum() * dv)
    mu = np.zeros(3)
    for k in range(3):
        weight = rho.sum(axis=tuple(j for j in range(3) if j != k)) * dv
        mu -= (weight * frac[k]).sum() * cell[k]
    scaled = atoms.get_scaled_positions(wrap=True)
    z = np.array([PSEUDO[s][2] for s in atoms.get_chemical_symbols()])
    mu += (z[:, None] * scaled).sum(axis=0) @ cell
    return mu, electrons, float(z.sum())


def dft(args):
    from ase.io import read, write

    src = Path(args.frames)
    names = args.names.split(",") if args.names else [p.name[:-len(".start.xyz")] for p in sorted(src.glob("*.start.xyz"))]
    for name in names:
        store = src / f"{name}.json"
        got = json.loads(store.read_text()) if store.exists() else {}
        a = read(src / f"{name}.start.xyz")
        work = src / "work" / name
        chain = a.info.get("triage_kind") == "chain"
        free = "z" if chain else "all"
        t0 = time.time()
        out = work / "relax.out"
        if not finished(out):
            begin = a
            if out.exists() and out.read_text(errors="replace").count("\n!") > 1:
                # a run that was stopped: go on from where it had got to, not from the start
                begin = read(out, format="espresso-out", index=-1)
                begin.info.update(a.info)
                begin.calc = None
                out.rename(work / f"relax.stopped{len(list(work.glob('relax.stopped*')))}.out")
            run(args.pw, args.np, pw_input(begin, args.pseudo, work / "scratch", "vc-relax", free), work, "relax")
        if not finished(out):
            print(f"{name}: the relaxation did not finish ({out})", flush=True)
            continue
        b, energy, _ = read_final(out)
        b.info.update(a.info)
        write(src / f"{name}.relaxed.xyz", b, format="extxyz")
        reps = int(a.info["triage_repeats"])
        if "triage_opposed_of" in a.info:
            got["opposed_of"] = a.info["triage_opposed_of"]
        got.update(name=name, chain=a.info["triage_chain"], state=a.info["triage_state"], kind=a.info["triage_kind"], repeats=reps,
                   atoms=len(b), energy_ev=energy, cell=b.cell.array.tolist(),
                   kmesh=" ".join(str(int(v)) for v in np.atleast_1d(a.info["kmesh"])) if not isinstance(a.info["kmesh"], str) else a.info["kmesh"],
                   rise_a=float(b.cell.lengths()[2] / reps), volume_a3=float(b.get_volume() / reps),
                   relax_s=got.get("relax_s", round(time.time() - t0)))
        if chain and "dipole_across_e_a" not in got:
            cube = work / "rho.cube"
            (work / "pp.in").write_text("\n".join([
                "&INPUTPP", "  prefix = 'frame'", f"  outdir = '{work / 'scratch'}'", f"  filplot = '{work / 'rho.dat'}'", "  plot_num = 0", "/",
                "&PLOT", "  iflag = 3", "  output_format = 6", f"  fileout = '{cube}'", "/", ""]))
            with open(work / "pp.out", "w") as fh:
                subprocess.run([str(Path(args.pw).with_name("pp.x")), "-in", "pp.in"], cwd=work, stdout=fh, stderr=subprocess.STDOUT,
                               check=False, env={**os.environ, "OMP_NUM_THREADS": "1", "PATH": f"{Path(args.pw).parent}:{os.environ.get('PATH', '/usr/bin:/bin')}"})
            if cube.exists():
                mu, electrons, ions = moment(cube, b)
                # across the chain only: along it the density crosses the cell's faces and the moment means nothing
                got.update(dipole_across_e_a=(mu[:2] / reps).tolist(), electrons=electrons, ion_charge=ions)
                cube.unlink()
                (work / "rho.dat").unlink(missing_ok=True)
        # an extended chain only: a helix held 1% off its length creeps for tens of steps along its soft
        # torsions (PVDF's TGTG', 27 steps and still moving), and its stiffness is not what is asked
        if chain and a.info["triage_state"] == "tt" and "stretch" not in got and not args.no_stretch:
            rows = []
            for strain in (-args.strain, args.strain):
                c = b.copy()
                c.info.update(a.info)
                cell = c.cell.array.copy()
                cell[2] *= 1 + strain
                c.set_cell(cell, scale_atoms=True)
                tag = f"stretch_{'m' if strain < 0 else 'p'}"
                if not finished(work / f"{tag}.out"):
                    run(args.pw, args.np, pw_input(c, args.pseudo, work / f"scratch_{tag}", "relax"), work, tag)
                if finished(work / f"{tag}.out"):
                    d, e, s = read_final(work / f"{tag}.out")
                    area = float(np.linalg.norm(np.cross(d.cell[0], d.cell[1])))
                    rows.append({"strain": strain, "energy_ev": e, "tension_ev_per_a": float(s[2, 2] * area)})
            if len(rows) == 2:
                c0 = b.cell.lengths()[2]
                by_energy = (rows[0]["energy_ev"] + rows[1]["energy_ev"] - 2 * energy) / args.strain**2 / c0
                by_stress = (rows[1]["tension_ev_per_a"] - rows[0]["tension_ev_per_a"]) / (2 * args.strain)
                got["stretch"] = {"rows": rows, "strain": args.strain, "force_per_strain_nn_by_energy": by_energy * EV_PER_A_TO_NN,
                                  "force_per_strain_nn_by_stress": by_stress * EV_PER_A_TO_NN}
        got["seconds"] = got.get("seconds", 0) + round(time.time() - t0)
        store.write_text(json.dumps(got, indent=1))
        debye = float(np.hypot(*got.get("dipole_across_e_a", [np.nan] * 2)) * E_ANGSTROM_TO_DEBYE)
        k = got.get("stretch", {})
        print(f"{name}: rise {got['rise_a']:.4f} A, {got['volume_a3']:.1f} A3 a repeat, dipole across {debye:.3f} D a repeat, "
              f"stretch {k.get('force_per_strain_nn_by_stress', float('nan')):.1f} nN (by energy {k.get('force_per_strain_nn_by_energy', float('nan')):.1f}), "
              f"{time.time() - t0:.0f} s", flush=True)


def read_berry(out: Path):
    """``(value, modulus, direction)`` in C/m2 of one Berry-phase run: its part of the polarisation."""
    text = out.read_text(errors="replace")
    line = [l for l in text.splitlines() if "C/m^2" in l and "mod" in l][-1]
    value = float(line.split("=")[1].split()[0])
    modulus = float(line.split("mod")[1].split(")")[0])
    where = [l for l in text.splitlines() if "polarization direction is" in l][-1]
    direction = np.array([float(v) for v in where.split("(")[1].split(")")[0].split(",")])
    return value, modulus, direction


def across(got: dict):
    """The lone chain's dipole across itself (e A a repeat), under either name it has been stored."""
    if "dipole_across_e_a" in got:
        return np.array(got["dipole_across_e_a"])
    return np.array(got["dipole_e_a"][:2]) if "dipole_e_a" in got else None


def mesh_of(atoms):
    m = atoms.info["kmesh"]
    return [int(v) for v in (m.split() if isinstance(m, str) else np.atleast_1d(m))]


def berry_parts(args, atoms, work: Path, tag: str, near=None):
    """The Berry-phase parts of a cell's polarisation, one for each periodic edge (C/m2).

    Each is known up to a whole number of its quantum.  It is taken on the branch nearest zero, or
    nearest ``near``'s part when a strained cell follows an unstrained one."""
    a = atoms.copy()
    a.info.update(atoms.info)
    if any(s in NOT_NORM_CONSERVING for s in a.get_chemical_symbols()):
        raise ValueError("no norm-conserving file for an element of this cell: its Berry phase would be wrong")
    mesh = mesh_of(a)
    chain = a.info.get("triage_kind") == "chain"
    a.info["kmesh"] = " ".join(str(v) for v in mesh)
    scratch = work / f"scratch_{tag}"
    scf = work / f"{tag}_scf.out"
    if not (scf.exists() and "JOB DONE" in scf.read_text(errors="replace")):
        run(args.pw, args.np, pw_input(a, args.pseudo_berry, scratch, "scf", conv_thr=1e-10, files=PSEUDO_BERRY), work, f"{tag}_scf")
    parts = []
    for g in ((1, 2) if chain else (1, 2, 3)):
        along = list(mesh)
        along[g - 1] = max(args.string, mesh[g - 1])
        a.info["kmesh"] = " ".join(str(v) for v in along)
        text = pw_input(a, args.pseudo_berry, scratch, "nscf", conv_thr=1e-10, files=PSEUDO_BERRY)
        text = text.replace("  prefix = 'frame'", f"  prefix = 'frame'\n  lberry = .true.\n  gdir = {g}\n  nppstr = {along[g - 1]}")
        out = work / f"{tag}_{g}.out"
        if not (out.exists() and "JOB DONE" in out.read_text(errors="replace")):
            run(args.pw, args.np, text, work, f"{tag}_{g}")
        value, modulus, direction = read_berry(out)
        centre = 0.0 if near is None else near[len(parts)]["c_per_m2"]
        value -= modulus * np.round((value - centre) / modulus)
        parts.append({"gdir": g, "c_per_m2": value, "quantum_c_per_m2": modulus, "direction": direction.tolist()})
    return parts


def total(parts) -> np.ndarray:
    return sum(np.array(r["direction"]) * r["c_per_m2"] for r in parts)


def berry(args):
    """The polarisation of a relaxed cell by the Berry phase.

    The branch nearest zero is right while a part is below half its quantum (about 0.4 C/m2 for
    these cells); the lone chain's own moment over the cell's volume is printed beside it."""
    from ase.io import read

    src = Path(args.frames)
    for name in args.names.split(","):
        store = src / f"{name}.json"
        got = json.loads(store.read_text())
        a = read(src / f"{name}.relaxed.xyz")
        t0 = time.time()
        parts = berry_parts(args, a, src / "work" / name, "berry")
        p = total(parts)
        got["berry"] = {"parts": parts, "polarisation_c_per_m2": p.tolist(), "magnitude_c_per_m2": float(np.linalg.norm(p))}
        if a.info.get("triage_kind") != "chain":
            along = a.cell.array[2] / np.linalg.norm(a.cell.array[2])
            got["berry"].update(along_chain_c_per_m2=float(p @ along), across_chain_c_per_m2=float(np.linalg.norm(p - (p @ along) * along)))
        store.write_text(json.dumps(got, indent=1))
        note = ""
        if across(got) is not None:
            rigid = across(got) * got["repeats"] / a.get_volume() * E_PER_A2_TO_C_PER_M2
            note = f"; the density's moment over the cell {rigid[0]:+.4f} {rigid[1]:+.4f}"
        print(f"{name}: Berry polarisation {p[0]:+.4f} {p[1]:+.4f} {p[2]:+.4f} C/m2 (magnitude {np.linalg.norm(p):.4f}){note}, {time.time() - t0:.0f} s", flush=True)


def piezo(args):
    """A crystal strained along its own polarisation, the cell held across it: what it gives back.

    ``c``: the stress along the axis per strain (GPa).  ``e``: the polarisation along the axis per
    strain (C/m2), which is the stress the clamped crystal makes per field, pN/nm2 per V/um.  For
    dipoles that do not change, e is minus the polarisation (the same moment in more volume).
    ``d = e / c``: strain per field with the faces across the axis held, pm/V."""
    from ase.io import read

    src = Path(args.frames)
    for name in args.names.split(","):
        store = src / f"{name}.json"
        got = json.loads(store.read_text())
        a = read(src / f"{name}.relaxed.xyz")
        work = src / "work" / name
        # the axis is the polarisation across the chains, the direction a field through a film of lying
        # chains has; along the chain a directed chain (an amide's N to C) has a part of its own
        p0 = np.array(got["berry"]["polarisation_c_per_m2"])
        along = a.cell.array[2] / np.linalg.norm(a.cell.array[2])
        p_across = p0 - (p0 @ along) * along
        axis = p_across / np.linalg.norm(p_across)
        t0 = time.time()
        rows = []
        for strain in (-args.strain, args.strain):
            tag = f"strain_{'m' if strain < 0 else 'p'}"
            b = a.copy()
            b.info.update(a.info)
            b.info["kmesh"] = " ".join(str(v) for v in mesh_of(a))
            b.set_cell(a.cell.array @ (np.eye(3) + strain * np.outer(axis, axis)).T, scale_atoms=True)
            if not finished(work / f"{tag}.out"):
                run(args.pw, args.np, pw_input(b, args.pseudo, work / f"scratch_{tag}_relax", "relax"), work, tag)
            if not finished(work / f"{tag}.out"):
                print(f"{name}: {tag} did not finish", flush=True)
                break
            c, energy, stress = read_final(work / f"{tag}.out")
            c.info.update(b.info)
            parts = berry_parts(args, c, work, f"{tag}_berry", near=got["berry"]["parts"])
            rows.append({"strain": strain, "energy_ev": energy, "stress_along_ev_a3": float(axis @ stress @ axis),
                         "polarisation_along_c_per_m2": float(total(parts) @ axis), "parts": parts})
        if len(rows) != 2:
            continue
        c_along = (rows[1]["stress_along_ev_a3"] - rows[0]["stress_along_ev_a3"]) / (2 * args.strain) * 160.21766
        e_along = (rows[1]["polarisation_along_c_per_m2"] - rows[0]["polarisation_along_c_per_m2"]) / (2 * args.strain)
        got["piezo"] = {"axis": axis.tolist(), "strain": args.strain, "rows": rows, "c_along_gpa": c_along, "e_along_c_per_m2": e_along,
                        "d_clamped_pm_per_v": e_along / c_along * 1000.0, "polarisation_across_c_per_m2": float(np.linalg.norm(p_across)),
                        "polarisation_along_chain_c_per_m2": float(p0 @ along)}
        store.write_text(json.dumps(got, indent=1))
        print(f"{name}: P across the chains {np.linalg.norm(p_across):.4f} C/m2, e {e_along:+.4f} C/m2, c {c_along:.1f} GPa, d (cell held across the axis) {e_along / c_along * 1000:+.1f} pm/V, "
              f"{time.time() - t0:.0f} s", flush=True)


def table(args):
    """Every record of a folder in one file, and the relaxed cells beside it."""
    from ase.io import read, write

    src = Path(args.frames)
    rows = [json.loads(p.read_text()) for p in sorted(src.glob("*.json")) if not p.name.endswith(".path.json")]
    paths = [json.loads(p.read_text()) for p in sorted(src.glob("*.path.json"))]
    out = {"what": "chains alone in a box and one chain to a cell, PBE-D3(BJ) at 90/360 Ry (examples/chain_triage.py); "
                   "polarisation by the Berry phase with norm-conserving files on the same geometry",
           "units": {"rise_a": "A a repeat along the chain", "volume_a3": "A3 a repeat (a lone chain's is its box)",
                     "dipole_across_e_a": "e A a repeat, across the chain, from the charge density",
                     "stretch": "nN per unit strain of the chain", "berry": "C/m2", "piezo": "c GPa, e C/m2 (pN/nm2 per V/um), d pm/V"},
           "rows": rows, "paths": paths}
    for q in paths:
        print(f"{q['name']:20s} turning every second chain of {q['of']}: top {q['top_kcal_per_turning_repeat']:.2f} kcal/mol a repeat that turns, "
              f"end {q['kcal_per_turning_repeat'][-1]:+.2f} (DFT on the potential's path)")
    if args.cells:
        cells = []
        for r in rows:
            c = read(src / f"{r['name']}.relaxed.xyz")
            c.info = {"name": r["name"], "kind": r["kind"], "repeats": r["repeats"], "energy_ev": r["energy_ev"]}
            c.calc = None
            cells.append(c)
        write(args.cells, cells, format="extxyz")
    for r in rows:
        mu = across(r) if across(r) is not None else [np.nan] * 2
        k = r.get("stretch", {})
        y = dict(r.get("berry", {}))
        if "polarisation_c_per_m2" in y and r["kind"] != "chain":  # a record made before the two parts were stored
            pol, along = np.array(y["polarisation_c_per_m2"]), np.array(r["cell"][2]) / np.linalg.norm(r["cell"][2])
            y.update(along_chain_c_per_m2=float(pol @ along), across_chain_c_per_m2=float(np.linalg.norm(pol - (pol @ along) * along)))
            r["berry"] = y
        z = r.get("piezo", {})
        print(f"{r['name']:14s} {r['kind']:7s} rise {r['rise_a']:.4f} A  {r['volume_a3']:6.1f} A3  E {r['energy_ev'] / r['repeats']:.5f} eV a repeat  "
              f"dipole across {np.hypot(*mu) * E_ANGSTROM_TO_DEBYE:.3f} D  stretch {k.get('force_per_strain_nn_by_stress', float('nan')):.1f} nN  "
              f"P across {y.get('across_chain_c_per_m2', y.get('magnitude_c_per_m2', float('nan'))):.4f}  "
              f"e {z.get('e_along_c_per_m2', float('nan')):+.4f}  c {z.get('c_along_gpa', float('nan')):.1f}  d {z.get('d_clamped_pm_per_v', float('nan')):+.1f}")
    by_name = {r["name"]: r for r in rows}
    for r in rows:
        if r.get("opposed_of") in by_name:
            polar = by_name[r["opposed_of"]]
            r["above_polar_kcal_per_repeat"] = (r["energy_ev"] / r["repeats"] - polar["energy_ev"] / polar["repeats"]) * KCAL
            print(f"{r['name']:14s} dipoles opposed: {r['above_polar_kcal_per_repeat']:+.3f} kcal/mol a repeat above {polar['name']} "
                  f"({r['above_polar_kcal_per_repeat'] / KCAL * 1000:+.2f} meV), {r['volume_a3']:.1f} against {polar['volume_a3']:.1f} A3")
    if args.out:
        args.out.write_text(json.dumps(out, indent=1) + "\n")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="command", required=True)
    b = sub.add_parser("build")
    b.add_argument("--chains", required=True)
    b.add_argument("--out", type=Path, required=True)
    b.add_argument("--model", type=Path, default=Path("/home/user/rowansci/egret-public/compiled_models/EGRET_1.model"))
    b.add_argument("--box", type=float, default=14.0)
    b.set_defaults(run=build)
    d = sub.add_parser("dft")
    d.add_argument("--frames", type=Path, required=True)
    d.add_argument("--names", default="")
    d.add_argument("--pw", default="/opt/mm/qe/bin/pw.x")
    d.add_argument("--pseudo", default="/opt/pseudo")
    d.add_argument("--np", type=int, default=4)
    d.add_argument("--strain", type=float, default=0.01)
    d.add_argument("--no-stretch", action="store_true")
    d.set_defaults(run=dft)
    k = sub.add_parser("pack")
    k.add_argument("--frames", type=Path, required=True)
    k.add_argument("--names", required=True)
    k.add_argument("--model", type=Path, default=Path("/home/user/rowansci/egret-public/compiled_models/EGRET_1.model"))
    k.add_argument("--starts", type=int, default=24)
    k.add_argument("--steps", type=int, default=300)
    k.add_argument("--keep", type=int, default=2)
    k.add_argument("--seed", type=int, default=11)
    k.set_defaults(run=pack)
    y = sub.add_parser("berry")
    y.add_argument("--frames", type=Path, required=True)
    y.add_argument("--names", required=True)
    y.add_argument("--pw", default="/opt/mm/qe/bin/pw.x")
    y.add_argument("--pseudo", default="/opt/pseudo")
    y.add_argument("--pseudo-berry", default="/opt/pseudo_nc", help="norm-conserving files, for the Berry phase")
    y.add_argument("--np", type=int, default=4)
    y.add_argument("--string", type=int, default=6)
    y.set_defaults(run=berry)
    z = sub.add_parser("piezo")
    z.add_argument("--frames", type=Path, required=True)
    z.add_argument("--names", required=True)
    z.add_argument("--pw", default="/opt/mm/qe/bin/pw.x")
    z.add_argument("--pseudo", default="/opt/pseudo")
    z.add_argument("--pseudo-berry", default="/opt/pseudo_nc", help="norm-conserving files, for the Berry phase")
    z.add_argument("--np", type=int, default=4)
    z.add_argument("--string", type=int, default=6)
    z.add_argument("--strain", type=float, default=0.01)
    z.set_defaults(run=piezo)
    o = sub.add_parser("opposed")
    o.add_argument("--frames", type=Path, required=True)
    o.add_argument("--names", required=True)
    o.add_argument("--model", type=Path, default=Path("/home/user/rowansci/egret-public/compiled_models/EGRET_1.model"))
    o.add_argument("--steps", type=int, default=200)
    o.set_defaults(run=opposed)
    u = sub.add_parser("turn")
    u.add_argument("--frames", type=Path, required=True)
    u.add_argument("--names", required=True, help="NAME:EDGE, the polar cell and the edge (1 or 2) it is doubled along")
    u.add_argument("--model", type=Path, default=Path("/home/user/rowansci/egret-public/compiled_models/EGRET_1.model"))
    u.add_argument("--images", type=int, default=5)
    u.add_argument("--fmax", type=float, default=0.15)
    u.add_argument("--steps", type=int, default=80)
    u.set_defaults(run=turn)
    v = sub.add_parser("wall")
    v.add_argument("--frames", type=Path, required=True)
    v.add_argument("--names", required=True)
    v.add_argument("--model", type=Path, default=Path("/home/user/rowansci/egret-public/compiled_models/EGRET_1.model"))
    v.add_argument("--images", type=int, default=5)
    v.add_argument("--fmax", type=float, default=0.15)
    v.add_argument("--steps", type=int, default=80)
    v.set_defaults(run=wall)
    x = sub.add_parser("nucleus")
    x.add_argument("--frames", type=Path, required=True)
    x.add_argument("--names", required=True)
    x.add_argument("--model", type=Path, default=Path("/home/user/rowansci/egret-public/compiled_models/EGRET_1.model"))
    x.add_argument("--methyl", action="store_true", help="make the N-H that holds the turning chain an N-CH3 first")
    x.add_argument("--images", type=int, default=5)
    x.add_argument("--fmax", type=float, default=0.15)
    x.add_argument("--steps", type=int, default=80)
    x.add_argument("--settle", type=int, default=150)
    x.add_argument("--gap", type=float, default=1.6, help="A of room made along the row for the methyl before the cell is settled")
    x.set_defaults(run=nucleus)
    w = sub.add_parser("points")
    w.add_argument("--frames", type=Path, required=True)
    w.add_argument("--names", required=True)
    w.add_argument("--pw", default="/opt/mm/qe/bin/pw.x")
    w.add_argument("--pseudo", default="/opt/pseudo")
    w.add_argument("--np", type=int, default=4)
    w.set_defaults(run=points)
    t = sub.add_parser("table")
    t.add_argument("--frames", type=Path, required=True)
    t.add_argument("--out", type=Path, default=None)
    t.add_argument("--cells", type=Path, default=None, help="write every relaxed cell to this extxyz file")
    t.set_defaults(run=table)
    args = p.parse_args()
    args.run(args)


if __name__ == "__main__":
    main()
