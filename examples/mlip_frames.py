"""DFT frames to train our own interatomic potential on, in the anchor's recipe.

`examples/mlip_second_opinion.py` showed that no pretrained potential resolves the one margin
this project holds from DFT. A potential of our own needs frames: cells with their DFT energy,
the force on every atom and the stress. The project held none (the stored relaxations keep
their end points only), so this script makes them.

The labels are computed exactly as the anchor was (`deliverables/polar_margin_dft/`):
Quantum ESPRESSO 7.5, PBE with two-body D3(BJ), the SSSP 1.3.0 PBE precision files, 90/360 Ry,
the anchor's k-point density. `check` reruns the anchor's own beta-PVDF input and must return
its energy to the last digit before any label is trusted.

Three steps; `label` is resumable:

    python examples/mlip_frames.py make    --frames DIR [--set pvdf | near | units]
    python examples/mlip_frames.py label   --frames DIR --pw PW.X --pseudo DIR [--np 4]
    python examples/mlip_frames.py collect --frames DIR --out docs/predictions/mlip/NAME.extxyz

and `cutoff`, which reruns labelled frames at a higher plane-wave cutoff, and `aim`, which
moves the units set's unlabelled frames to where the DFT stress on each cell's relaxed frame
points (label the relaxed frames first: `label --only units`), and `relax`, which relaxes a
cell in DFT itself and keeps every step as a frame.

The PVDF set, 76 frames from the four DFT-relaxed cells (helix alpha, polar beta, the two
antipolar betas):

- `end`: the relaxed cells themselves (4);
- `strain`: alpha, polar beta and antipolar beta under single strains of 2% along each axis,
  both signs, and a lateral shear (24). These carry the lattice's stiffness;
- `rattle`: the two-monomer-per-chain cells with a random strain up to 2% and every atom moved
  by a Gaussian of 0.05 A (18). These carry the stiff modes: bonds and angles;
- `md`: snapshots of molecular dynamics at 500 K in a pretrained potential, with a random
  strain (30). The pretrained potential is used here as a sampler only: it decides which
  geometries are labelled, not what their labels are. These carry the soft modes: torsions,
  the chains' rotation and sliding.

What the PVDF set does not hold, and so what a potential trained on it cannot be asked: any
unit but VDF, and any geometry on the way between all-trans and the helix (a barrier).

The near set, 54 frames, is PVDF again at the bottom of each well (strains of 1% and 4%, both
lateral axes at once, shakes of 0.01 A): what the first fits lacked to get the lattice right.

The units set, 168 frames, is for the first of those: TrFE, CFE, CTFE, acrylonitrile, the
reversed VDF unit and TFE put into the DFT-relaxed PVDF cells (a unit's groups replaced along
their own bonds), then relaxed and shaken in the pretrained sampler. Twelve compositions with
two monomers per chain (24 to 27 atoms: each unit alone, and pairs facing each other on
neighbouring chains) and six with four monomers per chain (48 to 51 atoms: the leading cells'
units together), each all-trans and helix.

Needs ase; `make` also needs torch and mace-torch for the sampler.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RELAXED = ROOT / "deliverables" / "polar_margin_dft" / "relaxed"
ANCHOR_INPUT = ROOT / "deliverables" / "polar_margin_dft" / "inputs" / "polar" / "scf_k1.in"
ANCHOR_ENERGY_RY = -275.80701219  # raw_stages.json, cells.polar.stages.scf_k1
ANCHOR_D3_RY = -0.05243146

PSEUDO = {
    "H": ("H_ONCV_PBE-1.0.oncvpsp.upf", 1.00794),
    "C": ("C.pbe-n-kjpaw_psl.1.0.0.UPF", 12.0107),
    "F": ("F.oncvpsp.upf", 18.998403163),
    "N": ("N.oncvpsp.upf", 14.0067),
    "Cl": ("Cl.pbe-n-rrkjus_psl.1.0.0.UPF", 35.453),
}
# name: (file, monomers in the cell, the anchor's k-point mesh)
CELLS = {
    "alpha": ("pvdf_alpha_pbe_d3bj.xyz", 4, (4, 2, 4)),
    "polar": ("pvdf_polar_pbe_d3bj.xyz", 2, (4, 2, 8)),
    "anti": ("pvdf_anti_free_pbe_d3bj.xyz", 2, (4, 2, 8)),
    "anti90": ("pvdf_anti90_pbe_d3bj.xyz", 2, (4, 2, 8)),
}
STRAIN = 0.02
RATTLE_A = 0.05
MD_K = 500.0
SEED = 20261116


def cell(name, repeat=(1, 1, 1)):
    from ase.io import read

    f, n, mesh = CELLS[name]
    a = read(RELAXED / f) * repeat
    a.info = {
        "source": name,
        "n_monomers": n * int(np.prod(repeat)),
        "kmesh": " ".join(str(max(1, m // r)) for m, r in zip(mesh, repeat)),
    }
    return a


def strained(atoms, e):
    """The cell and every atom moved by the small strain e (Voigt: xx yy zz yz xz xy)."""
    a = atoms.copy()
    eps = np.array([[e[0], e[5] / 2, e[4] / 2], [e[5] / 2, e[1], e[3] / 2], [e[4] / 2, e[3] / 2, e[2]]])
    a.set_cell(a.cell.array @ (np.eye(3) + eps), scale_atoms=True)
    return a


def md_snapshots(atoms, n, rng, model, temperature, every, threads):
    """Snapshots of Langevin dynamics at a fixed cell in a pretrained potential: a sampler."""
    import torch
    from ase import units
    from ase.md.langevin import Langevin
    from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
    from mace.calculators import mace_mp

    torch.set_num_threads(threads)
    a = atoms.copy()
    a.calc = mace_mp(model=model, dispersion=False, default_dtype="float32", device="cpu")
    MaxwellBoltzmannDistribution(a, temperature_K=temperature, rng=rng)
    dyn = Langevin(a, 1.0 * units.fs, temperature_K=temperature, friction=0.02 / units.fs, rng=rng)
    out = []
    for _ in range(n):
        dyn.run(every)
        s = a.copy()
        s.info = dict(atoms.info)
        out.append(s)
    return out


def pvdf_set(model, threads):
    rng = np.random.default_rng(SEED)
    frames = []

    def add(a, family, note):
        a = a.copy()
        a.calc = None
        a.info = {**a.info, "family": family, "note": note}
        frames.append(a)

    for name in CELLS:
        add(cell(name), "end", "the DFT-relaxed cell")
    for name in ("alpha", "polar", "anti"):
        base = cell(name)
        for axis, label in ((0, "xx"), (1, "yy"), (2, "zz"), (5, "xy")):
            for sign in (+1, -1):
                e = np.zeros(6)
                e[axis] = sign * STRAIN
                add(strained(base, e), "strain", f"{label} {sign * STRAIN:+.3f}")
    doubled = {"alpha": cell("alpha"), "polar": cell("polar", (1, 1, 2)), "anti": cell("anti", (1, 1, 2))}
    for name, base in doubled.items():
        for i in range(6):
            a = strained(base, rng.uniform(-STRAIN, STRAIN, 6))
            a.positions += rng.normal(0.0, RATTLE_A, a.positions.shape)
            add(a, "rattle", f"random strain to {STRAIN}, sigma {RATTLE_A} A, draw {i}")
    for name, base in doubled.items():
        for i, s in enumerate(md_snapshots(base, 10, rng, model, MD_K, 100, threads)):
            add(strained(s, rng.uniform(-STRAIN, STRAIN, 6)), "md",
                f"{MD_K:.0f} K in {model} (sampler only), 100 fs apart, random strain to {STRAIN}, snapshot {i}")
    return frames


def near_set(model, threads):
    """PVDF frames at the bottom of each well, where the margins and the lattice are decided.

    The first fits (docs/MLIP_SECOND_OPINION.md 3.3) had the atoms right and the lattice across
    the chains wrong: at the 2% strains of the PVDF set the energy is a quarter of a meV per
    atom above the minimum, under the fits' error. These are strains of 1% and 4% across the
    chains, 1% along them, both lateral axes at once, and shakes of 0.01 A.
    """
    rng = np.random.default_rng(SEED + 2)
    frames = []

    def add(a, family, note):
        a = a.copy()
        a.calc = None
        a.info = {**a.info, "family": family, "note": note}
        frames.append(a)

    for name in ("alpha", "polar", "anti"):
        base = cell(name)
        singles = [(0, "xx", e) for e in (-0.04, -0.01, 0.01, 0.04)] + [(1, "yy", e) for e in (-0.04, -0.01, 0.01, 0.04)]
        singles += [(2, "zz", e) for e in (-0.01, 0.01)]
        for axis, label, e in singles:
            strain = np.zeros(6)
            strain[axis] = e
            add(strained(base, strain), "n-strain", f"{label} {e:+.3f}")
        for e in (-0.04, -0.02, 0.02, 0.04):
            add(strained(base, np.array([e, e, 0.0, 0.0, 0.0, 0.0])), "n-strain", f"xx and yy {e:+.3f}")
        for i in range(4):
            a = strained(base, rng.uniform(-0.01, 0.01, 6))
            a.positions += rng.normal(0.0, 0.01, a.positions.shape)
            add(a, "n-rattle", f"random strain to 0.01, sigma 0.01 A, draw {i}")
    return frames


# What a VDF unit's two carbons carry after it is turned into another unit: (the CH2 carbon's
# two groups, the CF2 carbon's two groups). The letters are examples/candidate_screen.py's.
UNIT_GROUPS = {
    "V": (("H", "H"), ("F", "F")),
    "R": (("H", "F"), ("F", "F")),    # TrFE, CHF-CF2
    "C": (("H", "H"), ("F", "Cl")),   # CFE, CH2-CFCl
    "K": (("F", "F"), ("F", "Cl")),   # CTFE, CF2-CFCl
    "A": (("H", "H"), ("H", "CN")),   # acrylonitrile, CH2-CH(CN)
    "I": (("F", "F"), ("H", "H")),    # a reversed VDF unit, CF2-CH2 (the head-to-head defect)
    "F": (("F", "F"), ("F", "F")),    # TFE, CF2-CF2
}
BOND_A = {"H": 1.09, "F": 1.36, "Cl": 1.77, "CN": 1.47}
NITRILE_A = 1.16
# name: (the PVDF cell, its repeat along the chain, monomers per chain). Two chains each.
UNIT_BASES = {"TT2": ("polar", (1, 1, 2), 2), "TG2": ("alpha", (1, 1, 1), 2),
              "TT4": ("polar", (1, 1, 4), 4), "TG4": ("alpha", (1, 1, 2), 4)}
# One unit at a time in the smallest cell that holds both conformations (two monomers per
# chain), then pairs of units facing each other on neighbouring chains.
UNIT_SMALL = [("VR", "VV"), ("VR", "VR"), ("VC", "VV"), ("VK", "VV"), ("VA", "VV"), ("VI", "VV"), ("VF", "VV"),
              ("VR", "VC"), ("VR", "VA"), ("VI", "VA"), ("VR", "VK"), ("VI", "VK")]
# Four monomers per chain, after the cells the ranking is led by (docs/CANDIDATE_MATERIALS.md 5.6).
UNIT_LARGE = [
    (("VVRV", "VRVV"), "a TrFE matrix, 25 mol%"),
    (("VRVA", "IVRV"), "the nitrile with a reversed unit in a TrFE matrix (the leader's units)"),
    (("VRVC", "VVRV"), "CFE in a TrFE matrix (the film's units)"),
    (("VRVK", "IVRV"), "CTFE with a reversed unit in a TrFE matrix"),
    (("VIVV", "VVIV"), "PVDF with reversed units"),
    (("VAVF", "VRVV"), "the nitrile with TFE"),
]


def pvdf_chains(atoms):
    """Each chain's monomers along the chain: (CH2 carbon, its two H, CF2 carbon, its two F)."""
    from ase.neighborlist import natural_cutoffs, neighbor_list

    i, j = neighbor_list("ij", atoms, natural_cutoffs(atoms, mult=1.1))
    sym = atoms.get_chemical_symbols()
    bonded = {k: [int(b) for a, b in zip(i, j) if a == k] for k in range(len(atoms))}
    carbons = {k: [b for b in bonded[k] if sym[b] == "C"] for k in range(len(atoms)) if sym[k] == "C"}
    groups = {k: [b for b in bonded[k] if sym[b] != "C"] for k in carbons}
    if any(len(v) != 2 for v in groups.values()):
        raise ValueError("not a PVDF cell: a carbon without two groups")
    chains, seen = [], set()
    for start in sorted(k for k in carbons if sym[groups[k][0]] == "H"):
        if start in seen:
            continue
        if len(set(carbons[start])) != 2:
            raise ValueError("a chain of one monomer per cell: repeat the cell along the chain first")
        walk = [start, carbons[start][0]]  # either way along the chain: a sequence reads both ways
        while True:
            nxt = [b for b in carbons[walk[-1]] if b != walk[-2]][0]
            if nxt == start:
                break
            walk.append(nxt)
        seen.update(walk)
        if len(walk) % 2:
            raise ValueError("a chain with an odd number of carbons")
        chains.append([(walk[m], groups[walk[m]], walk[m + 1], groups[walk[m + 1]]) for m in range(0, len(walk), 2)])
    if len(seen) != len(carbons):
        raise ValueError("carbons left outside the chains")
    return chains


def with_units(atoms, sequences, rng):
    """The PVDF cell with each chain's monomers turned into the units of its sequence.

    A group is replaced along its own bond at the new bond's length; which of a carbon's two
    groups is the one replaced is drawn at random, so both hands of a chiral unit occur. The
    result is a start for a relaxation, not a geometry: the cell is still PVDF's.
    """
    from ase import Atom

    a = atoms.copy()
    chains = pvdf_chains(a)
    if len(chains) != len(sequences) or any(len(c) != len(s) for c, s in zip(chains, sequences)):
        raise ValueError(f"{[len(c) for c in chains]} monomers per chain against the sequences {sequences}")
    extra = []
    for chain, letters in zip(chains, sequences):
        for (c1, g1, c2, g2), letter in zip(chain, letters):
            for carbon, old, new in ((c1, g1, UNIT_GROUPS[letter][0]), (c2, g2, UNIT_GROUPS[letter][1])):
                new = list(new) if new[0] == new[1] else [new[k] for k in rng.permutation(2)]
                for atom, group in zip(old, new):
                    bond = a.get_distance(carbon, atom, mic=True, vector=True)
                    bond /= np.linalg.norm(bond)
                    a.positions[atom] = a.positions[carbon] + BOND_A[group] * bond
                    a[atom].symbol = "C" if group == "CN" else group
                    if group == "CN":
                        extra.append(Atom("N", a.positions[atom] + NITRILE_A * bond))
    for atom in extra:
        a.append(atom)
    return a


def sampler_relaxed(atoms, model, threads, fmax, steps):
    """Atoms and cell relaxed in the pretrained potential with D3: a sampler's geometry."""
    import torch
    from ase.filters import FrechetCellFilter
    from ase.optimize import BFGS
    from mace.calculators import mace_mp

    torch.set_num_threads(threads)
    a = atoms.copy()
    a.calc = mace_mp(model=model, dispersion=True, default_dtype="float64", device="cpu")
    # atoms first at PVDF's cell (a chlorine or a nitrile starts in contact), then the cell with them
    BFGS(a, logfile=None, maxstep=0.1).run(fmax=5 * fmax, steps=steps // 2)
    opt = BFGS(FrechetCellFilter(a), logfile=None, maxstep=0.1)
    opt.run(fmax=fmax, steps=steps)
    out = a.copy()
    out.info = {**atoms.info, "sampler_volume_change_percent": float(100 * (a.get_volume() / atoms.get_volume() - 1)),
                "sampler_largest_force": float(np.abs(a.get_forces()).max())}
    return out


def units_cells():
    """The cells of the units set: (PVDF base, each chain's sequence, what it is after)."""
    small = [(f"{conf}2", seqs, "") for seqs in UNIT_SMALL for conf in ("TT", "TG")]
    large = [(f"{conf}4", seqs, what) for seqs, what in UNIT_LARGE for conf in ("TT", "TG")]
    return small + large


def units_cell(index, model, threads):
    """One cell of the units set: relaxed in the sampler, a shaken copy, dynamics snapshots."""
    base_name, seqs, what = units_cells()[index]
    rng = np.random.default_rng([SEED + 1, index])  # a cell's frames do not depend on the others
    source, repeat, per_chain = UNIT_BASES[base_name]
    base = cell(source, repeat)
    tag = f"{base_name} {'|'.join(seqs)}"
    start = with_units(base, seqs, rng)
    start.info = {**base.info, "source": tag, "sequences": "|".join(seqs), "conformation": base_name[:2]}
    relaxed = sampler_relaxed(start, model, threads, 0.1, 200)
    frames = []

    def add(a, family, note):
        a = a.copy()
        a.calc = None
        a.info = {**relaxed.info, "family": family, "note": note}
        frames.append(a)

    add(relaxed, "units", f"{tag}{': ' + what if what else ''}; from the PVDF cell, relaxed in {model} with D3 "
                          f"(sampler only), volume {relaxed.info['sampler_volume_change_percent']:+.1f}%")
    r = strained(relaxed, rng.uniform(-STRAIN, STRAIN, 6))
    r.positions += rng.normal(0.0, RATTLE_A, r.positions.shape)
    add(r, "u-rattle", f"{tag}; random strain to {STRAIN}, sigma {RATTLE_A} A")
    for k, snap in enumerate(md_snapshots(relaxed, 3 if per_chain == 2 else 2, rng, model, MD_K, 100, threads)):
        add(strained(snap, rng.uniform(-STRAIN, STRAIN, 6)), "u-md",
            f"{tag}; {MD_K:.0f} K in {model} (sampler only), 100 fs apart, random strain to {STRAIN}, snapshot {k}")
    print(f"{tag}: {len(start)} atoms, sampler volume {relaxed.info['sampler_volume_change_percent']:+.1f}%, "
          f"largest force left {relaxed.info['sampler_largest_force']:.2f} eV/A", flush=True)
    return frames


def make(args):
    from ase.io import write

    args.frames.mkdir(parents=True, exist_ok=True)
    if args.set == "units":  # cell by cell, so a run can stop and the next one go on
        began, cells = time.time(), units_cells()
        for index in range(len(cells)):
            if (args.frames / f"units_{index:02d}_0.start.xyz").exists():
                continue
            if args.minutes and time.time() - began > 60 * args.minutes:
                break
            made = units_cell(index, args.sampler, args.threads)
            for k, a in reversed(list(enumerate(made))):  # the first frame last: it marks the cell as done
                a.info["frame"] = f"units_{index:02d}_{k}"
                write(args.frames / f"{a.info['frame']}.start.xyz", a, format="extxyz")
        done = sum((args.frames / f"units_{i:02d}_0.start.xyz").exists() for i in range(len(cells)))
        print(f"{done} of {len(cells)} cells made, {len(list(args.frames.glob('units_*.start.xyz')))} frames in {args.frames}")
        return
    frames = {"pvdf": pvdf_set, "near": near_set}[args.set](args.sampler, args.threads)
    for i, a in enumerate(frames):
        a.info["frame"] = f"{args.set}_{i:03d}"
        write(args.frames / f"{a.info['frame']}.start.xyz", a, format="extxyz")
    counts = {}
    for a in frames:
        counts[a.info["family"]] = counts.get(a.info["family"], 0) + 1
    print(f"{len(frames)} frames in {args.frames}: " + ", ".join(f"{k} {v}" for k, v in counts.items())
          + f"; atoms {min(len(a) for a in frames)} to {max(len(a) for a in frames)}")


def pw_input(atoms, pseudo_dir, outdir, conv_thr, ecut=(90, 360)):
    """The anchor's input (`deliverables/polar_margin_dft/inputs/*/scf_k1.in`) for this cell."""
    species = sorted(set(atoms.get_chemical_symbols()), key=lambda s: -atoms.get_chemical_symbols().count(s))
    lines = [
        "&CONTROL", "  calculation = 'scf'", "  prefix = 'frame'", f"  pseudo_dir = '{pseudo_dir}'",
        f"  outdir = '{outdir}'", "  tstress = .true.", "  tprnfor = .true.", "  disk_io = 'low'", "/", "",
        "&SYSTEM", "  ibrav = 0", f"  nat = {len(atoms)}", f"  ntyp = {len(species)}", "  input_dft = 'PBE'",
        "  occupations = 'fixed'", f"  ecutwfc = {ecut[0]}", f"  ecutrho = {ecut[1]}", "  vdw_corr = 'grimme-d3'",
        "  dftd3_version = 4", "  dftd3_threebody = .false.", "/", "",
        "&ELECTRONS", f"  conv_thr = {conv_thr:.1e}".replace("e", "d"), "  electron_maxstep = 200",
        "  mixing_mode = 'plain'", "  mixing_beta = 0.3", "  startingwfc = 'atomic+random'", "/", "",
        "ATOMIC_SPECIES",
    ]
    lines += [f"{s} {PSEUDO[s][1]} {PSEUDO[s][0]}" for s in species]
    lines += ["", "CELL_PARAMETERS angstrom"] + [" ".join(f"{v:.10f}" for v in row) for row in atoms.cell.array]
    lines += ["", "ATOMIC_POSITIONS angstrom"]
    lines += [f"{s} {p[0]:.10f} {p[1]:.10f} {p[2]:.10f}" for s, p in zip(atoms.get_chemical_symbols(), atoms.positions)]
    mesh = atoms.info["kmesh"]  # a string when made, an array once read back from the file
    mesh = mesh.split() if isinstance(mesh, str) else np.atleast_1d(mesh)
    lines += ["", "K_POINTS automatic", " ".join(str(int(m)) for m in mesh) + " 0 0 0", ""]
    return "\n".join(lines)


POOLS = 0  # k-point pools of a pw.x run; 0: one a process (--pools)


def run_pw(pw, np_, text, workdir, name):
    """Each k-point pool holds the whole cell's arrays: a chain in a wide box at one pool a process
    asks for more memory than the cell with four times its atoms, and fewer pools is the way out."""
    workdir.mkdir(parents=True, exist_ok=True)
    (workdir / f"{name}.in").write_text(text)
    mpirun = Path(pw).with_name("mpirun")
    cmd = [str(mpirun if mpirun.exists() else "mpirun"), "--allow-run-as-root", "--oversubscribe", "-np", str(np_),
           str(pw), "-nk", str(POOLS or np_), "-in", f"{name}.in"]
    env = {**os.environ, "OMP_NUM_THREADS": "1", "PATH": f"{Path(pw).parent}:{os.environ.get('PATH', '/usr/bin:/bin')}"}
    t0 = time.time()
    with open(workdir / f"{name}.out", "w") as out:
        subprocess.run(cmd, cwd=workdir, stdout=out, stderr=subprocess.STDOUT, check=False, env=env)
    return workdir / f"{name}.out", time.time() - t0


def read_pw(path):
    """Energy (eV), forces (eV/A), stress (eV/A3, ASE's sign) and the D3 energy (Ry) of a finished run."""
    from ase.io import read

    text = path.read_text(errors="replace")
    if "JOB DONE" not in text or "convergence has been achieved" not in text:
        return None
    a = read(path, format="espresso-out")
    d3 = [float(l.split("=")[1].split()[0]) for l in text.splitlines() if "DFT-D3 Dispersion         =" in l]
    ry = [float(l.split("=")[1].split()[0]) for l in text.splitlines() if l.startswith("!")]
    return {
        "energy": float(a.get_potential_energy()), "forces": a.get_forces(), "stress": a.get_stress(voigt=False),
        "energy_ry": ry[-1], "d3_ry": d3[-1] if d3 else None,
        "kpoints": int([l for l in text.splitlines() if "number of k points=" in l][0].split("=")[1].split()[0]),
    }


def check(args):
    """The anchor's own beta-PVDF input must return the anchor's energy."""
    work = args.frames / "check"
    text = ANCHOR_INPUT.read_text()
    text = "\n".join(
        f"  pseudo_dir = '{args.pseudo}'" if l.strip().startswith("pseudo_dir") else
        f"  outdir = '{work / 'scratch'}'" if l.strip().startswith("outdir") else l
        for l in text.splitlines()) + "\n"
    out, seconds = run_pw(args.pw, args.np, text, work, "anchor_beta")
    r = read_pw(out)
    shutil.rmtree(work / "scratch", ignore_errors=True)
    ok = r is not None and abs(r["energy_ry"] - ANCHOR_ENERGY_RY) < 2e-8 and abs(r["d3_ry"] - ANCHOR_D3_RY) < 2e-8
    print(f"anchor beta-PVDF: {r and r['energy_ry']} Ry (the anchor: {ANCHOR_ENERGY_RY}), D3 {r and r['d3_ry']} "
          f"(the anchor: {ANCHOR_D3_RY}), {seconds:.0f} s on {args.np} processes: {'reproduced' if ok else 'NOT REPRODUCED'}")
    (work / "check.json").write_text(json.dumps({"reproduced": ok, "energy_ry": r and r["energy_ry"], "seconds": seconds}))
    return ok


def label(args):
    from ase.io import read, write

    if not (args.frames / "check" / "check.json").exists() or not json.loads(
            (args.frames / "check" / "check.json").read_text())["reproduced"]:
        if not check(args):
            raise SystemExit("the anchor's energy is not reproduced: no label is computed")
    starts = sorted(args.frames.glob("*.start.xyz"))
    if args.only:
        starts = [s for s in starts if any(k in read(s).info["family"] for k in args.only.split(","))]
    if args.names:
        starts = [s for s in starts if s.name.replace(".start.xyz", "") in args.names.split(",")]
    began = time.time()
    for n, start in enumerate(starts):
        name = start.name.replace(".start.xyz", "")
        done = args.frames / f"{name}.dft.xyz"
        if done.exists():
            continue
        if args.minutes and time.time() - began > 60 * args.minutes:
            left = sum(not (args.frames / s.name.replace(".start.xyz", ".dft.xyz")).exists() for s in starts)
            print(f"the {args.minutes:g} minutes are used: {left} frames left to label", flush=True)
            return
        a = read(start)
        work = args.frames / "work" / name
        shutil.rmtree(work, ignore_errors=True)  # what an interrupted run left
        out, seconds = run_pw(args.pw, args.np, pw_input(a, args.pseudo, work / "scratch", args.conv_thr), work, name)
        r = read_pw(out)
        shutil.rmtree(work / "scratch", ignore_errors=True)
        if r is None:
            print(f"{name}: no converged result ({seconds:.0f} s), see {out}", flush=True)
            continue
        a.info.update({"REF_energy": r["energy"], "REF_stress": r["stress"], "dft_energy_ry": r["energy_ry"],
                       "dft_d3_ry": r["d3_ry"], "dft_kpoints": r["kpoints"], "dft_seconds": round(seconds, 1),
                       "dft_conv_thr_ry": args.conv_thr})
        a.arrays["REF_forces"] = r["forces"]
        write(done, a, format="extxyz")
        shutil.copy(out, args.frames / f"{name}.pw.out")
        shutil.rmtree(work, ignore_errors=True)
        print(f"{name} ({n + 1}/{len(starts)}) {a.info['family']:7s} {len(a)} atoms, {r['kpoints']} k: "
              f"{r['energy'] / a.info['n_monomers']:.5f} eV per monomer, largest force "
              f"{np.abs(r['forces']).max():.2f} eV/A, {seconds:.0f} s", flush=True)


def dft_relax(args):
    """Relax cells in DFT itself, atoms and lattice, every step kept as a labelled frame.

    The slow way to a DFT minimum (a 72-atom step is a quarter of an hour on four cores), for
    cells the potential cannot be trusted to relax yet, and a source of frames along the way.
    The named cells are stepped in turn; a run that is stopped goes on from each cell's last
    frame (the optimiser's memory of the curvature starts again, the geometry does not).
    Frames are NAME_dNN.dft.xyz, NN counting DFT evaluations.
    """
    from ase.calculators.calculator import Calculator, all_changes
    from ase.filters import FrechetCellFilter
    from ase.io import read, write
    from ase.optimize import BFGS

    class Dft(Calculator):
        implemented_properties = ["energy", "forces", "stress"]

        def __init__(self, name, count):
            super().__init__()
            self.cell, self.count = name, count  # `name` is the base class's own

        def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
            super().calculate(atoms, properties, system_changes)
            self.count += 1
            frame = f"{self.cell}_d{self.count:02d}"
            work = args.frames / "work" / frame
            shutil.rmtree(work, ignore_errors=True)
            a = self.atoms.copy()
            a.info = dict(start_info[self.cell])
            out, seconds = run_pw(args.pw, args.np, pw_input(a, args.pseudo, work / "scratch", args.conv_thr), work, frame)
            r = read_pw(out)
            if r is None:
                raise SystemExit(f"{frame}: no converged result, see {out}")
            a.info.update({"frame": frame, "family": "dft-relax", "REF_energy": r["energy"], "REF_stress": r["stress"],
                           "dft_energy_ry": r["energy_ry"], "dft_d3_ry": r["d3_ry"], "dft_kpoints": r["kpoints"],
                           "dft_seconds": round(seconds, 1), "dft_conv_thr_ry": args.conv_thr,
                           "note": f"DFT relaxation of {self.cell}, evaluation {self.count}"})
            a.arrays["REF_forces"] = r["forces"]
            write(args.frames / f"{frame}.dft.xyz", a, format="extxyz")
            shutil.copy(out, args.frames / f"{frame}.pw.out")
            shutil.rmtree(work, ignore_errors=True)
            stress = r["stress"]
            self.results = {"energy": r["energy"], "forces": r["forces"],
                            "stress": np.array([stress[0, 0], stress[1, 1], stress[2, 2], stress[1, 2], stress[0, 2], stress[0, 1]])}
            print(f"{frame}: {r['energy'] / a.info['n_monomers']:.5f} eV per monomer, largest force "
                  f"{np.abs(r['forces']).max():.3f} eV/A, stress " + " ".join(f"{v / 0.006241509:+.2f}" for v in np.diag(stress))
                  + f" GPa, {seconds:.0f} s", flush=True)

    began = time.time()
    names = args.names.split(",")
    start_info, optimisers = {}, {}
    for name in names:
        steps = sorted(args.frames.glob(f"{name}_d[0-9][0-9].dft.xyz"))
        last = read(steps[-1] if steps else args.frames / f"{name}.dft.xyz")
        start_info[name] = {k: last.info[k] for k in ("source", "n_monomers", "kmesh") if k in last.info}
        a = last.copy()
        a.calc = Dft(name, len(steps))
        # the last frame is labelled already: its DFT result is handed over, not computed again
        stress = np.asarray(last.info["REF_stress"]).reshape(3, 3)
        a.calc.atoms = a.copy()
        a.calc.results = {"energy": float(last.info["REF_energy"]), "forces": np.array(last.arrays["REF_forces"]),
                          "stress": np.array([stress[0, 0], stress[1, 1], stress[2, 2], stress[1, 2], stress[0, 2], stress[0, 1]])}
        optimisers[name] = (a, BFGS(FrechetCellFilter(a), logfile=None, maxstep=args.maxstep))
    done = {name: False for name in names}
    while not all(done.values()):
        for name in names:
            if done[name]:
                continue
            if args.minutes and time.time() - began > 60 * args.minutes:
                print(f"the {args.minutes:g} minutes are used; not yet relaxed: " + ", ".join(n for n in names if not done[n]), flush=True)
                return
            a, opt = optimisers[name]
            done[name] = bool(opt.run(fmax=args.fmax, steps=1))  # one more step each turn; the curvature is kept
            if done[name]:
                print(f"{name}: relaxed to {args.fmax} eV/A after {a.calc.count} DFT evaluations in all", flush=True)


# An assumed stiffness across the chains (C_11 = C_22, C_12, GPa), only to turn a stress into a
# strain to aim at: PVDF's is 23 to 25 and 6, the candidates' cells are softer.
AIM_C_GPA = (15.0, 5.0)
AIM_LIMIT = 0.12


def aim(args):
    """Move a cell's unlabelled frames to where DFT's stress on its relaxed frame points.

    The sampler's relaxed cells are not DFT's (it put many all-trans cells 20 to 37% above
    PVDF's volume). Once a cell's relaxed frame is labelled, its DFT stress across the chains
    and an assumed stiffness give the strain to DFT's minimum; the cell's other frames, not
    yet labelled, are strained by half of it, all of it and one and a half times it, so that
    the minimum is inside what is labelled whichever way the sampler was wrong.
    """
    from ase import units
    from ase.io import read, write

    c11, c12 = AIM_C_GPA
    moved = 0
    for first in sorted(args.frames.glob("units_*_0.dft.xyz")):
        cell_name = first.name[:8]
        a = read(first)
        s = np.asarray(a.info["REF_stress"]).reshape(3, 3) / units.GPa
        ex = -(c11 * s[0, 0] - c12 * s[1, 1]) / (c11 ** 2 - c12 ** 2)
        ey = -(c11 * s[1, 1] - c12 * s[0, 0]) / (c11 ** 2 - c12 ** 2)
        ex, ey = (float(np.clip(v, -AIM_LIMIT, AIM_LIMIT)) for v in (ex, ey))
        others = [f for f in sorted(args.frames.glob(f"{cell_name}_[1-9].start.xyz"))
                  if not (args.frames / f.name.replace(".start.xyz", ".dft.xyz")).exists()]
        todo = [(f, read(f)) for f in others]
        todo = [(f, b) for f, b in todo if "aimed_strain" not in b.info]
        for (f, b), step in zip(todo, np.linspace(0.5, 1.5, len(todo))):
            info = dict(b.info)
            b = strained(b, np.array([step * ex, step * ey, 0.0, 0.0, 0.0, 0.0]))
            b.info = {**info, "aimed_strain": f"{step * ex:+.4f} {step * ey:+.4f}",
                      "aimed_from_stress_GPa": f"{s[0, 0]:+.3f} {s[1, 1]:+.3f}"}
            write(f, b, format="extxyz")
            moved += 1
        print(f"{cell_name} {str(a.info['source']):16s} DFT stress xx {s[0, 0]:+.2f} yy {s[1, 1]:+.2f} zz {s[2, 2]:+.2f} GPa: "
              f"aim at strain xx {ex:+.3f} yy {ey:+.3f}; {len(todo)} frames moved", flush=True)
    print(f"{moved} frames moved")


def cutoff(args):
    """Labelled frames again at a higher plane-wave cutoff: what the anchor's 90/360 Ry costs.

    The anchor's cutoffs were set for H, C and F. The chlorine file is ultrasoft and the SSSP
    table asks 100/800 Ry of it, so the frames with chlorine (and nitrogen) are checked: the
    forces, the stress, and the energy of one frame against another of the same atoms, which
    is what a margin is. An energy's own shift with the cutoff is the same in both and is not
    an error of a margin.
    """
    from ase import units
    from ase.io import read

    ecut = tuple(int(v) for v in args.ecut.split(","))
    rows = {}
    for name in args.names.split(","):
        a = read(args.frames / f"{name}.dft.xyz")
        work = args.frames / "work" / f"{name}_ecut"
        out, seconds = run_pw(args.pw, args.np, pw_input(a, args.pseudo, work / "scratch", args.conv_thr, ecut), work, name)
        r = read_pw(out)
        shutil.rmtree(work / "scratch", ignore_errors=True)
        if r is None:
            raise SystemExit(f"{name}: no converged result at {ecut}, see {out}")
        n = a.info["n_monomers"]
        df = r["forces"] - a.arrays["REF_forces"]
        ds = (r["stress"] - np.asarray(a.info["REF_stress"]).reshape(3, 3)) / units.GPa
        rows[name] = {
            "source": str(a.info["source"]), "atoms": len(a), "seconds": round(seconds),
            "energy_eV_per_monomer": {"anchor_cutoffs": a.info["REF_energy"] / n, "higher": r["energy"] / n},
            "force_change_rms_meV_per_A": float(1e3 * np.sqrt((df ** 2).mean())),
            "force_change_largest_meV_per_A": float(1e3 * np.abs(df).max()),
            "stress_change_largest_GPa": float(np.abs(ds).max()),
        }
        print(f"{name} {rows[name]['source']}: forces change by {rows[name]['force_change_rms_meV_per_A']:.1f} meV/A rms "
              f"({rows[name]['force_change_largest_meV_per_A']:.1f} largest), stress by "
              f"{rows[name]['stress_change_largest_GPa']:.3f} GPa, {seconds:.0f} s", flush=True)
    names = list(rows)
    pairs = {}
    for a, b in zip(names[::2], names[1::2]):
        gap = {k: 1e3 * (rows[b]["energy_eV_per_monomer"][k] - rows[a]["energy_eV_per_monomer"][k])
               for k in ("anchor_cutoffs", "higher")}
        pairs[f"{b} minus {a}"] = {**{f"{k}_meV_per_monomer": v for k, v in gap.items()},
                                   "change_kcal_per_monomer": 0.023060548 * (gap["higher"] - gap["anchor_cutoffs"])}
        print(f"{rows[b]['source']} minus {rows[a]['source']}: {gap['anchor_cutoffs']:+.3f} meV per monomer at 90/360, "
              f"{gap['higher']:+.3f} at {ecut[0]}/{ecut[1]}: {pairs[f'{b} minus {a}']['change_kcal_per_monomer']:+.4f} kcal/mol",
              flush=True)
    args.out.write_text(json.dumps({"anchor_cutoffs_ry": [90, 360], "higher_cutoffs_ry": list(ecut),
                                    "frames": rows, "pairs": pairs}, indent=1) + "\n")
    print(f"wrote {args.out}")


def collect(args):
    from ase.io import read, write

    frames = [read(f) for f in sorted(args.frames.glob("*.dft.xyz"))]
    for a in frames:
        a.calc = None
    args.out.parent.mkdir(parents=True, exist_ok=True)
    write(args.out, frames, format="extxyz")
    counts = {}
    for a in frames:
        counts[a.info["family"]] = counts.get(a.info["family"], 0) + 1
    print(f"{len(frames)} labelled frames in {args.out}: " + ", ".join(f"{k} {v}" for k, v in counts.items())
          + f"; {sum(a.info['dft_seconds'] for a in frames) / 3600:.1f} hours of DFT")


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = p.add_subparsers(dest="step", required=True)
    for name in ("make", "check", "label", "aim", "cutoff", "relax", "collect"):
        s = sub.add_parser(name)
        s.add_argument("--frames", type=Path, required=True, help="working directory, outside the repository")
        if name == "make":
            s.add_argument("--set", default="pvdf", choices=["pvdf", "near", "units"])
            s.add_argument("--sampler", default="medium-0b3", help="pretrained MACE-MP model for the dynamics")
            s.add_argument("--threads", type=int, default=4)
            s.add_argument("--minutes", type=float, default=0.0,
                           help="units set: start no new cell after this long (the next run goes on)")
        if name in ("check", "label", "cutoff", "relax"):
            s.add_argument("--pw", required=True, help="path to pw.x (Quantum ESPRESSO 7.5)")
            s.add_argument("--pseudo", required=True, help="directory with the SSSP 1.3.0 PBE precision files")
            s.add_argument("--np", type=int, default=4)
            s.add_argument("--pools", type=int, default=0, help="k-point pools (default one a process); 1 for a frame whose arrays "
                                                                 "do not fit in memory once a pool")
        if name == "label":
            s.add_argument("--conv-thr", type=float, default=1e-9)
            s.add_argument("--only", default="", help="families to label, comma-separated")
            s.add_argument("--names", default="", help="frames to label, comma-separated")
            s.add_argument("--minutes", type=float, default=0.0,
                           help="start no new frame after this long, so a run ends inside a time limit")
        if name == "relax":
            s.add_argument("--names", required=True, help="labelled frames to relax, comma-separated; stepped in turn")
            s.add_argument("--fmax", type=float, default=0.05, help="eV/A, on the atoms and on the cell's gradient")
            s.add_argument("--maxstep", type=float, default=0.08)
            s.add_argument("--conv-thr", type=float, default=1e-9)
            s.add_argument("--minutes", type=float, default=0.0, help="start no new DFT evaluation after this long")
        if name == "cutoff":
            s.add_argument("--names", required=True, help="labelled frames, comma-separated; taken in pairs of the same atoms")
            s.add_argument("--ecut", default="100,800", help="wavefunction and density cutoffs, Ry")
            s.add_argument("--conv-thr", type=float, default=1e-9)
        if name in ("collect", "cutoff"):
            s.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    global POOLS
    POOLS = int(getattr(args, "pools", 0))
    {"make": make, "check": check, "label": label, "aim": aim, "cutoff": cutoff, "relax": dft_relax,
     "collect": collect}[args.step](args)


if __name__ == "__main__":
    main()
