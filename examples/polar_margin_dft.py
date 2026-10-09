"""Polar-versus-antipolar lattice energy margin of all-trans PVDF from periodic DFT (PBE-D3(BJ), Quantum ESPRESSO).

    python examples/polar_margin_dft.py starts [--alpha]          # polyfind's cells -> deliverables/polar_margin_dft/starts
    python examples/polar_margin_dft.py inputs --cell C --stage S # one pw.x input (run root, pseudo and scratch dirs are args)
    python examples/polar_margin_dft.py check --cell C --stage S  # exit 0 if the stage's last SCF meets the criteria
    python examples/polar_margin_dft.py collect                   # every output -> raw_stages.json, inputs/, relaxed/
    python examples/polar_margin_dft.py report                    # raw_stages.json -> polar_margin_dft.json

The margin the screen cannot resolve (``docs/SCREEN.md``, "Polar or antipolar"; ``docs/SCREEN_PHONONS.md``,
"What the gamma-free cells say about the screen's column") is computed here on the same Hamiltonian as the
provider's Quantum ESPRESSO references in ``sarco/materials/gpu_bundle/periodic_reference/beta_pvdf``
(``born/scf.in``): PBE, SSSP 1.3.0 PBE precision pseudopotentials (H_ONCV_PBE-1.0.oncvpsp.upf,
C.pbe-n-kjpaw_psl.1.0.0.UPF, F.oncvpsp.upf), 90/360 Ry, ``vdw_corr='grimme-d3'``, ``dftd3_version=4``
(zero-damping is 3, Becke-Johnson is 4), ``dftd3_threebody=.false.``, fixed occupations.

``starts`` regenerates the four cells of ``examples/screen_phonons.py --only pvdf`` exactly as that script
does (it calls its :func:`one_chemistry`) and records, for every cell, both the rigid-chain placement and
the all-atom-relaxed placement at fixed cell on the phonon packer (``phonon.relax_all_atom``, accepted only
when its own criterion is met).  The all-atom-relaxed placements are the DFT starts.  The cells:

* ``polar``      beta-PVDF, the polar all-trans cell (``polar90``; with gamma free it is the same cell);
* ``anti_free``  the gamma-free exactly-antipolar competitor (``anti_free``), degenerate with beta on the
                 classical potential (gap -0.006 kcal/mol per monomer after all-atom relaxation);
* ``anti90``     the gamma = 90 exactly-antipolar cell the screen compares (``anti90``), a saddle on the
                 classical potential (its chains slide 1.05 A when let go).

plus ``polar_dipole``: beta-PVDF from the provider's own geometry (``born/scf.in``, their CP2K PBE-D3(BJ)
zero-pressure cell), a cross-check that both starts reach the same DFT minimum; and ``alpha`` (TGTG',
24 atoms) for the alpha/beta polymorph-gap control.

``inputs`` writes the pw.x input for one stage of the protocol (``STAGES``), taking its geometry from the
start or from the stage it follows; ``deliverables/polar_margin_dft/runner/drive_cell.sh`` chains the stages
and runs pw.x.  ``collect`` reads every output under the run root; ``report`` derives the margins, their
convergence and D3 share, the lattices against experiment, and puts polyfind's margins beside them.  Paths
default to this machine's layout (WSL home for runs) and are arguments.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
OUT = os.path.join(REPO, "deliverables", "polar_margin_dft")
STARTS = os.path.join(OUT, "starts")
SARCO_BETA = r"E:\source\gh\sarco\materials\gpu_bundle\periodic_reference\beta_pvdf"

KCAL_PER_EV = 23.060547830619
RY_EV = 13.605693122994
RY_KCAL = RY_EV * KCAL_PER_EV
BOHR_A = 0.529177210903

# The provider's species card, byte for byte (born/scf.in).
SPECIES = [("H", 1.00794, "H_ONCV_PBE-1.0.oncvpsp.upf"),
           ("C", 12.0107, "C.pbe-n-kjpaw_psl.1.0.0.UPF"),
           ("F", 18.998403163, "F.oncvpsp.upf")]
MASS = {s: m for s, m, _ in SPECIES}

# polyfind's four cells -> the DFT cell names
CELLS = {"polar": "polar90", "anti_free": "anti_free", "anti90": "anti90"}


# ------------------------------------------------------------------------------------------- geometry I/O
def write_xyz(path: str, symbols, positions, lattice, comment: dict) -> None:
    lat = np.asarray(lattice, dtype=float)
    pos = np.asarray(positions, dtype=float)
    head = 'Lattice="' + " ".join(f"{v:.10f}" for v in lat.ravel()) + '" Properties=species:S:1:pos:R:3 pbc="T T T"'
    for k, v in comment.items():
        head += f" {k}={json.dumps(v) if isinstance(v, str) else v}"
    with open(path, "w", newline="\n") as f:
        f.write(f"{len(symbols)}\n{head}\n")
        for s, p in zip(symbols, pos):
            f.write(f"{s:<2s} {p[0]:16.10f} {p[1]:16.10f} {p[2]:16.10f}\n")


def read_xyz(path: str):
    with open(path) as f:
        lines = f.read().splitlines()
    n = int(lines[0])
    m = re.search(r'Lattice="([^"]+)"', lines[1])
    lat = np.array([float(v) for v in m.group(1).split()]).reshape(3, 3)
    sym, pos = [], []
    for ln in lines[2:2 + n]:
        t = ln.split()
        sym.append(t[0])
        pos.append([float(v) for v in t[1:4]])
    return sym, np.array(pos), lat


def sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        h.update(f.read())
    return h.hexdigest()


def reduce_inplane(lat: np.ndarray) -> tuple:
    """Lagrange-Gauss reduction of the two in-plane vectors (rows 0, 1); the chain axis (row 2) is kept.

    Cartesian positions are unchanged by a change of basis, so the reduced cell holds the same crystal;
    the reduced (shortest, most nearly orthogonal) basis gives the cleanest Monkhorst-Pack meshes.
    """
    a, b, c = (np.array(v, dtype=float) for v in lat)
    for _ in range(100):
        if np.dot(a, a) > np.dot(b, b):
            a, b = b, a
        mu = round(float(np.dot(a, b) / np.dot(a, a)))
        if mu == 0:
            break
        b = b - mu * a
    if np.cross(a, b)[2] < 0:
        b = -b
    # put the first vector along x, as polyfind does (a rotation about z, applied to the whole cell; the
    # caller applies the same R to the positions).  Returns ``(lattice, R)``.
    ang = math.atan2(a[1], a[0])
    R = np.array([[math.cos(-ang), -math.sin(-ang), 0.0], [math.sin(-ang), math.cos(-ang), 0.0], [0.0, 0.0, 1.0]])
    return np.array([R @ a, R @ b, R @ c]), R


def cell_params(lat: np.ndarray) -> dict:
    L = np.linalg.norm(lat, axis=1)
    def ang(u, v):
        return math.degrees(math.acos(float(np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v)))))
    return {"a": float(L[0]), "b": float(L[1]), "c": float(L[2]), "alpha": ang(lat[1], lat[2]),
            "beta": ang(lat[0], lat[2]), "gamma": ang(lat[0], lat[1]), "volume": float(abs(np.linalg.det(lat)))}


# ------------------------------------------------------------------------------------------- starts
def build_starts(out: str = STARTS) -> dict:
    """Run ``screen_phonons.one_chemistry`` for PVDF and capture every cell's rigid and relaxed placement."""
    sys.path.insert(0, HERE)
    import screen_phonons as SP  # noqa: E402
    from polyfind.polymers import get_polymer  # noqa: E402

    captured: dict = {}
    current = {"label": None}
    orig_stage, orig_relax = SP.phonon_stage, SP.relax_all_atom

    def stage(full, plain, params, label, h, n_report=12):
        current["label"] = label
        try:
            return orig_stage(full, plain, params, label, h, n_report)
        finally:
            current["label"] = None

    def relax(pk, params, Pn=None, latn=None, **kw):
        rel = orig_relax(pk, params, Pn=Pn, latn=latn, **kw)
        captured[current["label"]] = {"params": np.asarray(params, dtype=float).copy(), "P0": np.asarray(Pn).copy(),
                                      "latn": np.asarray(latn).copy(), "elements": list(pk.elements) * pk.n_chains,
                                      "relaxation": rel, "n_chains": pk.n_chains, "n_monomers": pk.chain.n_monomers}
        return rel

    SP.phonon_stage, SP.relax_all_atom = stage, relax
    t0 = time.time()
    try:
        rec = SP.one_chemistry(get_polymer("pvdf"), h=1e-3, anti_target=6000, anti_polish=8, anti_maxfev=900)
    finally:
        SP.phonon_stage, SP.relax_all_atom = orig_stage, orig_relax
    os.makedirs(out, exist_ok=True)
    n_mon = rec["n_monomers_per_cell"]
    summary = {"source": "examples/screen_phonons.py one_chemistry(pvdf), defaults of its main()",
               "preset": rec["preset"], "seconds": round(time.time() - t0, 1), "n_monomers_per_cell": n_mon,
               "screen": rec["screen"], "phonon_packer": rec["phonon_packer"], "cells": {}}
    for label, cap in captured.items():
        rel = cap["relaxation"]
        if not rel.converged:
            raise RuntimeError(f"{label}: all-atom relaxation not accepted ({rel.summary()}); no start written")
        for kind, P in (("rigid", cap["P0"]), ("relaxed", rel.Pn)):
            path = os.path.join(out, f"pvdf_{label}_{kind}.xyz")
            write_xyz(path, cap["elements"], P, cap["latn"],
                      {"polyfind_cell": label, "placement": kind, "params": json.dumps([float(v) for v in cap["params"]]),
                       "energy_kcal_per_cell": float(rel.energy) if kind == "relaxed" else float(rec["cells"][label]["rigid"]["energy"])})
        summary["cells"][label] = {
            "params": [float(v) for v in cap["params"]],
            "rigid_energy_kcal_per_cell": rec["cells"][label]["rigid"]["energy"],
            "relaxed_energy_kcal_per_cell": float(rel.energy),
            "relaxation": rel.summary(),
            "max_displacement": rec["cells"][label]["relaxed"]["max_displacement"],
            "chain_slides": rec["cells"][label]["relaxed"]["chain_slides"],
            "dipole_per_monomer_after": rec["cells"][label]["relaxed"]["dipole_per_monomer"],
            "lowest_optical_relaxed": rec["cells"][label]["relaxed"]["lowest_optical"],
            "files": {k: f"pvdf_{label}_{k}.xyz" for k in ("rigid", "relaxed")},
        }
    c = summary["cells"]

    def gap(a, p, kind):
        return (c[a][f"{kind}_energy_kcal_per_cell"] - c[p][f"{kind}_energy_kcal_per_cell"]) / n_mon

    summary["polyfind_margins_kcal_per_monomer"] = {
        "screen_gamma90_rigid (anti90 - best, screen packer)": rec["screen"]["gamma90"]["antipolar_gap"],
        "screen_gamma_free_rigid (anti_free - best, screen packer)": rec["screen"]["gamma_free"]["antipolar_gap"],
        "phonon_packer_rigid anti_free - polar": gap("anti_free", "polar90", "rigid"),
        "phonon_packer_rigid anti90 - polar": gap("anti90", "polar90", "rigid"),
        "phonon_packer_all_atom_relaxed anti_free - polar": gap("anti_free", "polar90", "relaxed"),
        "phonon_packer_all_atom_relaxed anti90 - polar": gap("anti90", "polar90", "relaxed"),
        "error_bar": rec["screen"]["error_bar"],
    }
    with open(os.path.join(out, "polyfind_starts.json"), "w", newline="\n") as f:
        json.dump(summary, f, indent=1)
    print(json.dumps(summary["polyfind_margins_kcal_per_monomer"], indent=1))
    return summary


def build_alpha(out: str = STARTS) -> dict:
    """alpha-PVDF (TGTG', 24 atoms) for the polymorph-gap control: the lowest antipolar two-chain cell.

    The chain is ``pvdf_polymorphs.py``'s (``periodic_chain(PVDF, [T, G+, T, G-], THREE_STATE)``, ideal
    angles) packed by ``pack()`` under the screen's fitted potential and Ewald sum, gamma free; among the
    returned cells the lowest one with ``|P| < 0.01 C/m^2`` is taken (alpha is the antipolar packing of
    this chain, delta its polar partner).  DFT relaxes everything, so this only has to be a sensible start.
    """
    sys.path.insert(0, HERE)
    import screen_electroactive as SE  # noqa: E402
    import polyfind.pack as pack_mod  # noqa: E402
    from polyfind.ewald import EwaldSpec  # noqa: E402
    from polyfind.polymers import PVDF, THREE_STATE  # noqa: E402

    T, GP, GM = 0, 1, 2
    chain = pack_mod.periodic_chain(PVDF, [T, GP, T, GM], THREE_STATE)
    FV, _, _ = SE._presets()
    with FV.applied():
        res = pack_mod.pack(chain, n_chains=2, table_cache_dir=None, coulomb="ewald", ewald=EwaldSpec(),
                            gamma_free=True, rng=np.random.default_rng(0))
    anti = [r for r in res if r.polarization_magnitude < 1e-2]
    if not anti:
        raise RuntimeError("pack() returned no antipolar alpha cell")
    r = anti[0]
    pk = pack_mod.CrystalPacker(chain, n_chains=2)
    P, lat = pk._place(r.params[None])
    from polyfind import backend as bk  # noqa: E402
    P = np.asarray(bk.to_numpy(P), dtype=float)[0]
    lat = np.asarray(bk.to_numpy(lat), dtype=float)[0]
    path = os.path.join(out, "pvdf_alpha_start.xyz")
    write_xyz(path, list(chain.elements) * 2, P, lat,
              {"polyfind_cell": "alpha TGTG' lowest antipolar, pack(gamma_free) under pvdf-dft-valence, Ewald",
               "row": r.row(), "params": json.dumps([float(v) for v in r.params])})
    best = res[0]
    info = {"alpha_row": r.row(), "alpha_e_per_monomer": float(r.energy_per_monomer),
            "best_row": best.row(), "best_e_per_monomer": float(best.energy_per_monomer), "file": os.path.basename(path)}
    print(json.dumps(info, indent=1))
    return info


def dipole_start(out: str = STARTS, sarco: str = SARCO_BETA) -> str:
    """The provider's beta cell (``born/scf.in``) as an extended xyz start; returns its path."""
    src = os.path.join(sarco, "born", "scf.in")
    txt = open(src).read()
    cell = re.search(r"CELL_PARAMETERS angstrom\n(.*?)\n\n", txt, re.S).group(1).split("\n")
    lat = np.array([[float(v) for v in ln.split()] for ln in cell])
    atoms = re.search(r"ATOMIC_POSITIONS angstrom\n(.*?)\n\n", txt, re.S).group(1).split("\n")
    sym = [ln.split()[0] for ln in atoms]
    pos = np.array([[float(v) for v in ln.split()[1:4]] for ln in atoms])
    path = os.path.join(out, "pvdf_polar_dipole_start.xyz")
    write_xyz(path, sym, pos, lat, {"source": "sarco materials/gpu_bundle/periodic_reference/beta_pvdf/born/scf.in",
                                    "source_sha256": sha256(src)})
    return path


# ------------------------------------------------------------------------------------------- QE inputs
def start_geometry(cell: str, starts: str = STARTS):
    """``(symbols, positions, lattice)`` of a DFT cell's start, in-plane basis reduced."""
    if cell == "polar_dipole":
        sym, pos, lat = read_xyz(os.path.join(starts, "pvdf_polar_dipole_start.xyz"))
        return sym, pos, lat
    name = "pvdf_alpha_start.xyz" if cell == "alpha" else f"pvdf_{CELLS[cell]}_relaxed.xyz"
    sym, pos, lat = read_xyz(os.path.join(starts, name))
    lat2, R = reduce_inplane(lat)
    pos2 = pos @ R.T
    return sym, pos2, lat2


def pw_input(calc: str, prefix: str, sym, pos, lat, kmesh, pseudo_dir: str, outdir: str,
             ecut=(90, 360), threebody: bool = False, conv_thr: float = 1e-10, extra_control: str = "",
             nosym: bool = False) -> str:
    species = [s for s in SPECIES if s[0] in set(sym)]
    idx = range(len(sym))
    ctrl = [f"  calculation = '{calc}'", f"  prefix = '{prefix}'", f"  pseudo_dir = '{pseudo_dir}'",
            f"  outdir = '{outdir}'", "  tstress = .true.", "  tprnfor = .true.", "  disk_io = 'low'"]
    if calc in ("relax", "vc-relax"):
        ctrl += ["  etot_conv_thr = 1.0d-6", "  forc_conv_thr = 1.0d-4", "  nstep = 300"]
    ctrl += [ln for ln in extra_control.splitlines() if ln.strip()]
    system = ["  ibrav = 0", f"  nat = {len(sym)}", f"  ntyp = {len(species)}", "  input_dft = 'PBE'",
              "  occupations = 'fixed'", f"  ecutwfc = {ecut[0]}", f"  ecutrho = {ecut[1]}",
              "  vdw_corr = 'grimme-d3'", "  dftd3_version = 4",
              f"  dftd3_threebody = {'.true.' if threebody else '.false.'}"]
    if nosym:
        system += ["  nosym = .true."]
    txt = "&CONTROL\n" + "\n".join(ctrl) + "\n/\n\n&SYSTEM\n" + "\n".join(system) + "\n/\n\n"
    thr = f"{conv_thr:.1e}".replace("e", "d")
    txt += (f"&ELECTRONS\n  conv_thr = {thr}\n  electron_maxstep = 200\n"
            "  mixing_mode = 'plain'\n  mixing_beta = 0.3\n  startingwfc = 'atomic+random'\n/\n\n")
    if calc in ("relax", "vc-relax"):
        txt += "&IONS\n  ion_dynamics = 'bfgs'\n/\n\n"
    if calc == "vc-relax":
        txt += ("&CELL\n  cell_dynamics = 'bfgs'\n  press = 0.0\n  press_conv_thr = 0.5\n"
                "  cell_dofree = 'all'\n/\n\n")
    txt += "ATOMIC_SPECIES\n" + "".join(f"{s} {m} {p}\n" for s, m, p in species) + "\n"
    txt += "CELL_PARAMETERS angstrom\n" + "".join(f"{v[0]:.10f} {v[1]:.10f} {v[2]:.10f}\n" for v in lat) + "\n"
    txt += "ATOMIC_POSITIONS angstrom\n" + "".join(
        f"{sym[i]:<2s} {pos[i][0]: .10f} {pos[i][1]: .10f} {pos[i][2]: .10f}\n" for i in idx) + "\n"
    txt += f"K_POINTS automatic\n{kmesh[0]} {kmesh[1]} {kmesh[2]} 0 0 0\n"
    return txt


# ------------------------------------------------------------------------------------------- QE output parsing
def parse_pw(path: str) -> dict:
    """Final-state observables of one pw.x output (the last SCF in it)."""
    txt = open(path, errors="replace").read()
    r: dict = {"file": os.path.basename(path), "job_done": "JOB DONE" in txt}
    e = re.findall(r"^!\s+total energy\s+=\s+(-?\d+\.\d+) Ry", txt, re.M)
    d3 = re.findall(r"DFT-D3 Dispersion\s+=\s+(-?\d+\.\d+) Ry", txt)
    r["energy_ry"] = float(e[-1]) if e else None
    r["d3_ry"] = float(d3[-1]) if d3 else None
    r["n_scf"] = len(e)
    k = re.findall(r"number of k points=\s+(\d+)", txt)
    r["n_kpoints_irreducible"] = int(k[-1]) if k else None
    s = re.findall(r"(\d+) Sym\. Ops\.", txt)
    r["sym_ops"] = int(s[-1]) if s else None
    p = re.findall(r"P=\s*(-?\d+\.\d+)", txt)
    r["pressure_kbar"] = float(p[-1]) if p else None
    st = re.findall(r"total   stress .*\n(.*)\n(.*)\n(.*)\n", txt)
    if st:
        r["stress_kbar"] = [[float(v) for v in ln.split()[3:6]] for ln in st[-1]]
    f = re.findall(r"Total force =\s+(\d+\.\d+)", txt)
    r["total_force_ry_bohr"] = float(f[-1]) if f else None
    blocks = re.findall(r"Forces acting on atoms.*?\n\n(.*?)\n\n", txt, re.S)
    if blocks:
        fr = [[float(v) for v in ln.split()[-3:]] for ln in blocks[-1].splitlines() if "force =" in ln]
        r["max_force_component_ry_bohr"] = float(np.abs(np.array(fr)).max()) if fr else None
    g = re.findall(r"highest occupied, lowest unoccupied level \(ev\):\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)", txt)
    if g:
        r["gap_ev"] = float(g[-1][1]) - float(g[-1][0])
    w = re.findall(r"PWSCF\s+:.*CPU\s+(.*) WALL", txt)
    r["wall"] = w[-1].strip() if w else None
    r["bfgs_converged"] = "bfgs converged" in txt
    r["bfgs_steps"] = len(re.findall(r"number of bfgs steps\s+=", txt)) or None
    r["has_final_coordinates"] = ("End final coordinates" in txt)
    return r


def final_geometry(path: str):
    """``(symbols, positions A, lattice A)`` of the last CELL_PARAMETERS / ATOMIC_POSITIONS block of a vc-relax."""
    txt = open(path, errors="replace").read()
    cells = re.findall(r"CELL_PARAMETERS \((alat=\s*([\d.]+)|angstrom)\)\n(.*?)\n(.*?)\n(.*?)\n", txt)
    if not cells:
        return None
    hdr, alat, *rows = cells[-1]
    lat = np.array([[float(v) for v in r.split()] for r in rows])
    if alat:
        lat = lat * float(alat) * BOHR_A
    blocks = re.findall(r"ATOMIC_POSITIONS \((angstrom|crystal|bohr|alat)\)\n(.*?)\n(?:End|\n)", txt, re.S)
    unit, body = blocks[-1]
    sym, pos = [], []
    for ln in body.strip().splitlines():
        t = ln.split()
        sym.append(t[0])
        pos.append([float(v) for v in t[1:4]])
    pos = np.array(pos)
    if unit == "crystal":
        pos = pos @ lat
    elif unit == "bohr":
        pos = pos * BOHR_A
    elif unit == "alat":
        pos = pos * float(alat) * BOHR_A
    return sym, pos, lat


# ------------------------------------------------------------------------------------------- the protocol
# Monkhorst-Pack meshes are given as (short in-plane axis, long in-plane axis, chain axis) and mapped onto
# each cell's own axes, so every cell gets the same sampling density (all four cells hold two chains in
# nearly the same area).  "k2" is the provider's 4x8x16 (their cell: 8.36 long x 4.73 short x 2.58 chain).
MESHES = {"k0": (2, 1, 4), "k1": (4, 2, 8), "k15": (6, 3, 12), "k2": (8, 4, 16), "k3": (12, 6, 24)}

# stage -> (calculation, mesh, (ecutwfc, ecutrho), three-body, geometry source)
STAGES = {
    "relax0": ("vc-relax", "k0", (90, 360), False, "start"),
    "relax1": ("vc-relax", "k1", (90, 360), False, "relax0"),
    "relax2": ("vc-relax", "k2", (90, 360), False, "final"),   # only if the k2 SCF fails the criteria
    "scf_k0": ("scf", "k0", (90, 360), False, "final"),
    "scf_k1": ("scf", "k1", (90, 360), False, "final"),
    "scf_k15": ("scf", "k15", (90, 360), False, "final"),
    "scf_k2": ("scf", "k2", (90, 360), False, "final"),
    "scf_k3": ("scf", "k3", (90, 360), False, "final"),
    "scf_k1_e110": ("scf", "k1", (110, 440), False, "final"),
    "scf_k1_e130": ("scf", "k1", (130, 520), False, "final"),
    "scf_k1_3body": ("scf", "k1", (90, 360), True, "final"),
    # symmetry-breaking check: the final geometry, every atom displaced (0.03 A rms, seeded) and the cell
    # strained (0.5 %), relaxed again with nosym; a symmetric stationary point that is not a minimum
    # would fall to a lower energy here.  Not a geometry source for any other stage.
    "relax_nosym": ("vc-relax", "k1", (90, 360), False, "final"),
}
PERTURB_SEED, PERTURB_RMS_A, PERTURB_STRAIN = 7, 0.03, 0.005
RELAX_ORDER = ["relax0", "relax1", "relax1_r1", "relax1_r2", "relax2", "relax2_r1"]
FORCE_THR = 1.0e-4   # Ry/bohr, largest Cartesian component
STRESS_THR = 0.5     # kbar, largest stress-tensor component
WSL_ROOT = "/home/niall_sweeny/polar_margin"


def mesh_for(lat: np.ndarray, name: str) -> tuple:
    short, long_, chain = MESHES[name]
    rec = 2 * np.pi * np.linalg.inv(lat).T
    d = 2 * np.pi / np.linalg.norm(rec, axis=1)  # interplanar spacings
    if d[2] > 3.5:  # a two-monomer repeat (alpha, c ~ 4.6 A): half the chain-axis divisions, same density
        chain = max(1, chain // 2)
    if d[0] <= d[1]:
        return (short, long_, chain)
    return (long_, short, chain)


def stage_spec(stage: str) -> tuple:
    base = re.sub(r"_r\d+$", "", stage)
    calc, mesh, ecut, three, src = STAGES[base]
    if stage != base:  # a repeat of a relaxation from its own predecessor's final geometry
        n = int(stage.rsplit("_r", 1)[1])
        src = base if n == 1 else f"{base}_r{n - 1}"
    return calc, mesh, ecut, three, src


def stage_out(root: str, cell: str, stage: str) -> str:
    return os.path.join(root, cell, stage, f"{stage}.out")


def last_relax(root: str, cell: str):
    """The latest relaxation stage of ``cell`` whose output finished, or ``None``."""
    done = [s for s in RELAX_ORDER if os.path.exists(stage_out(root, cell, s))
            and "JOB DONE" in open(stage_out(root, cell, s), errors="replace").read()]
    return done[-1] if done else None


def geometry_for(root: str, cell: str, src: str, starts: str):
    if src == "start":
        return start_geometry(cell, starts)
    if src == "final":
        src = last_relax(root, cell)
        if src is None:
            raise RuntimeError(f"{cell}: no finished relaxation to take a geometry from")
    g = final_geometry(stage_out(root, cell, src))
    if g is None:
        raise RuntimeError(f"{cell}: no geometry in {stage_out(root, cell, src)}")
    sym, pos, lat = g
    lat2, R = reduce_inplane(lat)
    return sym, pos @ R.T, lat2


def write_stage(root: str, cell: str, stage: str, starts: str, pseudo_dir: str, scratch: str) -> str:
    calc, mesh, ecut, three, src = stage_spec(stage)
    sym, pos, lat = geometry_for(root, cell, src, starts)
    k = mesh_for(lat, mesh)
    nosym = stage.startswith("relax_nosym")
    if nosym:
        rng = np.random.default_rng(PERTURB_SEED)
        eps = rng.normal(scale=PERTURB_STRAIN, size=(3, 3))
        eps = 0.5 * (eps + eps.T)
        pos = pos @ (np.eye(3) + eps).T + rng.normal(scale=PERTURB_RMS_A / math.sqrt(3), size=pos.shape)
        lat = lat @ (np.eye(3) + eps).T
    d = os.path.join(root, cell, stage)
    os.makedirs(d, exist_ok=True)
    txt = pw_input(calc, f"pvdf_{cell}", sym, pos, lat, k, pseudo_dir, f"{scratch}/{cell}_{stage}", ecut=ecut,
                   threebody=three, conv_thr=1e-10 if calc != "scf" else 1e-11, nosym=nosym)
    path = os.path.join(d, f"{stage}.in")
    with open(path, "w", newline="\n") as f:
        f.write(txt)
    print(f"{path}: {calc} mesh {k} ecut {ecut} three-body {three} geometry from {src}")
    return path


def check_stage(root: str, cell: str, stage: str) -> tuple:
    """``(ok, report)``: the last SCF of the stage meets the force and stress criteria (and BFGS converged)."""
    r = parse_pw(stage_out(root, cell, stage))
    smax = max(abs(v) for row in r.get("stress_kbar", [[float("inf")]]) for v in row)
    fmax = r.get("max_force_component_ry_bohr")
    fmax = float("inf") if fmax is None else fmax
    calc = stage_spec(stage)[0]
    ok = r["job_done"] and fmax < FORCE_THR and smax < STRESS_THR and (calc == "scf" or r["bfgs_converged"])
    return ok, f"{cell}/{stage}: done {r['job_done']} bfgs {r['bfgs_converged']} max force {fmax:.2e} Ry/bohr, max |stress| {smax:.3f} kbar"


# ------------------------------------------------------------------------------------------- analysis
def chain_orientation(sym, pos, lat) -> dict:
    """Per chain, the in-plane unit vector from each CF2 carbon toward its two fluorines' midpoint.

    The chain's dipole points opposite to it (F is the negative end).  Chains are told apart by
    connectivity (C-C bonds < 1.8 A across periodic images); a polar cell has every chain's vector
    parallel, an antipolar one has them cancelling.
    """
    sym = list(sym)
    pos = np.asarray(pos, dtype=float)
    lat = np.asarray(lat, dtype=float)
    inv = np.linalg.inv(lat)

    def mic(v):
        f = v @ inv
        f -= np.round(f)
        return f @ lat

    n = len(sym)
    parent = list(range(n))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i in range(n):
        for j in range(i + 1, n):
            lim = 1.8 if (sym[i], sym[j]) == ("C", "C") else (1.25 if "H" in (sym[i], sym[j]) else 1.6)
            if "C" in (sym[i], sym[j]) and np.linalg.norm(mic(pos[j] - pos[i])) < lim:
                parent[find(i)] = find(j)
    groups: dict = {}
    for i in range(n):
        groups.setdefault(find(i), []).append(i)
    chains = []
    centroids = []
    for g in groups.values():
        cs = [i for i in g if sym[i] == "C"]
        if cs:
            # in-plane: the carbon centroid; axial: the height of the first CH2 carbon (the centroid's z is
            # ambiguous by half a repeat for a one-monomer chain, depending on which image is nearest)
            ref = pos[cs[0]]
            cen = ref + np.mean([mic(pos[i] - ref) for i in cs], axis=0)
            ch2 = [i for i in cs if any(sym[j] == "H" and np.linalg.norm(mic(pos[j] - pos[i])) < 1.25 for j in g)]
            if ch2:
                cen[2] = pos[min(ch2, key=lambda i: (pos[i] @ inv)[2] % 1.0)][2]
            centroids.append(cen)
        vecs = []
        for i in g:
            if sym[i] != "C":
                continue
            fs = [j for j in g if sym[j] == "F" and np.linalg.norm(mic(pos[j] - pos[i])) < 1.6]
            if len(fs) == 2:
                v = (mic(pos[fs[0]] - pos[i]) + mic(pos[fs[1]] - pos[i])) / 2.0
                vecs.append(v)
        if vecs:
            v = np.mean(vecs, axis=0)
            vin = v.copy()
            vin[2] = 0.0
            chains.append({"atoms": len(g), "cf2_groups": len(vecs),
                           "cf2_to_f_unit_inplane": [round(float(x), 4) for x in vin / np.linalg.norm(vin)],
                           "azimuth_deg": round(math.degrees(math.atan2(vin[1], vin[0])), 2)})
    s = np.sum([c["cf2_to_f_unit_inplane"] for c in chains], axis=0)
    out = {"chains": chains, "sum_of_unit_vectors": [round(float(x), 4) for x in s],
           "arrangement": "polar" if np.linalg.norm(s) > 1.0 else "antipolar"}
    if len(centroids) == 2:
        # chain registry: the carbon centroid of chain 2 relative to chain 1, fractional, wrapped to [-0.5, 0.5)
        f = (centroids[1] - centroids[0]) @ inv
        f -= np.floor(f + 0.5)
        out["chain2_minus_chain1_fractional"] = [round(float(x), 4) for x in f]
        out["chain_separation_inplane_A"] = round(float(np.linalg.norm(((f @ lat) * np.array([1, 1, 0])))), 4)
    return out


def file_hashes(path: str) -> dict:
    data = open(path, "rb").read()
    return {"md5": hashlib.md5(data).hexdigest(), "sha256": hashlib.sha256(data).hexdigest(), "bytes": len(data)}


def collect(root: str, out_json: str, starts: str, pseudo_dir: str | None = None, out: str = OUT) -> dict:
    """Parse every stage output; copy every input into ``out/inputs``; write each final geometry to ``out/relaxed``."""
    res: dict = {"run_root": root, "cells": {}}
    if pseudo_dir and os.path.isdir(pseudo_dir):
        res["pseudopotentials"] = {p: file_hashes(os.path.join(pseudo_dir, p)) for _, _, p in SPECIES}
    cells = sorted(d for d in os.listdir(root) if os.path.isdir(os.path.join(root, d)) and d in (*CELLS, "polar_dipole", "alpha"))
    for cell in cells:
        c: dict = {"stages": {}}
        for stage in sorted(os.listdir(os.path.join(root, cell))):
            p = stage_out(root, cell, stage)
            if not os.path.exists(p):
                continue
            r = parse_pw(p)
            v = re.search(r"Program PWSCF (v\.\S+) starts", open(p, errors="replace").read())
            if v:
                res.setdefault("qe_versions", sorted(set(res.get("qe_versions", []) + [v.group(1)])))
            r["mesh"] = None
            inp = os.path.join(root, cell, stage, f"{stage}.in")
            if os.path.exists(inp):
                txt = open(inp).read()
                m = re.search(r"K_POINTS automatic\n(\d+) (\d+) (\d+)", txt)
                r["mesh"] = [int(x) for x in m.groups()] if m else None
                r["ecut"] = [int(x) for x in re.findall(r"ecut(?:wfc|rho) = (\d+)", txt)]
                r["threebody"] = ".true." in re.search(r"dftd3_threebody = (\S+)", txt).group(1)
                r["nosym"] = "nosym = .true." in txt
                os.makedirs(os.path.join(out, "inputs", cell), exist_ok=True)
                with open(os.path.join(out, "inputs", cell, f"{stage}.in"), "w", newline="\n") as f:
                    f.write(txt)
            r["meets_criteria"] = check_stage(root, cell, stage)[0] if r["job_done"] else False
            if stage.startswith("relax") and r["job_done"]:
                g = final_geometry(p)
                if g is not None:
                    sym, pos, lat = g
                    lat2, R = reduce_inplane(lat)
                    r["final_lattice"] = cell_params(lat2)
                    r["final_orientation"] = chain_orientation(sym, pos @ R.T, lat2)
                    if stage == "relax_nosym":
                        os.makedirs(os.path.join(out, "relaxed"), exist_ok=True)
                        write_xyz(os.path.join(out, "relaxed", f"pvdf_{cell}_nosym_pbe_d3bj.xyz"), sym, pos @ R.T, lat2,
                                  {"cell": cell, "relaxation": "relax_nosym (perturbed start, no symmetry, 4x2x8 mesh)",
                                   "hamiltonian": "PBE-D3(BJ) two-body, SSSP 1.3.0 PBE precision, 90/360 Ry"})
            c["stages"][stage] = r
        fr = last_relax(root, cell)
        c["final_relaxation"] = fr
        if fr:
            sym, pos, lat = geometry_for(root, cell, "final", starts)
            c["lattice"] = cell_params(lat)
            c["lattice_vectors_A"] = [[round(float(x), 8) for x in row] for row in lat]
            c["orientation"] = chain_orientation(sym, pos, lat)
            c["n_atoms"] = len(sym)
            c["n_monomers"] = sum(1 for s in sym if s == "C") // 2
            os.makedirs(os.path.join(out, "relaxed"), exist_ok=True)
            write_xyz(os.path.join(out, "relaxed", f"pvdf_{cell}_pbe_d3bj.xyz"), sym, pos, lat,
                      {"cell": cell, "relaxation": fr, "hamiltonian": "PBE-D3(BJ) two-body, SSSP 1.3.0 PBE precision, 90/360 Ry"})
        res["cells"][cell] = c
    with open(out_json, "w", newline="\n") as f:
        json.dump(res, f, indent=1)
    return res


# experiment: beta, Hasegawa et al. 1972 (sarco hasegawa_planar.cif); alpha, docs/REFERENCES.md 1.2
EXPERIMENT = {"beta": {"a_long": 8.58, "b_short_polar": 4.91, "c_chain": 2.56, "source": "Hasegawa et al. 1972"},
              "alpha": {"a_short": 4.96, "b_long": 9.64, "c_chain": 4.62, "source": "docs/REFERENCES.md 1.2 (Hasegawa et al. 1972)"}}
SCF_STAGES = ["scf_k0", "scf_k1", "scf_k15", "scf_k2", "scf_k3", "scf_k1_e110", "scf_k1_e130", "scf_k1_3body"]


def _per_mon(cell: dict, stage: str, key: str = "energy_ry"):
    st = cell["stages"].get(stage)
    if not st or not st.get("job_done") or st.get(key) is None:
        return None
    return st[key] / cell["n_monomers"]


def report(raw_json: str, starts: str, out_json: str) -> dict:
    raw = json.load(open(raw_json))
    pf = json.load(open(os.path.join(starts, "polyfind_starts.json")))
    cells = raw["cells"]
    meV, kcal = RY_EV * 1000.0, RY_KCAL
    rep: dict = {"hamiltonian": {"functional": "PBE", "dispersion": "DFT-D3(BJ) (dftd3_version=4), two-body (dftd3_threebody=.false.)",
                                 "pseudopotentials": "SSSP 1.3.0 PBE precision: " + ", ".join(p for _, _, p in SPECIES),
                                 "ecutwfc_ry": 90, "ecutrho_ry": 360, "occupations": "fixed"},
                 "qe_versions": raw.get("qe_versions"), "pseudopotentials": raw.get("pseudopotentials"),
                 "units": "per monomer (CH2-CF2); each beta-type cell holds 2, alpha holds 4", "cells": {}, "margins": {}}
    for name, c in cells.items():
        if "n_monomers" not in c:
            continue
        k2 = c["stages"].get("scf_k2", {})
        rep["cells"][name] = {
            "n_atoms": c["n_atoms"], "n_monomers": c["n_monomers"], "final_relaxation": c["final_relaxation"],
            "lattice_A_deg": {k: round(v, 4) for k, v in c["lattice"].items()},
            "orientation": c["orientation"],
            "scf_k2": {k: k2.get(k) for k in ("energy_ry", "d3_ry", "mesh", "n_kpoints_irreducible", "sym_ops",
                                               "max_force_component_ry_bohr", "stress_kbar", "pressure_kbar", "meets_criteria")},
            "energies_ry_per_monomer": {s: _per_mon(c, s) for s in SCF_STAGES},
            "d3_ry_per_monomer": {s: _per_mon(c, s, "d3_ry") for s in SCF_STAGES},
        }
        if "relax_nosym" in c["stages"]:
            ns = c["stages"]["relax_nosym"]
            e_ns = _per_mon(c, "relax_nosym")
            e_k1 = _per_mon(c, "scf_k1")
            rep["cells"][name]["nosym_check"] = {
                "job_done": ns.get("job_done"), "bfgs_converged": ns.get("bfgs_converged"),
                "final_energy_ry_per_monomer": e_ns,
                "minus_symmetric_k1_meV_per_monomer": None if (e_ns is None or e_k1 is None) else (e_ns - e_k1) * meV,
                "final_lattice": ns.get("final_lattice"), "final_orientation": ns.get("final_orientation")}
    pol = cells.get("polar")

    def diff(a, b, stage, key="energy_ry"):
        ea, eb = _per_mon(cells[a], stage, key), _per_mon(cells[b], stage, key)
        if ea is None or eb is None:
            return None
        d = ea - eb
        return {"ry": d, "meV": d * meV, "kcal_mol": d * kcal}

    for anti in ("anti_free", "anti90"):
        if pol is None or anti not in cells or "n_monomers" not in cells[anti]:
            continue
        m = {s: diff(anti, "polar", s) for s in SCF_STAGES}
        d3 = {s: diff(anti, "polar", s, "d3_ry") for s in SCF_STAGES}
        rep["margins"][f"{anti} - polar"] = {"total": m, "d3_part": d3,
                                             "pbe_part": {s: None if m[s] is None or d3[s] is None else
                                                          {k: m[s][k] - d3[s][k] for k in m[s]} for s in SCF_STAGES}}
    if pol is not None and "anti_free" in cells:
        # the lowest energies found on the 4x2x8 mesh, symmetric or not: the margin's lower and upper ends
        def best(cell):
            es = [x for x in (_per_mon(cells[cell], "scf_k1"), _per_mon(cells[cell], "relax_nosym")) if x is not None]
            return min(es) if es else None
        lo = [best(c) for c in ("anti_free", "anti90") if c in cells and "n_monomers" in cells[c]]
        lo = [x for x in lo if x is not None]
        ep = best("polar")
        if lo and ep is not None:
            d = min(lo) - ep
            rep["margins"]["lowest antipolar found - lowest polar found, 4x2x8 (any relaxation, symmetric or not)"] = {
                "ry": d, "meV": d * meV, "kcal_mol": d * kcal}
    if pol is not None and "polar_dipole" in cells and "n_monomers" in cells["polar_dipole"]:
        rep["margins"]["polar_dipole - polar (same minimum?)"] = {s: diff("polar_dipole", "polar", s) for s in ("scf_k1", "scf_k2")}
    if pol is not None and "alpha" in cells and "n_monomers" in cells["alpha"]:
        rep["margins"]["alpha - polar (beta)"] = {"total": {s: diff("alpha", "polar", s) for s in SCF_STAGES},
                                                  "d3_part": {s: diff("alpha", "polar", s, "d3_ry") for s in SCF_STAGES}}
    # lattice against experiment
    lat_cmp = {}
    for name in ("polar", "polar_dipole"):
        if name in rep["cells"]:
            L = rep["cells"][name]["lattice_A_deg"]
            long_, short = max(L["a"], L["b"]), min(L["a"], L["b"])
            e = EXPERIMENT["beta"]
            lat_cmp[name] = {"long (a_exp 8.58)": [long_, round(100 * (long_ / e["a_long"] - 1), 2)],
                             "short, polar (b_exp 4.91)": [short, round(100 * (short / e["b_short_polar"] - 1), 2)],
                             "chain (c_exp 2.56)": [L["c"], round(100 * (L["c"] / e["c_chain"] - 1), 2)],
                             "angles": [L["alpha"], L["beta"], L["gamma"]], "volume_per_monomer_A3": round(L["volume"] / 2, 3)}
    if "alpha" in rep["cells"]:
        L = rep["cells"]["alpha"]["lattice_A_deg"]
        long_, short = max(L["a"], L["b"]), min(L["a"], L["b"])
        e = EXPERIMENT["alpha"]
        lat_cmp["alpha"] = {"short (a_exp 4.96)": [short, round(100 * (short / e["a_short"] - 1), 2)],
                            "long (b_exp 9.64)": [long_, round(100 * (long_ / e["b_long"] - 1), 2)],
                            "chain (c_exp 4.62)": [L["c"], round(100 * (L["c"] / e["c_chain"] - 1), 2)],
                            "angles": [L["alpha"], L["beta"], L["gamma"]], "volume_per_monomer_A3": round(L["volume"] / 4, 3)}
    rep["lattice_vs_experiment"] = {"experiment": EXPERIMENT, "dft": lat_cmp}
    rep["polyfind_margins_kcal_per_monomer"] = pf["polyfind_margins_kcal_per_monomer"]
    with open(out_json, "w", newline="\n") as f:
        json.dump(rep, f, indent=1)
    # a compact printout
    for k, v in rep["margins"].items():
        print(k)
        if "meV" in v:
            print(f"  {'':<14s} {v['meV']:+10.4f} meV/mon  {v['kcal_mol']:+9.5f} kcal/mol/mon")
            continue
        tot = v.get("total", v)
        for s, d in tot.items():
            if d:
                extra = ""
                if "d3_part" in v and v["d3_part"].get(s):
                    extra = f"   D3 part {v['d3_part'][s]['meV']:+9.4f} meV"
                print(f"  {s:<14s} {d['meV']:+10.4f} meV/mon  {d['kcal_mol']:+9.5f} kcal/mol/mon{extra}")
    return rep


# ------------------------------------------------------------------------------------------- main
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("starts", help="regenerate polyfind's cells and the provider's beta start")
    s.add_argument("--out", default=STARTS)
    s.add_argument("--sarco", default=SARCO_BETA)
    s.add_argument("--skip-polyfind", action="store_true")
    s.add_argument("--alpha", action="store_true", help="also build the alpha-PVDF control start")
    i = sub.add_parser("inputs", help="write the pw.x input of one stage of one cell")
    i.add_argument("--cell", required=True)
    i.add_argument("--stage", required=True)
    i.add_argument("--root", default=f"{WSL_ROOT}/runs")
    i.add_argument("--starts", default=STARTS)
    i.add_argument("--pseudo-dir", default=f"{WSL_ROOT}/pseudo")
    i.add_argument("--scratch", default=f"{WSL_ROOT}/scratch")
    k = sub.add_parser("check", help="exit 0 when a stage's last SCF meets the force/stress criteria")
    k.add_argument("--cell", required=True)
    k.add_argument("--stage", required=True)
    k.add_argument("--root", default=f"{WSL_ROOT}/runs")
    c = sub.add_parser("collect", help="parse every output under the run root into one JSON")
    c.add_argument("--root", default=f"{WSL_ROOT}/runs")
    c.add_argument("--starts", default=STARTS)
    c.add_argument("--json", default=os.path.join(OUT, "raw_stages.json"))
    c.add_argument("--pseudo-dir", default=f"{WSL_ROOT}/pseudo")
    c.add_argument("--out", default=OUT, help="where inputs/ and relaxed/ are written")
    r = sub.add_parser("report", help="margins, convergence, lattices from raw_stages.json -> polar_margin_dft.json")
    r.add_argument("--raw", default=os.path.join(OUT, "raw_stages.json"))
    r.add_argument("--starts", default=STARTS)
    r.add_argument("--json", default=os.path.join(OUT, "polar_margin_dft.json"))
    args = ap.parse_args(argv)
    if args.cmd == "starts":
        if not args.skip_polyfind:
            build_starts(args.out)
        print("provider start:", dipole_start(args.out, args.sarco))
        if args.alpha:
            build_alpha(args.out)
    elif args.cmd == "inputs":
        write_stage(args.root, args.cell, args.stage, args.starts, args.pseudo_dir, args.scratch)
    elif args.cmd == "check":
        ok, rep = check_stage(args.root, args.cell, args.stage)
        print(("PASS " if ok else "FAIL ") + rep)
        return 0 if ok else 1
    elif args.cmd == "collect":
        collect(args.root, args.json, args.starts, args.pseudo_dir, args.out)
    elif args.cmd == "report":
        report(args.raw, args.starts, args.json)
    return 0


if __name__ == "__main__":
    sys.exit(main())
