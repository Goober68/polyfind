#!/usr/bin/env python3
"""A path between two states of one cell in a potential: a kink sheet forming from the all-trans cell.

The two states are frames of the same atoms in the same order (an all-trans cell and its kink
sheet, ``docs/predictions/mlip/cells14_frames_pbe_d3bj.extxyz``).  They differ in a few backbone
torsions (a kink pair is two of them a chain) and in the cell: the sheet's cell is sheared along
the chains.  A straight line between them in Cartesian or fractional coordinates turns a torsion
through its bonds, so the path is made the way the molecule would go:

``scan``   the torsions that differ are led from the start's values to the end's in equal steps by a
           restraint, and at each step everything else is relaxed in the potential, the cell and
           its shear included.  The images are the path; the energy stored with each is the
           potential's with D3 and without the restraint.  The last image is then released and
           relaxed, which says where the path arrived.
           ``--refine SCAN`` does the same in finer steps between the two images beside a scan's top.
``choose`` of a path's scans (one from each end), the one that arrives where it was led and has the
           lower top above the all-trans end: its top image and the two beside it, as the names of
           frames for a DFT single point each (``examples/gpu_jobs/dft_cells.sh``).
``release`` images of a scan relaxed with no restraint: an image before the top should fall back
           to the start and one after it on to the end, or the top is not what separates them.
``table``  scans (and the DFT labels of their images, where there are any) as one table.

A nudged elastic band is not used: ase's refuses images whose cells differ, and here the cell's
shear is part of the path.

Energies are kcal/mol per chain of the cell (two chains a cell) against the scan's first image.
What a path in a two-chain periodic cell is: every chain of the crystal forms its kink in step
with its neighbours.  It is the barrier of the sheet forming as a whole, per chain, and not the
barrier of one chain's kink forming among straight neighbours.

docs/MLIP_SECOND_OPINION.md 3.23.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import re
import warnings
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
KCAL = 23.0605
CHAINS = 2  # the rows' cells hold two chains


def _example(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "examples" / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def backbone_torsions(atoms, reach=1.85):
    """The backbone's torsions as index quadruples, chain by chain, in the order of the atoms.

    A backbone carbon has two carbon neighbours along the chain; a pendant carbon (a nitrile's) has
    one and is taken off first.  A chain closes on itself through the cell, so each has as many
    torsions as backbone carbons."""
    from ase.neighborlist import neighbor_list

    i, j = neighbor_list("ij", atoms, {("C", "C"): reach})
    nb = {}
    for x, y in zip(i, j):
        nb.setdefault(int(x), set()).add(int(y))
    back = set(nb)
    while True:
        off = {k for k in back if len(nb[k] & back) < 2}
        if not off:
            break
        back -= off
    if any(len(nb[k] & back) != 2 for k in back):
        raise SystemExit("the backbone is not a set of plain chains: a carbon has more than two backbone neighbours")
    quads, seen = [], set()
    for first in sorted(back):
        if first in seen:
            continue
        chain = [first, min(nb[first] & back)]   # the walk starts at the chain's lowest atom, toward its lower neighbour
        while True:
            nxt = (nb[chain[-1]] & back) - {chain[-2]}
            (y,) = nxt
            if y == first:
                break
            chain.append(y)
        seen.update(chain)
        n = len(chain)
        quads += [tuple(chain[(t + q) % n] for q in range(4)) for t in range(n)]
    return quads


def _unwrapped(atoms, quad):
    """The four atoms of a torsion as one piece: each placed at its nearest image to the one before."""
    from ase.geometry import find_mic

    p = [atoms.positions[quad[0]].copy()]
    for a, b in zip(quad[:-1], quad[1:]):
        step, _ = find_mic(atoms.positions[b] - atoms.positions[a], atoms.cell, pbc=True)
        p.append(p[-1] + step)
    return np.array(p)


def _torsion(p, xp=np):
    """The torsion angle of four points, radians in (-pi, pi]; trans is pi.  ``xp`` is numpy or torch."""
    b0, b1, b2 = p[0] - p[1], p[2] - p[1], p[3] - p[2]
    n = b1 / (b1 * b1).sum() ** 0.5
    v, w = b0 - (b0 * n).sum() * n, b2 - (b2 * n).sum() * n
    cross = xp.stack([n[1] * v[2] - n[2] * v[1], n[2] * v[0] - n[0] * v[2], n[0] * v[1] - n[1] * v[0]])
    return (np.arctan2 if xp is np else xp.atan2)((cross * w).sum(), (v * w).sum())


def torsions_deg(atoms, quads):
    return np.array([np.degrees(_torsion(_unwrapped(atoms, q))) for q in quads])


def arc(a, b):
    """From angle a to angle b the short way, degrees."""
    return (np.asarray(b) - np.asarray(a) + 180.0) % 360.0 - 180.0


def restraint_calculator(quads, targets_deg, k):
    """E = k * sum(1 - cos(phi - phi0)) over the torsions given (eV; k in eV, so k per rad2 near phi0),
    with its forces and its stress: a torsion changes when the cell is sheared."""
    import torch
    from ase.calculators.calculator import Calculator, all_changes

    class Restraint(Calculator):
        implemented_properties = ["energy", "forces", "stress"]

        def __init__(self):
            super().__init__()
            self.targets = np.radians(np.asarray(targets_deg, dtype=float))

        def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
            super().calculate(atoms, properties, system_changes)
            forces, virial, energy = np.zeros((len(atoms), 3)), np.zeros((3, 3)), 0.0
            for quad, phi0 in zip(quads, self.targets):
                p = torch.tensor(_unwrapped(atoms, quad), dtype=torch.float64, requires_grad=True)
                e = k * (1.0 - torch.cos(_torsion(p, torch) - phi0))
                (g,) = torch.autograd.grad(e, p)
                g = g.numpy()
                energy += float(e.detach())
                for a, ga, pa in zip(quad, g, p.detach().numpy()):
                    forces[a] -= ga
                    virial += np.outer(ga, pa)   # dE/d(strain): the four atoms move with the cell as one piece
            stress = 0.5 * (virial + virial.T) / atoms.get_volume()
            self.results = {"energy": energy, "forces": forces, "stress": stress[[0, 1, 2, 1, 0, 0], [0, 1, 2, 2, 2, 1]]}

    return Restraint()


def _potential(model):
    from ase.calculators.mixing import SumCalculator

    so = _example("mlip_second_opinion")
    return SumCalculator([so.file_calculator(model), so.d3_calculator()])


def _frame(files, name):
    from ase.io import read

    for f in files:
        for a in read(f, index=":"):
            if str(a.info.get("frame")) == name:
                return a
    raise SystemExit(f"no frame {name} in {[str(f) for f in files]}")


def _shear(atoms):
    """The cell's first two vectors along the third, in monomers (row_cells_frames.shear_along_chains)."""
    c = atoms.cell.array
    along = c[2] / np.linalg.norm(c[2])
    per = np.linalg.norm(c[2]) / (int(atoms.info["n_monomers"]) / CHAINS)
    return [float(c[0] @ along / per), float(c[1] @ along / per)]


def _bare(atoms, start, name, note, **more):
    b = atoms.copy()
    b.calc = None
    b.info = {k: v for k, v in start.info.items() if k in ("n_monomers", "kmesh", "family", "conformation", "source")}
    b.info.update(frame=name, family="path", note=note, **more)
    b.arrays = {k: v for k, v in b.arrays.items() if k in ("numbers", "positions")}
    return b


def scan(args):
    import torch
    from ase.calculators.mixing import SumCalculator
    from ase.filters import FrechetCellFilter
    from ase.io import write
    from ase.optimize import BFGS

    warnings.filterwarnings("ignore")
    torch.set_num_threads(args.threads)
    parent = json.loads(args.refine.read_text()) if args.refine else None
    if parent:   # the same two frames and the same torsions as the scan refined
        args.start, args.end = parent["start"], parent["end"]
    start, end = _frame(args.frames, args.start), _frame(args.frames, args.end)
    if (start.numbers != end.numbers).any():
        raise SystemExit("the two frames are not the same atoms in the same order")
    quads = backbone_torsions(start)
    if backbone_torsions(end) != quads:
        raise SystemExit("the two frames do not have the same backbone")
    a0, a1 = torsions_deg(start, quads), torsions_deg(end, quads)
    turn = arc(a0, a1)
    led = [q for q in range(len(quads)) if abs(turn[q]) > args.differ]
    if not led:
        raise SystemExit(f"no backbone torsion differs by more than {args.differ} degrees between the two frames")
    print(f"{args.name}: {len(quads)} backbone torsions, {len(led)} led: " + "; ".join(f"{quads[q]} {a0[q]:+.0f} to {a1[q]:+.0f}" for q in led), flush=True)
    base = _potential(args.model)
    n = int(start.info["n_monomers"])
    a = start.copy()
    a.calc = base
    e_start = a.get_potential_energy()
    lams = [point / args.points for point in range(1, args.points + 1)]
    lam0 = 0.0
    if parent:   # from the image before the scan's top to the image after it, in this scan's steps
        at, last = parent["top"]["index"], len(parent["images"]) - 2
        lo, hi = max(0, at - 1), min(last, at + 1)
        lam0, lam1 = parent["images"][lo]["lambda"], parent["images"][hi]["lambda"]
        lams = [lam0 + (lam1 - lam0) * point / args.points for point in range(1, args.points + 1)]
        e_start = parent["images"][0]["potential_eV_per_monomer"] * n   # energies stay against the path's first image
        a = _frame([args.refine.with_suffix(".extxyz")], f"{parent['name']}_i{lo:02d}")
        a.info = dict(start.info)
        a.calc = base
    b = end.copy()
    b.calc = base
    e_end = b.get_potential_energy()
    images, rows = [], []

    def keep(label, lam, steps, done, restraint_energy, note):
        a.calc = base
        e = a.get_potential_energy()
        now = torsions_deg(a, quads)
        row = {"image": label, "lambda": lam, "energy": (e - e_start) / CHAINS * KCAL, "potential_eV_per_monomer": e / n, "steps": steps, "converged": bool(done),
               "restraint_energy": restraint_energy / CHAINS * KCAL, "led_torsions_deg": [float(now[q]) for q in led],
               "largest_other_torsion_from_trans_deg": float(max(180.0 - abs(now[q]) for q in range(len(quads)) if q not in led)),
               "largest_force_eV_per_A": float(np.abs(a.get_forces()).max()), "shear_along_chains_monomers": _shear(a), "cell_A_deg": [float(v) for v in a.cell.cellpar()],
               "volume_A3_per_monomer": float(a.get_volume() / n)}
        rows.append(row)
        images.append(_bare(a, start, f"{args.name}_{label}", note, potential_eV_per_monomer=e / n, path=args.name, path_lambda=lam, potential=args.potential))
        write(args.out / f"{args.name}.extxyz", images, format="extxyz")
        (args.out / f"{args.name}.json").write_text(json.dumps(table, indent=1) + "\n")
        print(f"{args.name}_{label}: {row['energy']:+.3f} kcal/mol per chain, torsions " + " ".join(f"{v:+.0f}" for v in row["led_torsions_deg"])
              + f", shear {row['shear_along_chains_monomers'][0]:+.2f} {row['shear_along_chains_monomers'][1]:+.2f}, {steps} steps{'' if done else ' NOT converged'}", flush=True)

    table = {"what": "a scan from one labelled cell to another in a potential with D3: the backbone torsions that differ are led by a restraint in equal steps and "
                     "everything else, the cell included, is relaxed at each step. energy: kcal/mol per chain against the first image, without the restraint. The last "
                     "image, released, says where the path arrived (end_in_the_potential is the end frame's own energy there, not relaxed)",
             "name": args.name, "potential": args.potential, "model": Path(args.model).name, "start": args.start, "end": args.end, "points": args.points,
             "restraint_eV": args.k, "fmax_eV_per_A": args.fmax, "led_torsions": [list(quads[q]) for q in led], "start_torsions_deg": [float(a0[q]) for q in led],
             "end_torsions_deg": [float(a1[q]) for q in led], "end_in_the_potential": (e_end - e_start) / CHAINS * KCAL,
             "end_shear_along_chains_monomers": _shear(end), "images": rows}
    args.out.mkdir(parents=True, exist_ok=True)
    if parent:
        table["refines"] = {"scan": parent["name"], "from_image": f"i{lo:02d}", "to_image": f"i{hi:02d}", "lambda": [lam0, lam1]}
        keep("i00", lam0, 0, True, 0.0, f"image i{lo:02d} of {parent['name']}: the first image of its refinement")
    else:
        keep("i00", 0.0, 0, True, 0.0, f"{args.start}: the scan's first image, as labelled")
    for point, lam in enumerate(lams, start=1):
        restraint = restraint_calculator([quads[q] for q in led], [a0[q] + lam * turn[q] for q in led], args.k)
        a.calc = SumCalculator([base, restraint])
        opt = BFGS(FrechetCellFilter(a), logfile=None, maxstep=args.maxstep)
        done = opt.run(fmax=args.fmax, steps=args.steps)
        held = restraint.get_potential_energy(a)
        keep(f"i{point:02d}", lam, opt.get_number_of_steps(), done, held,
             f"{args.start} toward {args.end}, {len(led)} torsions led {100 * lam:.0f}% of the way by a restraint of {args.k} eV, the rest relaxed in {args.potential} with D3")
    on_the_scan = len(rows)
    if parent:
        table["arrived"] = parent.get("arrived")
    else:
        a.calc = base   # released: where did it arrive?
        opt = BFGS(FrechetCellFilter(a), logfile=None, maxstep=args.maxstep)
        done = opt.run(fmax=args.fmax, steps=args.release_steps)
        keep("released", 1.0, opt.get_number_of_steps(), done, 0.0, f"the last image of {args.name} released and relaxed in {args.potential} with D3")
        arrived = rows[-1]
        table["arrived"] = {"energy": arrived["energy"], "end_in_the_potential": table["end_in_the_potential"], "converged": arrived["converged"],
                            "torsions_from_the_end_deg": [float(v) for v in arc(a1[led], np.array(arrived["led_torsions_deg"]))],
                            "shear_along_chains_monomers": arrived["shear_along_chains_monomers"]}
    top = max(range(on_the_scan), key=lambda r: rows[r]["energy"])
    table["top"] = {"image": rows[top]["image"], "energy": rows[top]["energy"], "index": top}
    (args.out / f"{args.name}.json").write_text(json.dumps(table, indent=1) + "\n")
    print(f"{args.name}: top {rows[top]['energy']:+.3f} at {rows[top]['image']}"
          + (f"; arrived at {table['arrived']['energy']:+.3f} (the end frame itself {table['end_in_the_potential']:+.3f})" if table["arrived"] else ""), flush=True)


def _on_the_scan(table):
    return [r for r in table["images"] if r["image"] != "released"]


def _above_all_trans(table):
    """A scan's top against its all-trans end, and whether the scan arrived where it was led: the end frame's torsions,
    at an energy within a kcal/mol per chain of the end frame's own in the potential (which is not relaxed in it)."""
    arrived, top = table["arrived"], table["top"]["energy"]
    there = max(abs(v) for v in arrived["torsions_from_the_end_deg"]) < 30.0 and abs(arrived["energy"] - arrived["end_in_the_potential"]) < 1.0
    flat_first = max(180.0 - abs(v) for v in table["start_torsions_deg"]) < 30.0   # the scan starts at the all-trans cell
    return (top if flat_first else top - arrived["energy"]), bool(there)


def choose(args):
    """Of a path's scans, the one that arrived and has the lower top above all-trans; printed: the scan's file, then the
    frames to label, the top first, then the image before it and the image after it (a first image is labelled already)."""
    scans = {f: json.loads(Path(f).read_text()) for f in args.scans}
    marks = {f: _above_all_trans(t) for f, t in scans.items()}
    arrived = [f for f in scans if marks[f][1]] or list(scans)
    best = min(arrived, key=lambda f: marks[f][0])
    rows, name = _on_the_scan(scans[best]), scans[best]["name"]
    at = scans[best]["top"]["index"]
    names = [rows[at]["image"]] + [rows[r]["image"] for r in (at - 1, at + 1) if 0 < r < len(rows)]
    print(" ".join([best] + [f"{name}_{im}" for im in names if im != "i00"]))


def release(args):
    """Images of a scan relaxed with no restraint: where each falls."""
    import torch
    from ase.filters import FrechetCellFilter
    from ase.io import write
    from ase.optimize import BFGS

    warnings.filterwarnings("ignore")
    torch.set_num_threads(args.threads)
    table = json.loads(args.scan.read_text())
    name, quads = table["name"], [tuple(q) for q in table["led_torsions"]]
    base = _potential(args.model)
    images = args.images
    if not images:   # the image before the top and the image after it
        rows, at = _on_the_scan(table), table["top"]["index"]
        images = [rows[r]["image"] for r in (at - 1, at + 1) if 0 <= r < len(rows)]
    out = {"what": "images of a scan relaxed in the potential with no restraint: an image before the top should fall back to the scan's start and one after it on "
                   "to its end. energy: kcal/mol per chain against the path's first image", "scan": name, "start_torsions_deg": table["start_torsions_deg"],
           "end_torsions_deg": table["end_torsions_deg"], "top": table["top"], "images": {}}
    kept = []
    for im in images:
        a = _frame([args.scan.with_suffix(".extxyz")], f"{name}_{im}")
        n = int(a.info["n_monomers"])
        e_first = table["images"][0]["potential_eV_per_monomer"] * n
        a.calc = base
        e0 = a.get_potential_energy()
        opt = BFGS(FrechetCellFilter(a), logfile=None, maxstep=0.1)
        done = opt.run(fmax=args.fmax, steps=args.steps)
        e, now = a.get_potential_energy(), torsions_deg(a, quads)
        to_start = float(np.abs(arc(table["start_torsions_deg"], now)).max())
        to_end = float(np.abs(arc(table["end_torsions_deg"], now)).max())
        out["images"][im] = {"from": (e0 - e_first) / CHAINS * KCAL, "energy": (e - e_first) / CHAINS * KCAL, "steps": opt.get_number_of_steps(), "converged": bool(done),
                             "led_torsions_deg": [float(v) for v in now], "fell_to": "start" if to_start < 30.0 else "end" if to_end < 30.0 else "neither",
                             "shear_along_chains_monomers": _shear(a)}
        kept.append(_bare(a, a, f"{name}_{im}_free", f"image {im} of {name} relaxed with no restraint in {table['potential']} with D3 ({opt.get_number_of_steps()} steps"
                          f"{'' if done else ', NOT converged'})", potential_eV_per_monomer=e / n, path=name, potential=table["potential"]))
        args.out.mkdir(parents=True, exist_ok=True)
        write(args.out / f"{name}_free.extxyz", kept, format="extxyz")
        (args.out / f"{name}_free.json").write_text(json.dumps(out, indent=1) + "\n")
        r = out["images"][im]
        print(f"{name}_{im}: from {r['from']:+.3f} to {r['energy']:+.3f} kcal/mol per chain in {r['steps']} steps{'' if done else ' NOT converged'}: fell to the {r['fell_to']}", flush=True)


def table(args):
    """Scans, bands and the DFT labels of their images as one table (kcal/mol per chain against each scan's first image)."""
    from ase.io import read

    labelled = {}
    for f in args.labels:
        for a in read(f, index=":"):
            labelled[str(a.info["frame"])] = a
    out = {"what": "paths in a potential between labelled cells (examples/mlip_path.py scan), kcal/mol per chain of the two-chain cell against each path's first image; "
                   "dft: single points on the potential's images, against DFT's label of the first image. forward: the top against the first image; back: against the "
                   "released last image", "paths": {}}
    for f in args.scans:
        t = json.loads(Path(f).read_text())
        name, n = t["name"], None
        first = labelled.get(t["start"])
        rows = []
        for r in t["images"]:
            row = {k: r[k] for k in ("image", "lambda", "energy", "led_torsions_deg", "shear_along_chains_monomers", "converged", "steps", "volume_A3_per_monomer")}
            lab = labelled.get(f"{name}_{r['image']}")
            if r["image"] == "i00" and first is not None and "refines" not in t:
                row["dft"] = 0.0
            elif lab is not None and first is not None:
                n = int(lab.info["n_monomers"])
                row["dft"] = (float(lab.info["REF_energy"]) - float(first.info["REF_energy"])) / CHAINS * KCAL
                row["dft_largest_force_eV_per_A"] = float(np.abs(lab.arrays["REF_forces"]).max())
                row["potential_minus_dft"] = row["energy"] - row["dft"]
            rows.append(row)
        rec = {"potential": t["potential"], "start": t["start"], "end": t["end"], "led_torsions": t["led_torsions"], "points": t["points"],
               "top": t.get("top"), "arrived": t.get("arrived"), "end_in_the_potential": t["end_in_the_potential"], "images": rows}
        if t.get("top") and t.get("arrived"):
            rec["forward"] = t["top"]["energy"]
            rec["back"] = t["top"]["energy"] - t["arrived"]["energy"]
            rec["above_all_trans"], rec["arrived_where_led"] = _above_all_trans(t)
        with_dft = [r for r in rows if "dft" in r and r["image"] not in ("i00", "released")]
        if with_dft:
            hi = max(with_dft, key=lambda r: r["dft"])
            rec["dft_top"] = {"image": hi["image"], "energy": hi["dft"], "potential_there": hi["energy"], "images_labelled": [r["image"] for r in with_dft]}
        if "refines" in t:
            rec["refines"] = t["refines"]
        free = Path(f).with_name(f"{name}_free.json")
        if free.exists():
            rec["released"] = json.loads(free.read_text())["images"]
        out["paths"][name] = rec
        print(f"{name}: top {rec.get('forward', float('nan')):+.2f} above its first image, {rec.get('back', float('nan')):+.2f} above where it arrived"
              + (f"; DFT at the top image labelled {rec['dft_top']['energy']:+.2f} ({rec['dft_top']['image']}, the potential {rec['dft_top']['potential_there']:+.2f})" if with_dft else ""))
    Path(args.out).write_text(json.dumps(out, indent=1) + "\n")
    print(f"written {args.out}")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="step", required=True)
    s = sub.add_parser("scan")
    s.add_argument("--model", type=Path, required=True)
    s.add_argument("--potential", required=True, help="a name for the potential, kept with the images")
    s.add_argument("--frames", type=Path, nargs="+", required=True, help="files that hold the two frames")
    s.add_argument("--start", default="", help="the frame the scan starts from")
    s.add_argument("--end", default="", help="the frame whose torsions it is led to")
    s.add_argument("--refine", type=Path, help="a scan's .json: this scan goes from the image before that scan's top to the image after it, "
                   "between the same two frames (--start and --end are taken from it)")
    s.add_argument("--name", required=True, help="the path's name: its images are NAME_i00, NAME_i01, ...")
    s.add_argument("--out", type=Path, required=True)
    s.add_argument("--points", type=int, default=12)
    s.add_argument("--k", type=float, default=10.0, help="the restraint, eV (per rad2 near its target)")
    s.add_argument("--differ", type=float, default=60.0, help="degrees: a torsion that differs by more between the two frames is led")
    s.add_argument("--fmax", type=float, default=0.02)
    s.add_argument("--steps", type=int, default=400, help="the most steps at one point")
    s.add_argument("--release-steps", type=int, default=800)
    s.add_argument("--maxstep", type=float, default=0.1)
    s.add_argument("--threads", type=int, default=1)
    s = sub.add_parser("choose")
    s.add_argument("--scans", nargs="+", required=True, help="the .json of a path's scans, one from each end")
    s = sub.add_parser("release")
    s.add_argument("--model", type=Path, required=True)
    s.add_argument("--scan", type=Path, required=True, help="a scan's .json (its images beside it)")
    s.add_argument("--images", nargs="*", default=[], help="images to release, e.g. i05 i07 (default: the one before the top and the one after it)")
    s.add_argument("--out", type=Path, required=True)
    s.add_argument("--fmax", type=float, default=0.02)
    s.add_argument("--steps", type=int, default=600)
    s.add_argument("--threads", type=int, default=1)
    s = sub.add_parser("table")
    s.add_argument("--scans", nargs="+", required=True, help="scans' .json files")
    s.add_argument("--labels", nargs="*", default=[], help="files of DFT-labelled frames: the paths' first images and whichever images were labelled")
    s.add_argument("--out", required=True)
    args = p.parse_args()
    {"scan": scan, "choose": choose, "release": release, "table": table}[args.step](args)


if __name__ == "__main__":
    main()
