"""Can a pretrained machine-learned potential be a second opinion on the margins?

The test is the one margin this project holds from DFT: PVDF's all-trans (beta) minus its
helix (alpha), +0.726 kcal/mol per monomer in PBE-D3(BJ)
(`deliverables/polar_margin_dft/`). The switching window is 0.12 kcal/mol wide, so a second
opinion has to resolve about 0.1. A potential that misses this gap by more than that cannot
referee the candidates' margins, which are the same quantity in cells with pendants it knows
less well than PVDF.

Two readings per potential, both per monomer:

- single points on the three DFT-relaxed cells (alpha, polar beta, antipolar beta);
- each cell relaxed in the potential itself, atoms and lattice, from the DFT geometry, so the
  potential is not charged for a geometry it would not choose.

The potentials come in two kinds. The MACE-MP foundation models are fitted to PBE without
dispersion, so D3(BJ) is added the way the DFT anchor has it (two-body, PBE parameters); that
add-on alone is checked against the D3 part of the DFT margin. Models fitted to a
dispersion-corrected functional (Egret-1: wB97M-D3BJ) are given as files and run as they are;
their reference is then a different functional from the anchor's, which is a second reason a
difference can be there, and is said in the output.

Needs torch, mace-torch, ase and torch-dftd, none of which the package depends on:

    pip install torch --index-url https://download.pytorch.org/whl/cpu
    pip install mace-torch ase torch-dftd

    python examples/mlip_second_opinion.py --mace-mp medium-mpa-0,medium-0b3,medium \
        --model-file egret-1=/path/EGRET_1.model --relax

A model file is a pickled network: load one only from a source you trust.
"""

from __future__ import annotations

import argparse
import json
import time
import warnings
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RELAXED = ROOT / "deliverables" / "polar_margin_dft" / "relaxed"
ANCHOR = ROOT / "deliverables" / "polar_margin_dft" / "polar_margin_dft.json"
OUT = ROOT / "docs" / "predictions" / "mlip_second_opinion.json"

KCAL_PER_EV = 23.060548
CELLS = {
    "alpha": ("pvdf_alpha_pbe_d3bj.xyz", 4),
    "beta": ("pvdf_polar_pbe_d3bj.xyz", 2),
    "anti": ("pvdf_anti_free_pbe_d3bj.xyz", 2),
}
# PBE-D3(BJ), per monomer, from deliverables/polar_margin_dft (tests hold these to the file).
DFT = {
    "alpha_minus_beta_meV": -31.46,
    "alpha_minus_beta_pbe_part_meV": -51.55,
    "alpha_minus_beta_d3_part_meV": 20.09,
    "anti_minus_beta_meV": -1.11,
    "margin_kcal": 0.726,
}
NEEDED_KCAL = 0.1


def dft_cells():
    from ase.io import read

    return {k: (read(RELAXED / f), n) for k, (f, n) in CELLS.items()}


def d3_calculator():
    import torch
    from ase import units
    from torch_dftd.torch_dftd3_calculator import TorchDFTD3Calculator

    return TorchDFTD3Calculator(
        device="cpu", damping="bj", xc="pbe", dtype=torch.float64,
        cutoff=40.0 * units.Bohr, cnthr=40.0 * units.Bohr, abc=False,
    )


def single_points(calc, cells):
    from ase import units

    out = {}
    for k, (atoms, n) in cells.items():
        a = atoms.copy()
        a.calc = calc
        out[k] = {
            "energy_eV_per_monomer": float(a.get_potential_energy() / n),
            "largest_force_eV_per_A": float(np.abs(a.get_forces()).max()),
        }
        try:
            s = a.get_stress(voigt=True) / units.GPa
            out[k]["stress_xx_yy_zz_GPa"] = [float(v) for v in s[:3]]
        except Exception as err:  # a model without a stress head
            out[k]["stress_xx_yy_zz_GPa"] = None
            out[k]["stress_error"] = str(err)[:120]
    return out


def relaxed(calc, cells, fmax, steps, keep=None):
    """Each cell relaxed in the potential, atoms and lattice, from the DFT geometry.

    `keep`, a dict, receives the relaxed cells themselves (for a DFT single point on them).
    """
    from ase.filters import FrechetCellFilter
    from ase.optimize import BFGS

    out = {}
    for k, (atoms, n) in cells.items():
        a = atoms.copy()
        a.calc = calc
        try:
            opt = BFGS(FrechetCellFilter(a), logfile=None)
            done = bool(opt.run(fmax=fmax, steps=steps))
            cell0, cell1 = atoms.cell.cellpar(), a.cell.cellpar()
            out[k] = {
                "energy_eV_per_monomer": float(a.get_potential_energy() / n),
                "converged": done,
                "steps": int(opt.get_number_of_steps()),
                "largest_force_eV_per_A": float(np.abs(a.get_forces()).max()),
                "lattice_A_deg": [float(v) for v in cell1],
                "lattice_change_percent": [float(100 * (b / c - 1)) for b, c in zip(cell1[:3], cell0[:3])],
                "volume_change_percent": float(100 * (a.get_volume() / atoms.get_volume() - 1)),
            }
            if keep is not None:
                b = a.copy()
                b.calc = None
                keep[k] = b
        except Exception as err:
            out[k] = {"error": str(err)[:200]}
    return out


def gaps(points, add=None):
    """alpha - beta and anti - beta in meV per monomer, and the margin in kcal/mol."""
    if any("energy_eV_per_monomer" not in points[k] for k in CELLS):
        return None
    e = {k: points[k]["energy_eV_per_monomer"] + (add[k] if add else 0.0) for k in CELLS}
    ab = 1e3 * (e["alpha"] - e["beta"])
    margin = -ab * 1e-3 * KCAL_PER_EV
    return {
        "alpha_minus_beta_meV": ab,
        "anti_minus_beta_meV": 1e3 * (e["anti"] - e["beta"]),
        "margin_kcal": margin,
        "margin_error_kcal": margin - DFT["margin_kcal"],
        "times_the_needed_accuracy": abs(margin - DFT["margin_kcal"]) / NEEDED_KCAL,
    }


def run(models, relax, fmax, steps, threads, keep=None):
    import torch
    from ase.calculators.mixing import SumCalculator

    torch.set_num_threads(threads)
    cells = dft_cells()
    d3 = d3_calculator()
    d3_points = single_points(d3, cells)
    d3_e = {k: d3_points[k]["energy_eV_per_monomer"] for k in CELLS}
    d3_gap = 1e3 * (d3_e["alpha"] - d3_e["beta"])
    out = {
        "dft": DFT,
        "needed_accuracy_kcal": NEEDED_KCAL,
        "d3_bj_alone": {
            "alpha_minus_beta_meV": d3_gap,
            "dft_d3_part_meV": DFT["alpha_minus_beta_d3_part_meV"],
        },
        "models": {},
    }
    print(f"D3(BJ) alone: alpha - beta {d3_gap:+.2f} meV per monomer (the DFT margin's D3 part: "
          f"{DFT['alpha_minus_beta_d3_part_meV']:+.2f})", flush=True)
    for name, kind, source in models:
        t0 = time.time()
        if kind == "mace-mp":
            from mace.calculators import mace_mp

            bare = mace_mp(model=source, dispersion=False, default_dtype="float64", device="cpu")
            calc, add = SumCalculator([bare, d3]), d3_e
            fitted_to = "PBE without dispersion; D3(BJ) added as in the anchor"
        elif kind == "file+d3":
            from mace.calculators import MACECalculator

            bare = MACECalculator(model_paths=str(source), device="cpu", default_dtype="float64")
            calc, add = SumCalculator([bare, d3]), d3_e
            fitted_to = "the anchor's recipe less its D3(BJ) part; D3(BJ) added back"
        else:
            from mace.calculators import MACECalculator

            bare = MACECalculator(model_paths=str(source), device="cpu", default_dtype="float64")
            calc, add = bare, None
            fitted_to = "its own dispersion-corrected functional, not the anchor's PBE-D3(BJ)"
        rec = {"kind": kind, "fitted_to": fitted_to}
        try:
            rec["elements"] = sorted(int(z) for z in bare.models[0].atomic_numbers)
            rec["cutoff_A"] = float(bare.models[0].r_max)
        except Exception:
            pass
        points = single_points(bare, cells)
        rec["on_the_dft_geometries"] = {"cells": points, **(gaps(points, add) or {})}
        g = rec["on_the_dft_geometries"]
        line = (f"{name:16s} on the DFT geometries: margin {g['margin_kcal']:+.2f} kcal/mol "
                f"(DFT {DFT['margin_kcal']:+.3f}), off by {g['margin_error_kcal']:+.2f} = "
                f"{g['times_the_needed_accuracy']:.0f}x the needed {NEEDED_KCAL}; anti - beta "
                f"{g['anti_minus_beta_meV']:+.1f} meV (DFT {DFT['anti_minus_beta_meV']:+.2f}); largest force "
                f"{max(points[k]['largest_force_eV_per_A'] for k in CELLS):.2f} eV/A")
        if relax:
            own = relaxed(calc, cells, fmax, steps, keep)
            rec["relaxed_in_the_potential"] = {"cells": own, **(gaps(own) or {})}
            r = rec["relaxed_in_the_potential"]
            if "margin_kcal" in r:
                vol = ", ".join(f"{k} {own[k]['volume_change_percent']:+.1f}%" for k in CELLS)
                line += (f"\n{'':16s} relaxed in the potential: margin {r['margin_kcal']:+.2f}, off by "
                         f"{r['margin_error_kcal']:+.2f} = {r['times_the_needed_accuracy']:.0f}x; anti - beta "
                         f"{r['anti_minus_beta_meV']:+.1f} meV; volume against DFT's: {vol}; converged: "
                         + ", ".join(f"{k} {own[k]['converged']}" for k in CELLS))
            else:
                line += f"\n{'':16s} relaxed in the potential: failed, " + "; ".join(
                    f"{k}: {own[k].get('error', 'ok')}" for k in CELLS)
        rec["seconds"] = time.time() - t0
        out["models"][name] = rec
        print(line + f"  ({rec['seconds']:.0f} s)", flush=True)
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--mace-mp", default="", help="MACE-MP model names, comma-separated (PBE; D3(BJ) is added)")
    p.add_argument("--model-file", action="append", default=[], metavar="NAME=PATH",
                   help="a MACE-architecture model file with dispersion in its reference (run as it is)")
    p.add_argument("--model-file-d3", action="append", default=[], metavar="NAME=PATH",
                   help="a model file fitted to the anchor's recipe less D3(BJ) (ours); D3(BJ) is added back")
    p.add_argument("--relax", action="store_true", help="also relax each cell in the potential")
    p.add_argument("--fmax", type=float, default=0.01)
    p.add_argument("--steps", type=int, default=400)
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--out", type=Path, default=OUT)
    p.add_argument("--merge", action="store_true", help="keep the models already in the output file")
    args = p.parse_args()
    warnings.filterwarnings("ignore")
    models = [(m, "mace-mp", m) for m in args.mace_mp.split(",") if m]
    for kind, items in (("file", args.model_file), ("file+d3", args.model_file_d3)):
        for item in items:
            name, path = item.split("=", 1)
            models.append((name, kind, Path(path)))
    out = run(models, args.relax, args.fmax, args.steps, args.threads)
    if args.merge and args.out.exists():
        old = json.loads(args.out.read_text())
        old["models"].update(out["models"])
        out["models"] = old["models"]
    args.out.write_text(json.dumps(out, indent=1) + "\n")
    print(f"wrote {args.out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
