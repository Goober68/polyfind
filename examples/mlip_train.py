"""Fit a potential of our own to the DFT frames, and test it against the anchor.

The frames are `examples/mlip_frames.py`'s: cells with the energy, forces and stress of the
anchor's recipe (PBE with two-body D3(BJ)). Three steps:

    python examples/mlip_train.py prepare  --frames docs/predictions/mlip/pvdf_frames_pbe_d3bj.extxyz [MORE ...] --data DIR
    python examples/mlip_train.py train    --data DIR --name NAME [--foundation medium-0b3 | --scratch]
    python examples/mlip_train.py evaluate --data DIR --model DIR/NAME.model --name NAME

`prepare` takes the D3(BJ) part out of every label, so the network is fitted to the PBE part
alone and D3 is added back analytically when the potential is used. That is how the pretrained
MACE-MP models are built, and the reason is the same: the network sees 6 A, the dispersion tail
reaches further, and the add-on reproduces the DFT margin's dispersion part to 0.04 meV
(`examples/mlip_second_opinion.py`). The D3 energy the add-on gives for each frame is checked
against the one Quantum ESPRESSO printed for it: they differ by the tail between the add-on's
40 Bohr and the program's 95, which is 0.11 meV per atom, goes as the density and nothing else
(1.01 meV A3 per atom, the same in every frame to 2%), and is left in the labels for the
network. The four relaxed cells (`end`) are held out
of the fit whole; one frame in nine of the rest is held out for the force error.

`train` runs MACE's own training program. With `--foundation` the fit starts from a pretrained
MACE-MP model (MIT licence) and every weight is free; the reference energy of each element is
set so that the pretrained model's interaction energies match ours on average (the program's
`estimated`), which is what lets energies from a different code be fitted on. With `--scratch`
a small network starts from nothing.

`evaluate` asks the four things `docs/MLIP_SECOND_OPINION.md` 3.3 lists: beta minus alpha
against +0.726 kcal/mol per monomer, antipolar minus polar beta against -1.1 meV (both on the
DFT geometries and relaxed in the potential), beta-PVDF's elastic constants against the DFT
ones, and the errors on the held-out frames.

Needs torch, mace-torch, ase, torch-dftd (see `examples/mlip_second_opinion.py`).
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import subprocess
import sys
import time
import warnings
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "docs" / "predictions"
# Dipole's six +/-1% cells, relaxed ions (docs/CRYSTAL_STAGE.md); C_33 from the same record.
DFT_ELASTIC_GPA = {"C_11": 22.77, "C_22": 24.52, "C_12": 5.88, "C_33": 302.6}


def _example(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "examples" / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _second_opinion():
    return _example("mlip_second_opinion")


def prepare(args):
    from ase.io import read, write
    from ase.units import create_units

    ry = create_units("2006")["Ry"]  # the constant ase's reader used for the labels
    d3 = _second_opinion().d3_calculator()
    frames = [a for f in args.frames for a in read(f, index=":")]
    tail, train, test = [], [], []
    for i, a in enumerate(frames):
        b = a.copy()
        b.calc = d3
        e, f, s = b.get_potential_energy(), b.get_forces(), b.get_stress(voigt=False)
        # the add-on stops at 40 Bohr, the DFT program's D3 at 95: what is left is the far tail,
        # an energy per atom that goes as the density and nothing else (meV A3)
        tail.append(1e3 * (e - a.info["dft_d3_ry"] * ry) / len(a) * a.get_volume() / len(a))
        out = a.copy()
        out.calc = None
        for key in ("REF_energy", "REF_stress"):
            out.info.pop(key, None)
        out.info["REF_energy"] = float(a.info["REF_energy"] - e)
        out.info["REF_stress"] = np.asarray(a.info["REF_stress"]).reshape(3, 3) - s
        out.arrays["REF_forces"] = a.arrays["REF_forces"] - f
        out.info["d3_energy"] = float(e)
        (test if a.info["family"] == "end" or i % 9 == 5 else train).append(out)
    tail = np.array(tail)
    if tail.std() > 0.02 * abs(tail.mean()) or abs(tail.mean()) > 2.0:
        raise SystemExit(f"the D3 add-on and the DFT program's D3 differ by more than a far tail: "
                         f"{tail.mean():.3f} +/- {tail.std():.3f} meV A3 per atom")
    args.data.mkdir(parents=True, exist_ok=True)
    write(args.data / "train.xyz", train, format="extxyz")
    write(args.data / "test.xyz", test, format="extxyz")
    summary = {
        "frames": len(frames), "train": len(train), "test": len(test),
        "test_ends": sum(a.info["family"] == "end" for a in test),
        "d3_addon_less_dft_program_far_tail_meV_A3_per_atom": [float(tail.mean()), float(tail.std())],
        "largest_force_eV_per_A": float(max(np.abs(a.arrays["REF_forces"]).max() for a in frames)),
    }
    (args.data / "prepare.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary))


def train(args):
    common = [
        "mace_run_train", f"--name={args.name}", f"--work_dir={args.data}",
        f"--train_file={args.data / 'train.xyz'}", "--valid_fraction=0.1", f"--test_file={args.data / 'test.xyz'}",
        "--energy_key=REF_energy", "--forces_key=REF_forces", "--stress_key=REF_stress",
        "--loss=stress", "--compute_stress=True", "--error_table=PerAtomRMSEstressvirials",
        "--ema", "--ema_decay=0.99", f"--batch_size={args.batch}", f"--max_num_epochs={args.epochs}",
        f"--default_dtype={args.dtype}", "--device=" + args.device, f"--seed={args.seed}", "--eval_interval=5",
    ]
    if args.start:
        # one stage on from a model of ours, with the weights given: its reference energies are kept
        common += [f"--energy_weight={args.energy_weight}", "--forces_weight=100", f"--stress_weight={args.stress_weight}",
                   f"--foundation_model={args.start}", "--multiheads_finetuning=False", "--E0s=foundation",
                   f"--lr={args.lr / 10}"]
    else:
        common += ["--energy_weight=1", "--forces_weight=100", "--stress_weight=10",
                   "--swa", f"--start_swa={int(0.6 * args.epochs)}", "--swa_energy_weight=1000",
                   "--swa_forces_weight=100", "--swa_stress_weight=10", f"--swa_lr={args.lr / 10}"]
    if args.start:
        pass
    elif args.scratch:
        common += ["--model=ScaleShiftMACE", "--num_channels=64", "--max_L=1", "--r_max=6.0", "--E0s=average",
                   f"--lr={10 * args.lr}"]
    else:
        common += [f"--foundation_model={args.foundation}", "--multiheads_finetuning=False", "--E0s=estimated",
                   f"--lr={args.lr}"]
    t0 = time.time()
    with open(args.data / f"{args.name}.train.log", "w") as log:
        code = subprocess.run(common, stdout=log, stderr=subprocess.STDOUT, check=False).returncode
    print(f"{args.name}: exit {code}, {(time.time() - t0) / 60:.1f} min; log {args.data / (args.name + '.train.log')}")
    sys.exit(code)


def held_out(calc, frames):
    """Errors on the held-out frames that are not relaxed cells, against the PBE part."""
    from ase import units

    de, df, ds = [], [], []
    for a in frames:
        if a.info["family"] == "end":
            continue
        b = a.copy()
        b.calc = calc
        de.append((b.get_potential_energy() - a.info["REF_energy"]) / len(a))
        df.append((b.get_forces() - a.arrays["REF_forces"]).ravel())
        ds.append((b.get_stress(voigt=False) - np.asarray(a.info["REF_stress"]).reshape(3, 3)).ravel() / units.GPa)
    de = np.array(de)
    return {
        "frames": len(de),
        # a constant offset per atom is the reference energies' business, not the fit's
        "energy_rmse_meV_per_atom": float(1e3 * np.sqrt(np.mean((de - de.mean()) ** 2))),
        "energy_offset_meV_per_atom": float(1e3 * de.mean()),
        "force_rmse_meV_per_A": float(1e3 * np.sqrt(np.mean(np.concatenate(df) ** 2))),
        "force_rms_of_the_labels_meV_per_A": float(1e3 * np.sqrt(np.mean(np.concatenate(
            [a.arrays["REF_forces"].ravel() for a in frames if a.info["family"] != "end"]) ** 2))),
        "stress_rmse_GPa": float(np.sqrt(np.mean(np.concatenate(ds) ** 2))),
    }


def elastic(calc, atoms, h, fmax):
    """Relaxed-ion C_11, C_22, C_33, C_12 of a cell: its own minimum, then stress under +/-h."""
    from ase import units
    from ase.filters import FrechetCellFilter
    from ase.optimize import BFGS

    a = atoms.copy()
    a.calc = calc
    BFGS(FrechetCellFilter(a), logfile=None).run(fmax=fmax, steps=400)
    c = np.zeros((3, 3))
    for j in range(3):
        sigma = []
        for sign in (+1, -1):
            b = a.copy()
            eps = np.eye(3)
            eps[j, j] += sign * h
            b.set_cell(a.cell.array @ eps, scale_atoms=True)
            b.calc = calc
            BFGS(b, logfile=None).run(fmax=fmax, steps=400)
            sigma.append(b.get_stress(voigt=True)[:3])
        c[:, j] = (sigma[0] - sigma[1]) / (2 * h) / units.GPa
    return {
        "C_11": float(c[0, 0]), "C_22": float(c[1, 1]), "C_33": float(c[2, 2]),
        "C_12": float((c[0, 1] + c[1, 0]) / 2), "C_12_asymmetry": float(abs(c[0, 1] - c[1, 0])),
        "lattice_A": [float(v) for v in a.cell.cellpar()[:3]],
    }


def turned_chain(atoms, frames_module):
    """The helix cell with its second chain turned half a turn about its own axis.

    In alpha-PVDF the two chains' dipoles across the chain cancel; turned, they add. It is a
    polar packing of the helix (delta-PVDF's kind) that no training frame was made from.
    """
    a = atoms.copy()
    chain = frames_module.pvdf_chains(a)[1]
    members = sorted({i for c1, g1, c2, g2 in chain for i in (c1, *g1, c2, *g2)})
    first = chain[0][0]
    rel = np.array([a.get_distance(first, i, mic=True, vector=True) for i in members])
    carbons = [members.index(i) for c1, g1, c2, g2 in chain for i in (c1, c2)]
    axis = a.cell.array[2] / np.linalg.norm(a.cell.array[2])
    rel -= rel[carbons].mean(axis=0)
    centre = a.positions[first] - rel[members.index(first)]
    a.positions[members] = centre + (2.0 * np.outer(rel @ axis, axis) - rel)
    return a


def evaluate(args):
    import torch
    from ase.calculators.mixing import SumCalculator
    from ase.io import read, write
    from mace.calculators import MACECalculator

    warnings.filterwarnings("ignore")
    torch.set_num_threads(args.threads)
    so = _second_opinion()
    bare = MACECalculator(model_paths=str(args.model), device="cpu", default_dtype="float64")
    out = {"model": args.model.name, "prepare": json.loads((args.data / "prepare.json").read_text())}
    out["held_out"] = held_out(bare, read(args.data / "test.xyz", index=":"))
    print("held out:", json.dumps(out["held_out"]), flush=True)
    kept = {}
    margins = so.run([(args.name, "file+d3", args.model)], True, 0.01, 400, args.threads, kept)
    out["margins"] = {"dft": margins["dft"], **margins["models"][args.name]}
    if args.probe:
        # the potential's own minima, and a packing it was not shown, for a DFT single point each
        from ase.filters import FrechetCellFilter
        from ase.optimize import BFGS

        mf = _example("mlip_frames")
        args.probe.mkdir(parents=True, exist_ok=True)
        source = {"alpha": "alpha", "beta": "polar", "anti": "anti"}
        turned = turned_chain(mf.cell("alpha"), mf)
        turned.calc = SumCalculator([bare, so.d3_calculator()])
        opt = BFGS(FrechetCellFilter(turned), logfile=None, maxstep=0.1)
        opt.run(fmax=0.01, steps=600)
        kept["turned"] = turned
        for k, a in kept.items():
            b = a.copy()
            b.calc = None
            b.info = {**mf.cell(source.get(k, "alpha")).info, "family": "probe", "frame": f"{args.name}_{k}",
                      "note": f"{k}, relaxed in {args.name} with D3" + (
                          ": the helix cell with its second chain turned half a turn" if k == "turned" else "")}
            write(args.probe / f"{args.name}_{k}.start.xyz", b, format="extxyz")
        out["probe_cells_written"] = sorted(kept)
        print(f"wrote {len(kept)} cells relaxed in the potential to {args.probe} for DFT single points "
              f"(turned-chain helix: {opt.get_number_of_steps()} steps)", flush=True)
    beta = so.dft_cells()["beta"][0]
    out["beta_elastic_GPa"] = {"ours": elastic(SumCalculator([bare, so.d3_calculator()]), beta, 0.005, 0.005),
                               "dft": DFT_ELASTIC_GPA}
    e = out["beta_elastic_GPa"]["ours"]
    print("beta elastic constants, GPa: " + ", ".join(
        f"{k} {e[k]:.1f} (DFT {v})" for k, v in DFT_ELASTIC_GPA.items()), flush=True)
    path = OUT / f"mlip_own_{args.name}.json"
    path.write_text(json.dumps(out, indent=1) + "\n")
    print(f"wrote {path.relative_to(ROOT)}")


def closure(args):
    """DFT single points on the cells the potential relaxed: is its minimum the DFT's?

    For each cell: the DFT energy at the potential's geometry against DFT's own minimum (never
    below it; how far above is what the potential's geometry costs), the DFT force and stress
    left there, and the gaps between cells by DFT and by the potential on the same geometries.
    """
    import torch
    from ase import units
    from ase.calculators.mixing import SumCalculator
    from ase.io import read
    from mace.calculators import MACECalculator

    warnings.filterwarnings("ignore")
    torch.set_num_threads(args.threads)
    so = _second_opinion()
    calc = SumCalculator([MACECalculator(model_paths=str(args.model), device="cpu", default_dtype="float64"),
                          so.d3_calculator()])
    own_minimum = {a.info["source"]: a.info["REF_energy"] / a.info["n_monomers"]
                   for a in read(args.frames, index=":") if a.info["family"] == "end"}
    source = {"alpha": "alpha", "beta": "polar", "anti": "anti"}
    cells = {}
    for f in sorted(args.probe.glob(f"{args.name}_*.dft.xyz")):
        a = read(f)
        k = a.info["frame"].replace(f"{args.name}_", "")
        b = a.copy()
        b.calc = calc
        n = a.info["n_monomers"]
        cells[k] = {
            "dft_eV_per_monomer": a.info["REF_energy"] / n,
            "potential_eV_per_monomer": float(b.get_potential_energy() / n),
            "dft_largest_force_eV_per_A": float(np.abs(a.arrays["REF_forces"]).max()),
            "dft_stress_xx_yy_zz_GPa": [float(v) for v in np.diag(np.asarray(a.info["REF_stress"]).reshape(3, 3)) / units.GPa],
        }
        if k in source:
            cells[k]["dft_above_its_own_minimum_meV_per_monomer"] = 1e3 * (cells[k]["dft_eV_per_monomer"] - own_minimum[source[k]])
    gaps = {}
    for k in cells:
        if k == "beta" or "beta" not in cells:
            continue
        gaps[f"{k}_minus_beta_meV"] = {
            who: 1e3 * (cells[k][f"{who}_eV_per_monomer"] - cells["beta"][f"{who}_eV_per_monomer"])
            for who in ("dft", "potential")}
        gaps[f"{k}_minus_beta_meV"]["difference_kcal"] = so.KCAL_PER_EV * 1e-3 * (
            gaps[f"{k}_minus_beta_meV"]["potential"] - gaps[f"{k}_minus_beta_meV"]["dft"])
    path = OUT / f"mlip_own_{args.name}.json"
    out = json.loads(path.read_text())
    out["dft_on_the_potentials_cells"] = {"cells": cells, "gaps": gaps}
    path.write_text(json.dumps(out, indent=1) + "\n")
    for k, c in cells.items():
        print(f"{k:7s} DFT force left {c['dft_largest_force_eV_per_A']:.3f} eV/A, stress "
              + " ".join(f"{v:+.2f}" for v in c["dft_stress_xx_yy_zz_GPa"]) + " GPa"
              + (f", DFT energy {c['dft_above_its_own_minimum_meV_per_monomer']:+.2f} meV per monomer above its own minimum"
                 if "dft_above_its_own_minimum_meV_per_monomer" in c else ""))
    for k, g in gaps.items():
        print(f"{k}: DFT {g['dft']:+.2f}, the potential {g['potential']:+.2f} ({g['difference_kcal']:+.3f} kcal/mol)")
    print(f"added to {path.relative_to(ROOT)}")


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = p.add_subparsers(dest="step", required=True)
    s = sub.add_parser("prepare")
    s.add_argument("--frames", type=Path, required=True, nargs="+", help="one or more files of labelled frames")
    s.add_argument("--data", type=Path, required=True, help="working directory, outside the repository")
    s = sub.add_parser("train")
    s.add_argument("--data", type=Path, required=True)
    s.add_argument("--name", required=True)
    s.add_argument("--foundation", default="medium-0b3", help="pretrained MACE-MP model to start from (MIT)")
    s.add_argument("--scratch", action="store_true", help="a small network from nothing instead")
    s.add_argument("--start", type=Path, help="a model of ours to fit on from, in one stage with the weights below")
    s.add_argument("--energy-weight", type=float, default=3e5,
                   help="with --start; forces are 100. The first fit's 1000 left energy at 0.4%% of its loss")
    s.add_argument("--stress-weight", type=float, default=1e4, help="with --start; the first fit's 10 left stress at 0.05%%")
    s.add_argument("--epochs", type=int, default=60)
    s.add_argument("--batch", type=int, default=4)
    s.add_argument("--lr", type=float, default=0.001,
                   help="fitting on from a pretrained model; 0.005 undid the pretrained weights in one epoch")
    s.add_argument("--dtype", default="float64", choices=["float32", "float64"])
    s.add_argument("--seed", type=int, default=3)
    s.add_argument("--device", default="cpu")
    s = sub.add_parser("evaluate")
    s.add_argument("--data", type=Path, required=True)
    s.add_argument("--model", type=Path, required=True)
    s.add_argument("--name", required=True)
    s.add_argument("--threads", type=int, default=4)
    s.add_argument("--probe", type=Path, help="write the cells relaxed in the potential here, for DFT single points")
    s = sub.add_parser("closure")
    s.add_argument("--frames", type=Path, required=True, help="the labelled frames (for DFT's own minima)")
    s.add_argument("--probe", type=Path, required=True, help="the directory `evaluate --probe` wrote, labelled")
    s.add_argument("--model", type=Path, required=True)
    s.add_argument("--name", required=True)
    s.add_argument("--threads", type=int, default=4)
    args = p.parse_args()
    {"prepare": prepare, "train": train, "evaluate": evaluate, "closure": closure}[args.step](args)


if __name__ == "__main__":
    main()
