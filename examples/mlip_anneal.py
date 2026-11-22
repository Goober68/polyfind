"""Cells shaken at a temperature in the potential and then relaxed (docs/MLIP_SECOND_OPINION.md 3.24): does a cell that was
relaxed where it was built come to rest lower once its chains have been let move?

    python examples/mlip_anneal.py --model MODEL --name all_s1 --cells STARTS.extxyz [--pick REGEX] --out DIR \\
        --kelvin 350 --md-steps 1500 [--fs 1.0] [--seed 1] [--fmax 0.004] [--steps 700]

Each picked frame is run at --kelvin for --md-steps steps of --fs femtoseconds in its own cell (Langevin, the cell held),
then relaxed, atoms and lattice, as examples/mlip_train.py relax does.  Written to DIR/NAME_in_POTENTIAL.start.xyz, in the
form examples/registry_cells.py table reads; the frame keeps the start's registry_* entries.
"""

import argparse
import importlib.util
import re
import warnings
from pathlib import Path

import numpy as np


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--model", type=Path, required=True)
    p.add_argument("--name", required=True)
    p.add_argument("--cells", type=Path, required=True)
    p.add_argument("--pick", default="")
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--kelvin", type=float, default=350.0)
    p.add_argument("--md-steps", type=int, default=1500)
    p.add_argument("--fs", type=float, default=1.0)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--fmax", type=float, default=0.004)
    p.add_argument("--steps", type=int, default=700)
    p.add_argument("--threads", type=int, default=1)
    args = p.parse_args()

    import torch
    from ase import units
    from ase.calculators.mixing import SumCalculator
    from ase.filters import FrechetCellFilter
    from ase.io import read, write
    from ase.md.langevin import Langevin
    from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary
    from ase.optimize import BFGS

    warnings.filterwarnings("ignore")
    torch.set_num_threads(args.threads)
    spec = importlib.util.spec_from_file_location("mlip_train", Path(__file__).resolve().parent / "mlip_train.py")
    train = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(train)
    so = train._second_opinion()
    calc = SumCalculator([so.file_calculator(args.model), so.d3_calculator()])
    args.out.mkdir(parents=True, exist_ok=True)
    frames = [a for a in read(args.cells, index=":") if not args.pick or re.fullmatch(args.pick, str(a.info["frame"]))]
    print(f"{len(frames)} frames", flush=True)
    for start in frames:
        a = start.copy()
        a.calc = calc
        n = int(start.info["n_monomers"])
        e0 = a.get_potential_energy() / n
        rng = np.random.default_rng(args.seed)
        MaxwellBoltzmannDistribution(a, temperature_K=args.kelvin, rng=rng)
        Stationary(a)
        dyn = Langevin(a, timestep=args.fs * units.fs, temperature_K=args.kelvin, friction=0.01 / units.fs, rng=rng, logfile=None)
        hottest = 0.0
        for _ in range(max(1, args.md_steps // 100)):
            dyn.run(min(100, args.md_steps))
            hottest = max(hottest, a.get_temperature())
        e_hot = a.get_potential_energy() / n
        a.set_momenta(np.zeros((len(a), 3)))
        opt = BFGS(FrechetCellFilter(a), logfile=None, maxstep=0.1)
        done = bool(opt.run(fmax=args.fmax, steps=args.steps))
        e1 = a.get_potential_energy() / n
        frame = f"{start.info['frame']}_in_{args.name}"
        b = a.copy()
        b.calc = None
        b.info = {**{k: v for k, v in start.info.items() if k not in ("REF_energy", "REF_stress", "dft_energy_ry", "dft_d3_ry", "dft_kpoints", "dft_seconds", "dft_conv_thr_ry", "dft_run")},
                  "frame": frame, "family": "relaxed",
                  "note": f"{start.info.get('source', start.info['frame'])}, run {args.md_steps} steps of {args.fs} fs at {args.kelvin:.0f} K in {args.name} with D3 (seed {args.seed}), "
                          f"then relaxed to {args.fmax} eV/A ({opt.get_number_of_steps()} steps{'' if done else ', NOT converged'})",
                  "potential_eV_per_monomer": float(e1), "potential_start_eV_per_monomer": float(e0), "annealed_kelvin": args.kelvin, "annealed_md_steps": args.md_steps}
        for k in list(b.arrays):
            if k not in ("numbers", "positions"):
                del b.arrays[k]
        write(args.out / f"{frame}.start.xyz", b, format="extxyz")
        print(f"{frame}: {args.md_steps} steps at {args.kelvin:.0f} K (hottest {hottest:.0f} K, the potential's energy {1e3 * (e_hot - e0):+.1f} meV per monomer there), then "
              f"{opt.get_number_of_steps()} steps, converged {done}; the potential's energy {1e3 * (e1 - e0):+.2f} meV per monomer from the start", flush=True)


if __name__ == "__main__":
    main()
