"""The two forms of a mixed-cell evaluation on a row's own supercell: the loops against the batched form.

    python examples/mixed_forms_check.py --row vdf9-trfe4-cfe1-14 --row vdf10-trfe4-14 --nx 2 --ny 1

Every stored number was computed with the loop form (per chain and per chain pair, on the host).
The batched form is what runs on a device (``docs/PERFORMANCE_REVIEW.md`` 17.7); its tests are on
PVDF.  This evaluates the kinked crystal of a row of ``docs/predictions/kink_migration_belt4.json``
in the ``nx x ny`` tiling, at a random displacement of every variable, in both forms: the loops,
and the batched form on the array module the cell chooses by itself (the device where
``POLYFIND_EWALD_DEVICE=cuda`` and the cell is large enough; ``--host`` for NumPy whatever the
cell would choose).  Energy, gradient, closure and the closure's Jacobian must agree to rounding;
the exit code says whether they do.  A copolymer cell exercises what PVDF does not: sites that
differ along the chain, a 28-bond repeat, a chiral chain.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import sys
import time
import warnings
from pathlib import Path

import numpy as np

from polyfind import mixedcell as MC
from polyfind.fitting import PRESET_PARAMETERS

HERE = Path(__file__).resolve().parent
STORED = HERE.parent / "docs" / "predictions" / "kink_migration_belt4.json"


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def forms(row: str, nx: int, ny: int, kink_at: int = 0, table: dict | None = None, module=None, seed: int = 7, size: float = 0.4) -> dict:
    """Both forms at one displaced geometry of the row's kinked ``nx x ny`` supercell.

    ``module``: the batched form's array module (``None``: what the cell chooses).  Returns the
    energies, the largest differences, and the times of one evaluation each.  Call inside the
    preset the table was computed with."""
    DK = _load("dilute_kink")
    table = json.load(open(STORED)) if table is None else table
    letters, _ = DK.sequence(row)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # the isotactic-chain notice of a chiral repeat
        polymer = DK.C.make_polymer(row, letters)
        P = polymer.bonds_per_repeat
        seq = [int(v) for v in np.roll(np.array([1, 0, 2] + [0] * (P - 3)), kink_at % P)]
        d = DK.refined_cell(polymer, row, seq, f"{row} kink@{kink_at}", table["kinked_references"][f"{row}@{kink_at}"])
        model, x0, ch = DK.make_model(polymer, d, nx, ny)
    rng = np.random.default_rng(seed)
    x = x0.copy()
    x[model.index("internal")] += rng.uniform(-size, size, model.index("internal").size)
    x[model.index("place")] += rng.uniform(-0.05 * size, 0.05 * size, model.index("place").size)
    make = lambda xp: MC.MixedModel(polymer, "ref", model.pk, model.flips, [ch] * model.nc, angle_stiffness=DK.angle_stiffness(), xp=xp)  # noqa: E731
    loops, batched = make("loops"), make(module)
    chosen = batched.batched_module()
    if chosen is None:  # a cell too small for the device: the batched form on the host, so that there is something to compare
        batched = make(np)
        chosen = np
    batched.value_and_grad(x0)  # a device compiles its kernels on the first call
    t0 = time.perf_counter()
    e_l, g_l = loops.value_and_grad(x)
    t1 = time.perf_counter()
    e_b, g_b = batched.value_and_grad(x)
    t2 = time.perf_counter()
    (h_l, j_l), (h_b, j_b) = loops.closure(x), batched.closure(x)
    return {"row": row, "cell": [nx, ny], "chains": int(model.nc), "sites": int(model.pk.N), "variables": int(model.n_var),
            "module": getattr(chosen, "__name__", str(chosen)), "above_the_sheet": float(e_l - loops.value_and_grad(x0)[0]),
            "energy_loops": float(e_l), "energy_batched": float(e_b), "energy_difference": float(e_b - e_l),
            "gradient_difference": float(np.abs(g_b - g_l).max()), "gradient_largest": float(np.abs(g_l).max()),
            "closure_difference": float(np.abs(h_b - h_l).max()), "closure_jacobian_difference": float(np.abs(j_b - j_l).max()),
            "seconds_loops": t1 - t0, "seconds_batched": t2 - t1}


def agree(r: dict, rtol: float = 1e-10) -> bool:
    return (abs(r["energy_difference"]) <= rtol * max(1.0, abs(r["energy_loops"])) and r["gradient_difference"] <= rtol * max(1.0, r["gradient_largest"])
            and r["closure_difference"] <= 1e-9 and r["closure_jacobian_difference"] <= 1e-8)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--row", action="append", required=True, help="a row of the table's kinked references (repeatable)")
    ap.add_argument("--nx", type=int, default=2)
    ap.add_argument("--ny", type=int, default=1)
    ap.add_argument("--kink-at", type=int, default=0)
    ap.add_argument("--host", action="store_true", help="the batched form on NumPy, whatever the cell would choose")
    ap.add_argument("--table", default=str(STORED))
    ap.add_argument("--out", default=None, help="write the records here as JSON")
    args = ap.parse_args()
    table = json.load(open(args.table))
    out, ok = [], True
    with PRESET_PARAMETERS[table["preset"]].applied():
        for row in args.row:
            r = forms(row, args.nx, args.ny, args.kink_at, table, module=np if args.host else None)
            out.append(r)
            good = agree(r)
            ok = ok and good
            print(f"{row} {args.nx}x{args.ny}: {r['chains']} chains, {r['sites']} sites, {r['variables']} variables; batched on {r['module']}; "
                  f"{r['above_the_sheet']:+.1f} kcal/mol above the sheet; energy {r['energy_difference']:+.1e} of {r['energy_loops']:.3f}, "
                  f"gradient {r['gradient_difference']:.1e} of {r['gradient_largest']:.1f}, closure {r['closure_difference']:.1e}, "
                  f"its Jacobian {r['closure_jacobian_difference']:.1e}; loops {r['seconds_loops']:.2f} s, batched {r['seconds_batched']:.3f} s"
                  f"{'' if good else '   <-- THE FORMS DISAGREE'}", flush=True)
    if args.out:
        Path(args.out).write_text(json.dumps(out, indent=1))
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
