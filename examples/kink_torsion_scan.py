"""The energy of a lone kink pair against its own torsion: has the kinked state a rim?

examples/dilute_kink.py relaxes one kinked chain among all-trans neighbours in stages.  Whether
the result of a stage is a state with a lifetime or a shelf on the way back to all-trans cannot
be read off a relaxation: a quasi-Newton iteration is not bounded in its later steps, and a
damped descent with short steps can crawl for hours along a flat stretch without stopping or
arriving.  This script asks the question directly.  From the stored state of a stage of a stored
row, the two gauche bonds of the kink pair are set to ``+tau`` and ``-tau`` and held, and the
other variables of the stage are relaxed.  The energy against ``tau`` then shows the rim, if
there is one, as a rise between the kinked minimum (near 100 degrees for a strained kink) and
the eclipsed position at 120 degrees, and its height.

    python examples/kink_torsion_scan.py --row vdf10 --nx 2 --ny 1 --stage host_rigid_bodies --tau 97,99,101,102,103,104,106,112,120,130
    python examples/kink_torsion_scan.py --row vdf10 --nx 4 --ny 2 --stage host_rigid_bodies --tau 100,112,120,130

Every point starts from the same stored state, so a smooth curve (and a smoothly varying rest of
the chain, printed) is the sign that the points lie on one branch; a jump in either is not a
barrier but another branch.  Results go into the row as ``torsion_scan[stage][tau]``, kcal/mol
above the all-trans supercell, and are skipped if already there.
"""
import argparse
import importlib.util
import json
import time
from pathlib import Path

import numpy as np

from polyfind.fitting import PRESET_PARAMETERS

DS = importlib.util.spec_from_file_location("dk", Path(__file__).with_name("dilute_kink.py"))
DK = importlib.util.module_from_spec(DS)
DS.loader.exec_module(DK)


def held_kink(model, x_stage, free, closed, bonds, tau: float, maxiter: int = 600):
    """The stage's variables relaxed with the torsions ``bonds`` of chain 0 held at ``+tau, -tau``
    (signs as stored).  Returns ``(x, info)``."""
    x = np.asarray(x_stage, dtype=float).copy()
    t = model.index("tors", 0)
    held = np.array([t[bonds[0]], t[bonds[1]]])
    sign = np.sign(DK.wrap(x[held]))
    x[held] = sign * tau
    rest = np.setdiff1d(model.index("internal", 0), held)
    x = model.project(x, rest, [0])
    free2 = np.setdiff1d(free, held)
    xr, info = model.relax(x, free=free2, closed=closed, maxiter=maxiter)
    for _ in range(2):
        if "success" in info["message"].lower():
            break
        xr, more = model.relax(xr, free=free2, closed=closed, maxiter=maxiter)
        info = {**more, "iterations": info["iterations"] + more["iterations"], "evaluations": info["evaluations"] + more["evaluations"]}
    return xr, info


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--row", default="vdf10")
    ap.add_argument("--nx", type=int, default=2)
    ap.add_argument("--ny", type=int, default=1)
    ap.add_argument("--kink-at", type=int, default=0)
    ap.add_argument("--stage", default="host_rigid_bodies", help="the stage whose variables are relaxed and whose stored state every point starts from")
    ap.add_argument("--tau", default="97,99,101,102,103,104,106,112,120,130", help="comma list of the held gauche angle, degrees")
    ap.add_argument("--preset", default="pvdf-dft-valence-belt4")
    ap.add_argument("--maxiter", type=int, default=600)
    ap.add_argument("--out", default="docs/predictions/dilute_kink_belt4.json")
    args = ap.parse_args()
    log = lambda s: print(s, flush=True)  # noqa: E731
    out_path = Path(args.out)
    key = f"{args.row}@{args.kink_at} {args.nx}x{args.ny} kinked 0"
    disk = json.load(open(out_path))
    row = disk["rows"][key]
    if args.stage not in row.get("states", {}):
        raise SystemExit(f"no stored state of stage {args.stage} for {key}")
    letters, _ = DK.sequence(args.row)
    with PRESET_PARAMETERS[args.preset].applied():
        polymer = DK.C.make_polymer(args.row, letters)
        d = DK.all_trans_reference(polymer, args.row, disk["all_trans_references"][args.row])
        model, x_tt, _ = DK.make_model(polymer, d, args.nx, args.ny)
        e_tt = model.physical_energy(x_tt)
        x_stage = np.asarray(row["states"][args.stage], dtype=float)
        free, closed = DK.stage_sets(model, [0])[args.stage]
        P = model.P
        bonds = (args.kink_at % P, (args.kink_at + 2) % P)
        start = DK.wrap(x_stage[model.index("tors", 0)[list(bonds)]])
        log(f"{key}, stage {args.stage}: stored state {model.physical_energy(x_stage) - e_tt:+.4f} kcal/mol above all-trans, gauche bonds {np.round(start, 1).tolist()}")
        done = row.get("torsion_scan", {}).get(args.stage, {})
        for tau in [float(v) for v in args.tau.split(",") if v]:
            if f"{tau:g}" in done:
                log(f"   +/-{tau:g}: {done[f'{tau:g}']['energy']:+.4f} (stored)")
                continue
            t0 = time.time()
            xr, info = held_kink(model, x_stage, free, closed, bonds, tau, args.maxiter)
            tors = DK.wrap(model.state(xr).tors[0])
            others = float(np.abs(np.abs(np.delete(tors, list(bonds))) - 180.0).max())
            rec = {"energy": float(model.physical_energy(xr) - e_tt), "largest_other_torsion_off_trans_deg": others, "iterations": int(info["iterations"]),
                   "evaluations": int(info["evaluations"]), "message": info["message"], "seconds": time.time() - t0, "state": xr.tolist()}
            log(f"   +/-{tau:g}: {rec['energy']:+.4f}; largest other torsion off trans {others:.1f} deg; {info['iterations']} it, {info['evaluations']} ev, {time.time() - t0:.0f} s; {info['message']}")
            import fcntl
            with open(str(out_path) + ".lock", "w") as lock:
                fcntl.flock(lock, fcntl.LOCK_EX)
                fresh = json.load(open(out_path))
                fresh["rows"][key].setdefault("torsion_scan", {}).setdefault(args.stage, {})[f"{tau:g}"] = rec
                with open(out_path, "w") as fh:
                    json.dump(fresh, fh, indent=1)
    log(f"written {out_path}")


if __name__ == "__main__":
    main()
