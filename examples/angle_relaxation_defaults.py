"""Which chemistries need ``fit_ris`` to relax the backbone angles: the measurement behind
``polyfind.forcefield.ANGLE_RELAXATION_DEFAULTS``.

    python examples/angle_relaxation_defaults.py              # the registered polymers
    python examples/angle_relaxation_defaults.py --candidates # plus the screen's four candidates
    python examples/angle_relaxation_defaults.py --preset illustrative --json out.json

For each chemistry :func:`polyfind.forcefield.angle_relaxation_check` takes the one-bond scan
``fit_ris`` makes (step 10 deg, six monomers, every other bond trans) twice per bond type --
at the polymer's frozen backbone angles and with every backbone angle of the oligomer relaxed
per point -- and lists the local minima of each.  A rigid minimum with no relaxed minimum
within 30 deg is an *orphan*: a well the frozen geometry manufactures, which ``fit_ris``'s
basin assignment would fold into whichever state lies nearest.  A chemistry with an orphan on
any bond type gets ``angles="relaxed"`` by default; one without keeps the rigid scan, whose
numbers every earlier table depends on.  docs/NITRILE_LANDSCAPE.md is where the artefact was
found (VDCN's +/-120 deg rigid well, 8.7 kcal/mol below trans, which the relaxed chain does
not have).

The copolymer entry ``vdf11-vdcn1`` is not measured here: its 24-bond repeat makes the check
a 54-bond oligomer with 57 angle variables per point, and nothing fits its RIS model directly
(``polyfind.ris.transfer_ris`` assembles it from the homopolymer fits).  ``fit_ris`` on it
would measure on the fly.  Wall time under the fitted potential, one core: 3-5 s per
homopolymer, 1-2 s more for a candidate; FANOME's overlapping all-trans reference makes its
relaxations stall at the 135 deg bound and its row is reported for completeness only.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
import warnings

import numpy as np

import polyfind.fitting  # noqa: F401  (registers the fitted presets)
from polyfind.forcefield import ANGLE_RELAXATION_DEFAULTS, SimpleFF, angle_relaxation_check
from polyfind.polymers import POLYMERS

warnings.simplefilter("ignore")


def fmt_minima(ms) -> str:
    return "; ".join(f"{a:+.0f}: {e:+.2f}" for a, e, _ in ms) or "none"


def check_all(polymers, calc, step=10.0, n_monomers=6, verbose=True) -> dict:
    out = {}
    for p in polymers:
        t0 = time.time()
        chk = angle_relaxation_check(p, calc, step=step, n_monomers=n_monomers)
        rec = {"needs_relaxation": bool(chk["needs_relaxation"]), "seconds": round(time.time() - t0, 1),
               "restraint_k": chk["restraint_k"],
               "relaxed_reference_angles": np.round(chk["relaxed_reference_angles"], 1).tolist(),
               "bonds": {str(b): {k: d[k] for k in ("rigid_minima", "relaxed_minima", "orphans")}
                         for b, d in chk["bonds"].items()}}
        recorded = ANGLE_RELAXATION_DEFAULTS.get(p.name)
        rec["recorded"] = recorded
        out[p.name] = rec
        if verbose:
            flag = "" if recorded is None else ("  (matches the recorded default)" if recorded == rec["needs_relaxation"]
                                                 else "  ** DISAGREES with the recorded default **")
            print(f"\n{p.name}: needs relaxation = {rec['needs_relaxation']}  ({rec['seconds']} s){flag}")
            for b, d in chk["bonds"].items():
                print(f"  bond {b}: rigid   {fmt_minima(d['rigid_minima'])}")
                print(f"          relaxed {fmt_minima(d['relaxed_minima'])}")
                if d["orphans"]:
                    print(f"          orphans {fmt_minima(d['orphans'])}")
            sys.stdout.flush()
    return out


def table(results: dict) -> str:
    rows = ["| polymer | rigid minima (deg: kcal/mol) | relaxed minima | orphans | default |", "|---|---|---|---|---|"]
    for name, r in results.items():
        b0 = r["bonds"]["0"]
        rows.append(f"| {name} | {fmt_minima(b0['rigid_minima'])} | {fmt_minima(b0['relaxed_minima'])} | "
                    f"{fmt_minima(b0['orphans'])} | {'relaxed' if r['needs_relaxation'] else 'rigid'} |")
    return "\n".join(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--preset", default="pvdf-dft-valence", help="SimpleFF preset (illustrative for the unfitted potential)")
    ap.add_argument("--candidates", action="store_true", help="also the screen's four candidate chemistries")
    ap.add_argument("--only", nargs="*", help="restrict to these polymer names")
    ap.add_argument("--json", help="write every number here")
    args = ap.parse_args()
    calc = SimpleFF() if args.preset == "illustrative" else SimpleFF.from_preset(args.preset)
    polymers = [p for n, p in POLYMERS.items() if n != "vdf11-vdcn1"]
    if args.candidates:
        from screen_electroactive import new_candidates
        polymers += new_candidates()
    if args.only:
        polymers = [p for p in polymers if p.name in args.only]
    results = check_all(polymers, calc)
    print("\n" + table(results))
    if args.json:
        with open(args.json, "w") as f:
            json.dump(results, f, indent=1)
        print(f"\nwrote {args.json}")


if __name__ == "__main__":
    sys.path.insert(0, __file__.rsplit("/", 1)[0] if "/" in __file__ else __file__.rsplit("\\", 1)[0])
    main()
