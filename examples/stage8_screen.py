"""Screen candidate stage-8 terms against a fitted point's Jacobian: does a new term add a direction?

    python examples/stage8_screen.py --fit docs/predictions/crystal_stage_v7_trf2.json --sensitivity docs/predictions/crystal_stage_v7_trf2_sensitivity.json --out docs/predictions/stage8_screen_belt3.json

At the fitted point every candidate is perturbed both ways and the change of the residual vector
(tolerance units) per step is projected onto the column space of the Jacobian of the twelve fitted
variables (``examples/crystal_stage_sensitivity.py``).  The part outside that span is what the new
term could do that the fitted variables cannot already do by moving; what it would take out of the residual the
twelve variables leave is reported both for an unbounded move of the candidate (and how many
steps that needs) and within one step, which is the physically meaningful number.  The
candidates are the two fixed belt constants (the sites' distance from the fluorine and their
Lennard-Jones diameter) and pair-specific Lennard-Jones terms beside the fitted F-H one (F-F, H-H,
C-C, C-F, C-H, each at the mixing rule's value times a factor), and two electrostatic directions the form
carries as fixed quantities: the fluorine charge's offset from its nucleus and the polarizability
scale.  docs/CRYSTAL_STAGE.md 5.20.
"""
import argparse
import json
import time
from dataclasses import replace

import numpy as np

from polyfind import crystalfit as C
from polyfind import fitting as F


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fit", default="docs/predictions/crystal_stage_v7_trf2.json")
    ap.add_argument("--sensitivity", default="docs/predictions/crystal_stage_v7_trf2_sensitivity.json")
    ap.add_argument("--out", default="docs/predictions/stage8_screen_belt3.json")
    args = ap.parse_args()
    fit = json.load(open(args.fit))
    sens = json.load(open(args.sensitivity))
    x0 = np.asarray(fit["x"], dtype=float)
    C.BELT_GEOMETRY = fit.get("belt", {}).get("geometry", C.BELT_GEOMETRY)
    kw = sens.get("stage_kw") or dict(C.STAGE7)
    rows = sens["rows"]
    names = [n for n in sens["variables"] if n in sens["dz_dx"]]
    J = np.array([sens["dz_dx"][n] for n in names]).T
    z0 = np.asarray(sens["z0"], dtype=float)
    # the linear least-squares residual the twelve variables leave (bounds ignored): what a new column could still reduce
    coef0, *_ = np.linalg.lstsq(J, -z0, rcond=None)
    r_min = z0 + J @ coef0
    print(f"point: objective {fit['objective']:.3f}; anchors' sum of squares {float(z0 @ z0):.2f}, the linear model's attainable minimum with the twelve variables {float(r_min @ r_min):.2f}")
    base = C.parameters(x0)
    lj = {e: (xi, di) for e, xi, di in base.lj}
    orig = C.parameters

    def evaluate():
        e = C.evaluate(x0, **kw)
        return np.array([e.residuals[k][2] for k in rows]), float(e.objective)

    def pair_patch(a, b, xmul, dmul):
        xa, da = lj[a]
        xb, db = lj[b]
        pair = (a, b, float(np.sqrt(xa * xb) * xmul), float(np.sqrt(da * db) * dmul))

        def patched(x, base=F.FITTED_VALENCE):
            p = orig(x, base)
            return replace(p, lj_pairs=tuple(p.lj_pairs) + (pair,))
        return patched

    probes = []
    probes.append(("belt distance d, +-0.1 A", lambda s: setattr(C, "BELT_D", 0.7 + 0.1 * s), lambda: setattr(C, "BELT_D", 0.7)))
    probes.append(("belt site diameter, +-0.2 A", lambda s: setattr(F, "BELT_X_DIAMETER", 2.6 + 0.2 * s), lambda: setattr(F, "BELT_X_DIAMETER", 2.6)))
    for label, a, b, xm, dm in (("F-F pair diameter x1.02", "F", "F", 1.02, 1.0), ("F-F pair depth x1.5", "F", "F", 1.0, 1.5),
                                ("H-H pair diameter x1.02", "H", "H", 1.02, 1.0), ("C-C pair diameter x1.02", "C", "C", 1.02, 1.0),
                                ("C-F pair diameter x1.02", "C", "F", 1.02, 1.0), ("C-H pair diameter x1.02", "C", "H", 1.02, 1.0)):
        probes.append((label, (lambda s, a=a, b=b, xm=xm, dm=dm: setattr(C, "parameters", pair_patch(a, b, xm ** s, dm ** s))),
                       lambda: setattr(C, "parameters", orig)))
    # two electrostatic shape and response directions the form already has as fixed quantities
    import polyfind.polarizability as PZ
    orig_pol = PZ.Polarizable
    offs = dict(base.charge_offsets or ())
    f_off = float(offs.get("F", 0.0))

    def offset_patch(delta):
        def patched(x, b=F.FITTED_VALENCE):
            p = orig(x, b)
            return replace(p, charge_offsets=tuple((e, (d + delta) if e == "F" else d) for e, d in (p.charge_offsets or (("F", f_off),))))
        return patched

    def pol_scale(s):
        return lambda *a, **k: orig_pol(*a, **{**k, "scale": float(k.get("scale", 1.0)) * s})

    probes.append((f"F charge offset {f_off:.3f} +- 0.03 A", lambda s: setattr(C, "parameters", offset_patch(0.03 * s)), lambda: setattr(C, "parameters", orig)))
    probes.append(("polarizability scale x 1.05", lambda s: setattr(PZ, "Polarizable", pol_scale(1.05 ** s)), lambda: setattr(PZ, "Polarizable", orig_pol)))
    out = {"fit": args.fit, "rows": rows, "attainable_minimum_12": float(r_min @ r_min), "probes": {}}
    for label, apply, restore in probes:
        t0 = time.time()
        apply(+1)
        zp, op = evaluate()
        apply(-1)
        zm, om = evaluate()
        restore()
        dz = (zp - zm) / 2.0
        coef, *_ = np.linalg.lstsq(J, dz, rcond=None)
        resid = dz - J @ coef  # the new direction, per step of the candidate
        rr = float(resid @ resid)
        # moving the candidate by t steps (the twelve variables re-optimised) leaves |r_min - t resid|^2:
        # the full gain needs t* = (r_min . resid) / |resid|^2 steps; within one step the gain is bounded
        t_star = float(r_min @ resid) / rr if rr > 1e-12 else 0.0
        gain_full = float((r_min @ resid) ** 2 / rr) if rr > 1e-12 else 0.0
        t1 = float(np.clip(t_star, -1.0, 1.0))
        gain_one_step = float(2 * t1 * (r_min @ resid) - t1 * t1 * rr)
        out["probes"][label] = {"dz_per_step": dz.tolist(), "objective_minus": om, "objective_plus": op,
                                "norm": float(np.linalg.norm(dz)), "outside_span": float(np.sqrt(rr)),
                                "steps_for_full_gain": t_star, "gain_full": gain_full, "gain_within_one_step": gain_one_step}
        big = [(rows[i], dz[i]) for i in range(len(rows)) if abs(dz[i]) > 0.15]
        print(f"{label:30s} |dz| {np.linalg.norm(dz):5.2f} per step, outside the span {np.sqrt(rr):5.2f}; gain within one step "
              f"{gain_one_step:5.3f} of {float(r_min @ r_min):.2f} (full gain {gain_full:.2f} would need {t_star:+.1f} steps); "
              f"objective -/+ {om:.2f}/{op:.2f}  ({time.time() - t0:.0f} s)")
        print("    " + "; ".join(f"{r} {v:+.2f}" for r, v in big))
        with open(args.out, "w") as fh:
            json.dump(out, fh, indent=1)
    print(f"written {args.out}")


if __name__ == "__main__":
    main()
