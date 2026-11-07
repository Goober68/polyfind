"""Sensitivity of the crystal-stage residuals at a fitted point, and what the form can still move.

    python examples/crystal_stage_sensitivity.py --fit docs/predictions/crystal_stage_v7.json --out docs/predictions/crystal_stage_v7_sensitivity.json
    python examples/crystal_stage_sensitivity.py --analyze docs/predictions/crystal_stage_v7_sensitivity.json --verify

The first form takes central differences of every residual (in tolerance units, the fit's z) with
respect to every free variable at the fit's point, two evaluations per variable (a variable whose
bounds coincide, the carried charge penetration ``b_q``, is skipped), on the fit's own stage
(inferred from its residual rows) and belt geometry, and writes the Jacobian beside the fit.  The
second form reads it back and does the linear analysis of ``docs/CRYSTAL_STAGE.md`` 5.15: the
sensitivity per prior width (the ridge's sigma, so an entry of 1 means one width moves the
residual by one tolerance), the singular values of that matrix (how many directions of the
variable space the anchors fix to within a width), the minimum-norm step that moves each residual
alone, and the best values of the ``--target`` residuals reachable on the linear model with every
other residual held within ``--hold`` tolerances and the variables inside their bounds.
``--verify`` evaluates that point for real, which is the check the linear model needs (the
alpha-beta gap is the most nonlinear row), and records it in the file.
"""
import argparse
import json
import time

import numpy as np

from polyfind import crystalfit as C

STEPS = {"x_C": 0.005, "x_H": 0.005, "x_F": 0.005, "D_C": 0.02, "D_H": 0.02, "D_F": 0.02, "q_CH": 0.005,
         "q_CF": 0.005, "x_HF": 0.005, "D_HF": 0.02, "log10_D_X": 0.05, "b_q": 0.1, "log10_k_angle": 0.02}


def stage_kwargs(residual_rows) -> dict:
    """The evaluate() flags of the stage a saved fit was made on, read off its residual rows."""
    kw = dict(C.STAGE2) if "beta C_11 (polar)" in residual_rows else {}
    if "beta d|P|/d eta (polar)" in residual_rows:
        kw["slope"] = True
    if "beta backbone angle" in residual_rows:
        kw["angles"] = True
    return kw


def free_variables(x, variables, bounds):
    return [i for i, n in enumerate(variables) if bounds[i][1] - bounds[i][0] > 1e-6 and n in STEPS]


def jacobian(fit_path: str, out_path: str, geometry: str | None) -> dict:
    with open(fit_path) as fh:
        fit = json.load(fh)
    x0 = np.asarray(fit["x"], dtype=float)
    variables, _, bounds, _ = C._config(x0)
    kw = stage_kwargs(fit["residuals"])
    C.BELT_GEOMETRY = geometry or fit.get("belt", {}).get("geometry", C.BELT_GEOMETRY)
    t0 = time.time()
    e0 = C.evaluate(x0, **kw)
    rows = list(e0.residuals)
    out = {"fit": fit_path, "x0": x0.tolist(), "variables": list(variables), "rows": rows,
           "z0": [e0.residuals[k][2] for k in rows], "pred0": [e0.residuals[k][0] for k in rows],
           "target": [e0.residuals[k][1] for k in rows], "objective0": float(e0.objective), "ridge0": float(e0.ridge),
           "stage_kw": kw, "geometry": C.BELT_GEOMETRY, "steps": {}, "dz_dx": {}, "dpred_dx": {}}
    print(f"base objective {e0.objective:.3f} ({time.time() - t0:.0f} s); {len(rows)} rows", flush=True)
    for i in free_variables(x0, variables, bounds):
        n, h = variables[i], STEPS[variables[i]]
        xp, xm = x0.copy(), x0.copy()
        xp[i] += h
        xm[i] -= h
        ep, em = C.evaluate(xp, **kw), C.evaluate(xm, **kw)
        out["steps"][n] = h
        out["dz_dx"][n] = [(ep.residuals[k][2] - em.residuals[k][2]) / (2 * h) for k in rows]
        out["dpred_dx"][n] = [(ep.residuals[k][0] - em.residuals[k][0]) / (2 * h) for k in rows]
        with open(out_path, "w") as fh:
            json.dump(out, fh, indent=1)
        print(f"{n:14s} +-{h}: objective {ep.objective:.3f} / {em.objective:.3f}   ({time.time() - t0:.0f} s)", flush=True)
    return out


def prior_widths(x0, variables, sigma):
    """The ridge's widths in the units of x: the depth multipliers' widths are log widths."""
    return np.array([sigma[i] * (x0[i] if variables[i].startswith("D_") else 1.0) for i in range(len(variables))])


def analyze(path: str, targets, hold: float, verify: bool) -> None:
    from scipy.optimize import minimize
    with open(path) as fh:
        d = json.load(fh)
    x0 = np.asarray(d["x0"], dtype=float)
    variables, _, bounds, sigma = C._config(x0)
    names = [n for n in variables if n in d["dz_dx"]]
    idx = [list(variables).index(n) for n in names]
    rows, z0 = d["rows"], np.asarray(d["z0"], dtype=float)
    J = np.array([d["dz_dx"][n] for n in names]).T  # rows x free variables
    s = prior_widths(x0, variables, sigma)[idx]
    S = J * s[None, :]
    print(f"point: objective {d['objective0']:.3f}, ridge {d['ridge0']:.3f}; {len(rows)} residuals, {len(names)} free variables")
    print("\nsensitivity dz per prior width (|1| = one width moves the residual one tolerance):")
    print(f"{'':32s}" + "".join(f"{n:>9s}" for n in names) + "      z0")
    for i, r in enumerate(rows):
        print(f"{r:32s}" + "".join(f"{S[i, j]:9.2f}" for j in range(len(names))) + f"   {z0[i]:+.2f}")
    sv = np.linalg.svd(S, compute_uv=False)
    print("\nsingular values:", np.round(sv, 2).tolist(), "  (> 1:", int((sv > 1).sum()), "directions the anchors fix within a width)")
    Jp = np.linalg.pinv(S)
    print("\nminimum-norm step (prior widths) moving one residual by one tolerance with the others held:")
    for i, r in enumerate(rows):
        e = np.zeros(len(rows))
        e[i] = 1.0
        dx = Jp @ e
        k = int(np.argmax(np.abs(dx)))
        print(f"  {r:32s} |dx| {np.linalg.norm(dx):8.2f}   largest {names[k]:>13s} {dx[k]:+8.2f}   unreachable part {np.linalg.norm(S @ dx - e):.2f}")
    tgt = [rows.index(t) for t in targets]
    others = [i for i in range(len(rows)) if i not in tgt]
    lob = np.array([bounds[i][0] for i in idx]) - x0[idx]
    hib = np.array([bounds[i][1] for i in idx]) - x0[idx]

    def zlin(dx):
        return z0 + J @ dx
    cons = ([{"type": "ineq", "fun": lambda dx, i=i: hold - zlin(dx)[i]} for i in others]
            + [{"type": "ineq", "fun": lambda dx, i=i: hold + zlin(dx)[i]} for i in others])
    res = minimize(lambda dx: float((zlin(dx)[tgt] ** 2).sum()), np.zeros(len(names)), method="SLSQP",
                   bounds=list(zip(lob, hib)), constraints=cons, options={"maxiter": 500, "ftol": 1e-12})
    zl = zlin(res.x)
    print(f"\nlinear model: best {targets} with every other residual within {hold} and the variables in bounds (success {res.success}):")
    for t in tgt:
        print(f"  {rows[t]:32s} z {z0[t]:+.2f} -> {zl[t]:+.2f}")
    print("  step in prior widths:", {n: round(float(v), 2) for n, v in zip(names, res.x / s) if abs(v) > 0.05 * s[names.index(n)]})
    print("  rows at the limit:", [rows[i] for i in others if abs(zl[i]) > 0.98 * hold])
    record = {"targets": targets, "hold": hold, "x": None, "z_linear": zl.tolist(), "step_widths": (res.x / s).tolist()}
    x1 = x0.copy()
    x1[idx] += res.x
    record["x"] = x1.tolist()
    if verify:
        C.BELT_GEOMETRY = d.get("geometry", C.BELT_GEOMETRY)
        kw = d.get("stage_kw") or stage_kwargs(rows)
        t0 = time.time()
        ev = C.evaluate(x1, **kw)
        print(f"\nevaluated for real: objective {ev.objective:.2f}, ridge {ev.ridge:.2f} ({time.time() - t0:.0f} s)")
        print(f"  {'row':32s} {'z0':>7s} {'linear':>7s} {'true':>7s}")
        for i, r in enumerate(rows):
            print(f"  {r:32s} {z0[i]:+7.2f} {zl[i]:+7.2f} {ev.residuals[r][2]:+7.2f}")
        record.update({"z_true": [ev.residuals[r][2] for r in rows], "objective_true": float(ev.objective), "ridge_true": float(ev.ridge)})
    d["trade_off"] = record
    with open(path, "w") as fh:
        json.dump(d, fh, indent=1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fit", help="a crystal-stage result file to differentiate at")
    ap.add_argument("--out", help="where the Jacobian goes (with --fit)")
    ap.add_argument("--geometry", default=None, choices=(None, "textbook", "pbe"), help="belt chain geometry (default: the fit's)")
    ap.add_argument("--analyze", help="a Jacobian file to analyse")
    ap.add_argument("--target", action="append", default=None, help="residual rows to minimise on the linear model (default: beta C_12 and beta b)")
    ap.add_argument("--hold", type=float, default=1.0, help="tolerances every other residual is held within")
    ap.add_argument("--verify", action="store_true", help="evaluate the linear model's point for real")
    args = ap.parse_args()
    if args.fit:
        jacobian(args.fit, args.out or args.fit.replace(".json", "_sensitivity.json"), args.geometry)
    if args.analyze:
        analyze(args.analyze, args.target or ["beta C_12", "beta b (lateral)"], args.hold, args.verify)
    if not (args.fit or args.analyze):
        ap.error("give --fit (to differentiate) and/or --analyze (to analyse)")


if __name__ == "__main__":
    main()
