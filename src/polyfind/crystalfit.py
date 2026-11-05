"""The crystal stage of the fit: nonbonded terms against the DFT crystal anchors.

``docs/VALENCE_FIT.md`` fitted torsions, bonds, angles, Lennard-Jones multipliers and
bond-charge increments to 467 frames of PBE-D3 conformers.  That data cannot see a crystal:
relative conformer energies barely see the overall dispersion (so the Lennard-Jones terms
were bounded rather than determined), and energies and forces cannot separate a charge
scale from a dielectric constant at all (only a polarization can).  The result, measured
against Dipole's periodic PBE-D3(BJ) cells on 2026-11-05, is a potential that reproduces the
alpha-beta margin to 0.02 kcal/mol per monomer and the lateral and chain axes of beta to
0.1% and 0.7%, but packs beta's polar axis 4.5% too close, alpha's long axis 6% too close,
and carries a fixed-charge polarization of 0.115 C/m^2 against the DFT 0.18.

This stage moves only the terms the trainset left undetermined -- the Lennard-Jones
multipliers of C, H and F and the C-H and C-F bond-charge increments -- against the DFT
anchors below, with a ridge that keeps them near the trainset values, and leaves every
torsion, bond and angle term where the valence fit put it.  Every evaluation packs and
refines beta and the antipolar alpha cell exactly as the crystal fit of ``docs/SCREEN.md``
does (:func:`polyfind.fitting.predict`), in the fixed-charge Ewald setting the margins
were validated in.  The induced-dipole setting is evaluated afterwards, not fitted.
"""
from __future__ import annotations

import json
import time
from dataclasses import dataclass

import numpy as np

from . import fitting as F
from .ewald import EwaldSpec

# --- the anchors: Dipole's periodic PBE-D3(BJ), polyfind axis order --------------------
DFT_ANCHORS = {
    "beta_cell": (4.757, 8.528, 2.581),  # QE 7.6, pvdf_beta_ref_v1 (2026-11-04): polar, lateral, chain
    "alpha_cell": (4.945, 9.636, 4.660),  # QE, polyfind's antipolar TGTG' packing (docs/POLAR_MARGIN_DFT.md)
    "alpha_minus_beta": -0.726,  # kcal/mol per monomer, two-body D3
    "beta_antipolar_minus_polar": -0.026,  # kcal/mol per monomer; the antipolar beta packing is lower
    "p_beta": 0.18,  # C/m^2, published DFT Berry phase for beta (0.178-0.181)
}
TOL = {"cell": 0.01, "p": 0.01, "margin": 0.05, "gap": 0.05}  # relative for the cell, absolute otherwise

VARIABLES = ("x_C", "x_H", "x_F", "D_C", "D_H", "D_F", "q_CH", "q_CF")
# multipliers on the valence values for the Lennard-Jones terms, absolute shifts for the increments
X0 = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0])
BOUNDS = np.array([(0.90, 1.10)] * 3 + [(0.5, 2.0)] * 3 + [(-0.08, 0.08)] * 2)
SIGMA = np.array([0.05, 0.05, 0.05, 0.35, 0.35, 0.35, 0.05, 0.05])  # ridge widths (log scale for D)
RIDGE = 1.0


def parameters(x, base: F.FFParameters = F.FITTED_VALENCE) -> F.FFParameters:
    """The valence preset with its C, H, F Lennard-Jones terms and C-H, C-F increments moved
    (:func:`polyfind.fitting.crystal_stage_parameters`; the fitted result is the preset
    ``pvdf-dft-valence-crystal``)."""
    return F.crystal_stage_parameters(x, base)


def ridge(x) -> float:
    x = np.asarray(x, dtype=float)
    d = np.concatenate([(x[:3] - 1.0) / SIGMA[:3], np.log(x[3:6]) / SIGMA[3:6], x[6:] / SIGMA[6:]])
    return RIDGE * float(np.sum(d * d))


@dataclass
class Evaluation:
    x: np.ndarray
    residuals: dict
    objective: float
    ridge: float
    seconds: float
    beta: object = None
    alpha: object = None


def evaluate(x, base: F.FFParameters = F.FITTED_VALENCE, coulomb: str = "ewald", verbose: bool = False) -> Evaluation:
    """Pack and refine beta and alpha under ``x``; return the normalised residuals and the objective."""
    t0 = time.time()
    params = parameters(x, base)
    ewald = EwaldSpec() if coulomb == "ewald" else None
    beta_case, alpha_case = F.FITTED_CASES[1], F.FITTED_CASES[2]
    try:
        beta = F.predict(beta_case, params, coulomb=coulomb, ewald=ewald, always_gap=True)
        alpha = F.predict(alpha_case, params, coulomb=coulomb, ewald=ewald)
    except Exception as exc:  # a parameter set that cannot pack is simply bad
        if verbose:
            print(f"  evaluate: {type(exc).__name__}: {exc}", flush=True)
        return Evaluation(np.asarray(x, dtype=float), {}, F.BIG, ridge(x), time.time() - t0)
    res = {}
    ba, bb = beta.axes
    aa, ab = alpha.axes
    for name, pred, target in (("beta a (polar)", ba, DFT_ANCHORS["beta_cell"][0]),
                               ("beta b (lateral)", bb, DFT_ANCHORS["beta_cell"][1]),
                               ("beta c (chain)", beta.c, DFT_ANCHORS["beta_cell"][2]),
                               ("alpha a", aa, DFT_ANCHORS["alpha_cell"][0]),
                               ("alpha b", ab, DFT_ANCHORS["alpha_cell"][1]),
                               ("alpha c", alpha.c, DFT_ANCHORS["alpha_cell"][2])):
        res[name] = (pred, target, (pred - target) / (TOL["cell"] * target))
    res["beta |P|"] = (beta.polarization, DFT_ANCHORS["p_beta"], (beta.polarization - DFT_ANCHORS["p_beta"]) / TOL["p"])
    gap = alpha.energy_per_monomer - beta.energy_per_monomer
    res["E(alpha)-E(beta)"] = (gap, DFT_ANCHORS["alpha_minus_beta"], (gap - DFT_ANCHORS["alpha_minus_beta"]) / TOL["margin"])
    if beta.polar_gap is not None:
        res["beta antipolar-polar (rigid)"] = (beta.polar_gap, DFT_ANCHORS["beta_antipolar_minus_polar"],
                                             (beta.polar_gap - DFT_ANCHORS["beta_antipolar_minus_polar"]) / TOL["gap"])
    rg = ridge(x)
    obj = float(sum(r[2] ** 2 for r in res.values())) + rg
    return Evaluation(np.asarray(x, dtype=float), res, obj, rg, time.time() - t0, beta, alpha)


def table(ev: Evaluation) -> str:
    lines = [f"{'quantity':32s} {'polyfind':>10s} {'DFT':>10s} {'z':>7s}"]
    for name, (pred, target, z) in ev.residuals.items():
        lines.append(f"{name:32s} {pred:10.4f} {target:10.4f} {z:7.2f}")
    lines.append(f"{'ridge':32s} {'':>10s} {'':>10s} {ev.ridge:7.2f}   objective {ev.objective:.2f}   ({ev.seconds:.0f} s)")
    return "\n".join(lines)


def fit(maxfev: int = 240, coulomb: str = "ewald", log=None, x0=None) -> tuple[Evaluation, list[Evaluation]]:
    """Powell with bounds on the eight variables; every evaluation is logged."""
    from scipy.optimize import minimize
    history: list[Evaluation] = []

    def obj(x):
        xc = np.clip(x, BOUNDS[:, 0], BOUNDS[:, 1])
        ev = evaluate(xc, coulomb=coulomb)
        history.append(ev)
        if log is not None:
            print(f"[{len(history):3d}] obj {ev.objective:9.2f} ridge {ev.ridge:6.2f}  x = {np.array2string(xc, precision=4, suppress_small=True)}  ({ev.seconds:.0f} s)",
                  file=log, flush=True)
        return ev.objective

    start = X0 if x0 is None else np.asarray(x0, dtype=float)
    first = evaluate(start, coulomb=coulomb)
    history.append(first)
    if log is not None:
        print("start:\n" + table(first), file=log, flush=True)
    minimize(obj, start, method="Powell", bounds=[tuple(b) for b in BOUNDS],
             options={"maxfev": maxfev, "xtol": 1e-3, "ftol": 1e-3})
    best = min(history, key=lambda e: e.objective)
    return best, history


def save(ev: Evaluation, path: str) -> None:
    out = {"variables": list(VARIABLES), "x": ev.x.tolist(), "objective": ev.objective, "ridge": ev.ridge,
           "residuals": {k: {"polyfind": v[0], "dft": v[1], "z": v[2]} for k, v in ev.residuals.items()},
           "parameters": {"lj": parameters(ev.x).lj, "charge_increments": parameters(ev.x).charge_increments},
           "anchors": DFT_ANCHORS, "tolerances": TOL}
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1)
