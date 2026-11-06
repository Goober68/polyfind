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
    # C/m^2.  Stages 1-4 used the published 0.18 (other functionals and cells, 0.178-0.181); from
    # stage 5 the anchor is Dipole's own Berry phase at the pvdf_beta_ref_v1 cell, 0.1927 (three
    # repeats to 3e-10; docs/predictions/e1_dipole_existing.json, record beta_qe_berry_v1).
    "p_beta": 0.1927,
    "p_beta_published": 0.18,
    # relaxed-ion normal stiffness of beta, GPa, derived by polyfind from Dipole's six +/-1% axial QE
    # cells of 2026-11-05 (docs/predictions/e1_dipole_existing.json, record beta_qe_axial_v1):
    # stress and energy routes agree within 1.5%
    "beta_C": {"C11": 22.77, "C22": 24.52, "C12": 5.88},
}
TOL = {"cell": 0.01, "p": 0.01, "margin": 0.05, "gap": 0.05, "stiffness": 0.10}  # relative for cell and stiffness, absolute otherwise

VARIABLES = ("x_C", "x_H", "x_F", "D_C", "D_H", "D_F", "q_CH", "q_CF")
# multipliers on the valence values for the Lennard-Jones terms, absolute shifts for the increments
X0 = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0])
BOUNDS = np.array([(0.90, 1.10)] * 3 + [(0.5, 2.0)] * 3 + [(-0.08, 0.08)] * 2)
SIGMA = np.array([0.05, 0.05, 0.05, 0.35, 0.35, 0.35, 0.05, 0.05])  # ridge widths (log scale for D)
RIDGE = 1.0
# Stage 2 (docs/CRYSTAL_STAGE.md section 5): the same eight plus a pair-specific H-F term,
# multipliers on the geometric means of the mixed H and F diameter and well depth; fitted in the
# induced-dipole, force-shifted Lennard-Jones setting the E1/E2 predictions are made in, with
# the beta stiffness as anchors.
VARIABLES2 = VARIABLES + ("x_HF", "D_HF")
X0_2 = np.concatenate([X0, [1.0, 1.0]])
BOUNDS2 = np.vstack([BOUNDS, [(0.85, 1.15), (0.3, 3.0)]])
SIGMA2 = np.concatenate([SIGMA, [0.05, 0.5]])
STAGE2 = {"polarizable": True, "lj_cutoff": "force", "stiffness": True}
# Stage 3 (docs/CRYSTAL_STAGE.md section 5.2): stage 2 plus charge penetration, the permanent
# charges interacting as 1/r - (1 + b r/2) exp(-b r)/r in the lattice sum (EwaldSpec.charge_damping,
# b in 1/A).  b = 4 is point charges to 3e-4 at 2.5 A; the fit may bring it down to 1.
VARIABLES3 = VARIABLES2 + ("b_q",)
X0_3 = np.concatenate([X0_2, [4.0]])
BOUNDS3 = np.vstack([BOUNDS2, [(1.0, 6.0)]])
SIGMA3 = np.concatenate([SIGMA2, [2.0]])  # a weak pull towards point charges
# Stage 4 (docs/CRYSTAL_STAGE.md section 5.7): stage 2's ten plus the depth of the fluorine belt
# sites' wall against H (log10 of kcal/mol; -2 is the 0.01 that put beta's polar axis on the DFT
# value in the probe), with the charge penetration carried at its point-charge value.  The chains
# are polymers.pvdf_belt(BELT_D): four massless repulsion-only sites per CF2.
BELT_D = 0.7  # A, distance of the sites from the fluorine, perpendicular to the C-F bond
VARIABLES4 = VARIABLES2 + ("log10_D_X", "b_q")
X0_4 = np.concatenate([X0_2, [-2.0, 4.0]])
BOUNDS4 = np.vstack([BOUNDS2, [(-3.5, -0.3), (4.0, 4.0 + 1e-9)]])
SIGMA4 = np.concatenate([SIGMA2, [1.0, 2.0]])


def parameters(x, base: F.FFParameters = F.FITTED_VALENCE) -> F.FFParameters:
    """The valence preset with its C, H, F Lennard-Jones terms and C-H, C-F increments moved
    (:func:`polyfind.fitting.crystal_stage_parameters`; the fitted result is the preset
    ``pvdf-dft-valence-crystal``)."""
    return F.crystal_stage_parameters(x, base)


def _config(x):
    """``(variables, x0, bounds, sigma)`` for an 8- or 10-variable vector."""
    n = len(np.asarray(x))
    if n == 8:
        return VARIABLES, X0, BOUNDS, SIGMA
    if n == 10:
        return VARIABLES2, X0_2, BOUNDS2, SIGMA2
    if n == 11:
        return VARIABLES3, X0_3, BOUNDS3, SIGMA3
    if n == 12:
        return VARIABLES4, X0_4, BOUNDS4, SIGMA4
    raise ValueError(f"expected 8, 10, 11 or 12 variables, got {n}")


def ridge(x) -> float:
    x = np.asarray(x, dtype=float)
    _, _, _, sigma = _config(x)
    parts = [(x[:3] - 1.0) / sigma[:3], np.log(x[3:6]) / sigma[3:6], x[6:8] / sigma[6:8]]
    if len(x) >= 10:  # the pair term: a diameter multiplier and a well-depth multiplier
        parts += [(x[8:9] - 1.0) / sigma[8:9], np.log(x[9:10]) / sigma[9:10]]
    if len(x) == 11:  # charge penetration, measured from the point-charge start
        parts += [(x[10:11] - X0_3[10]) / sigma[10:11]]
    if len(x) == 12:  # the belt wall (log depth, from the probe's 0.01) and the carried penetration
        parts += [(x[10:11] - X0_4[10]) / sigma[10:11], (x[11:12] - X0_4[11]) / sigma[11:12]]
    d = np.concatenate(parts)
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


def reference_from_prediction(pred: F.Prediction, case, **packer_kw):
    """The :class:`polyfind.mechanics.Reference` of a refined prediction, built the way
    :func:`polyfind.mechanics.refined_reference` builds one (chain rebuilt from the refined torsions
    and angles in the refinement's own frame), for the response functions of :mod:`polyfind.mechanics`.
    Call inside ``params.applied()`` so the packer is the fitted one."""
    from .linegroup import repeat_chains
    from .mechanics import reference_from_chain
    rr, rigid, polymer = pred.refine_result, pred.rigid, case.polymer
    angles0 = np.array([polymer.backbone[k].backbone_angle for k in range(polymer.bonds_per_repeat)])
    frame = repeat_chains(polymer, rigid.chain, np.asarray(rigid.dihedrals, dtype=float)[None], angles0)[0]
    built = repeat_chains(polymer, rigid.chain, rr.torsions[None], rr.angles, align_to=frame)[0]
    return reference_from_chain(built, rr.result.params, n_chains=2, **packer_kw)


def cases_for(x):
    """The beta and alpha cases the variables ``x`` are evaluated on: the plain PVDF cases, or for
    stage 4 (twelve variables) the same sequences on :func:`polyfind.polymers.pvdf_belt`, under
    their own keys so the cached chains and screen tables stay separate."""
    from dataclasses import replace as _replace
    beta_case, alpha_case = F.FITTED_CASES[1], F.FITTED_CASES[2]
    if len(np.asarray(x)) == 12:
        from .polymers import pvdf_belt
        pol = pvdf_belt(BELT_D)
        return (_replace(beta_case, key=f"beta-belt-{BELT_D:.2f}", polymer=pol),
                _replace(alpha_case, key=f"alpha-belt-{BELT_D:.2f}", polymer=pol))
    return beta_case, alpha_case


def evaluate(x, base: F.FFParameters = F.FITTED_VALENCE, coulomb: str = "ewald", verbose: bool = False,
             polarizable: bool = False, lj_cutoff: str = "energy", stiffness: bool = False) -> Evaluation:
    """Pack and refine beta and alpha under ``x``; return the normalised residuals and the objective.

    ``polarizable`` adds the induced dipoles (:class:`polyfind.polarizability.Polarizable`),
    ``lj_cutoff="force"`` the force-shifted Lennard-Jones term, and ``stiffness`` the relaxed-ion
    C_11, C_22 and C_12 of beta against Dipole's strained-cell values: together the construction
    of the E1/E2 predictions (:data:`STAGE2`).
    """
    from .mechanics import elastic_constants, relax_reference
    from .polarizability import Polarizable
    t0 = time.time()
    params = parameters(x, base)
    ewald = EwaldSpec() if coulomb == "ewald" else None
    pol = Polarizable() if polarizable else None
    beta_case, alpha_case = cases_for(x)
    try:
        beta = F.predict(beta_case, params, coulomb=coulomb, ewald=ewald, always_gap=True,
                         polarizable=pol, lj_cutoff=lj_cutoff)
        alpha = F.predict(alpha_case, params, coulomb=coulomb, ewald=ewald, polarizable=pol, lj_cutoff=lj_cutoff)
        C = None
        if stiffness:
            with params.applied():
                ref = reference_from_prediction(beta, beta_case, coulomb=coulomb, ewald=ewald,
                                                polarizable=pol, lj_cutoff=lj_cutoff)
                ref = relax_reference(ref)
                C = elastic_constants(ref).block  # Voigt 1, 2, 6 in-plane block, GPa
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
    if C is not None:
        for name, pred, key in (("beta C_11 (polar)", float(C[0, 0]), "C11"), ("beta C_22 (lateral)", float(C[1, 1]), "C22"),
                                ("beta C_12", float(C[0, 1]), "C12")):
            target = DFT_ANCHORS["beta_C"][key]
            res[name] = (pred, target, (pred - target) / (TOL["stiffness"] * target))
    rg = ridge(x)
    obj = float(sum(r[2] ** 2 for r in res.values())) + rg
    return Evaluation(np.asarray(x, dtype=float), res, obj, rg, time.time() - t0, beta, alpha)


def table(ev: Evaluation) -> str:
    lines = [f"{'quantity':32s} {'polyfind':>10s} {'DFT':>10s} {'z':>7s}"]
    for name, (pred, target, z) in ev.residuals.items():
        lines.append(f"{name:32s} {pred:10.4f} {target:10.4f} {z:7.2f}")
    lines.append(f"{'ridge':32s} {'':>10s} {'':>10s} {ev.ridge:7.2f}   objective {ev.objective:.2f}   ({ev.seconds:.0f} s)")
    return "\n".join(lines)


def fit(maxfev: int = 240, coulomb: str = "ewald", log=None, x0=None, stage: int = 1,
        **eval_kw) -> tuple[Evaluation, list[Evaluation]]:
    """Powell with bounds on the eight (stage 1) or ten (stage 2) variables; every evaluation is
    logged.  ``stage=2`` uses :data:`X0_2`, :data:`BOUNDS2` and the :data:`STAGE2` evaluation
    settings unless ``eval_kw`` overrides them; ``x0`` resumes from an earlier best."""
    from scipy.optimize import minimize
    if stage in (2, 3, 4):
        eval_kw = {**STAGE2, **eval_kw}
    start = ({1: X0, 2: X0_2, 3: X0_3, 4: X0_4}[stage]) if x0 is None else np.asarray(x0, dtype=float)
    _, _, bounds, _ = _config(start)
    history: list[Evaluation] = []

    def obj(x):
        xc = np.clip(x, bounds[:, 0], bounds[:, 1])
        ev = evaluate(xc, coulomb=coulomb, **eval_kw)
        history.append(ev)
        if log is not None:
            print(f"[{len(history):3d}] obj {ev.objective:9.2f} ridge {ev.ridge:6.2f}  x = {np.array2string(xc, precision=4, suppress_small=True)}  ({ev.seconds:.0f} s)",
                  file=log, flush=True)
        return ev.objective

    first = evaluate(start, coulomb=coulomb, **eval_kw)
    history.append(first)
    if log is not None:
        print("start:\n" + table(first), file=log, flush=True)
    minimize(obj, start, method="Powell", bounds=[tuple(b) for b in bounds],
             options={"maxfev": maxfev, "xtol": 1e-3, "ftol": 1e-3})
    best = min(history, key=lambda e: e.objective)
    return best, history


def save(ev: Evaluation, path: str) -> None:
    variables = _config(ev.x)[0]
    prm = parameters(ev.x)
    out = {"variables": list(variables), "x": ev.x.tolist(), "objective": ev.objective, "ridge": ev.ridge,
           "residuals": {k: {"polyfind": v[0], "dft": v[1], "z": v[2]} for k, v in ev.residuals.items()},
           "parameters": {"lj": prm.lj, "lj_pairs": prm.lj_pairs, "charge_increments": prm.charge_increments,
                          "charge_damping": prm.charge_damping},
           "anchors": DFT_ANCHORS, "tolerances": TOL}
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1)
