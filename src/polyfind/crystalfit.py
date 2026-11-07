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
    # d|P|/d eta along the polar axis, C/m^2: the central difference of Dipole's Berry-phase
    # polarization at the -1%, 0, +1% polar-axis cells (ions relaxed at each fixed cell; 0.19600,
    # 0.19269, 0.18949), formed by polyfind (record beta_qe_berry_axial_v1).  For the longitudinal
    # component the proper and improper piezoelectric coefficients coincide, so on one branch this
    # is the relaxed-ion e_11; Dipole asserts no coefficient (its branch replay is running).
    "beta_dP_deta_polar": -0.3255,
    # Backbone C-C-C angles, degrees (docs/CRYSTAL_STAGE.md 5.12).  Beta: both carbons of Dipole's
    # native pvdf_beta_ref_v1, a planar zigzag (dihedral 180.00) with C-C 1.5331 A, so that
    # c = 2 r sin(theta/2) = 2.5807 is the cell's chain axis exactly.  Alpha: the CH2 and CF2
    # carbons of the PBE-D3(BJ) start geometry of pvdf_alpha_ref_v1, relaxed at the anchor cell
    # 4.945/9.636/4.660 (torsions 179.1, 59.9); Dipole's native alpha relaxation is still running.
    "beta_backbone_angle": 114.63,
    "alpha_backbone_angles": (118.23, 116.87),
}
TOL = {"cell": 0.01, "p": 0.01, "margin": 0.05, "gap": 0.05, "stiffness": 0.10, "slope": 0.03,
       "angle": 1.0}  # relative for cell and stiffness, absolute otherwise (angle in degrees)

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
BELT_GEOMETRY = "textbook"  # or "pbe": the chain built at Dipole's PBE-D3 bond lengths and angles (polymers.pvdf_belt)
VARIABLES4 = VARIABLES2 + ("log10_D_X", "b_q")
X0_4 = np.concatenate([X0_2, [-2.0, 4.0]])
BOUNDS4 = np.vstack([BOUNDS2, [(-3.5, -0.3), (4.0, 4.0 + 1e-9)]])
SIGMA4 = np.concatenate([SIGMA2, [1.0, 2.0]])
# Stage 6 (docs/CRYSTAL_STAGE.md section 5.11): the belt form of stage 4 under both of Dipole's
# Berry-phase anchors, P at the reference cell (0.1927, the anchor from stage 5 on) and d|P|/d eta
# along the polar axis from the -1/0/+1% cells (-0.3255 C/m^2); the strained states the stiffness
# block relaxes give the slope for nothing.
STAGE6 = {**STAGE2, "slope": True}
# Stage 7 (docs/CRYSTAL_STAGE.md section 5.13): stage 6 plus the stiffness of the refinement's
# harmonic backbone-angle restraint as a variable (log10 of its multiplier on
# refine.DEFAULT_ANGLE_STIFFNESS, 105 kcal/mol/rad^2, an invented constant until here), with the
# backbone angles of Dipole's beta and alpha cells as anchors.  Meant for the PBE-D3 chain
# geometry (BELT_GEOMETRY "pbe"), which builds the chain at beta's DFT angle so the restraint's
# equilibrium is the DFT value.
VARIABLES7 = VARIABLES4 + ("log10_k_angle",)
X0_7 = np.concatenate([X0_4, [0.0]])
BOUNDS7 = np.vstack([BOUNDS4, [(-0.7, 0.7)]])  # a fifth to five times the default
SIGMA7 = np.concatenate([SIGMA4, [0.3]])  # a factor of two is one width
STAGE7 = {**STAGE6, "angles": True}
STRAIN_STEP = 2e-3  # the strain of the stiffness and slope states (mechanics.strain_states' default)


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
    if n == 13:
        return VARIABLES7, X0_7, BOUNDS7, SIGMA7
    raise ValueError(f"expected 8, 10, 11, 12 or 13 variables, got {n}")


def ridge_parts(x) -> np.ndarray:
    """The ridge's components ``d_j`` (distance from the valence value in prior widths), one per
    variable, so that :func:`ridge` is ``RIDGE * sum(d^2)`` and the least-squares fitter can carry
    ``sqrt(RIDGE) * d`` as extra residual rows."""
    x = np.asarray(x, dtype=float)
    _, _, _, sigma = _config(x)
    parts = [(x[:3] - 1.0) / sigma[:3], np.log(x[3:6]) / sigma[3:6], x[6:8] / sigma[6:8]]
    if len(x) >= 10:  # the pair term: a diameter multiplier and a well-depth multiplier
        parts += [(x[8:9] - 1.0) / sigma[8:9], np.log(x[9:10]) / sigma[9:10]]
    if len(x) == 11:  # charge penetration, measured from the point-charge start
        parts += [(x[10:11] - X0_3[10]) / sigma[10:11]]
    if len(x) >= 12:  # the belt wall (log depth, from the probe's 0.01) and the carried penetration
        parts += [(x[10:11] - X0_4[10]) / sigma[10:11], (x[11:12] - X0_4[11]) / sigma[11:12]]
    if len(x) == 13:  # the backbone-angle restraint, log10 of its multiplier on the default
        parts += [x[12:13] / sigma[12:13]]
    return np.concatenate(parts)


def ridge(x) -> float:
    d = ridge_parts(x)
    return RIDGE * float(np.sum(d * d))


# Forward-difference steps of the least-squares fitter's Jacobian, per variable (the units of x).
FD_STEPS = {"x_C": 0.005, "x_H": 0.005, "x_F": 0.005, "D_C": 0.02, "D_H": 0.02, "D_F": 0.02, "q_CH": 0.005,
            "q_CF": 0.005, "x_HF": 0.005, "D_HF": 0.02, "log10_D_X": 0.05, "b_q": 0.1, "log10_k_angle": 0.02}


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
    if len(np.asarray(x)) >= 12:
        from .polymers import pvdf_belt
        pol = pvdf_belt(BELT_D, BELT_GEOMETRY)
        tag = f"belt-{BELT_D:.2f}" + ("" if BELT_GEOMETRY == "textbook" else f"-{BELT_GEOMETRY}")
        return (_replace(beta_case, key=f"beta-{tag}", polymer=pol),
                _replace(alpha_case, key=f"alpha-{tag}", polymer=pol))
    return beta_case, alpha_case


def polar_slope(states: dict, step: float = STRAIN_STEP) -> float:
    """``d|P_polar|/d eta`` along the polar axis from the ``(0, +1)`` and ``(0, -1)`` states of
    :func:`polyfind.mechanics.strain_states`: the central difference of the true polarization's
    polar component with the cell's orientation (the sign of P) taken out, so that it compares
    with the finite difference of a Berry-phase polarization whatever either side's sign
    convention.  Negative when stretching the polar axis lowers the polarization."""
    px_plus = float(states[(0, 1)].polarization[0])
    px_minus = float(states[(0, -1)].polarization[0])
    sgn = float(np.sign(px_plus + px_minus)) or 1.0
    return sgn * (px_plus - px_minus) / (2.0 * step)


def evaluate(x, base: F.FFParameters = F.FITTED_VALENCE, coulomb: str = "ewald", verbose: bool = False,
             polarizable: bool = False, lj_cutoff: str = "energy", stiffness: bool = False,
             slope: bool = False, angles: bool = False) -> Evaluation:
    """Pack and refine beta and alpha under ``x``; return the normalised residuals and the objective.

    ``polarizable`` adds the induced dipoles (:class:`polyfind.polarizability.Polarizable`),
    ``lj_cutoff="force"`` the force-shifted Lennard-Jones term, and ``stiffness`` the relaxed-ion
    C_11, C_22 and C_12 of beta against Dipole's strained-cell values: together the construction
    of the E1/E2 predictions (:data:`STAGE2`).  ``slope`` adds ``d|P|/d eta`` along the polar axis
    (:func:`polar_slope`) against Dipole's three-point Berry pair (:data:`STAGE6`), from the same
    strained states as the stiffness.  ``angles`` adds the refined backbone angles of both
    polymorphs against Dipole's (:data:`STAGE7`), which the refinement already has.
    """
    from .mechanics import elastic_constants, relax_reference, strain_states
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
        dp_slope = None
        if stiffness or slope:
            with params.applied():
                ref = reference_from_prediction(beta, beta_case, coulomb=coulomb, ewald=ewald,
                                                polarizable=pol, lj_cutoff=lj_cutoff)
                ref = relax_reference(ref)
                states = strain_states(ref, STRAIN_STEP)
                if stiffness:
                    C = elastic_constants(ref, STRAIN_STEP, states=states).block  # Voigt 1, 2, 6 in-plane block, GPa
                if slope:
                    dp_slope = polar_slope(states, STRAIN_STEP)
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
    if dp_slope is not None:
        target = DFT_ANCHORS["beta_dP_deta_polar"]
        res["beta d|P|/d eta (polar)"] = (dp_slope, target, (dp_slope - target) / TOL["slope"])
    if angles:
        b_ang = np.asarray(beta.refine_result.angles, dtype=float)  # [CH2 carbon, CF2 carbon]
        a_ang = np.asarray(alpha.refine_result.angles, dtype=float)
        t = DFT_ANCHORS["beta_backbone_angle"]
        res["beta backbone angle"] = (float(b_ang.mean()), t, (float(b_ang.mean()) - t) / TOL["angle"])
        for name, pred, t in (("alpha angle at CH2", float(a_ang[0]), DFT_ANCHORS["alpha_backbone_angles"][0]),
                              ("alpha angle at CF2", float(a_ang[1]), DFT_ANCHORS["alpha_backbone_angles"][1])):
            res[name] = (pred, t, (pred - t) / TOL["angle"])
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
        checkpoint: str | None = None, method: str = "powell", **eval_kw) -> tuple[Evaluation, list[Evaluation]]:
    """Powell with bounds on the eight (stage 1) or ten (stage 2) variables; every evaluation is
    logged.  ``stage=2`` uses :data:`X0_2`, :data:`BOUNDS2` and the :data:`STAGE2` evaluation
    settings unless ``eval_kw`` overrides them; ``x0`` resumes from an earlier best.  With
    ``checkpoint``, the best evaluation so far is :func:`save`\ d to that path (``complete``
    false) each time it improves, so a run cut short keeps its best point and resumes from it.

    ``method="trf"`` fits the same objective as a bounded least-squares problem instead
    (:func:`scipy.optimize.least_squares`, trust-region reflective): the residual vector is the
    ``z`` of every anchor row plus ``sqrt(RIDGE)`` times each ridge component (:func:`ridge_parts`),
    so its sum of squares is the objective, and the Jacobian is a forward difference of that vector
    with the per-variable steps :data:`FD_STEPS` (one evaluation per free variable, which can run in
    parallel where Powell's line searches cannot).  ``maxfev`` then bounds the residual evaluations
    least_squares asks for itself; the Jacobian's evaluations come on top and are in ``history``."""
    from scipy.optimize import least_squares, minimize
    if stage in (2, 3, 4):
        eval_kw = {**STAGE2, **eval_kw}
    elif stage == 6:
        eval_kw = {**STAGE6, **eval_kw}
    elif stage == 7:
        eval_kw = {**STAGE7, **eval_kw}
    start = ({1: X0, 2: X0_2, 3: X0_3, 4: X0_4, 6: X0_4, 7: X0_7}[stage]) if x0 is None else np.asarray(x0, dtype=float)
    if stage == 7 and len(start) == 12:  # from a stage-4/5/6 result: the restraint starts at the default
        start = np.concatenate([start, [0.0]])
    _, _, bounds, _ = _config(start)
    history: list[Evaluation] = []

    def record(ev):
        history.append(ev)
        if checkpoint is not None and ev.objective <= min(e.objective for e in history):
            save(ev, checkpoint, evaluations=len(history), complete=False)

    def obj(x):
        xc = np.clip(x, bounds[:, 0], bounds[:, 1])
        ev = evaluate(xc, coulomb=coulomb, **eval_kw)
        record(ev)
        if log is not None:
            print(f"[{len(history):3d}] obj {ev.objective:9.2f} ridge {ev.ridge:6.2f}  x = {np.array2string(xc, precision=4, suppress_small=True)}  ({ev.seconds:.0f} s)",
                  file=log, flush=True)
        return ev.objective

    first = evaluate(start, coulomb=coulomb, **eval_kw)
    record(first)
    if log is not None:
        print("start:\n" + table(first), file=log, flush=True)
    if method == "powell":
        minimize(obj, start, method="Powell", bounds=[tuple(b) for b in bounds],
                 options={"maxfev": maxfev, "xtol": 1e-3, "ftol": 1e-3})
    elif method == "trf":
        variables = _config(start)[0]
        rows = list(first.residuals)
        n_ridge = len(ridge_parts(start))
        m = len(rows) + n_ridge
        free = [i for i in range(len(start)) if bounds[i, 1] - bounds[i, 0] > 1e-6]

        def residual_vector(ev):
            if not ev.residuals:  # a parameter set that could not pack
                return np.full(m, np.sqrt(F.BIG / m))
            return np.concatenate([[ev.residuals[k][2] for k in rows], np.sqrt(RIDGE) * ridge_parts(ev.x)])

        def fun(x):
            xc = np.clip(np.asarray(x, dtype=float), bounds[:, 0], bounds[:, 1])
            ev = evaluate(xc, coulomb=coulomb, **eval_kw)
            record(ev)
            if log is not None:
                print(f"[{len(history):3d}] obj {ev.objective:9.2f} ridge {ev.ridge:6.2f}  x = {np.array2string(xc, precision=4, suppress_small=True)}  ({ev.seconds:.0f} s)",
                      file=log, flush=True)
            return residual_vector(ev)

        def jac(x):
            x = np.asarray(x, dtype=float)
            r0 = residual_vector(history[-1]) if np.allclose(history[-1].x, x, atol=1e-12) else fun(x)
            J = np.zeros((m, len(x)))
            for i in free:
                h = FD_STEPS.get(variables[i], 0.005)
                if x[i] + h > bounds[i, 1]:
                    h = -h
                xp = x.copy()
                xp[i] += h
                J[:, i] = (fun(xp) - r0) / h
            return J

        lo, hi = bounds[:, 0].copy(), bounds[:, 1].copy()
        x_start = np.clip(start, lo, hi)
        least_squares(fun, x_start, jac=jac, bounds=(lo, hi), method="trf", x_scale="jac",
                      ftol=1e-3, xtol=1e-3, gtol=1e-4, max_nfev=maxfev)
    else:
        raise ValueError(f"unknown method {method!r} (powell or trf)")
    best = min(history, key=lambda e: e.objective)
    return best, history


def save(ev: Evaluation, path: str, evaluations: int | None = None, complete: bool = True) -> None:
    """Write an evaluation as a result file (the format ``--x0`` resumes from).  ``evaluations``
    is how many the run had made when this was its best; ``complete`` false marks a checkpoint
    of a run still going (or cut short)."""
    variables = _config(ev.x)[0]
    prm = parameters(ev.x)
    out = {"variables": list(variables), "x": ev.x.tolist(), "objective": ev.objective, "ridge": ev.ridge,
           "complete": bool(complete), "evaluations": evaluations,
           "residuals": {k: {"polyfind": v[0], "dft": v[1], "z": v[2]} for k, v in ev.residuals.items()},
           "parameters": {"lj": prm.lj, "lj_pairs": prm.lj_pairs, "charge_increments": prm.charge_increments,
                          "charge_damping": prm.charge_damping},
           "anchors": DFT_ANCHORS, "tolerances": TOL}
    if len(ev.x) >= 12:
        out["belt"] = {"d": BELT_D, "geometry": BELT_GEOMETRY}
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1)
