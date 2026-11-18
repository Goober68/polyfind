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
# Stage 11 (docs/CRYSTAL_STAGE.md 5.24; stages 8 and 9 are that document's two designs that were
# never built, 5.21 and 5.23, and 10 is left free beside them): stage 7 plus what it costs to move
# neighbouring chains along each other.  Plain PVDF, two chains of one monomer, sheared rigidly
# along the chains in two packings (polyfind's beta cell, and the all-trans packing of the
# reversed-unit candidate with plain chains in it), against PBE-D3(BJ) single points on the same
# atoms (docs/predictions/shear_registry_belt4.json; examples/shear_registry.py).  Each packing is
# its mixed-cell state (lattice A, B, gamma, c; per chain fa, fb, setting angle, height; the
# backbone angles; which chain is flipped), fixed whatever the potential, so the atoms are the ones
# DFT was given.  Rows: the shear of a and of b in monomers and DFT's energy against the unsheared
# cell, kcal/mol per monomer.  The preset before this stage is 2.4 to 3.5 times DFT on them.
REGISTRY = {
    "beta": {"lattice": (4.754070390732875, 8.480849480989322, 90.00000001728047, 2.58768189713716),
             "place": ((-5.268608544477583e-10, 1.099112771305001e-10, 180.0000000520547, -3.962121860770875e-06),
                       (0.5, 0.5, 179.9999999825207, -0.15544595117388804)),
             "angles": (115.72082256222639, 115.7208225622264), "flips": (False, True),
             "rows": ((0.125, 0.0, 0.08472632240172957), (0.25, 0.0, 0.2963272064292811), (0.375, 0.0, 0.5329296346626006),
                      (0.5, 0.0, 0.6922574149016462), (0.0, 0.125, 0.018048739806792852), (0.0, 0.25, 0.06923154596166871),
                      (0.0, 0.375, 0.14694074422497957), (0.0, 0.5, 0.23980115753036843), (0.1, -0.3, 0.14764982987901032),
                      (0.2, -0.6, 0.4575218272023895), (0.303, -0.909, 0.6753695675395203), (0.4, -1.2, 0.7073160723060824),
                      (0.5, -1.5, 0.6913694890648986))},
    "candidate": {"lattice": (4.9666474417274875, 8.4095769676478, 87.74847749549396, 2.5820277490117967),
                  "place": ((-8.4433797942139e-09, -9.589605462579227e-10, 118.88700568626805, -0.00263537545088124),
                            (0.5, 0.5, 118.88700751806567, -0.0026353935304825765)),
                  "angles": (115.3233823251748, 115.32338232517456), "flips": (False, False),
                  "rows": ((0.125, 0.0, 0.03988136142748725), (0.25, 0.0, 0.1479008336393515), (0.375, 0.0, 0.2856862153557611),
                           (0.5, 0.0, 0.40511224011778885), (0.0, 0.125, 0.018728018754091302), (0.0, 0.25, 0.07573881867306073),
                           (0.0, 0.375, 0.16021256852688823), (0.0, 0.5, 0.2588413622091229), (0.1, -0.3, 0.157100121785401),
                           (0.2, -0.6, 0.4367857782869183), (0.303, -0.909, 0.5248112314772494), (0.4, -1.2, 0.39457165057252475),
                           (0.5, -1.5, 0.2935379270257595))},
}
# The rows the fit carries, per packing: a quarter and half a monomer along each lattice vector.
# Eight rows beside sixteen anchors; the other eighteen cells and the six E5 slides are the check.
REGISTRY_FIT_ROWS = ((0.25, 0.0), (0.5, 0.0), (0.0, 0.25), (0.0, 0.5))
TOL["registry"] = 0.2  # relative: a cost within a fifth of DFT's is one tolerance
STAGE11 = {**STAGE7, "registry": True}
# Stage 11's wide bounds on the diameters.  The earlier stages held every diameter within a tenth of
# its valence value; the registry cost is set by hydrogen's contacts and lowering it by a factor
# of three is outside that box, so the box is opened where the cost has its handles (the ridge
# still charges for the distance): hydrogen and the H-F pair may shrink to 0.70, carbon may grow
# to 1.30 (a smooth wall where hydrogen's is a bumpy one), fluorine to 0.85-1.15.
BOUNDS11_WIDE = BOUNDS7.copy()
BOUNDS11_WIDE[0] = (0.90, 1.30)  # x_C
BOUNDS11_WIDE[1] = (0.70, 1.10)  # x_H
BOUNDS11_WIDE[2] = (0.85, 1.15)  # x_F
BOUNDS11_WIDE[8] = (0.70, 1.15)  # x_HF


# Stage 12 (docs/CRYSTAL_STAGE.md 5.25): stage 11 with the pair term's form changed.  Stage 11
# showed the 12-6 wall cannot give the registry cost and the lattice together; against DFT on rigid
# chains moved along and towards each other it is too costly the closer the contact, by a factor
# that grows from 2 to 5 (examples/interchain_frames.py, interchain_wall.py).  The form here is
# Halgren's buffered 14-7 (polyfind.pack.wall_b147) on the same per-element diameters and depths,
# its first buffer the fourteenth variable (0.07 is his; larger is softer).  The belt sites stay
# what they were, a wall alone against hydrogen, in the new form (the trial fits are six times
# worse without them).  Two things follow from the form and are not choices of the fit: the
# diameters mean something else under a softer wall, so their box is wide and their ridge widths
# are three times the earlier stages'; and the torsion terms, fitted with 12-6 between a chain's
# own atoms, have to be fitted again afterwards.
VARIABLES12 = VARIABLES7 + ("wall_delta",)
X0_12 = np.concatenate([X0_7, [0.07]])
BOUNDS12 = np.vstack([BOUNDS7, [(0.0, 0.5)]])
BOUNDS12[0] = (0.90, 1.40)  # x_C
BOUNDS12[1] = (0.80, 1.40)  # x_H
BOUNDS12[2] = (0.85, 1.25)  # x_F
BOUNDS12[8] = (0.70, 1.20)  # x_HF
BOUNDS12[10] = (-5.0, -0.3)  # log10_D_X: the same wall-alone site, in a form where 1e-4 kcal/mol is still felt
SIGMA12 = np.concatenate([SIGMA7, [0.15]])
SIGMA12[[0, 1, 2, 8]] = 0.15
STAGE12 = dict(STAGE11)

# Stage 13 (docs/CRYSTAL_STAGE.md 5.26): stage 12 with two more terms of the charges.  DFT's
# electrostatic potential around single chains (examples/chain_potential.py) says polyfind's
# charges are right for plain PVDF and miss the same 0.4 kcal/mol per e around every chain, which
# is half of all there is around a chain without a dipole; a charge site on each C-H bond and a
# charge drawn along the backbone towards the fluorinated carbons repair most of it.  The
# fifteenth variable is the site's charge and the sixteenth what each fluorine draws
# (polyfind.fitting.FFParameters.bond_sites, backbone_pull); the chains are built with the sites.
# Their prior is what the potential alone gives (-0.318 e and 0.0155 e, with the two bond
# increments at 0.30 and 0.70 of belt4's), not zero: energies cannot set charges (fitted to the
# rigid cells' energies alone they end with no charge on the CF2 groups).
#
# The last three variables are the lattice's torsion triple, the CF2-CH2 type's.  The earlier stages
# held it at the valence fit's.  It cannot stay there: the charges inside a chain have changed, and
# refitted on the conformers under them the triple leaves the alpha-beta gap at -0.24 kcal/mol for
# DFT's -0.73, while along one direction of it (the one-fold term down, the two-fold up) the
# conformers lose 0.03 kcal/mol rms and the gap crosses DFT's.  So the triple is fitted here, with
# the conformers' own error among the rows (the rows from outside the module: evaluate(extra_rows=...),
# examples/stage13_rows.py), and it carries no ridge of its own.  ``VARIABLES13_CHARGES`` is the
# sixteen-variable form without it.
#
# The last two are the F-F pair's own diameter and depth, as multipliers on fluorine's (1 and 1 is
# the mixing rule).  Between rigid chains the reversed-unit candidate's cell turns on fluorine
# against fluorine at 2.5 A, and with the pair left to the mixing rule no wall holds that cell and
# the fully reversed chain's together (5.26).  ``VARIABLES13_TORSION`` is the nineteen-variable
# form without them.
VARIABLES13_CHARGES = VARIABLES12 + ("q_M", "pull_F")
VARIABLES13_TORSION = VARIABLES13_CHARGES + ("V1", "V2", "V3")
VARIABLES13 = VARIABLES13_TORSION + ("x_FF", "D_FF")
SITE_PRIOR = (-0.318, 0.0155)
X0_13_CHARGES = np.concatenate([X0_12, [0.0, 0.0]])
X0_13_TORSION = np.concatenate([X0_13_CHARGES, np.asarray(F.FITTED_VALENCE.torsion, dtype=float)])
X0_13 = np.concatenate([X0_13_TORSION, [1.0, 1.0]])
BOUNDS13_CHARGES = np.vstack([BOUNDS12, [(-0.6, 0.0), (-0.01, 0.05)]])
BOUNDS13_CHARGES[1] = (0.70, 1.40)  # x_H: with the H-F and F-F pairs on their own terms the trial fits want hydrogen itself smaller
BOUNDS13_CHARGES[8] = (0.70, 1.35)  # x_HF
BOUNDS13_TORSION = np.vstack([BOUNDS13_CHARGES, [(-5.0, 1.0), (-1.5, 3.0), (0.0, 3.0)]])
BOUNDS13 = np.vstack([BOUNDS13_TORSION, [(0.80, 1.40), (0.3, 4.0)]])
SIGMA13_CHARGES = np.concatenate([SIGMA12, [0.10, 0.010]])
SIGMA13_TORSION = np.concatenate([SIGMA13_CHARGES, [np.inf, np.inf, np.inf]])
SIGMA13 = np.concatenate([SIGMA13_TORSION, [0.15, 1.0]])
STAGE13 = dict(STAGE12)


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
    if n == 14:
        return VARIABLES12, X0_12, BOUNDS12, SIGMA12
    if n == 16:
        return VARIABLES13_CHARGES, X0_13_CHARGES, BOUNDS13_CHARGES, SIGMA13_CHARGES
    if n == 19:
        return VARIABLES13_TORSION, X0_13_TORSION, BOUNDS13_TORSION, SIGMA13_TORSION
    if n == 21:
        return VARIABLES13, X0_13, BOUNDS13, SIGMA13
    raise ValueError(f"expected 8, 10, 11, 12, 13, 14, 16, 19 or 21 variables, got {n}")


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
    if len(x) >= 13:  # the backbone-angle restraint, log10 of its multiplier on the default
        parts += [x[12:13] / sigma[12:13]]
    if len(x) >= 14:  # the wall's first buffer, from Halgren's 0.07
        parts += [(x[13:14] - X0_12[13]) / sigma[13:14]]
    if len(x) >= 16:  # the bond site's charge and the backbone pull, from what the potential alone gives
        parts += [(x[14:16] - np.array(SITE_PRIOR)) / sigma[14:16]]
    # (the torsion triple of the nineteen-variable form has no ridge: the conformers' rows hold it)
    if len(x) == 21:  # the F-F pair: a diameter multiplier and a depth multiplier, from the mixing rule
        parts += [(x[19:20] - 1.0) / sigma[19:20], np.log(x[20:21]) / sigma[20:21]]
    return np.concatenate(parts)


def ridge(x) -> float:
    d = ridge_parts(x)
    return RIDGE * float(np.sum(d * d))


# Forward-difference steps of the least-squares fitter's Jacobian, per variable (the units of x).
FD_STEPS = {"x_C": 0.005, "x_H": 0.005, "x_F": 0.005, "D_C": 0.02, "D_H": 0.02, "D_F": 0.02, "q_CH": 0.005,
            "q_CF": 0.005, "x_HF": 0.005, "D_HF": 0.02, "log10_D_X": 0.05, "b_q": 0.1, "log10_k_angle": 0.02,
            "wall_delta": 0.005, "q_M": 0.01, "pull_F": 0.002, "V1": 0.05, "V2": 0.05, "V3": 0.05,
            "x_FF": 0.005, "D_FF": 0.05}


@dataclass
class Evaluation:
    x: np.ndarray
    residuals: dict
    objective: float
    ridge: float
    seconds: float
    beta: object = None
    alpha: object = None
    extra: dict | None = None  # settings and side results that are not residuals (pendant tilts)


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
        if len(np.asarray(x)) >= 16:  # stage 13: the chains carry a charge site on every C-H bond
            from .polymers import with_bond_sites
            pol, tag = with_bond_sites(pol), tag + "-bs"
        return (_replace(beta_case, key=f"beta-{tag}", polymer=pol),
                _replace(alpha_case, key=f"alpha-{tag}", polymer=pol))
    return beta_case, alpha_case


def registry_costs(params: F.FFParameters, rows=None, coulomb: str = "ewald", polarizable: bool = True,
                   lj_cutoff: str = "force") -> dict:
    """``{(packing, shear of a, shear of b): (polyfind, DFT)}`` in kcal/mol per monomer: the cost
    under ``params`` of each rigid shear of :data:`REGISTRY` against the unsheared cell, beside
    DFT's on the same atoms.  ``rows`` restricts it to those ``(shear of a, shear of b)`` pairs
    (monomers) in every packing; ``None`` is all thirteen.  The chains are the candidate screen's
    VDF unit whatever :data:`BELT_GEOMETRY` says (the belt-site monomer at the PBE-D3 pendant
    geometry in a copolymer of :data:`polyfind.polymers.PVDF`'s backbone bond length,
    examples/candidate_screen.py ``make_polymer``): the cells DFT was given were built from it."""
    from . import mixedcell as MC
    from . import pack as pack_mod
    from .linegroup import repeat_chains
    from .polarizability import Polarizable
    from .polymers import PVDF, copolymer, monomer_of, pvdf_belt
    want = None if rows is None else {(float(a), float(b)) for a, b in rows}
    out = {}
    with params.applied():
        polymer = params.polymer(copolymer(name="vdf1-belt-pbe", formula="-(CH2-CF2)-", units=[monomer_of(pvdf_belt(BELT_D, "pbe"), "vdf")],
                                           bond_length=PVDF.bond_length))
        P = polymer.bonds_per_repeat
        for pack, data in REGISTRY.items():
            angs = np.asarray(data["angles"], dtype=float)
            ch = repeat_chains(polymer, "TT", np.full((1, P), 180.0), angs)[0]
            # the class is looked up on the module at the call: ``applied`` rebinds it there, and a name imported
            # beforehand would build a packer with the polymer's own charges and the default torsion term
            pk = pack_mod.CrystalPacker(ch, n_chains=2, lj_cutoff=lj_cutoff, coulomb=coulomb, ewald=EwaldSpec() if coulomb == "ewald" else None,
                               polarizable=Polarizable() if polarizable else None)
            st = MC.MixedState(np.array(data["lattice"], dtype=float), np.array(data["place"], dtype=float),
                               np.full((2, P), 180.0), np.tile(angs, (2, 1)), np.zeros(2))
            model = MC.MixedModel(polymer, "TT", pk, np.array(data["flips"], dtype=bool), [ch, ch], angle_stiffness=105.0, shear=True)
            x0 = st.vector()
            c = float(data["lattice"][3])
            e0 = model.physical_energy(x0)
            for fa, fb, dft in data["rows"]:
                if want is not None and (fa, fb) not in want:
                    continue
                x = x0.copy()
                x[model.index("shear")] = [fa * c, fb * c]
                out[(pack, fa, fb)] = ((model.physical_energy(x) - e0) / 2.0, float(dft))
    return out


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
             slope: bool = False, angles: bool = False, pendant_tilts: bool = False, registry: bool = False,
             extra_rows=None) -> Evaluation:
    """Pack and refine beta and alpha under ``x``; return the normalised residuals and the objective.

    ``polarizable`` adds the induced dipoles (:class:`polyfind.polarizability.Polarizable`),
    ``lj_cutoff="force"`` the force-shifted Lennard-Jones term, and ``stiffness`` the relaxed-ion
    C_11, C_22 and C_12 of beta against Dipole's strained-cell values: together the construction
    of the E1/E2 predictions (:data:`STAGE2`).  ``slope`` adds ``d|P|/d eta`` along the polar axis
    (:func:`polar_slope`) against Dipole's three-point Berry pair (:data:`STAGE6`), from the same
    strained states as the stiffness.  ``angles`` adds the refined backbone angles of both
    polymorphs against Dipole's (:data:`STAGE7`), which the refinement already has.
    ``pendant_tilts`` evaluates the alpha-beta gap with each polymorph's CH2/CF2 groups relaxed as
    rigid tilts about their carbons (:func:`polyfind.pendants.tilt_relaxation`, the relaxation
    DFT's geometry has and the rigid-pendant refinement lacks; docs/CRYSTAL_STAGE.md 5.22): the two
    relaxation energies are added to the refined energies before the difference is taken, and
    recorded in ``Evaluation.extra``.  ``registry`` adds the cost of moving neighbouring chains
    along each other (:func:`registry_costs`, the rows :data:`REGISTRY_FIT_ROWS` of each packing)
    against DFT's on the same atoms (:data:`STAGE11`).  ``extra_rows`` is a callable taking the
    parameter set and returning further rows ``{name: (polyfind, reference, z)}``, for anchors that
    live outside this module (the rigid cells and the single-chain potentials of
    examples/interchain_frames.py and examples/chain_potential.py); a failure in it is a parameter
    set that cannot be evaluated, like a failure to pack.
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
        tilt = None
        if pendant_tilts:
            from .pendants import tilt_relaxation
            kw = dict(coulomb=coulomb, ewald=ewald, polarizable=pol, lj_cutoff=lj_cutoff)
            with params.applied():
                tilt = {"alpha": tilt_relaxation(reference_from_prediction(alpha, alpha_case, **kw), alpha_case.polymer),
                        "beta": tilt_relaxation(reference_from_prediction(beta, beta_case, **kw), beta_case.polymer)}
        reg = registry_costs(params, REGISTRY_FIT_ROWS, coulomb=coulomb, polarizable=polarizable, lj_cutoff=lj_cutoff) if registry else None
        more = dict(extra_rows(params)) if extra_rows is not None else {}
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
    extra = None
    if tilt is not None:
        gap += tilt["alpha"].dE_per_monomer - tilt["beta"].dE_per_monomer
        extra = {"pendant_tilts": {s: {"dE_per_monomer": t.dE_per_monomer, "tilts_by_substituents_deg": t.tilts_by_substituents,
                                       "n_iterations": t.n_iterations} for s, t in tilt.items()}}
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
    if reg is not None:
        for (pack, fa, fb), (pred, t) in reg.items():
            res[f"registry {pack} ({fa:g}, {fb:g})"] = (pred, t, (pred - t) / (TOL["registry"] * t))
    for name, (pred, t, z) in more.items():
        if name in res:
            raise ValueError(f"extra row {name!r} has the name of a row of this module")
        res[name] = (float(pred), float(t), float(z))
    rg = ridge(x)
    obj = float(sum(r[2] ** 2 for r in res.values())) + rg
    return Evaluation(np.asarray(x, dtype=float), res, obj, rg, time.time() - t0, beta, alpha, extra)


def table(ev: Evaluation) -> str:
    lines = [f"{'quantity':32s} {'polyfind':>10s} {'DFT':>10s} {'z':>7s}"]
    for name, (pred, target, z) in ev.residuals.items():
        lines.append(f"{name:32s} {pred:10.4f} {target:10.4f} {z:7.2f}")
    lines.append(f"{'ridge':32s} {'':>10s} {'':>10s} {ev.ridge:7.2f}   objective {ev.objective:.2f}   ({ev.seconds:.0f} s)")
    return "\n".join(lines)


def fit(maxfev: int = 240, coulomb: str = "ewald", log=None, x0=None, stage: int = 1,
        checkpoint: str | None = None, method: str = "powell", bounds=None, **eval_kw) -> tuple[Evaluation, list[Evaluation]]:
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
    least_squares asks for itself; the Jacobian's evaluations come on top and are in ``history``.

    ``stage=11`` is stage 7 with the registry rows (:data:`STAGE11`); ``bounds`` replaces the
    stage's box (:data:`BOUNDS11_WIDE` opens the diameters)."""
    from scipy.optimize import least_squares, minimize
    if stage in (2, 3, 4):
        eval_kw = {**STAGE2, **eval_kw}
    elif stage == 6:
        eval_kw = {**STAGE6, **eval_kw}
    elif stage == 7:
        eval_kw = {**STAGE7, **eval_kw}
    elif stage == 11:
        eval_kw = {**STAGE11, **eval_kw}
    elif stage == 12:
        eval_kw = {**STAGE12, **eval_kw}
    elif stage == 13:
        eval_kw = {**STAGE13, **eval_kw}
    start = ({1: X0, 2: X0_2, 3: X0_3, 4: X0_4, 6: X0_4, 7: X0_7, 11: X0_7, 12: X0_12, 13: X0_13}[stage]) if x0 is None else np.asarray(x0, dtype=float)
    if stage in (7, 11, 12, 13) and len(start) == 12:  # from a stage-4/5/6 result: the restraint starts at the default
        start = np.concatenate([start, [0.0]])
    if stage in (12, 13) and len(start) == 13:  # from a 12-6 result: the wall starts at Halgren's buffer
        start = np.concatenate([start, [X0_12[13]]])
    if stage == 13 and len(start) == 14:  # from a stage-12 result: no charge on the sites, nothing drawn
        start = np.concatenate([start, [0.0, 0.0]])
    if stage == 13 and len(start) == 16:  # the torsion triple starts at the valence fit's
        start = np.concatenate([start, X0_13[16:19]])
    if stage == 13 and len(start) == 19:  # the F-F pair starts at the mixing rule
        start = np.concatenate([start, X0_13[19:]])
    bounds = (_config(start)[2] if bounds is None else np.asarray(bounds, dtype=float)).copy()
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
           "ridge_weight": RIDGE, "complete": bool(complete), "evaluations": evaluations,
           "residuals": {k: {"polyfind": v[0], "dft": v[1], "z": v[2]} for k, v in ev.residuals.items()},
           "parameters": {"lj": prm.lj, "lj_pairs": prm.lj_pairs, "charge_increments": prm.charge_increments,
                          "charge_damping": prm.charge_damping},
           "anchors": DFT_ANCHORS, "tolerances": TOL}
    if len(ev.x) >= 12:
        out["belt"] = {"d": BELT_D, "geometry": BELT_GEOMETRY}
    if prm.wall:
        out["parameters"]["wall"] = list(prm.wall)
    if prm.bond_sites or prm.backbone_pull:
        out["parameters"].update(bond_sites=prm.bond_sites, bond_site_fraction=prm.bond_site_fraction, backbone_pull=prm.backbone_pull)
    out["parameters"]["torsion"] = list(prm.torsion)
    if ev.extra:
        out.update(ev.extra)
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1)
