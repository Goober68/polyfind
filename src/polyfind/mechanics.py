"""Electromechanical response of a packed crystal: stiffness, piezoelectricity, actuator figures.

:mod:`polyfind.pack` answers *which* cell is stable and *how polar* it is.  This module
answers what mechanical work it does in a field: the elastic stiffness, the piezoelectric
coefficients, and the blocking stress / free strain / work density an actuator is judged
by.

**What a strain can be in this parametrisation, and what it cannot.**  A cell here is

    a_vec = (a, 0, 0)    b_vec = (b cos gamma, b sin gamma, 0)    c_vec = (0, 0, c)

so ``c_vec`` is pinned to z, which is also the chain axis.  Deforming by a homogeneous
strain ``eps`` sends the rows to ``r_i (1 + eps)``; the third row becomes
``c (eps_xz, eps_yz, 1 + eps_zz)``.  A rigid rotation about z can always put ``a_vec``
back on x -- it is absorbed exactly by the setting angles, and :func:`strained_cell`
folds it into ``phi1, phi2`` -- but nothing can put ``c_vec`` back on z once it has an x
or y component.  Hence

* ``eps_xx``, ``eps_yy``, ``eps_xy`` (Voigt 1, 2, 6) are expressible through
  ``(a, b, gamma)``: :data:`IN_PLANE`;
* ``eps_yz``, ``eps_xz`` (Voigt 4, 5), the two shears that involve the chain axis, are
  **not expressible at all**: :data:`UNREACHABLE`.  ``C_44``, ``C_55``, ``C_45`` and
  every ``C_I4``, ``C_I5`` are therefore not reported and no amount of relaxation would
  recover them -- the cell has no variable for them;
* ``eps_zz`` (Voigt 3) needs ``c``, which is not a cell variable at all: it is the rise
  of the chain's repeat transform, so it moves only if the chain's *shape* does.  With a
  rigid chain it cannot; with valence terms in the kernel and a :class:`Shape` it can, and
  the reachable set becomes :data:`WITH_AXIAL`.  See :func:`axial_report` for the rigid
  diagnosis and :func:`deformable_state` for the deformable measurement.

The reachable stiffness is therefore the 3x3 block ``{11, 22, 12, 16, 26, 66}`` with a
rigid chain, and the 4x4 block that adds ``{33, 13, 23, 36}`` with a deformable one.  Those
are the *true* constants and not a plane-stress reduction: ``C_11 = dsigma_1/deps_1`` is
defined at fixed ``eps_2 ... eps_6``, and fixed ``eps_3`` is exactly what a rigid chain
enforces (and what the constraint enforces when it does not).  What cannot be had from the
3x3 is the *compliance*: ``S = C^-1`` of that block is the compliance at clamped
``eps_zz``, so every rigid-path ``d`` below is an axially clamped piezoelectric strain
coefficient.  On the deformable path ``S`` is the inverse of the 4x4 and ``d`` is clamped
only in the two shears the cell cannot express at all.

**Axial strain: two paths, and only one of them is a measurement.**  ``c`` can always be
*changed* -- bond lengths are rigid but backbone bond angles are not, and opening every
backbone angle of PE by 8 degrees lengthens the repeat by 4.5% (PVDF alpha and gamma by
7.3%).  What decides whether that is a measurement is whether anything in the energy
resists it.

*Rigid path* (``shape=None``, the default, and the only thing available under the
illustrative potential).  The packing kernel is Lennard-Jones + damped-shifted-force
Coulomb + a Fourier torsion term and has **no bond or angle term**;
:meth:`polyfind.fitting.FFParameters.applied` does not forward the fitted valence terms
into packing, on the grounds that packing holds the chain rigid.  The only thing left
resisting an angle is :func:`polyfind.refine.refine_crystal`'s own harmonic restraint,
whose stiffness (``angle_stiffness = 105`` kcal/mol/rad^2, "UFF-like for sp3 carbon") is an
invented number deliberately excluded from the reported lattice energy.  A ``C_33``
obtained that way is a report of that constant: measured, it is exactly affine in it, 79%
of beta's and 84% of PE's coming from it, and the refined structure is stationary in ``c``
at ``k = 105`` and nowhere else.  :func:`axial_report` measures all of that, returns
``computable=False``, and :func:`electromechanical_response` leaves ``C_33`` as ``nan``.

*Deformable path* (``shape=`` a :class:`Shape`, which needs
``CrystalPacker(..., valence=...)``).  The fitted stretch and bend terms of
``docs/VALENCE_FIT.md`` are forwarded into the lattice kernel, which is opt-in and is the
one change that makes the axial direction real.  The chain's line-group parameters are then
relaxed at every strain state, subject to an equality constraint on the repeat, and the
axial stress is that constraint's Lagrange multiplier.  The evidence that the contamination
is gone is the same sweep run again: with valence terms, ``C_33`` does not move at all as
``angle_stiffness`` goes 0 -> 210 (six figures for beta, PE and alpha), because the
relaxation puts back whatever the refinement's restraint moved.  What remains is stated
rather than hidden: bond lengths are still rigid, so the stretch channel -- roughly half of
a real chain modulus -- is absent, and the line group restricts the chain to symmetric
conformations, so the number is an **upper bound** on both counts.

One thing is not fixed by any of this: a ``dE/dc`` taken at *fixed chain shape* is not an
axial stress in any physical sense, because the only way a rigid chain's ``c`` can grow is
for the repeats to slide apart across the periodic boundary and the nonbonded sum excludes
exactly the bonded pairs that cross it.  :func:`stress` returns it as ``sigma[2]`` on the
rigid path and says so; :func:`deformable_state` replaces it.

**Sign and unit conventions.**

* Engineering (Voigt) strain ``eps = (e_xx, e_yy, e_zz, 2 e_yz, 2 e_xz, 2 e_xy)``,
  dimensionless, tension positive.
* ``sigma_J = (1/V0) dE/de_J`` in GPa, tension positive.
* The electric enthalpy per unit reference volume is ``h(eps, E) = u(eps) - E . P``,
  which is what the kernel minimises (it adds ``-mu . E`` to the cell energy).  Then
  ``sigma_J = dh/de_J``, ``P_i = -dh/dE_i``, and

      e_iJ = (1/V0) dmu_i/de_J |_E  =  -dsigma_J/dE_i |_eps   (C/m^2)
      d_iJ = de_J/dE_i |_sigma=0    =  (e S)_iJ,  S = C^-1    (pC/N = pm/V)

  ``e`` is the **proper** piezoelectric constant: the dipole per *reference* volume, not
  ``dP/de`` with ``P = mu/V``.  The two differ by ``P_i`` on the three diagonal columns,
  and for a rigid dipole array that difference is the whole of the naive ``dP/de_xx``
  (dilating a box of fixed dipoles changes ``P`` only through ``V``, and does nothing
  piezoelectric at all).  :class:`Piezoelectric` reports both.  Both routes to ``d`` are
  computed independently here -- ``e`` from the dipole of strained, internally relaxed
  cells; ``d`` from the strain of field-relaxed cells -- and :class:`Response` reports how
  far apart they land.  One trap on the way: a shear leaves ``a_vec`` off the x axis, so
  the packer's canonical frame is rotated with respect to the reference by
  ``theta = -e_xy/2``; the dipole has to be rotated back out of it and the applied field
  rotated into it, and each of those terms is the same size as the coefficient itself.
* An actuator at field ``E``: free strain ``eps_free = d^T E`` (zero in-plane stress),
  blocking stress ``sigma_block = C eps_free = e^T E``, the stress the crystal exerts on
  a rigid clamp (the external stress needed to *hold* ``eps = 0`` is its negative), and
  work density ``w = (1/4) sigma_block . eps_free`` against a matched linear load; the
  full triangle ``(1/2) sigma_block . eps_free`` is reported alongside.

**"Relaxed-ion" means relaxed chains, not relaxed atoms.**  The internal degrees of
freedom of this parametrisation are the rigid-body ones: each chain's setting angle
``phi1, phi2`` and the axial offset ``dz``.  They are relaxed at every strain state, by
:func:`polyfind.pack.polish` on the analytic cell gradient, which is what makes a sweep
over strain states affordable.  With ``shape=None`` the chain conformation is held rigid,
so these are *rigid-chain* relaxed-ion constants: softer than clamped-ion, stiffer than
constants whose torsions also relax.  With a :class:`Shape` the conformation relaxes too,
and ``eps_zz`` stops being an identity and becomes an explicit equality constraint on the
repeat -- which is the only honest way to let the chain move and still say what strain the
calculation is at.

Nothing here is a prediction of experiment.  The illustrative potential is not quantitative
(DESIGN.md section 5.4) and the fitted presets are fitted to single-chain DFT energies and
forces, with no elastic or piezoelectric data in the objective at all.
``docs/ELECTROMECHANICS.md`` puts published PVDF and polyethylene numbers beside these, says
what could not be found to compare against, and records one further limit worth knowing
before reading any coefficient here.  With fixed point charges on a *rigid* chain the only
piezoelectric channel in the model is chain reorientation, so every diagonal column of ``e``
is zero.  Letting the chain deform does not by itself fix that: for a planar all-trans
zigzag with bond-charge-increment charges the cell dipole is a sum over bonds of fixed
length whose bisectors the zigzag's mirror symmetry pins perpendicular to the axis, so it is
*exactly* independent of the backbone angle -- beta-PVDF's and PE's diagonal columns stay
zero for that reason and not for the old one.  The helices are different: alpha and gamma
have three and five shape parameters and do show a diagonal response.

**The third path removes that last obstacle and brings its own.**
``CrystalPacker(charge_flux=...)`` lets the increments themselves depend on the backbone
angle (:class:`polyfind.forcefield.FluxTopology`), which is what a planar zigzag needs and
what gives beta a non-zero ``d_33`` and ``d_31`` -- with the measured signs and about an
order of magnitude short.  Nothing in this module branches on it: the flux reaches every
stress, every dipole derivative and every gradient through the kernel, and the agreement
between the two independent routes to ``d`` is what says so.  What it does *not* come with
is a trustworthy magnitude; ``polyfind.fitting.FITTED_VALENCE_FLUX`` and
``docs/ELECTROMECHANICS.md`` section 5.6 say how the two coefficients were fitted, how
poorly, and which of the resulting numbers are signs rather than values.
``examples/electromechanics.py`` runs the whole thing on the reference polymorphs.

**The fourth path gives the crystal a dielectric constant.**  ``CrystalPacker(polarizable=...)``
(:mod:`polyfind.polarizability`, Ewald only) adds self-consistent induced dipoles; nothing
in this module branches on it either -- the dipole the strain states differentiate is the
packer's total one and the stress carries the fixed-dipole gradient that the solve's
stationarity makes exact -- and :func:`dielectric_tensor` is the known-answer check that
must pass before any piezoelectric number with it is read.  ``docs/ELECTROMECHANICS.md``
section 5.8 has that check and what the dipoles do to ``d_33`` and ``d_31``:
a quarter of one shortfall and a sixteenth of the other, not a closure.
"""
from __future__ import annotations

import time
from dataclasses import dataclass, field as dc_field

import numpy as np

from .pack import EV_TO_KCAL, E_PER_A2_TO_C_PER_M2, PeriodicChain, polish

# --- units ---------------------------------------------------------------------------
# 1 kcal/mol / A^3 = (4184 / 6.02214076e23) J / 1e-30 m^3 = 6.9477 GPa.
KCAL_MOL_A3_TO_GPA = 4184.0 / 6.02214076e23 / 1e-30 / 1e9
V_PER_A_TO_V_PER_M = 1e10  # the kernel's field unit
# e (C/m^2) times S (1/GPa) is a d in C/N; 1 C/N = 1e12 pC/N.
C_PER_M2_PER_GPA_TO_PC_PER_N = 1e-9 * 1e12
GPA_TO_KJ_PER_M3 = 1e9 / 1e3  # sigma . eps in GPa -> kJ/m^3

# --- which strain components this cell can express ------------------------------------
VOIGT_LABELS = ("xx", "yy", "zz", "yz", "xz", "xy")
IN_PLANE = (0, 1, 5)  # Voigt 1, 2, 6: reachable through (a, b, gamma)
AXIAL = 2  # Voigt 3: needs c, the chain's own repeat
WITH_AXIAL = (0, 1, 2, 5)  # Voigt 1, 2, 3, 6: the reachable set once the chain can deform
UNREACHABLE = (3, 4)  # Voigt 4, 5: c_vec would have to leave z

# Voigt index -> the strain tensor it switches on (engineering convention).
_EPS_BASIS = np.zeros((6, 3, 3))
for _k, (_i, _j) in enumerate([(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)]):
    _w = 1.0 if _i == _j else 0.5
    _EPS_BASIS[_k, _i, _j] = _w
    _EPS_BASIS[_k, _j, _i] = _w
_EPS_BASIS.setflags(write=False)


def strain_tensor(eps6) -> np.ndarray:
    """(6,) engineering Voigt strain -> the symmetric (3, 3) tensor."""
    return np.einsum("k,kij->ij", np.asarray(eps6, dtype=float).reshape(6), _EPS_BASIS)


def lattice_matrix(a: float, b: float, gamma: float, c: float) -> np.ndarray:
    """The packer's lattice vectors as the rows of a (3, 3) matrix (A)."""
    g = np.deg2rad(gamma)
    return np.array([[a, 0.0, 0.0], [b * np.cos(g), b * np.sin(g), 0.0], [0.0, 0.0, c]])


def strain_from_cell(a0, b0, gamma0, a, b, gamma) -> np.ndarray:
    """The in-plane Voigt strain that takes cell ``(a0, b0, gamma0)`` to ``(a, b, gamma)``.

    Both cells are in the packer's canonical frame (``a_vec`` on x), so they differ by an
    unknown rigid rotation about z as well as by the strain.  With ``A`` the 2x2 matrix of
    in-plane lattice rows, ``A = A0 (1 + eps) R^T``, so ``M = A0^-1 A`` satisfies
    ``M M^T = (1 + eps)^2`` and the rotation drops out of the symmetric square root.  The
    inverse of :func:`strained_cell` restricted to :data:`IN_PLANE`.
    """
    def two_d(al, bl, gl):
        g = np.deg2rad(gl)
        return np.array([[al, 0.0], [bl * np.cos(g), bl * np.sin(g)]])

    M = np.linalg.solve(two_d(a0, b0, gamma0), two_d(a, b, gamma))
    w, V = np.linalg.eigh(M @ M.T)
    F = V @ np.diag(np.sqrt(np.maximum(w, 1e-300))) @ V.T  # the symmetric 1 + eps
    eps = np.zeros(6)
    eps[0], eps[1], eps[5] = F[0, 0] - 1.0, F[1, 1] - 1.0, F[0, 1] + F[1, 0]
    return eps


@dataclass
class StrainedCell:
    """A reference cell deformed by a Voigt strain, back in the packer's canonical frame.

    ``params`` is the 7-vector :meth:`polyfind.pack.CrystalPacker.energy` takes and ``c``
    the repeat; ``theta`` (deg) is the rigid rotation about z that put ``a_vec`` back on
    x, already folded into ``phi1`` and ``phi2``.  ``dparams`` (6, 7) and ``dc`` (6) are
    the derivatives with respect to the six Voigt components, ``nan`` for
    :data:`UNREACHABLE`; :func:`stress` contracts the kernel's analytic gradient with
    them.
    """

    params: np.ndarray
    c: float
    theta: float
    dparams: np.ndarray
    dc: np.ndarray

    def _rot(self, v, sign: float) -> np.ndarray:
        t = sign * np.deg2rad(self.theta)
        ct, st = np.cos(t), np.sin(t)
        v = np.asarray(v, dtype=float)
        return np.stack([ct * v[..., 0] - st * v[..., 1], st * v[..., 0] + ct * v[..., 1], v[..., 2]], axis=-1)

    def unrotate(self, v) -> np.ndarray:
        """A vector the packer reports in the canonical frame, back in the reference frame.

        The canonical frame is the reference frame turned by ``theta`` about z (that is what
        put ``a_vec`` back on x), so a polarization or a dipole read out of the packer has to
        be turned back by ``-theta`` before it can be differentiated against the strain --
        otherwise a pure shear looks as though it rotated the polarization, and the spurious
        term is first order in the strain.
        """
        return self._rot(v, -1.0)

    def rotate_field(self, e_lab) -> np.ndarray:
        """A field given in the reference (lab) frame, as the packer's canonical frame sees it.

        The mirror image of :meth:`unrotate`, and just as load bearing: the strain that
        defines a piezoelectric coefficient is the *symmetric* deformation, with the field
        components referred to the reference axes.  Handing the packer the lab field
        unrotated instead amounts to turning the crystal by ``-theta`` inside the field,
        which for a shear adds ``(dtheta/de_6) (z_hat x P) . E`` to the energy -- exactly
        ``P/2``, the same size as the coefficient being measured.
        """
        return self._rot(e_lab, +1.0)


def strained_cell(params0, c0: float, eps6) -> StrainedCell:
    """Apply a homogeneous Voigt strain to a reference cell.

    Raises for any strain with a non-zero :data:`UNREACHABLE` component: those deformations
    take ``c_vec`` off the z axis and this cell has no way to say so.
    """
    p0 = np.asarray(params0, dtype=float).reshape(7)
    e = np.asarray(eps6, dtype=float).reshape(6)
    bad = [k for k in UNREACHABLE if e[k] != 0.0]
    if bad:
        names = ", ".join(f"eps_{VOIGT_LABELS[k]} (Voigt {k + 1})" for k in bad)
        raise ValueError(
            f"{names} cannot be applied: this cell fixes c_vec along z (the chain axis), and a "
            "shear involving the chain axis would tilt it off z.  There is no cell variable for "
            "it, so C_44, C_55, C_45 and the mixed constants that use them are outside what this "
            "parametrisation can express at all (see the module docstring)."
        )
    a0, b0, gam0, phi1, phi2, dz0, flip = p0
    H0 = lattice_matrix(a0, b0, gam0, c0)
    H = H0 @ (np.eye(3) + strain_tensor(e))
    r1, r2, r3 = H
    if abs(r3[0]) > 1e-12 or abs(r3[1]) > 1e-12:  # unreachable was screened above
        raise AssertionError("c_vec left the z axis")
    a, b, c = (float(np.linalg.norm(r)) for r in (r1, r2, r3))
    cosg = float(r1 @ r2) / (a * b)
    gam = float(np.degrees(np.arccos(np.clip(cosg, -1.0, 1.0))))
    if float(np.cross(r1, r2)[2]) <= 0.0:
        raise ValueError("the strained cell is degenerate or left-handed; the strain is too large")
    theta = float(-np.arctan2(r1[1], r1[0]))  # rotate the crystal by theta to put a_vec on x
    th_deg = np.degrees(theta)
    params = np.array([a, b, gam, phi1 + th_deg, phi2 + th_deg, dz0 * c / c0, flip])

    dH = np.einsum("ij,kjl->kil", H0, _EPS_BASIS)  # (6, 3, 3)
    dr1, dr2, dr3 = dH[:, 0, :], dH[:, 1, :], dH[:, 2, :]
    da, db, dc = dr1 @ r1 / a, dr2 @ r2 / b, dr3 @ r3 / c
    dcos = (dr1 @ r2 + dr2 @ r1) / (a * b) - cosg * (da / a + db / b)
    sing = float(np.sqrt(max(1.0 - cosg * cosg, 1e-30)))
    dgam = np.degrees(-dcos / sing)
    dth = np.degrees(-(r1[0] * dr1[:, 1] - r1[1] * dr1[:, 0]) / (r1[0] ** 2 + r1[1] ** 2))
    dparams = np.zeros((6, 7))
    dparams[:, 0], dparams[:, 1], dparams[:, 2] = da, db, dgam
    dparams[:, 3] = dparams[:, 4] = dth
    dparams[:, 5] = dz0 * dc / c0
    for k in UNREACHABLE:
        dparams[k] = np.nan
        dc[k] = np.nan
    return StrainedCell(params=params, c=c, theta=float(th_deg), dparams=dparams, dc=dc)


# --- the reference state ---------------------------------------------------------------
def _packer_class():
    """``pack.CrystalPacker`` looked up late, so a :meth:`FFParameters.applied` block wins."""
    from . import pack as pack_mod

    return pack_mod.CrystalPacker


@dataclass
class Reference:
    """A packed, refined crystal ready to be strained.

    ``packer`` carries the potential and the (rigid) chain; ``params`` is its cell,
    relaxed to zero in-plane stress by :func:`relax_reference` before any strain is
    applied.  ``volume`` is the reference cell volume in A^3, the ``V0`` every stress and
    stiffness below is divided by.
    """

    packer: object
    params: np.ndarray
    c: float
    label: str = ""

    @property
    def volume(self) -> float:
        return float(self.packer.cell_volume(self.params[None], np.array([self.c]))[0])

    def polarization(self) -> np.ndarray:
        """Reference polarization in C/m^2 (zeros for a chain with no charges)."""
        return np.asarray(self.packer.polarization(self.params[None], c=np.array([self.c]))[0], dtype=float)


def reference_from_chain(chain: PeriodicChain, params, n_chains: int = 2, label: str = "", **packer_kw) -> Reference:
    """A :class:`Reference` for ``chain`` at cell ``params`` (7-vector, e.g. ``PackResult.params``)."""
    packer = _packer_class()(chain, n_chains=n_chains, **packer_kw)
    return Reference(packer=packer, params=np.asarray(params, dtype=float).reshape(7), c=float(chain.c), label=label)


def reference_from_result(chain: PeriodicChain, result, **kw) -> Reference:
    """A :class:`Reference` from a :class:`~polyfind.pack.PackResult` and the chain it describes."""
    return reference_from_chain(chain, result.params, n_chains=result.n_chains, label=result.chain, **kw)


def refined_reference(polymer, seq, states=None, label: str = "", n_chains: int = 2,
                      pack_kw: dict | None = None, refine_kw: dict | None = None, **packer_kw):
    """Pack and refine one polymorph, then wrap it as a :class:`Reference`.

    Returns ``(reference, refine_result)``.  The chain is rebuilt from the refined torsions
    and bond angles *aligned to the same frame the refinement used*, which is the one
    subtlety worth hiding behind a function: the frame
    :func:`polyfind.linegroup.repeat_chains` picks is a principal axis of the cross section
    and flips discontinuously, so a chain rebuilt without ``align_to`` has setting angles
    that no longer match the refined ``phi1, phi2``.
    """
    from .pack import pack, periodic_chain
    from .polymers import THREE_STATE
    from .refine import refine_crystal
    from .linegroup import repeat_chains

    states = THREE_STATE if states is None else states
    chain = periodic_chain(polymer, seq, states)
    res = pack(chain, n_chains=n_chains, **(pack_kw or {}))[0]
    rr = refine_crystal(polymer, res, **(refine_kw or {}))
    angles0 = np.array([polymer.backbone[k].backbone_angle for k in range(polymer.bonds_per_repeat)])
    frame = repeat_chains(polymer, res.chain, np.asarray(res.dihedrals, dtype=float)[None], angles0)[0]
    built = repeat_chains(polymer, res.chain, rr.torsions[None], rr.angles, align_to=frame)[0]
    return reference_from_chain(built, rr.result.params, n_chains=n_chains,
                                label=label or res.chain, **packer_kw), rr


_CELL_FREE = [0, 1, 2, 3, 4, 5]
_INTERNAL_FREE = [3, 4, 5]  # phi1, phi2, dz: the rigid-body internal coordinates


def _bounds(ref: Reference, params, free):
    """``lo, hi`` for :func:`polyfind.pack.polish`; the periodic variables get a full period."""
    p = np.asarray(params, dtype=float)
    lo = p[:6].copy()
    hi = p[:6].copy()
    if 0 in free:
        lo[0], hi[0] = 0.6 * p[0], 1.8 * p[0]
    if 1 in free:
        lo[1], hi[1] = 0.6 * p[1], 1.8 * p[1]
    if 2 in free:
        lo[2], hi[2] = 55.0, 125.0
    # The three periodic variables get a full period *centred on where they start*, so
    # polish's wrap brings them back next to their starting value instead of jumping to the
    # other end of the interval -- a difference of one period in dz is invisible in an
    # energy and fatal in a finite difference over strain.
    cz = float(ref.packer.chain.c)
    lo[3], hi[3] = p[3] - 180.0, p[3] + 180.0
    lo[4], hi[4] = p[4] - 180.0, p[4] + 180.0
    lo[5], hi[5] = p[5] - 0.5 * cz, p[5] + 0.5 * cz
    return lo, hi


def relax(ref: Reference, params, free=_INTERNAL_FREE, maxiter: int = 200) -> np.ndarray:
    """Minimise the cell energy over ``free`` (default: the internal coordinates only).

    :func:`polyfind.pack.polish` stops at ``gtol = 1e-5``, which is a loose tolerance to be
    taking *second* derivatives through -- a leftover gradient enters the stress difference
    divided by the strain step.  It was checked rather than assumed: restarting L-BFGS-B from
    its own answer, once or three more times, moves no constant of beta-PVDF or PE by 0.002
    GPa even at ``step = 5e-4``, where a residual would show most.  The finite difference here
    is limited by the nonbonded cutoff, not by the relaxation; see :class:`Elastic`.
    """
    p = np.asarray(params, dtype=float).reshape(7)
    free = list(free)
    lo, hi = _bounds(ref, p, free)
    return polish(ref.packer, p, lo, hi, free, method="lbfgs", gradient="analytic", maxiter=maxiter)


def relax_reference(ref: Reference, maxiter: int = 300) -> Reference:
    """Relax ``(a, b, gamma)`` and the internal coordinates, so the in-plane stress is zero.

    A second derivative is an elastic constant only about a stress-free state.  This does
    not touch ``c`` -- it cannot -- so the *axial* stress is whatever the refinement left
    (see :func:`axial_report`), and the constants below are the in-plane ones at fixed
    ``eps_zz``.
    """
    return Reference(packer=ref.packer, params=relax(ref, ref.params, _CELL_FREE, maxiter), c=ref.c, label=ref.label)


def stress(ref: Reference, sc: StrainedCell | None = None, params=None,
           rigid_chain: bool = True) -> tuple[float, np.ndarray]:
    """``(energy, sigma)`` for a strained cell: ``sigma_J = (1/V0) dE/de_J`` in GPa.

    The kernel's analytic ``dE/d(a, b, gamma, phi1, phi2, dz)`` and ``dE/dc`` are
    contracted with :attr:`StrainedCell.dparams` and :attr:`StrainedCell.dc`, so no extra
    energy evaluation is spent and the derivative is closed-form all the way from the
    strain.  ``sigma[3]``, ``sigma[4]`` are ``nan``: :data:`UNREACHABLE`.

    ``rigid_chain=True`` (the default) evaluates the packer's own chain at the strained
    ``c``, which is the rigid-chain reading: ``sigma[2]`` is then **not** an axial stress in
    any useful sense, because at a fixed chain shape the only way ``c`` can grow is to slide
    the repeats apart across the periodic boundary and the bonded exclusions remove exactly
    the pairs that cross it.  ``rigid_chain=False`` says the packer already carries the
    *deformed* chain whose own repeat is ``sc.c`` -- :func:`deformable_state` puts it there
    -- so the kernel is asked for that chain's energy with no override.  ``sigma[0]``,
    ``sigma[1]``, ``sigma[5]`` are then the true total stresses (the shape response drops out
    by the envelope theorem: the shape is relaxed and the constraint that fixes ``c`` does
    not move under an in-plane strain), while ``sigma[2]`` still is not, and
    :func:`deformable_state` replaces it with the constraint multiplier.
    """
    if sc is None:
        sc = strained_cell(ref.params, ref.c, np.zeros(6))
    p = sc.params if params is None else np.asarray(params, dtype=float).reshape(7)
    if rigid_chain:
        c_arr = None if abs(sc.c - ref.packer.chain.c) < 1e-12 else np.array([sc.c])
        E, g_cell, _, g_c = ref.packer.energy_and_grad(
            p, coords=None if c_arr is None else ref.packer.chain.coords, c=c_arr)
    else:
        E, g_cell, _, g_c = ref.packer.energy_and_grad(p)
    V0 = ref.volume
    sig = np.full(6, np.nan)
    for K in range(6):
        if K in UNREACHABLE:
            continue
        sig[K] = (float(g_cell @ sc.dparams[K, :6]) + float(g_c) * float(sc.dc[K])) / V0
    return float(E), sig * KCAL_MOL_A3_TO_GPA


@dataclass
class StrainState:
    """One strained, internally relaxed configuration."""

    cell: StrainedCell
    params: np.ndarray
    energy: float
    stress: np.ndarray  # (6,) GPa, nan where unreachable
    m: np.ndarray  # dipole per *reference* volume, C/m^2, in the reference frame
    polarization: np.ndarray  # dipole per current volume, C/m^2, in the reference frame
    relaxed: object = None  # the :class:`Relaxed` behind it, on the deformable path only


def _state_at_strain(ref: Reference, eps6, relax_internal: bool = True, field_lab=None) -> StrainState:
    """Strain the reference, relax the internal coordinates, and read off stress and dipole.

    ``field_lab`` is a uniform field given in the *reference* frame; it is rotated into the
    packer's canonical frame by :meth:`StrainedCell.rotate_field` before being applied, and
    the packer's field is restored afterwards.
    """
    sc = strained_cell(ref.params, ref.c, eps6)
    saved = ref.packer.field
    e_pk = None if field_lab is None else sc.rotate_field(np.asarray(field_lab, dtype=float))
    try:
        if e_pk is not None:
            ref.packer.set_field(e_pk)
        p = relax(ref, sc.params) if relax_internal else sc.params
        E, sig = stress(ref, sc, params=p)
        c_arr = np.array([sc.c])
        mu = np.asarray(ref.packer.dipole(p[None], c=c_arr)[0], dtype=float)
    finally:
        ref.packer.set_field(saved)
    if e_pk is not None:
        # The kernel's field is a *canonical-frame* vector, and the canonical frame turns
        # with the strain, so the strain derivative of ``-mu . E_pk`` has a second piece:
        # the field itself rotates, by ``dtheta/de_K`` about z.  It is not a correction --
        # for a shear it is ``P/2``, the size of the coefficient being measured -- and
        # without it the zero-stress strain below disagrees with ``e S`` by exactly that.
        zxE = np.array([-e_pk[1], e_pk[0], 0.0])
        extra = -EV_TO_KCAL * float(mu @ zxE) / ref.volume * KCAL_MOL_A3_TO_GPA
        for K in range(6):
            if K not in UNREACHABLE:
                sig[K] += extra * np.deg2rad(sc.dparams[K, 3])
    m = sc.unrotate(mu) / ref.volume * E_PER_A2_TO_C_PER_M2
    pol = sc.unrotate(mu) / float(ref.packer.cell_volume(p[None], c_arr)[0]) * E_PER_A2_TO_C_PER_M2
    return StrainState(cell=sc, params=p, energy=E, stress=sig, m=m, polarization=pol)


# --- the deformable chain -----------------------------------------------------------------
@dataclass
class Shape:
    """The chain's internal degrees of freedom: what makes ``eps_zz`` and a dilation real.

    Everything above this point holds the chain rigid, so ``c`` is a constant and a
    dilation cannot change the cell dipole.  A :class:`Shape` is the same parametrisation
    :func:`polyfind.refine.refine_crystal` uses -- the line group of the repeat
    (:mod:`polyfind.linegroup`), whose free parameters are the pattern's torsions and the
    backbone bond angles with the closure condition solved out, so the repeat stays exactly
    periodic and ``c`` is a smooth function of the parameters.

    Two things are needed before it means anything.  The **packer must carry valence
    terms** (``CrystalPacker(..., valence=...)``): without them nothing in the lattice
    energy resists opening an angle, and a relaxation over these parameters walks straight
    to a bound.  And the **frame** must be the one the reference was refined in, because
    ``phi1, phi2`` are measured against it; every chain built here is aligned to it.

    ``step`` is a *geometry-only* finite difference: the chain build (NeRF, closure solve,
    Kabsch alignment) is differenced, never the energy, so a displacement costs no kernel
    row -- the same trick :func:`polyfind.refine.refine_crystal` uses for its analytic
    gradient.
    """

    polymer: object
    lg: object  # polyfind.linegroup.LineGroup
    frame: PeriodicChain
    name: str
    step: float = 0.05  # deg
    penalty: float = 5.0  # on the commensurability residual, as in refine_crystal
    max_torsion_change: float = 40.0
    max_angle_change: float = 8.0

    @property
    def n(self) -> int:
        return int(self.lg.n_params)

    @property
    def x0(self) -> np.ndarray:
        return np.zeros(self.n)

    def limits(self) -> np.ndarray:
        return np.where(self.lg.is_angle(), self.max_angle_change, self.max_torsion_change)

    def labels(self) -> list:
        return list(self.lg.labels())

    def chains(self, X) -> list:
        """One chain per row of ``X`` (M, n), all aligned to :attr:`frame`."""
        return self.lg.chains(np.atleast_2d(np.asarray(X, dtype=float)), name=self.name, align_to=self.frame)

    def batch(self, x) -> list:
        """``[centre] + [+h, -h for every parameter]``, built in one pass."""
        x = np.asarray(x, dtype=float).reshape(self.n)
        rows = [x]
        for j in range(self.n):
            for sgn in (1.0, -1.0):
                y = x.copy()
                y[j] += sgn * self.step
                rows.append(y)
        return self.chains(np.array(rows))

    def conformation(self, x) -> tuple:
        """``(torsions, backbone angles)`` of one parameter vector."""
        tors, angs = self.lg.expand(np.asarray(x, dtype=float).reshape(1, self.n))
        return tors[0], angs[0]


def shape_of(polymer, ref: Reference, torsions=None, angles=None, frame=None) -> Shape:
    """The :class:`Shape` of a reference chain, at its own torsions and backbone angles.

    ``angles`` must be the **refined** backbone angles (``RefineResult.angles``), not the
    polymer's ideal ones: the reference chain was built from the refined ones, and a line
    group built at the wrong angles describes a different chain.  :func:`refined_reference`
    returns both, so ``shape_of(polymer, ref, angles=rr.angles)`` is the usual call.

    ``frame`` defaults to the reference's own chain, which is the frame ``phi1, phi2`` and
    ``dz`` are measured against.  Anything else silently rotates every chain built here
    relative to the cell the reference describes; the relaxation would absorb it (the
    setting angles are free) but the starting point would be wrong for no reason.
    """
    from .linegroup import line_group

    chain = ref.packer.chain
    tors = np.asarray(chain.dihedrals if torsions is None else torsions, dtype=float)
    B = polymer.bonds_per_repeat
    angs = (np.array([polymer.backbone[k].backbone_angle for k in range(B)], dtype=float)
            if angles is None else np.asarray(angles, dtype=float))
    lg = line_group(polymer, chain.name, tors, states=polymer.states, bond_angles=angs)
    return Shape(polymer=polymer, lg=lg, frame=chain if frame is None else frame, name=chain.name)


@dataclass
class Relaxed:
    """One relaxed configuration of cell *and* chain shape, and how well it converged.

    ``multiplier`` is ``dE/dc`` of the constrained minimum, in kcal/(mol A): the Lagrange
    multiplier of the constraint ``c(shape) = c_target``.  It is the axial stress, and it
    comes out of the stationarity condition rather than out of an extra relaxation --
    ``grad_shape = multiplier * dc_dshape`` at the solution, so ``residual`` (the part of
    the shape gradient that is *not* along the constraint, relative to the whole) measures
    how converged the constrained minimisation is.  ``c_error`` is how far the constraint
    itself is from satisfied and ``cell_residual`` the largest free cell gradient.

    **Read ``residual`` next to ``shape_gradient``, never alone.**  It is a *ratio*, so
    where the shape barely moves -- the reference itself, and any strain state whose
    symmetry forbids a conformational response -- both numerator and denominator are at the
    numerical floor and the ratio means nothing.  With one shape parameter it is identically
    zero for the opposite reason: the gradient and the constraint are parallel by
    construction.  It is informative exactly where there is a response to be converged, and
    there it runs at 1e-8 (alpha's axial column) to 1e-3 (alpha's in-plane columns).
    """

    params: np.ndarray
    x: np.ndarray
    chain: PeriodicChain
    energy: float
    multiplier: float
    residual: float
    c_error: float
    cell_residual: float
    iterations: int
    shape_gradient: float = 0.0  # |dE/d(shape)| at the solution, kcal/(mol deg)


def _shape_scalar(chain, gX, gc, packer, penalty: float) -> float:
    """``gX . coords + gc c + E_torsion + penalty`` for one displaced chain.

    The lattice energy *and* the valence energy depend on a conformation only through
    ``(coords, c)``, and the kernel returns both partial derivatives exactly, so this
    scalar has the same shape-gradient as the objective does at the point ``gX, gc`` were
    taken.  Differencing it needs the displaced *geometry* and not its energy, which is why
    a shape gradient costs no kernel row.
    """
    return (float((gX * chain.coords).sum()) + gc * chain.c
            + packer.torsion_energy(chain.dihedrals) * packer.n_chains
            + penalty * chain.rotation_error ** 2)


def relax_deformable(ref: Reference, shape: Shape, params, x=None, free=_INTERNAL_FREE,
                     c_target: float | None = None, maxiter: int = 120) -> Relaxed:
    """Minimise the energy over the free cell parameters *and* the chain's shape.

    ``c_target`` (A), when given, is imposed as an equality constraint on the repeat --
    SLSQP, with the analytic objective gradient and the constraint's own geometry-only
    Jacobian.  That is what makes a strain state well posed: an in-plane strain has to hold
    ``eps_zz = 0``, which with a deformable chain is a constraint and no longer an
    identity, and an axial strain *is* the constraint.  ``c_target=None`` leaves ``c`` free
    (L-BFGS-B), which is how the reference is brought to zero axial stress.

    Raises if the packer has no valence terms: without them the objective has no restoring
    force in the shape directions at all and the answer would be a report of the bounds.
    """
    from scipy.optimize import minimize

    packer = ref.packer
    if getattr(packer, "_valence", None) is None:
        raise ValueError(
            "relaxing the chain shape needs a packer with valence terms: the lattice energy on "
            "its own has no bond or angle term, so nothing resists opening an angle and the "
            "relaxation would run to its bounds.  Build the reference with "
            "reference_from_chain(..., valence=SimpleFF.from_preset('pvdf-dft-valence'))."
        )
    free = list(free)
    nf, ns = len(free), shape.n
    p_base = np.asarray(params, dtype=float).reshape(7).copy()
    x0 = shape.x0 if x is None else np.asarray(x, dtype=float).reshape(ns)
    cache: dict = {}

    def build(z):
        z = np.asarray(z, dtype=float)
        key = z.tobytes()
        hit = cache.get(key)
        if hit is not None:
            return hit
        p = p_base.copy()
        p[free] = z[:nf]
        chains = shape.batch(z[nf:])
        packer.update_chain(chains[0])
        E, g_cell, gX, gc = packer.energy_and_grad(p)
        E = float(E) + shape.penalty * chains[0].rotation_error ** 2
        grad = np.empty(nf + ns)
        grad[:nf] = g_cell[free]
        dc = np.empty(ns)
        for j in range(ns):
            plus, minus = chains[1 + 2 * j], chains[2 + 2 * j]
            grad[nf + j] = (_shape_scalar(plus, gX, gc, packer, shape.penalty)
                            - _shape_scalar(minus, gX, gc, packer, shape.penalty)) / (2.0 * shape.step)
            dc[j] = (plus.c - minus.c) / (2.0 * shape.step)
        out = (E, grad, chains[0], dc, p)
        if len(cache) > 96:
            cache.clear()
        cache[key] = out
        return out

    lo, hi = _bounds(ref, p_base, free)
    bounds = [(lo[i], hi[i]) for i in free] + [(-v, v) for v in shape.limits()]
    z0 = np.concatenate([p_base[free], x0])
    if c_target is None:
        res = minimize(lambda z: build(z)[:2], z0, method="L-BFGS-B", jac=True, bounds=bounds,
                       options={"ftol": 1e-12, "gtol": 1e-7, "maxiter": maxiter})
    else:
        cons = [{
            "type": "eq",
            "fun": lambda z: build(z)[2].c / c_target - 1.0,
            "jac": lambda z: np.concatenate([np.zeros(nf), build(z)[3] / c_target]),
        }]
        res = minimize(lambda z: build(z)[0], z0, method="SLSQP", jac=lambda z: build(z)[1],
                       bounds=bounds, constraints=cons, options={"ftol": 1e-12, "maxiter": maxiter})
    E, grad, chain, dc, p = build(res.x)
    packer.update_chain(chain)
    gs = grad[nf:]
    if c_target is None:
        lam, resid = 0.0, float(np.linalg.norm(gs))
    else:
        denom = float(dc @ dc)
        lam = float(gs @ dc) / denom if denom > 0 else 0.0
        scale = max(float(np.linalg.norm(gs)), 1e-12)
        resid = float(np.linalg.norm(gs - lam * dc)) / scale
    return Relaxed(params=p, x=np.asarray(res.x, dtype=float)[nf:], chain=chain, energy=E,
                   multiplier=lam, residual=resid,
                   c_error=0.0 if c_target is None else abs(chain.c / c_target - 1.0),
                   cell_residual=float(np.abs(grad[:nf]).max()), iterations=int(res.nit),
                   shape_gradient=float(np.linalg.norm(gs)))


def deformable_state(ref: Reference, shape: Shape, eps6, x=None, field_lab=None,
                     free=_INTERNAL_FREE, maxiter: int = 120) -> StrainState:
    """A strained state with the chain's conformation relaxed at fixed ``eps_zz``.

    The deformable twin of :func:`_state_at_strain`.  The repeat is held at the strained
    ``c`` by the constraint of :func:`relax_deformable`, so ``eps_zz`` is exactly what the
    strain says it is whether it is zero (an in-plane column) or not (the axial one) --
    the rigid path got that for free and this one has to ask for it.

    ``stress[2]`` is the constraint multiplier converted to GPa, which is the axial stress:
    ``dE*/deps_zz = c0 dE*/dc_target = c0 lambda``.  The other three come from the analytic
    kernel gradient exactly as in the rigid path.
    """
    sc = strained_cell(ref.params, ref.c, eps6)
    saved_field, saved_chain = ref.packer.field, ref.packer.chain
    e_pk = None if field_lab is None else sc.rotate_field(np.asarray(field_lab, dtype=float))
    try:
        if e_pk is not None:
            ref.packer.set_field(e_pk)
        rel = relax_deformable(ref, shape, sc.params, x, free=free, c_target=sc.c, maxiter=maxiter)
        E, sig = stress(ref, sc, params=rel.params, rigid_chain=False)
        sig[AXIAL] = rel.multiplier * ref.c / ref.volume * KCAL_MOL_A3_TO_GPA
        mu = np.asarray(ref.packer.dipole(rel.params[None])[0], dtype=float)
        V = float(ref.packer.cell_volume(rel.params[None])[0])
    finally:
        ref.packer.set_field(saved_field)
        ref.packer.update_chain(saved_chain)
    if e_pk is not None:
        zxE = np.array([-e_pk[1], e_pk[0], 0.0])
        extra = -EV_TO_KCAL * float(mu @ zxE) / ref.volume * KCAL_MOL_A3_TO_GPA
        for K in range(6):
            if K not in UNREACHABLE:
                sig[K] += extra * np.deg2rad(sc.dparams[K, 3])
    m = sc.unrotate(mu) / ref.volume * E_PER_A2_TO_C_PER_M2
    pol = sc.unrotate(mu) / V * E_PER_A2_TO_C_PER_M2
    return StrainState(cell=sc, params=rel.params, energy=E, stress=sig, m=m, polarization=pol,
                       relaxed=rel)


def relax_reference_deformable(ref: Reference, shape: Shape, maxiter: int = 300) -> tuple:
    """Relax the cell **and** the chain against the full energy; return ``(reference, shape)``.

    The reference of an elastic constant has to be a stationary point of the energy that is
    being differentiated twice, and with valence terms in the kernel the refined structure
    is not one: :func:`polyfind.refine.refine_crystal` held the angles with its own invented
    restraint, whose minimum is not the fitted bend terms' minimum.  This relaxes ``c`` too
    -- the only place in this module where ``c`` is free -- so the state that comes back has
    no axial stress *and* no in-plane stress, which is what makes the whole 4x4 block
    meaningful.  The returned :class:`Shape` is rebuilt at the relaxed conformation (so
    ``x = 0`` is the new reference) in the same frame, and the returned :class:`Reference`
    carries the relaxed repeat.
    """
    rel = relax_deformable(ref, shape, ref.params, shape.x0, free=_CELL_FREE, c_target=None,
                           maxiter=maxiter)
    ref2 = Reference(packer=ref.packer, params=rel.params, c=float(rel.chain.c), label=ref.label)
    tors, angs = shape.conformation(rel.x)
    shape2 = shape_of(shape.polymer, ref2, torsions=tors, angles=angs, frame=shape.frame)
    ref.packer.update_chain(rel.chain)
    return ref2, shape2


# --- elastic stiffness ------------------------------------------------------------------
@dataclass
class Elastic:
    """The reachable elastic stiffness, relaxed-ion in the rigid-chain sense.

    ``C`` is 6x6 in GPa with ``nan`` everywhere this parametrisation cannot go; ``block``
    is the 3x3 reachable part in the Voigt order ``(1, 2, 6)``, and ``S = block^-1`` its
    compliance *at clamped* ``eps_zz``.  ``asymmetry`` is the largest ``|C_IJ - C_JI|``
    before symmetrising, a check on the finite difference rather than on physics.

    **The nonbonded cutoff is the accuracy limit, and it has two separate faults.**  At the
    package default the kernel's Lennard-Jones term is energy-shifted but not force-shifted,
    so its force jumps at ``r = rc``
    (:meth:`polyfind.pack.CrystalPacker._pair_energy_and_dv` says so).  Straining the cell
    walks pairs across that jump, and with beta-PVDF's ``b = 8.59`` sitting right on
    ``cutoff=8.0``, ``C_22`` moves by 3.7 GPa between ``step = 5e-4`` and ``step = 4e-3``.
    ``CrystalPacker(lj_cutoff="force")`` removes that: the same spread is 0.43 GPa, and 0.05
    GPa at ``cutoff=12``.  It does **not** remove the other fault, which is truncation:
    beta's ``C_22`` is 37.8 GPa at ``rc=8`` and 42.9 at ``rc=20``, and the force-shifted form
    is if anything slightly further from converged at a given ``rc`` because it subtracts a
    tail as well as a value.  So the force shift fixes the *derivative* and the cutoff length
    fixes the *value*, and they are separate knobs; ``docs/ELECTROMECHANICS.md`` section 6
    measures both and says what adopting either as a default would cost.  The tables there
    are quoted at the package default with the caveat attached rather than silently at a
    different one.
    """

    C: np.ndarray
    block: np.ndarray
    S: np.ndarray
    residual_stress: np.ndarray  # sigma at the reference (GPa), nan where unreachable
    asymmetry: float
    step: float
    reachable: tuple = IN_PLANE
    unreachable: tuple = UNREACHABLE
    c33_from_energy: float = float("nan")  # the curvature route, as a check on C[2, 2]

    def table(self) -> str:
        names = {(0, 0): "C11", (1, 1): "C22", (5, 5): "C66", (0, 1): "C12", (0, 5): "C16", (1, 5): "C26"}
        if AXIAL in self.reachable:
            names.update({(2, 2): "C33", (0, 2): "C13", (1, 2): "C23", (2, 5): "C36"})
        return "  ".join(f"{n}={self.C[i, j]:7.2f}" for (i, j), n in sorted(names.items(), key=lambda kv: kv[1]))


def strain_states(ref: Reference, step: float = 2e-3, relax_internal: bool = True,
                  shape: Shape | None = None, reachable=None) -> dict:
    """``{(K, +1): state, (K, -1): state}`` for every reachable Voigt component.

    Both the stiffness and the piezoelectric constants are built from the same strained,
    relaxed configurations, and on the deformable path each one costs a constrained
    minimisation -- so they are computed once here and shared rather than twice over.
    """
    reachable = (IN_PLANE if shape is None else WITH_AXIAL) if reachable is None else tuple(reachable)
    out = {}
    for K in reachable:
        for sgn in (1, -1):
            s = np.zeros(6)
            s[K] = sgn * step
            out[(K, sgn)] = (_state_at_strain(ref, s, relax_internal) if shape is None
                             else deformable_state(ref, shape, s))
    return out


def elastic_constants(ref: Reference, step: float = 2e-3, relax_internal: bool = True,
                      shape: Shape | None = None, states: dict | None = None) -> Elastic:
    """Second derivatives of the energy density with respect to the reachable strains.

    Central differences of the *stress* (which is analytic, so this is one relaxation per
    displacement rather than one per pair): ``C_JK = dsigma_J/de_K``, two relaxations per
    reachable ``K``.  At the relaxed internal coordinates the partial derivative the kernel
    returns is the total one, so what comes back is the relaxed-ion constant.

    ``shape=None`` (the default) is the rigid-chain calculation: the reachable block is
    :data:`IN_PLANE` and ``C_33`` is ``nan``.  With a :class:`Shape` the block is
    :data:`WITH_AXIAL` -- the chain's conformation relaxes at every strain state, held at
    that state's ``c`` -- and ``C_33``, ``C_13``, ``C_23``, ``C_36`` come out with it.  The
    axial row and column are then measured *twice by different routes*: ``C_J3`` from the
    analytic cell gradient over the axial column, ``C_3J`` from the constraint multiplier
    over the in-plane columns, and :attr:`asymmetry` is their disagreement.  ``C_33``
    additionally gets :attr:`c33_from_energy`, the curvature of the relaxed energy itself.
    """
    reachable = IN_PLANE if shape is None else WITH_AXIAL
    if states is None:
        states = strain_states(ref, step, relax_internal, shape, reachable)
    C = np.full((6, 6), np.nan)
    for K in reachable:
        col = (states[(K, 1)].stress - states[(K, -1)].stress) / (2.0 * step)
        for J in reachable:
            C[J, K] = col[J]
    idx = np.ix_(reachable, reachable)
    blk = C[idx]
    asym = float(np.max(np.abs(blk - blk.T)))
    blk = 0.5 * (blk + blk.T)
    C[idx] = blk
    c33_energy = float("nan")
    _, sig0 = stress(ref)
    if shape is not None:
        # The reference in its own right: the same constrained relaxation at zero strain.
        # Its multiplier is the residual axial stress (which is what says the reference is
        # stress-free) and its energy is the middle point of the axial curvature.
        base = relax_deformable(ref, shape, ref.params, shape.x0, c_target=ref.c)
        sig0[AXIAL] = base.multiplier * ref.c / ref.volume * KCAL_MOL_A3_TO_GPA
        c33_energy = ((states[(AXIAL, 1)].energy + states[(AXIAL, -1)].energy - 2.0 * base.energy)
                      / step ** 2 / ref.volume * KCAL_MOL_A3_TO_GPA)
    return Elastic(C=C, block=blk, S=np.linalg.inv(blk), residual_stress=sig0, asymmetry=asym,
                   step=step, reachable=reachable, c33_from_energy=c33_energy)


# --- piezoelectric response --------------------------------------------------------------
@dataclass
class Piezoelectric:
    """Direct and converse piezoelectric coefficients, and the check that they agree.

    ``e`` (3, n) in C/m^2, rows ``x, y, z`` and columns the reachable Voigt order --
    ``(1, 2, 6)`` on the rigid path, ``(1, 2, 3, 6)`` with a :class:`Shape` -- is the
    **proper** piezoelectric stress constant ``(1/V0) dmu_i/de_J``, the dipole per
    *reference* volume.  ``e_improper`` is the naive ``dP_i/de_J``, the dipole per *current*
    volume, and the two differ by exactly ``P_i`` on the diagonal columns because a
    dilation changes ``P = mu/V`` through ``V`` alone.  Only the proper one is the
    ``-dsigma_J/dE_i`` that appears in the stress, so only the proper one satisfies
    ``d = e S``; keeping both is a cheap check that the difference really is ``P_i``, and a
    useful one, since for a rigid dipole array the naive ``dP/de`` is nothing *but* the
    volume derivative.  With a deformable chain it stops being nothing but that, and the
    difference between the two columns is exactly the argument the PVDF literature has
    about the dimensional effect -- see ``docs/ELECTROMECHANICS.md`` section 8.5.

    ``d_from_e = e S`` and ``d_direct`` (both pC/N) are the two routes to the converse
    coefficient; ``max_abs_difference`` and ``relative_difference`` say how far apart they
    land.  ``field`` is the magnitude (V/A) the converse route used and
    ``converse_iterations`` how many frame iterations it needed; ``d_direct`` is ``nan``
    when the converse route was not run (``converse=False``).

    ``d_improper = e_improper S`` is **not** a piezoelectric constant and no theorem
    connects it to the strain of a field-relaxed cell.  It is reported because on the
    diagonal columns it is exactly the *dimensional* (thickness) term the PVDF literature
    argues about -- ``-P_i`` times the sum of the compliance's diagonal-column entries, the
    charge an electrode sees when a crystal of fixed dipoles changes shape under stress.
    Where the proper ``d`` is zero on those columns, this is what a Broadhurst-Davis
    dimensional model would have said instead, evaluated with this crystal's own
    compliance.  See ``docs/ELECTROMECHANICS.md`` section 8.5.
    """

    e: np.ndarray
    e_improper: np.ndarray
    d_from_e: np.ndarray
    d_direct: np.ndarray
    d_improper: np.ndarray
    max_abs_difference: float
    relative_difference: float
    field: float
    step: float
    converse_iterations: int = 0


def free_strain(ref: Reference, field_lab, elastic: Elastic, iterations: int = 12,
                tol: float = 1e-8, shape: Shape | None = None,
                step_tol: float = 1e-12) -> tuple[np.ndarray, int]:
    """The strain at which the in-plane stress vanishes, with a field held in the reference frame.

    Newton on the analytic stress, whose Jacobian is the stiffness: ``eps -= S sigma``.  The
    converged point satisfies ``sigma = 0`` whatever Jacobian got it there, so using
    ``elastic.S`` here only makes the iteration short -- the answer is still an independent
    measurement of the zero-stress strain and not the ``S e^T`` it is compared against.  The
    internal coordinates are relaxed at every iterate, and with a :class:`Shape` so is the
    chain's conformation.  Returns ``(eps6, iterations used)``.

    ``tol`` should not be set below the reference's own residual stress: past that the
    iteration is chasing the relaxation's noise rather than a root.  :func:`piezoelectric`
    picks it from ``elastic.residual_stress`` for that reason; ``step_tol`` is the second
    stopping rule, on the size of the Newton step.
    """
    eps = np.zeros(6)
    idx = list(elastic.reachable)
    for it in range(1, iterations + 1):
        st = (_state_at_strain(ref, eps, field_lab=field_lab) if shape is None
              else deformable_state(ref, shape, eps, field_lab=field_lab))
        sig = st.stress[idx]
        if float(np.max(np.abs(sig))) < tol:
            return eps, it
        step = elastic.S @ sig
        eps[idx] = eps[idx] - step
        # A stress tolerance alone is not enough on the deformable path: the residual stress
        # of a constrained relaxation has a floor of its own (1e-8 GPa or so), and if ``tol``
        # sits below it the iteration runs to the cap doing nothing.  Stop when the *step*
        # stops mattering instead -- 1e-12 of strain is a coefficient of 5e-9 pC/N at the
        # fields used here, which is below every number this module reports.
        if float(np.max(np.abs(step))) < step_tol:
            return eps, it
    return eps, iterations


def piezoelectric(ref: Reference, elastic: Elastic, step: float = 2e-3, field: float = 0.02,
                  relax_internal: bool = True, converse_iterations: int = 12,
                  shape: Shape | None = None, states: dict | None = None,
                  converse: bool = True) -> Piezoelectric:
    """``e`` at fixed field, ``d = deps/dE`` at zero in-plane stress, and the check between them.

    The direct route strains the cell, relaxes the internal coordinates (and, with a
    :class:`Shape`, the chain's conformation at fixed ``eps_zz``) at each strain, and
    differentiates the cell dipole per reference volume.  The converse route hands each field
    axis in turn to :func:`free_strain`, which drives the analytic stress to zero over the
    reachable components (relaxing the same things at every iterate) and returns the strain
    it lands on.  The
    identity ``d = e S`` follows from ``d^2 h / de dE`` being symmetric; the two numbers come
    from disjoint machinery -- a dipole derivative against a stress root-find -- so agreeing
    is evidence and not bookkeeping.  It caught two sign-and-magnitude errors while this was
    being written: the improper-vs-proper volume term, and the frame the applied field is
    given in under a shear.
    """
    cols = elastic.reachable
    if states is None:
        states = strain_states(ref, step, relax_internal, shape, cols)
    e = np.zeros((3, len(cols)))
    e_improper = np.zeros((3, len(cols)))
    for n, K in enumerate(cols):
        sp, sm = states[(K, 1)], states[(K, -1)]
        e[:, n] = (sp.m - sm.m) / (2.0 * step)
        e_improper[:, n] = (sp.polarization - sm.polarization) / (2.0 * step)

    d_direct = np.full((3, len(cols)), np.nan)
    used = 0
    if converse:
        # The zero-stress root cannot be found more precisely than the reference's own
        # residual stress: below that the Newton is chasing the relaxation's noise, and on
        # the deformable path (where the floor is 1e-8 GPa rather than 1e-16) it would
        # otherwise run to the iteration cap for every field direction.
        floor = float(np.max(np.abs(elastic.residual_stress[list(cols)])))
        tol = max(1e-8, 20.0 * floor)
        d_direct = np.zeros((3, len(cols)))
        for i in range(3):
            strains = []
            for sgn in (1.0, -1.0):
                f = np.zeros(3)
                f[i] = sgn * field
                eps, it = free_strain(ref, f, elastic, iterations=converse_iterations,
                                      tol=tol, shape=shape)
                used = max(used, it)
                strains.append(eps)
            d_direct[i] = (strains[0] - strains[1])[list(cols)] / (2.0 * field * V_PER_A_TO_V_PER_M)
        d_direct *= 1e12  # C/N -> pC/N
    d_from_e = e @ elastic.S * C_PER_M2_PER_GPA_TO_PC_PER_N
    diff = float(np.max(np.abs(d_from_e - d_direct))) if converse else float("nan")
    scale = float(np.max(np.abs(d_from_e)))
    # Below a micro-pC/N there is no coefficient to be relatively wrong about: a crystal with
    # no dipoles lands at 1e-12 pC/N by both routes and a ratio of the two is meaningless.
    return Piezoelectric(e=e, e_improper=e_improper, d_from_e=d_from_e, d_direct=d_direct,
                         d_improper=e_improper @ elastic.S * C_PER_M2_PER_GPA_TO_PC_PER_N,
                         max_abs_difference=diff,
                         relative_difference=(diff / scale if scale > 1e-6 else 0.0) if converse else float("nan"),
                         field=field, step=step, converse_iterations=used)


# --- dielectric response ------------------------------------------------------------------
@dataclass
class Dielectric:
    """The static dielectric tensor of one crystal, at fixed cell, in the reference frame.

    ``clamped`` is the response at fixed geometry -- the induced dipoles alone, which is the
    electronic part ``eps_inf`` a frozen-phonon or DFPT calculation reports -- and is exactly
    1 for any packer without ``polarizable``.  ``relaxed`` lets the internal coordinates
    (and, with a :class:`Shape`, the chain's conformation) follow the field at fixed cell:
    the ionic part this parametrisation can express, which is chain *reorientation* and the
    line group's shape parameters, not the optical phonons of a full lattice dynamics.  Both
    are ``1 + (dP/dE) / eps_0`` from a central difference at ``+-field`` (V/A).
    """

    clamped: np.ndarray
    relaxed: np.ndarray | None
    field: float

    def table(self) -> str:
        def row(name, m):
            return f"{name}: diag({m[0, 0]:.4f}, {m[1, 1]:.4f}, {m[2, 2]:.4f})  max|offdiag| {np.abs(m - np.diag(np.diag(m))).max():.1e}"
        out = [row("eps_inf (clamped)", self.clamped)]
        if self.relaxed is not None:
            out.append(row("eps_0  (relaxed-ion, fixed cell)", self.relaxed))
        return "\n".join(out)


def dielectric_tensor(ref: Reference, field: float = 2e-3, relax: bool = True,
                      shape: Shape | None = None) -> Dielectric:
    r"""``eps_ij = delta_ij + (1/eps_0) dP_i/dE_j`` at fixed cell, by central differences in the field.

    The field is applied in the reference frame (at zero strain the packer's frame is the
    reference frame, so no rotation is needed) and the dipole read back is the packer's total
    one -- static charges plus induced dipoles, which :meth:`polyfind.pack.CrystalPacker.dipole`
    returns because the induced dipoles are stationary and ``mu = -dE/dE_applied`` holds for the
    whole.  Under tinfoil boundaries the applied field *is* the macroscopic field, so this is
    the dielectric constant and not a Clausius-Mossotti local-field quantity; the known-answer
    test for the machinery is the simple cubic lattice in ``tests/test_polarizable.py`` and the
    known answer for the physics is beta-PVDF's ``eps_inf`` (``docs/ELECTROMECHANICS.md``
    section 5.8).

    ``relax=True`` adds the tensor with the internal coordinates relaxed under the field at
    fixed cell (through :func:`_state_at_strain`, or :func:`deformable_state` with a
    :class:`Shape`).  The relaxation has a noise floor, so ``field`` should not be made very
    small for it; the electronic part is linear to machine precision at any small field.
    """
    from .polarizability import E_PER_A2_PER_V_PER_A_TO_CHI

    pk = ref.packer
    V = ref.volume
    c_arr = np.array([ref.c])
    clamped = np.eye(3)
    relaxed = np.eye(3) if relax else None
    saved = pk.field
    try:
        for j in range(3):
            mus, mus_r = [], []
            for sgn in (1.0, -1.0):
                f = np.zeros(3)
                f[j] = sgn * field
                pk.set_field(f)
                mus.append(np.asarray(pk.dipole(ref.params[None], c=c_arr)[0], dtype=float))
                pk.set_field(saved)
                if relax:
                    st = (_state_at_strain(ref, np.zeros(6), True, field_lab=f) if shape is None
                          else deformable_state(ref, shape, np.zeros(6), field_lab=f))
                    mus_r.append(st.m * V / E_PER_A2_TO_C_PER_M2)  # back to e.A
            clamped[:, j] += (mus[0] - mus[1]) / (2.0 * field * V) * E_PER_A2_PER_V_PER_A_TO_CHI
            if relax:
                relaxed[:, j] += (mus_r[0] - mus_r[1]) / (2.0 * field * V) * E_PER_A2_PER_V_PER_A_TO_CHI
    finally:
        pk.set_field(saved)
    return Dielectric(clamped=clamped, relaxed=relaxed, field=field)


# --- actuator figures ---------------------------------------------------------------------
@dataclass
class Actuator:
    """What the crystal does as an actuator at one field.

    ``free_strain = d^T E`` (dimensionless, zero in-plane stress), ``blocking_stress =
    C eps_free = e^T E`` in GPa (the stress the crystal exerts on a rigid clamp; hold it at
    ``eps = 0`` and the external stress is minus this), both in the Voigt order ``(1, 2, 6)``.
    ``work_density`` is ``(1/4) sigma_block . eps_free`` in kJ/m^3, the most work per cycle
    a matched linear load can take; ``elastic_energy`` is the full ``(1/2)`` triangle.
    """

    direction: np.ndarray
    field: float  # V/A
    free_strain: np.ndarray
    blocking_stress: np.ndarray
    work_density: float
    elastic_energy: float

    @property
    def field_v_per_m(self) -> float:
        return self.field * V_PER_A_TO_V_PER_M


def best_direction(elastic: Elastic, piezo: Piezoelectric) -> np.ndarray:
    """The poling direction that gets the most work out, as a unit vector.

    The work at field ``E0 n`` goes as ``n^T (d C d^T) n``, so the answer is the leading
    eigenvector of that symmetric positive-semidefinite 3x3 matrix.  Worth computing rather
    than assuming: the obvious choice, the polarization direction, is the *worst* one for a
    rigid dipole array -- a field along a dipole that is already aligned exerts no torque,
    and with fixed point charges on a rigid chain reorientation is the only piezoelectric
    channel there is.
    """
    d = piezo.d_from_e * 1e-12
    W = d @ elastic.block @ d.T
    w, V = np.linalg.eigh(W)
    return V[:, int(np.argmax(w))]


def actuator(elastic: Elastic, piezo: Piezoelectric, direction, field: float) -> Actuator:
    """Free strain, blocking stress and work density at field ``field`` (V/A) along ``direction``."""
    n = np.asarray(direction, dtype=float).reshape(3)
    norm = float(np.linalg.norm(n))
    n = n / norm if norm > 0 else n
    E = n * field * V_PER_A_TO_V_PER_M
    eps_free = (piezo.d_from_e.T @ E) * 1e-12  # pC/N -> C/N = m/V
    sig_block = elastic.block @ eps_free  # GPa
    w = float(sig_block @ eps_free)
    return Actuator(direction=n, field=field, free_strain=eps_free, blocking_stress=sig_block,
                    work_density=0.25 * w * GPA_TO_KJ_PER_M3, elastic_energy=0.5 * w * GPA_TO_KJ_PER_M3)


# --- axial strain -------------------------------------------------------------------------
@dataclass
class AxialReport:
    """Whether ``eps_zz`` is computable **on the rigid path**, and the evidence either way.

    This is the diagnosis, not the fix.  It measures what a ``C_33`` obtained without valence
    terms in the kernel would be made of, and its verdict (``computable=False``) is about
    that path only; :func:`deformable_state` and :func:`elastic_constants` with a
    :class:`Shape` are the path on which the same quantity *is* computable.

    ``geometrically_reachable`` is whether ``c`` can be moved at all (it can: through the
    backbone bond angles, bond lengths being rigid); ``strain_range`` is how far, within
    the refinement's own ``max_angle_change`` cap.  ``sigma_zz_rigid`` is the residual
    axial stress at a fixed chain shape, which is not a physical stress (see
    :func:`stress`).  ``c33`` and ``sigma_zz`` map a bend stiffness in kcal/mol/rad^2 to the
    ``C_33`` (GPa) and the axial stress (GPa) it produces along the angle path; the entries
    at ``0.0`` are the lattice energy alone, and the entry of ``sigma_zz`` that is nearest
    zero says which bend constant the refined structure is actually in equilibrium under.
    ``computable`` is the verdict this module acts on, and it is ``False``.
    """

    geometrically_reachable: bool
    strain_range: tuple
    angle_range: tuple
    sigma_zz_rigid: float
    c33: dict
    sigma_zz: dict
    computable: bool
    reason: str

    @property
    def bend_fraction(self) -> float:
        """The share of ``C_33`` at ``k = 105`` that the invented bend constant supplies."""
        full = self.c33.get(105.0)
        bare = self.c33.get(0.0)
        return float("nan") if not full else (full - bare) / full


def axial_report(polymer, ref: Reference, torsions=None, angles=None, max_angle_change: float = 8.0,
                 stiffnesses=(0.0, 52.5, 105.0, 210.0), n_points: int = 5,
                 strain: float = 0.01) -> AxialReport:
    """Measure what ``eps_zz`` costs, and where the cost comes from.

    ``c`` is moved by shifting every backbone bond angle of the repeat together (bond
    lengths are rigid, so that is the only handle), the chain is rebuilt, the internal
    coordinates are relaxed at fixed ``(a, b, gamma)``, and a parabola is fitted to the
    energy density.  The ``stiffnesses`` sweep adds
    :func:`polyfind.refine.refine_crystal`'s own harmonic bend restraint at several values
    of its constant, including zero, which makes two things visible at once: how much of
    ``C_33`` that invented constant supplies (79% for beta-PVDF, 84% for PE; less for the
    helical conformations, whose lattice energy resists the angle path more), and -- from
    ``sigma_zz`` -- that the refined structure is in axial equilibrium at ``k = 105`` and at
    no other value, the lattice energy alone having no axial minimum at all.  Hence
    ``computable=False``: this module will not report ``C_33`` as a constant of the potential.

    Two caveats on the path.  It is a *uniform* shift of all backbone angles, which is the
    direction the refinement relaxed along for an all-trans chain but only one direction of a
    larger shape manifold for a helix -- gamma-PVDF's ``sigma_zz`` along it is tens of GPa, so
    its ``C_33`` is partly the stiffness of that constraint.  And the torsions are frozen;
    letting them relax would soften every number.
    """
    from scipy.optimize import brentq

    from .linegroup import repeat_chains

    chain = ref.packer.chain
    B = polymer.bonds_per_repeat
    tors = np.asarray(chain.dihedrals if torsions is None else torsions, dtype=float)
    ang0 = (np.array([polymer.backbone[k].backbone_angle for k in range(B)], dtype=float)
            if angles is None else np.asarray(angles, dtype=float))
    c0, V0 = ref.c, ref.volume
    n_angles = len(tors)  # backbone angles per chain repeat, as refine.bend_energy counts them

    def c_of(delta: float) -> float:
        return float(repeat_chains(polymer, chain.name, tors[None], ang0 + delta)[0].c)

    lo, hi = -max_angle_change, max_angle_change
    c_lo, c_hi = c_of(lo), c_of(hi)
    reach = (c_lo / c0 - 1.0, c_hi / c0 - 1.0)
    target = min(strain, 0.5 * min(abs(reach[0]), abs(reach[1])))
    if not c_lo < c0 < c_hi:
        raise ValueError("the reference repeat is outside the range the bond angles can reach; "
                         "pass the refined backbone angles as `angles`")

    saved_chain = chain
    rows = []
    try:
        for eps in np.linspace(-target, target, n_points):
            ct = c0 * (1.0 + eps)
            # delta is measured from the *polymer's* equilibrium angles, so that the bend
            # restraint below is the one refine_crystal applies, but it is solved for the
            # target c rather than assumed zero at the reference: the reference chain has
            # already been refined and its angles are not the polymer's.
            delta = brentq(lambda d: c_of(d) - ct, lo, hi, xtol=1e-12)
            ch = repeat_chains(polymer, chain.name, tors[None], ang0 + delta, align_to=saved_chain)[0]
            ref.packer.update_chain(ch)
            p = ref.params.copy()
            p[5] = ref.params[5] * ch.c / c0
            p = relax(Reference(ref.packer, p, ch.c, ref.label), p)
            E = float(ref.packer.energy(p[None])[0])
            bend = 0.5 * np.deg2rad(delta) ** 2 * n_angles * ref.packer.n_chains
            rows.append((eps, E, bend))
    finally:
        ref.packer.update_chain(saved_chain)

    eps_v = np.array([r[0] for r in rows])
    E_v = np.array([r[1] for r in rows])
    bend_v = np.array([r[2] for r in rows])
    c33, sig_zz = {}, {}
    for k in stiffnesses:
        coef = np.polyfit(eps_v, (E_v + k * bend_v) / V0, 2)
        c33[float(k)] = float(2.0 * coef[0] * KCAL_MOL_A3_TO_GPA)
        sig_zz[float(k)] = float(coef[1] * KCAL_MOL_A3_TO_GPA)
    _, sig0 = stress(ref)
    share = ""
    if c33.get(105.0):
        share = f" (which supplies {(c33[105.0] - c33[0.0]) / c33[105.0] * 100:.0f}% of it at k = 105)"
    return AxialReport(
        geometrically_reachable=True,
        strain_range=reach,
        angle_range=(lo, hi),
        sigma_zz_rigid=float(sig0[AXIAL]),
        c33=c33,
        sigma_zz=sig_zz,
        computable=False,
        reason=(
            "eps_zz can be applied -- opening the backbone bond angles moves c over "
            f"{reach[0] * 100:+.1f}% to {reach[1] * 100:+.1f}% within the refinement's "
            f"+/-{max_angle_change:g} deg cap -- but on this path its energy is not the potential's.  "
            "A packer without valence terms has no bond or angle term at all (and "
            "FFParameters.applied does not forward the fitted ones), so the only thing resisting is "
            "refine_crystal's own harmonic bend restraint, an invented constant that is excluded from "
            f"the reported lattice energy.  C_33 is exactly affine in that constant{share}, and the "
            "entry of sigma_zz nearest zero shows the refined structure is in axial equilibrium under "
            "it and under nothing else.  Reported as nan.  Give the packer the fitted valence terms "
            "(CrystalPacker(valence=...)) and pass a Shape to elastic_constants and C_33 becomes a "
            "constant of the potential instead -- bond stretching still absent, because bond lengths "
            "are rigid, which makes even that an upper bound."
        ),
    )


# --- the whole response ---------------------------------------------------------------------
@dataclass
class Response:
    """Everything :func:`electromechanical_response` computes for one crystal."""

    label: str
    reference: Reference
    elastic: Elastic
    piezo: Piezoelectric
    actuator: Actuator
    polarization: np.ndarray
    volume: float
    seconds: float
    axial: AxialReport | None = None
    shape: Shape | None = None
    notes: list = dc_field(default_factory=list)

    def summary(self) -> str:
        el, pz, ac = self.elastic, self.piezo, self.actuator
        d = pz.d_from_e
        deformable = AXIAL in el.reachable
        cols = " ".join(f"eps_{VOIGT_LABELS[K]}" for K in el.reachable)
        lines = [
            f"{self.label}: V0={self.volume:.2f} A^3  |P|={np.linalg.norm(self.polarization):.4f} C/m^2  "
            f"({self.seconds:.1f} s)",
            f"  elastic (GPa, relaxed-ion, "
            f"{'chain deforms at fixed eps_zz' if deformable else 'rigid chain, eps_zz clamped'}): {el.table()}",
            (f"  C_33 also = {el.c33_from_energy:.2f} GPa by the energy-curvature route "
             f"(stress route {el.C[2, 2]:.2f}); axial asymmetry check in Elastic.asymmetry"
             if deformable else
             "  C_33 = nan (axial strain: see axial_report); C_44, C_55, C_45 and every mixed "
             "constant with 4 or 5: not expressible"),
            f"  e (C/m^2), rows x y z, cols {cols}:",
        ]
        for i, r in enumerate("xyz"):
            lines.append(f"    {r}  " + "  ".join(f"{v:+9.4f}" for v in pz.e[i]))
        lines.append("  d (pC/N) from e S | direct:")
        for i, r in enumerate("xyz"):
            lines.append(f"    {r}  " + "  ".join(f"{v:+9.3f}" for v in d[i])
                         + "   |   " + "  ".join(f"{v:+9.3f}" for v in pz.d_direct[i]))
        lines.append(f"  direct/converse agreement: max |dd| = {pz.max_abs_difference:.3g} pC/N "
                     f"({pz.relative_difference * 100:.2f}% of the largest coefficient)")
        lines.append(f"  actuator at E = {ac.field:g} V/A along ({', '.join(f'{v:+.3f}' for v in ac.direction)}):")
        lines.append("    free strain  " + "  ".join(f"{v:+.3e}" for v in ac.free_strain))
        lines.append("    blocking     " + "  ".join(f"{v:+.4f}" for v in ac.blocking_stress) + " GPa")
        lines.append(f"    work density {ac.work_density:.4g} kJ/m^3 (matched load; "
                     f"{ac.elastic_energy:.4g} kJ/m^3 for the full triangle)")
        return "\n".join(lines)


def electromechanical_response(ref: Reference, step: float = 2e-3, field: float = 0.02,
                               actuator_field: float = 0.01, direction=None,
                               relax_first: bool = True, polymer=None,
                               axial: bool = False, shape: Shape | None = None,
                               converse: bool = True) -> Response:
    """Elastic constants, piezoelectric coefficients and actuator figures for one crystal.

    ``actuator_field`` is in V/A; the default 0.01 V/A is 100 MV/m, a realistic drive field
    for a poled PVDF film and about twice its coercive field, and the figures it produces are
    a *linearised* response evaluated there.  ``direction`` is the poling direction; ``None``
    takes :func:`best_direction`.  ``axial=True`` additionally runs :func:`axial_report`,
    which needs ``polymer`` -- the *rigid* diagnosis of why ``eps_zz`` is not computable
    without valence terms.

    ``shape`` (a :class:`Shape`, which needs a packer built with ``valence=``) switches the
    whole calculation onto the deformable path: the reference is relaxed over ``c`` as well
    as the cell, the reachable block becomes :data:`WITH_AXIAL`, and the diagonal columns of
    ``e`` are no longer identically zero.  That is the one change that makes ``C_33`` and
    the ``d_31``/``d_33`` family reachable at all; ``shape=None`` leaves every number
    exactly what it was.
    """
    t0 = time.time()
    if relax_first:
        if shape is None:
            ref = relax_reference(ref)
        else:
            ref, shape = relax_reference_deformable(ref, shape)
    states = strain_states(ref, step, True, shape)
    el = elastic_constants(ref, step=step, shape=shape, states=states)
    pz = piezoelectric(ref, el, step=step, field=field, shape=shape, states=states, converse=converse)
    pol = ref.polarization()
    if direction is None:
        direction = best_direction(el, pz)
        if np.linalg.norm(direction) < 1e-12:
            direction = np.array([1.0, 0.0, 0.0])
    ac = actuator(el, pz, direction, actuator_field)
    ax = axial_report(polymer, ref) if axial else None
    rigid_notes = [
        "relaxed-ion means the rigid-body internal coordinates (phi1, phi2, dz); the chain "
        "conformation is held rigid",
        "C_33 is nan: eps_zz is reachable only through the bond angles, whose only restoring term "
        "anywhere in the pipeline is refine_crystal's invented bend constant (axial_report)",
        "d is the axially clamped strain coefficient: S is the inverse of the in-plane block only",
        "the charges are fixed point charges on a rigid chain, so chain *reorientation* is the only "
        "piezoelectric channel in the model: a dilation cannot change the dipole at all, and the "
        "proper e therefore vanishes on every diagonal column",
    ]
    deformable_notes = [
        "relaxed-ion means the rigid-body internal coordinates (phi1, phi2, dz) *and* the chain's "
        "line-group shape parameters, relaxed at every strain state against the valence terms",
        "eps_zz is a constraint on the repeat, not an identity: C_33 is the constraint multiplier's "
        "derivative and is checked against the curvature of the relaxed energy",
        "d is still not the free coefficient of a full 6x6: eps_yz and eps_xz remain inexpressible, "
        "so S is the inverse of the 4x4 (1, 2, 3, 6) block",
        "bond lengths are still rigid, so the stretch channel -- roughly half of a real chain "
        "modulus -- contributes a constant and nothing else",
    ]
    notes = (rigid_notes if shape is None else deformable_notes) + [
        "eps_yz and eps_xz are not expressible in this cell, so C_44, C_55, C_45 and every mixed "
        "constant involving them are absent, not merely unconverged",
        "not a prediction of experiment; see docs/ELECTROMECHANICS.md",
    ]
    return Response(label=ref.label, reference=ref, elastic=el, piezo=pz, actuator=ac, polarization=pol,
                    volume=ref.volume, seconds=time.time() - t0, axial=ax, shape=shape, notes=notes)
