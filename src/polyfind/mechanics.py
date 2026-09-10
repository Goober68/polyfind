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
  of the chain's repeat transform, fixed by the chain's geometry.  See
  :func:`axial_report`.

The reachable stiffness is therefore the 3x3 block ``{11, 22, 12, 16, 26, 66}``.  Those
are the *true* constants and not a plane-stress reduction: ``C_11 = dsigma_1/deps_1`` is
defined at fixed ``eps_2 ... eps_6``, and fixed ``eps_3`` is exactly what a rigid chain
enforces.  What cannot be had from them is the *compliance*: ``S = C^-1`` of the 3x3
block is the compliance at clamped ``eps_zz``, so every ``d`` below is an axially clamped
piezoelectric strain coefficient.

**Axial strain: reachable geometrically, not computable energetically.**  ``c`` can be
changed -- bond *lengths* are rigid but backbone bond *angles* are not, and opening every
backbone angle of PE by 8 degrees lengthens the repeat by 4.5% (PVDF alpha and gamma by
7.3%).  What is missing is the energy.  The packing kernel is Lennard-Jones +
damped-shifted-force Coulomb + a Fourier torsion term and has **no bond or angle term**;
:meth:`polyfind.fitting.FFParameters.applied` deliberately does not forward the fitted
valence terms into packing, on the grounds that packing holds the chain rigid.  The only
thing anywhere in the pipeline that resists opening an angle is
:func:`polyfind.refine.refine_crystal`'s own harmonic restraint, whose stiffness
(``angle_stiffness = 105`` kcal/mol/rad^2, "UFF-like for sp3 carbon") is an invented
number that is deliberately excluded from the reported lattice energy.  A ``C_33``
obtained along that path is a report of that constant; and it is missing the bond-stretch
contribution that dominates a real chain modulus, because bonds cannot stretch at all.
Worse, a ``dE/dc`` taken at *fixed chain shape* -- the other way the cell could change
``c`` -- is not an axial strain in any physical sense: it slides the repeats apart across
the periodic boundary, and the nonbonded sum excludes exactly the bonded pairs that cross
it, so nothing resists.  :func:`axial_report` measures all three of these.  Measured (see
``docs/ELECTROMECHANICS.md`` section 4): ``C_33`` is exactly affine in that invented
constant, and 79% of beta's and 84% of PE's comes from it; the refined structure is
stationary in ``c`` at ``k = 105`` and nowhere else, the lattice energy on its own wanting
to contract the chain at 2 to 5 GPa with no minimum at all.
:func:`electromechanical_response` therefore leaves ``C_33`` as ``nan``.

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
over strain states affordable.  The chain conformation is held rigid, so these are
*rigid-chain* relaxed-ion constants: softer than clamped-ion, stiffer than constants
whose torsions also relax.  Letting the conformation move is not offered for the in-plane
block on purpose -- it would change ``c``, and the calculation would no longer be at
fixed ``eps_zz``.

Nothing here is a prediction of experiment.  The illustrative potential is not quantitative
(DESIGN.md section 5.4) and the fitted presets are fitted to single-chain DFT energies and
forces, with no elastic or piezoelectric data in the objective at all.
``docs/ELECTROMECHANICS.md`` puts published PVDF and polyethylene numbers beside these, says
what could not be found to compare against, and records one further limit worth knowing
before reading any constant here: with fixed point charges on a rigid chain the *only*
piezoelectric channel in the model is chain reorientation, so every diagonal column of ``e``
is zero and neither of the two mechanisms the PVDF literature argues over is present.
``examples/electromechanics.py`` runs the whole thing on the reference polymorphs.
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
AXIAL = 2  # Voigt 3: needs c, the chain's own repeat (see axial_report)
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


def stress(ref: Reference, sc: StrainedCell | None = None, params=None) -> tuple[float, np.ndarray]:
    """``(energy, sigma)`` for a strained cell: ``sigma_J = (1/V0) dE/de_J`` in GPa.

    The kernel's analytic ``dE/d(a, b, gamma, phi1, phi2, dz)`` and ``dE/dc`` are
    contracted with :attr:`StrainedCell.dparams` and :attr:`StrainedCell.dc`, so no extra
    energy evaluation is spent and the derivative is closed-form all the way from the
    strain.  ``sigma[3]``, ``sigma[4]`` are ``nan``: :data:`UNREACHABLE`.

    ``sigma[2]`` is returned but is **not** an axial stress in any useful sense: at a fixed
    chain shape the only way ``c`` can grow is to slide the repeats apart across the
    periodic boundary, and the bonded exclusions remove exactly the pairs that cross it.
    :func:`axial_report` quotes it as a diagnostic and nothing else uses it.
    """
    if sc is None:
        sc = strained_cell(ref.params, ref.c, np.zeros(6))
    p = sc.params if params is None else np.asarray(params, dtype=float).reshape(7)
    c_arr = None if abs(sc.c - ref.packer.chain.c) < 1e-12 else np.array([sc.c])
    E, g_cell, _, g_c = ref.packer.energy_and_grad(p, coords=None if c_arr is None else ref.packer.chain.coords, c=c_arr)
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


# --- elastic stiffness ------------------------------------------------------------------
@dataclass
class Elastic:
    """The reachable elastic stiffness, relaxed-ion in the rigid-chain sense.

    ``C`` is 6x6 in GPa with ``nan`` everywhere this parametrisation cannot go; ``block``
    is the 3x3 reachable part in the Voigt order ``(1, 2, 6)``, and ``S = block^-1`` its
    compliance *at clamped* ``eps_zz``.  ``asymmetry`` is the largest ``|C_IJ - C_JI|``
    before symmetrising, a check on the finite difference rather than on physics.

    **The nonbonded cutoff is the accuracy limit, not the finite difference.**  The kernel's
    Lennard-Jones term is energy-shifted but not force-shifted, so its force jumps at
    ``r = rc`` (:meth:`polyfind.pack.CrystalPacker._pair_energy_and_dv` says so).  Straining
    the cell walks pairs across that jump, and at the package default ``cutoff=8.0`` -- with
    beta-PVDF's ``b = 8.59`` sitting right on it -- ``C_22`` still moves by 3.7 GPa between
    ``step = 5e-4`` and ``step = 4e-3``.  At ``cutoff=12`` the step dependence collapses to
    0.06 GPa and ``C_22`` settles 11% higher; at 16 it moves another 1.5%.  Elastic constants
    from this module therefore want a longer cutoff than a lattice-energy ranking does, and
    the ones in ``docs/ELECTROMECHANICS.md`` are quoted at the package default with that
    caveat attached rather than silently at a different one.
    """

    C: np.ndarray
    block: np.ndarray
    S: np.ndarray
    residual_stress: np.ndarray  # sigma at the reference (GPa), nan where unreachable
    asymmetry: float
    step: float
    reachable: tuple = IN_PLANE
    unreachable: tuple = UNREACHABLE

    def table(self) -> str:
        names = {(0, 0): "C11", (1, 1): "C22", (5, 5): "C66", (0, 1): "C12", (0, 5): "C16", (1, 5): "C26"}
        return "  ".join(f"{n}={self.C[i, j]:7.2f}" for (i, j), n in names.items())


def elastic_constants(ref: Reference, step: float = 2e-3, relax_internal: bool = True) -> Elastic:
    """Second derivatives of the energy density with respect to the reachable strains.

    Central differences of the *stress* (which is analytic, so this is one relaxation per
    displacement rather than one per pair): ``C_JK = dsigma_J/de_K``, two relaxations per
    ``K`` in :data:`IN_PLANE`, six in total for the whole block.  At the relaxed internal
    coordinates the partial derivative the kernel returns is the total one, so what comes
    back is the relaxed-ion constant.
    """
    C = np.full((6, 6), np.nan)
    cols = {}
    for K in IN_PLANE:
        s = np.zeros(6)
        s[K] = step
        cols[K] = (_state_at_strain(ref, s, relax_internal).stress
                   - _state_at_strain(ref, -s, relax_internal).stress) / (2.0 * step)
    for K in IN_PLANE:
        for J in IN_PLANE:
            C[J, K] = cols[K][J]
    idx = np.ix_(IN_PLANE, IN_PLANE)
    blk = C[idx]
    asym = float(np.max(np.abs(blk - blk.T)))
    blk = 0.5 * (blk + blk.T)
    C[idx] = blk
    _, sig0 = stress(ref)
    return Elastic(C=C, block=blk, S=np.linalg.inv(blk), residual_stress=sig0, asymmetry=asym, step=step)


# --- piezoelectric response --------------------------------------------------------------
@dataclass
class Piezoelectric:
    """Direct and converse piezoelectric coefficients, and the check that they agree.

    ``e`` (3, 3) in C/m^2, rows ``x, y, z`` and columns the Voigt order ``(1, 2, 6)``, is
    the **proper** piezoelectric stress constant ``(1/V0) dmu_i/de_J`` -- the dipole per
    *reference* volume.  ``e_improper`` is the naive ``dP_i/de_J``, the dipole per *current*
    volume, and the two differ by exactly ``P_i`` on the diagonal columns because a
    dilation changes ``P = mu/V`` through ``V`` alone.  Only the proper one is the
    ``-dsigma_J/dE_i`` that appears in the stress, so only the proper one satisfies
    ``d = e S``; keeping both is a cheap check that the difference really is ``P_i``, and a
    useful one, since for a rigid dipole array the naive ``dP/de`` is nothing but the
    volume derivative.

    ``d_from_e = e S`` and ``d_direct`` (both pC/N) are the two routes to the converse
    coefficient; ``max_abs_difference`` and ``relative_difference`` say how far apart they
    land.  ``field`` is the magnitude (V/A) the converse route used and
    ``converse_iterations`` how many frame iterations it needed.
    """

    e: np.ndarray
    e_improper: np.ndarray
    d_from_e: np.ndarray
    d_direct: np.ndarray
    max_abs_difference: float
    relative_difference: float
    field: float
    step: float
    converse_iterations: int = 0


def free_strain(ref: Reference, field_lab, elastic: Elastic, iterations: int = 12,
                tol: float = 1e-8) -> tuple[np.ndarray, int]:
    """The strain at which the in-plane stress vanishes, with a field held in the reference frame.

    Newton on the analytic stress, whose Jacobian is the stiffness: ``eps -= S sigma``.  The
    converged point satisfies ``sigma = 0`` whatever Jacobian got it there, so using
    ``elastic.S`` here only makes the iteration short -- the answer is still an independent
    measurement of the zero-stress strain and not the ``S e^T`` it is compared against.  The
    internal coordinates are relaxed at every iterate.  Returns ``(eps6, iterations used)``.
    """
    eps = np.zeros(6)
    idx = list(IN_PLANE)
    for it in range(1, iterations + 1):
        sig = _state_at_strain(ref, eps, field_lab=field_lab).stress[idx]
        if float(np.max(np.abs(sig))) < tol:
            return eps, it
        eps[idx] = eps[idx] - elastic.S @ sig
    return eps, iterations


def piezoelectric(ref: Reference, elastic: Elastic, step: float = 2e-3, field: float = 0.02,
                  relax_internal: bool = True, converse_iterations: int = 12) -> Piezoelectric:
    """``e`` at fixed field, ``d = deps/dE`` at zero in-plane stress, and the check between them.

    The direct route strains the cell, relaxes the internal coordinates at each strain, and
    differentiates the cell dipole per reference volume.  The converse route hands each field
    axis in turn to :func:`free_strain`, which drives the analytic in-plane stress to zero
    (internal coordinates relaxed at every iterate) and returns the strain it lands on.  The
    identity ``d = e S`` follows from ``d^2 h / de dE`` being symmetric; the two numbers come
    from disjoint machinery -- a dipole derivative against a stress root-find -- so agreeing
    is evidence and not bookkeeping.  It caught two sign-and-magnitude errors while this was
    being written: the improper-vs-proper volume term, and the frame the applied field is
    given in under a shear.
    """
    e = np.zeros((3, 3))
    e_improper = np.zeros((3, 3))
    for n, K in enumerate(IN_PLANE):
        s = np.zeros(6)
        s[K] = step
        sp, sm = _state_at_strain(ref, s, relax_internal), _state_at_strain(ref, -s, relax_internal)
        e[:, n] = (sp.m - sm.m) / (2.0 * step)
        e_improper[:, n] = (sp.polarization - sm.polarization) / (2.0 * step)

    d_direct = np.zeros((3, 3))
    used = 0
    for i in range(3):
        strains = []
        for sgn in (1.0, -1.0):
            f = np.zeros(3)
            f[i] = sgn * field
            eps, it = free_strain(ref, f, elastic, iterations=converse_iterations)
            used = max(used, it)
            strains.append(eps)
        d_direct[i] = (strains[0] - strains[1])[list(IN_PLANE)] / (2.0 * field * V_PER_A_TO_V_PER_M)
    d_direct *= 1e12  # C/N -> pC/N
    d_from_e = e @ elastic.S * C_PER_M2_PER_GPA_TO_PC_PER_N
    diff = float(np.max(np.abs(d_from_e - d_direct)))
    scale = float(np.max(np.abs(d_from_e)))
    # Below a micro-pC/N there is no coefficient to be relatively wrong about: a crystal with
    # no dipoles lands at 1e-12 pC/N by both routes and a ratio of the two is meaningless.
    return Piezoelectric(e=e, e_improper=e_improper, d_from_e=d_from_e, d_direct=d_direct,
                         max_abs_difference=diff, relative_difference=diff / scale if scale > 1e-6 else 0.0,
                         field=field, step=step, converse_iterations=used)


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
    """Whether ``eps_zz`` is computable here, and the evidence either way.

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
            f"+/-{max_angle_change:g} deg cap -- but its energy is not the potential's.  The packing "
            "kernel has no bond or angle term (and FFParameters.applied deliberately does not "
            "forward the fitted ones), so the only thing resisting is refine_crystal's own harmonic "
            "bend restraint, an invented constant that is excluded from the reported lattice energy.  "
            f"C_33 is exactly affine in that constant{share}, and the entry of sigma_zz nearest zero "
            "shows the refined structure is in axial equilibrium under it and under nothing else.  "
            "Bond stretching, which dominates a real chain modulus, is absent entirely because bond "
            "lengths are rigid.  Reported as nan."
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
    notes: list = dc_field(default_factory=list)

    def summary(self) -> str:
        el, pz, ac = self.elastic, self.piezo, self.actuator
        d = pz.d_from_e
        lines = [
            f"{self.label}: V0={self.volume:.2f} A^3  |P|={np.linalg.norm(self.polarization):.4f} C/m^2  "
            f"({self.seconds:.1f} s)",
            f"  elastic (GPa, relaxed-ion, rigid chain, eps_zz clamped): {el.table()}",
            "  C_33 = nan (axial strain: see axial_report); C_44, C_55, C_45 and every mixed "
            "constant with 4 or 5: not expressible",
            "  e (C/m^2), rows x y z, cols eps_xx eps_yy eps_xy:",
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
                               axial: bool = False) -> Response:
    """Elastic constants, piezoelectric coefficients and actuator figures for one crystal.

    ``actuator_field`` is in V/A; the default 0.01 V/A is 100 MV/m, a realistic drive field
    for a poled PVDF film and about twice its coercive field, and the figures it produces are
    a *linearised* response evaluated there.  ``direction`` is the poling direction; ``None``
    takes :func:`best_direction`.  ``axial=True`` additionally runs :func:`axial_report`,
    which needs ``polymer``.
    """
    t0 = time.time()
    if relax_first:
        ref = relax_reference(ref)
    el = elastic_constants(ref, step=step)
    pz = piezoelectric(ref, el, step=step, field=field)
    pol = ref.polarization()
    if direction is None:
        direction = best_direction(el, pz)
        if np.linalg.norm(direction) < 1e-12:
            direction = np.array([1.0, 0.0, 0.0])
    ac = actuator(el, pz, direction, actuator_field)
    ax = axial_report(polymer, ref) if axial else None
    notes = [
        "relaxed-ion means the rigid-body internal coordinates (phi1, phi2, dz); the chain "
        "conformation is held rigid",
        "eps_yz and eps_xz are not expressible in this cell, so C_44, C_55, C_45 and every mixed "
        "constant involving them are absent, not merely unconverged",
        "C_33 is nan: eps_zz is reachable only through the bond angles, whose only restoring term "
        "anywhere in the pipeline is refine_crystal's invented bend constant (axial_report)",
        "d is the axially clamped strain coefficient: S is the inverse of the in-plane block only",
        "the charges are fixed point charges on a rigid chain, so chain *reorientation* is the only "
        "piezoelectric channel in the model: a dilation cannot change the dipole at all, and the "
        "proper e therefore vanishes on every diagonal column",
        "not a prediction of experiment; see docs/ELECTROMECHANICS.md",
    ]
    return Response(label=ref.label, reference=ref, elastic=el, piezo=pz, actuator=ac, polarization=pol,
                    volume=ref.volume, seconds=time.time() - t0, axial=ax, notes=notes)
