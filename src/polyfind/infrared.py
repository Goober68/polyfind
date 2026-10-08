r"""Infrared intensities of Gamma-point phonons, from the model's own dipole derivative.

A Gamma-point mode ``k`` with mass-weighted eigenvector ``e_k`` (the columns of
:attr:`polyfind.phonon.Phonons.modes`) carries the cell dipole derivative

    ``d_k,a = d mu_a / d Q_k = sum_{i, b} Z*_{i, ab} e_{k, i b} / sqrt(m_i)``     (e / sqrt(amu)),

and its infrared intensity is proportional to ``|d_k|^2``; the direction of ``d_k`` is the
polarization of the light that excites it.  Everything here is that contraction and what it
needs, done on the same Hamiltonian the Hessian comes from:

* **The Born tensors** (:func:`born_charges_by_field`).  ``Z*_{i, ab} = d mu_a / d u_{i b}`` is
  the mixed second derivative of the cell energy with respect to a uniform applied field ``F``
  (V/A) and the displacement of atom ``i``.  The packer's field term is ``-mu . F`` in eV,
  so ``Z*_{i, ab} = -(1 / EV_TO_KCAL) d g_{i b} / d F_a``, with ``g`` the all-atom gradient of
  :func:`polyfind.phonon.cell_energy_and_grad` -- the very function the Hessian is a finite
  difference of.  The field enters that energy twice, as ``-(sum_j q_j r_j) . F`` with the
  fluxing charges and through the induced dipoles' permanent field ``E0``, and at the
  induced-dipole stationary point ``-dE/dF`` is the whole dipole ``sum_j q_j r_j + sum_j p_j``;
  so the static charge, the charge flux and the electronic screening are all in ``Z*``
  without a separate dipole function to keep consistent.  At fixed geometry the energy is
  exactly quadratic in ``F`` (the charges do not depend on the field; the induced dipoles are
  linear in it), so the central difference in ``F`` carries no truncation error, only
  rounding, and ``EV_TO_KCAL`` cancels between the field term and this division.
  :func:`born_charges_by_dipole` is the independent route -- central differences of any
  supplied cell-dipole function in the placed coordinates -- for the check.
* **The acoustic sum rule** ``sum_i Z*_i = 0``: a rigid translation moves no charge and
  induces no dipole in a neutral cell, so it holds identically and is a check on the
  arithmetic; when it holds, the three uniform translations carry no intensity.
* **The symmetry** (:func:`space_group`): every operation ``(R, t)`` of the placed cell,
  found by matching atoms of the same element modulo the lattice, and the character
  ``e_k . (g e_k)`` of every mode under each.  In a centrosymmetric cell the modes even under
  inversion (``g``) carry no dipole -- they are the Raman-active set -- and the odd ones
  (``u``) are the infrared-active set, Raman-inactive: the rule of mutual exclusion.  A
  computed intensity that does not vanish on a ``g`` mode is a defect in ``Z*`` or in the
  eigenvectors, and :attr:`InfraredMode.parity` lets a caller see it.  For the point group
  ``C2h`` (``P2_1/c``, the crystallographic alpha-PVDF cell) the modes are labelled
  ``Ag``/``Bg``/``Au``/``Bu``; ``Au`` dipoles lie along the two-fold axis and ``Bu`` dipoles
  perpendicular to it.  A cell without an inversion has no mutual exclusion (in ``C2v`` every
  mode is Raman-active and all but ``A2`` infrared-active), and nothing here assumes the cell
  is the group the crystallography says: the group is found from the coordinates.  Where every
  operation is along the cell's axes, :func:`allowed_axes` gives each mode's symmetry-allowed
  dipole directions, the selection rule the computed ``d_k`` must obey.

**Units.**  ``|d_k|^2`` is in ``e^2/amu``.  The integrated absorption of a mode, per mole of
cells, is ``A_k = N_A / (12 eps0 c^2) |d mu / d Q_k|^2`` (SI; the usual harmonic
double-harmonic expression), which is :data:`KM_PER_MOL_PER_E2_AMU` = 974.9 km/mol per
``e^2/amu``, derived below from CODATA 2018 values.  The absolute scale is the model's and is
given for completeness; a comparison with a measured spectrum should use the ratios and the
directions.

**What is not in it.**  The intensity is the harmonic, linear-in-``Q`` one (no electrical or
mechanical anharmonicity, no overtones or combinations, no temperature); the frequencies are
the transverse ones of the Gamma-point Hessian with tinfoil (zero macroscopic field)
electrostatics, i.e. no LO-TO splitting.  :func:`lo_shift` gives the first-order size of the
missing non-analytic term for one mode, with the model's own clamped-ion dielectric tensor
from :func:`electronic_dielectric`, so that the omission has a number rather than a word.
"""
from __future__ import annotations

import itertools
from dataclasses import dataclass, field

import numpy as np

from .forcefield import COULOMB
from .pack import EV_TO_KCAL

# --- the absolute intensity scale, from the stated constants ---------------------------------
_E_C = 1.602176634e-19  # C, exact (SI 2019)
_N_A = 6.02214076e23  # 1/mol, exact
_EPS0 = 8.8541878128e-12  # F/m, CODATA 2018
_C_M_PER_S = 299792458.0  # exact
_AMU_KG = 1.66053906660e-27  # CODATA 2018
#: integrated absorption ``N_A / (12 eps0 c^2) |dmu/dQ|^2`` in km/mol per (e^2/amu); 974.9
KM_PER_MOL_PER_E2_AMU = _N_A / (12.0 * _EPS0 * _C_M_PER_S ** 2) * _E_C ** 2 / _AMU_KG / 1e3
#: ``1 / eps0`` in (V/A) per (e/A^2): a polarization of 1 e/A^2 is a field of 180.95 V/A
INV_EPS0_V_A = _E_C / (_EPS0 * 1e-10)


# ------------------------------------------------------------------------- the Born tensors
def born_charges_by_field(packer, params, Pn, latn, step: float = 1e-3) -> np.ndarray:
    """``Z*`` (N, 3, 3) in e at placed coordinates ``Pn``: ``Z[i, a, b] = d mu_a / d u_{i b}``.

    Central differences of :func:`polyfind.phonon.cell_energy_and_grad`'s all-atom gradient
    in a uniform applied field of ``+-step`` V/A along each axis, divided by
    ``-EV_TO_KCAL`` (see the module docstring).  The packer's own field is restored on exit,
    whatever happens.  Six gradient evaluations in all.
    """
    from .phonon import cell_energy_and_grad

    params = np.asarray(params, dtype=float).reshape(7)
    Pn = np.asarray(Pn, dtype=float).reshape(packer.N, 3)
    latn = np.asarray(latn, dtype=float).reshape(3, 3)
    old = getattr(packer, "field", None)
    old = None if old is None else np.array(old, dtype=float)
    Z = np.zeros((packer.N, 3, 3))
    try:
        for a in range(3):
            F = np.zeros(3)
            F[a] = step
            packer.set_field(F)
            gp = cell_energy_and_grad(packer, params, Pn, latn)[1]
            packer.set_field(-F)
            gm = cell_energy_and_grad(packer, params, Pn, latn)[1]
            Z[:, a, :] = -(gp - gm) / (2.0 * step * EV_TO_KCAL)
    finally:
        packer.set_field(old)
    return Z


def born_charges_by_dipole(dipole, Pn, h: float = 1e-4) -> np.ndarray:
    """``Z*`` (N, 3, 3) in e by central differences of ``dipole(P) -> (3,)`` (e.A), one placed atom at a time."""
    Pn = np.asarray(Pn, dtype=float)
    N = Pn.shape[0]
    Z = np.zeros((N, 3, 3))
    for i in range(N):
        for b in range(3):
            Pp, Pm = Pn.copy(), Pn.copy()
            Pp[i, b] += h
            Pm[i, b] -= h
            Z[i, :, b] = (np.asarray(dipole(Pp), dtype=float) - np.asarray(dipole(Pm), dtype=float)) / (2.0 * h)
    return Z


def electronic_dielectric(packer, params, Pn, latn, step: float = 1e-2) -> np.ndarray:
    """The clamped-ion dielectric tensor ``eps_inf`` (3, 3) of the cell at fixed nuclei.

    ``alpha_ab = d mu_a / d F_b = -(1 / EV_TO_KCAL) d^2 E / dF_a dF_b`` from second differences
    of :func:`polyfind.phonon.cell_energy_and_grad`'s energy in the applied field (exactly
    quadratic at fixed nuclei), and ``eps_inf = 1 + alpha / (eps0 V)``.  Tinfoil boundaries make
    the applied field the macroscopic one.  Without induced dipoles it is the identity.
    """
    from .phonon import cell_energy_and_grad

    params = np.asarray(params, dtype=float).reshape(7)
    Pn = np.asarray(Pn, dtype=float).reshape(packer.N, 3)
    latn = np.asarray(latn, dtype=float).reshape(3, 3)
    old = getattr(packer, "field", None)
    old = None if old is None else np.array(old, dtype=float)

    def e(F):
        packer.set_field(np.asarray(F, dtype=float))
        return cell_energy_and_grad(packer, params, Pn, latn)[0]

    try:
        E0 = e(np.zeros(3))
        H = np.zeros((3, 3))
        I = np.eye(3) * step
        for a in range(3):
            H[a, a] = (e(I[a]) - 2.0 * E0 + e(-I[a])) / step ** 2
            for b in range(a):
                H[a, b] = H[b, a] = (e(I[a] + I[b]) - e(I[a] - I[b]) - e(-I[a] + I[b]) + e(-I[a] - I[b])) / (4.0 * step ** 2)
    finally:
        packer.set_field(old)
    alpha = -H / EV_TO_KCAL  # e.A per V/A
    V = abs(float(np.linalg.det(latn)))
    return np.eye(3) + alpha * INV_EPS0_V_A / V


def lo_shift(freq_cm1: float, d: np.ndarray, eps_inf: np.ndarray, volume: float) -> float:
    """First-order LO frequency (cm^-1) of an isolated mode, for a wavevector along its own dipole.

    ``omega_LO^2 = omega_TO^2 + (4 pi / V) COULOMB |d|^2 / (d_hat . eps_inf . d_hat)`` in
    kcal/(mol A^2 amu), ``d`` in e/sqrt(amu) and ``V`` in A^3.  Ignores the coupling between
    polar modes that the full non-analytic term would add; an estimate of the size of what the
    tinfoil Hessian leaves out, not a computed LO spectrum.
    """
    from .phonon import SQRT_KCAL_A2_AMU_TO_CM1

    d = np.asarray(d, dtype=float)
    n2 = float(d @ d)
    if n2 == 0.0:
        return float(freq_cm1)
    u = d / np.sqrt(n2)
    lam_to = np.sign(freq_cm1) * (freq_cm1 / SQRT_KCAL_A2_AMU_TO_CM1) ** 2
    lam_lo = lam_to + 4.0 * np.pi * COULOMB * n2 / (volume * float(u @ eps_inf @ u))
    return float(np.sign(lam_lo) * np.sqrt(abs(lam_lo)) * SQRT_KCAL_A2_AMU_TO_CM1)


# ------------------------------------------------------------------------- the symmetry
@dataclass
class SymmetryOperation:
    """``P -> P R^T + t`` maps atom ``i`` onto atom ``perm[i]`` (modulo the lattice) to within ``error`` A."""

    R: np.ndarray  # (3, 3) Cartesian
    t: np.ndarray  # (3,) Cartesian, reduced into the cell
    perm: np.ndarray  # (N,)
    error: float

    @property
    def kind(self) -> str:
        det, tr = float(np.linalg.det(self.R)), float(np.trace(self.R))
        if det > 0:
            return {3: "identity", -1: "2-fold"}.get(int(round(tr)), f"rotation (trace {tr:+.0f})")
        return {-3: "inversion", 1: "mirror"}.get(int(round(tr)), f"rotoinversion (trace {tr:+.0f})")

    @property
    def axis(self) -> np.ndarray | None:
        """The two-fold axis (eigenvalue +1 of ``R``) or the mirror's normal (eigenvalue -1); ``None`` otherwise."""
        k = self.kind
        if k not in ("2-fold", "mirror"):
            return None
        w, v = np.linalg.eigh(0.5 * (self.R + self.R.T))
        idx = int(np.argmax(w)) if k == "2-fold" else int(np.argmin(w))
        a = v[:, idx]
        return a * (np.sign(a[np.argmax(np.abs(a))]) or 1.0)

    def describe(self) -> str:
        ax = self.axis
        s = self.kind if ax is None else f"{self.kind} {'along' if self.kind == '2-fold' else 'normal'} " \
                                           f"({ax[0]:+.3f}, {ax[1]:+.3f}, {ax[2]:+.3f})"
        return f"{s}, t = ({self.t[0]:+.4f}, {self.t[1]:+.4f}, {self.t[2]:+.4f}) A, atom-match error {self.error:.1e} A"

    def apply(self, v: np.ndarray) -> np.ndarray:
        """The operation on a Gamma-point displacement field ``v`` (3N,): ``(g v)_{perm[i]} = R v_i``."""
        V = np.asarray(v, dtype=float).reshape(-1, 3)
        out = np.empty_like(V)
        out[self.perm] = V @ self.R.T
        return out.ravel()


def _lattice_rotations(latn: np.ndarray, tol: float) -> np.ndarray:
    """Orthogonal ``R`` with ``latn R^T = M latn`` for an integer ``M`` with entries in {-1, 0, 1}."""
    L = np.asarray(latn, dtype=float)
    Linv = np.linalg.inv(L)
    Ms = np.array(list(itertools.product((-1, 0, 1), repeat=9)), dtype=float).reshape(-1, 3, 3)
    Ms = Ms[np.abs(np.abs(np.linalg.det(Ms)) - 1.0) < 0.5]
    R = np.transpose(Linv[None] @ Ms @ L[None], (0, 2, 1))
    ok = np.abs(R @ np.transpose(R, (0, 2, 1)) - np.eye(3)[None]).max(axis=(1, 2)) < tol
    return R[ok]


def space_group(Pn, latn, elements, tol: float = 1e-3) -> list:
    """Every :class:`SymmetryOperation` of the placed cell with atom-match error below ``tol`` (A).

    The lattice's own rotations are enumerated (integer matrices on the lattice rows), each is
    tried with every translation that takes one atom onto another of its element, and an
    operation is kept when every atom lands on an atom of the same element, modulo the
    lattice, within ``tol``.  The identity comes first.  Gamma-point modes see only ``R`` and
    the permutation; ``t`` is reported for naming (screw and glide parts).
    """
    P = np.asarray(Pn, dtype=float)
    L = np.asarray(latn, dtype=float)
    Linv = np.linalg.inv(L)
    el = np.asarray(list(elements))
    N = P.shape[0]
    same = el[:, None] == el[None, :]
    counts = {e: int((el == e).sum()) for e in set(el)}
    a0 = int(np.argmin([counts[e] for e in el]))  # the rarest element's first atom anchors the translations
    ops, seen = [], set()
    for R in _lattice_rotations(L, 1e-6):
        PR = P @ R.T
        for j in np.flatnonzero(el == el[a0]):
            t = P[j] - PR[a0]
            D = (PR + t)[:, None, :] - P[None, :, :]
            F = D @ Linv
            F -= np.round(F)
            dist = np.linalg.norm(F @ L, axis=-1)
            dist[~same] = np.inf
            perm = np.argmin(dist, axis=1)
            err = float(dist[np.arange(N), perm].max())
            if err > tol or len(set(perm.tolist())) != N:
                continue
            key = (tuple(np.round(R, 6).ravel()), tuple(perm.tolist()))
            if key in seen:
                continue
            seen.add(key)
            tf = t @ Linv
            tf -= np.floor(tf + 1e-9)
            ops.append(SymmetryOperation(R=R, t=tf @ L, perm=perm, error=err))
    ops.sort(key=lambda o: (o.kind != "identity", o.kind, o.error))
    return ops


def point_group_name(ops) -> str:
    """The point group of the distinct ``R`` among ``ops``: ``C1``, ``Ci``, ``C2``, ``Cs``, ``C2h``, ``C2v``,
    ``D2``, ``D2h``, else ``other (...)``.  A centring translation (a second identity) does not count."""
    distinct = {}
    for o in ops:
        distinct.setdefault(tuple(np.round(o.R, 6).ravel()), o.kind)
    kinds = tuple(sorted(distinct.values()))
    table = {("identity",): "C1", ("identity", "inversion"): "Ci", ("2-fold", "identity"): "C2",
             ("identity", "mirror"): "Cs", ("2-fold", "identity", "inversion", "mirror"): "C2h",
             ("2-fold", "identity", "mirror", "mirror"): "C2v", ("2-fold", "2-fold", "2-fold", "identity"): "D2",
             ("2-fold", "2-fold", "2-fold", "identity", "inversion", "mirror", "mirror", "mirror"): "D2h"}
    return table.get(kinds, f"other ({len(kinds)} distinct rotations: {', '.join(kinds)})")


def characters(ops, modes) -> np.ndarray:
    """``chi[o, k] = e_k . (g_o e_k)`` (n_ops, M) for mass-weighted eigenvectors ``modes`` (3N, M).

    Symmetry-equivalent atoms have equal masses, so the operation acts on mass-weighted and
    Cartesian displacements alike.  For a non-degenerate mode every character is +-1.
    """
    E = np.asarray(modes, dtype=float)
    return np.array([[float(E[:, k] @ op.apply(E[:, k])) for k in range(E.shape[1])] for op in ops])


def allowed_axes(ops, chi) -> list | None:
    """Per mode, the Cartesian axes (``"x"``, ``"y"``, ``"z"``) its dipole may have by symmetry; ``""`` if none.

    A mode can carry a dipole along axis ``a`` only if it transforms as the coordinate ``a``
    does: ``chi[o, k] = R_o[a, a]`` for every operation.  That reading needs every ``R`` to be
    diagonal (operations along the cell's own axes, as for an orthogonal cell); otherwise
    ``None``.  The selection rule the computed intensities must obey.
    """
    if not all(np.allclose(o.R, np.diag(np.diag(o.R)), atol=1e-6) for o in ops):
        return None
    out = []
    for k in range(chi.shape[1]):
        out.append("".join(ax for a, ax in enumerate("xyz")
                           if all(abs(chi[o, k] - ops[o].R[a, a]) < 0.5 for o in range(len(ops)))))
    return out


def irrep_labels(ops, chi) -> list:
    """Per-mode labels from the characters.

    ``C2h``: ``Ag``/``Bg``/``Au``/``Bu``.  ``C2v``: ``A1``, ``A2``, or ``B`` followed by the axis
    it transforms like (``By``, ``Bz``, ... -- convention-free, unlike ``B1``/``B2``).  Any other
    group with an inversion: ``g``/``u``; without one: ``?``.
    """
    kinds = [o.kind for o in ops]
    inv = kinds.index("inversion") if "inversion" in kinds else None
    two = kinds.index("2-fold") if "2-fold" in kinds else None
    group = point_group_name(ops)
    axes = allowed_axes(ops, chi) if group == "C2v" else None
    out = []
    for k in range(chi.shape[1]):
        if group == "C2v":
            mirrors = [chi[o, k] for o in range(len(ops)) if kinds[o] == "mirror"]
            if chi[two, k] > 0:
                out.append("A1" if all(m > 0 for m in mirrors) else "A2")
            else:
                out.append("B" + (axes[k] if axes is not None else "?"))
            continue
        if inv is None:
            out.append("?")
            continue
        p = "g" if chi[inv, k] > 0 else "u"
        out.append((("A" if chi[two, k] > 0 else "B") + p) if group == "C2h" else p)
    return out


def born_symmetry_error(ops, Z) -> float:
    """``max |Z[perm i] - R Z[i] R^T|`` over the operations: the Born tensors' own symmetry check (e)."""
    worst = 0.0
    for op in ops:
        Zr = np.einsum("ab,nbc,dc->nad", op.R, Z, op.R)
        worst = max(worst, float(np.abs(Z[op.perm] - Zr).max()))
    return worst


# ------------------------------------------------------------------------- the intensities
def mode_dipoles(Z, modes, masses) -> np.ndarray:
    """``d[k] = sum_{i, b} Z[i, :, b] e_k[3i + b] / sqrt(m_i)`` (M, 3) in e/sqrt(amu).

    ``Z`` (N, 3, 3) in e, ``modes`` (3N, M) orthonormal mass-weighted eigenvectors as columns,
    ``masses`` (N,) in amu.
    """
    Z = np.asarray(Z, dtype=float)
    E = np.asarray(modes, dtype=float)
    N = Z.shape[0]
    U = E.reshape(N, 3, -1) / np.sqrt(np.asarray(masses, dtype=float))[:, None, None]
    return np.einsum("iab,ibk->ka", Z, U)


@dataclass
class InfraredMode:
    """One Gamma-point mode's infrared activity.

    ``dipole`` (3,) is ``d mu / d Q`` in e/sqrt(amu) in the packer's frame (x polar, y the long
    lateral axis, z the chain axis); ``intensity`` ``|dipole|^2`` in e^2/amu and ``km_per_mol``
    the same per mole of cells; ``shares`` the fraction of ``|dipole|^2`` along x, y, z;
    ``parity`` the character under inversion (``None`` without one), ``irrep`` the label from
    :func:`irrep_labels`.
    """

    index: int
    freq_cm1: float
    dipole: np.ndarray
    intensity: float
    km_per_mol: float
    parity: float | None = None
    irrep: str = "?"
    allowed: str | None = None  # axes the dipole may have by symmetry (allowed_axes), None if not determined

    @property
    def shares(self) -> np.ndarray:
        n2 = float(self.dipole @ self.dipole)
        return self.dipole ** 2 / n2 if n2 > 0 else np.zeros(3)

    @property
    def parallel(self) -> float:
        """Share of the intensity polarised along the chain axis."""
        return float(self.shares[2])


@dataclass
class Infrared:
    """Infrared activity of every mode of a :class:`polyfind.phonon.Phonons`, with the checks.

    ``asr`` is ``sum_i Z*_i`` (3, 3); ``ops`` the cell's symmetry operations; ``chi`` their
    characters on every mode; ``born_symmetry`` :func:`born_symmetry_error`.
    """

    modes: list
    Z: np.ndarray
    asr: np.ndarray
    ops: list = field(default_factory=list)
    chi: np.ndarray | None = None
    born_symmetry: float = float("nan")

    @property
    def point_group(self) -> str:
        return point_group_name(self.ops) if self.ops else "not determined"

    @property
    def centrosymmetric(self) -> bool:
        return any(o.kind == "inversion" for o in self.ops)

    def g_leak(self, acoustic=()) -> float:
        """Largest intensity on an inversion-even mode, relative to the largest intensity of all."""
        top = max(m.intensity for m in self.modes)
        g = [m.intensity for m in self.modes if m.parity is not None and m.parity > 0 and m.index not in acoustic]
        return (max(g) / top) if (g and top > 0) else 0.0

    def forbidden_leak(self, acoustic=()) -> float:
        """Largest share of any mode's ``|d|^2`` along an axis its symmetry forbids, weighted by that mode's
        intensity relative to the strongest: zero when the intensities obey the selection rule."""
        top = max(m.intensity for m in self.modes)
        worst = 0.0
        for m in self.modes:
            if m.allowed is None or m.index in acoustic or top <= 0:
                continue
            bad = sum(float(m.dipole[a] ** 2) for a, ax in enumerate("xyz") if ax not in m.allowed)
            worst = max(worst, bad / top)
        return worst

    def table(self, fmax: float = 200.0, acoustic=(), ref: float | None = None) -> str:
        """Rows of every mode up to ``fmax`` cm^-1: frequency, irrep, intensity (km/mol per cell and
        relative to ``ref`` km/mol, default the strongest optical mode in the window), the
        dipole direction and the share along the chain axis."""
        rows = [m for m in self.modes if m.freq_cm1 <= fmax and m.index not in acoustic]
        if ref is None:
            ref = max((m.km_per_mol for m in rows), default=1.0) or 1.0
        out = [f"{'cm^-1':>8s} {'irrep':>5s} {'allow':>5s} {'km/mol':>9s} {'rel':>6s}  {'x':>5s} {'y':>5s} {'z':>5s}  "
               "dipole d (e/sqrt(amu))"]
        for m in rows:
            s = m.shares
            allow = "?" if m.allowed is None else (m.allowed or "-")
            out.append(f"{m.freq_cm1:8.2f} {m.irrep:>5s} {allow:>5s} {m.km_per_mol:9.3f} {m.km_per_mol / ref:6.3f}  "
                       f"{s[0]:5.2f} {s[1]:5.2f} {s[2]:5.2f}  ({m.dipole[0]:+.4f}, {m.dipole[1]:+.4f}, {m.dipole[2]:+.4f})")
        return "\n".join(out)


def infrared(phonons, Z, ops=None) -> Infrared:
    """:class:`Infrared` for ``phonons`` (a :class:`polyfind.phonon.Phonons`) and Born tensors ``Z`` (N, 3, 3).

    ``ops`` (from :func:`space_group` at the Hessian's geometry) adds parities and labels; the
    intensities do not depend on it.
    """
    Z = np.asarray(Z, dtype=float)
    d = mode_dipoles(Z, phonons.modes, phonons.masses)
    chi = labels = axes = None
    inv = None
    if ops:
        chi = characters(ops, phonons.modes)
        labels = irrep_labels(ops, chi)
        axes = allowed_axes(ops, chi)
        kinds = [o.kind for o in ops]
        inv = kinds.index("inversion") if "inversion" in kinds else None
    modes = []
    for k in range(d.shape[0]):
        I2 = float(d[k] @ d[k])
        modes.append(InfraredMode(index=k, freq_cm1=float(phonons.freq_cm1[k]), dipole=d[k].copy(), intensity=I2,
                                  km_per_mol=I2 * KM_PER_MOL_PER_E2_AMU,
                                  parity=None if inv is None else float(chi[inv, k]),
                                  irrep="?" if labels is None else labels[k],
                                  allowed=None if axes is None else axes[k]))
    return Infrared(modes=modes, Z=Z, asr=Z.sum(axis=0), ops=list(ops or []), chi=chi,
                    born_symmetry=born_symmetry_error(ops, Z) if ops else float("nan"))
