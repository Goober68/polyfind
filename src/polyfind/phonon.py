r"""Gamma-point lattice dynamics through the packer's complete placed-cell owner.

Canonical packing/refinement, independent-atom phonon gradients and the
Cartesian Hessian consume one Hamiltonian. cell_energy_and_grad only
forwards the cell and discrete chain reversal to
CrystalPacker.placed_energy_and_grad; it does not assemble another potential
or reconstruct chain-local frames.

The owner includes short-range nonbonded interactions with complete triclinic
cutoff images/topological exclusions, valence, geometry-derived periodic
torsion, total charge-flux derivatives, Ewald/exclusions, stationary induced
dipoles and an applied field. Torsional forces and stiffness now enter this
Cartesian Hessian. The former value-only torsion treatment and the separate
potential assembly are removed.

Gamma-only support remains a physical sampling limitation. A two-backbone-site
primitive trans repeat cannot represent chain twisting; larger independent
repeats/BZ convergence are required. A Hessian at a nonstationary geometry
does not establish stable lattice modes. Model frequencies are not calibrated
material frequencies or general frequency bounds.


**Units.**  The Hessian is in kcal/(mol A^2); masses are :data:`polyfind.pack.MASS` in g/mol
(= amu).  With ``1 kcal = 4184 J`` (exact), ``N_A = 6.02214076e23 1/mol`` (exact, SI 2019),
``1 amu = 1.66053906660e-27 kg`` (CODATA 2018) and ``c = 2.99792458e10 cm/s`` (exact), an
eigenvalue ``lambda`` of the mass-weighted Hessian in kcal/(mol A^2 amu) is an angular
frequency squared ``omega^2 = lambda * 4184 / (N_A amu 1e-20)`` in rad^2/s^2, and the
wavenumber is ``omega / (2 pi c)``: ``sqrt(1 kcal/(mol A^2 amu)) = 108.59 cm^-1``
(:data:`SQRT_KCAL_A2_AMU_TO_CM1`).  A negative eigenvalue is reported as a negative
wavenumber, ``-sqrt(-lambda)``, the usual convention for an imaginary mode.

**Acoustic sum rule.**  Every term above is invariant under a uniform translation of the
whole cell (the Ewald sum under tinfoil boundaries included), so ``sum_j H[i a, j b] = 0`` at
any geometry -- it is a check on the arithmetic, and :attr:`Phonons.asr_residual` reports it
per atom.  ``asr="project"`` removes the three uniform translations from the dynamical matrix
exactly (the projector is applied on both sides), which pins the acoustic frequencies to zero
and leaves every optical mode where the raw matrix put it to within the residual.  The cell
is periodic in all three directions and has two chains, so exactly three zero modes are
expected at Gamma and no more; if more appear, or if any mode is imaginary, that is a
property of the geometry (not at a minimum over the all-atom coordinates) or of the
potential, and it is reported, not projected away.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from . import backend as bk
from .born import atom_labels
from .pack import MASS

# --- unit conversion, derived from the constants written down in the module docstring ------
KCAL_PER_MOL_TO_J = 4184.0 / 6.02214076e23  # J per (kcal/mol); 1 kcal = 4184 J, N_A exact
AMU_TO_KG = 1.66053906660e-27  # CODATA 2018
C_CM_PER_S = 2.99792458e10  # exact
# omega^2 in rad^2/s^2 per (kcal/mol/A^2/amu): energy / (mass length^2)
OMEGA2_PER_KCAL_A2_AMU = KCAL_PER_MOL_TO_J / (AMU_TO_KG * 1e-20)
SQRT_KCAL_A2_AMU_TO_CM1 = float(np.sqrt(OMEGA2_PER_KCAL_A2_AMU) / (2.0 * np.pi * C_CM_PER_S))  # 108.59
SQRT_KCAL_A2_AMU_TO_THZ = float(np.sqrt(OMEGA2_PER_KCAL_A2_AMU) / (2.0 * np.pi) * 1e-12)  # 3.2556


def _rotz_matrix(deg: float) -> np.ndarray:
    c, s = np.cos(np.deg2rad(deg)), np.sin(np.deg2rad(deg))
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def placing_frames(params, n_chains: int) -> list:
    """``R_s`` (3, 3) per chain with ``P_s = X_s @ R_s.T + offset_s``: the linear part of
    :meth:`polyfind.pack.CrystalPacker._place`, i.e. ``Rz(phi1)`` and ``Rz(phi2) . flip``."""
    a, b, gam, phi1, phi2, dz, flip = (float(v) for v in np.asarray(params, dtype=float).reshape(7))
    R = [_rotz_matrix(phi1)]
    if n_chains == 2:
        F = np.diag([1.0, -1.0, -1.0]) if flip > 0.5 else np.eye(3)
        R.append(_rotz_matrix(phi2) @ F)
    return R


def cell_energy_and_grad(packer, params, Pn, latn) -> tuple:
    """(E,gatoms) from the packer's complete independent placed-cell owner.

    Only the discrete chain reversal is supplied by canonical params.
    No setting-angle/offset reconstruction or parallel potential assembly.
    """
    params = np.asarray(params,dtype=float)
    if params.shape != (7,) or not np.isfinite(params).all():
        raise ValueError("phonon requires one finite canonical parameter row")
    E,gP,_ = packer.placed_energy_and_grad(Pn,latn,params[6])
    return E,gP



def cell_energy(packer, params, Pn, latn) -> float:
    """The energy half of :func:`cell_energy_and_grad`: kcal/mol per cell at placed coordinates ``Pn``."""
    return cell_energy_and_grad(packer, params, Pn, latn)[0]


def force_terms(packer, params, Pn, latn) -> dict:
    """Where the all-atom force sits: ``max |dE/dP|`` per term, kcal/(mol A), plus the bond and angle strain.

    ``bonds`` and ``angles`` are the valence stretch and bend gradients alone (each chain in its
    own frame, mapped to the placed one); ``other`` is everything else -- nonbonded, Ewald,
    exclusions, induced dipoles, flux and field -- as the total minus the two.  ``bond_strain``
    is ``max |r - r0|`` (A) and ``angle_strain`` ``max |theta - theta0|`` (deg) over the
    chains' valence terms.  The packer's own reference state relaxes cell, setting angles,
    torsions and backbone angles but never a bond length
    (:class:`polyfind.mechanics.Shape`), so a residual force here is ordinarily the stretch
    term of a bond built at the polymer's length and fitted to another ``r0``.
    """
    from dataclasses import replace as _replace

    params = np.asarray(params, dtype=float).reshape(7)
    n, nc = packer.n, packer.n_chains
    Pn = np.asarray(Pn, dtype=float).reshape(packer.N, 3)
    latn = np.asarray(latn, dtype=float).reshape(3, 3)
    _, g_all = cell_energy_and_grad(packer, params, Pn, latn)
    out = {"total": float(np.abs(g_all).max()), "bonds": 0.0, "angles": 0.0, "other": float(np.abs(g_all).max()),
           "bond_strain": 0.0, "angle_strain": 0.0}
    val = packer._valence
    if val is None:
        return out
    ints = lambda: np.zeros(0, dtype=int)  # noqa: E731
    reals = lambda: np.zeros(0)  # noqa: E731
    bonds_only = _replace(val, ang_i=ints(), ang_j=ints(), ang_k=ints(), ang_si=reals(), ang_sk=reals(), ang_kk=reals(), ang_t0=reals())
    angles_only = _replace(val, bond_i=ints(), bond_j=ints(), bond_s=reals(), bond_k=reals(), bond_r0=reals())
    gb, ga = np.zeros_like(Pn), np.zeros_like(Pn)
    for s in range(nc):
        sl = slice(s * n, (s + 1) * n)
        X = Pn[sl]
        repeat = latn[2]*(-1. if s == 1 and params[6] == 1. else 1.)
        gb[sl] = bonds_only.energy_and_repeat_grad(X, repeat, n_atoms=n)[1]
        ga[sl] = angles_only.energy_and_repeat_grad(X, repeat, n_atoms=n)[1]
        if val.n_bonds:
            r = val.bond_lengths(X, repeat=repeat)[0]
            out["bond_strain"] = max(out["bond_strain"], float(np.abs(r - val.bond_r0).max()))
        if val.n_angles:
            angles = val.angle_values(X, repeat=repeat)[0]
            out["angle_strain"] = max(out["angle_strain"], float(np.rad2deg(np.abs(angles - val.ang_t0)).max()))
    out["bonds"] = float(np.abs(gb).max())
    out["angles"] = float(np.abs(ga).max())
    out["other"] = float(np.abs(g_all - gb - ga).max())
    return out


def packer_with_built_bond_lengths(packer):
    """A shallow copy of ``packer`` whose valence stretch terms have ``r0`` pinned to the chain's own bond lengths.

    The packer's reference states are relaxed over cell, setting angles, torsions and backbone
    angles with every bond length rigid at the value :func:`polyfind.chain.build_chain` gave
    it, and the valence fit (``pvdf-dft-valence``) was made on such rigid-bond geometries: its
    stretch ``r0`` values are not bond potentials (C-F 1.467 A against a built 1.350 A, C-C
    1.458 A against 1.528 A, ``k(C-C)`` at the fit's bound) and the packer's own path never
    feels them.  Let every atom go against them and the C-F bond walks to 1.54 A.  Setting
    ``r0`` to the built lengths keeps the stiffness ``k`` and makes the reference geometry
    stationary over the stretches, so an all-atom relaxation moves the atoms by a few
    hundredths of an Angstrom rather than the bond lengths.  The energy is shifted by a constant
    (the stretch energy at the reference) and the copy is a consistent packer:
    :func:`check_known_answer` holds on it as on any other.  Nothing about ``packer`` changes.
    """
    import copy
    from dataclasses import replace as _replace

    val = packer._valence
    if val is None or not val.n_bonds:
        return packer
    X = np.asarray(packer.chain.coords, dtype=float)
    r = val.bond_lengths(X, float(packer.chain.c))[0]
    pk = copy.copy(packer)
    pk._valence = _replace(val, bond_r0=r.copy())
    pk.update_chain(packer.chain)
    return pk


def placed_coordinates(packer, params) -> tuple:
    """``(Pn (N, 3), latn (3, 3))`` of the packer's own chain at ``params``, host float64."""
    P, lat = packer._place(np.asarray(params, dtype=float).reshape(1, 7))
    return np.asarray(bk.to_numpy(P), dtype=float)[0], np.asarray(bk.to_numpy(lat), dtype=float)[0]


def fold_gradient(packer, params, gP) -> np.ndarray:
    """The all-atom gradient folded onto one repeat: ``sum_s R_s^T dE/dP_s``, the ``g_coords`` of
    :meth:`~polyfind.pack.CrystalPacker.energy_and_grad` (a repeat displacement moves every chain)."""
    n = packer.n
    R = placing_frames(params, packer.n_chains)
    gP = np.asarray(gP, dtype=float).reshape(packer.N, 3)
    return sum(gP[s * n:(s + 1) * n] @ R[s] for s in range(packer.n_chains))


def check_known_answer(packer, params) -> tuple:
    """``(E_cell, E_packer, |difference| / |E_packer|)`` at the undisplaced placed coordinates.

    The identity every number in this module rests on: the assembled energy must be the
    packer's, not a model of it.  Raises if the relative difference exceeds ``1e-9``.
    """
    params = np.asarray(params, dtype=float).reshape(7)
    Pn, latn = placed_coordinates(packer, params)
    e_cell = cell_energy(packer, params, Pn, latn)
    e_pack = float(packer.energy(params[None])[0])
    rel = abs(e_cell - e_pack) / max(abs(e_pack), 1e-300)
    if rel > 1e-9:
        raise RuntimeError(f"cell_energy {e_cell:.12g} differs from packer.energy {e_pack:.12g} "
                           f"by {rel:.2e} relative: a term is missing or miscounted")
    return e_cell, e_pack, rel


# ------------------------------------------------------------------------- the Hessian
@dataclass
class Hessian:
    """The Cartesian Gamma-point Hessian ``H`` (3N, 3N) in kcal/(mol A^2), symmetrised.

    ``asymmetry`` is ``|H - H^T|_F / |H|_F`` before symmetrisation -- the finite-difference
    diagnostic, since the exact Hessian is symmetric; ``gradient`` is the all-atom gradient at
    the geometry (kcal/(mol A)), whose size says how far from a stationary point the cell is;
    ``energy`` the known-answer energy.  Unpacks as ``(H, Pn, latn)``.
    """

    H: np.ndarray
    Pn: np.ndarray
    latn: np.ndarray
    asymmetry: float
    asymmetry_abs: float
    energy: float
    gradient: np.ndarray
    h: float

    def __iter__(self):
        return iter((self.H, self.Pn, self.latn))

    @property
    def max_force(self) -> float:
        return float(np.abs(self.gradient).max())


def gamma_hessian(packer, params, h: float = 1e-3, Pn=None, latn=None) -> Hessian:
    """Central differences of the analytic all-atom gradient, one placed coordinate at a time.

    ``H[:, k] = (g(P + h e_k) - g(P - h e_k)) / 2h``, then symmetrised.  Every displaced
    evaluation re-solves the fluxing charges of the displaced chain and the induced dipoles of
    the whole cell, so the Hessian is that of the *total* energy.  ``Pn``/``latn`` default to
    the packer's own chain placed at ``params``; pass them to take the Hessian at a geometry
    of your own (a relaxed one, say).
    """
    params = np.asarray(params, dtype=float).reshape(7)
    if Pn is None or latn is None:
        Pn, latn = placed_coordinates(packer, params)
    Pn = np.asarray(Pn, dtype=float).reshape(packer.N, 3)
    latn = np.asarray(latn, dtype=float).reshape(3, 3)
    N = packer.N
    E0, g0 = cell_energy_and_grad(packer, params, Pn, latn)
    H = np.zeros((3 * N, 3 * N))
    for k in range(3 * N):
        i, d = divmod(k, 3)
        Pp, Pm = Pn.copy(), Pn.copy()
        Pp[i, d] += h
        Pm[i, d] -= h
        gp = cell_energy_and_grad(packer, params, Pp, latn)[1]
        gm = cell_energy_and_grad(packer, params, Pm, latn)[1]
        H[:, k] = (gp - gm).ravel() / (2.0 * h)
    asym_abs = float(np.linalg.norm(H - H.T))
    asym = asym_abs / max(float(np.linalg.norm(H)), 1e-300)
    H = 0.5 * (H + H.T)
    return Hessian(H=H, Pn=Pn, latn=latn, asymmetry=asym, asymmetry_abs=asym_abs, energy=E0, gradient=g0, h=float(h))


def hessian_by_energy(packer, params, h: float = 1e-3, Pn=None, latn=None, pairs=None) -> np.ndarray:
    """Selected entries of the Hessian from second differences of the *energy* alone.

    ``H[k, l] = [E(+k+l) - E(+k-l) - E(-k+l) + E(-k-l)] / 4h^2`` (``k == l`` uses the three-point
    stencil).  Independent of the analytic gradient, so it is the check on it; ``pairs`` lists
    the ``(k, l)`` to evaluate (all of them by default, which is ``O((3N)^2)`` energies).
    Returns a (3N, 3N) array with ``nan`` where nothing was asked for.
    """
    params = np.asarray(params, dtype=float).reshape(7)
    if Pn is None or latn is None:
        Pn, latn = placed_coordinates(packer, params)
    Pn = np.asarray(Pn, dtype=float).reshape(packer.N, 3)
    M = 3 * packer.N
    out = np.full((M, M), np.nan)
    e = lambda P: cell_energy(packer, params, P, latn)  # noqa: E731
    E0 = None
    if pairs is None:
        pairs = [(k, l) for k in range(M) for l in range(k, M)]
    for k, l in pairs:
        if k == l:
            if E0 is None:
                E0 = e(Pn)
            Pp, Pm = Pn.copy(), Pn.copy()
            Pp.flat[k] += h
            Pm.flat[k] -= h
            out[k, k] = (e(Pp) - 2.0 * E0 + e(Pm)) / (h * h)
            continue
        vals = []
        for sk, sl_ in ((1, 1), (1, -1), (-1, 1), (-1, -1)):
            P = Pn.copy()
            P.flat[k] += sk * h
            P.flat[l] += sl_ * h
            vals.append(e(P))
        out[k, l] = out[l, k] = (vals[0] - vals[1] - vals[2] + vals[3]) / (4.0 * h * h)
    return out


# ------------------------------------------------------------------------- phonons
@dataclass
class Mode:
    """One Gamma-point mode: its frequency and where its (mass-weighted) amplitude sits."""

    index: int  # position in the sorted list
    freq_cm1: float
    freq_thz: float
    by_type: dict  # label -> share of |e|^2
    by_element: dict  # element -> share of |e|^2
    axis: np.ndarray  # (3,) share along x (polar), y (long lateral), z (chain axis)
    translation: float  # share along the three uniform translations
    rigid_chain: float  # share that is a rigid translation / rotation-about-axis of each chain
    by_chain: np.ndarray  # (n_chains,) share per chain

    @property
    def transverse(self) -> float:
        return float(self.axis[0] + self.axis[1])

    @property
    def axial(self) -> float:
        return float(self.axis[2])

    def describe(self) -> str:
        typ = " ".join(f"{k}:{v:.2f}" for k, v in self.by_type.items())
        return (f"{self.freq_cm1:8.2f} cm^-1 ({self.freq_thz:6.3f} THz)  x/y/z {self.axis[0]:.2f}/{self.axis[1]:.2f}/"
                f"{self.axis[2]:.2f}  transverse {self.transverse:.2f}  rigid-chain {self.rigid_chain:.2f}  "
                f"chains {'/'.join(f'{v:.2f}' for v in self.by_chain)}  {typ}")


@dataclass
class Phonons:
    """Gamma-point phonons of one packed cell.

    ``freq_cm1`` / ``freq_thz`` are sorted ascending, imaginary modes negative; ``modes`` holds
    the orthonormal eigenvectors of the mass-weighted Hessian as columns (Cartesian
    displacements are ``modes / sqrt(m)`` per atom); ``eigenvalues`` in kcal/(mol A^2 amu).
    ``asr_residual`` (N, 3, 3) is ``sum_j H[i a, j b]``, the force on atom ``i`` under a unit
    translation along ``b``; ``asr_max`` its largest entry, against ``hessian_scale`` the
    largest entry of ``H`` itself.  ``n_imaginary`` counts negative eigenvalues below
    ``-zero_tol``; ``n_zero`` the ones within ``+-zero_tol`` (three are the acoustic modes).
    ``acoustic`` are the indices of the three modes with the largest translational share, and
    ``lowest_optical`` the report on the lowest few of the rest.  With ``asr="project"``,
    ``raw`` carries the unprojected result.
    """

    freq_cm1: np.ndarray
    freq_thz: np.ndarray
    eigenvalues: np.ndarray
    modes: np.ndarray
    masses: np.ndarray
    elements: list
    labels: list
    hessian: Hessian
    asr_residual: np.ndarray
    asr_max: float
    hessian_scale: float
    n_imaginary: int
    n_zero: int
    zero_tol: float
    acoustic: list
    lowest_optical: list
    asr: str
    raw: "Phonons | None" = None
    notes: list = field(default_factory=list)
    _modes: list = field(default_factory=list, repr=False)
    _modes_optical: list = field(default_factory=list, repr=False)

    @property
    def optical(self) -> np.ndarray:
        """The sorted frequencies with the three acoustic modes removed (cm^-1)."""
        keep = [i for i in range(len(self.freq_cm1)) if i not in self.acoustic]
        return self.freq_cm1[keep]

    def mode(self, i: int) -> Mode:
        return self._modes[i]

    def report(self, n_lowest: int = 6) -> str:
        lines = [f"known-answer energy {self.hessian.energy:.6f} kcal/mol per cell; max |dE/dP| "
                 f"{self.hessian.max_force:.3e} kcal/(mol A); finite-difference step {self.hessian.h:g} A",
                 f"Hessian asymmetry before symmetrisation {self.hessian.asymmetry:.2e} (relative), "
                 f"{self.hessian.asymmetry_abs:.2e} kcal/(mol A^2) absolute; largest entry {self.hessian_scale:.4g}",
                 f"acoustic sum rule residual max |sum_j H_ij| = {self.asr_max:.3e} kcal/(mol A^2) "
                 f"({self.asr_max / self.hessian_scale:.1e} of the largest entry); asr='{self.asr}'",
                 f"imaginary modes: {self.n_imaginary}; modes within +-{self.zero_tol:g} cm^-1 of zero: {self.n_zero} "
                 f"(three acoustic expected); acoustic indices {self.acoustic}",
                 "frequencies (cm^-1): " + " ".join(f"{v:.1f}" for v in self.freq_cm1),
                 f"lowest {n_lowest} optical modes (shares of the mass-weighted eigenvector):"]
        for m in self._modes_optical[:n_lowest]:
            lines.append("  " + m.describe())
        lines.extend("note: " + s for s in self.notes)
        return "\n".join(lines)


def _mode_report(i, lam, e, masses, elements, labels, n, nc, Pn, T) -> Mode:
    w = (e.reshape(-1, 3) ** 2).sum(axis=1)  # per atom, sums to 1
    by_type, by_el = {}, {}
    for lab in dict.fromkeys(labels):
        by_type[lab] = float(sum(w[j] for j in range(len(labels)) if labels[j] == lab))
    for el in dict.fromkeys(elements):
        by_el[el] = float(sum(w[j] for j in range(len(elements)) if elements[j] == el))
    axis = (e.reshape(-1, 3) ** 2).sum(axis=0)
    trans = float(((T.T @ e) ** 2).sum())
    rigid = 0.0
    by_chain = np.zeros(nc)
    sm = np.sqrt(masses)
    for s in range(nc):
        sl = slice(s * n, (s + 1) * n)
        es = e.reshape(-1, 3)[sl].ravel()
        by_chain[s] = float(es @ es)
        # mass-weighted rigid motions of an infinite periodic chain: 3 translations and the
        # rotation about its own axis (through its mass centre)
        X = Pn[sl]
        com = (masses[sl][:, None] * X).sum(axis=0) / masses[sl].sum()
        basis = []
        for d in range(3):
            v = np.zeros((n, 3))
            v[:, d] = 1.0
            basis.append((v * sm[sl][:, None]).ravel())
        r = X - com
        v = np.stack([-r[:, 1], r[:, 0], np.zeros(n)], axis=1)
        basis.append((v * sm[sl][:, None]).ravel())
        Q, _ = np.linalg.qr(np.stack(basis, axis=1))
        rigid += float(((Q.T @ es) ** 2).sum())
    f = _freq(lam)
    return Mode(index=i, freq_cm1=f, freq_thz=f * SQRT_KCAL_A2_AMU_TO_THZ / SQRT_KCAL_A2_AMU_TO_CM1, by_type=by_type,
                by_element=by_el, axis=axis, translation=trans, rigid_chain=rigid, by_chain=by_chain)


def _freq(lam: float) -> float:
    return float(np.sign(lam) * np.sqrt(abs(lam)) * SQRT_KCAL_A2_AMU_TO_CM1)


def atom_masses(packer) -> np.ndarray:
    """(N,) masses in amu (:data:`polyfind.pack.MASS`), in the placed atom order."""
    return np.array([MASS[e] for e in packer.elements] * packer.n_chains, dtype=float)


def phonons_from_hessian(packer, hess: Hessian, asr: str = "none", zero_tol: float = 1.0, n_report: int = 6) -> Phonons:
    """:class:`Phonons` of a :class:`Hessian`; see :func:`phonons_gamma`."""
    if asr not in ("none", "project"):
        raise ValueError("asr must be 'none' or 'project'")
    N, n, nc = packer.N, packer.n, packer.n_chains
    masses = atom_masses(packer)
    elements = list(packer.elements) * nc
    labels = atom_labels(packer.chain) * nc
    H = hess.H
    scale = float(np.abs(H).max())
    asr_res = H.reshape(N, 3, N, 3).sum(axis=2)  # (N, 3, 3): force on i along a per unit translation along b
    inv = 1.0 / np.sqrt(masses)
    W = np.repeat(inv, 3)
    D = W[:, None] * H * W[None, :]
    # the three uniform translations, mass-weighted and orthonormal
    T = np.zeros((3 * N, 3))
    for d in range(3):
        T[d::3, d] = np.sqrt(masses)
    T /= np.linalg.norm(T, axis=0)[None, :]
    raw = None
    if asr == "project":
        raw = phonons_from_hessian(packer, hess, asr="none", zero_tol=zero_tol, n_report=n_report)
        Pj = np.eye(3 * N) - T @ T.T
        D = Pj @ D @ Pj
        D = 0.5 * (D + D.T)
    lam, V = np.linalg.eigh(D)
    freqs = np.array([_freq(v) for v in lam])
    thz = freqs * (SQRT_KCAL_A2_AMU_TO_THZ / SQRT_KCAL_A2_AMU_TO_CM1)
    modes = [_mode_report(i, lam[i], V[:, i], masses, elements, labels, n, nc, hess.Pn, T) for i in range(3 * N)]
    acoustic = sorted(sorted(range(3 * N), key=lambda i: -modes[i].translation)[:3])
    optical = [m for m in modes if m.index not in acoustic]
    lam_tol = (zero_tol / SQRT_KCAL_A2_AMU_TO_CM1) ** 2
    n_imag = int((lam < -lam_tol).sum())
    n_zero = int((np.abs(lam) <= lam_tol).sum())
    notes = ["the complete placed-cell model includes geometry-derived torsion; Gamma-only repeat support "
             "does not cover all chain-twisting branches, and model frequencies are not independently "
             "calibrated material frequencies (see polyfind.phonon)"]
    if hess.max_force > 1e-2:
        notes.append(f"the geometry is not a stationary point of the all-atom energy (max |dE/dP| = "
                     f"{hess.max_force:.2e} kcal/(mol A)); a Hessian there can have imaginary modes that a relaxed one would not")
    if n_zero != 3:
        notes.append(f"{n_zero} modes within +-{zero_tol:g} cm^-1 of zero where three (the acoustic modes) are expected")
    if n_imag:
        notes.append(f"{n_imag} imaginary mode(s): lowest {freqs[0]:.2f} cm^-1")
    ph = Phonons(freq_cm1=freqs, freq_thz=thz, eigenvalues=lam, modes=V, masses=masses, elements=elements,
                 labels=labels, hessian=hess, asr_residual=asr_res, asr_max=float(np.abs(asr_res).max()),
                 hessian_scale=scale, n_imaginary=n_imag, n_zero=n_zero, zero_tol=float(zero_tol),
                 acoustic=acoustic, lowest_optical=optical[:n_report], asr=asr, raw=raw, notes=notes,
                 _modes=modes, _modes_optical=optical)
    return ph


def phonons_gamma(packer, params, h: float = 1e-3, asr: str = "none", Pn=None, latn=None,
                  zero_tol: float = 1.0, n_report: int = 6) -> Phonons:
    """Gamma-point phonons of the packer's cell at ``params``.

    :func:`gamma_hessian` with step ``h`` (A), mass-weighted with :data:`polyfind.pack.MASS`,
    diagonalised; ``asr="project"`` removes the uniform translations first (the unprojected
    result is kept in :attr:`Phonons.raw`).  ``zero_tol`` (cm^-1) is the band around zero
    that counts a mode as acoustic-like.  The three lowest (and every) optical mode carry a
    per-type, per-axis and rigid-chain decomposition (:class:`Mode`).
    """
    hess = gamma_hessian(packer, params, h=h, Pn=Pn, latn=latn)
    return phonons_from_hessian(packer, hess, asr=asr, zero_tol=zero_tol, n_report=n_report)


def relax_all_atom(packer, params, Pn=None, latn=None, gtol: float = 1e-6, maxiter: int = 2000) -> tuple:
    """Independent atoms at fixed cell through the shared Cartesian relaxer.

    Returns (Pn,E,max absolute gradient component), preserving the legacy
    observation API. The optimizer uses a vector-norm force criterion gtol
    and an explicit mass-centre translation gauge, not post-hoc recentering.
    The tuple is not a stability or physical-trajectory claim.
    """
    from .cartesian_mechanics import CartesianChart,CartesianRelaxer,CartesianTolerance
    from .topology import build_cell

    params = np.asarray(params,dtype=float).reshape(7)
    cell = build_cell(packer.chain,params,packer=packer)
    if Pn is not None:
        cell.coords = np.asarray(Pn,dtype=float).reshape(packer.N,3)
    if latn is not None:
        cell.lattice = np.asarray(latn,dtype=float).reshape(3,3)
    chart = CartesianChart(cell,free_strain=np.zeros(6,dtype=bool))
    result = CartesianRelaxer(packer,chart,CartesianTolerance(gtol,1.)).run(maxiter=maxiter)
    geometry = result.geometry
    evaluation = packer.evaluate_chain_cell(geometry.coords,geometry.lattice,chart.layout)
    return geometry.coords.copy(),evaluation.terms.total,float(np.abs(evaluation.grad_coords).max())
