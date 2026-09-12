r"""Ewald summation of the lattice electrostatics.

The rest of this package evaluates the Coulomb energy of a crystal with a
damped-shifted-force sum truncated at 8 A (:meth:`polyfind.pack.CrystalPacker._pair_energy`).
That is fine for ranking neutral chains and it is *not* fine for the one question this
module exists to answer.  A lattice of dipoles interacts as ``1/r^3`` while the number of
pairs at distance ``r`` grows as ``r^2``, so the dipole-dipole part of the sum is only
conditionally convergent: its value depends on the order of summation, and **no cutoff
evaluates it**, whatever the cutoff is.  The polar/antipolar energy difference of a
ferroelectric is exactly that term.  ``docs/SCREEN.md``'s addendum is the measurement that
made this concrete.

**The decomposition.**  With ``E = (1/2) sum_{i,j,n}' q_i q_j / |r_ij + n|`` (the prime
skipping ``i = j`` at ``n = 0``), splitting ``1/r`` into ``erfc(a r)/r + erf(a r)/r`` gives

    E_real  = (1/2) sum_{i,j,n}' q_i q_j erfc(a r)/r                     (r < rc)
    E_recip = (2 pi / V) sum_{k != 0} exp(-k^2 / 4 a^2) / k^2 |S(k)|^2   (|k| < kmax)
    E_self  = -(a / sqrt(pi)) sum_i q_i^2
    E_bg    = -(pi / (2 V a^2)) Q^2                                       (Q = sum_i q_i)
    E_surf  = (2 pi / (3 V)) |M|^2                                        (M = sum_i q_i r_i)

with ``S(k) = sum_i q_i exp(i k . r_i)``.  The first two converge absolutely and the
splitting parameter ``a`` only moves work between them: the total is independent of it,
which is the check ``tests/test_ewald.py`` leans on hardest.  ``E_bg`` is the neutralising
background of a charged cell and vanishes for the neutral repeat units this package
builds -- it is implemented but has never been checked against anything, because nothing
here builds a charged cell.

**The surface term is the conditionally convergent part, and it is a physical choice.**
The reciprocal sum's ``k -> 0`` limit is what the truncated sum cannot evaluate, and its
value depends on the shape of the macroscopic sample and on what surrounds it -- not on
the crystal.  Two conventions:

``"tinfoil"`` (metallic, ``E_surf = 0``) -- **the default here.**  The infinite crystal is
    built up in a shape-independent way and the surface charge of the sample is screened
    by a conducting medium.  It is the bulk limit of a *short-circuited* or fully screened
    crystal, which is the condition under which a ferroelectric's spontaneous polarization
    is defined and measured (poled film between electrodes), and it is the convention
    every Berry-phase or Wannier polarization this package compares against is computed
    in.  It is also the only convention in which the energy is a property of the crystal
    rather than of the sample's outline.

``"vacuum"`` (spherical boundary, ``E_surf = 2 pi |M|^2 / 3 V``) -- the isolated sample in
    vacuum, including the depolarising field it makes for itself.  This is a real physical
    situation and it is not the bulk crystal: the term is a shape factor times ``P^2``,
    it penalises *every* polar cell, and a real unelectroded ferroelectric answers it by
    forming domains rather than by paying it.  Selectable because for a ferroelectric the
    choice is physics; not the default because a single polar unit cell with no domain
    structure and no screening is a slab that does not exist.

Which one is in force is carried on :class:`EwaldSpec` and reported by
:meth:`EwaldSpec.label`, so no energy from this module is quotable without it.

**Accuracy.**  ``alpha`` and ``kmax`` default to the standard accuracy-driven choice:
with ``p = sqrt(-ln(accuracy))``, ``alpha = p / rc`` (so ``erfc(alpha rc) ~ accuracy``)
and ``kmax = 2 alpha p`` (so ``exp(-kmax^2/4 alpha^2) = accuracy``).  Both are
overridable, and both *must* be overridable: the fact that the total does not move when
they are is the test that the implementation is right.

**Scope.**  NumPy and float64 only.  This is the deformable/direct path -- polish,
refinement, response -- not the batched screen; see :mod:`polyfind.lattice_table` for why
the screen cannot use it (the reciprocal sum is not pairwise, so it cannot be tabulated as
a chain-pair interaction).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.special import erfc

from .forcefield import COULOMB

BOUNDARIES = ("tinfoil", "vacuum")
_TWO_OVER_SQRT_PI = 2.0 / np.sqrt(np.pi)


@dataclass(frozen=True)
class EwaldSpec:
    """How the Ewald sum is split and what surrounds the crystal.

    ``boundary``
        ``"tinfoil"`` (default) or ``"vacuum"``; see the module docstring.  This is the
        one field that changes the *physics* rather than the numerics.
    ``accuracy``
        target relative truncation of both halves; sets ``alpha`` and ``kmax`` when they
        are not given.
    ``alpha``, ``rc``, ``kmax``
        explicit splitting parameter (1/A), real-space cutoff (A) and reciprocal cutoff
        (1/A).  ``rc`` defaults to the caller's own cutoff.  Setting ``alpha`` and
        checking the total does not move is the validation this module is built around.
    """

    boundary: str = "tinfoil"
    accuracy: float = 1e-8
    alpha: float | None = None
    rc: float | None = None
    kmax: float | None = None

    def __post_init__(self) -> None:
        if self.boundary not in BOUNDARIES:
            raise ValueError(f"unknown Ewald boundary {self.boundary!r} (expected one of {BOUNDARIES})")
        if not 0.0 < float(self.accuracy) < 1.0:
            raise ValueError(f"Ewald accuracy must be in (0, 1), got {self.accuracy!r}")
        for name in ("alpha", "rc", "kmax"):
            v = getattr(self, name)
            if v is not None and not float(v) > 0.0:
                raise ValueError(f"Ewald {name} must be positive, got {v!r}")

    def label(self) -> str:
        return (f"ewald[{self.boundary}, acc={self.accuracy:g}"
                + (f", alpha={self.alpha:g}" if self.alpha is not None else "")
                + (f", rc={self.rc:g}" if self.rc is not None else "")
                + (f", kmax={self.kmax:g}" if self.kmax is not None else "") + "]")


THOLE_FORMS = ("exp", "cubic")


@dataclass(frozen=True)
class Thole:
    r"""Thole damping of the dipole-dipole interaction at short range.

    Two point dipoles closer than about ``(4 alpha_i alpha_j)^(1/6)`` polarise each other
    without bound -- the polarization catastrophe -- and bonded neighbours are always that
    close.  Thole (Chem. Phys. **59**, 341 (1981)) replaces each point dipole by a smeared
    one, so the interaction tensor ``T = (3 r r^T / r^5 - I / r^3)`` becomes
    ``3 lambda5 r r^T / r^5 - lambda3 I / r^3`` with ``lambda -> 1`` beyond a few angstrom.
    The smearing is in the scaled distance ``u = r / (alpha_i alpha_j)^(1/6)``, so one
    dimensionless ``a`` serves every pair.

    ``form="exp"`` is the exponential density ``rho(u) = a^3/(8 pi) exp(-a u)``, giving
    ``lambda3 = 1 - (1 + a u + a^2 u^2 / 2) e^{-a u}`` and
    ``lambda5 = lambda3 - (a^3 u^3 / 6) e^{-a u}``; this is the form van Duijnen and Swart
    (J. Phys. Chem. A **102**, 2399 (1998)) refitted together with the atomic
    polarizabilities :mod:`polyfind.polarizability` ships, with ``a = 2.1304``.
    ``form="cubic"`` is ``rho(u) = 3a/(4 pi) exp(-a u^3)``, ``lambda3 = 1 - e^{-a u^3}``,
    ``lambda5 = 1 - (1 + a u^3) e^{-a u^3}``, the form AMOEBA uses with ``a = 0.39``
    (Ren and Ponder, J. Phys. Chem. B **107**, 5933 (2003)).

    :meth:`damping` returns the *undamped fractions* ``f3 = 1 - lambda3`` and
    ``f5 = 1 - lambda5`` and their ``u``-derivatives: the Ewald real-space sum adds
    ``-f3 / r^3`` and ``-3 f5 / r^5`` to its screened kernels, which is exactly the
    short-ranged difference between the damped and the bare interaction, so one cutoff
    serves both.
    """

    a: float = 2.1304
    form: str = "exp"

    def __post_init__(self) -> None:
        if self.form not in THOLE_FORMS:
            raise ValueError(f"unknown Thole form {self.form!r} (expected one of {THOLE_FORMS})")
        if not float(self.a) > 0.0:
            raise ValueError(f"Thole a must be positive, got {self.a!r}")

    def damping(self, u):
        """``(f3, f5, df3/du, df5/du)`` at scaled distance ``u`` (array)."""
        a = float(self.a)
        u = np.asarray(u, dtype=float)
        if self.form == "exp":
            au = a * u
            e = np.exp(-au)
            f3 = e * (1.0 + au + 0.5 * au * au)
            f5 = f3 + e * au * au * au / 6.0
            d3 = -0.5 * a ** 3 * u * u * e
            d5 = -(a ** 4) * u * u * u / 6.0 * e
        else:
            u3 = u * u * u
            e = np.exp(-a * u3)
            f3 = e
            f5 = (1.0 + a * u3) * e
            d3 = -3.0 * a * u * u * e
            d5 = -3.0 * a * a * u3 * u * u * e
        return f3, f5, d3, d5

    def label(self) -> str:
        return f"thole[{self.form}, a={self.a:g}]"


@dataclass
class EwaldTerms:
    """The five pieces of the sum, kept separate because that is what a check needs.

    With dipoles (see :meth:`Ewald.terms`) ``real``, ``recip``, ``self_`` and ``surface``
    include the charge-dipole and dipole-dipole parts; ``grad_dipoles`` is ``dE/dp_i``, i.e.
    minus the total field at ``i``; ``field`` is the field at every site from the *charges
    alone* (the permanent field an induced dipole responds to); and ``tensor`` is the
    ``(N, N, 3, 3)`` dipole interaction matrix ``T`` with ``field_i = sum_j T_ij p_j``,
    Ewald-summed, self-interaction with the images included on the diagonal blocks, and
    Thole-damped when a :class:`Thole` was given.
    """

    real: float
    recip: float
    self_: float
    background: float
    surface: float
    grad_coords: np.ndarray | None = None  # (N, 3) dE/dr_i at fixed lattice
    grad_lattice: np.ndarray | None = None  # (3, 3) dE/dh at fixed Cartesian r_i
    grad_charges: np.ndarray | None = None  # (N,) dE/dq_i
    grad_dipoles: np.ndarray | None = None  # (N, 3) dE/dp_i = -(field at i)
    field: np.ndarray | None = None  # (N, 3) field at i from the charges alone
    tensor: np.ndarray | None = None  # (N, N, 3, 3) dipole-dipole field matrix
    n_kvec: int = 0
    n_image: int = 0

    @property
    def total(self) -> float:
        return self.real + self.recip + self.self_ + self.background + self.surface


class Ewald:
    """Ewald sum of a periodic point-charge cell, with gradients.

    ``lat`` is always the (3, 3) matrix whose **rows** are the lattice vectors, matching
    :meth:`polyfind.pack.CrystalPacker._place`.  ``coords`` are Cartesian and are *not*
    wrapped into the cell: the real and reciprocal sums are invariant under wrapping any
    atom by a lattice vector, but the surface term is not, and the molecular branch (every
    chain kept whole, as :meth:`polyfind.pack.CrystalPacker.dipole` places it) is the one
    that makes ``M`` a property of the cell rather than of the wrapping.  For the same
    reason :attr:`EwaldTerms.grad_lattice` is the derivative at fixed *Cartesian*
    coordinates, which is the convention the packer's chain rule already uses for its
    truncated kernel.
    """

    def __init__(self, rc: float, spec: EwaldSpec | None = None, eps_r: float = 1.0,
                 coulomb: float = COULOMB, chunk_elems: int = 4_000_000):
        spec = EwaldSpec() if spec is None else spec
        self.spec = spec
        self.boundary = spec.boundary
        self.rc = float(rc if spec.rc is None else spec.rc)
        self.p = float(np.sqrt(-np.log(float(spec.accuracy))))
        self.alpha = float(self.p / self.rc if spec.alpha is None else spec.alpha)
        self.kmax = float(2.0 * self.alpha * self.p if spec.kmax is None else spec.kmax)
        self.eps_r = float(eps_r)
        self.pref = float(coulomb) / self.eps_r
        self.chunk_elems = int(chunk_elems)
        self._real_grid: dict = {}
        self._recip_grid: dict = {}

    def describe(self) -> str:
        return (f"Ewald(boundary={self.boundary}, alpha={self.alpha:.4f} 1/A, rc={self.rc:.3f} A, "
                f"kmax={self.kmax:.4f} 1/A, eps_r={self.eps_r:g})")

    # ------------------------------------------------------------------ grids
    @staticmethod
    def _grid(nmax, cache):
        key = tuple(int(v) for v in nmax)
        g = cache.get(key)
        if g is None:
            rngs = [np.arange(-k, k + 1, dtype=float) for k in key]
            g = np.stack(np.meshgrid(*rngs, indexing="ij"), axis=-1).reshape(-1, 3)
            cache[key] = g
        return g

    # ------------------------------------------------------------------ the sum
    def terms(self, coords, charges, lat, grad: bool = False, charge_grad: bool = False,
              dipoles=None, thole: Thole | None = None, thole_scale=None,
              dipole_grad: bool = False, field: bool = False, tensor: bool = False) -> EwaldTerms:
        r"""Every term of the sum for one configuration, optionally with gradients.

        With ``dipoles`` ``(N, 3)`` (e.A) every site carries a point dipole as well as a
        charge, and the sum is the Ewald sum of the charge-dipole system.  Writing
        ``B0 = erfc(a r)/r`` and ``B_{n+1} = [(2n+1) B_n + (2a^2)^{n+1}/(a sqrt(pi)) e^{-a^2 r^2}]/r^2``
        for the screened ``1/r``, ``1/r^3``, ``3/r^5`` kernels, the real-space pair term is

            q_i q_j B0 + (q_i p_j.d - q_j p_i.d) B1 + (p_i.p_j) B1 - (p_i.d)(p_j.d) B2

        with ``d = r_i - r_j + n``; the structure factor becomes
        ``S(k) = sum_j (q_j + i k.p_j) e^{i k.r_j}``; the dipole self energy is
        ``-(2 a^3 / 3 sqrt(pi)) sum |p_i|^2``; and the vacuum surface term uses
        ``M = sum q_i r_i + sum p_i``.  These are the standard expressions (Smith, CCP5
        Newsletter **46**, 18 (1998); Aguado and Madden, J. Chem. Phys. **119**, 7471
        (2003)), and the checks in ``tests/test_polarizable.py`` are the ones that matter:
        independence of the splitting parameter, the Lorentz field ``4 pi p / 3V`` of a
        simple cubic dipole lattice under tinfoil boundaries and zero under vacuum ones, and
        the Clausius-Mossotti dielectric constant of a cubic lattice of polarizable points.

        ``thole`` (a :class:`Thole`) with ``thole_scale`` ``(N, N)`` -- the per-pair
        ``(alpha_i alpha_j)^(1/6)`` -- damps the dipole-dipole part at short range; the
        charge-charge and charge-dipole parts are never damped.  ``dipole_grad`` asks for
        ``dE/dp_i``; ``field`` for the field at every site from the charges alone; ``tensor``
        for the dipole interaction matrix.  All three are what a self-consistent induced-dipole
        solve needs and none of them requires ``dipoles``: the field and the tensor are
        properties of the geometry, and ``dipole_grad`` with ``dipoles=None`` is minus the
        field.

        The charge-only path (``dipoles=None`` and none of the three flags) is the one every
        earlier number was measured with and its arithmetic is untouched, expression for
        expression.
        """
        X = np.asarray(coords, dtype=float).reshape(-1, 3)
        q = np.asarray(charges, dtype=float).reshape(-1)
        h = np.asarray(lat, dtype=float).reshape(3, 3)
        if X.shape[0] != q.shape[0]:
            raise ValueError(f"{X.shape[0]} coordinates but {q.shape[0]} charges")
        det = float(np.linalg.det(h))
        if abs(det) < 1e-12:
            raise ValueError("the lattice is singular: an Ewald sum needs three independent vectors")
        V = abs(det)
        hinv = np.linalg.inv(h)
        N = q.shape[0]
        p = None if dipoles is None else np.asarray(dipoles, dtype=float).reshape(N, 3)
        has_p = p is not None
        want_dip = has_p or field or tensor or dipole_grad
        if (thole is None) != (thole_scale is None):
            raise ValueError("Thole damping needs both `thole` and the (N, N) `thole_scale`")
        sth = None if thole is None else np.asarray(thole_scale, dtype=float).reshape(N, N)
        gX = np.zeros((N, 3)) if grad else None
        gq = np.zeros(N) if charge_grad else None
        gp = np.zeros((N, 3)) if dipole_grad else None
        Fq = np.zeros((N, 3)) if field else None
        Tm = np.zeros((N, N, 3, 3)) if tensor else None
        glat = np.zeros((3, 3)) if grad else None
        deps = np.zeros((3, 3)) if grad else None  # dE/d(strain), converted to dE/dh at the end
        I3 = np.eye(3)

        # --- real space ------------------------------------------------------------
        # Pair separations are taken in the minimum-image *fractional* frame and then
        # re-indexed: d = dX + (m - round(ds)) . h is the same set of images as the naive
        # sum, but the integer box needed to cover the cutoff is the tight one.  The
        # effective integers ``m - round(ds)`` are what the lattice derivative needs, since
        # ``round`` is locally constant and the atoms are held at fixed Cartesian positions.
        S = X @ hinv
        dS = S[:, None, :] - S[None, :, :]
        nint = np.round(dS)
        base = (dS - nint) @ h  # (N, N, 3)
        widths = V / np.linalg.norm(np.cross(np.roll(h, -1, axis=0), np.roll(h, -2, axis=0)), axis=1)
        nmax = np.ceil(self.rc / widths + 0.5).astype(int)
        m_all = self._grid(nmax, self._real_grid)
        e_real = 0.0
        gq_real = np.zeros(N) if charge_grad else None
        A1 = np.zeros((3, 3)) if grad else None
        Wij = np.zeros((N, N, 3)) if grad else None
        qq = q[:, None] * q[None, :]
        pp = None if not has_p else p @ p.T
        step = max(1, self.chunk_elems // max(1, N * N))
        for t in range(0, m_all.shape[0], step):
            m = m_all[t : t + step]
            D = base[None] + (m @ h)[:, None, None, :]  # (I, N, N, 3)
            r2 = np.einsum("mijc,mijc->mij", D, D)
            mask = (r2 < self.rc * self.rc) & (r2 > 1e-12)
            r2s = np.where(mask, r2, 1.0)
            r = np.sqrt(r2s)
            ar = self.alpha * r
            f = np.where(mask, erfc(ar) / r, 0.0)
            e_real += 0.5 * self.pref * float((qq[None] * f).sum())
            if charge_grad:
                gq_real = gq_real + self.pref * (f.sum(axis=0) @ q)
            if grad and not want_dip:
                # f'(r) = -erfc(a r)/r^2 - (2a/sqrt(pi)) exp(-a^2 r^2)/r
                fp = -(erfc(ar) / r2s + _TWO_OVER_SQRT_PI * self.alpha * np.exp(-ar * ar) / r)
                tt = np.where(mask, self.pref * qq[None] * fp / r, 0.0)
                W = tt[..., None] * D  # (I, N, N, 3);  dE/dd = W
                gX += W.sum(axis=(0, 2))
                A1 += m.T @ W.sum(axis=(1, 2))
                Wij += W.sum(axis=0)
            if want_dip:
                g = _TWO_OVER_SQRT_PI * self.alpha * np.exp(-ar * ar)
                B1 = np.where(mask, (f + g) / r2s, 0.0)
                B2 = np.where(mask, (3.0 * B1 + 2.0 * self.alpha ** 2 * g) / r2s, 0.0)
                C1, C2 = B1, B2  # the dipole-dipole kernels, Thole-damped when asked
                if thole is not None:
                    u = r / sth[None]
                    f3, f5, d3, d5 = thole.damping(u)
                    r3 = r2s * r
                    r5 = r3 * r2s
                    C1 = np.where(mask, B1 - f3 / r3, 0.0)
                    C2 = np.where(mask, B2 - 3.0 * f5 / r5, 0.0)
                if field:
                    Fq += self.pref * np.einsum("mij,j,mijc->ic", B1, q, D)
                if tensor:
                    Tm += self.pref * (np.einsum("mij,mija,mijb->ijab", C2, D, D)
                                       - C1.sum(axis=0)[:, :, None, None] * I3[None, None])
                if has_p:
                    pdj = np.einsum("mijc,jc->mij", D, p)  # p_j . d
                    pid = np.einsum("mijc,ic->mij", D, p)  # p_i . d
                    cross = q[:, None] * pdj - q[None, :] * pid  # q_i p_j.d - q_j p_i.d
                    e_real += 0.5 * self.pref * float((cross * B1).sum())
                    e_real += 0.5 * self.pref * float((pp[None] * C1 - pid * pdj * C2).sum())
                    if charge_grad:
                        gq_real = gq_real + self.pref * (pdj * B1).sum(axis=(0, 2))
                    if dipole_grad:
                        gp += self.pref * (np.einsum("mij,jc->ic", C1, p)
                                           - np.einsum("mij,mijc->ic", pdj * C2, D))
                    if grad:
                        B3 = np.where(mask, (5.0 * B2 + 4.0 * self.alpha ** 4 * g) / r2s, 0.0)
                        C1p, C2p = -B2, -B3  # (dC1/dr)/r and (dC2/dr)/r
                        if thole is not None:
                            C1p = np.where(mask, C1p + 3.0 * f3 / (r5) - d3 / (sth[None] * r2s * r2s), 0.0)
                            C2p = np.where(mask, C2p + 15.0 * f5 / (r5 * r2s)
                                           - 3.0 * d5 / (sth[None] * r5 * r), 0.0)
                        # dE/dd for the pair, both roles combined (the image set is closed
                        # under negation and every term is odd in d, so the 1/2 cancels)
                        scal = -qq[None] * B1 - cross * B2 + pp[None] * C1p - pid * pdj * C2p
                        p_i, p_j = p[None, :, None, :], p[None, None, :, :]  # (1, N, 1, 3), (1, 1, N, 3)
                        q_i, q_j = q[None, :, None, None], q[None, None, :, None]
                        W = self.pref * (scal[..., None] * D
                                         + B1[..., None] * (q_i * p_j - q_j * p_i)
                                         - C2[..., None] * (pdj[..., None] * p_i + pid[..., None] * p_j))
                        gX += W.sum(axis=(0, 2))
                        A1 += m.T @ W.sum(axis=(1, 2))
                        Wij += W.sum(axis=0)
                elif grad:
                    tt = -self.pref * qq[None] * B1
                    W = tt[..., None] * D
                    gX += W.sum(axis=(0, 2))
                    A1 += m.T @ W.sum(axis=(1, 2))
                    Wij += W.sum(axis=0)
                if dipole_grad:
                    gp += -self.pref * np.einsum("mij,j,mijc->ic", B1, q, D)
        if grad:
            glat = 0.5 * (A1 - np.einsum("ija,ijb->ab", nint, Wij))
            del A1, Wij
        if charge_grad:
            gq += gq_real

        # --- reciprocal space ------------------------------------------------------
        B = 2.0 * np.pi * hinv.T  # rows: reciprocal lattice vectors, b_i . a_j = 2 pi delta_ij
        nk = np.floor(self.kmax * np.linalg.norm(h, axis=1) / (2.0 * np.pi)).astype(int)
        mk = self._grid(nk, self._recip_grid)
        k = mk @ B
        k2 = np.einsum("kc,kc->k", k, k)
        keep = (k2 > 1e-12) & (k2 <= self.kmax * self.kmax)
        k, k2 = k[keep], k2[keep]
        Ak = np.exp(-k2 / (4.0 * self.alpha ** 2)) / k2
        kr = X @ k.T  # (N, K)
        cos, sin = np.cos(kr), np.sin(kr)
        ck, sk = q @ cos, q @ sin
        pre = 2.0 * np.pi / V * self.pref
        if has_p:
            kp = p @ k.T  # (N, K)  k . p_j
            # S(k) = sum_j (q_j + i k.p_j) e^{i k.r_j}
            ReS = ck - (kp * sin).sum(axis=0)
            ImS = sk + (kp * cos).sum(axis=0)
        else:
            ReS, ImS = ck, sk
        s2 = ReS * ReS + ImS * ImS
        e_recip = pre * float(Ak @ s2)
        if grad or charge_grad or dipole_grad:
            wgt = ImS[None, :] * cos - ReS[None, :] * sin  # (N, K)
            if has_p:
                wq = q[:, None] * wgt - kp * (ReS[None, :] * cos + ImS[None, :] * sin)
        if grad:
            if has_p:
                gX += 2.0 * pre * ((Ak[None, :] * wq) @ k)
            else:
                gX += 2.0 * pre * (q[:, None] * ((Ak[None, :] * wgt) @ k))
            # dE/dk for each k, then dE/d(strain) = -sum_k u_k (x) k  -  E_recip I
            u = 2.0 * pre * (Ak * (-1.0 / (4.0 * self.alpha ** 2) - 1.0 / k2) * s2 * k.T
                             + Ak * ((q[:, None] * X).T @ wgt))  # (3, K)
            if has_p:
                u = 2.0 * pre * (Ak * (-1.0 / (4.0 * self.alpha ** 2) - 1.0 / k2) * s2 * k.T
                                 + Ak * (X.T @ wq + p.T @ wgt))
            deps -= u @ k
            deps -= e_recip * np.eye(3)
        if charge_grad:
            gq += 2.0 * pre * ((Ak * ReS) @ cos.T + (Ak * ImS) @ sin.T)
        if dipole_grad:
            gp += 2.0 * pre * ((Ak[None, :] * wgt) @ k)
        if field:
            Fq -= 2.0 * pre * ((Ak[None, :] * (sk[None, :] * cos - ck[None, :] * sin)) @ k)
        if tensor:
            kk = Ak[:, None, None] * k[:, :, None] * k[:, None, :]  # (K, 3, 3)
            Tm -= 2.0 * pre * (np.einsum("jk,lk,kab->jlab", cos, cos, kk)
                               + np.einsum("jk,lk,kab->jlab", sin, sin, kk))

        # --- self, background, surface ---------------------------------------------
        e_self = -self.alpha / np.sqrt(np.pi) * self.pref * float(q @ q)
        Q = float(q.sum())
        e_bg = -np.pi / (2.0 * V * self.alpha ** 2) * self.pref * Q * Q
        Mq = q @ X
        M = Mq if not has_p else Mq + p.sum(axis=0)
        e_surf = 0.0
        if self.boundary == "vacuum":
            e_surf = 2.0 * np.pi / (3.0 * V) * self.pref * float(M @ M)
        if charge_grad:
            gq += -2.0 * self.alpha / np.sqrt(np.pi) * self.pref * q
            gq += -np.pi / (V * self.alpha ** 2) * self.pref * Q
            if self.boundary == "vacuum":
                gq += 4.0 * np.pi / (3.0 * V) * self.pref * (X @ M)
        if grad:
            deps -= (e_bg + e_surf) * np.eye(3)
            if self.boundary == "vacuum":
                gX += 4.0 * np.pi / (3.0 * V) * self.pref * q[:, None] * M[None, :]
            glat = glat + hinv.T @ deps
        dself = 2.0 * self.alpha ** 3 / (3.0 * np.sqrt(np.pi)) * self.pref
        if has_p:
            e_self += -dself * float((p * p).sum())
        if dipole_grad:
            if has_p:
                gp += -2.0 * dself * p
            if self.boundary == "vacuum":
                gp += 4.0 * np.pi / (3.0 * V) * self.pref * M[None, :]
        if field and self.boundary == "vacuum":
            Fq -= 4.0 * np.pi / (3.0 * V) * self.pref * Mq[None, :]
        if tensor:
            Tm[np.arange(N), np.arange(N)] += 2.0 * dself * I3
            if self.boundary == "vacuum":
                Tm -= 4.0 * np.pi / (3.0 * V) * self.pref * I3[None, None]

        return EwaldTerms(real=e_real, recip=e_recip, self_=e_self, background=e_bg, surface=e_surf,
                          grad_coords=gX, grad_lattice=glat, grad_charges=gq, grad_dipoles=gp,
                          field=Fq, tensor=Tm, n_kvec=int(k.shape[0]), n_image=int(m_all.shape[0]))

    def energy(self, coords, charges, lat) -> float:
        return self.terms(coords, charges, lat).total


# --------------------------------------------------------------- bonded exclusions
def exclusion_correction(coords, charges, c: float, scales, pref: float = COULOMB,
                         grad: bool = False, charge_grad: bool = False):
    r"""The full ``1/r`` an Ewald sum puts in that the bonded exclusions must take out.

    An Ewald sum evaluates every pair at full weight.  This package's potential does not:
    within a chain, pairs 1-2 and 1-3 apart are excluded and 1-4 pairs are scaled by
    ``scale14`` -- across the periodic boundary too, which is what ``scales`` (a
    ``(2K+1, n, n)`` stack indexed ``k = -K .. K``, the packer's
    :meth:`~polyfind.pack.CrystalPacker._scale_column`) records.  So

        E_correction = -(1/2) sum_{k,i,j} (1 - s_k[i,j]) q_i q_j / r_ijk

    with ``r_ijk = |r_i - r_j - k c z^|`` and the ``i = j, k = 0`` self pair left out (that
    one is the Ewald self energy).  Subtracting the *whole* ``1/r`` rather than only its
    ``erf`` part keeps this independent of the splitting parameter, which is what lets the
    alpha-independence check cover the correction as well as the sum.

    This is a property of one chain's own geometry -- a flip is an isometry, so both
    orientations of a two-chain cell give the same number -- so for a rigid chain it is a
    constant, exactly as :attr:`polyfind.pack.CrystalPacker.e_intra` is.

    Returns ``(E, dE/dcoords (n, 3), dE/dc, dE/dq (n,))`` for ONE chain; the trailing three
    are ``None`` unless asked for.
    """
    X = np.asarray(coords, dtype=float).reshape(-1, 3)
    q = np.asarray(charges, dtype=float).reshape(-1)
    Sk = np.asarray(scales, dtype=float)
    n = X.shape[0]
    K = (Sk.shape[0] - 1) // 2
    ks = np.arange(-K, K + 1, dtype=float)
    u = (1.0 - Sk) * (q[None, :, None] * q[None, None, :])  # (2K+1, n, n)
    D = X[:, None, :] - X[None, :, :]
    dz = D[None, :, :, 2] - ks[:, None, None] * float(c)
    r2 = (D[None, :, :, 0] ** 2 + D[None, :, :, 1] ** 2 + dz * dz)
    # keyed on the *scale*, not on ``u``: a pair whose charge product happens to be zero
    # still contributes to dE/dq, and keying on u would silently drop it
    live = (np.abs(1.0 - Sk) > 0.0) & (r2 > 1e-12)
    r2s = np.where(live, r2, 1.0)
    r = np.sqrt(r2s)
    e = -0.5 * pref * float(np.where(live, u / r, 0.0).sum())
    gX = gc = gq = None
    if grad:
        # phi(d) = -pref u / |d|;  dphi/dd = pref u d / r^3
        t = np.where(live, pref * u / (r2s * r), 0.0)
        fx, fy, fz = t * D[None, :, :, 0], t * D[None, :, :, 1], t * dz
        gX = 0.5 * np.stack([f.sum(axis=(0, 2)) - f.sum(axis=(0, 1)) for f in (fx, fy, fz)], axis=1)
        gc = -0.5 * float((fz * ks[:, None, None]).sum())
    if charge_grad:
        w = np.where(live, (1.0 - Sk) / r, 0.0)
        gq = -pref * (w.sum(axis=0) @ q)
    return e, gX, gc, gq


def charge_dipole_exclusion(coords, charges, c: float, scales, dipoles=None, pref: float = COULOMB,
                            grad: bool = False, charge_grad: bool = False):
    r"""The permanent field of bonded neighbours that an induced dipole must *not* see.

    The Ewald sum lets every charge polarise every site.  This package's potential excludes
    the 1-2 and 1-3 Coulomb pairs within a chain and scales the 1-4 ones by ``scale14``, and
    the induced dipoles follow the same rule -- the bonded charges' field is absorbed in the
    valence terms exactly as their energy is, and a C-F charge 1.35 A away would otherwise
    dominate every induced dipole on the carbon.  With ``w_k[i, j] = 1 - s_k[i, j]`` and
    ``r = r_i - r_j - k c z^`` (``scales`` as in :func:`exclusion_correction`), the correction is

        E = -(1/2) sum_{k,i,j} w_k[i,j] [ q_i (p_j . r) - q_j (p_i . r) ] / r^3

    for ONE chain in its *placed* coordinates (the induced dipoles of the two chains differ,
    so unlike :func:`exclusion_correction` this is not a per-chain constant; for a flipped
    chain pass ``scales[::-1]``, since the mirror ``z -> -z`` maps image ``k`` to ``-k``).
    The dipole-dipole interaction is not excluded -- every pair of induced dipoles
    interacts, Thole-damped, which is Thole's model -- and the permanent field is not
    Thole-damped either, which is also his model (the 1981 fit had no permanent charges).

    Returns ``(E, dE/dcoords (n, 3), dE/dc, dE/dq (n,), F (n, 3))``: ``F = -dE/dp`` is the
    field correction at every atom and is returned whatever ``dipoles`` is (it does not
    depend on them); the energy and the other gradients are ``0`` / ``None`` unless dipoles
    are given and asked for.
    """
    X = np.asarray(coords, dtype=float).reshape(-1, 3)
    q = np.asarray(charges, dtype=float).reshape(-1)
    Sk = np.asarray(scales, dtype=float)
    n = X.shape[0]
    K = (Sk.shape[0] - 1) // 2
    ks = np.arange(-K, K + 1, dtype=float)
    w = 1.0 - Sk  # (2K+1, n, n)
    D = X[:, None, :] - X[None, :, :]
    R = np.broadcast_to(D[None], (2 * K + 1, n, n, 3)).copy()
    R[..., 2] -= ks[:, None, None] * float(c)
    r2 = np.einsum("kijc,kijc->kij", R, R)
    live = (np.abs(w) > 0.0) & (r2 > 1e-12)
    r2s = np.where(live, r2, 1.0)
    inv3 = np.where(live, w / (r2s * np.sqrt(r2s)), 0.0)
    # F_a = -dE/dp_a, both roles of the sum written out (no symmetry of ``scales`` assumed)
    F = 0.5 * pref * (np.einsum("kia,i,kiac->ac", inv3, q, R) - np.einsum("kaj,j,kajc->ac", inv3, q, R))
    if dipoles is None:
        return 0.0, None, None, None, F
    p = np.asarray(dipoles, dtype=float).reshape(n, 3)
    pdj = np.einsum("kijc,jc->kij", R, p)
    pid = np.einsum("kijc,ic->kij", R, p)
    cross = q[None, :, None] * pdj - q[None, None, :] * pid
    e = -0.5 * pref * float((inv3 * cross).sum())
    gX = gc = gq = None
    if grad:
        inv5 = inv3 / r2s
        G = pref * (inv3[..., None] * (q[None, :, None, None] * p[None, None, :, :]
                                       - q[None, None, :, None] * p[None, :, None, :])
                    - 3.0 * (inv5 * cross)[..., None] * R)
        gX = -0.5 * (G.sum(axis=(0, 2)) - G.sum(axis=(0, 1)))
        gc = 0.5 * float((ks[:, None, None] * G[..., 2]).sum())
    if charge_grad:
        gq = -0.5 * pref * ((inv3 * pdj).sum(axis=(0, 2)) - (inv3 * pid).sum(axis=(0, 1)))
    return e, gX, gc, gq, F


# --------------------------------------------------------------------- validation
def madelung_nacl(spec: EwaldSpec | None = None, rc: float = 12.0, a: float = 1.0) -> float:
    """Madelung constant of the rock-salt lattice from this module's own sum.

    The eight-ion cubic cell of edge ``2a`` -- a simple cubic lattice of alternating unit
    charges at nearest-neighbour distance ``a``, which is what rock salt is -- evaluated
    with ``coulomb = 1`` and ``eps_r = 1``.  The Madelung constant is defined per *ion
    pair*: ``E = -M N_pairs q^2 / a`` with ``N_pairs = 4`` here, so ``M = -E a / 4``.

    Published value 1.7475645946331822.  This is the single best check on an Ewald
    implementation, because the lattice is the textbook conditionally convergent sum and
    the answer is known to more digits than double precision can hold -- no tolerance has
    to be invented for it.
    """
    ew = Ewald(rc, spec, eps_r=1.0, coulomb=1.0)
    idx = np.array([[i, j, k] for i in (0, 1) for j in (0, 1) for k in (0, 1)], dtype=float)
    coords = idx * a
    charges = np.array([1.0 if (i + j + k) % 2 == 0 else -1.0 for i, j, k in idx.astype(int)])
    lat = np.eye(3) * (2.0 * a)
    return -ew.energy(coords, charges, lat) * a / 4.0


MADELUNG_NACL = 1.7475645946331822  # rock salt, per ion pair (published)
MADELUNG_CSCL = 1.762674773  # caesium chloride, per ion pair (published)
MADELUNG_ZINCBLENDE = 1.6380550533  # zinc blende, per ion pair (published)
