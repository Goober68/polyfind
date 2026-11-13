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

**Scope.**  Float64 only.  This is the deformable/direct path -- polish, refinement,
response -- not the batched screen; see :mod:`polyfind.lattice_table` for why the screen
cannot use it (the reciprocal sum is not pairwise, so it cannot be tabulated as a
chain-pair interaction).

**Array module.**  The sum is written once against an array module ``xp``.  With NumPy,
the default, it is the arithmetic every stored number was measured with, operation for
operation.  With CuPy (``POLYFIND_EWALD_DEVICE=cuda``, or ``Ewald(..., xp=cupy)``) the
arrays that grow with the cell -- the pair separations, the neighbour list, the ``(N, K)``
phases and the ``(N, N, 3, 3)`` tensor -- are built and kept on the device, the 3 x 3
lattice algebra and the list of k vectors stay on the host, and what comes back is on the
host again except the tensor when the caller asks to keep it where it is (the
induced-dipole solve, :meth:`polyfind.pack.CrystalPacker._polarize`).  A cell of a few
dozen sites is a millisecond on the host and slower on a device, so cells below
``Ewald.device_min_sites`` stay on the host whatever the module.  A device sums in another
order (its matrix products, and the scatter-adds, whose order is not even fixed from run to
run), so its numbers agree with the host's to rounding, not bit for bit.
"""
from __future__ import annotations

import importlib
from dataclasses import dataclass

import numpy as np

from . import backend as bk
from .forcefield import COULOMB

BOUNDARIES = ("tinfoil", "vacuum")
_TWO_OVER_SQRT_PI = float(2.0 / np.sqrt(np.pi))  # a plain float: the same number on every array module


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
    # Charge penetration: the permanent charges interact as ``1/r - (1 + b r / 2) exp(-b r) / r``
    # (a point charge against a Slater-type density of exponent ``b``, 1/A) instead of ``1/r``,
    # so two charges that overlap at 2.5-3 A attract less than point charges do.  ``None`` is
    # point charges.  Applied to the charge-charge term of the real-space sum (it is short
    # ranged) and to the bonded-exclusion correction, so excluded pairs stay excluded; the
    # charge-dipole and dipole-dipole terms are left as they are (Thole damps the latter).
    charge_damping: float | None = None

    def __post_init__(self) -> None:
        if self.boundary not in BOUNDARIES:
            raise ValueError(f"unknown Ewald boundary {self.boundary!r} (expected one of {BOUNDARIES})")
        if self.charge_damping is not None and not float(self.charge_damping) > 0.0:
            raise ValueError(f"Ewald charge_damping must be positive (1/A), got {self.charge_damping!r}")
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
                + (f", kmax={self.kmax:g}" if self.kmax is not None else "")
                + (f", damping={self.charge_damping:g}" if self.charge_damping is not None else "") + "]")


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

    def damping(self, u, xp=np):
        """``(f3, f5, df3/du, df5/du)`` at scaled distance ``u`` (an array of the module ``xp``)."""
        a = float(self.a)
        u = xp.asarray(u, dtype=float)
        if self.form == "exp":
            au = a * u
            e = xp.exp(-au)
            f3 = e * (1.0 + au + 0.5 * au * au)
            f5 = f3 + e * au * au * au / 6.0
            d3 = -0.5 * a ** 3 * u * u * e
            d5 = -(a ** 4) * u * u * u / 6.0 * e
        else:
            u3 = u * u * u
            e = xp.exp(-a * u3)
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

    Every array is a NumPy array on the host, whatever array module the sum ran on, with one
    exception: ``tensor`` is left on the device when :meth:`Ewald.terms` was called with
    ``device_tensor=True``.
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
                 coulomb: float = COULOMB, chunk_elems: int = 4_000_000, xp=None):
        spec = EwaldSpec() if spec is None else spec
        # The array module of the sums: NumPy unless asked otherwise (see the module docstring
        # and polyfind.backend.ewald_backend).  Cells of fewer than ``device_min_sites`` sites
        # are summed on the host whatever the module.
        self.xp = bk.ewald_backend() if xp is None else xp
        self.device_min_sites = 256
        self.spec = spec
        self.boundary = spec.boundary
        self.rc = float(rc if spec.rc is None else spec.rc)
        self.p = float(np.sqrt(-np.log(float(spec.accuracy))))
        self.alpha = float(self.p / self.rc if spec.alpha is None else spec.alpha)
        self.kmax = float(2.0 * self.alpha * self.p if spec.kmax is None else spec.kmax)
        self.eps_r = float(eps_r)
        self.pref = float(coulomb) / self.eps_r
        self.charge_damping = None if spec.charge_damping is None else float(spec.charge_damping)
        self.chunk_elems = int(chunk_elems)
        self._real_grid: dict = {}
        self._recip_grid: dict = {}
        self._phase_cache: tuple | None = None
        self.prune_images = True  # False builds every image of the grid, the ones out of reach too (a test's switch)

    def describe(self) -> str:
        return (f"Ewald(boundary={self.boundary}, alpha={self.alpha:.4f} 1/A, rc={self.rc:.3f} A, "
                f"kmax={self.kmax:.4f} 1/A, eps_r={self.eps_r:g})")

    def array_module(self, n_sites: int):
        """The array module a cell of ``n_sites`` sites is summed with (NumPy for a small one)."""
        return np if self.xp is np or int(n_sites) < self.device_min_sites else self.xp

    def release(self) -> None:
        """Drop the phases kept from the last call (on a device they are the memory it still holds)."""
        self._phase_cache = None

    def __getstate__(self):
        # a module does not pickle and the kept phases are not worth sending: the module goes by
        # name and the phases are rebuilt by whoever needs them
        state = dict(self.__dict__)
        state["_phase_cache"] = None
        state["xp"] = None if self.xp is np else getattr(self.xp, "__name__", None)
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)
        self.xp = np if state.get("xp") is None else importlib.import_module(state["xp"])

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

    def _phases(self, X, k, xp=np, on_device=None):
        """``cos(k.r)`` and ``sin(k.r)``, ``(N, K)`` each, kept from the last call.

        An induced-dipole evaluation calls :meth:`terms` twice on one geometry (the field and the
        tensor for the solve, then the gradient at the solved dipoles), and on a large cell these
        two arrays are a tenth of each call.  The arrays of the last geometry are returned when
        the coordinates and the vectors are the same numbers, element for element; they are
        read-only so that no caller can change what the next call gets.  Small problems are not
        kept: there is nothing to save and nothing worth holding.

        ``X`` and ``k`` are host arrays (they are what the kept arrays are recognised by); the
        phases are arrays of ``xp`` and are kept for that module only.  ``on_device`` is the pair
        ``(X, k)`` as arrays of ``xp`` when the caller has them already.  A device array has no
        read-only flag, so there it is this module's own care that nothing writes into them.
        """
        hit = self._phase_cache
        if (hit is not None and hit[4] is xp and hit[0].shape == X.shape and hit[1].shape == k.shape
                and np.array_equal(hit[0], X) and np.array_equal(hit[1], k)):
            return hit[2], hit[3]
        self._phase_cache = hit = None  # the last geometry's arrays go before this one's are built
        Xd, kd = (xp.asarray(X), xp.asarray(k)) if on_device is None else on_device
        kr = Xd @ kd.T  # (N, K)
        cos, sin = xp.cos(kr), xp.sin(kr)
        if kr.size >= 1_000_000:
            if xp is np:
                cos.setflags(write=False)
                sin.setflags(write=False)
            self._phase_cache = (X.copy(), k.copy(), cos, sin, xp)
        return cos, sin

    # ------------------------------------------------------------------ the sum
    def terms(self, coords, charges, lat, grad: bool = False, charge_grad: bool = False,
              dipoles=None, thole: Thole | None = None, thole_scale=None,
              dipole_grad: bool = False, field: bool = False, tensor: bool = False,
              device_tensor: bool = False) -> EwaldTerms:
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

        ``coords``, ``charges``, ``lat`` and ``dipoles`` are host arrays.  On a device (see the
        module docstring) ``thole_scale`` may be handed over as a device array, which saves
        sending an ``(N, N)`` array at every call, and ``device_tensor=True`` leaves the tensor
        on the device for a solve there; every other result comes back as a host array.
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
        # Below, X, q, h, hinv, p and everything made of the k vectors alone are host arrays;
        # Xd, qd, hd, hinvd, pd, kd and every array that grows with the cell are arrays of ``xp``.
        # With NumPy the two are the same objects and each line is the expression it always was.
        xp = self.array_module(N)
        host = xp is np
        to_host = bk.to_numpy
        erfc = bk.erfc_of(xp)
        Xd, qd, hd, hinvd = xp.asarray(X), xp.asarray(q), xp.asarray(h), xp.asarray(hinv)
        pd = None if p is None else xp.asarray(p)
        sth = None if thole is None else xp.asarray(thole_scale, dtype=float).reshape(N, N)
        gX = xp.zeros((N, 3)) if grad else None
        gq = xp.zeros(N) if charge_grad else None
        gp = xp.zeros((N, 3)) if dipole_grad else None
        Fq = xp.zeros((N, 3)) if field else None
        Tm = xp.zeros((N, N, 3, 3)) if tensor else None
        glat = np.zeros((3, 3)) if grad else None
        deps = np.zeros((3, 3)) if grad else None  # dE/d(strain), converted to dE/dh at the end
        I3 = xp.eye(3)

        # --- real space ------------------------------------------------------------
        # Pair separations are taken in the minimum-image *fractional* frame and then
        # re-indexed: d = dX + (m - round(ds)) . h is the same set of images as the naive
        # sum, but the integer box needed to cover the cutoff is the tight one.  The
        # effective integers ``m - round(ds)`` are what the lattice derivative needs, since
        # ``round`` is locally constant and the atoms are held at fixed Cartesian positions.
        S = Xd @ hinvd
        dS = S[:, None, :] - S[None, :, :]
        nint = xp.round(dS)
        if host:
            base = (dS - nint) @ h  # (N, N, 3)
        else:  # one product of an (N^2, 3) array, not N stacked ones
            base = ((dS - nint).reshape(N * N, 3) @ hd).reshape(N, N, 3)
        widths = V / np.linalg.norm(np.cross(np.roll(h, -1, axis=0), np.roll(h, -2, axis=0)), axis=1)
        nmax = np.ceil(self.rc / widths + 0.5).astype(int)
        m_all = self._grid(nmax, self._real_grid)
        e_real = 0.0
        gq_real = xp.zeros(N) if charge_grad else None
        A1 = xp.zeros((3, 3)) if grad else None
        A2 = xp.zeros((3, 3)) if grad else None  # sum_k nint_k (x) W_k: the image-shift part of dE/dh
        qq = qd[:, None] * qd[None, :]
        pp = None if not has_p else pd @ pd.T
        rc2 = self.rc * self.rc
        step = max(1, self.chunk_elems // max(1, N * N))

        # The kernels are evaluated on the neighbour list only.  Every pair-image combination
        # within the cutoff is found from the full (I, N, N) distance array, which is one
        # cheap pass; everything after that runs on the K entries inside the sphere, which on
        # a 120-atom cell is about 2% of them (measured 2026-11-05: the dense form spent 83% of
        # a response calculation on masked-out arithmetic).  Sums over images and sites become
        # scatter-adds (bincount) and the lattice derivative's two moments are accumulated as
        # (3, 3) matrices, so no (N, N, 3) intermediate survives the loop.
        def scatter(idx, w):
            return xp.bincount(idx, weights=w, minlength=N)

        def scatter_vec(idx, W):
            return xp.stack([xp.bincount(idx, weights=W[:, c], minlength=N) for c in range(3)], axis=1)

        # An image whose shift alone puts every pair beyond the cutoff is not built: in the
        # minimum-image frame a pair's fractional separation along axis i is within 1/2 of the
        # image integer, so its distance is at least (|m_i| - 1/2) times the cell's width along
        # that axis.  In a cell wider than twice the cutoff that leaves the one image m = 0 of
        # the 27 of the grid (a sixteen-chain cell: 0.65 s of 0.68 spent on images that held
        # nothing).  The images that are kept are visited in the grid's order and in its chunks,
        # so every sum below runs over the same entries in the same order as before.
        can_hold = np.all((np.abs(m_all) - 0.5) * widths[None, :] < self.rc * (1.0 + 1e-9), axis=1)
        if not self.prune_images:
            can_hold[:] = True
        for t in range(0, m_all.shape[0], step):
            m = m_all[t : t + step]
            m = m[can_hold[t : t + step]]
            if m.shape[0] == 0:
                continue
            md = xp.asarray(m)
            D = base[None] + (md @ hd)[:, None, None, :]  # (I, N, N, 3)
            # (a device's einsum would run this as one tiny matrix product per pair)
            r2 = xp.einsum("mijc,mijc->mij", D, D) if host else (D * D).sum(axis=-1)
            mi, ii, jj = xp.nonzero((r2 < rc2) & (r2 > 1e-12))
            if mi.size == 0:
                continue
            Dk = D[mi, ii, jj]  # (K, 3)
            r2k = r2[mi, ii, jj]
            del D, r2
            rk = xp.sqrt(r2k)
            ark = self.alpha * rk
            erfck = erfc(ark)
            fk = erfck / rk
            qqk = qq[ii, jj]
            e_real += 0.5 * self.pref * float(qqk @ fk)
            if charge_grad:
                gq_real += self.pref * scatter(ii, fk * qd[jj])
            if grad:
                mk = md[mi]  # (K, 3) image integers
                nk = nint[ii, jj]  # (K, 3) the round() shifts of the minimum-image frame
            if grad and not want_dip:
                # f'(r) = -erfc(a r)/r^2 - (2a/sqrt(pi)) exp(-a^2 r^2)/r
                fp = -(erfck / r2k + _TWO_OVER_SQRT_PI * self.alpha * xp.exp(-ark * ark) / rk)
                W = (self.pref * qqk * fp / rk)[:, None] * Dk  # dE/dd per pair-image
                gX += scatter_vec(ii, W)
                A1 += mk.T @ W
                A2 += nk.T @ W
            if self.charge_damping is not None:
                # charge penetration: -q_i q_j g(r), g = (1 + b r/2) exp(-b r) / r, on the same list
                # (the bonded exclusions take their share back out in exclusion_correction)
                bq = self.charge_damping
                ebr = xp.exp(-bq * rk)
                gd = (1.0 + 0.5 * bq * rk) * ebr / rk
                e_real -= 0.5 * self.pref * float(qqk @ gd)
                if charge_grad:
                    gq_real -= self.pref * scatter(ii, gd * qd[jj])
                if grad:
                    # g'(r) = -exp(-b r) [ (b/2 + b^2 r/2) / r + (1 + b r/2) / r^2 ]
                    gdp = -ebr * ((0.5 * bq + 0.5 * bq * bq * rk) / rk + (1.0 + 0.5 * bq * rk) / r2k)
                    Wd = (-self.pref * qqk * gdp / rk)[:, None] * Dk  # the term's f' is -g'
                    gX += scatter_vec(ii, Wd)
                    A1 += mk.T @ Wd
                    A2 += nk.T @ Wd
            if want_dip:
                gk = _TWO_OVER_SQRT_PI * self.alpha * xp.exp(-ark * ark)
                B1 = (fk + gk) / r2k
                B2 = (3.0 * B1 + 2.0 * self.alpha ** 2 * gk) / r2k
                C1, C2 = B1, B2  # the dipole-dipole kernels, Thole-damped when asked
                if thole is not None:
                    sk = sth[ii, jj]
                    f3, f5, d3, d5 = thole.damping(rk / sk, xp)
                    r3 = r2k * rk
                    r5 = r3 * r2k
                    C1 = B1 - f3 / r3
                    C2 = B2 - 3.0 * f5 / r5
                if field:
                    Fq += self.pref * scatter_vec(ii, (B1 * qd[jj])[:, None] * Dk)
                if tensor:
                    pair = ii * N + jj
                    for a in range(3):
                        for b in range(a, 3):
                            tab = xp.bincount(pair, weights=C2 * Dk[:, a] * Dk[:, b], minlength=N * N).reshape(N, N)
                            Tm[:, :, a, b] += self.pref * tab
                            if b != a:
                                Tm[:, :, b, a] += self.pref * tab
                    c1 = xp.bincount(pair, weights=C1, minlength=N * N).reshape(N, N)
                    Tm -= self.pref * c1[:, :, None, None] * I3[None, None]
                if has_p:
                    pi_, pj_ = pd[ii], pd[jj]
                    pdj = (Dk * pj_).sum(axis=1)  # p_j . d
                    pid = (Dk * pi_).sum(axis=1)  # p_i . d
                    cross = qd[ii] * pdj - qd[jj] * pid  # q_i p_j.d - q_j p_i.d
                    ppk = pp[ii, jj]
                    e_real += 0.5 * self.pref * float(cross @ B1)
                    e_real += 0.5 * self.pref * float((ppk * C1 - pid * pdj * C2).sum())
                    if charge_grad:
                        gq_real += self.pref * scatter(ii, pdj * B1)
                    if dipole_grad:
                        gp += self.pref * (scatter_vec(ii, C1[:, None] * pj_) - scatter_vec(ii, (pdj * C2)[:, None] * Dk))
                    if grad:
                        B3 = (5.0 * B2 + 4.0 * self.alpha ** 4 * gk) / r2k
                        C1p, C2p = -B2, -B3  # (dC1/dr)/r and (dC2/dr)/r
                        if thole is not None:
                            C1p = C1p + 3.0 * f3 / r5 - d3 / (sk * r2k * r2k)
                            C2p = C2p + 15.0 * f5 / (r5 * r2k) - 3.0 * d5 / (sk * r5 * rk)
                        # dE/dd for the pair, both roles combined (the image set is closed
                        # under negation and every term is odd in d, so the 1/2 cancels)
                        scal = -qqk * B1 - cross * B2 + ppk * C1p - pid * pdj * C2p
                        W = self.pref * (scal[:, None] * Dk
                                         + B1[:, None] * (qd[ii][:, None] * pj_ - qd[jj][:, None] * pi_)
                                         - C2[:, None] * (pdj[:, None] * pi_ + pid[:, None] * pj_))
                        gX += scatter_vec(ii, W)
                        A1 += mk.T @ W
                        A2 += nk.T @ W
                elif grad:
                    W = (-self.pref * qqk * B1)[:, None] * Dk
                    gX += scatter_vec(ii, W)
                    A1 += mk.T @ W
                    A2 += nk.T @ W
                if dipole_grad:
                    gp += -self.pref * scatter_vec(ii, (B1 * qd[jj])[:, None] * Dk)
        del base, dS, nint, S
        if grad:
            glat = to_host(0.5 * (A1 - A2))
            del A1, A2
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
        mkept = mk[keep] if tensor else None
        Ak = np.exp(-k2 / (4.0 * self.alpha ** 2)) / k2
        kd, k2d, Akd = xp.asarray(k), xp.asarray(k2), xp.asarray(Ak)
        cos, sin = self._phases(X, k, xp, (Xd, kd))  # (N, K) each
        ck, sk = qd @ cos, qd @ sin
        pre = 2.0 * np.pi / V * self.pref
        if has_p:
            kp = pd @ kd.T  # (N, K)  k . p_j
            # S(k) = sum_j (q_j + i k.p_j) e^{i k.r_j}
            ReS = ck - (kp * sin).sum(axis=0)
            ImS = sk + (kp * cos).sum(axis=0)
        else:
            ReS, ImS = ck, sk
        s2 = ReS * ReS + ImS * ImS
        e_recip = pre * float(Akd @ s2)
        if grad or charge_grad or dipole_grad:
            wgt = ImS[None, :] * cos - ReS[None, :] * sin  # (N, K)
            if has_p:
                wq = qd[:, None] * wgt - kp * (ReS[None, :] * cos + ImS[None, :] * sin)
        if grad:
            if has_p:
                gX += 2.0 * pre * ((Akd[None, :] * wq) @ kd)
            else:
                gX += 2.0 * pre * (qd[:, None] * ((Akd[None, :] * wgt) @ kd))
            # dE/dk for each k, then dE/d(strain) = -sum_k u_k (x) k  -  E_recip I
            if has_p:
                u = 2.0 * pre * (Akd * (-1.0 / (4.0 * self.alpha ** 2) - 1.0 / k2d) * s2 * kd.T
                                 + Akd * (Xd.T @ wq + pd.T @ wgt))  # (3, K)
            else:
                u = 2.0 * pre * (Akd * (-1.0 / (4.0 * self.alpha ** 2) - 1.0 / k2d) * s2 * kd.T
                                 + Akd * ((qd[:, None] * Xd).T @ wgt))
            deps -= to_host(u @ kd)
            deps -= e_recip * np.eye(3)
        if charge_grad:
            gq += 2.0 * pre * ((Akd * ReS) @ cos.T + (Akd * ImS) @ sin.T)
        if dipole_grad:
            gp += 2.0 * pre * ((Akd[None, :] * wgt) @ kd)
        if field:
            Fq -= 2.0 * pre * ((Akd[None, :] * (sk[None, :] * cos - ck[None, :] * sin)) @ kd)
        if tensor:
            # T_jl,ab -= 2 pre sum_k A_k k_a k_b cos(k.(r_j - r_l)).  The summand is even in k, so
            # the sum runs over one of each pair (k, -k) at twice the weight; and with
            # M = [cos | sin] on that half, S(v) = sum_k 2 A_k (v.k)^2 cos(k.(r_j - r_l)) is the one
            # symmetric product (M s)(M s)^T with s = sqrt(2 A_k) |v.k|, which the BLAS forms at
            # half the cost of a general one.  The diagonal components are S(e_a) and the others
            # follow from k_a k_b = [(k_a + k_b)^2 - k_a^2 - k_b^2] / 2: six symmetric products of
            # width K where there were twelve general ones (cos and sin apart, k and -k both), a
            # quarter of the arithmetic.  On a sixteen-chain cell (960 sites, 14 482 vectors) this
            # was half of the whole energy evaluation (measured 2026-11-12: 4.8 s of 10.2).
            # (On a device the same line is a general product: CuPy 14 does not hand ``a @ a.T`` to
            # the symmetric update.  docs/PERFORMANCE_REVIEW.md 17.3.)
            half = (mkept[:, 0] > 0) | ((mkept[:, 0] == 0) & ((mkept[:, 1] > 0) | ((mkept[:, 1] == 0) & (mkept[:, 2] > 0))))
            kh = k[half]
            ih = xp.asarray(np.flatnonzero(half))  # the same columns as the mask, as positions
            Mh = xp.concatenate([cos[:, ih], sin[:, ih]], axis=1)  # (N, K)
            Mw = xp.empty_like(Mh)
            root = np.sqrt(2.0 * Ak[half])
            sym = {}
            for v in ((1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (1, 0, 1), (0, 1, 1)):
                sv = root * (kh @ np.array(v, dtype=float))
                xp.multiply(Mh, xp.asarray(np.concatenate([sv, sv]))[None, :], out=Mw)
                sym[v] = Mw @ Mw.T
            del Mh, Mw
            for a in range(3):
                for b in range(a, 3):
                    ea, eb = tuple(int(c == a) for c in range(3)), tuple(int(c == b) for c in range(3))
                    if a == b:
                        tab = sym[ea]
                    else:
                        tab = 0.5 * (sym[tuple(x + y for x, y in zip(ea, eb))] - sym[ea] - sym[eb])
                    Tm[:, :, a, b] -= 2.0 * pre * tab
                    if b != a:
                        Tm[:, :, b, a] -= 2.0 * pre * tab

        # --- self, background, surface ---------------------------------------------
        # From here on the per-site results are host arrays (on NumPy they already were, and
        # nothing is copied); only the tensor stays an array of ``xp``.
        del cos, sin
        gX, gq, gp, Fq = (None if arr is None else to_host(arr) for arr in (gX, gq, gp, Fq))
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
        dself = float(2.0 * self.alpha ** 3 / (3.0 * np.sqrt(np.pi)) * self.pref)
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
            diag = xp.arange(N)
            Tm[diag, diag] += 2.0 * dself * I3
            if self.boundary == "vacuum":
                Tm -= 4.0 * np.pi / (3.0 * V) * self.pref * I3[None, None]
            if not device_tensor:
                Tm = to_host(Tm)

        return EwaldTerms(real=e_real, recip=e_recip, self_=e_self, background=e_bg, surface=e_surf,
                          grad_coords=gX, grad_lattice=glat, grad_charges=gq, grad_dipoles=gp,
                          field=Fq, tensor=Tm, n_kvec=int(k.shape[0]), n_image=int(m_all.shape[0]))

    def energy(self, coords, charges, lat) -> float:
        return self.terms(coords, charges, lat).total


# --------------------------------------------------------------- bonded exclusions
def exclusion_correction(coords, charges, c: float, scales, pref: float = COULOMB,
                         grad: bool = False, charge_grad: bool = False, charge_damping: float | None = None):
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
    # h(r) is the interaction the Ewald sum put in for the pair: 1/r for point charges, and with
    # charge penetration (EwaldSpec.charge_damping) 1/r - (1 + b r/2) exp(-b r)/r
    h = 1.0 / r
    hp_over_r = -1.0 / (r2s * r)  # h'(r) / r
    if charge_damping is not None:
        bq = float(charge_damping)
        ebr = np.exp(-bq * r)
        h = h - (1.0 + 0.5 * bq * r) * ebr / r
        hp_over_r = hp_over_r + ebr * ((0.5 * bq + 0.5 * bq * bq * r) / r + (1.0 + 0.5 * bq * r) / r2s) / r
    e = -0.5 * pref * float(np.where(live, u * h, 0.0).sum())
    gX = gc = gq = None
    if grad:
        # phi(d) = -pref u h(|d|);  dphi/dd = -pref u h'(r) d / r
        t = np.where(live, -pref * u * hp_over_r, 0.0)
        fx, fy, fz = t * D[None, :, :, 0], t * D[None, :, :, 1], t * dz
        gX = 0.5 * np.stack([f.sum(axis=(0, 2)) - f.sum(axis=(0, 1)) for f in (fx, fy, fz)], axis=1)
        gc = -0.5 * float((fz * ks[:, None, None]).sum())
    if charge_grad:
        w = np.where(live, (1.0 - Sk) * h, 0.0)
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


# ------------------------------------------------ bonded exclusions, all chains of a cell at once
@dataclass(frozen=True)
class ExclusionPairs:
    """The bonded pairs of every chain of a cell as flat lists.

    One entry per ``(image k, site i, site j)`` of a chain whose scale differs from 1, for every
    chain: ``i`` and ``j`` index the sites of the *cell*, ``k`` is the image of ``j`` along z in
    the *placed* frame (a flipped chain's own image ``k`` is the placed ``-k``), and ``w`` is
    ``1 - scale``.  The dense routines above evaluate ``(2K + 1) n^2`` terms per chain of which a
    few hundred are not zero; these lists are those few hundred, for all chains in one array, so
    the corrections of a cell of any number of chains are a fixed number of array operations.
    """

    i: np.ndarray
    j: np.ndarray
    k: np.ndarray
    w: np.ndarray
    n_sites: int

    def __len__(self) -> int:
        return int(self.i.size)


def exclusion_pairs(scales, flipped, select=None) -> ExclusionPairs:
    """:class:`ExclusionPairs` of a cell whose chains all carry the stack ``scales`` (``(2K+1, n, n)``).

    ``flipped`` has one flag per chain (a mirrored chain reads the stack reversed, as
    :func:`charge_dipole_exclusion` is told to with ``scales[::-1]``); ``select`` (``(n,)`` bool)
    keeps the pairs both of whose sites are selected -- the electrostatically active sites, since
    a site with no charge and no dipole contributes nothing to either correction.
    """
    Sk = np.asarray(scales, dtype=float)
    n = Sk.shape[1]
    K = (Sk.shape[0] - 1) // 2
    w = 1.0 - Sk
    live = np.abs(w) > 0.0
    live[K, np.arange(n), np.arange(n)] = False  # the site with itself: the Ewald self term
    if select is not None:
        sel = np.asarray(select, dtype=bool)
        live &= sel[None, :, None] & sel[None, None, :]
    kk, ii, jj = np.nonzero(live)
    ww = w[kk, ii, jj]
    kz = (kk - K).astype(float)
    fl = np.asarray(flipped, dtype=bool).ravel()
    base = (np.arange(fl.size) * n)[:, None]
    sign = np.where(fl, -1.0, 1.0)[:, None]
    return ExclusionPairs(i=(base + ii[None, :]).ravel(), j=(base + jj[None, :]).ravel(),
                          k=(sign * kz[None, :]).ravel(), w=np.tile(ww, fl.size), n_sites=int(fl.size * n))


def _scatter(idx, W, n: int) -> np.ndarray:
    """``out[idx[e]] += W[e]`` for vectors: (L,) indices, (L, 3) values -> (n, 3)."""
    return np.stack([np.bincount(idx, weights=W[:, c], minlength=n) for c in range(3)], axis=1)


def exclusion_correction_pairs(coords, charges, c: float, pairs: ExclusionPairs, pref: float = COULOMB,
                               charge_damping: float | None = None):
    """:func:`exclusion_correction` of every chain of a cell at once: ``(E, dE/dcoords (N, 3), dE/dc)``.

    ``coords`` and ``charges`` are the cell's, in the placed frame; the sum of the per-chain
    corrections, which do not depend on the frame."""
    X = np.asarray(coords, dtype=float).reshape(-1, 3)
    q = np.asarray(charges, dtype=float).reshape(-1)
    d = X[pairs.i] - X[pairs.j]
    d[:, 2] -= pairs.k * float(c)
    r2 = np.einsum("ec,ec->e", d, d)
    r = np.sqrt(r2)
    u = pairs.w * q[pairs.i] * q[pairs.j]
    h = 1.0 / r
    hp_over_r = -1.0 / (r2 * r)
    if charge_damping is not None:
        bq = float(charge_damping)
        ebr = np.exp(-bq * r)
        h = h - (1.0 + 0.5 * bq * r) * ebr / r
        hp_over_r = hp_over_r + ebr * ((0.5 * bq + 0.5 * bq * bq * r) / r + (1.0 + 0.5 * bq * r) / r2) / r
    e = -0.5 * pref * float((u * h).sum())
    f = (-pref * u * hp_over_r)[:, None] * d
    N = X.shape[0]
    gX = 0.5 * (_scatter(pairs.i, f, N) - _scatter(pairs.j, f, N))
    gc = -0.5 * float((f[:, 2] * pairs.k).sum())
    return e, gX, gc


def charge_dipole_exclusion_pairs(coords, charges, c: float, pairs: ExclusionPairs, dipoles=None,
                                  pref: float = COULOMB, grad: bool = False):
    """:func:`charge_dipole_exclusion` of every chain of a cell at once.

    Returns ``(E, dE/dcoords (N, 3), dE/dc, F (N, 3))`` with the dense routine's conventions; the
    field ``F`` is returned whatever ``dipoles`` is.  A site outside the lists' selection gets no
    field: it has no dipole for one to act on."""
    X = np.asarray(coords, dtype=float).reshape(-1, 3)
    q = np.asarray(charges, dtype=float).reshape(-1)
    N = X.shape[0]
    i, j = pairs.i, pairs.j
    R = X[i] - X[j]
    R[:, 2] -= pairs.k * float(c)
    r2 = np.einsum("ec,ec->e", R, R)
    inv3 = pairs.w / (r2 * np.sqrt(r2))
    qi, qj = q[i], q[j]
    F = 0.5 * pref * (_scatter(j, (inv3 * qi)[:, None] * R, N) - _scatter(i, (inv3 * qj)[:, None] * R, N))
    if dipoles is None:
        return 0.0, None, None, F
    p = np.asarray(dipoles, dtype=float).reshape(N, 3)
    pi, pj = p[i], p[j]
    cross = qi * np.einsum("ec,ec->e", R, pj) - qj * np.einsum("ec,ec->e", R, pi)
    e = -0.5 * pref * float((inv3 * cross).sum())
    gX = gc = None
    if grad:
        G = pref * (inv3[:, None] * (qi[:, None] * pj - qj[:, None] * pi) - (3.0 * inv3 / r2 * cross)[:, None] * R)
        gX = -0.5 * (_scatter(i, G, N) - _scatter(j, G, N))
        gc = 0.5 * float((pairs.k * G[:, 2]).sum())
    return e, gX, gc, F


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
