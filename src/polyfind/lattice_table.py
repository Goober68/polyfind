"""Tabulated chain-pair interaction and an exhaustive packing screen.

The lattice energy of a cell of rigid 1-D-periodic chains is a sum of
*pair-of-chains* interactions.  For two copies of the same periodic chain the
interaction depends only on their relative geometry: the lateral offset vector
``rho``, the relative shift ``dz`` along the chain axis, the two setting angles
and the flip flag.  Rotating the whole pair about z is a symmetry, so

    W(|rho|, alpha1 = phi1 - theta_rho, alpha2 = phi2 - theta_rho, dz, flip)

is a function of four continuous variables and a flag (``theta_rho`` is the
direction of ``rho``).  With chain 2 of the two-chain cell at (1/2, 1/2),

    E(a, b, gamma, phi1, phi2, dz, flip) = E_intra(flip)
        + 1/2 sum_{s != 0}  [ Wself(|rho_s|, e1 (phi1 - theta_s))
                            + Wself(|rho_s|, e2 (phi2 - theta_s)) ]
        + sum_{centred s}   W(|rho_s|, phi1 - theta_s, phi2 - theta_s, dz, flip)

where ``rho_s = i a_vec + j b_vec`` runs over lattice sites within reach,
the centred sites are ``rho_s + (a_vec + b_vec)/2``, ``E_intra(flip)`` is the two
chains' own energy including their z-images (a constant per flip, see
:func:`intra_constants`), ``Wself(r, alpha) = W(r, alpha, alpha, 0, 0)`` and
``e1 = +1``; ``e2 = +1`` for a parallel cell but ``e2 = -1`` when chain 2 is
flipped, because a flipped chain interacting with its own (equally flipped)
lattice images is, after the global rotation ``(x, y, z) -> (x, -y, -z)`` that
undoes the flip, an unflipped chain at setting angle ``-phi2``.  (For a mirror-
symmetric chain such as all-trans PE the sign makes no difference; for a chiral
one, a TG+ helix say, it does.)  The decomposition is exact: it reproduces
``CrystalPacker(chain, 2).energy`` to 1e-10 kcal/mol when W is evaluated with
:func:`chain_pair_energy` instead of the table.

:class:`PairTable` tabulates ``W`` on a grid in ``(r, alpha1, alpha2, dz)`` for
both flips using exactly the pair potential of :class:`~polyfind.pack.CrystalPacker`
(LJ with an energy shift plus damped-shifted-force Coulomb, masked at the
cutoff), :func:`table_energy` evaluates cell energies as lattice sums of
multilinearly interpolated table values, and :func:`fft_screen` evaluates the
*whole* ``(phi1, phi2, dz)`` grid for each ``(a, b)`` with one inverse FFT,
using the fact that rotating a site by ``theta_s`` multiplies the Fourier
coefficient ``(m1, m2, n)`` of ``W`` by ``exp(-i (m1 + m2) theta_s)`` (here by
the Fourier image of the *linearly interpolated* rotation, so that the FFT
landscape is the same function as :func:`table_energy`, not merely close to it).

The build is the one expensive step, and it is a pure function of the chain geometry and
the potential parameters over a grid in which *radii do not interact*, so the radial axis
is split across a process pool (one worker per physical core; see
:func:`default_table_procs`) and the result is bit-identical to the serial build at any
worker count.  It is also memoised, in process and optionally on disk, by
:func:`pair_table`.
"""
from __future__ import annotations

import atexit
from dataclasses import dataclass, field

import numpy as np

from .forcefield import erfc_approx
from .pack import CrystalPacker, PeriodicChain, default_bounds

__all__ = [
    "chain_pair_energy",
    "intra_energy",
    "lattice_sites",
    "direct_lattice_energy",
    "PairTable",
    "table_energy",
    "fft_screen",
    "screen",
    "table_key",
    "pair_table",
    "pair_table_cache_info",
    "clear_pair_table_cache",
    "default_cache_dir",
    "default_table_procs",
    "shutdown_table_pool",
]


# ------------------------------------------------------------------ pair potential
def _pair_v(packer: CrystalPacker, r2: np.ndarray) -> np.ndarray:
    """The packer's exact pair potential (LJ + DSF Coulomb, shifted, masked) at r2.

    Pairwise by construction -- which is the whole basis of the tabulation -- so an Ewald
    packer is refused rather than quietly tabulated without its electrostatics.
    """
    packer._require_dsf("the tabulated chain-pair interaction")
    mask = (r2 < packer.rc ** 2) & (r2 > 1e-8)
    rr = np.sqrt(np.where(mask, r2, 1.0))
    s6 = (packer.lj_x / rr) ** 6
    v = packer.lj_d * (s6 * s6 - 2 * s6) - packer.lj_shift
    v = v + packer.qq * (erfc_approx(packer.alpha * rr, np) / rr - packer.dsf_shift + packer.dsf_force * (rr - packer.rc))
    return np.where(mask, v, 0.0)


def _rotz(X: np.ndarray, ang_deg: np.ndarray) -> np.ndarray:
    """Rotate (..., n, 3) about z by a per-configuration angle (deg)."""
    a = np.deg2rad(np.asarray(ang_deg, dtype=float))[..., None]
    ca, sa = np.cos(a), np.sin(a)
    return np.stack([ca * X[..., 0] - sa * X[..., 1], sa * X[..., 0] + ca * X[..., 1], X[..., 2]], axis=-1)


def chain_pair_energy(packer: CrystalPacker, rho_vec, dz, phi1, phi2, flip, chunk_elems: int = 4_000_000) -> np.ndarray:
    """Interaction energy (kcal/mol) between chain A at the origin with setting angle
    ``phi1`` and chain B at lateral offset ``rho_vec`` with setting angle ``phi2``,
    z-shift ``dz`` and flip flag ``flip``, summed over all z-images of B.

    ``packer`` must be a :class:`~polyfind.pack.CrystalPacker` built with
    ``n_chains=1`` (its per-pair matrices are then the (n, n) blocks used here);
    the pair potential, cutoff and image count are exactly the packer's.
    Vectorised over a batch of geometries: ``rho_vec`` is (G, 2), the rest are
    (G,) or scalars.
    """
    rho = np.atleast_2d(np.asarray(rho_vec, dtype=float))
    G = rho.shape[0]
    dz, phi1, phi2, flip = (np.broadcast_to(np.asarray(v, dtype=float), (G,)) for v in (dz, phi1, phi2, flip))
    X0 = np.asarray(packer.X0, dtype=float)
    n = X0.shape[0]
    ks = np.arange(-packer.K, packer.K + 1, dtype=float)
    nk = len(ks)
    out = np.empty(G)
    g_chunk = max(1, int(chunk_elems // (nk * n * n)))
    for s in range(0, G, g_chunk):
        sl = slice(s, min(s + g_chunk, G))
        m = sl.stop - sl.start
        f = flip[sl][:, None, None]
        Xf = np.where(f > 0.5, np.stack([X0[:, 0], -X0[:, 1], -X0[:, 2]], axis=-1)[None], X0[None])
        A = _rotz(np.broadcast_to(X0[None], (m, n, 3)), phi1[sl])
        B = _rotz(Xf, phi2[sl])
        B = B + np.concatenate([rho[sl], dz[sl][:, None]], axis=1)[:, None, :]
        D = A[:, :, None, :] - B[:, None, :, :]  # (m, n, n, 3)
        d2xy = (D[..., 0] ** 2 + D[..., 1] ** 2)[:, None, :, :]
        r2 = d2xy + (D[:, None, :, :, 2] - packer.chain.c * ks[None, :, None, None]) ** 2
        out[sl] = _pair_v(packer, r2).sum(axis=(1, 2, 3))
    return out


def intra_energy(packer1: CrystalPacker) -> float:
    """Energy of one chain with its own z-images (bonded exclusions applied) plus its
    torsion term: the constant part of the lattice energy, per chain."""
    packer1._require_dsf("intra_energy")
    big = packer1.reach + 20.0
    return float(packer1.energy(np.array([[big, big, 90.0, 0.0, 0.0, 0.0, 0.0]]))[0])


def intra_constants(packer2: CrystalPacker) -> np.ndarray:
    """``[E_intra(flip=0), E_intra(flip=1)]``: the constant part of a two-chain cell,
    isolated by evaluating the packer on a cell far larger than its reach.

    Refused for an Ewald packer: the "constant part" is isolated by evaluating a cell far
    larger than the cutoff, which for a conditionally convergent sum is not a limit at all
    (and for an Ewald sum is a very large image grid).

    The two entries are measured rather than assumed to be equal, because the
    intra-chain term of the *flipped* chain is whatever the packer says it is
    (with the current :mod:`polyfind.pack` the bonded-exclusion matrices of the
    (0,0,k) column are indexed by ``+k`` for the flipped chain as well, which
    makes ``E_intra(flip=1)`` differ from ``E_intra(flip=0)`` by a large constant;
    keeping it as a measured constant means the table reproduces the packer
    either way).  Both entries are independent of the cell and of the angles.
    """
    packer2._require_dsf("intra_constants")
    big = packer2.reach + 20.0
    p = np.array([[big, big, 90.0, 0.0, 0.0, 0.0, f] for f in (0.0, 1.0)])
    return packer2.energy(p)


# ------------------------------------------------------------------ lattice geometry
def lattice_sites(a: float, b: float, gamma: float, r_max: float, centred: bool = False):
    """Polar coordinates (r, theta in deg) of the lattice sites within ``r_max``.

    ``centred=False`` returns ``i a_vec + j b_vec`` with (i, j) != (0, 0);
    ``centred=True`` returns ``(i + 1/2) a_vec + (j + 1/2) b_vec``.
    """
    g = np.deg2rad(gamma)
    av = np.array([a, 0.0])
    bv = np.array([b * np.cos(g), b * np.sin(g)])
    h = min(a * np.sin(g), b * np.sin(g), a, b)
    R = int(np.ceil(r_max / max(h, 1e-3))) + 2
    i, j = np.meshgrid(np.arange(-R, R + 1), np.arange(-R, R + 1), indexing="ij")
    i = i.ravel().astype(float)
    j = j.ravel().astype(float)
    if centred:
        i, j = i + 0.5, j + 0.5
    else:
        keep = (i != 0) | (j != 0)
        i, j = i[keep], j[keep]
    pos = i[:, None] * av[None, :] + j[:, None] * bv[None, :]
    r = np.linalg.norm(pos, axis=1)
    sel = r < r_max
    r = r[sel]
    theta = np.degrees(np.arctan2(pos[sel, 1], pos[sel, 0]))
    return r, theta


def direct_lattice_energy(packer1: CrystalPacker, params, e_intra=None, r_max: float | None = None) -> np.ndarray:
    """Lattice energy of two-chain cells from the decomposition, with ``W`` evaluated
    *directly* (no table).  Comparable with ``CrystalPacker(chain, 2).energy``.

    ``e_intra`` is the pair ``[E_intra(flip=0), E_intra(flip=1)]`` from
    :func:`intra_constants`; it is measured from a two-chain packer if omitted.
    """
    params = np.atleast_2d(np.asarray(params, dtype=float))
    if e_intra is None:
        e_intra = intra_constants(CrystalPacker(packer1.chain, n_chains=2, cutoff=packer1.rc, alpha=packer1.alpha))
    e_intra = np.asarray(e_intra, dtype=float)
    if r_max is None:
        r_max = packer1.rc + 2 * packer1.chain.radius
    out = np.empty(params.shape[0])
    for m, (a, b, gam, phi1, phi2, dz, flip) in enumerate(params):
        e = float(e_intra[int(round(flip))])
        rs, ts = lattice_sites(a, b, gam, r_max, centred=False)
        rho = np.stack([rs, np.zeros_like(rs)], axis=1)  # canonical: rho along x, angles measured from it
        e2sign = -1.0 if flip > 0.5 else 1.0
        for phi, sgn in ((phi1, 1.0), (phi2, e2sign)):
            al = sgn * (phi - ts)
            e += 0.5 * chain_pair_energy(packer1, rho, 0.0, al, al, 0.0).sum()
        rc_, tc = lattice_sites(a, b, gam, r_max, centred=True)
        rhoc = np.stack([rc_, np.zeros_like(rc_)], axis=1)
        e += chain_pair_energy(packer1, rhoc, dz, phi1 - tc, phi2 - tc, flip).sum()
        out[m] = e
    return out


# ------------------------------------------------------------------ the table
def _pot_batch(lj_x6, lj_d, qq, lj_shift, r2, rc, alpha, dsf_shift, dsf_force, dtype=np.float32):
    """The same pair potential, for per-triple parameter vectors broadcast against r2.

    The reference form of the build's pair kernel (kept for tests and for reading; the
    build itself runs :class:`_PairKernel`, whose exact path evaluates the same expression
    with the same float32 rounding, in place).  ``lj_x6`` is ``lj_x ** 6`` (the LJ part is
    a function of ``r^2`` alone, so the sixth power is taken once per pair type instead of
    once per grid point).  ``r2`` is floored at (0.05 A)^2 so that ``(x/r)^12`` cannot
    overflow float32; that only affects pairs whose energy is astronomically large and
    therefore clipped to the table's cap.  Exactly coincident atoms give 0, as in the packer.
    """
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        mask = (r2 < dtype(rc * rc)) & (r2 > dtype(1e-8))
        r2s = np.where(mask, np.maximum(r2, dtype(2.5e-3)), dtype(1.0))
        ir2 = dtype(1.0) / r2s
        s6 = lj_x6 * (ir2 * ir2 * ir2)
        rr = np.sqrt(r2s)
        v = lj_d * (s6 * s6 - 2 * s6) - lj_shift
        v = v + qq * (erfc_approx(dtype(alpha) * rr, np) / rr - dtype(dsf_shift) + dtype(dsf_force) * (rr - dtype(rc)))
        return np.where(mask, v, dtype(0.0))


# Abramowitz-Stegun 7.1.26, the same constants as :func:`polyfind.forcefield.erfc_approx`.
_AS_P = 0.3275911
_AS_C = (0.254829592, -0.284496736, 1.421413741, -1.453152027, 1.061405429)
_R2_FLOOR = 2.5e-3  # (0.05 A)^2: keeps (x/r)^12 finite in float32; see _pot_batch
_SPLINE_R_LO = 1.0  # A: below this every pair type is >= 1e4 kcal/mol, far past any cap


def _pair_type_v_and_dv(lx6, ld, lq, ls, r, rc, alpha, dsf_shift, dsf_force):
    """``V(r)`` of one pair type and ``dV/dr``, in float64, for the spline knots.

    ``V = ld ((x/r)^12 - 2 (x/r)^6) - ls + lq (erfc(a r)/r - dsf_shift + dsf_force (r - rc))``
    with the *polynomial* erfc the kernel uses, differentiated as such (the same
    reasoning as :func:`polyfind.pack._derfc_pos`: the spline should reproduce the kernel,
    not the function the kernel approximates).
    """
    r = np.asarray(r, dtype=float)
    ir2 = 1.0 / (r * r)
    s6 = lx6 * ir2 ** 3
    v = ld * (s6 * s6 - 2 * s6) - ls
    dv = ld * (-12.0 * s6 * s6 + 12.0 * s6) / r
    x = alpha * r
    t = 1.0 / (1.0 + _AS_P * x)
    c1, c2, c3, c4, c5 = _AS_C
    poly = t * (c1 + t * (c2 + t * (c3 + t * (c4 + t * c5))))
    dpoly = c1 + t * (2.0 * c2 + t * (3.0 * c3 + t * (4.0 * c4 + t * 5.0 * c5)))
    ex = np.exp(-x * x)
    erfc = poly * ex
    derfc = (dpoly * (-_AS_P * t * t) - 2.0 * x * poly) * ex  # d/dx
    coul = erfc / r - dsf_shift + dsf_force * (r - rc)
    dcoul = alpha * derfc / r - erfc / (r * r) + dsf_force
    return v + lq * coul, dv + lq * dcoul


def _spline_coefficients(types: np.ndarray, h: float, rc: float, alpha: float, dsf_shift: float, dsf_force: float) -> np.ndarray:
    """Cubic Hermite coefficients of every pair type's ``V(r)`` on the knots ``r = k h``.

    ``(n_types, K, 4)`` float32 with ``V(r) = c0 + c1 f + c2 f^2 + c3 f^3``,
    ``f = r / h - k`` in ``[0, 1)``: the interpolant is C^1, local, and built from the
    analytic value and slope at the knots (error ``<= h^4 max|V''''| / 384`` on each
    interval).  Knots below :data:`_SPLINE_R_LO` hold the constant ``V(r_lo)``: every pair
    type is above 1e4 kcal/mol there, so any table entry containing such a pair is capped
    whichever value is used, and the wall's astronomical derivatives never enter a cubic.
    Knots run to ``rc`` and one beyond, so a clipped ``r <= rc`` always indexes inside.
    """
    K = int(np.floor(rc / h)) + 2
    rk = h * np.arange(K + 1, dtype=float)
    k_lo = int(np.ceil(_SPLINE_R_LO / h))
    out = np.empty((len(types), K, 4), dtype=np.float64)
    for t, (lx6, ld, lq, ls) in enumerate(np.asarray(types, dtype=float)):
        v, dv = _pair_type_v_and_dv(lx6, ld, lq, ls, np.maximum(rk, _SPLINE_R_LO), rc, alpha, dsf_shift, dsf_force)
        v0, v1, d0, d1 = v[:-1], v[1:], h * dv[:-1], h * dv[1:]
        out[t, :, 0] = v0
        out[t, :, 1] = d0
        out[t, :, 2] = 3.0 * (v1 - v0) - 2.0 * d0 - d1
        out[t, :, 3] = 2.0 * (v0 - v1) + d0 + d1
        out[t, :k_lo, 1:] = 0.0  # constant below r_lo
        out[t, :k_lo, 0] = v[k_lo] if k_lo < len(v) else v[-1]
    return out.astype(np.float32)


class _PairKernel:
    """The build's inner kernel: ``r^2`` on a ``(pairs, n_angle, n_angle)`` grid, the pair
    potential at every point and the sum over the pairs, with one set of preallocated
    float32 buffers reused for every chunk.

    Per chunk the terms all belong to one pair *type* (same ``lj_x``, ``lj_d``, charge
    product and shift), so the potential's constants are scalars rather than ``(B, 1, 1)``
    broadcasts.  The grid is assembled from three small per-pair tables:

        r2[b, i, j] = ea[b, i] + eb[b, j] + pe[b, n_angle - 1 + i - j]

    where ``ea`` and ``eb`` carry the terms in one angle each and ``pe[b, m]`` is
    ``-2 rp rq cos(psi_p - psi_q + 2 pi (m - n_angle + 1) / n_angle)``, the ``cos(u - v)``
    term that depends on the two angles only through their difference -- so it is a
    circulant, read through a negative-stride view instead of being formed from four
    ``(B, n_angle, n_angle)`` products.  Every later pass is in place.

    ``spline`` (knot spacing in A) switches the potential to a cubic Hermite interpolant per
    pair type in ``r`` (:func:`_spline_coefficients`): the elementwise work becomes a square
    root, an index, one gather and a cubic instead of a reciprocal, a square root, a
    division, the erfc polynomial with its exp and the LJ powers.  Off by default: the exact
    path reproduces :func:`_pot_batch` to float32 rounding, the spline to its interpolation
    error (measured below 1e-4 kcal/mol per grid point on the screen grid at 0.01 A, median
    1e-7 -- under the float32 accumulation noise of the build itself).  In NumPy it is *not*
    faster than the exact path (a gather costs about five arithmetic passes), which is why
    it is an option and not the default.
    """

    def __init__(self, spec: dict, batch: int, na: int):
        dt = np.dtype(spec["dtype"]).type
        self.dt, self.na = dt, na
        self.rc, self.alpha = dt(spec["rc"]), dt(spec["alpha"])
        self.rc2 = dt(spec["rc"] * spec["rc"])
        self.dsf_shift, self.dsf_force = dt(spec["dsf_shift"]), dt(spec["dsf_force"])
        shape = (batch, na, na)
        self.r2, self.b1, self.b2, self.b3, self.b4, self.b5 = (np.empty(shape, dtype=dt) for _ in range(6))
        self.m1, self.m2 = np.empty(shape, dtype=bool), np.empty(shape, dtype=bool)
        self.acc = np.empty((na, na), dtype=dt)
        self.ea, self.eb = np.empty((batch, na), dtype=dt), np.empty((batch, na), dtype=dt)
        # the circulant: pe[b, n_angle - 1 + i - j] read as a (batch, n_angle, n_angle) view
        self.pe = np.empty((batch, 2 * na - 1), dtype=dt)
        s0, s1 = self.pe.strides
        self.view = np.lib.stride_tricks.as_strided(self.pe[:, na - 1:], shape=shape, strides=(s0, s1, -s1))
        self.h = spec.get("spline")
        if self.h:
            self.idx = np.empty(shape, dtype=np.int64)
            self.inv_h = dt(1.0 / float(self.h))
            self.coef = None  # per flip, set by :meth:`set_pairs`

    def set_pairs(self, A: dict) -> None:
        """Point the kernel at one flip's pair tables (:func:`_pair_arrays`)."""
        self.types = np.asarray(A["types"], dtype=self.dt)
        self.cu, self.cv, self.cd = A["cu"], A["cv"], A["cd"]
        self.rp2q2, self.two_rp, self.two_rq = A["rp2q2"], A["two_rp"], A["two_rq"]
        if self.h:
            self.coef = _spline_coefficients(self.types, float(self.h), float(self.rc), float(self.alpha),
                                             float(self.dsf_shift), float(self.dsf_force))

    def _grid(self, b, zz, rd, B):
        """``r2[b, i, j]`` of the chunk's pairs ``b`` at radius ``rd`` and z-offsets ``zz``."""
        ea, eb, pe = self.ea[:B], self.eb[:B], self.pe[:B]
        # ea[b, i] = rp^2 + rq^2 + r^2 + zz^2 - 2 r rp cos(psi_p + ang_i);  eb[b, j] = 2 r rq cos(psi_q + ang_j)
        np.take(self.cu, b, axis=0, out=ea)
        np.multiply(ea, (rd * self.two_rp[b])[:, None], out=ea)
        np.subtract((self.rp2q2[b] + (rd * rd + zz * zz))[:, None], ea, out=ea)
        np.take(self.cv, b, axis=0, out=eb)
        np.multiply(eb, (rd * self.two_rq[b])[:, None], out=eb)
        np.take(self.cd, b, axis=0, out=pe)
        r2 = self.r2[:B]
        np.add(ea[:, :, None], eb[:, None, :], out=r2)
        np.add(r2, self.view[:B], out=r2)
        m = self.m1[:B]
        np.less(r2, self.rc2, out=m)
        np.greater(r2, self.dt(1e-8), out=self.m2[:B])
        np.logical_and(m, self.m2[:B], out=m)
        return r2, m

    def accumulate(self, b, zz, rd, t: int, plane: np.ndarray) -> None:
        """``plane += sum_b V(r2[b])`` for the chunk's pairs ``b``, all of type ``t``.

        Call inside ``np.errstate(over="ignore", invalid="ignore", divide="ignore")``: the
        floor of ``_R2_FLOOR`` keeps every value finite, but the masked-out points beyond
        the cutoff are evaluated too and can be anything.
        """
        B = len(b)
        r2, m = self._grid(b, zz, rd, B)
        v = self._spline(r2, B, t) if self.h else self._exact(r2, B, t)
        np.multiply(v, m, out=v)
        np.add.reduce(v, axis=0, out=self.acc)
        plane += self.acc

    def _exact(self, r2, B, t):
        dt = self.dt
        lx6, ld, lq, ls = self.types[t]
        np.maximum(r2, dt(_R2_FLOOR), out=r2)
        ir2, s6, rr, v, pl = self.b1[:B], self.b2[:B], self.b3[:B], self.b4[:B], self.b5[:B]
        np.divide(dt(1.0), r2, out=ir2)
        np.multiply(ir2, ir2, out=s6)
        np.multiply(s6, ir2, out=s6)
        np.multiply(s6, lx6, out=s6)  # (x/r)^6
        np.sqrt(r2, out=rr)
        np.subtract(s6, dt(2.0), out=v)
        np.multiply(v, s6, out=v)
        np.multiply(v, ld, out=v)
        np.subtract(v, ls, out=v)  # LJ, energy-shifted
        x, tt, ex = ir2, s6, r2  # buffers free again
        np.multiply(rr, self.alpha, out=x)
        np.multiply(x, dt(_AS_P), out=tt)
        np.add(tt, dt(1.0), out=tt)
        np.divide(dt(1.0), tt, out=tt)
        c1, c2, c3, c4, c5 = (dt(c) for c in _AS_C)
        np.multiply(tt, c5, out=pl)
        np.add(pl, c4, out=pl)
        np.multiply(pl, tt, out=pl)
        np.add(pl, c3, out=pl)
        np.multiply(pl, tt, out=pl)
        np.add(pl, c2, out=pl)
        np.multiply(pl, tt, out=pl)
        np.add(pl, c1, out=pl)
        np.multiply(pl, tt, out=pl)
        np.multiply(x, x, out=ex)
        np.negative(ex, out=ex)
        np.exp(ex, out=ex)
        np.multiply(pl, ex, out=pl)  # erfc(alpha r)
        np.divide(pl, rr, out=pl)
        np.subtract(pl, self.dsf_shift, out=pl)
        np.subtract(rr, self.rc, out=rr)
        np.multiply(rr, self.dsf_force, out=rr)
        np.add(pl, rr, out=pl)
        np.multiply(pl, lq, out=pl)
        np.add(v, pl, out=v)
        return v

    def _spline(self, r2, B, t):
        dt = self.dt
        np.clip(r2, dt(_R2_FLOOR), self.rc2, out=r2)
        g, fl, v = self.b1[:B], self.b2[:B], self.b3[:B]
        np.sqrt(r2, out=g)
        np.multiply(g, self.inv_h, out=g)  # grid coordinate r / h
        np.floor(g, out=fl)
        np.subtract(g, fl, out=g)  # f in [0, 1)
        idx = self.idx[:B]
        np.copyto(idx, fl, casting="unsafe")
        c = np.take(self.coef[t], idx, axis=0)  # (B, na, na, 4)
        np.multiply(c[..., 3], g, out=v)
        np.add(v, c[..., 2], out=v)
        np.multiply(v, g, out=v)
        np.add(v, c[..., 1], out=v)
        np.multiply(v, g, out=v)
        np.add(v, c[..., 0], out=v)
        return v


# Working set of one call of the build's inner kernel, in array elements.  The kernel
# (:class:`_PairKernel`) runs on six preallocated float32 buffers of this size plus two
# masks, so 65,536 elements is about 1.7 MB of live data per worker -- what the old kernel
# held at 16,384, when it built a dozen or so temporaries per call.  The working set, not
# the element width, is what decides whether the build is DRAM-bandwidth-bound: with the
# old kernel, six processes each running the whole build concurrently delivered 1.8x the
# throughput of one at 131,072 elements (8 MB per worker) and 3.9x at 16,384 (1 MB), and
# the serial build was faster at the smaller chunk too (gamma-PVDF screen table 30.4 s
# against 28.2 s).  With the buffered kernel the sweep moved: on a 4-core box, gamma-PVDF
# serial 19.5 / 14.6 / 12.1 / 12.0 s and on four processes 4.9 / 3.8 / 3.3 / 3.0 s at
# 8,192 / 16,384 / 32,768 / 65,536, flat from there to 262,144 -- the per-call overhead
# (~45 NumPy calls per chunk) is what the larger chunk amortises, and 65,536 is the smallest
# size on the plateau, leaving the most L3 headroom for more workers.  It is deliberately
# not part of :func:`table_key`, but it does fix the order of the float32 pair sums, so a
# table built at a different value differs from this one in its last bits (~1e-4 kcal/mol).
BUILD_CHUNK_ELEMS = 65_536


# ------------------------------------------------------------------ the build, by radius
# The build is a pure function of the chain geometry and the potential parameters over a
# grid, and *nothing in it couples one radius to another*: the (pair, z-image) selection,
# the angular planes, the two exchange symmetries and the cap all act within a single
# radius.  So the radial axis is the clean split, and a slice of it is a pure function --
# :func:`_build_block` -- of a small, picklable description of the problem
# (:func:`_build_spec`, ~30 kB) and a list of radius indices.  Assembling the slices in any
# order, from any number of threads or processes, reproduces the serial table bit for bit,
# which ``test_parallel_build_is_bit_identical_to_the_serial_one`` asserts as exact equality
# rather than as a tolerance.
def _build_spec(chain: PeriodicChain, cutoff: float, alpha: float, eps_r: float, r_min: float, dr: float,
                n_angle: int, n_z: int, cap: float, chunk_elems: int, symmetry: bool, dtype,
                spline: float | None = None) -> dict:
    """Everything a radius block of the build needs, as small picklable plain data.

    The four ``(n, n)`` pair-parameter matrices, the chain's coordinates and the radial
    grid: about 30 kB for a 24-atom chain, so a work item can carry the whole problem to a
    worker process without the serialisation eating the parallel gain -- and, in
    particular, without the worker having to rebuild a :class:`CrystalPacker`, which is
    the expensive thing to do per item.
    """
    pk1 = CrystalPacker(chain, n_chains=1, cutoff=cutoff, alpha=alpha, eps_r=eps_r)
    r_max = pk1.rc + 2 * chain.radius
    n_r = int(np.ceil((r_max - r_min) / dr)) + 1
    return {
        "key": f"{id(chain):x}-{_next_spec_serial()}",
        "coords": np.ascontiguousarray(chain.coords, dtype=float),
        "lj_x": np.asarray(pk1.lj_x, dtype=float),
        "lj_d": np.asarray(pk1.lj_d, dtype=float),
        "qq": np.asarray(pk1.qq, dtype=float),
        "lj_shift": np.asarray(pk1.lj_shift, dtype=float),
        "rc": float(pk1.rc),
        "alpha": float(alpha),
        "dsf_shift": float(pk1.dsf_shift),
        "dsf_force": float(pk1.dsf_force),
        "K": int(pk1.K),
        "c": float(chain.c),
        "r_values": r_min + dr * np.arange(n_r),
        "n_angle": int(n_angle),
        "n_z": int(n_z),
        "cap": float(cap),
        "chunk_elems": int(chunk_elems),
        "symmetry": bool(symmetry),
        "dtype": np.dtype(dtype).name,
        "spline": None if spline is None else float(spline),
    }


_SPEC_SERIAL = [0]
_PAIR_CACHE: dict = {}


def _next_spec_serial() -> int:
    _SPEC_SERIAL[0] += 1
    return _SPEC_SERIAL[0]


def _pair_arrays(spec: dict, flip: int) -> dict:
    """The per-flip atom-pair / z-image arrays of the build, memoised per process.

    A worker gets several radius blocks of the same build, so these are computed once per
    process and reused; they are also what would be far too large to send per work item
    (millions of entries for the fine grid, against 30 kB for the spec they come from).

    Besides the (pair, z-image) terms this holds what :class:`_PairKernel` reads per chunk:
    the pair *types* (the distinct ``(lj_x^6, lj_d, qq, lj_shift)`` rows, with the halved
    diagonal of the flip-1 exchange half-sum as types of their own) and, per pair, the
    three angular tables ``cos(psi_p + ang)``, ``cos(psi_q + ang)`` and the circulant
    ``cos(psi_p - psi_q + ang')`` over the ``2 n_angle - 1`` angle differences, so that no
    trigonometry happens per chunk.
    """
    ck = (spec["key"], int(flip))
    hit = _PAIR_CACHE.get(ck)
    if hit is not None:
        return hit
    dtype = np.dtype(spec["dtype"]).type
    X = spec["coords"]
    rad = np.hypot(X[:, 0], X[:, 1])
    psi = np.arctan2(X[:, 1], X[:, 0])
    zc = X[:, 2]
    n, n_z, K, c, na = len(rad), spec["n_z"], spec["K"], spec["c"], spec["n_angle"]
    if spec["symmetry"] and flip == 1:
        P, Q = np.triu_indices(n)  # exchange maps (p, q, t) -> (q, p, t) here
        half = P == Q
    else:
        P, Q = (v.ravel() for v in np.meshgrid(np.arange(n), np.arange(n), indexing="ij"))
        half = None
    pot = (spec["lj_x"], spec["lj_d"], spec["qq"], spec["lj_shift"])
    lxb, ldb, lqb, lsb = (a[P, Q] for a in pot)
    if half is not None:  # v is linear in (lj_d, qq, lj_shift): halve the diagonal terms
        ldb, lqb, lsb = (np.where(half, 0.5 * a, a) for a in (ldb, lqb, lsb))
    lx6 = (lxb ** 6).astype(dtype)
    ld, lq, ls = (a.astype(dtype) for a in (ldb, lqb, lsb))
    par = np.stack([lx6, ld, lq, ls], axis=1)
    types, type_of_pair = np.unique(par, axis=0, return_inverse=True)
    type_of_pair = np.asarray(type_of_pair).ravel()
    rp, rq = rad[P].astype(dtype), rad[Q].astype(dtype)
    psi_q = -psi if flip else psi
    z_q = -zc if flip else zc
    pp, pq = psi[P].astype(dtype), psi_q[Q].astype(dtype)
    ang = (2 * np.pi * np.arange(na) / na).astype(dtype)
    angd = (2 * np.pi * (np.arange(2 * na - 1) - (na - 1)) / na).astype(dtype)
    cu = np.cos(pp[:, None] + ang[None, :])
    cv = np.cos(pq[:, None] + ang[None, :])
    cd = np.cos((pp - pq)[:, None] + angd[None, :])
    cd *= (-2 * rp * rq)[:, None]  # the whole cos(u - v) term, ready to add
    dzpq = zc[P] - z_q[Q]
    # (pair, t) grid with t = j + k n_z, i.e. the z offset of the image is dz_j + k c
    t_all = np.arange(-K * n_z, (K + 1) * n_z)
    zt = (dzpq[:, None] - t_all[None, :] * (c / n_z)).ravel()
    pair_i = np.repeat(np.arange(len(P)), len(t_all))
    j_f = np.tile(np.mod(t_all, n_z), len(P))
    if spec["symmetry"] and flip == 0:
        sel0 = j_f <= n_z // 2
        pair_i, zt, j_f = pair_i[sel0], zt[sel0], j_f[sel0]
    gap = rp[pair_i].astype(float) + rq[pair_i].astype(float)  # closest possible xy approach is r - gap
    out = {"types": types, "cu": cu, "cv": cv, "cd": cd,
           "rp2q2": rp * rp + rq * rq, "two_rp": dtype(2.0) * rp, "two_rq": dtype(2.0) * rq,
           "pair_i": pair_i, "zt": zt, "zt2": zt ** 2, "j_f": j_f, "gap": gap,
           "key_t": j_f * len(types) + type_of_pair[pair_i]}
    if len(_PAIR_CACHE) >= 8:  # one build at a time; never let a long-lived worker grow
        _PAIR_CACHE.clear()
    _PAIR_CACHE[ck] = out
    return out


def _build_block(spec: dict, irs) -> np.ndarray:
    """``W`` for the radii ``irs``: ``(2, len(irs), n_angle, n_angle, n_z)`` float32, finished.

    "Finished" means both exchange symmetries have been applied and the cap taken, so the
    caller only has to write the block into place.  Only the (pair, z-image) terms that can
    reach within the cutoff at that radius are evaluated (the z-window of section 6 of the
    performance review: a pair whose closest possible lateral approach is
    ``r - rad_p - rad_q`` can only reach z-images with ``|z_p - z_q - k c| < sqrt(rc^2 - d^2)``),
    and with ``symmetry`` set the exchange relations

        flip 0:  W(r, a1, a2, dz) = W(r, a2 + 180, a1 + 180, -dz)
        flip 1:  W(r, a1, a2, dz) = W(r, 180 - a2, 180 - a1, dz)

    halve *both* flips.  Both are the same underlying symmetry -- swapping which chain sits
    at the origin -- but they act differently on the sum, so they are applied differently:
    exchange maps the term ``(p, q, t)`` (atom pair, z-image index) to ``(q, p, -t)`` for
    flip 0 and to ``(q, p, t)`` for flip 1.  For flip 0 that is a relation between the ``dz``
    planes ``j`` and ``-j``, so only planes up to ``c/2`` are summed and the rest are filled
    in; for flip 1 it acts inside each plane, so only the ``p <= q`` half of the atom pairs
    is summed (the diagonal at half weight) and the finished block is symmetrised.  Both act
    within one radius, which is why the radial axis is the split.

    Within a plane the surviving terms are grouped by pair type and handed to
    :class:`_PairKernel` in chunks of ``chunk_elems // n_angle^2`` pairs; the chunk sums
    are accumulated in float64 per plane.  The order is fixed by the radius alone, so the
    block is the same from any worker.
    """
    irs = np.atleast_1d(np.asarray(irs, dtype=np.int64))
    na, nz = spec["n_angle"], spec["n_z"]
    dtype = np.dtype(spec["dtype"]).type
    rc, cap = spec["rc"], spec["cap"]
    r_values = spec["r_values"]
    m_half = na // 2
    batch = max(1, spec["chunk_elems"] // (na * na))
    kern = _PairKernel(spec, batch, na)
    out = np.zeros((2, len(irs), na, na, nz), dtype=np.float32)
    for flip in (0, 1):
        A = _pair_arrays(spec, flip)
        pair_i, zt, zt2, gap, key_t = (A[k] for k in ("pair_i", "zt", "zt2", "gap", "key_t"))
        n_types = len(A["types"])
        kern.set_pairs(A)
        for m, ir in enumerate(irs):
            r = r_values[ir]
            keep = np.maximum(0.0, r - gap) ** 2 + zt2 < rc * rc
            if not keep.any():
                continue
            pi, zz, kk = pair_i[keep], zt[keep].astype(dtype), key_t[keep]
            order = np.argsort(kk, kind="stable")
            pi, zz, kk = pi[order], zz[order], kk[order]
            edges = np.searchsorted(kk, np.arange(nz * n_types + 1))
            rd = dtype(r)
            with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
                for j in range(nz):
                    if edges[(j + 1) * n_types] <= edges[j * n_types]:
                        continue
                    plane = np.zeros((na, na))
                    for t in range(n_types):
                        lo, hi = int(edges[j * n_types + t]), int(edges[j * n_types + t + 1])
                        for s in range(lo, hi, batch):
                            e = min(s + batch, hi)
                            kern.accumulate(pi[s:e], zz[s:e], rd, t, plane)
                    out[flip, m, :, :, j] += plane.astype(np.float32)
        if spec["symmetry"] and flip == 0:
            for j in range(nz // 2 + 1, nz):
                src = out[0, :, :, :, (-j) % nz]
                out[0, :, :, :, j] = np.roll(np.swapaxes(src, 1, 2), (m_half, m_half), axis=(1, 2))
        elif spec["symmetry"] and flip == 1:
            # W1 = V + S[V] with S[X][i1, i2] = X[m_half - i2, m_half - i1] and V the
            # p <= q half sum (diagonal at half weight)
            V = out[1]
            out[1] = V + np.roll(np.swapaxes(V, 1, 2)[:, ::-1, ::-1, :], (m_half + 1, m_half + 1), axis=(1, 2))
    np.minimum(out, np.float32(cap), out=out)
    return out


# ------------------------------------------------------------------ the process pool
# Threads could never have worked here, and not only because of the GIL: with the old
# 131,072-element inner chunk the build was DRAM-bandwidth-bound, so there was no headroom
# for *any* form of parallelism to use (see :data:`BUILD_CHUNK_ELEMS`).  With the chunk
# cache-sized, the radial split across processes gives 4.6x on six cores -- and saturates
# there, at the memory system rather than at the core count.  Windows spawns, so the worker
# function is module-level and its arguments are plain picklable data; the pool is kept
# alive between builds because spawning six interpreters that import NumPy costs about
# 2 s, which is more than a small table's whole build.
_POOL: dict = {"ex": None, "n": 0}
_POOL_LOCK = None


def _pool_lock():
    global _POOL_LOCK
    if _POOL_LOCK is None:
        import threading

        _POOL_LOCK = threading.Lock()
    return _POOL_LOCK


def physical_cores() -> int:
    """Number of *physical* cores, falling back to :func:`os.cpu_count` if unknown.

    Not the logical count: the build is DRAM-bandwidth-bound (see
    :data:`BUILD_CHUNK_ELEMS`), so an SMT sibling brings no extra load/store capacity and
    only adds pressure on the shared L3.  Measured on this 6-core / 12-thread i7-8700K,
    12 workers are 15-20% *slower* than 6 (gamma-PVDF screen grid: 7.7 s against 6.4 s).
    """
    import os
    import sys

    try:  # if psutil happens to be installed, it already knows
        import psutil

        n = psutil.cpu_count(logical=False)
        if n:
            return int(n)
    except Exception:
        pass
    try:
        if sys.platform == "win32":
            import ctypes
            from ctypes import wintypes

            k = ctypes.windll.kernel32
            size = wintypes.DWORD(0)
            k.GetLogicalProcessorInformation(None, ctypes.byref(size))
            buf = (ctypes.c_byte * size.value)()
            if k.GetLogicalProcessorInformation(buf, ctypes.byref(size)):
                raw = bytes(bytearray(buf))
                # SYSTEM_LOGICAL_PROCESSOR_INFORMATION on x64: 8-byte mask, 4-byte
                # relationship (0 = RelationProcessorCore), 16-byte union, 32 in all
                n = sum(int.from_bytes(raw[i + 8 : i + 12], "little") == 0 for i in range(0, size.value, 32))
                if n:
                    return n
        else:
            import glob

            sib = {open(p).read().strip() for p in
                   glob.glob("/sys/devices/system/cpu/cpu[0-9]*/topology/thread_siblings_list")}
            if sib:
                return len(sib)
    except Exception:
        pass
    return os.cpu_count() or 1


def default_table_procs(n_r: int | None = None) -> int:
    """Worker processes the table build should use: ``$POLYFIND_TABLE_PROCS``, or one per
    physical core (:func:`physical_cores`).

    Returns 1 -- "build it here, in this process" -- inside a worker of an outer process
    pool (``multiprocessing.parent_process()`` is not ``None``), because nesting the two
    levels only oversubscribes the machine: a parallel
    :func:`~polyfind.pipeline.run_pipeline` already has one process per candidate, and
    measured on the PVDF funnel, five candidate workers each building with two table
    workers ran past 600 s against 223 s unnested.  One worker per physical core is where
    the memory system saturates, wherever that worker comes from.  Set
    ``$POLYFIND_TABLE_PROCS`` to override in either direction.  Capped at ``n_r``, the
    number of radii, since that is the axis being split.
    """
    import multiprocessing
    import os

    env = os.environ.get("POLYFIND_TABLE_PROCS")
    if env:
        try:
            n = int(env)
        except ValueError:
            n = 1
    elif multiprocessing.parent_process() is not None:
        n = 1
    else:
        n = physical_cores()
    if n_r is not None:
        n = min(n, int(n_r))
    return max(1, n)


def _table_pool(n: int):
    """The shared build pool, created on demand and reused across builds."""
    from concurrent.futures import ProcessPoolExecutor

    with _pool_lock():
        if _POOL["ex"] is not None and _POOL["n"] == n:
            return _POOL["ex"]
        _shutdown_pool_locked()
        _POOL["ex"], _POOL["n"] = ProcessPoolExecutor(max_workers=n), n
        return _POOL["ex"]


def _shutdown_pool_locked() -> None:
    ex = _POOL["ex"]
    _POOL["ex"], _POOL["n"] = None, 0
    if ex is not None:
        ex.shutdown(wait=False, cancel_futures=True)


def shutdown_table_pool() -> None:
    """Stop the shared table-build worker processes (they are recreated on demand)."""
    with _pool_lock():
        _shutdown_pool_locked()


atexit.register(shutdown_table_pool)


def _radius_chunks(n_r: int, n_chunks: int) -> list:
    """Strided partitions of the radial axis.

    Strided, not contiguous: the cost of a radius falls steeply with ``r`` (the z-window
    admits far fewer image pairs as the chains separate), so contiguous blocks would leave
    the worker holding the small radii idle for most of the build.
    """
    return [c for c in (np.arange(i, n_r, n_chunks) for i in range(n_chunks)) if len(c)]


def _build_parallel(spec: dict, n_r: int, n_procs: int, W: np.ndarray) -> None:
    """Fill ``W`` from a process pool.  Raises if the pool cannot be created or used."""
    chunks = _radius_chunks(n_r, min(n_r, 4 * n_procs))
    ex = _table_pool(n_procs)
    try:
        for irs, block in zip(chunks, ex.map(_build_block, [spec] * len(chunks), chunks)):
            W[:, irs] = block
    except BaseException:
        shutdown_table_pool()  # a broken pool must not be handed to the next build
        raise


@dataclass
class PairTable:
    """``W(r, alpha1, alpha2, dz, flip)`` on a regular grid, stored as float32.

    Axes of :attr:`W`: ``(flip, r, alpha1, alpha2, dz)`` with
    ``r = r_min + dr i``, ``alpha = 360 i / n_angle`` (deg) and ``dz = c j / n_z``.
    Energies are clipped from above at :attr:`cap` so that interpolation stays
    finite on the repulsive wall.  ``W`` vanishes at and beyond
    ``r_max = cutoff + 2 radius`` (no atom pair can then be within the cutoff).
    """

    chain: PeriodicChain = field(repr=False)
    W: np.ndarray = field(repr=False)
    r_min: float
    dr: float
    n_angle: int
    n_z: int
    c: float
    rc: float
    cap: float
    e_intra: np.ndarray  # (2,) constant part of a two-chain cell, per flip
    build_time: float = 0.0
    _fourier: dict = field(default_factory=dict, repr=False)

    # --- construction ---------------------------------------------------------
    @classmethod
    def build(
        cls,
        chain: PeriodicChain,
        cutoff: float = 8.0,
        alpha: float = 0.2,
        eps_r: float = 1.0,
        r_min: float = 3.0,
        dr: float = 0.05,
        n_angle: int = 72,
        n_z: int = 16,
        cap: float = 50.0,
        chunk_elems: int = BUILD_CHUNK_ELEMS,
        symmetry: bool = True,
        n_threads: int | None = None,
        n_procs: int | None = None,
        dtype=np.float32,
        verbose: bool = False,
        spline: float | None = None,
    ) -> "PairTable":
        """Tabulate ``W`` with the exact pair potential of :class:`CrystalPacker`.

        Each grid point is an exact sum over all atom pairs and all z-images
        (``|k| <= K`` as in the packer); see :func:`_build_block` for the z-window and the
        two exchange symmetries that make one radius cheap.

        ``spline`` (a knot spacing in A, e.g. 0.01) replaces the pair potential by a cubic
        Hermite interpolant per pair type (:class:`_PairKernel`); off by default, since it
        changes the table's values by its interpolation error (below 1e-4 kcal/mol per entry
        at 0.01 A, median 1e-7) without being faster in NumPy; the table's only job is to
        choose starts for the exact polish.

        The radial axis is split across ``n_procs`` worker processes (default
        :func:`default_table_procs`: one per physical core, or 1 inside an outer pool's
        worker).  Radii are independent, so the table is **bit-identical** to the serial
        build whatever the worker count -- the split changes only which worker evaluates
        which radius, never what is summed inside one or in what order.  ``n_procs=1``
        builds in this process (with ``n_threads`` threads, 1 by default), and so does a
        machine or sandbox where a process pool cannot be created.
        """
        import time

        t_start = time.perf_counter()
        if n_angle % 2 or n_z % 2:
            raise ValueError("n_angle and n_z must be even (the exchange symmetries map grid to grid)")
        spec = _build_spec(chain, cutoff, alpha, eps_r, r_min, dr, n_angle, n_z, cap, chunk_elems, symmetry, dtype,
                           spline)
        n_r = len(spec["r_values"])
        W = np.empty((2, n_r, n_angle, n_angle, n_z), dtype=np.float32)
        if n_procs is None:
            n_procs = default_table_procs(n_r)
        n_procs = max(1, min(int(n_procs), n_r))
        fell_back = False
        if n_procs > 1:
            try:
                _build_parallel(spec, n_r, n_procs, W)
            except (OSError, RuntimeError, NotImplementedError, ImportError, EOFError) as exc:
                # e.g. a sandbox without process/spawn support, or a broken pool: fall back
                # to the in-process build rather than failing (as pipeline.py does).
                if verbose:
                    print(f"table {chain.name}: process pool unavailable ({exc!r}); building in this process")
                n_procs, fell_back = 1, True
        if n_procs == 1:
            from concurrent.futures import ThreadPoolExecutor

            # Threads used to be worth 1.3x here, with a 131,072-element inner chunk.  They
            # are now worse than useless: the cache-sized chunk of :data:`BUILD_CHUNK_ELEMS`
            # makes each NumPy call eight times shorter, so the loop is pure GIL contention
            # -- measured on the gamma-PVDF screen grid, 35.2 s on one thread against 61.3 s
            # on four.  Kept as an option because a future coarser grid could reverse that.
            if n_threads is None:
                n_threads = 1
            n_threads = max(1, min(int(n_threads), n_r))
            if n_threads > 1:
                chunks = _radius_chunks(n_r, n_threads)
                with ThreadPoolExecutor(max_workers=n_threads) as ex:
                    for irs, block in zip(chunks, ex.map(_build_block, [spec] * len(chunks), chunks)):
                        W[:, irs] = block
            else:
                W[...] = _build_block(spec, np.arange(n_r))
        dt = time.perf_counter() - t_start
        if verbose:
            how = (f"{n_procs} worker processes" if n_procs > 1 else
                   f"this process, {n_threads} thread" + ("" if n_threads == 1 else "s"))
            if fell_back:
                how += ", after the pool failed"
            print(f"table {chain.name}: {tuple(W.shape)} = {W.nbytes / 1e6:.0f} MB in {dt:.1f} s ({how})")
        pk2 = CrystalPacker(chain, n_chains=2, cutoff=cutoff, alpha=alpha, eps_r=eps_r)
        return cls(chain=chain, W=W, r_min=r_min, dr=dr, n_angle=n_angle, n_z=n_z, c=chain.c, rc=spec["rc"],
                   cap=cap, e_intra=intra_constants(pk2), build_time=dt)

    # --- lookup ---------------------------------------------------------------
    @property
    def n_r(self) -> int:
        return self.W.shape[1]

    @property
    def r_max(self) -> float:
        return self.r_min + self.dr * (self.n_r - 1)

    @property
    def nbytes(self) -> int:
        return int(self.W.nbytes)

    def interpolate(self, r, alpha1, alpha2, dz, flip) -> np.ndarray:
        """Multilinear interpolation of ``W``, periodic in the angles (deg) and in dz.

        ``r`` below ``r_min`` is clamped (the wall is capped anyway); beyond
        ``r_max`` the result is exactly zero.  Fully vectorised.

        The interpolation weight is taken from the *unwrapped* grid coordinate, not from
        the wrapped index: a tiny negative angle (``phi - theta_s`` when a lattice site
        sits a rounding error above a grid angle, which an oblique cell produces readily)
        gives ``a1 % 360.0 == 360.0`` exactly, whose wrapped index is 0 but whose fraction
        must be 0 and not ``n_angle``.  Getting that wrong put weights of -47 and +48 on
        two neighbouring table entries and made whole cells score hundreds of kcal/mol
        below anything physical.
        """
        r, a1, a2, z, f = np.broadcast_arrays(
            np.asarray(r, dtype=float), np.asarray(alpha1, dtype=float), np.asarray(alpha2, dtype=float),
            np.asarray(dz, dtype=float), np.asarray(flip, dtype=np.int64) if np.ndim(flip) else np.int64(flip))
        na, nz = self.n_angle, self.n_z
        fr = np.clip((r - self.r_min) / self.dr, 0.0, self.n_r - 1.0)
        i0 = np.minimum(fr.astype(np.int64), self.n_r - 2)
        tr = fr - i0
        ja, jb = (a1 % 360.0) * (na / 360.0), (a2 % 360.0) * (na / 360.0)
        jz = (z % self.c) * (nz / self.c)
        fa, fb, fz = ja.astype(np.int64), jb.astype(np.int64), jz.astype(np.int64)  # >= 0, so floor
        ta, tb, tz = ja - fa, jb - fb, jz - fz
        ia, ib, iz = fa % na, fb % na, fz % nz
        out = np.zeros(np.shape(r), dtype=float)
        for da in (0, 1):
            wa = 1 - ta if da == 0 else ta
            iia = (ia + da) % na
            for db in (0, 1):
                wb = (1 - tb if db == 0 else tb) * wa
                iib = (ib + db) % na
                for dj in (0, 1):
                    wz = (1 - tz if dj == 0 else tz) * wb
                    iiz = (iz + dj) % nz
                    for di in (0, 1):
                        w = (1 - tr if di == 0 else tr) * wz
                        out += w * self.W[f, i0 + di, iia, iib, iiz]
        return np.where(r > self.r_max, 0.0, out)

    # --- Fourier coefficients (for the FFT screen) ----------------------------
    def fourier(self, flip: int) -> np.ndarray:
        """FFT of ``W`` over ``(alpha1, alpha2, dz)`` per radius, ``(n_r, na, na, n_z//2+1)``."""
        key = int(flip)
        if key not in self._fourier:
            self._fourier[key] = np.fft.rfftn(self.W[key].astype(np.float64), axes=(1, 2, 3)).astype(np.complex64)
        return self._fourier[key]

    def radial_weights(self, r: np.ndarray):
        """Index and weight of the linear interpolation in r (the same rule as
        :meth:`interpolate`), so that Fourier coefficients can be interpolated in r."""
        fr = np.clip((np.asarray(r, dtype=float) - self.r_min) / self.dr, 0.0, self.n_r - 1.0)
        i0 = np.minimum(fr.astype(np.int64), self.n_r - 2)
        return i0, fr - i0


# ------------------------------------------------------------------ table cache
# A table costs 1.6-8 s to build on the screen grid of :data:`polyfind.pack.SCREEN_TABLE`
# and 18-87 s on this module's finer default (six worker processes on a 6-core i7-8700K;
# 5-25 s and 40-160 s before the build was parallelised), and depends only on the chain geometry
# (coordinates, elements, charges, bonded-exclusion topology, repeat length) and on the
# potential and grid parameters -- never on the cell.  So one build serves every
# ``pack()`` call on the same conformation, and (through the optional on-disk copy) every
# process of a multi-candidate run and every later run.
_KEY_PARAMS = ("cutoff", "alpha", "eps_r", "r_min", "dr", "n_angle", "n_z", "cap", "dtype")
_TABLE_CACHE: dict[str, PairTable] = {}
_TABLE_STATS = {"builds": 0, "memory_hits": 0, "disk_hits": 0}


def _build_defaults() -> dict:
    import inspect

    sig = inspect.signature(PairTable.build)
    return {k: sig.parameters[k].default for k in _KEY_PARAMS}


def table_key(chain: PeriodicChain, **kw) -> str:
    """Cache key of the ``W`` table for ``chain`` and the given potential/grid settings.

    Everything the table's *contents* depend on goes in: the chain's coordinates,
    elements, charges, repeat length and bonded-exclusion matrices, and the potential
    parameters (``cutoff``, ``alpha``, ``eps_r``) and grid (``r_min``, ``dr``,
    ``n_angle``, ``n_z``, ``cap``, ``dtype``, and ``spline`` when set).  ``symmetry``,
    ``chunk_elems``, ``n_threads``, ``n_procs`` and ``verbose`` only change how the same
    numbers are produced -- the parallel build is bit-identical to the serial one -- and are
    deliberately not part of the key.  The cell parameters are not part of it either --
    that is the point of the table.
    """
    import hashlib

    p = _build_defaults()
    for k, v in kw.items():
        if k in p:
            p[k] = v
    h = hashlib.blake2b(digest_size=12)
    h.update(np.ascontiguousarray(chain.coords, dtype=np.float64).tobytes())
    h.update(np.ascontiguousarray(chain.charges, dtype=np.float64).tobytes())
    h.update("|".join(chain.elements).encode())
    for k in sorted(chain.same_site_scale):
        h.update(f"{k}:".encode())
        h.update(np.ascontiguousarray(chain.same_site_scale[k], dtype=np.float64).tobytes())
    h.update(repr((float(chain.c),) + tuple(
        np.dtype(p[k]).name if k == "dtype" else float(p[k]) for k in _KEY_PARAMS)).encode())
    if kw.get("spline") is not None:  # the splined build is a different table; the exact keys are unchanged
        h.update(f"spline:{float(kw['spline'])!r}".encode())
    tag = "".join(ch if ch.isalnum() else "_" for ch in str(chain.name))[:24]
    return f"{tag}-{h.hexdigest()}"


def default_cache_dir() -> "str | None":
    """On-disk table cache directory, from ``$POLYFIND_TABLE_CACHE``; ``None`` if unset."""
    import os

    return os.environ.get("POLYFIND_TABLE_CACHE") or None


def _load_npz(path, chain: PeriodicChain) -> PairTable:
    with np.load(path) as z:
        m = z["meta"]
        return PairTable(chain=chain, W=z["W"], r_min=float(m[0]), dr=float(m[1]), n_angle=int(m[2]),
                         n_z=int(m[3]), c=float(m[4]), rc=float(m[5]), cap=float(m[6]),
                         e_intra=z["e_intra"], build_time=float(m[7]))


def _save_npz(path, table: PairTable) -> None:
    import os

    tmp = f"{path}.{os.getpid()}.tmp.npz"
    np.savez(tmp, W=table.W, e_intra=np.asarray(table.e_intra, dtype=float),
             meta=np.array([table.r_min, table.dr, table.n_angle, table.n_z, table.c, table.rc,
                            table.cap, table.build_time], dtype=float))
    os.replace(tmp, path)  # atomic: concurrent builders never see a half-written table


def pair_table(chain: PeriodicChain, cache_dir: "str | None" = "env", memory: bool = True,
               verbose: bool = False, **build_kw) -> PairTable:
    """:meth:`PairTable.build`, memoised on :func:`table_key`.

    Repeated calls for the same conformation and potential return *the same object*
    (``memory=False`` disables the in-process cache).  ``cache_dir`` adds an on-disk
    copy (``<key>.npz``, 7 MB on the screen grid, ~110 MB on the fine one) so that a
    fresh process -- one worker of a multi-candidate run, or the next run of the
    pipeline -- pays the build once too; the default ``"env"`` means
    "``$POLYFIND_TABLE_CACHE`` if it is set, else no disk cache" (opt-in, so nothing is
    written to a user's disk unasked), ``None`` disables it and any other value is used
    as the directory.  The write is atomic, so concurrent builders cannot read a
    half-written file, and an unreadable one is ignored rather than fatal.
    """
    import os

    key = table_key(chain, **build_kw)
    if memory and key in _TABLE_CACHE:
        _TABLE_STATS["memory_hits"] += 1
        return _TABLE_CACHE[key]
    if cache_dir == "env":
        cache_dir = default_cache_dir()
    path = os.path.join(cache_dir, key + ".npz") if cache_dir else None
    if path and os.path.exists(path):
        try:
            table = _load_npz(path, chain)
            _TABLE_STATS["disk_hits"] += 1
            if verbose:
                print(f"table {chain.name}: loaded {table.nbytes / 1e6:.0f} MB from {path}")
            if memory:
                _TABLE_CACHE[key] = table
            return table
        except Exception as exc:  # a truncated or stale file must not be fatal
            if verbose:
                print(f"table {chain.name}: ignoring unreadable cache file {path} ({exc!r})")
    table = PairTable.build(chain, verbose=verbose, **build_kw)
    _TABLE_STATS["builds"] += 1
    if memory:
        _TABLE_CACHE[key] = table
    if path:
        os.makedirs(cache_dir, exist_ok=True)
        try:
            _save_npz(path, table)
        except OSError as exc:
            if verbose:
                print(f"table {chain.name}: could not write {path} ({exc!r})")
    return table


def pair_table_cache_info() -> dict:
    """``{'entries', 'bytes', 'builds', 'memory_hits', 'disk_hits'}`` for the in-process cache."""
    return {"entries": len(_TABLE_CACHE), "bytes": sum(t.nbytes for t in _TABLE_CACHE.values()), **_TABLE_STATS}


def clear_pair_table_cache() -> None:
    """Drop the in-process tables (the on-disk copies, if any, are left alone)."""
    _TABLE_CACHE.clear()


# ------------------------------------------------------------------ lattice sums
def _self_sign(flip) -> float:
    """Sign of the setting angle in the self term of chain 2.

    A flipped chain interacting with its own (equally flipped) lattice images is,
    after the global rotation (x, y, z) -> (x, -y, -z), an unflipped chain whose
    setting angle and site angles are both negated.
    """
    return -1.0 if int(round(float(flip))) else 1.0


def table_energy(table: PairTable, params) -> np.ndarray:
    """Lattice energy per cell (kcal/mol) from interpolated table values.

    ``params`` is (M, 7) in the ``CrystalPacker.energy`` convention
    ``(a, b, gamma, phi1, phi2, dz, flip)``; the result includes the constant
    intra-chain part, so it is directly comparable with ``CrystalPacker.energy``.
    """
    params = np.atleast_2d(np.asarray(params, dtype=float))
    r_max = table.r_max
    out = np.empty(params.shape[0])
    for m, (a, b, gam, phi1, phi2, dz, flip) in enumerate(params):
        fl = int(round(flip))
        rs, ts = lattice_sites(a, b, gam, r_max, centred=False)
        e = float(table.e_intra[fl])
        for phi, sgn in ((phi1, 1.0), (phi2, _self_sign(fl))):
            al = sgn * (phi - ts)
            e += 0.5 * float(table.interpolate(rs, al, al, 0.0, 0).sum())
        rc_, tc = lattice_sites(a, b, gam, r_max, centred=True)
        e += float(table.interpolate(rc_, phi1 - tc, phi2 - tc, dz, fl).sum())
        out[m] = e
    return out


@dataclass
class ScreenResult:
    """Grid minima of a table screen over (a, b) at fixed gamma."""

    a_values: np.ndarray
    b_values: np.ndarray
    gamma: float
    flips: tuple
    energy: np.ndarray  # (n_a, n_b, n_flip) minimum over the (phi1, phi2, dz) grid
    phi1: np.ndarray
    phi2: np.ndarray
    dz: np.ndarray
    top: np.ndarray  # (n_top, 7) parameter rows, sorted by energy
    top_energy: np.ndarray
    time: float = 0.0

    def rows(self) -> np.ndarray:
        return self.top


def _shift_multiplier(theta_deg, n_angle: int) -> np.ndarray:
    """DFT multiplier of a *linearly interpolated* rotation of the angle axis by
    ``theta``: the exact Fourier image of what :meth:`PairTable.interpolate` does.

    ``theta_deg`` may be an array: the result has one more trailing axis of length
    ``n_angle`` (``(n_angle,)`` for a scalar), complex64.
    """
    theta = np.asarray(theta_deg, dtype=float)
    g = (theta % 360.0) * (n_angle / 360.0)
    q = np.floor(g)
    f = g - q
    m = np.fft.fftfreq(n_angle) * n_angle  # 0, 1, ..., -1
    ph = np.exp(-2j * np.pi * m * q[..., None] / n_angle)
    return (ph * ((1.0 - f)[..., None] + f[..., None] * np.exp(-2j * np.pi * m / n_angle))).astype(np.complex64)


def _lattice_sites_batch(a, b, gamma: float, r_max: float, centred: bool, chunk_elems: int = 2_000_000):
    """:func:`lattice_sites` for many cells at once.

    Returns flat arrays ``(cell, r, theta)`` sorted by cell index, the sites of each cell
    in exactly the order (and with exactly the values) that :func:`lattice_sites` gives
    them: one ``(i, j)`` grid large enough for the smallest cell serves every cell, and a
    site the per-cell grid would not have generated lies beyond ``r_max`` (``r >= h max(|i|,
    |j|)`` with ``h`` the shortest lattice height) and is dropped by the same test.
    """
    a = np.atleast_1d(np.asarray(a, dtype=float))
    b = np.atleast_1d(np.asarray(b, dtype=float))
    g = np.deg2rad(gamma)
    h = np.minimum(np.minimum(a * np.sin(g), b * np.sin(g)), np.minimum(a, b))
    R = int(np.ceil(r_max / max(float(h.min(initial=np.inf)), 1e-3))) + 2
    i, j = np.meshgrid(np.arange(-R, R + 1), np.arange(-R, R + 1), indexing="ij")
    i = i.ravel().astype(float)
    j = j.ravel().astype(float)
    if centred:
        i, j = i + 0.5, j + 0.5
    else:
        keep = (i != 0) | (j != 0)
        i, j = i[keep], j[keep]
    bx, by = b * np.cos(g), b * np.sin(g)
    n_cells, K = len(a), len(i)
    step = max(1, chunk_elems // max(K, 1))
    cells, rs, ths = [], [], []
    for c0 in range(0, n_cells, step):
        sl = slice(c0, min(c0 + step, n_cells))
        x = i[None, :] * a[sl, None] + j[None, :] * bx[sl, None]
        y = i[None, :] * 0.0 + j[None, :] * by[sl, None]
        r = np.sqrt(x * x + y * y)
        sel = r < r_max
        cells.append(np.broadcast_to(np.arange(c0, c0 + sel.shape[0])[:, None], sel.shape)[sel])
        rs.append(r[sel])
        ths.append(np.degrees(np.arctan2(y[sel], x[sel])))
    if not cells:
        return np.zeros(0, dtype=np.int64), np.zeros(0), np.zeros(0)
    return np.concatenate(cells), np.concatenate(rs), np.concatenate(ths)


def _cell_sums(v: np.ndarray, cell: np.ndarray, n_cells: int) -> np.ndarray:
    """Sum the rows of ``v`` (``(n_sites, k)``, the sites of a cell contiguous and ``cell``
    non-decreasing) over each cell.

    One ``np.sum`` per distinct site count, over a C-contiguous ``(cells, k, n_sites)``
    block, so that every cell gets exactly the pairwise summation NumPy would do for that
    cell alone (bit-identical to a per-cell ``values.sum(axis=1)``, which decides ties
    between symmetry-equivalent grid points).  Cells without sites sum to zero.
    """
    out = np.zeros((n_cells, v.shape[1]), dtype=v.dtype)
    counts = np.bincount(cell, minlength=n_cells)
    first = np.concatenate(([0], np.cumsum(counts)[:-1]))
    for n in np.unique(counts):
        if n == 0:
            continue
        cells_n = np.flatnonzero(counts == n)
        rows = first[cells_n][:, None] + np.arange(n)[None, :]
        out[cells_n] = np.ascontiguousarray(v[rows].transpose(0, 2, 1)).sum(axis=2)
    return out


def _interpolate_self(table: PairTable, r, alpha) -> np.ndarray:
    """``table.interpolate(r, alpha, alpha, 0.0, 0)``, bit for bit, for the self terms.

    Same multilinear rule in the same order of accumulation as :meth:`PairTable.interpolate`,
    minus the eight ``dz`` terms that carry an exactly-zero weight at ``dz = 0`` (adding
    ``0 * W`` changes nothing), so half the gathers.
    """
    r, a1 = np.broadcast_arrays(np.asarray(r, dtype=float), np.asarray(alpha, dtype=float))
    na = table.n_angle
    fr = np.clip((r - table.r_min) / table.dr, 0.0, table.n_r - 1.0)
    i0 = np.minimum(fr.astype(np.int64), table.n_r - 2)
    tr = fr - i0
    ja = (a1 % 360.0) * (na / 360.0)
    fa = ja.astype(np.int64)  # >= 0, so floor
    ta = ja - fa
    ia = fa % na
    W0 = table.W[0]
    out = np.zeros(np.shape(r), dtype=float)
    for da in (0, 1):
        wa = 1 - ta if da == 0 else ta
        iia = (ia + da) % na
        for db in (0, 1):
            wb = (1 - ta if db == 0 else ta) * wa
            iib = (ia + db) % na
            wz = (1 - 0.0) * wb  # the dz = 0 plane, weight (1 - tz) with tz = 0
            for di in (0, 1):
                w = (1 - tr if di == 0 else tr) * wz
                out += w * W0[i0 + di, iia, iib, 0]
    return np.where(r > table.r_max, 0.0, out)


def fft_screen(table: PairTable, a_values, b_values, gamma: float = 90.0, flips=(0, 1), n_top: int = 20,
               min_separation: float = 0.5, batch: int = 16, n_threads: int | None = None,
               ab_symmetry: bool = False) -> ScreenResult:
    """Evaluate the lattice energy on the *whole* ``(phi1, phi2, dz)`` grid for every
    ``(a, b)`` by one inverse FFT per cell.

    The cross term is a sum over the centred lattice sites of ``W`` rotated by
    ``theta_s`` in both angles; a rotation is a multiplication of the Fourier
    coefficients, and linear interpolation in ``r`` commutes with the transform, so
    the whole ``(n_angle, n_angle, n_z)`` landscape costs one inverse FFT.  The self
    terms depend on a single angle each and are interpolated directly on the phi
    grid.  The energies are identical (to float32) to :func:`table_energy` at the
    same grid points.

    Everything that depends only on the cell geometry is computed for all cells at once
    (the lattice sites, the self terms in one :meth:`PairTable.interpolate` call, the
    rotation multipliers and their per-radius outer products), and the Fourier
    accumulation, inverse FFT and grid minimum are done for ``batch`` cells per NumPy
    call; there is no per-cell Python work.  ``n_threads`` is accepted for compatibility
    and ignored: the screen is single-threaded (the earlier per-cell thread pool only
    contended for the GIL).

    ``ab_symmetry=True`` skips the cells with ``b < a``: at ``gamma = 90`` swapping the
    axes is a 90 deg rotation of the same lattice (``phi1, phi2 -> phi1 + 90, phi2 + 90``
    leaves the energy unchanged), so the minimum over the angle/dz grid at ``(b, a)``
    equals the one at ``(a, b)`` and half the FFTs are redundant.  Skipped cells are
    left at ``+inf`` in :attr:`ScreenResult.energy`.  It is only valid at ``gamma = 90``
    (elsewhere the swap is a reflection, which a chiral chain does not admit) and is
    rejected otherwise.

    A cell that puts two chain axes closer than ``table.r_min`` is left at ``+inf``:
    :meth:`PairTable.interpolate` *clamps* below ``r_min``, and W at ``r_min`` is a real
    (often attractive) value rather than the wall the true geometry would have, so such a
    cell would otherwise be scored far too low and win the screen.  With ``r_min = 3 A``
    and any real chain the two chains overlap there anyway, so nothing packable is lost;
    the effect is invisible at ``gamma = 90`` but dominates an oblique lattice, whose
    centred site ``(a_vec + b_vec)/2`` gets short quickly as gamma leaves 90.
    """
    import time

    t0 = time.perf_counter()
    a_values = np.atleast_1d(np.asarray(a_values, dtype=float))
    b_values = np.atleast_1d(np.asarray(b_values, dtype=float))
    flips = tuple(int(f) for f in flips)
    if ab_symmetry and abs(gamma - 90.0) > 1e-9:
        raise ValueError("ab_symmetry is only a symmetry at gamma = 90 deg")
    na, nz, c = table.n_angle, table.n_z, table.c
    nz2 = nz // 2 + 1
    phi = 360.0 * np.arange(na) / na
    r_max = table.r_max
    shape = (len(a_values), len(b_values), len(flips))
    energy = np.full(shape, np.inf)
    best = np.zeros(shape + (3,))

    # --- the cells to screen, and their lattice sites, all at once ---------------------
    IA, IB = np.meshgrid(np.arange(len(a_values)), np.arange(len(b_values)), indexing="ij")
    IA, IB = IA.ravel(), IB.ravel()
    if ab_symmetry:
        keep = b_values[IB] >= a_values[IA] - 1e-9
        IA, IB = IA[keep], IB[keep]
    A, B = a_values[IA], b_values[IB]
    n_cells = len(A)
    cs, rs, ts = _lattice_sites_batch(A, B, gamma, r_max, centred=False)
    cc, rc_, tc = _lattice_sites_batch(A, B, gamma, r_max, centred=True)
    r_min_cell = np.full(n_cells, np.inf)
    np.minimum.at(r_min_cell, cs, rs)
    np.minimum.at(r_min_cell, cc, rc_)
    ok = ~(r_min_cell < table.r_min)  # no table data below r_min: left at +inf (see the docstring)
    if not ok.all():
        new_index = np.cumsum(ok) - 1
        IA, IB = IA[ok], IB[ok]
        sel = ok[cs]
        cs, rs, ts = new_index[cs[sel]], rs[sel], ts[sel]
        sel = ok[cc]
        cc, rc_, tc = new_index[cc[sel]], rc_[sel], tc[sel]
        n_cells = len(IA)
    # --- self terms of every cell in one interpolation (they do not depend on flip) -------
    da = phi[None, :] - ts[:, None]  # (n_sites, na)
    self_pos = 0.5 * _cell_sums(_interpolate_self(table, rs[:, None], da), cs, n_cells)
    self_neg = None
    if any(_self_sign(f) < 0 for f in flips):
        self_neg = 0.5 * _cell_sums(_interpolate_self(table, rs[:, None], -da), cs, n_cells)

    # --- centred sites grouped by (cell, radius): sites at the same radius (a rectangular
    # lattice gives them in fours) share the coefficient array, so their rotation
    # multipliers are summed first, as (na, na) outer products -----------------------------
    keys = np.round(rc_, 9)
    order = np.lexsort((keys, cc))  # stable: sites keep their order inside a group
    cc, rc_, tc, keys = cc[order], rc_[order], tc[order], keys[order]
    if len(cc):
        new_group = np.concatenate(([True], (cc[1:] != cc[:-1]) | (keys[1:] != keys[:-1])))
    else:
        new_group = np.zeros(0, dtype=bool)
    g_start = np.flatnonzero(new_group)  # first site of each group
    n_groups = len(g_start)
    g_cell = cc[g_start]
    g_i0, g_w = table.radial_weights(rc_[g_start])
    w1 = (1.0 - g_w).astype(np.float32)[:, None, None, None]
    w2 = g_w.astype(np.float32)[:, None, None, None]
    # position of each site inside its group and of each group inside its cell: the sums
    # below add the p-th member of every group (cell) in one vectorised operation, p = 0,
    # 1, ..., so that each group (cell) is accumulated left to right in the same order as a
    # per-cell loop would, and the complex64 rounding is the same
    g_id = np.cumsum(new_group) - 1  # group of each site
    s_pos = np.arange(len(cc)) - g_start[g_id]
    g_counts = np.bincount(g_cell, minlength=n_cells)
    g_first_of_cell = np.concatenate(([0], np.cumsum(g_counts)))  # (n_cells + 1,): a batch of cells is a slice of groups
    g_pos = np.arange(n_groups) - g_first_of_cell[g_cell]
    s_counts = np.bincount(cc, minlength=n_cells)
    s_first_of_cell = np.concatenate(([0], np.cumsum(s_counts)))
    LL = np.empty((n_groups, na, na), dtype=np.complex64)
    for c0 in range(0, n_cells, 4 * batch):
        c1 = min(c0 + 4 * batch, n_cells)
        s0, s1 = s_first_of_cell[c0], s_first_of_cell[c1]
        if s1 == s0:
            continue
        L = _shift_multiplier(tc[s0:s1], na)  # (n_sites, na)
        outer = L[:, :, None] * L[:, None, :]
        pos, gid = s_pos[s0:s1], g_id[s0:s1]
        LL[gid[pos == 0]] = outer[pos == 0]
        for p in range(1, int(pos.max()) + 1):
            sel = pos == p
            LL[gid[sel]] += outer[sel]

    # --- per flip: Fourier accumulation, inverse FFT and grid minimum, batched over cells ---
    for fi, flip in enumerate(flips):
        Wh = table.fourier(flip)
        self2 = self_pos if _self_sign(flip) > 0 else self_neg
        e0 = float(table.e_intra[flip])
        for c0 in range(0, n_cells, batch):
            c1 = min(c0 + batch, n_cells)
            nb = c1 - c0
            q0, q1 = g_first_of_cell[c0], g_first_of_cell[c1]
            acc = np.zeros((nb, na, na, nz2), dtype=np.complex64)
            gsl = slice(q0, q1)
            pos, cell = g_pos[gsl], g_cell[gsl] - c0
            for p in range(int(pos.max(initial=-1)) + 1):
                sel = np.flatnonzero(pos == p) + q0
                i0 = g_i0[sel]
                # (1 - w) Wh[i0] + w Wh[i0 + 1], then LL * coef, each product in the operand
                # order of a per-cell loop (NumPy's complex64 product is not order symmetric)
                coef = Wh[i0]
                coef *= w1[sel]
                tmp = Wh[i0 + 1]
                tmp *= w2[sel]
                coef += tmp
                np.multiply(LL[sel, :, :, None], coef, out=coef)
                if len(sel) == nb:
                    acc += coef  # every cell of the batch has a p-th group
                else:
                    acc[cell[sel - q0]] += coef
            E = np.fft.irfftn(acc, s=(na, na, nz), axes=(1, 2, 3))
            E = E + self_pos[c0:c1, :, None, None]  # float64 from here
            E += self2[c0:c1, None, :, None]
            E += e0
            flat = E.reshape(nb, -1)
            k = flat.argmin(axis=1)
            i1, i2, j = np.unravel_index(k, (na, na, nz))
            ia, ib = IA[c0:c1], IB[c0:c1]
            energy[ia, ib, fi] = flat[np.arange(nb), k]
            best[ia, ib, fi, 0] = phi[i1]
            best[ia, ib, fi, 1] = phi[i2]
            best[ia, ib, fi, 2] = c * j / nz

    # top distinct cells: (a, b) and (b, a) are the same cell (a 90 deg rotation of it)
    order = np.argsort(energy, axis=None)
    rows, es, keys = [], [], []
    for idx in order:
        ia, ib, fi = np.unravel_index(idx, shape)
        if not np.isfinite(energy[ia, ib, fi]):
            break
        a, b = float(a_values[ia]), float(b_values[ib])
        p = np.array([a, b, gamma, best[ia, ib, fi, 0], best[ia, ib, fi, 1], best[ia, ib, fi, 2], flips[fi]])
        k = (min(a, b), max(a, b), flips[fi])
        if all(abs(k[0] - q[0]) + abs(k[1] - q[1]) > min_separation or k[2] != q[2] for q in keys):
            rows.append(p)
            es.append(energy[ia, ib, fi])
            keys.append(k)
        if len(rows) >= n_top:
            break
    return ScreenResult(a_values=a_values, b_values=b_values, gamma=gamma, flips=flips, energy=energy,
                        phi1=best[..., 0], phi2=best[..., 1], dz=best[..., 2],
                        top=np.array(rows), top_energy=np.array(es), time=time.perf_counter() - t0)


def screen(chain: PeriodicChain, n_top: int = 20, a_range=None, b_range=None, da: float = 0.1, db: float = 0.1,
           gamma: float = 90.0, flips=(0, 1), table: PairTable | None = None, min_separation: float = 0.5,
           ab_symmetry: bool = False, **build_kw):
    """Exhaustive table screen over ``(a, b, phi1, phi2, dz, flip)``.

    Returns ``(params, energies, result)``: ``params`` is an ``(n_top, 7)`` array of
    candidate cells in the :class:`CrystalPacker` convention, ``energies`` their
    table energies (comparable with ``CrystalPacker.energy``), and ``result`` the
    full :class:`ScreenResult`.  Without an explicit ``table`` the cached builder
    :func:`pair_table` is used, so a second call on the same conformation is free.
    """
    if table is None:
        table = pair_table(chain, **build_kw)
    bounds = default_bounds(chain)
    a_range = a_range or bounds["a"]
    b_range = b_range or bounds["b"]
    a_values = np.arange(a_range[0], a_range[1] + 1e-9, da)
    b_values = np.arange(b_range[0], b_range[1] + 1e-9, db)
    res = fft_screen(table, a_values, b_values, gamma=gamma, flips=flips, n_top=n_top,
                     min_separation=min_separation, ab_symmetry=ab_symmetry)
    return res.top, res.top_energy, res
