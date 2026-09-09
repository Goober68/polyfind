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
"""
from __future__ import annotations

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
]


# ------------------------------------------------------------------ pair potential
def _pair_v(packer: CrystalPacker, r2: np.ndarray) -> np.ndarray:
    """The packer's exact pair potential (LJ + DSF Coulomb, shifted, masked) at r2."""
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
    big = packer1.reach + 20.0
    return float(packer1.energy(np.array([[big, big, 90.0, 0.0, 0.0, 0.0, 0.0]]))[0])


def intra_constants(packer2: CrystalPacker) -> np.ndarray:
    """``[E_intra(flip=0), E_intra(flip=1)]``: the constant part of a two-chain cell,
    isolated by evaluating the packer on a cell far larger than its reach.

    The two entries are measured rather than assumed to be equal, because the
    intra-chain term of the *flipped* chain is whatever the packer says it is
    (with the current :mod:`polyfind.pack` the bonded-exclusion matrices of the
    (0,0,k) column are indexed by ``+k`` for the flipped chain as well, which
    makes ``E_intra(flip=1)`` differ from ``E_intra(flip=0)`` by a large constant;
    keeping it as a measured constant means the table reproduces the packer
    either way).  Both entries are independent of the cell and of the angles.
    """
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

    ``lj_x6`` is ``lj_x ** 6`` (the LJ part is a function of ``r^2`` alone, so the
    sixth power is taken once per pair type instead of once per grid point: ``**6``
    on the (pairs, n_angle, n_angle) array is `pow`, and removing it is a third of
    the builder's arithmetic).  ``r2`` is floored at (0.05 A)^2 so that ``(x/r)^12``
    cannot overflow float32; that only affects pairs whose energy is astronomically
    large and therefore clipped to the table's cap.  Exactly coincident atoms give 0,
    as in the packer.
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
        chunk_elems: int = 131_072,
        symmetry: bool = True,
        n_threads: int | None = None,
        dtype=np.float32,
        verbose: bool = False,
    ) -> "PairTable":
        """Tabulate ``W`` with the exact pair potential of :class:`CrystalPacker`.

        Each grid point is an exact sum over all atom pairs and all z-images
        (``|k| <= K`` as in the packer).  Only the (pair, z-image) terms that can
        reach within the cutoff at that radius are evaluated (the z-window of
        section 6 of the performance review: a pair whose closest possible lateral
        approach is ``r - rad_p - rad_q`` can only reach z-images with
        ``|z_p - z_q - k c| < sqrt(rc^2 - d^2)``), and with ``symmetry=True`` the
        exchange relations

            flip 0:  W(r, a1, a2, dz) = W(r, a2 + 180, a1 + 180, -dz)
            flip 1:  W(r, a1, a2, dz) = W(r, 180 - a2, 180 - a1, dz)

        halve *both* flips.  Both are the same underlying symmetry -- swapping which
        chain sits at the origin -- but they act differently on the sum, so they are
        applied differently: exchange maps the term ``(p, q, t)`` (atom pair, z-image
        index) to ``(q, p, -t)`` for flip 0 and to ``(q, p, t)`` for flip 1.  For
        flip 0 that is a relation between the ``dz`` planes ``j`` and ``-j``, so only
        planes up to ``c/2`` are summed and the rest are filled in; for flip 1 it acts
        inside each plane, so only the ``p <= q`` half of the atom pairs is summed
        (the diagonal at half weight) and the finished block is symmetrised.
        """
        import os
        import time
        from concurrent.futures import ThreadPoolExecutor

        t_start = time.perf_counter()
        if n_threads is None:  # NumPy releases the GIL inside the ufuncs, but the loop is
            n_threads = max(1, min(4, (os.cpu_count() or 1)))  # call-bound: 4 threads measured best
        if n_angle % 2 or n_z % 2:
            raise ValueError("n_angle and n_z must be even (the exchange symmetries map grid to grid)")
        pk1 = CrystalPacker(chain, n_chains=1, cutoff=cutoff, alpha=alpha, eps_r=eps_r)
        pk2 = CrystalPacker(chain, n_chains=2, cutoff=cutoff, alpha=alpha, eps_r=eps_r)
        rc, c, K = pk1.rc, chain.c, pk1.K
        r_max = rc + 2 * chain.radius
        n_r = int(np.ceil((r_max - r_min) / dr)) + 1
        r_values = r_min + dr * np.arange(n_r)
        W = np.zeros((2, n_r, n_angle, n_angle, n_z), dtype=np.float32)

        X = np.asarray(chain.coords, dtype=float)
        rad = np.hypot(X[:, 0], X[:, 1])
        psi = np.arctan2(X[:, 1], X[:, 0])
        zc = X[:, 2]
        n = len(rad)
        ang = (2 * np.pi * np.arange(n_angle) / n_angle).astype(dtype)
        m_half = n_angle // 2
        batch = max(1, chunk_elems // (n_angle * n_angle))
        pot = tuple(np.asarray(a, dtype=float) for a in (pk1.lj_x, pk1.lj_d, pk1.qq, pk1.lj_shift))

        for flip in (0, 1):
            if symmetry and flip == 1:
                P, Q = np.triu_indices(n)  # exchange maps (p, q, t) -> (q, p, t) here
                half = P == Q
            else:
                P, Q = (v.ravel() for v in np.meshgrid(np.arange(n), np.arange(n), indexing="ij"))
                half = None
            lxb, ldb, lqb, lsb = (a[P, Q] for a in pot)
            if half is not None:  # v is linear in (lj_d, qq, lj_shift): halve the diagonal terms
                ldb, lqb, lsb = (np.where(half, 0.5 * a, a) for a in (ldb, lqb, lsb))
            lx6 = (lxb ** 6).astype(dtype)
            ld, lq, ls = (a.astype(dtype) for a in (ldb, lqb, lsb))
            rp, rq = rad[P].astype(dtype), rad[Q].astype(dtype)
            psi_q = -psi if flip else psi
            z_q = -zc if flip else zc
            pp, pq = psi[P].astype(dtype), psi_q[Q].astype(dtype)
            dzpq = zc[P] - z_q[Q]
            # (pair, t) grid with t = j + k n_z, i.e. the z offset of the image is dz_j + k c
            t_all = np.arange(-K * n_z, (K + 1) * n_z)
            zt = (dzpq[:, None] - t_all[None, :] * (c / n_z)).ravel()
            pair_i = np.repeat(np.arange(len(P)), len(t_all))
            j_f = np.tile(np.mod(t_all, n_z), len(P))
            if symmetry and flip == 0:
                sel0 = j_f <= n_z // 2
                pair_i, zt, j_f = pair_i[sel0], zt[sel0], j_f[sel0]
            gap = rp[pair_i].astype(float) + rq[pair_i].astype(float)  # closest possible xy approach is r - gap
            zt2 = zt ** 2

            def radial_block(irs, flip=flip, pp=pp, pq=pq, pair_i=pair_i, zt=zt, zt2=zt2, j_f=j_f, gap=gap,
                             lx6=lx6, ld=ld, lq=lq, ls=ls, rp=rp, rq=rq):
                for ir in irs:
                    r = r_values[ir]
                    keep = np.maximum(0.0, r - gap) ** 2 + zt2 < rc * rc
                    if not keep.any():
                        continue
                    pi, zz, jj = pair_i[keep], zt[keep].astype(dtype), j_f[keep]
                    order = np.argsort(jj, kind="stable")
                    pi, zz, jj = pi[order], zz[order], jj[order]
                    edges = np.searchsorted(jj, np.arange(n_z + 1))
                    rd = dtype(r)
                    for j in range(n_z):
                        lo_j, hi_j = int(edges[j]), int(edges[j + 1])
                        if hi_j <= lo_j:
                            continue
                        plane = np.zeros((n_angle, n_angle))
                        for s in range(lo_j, hi_j, batch):
                            sl = slice(s, min(s + batch, hi_j))
                            b = pi[sl]
                            cu, su = np.cos(pp[b][:, None] + ang), np.sin(pp[b][:, None] + ang)
                            cv, sv = np.cos(pq[b][:, None] + ang), np.sin(pq[b][:, None] + ang)
                            rpb, rqb = rp[b][:, None], rq[b][:, None]
                            ea = rpb ** 2 + rqb ** 2 + rd * rd + zz[sl][:, None] ** 2 - 2 * rd * rpb * cu
                            eb = 2 * rd * rqb * cv
                            pref = (-2 * rp[b] * rq[b])[:, None, None]
                            r2 = ea[:, :, None] + eb[:, None, :] + pref * (cu[:, :, None] * cv[:, None, :] + su[:, :, None] * sv[:, None, :])
                            v = _pot_batch(lx6[b][:, None, None], ld[b][:, None, None], lq[b][:, None, None], ls[b][:, None, None],
                                           r2, rc, alpha, pk1.dsf_shift, pk1.dsf_force, dtype)
                            plane += v.sum(axis=0)
                        W[flip, ir, :, :, j] += plane.astype(np.float32)

            if n_threads > 1 and n_r > 1:
                blocks = [np.arange(n_r)[i::n_threads] for i in range(n_threads)]
                with ThreadPoolExecutor(max_workers=n_threads) as ex:
                    list(ex.map(radial_block, blocks))
            else:
                radial_block(range(n_r))
            if symmetry and flip == 0:
                for j in range(n_z // 2 + 1, n_z):
                    src = W[0, :, :, :, (-j) % n_z]
                    W[0, :, :, :, j] = np.roll(np.swapaxes(src, 1, 2), (m_half, m_half), axis=(1, 2))
            elif symmetry and flip == 1:
                # W1 = V + S[V] with S[X][i1, i2] = X[m_half - i2, m_half - i1] and V the
                # p <= q half sum (diagonal at half weight)
                V = W[1]
                W[1] = V + np.roll(np.swapaxes(V, 1, 2)[:, ::-1, ::-1, :], (m_half + 1, m_half + 1), axis=(1, 2))
        np.minimum(W, np.float32(cap), out=W)
        dt = time.perf_counter() - t_start
        if verbose:
            print(f"table {chain.name}: {tuple(W.shape)} = {W.nbytes / 1e6:.0f} MB in {dt:.1f} s")
        return cls(chain=chain, W=W, r_min=r_min, dr=dr, n_angle=n_angle, n_z=n_z, c=c, rc=rc, cap=cap,
                   e_intra=intra_constants(pk2), build_time=dt)

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
# A table costs 5-25 s to build on the screen grid of :data:`polyfind.pack.SCREEN_TABLE`
# and 35-150 s on this module's finer default, and depends only on the chain geometry
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
    ``n_angle``, ``n_z``, ``cap``, ``dtype``).  ``symmetry``, ``chunk_elems``,
    ``n_threads`` and ``verbose`` only change how the same numbers are produced and are
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


def _shift_multiplier(theta_deg: float, n_angle: int) -> np.ndarray:
    """DFT multiplier of a *linearly interpolated* rotation of the angle axis by
    ``theta``: the exact Fourier image of what :meth:`PairTable.interpolate` does."""
    g = (theta_deg % 360.0) * (n_angle / 360.0)
    q = int(np.floor(g))
    f = g - q
    m = np.fft.fftfreq(n_angle) * n_angle  # 0, 1, ..., -1
    ph = np.exp(-2j * np.pi * m * q / n_angle)
    return (ph * ((1.0 - f) + f * np.exp(-2j * np.pi * m / n_angle))).astype(np.complex64)


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
    import os
    import time
    from concurrent.futures import ThreadPoolExecutor

    t0 = time.perf_counter()
    if n_threads is None:
        n_threads = max(1, min(4, (os.cpu_count() or 1)))
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

    for fi, flip in enumerate(flips):
        Wh = table.fourier(flip)
        sgn2 = _self_sign(flip)
        e0 = float(table.e_intra[flip])

        def do_rows(ias, fi=fi, Wh=Wh, sgn2=sgn2, e0=e0):
            buf = np.empty((batch, na, na, nz2), dtype=np.complex64)
            meta = []

            def flush():
                if not meta:
                    return
                E = np.fft.irfftn(buf[: len(meta)], s=(na, na, nz), axes=(1, 2, 3))
                for k, (ia, ib, self1, self2) in enumerate(meta):
                    Ek = E[k] + self1[:, None, None] + self2[None, :, None] + e0
                    i1, i2, j = np.unravel_index(int(np.argmin(Ek)), Ek.shape)
                    energy[ia, ib, fi] = Ek[i1, i2, j]
                    best[ia, ib, fi] = (phi[i1], phi[i2], c * j / nz)
                meta.clear()

            for ia in ias:
                a = a_values[ia]
                for ib, b in enumerate(b_values):
                    if ab_symmetry and b < a - 1e-9:
                        continue
                    rs, ts = lattice_sites(a, b, gamma, r_max, centred=False)
                    rc_, tc = lattice_sites(a, b, gamma, r_max, centred=True)
                    if min(rs.min(initial=np.inf), rc_.min(initial=np.inf)) < table.r_min:
                        continue  # no table data there; left at +inf (see the docstring)
                    da = phi[:, None] - ts[None, :]
                    self1 = 0.5 * table.interpolate(rs[None, :], da, da, 0.0, 0).sum(axis=1)
                    self2 = self1 if sgn2 > 0 else 0.5 * table.interpolate(rs[None, :], -da, -da, 0.0, 0).sum(axis=1)
                    i0, tr = table.radial_weights(rc_)
                    acc = buf[len(meta)]
                    acc[...] = 0
                    # sites at the same radius (a rectangular lattice gives them in fours)
                    # share the coefficient array: sum their rotation multipliers first
                    keys = np.round(rc_, 9)
                    for key in np.unique(keys):
                        grp = keys == key
                        LL = np.zeros((na, na), dtype=np.complex64)
                        for t_i in tc[grp]:
                            L = _shift_multiplier(t_i, na)
                            LL += L[:, None] * L[None, :]
                        r_i, w_i = int(i0[grp][0]), float(tr[grp][0])
                        coef = Wh[r_i] if w_i == 0.0 else (1.0 - w_i) * Wh[r_i] + w_i * Wh[r_i + 1]
                        acc += LL[:, :, None] * coef
                    meta.append((ia, ib, self1, self2))
                    if len(meta) >= batch:
                        flush()
            flush()

        if n_threads > 1 and len(a_values) > 1:
            blocks = [range(i, len(a_values), n_threads) for i in range(n_threads)]
            with ThreadPoolExecutor(max_workers=n_threads) as ex:
                list(ex.map(do_rows, blocks))
        else:
            do_rows(range(len(a_values)))

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
