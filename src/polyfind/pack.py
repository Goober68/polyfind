"""Crystal packing of rigid helical chains with a handful of degrees of freedom.

Given the periodic chain produced by :func:`periodic_chain` (axis along z,
repeat ``c``), a two-chain cell is described by

    a, b, gamma      lattice in the plane perpendicular to the chains
    phi1, phi2       setting angles (rotation of each chain about its own axis,
                     measured from the a axis to the chain's principal plane)
    dz               relative shift of chain 2 along z
    flip             chain 2 parallel (0) or antiparallel (1)

with chain 2 centred at (1/2, 1/2) in the ab plane (as in PE, alpha-, beta- and
gamma-PVDF).  The lattice energy per cell is the sum of Lennard-Jones (UFF) and
damped-shifted-force Coulomb terms over all periodic images within a cutoff,
*including* the intra-chain energy (with bonded exclusions across the periodic
boundary) so that polymorphs with different chain conformations are comparable.

The energy kernel is batched over configurations on the active array backend:
the coarse search evaluates thousands of candidate cells in one call (a GPU
does this in a fraction of a second), and only the best few are polished with
a local optimiser.  Per configuration only the lateral images whose chain axes
can lie within the cutoff are evaluated.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import minimize

from . import backend as bk
from .chain import build_chain
from .forcefield import COULOMB, erfc_approx
from .helix import HelixParams, helix_parameters, kabsch, rotation_to_z
from .polymers import Polymer, RISStates, lj_params

MASS = {"C": 12.011, "H": 1.008, "F": 18.998, "Cl": 35.45}
_TOPO_CACHE: dict = {}


# ------------------------------------------------------------------ periodic chain
@dataclass
class PeriodicChain:
    polymer: Polymer
    name: str
    elements: list[str]
    coords: np.ndarray  # (n, 3), chain axis = z through the origin
    charges: np.ndarray
    c: float  # repeat along z (A)
    n_monomers: int
    helix: HelixParams
    dihedrals: np.ndarray  # the dihedrals of one crystallographic repeat (deg)
    same_site_scale: dict  # k (z-image shift) -> (n, n) nonbonded scale (0 excluded, 0.5 for 1-4, 1 otherwise)
    backbone: np.ndarray  # indices of backbone atoms within the block
    rotation_error: float = 0.0  # residual rotation of the repeat transform (deg); 0 for a true periodic chain

    @property
    def n_atoms(self) -> int:
        return len(self.elements)

    @property
    def mass(self) -> float:
        return float(sum(MASS[e] for e in self.elements))

    @property
    def radius(self) -> float:
        return float(np.linalg.norm(self.coords[:, :2], axis=1).max())


def _block_atoms(full, nb_block: int, kb: int):
    bb = [int(full.backbone[k]) for k in range(nb_block * (kb + 2), nb_block * (kb + 3))]
    atoms = [i for k in bb for i in [k] + list(full.subs_of[k])]
    return bb, atoms


def _same_site_scales(polymer: Polymer, full, nb_block: int, scale14: float) -> dict:
    """Nonbonded scale matrices between block 0 and blocks k = -2..2 (graph distances)."""
    key = (polymer.name, nb_block, full.n_dihedrals, scale14)
    if key in _TOPO_CACHE:
        return _TOPO_CACHE[key]
    _, atoms0 = _block_atoms(full, nb_block, 0)
    n = len(atoms0)
    adj = [[] for _ in range(full.n_atoms)]
    for a, b in full.bonds:
        adj[a].append(b)
        adj[b].append(a)
    dist = {}
    for s in atoms0:
        d = {s: 0}
        frontier = [s]
        for step in (1, 2, 3):
            nxt = []
            for u in frontier:
                for v in adj[u]:
                    if v not in d:
                        d[v] = step
                        nxt.append(v)
            frontier = nxt
        dist[s] = d
    out = {}
    for kshift in (-2, -1, 0, 1, 2):
        _, atomsk = _block_atoms(full, nb_block, kshift)
        S = np.ones((n, n))
        for p, ip in enumerate(atoms0):
            for q, iq in enumerate(atomsk):
                gd = dist[ip].get(iq, 99)
                if gd <= 2:
                    S[p, q] = 0.0
                elif gd == 3:
                    S[p, q] = scale14
        out[kshift] = S
    _TOPO_CACHE[key] = out
    return out


def _orient_block(coords: np.ndarray, bb_local: list[int]) -> np.ndarray:
    """z from 0; cross-section principal axis (zigzag plane for planar chains) along x."""
    coords = coords - np.array([0, 0, coords[:, 2].min()])
    xy = coords[bb_local][:, :2]
    if len(xy) > 1:
        w, v = np.linalg.eigh(np.cov(xy.T))
        ux = v[:, np.argmax(w)]
        ang = np.arctan2(ux[1], ux[0])
        ca, sa = np.cos(-ang), np.sin(-ang)
        coords = coords @ np.array([[ca, -sa, 0], [sa, ca, 0], [0, 0, 1]]).T
    return coords


def _assemble(polymer, name, dih_block, helix, axis, axis_point, c, scale14, rotation_error=0.0) -> PeriodicChain:
    nb_block = len(dih_block)
    full = build_chain(polymer, np.tile(dih_block, 5), cap=False)
    R = rotation_to_z(axis)
    X = (full.coords - axis_point) @ R.T
    bb0, atoms0 = _block_atoms(full, nb_block, 0)
    coords = X[atoms0]
    bb_local = [atoms0.index(k) for k in bb0]
    coords = _orient_block(coords, bb_local)
    return PeriodicChain(
        polymer=polymer,
        name=name,
        elements=[full.elements[i] for i in atoms0],
        coords=coords,
        charges=full.charges[atoms0],
        c=c,
        n_monomers=nb_block // polymer.bonds_per_repeat,
        helix=helix,
        dihedrals=np.asarray(dih_block, dtype=float),
        same_site_scale=_same_site_scales(polymer, full, nb_block, scale14),
        backbone=np.array(bb_local),
        rotation_error=rotation_error,
    )


def periodic_chain(polymer: Polymer, seq, states: RISStates, helix: HelixParams | None = None, scale14: float = 0.5) -> PeriodicChain:
    """One crystallographic repeat of the infinite chain generated by periodic ``seq``
    at the ideal state angles."""
    seq = [int(s) for s in seq]
    P = len(seq)
    if helix is None:
        helix = helix_parameters(polymer, np.array(seq), states)
    if helix.c is None:
        raise ValueError("chain is not commensurate; cannot build a periodic repeat")
    m = helix.periods_per_repeat
    per = np.array([states.angles[s] for s in seq])
    dih_block = np.tile(per, m)
    name = "".join(states.names[s] for s in seq)
    ch = _assemble(polymer, name, dih_block, helix, helix.axis, helix.axis_point, m * helix.rise_per_period, scale14)
    # sanity: block + c must coincide with the next block
    full = build_chain(polymer, np.tile(dih_block, 5), cap=False)
    R = rotation_to_z(helix.axis)
    X = (full.coords - helix.axis_point) @ R.T
    _, a0 = _block_atoms(full, len(dih_block), 0)
    _, a1 = _block_atoms(full, len(dih_block), 1)
    err = np.abs(X[a1] - X[a0] - np.array([0, 0, ch.c])).max()
    if err > 1e-3:
        raise RuntimeError(f"periodic repeat inconsistent (max deviation {err:.3g} A)")
    return ch


def periodic_chain_from_torsions(polymer: Polymer, name: str, torsions, states: RISStates | None = None, scale14: float = 0.5) -> PeriodicChain:
    """Periodic chain whose crystallographic repeat has the given (arbitrary) torsions.

    The transform between consecutive repeats is forced to a translation along its
    own axis; its residual rotation angle is reported in ``rotation_error`` (deg) so a
    refinement can penalise non-commensurate torsion sets.
    """
    tors = np.asarray(torsions, dtype=float)
    nb = len(tors)
    from .chain import build_backbone

    bb = build_backbone(polymer, np.tile(tors, 5), xp=np)
    X, Y = bb[2 * nb : 3 * nb], bb[3 * nb : 4 * nb]
    if nb < 3:  # need >= 3 points for a rigid fit; use two blocks
        X, Y = bb[2 * nb : 2 * nb + 4], bb[3 * nb : 3 * nb + 4]
    R, t = kabsch(X, Y)
    ang = float(np.degrees(np.arccos(np.clip((np.trace(R) - 1) / 2, -1, 1))))
    ang = min(ang, 360.0 - ang)
    c = float(np.linalg.norm(t))
    axis = t / c
    full = build_chain(polymer, np.tile(tors, 5), cap=False)
    _, a0 = _block_atoms(full, nb, 0)
    point = full.coords[a0].mean(axis=0)
    helix = HelixParams(
        sequence=name, period_bonds=nb, monomers_per_period=nb // polymer.bonds_per_repeat,
        rotation_per_period=ang, rise_per_period=c, rise_per_bond=c / nb, periods_per_repeat=1,
        turns_per_repeat=0, c=c, radius_backbone=0.0, radius_all=0.0, axis=axis, axis_point=point, label="refined",
    )
    ch = _assemble(polymer, name, tors, helix, axis, point, c, scale14, rotation_error=ang)
    ch.helix.radius_all = ch.radius
    ch.helix.radius_backbone = float(np.linalg.norm(ch.coords[ch.backbone][:, :2], axis=1).max())
    return ch


# ------------------------------------------------------------------ energy kernel
@dataclass
class PackResult:
    chain: str
    energy_per_cell: float
    energy_per_monomer: float
    a: float
    b: float
    gamma: float
    c: float
    phi1: float
    phi2: float
    dz: float
    flip: int
    n_chains: int
    density: float
    dihedrals: np.ndarray | None = None
    cell_coords: np.ndarray | None = None
    cell_elements: list[str] | None = None

    @property
    def params(self) -> np.ndarray:
        return np.array([self.a, self.b, self.gamma, self.phi1, self.phi2, self.dz, self.flip], dtype=float)

    def row(self) -> str:
        return (
            f"{self.chain:<14s} E/mon={self.energy_per_monomer:8.3f}  a={self.a:5.2f} b={self.b:5.2f} c={self.c:5.2f} "
            f"gamma={self.gamma:5.1f}  phi=({self.phi1:5.1f},{self.phi2:5.1f}) dz={self.dz:4.2f} {'anti' if self.flip else 'para'}  rho={self.density:5.3f}"
        )


class CrystalPacker:
    """Batched lattice-energy evaluator for a rigid periodic chain."""

    PARAMS = ("a", "b", "gamma", "phi1", "phi2", "dz", "flip")

    def __init__(self, chain: PeriodicChain, n_chains: int = 2, cutoff: float = 8.0, alpha: float = 0.2, eps_r: float = 1.0, torsion=(1.3, -0.05, 2.5), xp=None):
        self.chain = chain
        self.n_chains = n_chains
        self.rc = cutoff
        self.alpha = alpha
        self.xp = xp or bk.get_backend()
        xp = self.xp
        dt = bk.float_dtype()
        n = chain.n_atoms
        x, d = lj_params(chain.elements)
        xs = np.sqrt(np.outer(x, x))
        ds = np.sqrt(np.outer(d, d))
        qq = np.outer(chain.charges, chain.charges) * COULOMB / eps_r
        self.lj_x = xp.asarray(np.tile(xs, (n_chains, n_chains)), dtype=dt)
        self.lj_d = xp.asarray(np.tile(ds, (n_chains, n_chains)), dtype=dt)
        self.qq = xp.asarray(np.tile(qq, (n_chains, n_chains)), dtype=dt)
        rc = cutoff
        erc = float(erfc_approx(np.array(alpha * rc)))
        self.dsf_shift = erc / rc
        self.dsf_force = erc / rc ** 2 + 2 * alpha / np.sqrt(np.pi) * np.exp(-(alpha * rc) ** 2) / rc
        s6 = (xs / rc) ** 6
        self.lj_shift = xp.asarray(np.tile(ds * (s6 * s6 - 2 * s6), (n_chains, n_chains)), dtype=dt)
        V1, V2, V3 = torsion
        phi = np.deg2rad(chain.dihedrals)
        self.e_torsion = float((0.5 * (V1 * (1 + np.cos(phi)) + V2 * (1 - np.cos(2 * phi)) + V3 * (1 + np.cos(3 * phi)))).sum()) * n_chains
        N = n * n_chains
        self.same_site = {}
        for k, S in chain.same_site_scale.items():
            M = np.ones((N, N))
            for s in range(n_chains):
                M[s * n : (s + 1) * n, s * n : (s + 1) * n] = S
            self.same_site[k] = M
        self.n, self.N = n, N
        self.K = int(np.ceil(self.rc / self.chain.c)) + 1
        self.X0 = xp.asarray(chain.coords, dtype=dt)
        self.reach = self.rc + 2 * chain.radius + 0.5  # axis-axis distance beyond which no atom pair is within the cutoff
        # scale tensor for the (0,0,k) column (index 0..2K) then ones for the other lateral images
        ks = list(range(-self.K, self.K + 1))
        col = np.stack([self.same_site.get(k, np.ones((N, N))) for k in ks])
        self._scale_col = xp.asarray(col, dtype=dt)  # (2K+1, N, N)
        self._dt = dt

    # --- geometry --------------------------------------------------------------
    def _place(self, params):
        """params (M, 7) -> primary coordinates (M, N, 3), lattice vectors (M, 3, 3)."""
        xp = self.xp
        p = xp.asarray(params, dtype=self._dt)
        M = p.shape[0]
        a, b, gam, phi1, phi2, dz, flip = (p[:, i] for i in range(7))
        g = xp.deg2rad(gam)
        z = xp.zeros(M, dtype=self._dt)
        avec = xp.stack([a, z, z], axis=1)
        bvec = xp.stack([b * xp.cos(g), b * xp.sin(g), z], axis=1)
        cvec = xp.stack([z, z, xp.full(M, self.chain.c, dtype=self._dt)], axis=1)

        def rotz(X, ang):
            ca, sa = xp.cos(xp.deg2rad(ang))[:, None], xp.sin(xp.deg2rad(ang))[:, None]
            return xp.stack([ca * X[..., 0] - sa * X[..., 1], sa * X[..., 0] + ca * X[..., 1], X[..., 2]], axis=-1)

        X = xp.broadcast_to(self.X0, (M,) + self.X0.shape)
        chains = [rotz(X, phi1)]
        if self.n_chains == 2:
            Xf = xp.where(flip[:, None, None] > 0.5, xp.stack([X[..., 0], -X[..., 1], -X[..., 2]], axis=-1), X)
            X2 = rotz(Xf, phi2) + (0.5 * avec + 0.5 * bvec)[:, None, :] + xp.stack([z, z, dz], axis=1)[:, None, :]
            chains.append(X2)
        return xp.concatenate(chains, axis=1), xp.stack([avec, bvec, cvec], axis=1)

    def _select_images(self, params: np.ndarray):
        """Per configuration, the lateral images (i, j) whose chain axes can be within reach.

        Returns ijk (M, I, 3) with the (0,0,k) column first (so a shared scale tensor
        applies), padded with far images that fall outside the cutoff.
        """
        M = params.shape[0]
        a, b, g = params[:, 0], params[:, 1], np.deg2rad(params[:, 2])
        a_min = float(min(a.min(), (b * np.sin(g)).min()))
        R = int(np.ceil(self.reach / max(a_min, 1e-3))) + 1
        grid = [(0, 0)] + [(i, j) for i in range(-R, R + 1) for j in range(-R, R + 1) if (i, j) != (0, 0)]
        G = np.array(grid, dtype=float)  # (Gn, 2)
        av = np.stack([a, np.zeros(M)], axis=1)
        bv = np.stack([b * np.cos(g), b * np.sin(g)], axis=1)
        pos = G[None, :, 0, None] * av[:, None, :] + G[None, :, 1, None] * bv[:, None, :]  # (M, Gn, 2)
        deltas = [np.zeros((M, 2))]
        if self.n_chains == 2:
            o = 0.5 * av + 0.5 * bv
            deltas += [o, -o]
        dmin = np.min(np.stack([np.linalg.norm(pos + d[:, None, :], axis=2) for d in deltas]), axis=0)  # (M, Gn)
        keep = dmin < self.reach
        keep[:, 0] = True
        counts = keep.sum(axis=1)
        L = int(counts.max())
        far = np.array([3 * R + 5, 3 * R + 5], dtype=float)
        lat = np.empty((M, L, 2))
        for m in range(M):
            sel = G[keep[m]]
            lat[m, : len(sel)] = sel
            lat[m, len(sel) :] = far
        ks = np.arange(-self.K, self.K + 1, dtype=float)
        ijk = np.concatenate([np.repeat(lat[:, :, None, :], len(ks), axis=2), np.broadcast_to(ks[None, None, :, None], (M, L, len(ks), 1))], axis=3)
        return ijk.reshape(M, L * len(ks), 3), L

    # --- energy ------------------------------------------------------------------
    def energy(self, params, chunk_elems: int | None = None) -> np.ndarray:
        """Lattice energy per cell (kcal/mol) for each row of params (M, 7).

        Configurations are sorted by cell area (so chunks share similar image counts and
        little padding is wasted) and processed in chunks of ``chunk_elems`` array
        elements: small on CPU (cache-resident), large on GPU (fewer kernel launches).
        """
        xp = self.xp
        params = np.atleast_2d(np.asarray(params, dtype=float))
        if chunk_elems is None:
            chunk_elems = 60_000_000 if bk.device_name() == "cuda" else 500_000
        order = np.argsort(params[:, 0] * params[:, 1] * np.sin(np.deg2rad(params[:, 2])))
        params = params[order]
        out = np.empty(params.shape[0])
        N = self.N
        nk = 2 * self.K + 1
        est_I = 30 * nk
        m_chunk = max(1, chunk_elems // (est_I * N * N))
        for s in range(0, params.shape[0], m_chunk):
            sub = params[s : s + m_chunk]
            ijk_np, L = self._select_images(sub)
            ijk = xp.asarray(ijk_np, dtype=self._dt)
            P, lat = self._place(sub)
            D = P[:, :, None, :] - P[:, None, :, :]  # (M, N, N, 3)
            shift = xp.einsum("mic,mcd->mid", ijk, lat)  # (M, I, 3)
            dx = D[:, None, :, :, 0] - shift[:, :, None, None, 0]
            dy = D[:, None, :, :, 1] - shift[:, :, None, None, 1]
            dzz = D[:, None, :, :, 2] - shift[:, :, None, None, 2]
            r2 = dx * dx + dy * dy + dzz * dzz  # (M, I, N, N)
            mask = (r2 < self.rc ** 2) & (r2 > 1e-8)
            rr = xp.sqrt(xp.where(mask, r2, 1.0))
            s6 = (self.lj_x / rr) ** 6
            v = self.lj_d * (s6 * s6 - 2 * s6) - self.lj_shift
            v = v + self.qq * (erfc_approx(self.alpha * rr, xp) / rr - self.dsf_shift + self.dsf_force * (rr - self.rc))
            v = xp.where(mask, v, 0.0)
            # scale: first nk images are the (0,0,k) column
            e0 = (v[:, :nk] * self._scale_col[None]).sum(axis=(1, 2, 3))
            e1 = v[:, nk:].sum(axis=(1, 2, 3))
            out[s : s + m_chunk] = bk.to_numpy(0.5 * (e0 + e1))
        res = np.empty_like(out)
        res[order] = out
        return res + self.e_torsion

    def energy_per_monomer(self, params) -> np.ndarray:
        return self.energy(params) / (self.n_chains * self.chain.n_monomers)

    # --- helpers ---------------------------------------------------------------
    def density(self, a, b, gamma) -> float:
        V = a * b * np.sin(np.deg2rad(gamma)) * self.chain.c
        return self.n_chains * self.chain.mass / (V * 0.602214)

    def result(self, params) -> PackResult:
        params = np.asarray(params, dtype=float)
        e = float(self.energy(params[None])[0])
        P, _ = self._place(params[None])
        a, b, gam, phi1, phi2, dz, flip = params
        return PackResult(
            chain=self.chain.name,
            energy_per_cell=e,
            energy_per_monomer=e / (self.n_chains * self.chain.n_monomers),
            a=float(a), b=float(b), gamma=float(gam), c=self.chain.c,
            phi1=float(phi1 % 360), phi2=float(phi2 % 360), dz=float(dz % self.chain.c), flip=int(round(flip)),
            n_chains=self.n_chains,
            density=self.density(a, b, gam),
            dihedrals=self.chain.dihedrals.copy(),
            cell_coords=bk.to_numpy(P[0]),
            cell_elements=list(self.chain.elements) * self.n_chains,
        )


# ------------------------------------------------------------------ search
def default_bounds(chain: PeriodicChain, gamma_free: bool = False) -> dict:
    r = chain.radius + 1.5
    return {
        "a": (1.4 * r, 3.2 * r),
        "b": (1.4 * r, 3.2 * r),
        "gamma": (60.0, 120.0) if gamma_free else (90.0, 90.0),
        "phi1": (0.0, 360.0),
        "phi2": (0.0, 360.0),
        "dz": (0.0, chain.c),
    }


def polish(packer: CrystalPacker, params: np.ndarray, lo: np.ndarray, hi: np.ndarray, free: list[int], maxfev: int = 1500) -> np.ndarray:
    """Nelder-Mead on the continuous cell parameters (flip fixed), bounds by penalty."""
    x0 = np.asarray(params, dtype=float).copy()

    def objective(xf):
        x = x0.copy()
        x[free] = xf
        pen = 0.0
        for i in free:
            if x[i] < lo[i]:
                pen += 100 * (lo[i] - x[i]) ** 2
                x[i] = lo[i]
            elif x[i] > hi[i]:
                pen += 100 * (x[i] - hi[i]) ** 2
                x[i] = hi[i]
        return float(packer.energy(x[None])[0]) + pen

    res = minimize(objective, x0[free], method="Nelder-Mead", options={"xatol": 1e-3, "fatol": 1e-4, "maxfev": maxfev})
    x = x0.copy()
    x[free] = res.x
    x[:6] = np.clip(x[:6], lo, hi)
    return x


def pack(
    chain: PeriodicChain,
    n_chains: int = 2,
    n_random: int = 3000,
    n_refine: int = 5,
    bounds: dict | None = None,
    gamma_free: bool = False,
    flips=(0, 1),
    rng=None,
    cutoff: float = 8.0,
    eps_r: float = 1.0,
    maxfev: int = 1500,
    verbose: bool = False,
) -> list[PackResult]:
    """Coarse batched random search followed by Nelder-Mead polishing of the best cells."""
    rng = rng or np.random.default_rng(0)
    packer = CrystalPacker(chain, n_chains=n_chains, cutoff=cutoff, eps_r=eps_r)
    bounds = bounds or default_bounds(chain, gamma_free)
    keys = ["a", "b", "gamma", "phi1", "phi2", "dz"]
    lo = np.array([bounds[k][0] for k in keys])
    hi = np.array([bounds[k][1] for k in keys])
    flips = list(flips) if n_chains == 2 else [0]
    cont = lo + (hi - lo) * rng.random((n_random, 6))
    if n_chains == 1:
        cont[:, 4] = 0.0
        cont[:, 5] = 0.0
    Es, Ps = [], []
    for f in flips:
        params = np.concatenate([cont, np.full((n_random, 1), float(f))], axis=1)
        Es.append(packer.energy(params))
        Ps.append(params)
    E = np.concatenate(Es)
    params = np.concatenate(Ps)
    order = np.argsort(E)
    starts = []
    for idx in order:
        p = params[idx]
        if all(abs(p[0] - s[0]) + abs(p[1] - s[1]) > 0.5 or p[6] != s[6] for s in starts):
            starts.append(p)
        if len(starts) >= n_refine:
            break
    if verbose:
        print(f"coarse search: {len(E)} cells, best E/cell = {E[order[0]]:.3f}")
    free = [i for i in range(6) if hi[i] > lo[i] and not (n_chains == 1 and i in (4, 5))]
    results = [packer.result(polish(packer, p, lo, hi, free, maxfev)) for p in starts]
    results.sort(key=lambda r: r.energy_per_cell)
    uniq = []
    for r in results:
        if all(abs(r.energy_per_cell - u.energy_per_cell) > 1e-3 for u in uniq):
            uniq.append(r)
    return uniq


def to_cif(res: PackResult, title: str = "polyfind") -> str:
    """P1 CIF of the packed cell (fractional coordinates)."""
    a, b, c, gam = res.a, res.b, res.c, res.gamma
    g = np.deg2rad(gam)
    lat = np.array([[a, 0, 0], [b * np.cos(g), b * np.sin(g), 0], [0, 0, c]])
    frac = np.linalg.solve(lat.T, res.cell_coords.T).T % 1.0
    lines = [
        f"data_{title}",
        f"_cell_length_a {a:.4f}",
        f"_cell_length_b {b:.4f}",
        f"_cell_length_c {c:.4f}",
        "_cell_angle_alpha 90",
        "_cell_angle_beta 90",
        f"_cell_angle_gamma {gam:.3f}",
        "_symmetry_space_group_name_H-M 'P 1'",
        "loop_",
        "_atom_site_label",
        "_atom_site_type_symbol",
        "_atom_site_fract_x",
        "_atom_site_fract_y",
        "_atom_site_fract_z",
    ]
    for i, (e, f) in enumerate(zip(res.cell_elements, frac)):
        lines.append(f"{e}{i + 1} {e} {f[0]:.5f} {f[1]:.5f} {f[2]:.5f}")
    return "\n".join(lines) + "\n"
