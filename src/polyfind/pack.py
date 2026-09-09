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

# Working-set target for one chunk of the pair kernel on the CPU.  The kernel is
# memory-bound, so the fastest chunk is the one whose temporaries stay resident:
# measured on PE/alpha/gamma, one configuration per chunk beats ten by 2-4x.  (On a
# GPU the opposite holds and the default is 60M elements.)
CPU_CHUNK_ELEMS = 60_000


def _erfc_pos(x, xp=np):
    """:func:`polyfind.forcefield.erfc_approx` restricted to x >= 0 (the only case here).

    Identical values for non-negative arguments; skips the ``abs`` and the ``where``,
    which on the kernel's largest array is a measurable fraction of the work.
    """
    t = 1.0 / (1.0 + 0.3275911 * x)
    poly = t * (0.254829592 + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))))
    return poly * xp.exp(-x * x)


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


def _align_about_z(coords: np.ndarray, ref: np.ndarray) -> np.ndarray:
    """Rotate ``coords`` about z and shift it along z to best superimpose it on ``ref``.

    The frame chosen by :func:`_orient_block` is a principal axis of the cross section,
    whose sign (and, for a nearly isotropic cross section, whose choice) flips
    discontinuously as the torsions vary.  A rotation about z combined with a shift
    along z is an exact symmetry of the lattice energy when absorbed into
    ``phi1, phi2, dz``, so removing it makes the map torsions -> coordinates smooth,
    which is what a finite-difference gradient needs.
    """
    cw = float((coords[:, 0] * ref[:, 0] + coords[:, 1] * ref[:, 1]).sum())
    sw = float((coords[:, 0] * ref[:, 1] - coords[:, 1] * ref[:, 0]).sum())
    ang = np.arctan2(sw, cw)
    ca, sa = np.cos(ang), np.sin(ang)
    out = np.empty_like(coords)
    out[:, 0] = ca * coords[:, 0] - sa * coords[:, 1]
    out[:, 1] = sa * coords[:, 0] + ca * coords[:, 1]
    out[:, 2] = coords[:, 2] + (ref[:, 2] - coords[:, 2]).mean()
    return out


def _assemble(polymer, name, dih_block, helix, axis, axis_point, c, scale14, rotation_error=0.0, full=None) -> PeriodicChain:
    nb_block = len(dih_block)
    if full is None:
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
    full = build_chain(polymer, np.tile(dih_block, 5), cap=False)
    ch = _assemble(polymer, name, dih_block, helix, helix.axis, helix.axis_point, m * helix.rise_per_period, scale14, full=full)
    # sanity: block + c must coincide with the next block
    R = rotation_to_z(helix.axis)
    X = (full.coords - helix.axis_point) @ R.T
    _, a0 = _block_atoms(full, len(dih_block), 0)
    _, a1 = _block_atoms(full, len(dih_block), 1)
    err = np.abs(X[a1] - X[a0] - np.array([0, 0, ch.c])).max()
    if err > 1e-3:
        raise RuntimeError(f"periodic repeat inconsistent (max deviation {err:.3g} A)")
    return ch


_TEMPLATE_CACHE: dict = {}


def _template(polymer: Polymer, nb: int):
    """A five-block oligomer of the right topology; only its bonding is used."""
    key = (polymer.name, nb)
    tpl = _TEMPLATE_CACHE.get(key)
    if tpl is None:
        tpl = _TEMPLATE_CACHE[key] = build_chain(polymer, np.zeros(5 * nb), cap=False)
    return tpl


def _batch_block_coords(polymer: Polymer, tors_batch: np.ndarray) -> np.ndarray:
    """(M, nb) repeat torsions -> (M, n_full, 3), atom order as ``build_chain(cap=False)``.

    One batched NeRF pass for the whole batch: the per-atom Python loop is paid once for
    the batch instead of once per chain, which is what makes a gradient over torsions
    affordable.
    """
    from .chain import build_backbone, nerf, substituent_positions

    M, nb = tors_batch.shape
    dih = np.tile(tors_batch, (1, 5))
    N = dih.shape[1]
    B, L = polymer.bonds_per_repeat, polymer.bond_length
    bb = build_backbone(polymer, dih, xp=np)  # (M, N+3, 3)
    v0 = nerf(bb[:, 2], bb[:, 1], bb[:, 0], L, polymer.backbone[0].backbone_angle, 180.0, xp=np)
    v1 = nerf(bb[:, N], bb[:, N + 1], bb[:, N + 2], L, polymer.backbone[(N + 2) % B].backbone_angle, 180.0, xp=np)
    bb_ext = np.concatenate([v0[:, None], bb, v1[:, None]], axis=1)  # (M, N+5, 3)
    out = np.empty((M, 3 * (N + 3), 3))
    for k in range(N + 3):
        spec = polymer.backbone[k % B]
        x = bb_ext[:, k + 1]
        s1, s2 = substituent_positions(bb_ext[:, k], x, bb_ext[:, k + 2], spec.sub_bond, spec.sub_angle, xp=np)
        out[:, 3 * k], out[:, 3 * k + 1], out[:, 3 * k + 2] = x, s1, s2
    return out


def repeat_chains_from_torsions(polymer: Polymer, name: str, torsions_batch, scale14: float = 0.5, align_to: "PeriodicChain | np.ndarray | None" = None) -> list[PeriodicChain]:
    """:func:`periodic_chain_from_torsions` for a whole batch of torsion sets at once."""
    tors_batch = np.atleast_2d(np.asarray(torsions_batch, dtype=float))
    M, nb = tors_batch.shape
    tpl = _template(polymer, nb)
    coords_full = _batch_block_coords(polymer, tors_batch)
    bb0, atoms0 = _block_atoms(tpl, nb, 0)
    bb_local = [atoms0.index(k) for k in bb0]
    scales = _same_site_scales(polymer, tpl, nb, scale14)
    elements = [tpl.elements[i] for i in atoms0]
    charges = tpl.charges[atoms0]
    ref = None
    if align_to is not None:
        ref = align_to.coords if isinstance(align_to, PeriodicChain) else np.asarray(align_to, dtype=float)
    out = []
    for m in range(M):
        C = coords_full[m]
        bbc = C[tpl.backbone]
        X, Y = bbc[2 * nb : 3 * nb], bbc[3 * nb : 4 * nb]
        if nb < 3:  # need >= 3 points for a rigid fit; use two blocks
            X, Y = bbc[2 * nb : 2 * nb + 4], bbc[3 * nb : 3 * nb + 4]
        R, t = kabsch(X, Y)
        ang = float(np.degrees(np.arccos(np.clip((np.trace(R) - 1) / 2, -1, 1))))
        ang = min(ang, 360.0 - ang)
        c = float(np.linalg.norm(t))
        axis = t / c
        point = C[atoms0].mean(axis=0)
        coords = _orient_block((C[atoms0] - point) @ rotation_to_z(axis).T, bb_local)
        if ref is not None:
            coords = _align_about_z(coords, ref)
        tors = tors_batch[m]
        helix = HelixParams(
            sequence=name, period_bonds=nb, monomers_per_period=nb // polymer.bonds_per_repeat,
            rotation_per_period=ang, rise_per_period=c, rise_per_bond=c / nb, periods_per_repeat=1,
            turns_per_repeat=0, c=c, radius_backbone=0.0, radius_all=0.0, axis=axis, axis_point=point, label="refined",
        )
        ch = PeriodicChain(
            polymer=polymer, name=name, elements=list(elements), coords=coords, charges=charges.copy(),
            c=c, n_monomers=nb // polymer.bonds_per_repeat, helix=helix, dihedrals=tors.copy(),
            same_site_scale=scales, backbone=np.array(bb_local), rotation_error=ang,
        )
        ch.helix.radius_all = ch.radius
        ch.helix.radius_backbone = float(np.linalg.norm(coords[ch.backbone][:, :2], axis=1).max())
        out.append(ch)
    return out


def periodic_chain_from_torsions(polymer: Polymer, name: str, torsions, states: RISStates | None = None, scale14: float = 0.5, align_to: "PeriodicChain | np.ndarray | None" = None) -> PeriodicChain:
    """Periodic chain whose crystallographic repeat has the given (arbitrary) torsions.

    The transform between consecutive repeats is forced to a translation along its
    own axis; its residual rotation angle is reported in ``rotation_error`` (deg) so a
    refinement can penalise non-commensurate torsion sets.

    ``align_to`` (a :class:`PeriodicChain` or an ``(n, 3)`` coordinate array with the
    same atom order) rotates the finished repeat about z and shifts it along z onto
    that reference.  Both are exact symmetries of the lattice energy (absorbed by
    ``phi1, phi2, dz``), and using them makes the chain vary smoothly with the
    torsions, which a finite-difference gradient requires; see :func:`_align_about_z`.
    """
    return repeat_chains_from_torsions(polymer, name, np.asarray(torsions, dtype=float)[None], scale14=scale14, align_to=align_to)[0]


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
        self.n_chains = n_chains
        self.rc = cutoff
        self.rc2 = cutoff ** 2
        self.alpha = alpha
        self.eps_r = eps_r
        self.torsion = tuple(torsion)
        self.xp = xp or bk.get_backend()
        self._dt = bk.float_dtype()
        self.n_energy_calls = 0  # batched kernel calls (one per ``energy()`` invocation)
        self.n_energy_rows = 0  # configurations evaluated in total
        self._setup_topology(chain)
        self.update_chain(chain)

    # --- setup -----------------------------------------------------------------
    def _setup_topology(self, chain: PeriodicChain) -> None:
        """Everything that depends only on the elements, charges and bonded topology.

        Survives :meth:`update_chain`: during a torsion refinement the chain's geometry
        changes every iteration but its topology never does.
        """
        xp, dt = self.xp, self._dt
        n, n_chains, rc, alpha = chain.n_atoms, self.n_chains, self.rc, self.alpha
        self.elements = list(chain.elements)
        x, d = lj_params(chain.elements)
        xs = np.sqrt(np.outer(x, x))
        ds = np.sqrt(np.outer(d, d))
        qq = np.outer(chain.charges, chain.charges) * COULOMB / self.eps_r
        erc = float(erfc_approx(np.array(alpha * rc)))
        self.dsf_shift = erc / rc
        self.dsf_force = erc / rc ** 2 + 2 * alpha / np.sqrt(np.pi) * np.exp(-(alpha * rc) ** 2) / rc
        s6 = (xs / rc) ** 6
        lj_shift = ds * (s6 * s6 - 2 * s6)
        # per-pair constants: v = ir6 (A ir6 - B) + qq erfc(alpha r)/r + qq_f r + const
        A = ds * xs ** 12
        B = 2.0 * ds * xs ** 6
        qq_f = qq * self.dsf_force
        const = -lj_shift - qq * (self.dsf_shift + self.dsf_force * rc)
        tile = lambda Z: np.tile(Z, (n_chains, n_chains))  # noqa: E731
        self._A, self._B = xp.asarray(tile(A), dtype=dt), xp.asarray(tile(B), dtype=dt)
        self._qq, self._qqf = xp.asarray(tile(qq), dtype=dt), xp.asarray(tile(qq_f), dtype=dt)
        self._const = xp.asarray(tile(const), dtype=dt)
        # the same n x n block also serves the (0,0,k) chain1-chain2 column and the intra sum
        self._A_nn, self._B_nn = xp.asarray(A, dtype=dt), xp.asarray(B, dtype=dt)
        self._qq_nn, self._qqf_nn = xp.asarray(qq, dtype=dt), xp.asarray(qq_f, dtype=dt)
        self._const_nn = xp.asarray(const, dtype=dt)
        # kept for backwards compatibility / introspection
        self.lj_x = xp.asarray(tile(xs), dtype=dt)
        self.lj_d = xp.asarray(tile(ds), dtype=dt)
        self.qq = self._qq
        self.lj_shift = xp.asarray(tile(lj_shift), dtype=dt)
        N = n * n_chains
        self.same_site = {}
        for k, S in chain.same_site_scale.items():
            Mk = np.ones((N, N))
            for s in range(n_chains):
                Mk[s * n : (s + 1) * n, s * n : (s + 1) * n] = S
            self.same_site[k] = Mk
        self._scale_nn_base = dict(chain.same_site_scale)
        self.n, self.N = n, N
        self._scale_cache: dict[int, object] = {}
        self._image_cache: dict = {}

    def update_chain(self, chain: PeriodicChain) -> None:
        """Replace the chain-dependent state, keeping the topology-dependent tables.

        Used by the refinement, where only the torsions (hence the coordinates, the
        repeat ``c``, the torsion energy, the intra-chain constant and the image
        selection) change from one evaluation to the next.
        """
        if list(chain.elements) != getattr(self, "elements", list(chain.elements)):
            raise ValueError("update_chain requires a chain with the same element list")
        self.chain = chain
        self.X0 = self.xp.asarray(chain.coords, dtype=self._dt)
        self.K = int(np.ceil(self.rc / chain.c)) + 1
        self.reach = self.rc + 2 * chain.radius + 0.5  # axis-axis distance beyond which no atom pair is within the cutoff
        self.e_torsion = self.torsion_energy(chain.dihedrals) * self.n_chains
        self._image_cache.clear()
        # The intra-chain (same-site, (0,0,k)) sum is a constant for a rigid chain, and a
        # flip is an isometry of the chain, so one value serves both orientations.
        self.e_intra = 0.5 * self.n_chains * float(bk.to_numpy(self._intra_column(self.X0[None], np.array([chain.c]), self.K))[0])

    def torsion_energy(self, dihedrals) -> float:
        """Fourier torsion energy of one chain's repeat (kcal/mol)."""
        V1, V2, V3 = self.torsion
        phi = np.deg2rad(np.asarray(dihedrals, dtype=float))
        return float((0.5 * (V1 * (1 + np.cos(phi)) + V2 * (1 - np.cos(2 * phi)) + V3 * (1 + np.cos(3 * phi)))).sum())

    def chain_batch(self, chains) -> dict:
        """``coords``/``c``/``e_torsion`` keyword arguments of :meth:`energy` for a list of chains."""
        chains = list(chains)
        return {
            "coords": np.stack([ch.coords for ch in chains]),
            "c": np.array([ch.c for ch in chains], dtype=float),
            "e_torsion": np.array([self.torsion_energy(ch.dihedrals) * self.n_chains for ch in chains]),
        }

    def _scale_column(self, K: int):
        """(2K+1, n, n) bonded-exclusion scales for the same-site (0,0,k) images."""
        S = self._scale_cache.get(K)
        if S is None:
            ones = np.ones((self.n, self.n))
            S = self.xp.asarray(np.stack([self._scale_nn_base.get(k, ones) for k in range(-K, K + 1)]), dtype=self._dt)
            self._scale_cache[K] = S
        return S

    # --- geometry --------------------------------------------------------------
    def _place(self, params, coords=None, c=None):
        """params (M, 7) -> primary coordinates (M, N, 3), lattice vectors (M, 3, 3).

        ``coords`` (M, n, 3) and ``c`` (M,) override the packer's own chain, so that one
        batched call can carry a different chain geometry per row.
        """
        xp = self.xp
        p = xp.asarray(params, dtype=self._dt)
        M = p.shape[0]
        a, b, gam, phi1, phi2, dz, flip = (p[:, i] for i in range(7))
        g = xp.deg2rad(gam)
        z = xp.zeros(M, dtype=self._dt)
        cz = xp.full(M, self.chain.c, dtype=self._dt) if c is None else xp.asarray(c, dtype=self._dt)
        avec = xp.stack([a, z, z], axis=1)
        bvec = xp.stack([b * xp.cos(g), b * xp.sin(g), z], axis=1)
        cvec = xp.stack([z, z, cz], axis=1)

        def rotz(X, ang):
            ca, sa = xp.cos(xp.deg2rad(ang))[:, None], xp.sin(xp.deg2rad(ang))[:, None]
            return xp.stack([ca * X[..., 0] - sa * X[..., 1], sa * X[..., 0] + ca * X[..., 1], X[..., 2]], axis=-1)

        X = xp.broadcast_to(self.X0, (M,) + self.X0.shape) if coords is None else xp.asarray(coords, dtype=self._dt)
        chains = [rotz(X, phi1)]
        if self.n_chains == 2:
            Xf = xp.where(flip[:, None, None] > 0.5, xp.stack([X[..., 0], -X[..., 1], -X[..., 2]], axis=-1), X)
            X2 = rotz(Xf, phi2) + (0.5 * avec + 0.5 * bvec)[:, None, :] + xp.stack([z, z, dz], axis=1)[:, None, :]
            chains.append(X2)
        return xp.concatenate(chains, axis=1), xp.stack([avec, bvec, cvec], axis=1)

    def _select_images(self, params: np.ndarray, K: int, reach: float):
        """Per configuration, the lateral images (i, j) != (0, 0) whose chain axes can be within reach.

        Returns ijk (M, I, 3), padded with far images that fall outside the cutoff.  The
        (0, 0) column is handled separately: its same-site part is the constant
        :attr:`e_intra` and its chain1-chain2 part is evaluated on an n x n block.
        """
        M = params.shape[0]
        a, b, g = params[:, 0], params[:, 1], np.deg2rad(params[:, 2])
        a_min = float(min(a.min(), (b * np.sin(g)).min()))
        R = int(np.ceil(reach / max(a_min, 1e-3))) + 1
        G = self._image_cache.get(R)
        if G is None:
            ij = np.stack(np.meshgrid(np.arange(-R, R + 1.0), np.arange(-R, R + 1.0), indexing="ij"), axis=-1).reshape(-1, 2)
            G = self._image_cache[R] = ij[np.any(ij != 0, axis=1)]
        av = np.stack([a, np.zeros(M)], axis=1)
        bv = np.stack([b * np.cos(g), b * np.sin(g)], axis=1)
        pos = G[None, :, 0, None] * av[:, None, :] + G[None, :, 1, None] * bv[:, None, :]  # (M, Gn, 2)
        deltas = [np.zeros((M, 2))]
        if self.n_chains == 2:
            o = 0.5 * av + 0.5 * bv
            deltas += [o, -o]
        dmin = np.min(np.stack([np.linalg.norm(pos + d[:, None, :], axis=2) for d in deltas]), axis=0)  # (M, Gn)
        keep = dmin < reach
        counts = keep.sum(axis=1)
        L = int(counts.max())
        if L == 0:
            return np.zeros((M, 0, 3)), 0
        # vectorised gather of the kept sites, padded with a far image (stable argsort
        # puts the kept columns first, in grid order)
        rank = np.argsort(~keep, axis=1, kind="stable")[:, :L]
        lat = G[rank]  # (M, L, 2)
        lat[np.arange(L)[None, :] >= counts[:, None]] = 3 * R + 5
        ks = np.arange(-K, K + 1, dtype=float)
        ijk = np.concatenate([np.repeat(lat[:, :, None, :], len(ks), axis=2), np.broadcast_to(ks[None, None, :, None], (M, L, len(ks), 1))], axis=3)
        return ijk.reshape(M, L * len(ks), 3), L

    # --- energy ------------------------------------------------------------------
    def _pair_energy(self, r2, A, B, qq, qqf, const):
        """Pair potential (LJ + damped-shifted-force Coulomb) on squared distances."""
        xp = self.xp
        mask = (r2 < self.rc2) & (r2 > 1e-8)
        r2s = xp.where(mask, r2, 1.0)
        rr = xp.sqrt(r2s)
        ir6 = 1.0 / (r2s * r2s * r2s)
        v = ir6 * (A * ir6 - B) + qq * (_erfc_pos(self.alpha * rr, xp) / rr) + qqf * rr + const
        return xp.where(mask, v, 0.0)

    def _intra_column(self, coords, c, K: int):
        """Same-site (0,0,k) energy of ONE chain, per row of ``coords`` (M, n, 3).

        This is the only place the bonded-exclusion scales are needed.  Evaluating it
        once per chain rather than once per configuration is what removes them from the
        kernel -- and it is also a correctness fix: a flipped chain's coordinates are
        mirrored to (x, -y, -z), which maps its z-image ``k`` onto image ``-k``, so
        scoring it with ``S[k]`` (as a per-configuration (0,0,k) column must, having one
        scale tensor for both chains) fails to exclude its bonded pairs and adds
        thousands of kcal/mol to every antiparallel cell.  A flip is an isometry of the
        chain, so the constant computed here is the same for either orientation.
        """
        xp = self.xp
        ks = xp.asarray(np.arange(-K, K + 1, dtype=float), dtype=self._dt)
        cz = xp.asarray(np.asarray(c, dtype=float), dtype=self._dt)
        D = coords[:, :, None, :] - coords[:, None, :, :]  # (M, n, n, 3)
        dz = D[:, None, :, :, 2] - ks[None, :, None, None] * cz[:, None, None, None]
        r2 = D[:, None, :, :, 0] ** 2 + D[:, None, :, :, 1] ** 2 + dz * dz
        v = self._pair_energy(r2, self._A_nn, self._B_nn, self._qq_nn, self._qqf_nn, self._const_nn)
        return (v * self._scale_column(K)[None]).sum(axis=(1, 2, 3))

    def energy(self, params, chunk_elems: int | None = None, coords=None, c=None, e_torsion=None) -> np.ndarray:
        """Lattice energy per cell (kcal/mol) for each row of params (M, 7).

        Configurations are sorted by cell area (so chunks share similar image counts and
        little padding is wasted) and processed in chunks of ``chunk_elems`` array
        elements: small on CPU (cache-resident), large on GPU (fewer kernel launches).

        ``coords`` (M, n, 3), ``c`` (M,) and ``e_torsion`` (M,) optionally give a
        *different chain geometry per row*, so that a gradient batch over torsions is a
        single call; without them the packer's own chain is used for every row and its
        intra-chain energy is the constant computed in :meth:`update_chain`.
        """
        xp = self.xp
        params = np.atleast_2d(np.asarray(params, dtype=float))
        M = params.shape[0]
        self.n_energy_calls += 1
        self.n_energy_rows += M
        if chunk_elems is None:
            chunk_elems = 60_000_000 if bk.device_name() == "cuda" else CPU_CHUNK_ELEMS
        per_row = coords is not None
        if per_row:
            coords = np.asarray(coords, dtype=float)
            if coords.ndim == 2:
                coords = coords[None]
            if coords.shape[0] == 1 and M > 1:
                coords = np.broadcast_to(coords, (M, self.n, 3))
            c_arr = np.full(M, self.chain.c) if c is None else np.broadcast_to(np.asarray(c, dtype=float).ravel(), (M,))
            K = int(np.ceil(self.rc / float(c_arr.min()))) + 1
            reach = self.rc + 2 * float(np.linalg.norm(coords[:, :, :2], axis=2).max()) + 0.5
            e_add = np.full(M, self.e_torsion) if e_torsion is None else np.broadcast_to(np.asarray(e_torsion, dtype=float).ravel(), (M,)).copy()
            i_chunk = max(1, chunk_elems // ((2 * K + 1) * self.n * self.n))
            intra = np.concatenate([
                bk.to_numpy(self._intra_column(xp.asarray(coords[t : t + i_chunk], dtype=self._dt), c_arr[t : t + i_chunk], K))
                for t in range(0, M, i_chunk)
            ])
            e_add = e_add + 0.5 * self.n_chains * intra
        else:
            c_arr = None
            K, reach = self.K, self.reach
            e_add = np.full(M, self.e_torsion + self.e_intra)
        order = np.argsort(params[:, 0] * params[:, 1] * np.sin(np.deg2rad(params[:, 2])))
        params = params[order]
        out = np.empty(M)
        N, n, nk = self.N, self.n, 2 * K + 1
        m_chunk = max(1, chunk_elems // (30 * nk * N * N))
        for s in range(0, M, m_chunk):
            sl = order[s : s + m_chunk]
            sub = params[s : s + m_chunk]
            sub_coords = coords[sl] if per_row else None
            sub_c = c_arr[sl] if per_row else None
            P, lat = self._place(sub, sub_coords, sub_c)
            e = xp.zeros(sub.shape[0], dtype=self._dt)
            if self.n_chains == 2:
                # (0, 0, k) column: chain1-chain2 only; the (2,1) block equals the (1,2) block
                Dh = P[:, :n, None, :] - P[:, None, n:, :]  # (m, n, n, 3)
                ks = xp.asarray(np.arange(-K, K + 1, dtype=float), dtype=self._dt)
                cz = lat[:, 2, 2]
                dzh = Dh[:, None, :, :, 2] - ks[None, :, None, None] * cz[:, None, None, None]
                r2h = Dh[:, None, :, :, 0] ** 2 + Dh[:, None, :, :, 1] ** 2 + dzh * dzh
                vh = self._pair_energy(r2h, self._A_nn, self._B_nn, self._qq_nn, self._qqf_nn, self._const_nn)
                e = e + 2.0 * vh.sum(axis=(1, 2, 3))
            ijk_np, L = self._select_images(sub, K, reach)
            if L:
                ijk = xp.asarray(ijk_np, dtype=self._dt)
                D = P[:, :, None, :] - P[:, None, :, :]  # (m, N, N, 3)
                shift = xp.einsum("mic,mcd->mid", ijk, lat)  # (m, I, 3)
                dx = D[:, None, :, :, 0] - shift[:, :, None, None, 0]
                dy = D[:, None, :, :, 1] - shift[:, :, None, None, 1]
                dzz = D[:, None, :, :, 2] - shift[:, :, None, None, 2]
                r2 = dx * dx + dy * dy + dzz * dzz  # (m, I, N, N)
                v = self._pair_energy(r2, self._A, self._B, self._qq, self._qqf, self._const)
                e = e + v.sum(axis=(1, 2, 3))
            out[s : s + m_chunk] = bk.to_numpy(0.5 * e)
        res = np.empty_like(out)
        res[order] = out
        return res + e_add

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


# finite-difference steps for the cell variables (a, b, gamma, phi1, phi2, dz)
FD_STEPS = np.array([1e-3, 1e-3, 0.05, 0.05, 0.05, 1e-3])


def cell_value_and_grad(packer: CrystalPacker, x: np.ndarray, free, lo=None, hi=None, steps=None):
    """Energy and its gradient w.r.t. the free cell parameters from ONE batched call.

    Builds the (1 + 2 n_free, 7) matrix of the centre and of central-difference
    displacements (clipped to the bounds, which keeps the difference one-sided rather
    than stepping outside a bound) and evaluates it in a single kernel call.
    """
    x = np.asarray(x, dtype=float)
    free = list(free)
    steps = FD_STEPS if steps is None else steps
    P = np.repeat(x[None], 1 + 2 * len(free), axis=0)
    for i, idx in enumerate(free):
        h = steps[idx]
        xp_, xm = x[idx] + h, x[idx] - h
        if lo is not None:
            xp_, xm = min(xp_, hi[idx]), max(xm, lo[idx])
        P[1 + 2 * i, idx], P[2 + 2 * i, idx] = xp_, xm
    E = packer.energy(P)
    grad = np.zeros(len(free))
    for i, idx in enumerate(free):
        d = P[1 + 2 * i, idx] - P[2 + 2 * i, idx]
        grad[i] = (E[1 + 2 * i] - E[2 + 2 * i]) / d if d else 0.0
    return float(E[0]), grad


def polish(
    packer: CrystalPacker,
    params: np.ndarray,
    lo: np.ndarray,
    hi: np.ndarray,
    free: list[int],
    maxfev: int = 1500,
    method: str = "lbfgs",
    maxiter: int = 200,
) -> np.ndarray:
    """Local minimisation of the continuous cell parameters (flip fixed).

    ``method="lbfgs"`` (default) runs L-BFGS-B on batched central-difference gradients:
    one kernel call of ``1 + 2 n_free`` configurations per function evaluation instead of
    the ~10-50x more sequential single-cell evaluations Nelder-Mead needs.
    ``method="nelder-mead"`` is the original simplex search with bounds by penalty.
    """
    x0 = np.asarray(params, dtype=float).copy()
    if method not in ("lbfgs", "l-bfgs-b", "nelder-mead", "nelder_mead"):
        raise ValueError(f"unknown polish method {method!r}")
    if method.startswith("nelder"):
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

    # phi1, phi2 and dz are periodic: leave them unbounded and wrap afterwards, so the
    # optimiser is never stopped by an artificial boundary at 0 or 360 deg.
    period = {3: 360.0, 4: 360.0, 5: float(packer.chain.c)}
    wrap = [i for i in free if i in period and hi[i] - lo[i] >= period[i] - 1e-9]
    bounds = [(None, None) if i in wrap else (lo[i], hi[i]) for i in free]
    lo_fd, hi_fd = np.asarray(lo, dtype=float).copy(), np.asarray(hi, dtype=float).copy()
    lo_fd[wrap], hi_fd[wrap] = -np.inf, np.inf  # no clipping of the periodic variables

    def fun(xf):
        x = x0.copy()
        x[free] = xf
        return cell_value_and_grad(packer, x, free, lo_fd, hi_fd)

    res = minimize(
        fun, x0[free], method="L-BFGS-B", jac=True, bounds=bounds,
        options={"ftol": 1e-9, "gtol": 1e-5, "maxiter": maxiter, "maxfun": max(maxfev, 2 * maxiter)},
    )
    x = x0.copy()
    x[free] = res.x
    for i in wrap:
        x[i] = lo[i] + (x[i] - lo[i]) % period[i]
    bounded = [i for i in range(6) if i not in wrap]
    x[bounded] = np.clip(x[bounded], lo[bounded], hi[bounded])
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
    method: str = "lbfgs",
) -> list[PackResult]:
    """Coarse batched random search followed by local polishing of the best cells.

    ``method`` selects the polisher: ``"lbfgs"`` (default, batched-gradient L-BFGS-B) or
    ``"nelder-mead"``.
    """
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
    results = [packer.result(polish(packer, p, lo, hi, free, maxfev, method=method)) for p in starts]
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
