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

:func:`pack` offers two coarse screens.  ``screen="random"`` samples cells
uniformly and scores them with this kernel.  ``screen="table"`` uses
:mod:`polyfind.lattice_table`: a chain-pair interaction tabulated once per
conformation makes the whole ``(phi1, phi2, dz)`` landscape of an ``(a, b)``
cell one inverse FFT, so the screen is exhaustive rather than sampled.  Both
hand their best cells to the same exact-kernel polish, so the energies that
come back are exact either way and only the choice of starts differs.

**Rigid table, deformable kernel.**  That tabulated interaction is built for one
fixed set of chain coordinates, one set of point charges and one energy-shifted pair
potential: it is a *screen*, and it stays rigid.  Everything that lets the chain deform
-- :class:`ChainValence`'s bond and angle terms, the force-shifted Lennard-Jones
form, and the geometry-dependent charges of :func:`chain_flux` -- is opt-in on
:class:`CrystalPacker` and lives on the direct-kernel path
(:meth:`CrystalPacker.energy` and :meth:`CrystalPacker.energy_and_grad`), which is
what :mod:`polyfind.refine` and :mod:`polyfind.mechanics` call.  The two are kept
apart by :func:`_table_starts`, which refuses a packer carrying either and refuses
a table whose chain is not the chain being packed, so a deformed chain can never
reach a rigid tabulation by accident.

A uniform applied electric field is optional (``field=(Ex, Ey, Ez)`` in V/A).  It
adds ``-mu_cell . E`` to the cell energy, where ``mu_cell = sum_i q_i r_i`` is the
dipole moment of the placed cell; the dipole depends on the setting angles and the
flip, so it is recomputed per configuration inside the kernel.  This is the cheap
screening heuristic of ``docs/CHEMISTRY_EXTENSION.md`` section 3, not a
Berry-phase or DFT treatment: the charges never depend on the field, so the cell has
no polarizability and no depolarisation field, and only the *orientational* and
packing response to the field is captured.  With ``charge_flux`` they do depend on
the chain's own geometry, which is a different thing and is what gives a planar
zigzag a piezoelectric response at all (:func:`chain_flux`).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import minimize

from . import backend as bk
from .chain import build_chain
from .forcefield import COULOMB, _valence_lookup, erfc_approx, valence_topology
from .helix import HelixParams, helix_parameters, kabsch, rotation_to_z
from .polymers import Polymer, RISStates, lj_params

MASS = {"C": 12.011, "H": 1.008, "N": 14.007, "O": 15.999, "F": 18.998, "Cl": 35.45}
_TOPO_CACHE: dict = {}

# 1 e/A^2 in C/m^2:  1.602176634e-19 C / (1e-10 m)^2.
E_PER_A2_TO_C_PER_M2 = 16.0218
# 1 eV in kcal/mol; converts (mu in e.A) . (E in V/A) to the kernel's energy unit.
EV_TO_KCAL = 23.0605

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


def _derfc_pos(x, xp=np):
    """d/dx of :func:`_erfc_pos`, for x >= 0.

    The *approximation* is differentiated, not ``erfc`` itself: the kernel's Coulomb
    term is ``_erfc_pos``, so this is what makes the analytic gradient the exact
    gradient of the function the kernel evaluates (to within rounding) rather than of
    the function it approximates.  It matters more than the 1.5e-7 accuracy of the
    polynomial's *value* suggests: its derivative differs from the true
    ``-2/sqrt(pi) exp(-x^2)`` by 1.2e-5 relative over the range used here, and
    substituting the true one costs an order of magnitude on the measured gradient
    error (see docs/PERFORMANCE_REVIEW.md section 13).
    """
    c1, c2, c3, c4, c5 = 0.254829592, -0.284496736, 1.421413741, -1.453152027, 1.061405429
    t = 1.0 / (1.0 + 0.3275911 * x)
    poly = t * (c1 + t * (c2 + t * (c3 + t * (c4 + t * c5))))
    dpoly = c1 + t * (2.0 * c2 + t * (3.0 * c3 + t * (4.0 * c4 + t * 5.0 * c5)))
    return (dpoly * (-0.3275911 * t * t) - 2.0 * x * poly) * xp.exp(-x * x)


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
    from .chain import build_backbone, nerf, pendant_positions

    M, nb = tors_batch.shape
    dih = np.tile(tors_batch, (1, 5))
    N = dih.shape[1]
    B, L = polymer.bonds_per_repeat, polymer.bond_length
    bb = build_backbone(polymer, dih, xp=np)  # (M, N+3, 3)
    v0 = nerf(bb[:, 2], bb[:, 1], bb[:, 0], L, polymer.backbone[0].backbone_angle, 180.0, xp=np)
    v1 = nerf(bb[:, N], bb[:, N + 1], bb[:, N + 2], L, polymer.backbone[(N + 2) % B].backbone_angle, 180.0, xp=np)
    bb_ext = np.concatenate([v0[:, None], bb, v1[:, None]], axis=1)  # (M, N+5, 3)
    # A pendant is a group of one or more atoms (phase 3), so the stride is per backbone
    # atom rather than a fixed 3; the layout matches build_chain's exactly (backbone atom,
    # then each pendant's atoms in order).
    widths = [1 + polymer.backbone[k % B].n_pendant_atoms for k in range(N + 3)]
    starts = np.concatenate([[0], np.cumsum(widths)])
    out = np.empty((M, int(starts[-1]), 3))
    for k in range(N + 3):
        spec = polymer.backbone[k % B]
        x = bb_ext[:, k + 1]
        at = int(starts[k])
        out[:, at] = x
        groups = pendant_positions(bb_ext[:, k], x, bb_ext[:, k + 2], spec, xp=np)
        at += 1
        for positions in groups:
            for pos in positions:
                out[:, at] = pos
                at += 1
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


# ------------------------------------------------------------------ valence terms
@dataclass
class ChainValence:
    """Harmonic bond and angle terms of **one repeat** of an infinite periodic chain.

    The chain the packer carries is one crystallographic repeat, and the bonds and angles
    of the infinite chain do not stop at its edges: a backbone bond joins the last atom of
    the repeat to the first atom of the *next* one, which is the same atom translated by
    ``c`` along z.  Every term here therefore carries the image index of each of its atoms
    (``0`` for the repeat itself, ``+/-1`` for a neighbouring one), and the energy depends
    on ``c`` as well as on the coordinates -- which is exactly the dependence that gives an
    axial strain something to push against.

    One representative per term per repeat: a bond is kept when the lower of its two image
    indices is 0, an angle when its *central* atom is in the repeat.  Both rules pick each
    term of the infinite chain exactly once, so :meth:`energy` is the valence energy *per
    repeat* and adding it to a lattice energy per cell needs one factor of ``n_chains`` and
    nothing else.

    **Bond stretching is inert here, and that is a property of the geometry, not of this
    class.**  :func:`polyfind.chain.build_chain` places every atom at the polymer's own bond
    length, so ``r - r0`` is fixed by the chemistry and the stretch block contributes a
    constant and a zero gradient.  It is evaluated anyway, because leaving it out would make
    the reported energy not the potential's, and because a future flexible builder would
    need nothing changed here.
    """

    bond_i: np.ndarray  # (nb,) atom in the repeat
    bond_j: np.ndarray  # (nb,) the other atom, in image ``bond_s``
    bond_s: np.ndarray  # (nb,) image of ``bond_j``, in units of c along z
    bond_k: np.ndarray  # (nb,) kcal/(mol A^2)
    bond_r0: np.ndarray  # (nb,) A
    ang_i: np.ndarray  # (na,)
    ang_j: np.ndarray  # (na,) the central atom, always in the repeat
    ang_k: np.ndarray  # (na,)
    ang_si: np.ndarray  # (na,) image of ``ang_i``
    ang_sk: np.ndarray  # (na,) image of ``ang_k``
    ang_kk: np.ndarray  # (na,) kcal/(mol rad^2)
    ang_t0: np.ndarray  # (na,) radians

    @property
    def n_bonds(self) -> int:
        return int(self.bond_i.size)

    @property
    def n_angles(self) -> int:
        return int(self.ang_i.size)

    def _bond_vectors(self, X: np.ndarray, cz: np.ndarray) -> np.ndarray:
        d = X[:, self.bond_i] - X[:, self.bond_j]
        d[..., 2] -= self.bond_s[None, :] * cz[:, None]
        return d

    def _angle_vectors(self, X: np.ndarray, cz: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        u = X[:, self.ang_i] - X[:, self.ang_j]
        v = X[:, self.ang_k] - X[:, self.ang_j]
        u[..., 2] += self.ang_si[None, :] * cz[:, None]
        v[..., 2] += self.ang_sk[None, :] * cz[:, None]
        return u, v

    def energy(self, coords, c) -> np.ndarray:
        """Valence energy per repeat (kcal/mol), one entry per row of ``coords`` (M, n, 3)."""
        X = np.array(coords, dtype=float)
        if X.ndim == 2:
            X = X[None]
        M = X.shape[0]
        cz = np.broadcast_to(np.asarray(c, dtype=float).ravel(), (M,))
        out = np.zeros(M)
        if self.n_bonds:
            r = np.linalg.norm(self._bond_vectors(X, cz), axis=-1)
            out += 0.5 * (self.bond_k * (r - self.bond_r0) ** 2).sum(axis=1)
        if self.n_angles:
            u, v = self._angle_vectors(X, cz)
            nu, nv = np.linalg.norm(u, axis=-1), np.linalg.norm(v, axis=-1)
            cos = np.clip((u * v).sum(-1) / (nu * nv), -1.0, 1.0)
            da = np.arccos(cos) - self.ang_t0
            out += 0.5 * (self.ang_kk * da * da).sum(axis=1)
        return out

    def energy_and_grad(self, coords, c, n_atoms: int | None = None):
        """``(E, dE/dcoords (n, 3), dE/dc)`` for ONE repeat geometry.

        Closed form throughout: a stretch contributes ``k (r - r0) rhat`` along its own
        separation, a bend the usual ``dtheta/dx`` of the two arms and minus their sum at
        the centre.  ``c`` enters through the image offsets alone, so its derivative is the
        z components of those forces weighted by the image indices.
        """
        X = np.asarray(coords, dtype=float).reshape(1, -1, 3)
        n = X.shape[1] if n_atoms is None else int(n_atoms)
        cz = np.asarray(c, dtype=float).reshape(1)
        E = 0.0
        gX = np.zeros((n, 3))
        gc = 0.0
        if self.n_bonds:
            d = self._bond_vectors(X, cz)[0]
            r = np.linalg.norm(d, axis=-1)
            dr = r - self.bond_r0
            E += 0.5 * float((self.bond_k * dr * dr).sum())
            f = (self.bond_k * dr / r)[:, None] * d  # dE/d(separation vector)
            np.add.at(gX, self.bond_i, f)
            np.add.at(gX, self.bond_j, -f)
            gc -= float((self.bond_s * f[:, 2]).sum())
        if self.n_angles:
            u, v = self._angle_vectors(X, cz)
            u, v = u[0], v[0]
            nu, nv = np.linalg.norm(u, axis=-1), np.linalg.norm(v, axis=-1)
            uh, vh = u / nu[:, None], v / nv[:, None]
            cos = np.clip((uh * vh).sum(-1), -1.0, 1.0)
            theta = np.arccos(cos)
            da = theta - self.ang_t0
            E += 0.5 * float((self.ang_kk * da * da).sum())
            s = np.maximum(np.sqrt(1.0 - cos * cos), 1e-9)
            gi = (cos[:, None] * uh - vh) / (nu * s)[:, None]
            gk = (cos[:, None] * vh - uh) / (nv * s)[:, None]
            w = (self.ang_kk * da)[:, None]
            fi, fk = w * gi, w * gk
            np.add.at(gX, self.ang_i, fi)
            np.add.at(gX, self.ang_k, fk)
            np.add.at(gX, self.ang_j, -(fi + fk))
            gc += float((self.ang_si * fi[:, 2]).sum() + (self.ang_sk * fk[:, 2]).sum())
        return E, gX, gc


def chain_valence(chain: PeriodicChain, bond_table: dict, angle_table: dict) -> ChainValence | None:
    """The :class:`ChainValence` of ``chain`` for ``{type: (k, r0)}`` / ``{type: (k, theta0_deg)}``.

    ``None`` when both tables are empty, which is what keeps a packer without valence terms
    computing exactly the expression it always did.  Only the *topology* of the chain is
    used, so the result survives any change of conformation: the same object serves every
    geometry of the same repeat, which is what a refinement needs.
    """
    if not bond_table and not angle_table:
        return None
    nb_block = len(chain.dihedrals)
    tpl = _template(chain.polymer, nb_block)
    where: dict[int, tuple[int, int]] = {}
    for k in (-2, -1, 0, 1, 2):
        for p, idx in enumerate(_block_atoms(tpl, nb_block, k)[1]):
            where[int(idx)] = (p, k)
    bl, al = valence_topology(tpl.elements, tpl.bonds)
    bi, bj, bs, bk, br = [], [], [], [], []
    for i, j, name in bl:
        wi, wj = where.get(i), where.get(j)
        if wi is None or wj is None or min(wi[1], wj[1]) != 0:
            continue  # the representative of this bond sits in another repeat
        if not bond_table:
            continue
        (p, _), (q, s) = (wi, wj) if wi[1] == 0 else (wj, wi)
        k, r0 = _valence_lookup(bond_table, name, "bond")
        bi.append(p), bj.append(q), bs.append(s), bk.append(k), br.append(r0)
    ai, aj, ak, asi, ask, akk, at = [], [], [], [], [], [], []
    for i, j, k, name in al:
        wj = where.get(j)
        if wj is None or wj[1] != 0 or not angle_table:
            continue
        wi, wk = where.get(i), where.get(k)
        if wi is None or wk is None:
            continue  # pragma: no cover - a neighbour of a block-0 atom is always in -1..1
        kk, t0 = _valence_lookup(angle_table, name, "angle")
        ai.append(wi[0]), aj.append(wj[0]), ak.append(wk[0])
        asi.append(wi[1]), ask.append(wk[1]), akk.append(kk), at.append(t0)
    ints = lambda v: np.array(v, dtype=int)  # noqa: E731
    reals = lambda v: np.array(v, dtype=float)  # noqa: E731
    return ChainValence(
        bond_i=ints(bi), bond_j=ints(bj), bond_s=reals(bs), bond_k=reals(bk), bond_r0=reals(br),
        ang_i=ints(ai), ang_j=ints(aj), ang_k=ints(ak), ang_si=reals(asi), ang_sk=reals(ask),
        ang_kk=reals(akk), ang_t0=np.deg2rad(reals(at)),
    )


# ------------------------------------------------------------------ charge flux
def _repeat_images(chain: PeriodicChain) -> tuple:
    """``(template, {template atom: (position in the repeat, image)})`` for ``chain``.

    The same five-block template :func:`chain_valence` walks, mapped so that a term of
    the infinite chain can be written in the repeat's own atom order with an image index
    for each atom.  Blocks -2..2 are enough: every neighbour of an atom of block 0 is in
    block -1, 0 or 1, and every neighbour of *those* in -2..2.
    """
    tpl = _template(chain.polymer, len(chain.dihedrals))
    where: dict[int, tuple[int, int]] = {}
    for k in (-2, -1, 0, 1, 2):
        for p, idx in enumerate(_block_atoms(tpl, len(chain.dihedrals), k)[1]):
            where[int(idx)] = (p, k)
    return tpl, where


def chain_flux(chain: PeriodicChain, ff):
    """The :class:`~polyfind.forcefield.FluxTopology` of one repeat of ``chain``.

    ``ff`` is anything with :meth:`~polyfind.forcefield.SimpleFF.flux_topology` -- a
    :class:`~polyfind.forcefield.SimpleFF` carrying ``charge_increments`` and
    ``charge_flux``.  ``None`` when the potential has no non-zero flux coefficient, which
    is what keeps a packer without flux computing exactly the expression it always did.

    Only the *topology* of the chain is used, so the result survives any change of
    conformation, exactly as :func:`chain_valence` does.  The flux's drivers reach across
    the repeat boundary -- a pendant bond's angle driver involves the backbone neighbours
    of its own carbon, one of which belongs to the next repeat -- so the charges depend on
    ``c`` as well as on the coordinates, and both derivatives are carried.
    """
    if ff is None or not ff.has_flux():
        return None
    tpl, where = _repeat_images(chain)
    return ff.flux_topology(tpl.elements, tpl.bonds, images=where)


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
    dipole: np.ndarray | None = None  # cell dipole mu = sum_i q_i r_i (e.A); None if the cell is not neutral
    polarization: np.ndarray | None = None  # P = mu / V (C/m^2)
    field: tuple[float, float, float] = (0.0, 0.0, 0.0)  # applied uniform field (V/A)
    field_energy: float = 0.0  # the -mu.E term already included in energy_per_cell (kcal/mol)

    @property
    def params(self) -> np.ndarray:
        return np.array([self.a, self.b, self.gamma, self.phi1, self.phi2, self.dz, self.flip], dtype=float)

    @property
    def polarization_magnitude(self) -> float:
        """|P| in C/m^2 (0 when the cell charge made the dipole undefined)."""
        return 0.0 if self.polarization is None else float(np.linalg.norm(self.polarization))

    def row(self) -> str:
        pol = "" if self.polarization is None else f"  |P|={self.polarization_magnitude:5.3f}"
        return (
            f"{self.chain:<14s} E/mon={self.energy_per_monomer:8.3f}  a={self.a:5.2f} b={self.b:5.2f} c={self.c:5.2f} "
            f"gamma={self.gamma:5.1f}  phi=({self.phi1:5.1f},{self.phi2:5.1f}) dz={self.dz:4.2f} {'anti' if self.flip else 'para'}  rho={self.density:5.3f}"
            + pol
        )


class CrystalPacker:
    """Batched lattice-energy evaluator for a rigid periodic chain.

    ``field`` is an optional uniform applied electric field ``(Ex, Ey, Ez)`` in V/A,
    which adds ``-mu_cell . E`` to every configuration's energy; see :meth:`dipole`
    and :meth:`field_energy`.  ``field=None`` (the default) leaves the kernel exactly
    as it was -- the field term is not evaluated at all -- so zero-field energies are
    bit-for-bit those of the fieldless code.

    **Two opt-in changes to the potential, both off by default.**  A default-constructed
    packer is bit-for-bit the one every number in this package was measured with
    (``tests/test_pack.py`` and ``tests/test_mechanics.py`` assert it with ``==``), and
    each of these has to be asked for by name:

    ``valence``
        anything with ``bond_table()`` and ``angle_table()`` methods -- a
        :class:`polyfind.forcefield.SimpleFF`, typically
        ``SimpleFF.from_preset("pvdf-dft-valence")``.  Its harmonic stretch and bend terms
        are then evaluated on the chain *as a periodic object* (:class:`ChainValence`) and
        added to every configuration, gradient included.  Within a rigid packing that is a
        constant per chain, which is why :meth:`polyfind.fitting.FFParameters.applied`
        does not forward it; it stops being a constant the moment the chain is allowed to
        deform, which is what :func:`polyfind.mechanics.axial_response` needs and what
        makes ``dE/dc`` a restoring force rather than a report of an invented restraint.
    ``lj_cutoff``
        ``"energy"`` (default) shifts the Lennard-Jones term by ``V(rc)``, leaving its
        *force* discontinuous at the cutoff; ``"force"`` subtracts the tangent instead,
        ``V(r) - V(rc) - (r - rc) V'(rc)``, so value and force both reach zero at ``rc``.
        The damped-shifted-force Coulomb term is already force-shifted, so with
        ``"force"`` the whole pair potential is C1 and a strain no longer walks pairs
        across a jump in the force (see :meth:`_pair_energy_and_dv`).

    ``charge_flux``
        a :class:`~polyfind.forcefield.SimpleFF` carrying ``charge_increments`` **and**
        ``charge_flux``, typically ``SimpleFF.from_preset("pvdf-dft-valence-flux")``.  The
        chain's point charges then stop being constants: they are recomputed from the
        chain's own geometry at every :meth:`update_chain` and per row wherever a row
        carries its own ``coords``, and every derivative of the energy and of the dipole
        carries the extra ``dE/dq . dq/dgeometry`` term
        (:class:`~polyfind.forcefield.FluxTopology`).  Without it a planar all-trans
        zigzag's dipole is *exactly* independent of its backbone angle and beta-PVDF's
        ``d_33`` and ``d_31`` are identically zero; see ``docs/ELECTROMECHANICS.md``.  The
        charges it starts from are the potential's own bond-charge increments, so a packer
        given flux **replaces** whatever charges the chain carried, and a mismatch between
        the two at zero flux is refused rather than absorbed.

    All three change the energy, so all three are refused by the tabulated screen
    (:func:`_table_starts`), which is a rigid-chain, energy-shifted, fixed-charge object by
    construction.
    """

    PARAMS = ("a", "b", "gamma", "phi1", "phi2", "dz", "flip")
    NEUTRAL_TOL = 1e-6  # |cell charge| (e) above which the dipole is origin-dependent
    LJ_CUTOFFS = ("energy", "force")

    def __init__(self, chain: PeriodicChain, n_chains: int = 2, cutoff: float = 8.0, alpha: float = 0.2, eps_r: float = 1.0, torsion=(1.3, -0.05, 2.5), xp=None, field=None, valence=None, lj_cutoff: str = "energy", charge_flux=None):
        self.n_chains = n_chains
        self.rc = cutoff
        self.rc2 = cutoff ** 2
        self.alpha = alpha
        self.eps_r = eps_r
        self.torsion = tuple(torsion)
        if lj_cutoff not in self.LJ_CUTOFFS:
            raise ValueError(f"unknown lj_cutoff {lj_cutoff!r} (expected one of {self.LJ_CUTOFFS})")
        self.lj_cutoff = lj_cutoff
        self.valence = valence
        self.charge_flux = charge_flux
        self.xp = xp or bk.get_backend()
        self._dt = bk.float_dtype()
        self.n_energy_calls = 0  # batched kernel calls (one per ``energy()`` invocation)
        self.n_energy_rows = 0  # configurations evaluated in total
        self.n_grad_calls = 0  # :meth:`energy_and_grad` invocations (one configuration each)
        self._setup_topology(chain)
        self.update_chain(chain)
        self.set_field(field)

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
        erc = float(erfc_approx(np.array(alpha * rc)))
        self.dsf_shift = erc / rc
        self.dsf_force = erc / rc ** 2 + 2 * alpha / np.sqrt(np.pi) * np.exp(-(alpha * rc) ** 2) / rc
        s6 = (xs / rc) ** 6
        lj_shift = ds * (s6 * s6 - 2 * s6)
        # per-pair constants: v = ir6 (A ir6 - B) + qq erfc(alpha r)/r + qq_f r + const
        A = ds * xs ** 12
        B = 2.0 * ds * xs ** 6
        # Force-shifted Lennard-Jones: subtract the tangent at rc, not just the value.
        # ``qq_f`` is already the linear coefficient the damped-shifted-force Coulomb term
        # needs, so the LJ force shift rides in the same array at no extra cost -- after
        # this it is no longer "qq times something", which is why the kernel calls it a
        # linear coefficient and not a Coulomb term.
        self._lj_shift_np = lj_shift
        self._lj_f_np = (12.0 * A / rc ** 13 - 6.0 * B / rc ** 7) if self.lj_cutoff == "force" else None
        tile = lambda Z: np.tile(Z, (n_chains, n_chains))  # noqa: E731
        self._A, self._B = xp.asarray(tile(A), dtype=dt), xp.asarray(tile(B), dtype=dt)
        # the same n x n block also serves the (0,0,k) chain1-chain2 column and the intra sum
        self._A_nn, self._B_nn = xp.asarray(A, dtype=dt), xp.asarray(B, dtype=dt)
        self._set_charge_tables(chain.charges)
        # kept for backwards compatibility / introspection
        self.lj_x = xp.asarray(tile(xs), dtype=dt)
        self.lj_d = xp.asarray(tile(ds), dtype=dt)
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
        # Valence terms depend on the topology alone, so they are resolved once here and
        # re-evaluated (not rebuilt) by every :meth:`update_chain`.
        self._valence = None if self.valence is None else chain_valence(
            chain, self.valence.bond_table(), self.valence.angle_table())
        # Charge flux likewise: the topology is resolved once, the charges are recomputed.
        self._flux = chain_flux(chain, self.charge_flux)
        if self._flux is not None:
            self._check_flux_base(chain)

    def _set_charge_tables(self, q) -> None:
        """(Re)build every charge-dependent pair table from the repeat's charges ``q`` (n,).

        Called once from :meth:`_setup_topology` with the chain's own charges -- which is the
        only thing that happens without charge flux, so the tables are bit-for-bit what they
        always were -- and again from :meth:`update_chain` whenever the flux has moved them.
        """
        xp, dt, rc = self.xp, self._dt, self.rc
        qq = np.outer(q, q) * COULOMB / self.eps_r
        qq_f = qq * self.dsf_force
        const = -self._lj_shift_np - qq * (self.dsf_shift + self.dsf_force * rc)
        if self._lj_f_np is not None:
            qq_f = qq_f + self._lj_f_np
            const = const - self._lj_f_np * rc
        tile = lambda Z: np.tile(Z, (self.n_chains, self.n_chains))  # noqa: E731
        self._qq, self._qqf = xp.asarray(tile(qq), dtype=dt), xp.asarray(tile(qq_f), dtype=dt)
        self._const = xp.asarray(tile(const), dtype=dt)
        self._qq_nn, self._qqf_nn = xp.asarray(qq, dtype=dt), xp.asarray(qq_f, dtype=dt)
        self._const_nn = xp.asarray(const, dtype=dt)
        self.qq = self._qq  # kept for backwards compatibility / introspection
        # charges of the whole cell, in the atom order :meth:`_place` produces
        self._q_cell = np.tile(np.asarray(q, dtype=float), self.n_chains)
        self._q_cell_xp = xp.asarray(self._q_cell, dtype=dt)

    def _row_charge_tables(self, q, cell: bool):
        """``(qq, qq_f, const)`` for per-row charges ``q`` (M, n); shapes (M, 1, n, n).

        The kernel's pair arrays are indexed ``(row, image, atom, atom)``, so a leading row
        axis and a broadcast image axis is all a per-row charge needs: every expression in
        :meth:`_pair_energy` is already elementwise.  ``cell=True`` tiles the block up to
        the whole cell, as :meth:`_set_charge_tables` does for the fixed case.
        """
        xp, dt, rc = self.xp, self._dt, self.rc
        q = np.asarray(q, dtype=float)
        qq = q[:, :, None] * q[:, None, :] * COULOMB / self.eps_r
        if cell:
            qq = np.tile(qq, (1, self.n_chains, self.n_chains))
        qq_f = qq * self.dsf_force
        shift = self._lj_shift_np if not cell else np.tile(self._lj_shift_np, (self.n_chains, self.n_chains))
        const = -shift[None] - qq * (self.dsf_shift + self.dsf_force * rc)
        if self._lj_f_np is not None:
            f = self._lj_f_np if not cell else np.tile(self._lj_f_np, (self.n_chains, self.n_chains))
            qq_f = qq_f + f[None]
            const = const - f[None] * rc
        to = lambda Z: xp.asarray(Z[:, None], dtype=dt)  # noqa: E731 -- the broadcast image axis
        return to(qq), to(qq_f), to(const)

    def _check_flux_base(self, chain: PeriodicChain) -> None:
        """Refuse a flux whose zero-coefficient charges are not the chain's own.

        The flux is a perturbation of the potential's bond-charge increments, so at zero
        coefficient it must reproduce exactly the charges the chain already carries.  If it
        does not, the packer has been given a chain built for a different charge model and
        switching it on would quietly change the electrostatics as well as add the flux.
        """
        from dataclasses import replace as _replace

        base = _replace(self._flux, ka=np.zeros_like(self._flux.ka), kb=np.zeros_like(self._flux.kb),
                        ea=self._flux.ea[:0], ek=self._flux.ek[:0], eks=self._flux.eks[:0])
        q0 = base.charges(chain.coords, chain.c)[0]
        err = float(np.abs(q0 - np.asarray(chain.charges, dtype=float)).max())
        if err > 1e-9:
            raise ValueError(
                f"charge flux: at zero coefficients the potential's increments give charges that "
                f"differ from the chain's by {err:.3e} e.  The flux perturbs the increments, so "
                "the chain must already carry them -- build it inside FFParameters.applied() (or "
                "pass the same preset's charges), rather than letting the flux silently replace a "
                "different charge model.")

    def update_chain(self, chain: PeriodicChain) -> None:
        """Replace the chain-dependent state, keeping the topology-dependent tables.

        Used by the refinement, where only the torsions (hence the coordinates, the
        repeat ``c``, the torsion energy, the intra-chain constant and the image
        selection) change from one evaluation to the next.

        With charge flux the *charges* move too -- they are a function of the conformation
        -- so the charge-dependent pair tables are rebuilt here as well.
        """
        if list(chain.elements) != getattr(self, "elements", list(chain.elements)):
            raise ValueError("update_chain requires a chain with the same element list")
        self.chain = chain
        if getattr(self, "_flux", None) is not None:
            self._set_charge_tables(self._flux.charges(chain.coords, chain.c)[0])
        self.X0 = self.xp.asarray(chain.coords, dtype=self._dt)
        self.K = int(np.ceil(self.rc / chain.c)) + 1
        self.reach = self.rc + 2 * chain.radius + 0.5  # axis-axis distance beyond which no atom pair is within the cutoff
        self.e_torsion = self.torsion_energy(chain.dihedrals) * self.n_chains
        self._image_cache.clear()
        # The intra-chain (same-site, (0,0,k)) sum is a constant for a rigid chain, and a
        # flip is an isometry of the chain, so one value serves both orientations.
        self.e_intra = 0.5 * self.n_chains * float(bk.to_numpy(self._intra_column(self.X0[None], np.array([chain.c]), self.K))[0])
        # The valence energy of the *whole cell*: one repeat's worth per chain (a flip is an
        # isometry, so both orientations cost the same).  Unlike ``e_intra`` this is a
        # constant only while the chain is; it moves as soon as an angle does.
        self.e_valence = 0.0 if self._valence is None else self.n_chains * float(
            self._valence.energy(chain.coords[None], np.array([chain.c]))[0])

    def torsion_energy(self, dihedrals) -> float:
        """Fourier torsion energy of one chain's repeat (kcal/mol)."""
        V1, V2, V3 = self.torsion
        phi = np.deg2rad(np.asarray(dihedrals, dtype=float))
        return float((0.5 * (V1 * (1 + np.cos(phi)) + V2 * (1 - np.cos(2 * phi)) + V3 * (1 + np.cos(3 * phi)))).sum())

    # --- dipole, polarization and the applied field --------------------------------
    @property
    def cell_charge(self) -> float:
        """Total charge of the cell (e).  Zero for chains built from neutral repeat units."""
        return float(self._q_cell.sum())

    @property
    def is_neutral(self) -> bool:
        return abs(self.cell_charge) <= self.NEUTRAL_TOL

    def _require_neutral(self, what: str) -> None:
        if not self.is_neutral:
            raise ValueError(
                f"{what} needs a neutral cell: sum_i q_i = {self.cell_charge:+.3e} e for this chain, and for a "
                "charged cell sum_i q_i r_i depends on where the origin is put, so it is not a dipole moment at all"
            )

    def set_field(self, field) -> None:
        """Set (``(Ex, Ey, Ez)`` in V/A) or clear (``None``) the uniform applied field."""
        if field is None:
            self.field, self._field_on, self._field_xp = None, False, None
            return
        f = np.asarray(field, dtype=float).reshape(3)
        self.field = f
        self._field_on = bool(np.any(f != 0.0))
        self._field_xp = self.xp.asarray(f, dtype=self._dt) if self._field_on else None
        if self._field_on:
            self._require_neutral("an applied field")

    def dipole(self, params, coords=None, c=None) -> np.ndarray:
        """Cell dipole moment ``mu = sum_i q_i r_i`` in e.A, one (3,) row per configuration.

        ``params`` is (M, 7) as for :meth:`energy`; ``coords`` (M, n, 3) and ``c`` (M,)
        optionally give a different chain geometry per row.  The sum runs over the atoms
        of the cell *as placed*: each chain's crystallographic repeat is kept whole (the
        pendant atoms that overhang the repeat are not wrapped back into it), which is
        the molecular choice of branch for the polarization -- wrapping an atom by a
        lattice vector would shift ``mu`` by ``q`` times that vector, the usual
        polarization quantum.

        **This is only well defined because every repeat unit of these polymers is
        constructed neutral.**  ``sum_i q_i r_i`` changes by ``Q * d`` when the origin
        moves by ``d``, so for a cell of total charge ``Q != 0`` it is not a property of
        the cell at all; the method therefore refuses to compute it (see
        :attr:`cell_charge`).  Being neutral, the result is independent of the origin,
        of which periodic image of each chain is used, and of ``dz`` (translating a
        neutral chain does not change the cell dipole) -- but *not* of the setting
        angles or the flip, which rotate and mirror the charge distribution.

        With charge flux the charges are a function of the chain's geometry, so a row that
        carries its own ``coords`` carries its own charges too, and that is where
        beta-PVDF's piezoelectric response comes from: the fixed-increment dipole of a
        planar zigzag is exactly invariant under the one internal coordinate an axial
        strain moves.
        """
        self._require_neutral("the cell dipole")
        params = np.atleast_2d(np.asarray(params, dtype=float))
        P, _ = self._place(params, coords, c)
        Pn = np.asarray(bk.to_numpy(P), dtype=float)
        if self._flux is not None and coords is not None:
            q = self._row_charges(coords, c, Pn.shape[0])
            return np.einsum("mn,mnc->mc", np.tile(q, (1, self.n_chains)), Pn)
        return np.einsum("n,mnc->mc", self._q_cell, Pn)

    def _row_charges(self, coords, c, M: int) -> np.ndarray:
        """The repeat's fluxed charges for each of ``M`` rows, ``(M, n)``."""
        X = np.asarray(coords, dtype=float)
        if X.ndim == 2:
            X = X[None]
        if X.shape[0] == 1 and M > 1:
            X = np.broadcast_to(X, (M, self.n, 3))
        cz = np.full(M, self.chain.c) if c is None else np.broadcast_to(
            np.asarray(c, dtype=float).ravel(), (M,))
        return self._flux.charges(X, cz)

    def cell_volume(self, params, c=None) -> np.ndarray:
        """Cell volume ``a b sin(gamma) c`` in A^3, one entry per row of ``params``."""
        params = np.atleast_2d(np.asarray(params, dtype=float))
        cz = np.full(params.shape[0], self.chain.c) if c is None else np.broadcast_to(
            np.asarray(c, dtype=float).ravel(), (params.shape[0],))
        return params[:, 0] * params[:, 1] * np.sin(np.deg2rad(params[:, 2])) * cz

    def polarization(self, params, coords=None, c=None) -> np.ndarray:
        """Cell polarization ``P = mu / V`` in C/m^2, one (3,) row per configuration.

        The conversion is ``1 e/A^2 = 16.0218 C/m^2``.  The caveats of :meth:`dipole`
        apply, and one more: these are fixed point charges, so this is the *rigid-ion*
        polarization of the illustrative charge model, not a Berry-phase polarization,
        and it carries no electronic contribution.
        """
        return self.dipole(params, coords, c) / self.cell_volume(params, c)[:, None] * E_PER_A2_TO_C_PER_M2

    def field_energy(self, params, coords=None, c=None) -> np.ndarray:
        """The ``-mu . E`` coupling of each configuration, in kcal/mol per cell.

        ``mu`` in e.A dotted with a field in V/A is an energy in eV, hence the factor
        :data:`EV_TO_KCAL`.  Zeros (without touching the coordinates) when no field is
        set: this is the term :meth:`energy` adds, so it is exactly what the field
        contributes to a reported lattice energy.
        """
        M = np.atleast_2d(np.asarray(params, dtype=float)).shape[0]
        if not self._field_on:
            return np.zeros(M)
        return -(self.dipole(params, coords, c) @ self.field) * EV_TO_KCAL

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

    def _pair_energy_and_dv(self, r2, A, B, qq, qqf, const):
        """:meth:`_pair_energy` together with ``dV/d(r^2)``, on the same mask.

        The value expression is character for character the one in :meth:`_pair_energy`,
        so the energies built from this routine are bit-identical to the kernel's.  With
        ``V = A u^-6 - B u^-3 + qq erfc(a r)/r + qq_f r + const`` and ``u = r^2``,

            dV/du = 3 u^-3 (B - 2 A u^-3) / u + (qq f'(r) + qq_f) / (2 r),
            f(r)  = erfc(a r) / r,   f'(r) = a erfc'(a r) / r - erfc(a r) / r^2.

        Outside the cutoff both are zero.  The shifted-force Coulomb term is smooth there by
        construction; the Lennard-Jones term is smooth there only with ``lj_cutoff="force"``.
        At the default ``lj_cutoff="energy"`` it is energy-shifted alone, so its force has a
        jump at ``r = rc`` that a finite difference straddling the cutoff sees and this
        derivative does not -- measurably, in an elastic constant: see
        :class:`polyfind.mechanics.Elastic`.
        """
        xp = self.xp
        mask = (r2 < self.rc2) & (r2 > 1e-8)
        r2s = xp.where(mask, r2, 1.0)
        rr = xp.sqrt(r2s)
        ir6 = 1.0 / (r2s * r2s * r2s)
        ar = self.alpha * rr
        ef = _erfc_pos(ar, xp)
        v = ir6 * (A * ir6 - B) + qq * (ef / rr) + qqf * rr + const
        dcoul = self.alpha * _derfc_pos(ar, xp) / rr - ef / r2s
        dv = 3.0 * ir6 * (B - 2.0 * A * ir6) / r2s + (qq * dcoul + qqf) / (2.0 * rr)
        return xp.where(mask, v, 0.0), xp.where(mask, dv, 0.0)

    def _pair_phi(self, r2):
        """The charge-bilinear factor of the pair potential: ``v_coulomb = qq * phi(r)``.

        Everything in :meth:`_pair_energy` that is proportional to ``qq`` collected into one
        expression -- the damped complementary error function, the shifted-force linear term,
        and the two constants that live in ``qq_f`` and ``const``.  It is what
        ``dE/dq`` needs: with charge flux the charges are functions of the geometry, so every
        derivative of the energy picks up ``sum_a (dE/dq_a)(dq_a/dgeometry)`` and
        ``dE/dq_a = (COULOMB/eps_r) sum_b q_b phi(r_ab)`` is the electrostatic potential at
        ``a``.  Evaluated only when the flux is on; the fixed-charge kernel never calls it.
        """
        xp = self.xp
        mask = (r2 < self.rc2) & (r2 > 1e-8)
        r2s = xp.where(mask, r2, 1.0)
        rr = xp.sqrt(r2s)
        v = (_erfc_pos(self.alpha * rr, xp) / rr + self.dsf_force * rr
             - (self.dsf_shift + self.dsf_force * self.rc))
        return xp.where(mask, v, 0.0)

    def _intra_column_and_grad(self, coords, c, K: int, tables=None, q_repeat=None):
        """:meth:`_intra_column` for ONE row, plus its gradient w.r.t. the coordinates and ``c``.

        Returns ``(intra, dintra/dcoords (n, 3), dintra/dc, dintra/dq (n,) or None)``; the
        ``0.5 n_chains`` prefactor of :attr:`e_intra` is the caller's business.  ``tables``
        overrides the packer's ``(qq, qq_f, const)`` for a row carrying its own charges, and
        ``q_repeat`` (the repeat's charges) asks for the ``dq`` term; ``None`` skips it, which
        is the fixed-charge path and costs nothing.
        """
        xp = self.xp
        qq, qqf, const = tables if tables is not None else (self._qq_nn, self._qqf_nn, self._const_nn)
        kz = np.arange(-K, K + 1, dtype=float)
        ks = xp.asarray(kz, dtype=self._dt)
        cz = xp.asarray(np.asarray(c, dtype=float), dtype=self._dt)
        D = coords[:, :, None, :] - coords[:, None, :, :]  # (1, n, n, 3)
        dz = D[:, None, :, :, 2] - ks[None, :, None, None] * cz[:, None, None, None]
        r2 = D[:, None, :, :, 0] ** 2 + D[:, None, :, :, 1] ** 2 + dz * dz
        v, gv = self._pair_energy_and_dv(r2, self._A_nn, self._B_nn, qq, qqf, const)
        S = self._scale_column(K)[None]
        intra = (v * S).sum(axis=(1, 2, 3))
        gq = None
        if q_repeat is not None:
            W = np.asarray(bk.to_numpy(self._pair_phi(r2) * S), dtype=float).sum(axis=(0, 1))
            qc = np.asarray(q_repeat, dtype=float) * (COULOMB / self.eps_r)
            gq = W @ qc + W.T @ qc
        t = 2.0 * gv * S  # dV/d(separation vector) = 2 dV/dr2 * (separation vector)
        f = [np.asarray(bk.to_numpy(t * w), dtype=float) for w in (D[:, None, :, :, 0], D[:, None, :, :, 1], dz)]
        gX = np.stack([q.sum(axis=(0, 1, 3)) - q.sum(axis=(0, 1, 2)) for q in f], axis=1)
        gc = -float((f[2] * kz[None, :, None, None]).sum())
        return intra, gX, gc, gq

    def _intra_column(self, coords, c, K: int, tables=None):
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
        qq, qqf, const = tables if tables is not None else (self._qq_nn, self._qqf_nn, self._const_nn)
        v = self._pair_energy(r2, self._A_nn, self._B_nn, qq, qqf, const)
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

        When a ``field`` is set, ``-mu . E`` (:meth:`field_energy`) is added per row from
        that row's placed coordinates.  With no field the kernel is untouched.
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
            # With charge flux each row's charges follow that row's conformation, so the
            # pair tables become per-row arrays with a broadcast image axis.
            q_row = self._row_charges(coords, c_arr, M) if self._flux is not None else None
            intra = np.concatenate([
                bk.to_numpy(self._intra_column(
                    xp.asarray(coords[t : t + i_chunk], dtype=self._dt), c_arr[t : t + i_chunk], K,
                    tables=None if q_row is None else self._row_charge_tables(q_row[t : t + i_chunk], False)))
                for t in range(0, M, i_chunk)
            ])
            e_add = e_add + 0.5 * self.n_chains * intra
            if self._valence is not None:
                e_add = e_add + self.n_chains * self._valence.energy(coords, c_arr)
        else:
            c_arr = None
            q_row = None
            K, reach = self.K, self.reach
            e_add = np.full(M, self.e_torsion + self.e_intra)
            if self._valence is not None:
                e_add = e_add + self.e_valence
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
            t_nn = self._row_charge_tables(q_row[sl], False) if q_row is not None else None
            t_cell = self._row_charge_tables(q_row[sl], True) if q_row is not None else None
            qc_xp = (xp.asarray(np.tile(q_row[sl], (1, self.n_chains)), dtype=self._dt)
                     if q_row is not None else None)
            e = xp.zeros(sub.shape[0], dtype=self._dt)
            if self.n_chains == 2:
                # (0, 0, k) column: chain1-chain2 only; the (2,1) block equals the (1,2) block
                Dh = P[:, :n, None, :] - P[:, None, n:, :]  # (m, n, n, 3)
                ks = xp.asarray(np.arange(-K, K + 1, dtype=float), dtype=self._dt)
                cz = lat[:, 2, 2]
                dzh = Dh[:, None, :, :, 2] - ks[None, :, None, None] * cz[:, None, None, None]
                r2h = Dh[:, None, :, :, 0] ** 2 + Dh[:, None, :, :, 1] ** 2 + dzh * dzh
                qq, qqf, cst = t_nn if t_nn is not None else (self._qq_nn, self._qqf_nn, self._const_nn)
                vh = self._pair_energy(r2h, self._A_nn, self._B_nn, qq, qqf, cst)
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
                qq, qqf, cst = t_cell if t_cell is not None else (self._qq, self._qqf, self._const)
                v = self._pair_energy(r2, self._A, self._B, qq, qqf, cst)
                e = e + v.sum(axis=(1, 2, 3))
            if self._field_on:
                # -mu . E, from the placed coordinates of this chunk: the dipole follows
                # the setting angles and the flip, so it is a per-row quantity.
                q_use = self._q_cell_xp[None, :, None] if qc_xp is None else qc_xp[:, :, None]
                mu = (P * q_use).sum(axis=1)  # (m, 3), e.A
                out[s : s + m_chunk] = bk.to_numpy(0.5 * e - EV_TO_KCAL * (mu * self._field_xp[None, :]).sum(axis=1))
            else:
                out[s : s + m_chunk] = bk.to_numpy(0.5 * e)
        res = np.empty_like(out)
        res[order] = out
        return res + e_add

    def energy_and_grad(self, params, coords=None, c=None, e_torsion=None):
        """Energy of ONE configuration and its analytic gradient, in a single pass.

        Returns ``(E, g_cell, g_coords, g_c)``:

        ``E``
            the lattice energy of this configuration -- bit for bit the number
            :meth:`energy` returns for the same row, because the value is built from the
            same expressions in the same order (``tests/test_pack.py`` asserts equality
            with ``==``, not a tolerance).
        ``g_cell``
            ``dE/d(a, b, gamma, phi1, phi2, dz)``, per A for the lengths and the shift and
            per *degree* for the angles (the units the optimiser's variables are in).
            ``flip`` is discrete and has no entry.
        ``g_coords``
            ``dE/d(chain coordinates)``, ``(n, 3)``, including the intra-chain
            ``(0, 0, k)`` term that :attr:`e_intra` carries.  This is the quantity a
            refinement over torsions needs: the chain rule from here to any
            conformational parametrisation involves no further kernel evaluation.
        ``g_c``
            ``dE/dc``, the repeat along z, which enters both the lattice vector and the
            ``(0, 0, k)`` image shifts.

        The derivatives are closed-form throughout.  Each pair term contributes
        ``2 (dV/dr^2) * (separation vector)`` to the two atoms it joins (and, for an
        image, minus that to the image shift); from there

        * ``a``, ``b`` and ``gamma`` enter through the lattice vectors and through
          chain 2's ``(avec + bvec) / 2`` offset,
        * ``dz`` through chain 2's offset alone,
        * ``phi1`` and ``phi2`` through ``dr/dphi = z_hat x r``, r measured from each
          chain's own axis,
        * the applied field through ``-mu . E`` with ``mu = sum_i q_i r_i``, whose
          coordinate gradient is just ``-q_i E``.

        ``coords`` / ``c`` / ``e_torsion`` override the packer's chain exactly as in
        :meth:`energy`; ``e_torsion`` only shifts the value, never the gradient.
        """
        xp, dt = self.xp, self._dt
        p = np.asarray(params, dtype=float).reshape(-1)
        if p.shape[0] != 7:
            raise ValueError("energy_and_grad evaluates one configuration: params must have 7 entries")
        self.n_grad_calls += 1
        self.n_energy_calls += 1
        self.n_energy_rows += 1
        n, N, nc = self.n, self.N, self.n_chains
        sub = p[None]
        per_row = coords is not None
        if per_row:
            Xn = np.asarray(coords, dtype=float).reshape(1, n, 3)
            cz = float(self.chain.c if c is None else np.asarray(c, dtype=float).ravel()[0])
            K = int(np.ceil(self.rc / cz)) + 1
            reach = self.rc + 2.0 * float(np.linalg.norm(Xn[:, :, :2], axis=2).max()) + 0.5
            e_t = self.e_torsion if e_torsion is None else float(np.asarray(e_torsion, dtype=float).ravel()[0])
            c_arr = np.array([cz])
            X, intra_in = Xn, xp.asarray(Xn, dtype=dt)
        else:
            X, c_arr = None, None
            cz = float(self.chain.c)
            K, reach, e_t = self.K, self.reach, self.e_torsion
            intra_in = self.X0[None]
        # Charge flux: this row's charges and their exact derivatives with respect to the
        # repeat's coordinates and to c.  ``gq`` accumulates dE/dq alongside the usual
        # dE/dX, and the two are combined at the end, so ``g_coords`` and ``g_c`` come back
        # as *total* derivatives -- which is what lets every consumer's chain rule through
        # (coords, c) stay exactly as it was.
        flux_q = flux_dq = flux_dqc = None
        t_nn = t_cell = None
        q_cell = self._q_cell
        if self._flux is not None:
            Xf = Xn[0] if per_row else np.asarray(self.chain.coords, dtype=float)
            flux_q, flux_dq, flux_dqc = self._flux.charges_and_grad(Xf, cz)
            if per_row:
                t_nn = self._row_charge_tables(flux_q[None], False)
                t_cell = self._row_charge_tables(flux_q[None], True)
                q_cell = np.tile(flux_q, nc)
        intra, gX, gc, gq = self._intra_column_and_grad(
            intra_in, np.array([cz]), K, tables=t_nn, q_repeat=flux_q)
        pref = 0.5 * nc
        # identical arithmetic to ``energy``'s two branches, so the value matches bit for bit
        e_add = (e_t + pref * float(bk.to_numpy(intra)[0])) if per_row else (self.e_torsion + self.e_intra)
        gX, gc = pref * gX, pref * gc
        gq_chain = None if gq is None else pref * gq
        if self._valence is not None:
            # One repeat's valence energy per chain, with its exact derivatives with respect
            # to the repeat's coordinates and to c.  This is the whole of the plumbing the
            # deformable path needs: every consumer of ``g_coords``/``g_c`` -- the
            # refinement's chain rule, the axial relaxation -- now sees a restoring force.
            Ev, gXv, gcv = self._valence.energy_and_grad(
                Xn if per_row else self.chain.coords, cz, n_atoms=n)
            # ``self.e_valence`` rather than ``nc * Ev`` off the per-row path, so that the
            # value stays bit-for-bit the one :meth:`energy` returns for the same row.
            e_add = e_add + (nc * Ev if per_row else self.e_valence)
            gX, gc = gX + nc * gXv, gc + nc * gcv

        P, lat = self._place(sub, X, c_arr)
        latn = np.asarray(bk.to_numpy(lat), dtype=float)[0]
        Pn = np.asarray(bk.to_numpy(P), dtype=float)[0]
        gP = np.zeros((N, 3))
        e = xp.zeros(1, dtype=dt)
        kz = np.arange(-K, K + 1, dtype=float)
        if nc == 2:
            # (0, 0, k) column: chain1-chain2 only, counted twice (the (2,1) block equals (1,2))
            Dh = P[:, :n, None, :] - P[:, None, n:, :]
            ks = xp.asarray(kz, dtype=dt)
            czl = lat[:, 2, 2]
            dzh = Dh[:, None, :, :, 2] - ks[None, :, None, None] * czl[:, None, None, None]
            r2h = Dh[:, None, :, :, 0] ** 2 + Dh[:, None, :, :, 1] ** 2 + dzh * dzh
            qq, qqf, cst = t_nn if t_nn is not None else (self._qq_nn, self._qqf_nn, self._const_nn)
            vh, gh = self._pair_energy_and_dv(r2h, self._A_nn, self._B_nn, qq, qqf, cst)
            e = e + 2.0 * vh.sum(axis=(1, 2, 3))
            if gq_chain is not None:
                # Chain 1 sees this block as rows, chain 2 as columns; both carry the *same*
                # repeat charges, so both roles add into the one dE/dq per repeat atom.
                W = np.asarray(bk.to_numpy(self._pair_phi(r2h)), dtype=float).sum(axis=(0, 1))
                qc = flux_q * (COULOMB / self.eps_r)
                gq_chain = gq_chain + W @ qc + W.T @ qc
            # E gets 0.5 * (2 * sum vh) = sum vh, so the prefactor on dV/dw is exactly 1
            t = 2.0 * gh
            f = [np.asarray(bk.to_numpy(t * w), dtype=float)
                 for w in (Dh[:, None, :, :, 0], Dh[:, None, :, :, 1], dzh)]
            for d, q in enumerate(f):
                gP[:n, d] += q.sum(axis=(0, 1, 3))
                gP[n:, d] -= q.sum(axis=(0, 1, 2))
            gc += -float((f[2] * kz[None, :, None, None]).sum())
            del f, t, vh, gh, r2h, dzh, Dh
        glat = np.zeros((3, 3))
        ijk_np, L = self._select_images(sub, K, reach)
        if L:
            ijk = xp.asarray(ijk_np, dtype=dt)
            D = P[:, :, None, :] - P[:, None, :, :]
            shift = xp.einsum("mic,mcd->mid", ijk, lat)
            dx = D[:, None, :, :, 0] - shift[:, :, None, None, 0]
            dy = D[:, None, :, :, 1] - shift[:, :, None, None, 1]
            dzz = D[:, None, :, :, 2] - shift[:, :, None, None, 2]
            r2 = dx * dx + dy * dy + dzz * dzz
            qq, qqf, cst = t_cell if t_cell is not None else (self._qq, self._qqf, self._const)
            v, gv = self._pair_energy_and_dv(r2, self._A, self._B, qq, qqf, cst)
            e = e + v.sum(axis=(1, 2, 3))
            if gq_chain is not None:
                # E gets 0.5 * sum v; the image set is closed under negation, so the (a, b)
                # and (b, a) sums are equal and the 0.5 cancels against counting both.
                Wc = np.asarray(bk.to_numpy(self._pair_phi(r2)), dtype=float).sum(axis=(0, 1))
                qq_c = q_cell * (COULOMB / self.eps_r)
                gq_cell = 0.5 * (Wc @ qq_c + Wc.T @ qq_c)
                gq_chain = gq_chain + gq_cell[:n] + (gq_cell[n:] if nc == 2 else 0.0)
            del v, r2
            # E gets 0.5 * sum v, so dE/d(separation) = 0.5 * 2 * gv * separation = gv * separation
            gsh = np.zeros((ijk_np.shape[1], 3))
            for d, w in enumerate((dx, dy, dzz)):
                q = np.asarray(bk.to_numpy(gv * w), dtype=float)
                gP[:, d] += q.sum(axis=(0, 1, 3)) - q.sum(axis=(0, 1, 2))
                gsh[:, d] = -q.sum(axis=(0, 2, 3))
            glat = ijk_np[0].T @ gsh
            del dx, dy, dzz, gv, q, gsh

        if self._field_on:
            q_use = self._q_cell_xp if t_cell is None else xp.asarray(q_cell, dtype=dt)
            mu = (P * q_use[None, :, None]).sum(axis=1)
            out = float(bk.to_numpy(0.5 * e - EV_TO_KCAL * (mu * self._field_xp[None, :]).sum(axis=1))[0])
            gP += -EV_TO_KCAL * q_cell[:, None] * np.asarray(self.field, dtype=float)[None, :]
            if gq_chain is not None:
                # d(-mu . E)/dq_a = -(r_a . E): the field term's charge derivative, which is
                # what carries the *intrinsic* piezoelectric response of a fluxing chain.
                w = -EV_TO_KCAL * (Pn @ np.asarray(self.field, dtype=float))
                gq_chain = gq_chain + w[:n] + (w[n:] if nc == 2 else 0.0)
        else:
            out = float(bk.to_numpy(0.5 * e)[0])
        if gq_chain is not None:
            gX = gX + np.einsum("a,abd->bd", gq_chain, flux_dq)
            gc = gc + float(gq_chain @ flux_dqc)

        # --- chain rule: placed coordinates -> cell parameters and chain coordinates ---
        a, b, gam, phi1, phi2, dz, flip = p
        g = np.deg2rad(gam)
        cg, sg = np.cos(g), np.sin(g)
        g1, g2 = gP[:n], gP[n:]
        P1 = Pn[:n]
        rad = np.pi / 180.0
        g_cell = np.zeros(6)
        g_cell[0] = glat[0, 0]
        g_cell[1] = cg * glat[1, 0] + sg * glat[1, 1]
        g_cell[2] = rad * b * (-sg * glat[1, 0] + cg * glat[1, 1])
        # dr/dphi = z_hat x r about each chain's own axis (which passes through the origin
        # of that chain's own coordinates, i.e. before the cell offset is added)
        g_cell[3] = rad * float((g1[:, 0] * -P1[:, 1] + g1[:, 1] * P1[:, 0]).sum())
        g_c = gc + glat[2, 2]
        # rotation back out of the placed frame: Rz(phi)^T, and the flip's mirror
        ca, sa = np.cos(np.deg2rad(phi1)), np.sin(np.deg2rad(phi1))
        gX[:, 0] += ca * g1[:, 0] + sa * g1[:, 1]
        gX[:, 1] += -sa * g1[:, 0] + ca * g1[:, 1]
        gX[:, 2] += g1[:, 2]
        if nc == 2:
            off = 0.5 * (latn[0] + latn[1]) + np.array([0.0, 0.0, dz])
            Q2 = Pn[n:] - off[None, :]
            G2 = g2.sum(axis=0)
            g_cell[0] += 0.5 * G2[0]
            g_cell[1] += 0.5 * (cg * G2[0] + sg * G2[1])
            g_cell[2] += rad * 0.5 * b * (-sg * G2[0] + cg * G2[1])
            g_cell[4] = rad * float((g2[:, 0] * -Q2[:, 1] + g2[:, 1] * Q2[:, 0]).sum())
            g_cell[5] = G2[2]
            ca, sa = np.cos(np.deg2rad(phi2)), np.sin(np.deg2rad(phi2))
            gf = np.stack([ca * g2[:, 0] + sa * g2[:, 1], -sa * g2[:, 0] + ca * g2[:, 1], g2[:, 2]], axis=1)
            if flip > 0.5:  # chain 2 carries (x, -y, -z) of the chain, an involution
                gf[:, 1] *= -1.0
                gf[:, 2] *= -1.0
            gX += gf
        return out + e_add, g_cell, gX, g_c

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
        cell = np.asarray(bk.to_numpy(P[0]), dtype=float)
        mu = pol = None
        e_field = 0.0
        if self.is_neutral:
            mu = np.einsum("n,nc->c", self._q_cell, cell)
            pol = mu / float(self.cell_volume(params[None])[0]) * E_PER_A2_TO_C_PER_M2
            if self._field_on:
                e_field = float(-(mu @ self.field) * EV_TO_KCAL)
        return PackResult(
            chain=self.chain.name,
            energy_per_cell=e,
            energy_per_monomer=e / (self.n_chains * self.chain.n_monomers),
            a=float(a), b=float(b), gamma=float(gam), c=self.chain.c,
            phi1=float(phi1 % 360), phi2=float(phi2 % 360), dz=float(dz % self.chain.c), flip=int(round(flip)),
            n_chains=self.n_chains,
            density=self.density(a, b, gam),
            dihedrals=self.chain.dihedrals.copy(),
            cell_coords=cell,
            cell_elements=list(self.chain.elements) * self.n_chains,
            dipole=mu,
            polarization=pol,
            field=(0.0, 0.0, 0.0) if self.field is None else tuple(float(x) for x in self.field),
            field_energy=e_field,
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

# Grid of the chain-pair table built for the ``screen="table"`` path.  It is coarser than
# :meth:`polyfind.lattice_table.PairTable.build`'s own default (5 deg, c/16, 0.05 A),
# which is sized for *evaluating* lattice energies; here the table only has to say which
# cells are worth polishing, and the polish is exact.  Measured on PE all-trans and beta-,
# alpha- and gamma-PVDF: this grid selects starts that polish to the same global minimum
# as the fine one, for about 1/10 of the build cost (on six worker processes: PE 1.9 s
# against 18 s, gamma 8.2 s against 87 s) and 1/9 of the memory.  Pass ``table_kw`` to
# override it.
SCREEN_TABLE = {"n_angle": 48, "n_z": 8, "dr": 0.1}


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


def cell_value_and_grad_analytic(packer: CrystalPacker, x: np.ndarray, free):
    """Energy and the analytic gradient w.r.t. the free cell parameters, from ONE row.

    The same interface as :func:`cell_value_and_grad` (minus the bounds, which only
    existed to keep a finite difference from stepping outside them), backed by
    :meth:`CrystalPacker.energy_and_grad`: one configuration per function evaluation
    instead of ``1 + 2 n_free``.
    """
    x = np.asarray(x, dtype=float)
    E, g_cell, _, _ = packer.energy_and_grad(x)
    return float(E), g_cell[np.asarray(list(free), dtype=int)]


def polish(
    packer: CrystalPacker,
    params: np.ndarray,
    lo: np.ndarray,
    hi: np.ndarray,
    free: list[int],
    maxfev: int = 1500,
    method: str = "lbfgs",
    maxiter: int = 200,
    gradient: str = "analytic",
) -> np.ndarray:
    """Local minimisation of the continuous cell parameters (flip fixed).

    ``method="lbfgs"`` (default) runs L-BFGS-B with ``gradient="analytic"``
    (:meth:`CrystalPacker.energy_and_grad`: one kernel row per function evaluation) or
    ``gradient="fd"`` (the original batched central differences: one call of
    ``1 + 2 n_free`` rows).  Either way that is far fewer sequential evaluations than
    Nelder-Mead needs; ``method="nelder-mead"`` is the original simplex search with
    bounds by penalty.
    """
    x0 = np.asarray(params, dtype=float).copy()
    if method not in ("lbfgs", "l-bfgs-b", "nelder-mead", "nelder_mead"):
        raise ValueError(f"unknown polish method {method!r}")
    if gradient not in ("analytic", "fd"):
        raise ValueError(f"unknown polish gradient {gradient!r} (expected 'analytic' or 'fd')")
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
        if gradient == "analytic":
            return cell_value_and_grad_analytic(packer, x, free)
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


def _distinct_starts(params: np.ndarray, energies: np.ndarray, n: int, min_separation: float = 0.5,
                     sort_ab: bool = False, tags=None) -> list[np.ndarray]:
    """The ``n`` lowest rows of ``params`` that are not near-duplicates of one another.

    Two cells count as duplicates when ``|da| + |db| <= min_separation`` and they share a
    tag (the flip, plus gamma when several are screened).  ``sort_ab`` compares the
    *sorted* pair, so that ``(a, b)`` and ``(b, a)`` -- the same lattice rotated by
    90 deg -- are one cell.
    """
    order = np.argsort(energies)
    starts, keys = [], []
    for idx in order:
        p = params[idx]
        ab = (min(p[0], p[1]), max(p[0], p[1])) if sort_ab else (p[0], p[1])
        k = ab + (p[6] if tags is None else tags[idx],)
        if all(abs(k[0] - q[0]) + abs(k[1] - q[1]) > min_separation or k[2] != q[2] for q in keys):
            starts.append(p)
            keys.append(k)
        if len(starts) >= n:
            break
    return starts


def _random_starts(packer, lo, hi, n_random, n_refine, flips, n_chains, rng, verbose):
    """The original screen: ``n_random`` uniform cells per flip through the exact kernel."""
    cont = lo + (hi - lo) * rng.random((n_random, 6))
    if n_chains == 1:
        cont[:, 4] = 0.0
        cont[:, 5] = 0.0
    Es, Ps = [], []
    for f in flips:
        params = np.concatenate([cont, np.full((n_random, 1), float(f))], axis=1)
        Es.append(packer.energy(params))
        Ps.append(params)
    E, params = np.concatenate(Es), np.concatenate(Ps)
    if verbose:
        print(f"coarse search: {len(E)} cells, best E/cell = {E.min():.3f}")
    return _distinct_starts(params, E, n_refine)


def _table_starts(packer, chain, lo, hi, n_refine, flips, step, gammas, table, cutoff, alpha, eps_r,
                  cache_dir, table_kw, verbose):
    """The exhaustive screen: every ``(phi1, phi2, dz)`` of a tabulated ``W``, per ``(a, b)``.

    One inverse FFT per ``(a, b, flip)`` covers the whole angle/z landscape, so the grid
    the starts come from is the full search space at the table's resolution rather than a
    random sample of it.  The starts are the lowest grid cells that are not near-duplicates
    of one another, by the same rule the random screen uses (an (a, b) landscape minimum
    rather than plain energy order was tried and is worse: it promotes shallow far-away
    basins ahead of the near-degenerate deep ones that matter).  No exact-kernel evaluation
    happens here: the selected cells are handed to the same :func:`polish` as the random
    screen, so the returned energies are exact either way.

    The tabulated interaction knows nothing about an applied field, so with one set the
    ``(phi1, phi2, dz)`` minimisation inside each ``(a, b, flip)`` cell is still the
    *field-free* one.  Four times as many candidates are then taken from the screen and
    re-ranked with their exact ``-mu . E`` before the starts are chosen, and the polish
    that follows is fully field-aware -- but a field large enough to reorder the setting
    angles within a cell can still hide a basin from this screen, and ``screen="random"``
    (whose every evaluation goes through the field-aware kernel) is the honest choice
    there.

    **The table is a rigid-chain, energy-shifted object, and this is where that is
    enforced.**  ``W(r, alpha1, alpha2, dz, flip)`` is tabulated for *one* set of chain
    coordinates and re-used at every ``(a, b)``, and
    :func:`polyfind.lattice_table.chain_pair_energy` rebuilds the pair potential from
    ``packer.lj_shift`` and the two DSF constants -- it knows nothing about a force-shifted
    Lennard-Jones term and nothing about valence energy.  So a packer carrying either is
    refused outright, and so is a table built for coordinates other than the chain's.  Both
    are the same rule: **the table serves the screen and the screen is rigid; deformation
    belongs to the direct kernel** (:meth:`CrystalPacker.energy_and_grad`), which is what
    :mod:`polyfind.refine` and :mod:`polyfind.mechanics` use.  Nothing here quietly averages
    a deformed chain into a rigid tabulation.
    """
    from .lattice_table import fft_screen, pair_table

    if packer.valence is not None or packer.lj_cutoff != "energy" or packer._flux is not None:
        what = ("valence terms" if packer.valence is not None else
                "charge flux" if packer._flux is not None else f"lj_cutoff={packer.lj_cutoff!r}")
        raise ValueError(
            f"screen='table' cannot be used with a packer carrying {what}: the tabulated chain-pair "
            "interaction is built for a rigid chain with fixed charges from the energy-shifted pair potential, so it "
            "would not be the potential being polished.  Use screen='random' (every evaluation goes "
            "through the exact kernel), or screen a rigid, energy-shifted packer and hand the result "
            "to the direct-kernel path (refine_crystal, polyfind.mechanics) for the deformable part."
        )
    if table is not None and (table.chain.coords.shape != chain.coords.shape
                              or float(np.abs(table.chain.coords - chain.coords).max()) > 1e-6):
        raise ValueError(
            "the given table was built for different chain coordinates than the chain being packed: "
            "W is tabulated for one rigid conformation and says nothing about a deformed one"
        )
    if table is None:
        table = pair_table(chain, cutoff=cutoff, alpha=alpha, eps_r=eps_r, cache_dir=cache_dir,
                           verbose=verbose, **{**SCREEN_TABLE, **(table_kw or {})})
    elif abs(table.c - chain.c) > 1e-6 or abs(table.rc - cutoff) > 1e-9:
        raise ValueError(f"the given table is for c={table.c:.4f}, cutoff={table.rc}, "
                         f"not c={chain.c:.4f}, cutoff={cutoff}")
    a_values = np.arange(lo[0], hi[0] + 1e-9, step)
    b_values = np.arange(lo[1], hi[1] + 1e-9, step)
    # (a, b) and (b, a) are the same lattice rotated by 90 deg, so half the grid is
    # redundant -- but only when both axes are screened over the same values, otherwise
    # skipping b < a would drop cells whose mirror image is not on the grid at all.
    same_grid = a_values.shape == b_values.shape and np.allclose(a_values, b_values)
    n_top = 4 * n_refine if packer._field_on else n_refine
    rows, es, tags = [], [], []
    for g in gammas:
        res = fft_screen(table, a_values, b_values, gamma=float(g), flips=tuple(flips), n_top=n_top,
                         min_separation=0.5, ab_symmetry=same_grid and abs(float(g) - 90.0) < 1e-9)
        rows.append(res.top)
        es.append(res.top_energy + packer.field_energy(res.top) if packer._field_on else res.top_energy)
        tags += [(float(p[6]), round(float(g), 6)) for p in res.top]
        if verbose:
            n_pts = len(a_values) * len(b_values) * len(flips) * table.n_angle ** 2 * table.n_z
            print(f"table screen gamma={g:.1f}: {len(a_values)}x{len(b_values)} cells, {n_pts / 1e6:.0f} M "
                  f"landscape points in {res.time:.2f} s, best E/cell = {res.top_energy[0]:.3f}")
    params, E = np.concatenate(rows), np.concatenate(es)
    return _distinct_starts(params, E, n_refine, sort_ab=True, tags=tags)


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
    screen: str = "table",
    alpha: float = 0.2,
    screen_step: float = 0.25,
    screen_gammas=None,
    table=None,
    table_cache_dir: str | None = "env",
    table_kw: dict | None = None,
    field=None,
    gradient: str = "analytic",
) -> list[PackResult]:
    """Coarse screen of the cell parameters followed by local polishing of the best cells.

    ``screen`` selects how the starts are found:

    * ``"table"`` (default) -- the exhaustive screen of :mod:`polyfind.lattice_table`.
      A tabulated chain-pair interaction ``W`` (built once per conformation and cached by
      :func:`~polyfind.lattice_table.pair_table`) turns each ``(a, b, flip)`` into one
      inverse FFT that evaluates the *entire* ``(phi1, phi2, dz)`` landscape, so the screen
      covers 10^7-10^8 cells in a second or two instead of sampling ``2 n_random`` of them
      with the exact kernel, and no basin can be missed by bad luck.  ``screen_step`` is the
      ``(a, b)`` grid spacing in A; ``screen_gammas`` (default: the gamma bounds, five
      values when they are free) the gammas screened; ``table`` accepts a prebuilt
      :class:`~polyfind.lattice_table.PairTable`; ``table_cache_dir`` and ``table_kw``
      (defaulting to :data:`SCREEN_TABLE`) are passed to the cached builder.  The build is
      the largest single cost of this path -- 1.6-8 s per conformation, on as many worker
      processes as the machine has physical cores -- and it is paid once: the table depends
      only on the chain and the potential, so a second ``pack()`` of the same conformation
      is free, and setting ``$POLYFIND_TABLE_CACHE`` (or ``table_cache_dir``) shares it
      with other processes and later runs.  Only for ``n_chains == 2``; a one-chain cell
      falls back to ``"random"``.
    * ``"random"`` -- ``n_random`` uniform random cells per flip through the exact kernel.

    Either way the starts are polished with the exact kernel, so the returned energies,
    the return type, the ordering and the deduplication of near-identical minima are the
    same.  ``method`` selects the polisher: ``"lbfgs"`` (default, gradient L-BFGS-B) or
    ``"nelder-mead"``; ``gradient`` selects where an L-BFGS-B polish gets its gradient,
    ``"analytic"`` (default, one kernel row per evaluation) or ``"fd"`` (batched central
    differences, ``1 + 2 n_free`` rows).

    ``field=(Ex, Ey, Ez)`` (V/A) applies a uniform electric field: the cell is screened
    and polished against the lattice energy *plus* ``-mu_cell . E``, so what comes back
    is the cell the field selects, and every :class:`PackResult` carries its dipole and
    polarization.  ``screen="random"`` is field-aware throughout; ``screen="table"``
    screens field-free and re-ranks (see :func:`_table_starts`).
    """
    if screen not in ("table", "random"):
        raise ValueError(f"unknown screen {screen!r} (expected 'table' or 'random')")
    rng = rng or np.random.default_rng(0)
    packer = CrystalPacker(chain, n_chains=n_chains, cutoff=cutoff, alpha=alpha, eps_r=eps_r, field=field)
    bounds = bounds or default_bounds(chain, gamma_free)
    keys = ["a", "b", "gamma", "phi1", "phi2", "dz"]
    lo = np.array([bounds[k][0] for k in keys])
    hi = np.array([bounds[k][1] for k in keys])
    flips = list(flips) if n_chains == 2 else [0]
    if screen == "table" and n_chains == 2:
        if screen_gammas is None:
            screen_gammas = [lo[2]] if hi[2] <= lo[2] + 1e-9 else np.linspace(lo[2], hi[2], 5)
        starts = _table_starts(packer, chain, lo, hi, n_refine, flips, screen_step, screen_gammas, table,
                               cutoff, alpha, eps_r, table_cache_dir, table_kw, verbose)
    else:
        starts = _random_starts(packer, lo, hi, n_random, n_refine, flips, n_chains, rng, verbose)
    free = [i for i in range(6) if hi[i] > lo[i] and not (n_chains == 1 and i in (4, 5))]
    results = [packer.result(polish(packer, p, lo, hi, free, maxfev, method=method, gradient=gradient)) for p in starts]
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
