"""Energy backends and RIS parameter fitting.

Two roles:

1. A pluggable :class:`Calculator` protocol.  :class:`SimpleFF` is a small,
   transparent intramolecular potential (Fourier torsion on the backbone
   dihedral + UFF Lennard-Jones + Coulomb for 1-4 and beyond) that exists so the
   whole pipeline runs and can be tested without heavy dependencies.  Its
   parameters are *illustrative*.  :class:`ASECalculator` wraps any ASE
   calculator (e.g. a MACE machine-learned potential) behind the same
   interface, which is the intended production path.

2. :func:`fit_ris`: derive the RIS first- and second-order energies from a
   calculator by scanning one and two consecutive backbone dihedrals of a short
   oligomer.  This is the *only* place the (possibly expensive) potential is
   evaluated densely: B*(360/step) + B*(360/step)^2 single points, about 2.6k
   for step = 10 deg and B = 2.  The scan is embarrassingly parallel and goes
   through ``energy_batch`` so a GPU-backed calculator can evaluate all
   conformers at once.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Protocol, Sequence

import numpy as np

from . import backend as bk
from .chain import Structure, build_chain, build_chain_batch
from .polymers import Polymer, RISStates, THREE_STATE, lj_params
from .ris import RISModel

COULOMB = 332.0637  # kcal A / (mol e^2)

# A batch of conformers sharing one topology: (template, coords) with coords shape (M, n, 3),
# as produced by :func:`polyfind.chain.build_chain_batch`.
ConformerBatch = tuple[Structure, np.ndarray]


class Calculator(Protocol):
    """Energies in kcal/mol for one or many rigid-geometry oligomers.

    ``energy_batch`` accepts either a list of :class:`~polyfind.chain.Structure`
    (independent topologies) or a single :data:`ConformerBatch` ``(template, coords)``
    sharing one topology -- the form produced by :func:`~polyfind.chain.build_chain_batch`.
    The tuple form lets a GPU-resident or otherwise batched calculator (an MLIP, say)
    evaluate every row as one forward pass instead of looping; :class:`SimpleFF` does
    this via :meth:`SimpleFF.energy_coords`, and :class:`ASECalculator` loops but still
    accepts the tuple form so callers need not branch on it.
    """

    def energy(self, struct: Structure) -> float: ...

    def energy_batch(self, structs: Sequence[Structure] | ConformerBatch) -> np.ndarray: ...


def erfc_approx(x, xp=np):
    """Abramowitz-Stegun 7.1.26 (|err| < 1.5e-7), elementary ops only (runs on CuPy)."""
    ax = xp.abs(x)
    t = 1.0 / (1.0 + 0.3275911 * ax)
    poly = t * (0.254829592 + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))))
    y = poly * xp.exp(-ax * ax)
    return xp.where(x >= 0, y, 2.0 - y)


@dataclass
class _Topology:
    pairs_i: np.ndarray
    pairs_j: np.ndarray
    lj_x: np.ndarray  # per pair
    lj_d: np.ndarray
    qq: np.ndarray  # per pair q_i q_j * scale
    torsions: np.ndarray  # (n_tors, 4) backbone quadruples


@dataclass
class SimpleFF:
    """Illustrative intramolecular potential (kcal/mol).

    torsion = (V1, V2, V3) Fourier terms on each backbone C-C-C-C dihedral:
        V(phi) = 1/2 [V1 (1 + cos phi) + V2 (1 - cos 2 phi) + V3 (1 + cos 3 phi)]
    Nonbonded (LJ 12-6 with UFF radii/well depths, geometric combining; Coulomb with
    relative permittivity eps_r) for all pairs separated by >= 3 bonds, with 1-4 pairs
    scaled by ``scale14``.
    """

    torsion: tuple[float, float, float] = (1.3, -0.05, 2.5)
    scale14: float = 0.5
    eps_r: float = 1.0
    _cache: dict = field(default_factory=dict, repr=False)

    # --- topology -------------------------------------------------------------
    def _topology(self, struct: Structure) -> _Topology:
        key = (struct.polymer.name, struct.n_dihedrals, len(struct.elements))
        top = self._cache.get(key)
        if top is not None:
            return top
        n = struct.n_atoms
        adj = [[] for _ in range(n)]
        for a, b in struct.bonds:
            adj[a].append(b)
            adj[b].append(a)
        # graph distances up to 3 via BFS from each atom
        dist = np.full((n, n), 99, dtype=int)
        for s in range(n):
            dist[s, s] = 0
            frontier = [s]
            for d in (1, 2, 3):
                nxt = []
                for u in frontier:
                    for v in adj[u]:
                        if dist[s, v] == 99:
                            dist[s, v] = d
                            nxt.append(v)
                frontier = nxt
        iu, ju = np.triu_indices(n, 1)
        keep = dist[iu, ju] >= 3
        iu, ju = iu[keep], ju[keep]
        scale = np.where(dist[iu, ju] == 3, self.scale14, 1.0)
        x, d = lj_params(struct.elements)
        lj_x = np.sqrt(x[iu] * x[ju])
        lj_d = np.sqrt(d[iu] * d[ju]) * scale
        qq = struct.charges[iu] * struct.charges[ju] * scale * COULOMB / self.eps_r
        tors = np.array([struct.dihedral_atoms(j) for j in range(struct.n_dihedrals)], dtype=int).reshape(-1, 4)
        top = _Topology(iu, ju, lj_x, lj_d, qq, tors)
        self._cache[key] = top
        return top

    # --- energies -------------------------------------------------------------
    def _energy_from_coords(self, top: _Topology, coords, xp):
        """coords: (M, n, 3) on the active backend -> (M,) energies."""
        pi, pj = xp.asarray(top.pairs_i), xp.asarray(top.pairs_j)
        ri = coords[:, pi]
        rj = coords[:, pj]
        r = xp.sqrt(((ri - rj) ** 2).sum(axis=-1))
        lj_x = xp.asarray(top.lj_x)
        lj_d = xp.asarray(top.lj_d)
        qq = xp.asarray(top.qq)
        s6 = (lj_x / r) ** 6
        e_lj = (lj_d * (s6 * s6 - 2.0 * s6)).sum(axis=1)
        e_c = (qq / r).sum(axis=1)
        e_t = xp.zeros(coords.shape[0])
        if top.torsions.size:
            tors = xp.asarray(top.torsions)
            a, b, c, d = (coords[:, tors[:, k]] for k in range(4))
            b0, b1, b2 = b - a, c - b, d - c
            n1, n2 = xp.cross(b0, b1), xp.cross(b1, b2)
            x = (n1 * n2).sum(-1)
            y = xp.sqrt((b1 * b1).sum(-1)) * (b0 * n2).sum(-1)
            phi = xp.arctan2(y, x)
            V1, V2, V3 = self.torsion
            e_t = (0.5 * (V1 * (1 + xp.cos(phi)) + V2 * (1 - xp.cos(2 * phi)) + V3 * (1 + xp.cos(3 * phi)))).sum(axis=1)
        return e_lj + e_c + e_t

    def energy(self, struct: Structure) -> float:
        top = self._topology(struct)
        return float(self._energy_from_coords(top, struct.coords[None], np)[0])

    def energy_coords(self, template: Structure, coords: np.ndarray) -> np.ndarray:
        """Energies for an (M, n, 3) coordinate array using ``template``'s cached topology.

        Same maths as :meth:`_energy_from_coords`/:meth:`energy_batch`, but skips building
        one :class:`Structure` per conformer -- the array typically comes straight from
        :func:`polyfind.chain.build_chain_batch`. Runs on GPU if available.
        """
        coords = np.asarray(coords)
        if coords.shape[0] == 0:
            return np.empty(0)
        xp = bk.get_backend()
        top = self._topology(template)
        coords_xp = xp.asarray(coords, dtype=bk.float_dtype())
        out = np.empty(coords.shape[0])
        chunk = 2048
        for k in range(0, coords.shape[0], chunk):
            out[k : k + chunk] = bk.to_numpy(self._energy_from_coords(top, coords_xp[k : k + chunk], xp))
        return out

    def energy_batch(self, structs: Sequence[Structure] | ConformerBatch) -> np.ndarray:
        """Vectorised over conformers (same topology); runs on GPU if available.

        Accepts either a list of :class:`Structure` or a ``(template, coords)`` tuple
        (see :class:`Calculator`); the latter skips restacking coordinates that are
        already a batched array.
        """
        if isinstance(structs, tuple):
            template, coords = structs
            return self.energy_coords(template, coords)
        structs = list(structs)
        if not structs:
            return np.empty(0)
        coords = np.stack([s.coords for s in structs])
        return self.energy_coords(structs[0], coords)

    def components(self, struct: Structure) -> dict:
        top = self._topology(struct)
        c = struct.coords
        r = np.linalg.norm(c[top.pairs_i] - c[top.pairs_j], axis=1)
        s6 = (top.lj_x / r) ** 6
        return {
            "lj": float((top.lj_d * (s6 * s6 - 2 * s6)).sum()),
            "coulomb": float((top.qq / r).sum()),
            "torsion": self.energy(struct) - float((top.lj_d * (s6 * s6 - 2 * s6)).sum()) - float((top.qq / r).sum()),
        }


class ASECalculator:
    """Adapter: any ASE calculator (MACE, CHGNet, DFT codes...) as a polyfind Calculator.

    Energies are converted from eV to kcal/mol.  ``energy_batch`` loops by default;
    override for calculators with native batching.
    """

    EV_TO_KCAL = 23.0605

    def __init__(self, ase_calc):
        try:
            from ase import Atoms  # noqa: F401
        except ImportError as e:  # pragma: no cover
            raise ImportError("ASECalculator requires the 'ase' package") from e
        self.calc = ase_calc

    def _atoms_from(self, elements, coords):
        from ase import Atoms

        atoms = Atoms(symbols=elements, positions=coords)
        atoms.calc = self.calc
        return atoms

    def _atoms(self, struct: Structure):
        return self._atoms_from(struct.elements, struct.coords)

    def energy(self, struct: Structure) -> float:
        return float(self._atoms(struct).get_potential_energy()) * self.EV_TO_KCAL

    def energy_batch(self, structs: Sequence[Structure] | ConformerBatch) -> np.ndarray:
        """Loops (no native batching here); accepts the ``(template, coords)`` tuple form
        too, building one ``Atoms`` per row from the shared element list."""
        if isinstance(structs, tuple):
            template, coords = structs
            out = np.empty(len(coords))
            for i, c in enumerate(coords):
                atoms = self._atoms_from(template.elements, c)
                out[i] = float(atoms.get_potential_energy()) * self.EV_TO_KCAL
            return out
        return np.array([self.energy(s) for s in structs])


# ---------------------------------------------------------------------- fitting
def _basins(states: RISStates, grid: np.ndarray) -> list[np.ndarray]:
    """Grid indices belonging to each state's basin (nearest ideal angle on the circle)."""
    diff = np.abs((grid[:, None] - np.array(states.angles)[None, :] + 180.0) % 360.0 - 180.0)
    owner = diff.argmin(axis=1)
    return [np.where(owner == s)[0] for s in range(states.n)]


def _wrap180(x):
    """Wrap an angle (deg, array or scalar) to (-180, 180]."""
    return ((np.asarray(x, dtype=float) + 180.0) % 360.0) - 180.0


def _basin_bounds(states: RISStates, eps: float = 1e-6) -> list[tuple[float, float]]:
    """Per-state ``(lo, hi)`` offsets (deg) from the state's own ideal angle to the
    boundary with its nearest neighbour on each side of the circle, matching
    :func:`_basins`'s tie-break EXACTLY: ``argmin`` over ``|grid - ideal|`` hands an exact
    boundary point to the lowest-index state. So a boundary shared with a lower-indexed
    neighbour belongs to that neighbour, not this state -- this state's interval is
    nudged open on that side by ``eps``; a boundary shared with a higher-indexed
    neighbour belongs to this state, so that side stays closed. Without this, a
    continuous local search can sit exactly on a boundary the discrete (dense-grid)
    bucketing would have assigned to the other state, comparing two different basins.
    """
    angs = np.array(states.angles, dtype=float)
    bounds = []
    for s in range(states.n):
        hi, hi_sp = 180.0, None
        lo, lo_sp = -180.0, None
        for sp in range(states.n):
            if sp == s:
                continue
            d = float(_wrap180(angs[sp] - angs[s]))  # signed shortest distance a_s -> a_sp
            half = d / 2.0
            if d > 0 and half < hi:
                hi, hi_sp = half, sp
            elif d < 0 and half > lo:
                lo, lo_sp = half, sp
        if hi_sp is not None and hi_sp < s:
            hi -= eps
        if lo_sp is not None and lo_sp < s:
            lo += eps
        bounds.append((lo, hi))
    return bounds


def _clip_to_basin(x, ideal: float, lo: float, hi: float):
    """Clip angle(s) ``x`` to ``ideal``'s basin, expressed via the ``(lo, hi)`` offsets
    from :func:`_basin_bounds`; wraps back to a canonical angle."""
    off = np.clip(_wrap180(np.asarray(x, dtype=float) - ideal), lo, hi)
    return _wrap180(ideal + off)


def _quad_vertex(x0, x1, x2, y0, y1, y2, fallback):
    """Vertex of the parabola through three (possibly unequally spaced) points; falls
    back to ``fallback`` if the points are degenerate or the fit is not convex (no
    interior minimum, e.g. a flat or concave stencil)."""
    xs, ys = np.array([x0, x1, x2], dtype=float), np.array([y0, y1, y2], dtype=float)
    if len(set(np.round(xs, 9))) < 3:
        return fallback
    a, b, _ = np.polyfit(xs, ys, 2)
    if not np.isfinite(a) or a <= 1e-12:
        return fallback
    return float(-b / (2.0 * a))


def _refine_coord(polymer, calc, base, ref, j, other_j, other_val, states, bounds, s_idx, center, energy, h):
    """One coarse-to-fine refinement step along dihedral index ``j`` for a whole batch of
    1-D basins at once (one state each, flattened): evaluate ``center +/- h`` (clipped to
    each point's basin), fit a parabola per point and evaluate its vertex too (also
    clipped), and move to whichever of the four points evaluated has the lowest energy --
    so a bad fit (flat/concave stencil, or a vertex estimate that turns out worse) can
    never make things worse, only fail to help. ``other_j``/``other_val`` hold the other
    dihedral's index and (fixed, per-point) value when this is used for one direction of a
    2-D scan; ``None`` for a plain 1-D scan. Two batched ``energy_batch`` calls cover every
    basin regardless of how many there are. Returns ``(new_center, new_energy,
    n_evaluations)`` with the same shape as ``center``.
    """
    shape = center.shape
    c = center.reshape(-1)
    e = energy.reshape(-1)
    sidx = s_idx.reshape(-1)
    other = None if other_val is None else other_val.reshape(-1)
    n = c.shape[0]
    lo = np.array([_clip_to_basin(c[i] - h, states.angles[sidx[i]], *bounds[sidx[i]]) for i in range(n)])
    hi = np.array([_clip_to_basin(c[i] + h, states.angles[sidx[i]], *bounds[sidx[i]]) for i in range(n)])

    dihs = np.tile(base, (2 * n, 1))
    dihs[:n, j] = lo
    dihs[n:, j] = hi
    if other_j is not None:
        dihs[:n, other_j] = other
        dihs[n:, other_j] = other
    template, coords = build_chain_batch(polymer, dihs)
    E = calc.energy_batch((template, coords)) - ref
    e_lo, e_hi = E[:n], E[n:]

    vtx = np.array([_quad_vertex(lo[i], c[i], hi[i], e_lo[i], e[i], e_hi[i], c[i]) for i in range(n)])
    vtx = np.array([_clip_to_basin(vtx[i], states.angles[sidx[i]], *bounds[sidx[i]]) for i in range(n)])
    dihs_v = np.tile(base, (n, 1))
    dihs_v[:, j] = vtx
    if other_j is not None:
        dihs_v[:, other_j] = other
    template_v, coords_v = build_chain_batch(polymer, dihs_v)
    e_v = calc.energy_batch((template_v, coords_v)) - ref

    new_c, new_e = np.empty(n), np.empty(n)
    for i in range(n):
        candidates = [(e[i], c[i]), (e_lo[i], lo[i]), (e_hi[i], hi[i]), (e_v[i], vtx[i])]
        best_e, best_x = min(candidates, key=lambda t: t[0])
        new_c[i], new_e[i] = best_x, best_e
    return new_c.reshape(shape), new_e.reshape(shape), 3 * n


_STENCIL2D = [(dx, dy) for dx in (-1, 0, 1) for dy in (-1, 0, 1) if (dx, dy) != (0, 0)]  # 3x3 minus centre, incl. diagonals


def _refine_pair(polymer, calc, base, ref, j0, j1, states, bounds, s_idx, sp_idx, cx, cy, ce, h):
    """One coarse-to-fine pattern-search step for a whole batch of 2-D basins (state
    pairs) at once: evaluate the full 3x3 stencil (8 neighbours; the centre is already
    known) and move to whichever of the 9 points has the lowest energy.  The diagonal
    stencil directions let this track a tilted or curved valley (e.g. the pentane-effect
    G+G- basin, where the two dihedrals must move together to avoid a steric clash) in one
    step, unlike two independent 1-D scans; the plain greedy choice (no curve fit) is
    robust to the near-singular repulsive wall right next to that basin's minimum, which
    wrecked a quadratic/Newton fit (huge outlier values, degenerate points when a stencil
    arm clips onto the basin edge). One batched ``energy_batch`` call covers every basin
    regardless of how many there are. Returns ``(new_cx, new_cy, new_ce,
    n_evaluations)``.
    """
    shape = cx.shape
    x0, y0, e0 = cx.reshape(-1), cy.reshape(-1), ce.reshape(-1)
    si, pi = s_idx.reshape(-1), sp_idx.reshape(-1)
    n = x0.shape[0]
    n_off = len(_STENCIL2D)
    X = np.empty((n_off, n))
    Y = np.empty((n_off, n))
    for oi, (dx, dy) in enumerate(_STENCIL2D):
        X[oi] = [_clip_to_basin(x0[i] + dx * h, states.angles[si[i]], *bounds[si[i]]) for i in range(n)]
        Y[oi] = [_clip_to_basin(y0[i] + dy * h, states.angles[pi[i]], *bounds[pi[i]]) for i in range(n)]
    dihs = np.tile(base, (n_off * n, 1))
    for oi in range(n_off):
        dihs[oi * n : (oi + 1) * n, j0] = X[oi]
        dihs[oi * n : (oi + 1) * n, j1] = Y[oi]
    template, coords = build_chain_batch(polymer, dihs)
    E = (calc.energy_batch((template, coords)) - ref).reshape(n_off, n)

    out_x, out_y, out_e = x0.copy(), y0.copy(), e0.copy()
    for oi in range(n_off):
        better = E[oi] < out_e
        out_x[better] = X[oi][better]
        out_y[better] = Y[oi][better]
        out_e[better] = E[oi][better]
    n_eval = n_off * n
    return out_x.reshape(shape), out_y.reshape(shape), out_e.reshape(shape), n_eval


def _polish_pair(polymer, calc, base, ref, j0, j1, states, bounds, s_idx, sp_idx, cx, cy, ce, h):
    """A single quadratic-surface (Newton) polish on top of :func:`_refine_pair`'s
    pattern-search result: fit ``f0 + g.d + 1/2 d^T H d`` to the same 3x3 stencil, step to
    its stationary point if the fit is locally convex (else stay put), and keep whichever
    of the 10 points evaluated (9 stencil, already-known centre excepted, plus the Newton
    point) is lowest -- so this can only match or improve on the pattern-search centre,
    never lose it. Doing this once, after pattern search has already found the basin's
    true sub-region, is safe where iterating it from a coarse start was not (see
    :func:`_refine_pair`'s docstring): the neighbourhood is by now locally well-behaved, so
    the quadratic fit is no longer at the mercy of the repulsive wall or a boundary-clipped
    duplicate stencil point corrupting its curvature. Returns ``(cx, cy, ce,
    n_evaluations)``.
    """
    shape = cx.shape
    x0, y0, e0 = cx.reshape(-1), cy.reshape(-1), ce.reshape(-1)
    si, pi = s_idx.reshape(-1), sp_idx.reshape(-1)
    n = x0.shape[0]
    n_off = len(_STENCIL2D)
    X = np.empty((n_off, n))
    Y = np.empty((n_off, n))
    for oi, (dx, dy) in enumerate(_STENCIL2D):
        X[oi] = [_clip_to_basin(x0[i] + dx * h, states.angles[si[i]], *bounds[si[i]]) for i in range(n)]
        Y[oi] = [_clip_to_basin(y0[i] + dy * h, states.angles[pi[i]], *bounds[pi[i]]) for i in range(n)]
    dihs = np.tile(base, (n_off * n, 1))
    for oi in range(n_off):
        dihs[oi * n : (oi + 1) * n, j0] = X[oi]
        dihs[oi * n : (oi + 1) * n, j1] = Y[oi]
    template, coords = build_chain_batch(polymer, dihs)
    E = (calc.energy_batch((template, coords)) - ref).reshape(n_off, n)

    new_x, new_y = x0.copy(), y0.copy()
    for i in range(n):
        dx = np.concatenate([[0.0], X[:, i] - x0[i]])
        dy = np.concatenate([[0.0], Y[:, i] - y0[i]])
        de = np.concatenate([[e0[i]], E[:, i]])
        A = np.column_stack([np.ones_like(dx), dx, dy, dx**2, dx * dy, dy**2])
        try:
            coef, *_ = np.linalg.lstsq(A, de, rcond=None)
            _, gx, gy, hxx, hxy, hyy = coef
            H = np.array([[2 * hxx, hxy], [hxy, 2 * hyy]])
            if np.all(np.linalg.eigvalsh(H) > 1e-9):
                delta = np.linalg.solve(-H, [gx, gy])
                new_x[i] = _clip_to_basin(x0[i] + delta[0], states.angles[si[i]], *bounds[si[i]])
                new_y[i] = _clip_to_basin(y0[i] + delta[1], states.angles[pi[i]], *bounds[pi[i]])
        except np.linalg.LinAlgError:
            pass

    dihs_v = np.tile(base, (n, 1))
    dihs_v[:, j0] = new_x
    dihs_v[:, j1] = new_y
    template_v, coords_v = build_chain_batch(polymer, dihs_v)
    e_v = calc.energy_batch((template_v, coords_v)) - ref

    out_x, out_y, out_e = np.empty(n), np.empty(n), np.empty(n)
    for i in range(n):
        candidates = [(e0[i], x0[i], y0[i])]
        candidates += [(E[oi, i], X[oi, i], Y[oi, i]) for oi in range(n_off)]
        candidates.append((e_v[i], new_x[i], new_y[i]))
        best_e, best_x, best_y = min(candidates, key=lambda t: t[0])
        out_x[i], out_y[i], out_e[i] = best_x, best_y, best_e
    n_eval = n_off * n + n
    return out_x.reshape(shape), out_y.reshape(shape), out_e.reshape(shape), n_eval


def _dense_scan_bond(polymer, calc, states, base, j0, grid, basins, ref):
    """Bond-type ``b``'s share of the dense scan (batched conformer building): full
    ``len(grid)`` 1-D scan and ``len(grid)^2`` 2-D scan, basin minima on the grid exactly
    as before. Returns ``(E1, E2, e1, arg1, e2, arg2, n_eval)``."""
    S = states.n
    j1 = j0 + 1
    dihs1 = np.tile(base, (len(grid), 1))
    dihs1[:, j0] = grid
    template1, coords1 = build_chain_batch(polymer, dihs1)
    E1 = calc.energy_batch((template1, coords1)) - ref
    n_eval = len(grid)
    e1, arg1 = np.zeros(S), {}
    for s, idx in enumerate(basins):
        k = idx[np.argmin(E1[idx])]
        e1[s] = E1[k]
        arg1[states.names[s]] = float(grid[k])

    phis, psis = np.meshgrid(grid, grid, indexing="ij")
    dihs2 = np.tile(base, (len(grid) ** 2, 1))
    dihs2[:, j0] = phis.ravel()
    dihs2[:, j1] = psis.ravel()
    template2, coords2 = build_chain_batch(polymer, dihs2)
    E2 = (calc.energy_batch((template2, coords2)) - ref).reshape(len(grid), len(grid))
    n_eval += len(grid) ** 2
    e2, arg2 = np.zeros((S, S)), {}
    for s, idx in enumerate(basins):
        for sp, idxp in enumerate(basins):
            sub = E2[np.ix_(idx, idxp)]
            k = np.unravel_index(np.argmin(sub), sub.shape)
            arg2[(states.names[s], states.names[sp])] = (float(grid[idx[k[0]]]), float(grid[idxp[k[1]]]))
            e2[s, sp] = sub[k]
    return E1, E2, e1, arg1, e2, arg2, n_eval


def _refine_schedule(coarse_step: float, step: float, rounds: int) -> list[float]:
    """Half-widths for the refinement rounds: the first bridges the coarse grid's spacing
    (a coarse-grid pick can be up to ``coarse_step/2`` away from the true basin minimum,
    e.g. when a steric-clash basin's minimum sits near the basin edge -- a stencil narrower
    than that would have to extrapolate rather than interpolate to find it), then each
    round halves down to ``step`` and stays there.
    """
    return [max(step, coarse_step / (2.0 * (r + 1))) for r in range(rounds)]


def _adaptive_scan_bond(polymer, calc, states, base, j0, coarse_grid, basins_coarse, bounds, ref, step, coarse_step, rounds=2):
    """Bond-type ``b``'s share of the coarse-to-fine scan: a coarse (``coarse_step``) grid
    to find each basin's approximate minimum, batched across every state (1-D) / state
    pair (2-D) at once, then ``rounds`` of refinement (see :func:`_refine_schedule` for the
    half-width schedule): :func:`_refine_coord` (quadratic-interpolation, safe by
    construction) for the 1-D basins; for the 2-D basins the pattern-search
    :func:`_refine_pair` (robust to the steep repulsive wall right next to a
    steric-clash basin's minimum) followed by one :func:`_polish_pair` quadratic step for
    sub-step precision. A couple of batched calculator calls per round cover every basin
    simultaneously. Returns the same shape as :func:`_dense_scan_bond`, with ``E1``/``E2``
    the coarse-grid energies (for diagnostics) and ``e1``/``e2``/``arg1``/``arg2`` the
    refined basin minima."""
    S = states.n
    j1 = j0 + 1
    n_eval = 0
    schedule = _refine_schedule(coarse_step, step, rounds)

    # ---- 1-D: coarse grid, then refine every state's basin minimum
    dihs1 = np.tile(base, (len(coarse_grid), 1))
    dihs1[:, j0] = coarse_grid
    template1, coords1 = build_chain_batch(polymer, dihs1)
    E1 = calc.energy_batch((template1, coords1)) - ref
    n_eval += len(coarse_grid)
    center1, energy1, s_idx1 = np.empty(S), np.empty(S), np.arange(S)
    for s, idx in enumerate(basins_coarse):
        k = idx[np.argmin(E1[idx])]
        center1[s], energy1[s] = coarse_grid[k], E1[k]
    for h in schedule:
        center1, energy1, ne = _refine_coord(polymer, calc, base, ref, j0, None, None, states, bounds, s_idx1, center1, energy1, h)
        n_eval += ne
    arg1 = {states.names[s]: float(center1[s]) for s in range(S)}
    e1 = energy1

    # ---- 2-D: coarse grid, then refine every state-pair's basin minimum (phi then psi)
    phis, psis = np.meshgrid(coarse_grid, coarse_grid, indexing="ij")
    dihs2 = np.tile(base, (len(coarse_grid) ** 2, 1))
    dihs2[:, j0] = phis.ravel()
    dihs2[:, j1] = psis.ravel()
    template2, coords2 = build_chain_batch(polymer, dihs2)
    E2 = (calc.energy_batch((template2, coords2)) - ref).reshape(len(coarse_grid), len(coarse_grid))
    n_eval += len(coarse_grid) ** 2
    cx, cy, ce = np.empty((S, S)), np.empty((S, S)), np.empty((S, S))
    s_idx_x, s_idx_y = np.empty((S, S), dtype=int), np.empty((S, S), dtype=int)
    for s, idx in enumerate(basins_coarse):
        for sp, idxp in enumerate(basins_coarse):
            sub = E2[np.ix_(idx, idxp)]
            k = np.unravel_index(np.argmin(sub), sub.shape)
            cx[s, sp], cy[s, sp], ce[s, sp] = coarse_grid[idx[k[0]]], coarse_grid[idxp[k[1]]], sub[k]
            s_idx_x[s, sp], s_idx_y[s, sp] = s, sp
    for h in schedule:
        cx, cy, ce, ne = _refine_pair(polymer, calc, base, ref, j0, j1, states, bounds, s_idx_x, s_idx_y, cx, cy, ce, h)
        n_eval += ne
    # one quadratic polish for sub-step precision, now that pattern search has found the
    # right sub-region (safe: see _polish_pair)
    cx, cy, ce, ne = _polish_pair(polymer, calc, base, ref, j0, j1, states, bounds, s_idx_x, s_idx_y, cx, cy, ce, schedule[-1] / 2.0)
    n_eval += ne
    arg2 = {(states.names[s], states.names[sp]): (float(cx[s, sp]), float(cy[s, sp])) for s in range(S) for sp in range(S)}
    e2 = ce
    return E1, E2, e1, arg1, e2, arg2, n_eval


@dataclass
class FitReport:
    grid: np.ndarray
    scan1: dict  # bond type -> (n_grid,) energies relative to all-trans
    scan2: dict  # bond type -> (n_grid, n_grid) energies relative to all-trans
    n_evaluations: int
    model: RISModel
    argmin1: dict  # bond type -> {state: angle at basin minimum}
    argmin2: dict  # bond type -> {(s, s'): (angle, angle)}


def fit_ris(
    polymer: Polymer,
    calc: Calculator,
    states: RISStates = THREE_STATE,
    step: float = 10.0,
    n_monomers: int = 5,
    name: str | None = None,
    symmetrize: bool = True,
    adapt_angles: bool = True,
    third_order: bool = False,
    cap: float = 50.0,
    scan: str = "adaptive",
    coarse_step: float = 30.0,
) -> FitReport:
    """Derive an RIS model from dihedral scans of a short oligomer.

    First-order energies come from a 1D scan of one bond of each type (all other
    bonds trans); pair energies from a 2D scan of two consecutive bonds, minus the
    first-order terms.  Energies are basin minima, relative to all-trans.  With
    ``symmetrize`` the model is averaged with its mirror image (G+ <-> G-), which is
    exact for achiral chains and removes grid/refinement noise.  With ``adapt_angles``
    the state dihedral angles are moved to the 1D-scan basin minima (averaged over
    bond types and mirror pairs), as in classical RIS parametrisations.  With
    ``third_order`` triplet corrections are added (see :func:`_fit_third_order`),
    which is what distinguishes e.g. TG+TG+ (3/1 helix) from TG+TG- (alpha-PVDF).
    Energies above ``cap`` (steric overlap) are clipped to ``cap``.

    Every conformer of every scan is built in one batched :func:`~polyfind.chain.build_chain_batch`
    call and evaluated in one ``calc.energy_batch`` call per scan (or refinement round), so an
    expensive (MLIP/DFT) calculator sees few, large batches rather than many single points.

    ``scan`` selects the search strategy:

    * ``"dense"``: the original ``step``-degree grid over the full 360 deg, basin minima
      read off the grid -- exact and reproducible, but most points are far from any basin.
    * ``"adaptive"`` (default): a coarse (``coarse_step``, default 30 deg) grid locates each
      basin, then a bounded quadratic-interpolation refinement (not restricted to a grid, so
      it can find basin minima at least as good as the dense scan's, usually better) polishes
      each basin's minimum at ``step`` resolution. Several times fewer evaluations than the
      dense scan for comparable or better accuracy.
    """
    B = polymer.bonds_per_repeat
    N = max(n_monomers * B, 2 * B + 6)
    base = np.full(N, 180.0)
    ref = calc.energy(build_chain(polymer, base))
    n_eval = 1
    e1 = np.zeros((B, states.n))
    e2 = np.zeros((B, states.n, states.n))
    scan1, scan2, arg1, arg2 = {}, {}, {}, {}
    if scan == "dense":
        report_grid = np.arange(-180.0, 180.0, step)
        basins = _basins(states, report_grid)
        for b in range(B):
            # scanned bond j0 of type b, nearest the middle of the oligomer
            j0 = (N // 2 - 1) - ((N // 2 - 1) - b) % B
            E1, E2, e1_b, arg1_b, e2_b, arg2_b, ne = _dense_scan_bond(polymer, calc, states, base, j0, report_grid, basins, ref)
            scan1[b], scan2[b] = E1, E2
            e1[b], arg1[b] = e1_b, arg1_b
            e2[b], arg2[b] = e2_b, arg2_b
            n_eval += ne
    elif scan == "adaptive":
        # offset by half a step so grid points never land exactly on a basin boundary
        # (they would if coarse_step divides the state spacing evenly, e.g. 30 deg into
        # a 3-state model's 120 deg spacing); a boundary tie excludes the outer edge of
        # every basin from the coarse sample, which is exactly where a steric-clash
        # basin's true minimum can sit, right next to the repulsive wall.
        report_grid = np.arange(-180.0 + coarse_step / 2.0, 180.0, coarse_step)
        basins_coarse = _basins(states, report_grid)
        bounds = _basin_bounds(states)
        for b in range(B):
            j0 = (N // 2 - 1) - ((N // 2 - 1) - b) % B
            E1, E2, e1_b, arg1_b, e2_b, arg2_b, ne = _adaptive_scan_bond(
                polymer, calc, states, base, j0, report_grid, basins_coarse, bounds, ref, step, coarse_step
            )
            scan1[b], scan2[b] = E1, E2
            e1[b], arg1[b] = e1_b, arg1_b
            e2[b], arg2[b] = e2_b, arg2_b
            n_eval += ne
    else:
        raise ValueError(f"unknown scan mode {scan!r}; expected 'dense' or 'adaptive'")
    # remove first-order contributions from the pair scan (the 1D scans were done with the
    # partner trans, so e1 already includes the state's interaction with a trans partner)
    for b in range(B):
        bn = (b + 1) % B
        e2[b] = e2[b] - e1[b][:, None] - e1[bn][None, :]
    if symmetrize:
        m = np.array(states.mirror)
        e1 = 0.5 * (e1 + e1[:, m])
        e2 = 0.5 * (e2 + e2[:, m][:, :, m])
    e1 = np.minimum(e1, cap)
    e2 = np.minimum(e2, cap)
    if adapt_angles:
        angs = []
        for s, nm in enumerate(states.names):
            vals = np.array([arg1[b][nm] for b in range(B)])
            # circular mean of the basin minima across bond types
            ang = np.degrees(np.arctan2(np.sin(np.deg2rad(vals)).mean(), np.cos(np.deg2rad(vals)).mean()))
            angs.append(ang)
        angs = np.array(angs)
        if symmetrize:  # mirror pairs get equal magnitude and opposite sign
            m = np.array(states.mirror)
            self_mirror = m == np.arange(states.n)
            angs = np.where(self_mirror, angs, np.sign(angs) * 0.5 * (np.abs(angs) + np.abs(angs[m])))
        # keep trans at exactly 180 (a basin minimum at -180 == 180)
        angs = np.where(np.abs(np.abs(angs) - 180.0) < 1e-6, 180.0, angs)
        states = RISStates(states.names, tuple(float(a) for a in angs), states.mirror)
    e3 = None
    if third_order:
        e3 = np.clip(_fit_third_order(polymer, calc, states, N, base), -cap, cap)
        n_eval += B * states.n ** 3 * 4
        if symmetrize:
            m = np.array(states.mirror)
            e3 = 0.5 * (e3 + e3[:, m][:, :, m][:, :, :, m])
    model = RISModel(states, B, e1, e2, e3, name=name or f"{polymer.name}-fit")
    return FitReport(report_grid, scan1, scan2, n_eval, model, arg1, arg2)


def _fit_third_order(polymer: Polymer, calc: Calculator, states: RISStates, N: int, base: np.ndarray) -> np.ndarray:
    """Triplet corrections by inclusion-exclusion at the state angles:

        e3[b, a, c, d] = E(a c d) - E(a c T) - E(T c d) + E(T c T)

    i.e. the part of the energy of three consecutive states not captured by the
    pairs (dominated by 1,6-type contacts).  4 * S^3 evaluations per bond type, built
    and evaluated as one batch per bond type.
    """
    B, S = polymer.bonds_per_repeat, states.n
    t_idx = states.index("T") if "T" in states.names else 0
    e3 = np.zeros((B, S, S, S))
    for b in range(B):
        j0 = (N // 2 - 1) - ((N // 2 - 1) - b) % B
        rows, keys = [], []
        for a in range(S):
            for c in range(S):
                for d in range(S):
                    for (x, y, z) in ((a, c, d), (a, c, t_idx), (t_idx, c, d), (t_idx, c, t_idx)):
                        dd = base.copy()
                        dd[j0], dd[j0 + 1], dd[j0 + 2] = states.angles[x], states.angles[y], states.angles[z]
                        rows.append(dd)
                    keys.append((a, c, d))
        dihs = np.array(rows)
        template, coords = build_chain_batch(polymer, dihs)
        E = calc.energy_batch((template, coords)).reshape(-1, 4)
        for (a, c, d), row in zip(keys, E):
            e3[b, a, c, d] = row[0] - row[1] - row[2] + row[3]
    return e3


def oligomer_energy_of_sequence(polymer: Polymer, calc: Calculator, seq, states: RISStates = THREE_STATE, n_periods: int = 6) -> float:
    """Energy per monomer of a finite oligomer built from a periodic sequence (sanity check)."""
    seq = list(seq)
    dih = np.array([states.angles[s] for s in seq * n_periods])
    return calc.energy(build_chain(polymer, dih)) / (len(dih) / polymer.bonds_per_repeat)
