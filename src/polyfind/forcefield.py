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
from .chain import Structure, build_chain
from .polymers import Polymer, RISStates, THREE_STATE, lj_params
from .ris import RISModel

COULOMB = 332.0637  # kcal A / (mol e^2)


class Calculator(Protocol):
    def energy(self, struct: Structure) -> float: ...

    def energy_batch(self, structs: Sequence[Structure]) -> np.ndarray: ...


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

    def energy_batch(self, structs: Sequence[Structure]) -> np.ndarray:
        """Vectorised over conformers (same topology); runs on GPU if available."""
        if not structs:
            return np.empty(0)
        xp = bk.get_backend()
        top = self._topology(structs[0])
        coords = xp.asarray(np.stack([s.coords for s in structs]), dtype=bk.float_dtype())
        out = np.empty(len(structs))
        chunk = 2048
        for k in range(0, len(structs), chunk):
            out[k : k + chunk] = bk.to_numpy(self._energy_from_coords(top, coords[k : k + chunk], xp))
        return out

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

    def _atoms(self, struct: Structure):
        from ase import Atoms

        atoms = Atoms(symbols=struct.elements, positions=struct.coords)
        atoms.calc = self.calc
        return atoms

    def energy(self, struct: Structure) -> float:
        return float(self._atoms(struct).get_potential_energy()) * self.EV_TO_KCAL

    def energy_batch(self, structs: Sequence[Structure]) -> np.ndarray:
        return np.array([self.energy(s) for s in structs])


# ---------------------------------------------------------------------- fitting
def _basins(states: RISStates, grid: np.ndarray) -> list[np.ndarray]:
    """Grid indices belonging to each state's basin (nearest ideal angle on the circle)."""
    diff = np.abs((grid[:, None] - np.array(states.angles)[None, :] + 180.0) % 360.0 - 180.0)
    owner = diff.argmin(axis=1)
    return [np.where(owner == s)[0] for s in range(states.n)]


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
) -> FitReport:
    """Derive an RIS model from dihedral scans of a short oligomer.

    First-order energies come from a 1D scan of one bond of each type (all other
    bonds trans); pair energies from a 2D scan of two consecutive bonds, minus the
    first-order terms.  Energies are basin minima on the grid, relative to all-trans.
    With ``symmetrize`` the model is averaged with its mirror image (G+ <-> G-),
    which is exact for achiral chains and removes grid noise.  With ``adapt_angles``
    the state dihedral angles are moved to the 1D-scan basin minima (averaged over
    bond types and mirror pairs), as in classical RIS parametrisations.  With
    ``third_order`` triplet corrections are added (see :func:`_fit_third_order`),
    which is what distinguishes e.g. TG+TG+ (3/1 helix) from TG+TG- (alpha-PVDF).
    Energies above ``cap`` (steric overlap) are clipped to ``cap``.
    """
    B = polymer.bonds_per_repeat
    N = max(n_monomers * B, 2 * B + 6)
    grid = np.arange(-180.0, 180.0, step)
    basins = _basins(states, grid)
    base = np.full(N, 180.0)
    ref = calc.energy(build_chain(polymer, base))
    n_eval = 1
    e1 = np.zeros((B, states.n))
    e2 = np.zeros((B, states.n, states.n))
    scan1, scan2, arg1, arg2 = {}, {}, {}, {}
    for b in range(B):
        # scanned bond j0 of type b, nearest the middle of the oligomer
        j0 = (N // 2 - 1) - ((N // 2 - 1) - b) % B
        # ---- first order
        structs = []
        for phi in grid:
            d = base.copy()
            d[j0] = phi
            structs.append(build_chain(polymer, d))
        E = calc.energy_batch(structs) - ref
        n_eval += len(structs)
        scan1[b] = E
        arg1[b] = {}
        for s, idx in enumerate(basins):
            k = idx[np.argmin(E[idx])]
            e1[b, s] = E[k]
            arg1[b][states.names[s]] = float(grid[k])
        # ---- second order: bonds j0 (type b) and j0+1 (type b+1)
        structs = []
        for phi in grid:
            for psi in grid:
                d = base.copy()
                d[j0], d[j0 + 1] = phi, psi
                structs.append(build_chain(polymer, d))
        E2 = (calc.energy_batch(structs) - ref).reshape(len(grid), len(grid))
        n_eval += len(structs)
        scan2[b] = E2
        arg2[b] = {}
        for s, idx in enumerate(basins):
            for sp, idxp in enumerate(basins):
                sub = E2[np.ix_(idx, idxp)]
                k = np.unravel_index(np.argmin(sub), sub.shape)
                arg2[b][(states.names[s], states.names[sp])] = (float(grid[idx[k[0]]]), float(grid[idxp[k[1]]]))
                e2[b, s, sp] = sub[k]
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
    return FitReport(grid, scan1, scan2, n_eval, model, arg1, arg2)


def _fit_third_order(polymer: Polymer, calc: Calculator, states: RISStates, N: int, base: np.ndarray) -> np.ndarray:
    """Triplet corrections by inclusion-exclusion at the state angles:

        e3[b, a, c, d] = E(a c d) - E(a c T) - E(T c d) + E(T c T)

    i.e. the part of the energy of three consecutive states not captured by the
    pairs (dominated by 1,6-type contacts).  4 * S^3 evaluations per bond type.
    """
    B, S = polymer.bonds_per_repeat, states.n
    t_idx = states.index("T") if "T" in states.names else 0
    e3 = np.zeros((B, S, S, S))
    for b in range(B):
        j0 = (N // 2 - 1) - ((N // 2 - 1) - b) % B
        structs, keys = [], []
        for a in range(S):
            for c in range(S):
                for d in range(S):
                    for (x, y, z) in ((a, c, d), (a, c, t_idx), (t_idx, c, d), (t_idx, c, t_idx)):
                        dd = base.copy()
                        dd[j0], dd[j0 + 1], dd[j0 + 2] = states.angles[x], states.angles[y], states.angles[z]
                        structs.append(build_chain(polymer, dd))
                    keys.append((a, c, d))
        E = calc.energy_batch(structs).reshape(-1, 4)
        for (a, c, d), row in zip(keys, E):
            e3[b, a, c, d] = row[0] - row[1] - row[2] + row[3]
    return e3


def oligomer_energy_of_sequence(polymer: Polymer, calc: Calculator, seq, states: RISStates = THREE_STATE, n_periods: int = 6) -> float:
    """Energy per monomer of a finite oligomer built from a periodic sequence (sanity check)."""
    seq = list(seq)
    dih = np.array([states.angles[s] for s in seq * n_periods])
    return calc.energy(build_chain(polymer, dih)) / (len(dih) / polymer.bonds_per_repeat)
