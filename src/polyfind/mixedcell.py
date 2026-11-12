r"""Cells of several chains of ONE sequence, each with its own conformation and placement.

:class:`polyfind.pack.CrystalPacker` places one chain, or two copies of one chain, per cell:
every chain of the crystal is the same rigid repeat.  That is the right object for a polymorph
and the wrong one for a defect.  A kink pair nucleated in one chain of an all-trans crystal has
neighbours that are still all-trans; modelling it with the packer's cell kinks every chain in
register, which is a different crystal (docs/CANDIDATE_MATERIALS.md section 8).  This module is
the cell in which the chains differ:

* :func:`placed_energy_and_grad` is :func:`polyfind.phonon.cell_energy_and_grad` for any number
  of chains, each with its own placing frame, plus the derivative with respect to the lattice
  vectors.  It assembles the packer's own energy term for term -- the same-site column of each
  chain with the bonded exclusions, the pair kernel between every pair of chains over every
  lattice image in reach, the Ewald sum with its per-chain exclusion correction and the induced
  dipoles -- from the tables of a packer built with ``n_chains`` equal to the number of chains
  in the cell.  The chains share one sequence, hence one set of elements, charges and bonded
  topology; only their coordinates differ.  At the packer's own two-chain placement it is the
  packer's energy to rounding, and a tiling of that cell gives the same energy per chain
  (``tests/test_mixedcell.py``).
* :class:`MixedModel` puts the chain builder's internal coordinates on top: each chain is built
  from its own torsions and backbone angles (:func:`polyfind.linegroup.repeat_chains_and_jacobian`,
  the builder and Jacobian the crystal refinement uses), placed by a fractional lateral position,
  a setting angle and an axial shift, and the energy carries what the refinement's does: the
  lattice energy, the Fourier torsion term and the backbone-angle restraint of every chain.  All
  derivatives are analytic.

**One repeat for all chains, held exactly.**  The cell has one ``c``, and a chain built from its
internal coordinates has its own natural repeat ``c_s`` and, off a commensurate conformation, a
residual rotation ``w_s`` of its repeat transform.  Periodicity is the four constraints
``w_s = 0``, ``c_s = c`` per chain (:meth:`MixedModel.closure`, with their analytic Jacobian), and
:meth:`MixedModel.relax` imposes them as equality constraints of a sequential quadratic programme;
:func:`neb` keeps every image on the constraint surface by projection.  A penalty is not enough,
and it is also the wrong tool: with every torsion free a chain can buy lattice energy with a tenth
of a degree of residual rotation (a few tenths of a kcal/mol per chain at the refinement's default
penalty, the size of the differences this module is for), and a penalty stiff enough to stop that
makes the torsion space too ill-conditioned for a quasi-Newton method to finish.  The penalty and
multiplier terms remain available (``pen_rot``, ``k_c``, off by default) for a caller who wants
the refinement's own objective.  A chain shorter than its neighbours, which is what a kink pair
makes, reaches the cell's repeat through its own bends and torsions and carries the tension that
costs -- the periodic model's finite-length effect, which falls as the inverse repeat length.

**What is not here.**  Charge flux and an applied field (refused), and the pendant tilts, which
:func:`tilt_relaxation` adds afterwards exactly as :func:`polyfind.pendants.tilt_relaxation` does
for the packer's cell.
"""
from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np
from scipy.optimize import minimize

from . import backend as bk
from . import phonon as PH
from .linegroup import repeat_chains_and_jacobian
from .pack import PeriodicChain

__all__ = [
    "MixedModel",
    "MixedState",
    "lattice_rows",
    "neb",
    "placed_energy_and_grad",
    "placing_frame",
    "supercell_state",
    "tile_state",
    "tilt_relaxation",
]

_FLIP = np.diag([1.0, -1.0, -1.0])
_RAD = np.pi / 180.0
CLOSURE_C = 10.0  # the repeat-mismatch constraint in tenths of an angstrom, the size of a degree of residual rotation


def placing_frame(phi_deg: float, flip: bool) -> np.ndarray:
    """``R`` (3, 3) with ``P = X @ R.T + offset``: ``Rz(phi)``, after the ``(x, -y, -z)`` flip if set.

    The linear part of :meth:`polyfind.pack.CrystalPacker._place` for one chain."""
    c, s = np.cos(phi_deg * _RAD), np.sin(phi_deg * _RAD)
    R = np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
    return R @ _FLIP if flip else R


def lattice_rows(A: float, B: float, gamma_deg: float, c: float) -> np.ndarray:
    """Lattice vectors as rows: ``A`` along x, ``B`` in the plane at ``gamma``, ``c`` along z."""
    g = gamma_deg * _RAD
    return np.array([[A, 0.0, 0.0], [B * np.cos(g), B * np.sin(g), 0.0], [0.0, 0.0, c]])


# ----------------------------------------------------------------------- the energy
def placed_energy_and_grad(packer, frames, Pn, latn, valences=None):
    """``(E, dE/dPn (N, 3), dE/dlattice (3, 3))`` of ``n_chains`` placed chains of one sequence.

    ``packer`` is a :class:`~polyfind.pack.CrystalPacker` built with ``n_chains`` equal to the
    number of chains in this cell (its pair, charge and polarizability tables are the repeat's,
    tiled); ``frames`` (n_chains, 3, 3) are the placing frames of :func:`placing_frame`, which
    take each chain's block of ``Pn`` back to that chain's own frame, where its same-site column
    and its exclusion correction are evaluated; ``latn`` has ``c`` along z and the other two
    vectors in the plane.  ``valences`` is ``None`` (the packer's own valence terms on every
    chain, if it has any) or one :class:`~polyfind.pack.ChainValence` (or ``None``) per chain.

    ``E`` is the packer's lattice energy in kcal/mol per cell **without** the Fourier torsion
    term, which is a function of the dihedrals and is the caller's to add per chain.  The
    lattice derivative is taken at fixed Cartesian coordinates, the packer's own convention; its
    ``[2, 2]`` entry is ``dE/dc`` and carries every dependence on the repeat (the image shifts
    of the same-site columns and of the exclusions included).
    """
    n, N, nc = packer.n, packer.N, packer.n_chains
    Pn = np.asarray(Pn, dtype=float).reshape(N, 3)
    latn = np.asarray(latn, dtype=float).reshape(3, 3)
    frames = np.asarray(frames, dtype=float).reshape(nc, 3, 3)
    if getattr(packer, "_flux", None) is not None:
        raise ValueError("charge flux is not carried by the mixed cell: every chain has the repeat's fixed charges")
    if getattr(packer, "_field_on", False):
        raise ValueError("an applied field is not carried by the mixed cell; clear it with set_field(None)")
    if abs(latn[0, 2]) + abs(latn[1, 2]) + abs(latn[2, 0]) + abs(latn[2, 1]) > 1e-12:
        raise ValueError("the mixed cell wants c along z and the other two lattice vectors in the plane")
    cz = float(latn[2, 2])
    rc, rc2 = float(packer.rc), float(packer.rc2)
    K = int(np.ceil(rc / cz)) + 1
    flips = frames[:, 2, 2] < 0.0
    sl = [slice(s * n, (s + 1) * n) for s in range(nc)]
    Xs = [Pn[sl[s]] @ frames[s] for s in range(nc)]  # each chain's own frame, up to a translation
    q = np.asarray(packer._q_cell[:n], dtype=float)
    q_cell = np.asarray(packer._q_cell, dtype=float)
    xp, dt = packer.xp, packer._dt
    E = 0.0
    gP = np.zeros((N, 3))
    glat = np.zeros((3, 3))
    gX = [np.zeros((n, 3)) for _ in range(nc)]

    # --- per chain, in its own frame: same-site column and valence ------------------------------
    for s in range(nc):
        intra, gXi, gci, _ = packer._intra_column_and_grad(xp.asarray(Xs[s][None], dtype=dt), np.array([cz]), K)
        E += 0.5 * float(bk.to_numpy(intra)[0])
        gX[s] += 0.5 * gXi
        glat[2, 2] += 0.5 * gci
        val = packer._valence if valences is None else valences[s]
        if val is not None:
            Ev, gXv, gcv = val.energy_and_grad(Xs[s], cz, n_atoms=n)
            E += float(Ev)
            gX[s] += gXv
            glat[2, 2] += float(gcv)

    # --- every pair of chains over every lattice image in reach, on the placed coordinates ------
    # E = 1/2 sum over ordered (s, t, image) except a chain with itself at no lateral shift, and
    # (s, t, image) equals (t, s, -image): so each unordered block is evaluated once at weight 1.
    A_nn = np.asarray(bk.to_numpy(packer._A_nn), dtype=float)
    B_nn = np.asarray(bk.to_numpy(packer._B_nn), dtype=float)
    qq, qqf, cst = PH._charge_tables(packer, q, q)
    A2, B2 = latn[0, :2], latn[1, :2]
    area = abs(A2[0] * B2[1] - A2[1] * B2[0])
    wa, wb = area / float(np.linalg.norm(B2)), area / float(np.linalg.norm(A2))
    Minv = np.linalg.inv(np.array([A2, B2]))  # lateral Cartesian -> fractional (row vectors)
    cen = np.array([Pn[sl[s]][:, :2].mean(axis=0) for s in range(nc)])
    rad = np.array([float(np.linalg.norm(Pn[sl[s]][:, :2] - cen[s], axis=1).max()) for s in range(nc)])
    for s in range(nc):
        Ps = Pn[sl[s]]
        for t in range(s, nc):
            Pt = Pn[sl[t]]
            reach = rc + rad[s] + rad[t]
            d0 = cen[s] - cen[t]
            f0 = d0 @ Minv
            ri, rj = int(np.ceil(reach / wa)) + 1, int(np.ceil(reach / wb)) + 1
            i0, j0 = int(np.round(f0[0])), int(np.round(f0[1]))
            dz0 = Ps[:, None, 2] - Pt[None, :, 2]
            ks = np.arange(int(np.floor((float(dz0.min()) - rc) / cz)), int(np.ceil((float(dz0.max()) + rc) / cz)) + 1, dtype=float)
            dz = dz0[None, :, :] - ks[:, None, None] * cz
            dz2 = dz * dz
            for i in range(i0 - ri, i0 + ri + 1):
                for j in range(j0 - rj, j0 + rj + 1):
                    if s == t and (i < 0 or (i == 0 and j <= 0)):
                        continue  # the chain itself, and the half of its images counted by the other half
                    sh = i * A2 + j * B2
                    if float(np.linalg.norm(d0 - sh)) >= reach:
                        continue
                    dx = Ps[:, None, 0] - Pt[None, :, 0] - sh[0]
                    dy = Ps[:, None, 1] - Pt[None, :, 1] - sh[1]
                    dxy2 = dx * dx + dy * dy
                    if float(dxy2.min()) >= rc2:
                        continue
                    v, dv = PH._pair_and_dv(packer, dxy2[None, :, :] + dz2, A_nn, B_nn, qq, qqf, cst)
                    E += float(v.sum())
                    t2 = 2.0 * dv  # dV/d(separation vector) = 2 dV/dr^2 * (separation vector)
                    tk = t2.sum(axis=0)
                    fx, fy = tk * dx, tk * dy
                    fzk = t2 * dz
                    fz = fzk.sum(axis=0)
                    for d, f in enumerate((fx, fy, fz)):
                        gP[sl[s], d] += f.sum(axis=1)
                        gP[sl[t], d] -= f.sum(axis=0)
                    # the image shift is i A + j B + k C, and dE/d(shift) = -sum of the pair forces
                    G = -np.array([float(fx.sum()), float(fy.sum()), float(fz.sum())])
                    glat[0] += i * G
                    glat[1] += j * G
                    tkk = (t2 * ks[:, None, None]).sum(axis=0)
                    glat[2, 0] -= float((tkk * dx).sum())
                    glat[2, 1] -= float((tkk * dy).sum())
                    glat[2, 2] -= float((fzk * ks[:, None, None]).sum())

    # --- electrostatics of an Ewald packer --------------------------------------------------------
    if packer._ewald is not None:
        if packer.polarizable is not None:
            e_es, _, _, (gPe, glate, _, gce) = packer._polarize(Pn, q_cell, latn, flips.astype(float), grad=True)
            glat[2, 2] += float(gce)
        else:
            idx = packer._es_idx
            if idx is None:
                terms = packer._ewald.terms(Pn, q_cell, latn, grad=True)
                gPe = terms.grad_coords
            else:
                terms = packer._ewald.terms(Pn[idx], q_cell[idx], latn, grad=True)
                gPe = np.zeros((N, 3))
                gPe[idx] = terms.grad_coords
            e_es, glate = terms.total, terms.grad_lattice
        E += float(e_es)
        gP += gPe
        glat += glate
        for s in range(nc):
            ex_e, ex_gX, ex_gc, _ = packer._exclusion(Xs[s], q, cz, K, grad=True)
            E += float(ex_e)
            gX[s] += ex_gX
            glat[2, 2] += float(ex_gc)

    for s in range(nc):
        gP[sl[s]] += gX[s] @ frames[s].T
    return E, gP, glat


# ------------------------------------------------------------------------- the model
@dataclass
class MixedState:
    """The variables of a mixed cell, in the units the optimiser sees.

    ``lattice`` is ``(A, B, gamma, c)`` in A and degrees; ``place`` (n_chains, 4) is per chain the
    fractional lateral position of its axis ``(fa, fb)``, its setting angle in degrees and its
    axial shift in A; ``tors`` (n_chains, P) and ``angs`` (n_chains, R) are the torsions and
    backbone angles of each chain's repeat in degrees."""

    lattice: np.ndarray
    place: np.ndarray
    tors: np.ndarray
    angs: np.ndarray

    def vector(self) -> np.ndarray:
        nc = self.place.shape[0]
        return np.concatenate([np.asarray(self.lattice, dtype=float).ravel(), np.asarray(self.place, dtype=float).ravel(),
                               np.concatenate([np.asarray(self.tors, dtype=float), np.asarray(self.angs, dtype=float)], axis=1).reshape(nc, -1).ravel()])

    def copy(self) -> "MixedState":
        return MixedState(np.array(self.lattice, dtype=float), np.array(self.place, dtype=float),
                          np.array(self.tors, dtype=float), np.array(self.angs, dtype=float))

    @staticmethod
    def from_vector(x, n_chains: int, n_torsions: int) -> "MixedState":
        """The state a variable vector of ``n_chains`` chains of ``n_torsions`` torsions each
        stands for (the inverse of :meth:`vector`)."""
        x = np.asarray(x, dtype=float).ravel()
        nc, P = int(n_chains), int(n_torsions)
        internal = x[4 + 4 * nc:].reshape(nc, -1)
        return MixedState(x[:4].copy(), x[4:4 + 4 * nc].reshape(nc, 4).copy(), internal[:, :P].copy(), internal[:, P:].copy())


def supercell_state(params, c: float, tors, angs, nx: int = 2, ny: int = 1):
    """``(state, flips)`` of the ``nx x ny`` tiling of the packer's two-chain cell at ``params``.

    Chain order is ``(ix, iy, k)`` with ``k`` the chain within the packer's cell: ``k = 0`` at the
    cell origin with setting angle ``phi1``, ``k = 1`` at ``(a + b) / 2``, shifted by ``dz`` and
    flipped if the cell is.  Every chain starts with the same ``tors`` (P,) and ``angs`` (R,)."""
    a, b, gam, phi1, phi2, dz, flip = (float(v) for v in np.asarray(params, dtype=float).reshape(7))
    place, flips = [], []
    for ix in range(nx):
        for iy in range(ny):
            for k in (0, 1):
                place.append([(ix + 0.5 * k) / nx, (iy + 0.5 * k) / ny, phi2 if k else phi1, dz if k else 0.0])
                flips.append(bool(k and flip > 0.5))
    nc = len(place)
    tors, angs = np.asarray(tors, dtype=float).ravel(), np.asarray(angs, dtype=float).ravel()
    return MixedState(np.array([nx * a, ny * b, gam, float(c)]), np.array(place),
                      np.tile(tors, (nc, 1)), np.tile(angs, (nc, 1))), np.array(flips)


def tile_state(state: MixedState, flips, nx: int = 2, ny: int = 1):
    """``(state, flips)`` of the ``nx x ny`` tiling of a mixed cell: every chain of ``state``
    repeated at each lattice point of the tiling with its own conformation, setting angle and
    axial shift.  Chain order is ``(ix, iy, k)`` with ``k`` the chain within ``state``, which for a
    two-chain cell is the order of :func:`supercell_state`."""
    A, B, gam, c = (float(v) for v in state.lattice)
    flips = np.asarray(flips, dtype=bool).ravel()
    place, tors, angs, fl = [], [], [], []
    for ix in range(nx):
        for iy in range(ny):
            for k in range(state.place.shape[0]):
                fa, fb, phi, z = state.place[k]
                place.append([(ix + fa) / nx, (iy + fb) / ny, phi, z])
                tors.append(state.tors[k])
                angs.append(state.angs[k])
                fl.append(bool(flips[k]))
    return MixedState(np.array([nx * A, ny * B, gam, c]), np.array(place), np.array(tors, dtype=float), np.array(angs, dtype=float)), np.array(fl)


class MixedModel:
    """Energy and analytic gradient of a mixed cell in the chain builder's internal coordinates.

    ``packer`` is built with ``n_chains`` equal to the number of chains here, from any chain of
    the sequence (the tables depend on the topology alone); ``refs`` is one
    :class:`~polyfind.pack.PeriodicChain` per chain, the frame each chain's coordinates are
    aligned to (the builder's own frame is a principal axis of the cross section and flips
    discontinuously; see :func:`polyfind.pack._align_about_z`), fixed for the life of the model so
    a setting angle means the same thing at every point of a path.  ``angle_stiffness`` is the
    refinement's backbone-angle restraint (kcal/mol/rad^2,
    :func:`polyfind.refine.resolve_angle_stiffness`).  ``pen_rot`` (the refinement's
    commensurability penalty, per degree^2) and ``k_c`` (a penalty tying each chain's natural repeat
    to the cell's, kcal/mol/A^2) are zero by default: :meth:`relax` holds the closure exactly.
    """

    def __init__(self, polymer, name: str, packer, flips, refs, angle_stiffness: float,
                 pen_rot: float = 0.0, k_c: float = 0.0, valences=None):
        self.polymer, self.name, self.pk = polymer, name, packer
        self.nc, self.n = int(packer.n_chains), int(packer.n)
        self.flips = np.asarray(flips, dtype=bool).reshape(self.nc)
        if len(refs) != self.nc:
            raise ValueError(f"{len(refs)} reference chains for a cell of {self.nc}")
        self.refs = [np.asarray(r.coords if isinstance(r, PeriodicChain) else r, dtype=float) for r in refs]
        first = next((r for r in refs if isinstance(r, PeriodicChain)), None)
        if first is None:
            raise ValueError("at least one reference must be a PeriodicChain (it carries the repeat's torsion count)")
        self.P = int(len(first.dihedrals))
        self.angles0 = np.array([polymer.backbone[k].backbone_angle for k in range(polymer.bonds_per_repeat)], dtype=float)
        self.R = int(self.angles0.size)
        self.n_repeats = self.P / self.R
        self.n_monomers = int(first.n_monomers)
        self.angle_stiffness = float(angle_stiffness)
        self.pen_rot, self.k_c = float(pen_rot), float(k_c)
        self.valences = valences
        self.n_evaluations = 0
        self.reset_multipliers()

    # --- closure multipliers -----------------------------------------------------------------
    def reset_multipliers(self) -> None:
        """Zero the augmented-Lagrangian multipliers of the closure constraints."""
        self.lam_w = np.zeros((self.nc, 3))  # kcal/mol per degree, on each chain's rotation vector
        self.lam_c = np.zeros(self.nc)       # kcal/mol per A, on each chain's repeat mismatch

    def multipliers(self) -> tuple:
        return self.lam_w.copy(), self.lam_c.copy()

    def set_multipliers(self, lam) -> None:
        self.lam_w, self.lam_c = np.array(lam[0], dtype=float).reshape(self.nc, 3), np.array(lam[1], dtype=float).reshape(self.nc)

    def update_multipliers(self, x) -> tuple:
        """The first-order multiplier update at ``x``; returns the residuals it used,
        ``(max |w| in degrees, max |c_s - c| in A)``."""
        c = float(np.asarray(x, dtype=float).ravel()[3])
        built = self.build(x)[4]
        w = np.array([b[3] for b in built])
        dc = np.array([float(b[0].c) - c for b in built])
        self.lam_w = self.lam_w + 2.0 * self.pen_rot * w
        self.lam_c = self.lam_c + 2.0 * self.k_c * dc
        return float(np.linalg.norm(w, axis=1).max()), float(np.abs(dc).max())

    # --- the variable vector -------------------------------------------------------------
    @property
    def n_var(self) -> int:
        return 4 + self.nc * (4 + self.P + self.R)

    def state(self, x) -> MixedState:
        x = np.asarray(x, dtype=float).ravel()
        nc, P, R = self.nc, self.P, self.R
        internal = x[4 + 4 * nc:].reshape(nc, P + R)
        return MixedState(x[:4].copy(), x[4:4 + 4 * nc].reshape(nc, 4).copy(), internal[:, :P].copy(), internal[:, P:].copy())

    def index(self, what: str, chains=None) -> np.ndarray:
        """Indices into the variable vector: ``what`` is ``"lattice"`` (A, B, gamma), ``"c"``,
        ``"place"`` (fa, fb, phi, z of ``chains``), ``"xy"``, ``"phi"``, ``"z"``, ``"tors"``,
        ``"angs"`` or ``"internal"``; ``chains`` defaults to all."""
        nc, P, R = self.nc, self.P, self.R
        ch = range(nc) if chains is None else [int(s) for s in np.atleast_1d(chains)]
        if what == "lattice":
            return np.array([0, 1, 2])
        if what == "c":
            return np.array([3])
        cols = {"place": [0, 1, 2, 3], "xy": [0, 1], "phi": [2], "z": [3]}
        if what in cols:
            return np.array([4 + 4 * s + k for s in ch for k in cols[what]], dtype=int)
        base = 4 + 4 * nc
        span = {"tors": range(0, P), "angs": range(P, P + R), "internal": range(0, P + R)}[what]
        return np.array([base + s * (P + R) + k for s in ch for k in span], dtype=int)

    def scale(self, x) -> np.ndarray:
        """A length per variable that makes one unit of each comparable (about 1 A of atom motion):
        the preconditioner of :meth:`relax` and the metric of :func:`neb`."""
        st = self.state(x)
        sc = np.ones(self.n_var)
        sc[2] = 2.0
        sc[self.index("xy")] = np.tile([1.0 / st.lattice[0], 1.0 / st.lattice[1]], self.nc)
        sc[self.index("phi")] = 10.0
        sc[self.index("tors")] = 10.0
        sc[self.index("angs")] = 2.0
        return sc

    # --- geometry --------------------------------------------------------------------------
    def build(self, x):
        """``(latn, frames, offsets, Pn, built)``: the placed cell and each chain's builder output
        ``(chain, dcoords (P + R, n, 3), dc (P + R,), w (3,), dw (P + R, 3))``."""
        st = self.state(x)
        A, B, gam, c = st.lattice
        latn = lattice_rows(A, B, gam, c)
        nc, n = self.nc, self.n
        Pn = np.empty((nc * n, 3))
        frames = np.empty((nc, 3, 3))
        offs = np.empty((nc, 3))
        built = []
        for s in range(nc):
            chains, dX, dc, w, dw = repeat_chains_and_jacobian(self.polymer, self.name, st.tors[s][None], st.angs[s], self.refs[s])
            built.append((chains[0], dX[0], dc[0], w[0], dw[0]))
            fa, fb, phi, z = st.place[s]
            frames[s] = placing_frame(phi, bool(self.flips[s]))
            offs[s] = fa * latn[0] + fb * latn[1]
            offs[s, 2] += z
            Pn[s * n:(s + 1) * n] = chains[0].coords @ frames[s].T + offs[s]
        return latn, frames, offs, Pn, built

    def chains(self, x) -> list:
        """The :class:`~polyfind.pack.PeriodicChain` of every chain at ``x`` (own frame)."""
        return [b[0] for b in self.build(x)[4]]

    # --- energy ------------------------------------------------------------------------------
    def value_and_grad(self, x, parts: dict | None = None):
        """Energy with the closure penalties, and its gradient with respect to every variable.

        ``parts``, if a dict is given, receives the pieces: ``lattice`` (the packer's energy of the
        placed cell without torsions), ``torsion``, ``bend``, ``penalty_rotation``, ``penalty_c``
        (the closure terms, multipliers included), ``physical`` (the first three), and per chain
        the natural repeats, their mismatch to the cell's and the residual rotations."""
        x = np.asarray(x, dtype=float).ravel()
        self.n_evaluations += 1
        st = self.state(x)
        A, B, gam, c = st.lattice
        latn, frames, offs, Pn, built = self.build(x)
        E_lat, gP, glat = placed_energy_and_grad(self.pk, frames, Pn, latn, valences=self.valences)
        nc, n, P, R = self.nc, self.n, self.P, self.R
        cg, sg = np.cos(gam * _RAD), np.sin(gam * _RAD)
        g = np.zeros(self.n_var)
        gA = glat[0, 0]
        gB = cg * glat[1, 0] + sg * glat[1, 1]
        gG = _RAD * B * (-sg * glat[1, 0] + cg * glat[1, 1])
        gc = glat[2, 2]
        e_tors = e_bend = p_rot = p_c = 0.0
        base = 4 + 4 * nc
        k_ang = self.angle_stiffness
        for s in range(nc):
            chain, dX, dc, w, dw = built[s]
            sl = slice(s * n, (s + 1) * n)
            gPs = gP[sl]
            F = gPs.sum(axis=0)
            fa, fb = st.place[s, 0], st.place[s, 1]
            gA += fa * F[0]
            gB += fb * (cg * F[0] + sg * F[1])
            gG += _RAD * B * fb * (-sg * F[0] + cg * F[1])
            Q = Pn[sl] - offs[s]
            g[4 + 4 * s:8 + 4 * s] = [F @ latn[0], F @ latn[1],
                                      _RAD * float((gPs[:, 1] * Q[:, 0] - gPs[:, 0] * Q[:, 1]).sum()), F[2]]
            gv = np.einsum("ni,vni->v", gPs @ frames[s], dX)
            e_tors += self.pk.torsion_energy(st.tors[s])
            gv[:P] += self.pk.torsion_energy_grad(st.tors[s])
            d = (st.angs[s] - self.angles0) * _RAD
            e_bend += 0.5 * k_ang * float((d * d).sum()) * self.n_repeats
            gv[P:] += k_ang * d * _RAD * self.n_repeats
            lw, lc = self.lam_w[s], float(self.lam_c[s])
            p_rot += self.pen_rot * float(w @ w) + float(lw @ w)
            gv += dw @ (2.0 * self.pen_rot * w + lw)
            dcz = float(chain.c) - c
            p_c += self.k_c * dcz * dcz + lc * dcz
            tension = 2.0 * self.k_c * dcz + lc
            gv += tension * dc
            gc -= tension
            g[base + s * (P + R):base + (s + 1) * (P + R)] = gv
        g[:4] = [gA, gB, gG, gc]
        if parts is not None:
            parts.update(lattice=E_lat, torsion=e_tors, bend=e_bend, penalty_rotation=p_rot, penalty_c=p_c,
                         physical=E_lat + e_tors + e_bend,
                         repeats=np.array([float(b[0].c) for b in built]),
                         c_mismatch=np.array([float(b[0].c) - c for b in built]),
                         rotation_error=np.array([float(np.linalg.norm(b[3])) for b in built]))
        return E_lat + e_tors + e_bend + p_rot + p_c, g

    def terms(self, x) -> dict:
        """The pieces of the energy at ``x`` (see :meth:`value_and_grad`), plus the total."""
        parts: dict = {}
        total, g = self.value_and_grad(x, parts)
        parts["total"] = total
        parts["grad_max"] = float(np.abs(g).max())
        return parts

    def physical_energy(self, x) -> float:
        """Lattice energy plus every chain's torsion and bend energy: the number the crystal
        refinement reports for the packer's cell, kcal/mol per cell.  The closure terms are left
        out; at a closed structure (:meth:`relax_closed`) they are zero."""
        return float(self.terms(x)["physical"])

    # --- closure ---------------------------------------------------------------------------------
    def closure(self, x, chains=None):
        """``(h, J)``: the closure residuals of ``chains`` (default all) and their Jacobian.

        Four per chain, in the order ``(w_x, w_y, w_z, CLOSURE_C * (c_s - c))``: the rotation
        vector of the chain's repeat transform in degrees and its repeat mismatch to the cell's in
        A, the second scaled so that both are of one size.  ``J`` is ``(4 m, n_var)``."""
        x = np.asarray(x, dtype=float).ravel()
        st = self.state(x)
        ch = range(self.nc) if chains is None else [int(s) for s in np.atleast_1d(chains)]
        P, R = self.P, self.R
        base = 4 + 4 * self.nc
        h = np.zeros(4 * len(ch))
        J = np.zeros((4 * len(ch), self.n_var))
        for k, s in enumerate(ch):
            chains_, _, dc, w, dw = repeat_chains_and_jacobian(self.polymer, self.name, st.tors[s][None], st.angs[s], self.refs[s])
            cols = slice(base + s * (P + R), base + (s + 1) * (P + R))
            h[4 * k:4 * k + 3] = w[0]
            J[4 * k:4 * k + 3, cols] = dw[0].T
            h[4 * k + 3] = CLOSURE_C * (float(chains_[0].c) - st.lattice[3])
            J[4 * k + 3, cols] = CLOSURE_C * dc[0]
            J[4 * k + 3, 3] = -CLOSURE_C
        return h, J

    def closure_residual(self, x) -> tuple:
        """``(largest residual rotation in degrees, largest repeat mismatch in A)`` over the chains."""
        h, _ = self.closure(x)
        h = h.reshape(-1, 4)
        return float(np.linalg.norm(h[:, :3], axis=1).max()), float(np.abs(h[:, 3]).max() / CLOSURE_C)

    def project(self, x, free, chains=None, tol: float = 1e-9, maxiter: int = 30):
        """``x`` moved onto the closure surface of ``chains`` by the smallest change of the ``free``
        variables (in the scaled metric of :meth:`scale`): Newton steps ``-J^+ h``."""
        x = np.asarray(x, dtype=float).ravel().copy()
        free = np.unique(np.asarray(free, dtype=int))
        sc = self.scale(x)[free]
        for _ in range(maxiter):
            h, J = self.closure(x, chains)
            if float(np.abs(h).max()) < tol:
                break
            Jf = J[:, free] * sc
            x[free] -= sc * (Jf.T @ np.linalg.lstsq(Jf @ Jf.T, h, rcond=None)[0])
        return x

    # --- relaxation --------------------------------------------------------------------------
    def _bounds(self, x0, free, sc, angle_window, cell_window, c_window) -> list:
        """Per free variable ``(low, high)`` in scaled units: the backbone-angle window, and the
        cell edges, cell angle and repeat about their values at ``x0``."""
        bounds = [(None, None)] * len(free)
        ang = {int(i): k for k, i in enumerate(self.index("angs"))} if angle_window is not None else {}
        for j, i in enumerate(free):
            i = int(i)
            if i in ang:
                a0 = self.angles0[ang[i] % self.R]
                bounds[j] = ((a0 - angle_window) / sc[j], (a0 + angle_window) / sc[j])
            elif i in (0, 1):
                lo, hi = (0.0, None) if cell_window is None else (cell_window[0] * x0[i], cell_window[1] * x0[i])
                bounds[j] = (max(2.0, lo) / sc[j], None if hi is None else hi / sc[j])
            elif i == 2:
                bounds[j] = (55.0 / sc[j], 125.0 / sc[j])
            elif i == 3 and c_window is not None:
                bounds[j] = (c_window[0] * x0[i] / sc[j], c_window[1] * x0[i] / sc[j])
        return bounds

    def descend(self, x0, free=None, closed=None, steps: int = 200, max_step: float = 0.1, gtol: float = 1.5,
                angle_window: float | None = 8.0, cell_window=(0.6, 1.6), c_window=(0.85, 1.15), progress=None):
        """Damped descent (FIRE) on the closure surface of ``closed``, no step longer than
        ``max_step`` scaled units in any variable: what takes a strained start down to where a
        quasi-Newton method is safe, without leaving the basin the start is in.  Stops when the
        largest component of the gradient along the surface (kcal/mol per scaled unit) is under
        ``gtol`` or after ``steps`` evaluations.  ``progress(evaluation, energy, x)``, if given, is
        called after every evaluation.  Returns ``(x, largest gradient component, evaluations)``."""
        x0 = np.asarray(x0, dtype=float).ravel().copy()
        free = np.arange(self.n_var) if free is None else np.unique(np.asarray(free, dtype=int))
        closed = [] if closed is None else [int(s) for s in closed]
        sc = self.scale(x0)[free]
        bounds = self._bounds(x0, free, sc, angle_window, cell_window, c_window)
        lo = np.array([-np.inf if b[0] is None else b[0] for b in bounds])
        hi = np.array([np.inf if b[1] is None else b[1] for b in bounds])
        y = np.clip(x0[free] / sc, lo, hi)
        v = np.zeros_like(y)
        dt, dt_max, a_fire, n_pos = 0.05, 0.5, 0.1, 0
        gm, it = float("inf"), 0

        def full(yy):
            x = x0.copy()
            x[free] = yy * sc
            return x

        for it in range(1, steps + 1):
            x = full(y)
            e, g = self.value_and_grad(x)
            if progress is not None:
                progress(it, float(e), x)
            f = -g[free] * sc
            f[((y <= lo) & (f < 0.0)) | ((y >= hi) & (f > 0.0))] = 0.0  # a variable at its bound does not push through it
            if closed:
                Jf = self.closure(x, closed)[1][:, free] * sc
                _, sv, Vt = np.linalg.svd(Jf, full_matrices=False)
                Q = Vt[sv > 1e-10 * max(float(sv.max()), 1e-300)].T
                f -= Q @ (Q.T @ f)
                v -= Q @ (Q.T @ v)
            gm = float(np.abs(f).max())
            if gm < gtol:
                break
            if float(f @ v) > 0.0:
                n_pos += 1
                v = (1.0 - a_fire) * v + a_fire * (float(np.linalg.norm(v)) / max(float(np.linalg.norm(f)), 1e-300)) * f
                if n_pos > 5:
                    dt, a_fire = min(dt * 1.1, dt_max), a_fire * 0.99
            else:
                n_pos, v, dt, a_fire = 0, np.zeros_like(v), dt * 0.5, 0.1
            v = v + dt * f
            step = dt * v
            big = float(np.abs(step).max())
            if big > max_step:
                step *= max_step / big
            y = np.clip(y + step, lo, hi)
            if closed:
                for _ in range(12):
                    h, J = self.closure(full(y), closed)
                    if float(np.abs(h).max()) < 1e-9:
                        break
                    Jf = J[:, free] * sc
                    y = np.clip(y - Jf.T @ np.linalg.lstsq(Jf @ Jf.T, h, rcond=None)[0], lo, hi)
        return full(y), gm, it

    def relax(self, x0, free=None, closed=None, ftol: float = 1e-9, maxiter: int = 400, angle_window: float | None = 8.0,
              first_step: float | None = 0.3, presteps: int = 300, cell_window=(0.6, 1.6), c_window=(0.85, 1.15),
              trust: float | None = None, progress=None, patience: int | None = None, patience_tol: float = 1e-6):
        """Constrained minimisation (SLSQP) over the ``free`` variables (indices; default all), in
        the scaled variables of :meth:`scale`, with the closure of the chains in ``closed`` held as
        equality constraints.  ``closed`` defaults to every chain with a free internal coordinate,
        or to every chain when the repeat ``c`` is free.  ``angle_window`` bounds each backbone
        angle to the polymer's own value +/- that many degrees, the refinement's bound; ``None``
        leaves them free.  ``cell_window`` and ``c_window`` bound the cell edges and the repeat to
        those multiples of their starting values.  Returns ``(x, info)`` with the optimiser's
        message, its iteration and evaluation counts and the closure residuals at the result.

        SLSQP starts from a unit Hessian, so its first step is the (constrained) steepest-descent
        step at full length: from a strained start that is hundreds of degrees of torsion, enough
        to leave the basin the caller put the structure in.  With ``first_step`` set (scaled
        units, about 1 A of atom motion; ``None`` for the raw behaviour) a strained start is first
        taken down by :meth:`descend` (at most ``presteps`` evaluations, no step longer than a
        third of ``first_step``) until the gradient along the surface is under five times
        ``first_step``, and the objective handed to SLSQP is scaled so that its first step is no
        longer than ``first_step``: a factor of at most five, which leaves its unit Hessian
        inside the range of the real one.

        Neither of those bounds the later steps, so the basin this method ends in is not
        guaranteed to be the one it started in, and a structure it returns to another basin is not
        by itself shown to have had no rim (nor is a slow damped descent shown to have found one:
        examples/kink_torsion_scan.py is the test for a kink).  ``trust`` (scaled units) confines
        every free variable to that distance from where the optimiser starts, so one call cannot
        move a torsion by more than ``10 * trust`` degrees or a chain by more than ``trust`` A;
        ``info["at_trust_bound"]`` says whether a variable ended on the box.
        :meth:`relax_in_basin` repeats such calls until none does.

        ``progress(stage, evaluation, energy, x)``, if given, is called after every energy
        evaluation, ``stage`` being ``"descent"`` or ``"slsqp"``: a log line and a checkpoint for a
        relaxation that runs for hours (a cell of sixteen ten-monomer chains is seconds per
        evaluation and hundreds of evaluations).

        ``patience`` (a number of evaluations) with ``patience_tol`` (kcal/mol per cell) is a second
        way to stop: when the lowest energy found on the closure surface (residuals under 1e-3)
        has not fallen by ``patience_tol`` over the last ``patience`` evaluations made on it, the
        lowest point is returned with ``info["stalled"]`` true.  SLSQP's own test asks, from one
        iteration to the next, for an energy change and a closure residual both under ``ftol``
        (1e-9), and a large cell can go on failing it long after the energy has settled: the two
        end states of the sixteen-chain step were still iterating 220 and 80 evaluations after
        their lowest energy had stopped changing in the fourth decimal, at gradients of 7e-5 and
        1e-3 kcal/mol per scaled unit.  It is not rounding (the energy of that cell at points 1e-9
        scaled units apart scatters by 8e-12 kcal/mol about a straight line); which of the two
        conditions held it was not looked into.  Either way ``info["gradient_max"]`` says how
        stationary the result is: the largest component, in kcal/mol per scaled unit, of the
        gradient along the closure surface (a variable held by a bound it pushes against left out).

        The variables' units are those of :meth:`scale`, equal amounts of atom motion.  Units that
        even out the curvature instead (one factor per kind of variable, from the mean of
        ``z . H z`` over random signs ``z`` on that kind) were tried on the four-chain step and are
        worse, 1072 evaluations against 429 to the same minimum: one torsion turned alone swings
        the rest of its chain through the neighbours and is the stiffest direction there is
        (2000 kcal/mol per (10 deg)^2), while what a relaxation has to find is concerted turns of
        several, which no diagonal scaling can favour."""
        x0 = np.asarray(x0, dtype=float).ravel().copy()
        free = np.arange(self.n_var) if free is None else np.unique(np.asarray(free, dtype=int))
        if closed is None:
            fs = set(free.tolist())
            closed = [s for s in range(self.nc) if 3 in fs or fs.intersection(self.index("internal", s).tolist())]
        closed = [int(s) for s in closed]
        sc = self.scale(x0)[free]
        bounds = self._bounds(x0, free, sc, angle_window, cell_window, c_window)
        pre = 0
        if first_step is not None and presteps > 0:
            xd, _, pre = self.descend(x0, free, closed, steps=presteps, max_step=first_step / 3.0, gtol=5.0 * first_step,
                                      angle_window=angle_window, cell_window=cell_window, c_window=c_window,
                                      progress=None if progress is None else (lambda k, e, x: progress("descent", k, e, x)))
            start = xd
        else:
            start = x0

        def full(y):
            x = x0.copy()
            x[free] = y * sc
            return x

        cache: dict = {}
        ccache: dict = {}
        low = {"e": np.inf, "y": None, "mark": np.inf, "since": 0}

        class _Stalled(Exception):
            pass

        def hj(y):
            key = y.tobytes()
            if ccache.get("key") != key:
                h, J = self.closure(full(y), closed)
                ccache.update(key=key, h=h, J=J[:, free] * sc)
            return ccache["h"], ccache["J"]

        def fg(y):
            key = y.tobytes()
            if cache.get("key") != key:
                x = full(y)
                e, g = self.value_and_grad(x)
                cache.update(key=key, e=e, g=g[free] * sc, n=cache.get("n", 0) + 1)
                if progress is not None:
                    progress("slsqp", cache["n"], float(e), x)
                if patience is not None and (not closed or float(np.abs(hj(y)[0]).max()) < 1e-3):
                    # only points on the closure surface count, as a new lowest or as an evaluation without one:
                    # the optimiser's early iterates leave the surface, and their energies say nothing either way
                    if e < low["e"]:
                        low.update(e=float(e), y=y.copy())
                    if low["e"] < low["mark"] - patience_tol:
                        low.update(mark=low["e"], since=0)
                    else:
                        low["since"] += 1
                    if low["since"] >= patience:
                        raise _Stalled
            return cache["e"], cache["g"]

        cons = [{"type": "eq", "fun": lambda y: hj(y)[0], "jac": lambda y: hj(y)[1]}] if closed else []
        y0 = start[free] / sc
        if trust is not None:
            bounds = [(y - trust if lo is None else max(lo, y - trust), y + trust if hi is None else min(hi, y + trust)) for (lo, hi), y in zip(bounds, y0)]
        fs = 1.0
        if first_step is not None:
            g0 = fg(y0)[1].copy()
            if closed:  # the part of the gradient the constraints do not carry
                J0 = hj(y0)[1]
                g0 -= J0.T @ np.linalg.lstsq(J0 @ J0.T, J0 @ g0, rcond=None)[0]
            fs = max(1.0, float(np.abs(g0).max()) / float(first_step))
        its = [0]
        try:
            res = minimize(lambda y: fg(y)[0] / fs, y0, jac=lambda y: fg(y)[1] / fs, method="SLSQP", bounds=bounds,
                           constraints=cons, options={"ftol": ftol / fs, "maxiter": maxiter}, callback=lambda yk: its.__setitem__(0, its[0] + 1))
            y_end, message, nit, nfev, stalled = res.x, str(res.message), int(res.nit), int(res.nfev), False
        except _Stalled:
            y_end, nit, nfev, stalled = low["y"], its[0], int(cache["n"]), True
            message = f"stalled: the lowest energy on the closure surface has not fallen by {patience_tol:g} kcal/mol in {patience} evaluations"
        x = full(y_end)
        if closed:  # an iterate that stopped short of convergence is still returned periodic
            on_surface = np.intersect1d(free, np.concatenate([self.index("internal", closed), self.index("c")]))
            if on_surface.size:
                x = self.project(x, on_surface, closed)
        h = self.closure(x, closed)[0].reshape(-1, 4) if closed else np.zeros((1, 4))
        rot, dc = float(np.linalg.norm(h[:, :3], axis=1).max()), float(np.abs(h[:, 3]).max() / CLOSURE_C)
        at_bound = bool(trust is not None and float(np.abs(y_end - y0).max()) > trust * (1.0 - 1e-6))
        # how stationary the result is: the gradient along the closure surface, bounds that hold a variable left out
        parts: dict = {}
        _, g = self.value_and_grad(x, parts)
        yf, gf = x[free] / sc, g[free] * sc
        lo = np.array([-np.inf if b[0] is None else b[0] for b in bounds])
        hi = np.array([np.inf if b[1] is None else b[1] for b in bounds])
        gf[((yf <= lo + 1e-9) & (gf > 0.0)) | ((yf >= hi - 1e-9) & (gf < 0.0))] = 0.0
        if closed:
            Jf = self.closure(x, closed)[1][:, free] * sc
            _, sv, Vt = np.linalg.svd(Jf, full_matrices=False)
            Q = Vt[sv > 1e-10 * max(float(sv.max()), 1e-300)].T
            gf = gf - Q @ (Q.T @ gf)
        return x, {"message": message, "iterations": nit, "evaluations": nfev + pre,
                   "energy": float(parts["physical"]), "rotation_error_max": rot, "c_mismatch_max": dc, "closed_chains": closed,
                   "objective_scale": fs, "descent_evaluations": pre, "at_trust_bound": at_bound, "stalled": stalled,
                   "gradient_max": float(np.abs(gf).max()) if gf.size else 0.0}

    def relax_in_basin(self, x0, free=None, closed=None, trust: float = 0.3, rounds: int = 100, callback=None, **kw):
        """:meth:`relax` in a trust box of half-width ``trust`` (scaled units), repeated from each
        result until no variable ends on its box and the optimiser reports success: quasi-Newton
        steps that cannot leave the basin they start in unless it has no rim.  Each round ends at
        the lowest point of its box that the optimiser finds, so the energy never rises from one
        round to the next, and a minimum behind a barrier wider than the box is out of reach.
        ``callback(round, x, info)`` is called after every round.  Returns ``(x, info)``; ``info``
        has the last round's fields, ``rounds``, the summed ``iterations`` and ``evaluations``,
        and ``converged`` (false if the rounds ran out with a variable still on its box)."""
        x = np.asarray(x0, dtype=float).ravel().copy()
        kw = {"presteps": 0, "first_step": None, **kw}
        info, its, evs, done = {}, 0, 0, False
        for k in range(int(rounds)):
            x, info = self.relax(x, free=free, closed=closed, trust=trust, **kw)
            its += info["iterations"]
            evs += info["evaluations"]
            if callback is not None:
                callback(k, x, info)
            if not info["at_trust_bound"] and "success" in info["message"].lower():
                done = True
                break
        return x, {**info, "rounds": k + 1, "iterations": its, "evaluations": evs, "converged": done}


# ------------------------------------------------------------------- pendant tilts
def tilt_relaxation(model: MixedModel, x, polymer=None, k_backbone: float | None = None, gtol: float = 1e-7,
                    perturb: float = 0.0, seed: int = 0):
    """:func:`polyfind.pendants.tilt_relaxation` for a mixed cell: ``(dE_per_cell, tilts_deg, groups)``.

    The search starts from zero tilts, where a structure with a symmetry has no gradient along a
    tilt that would break it, so a symmetric arrangement is returned untilted whether or not it is
    stable.  ``perturb`` (degrees) starts instead from random tilts of up to that size (``seed``):
    the same result from both starts is the evidence that the untilted groups are a minimum.

    Every CH2/CF2/CHF group of every chain may turn rigidly about its carbon, with the carbons,
    the cell, the bond lengths and the X-C-X angles fixed, the belt sites following their groups
    and each chain's valence bends pinned at that chain's own geometry.  ``dE_per_cell`` is the
    energy released (<= 0, kcal/mol per cell)."""
    from . import pendants as PD

    polymer = model.polymer if polymer is None else polymer
    latn, frames, offs, P0, built = model.build(x)
    pk = model.pk
    cz = float(latn[2, 2])
    vals = [PD.pinned_chain_valence(b[0], polymer, k_backbone, cz=cz)[0] for b in built]
    J, atoms = PD.belt_site_constraint(pk, P0, cz=float(latn[2, 2]))
    el = np.asarray(list(pk.elements) * pk.n_chains)
    sym = [str(e) for e in el[atoms]]
    groups = PD.pendant_groups(latn, P0[atoms], sym)
    if not groups:
        raise ValueError("no backbone atom with two backbone neighbours and two single-atom pendants")
    E0, _, _ = placed_energy_and_grad(pk, frames, P0, latn, valences=vals)
    h = 1e-6

    def eg(t):
        dX = PD.group_displacements(groups, t, len(atoms))
        E, g, _ = placed_energy_and_grad(pk, frames, P0 + (J @ dX).reshape(-1, 3), latn, valences=vals)
        gX = J.T @ np.asarray(g, dtype=float).ravel()
        grad = np.empty_like(t)
        for i in range(len(t)):
            tp, tm = t.copy(), t.copy()
            tp[i] += h
            tm[i] -= h
            grad[i] = gX @ (PD.group_displacements(groups, tp, len(atoms)) - PD.group_displacements(groups, tm, len(atoms))) / (2 * h)
        return float(E), grad

    t0 = np.zeros(3 * len(groups))
    if perturb:
        t0 = np.radians(np.random.default_rng(seed).uniform(-float(perturb), float(perturb), t0.size))
    res = minimize(eg, t0, jac=True, method="L-BFGS-B",
                   options={"maxiter": 2000, "gtol": gtol, "ftol": 1e-15})
    if float(res.fun) > float(E0) + 1e-9:  # the start is admissible; never report a rise
        return 0.0, np.zeros((len(groups), 3)), groups
    return float(res.fun) - float(E0), np.degrees(res.x.reshape(len(groups), 3)), groups


# ----------------------------------------------------------------------------- NEB
@dataclass
class NEBResult:
    images: np.ndarray        # (M, n_var), end points included
    energies: np.ndarray      # (M,)
    force_max: float          # largest scaled NEB force component on the moving images
    climbing: int             # index of the climbing image, -1 before it is switched on
    iterations: int
    converged: bool
    history: list


def neb(evaluate, images, scale=None, free=None, k_spring: float = 1.0, climb_after: int = 30,
        fmax: float = 0.05, maxiter: int = 300, max_step: float = 0.15, callback=None,
        constraint=None, constraint_tol: float = 1e-8, climb_force: float | None = None, snapshot=None) -> NEBResult:
    """Nudged elastic band with the improved tangent and a climbing image, relaxed by FIRE.

    ``evaluate(X)`` takes the moving images ``(m, n_var)`` and returns ``(energies (m,), gradients
    (m, n_var))``, so the caller can evaluate them in parallel.  ``images`` (M, n_var) is the
    initial path with both end points, which stay fixed and must be relaxed already; ``scale``
    (n_var,) converts each variable to a common length (the band's metric), ``free`` restricts
    which variables move (the others are held at their per-image values).  ``k_spring`` is in
    kcal/mol per scaled length squared.  The band converges when the largest scaled force
    component (perpendicular true force plus parallel spring force; the full reversed force on
    the climbing image) falls under ``fmax``.  The climbing image is switched on after
    ``climb_after`` iterations and, if ``climb_force`` is given, not before the band's largest force
    has fallen under it: an image of a straight initial path can start on a steric clash far above
    the saddle, and one switched to climbing there climbs the clash.

    ``constraint(x)``, if given, returns ``(h (k,), J (k, n_var))`` for one image: equality
    constraints every image keeps (the chain closure, :meth:`MixedModel.closure`).  The band then
    lives on the constraint surface: gradients, tangents and velocities are projected onto its
    tangent space at each image, and after every step each image is returned to the surface by
    Newton steps along the constraint normals (the smallest move in the scaled metric).  The end
    points and the initial images must lie on it already.

    ``callback(iteration, energies, force_max, climbing)`` is called once per iteration, and
    ``snapshot(iteration, images)`` with the whole band as it stands (end points included): what a
    caller writes to disk to resume a band that runs for hours.
    """
    X = np.array(images, dtype=float)
    M, nv = X.shape
    free = np.arange(nv) if free is None else np.unique(np.asarray(free, dtype=int))
    sc = np.ones(nv) if scale is None else np.asarray(scale, dtype=float)
    scf = sc[free]
    Y = X[:, free] / scf  # scaled, free

    def full(i):
        x = X[i].copy()
        x[free] = Y[i] * scf
        return x

    def normals(i):
        """Orthonormal columns spanning the constraint normals of image ``i`` in the scaled free
        variables, and the residual."""
        h, J = constraint(full(i))
        Jf = np.asarray(J, dtype=float)[:, free] * scf
        _, sv, Vt = np.linalg.svd(Jf, full_matrices=False)
        keep = sv > 1e-10 * max(float(sv.max()), 1e-300)
        return Vt[keep].T, np.asarray(h, dtype=float), Jf

    def retract(i):
        for _ in range(12):
            h, J = constraint(full(i))
            h = np.asarray(h, dtype=float)
            if float(np.abs(h).max()) < constraint_tol:
                return
            Jf = np.asarray(J, dtype=float)[:, free] * scf
            Y[i] -= Jf.T @ np.linalg.lstsq(Jf @ Jf.T, h, rcond=None)[0]

    e_end, _ = evaluate(X[[0, M - 1]])
    E = np.empty(M)
    E[0], E[-1] = float(e_end[0]), float(e_end[1])
    V = np.zeros((M - 2, len(free)))
    dt, dt_max, a_fire, n_pos = 0.05, 0.5, 0.1, 0
    climbing, history = -1, []
    converged, it, fm = False, 0, float("inf")
    for it in range(1, maxiter + 1):
        Xm = X[1:-1].copy()
        Xm[:, free] = Y[1:-1] * scf
        e, g = evaluate(Xm)
        E[1:-1] = np.asarray(e, dtype=float)
        G = np.asarray(g, dtype=float)[:, free] * scf  # scaled gradient of the moving images
        Q = [None] * M
        if constraint is not None:
            for i in range(1, M - 1):
                Q[i] = normals(i)[0]
                G[i - 1] -= Q[i] @ (Q[i].T @ G[i - 1])
                V[i - 1] -= Q[i] @ (Q[i].T @ V[i - 1])
        if climbing < 0 and it > climb_after and (climb_force is None or fm < climb_force):
            climbing = int(np.argmax(E[1:-1])) + 1
        Fn = np.empty_like(G)
        for i in range(1, M - 1):
            tp, tm = Y[i + 1] - Y[i], Y[i] - Y[i - 1]
            ep, e0, em = E[i + 1], E[i], E[i - 1]
            if ep > e0 > em:
                tau = tp
            elif ep < e0 < em:
                tau = tm
            else:
                dmax, dmin = max(abs(ep - e0), abs(em - e0)), min(abs(ep - e0), abs(em - e0))
                tau = tp * dmax + tm * dmin if ep > em else tp * dmin + tm * dmax
            if Q[i] is not None:
                tau = tau - Q[i] @ (Q[i].T @ tau)
            tau = tau / max(float(np.linalg.norm(tau)), 1e-300)
            f_true = -G[i - 1]
            if i == climbing:
                Fn[i - 1] = f_true - 2.0 * float(f_true @ tau) * tau
            else:
                f_spring = k_spring * (float(np.linalg.norm(tp)) - float(np.linalg.norm(tm))) * tau
                Fn[i - 1] = f_true - float(f_true @ tau) * tau + f_spring
        fm = float(np.abs(Fn).max())
        history.append({"iteration": it, "force_max": fm, "barrier": float(E.max() - E[0]), "climbing": climbing})
        if callback is not None:
            callback(it, E.copy(), fm, climbing)
        if snapshot is not None:
            Xs = X.copy()
            Xs[:, free] = Y * scf
            snapshot(it, Xs)
        if fm < fmax and (climbing > 0 or climb_after >= maxiter):
            converged = True
            break
        # FIRE on the whole band
        power = float((Fn * V).sum())
        if power > 0.0:
            n_pos += 1
            fn, vn = float(np.linalg.norm(Fn)), float(np.linalg.norm(V))
            V = (1.0 - a_fire) * V + a_fire * (vn / max(fn, 1e-300)) * Fn
            if n_pos > 5:
                dt, a_fire = min(dt * 1.1, dt_max), a_fire * 0.99
        else:
            n_pos, V, dt, a_fire = 0, np.zeros_like(V), dt * 0.5, 0.1
        V = V + dt * Fn
        step = dt * V
        big = float(np.abs(step).max())
        if big > max_step:
            step *= max_step / big
        Y[1:-1] += step
        if constraint is not None:
            for i in range(1, M - 1):
                retract(i)
    X[:, free] = Y * scf
    return NEBResult(images=X, energies=E, force_max=fm, climbing=climbing, iterations=it, converged=converged, history=history)
