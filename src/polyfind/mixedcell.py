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

**One repeat for all chains.**  The cell has one ``c``, and a chain built from its internal
coordinates has its own natural repeat ``c_s`` and, off a commensurate conformation, a residual
rotation ``w_s`` of its repeat transform.  Both are held by penalties, ``k_c (c_s - c)^2`` and
``pen_rot |w_s|^2`` (the second is the refinement's own, with the same default); at a relaxed
structure they are reported (:meth:`MixedModel.terms`) and are not part of the physical energy.
A chain shorter than its neighbours, which is what a kink pair makes, therefore reaches the
cell's repeat through its own bends and torsions and carries the tension that costs -- the
periodic model's finite-length effect, which falls as the inverse repeat length.

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
    "tilt_relaxation",
]

_FLIP = np.diag([1.0, -1.0, -1.0])
_RAD = np.pi / 180.0


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


class MixedModel:
    """Energy and analytic gradient of a mixed cell in the chain builder's internal coordinates.

    ``packer`` is built with ``n_chains`` equal to the number of chains here, from any chain of
    the sequence (the tables depend on the topology alone); ``refs`` is one
    :class:`~polyfind.pack.PeriodicChain` per chain, the frame each chain's coordinates are
    aligned to (the builder's own frame is a principal axis of the cross section and flips
    discontinuously; see :func:`polyfind.pack._align_about_z`), fixed for the life of the model so
    a setting angle means the same thing at every point of a path.  ``angle_stiffness`` is the
    refinement's backbone-angle restraint (kcal/mol/rad^2,
    :func:`polyfind.refine.resolve_angle_stiffness`), ``pen_rot`` its commensurability penalty per
    degree^2 and ``k_c`` the penalty tying each chain's natural repeat to the cell's (kcal/mol/A^2).
    """

    def __init__(self, polymer, name: str, packer, flips, refs, angle_stiffness: float,
                 pen_rot: float = 5.0, k_c: float = 2000.0, valences=None):
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
        placed cell without torsions), ``torsion``, ``bend``, ``penalty_rotation``, ``penalty_c``,
        ``physical`` (the first three), and per chain the natural repeats, their mismatch to the
        cell's and the residual rotations."""
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
            p_rot += self.pen_rot * float(w @ w)
            gv += 2.0 * self.pen_rot * (dw @ w)
            dcz = float(chain.c) - c
            p_c += self.k_c * dcz * dcz
            gv += 2.0 * self.k_c * dcz * dc
            gc -= 2.0 * self.k_c * dcz
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
        refinement reports for the packer's cell, kcal/mol per cell.  The closure penalties are
        left out; at a relaxed structure they are small and :meth:`terms` says how small."""
        return float(self.terms(x)["physical"])

    # --- relaxation --------------------------------------------------------------------------
    def relax(self, x0, free=None, gtol: float = 2e-3, maxiter: int = 400, angle_window: float | None = 8.0):
        """L-BFGS-B over the ``free`` variables (indices; default all), in the scaled variables of
        :meth:`scale`.  ``angle_window`` bounds each backbone angle to the polymer's own value +/-
        that many degrees, the refinement's bound; ``None`` leaves them free.  Returns
        ``(x, info)`` with ``info`` the optimiser's message, iteration count and final gradient."""
        x0 = np.asarray(x0, dtype=float).ravel().copy()
        free = np.arange(self.n_var) if free is None else np.unique(np.asarray(free, dtype=int))
        sc = self.scale(x0)[free]
        bounds = [(None, None)] * len(free)
        if angle_window is not None:
            ang = {int(i): k for k, i in enumerate(self.index("angs"))}
            for j, i in enumerate(free):
                if int(i) in ang:
                    a0 = self.angles0[ang[int(i)] % self.R]
                    bounds[j] = ((a0 - angle_window) / sc[j], (a0 + angle_window) / sc[j])
        for j, i in enumerate(free):
            if int(i) in (0, 1):
                bounds[j] = (2.0 / sc[j], None)
            elif int(i) == 2:
                bounds[j] = (55.0 / sc[j], 125.0 / sc[j])

        def fg(y):
            x = x0.copy()
            x[free] = y * sc
            e, g = self.value_and_grad(x)
            return e, g[free] * sc

        res = minimize(fg, x0[free] / sc, jac=True, method="L-BFGS-B", bounds=bounds,
                       options={"gtol": gtol, "ftol": 1e-12, "maxiter": maxiter, "maxfun": 3 * maxiter})
        x = x0.copy()
        x[free] = res.x * sc
        return x, {"message": str(res.message), "iterations": int(res.nit), "evaluations": int(res.nfev),
                   "grad_max_scaled": float(np.abs(res.jac).max()) if res.jac is not None else float("nan"),
                   "energy": float(res.fun)}


# ------------------------------------------------------------------- pendant tilts
def tilt_relaxation(model: MixedModel, x, polymer=None, k_backbone: float | None = None, gtol: float = 1e-7):
    """:func:`polyfind.pendants.tilt_relaxation` for a mixed cell: ``(dE_per_cell, tilts_deg, groups)``.

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

    res = minimize(eg, np.zeros(3 * len(groups)), jac=True, method="L-BFGS-B",
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
        fmax: float = 0.05, maxiter: int = 300, max_step: float = 0.15, callback=None) -> NEBResult:
    """Nudged elastic band with the improved tangent and a climbing image, relaxed by FIRE.

    ``evaluate(X)`` takes the moving images ``(m, n_var)`` and returns ``(energies (m,), gradients
    (m, n_var))``, so the caller can evaluate them in parallel.  ``images`` (M, n_var) is the
    initial path with both end points, which stay fixed and must be relaxed already; ``scale``
    (n_var,) converts each variable to a common length (the band's metric), ``free`` restricts
    which variables move (the others are held at their per-image values).  ``k_spring`` is in
    kcal/mol per scaled length squared.  The band converges when the largest scaled force
    component (perpendicular true force plus parallel spring force; the full reversed force on
    the climbing image) falls under ``fmax``.
    """
    X = np.array(images, dtype=float)
    M, nv = X.shape
    free = np.arange(nv) if free is None else np.unique(np.asarray(free, dtype=int))
    sc = np.ones(nv) if scale is None else np.asarray(scale, dtype=float)
    scf = sc[free]
    Y = X[:, free] / scf  # scaled, free
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
        if climbing < 0 and it > climb_after:
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
    X[:, free] = Y * scf
    return NEBResult(images=X, energies=E, force_max=fm, climbing=climbing, iterations=it, converged=converged, history=history)
