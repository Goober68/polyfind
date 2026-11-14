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

**An applied field.**  A uniform field set on the packer (:meth:`polyfind.pack.CrystalPacker.set_field`,
V/A) is carried as it is in the packer's own cell: ``-mu . E`` on the placed charges, and the field
added to the permanent field the induced dipoles answer, so the energy is the electric enthalpy
``U - mu_total . E`` of the cell at that macroscopic field and its derivative with respect to the
field is minus the cell's dipole, static plus induced (:func:`cell_dipole`).  The dipole is the
sum over the atoms as placed, every chain whole: between two states joined by a continuous path
its difference is defined without a polarization quantum.

**What is not here.**  Charge flux (refused), and the pendant tilts, which
:func:`tilt_relaxation` adds afterwards exactly as :func:`polyfind.pendants.tilt_relaxation` does
for the packer's cell.
"""
from __future__ import annotations

import os
import time
from dataclasses import dataclass, replace

import numpy as np
from scipy.optimize import minimize

from . import backend as bk
from . import chainjac as CJ
from . import phonon as PH
from .ewald import exclusion_correction_pairs
from .forcefield import COULOMB
from .linegroup import _windows, repeat_chains_and_jacobian
from .pack import (
    EV_TO_KCAL,
    PeriodicChain,
    _block_atoms,
    _derfc_pos,
    _erfc_pos,
    _finish_repeat_chains,
    _template,
)

__all__ = [
    "MixedModel",
    "MixedState",
    "cell_dipole",
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
def placed_energy_and_grad(packer, frames, Pn, latn, valences=None, xp=None, timings=None):
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

    ``xp`` chooses how the sums are organised.  ``None`` is the loop over chains and chain pairs
    below, on the host, block by block: the form every stored number was computed with.  An
    array module -- NumPy itself, or a device's -- evaluates the same terms as a fixed number of
    array operations whatever the number of chains (:func:`_placed_batched`): the pair sums of
    all chain pairs as one stack of blocks on that module, the bonded exclusions and the
    valence terms of all chains from flat lists.  The two agree to rounding, not bit for bit
    (the sums are taken in another order).  ``timings``, a dict, receives the seconds each part
    of the batched form took.
    """
    n, N, nc = packer.n, packer.N, packer.n_chains
    Pn = np.asarray(Pn, dtype=float).reshape(N, 3)
    latn = np.asarray(latn, dtype=float).reshape(3, 3)
    frames = np.asarray(frames, dtype=float).reshape(nc, 3, 3)
    if getattr(packer, "_flux", None) is not None:
        raise ValueError("charge flux is not carried by the mixed cell: every chain has the repeat's fixed charges")
    if abs(latn[0, 2]) + abs(latn[1, 2]) + abs(latn[2, 0]) + abs(latn[2, 1]) > 1e-12:
        raise ValueError("the mixed cell wants c along z and the other two lattice vectors in the plane")
    if xp is not None:
        return _placed_batched(packer, frames, Pn, latn, valences, xp, timings)
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

    if getattr(packer, "_field_on", False):
        e_f, g_f = _applied_field(packer, Pn)
        E += e_f
        gP += g_f
    for s in range(nc):
        gP[sl[s]] += gX[s] @ frames[s].T
    return E, gP, glat


def _applied_field(packer, Pn) -> tuple:
    """``(-mu_static . E, its derivative with respect to the placed coordinates)`` in kcal/mol, for
    the packer's set field (V/A).  The induced dipoles' part is in the dipole solve, whose permanent
    field carries the applied one (:meth:`polyfind.pack.CrystalPacker._polarize`).  The lattice
    derivative at fixed Cartesian coordinates is zero, and so is the net force on a neutral chain."""
    f = np.asarray(packer.field, dtype=float)
    q = np.asarray(packer._q_cell, dtype=float)
    return -EV_TO_KCAL * float(q @ (Pn @ f)), -EV_TO_KCAL * q[:, None] * f[None, :]


def cell_dipole(packer, frames, Pn, latn) -> tuple:
    """``(static, induced)`` dipole of the placed cell in e.A, at the packer's set field: the sum
    of ``q r`` over the atoms as placed (every chain whole), and the sum of the induced dipoles
    (zeros without a polarizable model).  Their sum is ``-dE/dE_applied`` of
    :func:`placed_energy_and_grad`, in eV per V/A."""
    Pn = np.asarray(Pn, dtype=float).reshape(packer.N, 3)
    q = np.asarray(packer._q_cell, dtype=float)
    static = q @ Pn
    induced = np.zeros(3)
    if packer.polarizable is not None and packer._ewald is not None:
        flips = np.asarray(frames, dtype=float).reshape(packer.n_chains, 3, 3)[:, 2, 2] < 0.0
        induced = packer._polarize(Pn, q, np.asarray(latn, dtype=float).reshape(3, 3), flips.astype(float))[1].sum(axis=0)
    return static, induced


# ------------------------------------------- the same energy as a fixed number of array operations
HOST_BLOCK_ELEMS = 120_000    # pair elements per stack of blocks on the host: what stays in cache
DEVICE_BLOCK_ELEMS = 8_000_000  # and on a device, where one large stack is one launch


class _Stopwatch:
    """Adds the time since the last mark to ``timings[name]``.  On a device the mark waits for the
    device first: a device call returns when the work is queued, not when it is done."""

    def __init__(self, timings, xp=np):
        self.timings, self.xp = timings, xp
        self.last = time.perf_counter() if timings is not None else 0.0

    def __call__(self, name: str) -> None:
        if self.timings is None:
            return
        if self.xp is not np:
            bk.synchronize(self.xp)
        now = time.perf_counter()
        self.timings[name] = self.timings.get(name, 0.0) + (now - self.last)
        self.last = now


def _pair_kernel(xp, r2, tab, rc2: float, alpha: float):
    """``(V, 2 dV/dr^2)`` of the packer's pair potential on squared distances, on any array module.

    The expressions of :meth:`polyfind.pack.CrystalPacker._pair_energy_and_dv` with its mask;
    ``tab`` is ``(A, B, qq or None, qq_f, const)`` on that module (``qq`` is ``None`` under Ewald,
    where the kernel carries no Coulomb term).  The factor two of ``dV/d(separation vector) =
    2 dV/dr^2 * (separation vector)`` is folded in."""
    A, B, qq, qqf, cst = tab
    mask = (r2 < rc2) & (r2 > 1e-8)
    r2s = xp.where(mask, r2, 1.0)
    rr = xp.sqrt(r2s)
    ir6 = 1.0 / (r2s * r2s * r2s)
    y = A * ir6
    v = ir6 * (y - B) + qqf * rr + cst
    t2 = 6.0 * ir6 * (B - 2.0 * y) / r2s + qqf / rr
    if qq is not None:
        ar = alpha * rr
        ef = _erfc_pos(ar, xp)
        v = v + qq * (ef / rr)
        t2 = t2 + qq * (alpha * _derfc_pos(ar, xp) / rr - ef / r2s) / rr
    return xp.where(mask, v, 0.0), xp.where(mask, t2, 0.0)


def _module_key(xp) -> str:
    return getattr(xp, "__name__", repr(xp))


def _pair_tables(packer, xp) -> tuple:
    """The repeat's pair tables ``(A, B, qq or None, qq_f, const)``, each ``(n, n)``, on ``xp``.
    Kept on the packer per array module: they depend on the topology and the charges alone."""
    cache = packer.__dict__.setdefault("_mixed_pair_tables", {})
    key = (_module_key(xp), id(packer._q_cell))
    hit = cache.get(key)
    if hit is None:
        n = packer.n
        q = np.asarray(packer._q_cell[:n], dtype=float)
        qq, qqf, cst = PH._charge_tables(packer, q, q)
        A = np.asarray(bk.to_numpy(packer._A_nn), dtype=float)
        B = np.asarray(bk.to_numpy(packer._B_nn), dtype=float)
        cache.clear()
        hit = cache[key] = (xp.asarray(A), xp.asarray(B), None if packer._ewald is not None else xp.asarray(qq),
                            xp.asarray(qqf), xp.asarray(cst))
    return hit


def _scale_stack(packer, K: int, xp) -> tuple:
    """The bonded-scale stack ``(2K+1, n, n)`` on ``xp`` with the two index grids that gather from it."""
    cache = packer.__dict__.setdefault("_mixed_scale_stacks", {})
    key = (_module_key(xp), int(K))
    hit = cache.get(key)
    if hit is None:
        n = packer.n
        ar = xp.arange(n)
        if len(cache) > 4:
            cache.clear()
        hit = cache[key] = (xp.asarray(packer._scale_column_np(K)), ar[None, :, None], ar[None, None, :])
    return hit


def _chain_pair_blocks(P3, latn, rc: float) -> tuple:
    """Every (chain, chain, lateral image) that can hold a pair of sites inside the cutoff.

    ``(s, t, i, j, shift (nb, 2))``: chain ``s`` against chain ``t`` moved by ``i A + j B``, each
    unordered block once -- ``s < t`` with any image, and a chain with its own images in one
    half plane (the other half is the same blocks seen from the other side).  The test of the
    loop form, then one more that needs no pair distances: the nearest site of one chain to the
    moved centre of the other is at least the other's radius nearer than any of its sites.
    """
    nc = P3.shape[0]
    A2, B2 = latn[0, :2], latn[1, :2]
    area = abs(A2[0] * B2[1] - A2[1] * B2[0])
    wa, wb = area / float(np.linalg.norm(B2)), area / float(np.linalg.norm(A2))
    Minv = np.linalg.inv(np.array([A2, B2]))  # lateral Cartesian -> fractional (row vectors)
    xy = P3[:, :, :2]
    cen = xy.mean(axis=1)
    rad = np.sqrt(((xy - cen[:, None, :]) ** 2).sum(axis=-1)).max(axis=1)
    s, t = np.triu_indices(nc)
    reach = rc + rad[s] + rad[t]
    d0 = cen[s] - cen[t]
    f0 = d0 @ Minv
    ri, rj = int(np.ceil(float(reach.max()) / wa)) + 1, int(np.ceil(float(reach.max()) / wb)) + 1
    di, dj = (g.ravel() for g in np.meshgrid(np.arange(-ri, ri + 1), np.arange(-rj, rj + 1), indexing="ij"))
    I = np.round(f0[:, 0]).astype(int)[:, None] + di[None, :]
    J = np.round(f0[:, 1]).astype(int)[:, None] + dj[None, :]
    shx, shy = I * A2[0] + J * B2[0], I * A2[1] + J * B2[1]
    keep = np.hypot(d0[:, 0:1] - shx, d0[:, 1:2] - shy) < reach[:, None]
    keep &= (s != t)[:, None] | (I > 0) | ((I == 0) & (J > 0))
    pi, im = np.nonzero(keep)
    bs, bt, bi, bj = s[pi], t[pi], I[pi, im], J[pi, im]
    sh = np.stack([shx[pi, im], shy[pi, im]], axis=1)
    da = np.sqrt(((xy[bs] - (cen[bt] + sh)[:, None, :]) ** 2).sum(axis=-1)).min(axis=1) - rad[bt]
    db = np.sqrt(((xy[bt] - (cen[bs] - sh)[:, None, :]) ** 2).sum(axis=-1)).min(axis=1) - rad[bs]
    near = (da < rc) & (db < rc)
    return bs[near], bt[near], bi[near], bj[near], sh[near]


def _pair_sums(packer, P3, latn, flips, K: int, xp, block_elems: int | None = None):
    """The pair energy of the placed cell, ``(E, dE/dP (nc, n, 3), dE/dlattice (3, 3))``, on ``xp``.

    Two stacks of ``(n, n)`` blocks through one kernel.  The blocks of
    :func:`_chain_pair_blocks`, at weight one, each site pair at its nearest image along z and
    as many further images either side as the cutoff can reach (none once ``c`` is at least
    twice the cutoff); and every chain against itself at no lateral shift, at half weight with
    the bonded scales gathered by image -- for a flipped chain by minus the image, since the
    stack is the chain's own frame's.  The lattice derivative follows the loop form's
    convention: the same-site column gives its ``[2, 2]`` entry alone.
    """
    nc, n = P3.shape[0], P3.shape[1]
    cz, rc, rc2 = float(latn[2, 2]), float(packer.rc), float(packer.rc2)
    alpha = float(packer.alpha)
    if block_elems is None:
        block_elems = HOST_BLOCK_ELEMS if xp is np else DEVICE_BLOCK_ELEMS
    per = max(1, int(block_elems) // (n * n))
    bs, bt, bi, bj, sh = _chain_pair_blocks(P3, latn, rc)
    tab = _pair_tables(packer, xp)
    Pd = xp.asarray(P3)
    reach_k = max(0, int(np.ceil(rc / cz - 0.5)))  # images either side of the nearest that the cutoff reaches
    offsets = [float(o) for o in range(-reach_k, reach_k + 1)]
    inv_cz = 1.0 / cz
    E = 0.0
    g3 = xp.zeros((nc, n * 3))
    lat01 = xp.zeros((2, 3))
    lat2 = xp.zeros(3)
    chain_ids = xp.arange(nc)
    ij = np.stack([bi, bj]).astype(float)

    def forces(t2, dx, dy, dz):
        fx, fy, fz = t2 * dx, t2 * dy, t2 * dz
        on_first = xp.stack([fx.sum(axis=2), fy.sum(axis=2), fz.sum(axis=2)], axis=2)
        on_second = xp.stack([fx.sum(axis=1), fy.sum(axis=1), fz.sum(axis=1)], axis=2)
        return fx, fy, fz, on_first, on_second

    for lo in range(0, bs.size, per):
        hi = min(lo + per, bs.size)
        sd, td, shd = xp.asarray(bs[lo:hi]), xp.asarray(bt[lo:hi]), xp.asarray(sh[lo:hi])
        Pa, Pb = Pd[sd], Pd[td]
        dx = Pa[:, :, None, 0] - (Pb[:, :, 0] + shd[:, 0][:, None])[:, None, :]
        dy = Pa[:, :, None, 1] - (Pb[:, :, 1] + shd[:, 1][:, None])[:, None, :]
        dz0 = Pa[:, :, None, 2] - Pb[:, None, :, 2]
        dxy2 = dx * dx + dy * dy
        kn = xp.round(dz0 * inv_cz)
        gs = gt = None
        for off in offsets:
            kel = kn + off if off else kn
            dz = dz0 - kel * cz
            v, t2 = _pair_kernel(xp, dxy2 + dz * dz, tab, rc2, alpha)
            E = E + v.sum()
            fx, fy, fz, a, b = forces(t2, dx, dy, dz)
            gs, gt = (a, b) if gs is None else (gs + a, gt + b)
            # the image shift is i A + j B + k C, and dE/d(shift) = -sum of the pair forces
            lat2 = lat2 - xp.stack([(kel * fx).sum(), (kel * fy).sum(), (kel * fz).sum()])
        nb = hi - lo
        hot_s = (chain_ids[:, None] == sd[None, :]).astype(float)
        hot_t = (chain_ids[:, None] == td[None, :]).astype(float)
        g3 = g3 + hot_s @ gs.reshape(nb, n * 3) - hot_t @ gt.reshape(nb, n * 3)
        lat01 = lat01 - xp.asarray(np.ascontiguousarray(ij[:, lo:hi])) @ gs.sum(axis=1)

    Sd, row, col = _scale_stack(packer, K, xp)
    sign = xp.asarray(np.where(np.asarray(flips, dtype=bool), -1.0, 1.0))
    lat22 = 0.0
    for lo in range(0, nc, per):
        hi = min(lo + per, nc)
        Pc = Pd[lo:hi]
        dx = Pc[:, :, None, 0] - Pc[:, None, :, 0]
        dy = Pc[:, :, None, 1] - Pc[:, None, :, 1]
        dz0 = Pc[:, :, None, 2] - Pc[:, None, :, 2]
        dxy2 = dx * dx + dy * dy
        kn = xp.round(dz0 * inv_cz)
        sg = sign[lo:hi][:, None, None]
        for off in offsets:
            kel = kn + off if off else kn
            dz = dz0 - kel * cz
            v, t2 = _pair_kernel(xp, dxy2 + dz * dz, tab, rc2, alpha)
            inside = abs(kel) <= K  # the loop form's own range of images
            kidx = (xp.where(inside, sg * kel, 0.0) + K).astype(np.int64)
            wgt = xp.where(inside, Sd[kidx, row, col], 0.0)
            E = E + 0.5 * (v * wgt).sum()
            fx, fy, fz, a, b = forces(t2 * wgt, dx, dy, dz)
            g3[lo:hi] += 0.5 * (a - b).reshape(hi - lo, n * 3)
            lat22 = lat22 - 0.5 * (kel * fz).sum()

    glat = np.zeros((3, 3))
    glat[:2] = bk.to_numpy(lat01)
    glat[2] = bk.to_numpy(lat2)
    glat[2, 2] += float(lat22)
    return float(E), np.asarray(bk.to_numpy(g3), dtype=float).reshape(nc, n, 3), glat


def _valence_groups(packer, valences, nc: int) -> list:
    """``[(ChainValence, chain indices)]``: the chains that share one valence object, together."""
    if valences is None:
        return [] if packer._valence is None else [(packer._valence, np.arange(nc))]
    groups: dict = {}
    for s, val in enumerate(valences):
        if val is not None:
            groups.setdefault(id(val), (val, []))[1].append(s)
    return [(val, np.array(chains)) for val, chains in groups.values()]


def _placed_batched(packer, frames, Pn, latn, valences, xp, timings=None):
    """:func:`placed_energy_and_grad` as a fixed number of array operations: see there."""
    n, N, nc = packer.n, packer.N, packer.n_chains
    cz, rc = float(latn[2, 2]), float(packer.rc)
    K = int(np.ceil(rc / cz)) + 1
    flips = frames[:, 2, 2] < 0.0
    P3 = Pn.reshape(nc, n, 3)
    mark = _Stopwatch(timings, xp)
    E = 0.0
    gP = np.zeros((N, 3))
    glat = np.zeros((3, 3))
    g3 = gP.reshape(nc, n, 3)

    # --- valence, on the placed coordinates: a flipped chain's neighbouring repeats lie the other way
    sign = np.where(flips, -1.0, 1.0)
    for val, chains in _valence_groups(packer, valences, nc):
        Ev, gv, gcv = val.energy_and_grad_batch(P3[chains], cz, image_sign=sign[chains])
        E += Ev
        g3[chains] += gv
        glat[2, 2] += gcv
    mark("valence")

    # --- every pair of sites inside the cutoff, same-site columns included
    Ep, gp, glp = _pair_sums(packer, P3, latn, flips, K, xp)
    E += Ep
    g3 += gp
    glat += glp
    mark("pair sums")

    # --- electrostatics of an Ewald packer
    if packer._ewald is not None:
        q_cell = np.asarray(packer._q_cell, dtype=float)
        pairs = packer._exclusion_pairs(K, flips)
        if packer.polarizable is not None:
            e_es, _, _, (gPe, glate, _, gce) = packer._polarize(Pn, q_cell, latn, flips.astype(float), grad=True, pairs=pairs)
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
        mark("electrostatics")
        ex_e, ex_g, ex_c = exclusion_correction_pairs(Pn, q_cell, cz, pairs, pref=COULOMB / packer.eps_r,
                                                      charge_damping=packer._ewald.charge_damping)
        E += ex_e
        gP += ex_g
        glat[2, 2] += ex_c
        mark("exclusions")
    if getattr(packer, "_field_on", False):
        e_f, g_f = _applied_field(packer, Pn)
        E += e_f
        gP += g_f
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

    **Two forms of one evaluation.**  :meth:`value_and_grad` is either the loop over chains and
    chain pairs every stored number was computed with, or the batched form
    (:meth:`_value_and_grad_batched`): the same energy to rounding as a fixed number of array
    operations, with the pair sums on an array module.  ``xp`` chooses: ``None`` (the default)
    follows the packer's Ewald sums -- on a device when they are (``POLYFIND_EWALD_DEVICE=cuda``
    and a cell of at least ``Ewald.device_min_sites`` sites), where the loops would leave the
    device idle most of the evaluation, otherwise the loops, or the batched form on the host if
    ``POLYFIND_MIXED_BATCHED=1``; an array module asks for the batched form on it (NumPy
    included), and ``"loops"`` for the loops wherever the Ewald sums are.  ``timings``, set to a
    dict, collects the seconds each part of the batched form takes.
    """

    def __init__(self, polymer, name: str, packer, flips, refs, angle_stiffness: float,
                 pen_rot: float = 0.0, k_c: float = 0.0, valences=None, xp=None):
        self.polymer, self.name, self.pk = polymer, name, packer
        self.xp = xp
        self.timings = None
        self.nc, self.n = int(packer.n_chains), int(packer.n)
        # The pair sums of a mixed cell are host loops in float64 that borrow the packer's kernel.
        # On a machine with a GPU the packer's batched kernels move to the device in float32 by
        # themselves (polyfind.backend), and this cell would then be evaluated block by block
        # across the bus at seven digits: wrong for a relaxation and slow.  Say so instead.
        if getattr(packer, "xp", np) is not np or np.dtype(getattr(packer, "_dt", np.float64)) != np.float64:
            raise ValueError(
                "a mixed cell evaluates its pair sums on the host in float64, and this packer was built for the "
                "batched device backend (float32): set POLYFIND_DEVICE=cpu, or call polyfind.backend.set_backend('cpu'), "
                "before building the packer.  The Ewald sums and the dipole solve go to a device separately, with "
                "POLYFIND_EWALD_DEVICE=cuda.")
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
        st = self.state(x)
        chains, tr, _ = self._build_light(st.tors, st.angs, np.arange(self.nc))
        w = np.degrees(tr["w"])
        dc = np.array([float(ch.c) for ch in chains]) - float(st.lattice[3])
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

    def dipole(self, x) -> dict:
        """The cell's dipole at ``x`` in e.A, at the packer's set field: ``static`` (the placed
        charges), ``induced`` and ``total`` (:func:`cell_dipole`)."""
        latn, frames, _, Pn, _ = self.build(np.asarray(x, dtype=float).ravel())
        static, induced = cell_dipole(self.pk, frames, Pn, latn)
        return {"static": static, "induced": induced, "total": static + induced}

    def _topology(self) -> dict:
        """What the light build reads from the repeat's topology, resolved once."""
        top = self.__dict__.get("_light_topology")
        if top is None:
            tpl = _template(self.polymer, self.P)
            _, atoms0 = _block_atoms(tpl, self.P, 0)
            sx, sy = _windows(self.P)
            back = np.asarray(tpl.backbone)
            top = self._light_topology = {
                "atoms0": np.asarray(atoms0, dtype=int), "sx": sx, "sy": sy, "n_window": int(back[sx].size),
                "window": np.concatenate([back[sx], back[sy]]), "refs": np.stack(self.refs)}
        return top

    def _build_light(self, tors, angs, chains, sums: bool = False):
        """``(chains, transform, oligomer)`` of the chains ``chains`` built together.

        The :class:`~polyfind.pack.PeriodicChain` objects are those of :meth:`build` (the same
        oligomer, the same batched tail), and ``transform`` is the repeat transform of the two
        backbone windows -- ``c``, the rotation vector ``w`` and the repeat axis with their
        derivatives, in radians.  Its derivative comes from the Jacobian of the window atoms
        alone (:func:`polyfind.chainjac.window_transform_jacobian`: bit for bit what the whole
        chain's Jacobian gives), or with ``sums`` from running sums along the backbone
        (:func:`polyfind.chainjac.backbone_transform_jacobian`: the same to rounding, a fifth of
        the time), which is what the batched evaluation uses.  The derivative of the coordinates
        is not formed: a gradient is carried back through ``oligomer`` ``(coords, bb, owner,
        is_pendant)`` by :func:`polyfind.chainjac.oligomer_vjp`.
        """
        top = self._topology()
        tors, angs = np.atleast_2d(np.asarray(tors, dtype=float)), np.atleast_2d(np.asarray(angs, dtype=float))
        coords, bb, owner, is_pendant = CJ.oligomer(self.polymer, tors, angs)
        built = _finish_repeat_chains(self.polymer, self.name, tors, coords, 0.5, top["refs"][np.asarray(chains, dtype=int)])
        if sums:
            tr = CJ.backbone_transform_jacobian(bb, top["sx"], top["sy"], self.P, self.R)
        else:
            jw = CJ.oligomer_jacobian(coords, bb, owner, is_pendant, self.P, self.R, atoms=top["window"])
            nw = top["n_window"]
            tr = CJ.window_transform_jacobian(bb[:, top["sx"]], bb[:, top["sy"]], jw[:, :, :nw], jw[:, :, nw:])
        return built, tr, (coords, bb, owner, is_pendant)

    def batched_module(self):
        """The array module of the batched evaluation, or ``None`` for the loops (see the class)."""
        if isinstance(self.xp, str):
            if self.xp != "loops":
                raise ValueError(f"MixedModel.xp is an array module, None or 'loops', not {self.xp!r}")
            return None
        if self.xp is not None:
            return self.xp
        ew = self.pk._ewald
        if ew is not None:
            n_sites = self.pk.N if self.pk._es_idx is None else len(self.pk._es_idx)
            xp = ew.array_module(n_sites)
            if xp is not np:
                return xp
        return np if os.environ.get("POLYFIND_MIXED_BATCHED", "") == "1" else None

    def _value_and_grad_batched(self, x, parts, xp):
        """:meth:`value_and_grad` as a fixed number of array operations, the pair sums on ``xp``.

        Every chain is built in one batch; the placed cell is evaluated by
        :func:`placed_energy_and_grad` with ``xp``; and the gradient with respect to the torsions
        and angles is the torque of the Cartesian gradient about each bond and each angle's
        normal (:func:`polyfind.chainjac.oligomer_vjp`), not a product with the Jacobian of the
        coordinates, which is forty variables times a hundred sites per chain and was a third of
        the host's work.
        """
        mark = _Stopwatch(self.timings)
        st = self.state(x)
        B, gam, c = st.lattice[1:]
        latn = lattice_rows(*st.lattice)
        nc, n, P, R = self.nc, self.n, self.P, self.R
        top = self._topology()
        built, tr, (coords, bb, owner, is_pendant) = self._build_light(st.tors, st.angs, np.arange(nc), sums=True)
        Xa = np.stack([ch.coords for ch in built])
        cs = np.array([float(ch.c) for ch in built])
        phi = st.place[:, 2] * _RAD
        cp, sp = np.cos(phi), np.sin(phi)
        frames = np.zeros((nc, 3, 3))
        frames[:, 0, 0], frames[:, 0, 1], frames[:, 1, 0], frames[:, 1, 1], frames[:, 2, 2] = cp, -sp, sp, cp, 1.0
        frames[self.flips, :, 1:] *= -1.0  # Rz(phi) after the (x, -y, -z) flip
        offs = st.place[:, 0:1] * latn[0] + st.place[:, 1:2] * latn[1]
        offs[:, 2] += st.place[:, 3]
        Q = np.einsum("sni,sji->snj", Xa, frames)
        Pn = (Q + offs[:, None, :]).reshape(nc * n, 3)
        mark("build")
        E_lat, gP, glat = placed_energy_and_grad(self.pk, frames, Pn, latn, valences=self.valences, xp=xp, timings=self.timings)
        mark = _Stopwatch(self.timings)
        cg, sg = np.cos(gam * _RAD), np.sin(gam * _RAD)
        g3 = gP.reshape(nc, n, 3)
        F = g3.sum(axis=1)
        fa, fb = st.place[:, 0], st.place[:, 1]
        g = np.zeros(self.n_var)
        gA = glat[0, 0] + float((fa * F[:, 0]).sum())
        gB = cg * glat[1, 0] + sg * glat[1, 1] + float((fb * (cg * F[:, 0] + sg * F[:, 1])).sum())
        gG = _RAD * B * (-sg * glat[1, 0] + cg * glat[1, 1]) + _RAD * B * float((fb * (-sg * F[:, 0] + cg * F[:, 1])).sum())
        gc = glat[2, 2]
        g[4:4 + 4 * nc] = np.stack([F @ latn[0], F @ latn[1],
                                    _RAD * (g3[:, :, 1] * Q[:, :, 0] - g3[:, :, 0] * Q[:, :, 1]).sum(axis=1), F[:, 2]], axis=1).ravel()
        # the Cartesian gradient in each chain's own frame, carried back to its torsions and angles
        G = np.einsum("sni,sij->snj", g3, frames)
        gv, Gp = CJ.aligned_block_vjp(coords[:, top["atoms0"]], tr["axis"], tr["daxis"], top["refs"], G)
        gv += CJ.oligomer_vjp(coords, bb, owner, is_pendant, P, R, top["atoms0"], Gp)
        V1, V2, V3 = self.pk.torsion
        ph = np.deg2rad(st.tors)
        e_tors = float((0.5 * (V1 * (1 + np.cos(ph)) + V2 * (1 - np.cos(2 * ph)) + V3 * (1 + np.cos(3 * ph)))).sum())
        gv[:, :P] += 0.5 * (-V1 * np.sin(ph) + 2.0 * V2 * np.sin(2 * ph) - 3.0 * V3 * np.sin(3 * ph)) * (np.pi / 180.0)
        d = (st.angs - self.angles0) * _RAD
        e_bend = 0.5 * self.angle_stiffness * float((d * d).sum()) * self.n_repeats
        gv[:, P:] += self.angle_stiffness * d * _RAD * self.n_repeats
        w, dw = np.degrees(tr["w"]), np.degrees(tr["dw"])
        p_rot = self.pen_rot * float((w * w).sum()) + float((self.lam_w * w).sum())
        gv += np.einsum("svi,si->sv", dw, 2.0 * self.pen_rot * w + self.lam_w)
        dcz = cs - c
        p_c = float((self.k_c * dcz * dcz + self.lam_c * dcz).sum())
        tension = 2.0 * self.k_c * dcz + self.lam_c
        gv += tension[:, None] * tr["dc"]
        gc -= float(tension.sum())
        g[4 + 4 * nc:] = gv.ravel()
        g[:4] = [gA, gB, gG, gc]
        if parts is not None:
            parts.update(lattice=E_lat, torsion=e_tors, bend=e_bend, penalty_rotation=p_rot, penalty_c=p_c,
                         physical=E_lat + e_tors + e_bend, repeats=cs.copy(), c_mismatch=dcz,
                         rotation_error=np.linalg.norm(w, axis=1))
        mark("assembly")
        return E_lat + e_tors + e_bend + p_rot + p_c, g

    # --- energy ------------------------------------------------------------------------------
    def value_and_grad(self, x, parts: dict | None = None):
        """Energy with the closure penalties, and its gradient with respect to every variable.

        ``parts``, if a dict is given, receives the pieces: ``lattice`` (the packer's energy of the
        placed cell without torsions), ``torsion``, ``bend``, ``penalty_rotation``, ``penalty_c``
        (the closure terms, multipliers included), ``physical`` (the first three), and per chain
        the natural repeats, their mismatch to the cell's and the residual rotations."""
        x = np.asarray(x, dtype=float).ravel()
        self.n_evaluations += 1
        xp = self.batched_module()
        if xp is not None:
            return self._value_and_grad_batched(x, parts, xp)
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
        ch = np.arange(self.nc) if chains is None else np.array([int(s) for s in np.atleast_1d(chains)], dtype=int)
        P, R = self.P, self.R
        base = 4 + 4 * self.nc
        h = np.zeros(4 * len(ch))
        J = np.zeros((4 * len(ch), self.n_var))
        if len(ch) == 0:
            return h, J
        # the repeat transform reads two windows of the backbone: all chains in one build, and the
        # derivative of those atoms alone (the same numbers as the whole chain's Jacobian gives);
        # beside the batched evaluation, from running sums along the backbone (to rounding)
        built, tr, _ = self._build_light(st.tors[ch], st.angs[ch], ch, sums=self.batched_module() is not None)
        w, dw = np.degrees(tr["w"]), np.degrees(tr["dw"])
        for k, s in enumerate(ch):
            cols = slice(base + s * (P + R), base + (s + 1) * (P + R))
            h[4 * k:4 * k + 3] = w[k]
            J[4 * k:4 * k + 3, cols] = dw[k].T
            h[4 * k + 3] = CLOSURE_C * (float(built[k].c) - st.lattice[3])
            J[4 * k + 3, cols] = CLOSURE_C * tr["dc"][k]
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
              trust: float | None = None, progress=None, patience: int | None = None, patience_tol: float = 1e-6, linear=None):
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

        ``linear`` = ``(a, b)`` holds ``a . y = b`` as one more equality, ``y`` being the free
        variables in scaled units and ``a`` a vector of their number: the relaxation then finds the
        minimum in a hyperplane (:func:`saddle_polish`).  The damped pre-descent, which knows only
        the closure, is skipped; ``gradient_max`` leaves the hyperplane's normal out as it does the
        closure's; and ``info["linear_force"]`` is the multiplier of the new constraint, the
        derivative of the constrained minimum's energy with respect to ``b``.

        In a cell of many chains this method is slow and stops short, for two reasons that
        :meth:`relax_on_surface` removes (trial points off the closure surface, and the
        conditioning): use that one there.

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
        a_lin = None if linear is None else np.asarray(linear[0], dtype=float).ravel()
        if a_lin is not None and a_lin.size != free.size:
            raise ValueError(f"linear: {a_lin.size} coefficients for {free.size} free variables")
        if first_step is not None and presteps > 0 and linear is None:
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
        if a_lin is not None:
            b_lin = float(linear[1])
            cons.append({"type": "eq", "fun": lambda y: np.array([float(a_lin @ y) - b_lin]), "jac": lambda y: a_lin[None, :]})
        y0 = start[free] / sc
        if trust is not None:
            bounds = [(y - trust if lo is None else max(lo, y - trust), y + trust if hi is None else min(hi, y + trust)) for (lo, hi), y in zip(bounds, y0)]
        fs = 1.0
        if first_step is not None:
            g0 = fg(y0)[1].copy()
            if closed or a_lin is not None:  # the part of the gradient the constraints do not carry
                J0 = np.vstack(([hj(y0)[1]] if closed else []) + ([a_lin[None, :]] if a_lin is not None else []))
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
        normals = ([self.closure(x, closed)[1][:, free] * sc] if closed else []) + ([a_lin[None, :]] if a_lin is not None else [])
        multiplier = None
        if normals:
            Jf = np.vstack(normals)
            if a_lin is not None:  # gradient = (closure normals) mu + a lambda: lambda is dE/db of the constrained minimum
                multiplier = float(np.linalg.lstsq(Jf.T, gf, rcond=None)[0][-1])
            _, sv, Vt = np.linalg.svd(Jf, full_matrices=False)
            Q = Vt[sv > 1e-10 * max(float(sv.max()), 1e-300)].T
            gf = gf - Q @ (Q.T @ gf)
        return x, {"message": message, "iterations": nit, "evaluations": nfev + pre,
                   "energy": float(parts["physical"]), "rotation_error_max": rot, "c_mismatch_max": dc, "closed_chains": closed,
                   "objective_scale": fs, "descent_evaluations": pre, "at_trust_bound": at_bound, "stalled": stalled,
                   "gradient_max": float(np.abs(gf).max()) if gf.size else 0.0, "linear_force": multiplier}

    def relax_on_surface(self, x0, free=None, closed=None, metric=None, linear=None, gtol: float = 1e-4, ftol: float = 1e-9,
                         maxiter: int = 500, memory: int = 30, longest: float = 0.5, first_step: float | None = 0.3,
                         presteps: int = 0, angle_window: float | None = 8.0, cell_window=(0.6, 1.6), c_window=(0.85, 1.15),
                         progress=None, patience: int | None = None, patience_tol: float = 1e-6):
        """The minimisation of :meth:`relax` done on the closure surface: a limited-memory
        quasi-Newton method in the surface's tangent space, every trial point returned to the
        surface (and to the hyperplane of ``linear``) before its energy is evaluated.

        :meth:`relax` hands the closure to SLSQP as equality constraints, and SLSQP tries straight
        steps: a step along the tangent leaves the surface at second order, a chain that does not
        close meets its own image, and the energy there says nothing about the energy on the
        surface (a step of Newton's length in the first-saddle hyperplane of the four-chain step:
        10.331 off the surface at a closure residual of 0.017 degrees, 10.279 once returned to it,
        from 10.290).  So its steps stay short and most of its evaluations are spent on trial
        points it then rejects (2.2 to 2.8 per iteration in cells of four to sixteen chains).  Here
        nothing is evaluated off the surface, a step is judged by the energy alone, and nearly
        every step is accepted.

        ``metric`` is a matrix ``T`` (free x free): the method works in ``w`` with ``dy = T dw``,
        ``y`` being the free variables in the units of :meth:`scale`.  With ``T`` an inverse square
        root of the Hessian the first step is Newton's; :func:`chain_metric` builds one from one
        block per chain, which is what the conditioning needs (the constrained Hessian of the
        sixteen-chain sheet has eigenvalues from 0.03 to 420 in the scaled variables, 0.02 to 2
        in those).  Without it the method is the same one in the scaled variables.

        ``gtol`` is the largest component of the gradient along the surface at which it stops
        (kcal/mol per scaled unit; what :meth:`relax` reports as ``gradient_max``), ``ftol`` a
        fall of the energy (kcal/mol per cell) under which two successive steps also end it,
        ``longest`` the longest move of any variable in one step (scaled units), ``first_step``
        that of a step taken without curvature information, ``memory`` the number of steps
        remembered.  ``presteps`` > 0 runs :meth:`descend` first, as :meth:`relax` does.
        ``progress(stage, evaluation, energy, x)``, ``patience`` and ``patience_tol`` are those of
        :meth:`relax`; every evaluation here is on the surface.  The windows on the backbone
        angles and the cell are not enforced (near a minimum none is active):
        ``info["outside_windows"]`` counts the variables that ended outside theirs.  Returns
        ``(x, info)`` with the fields of :meth:`relax`."""
        from collections import deque

        x0 = np.asarray(x0, dtype=float).ravel().copy()
        free = np.arange(self.n_var) if free is None else np.unique(np.asarray(free, dtype=int))
        if closed is None:
            fset = set(free.tolist())
            closed = [s for s in range(self.nc) if 3 in fset or fset.intersection(self.index("internal", s).tolist())]
        closed = [int(s) for s in closed]
        sc = self.scale(x0)[free]
        n = free.size
        bounds = self._bounds(x0, free, sc, angle_window, cell_window, c_window)
        a_lin = None if linear is None else np.asarray(linear[0], dtype=float).ravel()
        if a_lin is not None and a_lin.size != n:
            raise ValueError(f"linear: {a_lin.size} coefficients for {n} free variables")
        b_lin = None if linear is None else float(linear[1])
        T = None if metric is None else np.asarray(metric, dtype=float)
        if T is not None and T.shape != (n, n):
            raise ValueError(f"metric: shape {T.shape} for {n} free variables")
        Tinv = None if T is None else np.linalg.inv(T)
        pre, start = 0, x0
        if presteps > 0 and linear is None:
            fstep = 0.3 if first_step is None else float(first_step)
            start, _, pre = self.descend(x0, free, closed, steps=presteps, max_step=fstep / 3.0, gtol=5.0 * fstep,
                                         angle_window=angle_window, cell_window=cell_window, c_window=c_window,
                                         progress=None if progress is None else (lambda k, e, x: progress("descent", k, e, x)))

        def full(y):
            x = x0.copy()
            x[free] = y * sc
            return x

        def constraint(y):
            vals, rows = [], []
            if closed:
                h, J = self.closure(full(y), closed)
                vals.append(h)
                rows.append(J[:, free] * sc)
            if a_lin is not None:
                vals.append(np.array([float(a_lin @ y) - b_lin]))
                rows.append(a_lin[None, :])
            if not rows:
                return np.zeros(0), np.zeros((0, n))
            return np.concatenate(vals), np.vstack(rows)

        def retract(y):  # back to the constraints by Newton steps, each the shortest in the metric; (None, None) if they do not converge
            for _ in range(40):
                c, N = constraint(y)
                if c.size == 0 or float(np.abs(c).max()) < 1e-9:
                    return y, N
                if not np.all(np.isfinite(c)):
                    break
                Nw = N if T is None else N @ T
                dw = Nw.T @ np.linalg.lstsq(Nw @ Nw.T, c, rcond=None)[0]
                y = y - (dw if T is None else T @ dw)
            return None, None

        def normals(N, in_metric: bool):  # orthonormal columns spanning the constraints' normals, in w or in y
            if N.shape[0] == 0:
                return np.zeros((n, 0))
            M = N @ T if (in_metric and T is not None) else N
            _, sv, Vt = np.linalg.svd(M, full_matrices=False)
            return Vt[sv > 1e-10 * max(float(sv.max()), 1e-300)].T

        def along(gy, N):  # the gradient along the surface, in the scaled variables
            Qy = normals(N, False)
            return gy - Qy @ (Qy.T @ gy)

        count = {"n": 0}
        low = {"e": np.inf, "y": None, "mark": np.inf, "since": 0}

        class _Stalled(Exception):
            pass

        def energy(y):
            x = full(y)
            e, g = self.value_and_grad(x)
            count["n"] += 1
            if progress is not None:
                progress("surface", count["n"], float(e), x)
            if e < low["e"]:
                low.update(e=float(e), y=y.copy())
            if patience is not None:
                if low["e"] < low["mark"] - patience_tol:
                    low.update(mark=low["e"], since=0)
                else:
                    low["since"] += 1
                if low["since"] >= patience:
                    raise _Stalled
            return float(e), g[free] * sc

        pairs: deque = deque(maxlen=int(memory))

        def two_loop(v):
            q, hist = v.copy(), []
            for s_, y_ in reversed(pairs):
                rho = 1.0 / float(s_ @ y_)
                al = rho * float(s_ @ q)
                hist.append((al, rho, s_, y_))
                q -= al * y_
            if pairs:
                s_, y_ = pairs[-1]
                q *= float(s_ @ y_) / float(y_ @ y_)
            for al, rho, s_, y_ in reversed(hist):
                q += (al - rho * float(y_ @ q)) * s_
            return q

        y, N = retract(start[free] / sc)
        if y is None:
            raise RuntimeError("the start cannot be returned to the closure surface")
        message, stalled, its, small = "iteration limit reached", False, 0, 0
        try:
            e, gy = energy(y)
            Q = normals(N, True)
            r = gy if T is None else T.T @ gy
            r = r - Q @ (Q.T @ r)
            for its in range(1, int(maxiter) + 1):
                if float(np.abs(along(gy, N)).max()) < gtol:
                    message = "Optimization terminated successfully: the gradient along the surface is under gtol"
                    break
                d = -two_loop(r)
                d -= Q @ (Q.T @ d)
                slope = float(r @ d)
                if not slope < -1e-12 * float(np.linalg.norm(r) * np.linalg.norm(d)):
                    pairs.clear()
                    d, slope = -r, -float(r @ r)
                step = d if T is None else T @ d
                big = max(float(np.abs(step).max()), 1e-300)
                cap = longest if (pairs or first_step is None) else min(longest, float(first_step))
                alpha = min(1.0, cap / big)
                accepted = False
                for _ in range(30):
                    yt, Nt = retract(y + alpha * step)
                    if yt is None:
                        alpha *= 0.25
                    else:
                        et, gyt = energy(yt)
                        if et <= e + 1e-4 * alpha * slope:
                            accepted = True
                            break
                        denom = 2.0 * (et - e - alpha * slope)
                        new = -slope * alpha * alpha / denom if denom > 0.0 and np.isfinite(denom) else 0.25 * alpha
                        alpha = min(max(new, 0.1 * alpha), 0.5 * alpha)
                    if alpha * big < 1e-10:
                        break
                if not accepted:
                    if pairs:  # the curvature remembered is no guide here: once more along the steepest descent
                        pairs.clear()
                        continue
                    message = "no lower point along the steepest descent"
                    break
                Qt = normals(Nt, True)
                rt = gyt if T is None else T.T @ gyt
                rt = rt - Qt @ (Qt.T @ rt)
                sw = (yt - y) if T is None else Tinv @ (yt - y)
                sw = sw - Qt @ (Qt.T @ sw)
                yw = rt - (r - Qt @ (Qt.T @ r))
                if float(sw @ yw) > 1e-10 * float(np.linalg.norm(sw) * np.linalg.norm(yw)):
                    pairs.append((sw, yw))
                small = small + 1 if e - et < ftol else 0
                y, N, Q, e, gy, r = yt, Nt, Qt, et, gyt, rt
                if small >= 2:
                    message = "Optimization terminated successfully: the energy has stopped falling (ftol)"
                    break
        except _Stalled:
            y, stalled = low["y"], True
            N = constraint(y)[1]
            message = f"stalled: the lowest energy on the closure surface has not fallen by {patience_tol:g} kcal/mol in {patience} evaluations"
        x = full(y)
        h = self.closure(x, closed)[0].reshape(-1, 4) if closed else np.zeros((1, 4))
        rot, dc = float(np.linalg.norm(h[:, :3], axis=1).max()), float(np.abs(h[:, 3]).max() / CLOSURE_C)
        parts: dict = {}
        _, g = self.value_and_grad(x, parts)
        yf, gf = x[free] / sc, g[free] * sc
        lo = np.array([-np.inf if b[0] is None else b[0] for b in bounds])
        hi = np.array([np.inf if b[1] is None else b[1] for b in bounds])
        outside = int(np.count_nonzero((yf < lo - 1e-9) | (yf > hi + 1e-9)))
        multiplier = None
        if N.shape[0]:
            if a_lin is not None:
                multiplier = float(np.linalg.lstsq(N.T, gf, rcond=None)[0][-1])
            gf = along(gf, N)
        return x, {"message": message, "iterations": its, "evaluations": count["n"] + pre,
                   "energy": float(parts["physical"]), "rotation_error_max": rot, "c_mismatch_max": dc, "closed_chains": closed,
                   "objective_scale": 1.0, "descent_evaluations": pre, "at_trust_bound": False, "stalled": stalled,
                   "gradient_max": float(np.abs(gf).max()) if gf.size else 0.0, "linear_force": multiplier, "outside_windows": outside}

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


# --------------------------------------------------------------------- a saddle, settled
def saddle_polish(model: MixedModel, x, tangent, free=None, closed=None, tol: float = 0.01, first: float = 0.05,
                  longest: float = 0.3, rounds: int = 10, log=None, on_surface: bool = False, **relax_kw):
    """The stationary point near ``x`` that is a maximum along ``tangent`` and a minimum across it.

    A band's climbing image finds where a pass is; this settles it.  With ``a`` the unit vector of
    ``tangent`` (a difference of two states, in the variables' own units) in the scaled free
    variables, taken along the closure surface at ``x``, and ``s = a . (y - y(x))``: the energy is
    minimised in the hyperplane of fixed ``s`` by :meth:`MixedModel.relax` (every direction
    across the path at once, soft ones included, by the quasi-Newton method the end states are
    relaxed with), which also gives ``dE_min/ds``, the multiplier of the hyperplane; and ``s`` is
    moved by secant steps to where that derivative vanishes.  At that point the gradient has no
    component in the hyperplane (a minimum there) and none along its normal: a stationary point,
    and a first-order saddle if the derivative falls through zero as ``s`` grows, which the secant
    slope says.  The normal need not be the unstable direction exactly, only close enough for the
    constrained minima to form a continuous curve.

    ``first`` is the first move of ``s`` (scaled units, uphill), ``longest`` the longest later one,
    ``tol`` the largest ``|dE_min/ds|`` accepted (kcal/mol per scaled unit), ``rounds`` the most
    hyperplanes tried; ``relax_kw`` goes to :meth:`MixedModel.relax` (``patience`` is the one a
    large cell needs).  With ``trust`` among them each minimisation is confined to that distance
    (scaled units) from its start, which keeps the optimiser's trial steps out of the far field,
    and one that ends on its box is continued from there in the same hyperplane, up to five
    times, before its energy and force are used.  ``on_surface`` minimises with
    :meth:`MixedModel.relax_on_surface` instead (``relax_kw`` then goes to it: ``metric``, ``gtol``),
    which is what a cell of sixteen chains needs: the first-saddle hyperplane of the sixteen-chain
    step is at +10.6596 after 64 evaluations that way, and was at +10.6688 after 246 and +10.6637
    after 1435, still falling, with :meth:`MixedModel.relax`.  ``log(text)`` gets a line per
    hyperplane.  Returns ``(x, info)``: ``energy``,
    ``force`` (the last ``dE_min/ds``), ``gradient_max`` (across the path), ``slope``
    (``d2E_min/ds2`` from the last two hyperplanes, negative at a saddle; ``None`` if only one was
    needed), ``displacement`` (the final ``s``), ``converged``, ``evaluations`` and ``rounds``, the
    list of hyperplanes tried."""
    x = np.asarray(x, dtype=float).ravel().copy()
    free = np.arange(model.n_var) if free is None else np.unique(np.asarray(free, dtype=int))
    closed = [] if closed is None else [int(c) for c in closed]
    sc = model.scale(x)[free]
    a = np.asarray(tangent, dtype=float).ravel()[free] / sc
    if closed:  # along the closure surface: the part of the tangent the closure's normals do not carry
        Jf = model.closure(x, closed)[1][:, free] * sc
        _, sv, Vt = np.linalg.svd(Jf, full_matrices=False)
        Q = Vt[sv > 1e-10 * max(float(sv.max()), 1e-300)].T
        a = a - Q @ (Q.T @ a)
    norm = float(np.linalg.norm(a))
    if norm < 1e-12:
        raise ValueError("the tangent has no component along the closure surface in the free variables")
    a /= norm
    b0 = float(a @ (x[free] / sc))
    inner = np.intersect1d(free, model.index("internal", closed)) if closed else np.zeros(0, dtype=int)  # what returns a moved point to the surface
    s, tried, slope, done, total = 0.0, [], None, False, 0
    kw = {"presteps": 0, **relax_kw}
    minimise = model.relax_on_surface if on_surface else model.relax
    info: dict = {}
    for k in range(int(rounds)):
        x, info = minimise(x, free=free, closed=closed, linear=(a, b0 + s), **kw)
        total += info["evaluations"]
        for _ in range(5):  # a boxed minimisation that stopped on its box has not found the hyperplane's minimum yet
            if not info.get("at_trust_bound"):
                break
            x, more = minimise(x, free=free, closed=closed, linear=(a, b0 + s), **kw)
            total += more["evaluations"]
            info = {**more, "evaluations": info["evaluations"] + more["evaluations"]}
        tried.append({"displacement": s, "energy": info["energy"], "force": info["linear_force"], "gradient_max": info["gradient_max"],
                      "evaluations": info["evaluations"], "message": info["message"]})
        if log is not None:
            log(f"hyperplane {k + 1}: s {s:+.4f}, energy {info['energy']:.5f}, dE/ds {info['linear_force']:+.4f}, across {info['gradient_max']:.1e}, "
                f"{info['evaluations']} evaluations; {info['message'][:40]}")
        g = float(info["linear_force"])
        if len(tried) > 1:
            ds_last = s - tried[-2]["displacement"]
            slope = (g - tried[-2]["force"]) / ds_last if abs(ds_last) > 1e-12 else slope
        if abs(g) < tol:
            done = True
            break
        step = -g / slope if slope is not None and slope < 0.0 else float(np.sign(g)) * first  # uphill until the curve is seen to bend down
        step = float(np.clip(step, -longest, longest))
        s += step
        x = x.copy()
        x[free] += sc * a * step
        if inner.size:
            x = model.project(x, inner, closed)
    return x, {"energy": info["energy"], "force": float(info["linear_force"]), "gradient_max": info["gradient_max"], "slope": slope, "displacement": s,
               "converged": bool(done), "evaluations": total, "rounds": tried}


# ------------------------------------------------------------- curvature and a metric
def hessian_columns(model: MixedModel, x, columns, step: float = 1e-3, progress=None) -> np.ndarray:
    """``H[:, columns]`` of the energy at ``x`` in the scaled variables of :meth:`MixedModel.scale`
    (rows and columns both): central differences of the analytic gradient, two evaluations per
    column, ``step`` in scaled units.  ``progress(done, of)`` is called after each column."""
    x = np.asarray(x, dtype=float).ravel()
    sc = model.scale(x)
    columns = [int(j) for j in np.atleast_1d(columns)]
    out = np.zeros((model.n_var, len(columns)))
    for k, j in enumerate(columns):
        xp, xm = x.copy(), x.copy()
        xp[j] += step * sc[j]
        xm[j] -= step * sc[j]
        out[:, k] = (model.value_and_grad(xp)[1] - model.value_and_grad(xm)[1]) / (2.0 * step) * sc
        if progress is not None:
            progress(k + 1, len(columns))
    return out


def tiled_hessian(model: MixedModel, x, nx: int, ny: int, step: float = 1e-3, progress=None) -> np.ndarray:
    """The whole Hessian (scaled variables) at a state that is the ``nx x ny`` tiling of a smaller
    cell, chains in the order of :func:`tile_state`: the columns of the first cell's chains and of
    the four cell variables by :func:`hessian_columns`, every other chain's by translation (the
    column of a chain is the column of its translate with the rows moved along).  For a
    sixteen-chain tiling of a two-chain cell that is 92 columns instead of 708.  Symmetrised."""
    x = np.asarray(x, dtype=float).ravel()
    nc, n = model.nc, model.n_var
    n0 = nc // (nx * ny)
    if n0 * nx * ny != nc:
        raise ValueError(f"{nc} chains are not an {nx} x {ny} tiling")
    slots = lambda s: np.concatenate([model.index("place", [s]), model.index("internal", [s])])  # noqa: E731
    chain = lambda ix, iy, k: ((ix % nx) * ny + (iy % ny)) * n0 + k  # noqa: E731
    first = np.concatenate([np.arange(4)] + [slots(k) for k in range(n0)])
    H = np.zeros((n, n))
    H[:, first] = hessian_columns(model, x, first, step, progress)
    for dx in range(nx):
        for dy in range(ny):
            if dx == 0 and dy == 0:
                continue
            moved = np.arange(n)  # where each row goes under the translation by (dx, dy)
            for ix in range(nx):
                for iy in range(ny):
                    for k in range(n0):
                        moved[slots(chain(ix, iy, k))] = slots(chain(ix + dx, iy + dy, k))
            for k in range(n0):
                H[np.ix_(moved, slots(chain(dx, dy, k)))] = H[:, slots(k)]
    return 0.5 * (H + H.T)


def closure_multipliers(model: MixedModel, x, free, closed) -> np.ndarray:
    """``mu`` (chains of ``closed`` x 4) with ``gradient + sum_k mu_k dh_k/dy = 0`` in the least
    squares over the ``free`` variables: the closure's multipliers at a constrained minimum, the
    tension and the torques each chain's periodicity holds."""
    x = np.asarray(x, dtype=float).ravel()
    free = np.unique(np.asarray(free, dtype=int))
    closed = [int(s) for s in closed]
    sc = model.scale(x)[free]
    g = model.value_and_grad(x)[1][free] * sc
    J = model.closure(x, closed)[1][:, free] * sc
    return -np.linalg.lstsq(J.T, g, rcond=None)[0].reshape(len(closed), 4)


def closure_curvature(model: MixedModel, x, closed, multipliers, step: float = 1e-4) -> np.ndarray:
    """``sum_k mu_k d2h_k/dy dy`` (scaled variables, the full ``n_var x n_var`` matrix): what the
    closure adds to the Hessian of the energy to make the Hessian of the Lagrangian, whose values
    along the closure surface are the curvatures of the energy on it.  At the kinked sheet the
    energy's own Hessian has negative values along the surface (-3.1 the lowest in the four-chain
    cell) where the energy on the surface rises (+2.0 by second differences): the chains are under
    tension, and the tension's part is not small.  With it the two agree to three figures.
    Differences of the analytic closure Jacobian, which costs no energy evaluation."""
    x = np.asarray(x, dtype=float).ravel()
    sc = model.scale(x)
    out = np.zeros((model.n_var, model.n_var))
    mu = np.asarray(multipliers, dtype=float).reshape(-1, 4)
    for q, s in enumerate([int(v) for v in closed]):
        for j in np.concatenate([model.index("internal", [s]), model.index("c")]):
            xp, xm = x.copy(), x.copy()
            xp[j] += step * sc[j]
            xm[j] -= step * sc[j]
            out[:, j] += (mu[q] @ (model.closure(xp, [s])[1] - model.closure(xm, [s])[1])) / (2.0 * step) * sc
    return 0.5 * (out + out.T)


def chain_blocks(model: MixedModel, x, chains, step: float = 1e-3, progress=None) -> dict:
    """``{chain: block}``: the Hessian of the energy in each listed chain's own 44-odd variables
    (its placement, then its torsions and backbone angles; scaled) with every other chain held,
    by :func:`hessian_columns`.  The block of a chain in a sheet hardly depends on the size of the
    cell around it (four chains against sixteen: 3.6 in elements of up to 1400), so the blocks of a
    large cell can be taken from a small one."""
    out = {}
    for s in [int(v) for v in np.atleast_1d(chains)]:
        own = np.concatenate([model.index("place", [s]), model.index("internal", [s])])
        out[s] = hessian_columns(model, x, own, step, progress)[own]
    return out


def chain_metric(model: MixedModel, x, free, closed, blocks: dict, floor: float = 0.1, across: float = 10.0) -> np.ndarray:
    """The metric of :meth:`MixedModel.relax_on_surface` from one block per chain.

    For each chain the curvature of the energy on its own closure surface (its block of
    ``blocks``, plus the closure's curvature with the multipliers at ``x``, projected onto the
    tangent space of its four closure conditions; ``across`` along their normals) and its inverse
    square root, eigenvalues under ``floor`` raised to it; the blocks side by side are ``T``, with
    ones for a free variable that belongs to no chain.  ``blocks`` holds the blocks of some chains,
    in the order :func:`chain_blocks` uses; chain ``s`` takes that of the chain ``s`` modulo their
    number, which for a tiling of the cell they were computed in is its own translate.

    What it does: a chain's stiff directions are single torsions, each of which swings the rest
    of the chain through its neighbours unless the chain's placement and the other torsions
    follow, and its soft ones are those concerted moves.  No diagonal scaling tells the two apart;
    the chain's own block does.  The constrained Hessian of the kinked sheet, sixteen chains, has
    eigenvalues from 0.031 to 421 in the scaled variables and from 0.018 to 2.1 in these, half of
    them within 0.13 of one; what is left is the coupling between chains.  (Blocks
    of the internal coordinates alone do nothing: 0.08 to 415.)  The sheet's blocks serve for a
    state with one chain's kink elsewhere: the four-chain first-saddle hyperplane takes 29
    evaluations with them and 26 with the moving chain's block recomputed at the start."""
    x = np.asarray(x, dtype=float).ravel()
    free = np.unique(np.asarray(free, dtype=int))
    closed = [int(s) for s in closed]
    sc = model.scale(x)
    pos = {int(v): i for i, v in enumerate(free)}
    keys = sorted(blocks)
    mu = closure_multipliers(model, x, free, closed) if closed else np.zeros((0, 4))
    T = np.eye(free.size)
    for s in range(model.nc):
        own = np.concatenate([model.index("place", [s]), model.index("internal", [s])])
        keep = np.array([i for i, v in enumerate(own) if int(v) in pos], dtype=int)
        if keep.size == 0:
            continue
        Hs = np.asarray(blocks[keys[s % len(keys)]], dtype=float)
        if Hs.shape != (own.size, own.size):
            raise ValueError(f"block of chain {keys[s % len(keys)]}: shape {Hs.shape} for {own.size} variables")
        N = np.zeros((keep.size, 0))
        if s in closed:
            Hs = Hs + closure_curvature(model, x, [s], mu[closed.index(s)])[np.ix_(own, own)]
            Js = model.closure(x, [s])[1][:, own[keep]] * sc[own[keep]]
            _, sv, Vt = np.linalg.svd(Js, full_matrices=False)
            N = Vt[sv > 1e-10 * max(float(sv.max()), 1e-300)].T
        P = np.eye(keep.size) - N @ N.T
        M = P @ Hs[np.ix_(keep, keep)] @ P + across * (N @ N.T)
        w, V = np.linalg.eigh(0.5 * (M + M.T))
        loc = np.array([pos[int(own[i])] for i in keep])
        T[np.ix_(loc, loc)] = (V / np.sqrt(np.maximum(w, floor))) @ V.T
    return T


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
    xp = model.batched_module()  # the form the model's own evaluations take: the loops, or batched on a module
    E0, _, _ = placed_energy_and_grad(pk, frames, P0, latn, valences=vals, xp=xp)
    h = 1e-6

    def eg(t):
        dX = PD.group_displacements(groups, t, len(atoms))
        E, g, _ = placed_energy_and_grad(pk, frames, P0 + (J @ dX).reshape(-1, 3), latn, valences=vals, xp=xp)
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
        constraint=None, constraint_tol: float = 1e-8, climb_force: float | None = None, snapshot=None,
        metric=None) -> NEBResult:
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

    The band is good for finding where the top is and poor at settling it: where the energy has a
    soft direction across the path the climbing image creeps along it for hundreds of iterations at
    forces already under the criterion (the second top of the four-chain single-chain step: 12.66
    at a force of 0.16, 12.51 after 320 more iterations at 0.11, still falling).  A limited-memory
    quasi-Newton step on the band force was tried for that and is no better (forces of 0.15 to 0.8
    for 130 iterations, the images bunching towards one end).  :func:`saddle_polish` takes the
    climbing image and the band's tangent there and finds the stationary point itself.

    ``callback(iteration, energies, force_max, climbing)`` is called once per iteration, and
    ``snapshot(iteration, images)`` with the whole band as it stands (end points included): what a
    caller writes to disk to resume a band that runs for hours.

    ``metric`` is a matrix ``T`` (free x free), as :meth:`MixedModel.relax_on_surface` takes: the
    band then lives in ``w`` with ``y = T w``, ``y`` being the scaled free variables.  Tangents,
    springs, the projections onto the constraint surface and FIRE's velocities are all taken in
    ``w``; with ``T`` an inverse square root of the Hessian's chain blocks
    (:func:`chain_metric`) a chain's stiff single-torsion directions and its soft concerted ones
    have curvatures of order one alike, which is what lets FIRE's one time step serve both (the
    soft ones are where a band in per-variable units creeps for hundreds of iterations).  The
    force the band reports and converges on is carried back to the scaled variables
    (``T^-T F``), and ``max_step`` still limits a move in them, so both keep their meaning.
    Without ``metric`` nothing changes.
    """
    X = np.array(images, dtype=float)
    M, nv = X.shape
    free = np.arange(nv) if free is None else np.unique(np.asarray(free, dtype=int))
    sc = np.ones(nv) if scale is None else np.asarray(scale, dtype=float)
    scf = sc[free]
    T = None if metric is None else np.asarray(metric, dtype=float)
    if T is not None and T.shape != (len(free), len(free)):
        raise ValueError(f"metric: shape {T.shape} for {len(free)} free variables")
    Tinv = None if T is None else np.linalg.inv(T)
    Y = X[:, free] / scf  # scaled, free; in the metric's variables w (y = T w) when there is one
    if T is not None:
        Y = Y @ Tinv.T

    def scaled(w):  # the scaled free variables of a point of the band
        return w if T is None else T @ w

    def full(i):
        x = X[i].copy()
        x[free] = scaled(Y[i]) * scf
        return x

    def jacobian(J):  # the constraints' Jacobian in the band's variables
        Jf = np.asarray(J, dtype=float)[:, free] * scf
        return Jf if T is None else Jf @ T

    def normals(i):
        """Orthonormal columns spanning the constraint normals of image ``i`` in the band's
        variables, and the residual."""
        h, J = constraint(full(i))
        Jf = jacobian(J)
        _, sv, Vt = np.linalg.svd(Jf, full_matrices=False)
        keep = sv > 1e-10 * max(float(sv.max()), 1e-300)
        return Vt[keep].T, np.asarray(h, dtype=float), Jf

    def retract(i):
        for _ in range(12):
            h, J = constraint(full(i))
            h = np.asarray(h, dtype=float)
            if float(np.abs(h).max()) < constraint_tol:
                return
            Jf = jacobian(J)
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
        Xm[:, free] = (Y[1:-1] if T is None else Y[1:-1] @ T.T) * scf
        e, g = evaluate(Xm)
        E[1:-1] = np.asarray(e, dtype=float)
        G = np.asarray(g, dtype=float)[:, free] * scf  # scaled gradient of the moving images
        if T is not None:
            G = G @ T  # the gradient in the metric's variables (T^T g)
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
        fm = float(np.abs(Fn if T is None else Fn @ Tinv).max())  # in the scaled variables, with or without a metric
        history.append({"iteration": it, "force_max": fm, "barrier": float(E.max() - E[0]), "climbing": climbing})
        if callback is not None:
            callback(it, E.copy(), fm, climbing)
        if snapshot is not None:
            Xs = X.copy()
            if T is None:
                Xs[:, free] = Y * scf
            else:  # the end points as they were given, not through the metric and back
                Xs[1:-1, free] = (Y[1:-1] @ T.T) * scf
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
        big = float(np.abs(step if T is None else step @ T.T).max())  # the longest move of a scaled variable
        if big > max_step:
            step *= max_step / big
        Y[1:-1] += step
        if constraint is not None:
            for i in range(1, M - 1):
                retract(i)
    if T is None:
        X[:, free] = Y * scf
    else:
        X[1:-1, free] = (Y[1:-1] @ T.T) * scf
    return NEBResult(images=X, energies=E, force_max=fm, climbing=climbing, iterations=it, converged=converged, history=history)
