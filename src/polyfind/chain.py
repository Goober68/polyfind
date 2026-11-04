"""Torsion sequence -> Cartesian coordinates (rigid bond geometry).

Backbone atoms are placed with the NeRF (natural extension reference frame)
recursion; substituents are placed in the plane bisecting the backbone angle.
All routines accept a leading batch dimension so that thousands of chains can
be built at once on the active array backend (see :mod:`polyfind.backend`).

Conventions
-----------
* Backbone atom ``k`` has type ``k mod B`` (B = bonds per repeat); for PVDF
  even atoms are CH2 and odd atoms are CF2.
* RIS bond ``j`` is the bond between backbone atoms ``j+1`` and ``j+2``; its
  dihedral is defined by atoms ``j, j+1, j+2, j+3``.  A chain with ``N`` RIS
  bonds therefore has ``N+3`` backbone atoms.
* Consecutive RIS bonds ``j, j+1`` share atom ``j+2``; so pair type ``b``
  in :class:`polyfind.ris.RISModel` is centred on a backbone atom of type
  ``(b + 2) mod B``.
* A backbone atom's two pendants may differ (CFCl, CHCl); each then has its own
  element, bond length and charge, and pendant 1 always goes on the ``+w`` side of
  the local frame, which makes every such chain isotactic.  See
  :func:`substituent_positions` and :attr:`polyfind.polymers.Polymer.is_chiral`.
* A pendant is either a single atom or a small rigid fragment (nitrile, methoxy);
  the fragment's atoms are placed at fixed coordinates in the pendant frame built by
  :func:`pendant_frames`, of which the single atom is the degenerate case -- it sits
  on the frame's first axis at the bond length, which is exactly where
  :func:`substituent_positions` has always put it.  Atom order within a backbone atom
  is: the backbone atom, then pendant 1's atoms in order, then pendant 2's.

One thing elsewhere has not caught up with multi-atom pendants
-------------------------------------------------------------
``pack.MASS`` has no entry for N or O, so ``PeriodicChain.mass`` (and the densities
derived from it) raises for the phase-3 chemistries.  Everything that builds
coordinates -- :func:`build_chain`, :func:`build_chain_batch`, and the refinement's
``pack._batch_block_coords`` / ``linegroup._block_coords``, which are
:func:`build_oligomer_coords` -- places pendants through the same
:func:`pendant_frames` and is general.

Performance notes
-----------------
The builders are dominated by per-call overhead on small arrays, not arithmetic.
:func:`build_backbone` therefore runs the NeRF recursion on component arrays
``X, Y, Z`` of shape ``(N+3, M)`` with the per-atom-type scalars hoisted and the
dihedral trigonometry done for all steps at once (and in Python floats for a single
chain on NumPy); :func:`place_pendants` builds the pendant frames of every backbone
atom in one call and places each pendant atom of each atom type with one broadcast.
Both are bit for bit the per-step / per-atom forms (:func:`nerf`,
:func:`pendant_positions`), which stay as the reference implementations.
"""
from __future__ import annotations

import math
import warnings
from dataclasses import dataclass, field

import numpy as np

from . import backend as bk
from .polymers import Polymer, pendant_pair
from .chemical_graph import ChemicalBondGraph


def _unit(v, xp):
    return v / xp.linalg.norm(v, axis=-1, keepdims=True)


def _cross(u, v, xp):
    """Component cross product of (..., 3) arrays.

    The same arithmetic, in the same order, as ``xp.cross`` for 3-vectors, so the
    result is bit-for-bit the same on NumPy -- minus the ~20x call overhead of
    ``xp.cross`` on small arrays, which is what the per-step NeRF recursion pays.
    """
    ux, uy, uz = u[..., 0], u[..., 1], u[..., 2]
    vx, vy, vz = v[..., 0], v[..., 1], v[..., 2]
    return xp.stack([uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx], axis=-1)


def nerf(a, b, c, bond, angle_deg, dihedral_deg, xp=np):
    """Place atom d given a, b, c (..., 3), |cd| = bond, angle bcd, dihedral abcd."""
    ang = xp.deg2rad(xp.asarray(angle_deg, dtype=float))
    dih = xp.deg2rad(xp.asarray(dihedral_deg, dtype=float))
    bond = xp.asarray(bond, dtype=float)
    ang, dih, bond = xp.broadcast_arrays(ang, dih, bond)
    bc = _unit(c - b, xp)
    n = _unit(xp.cross(b - a, bc), xp)
    m = xp.cross(n, bc)
    # local coordinates of d in the frame (bc, m, n) centred on c
    d2 = xp.stack([-bond * xp.cos(ang), bond * xp.sin(ang) * xp.cos(dih), bond * xp.sin(ang) * xp.sin(dih)], axis=-1)
    return c + d2[..., 0:1] * bc + d2[..., 1:2] * m + d2[..., 2:3] * n


def backbone_angle_lookup(polymer: Polymer, bond_angles=None, xp=np):
    """``k -> backbone angle at backbone atom k`` for an optional per-repeat override.

    ``bond_angles`` is ``None`` (use the polymer's own frozen angles), an array of
    shape ``(R,)``, or a batch of shape ``(M, R)``; entry ``r`` is the angle of every
    backbone atom ``k`` with ``k % R == r``, i.e. the override is broadcast along the
    chain exactly as the polymer's own repeat is.  ``R`` need not equal
    ``polymer.bonds_per_repeat``: a crystallographic repeat of several monomers can
    carry one angle per backbone atom of the repeat.

    The returned callable gives a Python float when no override is in force (so the
    frozen :class:`~polyfind.polymers.Polymer` path is bit-for-bit unchanged) and an
    array of shape ``(M,)`` or ``(1,)`` otherwise, which every consumer here
    broadcasts against the batch dimension.
    """
    B = polymer.bonds_per_repeat
    if bond_angles is None:
        return lambda k: polymer.backbone[k % B].backbone_angle
    ba = xp.asarray(bond_angles, dtype=float)
    if ba.ndim == 1:
        ba = ba[None, :]
    if ba.ndim != 2 or ba.shape[1] == 0:
        raise ValueError("bond_angles must have shape (R,) or (M, R) with R >= 1")
    R = ba.shape[1]
    return lambda k: ba[:, k % R]


def _backbone_scalar(N: int, L: float, cos_c, sin_c, P, Q) -> np.ndarray:
    """The NeRF recursion of :func:`build_backbone` for one chain in Python floats.

    Same formulas and operation order as the array loop there; ``cos_c[k]``,
    ``sin_c[k]`` are ``-L cos``, ``L sin`` of the backbone angle at atom ``k`` and
    ``P[j]``, ``Q[j]`` the dihedral coefficients of RIS bond ``j`` (all plain floats).
    Returns ``(N+3, 3)``.  Bit-identical to the array loop (IEEE doubles either way),
    and ~20x cheaper than NumPy dispatch on 1-element arrays.
    """
    sqrt = math.sqrt
    n = N + 3
    xs, ys, zs = [0.0] * n, [0.0] * n, [0.0] * n
    xs[1] = L
    xs[2] = L + cos_c[1]  # L - L cos(theta_1)
    ys[2] = sin_c[1]
    for k in range(3, n):
        ax, ay, az = xs[k - 3], ys[k - 3], zs[k - 3]
        bx, by, bz = xs[k - 2], ys[k - 2], zs[k - 2]
        cx, cy, cz = xs[k - 1], ys[k - 1], zs[k - 1]
        cc, p, q = cos_c[k - 1], P[k - 3], Q[k - 3]
        vx, vy, vz = cx - bx, cy - by, cz - bz
        nv = sqrt(vx * vx + vy * vy + vz * vz)
        vx, vy, vz = vx / nv, vy / nv, vz / nv  # bc
        ux, uy, uz = bx - ax, by - ay, bz - az
        nx, ny, nz = uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx
        nn = sqrt(nx * nx + ny * ny + nz * nz)
        nx, ny, nz = nx / nn, ny / nn, nz / nn  # n
        mx, my, mz = ny * vz - nz * vy, nz * vx - nx * vz, nx * vy - ny * vx  # m = n x bc
        xs[k] = cx + cc * vx + p * mx + q * nx
        ys[k] = cy + cc * vy + p * my + q * ny
        zs[k] = cz + cc * vz + p * mz + q * nz
    return np.array([xs, ys, zs]).T.copy()


def build_backbone(polymer: Polymer, dihedrals_deg, xp=None, bond_angles=None):
    """Backbone coordinates for a dihedral array of shape (N,) or (M, N).

    Returns (N+3, 3) or (M, N+3, 3).  ``bond_angles`` optionally overrides the
    polymer's frozen backbone angles; see :func:`backbone_angle_lookup`.

    The recursion is the NeRF step of :func:`nerf` -- atom ``k`` from atoms
    ``k-3, k-2, k-1``, the bond length, the backbone angle at ``k-1`` and RIS dihedral
    ``k-3`` -- but written out on component arrays ``X, Y, Z`` of shape ``(N+3, M)``,
    with ``bc = unit(c - b)``, ``n = unit((b - a) x bc)`` and ``m = n x bc``:

        d = c - L cos(theta) bc + L sin(theta) cos(phi) m + L sin(theta) sin(phi) n

    The per-atom-type scalars ``-L cos(theta)``, ``L sin(theta)`` are hoisted out of the
    loop and the dihedral trigonometry is done for all steps at once, so each step is
    ~35 element-wise operations on ``(M,)`` rows instead of ~15 NumPy calls on ``(M, 3)``
    blocks (``xp.cross`` alone costs ~20x an explicit component product).  Every
    operation is the one ``nerf`` performs, in its order (division by the norm, the
    same trigonometric ufuncs, the same left-to-right sum), so the coordinates are bit
    for bit those of the per-step form.  A single chain on NumPy runs the same recursion
    in Python floats (:func:`_backbone_scalar`), which is again bit-identical.
    """
    xp = xp or bk.get_backend()
    dih = xp.asarray(dihedrals_deg, dtype=float)
    batched = dih.ndim == 2
    if not batched:
        dih = dih[None, :]
    M, N = dih.shape
    n_atoms = N + 3
    L = float(polymer.bond_length)
    B = polymer.bonds_per_repeat
    # Backbone angle at every atom (the same ufuncs, on the same values, as ``nerf`` uses,
    # so the coordinates are bit for bit those of the per-step form): ``cos_c[k]`` is
    # ``-L cos(theta_k)`` (the coefficient of ``bc``) and ``sin_c[k]`` is ``L sin(theta_k)``.
    if bond_angles is None:
        # the polymer's own frozen angles: one value per atom type, as Python floats
        rad = xp.deg2rad(xp.asarray([polymer.backbone[k % B].backbone_angle for k in range(n_atoms)], dtype=float))
        cos_c = (-L * xp.cos(rad)).tolist()
        sin_c = (L * xp.sin(rad)).tolist()
        fac = xp.asarray(sin_c[2 : N + 2], dtype=float)[:, None]  # (N, 1): step k uses the angle at atom k-1
        scalar_ok = True
    else:
        ba = xp.asarray(bond_angles, dtype=float)
        if ba.ndim == 1:
            ba = ba[None, :]
        if ba.ndim != 2 or ba.shape[1] == 0:
            raise ValueError("bond_angles must have shape (R,) or (M, R) with R >= 1")
        Mb, R = ba.shape
        rad = xp.deg2rad(ba)
        idx = np.arange(n_atoms) % R
        cos_c = (-L * xp.cos(rad))[:, idx].T  # (n_atoms, Mb); Mb is M or 1, broadcast along the batch
        sin_c = (L * xp.sin(rad))[:, idx].T
        fac = sin_c[2 : N + 2]  # (N, Mb)
        scalar_ok = Mb == 1
        if M == 1 and xp is np:
            cos_c, sin_c = cos_c[:, 0].tolist(), sin_c[:, 0].tolist()
    # dihedral coefficients for every step at once: P multiplies (n x v)^, Q multiplies n^
    drad = xp.deg2rad(dih)
    P = xp.cos(drad).T * fac  # (N, M)
    Q = xp.sin(drad).T * fac
    if M == 1 and xp is np and scalar_ok:
        coords = _backbone_scalar(N, L, cos_c, sin_c, P[:, 0].tolist(), Q[:, 0].tolist())
        return coords[None] if batched else coords
    # component rows, atom-major: every step touches three contiguous (M,) rows
    X = xp.zeros((n_atoms, M), dtype=float)
    Y = xp.zeros((n_atoms, M), dtype=float)
    Z = xp.zeros((n_atoms, M), dtype=float)
    # first three atoms in the xy plane
    X[1] = L
    X[2] = L + cos_c[1]  # L - L cos(theta_1)
    Y[2] = sin_c[1]
    for k in range(3, n_atoms):
        ax, ay, az = X[k - 3], Y[k - 3], Z[k - 3]
        bx, by, bz = X[k - 2], Y[k - 2], Z[k - 2]
        cx, cy, cz = X[k - 1], Y[k - 1], Z[k - 1]
        cc, p, q = cos_c[k - 1], P[k - 3], Q[k - 3]
        vx, vy, vz = cx - bx, cy - by, cz - bz
        nv = xp.sqrt(vx * vx + vy * vy + vz * vz)
        vx, vy, vz = vx / nv, vy / nv, vz / nv  # bc
        ux, uy, uz = bx - ax, by - ay, bz - az
        nx, ny, nz = uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx
        nn = xp.sqrt(nx * nx + ny * ny + nz * nz)
        nx, ny, nz = nx / nn, ny / nn, nz / nn  # n
        mx, my, mz = ny * vz - nz * vy, nz * vx - nx * vz, nx * vy - ny * vx  # m = n x bc
        X[k] = cx + cc * vx + p * mx + q * nx
        Y[k] = cy + cc * vy + p * my + q * ny
        Z[k] = cz + cc * vz + p * mz + q * nz
    coords = xp.ascontiguousarray(xp.stack([X, Y, Z], axis=-1).transpose(1, 0, 2))  # (M, N+3, 3)
    return coords if batched else coords[0]


def substituent_positions(x_prev, x, x_next, bond, sub_angle_deg, xp=np):
    """Two substituents on atom x with backbone neighbours x_prev, x_next; shapes (..., 3).

    Both pendants lie in the plane through ``x`` spanned by the backbone-angle bisector
    ``u`` and the backbone-plane normal ``w`` (i.e. the plane perpendicular to the
    backbone plane), one on each side of ``u``, and the angle between them is exactly
    ``sub_angle_deg``.

    ``bond`` is either one length, shared by both pendants (the symmetric case, and
    what every existing caller passes as ``spec.sub_bond``), or a ``(first, second)``
    pair for a backbone atom bearing two different pendants -- so a caller that simply
    forwards ``spec.sub_bond`` stays correct for asymmetric monomers too.

    Angle-splitting convention: **equal split**, i.e. each pendant sits at
    ``sub_angle_deg / 2`` from the bisector, whatever its bond length.  This is
    defensible and is what the model can actually support: the pendants' bond lengths
    do not by themselves determine how the angle should be shared, and a physically
    motivated unequal split (VSEPR: the bulkier pendant pushes the smaller one, so the
    bisector tilts towards the smaller one) would need a further per-monomer parameter
    -- the two backbone-substituent angles separately -- that :class:`BackboneAtom`
    does not carry.  Any split preserves the specified inter-substituent angle, since
    the two directions are unit vectors at +/- half from ``u`` regardless of ``bond``;
    only the placement relative to the backbone changes.  The equal split also keeps
    the symmetric case bit-for-bit as it was.

    Handedness, and hence tacticity: pendant 1 always goes on the ``+w`` side, with
    ``w = unit((x_prev - x) x (x_next - x))`` fixed by the chain direction, so every
    stereocentre of a chain gets the same configuration -- the chain built is the
    isotactic one.  See :attr:`polyfind.polymers.Polymer.is_chiral`.

    This places one atom per pendant, which is all a pendant used to be.  For a pendant
    that is a rigid *fragment* it places the fragment's first atom -- the one bonded to
    the backbone -- and :func:`pendant_positions` places the rest, off the frame that
    :func:`pendant_frames` builds around these same two directions.
    """
    n1 = _unit(x_prev - x, xp)
    n2 = _unit(x_next - x, xp)
    u = -_unit(n1 + n2, xp)  # bisector, pointing away from the backbone
    w = _unit(xp.cross(n1, n2), xp)  # normal to the backbone plane
    b1, b2 = pendant_pair(bond, "bond")
    half = np.deg2rad(sub_angle_deg) / 2.0
    s1 = x + b1 * (np.cos(half) * u + np.sin(half) * w)
    s2 = x + b2 * (np.cos(half) * u - np.sin(half) * w)
    return s1, s2


def pendant_frames(x_prev, x, x_next, sub_angle_deg, xp=np):
    """Orthonormal frame ``(e1, e2, e3)`` per pendant of atom x; shapes (..., 3).

    Returns ``((e1, e2, e3), (e1, e2, e3))``, one triple per pendant, in which a
    fragment's local coordinates are read: an atom at local ``(lx, ly, lz)`` sits at
    ``x + lx e1 + ly e2 + lz e3``.  The axes are

    * ``e1``: the pendant bond direction, ``cos(a/2) u +- sin(a/2) w`` with ``u`` the
      backbone-angle bisector and ``w`` the backbone-plane normal, exactly the direction
      :func:`substituent_positions` uses -- so an atom at ``(bond, 0, 0)`` lands
      bit-for-bit where a single-atom substituent always did;
    * ``e2``: perpendicular to ``e1`` in the ``(u, w)`` plane, on the ``u`` side;
    * ``e3``: ``+e1 x e2`` for pendant 1 and ``-e1 x e2`` for pendant 2.

    That sign on ``e3`` is the one real choice here, and it is made so that the two
    pendant frames are **mirror images of one another** through the local backbone plane
    (the plane through ``x`` containing its two backbone neighbours).  The consequence is
    that two *identical* fragments leave the backbone atom locally mirror-symmetric,
    exactly as two identical single atoms do -- so VDCN, with two nitriles, is achiral
    for the same reason PVDF is, and
    :attr:`polyfind.polymers.BackboneAtom.is_stereocentre` can go on
    being a plain comparison of the two pendants.  With the other sign, two identical
    fragments would be placed as a *rotation* of one another and a symmetric monomer
    would come out spuriously chiral.

    The price is that pendant 2's frame is left-handed, so a fragment placed there is the
    *reflection* of its local specification.  That is invisible for a fragment with a
    mirror plane of its own, which nitrile and methoxy both have, but it means a genuinely
    chiral fragment would appear as its enantiomer on pendant 2 -- and two identical
    chiral fragments would be placed as a meso pair rather than as two of the same hand.
    Getting that case right needs the fragment to declare its own handedness; nothing
    here does yet.
    """
    n1 = _unit(x_prev - x, xp)
    n2 = _unit(x_next - x, xp)
    u = -_unit(n1 + n2, xp)
    w = _unit(xp.cross(n1, n2), xp)
    # ``sub_angle_deg`` is one angle, or an array broadcastable against the leading
    # shape (one per backbone atom, as :func:`place_pendants` passes it)
    half = np.deg2rad(sub_angle_deg) / 2.0
    c, s = np.cos(half), np.sin(half)
    if xp is not np:
        c, s = xp.asarray(c), xp.asarray(s)
    e1_a, e2_a = c * u + s * w, s * u - c * w
    e1_b, e2_b = c * u - s * w, s * u + c * w
    return (e1_a, e2_a, xp.cross(e1_a, e2_a)), (e1_b, e2_b, -xp.cross(e1_b, e2_b))


def pendant_positions(x_prev, x, x_next, spec, xp=np):
    """Positions of every pendant atom of backbone atom ``x``; shapes (..., 3).

    Returns one list of positions per pendant, in the pendant's own atom order, from the
    frames of :func:`pendant_frames`.  A single-atom pendant gives a one-element list
    holding exactly ``substituent_positions``' answer, bit for bit: its only atom is at
    local ``(bond, 0, 0)`` and the zero components are skipped rather than added.
    """
    frames = pendant_frames(x_prev, x, x_next, spec.sub_angle, xp=xp)
    out = []
    for pendant, (e1, e2, e3) in zip(spec.pendants, frames):
        positions = []
        for atom in pendant.atoms:
            lx, ly, lz = atom.offset
            p = x + lx * e1
            if ly:
                p = p + ly * e2
            if lz:
                p = p + lz * e3
            positions.append(p)
        out.append(positions)
    return tuple(out)


def extended_backbone(polymer: Polymer, dihedrals_deg, bond_angles=None, xp=np, backbone=None):
    """Backbone of ``(M, N)`` dihedrals plus one virtual trans atom beyond each end.

    Returns ``(M, N+5, 3)``; real backbone atom ``k`` sits at ``[:, k+1]``.  The virtual
    atoms are what the terminal backbone atoms' pendants are placed from, exactly as
    :func:`build_chain` does it.

    ``backbone``, if given, is the ``(M, N+3, 3)`` result of :func:`build_backbone` for
    these very dihedrals and angles, and is used instead of building it again (the
    analytic Jacobian of :mod:`polyfind.chainjac` differentiates the backbone it has
    already built, and places the pendants on it through here, so its coordinates are
    those of :func:`build_oligomer_coords` by construction).
    """
    dih = xp.asarray(dihedrals_deg, dtype=float)
    M, N = dih.shape
    L = polymer.bond_length
    ang_of = backbone_angle_lookup(polymer, bond_angles, xp=xp)
    bb = build_backbone(polymer, dih, xp=xp, bond_angles=bond_angles) if backbone is None else backbone  # (M, N+3, 3)
    v0 = nerf(bb[:, 2], bb[:, 1], bb[:, 0], L, ang_of(0), 180.0, xp=xp)
    v1 = nerf(bb[:, N], bb[:, N + 1], bb[:, N + 2], L, ang_of(N + 2), 180.0, xp=xp)
    return xp.concatenate([xp.broadcast_to(v0, (M, 3))[:, None], bb, xp.broadcast_to(v1, (M, 3))[:, None]], axis=1)


def atom_starts(polymer: Polymer, n_backbone: int) -> np.ndarray:
    """``starts[k]`` = index of backbone atom ``k`` in :func:`build_chain`'s atom order.

    Backbone atom ``k`` is followed by its pendant atoms (pendant 1's, then pendant 2's),
    so ``starts[k+1] - starts[k] = 1 + n_pendant_atoms``; ``starts[-1]`` is the number
    of atoms before any end cap.
    """
    B = polymer.bonds_per_repeat
    widths = [1 + polymer.backbone[k % B].n_pendant_atoms for k in range(n_backbone)]
    return np.concatenate([[0], np.cumsum(widths)]).astype(int)


def place_pendants(polymer: Polymer, bb_ext, xp=np):
    """All-atom coordinates ``(M, n, 3)`` in :func:`build_chain` order from an extended backbone.

    ``bb_ext`` is ``(M, n_backbone + 2, 3)`` as returned by :func:`extended_backbone`.
    Every backbone atom's pendant frames are built in **one** :func:`pendant_frames` call
    over the whole ``(M, n_backbone)`` block (the pendant angle is broadcast per atom),
    and then every pendant atom of every backbone-atom *type* is placed with one
    broadcast: the gather of the frames of that type's atoms, ``x + lx e1 [+ ly e2]
    [+ lz e3]``, and one scatter into the output.  The arithmetic per atom is the same
    as :func:`pendant_positions` on the same frames, so the coordinates are bit for bit
    those of the per-atom loop (the zero components of a single-atom pendant are skipped
    rather than added, as there).  For PVDF that is four placements for a 43-atom
    backbone instead of 43 calls of the per-atom routine.
    """
    M, n_ext, _ = bb_ext.shape
    n_bb = n_ext - 2
    B = polymer.bonds_per_repeat
    specs = [polymer.backbone[k % B] for k in range(n_bb)]
    starts = atom_starts(polymer, n_bb)
    out = xp.empty((M, int(starts[-1]), 3), dtype=float)
    x_prev, x, x_next = bb_ext[:, :-2], bb_ext[:, 1:-1], bb_ext[:, 2:]
    out[:, starts[:-1]] = x
    sub_angle = np.array([s.sub_angle for s in specs], dtype=float)[:, None]  # (n_bb, 1), per atom
    frames = pendant_frames(x_prev, x, x_next, sub_angle, xp=xp)  # each axis (M, n_bb, 3)
    for r in range(min(B, n_bb)):  # backbone-atom type: same pendants, same local offsets
        ks = np.arange(r, n_bb, B)
        xr = x[:, ks]
        at = starts[ks] + 1
        for pendant, (e1, e2, e3) in zip(specs[r].pendants, frames):
            e1r, e2r, e3r = e1[:, ks], None, None
            for atom in pendant.atoms:
                lx, ly, lz = atom.offset
                p = xr + lx * e1r
                if ly:
                    e2r = e2[:, ks] if e2r is None else e2r
                    p = p + ly * e2r
                if lz:
                    e3r = e3[:, ks] if e3r is None else e3r
                    p = p + lz * e3r
                out[:, at] = p
                at = at + 1
    return out


def build_oligomer_coords(polymer: Polymer, dihedrals_deg, bond_angles=None, xp=np):
    """``(M, N)`` dihedrals -> ``(M, n, 3)`` all-atom coordinates, ``build_chain(cap=False)`` order.

    :func:`extended_backbone` followed by :func:`place_pendants`: the whole oligomer in a
    handful of vectorised calls, which is what the refinement's chain builders
    (:func:`polyfind.pack._batch_block_coords`, :func:`polyfind.linegroup._block_coords`)
    and :func:`build_chain_batch` are built on.
    """
    return place_pendants(polymer, extended_backbone(polymer, dihedrals_deg, bond_angles, xp=xp), xp=xp)


_CHIRAL_WARNED: set[str] = set()


def warn_if_chiral(polymer: Polymer) -> None:
    """Warn once per polymer that a chain with stereocentres has a tacticity.

    Emitted from :func:`build_chain`, which every path into the funnel goes through
    (the batched builder, ``pack.periodic_chain``, ``pack._template`` and hence the
    packing search all build a template with it), so a user who fits or packs one of
    these chemistries is told before they get numbers back.  What the tacticity costs
    the user is the choice of enantiomer, which the model cannot express; the places
    that used to assume an achiral repeat (``fit_ris``'s ``symmetrize``,
    ``enumerate_periodic``'s deduplication and ``linegroup``'s glide) now take it from
    :attr:`~polyfind.polymers.Polymer.is_chiral` themselves.
    """
    if not polymer.is_chiral or polymer.name in _CHIRAL_WARNED:
        return
    _CHIRAL_WARNED.add(polymer.name)
    warnings.warn(
        f"polymer {polymer.name!r} has backbone atoms with two different substituents, i.e. "
        "stereocentres. The chain built here is the ISOTACTIC one (pendant 1 on a fixed side "
        "of the local frame at every backbone atom); syndiotactic and atactic chains are not "
        "expressible in this model. Because the repeat is chiral, G+ and G- conformers are no "
        "longer mirror-equivalent: they are the two handednesses of a one-handed helix and "
        "have genuinely different energies. The relation that does hold is mirror composed "
        "with chain reversal, E(phi_1..phi_N) == E(-phi_N..-phi_1). fit_ris (symmetrize="
        "'auto'), enumerate_periodic and linegroup all take that from Polymer.is_chiral, so "
        "forcing symmetrize=True is the one way left to average the asymmetry away. "
        "See Polymer.is_chiral and docs/CHEMISTRY_EXTENSION.md.",
        UserWarning,
        stacklevel=3,
    )


@dataclass
class Structure:
    """All-atom oligomer with fixed bond geometry; only dihedrals vary."""

    polymer: Polymer
    elements: list[str]
    coords: np.ndarray  # (n, 3)
    charges: np.ndarray  # (n,)
    chemical_graph: ChemicalBondGraph
    backbone: np.ndarray  # indices of the real backbone atoms in chain order
    n_dihedrals: int
    dihedrals: np.ndarray  # (N,) the RIS dihedrals used to build it
    # backbone atom index -> its pendant atom indices: pendant 1's atoms then pendant 2's,
    # each in the pendant's own order (one index per pendant unless a pendant is a fragment)
    subs_of: dict = field(default_factory=dict)

    def __post_init__(self):
        if self.chemical_graph.atom_count!=len(self.elements):
            raise ValueError('structure and chemical graph atom counts differ')

    @property
    def bonds(self):
        return self.chemical_graph.bonds

    def dihedral_atoms(self, j: int) -> tuple[int, int, int, int]:
        bb = self.backbone
        return int(bb[j]), int(bb[j + 1]), int(bb[j + 2]), int(bb[j + 3])

    @property
    def n_atoms(self) -> int:
        return len(self.elements)

    def to_xyz(self, comment: str = "") -> str:
        lines = [str(self.n_atoms), comment]
        for e, (x, y, z) in zip(self.elements, self.coords):
            lines.append(f"{e} {x:.6f} {y:.6f} {z:.6f}")
        return "\n".join(lines) + "\n"


def build_chain(polymer: Polymer, dihedrals_deg, cap: bool = True, bond_angles=None) -> Structure:
    """All-atom oligomer for ``N`` RIS dihedrals (N+3 backbone atoms), H-capped by default.

    The chain is built with one virtual backbone atom beyond each end (trans) so
    that terminal substituents are placed consistently; with ``cap=True`` an H atom
    is placed along each virtual bond, i.e. CH3-/CH2F-type end groups.

    ``bond_angles`` (shape ``(R,)``, see :func:`backbone_angle_lookup`) overrides the
    polymer's frozen backbone angles; the substituents are placed from the *actual*
    backbone neighbours, so they follow the changed geometry automatically.

    For a polymer with stereocentres the chain built is the isotactic one; a warning
    says so once (see :func:`warn_if_chiral`).
    """
    warn_if_chiral(polymer)
    dih = np.asarray(dihedrals_deg, dtype=float)
    N = len(dih)
    B = polymer.bonds_per_repeat
    L = polymer.bond_length
    ang_of = backbone_angle_lookup(polymer, bond_angles, xp=np)
    bb = build_backbone(polymer, dih, xp=np, bond_angles=bond_angles)  # (N+3, 3) real atoms, canonical frame
    # virtual (trans) backbone atoms beyond each end, in the *same* frame as ``bb`` so that
    # coordinates agree with build_backbone (helix analysis relies on this)
    v0 = nerf(bb[2], bb[1], bb[0], L, ang_of(0), 180.0, xp=np)
    v1 = nerf(bb[N], bb[N + 1], bb[N + 2], L, ang_of(N + 2), 180.0, xp=np)
    bb_ext = np.vstack([np.reshape(v0, 3)[None], bb, np.reshape(v1, 3)[None]])  # (N+5, 3); real atom k sits at bb_ext[k+1]
    elements, coords, charges, bonds = [], [], [], []
    backbone_idx, subs_of = [], {}
    for k in range(N + 3):
        spec = polymer.backbone[k % B]
        x = bb_ext[k + 1]
        idx = len(elements)
        elements.append(spec.element)
        coords.append(x)
        charges.append(spec.charge)
        backbone_idx.append(idx)
        if k > 0:
            bonds.append((backbone_idx[k - 1], idx,1))
        groups = pendant_positions(bb_ext[k], x, bb_ext[k + 2], spec, xp=np)
        subs = []
        for pendant, positions in zip(spec.pendants, groups):
            first = len(elements)
            for atom, pos in zip(pendant.atoms, positions):
                elements.append(atom.element)
                coords.append(pos)
                charges.append(atom.charge)
                subs.append(len(elements) - 1)
            bonds.append((idx,first,1))
            bonds.extend((first+i,first+j,order) for i,j,order in pendant.chemical_graph.edges)
        subs_of[idx] = subs
    if cap:
        for k, virt in ((0, bb_ext[0]), (N + 2, bb_ext[N + 4])):
            idx = backbone_idx[k]
            x = coords[idx]
            v = virt - x
            pos = x + 1.09 * v / np.linalg.norm(v)
            elements.append("H")
            coords.append(pos)
            charges.append(-float(np.sum(charges[idx:idx + 1])) * 0.0)  # neutral cap
            bonds.append((idx,len(elements)-1,1))
    return Structure(
        polymer=polymer,
        elements=elements,
        coords=np.array(coords),
        charges=np.array(charges),
        chemical_graph=ChemicalBondGraph(len(elements),tuple(bonds)),
        backbone=np.array(backbone_idx),
        n_dihedrals=N,
        dihedrals=dih,
        subs_of=subs_of,
    )


def build_chain_batch(polymer: Polymer, dihedrals_deg, cap: bool = True, bond_angles=None) -> tuple[Structure, np.ndarray]:
    """Batched version of :func:`build_chain`.

    ``dihedrals_deg`` has shape ``(M, N)``.  Returns ``(template, coords)`` where
    ``template`` is ``build_chain(polymer, dihedrals_deg[0], cap)`` (topology,
    elements, bonds -- identical for every row since only the dihedrals vary)
    and ``coords`` is an ``(M, n_atoms, 3)`` array whose row ``m`` equals
    ``build_chain(polymer, dihedrals_deg[m], cap).coords`` to numerical precision.

    Implemented with the batched :func:`build_backbone`, the block-wide pendant
    placement of :func:`place_pendants` and :func:`nerf` for the virtual end atoms, so
    the whole oligomer (backbone, substituents, end caps) is placed in a handful of
    vectorised calls with no Python loop over conformers or over backbone atoms.

    ``bond_angles`` optionally overrides the polymer's frozen backbone angles, either
    once for the whole batch (shape ``(R,)``) or per row (shape ``(M, R)``), which is
    what makes a finite-difference gradient over the bond angles one batched call.
    """
    dih = np.asarray(dihedrals_deg, dtype=float)
    M, N = dih.shape
    row0 = None if bond_angles is None else np.atleast_2d(np.asarray(bond_angles, dtype=float))[0]
    template = build_chain(polymer, dih[0], cap, bond_angles=row0)
    bb_ext = extended_backbone(polymer, dih, bond_angles, xp=np)  # (M, N+5, 3)
    body = place_pendants(polymer, bb_ext, xp=np)  # backbone atom k at template.backbone[k], then its pendants
    if not cap:
        return template, body
    coords = np.empty((M, template.n_atoms, 3), dtype=float)
    n_body = body.shape[1]
    coords[:, :n_body] = body
    for i, (k, virt) in enumerate(((0, bb_ext[:, 0]), (N + 2, bb_ext[:, N + 4]))):
        x = bb_ext[:, k + 1]
        v = virt - x
        coords[:, n_body + i] = x + 1.09 * v / np.linalg.norm(v, axis=-1, keepdims=True)
    return template, coords


def set_dihedrals(struct: Structure, dihedrals_deg) -> Structure:
    """Rebuild the same oligomer with new dihedrals (cheap: rigid geometry)."""
    # backbone atoms plus every pendant atom; the only extras are the two end caps, so
    # anything beyond that count means the structure was capped.  Counted from the
    # structure itself rather than as ``3 * (N + 3)``, which held only while every
    # pendant was a single atom.
    uncapped = len(struct.backbone) + sum(len(v) for v in struct.subs_of.values())
    cap = len(struct.elements) > uncapped
    return build_chain(struct.polymer, dihedrals_deg, cap=cap)


# ------------------------------------------------------------------ measurements
def distance(c, i, j):
    return float(np.linalg.norm(c[i] - c[j]))


def angle(c, i, j, k):
    a = c[i] - c[j]
    b = c[k] - c[j]
    return float(np.degrees(np.arccos(np.clip(a @ b / np.linalg.norm(a) / np.linalg.norm(b), -1, 1))))


def dihedral(c, i, j, k, l):
    b0, b1, b2 = c[j] - c[i], c[k] - c[j], c[l] - c[k]
    n1, n2 = np.cross(b0, b1), np.cross(b1, b2)
    x = n1 @ n2
    y = np.linalg.norm(b1) * (b0 @ n2)
    return float(np.degrees(np.arctan2(y, x)))
