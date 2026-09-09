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
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from . import backend as bk
from .polymers import Polymer


def _unit(v, xp):
    return v / xp.linalg.norm(v, axis=-1, keepdims=True)


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


def build_backbone(polymer: Polymer, dihedrals_deg, xp=None):
    """Backbone coordinates for a dihedral array of shape (N,) or (M, N).

    Returns (N+3, 3) or (M, N+3, 3).
    """
    xp = xp or bk.get_backend()
    dih = xp.asarray(dihedrals_deg, dtype=float)
    batched = dih.ndim == 2
    if not batched:
        dih = dih[None, :]
    M, N = dih.shape
    B = polymer.bonds_per_repeat
    L = polymer.bond_length
    angles = [polymer.backbone[k % B].backbone_angle for k in range(N + 3)]
    coords = xp.zeros((M, N + 3, 3), dtype=float)
    # first three atoms in the xy plane
    coords[:, 1, 0] = L
    a1 = np.deg2rad(angles[1])
    coords[:, 2, 0] = L - L * np.cos(a1)
    coords[:, 2, 1] = L * np.sin(a1)
    for k in range(3, N + 3):
        coords[:, k] = nerf(coords[:, k - 3], coords[:, k - 2], coords[:, k - 1], L, angles[k - 1], dih[:, k - 3], xp=xp)
    return coords if batched else coords[0]


def substituent_positions(x_prev, x, x_next, bond, sub_angle_deg, xp=np):
    """Two substituents on atom x with backbone neighbours x_prev, x_next; shapes (..., 3)."""
    n1 = _unit(x_prev - x, xp)
    n2 = _unit(x_next - x, xp)
    u = -_unit(n1 + n2, xp)  # bisector, pointing away from the backbone
    w = _unit(xp.cross(n1, n2), xp)  # normal to the backbone plane
    half = np.deg2rad(sub_angle_deg) / 2.0
    s1 = x + bond * (np.cos(half) * u + np.sin(half) * w)
    s2 = x + bond * (np.cos(half) * u - np.sin(half) * w)
    return s1, s2


@dataclass
class Structure:
    """All-atom oligomer with fixed bond geometry; only dihedrals vary."""

    polymer: Polymer
    elements: list[str]
    coords: np.ndarray  # (n, 3)
    charges: np.ndarray  # (n,)
    bonds: list[tuple[int, int]]
    backbone: np.ndarray  # indices of the real backbone atoms in chain order
    n_dihedrals: int
    dihedrals: np.ndarray  # (N,) the RIS dihedrals used to build it
    subs_of: dict = field(default_factory=dict)  # backbone atom index -> substituent atom indices

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


def build_chain(polymer: Polymer, dihedrals_deg, cap: bool = True) -> Structure:
    """All-atom oligomer for ``N`` RIS dihedrals (N+3 backbone atoms), H-capped by default.

    The chain is built with one virtual backbone atom beyond each end (trans) so
    that terminal substituents are placed consistently; with ``cap=True`` an H atom
    is placed along each virtual bond, i.e. CH3-/CH2F-type end groups.
    """
    dih = np.asarray(dihedrals_deg, dtype=float)
    N = len(dih)
    B = polymer.bonds_per_repeat
    L = polymer.bond_length
    bb = build_backbone(polymer, dih, xp=np)  # (N+3, 3) real atoms, canonical frame
    # virtual (trans) backbone atoms beyond each end, in the *same* frame as ``bb`` so that
    # coordinates agree with build_backbone (helix analysis relies on this)
    v0 = nerf(bb[2], bb[1], bb[0], L, polymer.backbone[0].backbone_angle, 180.0, xp=np)
    v1 = nerf(bb[N], bb[N + 1], bb[N + 2], L, polymer.backbone[(N + 2) % B].backbone_angle, 180.0, xp=np)
    bb_ext = np.vstack([v0[None], bb, v1[None]])  # (N+5, 3); real atom k sits at bb_ext[k+1]
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
            bonds.append((backbone_idx[k - 1], idx))
        s1, s2 = substituent_positions(bb_ext[k], x, bb_ext[k + 2], spec.sub_bond, spec.sub_angle, xp=np)
        subs = []
        for s in (s1, s2):
            elements.append(spec.substituent)
            coords.append(s)
            charges.append(spec.sub_charge)
            bonds.append((idx, len(elements) - 1))
            subs.append(len(elements) - 1)
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
            bonds.append((idx, len(elements) - 1))
    return Structure(
        polymer=polymer,
        elements=elements,
        coords=np.array(coords),
        charges=np.array(charges),
        bonds=bonds,
        backbone=np.array(backbone_idx),
        n_dihedrals=N,
        dihedrals=dih,
        subs_of=subs_of,
    )


def build_chain_batch(polymer: Polymer, dihedrals_deg, cap: bool = True) -> tuple[Structure, np.ndarray]:
    """Batched version of :func:`build_chain`.

    ``dihedrals_deg`` has shape ``(M, N)``.  Returns ``(template, coords)`` where
    ``template`` is ``build_chain(polymer, dihedrals_deg[0], cap)`` (topology,
    elements, bonds -- identical for every row since only the dihedrals vary)
    and ``coords`` is an ``(M, n_atoms, 3)`` array whose row ``m`` equals
    ``build_chain(polymer, dihedrals_deg[m], cap).coords`` to numerical precision.

    Implemented with the batched :func:`build_backbone` plus :func:`substituent_positions`
    and :func:`nerf`, both of which already accept a leading batch dimension, so the whole
    oligomer (backbone, substituents, end caps) is placed in a handful of vectorised calls
    instead of one Python loop per conformer.
    """
    dih = np.asarray(dihedrals_deg, dtype=float)
    M, N = dih.shape
    B = polymer.bonds_per_repeat
    L = polymer.bond_length
    template = build_chain(polymer, dih[0], cap)
    bb = build_backbone(polymer, dih, xp=np)  # (M, N+3, 3)
    v0 = nerf(bb[:, 2], bb[:, 1], bb[:, 0], L, polymer.backbone[0].backbone_angle, 180.0, xp=np)
    v1 = nerf(bb[:, N], bb[:, N + 1], bb[:, N + 2], L, polymer.backbone[(N + 2) % B].backbone_angle, 180.0, xp=np)
    bb_ext = np.concatenate([v0[:, None, :], bb, v1[:, None, :]], axis=1)  # (M, N+5, 3)
    coords = np.zeros((M, template.n_atoms, 3), dtype=float)
    for k in range(N + 3):
        spec = polymer.backbone[k % B]
        idx = int(template.backbone[k])
        x = bb_ext[:, k + 1]
        coords[:, idx] = x
        s1, s2 = substituent_positions(bb_ext[:, k], x, bb_ext[:, k + 2], spec.sub_bond, spec.sub_angle, xp=np)
        s1_idx, s2_idx = template.subs_of[idx]
        coords[:, s1_idx] = s1
        coords[:, s2_idx] = s2
    if cap:
        cap0 = template.n_atoms - 2
        for i, (k, virt) in enumerate(((0, bb_ext[:, 0]), (N + 2, bb_ext[:, N + 4]))):
            idx = int(template.backbone[k])
            x = coords[:, idx]
            v = virt - x
            coords[:, cap0 + i] = x + 1.09 * v / np.linalg.norm(v, axis=-1, keepdims=True)
    return template, coords


def set_dihedrals(struct: Structure, dihedrals_deg) -> Structure:
    """Rebuild the same oligomer with new dihedrals (cheap: rigid geometry)."""
    cap = len(struct.elements) > 3 * (struct.n_dihedrals + 3)
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
