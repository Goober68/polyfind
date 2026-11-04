"""Helical symmetry of periodic torsion sequences.

A periodic sequence of P rotational states builds an infinite chain that is
mapped onto itself by a rigid *screw* motion (rotation ``theta`` about an axis
plus translation ``d`` along it).  We obtain that motion by superimposing one
period onto the next (Kabsch), then decompose it.  This gives the chain axis,
the rise per period, the crystallographic repeat ``c`` and the chain radius
without any fitting or optimisation, and it is what lets crystal packing be
searched in a handful of variables.

Also here: canonicalisation of periodic sequences under the symmetry
operations of the isolated chain (cyclic shift by a whole repeat unit, chain
reversal, mirror image), used to deduplicate candidates.  Which of those are
symmetries depends on the repeat: see :func:`sequence_images`.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .chain import build_backbone, pendant_positions, warn_if_chiral
from .polymers import Polymer, RISStates


# ------------------------------------------------------------------ rigid fits
def kabsch(P: np.ndarray, Q: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Rotation R and translation t with Q ~= P @ R.T + t (least squares)."""
    pc, qc = P.mean(axis=0), Q.mean(axis=0)
    H = (P - pc).T @ (Q - qc)
    U, _, Vt = np.linalg.svd(H)
    dsign = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1.0, 1.0, dsign])
    R = Vt.T @ D @ U.T
    t = qc - pc @ R.T
    return R, t


def screw_decompose(R: np.ndarray, t: np.ndarray, points: np.ndarray | None = None, tol_deg: float = 0.05):
    """Decompose x -> R x + t into a screw motion.

    Returns (axis n, angle theta in degrees, rise d, point p on the axis).
    ``n`` is oriented so that the rise ``d`` is positive; theta's sign then
    encodes handedness (right-hand rule about n).  For a pure translation the
    lateral position of the axis is undefined and ``p`` is set to the centroid
    of ``points``.
    """
    cos_t = np.clip((np.trace(R) - 1.0) / 2.0, -1.0, 1.0)
    theta = np.degrees(np.arccos(cos_t))
    if theta < tol_deg:  # pure translation
        d = np.linalg.norm(t)
        n = t / d if d > 0 else np.array([0.0, 0.0, 1.0])
        p = points.mean(axis=0) if points is not None else np.zeros(3)
        return n, 0.0, float(d), p
    if 180.0 - theta < tol_deg:
        w, v = np.linalg.eigh((R + R.T) / 2.0)
        n = v[:, np.argmax(w)]
    else:
        n = np.array([R[2, 1] - R[1, 2], R[0, 2] - R[2, 0], R[1, 0] - R[0, 1]])
        n /= np.linalg.norm(n)
    d = float(t @ n)
    if d < 0:
        n, d = -n, -d
    # signed angle about n: compare R applied to a vector perpendicular to n
    a = np.cross(n, [1.0, 0.0, 0.0])
    if np.linalg.norm(a) < 1e-6:
        a = np.cross(n, [0.0, 1.0, 0.0])
    a /= np.linalg.norm(a)
    b = R @ a
    theta = float(np.degrees(np.arctan2(np.cross(a, b) @ n, a @ b)))
    # axis point: (I - R) p = t - d n, with p . n = 0
    A = np.vstack([np.eye(3) - R, n[None, :]])
    rhs = np.concatenate([t - d * n, [0.0]])
    p, *_ = np.linalg.lstsq(A, rhs, rcond=None)
    return n, theta, d, p


# ------------------------------------------------------------------ helix params
@dataclass
class HelixParams:
    sequence: str
    period_bonds: int
    monomers_per_period: int
    rotation_per_period: float  # degrees, signed (handedness)
    rise_per_period: float  # A
    rise_per_bond: float  # A
    periods_per_repeat: int | None  # m such that m*theta = 0 mod 360 (None if incommensurate)
    turns_per_repeat: int | None
    c: float | None  # crystallographic repeat along the chain axis (A)
    radius_backbone: float  # max distance of backbone atoms from the axis (A)
    radius_all: float  # max distance of any atom from the axis (A)
    axis: np.ndarray
    axis_point: np.ndarray
    label: str

    def summary(self) -> str:
        c = f"{self.c:.2f}" if self.c is not None else "incommensurate"
        return (
            f"{self.sequence:<16s} P={self.period_bonds:2d} theta/period={self.rotation_per_period:7.1f} deg "
            f"rise/period={self.rise_per_period:5.2f} A  c={c:>6s} A  r_bb={self.radius_backbone:.2f} "
            f"r_all={self.radius_all:.2f}  {self.label}"
        )


def _commensurate(theta: float, max_m: int = 24, tol: float = 1.0):
    for m in range(1, max_m + 1):
        turns = m * theta / 360.0
        if abs(turns - round(turns)) * 360.0 / 1.0 <= tol * m ** 0.5:  # tolerance grows mildly with m
            return m, int(abs(round(turns)))
    return None, None


def kabsch_batch(P: np.ndarray, Q: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """:func:`kabsch` over a leading batch axis: ``P``, ``Q`` are ``(M, L, 3)``.

    Returns ``R`` ``(M, 3, 3)`` and ``t`` ``(M, 3)``.  Row ``m`` is ``kabsch(P[m], Q[m])``
    to the last bit: the centring, the ``3 x 3`` cross-covariance, the SVD, the sign
    fix and the translation are the same operations applied per matrix by NumPy's
    stacked linear algebra, not a differently ordered reduction.
    """
    pc, qc = P.mean(axis=1, keepdims=True), Q.mean(axis=1, keepdims=True)
    H = np.matmul(np.transpose(P - pc, (0, 2, 1)), Q - qc)
    U, _, Vt = np.linalg.svd(H)
    Ut, V = np.transpose(U, (0, 2, 1)), np.transpose(Vt, (0, 2, 1))
    dsign = np.sign(np.linalg.det(np.matmul(V, Ut)))
    VD = V.copy()
    VD[:, :, 2] *= dsign[:, None]
    R = np.matmul(VD, Ut)
    t = qc[:, 0] - np.matmul(pc, np.transpose(R, (0, 2, 1)))[:, 0]
    return R, t


def _period_atoms(polymer: Polymer, bb: np.ndarray, k0: int, k1: int) -> np.ndarray:
    """All-atom coordinates of backbone atoms ``k0 .. k1-1`` of the backbones ``bb``.

    ``bb`` is ``(M, N+3, 3)`` from a batched :func:`polyfind.chain.build_backbone`; the
    result is ``(M, n_atoms, 3)`` in :func:`polyfind.chain.build_chain`'s atom order --
    backbone atom, then pendant 1's atoms, then pendant 2's, for each ``k`` in turn -- so
    that reductions over it (the centroid of a pure-translation repeat) sum in the same
    order as over ``Structure.coords[atoms]`` and give the same bits.  The pendants are
    placed by :func:`polyfind.chain.pendant_positions` from the *real* backbone
    neighbours, which is what ``build_chain`` does for every atom that is not at a chain
    end; the callers here only ever ask for interior atoms (``k0 >= 1``,
    ``k1 <= N + 2``), where no virtual end atom is involved.
    """
    B = polymer.bonds_per_repeat
    if k0 < 1 or k1 > bb.shape[1] - 1:
        raise ValueError("_period_atoms needs interior backbone atoms (both neighbours real)")
    pieces = []
    for k in range(k0, k1):
        spec = polymer.backbone[k % B]
        x = bb[:, k]
        pieces.append(x)
        groups = pendant_positions(bb[:, k - 1], x, bb[:, k + 1], spec, xp=np)
        for g in groups:
            pieces.extend(g)
    return np.stack(pieces, axis=1)


def helix_parameters_batch(polymer: Polymer, seqs, states: RISStates, max_m: int = 24) -> list[HelixParams]:
    """:func:`helix_parameters` for several periodic sequences of the *same* period at once.

    ``seqs`` is ``(M, P)`` (state indices).  All ``M`` chains are built by one batched
    :func:`polyfind.chain.build_backbone` call and superimposed by one batched Kabsch;
    only the screw decomposition, the commensurability search and the radii (a few
    scalar operations each) are done per chain.  The descriptors are identical to
    ``M`` separate :func:`helix_parameters` calls: the batched NeRF, pendant placement
    and Kabsch apply the same floating-point operations per row.

    This is the routine :func:`polyfind.enumerate.enumerate_periodic` uses, once per
    period, for every new canonical sequence of that period; it is what turned the
    enumeration from a per-candidate rebuild of the chain (three times over) into a
    handful of batched builds.
    """
    seqs = np.asarray(seqs, dtype=int)
    if seqs.ndim == 1:
        seqs = seqs[None, :]
    M, P = seqs.shape
    B = polymer.bonds_per_repeat
    if P % B:
        raise ValueError("period must be a multiple of bonds_per_repeat")
    if M == 0:
        return []
    warn_if_chiral(polymer)
    L = max(P, 4)
    n_periods = int(np.ceil((2 * P + L + 3) / P)) + 1
    angles = np.array([float(a) for a in states.angles])
    dih = angles[np.tile(seqs, (1, n_periods))]  # (M, n_periods * P)
    bb = build_backbone(polymer, dih, xp=np)  # (M, N+3, 3), one NeRF recursion for the batch
    X = bb[:, P : P + L]
    Y = bb[:, 2 * P : 2 * P + L]
    R, t = kabsch_batch(X, Y)
    assert np.allclose(np.matmul(X, np.transpose(R, (0, 2, 1))) + t[:, None, :], Y, atol=1e-6), "periodic transform is not rigid?"
    # every atom of the period P .. 2P-1 (backbone + pendants), all rows at once
    atoms = _period_atoms(polymer, bb, P, 2 * P)  # (M, n_atoms, 3)
    names = states.names
    out = []
    for j in range(M):
        seq = seqs[j]
        block = bb[j, P : 2 * P]
        n, theta, d, p = screw_decompose(R[j], t[j], points=block)
        if abs(theta) < 1e-6:
            # pure translation: define the axis through the centroid of all atoms of one period
            p = atoms[j].mean(axis=0)
        m, turns = _commensurate(theta, max_m=max_m)

        def radius(pts):
            v = pts - p
            v -= np.outer(v @ n, n)
            return float(np.linalg.norm(v, axis=1).max())

        seq_name = "".join(names[s] for s in seq)
        all_trans = all(names[s] == "T" for s in seq)
        if all_trans:
            label = "planar zigzag"
        elif abs(theta) < 1e-6:
            label = f"glide/translation repeat of {P // B} monomer(s)"
        elif m is not None:
            label = f"{m}/{turns} helix (periods/turns), {'right' if theta > 0 else 'left'}-handed"
        else:
            label = f"helix, >{max_m} periods per repeat (not packed)"
        out.append(HelixParams(
            sequence=seq_name,
            period_bonds=P,
            monomers_per_period=P // B,
            rotation_per_period=theta,
            rise_per_period=d,
            rise_per_bond=d / P,
            periods_per_repeat=m,
            turns_per_repeat=turns,
            c=(m * d) if m is not None else None,
            radius_backbone=radius(block),
            radius_all=radius(atoms[j]),
            axis=n,
            axis_point=p,
            label=label,
        ))
    return out


def helix_parameters(polymer: Polymer, seq: np.ndarray, states: RISStates, max_m: int = 24) -> HelixParams:
    """Helix descriptors of the infinite chain generated by periodic sequence ``seq`` (state indices).

    The chain is built once (backbone by NeRF, pendants from the backbone) and the
    descriptors -- screw axis, rotation and rise per period, commensurability, radii --
    all come from that single build; see :func:`helix_parameters_batch` for many
    sequences of one period.
    """
    seq = np.asarray(seq, dtype=int)
    if seq.ndim != 1:
        raise ValueError("helix_parameters takes one sequence; use helix_parameters_batch for several")
    return helix_parameters_batch(polymer, seq[None, :], states, max_m=max_m)[0]


def rotation_to_z(n: np.ndarray) -> np.ndarray:
    """Rotation matrix taking unit vector n to +z."""
    n = n / np.linalg.norm(n)
    z = np.array([0.0, 0.0, 1.0])
    v = np.cross(n, z)
    s, c = np.linalg.norm(v), n @ z
    if s < 1e-12:
        return np.eye(3) if c > 0 else np.diag([1.0, -1.0, -1.0])
    vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + vx + vx @ vx * ((1 - c) / s ** 2)


# ---------------------------------------------------------------- symmetry
def sequence_images(
    seq,
    states: RISStates,
    bonds_per_repeat: int,
    reversal: bool = True,
    mirror: bool = True,
    chiral: bool = False,
):
    """All images of a periodic sequence under the symmetries of the isolated chain.

    A shift by a whole repeat unit is always one of them.  Which of chain reversal
    and reflection (the G+/G- state mirror) are is a property of the *repeat*:

    * an **achiral** repeat (PE, PVDF, PVDC) has both separately, so the images are
      shifts of ``{seq, reversed, mirrored, mirrored-and-reversed}``;
    * a **chiral** repeat (CFE, CDFE: a backbone atom with two different pendants is
      a stereocentre, and the builder makes the isotactic chain) has *neither*
      separately, only their composition.  Reflection maps the chain to its
      enantiomer; and reading an isotactic chain backwards swaps ``prev`` and
      ``next`` at every stereocentre, which flips each configuration relative to the
      chain direction, so plain reversal lands on the enantiomer too --
      ``E(phi_N..phi_1) == E(-phi_1..-phi_N)`` exactly, and both differ from
      ``E(phi_1..phi_N)`` (measured at up to 750 kcal/mol on a 24-bond CFE oligomer,
      against 0 for PVDF).  What survives is reflection *composed with* reversal,
      ``E(phi_1..phi_N) == E(-phi_N..-phi_1)``, the same operation that
      :func:`polyfind.forcefield.fit_ris` symmetrises a chiral fit over.  Pass
      ``chiral=True`` for that group: a sequence and its mirror are then distinct
      candidates, as they must be, since they are the right- and left-handed helices
      of a one-handed chain and have genuinely different energies.

    ``reversal`` and ``mirror`` switch off the corresponding generator; with
    ``chiral=True`` the surviving operation needs both, so switching off either
    leaves only the shifts.
    """
    seq = tuple(int(s) for s in seq)
    P = len(seq)
    B = bonds_per_repeat
    mirrored = lambda t: tuple(states.mirror[s] for s in t)  # noqa: E731
    base = [seq]
    if chiral:
        if reversal and mirror:
            base.append(mirrored(reversed(seq)))
    else:
        if reversal:
            base.append(tuple(reversed(seq)))
        if mirror:
            base += [mirrored(b) for b in list(base)]
    out = set()
    for b in base:
        for r in range(0, P, B):
            out.add(b[r:] + b[:r])
    return out


def canonical_sequence(
    seq,
    states: RISStates,
    bonds_per_repeat: int,
    reversal: bool = True,
    mirror: bool = True,
    chiral: bool = False,
) -> tuple:
    """Lexicographically smallest image (a hashable canonical form).

    ``chiral`` selects the smaller symmetry group of a chiral repeat; see
    :func:`sequence_images`.
    """
    return min(sequence_images(seq, states, bonds_per_repeat, reversal, mirror, chiral))


def is_primitive(seq, bonds_per_repeat: int) -> bool:
    """True if the sequence is not a repetition of a shorter period (multiple of B)."""
    seq = tuple(seq)
    P = len(seq)
    for q in range(bonds_per_repeat, P, bonds_per_repeat):
        if P % q == 0 and seq == seq[:q] * (P // q):
            return False
    return True
