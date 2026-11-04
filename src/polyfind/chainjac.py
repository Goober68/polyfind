"""Analytic Jacobian of the chain build: repeat torsions and bond angles -> coordinates.

The continuous refinement (:mod:`polyfind.refine`) and the closure solve of the
line-group parametrisation (:mod:`polyfind.linegroup`) both need the derivative of
the *geometry* of one crystallographic repeat with respect to the repeat's torsions
and backbone angles: the aligned coordinates of the repeat, its length ``c`` and the
rotation vector of the repeat transform.  Until now that Jacobian was taken by
central differences, which costs ``2 (P + R)`` chain builds per gradient (and per
Newton iteration of the closure solve).  Here it is closed-form, from the one build
already being made, in three layers.

**The oligomer.**  Rigid bond geometry means every internal coordinate is a
rotation generator.  Rotating RIS bond ``j`` (between backbone atoms ``j+1`` and
``j+2``) by ``d phi`` moves every atom downstream of it by ``u_j x (r - r_{j+2})``,
``u_j`` the unit bond vector; a pendant atom is rigid in the frame of its backbone
atom and its two neighbours, and for the atoms whose frames touch the bond the frame
turns rigidly with the downstream part because the bond's two atoms lie on the axis.
Opening the backbone angle at atom ``k`` by ``d theta`` rotates everything downstream
of ``k`` about the normal ``n_k`` of the plane ``(k-1, k, k+1)`` through ``r_k``, and
the pendants *of atom k* by half that angle about the same axis (their frame is built
on the bisector, which turns by ``d theta / 2``).  A repeat torsion or angle acts on
every one of its instances along the five-block oligomer, so its generator is the sum
over instances.

The coordinates being differentiated are the ones the refinement builds: the backbone
comes from :func:`polyfind.chain.build_backbone` and the pendants are placed on it by
:func:`polyfind.chain.extended_backbone` and :func:`polyfind.chain.place_pendants`,
which is :func:`polyfind.chain.build_oligomer_coords` step by step, so the oligomer
here is bit for bit the one :func:`polyfind.pack._batch_block_coords` and
:func:`polyfind.linegroup._block_coords` build.

**The repeat transform.**  The rigid fit ``Y ~= R X + t`` of two backbone windows
(:func:`polyfind.helix.kabsch`) maximises ``tr(R H)``, ``H = Xc^T Yc``, so at the
optimum ``A = R H`` is symmetric, and a perturbation ``dH`` moves ``R`` to
``(I + [delta]_x) R`` with ``([delta]_x A + A [delta]_x) = B^T - B``, ``B = R dH``;
for symmetric ``A`` the left side is ``[(tr A I - A) delta]_x``, a 3x3 solve.  From
``delta`` follow ``dt``, ``dc = t . dt / c``, the axis, and the rotation vector through
the inverse left Jacobian of SO(3).

**The frame.**  The repeat is expressed with its axis along z, oriented about z and
shifted along z to superimpose on a reference (:func:`polyfind.pack._orient_blocks`
and :func:`polyfind.pack._align_about_z_batch`).  The composite of those two about-z
steps is the *optimal* about-z rotation onto the reference, whatever the principal-axis
step chose, so the map is smooth and its derivative needs no eigendecomposition: any
rotation taking the axis to z serves, differing from ``rotation_to_z`` by an about-z
rotation the alignment absorbs, and the aligning angle is ``atan2(sw, cw)`` of the
correlation with the reference, differentiated directly.

Every array here carries a leading batch dimension ``M`` and a variable dimension
``nv = P + R`` (torsions first, then angles), and every derivative is per **degree**,
the unit the optimiser's variables are in.
"""
from __future__ import annotations

import numpy as np

from .chain import atom_starts, build_backbone, extended_backbone, place_pendants
from .pack import _rotations_to_z
from .polymers import Polymer

DEG = np.pi / 180.0


def _unit(v):
    return v / np.linalg.norm(v, axis=-1, keepdims=True)


_BB_CACHE: dict = {}


def _backbone_cached(polymer: Polymer, dih: np.ndarray, angles_batch):
    """``build_backbone`` with a one-entry memo on its inputs.

    The closure solve's last accepted Newton trial and the chain build that follows it
    ask for the very same backbone; the memo hands the second caller the first one's
    array (bit for bit: the function is deterministic in its inputs).
    """
    # keyed on the polymer's identity *and* the numbers the backbone build reads from
    # it, so two definitions sharing a name (a fitted variant, say) cannot collide
    key = (id(polymer), polymer.name, polymer.bond_length, tuple(b.backbone_angle for b in polymer.backbone),
           dih.tobytes(), dih.shape,
           None if angles_batch is None else (np.asarray(angles_batch, dtype=float).tobytes(), np.shape(angles_batch)))
    hit = _BB_CACHE.get(key)
    if hit is None:
        _BB_CACHE.clear()
        hit = _BB_CACHE[key] = build_backbone(polymer, dih, xp=np, bond_angles=angles_batch)
    return hit


def oligomer_and_jacobian(polymer: Polymer, tors_batch, angles_batch=None, pendants: bool = True):
    """Five-block oligomer coordinates and their Jacobian w.r.t. the repeat's variables.

    ``tors_batch`` is ``(M, P)``; ``angles_batch`` is ``None`` (the polymer's own angles,
    ``R = bonds_per_repeat``), ``(R,)`` or ``(M, R)``.  With ``pendants=False`` only the
    real backbone atoms are placed and differentiated (``coords is bb``, ``owner`` is
    ``arange``), which is all the repeat transform needs.

    Returns ``(coords, bb, jac, owner)``: ``coords`` ``(M, n_full, 3)`` in the atom order
    of ``build_chain(..., cap=False)`` (bit for bit what
    :func:`polyfind.chain.build_oligomer_coords` builds: the same backbone, the same
    :func:`polyfind.chain.place_pendants` on it), ``bb`` ``(M, N+3, 3)`` the real
    backbone atoms, ``jac`` ``(M, P + R, n_full, 3)`` = ``d coords / d (torsions, angles)``
    per degree, and ``owner`` ``(n_full,)`` the backbone atom each atom belongs to.

    The derivative is exact for every atom except the pendants of the two *end* backbone
    atoms (0 and N+2), whose frames hang off the virtual trans atoms and which no consumer
    of the middle block ever reads; their rows are left as if rigid.
    """
    tors_batch = np.atleast_2d(np.asarray(tors_batch, dtype=float))
    M, nb = tors_batch.shape
    dih = np.tile(tors_batch, (1, 5))
    N = dih.shape[1]
    B = polymer.bonds_per_repeat
    if angles_batch is None:
        R = B
    else:
        ab = np.asarray(angles_batch, dtype=float)
        R = ab.shape[-1]
    bb = _backbone_cached(polymer, dih, angles_batch)  # (M, N+3, 3)
    if pendants:
        # the pendants of the very builder the refinement uses, on this backbone
        bb_ext = extended_backbone(polymer, dih, angles_batch, xp=np, backbone=bb)
        coords = place_pendants(polymer, bb_ext, xp=np)
        starts = atom_starts(polymer, N + 3)
        n_full = int(starts[-1])
        owner = np.repeat(np.arange(N + 3), np.diff(starts))
        is_pendant = np.ones(n_full, dtype=bool)
        is_pendant[starts[:-1]] = False
    else:
        coords, n_full = bb, N + 3
        owner = np.arange(N + 3)
        is_pendant = np.zeros(n_full, dtype=bool)

    # --- torsions: bond j turns everything owned by atoms >= j + 2 about (r_{j+1}, r_{j+2})
    a = bb[:, 2 : N + 2]  # (M, N, 3) a point on each bond
    u = _unit(bb[:, 2 : N + 2] - bb[:, 1 : N + 1])
    mask_t = (owner[None, :] >= (np.arange(N) + 2)[:, None]).astype(float)  # (N, n_full)
    g = np.cross(u[:, :, None, :], coords[:, None, :, :] - a[:, :, None, :]) * mask_t[None, :, :, None]
    jac_t = g.reshape(M, 5, nb, n_full, 3).sum(axis=1)  # one repeat torsion = five instances

    # --- angles: the angle at atom k turns everything owned by atoms > k about n_k through
    # r_k, and the pendants of atom k itself by half that
    ks = np.arange(1, N + 2)
    n = _unit(np.cross(bb[:, ks - 1] - bb[:, ks], bb[:, ks + 1] - bb[:, ks]))  # (M, N+1, 3)
    w = np.where(owner[None, :] > ks[:, None], 1.0,
                 np.where((owner[None, :] == ks[:, None]) & is_pendant[None, :], 0.5, 0.0))
    ga = np.cross(n[:, :, None, :], coords[:, None, :, :] - bb[:, ks][:, :, None, :]) * w[None, :, :, None]
    jac_a = np.zeros((M, R, n_full, 3))
    types = ks % R
    for r in range(R):
        jac_a[:, r] = ga[:, types == r].sum(axis=1)
    jac = np.concatenate([jac_t, jac_a], axis=1) * DEG
    return coords, bb, jac, owner


def _vee(S):
    """Vector of a skew-symmetric matrix, batched over leading axes."""
    return np.stack([S[..., 2, 1], S[..., 0, 2], S[..., 1, 0]], axis=-1)


def _left_jacobian_inverse(w):
    """Inverse left Jacobian of SO(3) at rotation vector ``w`` (rad), batched ``(M, 3)`` -> ``(M, 3, 3)``."""
    th = np.linalg.norm(w, axis=-1)
    M = w.shape[0]
    K = np.zeros((M, 3, 3))
    K[:, 0, 1], K[:, 0, 2], K[:, 1, 0] = -w[:, 2], w[:, 1], w[:, 2]
    K[:, 1, 2], K[:, 2, 0], K[:, 2, 1] = -w[:, 0], -w[:, 1], w[:, 0]
    small = th < 1e-6
    ths = np.where(small, 1.0, th)
    coef = np.where(small, 1.0 / 12.0, 1.0 / ths ** 2 - (1.0 + np.cos(ths)) / (2.0 * ths * np.sin(ths)))
    return np.eye(3)[None] - 0.5 * K + coef[:, None, None] * (K @ K)


def repeat_transform_jacobian(bb, jac_bb, sx, sy):
    """The rigid fit between two backbone windows and its derivative.

    ``bb`` is ``(M, n_bb, 3)`` and ``jac_bb`` ``(M, nv, n_bb, 3)``; ``sx``, ``sy`` the two
    window slices.  Returns a dict with ``R`` ``(M, 3, 3)``, ``t`` ``(M, 3)``, ``c`` ``(M,)``,
    ``axis`` ``(M, 3)``, ``w`` ``(M, 3)`` the rotation vector (rad) and, per variable,
    ``delta`` ``(M, nv, 3)`` with ``dR = [delta]_x R``, ``dt`` ``(M, nv, 3)``, ``dc`` ``(M, nv)``,
    ``daxis`` ``(M, nv, 3)`` and ``dw`` ``(M, nv, 3)`` (rad per unit of the variable).

    ``R``, ``t`` and ``w`` are computed exactly as :func:`polyfind.linegroup.repeat_rotvec`
    and :func:`polyfind.helix.kabsch` compute them.
    """
    X, Y = bb[:, sx], bb[:, sy]
    pc, qc = X.mean(axis=1), Y.mean(axis=1)
    Xc, Yc = X - pc[:, None], Y - qc[:, None]
    H = np.einsum("mij,mik->mjk", Xc, Yc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(np.einsum("mji,mkj->mik", Vt, U)))
    V = np.transpose(Vt, (0, 2, 1)).copy()
    V[:, :, 2] *= d[:, None]
    R = np.einsum("mij,mkj->mik", V, U)  # V D U^T
    t = qc - np.einsum("mij,mj->mi", R, pc)
    c = np.linalg.norm(t, axis=1)
    axis = t / c[:, None]
    wv = 0.5 * np.stack([R[:, 2, 1] - R[:, 1, 2], R[:, 0, 2] - R[:, 2, 0], R[:, 1, 0] - R[:, 0, 1]], axis=1)
    s = np.linalg.norm(wv, axis=1)
    cth = 0.5 * (np.trace(R, axis1=1, axis2=2) - 1.0)
    theta = np.arctan2(s, cth)
    scale = np.where(s > 1e-12, theta / np.where(s > 1e-12, s, 1.0), 1.0)
    w = scale[:, None] * wv
    # --- derivative of the fit
    dX, dY = jac_bb[:, :, sx], jac_bb[:, :, sy]
    dpc, dqc = dX.mean(axis=2), dY.mean(axis=2)
    dXc, dYc = dX - dpc[:, :, None], dY - dqc[:, :, None]
    dH = np.einsum("mvij,mik->mvjk", dXc, Yc) + np.einsum("mij,mvik->mvjk", Xc, dYc)
    Bm = np.einsum("mij,mvjk->mvik", R, dH)
    rhs = _vee(np.swapaxes(Bm, -1, -2) - Bm)  # (M, nv, 3)
    A = R @ H
    Mm = np.trace(A, axis1=1, axis2=2)[:, None, None] * np.eye(3)[None] - A
    delta = np.swapaxes(np.linalg.solve(Mm, np.swapaxes(rhs, -1, -2)), -1, -2)  # (M, nv, 3)
    Rpc = np.einsum("mij,mj->mi", R, pc)
    dt = dqc - np.cross(delta, Rpc[:, None, :]) - np.einsum("mij,mvj->mvi", R, dpc)
    dc = np.einsum("mi,mvi->mv", t, dt) / c[:, None]
    daxis = (dt - axis[:, None, :] * dc[:, :, None]) / c[:, None, None]
    dw = np.einsum("mij,mvj->mvi", _left_jacobian_inverse(w), delta)
    return dict(R=R, t=t, c=c, axis=axis, w=w, delta=delta, dt=dt, dc=dc, daxis=daxis, dw=dw)


def aligned_block_jacobian(coords_block, jac_block, axis, daxis, ref):
    """Derivative of the aligned repeat coordinates.

    ``coords_block`` ``(M, n, 3)`` are the block's raw coordinates, ``jac_block``
    ``(M, nv, n, 3)`` their derivative, ``axis``/``daxis`` the repeat axis and its
    derivative, ``ref`` ``(n, 3)`` the reference frame the block is aligned to.  Returns
    ``(F, dF)``: the aligned coordinates -- equal to what
    ``_align_about_z_batch(_orient_blocks(...))`` produces up to rounding -- and their
    derivative ``(M, nv, n, 3)``.
    """
    M, n, _ = coords_block.shape
    point = coords_block.mean(axis=1)
    dpoint = jac_block.mean(axis=2)
    Rot = _rotations_to_z(axis)  # (M, 3, 3)
    W = np.einsum("mij,mnj->mni", Rot, coords_block - point[:, None])
    eps = np.cross(np.einsum("mij,mvj->mvi", Rot, daxis), np.array([0.0, 0.0, 1.0]))  # (M, nv, 3)
    dW = np.cross(eps[:, :, None, :], W[:, None]) + np.einsum("mij,mvnj->mvni", Rot, jac_block - dpoint[:, :, None])
    xr, yr = ref[:, 0], ref[:, 1]
    cw = (W[:, :, 0] * xr + W[:, :, 1] * yr).sum(axis=1)
    sw = (W[:, :, 0] * yr - W[:, :, 1] * xr).sum(axis=1)
    th = np.arctan2(sw, cw)
    dcw = (dW[..., 0] * xr + dW[..., 1] * yr).sum(axis=2)
    dsw = (dW[..., 0] * yr - dW[..., 1] * xr).sum(axis=2)
    dth = (cw[:, None] * dsw - sw[:, None] * dcw) / (cw ** 2 + sw ** 2)[:, None]
    ca, sa = np.cos(th), np.sin(th)
    F = np.empty_like(W)
    F[:, :, 0] = ca[:, None] * W[:, :, 0] - sa[:, None] * W[:, :, 1]
    F[:, :, 1] = sa[:, None] * W[:, :, 0] + ca[:, None] * W[:, :, 1]
    F[:, :, 2] = W[:, :, 2] - W[:, :, 2].mean(axis=1, keepdims=True) + ref[:, 2].mean()
    dF = np.empty_like(dW)
    rx = ca[:, None, None] * dW[..., 0] - sa[:, None, None] * dW[..., 1]
    ry = sa[:, None, None] * dW[..., 0] + ca[:, None, None] * dW[..., 1]
    dF[..., 0] = rx - dth[:, :, None] * F[:, None, :, 1]
    dF[..., 1] = ry + dth[:, :, None] * F[:, None, :, 0]
    dF[..., 2] = dW[..., 2] - dW[..., 2].mean(axis=2, keepdims=True)
    return F, dF
