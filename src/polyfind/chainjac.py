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


def oligomer(polymer: Polymer, tors_batch, angles_batch=None, pendants: bool = True):
    """The five-block oligomer alone: ``(coords, bb, owner, is_pendant)``.

    The build half of :func:`oligomer_and_jacobian` (same arguments, same arrays), for a caller
    that differentiates only some of the atoms (:func:`oligomer_jacobian`) or contracts the
    derivative with a gradient instead of forming it (:func:`oligomer_vjp`).
    """
    tors_batch = np.atleast_2d(np.asarray(tors_batch, dtype=float))
    dih = np.tile(tors_batch, (1, 5))
    N = dih.shape[1]
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
    return coords, bb, owner, is_pendant


def oligomer_jacobian(coords, bb, owner, is_pendant, nb: int, R: int, atoms=None):
    """``d coords[:, atoms] / d (torsions, angles)`` per degree, ``(M, nb + R, len(atoms), 3)``.

    ``atoms=None`` differentiates every atom.  A column of the result does not depend on which
    other columns were asked for (every operation is per atom), so a subset is bit for bit the
    same columns of the full Jacobian at a fraction of the work: the full array is
    ``(repeat variables) x (five blocks of atoms)`` when a consumer reads one block.
    """
    M = coords.shape[0]
    N = 5 * nb
    if atoms is not None:
        atoms = np.asarray(atoms, dtype=int)
        coords, owner, is_pendant = coords[:, atoms], owner[atoms], is_pendant[atoms]
    n_cols = coords.shape[1]
    # --- torsions: bond j turns everything owned by atoms >= j + 2 about (r_{j+1}, r_{j+2})
    a = bb[:, 2 : N + 2]  # (M, N, 3) a point on each bond
    u = _unit(bb[:, 2 : N + 2] - bb[:, 1 : N + 1])
    mask_t = (owner[None, :] >= (np.arange(N) + 2)[:, None]).astype(float)  # (N, n_cols)
    g = np.cross(u[:, :, None, :], coords[:, None, :, :] - a[:, :, None, :]) * mask_t[None, :, :, None]
    jac_t = g.reshape(M, 5, nb, n_cols, 3).sum(axis=1)  # one repeat torsion = five instances

    # --- angles: the angle at atom k turns everything owned by atoms > k about n_k through
    # r_k, and the pendants of atom k itself by half that
    ks = np.arange(1, N + 2)
    n = _unit(np.cross(bb[:, ks - 1] - bb[:, ks], bb[:, ks + 1] - bb[:, ks]))  # (M, N+1, 3)
    w = np.where(owner[None, :] > ks[:, None], 1.0,
                 np.where((owner[None, :] == ks[:, None]) & is_pendant[None, :], 0.5, 0.0))
    ga = np.cross(n[:, :, None, :], coords[:, None, :, :] - bb[:, ks][:, :, None, :]) * w[None, :, :, None]
    jac_a = np.zeros((M, R, n_cols, 3))
    types = ks % R
    for r in range(R):
        jac_a[:, r] = ga[:, types == r].sum(axis=1)
    return np.concatenate([jac_t, jac_a], axis=1) * DEG


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
    nb = tors_batch.shape[1]
    R = polymer.bonds_per_repeat if angles_batch is None else np.asarray(angles_batch, dtype=float).shape[-1]
    coords, bb, owner, is_pendant = oligomer(polymer, tors_batch, angles_batch, pendants)
    return coords, bb, oligomer_jacobian(coords, bb, owner, is_pendant, nb, R), owner


def oligomer_vjp(coords, bb, owner, is_pendant, nb: int, R: int, atoms, G):
    """``sum_n G[m, n] . d coords[m, atoms[n]] / d (torsions, angles)``, ``(M, nb + R)`` per degree.

    The contraction of :func:`oligomer_jacobian` with a gradient ``G`` ``(M, len(atoms), 3)``,
    without the Jacobian.  Every variable is a rotation generator, so its derivative of an
    energy is a torque: turning bond ``j`` moves the atoms downstream of it by
    ``u_j x (r - a_j)``, and ``sum_n G_n . (u_j x (r_n - a_j)) = u_j . (T - a_j x F)`` with
    ``F`` and ``T`` the total force and the total torque about the origin of those atoms.
    ``F`` and ``T`` for every bond at once are suffix sums along the chain, so the cost is that
    of one pass over the atoms where the Jacobian costs atoms times variables.
    """
    atoms = np.asarray(atoms, dtype=int)
    M = coords.shape[0]
    N = 5 * nb
    C = coords[:, atoms]
    own = owner[atoms]
    pend = is_pendant[atoms].astype(float)
    hot = np.zeros((N + 3, atoms.size))
    hot[own, np.arange(atoms.size)] = 1.0
    T = np.cross(C, G)
    F_own, T_own = hot @ G, hot @ T  # (M, N+3, 3): per owning backbone atom
    hot_p = hot * pend[None, :]
    Fp_own, Tp_own = hot_p @ G, hot_p @ T  # the pendants alone
    Fs = np.cumsum(F_own[:, ::-1], axis=1)[:, ::-1]  # owner >= o
    Ts = np.cumsum(T_own[:, ::-1], axis=1)[:, ::-1]
    # --- torsions: bond j turns everything owned by atoms >= j + 2
    a = bb[:, 2 : N + 2]
    u = _unit(bb[:, 2 : N + 2] - bb[:, 1 : N + 1])
    gt = (u * (Ts[:, 2 : N + 2] - np.cross(a, Fs[:, 2 : N + 2]))).sum(axis=-1)  # (M, N)
    gt = gt.reshape(M, 5, nb).sum(axis=1)
    # --- angles: atom k turns everything owned by atoms > k, and its own pendants by half
    ks = np.arange(1, N + 2)
    bk = bb[:, ks]
    nk = _unit(np.cross(bb[:, ks - 1] - bk, bb[:, ks + 1] - bk))
    down = Ts[:, ks + 1] - np.cross(bk, Fs[:, ks + 1])
    half = Tp_own[:, ks] - np.cross(bk, Fp_own[:, ks])
    ga_all = (nk * (down + 0.5 * half)).sum(axis=-1)  # (M, N+1)
    types = ks % R
    ga = np.zeros((M, R))
    for r in range(R):
        ga[:, r] = ga_all[:, types == r].sum(axis=1)
    return np.concatenate([gt, ga], axis=1) * DEG


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
    return window_transform_jacobian(bb[:, sx], bb[:, sy], jac_bb[:, :, sx], jac_bb[:, :, sy])


def _window_fit(X, Y):
    """The rigid fit of two windows: ``(pc, qc, Xc, Yc, H, R, t, c, axis, w)``."""
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
    return pc, qc, Xc, Yc, H, R, t, c, axis, w


def _fit_derivative(fit, dH, dpc, dqc) -> dict:
    """The derivative of the fit from that of ``H`` and of the two centroids, per variable."""
    pc, _, _, _, H, R, t, c, axis, w = fit
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


def window_transform_jacobian(X, Y, dX, dY):
    """:func:`repeat_transform_jacobian` on the two windows themselves: ``X``, ``Y`` ``(M, n_w, 3)``
    and their derivatives ``dX``, ``dY`` ``(M, nv, n_w, 3)``.  For a caller that has differentiated
    the window atoms alone (:func:`oligomer_jacobian` with ``atoms``)."""
    fit = _window_fit(X, Y)
    Xc, Yc = fit[2], fit[3]
    # --- derivative of the fit
    dpc, dqc = dX.mean(axis=2), dY.mean(axis=2)
    dXc, dYc = dX - dpc[:, :, None], dY - dqc[:, :, None]
    dH = np.einsum("mvij,mik->mvjk", dXc, Yc) + np.einsum("mij,mvik->mvjk", Xc, dYc)
    return _fit_derivative(fit, dH, dpc, dqc)


def _suffix(a):
    """``out[:, i] = a[:, i:].sum(axis=1)`` along the window axis, with a row of zeros at the end."""
    out = np.zeros((a.shape[0], a.shape[1] + 1) + a.shape[2:])
    out[:, :-1] = np.cumsum(a[:, ::-1], axis=1)[:, ::-1]
    return out


def _skew(e):
    """``[e]_x`` with ``[e]_x v = e x v``, batched over leading axes."""
    out = np.zeros(e.shape + (3,))
    out[..., 0, 1], out[..., 0, 2], out[..., 1, 0] = -e[..., 2], e[..., 1], e[..., 2]
    out[..., 1, 2], out[..., 2, 0], out[..., 2, 1] = -e[..., 0], -e[..., 1], e[..., 0]
    return out


def backbone_transform_jacobian(bb, sx, sy, nb: int, R: int):
    """:func:`repeat_transform_jacobian` from the backbone alone, without the windows' Jacobian.

    The fit's derivative needs three sums over the window atoms of their derivatives: that of
    ``H = Xc^T Yc`` and of the two centroids.  Every variable is a rotation generator that turns
    the atoms downstream of it -- axis ``e`` through ``a``, ``dr = e x (r - a)`` -- so each of
    those sums is ``[e]_x`` applied to a sum over the downstream window atoms of ``(X - a) Yc^T``
    and the like, and the downstream sums for every generator are running sums along the
    window.  A hundred generators cost a hundred 3x3 products, where the Jacobian of forty window
    atoms with respect to forty variables costs a five-block pass per atom.  The same numbers to
    rounding (not bit for bit: the sums are taken in another order).
    """
    M = bb.shape[0]
    N = 5 * nb
    X, Y = bb[:, sx], bb[:, sy]
    nw = X.shape[1]
    fit = _window_fit(X, Y)
    Xc, Yc = fit[2], fit[3]
    pc, qc = fit[0], fit[1]
    x0, y0 = sx.start, sy.start
    # running sums over the window atoms i >= a given one (row nw: none), of the centred windows:
    # the generators far up the chain turn both windows almost as one body, and sums about the
    # windows' own centroids keep the small difference that matters out of large numbers
    SH = _suffix(Xc[:, :, :, None] * Yc[:, :, None, :])  # sum Xc_i Yc_i^T
    SXc = _suffix(Xc)
    SYc = _suffix(Yc)
    count = (nw - np.arange(nw + 1)).astype(float)

    def generator(e, a, first):
        """Sums over the window atoms with backbone index >= ``first`` of their derivatives."""
        ix = np.clip(first - x0, 0, nw)
        iy = np.clip(first - y0, 0, nw)
        E = _skew(e)
        ax, ay = a - pc[:, None, :], a - qc[:, None, :]
        TX = SH[:, ix] - ax[:, :, :, None] * SYc[:, ix][:, :, None, :]   # sum (X_i - a) Yc_i^T
        TY = SH[:, iy] - SXc[:, iy][:, :, :, None] * ay[:, :, None, :]   # sum Xc_i (Y_i - a)^T
        dH = E @ TX - TY @ E
        dpc = np.cross(e, SXc[:, ix] - count[ix][None, :, None] * ax) / nw
        dqc = np.cross(e, SYc[:, iy] - count[iy][None, :, None] * ay) / nw
        return dH, dpc, dqc

    # torsions: bond j turns the backbone atoms >= j + 2 about (r_{j+1}, r_{j+2})
    j = np.arange(N)
    u = _unit(bb[:, 2 : N + 2] - bb[:, 1 : N + 1])
    tH, tp, tq = generator(u, bb[:, 2 : N + 2], j + 2)
    tH = tH.reshape(M, 5, nb, 3, 3).sum(axis=1)
    tp, tq = tp.reshape(M, 5, nb, 3).sum(axis=1), tq.reshape(M, 5, nb, 3).sum(axis=1)
    # angles: atom k turns the backbone atoms > k about its normal through itself
    ks = np.arange(1, N + 2)
    bk = bb[:, ks]
    n = _unit(np.cross(bb[:, ks - 1] - bk, bb[:, ks + 1] - bk))
    aH, ap, aq = generator(n, bk, ks + 1)
    kinds = (ks[None, :] % R == np.arange(R)[:, None]).astype(float)  # (R, N+1)
    aH = np.einsum("rk,mkij->mrij", kinds, aH)
    ap, aq = np.einsum("rk,mki->mri", kinds, ap), np.einsum("rk,mki->mri", kinds, aq)
    dH = np.concatenate([tH, aH], axis=1) * DEG
    dpc = np.concatenate([tp, ap], axis=1) * DEG
    dqc = np.concatenate([tq, aq], axis=1) * DEG
    return _fit_derivative(fit, dH, dpc, dqc)


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


def aligned_block_vjp(coords_block, axis, daxis, ref, G):
    """The contraction of :func:`aligned_block_jacobian` with a gradient, without the Jacobian.

    ``G`` ``(M, n, 3)`` is the derivative of an energy with respect to the aligned coordinates;
    ``ref`` is ``(n, 3)`` or one reference per row, ``(M, n, 3)``.  Returns ``(gv, Gp)``: ``gv``
    ``(M, nv)`` is the part of ``sum_n G_n . dF_n / dv`` that comes through the repeat axis
    (``daxis``), and ``Gp`` ``(M, n, 3)`` is the gradient carried back to the block's raw
    coordinates, to be contracted with their derivative (:func:`oligomer_vjp`):
    ``sum_n G_n . dF_n/dv = gv + sum_n Gp_n . d coords_block_n / dv``.

    The chain rule of :func:`aligned_block_jacobian` read backwards.  With ``W`` the block turned
    onto z and ``F`` the aligned block, ``dF = Rz(th) dW + dth (z x F)`` less the mean of
    ``dW_z``, and ``dth`` is itself a sum over the atoms of ``dW``; collecting the factor of each
    ``dW_n`` gives ``H_n`` below.  Then ``dW_n = eps x W_n + Rot (dC_n - mean dC)`` with ``eps``
    the turn of the axis, whose factor is the torque ``sum_n W_n x H_n``.
    """
    point = coords_block.mean(axis=1)
    Rot = _rotations_to_z(axis)  # (M, 3, 3)
    W = np.einsum("mij,mnj->mni", Rot, coords_block - point[:, None])
    xr, yr = ref[..., 0], ref[..., 1]
    cw = (W[:, :, 0] * xr + W[:, :, 1] * yr).sum(axis=1)
    sw = (W[:, :, 0] * yr - W[:, :, 1] * xr).sum(axis=1)
    th = np.arctan2(sw, cw)
    ca, sa = np.cos(th)[:, None], np.sin(th)[:, None]
    Fx = ca * W[:, :, 0] - sa * W[:, :, 1]
    Fy = sa * W[:, :, 0] + ca * W[:, :, 1]
    tau = ((G[:, :, 1] * Fx - G[:, :, 0] * Fy).sum(axis=1) / (cw ** 2 + sw ** 2))[:, None]
    H = np.empty_like(W)
    H[:, :, 0] = ca * G[:, :, 0] + sa * G[:, :, 1] + tau * (cw[:, None] * yr - sw[:, None] * xr)
    H[:, :, 1] = -sa * G[:, :, 0] + ca * G[:, :, 1] - tau * (cw[:, None] * xr + sw[:, None] * yr)
    H[:, :, 2] = G[:, :, 2] - G[:, :, 2].mean(axis=1, keepdims=True)
    L = np.cross(W, H).sum(axis=1)  # (M, 3)
    eps = np.cross(np.einsum("mij,mvj->mvi", Rot, daxis), np.array([0.0, 0.0, 1.0]))  # (M, nv, 3)
    gv = np.einsum("mvi,mi->mv", eps, L)
    Gp = np.einsum("mij,mni->mnj", Rot, H)
    Gp -= Gp.mean(axis=1, keepdims=True)
    return gv, Gp
