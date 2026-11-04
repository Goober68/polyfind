"""Line-group parametrisation of one crystallographic repeat.

The refinement of :mod:`polyfind.refine` may only move a chain inside the set of
conformations whose repeat-to-repeat transform is a *pure translation*: anything
else is not a chain of an infinite crystal.  Writing ``R(y)`` for the rotation
part of that transform as a function of the repeat's torsions and backbone
angles ``y``, the admissible set is ``{y : R(y) = I}``.  The penalty method of
the original refinement pushes towards that set from outside; this module
parametrises it directly, which removes the stiff term from the objective and
cuts the variable count.

Two ingredients.

**The symmetry pattern.**  A periodic state sequence carries its own line-group
symmetry, and the symmetry is readable from the sequence itself.  Let ``s`` be
the state sequence of the repeat (length ``P``) and ``m`` the mirror map of the
RIS states (``G+ <-> G-``, ``T -> T``).  Then

* if ``s[j+Q] == m(s[j])`` for all ``j``, for some ``Q`` with ``P / Q`` even,
  the chain has a *glide*-type symmetry: mirror composed with a shift by ``Q``
  bonds.  The torsions inherit ``phi[j+Q] = -phi[j]``, so only ``Q`` of them are
  free.  ``T`` maps to itself under the mirror and ``-180 == 180``, so a run of
  trans states is free to deflect, ``(180-d, ..., 180+d, ...)``;
* otherwise, if ``s[j+Q] == s[j]`` for some ``Q`` dividing ``P``, the repeat is
  ``P/Q`` turns of a *screw* and the torsions inherit ``phi[j+Q] = phi[j]``;
* otherwise no pattern is recognised and every torsion stays free.

A glide contains a reflection, and a chiral object admits no improper isometry:
a chain whose repeat carries a stereocentre (CFE, CDFE -- see
:attr:`polyfind.polymers.Polymer.is_chiral`) therefore does *not* have one, however
glide-like its state sequence reads, because the pendants break it even where the
torsions do not.  Imposing ``phi[j+Q] = -phi[j]`` there would restrict the
refinement to a subspace the true minimum need not lie in, so glide detection is
conditional on the repeat being achiral (``torsion_pattern(..., chiral=True)``,
which :func:`line_group` passes for a chiral polymer).  A screw is a proper
isometry -- a rotation and a translation -- and stays valid for both.  A chiral
sequence that would have matched a glide falls through to the screw test and then
to ``"free"``, i.e. to the penalty method, which is the documented behaviour for
any unrecognised sequence.

For PVDF this gives: ``TT`` -> glide with ``Q = 1``, one free deflection;
``TG+TG-`` (alpha) -> glide with ``Q = 2``, free ``(t, g)`` expanding to
``(t, g, -t, -g)``; ``TTTG+TTTG-`` (gamma) -> glide with ``Q = 4``, free
``(t1, t2, t3, g)``; ``TG+`` repeated -> screw with ``Q = 2``, free ``(t, g)``.

**The closure condition.**  Imposing the symmetry alone does *not* make the
repeat a translation: with a glide, ``R`` is the square of an improper isometry
``S(omega, n)``, i.e. a rotation by ``2 omega``, so ``R = I`` collapses from
three scalar conditions to one, but that one still has to be solved.  (This is
the geometry that makes an all-trans chain with unequal backbone angles curve:
its repeat is a rotation by ``theta_1 - theta_2`` per monomer.)  So the free
parameters here are the symmetry pattern's parameters *plus* the backbone
angles, minus the rank of the closure Jacobian: a few of them are solved for by
Newton's method at every evaluation, and the resulting chain is periodic to
1e-10 deg without any penalty term.

Everything is batched: one gradient batch of chains is built, closed and
assembled in a handful of vectorised calls.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .chain import build_backbone, build_oligomer_coords
from .chainjac import aligned_block_jacobian, oligomer_and_jacobian, repeat_transform_jacobian
from .helix import canonical_sequence
from .pack import PeriodicChain, _block_atoms, _finish_repeat_chains, _template
from .polymers import Polymer, RISStates


class LineGroupError(ValueError):
    """The sequence has no usable line-group parametrisation."""


def wrap180(x):
    """Angle(s) folded into (-180, 180]."""
    w = (np.asarray(x, dtype=float) + 180.0) % 360.0 - 180.0
    return np.where(w == -180.0, 180.0, w)


# ------------------------------------------------------------------ the pattern
@dataclass(frozen=True)
class TorsionPattern:
    """Linear map from a few free torsion parameters to the repeat's ``P`` torsions.

    ``torsions = signs * params[group]``; ``kind`` is ``"glide"``, ``"screw"`` or
    ``"free"`` and ``period`` is the ``Q`` of the docstring above.
    """

    kind: str
    period: int
    group: tuple[int, ...]
    signs: tuple[float, ...]
    sequence: str
    canonical: str

    @property
    def n_torsions(self) -> int:
        return len(self.group)

    @property
    def n_params(self) -> int:
        return max(self.group) + 1

    @property
    def labels(self) -> list[str]:
        return [f"t{i}" for i in range(self.n_params)]

    def torsions(self, params) -> np.ndarray:
        """(..., n_params) -> (..., P)."""
        p = np.asarray(params, dtype=float)
        return np.asarray(self.signs) * p[..., np.asarray(self.group)]

    def params_from_torsions(self, torsions) -> tuple[np.ndarray, float]:
        """Least-squares (circular) inverse and the residual it leaves behind."""
        t = np.asarray(torsions, dtype=float)
        group, signs = np.asarray(self.group), np.asarray(self.signs)
        p = np.empty(self.n_params)
        for i in range(self.n_params):
            v = signs[group == i] * t[group == i]
            p[i] = v[0] + wrap180(v - v[0]).mean()
        err = float(np.abs(wrap180(self.torsions(p) - t)).max()) if len(t) else 0.0
        return p, err


def state_sequence(name: str, states: RISStates, length: int | None = None) -> list[int]:
    """Parse a sequence name (``"TG+TG-"``) into state indices, tiled to ``length``."""
    order = sorted(range(states.n), key=lambda i: -len(states.names[i]))
    seq, i = [], 0
    while i < len(name):
        for s in order:
            nm = states.names[s]
            if name.startswith(nm, i):
                seq.append(s)
                i += len(nm)
                break
        else:
            raise LineGroupError(f"cannot parse sequence name {name!r}")
    if length is not None:
        if not seq or length % len(seq):
            raise LineGroupError(f"sequence {name!r} does not tile a repeat of {length} bonds")
        seq = seq * (length // len(seq))
    return seq


def torsion_pattern(seq, states: RISStates, bonds_per_repeat: int = 1, chiral: bool = False) -> TorsionPattern:
    """Symmetry pattern of a periodic state sequence (see the module docstring).

    ``chiral`` marks a repeat with a stereocentre, whose chain has no improper
    symmetry: glide detection is then skipped and such a sequence falls back to the
    screw test and, failing that, to ``"free"``.
    """
    seq = [int(s) for s in seq]
    P = len(seq)
    name = "".join(states.names[s] for s in seq)
    canon = "".join(states.names[s] for s in canonical_sequence(seq, states, bonds_per_repeat, chiral=chiral)) if P else ""
    if not chiral:  # a chiral chain has no improper symmetry, so no glide
        for Q in range(1, P):  # glide: mirror composed with a shift by Q bonds
            if P % (2 * Q) == 0 and all(seq[(j + Q) % P] == states.mirror[seq[j]] for j in range(P)):
                group = tuple(j % Q for j in range(P))
                signs = tuple(1.0 if (j // Q) % 2 == 0 else -1.0 for j in range(P))
                return TorsionPattern("glide", Q, group, signs, name, canon)
    for Q in range(1, P):  # screw: a plain shift by Q bonds
        if P % Q == 0 and all(seq[(j + Q) % P] == seq[j] for j in range(P)):
            return TorsionPattern("screw", Q, tuple(j % Q for j in range(P)), (1.0,) * P, name, canon)
    return TorsionPattern("free", P, tuple(range(P)), (1.0,) * P, name, canon)


# ------------------------------------------------------------------ geometry
def _windows(nb: int):
    """The two backbone windows whose rigid fit is the repeat transform.

    Matches :func:`polyfind.pack.repeat_chains_from_torsions` exactly, so that the
    residual solved here is the ``rotation_error`` reported there.
    """
    if nb < 3:
        return slice(2 * nb, 2 * nb + 4), slice(3 * nb, 3 * nb + 4)
    return slice(2 * nb, 3 * nb), slice(3 * nb, 4 * nb)


def repeat_rotvec(polymer: Polymer, torsions, bond_angles=None) -> np.ndarray:
    """Rotation vector (deg) of the repeat transform, batched over rows.

    ``torsions`` is ``(M, P)`` and ``bond_angles`` ``None``, ``(R,)`` or ``(M, R)``.
    The norm of the returned vector is the ``rotation_error`` of the chain that
    :mod:`polyfind.pack` would build from the same torsions.
    """
    tors = np.atleast_2d(np.asarray(torsions, dtype=float))
    nb = tors.shape[1]
    bb = build_backbone(polymer, np.tile(tors, (1, 5)), xp=np, bond_angles=bond_angles)
    sx, sy = _windows(nb)
    X, Y = bb[:, sx], bb[:, sy]
    Xc = X - X.mean(axis=1, keepdims=True)
    Yc = Y - Y.mean(axis=1, keepdims=True)
    H = np.einsum("mij,mik->mjk", Xc, Yc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(np.einsum("mji,mkj->mik", Vt, U)))  # det(V U^T)
    V = np.transpose(Vt, (0, 2, 1)).copy()
    V[:, :, 2] *= d[:, None]
    R = np.einsum("mij,mkj->mik", V, U)  # V D U^T
    w = 0.5 * np.stack([R[:, 2, 1] - R[:, 1, 2], R[:, 0, 2] - R[:, 2, 0], R[:, 1, 0] - R[:, 0, 1]], axis=1)
    s = np.linalg.norm(w, axis=1)
    c = 0.5 * (np.trace(R, axis1=1, axis2=2) - 1.0)
    theta = np.arctan2(s, c)
    scale = np.where(s > 1e-12, theta / np.where(s > 1e-12, s, 1.0), 1.0)
    return np.degrees(scale[:, None] * w)


def _block_coords(polymer: Polymer, tors_batch: np.ndarray, angles_batch=None) -> np.ndarray:
    """(M, nb) torsions [+ (M, R) angles] -> (M, n_full, 3) five-block oligomer coordinates.

    :func:`polyfind.pack._batch_block_coords` with the backbone angles opened up; the
    atom order is that of ``build_chain(..., cap=False)``.
    """
    return build_oligomer_coords(polymer, np.tile(tors_batch, (1, 5)), angles_batch, xp=np)


def repeat_chains(
    polymer: Polymer,
    name: str,
    torsions_batch,
    bond_angles=None,
    scale14: float = 0.5,
    align_to: "PeriodicChain | np.ndarray | None" = None,
) -> list[PeriodicChain]:
    """:func:`polyfind.pack.repeat_chains_from_torsions` with per-row backbone angles.

    ``bond_angles`` is ``None`` (the polymer's own angles), ``(R,)`` for the whole
    batch or ``(M, R)`` per row.  With ``bond_angles=None`` the result is identical
    to ``pack.repeat_chains_from_torsions``.
    """
    tors_batch = np.atleast_2d(np.asarray(torsions_batch, dtype=float))
    coords_full = _block_coords(polymer, tors_batch, bond_angles)
    # the rigid fit, orientation and alignment of every row are batched matrix
    # operations in pack (``_finish_repeat_chains``, windows as in ``_windows``)
    return _finish_repeat_chains(polymer, name, tors_batch, coords_full, scale14, align_to)


def repeat_chains_and_jacobian(
    polymer: Polymer,
    name: str,
    torsions_batch,
    bond_angles,
    align_to,
    scale14: float = 0.5,
):
    """:func:`repeat_chains` plus the analytic Jacobian of what the packer reads from it.

    ``bond_angles`` is ``(R,)`` or ``(M, R)`` and ``align_to`` is required: the derivative
    is that of the *aligned* coordinates, which is what the refinement differentiates.
    Returns ``(chains, dcoords, dc, w, dw)`` with ``chains`` exactly what
    :func:`repeat_chains` builds (the same oligomer from
    :func:`polyfind.chainjac.oligomer_and_jacobian`, the same batched tail
    :func:`polyfind.pack._finish_repeat_chains`), ``dcoords`` ``(M, P + R, n, 3)`` =
    ``d chain.coords / d (torsions, angles)``, ``dc`` ``(M, P + R)``, ``w`` ``(M, 3)`` the
    rotation vector of the repeat transform in degrees (its norm is
    ``chain.rotation_error``) and ``dw`` ``(M, P + R, 3)`` its derivative -- all per
    degree.  See :mod:`polyfind.chainjac`.
    """
    tors_batch = np.atleast_2d(np.asarray(torsions_batch, dtype=float))
    M, nb = tors_batch.shape
    coords_full, bb, jac, _ = oligomer_and_jacobian(polymer, tors_batch, bond_angles, pendants=True)
    chains = _finish_repeat_chains(polymer, name, tors_batch, coords_full, scale14, align_to)
    tpl = _template(polymer, nb)
    _, atoms0 = _block_atoms(tpl, nb, 0)
    sx, sy = _windows(nb)
    tr = repeat_transform_jacobian(bb, jac[:, :, tpl.backbone], sx, sy)
    ref = align_to.coords if isinstance(align_to, PeriodicChain) else np.asarray(align_to, dtype=float)
    _, dF = aligned_block_jacobian(coords_full[:, atoms0], jac[:, :, atoms0], tr["axis"], tr["daxis"], ref)
    return chains, dF, tr["dc"], np.degrees(tr["w"]), np.degrees(tr["dw"])


# ------------------------------------------------------------------ parametrisation
@dataclass
class LineGroup:
    """Free-parameter vector of an exactly periodic repeat.

    The full parameter vector is ``y = (pattern torsion parameters, backbone angles)``.
    ``dep`` of its entries are solved for so that the repeat transform is a pure
    translation; the remaining ``free`` entries are what the optimiser sees, and
    ``expand`` turns them into torsions and angles.
    """

    polymer: Polymer
    name: str
    pattern: TorsionPattern
    y0: np.ndarray  # reference parameters (the start conformation)
    free: np.ndarray  # indices of y that stay free
    dep: np.ndarray  # indices of y solved for by the closure condition
    jac0: np.ndarray  # closure Jacobian w.r.t. y[dep] at the reference
    tors0: np.ndarray  # reference torsions (for the wrap convention)
    tol: float = 1e-9
    max_newton: int = 12
    last_residual: float = 0.0
    n_builds: int = 0  # backbone builds spent in closure solves (diagnostics)

    # --- construction -------------------------------------------------------
    @property
    def n_params(self) -> int:
        return len(self.free)

    @property
    def n_torsions(self) -> int:
        return self.pattern.n_params

    @property
    def kind(self) -> str:
        return self.pattern.kind

    def labels(self) -> list[str]:
        names = self.pattern.labels + [f"angle{i}" for i in range(len(self.y0) - self.pattern.n_params)]
        return [names[i] for i in self.free]

    def is_angle(self) -> np.ndarray:
        """Boolean mask: which free parameters are backbone angles."""
        return self.free >= self.pattern.n_params

    @property
    def x0(self) -> np.ndarray:
        return np.zeros(len(self.free))

    # --- evaluation ---------------------------------------------------------
    def _split(self, Y):
        nt = self.pattern.n_params
        tors = self.pattern.torsions(Y[:, :nt])
        tors = self.tors0 + wrap180(tors - self.tors0)
        return tors, Y[:, nt:]

    def _residual(self, Y):
        tors, angles = self._split(Y)
        return repeat_rotvec(self.polymer, tors, angles)

    def _jacobian(self, Y, cols=None, h: float = 1e-3):
        """d(residual)/d(y[cols]) for every row, by central differences in one batched build."""
        cols = self.dep if cols is None else np.asarray(cols, dtype=int)
        M, k = Y.shape[0], len(cols)
        big = np.repeat(Y, 2 * k, axis=0)
        for i, j in enumerate(cols):
            big[2 * i :: 2 * k, j] += h
            big[2 * i + 1 :: 2 * k, j] -= h
        r = self._residual(big).reshape(M, k, 2, 3)
        return np.transpose((r[:, :, 0] - r[:, :, 1]) / (2 * h), (0, 2, 1))

    @property
    def tmat(self) -> np.ndarray:
        """``d torsions / d (pattern parameters)``, ``(P, n_t)``: ``torsions = tmat @ params``."""
        T = np.zeros((len(self.pattern.group), self.pattern.n_params))
        T[np.arange(len(self.pattern.group)), np.asarray(self.pattern.group)] = np.asarray(self.pattern.signs)
        return T

    def _residual_and_jac(self, Y):
        """The residual (M, 3) and its closed-form Jacobian ``d r / d y`` (M, 3, n_y), deg per deg.

        One backbone build per call, the same one :func:`repeat_rotvec` makes; the
        residual is bit for bit :meth:`_residual`'s.  See :mod:`polyfind.chainjac`.
        """
        tors, angles = self._split(Y)
        nb = tors.shape[1]
        _, bb, jac, _ = oligomer_and_jacobian(self.polymer, tors, angles, pendants=False)
        sx, sy = _windows(nb)
        tr = repeat_transform_jacobian(bb, jac, sx, sy)
        dw = np.degrees(tr["dw"])  # (M, P + R, 3), per degree of torsion / angle
        J = np.concatenate([np.einsum("mpi,pt->mit", dw[:, :nb], self.tmat), np.transpose(dw[:, nb:], (0, 2, 1))], axis=2)
        return np.degrees(tr["w"]), J

    def _close(self, Y, max_step: float = 20.0):
        """Newton on y[dep] until the repeat transform is a pure translation.

        The Jacobian is closed-form and comes with every residual (one backbone build
        each, no finite differences), so an iteration costs one build when the full
        step is accepted.  Damped: a step that makes a row's residual worse is halved
        for that row, which keeps far-from-reference parameters (a helix pushed away
        from commensurability, say) from running away.  A row that is already inside
        tolerance is frozen (its step is never taken), so where a row ends up does not
        depend on which other rows share its batch.  A row that still cannot be
        closed keeps its residual, and the caller's penalty term -- inactive at 1e-6
        deg, enormous at 10 deg -- pushes the optimiser back out of that region.

        Returns ``(Y, iterations, J)`` with ``J`` the closure Jacobian ``(M, 3, n_y)`` at
        the returned ``Y``.
        """
        Y = np.array(Y, dtype=float)
        r, J = self._residual_and_jac(Y)
        self.n_builds += 1
        for it in range(self.max_newton):
            norm = np.linalg.norm(r, axis=1)
            if norm.max() < self.tol:
                self.last_residual = float(norm.max())
                return Y, it, J
            step = np.clip(-np.einsum("mij,mj->mi", np.linalg.pinv(J[:, :, self.dep]), r), -max_step, max_step)
            Ybest, rbest, Jbest = Y.copy(), r.copy(), J.copy()
            for _ in range(4):
                Yt = Y.copy()
                Yt[:, self.dep] += step
                rt, Jt = self._residual_and_jac(Yt)
                self.n_builds += 1
                # a row that is already closed stays put, so that its result does not
                # depend on which other rows share its batch
                take = (np.linalg.norm(rt, axis=1) < np.linalg.norm(rbest, axis=1)) & (norm >= self.tol)
                Ybest[take], rbest[take], Jbest[take] = Yt[take], rt[take], Jt[take]
                bad = (np.linalg.norm(rbest, axis=1) >= norm) & (norm > self.tol)
                if not bad.any():
                    break
                step = np.where(bad[:, None], 0.5 * step, step)
            Y, r, J = Ybest, rbest, Jbest
        self.last_residual = float(np.linalg.norm(r, axis=1).max())
        return Y, self.max_newton, J

    def expand(self, x) -> tuple[np.ndarray, np.ndarray]:
        """Free parameters (n,) or (M, n) -> torsions (M, P) and backbone angles (M, B)."""
        X = np.atleast_2d(np.asarray(x, dtype=float))
        Y = np.repeat(self.y0[None], X.shape[0], axis=0)
        Y[:, self.free] += X
        Y, _, _ = self._close(Y)
        return self._split(Y)

    def expand_jac(self, x) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """:meth:`expand` plus ``d (torsions, angles) / d x``, ``(M, P + R, n_params)``.

        The closure is an implicit function: ``r(y) = 0`` with ``y[dep]`` solved for, so
        ``d y[dep] / d x = -pinv(J[:, dep]) J[:, free]`` at the solution, with ``J`` the
        closure Jacobian the solve ends on.  Torsions follow the pattern's linear map.
        """
        X = np.atleast_2d(np.asarray(x, dtype=float))
        M = X.shape[0]
        Y = np.repeat(self.y0[None], M, axis=0)
        Y[:, self.free] += X
        Y, _, J = self._close(Y)
        tors, angles = self._split(Y)
        n_y, nt = len(self.y0), self.pattern.n_params
        dYdx = np.zeros((M, n_y, len(self.free)))
        dYdx[:, self.free, np.arange(len(self.free))] = 1.0
        if len(self.dep):
            dYdx[:, self.dep] = -np.einsum("mij,mjk->mik", np.linalg.pinv(J[:, :, self.dep]), J[:, :, self.free])
        dVdy = np.zeros((tors.shape[1] + angles.shape[1], n_y))
        dVdy[: tors.shape[1], :nt] = self.tmat
        dVdy[tors.shape[1] :, nt:] = np.eye(n_y - nt)
        return tors, angles, np.einsum("vy,myk->mvk", dVdy, dYdx)

    def chains(self, x, name: str | None = None, align_to=None, scale14: float = 0.5) -> list[PeriodicChain]:
        tors, angles = self.expand(x)
        return repeat_chains(self.polymer, name or self.name, tors, angles, scale14=scale14, align_to=align_to)


def line_group(
    polymer: Polymer,
    name: str,
    torsions,
    states: RISStates | None = None,
    bond_angles=None,
    tol: float = 1e-9,
) -> LineGroup:
    """Build the :class:`LineGroup` of the repeat ``torsions`` of sequence ``name``.

    Raises :class:`LineGroupError` when the sequence has no recognised pattern, when
    the given torsions are not of that pattern, or when nothing is left free after
    the closure condition -- the caller then falls back to the penalty method.  A
    chiral repeat has no glide (see the module docstring), so a glide-patterned
    sequence of e.g. CFE takes that fallback rather than a symmetry its chain does
    not have.
    """
    states = states or polymer.states
    tors0 = np.asarray(torsions, dtype=float)
    P = len(tors0)
    seq = state_sequence(name, states, length=P)
    pattern = torsion_pattern(seq, states, polymer.bonds_per_repeat, chiral=bool(polymer.is_chiral))
    if pattern.kind == "free":
        raise LineGroupError(f"no line-group pattern for sequence {name!r}")
    p0, err = pattern.params_from_torsions(tors0)
    if err > 1e-6:
        raise LineGroupError(f"torsions do not follow the {pattern.kind} pattern of {name!r} ({err:.3g} deg)")
    B = polymer.bonds_per_repeat
    a0 = np.array([polymer.backbone[k].backbone_angle for k in range(B)]) if bond_angles is None else np.asarray(bond_angles, dtype=float)
    y0 = np.concatenate([p0, a0])
    lg = LineGroup(polymer, name, pattern, y0, np.arange(len(y0)), np.zeros(0, dtype=int),
                   np.zeros((3, 0)), tors0, tol=tol)
    r0 = lg._residual(y0[None])[0]
    if np.abs(r0).max() > 1e-6:
        raise LineGroupError(f"start conformation of {name!r} is not periodic ({np.abs(r0).max():.3g} deg)")
    # Closure Jacobian at the reference.  Its rank is the number of parameters that have
    # to be solved for; they are picked greedily, torsion parameters before backbone
    # angles (so the angles stay free wherever the symmetry allows it) and within each
    # group by how much a column adds to the span already selected.
    J = lg._residual_and_jac(y0[None])[1][0]  # (3, n_y), closed-form
    sv = np.linalg.svd(J, compute_uv=False)
    k = int((sv > 1e-4 * max(sv.max(), 1e-12)).sum())
    dep, basis = [], np.zeros((3, 0))
    for j in list(range(pattern.n_params)) + list(range(pattern.n_params, len(y0))):
        if len(dep) >= k:
            break
        v = J[:, j] - basis @ (basis.T @ J[:, j])
        if np.linalg.norm(v) > 1e-3 * max(np.linalg.norm(J[:, j]), 1e-12):
            dep.append(j)
            basis = np.column_stack([basis, v / np.linalg.norm(v)])
    if len(dep) < k or len(y0) - len(dep) < 1:
        raise LineGroupError(f"closure leaves no usable free parameters for {name!r}")
    lg.dep = np.array(sorted(dep))
    lg.free = np.array([j for j in range(len(y0)) if j not in dep])
    lg.jac0 = J[:, lg.dep]
    Y, iters, _ = lg._close(y0[None])
    if iters >= lg.max_newton or np.abs(Y[0] - y0).max() > 1e-6:
        raise LineGroupError(f"closure solve did not reproduce the start conformation of {name!r}")
    return lg
