"""Continuous refinement of a packed crystal: chain conformation + cell together.

The discrete RIS stage uses ideal state angles and rigid geometry; real chains
deflect (beta-PVDF dihedrals ~ +/-172 deg, alpha-PVDF gauche ~ +/-45 deg).  Here
the conformation of one crystallographic repeat and the cell parameters are
relaxed together against the lattice energy, subject to the chain staying
periodic: the transform between consecutive repeats must be a pure translation.

Two parametrisations of "the chain stays periodic":

``parametrisation="linegroup"`` (default)
    The repeat is described by its line-group symmetry (:mod:`polyfind.linegroup`):
    the state sequence fixes a pattern for the torsions (a glide's
    ``phi[j+Q] = -phi[j]``, a screw's ``phi[j+Q] = phi[j]``), the backbone bond
    angles of the repeat join as variables, and the closure condition that makes
    the repeat transform a pure translation is solved for a couple of those
    parameters at every evaluation.  Periodicity then holds to ~1e-10 deg by
    construction, the penalty term is switched off, and the variable count drops
    from ``6 + P`` to ``6 + (pattern parameters + bond angles - closure rank)``.
``parametrisation="free"``
    The original method: all ``P`` torsions of the repeat are free (optionally
    plus the bond angles) and periodicity is enforced by a quadratic penalty on
    the residual rotation angle of the repeat transform (the linked-atom /
    helical-constraint idea of fibre-diffraction refinement).  Used automatically
    when the sequence has no recognised line-group pattern.

The objective is smooth, so the default method is L-BFGS-B with gradients, and
``gradient`` chooses where they come from.

``gradient="analytic"`` (default)
    :meth:`polyfind.pack.CrystalPacker.energy_and_grad` returns the energy of the
    current configuration together with closed-form derivatives with respect to the
    cell parameters, to the repeat's atom coordinates and to ``c``, from **one**
    kernel row.  The cell variables are then analytic outright.  For the
    conformational variables the remaining link is ``d(coords, c)/d(shape)``, which
    runs through the NeRF build, the closure solve and the Kabsch alignment; that
    Jacobian is taken by central differences of the *geometry only*, contracted
    against the exact ``dE/d(coords)`` and ``dE/dc``.  No kernel row is spent on it,
    so an iteration costs one row instead of ``1 + 2 n_vars`` (29 for the gamma
    chain), and the chains of the displacement batch are built in the single NeRF
    pass they were built in before.
``gradient="fd"``
    The original route: central differences of the *energy*, evaluated as one batched
    kernel call of ``1 + 2 n_vars`` configurations per iteration -- the cell
    displacements reusing the current chain and each conformational displacement a
    different chain, all carried as per-row coordinates by
    :meth:`polyfind.pack.CrystalPacker.energy`.

Either way the packer's topology-dependent tables (tiled force-field parameters,
exclusion scales, shift tables) are built once and kept across evaluations by
:meth:`polyfind.pack.CrystalPacker.update_chain`.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import minimize

from .linegroup import LineGroupError, line_group, repeat_chains, wrap180
from .pack import FD_STEPS, CrystalPacker, PackResult, PeriodicChain
from .polymers import Polymer


@dataclass
class RefineResult:
    start: PackResult
    result: PackResult
    torsions: np.ndarray
    rotation_error: float  # residual rotation per repeat (deg)
    n_evaluations: int
    angles: np.ndarray = field(default_factory=lambda: np.zeros(0))  # backbone angles of the repeat (deg)
    parametrisation: str = "free"
    n_variables: int = 0
    angle_energy: float = 0.0  # bond-angle strain per monomer (kcal/mol), 0 when the angles are frozen

    def summary(self) -> str:
        t = ", ".join(f"{x:6.1f}" for x in self.torsions)
        a = ", ".join(f"{x:5.1f}" for x in self.angles)
        return (
            f"{self.result.row()}\n    torsions: [{t}] deg; angles: [{a}] deg (strain {self.angle_energy:+.3f}); "
            f"commensurability residual {self.rotation_error:.3g} deg; "
            f"dE/mon = {self.result.energy_per_monomer - self.start.energy_per_monomer:+.3f} kcal/mol "
            f"({self.n_evaluations} evaluations of {self.n_variables} {self.parametrisation} variables)"
        )


def refine_crystal(
    polymer: Polymer,
    start: PackResult,
    torsions: np.ndarray | None = None,
    penalty: float = 5.0,
    maxfev: int = 3000,
    refine_torsions: bool = True,
    max_torsion_change: float = 40.0,
    cutoff: float = 8.0,
    eps_r: float = 1.0,
    gamma_free: bool = False,
    method: str = "lbfgs",
    maxiter: int = 200,
    h_torsion: float = 0.05,
    parametrisation: str = "linegroup",
    refine_angles: bool = True,
    max_angle_change: float = 8.0,
    angle_stiffness: float = 105.0,
    field=None,
    gradient: str = "analytic",
    probe: dict | None = None,
) -> RefineResult:
    """Relax cell (a, b, [gamma], phi1, phi2, dz) and the repeat's conformation.

    ``parametrisation`` is ``"linegroup"`` (default; falls back to ``"free"`` when the
    sequence has no recognised pattern) or ``"free"``.  ``refine_angles`` adds the
    repeat's backbone bond angles to the variables, bounded to the polymer's own
    values +/- ``max_angle_change``.  ``refine_torsions=False`` freezes the chain
    entirely and moves only the cell.

    ``angle_stiffness`` (kcal/mol/rad^2, UFF-like for sp3 carbon) is the harmonic
    bond-angle term ``0.5 k (theta - theta_0)^2`` summed over the repeat's backbone
    angles.  The packing kernel has no valence terms at all, so without it nothing
    resists opening the angles and they run straight to their bounds; it plays the
    same role for the angles as the Fourier torsion term does for the torsions.  It
    is reported separately as :attr:`RefineResult.angle_energy` and is *not* part of
    the lattice energy in :attr:`RefineResult.result`.

    ``method="lbfgs"`` (default) is L-BFGS-B with gradients; ``gradient="analytic"``
    (default) takes them from :meth:`polyfind.pack.CrystalPacker.energy_and_grad` in one
    kernel row per iteration, ``gradient="fd"`` from a batched central difference of
    ``1 + 2 n_vars`` rows (see the module docstring).  ``maxfev`` caps Nelder-Mead's
    function evaluations and is mapped onto L-BFGS-B's iteration limit.
    ``method="nelder-mead"`` is the original simplex search.

    ``field=(Ex, Ey, Ez)`` (V/A) relaxes the cell *and* the conformation in the presence
    of a uniform applied field, the lattice energy carrying the ``-mu_cell . E`` term of
    :meth:`polyfind.pack.CrystalPacker.field_energy`; the dipole follows the torsions and
    the setting angles, so the field acts on both.  ``field=None`` (the default) inherits
    the field ``start`` was packed under, so a refinement never silently drops it; pass
    ``field=(0, 0, 0)`` to refine a field-packed cell at zero field.

    ``probe``, if a dict is given, is filled before the optimiser runs with the objective
    (``value``), both gradient routes (``value_and_grad_analytic``,
    ``value_and_grad_fd``), the starting vector and the variable labels.  That is how the
    analytic gradient is checked component by component against a central difference of
    the very function the optimiser minimises, at an arbitrary point, rather than against
    a re-implementation of it; see ``tests/test_refine.py``.
    """
    if method not in ("lbfgs", "l-bfgs-b", "nelder-mead", "nelder_mead"):
        raise ValueError(f"unknown refine method {method!r}")
    if parametrisation not in ("linegroup", "free"):
        raise ValueError(f"unknown parametrisation {parametrisation!r}")
    if gradient not in ("analytic", "fd"):
        raise ValueError(f"unknown refine gradient {gradient!r} (expected 'analytic' or 'fd')")
    tors0 = np.asarray(start.dihedrals if torsions is None else torsions, dtype=float)
    B = polymer.bonds_per_repeat
    angles0 = np.array([polymer.backbone[k].backbone_angle for k in range(B)])
    p0 = start.params
    cell_idx = [0, 1, 3, 4, 5] + ([2] if gamma_free else [])
    n_cell = len(cell_idx)

    # --- the conformational block of the variable vector ------------------------
    lg = None
    if parametrisation == "linegroup" and refine_torsions:
        try:
            lg = line_group(polymer, start.chain, tors0, states=polymer.states, bond_angles=angles0)
        except LineGroupError:
            lg = None  # no pattern for this sequence: fall back to the penalty method
    mode = "linegroup" if lg is not None else "free"
    if not refine_torsions:
        n_shape = 0
        s0 = np.zeros(0)
        s_lo, s_hi, s_step = np.zeros(0), np.zeros(0), np.zeros(0)
    elif mode == "linegroup":
        n_shape = lg.n_params
        s0 = lg.x0
        lim = np.where(lg.is_angle(), max_angle_change, max_torsion_change)
        s_lo, s_hi = -lim, lim
        s_step = np.full(n_shape, h_torsion)
    else:
        n_t = len(tors0)
        s0 = np.concatenate([tors0, angles0]) if refine_angles else tors0.copy()
        n_shape = len(s0)
        s_lo = np.concatenate([np.full(n_t, -np.inf), angles0 - max_angle_change]) if refine_angles else np.full(n_t, -np.inf)
        s_hi = np.concatenate([np.full(n_t, np.inf), angles0 + max_angle_change]) if refine_angles else np.full(n_t, np.inf)
        s_step = np.full(n_shape, h_torsion)
    # The penalty stays in the objective for both paths, but with the line-group
    # parametrisation it does no work: the closure solve leaves a residual of ~1e-11 deg
    # (and ~2e-6 deg as pack measures it, where arccos near 1 loses half the mantissa),
    # so the term is below 1e-10 kcal/mol.  It only bites if a closure solve fails --
    # a helix pushed far off commensurability, say -- and then it correctly pushes the
    # optimiser back into the region where the chain really is periodic.
    pen_rot = penalty

    def conformations(S) -> tuple[np.ndarray, np.ndarray]:
        """(M, n_shape) -> torsions (M, P) and backbone angles (M, B)."""
        S = np.atleast_2d(np.asarray(S, dtype=float))
        M = S.shape[0]
        if n_shape == 0:
            return np.repeat(tors0[None], M, axis=0), np.repeat(angles0[None], M, axis=0)
        if mode == "linegroup":
            return lg.expand(S)
        if refine_angles:
            return S[:, : len(tors0)].copy(), S[:, len(tors0) :].copy()
        return S.copy(), np.repeat(angles0[None], M, axis=0)

    # The reference frame every chain of the refinement is expressed in: the setting
    # angles phi1, phi2 and the shift dz are measured against it, and aligning to it
    # makes the map conformation -> coordinates continuous (the raw principal-axis
    # frame flips by 180 deg as the torsions vary, which no gradient can follow).
    ref = repeat_chains(polymer, start.chain, tors0[None], angles0)[0]
    if field is None:
        field = getattr(start, "field", (0.0, 0.0, 0.0))
    field = None if not np.any(np.asarray(field, dtype=float)) else field
    packer = CrystalPacker(ref, n_chains=start.n_chains, cutoff=cutoff, eps_r=eps_r, field=field)
    count = {"n": 0}

    def chains_for(S) -> tuple[list[PeriodicChain], np.ndarray, np.ndarray]:
        tors, angs = conformations(S)
        chains = repeat_chains(polymer, start.chain, tors, angs, align_to=ref)
        return chains, tors, angs

    n_repeats = len(tors0) / B  # backbone angles of one crystallographic repeat, per angle type

    def bend_energy(angles) -> float:
        """Harmonic bond-angle strain of the whole cell (kcal/mol)."""
        d = np.deg2rad(np.asarray(angles, dtype=float) - angles0)
        return 0.5 * angle_stiffness * float(np.sum(d * d)) * n_repeats * start.n_chains

    def restraints(tors, angles, chain) -> float:
        """Bond-angle strain, commensurability penalty, and the caps on how far the chain may move."""
        pen = pen_rot * chain.rotation_error ** 2 + bend_energy(angles)
        dev = np.abs(wrap180(np.asarray(tors) - tors0))
        pen += 10.0 * float(np.sum(np.maximum(dev - max_torsion_change, 0) ** 2))
        dang = np.abs(np.asarray(angles) - angles0)
        pen += 10.0 * float(np.sum(np.maximum(dang - max_angle_change, 0) ** 2))
        return pen

    def unpack(x):
        p = p0.copy()
        p[cell_idx] = x[:n_cell]
        return x[n_cell:], p

    x0 = np.concatenate([p0[cell_idx], s0])
    if method.startswith("nelder"):
        def objective(x):
            count["n"] += 1
            s, p = unpack(x)
            if x[0] < 2.0 or x[1] < 2.0:
                return 1e6
            pen = 0.0
            if n_shape:
                lo_pen = np.maximum(s_lo - s, 0.0)
                hi_pen = np.maximum(s - s_hi, 0.0)
                pen += 10.0 * float(np.sum(np.where(np.isfinite(lo_pen), lo_pen, 0.0) ** 2 + np.where(np.isfinite(hi_pen), hi_pen, 0.0) ** 2))
            chains, tors, angs = chains_for(s[None])
            packer.update_chain(chains[0])
            return float(packer.energy(p[None])[0]) + restraints(tors[0], angs[0], chains[0]) + pen

        res = minimize(objective, x0, method="Nelder-Mead", options={"xatol": 1e-3, "fatol": 1e-4, "maxfev": maxfev, "adaptive": True})
        xbest = res.x
    else:
        def shape_batch(s):
            """The centre conformation and both displacements of every shape variable.

            One NeRF pass (and one closure solve) for the whole batch, as before; what
            changes between the two gradient routes is only what the chains are used for.
            """
            batch = [s]
            for j in range(n_shape):
                for sgn in (1.0, -1.0):
                    sj = s.copy()
                    sj[j] += sgn * s_step[j]
                    batch.append(sj)
            return chains_for(np.array(batch) if n_shape else s[None])

        def value_and_grad_analytic(x):
            count["n"] += 1
            s, p = unpack(x)
            chains, tors, angs = shape_batch(s)
            packer.update_chain(chains[0])
            # one kernel row: the energy and its closed-form derivatives w.r.t. the cell
            # parameters, the repeat's coordinates and the repeat length c
            E, g_cell, gX, gc = packer.energy_and_grad(p)
            grad = np.empty(n_cell + n_shape)
            grad[:n_cell] = g_cell[cell_idx]
            if n_shape:
                # The lattice energy sees a conformation only through (coords, c), and
                # dE/dcoords, dE/dc are exact here, so
                #     S(s) = gX . coords(s) + gc c(s) + E_torsion(s) + restraints(s)
                # has the same gradient in s as the objective does at this point.  Its
                # central difference therefore needs the geometry of each displacement
                # but not its energy: the remaining numerical link is the Jacobian of the
                # chain build (NeRF + closure solve + Kabsch alignment), not the kernel.
                S = np.array([
                    float((gX * ch.coords).sum()) + gc * ch.c
                    + packer.torsion_energy(ch.dihedrals) * packer.n_chains
                    + restraints(tors[j], angs[j], ch)
                    for j, ch in enumerate(chains)
                ])
                for j in range(n_shape):
                    grad[n_cell + j] = (S[1 + 2 * j] - S[2 + 2 * j]) / (2 * s_step[j])
            return float(E) + restraints(tors[0], angs[0], chains[0]), grad

        def value_and_grad_fd(x):
            count["n"] += 1
            s, p = unpack(x)
            # one NeRF pass (and one closure solve) for the centre chain and both
            # displacements of every conformational variable, then one kernel call
            chains, tors, angs = shape_batch(s)
            packer.update_chain(chains[0])
            rows, which = [p], [0]
            for idx in cell_idx:
                for sgn in (1.0, -1.0):
                    q = p.copy()
                    q[idx] += sgn * FD_STEPS[idx]
                    rows.append(q)
                    which.append(0)
            for j in range(1, len(chains)):
                rows.append(p)
                which.append(j)
            sel = [chains[i] for i in which]
            E = packer.energy(np.array(rows), **packer.chain_batch(sel))
            # row j of the batch carries conformation batch[j] and chain chains[j]; the
            # restraints are identical for the cell rows and cancel in their gradient
            F = E + np.array([restraints(tors[j], angs[j], chains[j]) for j in which])
            grad = np.empty(n_cell + n_shape)
            for i, idx in enumerate(cell_idx):
                grad[i] = (F[1 + 2 * i] - F[2 + 2 * i]) / (2 * FD_STEPS[idx])
            for j in range(n_shape):
                grad[n_cell + j] = (F[1 + 2 * n_cell + 2 * j] - F[2 + 2 * n_cell + 2 * j]) / (2 * s_step[j])
            return float(F[0]), grad

        def value_only(x):
            """The objective alone -- what both routes above return as their value."""
            s, p = unpack(x)
            chains, tors, angs = chains_for(s[None])
            packer.update_chain(chains[0])
            return float(packer.energy(p[None])[0]) + restraints(tors[0], angs[0], chains[0])

        value_and_grad = value_and_grad_analytic if gradient == "analytic" else value_and_grad_fd
        if probe is not None:
            shape_labels = [] if not n_shape else (
                lg.labels() if mode == "linegroup" else
                [f"t{j}" for j in range(len(tors0))] + ([f"angle{j}" for j in range(B)] if refine_angles else []))
            probe.update(
                x0=x0.copy(), value=value_only, value_and_grad_analytic=value_and_grad_analytic,
                value_and_grad_fd=value_and_grad_fd, n_cell=n_cell, n_shape=n_shape, mode=mode,
                cell_idx=list(cell_idx), s_step=np.asarray(s_step, dtype=float), packer=packer,
                labels=[CrystalPacker.PARAMS[i] for i in cell_idx] + shape_labels,
            )
        bounds = []
        for idx in cell_idx:
            bounds.append((2.0, None) if idx in (0, 1) else ((60.0, 120.0) if idx == 2 else (None, None)))
        bounds += [(lo if np.isfinite(lo) else None, hi if np.isfinite(hi) else None) for lo, hi in zip(s_lo, s_hi)]
        res = minimize(
            value_and_grad, x0, method="L-BFGS-B", jac=True, bounds=bounds,
            options={"ftol": 1e-9, "gtol": 1e-5, "maxiter": max(1, min(maxiter, maxfev)), "maxfun": maxfev},
        )
        xbest = res.x

    s, p = unpack(xbest)
    chains, tors, angs = chains_for(s[None])
    packer.update_chain(chains[0])
    out = packer.result(p)
    return RefineResult(
        start=start, result=out, torsions=chains[0].dihedrals.copy(), rotation_error=chains[0].rotation_error,
        n_evaluations=count["n"], angles=angs[0].copy(), parametrisation=mode, n_variables=n_cell + n_shape,
        angle_energy=bend_energy(angs[0]) / (start.n_chains * chains[0].n_monomers),
    )
