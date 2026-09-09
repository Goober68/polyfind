"""Continuous refinement of a packed crystal: torsions + cell together.

The discrete RIS stage uses ideal state angles and rigid geometry; real chains
deflect (beta-PVDF dihedrals ~ +/-172 deg, alpha-PVDF gauche ~ +/-45 deg).  Here
the torsions of one crystallographic repeat and the cell parameters are relaxed
together against the lattice energy, subject to the chain staying periodic: the
transform between consecutive repeats must be a pure translation, enforced by a
quadratic penalty on its residual rotation angle (the linked-atom / helical-
constraint idea of fibre-diffraction refinement).

The objective is smooth, so the default method is L-BFGS-B on finite-difference
gradients obtained from *one* batched kernel call per iteration: the cell
displacements reuse the current chain, and each torsion displacement is a different
chain, all of them carried as per-row coordinates by
:meth:`polyfind.pack.CrystalPacker.energy`.  The chains of one gradient batch are
built in a single NeRF pass, and the packer's topology-dependent tables (tiled
force-field parameters, exclusion scales, shift tables) are built once and kept
across evaluations by :meth:`polyfind.pack.CrystalPacker.update_chain`.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import minimize

from .pack import (
    FD_STEPS,
    CrystalPacker,
    PackResult,
    PeriodicChain,
    periodic_chain_from_torsions,
    repeat_chains_from_torsions,
)
from .polymers import Polymer


@dataclass
class RefineResult:
    start: PackResult
    result: PackResult
    torsions: np.ndarray
    rotation_error: float  # residual rotation per repeat (deg)
    n_evaluations: int

    def summary(self) -> str:
        t = ", ".join(f"{x:6.1f}" for x in self.torsions)
        return (
            f"{self.result.row()}\n    torsions: [{t}] deg; commensurability residual {self.rotation_error:.3f} deg; "
            f"dE/mon = {self.result.energy_per_monomer - self.start.energy_per_monomer:+.3f} kcal/mol ({self.n_evaluations} evaluations)"
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
) -> RefineResult:
    """Relax cell (a, b, [gamma], phi1, phi2, dz) and the repeat's torsions from a packed start.

    ``method="lbfgs"`` (default) uses batched finite-difference gradients; ``maxfev``
    caps Nelder-Mead's function evaluations and is mapped onto L-BFGS-B's iteration
    limit.  ``method="nelder-mead"`` is the original simplex search.
    """
    if method not in ("lbfgs", "l-bfgs-b", "nelder-mead", "nelder_mead"):
        raise ValueError(f"unknown refine method {method!r}")
    tors0 = np.asarray(start.dihedrals if torsions is None else torsions, dtype=float)
    p0 = start.params
    cell_idx = [0, 1, 3, 4, 5] + ([2] if gamma_free else [])
    n_cell = len(cell_idx)
    n_t = len(tors0) if refine_torsions else 0
    x0 = np.concatenate([p0[cell_idx], tors0[:n_t]])
    count = {"n": 0}
    # the reference frame every chain of the refinement is expressed in: the setting
    # angles phi1, phi2 and the shift dz are measured against it, and aligning to it
    # makes the map torsions -> coordinates continuous (the raw principal-axis frame
    # flips by 180 deg as the torsions vary, which no gradient can follow)
    ref = periodic_chain_from_torsions(polymer, start.chain, tors0)
    packer = CrystalPacker(ref, n_chains=start.n_chains, cutoff=cutoff, eps_r=eps_r)

    def chains_for(tors_batch) -> list[PeriodicChain]:
        return repeat_chains_from_torsions(polymer, start.chain, np.atleast_2d(tors_batch), align_to=ref)

    def restraints(tors, chain) -> float:
        """Commensurability penalty plus the cap on how far the torsions may move."""
        pen = penalty * chain.rotation_error ** 2
        if n_t:
            dev = np.abs(np.asarray(tors)[:n_t] - tors0[:n_t])
            pen += 10.0 * float(np.sum(np.maximum(dev - max_torsion_change, 0) ** 2))
        return pen

    def unpack(x):
        tors = tors0.copy()
        if n_t:
            tors[:n_t] = x[n_cell:]
        p = p0.copy()
        p[cell_idx] = x[:n_cell]
        return tors, p

    if method.startswith("nelder"):
        def objective(x):
            count["n"] += 1
            tors, p = unpack(x)
            if x[0] < 2.0 or x[1] < 2.0:
                return 1e6
            chain = chains_for(tors)[0]
            packer.update_chain(chain)
            return float(packer.energy(p[None])[0]) + restraints(tors, chain)

        res = minimize(objective, x0, method="Nelder-Mead", options={"xatol": 1e-3, "fatol": 1e-4, "maxfev": maxfev, "adaptive": True})
        xbest = res.x
    else:
        def value_and_grad(x):
            count["n"] += 1
            tors, p = unpack(x)
            # one NeRF pass for the centre chain and both torsion displacements of each
            # free torsion, then one kernel call for the whole gradient
            batch = [tors]
            for j in range(n_t):
                for sgn in (1.0, -1.0):
                    tj = tors.copy()
                    tj[j] += sgn * h_torsion
                    batch.append(tj)
            chains = chains_for(np.array(batch))
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
            # row j of the batch carries torsion set batch[j] and chain chains[j]; the
            # restraints are identical for the cell rows and cancel in their gradient
            F = E + np.array([restraints(batch[j], chains[j]) for j in which])
            grad = np.empty(n_cell + n_t)
            for i, idx in enumerate(cell_idx):
                grad[i] = (F[1 + 2 * i] - F[2 + 2 * i]) / (2 * FD_STEPS[idx])
            for j in range(n_t):
                grad[n_cell + j] = (F[1 + 2 * n_cell + 2 * j] - F[2 + 2 * n_cell + 2 * j]) / (2 * h_torsion)
            return float(F[0]), grad

        bounds = []
        for idx in cell_idx:
            bounds.append((2.0, None) if idx in (0, 1) else ((60.0, 120.0) if idx == 2 else (None, None)))
        bounds += [(None, None)] * n_t
        res = minimize(
            value_and_grad, x0, method="L-BFGS-B", jac=True, bounds=bounds,
            options={"ftol": 1e-9, "gtol": 1e-5, "maxiter": max(1, min(maxiter, maxfev)), "maxfun": maxfev},
        )
        xbest = res.x

    tors, p = unpack(xbest)
    chain = chains_for(tors)[0]
    packer.update_chain(chain)
    out = packer.result(p)
    return RefineResult(start=start, result=out, torsions=chain.dihedrals.copy(), rotation_error=chain.rotation_error, n_evaluations=count["n"])
