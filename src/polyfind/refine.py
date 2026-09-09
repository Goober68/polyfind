"""Continuous refinement of a packed crystal: torsions + cell together.

The discrete RIS stage uses ideal state angles and rigid geometry; real chains
deflect (beta-PVDF dihedrals ~ +/-172 deg, alpha-PVDF gauche ~ +/-45 deg).  Here
the torsions of one crystallographic repeat and the cell parameters are relaxed
together against the lattice energy, subject to the chain staying periodic: the
transform between consecutive repeats must be a pure translation, enforced by a
quadratic penalty on its residual rotation angle (the linked-atom / helical-
constraint idea of fibre-diffraction refinement).

Cost: a few thousand lattice-energy evaluations per candidate, each on a
handful of atoms; it is the last cheap stage before an optional high-fidelity
relaxation with an external calculator.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import minimize

from .pack import CrystalPacker, PackResult, PeriodicChain, periodic_chain_from_torsions
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
) -> RefineResult:
    """Relax cell (a, b, [gamma], phi1, phi2, dz) and the repeat's torsions from a packed start."""
    tors0 = np.asarray(start.dihedrals if torsions is None else torsions, dtype=float)
    p0 = start.params
    cell_idx = [0, 1, 3, 4, 5] + ([2] if gamma_free else [])
    n_cell = len(cell_idx)
    n_t = len(tors0) if refine_torsions else 0
    x0 = np.concatenate([p0[cell_idx], tors0[:n_t]])
    count = {"n": 0}

    def build(x):
        tors = tors0.copy()
        if n_t:
            tors[:n_t] = x[n_cell:]
        chain = periodic_chain_from_torsions(polymer, start.chain, tors)
        packer = CrystalPacker(chain, n_chains=start.n_chains, cutoff=cutoff, eps_r=eps_r)
        p = p0.copy()
        p[cell_idx] = x[:n_cell]
        return chain, packer, p

    def objective(x):
        count["n"] += 1
        pen = 0.0
        if n_t:
            dev = np.abs(x[n_cell:] - tors0[:n_t])
            pen += 10.0 * np.sum(np.maximum(dev - max_torsion_change, 0) ** 2)
        if x[0] < 2.0 or x[1] < 2.0:
            return 1e6
        chain, packer, p = build(x)
        e = float(packer.energy(p[None])[0])
        return e + penalty * chain.rotation_error ** 2 + pen

    res = minimize(objective, x0, method="Nelder-Mead", options={"xatol": 1e-3, "fatol": 1e-4, "maxfev": maxfev, "adaptive": True})
    chain, packer, p = build(res.x)
    out = packer.result(p)
    return RefineResult(start=start, result=out, torsions=chain.dihedrals.copy(), rotation_error=chain.rotation_error, n_evaluations=count["n"])
