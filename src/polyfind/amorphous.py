"""Amorphous and semi-crystalline ensembles from the transfer matrix.

No molecular dynamics or Monte Carlo equilibration: chains are drawn *exactly*
from the RIS Boltzmann distribution (:meth:`RISModel.sample`), batched on the
active backend, and their backbones are built in batch.  Clamping a segment to
a fixed state models a crystalline stem (chains entering/leaving a lamella).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import backend as bk
from .chain import build_backbone
from .polymers import Polymer
from .ris import RISModel


@dataclass
class EnsembleStats:
    temperature: float
    n_chains: int
    n_bonds: int
    state_fractions: dict
    diad_fractions: dict  # e.g. 'TT': 0.41
    trans_run_mean: float
    trans_run_p_ge: dict  # k -> probability that a bond lies in a trans run of length >= k
    tetrad_TTTT_fraction: float  # fraction of bond quadruples that are all trans (beta-like nuclei)
    mean_sq_end_to_end: float | None
    characteristic_ratio: float | None
    free_energy_per_monomer: float

    def summary(self) -> str:
        lines = [
            f"T = {self.temperature:.0f} K, {self.n_chains} chains x {self.n_bonds} bonds",
            "state fractions: " + ", ".join(f"{k}={v:.3f}" for k, v in self.state_fractions.items()),
            "diad fractions:  " + ", ".join(f"{k}={v:.3f}" for k, v in sorted(self.diad_fractions.items(), key=lambda kv: -kv[1])[:6]),
            f"mean trans-run length: {self.trans_run_mean:.2f} bonds; TTTT tetrad fraction: {self.tetrad_TTTT_fraction:.3f}",
            "P(bond in trans run >= k): " + ", ".join(f"k={k}:{p:.3f}" for k, p in self.trans_run_p_ge.items()),
            f"conformational free energy per monomer: {self.free_energy_per_monomer:.3f} kcal/mol",
        ]
        if self.characteristic_ratio is not None:
            lines.append(f"<R^2> = {self.mean_sq_end_to_end:.1f} A^2, characteristic ratio C_n = {self.characteristic_ratio:.2f}")
        return "\n".join(lines)


def sample_ensemble(polymer: Polymer, model: RISModel, n_bonds: int, n_chains: int, T: float, rng=None, clamp=None, build: bool = True):
    """Exact Boltzmann sample of ``n_chains`` chains; optionally build their backbones (batched)."""
    rng = rng or bk.default_rng(0)
    states = model.sample(n_bonds, T, n_chains, rng=rng, clamp=clamp)
    coords = None
    if build:
        xp = bk.get_backend()
        angles = xp.asarray(np.array(model.states.angles))
        dih = angles[states]
        coords = build_backbone(polymer, dih, xp=xp)
    return states, coords


def run_lengths(states: np.ndarray, s: int) -> np.ndarray:
    """Lengths of all maximal runs of state ``s`` across the rows of ``states``."""
    out = []
    for row in states:
        x = np.concatenate([[False], row == s, [False]]).astype(int)
        d = np.diff(x)
        starts, ends = np.where(d == 1)[0], np.where(d == -1)[0]
        out.append(ends - starts)
    return np.concatenate(out) if out else np.empty(0, dtype=int)


def ensemble_stats(polymer: Polymer, model: RISModel, states, coords, T: float, ks=(4, 6, 8, 12)) -> EnsembleStats:
    states = bk.to_numpy(states)
    n_chains, n_bonds = states.shape
    names = model.states.names
    S = model.S
    frac = {names[s]: float((states == s).mean()) for s in range(S)}
    diads = {}
    for a in range(S):
        for b in range(S):
            diads[names[a] + names[b]] = float(((states[:, :-1] == a) & (states[:, 1:] == b)).mean())
    t = model.states.index("T") if "T" in names else 0
    runs = run_lengths(states, t)
    trans_run_mean = float(runs.mean()) if runs.size else 0.0
    p_ge = {}
    for k in ks:
        # fraction of bonds lying in a trans run of length >= k
        p_ge[k] = float((runs[runs >= k]).sum() / (n_chains * n_bonds))
    tt = states == t
    tetrad = float((tt[:, :-3] & tt[:, 1:-2] & tt[:, 2:-1] & tt[:, 3:]).mean())
    msq, cn = None, None
    if coords is not None:
        c = bk.to_numpy(coords)
        r = c[:, -1] - c[:, 0]
        msq = float((r ** 2).sum(axis=1).mean())
        n_backbone_bonds = c.shape[1] - 1
        cn = msq / (n_backbone_bonds * polymer.bond_length ** 2)
    f = model.free_energy(n_bonds, T) / (n_bonds / polymer.bonds_per_repeat)
    return EnsembleStats(T, n_chains, n_bonds, frac, diads, trans_run_mean, p_ge, tetrad, msq, cn, f)


def lamella_interface(polymer: Polymer, model: RISModel, stem_bonds: int, tail_bonds: int, T: float, n_chains: int = 2000, rng=None):
    """Chains leaving a crystalline stem: first ``stem_bonds`` clamped trans, tail free.

    Returns the tail states, tail backbone coordinates (chain frame) and a dict of
    interface statistics: how far the tail end lies from the stem axis, the
    probability the tail re-enters the plane of the lamella (loop-like, z-turnback),
    and trans-run statistics in the tail relative to the bulk.
    """
    t = model.states.index("T")
    n_bonds = stem_bonds + tail_bonds
    clamp = {i: t for i in range(stem_bonds)}
    states, coords = sample_ensemble(polymer, model, n_bonds, n_chains, T, rng=rng, clamp=clamp, build=True)
    states = bk.to_numpy(states)
    c = bk.to_numpy(coords)
    # stem axis direction: from first to last stem backbone atom
    stem_end = stem_bonds + 2
    axis = c[:, stem_end] - c[:, 0]
    axis /= np.linalg.norm(axis, axis=1, keepdims=True)
    tail_vec = c[:, -1] - c[:, stem_end]
    along = (tail_vec * axis).sum(axis=1)
    perp = np.linalg.norm(tail_vec - along[:, None] * axis, axis=1)
    tail_states = states[:, stem_bonds:]
    bulk_states = model.sample(tail_bonds, T, n_chains, rng=bk.default_rng(1))
    stats = {
        "mean_advance_along_stem_axis": float(along.mean()),
        "p_turnback": float((along < 0).mean()),
        "mean_lateral_excursion": float(perp.mean()),
        "tail_trans_fraction": float((tail_states == t).mean()),
        "bulk_trans_fraction": float((bk.to_numpy(bulk_states) == t).mean()),
        "tail_first_bond_trans_fraction": float((tail_states[:, 0] == t).mean()),
    }
    return tail_states, c, stats
