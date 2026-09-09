"""The multi-fidelity funnel, end to end.

    fit RIS  ->  enumerate periodic chains (DP)  ->  pack top-k (batched)  ->
    refine torsions + cell  ->  amorphous / interface statistics  ->  report

The expensive calculator is used only in the first stage (a few thousand
single points on a short oligomer) and, optionally, to re-score the final
candidates; every other stage runs on the cheap RIS / rigid-chain models.
"""
from __future__ import annotations

import time
from dataclasses import dataclass, field

import numpy as np

from .amorphous import EnsembleStats, ensemble_stats, lamella_interface, sample_ensemble
from .enumerate import Candidate, KNOWN_CHAINS, enumerate_periodic, table
from .forcefield import Calculator, FitReport, SimpleFF, fit_ris
from .helix import canonical_sequence
from .pack import PackResult, pack, periodic_chain
from .polymers import Polymer, get_polymer
from .refine import RefineResult, refine_crystal
from .ris import RISModel

# Reference unit cells (a, b, c in A, density g/cm^3), literature values, approximate.
EXPERIMENTAL_CELLS = {
    "pvdf": {
        "beta (TTTT)": (8.58, 4.91, 2.56, 1.97),
        "alpha/delta (TGTG')": (4.96, 9.64, 4.62, 1.92),
        "gamma/epsilon (T3GT3G')": (4.96, 9.67, 9.20, 1.94),
    },
    "pe": {"orthorhombic (all-trans)": (7.42, 4.95, 2.55, 1.00)},
}


@dataclass
class PipelineConfig:
    polymer: str = "pvdf"
    calculator: Calculator | None = None  # default SimpleFF()
    ris_model: RISModel | None = None  # skip fitting if given
    fit_step: float = 10.0
    fit_monomers: int = 6
    third_order: bool = True
    max_period: int = 8
    k_per_period: int = 60
    top_k_pack: int = 3
    pack_known: bool = True  # always pack the known polymorph chains too (for validation)
    max_atoms_per_chain: int = 48  # skip candidates whose crystallographic repeat is larger (cost ~ atoms^2)
    n_random: int = 3000
    n_refine: int = 4
    refine: bool = True
    refine_maxfev: int = 2000
    temperature: float = 450.0  # for the amorphous ensemble (PVDF melt ~ 450 K)
    n_chains: int = 2000
    n_bonds: int = 200
    seed: int = 0


@dataclass
class PipelineResult:
    polymer: Polymer
    model: RISModel
    fit: FitReport | None
    candidates: list[Candidate]
    packed: list[tuple[Candidate, PackResult]]
    refined: list[RefineResult]
    amorphous: EnsembleStats | None
    interface: dict | None
    timings: dict = field(default_factory=dict)

    def report(self) -> str:
        p = self.polymer
        out = [f"polyfind report for {p.name} {p.formula}", "=" * 72]
        out.append(self.model.describe())
        if self.fit:
            out.append(f"(fitted from {self.fit.n_evaluations} calculator evaluations in {self.timings.get('fit', 0):.1f} s)")
        out.append("")
        out.append(f"Periodic chain conformations (RIS energy per monomer, {self.timings.get('enumerate', 0):.2f} s):")
        out.append(table(self.candidates, top=12))
        out.append("")
        out.append(f"Crystal packing of {len(self.packed)} chain conformations ({self.timings.get('pack', 0):.1f} s):")
        expt = EXPERIMENTAL_CELLS.get(p.name, {})
        for cand, res in self.packed:
            out.append(f"  {res.row()}" + (f"   [{cand.known_as}]" if cand.known_as else ""))
            if cand.known_as in expt:
                a, b, c, rho = expt[cand.known_as]
                out.append(f"  {'':14s} experiment:      a={a:5.2f} b={b:5.2f} c={c:5.2f}  rho={rho:5.3f}")
        if self.refined:
            out.append("")
            out.append(f"Continuous refinement (torsions + cell, {self.timings.get('refine', 0):.1f} s):")
            for r in self.refined:
                out.append("  " + r.summary().replace("\n", "\n  "))
        if self.packed:
            out.append("")
            best = min(self.packed, key=lambda cr: cr[1].energy_per_monomer)
            out.append("Lattice energy ranking (kcal/mol per monomer, relative to the best):")
            for cand, res in sorted(self.packed, key=lambda cr: cr[1].energy_per_monomer):
                out.append(f"  {res.chain:<14s} {res.energy_per_monomer - best[1].energy_per_monomer:+7.3f}" + (f"   {cand.known_as}" if cand.known_as else ""))
        if self.amorphous:
            out.append("")
            out.append(f"Amorphous ensemble ({self.timings.get('amorphous', 0):.2f} s, exact Boltzmann sampling):")
            out.append("  " + self.amorphous.summary().replace("\n", "\n  "))
        if self.interface:
            out.append("  chains leaving a crystalline stem: " + ", ".join(f"{k}={v:.3f}" for k, v in self.interface.items()))
        return "\n".join(out)


def run_pipeline(cfg: PipelineConfig, verbose: bool = True) -> PipelineResult:
    log = print if verbose else (lambda *a, **k: None)
    polymer = get_polymer(cfg.polymer)
    calc = cfg.calculator or SimpleFF()
    rng = np.random.default_rng(cfg.seed)
    timings = {}

    # 1. RIS parameters
    t0 = time.time()
    fit = None
    if cfg.ris_model is not None:
        model = cfg.ris_model
    else:
        fit = fit_ris(polymer, calc, step=cfg.fit_step, n_monomers=cfg.fit_monomers, third_order=cfg.third_order)
        model = fit.model
    timings["fit"] = time.time() - t0
    log(f"[1/5] RIS model ready ({timings['fit']:.1f} s)")

    # 2. periodic chain conformations
    t0 = time.time()
    cands = enumerate_periodic(polymer, model, max_period=cfg.max_period, k_per_period=cfg.k_per_period)
    timings["enumerate"] = time.time() - t0
    log(f"[2/5] {len(cands)} distinct periodic chain conformations up to period {cfg.max_period} ({timings['enumerate']:.2f} s)")

    # 3. packing
    t0 = time.time()
    def n_atoms(c):
        return 3 * c.period * (c.helix.periods_per_repeat or 0)

    packable = [c for c in cands if c.helix.c is not None and n_atoms(c) <= cfg.max_atoms_per_chain]
    skipped = [c for c in cands if c.helix.c is not None and n_atoms(c) > cfg.max_atoms_per_chain]
    for c in skipped[:5]:
        log(f"      not packing {c.name} ({n_atoms(c)} atoms per chain repeat > {cfg.max_atoms_per_chain})")
    to_pack = list(packable[: cfg.top_k_pack])
    if cfg.pack_known:
        for c in packable:
            if c.known_as and c not in to_pack:
                to_pack.append(c)
    packed = []
    for cand in to_pack:
        try:
            chain = periodic_chain(polymer, cand.seq, model.states, helix=cand.helix)
        except (ValueError, RuntimeError) as e:
            log(f"      skip {cand.name}: {e}")
            continue
        res = pack(chain, n_chains=2, n_random=cfg.n_random, n_refine=cfg.n_refine, rng=rng)
        packed.append((cand, res[0]))
        log(f"      packed {cand.name:<14s} E/mon={res[0].energy_per_monomer:8.3f}  a={res[0].a:.2f} b={res[0].b:.2f} c={res[0].c:.2f}")
    timings["pack"] = time.time() - t0
    log(f"[3/5] packed {len(packed)} chain conformations ({timings['pack']:.1f} s)")

    # 4. refinement
    refined = []
    t0 = time.time()
    if cfg.refine:
        for cand, res in packed:
            refined.append(refine_crystal(polymer, res, maxfev=cfg.refine_maxfev))
        # replace packed energies by refined ones for the ranking
        packed = [(cand, r.result) for (cand, _), r in zip(packed, refined)]
    timings["refine"] = time.time() - t0
    log(f"[4/5] refinement done ({timings['refine']:.1f} s)")

    # 5. amorphous ensemble
    t0 = time.time()
    states, coords = sample_ensemble(polymer, model, cfg.n_bonds, cfg.n_chains, cfg.temperature, rng=rng)
    amorphous = ensemble_stats(polymer, model, states, coords, cfg.temperature)
    _, _, interface = lamella_interface(polymer, model, stem_bonds=12, tail_bonds=60, T=cfg.temperature, n_chains=min(cfg.n_chains, 2000), rng=rng)
    timings["amorphous"] = time.time() - t0
    log(f"[5/5] amorphous statistics ({timings['amorphous']:.2f} s)")
    return PipelineResult(polymer, model, fit, cands, packed, refined, amorphous, interface, timings)
