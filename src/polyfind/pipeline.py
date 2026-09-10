"""The multi-fidelity funnel, end to end.

    fit RIS  ->  enumerate periodic chains (DP)  ->  pack top-k (batched)  ->
    refine torsions + cell  ->  amorphous / interface statistics  ->  report

The expensive calculator is used only in the first stage (a few thousand
single points on a short oligomer) and, optionally, to re-score the final
candidates; every other stage runs on the cheap RIS / rigid-chain models.
"""
from __future__ import annotations

import os
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, field

import numpy as np

from .amorphous import EnsembleStats, ensemble_stats, lamella_interface, sample_ensemble
from .enumerate import Candidate, KNOWN_CHAINS, enumerate_periodic, table
from .forcefield import Calculator, FitReport, SimpleFF, fit_ris
from .helix import canonical_sequence
from .pack import PackResult, pack, periodic_chain
from .polymers import Polymer, get_polymer
from .lattice_table import physical_cores as _physical_cores
from .refine import RefineResult, refine_crystal
from .ris import RISModel

# Reference unit cells (a, b, c in A, density g/cm^3), literature values, approximate.
# Every entry is sourced and given a verdict in docs/REFERENCES.md section 1; consult
# that before treating any of these as ground truth.  Two caveats recorded there:
#   * the gamma cell is monoclinic with beta ~ 93 deg, which this table does not carry
#     (the search fits the cell angle freely, so nothing computed depends on it);
#   * the gamma density of 1.94 does not follow from the gamma cell beside it, which
#     gives 1.93.  Left alone deliberately -- changing it shifts a reported table.
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
    workers: int | None = None  # candidates packed/refined in parallel processes; None = a bandwidth-aware default, see run_pipeline; 1 = serial


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


def _pack_and_refine_one(
    polymer_name: str,
    model_dict: dict,
    seq: tuple,
    seed: int,
    n_random: int,
    n_refine: int,
    refine: bool,
    refine_maxfev: int,
    n_chains: int = 2,
) -> tuple[str | None, PackResult | None, RefineResult | None]:
    """Pack (and, if requested, refine) one candidate chain conformation.

    Top-level and picklable (importable as ``polyfind.pipeline._pack_and_refine_one``) so
    it can run in a worker process under Windows' spawn model: it takes only plain,
    picklable data -- no ``Candidate``/``RISModel``/``Polymer`` objects -- and rebuilds the
    ``RISModel`` and the ``PeriodicChain`` itself. Doing both stages for one candidate here
    means a candidate's refinement never has to wait for another candidate's packing.

    Returns ``(error, pack_result, refine_result)``; on a chain-build failure ``error`` is
    the message and the other two are ``None`` (mirrors the try/except the serial code used
    to have around ``periodic_chain``).
    """
    polymer = get_polymer(polymer_name)
    model = RISModel.from_dict(model_dict)
    try:
        chain = periodic_chain(polymer, seq, model.states)
    except (ValueError, RuntimeError) as e:
        return str(e), None, None
    rng = np.random.default_rng(seed)
    pack_res = pack(chain, n_chains=n_chains, n_random=n_random, n_refine=n_refine, rng=rng)[0]
    refine_res = refine_crystal(polymer, pack_res, maxfev=refine_maxfev) if refine else None
    return None, pack_res, refine_res


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

    # 3+4. packing and refinement -- candidates are independent, so each candidate's pack
    # and refine run together in one worker (process or, for workers<=1, this process),
    # so a candidate's refinement never waits on another candidate's packing.
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

    model_dict = model.to_dict()
    jobs = [
        dict(
            polymer_name=polymer.name,
            model_dict=model_dict,
            seq=cand.seq,
            # per-candidate seed, independent of worker count: workers=1 and workers=N
            # must produce byte-identical results.
            seed=cfg.seed + 1000 * i,
            n_random=cfg.n_random,
            n_refine=cfg.n_refine,
            refine=cfg.refine,
            refine_maxfev=cfg.refine_maxfev,
        )
        for i, cand in enumerate(to_pack)
    ]
    n_workers = cfg.workers
    if n_workers is None:
        # Physical cores, not logical, and not cpu_count-1.  Measured on a 6-physical /
        # 12-logical machine: 3 workers beat 5 by 1.24x, because each worker's packing
        # kernel is memory-bandwidth bound rather than core bound (the same effect that
        # made the table build saturate at ~4 workers; see lattice_table).  Over-
        # subscribing logical cores makes them contend for bandwidth they already lack.
        n_workers = min(len(to_pack), max(1, _physical_cores() // 2)) if to_pack else 1
    n_workers = max(1, n_workers)

    outcomes: list[tuple] = [None] * len(jobs)

    def log_outcome(i, err, pr, rr):
        cand = to_pack[i]
        if err is not None:
            log(f"      skip {cand.name}: {err}")
        else:
            log(f"      packed {pr.chain:<14s} E/mon={pr.energy_per_monomer:8.3f}  a={pr.a:.2f} b={pr.b:.2f} c={pr.c:.2f}")

    def run_serial():
        for i, job in enumerate(jobs):
            err, pr, rr = _pack_and_refine_one(**job)
            outcomes[i] = (err, pr, rr)
            log_outcome(i, err, pr, rr)

    if not jobs:
        pass
    elif n_workers <= 1:
        run_serial()
    else:
        try:
            with ProcessPoolExecutor(max_workers=n_workers) as ex:
                futs = {ex.submit(_pack_and_refine_one, **job): i for i, job in enumerate(jobs)}
                for fut in as_completed(futs):
                    i = futs[fut]
                    err, pr, rr = fut.result()
                    outcomes[i] = (err, pr, rr)
                    log_outcome(i, err, pr, rr)
        except (OSError, RuntimeError, NotImplementedError) as e:
            # e.g. a sandbox without process/spawn support, or a broken process pool:
            # fall back to serial rather than failing the whole pipeline.
            log(f"      warning: process pool unavailable ({e!r}); falling back to serial packing/refinement")
            outcomes = [None] * len(jobs)
            run_serial()

    # keep candidate order (not completion order) in the report
    packed = []
    refined = []
    for i in range(len(jobs)):
        outc = outcomes[i]
        if outc is None or outc[0] is not None:
            continue
        _, pr, rr = outc
        packed.append((to_pack[i], pr))
        if rr is not None:
            refined.append(rr)
    if cfg.refine:
        # replace packed energies by refined ones for the ranking
        packed = [(cand, r.result) for (cand, _), r in zip(packed, refined)]
    timings["pack"] = time.time() - t0
    log(f"[3/5] packed {len(packed)} chain conformations ({timings['pack']:.1f} s)")
    # packing and refinement are fused per candidate (see above), so both stages share
    # the same measured wall time.
    timings["refine"] = timings["pack"] if cfg.refine else 0.0
    log(f"[4/5] refinement done ({timings['refine']:.1f} s)")

    # 5. amorphous ensemble
    t0 = time.time()
    states, coords = sample_ensemble(polymer, model, cfg.n_bonds, cfg.n_chains, cfg.temperature, rng=rng)
    amorphous = ensemble_stats(polymer, model, states, coords, cfg.temperature)
    _, _, interface = lamella_interface(polymer, model, stem_bonds=12, tail_bonds=60, T=cfg.temperature, n_chains=min(cfg.n_chains, 2000), rng=rng)
    timings["amorphous"] = time.time() - t0
    log(f"[5/5] amorphous statistics ({timings['amorphous']:.2f} s)")
    return PipelineResult(polymer, model, fit, cands, packed, refined, amorphous, interface, timings)
