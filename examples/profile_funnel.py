"""Stage-by-stage CPU profile of the default PVDF funnel in one process, with cProfile
splits of every candidate's refinement (geometry vs kernel).  Mirrors PipelineConfig
defaults; no table cache.  Run: python examples/profile_funnel.py"""
import cProfile, os, pstats, time
import numpy as np
from polyfind.polymers import PVDF
from polyfind.forcefield import SimpleFF, fit_ris
from polyfind.enumerate import enumerate_periodic
from polyfind.pack import periodic_chain, CrystalPacker, default_bounds, polish, SCREEN_TABLE, _table_starts
from polyfind.lattice_table import pair_table, clear_pair_table_cache, shutdown_table_pool, default_table_procs
from polyfind.refine import refine_crystal
from polyfind.pipeline import PipelineConfig

os.environ.pop("POLYFIND_TABLE_CACHE", None)
cfg = PipelineConfig()
tic = time.perf_counter
T = {}
t = tic(); rep = fit_ris(PVDF, SimpleFF(), step=cfg.fit_step, n_monomers=cfg.fit_monomers, third_order=cfg.third_order); T["fit"] = tic() - t
model = rep.model
t = tic(); cands = enumerate_periodic(PVDF, model, max_period=cfg.max_period, k_per_period=cfg.k_per_period); T["enumerate"] = tic() - t
def n_atoms(c): return 3 * c.period * (c.helix.periods_per_repeat or 0)
packable = [c for c in cands if c.helix.c is not None and n_atoms(c) <= cfg.max_atoms_per_chain]
to_pack = list(packable[: cfg.top_k_pack])
for c in packable:
    if c.known_as and c not in to_pack: to_pack.append(c)
print(f"fit {T['fit']:.2f} s ({rep.n_evaluations} evals); enumerate {T['enumerate']:.2f} s ({len(cands)} cands); packing {[c.name for c in to_pack]}; table procs {default_table_procs()}", flush=True)
rows = []
prof_keys = ["build_backbone", "nerf", "_block_coords", "repeat_chains", "repeat_rotvec", "_close", "_jacobian", "_residual",
             "energy_and_grad", "kabsch", "pendant_positions", "_orient_block", "_align_about_z", "shape_batch", "value_and_grad_analytic", "svd", "minimize"]
for cand in to_pack:
    r = {"name": cand.name}
    t = tic(); chain = periodic_chain(PVDF, cand.seq, model.states); r["chain"] = tic() - t
    r["atoms"] = chain.n_atoms
    clear_pair_table_cache()
    t = tic(); tbl = pair_table(chain, cache_dir=None, memory=True, **SCREEN_TABLE); r["table"] = tic() - t
    packer = CrystalPacker(chain, n_chains=2)
    b = default_bounds(chain); keys = ["a", "b", "gamma", "phi1", "phi2", "dz"]
    lo = np.array([b[k][0] for k in keys]); hi = np.array([b[k][1] for k in keys])
    gammas = [lo[2]]
    t = tic(); starts = _table_starts(packer, chain, lo, hi, cfg.n_refine, [0, 1], 0.25, gammas, tbl, 8.0, 0.2, 1.0, None, SCREEN_TABLE, False); r["screen"] = tic() - t
    free = [i for i in range(6) if hi[i] > lo[i]]
    packer.n_energy_rows = 0
    t = tic(); results = [packer.result(polish(packer, p, lo, hi, free, 1500)) for p in starts]; r["polish"] = tic() - t
    r["polish_rows"] = packer.n_energy_rows; r["n_starts"] = len(starts)
    results.sort(key=lambda x: x.energy_per_cell); best = results[0]
    pr = cProfile.Profile()
    t = tic(); pr.enable(); rr = refine_crystal(PVDF, best, maxfev=cfg.refine_maxfev); pr.disable(); r["refine"] = tic() - t
    r["refine_evals"] = rr.n_evaluations; r["refine_vars"] = rr.n_variables; r["param"] = rr.parametrisation
    st = pstats.Stats(pr); st.sort_stats("cumulative")
    tot = st.total_tt
    split = {}
    for (file, line, fn), (cc, nc, tt, ct, callers) in st.stats.items():
        for k in prof_keys:
            if fn == k or fn.endswith(k):
                split[k] = max(split.get(k, 0.0), ct)
    r["prof_total"] = tot; r["split"] = {k: round(v, 3) for k, v in split.items()}
    r["E"] = rr.result.energy_per_monomer
    rows.append(r)
    print(f"{r['name']:<14s} atoms={r['atoms']:2d} chain {r['chain']*1e3:5.0f} ms | table {r['table']:5.2f} s | screen {r['screen']:5.2f} s | polish {r['polish']:5.2f} s ({r['polish_rows']} rows, {r['n_starts']} starts) | refine {r['refine']:5.2f} s ({r['refine_evals']} evals, {r['refine_vars']} vars, {r['param']}) E/mon {r['E']:.4f}", flush=True)
    print("   refine cProfile cumulative (s):", r["split"], flush=True)
shutdown_table_pool()
tot = {k: sum(r[k] for r in rows) for k in ("chain", "table", "screen", "polish", "refine")}
grand = T["fit"] + T["enumerate"] + sum(tot.values())
print("\nSTAGE TOTALS (single process):")
print(f"  fit        {T['fit']:6.2f} s  {100*T['fit']/grand:4.1f}%")
print(f"  enumerate  {T['enumerate']:6.2f} s  {100*T['enumerate']/grand:4.1f}%")
for k in ("chain", "table", "screen", "polish", "refine"):
    print(f"  {k:<10s} {tot[k]:6.2f} s  {100*tot[k]/grand:4.1f}%")
print(f"  total      {grand:6.2f} s")
