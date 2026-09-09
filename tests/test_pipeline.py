import numpy as np

from polyfind.pipeline import PipelineConfig, run_pipeline
from polyfind.ris import polyethylene_like_model


def test_pipeline_smoke_pe():
    cfg = PipelineConfig(polymer="pe", ris_model=polyethylene_like_model(), max_period=4, top_k_pack=1, pack_known=False,
                         n_random=200, n_refine=1, refine=True, refine_maxfev=100, n_chains=100, n_bonds=40, temperature=400.0)
    res = run_pipeline(cfg, verbose=False)
    assert res.candidates[0].name == "T"
    assert len(res.packed) == 1 and res.packed[0][1].chain == "T"
    assert res.amorphous is not None and 0 < res.amorphous.state_fractions["T"] < 1
    rep = res.report()
    assert "Crystal packing" in rep and "Amorphous ensemble" in rep and "experiment" in rep


def test_pipeline_workers_deterministic():
    """Independent candidates run in parallel processes (stages 3+4); workers=1 (serial)
    and workers=2 (ProcessPoolExecutor) must give byte-identical packing and refinement
    results, since each candidate's RNG seed is derived from its index, not the number
    of workers."""
    def make_cfg(workers):
        return PipelineConfig(
            polymer="pe", ris_model=polyethylene_like_model(), max_period=2, top_k_pack=2, pack_known=False,
            n_random=50, n_refine=1, refine=True, refine_maxfev=50, n_chains=50, n_bonds=20, temperature=400.0,
            seed=0, workers=workers,
        )

    res_serial = run_pipeline(make_cfg(1), verbose=False)
    res_parallel = run_pipeline(make_cfg(2), verbose=False)

    assert len(res_serial.packed) == 2 and len(res_parallel.packed) == 2

    names_serial = [r.chain for _, r in res_serial.packed]
    names_parallel = [r.chain for _, r in res_parallel.packed]
    assert names_serial == names_parallel == ["T", "TG+"]

    for (_, p1), (_, p2) in zip(res_serial.packed, res_parallel.packed):
        assert p1.energy_per_monomer == p2.energy_per_monomer
        assert (p1.a, p1.b, p1.c, p1.gamma) == (p2.a, p2.b, p2.c, p2.gamma)

    assert len(res_serial.refined) == len(res_parallel.refined) == 2
    for r1, r2 in zip(res_serial.refined, res_parallel.refined):
        assert np.array_equal(r1.torsions, r2.torsions)
        assert r1.result.energy_per_monomer == r2.result.energy_per_monomer
        assert (r1.result.a, r1.result.b, r1.result.c) == (r2.result.a, r2.result.b, r2.result.c)


def test_cli_enumerate_runs(capsys):
    from polyfind.cli import main
    import json, tempfile, os
    m = polyethylene_like_model()
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "m.json")
        m.save(p)
        main(["enumerate", "--polymer", "pe", "--model", p, "--max-period", "4", "--top", "5"])
    out = capsys.readouterr().out
    assert "planar zigzag" in out
