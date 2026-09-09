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
