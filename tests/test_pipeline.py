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


def _a_packed_cell(chain):
    from polyfind.pack import CrystalPacker

    return CrystalPacker(chain, n_chains=2).result(
        np.array([5.0, 9.0, 87.0, 10.0, 20.0, 0.3 * chain.c, 0.0]))


def test_reference_cell_angles_record_what_the_model_cannot_express():
    """gamma-PVDF's reference cell is monoclinic and this parametrisation is not.

    ``pack`` builds a = (a, 0, 0), b = (b cos g, b sin g, 0), c = (0, 0, c) with the chain
    axis along z, so the crystallographic *gamma* (a to b) is a search variable while
    *alpha* and *beta* are 90 deg by construction.  gamma-PVDF's unique angle is
    beta = 93 deg, between a and the chain axis, so it is precisely the one that cannot be
    represented; the package approximates that cell as orthorhombic and
    ``REFERENCE_CELL_ANGLES`` records both the real angles and that limitation.  This test
    pins all three halves of the statement: the record, the pin at 90 deg, and the one
    angle that really is free.
    """
    from polyfind.pack import default_bounds, periodic_chain, to_cif
    from polyfind.pipeline import EXPERIMENTAL_CELLS, REFERENCE_CELL_ANGLES
    from polyfind.polymers import PVDF, THREE_STATE

    assert REFERENCE_CELL_ANGLES.keys() == EXPERIMENTAL_CELLS.keys()
    for name, cells in EXPERIMENTAL_CELLS.items():
        assert REFERENCE_CELL_ANGLES[name].keys() == cells.keys()
    non_right = {(p, lab) for p, d in REFERENCE_CELL_ANGLES.items()
                 for lab, angs in d.items() if any(a != 90.0 for a in angs)}
    assert non_right == {("pvdf", "gamma/epsilon (T3GT3G')")}
    assert REFERENCE_CELL_ANGLES["pvdf"]["gamma/epsilon (T3GT3G')"] == (90.0, 93.0, 90.0)

    ch = periodic_chain(PVDF, [0, 0], THREE_STATE)
    assert default_bounds(ch)["gamma"] == (90.0, 90.0)  # pinned unless asked
    assert default_bounds(ch, gamma_free=True)["gamma"] == (60.0, 120.0)
    # ... while the other two angles are not variables at all: the CIF writer emits them as
    # literal 90, which is the parametrisation showing through.
    cif = to_cif(_a_packed_cell(ch))
    assert "_cell_angle_alpha 90\n" in cif and "_cell_angle_beta 90\n" in cif
    assert "_cell_angle_gamma 87.000" in cif  # the one angle that is a real variable
