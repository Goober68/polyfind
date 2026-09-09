import numpy as np

from polyfind.enumerate import enumerate_periodic, table
from polyfind.polymers import PE, PVDF
from polyfind.ris import polyethylene_like_model, RISModel
from polyfind.polymers import THREE_STATE


def test_pe_textbook_model_ranks_all_trans_first():
    m = polyethylene_like_model()
    cands = enumerate_periodic(PE, m, max_period=6)
    assert cands[0].name == "T" and cands[0].known_as is not None
    names = [c.name for c in cands]
    assert len(names) == len(set(names))
    # the TG helix should be there and cost 0.5 kcal/mol per monomer (period 2 -> 0.25 per bond... per monomer = per bond for PE)
    tg = next(c for c in cands if c.name in ("TG+", "TG-"))
    assert abs(tg.energy_per_monomer - 0.25) < 1e-9
    assert tg.helix.periods_per_repeat == 3


def test_pvdf_enumeration_contains_polymorph_chains_once():
    rng = np.random.default_rng(0)
    e1 = np.abs(rng.normal(size=(2, 3))) * 0.3
    e1[:, 0] = 0
    e2 = rng.normal(size=(2, 3, 3)) * 0.5
    e2[:, 0, 0] = 0
    mir = np.array(THREE_STATE.mirror)
    e1 = 0.5 * (e1 + e1[:, mir])
    e2 = 0.5 * (e2 + e2[:, mir][:, :, mir])
    m = RISModel(THREE_STATE, 2, e1, e2)
    cands = enumerate_periodic(PVDF, m, max_period=8, k_per_period=80)
    known = [c.known_as for c in cands if c.known_as]
    assert "beta (TTTT)" in known and "alpha/delta (TGTG')" in known and "gamma/epsilon (T3GT3G')" in known
    assert len(known) == len(set(known))
    # all candidates propagate and are commensurate
    assert all(c.helix.rise_per_bond >= 0.3 for c in cands)
    # energies are per monomer and sorted
    ems = [c.energy_per_monomer for c in cands]
    assert ems == sorted(ems)
    assert "beta" in table(cands)
