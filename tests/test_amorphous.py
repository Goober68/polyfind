import numpy as np

from polyfind.amorphous import sample_ensemble, ensemble_stats, run_lengths, lamella_interface
from polyfind.polymers import PE, PVDF
from polyfind.ris import polyethylene_like_model


def test_run_lengths():
    s = np.array([[0, 0, 1, 0, 0, 0, 2, 0]])
    assert sorted(run_lengths(s, 0).tolist()) == [1, 2, 3]


def test_pe_ensemble_stats_and_characteristic_ratio():
    m = polyethylene_like_model(e_gauche=0.5, e_pentane=2.0)
    states, coords = sample_ensemble(PE, m, n_bonds=200, n_chains=400, T=400.0, rng=np.random.default_rng(0))
    st = ensemble_stats(PE, m, states, coords, 400.0)
    assert 0.5 < st.state_fractions["T"] < 0.8
    assert abs(st.state_fractions["G+"] - st.state_fractions["G-"]) < 0.05
    # textbook PE characteristic ratio is ~6-7 at ~400 K for this class of model
    assert 4.0 < st.characteristic_ratio < 9.0
    assert st.diad_fractions["TT"] > st.diad_fractions["G+G-"]
    assert 0 < st.tetrad_TTTT_fraction < st.state_fractions["T"]
    assert "characteristic ratio" in st.summary()


def test_lamella_interface_clamps_stem():
    m = polyethylene_like_model()
    tail, coords, stats = lamella_interface(PE, m, stem_bonds=10, tail_bonds=30, T=400.0, n_chains=300, rng=np.random.default_rng(1))
    assert tail.shape == (300, 30)
    assert 0 <= stats["p_turnback"] <= 1
    assert stats["mean_advance_along_stem_axis"] > 0  # on average the tail continues forward
