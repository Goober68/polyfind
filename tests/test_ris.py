import numpy as np
import pytest

from polyfind.ris import RISModel, polyethylene_like_model, KB
from polyfind.polymers import THREE_STATE


def pvdf_like_model():
    rng = np.random.default_rng(0)
    e1 = rng.normal(size=(2, 3)) * 0.5
    e1[:, 0] = 0.0
    e2 = rng.normal(size=(2, 3, 3))
    e2[:, 0, 0] = 0.0
    return RISModel(THREE_STATE, 2, e1, e2, name="random-2type")


@pytest.mark.parametrize("model", [polyethylene_like_model(), pvdf_like_model()])
@pytest.mark.parametrize("n", [1, 2, 5, 8])
def test_viterbi_matches_brute_force(model, n):
    e_min, seq = model.minimum(n)
    brute = model.enumerate_all(n)
    assert e_min == pytest.approx(brute[0][0])
    assert model.energy(seq) == pytest.approx(e_min)


@pytest.mark.parametrize("model", [polyethylene_like_model(), pvdf_like_model()])
def test_k_best_matches_brute_force(model):
    n, k = 7, 25
    kb = model.k_best(n, k)
    brute = model.enumerate_all(n)
    assert [round(e, 9) for e, _ in kb] == [round(e, 9) for e, _ in brute[:k]]
    for e, s in kb:
        assert model.energy(s) == pytest.approx(e)


@pytest.mark.parametrize("model", [polyethylene_like_model(), pvdf_like_model()])
@pytest.mark.parametrize("period", [2, 4, 6])
def test_cyclic_k_best_matches_brute_force(model, period):
    ck = model.cyclic_k_best(period, k=30)
    brute = model.enumerate_all(period, periodic=True)
    # every reported cycle has the right energy, and the best energy matches brute force
    for e, s in ck:
        assert model.energy(s, periodic=True) == pytest.approx(e)
    assert ck[0][0] == pytest.approx(brute[0][0])
    # the first ten brute-force energies are all present among the reported cycles
    reported = sorted(round(e, 9) for e, _ in ck)
    for e, _ in brute[:10]:
        assert round(e, 9) in reported


@pytest.mark.parametrize("model", [polyethylene_like_model(), pvdf_like_model()])
def test_partition_function_matches_brute_force(model):
    n, T = 6, 300.0
    brute = model.enumerate_all(n)
    z = sum(np.exp(-e / (KB * T)) for e, _ in brute)
    assert model.log_partition(n, T) == pytest.approx(np.log(z), rel=1e-10)
    # marginals from brute force
    p = np.zeros((n, 3))
    for e, s in brute:
        p[np.arange(n), s] += np.exp(-e / (KB * T)) / z
    assert np.allclose(model.marginals(n, T), p, atol=1e-10)


def test_clamped_partition_and_marginals():
    model = pvdf_like_model()
    n, T = 6, 300.0
    clamp = {2: 1, 3: 0}
    brute = [(e, s) for e, s in model.enumerate_all(n) if s[2] == 1 and s[3] == 0]
    z = sum(np.exp(-e / (KB * T)) for e, _ in brute)
    assert model.log_partition(n, T, clamp) == pytest.approx(np.log(z), rel=1e-10)
    m = model.marginals(n, T, clamp)
    assert m[2, 1] == pytest.approx(1.0) and m[3, 0] == pytest.approx(1.0)


def test_sampling_frequencies_match_marginals():
    model = pvdf_like_model()
    n, T = 8, 400.0
    rng = np.random.default_rng(1)
    samp = model.sample(n, T, 40000, rng=rng)
    freq = np.stack([(samp == s).mean(axis=0) for s in range(3)], axis=1)
    assert np.allclose(freq, model.marginals(n, T), atol=0.01)
    # pair statistics too (joint of bonds 3,4) vs brute force
    brute = model.enumerate_all(n)
    z = sum(np.exp(-e / (KB * T)) for e, _ in brute)
    joint = np.zeros((3, 3))
    for e, s in brute:
        joint[s[3], s[4]] += np.exp(-e / (KB * T)) / z
    emp = np.zeros((3, 3))
    for a in range(3):
        for b in range(3):
            emp[a, b] = np.mean((samp[:, 3] == a) & (samp[:, 4] == b))
    assert np.allclose(emp, joint, atol=0.01)


def test_clamped_sampling_respects_clamp():
    model = pvdf_like_model()
    samp = model.sample(10, 300.0, 500, rng=np.random.default_rng(2), clamp={0: 2, 5: 1, 9: 0})
    assert (samp[:, 0] == 2).all() and (samp[:, 5] == 1).all() and (samp[:, 9] == 0).all()


def test_infinite_chain_free_energy_limit():
    model = polyethylene_like_model()
    T = 300.0
    f_inf = model.infinite_chain_free_energy_per_repeat(T)
    f_n = [model.free_energy(n, T) for n in (200, 400)]
    # increments converge to the per-repeat value
    assert (f_n[1] - f_n[0]) / 200 == pytest.approx(f_inf, rel=1e-6)


def test_parse_and_names_roundtrip():
    model = pvdf_like_model()
    seq = model.parse("TG+TG-")
    assert seq == [0, 1, 0, 2]
    assert model.names(seq) == "TG+TG-"


def test_json_roundtrip(tmp_path):
    model = pvdf_like_model()
    p = tmp_path / "m.json"
    model.save(str(p))
    m2 = RISModel.load(str(p))
    assert np.allclose(m2.first_order, model.first_order)
    assert np.allclose(m2.second_order, model.second_order)
    assert m2.states == model.states
