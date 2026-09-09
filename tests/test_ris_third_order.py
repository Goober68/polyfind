import numpy as np
import pytest

from polyfind.polymers import THREE_STATE
from polyfind.ris import RISModel, KB


def model3(seed=0, B=2):
    rng = np.random.default_rng(seed)
    e1 = rng.normal(size=(B, 3)) * 0.4
    e1[:, 0] = 0
    e2 = rng.normal(size=(B, 3, 3)) * 0.6
    e2[:, 0, 0] = 0
    e3 = rng.normal(size=(B, 3, 3, 3)) * 0.5
    return RISModel(THREE_STATE, B, e1, e2, e3, name="rand3")


@pytest.mark.parametrize("B", [1, 2])
@pytest.mark.parametrize("n", [2, 3, 5, 7])
def test_third_order_minimum_and_kbest_vs_brute_force(B, n):
    m = model3(B=B)
    brute = m.enumerate_all(n)
    e, s = m.minimum(n)
    assert e == pytest.approx(brute[0][0]) and m.energy(s) == pytest.approx(e)
    kb = m.k_best(n, 12)
    assert [round(x, 9) for x, _ in kb] == [round(x, 9) for x, _ in brute[:12]]
    for e, s in kb:
        assert m.energy(s) == pytest.approx(e)


@pytest.mark.parametrize("period", [2, 4, 6])
def test_third_order_cyclic_vs_brute_force(period):
    m = model3(B=2)
    ck = m.cyclic_k_best(period, 40)
    brute = m.enumerate_all(period, periodic=True)
    for e, s in ck:
        assert m.energy(s, periodic=True) == pytest.approx(e)
    assert ck[0][0] == pytest.approx(brute[0][0])
    reported = sorted(round(e, 9) for e, _ in ck)
    for e, _ in brute[:8]:
        assert round(e, 9) in reported


def test_third_order_partition_marginals_sampling():
    m = model3(B=2)
    n, T = 6, 350.0
    brute = m.enumerate_all(n)
    w = np.array([np.exp(-e / (KB * T)) for e, _ in brute])
    z = w.sum()
    assert m.log_partition(n, T) == pytest.approx(np.log(z), rel=1e-10)
    p = np.zeros((n, 3))
    for wi, (e, s) in zip(w, brute):
        p[np.arange(n), s] += wi / z
    assert np.allclose(m.marginals(n, T), p, atol=1e-10)
    samp = m.sample(n, T, 30000, rng=np.random.default_rng(0))
    freq = np.stack([(samp == s).mean(axis=0) for s in range(3)], axis=1)
    assert np.allclose(freq, p, atol=0.012)
    # clamped
    clamp = {1: 2, 4: 0}
    sel = [(e, s) for e, s in brute if s[1] == 2 and s[4] == 0]
    zc = sum(np.exp(-e / (KB * T)) for e, _ in sel)
    assert m.log_partition(n, T, clamp) == pytest.approx(np.log(zc), rel=1e-10)
    sc = m.sample(n, T, 200, rng=np.random.default_rng(1), clamp=clamp)
    assert (sc[:, 1] == 2).all() and (sc[:, 4] == 0).all()


def test_third_order_breaks_TG_degeneracy():
    """A pairwise model cannot tell TG+TG+ (helix) from TG+TG- (glide); a triplet term can."""
    e1 = np.zeros((2, 3))
    e2 = np.zeros((2, 3, 3))
    m2 = RISModel(THREE_STATE, 2, e1, e2)
    tgtg_p = m2.parse("TG+TG+")
    tgtg_m = m2.parse("TG+TG-")
    assert m2.energy(tgtg_p, periodic=True) == m2.energy(tgtg_m, periodic=True)
    e3 = np.zeros((2, 3, 3, 3))
    e3[1, 1, 0, 1] = e3[1, 2, 0, 2] = 0.7  # G T G same-sign penalty (type-1 triple)
    m3 = RISModel(THREE_STATE, 2, e1, e2, e3)
    assert m3.energy(tgtg_m, periodic=True) < m3.energy(tgtg_p, periodic=True)
    best = m3.cyclic_k_best(4, 5)[0][1]
    assert m3.energy(best, periodic=True) <= m3.energy(tgtg_m, periodic=True)


def test_infinite_chain_free_energy_third_order():
    m = model3(B=2)
    f = m.infinite_chain_free_energy_per_repeat(300.0)
    fa, fb = m.free_energy(200, 300.0), m.free_energy(400, 300.0)
    assert (fb - fa) / 100 == pytest.approx(f, rel=1e-6)


def test_json_roundtrip_third_order(tmp_path):
    m = model3()
    p = tmp_path / "m3.json"
    m.save(str(p))
    m2 = RISModel.load(str(p))
    assert np.allclose(m2.third_order, m.third_order)
