import numpy as np
import pytest

from polyfind.chain import build_chain, build_chain_batch
from polyfind.forcefield import SimpleFF, fit_ris, erfc_approx
from polyfind.polymers import PE, PVDF, THREE_STATE
from scipy.special import erfc


def test_erfc_approx():
    x = np.linspace(-3, 3, 61)
    assert np.allclose(erfc_approx(x), erfc(x), atol=2e-7)


def test_energy_batch_matches_single():
    ff = SimpleFF()
    rng = np.random.default_rng(0)
    structs = [build_chain(PVDF, rng.uniform(-180, 180, 8)) for _ in range(6)]
    eb = ff.energy_batch(structs)
    for s, e in zip(structs, eb):
        assert ff.energy(s) == pytest.approx(e, rel=1e-8)


def test_pe_trans_below_gauche_and_pentane_effect():
    ff = SimpleFF()
    rep = fit_ris(PE, ff, step=15.0, n_monomers=8, scan="dense")
    m = rep.model
    T, GP, GM = 0, 1, 2
    assert abs(m.first_order[0, T]) < 1e-9
    assert 0.2 < m.first_order[0, GP] < 2.0  # gauche costs a fraction of a kcal/mol
    assert m.first_order[0, GP] == pytest.approx(m.first_order[0, GM])  # mirror symmetry
    assert m.second_order[0, GP, GM] > 1.0  # pentane effect: G+G- strongly disfavoured
    assert abs(m.second_order[0, T, T]) < 1e-9
    assert rep.n_evaluations == 1 + 24 + 24 * 24


def test_pvdf_fit_has_two_pair_types_and_is_mirror_symmetric():
    ff = SimpleFF()
    rep = fit_ris(PVDF, ff, step=20.0, n_monomers=5)
    m = rep.model
    assert m.B == 2 and m.second_order.shape == (2, 3, 3)
    mir = np.array(THREE_STATE.mirror)
    assert np.allclose(m.second_order, m.second_order[:, mir][:, :, mir])
    # the two pair matrices differ (CH2- vs CF2-centred)
    assert not np.allclose(m.second_order[0], m.second_order[1])


def test_energy_coords_matches_energy_batch_list_and_single():
    ff = SimpleFF()
    rng = np.random.default_rng(3)
    dih = rng.uniform(-180, 180, size=(5, 8))
    template, coords = build_chain_batch(PVDF, dih)
    e_tuple = ff.energy_batch((template, coords))
    structs = [build_chain(PVDF, dih[m]) for m in range(5)]
    e_list = ff.energy_batch(structs)
    np.testing.assert_allclose(e_tuple, e_list, atol=1e-8)
    for m, e in enumerate(e_tuple):
        assert ff.energy(structs[m]) == pytest.approx(e, rel=1e-8)


def test_fit_ris_unknown_scan_raises():
    with pytest.raises(ValueError):
        fit_ris(PE, SimpleFF(), scan="bogus")


def test_dense_scan_batched_matches_unbatched_reference():
    """The batched dense scan (build_chain_batch + energy_batch) must reproduce the
    original per-point build_chain/calc.energy_batch(list) path exactly."""
    ff = SimpleFF()
    rep = fit_ris(PVDF, ff, step=20.0, n_monomers=5, scan="dense")
    grid = rep.grid
    N = max(5 * 2, 2 * 2 + 6)
    base = np.full(N, 180.0)
    ref = ff.energy(build_chain(PVDF, base))
    j0 = (N // 2 - 1) - ((N // 2 - 1) - 0) % 2
    structs = []
    for phi in grid:
        d = base.copy()
        d[j0] = phi
        structs.append(build_chain(PVDF, d))
    expected = ff.energy_batch(structs) - ref
    np.testing.assert_allclose(rep.scan1[0], expected, atol=1e-9)


@pytest.mark.parametrize("polymer,n_monomers", [(PE, 8), (PVDF, 6)])
def test_adaptive_matches_dense_10deg_fit(polymer, n_monomers):
    """The adaptive (coarse-to-fine) fit must reproduce the dense 10-degree scan's first-
    and second-order energies: never meaningfully worse, often better since it is not
    confined to a grid, with the same mirror symmetry, the same basin-angle signs, and
    several times fewer calculator evaluations.

    The lower bound on second_order is wider than first_order's, for one documented case:
    PVDF's CF2-centred G+/G- pair (bond type 1). Its energy has no interior critical point
    near the G+/T basin boundary -- scanning phi from 100 to 140 deg (crossing straight
    through the boundary at 120) at each point's own best psi gives a smooth, monotonically
    *decreasing* profile (E ~= 7.8, 5.1, 4.88 (at phi=120), 4.7, 4.04 (at phi=140), ...; no
    kink, no minimum). So the true infimum of a *correctly* bounded search of the G+/G-
    basin (open at 120, matching ``_basins``' tie-break -- see ``_basin_bounds``) sits
    right at that basin edge, at ~4.85 kcal/mol: excluding the single boundary point cannot
    raise the infimum of a continuous, monotonic function approaching it. The dense scan
    reports a much higher 8.17 kcal/mol only because step=10 deg is too coarse to sample
    anywhere near that edge (its adjacent grid points are 110 deg, E ~= 10.2, and 120 deg
    itself, excluded by the same tie-break). Adaptive is not confined to that resolution
    limit and correctly finds the lower, edge-of-basin value; everywhere else the two agree
    closely.
    """
    dense = fit_ris(polymer, SimpleFF(), step=10.0, n_monomers=n_monomers, third_order=True, scan="dense")
    adap = fit_ris(polymer, SimpleFF(), step=10.0, n_monomers=n_monomers, third_order=True, scan="adaptive")

    d1, a1 = dense.model.first_order, adap.model.first_order
    d2, a2 = dense.model.second_order, adap.model.second_order

    assert np.all(a1 <= d1 + 0.1)
    assert np.all(a1 >= d1 - 0.5)
    assert np.all(a2 <= d2 + 0.1)
    assert np.all(a2 >= d2 - 5.0)

    mir = np.array(polymer.states.mirror)
    assert np.allclose(a1, a1[:, mir], atol=1e-8)
    assert np.allclose(a2, a2[:, mir][:, :, mir], atol=1e-8)

    for b in range(polymer.bonds_per_repeat):
        assert np.sign(adap.argmin1[b]["G+"]) == np.sign(dense.argmin1[b]["G+"])
        assert np.sign(adap.argmin1[b]["G-"]) == np.sign(dense.argmin1[b]["G-"])

    assert dense.n_evaluations >= 2.5 * adap.n_evaluations
