import warnings

import numpy as np
import pytest

from polyfind.chain import build_chain, build_chain_batch
from polyfind.forcefield import (
    SimpleFF,
    _reversal_angle_residual,
    _reversal_image3,
    _reversal_images,
    erfc_approx,
    fit_ris,
)
from polyfind.polymers import PE, PVDF, THREE_STATE, get_polymer
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


# ------------------------------------------------ reflection-with-reversal, orders 1-3
def _raw_fit(name, third_order=True, step=20.0):
    """An unsymmetrised fit, so the relations can be measured rather than assumed."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # the chiral polymers warn once on build_chain
        return fit_ris(get_polymer(name), SimpleFF(), step=step, n_monomers=3,
                       third_order=third_order, symmetrize=False)


@pytest.mark.parametrize("name", ["pvdf", "pvdc", "cfe", "cdfe"])
def test_third_order_terms_obey_reflection_with_reversal(name):
    """The third-order image of E(phi_1..phi_N) == E(-phi_N..-phi_1), measured.

    Reversal reverses the order of a triple as well as mirroring its states, and shifts
    the bond type by one more than the pair term does, so

        e3[b, x, y, z] == e3[(c - 2 - b) % B, m(z), m(y), m(x)].

    PVDF and PVDC are achiral controls: any true symmetry of the model has to hold there
    too, and does.  The near misses -- mirroring without reversing the triple, reversing
    without mirroring, and the right form at the wrong bond-type shift -- are all wrong
    for the chiral pair, which is what makes this a measurement and not a relabelling.
    """
    p = get_polymer(name)
    rep = _raw_fit(name)
    e1, e2, e3 = rep.model.first_order, rep.model.second_order, rep.model.third_order
    B, m = p.bonds_per_repeat, np.array(THREE_STATE.mirror)
    *_, c = _reversal_images(e1, e2, THREE_STATE.mirror, B)
    assert c == 1  # for B = 2 the shift is decided by the first- and second-order fit

    _, resid = _reversal_image3(e3, THREE_STATE.mirror, B, c)
    assert resid < 0.5, resid  # grid noise: 0.00 (pvdf, cfe, cdfe), 0.02 (pvdc)

    plain = float(np.abs(e3 - e3[:, m][:, :, m][:, :, :, m]).max())
    wrong_shift = _reversal_image3(e3, THREE_STATE.mirror, B, c + 1)[1]
    no_reverse = float(np.abs(e3 - e3[[(c - 2 - b) % B for b in range(B)]][:, m][:, :, m][:, :, :, m]).max())
    no_mirror = float(np.abs(e3 - np.transpose(e3[[(c - 2 - b) % B for b in range(B)]], (0, 3, 2, 1))).max())
    if p.is_chiral:
        assert plain > 5.0 and wrong_shift > 5.0 and no_reverse > 5.0 and no_mirror > 5.0
    else:
        assert plain < 0.5  # an achiral chain has the plain mirror as well


@pytest.mark.parametrize("name", ["cfe", "cdfe"])
def test_a_chiral_fit_symmetrises_all_three_orders(name):
    """The default fit averages a chiral model over the operation it just validated."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        model = fit_ris(get_polymer(name), SimpleFF(), step=20.0, n_monomers=3, third_order=True).model
    B, m = model.B, np.array(model.states.mirror)
    *_, c = _reversal_images(model.first_order, model.second_order, model.states.mirror, B)
    assert _reversal_image3(model.third_order, model.states.mirror, B, c)[1] < 1e-9
    # and the real G+/G- asymmetry is still there, in the third order too
    assert np.abs(model.third_order - model.third_order[:, m][:, :, m][:, :, :, m]).max() > 5.0


@pytest.mark.parametrize("name", ["pvdf", "pvdc", "cfe", "cdfe"])
def test_state_angles_obey_reflection_with_reversal(name):
    """``adapt_angles``' mirror pairing is valid for a chiral fit too, and is checked.

    The 1-D scan of a bond of type ``b`` is the scan of type ``(c - b) % B`` read at the
    negated angle, so ``arg1[b, s] == -arg1[(c - b) % B, m(s)]``.  Averaging that over the
    bond types (a bijection) is what makes the adapted angles of a mirror pair equal and
    opposite -- for a chiral chain as much as an achiral one, even though the *per bond
    type* minima are wildly asymmetric there (80 deg apart for CFE and CDFE).
    """
    p = get_polymer(name)
    rep = _raw_fit(name, third_order=False)
    B = p.bonds_per_repeat
    *_, c = _reversal_images(rep.model.first_order, rep.model.second_order, THREE_STATE.mirror, B)
    assert _reversal_angle_residual(rep.argmin1, THREE_STATE, B, c) == 0.0  # on a dense grid, exactly

    gp, gm = THREE_STATE.names[1], THREE_STATE.names[2]
    per_type = max(abs(rep.argmin1[b][gp] + rep.argmin1[b][gm]) for b in range(B))
    assert (per_type > 10.0) if p.is_chiral else (per_type == 0.0)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        model = fit_ris(get_polymer(name), SimpleFF(), step=20.0, n_monomers=3).model
    assert model.states.angles[1] == -model.states.angles[2]  # aggregated: exactly paired


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
    """The adaptive (coarse-to-fine) fit must reproduce the dense (default) 10-degree
    scan's first- and second-order energies: never meaningfully worse, often better since
    it is not confined to a grid, with the same mirror symmetry, the same basin-angle
    signs, and several times fewer calculator evaluations.

    One matrix entry is excluded from the ``second_order`` comparison here: PVDF's
    CF2-centred G+/G- pair (bond type 1), whose basin has no interior minimum -- see
    ``test_adaptive_pvdf_g_plus_g_minus_has_no_interior_minimum``, which checks it
    explicitly instead. Everywhere else -- including PE's own G+/G- pair (whose basin
    *does* have an interior minimum, just close to the basin edge) and PVDF's CH2-centred
    G+/G+ pair (a real, fairly sharp bowl around (86, 86) deg that the dense grid samples
    only at (90, 80), verified by a 2-degree brute-force scan) -- adaptive matches dense
    within [-0.6, +0.1] kcal/mol: it is never meaningfully worse, and where it is better
    by more than refinement noise, that reflects genuine sub-grid curvature dense's grid
    is too coarse to resolve, not an error.
    """
    dense = fit_ris(polymer, SimpleFF(), step=10.0, n_monomers=n_monomers, third_order=True, scan="dense")
    adap = fit_ris(polymer, SimpleFF(), step=10.0, n_monomers=n_monomers, third_order=True, scan="adaptive")

    d1, a1 = dense.model.first_order, adap.model.first_order
    d2, a2 = dense.model.second_order, adap.model.second_order

    mask = np.ones_like(d2, dtype=bool)
    if polymer is PVDF:
        gp, gm = THREE_STATE.index("G+"), THREE_STATE.index("G-")
        mask[1, gp, gm] = mask[1, gm, gp] = False

    assert np.all(a1 <= d1 + 0.1)
    assert np.all(a1 >= d1 - 0.5)
    assert np.all(a2[mask] <= d2[mask] + 0.1)
    assert np.all(a2[mask] >= d2[mask] - 0.6)

    mir = np.array(polymer.states.mirror)
    assert np.allclose(a1, a1[:, mir], atol=1e-8)
    assert np.allclose(a2, a2[:, mir][:, :, mir], atol=1e-8)

    for b in range(polymer.bonds_per_repeat):
        assert np.sign(adap.argmin1[b]["G+"]) == np.sign(dense.argmin1[b]["G+"])
        assert np.sign(adap.argmin1[b]["G-"]) == np.sign(dense.argmin1[b]["G-"])

    assert dense.n_evaluations >= 2.5 * adap.n_evaluations


def test_adaptive_pvdf_g_plus_g_minus_has_no_interior_minimum():
    """PVDF's CF2-centred (bond type 1) G+/G- pair is the one place adaptive and dense
    disagree by more than grid/refinement noise, and it is expected, not a bug: this
    basin has no interior critical point near the G+/T boundary. Scanning phi from 100 to
    140 deg (crossing straight through the nominal boundary at 120) at each point's own
    best psi gives a smooth, monotonically *decreasing* profile (E ~= 7.8, 5.1, 4.88 (at
    phi=120), 4.7, 4.04 (at phi=140), ...; no kink, no minimum). So the true infimum of a
    correctly bounded search of the open G+/G- basin (matching ``_basins``' tie-break --
    see ``_basin_bounds``) sits at its own edge, at ~4.85 kcal/mol; excluding the single
    boundary point cannot raise the infimum of a continuous, monotonic function
    approaching it. The dense step=10 deg scan reports a much higher ~8.17 kcal/mol only
    because it is too coarse to sample anywhere near that edge (its neighbouring grid
    points are 110 deg, E ~= 10.2, and 120 deg itself, excluded by the same tie-break).
    This is asserted explicitly, as expected behaviour, rather than folded into a widened
    tolerance in the general accuracy test above.
    """
    dense = fit_ris(PVDF, SimpleFF(), step=10.0, n_monomers=6, scan="dense")
    adap = fit_ris(PVDF, SimpleFF(), step=10.0, n_monomers=6, scan="adaptive")
    gp, gm = THREE_STATE.index("G+"), THREE_STATE.index("G-")

    d = dense.model.second_order[1, gp, gm]
    a = adap.model.second_order[1, gp, gm]
    assert d == pytest.approx(8.17, abs=0.2)
    assert a == pytest.approx(4.85, abs=0.2)
    assert a - d < -3.0  # a large, expected improvement -- not refinement noise

    # neither coordinate sits exactly on the G+/T (or G-/T) tie boundary that _basins
    # would reassign to T, but adaptive does ride right up against one of them
    phi, psi = adap.argmin2[1][("G+", "G-")]
    assert phi != 120.0 and psi != -120.0
    assert min(abs(phi - 120.0), abs(psi + 120.0)) < 1e-3
