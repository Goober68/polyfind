"""The intrinsic electrostriction by two routes, and the two-state switch."""
import numpy as np
import pytest

from polyfind import electrostriction as ES
from polyfind import mechanics as M
from polyfind.polymers import PE, PVDF, THREE_STATE

T, GP, GM = 0, 1, 2


# Force-shifted Lennard-Jones: the energy-shifted default has a force jump at the cutoff, which
# beta's b axis straddles, and the polarization route lands its stress-targeted states on
# either side of it (see the electrostriction() docstring for the measured effect).
@pytest.fixture(scope="module")
def beta():
    ref, _ = M.refined_reference(PVDF, [T, T], label="beta", lj_cutoff="force")
    ref = M.relax_reference(ref)
    el = M.elastic_constants(ref)
    return ref, el


@pytest.fixture(scope="module")
def pe():
    ref, _ = M.refined_reference(PE, [T], label="PE", lj_cutoff="force")
    ref = M.relax_reference(ref)
    el = M.elastic_constants(ref)
    return ref, el


@pytest.fixture(scope="module")
def beta_es(beta):
    ref, el = beta
    return ES.electrostriction(ref, el, field=0.005)


def test_polyethylene_with_no_charges_has_no_electrostriction(pe):
    ref, el = pe
    es = ES.electrostriction(ref, el, field=0.005)
    idx = list(es.reachable)
    assert np.all(np.abs(es.M_strain[:, idx]) < 1e-30)
    assert np.all(np.abs(es.M_polarization[:, idx]) < 1e-30)
    assert np.all(np.abs(es.linear[:, idx]) < 1e-12)
    assert np.all(np.isnan(es.Q[:, idx]))
    assert np.allclose(es.chi_free, 0.0) and np.allclose(es.chi_clamped, 0.0)


def test_the_two_routes_agree_on_rigid_beta(beta_es):
    es = beta_es
    idx = list(es.reachable)
    a, b = es.M_strain[:, idx], es.M_polarization[:, idx]
    assert np.isfinite(a).all() and np.isfinite(b).all()
    assert float(np.max(np.abs(a))) > 1e-22, "beta-PVDF with fixed charges still turns in a field"
    assert es.max_relative_disagreement < 0.01
    # the free and clamped susceptibilities differ by the piezoelectric term, which is not small here
    assert es.chi_free[1] > es.chi_clamped[1] > 0.0


def test_the_even_part_is_quadratic_in_the_field(beta):
    ref, el = beta
    es1 = ES.electrostriction(ref, el, field=0.005, polarization_route=False)
    es2 = ES.electrostriction(ref, el, field=0.0025, polarization_route=False)
    idx = list(es1.reachable)
    scale = float(np.max(np.abs(es1.M_strain[:, idx])))
    assert np.max(np.abs(es1.M_strain[:, idx] - es2.M_strain[:, idx])) < 0.03 * scale


def test_the_odd_part_is_the_converse_piezoelectric_coefficient(beta, beta_es):
    ref, el = beta
    pz = M.piezoelectric(ref, el, field=0.005)
    d = pz.d_direct  # (3, n) pC/N == pm/V, columns in elastic.reachable order
    lin = beta_es.linear[:, list(beta_es.reachable)]
    assert np.isfinite(d).all()
    assert np.max(np.abs(lin - d)) < 0.02 * max(float(np.max(np.abs(d))), 1e-6)


def test_q_is_m_over_the_squared_susceptibility(beta_es):
    es = beta_es
    for i in range(3):
        if es.chi_free[i] > 1e-9:
            J = es.reachable[0]
            assert es.Q[i, J] == pytest.approx(es.M_strain[i, J] / (ES.EPS0 * es.chi_free[i]) ** 2, rel=1e-12)
        else:
            assert np.isnan(es.Q[i, es.reachable[0]])


def test_threshold_and_two_state_fraction():
    # 1.08 kcal/mol per monomer (the alpha-beta margin of the DFT fit) against 0.54 e.A
    e_star = ES.threshold_field(1.08, 0.54)
    assert e_star == pytest.approx(1.08 / ES.KCAL_PER_EV / 0.54)
    assert e_star * 1e4 > 500  # MV/m: a perfect PVDF crystal does not switch at working fields
    ts = ES.TwoState(delta_e=1.08, delta_mu=0.54, strain=0.11, n=1)
    assert ts.threshold == pytest.approx(e_star)
    assert ts.fraction(e_star, 300.0) == pytest.approx(0.5)
    # a lone monomer 47 meV up is 14% populated at 300 K; a 20-monomer cooperative unit is not
    assert 0.1 < ts.fraction(0.0, 300.0) < 0.2
    assert ES.TwoState(delta_e=1.08, delta_mu=0.54, strain=0.11, n=20).fraction(0.0, 300.0) < 1e-12
    assert ts.fraction(10.0 * e_star, 300.0) == pytest.approx(1.0)
    assert ts.fraction(0.0, 0.0) == 0.0 and ts.fraction(2 * e_star, 0.0) == 1.0
    assert ts.strain_at(10.0 * e_star) == pytest.approx(0.11)
    # a cooperative unit sharpens the switch about the same threshold
    coop = ES.TwoState(delta_e=1.08, delta_mu=0.54, strain=0.11, n=10)
    assert coop.fraction(e_star) == pytest.approx(0.5)
    assert coop.slope(e_star) > 5 * ts.slope(e_star)
    # a near-degenerate pair switches at working fields: 0.05 kcal/mol over 0.54 e.A
    assert ES.threshold_field(0.05, 0.54) * 1e4 == pytest.approx(40.1, rel=0.01)
    with pytest.raises(ValueError):
        ES.threshold_field(1.0, 0.0)


def test_dipole_per_monomer_units():
    # 0.16 C/m^2 over 108 A^3 and two monomers: 0.54 e.A each
    assert ES.dipole_per_monomer(0.16, 108.0, 2) == pytest.approx(0.539, rel=0.01)


def test_pvdf_alpha_to_beta_transition_strain_is_about_eleven_percent():
    ts = ES.transition_strain(PVDF, [T, GP, T, GM], [T, T], THREE_STATE)
    assert 2.2 < ts["rise_ground"] < 2.45
    assert 2.5 < ts["rise_polar"] < 2.65
    assert 0.07 < ts["axial"] < 0.16
    with_area = ES.transition_strain(PVDF, [T, GP, T, GM], [T, T], THREE_STATE, area_ground=23.9, area_polar=21.1)
    assert with_area["areal"] == pytest.approx(21.1 / 23.9 - 1.0)
    assert with_area["transverse_isotropic"] < 0.0 < with_area["axial"]
