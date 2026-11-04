"""Level 1 of the loss model: the derived small-signal tan delta."""
import numpy as np
import pytest

from polyfind import loss as L


def test_vogel_fulcher_reproduces_the_terpolymers_room_temperature_peak():
    tau = float(L.vogel_fulcher_tau(300.0, 8.08e-13, 0.124, 191.0))
    assert 3e-7 < tau < 6e-7
    f = 1.0 / (2 * np.pi * tau)
    assert 2.5e5 < f < 5e5  # about 0.4 MHz
    assert float(L.vogel_fulcher_tau(300.0, 1e-13, 0.3, 0.0)) == pytest.approx(1e-13 * np.exp(0.3 / (L.K_B_EV * 300.0)))


def test_the_vogel_fulcher_fit_recovers_its_parameters():
    T = np.array([250.0, 265.0, 280.0, 300.0, 320.0, 340.0])
    f = 1.0 / (2 * np.pi * L.vogel_fulcher_tau(T, 8.08e-13, 0.124, 191.0))
    fit = L.fit_vogel_fulcher(T, f)
    assert fit["tf"] == pytest.approx(191.0, abs=1.0)
    assert fit["ea_ev"] == pytest.approx(0.124, rel=0.03)
    assert fit["tau0"] == pytest.approx(8.08e-13, rel=0.2)
    assert fit["residual_ln_f"] < 1e-3


def test_bias_suppression_inverts_to_the_unit_dipole():
    dmu = L.unit_dipole_from_bias_suppression(0.5, 0.006, 300.0)
    ch = L.SwitchingChannel(dmu_mon=dmu, de_mon=0.0, n=1, fraction=1.0, rho_mon=1e27, tau0=1e-12, ea_ev=0.1)
    assert ch.delta_eps(0.006, 300.0) / ch.delta_eps(0.0, 300.0) == pytest.approx(0.5, rel=1e-9)
    with pytest.raises(ValueError):
        L.unit_dipole_from_bias_suppression(1.0, 0.006)


def test_degenerate_channel_has_the_textbook_strength_and_debye_limits():
    ch = L.SwitchingChannel(dmu_mon=1.0, de_mon=0.0, n=4, fraction=1.0, rho_mon=1e27, tau0=1e-12, ea_ev=0.2,
                            eps_inf=3.0)
    mu = 4.0 * L.E_A_TO_C_M
    kT = L.K_B_EV * 300.0 * 1.602176634e-19
    assert ch.delta_eps(0.0, 300.0) == pytest.approx((1e27 / 4) * mu * mu / (4 * L.EPS0 * kT), rel=1e-9)
    assert ch.polar_fraction(0.0, 300.0) == pytest.approx(0.5)
    f_lo = ch.f_relax(300.0) * 1e-3
    f_hi = ch.f_relax(300.0) * 1e3
    de = ch.delta_eps()
    # far below the relaxation: eps'' ~ de * w tau, eps' ~ eps_inf + de
    e = ch.permittivity(f_lo, 300.0)
    assert e.real == pytest.approx(3.0 + de, rel=1e-5)
    assert -e.imag == pytest.approx(de * 2 * np.pi * f_lo * ch.tau(300.0), rel=1e-4)
    # at the peak, tan delta of the channel alone is de / (2 eps_inf + de) and eps'' = de / 2
    e = ch.permittivity(ch.f_relax(300.0), 300.0)
    assert -e.imag == pytest.approx(de / 2, rel=1e-9)
    assert float(ch.tan_delta(ch.f_relax(300.0), 300.0)) == pytest.approx(de / (2 * 3.0 + de), rel=1e-9)
    # far above: nothing relaxes
    e = ch.permittivity(f_hi, 300.0)
    assert e.real == pytest.approx(3.0, rel=1e-4)


def test_bias_both_suppresses_the_strength_and_shortens_the_relaxation():
    ch = L.literature_terpolymer()
    assert ch.de_mon == 0.0 and ch.n >= 10 and 0.2 < ch.fraction < 0.8
    assert ch.delta_eps(0.006) == pytest.approx(0.5 * ch.delta_eps(0.0), rel=0.05)  # n is rounded to a whole monomer count
    assert ch.tau(300.0, 0.006) < ch.tau(300.0, 0.0)
    assert ch.f_relax(300.0, 0.0) == pytest.approx(1 / (2 * np.pi * float(L.vogel_fulcher_tau(300.0, 8.08e-13, 0.124, 191.0))))
    td0 = float(ch.tan_delta(5e4, 300.0, 0.0))
    tdb = float(ch.tan_delta(5e4, 300.0, 0.0055))
    assert 0.02 < td0 < 0.3 and tdb < td0


def test_cole_cole_broadening_lowers_the_peak_and_spreads_the_loss():
    deb = L.literature_terpolymer(alpha=0.0)
    cc = L.literature_terpolymer(alpha=0.4)
    fr = deb.f_relax(300.0)
    assert float(cc.tan_delta(fr)) < float(deb.tan_delta(fr))
    assert float(cc.tan_delta(fr * 1e-3)) > float(deb.tan_delta(fr * 1e-3))


def test_a_margin_above_the_ground_state_kills_the_zero_bias_channel():
    base = L.literature_terpolymer()
    up = base.scaled(de_mon=0.3, label="0.3 kcal/mol up")  # 0.3 kcal/mol per monomer, n monomers per unit
    assert up.delta_eps(0.0) < 1e-3 * base.delta_eps(0.0)
    # and the threshold field brings it back: at E* the polar fraction is one half
    e_star = up.de_unit_ev / up.dmu_unit
    assert up.polar_fraction(e_star) == pytest.approx(0.5)
    assert up.delta_eps(e_star) == pytest.approx(base.delta_eps(0.0), rel=1e-6)


def test_operating_point_bookkeeping():
    ch = L.literature_terpolymer()
    op = ch.operating_point(freq=5e4, field_bias=0.0055, field_ac=0.00018, young_modulus_pa=2e8,
                            tan_delta_mech=0.05)
    assert op["tan_delta"] == pytest.approx(op["eps_imag"] / op["eps_real"])
    E_ac = 0.00018 * 1e10
    assert op["loss_electrical_j_per_m3"] == pytest.approx(np.pi * L.EPS0 * op["eps_imag"] * E_ac ** 2)
    assert op["strain_ac"] == pytest.approx(op["d_eff_m_per_v"] * E_ac)
    assert op["work_j_per_m3"] == pytest.approx(0.25 * 2e8 * op["strain_ac"] ** 2)
    assert op["loss_mechanical_j_per_m3"] == pytest.approx(np.pi * 2e8 * 0.05 * op["strain_ac"] ** 2)
    total = op["work_j_per_m3"] + op["loss_electrical_j_per_m3"] + op["loss_mechanical_j_per_m3"]
    assert op["efficiency"] == pytest.approx(op["work_j_per_m3"] / total)
    assert op["heat_w_per_m3"] == pytest.approx(5e4 * (op["loss_electrical_j_per_m3"] + op["loss_mechanical_j_per_m3"]))
    assert 0.0 < op["efficiency"] < 1.0
    assert -90.0 < op["strain_phase_deg"] <= 0.0


def test_three_state_unit_has_the_full_strength_but_no_linear_strain_at_zero_bias():
    import dataclasses

    two = L.SwitchingChannel(dmu_mon=0.45, de_mon=0.0, n=17, fraction=0.5, rho_mon=1.5e28, tau0=8.08e-13, ea_ev=0.124,
                             tf=191.0)
    three = dataclasses.replace(two, states=3)
    assert three.polar_fraction(0.0) == pytest.approx(2.0 / 3.0)
    assert three.net_polar_fraction(0.0) == 0.0
    # both polar states respond: the strength is x_polar / (x (1 - x)) times the two-state one
    assert three.delta_eps(0.0) == pytest.approx(two.delta_eps(0.0) * (2.0 / 3.0) / 0.25)
    assert three.dx_de(0.0) == 0.0
    assert three.tau(300.0, 0.0) == pytest.approx(two.tau(300.0, 0.0))
    kT = L.K_B_EV * 300.0
    for E in (0.001, 0.003, 0.0055):
        x, m = three.occupancies(E)
        assert 0.0 < m < x < 1.0
        # the strain per field is a two-state unit's at the same polar fraction times tanh(dmu_u E / kT)
        ref = three.dmu_unit * x * (1.0 - x) / kT / L.V_PER_A_TO_V_PER_M
        assert three.dx_de(E) / ref == pytest.approx(np.tanh(three.dmu_unit * E / kT))
        assert three.delta_eps(E) / three.delta_eps(0.0) == pytest.approx(L.bias_suppression(three.dmu_unit * E / kT, 3))
    assert three.tau(300.0, 0.0055) < three.tau(300.0, 0.0)  # the tilt speeds the slow mode too
    # a margin: the polar states are rarer, the strain per field smaller, the slow mode faster
    up = dataclasses.replace(three, de_mon=0.1)
    assert up.polar_fraction(0.0) < three.polar_fraction(0.0)
    assert up.dx_de(0.003) < three.dx_de(0.003)
    assert up.tau(300.0, 0.0) < three.tau(300.0, 0.0)


def test_three_state_calibration_needs_a_smaller_unit_for_the_same_suppression():
    d2 = L.unit_dipole_from_bias_suppression(0.5, 0.006, 300.0)
    d3 = L.unit_dipole_from_bias_suppression(0.5, 0.006, 300.0, states=3)
    assert d3 < d2
    kT = L.K_B_EV * 300.0
    assert L.bias_suppression(d3 * 0.006 / kT, 3) == pytest.approx(0.5, rel=1e-6)
    assert L.bias_suppression(d2 * 0.006 / kT, 2) == pytest.approx(0.5, rel=1e-6)
    assert L.bias_suppression(0.0, 3) == 1.0
    ch3 = L.literature_terpolymer(states=3)
    ch2 = L.literature_terpolymer()
    assert ch3.states == 3 and ch3.n < ch2.n
    assert ch3.delta_eps(0.006) / ch3.delta_eps(0.0) == pytest.approx(0.5, abs=0.05)  # n is an integer
    assert ch3.delta_eps(0.0) == pytest.approx(45.0)
    op = ch3.operating_point(50e3, 0.0055, 0.00018, 2e8)
    assert 0.0 < op["net_polar_fraction"] < op["polar_fraction"] < 1.0
