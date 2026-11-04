"""Level 2 of the loss model: the driven kinetic Monte Carlo against what it must reproduce."""
import numpy as np
import pytest

from polyfind import kmc as KM
from polyfind import loss as L

T = 300.0
KT = L.K_B_EV * T


def _independent(L_side=8, dmu=4.0, margin=0.0, ea=0.15, nu0=1e12, seed=1):
    """Non-interacting units: an exact two-state reference exists for everything."""
    return KM.IsingKMC(L=L_side, dmu_unit=dmu, margins_ev=margin, coupling_ev=0.0, barrier_ev=ea, attempt_rate=nu0,
                       temperature=T, n_units=1e27, s_switch=0.1, seed=seed)


def test_independent_units_reach_the_two_state_occupancy():
    m = _independent(L_side=12, margin=0.02)
    field = 0.004  # V/A: dmu E = 0.016 eV against a 0.02 eV margin
    x_exact = 1.0 / (1.0 + np.exp((0.02 - 4.0 * field) / KT))
    r1 = m.nu0 * np.exp(-m.ea / m.kT)
    m.equilibrate(field, 200.0 / r1)
    xs = []
    for _ in range(40):
        m.equilibrate(field, 20.0 / r1)
        xs.append(m.polar_fraction)
    assert np.mean(xs) == pytest.approx(x_exact, abs=0.03)


def test_independent_units_give_debye_loss_with_tan_delta_omega_tau():
    """Two states, no interactions, zero margin: tau = 1 / (2 r1), eps_inf = 0, so the
    fundamental's tan delta is exactly omega tau and the loop loss is pi eps0 eps'' E^2."""
    m = _independent(L_side=16, ea=0.20)
    r1 = m.nu0 * np.exp(-m.ea / m.kT)
    tau = 1.0 / (2.0 * r1)
    f = 0.3 / (2.0 * np.pi * tau)  # omega tau = 0.3
    # 0.0015 V/A keeps dmu E / 2kT at 0.12 (linear to half a percent) and the response above the
    # thermal noise of 256 units; 24 periods average it down
    res = m.run(freq=f, field_dc=0.0, field_ac=0.0015, periods=24, equilibrate_periods=4, samples_per_period=48)
    assert res.eps_real > 0 and res.eps_imag > 0
    assert res.tan_delta == pytest.approx(0.3, rel=0.3)
    assert res.s_phase_deg == pytest.approx(-16.7, abs=8.0)
    E_ac = 0.0015 * 1e10
    assert res.loop_loss_j_per_m3 == pytest.approx(np.pi * L.EPS0 * res.eps_imag * E_ac ** 2, rel=0.1)
    # the strain is the polar fraction's fundamental times the switch strain
    assert res.s_ac == pytest.approx(0.1 * res.eps_real * L.EPS0 * E_ac / np.cos(np.radians(res.s_phase_deg)) / (m.n_units * m.dmu * L.E_A_TO_C_M), rel=0.05)
    # and the level-1 channel built from the same numbers agrees on the strength
    ch = L.SwitchingChannel(dmu_mon=4.0, de_mon=0.0, n=1, fraction=1.0, rho_mon=1e27, tau0=1 / (2 * m.nu0), ea_ev=0.20,
                            eps_inf=0.0)
    assert res.eps_real == pytest.approx(ch.delta_eps(0.0, T) / (1 + 0.3 ** 2), rel=0.25)


def test_the_autocorrelation_time_of_independent_units_is_the_two_state_time():
    m = _independent(L_side=10, ea=0.18)
    r1 = m.nu0 * np.exp(-m.ea / m.kT)
    tau = m.autocorrelation_time(0.0, samples=6000)
    assert tau == pytest.approx(1.0 / (2.0 * r1), rel=0.25)


def test_coupling_slows_the_relaxation_and_the_calibration_accounts_for_it():
    ch = L.SwitchingChannel(dmu_mon=0.45, de_mon=0.0, n=17, fraction=0.5, rho_mon=1.5e28, tau0=8.08e-13, ea_ev=0.124,
                            tf=191.0)
    tau_vf = ch.tau(T)
    free = KM.IsingKMC.from_channel(ch, L=10, coupling_ev=0.0, temperature=T, seed=3, calibrate=False)
    coupled = KM.IsingKMC.from_channel(ch, L=10, coupling_ev=0.03, temperature=T, seed=3, calibrate=False)
    assert free.ea > coupled.ea  # the barrier makes room for the neighbours' hold
    assert coupled.ea == pytest.approx(free.ea - 4 * 0.03)
    # the single-flip time of the free lattice is the channel's time (that is the formula)
    assert 1.0 / (2.0 * free.nu0 * np.exp(-free.ea / free.kT)) == pytest.approx(tau_vf / 2.0, rel=1e-6)
    with pytest.raises(ValueError):
        KM.IsingKMC.from_channel(ch, L=6, coupling_ev=0.2, temperature=T, calibrate=False)
    # the measured calibration: independent units relax on 1/(2 r1) = tau/2, so Ea moves up by kT ln 2
    cal = KM.IsingKMC.from_channel(ch, L=10, coupling_ev=0.0, temperature=T, seed=3, calibrate=True, samples=2000)
    assert cal.calibration["barrier_shift_ev"] == pytest.approx(cal.kT * np.log(2.0), abs=0.6 * cal.kT)
    assert cal.calibration["tau_after"] == pytest.approx(tau_vf, rel=0.5)


def test_large_amplitude_loss_is_not_the_small_signal_one():
    """Driving a degenerate interacting lattice through most of its switching range makes
    the loop loss grow faster than E^2 and the third harmonic appear."""
    ch = L.SwitchingChannel(dmu_mon=0.45, de_mon=0.0, n=17, fraction=0.5, rho_mon=1.5e28, tau0=8.08e-13, ea_ev=0.124,
                            tf=191.0)
    m = KM.IsingKMC.from_channel(ch, L=12, coupling_ev=0.01, temperature=T, seed=5, calibrate=False)
    f = 0.2 * m.nu0 * np.exp(-m.ea / m.kT) / (2 * np.pi)
    small = m.run(f, 0.0, 0.0005, periods=10, equilibrate_periods=2, samples_per_period=32)
    big = m.run(f, 0.0, 0.01, periods=10, equilibrate_periods=2, samples_per_period=32)
    # dmu_u E = 0.077 eV = 3 kT at the large drive: the units saturate, the strain amplitude grows
    # far less than the field, the loss per E^2 falls, and the third harmonic is real (at the small
    # drive it is thermal noise on 144 units, so it is not compared)
    assert big.s_ac < 0.5 * (0.01 / 0.0005) * small.s_ac
    assert big.loop_loss_j_per_m3 < 0.5 * small.loop_loss_j_per_m3 * (0.01 / 0.0005) ** 2
    assert big.harmonic_ratio > 0.05
    assert abs(big.drift_per_cycle) < 0.2
    op = big.operating_point(2e8, 0.05)
    assert 0.0 < op["efficiency"] < 1.0 and op["heat_w_per_m3"] == pytest.approx(
        f * (op["loss_electrical_j_per_m3"] + op["loss_mechanical_j_per_m3"]))


def test_fluctuation_dissipation_reproduces_debye_for_independent_units():
    """At a constant field the lattice's own fluctuations give the small-signal permittivity:
    static strength N dmu^2 x(1-x) / (eps0 kT) and tan delta = omega tau, tau = 1/(2 r1)."""
    m = _independent(L_side=16, ea=0.20, margin=0.01)
    r1 = m.nu0 * np.exp(-m.ea / m.kT)
    tau = 1.0 / (2.0 * r1)
    field = 0.002
    x = 1.0 / (1.0 + np.exp((0.01 - 4.0 * field) / KT))
    lam = 2.0 * r1 * np.cosh((4.0 * field - 0.01) / (2.0 * KT))  # the tilted two-state rate
    tau_e = 1.0 / lam
    freqs = np.array([0.1, 0.3, 1.0]) / (2.0 * np.pi * tau_e)  # at and below the relaxation (see the docstring)
    lr = m.linear_response(freqs, field_dc=field, samples=20000, dt=tau_e / 20.0, relaxation=100.0 * tau_e)
    ch = L.SwitchingChannel(dmu_mon=4.0, de_mon=0.01 * 23.060548, n=1, fraction=1.0,
                            rho_mon=1e27, tau0=1.0 / (2.0 * m.nu0), ea_ev=0.20, eps_inf=0.0)
    assert lr["x_mean"] == pytest.approx(x, abs=0.03)
    assert lr["chi_static"] == pytest.approx(ch.delta_eps(field, T), rel=0.2)
    assert lr["tau_integral"] == pytest.approx(tau_e, rel=0.25)
    for wt, td in zip((0.1, 0.3, 1.0), lr["tan_delta"]):
        assert td == pytest.approx(wt, rel=0.3)


def _three(L_side=10, dmu=4.0, margin=0.0, ea=0.15, nu0=1e12, seed=1, **kw):
    """Non-interacting three-state units: ground, polar up, polar down."""
    return KM.IsingKMC(L=L_side, dmu_unit=dmu, margins_ev=margin, coupling_ev=0.0, barrier_ev=ea, attempt_rate=nu0,
                       temperature=T, n_units=1e27, s_switch=0.1, seed=seed, three_state=True, **kw)


def test_three_state_units_reach_the_three_state_occupancies():
    m = _three(L_side=12, margin=0.02)
    field = 0.004
    a = np.exp(-(0.02 - 4.0 * field) / KT)
    b = np.exp(-(0.02 + 4.0 * field) / KT)
    z = 1.0 + a + b
    r = m.nu0 * np.exp(-m.ea / m.kT)
    m.equilibrate(field, 200.0 / r)
    xs, ms = [], []
    for _ in range(40):
        m.equilibrate(field, 20.0 / r)
        xs.append(m.polar_fraction)
        ms.append(m.net_polar_fraction)
    assert np.mean(xs) == pytest.approx((a + b) / z, abs=0.03)
    assert np.mean(ms) == pytest.approx((a - b) / z, abs=0.03)


def test_three_state_fluctuations_give_the_strength_and_a_strain_that_vanishes_at_zero_bias():
    m = _three(L_side=16, ea=0.20, margin=0.01)
    three = L.SwitchingChannel(dmu_mon=4.0, de_mon=0.01 * 23.060548, n=1, fraction=1.0, rho_mon=1e27,
                               tau0=1.0 / m.nu0, ea_ev=0.20, eps_inf=0.0, states=3)
    two = L.SwitchingChannel(dmu_mon=4.0, de_mon=0.01 * 23.060548, n=1, fraction=1.0, rho_mon=1e27,
                             tau0=1.0 / m.nu0, ea_ev=0.20, eps_inf=0.0)
    for field in (0.0, 0.003):
        tau_e = three.tau(T, field)  # the slow mode of the independent three-state unit
        lr = m.linear_response(np.array([0.2]) / (2.0 * np.pi * tau_e), field_dc=field, samples=20000, dt=tau_e / 20.0,
                               relaxation=100.0 * tau_e)
        x, mm = three.occupancies(field, T)
        assert lr["x_mean"] == pytest.approx(x, abs=0.03)
        assert lr["m_mean"] == pytest.approx(mm, abs=0.03)
        assert lr["chi_static"] == pytest.approx(three.delta_eps(field, T), rel=0.2)
        if field == 0.0:
            assert lr["tau_integral"] == pytest.approx(tau_e, rel=0.3)  # at zero bias the polarization is the slow mode
            assert abs(lr["dx_de"]) < 0.15 * two.dx_de(0.0, T)  # no linear strain: noise about zero
        else:
            assert lr["dx_de"] == pytest.approx(three.dx_de(field, T), rel=0.3)
            assert abs(lr["d_eff"][0]) == pytest.approx(0.1 * three.dx_de(field, T), rel=0.35)


def test_three_state_strain_is_at_twice_the_drive_frequency_at_zero_bias():
    m = _three(L_side=16, ea=0.20)
    r = m.nu0 * np.exp(-m.ea / m.kT)
    f = 0.1 * r / (2.0 * np.pi)
    res = m.run(f, 0.0, 0.01, periods=12, equilibrate_periods=2, samples_per_period=48)  # dmu E = 1.5 kT
    assert res.s_2f > 3.0 * res.s_ac
    assert res.eps_real > 0.0  # the polarization still responds at the drive frequency
    assert abs(res.p_dc) < 0.2 * res.eps_real * L.EPS0 * 0.01 * 1e10


def test_a_direct_reversal_path_speeds_the_polarization_relaxation():
    slow = _three(L_side=10, ea=0.18, seed=2)
    fast = _three(L_side=10, ea=0.18, seed=2, flip_barrier_ev=0.18)
    r = slow.nu0 * np.exp(-slow.ea / slow.kT)
    t_slow = slow.autocorrelation_time(0.0, samples=6000, dt=0.05 / r, relaxation=50.0 / r)
    t_fast = fast.autocorrelation_time(0.0, samples=6000, dt=0.05 / r, relaxation=50.0 / r)
    assert t_slow == pytest.approx(1.0 / r, rel=0.25)  # the polarization mode of the degenerate unit: 1/r, not 1/(2r)
    assert t_fast < 0.7 * t_slow  # with the reversal at the same barrier the mode relaxes three times faster


def test_from_channel_follows_the_channels_states_and_clones_at_temperature():
    ch = L.literature_terpolymer(states=3)
    m = KM.IsingKMC.from_channel(ch, L=8, coupling_ev=0.01, temperature=T, seed=3, calibrate=False)
    assert m.three_state and m.K == 2 and m.dmu == pytest.approx(ch.dmu_unit)
    hot = m.at_temperature(330.0)
    assert hot.three_state and hot.T == 330.0 and hot.ea == m.ea and np.array_equal(hot.margins, m.margins)
    two = KM.IsingKMC.from_channel(L.literature_terpolymer(), L=8, coupling_ev=0.01, temperature=T, calibrate=False)
    assert not two.three_state and two.K == 1


def test_the_loss_peak_of_independent_units_is_the_debye_peak():
    m = _independent(L_side=16, ea=0.20)
    r1 = m.nu0 * np.exp(-m.ea / m.kT)
    tau = 1.0 / (2.0 * r1)
    pk = m.loss_peak(0.0, samples=12000, tau_center=tau)
    assert not pk["at_edge"]
    assert pk["tau_peak"] == pytest.approx(tau, rel=0.25)
    assert pk["tau_centroid"] == pytest.approx(tau, rel=0.2)  # symmetric in log f: the centroid is the peak
    assert pk["width_decades"] == pytest.approx(1.14, abs=0.3)  # log10((2 + sqrt 3) / (2 - sqrt 3))
    assert pk["eps_peak"] == pytest.approx(0.5 * pk["chi_static"], rel=0.3)  # a Debye peak is half the strength
    # a centre 1.5 decades above the peak with a 2-decade grid lands the maximum on the grid's edge and says so
    off = m.loss_peak(0.0, samples=3000, tau_center=tau / 30.0, decades=2.0)
    assert off["at_edge"]
