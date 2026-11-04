"""Level 2 of the loss model: kinetic Monte Carlo of switching units on a lattice, driven.

Level 1 (:mod:`polyfind.loss`) is a small-signal relaxation: one dielectric strength, one
relaxation time, a linear response about the bias point.  The program's operating point
is a large standing field with a small AC drive, but the question "tan delta at each
performance level" needs the response at any amplitude, where the switching units no
longer respond linearly, the loss is a minor-loop area rather than ``pi eps0 eps'' E^2``,
and the cooperativity that the Vogel-Fulcher ``Tf`` only parametrises has to come from the
units' interactions.  This module simulates that directly.

**Model.**  Units of ``n`` monomers (the cooperative units of level 1) sit on a square
lattice; along ``x`` they are consecutive segments of one chain, along ``y`` neighbouring
chains.  Each unit is in its ground (``c = 0``) or polar (``c = 1``) conformation; the
polar state carries ``dmu_u`` more dipole along the drive axis and costs ``dE_i`` more
energy (the margin of the sequence that unit sits on; a spread of margins is the
relaxor's disorder).  With the conformational coupling in Ising form,
``H = -sum_<ij> J_ij (2c_i - 1)(2c_j - 1) + sum_i dE_i c_i - dmu_u E(t) sum_i p_i`` with
``p_i = c_i`` (two states), the energy difference between the two states of an isolated
unit is ``dE_i - dmu_u E``, exactly level 1's two-state unit, and ``J`` is what makes units
switch together.

**Three states.**  With ``three_state=True`` the polar conformation carries its dipole up
or down along the drive axis, ``p_i`` in ``{-1, 0, +1}``, the strain the same for either
(it is conformational), and a dipolar coupling ``-sum_<ij> J_d,ij p_i p_j`` between polar
neighbours.  At zero bias the two polar states are equally occupied, so the linear strain
vanishes (the strain is at twice the drive frequency) while the dielectric strength does
not; the linear strain grows with the bias as the down state freezes out.  Transitions
are ground to either polar state and back; a direct reversal (up to down) is allowed with
its own barrier ``flip_barrier_ev`` or, by default, goes through the ground state.

**Dynamics.**  Single transitions with Arrhenius rates ``nu0 exp(-max(Ea + dH/2, 0) / kT)``
(``dH`` the energy change of the transition; the barrier top midway, a
Bronsted-Evans-Polanyi slope of one half, the same rule as level 1's tilted well), which
satisfy detailed balance wherever the floor at zero is not reached.  Rejection-free
(Gillespie) kinetic Monte Carlo with the field ``E(t) = E_dc + E_ac sin(w t)`` evaluated at
each event; the events are nanoseconds apart against drive periods of microseconds, so
holding the field fixed between events is a safe approximation at every frequency this
module is used for.

**Outputs.**  The polarization ``P(t) = N dmu_u m`` (``m`` the mean of ``p``, the net polar
fraction) and the strain ``S(t) = s_switch x`` (``x`` the mean of ``c``, the polar fraction),
sampled over the measured periods; the fundamental Fourier components give an effective
``eps'``, ``eps''`` and ``tan delta`` at that amplitude, the loop area ``oint E dP`` gives the
loss per cycle directly, the strain's second harmonic is reported beside its fundamental,
and the same per-cycle bookkeeping as level 1 gives work, efficiency and heat.  At a
constant field the lattice's own fluctuations give the small-signal permittivity and the
strain per field (:meth:`IsingKMC.linear_response`: the variance of ``P`` and the
cross-correlation of ``x`` with ``P``, with their time dependence).

**Calibration.**  :meth:`IsingKMC.from_channel` takes a level-1 channel (its unit size,
unit dipole, density, margin, attempt time and relaxation time at the calibration
temperature) and a coupling ``J``, and sets the barrier so that a unit flipping against its
``z`` aligned neighbours relaxes on the channel's time scale: ``Ea = kT ln(nu0 tau) - z J``.
The lattice then supplies its own temperature and amplitude dependence, which is what
level 1 could not.  What the lattice still is not: the crystal-amorphous interface, where
the neutron work puts the dispersion; that is level 3.

Cost: events per drive period is about ``N / (f tau)``, each a few tens of microseconds
in NumPy, so a 24 x 24 lattice at 50 kHz is a few seconds per period.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .loss import E_A_TO_C_M, EPS0, K_B_EV, V_PER_A_TO_V_PER_M, SwitchingChannel, vogel_fulcher_tau


@dataclass
class KMCResult:
    """One driven run: the sampled trajectory and what the fundamental harmonic says."""

    freq: float
    field_dc: float
    field_ac: float
    temperature: float
    t: np.ndarray
    field: np.ndarray
    polar_fraction: np.ndarray
    net_polar_fraction: np.ndarray
    polarization: np.ndarray  # C/m^2, above the all-ground state
    strain: np.ndarray
    events: int
    periods: int
    eps_real: float
    eps_imag: float
    tan_delta: float
    loop_loss_j_per_m3: float  # oint E dP per cycle
    p_dc: float
    s_dc: float
    s_ac: float
    s_phase_deg: float
    s_2f: float  # the strain's second-harmonic amplitude (the electrostrictive response at zero bias)
    harmonic_ratio: float  # |third harmonic| / |fundamental| of P: the nonlinearity
    drift_per_cycle: float  # change of the polar fraction per cycle over the window, over its AC amplitude

    def operating_point(self, young_modulus_pa: float, tan_delta_mech: float = 0.0) -> dict:
        """Per-cycle bookkeeping with the loop loss as the electrical loss."""
        loss_e = self.loop_loss_j_per_m3
        loss_m = np.pi * young_modulus_pa * tan_delta_mech * self.s_ac ** 2
        work = 0.25 * young_modulus_pa * self.s_ac ** 2
        total = work + loss_e + loss_m
        return {
            "frequency_hz": self.freq, "tan_delta": self.tan_delta, "eps_real": self.eps_real,
            "eps_imag": self.eps_imag, "strain_ac": self.s_ac, "strain_phase_deg": self.s_phase_deg,
            "strain_2f": self.s_2f,
            "blocking_stress_pa": young_modulus_pa * self.s_ac, "work_j_per_m3": float(work),
            "loss_electrical_j_per_m3": float(loss_e), "loss_mechanical_j_per_m3": float(loss_m),
            "efficiency": float(work / total) if total > 0 else 0.0,
            "heat_w_per_m3": float(self.freq * (loss_e + loss_m)),
            "stroke_rate_per_lambda": self.s_ac * self.freq, "harmonic_ratio": self.harmonic_ratio,
        }


class IsingKMC:
    """Switching units on an ``L x L`` periodic square lattice (see the module docstring).

    ``margins_ev`` is the polar-state cost per unit, a scalar or an ``(L*L,)`` array;
    ``coupling_ev`` is ``J`` (scalar) or ``(J_along_chain, J_across_chains)``, the
    conformational coupling (like conformations on neighbouring units are favoured);
    ``barrier_ev`` the switching barrier ``Ea``; ``attempt_rate`` ``nu0`` in 1/s; ``n_units``
    the unit density per m^3 (for ``P`` in C/m^2); ``s_switch`` the strain of a switched
    unit's monomers along the drive axis.  ``three_state=True`` gives the polar state a
    dipole sign (up or down along the drive axis; the strain is the same for either), with
    ``dipolar_coupling_ev`` (``J_d``, scalar or a pair) between the dipoles of polar
    neighbours and ``flip_barrier_ev`` the barrier of a direct up-down reversal (``None``:
    reversal goes through the ground state).
    """

    def __init__(self, L: int, dmu_unit: float, margins_ev, coupling_ev, barrier_ev: float, attempt_rate: float,
                 temperature: float, n_units: float, s_switch: float = 0.11, seed: int | None = 0,
                 three_state: bool = False, dipolar_coupling_ev=0.0, flip_barrier_ev: float | None = None):
        self.L = int(L)
        self.N = self.L * self.L
        self.dmu = float(dmu_unit)
        self.margins = np.broadcast_to(np.asarray(margins_ev, dtype=float), (self.N,)).copy()
        jx, jy = (coupling_ev, coupling_ev) if np.isscalar(coupling_ev) else coupling_ev
        self.jx, self.jy = float(jx), float(jy)
        jdx, jdy = (dipolar_coupling_ev, dipolar_coupling_ev) if np.isscalar(dipolar_coupling_ev) else dipolar_coupling_ev
        self.jdx, self.jdy = float(jdx), float(jdy)
        self.ea = float(barrier_ev)
        self.ea_flip = None if flip_barrier_ev is None else float(flip_barrier_ev)
        self.nu0 = float(attempt_rate)
        self.T = float(temperature)
        self.kT = K_B_EV * self.T
        self.n_units = float(n_units)
        self.s_switch = float(s_switch)
        self.three_state = bool(three_state)
        self.K = 2 if self.three_state else 1  # candidate transitions per unit
        self.rng = np.random.default_rng(seed)
        idx = np.arange(self.N).reshape(self.L, self.L)
        self.right = np.roll(idx, -1, axis=1).ravel()
        self.left = np.roll(idx, 1, axis=1).ravel()
        self.up = np.roll(idx, -1, axis=0).ravel()
        self.down = np.roll(idx, 1, axis=0).ravel()
        self.c = np.zeros(self.N, dtype=float)  # conformation: 1 polar, 0 ground
        self.p = np.zeros(self.N, dtype=float)  # dipole sign of the polar state, 0 in the ground state
        self.time = 0.0
        self.tau_ref = 1.0 / (2.0 * self.nu0 * np.exp(-self.ea / self.kT))  # the two-state time; from_channel replaces it

    @property
    def sigma(self) -> np.ndarray:
        """The conformation as an Ising spin, ``2c - 1``."""
        return 2.0 * self.c - 1.0

    # --- construction from level 1 ------------------------------------------------------
    @classmethod
    def from_channel(cls, channel: SwitchingChannel, L: int = 24, coupling_ev=0.02, temperature: float = 300.0,
                     margin_spread_ev: float = 0.0, barrier_ev: float | None = None, seed: int | None = 0,
                     calibrate: bool = True, samples: int = 8000, three_state: bool | None = None,
                     dipolar_coupling_ev=0.0, flip_barrier_ev: float | None = None) -> "IsingKMC":
        """A lattice whose units are the channel's.

        ``coupling_ev`` is ``J`` (or ``(J_x, J_y)``); ``margin_spread_ev`` the standard
        deviation of the units' margins about the channel's (the random field that keeps an
        interacting lattice from ordering into one domain, which is what makes it a relaxor
        and not a ferroelectric).  ``three_state`` defaults to the channel's ``states == 3``.
        The barrier starts at ``kT ln(nu0 tau) - z J``, the value at which one unit flipping
        against ``z`` aligned neighbours relaxes on the channel's Vogel-Fulcher time, and with
        ``calibrate=True`` is then corrected by measurement: every rate carries the same
        ``exp(-Ea / kT)``, so the lattice's own zero-field autocorrelation time scales exactly
        with it and one measurement fixes ``Ea`` to the channel's ``tau`` at ``temperature``.
        The lattice then supplies the temperature and amplitude dependence on its own.
        """
        jx, jy = (coupling_ev, coupling_ev) if np.isscalar(coupling_ev) else coupling_ev
        zj = 2.0 * float(jx) + 2.0 * float(jy)
        nu0 = 1.0 / channel.tau0
        kT = K_B_EV * temperature
        tau = float(vogel_fulcher_tau(temperature, channel.tau0, channel.ea_ev, channel.tf))
        if barrier_ev is None:
            barrier_ev = kT * np.log(nu0 * tau) - zj
            if barrier_ev < 0.0:
                raise ValueError(f"coupling too strong for the channel's relaxation time: Ea would be {barrier_ev:.3f} eV")
        if three_state is None:
            three_state = getattr(channel, "states", 2) == 3
        rng = np.random.default_rng(seed)
        margins = channel.de_unit_ev + (margin_spread_ev * rng.standard_normal(L * L) if margin_spread_ev > 0 else 0.0)
        model = cls(L=L, dmu_unit=channel.dmu_unit, margins_ev=margins, coupling_ev=(jx, jy), barrier_ev=float(barrier_ev),
                    attempt_rate=nu0, temperature=temperature, n_units=channel.n_units, s_switch=channel.s_switch, seed=seed,
                    three_state=three_state, dipolar_coupling_ev=dipolar_coupling_ev, flip_barrier_ev=flip_barrier_ev)
        model.calibration = {"tau_target": tau, "barrier_formula": float(barrier_ev)}
        model.tau_ref = tau
        if calibrate:
            model.calibrate(tau, samples=samples)
        return model

    def at_temperature(self, temperature: float, seed: int | None = 1, tau_target: float | None = None,
                       samples: int = 8000) -> "IsingKMC":
        """The same lattice (units, margins, couplings, barriers, density) at another temperature.

        With ``tau_target`` (the bench fit's loss-peak time at that temperature) the clone's
        barrier is recalibrated to it: the lattice's own temperature dependence is the
        barrier plus a few hundredths of an eV of cooperativity whatever the coupling (the
        coupling x spread grid of plan section 12), so the Vogel-Fulcher slope is taken as
        an input this way, and the lattice supplies the bias and amplitude dependence at
        each temperature.
        """
        clone = IsingKMC(L=self.L, dmu_unit=self.dmu, margins_ev=self.margins, coupling_ev=(self.jx, self.jy),
                         barrier_ev=self.ea, attempt_rate=self.nu0, temperature=temperature, n_units=self.n_units,
                         s_switch=self.s_switch, seed=seed, three_state=self.three_state,
                         dipolar_coupling_ev=(self.jdx, self.jdy), flip_barrier_ev=self.ea_flip)
        if tau_target is not None:
            clone.calibration = {"tau_target": float(tau_target), "barrier_formula": self.ea}
            clone.calibrate(float(tau_target), samples=samples)
        return clone

    def calibrate(self, tau_target: float, samples: int = 8000, field_dc: float = 0.0, iterations: int = 2,
                  tolerance: float = 0.05, observable: str = "centroid") -> dict:
        """Shift ``Ea`` so that the lattice's relaxation time equals ``tau_target`` in the chosen observable.

        ``observable="centroid"`` (default): the loss peak sits at ``1 / (2 pi tau_target)``.
        ``observable="integral"``: the integral (mean) relaxation time of the polarization,
        ``int phi dt``, equals ``tau_target``; this is what a meter that sees only the
        low-frequency flank of the peak delivers, as ``eps'' / (delta_eps omega)`` at
        ``omega tau << 1`` (the LCR-6100's 10 Hz to 100 kHz at room temperature against a
        peak near 0.5 MHz), so the bench's flank time and the lattice's are the same quantity.

        The position of the small-signal loss peak (:meth:`loss_peak`, from the
        fluctuations at ``field_dc``; its centroid in ``log f``, which is the peak for a
        symmetric peak and far less noisy than the maximum of a broad one) is the quantity
        the Vogel-Fulcher fit of the bench sweep is made on, and it scales exactly with
        ``exp(-Ea / kT)`` because that factor multiplies every rate (up to the floor on
        negative barriers; the direct-reversal barrier, when there is one, is shifted by the
        same amount).  A fitted time constant does not serve here: the polarization
        autocorrelation of a coupled lattice is multi-exponential (units at a conformational
        interface reverse fast, units inside a polar region slowly), so its fit depends on
        the sampling step.  The grid is centred on the target, and the centroid of a broad
        peak that sits off-centre is pulled toward the centre by the grid's truncation, so a
        single shift under-corrects; the measurement and shift are repeated ``iterations``
        times (stopping early once within ``tolerance`` of the target), and the bias vanishes
        as the peak approaches the centre.  Returns and stores the centroid time before and
        after, the history, the peak time after, and the peak's width in decades.
        """
        history = []
        shift_total = 0.0
        at_edge = False
        first = float("nan")
        if observable not in ("centroid", "integral"):
            raise ValueError("observable must be 'centroid' or 'integral'")
        for _ in range(max(1, int(iterations))):
            if observable == "centroid":
                pk = self.loss_peak(field_dc, samples=samples, tau_center=tau_target)
                measured = pk["tau_centroid"]
                at_edge = at_edge or bool(pk["at_edge"])
            else:
                lr = self.linear_response(np.array([1.0]), field_dc=field_dc, samples=samples, dt=tau_target / 20.0,
                                          relaxation=100.0 * tau_target)
                measured = lr["tau_integral"]
            history.append(float(measured))
            if not np.isfinite(first):
                first = float(measured)
            if not (np.isfinite(measured) and measured > 0):
                break
            shift = self.kT * np.log(tau_target / measured)
            if abs(shift) < tolerance * self.kT:
                break
            self.ea = float(self.ea + shift)
            if self.ea_flip is not None:
                self.ea_flip = float(self.ea_flip + shift)
            shift_total += shift
        self.tau_ref = float(tau_target)
        after = self.loss_peak(field_dc, samples=samples, tau_center=tau_target)
        tau_after = float(after["tau_centroid"]) if observable == "centroid" else float(after["tau_integral"])
        info = {"tau_target": float(tau_target), "tau_before": first, "tau_history": history, "observable": observable,
                "barrier_shift_ev": float(shift_total), "barrier_ev": self.ea, "tau_after": tau_after,
                "tau_centroid_after": float(after["tau_centroid"]), "tau_integral_after": float(after["tau_integral"]),
                "tau_peak_after": float(after["tau_peak"]), "width_decades": float(after["width_decades"]),
                "peak_at_edge": bool(at_edge or after["at_edge"])}
        getattr(self, "calibration", {}).update(info) if hasattr(self, "calibration") else setattr(self, "calibration", info)
        return info

    # --- energetics ------------------------------------------------------------------------
    def local_field(self, field: float) -> np.ndarray:
        """Two states: ``h_i = (dmu_u E - dE_i) / 2``, the field on the Ising spin."""
        return 0.5 * (self.dmu * field - self.margins)

    def neighbour_sum(self) -> np.ndarray:
        """``sum_j J_ij (2 c_j - 1)``: the conformational field of the neighbours."""
        s = 2.0 * self.c - 1.0
        return self.jx * (s[self.right] + s[self.left]) + self.jy * (s[self.up] + s[self.down])

    def dipolar_sum(self) -> np.ndarray:
        """``sum_j J_d,ij p_j``: the dipolar field of the neighbours (three states)."""
        p = self.p
        return self.jdx * (p[self.right] + p[self.left]) + self.jdy * (p[self.up] + p[self.down])

    def transition_energies(self, field: float) -> np.ndarray:
        """``dH`` of each candidate transition at the given field, shape ``(N, K)``.

        Two states: the flip.  Three states: column 0 is ground to up or polar to ground,
        column 1 is ground to down or polar to the reversed polar state.
        """
        sc = self.neighbour_sum()
        if not self.three_state:
            return (2.0 * (2.0 * self.c - 1.0) * (sc + self.local_field(field)))[:, None]
        sp = self.dipolar_sum()
        polar = self.c > 0.5
        p = self.p
        dc0 = np.where(polar, -1.0, 1.0)
        dp0 = np.where(polar, -p, 1.0)
        dc1 = np.where(polar, 0.0, 1.0)
        dp1 = np.where(polar, -2.0 * p, -1.0)

        def dh(dc, dp):
            return -2.0 * dc * sc - dp * sp + self.margins * dc - self.dmu * field * dp

        return np.stack([dh(dc0, dp0), dh(dc1, dp1)], axis=1)

    def flip_energies(self, field: float) -> np.ndarray:
        """``dH_i`` for the first candidate transition of each unit (two states: the flip)."""
        return self.transition_energies(field)[:, 0]

    def rates(self, field: float) -> np.ndarray:
        """The rate of every candidate transition, flat ``(N * K,)``."""
        dh = self.transition_energies(field)
        ea = np.full(dh.shape, self.ea)
        if self.three_state:
            polar = self.c > 0.5
            ea[polar, 1] = np.inf if self.ea_flip is None else self.ea_flip
        barrier = np.maximum(ea + 0.5 * dh, 0.0)
        return (self.nu0 * np.exp(-barrier / self.kT)).ravel()

    @property
    def polar_fraction(self) -> float:
        return float(self.c.mean())

    @property
    def net_polar_fraction(self) -> float:
        return float(self.p.mean())

    def polarization(self) -> float:
        return self.n_units * self.dmu * E_A_TO_C_M * self.net_polar_fraction

    # --- the simulation --------------------------------------------------------------------
    def _apply(self, idx: int) -> None:
        i, k = divmod(int(idx), self.K)
        if not self.three_state:
            self.c[i] = 1.0 - self.c[i]
            self.p[i] = self.c[i]
            return
        if self.c[i] > 0.5:
            if k == 0:
                self.c[i] = 0.0
                self.p[i] = 0.0
            else:
                self.p[i] = -self.p[i]
        else:
            self.c[i] = 1.0
            self.p[i] = 1.0 if k == 0 else -1.0

    def _advance(self, field_of_t, t_end: float, sample_times=None, record=None) -> int:
        """Run events until ``t_end``; record the state at each sample time passed."""
        events = 0
        k = 0
        n_samples = 0 if sample_times is None else len(sample_times)
        last = self.N * self.K - 1
        while self.time < t_end:
            r = self.rates(field_of_t(self.time))
            R = float(r.sum())
            if R <= 0.0:
                self.time = t_end
                break
            dt = -np.log(self.rng.random()) / R
            t_new = self.time + dt
            while k < n_samples and sample_times[k] <= t_new:
                if sample_times[k] <= t_end:
                    record(sample_times[k])
                k += 1
            if t_new > t_end:
                self.time = t_end
                break
            self.time = t_new
            i = int(np.searchsorted(np.cumsum(r), self.rng.random() * R))
            self._apply(min(i, last))
            events += 1
        while k < n_samples and sample_times[k] <= t_end:
            record(sample_times[k])
            k += 1
        return events

    def equilibrate(self, field_dc: float, duration: float) -> int:
        """Hold the DC field for ``duration`` seconds (events only)."""
        t_end = self.time + duration
        return self._advance(lambda t: field_dc, t_end)

    def run(self, freq: float, field_dc: float, field_ac: float, periods: int = 6, equilibrate_periods: int = 2,
            samples_per_period: int = 64, dc_relaxation: float | None = None, detrend: bool = True) -> KMCResult:
        """Drive with ``E(t) = E_dc + E_ac sin(2 pi f t)`` (V/A) and measure the last ``periods``.

        ``dc_relaxation`` seconds at the DC field first (default: two hundred reference
        relaxation times, at most a millisecond; a relaxor settles slowly under a bias, and
        the result's ``drift_per_cycle`` says whether it had), then ``equilibrate_periods``
        of drive unrecorded, then the measured periods.  ``detrend`` removes a linear drift
        of the polar fractions over the window before the harmonic analysis, so that a slow
        residual settling does not masquerade as an in-phase response.
        """
        w = 2.0 * np.pi * freq
        period = 1.0 / freq
        if dc_relaxation is None:
            dc_relaxation = min(200.0 * self.tau_ref, 1e-3)
        self.equilibrate(field_dc, dc_relaxation)
        t0 = self.time
        field_of_t = lambda t: field_dc + field_ac * np.sin(w * (t - t0))  # noqa: E731
        self._advance(field_of_t, t0 + equilibrate_periods * period)
        start = t0 + equilibrate_periods * period
        samples = start + np.arange(periods * samples_per_period) * period / samples_per_period
        rec_t, rec_x, rec_m = [], [], []

        def record(ts):
            rec_t.append(ts)
            rec_x.append(self.polar_fraction)
            rec_m.append(self.net_polar_fraction)

        events = self._advance(field_of_t, start + periods * period, samples, record)
        t = np.asarray(rec_t)
        x = np.asarray(rec_x)
        m = np.asarray(rec_m)
        E = field_dc + field_ac * np.sin(w * (t - t0))
        return self._analyse(freq, field_dc, field_ac, t - t0, E, x, m, events, periods, detrend)

    def _analyse(self, freq, field_dc, field_ac, t, E, x, m, events, periods, detrend=True) -> KMCResult:
        w = 2.0 * np.pi * freq
        n = len(t)
        n_per = max(1, n // periods)
        drift_x = float(np.mean(x[-n_per:]) - np.mean(x[:n_per])) / max(periods - 1, 1)
        x_used, m_used = x.copy(), m.copy()
        if detrend and n > 2:
            for arr in (x_used, m_used):
                slope, _ = np.polyfit(t, arr, 1)
                arr -= slope * (t - t.mean())
        p_unit = self.n_units * self.dmu * E_A_TO_C_M
        P = p_unit * m_used
        S = self.s_switch * x_used
        sn, cs = np.sin(w * t), np.cos(w * t)
        p_s, p_c = 2.0 * np.mean(P * sn), 2.0 * np.mean(P * cs)
        s_s, s_c = 2.0 * np.mean(S * sn), 2.0 * np.mean(S * cs)
        s_2f = float(np.hypot(2.0 * np.mean(S * np.sin(2 * w * t)), 2.0 * np.mean(S * np.cos(2 * w * t))))
        p3 = np.hypot(2.0 * np.mean(P * np.sin(3 * w * t)), 2.0 * np.mean(P * np.cos(3 * w * t)))
        E_ac_si = field_ac * V_PER_A_TO_V_PER_M
        eps_real = p_s / (EPS0 * E_ac_si) if field_ac > 0 else float("nan")
        eps_imag = -p_c / (EPS0 * E_ac_si) if field_ac > 0 else float("nan")
        # loop area oint E_ac dP over the measured periods, per cycle (trapezoid in P): the AC
        # field only, so that a residual settling of P under the DC bias does not enter as E_dc dP
        E_si = (E - field_dc) * V_PER_A_TO_V_PER_M
        loop = float(np.sum(0.5 * (E_si[1:] + E_si[:-1]) * np.diff(P))) / periods if n > 1 else 0.0
        s_amp = float(np.hypot(s_s, s_c))
        s_phase = float(np.degrees(np.arctan2(s_c, s_s)))
        x_amp = s_amp / self.s_switch if self.s_switch else float("nan")
        return KMCResult(freq=float(freq), field_dc=float(field_dc), field_ac=float(field_ac), temperature=self.T,
                         t=t, field=E, polar_fraction=x, net_polar_fraction=m, polarization=p_unit * m,
                         strain=self.s_switch * x, events=int(events), periods=int(periods),
                         eps_real=float(eps_real), eps_imag=float(eps_imag),
                         tan_delta=float(eps_imag / eps_real) if field_ac > 0 and eps_real != 0 else float("nan"),
                         loop_loss_j_per_m3=loop, p_dc=float(P.mean()), s_dc=float(S.mean()), s_ac=s_amp,
                         s_phase_deg=s_phase, s_2f=s_2f,
                         harmonic_ratio=float(p3 / np.hypot(p_s, p_c)) if field_ac > 0 else float("nan"),
                         drift_per_cycle=float(drift_x / x_amp) if x_amp > 0 else float("nan"))

    def calibrate_strength(self, delta_eps_target: float, field_dc: float = 0.0, samples: int = 12000) -> dict:
        """Rescale the unit density so that the lattice's static susceptibility at ``field_dc``
        equals the measured dielectric strength.

        The dynamics do not depend on ``n_units`` (it only converts a polar fraction to a
        polarization), so this calibration is independent of :meth:`calibrate`; it absorbs
        what the coupling does to the strength (a Curie-Weiss enhancement at zero bias) into
        the switchable fraction, which is the quantity level 1 calibrated the same way.
        """
        lr = self.linear_response(np.array([1.0]), field_dc=field_dc, samples=samples)
        chi = float(lr["chi_static"])
        factor = delta_eps_target / chi if chi > 0 else 1.0
        self.n_units *= factor
        info = {"chi_before": chi, "density_factor": float(factor), "n_units": float(self.n_units)}
        if hasattr(self, "calibration"):
            self.calibration.update(info)
        return info

    def trajectory(self, field_dc: float, samples: int, dt: float, relaxation: float) -> tuple[np.ndarray, np.ndarray]:
        """The polar fraction and the net polar fraction sampled every ``dt`` at a constant
        field after ``relaxation`` seconds."""
        self.equilibrate(field_dc, relaxation)
        ts = self.time + np.arange(int(samples)) * dt
        xs: list[float] = []
        ms: list[float] = []

        def record(t_s):
            xs.append(self.polar_fraction)
            ms.append(self.net_polar_fraction)

        self._advance(lambda t: field_dc, ts[-1] + dt, ts, record)
        return np.asarray(xs), np.asarray(ms)

    @staticmethod
    def _response_kernel(corr: np.ndarray, dt: float, w: np.ndarray) -> tuple[np.ndarray, float, np.ndarray]:
        """``1 - i w int phi(t) exp(-i w t) dt`` for the normalised correlation ``phi = corr / corr[0]``.

        The measured ``phi`` is used where it is significant (positive and above twice its
        standard error by Bartlett's estimate, ``se_k^2 = (1 + 2 sum_{j<k} phi_j^2) / n``) and
        continued beyond that point by the exponential tail with the time constant of its
        initial decay (a weighted log-linear fit over ``phi > 0.3``).  The tail starts below a
        few percent of ``phi(0)``, so what it assumes (a Debye tail) bounds the error at a few
        percent of the integral time; integrating a noisy tail to its first zero crossing was
        found to inflate the time by up to half on a few hundred units.  Returns the kernel at
        each ``w``, the integral time (the same construction) and the ``phi`` used.
        """
        phi = corr / corr[0]
        n = len(phi)
        se = np.sqrt((1.0 + 2.0 * np.concatenate(([0.0, 0.0], np.cumsum(phi[1:] ** 2)[:-1]))) / n)  # sum_{j<k}, k >= 2
        weak = (phi <= 0.0) | (phi < 2.0 * se)
        weak[0] = False
        cut = int(np.argmax(weak)) if np.any(weak) else n
        cut = max(cut, 2)
        upto = int(np.argmax(phi < 0.3)) if np.any(phi < 0.3) else n
        upto = max(min(upto, cut), 3)
        tt = np.arange(upto) * dt
        wts = phi[:upto]
        denom = np.sum(wts * tt * tt)
        slope = np.sum(wts * tt * np.log(np.maximum(phi[:upto], 1e-300))) / denom if denom > 0 else 0.0
        tau_fit = -1.0 / slope if slope < 0 else 0.0
        phi = phi[:cut]
        t = np.arange(cut) * dt
        t_c, phi_c = t[-1], max(float(phi[-1]), 0.0)
        tau_int = float(np.trapezoid(phi, t) + phi_c * tau_fit)
        kernel = np.exp(-1j * np.outer(w, t))  # (n_freq, n_t)
        integral = np.trapezoid(phi[None, :] * kernel, t, axis=1)
        integral = integral + phi_c * np.exp(-1j * w * t_c) * tau_fit / (1.0 + 1j * w * tau_fit)
        return 1.0 - 1j * w * integral, tau_int, phi

    def linear_response(self, freqs, field_dc: float = 0.0, samples: int = 20000, dt: float | None = None,
                        relaxation: float | None = None, eps_inf: float = 0.0) -> dict:
        """The small-signal permittivity and strain per field at a bias from the equilibrium
        fluctuations (fluctuation-dissipation).

        No drive is applied: the lattice sits at ``field_dc`` and its polar fractions are
        sampled every ``dt`` (default ``tau_ref / 20``) for ``samples`` steps after
        ``relaxation`` (default ``100 tau_ref``).  The static susceptibility is
        ``V <dP^2> / (eps0 kT)`` with ``V`` the lattice's volume, and the frequency
        dependence is ``chi(w) / chi(0) = 1 - i w int_0^inf phi(t) exp(-i w t) dt`` with ``phi``
        the normalised autocorrelation of ``P`` (where it is significant, with a fitted
        exponential tail beyond; see :meth:`_response_kernel`).  The strain per field is the same construction on the cross-correlation of the polar
        fraction with the polarization (Kubo): ``dx/dE = dmu_u N <dx dm> / kT``, with the
        normalised cross-correlation for its time dependence; for a three-state lattice at
        zero bias it is zero within noise.  This uses every fluctuation of the trajectory,
        which a small drive cannot, so it is the right estimator for the small-signal column;
        :meth:`run` is for the large-signal one.

        Returns ``freqs``, ``eps_real``, ``eps_imag``, ``tan_delta`` (with ``eps_inf`` added
        to the real part), ``chi_static``, ``tau_integral`` (the integral of ``phi``, the
        Debye-equivalent time), ``x_mean``, ``m_mean``, ``dx_de`` (static, per V/m), ``d_eff``
        (the complex strain per V/m at each frequency), and the sampled ``phi`` and
        ``phi_strain`` with their ``dt``.

        Measured against the exact two-state reference: the strength to 2%, the time to
        12%, tan delta to 10% at ``omega tau`` of 0.3 and 1.  Above the relaxation the real
        part is the small difference ``1 - w int phi sin(wt)`` and a few percent of noise in
        ``phi`` becomes tens of percent in ``eps'``, so the estimator is for frequencies at
        and below the relaxation, which is where the operating points are.
        """
        freqs = np.atleast_1d(np.asarray(freqs, dtype=float))
        if dt is None:
            dt = self.tau_ref / 20.0
        if relaxation is None:
            relaxation = 100.0 * self.tau_ref
        x, m = self.trajectory(field_dc, samples, dt, relaxation)
        x_mean, m_mean = float(x.mean()), float(m.mean())
        dx, dm = x - x_mean, m - m_mean
        volume = self.N / self.n_units  # m^3
        p_unit = self.n_units * self.dmu * E_A_TO_C_M  # C/m^2 per unit of net polar fraction
        var_p = float(np.var(dm)) * p_unit ** 2
        kT_j = self.kT * 1.602176634e-19
        chi_static = volume * var_p / (EPS0 * kT_j)
        n = len(dm)
        norm = np.arange(n, 0, -1)
        c = np.correlate(dm, dm, mode="full")[n - 1:] / norm
        cxm = np.correlate(dx, dm, mode="full")[n - 1:] / norm  # <dx(t) dm(0)>
        dx_de = float(self.dmu * self.N * cxm[0] / self.kT / V_PER_A_TO_V_PER_M)  # per V/m, static
        w = 2.0 * np.pi * freqs
        if c[0] <= 0.0:
            nanv = np.full(len(freqs), np.nan)
            return {"freqs": freqs, "eps_real": nanv, "eps_imag": nanv, "tan_delta": nanv, "chi_static": 0.0,
                    "tau_integral": float("nan"), "x_mean": x_mean, "m_mean": m_mean, "dx_de": dx_de,
                    "d_eff": np.zeros(len(freqs), dtype=complex), "phi": c, "phi_strain": cxm, "dt": dt}
        kernel, tau_int, phi = self._response_kernel(c, dt, w)
        chi_w = chi_static * kernel
        eps_real = eps_inf + chi_w.real
        eps_imag = -chi_w.imag
        if cxm[0] > 0.0:
            kernel_s, _, phi_s = self._response_kernel(cxm, dt, w)
        else:  # no resolvable strain response (three states at zero bias): static value, which is noise about zero
            kernel_s, phi_s = np.ones(len(freqs), dtype=complex), cxm
        d_eff = self.s_switch * dx_de * kernel_s
        return {"freqs": freqs, "eps_real": eps_real, "eps_imag": eps_imag, "tan_delta": eps_imag / eps_real,
                "chi_static": float(chi_static), "tau_integral": tau_int, "x_mean": x_mean, "m_mean": m_mean,
                "dx_de": dx_de, "d_eff": d_eff, "phi": phi, "phi_strain": phi_s, "dt": dt}

    def loss_peak(self, field_dc: float = 0.0, samples: int = 3000, dt: float | None = None, relaxation: float | None = None,
                  decades: float = 3.0, points: int = 61, tau_center: float | None = None) -> dict:
        """The loss peak of the small-signal spectrum at a bias, from the lattice's fluctuations.

        ``eps''(f)`` by :meth:`linear_response` on ``points`` logarithmically spaced
        frequencies over ``decades`` centred on ``1 / (2 pi tau_center)`` (default
        ``tau_ref``; the sampling step defaults to a twentieth of it), its maximum located by
        a parabola in ``log f`` through the three points about it.  Returns ``f_peak``,
        ``tau_peak = 1 / (2 pi f_peak)`` (the quantity a Vogel-Fulcher fit of a sweep is made
        on), ``f_centroid`` and ``tau_centroid`` (the first moment of ``eps''`` in ``log f``
        over the grid: the same as the peak for a symmetric peak on a centred grid, and far
        less noisy for a broad one, whose maximum wanders over its flat top between
        realisations), ``eps_peak``, ``width_decades`` (the full width at half maximum in
        decades of frequency; 1.14 for a Debye relaxation, more for a broad one; ``nan`` when
        a half maximum lies outside the grid), ``at_edge`` (the maximum on the grid's end, so
        the centre was off), the static ``chi_static``, ``dx_de``, ``x_mean``, ``m_mean`` and
        the spectrum itself.
        """
        tau_c = self.tau_ref if tau_center is None else float(tau_center)
        if dt is None:
            dt = tau_c / 20.0
        if relaxation is None:
            relaxation = 100.0 * tau_c
        f0 = 1.0 / (2.0 * np.pi * tau_c)
        freqs = f0 * np.logspace(-0.5 * decades, 0.5 * decades, int(points))
        lr = self.linear_response(freqs, field_dc=field_dc, samples=samples, dt=dt, relaxation=relaxation)
        ei = np.asarray(lr["eps_imag"], dtype=float)
        out = {"freqs": freqs, "eps_imag": ei, "eps_real": lr["eps_real"], "chi_static": lr["chi_static"],
               "dx_de": lr["dx_de"], "x_mean": lr["x_mean"], "m_mean": lr["m_mean"], "tau_integral": lr["tau_integral"]}
        if not np.all(np.isfinite(ei)) or np.max(ei) <= 0.0:
            out.update({"f_peak": float("nan"), "tau_peak": float("nan"), "eps_peak": float("nan"),
                        "width_decades": float("nan"), "at_edge": True})
            return out
        lf = np.log10(freqs)
        step = lf[1] - lf[0]
        k = int(np.argmax(ei))
        at_edge = k in (0, len(ei) - 1)
        if not at_edge:
            y0, y1, y2 = ei[k - 1], ei[k], ei[k + 1]
            denom = y0 - 2.0 * y1 + y2
            frac = 0.5 * (y0 - y2) / denom if denom < 0 else 0.0
            lf_peak = lf[k] + frac * step
            e_peak = y1 - 0.25 * (y0 - y2) * frac
        else:
            lf_peak, e_peak = lf[k], ei[k]
        half = 0.5 * e_peak
        above = np.flatnonzero(ei >= half)
        width = float("nan")
        if len(above) and above[0] > 0 and above[-1] < len(ei) - 1 and not at_edge:
            lo, hi = above[0], above[-1]
            lf_lo = lf[lo - 1] + (half - ei[lo - 1]) / (ei[lo] - ei[lo - 1]) * step
            lf_hi = lf[hi] + (ei[hi] - half) / (ei[hi] - ei[hi + 1]) * step
            width = float(lf_hi - lf_lo)
        f_peak = float(10.0 ** lf_peak)
        pos = np.clip(ei, 0.0, None)
        lf_cen = float(np.sum(pos * lf) / np.sum(pos))
        f_cen = float(10.0 ** lf_cen)
        out.update({"f_peak": f_peak, "tau_peak": 1.0 / (2.0 * np.pi * f_peak), "eps_peak": float(e_peak),
                    "f_centroid": f_cen, "tau_centroid": 1.0 / (2.0 * np.pi * f_cen),
                    "width_decades": width, "at_edge": at_edge})
        return out

    def autocorrelation_time(self, field_dc: float = 0.0, samples: int = 4000, dt: float | None = None,
                             relaxation: float | None = None) -> float:
        """Integral relaxation time of the polarization at a constant field, in seconds.

        Samples ``P`` every ``dt`` (default: a twentieth of the single-flip time
        ``1 / (2 nu0 exp(-Ea/kT))``, right for independent units; an interacting lattice
        relaxes more slowly and wants ``dt`` set from its expected time, as
        :meth:`calibrate` does) after ``relaxation`` seconds, and fits the initial decay of
        the normalised autocorrelation.  A measurement of the lattice's own collective time,
        to be held against the channel's Vogel-Fulcher value.
        """
        r1 = self.nu0 * np.exp(-self.ea / self.kT)
        if dt is None:
            dt = 0.05 / (2.0 * r1)
        if relaxation is None:
            relaxation = 50.0 / r1
        _, m = self.trajectory(field_dc, samples, dt, relaxation)
        x = m - np.mean(m)
        if np.allclose(x, 0.0):
            return float("nan")
        n = len(x)
        c = np.correlate(x, x, mode="full")[n - 1:] / (np.arange(n, 0, -1))
        c = c / c[0]
        # the time constant of the initial decay: a weighted log-linear fit over the part of the
        # autocorrelation above 0.3, which is far less noisy than the integral to the first zero
        # crossing and, for a stretched decay, is its short-time scale
        upto = int(np.argmax(c < 0.3)) if np.any(c < 0.3) else n
        upto = max(upto, 3)
        tt = np.arange(upto) * dt
        wts = c[:upto]
        slope = np.sum(wts * tt * np.log(c[:upto])) / np.sum(wts * tt * tt)
        return float(-1.0 / slope) if slope < 0 else float("nan")


def lattice_observables(model: IsingKMC, temperatures=(270.0, 330.0), bias: float = 0.006, samples: int = 8000,
                        fdt_samples: int = 8000, ea_guess: float = 0.5, seed: int = 1) -> dict:
    """The bench's calibration numbers, measured on the lattice.

    ``apparent_ea_ev`` is the slope of ``ln tau`` against ``1 / kT`` between the two
    temperatures (the lattice's own loss-peak centroid times at each, :meth:`IsingKMC.loss_peak`
    centred on an Arrhenius guess at ``ea_guess`` from the reference time and re-centred
    when the peak lands on the grid's edge or more than 3x from the guess); ``suppression``
    the static strength at ``bias`` (V/A) over its zero-bias value (fluctuation-dissipation);
    ``width_decades`` the loss peak's full width at half maximum at the reference
    temperature (1.14 for Debye); ``chi_zero`` the zero-bias strength before any density
    scaling.  These are the quantities the bench's bias-suppression ratio, Vogel-Fulcher
    slope and loss-peak width calibrate ``n``, ``J`` and the margin spread against (the
    dielectric strength then sets the density, which the dynamics do not depend on).
    """
    t_lo, t_hi = sorted(float(x) for x in temperatures)
    taus = {}
    for tk in (t_lo, t_hi):
        mt = model.at_temperature(tk, seed=seed)
        guess = float(model.tau_ref * np.exp(ea_guess * (1.0 / mt.kT - 1.0 / model.kT)))
        tau = float("nan")
        for _ in range(4):
            pk = mt.loss_peak(0.0, samples=samples, tau_center=guess)
            tau = pk["tau_centroid"]
            if np.isfinite(tau) and tau > 0 and not pk["at_edge"] and 1.0 / 3.0 < tau / guess < 3.0:
                break
            guess = tau if np.isfinite(tau) and tau > 0 else guess * 3.0
        taus[tk] = tau
    k_lo, k_hi = K_B_EV * t_lo, K_B_EV * t_hi
    ea_app = float(np.log(taus[t_lo] / taus[t_hi]) / (1.0 / k_lo - 1.0 / k_hi))
    pk0 = model.loss_peak(0.0, samples=fdt_samples)
    lrb = model.linear_response(np.array([1.0]), field_dc=bias, samples=fdt_samples)
    chi0 = float(pk0["chi_static"])
    return {"tau_lo": taus[t_lo], "tau_hi": taus[t_hi], "temperatures": (t_lo, t_hi), "apparent_ea_ev": ea_app,
            "tau_peak": float(pk0["tau_peak"]), "tau_centroid": float(pk0["tau_centroid"]),
            "width_decades": float(pk0["width_decades"]),
            "chi_zero": chi0, "chi_bias": float(lrb["chi_static"]),
            "suppression": float(lrb["chi_static"] / chi0) if chi0 > 0 else float("nan"),
            "x_zero": float(pk0["x_mean"]), "x_bias": float(lrb["x_mean"]), "m_bias": float(lrb["m_mean"]),
            "dx_de_bias": float(lrb["dx_de"])}
