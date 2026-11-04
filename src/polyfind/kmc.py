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
chains.  Each unit is in its ground (``sigma = -1``) or polar (``sigma = +1``) state; the
polar state carries ``dmu_u`` more dipole along the drive axis and costs ``dE_i`` more
energy (the margin of the sequence that unit sits on; a spread of margins is the
relaxor's disorder).  With the Ising form ``H = -sum_<ij> J_ij sigma_i sigma_j - sum_i h_i(t)
sigma_i`` and ``h_i(t) = (dmu_u E(t) - dE_i) / 2``, the energy difference between the two
states of an isolated unit is ``dE_i - dmu_u E``, exactly level 1's two-state unit, and
``J`` is what makes units switch together.

**Dynamics.**  Single flips with Arrhenius rates ``nu0 exp(-max(Ea + dH/2, 0) / kT)``
(``dH`` the energy change of the flip; the barrier top midway, a Bronsted-Evans-Polanyi
slope of one half, the same rule as level 1's tilted well), which satisfy detailed balance
wherever the floor at zero is not reached.  Rejection-free (Gillespie) kinetic Monte Carlo
with the field ``E(t) = E_dc + E_ac sin(w t)`` evaluated at each event; the events are
nanoseconds apart against drive periods of microseconds, so holding the field fixed
between events is a safe approximation at every frequency this module is used for.

**Outputs.**  The polarization ``P(t) = N dmu_u <s>`` (``s = (sigma + 1) / 2`` the polar
fraction) and the strain ``S(t) = s_switch <s>``, sampled over the measured periods; the
fundamental Fourier components give an effective ``eps'``, ``eps''`` and ``tan delta`` at
that amplitude, the loop area ``oint E dP`` gives the loss per cycle directly, and the
same per-cycle bookkeeping as level 1 gives work, efficiency and heat.

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
            "blocking_stress_pa": young_modulus_pa * self.s_ac, "work_j_per_m3": float(work),
            "loss_electrical_j_per_m3": float(loss_e), "loss_mechanical_j_per_m3": float(loss_m),
            "efficiency": float(work / total) if total > 0 else 0.0,
            "heat_w_per_m3": float(self.freq * (loss_e + loss_m)),
            "stroke_rate_per_lambda": self.s_ac * self.freq, "harmonic_ratio": self.harmonic_ratio,
        }


class IsingKMC:
    """Switching units on an ``L x L`` periodic square lattice (see the module docstring).

    ``margins_ev`` is the polar-state cost per unit, a scalar or an ``(L*L,)`` array;
    ``coupling_ev`` is ``J`` (scalar) or ``(J_along_chain, J_across_chains)``;
    ``barrier_ev`` the flip barrier ``Ea``; ``attempt_rate`` ``nu0`` in 1/s; ``n_units`` the
    unit density per m^3 (for ``P`` in C/m^2); ``s_switch`` the strain of a switched unit's
    monomers along the drive axis.
    """

    def __init__(self, L: int, dmu_unit: float, margins_ev, coupling_ev, barrier_ev: float, attempt_rate: float,
                 temperature: float, n_units: float, s_switch: float = 0.11, seed: int | None = 0):
        self.L = int(L)
        self.N = self.L * self.L
        self.dmu = float(dmu_unit)
        self.margins = np.broadcast_to(np.asarray(margins_ev, dtype=float), (self.N,)).copy()
        jx, jy = (coupling_ev, coupling_ev) if np.isscalar(coupling_ev) else coupling_ev
        self.jx, self.jy = float(jx), float(jy)
        self.ea = float(barrier_ev)
        self.nu0 = float(attempt_rate)
        self.T = float(temperature)
        self.kT = K_B_EV * self.T
        self.n_units = float(n_units)
        self.s_switch = float(s_switch)
        self.rng = np.random.default_rng(seed)
        idx = np.arange(self.N).reshape(self.L, self.L)
        self.right = np.roll(idx, -1, axis=1).ravel()
        self.left = np.roll(idx, 1, axis=1).ravel()
        self.up = np.roll(idx, -1, axis=0).ravel()
        self.down = np.roll(idx, 1, axis=0).ravel()
        self.sigma = -np.ones(self.N, dtype=float)
        self.time = 0.0
        self.tau_ref = 1.0 / (2.0 * self.nu0 * np.exp(-self.ea / self.kT))  # the two-state time; from_channel replaces it

    # --- construction from level 1 ------------------------------------------------------
    @classmethod
    def from_channel(cls, channel: SwitchingChannel, L: int = 24, coupling_ev=0.02, temperature: float = 300.0,
                     margin_spread_ev: float = 0.0, barrier_ev: float | None = None, seed: int | None = 0,
                     calibrate: bool = True, samples: int = 3000) -> "IsingKMC":
        """A lattice whose units are the channel's.

        ``coupling_ev`` is ``J`` (or ``(J_x, J_y)``); ``margin_spread_ev`` the standard
        deviation of the units' margins about the channel's (the random field that keeps an
        interacting lattice from ordering into one domain, which is what makes it a relaxor
        and not a ferroelectric).  The barrier starts at ``kT ln(nu0 tau) - z J``, the value
        at which one unit flipping against ``z`` aligned neighbours relaxes on the channel's
        Vogel-Fulcher time, and with ``calibrate=True`` is then corrected by measurement: every
        rate carries the same ``exp(-Ea / kT)``, so the lattice's own zero-field autocorrelation
        time scales exactly with it and one measurement fixes ``Ea`` to the channel's ``tau``
        at ``temperature``.  The lattice then supplies the temperature and amplitude
        dependence on its own.
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
        rng = np.random.default_rng(seed)
        margins = channel.de_unit_ev + (margin_spread_ev * rng.standard_normal(L * L) if margin_spread_ev > 0 else 0.0)
        model = cls(L=L, dmu_unit=channel.dmu_unit, margins_ev=margins, coupling_ev=(jx, jy), barrier_ev=float(barrier_ev),
                    attempt_rate=nu0, temperature=temperature, n_units=channel.n_units, s_switch=channel.s_switch, seed=seed)
        model.calibration = {"tau_target": tau, "barrier_formula": float(barrier_ev)}
        model.tau_ref = tau
        if calibrate:
            model.calibrate(tau, samples=samples)
        return model

    def calibrate(self, tau_target: float, samples: int = 3000, field_dc: float = 0.0) -> dict:
        """Shift ``Ea`` so that the measured zero-field autocorrelation time equals ``tau_target``.

        Exact up to the floor on negative barriers, because ``exp(-Ea / kT)`` multiplies
        every rate.  Returns and stores the measurement before and after.
        """
        measured = self.autocorrelation_time(field_dc, samples=samples, dt=tau_target / 20.0, relaxation=20.0 * tau_target)
        shift = self.kT * np.log(tau_target / measured) if np.isfinite(measured) and measured > 0 else 0.0
        self.ea = float(self.ea + shift)
        self.tau_ref = float(tau_target)
        check = self.autocorrelation_time(field_dc, samples=samples, dt=tau_target / 20.0, relaxation=20.0 * tau_target)
        info = {"tau_target": float(tau_target), "tau_before": float(measured), "barrier_shift_ev": float(shift),
                "barrier_ev": self.ea, "tau_after": float(check)}
        getattr(self, "calibration", {}).update(info) if hasattr(self, "calibration") else setattr(self, "calibration", info)
        return info

    # --- energetics ------------------------------------------------------------------------
    def local_field(self, field: float) -> np.ndarray:
        return 0.5 * (self.dmu * field - self.margins)

    def neighbour_sum(self) -> np.ndarray:
        s = self.sigma
        return self.jx * (s[self.right] + s[self.left]) + self.jy * (s[self.up] + s[self.down])

    def flip_energies(self, field: float) -> np.ndarray:
        """``dH_i`` for flipping each unit at the given field."""
        return 2.0 * self.sigma * (self.neighbour_sum() + self.local_field(field))

    def rates(self, field: float) -> np.ndarray:
        barrier = np.maximum(self.ea + 0.5 * self.flip_energies(field), 0.0)
        return self.nu0 * np.exp(-barrier / self.kT)

    @property
    def polar_fraction(self) -> float:
        return float(0.5 * (self.sigma.mean() + 1.0))

    def polarization(self) -> float:
        return self.n_units * self.dmu * E_A_TO_C_M * self.polar_fraction

    # --- the simulation --------------------------------------------------------------------
    def _advance(self, field_of_t, t_end: float, sample_times=None, record=None) -> int:
        """Run events until ``t_end``; record the state at each sample time passed."""
        events = 0
        k = 0
        n_samples = 0 if sample_times is None else len(sample_times)
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
            self.sigma[min(i, self.N - 1)] *= -1.0
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
        of the polar fraction over the window before the harmonic analysis, so that a slow
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
        rec_t, rec_x = [], []

        def record(ts):
            rec_t.append(ts)
            rec_x.append(self.polar_fraction)

        events = self._advance(field_of_t, start + periods * period, samples, record)
        t = np.asarray(rec_t)
        x = np.asarray(rec_x)
        E = field_dc + field_ac * np.sin(w * (t - t0))
        return self._analyse(freq, field_dc, field_ac, t - t0, E, x, events, periods, detrend)

    def _analyse(self, freq, field_dc, field_ac, t, E, x, events, periods, detrend=True) -> KMCResult:
        w = 2.0 * np.pi * freq
        m = len(t)
        n_per = max(1, m // periods)
        drift_x = float(np.mean(x[-n_per:]) - np.mean(x[:n_per])) / max(periods - 1, 1)
        x_used = x.copy()
        if detrend and m > 2:
            slope, intercept = np.polyfit(t, x, 1)
            x_used = x - slope * (t - t.mean())
        P = self.n_units * self.dmu * E_A_TO_C_M * x_used
        S = self.s_switch * x_used
        sn, cs = np.sin(w * t), np.cos(w * t)
        p_s, p_c = 2.0 * np.mean(P * sn), 2.0 * np.mean(P * cs)
        s_s, s_c = 2.0 * np.mean(S * sn), 2.0 * np.mean(S * cs)
        p3 = np.hypot(2.0 * np.mean(P * np.sin(3 * w * t)), 2.0 * np.mean(P * np.cos(3 * w * t)))
        E_ac_si = field_ac * V_PER_A_TO_V_PER_M
        eps_real = p_s / (EPS0 * E_ac_si) if field_ac > 0 else float("nan")
        eps_imag = -p_c / (EPS0 * E_ac_si) if field_ac > 0 else float("nan")
        # loop area oint E_ac dP over the measured periods, per cycle (trapezoid in P): the AC
        # field only, so that a residual settling of P under the DC bias does not enter as E_dc dP
        E_si = (E - field_dc) * V_PER_A_TO_V_PER_M
        loop = float(np.sum(0.5 * (E_si[1:] + E_si[:-1]) * np.diff(P))) / periods if m > 1 else 0.0
        s_amp = float(np.hypot(s_s, s_c))
        s_phase = float(np.degrees(np.arctan2(s_c, s_s)))
        x_amp = s_amp / self.s_switch if self.s_switch else float("nan")
        return KMCResult(freq=float(freq), field_dc=float(field_dc), field_ac=float(field_ac), temperature=self.T,
                         t=t, field=E, polar_fraction=x, polarization=self.n_units * self.dmu * E_A_TO_C_M * x,
                         strain=self.s_switch * x, events=int(events), periods=int(periods),
                         eps_real=float(eps_real), eps_imag=float(eps_imag),
                         tan_delta=float(eps_imag / eps_real) if field_ac > 0 and eps_real != 0 else float("nan"),
                         loop_loss_j_per_m3=loop, p_dc=float(P.mean()), s_dc=float(S.mean()), s_ac=s_amp,
                         s_phase_deg=s_phase, harmonic_ratio=float(p3 / np.hypot(p_s, p_c)) if field_ac > 0 else float("nan"),
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

    def trajectory(self, field_dc: float, samples: int, dt: float, relaxation: float) -> np.ndarray:
        """The polar fraction sampled every ``dt`` at a constant field after ``relaxation`` seconds."""
        self.equilibrate(field_dc, relaxation)
        ts = self.time + np.arange(int(samples)) * dt
        xs: list[float] = []
        self._advance(lambda t: field_dc, ts[-1] + dt, ts, lambda t_s: xs.append(self.polar_fraction))
        return np.asarray(xs)

    def linear_response(self, freqs, field_dc: float = 0.0, samples: int = 20000, dt: float | None = None,
                        relaxation: float | None = None, eps_inf: float = 0.0) -> dict:
        """The small-signal permittivity at a bias from the equilibrium fluctuations (fluctuation-dissipation).

        No drive is applied: the lattice sits at ``field_dc`` and its polar fraction is
        sampled every ``dt`` (default ``tau_ref / 20``) for ``samples`` steps after
        ``relaxation`` (default ``100 tau_ref``).  The static susceptibility is
        ``V <dP^2> / (eps0 kT)`` with ``V`` the lattice's volume, and the frequency
        dependence is ``chi(w) / chi(0) = 1 - i w int_0^inf phi(t) exp(-i w t) dt`` with ``phi``
        the normalised autocorrelation of ``P``, integrated to its first zero crossing.  This
        uses every fluctuation of the trajectory, which a small drive cannot, so it is the
        right estimator for the small-signal column; :meth:`run` is for the large-signal one.

        Returns ``freqs``, ``eps_real``, ``eps_imag``, ``tan_delta`` (with ``eps_inf`` added
        to the real part), ``chi_static``, ``tau_integral`` (the integral of ``phi``, the
        Debye-equivalent time), ``x_mean`` and the sampled ``phi`` with its ``dt``.

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
        x = self.trajectory(field_dc, samples, dt, relaxation)
        x_mean = float(x.mean())
        dx = x - x_mean
        volume = self.N / self.n_units  # m^3
        p_unit = self.n_units * self.dmu * E_A_TO_C_M  # C/m^2 per unit of polar fraction
        var_p = float(np.var(dx)) * p_unit ** 2
        kT_j = self.kT * 1.602176634e-19
        chi_static = volume * var_p / (EPS0 * kT_j)
        n = len(dx)
        c = np.correlate(dx, dx, mode="full")[n - 1:] / np.arange(n, 0, -1)
        if c[0] <= 0.0:
            nanv = np.full(len(freqs), np.nan)
            return {"freqs": freqs, "eps_real": nanv, "eps_imag": nanv, "tan_delta": nanv, "chi_static": 0.0,
                    "tau_integral": float("nan"), "x_mean": x_mean, "phi": c, "dt": dt}
        phi = c / c[0]
        zero = int(np.argmax(phi <= 0.0)) if np.any(phi <= 0.0) else n
        phi = phi[:max(zero, 2)]
        t = np.arange(len(phi)) * dt
        tau_int = float(np.trapezoid(phi, t))
        w = 2.0 * np.pi * freqs
        kernel = np.exp(-1j * np.outer(w, t))  # (n_freq, n_t)
        integral = np.trapezoid(phi[None, :] * kernel, t, axis=1)
        chi_w = chi_static * (1.0 - 1j * w * integral)
        eps_real = eps_inf + chi_w.real
        eps_imag = -chi_w.imag
        return {"freqs": freqs, "eps_real": eps_real, "eps_imag": eps_imag, "tan_delta": eps_imag / eps_real,
                "chi_static": float(chi_static), "tau_integral": tau_int, "x_mean": x_mean, "phi": phi, "dt": dt}

    def autocorrelation_time(self, field_dc: float = 0.0, samples: int = 4000, dt: float | None = None,
                             relaxation: float | None = None) -> float:
        """Integral relaxation time of the polarization at a constant field, in seconds.

        Samples ``P`` every ``dt`` (default: a twentieth of the single-flip time
        ``1 / (2 nu0 exp(-Ea/kT))``, right for independent units; an interacting lattice
        relaxes more slowly and wants ``dt`` set from its expected time, as
        :meth:`calibrate` does) after ``relaxation`` seconds, and integrates the normalised
        autocorrelation to its first zero crossing.  A measurement of the lattice's own
        collective time, to be held against the channel's Vogel-Fulcher value.
        """
        r1 = self.nu0 * np.exp(-self.ea / self.kT)
        if dt is None:
            dt = 0.05 / (2.0 * r1)
        if relaxation is None:
            relaxation = 50.0 / r1
        self.equilibrate(field_dc, relaxation)
        ts = self.time + np.arange(samples) * dt
        xs = []
        self._advance(lambda t: field_dc, ts[-1] + dt, ts, lambda t_s: xs.append(self.polar_fraction))
        x = np.asarray(xs) - np.mean(xs)
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
