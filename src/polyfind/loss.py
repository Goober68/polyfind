"""Level 1 of the loss model: a derived small-signal tan delta for the switching channel.

The relaxor terpolymers' loss, strain and frequency limit are one relaxation
(``docs/ELECTROSTRICTION_RESEARCH.md`` sections 2 and 3): units of ``n`` monomers switch
between a ground and a polar conformation, their switchable polarization relaxes with a
Vogel-Fulcher time, and a DC bias both suppresses the switchable fraction and speeds the
relaxation.  This module writes that down in closed form so that tan delta, the strain per
AC field and the cycle efficiency can be evaluated at any frequency, temperature and bias
from inputs that are computed (``dmu`` per monomer, the margin ``dE``, the attempt time from
the libration frequencies, the barrier from a NEB) or calibrated on the bench sweep
(``tau0, Ea, Tf``; the bias suppression; the dielectric strength).  Everything is a
single-channel estimate: a Debye relaxation, broadened to Cole-Cole when a distribution of
barriers is wanted, plus a frequency-independent background.  Level 2 (kinetic Monte Carlo
on the sequence) replaces the small-signal assumption; level 3 (molecular dynamics of the
interface) replaces the single channel.

**Switching channel.**  ``N`` units per m^3, each with dipole gain ``dmu_u`` (e.A) and
energy cost ``dE_u`` (eV) for the polar state, in nominal field ``E`` (V/A): the polar
fraction is ``x = 1 / (1 + exp((dE_u - dmu_u E) / kT))`` and the differential
susceptibility ``chi = N dmu_u^2 x (1 - x) / (eps0 kT)``.  A degenerate unit (``dE_u = 0``)
gives the familiar ``chi(0) = N dmu_u^2 / (4 eps0 kT)`` and the bias suppression
``sech^2(dmu_u E / 2kT)``; a unit far above its ground state contributes nothing at zero
field, which is why the margin decides the strain at working fields.

**Three states.**  With ``states=3`` the polar conformation carries its dipole up or down
along the drive axis (ground, polar up, polar down): ``x_polar = (a + b) / Z`` and the net
polar fraction ``m = (a - b) / Z`` with ``a, b = exp(-(dE_u -/+ dmu_u E) / kT)``,
``Z = 1 + a + b``; the susceptibility is ``N dmu_u^2 (x_polar - m^2) / (eps0 kT)`` and the
strain per field ``S_switch (dmu_u / kT) m (1 - x_polar)``, which vanishes at zero bias (the
strain is then at twice the drive frequency, the ``Q P^2`` of the literature) and rises with
the bias as the down state freezes out; the two-state forms are the same expressions with
``m = x``.  The bias suppression of a degenerate three-state unit is
``3 (cosh y + 2) / (1 + 2 cosh y)^2`` with ``y = dmu_u E / kT`` (two states:
``sech^2(y / 2)``), so the same measured suppression calibrates a smaller unit.

**Relaxation.**  ``tau(T) = tau0 exp(Ea / k (T - Tf))`` (``Tf = 0`` is Arrhenius); under
bias the double well tilts and ``tau(T, E) = tau(T) / cosh((dmu_u E - dE_u) / 2kT)`` with
the barrier top taken midway (a Bronsted-Evans-Polanyi slope of one half).  Then
``eps*(w) = eps_inf + chi / (1 + (i w tau)^(1 - alpha))`` (Cole-Cole; ``alpha = 0`` is
Debye) and the strain of the channel follows the same spectral function:
``S*(w) = S_switch (dx/dE) E_ac / (1 + (i w tau)^(1 - alpha))``.

**Operating point.**  Per cycle and per m^3: electrical loss ``pi eps0 eps'' E_ac^2`` (plus
a background ``pi eps0 eps' tan_bg E_ac^2``), mechanical loss ``pi Y tan_m S_ac^2``, work at
a matched load ``Y S_ac^2 / 4``; efficiency is work over work plus losses, and the heat is
the losses times the frequency, to be held under the program's budget.
"""
from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np

K_B_EV = 8.617333262e-5  # eV/K
EPS0 = 8.8541878128e-12  # F/m
E_A_TO_C_M = 1.602176634e-29  # C.m per e.A
V_PER_A_TO_V_PER_M = 1e10
AVOGADRO = 6.02214076e23


# --- relaxation time ---------------------------------------------------------------------
def vogel_fulcher_tau(temperature, tau0: float, ea_ev: float, tf: float = 0.0):
    """``tau0 exp(Ea / k (T - Tf))`` in seconds; ``Tf = 0`` is the Arrhenius form."""
    T = np.asarray(temperature, dtype=float)
    return tau0 * np.exp(ea_ev / (K_B_EV * (T - tf)))


def fit_vogel_fulcher(t_max, freq, tf_grid=None) -> dict:
    """Fit ``f = f0 exp(-Ea / k (T_max - Tf))`` to the temperature of the loss peak at each
    frequency, as the bench sweep delivers it.

    Linear in ``ln f0`` and ``Ea`` at each trial ``Tf``, the trial with the least residual
    kept.  Returns ``tau0 = 1 / (2 pi f0)``, ``ea_ev``, ``tf`` and the rms residual in
    ``ln f``.  ``tf_grid`` defaults to 0 .. min(T_max) - 5 K in 0.5 K steps.
    """
    T = np.asarray(t_max, dtype=float)
    lf = np.log(np.asarray(freq, dtype=float))
    if tf_grid is None:
        tf_grid = np.arange(0.0, float(T.min()) - 5.0, 0.5)
    best = None
    for tf in np.asarray(tf_grid, dtype=float):
        xinv = 1.0 / (T - tf)
        A = np.c_[np.ones_like(xinv), -xinv]
        coef, *_ = np.linalg.lstsq(A, lf, rcond=None)
        res = float(np.sqrt(np.mean((A @ coef - lf) ** 2)))
        if best is None or res < best[0]:
            best = (res, float(tf), float(coef[0]), float(coef[1]))
    res, tf, ln_f0, ea_over_k = best
    f0 = float(np.exp(ln_f0))
    return {"f0": f0, "tau0": 1.0 / (2.0 * np.pi * f0), "ea_ev": ea_over_k * K_B_EV, "tf": tf, "residual_ln_f": res}


# --- calibration helpers -----------------------------------------------------------------
def bias_suppression(y: float, states: int = 2) -> float:
    """The small-signal susceptibility of a degenerate unit at ``y = dmu_u E / kT`` over its
    zero-bias value: ``sech^2(y / 2)`` for two states, ``3 (cosh y + 2) / (1 + 2 cosh y)^2`` for three."""
    if states == 2:
        return float(1.0 / np.cosh(0.5 * y) ** 2)
    if states == 3:
        c = np.cosh(y)
        return float(3.0 * (c + 2.0) / (1.0 + 2.0 * c) ** 2)
    raise ValueError("states must be 2 or 3")


def unit_dipole_from_bias_suppression(ratio: float, field: float, temperature: float = 300.0, states: int = 2) -> float:
    """``dmu_u`` (e.A) of a degenerate unit whose small-signal susceptibility falls to
    ``ratio`` of its zero-bias value at nominal field ``field`` (V/A): the inverse of
    :func:`bias_suppression` (closed form for two states, bisection for three)."""
    if not 0.0 < ratio < 1.0:
        raise ValueError("ratio must lie strictly between 0 and 1")
    if states == 2:
        y = 2.0 * float(np.arccosh(1.0 / np.sqrt(ratio)))
    else:
        lo, hi = 0.0, 1.0
        while bias_suppression(hi, states) > ratio:
            hi *= 2.0
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if bias_suppression(mid, states) > ratio:
                lo = mid
            else:
                hi = mid
        y = 0.5 * (lo + hi)
    return y * K_B_EV * temperature / field


def monomer_density(density_g_cm3: float, monomer_mass_g_mol: float) -> float:
    """Monomers per m^3."""
    return density_g_cm3 / monomer_mass_g_mol * AVOGADRO * 1e6


# --- the channel ---------------------------------------------------------------------------
@dataclass
class SwitchingChannel:
    """One relaxing population of switching units, plus a background.

    ``dmu_mon`` (e.A per monomer, the dipole gained along the field by the polar
    conformation; for a flip between +mu and -mu use ``2 mu``), ``de_mon`` (kcal/mol per
    monomer, polar minus ground), ``n`` (monomers per cooperatively switching unit),
    ``fraction`` (the fraction of monomers that belong to switchable units),
    ``rho_mon`` (monomers per m^3), ``tau0`` (s), ``ea_ev``, ``tf`` (K), ``alpha``
    (Cole-Cole breadth, 0 is Debye), ``eps_inf`` (everything that does not relax in the
    window), ``tan_bg`` (a frequency-independent background loss on ``eps'``),
    ``s_switch`` (the strain of the switch along the drive axis, dimensionless; +0.11 to
    +0.13 for TGTG' to TTTT along the chain, about -0.06 across it).
    """

    dmu_mon: float
    de_mon: float
    n: int
    fraction: float
    rho_mon: float
    tau0: float
    ea_ev: float
    tf: float = 0.0
    alpha: float = 0.0
    eps_inf: float = 4.0
    tan_bg: float = 0.0
    s_switch: float = 0.11
    label: str = ""
    states: int = 2  # 2: ground / polar along the field; 3: ground / polar up / polar down

    KCAL_PER_EV = 23.060548

    # unit-level quantities
    @property
    def dmu_unit(self) -> float:
        return self.n * self.dmu_mon

    @property
    def de_unit_ev(self) -> float:
        return self.n * self.de_mon / self.KCAL_PER_EV

    @property
    def n_units(self) -> float:
        return self.fraction * self.rho_mon / self.n

    def occupancies(self, field: float = 0.0, temperature: float = 300.0) -> tuple[float, float]:
        """``(x_polar, m)``: the polar fraction and the net polar fraction (the polarization over
        ``N dmu_u``) at bias ``field`` (V/A).  Two states: ``m = x``.  Three: ``x = (a + b) / Z``,
        ``m = (a - b) / Z`` with ``a, b = exp(-(dE_u -/+ dmu_u E) / kT)``, ``Z = 1 + a + b``."""
        kT = K_B_EV * temperature
        if self.states == 2:
            g = (self.de_unit_ev - self.dmu_unit * field) / kT
            x = float(1.0 / (1.0 + np.exp(np.clip(g, -700.0, 700.0))))
            return x, x
        if self.states != 3:
            raise ValueError("states must be 2 or 3")
        a = np.exp(np.clip(-(self.de_unit_ev - self.dmu_unit * field) / kT, -700.0, 700.0))
        b = np.exp(np.clip(-(self.de_unit_ev + self.dmu_unit * field) / kT, -700.0, 700.0))
        z = 1.0 + a + b
        return float((a + b) / z), float((a - b) / z)

    def polar_fraction(self, field: float = 0.0, temperature: float = 300.0) -> float:
        return self.occupancies(field, temperature)[0]

    def net_polar_fraction(self, field: float = 0.0, temperature: float = 300.0) -> float:
        return self.occupancies(field, temperature)[1]

    def delta_eps(self, field: float = 0.0, temperature: float = 300.0) -> float:
        """Dielectric strength of the channel at bias ``field`` (V/A): ``N dmu^2 (x - m^2) / (eps0 kT)``,
        which is ``x (1 - x)`` for two states."""
        x, m = self.occupancies(field, temperature)
        mu = self.dmu_unit * E_A_TO_C_M
        kT = K_B_EV * temperature * 1.602176634e-19
        return float(self.n_units * mu * mu * (x - m * m) / (EPS0 * kT))

    def dx_de(self, field: float = 0.0, temperature: float = 300.0) -> float:
        """``dx/dE`` of the polar fraction (the strain's variable) per V/m at the bias point:
        ``(dmu_u / kT) m (1 - x)``, which is ``x (1 - x)`` for two states and zero at zero bias for three."""
        x, m = self.occupancies(field, temperature)
        return float(self.dmu_unit * m * (1.0 - x) / (K_B_EV * temperature) / V_PER_A_TO_V_PER_M)

    def tau(self, temperature: float = 300.0, field: float = 0.0) -> float:
        """The relaxation time at the bias: Vogel-Fulcher, tilted.  Two states: ``tau / cosh((dmu_u E
        - dE_u) / 2kT)``.  Three: ``tau / lambda`` with ``lambda`` the slow mode of the three-state
        master equation (rates ``r g a`` and ``r b / g`` between ground and up, ``r g b`` and
        ``r a / g`` between ground and down; ``a, b = exp(+/- dmu_u E / 2kT)``, ``g = exp(-dE_u /
        2kT)``), normalised to the polarization mode of the degenerate unit at zero bias."""
        t = float(vogel_fulcher_tau(temperature, self.tau0, self.ea_ev, self.tf))
        kT = K_B_EV * temperature
        if self.states == 2:
            tilt = (self.dmu_unit * field - self.de_unit_ev) / (2.0 * kT)
            return t / float(np.cosh(np.clip(tilt, -700.0, 700.0)))
        y2 = float(np.clip(0.5 * self.dmu_unit * field / kT, -350.0, 350.0))
        a, b = np.exp(y2), np.exp(-y2)
        g = float(np.exp(np.clip(-0.5 * self.de_unit_ev / kT, -350.0, 350.0)))
        tr = (a + b) * (g + 1.0 / g)
        det = a * a + b * b + 1.0 / (g * g)
        lam = 0.5 * (tr - np.sqrt(max(tr * tr - 4.0 * det, 0.0)))
        return t / float(lam)

    def f_relax(self, temperature: float = 300.0, field: float = 0.0) -> float:
        return 1.0 / (2.0 * np.pi * self.tau(temperature, field))

    def spectral(self, freq, temperature: float = 300.0, field: float = 0.0):
        """``1 / (1 + (i w tau)^(1 - alpha))`` at the given frequencies (Hz)."""
        w = 2.0 * np.pi * np.asarray(freq, dtype=float)
        iwt = 1j * w * self.tau(temperature, field)
        return 1.0 / (1.0 + iwt ** (1.0 - self.alpha))

    def permittivity(self, freq, temperature: float = 300.0, field: float = 0.0):
        """Complex relative permittivity ``eps' - i eps''`` (loss positive)."""
        e = self.eps_inf + self.delta_eps(field, temperature) * self.spectral(freq, temperature, field)
        return e.real - 1j * (-e.imag + self.tan_bg * e.real)

    def tan_delta(self, freq, temperature: float = 300.0, field: float = 0.0):
        e = self.permittivity(freq, temperature, field)
        return -e.imag / e.real

    def strain_coefficient(self, freq, temperature: float = 300.0, field: float = 0.0):
        """Complex ``d_eff`` (strain per V/m) of the switching channel at the bias point."""
        return self.s_switch * self.dx_de(field, temperature) * self.spectral(freq, temperature, field)

    def operating_point(self, freq: float, field_bias: float, field_ac: float, young_modulus_pa: float,
                        temperature: float = 300.0, tan_delta_mech: float = 0.0) -> dict:
        """Loss, strain, work, efficiency and heat at one drive condition.

        ``field_bias`` and ``field_ac`` in V/A (0.0055 and 0.00018 are the program's 55 and
        1.8 V/um), ``young_modulus_pa`` in Pa.  Per cycle and per m^3: ``loss_electrical =
        pi eps0 eps'' E_ac^2``, ``loss_mechanical = pi Y tan_m S_ac^2``, ``work = Y S_ac^2 /
        4`` (matched load); ``efficiency = work / (work + losses)``; ``heat_w_per_m3 = f *
        losses``; ``stroke_rate = S_ac * f`` (per unit of the clutch's lever ``lambda``).
        """
        e = self.permittivity(freq, temperature, field_bias)
        eps_r, eps_i = float(e.real), float(-e.imag)
        E_ac = field_ac * V_PER_A_TO_V_PER_M
        d = self.strain_coefficient(freq, temperature, field_bias)
        s_ac = float(abs(d)) * E_ac
        loss_e = np.pi * EPS0 * eps_i * E_ac ** 2
        loss_m = np.pi * young_modulus_pa * tan_delta_mech * s_ac ** 2
        work = 0.25 * young_modulus_pa * s_ac ** 2
        total = work + loss_e + loss_m
        return {
            "frequency_hz": float(freq), "temperature_k": float(temperature),
            "eps_real": eps_r, "eps_imag": eps_i, "tan_delta": eps_i / eps_r,
            "f_relax_hz": self.f_relax(temperature, field_bias),
            "polar_fraction": self.polar_fraction(field_bias, temperature),
            "net_polar_fraction": self.net_polar_fraction(field_bias, temperature),
            "delta_eps_at_bias": self.delta_eps(field_bias, temperature),
            "d_eff_m_per_v": float(abs(d)), "strain_phase_deg": float(np.degrees(np.angle(d))),
            "strain_ac": s_ac, "blocking_stress_pa": young_modulus_pa * s_ac,
            "work_j_per_m3": float(work), "loss_electrical_j_per_m3": float(loss_e),
            "loss_mechanical_j_per_m3": float(loss_m),
            "efficiency": float(work / total) if total > 0 else 0.0,
            "heat_w_per_m3": float(freq * (loss_e + loss_m)),
            "stroke_rate_per_lambda": s_ac * float(freq),
        }

    def scaled(self, dmu_factor: float = 1.0, de_mon: float | None = None, ea_shift_ev: float = 0.0,
               tf: float | None = None, fraction: float | None = None, n: int | None = None,
               label: str = "") -> "SwitchingChannel":
        """A candidate relative to this channel: its monomer dipole scaled, its margin set,
        its barrier shifted, its cooperativity or switchable fraction replaced."""
        return replace(self, dmu_mon=self.dmu_mon * dmu_factor,
                       de_mon=self.de_mon if de_mon is None else de_mon,
                       ea_ev=self.ea_ev + ea_shift_ev, tf=self.tf if tf is None else tf,
                       fraction=self.fraction if fraction is None else fraction,
                       n=self.n if n is None else n, label=label or self.label)


def calibrate_degenerate_channel(delta_eps: float, suppression_ratio: float, field: float, dmu_mon: float,
                                 rho_mon: float, tau0: float, ea_ev: float, tf: float = 0.0,
                                 temperature: float = 300.0, states: int = 2, **kw) -> SwitchingChannel:
    """A degenerate (``dE = 0``) channel from three measured numbers and one computed one.

    ``delta_eps`` is the dielectric strength (peak ``eps'`` minus ``eps_inf``),
    ``suppression_ratio`` the small-signal ``eps`` under bias ``field`` (V/A) over its
    zero-bias value, ``dmu_mon`` the computed dipole gain per monomer, ``rho_mon`` the
    monomer density.  The bias suppression sets the unit dipole, hence the cooperative size
    ``n``; the dielectric strength then sets the switchable fraction.  ``states=3`` uses the
    three-state suppression, which gives a smaller unit for the same measurement.  ``kw``
    passes ``alpha``, ``eps_inf``, ``tan_bg``, ``s_switch``, ``label``.
    """
    dmu_u = unit_dipole_from_bias_suppression(suppression_ratio, field, temperature, states)
    n = max(1, int(round(dmu_u / dmu_mon)))
    probe = SwitchingChannel(dmu_mon=dmu_mon, de_mon=0.0, n=n, fraction=1.0, rho_mon=rho_mon,
                             tau0=tau0, ea_ev=ea_ev, tf=tf, states=states, **kw)
    fraction = delta_eps / probe.delta_eps(0.0, temperature)
    return replace(probe, fraction=float(fraction))


def literature_terpolymer(dmu_mon: float = 0.45, **kw) -> SwitchingChannel:
    """The 63/30/7 P(VDF-TrFE-CFE) terpolymer from published numbers, for calibration checks.

    Vogel-Fulcher ``tau0 = 8.08e-13 s``, ``Ea = 0.124 eV``, ``Tf = 191 K`` and a small-signal
    ``eps`` that halves under a DC bias of about 60 MV/m (Scientific Reports 2025, 330 nm
    spin-coated films, 100 Hz to 1 MHz); dielectric strength about 45 (peak ``eps'`` of
    50-60 for the terpolymers, ``eps_inf`` about 4); density 1.78 g/cm^3 and a mean monomer
    mass of 70.5 g/mol.  ``dmu_mon`` is polyfind's beta-PVDF value.  Sources in
    ``sarco/materials/computational/papers/2026-11-04_electrostriction_literature_pass.md``.
    """
    kw.setdefault("label", "P(VDF-TrFE-CFE) 63/30/7, literature")
    return calibrate_degenerate_channel(delta_eps=45.0, suppression_ratio=0.5, field=0.006, dmu_mon=dmu_mon,
                                        rho_mon=monomer_density(1.78, 70.5), tau0=8.08e-13, ea_ev=0.124, tf=191.0,
                                        temperature=300.0, **kw)
