"""Electrostriction: the quadratic strain of a crystal in a field, and the conformational
switch that gives relaxor polymers their giant one.

Two things live here, and they are different in kind.

**The intrinsic (lattice) coefficient.**  ``eta_J = M_iiJ Et_i^2`` for a nominal field along
reference axis ``i``: the even part of the zero-stress strain, which
:func:`polyfind.mechanics.free_strain` measures at ``+E, 0, -E``.  The same number is a mixed
derivative of the nominal polarization.  The fixed-voltage enthalpy ``h(eta, Et)`` that
``free_strain`` minimises has the nominal polarization ``Pt = F^-1 mu / V0`` as the conjugate
of the nominal field, so in the variables ``(sigma, Et)`` the strain is ``eta_J = S_JK sigma_K
+ d_iJ Et_i + M_iiJ Et_i^2 + ...`` and the polarization ``Pt_i = P_i + d_iJ sigma_J + eps0
chi_ii Et_i + 2 M_iiJ sigma_J Et_i + ...``: the coefficient of ``sigma_J Et_i`` in ``Pt_i`` is
``2 M_iiJ``.  :func:`electrostriction` computes both -- the strain route and the polarization
route, from the state at ``(+/- stress, +/- field)`` -- and reports them side by side.  They
agree because they are two second derivatives of one function, and when they do not the finite
differences are the thing to look at.  (The textbook form ``2M = d(eps0 eps)/d sigma`` is this
relation with the *fixed-stress* permittivity; the fixed-strain one that a DFPT run at a
strained cell gives differs from it by the piezoelectric term ``e d / eps0``, which for a
rigid fixed-charge beta-PVDF is half of ``chi`` and not a correction.)  ``Q = M / (eps0
chi^sigma)^2`` converts to the polarization form with the free (zero-stress) susceptibility,
the one a free film's dielectric constant measures.  Rigid or deformable path, fixed charges
or induced dipoles, as the :class:`~polyfind.mechanics.Reference` was built: with fixed charges
and a rigid chain the only response is the setting angles and the offset turning in the
field, so ``chi`` is the orientational part alone and ``Q`` is correspondingly large for a
small ``M``.

This is the floor of the response, what is left above the relaxor's relaxation frequency.
Dipole rotation in the ferroelectric phase gives about 1% strain and ``Q33 < 3 m^4/C^2`` in
Bernholc's DFT (Chen et al., Science 375, 2022), and that is the order this computes.

**The switch.**  The giant strain of the relaxor terpolymers is a field-induced change of
conformation, TG or T3G to all-trans, shared by the crystals and the oriented amorphous
fraction (``docs/ELECTROSTRICTION_RESEARCH.md`` section 2).  :class:`TwoState` is the smallest
model of it: a unit of ``n`` monomers that costs ``n dE`` to switch and gains ``n dmu`` of dipole
along the field, so the threshold field is ``E* = dE / dmu`` and the switched fraction is a
Boltzmann two-state.  :func:`transition_strain` is the geometry of the switch, the rise per
monomer of one periodic chain over another, which the chain builder knows exactly.  The
numbers that feed it (``dE``, ``dmu``) are the lattice energies and polarizations
``polyfind.pack`` and ``polyfind.refine`` already compute, and the plan says to what resolution
they are needed and how they are to be verified.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import mechanics as M
from .mechanics import Elastic, Reference, Shape, _state_at_strain, deformable_state
from .polarizability import EPS0

K_B_EV = 8.617333262e-5  # eV/K
KCAL_PER_EV = 23.060548
V_PER_A_TO_V_PER_M = M.V_PER_A_TO_V_PER_M  # 1e10
E_A_PER_C_M = 1.0 / 1.602176634e-29  # e.A per C.m


# --- the intrinsic coefficient -------------------------------------------------------------
@dataclass
class Electrostriction:
    """The field electrostriction tensor of one crystal, by two routes.

    ``M_strain`` and ``M_polarization`` are (3, 6) in m^2/V^2: row ``i`` is the nominal field
    along reference axis ``i`` (x polar for beta-PVDF, y lateral, z chain), column ``J`` the
    Voigt strain; ``nan`` where the parametrisation cannot strain.  ``M_strain`` is the even
    part of the zero-stress strain over ``E^2``; ``M_polarization`` is half the mixed
    derivative of the nominal polarization with respect to stress and field.  ``Q`` (3, 6) in
    m^4/C^2 is ``M_strain / (eps0 chi_free_ii)^2``, ``nan`` when there is no polarization
    response to divide by.  ``chi_free`` and ``chi_clamped`` (3,) are the relaxed-ion nominal
    susceptibilities (``eps - 1``) at zero stress and at zero strain; their difference is the
    piezoelectric term.  ``linear`` (3, 6) in pm/V is the odd part of the same strains,
    ``d eta_J / d Et_i``: the converse piezoelectric coefficient, reported so that it can be
    held against :func:`polyfind.mechanics.piezoelectric`.  ``field`` (V/A) and ``stress``
    (GPa) are the steps used.
    """

    M_strain: np.ndarray
    M_polarization: np.ndarray
    Q: np.ndarray
    chi_free: np.ndarray
    chi_clamped: np.ndarray
    linear: np.ndarray
    field: float
    stress: float
    reachable: tuple

    @property
    def max_relative_disagreement(self) -> float:
        """Largest ``|M_strain - M_polarization|`` over ``max |M_strain|`` on the reachable columns."""
        idx = list(self.reachable)
        a, b = self.M_strain[:, idx], self.M_polarization[:, idx]
        if not np.isfinite(b).any():
            return float("nan")
        scale = float(np.nanmax(np.abs(a)))
        if scale < 1e-28:  # no electrostriction to speak of (PE with no charges): both routes are noise
            return 0.0 if float(np.nanmax(np.abs(b))) < 1e-28 else float("inf")
        return float(np.nanmax(np.abs(a - b)) / scale)

    def table(self) -> str:
        cols = [f"{J + 1}" for J in self.reachable]
        out = [f"field {self.field:g} V/A ({self.field * V_PER_A_TO_V_PER_M / 1e6:.0f} MV/m); "
               f"stress step {self.stress:g} GPa; columns Voigt {' '.join(cols)}",
               "M (1e-18 m^2/V^2), strain route / polarization route:"]
        for i, r in enumerate("xyz"):
            ff = "  ".join(f"{self.M_strain[i, J] * 1e18:+10.5f}" for J in self.reachable)
            mx = "  ".join(f"{self.M_polarization[i, J] * 1e18:+10.5f}" for J in self.reachable)
            out.append(f"  E_{r}: {ff}   |   {mx}")
        out.append(f"chi (relaxed-ion, eps - 1), free: diag({self.chi_free[0]:.4f}, {self.chi_free[1]:.4f}, "
                   f"{self.chi_free[2]:.4f}); clamped: diag({self.chi_clamped[0]:.4f}, {self.chi_clamped[1]:.4f}, "
                   f"{self.chi_clamped[2]:.4f})")
        out.append("Q (m^4/C^2) from M_strain and chi_free:")
        for i, r in enumerate("xyz"):
            out.append(f"  E_{r}: " + "  ".join(f"{self.Q[i, J]:+9.3f}" for J in self.reachable))
        out.append(f"max relative disagreement between the routes: {self.max_relative_disagreement:.3g}")
        return "\n".join(out)


TIGHT = {"gtol": 1e-10, "ftol": 1e-15}  # the rigid path's internal relaxation, for the polarization route


def _state(ref: Reference, eps6, field, shape, tight: dict | None = None):
    if shape is None:
        return _state_at_strain(ref, eps6, True, field=field, **(tight or {}))
    return deformable_state(ref, shape, eps6, field=field)


def state_at_stress(ref: Reference, elastic: Elastic, target6, field=None, shape: Shape | None = None,
                    iterations: int = 12, tol: float = 1e-8, step_tol: float = 1e-12, tight: dict | None = None):
    """The strained, internally relaxed state at applied stress ``target6`` (GPa) and nominal field.

    :func:`polyfind.mechanics.free_strain` with a stress target: Newton on the analytic stress
    with the compliance as Jacobian, on the reachable components only.  Returns ``(eps6,
    state)``; the state's ``nominal`` is the polarization at that stress and field.  ``tight``
    is passed to the internal relaxation on the rigid path (:data:`TIGHT` by default): the
    stress is quadratic in a leftover gradient of that relaxation but the polarization is
    linear in it, so a stress target converges at the default tolerance while a polarization
    read through the same state does not.
    """
    tight = TIGHT if tight is None else tight
    idx = list(elastic.reachable)
    target = np.asarray(target6, dtype=float).reshape(6)
    eps = np.zeros(6)
    st = None
    for _ in range(iterations):
        st = _state(ref, eps, field, shape, tight)
        sig = st.stress[idx] - target[idx]
        if float(np.max(np.abs(sig))) < tol:
            return eps, st
        step = elastic.S @ sig
        eps[idx] = eps[idx] - step
        if float(np.max(np.abs(step))) < step_tol:
            return eps, st
    return eps, _state(ref, eps, field, shape, tight)


def nominal_susceptibility(ref: Reference, eps6, field: float = 2e-3, shape: Shape | None = None) -> np.ndarray:
    """``chi^eta_ij = (1/eps0) d Pt_i / d Et_j`` at fixed strain ``eps6``, relaxed-ion, dimensionless.

    ``Pt = F^-1 mu / V0`` is the nominal polarization (the charge per unit *reference*
    electrode area) and ``Et`` the nominal field, the conjugate pair of the fixed-voltage
    enthalpy; at zero strain this is the relaxed-ion susceptibility of
    :func:`polyfind.mechanics.dielectric_tensor`.  Central differences in the field; the
    internal coordinates (and the chain's shape, with a :class:`~polyfind.mechanics.Shape`)
    are relaxed at every evaluation.
    """
    eps6 = np.asarray(eps6, dtype=float).reshape(6)
    chi = np.zeros((3, 3))
    for j in range(3):
        ms = []
        for sgn in (1.0, -1.0):
            f = np.zeros(3)
            f[j] = sgn * field
            ms.append(np.asarray(_state(ref, eps6, f, shape, TIGHT).nominal, dtype=float))
        chi[:, j] = (ms[0] - ms[1]) / (2.0 * field)  # C/m^2 per V/A
    return chi / V_PER_A_TO_V_PER_M / EPS0


def electrostriction(ref: Reference, elastic: Elastic | None = None, field: float = 0.005,
                     stress: float = 0.005, shape: Shape | None = None, chi_field: float = 2e-3,
                     polarization_route: bool = True) -> Electrostriction:
    """The intrinsic electrostriction of ``ref`` by the strain and the polarization routes.

    ``field`` (V/A; 0.005 is 50 MV/m) is the magnitude of the nominal field applied along
    each reference axis.  Strain route: the even part of the zero-stress strain over ``E^2``;
    the quadratic strain at 50 MV/m is of order ``1e-4`` for a real crystal, well above the
    strain noise of the zero-stress search, and the odd part cancels exactly, so the field
    need not be small; halving it is the convergence check.  Polarization route: the nominal
    polarization at the four states ``(+/- stress, +/- field)`` for each reachable stress
    component (``stress`` in GPa), mixed central difference, halved.  ``chi_field`` is the
    field step of the susceptibilities.  ``elastic`` is computed if not given.

    Measured on rigid beta-PVDF under the illustrative potential with
    ``CrystalPacker(lj_cutoff="force")``: the two routes agree to 5e-5 at ``stress=0.005``,
    8e-5 at 0.01, 2e-4 at 0.02 and 2% at 0.04 GPa (the stress-squared term), and the strain
    route moves by 2.5e-3 between ``field`` 0.005 and 0.0025.  With the package default
    ``lj_cutoff="energy"`` the Lennard-Jones force jumps at the cutoff, which beta's ``b``
    axis straddles (:class:`polyfind.mechanics.Elastic` measures the same fault on ``C_22``),
    and the stress-targeted states land on either side of the jump: the polarization route
    then scatters by up to 3x between stress steps while the strain route barely moves.
    Build the reference force-shifted before trusting either number to better than 20%.
    """
    if elastic is None:
        elastic = M.elastic_constants(ref, shape=shape)
    reach = tuple(elastic.reachable)
    idx = list(reach)
    E_si = field * V_PER_A_TO_V_PER_M
    S_si = stress * 1e9

    zero = np.zeros(6)
    eps0, _ = state_at_stress(ref, elastic, zero, np.zeros(3), shape)
    M_strain = np.full((3, 6), np.nan)
    linear = np.full((3, 6), np.nan)
    chi_free = np.zeros(3)
    for i in range(3):
        f = np.zeros(3)
        f[i] = field
        ep, _ = state_at_stress(ref, elastic, zero, f, shape)
        em, _ = state_at_stress(ref, elastic, zero, -f, shape)
        M_strain[i, idx] = (ep[idx] + em[idx] - 2.0 * eps0[idx]) / (2.0 * E_si ** 2)
        linear[i, idx] = (ep[idx] - em[idx]) / (2.0 * E_si) * 1e12
        fc = np.zeros(3)
        fc[i] = chi_field
        _, sp = state_at_stress(ref, elastic, zero, fc, shape)
        _, sm = state_at_stress(ref, elastic, zero, -fc, shape)
        chi_free[i] = (sp.nominal[i] - sm.nominal[i]) / (2.0 * chi_field) / V_PER_A_TO_V_PER_M / EPS0

    chi_clamped = np.diag(nominal_susceptibility(ref, np.zeros(6), chi_field, shape)).copy()

    M_pol = np.full((3, 6), np.nan)
    if polarization_route:
        for i in range(3):
            f = np.zeros(3)
            f[i] = field
            for J in idx:
                tgt = np.zeros(6)
                tgt[J] = stress
                _, spp = state_at_stress(ref, elastic, tgt, f, shape)
                _, spm = state_at_stress(ref, elastic, tgt, -f, shape)
                _, smp = state_at_stress(ref, elastic, -tgt, f, shape)
                _, smm = state_at_stress(ref, elastic, -tgt, -f, shape)
                mixed = (spp.nominal[i] - spm.nominal[i] - smp.nominal[i] + smm.nominal[i]) / (4.0 * S_si * E_si)
                M_pol[i, J] = 0.5 * mixed

    Q = np.full((3, 6), np.nan)
    for i in range(3):
        denom = (EPS0 * chi_free[i]) ** 2
        if denom > 0.0 and chi_free[i] > 1e-9:
            Q[i, idx] = M_strain[i, idx] / denom
    return Electrostriction(M_strain=M_strain, M_polarization=M_pol, Q=Q, chi_free=chi_free,
                            chi_clamped=chi_clamped, linear=linear, field=float(field), stress=float(stress),
                            reachable=reach)


# --- the switch ----------------------------------------------------------------------------
def threshold_field(delta_e_kcal: float, delta_mu_ea: float) -> float:
    """``E* = dE / dmu`` in V/A: the field at which a unit's polar conformation breaks even.

    ``delta_e_kcal`` is the energy of the polar conformation above the ground one in kcal/mol
    per monomer; ``delta_mu_ea`` the dipole gained along the field in e.A per monomer.  The
    ratio is intensive, so the size of the cooperatively switching unit does not enter.
    Multiply by 1e4 for MV/m.
    """
    if delta_mu_ea <= 0.0:
        raise ValueError("delta_mu must be positive: the polar conformation must gain dipole along the field")
    return delta_e_kcal / KCAL_PER_EV / delta_mu_ea


@dataclass
class TwoState:
    """A unit of ``n`` monomers switching between a non-polar and a polar conformation.

    ``delta_e`` kcal/mol per monomer (polar minus ground, zero field), ``delta_mu`` e.A per
    monomer gained along the field, ``strain`` the strain of the switch (any shape; axial,
    transverse or a Voigt vector), ``n`` the cooperative unit.  :meth:`fraction` is the
    Boltzmann occupancy of the polar state at field ``E`` (V/A) and temperature ``T``;
    :meth:`strain_at` is that fraction times the switch strain.  Zero kelvin is the step at
    :attr:`threshold`.  Nothing here knows about the barrier: this is the equilibrium curve,
    the one a slow measurement sees.
    """

    delta_e: float
    delta_mu: float
    strain: np.ndarray | float
    n: int = 1

    @property
    def threshold(self) -> float:
        return threshold_field(self.delta_e, self.delta_mu)

    def fraction(self, field: float, temperature: float = 300.0) -> float:
        """Occupancy of the polar state at nominal field ``field`` (V/A) along the dipole."""
        gap_ev = self.n * (self.delta_e / KCAL_PER_EV - self.delta_mu * field)
        if temperature <= 0.0:
            return 1.0 if gap_ev < 0.0 else (0.5 if gap_ev == 0.0 else 0.0)
        x = gap_ev / (K_B_EV * temperature)
        if x > 700.0:
            return 0.0
        return float(1.0 / (1.0 + np.exp(x)))

    def strain_at(self, field: float, temperature: float = 300.0):
        return self.fraction(field, temperature) * np.asarray(self.strain, dtype=float)

    def slope(self, field: float, temperature: float = 300.0, h: float = 1e-5):
        """``d strain / d E`` at ``field`` (per V/A): the bias-hold small-signal response."""
        return (self.strain_at(field + h, temperature) - self.strain_at(field - h, temperature)) / (2.0 * h)


def dipole_per_monomer(polarization_c_m2: float, volume_a3: float, monomers: int) -> float:
    """``P * V / n`` in e.A: the dipole of one monomer from a cell polarization and volume."""
    return float(polarization_c_m2) * float(volume_a3) * 1e-30 * E_A_PER_C_M / int(monomers)


def transition_strain(polymer, seq_ground, seq_polar, states, area_ground: float | None = None,
                      area_polar: float | None = None) -> dict:
    """The strain of switching one chain from ``seq_ground`` to ``seq_polar``.

    Axial: the rise per monomer of the polar periodic chain over the ground one, minus one
    (+0.11 to +0.13 for PVDF TGTG' to TTTT: the chain extends).  With the cross-section per
    chain of the two packed cells (A^2), the areal strain and its isotropic transverse
    equivalent ``sqrt(1 + areal) - 1`` are added.  Both chains are built at the ideal state
    angles of ``states``; refined torsions move the number by a few percent.
    """
    from .helix import helix_parameters

    hg = helix_parameters(polymer, np.array([int(s) for s in seq_ground]), states)
    hp = helix_parameters(polymer, np.array([int(s) for s in seq_polar]), states)
    rise_g = hg.rise_per_period / hg.monomers_per_period
    rise_p = hp.rise_per_period / hp.monomers_per_period
    out = {"rise_ground": float(rise_g), "rise_polar": float(rise_p), "axial": float(rise_p / rise_g - 1.0)}
    if area_ground is not None and area_polar is not None:
        areal = float(area_polar) / float(area_ground) - 1.0
        out["areal"] = areal
        out["transverse_isotropic"] = float(np.sqrt(1.0 + areal) - 1.0)
        out["volumetric"] = float((1.0 + areal) * (1.0 + out["axial"]) - 1.0)
    return out
