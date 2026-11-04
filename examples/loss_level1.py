"""Level 1 of the loss model on the baseline terpolymer: tan delta against frequency and bias,
and the efficiency at the program's bias-hold point.

    python examples/loss_level1.py

The channel is calibrated on published numbers (Vogel-Fulcher parameters, the dielectric
strength, the halving of eps under 60 MV/m) with polyfind's beta-PVDF monomer dipole; the
bench sweep (E4 of the ask) replaces every one of those numbers for the actual film.  The
second table compares the VF relaxation with an Arrhenius one carrying polyfind's computed
kink-row barrier (0.357 eV): both put the room-temperature relaxation in the same decade
and they differ in temperature, which is what the sweep in T decides.  The third is the
bias-hold operating point of ``sarco/materials/computational/objective_function_v1.md``.
"""
import numpy as np

from polyfind import loss as L

T = 300.0
base = L.literature_terpolymer()
print(f"baseline channel: n = {base.n} monomers per unit, switchable fraction {base.fraction:.2f}, "
      f"dmu_unit {base.dmu_unit:.1f} e.A, f_relax(300 K) = {base.f_relax(T) / 1e3:.0f} kHz unbiased, "
      f"{base.f_relax(T, 0.0055) / 1e3:.0f} kHz at 55 MV/m")

freqs = [1e2, 1e3, 1e4, 3e4, 5e4, 1e5, 3e5, 1e6]
print("\ntan delta at 300 K (Debye / Cole-Cole alpha = 0.3), unbiased and at 55 MV/m DC:")
cc = L.literature_terpolymer(alpha=0.3)
print(f"{'f (Hz)':>8s} {'Debye, 0':>10s} {'Debye, 55':>10s} {'CC, 0':>10s} {'CC, 55':>10s}")
for f in freqs:
    print(f"{f:8.0f} {float(base.tan_delta(f, T, 0.0)):10.4f} {float(base.tan_delta(f, T, 0.0055)):10.4f} "
          f"{float(cc.tan_delta(f, T, 0.0)):10.4f} {float(cc.tan_delta(f, T, 0.0055)):10.4f}")

print("\nf_relax against temperature: the literature Vogel-Fulcher fit and an Arrhenius channel with the computed 0.357 eV kink hop")
arr = L.SwitchingChannel(dmu_mon=base.dmu_mon, de_mon=0.0, n=base.n, fraction=base.fraction, rho_mon=base.rho_mon,
                         tau0=base.tau0, ea_ev=0.357, tf=0.0, eps_inf=base.eps_inf, label="Arrhenius, kink hop")
print(f"{'T (K)':>6s} {'VF f_relax (kHz)':>18s} {'Arrhenius (kHz)':>16s} {'tan delta at 50 kHz, VF / Arr':>30s}")
for Tk in (260.0, 280.0, 300.0, 320.0, 340.0):
    print(f"{Tk:6.0f} {base.f_relax(Tk) / 1e3:18.1f} {arr.f_relax(Tk) / 1e3:16.1f} "
          f"{float(base.tan_delta(5e4, Tk)):15.4f} / {float(arr.tan_delta(5e4, Tk)):.4f}")

print("\nbias-hold operating point (55 V/um DC, 1.8 V/um AC, Y = 0.2 GPa, tan delta_mech 0.05):")
print(f"{'f (kHz)':>8s} {'tan d':>7s} {'S_ac':>9s} {'sigma_b (MPa)':>14s} {'work (J/m3)':>12s} {'loss_e':>10s} {'eff.':>6s} {'heat (MW/m3)':>13s}")
for f in (1e3, 1e4, 3.4e4, 5e4, 6.6e4, 1e5):
    op = base.operating_point(f, 0.0055, 0.00018, 2e8, T, tan_delta_mech=0.05)
    print(f"{f / 1e3:8.1f} {op['tan_delta']:7.4f} {op['strain_ac']:9.2e} {op['blocking_stress_pa'] / 1e6:14.3f} "
          f"{op['work_j_per_m3']:12.1f} {op['loss_electrical_j_per_m3']:10.1f} {op['efficiency']:6.3f} "
          f"{op['heat_w_per_m3'] / 1e6:13.2f}")
print("the program's optimizer assumed tan delta_0 = 0.0125 (or 0.05) at 34-66 kHz and a 30 kHz f_relax; "
      "read the table against those and against the 4 MW/m3 heat budget.")
