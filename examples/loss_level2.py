"""Level 2 of the loss model on the baseline terpolymer: the driven lattice against level 1.

    python examples/loss_level2.py            # 32 x 32 units, the full sweep (an hour or two)
    python examples/loss_level2.py --quick    # 20 x 20, fewer periods and points (minutes)
    python examples/loss_level2.py --states 2 # the two-state unit of the first runs

The lattice is built from the level-1 channel (its unit size, dipole, density, attempt time
and Vogel-Fulcher time at 300 K) with a coupling J and a spread of margins (the random field
that keeps it a relaxor), and its barrier is calibrated by measurement to the channel's
relaxation time.  With the three-state unit (the default) the polar state's dipole points
up or down along the field, so the linear strain vanishes at zero bias and grows with it.
It then reports what level 1 could not: its own temperature dependence of the relaxation
time, the small-signal permittivity and strain per field along the bias axis from its own
fluctuations, and the driven response at the program's bias-hold point for a small and a
large AC drive, with the loop loss, the effective tan delta, the strain amplitude, phase
and second harmonic, the third-harmonic content and the per-cycle efficiency, beside level
1's Debye numbers.
"""
import argparse
import time

import numpy as np

from polyfind import kmc as KM
from polyfind import loss as L


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--L", type=int, default=None)
    ap.add_argument("--states", type=int, default=3, choices=(2, 3))
    ap.add_argument("--coupling", type=float, default=0.02, help="J per bond, eV")
    ap.add_argument("--spread", type=float, default=0.08, help="margin spread per unit, eV")
    ap.add_argument("--dipolar", type=float, default=0.0, help="J_d per bond, eV (three states)")
    args = ap.parse_args()
    L_side = args.L or (20 if args.quick else 32)
    periods, eq = (12, 3) if args.quick else (24, 4)
    T = 300.0
    ch = L.literature_terpolymer(states=args.states)
    t0 = time.time()
    m = KM.IsingKMC.from_channel(ch, L=L_side, coupling_ev=args.coupling, temperature=T, margin_spread_ev=args.spread,
                                 seed=0, samples=2000 if args.quick else 3000, dipolar_coupling_ev=args.dipolar)
    cal = m.calibration
    print(f"lattice {L_side}x{L_side}, {args.states} states, J = {args.coupling} eV, J_d = {args.dipolar} eV, margin spread "
          f"{args.spread} eV, units of {ch.n} monomers, dmu_unit {ch.dmu_unit:.1f} e.A (kT/dmu_unit = "
          f"{L.K_B_EV * T / ch.dmu_unit * 1e4:.0f} V/um); barrier: formula {cal['barrier_formula']:.3f} eV -> calibrated "
          f"{cal['barrier_ev']:.3f} eV (tau {cal['tau_before']:.2e} -> {cal['tau_after']:.2e} s, target {cal['tau_target']:.2e}; "
          f"{time.time() - t0:.0f} s)")
    print(f"zero-field polar fraction {m.polar_fraction:.3f}")
    st = m.calibrate_strength(45.0, 0.0, samples=8000 if args.quick else 16000)
    print(f"static susceptibility at zero bias {st['chi_before']:.1f} against the measured strength 45: unit density scaled by "
          f"{st['density_factor']:.3f} (the coupling's Curie-Weiss enhancement absorbed into the switchable fraction)")

    print("\nrelaxation time against temperature: the lattice's own, the channel's Vogel-Fulcher, and Arrhenius at the lattice barrier")
    temps = (280.0, 300.0, 320.0) if args.quick else (270.0, 285.0, 300.0, 315.0, 330.0)
    for Tk in temps:
        mt = m.at_temperature(Tk, seed=1)
        guess = float(L.vogel_fulcher_tau(Tk, ch.tau0, ch.ea_ev, ch.tf))
        t1 = time.time()
        tau_k = mt.autocorrelation_time(0.0, samples=6000 if args.quick else 10000, dt=guess / 10.0, relaxation=50.0 * guess)
        arr = (1.0 / m.nu0) * np.exp(m.ea / (L.K_B_EV * Tk))
        print(f"  T = {Tk:5.0f} K: lattice {tau_k:9.2e} s, Vogel-Fulcher {guess:9.2e} s, Arrhenius(Ea = {m.ea:.3f}) {arr:9.2e} s   ({time.time() - t1:.0f} s)")

    freqs = np.array((20e3, 50e3, 100e3) if args.quick else (10e3, 20e3, 34e3, 50e3, 66e3, 100e3))
    print("\nsmall-signal permittivity from the lattice's own fluctuations (fluctuation-dissipation), eps_inf = 4, against level 1:")
    for bias in (0.0, 0.0055):
        t1 = time.time()
        lr = m.linear_response(freqs, field_dc=bias, samples=12000 if args.quick else 24000, eps_inf=ch.eps_inf)
        print(f"  bias {bias * 1e4:4.1f} V/um: polar fraction {lr['x_mean']:.3f} (net {lr['m_mean']:.3f}), chi_static {lr['chi_static']:.1f} "
              f"(level 1 {ch.delta_eps(bias, T):.1f}), tau_integral {lr['tau_integral']:.2e} s (level 1 {ch.tau(T, bias):.2e}), "
              f"dx/dE {lr['dx_de']:.2e} per V/m (level 1 {ch.dx_de(bias, T):.2e})   ({time.time() - t1:.0f} s)")
        print(f"    {'f (kHz)':>8s} {'eps_r':>7s} {'eps_i':>7s} {'tan d':>8s} {'tan d (L1)':>11s}")
        for f, er, ei, td in zip(freqs, lr["eps_real"], lr["eps_imag"], lr["tan_delta"]):
            print(f"    {f / 1e3:8.0f} {er:7.1f} {ei:7.2f} {td:8.4f} {float(ch.tan_delta(f, T, bias)):11.4f}")

    print("\nbias sweep (fluctuation-dissipation): where the strain per field and the loss sit along the DC axis; 1.8 V/um AC, "
          "Y = 0.2 GPa, mechanical tan delta 0.05; S_budget is the strain amplitude at the 4 MW/m^3 heat budget (linear scaling)")
    print(f"{'bias (V/um)':>12s} {'x':>6s} {'m':>6s} {'chi_s':>7s} {'tau (s)':>9s} {'f (kHz)':>8s} {'tan d':>8s} {'S_ac':>9s} "
          f"{'work':>7s} {'loss_e':>8s} {'eff.':>6s} {'heat (MW/m3)':>13s} {'S_budget':>9s} {'S_budget*f':>11s}")
    e_ac = 0.00018
    E_ac_si = e_ac * 1e10
    for bias in (0.0, 0.001, 0.002, 0.003, 0.0045, 0.0055, 0.0065):
        lr = m.linear_response(freqs, field_dc=bias, samples=8000 if args.quick else 16000, eps_inf=ch.eps_inf)
        for f, ei, td, d in zip(freqs, lr["eps_imag"], lr["tan_delta"], lr["d_eff"]):
            s_ac = abs(d) * E_ac_si
            loss_e = np.pi * L.EPS0 * ei * E_ac_si ** 2
            work = 0.25 * 2e8 * s_ac ** 2
            loss_m = np.pi * 2e8 * 0.05 * s_ac ** 2
            eff = work / (work + loss_e + loss_m) if work + loss_e + loss_m > 0 else 0.0
            heat = f * (loss_e + loss_m)
            s_budget = s_ac * np.sqrt(4e6 / heat) if heat > 0 else float("nan")
            print(f"{bias * 1e4:12.1f} {lr['x_mean']:6.3f} {lr['m_mean']:6.3f} {lr['chi_static']:7.1f} {lr['tau_integral']:9.2e} {f / 1e3:8.0f} "
                  f"{td:8.4f} {s_ac:9.2e} {work:7.1f} {loss_e:8.1f} {eff:6.3f} {heat / 1e6:13.2f} {s_budget:9.2e} {s_budget * f:11.1f}")

    print("\nbias-hold point, 55 V/um DC, driven: the large-signal response (10 V/um AC) and the program's 1.8 V/um, against level 1")
    print(f"{'f (kHz)':>8s} {'E_ac (V/um)':>12s} {'tan d (lattice)':>16s} {'tan d (L1)':>11s} {'S_ac':>9s} {'S_ac (L1)':>10s} "
          f"{'phase':>7s} {'2f/1f':>6s} {'3rd/1st':>8s} {'drift':>7s} {'loss (J/m3)':>12s} {'eff.':>6s} {'heat (MW/m3)':>13s} {'events':>8s} {'s':>5s}")
    for f in freqs:
        for e_ac in (0.001, 0.00018):
            t1 = time.time()
            res = m.run(f, 0.0055, e_ac, periods=periods, equilibrate_periods=eq, samples_per_period=64)
            op = res.operating_point(2e8, 0.05)
            l1 = ch.operating_point(f, 0.0055, e_ac, 2e8, T, tan_delta_mech=0.05)
            print(f"{f / 1e3:8.0f} {e_ac * 1e4:12.1f} {res.tan_delta:16.4f} {l1['tan_delta']:11.4f} {res.s_ac:9.2e} {l1['strain_ac']:10.2e} "
                  f"{res.s_phase_deg:7.1f} {res.s_2f / res.s_ac if res.s_ac else float('nan'):6.2f} {res.harmonic_ratio:8.3f} "
                  f"{res.drift_per_cycle:7.3f} {res.loop_loss_j_per_m3:12.1f} {op['efficiency']:6.3f} "
                  f"{op['heat_w_per_m3'] / 1e6:13.2f} {res.events:8d} {time.time() - t1:5.0f}")
    print("level 1 is linear about the bias point; the driven lattice's tan delta at the large drive is the minor-loop loss over\n"
          "the stored energy at that amplitude; at 1.8 V/um the driven rows are noise-limited and the fluctuation-dissipation\n"
          "table above is the small-signal answer.")

    if args.states == 3:
        print("\nzero bias, driven at 10 V/um AC: the three-state strain is at twice the drive frequency")
        for f in freqs[:2]:
            res = m.run(f, 0.0, 0.001, periods=periods, equilibrate_periods=eq, samples_per_period=64)
            print(f"  {f / 1e3:4.0f} kHz: S at f {res.s_ac:.2e}, S at 2f {res.s_2f:.2e}, eps_r {res.eps_real:.1f}, tan d {res.tan_delta:.4f}")


if __name__ == "__main__":
    main()
