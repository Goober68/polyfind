"""Level 2 of the loss model on the baseline terpolymer: the driven lattice against level 1.

    python examples/loss_level2.py            # 24 x 24 units, the full sweep (tens of minutes)
    python examples/loss_level2.py --quick    # 16 x 16, fewer periods and points (a few minutes)

The lattice is built from the level-1 channel (its unit size, dipole, density, attempt time
and Vogel-Fulcher time at 300 K) with a coupling J and a spread of margins (the random field
that keeps it a relaxor), and its barrier is calibrated by measurement to the channel's
relaxation time.  It then reports what level 1 could not: its own temperature dependence of
the relaxation time, and the response at the program's bias-hold point for a small and a
large AC drive, with the loop loss, the effective tan delta, the strain amplitude and phase,
the third-harmonic content and the per-cycle efficiency, beside level 1's Debye numbers.
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
    ap.add_argument("--coupling", type=float, default=0.02, help="J per bond, eV")
    ap.add_argument("--spread", type=float, default=0.08, help="margin spread per unit, eV")
    args = ap.parse_args()
    L_side = args.L or (20 if args.quick else 32)
    periods, eq = (12, 3) if args.quick else (24, 4)
    T = 300.0
    ch = L.literature_terpolymer()
    t0 = time.time()
    m = KM.IsingKMC.from_channel(ch, L=L_side, coupling_ev=args.coupling, temperature=T, margin_spread_ev=args.spread,
                                 seed=0, samples=2000 if args.quick else 3000)
    cal = m.calibration
    print(f"lattice {L_side}x{L_side}, J = {args.coupling} eV, margin spread {args.spread} eV, units of {ch.n} monomers, "
          f"dmu_unit {ch.dmu_unit:.1f} e.A; barrier: formula {cal['barrier_formula']:.3f} eV -> calibrated {cal['barrier_ev']:.3f} eV "
          f"(tau {cal['tau_before']:.2e} -> {cal['tau_after']:.2e} s, target {cal['tau_target']:.2e}; {time.time() - t0:.0f} s)")
    print(f"zero-field polar fraction {m.polar_fraction:.3f}")
    st = m.calibrate_strength(45.0, 0.0, samples=8000 if args.quick else 16000)
    print(f"static susceptibility at zero bias {st['chi_before']:.1f} against the measured strength 45: unit density scaled by "
          f"{st['density_factor']:.3f} (the coupling's Curie-Weiss enhancement absorbed into the switchable fraction)")

    print("\nrelaxation time against temperature: the lattice's own, the channel's Vogel-Fulcher, and Arrhenius at the lattice barrier")
    temps = (280.0, 300.0, 320.0) if args.quick else (270.0, 285.0, 300.0, 315.0, 330.0)
    for Tk in temps:
        mt = KM.IsingKMC(L=L_side, dmu_unit=m.dmu, margins_ev=m.margins, coupling_ev=(m.jx, m.jy), barrier_ev=m.ea,
                         attempt_rate=m.nu0, temperature=Tk, n_units=m.n_units, s_switch=m.s_switch, seed=1)
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
        print(f"  bias {bias * 1e4:4.1f} V/um: polar fraction {lr['x_mean']:.3f}, chi_static {lr['chi_static']:.1f} "
              f"(level 1 {ch.delta_eps(bias, T):.1f}), tau_integral {lr['tau_integral']:.2e} s (level 1 {ch.tau(T, bias):.2e})   ({time.time() - t1:.0f} s)")
        print(f"    {'f (kHz)':>8s} {'eps_r':>7s} {'eps_i':>7s} {'tan d':>8s} {'tan d (L1)':>11s}")
        for f, er, ei, td in zip(freqs, lr["eps_real"], lr["eps_imag"], lr["tan_delta"]):
            print(f"    {f / 1e3:8.0f} {er:7.1f} {ei:7.2f} {td:8.4f} {float(ch.tan_delta(f, T, bias)):11.4f}")

    print("\nbias sweep (fluctuation-dissipation): where the strain per field and the loss sit along the DC axis; 1.8 V/um AC, "
          "Y = 0.2 GPa, mechanical tan delta 0.05")
    print(f"{'bias (V/um)':>12s} {'x':>6s} {'chi_s':>7s} {'tau (s)':>9s} {'f (kHz)':>8s} {'tan d':>8s} {'S_ac':>9s} {'work':>7s} "
          f"{'loss_e':>8s} {'eff.':>6s} {'heat (MW/m3)':>13s}")
    e_ac = 0.00018
    E_ac_si = e_ac * 1e10
    for bias in (0.0, 0.001, 0.002, 0.003, 0.0045, 0.0055):
        lr = m.linear_response(freqs, field_dc=bias, samples=8000 if args.quick else 16000, eps_inf=ch.eps_inf)
        dx_de = L.EPS0 * lr["chi_static"] / (m.n_units * m.dmu * L.E_A_TO_C_M)  # per V/m, static
        for f, er, ei, td in zip(freqs, lr["eps_real"], lr["eps_imag"], lr["tan_delta"]):
            roll = (er - ch.eps_inf) / lr["chi_static"] if lr["chi_static"] > 0 else 0.0  # the channel's in-phase fraction at f
            s_ac = m.s_switch * dx_de * E_ac_si * max(roll, 0.0)
            loss_e = np.pi * L.EPS0 * ei * E_ac_si ** 2
            work = 0.25 * 2e8 * s_ac ** 2
            loss_m = np.pi * 2e8 * 0.05 * s_ac ** 2
            eff = work / (work + loss_e + loss_m) if work + loss_e + loss_m > 0 else 0.0
            print(f"{bias * 1e4:12.1f} {lr['x_mean']:6.3f} {lr['chi_static']:7.1f} {lr['tau_integral']:9.2e} {f / 1e3:8.0f} {td:8.4f} "
                  f"{s_ac:9.2e} {work:7.1f} {loss_e:8.1f} {eff:6.3f} {f * (loss_e + loss_m) / 1e6:13.2f}")

    print("\nbias-hold point, 55 V/um DC, driven: the large-signal response (10 V/um AC) and the program's 1.8 V/um, against level 1")
    print(f"{'f (kHz)':>8s} {'E_ac (V/um)':>12s} {'tan d (lattice)':>16s} {'tan d (L1)':>11s} {'S_ac':>9s} {'S_ac (L1)':>10s} "
          f"{'phase':>7s} {'3rd/1st':>8s} {'drift':>7s} {'loss (J/m3)':>12s} {'eff.':>6s} {'heat (MW/m3)':>13s} {'events':>8s} {'s':>5s}")
    for f in freqs:
        for e_ac in (0.001, 0.00018):
            t1 = time.time()
            res = m.run(f, 0.0055, e_ac, periods=periods, equilibrate_periods=eq, samples_per_period=64)
            op = res.operating_point(2e8, 0.05)
            l1 = ch.operating_point(f, 0.0055, e_ac, 2e8, T, tan_delta_mech=0.05)
            print(f"{f / 1e3:8.0f} {e_ac * 1e4:12.1f} {res.tan_delta:16.4f} {l1['tan_delta']:11.4f} {res.s_ac:9.2e} {l1['strain_ac']:10.2e} "
                  f"{res.s_phase_deg:7.1f} {res.harmonic_ratio:8.3f} {res.drift_per_cycle:7.3f} {res.loop_loss_j_per_m3:12.1f} {op['efficiency']:6.3f} "
                  f"{op['heat_w_per_m3'] / 1e6:13.2f} {res.events:8d} {time.time() - t1:5.0f}")
    print("level 1 is linear about the bias point; the driven lattice's tan delta at the large drive is the minor-loop loss over\n"
          "the stored energy at that amplitude; at 1.8 V/um the driven rows are noise-limited and the fluctuation-dissipation\n"
          "table above is the small-signal answer.")


if __name__ == "__main__":
    main()
