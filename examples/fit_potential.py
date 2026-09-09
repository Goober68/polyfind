"""Fit the built-in potential to the experimental crystal data, and check whether it generalises.

Reproduces the fit of :mod:`polyfind.fitting`: five parameters of ``SimpleFF``
(eps_r, charge_scale, the LJ minimum distances of H and F, and the one-fold torsion
coefficient) against the cells and densities of PE, beta- and alpha-PVDF, beta's
polarization, alpha's antipolar ground state and the alpha/beta near-degeneracy --
with gamma-PVDF held out of the fit entirely.

    python examples/fit_potential.py            # the full fit, ~20 minutes on 4 CPU cores
    python examples/fit_potential.py --quick    # a much smaller budget, a few minutes
    python examples/fit_potential.py --preset   # no fitting: score the recorded preset

Set ``POLYFIND_TABLE_CACHE`` to a directory to keep the screen tables between runs;
the first run pays 5-25 s per conformation to build them.

The last table it prints is the one that matters: every target's error before and
after, with the held-out row marked.  A fit that improves the fitted rows and not the
held-out one has told you that the parameterisation is too flexible for the data, which
is a result and not a failure.
"""
import argparse
import time

from polyfind import fitting as F


def rms_error(preds, keys):
    """RMS relative error (%) of the cell edges and the repeat, over the given structures."""
    errs = []
    for case in F.ALL_CASES:
        if case.key not in keys or case.key not in preds:
            continue
        p = preds[case.key]
        ea, eb, ec, _ = case.cell
        ea, eb = sorted((ea, eb))
        errs += [(p.axes[0] - ea) / ea, (p.axes[1] - eb) / eb, (p.c - ec) / ec]
    return 100.0 * (sum(e * e for e in errs) / len(errs)) ** 0.5


def show_chain_check(label, params):
    """The isolated-chain ranking, which the crystal objective never sees."""
    chk = F.chain_check(params)
    print(f"  {label:<12s} gauche E/bond: CH2 {chk['gauche_CH2']:+6.2f}, CF2 {chk['gauche_CF2']:+6.2f} kcal/mol; "
          f"RIS ranking: " + ", ".join(f"{n} ({e:+.2f})" for n, e, _ in chk["ranking"]))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--quick", action="store_true", help="small budget (a few minutes)")
    ap.add_argument("--preset", action="store_true", help="score the recorded preset, do not fit")
    ap.add_argument("--ablate", action="store_true",
                    help="score the two halves of the fitted vector separately and stop")
    ap.add_argument("--screen", type=int, default=16)
    ap.add_argument("--maxfev", type=int, default=70)
    args = ap.parse_args()

    print(__doc__.split("\n\n")[1].strip(), "\n")
    print("fitted targets:  " + ", ".join(f"{c.key} ({c.label})" for c in F.FITTED_CASES))
    print("held out:        " + ", ".join(f"{c.key} ({c.label})" for c in F.HELD_OUT_CASES))
    print("parameters:      " + ", ".join(F.VARIABLES))
    print()

    if args.ablate:
        print("Which half of the fitted vector transferred?  Watch the held-out column.\n")
        print(f"{'parameter set':<26s} {'objective':>9s} {'fitted RMS':>11s} {'gamma RMS':>10s} "
              f"{'gamma rho':>10s} {'beta |P|':>9s}")
        sets = {"illustrative": F.X0, **F.ABLATIONS, "full fit": F.vector_from_parameters(F.FITTED_PVDF)}
        for name, x in sets.items():
            preds = F.predict_all(F.ALL_CASES, F.parameters_from_vector(x))
            obj = sum(r.contribution for r in F.all_residuals(preds))
            print(f"{name:<26s} {obj:9.2f} {rms_error(preds, {c.key for c in F.FITTED_CASES}):10.2f}% "
                  f"{rms_error(preds, {c.key for c in F.HELD_OUT_CASES}):9.2f}% "
                  f"{preds['gamma'].density:10.3f} {preds['beta'].polarization:9.3f}")
        print("\nexperiment: gamma rho = 1.940 g/cm3, beta |P| = 0.13 C/m2")
        return

    if args.preset:
        t0 = time.time()
        fitted_params = F.FITTED_PVDF
        before = F.predict_all(F.ALL_CASES, F.ILLUSTRATIVE)
        after = F.predict_all(F.ALL_CASES, fitted_params)
        print(f"scored the recorded preset in {time.time() - t0:.0f} s")
        print(f"illustrative: {F.ILLUSTRATIVE.describe()}")
        print(f"fitted:       {F.FITTED_PVDF.describe()}")
        print(F.comparison_table(before, after))
        print("\nfitted targets only, before:")
        print(F.residual_table(before))
        print("\nfitted targets only, after:")
        print(F.residual_table(after))
    else:
        n_screen, maxfev = (4, 20) if args.quick else (args.screen, args.maxfev)
        res = F.fit(n_screen=n_screen, maxfev=maxfev, verbose=True)
        fitted_params = res.parameters
        print()
        print(f"illustrative: {F.ILLUSTRATIVE.describe()}")
        print(f"fitted:       {res.parameters.describe()}")
        print(f"objective {res.objective_before:.2f} -> {res.objective_after:.2f} "
              f"in {res.n_evaluations} objective evaluations, {res.seconds:.0f} s")
        print(res.table())
        print("\nfitted targets only, before:")
        print(F.residual_table(res.before))
        print("\nfitted targets only, after:")
        print(F.residual_table(res.after))
        before, after = res.before, res.after

    fitted = {c.key for c in F.FITTED_CASES}
    held = {c.key for c in F.HELD_OUT_CASES}
    print("\nGeneralisation (RMS relative error on the cell edges and the repeat):")
    print(f"  fitted structures  {rms_error(before, fitted):5.2f}%  ->  {rms_error(after, fitted):5.2f}%")
    print(f"  held out (gamma)   {rms_error(before, held):5.2f}%  ->  {rms_error(after, held):5.2f}%")
    print(f"  held-out density   {before['gamma'].density:5.3f}   ->  {after['gamma'].density:5.3f}"
          f"   (experiment 1.940)")
    print("\nAnd the isolated-chain physics, which no term of the objective constrains:")
    show_chain_check("illustrative", F.ILLUSTRATIVE)
    show_chain_check("fitted", fitted_params)


if __name__ == "__main__":
    main()
