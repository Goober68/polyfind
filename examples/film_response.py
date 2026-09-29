"""Film d33, d31, d32 of drawn, poled beta-PVDF from the crystal, through an exact laminate.

The crystal's compliance and proper (Vanderbilt) ``d`` come from the deformable path under the
Born-fitted potential with induced dipoles (``docs/ELECTROMECHANICS.md`` 5.10); the film is
:func:`polyfind.film.laminate` of crystal lamellae and passive isotropic amorphous layers
stacked along the draw axis.  The three film inputs -- crystallinity, amorphous Young's modulus
and Poisson ratio -- are *swept*, not fitted: the measured coefficients (Nix and Ward 1986,
``d33 = -32``, ``d31 = +20``, ``d32 = +1.5`` pC/N) are acceptance quantities and never enter.
``docs/ELECTROMECHANICS.md`` section 5.12 reads the output.

Run: ``PYTHONPATH=src python examples/film_response.py [--json out.json]`` (about a minute).
"""
import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

MEASURED = {"d33": -32.0, "d31": 20.0, "d32": 1.5}
CRYSTALLINITY = (0.4, 0.5, 0.6)
AMORPHOUS_YOUNG = (0.1, 0.3, 1.0, 3.0)  # GPa
AMORPHOUS_POISSON = (0.35, 0.45, 0.49, 0.499)
SINGLE_DEFICIT_CASES = ((0.4, 0.1, 0.49), (0.4, 0.3, 0.49), (0.5, 0.3, 0.49), (0.5, 1.0, 0.45), (0.6, 1.0, 0.49))


def crystal(preset: str):
    """``(S, d, e_polar_row, |P|)`` for deformable beta-PVDF, the polar rows oriented along +P."""
    import fit_born_flux as F
    from polyfind import mechanics as M
    from polyfind.fitting import FITTED_VALENCE
    from polyfind.forcefield import SimpleFF
    from polyfind.polarizability import Polarizable

    packer, params, shape, _ = F.beta_reference(SimpleFF.from_preset(preset), Polarizable())
    ref = M.Reference(packer=packer, params=params, c=float(packer.chain.c), label="beta")
    with FITTED_VALENCE.applied():
        resp = M.electromechanical_response(ref, polymer=F.BETA[0], shape=shape)
    el, pz, P = resp.elastic, resp.piezo, resp.polarization
    if tuple(el.reachable) != (0, 1, 2, 5):
        raise RuntimeError(f"expected the deformable (xx, yy, zz, xy) block, got {el.reachable}")
    if abs(P[1]) > 1e-6 * abs(P[0]):
        raise RuntimeError("expected beta's polarization along the packer's x")
    sgn = 1.0 if P[0] >= 0 else -1.0
    d = pz.d_from_e.copy()
    d[0] *= sgn  # the polar row, quoted with the poling axis along +P
    return el.S, d, sgn * pz.e[0], float(np.linalg.norm(P))


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--preset", default="pvdf-dft-valence-flux-born")
    ap.add_argument("--json", default=None)
    args = ap.parse_args()

    from polyfind.film import laminate

    S, d, e_row, absP = crystal(args.preset)
    xtal = {"d33": d[0, 0], "d32": d[0, 1], "d31": d[0, 2], "young_draw": 1.0 / S[2, 2]}
    print(f"crystal ({args.preset}): d33 {xtal['d33']:+.2f}  d32 {xtal['d32']:+.2f}  d31 {xtal['d31']:+.2f} pC/N; "
          f"E_draw {xtal['young_draw']:.0f} GPa; |P| {absP:.4f} C/m^2")
    print(f"measured film: d33 {MEASURED['d33']:+.0f}  d32 {MEASURED['d32']:+.1f}  d31 {MEASURED['d31']:+.0f} pC/N\n")
    print(f"{'X_c':>5} {'E_a GPa':>8} {'nu_a':>6} | {'d33':>7} {'d32':>7} {'d31':>7} | {'d31/-d33':>8} | {'E_draw GPa':>10}")
    rows = []
    for phi in CRYSTALLINITY:
        for Ea in AMORPHOUS_YOUNG:
            for nu in AMORPHOUS_POISSON:
                f = laminate(S, d, phi, Ea, nu)
                c = f.film_coefficients(polar=0)
                ratio = c["d31"] / -c["d33"]
                row = {"crystallinity": phi, "amorphous_young": Ea, "amorphous_poisson": nu, **c,
                       "d31_over_minus_d33": ratio, "young_draw": f.young_draw}
                rows.append(row)
                print(f"{phi:5.2f} {Ea:8.2f} {nu:6.3f} | {c['d33']:+7.2f} {c['d32']:+7.2f} {c['d31']:+7.2f} | "
                      f"{ratio:8.3f} | {f.young_draw:10.3f}")
        print()
    d31 = [r["d31"] for r in rows]
    print(f"d31 over the sweep: {min(d31):+.2f} to {max(d31):+.2f} pC/N (crystal {xtal['d31']:+.2f}); "
          f"d33 {min(r['d33'] for r in rows):+.2f} to {max(r['d33'] for r in rows):+.2f}")
    ratios = [r["d31_over_minus_d33"] for r in rows]
    print(f"d31/-d33 over the sweep: {min(ratios):.3f} to {max(ratios):.3f}; measured "
          f"{MEASURED['d31'] / -MEASURED['d33']:.3f} (the ratio does not depend on the crystal's magnitude)")

    # One deficit or two?  Raise only the crystal's polar column e_x,xx until the film d33 is the
    # measured one, and read off what the film's d31 and d32 become.  Film d is linear in e_x,xx,
    # whose unit change adds S[0] (times the unit conversion) to the crystal's polar d row.
    print("\nraising only the crystal polar column e_x,xx until the film d33 is the measured -32:")
    print(f"{'X_c':>5} {'E_a GPa':>8} {'nu_a':>6} | {'e_x,xx':>7} {'x ours':>6} | {'d33':>6} {'d32':>6} {'d31':>6} | {'E_draw GPa':>10}")
    k = 1e3  # C/m^2 / GPa -> pC/N
    single = []
    for phi, Ea, nu in SINGLE_DEFICIT_CASES:
        base = laminate(S, d, phi, Ea, nu).film_coefficients(polar=0)
        d_unit = d.copy()
        d_unit[0] += S[0] * k
        slope = laminate(S, d_unit, phi, Ea, nu).film_coefficients(polar=0)["d33"] - base["d33"]
        de = (MEASURED["d33"] - base["d33"]) / slope
        d_new = d.copy()
        d_new[0] += de * S[0] * k
        f = laminate(S, d_new, phi, Ea, nu)
        c = f.film_coefficients(polar=0)
        exx = float(e_row[0] + de)
        single.append({"crystallinity": phi, "amorphous_young": Ea, "amorphous_poisson": nu,
                       "e_xx_required": exx, "factor": exx / float(e_row[0]), **c, "young_draw": f.young_draw})
        print(f"{phi:5.2f} {Ea:8.2f} {nu:6.3f} | {exx:+7.3f} {exx / e_row[0]:6.1f} | {c['d33']:+6.1f} {c['d32']:+6.2f} "
              f"{c['d31']:+6.1f} | {f.young_draw:10.2f}")
    print(f"measured{'':>16} | {'':>14} | {MEASURED['d33']:+6.1f} {MEASURED['d32']:+6.2f} {MEASURED['d31']:+6.1f} |")
    if args.json:
        with open(args.json, "w") as fh:
            json.dump({"preset": args.preset, "crystal": xtal, "crystal_e_polar_row": e_row.tolist(), "absP": absP,
                       "crystal_S": S.tolist(), "crystal_d": d.tolist(), "measured": MEASURED, "rows": rows,
                       "single_deficit": single}, fh, indent=1)


if __name__ == "__main__":
    main()
