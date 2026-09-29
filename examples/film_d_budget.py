"""What crystal ``e`` would reproduce the measured film ``d`` with this crystal's own stiffness?

``d = e S`` inverts to ``e = d C``: given the measured film coefficients (Nix and Ward 1986,
poled drawn beta film: d33 = -32, d32 = +1.5, d31 = +20 pC/N; 3 = poling = the packer's polar
x, 1 = draw = the chain axis z) and the model's stiffness block, this is the proper
(Vanderbilt) polar row a crystal would need.  It separates what a better crystal ``e`` could
close from what only a different compliance could, and prices the clamped-ion transverse
differences reported by the provider (``docs/REFERENCE_DATA_REQUEST.md``, 2026-09-13)
through the same compliance.  ``docs/ELECTROMECHANICS.md`` section 5.11 reads the output.

Run: ``PYTHONPATH=src python examples/film_d_budget.py [--json out.json]`` (about a minute).
"""
import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

MEASURED_FILM = {"d33": -32.0, "d32": 1.5, "d31": 20.0}  # pC/N, Nix and Ward 1986
# Provider's clamped-ion e_y,jj minus ours, their frame (P_y < 0), Vanderbilt convention,
# 2026-09-13: lateral xx +0.129, chain zz +0.123 C/m^2.  Their y is our polar x; the sign
# flips on taking the poling axis along +P.
CLAMPED_ION_DIFF_THEIR_FRAME = {"lateral": 0.129, "chain": 0.123}


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--preset", default="pvdf-dft-valence-flux-born")
    ap.add_argument("--json", default=None)
    args = ap.parse_args()

    import fit_born_flux as F
    from polyfind import mechanics as M
    from polyfind.fitting import FITTED_VALENCE
    from polyfind.forcefield import SimpleFF
    from polyfind.polarizability import Polarizable

    packer, params, shape, _ = F.beta_reference(SimpleFF.from_preset(args.preset), Polarizable())
    ref = M.Reference(packer=packer, params=params, c=float(packer.chain.c), label="beta")
    with FITTED_VALENCE.applied():
        resp = M.electromechanical_response(ref, polymer=F.BETA[0], shape=shape)
    P, el, pz = resp.polarization, resp.elastic, resp.piezo
    sgn = 1.0 if P[0] >= 0 else -1.0  # poling axis along +P
    cols = list(el.reachable)
    ix, iy, iz = cols.index(0), cols.index(1), cols.index(2)
    C = el.C[np.ix_(cols, cols)]
    S = el.S
    k = M.C_PER_M2_PER_GPA_TO_PC_PER_N

    e_ours = sgn * pz.e[0]
    d_ours = sgn * pz.d_from_e[0]
    d_meas = np.zeros(len(cols))
    d_meas[ix], d_meas[iy], d_meas[iz] = MEASURED_FILM["d33"], MEASURED_FILM["d32"], MEASURED_FILM["d31"]
    e_req = d_meas @ C / k

    # the polar column alone: d33 = -32 with every other column of e left as ours
    e_xx_d33_only = (d_meas[ix] / k - sum(e_ours[j] * S[j, ix] for j in range(len(cols)) if j != ix)) / S[ix, ix]

    de = np.zeros(len(cols))  # the provider-minus-ours clamped-ion difference, poling along +P
    de[iy] = -CLAMPED_ION_DIFF_THEIR_FRAME["lateral"]
    de[iz] = -CLAMPED_ION_DIFF_THEIR_FRAME["chain"]
    dd = de @ S * k

    out = {"preset": args.preset, "reachable_voigt": [c + 1 for c in cols], "absP": float(np.linalg.norm(P)),
           "C_GPa": C.tolist(), "S_per_GPa": S.tolist(),
           "e_polar_row_ours": e_ours.tolist(), "d_polar_row_ours": d_ours.tolist(),
           "d_polar_row_measured": d_meas.tolist(), "e_polar_row_required": e_req.tolist(),
           "e_xx_required_for_d33_alone": float(e_xx_d33_only),
           "clamped_ion_difference_e": de.tolist(), "clamped_ion_difference_d": dd.tolist()}
    names = ["xx", "yy", "zz", "xy"]
    print(f"{args.preset}: |P| = {out['absP']:.4f} C/m^2; columns {[names[cols.index(c)] for c in cols]} "
          "(x polar, z chain), poling along +P")
    print("  e ours      (C/m^2) " + "  ".join(f"{v:+8.4f}" for v in e_ours))
    print("  e required  (C/m^2) " + "  ".join(f"{v:+8.4f}" for v in e_req))
    print("  d ours       (pC/N) " + "  ".join(f"{v:+8.3f}" for v in d_ours))
    print("  d measured   (pC/N) " + "  ".join(f"{v:+8.3f}" for v in d_meas))
    print(f"  e_xx for d33 = {d_meas[ix]:+.0f} alone, other columns ours: {e_xx_d33_only:+.4f} C/m^2")
    print(f"  clamped-ion transverse difference through S: d33 {dd[ix]:+.3f}, d31 {dd[iz]:+.3f} pC/N")
    print(f"  required chain-axis e is {e_req[iz] / out['absP']:.0f} x |P|; C_33 = {C[iz, iz]:.0f} GPa")
    if args.json:
        with open(args.json, "w") as fh:
            json.dump(out, fh, indent=1)


if __name__ == "__main__":
    main()
