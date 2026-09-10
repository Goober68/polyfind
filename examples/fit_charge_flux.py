"""Fit the charge-flux coefficients to a reference dipole-strain response, and audit the fit.

    python examples/fit_charge_flux.py --data <dir>          # the fit and its diagnostics
    python examples/fit_charge_flux.py --data <dir> --beta   # and beta-PVDF's d33/d31 with it

**What the reference data is, and what it is not.**  ``--data`` points at
``results/field_neighborhood_refined`` of the sibling ``sarco/materials/gpu_bundle``
project: for each of four chemistries (PVDF, VDCN, AN, CNEPO) a relaxed geometry at zero
axial strain and one at +2%, both at zero field, each carrying an absolute dipole in e.A.
Differencing the two gives ``dmu/d(strain)``, twelve numbers in all (four systems, three
Cartesian components).  Its own ``interpretation`` field reads:

    "Exploratory GFN2 prediction; finite-size, packing, stereochemistry and higher-level
    DFT validation outstanding.  Strain is relative to fixed seed span, not a stress-free
    bulk lattice."

It is a **finite two-chain pair in vacuum with the terminal backbone atoms pinned**, and
the method is GFN2-xTB, not the PBE-D3 the rest of this package's potential is fitted to.
So it calibrates a *mechanism* and an order of magnitude, not a bulk coefficient, and every
number fitted to it inherits that.  Say so wherever one is quoted.

**How the fit works.**  The flux is evaluated on the reference's *own* geometries, so the
only thing being asked of it is "do geometry-dependent charges reproduce the dipole change
these geometries produce?"  The fixed-increment model already answers part of that -- those
geometries relax, and moving fixed charges moves a dipole -- so the fit targets the
*residual*, and the residual is what the flux is fitted to.  At fixed geometry the dipole is
exactly linear in every flux coefficient, so this is an ordinary linear least squares and
the columns below are exact sensitivities rather than finite differences.

**Two channels, and only one of them is usable by a crystal.**  ``k_angle`` drives the
increment from the angles the bond makes at its own atom; ``k_bond`` drives it from the
bond's own length.  ``polyfind.chain.build_chain`` places every atom at the polymer's own
bond length, so ``k_bond`` contributes a constant charge shift with **zero gradient** inside
a packing -- measurably: with ``k_bond = 1`` on either bond type, beta's dipole holds every
digit over the whole shape sweep, exactly as the fixed model does.  It is fitted anyway,
and that is the point: without it the angle channel absorbs a stretch-driven signal it
cannot reproduce, and the fitted ``k_angle`` comes out five times larger.  Both fits are
reported below because the difference between them is the honest measure of how well this
data determines the channel a crystal can actually use.
"""
import argparse
import json
import os
import sys

import numpy as np

from polyfind.forcefield import SimpleFF, infer_bonds

SYSTEMS = ("pvdf", "vdcn", "an", "cnepo")
TAGS = ("strain+0.000_field+0", "strain+0.020_field+0")
# The ordered element pairs the four reference chemistries need.  C-C is absent on purpose:
# a homonuclear increment has no orientation (see forcefield.flux_topology), so the backbone
# channel -- the one most likely to carry an *axial* response -- is outside this model.
PAIRS = (("C", "H"), ("C", "F"), ("C", "N"), ("C", "O"))


def load(root: str) -> dict:
    """``{system: [(coords, elements, bonds, strain, dipole), ...]}`` for the two strain states."""
    out = {}
    for s in SYSTEMS:
        rows = []
        for t in TAGS:
            path = os.path.join(root, s, t + ".json")
            with open(path) as fh:
                d = json.load(fh)
            X = np.array(d["positions_A"], dtype=float)
            el = list(d["symbols"])
            bonds = infer_bonds(el, X)
            recorded = {tuple(sorted(b)) for b in d["bonds"]}
            if {tuple(sorted(b)) for b in bonds} != recorded:
                raise ValueError(f"{path}: inferred bonds disagree with the file's own list")
            rows.append((el, X, bonds, float(d["axial_strain"]), np.array(d["dipole_eA"], dtype=float)))
        out[s] = rows
    return out


def flux_ff(base: SimpleFF, coeffs: dict) -> SimpleFF:
    """``base`` with the given ``{(e_i, e_j): (k_angle, k_bond)}`` flux coefficients."""
    kw = {k: v for k, v in base.__dict__.items() if k != "_cache"}
    kw["charge_flux"] = tuple((a, b, float(ka), float(kb)) for (a, b), (ka, kb) in coeffs.items())
    return SimpleFF(**kw)


def response(data: dict, base: SimpleFF, coeffs: dict) -> np.ndarray:
    """``dmu/d(strain)`` predicted for every system, flattened to (12,)."""
    ff = base if not coeffs else flux_ff(base, coeffs)
    out = []
    for s in SYSTEMS:
        (el, X0, b0, e0, _), (_, X1, b1, e1, _) = data[s]
        m0 = ff.charges_at(el, b0, X0) @ X0
        m1 = ff.charges_at(el, b1, X1) @ X1
        out.append((m1 - m0) / (e1 - e0))
    return np.concatenate(out)


def target(data: dict) -> np.ndarray:
    return np.concatenate([(data[s][1][4] - data[s][0][4]) / (data[s][1][3] - data[s][0][3])
                           for s in SYSTEMS])


def columns(data: dict, base: SimpleFF, keys) -> tuple:
    """``(base response, sensitivity matrix (12, len(keys)))``; exact, the model being linear."""
    b = response(data, base, {})
    cols = []
    for k in keys:
        coeffs = {p: (1.0 if k == ("angle", p) else 0.0, 1.0 if k == ("bond", p) else 0.0)
                  for p in PAIRS}
        cols.append(response(data, base, coeffs) - b)
    return b, (np.stack(cols, axis=1) if cols else np.zeros((12, 0)))


def fit(A: np.ndarray, r: np.ndarray) -> tuple:
    k, *_ = np.linalg.lstsq(A, r, rcond=None)
    res = r - A @ k
    ss = float((r ** 2).sum())
    return k, float(np.sqrt((res ** 2).mean())), (1.0 - float((res ** 2).sum()) / ss if ss else 0.0)


VARIANTS = {
    "angle only (C-H, C-F)": [("angle", PAIRS[0]), ("angle", PAIRS[1])],
    "angle + bond (C-H, C-F)": [("angle", p) for p in PAIRS[:2]] + [("bond", p) for p in PAIRS[:2]],
    "angle only, four pairs": [("angle", p) for p in PAIRS],
    "angle + bond, four pairs": [("angle", p) for p in PAIRS] + [("bond", p) for p in PAIRS],
}


def report(data: dict, base: SimpleFF) -> dict:
    t = target(data)
    b0, _ = columns(data, base, [])
    r = t - b0
    print("Reference dmu/d(strain), e.A per unit strain (four systems x xyz):")
    for n, s in enumerate(SYSTEMS):
        sl = slice(3 * n, 3 * n + 3)
        print(f"  {s:6s} reference {np.round(t[sl], 3)}   fixed increments {np.round(b0[sl], 3)}"
              f"   residual {np.round(r[sl], 3)}")
    print(f"  rms residual with fixed increments: {np.sqrt((r ** 2).mean()):.4f}")

    print("\nFits (12 observations; the parameter count is the honest denominator):")
    out = {}
    for name, keys in VARIANTS.items():
        _, A = columns(data, base, keys)
        k, rms, r2 = fit(A, r)
        out[name] = dict(zip(keys, k))
        print(f"  {name:26s} {len(keys)} params ({12 / len(keys):.0f}:1)  rms {rms:6.3f}  R2 {r2:+.3f}")
        for key, v in zip(keys, k):
            print(f"      k_{key[0]:<5s} {key[1][0]}-{key[1][1]}  = {v:+.4f}")

    print("\nLeave one system out (the same fit with each chemistry held back):")
    for name in ("angle only (C-H, C-F)", "angle + bond (C-H, C-F)"):
        keys = VARIANTS[name]
        _, A = columns(data, base, keys)
        print(f"  {name}")
        for n, s in enumerate(SYSTEMS):
            m = np.ones(12, dtype=bool)
            m[3 * n : 3 * n + 3] = False
            k, _, _ = fit(A[m], r[m])
            out[f"{name} / without {s}"] = dict(zip(keys, k))
            print(f"    without {s:6s}: " + "  ".join(
                f"k_{key[0]}({key[1][0]}-{key[1][1]})={v:+.4f}" for key, v in zip(keys, k)))

    print("\nWhat else would explain the same residual, for scale:")
    scaled = SimpleFF(**{**{kk: vv for kk, vv in base.__dict__.items() if kk != "_cache"},
                         "charge_scale": base.charge_scale * 2.0})
    scale_col = response(data, scaled, {}) - b0
    k, rms, r2 = fit(scale_col[:, None], r)
    print(f"  a plain charge-magnitude scale, 1 param: charge_scale x {1 + k[0]:.3f}"
          f"  rms {rms:6.3f}  R2 {r2:+.3f}")
    keys = [("bond", p) for p in PAIRS[:2]]
    _, A = columns(data, base, keys)
    k, rms, r2 = fit(A, r)
    print(f"  the bond channel alone, 2 params:        rms {rms:6.3f}  R2 {r2:+.3f}"
          f"  (inert in a rigid-bonded crystal)")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True, help="path to results/field_neighborhood_refined")
    ap.add_argument("--preset", default="pvdf-dft-valence")
    ap.add_argument("--beta", action="store_true", help="also run beta-PVDF with each fit")
    args = ap.parse_args()

    import polyfind.fitting  # noqa: F401 - registers the presets

    base = SimpleFF.from_preset(args.preset)
    data = load(args.data)
    fits = report(data, base)
    if not args.beta:
        return

    from polyfind import mechanics as M
    from polyfind.fitting import FITTED_VALENCE
    from polyfind.polymers import PVDF

    print("\n" + "=" * 104)
    print("beta-PVDF with each fit.  Film axes: 3 = poling = the packer's x, 1 = draw = the chain")
    print("axis z, so the film's d_33 is d_x,xx and its d_31 is d_x,zz.  Measured d_33 = -32,")
    print("d_31 = +20 pC/N (Nix and Ward 1986).  Signs are quoted with the poling axis along +P.")
    print("'proper' is e S, the piezoelectric constant; 'total' adds the dimensional (thickness)")
    print("term a Broadhurst-Davis reading would add.  Neither is a prediction of experiment: the")
    print("flux is calibrated to an exploratory GFN2 finite-oligomer response.")
    print("=" * 104)
    # Pack and refine once and reuse: only the packer's charge model differs between rows.
    with FITTED_VALENCE.applied():
        ref0, rr = M.refined_reference(PVDF, [0, 0], label="beta", valence=base,
                                       refine_kw={"valence": base}, pack_kw={"table_cache_dir": None})
    chain0, params0 = ref0.packer.chain, ref0.params
    print(f"{'fit':38s} {'c':>7s} {'P_x':>8s} {'e_x,zz':>9s} {'C_33':>7s} {'sig_zz':>8s} "
          f"{'d33':>8s} {'d31':>8s} {'d33 tot':>8s} {'d31 tot':>8s} {'eS/direct':>10s}")
    print("  (sig_zz is the residual axial stress at the relaxed reference: a row whose value is "
          "not ~0 is\n   not a measurement -- the line group's +/-8 deg cap bound and the "
          "reference is not stress-free.)")
    rows = [("no flux", None)]
    # A scale sweep on the shipped fit, which is the check that the *sign* is a property of
    # the mechanism and the magnitude is just proportional to a poorly determined number.
    shipped = fits["angle only (C-H, C-F)"]
    for f in (0.25, 0.5, 2.0):
        rows.append((f"angle only x {f:g}", {k: f * v for k, v in shipped.items()}))
    rows += list(fits.items())
    for name, coeffs in rows:
        cf = None
        if coeffs is not None:
            per_pair = {p: (coeffs.get(("angle", p), 0.0), coeffs.get(("bond", p), 0.0)) for p in PAIRS}
            cf = flux_ff(base, per_pair)
        with FITTED_VALENCE.applied():
            ref = M.reference_from_chain(chain0, params0, n_chains=2, label="beta",
                                         valence=base, charge_flux=cf)
            shape = M.shape_of(PVDF, ref, angles=rr.angles)
            resp = M.electromechanical_response(ref, shape=shape)
        cols = list(resp.elastic.reachable)
        ix, ixx, izz = 0, cols.index(0), cols.index(2)
        sgn = 1.0 if resp.polarization[0] >= 0 else -1.0
        d, di = resp.piezo.d_from_e, resp.piezo.d_improper
        print(f"{name:38s} {resp.reference.c:7.4f} {resp.polarization[0]:+8.4f} "
              f"{resp.piezo.e[ix, izz]:+9.5f} {resp.elastic.C[2, 2]:7.1f} "
              f"{resp.elastic.residual_stress[2]:+8.2e} "
              f"{sgn * d[ix, ixx]:+8.2f} {sgn * d[ix, izz]:+8.2f} "
              f"{sgn * di[ix, ixx]:+8.2f} {sgn * di[ix, izz]:+8.2f} "
              f"{resp.piezo.relative_difference * 100:9.2f}%")
        sys.stdout.flush()


if __name__ == "__main__":
    main()
