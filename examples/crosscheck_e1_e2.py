"""Cross-check polyfind's pre-registered E1/E2 predictions against Dipole's records, quantity by quantity.

    python examples/crosscheck_e1_e2.py --polyfind docs/predictions/e1_polyfind.json --dipole docs/predictions/e1_dipole_existing.json
    python examples/crosscheck_e1_e2.py --e2 docs/predictions/e2_polyfind.json --e2-dipole <file>

Dipole's values are read in Dipole's frame (x = a lateral, y = b polar, z = chain) and mapped
to polyfind's (x polar, y lateral, z chain) by swapping x and y.  The acceptance rules are the
ask's (sarco/materials/ask-polyfind-dipole-electrostriction-verification-v1.md): E1, sign on
every compared component and magnitude within a factor 1.5 on the polar-polar components;
E2, margins within 0.1 kcal/mol per monomer, the ordering of conformations preserved, the
conformational strain within 10%.  Components whose magnitude is below ``--floor`` on both
sides are reported as "small" and not sign-tested (PE's lattice response is the case).

Dipole-side E1 file, per system: ``C_GPa.block`` (3x3 normal block, x y z), ``eps_inf.diag``,
``eps_relaxed_ion.diag`` (zero strain), and optionally ``strained``: ``{axis: {"+0.010": {"eps_diag":
[...], "P_C_m2": [...]}, "-0.010": {...}}}`` from which the derivatives are formed here.  E2 file,
per system: ``conformations: {name: {"energy_per_monomer_kcal": ..., "P_C_m2": [...],
"c_per_monomer_A": ...}}``.
"""
import argparse
import json

import numpy as np

SWAP = np.array([1, 0, 2])  # Dipole (x lateral, y polar, z chain) -> polyfind (x polar, y lateral, z chain)
AXIS_D2P = {"x": "yy (lateral)", "y": "xx (polar)", "z": "zz (chain)"}  # Dipole's strain axis -> polyfind's strained key


def to_pf_vec(v):
    return np.asarray(v, dtype=float)[SWAP]


def to_pf_mat(m):
    m = np.asarray(m, dtype=float)
    return m[np.ix_(SWAP, SWAP)]


def verdict(pf, dp, polar_polar: bool, floor: float) -> str:
    if not (np.isfinite(pf) and np.isfinite(dp)):
        return "n/a"
    if abs(pf) < floor and abs(dp) < floor:
        return "small on both sides"
    sign_ok = np.sign(pf) == np.sign(dp)
    if not sign_ok:
        return "SIGN DISAGREES"
    if polar_polar:
        ratio = abs(pf) / abs(dp) if dp != 0 else float("inf")
        return f"ratio {ratio:.2f}: " + ("within 1.5" if 1 / 1.5 <= ratio <= 1.5 else "OUTSIDE 1.5")
    return "sign agrees"


def line(name, pf, dp, pp, floor):
    print(f"  {name:36s} polyfind {pf:12.4f}   Dipole {dp:12.4f}   {verdict(pf, dp, pp, floor)}")


def crosscheck_e1(pf_all, dp_all, floor):
    global to_pf_vec, to_pf_mat
    for key, dp in dp_all.get("systems", {}).items():
        pf = pf_all.get("systems", {}).get(dp.get("polyfind_system", key))
        if pf is None:
            print(f"\n[{key}] no polyfind prediction yet")
            continue
        # a record in polyfind's own frame needs no axis swap
        same_frame = str(dp.get("frame", dp_all.get("frame", ""))).lower().startswith("polyfind")
        to_pf_vec = (lambda v: np.asarray(v, dtype=float)) if same_frame else (lambda v: np.asarray(v, dtype=float)[SWAP])
        to_pf_mat = (lambda m: np.asarray(m, dtype=float)) if same_frame else (lambda m: np.asarray(m, dtype=float)[np.ix_(SWAP, SWAP)])
        print(f"\n=== E1 {pf['label']}  <->  {dp.get('label', key)} ===")
        if dp.get("status"):
            print(f"  Dipole status: {dp['status']}")
        if "cell" in dp:
            c = dp["cell"]
            if same_frame:
                pa, pb, pc = pf["cell"]["a"], pf["cell"]["b"], pf["cell"]["c"]
                da, db, dc = c.get("a", float("nan")), c.get("b", float("nan")), c.get("c", float("nan"))
                print(f"  cell (A), polyfind / Dipole, same axis order (x polar, y lateral, z chain):")
                for name, pv, dv in (("a (polar)", pa, da), ("b (lateral)", pb, db), ("c (chain)", pc, dc)):
                    print(f"  {name:36s} polyfind {pv:12.4f}   Dipole {dv:12.4f}   {100 * (pv / dv - 1):+.1f}%")
                if "volume_A3" in c:
                    print(f"  {'volume (A^3)':36s} polyfind {pf['cell']['volume_A3']:12.3f}   Dipole {c['volume_A3']:12.3f}   {100 * (pf['cell']['volume_A3'] / c['volume_A3'] - 1):+.1f}%")
            else:
                print(f"  cell (polyfind a b c / Dipole a b c, A): {pf['cell']['a']:.3f} {pf['cell']['b']:.3f} {pf['cell']['c']:.3f} / "
                      f"{c.get('a', float('nan')):.3f} {c.get('b', float('nan')):.3f} {c.get('c', float('nan')):.3f}  "
                      "(polyfind's a is along x polar; Dipole's a is the lateral axis)")
        pf_voigt = pf["C_GPa"]["voigt"]
        pf_block = np.asarray(pf["C_GPa"]["block"], dtype=float)
        if "energy_Eh_per_cell" in dp:
            print(f"  Dipole total energy {dp['energy_Eh_per_cell']:.6f} Eh per {dp.get('monomers_per_cell', '?')} monomers (differences only are comparable; E2)")
        if "C_GPa" in dp:
            dpc = to_pf_mat(dp["C_GPa"]["block"])
            print(f"  C (GPa), Dipole status: {dp['C_GPa'].get('status', '')}")
            for I in range(2):
                for J in range(2):
                    if I + 1 in pf_voigt and J + 1 in pf_voigt:
                        pi, pj = pf_voigt.index(I + 1), pf_voigt.index(J + 1)
                        line(f"C_{I + 1}{J + 1} ({'polar' if I == 0 else 'lateral'},{'polar' if J == 0 else 'lateral'})",
                             pf_block[pi, pj], dpc[I, J], I == 0 and J == 0, floor)
            print(f"  C_33 (chain): polyfind rigid chain, none; Dipole {dpc[2, 2]:.1f}")
        for name, pfk in (("eps_inf", "eps_inf"), ("eps_relaxed_ion", "eps_relaxed_ion")):
            if name in dp:
                d = to_pf_vec(dp[name]["diag"])
                p = np.diag(np.asarray(pf[pfk], dtype=float))
                print(f"  {name}, Dipole status: {dp[name].get('status', '')}")
                for i, ax in enumerate(("x polar", "y lateral", "z chain")):
                    line(f"{name}_{ax}", p[i], d[i], i == 0, floor)
        if "strained" in dp:
            for axis, cells in dp["strained"].items():
                # a record in polyfind's frame may key the strain by polyfind's own column name
                key_pf = axis if axis in pf["strained"] else AXIS_D2P.get(axis, axis)
                if key_pf not in pf["strained"]:
                    print(f"  strain along Dipole's {axis} ({key_pf}): polyfind has no column (chain rigid)")
                    continue
                plus, minus = cells["+0.010"], cells["-0.010"]
                if "eps_diag" in plus:
                    deps_d = (to_pf_vec(plus["eps_diag"]) - to_pf_vec(minus["eps_diag"])) / 0.02
                    deps_p = np.diag(np.asarray(pf["strained"][key_pf]["d_eps_d_eta"], dtype=float))
                    print(f"  d eps_ii / d eta along {key_pf}:")
                    for i, ax in enumerate(("x polar", "y lateral", "z chain")):
                        line(f"d eps_{ax} / d eta", deps_p[i], deps_d[i], key_pf.startswith("xx") and i == 0, floor)
                if "P_C_m2" in plus:
                    dp_d = (to_pf_vec(plus["P_C_m2"]) - to_pf_vec(minus["P_C_m2"])) / 0.02
                    dp_p = np.asarray(pf["strained"][key_pf]["dP_true_d_eta_C_m2"], dtype=float)
                    # the sign of P is the orientation of each side's cell; d|P|/d eta of the polar component is the comparable number
                    p0_d = to_pf_vec(dp["P0_C_m2"]) if "P0_C_m2" in dp else (to_pf_vec(plus["P_C_m2"]) + to_pf_vec(minus["P_C_m2"])) / 2
                    s_pf = np.sign(pf["P0_C_m2"][0]) or 1.0
                    s_dp = np.sign(p0_d[0]) or 1.0
                    print(f"  polarization against strain along {key_pf} (ions relaxed at the fixed strained cell on both sides):")
                    line("d|P_polar| / d eta (true, C/m^2)", s_pf * dp_p[0], s_dp * dp_d[0], key_pf.startswith("xx"), floor)
                    if abs(pf["P0_C_m2"][0]) > floor and abs(p0_d[0]) > floor:
                        line("(d|P_polar| / d eta) / |P_polar|", s_pf * dp_p[0] / abs(pf["P0_C_m2"][0]), s_dp * dp_d[0] / abs(p0_d[0]), key_pf.startswith("xx"), floor)
                    print(f"    raw dP/d eta, polyfind {np.round(dp_p, 4).tolist()} / Dipole {np.round(dp_d, 4).tolist()} C/m^2 (sign is the orientation of the cell's dipoles)")
        if "P_C_m2" in dp:
            dpv = to_pf_vec(dp["P_C_m2"]); pfv = np.asarray(pf["P0_C_m2"], dtype=float)
            print(f"  polarization (C/m^2), Dipole status: {dp.get('status', '')}")
            line("|P|", float(np.linalg.norm(pfv)), float(np.linalg.norm(dpv)), True, floor)
            print(f"    components polyfind {np.round(pfv, 4).tolist()} / Dipole {np.round(dpv, 4).tolist()} (sign is the orientation of the cell's dipoles)")
        if "berry_slope_C_m2" in dp:
            b = dp["berry_slope_C_m2"]
            print(f"  Berry slope dP_polar/d eta_lateral (Dipole, {b.get('status', '')}): "
                  + ", ".join(f"{k} {v:+.3f}" for k, v in b["dP_y_d_eta_x"].items()))
            if "yy (lateral)" in pf["strained"]:
                dp_p = np.asarray(pf["strained"]["yy (lateral)"]["dP_true_d_eta_C_m2"], dtype=float)
                print(f"    polyfind dP_polar/d eta_lateral (true) {dp_p[0]:+.3f} C/m^2; proper e_12 "
                      f"{np.asarray(pf['e_proper_C_m2']['values'])[0][pf['e_proper_C_m2']['cols_voigt'].index(2)]:+.3f}")


def crosscheck_e2(pf_all, dp_all):
    for key, dp in dp_all.get("systems", {}).items():
        pf = pf_all.get("systems", {}).get(key)
        if pf is None:
            print(f"\n[{key}] no polyfind prediction yet")
            continue
        print(f"\n=== E2 {pf['label']}  <->  {dp.get('label', key)} ===")
        pfc, dpc = pf["conformations"], dp["conformations"]
        common = [c for c in pfc if c in dpc and "failed" not in pfc[c]]
        if len(common) < 2:
            print("  fewer than two conformations on both sides")
            continue
        g_pf = min(common, key=lambda c: pfc[c]["energy_per_monomer_kcal"])
        g_dp = min(common, key=lambda c: dpc[c]["energy_per_monomer_kcal"])
        print(f"  ground conformation: polyfind {g_pf}, Dipole {g_dp}  {'(ordering of the ground state agrees)' if g_pf == g_dp else '(GROUND STATE DISAGREES)'}")
        order_pf = sorted(common, key=lambda c: pfc[c]["energy_per_monomer_kcal"])
        order_dp = sorted(common, key=lambda c: dpc[c]["energy_per_monomer_kcal"])
        print(f"  ordering: polyfind {' < '.join(order_pf)}; Dipole {' < '.join(order_dp)}  {'agrees' if order_pf == order_dp else 'DIFFERS'}")
        for c in common:
            m_pf = pfc[c]["energy_per_monomer_kcal"] - pfc[g_dp]["energy_per_monomer_kcal"]
            m_dp = dpc[c]["energy_per_monomer_kcal"] - dpc[g_dp]["energy_per_monomer_kcal"]
            print(f"  margin {c:8s} vs {g_dp}: polyfind {m_pf:+.3f}  Dipole {m_dp:+.3f} kcal/mol per monomer  "
                  f"{'within 0.1' if abs(m_pf - m_dp) <= 0.1 else 'OUTSIDE 0.1'}")
            s_pf = pfc[c]["c_per_monomer_A"] / pfc[g_dp]["c_per_monomer_A"] - 1.0
            s_dp = dpc[c]["c_per_monomer_A"] / dpc[g_dp]["c_per_monomer_A"] - 1.0
            if abs(s_dp) > 1e-6:
                print(f"    conformational strain vs {g_dp}: polyfind {s_pf * 100:+.2f}%  Dipole {s_dp * 100:+.2f}%  "
                      f"{'within 10%' if abs(s_pf / s_dp - 1.0) <= 0.1 else 'OUTSIDE 10%'}")
            if "P_C_m2" in dpc[c]:
                print(f"    |P|: polyfind {np.linalg.norm(pfc[c]['P_C_m2']):.3f}  Dipole {np.linalg.norm(dpc[c]['P_C_m2']):.3f} C/m^2")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--polyfind", default="docs/predictions/e1_polyfind.json")
    ap.add_argument("--dipole", default="docs/predictions/e1_dipole_existing.json")
    ap.add_argument("--e2", default=None)
    ap.add_argument("--e2-dipole", default=None)
    ap.add_argument("--floor", type=float, default=0.05)
    args = ap.parse_args()
    if args.polyfind and args.dipole:
        with open(args.polyfind) as fh:
            pf = json.load(fh)
        with open(args.dipole) as fh:
            dp = json.load(fh)
        print(f"E1: polyfind {pf.get('potential', '')}; Dipole frame: {dp.get('frame', '')}")
        crosscheck_e1(pf, dp, args.floor)
    if args.e2 and args.e2_dipole:
        with open(args.e2) as fh:
            pf2 = json.load(fh)
        with open(args.e2_dipole) as fh:
            dp2 = json.load(fh)
        crosscheck_e2(pf2, dp2)


if __name__ == "__main__":
    main()
