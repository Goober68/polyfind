"""Render ``docs/E1_E2_PREDICTIONS.md`` from the prediction files (tables from the JSON, never typed by hand).

    python examples/e1_e2_report.py [--notes docs/predictions/e1_e2_notes.md] [--out docs/E1_E2_PREDICTIONS.md]

The notes file (hand-written prose: status, definitions, flags) goes first; the tables follow; the
cross-check against Dipole's existing numbers is appended as it prints.
"""
import argparse
import io
import json
import subprocess
import sys

import numpy as np


def fmt(x, nd=3):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return str(x)
    if not np.isfinite(v):
        return "-"
    if v != 0 and abs(v) < 10 ** (-nd):
        return f"{v:.2e}"
    return f"{v:.{nd}f}"


def e1_tables(pf):
    out = [f"\nPotential: {pf['potential']}. Frame: {pf['frame']}. Strain step {pf['strain']:+.3f} on each side.\n"]
    out.append("\n### E1, cell and stiffness\n")
    out.append("| system | a (A, polar) | b (A, lateral) | gamma | c (A, chain) | V (A^3) | C_11 | C_12 | C_22 | C_66 (GPa) |")
    out.append("|---|---|---|---|---|---|---|---|---|---|")
    for key, s in pf["systems"].items():
        c = s["cell"]
        b = np.asarray(s["C_GPa"]["block"], dtype=float)
        v = s["C_GPa"]["voigt"]
        def C(I, J):
            return fmt(b[v.index(I), v.index(J)], 2) if I in v and J in v else "-"
        out.append(f"| {s['label']} | {fmt(c['a'])} | {fmt(c['b'])} | {fmt(c['gamma_deg'], 2)} | {fmt(c['c'], 4)} | {fmt(c['volume_A3'], 2)} | "
                   f"{C(1, 1)} | {C(1, 2)} | {C(2, 2)} | {C(6, 6)} |")
    out.append("\n### E1, dielectric tensors at the stress-free cell, proper piezoelectric constants, polarization\n")
    out.append("| system | eps_inf (x, y, z) | eps_0 relaxed-ion (x, y, z) | e_11, e_12, e_26 (C/m^2) | d_11, d_12 (pm/V) | P0 (C/m^2) |")
    out.append("|---|---|---|---|---|---|")
    for key, s in pf["systems"].items():
        ei = np.diag(np.asarray(s["eps_inf"], dtype=float))
        e0 = np.diag(np.asarray(s["eps_relaxed_ion"], dtype=float))
        e = np.asarray(s["e_proper_C_m2"]["values"], dtype=float)
        cols = s["e_proper_C_m2"]["cols_voigt"]
        d = np.asarray(s["d_direct_pm_V"], dtype=float)
        def E(i, J):
            return fmt(e[i, cols.index(J)], 4) if J in cols else "-"
        def D(i, J):
            return fmt(d[i, cols.index(J)], 2) if J in cols else "-"
        p0 = np.asarray(s["P0_C_m2"], dtype=float)
        out.append(f"| {s['label']} | {fmt(ei[0], 4)}, {fmt(ei[1], 4)}, {fmt(ei[2], 4)} | {fmt(e0[0], 4)}, {fmt(e0[1], 4)}, {fmt(e0[2], 4)} | "
                   f"{E(0, 1)}, {E(0, 2)}, {E(1, 6)} | {D(0, 1)}, {D(0, 2)} | {fmt(p0[0], 4)}, {fmt(p0[1], 4)}, {fmt(p0[2], 4)} |")
    out.append("\n### E1, the true relaxed-ion dielectric tensor against normal strain (central differences at +/-1%)\n")
    out.append("| system | strain axis | d eps_xx / d eta | d eps_yy / d eta | d eps_zz / d eta | dP_x / d eta (true, C/m^2) | eps diag at +1% | eps diag at -1% |")
    out.append("|---|---|---|---|---|---|---|---|")
    for key, s in pf["systems"].items():
        for axis, r in s["strained"].items():
            de = np.diag(np.asarray(r["d_eps_d_eta"], dtype=float))
            dp = np.asarray(r["dP_true_d_eta_C_m2"], dtype=float)
            plus = np.diag(np.asarray(r["cells"][f"{pf['strain']:+.3f}"]["eps"], dtype=float))
            minus = np.diag(np.asarray(r["cells"][f"{-pf['strain']:+.3f}"]["eps"], dtype=float))
            out.append(f"| {s['label']} | {axis} | {fmt(de[0])} | {fmt(de[1])} | {fmt(de[2])} | {fmt(dp[0], 4)} | "
                       f"{', '.join(fmt(v, 4) for v in plus)} | {', '.join(fmt(v, 4) for v in minus)} |")
    out.append("\n### E1, electrostriction (field along reference axis i, Voigt strain J), both routes\n")
    out.append("| system | field | M_strain (1e-18 m^2/V^2), J = 1, 2 | M_polarization, J = 1, 2 | chi_free | Q (m^4/C^2), J = 1, 2 |")
    out.append("|---|---|---|---|---|---|")
    for key, s in pf["systems"].items():
        if "M_strain_m2_V2" not in s:
            out.append(f"| {s['label']} | - | not computed on this cell (cost) | | | |")
            continue
        ms = np.asarray(s["M_strain_m2_V2"], dtype=float) * 1e18
        mp = np.asarray(s["M_polarization_m2_V2"], dtype=float) * 1e18
        q = np.asarray(s["Q_m4_C2"], dtype=float)
        chi = np.asarray(s["chi_free"], dtype=float)
        chi = np.diag(chi) if chi.ndim == 2 else chi.reshape(-1)  # stored as the diagonal
        for i, ax in enumerate(("x (polar)", "y (lateral)", "z (chain)")):
            out.append(f"| {s['label'] if i == 0 else ''} | {ax} | {fmt(ms[i, 0], 5)}, {fmt(ms[i, 1], 5)} | {fmt(mp[i, 0], 5)}, {fmt(mp[i, 1], 5)} | "
                       f"{fmt(chi[i])} | {fmt(q[i, 0], 2)}, {fmt(q[i, 1], 2)} |")
    return "\n".join(out)


def e2_tables(pf):
    out = [f"\nPotential: {pf['potential']}. Energies per monomer on the same potential, so only differences carry meaning.\n"]
    out.append("\n### E2, conformations\n")
    out.append("| system | conformation | E per monomer (kcal/mol) | margin vs ground | |P| (C/m^2) | c per monomer (A) | conformational strain vs ground | V per monomer (A^3) | cell a, b, gamma, c |")
    out.append("|---|---|---|---|---|---|---|---|---|")
    for key, s in pf["systems"].items():
        for conf, r in s["conformations"].items():
            if "failed" in r:
                out.append(f"| {s['label']} | {conf} | failed to pack: {r['failed'][:60]} | | | | | | |")
                continue
            c = r["cell"]
            out.append(f"| {s['label']} | {conf} | {fmt(r['energy_per_monomer_kcal'], 4)} | {fmt(r.get('margin_vs_ground_kcal', float('nan')), 3)} | "
                       f"{fmt(np.linalg.norm(r['P_C_m2']), 4)} | {fmt(r['c_per_monomer_A'], 4)} | {fmt(100 * r.get('conformational_strain_vs_ground', float('nan')), 1)}% | "
                       f"{fmt(r['volume_per_monomer_A3'], 2)} | {fmt(c['a'])}, {fmt(c['b'])}, {fmt(c['gamma_deg'], 1)}, {fmt(c['c'])} |")
    return "\n".join(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--e1", default="docs/predictions/e1_polyfind.json")
    ap.add_argument("--e2", default="docs/predictions/e2_polyfind.json")
    ap.add_argument("--dipole", default="docs/predictions/e1_dipole_existing.json")
    ap.add_argument("--notes", default="docs/predictions/e1_e2_notes.md")
    ap.add_argument("--out", default="docs/E1_E2_PREDICTIONS.md")
    args = ap.parse_args()
    parts = []
    try:
        with open(args.notes) as fh:
            parts.append(fh.read())
    except FileNotFoundError:
        parts.append("# polyfind's pre-registered predictions for E1 and E2\n")
    with open(args.e1) as fh:
        parts.append("\n## E1 predictions\n" + e1_tables(json.load(fh)))
    try:
        with open(args.e2) as fh:
            parts.append("\n## E2 predictions\n" + e2_tables(json.load(fh)))
    except FileNotFoundError:
        pass
    try:
        res = subprocess.run([sys.executable, "examples/crosscheck_e1_e2.py", "--polyfind", args.e1, "--dipole", args.dipole],
                             capture_output=True, text=True, check=False)
        parts.append("\n## First cross-check against what Dipole already has\n\nOutput of `examples/crosscheck_e1_e2.py` on "
                     f"`{args.dipole}` (Dipole's statuses quoted inline):\n\n```\n{res.stdout.strip()}\n```\n")
    except Exception as exc:  # the report is still useful without it
        parts.append(f"\n(cross-check not run: {exc!r})\n")
    text = "\n".join(parts) + f"\n\n*Rendered by `examples/e1_e2_report.py` from `{args.e1}`, `{args.e2}` and `{args.dipole}`; re-run it when a file changes.*\n"
    with open(args.out, "w") as fh:
        fh.write(text)
    print(f"written {args.out} ({len(text)} chars)")


if __name__ == "__main__":
    main()
