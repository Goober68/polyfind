"""The frequency and loss column: candidates through the loss channel with their own relaxation rate.

The relaxation rate of the switching channel is taken as nucleation-limited and scaled between
candidates by the kink-pair formation energy of examples/kink_formation_screen.py:
``rate / rate_ref = exp(-(E_form - E_form_ref) / kT)`` at 300 K, the reference being the CFE 70/20/10
row, whose channel is anchored on the measured Vogel-Fulcher relaxation of the 63/30/7 terpolymer
film (f_relax 364 kHz at 300 K; polyfind.loss.literature_terpolymer).  The migration barrier is
assumed the same for every row (not computed), the margin, dipole gain and switch strain come from
the candidate screen, and the operating point is the program's (55 V/um bias, 1.8 V/um AC, 34 and
66 kHz, Y 0.2 GPa, mechanical tan delta 0.05).  Everything is a model proxy and is labelled so in
docs/CANDIDATE_MATERIALS.md section 8.
"""
import argparse
import json

import numpy as np

from polyfind import loss as L
from polyfind.electrostriction import dipole_per_monomer

KT = 0.596  # kcal/mol at 300 K
V_UM = 1e-4
REF = "vdf7-trfe2-cfe1"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--screen", default="docs/predictions/candidate_screen_belt4.json", help="comma list of screen files; rows merge by name")
    ap.add_argument("--kinks", default="docs/predictions/kink_formation_belt4.json", help="comma list of kink files; rows merge by name")
    ap.add_argument("--ref", default=REF, help="the row anchored on the measured 364 kHz channel (the CFE terpolymer cell)")
    ap.add_argument("--e2", default="docs/predictions/e2_polyfind_belt4.json")
    ap.add_argument("--out", default="docs/predictions/candidate_frequency_map_belt4.json")
    ap.add_argument("--bias", type=float, default=55.0)
    ap.add_argument("--ac", type=float, default=1.8)
    ap.add_argument("--steps", default="", help="kink_migration result file: scale each row's rate by the barrier of the sheet's one-monomer step "
                    "(plus the climb along the sequence from --sites) instead of the kink-pair formation energy; rows are those with a band")
    ap.add_argument("--sites", default="docs/predictions/kink_site_scan_belt4.json", help="kink_site_scan result file (with --steps)")
    ap.add_argument("--margins", default="docs/predictions/polished_margins_belt4.json", help="polished margins; a row found there takes its margin and switch strain from it")
    ap.add_argument("--step-cell", default="1x1", help="which supercell's band to read with --steps")
    ap.add_argument("--step-mover", default="all", help="which mover's band to read with --steps: 'all' (the sheet in register) or a chain index")
    args = ap.parse_args()
    screen = {}
    for path in args.screen.split(","):
        screen.update(json.load(open(path))["sequences"])
    kinks = {}
    for path in args.kinks.split(","):
        kinks.update(json.load(open(path))["rows"])
    e2 = json.load(open(args.e2))["systems"]
    base = L.literature_terpolymer()
    T = 300.0
    rows = {}
    for name, k in kinks.items():
        if "formation_kcal_per_kink_pair" not in k or "@" in name or name.endswith("~"):
            continue
        # statics: from the candidate screen (vdf10 and vdf7-trfe3 are the E2 set's PVDF and 70/30 rows)
        if name in screen and screen[name].get("derived"):
            d = screen[name]["derived"]; tt = screen[name]["conformations"]["TTTT"]
            margin, mu, axial = d["margin_TT_minus_TG_kcal_per_monomer"], tt["dipole_per_monomer_eA"], d["axial_transition_strain"]
        elif name not in ("vdf10", "vdf7-trfe3"):
            continue  # a sequence with a kink energy and no statics row (a site-scan sequence)
        else:
            key = {"vdf10": "pvdf", "vdf7-trfe3": "vdftrfe"}[name]
            c = e2[key]["conformations"]; tt, tg = c["TTTT"], c["TGTG'"]
            margin = tt["energy_per_monomer_kcal"] - tg["energy_per_monomer_kcal"]
            mu = dipole_per_monomer(float(np.linalg.norm(tt["P_C_m2"])), tt["cell"]["volume_A3"], int(round(tt["monomers_per_cell"])))
            axial = tt["c_per_monomer_A"] / tg["c_per_monomer_A"] - 1.0
        rows[name] = {"E_form": k["formation_kcal_per_kink_pair"], "margin": margin, "mu_TT": mu, "axial": axial, "validated_chemistry": k["validated_chemistry"]}
    try:
        polished = json.load(open(args.margins)).get("rows", {})
    except FileNotFoundError:
        polished = {}
    for name, r in rows.items():  # a polished margin supersedes the screen's (docs/CANDIDATE_MATERIALS.md 2.1)
        if name in polished and "margin" in polished[name]:
            r["margin_screen"] = r["margin"]
            r["margin"] = polished[name]["margin"]
            r["axial"] = polished[name].get("axial_transition_strain", r["axial"])
    if args.steps:  # the rate is the sheet's step, not the kink pair's formation
        steps = json.load(open(args.steps))["rows"]
        try:
            sites = json.load(open(args.sites)).get("rows", {})
        except FileNotFoundError:
            sites = {}
        stepped = {}
        for key, e in steps.items():
            if ("barrier_forward" not in e or e.get("host") != "flexible" or f"{e['supercell'][0]}x{e['supercell'][1]}" != args.step_cell
                    or e.get("direction") != "ahead" or str(e.get("mover")) != args.step_mover):
                continue
            name = e["row"]
            if name not in rows:
                continue
            scan = sites.get(f"{name}@{e['kink_offset_bonds']}", {}).get("sites", {})
            # the climb along the matrix: the highest sheet energy on the way between two of its low sites, pins left out (a sheet does not cross them)
            vals = sorted(v["relative"] for v in scan.values() if v.get("converged", True) and v["relative"] < 4.0)
            climb = (max(vals) - min(vals)) if len(vals) > 2 else 0.0
            barrier = e.get("barrier_forward_polished", e["barrier_forward"])  # the settled saddles where both were polished (kink_migration.py --polish), else the bands' tops
            stepped[name] = dict(rows[name], E_form_sheet=rows[name]["E_form"], E_step=barrier, climb=climb, E_form=barrier, step_saddles="polished" if "barrier_forward_polished" in e else "band tops")
        rows = stepped
    if args.ref not in rows:
        raise SystemExit(f"the anchor row {args.ref!r} is not among the rows with a kink energy: {sorted(rows)}")
    e_ref = rows[args.ref]["E_form"]
    tf_factor = (T - base.tf) / T  # an Arrhenius rate change written as a shift of the Vogel-Fulcher activation
    first = "step barrier, kcal/mol" if args.steps else "kink pair, kcal/mol"
    md = [f"| sequence | {first} | rate vs anchor row | f_relax proxy at 300 K | margin | polar fraction at bias | S_ac 34 / 66 kHz | tan delta 34 / 66 kHz | tan delta per % strain, 34 kHz | efficiency 34 / 66 kHz |",
          "|---|---|---|---|---|---|---|---|---|---|"]
    out = {"reference": args.ref, "operating_point": {"bias_V_um": args.bias, "ac_V_um": args.ac}, "rows": {}}
    for name, r in rows.items():
        rate = float(np.exp(-(r["E_form"] - e_ref) / KT))
        ea_shift = ((r["E_form"] - e_ref) / L.SwitchingChannel.KCAL_PER_EV) * tf_factor
        ch = base.scaled(dmu_factor=r["mu_TT"] / base.dmu_mon, de_mon=max(r["margin"], 0.0), ea_shift_ev=ea_shift, label=name)
        ch = L.replace(ch, s_switch=r["axial"])
        f_relax = base.f_relax(T, 0.0) * rate  # the degenerate channel's relaxation frequency, scaled by the nucleation factor
        ops = {f: ch.operating_point(f, args.bias * V_UM, args.ac * V_UM, 0.2e9, T, 0.05) for f in (34e3, 66e3)}
        o34, o66 = ops[34e3], ops[66e3]
        per_strain = o34["tan_delta"] / (100.0 * o34["strain_ac"]) if o34["strain_ac"] > 1e-12 else float("inf")
        out["rows"][name] = dict(r, rate_vs_ref=rate, f_relax_hz=f_relax, operating={"34kHz": o34, "66kHz": o66})
        fr = f"{f_relax / 1e3:.0f} kHz" if f_relax >= 1e3 else f"{f_relax:.0f} Hz"
        md.append(f"| {name} | {r['E_form']:+.2f} | {rate:.3g} | {fr} | {r['margin']:+.3f} | {o34['polar_fraction']:.2f} | {o34['strain_ac'] * 100:.3f}% / {o66['strain_ac'] * 100:.3f}% | "
                  f"{o34['tan_delta']:.3f} / {o66['tan_delta']:.3f} | {per_strain:.2f} | {o34['efficiency']:.2f} / {o66['efficiency']:.2f} |")
        print(md[-1])
    out["table_md"] = "\n".join(md)
    with open(args.out, "w") as fh:
        json.dump(out, fh, indent=1)
    with open(args.out.replace(".json", ".md"), "w") as fh:
        fh.write("\n".join(md) + "\n")
    print(f"written {args.out}")


if __name__ == "__main__":
    main()
