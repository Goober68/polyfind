"""Map the candidate screen's margins, dipole gains and switch strains onto the loss model.

Each candidate is the literature terpolymer channel (polyfind.loss.literature_terpolymer: the
Vogel-Fulcher relaxation of the 63/30/7 P(VDF-TrFE-CFE) film, which no chemistry here changes)
with its margin, monomer dipole gain and switch strain replaced by the screen's numbers, evaluated
at the program's operating point (34 and 66 kHz, 55 V/um bias, 1.8 V/um AC, Y 0.2 GPa, mechanical
tan delta 0.05; docs/ELECTROSTRICTION_RESEARCH.md 11-12).  What moves between candidates is the
polar fraction at bias, the strain per field and the heat; the relaxation time is the same for
all, so the loss per unit strain at the drive frequency is NOT a chemistry result here -- it is the
bench's VF fit.  The threshold field E* = dE / dmu is the equilibrium number (zero kelvin step).
"""
import argparse
import json

import numpy as np

from polyfind import loss as L

V_UM = 1e-4  # V/um in V/A


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--screen", default="docs/predictions/candidate_screen_belt4.json")
    ap.add_argument("--e2", default="docs/predictions/e2_polyfind_belt4.json", help="the PVDF and 70/30 rows come from the E2 set")
    ap.add_argument("--bias", type=float, default=55.0, help="V/um")
    ap.add_argument("--ac", type=float, default=1.8, help="V/um")
    ap.add_argument("--young", type=float, default=0.2e9)
    ap.add_argument("--out", default="docs/predictions/candidate_loss_map_belt4.json")
    args = ap.parse_args()
    screen = json.load(open(args.screen))
    rows = {}
    # the two E2 rows: PVDF and the 70/30 model, in the same quantities
    e2 = json.load(open(args.e2))["systems"]
    for key, label in (("pvdf", "pvdf"), ("vdftrfe", "vdf7-trfe3")):
        c = e2[key]["conformations"]
        tt, tg = c["TTTT"], c["TGTG'"]
        from polyfind.electrostriction import dipole_per_monomer
        mu = lambda r: dipole_per_monomer(float(np.linalg.norm(r["P_C_m2"])), r["cell"]["volume_A3"], int(round(r["monomers_per_cell"])))
        rows[label] = {"validated_chemistry": True, "note": e2[key]["label"],
                       "margin": tt["energy_per_monomer_kcal"] - tg["energy_per_monomer_kcal"],
                       "margin_rigid": tt.get("energy_rigid_pendants_per_monomer_kcal", tt["energy_per_monomer_kcal"]) - tg.get("energy_rigid_pendants_per_monomer_kcal", tg["energy_per_monomer_kcal"]),
                       "mu_TT": mu(tt), "dmu_diff": mu(tt) - mu(tg), "axial": tt["c_per_monomer_A"] / tg["c_per_monomer_A"] - 1.0,
                       "areal": (tt["volume_per_monomer_A3"] / tt["c_per_monomer_A"]) / (tg["volume_per_monomer_A3"] / tg["c_per_monomer_A"]) - 1.0,
                       "P_TT": float(np.linalg.norm(tt["P_C_m2"])), "P_TG": float(np.linalg.norm(tg["P_C_m2"])),
                       "tilt_TG": tg.get("pendant_tilt_relaxation_per_monomer_kcal", 0.0) or 0.0}
    for name, e in screen["sequences"].items():
        d = e.get("derived") or {}
        if not d:
            continue
        c = e["conformations"]
        rows[name] = {"validated_chemistry": e["validated_chemistry"], "note": e["note"], "margin": d["margin_TT_minus_TG_kcal_per_monomer"],
                      "margin_rigid": d["margin_rigid_kcal_per_monomer"], "mu_TT": c["TTTT"]["dipole_per_monomer_eA"], "dmu_diff": d["dipole_gain_eA_per_monomer"],
                      "axial": d["axial_transition_strain"], "areal": d["areal_transition_strain"], "P_TT": d["P_TT_C_m2"], "P_TG": d["P_TG_C_m2"],
                      "tilt_TG": c["TGTG'"].get("pendant_tilt_relaxation_per_monomer_kcal", 0.0) or 0.0}
    base = L.literature_terpolymer()
    out = {"operating_point": {"bias_V_um": args.bias, "ac_V_um": args.ac, "young_Pa": args.young, "frequencies_hz": [34e3, 66e3]},
           "baseline": {"dmu_mon": base.dmu_mon, "de_mon": base.de_mon, "n": base.n, "fraction": base.fraction, "f_relax_hz_300K": base.f_relax(300.0)},
           "rows": {}}
    print(f"baseline channel: dmu {base.dmu_mon:.3f} e.A, margin {base.de_mon:.3f} kcal/mol per monomer, n {base.n}, fraction {base.fraction:.2f}, f_relax {base.f_relax(300.0) / 1e3:.0f} kHz")
    print(f"{'candidate':24s} {'chem':5s} {'margin':>7s} {'dmu':>6s} {'E*':>6s} {'axial':>6s} | {'x(bias)':>7s} {'S_ac 34k':>9s} {'tan d 34k':>9s} {'heat 34k':>8s} {'eff 34k':>7s} | {'tan d 66k':>9s} {'S_ac 66k':>9s}")
    def threshold(margin, dmu):
        """E* in MV/m: zero when TTTT is already the ground conformation."""
        return 0.0 if margin <= 0 else margin / L.SwitchingChannel.KCAL_PER_EV / dmu * 1e4

    md = ["| sequence | chemistry | margin TT-TG (rigid) kcal/mol per monomer | TGTG' tilt relaxation | axial, areal switch strain | P(TTTT) C/m^2 | mu(TTTT) e.A | E* (mu_TT) MV/m | E* (mu_TT - mu_TG) MV/m |",
          "|---|---|---|---|---|---|---|---|---|"]
    for name, r in rows.items():
        dmu = r["mu_TT"]  # the gain along the field: the TGTG' state's dipole has no preferred direction in a relaxor
        est = threshold(r["margin"], dmu)
        est_lo = threshold(r["margin"], r["dmu_diff"]) if r["dmu_diff"] > 0 else float("inf")
        md.append(f"| {name} | {'VDF/TrFE (validated)' if r['validated_chemistry'] else 'chlorine (extrapolated)'} | {r['margin']:+.3f} ({r['margin_rigid']:+.3f}) | {r['tilt_TG']:+.3f} | "
                  f"{r['axial'] * 100:+.1f}%, {r['areal'] * 100:+.1f}% | {r['P_TT']:.3f} | {r['mu_TT']:.3f} | {est:.0f} | {est_lo:.0f} |")
        ch = base.scaled(dmu_factor=dmu / base.dmu_mon, de_mon=max(r["margin"], 0.0), label=name)
        ch = L.replace(ch, s_switch=r["axial"])
        ops = {}
        for f in (34e3, 66e3):
            o = ch.operating_point(f, args.bias * V_UM, args.ac * V_UM, args.young, 300.0, 0.05)
            ops[f"{int(f / 1e3)}kHz"] = o
        r2 = dict(r, threshold_MV_m=est, threshold_lower_dmu_MV_m=est_lo, channel={"dmu_mon": ch.dmu_mon, "de_mon": ch.de_mon, "n": ch.n, "s_switch": ch.s_switch}, operating=ops)
        out["rows"][name] = r2
        o34, o66 = ops["34kHz"], ops["66kHz"]
        print(f"{name:24s} {'yes' if r['validated_chemistry'] else 'Cl':5s} {r['margin']:+7.3f} {dmu:6.3f} {est:6.0f} {r['axial'] * 100:+5.1f}% | {o34['polar_fraction']:7.3f} {o34['strain_ac']:9.2e} {o34['tan_delta']:9.3f} {o34['heat_w_per_m3'] / 1e6:8.1f} {o34['efficiency']:7.3f} | {o66['tan_delta']:9.3f} {o66['strain_ac']:9.2e}")
        out["rows"][name]["loss_row_md"] = (f"| {name} | {r['margin']:+.3f} | {est:.0f} | {o34['polar_fraction']:.2f} | {o34['strain_ac'] * 100:.3f}% | {o34['tan_delta']:.3f} / {o66['tan_delta']:.3f} | "
                                            f"{o34['heat_w_per_m3'] / 1e6:.1f} / {o66['heat_w_per_m3'] / 1e6:.1f} | {o34['efficiency']:.2f} / {o66['efficiency']:.2f} | {o34['stroke_rate_per_lambda']:.0f} / {o66['stroke_rate_per_lambda']:.0f} |")
    md_loss = ["| sequence | margin kcal/mol per monomer | E* MV/m | polar fraction at 55 V/um | S_ac at 1.8 V/um, 34 kHz | tan delta 34 / 66 kHz | heat MW/m^3 34 / 66 kHz | efficiency 34 / 66 kHz | stroke per lever 34 / 66 kHz |",
               "|---|---|---|---|---|---|---|---|---|"] + [out["rows"][n]["loss_row_md"] for n in out["rows"]]
    out["tables_md"] = {"screen": "\n".join(md), "loss": "\n".join(md_loss)}
    with open(args.out, "w") as fh:
        json.dump(out, fh, indent=1)
    with open(args.out.replace(".json", ".md"), "w") as fh:
        fh.write("Screen\n\n" + "\n".join(md) + "\n\nLoss model\n\n" + "\n".join(md_loss) + "\n")
    print(f"written {args.out} and the tables beside it")


if __name__ == "__main__":
    main()
