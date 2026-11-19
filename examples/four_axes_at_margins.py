"""The four axes of the cells whose margin another evaluator has, at that margin: strain, tan delta, frequency, blocked force per area.

docs/CANDIDATE_MATERIALS.md 5.6 ranks cells on the four axes with belt4's margins, and DFT has since
said where those are wrong (docs/MLIP_SECOND_OPINION.md 3.10 to 3.15).  This puts another
evaluator's margins through the same tables (examples/pin_design.py) and changes nothing else: the
dipole and the switch strain are belt4's cells', the kink sheets and the single-chain moves are
belt4's energies, the crystal stiffness is belt4's.  So it says what the margin does to the four
axes.  It does not say what the other inputs would be under an evaluator that has been checked:
none of them has a DFT test.

The margins come from ``--margins-from`` (``docs/predictions/rows14_potential.json``: the
potential's own, on the fourteen-monomer cells).  Each is also shifted by every value of
``--shifts``, the evaluator's error bar, so the table shows which conclusions survive it.

The loss model takes a margin below zero as zero (``polyfind.loss``: the channel is then degenerate
at no field, as the measured film's is calibrated).  A cell with all-trans the lower one is a
ferroelectric with a remanent polarization, which the model does not hold; such rows are marked.

    python examples/four_axes_at_margins.py --out docs/predictions/four_axes_potential_margins.json
"""
import argparse
import json
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
PRED = HERE.parent / "docs" / "predictions"
PN_NM2 = 1e6


def axes_of(pt):
    hb = pt["at_the_heat_budget"]
    through = hb.get("blocked_stress_through_the_crystal_pa")
    return {"margin_the_model_used": max(pt["margin"], 0.0), "polar_fraction_at_the_bias": pt["polar_fraction"], "f_relax_hz": pt["f_relax_hz"],
            "tan_delta_34kHz": pt["tan_delta"][0], "strain_ac_at_the_heat_budget_34kHz": hb["strain_ac"][0], "field_ac_V_um": hb["field_ac_V_um"][0], "binds": hb["binds"][0],
            "blocked_pn_nm2_by_the_films_assumed_modulus": hb["blocked_stress_pa"][0] / PN_NM2,
            "blocked_pn_nm2_through_the_crystal": None if through is None else through[0] / PN_NM2}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--margins-from", type=Path, default=PRED / "rows14_potential.json", help="rows -> margin, from the evaluator that replaces the tables' preset")
    ap.add_argument("--base", type=Path, default=PRED / "polished_margins_belt4.json", help="the tables' own margins file: everything but the margin is kept")
    ap.add_argument("--shifts", default="-0.1,0,0.05,0.1", help="kcal/mol per monomer added to every replaced margin, one table each")
    ap.add_argument("--sheet-shift", nargs="*", default=[], help="ROW=DELTA, kcal/mol per kink pair: move every scanned kink site of ROW (both hands) by DELTA, the difference "
                                                                  "between another evaluator's kink sheet and the tables' at one site, taken as common to the row's sites")
    ap.add_argument("--sites", type=Path, default=PRED / "kink_site_scan_belt4.json", help="the tables' kink-site scans (examples/kink_site_scan.py)")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    shifts = {item.split("=")[0]: float(item.split("=")[1]) for item in args.sheet_shift}
    theirs = {k: float(v["margin"]) for k, v in json.loads(args.margins_from.read_text())["rows"].items()}
    base = json.loads(args.base.read_text())
    out = {"what": "the four axes (examples/pin_design.py) with the margins of some rows replaced and everything else as the tables have it: dipole, switch strain, kink "
                   "energies and crystal stiffness are belt4's and have no DFT test. Frequency is the relaxation frequency of the small-signal permittivity; strain is the AC "
                   "strain at the heat budget at 34 kHz; blocked force is per area, pN/nm2 (= MPa)",
           "margins_from": str(args.margins_from.name), "base": str(args.base.name), "replaced": theirs, "belt4_margins": {k: base["rows"][k]["margin"] for k in theirs},
           "shifts": {}}
    if shifts:
        out["sheet_shifts_kcal_per_kink_pair"] = shifts
        out["what"] += ". The kink sites of the rows in sheet_shifts are moved by those amounts (both hands, every site alike): only the rate model that sums over kink sites sees it"
    with tempfile.TemporaryDirectory() as tmp:
        extra = []
        if shifts:
            scans = json.loads(args.sites.read_text())
            for key, row in scans["rows"].items():
                if key.split("@")[0] in shifts and "fixed cell" not in key:
                    for site in row["sites"].values():
                        site["formation"] += shifts[key.split("@")[0]]
            (Path(tmp) / "sites.json").write_text(json.dumps(scans))
            extra = ["--sites", str(Path(tmp) / "sites.json")]
        for shift in [float(x) for x in args.shifts.split(",")]:
            rows = json.loads(json.dumps(base))
            for k, m in theirs.items():
                rows["rows"][k]["margin"] = m + shift
            margins, result = Path(tmp) / "margins.json", Path(tmp) / "out.json"
            margins.write_text(json.dumps(rows))
            subprocess.run([sys.executable, "-W", "ignore", str(HERE / "pin_design.py"), "--margins", str(margins), "--out", str(result), *extra], stdout=subprocess.DEVNULL, check=True)
            res = json.loads(result.read_text())
            entry = {"film_channel_as_calibrated": axes_of(res["three_axes"]["rows"]["calibrated_film"]), "rows": {}}
            for k, m in theirs.items():
                rec = {"margin": m + shift, "all_trans_is_lower": m + shift < 0.0}
                for model, key in (("by_the_sum_over_kink_sites", "sites"),) + (() if shifts else (("by_one_chains_move", "moves"),)):
                    if k in res.get(key, {}).get("rows", {}):
                        rec[model] = axes_of(res[key]["rows"][k]["three_axes"])
                entry["rows"][k] = rec
            out["shifts"][f"{shift:+.2f}"] = entry
    args.out.write_text(json.dumps(out, indent=1))
    for shift, entry in out["shifts"].items():
        film = entry["film_channel_as_calibrated"]
        print(f"\nmargins shifted by {shift}; the measured film's channel: {film['f_relax_hz'] / 1e3:.0f} kHz, tan delta {film['tan_delta_34kHz']:.3f}, "
              f"strain {100 * film['strain_ac_at_the_heat_budget_34kHz']:.3f}%, blocked {film['blocked_pn_nm2_by_the_films_assumed_modulus']:.2f}")
        for k, rec in entry["rows"].items():
            for model in ("by_the_sum_over_kink_sites", "by_one_chains_move"):
                if model not in rec:
                    continue
                a = rec[model]
                through = "-" if a["blocked_pn_nm2_through_the_crystal"] is None else f"{a['blocked_pn_nm2_through_the_crystal']:.1f}"
                print(f"  {k:26s} {rec['margin']:+.3f}{' (below zero: taken as 0)' if rec['all_trans_is_lower'] else '':26s} {model:27s} polar {a['polar_fraction_at_the_bias']:.2f}  "
                      f"{a['f_relax_hz'] / 1e3:8.0f} kHz  tan delta {a['tan_delta_34kHz']:.4f}  strain {100 * a['strain_ac_at_the_heat_budget_34kHz']:.3f}%  "
                      f"blocked {a['blocked_pn_nm2_by_the_films_assumed_modulus']:.2f}; {through}")
    print(f"written {args.out}")


if __name__ == "__main__":
    main()
