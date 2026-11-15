"""Pin strength per unit and what a pin density does to the switching curve, on the polished margins.

Three tables for docs/CANDIDATE_MATERIALS.md sections 4 and 5, from
docs/predictions/polished_margins_belt4.json (examples/polish_margins.py) and the loss model's
channel (polyfind.loss.literature_terpolymer, as examples/candidate_loss_map.py uses it):

* **pins**: for each sequence that is a matrix with one unit (or two) replaced, the cost of the
  unit in the all-trans and in the helical lattice (the cell's energy against the matrix's, per
  unit) and their difference, the pin strength: how much the unit prefers the helix.  Only the
  difference is a physical quantity (the two costs each carry the unit's own valence energy).
* **uniform**: a matrix of margin ``m0`` with one pin of strength ``s`` every ``N`` monomers has
  the margin ``m0 + s / N`` (the law the fourteen- and twenty-monomer cells follow, 2.1); where
  that puts the channel at the operating point.
* **units** and **three axes**: what a unit costs (its pin) against what it buys (the kink sheet
  made cheaper), and where a unit at a given spacing sits on strain, tan delta and frequency
  together, the frequency on the rate model in which the kink sheets' population sets it,
  anchored on the 64/29/7 cell at the measured film's 364 kHz.  The pin and the cut are both
  taken to fall as one over the spacing (measured for CFE at fourteen and twenty).
* **sites**: where a cell's kink sites have been scanned (examples/kink_site_scan.py), the rate
  goes by all of them: the sum of ``exp(-E / kT)`` over the sites of the repeat and both hands
  of the kink pair, against the same sum for the 64/29/7 cell, in place of one site's energy.
* **statistical**: the same law read locally, each pin-free run of ``L`` monomers at the margin
  ``m0 + s / (L + 1)``, averaged over the runs of a random sequence with pin fraction ``p``
  (run lengths geometric, each monomer weighted by its run).  A hypothesis, not a result: the law
  was measured on periodic cells whose chains are in register.
"""
import argparse
import json

import numpy as np

from polyfind import loss as L

V_UM = 1e-4

# sequence -> (its matrix, monomers in the cell, units replaced, what the unit is)
PINS = {
    "vdf9-trfe4-cfe1-14": ("vdf10-trfe4-14", 14, 1, "CFE, one per fourteen in 29 mol% TrFE"),
    "vdf9-trfe4-ctfe1-14": ("vdf10-trfe4-14", 14, 1, "CTFE, one per fourteen in 29 mol% TrFE"),
    "vdf9-trfe4-an1-14": ("vdf10-trfe4-14", 14, 1, "AN, one per fourteen in 29 mol% TrFE"),
    "vdf9-trfe4-cdfe1-14": ("vdf10-trfe4-14", 14, 1, "CDFE, one per fourteen in 29 mol% TrFE"),
    "vdf8-trfe4-cfe1-cdfe1-14": ("vdf10-trfe4-14", 14, 2, "CFE and CDFE, one of each per fourteen (mean per unit)"),
    # the unit screen of 2026-11-14: the tested groups in another order
    "vdf9-trfe4-tfe1-14": ("vdf10-trfe4-14", 14, 1, "TFE (CF2-CF2), one per fourteen in 29 mol% TrFE"),
    "vdf9-trfe4-rvdf1-14": ("vdf10-trfe4-14", 14, 1, "reversed VDF (CF2-CH2), one per fourteen in 29 mol% TrFE"),
    "vdf9-trfe4-vf1-14": ("vdf10-trfe4-14", 14, 1, "VF (CH2-CHF), one per fourteen in 29 mol% TrFE"),
    "vdf9-trfe4-e1-14": ("vdf10-trfe4-14", 14, 1, "ethylene (CH2-CH2), one per fourteen in 29 mol% TrFE"),
    "vdf9-trfe4-dfe1-14": ("vdf10-trfe4-14", 14, 1, "DFE (CHF-CHF), one per fourteen in 29 mol% TrFE"),
    "vdf7-trfe2-cfe1": ("vdf8-trfe2", 10, 1, "CFE, one per ten beside 20 mol% TrFE"),
    "vdf7-trfe2-cfe1-adj": ("vdf8-trfe2", 10, 1, "CFE next to a TrFE, one per ten"),
    "vdf7-trfe2-ctfe1": ("vdf8-trfe2", 10, 1, "CTFE, one per ten beside 20 mol% TrFE"),
    "vdf7-trfe2-cdfe1": ("vdf8-trfe2", 10, 1, "CDFE, one per ten beside 20 mol% TrFE"),
    "vdf7-trfe2-an1": ("vdf8-trfe2", 10, 1, "AN, one per ten beside 20 mol% TrFE"),
}
# two units in one cell, seven monomers apart: (the cell, its first unit's cell, its second unit's cell, what it is)
COMBINATIONS = {
    "vdf8-trfe4-an1-tfe1-14": ("vdf9-trfe4-an1-14", "vdf9-trfe4-tfe1-14", "AN with TFE"),
    "vdf8-trfe4-an1-e1-14": ("vdf9-trfe4-an1-14", "vdf9-trfe4-e1-14", "AN with ethylene"),
    "vdf8-trfe4-ctfe1-rvdf1-14": ("vdf9-trfe4-ctfe1-14", "vdf9-trfe4-rvdf1-14", "CTFE with reversed VDF"),
    "vdf8-trfe4-an1-rvdf1-14": ("vdf9-trfe4-an1-14", "vdf9-trfe4-rvdf1-14", "AN with reversed VDF"),
    "vdf8-trfe4-cfe1-rvdf1-14": ("vdf9-trfe4-cfe1-14", "vdf9-trfe4-rvdf1-14", "CFE with reversed VDF"),
    "vdf8-trfe4-cfe1-cdfe1-14": ("vdf9-trfe4-cfe1-14", "vdf9-trfe4-cdfe1-14", "CFE with CDFE"),
}
# a matrix without TrFE: (the cell, reversed units per fourteen)
REGIO_MATRICES = {"vdf12-rvdf2-14": 2}
# designed cells on another matrix: given as computed, with no sum of single units beside them
OTHER_CELLS = {"vdf9-trfe3-an1-rvdf1-14": "AN with reversed VDF, one of each per fourteen, on three TrFE in fourteen (21 mol%)"}
# a pin in PVDF: the matrix's margin is the E2 set's (polishing changes nothing in PVDF's cells)
PINS_IN_PVDF = {"vdf9-cfe1": (10, 1, "CFE, one per ten VDF"), "vdf8-cfe2": (10, 2, "CFE, two per ten VDF (per unit)"), "vdf9-an1": (10, 1, "AN, one per ten VDF")}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--margins", default="docs/predictions/polished_margins_belt4.json")
    ap.add_argument("--loss", default="docs/predictions/candidate_loss_map_belt4.json", help="the loss map: the matrix's dipole, and PVDF's margin")
    ap.add_argument("--matrix", default="vdf10-trfe4-14")
    ap.add_argument("--sheets", default="docs/predictions/kink_formation_belt4_polished.json", help="the kink sheets' energies like for like (examples/kink_summary.py)")
    ap.add_argument("--bias", type=float, default=55.0)
    ap.add_argument("--ac", type=float, default=1.8)
    ap.add_argument("--heat-budget", type=float, default=4e6, help="W/m^3")
    ap.add_argument("--ac-max", type=float, default=10.0, help="V/um: the largest AC field the small-signal response is used at, at the heat budget")
    ap.add_argument("--sites", default="docs/predictions/kink_site_scan_belt4.json", help="the kink-site scans (examples/kink_site_scan.py)")
    ap.add_argument("--steps", default="docs/predictions/kink_migration_belt4.json", help="the single-chain moves in sixteen chains (rows '<row>@0 4x2 mover 0 ahead host flexible')")
    ap.add_argument("--out", default="docs/predictions/pin_design_belt4.json")
    args = ap.parse_args()
    rows = json.load(open(args.margins))["rows"]
    loss = json.load(open(args.loss))["rows"]

    def cells(name):
        c, tg = rows[name]["cells"], rows[name]["tg_conformer"]
        return c["TTTT"]["energy_per_monomer"], c[tg]["energy_per_monomer"]

    out = {"margins": args.margins, "operating_point": {"bias_V_um": args.bias, "ac_V_um": args.ac}, "pins": {}}
    md = ["| unit and dilution | margin | matrix's margin | all-trans cost per unit | helix cost per unit | pin strength per unit (TT - TG) |", "|---|---|---|---|---|---|"]
    for name, (matrix, n, units, what) in PINS.items():
        if "margin" not in rows.get(name, {}) or "margin" not in rows.get(matrix, {}):
            continue  # not polished (or its polish is still running)
        (t, g), (t0, g0) = cells(name), cells(matrix)
        rec = {"what": what, "matrix": matrix, "margin": rows[name]["margin"], "matrix_margin": rows[matrix]["margin"], "tt_cost_per_unit": n * (t - t0) / units,
               "tg_cost_per_unit": n * (g - g0) / units, "pin_per_unit": n * (rows[name]["margin"] - rows[matrix]["margin"]) / units, "tg_conformer": rows[name]["tg_conformer"]}
        out["pins"][name] = rec
        md.append(f"| {what} | {rec['margin']:+.3f} | {rec['matrix_margin']:+.3f} | {rec['tt_cost_per_unit']:+.1f} | {rec['tg_cost_per_unit']:+.1f} | **{rec['pin_per_unit']:+.1f}** |")
    m_pvdf = loss["pvdf"]["margin"]
    for name, (n, units, what) in PINS_IN_PVDF.items():
        if "margin" not in rows.get(name, {}):
            continue
        rec = {"what": what, "matrix": "pvdf", "margin": rows[name]["margin"], "matrix_margin": m_pvdf, "pin_per_unit": n * (rows[name]["margin"] - m_pvdf) / units,
               "tg_conformer": rows[name]["tg_conformer"]}
        out["pins"][name] = rec
        md.append(f"| {what} | {rec['margin']:+.3f} | {m_pvdf:+.3f} | - | - | **{rec['pin_per_unit']:+.1f}** |")

    # what a unit costs (the pin) against what it buys (the kink sheet made cheaper), at one per fourteen in the matrix
    try:
        sheets = json.load(open(args.sheets))["rows"]
    except FileNotFoundError:
        sheets = {}
    md_k = ["| unit, one per fourteen | pin per unit | kink sheet, kcal/mol per kink pair | cut against the matrix | cut per kcal/mol of pin | decades of frequency per 0.1 kcal/mol of margin | sheets' population against the CFE cell's | all-trans dipole per monomer, e.A | axial switch strain |",
            "|---|---|---|---|---|---|---|---|---|"]
    out["units"] = {}
    if args.matrix in sheets:
        e0 = sheets[args.matrix]["formation_kcal_per_kink_pair"]
        e_cfe = sheets.get("vdf9-trfe4-cfe1-14", {}).get("formation_kcal_per_kink_pair")
        md_k.append(f"| none: the matrix ({args.matrix}) | 0 | {e0:.2f} | 0 | - | - | {np.exp(-(e0 - e_cfe) / 0.596):.3f} | {loss[args.matrix]['mu_TT']:.3f} | {loss[args.matrix]['axial'] * 100:+.1f}% |")
        for name, rec in sorted(out["pins"].items(), key=lambda kv: kv[1]["pin_per_unit"]):
            if rec["matrix"] != args.matrix or name not in sheets or "one of each" in rec["what"]:
                continue
            e = sheets[name]["formation_kcal_per_kink_pair"]
            pin, cut = rec["pin_per_unit"], e0 - e
            # the exchange rate between the two: the pin moves the margin by pin / N and the cut moves ln(rate) by (14 cut / N) / kT
            pays = pin > 0.05 and cut > 0
            u = {"pin_per_unit": pin, "kink_sheet": e, "cut": cut, "cut_per_pin": cut / pin if pays else None,
                 "decades_per_tenth_of_margin": 14.0 * cut / pin * 0.1 / (0.596 * np.log(10.0)) if pays else None,
                 "population_against_cfe": float(np.exp(-(e - e_cfe) / 0.596)), "mu_TT": loss[name]["mu_TT"], "axial": loss[name]["axial"], "what": rec["what"].split(",")[0]}
            out["units"][name] = u
            ratio = f"**{u['cut_per_pin']:.2f}** | {u['decades_per_tenth_of_margin']:.1f}" if pays else ("no cost: it lowers the margin | -" if cut > 0 else "it stiffens the lattice | -")
            md_k.append(f"| {u['what']} | {pin:+.1f} | {e:.2f} | {cut:+.2f} | {ratio} | {u['population_against_cfe']:.2g} | {u['mu_TT']:.3f} | {u['axial'] * 100:+.1f}% |")

    base = L.literature_terpolymer()
    mat = loss[args.matrix]
    m0 = rows[args.matrix]["margin"]

    # strain, tan delta and frequency together: a unit at one per N, the rate from the kink sheets' population
    md_3 = [f"| unit, one per N | margin | relaxation frequency | AC strain at {args.ac:g} V/um, 34 / 66 kHz | tan delta 34 / 66 kHz | AC strain at the heat budget, 34 / 66 kHz | AC field there, V/um, and what binds |",
            "|---|---|---|---|---|---|---|"]
    out["three_axes"] = {"heat_budget_w_per_m3": args.heat_budget, "ac_max_V_um": args.ac_max,
                         "frequency_is": "the relaxation frequency of the small-signal permittivity (a Cole-Cole fit to a thin film's permittivity sets the anchor): how fast the "
                                         "polarization follows, not a measured strain response and not an actuator's speed",
                         "rate_model": "the kink sheets' population sets the rate: exp(-(E_sheet - E_sheet of the 64/29/7 cell) / kT) times the measured film's 364 kHz",
                         "scaling": "pin and cut as one over the spacing; dipole and switch strain between the matrix's and the unit's cell the same way", "rows": {}}
    if out["units"] and args.matrix in sheets:
        tf_factor = (300.0 - base.tf) / 300.0
        f_anchor = base.f_relax(300.0, 0.0)

        def point(margin, e_sheet, mu, axial):
            shift = ((e_sheet - e_cfe) / L.SwitchingChannel.KCAL_PER_EV) * tf_factor
            ch = L.replace(base.scaled(dmu_factor=mu / base.dmu_mon, de_mon=max(margin, 0.0), ea_shift_ev=shift, label="x"), s_switch=axial)
            return dict(axes(ch), margin=margin, kink_sheet=e_sheet, f_relax_hz=f_anchor * float(np.exp(-(e_sheet - e_cfe) / 0.596)))

        def axes(ch):
            """The three axes of a channel: the point at the program's AC field and the point at the heat budget
            (the field solved from the budget, capped at --ac-max), at both drive frequencies."""
            o = {f: ch.operating_point(f, args.bias * V_UM, args.ac * V_UM, 0.2e9, 300.0, 0.05) for f in (34e3, 66e3)}
            b = {f: ch.budget_point(f, args.bias * V_UM, 0.2e9, 300.0, 0.05, args.heat_budget, args.ac_max * V_UM) for f in (34e3, 66e3)}
            return {"strain_ac": [o[34e3]["strain_ac"], o[66e3]["strain_ac"]], "tan_delta": [o[34e3]["tan_delta"], o[66e3]["tan_delta"]],
                    "efficiency": [o[34e3]["efficiency"], o[66e3]["efficiency"]], "polar_fraction": o[34e3]["polar_fraction"],
                    "at_the_heat_budget": {"strain_ac": [b[34e3]["strain_ac"], b[66e3]["strain_ac"]], "field_ac_V_um": [b[34e3]["field_ac"] / V_UM, b[66e3]["field_ac"] / V_UM],
                                           "binds": [b[34e3]["binds"], b[66e3]["binds"]], "heat_w_per_m3": [b[34e3]["heat_w_per_m3"], b[66e3]["heat_w_per_m3"]]}}

        def line(label, pt):
            fr = pt["f_relax_hz"]
            fr = f"{fr / 1e6:.1f} MHz" if fr >= 1e6 else (f"{fr / 1e3:.0f} kHz" if fr >= 1e3 else f"{fr:.0f} Hz")
            hb = pt["at_the_heat_budget"]
            return (f"| {label} | {pt['margin']:+.3f} | {fr} | {pt['strain_ac'][0] * 100:.3f}% / {pt['strain_ac'][1] * 100:.3f}% | {pt['tan_delta'][0]:.3f} / {pt['tan_delta'][1]:.3f} | "
                    f"**{hb['strain_ac'][0] * 100:.3f}%** / {hb['strain_ac'][1] * 100:.3f}% | {hb['field_ac_V_um'][0]:.1f} ({hb['binds'][0]}) / {hb['field_ac_V_um'][1]:.1f} ({hb['binds'][1]}) |")

        film = dict(axes(base), margin=0.0, f_relax_hz=f_anchor)
        out["three_axes"]["rows"]["calibrated_film"] = film
        md_3.append(line("the measured 63/30/7 film's channel as calibrated", film))

        pt = point(m0, e0, mat["mu_TT"], mat["axial"])
        out["three_axes"]["rows"]["matrix"] = pt
        md_3.append(line("none: the matrix", pt))
        # designed cells: two units together, the margin and the kink sheet computed, beside what the single units would add up to
        md_c = ["| cell | margin: computed, if the pins added | kink sheet: computed, if the cuts added | relaxation frequency | AC strain at "
                f"{args.ac:g} V/um, 34 / 66 kHz | tan delta 34 / 66 kHz | AC strain at the heat budget, 34 / 66 kHz | AC field there, V/um, and what binds |", "|---|---|---|---|---|---|---|---|"]
        out["combinations"] = {}
        for name, (first, second, what) in COMBINATIONS.items():
            if "margin" not in rows.get(name, {}) or first not in out["pins"] or second not in out["pins"]:
                continue
            added_m = m0 + (out["pins"][first]["pin_per_unit"] + out["pins"][second]["pin_per_unit"]) / 14.0
            rec = {"what": what, "margin": rows[name]["margin"], "margin_if_the_pins_added": added_m}
            cells = f"| {what}, one of each per fourteen | {rows[name]['margin']:+.3f}, {added_m:+.3f} | "
            if name in sheets and first in out["units"] and second in out["units"]:
                e = sheets[name]["formation_kcal_per_kink_pair"]
                added_e = e0 - out["units"][first]["cut"] - out["units"][second]["cut"]
                pt = point(rows[name]["margin"], e, loss[name]["mu_TT"], loss[name]["axial"])
                rec.update(kink_sheet=e, kink_sheet_if_the_cuts_added=added_e, three_axes=pt)
                md_c.append(cells + f"{e:.2f}, {added_e:.2f} | " + line("", pt).split("|", 3)[3].lstrip())
            else:
                md_c.append(cells + "not computed | - | - | - | - | - |")
            out["combinations"][name] = rec
        for name, n_rev in REGIO_MATRICES.items():  # PVDF with reversed units and no TrFE, against PVDF itself
            if "margin" not in rows.get(name, {}) or name not in sheets or "vdf10" not in sheets:
                continue
            e = sheets[name]["formation_kcal_per_kink_pair"]
            pt = point(rows[name]["margin"], e, loss[name]["mu_TT"], loss[name]["axial"])
            out["combinations"][name] = {"what": f"PVDF with {n_rev} reversed units per fourteen, no TrFE", "margin": rows[name]["margin"], "pvdf_margin": m_pvdf,
                                         "per_reversed_unit": 14.0 * (rows[name]["margin"] - m_pvdf) / n_rev, "kink_sheet": e,
                                         "pvdf_kink_sheet": sheets["vdf10"]["formation_kcal_per_kink_pair"], "mu_TT": loss[name]["mu_TT"], "three_axes": pt}
            md_c.append(f"| PVDF with {n_rev} reversed units per fourteen, no TrFE (PVDF itself: {m_pvdf:+.3f}, {sheets['vdf10']['formation_kcal_per_kink_pair']:.2f}) | {rows[name]['margin']:+.3f}, - | {e:.2f}, - | "
                        + line("", pt).split("|", 3)[3].lstrip())
        for name, what in OTHER_CELLS.items():
            if "margin" not in rows.get(name, {}) or name not in sheets or name not in loss:
                continue
            e = sheets[name]["formation_kcal_per_kink_pair"]
            pt = point(rows[name]["margin"], e, loss[name]["mu_TT"], loss[name]["axial"])
            out["combinations"][name] = {"what": what, "margin": rows[name]["margin"], "kink_sheet": e, "mu_TT": loss[name]["mu_TT"], "three_axes": pt}
            md_c.append(f"| {what} | {rows[name]['margin']:+.3f}, - | {e:.2f}, - | " + line("", pt).split("|", 3)[3].lstrip())
        # the rate by every scanned kink site: sum over the repeat's sites and both hands of exp(-E / kT), against the anchor cell's sum
        try:
            scans = json.load(open(args.sites)).get("rows", {})
        except FileNotFoundError:
            scans = {}

        def site_sum(row):
            hands = [scans[k] for k in (f"{row}@0", f"{row}@0~") if k in scans]
            if not hands:
                return None
            zs, low = [], None
            for h in hands:
                es = [(v["formation"], v["units"], int(o)) for o, v in h["sites"].items() if v.get("converged", True) and v.get("kink_kept", True)]
                zs.append(float(sum(np.exp(-e / 0.596) for e, _, _ in es)) * (len(h["letters"]) * 1.0 / len(es)) / len(h["letters"]))  # per site, so that a half-scanned repeat compares
                best = min(es)
                low = best if low is None or best < low else low
            return {"z_per_site": float(np.mean(zs)), "hands_scanned": len(hands), "sites_scanned": [len(h["sites"]) for h in hands],
                    "cheapest": {"energy": low[0], "units": low[1], "offset_bonds": low[2]}, "start_site": hands[0]["sites"]["0"]["formation"]}

        md_z = ["| cell | kink sheet at the start site | cheapest site scanned (the two units it sits on) | sites' sum against the 64/29/7 cell's | relaxation frequency by the sum | AC strain at "
                f"{args.ac:g} V/um, 34 kHz | tan delta, 34 kHz | AC strain at the heat budget, 34 kHz |", "|---|---|---|---|---|---|---|---|"]
        out["sites"] = {"note": "kcal/mol per kink pair against all-trans, rigid pendants; the sum is over the scanned sites of the repeat and the scanned hands of the kink pair, "
                                "per site; the frequency is the measured film's 364 kHz times the ratio of sums", "rows": {}}
        z_ref = site_sum("vdf9-trfe4-cfe1-14")
        if z_ref is not None:
            cells_z = [("vdf9-trfe4-cfe1-14", "64/29/7 (one CFE per fourteen): the anchor", rows["vdf9-trfe4-cfe1-14"]["margin"]),
                       (args.matrix, "the 71/29 matrix", m0)]
            cells_z += [(n, out["units"][n]["what"] + ", one per fourteen", rows[n]["margin"]) for n in out["units"] if n != "vdf9-trfe4-cfe1-14"]
            cells_z += [(n, c["what"] + ("" if n in REGIO_MATRICES or n in OTHER_CELLS else ", one of each per fourteen"), c["margin"]) for n, c in out["combinations"].items()]
            for name, label, margin in cells_z:
                z = site_sum(name)
                if z is None or name not in loss:
                    continue
                ratio = z["z_per_site"] / z_ref["z_per_site"]
                pt = point(margin, e_cfe - 0.596 * float(np.log(ratio)), loss[name]["mu_TT"], loss[name]["axial"])
                out["sites"]["rows"][name] = dict(z, what=label, ratio_to_the_anchor=ratio, three_axes=pt)
                hb = pt["at_the_heat_budget"]
                fr = pt["f_relax_hz"]
                fr = f"{fr / 1e6:.1f} MHz" if fr >= 1e6 else (f"{fr / 1e3:.0f} kHz" if fr >= 1e3 else f"{fr:.0f} Hz")
                md_z.append(f"| {label} | {z['start_site']:.2f} | {z['cheapest']['energy']:.2f} ({z['cheapest']['units']}) | {ratio:.3g} | {fr} | {pt['strain_ac'][0] * 100:.3f}% | {pt['tan_delta'][0]:.3f} | "
                            f"**{hb['strain_ac'][0] * 100:.3f}%** |")
        # the other rate model: one chain's cheapest move out of the sheet in sixteen chains (the step to B, or straight to N),
        # against the same move in the anchor cell
        try:
            steps = json.load(open(args.steps)).get("rows", {})
        except FileNotFoundError:
            steps = {}

        def cheapest_move(row):
            e = steps.get(f"{row}@0 4x2 mover 0 ahead host flexible")
            if e is None:
                return None
            step = e.get("barrier_forward_polished", e.get("barrier_forward"))
            moves = {} if step is None else {"the step to B": step}
            if "to_N" in e:
                moves["to N"] = e["to_N"]["barrier_forward"]
            if not moves:
                return None
            name = min(moves, key=moves.get)
            return {"move": name, "barrier": moves[name], "step_to_B": step, "to_N": moves.get("to N"), "start_sheet": e["kinked_crystal"]["formation_per_kink_pair"],
                    "N_back_to_the_sheet": e["to_N"]["barrier_back"] if "to_N" in e else None,
                    "L": e["to_L"]["L"] if "to_L" in e else None, "N_to_L": e["to_L"]["barrier_from_N"] if "to_L" in e else None}

        md_m = ["| cell | the cheapest move, barrier from the sheet | the step to B | rate against the 64/29/7 cell's | relaxation frequency by the move | by the sites' sum | AC strain at "
                f"{args.ac:g} V/um, 34 kHz | tan delta, 34 kHz | AC strain at the heat budget, 34 kHz |", "|---|---|---|---|---|---|---|---|---|"]
        out["moves"] = {"note": "kcal/mol against the sheet at the start site (kink offset 0), sixteen chains, neighbours flexible, cell fixed; the frequency is the measured film's "
                                "364 kHz times exp(-(barrier - the 64/29/7 cell's) / kT)", "rows": {}}
        m_ref = cheapest_move("vdf9-trfe4-cfe1-14")
        if m_ref is not None:
            for name, label, margin in cells_z if z_ref is not None else []:
                mv = cheapest_move(name)
                if mv is None or name not in loss:
                    continue
                ratio = float(np.exp(-(mv["barrier"] - m_ref["barrier"]) / 0.596))
                pt = point(margin, e_cfe - 0.596 * float(np.log(ratio)), loss[name]["mu_TT"], loss[name]["axial"])
                by_sum = out["sites"]["rows"].get(name, {}).get("three_axes", {}).get("f_relax_hz")
                out["moves"]["rows"][name] = dict(mv, what=label, ratio_to_the_anchor=ratio, three_axes=pt, f_relax_by_the_sites_sum_hz=by_sum)
                hb = pt["at_the_heat_budget"]

                def hz(fr):
                    return "-" if fr is None else (f"{fr / 1e6:.1f} MHz" if fr >= 1e6 else (f"{fr / 1e3:.0f} kHz" if fr >= 1e3 else f"{fr:.0f} Hz"))
                to_b = "no path found" if mv["step_to_B"] is None else f"{mv['step_to_B']:.2f}"
                md_m.append(f"| {label} | {mv['move']}, {mv['barrier']:.2f} | {to_b} | {ratio:.3g} | {hz(pt['f_relax_hz'])} | {hz(by_sum)} | {pt['strain_ac'][0] * 100:.3f}% | {pt['tan_delta'][0]:.3f} | "
                            f"**{hb['strain_ac'][0] * 100:.3f}%** |")
        for name, u in out["units"].items():
            for n in (14, 20, 28, 40):
                w = 14.0 / n
                pt = point(m0 + u["pin_per_unit"] / n, e0 - u["cut"] * w, mat["mu_TT"] + (u["mu_TT"] - mat["mu_TT"]) * w, mat["axial"] + (u["axial"] - mat["axial"]) * w)
                out["three_axes"]["rows"].setdefault(name, {})[str(n)] = pt
                md_3.append(line(f"{u['what']}, {n}", pt))

    def op(m):
        ch = L.replace(base.scaled(dmu_factor=mat["mu_TT"] / base.dmu_mon, de_mon=max(m, 0.0), label="x"), s_switch=mat["axial"])
        return ch.operating_point(34e3, args.bias * V_UM, args.ac * V_UM, 0.2e9, 300.0, 0.05)

    pins = {k: out["pins"][n]["pin_per_unit"] for k, n in (("AN", "vdf9-trfe4-an1-14"), ("CTFE", "vdf9-trfe4-ctfe1-14"), ("CFE", "vdf9-trfe4-cfe1-14"), ("CDFE", "vdf9-trfe4-cdfe1-14"))
            if n in out["pins"]}
    spacings = (10, 14, 20, 28, 40, 56, 80, 120)
    uni = {}
    md_u = ["| one pin per | mol% | " + " | ".join(f"{k} ({s:+.1f}): margin, polar fraction, S_ac, tan delta 34 kHz" for k, s in pins.items()) + " |", "|---|---|" + "---|" * len(pins)]
    for n in spacings:
        cellsmd = []
        for k, s in pins.items():
            o = op(m0 + s / n)
            uni.setdefault(k, {})[str(n)] = {"margin": m0 + s / n, "polar_fraction": o["polar_fraction"], "strain_ac": o["strain_ac"], "tan_delta": o["tan_delta"]}
            cellsmd.append(f"{m0 + s / n:+.3f}, {o['polar_fraction']:.2f}, {o['strain_ac'] * 100:.3f}%, {o['tan_delta']:.3f}")
        md_u.append(f"| {n} | {100 / n:.1f} | " + " | ".join(cellsmd) + " |")
    out["uniform"] = {"matrix": args.matrix, "matrix_margin": m0, "pins": pins, "rows": uni}

    stat = {}
    md_s = ["| pin, mol% | uniform spacing: margin, polar fraction, S_ac | statistical sequence (local law): polar fraction, S_ac | monomers in pin-free runs of forty or more |", "|---|---|---|---|"]
    runs = np.arange(1, 2000)
    for k, s in pins.items():
        for p in (0.03, 0.05, 0.07, 0.10):
            w = runs * p * (1.0 - p) ** runs
            w /= w.sum()
            ops = [op(m0 + s / (r + 1)) for r in runs[:400]]  # beyond four hundred monomers a run is the matrix
            tail = op(m0)
            x = float(w[:400] @ np.array([o["polar_fraction"] for o in ops]) + w[400:].sum() * tail["polar_fraction"])
            sa = float(w[:400] @ np.array([o["strain_ac"] for o in ops]) + w[400:].sum() * tail["strain_ac"])
            u = op(m0 + s * p)
            stat.setdefault(k, {})[f"{p:g}"] = {"statistical": {"polar_fraction": x, "strain_ac": sa}, "share_in_runs_of_forty_or_more": float(w[runs >= 40].sum()),
                                                "uniform": {"margin": m0 + s * p, "polar_fraction": u["polar_fraction"], "strain_ac": u["strain_ac"]}}
            md_s.append(f"| {k}, {100 * p:g} | {m0 + s * p:+.3f}, {u['polar_fraction']:.2f}, {u['strain_ac'] * 100:.3f}% | {x:.2f}, {sa * 100:.3f}% | {float(w[runs >= 40].sum()):.2f} |")
    out["statistical"] = {"hypothesis": "each pin-free run of L monomers at the margin m0 + s / (L + 1); runs geometric; every monomer weighted by its run", "rows": stat}
    film = base.operating_point(34e3, args.bias * V_UM, args.ac * V_UM, 0.2e9, 300.0, 0.05)
    out["calibrated_film"] = {"polar_fraction": film["polar_fraction"], "strain_ac": film["strain_ac"], "tan_delta": film["tan_delta"],
                              "note": "the literature channel as calibrated on the 63/30/7 film (margin zero, n 17)"}
    md_c = md_c if out.get("combinations") else []
    md_z = md_z if out.get("sites", {}).get("rows") else []
    md_m = md_m if out.get("moves", {}).get("rows") else []
    out["tables_md"] = {"pins": "\n".join(md), "units": "\n".join(md_k), "three_axes": "\n".join(md_3), "combinations": "\n".join(md_c), "sites": "\n".join(md_z), "moves": "\n".join(md_m),
                        "uniform": "\n".join(md_u), "statistical": "\n".join(md_s)}
    with open(args.out, "w") as fh:
        json.dump(out, fh, indent=1)
    text = ("Pin strength per unit\n\n" + "\n".join(md) + "\n\nWhat a unit costs and what it buys\n\n" + "\n".join(md_k) + "\n\nStrain, tan delta and frequency together (the rate from the kink sheets' population)\n\n" + "\n".join(md_3)
            + "\n\nDesigned cells: two units together\n\n" + "\n".join(md_c) + "\n\nThe rate by every scanned kink site\n\n" + "\n".join(md_z) + "\n\nThe rate by one chain's cheapest move in sixteen chains\n\n" + "\n".join(md_m) + f"\n\nUniform spacing on the {args.matrix} matrix (margin {m0:+.3f})\n\n" + "\n".join(md_u)
            + "\n\nStatistical against uniform, the margin law read locally\n\n" + "\n".join(md_s)
            + f"\n\nThe calibrated film's channel at the operating point: polar fraction {film['polar_fraction']:.2f}, S_ac {film['strain_ac'] * 100:.3f}%, tan delta {film['tan_delta']:.3f}.\n")
    with open(args.out.replace(".json", ".md"), "w") as fh:
        fh.write(text)
    print(text)
    print(f"written {args.out}")


if __name__ == "__main__":
    main()
