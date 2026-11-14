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
    "vdf7-trfe2-cfe1": ("vdf8-trfe2", 10, 1, "CFE, one per ten beside 20 mol% TrFE"),
    "vdf7-trfe2-cfe1-adj": ("vdf8-trfe2", 10, 1, "CFE next to a TrFE, one per ten"),
    "vdf7-trfe2-ctfe1": ("vdf8-trfe2", 10, 1, "CTFE, one per ten beside 20 mol% TrFE"),
    "vdf7-trfe2-cdfe1": ("vdf8-trfe2", 10, 1, "CDFE, one per ten beside 20 mol% TrFE"),
    "vdf7-trfe2-an1": ("vdf8-trfe2", 10, 1, "AN, one per ten beside 20 mol% TrFE"),
}
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
    md_k = ["| unit, one per fourteen | pin per unit | kink sheet, kcal/mol per kink pair | cut against the matrix | cut per kcal/mol of pin | sheets' population against the CFE cell's | all-trans dipole per monomer, e.A | axial switch strain |",
            "|---|---|---|---|---|---|---|---|"]
    out["units"] = {}
    if args.matrix in sheets:
        e0 = sheets[args.matrix]["formation_kcal_per_kink_pair"]
        e_cfe = sheets.get("vdf9-trfe4-cfe1-14", {}).get("formation_kcal_per_kink_pair")
        md_k.append(f"| none: the matrix ({args.matrix}) | 0 | {e0:.2f} | 0 | - | {np.exp(-(e0 - e_cfe) / 0.596):.3f} | {loss[args.matrix]['mu_TT']:.3f} | {loss[args.matrix]['axial'] * 100:+.1f}% |")
        for name, rec in sorted(out["pins"].items(), key=lambda kv: kv[1]["pin_per_unit"]):
            if rec["matrix"] != args.matrix or name not in sheets or "one of each" in rec["what"]:
                continue
            e = sheets[name]["formation_kcal_per_kink_pair"]
            u = {"pin_per_unit": rec["pin_per_unit"], "kink_sheet": e, "cut": e0 - e, "cut_per_pin": (e0 - e) / rec["pin_per_unit"],
                 "population_against_cfe": float(np.exp(-(e - e_cfe) / 0.596)), "mu_TT": loss[name]["mu_TT"], "axial": loss[name]["axial"]}
            out["units"][name] = u
            md_k.append(f"| {rec['what'].split(',')[0]} | {u['pin_per_unit']:+.1f} | {e:.2f} | {u['cut']:.2f} | **{u['cut_per_pin']:.2f}** | {u['population_against_cfe']:.1f} | {u['mu_TT']:.3f} | {u['axial'] * 100:+.1f}% |")

    base = L.literature_terpolymer()
    mat = loss[args.matrix]
    m0 = rows[args.matrix]["margin"]

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
    out["tables_md"] = {"pins": "\n".join(md), "units": "\n".join(md_k), "uniform": "\n".join(md_u), "statistical": "\n".join(md_s)}
    with open(args.out, "w") as fh:
        json.dump(out, fh, indent=1)
    text = ("Pin strength per unit\n\n" + "\n".join(md) + "\n\nWhat a unit costs and what it buys\n\n" + "\n".join(md_k) + f"\n\nUniform spacing on the {args.matrix} matrix (margin {m0:+.3f})\n\n" + "\n".join(md_u)
            + "\n\nStatistical against uniform, the margin law read locally\n\n" + "\n".join(md_s)
            + f"\n\nThe calibrated film's channel at the operating point: polar fraction {film['polar_fraction']:.2f}, S_ac {film['strain_ac'] * 100:.3f}%, tan delta {film['tan_delta']:.3f}.\n")
    with open(args.out.replace(".json", ".md"), "w") as fh:
        fh.write(text)
    print(text)
    print(f"written {args.out}")


if __name__ == "__main__":
    main()
