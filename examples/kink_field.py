"""What a bias field does to the kink states and the barriers between them.

The mixed cell carries a uniform applied field as the packer does (polyfind.mixedcell, "An
applied field"): the energy is the electric enthalpy ``U - mu . E`` at the macroscopic field
``E``, and its derivative with respect to the field is minus the cell's dipole, static plus
induced.  So to first order in the field the energy of a state or a saddle moves by
``-(mu - mu_sheet) . E`` against the sheet, and a barrier by ``-(mu_saddle - mu_start) . E``:
the dipole difference is the move's activation dipole.

``dipoles`` evaluates the cell dipole of every state a row holds (the sheet, its end states, its
settled saddles or band tops, and the states and bands of the states file) at zero field and
writes them with their differences from the sheet, resolved along the sheet's own polar axis,
across it in the plane, and along the chains.

``relaxed`` is the check of the first order: with the field set along (or against) the sheet's
polar axis it relaxes the sheet and the named states again and settles the named saddles again,
from their zero-field states, and compares the barriers with the dipoles' prediction.  At a
fixed cell, as the bands were run.

``table`` writes the moves' activation dipoles and their barrier shifts at a bias.

Units: dipoles in e.A (1 e.A = 4.803 D), fields in V/um (1 V/um = 1e-4 V/A), energies in
kcal/mol (1 e.A in 1 V/A is 23.06 kcal/mol; in 55 V/um it is 0.127).
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
V_UM = 1e-4  # V/A per V/um
EV_TO_KCAL = 23.0605
KT = 0.596  # kcal/mol at 300 K
SEGMENTS = ("A-I", "I-B")


def _km():
    spec = importlib.util.spec_from_file_location("kink_migration", ROOT / "examples" / "kink_migration.py")
    km = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(km)
    return km


def load(path) -> dict:
    return json.loads(Path(path).read_text())


def saddle_of(band: dict) -> tuple:
    """``(state, energy, how)`` of a band record: its settled saddle if the search converged and
    was accepted, else the band's top image."""
    pol = band.get("polished", {})
    i = int(band["saddle_image"])
    top = float(band["energies"][i])
    # A band converged in the chains' metric has its top on the saddle to 0.01; a search from such a top that ends well below it has slid to
    # another saddle beside an end state.  (A band in per-variable units overstates its top by up to 0.1: there the settled saddle is the number.)
    slid = bool(band.get("band_variables")) and bool(band.get("converged")) and top > float(pol.get("saddle", top)) + 0.05
    if pol.get("converged") and pol.get("a_saddle_between_the_two", True) and "state" in pol and not slid:
        return np.asarray(pol["state"], dtype=float), float(pol["saddle"]), "settled"
    return np.asarray(band["images"][i], dtype=float), top, "band top" + ("" if band.get("converged") else ", band not converged")


def held(row: dict, key: str, states: dict) -> dict:
    """``{name: (state, stored energy against the sheet, kind)}`` of everything the row holds."""
    out = {"A": (np.asarray(row["kinked_crystal"]["state"], dtype=float), 0.0, "the sheet")}
    for name in ("I", "B"):
        if "state" in row.get(name, {}):
            out[name] = (np.asarray(row[name]["state"], dtype=float), float(row[name]["energy"]), "state")
    for seg in SEGMENTS:
        if "images" in row.get(seg, {}):
            x, e, how = saddle_of(row[seg])
            out[f"{seg} saddle"] = (x, e, f"saddle ({how})")
    for name, rec in states.get("settled", {}).get(key, {}).items():
        out[name] = (np.asarray(rec["state"], dtype=float), float(rec["energy"]), "state")
    for name, band in states.get("bands", {}).get(key, {}).items():
        x, e, how = saddle_of(band)
        out[f"{name} saddle"] = (x, e, f"saddle ({how})")
    return out


def build(src: dict, row: dict, preset_applied, km):
    """``(model, sheet state, free, closed)`` of a row, inside an applied preset."""
    DK = km.DK
    row_name, kink = row["row"], int(row["kink_offset_bonds"])
    kkey = f"{row_name}@{kink}{'~' if row.get('kink_mirrored') else ''}"
    letters, _ = DK.sequence(row_name)
    polymer = DK.C.make_polymer(row_name, letters)
    nx, ny = row["supercell"]
    model, _, _ = DK.make_model(polymer, src["kinked_references"][kkey], nx, ny)
    movers = [int(m) for m in str(row["mover"]).split(",")] if str(row["mover"]) != "all" else list(range(model.nc))
    free, closed = km.free_sets(model, movers, row["host"])
    return model, polymer, kkey, np.unique(free), closed


def axes(mu_sheet: np.ndarray) -> np.ndarray:
    """Rows: the sheet's polar axis (its dipole's direction in the plane), the in-plane normal to
    it, and the chain axis."""
    u = np.array([mu_sheet[0], mu_sheet[1], 0.0])
    u /= np.linalg.norm(u)
    return np.array([u, [-u[1], u[0], 0.0], [0.0, 0.0, 1.0]])


def dipoles(args) -> None:
    from polyfind.fitting import PRESET_PARAMETERS

    km = _km()
    src, states = load(args.source), (load(args.states) if args.states and Path(args.states).exists() else {})
    out = Path(args.out)
    disk = load(out) if out.exists() else {"what": __doc__.split("\n\n")[0], "preset": args.preset, "rows": {}}
    for key in args.row:
        row = src["rows"][key]
        rec = disk["rows"].setdefault(key, {})
        rec.update(row=row["row"], supercell=row["supercell"], chains=row["chains"], monomers_per_chain=row["monomers_per_chain"])
        got = rec.setdefault("states", {})
        with PRESET_PARAMETERS[args.preset].applied():
            model, _, _, _, _ = build(src, row, None, km)
            model.pk.set_field(None)
            things = held(row, key, states)
            eA = None
            for name, (x, e_stored, kind) in things.items():
                if args.only and name not in args.only and name != "A":
                    continue
                if name in got and not args.again and np.allclose(got[name].get("stored_energy", 1e9), e_stored, atol=1e-9):
                    continue
                if x.size != model.n_var:
                    print(f"   {name}: {x.size} numbers, the cell has {model.n_var}; left out", flush=True)
                    continue
                t0 = time.time()
                mu = model.dipole(x)
                e = model.physical_energy(x) if (args.energies or name == "A") else None
                if name == "A":
                    eA = e
                got[name] = {"kind": kind, "stored_energy": e_stored, "static": mu["static"].tolist(), "induced": mu["induced"].tolist(),
                             "total": mu["total"].tolist()}
                if e is not None and eA is not None:
                    got[name]["energy"] = e - eA
                print(f"   {key} | {name:14s} {kind:38s} mu = {np.round(mu['total'], 4)} e.A  ({time.time() - t0:.0f} s)", flush=True)
                out.write_text(json.dumps(disk, indent=1))
        mu_a = np.asarray(got["A"]["total"])
        per_chain = mu_a / row["chains"]
        refs = load(args.all_trans).get("all_trans_references", {}) if args.all_trans and Path(args.all_trans).exists() else {}
        if row["row"] in refs:  # the polar axis is the all-trans crystal's; the sheet's own dipole is turned from it by its kinks
            if "all_trans" not in rec or args.again:
                with PRESET_PARAMETERS[args.preset].applied():
                    DK = km.DK
                    polymer = DK.C.make_polymer(row["row"], DK.sequence(row["row"])[0])
                    tt, x_tt, _ = DK.make_model(polymer, refs[row["row"]], 1, 1)
                    mu_tt = tt.dipole(x_tt)
                rec["all_trans"] = {"static": mu_tt["static"].tolist(), "induced": mu_tt["induced"].tolist(), "total": mu_tt["total"].tolist(), "chains": 2}
            mu_tt = np.asarray(rec["all_trans"]["total"]) / 2.0
            turned = 1.0 if per_chain[:2] @ mu_tt[:2] > 0 else -1.0  # the two cells were polished apart: the same crystal, possibly turned half round about the chains
            mu_tt = turned * mu_tt
            ax = axes(mu_tt)
            rec["polar_axis_from"] = "the all-trans crystal's dipole in the lattice's frame" + ("" if turned > 0 else " (that cell turned half round about the chain axis)")
            rec["all_trans_dipole_per_monomer"] = float(np.linalg.norm(mu_tt[:2]) / row["monomers_per_chain"])
            rec["all_trans_induced_fraction"] = float(turned * (np.asarray(rec["all_trans"]["induced"]) / 2.0) @ ax[0] / (mu_tt @ ax[0]))
            d = per_chain - mu_tt
            rec["kink_pair_in_the_sheet"] = {"dipole_against_all_trans_per_kink_pair": {"along_the_polar_axis": float(d @ ax[0]), "across_it": float(d @ ax[1])},
                                             "note": "the sheet's dipole per chain less the all-trans crystal's; the component along the chains depends on where a repeat is cut and is left out"}
        else:
            ax = axes(mu_a)
            rec["polar_axis_from"] = "the sheet's own dipole"
        rec["polar_axis"] = ax[0].tolist()
        rec["sheet_dipole_per_monomer"] = float(per_chain @ ax[0] / row["monomers_per_chain"])
        rec["sheet_dipole_turned_from_the_axis_deg"] = float(np.degrees(np.arctan2(per_chain @ ax[1], per_chain @ ax[0])))
        for name, g in got.items():
            d = np.asarray(g["total"]) - mu_a
            g["against_the_sheet"] = {"along_the_polar_axis": float(d @ ax[0]), "across_it": float(d @ ax[1]), "along_the_chains": float(d @ ax[2]),
                                      "induced_part_along_the_polar_axis": float((np.asarray(g["induced"]) - np.asarray(got["A"]["induced"])) @ ax[0])}
        out.write_text(json.dumps(disk, indent=1))
        print(f"{key}: {len(got)} dipoles; polar axis {np.round(ax[0], 4)} ({rec['polar_axis_from']}); the sheet's dipole along it {rec['sheet_dipole_per_monomer']:.4f} e.A per monomer, "
              f"turned {rec['sheet_dipole_turned_from_the_axis_deg']:+.1f} deg from it" + (f"; all-trans {rec['all_trans_dipole_per_monomer']:.4f} per monomer, a kink pair in the sheet "
              f"{rec['kink_pair_in_the_sheet']['dipole_against_all_trans_per_kink_pair']['along_the_polar_axis']:+.3f} along and "
              f"{rec['kink_pair_in_the_sheet']['dipole_against_all_trans_per_kink_pair']['across_it']:+.3f} across" if "all_trans" in rec else ""), flush=True)
    print(f"written {out}")


def relaxed(args) -> None:
    from polyfind import mixedcell as MC
    from polyfind.fitting import PRESET_PARAMETERS

    km = _km()
    DK = km.DK
    src, states = load(args.source), (load(args.states) if args.states and Path(args.states).exists() else {})
    out = Path(args.out)
    disk = load(out)
    key = args.row
    row, rec = src["rows"][key], disk["rows"][key]
    u = np.asarray(rec["polar_axis"])
    log = lambda s: print(s, flush=True)  # noqa: E731
    with PRESET_PARAMETERS[args.preset].applied():
        model, polymer, kkey, free, closed = build(src, row, None, km)
        things = held(row, key, states)
        sheet = things["A"][0]
        small, xs_, _ = DK.make_model(polymer, src["kinked_references"][kkey], 2, 1)
        blocks = MC.chain_blocks(small, xs_, range(2))
        metric = MC.chain_metric(model, sheet, free, closed, blocks)
        kw = dict(free=free, closed=closed, metric=metric, gtol=args.gtol, maxiter=args.maxiter, patience=args.patience, patience_tol=1e-6)
        for f_vum in args.field:
            field = f_vum * V_UM * u
            tag = f"{f_vum:+g} V/um"
            res = rec.setdefault("relaxed_under_the_field", {}).setdefault(tag, {})
            model.pk.set_field(None)
            e0A = model.physical_energy(sheet)
            model.pk.set_field(field)
            t0 = time.time()
            xA, info = model.relax_on_surface(sheet, **kw)
            eA = model.physical_energy(xA)
            res["A"] = {"enthalpy_change": eA - e0A, "first_order": -float(np.asarray(rec["states"]["A"]["total"]) @ field) * EV_TO_KCAL,
                        "gradient_max": info["gradient_max"], "evaluations": info["evaluations"]}
            log(f"{tag}: the sheet {eA - e0A:+.5f} (first order {res['A']['first_order']:+.5f}); gradient {info['gradient_max']:.1e}; {time.time() - t0:.0f} s")
            for name in args.states_named:
                x0, e_stored, _ = things[name]
                first = e_stored - float(np.asarray(rec["states"][name]["against_the_sheet"]["along_the_polar_axis"])) * f_vum * V_UM * EV_TO_KCAL
                if " saddle" in name:
                    band = (row if name[:-7] in SEGMENTS else states["bands"][key])[name[:-7]]
                    imgs, i = np.asarray(band["images"], dtype=float), int(band["saddle_image"])
                    tangent = imgs[min(i + 1, len(imgs) - 1)] - imgs[max(i - 1, 0)]
                    turning = np.concatenate([model.index("tors"), model.index("phi")])
                    tangent[turning] = km.wrap(tangent[turning])
                    x, info = MC.saddle_polish(model, x0, tangent, free, closed, tol=args.tol, on_surface=True, metric=metric, gtol=args.gtol,
                                               patience=args.patience, patience_tol=1e-6)
                    e = info["energy"] - eA
                    res[name] = {"energy": e, "zero_field": e_stored, "first_order": first, "force_along": info["force"], "gradient_across": info["gradient_max"],
                                 "curvature_along": info["slope"], "converged": bool(info["converged"]), "evaluations": info["evaluations"]}
                else:
                    x, info = model.relax_on_surface(x0, **kw)
                    e = model.physical_energy(x) - eA
                    res[name] = {"energy": e, "zero_field": e_stored, "first_order": first, "gradient_max": info["gradient_max"], "evaluations": info["evaluations"]}
                res[name]["beyond_first_order"] = e - first
                log(f"{tag}: {name:12s} {e:+.5f}  zero field {e_stored:+.5f}  first order {first:+.5f}  beyond it {e - first:+.5f}  ({time.time() - t0:.0f} s)")
                out.write_text(json.dumps(disk, indent=1))
        model.pk.set_field(None)
    print(f"written {out}")


# The moves read from a row's dipoles: (label, start, saddle, end).  A saddle that is not held is skipped.
MOVES = (
    ("the sheet to the wide kink (I)", "A", "A-I saddle", "I"),
    ("I to the completed step (B)", "I", "I-B saddle", "B"),
    ("I to M", "I", "I-M saddle", "M"),
    ("M to B", "M", "M-B saddle", "B"),
    ("the sheet to N", "A", "A-N saddle", "N"),
    ("N to L", "N", "N-L saddle", "L"),
    ("the sheet straight to L", "A", "A-L saddle", "L"),
    ("N to N2", "N", "N-N2 saddle", "N2"),
    ("N: a neighbour joins (chain 4)", "N", "N-NN4 saddle", "NN4"),
    ("N: a neighbour joins (chain 12)", "N", "N-NN12 saddle", "NN12"),
)


def table(args) -> None:
    disk = load(args.source)
    f = args.bias * V_UM * EV_TO_KCAL  # kcal/mol per e.A along the polar axis
    lines = [f"Bias {args.bias:g} V/um along the sheet's polar axis: 1 e.A of dipole is {f:.4f} kcal/mol; kT at 300 K is {KT}.", ""]
    summary = {}
    for key, rec in disk["rows"].items():
        st = rec["states"]
        par = lambda n: st[n]["against_the_sheet"]["along_the_polar_axis"]  # noqa: E731
        lines += [f"### {key}", "",
                  f"{rec['chains']} chains of {rec['monomers_per_chain']} monomers; polar axis from {rec['polar_axis_from']}; the sheet's dipole along it is "
                  f"{rec['sheet_dipole_per_monomer']:.3f} e.A per monomer" + (f" (all-trans {rec['all_trans_dipole_per_monomer']:.3f}; a kink pair in the sheet "
                  f"{rec['kink_pair_in_the_sheet']['dipole_against_all_trans_per_kink_pair']['along_the_polar_axis']:+.2f} e.A along the axis)" if "all_trans" in rec else "") + ".", "",
                  "| state | energy, zero field | dipole against the sheet, along the polar axis (e.A) | across it | along the chains | shift at the bias |",
                  "|---|---|---|---|---|---|"]
        for name, g in st.items():
            if " saddle" in name or name == "A" or abs(g["stored_energy"]) < 1e-3:  # a state that settled back on the sheet is the sheet
                continue
            a = g["against_the_sheet"]
            lines.append(f"| {name} | {g['stored_energy']:+.2f} | {a['along_the_polar_axis']:+.3f} | {a['across_it']:+.3f} | {a['along_the_chains']:+.3f} | {-a['along_the_polar_axis'] * f:+.3f} |")
        lines += ["", "| move | barrier, zero field | activation dipole along the polar axis (e.A) | barrier at the bias | rate factor | back: barrier | activation dipole | rate factor |",
                  "|---|---|---|---|---|---|---|---|"]
        moves = {}
        for label, a, s, b in MOVES:
            if not all(n in st for n in (a, s, b)):
                continue
            ea, es, eb = (st[n]["stored_energy"] for n in (a, s, b))
            da, db = par(s) - par(a), par(s) - par(b)
            moves[label] = {"from": a, "saddle": s, "to": b, "barrier": es - ea, "activation_dipole": da, "barrier_shift_at_bias": -da * f,
                            "rate_factor_at_bias": float(np.exp(da * f / KT)), "barrier_back": es - eb, "activation_dipole_back": db,
                            "rate_factor_back_at_bias": float(np.exp(db * f / KT)), "saddle_kind": st[s]["kind"]}
            m = moves[label]
            lines.append(f"| {label} | {m['barrier']:.2f} | {da:+.3f} | {m['barrier'] - da * f:.2f} | {m['rate_factor_at_bias']:.3f} | {m['barrier_back']:.2f} | {db:+.3f} | {m['rate_factor_back_at_bias']:.3f} |")
        summary[key] = moves
        rel = rec.get("relaxed_under_the_field")
        if rel:
            lines += ["", "Relaxed and settled again under the field (energies against the sheet relaxed under the same field):", "",
                      "| field | state or saddle | zero field | first order | relaxed | beyond first order |", "|---|---|---|---|---|---|"]
            for tag, r in rel.items():
                for name, g in r.items():
                    if name == "A":
                        continue
                    lines.append(f"| {tag} | {name} | {g['zero_field']:+.4f} | {g['first_order']:+.4f} | {g['energy']:+.4f} | {g['beyond_first_order']:+.4f} |")
        lines.append("")
    disk["bias_V_um"] = args.bias
    disk["moves"] = summary
    Path(args.source).write_text(json.dumps(disk, indent=1))
    Path(args.out).write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"written {args.out}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    d = sub.add_parser("dipoles", help="the cell dipole of every state and saddle a row holds")
    d.add_argument("--source", default="docs/predictions/kink_migration_belt4.json")
    d.add_argument("--states", default="docs/predictions/kink_states_belt4.json")
    d.add_argument("--row", action="append", required=True)
    d.add_argument("--only", nargs="*", default=None, help="names of the states to do (the sheet is always done)")
    d.add_argument("--again", action="store_true")
    d.add_argument("--energies", action="store_true", help="also evaluate each state's energy (a check against the stored one)")
    d.add_argument("--all-trans", default="docs/predictions/dilute_kink_belt4.json", help="the file holding the all-trans reference cells: their dipole is the polar axis")
    d.add_argument("--preset", default="pvdf-dft-valence-belt4")
    d.add_argument("--out", default="docs/predictions/kink_field_belt4.json")
    d.set_defaults(func=dipoles)
    r = sub.add_parser("relaxed", help="states relaxed and saddles settled again under the field: the check of the first order")
    r.add_argument("--source", default="docs/predictions/kink_migration_belt4.json")
    r.add_argument("--states", default="docs/predictions/kink_states_belt4.json")
    r.add_argument("--row", required=True)
    r.add_argument("--field", type=float, nargs="+", default=[55.0, -55.0], help="V/um along the sheet's polar axis")
    r.add_argument("--states-named", nargs="+", default=["I", "B", "A-I saddle", "I-B saddle"])
    r.add_argument("--gtol", type=float, default=1e-4)
    r.add_argument("--tol", type=float, default=0.003)
    r.add_argument("--maxiter", type=int, default=600)
    r.add_argument("--patience", type=int, default=40)
    r.add_argument("--preset", default="pvdf-dft-valence-belt4")
    r.add_argument("--out", default="docs/predictions/kink_field_belt4.json")
    r.set_defaults(func=relaxed)
    t = sub.add_parser("table", help="the moves' activation dipoles and barrier shifts at a bias")
    t.add_argument("--source", default="docs/predictions/kink_field_belt4.json")
    t.add_argument("--bias", type=float, default=55.0)
    t.add_argument("--out", default="docs/predictions/kink_field_belt4.md")
    t.set_defaults(func=table)
    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
