"""Tables of the lone kink pair, the kink sheet and its one-monomer step, from the result files.

Reads ``docs/predictions/dilute_kink_belt4.json`` (examples/dilute_kink.py),
``kink_migration_belt4.json`` (examples/kink_migration.py) and ``polished_margins_belt4.json``
(examples/polish_margins.py), adds the pendant tilts of each polished all-trans and kinked
two-chain cell where they have not been computed yet (cached in the files), and writes the tables
of docs/CANDIDATE_MATERIALS.md section 8 as markdown: ``docs/predictions/kink_summary_belt4.md``.
"""
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np

from polyfind import mixedcell as MC
from polyfind.fitting import PRESET_PARAMETERS

DS = importlib.util.spec_from_file_location("dk", Path(__file__).with_name("dilute_kink.py"))
DK = importlib.util.module_from_spec(DS)
DS.loader.exec_module(DK)
KT = 0.596  # kcal/mol at 300 K


def locked_update(path: Path, fn):
    import fcntl
    with open(str(path) + ".lock", "w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        disk = json.load(open(path)) if path.exists() else {}
        fn(disk)
        with open(path, "w") as fh:
            json.dump(disk, fh, indent=1)


def tilt_of(polymer, d: dict) -> float:
    """Pendant-tilt relaxation of a polished two-chain cell, kcal/mol per cell."""
    model, x, _ = DK.make_model(polymer, {k: (np.asarray(v, dtype=float) if isinstance(v, list) else v) for k, v in d.items()}, 1, 1)
    return float(MC.tilt_relaxation(model, x)[0])


def step_response(row: dict, name: str) -> dict:
    """How each chain of a step row's cell sits in the state ``name`` (I or B) against the kinked
    crystal: pure geometry from the stored state vectors.  The lateral and axial displacements are
    taken against the median of the chains that do not move their kink pair (one chain of the cell
    is held, so the raw numbers carry that chain's share).  Returns the mover's own record and one
    record per shell of neighbours (chains at the same distance from the mover to 0.2 A)."""
    wrap = lambda a: (np.asarray(a, dtype=float) + 180.0) % 360.0 - 180.0  # noqa: E731
    nc, P = int(row["chains"]), 2 * int(row["monomers_per_chain"])
    xA, x = np.asarray(row["kinked_crystal"]["state"], dtype=float), np.asarray(row[name]["state"], dtype=float)
    A, B, gam = xA[:3]
    lat = np.array([[A, 0.0], [B * np.cos(np.deg2rad(gam)), B * np.sin(np.deg2rad(gam))]])
    pA, px = xA[4:4 + 4 * nc].reshape(nc, 4), x[4:4 + 4 * nc].reshape(nc, 4)
    iA, ix = xA[4 + 4 * nc:].reshape(nc, -1), x[4 + 4 * nc:].reshape(nc, -1)
    m = int(row["mover"])
    d = pA[:, :2] - pA[m, :2]
    dist = np.linalg.norm((d - np.round(d)) @ lat, axis=1)
    others = [s for s in range(nc) if s != m]
    shift = (px[:, :2] - pA[:, :2]) @ lat
    shift -= np.median(shift[others], axis=0)
    axial = px[:, 3] - pA[:, 3]
    axial = axial - np.median(axial[others])
    turn = wrap(px[:, 2] - pA[:, 2])
    dt = np.abs(wrap(ix[:, :P] - iA[:, :P]))
    one = lambda s: {"move_A": float(np.linalg.norm(shift[s])), "turn_deg": float(turn[s]), "axial_A": float(axial[s]), "torsion_deg": float(dt[s].max())}  # noqa: E731
    shells: dict = {}
    for s in sorted(others, key=lambda s: dist[s]):
        key = next((k for k in shells if abs(k - dist[s]) < 0.2), float(dist[s]))
        shells.setdefault(key, []).append(one(s))
    return {"mover": one(m), "mover_gauche_deg": [float(v) for v in wrap(ix[m, :P]) if abs(v) < 120.0],
            "shells": [{"distance_A": k, "chains": len(v), "move_A": max(c["move_A"] for c in v), "turn_deg": max(abs(c["turn_deg"]) for c in v),
                        "axial_A": max((c["axial_A"] for c in v), key=abs), "torsion_deg": max(c["torsion_deg"] for c in v)} for k, v in shells.items()]}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--dilute", default="docs/predictions/dilute_kink_belt4.json")
    ap.add_argument("--migration", default="docs/predictions/kink_migration_belt4.json")
    ap.add_argument("--margins", default="docs/predictions/polished_margins_belt4.json")
    ap.add_argument("--sites", default="docs/predictions/kink_site_scan_belt4.json")
    ap.add_argument("--slide", default="docs/predictions/chain_slide_belt4.json")
    ap.add_argument("--preset", default="pvdf-dft-valence-belt4")
    ap.add_argument("--no-tilts", action="store_true", help="do not compute missing pendant tilts")
    ap.add_argument("--out", default="docs/predictions/kink_summary_belt4.md")
    ap.add_argument("--sheets", default="docs/predictions/kink_formation_belt4_polished.json")
    args = ap.parse_args()
    dil_p, mig_p, mar_p = Path(args.dilute), Path(args.migration), Path(args.margins)
    dil = json.load(open(dil_p)) if dil_p.exists() else {}
    mig = json.load(open(mig_p)) if mig_p.exists() else {}
    mar = json.load(open(mar_p)) if mar_p.exists() else {}
    sit = json.load(open(args.sites)) if Path(args.sites).exists() else {}
    md = []

    # --- pendant tilts of the polished cells the sheet energies need -------------------------------
    if not args.no_tilts:
        with PRESET_PARAMETERS[args.preset].applied():
            for kkey, d in mig.get("kinked_references", {}).items():
                row = kkey.split("@")[0]
                letters, _ = DK.sequence(row)
                polymer = None
                if "tilt_per_cell" not in d and "state" in d:
                    polymer = DK.C.make_polymer(row, letters)
                    d["tilt_per_cell"] = tilt_of(polymer, d)
                    locked_update(mig_p, lambda disk, kkey=kkey, v=d["tilt_per_cell"]: disk["kinked_references"][kkey].__setitem__("tilt_per_cell", v))
                    print(f"tilts of the kinked cell {kkey}: {d['tilt_per_cell']:+.4f} kcal/mol per cell", flush=True)
                tt = dil.get("all_trans_references", {}).get(row)
                if tt is not None and "tilt_per_cell" not in tt and "state" in tt:
                    polymer = polymer or DK.C.make_polymer(row, letters)
                    tt["tilt_per_cell"] = tilt_of(polymer, tt)
                    locked_update(dil_p, lambda disk, row=row, v=tt["tilt_per_cell"]: disk["all_trans_references"][row].__setitem__("tilt_per_cell", v))
                    print(f"tilts of the all-trans cell {row}: {tt['tilt_per_cell']:+.4f} kcal/mol per cell", flush=True)

    # --- the lone kink pair -------------------------------------------------------------------------
    md += ["**Lone kink pair in the all-trans crystal** (kcal/mol per kink pair; one chain of the supercell kinked, cell held at the all-trans values)", "",
           "| sequence | chains x monomers | axial registry: barrier, at 1 A (per monomer, rigid slide) | neighbours fixed | neighbours rigid bodies | neighbours flexible |", "|---|---|---|---|---|---|"]
    for key, e in sorted(dil.get("rows", {}).items(), key=lambda kv: (kv[1].get("row", ""), kv[1].get("monomers_per_chain", 0), kv[0])):
        st = e.get("stages", {})
        if not st or len(e.get("kinked_chains", [])) != 1:
            continue
        sl = e.get("slide")
        reg = f"{sl['barrier_per_monomer']:.2f}, {sl['at_1A_per_monomer']:.2f}" if sl else "-"

        def cell(name):
            r = st.get(name)
            if r is None:
                return "-"
            if r.get("collapsed"):
                return f"reverts to all-trans ({r['formation']:+.2f})"
            g = np.round(r["kinks"][0]["gauche_deg"], 0)
            capped = f"; {r['angles_at_window']} angles at the window" if r.get("angles_at_window") else ""
            if r.get("descent") and not r["descent"]["converged"]:
                capped += f"; a damped descent stopped at its budget, NOT a minimum: gradient {r['descent']['largest_projected_gradient']:.2f}"
            return f"{r['formation']:+.2f} (gauche bonds {g[0]:.0f}, {g[1]:.0f}{capped})"

        extra = ""
        if abs(e.get("angle_stiffness", 233.037) - 233.037) > 0.5:
            extra += f", angle restraint {e['angle_stiffness']:g} (chain stiffness matched to DFT)"
        if abs(e.get("angle_window_deg", 8.0) - 8.0) > 1e-9:
            extra += f", angle window {e['angle_window_deg']:g} deg"
        md.append(f"| {key.split(' ')[0]}{extra} | {e['chains']} x {e['monomers_per_chain']} | {reg} | {cell('rigid_host')} | {cell('host_rigid_bodies')} | {cell('flexible_host')} |")

    held = [(key, e) for key, e in sorted(dil.get("rows", {}).items(), key=lambda kv: (kv[1].get("row", ""), kv[1].get("chains", 0), kv[0])) if e.get("torsion_scan")]
    if held:
        md += ["", "Energy against the kink's own torsion (the two gauche bonds held at +tau and -tau degrees, the rest of the stage relaxed from the stage's stored state; "
               "kcal/mol above all-trans; examples/kink_torsion_scan.py). A rise before 120 degrees is the rim that holds the kink; a break in the curve is another branch.", ""]
        for key, e in held:
            for stage in [st for st in ("rigid_host", "host_rigid_bodies", "flexible_host") if st in e["torsion_scan"]]:
                pts = sorted(e["torsion_scan"][stage].items(), key=lambda kv: float(kv[0]))
                stored = e.get("stages", {}).get(stage, {})
                note = f"stored state {stored['formation']:+.2f}" + (", a descent's last point and not a minimum" if stored.get("descent") and not stored["descent"]["converged"] else "") if stored else ""
                md.append(f"* `{key}`, {stage.replace('_', ' ')} ({e['chains']} chains; {note}): " + "; ".join(f"{float(k):g}: {v['energy']:.2f}" for k, v in pts))

    # --- one chain slid, the lattice responding ----------------------------------------------------
    sld = json.load(open(args.slide)) if Path(args.slide).exists() else {}
    rows = {k: r for k, r in sld.get("rows", {}).items() if r.get("chains", 0) >= 4}
    if rows:
        def val(r, name):
            if name not in r:
                return "-"
            moving = abs(r[name].get("axial_force_on_slid_chain", 0.0)) > 0.01 * r.get("repeat_monomers", 1)
            return f"{r[name]['barrier']:.3f}" + (" (not stationary along the slide)" if moving else "")

        def moves(r):
            sh = list(r.get("shift", {}).get("shells_A", {}).items())[:2]
            return "; ".join(f"{float(k):.2f} A: {v['largest_move_A']:.3f} A, {v['largest_turn_deg']:.1f} deg" for k, v in sh) or "-"

        def order(kv):
            name, cell = kv[0].split(" ")[:2]
            return (not name.startswith("vdf"), len(kv[0].split(" ")) > 2, len(name), name, kv[1]["chains"], cell)

        md += ["", "**One chain slid half a monomer along its axis, the lattice responding** (kcal/mol per monomer of the slid chain against the unslid lattice; "
               "cell fixed unless stated; examples/chain_slide.py)", "",
               "| repeat and cell | chains | rigid lattice | neighbours turn | turn and move aside | and flex | with pendant tilts | turn and move aside, lateral cell free "
               "| the two nearest shells at 'turn and move aside': largest move, turn |", "|---|---|---|---|---|---|---|---|---|"]
        for key, r in sorted(rows.items(), key=order):
            if key.split(" ")[0] not in ("vdf1", "vdf2", "vdf10"):
                continue
            free = r.get("shift_cell_free")
            cf = f"{free['barrier']:.3f} (a {free['cell_change']['a']:+.2%}, b {free['cell_change']['b']:+.2%})" if free else "-"
            tl = r.get("flexible_with_tilts")
            md.append(f"| {key} | {r['chains']} | {r['rigid']:.3f} | {val(r, 'turn')} | {val(r, 'shift')} | {val(r, 'flexible')} | {tl['barrier']:.3f} | {cf} | {moves(r)} |" if tl else
                      f"| {key} | {r['chains']} | {r['rigid']:.3f} | {val(r, 'turn')} | {val(r, 'shift')} | {val(r, 'flexible')} | - | {cf} | {moves(r)} |")
        fam = {r["chains"]: r for k, r in rows.items() if k.split(" ")[0] == "vdf1" and len(k.split(" ")) == 2 and k.split(" ")[1] in ("2x1", "4x2", "6x3", "8x4")}

        def limit(name):
            ser = {n: r[name]["barrier"] for n, r in fam.items() if name in r}
            if len(ser) < 2:
                return None
            n1, n2 = sorted(ser)[-2:]
            slope = (ser[n1] - ser[n2]) / (1 / n1 - 1 / n2)
            return ser[n2] - slope / n2, n1, n2

        lf, lc = limit("shift"), limit("shift_cell_free")
        if lf and lc:
            md += ["", f"Dilute limit, linear in one over the number of chains from the two largest cells of the 2:1 family: fixed cell {lf[0]:.3f} (from {lf[1]} and {lf[2]} chains), "
                       f"lateral cell free {lc[0]:.3f} (from {lc[1]} and {lc[2]}): **{0.5 * (lf[0] + lc[0]):.2f}**, {0.5 * (lf[0] + lc[0]) / fam[max(fam)]['rigid']:.0%} of the rigid lattice's "
                       f"{fam[max(fam)]['rigid']:.2f}. The fixed cell bounds it from above and the free cell from below at every size."]
        prof = [(k, r, name) for k, r in sorted(rows.items(), key=order) for name in r if name.startswith("profile_")]
        if prof:
            md += ["", "Along the slide (slip in monomer periods: rigid lattice, relaxed):", ""]
            for key, r, name in prof:
                pts = sorted(r[name].items(), key=lambda kv: float(kv[0]))
                top = max(pts, key=lambda kv: kv[1]["relaxed"])
                rtop = max(pts, key=lambda kv: kv[1]["rigid"])
                hyst = max(abs(e.get("up", e["relaxed"]) - e.get("down", e["relaxed"])) for _, e in pts)
                md.append(f"* `{key}`, {name[8:].replace('_', ' ')}: " + "; ".join(f"{float(k):.3g}: {e['rigid']:.2f}, {e['relaxed']:.2f}" for k, e in pts)
                          + f". Highest relaxed {top[1]['relaxed']:.2f} at {float(top[0]):.3g} (rigid {rtop[1]['rigid']:.2f} at {float(rtop[0]):.3g}); the two sweeps differ by {hyst:.3f} at most.")

    # --- the sheet ------------------------------------------------------------------------------------
    md += ["", "**Kink sheet: every chain kinked in register, own cell** (kcal/mol per kink pair against the all-trans crystal)", "",
           "| sequence | the search's cells (line-group all-trans, soft closure) | both cells polished, rigid pendants | with pendant tilts | kinked cell a, b, gamma, c | all-trans cell |", "|---|---|---|---|---|---|"]
    sheets = {}
    for kkey, d in sorted(mig.get("kinked_references", {}).items()):
        row = kkey.split("@")[0]
        tt = dil.get("all_trans_references", {}).get(row)
        if tt is None or "state" not in d or "state" not in tt:
            continue
        search = (d["search_energy"] - tt["search_energy"]) / 2
        pol = (d["energy"] - tt["energy"]) / 2
        tilt = f"{pol + (d['tilt_per_cell'] - tt['tilt_per_cell']) / 2:+.2f}" if "tilt_per_cell" in d and "tilt_per_cell" in tt else "-"
        if tilt != "-" and kkey.endswith("@0"):
            letters = DK.sequence(row)[0]
            sheets[row] = {"letters": letters, "validated_chemistry": set(letters) <= {"V", "R"}, "kink_offset_bonds": 0,
                           "formation_kcal_per_kink_pair": pol + (d["tilt_per_cell"] - tt["tilt_per_cell"]) / 2, "formation_rigid_kcal_per_kink_pair": pol,
                           "formation_search_rigid_kcal_per_kink_pair": search}
        md.append(f"| {kkey} | {search:+.2f} | {pol:+.2f} | {tilt} | {', '.join(f'{v:.3f}' for v in d['state'][:4])} | {', '.join(f'{v:.3f}' for v in tt['state'][:4])} |")

    if sheets:  # the sheet energies in the kink-formation files' format, for examples/candidate_frequency_map.py --kinks
        Path(args.sheets).write_text(json.dumps({"preset": args.preset, "protocol": "the kink sheet (every chain kinked in register, own cell) against the all-trans crystal, "
                                                 "both cells polished in the mixed cell (every torsion free, closure exact), pendant tilts included",
                                                 "energy": "kcal/mol per kink pair", "rows": sheets}, indent=1))

    # --- the one-monomer step -----------------------------------------------------------------------
    md += ["", "**One chain's kink pair one monomer ahead of its neighbours'** (kcal/mol per moving chain against the sheet; cell fixed)", "",
           "| sequence | cell, moving chains | host | wide kink I | saddle A-I | saddle I-B | B (one monomer ahead) | barrier ahead / back | rate at 300 K, 3e12 /s attempt | bands converged |", "|---|---|---|---|---|---|---|---|---|---|"]
    for key, e in sorted(mig.get("rows", {}).items()):
        if "barrier_forward" not in e:
            continue
        nc = e["chains"]
        rate = 3e12 * np.exp(-e["barrier_forward"] * e.get("moving_chains", 1) / KT)  # what moves together must gather its barrier together
        conv = "yes" if e["A-I"]["converged"] and e["I-B"]["converged"] else f"A-I force {e['A-I']['force_max']:.2f}, I-B force {e['I-B']['force_max']:.2f}"
        who = f"all {nc} (the sheet in register)" if str(e.get("mover")) == "all" else f"1 of {nc}"
        md.append(f"| {e['row']}@{e['kink_offset_bonds']} {e['direction']} | {e['supercell'][0]}x{e['supercell'][1]}, {who} | {e['host']} | {e['I']['energy']:+.2f} | {e['A-I']['saddle']:+.2f} | "
                  f"{e['I-B']['saddle']:+.2f} | {e['B']['energy']:+.2f} | {e['barrier_forward']:.2f} / {e['barrier_back']:.2f} | {rate:.2g} /s | {conv} |")

    ends = [(k, e) for k, e in sorted(mig.get("rows", {}).items(), key=lambda kv: (kv[1].get("row", ""), kv[1].get("host", ""), kv[1].get("chains", 0), str(kv[1].get("mover"))))
            if "barrier_forward" not in e and "I" in e and "B" in e and str(e.get("mover")) != "all"]
    if ends:
        md += ["", "End points only (no band)", "",
               "| sequence | cell, moving chains | host | wide kink I | B (one monomer ahead) |", "|---|---|---|---|---|"]
        for key, e in ends:
            who = f"1 of {e['chains']}" if "," not in str(e.get("mover")) else f"chains {e['mover']} of {e['chains']} together, per chain"
            md.append(f"| {e['row']}@{e['kink_offset_bonds']} {e['direction']} | {e['supercell'][0]}x{e['supercell'][1]}, {who} | {e['host']} | {e['I']['energy']:+.2f} | {e['B']['energy']:+.2f} |")

    # --- the single-chain step against the size of the cell ------------------------------------------
    single = {}
    for key, e in mig.get("rows", {}).items():
        if e.get("row") == "vdf10" and str(e.get("mover")) == "0" and e.get("direction") == "ahead" and e.get("kink_offset_bonds") == 0 and not e.get("kink_mirrored") and "B" in e and "I" in e:
            single.setdefault(tuple(e["supercell"]), {})[e["host"]] = e
    if single:
        def band(e):
            if e is None or "A-I" not in e or "I-B" not in e:
                return "- | - | -"
            unc = "" if e["A-I"].get("converged") and e["I-B"].get("converged") else f" (forces {e['A-I']['force_max']:.2f}, {e['I-B']['force_max']:.2f})"
            top = max(0.0, e["A-I"]["saddle"], e["I-B"]["saddle"])
            return f"{e['A-I']['saddle']:.2f} | {e['I-B']['saddle']:.2f} | {top:.2f} / {top - e['B']['energy']:.2f}{unc}"

        md += ["", "**The single-chain step against the size of the cell** (PVDF, chain 0's kink pair one monomer ahead of every other chain's; kcal/mol against the sheet; cell fixed)", "",
               "| cell (a x b tiling of the two-chain cell) | one chain in | neighbours rigid chains: I, B | neighbours flexible: I | saddle A-I | saddle I-B | barrier ahead / back | B |",
               "|---|---|---|---|---|---|---|---|"]
        for cell in sorted(single, key=lambda c: (2 * c[0] * c[1], -c[0])):
            fl, bd = single[cell].get("flexible"), single[cell].get("bodies")
            n = 2 * cell[0] * cell[1]
            note = " (its a-neighbours move too: a stack)" if cell == (1, 1) else ""
            md.append(f"| {cell[0]}x{cell[1]} | {n}{note} | {f'{bd['I']['energy']:+.2f}, {bd['B']['energy']:+.2f}' if bd else '-'} | {f'{fl['I']['energy']:+.2f}' if fl else '-'} | "
                      f"{band(fl)} | {f'{fl['B']['energy']:+.3f}' if fl else '-'} |")
        big = max((c for c in single if "flexible" in single[c]), key=lambda c: c[0] * c[1])
        e = single[big]["flexible"]
        md += ["", f"How the lattice answers in the {big[0]}x{big[1]} cell ({e['chains']} chains, neighbours flexible): displacements against the kinked crystal, "
               "relative to the median of the chains whose kink pair stays", "",
               "| state | chain or shell (distance from the moving chain) | chains | moves aside by, at most | turns, at most | shifts along the axis, at most | largest torsion change |", "|---|---|---|---|---|---|---|"]
        for name, label in (("I", "wide kink I"), ("B", "B, one monomer ahead")):
            r = step_response(e, name)
            g = ", ".join(f"{v:+.0f}" for v in r["mover_gauche_deg"])
            md.append(f"| {label} | the moving chain (gauche bonds at {g} deg) | 1 | {r['mover']['move_A']:.2f} A | {r['mover']['turn_deg']:+.1f} deg | {r['mover']['axial_A']:+.2f} A | - |")
            for sh in r["shells"]:
                md.append(f"| | {sh['distance_A']:.2f} A | {sh['chains']} | {sh['move_A']:.2f} A | {sh['turn_deg']:.1f} deg | {sh['axial_A']:+.2f} A | {sh['torsion_deg']:.1f} deg |")

    # --- the sheet along the sequence -----------------------------------------------------------------
    scans = {k: v for k, v in sit.get("rows", {}).items() if not v.get("fixed_cell")}
    if scans:
        md += ["", "**The kink sheet along the sequence** (kcal/mol per kink pair against the all-trans crystal, rigid pendants; the kink pair's first unit is the monomer named, "
               "gauche bonds on the even bonds; `~` is the mirrored kink pair G- T G+)", ""]
        stats = {}
        for key, e in sorted(scans.items()):
            letters = e["letters"]
            n = len(letters)
            row = key.split("@")[0]
            by = {int(o) // 2: v for o, v in e["sites"].items() if v.get("converged", True)}
            cells = " ".join(f"{letters[m]}{letters[(m + 1) % n]} {by[m]['formation']:.2f}" if m in by else f"{letters[m]}{letters[(m + 1) % n]} -" for m in range(n))
            md.append(f"* `{letters}` {'G- T G+' if e.get('mirrored') else 'G+ T G-'} ({len(by)} of {n} positions): {cells}")
            st = stats.setdefault(row, {"letters": letters, "E": [], "scans": 0})
            st["E"] += [v["formation"] for v in by.values()]
            st["scans"] += 1
        md += ["", "| sequence | scans (handednesses) | lowest site | highest site | sheets per repeat at 300 K, sum over sites of exp(-E/kT) | against 64/29/7 |", "|---|---|---|---|---|---|"]
        ref = stats.get("vdf9-trfe4-cfe1-14")
        zref = float(np.exp(-np.array(ref["E"]) / KT).sum()) / ref["scans"] if ref else None
        for row, st in stats.items():
            E = np.array(st["E"])
            z = float(np.exp(-E / KT).sum()) / st["scans"]  # per handedness scanned, so that rows with one and with two scans compare
            md.append(f"| {row} ({st['letters']}) | {st['scans']} | {E.min():.2f} | {E.max():.2f} | {z:.2e} per handedness | {z / zref:.3f} |" if zref else
                      f"| {row} ({st['letters']}) | {st['scans']} | {E.min():.2f} | {E.max():.2f} | {z:.2e} per handedness | - |")

    # --- the margins ------------------------------------------------------------------------------------
    md += ["", "**Margins with both cells polished** (kcal/mol per monomer, all-trans minus TG; positive when TG is lower)", "",
           "| sequence | the screen's margin | polished | all-trans cell lower by | TG cell lower by | TG conformer | torsions off the ideal pattern, TT / TG | axial transition strain |", "|---|---|---|---|---|---|---|---|"]
    for row, e in mar.get("rows", {}).items():
        if "margin" not in e:
            continue
        tt, tg = e["cells"]["TTTT"], e["cells"][e["tg_conformer"]]
        md.append(f"| {row} ({e['letters']}) | {e['screen']['margin']:+.3f} | {e['margin']:+.3f} | {tt['polish_gain_per_monomer']:.3f} | {tg['polish_gain_per_monomer']:.3f} | {e['tg_conformer']} | "
                  f"{tt['largest_torsion_off_ideal_deg']:.0f} / {tg['largest_torsion_off_ideal_deg']:.0f} deg | {100 * e.get('axial_transition_strain', tt['c_per_monomer'] / tg['c_per_monomer'] - 1):+.1f}% |")
    text = "\n".join(md) + "\n"
    Path(args.out).write_text(text)
    print(text)
    print(f"written {args.out}")


if __name__ == "__main__":
    main()
