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


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--dilute", default="docs/predictions/dilute_kink_belt4.json")
    ap.add_argument("--migration", default="docs/predictions/kink_migration_belt4.json")
    ap.add_argument("--margins", default="docs/predictions/polished_margins_belt4.json")
    ap.add_argument("--preset", default="pvdf-dft-valence-belt4")
    ap.add_argument("--no-tilts", action="store_true", help="do not compute missing pendant tilts")
    ap.add_argument("--out", default="docs/predictions/kink_summary_belt4.md")
    args = ap.parse_args()
    dil_p, mig_p, mar_p = Path(args.dilute), Path(args.migration), Path(args.margins)
    dil = json.load(open(dil_p)) if dil_p.exists() else {}
    mig = json.load(open(mig_p)) if mig_p.exists() else {}
    mar = json.load(open(mar_p)) if mar_p.exists() else {}
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
            return f"{r['formation']:+.2f} (gauche bonds {g[0]:.0f}, {g[1]:.0f}{capped})"

        extra = ""
        if abs(e.get("angle_stiffness", 233.037) - 233.037) > 0.5:
            extra += f", angle restraint {e['angle_stiffness']:g} (chain stiffness matched to DFT)"
        if abs(e.get("angle_window_deg", 8.0) - 8.0) > 1e-9:
            extra += f", angle window {e['angle_window_deg']:g} deg"
        md.append(f"| {key.split(' ')[0]}{extra} | {e['chains']} x {e['monomers_per_chain']} | {reg} | {cell('rigid_host')} | {cell('host_rigid_bodies')} | {cell('flexible_host')} |")

    # --- the sheet ------------------------------------------------------------------------------------
    md += ["", "**Kink sheet: every chain kinked in register, own cell** (kcal/mol per kink pair against the all-trans crystal)", "",
           "| sequence | the search's cells (line-group all-trans, soft closure) | both cells polished, rigid pendants | with pendant tilts | kinked cell a, b, gamma, c | all-trans cell |", "|---|---|---|---|---|---|"]
    for kkey, d in sorted(mig.get("kinked_references", {}).items()):
        row = kkey.split("@")[0]
        tt = dil.get("all_trans_references", {}).get(row)
        if tt is None or "state" not in d or "state" not in tt:
            continue
        search = (d["search_energy"] - tt["search_energy"]) / 2
        pol = (d["energy"] - tt["energy"]) / 2
        tilt = f"{pol + (d['tilt_per_cell'] - tt['tilt_per_cell']) / 2:+.2f}" if "tilt_per_cell" in d and "tilt_per_cell" in tt else "-"
        md.append(f"| {kkey} | {search:+.2f} | {pol:+.2f} | {tilt} | {', '.join(f'{v:.3f}' for v in d['state'][:4])} | {', '.join(f'{v:.3f}' for v in tt['state'][:4])} |")

    # --- the one-monomer step -----------------------------------------------------------------------
    md += ["", "**One chain's kink pair one monomer ahead of its neighbours'** (kcal/mol per moving chain against the sheet; cell fixed)", "",
           "| sequence | cell, moving chains | host | wide kink I | saddle A-I | saddle I-B | B (one monomer ahead) | barrier ahead / back | rate at 300 K, 3e12 /s attempt | bands converged |", "|---|---|---|---|---|---|---|---|---|---|"]
    for key, e in sorted(mig.get("rows", {}).items()):
        if "barrier_forward" not in e:
            continue
        nc = e["chains"]
        rate = 3e12 * np.exp(-e["barrier_forward"] / KT)
        conv = "yes" if e["A-I"]["converged"] and e["I-B"]["converged"] else f"A-I force {e['A-I']['force_max']:.2f}, I-B force {e['I-B']['force_max']:.2f}"
        md.append(f"| {e['row']}@{e['kink_offset_bonds']} {e['direction']} | {e['supercell'][0]}x{e['supercell'][1]}, 1 of {nc} | {e['host']} | {e['I']['energy']:+.2f} | {e['A-I']['saddle']:+.2f} | "
                  f"{e['I-B']['saddle']:+.2f} | {e['B']['energy']:+.2f} | {e['barrier_forward']:.2f} / {e['barrier_back']:.2f} | {rate:.2g} /s | {conv} |")

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
