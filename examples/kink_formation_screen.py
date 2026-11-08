"""Tier 3, first number: the kink-pair formation energy in the packed crystal, per sequence.

The screen of examples/candidate_screen.py says whether a chemistry switches at the operating
field and how much it strains; it says nothing about how fast, which is what tan delta at the
drive frequency is.  Both seats agree the switch is collective (Dipole, field_barrier_v1.md:
the field moves a single bond's barrier by 1-4%), so the rate is set by nucleating a kink pair
(T G+ T G- T inside a trans run) in the lattice and moving it.  This screen computes the first
of the two: the energy of a kink pair embedded in the all-trans two-chain crystal, relaxed with
the cell, against the all-trans cell, under the same preset and pendant relaxation as the
candidate screen.  The pair sits in the first two VDF monomers of each ten-monomer sequence
(bonds 0-2), away from the TrFE and the defect, so the rows compare the lattice's cost of a
kink in a VDF run as the chemistry around it changes.  The migration barrier is not here; a
rate proxy f ~ exp(-E_form / kT) relative to the CFE terpolymer row is reported, labelled as
such.  The kinks of the two chains are in register (one per chain per repeat), which is the
model cell's limitation.
"""
import argparse
import importlib.util
import json
import time
from pathlib import Path

import numpy as np

from polyfind.fitting import PRESET_PARAMETERS

CS = importlib.util.spec_from_file_location("cs", Path(__file__).with_name("candidate_screen.py"))
C = importlib.util.module_from_spec(CS)
CS.loader.exec_module(C)
T, GP, GM = 0, 1, 2
KT = 0.596  # kcal/mol at 300 K
KINK = [GP, T, GM] + [T] * 17  # a kink pair in the first two monomers; the chain direction is restored after it
ROWS = ["vdf10", "vdf8-trfe2", "vdf7-trfe3", "vdf7-trfe2-cfe1", "vdf7-trfe2-cfe1-adj", "vdf7-trfe2-ctfe1", "vdf7-trfe2-cdfe1", "vdf7-trfe2-an1", "vdf9-an1"]
EXTRA = {"vdf10": ("VVVVVVVVVV", "PVDF as a ten-monomer sequence"), "vdf7-trfe3": ("VVRVVRVVVR", "the 70/30 model of the E2 set as a ten-monomer sequence")}


def build_from(polymer, seq, label, start_params):
    """E2's build, started from another reference's packing parameters (the same basin)."""
    from polyfind import mechanics as M
    from polyfind.ewald import EwaldSpec
    from polyfind.polarizability import Polarizable
    kw = {"pack_kw": {"table_cache_dir": None}, "lj_cutoff": "force", "coulomb": "ewald", "ewald": EwaldSpec(), "polarizable": Polarizable()}
    refine_kw = {k: kw[k] for k in ("lj_cutoff", "coulomb", "ewald", "polarizable")}
    ref, rr = M.refined_reference(polymer, seq, label=label, refine_kw=refine_kw, start_params=start_params, **kw)
    return M.relax_reference(ref), rr


def one_cell_from(polymer, seq, label, start_params):
    """candidate_screen.one_cell with the build seeded from ``start_params``."""
    import types
    saved = C.E2.build
    C.E2.build = lambda pol, sq, lab: build_from(pol, sq, lab, start_params)
    try:
        return C.one_cell(polymer, seq, label)
    finally:
        C.E2.build = saved


def reference_params(polymer, label):
    """The all-trans cell's refined packing parameters (a, b, gamma, phi1, phi2, dz, flip)."""
    ref, rr = C.E2.build(polymer, [T] * 20, label)
    return np.asarray(ref.params, dtype=float)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--same-basin", action="store_true", help="refine the kinked chain from the all-trans cell's own packing parameters instead of its own packing screen (the kink's formation energy in the trans crystal's basin)")
    ap.add_argument("--preset", default="pvdf-dft-valence-belt4")
    ap.add_argument("--screen", default="docs/predictions/candidate_screen_belt4.json", help="the candidate screen: its TTTT rows are the references")
    ap.add_argument("--rows", default=",".join(ROWS))
    ap.add_argument("--out", default="docs/predictions/kink_formation_belt4.json")
    args = ap.parse_args()
    screen = json.load(open(args.screen))["sequences"]
    out_path = Path(args.out)
    results = json.load(open(out_path)) if out_path.exists() else {}
    results.update({"preset": args.preset, "protocol": "same basin: kinked chain refined from the all-trans cell's packing parameters" if args.same_basin else "own packing screen", "kink": "T G+ T G- T pair at bonds 0-2 of the twenty-bond repeat, one per chain, relaxed with the cell; pendant tilts as in the candidate screen",
                    "energy": "kcal/mol per kink pair: (E_kinked - E_TTTT) per cell / n_chains", "rate_proxy": "exp(-(E_form - E_form_ref) / kT) at 300 K, reference the CFE 70/20/10 row"})
    results.setdefault("rows", {})
    with PRESET_PARAMETERS[args.preset].applied():
        for name in args.rows.split(","):
            letters, note = (C.SEQUENCES.get(name) or EXTRA[name])
            polymer = C.make_polymer(name, letters)
            entry = results["rows"].setdefault(name, {"letters": letters, "note": note, "validated_chemistry": set(letters) <= C.VALIDATED})
            if "kinked" not in entry or "failed" in entry["kinked"]:
                try:
                    if args.same_basin:
                        p0 = reference_params(polymer, f"{name} TTTT")
                        entry["start_params"] = p0.tolist()
                        entry["kinked"] = one_cell_from(polymer, KINK, f"{name} kink", p0)
                    else:
                        entry["kinked"] = C.one_cell(polymer, KINK, f"{name} kink")
                except Exception as exc:
                    entry["kinked"] = {"failed": repr(exc)}
                    print(f"{name} kink: failed: {exc!r}", flush=True)
            if "TTTT" not in entry or "failed" in entry["TTTT"]:
                ref = screen.get(name, {}).get("conformations", {}).get("TTTT")
                if ref is None or "failed" in ref:
                    try:
                        ref = C.one_cell(polymer, [T] * 20, f"{name} TTTT")
                    except Exception as exc:
                        ref = {"failed": repr(exc)}
                entry["TTTT"] = ref
            k, r = entry["kinked"], entry["TTTT"]
            if "failed" not in k and "failed" not in r:
                n_mon = k["monomers_per_cell"]
                n_chains = 2
                e_form = (k["energy_per_monomer_kcal"] - r["energy_per_monomer_kcal"]) * n_mon / n_chains
                e_form_rigid = (k["energy_rigid_per_monomer_kcal"] - r["energy_rigid_per_monomer_kcal"]) * n_mon / n_chains
                entry["formation_kcal_per_kink_pair"] = e_form
                entry["formation_rigid_kcal_per_kink_pair"] = e_form_rigid
                entry["kink_tilt_release_kcal_per_cell"] = k["pendant_tilt_relaxation_per_monomer_kcal"] * n_mon
                entry["c_change_pct"] = 100.0 * (k["c_per_monomer_A"] / r["c_per_monomer_A"] - 1.0)
                print(f"{name:22s}: kink pair {e_form:+.2f} kcal/mol (rigid {e_form_rigid:+.2f}); chain shortens {entry['c_change_pct']:+.2f}%; |P| {k['P_abs_C_m2']:.3f} vs {r['P_abs_C_m2']:.3f}  ({k['seconds']:.0f} s)", flush=True)
            with open(out_path, "w") as fh:
                json.dump(results, fh, indent=1)
    ref_row = results["rows"].get("vdf7-trfe2-cfe1", {})
    e_ref = ref_row.get("formation_kcal_per_kink_pair")
    if e_ref is not None:
        for name, e in results["rows"].items():
            if "formation_kcal_per_kink_pair" in e:
                e["rate_factor_vs_cfe"] = float(np.exp(-(e["formation_kcal_per_kink_pair"] - e_ref) / KT))
        with open(out_path, "w") as fh:
            json.dump(results, fh, indent=1)
        print("rate proxy vs the CFE row:", {n: round(e["rate_factor_vs_cfe"], 3) for n, e in results["rows"].items() if "rate_factor_vs_cfe" in e})
    print(f"written {out_path}")


if __name__ == "__main__":
    main()
