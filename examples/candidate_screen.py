"""Candidate chains and defect placements for the loss / frequency / strain optimum (first pass).

Ten-monomer model sequences packed as two-chain crystals in the TTTT (polar, beta-like) and TGTG'
(non-polar ground of PVDF) conformations under a fitted preset, each conformation's pendant groups
relaxed as rigid tilts (docs/CRYSTAL_STAGE.md 5.22), and from the pair: the switching margin
E(TTTT) - E(TGTG') per monomer, the transition strain of the switch (axial and areal), the dipole
gained per monomer and the threshold field E* = dE / dmu (polyfind.electrostriction.threshold_field).
Units: VDF and TrFE on the belt-site PBE-D3 geometry the preset was fitted on; CFE, CTFE and CDFE
with belt sites on their fluorines and the chlorine on the valence fit's parameters, which no
crystal anchor has tested -- those rows are extrapolations and are labelled so.
"""
import argparse
import importlib.util
import json
import time
from pathlib import Path

import numpy as np

from polyfind import electrostriction as ES
from polyfind.fitting import PRESET_PARAMETERS
from polyfind.pendants import tilt_relaxation
from polyfind.polymers import (NITRILE, PVDF, BackboneAtom, Polymer, belt_fluorine, copolymer, monomer_of, pvdf_belt, trfe_belt)

SPEC = importlib.util.spec_from_file_location("e2", Path(__file__).with_name("e2_predictions.py"))
E2 = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(E2)
T, GP, GM = 0, 1, 2
D = 0.7  # the belt-site distance of every fitted belt preset (crystalfit.BELT_D)


def cfe_belt(d: float = D) -> Polymer:
    """CH2-CFCl: PVDF's PBE-D3 C-H and C-F, the chlorine at 1.77 A; belt sites on the fluorine."""
    return Polymer(name="cfe-belt-pbe", formula="-(CH2-CFCl)n-", bond_length=1.5331,
                   backbone=(BackboneAtom("C", "H", 1.0960, 114.6, 107.64, -0.20, +0.10),
                             BackboneAtom("C", (belt_fluorine(1.3863, -0.20, d), "Cl"), (None, 1.77), 114.6, 108.0, +0.30, (None, -0.10))))


def ctfe_belt(d: float = D) -> Polymer:
    """CF2-CFCl."""
    return Polymer(name="ctfe-belt-pbe", formula="-(CF2-CFCl)n-", bond_length=1.5331,
                   backbone=(BackboneAtom("C", belt_fluorine(1.3863, -0.20, d), None, 114.6, 105.12, +0.40, None),
                             BackboneAtom("C", (belt_fluorine(1.3863, -0.20, d), "Cl"), (None, 1.77), 114.6, 108.0, +0.30, (None, -0.10))))


def cdfe_belt(d: float = D) -> Polymer:
    """CHCl-CF2 (the CFE isomer with the chlorine on the hydrogen-bearing carbon)."""
    return Polymer(name="cdfe-belt-pbe", formula="-(CHCl-CF2)n-", bond_length=1.5331,
                   backbone=(BackboneAtom("C", ("H", "Cl"), (1.0960, 1.77), 114.6, 109.0, 0.00, (+0.10, -0.10)),
                             BackboneAtom("C", belt_fluorine(1.3863, -0.20, d), None, 114.6, 105.12, +0.40, None)))


def an_belt() -> Polymer:
    """CH2-CH(CN), acrylonitrile, at the belt chain's PBE-D3 backbone geometry (the nitrile is a rigid
    fragment the pendant-tilt relaxation leaves alone; its parameters are the valence fit's, held out at
    0.9 kcal/mol)."""
    return Polymer(name="an-belt-pbe", formula="-(CH2-CH(CN))n-", bond_length=1.5331,
                   backbone=(BackboneAtom("C", "H", 1.0960, 114.6, 107.64, -0.20, +0.10),
                             BackboneAtom("C", ("H", NITRILE), (1.0960, None), 114.6, 109.5, +0.05, (+0.10, None))))


def fan_belt(d: float = D) -> Polymer:
    """CH2-CF(CN), 2-fluoroacrylonitrile (the gen-1 slate's highest dipole density), belt sites on its fluorine."""
    return Polymer(name="fan-belt-pbe", formula="-(CH2-CF(CN))n-", bond_length=1.5331,
                   backbone=(BackboneAtom("C", "H", 1.0960, 114.6, 107.64, -0.20, +0.10),
                             BackboneAtom("C", (belt_fluorine(1.3863, -0.20, d), NITRILE), (None, None), 114.6, 108.0, +0.30, (None, None))))


UNITS = {"V": lambda: monomer_of(pvdf_belt(D, "pbe"), "vdf"), "R": lambda: monomer_of(trfe_belt(D, "pbe"), "trfe"),
         "C": lambda: monomer_of(cfe_belt(), "cfe"), "K": lambda: monomer_of(ctfe_belt(), "ctfe"), "D": lambda: monomer_of(cdfe_belt(), "cdfe"),
         "A": lambda: monomer_of(an_belt(), "an"), "N": lambda: monomer_of(fan_belt(), "fan")}
VALIDATED = set("VR")  # chemistries the crystal anchors have tested

# name: (sequence letters, what the row asks)
SEQUENCES = {
    "vdf8-trfe2": ("VVRVVVVRVV", "P(VDF-TrFE) 80/20: less TrFE than the 70/30 model"),
    "vdf6-trfe4": ("VRVRVVRVRV", "P(VDF-TrFE) 60/40"),
    "vdf5-trfe5": ("VRVRVRVRVR", "P(VDF-TrFE) 50/50, alternating"),
    "vdf7-trfe2-cfe1": ("VVRVVCVVRV", "terpolymer model 70/20/10, the CFE two VDF away from either TrFE"),
    "vdf7-trfe2-cfe1-adj": ("VVRCVVVVRV", "the same composition with the CFE next to a TrFE"),
    "vdf7-trfe2-ctfe1": ("VVRVVKVVRV", "CTFE in place of CFE"),
    "vdf7-trfe2-cdfe1": ("VVRVVDVVRV", "CDFE (CHCl-CF2) in place of CFE"),
    "vdf9-cfe1": ("VVVVCVVVVV", "one CFE per ten VDF, no TrFE"),
    "vdf8-cfe2": ("VVCVVVVCVV", "two CFE per ten, no TrFE"),
    "vdf7-trfe2-an1": ("VVRVVAVVRV", "the slate's nitrile: one acrylonitrile per ten beside 20 mol% TrFE (Dipole's makeable track)"),
    "vdf7-trfe2-fan1": ("VVRVVNVVRV", "2-fluoroacrylonitrile in place of CFE (the slate's dipole-density anomaly)"),
    "vdf9-an1": ("VVVVAVVVVV", "one acrylonitrile per ten VDF, no TrFE: a trans-favouring defect pulling PVDF's margin toward zero from the other side"),
}
CONFORMATIONS = {"TTTT": [T] * 20, "TGTG'": [T, GP, T, GM] * 5}


def make_polymer(name: str, letters: str) -> Polymer:
    units = [UNITS[c]() for c in letters]
    formula = "-" + "-".join({"V": "(CH2-CF2)", "R": "(CHF-CF2)", "C": "(CH2-CFCl)", "K": "(CF2-CFCl)", "D": "(CHCl-CF2)", "A": "(CH2-CH(CN))", "N": "(CH2-CF(CN))"}[c] for c in letters) + "-"
    return copolymer(name=f"{name}-belt-pbe", formula=formula, units=units, bond_length=PVDF.bond_length)


def one_cell(polymer, seq, label):
    t0 = time.time()
    ref, rr = E2.build(polymer, seq, label)
    st = ES._state(ref, np.zeros(6), None, None, ES.TIGHT)
    n_chains = int(getattr(ref.packer, "n_chains", 2))
    n_mon = n_chains * int(getattr(ref.packer.chain, "n_monomers", len(seq) // 2))  # vinyl monomers; a carbon count would include nitrile carbons
    p = np.asarray(ref.polarization(), dtype=float)
    bend_cell = float(rr.angle_energy) * n_chains * float(getattr(ref.packer.chain, "n_monomers", len(seq) / 2.0))
    row = {"energy_rigid_per_monomer_kcal": (float(st.energy) + bend_cell) / n_mon, "P_C_m2": p.tolist(), "P_abs_C_m2": float(np.linalg.norm(p)),
           "c_per_monomer_A": float(ref.c) * n_chains / n_mon, "volume_per_monomer_A3": float(ref.volume) / n_mon,
           "cell": {"a": float(ref.params[0]), "b": float(ref.params[1]), "gamma_deg": float(ref.params[2]), "c": float(ref.c), "volume_A3": float(ref.volume)},
           "monomers_per_cell": n_mon, "dipole_per_monomer_eA": ES.dipole_per_monomer(float(np.linalg.norm(p)), float(ref.volume), int(round(n_mon)))}
    try:
        tr = tilt_relaxation(ref, polymer)
        row["pendant_tilt_relaxation_per_monomer_kcal"] = float(tr.dE_per_cell) / n_mon
        row["pendant_tilts_deg"] = tr.tilts_by_substituents
    except Exception as exc:
        row["pendant_tilt_relaxation_per_monomer_kcal"] = 0.0
        row["pendant_tilts_deg"] = f"unavailable: {exc!r}"
    row["energy_per_monomer_kcal"] = row["energy_rigid_per_monomer_kcal"] + row["pendant_tilt_relaxation_per_monomer_kcal"]
    row["seconds"] = time.time() - t0
    return row


def derive(rows: dict) -> dict:
    tt, tg = rows.get("TTTT"), rows.get("TGTG'")
    if not tt or not tg or "failed" in tt or "failed" in tg:
        return {}
    de = tt["energy_per_monomer_kcal"] - tg["energy_per_monomer_kcal"]
    dmu = tt["dipole_per_monomer_eA"] - tg["dipole_per_monomer_eA"]
    out = {"margin_TT_minus_TG_kcal_per_monomer": de, "margin_rigid_kcal_per_monomer": tt["energy_rigid_per_monomer_kcal"] - tg["energy_rigid_per_monomer_kcal"],
           "axial_transition_strain": tt["c_per_monomer_A"] / tg["c_per_monomer_A"] - 1.0,
           "areal_transition_strain": (tt["volume_per_monomer_A3"] / tt["c_per_monomer_A"]) / (tg["volume_per_monomer_A3"] / tg["c_per_monomer_A"]) - 1.0,
           "dipole_gain_eA_per_monomer": dmu, "P_TT_C_m2": tt["P_abs_C_m2"], "P_TG_C_m2": tg["P_abs_C_m2"]}
    if de > 0 and dmu > 0:
        out["threshold_field_MV_m"] = ES.threshold_field(de, dmu) * 1e4
    elif de <= 0:
        out["threshold_field_MV_m"] = 0.0  # TTTT is already the ground conformation
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--preset", default="pvdf-dft-valence-belt4")
    ap.add_argument("--sequences", default=",".join(SEQUENCES))
    ap.add_argument("--conformations", default=",".join(CONFORMATIONS))
    ap.add_argument("--out", default="docs/predictions/candidate_screen_belt4.json")
    args = ap.parse_args()
    params = PRESET_PARAMETERS[args.preset]
    out_path = Path(args.out)
    results = json.load(open(out_path)) if out_path.exists() else {}
    results.update({"preset": args.preset, "potential": f"{args.preset}, induced dipoles, Ewald, force-shifted LJ, belt sites on every fluorine, PBE-D3 chain geometry, pendant groups relaxed as rigid tilts, 0 K",
                    "energy": "lattice energy plus the refined backbone angles' strain plus the pendant-tilt relaxation, per vinyl monomer",
                    "frame": "x polar, y lateral, z chain"})
    results.setdefault("sequences", {})
    with params.applied():
        for name in args.sequences.split(","):
            letters, note = SEQUENCES[name]
            polymer = make_polymer(name, letters)
            entry = results["sequences"].setdefault(name, {"letters": letters, "note": note, "validated_chemistry": set(letters) <= VALIDATED, "conformations": {}})
            entry["validated_chemistry"] = set(letters) <= VALIDATED
            for conf in args.conformations.split(","):
                if conf in entry["conformations"] and "failed" not in entry["conformations"][conf]:
                    continue
                try:
                    row = one_cell(polymer, CONFORMATIONS[conf], f"{name} {conf}")
                    tilt = row["pendant_tilt_relaxation_per_monomer_kcal"]
                    print(f"{name:22s} {conf:6s}: E/monomer {row['energy_per_monomer_kcal']:.4f} (tilts {tilt:+.4f}) kcal/mol, |P| {row['P_abs_C_m2']:.4f} C/m^2, "
                          f"c/monomer {row['c_per_monomer_A']:.4f} A, V/monomer {row['volume_per_monomer_A3']:.2f} A^3  ({row['seconds']:.0f} s)", flush=True)
                except Exception as exc:
                    row = {"failed": repr(exc)}
                    print(f"{name} {conf}: failed: {exc!r}", flush=True)
                entry["conformations"][conf] = row
                entry["derived"] = derive(entry["conformations"])
                with open(out_path, "w") as fh:
                    json.dump(results, fh, indent=1)
            d = entry["derived"]
            if d:
                print(f"  -> margin TT-TG {d['margin_TT_minus_TG_kcal_per_monomer']:+.3f} kcal/mol per monomer (rigid {d['margin_rigid_kcal_per_monomer']:+.3f}); "
                      f"axial strain {d['axial_transition_strain'] * 100:+.1f}%, areal {d['areal_transition_strain'] * 100:+.1f}%; dmu {d['dipole_gain_eA_per_monomer']:+.3f} e.A; "
                      f"E* {d.get('threshold_field_MV_m', float('nan')):.0f} MV/m", flush=True)
    print(f"written {out_path}")


if __name__ == "__main__":
    main()
