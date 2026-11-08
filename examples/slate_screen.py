"""The slate's non-PVDF chemistries through the same crystal screen as examples/candidate_screen.py.

PVDC, polyacrylonitrile (AN), poly(vinylidene cyanide) (VDCN) and the 11:1 VDF/VDCN copolymer, with
PVDF as the control, under the plain valence preset ``pvdf-dft-valence`` (no belt sites exist for these
chemistries, so the fitted crystal terms of the belt presets do not apply; the nonbonded terms are the
valence fit's).  TTTT against TGTG': margin, switch strain, polarization of the TTTT crystal, dipole
gained per monomer and the threshold field.  Nitrile pendants are multi-atom fragments, which the
pendant-tilt relaxation does not treat, so those rows are rigid-pendant numbers (recorded as such).
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
from polyfind.polymers import POLYMERS

SPEC = importlib.util.spec_from_file_location("e2", Path(__file__).with_name("e2_predictions.py"))
E2 = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(E2)
CS = importlib.util.spec_from_file_location("cs", Path(__file__).with_name("candidate_screen.py"))
C = importlib.util.module_from_spec(CS)
CS.loader.exec_module(C)
T, GP, GM = 0, 1, 2

SLATE = {"pvdf": "the control, plain chain under the valence preset",
         "pvdc": "poly(vinylidene chloride), CH2-CCl2",
         "an": "polyacrylonitrile, CH2-CH(CN): the gen-1 slate's leader on dipole and degeneracy",
         "vdcn": "poly(vinylidene cyanide), CH2-C(CN)2: best degeneracy, dipoles largely cancelling",
         "vdf11-vdcn1": "the 11:1 VDF/VDCN copolymer (8.3 mol% VDCN)"}


def conformations(polymer):
    b = polymer.bonds_per_repeat
    tg = [T, GP, T, GM] * (b // 4) if b % 4 == 0 else [T, GP, T, GM] * (b // 2)
    return {"TTTT": [T] * max(b, 2), "TGTG'": tg}


def one_cell(polymer, seq, label):
    t0 = time.time()
    ref, rr = E2.build(polymer, seq, label)
    st = ES._state(ref, np.zeros(6), None, None, ES.TIGHT)
    pk = ref.packer
    n_chains = int(getattr(pk, "n_chains", 2))
    n_mon = n_chains * int(getattr(pk.chain, "n_monomers", 1))  # the chain counts vinyl monomers for homopolymers and copolymers alike
    p = np.asarray(ref.polarization(), dtype=float)
    bend_cell = float(rr.angle_energy) * n_chains * int(getattr(pk.chain, "n_monomers", 1))  # angle_energy is per vinyl monomer
    row = {"energy_rigid_per_monomer_kcal": (float(st.energy) + bend_cell) / n_mon, "P_C_m2": p.tolist(), "P_abs_C_m2": float(np.linalg.norm(p)),
           "c_per_monomer_A": float(ref.c) * n_chains / n_mon, "volume_per_monomer_A3": float(ref.volume) / n_mon,
           "cell": {"a": float(ref.params[0]), "b": float(ref.params[1]), "gamma_deg": float(ref.params[2]), "c": float(ref.c), "volume_A3": float(ref.volume)},
           "monomers_per_cell": n_mon, "dipole_per_monomer_eA": ES.dipole_per_monomer(float(np.linalg.norm(p)), float(ref.volume), n_mon)}
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--preset", default="pvdf-dft-valence")
    ap.add_argument("--polymers", default=",".join(SLATE))
    ap.add_argument("--out", default="docs/predictions/slate_screen_valence.json")
    args = ap.parse_args()
    out_path = Path(args.out)
    results = json.load(open(out_path)) if out_path.exists() else {}
    results.update({"preset": args.preset, "potential": f"{args.preset}, induced dipoles, Ewald, force-shifted LJ, plain chains (no belt sites), textbook geometry, pendant tilts where the pendants are single atoms, 0 K",
                    "frame": "x polar, y lateral, z chain"})
    results.setdefault("polymers", {})
    with PRESET_PARAMETERS[args.preset].applied():
        for name in args.polymers.split(","):
            polymer = POLYMERS[name]
            entry = results["polymers"].setdefault(name, {"note": SLATE.get(name, ""), "formula": polymer.formula, "conformations": {}})
            for conf, seq in conformations(polymer).items():
                if conf in entry["conformations"] and "failed" not in entry["conformations"][conf]:
                    continue
                try:
                    row = one_cell(polymer, seq, f"{name} {conf}")
                    print(f"{name:12s} {conf:6s}: E/monomer {row['energy_per_monomer_kcal']:.4f} (tilts {row['pendant_tilt_relaxation_per_monomer_kcal']:+.4f}) kcal/mol, |P| {row['P_abs_C_m2']:.4f} C/m^2, "
                          f"c/monomer {row['c_per_monomer_A']:.4f} A, V/monomer {row['volume_per_monomer_A3']:.2f} A^3, {row['monomers_per_cell']} monomers  ({row['seconds']:.0f} s)", flush=True)
                except Exception as exc:
                    row = {"failed": repr(exc)}
                    print(f"{name} {conf}: failed: {exc!r}", flush=True)
                entry["conformations"][conf] = row
                entry["derived"] = C.derive(entry["conformations"])
                with open(out_path, "w") as fh:
                    json.dump(results, fh, indent=1)
            d = entry["derived"]
            if d:
                mu = entry["conformations"]["TTTT"]["dipole_per_monomer_eA"]
                est = 0.0 if d["margin_TT_minus_TG_kcal_per_monomer"] <= 0 else ES.threshold_field(d["margin_TT_minus_TG_kcal_per_monomer"], mu) * 1e4
                print(f"  -> margin TT-TG {d['margin_TT_minus_TG_kcal_per_monomer']:+.3f} kcal/mol per monomer (rigid {d['margin_rigid_kcal_per_monomer']:+.3f}); axial {d['axial_transition_strain'] * 100:+.1f}%, "
                      f"areal {d['areal_transition_strain'] * 100:+.1f}%; P_TT {d['P_TT_C_m2']:.3f}; mu_TT {mu:.3f} e.A; E* (mu_TT) {est:.0f} MV/m", flush=True)
    print(f"written {out_path}")


if __name__ == "__main__":
    main()
