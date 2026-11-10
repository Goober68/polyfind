"""The polished cells of the E2 addendum, in the schema of the E2 predictions.

examples/e2_predictions.py registered E2 on cells refined under the line group of the state
pattern, which a chiral sequence does not have.  examples/polish_margins.py relaxed them again
with every torsion and angle free (closure exact); those cells are the ones exported for Dipole
(``docs/predictions/polished_<row>_<cell>.xyz``) in the E2 addendum of
sarco/materials/ask-polyfind-dipole-registry-kink-v1.md.  This script writes their energies,
repeats and margins in the layout examples/crosscheck_e1_e2.py reads (``--e2``), with the hash
of each exported file, so the comparison needs nothing but Dipole's energies.  Each cell is
exported twice (examples/export_mixed_cell.py): with rigid pendant groups, whose model energy is
``energy_rigid_pendants_per_monomer_kcal`` (what a DFT single point on that file is compared
with), and with the groups tilt-relaxed (``--tilts``), whose model energy is
``energy_per_monomer_kcal`` (a single point on the tilted file, or a relaxation from either):

    python examples/e2_polished_predictions.py
    python examples/crosscheck_e1_e2.py --polyfind "" --dipole "" --e2 docs/predictions/e2_polished_polyfind_belt4.json --e2-dipole <file>

No polarization is registered for these cells (the mixed cell has no polarization routine yet).
"""
import argparse
import hashlib
import json
from pathlib import Path

SYSTEMS = {"vdftrfe": ("vdf7-trfe3", "VDF/TrFE 70/30 model, polished cells"),
           "vdftrfecfe_64_29_7": ("vdf9-trfe4-cfe1-14", "VDF/TrFE/CFE 64/29/7 (one CFE in fourteen), polished cells")}
FILE_TAG = {"TTTT": "TTTT", "TGTG'": "TGTGp", "TG'TG": "TGpTG"}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--margins", default="docs/predictions/polished_margins_belt4.json")
    ap.add_argument("--out", default="docs/predictions/e2_polished_polyfind_belt4.json")
    args = ap.parse_args()
    src = json.load(open(args.margins))
    base = Path(args.margins).parent
    out = {"potential": f"{src['preset']}, induced dipoles, Ewald, force-shifted LJ, fluorine belt sites, PBE-D3 chain geometry, rigid bonds, 0 K; "
                        "every torsion and backbone angle of both chains free, closure exact (examples/polish_margins.py)",
           "frame": "x polar, y lateral, z chain",
           "energy": "kcal/mol per monomer, with the pendant-tilt relaxation (the atoms of file_tilted); energy_rigid_pendants_per_monomer_kcal without it "
                     "(the atoms of file_rigid_pendants)",
           "systems": {}}
    for key, (row, label) in SYSTEMS.items():
        r = src["rows"][row]
        ground = r["tg_conformer"]
        cells = r["cells"]
        conf = {}
        for name, c in cells.items():
            a, b, gamma, cc = (float(v) for v in c["cell"])
            entry = {"energy_per_monomer_kcal": float(c["energy_per_monomer"]), "energy_rigid_pendants_per_monomer_kcal": float(c["rigid_per_monomer"]),
                     "pendant_tilt_relaxation_per_monomer_kcal": float(c["tilt_per_monomer"]), "c_per_monomer_A": float(c["c_per_monomer"]),
                     "cell": {"a": a, "b": b, "gamma_deg": gamma, "c": cc}, "n_chains": 2, "monomers_per_cell": 2 * len(r["letters"]),
                     "largest_torsion_off_ideal_deg": float(c["largest_torsion_off_ideal_deg"]), "gain_over_the_registered_cell_per_monomer_kcal": float(c["polish_gain_per_monomer"]),
                     "margin_vs_ground_kcal": float(c["energy_per_monomer"] - cells[ground]["energy_per_monomer"]),
                     "margin_rigid_pendants_vs_ground_kcal": float(c["rigid_per_monomer"] - cells[ground]["rigid_per_monomer"]),
                     "conformational_strain_vs_ground": float(c["c_per_monomer"] / cells[ground]["c_per_monomer"] - 1.0)}
            # the rigid-pendant atoms (the energy_rigid_pendants number) and the tilt-relaxed ones (the energy number)
            for field, suffix in (("file_rigid_pendants", ""), ("file_tilted", "_tilted")):
                path = base / f"polished_{row}_{FILE_TAG.get(name, name)}{suffix}.xyz"
                if path.exists():
                    entry[field] = {"name": path.name, "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
            conf[name] = entry
        out["systems"][key] = {"label": label, "sequence": r["letters"], "row": row, "ground_helix": ground, "margin_TT_minus_TG_kcal": float(r["margin"]),
                               "registered_margin_kcal": float(r["screen"]["margin"]), "conformations": conf}
        print(f"{key} ({r['letters']}): " + "; ".join(f"{n} {v['energy_per_monomer_kcal']:.4f} (rigid {v['energy_rigid_pendants_per_monomer_kcal']:.4f}), c {v['c_per_monomer_A']:.4f}"
                                                      + (f", {v['file_rigid_pendants']['name']}" + (" and _tilted" if "file_tilted" in v else "") if "file_rigid_pendants" in v else ", no exported file")
                                                      for n, v in conf.items())
              + f"; margin TTTT - {ground} {r['margin']:+.3f} (registered cells {r['screen']['margin']:+.3f})")
    Path(args.out).write_text(json.dumps(out, indent=1) + "\n")
    print(f"written {args.out}")


if __name__ == "__main__":
    main()
