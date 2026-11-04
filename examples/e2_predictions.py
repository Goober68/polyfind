"""polyfind's pre-registered predictions for E2, the two-phase energetics per chemistry.

    python examples/e2_predictions.py [--out docs/predictions/e2_polyfind.json] [--systems pvdf,vdftrfe]

Fully relaxed cells of the three conformations TTTT, TGTG' and T3GT3G' for PVDF and for the
70/30 VDF/TrFE model, on the fitted polarizable potential: energy per monomer and the
margins between conformations (kcal/mol per monomer), the polarization, ``c`` per monomer
along the chain and the conformational strain ``c_polar / c_ground - 1``, the cell.  The
defect-unit cells of the ask (CFE, CTFE, CDFE, FA, DB in a beta host) are not here: the
potential's chlorine parameters are not validated, which is one of the things E2 tests.
"""
import argparse
import json
import time

import numpy as np

from polyfind import electrostriction as ES
from polyfind import mechanics as M
from polyfind.polymers import PVDF, THREE_STATE, VDF_TRFE_7_3

T, GP, GM = 0, 1, 2
SYSTEMS = {
    "pvdf": ("PVDF", PVDF, {"TTTT": [T, T], "TGTG'": [T, GP, T, GM], "T3GT3G'": [T, T, T, GP, T, T, T, GM]}),
    "vdftrfe": ("VDF/TrFE 70/30 model", VDF_TRFE_7_3,
                {"TTTT": [T] * 20, "TGTG'": [T, GP, T, GM] * 5, "T3GT3G'": [T, T, T, GP, T, T, T, GM] * 5}),
}


def build(polymer, seq, label):
    from polyfind.ewald import EwaldSpec
    from polyfind.polarizability import Polarizable
    kw = {"pack_kw": {"table_cache_dir": None}, "lj_cutoff": "force", "coulomb": "ewald", "ewald": EwaldSpec(),
          "polarizable": Polarizable()}
    ref, _ = M.refined_reference(polymer, seq, label=label, **kw)
    return M.relax_reference(ref)


def monomers_in_cell(polymer, seq, ref):
    """Vinyl monomers in the packed cell, counted from the backbone carbons (two per monomer);
    the helix routine's "monomer" is the chain's repeat unit, which for a copolymer is the
    whole 10-unit sequence, so it is not used here."""
    pk = ref.packer
    n_chains = int(getattr(pk, "n_chains", 2))
    els = [str(e) for e in np.ravel(np.asarray(getattr(pk, "elements", []), dtype=object))]
    n_c = sum(1 for e in els if e.upper() == "C")
    per_chain_from_seq = len(seq) / 2.0  # two backbone bonds per vinyl monomer, one chain period
    if n_c > 0:
        mon = n_c / 2.0
        if abs(mon - per_chain_from_seq) < 0.5:  # the packer lists one chain's atoms
            mon *= n_chains
    else:
        mon = per_chain_from_seq * n_chains
    return float(mon), n_chains


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="docs/predictions/e2_polyfind.json")
    ap.add_argument("--systems", default="pvdf,vdftrfe")
    args = ap.parse_args()
    from polyfind.fitting import FITTED_VALENCE
    results = {"potential": "pvdf-dft-valence, induced dipoles, Ewald, rigid chain, force-shifted LJ, 0 K",
               "frame": "x polar, y lateral, z chain", "systems": {}}
    with FITTED_VALENCE.applied():
        for key in args.systems.split(","):
            label, polymer, confs = SYSTEMS[key]
            rows = {}
            for conf, seq in confs.items():
                t0 = time.time()
                try:
                    ref = build(polymer, seq, f"{label} {conf}")
                    st = ES._state(ref, np.zeros(6), None, None, ES.TIGHT)
                    n_mon, n_chains = monomers_in_cell(polymer, seq, ref)
                    p = np.asarray(ref.polarization(), dtype=float)
                    rows[conf] = {"sequence": seq, "n_chains": n_chains, "monomers_per_cell": n_mon,
                                  "energy_per_monomer_kcal": float(st.energy) / n_mon,
                                  "P_C_m2": p.tolist(), "c_per_monomer_A": float(ref.c) * n_chains / n_mon,
                                  "cell": {"a": float(ref.params[0]), "b": float(ref.params[1]), "gamma_deg": float(ref.params[2]),
                                           "c": float(ref.c), "volume_A3": float(ref.volume)},
                                  "volume_per_monomer_A3": float(ref.volume) / n_mon, "seconds": time.time() - t0}
                    print(f"{label} {conf:8s}: E/monomer {rows[conf]['energy_per_monomer_kcal']:.4f} kcal/mol, |P| {np.linalg.norm(p):.4f} C/m^2, "
                          f"c/monomer {rows[conf]['c_per_monomer_A']:.4f} A, V/monomer {rows[conf]['volume_per_monomer_A3']:.2f} A^3, "
                          f"cell a {ref.params[0]:.3f} b {ref.params[1]:.3f} gamma {ref.params[2]:.1f} c {ref.c:.3f}   ({time.time() - t0:.0f} s)", flush=True)
                except Exception as exc:  # a conformation that does not pack is a result too
                    rows[conf] = {"sequence": seq, "failed": repr(exc)}
                    print(f"{label} {conf}: failed: {exc!r}", flush=True)
            ok = {k: v for k, v in rows.items() if "failed" not in v}
            if ok:
                ground = min(ok, key=lambda k: ok[k]["energy_per_monomer_kcal"])
                for conf, r in ok.items():
                    r["margin_vs_ground_kcal"] = r["energy_per_monomer_kcal"] - ok[ground]["energy_per_monomer_kcal"]
                    r["conformational_strain_vs_ground"] = r["c_per_monomer_A"] / ok[ground]["c_per_monomer_A"] - 1.0
                print(f"  ground conformation {ground}; margins (kcal/mol per monomer): "
                      + ", ".join(f"{c} {r['margin_vs_ground_kcal']:+.3f} (strain {r['conformational_strain_vs_ground'] * 100:+.1f}%)" for c, r in ok.items()), flush=True)
            results["systems"][key] = {"label": label, "conformations": rows}
            with open(args.out, "w") as fh:
                json.dump(results, fh, indent=1)


if __name__ == "__main__":
    main()
