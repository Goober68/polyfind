"""The packing search's candidate basins for a cell, each refined under the full potential.

    python examples/copolymer_basins.py --preset pvdf-dft-valence-belt2 --variant belt --geometry pbe
    python examples/copolymer_basins.py --preset pvdf-dft-valence --system beta --top 3

``mechanics.refined_reference`` refines the first of the distinct minima the packing screen returns
(``pack()``: an exhaustive table screen under the truncated Coulomb kernel, the best few polished and
deduplicated), and the E1/E2 sets stand on that one basin.  The 70/30 copolymer model cell landed in
a different arrangement under each preset (docs/CRYSTAL_STAGE.md 5.10, 5.14).  This script takes the
top ``--top`` candidates of the screen, refines each under the same full potential the sets use
(induced dipoles, Ewald, force-shifted Lennard-Jones, the chain's conformation and the cell), relaxes
the cell, and lists them by their full-potential energy per monomer with their cells and dipole
directions, so that the basin a set sits in can be named and the question "is it the lowest?" has an
answer.  Output: ``docs/predictions/basins_<preset>_<system>.json``.
"""
import argparse
import json
import time

import numpy as np

from polyfind import mechanics as M
from polyfind.ewald import EwaldSpec
from polyfind.linegroup import repeat_chains
from polyfind.pack import pack, periodic_chain
from polyfind.polarizability import Polarizable
from polyfind.polymers import PVDF, THREE_STATE, VDF_TRFE_7_3, belt_variant
from polyfind.refine import refine_crystal

SYSTEMS = {"vdftrfe": ("VDF/TrFE 70/30 model (TTTT)", VDF_TRFE_7_3, [0] * 20), "beta": ("beta-PVDF (TTTT)", PVDF, [0, 0]),
           "alpha": ("alpha-PVDF (TGTG')", PVDF, [0, 1, 0, 2])}
KW = {"lj_cutoff": "force", "coulomb": "ewald", "ewald": EwaldSpec(), "polarizable": Polarizable()}


def refine_candidate(polymer, res, label):
    """``mechanics.refined_reference``'s tail for one PackResult: refine, rebuild aligned, relax the cell."""
    rr = refine_crystal(polymer, res, **KW)
    angles0 = np.array([polymer.backbone[k].backbone_angle for k in range(polymer.bonds_per_repeat)])
    frame = repeat_chains(polymer, res.chain, np.asarray(res.dihedrals, dtype=float)[None], angles0)[0]
    built = repeat_chains(polymer, res.chain, rr.torsions[None], rr.angles, align_to=frame)[0]
    ref = M.reference_from_chain(built, rr.result.params, n_chains=res.n_chains, label=label, **KW)
    ref = M.relax_reference(ref)
    return ref, rr


def write_extxyz(ref, path: str, comment: str) -> None:
    """The cell's atoms (belt sites left out) in extended XYZ with the lattice, polyfind's frame."""
    from polyfind.phonon import placed_coordinates
    pk = ref.packer
    Pn, latn = placed_coordinates(pk, np.asarray(ref.params, dtype=float))
    Pn, latn = np.asarray(Pn, dtype=float), np.asarray(latn, dtype=float)
    el = np.asarray(list(pk.elements) * pk.n_chains)
    keep = el != "X"
    lat = " ".join(f"{v:.6f}" for v in latn.reshape(-1))
    with open(path, "w") as fh:
        fh.write(f"{int(keep.sum())}\n")
        fh.write(f'Lattice="{lat}" Properties=species:S:1:pos:R:3 pbc="T T T" {comment}\n')
        for e, r in zip(el[keep], Pn[keep]):
            fh.write(f"{e:2s} {r[0]:14.8f} {r[1]:14.8f} {r[2]:14.8f}\n")


def describe(ref, rr, res):
    pk = ref.packer
    n_mon = pk.n_chains * pk.chain.n_monomers
    e_cell = float(pk.energy(ref.params[None], c=np.array([ref.c]))[0])
    P = np.asarray(ref.polarization(), dtype=float)
    tilt = float(np.degrees(np.arctan2(P[1], P[0])))
    return {"screen": {"energy_per_monomer_dsf": float(res.energy_per_monomer), "a": float(res.a), "b": float(res.b), "gamma": float(res.gamma),
                       "phi1": float(res.phi1), "phi2": float(res.phi2), "dz": float(res.dz), "flip": int(res.flip)},
            "refined": {"energy_per_monomer": e_cell / n_mon, "angle_strain_per_monomer": float(rr.angle_energy),
                        "total_per_monomer": e_cell / n_mon + float(rr.angle_energy),
                        "a": float(ref.params[0]), "b": float(ref.params[1]), "gamma": float(ref.params[2]), "c": float(ref.c),
                        "phi1": float(ref.params[3]), "phi2": float(ref.params[4]), "dz": float(ref.params[5]), "flip": int(round(ref.params[6])),
                        "P_C_m2": P.tolist(), "P_abs": float(np.linalg.norm(P)), "tilt_from_a_deg": tilt,
                        "backbone_angles_deg": np.asarray(rr.angles, dtype=float).tolist()}}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--preset", default="pvdf-dft-valence-belt2")
    ap.add_argument("--variant", default="plain", choices=("plain", "belt"))
    ap.add_argument("--geometry", default="textbook", choices=("textbook", "pbe"))
    ap.add_argument("--system", default="vdftrfe", choices=sorted(SYSTEMS))
    ap.add_argument("--top", type=int, default=4, help="how many of the screen's distinct minima to refine")
    ap.add_argument("--out", default=None)
    ap.add_argument("--xyz", action="store_true", help="also write each refined basin's cell as extended XYZ beside the output")
    args = ap.parse_args()
    from polyfind.fitting import PRESET_PARAMETERS
    label, polymer, seq = SYSTEMS[args.system]
    if args.variant == "belt":
        polymer = belt_variant(polymer, geometry=args.geometry)
    out_path = args.out or f"docs/predictions/basins_{args.preset}_{args.system}.json"
    results = {"preset": args.preset, "variant": args.variant, "geometry": args.geometry, "system": label,
               "frame": "x polar (a), y lateral (b), z chain", "candidates": []}
    with PRESET_PARAMETERS[args.preset].applied():
        t0 = time.time()
        chain = periodic_chain(polymer, seq, THREE_STATE)
        cands = pack(chain, n_chains=2, n_refine=max(args.top, 5), table_cache_dir=None)
        print(f"screen: {len(cands)} distinct minima under the truncated kernel ({time.time() - t0:.0f} s)")
        for r in cands:
            print("  " + r.row())
        for k, res in enumerate(cands[:args.top]):
            t1 = time.time()
            ref, rr = refine_candidate(polymer, res, f"{label} basin {k}")
            d = describe(ref, rr, res)
            d["index"] = k
            d["seconds"] = time.time() - t1
            results["candidates"].append(d)
            f = d["refined"]
            print(f"basin {k}: refined E/monomer {f['energy_per_monomer']:.4f} (+ angle strain {f['angle_strain_per_monomer']:.4f} = {f['total_per_monomer']:.4f}) kcal/mol; "
                  f"cell a {f['a']:.3f} b {f['b']:.3f} gamma {f['gamma']:.1f} c {f['c']:.3f}; |P| {f['P_abs']:.3f} C/m^2 at {f['tilt_from_a_deg']:.0f} deg from a; "
                  f"phi ({f['phi1']:.1f}, {f['phi2']:.1f}) dz {f['dz']:.2f} {'anti' if f['flip'] else 'para'}   ({d['seconds']:.0f} s)", flush=True)
            if args.xyz:
                xyz = out_path.replace(".json", f"_basin{k}.xyz")
                write_extxyz(ref, xyz, f"polyfind_preset={args.preset} system=\"{label}\" basin={k} "
                                       f"E_per_monomer_kcal={f['total_per_monomer']:.5f} P_C_m2=\"{f['P_C_m2'][0]:.4f} {f['P_C_m2'][1]:.4f} {f['P_C_m2'][2]:.4f}\" "
                                       f"frame=\"x polar, y lateral, z chain\"")
                d["xyz"] = xyz
            with open(out_path, "w") as fh:
                json.dump(results, fh, indent=1)
    order = sorted(results["candidates"], key=lambda d: d["refined"]["total_per_monomer"])
    print("\nby full-potential energy (with the angle strain):")
    for d in order:
        f = d["refined"]
        print(f"  basin {d['index']}: {f['total_per_monomer']:.4f} kcal/mol per monomer, gamma {f['gamma']:.1f}, tilt {f['tilt_from_a_deg']:.0f} deg, |P| {f['P_abs']:.3f}")
    print(f"written {out_path}")


if __name__ == "__main__":
    main()
