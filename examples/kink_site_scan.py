"""The energy of a kink sheet as a function of where along the sequence it sits.

In PVDF every monomer is the same and a kink pair costs the same anywhere.  In a copolymer the
kink pair T G+ T G- T sits on particular units (VDF, TrFE, a chlorinated defect), and the sheet of
kink pairs (every chain kinked in register, examples/kink_migration.py) has a different energy at
each position: that profile is what pins the sheet and what it must climb to move.

Starting from the packing search's kinked crystal at ``--kink-at`` (polished: every torsion free,
closure exact), the kink pair of every chain is moved along the chain two bonds (one monomer) at
a time by a local edit of its three torsions, and the cell is relaxed again with everything free.
The lattice is therefore the one the sheet carries with it, not a new packing search at each
position: the nearest minimum, which is what a moving sheet sees.  A chain packed head-to-tail
against its neighbour (a flipped chain) has its kink moved the other way along its own sequence,
so that the sheet stays one plane in space.

Energies are in kcal/mol per kink pair against the polished all-trans crystal, rigid pendants.
Relaxing the cell at each position could mix two things, the site's own energy and the packing
basin the relaxation falls in.  In the 71/29 cell, where neighbouring VDF-VDF positions come out
2 kcal/mol apart in cells 0.15 A apart in a, it does not: four sites re-relaxed from another
site's lattice returned to the same energies and the same cells, so the lattice follows the site.
``--fixed-cell`` is there for a lattice where that check fails: it holds the cell at the starting
kinked crystal's and relaxes only the chains.
``--parity 1`` puts the gauche bonds on the other kind of backbone bond (offsets ``--kink-at`` + 1,
+ 3, ...), which is a different kink pair.
"""
import argparse
import importlib.util
import json
import time
from pathlib import Path

import numpy as np

from polyfind.fitting import PRESET_PARAMETERS

MS = importlib.util.spec_from_file_location("km", Path(__file__).with_name("kink_migration.py"))
KM = importlib.util.module_from_spec(MS)
MS.loader.exec_module(KM)
DK = KM.DK


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--row", default="vdf9-trfe4-cfe1-14")
    ap.add_argument("--kink-at", type=int, default=0, help="the offset of the stored kinked crystal the scan starts from")
    ap.add_argument("--mirror-kink", action="store_true")
    ap.add_argument("--parity", type=int, default=0, choices=(0, 1), help="0: offsets kink-at + 0, 2, 4, ...; 1: kink-at + 1, 3, ...")
    ap.add_argument("--offsets", default="", help="comma list of offsets (bonds) to visit; default every position of the parity, walked outward from the start")
    ap.add_argument("--preset", default="pvdf-dft-valence-belt4")
    ap.add_argument("--maxiter", type=int, default=600)
    ap.add_argument("--fixed-cell", action="store_true", help="hold the cell at the starting kinked crystal's; relax only the chains")
    ap.add_argument("--references", default="docs/predictions/dilute_kink_belt4.json")
    ap.add_argument("--kinked-references", default="docs/predictions/kink_migration_belt4.json")
    ap.add_argument("--out", default="docs/predictions/kink_site_scan_belt4.json")
    args = ap.parse_args()
    log = lambda s: print(s, flush=True)  # noqa: E731
    out_path = Path(args.out)
    letters, note = DK.sequence(args.row)
    key = f"{args.row}@{args.kink_at}{'~' if args.mirror_kink else ''}"
    out_key = key + (" fixed cell" if args.fixed_cell else "")  # the stored kinked crystal is under ``key``; the scan's rows under ``out_key``

    def locked_update(path: Path, fn):
        import fcntl
        with open(str(path) + ".lock", "w") as lock:
            fcntl.flock(lock, fcntl.LOCK_EX)
            disk = json.load(open(path)) if path.exists() else {}
            fn(disk)
            with open(path, "w") as fh:
                json.dump(disk, fh, indent=1)

    with PRESET_PARAMETERS[args.preset].applied():
        polymer = DK.C.make_polymer(args.row, letters)
        P = polymer.bonds_per_repeat
        ref_path, kref_path = Path(args.references), Path(args.kinked_references)
        stored = (json.load(open(ref_path)) if ref_path.exists() else {}).get("all_trans_references", {}).get(args.row)
        d_tt = DK.all_trans_reference(polymer, args.row, stored)
        if stored is None or "state" not in stored:
            locked_update(ref_path, lambda disk: disk.setdefault("all_trans_references", {}).__setitem__(args.row, DK.as_json(d_tt)))
        stored = (json.load(open(kref_path)) if kref_path.exists() else {}).get("kinked_references", {}).get(key)
        g = 2 if args.mirror_kink else 1
        seq = [int(v) for v in np.roll(np.array([g, 0, 3 - g] + [0] * (P - 3)), args.kink_at % P)]
        d = DK.refined_cell(polymer, args.row, seq, f"{args.row} kink@{args.kink_at}{'~' if args.mirror_kink else ''}", stored)
        if stored is None or "state" not in stored:
            locked_update(kref_path, lambda disk: disk.setdefault("kinked_references", {}).__setitem__(key, DK.as_json(d)))
        model, xA, ch = DK.make_model(polymer, d, 1, 1)
        nc, n_mon = model.nc, ch.n_monomers
        e_tt = float(d_tt["energy"])
        held = [model.index("xy", nc - 1), model.index("z", nc - 1)] + ([model.index("lattice"), model.index("c")] if args.fixed_cell else [])
        free = np.setdiff1d(np.arange(model.n_var), np.concatenate(held))
        everyone = list(range(nc))
        sense = np.where(model.flips, -1, 1)  # a flipped chain runs the other way along the cell's axis
        prior = (json.load(open(out_path)) if out_path.exists() else {}).get("rows", {}).get(out_key, {}).get("sites", {})
        sites = dict(prior)
        e0 = (model.physical_energy(xA) - e_tt) / nc
        sites[str(args.kink_at)] = {"formation": e0, "relative": 0.0, "units": letters[(args.kink_at // 2) % n_mon] + letters[(args.kink_at // 2 + 1) % n_mon],
                                    "torsions": KM.letters_of(model.state(xA).tors[0]), "cell": model.state(xA).lattice.tolist(), "converged": True, "state": xA.tolist()}
        log(f"{out_key}: {letters}; kink sheet at offset {args.kink_at} (units {sites[str(args.kink_at)]['units']}): {e0:+.3f} kcal/mol per kink pair against all-trans")

        def save():
            locked_update(out_path, lambda disk: (disk.update(preset=args.preset, energy="kcal/mol per kink pair against the polished all-trans crystal, rigid pendants; "
                                                              "the sheet moved along the sequence from the stored kinked crystal, cell relaxed at each position"),
                                                  disk.setdefault("rows", {}).__setitem__(out_key, {"letters": letters, "note": note, "start_offset": args.kink_at,
                                                                                                   "mirrored": bool(args.mirror_kink), "fixed_cell": bool(args.fixed_cell), "sites": sites})))
        save()
        if args.offsets:
            targets = [int(v) for v in args.offsets.split(",")]
        else:  # outward from the start, one monomer at a time in each direction, so each start is the neighbouring site's relaxed cell
            half = P // 2
            steps = [s for k in range(1, half // 2 + 1) for s in (2 * k, -2 * k)]
            targets = [(args.kink_at + args.parity + s) % P for s in ([0] if args.parity else []) + steps]
            targets = list(dict.fromkeys(t for t in targets if t != args.kink_at % P or args.parity))
        for p in targets:
            if str(p) in sites and sites[str(p)].get("converged"):
                continue
            # start from the nearest site already relaxed (along the ring of positions)
            done = {int(k): v for k, v in sites.items() if v.get("converged") and "state" in v}
            near = min(done, key=lambda q: min((p - q) % P, (q - p) % P))
            step = (p - near + P // 2) % P - P // 2
            x = np.asarray(done[near]["state"], dtype=float).copy()
            for s in range(nc):
                idx = model.index("tors", s)
                x[idx] = KM.move_kink(x[idx], near if sense[s] > 0 else (2 * args.kink_at - near) % P, int(step * sense[s]))
            t0 = time.time()
            x = model.project(x, model.index("internal"), everyone)
            start = (model.physical_energy(x) - e_tt) / nc
            xr, info = model.relax(x, free=free, closed=everyone, maxiter=args.maxiter)
            for _ in range(3):
                if "success" in info["message"].lower():
                    break
                xr, more = model.relax(xr, free=free, closed=everyone, maxiter=args.maxiter)
                info = {**more, "iterations": info["iterations"] + more["iterations"]}
            st = model.state(xr)
            e = (model.physical_energy(xr) - e_tt) / nc
            lets = KM.letters_of(st.tors[0])
            kept = lets.split()[p % P] != "T" and lets.split()[(p + 2) % P] != "T" and lets.count("G") == 2
            sites[str(p)] = {"formation": e, "relative": e - e0, "units": letters[(p // 2) % n_mon] + letters[(p // 2 + 1) % n_mon], "torsions": lets,
                             "kink_kept": bool(kept), "cell": st.lattice.tolist(), "start_energy": start, "from_offset": int(near), "iterations": info["iterations"],
                             "message": info["message"], "converged": "success" in info["message"].lower(), "seconds": time.time() - t0, "state": xr.tolist()}
            log(f"   offset {p:3d} (units {sites[str(p)]['units']}, bonds {p}-{(p + 2) % P}): {e:+.3f} per kink pair ({e - e0:+.3f} from the start site){'' if kept else '  KINK NOT KEPT: ' + lets}; "
                f"cell {np.round(st.lattice, 3).tolist()}; from offset {near}, {info['iterations']} it, {time.time() - t0:.0f} s; {info['message']}")
            save()
    log(f"written {out_path}")


if __name__ == "__main__":
    main()
