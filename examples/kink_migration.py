"""The barrier for one chain's kink pair to move one monomer along the chain, by nudged elastic band.

A kink pair T G+ T G- T shortens its chain by about 1 A, and a chain is held in axial registry
with its neighbours (examples/dilute_kink.py: a lone kink pair in an all-trans crystal carries a
strain field that costs more than the kink).  The kink pairs of a crystal therefore sit in
register, on every chain at the same height, and that kinked crystal is the host here: every
chain of the supercell carries one kink pair per repeat, so no chain is shorter than another, and
ONE chain's kink pair is moved one monomer ahead of its neighbours'.

The move is two gauche migrations (counter-rotations of two second-neighbour bonds, the only
torsional move that leaves both stems of the chain parallel), through the wide kink:

    A   G+ T G- T T     the kink pair in register with the neighbours' (the kinked crystal)
    I   G+ T T T G-     the leading gauche bond moved two bonds on
    B   T T G+ T G-     the trailing one followed: the kink pair one monomer ahead

The end points A, I, B are relaxed at the kinked crystal's cell with every chain periodic exactly
(polyfind.mixedcell), and the two segments A -> I and I -> B are relaxed as bands on the closure
surface with a climbing image.  A is the kinked crystal of the packing search
(examples/kink_formation_screen.py's own-basin cell: a relaxation started from the all-trans
lattice stops in a symmetric cell that lies higher), polished with every torsion free and the
closure exact, as the all-trans reference it is measured against is (examples/dilute_kink.py).  Reported per row: the energies of I and B above A (the price of
leaving the sheet by one monomer), the two saddle energies, and so the forward barrier of the
one-monomer step, in kcal/mol per chain that moves.

``--host`` sets what may respond: ``rigid`` (the neighbours fixed as they are in the kinked
crystal), ``bodies`` (each neighbour free to shift and turn as a rigid chain) or ``flexible``
(every torsion and backbone angle of every chain).
"""
import argparse
import importlib.util
import json
import multiprocessing as mp
import time
from pathlib import Path

import numpy as np

from polyfind import mixedcell as MC
from polyfind.fitting import PRESET_PARAMETERS

DS = importlib.util.spec_from_file_location("dk", Path(__file__).with_name("dilute_kink.py"))
DK = importlib.util.module_from_spec(DS)
DS.loader.exec_module(DK)

_MODEL = None


def _eval_one(x):
    return _MODEL.value_and_grad(x)


wrap = DK.wrap


def letters_of(t) -> str:
    """G+ / G- / T per torsion (a gauche bond is within 60 degrees of +-60)."""
    return " ".join("G+" if 0.0 < v < 120.0 else ("G-" if -120.0 < v < 0.0 else "T") for v in wrap(t))


def move_kink(t, k0: int, step: int):
    """Torsions ``t`` with the kink pair at bonds ``k0, k0 + 2`` moved ``step`` bonds along: the two
    gauche values carried to their new bonds, the bonds they leave and the bond between set
    trans.  Every other torsion keeps its value, so a copolymer chain keeps the deflections that
    belong to its sequence."""
    t = np.array(t, dtype=float)
    P = len(t)
    g1, g2 = t[k0 % P], t[(k0 + 2) % P]
    for k in (k0, k0 + 1, k0 + 2):
        t[k % P] = 180.0
    t[(k0 + step) % P], t[(k0 + step + 1) % P], t[(k0 + step + 2) % P] = g1, 180.0, g2
    return t


def free_sets(model, movers, host: str):
    """``(free indices, chains whose closure is held)`` at a fixed cell for the chosen host freedom.
    One chain's lateral position and axial shift stay fixed (the last chain that does not move, or
    the last chain when all do), which removes the cell's rigid translation."""
    nc = model.nc
    movers = [int(m) for m in movers]
    others = [s for s in range(nc) if s not in movers]
    anchor = others[-1] if others else movers[-1]
    place = lambda chains: np.concatenate([model.index("place", [s for s in chains if s != anchor]),  # noqa: E731
                                           model.index("phi", [s for s in chains if s == anchor])]).astype(int)
    m_int = model.index("internal", movers)
    if host == "rigid" and others:
        return np.concatenate([m_int, place(movers)]).astype(int), movers
    if host == "bodies" and others:
        return np.concatenate([m_int, place(movers), place(others)]).astype(int), movers
    o_int = model.index("internal", others) if others else np.zeros(0, dtype=int)
    o_place = place(others) if others else np.zeros(0, dtype=int)
    return np.concatenate([m_int, place(movers), o_place, o_int]).astype(int), list(range(nc))


def band(model, xa, xb, n_images: int, free, closed):
    """The straight path from ``xa`` to ``xb`` (each torsion and setting angle by its shorter way
    round), every image returned to the closure surface."""
    d = xb - xa
    per = np.concatenate([model.index("tors"), model.index("phi")])
    d[per] = wrap(d[per])
    imgs = np.array([xa + d * t for t in np.linspace(0.0, 1.0, n_images)])
    internal = np.intersect1d(free, model.index("internal", closed))
    for i in range(1, n_images - 1):
        imgs[i] = model.project(imgs[i], internal, closed)
    return imgs


def main():
    global _MODEL
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--row", default="vdf10")
    ap.add_argument("--nx", type=int, default=1)
    ap.add_argument("--ny", type=int, default=1)
    ap.add_argument("--mover", default="0", help="the chain whose kink pair moves, or 'all': every chain's moves together, the sheet in register "
                    "(the step without the cost of a chain leaving its neighbours; energies are then per chain)")
    ap.add_argument("--kink-at", type=int, default=0, help="bond offset of every chain's kink pair in the kinked crystal")
    ap.add_argument("--mirror-kink", action="store_true", help="G- T G+ instead of G+ T G-")
    ap.add_argument("--backward", action="store_true", help="move the kink pair one monomer back instead of ahead")
    ap.add_argument("--host", default="flexible", choices=("rigid", "bodies", "flexible"))
    ap.add_argument("--images", type=int, default=9, help="images per segment, end points included")
    ap.add_argument("--fmax", type=float, default=0.1, help="band convergence, kcal/mol per scaled unit (10 degrees of torsion, 1 A)")
    ap.add_argument("--maxiter", type=int, default=400)
    ap.add_argument("--k-spring", type=float, default=1.0)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--preset", default="pvdf-dft-valence-belt4")
    ap.add_argument("--angle-window", type=float, default=8.0)
    ap.add_argument("--references", default="docs/predictions/dilute_kink_belt4.json", help="the all-trans reference cells (shared with dilute_kink.py)")
    ap.add_argument("--out", default="docs/predictions/kink_migration_belt4.json")
    args = ap.parse_args()
    log = lambda s: print(s, flush=True)  # noqa: E731
    out_path, ref_path = Path(args.out), Path(args.references)
    letters, note = DK.sequence(args.row)
    key = (f"{args.row}@{args.kink_at}{'~' if args.mirror_kink else ''} {args.nx}x{args.ny} mover {args.mover} "
           f"{'back' if args.backward else 'ahead'} host {args.host}")
    entry: dict = {}
    prior = (json.load(open(out_path)) if out_path.exists() else {}).get("rows", {}).get(key, {})  # a stored run of this row: its end points and bands are reused

    def locked_update(path: Path, fn):
        import fcntl
        with open(str(path) + ".lock", "w") as lock:
            fcntl.flock(lock, fcntl.LOCK_EX)
            disk = json.load(open(path)) if path.exists() else {}
            fn(disk)
            with open(path, "w") as fh:
                json.dump(disk, fh, indent=1)

    def save():
        def put(disk):
            disk.update(preset=args.preset, energy="kcal/mol per chain whose kink pair moves, against the kinked crystal with every kink pair in register; "
                        "cell fixed at the kinked crystal's, every chain periodic exactly")
            disk.setdefault("rows", {})[key] = entry
        locked_update(out_path, put)

    with PRESET_PARAMETERS[args.preset].applied():
        polymer = DK.C.make_polymer(args.row, letters)
        stored = (json.load(open(ref_path)) if ref_path.exists() else {}).get("all_trans_references", {}).get(args.row)
        d_tt = DK.all_trans_reference(polymer, args.row, stored)
        if stored is None or "state" not in stored:
            locked_update(ref_path, lambda disk: disk.setdefault("all_trans_references", {}).__setitem__(args.row, DK.as_json(d_tt)))
        kkey = f"{args.row}@{args.kink_at}{'~' if args.mirror_kink else ''}"
        stored = (json.load(open(out_path)) if out_path.exists() else {}).get("kinked_references", {}).get(kkey)
        P = polymer.bonds_per_repeat
        g = 2 if args.mirror_kink else 1
        seq = [int(v) for v in np.roll(np.array([g, 0, 3 - g] + [0] * (P - 3)), args.kink_at % P)]
        d = DK.refined_cell(polymer, args.row, seq, f"{args.row} kink@{args.kink_at}{'~' if args.mirror_kink else ''}", stored)
        if stored is None or "state" not in stored:
            locked_update(out_path, lambda disk: disk.setdefault("kinked_references", {}).__setitem__(kkey, DK.as_json(d)))
        model, x0, ch = DK.make_model(polymer, d, args.nx, args.ny)
        _MODEL = model
        nc, n_mon, P = model.nc, ch.n_monomers, model.P
        movers = list(range(nc)) if args.mover == "all" else [int(args.mover)]
        nm, m0 = len(movers), movers[0]
        e_tt = float(d_tt["energy"]) * nc / 2  # the polished all-trans crystal, per supercell
        entry.update(row=args.row, letters=letters, note=note, supercell=[args.nx, args.ny], chains=nc, sites=int(model.pk.N), monomers_per_chain=int(n_mon),
                     mover=args.mover, kink_offset_bonds=args.kink_at, kink_mirrored=bool(args.mirror_kink), direction="back" if args.backward else "ahead",
                     host=args.host, images_per_segment=args.images, fmax=args.fmax)
        search = (float(d["search_energy"]) - float(d_tt["search_energy"])) / 2
        log(f"{key}: {nc} chains of {n_mon} monomers, {model.pk.N} sites; kink pair in the search's cells {search:+.3f} kcal/mol "
            f"(kinked cell's closure residual {d['search_closure_residual'][0]:.3f} deg), polished {(float(d['energy']) - float(d_tt['energy'])) / 2:+.3f}")

        # --- A: the kinked crystal, every chain's kink pair at the same offset ----------------------
        t0 = time.time()
        everyone = list(range(nc))
        xA, eA, info = x0, model.physical_energy(x0), {"iterations": 0, "message": "the polished two-chain cell, tiled"}
        if nc > 2:  # in the supercell the chains may differ: is the tiled cell still the minimum?
            cell_free = np.setdiff1d(np.arange(model.n_var), np.concatenate([model.index("xy", nc - 1), model.index("z", nc - 1)]))
            xp = x0.copy()
            xp[model.index("tors")] += np.random.default_rng(1).uniform(-0.5, 0.5, nc * P)
            xp = model.project(xp, model.index("internal"), everyone)
            xr, info = model.relax(xp, free=cell_free, closed=everyone, maxiter=600, angle_window=args.angle_window)
            if model.physical_energy(xr) < eA - 1e-4:
                xA, eA = xr, model.physical_energy(xr)
        stA = model.state(xA)
        form = (eA - e_tt) / nc
        entry["kinked_crystal"] = {"formation_per_kink_pair": form, "formation_per_kink_pair_search": search,
                                   "formation_per_kink_pair_two_chain_cell": (float(d["energy"]) - float(d_tt["energy"])) / 2, "cell": stA.lattice.tolist(),
                                   "all_trans_cell": (np.asarray(d_tt["state"][:4]) * np.array([args.nx, args.ny, 1.0, 1.0])).tolist(),
                                   "torsions": [letters_of(stA.tors[s]) for s in range(nc)], "iterations": info["iterations"], "message": info["message"],
                                   "seconds": time.time() - t0, "state": xA.tolist()}
        log(f"   kinked crystal: {form:+.3f} kcal/mol per kink pair against all-trans; cell {np.round(stA.lattice, 3).tolist()} "
            f"(all-trans {np.round(entry['kinked_crystal']['all_trans_cell'], 3).tolist()}); {info['iterations']} it, {time.time() - t0:.0f} s; {info['message']}")
        log(f"      chain {m0}: {letters_of(stA.tors[m0])}")
        save()

        # --- I and B: the mover's kink pair widened, then one monomer on ---------------------------
        free, closed = free_sets(model, movers, args.host)
        k0, step = args.kink_at % P, (-2 if args.backward else 2)
        lead = (k0 + 2) % P if not args.backward else k0  # the gauche bond that moves first
        xI0, xB0 = xA.copy(), xA.copy()
        for m in movers:
            m_t = model.index("tors", m)
            tA = xA[m_t].copy()
            tI = tA.copy()
            tI[(lead + step) % P], tI[lead] = tA[lead], 180.0
            xI0[m_t] = tI
            xB0[m_t] = move_kink(tA, k0, step)
        ends = {"A": xA}
        for name, x0 in (("I", xI0), ("B", xB0)):
            t0 = time.time()
            old = prior.get(name, {})
            if "state" in old and old.get("converged") and len(old["state"]) == model.n_var and abs(prior.get("kinked_crystal", {}).get("formation_per_kink_pair", 1e9) - form) < 1e-6:
                xr, rec = np.asarray(old["state"], dtype=float), old
                log(f"   {name}: {rec['energy']:+.3f} kcal/mol above the kinked crystal (from the stored run); {rec['torsions']}")
            else:
                x0 = model.project(x0, model.index("internal", movers), movers)
                xr, info = model.relax(x0, free=free, closed=closed, maxiter=600, angle_window=args.angle_window)
                for _ in range(3):  # a strained start can take more than one round
                    if "success" in info["message"].lower():
                        break
                    xr, more = model.relax(xr, free=free, closed=closed, maxiter=600, angle_window=args.angle_window)
                    info = {**more, "iterations": info["iterations"] + more["iterations"]}
                st = model.state(xr)
                rec = {"energy": (model.physical_energy(xr) - eA) / nm, "torsions": letters_of(st.tors[m0]), "torsions_deg": wrap(st.tors[m0]).tolist(),
                       "start_energy": (model.physical_energy(x0) - eA) / nm, "iterations": info["iterations"], "message": info["message"],
                       "converged": "success" in info["message"].lower(),
                       "angles_at_window": int((np.abs(st.angs - model.angles0[None]) > args.angle_window - 1e-3).sum()),
                       "seconds": time.time() - t0, "state": xr.tolist()}
                log(f"   {name}: {rec['energy']:+.3f} kcal/mol above the kinked crystal (start {rec['start_energy']:+.3f}); {letters_of(st.tors[m0])}; "
                    f"{info['iterations']} it, {time.time() - t0:.0f} s; {rec['angles_at_window']} angles at the window; {info['message']}")
            ends[name] = xr
            entry[name] = rec
            save()

        # --- the two bands ----------------------------------------------------------------------------
        pool = mp.get_context("fork").Pool(args.workers) if args.workers > 1 else None

        def evaluate(X):
            res = pool.map(_eval_one, list(X)) if pool is not None else [_eval_one(x) for x in X]
            return np.array([r[0] for r in res]), np.array([r[1] for r in res])

        constraint = lambda x: model.closure(x, closed)  # noqa: E731
        sc = model.scale(xA)
        top = eA
        for seg, (a, b) in (("A-I", ("A", "I")), ("I-B", ("I", "B"))):
            t0 = time.time()
            old = prior.get(seg, {})
            same_ends = "images" in old and len(old["images"]) == args.images and np.allclose(old["images"][0], ends[a], atol=1e-9) and np.allclose(old["images"][-1], ends[b], atol=1e-9)
            if same_ends and old.get("converged"):
                entry[seg] = old
                top = max(top, eA + old["saddle"] * nm)
                log(f"   {seg}: saddle {old['saddle']:+.3f} (from the stored run)")
                continue
            imgs = np.asarray(old["images"], dtype=float) if same_ends else band(model, ends[a], ends[b], args.images, free, closed)

            def progress(it, E, fm, climbing, seg=seg, t0=t0):
                if it == 1 or it % 10 == 0:
                    log(f"      {seg} it {it:4d}: force {fm:8.4f}, top {(E.max() - eA) / nm:+8.4f}, climbing {climbing}, energies {np.round((E - eA) / nm, 3).tolist()}  ({time.time() - t0:.0f} s)")

            res = MC.neb(evaluate, imgs, scale=sc, free=free, k_spring=args.k_spring, climb_after=40, climb_force=1.0, fmax=args.fmax, maxiter=args.maxiter,
                         constraint=constraint, callback=progress)
            i_top = int(np.argmax(res.energies))
            st = model.state(res.images[i_top])
            dips = [i for i in range(1, args.images - 1) if res.energies[i] < min(res.energies[i - 1], res.energies[i + 1]) - 1e-3]
            entry[seg] = {"energies": ((res.energies - eA) / nm).tolist(), "saddle": float(res.energies.max() - eA) / nm, "saddle_image": i_top, "converged": bool(res.converged),
                          "force_max": res.force_max, "iterations": res.iterations, "interior_minima": dips,
                          "saddle_torsions_deg": wrap(st.tors[m0]).tolist(), "seconds": time.time() - t0, "images": res.images.tolist()}
            top = max(top, float(res.energies.max()))
            log(f"   {seg}: saddle {(res.energies.max() - eA) / nm:+.3f} at image {i_top}{'' if res.converged else ' (NOT converged)'}; force {res.force_max:.4f}; "
                f"{res.iterations} it, {time.time() - t0:.0f} s{'; interior minima at ' + str(dips) if dips else ''}")
            save()
        if pool is not None:
            pool.close()
        entry["barrier_forward"] = (top - eA) / nm
        entry["barrier_back"] = (top - model.physical_energy(ends["B"])) / nm
        entry["moving_chains"] = nm
        log(f"   one-monomer step: barrier {entry['barrier_forward']:+.3f} ahead, {entry['barrier_back']:+.3f} back; B lies {entry['B']['energy']:+.3f} above A")
        save()
    log(f"written {out_path}")


if __name__ == "__main__":
    main()
