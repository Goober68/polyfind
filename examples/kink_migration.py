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
surface with a climbing image.  Reported per row: the energies of I and B above A (the price of
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


def wrap(a):
    return (np.asarray(a, dtype=float) + 180.0) % 360.0 - 180.0


def letters_of(t) -> str:
    """G+ / G- / T per torsion (a gauche bond is within 60 degrees of +-60)."""
    return " ".join("G+" if 0.0 < v < 120.0 else ("G-" if -120.0 < v < 0.0 else "T") for v in wrap(t))


def shifted(v, k: int):
    """A per-bond array moved ``k`` bonds along the chain (cyclically)."""
    return np.roll(np.asarray(v, dtype=float), k)


def free_sets(model, mover: int, host: str):
    """``(free indices, chains whose closure is held)`` at a fixed cell for the chosen host freedom."""
    nc = model.nc
    others = [s for s in range(nc) if s != mover]
    anchor = others[-1]
    place = lambda chains: np.concatenate([model.index("place", [s for s in chains if s != anchor]),  # noqa: E731
                                           model.index("phi", [s for s in chains if s == anchor])]).astype(int)
    m_int = model.index("internal", mover)
    if host == "rigid":
        return np.concatenate([m_int, model.index("place", mover)]).astype(int), [mover]
    if host == "bodies":
        return np.concatenate([m_int, model.index("place", mover), place(others)]).astype(int), [mover]
    return np.concatenate([m_int, model.index("place", mover), place(others), model.index("internal", others)]).astype(int), list(range(nc))


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
    ap.add_argument("--mover", type=int, default=0, help="the chain whose kink pair moves")
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
        d = DK.all_trans_reference(polymer, args.row, stored)
        if stored is None:
            locked_update(ref_path, lambda disk: disk.setdefault("all_trans_references", {}).__setitem__(
                args.row, {k: (v.tolist() if isinstance(v, np.ndarray) else v) for k, v in d.items()}))
        model, x_tt, tt = DK.make_model(polymer, d, args.nx, args.ny)
        _MODEL = model
        nc, n_mon, P = model.nc, tt.n_monomers, model.P
        e_tt = model.physical_energy(x_tt)
        entry.update(row=args.row, letters=letters, note=note, supercell=[args.nx, args.ny], chains=nc, sites=int(model.pk.N), monomers_per_chain=int(n_mon),
                     mover=args.mover, kink_offset_bonds=args.kink_at, kink_mirrored=bool(args.mirror_kink), direction="back" if args.backward else "ahead",
                     host=args.host, images_per_segment=args.images, fmax=args.fmax)
        log(f"{key}: {nc} chains of {n_mon} monomers, {model.pk.N} sites")

        # --- A: the kinked crystal, every chain's kink pair at the same offset, cell free -----------
        t0 = time.time()
        x = x_tt.copy()
        for s in range(nc):
            x[model.index("tors", s)] = DK.kinked_torsions(d["tors"], args.kink_at, args.mirror_kink)
        x[3] = float(model.terms(x)["repeats"][0])  # the kinked chain's own repeat: no chain starts under tension
        everyone = list(range(nc))
        x = model.project(x, model.index("internal"), everyone)
        cell_free = np.setdiff1d(np.arange(model.n_var), np.concatenate([model.index("xy", nc - 1), model.index("z", nc - 1)]))
        xA, info = model.relax(x, free=cell_free, closed=everyone, maxiter=600, angle_window=args.angle_window)
        eA = model.physical_energy(xA)
        stA = model.state(xA)
        form = (eA - e_tt) / nc
        entry["kinked_crystal"] = {"formation_per_kink_pair": form, "cell": stA.lattice.tolist(), "all_trans_cell": model.state(x_tt).lattice.tolist(),
                                   "torsions": [letters_of(stA.tors[s]) for s in range(nc)], "iterations": info["iterations"], "message": info["message"],
                                   "seconds": time.time() - t0, "state": xA.tolist()}
        log(f"   kinked crystal: {form:+.3f} kcal/mol per kink pair against all-trans; cell {np.round(stA.lattice, 3).tolist()} (all-trans {np.round(model.state(x_tt).lattice, 3).tolist()}); "
            f"{info['iterations']} it, {time.time() - t0:.0f} s; {info['message']}")
        log(f"      chain {args.mover}: {letters_of(stA.tors[args.mover])}")
        save()

        # --- I and B: the mover's kink pair widened, then one monomer on ---------------------------
        free, closed = free_sets(model, args.mover, args.host)
        m_t, m_a = model.index("tors", args.mover), model.index("angs", args.mover)
        k0, step = args.kink_at % P, (-2 if args.backward else 2)
        tA, aA = xA[m_t].copy(), xA[m_a].copy()
        lead = (k0 + 2) % P if not args.backward else k0  # the gauche bond that moves first
        tI = tA.copy()
        tI[(lead + step) % P], tI[lead] = tA[lead], 180.0
        xI0, xB0 = xA.copy(), xA.copy()
        xI0[m_t] = tI
        xB0[m_t], xB0[m_a] = shifted(tA, step), shifted(aA, step)
        ends = {"A": xA}
        for name, x0 in (("I", xI0), ("B", xB0)):
            t0 = time.time()
            x0 = model.project(x0, model.index("internal", args.mover), [args.mover])
            xr, info = model.relax(x0, free=free, closed=closed, maxiter=600, angle_window=args.angle_window)
            ends[name] = xr
            st = model.state(xr)
            entry[name] = {"energy": model.physical_energy(xr) - eA, "torsions": letters_of(st.tors[args.mover]), "torsions_deg": wrap(st.tors[args.mover]).tolist(),
                           "start_energy": model.physical_energy(x0) - eA, "iterations": info["iterations"], "message": info["message"],
                           "angles_at_window": int((np.abs(st.angs - model.angles0[None]) > args.angle_window - 1e-3).sum()),
                           "seconds": time.time() - t0, "state": xr.tolist()}
            log(f"   {name}: {entry[name]['energy']:+.3f} kcal/mol above the kinked crystal (start {entry[name]['start_energy']:+.3f}); {letters_of(st.tors[args.mover])}; "
                f"{info['iterations']} it, {time.time() - t0:.0f} s; {entry[name]['angles_at_window']} angles at the window; {info['message']}")
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
            imgs = band(model, ends[a], ends[b], args.images, free, closed)

            def progress(it, E, fm, climbing, seg=seg, t0=t0):
                if it == 1 or it % 10 == 0:
                    log(f"      {seg} it {it:4d}: force {fm:8.4f}, top {E.max() - eA:+8.4f}, climbing {climbing}, energies {np.round(E - eA, 3).tolist()}  ({time.time() - t0:.0f} s)")

            res = MC.neb(evaluate, imgs, scale=sc, free=free, k_spring=args.k_spring, climb_after=40, fmax=args.fmax, maxiter=args.maxiter,
                         constraint=constraint, callback=progress)
            i_top = int(np.argmax(res.energies))
            st = model.state(res.images[i_top])
            dips = [i for i in range(1, args.images - 1) if res.energies[i] < min(res.energies[i - 1], res.energies[i + 1]) - 1e-3]
            entry[seg] = {"energies": (res.energies - eA).tolist(), "saddle": float(res.energies.max() - eA), "saddle_image": i_top, "converged": bool(res.converged),
                          "force_max": res.force_max, "iterations": res.iterations, "interior_minima": dips,
                          "saddle_torsions_deg": wrap(st.tors[args.mover]).tolist(), "seconds": time.time() - t0, "images": res.images.tolist()}
            top = max(top, float(res.energies.max()))
            log(f"   {seg}: saddle {res.energies.max() - eA:+.3f} at image {i_top}{'' if res.converged else ' (NOT converged)'}; force {res.force_max:.4f}; "
                f"{res.iterations} it, {time.time() - t0:.0f} s{'; interior minima at ' + str(dips) if dips else ''}")
            save()
        if pool is not None:
            pool.close()
        entry["barrier_forward"] = top - eA
        entry["barrier_back"] = top - model.physical_energy(ends["B"])
        log(f"   one-monomer step: barrier {entry['barrier_forward']:+.3f} ahead, {entry['barrier_back']:+.3f} back; B lies {entry['B']['energy']:+.3f} above A")
        save()
    log(f"written {out_path}")


if __name__ == "__main__":
    main()
