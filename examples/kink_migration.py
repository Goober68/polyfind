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
    ap.add_argument("--ends-only", action="store_true", help="relax the end points A, I, B and stop: no bands")
    ap.add_argument("--trust-tiled", action="store_true",
                    help="take the tiled polished kinked cell as A without relaxing the supercell with every chain free (in the four-chain PVDF cell that check returns the tiled cell)")
    ap.add_argument("--only", default=None, metavar="PIECES",
                    help="comma-separated pieces to compute in this process, of A (the supercell check of the kinked crystal), I, B, A-I, I-B; the others are "
                    "read from the stored row.  Each process writes the fields it computed, so the pieces of one row can run side by side")
    ap.add_argument("--seed-from", default=None, metavar="ROW",
                    help="start I and B from the stored I and B of another row of the same supercell (the rigid-chain-neighbour row, say) "
                    "instead of from the kinked crystal with the mover's torsions reset")
    ap.add_argument("--band-from", default=None, metavar="ROW",
                    help="start each band from the converged band of another row (a smaller supercell, say): the mover's torsions and backbone angles are "
                    "taken from that row's images, everything else runs straight between this row's end points")
    ap.add_argument("--end-maxiter", type=int, default=600, help="iterations of one relaxation round of an end point (three rounds at most)")
    ap.add_argument("--patience", type=int, default=None,
                    help="stop a relaxation when the lowest energy has not fallen by --patience-tol over this many evaluations (the optimiser's own test, an energy change "
                    "under 1e-9 kcal/mol, is at the rounding of a 1600-site cell's energy and is not met there); the result counts as converged if its gradient along the "
                    "closure surface is under --end-gradient")
    ap.add_argument("--patience-tol", type=float, default=1e-5, help="kcal/mol per supercell")
    ap.add_argument("--end-gradient", type=float, default=0.01, help="kcal/mol per scaled unit (a tenth of the band's default convergence)")
    ap.add_argument("--checkpoints", default=None, metavar="DIR",
                    help="write the running state of every relaxation and band into DIR every ten evaluations, and resume from it when it is there")
    ap.add_argument("--references", default="docs/predictions/dilute_kink_belt4.json", help="the all-trans reference cells (shared with dilute_kink.py)")
    ap.add_argument("--out", default="docs/predictions/kink_migration_belt4.json")
    args = ap.parse_args()
    log = lambda s: print(s, flush=True)  # noqa: E731
    out_path, ref_path = Path(args.out), Path(args.references)
    letters, note = DK.sequence(args.row)
    key = (f"{args.row}@{args.kink_at}{'~' if args.mirror_kink else ''} {args.nx}x{args.ny} mover {args.mover} "
           f"{'back' if args.backward else 'ahead'} host {args.host}")
    entry: dict = {}
    mine: set = set()  # the fields this process computed and has not written yet
    stored_rows = (json.load(open(out_path)) if out_path.exists() else {}).get("rows", {})
    prior = stored_rows.get(key, {})  # a stored run of this row: its end points and bands are reused
    pieces = ("A", "I", "B", "A-I", "I-B")
    only = None if args.only is None else {v.strip() for v in args.only.split(",")}
    if only is not None and not only <= set(pieces):
        raise SystemExit(f"--only takes {', '.join(pieces)}; got {args.only}")
    do = lambda piece: only is None or piece in only  # noqa: E731
    ck_dir = None if args.checkpoints is None else Path(args.checkpoints)
    if ck_dir is not None:
        ck_dir.mkdir(parents=True, exist_ok=True)
    slug = "".join(ch if ch.isalnum() else "_" for ch in key)

    def checkpoint(piece: str):
        return None if ck_dir is None else ck_dir / f"{slug}.{piece}.npy"

    def write_checkpoint(path, arr):
        tmp = path.with_suffix(".tmp.npy")
        np.save(tmp, np.asarray(arr, dtype=float))
        tmp.replace(path)

    def locked_update(path: Path, fn):
        import fcntl
        with open(str(path) + ".lock", "w") as lock:
            fcntl.flock(lock, fcntl.LOCK_EX)
            disk = json.load(open(path)) if path.exists() else {}
            fn(disk)
            with open(path, "w") as fh:
                json.dump(disk, fh, indent=1)

    def save():
        """Merge what this process computed into the stored row (field by field: another process may hold the other pieces of the row)."""
        def put(disk):
            disk.update(preset=args.preset, energy="kcal/mol per chain whose kink pair moves, against the kinked crystal with every kink pair in register; "
                        "cell fixed at the kinked crystal's, every chain periodic exactly")
            row = disk.setdefault("rows", {}).setdefault(key, {})
            row.update({k: v for k, v in entry.items() if k in mine or k not in pieces_fields})
            mine.clear()
        locked_update(out_path, put)

    pieces_fields = ("kinked_crystal", "I", "B", "A-I", "I-B")

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

        def watch(piece, e_ref, per, t_start):
            """A progress hook for model.relax: every ten evaluations a log line, and the lowest point
            seen so far as the checkpoint (the latest evaluation may be a trial step of a line search)."""
            path = checkpoint(piece)
            best = {"e": np.inf, "x": None, "n": 0}

            def hook(stage, k, e, x):
                best["n"] += 1
                if e < best["e"]:
                    best.update(e=e, x=x.copy())
                if k == 1 or k % 10 == 0:
                    log(f"      {piece} {stage} {k:5d}: {(e - e_ref) / per:+10.4f}, lowest so far {(best['e'] - e_ref) / per:+10.4f}  "
                        f"({best['n']} evaluations, {time.time() - t_start:.0f} s)")
                    if path is not None:
                        write_checkpoint(path, best["x"])
            return hook

        def resumed(piece, x_start, free_, closed_):
            """The checkpoint of ``piece`` if there is one, returned to the closure surface; else ``x_start``."""
            ck = checkpoint(piece)
            if ck is None or not ck.exists() or np.load(ck).shape != (model.n_var,):
                return x_start
            log(f"   {piece}: resumed from {ck}")
            return model.project(np.load(ck), np.intersect1d(free_, model.index("internal", closed_)), closed_)

        stop = dict(patience=args.patience, patience_tol=args.patience_tol)
        ended = lambda info: "success" in info["message"].lower() or (info["stalled"] and info["gradient_max"] < args.end_gradient)  # noqa: E731

        old = prior.get("kinked_crystal", {})
        check = nc > 2 and not args.trust_tiled and do("A")  # in the supercell the chains may differ: is the tiled cell still the minimum?
        if not check and "state" in old and len(old["state"]) == model.n_var:
            xA = np.asarray(old["state"], dtype=float)
            eA = model.physical_energy(xA)
            entry["kinked_crystal"] = old
            stA = model.state(xA)
            form = (eA - e_tt) / nc
            log(f"   kinked crystal: {form:+.3f} kcal/mol per kink pair against all-trans (from the stored run: {old.get('message')})")
        else:
            xA, eA, info = x0, model.physical_energy(x0), {"iterations": 0, "message": "the polished two-chain cell, tiled"}
            checked = None
            if check:
                cell_free = np.setdiff1d(np.arange(model.n_var), np.concatenate([model.index("xy", nc - 1), model.index("z", nc - 1)]))
                xp = x0.copy()
                xp[model.index("tors")] += np.random.default_rng(1).uniform(-0.5, 0.5, nc * P)
                xp = model.project(xp, model.index("internal"), everyone)
                e_start = model.physical_energy(xp)
                xr, info = model.relax(resumed("A", xp, cell_free, everyone), free=cell_free, closed=everyone, maxiter=args.end_maxiter,
                                       angle_window=args.angle_window, progress=watch("A", eA, nc, t0), **stop)
                e_r = model.physical_energy(xr)
                checked = {"start_above_tiled_per_chain": (e_start - eA) / nc, "relaxed_above_tiled_per_chain": (e_r - eA) / nc,
                           "cell_relaxed": model.state(xr).lattice.tolist(), "converged": bool(ended(info)), "gradient_max": info["gradient_max"],
                           "largest_torsion_change_deg": float(np.abs(wrap(model.state(xr).tors - model.state(x0).tors)).max()),
                           "state_relaxed": xr.tolist()}
                log(f"   A check: the perturbed supercell, every chain and the cell free, relaxed to {(e_r - eA) / nc:+.6f} kcal/mol per chain against the tiled cell "
                    f"(start {(e_start - eA) / nc:+.4f}); largest torsion change {checked['largest_torsion_change_deg']:.3f} deg; {info['message']}")
                if e_r < eA - 1e-4:
                    xA, eA = xr, e_r
            stA = model.state(xA)
            form = (eA - e_tt) / nc
            entry["kinked_crystal"] = {"formation_per_kink_pair": form, "formation_per_kink_pair_search": search,
                                       "formation_per_kink_pair_two_chain_cell": (float(d["energy"]) - float(d_tt["energy"])) / 2, "cell": stA.lattice.tolist(),
                                       "all_trans_cell": (np.asarray(d_tt["state"][:4]) * np.array([args.nx, args.ny, 1.0, 1.0])).tolist(),
                                       "torsions": [letters_of(stA.tors[s]) for s in range(nc)], "iterations": info["iterations"], "message": info["message"],
                                       "seconds": time.time() - t0, "state": xA.tolist()}
            if checked is not None:
                entry["kinked_crystal"]["supercell_check"] = checked
            mine.add("kinked_crystal")
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
        seed = None
        if args.seed_from is not None:
            seed = stored_rows.get(args.seed_from)
            if seed is None or abs(seed.get("kinked_crystal", {}).get("formation_per_kink_pair", 1e9) - form) > 1e-6:
                raise SystemExit(f"--seed-from: no stored row {args.seed_from!r} on this kinked crystal")
        for name, x0 in (("I", xI0), ("B", xB0)):
            t0 = time.time()
            old = prior.get(name, {})
            if "state" in old and old.get("converged") and len(old["state"]) == model.n_var and abs(prior.get("kinked_crystal", {}).get("formation_per_kink_pair", 1e9) - form) < 1e-6:
                xr, rec = np.asarray(old["state"], dtype=float), old
                log(f"   {name}: {rec['energy']:+.3f} kcal/mol above the kinked crystal (from the stored run); {rec['torsions']}")
            elif not do(name):
                continue
            else:
                start = "the kinked crystal with the mover's torsions reset"
                if seed is not None:
                    if len(seed.get(name, {}).get("state", [])) != model.n_var:
                        raise SystemExit(f"--seed-from: row {args.seed_from!r} holds no {name} of this supercell")
                    x0, start = np.asarray(seed[name]["state"], dtype=float), f"{name} of the row {args.seed_from!r}"
                x0 = model.project(x0, model.index("internal", movers), movers)
                e_start = model.physical_energy(x0)
                x0 = resumed(name, x0, free, closed)
                hook = watch(name, eA, nm, t0)
                xr, info = model.relax(x0, free=free, closed=closed, maxiter=args.end_maxiter, angle_window=args.angle_window, progress=hook, **stop)
                for _ in range(3):  # a strained start can take more than one round
                    if ended(info):
                        break
                    xr, more = model.relax(xr, free=free, closed=closed, maxiter=args.end_maxiter, angle_window=args.angle_window, progress=hook, **stop)
                    info = {**more, "iterations": info["iterations"] + more["iterations"], "evaluations": info["evaluations"] + more["evaluations"]}
                st = model.state(xr)
                rec = {"energy": (model.physical_energy(xr) - eA) / nm, "torsions": letters_of(st.tors[m0]), "torsions_deg": wrap(st.tors[m0]).tolist(),
                       "start": start, "start_energy": (e_start - eA) / nm, "iterations": info["iterations"], "evaluations": info["evaluations"],
                       "message": info["message"], "converged": bool(ended(info)), "gradient_max": info["gradient_max"],
                       "angles_at_window": int((np.abs(st.angs - model.angles0[None]) > args.angle_window - 1e-3).sum()),
                       "seconds": time.time() - t0, "state": xr.tolist()}
                mine.add(name)
                log(f"   {name}: {rec['energy']:+.3f} kcal/mol above the kinked crystal (start {rec['start_energy']:+.3f}); {letters_of(st.tors[m0])}; "
                    f"{info['iterations']} it, {info['evaluations']} evaluations, {time.time() - t0:.0f} s; gradient {info['gradient_max']:.1e}; "
                    f"{rec['angles_at_window']} angles at the window; {info['message']}")
            ends[name] = xr
            entry[name] = rec
            save()

        if args.ends_only or (only is not None and not (only & {"A-I", "I-B"})):
            log(f"   end points only; written {out_path}")
            return

        # --- the two bands ----------------------------------------------------------------------------
        pool = None  # made when a band is run: the workers are forks of this process and its model

        def evaluate(X):
            res = pool.map(_eval_one, list(X)) if pool is not None else [_eval_one(x) for x in X]
            return np.array([r[0] for r in res]), np.array([r[1] for r in res])

        constraint = lambda x: model.closure(x, closed)  # noqa: E731
        sc = model.scale(xA)
        per = np.concatenate([model.index("tors"), model.index("phi")])

        def same(u, v):  # the same structure: equal up to whole turns of the torsions and setting angles
            dv = np.asarray(u, dtype=float) - np.asarray(v, dtype=float)
            dv[per] = wrap(dv[per])
            return bool(np.abs(dv).max() < 1e-8)

        guide = None
        if args.band_from is not None:
            guide = stored_rows.get(args.band_from)
            if guide is None:
                raise SystemExit(f"--band-from: no stored row {args.band_from!r}")

        def guided(seg, a, b):
            """The straight band between this row's end points with the movers' internal coordinates
            taken from the images of the guide row's band (the same chain in another supercell)."""
            g_imgs = np.asarray(guide[seg]["images"], dtype=float)
            g_nc = int(guide["chains"])
            g_base, width = 4 + 4 * g_nc, model.P + model.R
            if len(g_imgs) != args.images or g_imgs.shape[1] != g_base + g_nc * width or str(guide["mover"]) != str(args.mover) or len(movers) != 1:
                raise SystemExit(f"--band-from: the band {seg} of {args.band_from!r} does not fit this row (images, chain length or mover)")
            imgs = band(model, ends[a], ends[b], args.images, free, closed)
            m_int = model.index("internal", m0)
            g_int = np.arange(g_base + m0 * width, g_base + (m0 + 1) * width)
            d0 = wrap(g_imgs[0, g_int] - ends[a][m_int])  # the guide's end points are not this row's: carry the difference along the band
            d1 = wrap(g_imgs[-1, g_int] - ends[b][m_int])
            t_tors = np.arange(width) < model.P
            for i in range(1, args.images - 1):
                w = i / (args.images - 1.0)
                v = g_imgs[i, g_int] - ((1.0 - w) * d0 + w * d1)
                ref = imgs[i][m_int]
                imgs[i][m_int] = np.where(t_tors, ref + wrap(v - ref), v)
                imgs[i] = model.project(imgs[i], np.intersect1d(free, model.index("internal", closed)), closed)
            return imgs

        for seg, (a, b) in (("A-I", ("A", "I")), ("I-B", ("I", "B"))):
            t0 = time.time()
            old = prior.get(seg, {})
            if a not in ends or b not in ends:
                if do(seg):
                    log(f"   {seg}: the end points {a} and {b} are not both stored for this row; nothing done")
                continue
            same_ends = "images" in old and len(old["images"]) == args.images and same(old["images"][0], ends[a]) and same(old["images"][-1], ends[b])
            if same_ends and old.get("converged"):
                entry[seg] = old
                log(f"   {seg}: saddle {old['saddle']:+.3f} (from the stored run)")
                continue
            if not do(seg):
                continue
            ck = checkpoint(seg)
            start = "the stored band, continued" if same_ends else ("straight between the end points" if guide is None else f"the movers' internal coordinates from the band of {args.band_from!r}")
            if ck is not None and ck.exists() and np.load(ck).shape == (args.images, model.n_var) and same(np.load(ck)[0], ends[a]) and same(np.load(ck)[-1], ends[b]):
                imgs, start = np.load(ck), start + f"; resumed from {ck}"
            elif same_ends:
                imgs = np.asarray(old["images"], dtype=float)
            else:
                imgs = band(model, ends[a], ends[b], args.images, free, closed) if guide is None else guided(seg, a, b)
            log(f"   {seg}: start: {start}")
            if pool is None and args.workers > 1:
                pool = mp.get_context("fork").Pool(args.workers)

            def progress(it, E, fm, climbing, seg=seg, t0=t0):
                if it == 1 or it % 10 == 0:
                    log(f"      {seg} it {it:4d}: force {fm:8.4f}, top {(E.max() - eA) / nm:+8.4f}, climbing {climbing}, energies {np.round((E - eA) / nm, 3).tolist()}  ({time.time() - t0:.0f} s)")

            def snapshot(it, X, ck=ck):
                if ck is not None and it % 10 == 0:
                    write_checkpoint(ck, X)

            res = MC.neb(evaluate, imgs, scale=sc, free=free, k_spring=args.k_spring, climb_after=40, climb_force=1.0, fmax=args.fmax, maxiter=args.maxiter,
                         constraint=constraint, callback=progress, snapshot=snapshot)
            if ck is not None:
                write_checkpoint(ck, res.images)
            i_top = int(np.argmax(res.energies))
            st = model.state(res.images[i_top])
            dips = [i for i in range(1, args.images - 1) if res.energies[i] < min(res.energies[i - 1], res.energies[i + 1]) - 1e-3]
            entry[seg] = {"energies": ((res.energies - eA) / nm).tolist(), "saddle": float(res.energies.max() - eA) / nm, "saddle_image": i_top, "converged": bool(res.converged),
                          "force_max": res.force_max, "fmax": args.fmax, "iterations": res.iterations + (old.get("iterations", 0) if same_ends else 0), "interior_minima": dips,
                          "start": start, "saddle_torsions_deg": wrap(st.tors[m0]).tolist(), "seconds": time.time() - t0, "images": res.images.tolist()}
            mine.add(seg)
            log(f"   {seg}: saddle {(res.energies.max() - eA) / nm:+.3f} at image {i_top}{'' if res.converged else ' (NOT converged)'}; force {res.force_max:.4f}; "
                f"{res.iterations} it, {time.time() - t0:.0f} s{'; interior minima at ' + str(dips) if dips else ''}")
            save()
        if pool is not None:
            pool.close()

        def summarise(disk):  # the barriers, from the row as it now stands on disk (another process may have written the other band)
            row = disk.get("rows", {}).get(key, {})
            if all(k in row for k in ("A-I", "I-B", "B")):
                top = max(0.0, float(row["A-I"]["saddle"]), float(row["I-B"]["saddle"]))
                row.update(barrier_forward=top, barrier_back=top - float(row["B"]["energy"]), moving_chains=nm,
                           bands_converged=bool(row["A-I"].get("converged") and row["I-B"].get("converged")))
                log(f"   one-monomer step: barrier {top:+.3f} ahead, {top - float(row['B']['energy']):+.3f} back; B lies {float(row['B']['energy']):+.3f} above A"
                    f"{'' if row['bands_converged'] else ' (a band is NOT converged)'}")
        locked_update(out_path, summarise)
    log(f"written {out_path}")


if __name__ == "__main__":
    main()
