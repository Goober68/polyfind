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

A band finds where a pass is and does not settle it: ``--polish`` takes each stored band's top to
the stationary point (polyfind.mixedcell.saddle_polish).  In a cell of sixteen chains and more the
minimisations (the end states, the hyperplanes of ``--polish``) need ``--surface``: on the closure
surface with one curvature block per chain as the metric (MixedModel.relax_on_surface,
chain_metric), tens of evaluations where SLSQP with the closure as constraints takes hundreds
and stops short of the minimum.  ``--refine`` relaxes a stored end state further as a check.
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


def settled_barrier(row: dict, log=None) -> None:
    """Put the row's barriers on its settled saddles (``barrier_forward_polished``, ``barrier_back_polished``).

    Both saddles settled: the higher of the two.  One settled and the other not: still given when the
    settled one lies above the other band's top and that band had converged (a converged band's top has
    been at or above its saddle in every case settled so far), with ``barrier_polished_note`` saying
    which top is a band's; otherwise the fields are removed."""
    for k in ("barrier_forward_polished", "barrier_back_polished", "barrier_polished_note"):
        row.pop(k, None)
    done = {k: row[k].get("polished", {}).get("converged") for k in ("A-I", "I-B")}
    if not any(done.values()):
        return
    tops = {k: float(row[k]["polished"]["saddle"] if done[k] else row[k]["saddle"]) for k in done}
    high = max(tops, key=tops.get)
    if not done[high]:
        return  # the higher top is a band's: nothing settled about the barrier
    other = "A-I" if high == "I-B" else "I-B"
    if not done[other] and not (row[other].get("force_max", 1e9) < 0.11):
        return
    ptop = max(0.0, tops[high])
    row.update(barrier_forward_polished=ptop, barrier_back_polished=ptop - float(row["B"]["energy"]))
    if not done[other]:
        row["barrier_polished_note"] = f"the saddle {high} is settled; {other} is its band's top ({tops[other]:+.3f} at a force of {row[other]['force_max']:.2f}), which lies below it"
    if log is not None:
        log(f"   on the settled saddle{'s' if all(done.values()) else ''}: barrier {ptop:+.3f} ahead, {ptop - float(row['B']['energy']):+.3f} back"
            f"{'' if all(done.values()) else ' (' + row['barrier_polished_note'] + ')'}")


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
    ap.add_argument("--mover", default="0", help="the chain whose kink pair moves; several chains separated by commas (they move together; energies are per moving chain); "
                    "or 'all': every chain's moves together, the sheet in register (the step without the cost of a chain leaving its neighbours)")
    ap.add_argument("--kink-at", type=int, default=0, help="bond offset of every chain's kink pair in the kinked crystal")
    ap.add_argument("--mirror-kink", action="store_true", help="G- T G+ instead of G+ T G-")
    ap.add_argument("--backward", action="store_true", help="move the kink pair one monomer back instead of ahead")
    ap.add_argument("--host", default="flexible", choices=("rigid", "bodies", "flexible"))
    ap.add_argument("--images", type=int, default=9, help="images per segment, end points included")
    ap.add_argument("--fmax", type=float, default=0.1, help="band convergence, kcal/mol per scaled unit (10 degrees of torsion, 1 A)")
    ap.add_argument("--maxiter", type=int, default=400)
    ap.add_argument("--k-spring", type=float, default=1.0)
    ap.add_argument("--polish", action="store_true",
                    help="settle the top of each stored band instead of relaxing the band further: from its highest image, the energy minimised in the hyperplane normal "
                    "to the band's tangent and the hyperplane moved to where the force along the tangent vanishes (polyfind.mixedcell.saddle_polish); "
                    "the result is stored in the band's record as 'polished'")
    ap.add_argument("--polish-trust", type=float, default=None, help="confine each hyperplane's minimisation to this distance (scaled units) from its start; "
                    "one that ends on its box is continued in the same hyperplane")
    ap.add_argument("--polish-tol", type=float, default=0.01, help="the force along the tangent at which a polished saddle is accepted, kcal/mol per scaled unit")
    ap.add_argument("--surface", action="store_true",
                    help="minimise on the closure surface with one curvature block per chain as the metric (MixedModel.relax_on_surface, chain_metric) instead of "
                    "handing the closure to SLSQP: the end states I and B and the hyperplanes of --polish.  What a cell of sixteen chains needs (tens of evaluations "
                    "where SLSQP takes hundreds and stops short); the blocks come from the four-chain tiling of the same sheet, about a minute")
    ap.add_argument("--surface-gtol", type=float, default=1e-4, help="with --surface: the gradient along the surface at which a minimisation ends, kcal/mol per scaled unit")
    ap.add_argument("--refine", action="store_true",
                    help="relax the stored end states named by --only further instead of reusing them (with --surface, on the closure surface); a result is stored "
                    "only if it is lower, with the earlier energy kept beside it.  The stored bands end on the earlier states and are kept")
    ap.add_argument("--band-step", type=float, default=0.15, help="the longest move of any variable in one band iteration, scaled units")
    ap.add_argument("--climb-after", type=int, default=40, help="band iterations before the climbing image may be switched on (it also waits for the band's largest force "
                    "to fall under 1); 0 when a band that was already climbing is continued")
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
                    help="stop a relaxation when the lowest energy has not fallen by --patience-tol over this many evaluations (in a sixteen-chain cell the optimiser's own "
                    "test, an energy change and a closure residual under 1e-9 between iterations, went on being failed for hundreds of evaluations after the energy had settled to 1e-4); "
                    "the result counts as converged if its gradient along the closure surface is under --end-gradient")
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
            tmp = Path(str(path) + ".tmp")  # written whole, then moved: a reader without the lock never sees half a file
            with open(tmp, "w") as fh:
                json.dump(disk, fh, indent=1)
            tmp.replace(path)

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
        movers = list(range(nc)) if args.mover == "all" else [int(v) for v in str(args.mover).split(",")]
        if len(set(movers)) != len(movers) or min(movers) < 0 or max(movers) >= nc:
            raise SystemExit(f"--mover: chains {args.mover} in a cell of {nc}")
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
        metric = None
        if args.surface:
            t_m = time.time()
            mfile = None if ck_dir is None else ck_dir / f"{''.join(c if c.isalnum() else '_' for c in kkey)}.chain_blocks.npz"
            if mfile is not None and mfile.exists():
                blocks = {int(k): v for k, v in np.load(mfile).items()}
            else:  # one chain's block hardly depends on the cell around it: take them from the four-chain tiling
                small, xs, _ = (model, x0, ch) if (args.nx, args.ny) == (2, 1) else DK.make_model(polymer, d, 2, 1)
                blocks = MC.chain_blocks(small, xs, range(2))
                if mfile is not None:
                    np.savez(mfile, **{str(k): v for k, v in blocks.items()})
            metric = MC.chain_metric(model, xA, free, closed, blocks)
            log(f"   metric: one curvature block per chain, from the four-chain tiling of the sheet ({time.time() - t_m:.0f} s)")

        def minimise(x, hook):
            """One end state's relaxation: SLSQP with the closure as constraints, or with --surface on the closure surface."""
            if args.surface:
                return model.relax_on_surface(x, free=free, closed=closed, metric=metric, gtol=args.surface_gtol, maxiter=args.end_maxiter,
                                              presteps=300, angle_window=args.angle_window, progress=hook, **stop)
            return model.relax(x, free=free, closed=closed, maxiter=args.end_maxiter, angle_window=args.angle_window, progress=hook, **stop)

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
            stored_end = ("state" in old and old.get("converged") and len(old["state"]) == model.n_var
                          and abs(prior.get("kinked_crystal", {}).get("formation_per_kink_pair", 1e9) - form) < 1e-6)
            if stored_end and not (args.refine and only is not None and name in only):
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
                if stored_end:  # --refine
                    x0, start = np.asarray(old["state"], dtype=float), f"the stored {name} ({old['energy']:+.6f}; {old.get('message', '')[:40]}), relaxed further"
                x0 = model.project(x0, model.index("internal", movers), movers)
                e_start = model.physical_energy(x0)
                x0 = resumed(name, x0, free, closed)
                hook = watch(name, eA, nm, t0)
                xr, info = minimise(x0, hook)
                for _ in range(3):  # a strained start can take more than one round
                    if ended(info):
                        break
                    xr, more = minimise(xr, hook)
                    info = {**more, "iterations": info["iterations"] + more["iterations"], "evaluations": info["evaluations"] + more["evaluations"]}
                st = model.state(xr)
                rec = {"energy": (model.physical_energy(xr) - eA) / nm, "torsions": letters_of(st.tors[m0]), "torsions_deg": wrap(st.tors[m0]).tolist(),
                       **({"torsions_movers": {str(m): letters_of(st.tors[m]) for m in movers}} if 1 < nm < nc else {}),
                       "start": start, "start_energy": (e_start - eA) / nm, "iterations": info["iterations"], "evaluations": info["evaluations"],
                       "message": info["message"], "converged": bool(ended(info)), "gradient_max": info["gradient_max"],
                       "method": "on the closure surface, one curvature block per chain as the metric" if args.surface else "SLSQP, the closure as equality constraints",
                       "angles_at_window": int((np.abs(st.angs - model.angles0[None]) > args.angle_window - 1e-3).sum()),
                       "seconds": time.time() - t0, "state": xr.tolist()}
                log(f"   {name}: {rec['energy']:+.3f} kcal/mol above the kinked crystal (start {rec['start_energy']:+.3f}); {letters_of(st.tors[m0])}; "
                    f"{info['iterations']} it, {info['evaluations']} evaluations, {time.time() - t0:.0f} s; gradient {info['gradient_max']:.1e}; "
                    f"{rec['angles_at_window']} angles at the window; {info['message']}")
                if stored_end:
                    gain = float(old["energy"]) - rec["energy"]
                    log(f"   {name}: the stored state relaxed further falls by {gain:.6f} kcal/mol to {rec['energy']:+.6f} (gradient {old.get('gradient_max', float('nan')):.1e} before, {info['gradient_max']:.1e} now)")
                    check = {"energy": rec["energy"], "fall": gain, "gradient_max": info["gradient_max"], "evaluations": info["evaluations"], "message": info["message"], "method": rec["method"],
                             "largest_torsion_change_deg": float(np.abs(wrap(st.tors - model.state(np.asarray(old["state"], dtype=float)).tors)).max()), "seconds": time.time() - t0}
                    xr, rec = np.asarray(old["state"], dtype=float), {**old, "relaxed_further": check}  # the bands end on the stored state: it stays, with the check beside it
                mine.add(name)
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
            if args.polish:
                if not do(seg):
                    continue
                if not same_ends:
                    log(f"   {seg}: no stored band between this row's end points to polish")
                    continue
                imgs, i_top = np.asarray(old["images"], dtype=float), int(old["saddle_image"])
                how = dict(on_surface=True, metric=metric, gtol=args.surface_gtol) if args.surface else ({} if args.polish_trust is None else {"trust": args.polish_trust})
                xs, info = MC.saddle_polish(model, imgs[i_top], imgs[i_top + 1] - imgs[i_top - 1], free, closed, tol=args.polish_tol,
                                            log=lambda text, seg=seg, t0=t0: log(f"      {seg} polish, {text}  ({time.time() - t0:.0f} s)"),
                                            patience=args.patience if args.patience is not None else 40, patience_tol=args.patience_tol, angle_window=args.angle_window,
                                            progress=watch(f"{seg} polish", eA, nm, t0), **how)
                st = model.state(xs)
                earlier = [{k: v for k, v in old["polished"].items() if k not in ("state", "rounds", "torsions_deg", "earlier")}] if "polished" in old else []
                entry[seg] = {**old, "polished": {
                    "method": "on the closure surface, one curvature block per chain as the metric" if args.surface else "SLSQP, the closure as equality constraints",
                    **({"earlier": old["polished"].get("earlier", []) + earlier} if earlier else {}),
                    "saddle": (info["energy"] - eA) / nm, "band_top": old["saddle"], "band_force": old["force_max"], "force_along": info["force"], "gradient_across": info["gradient_max"],
                    "curvature_along": info["slope"], "displacement": info["displacement"], "converged": info["converged"], "evaluations": info["evaluations"],
                    "rounds": info["rounds"], "torsions_deg": wrap(st.tors[m0]).tolist(), "seconds": time.time() - t0, "state": xs.tolist()}}
                mine.add(seg)
                log(f"   {seg}: polished saddle {entry[seg]['polished']['saddle']:+.4f} (the band's top {old['saddle']:+.4f} at force {old['force_max']:.3f}); force along the tangent "
                    f"{info['force']:+.4f}, across {info['gradient_max']:.1e}, curvature along {'not measured' if info['slope'] is None else format(info['slope'], '+.3f')}; "
                    f"{info['evaluations']} evaluations, {time.time() - t0:.0f} s{'' if info['converged'] else ' (NOT converged)'}")
                save()
                continue
            if same_ends and old.get("converged") and old.get("fmax", prior.get("fmax", args.fmax)) <= args.fmax:  # a tighter --fmax continues a stored band
                entry[seg] = old
                log(f"   {seg}: saddle {old['saddle']:+.3f} (from the stored run)")
                continue
            if not do(seg):
                continue
            ck = checkpoint(seg)
            start = (f"the stored band (saddle {old['saddle']:+.4f} at force {old['force_max']:.3f} after {old.get('iterations', 0)} iterations), continued"
                     if same_ends else ("straight between the end points" if guide is None else f"the movers' internal coordinates from the band of {args.band_from!r}"))
            if ck is not None and ck.exists() and np.load(ck).shape == (args.images, model.n_var) and same(np.load(ck)[0], ends[a]) and same(np.load(ck)[-1], ends[b]):
                imgs, start = np.load(ck), start + f"; resumed from {ck}"  # a checkpoint is removed when its band is stored, so one that is there is newer than the record
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

            res = MC.neb(evaluate, imgs, scale=sc, free=free, k_spring=args.k_spring, climb_after=args.climb_after, climb_force=1.0, fmax=args.fmax, maxiter=args.maxiter,
                         constraint=constraint, callback=progress, snapshot=snapshot, max_step=args.band_step)
            i_top = int(np.argmax(res.energies))
            st = model.state(res.images[i_top])
            dips = [i for i in range(1, args.images - 1) if res.energies[i] < min(res.energies[i - 1], res.energies[i + 1]) - 1e-3]
            entry[seg] = {"energies": ((res.energies - eA) / nm).tolist(), "saddle": float(res.energies.max() - eA) / nm, "saddle_image": i_top, "converged": bool(res.converged),
                          "force_max": res.force_max, "fmax": args.fmax, "iterations": res.iterations + (old.get("iterations", 0) if same_ends else 0), "interior_minima": dips,
                          "start": old.get("start", start) if same_ends else start, "saddle_history": (old.get("saddle_history", [[old.get("iterations", 0), old["saddle"], old["force_max"]]]) if same_ends else [])
                          + [[res.iterations + (old.get("iterations", 0) if same_ends else 0), float(res.energies.max() - eA) / nm, res.force_max]],
                          "saddle_torsions_deg": wrap(st.tors[m0]).tolist(), "seconds": time.time() - t0, "images": res.images.tolist()}
            mine.add(seg)
            log(f"   {seg}: saddle {(res.energies.max() - eA) / nm:+.3f} at image {i_top}{'' if res.converged else ' (NOT converged)'}; force {res.force_max:.4f}; "
                f"{res.iterations} it, {time.time() - t0:.0f} s{'; interior minima at ' + str(dips) if dips else ''}")
            save()
            if ck is not None and ck.exists():
                ck.unlink()  # the stored record is the band now
        if pool is not None:
            pool.close()

        def summarise(disk):  # the barriers, from the row as it now stands on disk (another process may have written the other band)
            row = disk.get("rows", {}).get(key, {})
            if all(k in row for k in ("A-I", "I-B", "B")):
                top = max(0.0, float(row["A-I"]["saddle"]), float(row["I-B"]["saddle"]))
                row.update(barrier_forward=top, barrier_back=top - float(row["B"]["energy"]), moving_chains=nm,
                           bands_converged=bool(row["A-I"].get("converged") and row["I-B"].get("converged")))
                settled_barrier(row, log)
                log(f"   one-monomer step: barrier {top:+.3f} ahead, {top - float(row['B']['energy']):+.3f} back; B lies {float(row['B']['energy']):+.3f} above A"
                    f"{'' if row['bands_converged'] else ' (a band is NOT converged)'}")
        locked_update(out_path, summarise)
    log(f"written {out_path}")


if __name__ == "__main__":
    main()
