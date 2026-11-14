"""A step of the kink pair in a wider cell, started from a narrower one.

    # the stored states of a row, let into the cell twice as long along b, as a seed row for
    # kink_migration.py --seed-from (written to its own file: the stored table is not touched)
    python examples/kink_wide_cell.py seed --row "vdf10@0 8x2 mover 0 ahead host flexible" --along b \
        --name "seed 8x4" --out wide64.json
    python examples/kink_migration.py --row vdf10 --nx 8 --ny 4 --mover 0 --host flexible --surface \
        --trust-tiled --ends-only --seed-from "seed 8x4" --checkpoints ck --out wide64.json

    # a saddle of the wider cell from a narrower cell's: its state let in as the start, the tangent
    # of the band it came from let in through the same doublings, settled on the closure surface
    python examples/kink_wide_cell.py polish --segment I-B --nx 8 --ny 4 \
        --band-row "vdf10@0 4x2 mover 0 ahead host flexible" --band-path a b \
        --start-row "vdf10@0 8x2 mover 0 ahead host flexible" --start-path b \
        --wide-file wide64.json --blocks ck/vdf10_0.chain_blocks.npz

The cell-size series of the single-chain step (docs/CANDIDATE_MATERIALS.md 8.3) is made this way:
the one moving chain keeps the relaxed neighbours it has in the narrower cell, on both sides, and
the sheet is let in where that cell is farthest from it.  A relaxation then has only the new
room to use, tens of evaluations where a start from the sheet takes hundreds.  ``--seed-from``
of ``kink_migration.py`` takes only a row of the same cell; this makes one.

Chain order in a state vector is (layer along a, row along b, chain of the two-chain cell).
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
STORED = ROOT / "docs" / "predictions" / "kink_migration_belt4.json"
PER_CELL = 2  # chains of the two-chain cell


def let_in(x_small, x_sheet, nx: int, ny: int, along: str) -> np.ndarray:
    """The state of the cell twice as long along ``along`` ('a' or 'b').

    ``x_small`` is a state of the ``nx`` x ``ny`` cell in which chain 0 (layer 0, row 0) is the one
    that differs, ``x_sheet`` the sheet in that cell.  Along the doubled direction the layers on
    chain 0's side keep ``x_small`` -- going forward the first half of the cell (and always chain
    0's direct neighbour), going back through the boundary the other half -- and the layers in
    between are the sheet's.  The sheet let into itself is the sheet of the wider cell.
    """
    ax = {"a": 0, "b": 1}[along]
    n = (nx, ny)[ax]
    if n < 2:
        raise ValueError(f"a cell of one layer along {along} has no far side to let the sheet into")
    nc = PER_CELL * nx * ny
    x_small, x_sheet = np.asarray(x_small, dtype=float), np.asarray(x_sheet, dtype=float)
    width = (x_sheet.size - 4 - 4 * nc) // nc
    if x_sheet.size != 4 + nc * (4 + width) or x_small.size != x_sheet.size:
        raise ValueError(f"states of {x_small.size} and {x_sheet.size} variables for a {nx}x{ny} cell of {nc} chains")

    def parts(x):
        return x[4:4 + 4 * nc].reshape(nc, 4), x[4 + 4 * nc:].reshape(nc, width)

    (pS, iS), (pA, iA) = parts(x_small), parts(x_sheet)
    nx2, ny2 = (2 * nx, ny) if ax == 0 else (nx, 2 * ny)
    place, inner = np.zeros((2 * nc, 4)), np.zeros((2 * nc, width))
    for ix2 in range(nx2):
        for iy2 in range(ny2):
            j2 = (ix2, iy2)[ax]
            if j2 < max(n // 2, 2):         # chain 0's side of the narrower cell
                src_p, src_i, j, shift = pS, iS, j2, 0.0
            elif j2 >= 2 * n - n // 2:       # its far side, which is the near side through the boundary
                src_p, src_i, j, shift = pS, iS, j2 - n, float(n)
            else:                            # the sheet, let in between
                src_p, src_i, j, shift = pA, iA, j2 % n, float(n * (j2 // n))
            ix, iy = (j, iy2) if ax == 0 else (ix2, j)
            for k in range(PER_CELL):
                new, old = (ix2 * ny2 + iy2) * PER_CELL + k, (ix * ny + iy) * PER_CELL + k
                pl = src_p[old].copy()
                pl[ax] = (pl[ax] * n + shift) / (2.0 * n)
                place[new], inner[new] = pl, src_i[old]
    cell = x_sheet[:4].copy()
    cell[ax] *= 2.0
    return np.concatenate([cell, place.ravel(), inner.ravel()])


def through(x, x_sheet, nx: int, ny: int, path) -> tuple:
    """``x`` let in through the doublings of ``path`` (a list of 'a'/'b'), the sheet with it:
    ``(x_wide, sheet_wide, nx_wide, ny_wide)``."""
    x, x_sheet = np.asarray(x, dtype=float), np.asarray(x_sheet, dtype=float)
    for along in path:
        x, x_sheet = let_in(x, x_sheet, nx, ny, along), let_in(x_sheet, x_sheet, nx, ny, along)
        nx, ny = (2 * nx, ny) if along == "a" else (nx, 2 * ny)
    return x, x_sheet, nx, ny


def load_rows(path) -> dict:
    return json.loads(Path(path).read_text())


def seed(args) -> None:
    src = load_rows(args.source)
    row = src["rows"][args.row]
    nx, ny = row["supercell"]
    sheet = np.asarray(row["kinked_crystal"]["state"], dtype=float)
    path = list(args.along)
    new = {"kinked_crystal": {**{k: v for k, v in row["kinked_crystal"].items() if k not in ("state", "supercell_check")},
                              "state": through(sheet, sheet, nx, ny, path)[0].tolist()},
           "chains": int(row["chains"]) * 2 ** len(path), "seed_of": args.row, "let_in_along": path}
    for name in ("I", "B"):
        if "state" in row.get(name, {}):
            new[name] = {"state": through(row[name]["state"], sheet, nx, ny, path)[0].tolist()}
    out = Path(args.out)
    disk = load_rows(out) if out.exists() else {k: v for k, v in src.items() if k != "rows"}
    disk.setdefault("rows", {})[args.name] = new
    out.write_text(json.dumps(disk, indent=1))
    wide = through(sheet, sheet, nx, ny, path)[2:]
    print(f"written {out}: row {args.name!r}, {new['chains']} chains ({wide[0]}x{wide[1]}), states {[k for k in ('I', 'B') if k in new]} "
          f"of {args.row!r} let in along {' then '.join(path)}", flush=True)


def polish(args) -> None:
    from polyfind import mixedcell as MC
    from polyfind.fitting import PRESET_PARAMETERS

    spec = importlib.util.spec_from_file_location("kink_migration", ROOT / "examples" / "kink_migration.py")
    KM = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(KM)
    DK = KM.DK
    log = lambda s: print(s, flush=True)
    seg = args.segment
    stored = load_rows(args.source)
    band_row = stored["rows"][args.band_row]
    start_row = (load_rows(args.start_file) if args.start_file else stored)["rows"][args.start_row]
    bx, by = band_row["supercell"]
    sx, sy = start_row["supercell"]
    band_sheet = np.asarray(band_row["kinked_crystal"]["state"], dtype=float)
    start_sheet = np.asarray(start_row["kinked_crystal"]["state"], dtype=float)
    row_name, kink = band_row["row"], int(band_row["kink_offset_bonds"])
    kkey = f"{row_name}@{kink}{'~' if band_row.get('kink_mirrored') else ''}"
    letters, _ = DK.sequence(row_name)
    with PRESET_PARAMETERS[args.preset].applied():
        polymer = DK.C.make_polymer(row_name, letters)
        model, x_tiled, _ = DK.make_model(polymer, stored["kinked_references"][kkey], args.nx, args.ny)
        free, closed = KM.free_sets(model, [0], "flexible")
        free = np.unique(free)
        # the start: the narrower cell's saddle (its settled state if it has one, else its band's top image)
        if "state" in start_row.get("saddles", {}).get(seg, {}):
            x_from, what = start_row["saddles"][seg]["state"], f"the settled saddle of the {sx}x{sy} cell"
            known = float(start_row["saddles"][seg]["saddle"])
        elif "state" in start_row.get(seg, {}).get("polished", {}):
            x_from, what = start_row[seg]["polished"]["state"], f"the settled saddle of the {sx}x{sy} cell"
            known = float(start_row[seg]["polished"]["saddle"])
        else:
            x_from, what = start_row[seg]["images"][int(start_row[seg]["saddle_image"])], f"the top image of the {sx}x{sy} band"
            known = float(start_row[seg].get("polished", {}).get("saddle", start_row[seg]["saddle"]))
        x_top, sheet, wx, wy = through(x_from, start_sheet, sx, sy, list(args.start_path))
        if (wx, wy) != (args.nx, args.ny) or sheet.size != model.n_var:
            raise SystemExit(f"--start-path takes the {sx}x{sy} cell to {wx}x{wy}, not to {args.nx}x{args.ny}")
        eA = model.physical_energy(sheet)
        tiled = eA - model.physical_energy(x_tiled)
        log(f"{model.nc} chains ({args.nx}x{args.ny}); the sheet let in against the tiled two-chain cell: {tiled:+.2e} kcal/mol")
        if abs(tiled) > 1e-6:
            raise SystemExit("the sheet let into the wider cell is not the tiled cell: the chain order or the doublings are wrong")
        # the tangent: the band's, let in through the doublings from its own cell
        imgs, i = np.asarray(band_row[seg]["images"], dtype=float), int(band_row[seg]["saddle_image"])
        ahead, _, tx, ty = through(imgs[i + 1], band_sheet, bx, by, list(args.band_path))
        if (tx, ty) != (args.nx, args.ny):
            raise SystemExit(f"--band-path takes the {bx}x{by} cell to {tx}x{ty}, not to {args.nx}x{args.ny}")
        tangent = ahead - through(imgs[i - 1], band_sheet, bx, by, list(args.band_path))[0]
        turning = np.concatenate([model.index("tors"), model.index("phi")])
        tangent[turning] = KM.wrap(tangent[turning])
        if args.blocks and Path(args.blocks).exists():
            blocks = {int(k): v for k, v in np.load(args.blocks).items()}
        else:  # one chain's block hardly depends on the cell around it: from the four-chain tiling, as the driver takes them
            small, xs, _ = DK.make_model(polymer, stored["kinked_references"][kkey], 2, 1)
            blocks = MC.chain_blocks(small, xs, range(2))
            if args.blocks:
                np.savez(args.blocks, **{str(k): v for k, v in blocks.items()})
        metric = MC.chain_metric(model, sheet, free, closed, blocks)
        x_top = model.project(x_top, np.intersect1d(free, model.index("internal", closed)), closed)
        e_start = model.physical_energy(x_top) - eA
        log(f"{seg}: {what} let in: {e_start:+.4f} (in its own cell {known:+.4f}); closure residual {max(model.closure_residual(x_top)):.1e}")
        t0 = time.time()
        best = {"e": np.inf, "n": 0}
        ck = None if args.checkpoint is None else Path(args.checkpoint)

        def progress(stage, k, e, x):
            best["n"] += 1
            if e < best["e"]:
                best.update(e=e)
            if k == 1 or k % 10 == 0:
                log(f"      {stage} {k:4d}: {e - eA:+.4f}, lowest so far {best['e'] - eA:+.4f}  ({best['n']} evaluations, {time.time() - t0:.0f} s)")
                if ck is not None:
                    np.save(ck, x)

        xs, info = MC.saddle_polish(model, x_top, tangent, free, closed, tol=args.tol, on_surface=True, metric=metric, gtol=args.gtol,
                                    patience=args.patience, patience_tol=1e-5, progress=progress,
                                    log=lambda text: log(f"   {seg} polish, {text}  ({time.time() - t0:.0f} s)"))
        st = model.state(xs)
        rec = {"saddle": info["energy"] - eA, "start": f"{what} let into the {args.nx}x{args.ny} cell ({e_start:+.4f}), the tangent of the {bx}x{by} band let in with it",
               f"in_the_{sx}x{sy}_cell": known, "force_along": info["force"], "gradient_across": info["gradient_max"], "curvature_along": info["slope"],
               "displacement": info["displacement"], "converged": bool(info["converged"]), "evaluations": info["evaluations"], "rounds": info["rounds"],
               "method": "on the closure surface, one curvature block per chain as the metric", "torsions_deg": KM.wrap(st.tors[0]).tolist(),
               "torsions": KM.letters_of(st.tors[0]), "seconds": time.time() - t0, "state": xs.tolist(), "from_cell": f"{sx}x{sy}",
               "evaluated": "batched, on " + getattr(model.batched_module(), "__name__", "?") if model.batched_module() is not None else "loops, host"}
        log(f"{seg}: saddle {rec['saddle']:+.4f} in {model.nc} chains ({known:+.4f} in the {sx}x{sy} cell); along {info['force']:+.4f}, across {info['gradient_max']:.1e}, "
            f"curvature {info['slope']}; {info['evaluations']} evaluations, {time.time() - t0:.0f} s{'' if info['converged'] else ' (NOT converged)'}")
    out = Path(args.out)
    disk = load_rows(out) if out.exists() else {}
    disk.setdefault("saddles", {})[seg] = rec
    disk["cell"] = [args.nx, args.ny]
    out.write_text(json.dumps(disk, indent=1))
    log(f"written {out}")


def store(args) -> None:
    """A wider cell's row (and its saddles) from a job's files into the stored table; or, with
    ``--beside``, its energies alone beside the row the table already holds, as a second route."""
    wide = load_rows(args.wide_file)["rows"][args.row]
    table_path = Path(args.table)
    table = load_rows(table_path)
    saddles = load_rows(args.saddles)["saddles"] if args.saddles else {}
    if args.beside:
        have = table["rows"][args.row]
        again = [name for name in ("I", "B") if name in wide and not (name in have and have[name]["energy"] == wide[name]["energy"])]  # not the stored run itself
        rec = {"note": args.note, **{name: {k: v for k, v in wide[name].items() if k not in ("state", "torsions_deg")} for name in again}}
        if saddles:
            rec["saddles"] = {seg: {k: v for k, v in sd.items() if k not in ("state", "torsions_deg", "rounds")} for seg, sd in saddles.items()}
        have[args.beside] = rec
        for name in again:
            if name in have:
                print(f"   {name}: {wide[name]['energy']:+.4f} again, {have[name]['energy']:+.4f} stored ({wide[name]['energy'] - have[name]['energy']:+.4f})")
        for seg, sd in saddles.items():
            if seg in have.get("saddles", {}):
                print(f"   saddle {seg}: {sd['saddle']:+.4f} again, {have['saddles'][seg]['saddle']:+.4f} stored ({sd['saddle'] - have['saddles'][seg]['saddle']:+.4f})")
    else:
        row = {k: v for k, v in wide.items() if k not in ("seed_of", "let_in_along")}
        if args.note:
            row["evaluated"] = args.note
        before = table["rows"].get(args.row)
        if before is not None and args.keep_old:  # the row it replaces: its energies stay, as the other route
            row[args.keep_old] = {name: {k: v for k, v in before[name].items() if k not in ("state", "torsions_deg")} for name in ("I", "B") if name in before}
            for name in ("I", "B"):
                if name in row and name in before:
                    print(f"   {name}: {row[name]['energy']:+.4f}, {before[name]['energy']:+.4f} in the row it replaces ({row[name]['energy'] - before[name]['energy']:+.1e})")
        if saddles:
            row["saddles"] = saddles
            row["saddles_note"] = ("no band in this cell: each saddle settled on the closure surface from a narrower cell's saddle let into this one, "
                                   "with the tangent of the sixteen-chain band let in the same way (examples/kink_wide_cell.py polish)")
            done = [sd["saddle"] for sd in saddles.values() if sd.get("converged")]
            if len(done) == 2 and "B" in row:
                row["barrier_forward_polished"], row["barrier_back_polished"] = max(done), max(done) - row["B"]["energy"]
        table["rows"][args.row] = row
        print(f"   stored {args.row!r}: " + ", ".join(f"{name} {row[name]['energy']:+.4f}" for name in ("I", "B") if name in row)
              + "".join(f", saddle {seg} {sd['saddle']:+.4f}" for seg, sd in saddles.items()))
    tmp = table_path.with_suffix(".tmp")
    tmp.write_text(json.dumps(table, indent=1))
    tmp.replace(table_path)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="what", required=True)
    s = sub.add_parser("seed", help="a seed row for kink_migration.py --seed-from: a stored row's states let into a wider cell")
    s.add_argument("--source", default=str(STORED), help="the table holding the row")
    s.add_argument("--row", required=True)
    s.add_argument("--along", nargs="+", choices=("a", "b"), required=True, help="the doublings, in order")
    s.add_argument("--name", required=True, help="the seed row's name in --out")
    s.add_argument("--out", required=True, help="the file kink_migration.py will be given as --out (made from the table's references if it is not there)")
    s.set_defaults(run=seed)
    p = sub.add_parser("polish", help="a saddle of a wider cell from a narrower cell's saddle and band")
    p.add_argument("--segment", required=True, choices=("A-I", "I-B"))
    p.add_argument("--nx", type=int, required=True)
    p.add_argument("--ny", type=int, required=True)
    p.add_argument("--source", default=str(STORED), help="the stored table (bands, references)")
    p.add_argument("--band-row", required=True, help="the row whose band gives the tangent")
    p.add_argument("--band-path", nargs="+", choices=("a", "b"), required=True, help="the doublings from the band's cell to this one")
    p.add_argument("--start-row", required=True, help="the row whose saddle is the start")
    p.add_argument("--start-file", default=None, help="the file holding --start-row, if not the stored table")
    p.add_argument("--start-path", nargs="+", choices=("a", "b"), required=True, help="the doublings from the start row's cell to this one")
    p.add_argument("--blocks", default=None, help="the chains' curvature blocks (kink_migration.py --surface writes them into its --checkpoints)")
    p.add_argument("--checkpoint", default=None, help="a file for the running state")
    p.add_argument("--tol", type=float, default=0.01)
    p.add_argument("--gtol", type=float, default=1e-4)
    p.add_argument("--patience", type=int, default=40)
    p.add_argument("--preset", default="pvdf-dft-valence-belt4")
    p.add_argument("--out", required=True)
    p.set_defaults(run=polish)
    t = sub.add_parser("store", help="a job's row of a wider cell, and its saddles, into the stored table")
    t.add_argument("--wide-file", required=True, help="the job's table (jobs/shared/wide_steps.json)")
    t.add_argument("--row", required=True)
    t.add_argument("--saddles", default=None, help="the job's saddles of that cell (jobs/shared/saddles_NXxNY.json)")
    t.add_argument("--beside", default=None, metavar="FIELD",
                   help="the table holds this row already: keep it, and put this run's energies (no states) beside it under FIELD")
    t.add_argument("--keep-old", default=None, metavar="FIELD",
                   help="the table holds this row already: replace it, and keep the old row's energies (no states) in the new one under FIELD")
    t.add_argument("--note", default=None, help="where and how it was evaluated")
    t.add_argument("--table", default=str(STORED))
    t.set_defaults(run=store)
    args = ap.parse_args()
    args.run(args)


if __name__ == "__main__":
    main()
