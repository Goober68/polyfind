"""A stored kink sheet relaxed again with the cell free to shear along its chains.

Every kink cell of examples/dilute_kink.py, kink_migration.py and kink_site_scan.py is upright: the
kinks of all the chains sit at one height.  DFT's stress on such a sheet, and a relaxation of it in
our own potential to a tight stop, say it is not a minimum: the cell shears until each chain's kink
is one monomer along from its neighbour's, an inclined sheet, and the energy falls by kcal/mol per
kink pair (docs/MLIP_SECOND_OPINION.md 3.19).  This gives polyfind's cell the same freedom
(``polyfind.mixedcell.MixedModel(shear=True)``, as examples/shear_polish.py does for the margins'
cells) and relaxes a stored sheet from the upright cell and from starts sheared by whole monomers
along ``a`` and ``b``, with the row's all-trans cell beside it.

    python examples/kink_sheet_shear.py --key "vdf10@0" --out docs/predictions/kink_sheet_shear_belt4.json
    python examples/kink_sheet_shear.py --key "vdf10-trfe4-14@0" --sites 12 --out docs/predictions/kink_sheet_shear_belt4.json

The first is a kinked crystal of the packing search (``kink_migration_belt4.json``); the second a
site of examples/kink_site_scan.py.  The energies are with the pendants as built (rigid), as the
scans' are.  The file is added to: a key that is there is run again only with ``--redo``.
"""
import argparse
import importlib.util
import json
import time
from pathlib import Path

import numpy as np

from polyfind.fitting import parameters_of

HERE = Path(__file__).resolve().parent
PRED = HERE.parent / "docs" / "predictions"
SPEC = importlib.util.spec_from_file_location("shear_polish", HERE / "shear_polish.py")
SP = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(SP)
DK = SP.DK
WHAT = ("a stored kink sheet (a kink pair on every chain) relaxed with the cell free to shear along its chains, from the upright cell and from starts sheared by whole monomers "
        "along a and b; kcal/mol per kink pair against the row's all-trans cell relaxed the same way, pendants as built. shear_monomers: what a and b carry along the chains at "
        "the end, in monomers; kink_kept: every chain still has two gauche bonds of opposite sense about a trans one (examples/kink_sheet_shear.py)")


def arrays(d):
    return {k: (np.asarray(v, dtype=float) if isinstance(v, list) else v) for k, v in d.items()}


def kink_kept(tors) -> bool:
    """Whether a chain's torsions (degrees) hold one kink pair: two gauche bonds of opposite sense two bonds apart, the rest trans."""
    t = (np.asarray(tors, dtype=float) + 180.0) % 360.0 - 180.0
    g = [(k, v) for k, v in enumerate(t) if abs(v) < 120.0]
    return len(g) == 2 and g[0][1] * g[1][1] < 0 and (g[1][0] - g[0][0]) in (2, len(t) - 2)


def upright_minimum(model, x, maxiter):
    """``x`` relaxed with the shear held at none: for cells stored under another preset."""
    everything, _ = SP.free_sets(model)
    x1, _ = model.relax(x, free=np.setdiff1d(everything, model.index("shear")), closed=[0, 1], maxiter=maxiter, patience=60, patience_tol=1e-6)
    return x1


def relax_states(model, states, e_flat, c0, maxiter, t0, say=print):
    """Each of ``states`` (full states of this cell, their shears with them) relaxed as it stands: sheets another preset found, asked of this one."""
    everything, _ = SP.free_sets(model)
    out = []
    for x in states:
        x = np.asarray(x, dtype=float)
        start = [float(v) / c0 for v in model.state(x).shear]
        try:
            x1, info = model.relax(x, free=everything, closed=[0, 1], maxiter=maxiter, patience=60, patience_tol=1e-6)
            if info["stalled"] or info["gradient_max"] > 1e-3:
                x1, last = model.relax_on_surface(x1, free=everything, closed=[0, 1], gtol=1e-4, maxiter=300, patience=60)
                info = {**info, "gradient_max": last["gradient_max"], "stalled": last["stalled"]}
        except RuntimeError as err:
            out.append({"start_monomers": start, "from_a_stored_sheet": True, "failed": str(err)})
            say(f"  stored sheet at {start[0]:+.2f}, {start[1]:+.2f} monomer: FAILED, {err} ({time.time() - t0:.0f} s)")
            continue
        e1 = float(model.physical_energy(x1))
        if abs(e1 - e_flat) > 1e3:
            out.append({"start_monomers": start, "from_a_stored_sheet": True, "failed": "collapsed"})
            say(f"  stored sheet at {start[0]:+.2f}, {start[1]:+.2f} monomer: FAILED, collapsed ({time.time() - t0:.0f} s)")
            continue
        s1 = model.state(x1)
        tors = (np.asarray(s1.tors, dtype=float).reshape(model.nc, -1) + 180.0) % 360.0 - 180.0
        out.append({"start_monomers": start, "from_a_stored_sheet": True, "shear_monomers": [float(v) / c0 for v in s1.shear], "formation": (e1 - e_flat) / 2.0,
                    "kink_kept": bool(all(kink_kept(row) for row in tors)), "gauche_bonds_deg": [[round(float(v), 1) for v in row if abs(v) < 120.0] for row in tors],
                    "gradient_max": float(info["gradient_max"]), "stalled": bool(info["stalled"]), "lattice": [float(v) for v in s1.lattice], "state": [float(v) for v in x1]})
        st = out[-1]
        say(f"  stored sheet at {start[0]:+.2f}, {start[1]:+.2f} monomer: ends {st['shear_monomers'][0]:+.2f}, {st['shear_monomers'][1]:+.2f}; "
            f"{st['formation']:+.3f} kcal/mol per kink pair{'' if st['kink_kept'] else ', KINK NOT KEPT'} ({time.time() - t0:.0f} s)")
    return out


def stored_sheets(rec) -> list:
    """The distinct sheets of a site's record that kept the kink (one per cell, cells told apart to half a monomer of shear), lowest first."""
    seen, out = set(), []
    for s in sorted(held_starts(rec), key=lambda s: s["formation"]):
        cell = tuple(round(2.0 * v) for v in s["shear_monomers"])
        if cell not in seen:
            seen.add(cell)
            out.append(s["state"])
    return out


def relax_sheet(model, x, e_flat, c0, starts, maxiter, t0, say=print, upright=True):
    """The sheet at ``x`` (upright) relaxed from the kicked upright cell (unless ``upright`` is off) and from each start ``(ma, mb)`` in monomers."""
    out = []
    for start in ([None] if upright else []) + list(starts):
        try:
            x1, info, (t,) = SP.relax_sheared(model, x, starts=() if start is None else ((start[0] * c0, start[1] * c0),), maxiter=maxiter, upright=start is None)
        except RuntimeError as err:  # a start the chains cannot be closed from: recorded, not a result
            out.append({"start_monomers": [0.0, 0.0] if start is None else list(start), "failed": str(err)})
            say(f"  start {out[-1]['start_monomers'][0]:+.2f}, {out[-1]['start_monomers'][1]:+.2f} monomer: FAILED, {err} ({time.time() - t0:.0f} s)")
            continue
        if abs(t["energy"] - e_flat) > 1e3:  # atoms driven through each other from a start this preset has no cell near: not a structure
            out.append({"start_monomers": [v / c0 for v in t["start_shear_A"]], "failed": "collapsed"})
            say(f"  start {out[-1]['start_monomers'][0]:+.2f}, {out[-1]['start_monomers'][1]:+.2f} monomer: FAILED, collapsed ({time.time() - t0:.0f} s)")
            continue
        s1 = model.state(x1)
        tors = (np.asarray(s1.tors, dtype=float).reshape(model.nc, -1) + 180.0) % 360.0 - 180.0
        out.append({"start_monomers": [v / c0 for v in t["start_shear_A"]], "shear_monomers": [v / c0 for v in t["shear_A"]], "formation": (t["energy"] - e_flat) / 2.0,
                    "kink_kept": bool(all(kink_kept(row) for row in tors)), "gauche_bonds_deg": [[round(float(v), 1) for v in row if abs(v) < 120.0] for row in tors],
                    "gradient_max": float(info["gradient_max"]), "stalled": bool(info["stalled"]), "lattice": [float(v) for v in s1.lattice], "state": [float(v) for v in x1]})
        st = out[-1]
        say(f"  start {st['start_monomers'][0]:+.2f}, {st['start_monomers'][1]:+.2f} monomer: ends {st['shear_monomers'][0]:+.2f}, {st['shear_monomers'][1]:+.2f}; "
            f"{st['formation']:+.3f} kcal/mol per kink pair{'' if st['kink_kept'] else ', KINK NOT KEPT'} ({time.time() - t0:.0f} s)")
    return out


def held_starts(rec) -> list:
    """The starts of a record that ended in a sheet with its kink kept."""
    return [s for s in rec["starts"] if s.get("kink_kept") and "failed" not in s]


def gathered(files) -> dict:
    """The site records of several result files; a site in more than one has its starts put together."""
    done = {}
    for f in files:
        for key, rec in json.loads(Path(f).read_text())["rows"].items():
            if " site " not in key:
                continue
            if key in done:
                seen = {tuple(round(v, 2) for v in s["start_monomers"]) for s in done[key]["starts"]}
                done[key]["starts"] += [s for s in rec["starts"] if tuple(round(v, 2) for v in s["start_monomers"]) not in seen]
            else:
                done[key] = rec
            held = held_starts(done[key])
            done[key]["lowest"] = min(s["formation"] for s in held) if held else None
    return done


def as_scan(args):
    """The scan of ``--scan`` with each site's energy replaced by the lowest sheet found here with the shear free: the file
    ``examples/pin_design.py --sites`` reads.  A site not run here keeps the scan's upright energy and is listed."""
    scan = json.loads(args.scan.read_text())
    done = gathered(args.merge)
    rows, left = {}, []
    for key, row in scan["rows"].items():
        mine = {k.split(" site ")[1].split(" ")[0]: v for k, v in done.items() if k.split(" site ")[0] == key}
        if not mine:
            continue
        row = json.loads(json.dumps(row))
        for site, v in row["sites"].items():
            if site in mine and mine[site]["lowest"] is not None:
                best = min(held_starts(mine[site]), key=lambda s: s["formation"])
                v.update(upright_formation=v["formation"], formation=best["formation"], shear_monomers=best["shear_monomers"], kink_kept=True, converged=True)
                v.pop("state", None)
            elif site in mine:
                v.update(upright_formation=v["formation"], kink_kept=False)
            elif v.get("converged", True) and v.get("kink_kept", True):
                left.append(f"{key} site {site}")
        rows[key] = row
    out = {**{k: v for k, v in scan.items() if k != "rows"}, "rows": rows,
           "shear_free": "each site's formation is the lowest kink sheet found with the cell free to shear along its chains (examples/kink_sheet_shear.py), against the "
                         "all-trans cell relaxed the same way; upright_formation is the scan's own", "sites_left_upright": left}
    args.as_scan.write_text(json.dumps(out, indent=1))
    free = sum(1 for r in rows.values() for v in r["sites"].values() if "shear_monomers" in v)
    print(f"{len(rows)} rows, {free} sites with the shear free, {len(left)} left upright; written {args.as_scan}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--key", default="", help='a kinked crystal of --kinked-references, as "vdf10@0" (with --sites: the scan\'s row, as "vdf10-trfe4-14@0~")')
    ap.add_argument("--merge", nargs="*", default=[], help="with --as-scan: result files of this script")
    ap.add_argument("--as-scan", type=Path, default=None, help="write --scan again with the sites of --merge at their shear-free energies, and stop")
    ap.add_argument("--gather", type=Path, default=None, help="write the site records of --merge as one result file (one site's starts put together), and stop")
    ap.add_argument("--sites", default="", help="comma list of sites of --scan (offsets in bonds), or 'all' (every site the scan kept a kink at); without it the kinked crystal itself")
    ap.add_argument("--scan", type=Path, default=PRED / "kink_site_scan_belt4.json")
    ap.add_argument("--kinked-references", type=Path, default=PRED / "kink_migration_belt4.json")
    ap.add_argument("--dilute", type=Path, default=PRED / "dilute_kink_belt4.json", help="the store with the row's all-trans cell and the preset")
    ap.add_argument("--preset", default="", help="evaluate under this preset instead of the store's (each stored cell is first relaxed upright under it)")
    ap.add_argument("--starts", default="-1:-1,1:1,-1:1,1:-1", help="starting shears in monomers along a and b, 'ma:mb,...'; the upright cell is always tried")
    ap.add_argument("--from-sheets", type=Path, default=None, help="with --sites and --preset: a result file of this script under the store's preset; each site is also relaxed "
                                                                   "from every distinct sheet stored for it there (a preset that slides a rigid start back does not find them itself)")
    ap.add_argument("--maxiter", type=int, default=800)
    ap.add_argument("--redo", action="store_true")
    ap.add_argument("--no-upright", action="store_true", help="leave the kicked upright start out: more --starts for sites stored in another file (--merge puts one site's starts together)")
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()
    if args.gather:
        rows = gathered(args.merge)
        args.gather.write_text(json.dumps({"what": WHAT, "rows": {k: rows[k] for k in sorted(rows, key=lambda k: (k.split(" site ")[0], int(k.split(" site ")[1].split(" ")[0])))}}, indent=1))
        print(f"{len(rows)} sites; written {args.gather}")
    if args.as_scan:
        return as_scan(args)
    if args.gather:
        return None
    if not args.key or args.out is None:
        ap.error("--key and --out are needed")
    row = args.key.split("@")[0]
    store = json.loads(args.dilute.read_text())
    flat = arrays(store["all_trans_references"][row])
    kinked = arrays(json.loads(args.kinked_references.read_text())["kinked_references"][args.key])
    preset = args.preset or store["preset"]
    own = preset == store["preset"]
    starts = [tuple(float(v) for v in item.split(":")) for item in args.starts.split(",") if item]
    out = json.loads(args.out.read_text()) if args.out.exists() else {"what": WHAT, "rows": {}}
    theirs = gathered([args.from_sheets]) if args.from_sheets else None
    t0 = time.time()
    say = lambda text: print(text, flush=True)  # noqa: E731
    with parameters_of(preset).applied():
        model, x, ch = SP.sheared_model(DK.C.make_polymer(row, DK.sequence(row)[0]), flat, check=own)
        c0 = float(model.state(x).lattice[3]) / ch.n_monomers
        if not own:
            x = upright_minimum(model, x, args.maxiter)
        e_flat_upright = float(model.physical_energy(x))
        x_flat, _, _ = SP.relax_sheared(model, x, maxiter=args.maxiter)
        e_flat = float(model.physical_energy(x_flat))
        flat_shear = [float(v) / c0 for v in model.state(x_flat).shear]
        say(f"{row} all-trans under {preset}: upright {e_flat_upright:.4f}, with the shear free {e_flat:.4f} kcal/mol per cell, shear {flat_shear[0]:+.2f}, {flat_shear[1]:+.2f} monomer "
            f"({time.time() - t0:.0f} s)")
        if args.sites:
            scan = json.loads(args.scan.read_text())["rows"][args.key]
            polymer = DK.C.make_polymer(row, scan["letters"])
            wanted = [o for o, v in scan["sites"].items() if v.get("converged", True) and v.get("kink_kept", True)] if args.sites == "all" else [o for o in args.sites.split(",") if o]
            cells = [(f"{args.key} site {site} {preset}", np.concatenate([np.asarray(scan["sites"][site]["state"], dtype=float), np.zeros(2)]), float(scan["sites"][site]["formation"]),
                      {"site": site, "units": scan["sites"][site]["units"]}) for site in wanted]
        else:
            polymer = DK.C.make_polymer(row, DK.sequence(row)[0])
            cells = [(f"{args.key} {preset}", None, (float(kinked["energy"]) - float(flat["energy"])) / 2.0, {})]
        model, x_ref, ch = SP.sheared_model(polymer, kinked, check=own)
        for key, x, stored, more in cells:
            if key in out["rows"] and not args.redo:
                say(f"{key}: stored, skipped")
                continue
            x = x_ref if x is None else x
            if not own:
                x = upright_minimum(model, x, args.maxiter)
            upright = (float(model.physical_energy(x)) - e_flat_upright) / 2.0
            if own and abs(upright - stored) > 1e-3:
                raise RuntimeError(f"{key}: the upright sheet is {upright:+.4f} here and {stored:+.4f} in the store")
            say(f"{key} upright: {upright:+.3f} kcal/mol per kink pair" + ("" if own else f" (the store's preset: {stored:+.3f})"))
            rec = {"preset": preset, **more, "monomer_repeat_A": c0, "all_trans": {"upright": e_flat_upright, "shear_free": e_flat, "shear_monomers": flat_shear}, "upright": upright,
                   "starts": relax_sheet(model, x, e_flat, c0, starts, args.maxiter, t0, say, upright=not args.no_upright)}
            if theirs is not None and f"{key.rsplit(' ', 1)[0]} {store['preset']}" in theirs:
                rec["starts"] += relax_states(model, stored_sheets(theirs[f"{key.rsplit(' ', 1)[0]} {store['preset']}"]), e_flat, c0, args.maxiter, t0, say)
            held = held_starts(rec)
            rec["lowest"] = min(s["formation"] for s in held) if held else None
            out = json.loads(args.out.read_text()) if args.out.exists() else out   # another run may have added to the file meanwhile
            out["rows"][key] = rec
            args.out.write_text(json.dumps(out, indent=1))
            say(f"{key}: upright {upright:+.3f}, lowest with the shear free " + ("no kink kept" if rec["lowest"] is None else f"{rec['lowest']:+.3f}") + f" kcal/mol per kink pair; written {args.out}")


if __name__ == "__main__":
    main()
