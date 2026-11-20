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


def relax_sheet(model, x, e_flat, c0, starts, maxiter, t0, say=print):
    """The sheet at ``x`` (upright) relaxed from the kicked upright cell and from each start ``(ma, mb)`` in monomers."""
    out = []
    for start in [None] + list(starts):
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


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--key", required=True, help='a kinked crystal of --kinked-references, as "vdf10@0" (with --sites: the scan\'s row, as "vdf10-trfe4-14@0~")')
    ap.add_argument("--sites", default="", help="comma list of sites of --scan (offsets in bonds); without it the kinked crystal itself")
    ap.add_argument("--scan", type=Path, default=PRED / "kink_site_scan_belt4.json")
    ap.add_argument("--kinked-references", type=Path, default=PRED / "kink_migration_belt4.json")
    ap.add_argument("--dilute", type=Path, default=PRED / "dilute_kink_belt4.json", help="the store with the row's all-trans cell and the preset")
    ap.add_argument("--preset", default="", help="evaluate under this preset instead of the store's (each stored cell is first relaxed upright under it)")
    ap.add_argument("--starts", default="-1:-1,1:1,-1:1,1:-1", help="starting shears in monomers along a and b, 'ma:mb,...'; the upright cell is always tried")
    ap.add_argument("--maxiter", type=int, default=800)
    ap.add_argument("--redo", action="store_true")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    row = args.key.split("@")[0]
    store = json.loads(args.dilute.read_text())
    flat = arrays(store["all_trans_references"][row])
    kinked = arrays(json.loads(args.kinked_references.read_text())["kinked_references"][args.key])
    preset = args.preset or store["preset"]
    own = preset == store["preset"]
    starts = [tuple(float(v) for v in item.split(":")) for item in args.starts.split(",") if item]
    out = json.loads(args.out.read_text()) if args.out.exists() else {"what": WHAT, "rows": {}}
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
            cells = [(f"{args.key} site {site} {preset}", np.concatenate([np.asarray(scan["sites"][site]["state"], dtype=float), np.zeros(2)]), float(scan["sites"][site]["formation"]),
                      {"site": site, "units": scan["sites"][site]["units"]}) for site in args.sites.split(",") if site]
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
                   "starts": relax_sheet(model, x, e_flat, c0, starts, args.maxiter, t0, say)}
            held = [s for s in rec["starts"] if s.get("kink_kept") and "failed" not in s]
            rec["lowest"] = min(s["formation"] for s in held) if held else None
            out = json.loads(args.out.read_text()) if args.out.exists() else out   # another run may have added to the file meanwhile
            out["rows"][key] = rec
            args.out.write_text(json.dumps(out, indent=1))
            say(f"{key}: upright {upright:+.3f}, lowest with the shear free " + ("no kink kept" if rec["lowest"] is None else f"{rec['lowest']:+.3f}") + f" kcal/mol per kink pair; written {args.out}")


if __name__ == "__main__":
    main()
