"""The barrier to sliding one chain of beta-PVDF along its axis, with the lattice around it responding.

docs/CANDIDATE_MATERIALS.md 8.1 quotes 6.80 kcal/mol per monomer for one all-trans chain slid
rigidly by half a monomer in a rigid lattice.  That is an upper bound: the neighbours of a
half-slipped chain can turn about their axes, move aside, and flex.  This script holds one chain
of an ``nx x ny`` tiling of the polished two-chain cell at a given slip and relaxes the rest, at
fixed cell (one chain cannot change the cell of a crystal; ``--cell-free`` gives the other bound,
one slid chain per supercell at zero lateral stress), at three levels:

``turn``      every chain turns about its own axis (setting angles);
``shift``     and moves laterally (the slid chain too);
``flexible``  and every torsion and backbone angle is free, each chain periodic exactly.

**Why the half-slip is the saddle.**  The all-trans chain has a mirror plane normal to its axis,
so at a slip of half a monomer the axial force on every chain vanishes by symmetry whatever the
chains do laterally, and the state is stationary along the slide.  If it is also a minimum in
every other coordinate with the slip held, the full Hessian has exactly one negative eigenvalue
(interlacing) and the state is the transition state of the slide: the barrier needs no path.
A relaxation started exactly on the symmetric half-slipped lattice stays on it, a saddle held by
symmetry, so every relaxation is also started from ``--seeds`` perturbed copies and the lowest
result is kept; the axial forces at the result are printed as the check of stationarity.  The
other chains' axial shifts are held during the relaxation and released afterwards (``drag``:
all but one far chain free) to show the state does not move.

A one-monomer repeat has no internal freedom (two torsions and two angles against four closure
conditions), so ``flexible`` needs ``--repeat 2`` at least, which admits the alternately deflected
chain; ``--repeat 10`` in a small cell checks that nothing of longer wavelength along the chain
is softer.  The rigid numbers are independent of the repeat (tests/test_chain_slide.py).

    python examples/chain_slide.py --repeat 1 --cells 2x1,4x2,6x3,8x4 --levels turn,shift
    python examples/chain_slide.py --repeat 2 --cells 2x1,4x2 --levels turn,shift,flexible
    python examples/chain_slide.py --repeat 1 --cells 4x2 --levels shift --profile 8

Energies are kcal/mol per monomer of the slid chain, rigid pendants unless stated, against the
unslid lattice (which is the minimum at every level: the polished cell).  ``--profile K`` adds
the relaxed energy at K slips up to half a monomer, each continued from the last and from
perturbed copies, there and back.
"""
import argparse
import importlib.util
import json
import time
from pathlib import Path

import numpy as np

from polyfind.fitting import PRESET_PARAMETERS

DS = importlib.util.spec_from_file_location("dk", Path(__file__).with_name("dilute_kink.py"))
DK = importlib.util.module_from_spec(DS)
DS.loader.exec_module(DK)
MC = DK.MC
LEVELS = ("turn", "shift", "flexible")


def reference(letters: str, stored: dict | None = None):
    """``(polymer, name, cell dict)`` of the polished all-trans two-chain cell of a PVDF repeat of
    ``len(letters)`` monomers (the packing search and its polish, or the stored cell)."""
    name = f"vdf{len(letters)}"
    polymer = DK.C.make_polymer(name, letters)
    return polymer, name, DK.all_trans_reference(polymer, name, stored)


def axis_positions(model, x) -> np.ndarray:
    """Cartesian lateral position (A) of every chain axis."""
    st = model.state(x)
    A, B, gam = (float(v) for v in st.lattice[:3])
    a1 = np.array([A, 0.0])
    a2 = B * np.array([np.cos(np.radians(gam)), np.sin(np.radians(gam))])
    return st.place[:, :1] * a1 + st.place[:, 1:2] * a2


def distances_from(model, x, chain: int = 0) -> np.ndarray:
    """Minimum-image lateral distance (A) of every chain axis from chain ``chain``."""
    st = model.state(x)
    A, B, gam = (float(v) for v in st.lattice[:3])
    a1 = np.array([A, 0.0])
    a2 = B * np.array([np.cos(np.radians(gam)), np.sin(np.radians(gam))])
    f = st.place[:, :2] - st.place[chain, :2]
    f -= np.round(f)
    best = np.full(len(f), np.inf)
    for i in (-1, 0, 1):
        for j in (-1, 0, 1):
            d = (f[:, :1] + i) * a1 + (f[:, 1:2] + j) * a2
            best = np.minimum(best, np.linalg.norm(d, axis=1))
    return best


def free_set(model, level: str, anchor: int, drag: bool = False, cell: bool = False) -> np.ndarray:
    """The variables a level relaxes.  The slip (chain 0's axial shift) is never among them;
    ``anchor`` is the far chain whose lateral position fixes the lattice's.  ``drag`` frees the
    axial shifts of every chain but the slid one and the anchor; ``cell`` frees the lateral cell
    (a, b, gamma), never the repeat."""
    nc = model.nc
    others = [s for s in range(nc) if s != anchor]
    idx = [model.index("phi")]
    if level in ("shift", "flexible"):
        idx.append(model.index("xy", others))
    if level == "flexible":
        idx.append(model.index("internal"))
    if drag:
        idx.append(model.index("z", [s for s in range(nc) if s not in (0, anchor)]))
    if cell:
        idx.append(model.index("lattice"))
    return np.unique(np.concatenate(idx).astype(int))


def perturbed(model, x, level: str, anchor: int, seed: int, turn: float = 3.0, shift: float = 0.05, torsion: float = 1.0):
    """``x`` with its symmetry broken in the variables of ``level``: setting angles within ``turn``
    degrees, lateral positions within ``shift`` A, torsions within ``torsion`` degrees (closure restored)."""
    rng = np.random.default_rng(seed)
    y = np.asarray(x, dtype=float).copy()
    st = model.state(y)
    y[model.index("phi")] += rng.uniform(-turn, turn, model.nc)
    if level in ("shift", "flexible"):
        others = [s for s in range(model.nc) if s != anchor]
        y[model.index("xy", others)] += (rng.uniform(-shift, shift, (len(others), 2)) / np.array([st.lattice[0], st.lattice[1]])).ravel()
    if level == "flexible":
        y[model.index("tors")] += rng.uniform(-torsion, torsion, model.nc * model.P)
        y = model.project(y, model.index("internal"), list(range(model.nc)))
    return y


def relax(model, x, free, maxiter: int):
    """SLSQP over ``free`` until it reports success (three restarts at most)."""
    xr, info = model.relax(x, free=free, maxiter=maxiter)
    for _ in range(3):
        if "success" in info["message"].lower():
            break
        xr, more = model.relax(xr, free=free, maxiter=maxiter)
        info = {**more, "iterations": info["iterations"] + more["iterations"]}
    return xr, info


def describe(model, x, x_ref, anchor: int) -> dict:
    """What the lattice did: turns and lateral moves by distance from the slid chain, the torsions,
    and the axial forces (the gradient along each chain's axial shift), which vanish at a saddle."""
    st, s0 = model.state(x), model.state(x_ref)
    turn = DK.wrap(st.place[:, 2] - s0.place[:, 2])
    move = np.linalg.norm(axis_positions(model, x) - axis_positions(model, x_ref), axis=1)
    dist = distances_from(model, x_ref, 0)
    _, g = model.value_and_grad(x)
    gz = g[model.index("z")]
    shells = {}
    for s in range(1, model.nc):
        key = f"{dist[s]:.2f}"
        e = shells.setdefault(key, {"chains": 0, "turn": 0.0, "move": 0.0})
        e["chains"] += 1
        e["turn"] = max(e["turn"], float(abs(turn[s])))
        e["move"] = max(e["move"], float(move[s]))
    near = sorted(shells, key=float)[:4]
    return {"slid_chain_turn_deg": float(turn[0]), "slid_chain_move_A": float(move[0]),
            "largest_turn_deg": float(np.abs(turn).max()), "largest_move_A": float(move.max()),
            "largest_torsion_off_trans_deg": float(np.abs(DK.wrap(st.tors - 180.0)).max()),
            "axial_force_on_slid_chain": float(gz[0]), "largest_axial_force_on_another": float(np.abs(gz[1:]).max()),
            "shells_A": {k: {"chains": shells[k]["chains"], "largest_turn_deg": shells[k]["turn"], "largest_move_A": shells[k]["move"]} for k in near}}


def held_slip(model, x_ref, x_start, level: str, anchor: int, seeds, maxiter: int, n_mon: int, log, extra_starts=(), cell: bool = False):
    """The lowest relaxed state at the slip of ``x_start`` over the symmetric start, ``extra_starts``
    (states of a lower level or a neighbouring slip, moved to this slip by the caller) and perturbed
    copies of each.  Returns ``(x, row)``, energies per monomer of the slid chain."""
    free = free_set(model, level, anchor, cell=cell)
    e_ref = model.physical_energy(x_ref)
    starts = [("as started", x_start)] + [(f"continued {k}", y) for k, y in enumerate(extra_starts)]
    starts += [(f"{lab}, seed {s}", perturbed(model, y, level, anchor, s)) for lab, y in list(starts) for s in seeds]
    found = []
    for label, y in starts:
        t0 = time.time()
        xr, info = relax(model, y, free, maxiter)
        e = (model.physical_energy(xr) - e_ref) / n_mon
        found.append((e, label, xr, info))
        log(f"         {label:28s} {e:+9.4f}   {info['iterations']:4d} it, {info['message']}, {time.time() - t0:.0f} s")
    found.sort(key=lambda f: f[0])
    e, label, xb, info = found[0]
    row = {"barrier": float(e), "from": label, "message": info["message"], "starts": {lab: float(v) for v, lab, _, _ in found},
           "distinct_minima": [float(v) for v in sorted({round(v, 3) for v, _, _, _ in found})], **describe(model, xb, x_ref, anchor)}
    return xb, row


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--repeat", type=int, default=1, help="monomers in the chain repeat (1: rigid-body levels only; 2 or more for 'flexible')")
    ap.add_argument("--row", default="", help="a sequence row of examples/dilute_kink.py in place of PVDF (its polished all-trans cell must be stored there)")
    ap.add_argument("--c-per-monomer", type=float, default=None,
                    help="set the repeat to this many A per monomer and relax the two-chain reference laterally at it (for a comparison at another seat's chain repeat); rows are tagged")
    ap.add_argument("--period", action="store_true",
                    help="--profile over a whole monomer period, for a chain without the mirror plane (a copolymer: its slide has no symmetric saddle)")
    ap.add_argument("--cells", default="2x1,4x2", help="comma list of tilings nx x ny of the two-chain cell")
    ap.add_argument("--levels", default="turn,shift")
    ap.add_argument("--seeds", type=int, default=2, help="perturbed starts per start")
    ap.add_argument("--profile", type=int, default=0, help="K: also the relaxed energy at K slips up to half a monomer, at the last level given")
    ap.add_argument("--tilts", action="store_true", help="add the pendant-tilt relaxation at the last level's half-slip")
    ap.add_argument("--cell-free", action="store_true",
                    help="also the last level with the lateral cell free: one slid chain per cell at zero lateral stress, which bounds the dilute limit from below as the fixed cell does from above")
    ap.add_argument("--preset", default="pvdf-dft-valence-belt4")
    ap.add_argument("--maxiter", type=int, default=600)
    ap.add_argument("--kink-references", default="docs/predictions/dilute_kink_belt4.json", help="polished all-trans cells of longer repeats, used if present")
    ap.add_argument("--out", default="docs/predictions/chain_slide_belt4.json")
    args = ap.parse_args()
    log = lambda s: print(s, flush=True)  # noqa: E731
    levels = [v for v in args.levels.split(",") if v]
    if any(v not in LEVELS for v in levels):
        raise SystemExit(f"levels are {LEVELS}")
    if "flexible" in levels and args.repeat < 2 and not args.row:
        raise SystemExit("a one-monomer repeat has no internal freedom: use --repeat 2 or more for 'flexible'")
    out_path = Path(args.out)
    seeds = list(range(args.seeds))

    def save(key, row):
        import fcntl
        with open(str(out_path) + ".lock", "w") as lock:
            fcntl.flock(lock, fcntl.LOCK_EX)
            disk = json.load(open(out_path)) if out_path.exists() else {}
            disk.update(preset=args.preset, energy="kcal/mol per monomer of the slid chain against the unslid lattice; fixed cell; rigid pendants unless stated")
            disk.setdefault("rows", {}).setdefault(key, {}).update(row)
            disk["rows"] = dict(sorted(disk["rows"].items()))
            with open(out_path, "w") as fh:
                json.dump(disk, fh, indent=1)

    with PRESET_PARAMETERS[args.preset].applied():
        kinks = Path(args.kink_references)
        if args.row:
            letters, _ = DK.sequence(args.row)
            stored = json.load(open(kinks)).get("all_trans_references", {}).get(args.row) if kinks.exists() else None
            if stored is None:
                raise SystemExit(f"no polished all-trans cell of {args.row} in {kinks}")
            polymer, name = DK.C.make_polymer(args.row, letters), args.row
            d = DK.all_trans_reference(polymer, name, stored)
            n_mon = len(letters)
        else:
            stored = (json.load(open(out_path)) if out_path.exists() else {}).get("references", {}).get(f"vdf{args.repeat}")
            from_kinks = stored is None and kinks.exists() and f"vdf{args.repeat}" in json.load(open(kinks)).get("all_trans_references", {})
            if from_kinks:  # the long repeats of examples/dilute_kink.py are polished already
                stored = json.load(open(kinks))["all_trans_references"][f"vdf{args.repeat}"]
            polymer, name, d = reference("V" * args.repeat, stored)
            if stored is None:
                import fcntl
                with open(str(out_path) + ".lock", "w") as lock:
                    fcntl.flock(lock, fcntl.LOCK_EX)
                    disk = json.load(open(out_path)) if out_path.exists() else {}
                    disk.setdefault("references", {})[name] = DK.as_json(d)
                    with open(out_path, "w") as fh:
                        json.dump(disk, fh, indent=1)
            n_mon = args.repeat
        tag = ""
        if args.c_per_monomer:
            # the two-chain cell at the given repeat: lateral cell and placements relaxed, the repeat held
            m2, x2, _ = DK.make_model(polymer, d, 1, 1)
            x2 = x2.copy()
            x2[3] = args.c_per_monomer * n_mon
            x2 = m2.project(x2, m2.index("internal"), [0, 1])
            free2 = np.setdiff1d(np.arange(m2.n_var), np.concatenate([m2.index("c"), m2.index("xy", 1), m2.index("z", 1)]))
            x2, info2 = relax(m2, x2, free2, args.maxiter)
            log(f"   repeat set to {args.c_per_monomer:.6f} A per monomer: lateral cell {np.round(x2[:3], 5).tolist()} (was {np.round(np.asarray(d['state'][:3], dtype=float), 5).tolist()}), "
                f"E {m2.physical_energy(x2) / (2 * n_mon):.6f} kcal/mol per monomer; {info2['message']}")
            d = dict(d, state=np.asarray(x2, dtype=float), energy=float(m2.physical_energy(x2)))
            tag = f" c{args.c_per_monomer:.4f}"
        cm = float(d["state"][3]) / n_mon
        log(f"{name}: cell {np.round(np.asarray(d['state'][:3], dtype=float), 5).tolist()}, monomer period {cm:.5f} A, E {float(d['energy']) / (2 * n_mon):.6f} kcal/mol per monomer")
        for cell in args.cells.split(","):
            nx, ny = (int(v) for v in cell.split("x"))
            key = f"{name} {nx}x{ny}{tag}"
            t_cell = time.time()
            model, x0, ch = DK.make_model(polymer, d, nx, ny)
            dist = distances_from(model, x0, 0)
            anchor = int(np.argmax(dist))
            iz = model.index("z", 0)[0]
            e0 = model.physical_energy(x0)
            xs = x0.copy()
            xs[iz] += 0.5 * cm
            rigid = (model.physical_energy(xs) - e0) / n_mon
            row = {"chains": int(model.nc), "sites": int(model.pk.N), "repeat_monomers": n_mon, "monomer_period_A": cm, "anchor_chain_distance_A": float(dist[anchor]),
                   "nearest_neighbours_A": [float(v) for v in np.sort(dist[1:])[:6]], "rigid": float(rigid)}
            log(f"{key}: {model.nc} chains, {model.pk.N} sites; neighbours at {np.round(np.sort(dist[1:])[:6], 2).tolist()} A, anchor at {dist[anchor]:.1f} A; rigid half-slip {rigid:+.4f}")
            save(key, row)
            carried, xb = [], xs
            for level in levels:
                log(f"   {level}:")
                xb, r = held_slip(model, x0, xs, level, anchor, seeds, args.maxiter, n_mon, log, extra_starts=carried)
                # the other chains' axial shifts released: a saddle does not move
                t0 = time.time()
                xd = xb.copy()
                free_d = free_set(model, level, anchor, drag=True)
                zs = model.index("z", [s for s in range(model.nc) if s not in (0, anchor)])
                if len(zs) == 0:
                    raise SystemExit("a two-chain cell has no single chain to slide: use 2x1 or larger")
                xd[zs] += np.random.default_rng(99).uniform(-0.01, 0.01, len(zs))
                xd, info = relax(model, xd, free_d, args.maxiter)
                r["drag_released"] = {"energy_change": float((model.physical_energy(xd) - model.physical_energy(xb)) / n_mon),
                                      "largest_axial_shift_A": float(np.abs(xd[zs] - x0[zs]).max()), "message": info["message"]}
                log(f"      {level}: barrier {r['barrier']:+.4f} ({r['barrier'] / rigid:.2f} of rigid; from '{r['from']}'; minima seen {r['distinct_minima']}); slid chain turned "
                    f"{r['slid_chain_turn_deg']:+.2f} deg, moved {r['slid_chain_move_A']:.3f} A; largest turn {r['largest_turn_deg']:.2f} deg, move {r['largest_move_A']:.3f} A, torsion "
                    f"{r['largest_torsion_off_trans_deg']:.2f} deg off trans; axial force on the slid chain {r['axial_force_on_slid_chain']:+.1e}, on another {r['largest_axial_force_on_another']:.1e}; "
                    f"axial shifts released: {r['drag_released']['energy_change']:+.1e}, largest shift {r['drag_released']['largest_axial_shift_A']:.1e} A ({time.time() - t0:.0f} s)")
                log("         shells (distance A: chains, largest turn deg, largest move A): "
                    + "; ".join(f"{k}: {v['chains']}, {v['largest_turn_deg']:.2f}, {v['largest_move_A']:.3f}" for k, v in r["shells_A"].items()))
                row[level] = r
                save(key, {level: r})
                carried = [xb]
            last = levels[-1] if levels else None
            if args.cell_free and last:
                log(f"   {last}, lateral cell free:")
                xc, r = held_slip(model, x0, xs, last, anchor, seeds[:1], args.maxiter, n_mon, log, extra_starts=[xb], cell=True)
                sc, s0 = model.state(xc), model.state(x0)
                r["cell_change"] = {"a": float(sc.lattice[0] / s0.lattice[0] - 1), "b": float(sc.lattice[1] / s0.lattice[1] - 1), "gamma_deg": float(sc.lattice[2] - s0.lattice[2])}
                log(f"      {last}, lateral cell free: barrier {r['barrier']:+.4f} ({r['barrier'] / rigid:.2f} of rigid); a {r['cell_change']['a']:+.3%}, b {r['cell_change']['b']:+.3%}, "
                    f"gamma {r['cell_change']['gamma_deg']:+.3f} deg; axial force on the slid chain {r['axial_force_on_slid_chain']:+.1e}")
                save(key, {f"{last}_cell_free": r})
            if args.tilts and last:
                t0 = time.time()
                dE, tilts, _ = MC.tilt_relaxation(model, xb)
                dE0, _, _ = MC.tilt_relaxation(model, x0)
                r = {"barrier": float(row[last]["barrier"] + (dE - dE0) / n_mon), "tilt_gain": float((dE - dE0) / n_mon), "largest_tilt_deg": float(np.abs(tilts).max())}
                log(f"   {last} with pendant tilts: {r['barrier']:+.4f} (tilts {r['tilt_gain']:+.4f}, largest {r['largest_tilt_deg']:.2f} deg); {time.time() - t0:.0f} s")
                save(key, {f"{last}_with_tilts": r})
            if args.profile and last:
                what = last + ("_cell_free" if args.cell_free else "")
                log(f"   profile at '{what}' (slip in monomer periods: relaxed, rigid):")
                fr = [(1.0 if args.period else 0.5) * (k + 1) / args.profile for k in range(args.profile)]
                prof, prev = {}, None
                for sweep, order in (("up", fr), ("down", fr[::-1])):
                    prev = None
                    for f in order:
                        xf = x0.copy()
                        xf[iz] += f * cm
                        rig = (model.physical_energy(xf) - e0) / n_mon
                        extra = []
                        if prev is not None:
                            y = prev.copy()
                            y[iz] = x0[iz] + f * cm
                            extra = [y]
                        xr, r = held_slip(model, x0, xf, last, anchor, seeds[:1], args.maxiter, n_mon, lambda s: None, extra_starts=extra, cell=args.cell_free)
                        prev = xr
                        e = prof.setdefault(f"{f:.4f}", {"rigid": float(rig)})
                        e[sweep] = float(r["barrier"])
                        e[f"{sweep}_axial_force"] = float(r["axial_force_on_slid_chain"])
                        log(f"      {sweep:4s} {f:.4f}: {r['barrier']:+.4f} (rigid {rig:+.4f}); axial force on the slid chain {r['axial_force_on_slid_chain']:+.3f} kcal/mol/A per repeat")
                for e in prof.values():
                    e["relaxed"] = float(min(e.get("up", np.inf), e.get("down", np.inf)))
                top = max(prof, key=lambda k: prof[k]["relaxed"])
                rtop = max(prof, key=lambda k: prof[k]["rigid"])
                log(f"      highest relaxed {prof[top]['relaxed']:+.4f} at {top} of a period (rigid {prof[rtop]['rigid']:+.4f} at {rtop}); "
                    f"largest difference between the sweeps {max(abs(e.get('up', e['relaxed']) - e.get('down', e['relaxed'])) for e in prof.values()):.4f}")
                save(key, {f"profile_{what}": prof})
            log(f"   {key} done in {time.time() - t_cell:.0f} s")
    log(f"written {out_path}")


if __name__ == "__main__":
    main()
