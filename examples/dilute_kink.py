"""Kink pairs in a cell whose chains may differ: one kinked chain among all-trans neighbours, or several.

examples/kink_formation_screen.py kinks every chain of the packer's two-chain cell in register,
which is a kinked crystal.  polyfind.mixedcell is the cell in which the chains differ, and this
script asks it the questions the kinked crystal cannot answer:

* **a lone kink pair** (``--kinked 0``): the all-trans two-chain cell tiled ``nx x ny`` (four
  chains at 2x1, eight at 2x2), chain 0 given T G+ T G- T at ``--kink-at`` (``--mirror-kink`` for
  G- T G+), relaxed in stages that are each an upper bound on the next: ``rigid_host`` (the
  all-trans neighbours fixed), ``host_rigid_bodies`` (each neighbour free to shift and turn as a
  rigid chain), ``flexible_host`` (every chain's torsions and backbone angles free), all three
  with the cell held at the all-trans values, which is the boundary condition of a kink far from
  any other; then ``c_free`` (the repeat free) and ``cell_free`` (the whole cell), which are a kink
  density of one per supercell;
* **several chains kinked in register** (``--kinked 0,1`` ... ``--kinked all``): the same stages,
  so the energy per kink pair as a function of how many neighbours kink with it is read off one
  model;
* **the axial registry** (stage ``slide``): the energy of sliding chain 0 rigidly along its axis
  through one monomer period in the all-trans supercell, which is what holds a kinked chain's
  1 A length deficit in place.

**The all-trans reference** is the packing search's refined two-chain cell, polished in the mixed
cell: every torsion and backbone angle of both chains free, the cell free, each chain periodic
exactly.  The search refines an all-trans chain under the line group of its state sequence, which
for a chiral (copolymer) sequence is not a symmetry of the molecule and holds the chain planar:
the polished chain is 0.1 kcal/mol per monomer lower for a VDF/TrFE sequence, and a kinked chain,
which the search refines with every torsion free, would otherwise be compared with a reference
that was not.  Stage ``host_check`` then relaxes the tiled all-trans supercell with every chain
free and takes that as the reference if it is lower (chains that want to differ along a or b).

Every chain is periodic exactly (the closure constraints of polyfind.mixedcell), so a kinked chain,
which is about 1 A shorter per kink pair than the all-trans stem, reaches the cell's repeat by
flattening its kink and opening its angles.  Formation energies are in kcal/mol per kink pair
against the all-trans supercell.  A kink whose gauche bonds have returned to within 30 degrees of
trans is reported as collapsed: the structure is then the all-trans one.
"""
import argparse
import importlib.util
import json
import time
from pathlib import Path

import numpy as np

from polyfind import mixedcell as MC
from polyfind import pack as pack_mod
from polyfind.ewald import EwaldSpec
from polyfind.fitting import PRESET_PARAMETERS
from polyfind.linegroup import repeat_chains
from polyfind.polarizability import Polarizable
from polyfind.refine import resolve_angle_stiffness

CS = importlib.util.spec_from_file_location("cs", Path(__file__).with_name("candidate_screen.py"))
C = importlib.util.module_from_spec(CS)
CS.loader.exec_module(C)
T = 0
EXTRA = {"vdf10": ("VVVVVVVVVV", "PVDF as a ten-monomer sequence"), "vdf7-trfe3": ("VVRVVRVVVR", "the 70/30 model of the E2 set"),
         "vdf14": ("V" * 14, "PVDF as a fourteen-monomer sequence"), "vdf20": ("V" * 20, "PVDF as a twenty-monomer sequence")}
KW = dict(lj_cutoff="force", coulomb="ewald", ewald=EwaldSpec(), polarizable=Polarizable())
GAUCHE = 60.0
K_ANGLE = None  # backbone-angle restraint, kcal/mol/rad^2; None: the preset's (--angle-stiffness sets it for a sensitivity run)


def angle_stiffness() -> float:
    return resolve_angle_stiffness(K_ANGLE, None)
STAGES = ("host_check", "slide", "rigid_host", "host_rigid_bodies", "flexible_host", "c_free", "cell_free")


def sequence(name: str):
    return C.SEQUENCES.get(name) or EXTRA[name]


def wrap(a):
    return (np.asarray(a, dtype=float) + 180.0) % 360.0 - 180.0


def reference_chain(polymer, d: dict):
    """The chain of a refined cell in the frame its placement refers to."""
    P = len(d["tors"])
    if d.get("frame") is not None:
        frame = np.asarray(d["frame"], dtype=float)
    else:  # an all-trans cell stored without its frame: the ideal all-trans chain's, which is the same
        ang0 = np.array([polymer.backbone[k].backbone_angle for k in range(polymer.bonds_per_repeat)], dtype=float)
        frame = repeat_chains(polymer, "TT", np.full((1, P), 180.0), ang0)[0]
    return repeat_chains(polymer, "ref", np.asarray(d["tors"], dtype=float)[None], np.asarray(d["angs"], dtype=float), align_to=frame)[0]


def polish(polymer, d: dict, seed: int = 0, maxiter: int = 800) -> dict:
    """The search's two-chain cell relaxed in the mixed cell: every torsion and angle of both
    chains free, the cell free, the closure exact.  The torsions start half a degree off the
    search's so that a stationary point held by a symmetry the molecule does not have is left."""
    ch = reference_chain(polymer, d)
    pk = pack_mod.CrystalPacker(ch, n_chains=2, **KW)
    st, flips = MC.supercell_state(d["params"], ch.c, d["tors"], d["angs"], 1, 1)
    model = MC.MixedModel(polymer, "ref", pk, flips, [ch, ch], angle_stiffness=angle_stiffness())
    x0 = st.vector()
    # the search's bend energy was computed with the preset's restraint; it is linear in the restraint
    e_search = float(d["E_cell"]) + float(d["bend_per_monomer"]) * (angle_stiffness() / resolve_angle_stiffness(None, None)) * 2 * ch.n_monomers
    if abs(model.physical_energy(x0) - e_search) > 1e-6:
        raise RuntimeError("the mixed cell does not reproduce the packing search's two-chain cell")
    residual = model.closure_residual(x0)
    x = x0.copy()
    x[model.index("tors")] += np.random.default_rng(seed).uniform(-0.5, 0.5, 2 * model.P)
    x = model.project(x, model.index("internal"), [0, 1])
    free = np.setdiff1d(np.arange(model.n_var), np.concatenate([model.index("xy", 1), model.index("z", 1)]))
    t0 = time.time()
    x1, info = model.relax(x, free=free, closed=[0, 1], maxiter=maxiter)
    return {"state": x1, "flips": flips.astype(float), "frame": np.asarray(ch.coords, dtype=float), "energy": model.physical_energy(x1),
            "search_energy": e_search, "search_closure_residual": [float(residual[0]), float(residual[1])],
            "polish": {"iterations": info["iterations"], "evaluations": info["evaluations"], "message": info["message"], "seconds": time.time() - t0}}


def refined_cell(polymer, row: str, seq, label: str, stored: dict | None = None) -> dict:
    """The packing search's refined two-chain cell of the state sequence ``seq`` (cell parameters,
    repeat, torsions, angles, energy, frame), and its polish in the mixed cell (``state``, the
    variable vector of the two-chain mixed cell; ``energy``, kcal/mol per two-chain cell)."""
    d = {k: (np.asarray(v, dtype=float) if isinstance(v, list) else v) for k, v in stored.items()} if stored else None
    if d is None:
        t0 = time.time()
        ref, rr = C.E2.build(polymer, seq, label)
        d = {"params": np.asarray(ref.params, dtype=float), "c": float(ref.c), "tors": np.asarray(rr.torsions, dtype=float),
             "angs": np.asarray(rr.angles, dtype=float), "E_cell": float(ref.packer.energy(np.asarray(ref.params)[None])[0]),
             "bend_per_monomer": float(rr.angle_energy), "frame": np.asarray(ref.packer.chain.coords, dtype=float), "seconds": time.time() - t0}
    if "state" not in d:
        d.update(polish(polymer, d))
    return d


def all_trans_reference(polymer, name: str, stored: dict | None = None) -> dict:
    """The refined, polished all-trans two-chain cell."""
    return refined_cell(polymer, name, [T] * polymer.bonds_per_repeat, f"{name} TTTT", stored)


def as_json(d: dict) -> dict:
    return {k: (v.tolist() if isinstance(v, np.ndarray) else v) for k, v in d.items()}


def make_model(polymer, d: dict, nx: int, ny: int):
    """``(model, x, chain)``: the ``nx x ny`` tiling of the polished two-chain cell ``d``."""
    ch = reference_chain(polymer, d)
    nc = 2 * nx * ny
    pk = pack_mod.CrystalPacker(ch, n_chains=nc, **KW)
    st, flips = MC.tile_state(MC.MixedState.from_vector(d["state"], 2, len(d["tors"])), np.asarray(d["flips"]) > 0.5, nx, ny)
    model = MC.MixedModel(polymer, "ref", pk, flips, [ch] * nc, angle_stiffness=angle_stiffness())
    x = st.vector()
    if abs(model.physical_energy(x) * 2 / nc - float(d["energy"])) > 1e-6:
        raise RuntimeError("the tiled supercell does not reproduce the polished two-chain cell")
    return model, x, ch


def kinked_torsions(tors, offset: int, mirror: bool) -> np.ndarray:
    """``tors`` with G+ T G- (G- T G+ if ``mirror``) on bonds ``offset .. offset + 2``."""
    t = np.array(tors, dtype=float)
    P = len(t)
    g = -GAUCHE if mirror else GAUCHE
    t[offset % P], t[(offset + 1) % P], t[(offset + 2) % P] = g, 180.0, -g
    return t


def stage_sets(model, kinked) -> dict:
    """``{stage: (free indices, chains whose closure is held)}`` for the kinked chains ``kinked``.

    One chain's lateral position and axial shift stay fixed (the last host chain, or the last
    chain when every chain is kinked), which removes the cell's rigid translation."""
    nc = model.nc
    kinked = sorted(int(s) for s in kinked)
    host = [s for s in range(nc) if s not in kinked]
    anchor = host[-1] if host else kinked[-1]
    place = lambda chains: np.concatenate([model.index("place", [s for s in chains if s != anchor]),  # noqa: E731
                                           model.index("phi", [s for s in chains if s == anchor])]).astype(int)
    k_int, h_int = model.index("internal", kinked), (model.index("internal", host) if host else np.zeros(0, dtype=int))
    k_place, h_place = place(kinked), (place(host) if host else np.zeros(0, dtype=int))
    everyone = list(range(nc))
    sets = {}
    if host:
        sets["rigid_host"] = (np.concatenate([k_int, k_place]), kinked)
        sets["host_rigid_bodies"] = (np.concatenate([k_int, k_place, h_place]), kinked)
    flexible = np.concatenate([k_int, k_place, h_place, h_int]).astype(int)
    sets["flexible_host"] = (flexible, everyone)
    sets["c_free"] = (np.concatenate([flexible, model.index("c")]), everyone)
    sets["cell_free"] = (np.concatenate([flexible, model.index("c"), model.index("lattice")]), everyone)
    return sets


def kink_shape(model, x, chain: int, offset: int, mirror: bool) -> dict:
    """The kink's torsions (bonds offset-2 .. offset+4, signed, deg), whether its gauche bonds are
    still gauche, and how far the rest of the chain and its backbone angles have moved."""
    st = model.state(x)
    P = model.P
    t = [float(v) for v in wrap(st.tors[chain])]
    g1, g2 = t[offset % P], t[(offset + 2) % P]
    rest = [abs(abs(t[(offset + k) % P]) - 180.0) for k in range(5, P - 2)]
    return {"kink_torsions_deg": [t[(offset + k) % P] for k in range(-2, 5)], "gauche_deg": [g1, g2],
            "collapsed": bool(abs(g1) > 150.0 and abs(g2) > 150.0),
            "largest_trans_deviation_elsewhere_deg": float(max(rest)) if rest else 0.0,
            "backbone_angle_range_deg": [float(st.angs[chain].min()), float(st.angs[chain].max())]}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--row", default="vdf10")
    ap.add_argument("--nx", type=int, default=2)
    ap.add_argument("--ny", type=int, default=1)
    ap.add_argument("--kinked", default="0", help="comma list of the chains that carry the kink pair, or 'all'")
    ap.add_argument("--kink-at", type=int, default=0, help="bond offset of the kink pair along the repeat")
    ap.add_argument("--mirror-kink", action="store_true", help="G- T G+ instead of G+ T G-")
    ap.add_argument("--preset", default="pvdf-dft-valence-belt4")
    ap.add_argument("--maxiter", type=int, default=600)
    ap.add_argument("--angle-window", type=float, default=8.0, help="bound on each backbone angle about the polymer's own, degrees (the refinement's)")
    ap.add_argument("--stages", default=",".join(STAGES))
    ap.add_argument("--tilt-stages", default="", help="comma list of stages after which the pendant tilts of the cell are relaxed")
    ap.add_argument("--angle-stiffness", type=float, default=None,
                    help="backbone-angle restraint in kcal/mol/rad^2 for a sensitivity run (default: the preset's 233). The mixed cell's chains have rigid bonds, "
                         "so with the preset's restraint their axial stiffness is about twice the full valence form's; about 114 matches that form's C33")
    ap.add_argument("--out", default="docs/predictions/dilute_kink_belt4.json")
    args = ap.parse_args()
    global K_ANGLE
    K_ANGLE = args.angle_stiffness
    tag = "" if K_ANGLE is None else f" k{K_ANGLE:g}"  # a sensitivity run keeps its own reference cell and its own rows
    wtag = "" if abs(args.angle_window - 8.0) < 1e-9 else f" w{args.angle_window:g}"  # and so does a run with another angle window
    out_path = Path(args.out)
    letters, note = sequence(args.row)
    log = lambda s: print(s, flush=True)  # noqa: E731
    entry: dict = {}
    key = None

    def locked_update(fn):
        import fcntl
        with open(str(out_path) + ".lock", "w") as lock:
            fcntl.flock(lock, fcntl.LOCK_EX)
            disk = json.load(open(out_path)) if out_path.exists() else {}
            disk.setdefault("rows", {})
            disk.setdefault("all_trans_references", {})
            disk.update(preset=args.preset, energy="kcal/mol per kink pair: the cell with the kinked chains minus the all-trans supercell, divided by the number of kinked chains; "
                        "every chain periodic exactly (polyfind.mixedcell closure constraints)")
            fn(disk)
            with open(out_path, "w") as fh:
                json.dump(disk, fh, indent=1)

    def save():
        locked_update(lambda disk: disk["rows"].__setitem__(key, entry))

    with PRESET_PARAMETERS[args.preset].applied():
        polymer = C.make_polymer(args.row, letters)
        refs = (json.load(open(out_path)) if out_path.exists() else {}).get("all_trans_references", {})
        stored = refs.get(args.row + tag)
        if stored is None and tag and args.row in refs:  # the search's cell is the same; only its polish depends on the restraint
            stored = {k: v for k, v in refs[args.row].items() if k in ("params", "c", "tors", "angs", "E_cell", "bend_per_monomer", "frame", "seconds")}
            if "frame" not in stored and "state" in refs[args.row]:
                stored["frame"] = refs[args.row]["frame"]
        d = all_trans_reference(polymer, args.row, stored)
        if stored is None or "state" not in stored:
            locked_update(lambda disk: disk["all_trans_references"].__setitem__(args.row + tag, as_json(d)))
        model, x_tt, tt = make_model(polymer, d, args.nx, args.ny)
        nc, n_mon, P = model.nc, tt.n_monomers, model.P
        kinked = list(range(nc)) if args.kinked == "all" else sorted(int(v) for v in args.kinked.split(","))
        key = f"{args.row}@{args.kink_at}{'~' if args.mirror_kink else ''} {args.nx}x{args.ny} kinked {','.join(map(str, kinked))}{tag}{wtag}"
        e_tt = model.physical_energy(x_tt)
        st_tt = model.state(x_tt)
        entry.update(row=args.row, letters=letters, note=note, supercell=[args.nx, args.ny], chains=nc, sites=int(model.pk.N), monomers_per_chain=int(n_mon),
                     kinked_chains=kinked, kink_offset_bonds=args.kink_at, kink_mirrored=bool(args.mirror_kink), angle_window_deg=args.angle_window,
                     angle_stiffness=angle_stiffness(),
                     all_trans_cell={"lattice": st_tt.lattice.tolist(), "energy_per_monomer": e_tt / (nc * n_mon),
                                     "search_energy_per_monomer": float(d["search_energy"]) / (2 * n_mon),
                                     "polish_gain_per_monomer": (float(d["search_energy"]) - float(d["energy"])) / (2 * n_mon),
                                     "largest_torsion_off_trans_deg": float(np.abs(np.abs(wrap(st_tt.tors)) - 180.0).max())},
                     stages={}, states={})
        log(f"{key}: {nc} chains of {n_mon} monomers, {model.pk.N} sites; all-trans cell {np.round(st_tt.lattice, 4).tolist()}, {e_tt / (nc * n_mon):.5f} kcal/mol per monomer "
            f"({entry['all_trans_cell']['polish_gain_per_monomer']:.4f} below the search's line-group cell; torsions up to {entry['all_trans_cell']['largest_torsion_off_trans_deg']:.1f} deg off trans)")
        stages = [s for s in args.stages.split(",") if s]
        sets = stage_sets(model, kinked)
        everyone = list(range(nc))

        if "host_check" in stages:  # the tiled all-trans supercell with every chain free: is it the reference?
            t0 = time.time()
            free, closed = stage_sets(model, [0])["flexible_host"]
            xp = x_tt.copy()
            xp[model.index("tors")] += np.random.default_rng(1).uniform(-0.5, 0.5, nc * P)
            xp = model.project(xp, model.index("internal"), everyone)
            xr, info = model.relax(xp, free=free, closed=closed, maxiter=args.maxiter, angle_window=args.angle_window)
            de = model.physical_energy(xr) - e_tt
            entry["host_check"] = {"energy_below_tiled": de, "per_monomer": de / (nc * n_mon), "iterations": info["iterations"], "message": info["message"],
                                   "seconds": time.time() - t0, "adopted": bool(de < -0.01)}
            log(f"   host check        : the all-trans supercell with every chain free is {de:+.4f} kcal/mol from the tiled cell ({de / (nc * n_mon):+.5f} per monomer); "
                f"{info['iterations']} it, {time.time() - t0:.0f} s; {info['message']}")
            if de < -0.01:
                x_tt, e_tt = xr, model.physical_energy(xr)
                entry["states"]["all_trans"] = x_tt.tolist()
                log("   the relaxed supercell is the reference from here on")
            save()
        entry["all_trans_supercell_energy"] = e_tt

        if "slide" in stages:  # the axial registry: chain 0 slid rigidly through one monomer period
            a_mon = st_tt.lattice[3] / n_mon
            iz = model.index("z", 0)[0]
            us = np.linspace(0.0, a_mon, 21)
            es = []
            for u in us:
                xs = x_tt.copy()
                xs[iz] += u
                es.append(model.physical_energy(xs) - e_tt)
            es = np.array(es)
            entry["slide"] = {"shift_A": us.tolist(), "energy": es.tolist(), "monomer_period_A": float(a_mon), "barrier_per_monomer": float(es.max() / n_mon),
                              "curvature_per_monomer": float(2.0 * es[1] / us[1] ** 2 / n_mon), "at_1A_per_monomer": float(np.interp(1.0, us, es) / n_mon)}
            log(f"   axial registry    : sliding chain 0 rigidly costs up to {es.max() / n_mon:.3f} kcal/mol per monomer over a period of {a_mon:.3f} A "
                f"({np.interp(1.0, us, es) / n_mon:.3f} at 1 A; curvature {entry['slide']['curvature_per_monomer']:.2f} kcal/mol/A^2 per monomer)")
            save()

        x = x_tt.copy()
        for s in kinked:
            x[model.index("tors", s)] = kinked_torsions(x_tt[model.index("tors", s)], args.kink_at, args.mirror_kink)
        start = model.terms(x)
        c_cell = float(st_tt.lattice[3])
        entry["ideal_kink"] = {"natural_repeat_A": float(start["repeats"][kinked[0]]), "shortening_A": float(start["repeats"][kinked[0]] - c_cell),
                               "rotation_error_deg": float(start["rotation_error"][kinked[0]])}
        log(f"   ideal kink        : natural repeat {start['repeats'][kinked[0]]:.4f} A against the cell's {c_cell:.4f} ({start['repeats'][kinked[0]] - c_cell:+.3f}), rotation error {start['rotation_error'][kinked[0]]:.3f} deg")
        x = model.project(x, model.index("internal", kinked), kinked)  # the kinked chains made periodic with the cell's repeat
        entry["projected_start"] = {"energy_above_all_trans": model.physical_energy(x) - e_tt, **kink_shape(model, x, kinked[0], args.kink_at, args.mirror_kink)}
        tilt_stages = [s for s in args.tilt_stages.split(",") if s]
        tilt_ref = None
        for stage in stages:
            if stage not in sets:
                continue
            free, closed = sets[stage]
            t0 = time.time()
            x, info = model.relax(x, free=free, closed=closed, maxiter=args.maxiter, angle_window=args.angle_window)
            e = model.physical_energy(x)
            st = model.state(x)
            shapes = [kink_shape(model, x, s, args.kink_at, args.mirror_kink) for s in kinked]
            rec = {"formation": (e - e_tt) / len(kinked), "energy_above_all_trans": e - e_tt, "iterations": info["iterations"], "evaluations": info["evaluations"],
                   "message": info["message"], "rotation_error_max_deg": info["rotation_error_max"], "c_mismatch_max_A": info["c_mismatch_max"],
                   "cell": st.lattice.tolist(), "seconds": time.time() - t0, "kinks": shapes, "collapsed": bool(all(sh["collapsed"] for sh in shapes)),
                   "host_setting_angles_deg": [float(st.place[s, 2]) for s in range(nc) if s not in kinked],
                   "angles_at_window": int(sum(int((np.abs(st.angs[s] - model.angles0) > args.angle_window - 1e-3).sum()) for s in range(nc)))}
            entry["stages"][stage] = rec
            entry["states"][stage] = x.tolist()
            log(f"   {stage:18s}: {rec['formation']:+8.3f} kcal/mol per kink pair{' (COLLAPSED to all-trans)' if rec['collapsed'] else ''}; gauche bonds {np.round(shapes[0]['gauche_deg'], 1).tolist()}, "
                f"angles {np.round(shapes[0]['backbone_angle_range_deg'], 1).tolist()} ({rec['angles_at_window']} at the window); cell {np.round(st.lattice, 3).tolist()}; "
                f"{info['iterations']} it, {info['evaluations']} ev, {rec['seconds']:.0f} s; {info['message']}")
            if stage in tilt_stages and not rec["collapsed"]:
                t0 = time.time()
                dE_k, tilts_k, _ = MC.tilt_relaxation(model, x)
                if tilt_ref is None:
                    tilt_ref = MC.tilt_relaxation(model, x_tt)[0]
                rec["tilts"] = {"kinked_cell": dE_k, "all_trans_supercell": tilt_ref, "per_kink_pair": (dE_k - tilt_ref) / len(kinked),
                                "largest_tilt_deg": float(np.abs(tilts_k).max()), "seconds": time.time() - t0}
                rec["formation_with_tilts"] = rec["formation"] + rec["tilts"]["per_kink_pair"]
                log(f"   {'':18s}  pendant tilts {rec['tilts']['per_kink_pair']:+.3f} per kink pair: {rec['formation_with_tilts']:+.3f} with them  ({time.time() - t0:.0f} s)")
            save()
            if rec["collapsed"]:
                log("   the kink has returned to all-trans; later stages would only repeat the all-trans cell")
                break
    log(f"written {out_path}")


if __name__ == "__main__":
    main()
