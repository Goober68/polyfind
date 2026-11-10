"""Kink pairs in a cell whose chains may differ: one kinked chain among all-trans neighbours, or several.

examples/kink_formation_screen.py kinks every chain of the packer's two-chain cell in register,
which is a kinked crystal.  polyfind.mixedcell is the cell in which the chains differ, and this
script asks it the questions the kinked crystal cannot answer:

* **a lone kink pair** (``--kinked 0``): the refined all-trans two-chain cell tiled ``nx x ny``
  (four chains at 2x1, eight at 2x2), chain 0 given T G+ T G- T at ``--kink-at``
  (``--mirror-kink`` for G- T G+), relaxed in stages that are each an upper bound on the next:
  ``rigid_host`` (the all-trans neighbours fixed), ``host_rigid_bodies`` (each neighbour free to
  shift and turn as a rigid chain), ``flexible_host`` (every chain's torsions and backbone angles
  free), all three with the cell held at the all-trans values, which is the boundary condition of
  a kink far from any other; then ``c_free`` (the repeat free) and ``cell_free`` (the whole cell),
  which are a kink density of one per supercell;
* **several chains kinked in register** (``--kinked 0,1`` ... ``--kinked all``): the same stages,
  so the energy per kink pair as a function of how many neighbours kink with it is read off one
  model, and ``--kinked all`` at ``cell_free`` is the kinked crystal of the screen, here with the
  chain closure held exactly.

Every chain is periodic exactly (the closure constraints of polyfind.mixedcell), so a kinked chain,
which is about 1 A shorter per kink pair than the all-trans stem, reaches the cell's repeat by
flattening its kink and opening its angles.  Formation energies are in kcal/mol per kink pair
against the all-trans supercell, which is the tiled two-chain cell (checked to be stationary and
stable under every freedom the stages release).  A kink whose gauche bonds have returned to within
30 degrees of trans is reported as collapsed: the structure is then the all-trans one.
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
STAGES = ("rigid_host", "host_rigid_bodies", "flexible_host", "c_free", "cell_free")


def sequence(name: str):
    return C.SEQUENCES.get(name) or EXTRA[name]


def all_trans_reference(polymer, name: str, stored: dict | None = None) -> dict:
    """The refined, relaxed all-trans two-chain cell: cell parameters, repeat, torsions, angles."""
    if stored:
        return {k: (np.asarray(v, dtype=float) if isinstance(v, list) else v) for k, v in stored.items()}
    t0 = time.time()
    ref, rr = C.E2.build(polymer, [T] * polymer.bonds_per_repeat, f"{name} TTTT")
    return {"params": np.asarray(ref.params, dtype=float), "c": float(ref.c), "tors": np.asarray(rr.torsions, dtype=float),
            "angs": np.asarray(rr.angles, dtype=float), "E_cell": float(ref.packer.energy(np.asarray(ref.params)[None])[0]),
            "bend_per_monomer": float(rr.angle_energy), "seconds": time.time() - t0}


def kinked_torsions(tors, offset: int, mirror: bool) -> np.ndarray:
    """All-trans torsions with G+ T G- (G- T G+ if ``mirror``) on bonds ``offset .. offset + 2``."""
    t = np.array(tors, dtype=float)
    P = len(t)
    g = -GAUCHE if mirror else GAUCHE
    t[offset % P], t[(offset + 1) % P], t[(offset + 2) % P] = g, 180.0, -g
    return t


def make_model(polymer, d: dict, nx: int, ny: int):
    """``(model, x_all_trans, all-trans chain)`` for the ``nx x ny`` tiling of the reference cell."""
    P = polymer.bonds_per_repeat
    ang0 = np.array([polymer.backbone[k].backbone_angle for k in range(P)], dtype=float)
    frame = repeat_chains(polymer, "TT", np.full((1, P), 180.0), ang0)[0]
    tt = repeat_chains(polymer, "TT", d["tors"][None], d["angs"], align_to=frame)[0]
    nc = 2 * nx * ny
    pk = pack_mod.CrystalPacker(tt, n_chains=nc, **KW)
    st, flips = MC.supercell_state(d["params"], tt.c, d["tors"], d["angs"], nx, ny)
    model = MC.MixedModel(polymer, "TT", pk, flips, [tt] * nc, angle_stiffness=resolve_angle_stiffness(None, None))
    return model, st.vector(), tt


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
    w = lambda a: float((a + 180.0) % 360.0 - 180.0)  # noqa: E731
    t = [w(v) for v in st.tors[chain]]
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
    ap.add_argument("--maxiter", type=int, default=400)
    ap.add_argument("--angle-window", type=float, default=8.0, help="bound on each backbone angle about the polymer's own, degrees (the refinement's)")
    ap.add_argument("--stages", default=",".join(STAGES))
    ap.add_argument("--tilt-stages", default="", help="comma list of stages after which the pendant tilts of the cell are relaxed")
    ap.add_argument("--out", default="docs/predictions/dilute_kink_belt4.json")
    args = ap.parse_args()
    out_path = Path(args.out)
    results = json.load(open(out_path)) if out_path.exists() else {}
    results.setdefault("rows", {})
    results.update(preset=args.preset, energy="kcal/mol per kink pair: the cell with the kinked chains minus the tiled all-trans supercell, divided by the number of kinked chains; "
                   "every chain periodic exactly (polyfind.mixedcell closure constraints)")
    letters, note = sequence(args.row)
    log = lambda s: print(s, flush=True)  # noqa: E731
    mine: dict = {"rows": set(), "refs": set()}

    def save():
        """Write this run's own rows and references into the file as it stands now, under a lock,
        so several runs can share one output file."""
        import fcntl
        with open(str(out_path) + ".lock", "w") as lock:
            fcntl.flock(lock, fcntl.LOCK_EX)
            disk = json.load(open(out_path)) if out_path.exists() else {}
            disk.setdefault("rows", {})
            disk.setdefault("all_trans_references", {})
            disk.update(preset=results["preset"], energy=results["energy"])
            for k in mine["rows"]:
                disk["rows"][k] = results["rows"][k]
            for k in mine["refs"]:
                disk["all_trans_references"][k] = results["all_trans_references"][k]
            with open(out_path, "w") as fh:
                json.dump(disk, fh, indent=1)

    with PRESET_PARAMETERS[args.preset].applied():
        polymer = C.make_polymer(args.row, letters)
        refs = results.setdefault("all_trans_references", {})
        d = all_trans_reference(polymer, args.row, refs.get(args.row))
        if args.row not in refs:
            mine["refs"].add(args.row)
        refs[args.row] = {k: (v.tolist() if isinstance(v, np.ndarray) else v) for k, v in d.items()}
        save()
        model, x_tt, tt = make_model(polymer, d, args.nx, args.ny)
        nc, n_mon = model.nc, tt.n_monomers
        kinked = list(range(nc)) if args.kinked == "all" else sorted(int(v) for v in args.kinked.split(","))
        key = f"{args.row}@{args.kink_at}{'~' if args.mirror_kink else ''} {args.nx}x{args.ny} kinked {','.join(map(str, kinked))}"
        e_cell = d["E_cell"] + d["bend_per_monomer"] * 2 * n_mon
        e_tt = model.physical_energy(x_tt)
        if abs(e_tt * 2 / nc - e_cell) > 1e-6:
            raise RuntimeError("the tiled all-trans supercell does not reproduce the two-chain cell")
        entry = results["rows"].setdefault(key, {})
        mine["rows"].add(key)
        entry.update(row=args.row, letters=letters, note=note, supercell=[args.nx, args.ny], chains=nc, sites=int(model.pk.N), monomers_per_chain=int(n_mon),
                     kinked_chains=kinked, kink_offset_bonds=args.kink_at, kink_mirrored=bool(args.mirror_kink), angle_window_deg=args.angle_window,
                     all_trans_cell={"params": d["params"].tolist(), "c": d["c"], "energy_per_monomer": e_cell / (2 * n_mon)},
                     all_trans_supercell_energy=e_tt)
        log(f"{key}: {nc} chains of {n_mon} monomers, {model.pk.N} sites; all-trans cell {np.round(d['params'], 4).tolist()} c {d['c']:.4f}, {e_cell / (2 * n_mon):.5f} kcal/mol per monomer")
        sets = stage_sets(model, kinked)
        x = x_tt.copy()
        for s in kinked:
            x[model.index("tors", s)] = kinked_torsions(d["tors"], args.kink_at, args.mirror_kink)
        start = model.terms(x)
        entry["ideal_kink"] = {"natural_repeat_A": float(start["repeats"][kinked[0]]), "shortening_A": float(start["repeats"][kinked[0]] - d["c"]),
                               "rotation_error_deg": float(start["rotation_error"][kinked[0]])}
        log(f"   ideal kink: natural repeat {start['repeats'][kinked[0]]:.4f} A against the cell's {d['c']:.4f} ({start['repeats'][kinked[0]] - d['c']:+.3f}), rotation error {start['rotation_error'][kinked[0]]:.3f} deg")
        x = model.project(x, model.index("internal", kinked), kinked)  # the kinked chains made periodic with the cell's repeat
        entry["projected_start"] = {"energy_above_all_trans": model.physical_energy(x) - e_tt, **kink_shape(model, x, kinked[0], args.kink_at, args.mirror_kink)}
        entry.setdefault("stages", {})
        states = entry.setdefault("states", {})
        tilt_stages = [s for s in args.tilt_stages.split(",") if s]
        tilt_ref = None
        for stage in [s for s in args.stages.split(",") if s]:
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
            states[stage] = x.tolist()
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
