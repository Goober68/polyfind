"""The dilute kink: one chain carrying a kink pair among all-trans neighbours.

examples/kink_formation_screen.py kinks every chain of the packer's two-chain cell in register,
which is a kinked crystal, and brackets the number a rate needs between that crystal relaxed and
the same chains held in the all-trans lattice.  The object the rate needs is neither: one kinked
chain whose neighbours are still all-trans and free to give way.  polyfind.mixedcell is the cell
for it, and this script is the calculation:

* the refined all-trans two-chain cell of the sequence, tiled ``nx x ny`` (four chains at 2x1,
  eight at 2x2), every chain with its own torsions, backbone angles, position, setting angle and
  axial shift;
* chain 0 given the kink pair (T G+ T G- T at ``--kink-at``, ``--mirror-kink`` for the other
  handedness), fitted into the cell's repeat, then relaxed in stages that are each an upper bound
  on the next: the host rigid and fixed, the host chains free as rigid bodies, every chain's
  conformation free with the cell held at the all-trans values (the dilute boundary condition),
  then the repeat ``c`` free and last the whole cell (a kink density of one per supercell);
* the pendant tilts of both cells relaxed afterwards, as in the screen;
* optionally (``--c-scan``) the cell's repeat scanned, which separates the tension the kinked
  chain carries for being shorter than its neighbours from the kink itself.

Formation energies are per kink pair against the all-trans supercell relaxed under the same
freedoms, in kcal/mol.  The closure penalties of the mixed cell are left out of the energies and a
first-order correction for the slack they allow (twice the penalty) is added; both are reported.
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


def make_model(polymer, d: dict, nx: int, ny: int, k_c: float, pen_rot: float):
    """``(model, x_all_trans, all-trans chain)`` for the ``nx x ny`` tiling of the reference cell."""
    P = polymer.bonds_per_repeat
    ang0 = np.array([polymer.backbone[k].backbone_angle for k in range(P)], dtype=float)
    frame = repeat_chains(polymer, "TT", np.full((1, P), 180.0), ang0)[0]
    tt = repeat_chains(polymer, "TT", d["tors"][None], d["angs"], align_to=frame)[0]
    nc = 2 * nx * ny
    pk = pack_mod.CrystalPacker(tt, n_chains=nc, **KW)
    st, flips = MC.supercell_state(d["params"], tt.c, d["tors"], d["angs"], nx, ny)
    model = MC.MixedModel(polymer, "TT", pk, flips, [tt] * nc, angle_stiffness=resolve_angle_stiffness(None, None),
                          pen_rot=pen_rot, k_c=k_c)
    return model, st.vector(), tt


def summary(model, x) -> dict:
    """Energies of ``x``: physical (lattice + torsion + bend), the closure penalties, and the
    physical energy corrected to first order for the slack the penalties allow."""
    p = model.terms(x)
    pen = float(p["penalty_c"] + p["penalty_rotation"])
    return {"physical": float(p["physical"]), "penalty_c": float(p["penalty_c"]), "penalty_rotation": float(p["penalty_rotation"]),
            "corrected": float(p["physical"]) + 2.0 * pen, "lattice": float(p["lattice"]), "torsion": float(p["torsion"]), "bend": float(p["bend"]),
            "c_mismatch_max_A": float(np.abs(p["c_mismatch"]).max()), "rotation_error_max_deg": float(p["rotation_error"].max()),
            "grad_max": float(p["grad_max"])}


def relax_stage(model, x, free, label: str, maxiter: int, log) -> tuple:
    t0 = time.time()
    x1, info = model.relax(x, free=free, maxiter=maxiter)
    s = summary(model, x1)
    s.update(iterations=info["iterations"], evaluations=info["evaluations"], message=info["message"], seconds=time.time() - t0)
    log(f"   {label:34s}: E {s['corrected']:.4f} (physical {s['physical']:.4f}, penalties c {s['penalty_c']:.4f} rot {s['penalty_rotation']:.4f}); "
        f"c mismatch {s['c_mismatch_max_A']:.4f} A, rotation {s['rotation_error_max_deg']:.3f} deg; {info['iterations']} it, {s['seconds']:.0f} s")
    return x1, s


def free_sets(model, kink: int = 0, anchor: int | None = None) -> dict:
    """The variables each stage frees.  ``anchor`` is a host chain whose position stays fixed, which
    removes the cell's rigid translation (default the last chain)."""
    nc = model.nc
    anchor = nc - 1 if anchor is None else anchor
    host = [s for s in range(nc) if s != kink]
    moving_host = [s for s in host if s != anchor]
    kink_int, kink_place = model.index("internal", kink), model.index("place", kink)
    host_place = np.concatenate([model.index("place", moving_host), model.index("phi", anchor)])
    host_int = model.index("internal", host)
    return {"fit": kink_int,
            "rigid_host": np.concatenate([kink_int, kink_place]),
            "host_rigid_bodies": np.concatenate([kink_int, kink_place, host_place]),
            "pinned_cell": np.concatenate([kink_int, kink_place, host_place, host_int]),
            "c_free": np.concatenate([kink_int, kink_place, host_place, host_int, model.index("c")]),
            "cell_free": np.concatenate([kink_int, kink_place, host_place, host_int, model.index("c"), model.index("lattice")])}


def kink_shape(model, x, kink: int, offset: int) -> dict:
    st = model.state(x)
    P = model.P
    t = st.tors[kink]
    w = lambda a: float((a + 180.0) % 360.0 - 180.0)  # noqa: E731
    return {"kink_torsions_deg": [w(t[(offset + k) % P]) for k in range(-2, 5)],
            "largest_trans_deviation_elsewhere_deg": float(max(abs(abs(w(t[(offset + k) % P])) - 180.0) for k in range(5, P - 2))),
            "backbone_angle_range_deg": [float(st.angs[kink].min()), float(st.angs[kink].max())],
            "host_torsion_deviation_max_deg": float(max(np.abs(np.abs((st.tors[s] + 180.0) % 360.0 - 180.0) - 180.0).max() for s in range(model.nc) if s != kink))}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--row", default="vdf10")
    ap.add_argument("--nx", type=int, default=2)
    ap.add_argument("--ny", type=int, default=1)
    ap.add_argument("--kink-at", type=int, default=0, help="bond offset of the kink pair along the repeat")
    ap.add_argument("--mirror-kink", action="store_true", help="G- T G+ instead of G+ T G-")
    ap.add_argument("--preset", default="pvdf-dft-valence-belt4")
    ap.add_argument("--k-c", type=float, default=2000.0, help="penalty tying each chain's repeat to the cell's, kcal/mol/A^2")
    ap.add_argument("--pen-rot", type=float, default=5.0)
    ap.add_argument("--maxiter", type=int, default=400)
    ap.add_argument("--stages", default="rigid_host,host_rigid_bodies,pinned_cell,c_free,cell_free")
    ap.add_argument("--c-scan", default="", help="comma list of changes of the repeat c (A) at which the pinned-cell stage is repeated for both cells, e.g. -0.6,-0.4,-0.2")
    ap.add_argument("--no-tilts", action="store_true")
    ap.add_argument("--out", default="docs/predictions/dilute_kink_belt4.json")
    args = ap.parse_args()
    out_path = Path(args.out)
    results = json.load(open(out_path)) if out_path.exists() else {}
    results.setdefault("rows", {})
    results.update(preset=args.preset, energy="kcal/mol per kink pair: the mixed cell with one kinked chain minus the all-trans supercell relaxed under the same freedoms; "
                   "'corrected' adds twice the closure penalties (first order in the slack they allow)")
    key = f"{args.row}@{args.kink_at}{'~' if args.mirror_kink else ''} {args.nx}x{args.ny}"
    letters, note = sequence(args.row)
    log = lambda s: print(s, flush=True)  # noqa: E731

    def save():
        with open(out_path, "w") as fh:
            json.dump(results, fh, indent=1)

    with PRESET_PARAMETERS[args.preset].applied():
        polymer = C.make_polymer(args.row, letters)
        refs = results.setdefault("all_trans_references", {})
        d = all_trans_reference(polymer, args.row, refs.get(args.row))
        refs[args.row] = {k: (v.tolist() if isinstance(v, np.ndarray) else v) for k, v in d.items()}
        save()
        model, x_tt, tt = make_model(polymer, d, args.nx, args.ny, args.k_c, args.pen_rot)
        nc, n_mon = model.nc, tt.n_monomers
        e_cell = d["E_cell"] + d["bend_per_monomer"] * 2 * n_mon
        entry = results["rows"].setdefault(key, {})
        entry.update(row=args.row, letters=letters, note=note, supercell=[args.nx, args.ny], chains=nc, sites=int(model.pk.N), kink_offset_bonds=args.kink_at,
                     kink_mirrored=bool(args.mirror_kink), k_c=args.k_c, pen_rot=args.pen_rot,
                     all_trans_cell={"params": d["params"].tolist(), "c": d["c"], "energy_per_monomer": e_cell / (2 * n_mon)})
        log(f"{key}: {nc} chains, {model.pk.N} sites, {model.n_var} variables; all-trans cell {np.round(d['params'], 4).tolist()} c {d['c']:.4f}, {e_cell / (2 * n_mon):.5f} kcal/mol per monomer")
        s0 = summary(model, x_tt)
        log(f"   tiling check: physical energy per two chains {s0['physical'] * 2 / nc:.6f} against the two-chain cell's {e_cell:.6f}")
        if abs(s0["physical"] * 2 / nc - e_cell) > 1e-6:
            raise RuntimeError("the tiled all-trans supercell does not reproduce the two-chain cell")
        sets = free_sets(model)
        stages = [s for s in args.stages.split(",") if s]

        # --- the all-trans supercell under the same freedoms (it barely moves; the closure
        # penalties and the independent repeat are what it settles) ---
        ref_e, x_ref = {}, x_tt.copy()
        rng = np.random.default_rng(0)
        x_ref[model.index("phi")] += rng.uniform(-0.3, 0.3, nc)  # off the symmetric point, so an instability of the host would show
        x_ref[model.index("tors")] += rng.uniform(-0.3, 0.3, nc * model.P)
        log("all-trans supercell:")
        for stage in ("pinned_cell", "c_free", "cell_free"):
            if stage == "pinned_cell" or stage in stages:
                x_ref, s = relax_stage(model, x_ref, sets[stage], stage, args.maxiter, log)
                if stage == "pinned_cell":
                    x_ref_pinned = x_ref.copy()
                ref_e[stage] = s
                entry.setdefault("all_trans_supercell", {})[stage] = s
        for stage in ("rigid_host", "host_rigid_bodies"):
            ref_e[stage] = ref_e["pinned_cell"]
        drift = (ref_e["pinned_cell"]["corrected"] - s0["physical"]) * 2 / nc
        log(f"   the relaxed all-trans supercell is {drift:+.4f} kcal/mol per two chains from the tiled cell")
        entry["all_trans_supercell"]["drift_per_two_chains"] = drift
        save()

        # --- the kinked chain ---
        log("kinked chain 0:")
        x = x_tt.copy()
        x[model.index("tors", 0)] = kinked_torsions(d["tors"], args.kink_at, args.mirror_kink)
        start = model.terms(x)
        log(f"   ideal kink in the all-trans site: natural repeat {start['repeats'][0]:.4f} against the cell's {d['c']:.4f}, rotation error {start['rotation_error'][0]:.2f} deg")
        entry["ideal_kink"] = {"natural_repeat_A": float(start["repeats"][0]), "rotation_error_deg": float(start["rotation_error"][0])}
        for kc in (20.0, 200.0, args.k_c):  # fit the chain into the cell's repeat, the tie stiffened in steps
            model.k_c = min(kc, args.k_c)
            x, s = relax_stage(model, x, sets["fit"], f"fit (k_c {model.k_c:g})", args.maxiter, log)
        model.k_c = args.k_c
        entry["stages"] = entry.get("stages", {})
        states = entry.setdefault("states", {})
        for stage in stages:
            x, s = relax_stage(model, x, sets[stage], stage, args.maxiter, log)
            s["formation"] = s["corrected"] - ref_e[stage]["corrected"]
            s["formation_uncorrected"] = s["physical"] - ref_e[stage]["physical"]
            s.update(kink_shape(model, x, 0, args.kink_at))
            st = model.state(x)
            s["cell"] = st.lattice.tolist()
            entry["stages"][stage] = s
            states[stage] = x.tolist()
            log(f"   -> {stage}: kink pair {s['formation']:+.3f} kcal/mol (uncorrected {s['formation_uncorrected']:+.3f}); kink torsions {np.round(s['kink_torsions_deg'], 1).tolist()}")
            if stage == "pinned_cell":
                x_pinned = x.copy()
            save()

        # --- pendant tilts, at the pinned-cell structures ---
        if not args.no_tilts and "pinned_cell" in entry["stages"]:
            t0 = time.time()
            dE_k, tilts_k, groups = MC.tilt_relaxation(model, x_pinned)
            dE_r, _, _ = MC.tilt_relaxation(model, x_ref_pinned)
            entry["tilts"] = {"kinked_cell": dE_k, "all_trans_supercell": dE_r, "difference": dE_k - dE_r,
                              "largest_tilt_deg": float(np.abs(tilts_k).max()), "seconds": time.time() - t0}
            f = entry["stages"]["pinned_cell"]["formation"]
            entry["formation_pinned_cell_with_tilts"] = f + dE_k - dE_r
            log(f"   pendant tilts: kinked cell {dE_k:+.3f}, all-trans supercell {dE_r:+.3f} kcal/mol per cell; kink pair with tilts {f + dE_k - dE_r:+.3f}  ({time.time() - t0:.0f} s)")
            save()

        # --- the repeat scanned: the kinked chain's tension apart from the kink ---
        if args.c_scan and "pinned_cell" in entry["stages"]:
            scan = entry.setdefault("c_scan", {})
            free_k = free_r = sets["pinned_cell"]
            xk, xr = x_pinned.copy(), x_ref_pinned.copy()
            for dc in [0.0] + [float(v) for v in args.c_scan.split(",")]:
                xk[3] = xr[3] = d["c"] + dc
                xk, sk = relax_stage(model, xk, free_k, f"c {dc:+.2f}: kinked", args.maxiter, log)
                xr, sr = relax_stage(model, xr, free_r, f"c {dc:+.2f}: all-trans", args.maxiter, log)
                scan[f"{dc:+.3f}"] = {"kinked": sk["corrected"], "all_trans": sr["corrected"],
                                      "kinked_chain": sk["corrected"] - (nc - 1) / nc * sr["corrected"]}
                save()
            base = ref_e["pinned_cell"]["corrected"] / nc
            vals = {k: v["kinked_chain"] - base for k, v in scan.items()}
            entry["c_scan_formation"] = vals
            entry["formation_tension_free_estimate"] = float(min(vals.values()))
            log(f"   kinked chain against one all-trans chain at the cell's repeat, by change of c: { {k: round(v, 3) for k, v in vals.items()} }")
            save()
    log(f"written {out_path}")


if __name__ == "__main__":
    main()
