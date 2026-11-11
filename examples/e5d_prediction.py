"""E5 d: the half-slip of one chain in four with the lattice giving way, as the DFT seat will run it.

E5 c is a rigid lattice: one chain of the doubled beta cell slid by half a monomer, a single
point.  examples/chain_slide.py finds that the lattice around a half-slipped chain gives way,
the neighbours moving aside, and that this more than halves the cost (6.80 to 3.15 kcal/mol per
monomer in the dilute limit).  The cheapest DFT test of that is the cell the DFT seat already
holds for E5 c, ``beta_2a_slide_c2.xyz`` (24 atoms), relaxed with the lateral cell free and the
chain repeat held.  The slip is held by symmetry, as in E5 b: every atom of that cell sits on a
mirror plane normal to c.  The PI dispatched it on 2026-11-11 ("Go on the DFT check").

This script registers polyfind's number for exactly that protocol:

* the chain repeat set to the DFT cell's (``--c-per-monomer``; Dipole's is 2.5807 A) and held;
* the zero: the unslid cell at that repeat with its lateral cell relaxed (the DFT's zero is its
  accepted cell, which is that for the DFT);
* one chain of four at the half-slip, every chain free to turn and move aside, the lateral cell
  (a, b, gamma) free, the other chains' axial shifts held (in the DFT the mirror plane holds them);
* started on the symmetric lattice, as a DFT code that keeps the symmetry it finds starts, and
  from perturbed copies, to show the symmetric start reaches the lowest state;
* one monomer per chain, as in the DFT's 24-atom cell: at a held repeat such a chain has no
  internal freedom in the model (two torsions and two angles against four closure conditions),
  and in the DFT only its bond lengths;
* then with the pendant groups tilt-relaxed, which a DFT ion relaxation does without being asked.

``--period-two`` adds, for information and not for the comparison, the same with two monomers per
chain and every torsion and angle free: a period-doubling the DFT's cell cannot have.

Energies are kcal/mol: per monomer of the 24-atom cell (four monomers; the unit of E5), per
monomer of the moved chain, and per cell in mEh.  Also written: the same at polyfind's own
repeat, and the rigid half-slip at each repeat, whose ratio to the relaxed number is the size of
the effect and can be formed on the DFT side from E5 c and E5 d alone.
"""
import argparse
import hashlib
import importlib.util
import json
import time
from pathlib import Path

import numpy as np

from polyfind.fitting import PRESET_PARAMETERS

SS = importlib.util.spec_from_file_location("cs", Path(__file__).with_name("chain_slide.py"))
CS = importlib.util.module_from_spec(SS)
SS.loader.exec_module(CS)
DK, MC = CS.DK, CS.MC
KCAL_PER_EH = 627.5094740631
REPEAT = 1  # monomers per chain in the model: the DFT cell's


def energies(dE_cell: float) -> dict:
    """``dE_cell`` is per 24-atom cell (one monomer of each of four chains)."""
    return {"dE_cell_kcal": float(dE_cell), "dE_cell_mEh": float(dE_cell / KCAL_PER_EH * 1e3), "per_cell_monomer_kcal": float(dE_cell / 4), "per_moved_monomer_kcal": float(dE_cell)}


def block(polymer, d, log, seeds, maxiter: int, n_mon: int = REPEAT, levels=("shift",)) -> dict:
    """Everything at the chain repeat of the two-chain cell ``d`` (``n_mon`` monomers per chain)."""
    model, x0, _ = DK.make_model(polymer, d, 2, 1)
    cm = float(d["state"][3]) / n_mon
    dist = CS.distances_from(model, x0, 0)
    anchor = int(np.argmax(dist))
    iz = model.index("z", 0)[0]
    e0 = model.physical_energy(x0)
    xs = x0.copy()
    xs[iz] += 0.5 * cm
    rigid = (model.physical_energy(xs) - e0) / n_mon
    s0 = model.state(x0)
    out = {"chain_repeat_A": cm, "monomers_per_chain_in_the_model": n_mon,
           "zero_cell": {"a": float(s0.lattice[0]) / 2, "b": float(s0.lattice[1]), "gamma_deg": float(s0.lattice[2])},
           "zero_energy_per_monomer_kcal": float(e0 / (4 * n_mon)), "rigid_half_slip": energies(rigid)}
    log(f"   repeat {cm:.5f} A: unslid two-chain cell a {s0.lattice[0] / 2:.5f}, b {s0.lattice[1]:.5f}, gamma {s0.lattice[2]:.3f}; rigid half-slip {rigid / 4:+.4f} per monomer of the cell ({rigid:+.3f} per moved monomer)")
    carried, xb, last = [], xs, levels[-1]
    for level in levels:
        xb, r = CS.held_slip(model, x0, xs, level, anchor, seeds, maxiter, n_mon, lambda s: None, extra_starts=carried, cell=True)
        sb = model.state(xb)
        change = {"a": float(sb.lattice[0] / s0.lattice[0] - 1), "b": float(sb.lattice[1] / s0.lattice[1] - 1), "gamma_deg": float(sb.lattice[2] - s0.lattice[2])}
        frac = sb.place[:, :2] - s0.place[:, :2]
        turn = DK.wrap(sb.place[:, 2] - s0.place[:, 2])
        e = energies(r["barrier"])
        sym = r["starts"]["as started"]
        out[level] = {**e, "ratio_to_rigid": float(r["barrier"] / rigid), "cell_change": change,
                      "relaxed_cell": {"a": float(sb.lattice[0]) / 2, "b": float(sb.lattice[1]), "gamma_deg": float(sb.lattice[2])},
                      "symmetric_start_per_moved_monomer": float(sym), "lowest_start": r["from"], "starts_per_moved_monomer": r["starts"],
                      "symmetric_start_is_lowest": bool(sym - r["barrier"] < 1e-4),
                      "chain_turns_deg": [float(v) for v in turn], "chain_fractional_shifts": [[float(a), float(b)] for a, b in frac],
                      "largest_torsion_off_trans_deg": r["largest_torsion_off_trans_deg"], "axial_force_on_slid_chain": r["axial_force_on_slid_chain"],
                      "largest_axial_force_on_another": r["largest_axial_force_on_another"], "message": r["message"]}
        log(f"      {level:8s}: {e['per_cell_monomer_kcal']:+.4f} per monomer of the cell ({e['per_moved_monomer_kcal']:+.4f} per moved monomer, {e['dE_cell_mEh']:+.3f} mEh per cell; "
            f"{r['barrier'] / rigid:.3f} of rigid); a {change['a']:+.2%}, b {change['b']:+.2%}, gamma {change['gamma_deg']:+.3f} deg; symmetric start {sym:+.4f}, lowest of "
            f"{len(r['starts'])} starts {r['barrier']:+.4f}; turns {np.round(turn, 2).tolist()} deg; torsions {r['largest_torsion_off_trans_deg']:.2f} deg off trans")
        carried = [xb]
    dE, tilts, _ = MC.tilt_relaxation(model, xb)
    dE0, _, _ = MC.tilt_relaxation(model, x0)
    gain = (dE - dE0) / n_mon
    # the same from tilted starts: untilted groups are a minimum only if these return to them
    broken = [MC.tilt_relaxation(model, xb, perturb=2.0, seed=k) for k in range(2)]
    gain_b = min((b[0] - dE0) / n_mon for b in broken)
    out["with_pendant_tilts"] = {**energies(out[last]["dE_cell_kcal"] + gain), "tilt_gain_per_cell_kcal": float(gain), "largest_tilt_deg": float(np.abs(tilts).max()),
                                 "ratio_to_rigid": float((out[last]["dE_cell_kcal"] + gain) / rigid),
                                 "from_tilted_starts": {"tilt_gain_per_cell_kcal": float(gain_b), "largest_tilt_deg": float(max(np.abs(b[1]).max() for b in broken)),
                                                        "per_cell_monomer_kcal": float((out[last]["dE_cell_kcal"] + gain_b) / 4)}}
    log(f"      with pendant tilts: {out['with_pendant_tilts']['per_cell_monomer_kcal']:+.4f} per monomer of the cell (tilts {gain:+.4f} per cell, largest {np.abs(tilts).max():.2f} deg; "
        f"{out['with_pendant_tilts']['ratio_to_rigid']:.3f} of rigid); from tilted starts {out['with_pendant_tilts']['from_tilted_starts']['per_cell_monomer_kcal']:+.4f} "
        f"(tilts {gain_b:+.4f}, largest {out['with_pendant_tilts']['from_tilted_starts']['largest_tilt_deg']:.2f} deg)")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--c-per-monomer", type=float, default=2.5807377932516, help="the DFT cell's chain repeat, A per monomer (Dipole's pvdf_beta_ref_v1)")
    ap.add_argument("--seeds", type=int, default=3)
    ap.add_argument("--period-two", action="store_true", help="also two monomers per chain with every torsion and angle free (information only)")
    ap.add_argument("--maxiter", type=int, default=600)
    ap.add_argument("--preset", default="pvdf-dft-valence-belt4")
    ap.add_argument("--references", default="docs/predictions/chain_slide_belt4.json", help="the stored short-repeat cells of examples/chain_slide.py")
    ap.add_argument("--structures", default="docs/predictions/e5_structures/e5_structures.json")
    ap.add_argument("--out", default="docs/predictions/e5d_polyfind_belt4.json")
    args = ap.parse_args()
    log = lambda s: print(s, flush=True)  # noqa: E731
    t0 = time.time()
    idx = json.load(open(args.structures))
    start = idx["structures"]["beta_2a_slide_c2"]
    result = {"schema": "polyfind_e5d_prediction_v1",
              "potential": f"{args.preset}: induced dipoles, Ewald, force-shifted LJ, fluorine belt sites, PBE-D3 chain geometry, rigid bonds, 0 K",
              "protocol": "one chain of the four-chain cell (the beta cell doubled along its short axis) at the half-monomer slip; the lateral cell (a, b, gamma) and every chain's "
                          "lateral position and setting angle free; the chain repeat held; the other chains' axial shifts held; zero: the unslid cell at the same repeat with its "
                          "lateral cell relaxed",
              "dft_start": {"file": start["file"], "sha256": start["sha256"], "atoms": start["atoms"], "mesh": start["mesh"], "zero": "beta_2a_reference (or twice the accepted energy)"},
              "energy": "kcal/mol; per_cell_monomer: over the four monomers of the 24-atom cell (E5's unit); per_moved_monomer: per monomer of the slid chain; dE_cell: per 24-atom cell",
              "model_cell": f"four chains of {REPEAT} VDF monomer, as in the DFT's 24-atom cell"}
    refs = json.load(open(args.references))["references"]
    with PRESET_PARAMETERS[args.preset].applied():
        polymer, name, d = CS.reference("V" * REPEAT, refs[f"vdf{REPEAT}"])
        log(f"E5 d, {name}, at the DFT seat's repeat:")
        d_dft, info = CS.at_repeat(polymer, d, args.c_per_monomer, REPEAT, args.maxiter)
        result["dft_repeat"] = block(polymer, d_dft, log, list(range(args.seeds)), args.maxiter)
        log("at polyfind's own repeat:")
        result["own_repeat"] = block(polymer, d, log, list(range(args.seeds)), args.maxiter)
        if args.period_two:
            polymer2, name2, d2 = CS.reference("VV", refs["vdf2"])
            log(f"{name2}, every torsion and angle free, at the DFT seat's repeat (information: the DFT's cell cannot double its period):")
            d2_dft, info = CS.at_repeat(polymer2, d2, args.c_per_monomer, 2, args.maxiter)
            result["period_two_flexible_dft_repeat"] = block(polymer2, d2_dft, log, [0], args.maxiter, n_mon=2, levels=("shift", "flexible"))
    reg = result["dft_repeat"]
    result["registered"] = {"what": "the number to compare with the DFT's relaxed energy of beta_2a_slide_c2 against its zero: pendant groups tilt-relaxed, at the DFT seat's repeat",
                            **reg["with_pendant_tilts"], "rigid_pendants": {k: reg["shift"][k] for k in ("dE_cell_kcal", "dE_cell_mEh", "per_cell_monomer_kcal", "per_moved_monomer_kcal")},
                            "cell_change": reg["shift"]["cell_change"], "rigid_half_slip_at_this_repeat": reg["rigid_half_slip"],
                            "ratio_to_rigid": reg["with_pendant_tilts"]["ratio_to_rigid"]}
    result["acceptance"] = "the sign; the relaxed energy within a factor 1.5; a and b both opening; the relaxed number under two thirds of the DFT's own rigid one (E5 c at c/2)"
    result["seconds"] = time.time() - t0
    Path(args.out).write_text(json.dumps(result, indent=1) + "\n")
    log(f"registered: {reg['with_pendant_tilts']['per_cell_monomer_kcal']:+.4f} kcal/mol per monomer of the cell ({reg['with_pendant_tilts']['dE_cell_mEh']:+.3f} mEh per cell), "
        f"a {reg['shift']['cell_change']['a']:+.2%}, b {reg['shift']['cell_change']['b']:+.2%}; rigid {reg['rigid_half_slip']['per_cell_monomer_kcal']:+.4f}; "
        f"sha256 of the start {hashlib.sha256((Path(args.structures).parent / start['file']).read_bytes()).hexdigest()[:16]}...")
    log(f"written {args.out} ({time.time() - t0:.0f} s)")


if __name__ == "__main__":
    main()
