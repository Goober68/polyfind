"""Cross-check E5 (the axial registry of beta-PVDF) against Dipole's records, structure by structure.

    python examples/crosscheck_e5.py --dipole <result.json> [<result.json> ...]

Each Dipole record (``payload.records.<case>.reference`` with ``cell_A``, ``positions_A``,
``symbols``, ``energy_Eh``, as examples/dft_cell_geometry.py reads them) is recognised by its
geometry, not by its name: the chains are found from the bonds, the height of each chain's CH2
carbons along c gives its slip in units of c, and the cell says which part it is:

    12 atoms, slip 0                      the zero of (a) and (b)
    12 atoms, slip c/8 ... c/2, the reference cell      part a
    12 atoms, slip c/2, another cell                    part b (relaxed with the slip held)
    24 atoms, one chain of four slipped by 0, c/4, c/2  part c (0: its zero)
    24 atoms, one chain of four at c/2, another lateral cell, the same repeat
                                                        part d (relaxed with the lateral cell free)

A plain table is read as well: ``{"structures": {"beta_a_slide_c2": {"energy_Eh": ...}, ...}}``
with the names of docs/predictions/e5_structures/e5_structures.json (``beta_b_relaxed`` for part b
and ``beta_d_relaxed`` for part d, each with ``cell_A`` if the cell change is to be compared).

Part d (dispatched 2026-11-11) is the c/2 cell of part c relaxed with the lateral cell free and the
repeat held: polyfind's number for it is docs/predictions/e5d_polyfind_belt4.json
(examples/e5d_prediction.py), read with ``--part-d``.  Its zero is part c's; its acceptance, the
sign, a factor 1.5, a and b both opening, and the relaxed energy under two thirds of the DFT's own
rigid one (part c at c/2).

Energies are compared per monomer of the cell (kcal/mol), against polyfind's two registered
sets (docs/predictions/e5_polyfind_belt4.json): ``own_cell``, polyfind's own minimum slid the same
way, and ``dipole_cell``, polyfind in Dipole's cell.  Acceptance (the ask's): the sign, the
ordering along the slide, one chain of four above one sublattice per moved monomer, and
magnitudes within a factor 1.5.
"""
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np

GS = importlib.util.spec_from_file_location("dg", Path(__file__).with_name("dft_cell_geometry.py"))
DG = importlib.util.module_from_spec(GS)
GS.loader.exec_module(DG)

KCAL_PER_EH = 627.5094740631
A_KEYS = {0.125: "c/8", 0.25: "c/4", 0.375: "3c/8", 0.5: "c/2"}
C_KEYS = {0.25: "c/4", 0.5: "c/2"}


def chains_of(cell, pos, sym) -> list:
    n = len(sym)
    parent = list(range(n))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i, j, _, _ in DG.bonds(cell, pos, sym):
        parent[find(i)] = find(j)
    groups = {}
    for i in range(n):
        groups.setdefault(find(i), []).append(i)
    return sorted(groups.values(), key=lambda g: g[0])


def chain_heights(cell, pos, sym) -> np.ndarray:
    """Fractional height along c of each chain's CH2 carbon (modulo one), from the bonds."""
    cell, pos = np.asarray(cell, dtype=float), np.asarray(pos, dtype=float)
    bonded = {}
    for i, j, _, _ in DG.bonds(cell, pos, sym):
        bonded.setdefault(i, []).append(j)
        bonded.setdefault(j, []).append(i)
    frac = pos @ np.linalg.inv(cell)
    out = []
    for g in chains_of(cell, pos, sym):
        ch2 = [i for i in g if sym[i] == "C" and sum(sym[j] == "H" for j in bonded.get(i, [])) == 2]
        if not ch2:
            raise ValueError("a chain without a CH2 carbon")
        out.append(float(np.angle(np.exp(2j * np.pi * frac[ch2, 2]).mean()) / (2 * np.pi)) % 1.0)
    return np.array(out)


def fold(s: float) -> float:
    """A slip in units of c folded to [0, 1/2] (the registry is even and periodic)."""
    s = s % 1.0
    return min(s, 1.0 - s)


def classify(rec: dict, ref_cell, tol_cell: float = 2e-3, tol_slip: float = 0.01):
    """``(part, slide key or None)`` of a record: part is ``zero``, ``a``, ``b``, ``c0``, ``c`` or ``None``."""
    cell, pos, sym = np.asarray(rec["cell_A"], dtype=float), np.asarray(rec["positions_A"], dtype=float), list(rec["symbols"])
    h = chain_heights(cell, pos, sym)
    ref_cell = np.asarray(ref_cell, dtype=float)

    def near(values, v):
        k = min(values, key=lambda q: abs(q - v))
        return k if abs(k - v) < tol_slip else None

    if len(sym) == 12 and len(h) == 2:
        s = fold(h[1] - h[0])
        same_cell = np.allclose(np.abs(cell), np.abs(ref_cell), atol=tol_cell * np.abs(ref_cell).max())
        if s < tol_slip:
            return ("zero", None) if same_cell else (None, None)
        k = near(A_KEYS, s)
        if k is None:
            return None, None
        if same_cell:
            return "a", A_KEYS[k]
        return ("b", "c/2") if k == 0.5 else (None, None)
    if len(sym) == 24 and len(h) == 4:
        ref2 = ref_cell.copy()
        ref2[0] *= 2.0
        same_cell = np.allclose(np.abs(cell), np.abs(ref2), atol=tol_cell * np.abs(ref2).max())
        same_repeat = abs(np.linalg.norm(cell[2]) - np.linalg.norm(ref2[2])) < tol_cell * np.linalg.norm(ref2[2])
        if not same_cell and not same_repeat:
            return None, None
        # the odd chain out: the one whose height differs from the other three
        d = np.array([[fold(h[i] - h[j]) for j in range(4)] for i in range(4)])
        if d.max() < tol_slip:
            return ("c0", None) if same_cell else (None, None)
        odd = int(np.argmax(d.sum(axis=1)))
        rest = [i for i in range(4) if i != odd]
        if max(d[i, j] for i in rest for j in rest) > tol_slip:
            return None, None
        k = near(C_KEYS, float(np.mean([d[odd, i] for i in rest])))
        if k is None:
            return None, None
        if same_cell:
            return "c", C_KEYS[k]
        return ("d", "c/2") if k == 0.5 else (None, None)  # the half-slipped cell with its lateral cell relaxed, the repeat held
    return None, None


def read_records(paths, structures: dict) -> list:
    """Every record of every file as ``{name, energy_Eh, cell_A, positions_A, symbols}``; a plain
    table's entries take their geometry from the structure index."""
    out = []
    for path in paths:
        d = json.load(open(path))
        recs = d["payload"]["records"] if "payload" in d else d.get("records", {})
        for name, r in recs.items():
            ref = r.get("reference") or r
            if all(k in ref for k in ("cell_A", "positions_A", "symbols", "energy_Eh")):
                out.append({"name": name, "energy_Eh": float(ref["energy_Eh"]), "cell_A": ref["cell_A"], "positions_A": ref["positions_A"], "symbols": ref["symbols"]})
        for name, r in d.get("structures", {}).items():
            if "energy_Eh" not in r:
                continue
            init = structures.get(name, {}).get("initial", {})
            entry = {"name": name, "energy_Eh": float(r["energy_Eh"]), "cell_A": r.get("cell_A", init.get("cell_A")), "positions_A": r.get("positions_A", init.get("positions_A")),
                     "symbols": r.get("symbols", init.get("symbols")), "table_part": {"beta_b_relaxed": "b", "beta_d_relaxed": "d"}.get(name)}
            out.append(entry)
    return out


def verdict(pf: float, dp: float) -> str:
    if np.sign(pf) != np.sign(dp):
        return "SIGN DISAGREES"
    ratio = abs(pf) / abs(dp) if dp else float("inf")
    return f"ratio {ratio:.2f}: " + ("within 1.5" if 1 / 1.5 <= ratio <= 1.5 else "OUTSIDE 1.5")


def crosscheck(pred: dict, structures: dict, records: list, log=print, pred_d: dict | None = None) -> dict:
    """Compare; returns ``{part: {key: {dipole, own_cell, dipole_cell, ...}}}`` for a caller or a test.
    ``pred_d`` is the part-d registration (docs/predictions/e5d_polyfind_belt4.json), if there is one."""
    ref_cell = np.diag([pred["dipole_reference"]["cell_A"][0], pred["dipole_reference"]["cell_A"][1], pred["dipole_reference"]["cell_A"][3]])
    e_ref = float(pred["dipole_reference"]["energy_Eh"])
    found = {"a": {}, "c": {}}
    zero, c0, b, dd = None, None, None, None
    for r in records:
        if r.get("table_part") == "b":
            b = r
            continue
        if r.get("table_part") == "d":
            dd = r
            continue
        if r["cell_A"] is None or r["positions_A"] is None:
            log(f"  {r['name']}: no geometry, skipped")
            continue
        part, key = classify(r, ref_cell)
        if part == "zero":
            zero = r
        elif part == "c0":
            c0 = r
        elif part == "b":
            b = r
        elif part == "d":
            dd = r
        elif part in ("a", "c"):
            found[part][key] = r
        else:
            log(f"  {r['name']}: not an E5 structure, skipped")
            continue
        log(f"  {r['name']}: part {part}{' ' + key if key else ''}, E {r['energy_Eh']:.10f} Eh")
    e0 = zero["energy_Eh"] if zero else e_ref
    log(f"zero of a and b: {'the unslid record ' + zero['name'] if zero else 'the accepted reference energy'} ({e0:.10f} Eh)"
        + (f"; differs from the accepted reference by {(e0 - e_ref) * KCAL_PER_EH:+.5f} kcal/mol per cell" if zero else ""))
    out = {"a": {}, "b": {}, "c": {}, "d": {}}
    own, dip = pred["own_cell"], pred.get("dipole_cell", {})
    if found["a"]:
        log("\npart a: one sublattice slid rigidly; kcal/mol per monomer of the cell")
        prev = 0.0
        for key in ("c/8", "c/4", "3c/8", "c/2"):
            if key not in found["a"]:
                continue
            dE = (found["a"][key]["energy_Eh"] - e0) * KCAL_PER_EH / 2
            p_own = own["a"]["slides"][key]["per_cell_monomer_kcal"]
            p_dip = dip["a"]["slides"][key]["per_cell_monomer_kcal"] if dip else float("nan")
            rising = dE > prev
            prev = dE
            out["a"][key] = {"dipole": dE, "own_cell": p_own, "dipole_cell": p_dip, "rising": bool(rising)}
            log(f"  {key:5s} Dipole {dE:+8.3f}   polyfind own cell {p_own:+8.3f} ({verdict(p_own, dE)})   in Dipole's cell {p_dip:+8.3f} ({verdict(p_dip, dE)})"
                + ("" if rising else "   NOT RISING along the slide"))
    if b is not None:
        dE = (b["energy_Eh"] - e0) * KCAL_PER_EH / 2
        pb = own["b"]
        out["b"] = {"dipole": dE, "own_cell": pb["per_cell_monomer_kcal"], "own_cell_symmetry_broken": pb.get("symmetry_broken", {}).get("per_cell_monomer_kcal")}
        log(f"\npart b: c/2 relaxed with the slip held; kcal/mol per monomer of the cell\n  Dipole {dE:+8.3f}   polyfind {pb['per_cell_monomer_kcal']:+8.3f} ({verdict(pb['per_cell_monomer_kcal'], dE)})"
            + (f"   polyfind from a symmetry-broken start {pb['symmetry_broken']['per_cell_monomer_kcal']:+8.3f}" if "symmetry_broken" in pb else ""))
        if b.get("cell_A") is not None:
            cell = np.asarray(b["cell_A"], dtype=float)
            ch = {"a": float(np.linalg.norm(cell[0]) / ref_cell[0, 0] - 1), "b": float(np.linalg.norm(cell[1]) / ref_cell[1, 1] - 1), "c": float(np.linalg.norm(cell[2]) / ref_cell[2, 2] - 1)}
            out["b"]["cell_change"] = ch
            pc = pb["cell_change"]
            log(f"  cell change: Dipole a {ch['a']:+.2%}, b {ch['b']:+.2%}, c {ch['c']:+.2%}   polyfind a {pc['a']:+.2%}, b {pc['b']:+.2%}, c {pc['c']:+.2%}   "
                + ("b opens and a closes on both sides" if ch["b"] > 0 and ch["a"] < 0 else "THE CELL CHANGE DIFFERS IN SIGN"))
    if found["c"]:
        e0c = c0["energy_Eh"] if c0 else 2.0 * e0
        log(f"\npart c: one chain of four slid rigidly; zero: {'the unslid doubled record ' + c0['name'] if c0 else 'twice the zero of a (no doubled reference record)'}"
            + (f", which differs from twice the zero of a by {(c0['energy_Eh'] - 2 * e0) * KCAL_PER_EH:+.5f} kcal/mol per cell" if c0 else ""))
        for key in ("c/4", "c/2"):
            if key not in found["c"]:
                continue
            dE = (found["c"][key]["energy_Eh"] - e0c) * KCAL_PER_EH / 4
            p_own = own["c"]["slides"][key]["per_cell_monomer_kcal"]
            p_dip = dip["c"]["slides"][key]["per_cell_monomer_kcal"] if dip else float("nan")
            out["c"][key] = {"dipole": dE, "own_cell": p_own, "dipole_cell": p_dip}
            log(f"  {key:5s} Dipole {dE:+8.3f} ({4 * dE:+.2f} per monomer of the moved chain)   polyfind own cell {p_own:+8.3f} ({verdict(p_own, dE)})   "
                f"in Dipole's cell {p_dip:+8.3f} ({verdict(p_dip, dE)})")
        if "c/2" in out["c"] and "c/2" in out["a"]:
            ratio_d = (4 * out["c"]["c/2"]["dipole"]) / (2 * out["a"]["c/2"]["dipole"])
            ratio_p = own["anisotropy_c_over_a_per_moved_monomer_at_half"]
            out["anisotropy"] = {"dipole": float(ratio_d), "own_cell": float(ratio_p), "dipole_cell": float(dip.get("anisotropy_c_over_a_per_moved_monomer_at_half", float("nan")))}
            log(f"\none chain of four against one sublattice, per moved monomer at c/2: Dipole {ratio_d:.2f}   polyfind {ratio_p:.2f} (own cell), "
                f"{out['anisotropy']['dipole_cell']:.2f} (Dipole's cell)   {'one chain costs more on both sides' if ratio_d > 1 else 'THE ANISOTROPY IS REVERSED'}")
    if dd is not None:
        e0c = c0["energy_Eh"] if c0 else 2.0 * e0
        dE = (dd["energy_Eh"] - e0c) * KCAL_PER_EH / 4
        out["d"] = {"dipole": dE}
        log(f"\npart d: one chain of four at c/2, the lateral cell relaxed, the repeat held; kcal/mol per monomer of the cell; zero: "
            f"{'the unslid doubled record ' + c0['name'] if c0 else 'twice the zero of a (no doubled reference record)'}")
        if pred_d is None:
            log(f"  Dipole {dE:+8.3f}   (no part-d registration given: --part-d)")
        else:
            reg = pred_d["registered"]
            p_d = reg["per_cell_monomer_kcal"]
            out["d"].update(polyfind=p_d, polyfind_rigid_pendants=reg["rigid_pendants"]["per_cell_monomer_kcal"])
            log(f"  Dipole {dE:+8.3f}   polyfind {p_d:+8.3f} ({verdict(p_d, dE)})")
            if dd.get("cell_A") is not None:
                cell = np.asarray(dd["cell_A"], dtype=float)
                ch = {"a": float(np.linalg.norm(cell[0]) / (2 * ref_cell[0, 0]) - 1), "b": float(np.linalg.norm(cell[1]) / ref_cell[1, 1] - 1),
                      "c": float(np.linalg.norm(cell[2]) / ref_cell[2, 2] - 1)}
                out["d"]["cell_change"] = ch
                pc = reg["cell_change"]
                log(f"  cell change: Dipole a {ch['a']:+.2%}, b {ch['b']:+.2%}, c {ch['c']:+.2%} (held)   polyfind a {pc['a']:+.2%}, b {pc['b']:+.2%}   "
                    + ("a and b open on both sides" if ch["a"] > 0 and ch["b"] > 0 else "THE CELL DOES NOT OPEN IN BOTH DIRECTIONS"))
            if "c/2" in out["c"]:
                ratio_d = dE / out["c"]["c/2"]["dipole"]
                ratio_p = p_d / out["c"]["c/2"]["dipole_cell"] if np.isfinite(out["c"]["c/2"]["dipole_cell"]) else float("nan")
                out["d"]["relaxed_over_rigid"] = {"dipole": float(ratio_d), "polyfind": float(ratio_p)}
                log(f"  relaxed over rigid (part d over part c at c/2): Dipole {ratio_d:.2f}   polyfind {ratio_p:.2f}   "
                    + ("the lattice's give takes more than a third off on both sides" if ratio_d < 2 / 3 else "THE DFT LATTICE GIVES WAY BY LESS THAN A THIRD"))
    if not (found["a"] or found["c"] or b is not None or dd is not None):
        log("no E5 structure among the records")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--polyfind", default="docs/predictions/e5_polyfind_belt4.json")
    ap.add_argument("--structures", default="docs/predictions/e5_structures/e5_structures.json")
    ap.add_argument("--part-d", default="docs/predictions/e5d_polyfind_belt4.json", help="the part-d registration, used if the file exists")
    ap.add_argument("--dipole", nargs="+", required=True, help="Dipole result.json files, or a plain table of energies by structure name")
    args = ap.parse_args()
    pred = json.load(open(args.polyfind))
    structures = json.load(open(args.structures))["structures"] if Path(args.structures).exists() else {}
    pred_d = json.load(open(args.part_d)) if args.part_d and Path(args.part_d).exists() else None
    print(f"E5: polyfind {pred['potential']}")
    crosscheck(pred, structures, read_records(args.dipole, structures), pred_d=pred_d)


if __name__ == "__main__":
    main()
