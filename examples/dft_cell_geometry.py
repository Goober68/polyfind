"""Backbone angles, torsions and bond lengths of a periodic DFT cell, as the fit's angle anchors need them.

    python examples/dft_cell_geometry.py /path/to/result.json            # a Dipole qe_cell record (every record in it)
    python examples/dft_cell_geometry.py cell.xyz                        # an extended XYZ with a Lattice line
    python examples/dft_cell_geometry.py result.json --out angles.json

Bonds are found by distance over the 27 neighbouring images (C-C below 1.75 A, C-H below 1.25,
C-F below 1.55, C-Cl below 1.95); the backbone is the carbons; at every carbon the C-C-C angle is
reported with the carbon's substituents (CH2 or CF2), and along the chain the C-C-C-C torsions.
This is what turned Dipole's beta reference into the anchor 114.63 deg and its alpha start
geometry into 118.23 / 116.87 (docs/CRYSTAL_STAGE.md 5.12), and what the alpha successor's
positions go through when they land.
"""
import argparse
import itertools
import json

import numpy as np

BOND_MAX = {("C", "C"): 1.75, ("C", "H"): 1.25, ("C", "F"): 1.55, ("C", "Cl"): 1.95}
IMAGES = np.array(list(itertools.product((-1, 0, 1), repeat=3)), dtype=float)


def read_extxyz(path):
    lines = open(path).read().splitlines()
    n = int(lines[0])
    comment = lines[1]
    lat = np.array(comment.split('Lattice="')[1].split('"')[0].split(), dtype=float).reshape(3, 3)
    sym, pos = [], []
    for ln in lines[2:2 + n]:
        parts = ln.split()
        sym.append(parts[0])
        pos.append([float(v) for v in parts[1:4]])
    return {"extxyz": {"cell_A": lat.tolist(), "positions_A": pos, "symbols": sym}}


def read_dipole(path):
    d = json.load(open(path))
    recs = d["payload"]["records"] if "payload" in d else d.get("records", {})
    out = {}
    for name, r in recs.items():
        ref = r.get("reference") or r
        if "cell_A" in ref and "positions_A" in ref:
            out[name] = {"cell_A": ref["cell_A"], "positions_A": ref["positions_A"], "symbols": ref["symbols"]}
    if not out and "cell_A" in d:
        out["cell"] = {"cell_A": d["cell_A"], "positions_A": d["positions_A"], "symbols": d["symbols"]}
    return out


def bonds(cell, pos, sym):
    """(i, j, shift, d) for every bonded pair once, the shift (3,) the lattice translation of j."""
    cell = np.asarray(cell, dtype=float)
    pos = np.asarray(pos, dtype=float)
    out = []
    for i in range(len(sym)):
        for j in range(i, len(sym)):
            key = tuple(sorted((sym[i], sym[j]), key=lambda e: (e != "C", e)))
            if key not in BOND_MAX:
                continue
            for img in IMAGES:
                if i == j and (not img.any() or tuple(img) < tuple(-img)):
                    continue  # an atom's own image: each +-shift pair once
                d = float(np.linalg.norm(pos[j] + img @ cell - pos[i]))
                if d < BOND_MAX[key]:
                    out.append((i, j, img.copy(), d))
    return out


def analyse(cell, pos, sym):
    cell = np.asarray(cell, dtype=float)
    pos = np.asarray(pos, dtype=float)
    bl = bonds(cell, pos, sym)
    nb = {i: [] for i in range(len(sym))}
    for i, j, img, d in bl:
        nb[i].append((j, img, d))
        nb[j].append((i, -img, d))
    lengths = {}
    for i, j, img, d in bl:
        key = "-".join(sorted((sym[i], sym[j]), key=lambda e: (e != "C", e)))
        lengths.setdefault(key, []).append(d)
    angles = []
    for c in [i for i, e in enumerate(sym) if e == "C"]:
        cc = [(j, img) for j, img, _ in nb[c] if sym[j] == "C"]
        subs = sorted(sym[j] for j, _, _ in nb[c] if sym[j] != "C")
        if len(cc) == 2:
            (j1, s1), (j2, s2) = cc
            u = pos[j1] + s1 @ cell - pos[c]
            v = pos[j2] + s2 @ cell - pos[c]
            ang = float(np.degrees(np.arccos(np.clip(u @ v / np.linalg.norm(u) / np.linalg.norm(v), -1, 1))))
            angles.append({"carbon": c, "substituents": "".join(subs), "angle_deg": ang})
    torsions = []
    for i, j, img, d in bl:
        if sym[i] != "C" or sym[j] != "C":
            continue
        pi, pj = pos[i], pos[j] + img @ cell
        for k, sk, _ in nb[i]:
            if sym[k] != "C":
                continue
            pk = pos[k] + sk @ cell
            if np.linalg.norm(pk - pj) < 1e-6:
                continue
            for l, sl, _ in nb[j]:
                if sym[l] != "C":
                    continue
                pl = pos[l] + (sl + img) @ cell
                if np.linalg.norm(pl - pi) < 1e-6:
                    continue
                b1, b2, b3 = pi - pk, pj - pi, pl - pj
                n1, n2 = np.cross(b1, b2), np.cross(b2, b3)
                m1 = np.cross(n1, b2 / np.linalg.norm(b2))
                tor = float(np.degrees(np.arctan2(m1 @ n2, n1 @ n2)))
                torsions.append({"atoms": [int(k), int(i), int(j), int(l)], "torsion_deg": tor})
    lengths = {k: {"mean": float(np.mean(v)), "min": float(np.min(v)), "max": float(np.max(v)), "n": len(v)} for k, v in lengths.items()}
    a, b, c = (float(np.linalg.norm(cell[i])) for i in range(3))
    gamma = float(np.degrees(np.arccos(cell[0] @ cell[1] / a / b)))
    return {"cell": {"a": a, "b": b, "c": c, "gamma_deg": gamma, "volume_A3": float(abs(np.linalg.det(cell)))},
            "n_atoms": len(sym), "bond_lengths_A": lengths, "backbone_angles": angles,
            "backbone_angle_mean_deg": float(np.mean([x["angle_deg"] for x in angles])) if angles else None,
            "backbone_angle_by_substituents_deg": {s: float(np.mean([x["angle_deg"] for x in angles if x["substituents"] == s]))
                                                   for s in sorted({x["substituents"] for x in angles})},
            "torsions": torsions, "torsion_abs_values_deg": sorted({round(abs(t["torsion_deg"]), 1) for t in torsions})}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("path")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    cells = read_extxyz(args.path) if args.path.endswith(".xyz") else read_dipole(args.path)
    if not cells:
        raise SystemExit("no cell with positions found")
    results = {}
    for name, c in cells.items():
        r = analyse(c["cell_A"], c["positions_A"], c["symbols"])
        results[name] = r
        print(f"{name}: {r['n_atoms']} atoms; cell a {r['cell']['a']:.4f} b {r['cell']['b']:.4f} c {r['cell']['c']:.4f} A gamma {r['cell']['gamma_deg']:.2f}")
        for k, v in r["bond_lengths_A"].items():
            print(f"  {k}: {v['mean']:.4f} A ({v['min']:.4f} to {v['max']:.4f}, {v['n']} bonds)")
        for s, v in r["backbone_angle_by_substituents_deg"].items():
            print(f"  backbone angle at C{s}: {v:.2f} deg")
        print(f"  backbone angle mean {r['backbone_angle_mean_deg']:.2f} deg; torsions |phi| {r['torsion_abs_values_deg']}")
    if args.out:
        with open(args.out, "w") as fh:
            json.dump(results, fh, indent=1)
        print(f"written {args.out}")


if __name__ == "__main__":
    main()
