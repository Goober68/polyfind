"""The stiffness of a polished two-chain crystal cell: elastic constants from strained, re-relaxed cells.

The fourth part of the objective is the blocked force per area, stiffness times free strain
(docs/CANDIDATE_MATERIALS.md 5.6).  This script gives the crystal's share of the stiffness, for
the cells the candidate ranking rests on: the polished cells of examples/polish_margins.py (every
torsion and backbone angle of both chains free, each chain periodic exactly).

A cell is strained homogeneously (xx, yy, zz and the in-plane shear xy: the four strains its
lattice ``(A, B, gamma, c)`` can carry; x is the polar axis of an all-trans cell, z the chain),
every internal coordinate is relaxed at the fixed cell on the closure surface, and the energy's
second derivatives are taken by central differences over thirty-three cells.  Over the cell's volume they are the elastic
constants ``C_11, C_22, C_33, C_66`` and their couplings (Voigt; 6 is the engineering shear xy),
relaxed-ion and at zero field.  The compliance's diagonal gives the Young's moduli along x, y, z.

What the model does not hold: the bonds are rigid, so the chain-axis constant ``C_33`` is the
bending and twisting of the backbone alone and an upper bound; and the pendant groups are the
builder's (their tilts are not re-relaxed at each strain), so the lateral constants are upper
bounds by a few per cent.  The lateral constants of beta-PVDF are anchored on DFT
(docs/E1_E2_PREDICTIONS.md); nothing here is for the copolymers.

    python examples/crystal_stiffness.py --rows vdf10 --cells TTTT
    python examples/crystal_stiffness.py --rows vdf9-trfe3-an1-rvdf1-14,vdf9-trfe4-an1-14 --cells TTTT,helix
"""
import argparse
import fcntl
import importlib.util
import json
import time
from pathlib import Path

import numpy as np

from polyfind.fitting import PRESET_PARAMETERS

DS = importlib.util.spec_from_file_location("dk", Path(__file__).with_name("dilute_kink.py"))
DK = importlib.util.module_from_spec(DS)
DS.loader.exec_module(DK)
C = DK.C
KCAL_A3_TO_GPA = 4184.0 / 6.02214076e23 / 1e-30 / 1e9  # kcal/mol per A^3 in GPa (6.9477)
VOIGT = ("1", "2", "3", "6")  # xx, yy, zz, the engineering shear xy


def strained_lattice(lattice, e) -> np.ndarray:
    """``(A, B, gamma, c)`` of the cell under the strain ``e = (e1, e2, e3, e6)``: the symmetric
    deformation ``[[1 + e1, e6 / 2], [e6 / 2, 1 + e2]]`` of the lateral cell, and ``c (1 + e3)``.
    The energy depends on the cell's shape alone, so the rotation the deformation leaves is dropped."""
    A, B, gamma, c = (float(v) for v in lattice)
    g = np.radians(gamma)
    F = np.array([[1.0 + e[0], 0.5 * e[3]], [0.5 * e[3], 1.0 + e[1]]])
    a, b = F @ np.array([A, 0.0]), F @ np.array([B * np.cos(g), B * np.sin(g)])
    A2, B2 = float(np.linalg.norm(a)), float(np.linalg.norm(b))
    return np.array([A2, B2, float(np.degrees(np.arccos(a @ b / (A2 * B2)))), c * (1.0 + e[2])])


def strain_set(h: float, h_shear: float) -> list:
    """The origin, plus and minus along each of the four strains, and for each pair of strains its
    four sign combinations: thirty-three cells, so that every second derivative is a central
    difference and the cubic terms of the energy cancel."""
    steps = np.array([h, h, h, h_shear])
    pts = [np.zeros(4)]
    for i in range(4):
        for s in (1.0, -1.0):
            e = np.zeros(4)
            e[i] = s * steps[i]
            pts.append(e)
    for i in range(4):
        for j in range(i + 1, 4):
            for si in (1.0, -1.0):
                for sj in (1.0, -1.0):
                    e = np.zeros(4)
                    e[i], e[j] = si * steps[i], sj * steps[j]
                    pts.append(e)
    return pts


def central_differences(points, energies, steps) -> tuple:
    """``(g, H, asymmetry)`` of ``E = e0 + g . e + e . H . e / 2`` from the cells of
    :func:`strain_set`: the first and second derivatives by central differences, and the largest
    ``|E(+h) - E(-h)|`` along one strain (the cubic term the differences cancel)."""
    E = {tuple(np.sign(p).astype(int)): float(v) for p, v in zip(points, energies)}
    unit = lambda i, s: tuple(s if k == i else 0 for k in range(4))  # noqa: E731
    pair = lambda i, j, si, sj: tuple(si if k == i else (sj if k == j else 0) for k in range(4))  # noqa: E731
    e0 = E[(0, 0, 0, 0)]
    g, H, asym = np.zeros(4), np.zeros((4, 4)), 0.0
    for i in range(4):
        up, dn = E[unit(i, 1)], E[unit(i, -1)]
        g[i] = (up - dn) / (2.0 * steps[i])
        H[i, i] = (up + dn - 2.0 * e0) / steps[i] ** 2
        asym = max(asym, abs(up - dn))
    for i in range(4):
        for j in range(i + 1, 4):
            H[i, j] = H[j, i] = (E[pair(i, j, 1, 1)] + E[pair(i, j, -1, -1)] - E[pair(i, j, 1, -1)] - E[pair(i, j, -1, 1)]) / (4.0 * steps[i] * steps[j])
    return g, H, asym


def stored_cell(row: str, cell: str, margins: dict, references: dict) -> tuple:
    """``(name of the conformer, its stored polished cell)``: ``cell`` is TTTT, ``helix`` (the lowest
    polished helix of the row) or a helix conformer's name."""
    if cell == "TTTT":
        d = references["all_trans_references"][row]
        return "TTTT", d
    name = margins["rows"][row]["tg_conformer"] if cell == "helix" else cell
    return name, margins["references"][f"{row} {name}"]


def stiffness(polymer, d: dict, h: float, h_shear: float, gtol: float, log=print) -> dict:
    d = {k: (np.asarray(v, dtype=float) if isinstance(v, list) else v) for k, v in d.items()}
    model, x0, ch = DK.make_model(polymer, d, 1, 1)
    # one chain's lateral position and axial shift stay put (the cell's translations), the lattice is held at each strain
    held = np.concatenate([model.index("lattice"), model.index("c"), model.index("xy", 1), model.index("z", 1)])
    free = np.setdiff1d(np.arange(model.n_var), held)
    lat0 = model.state(x0).lattice
    volume = float(lat0[0] * lat0[1] * np.sin(np.radians(lat0[2])) * lat0[3])
    pts, energies, worst = strain_set(h, h_shear), [], 0.0
    t0 = time.time()
    for k, e in enumerate(pts):
        x = x0.copy()
        x[:4] = strained_lattice(lat0, e)
        x = model.project(x, model.index("internal"), [0, 1])
        x1, info = model.relax_on_surface(x, free=free, closed=[0, 1], gtol=gtol, maxiter=800)
        energies.append(model.physical_energy(x1))
        worst = max(worst, float(info["gradient_max"]))
        log(f"      strain {k:2d} ({', '.join(f'{v:+.4f}' for v in e)}): {energies[-1] - energies[0]:+.5f} kcal/mol per cell; "
            f"{info['evaluations']} evaluations, gradient {info['gradient_max']:.1e}  ({time.time() - t0:.0f} s)")
    g, H, resid = central_differences(pts, energies, np.array([h, h, h, h_shear]))
    Cm = H / volume * KCAL_A3_TO_GPA
    S = np.linalg.inv(Cm)
    n_mon = 2 * ch.n_monomers
    return {"cell": [float(v) for v in lat0], "volume_A3": volume, "volume_per_monomer_A3": volume / n_mon, "monomers_per_cell": n_mon,
            "cross_section_per_chain_nm2": float(lat0[0] * lat0[1] * np.sin(np.radians(lat0[2])) / 2.0 / 100.0),
            "C_GPa": {f"{VOIGT[i]}{VOIGT[j]}": float(Cm[i, j]) for i in range(4) for j in range(i, 4)},
            "young_GPa": {"x": float(1.0 / S[0, 0]), "y": float(1.0 / S[1, 1]), "z": float(1.0 / S[2, 2])},
            "shear_xy_GPa": float(1.0 / S[3, 3]),
            "poisson": {"zx": float(-S[0, 2] / S[2, 2]), "zy": float(-S[1, 2] / S[2, 2]), "xy": float(-S[1, 0] / S[0, 0])},
            "stable": bool(np.all(np.linalg.eigvalsh(Cm) > 0.0)),
            "residual_stress_GPa": [float(v) for v in g / volume * KCAL_A3_TO_GPA],
            "relaxed_base_below_stored_kcal": float(energies[0] - model.physical_energy(x0)),
            "plus_minus_asymmetry_kcal": resid, "largest_gradient": worst, "strain_step": h, "shear_step": h_shear,
            "energies_kcal": [float(v) for v in energies], "strains": [[float(v) for v in e] for e in pts], "seconds": time.time() - t0}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--rows", required=True, help="comma list of sequence rows with polished cells")
    ap.add_argument("--cells", default="TTTT,helix", help="comma list of TTTT, helix (the row's lowest polished helix) or a helix conformer's name")
    ap.add_argument("--step", type=float, default=0.005, help="the strain step along xx, yy, zz")
    ap.add_argument("--shear-step", type=float, default=0.01, help="the engineering shear step")
    ap.add_argument("--gtol", type=float, default=1e-5)
    ap.add_argument("--preset", default="pvdf-dft-valence-belt4")
    ap.add_argument("--margins", default="docs/predictions/polished_margins_belt4.json")
    ap.add_argument("--references", default="docs/predictions/dilute_kink_belt4.json")
    ap.add_argument("--out", default="docs/predictions/crystal_stiffness_belt4.json")
    args = ap.parse_args()
    out_path = Path(args.out)
    margins, references = json.load(open(args.margins)), json.load(open(args.references))

    def store(key, rec):
        with open(str(out_path) + ".lock", "w") as lock:
            fcntl.flock(lock, fcntl.LOCK_EX)
            disk = json.load(open(out_path)) if out_path.exists() else {}
            disk.update(preset=args.preset,
                        what="elastic constants of polished two-chain cells, relaxed-ion, zero field: energies of thirty-three strained cells with every "
                             "internal coordinate relaxed at the fixed cell, second derivatives by central differences (examples/crystal_stiffness.py)",
                        frame="Voigt 1 xx (the all-trans cell's polar axis), 2 yy (lateral), 3 zz (the chain), 6 the engineering shear xy",
                        limits="rigid bonds: C_33 and the Young's modulus along z hold the backbone's bending and twisting alone and are upper bounds; "
                               "the pendant groups' tilts are not re-relaxed: the lateral constants are upper bounds by a few per cent")
            disk.setdefault("cells", {})[key] = rec
            out_path.write_text(json.dumps(disk, indent=1))

    with PRESET_PARAMETERS[args.preset].applied():
        for row in [r for r in args.rows.split(",") if r]:
            letters, note = DK.sequence(row)
            polymer = C.make_polymer(row, letters)
            for cell in [c for c in args.cells.split(",") if c]:
                try:
                    name, d = stored_cell(row, cell, margins, references)
                except KeyError:
                    print(f"{row} {cell}: no polished cell stored; left out", flush=True)
                    continue
                print(f"{row} {name} ({letters}):", flush=True)
                rec = stiffness(polymer, d, args.step, args.shear_step, args.gtol, log=lambda m: print(m, flush=True))
                rec.update(row=row, letters=letters, conformer=name)
                store(f"{row} {name}", rec)
                c = rec["C_GPa"]
                print(f"   C_11 {c['11']:.2f}, C_22 {c['22']:.2f}, C_12 {c['12']:.2f}, C_66 {c['66']:.2f}, C_33 {c['33']:.1f}, C_13 {c['13']:.2f}, C_23 {c['23']:.2f} GPa; "
                      f"Young's x {rec['young_GPa']['x']:.2f}, y {rec['young_GPa']['y']:.2f}, z {rec['young_GPa']['z']:.1f} GPa; "
                      f"plus/minus asymmetry {rec['plus_minus_asymmetry_kcal']:.1e} kcal/mol, base {rec['relaxed_base_below_stored_kcal']:+.1e}; {rec['seconds']:.0f} s", flush=True)
    print(f"written {out_path}", flush=True)


if __name__ == "__main__":
    main()
