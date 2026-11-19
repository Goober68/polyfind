"""The cells a preset gives a row, written as atoms for DFT.

examples/rows_under_preset.py stores each relaxed cell's state.  This writes those cells as
``NAME.start.xyz`` for ``examples/mlip_frames.py label``: the first step of testing a row's
margin against DFT (the cells as polyfind has them; then relaxed with the potential tilted to
DFT's forces, ``examples/mlip_train.py relax --dft-gradient``, and labelled again).  The pendant
groups are as the chain builder puts them: the rigid tilts polyfind relaxes them by are an
energy in the stored record, not a geometry, and DFT's relaxation moves them itself.

    python examples/row_cells_frames.py write --margins docs/predictions/margins_chain1_dft_rows.json --row vdf3-trfe1-4 --kmesh "4 2 2" --out DIR
    python examples/row_cells_frames.py collect --frames DIR --row vdf3-trfe1-4 --pairs free=TT_LINE,TG_LINE --margins chain1=FILE belt4=FILE --out JSON

``collect`` reads the labelled frames back: each line (a start and the cells relaxed on from it)
with its DFT energy and the force left, and the margins between the last frames of chosen lines.
"""
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np

from polyfind.fitting import PRESET_PARAMETERS

SPEC = importlib.util.spec_from_file_location("shear_polish", Path(__file__).with_name("shear_polish.py"))
SP = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(SP)
ATOMS = ("C", "H", "F", "Cl", "N", "O")  # the sites that are atoms; the others are the potential's own (belt sites)


def safe(conf: str) -> str:
    """A conformation's name as part of a file name: TG'TG -> TGpTG."""
    return conf.replace("'", "p")


def atoms_of(model, x):
    """The cell of ``model`` at ``x`` as an ase.Atoms of its real atoms, chain after chain in the builder's order."""
    from ase import Atoms

    latn, _, _, Pn, _ = model.build(np.asarray(x, dtype=float))
    el = np.array([str(e) for e in list(model.pk.elements) * model.nc])
    keep = np.array([e in ATOMS for e in el])
    a = Atoms(symbols=list(el[keep]), positions=np.asarray(Pn, dtype=float)[keep], cell=np.asarray(latn, dtype=float), pbc=True)
    a.wrap()
    return a


def write_cells(args):
    from ase.io import write

    disk = json.loads(Path(args.margins).read_text())
    entry = disk["rows"][args.row]
    preset = disk["preset"]
    tag = args.tag or preset.split("-")[-1]
    pol, ref = json.load(open(disk.get("polished", SP.POLISHED))), json.load(open(disk.get("references", SP.REFERENCES)))
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    prm = PRESET_PARAMETERS[preset]
    SP.C.SITES_FOR = prm
    with prm.applied():
        polymer = SP.C.make_polymer(args.row, entry["letters"])
        for conf in [c for c in args.confs.split(",") if c] or list(entry["cells"]):
            rec = entry["cells"][conf]
            model, _, ch = SP.sheared_model(polymer, SP.stored_cell(args.row, conf, pol, ref), check=False)
            x = np.asarray(rec[f"state_{args.which}"], dtype=float)
            a = atoms_of(model, x)
            energy = rec[args.which] + rec.get(f"tilt_{args.which}", 0.0)
            name = f"row_{args.row}_{safe(conf)}_{tag}"
            a.info.update(family="rowcell", frame=name, n_monomers=model.nc * ch.n_monomers, kmesh=args.kmesh, conformation=conf,
                          source=f"polyfind {preset}, {args.row} {conf}, cell {args.which}, pendants as built",
                          polyfind_kcal_per_monomer=float(energy))
            write(out / f"{name}.start.xyz", a, format="extxyz")
            print(f"{name}: {len(a)} atoms, cell {np.round(a.cell.cellpar(), 3).tolist()}, polyfind {energy:+.4f} kcal/mol per monomer")
    print(f"written {out}; label with examples/mlip_frames.py label --frames {out}")


KCAL = 23.0605  # kcal/mol in one eV
GPA = 160.21766  # GPa in one eV/A^3


def lines_of(frames: list) -> dict:
    """The labelled frames grouped into relaxation lines: a start and the cells relaxed on from it
    (``START_in_POTENTIAL``, ``START_in_POTENTIAL_n2``, ...), each line in the order it was walked."""
    by_name = {a.info["frame"]: a for a in frames}
    lines = {}
    for name in sorted(by_name):
        root = name.split("_in_")[0]
        lines.setdefault(root, []).append(name)
    out = {}
    for root, names in lines.items():
        order = sorted(names, key=lambda n: (n != root, len(n), n))
        out[root] = [by_name[n] for n in order]
    return out


def collect(args):
    """The DFT energies of the labelled lines of one row, the margins they give, polyfind's beside them."""
    from ase.io import read, write

    frames = [read(f) for f in sorted(Path(args.frames).glob("*.dft.xyz"))]
    lines = lines_of(frames)
    table = {}
    for root, steps in lines.items():
        rows = []
        for a in steps:
            n = a.info["n_monomers"]
            stress = np.asarray(a.info["REF_stress"], dtype=float).reshape(3, 3)
            rows.append({"frame": a.info["frame"], "eV_per_monomer": float(a.info["REF_energy"]) / n, "largest_force_eV_per_A": float(np.abs(a.arrays["REF_forces"]).max()),
                         "stress_GPa_xx_yy_zz_yz_xz_xy": [float(GPA * stress[i, j]) for i, j in ((0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1))],
                         "cell_A_deg": [float(v) for v in a.cell.cellpar()], "volume_A3_per_monomer": float(a.get_volume()) / n})
        conf = steps[0].info.get("conformation", "")
        table[root] = {"all_trans": conf in ("TTTT", "TT"), "conformation": conf, "source": str(steps[0].info.get("source", steps[0].info.get("note", ""))), "steps": rows}
    out = {"row": args.row, "recipe": "PBE-D3(BJ), 90/360 Ry, the anchor's (examples/mlip_frames.py label)", "energy": "eV per monomer; margins kcal/mol per monomer, all-trans minus helix",
           "lines": table, "margins": {}, "polyfind": {}}
    last = {root: t["steps"][-1] for root, t in table.items()}
    for pair in args.pairs:
        label, tt, tg = pair.split("=")[0], *pair.split("=")[1].split(",")
        out["margins"][label] = {"all_trans": tt, "helix": tg, "dft": KCAL * (last[tt]["eV_per_monomer"] - last[tg]["eV_per_monomer"]),
                                 "dft_at_the_starts": KCAL * (table[tt]["steps"][0]["eV_per_monomer"] - table[tg]["steps"][0]["eV_per_monomer"]),
                                 "largest_force_left_eV_per_A": max(last[tt]["largest_force_eV_per_A"], last[tg]["largest_force_eV_per_A"])}
    for item in args.margins:
        name, path = item.split("=")
        e = json.loads(Path(path).read_text())["rows"][args.row]
        out["polyfind"][name] = {"margin": e["margin"], "margin_upright": e["margin_upright"], "helix_sheared": e["helix_sheared"], "helix_upright": e["helix_upright"],
                                 "all_trans_shear_gain": e["all_trans_shear_gain"], "helix_shear_gain": e["helix_shear_gain"]}
    Path(args.out).write_text(json.dumps(out, indent=1))
    if args.extxyz:
        write(args.extxyz, [a for steps in lines.values() for a in steps], format="extxyz")
    for root, t in table.items():
        print(f"{root} ({t['conformation']}): " + "; ".join(f"{s['eV_per_monomer']:.5f} eV, {s['largest_force_eV_per_A']:.2f} eV/A" for s in t["steps"]))
    for label, m in out["margins"].items():
        print(f"  {label}: DFT {m['dft']:+.3f} kcal/mol per monomer (at the starts {m['dft_at_the_starts']:+.3f}; largest force left {m['largest_force_left_eV_per_A']:.2f} eV/A)")
    for name, e in out["polyfind"].items():
        print(f"  polyfind {name}: margin {e['margin']:+.3f}, upright {e['margin_upright']:+.3f}; shear gain all-trans {e['all_trans_shear_gain']:+.3f}, helix {e['helix_shear_gain']:+.3f}")
    print(f"written {args.out}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("write", help="a row's stored cells as NAME.start.xyz")
    s.add_argument("--margins", required=True, help="a file of examples/rows_under_preset.py")
    s.add_argument("--row", required=True)
    s.add_argument("--confs", default="", help="comma list; default every cell of the row")
    s.add_argument("--which", default="sheared", choices=("sheared", "upright"), help="the cell free to shear along its chains, or held upright")
    s.add_argument("--kmesh", required=True, help='the k-point mesh of the DFT, e.g. "4 2 2"')
    s.add_argument("--tag", default="", help="a word for the frames' names; default the preset's last part")
    s.add_argument("--out", required=True, help="working directory for the DFT, outside the repository")
    s = sub.add_parser("collect", help="the labelled lines' DFT energies and the margins between them")
    s.add_argument("--frames", required=True)
    s.add_argument("--row", required=True)
    s.add_argument("--pairs", nargs="+", required=True, metavar="LABEL=TT_LINE,TG_LINE", help="margins to take, each between the last frames of two lines")
    s.add_argument("--margins", nargs="*", default=[], metavar="NAME=FILE", help="polyfind's margins of the row (examples/rows_under_preset.py)")
    s.add_argument("--out", required=True)
    s.add_argument("--extxyz", default="", help="store the labelled frames here")
    args = ap.parse_args()
    {"write": write_cells, "collect": collect}[args.cmd](args)


if __name__ == "__main__":
    main()
