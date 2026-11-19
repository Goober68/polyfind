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

from polyfind.fitting import parameters_of

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
    tag = args.tag or Path(preset).stem.split("-")[-1]
    pol, ref = json.load(open(disk.get("polished", SP.POLISHED))), json.load(open(disk.get("references", SP.REFERENCES)))
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    prm = parameters_of(preset)
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


def write_kink(args):
    """One run of examples/dilute_kink.py as atoms: its all-trans supercell and the kinked cell at each of ``--stages``.

    The kink's energy in another evaluator is then the difference of the two cells there, per kinked
    chain.  A stage that holds the lattice (``flexible_host``) is compared at that lattice: relax the
    atoms alone (``examples/mlip_train.py relax --mask "0 0 0 0 0 0"``), or the kink may leave."""
    from ase.io import write

    spec = importlib.util.spec_from_file_location("dilute_kink", Path(__file__).with_name("dilute_kink.py"))
    DK = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(DK)
    arrays = lambda d: {k: (np.asarray(v, dtype=float) if isinstance(v, list) else v) for k, v in d.items()}  # noqa: E731
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    if args.scan:  # the kink sheet at chosen sites of examples/kink_site_scan.py: every chain kinked in register, in the cell the sheet carries
        scan = json.loads(Path(args.scan).read_text())
        row, sites = args.key.split("@")[0], scan["rows"][args.key]["sites"]
        kinked = arrays(json.loads(Path(args.kinked_references).read_text())["kinked_references"][args.key])
        tag = args.tag or Path(scan["preset"]).stem.split("-")[-1]
        hand = "m" if args.key.endswith("~") else ""
        with parameters_of(scan["preset"]).applied():
            polymer = DK.C.make_polymer(row, scan["rows"][args.key]["letters"])
            model, _, ch = DK.make_model(polymer, kinked, 1, 1)
            for site in [s for s in args.sites.split(",") if s]:
                a = atoms_of(model, np.asarray(sites[site]["state"], dtype=float))
                name = f"sheet_{row}_{tag}_at{site}{hand}"
                a.info.update(family="kinkcell", frame=name, n_monomers=model.nc * ch.n_monomers, kmesh=args.kmesh, conformation=f"kink sheet at bond {site}", kinked_chains=model.nc,
                              source=f"polyfind {scan['preset']}, {args.key}, kink sheet at offset {site} (units {sites[site]['units']}), pendants as built",
                              polyfind_kcal_per_kink_pair=float(sites[site]["formation"]))
                write(out / f"{name}.start.xyz", a, format="extxyz")
                print(f"{name}: {len(a)} atoms, cell {np.round(a.cell.cellpar(), 3).tolist()}, polyfind {sites[site]['formation']:+.3f} kcal/mol per kink pair")
        print(f"written {out}")
        return
    disk = json.loads(Path(args.dilute).read_text())
    entry = disk["rows"][args.key]
    row, (nx, ny), kinked = entry["row"], entry["supercell"], entry["kinked_chains"]
    ref = arrays(disk["all_trans_references"][row])
    tag = args.tag or Path(disk["preset"]).stem.split("-")[-1]
    with parameters_of(disk["preset"]).applied():
        polymer = DK.C.make_polymer(row, entry["letters"])
        model, x_tt, ch = DK.make_model(polymer, ref, nx, ny)
        e_tt = model.physical_energy(x_tt)
        for label, x in [("TT", x_tt)] + [(stage, np.asarray(entry["states"][stage], dtype=float)) for stage in args.stages.split(",") if stage]:
            a = atoms_of(model, x)
            formation = (model.physical_energy(x) - e_tt) / len(kinked)
            name = f"kink_{row}_{tag}_{label}"
            a.info.update(family="kinkcell", frame=name, n_monomers=model.nc * ch.n_monomers, kmesh=args.kmesh, conformation=label, kinked_chains=len(kinked),
                          source=f"polyfind {disk['preset']}, {args.key}, {'the all-trans supercell' if label == 'TT' else 'stage ' + label}, pendants as built",
                          polyfind_kcal_per_kink_pair=float(formation))
            write(out / f"{name}.start.xyz", a, format="extxyz")
            print(f"{name}: {len(a)} atoms, cell {np.round(a.cell.cellpar(), 3).tolist()}, polyfind {formation:+.3f} kcal/mol per kink pair")
    print(f"written {out}")


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
                         "cell_A_deg": [float(v) for v in a.cell.cellpar()], "volume_A3_per_monomer": float(a.get_volume()) / n,
                         "potential_eV_per_monomer": a.info.get("potential_eV_per_monomer")})
        conf = steps[0].info.get("conformation", "")
        table[root] = {"all_trans": conf in ("TTTT", "TT"), "conformation": conf, "source": str(steps[0].info.get("source", steps[0].info.get("note", ""))), "steps": rows}
    out = {"row": args.row, "recipe": "PBE-D3(BJ), 90/360 Ry, the anchor's (examples/mlip_frames.py label)", "energy": "eV per monomer; margins kcal/mol per monomer, all-trans minus helix",
           "lines": table, "margins": {}, "polyfind": {}}
    last = {root: t["steps"][-1] for root, t in table.items()}
    for pair in args.pairs:
        label, tt, tg = pair.split("=")[0], *pair.split("=")[1].split(",")
        first = [table[k]["steps"][0] for k in (tt, tg)]
        labelled_starts = first[0]["frame"] == tt and first[1]["frame"] == tg  # a line whose start was not labelled begins at its first relaxed cell
        both = last[tt]["potential_eV_per_monomer"] is not None and last[tg]["potential_eV_per_monomer"] is not None
        out["margins"][label] = {"all_trans": tt, "helix": tg, "dft": KCAL * (last[tt]["eV_per_monomer"] - last[tg]["eV_per_monomer"]),
                                 "potential": KCAL * (last[tt]["potential_eV_per_monomer"] - last[tg]["potential_eV_per_monomer"]) if both else None,
                                 "dft_at_the_starts": KCAL * (first[0]["eV_per_monomer"] - first[1]["eV_per_monomer"]) if labelled_starts else None,
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
        starts = "the starts not labelled" if m["dft_at_the_starts"] is None else f"at the starts {m['dft_at_the_starts']:+.3f}"
        own = "" if m["potential"] is None else f", the potential's own {m['potential']:+.3f}"
        print(f"  {label}: DFT {m['dft']:+.3f} kcal/mol per monomer{own} ({starts}; largest force left {m['largest_force_left_eV_per_A']:.2f} eV/A)")
    for name, e in out["polyfind"].items():
        print(f"  polyfind {name}: margin {e['margin']:+.3f}, upright {e['margin_upright']:+.3f}; shear gain all-trans {e['all_trans_shear_gain']:+.3f}, helix {e['helix_shear_gain']:+.3f}")
    print(f"written {args.out}")


def potential(args):
    """Rows' cells relaxed in a potential and never labelled (examples/mlip_train.py relax): the
    potential's own energy of each, and all-trans minus each helix placement.  For cells too large
    for DFT here; what the potential's own margin is worth is in docs/MLIP_SECOND_OPINION.md 3.13."""
    from ase.io import read

    out = {"what": f"cells relaxed in the potential {args.name} (examples/mlip_train.py relax) from a preset's cells, the potential's own energies: "
                   "eV per monomer; margins kcal/mol per monomer, all-trans minus helix; no DFT on these cells", "potential": args.name, "rows": {}}
    for row in args.rows:
        cells = {}
        for d in args.frames:
            for f in sorted(Path(d).glob(f"row_{row}_*_in_{args.name}.start.xyz")):
                a = read(f)
                conf = a.info.get("conformation") or f.name[len(f"row_{row}_"):].split("_")[0].replace("p", "'")
                cells[conf] = {"frame": f.name[:-len(".start.xyz")], "eV_per_monomer": float(a.info["potential_eV_per_monomer"]), "cell_A_deg": [float(v) for v in a.cell.cellpar()],
                               "atoms": len(a), "n_monomers": int(a.info["n_monomers"]), "source": str(a.info.get("source", ""))}
        if "TTTT" not in cells or len(cells) < 2:
            print(f"{row}: {sorted(cells)} only, no margin")
            continue
        tt = cells["TTTT"]["eV_per_monomer"]
        minus = {c: KCAL * (tt - v["eV_per_monomer"]) for c, v in cells.items() if c != "TTTT"}
        low = max(minus, key=minus.get)
        out["rows"][row] = {"cells": cells, "all_trans_minus_helix": minus, "margin": minus[low], "helix": low, "helix_placements": len(minus),
                            "strain_along_chain": cells["TTTT"]["cell_A_deg"][2] / cells[low]["cell_A_deg"][2] - 1.0}
        print(f"{row}: margin {minus[low]:+.3f} (helix {low}; " + ", ".join(f"{c} {v:+.3f}" for c, v in minus.items()) + ")")
    Path(args.out).write_text(json.dumps(out, indent=1))
    print(f"written {args.out}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("kink", help="a dilute_kink.py run's all-trans and kinked cells as NAME.start.xyz")
    s.add_argument("--dilute", default="", help="the file examples/dilute_kink.py wrote")
    s.add_argument("--scan", default="", help="instead of --dilute: the file examples/kink_site_scan.py wrote (its key is ROW@0, or ROW@0~ for the other hand)")
    s.add_argument("--kinked-references", default="docs/predictions/kink_migration_belt4.json", help="with --scan: the kinked crystals the scans started from")
    s.add_argument("--sites", default="0", help="with --scan: comma list of the sites (bond offsets) whose cell is written")
    s.add_argument("--key", required=True, help="the run's key there, as 'vdf4-trfe2-6@0 1x1 kinked 0,1'")
    s.add_argument("--stages", default="flexible_host", help="comma list of the stages whose kinked cell is written")
    s.add_argument("--kmesh", required=True)
    s.add_argument("--tag", default="")
    s.add_argument("--out", required=True)
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
    s = sub.add_parser("potential", help="cells relaxed in a potential and not labelled: its own energies and margins")
    s.add_argument("--frames", nargs="+", required=True, help="folders of NAME_in_POTENTIAL.start.xyz")
    s.add_argument("--rows", nargs="+", required=True)
    s.add_argument("--name", default="all_s1", help="the potential's name in the frames' names")
    s.add_argument("--out", required=True)
    args = ap.parse_args()
    {"write": write_cells, "collect": collect, "potential": potential, "kink": write_kink}[args.cmd](args)


if __name__ == "__main__":
    main()
