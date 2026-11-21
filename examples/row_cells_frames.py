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
    if args.sheet_of:  # a stored kinked crystal itself (the packing search's, polished) and the row's all-trans cell: for a row with no site scan (PVDF)
        row = args.sheet_of.split("@")[0]
        kinked = arrays(json.loads(Path(args.kinked_references).read_text())["kinked_references"][args.sheet_of])
        store = json.loads(Path(args.dilute).read_text())
        flat = arrays(store["all_trans_references"][row])
        tag = args.tag or Path(store["preset"]).stem.split("-")[-1]
        with parameters_of(store["preset"]).applied():
            polymer = DK.C.make_polymer(row, DK.sequence(row)[0])
            formation = (float(kinked["energy"]) - float(flat["energy"])) / 2.0
            for label, d, e in (("TT", flat, 0.0), ("ref", kinked, formation)):
                model, x, ch = DK.make_model(polymer, d, 1, 1)
                a = atoms_of(model, x)
                name = f"sheet_{row}_{tag}_{label}"
                a.info.update(family="kinkcell", frame=name, n_monomers=model.nc * ch.n_monomers, kmesh=args.kmesh, conformation="all-trans" if label == "TT" else "kink sheet",
                              kinked_chains=model.nc, source=f"polyfind {store['preset']}, {args.sheet_of}, " + ("the all-trans cell" if label == "TT" else "the kinked crystal") + ", pendants as built",
                              polyfind_kcal_per_kink_pair=float(e))
                write(out / f"{name}.start.xyz", a, format="extxyz")
                print(f"{name}: {len(a)} atoms, cell {np.round(a.cell.cellpar(), 3).tolist()}, polyfind {e:+.3f} kcal/mol per kink pair")
        print(f"written {out}")
        return
    if args.scan:  # the kink sheet at chosen sites of examples/kink_site_scan.py: every chain kinked in register, in the cell the sheet carries
        scan = json.loads(Path(args.scan).read_text())
        row, sites = args.key.split("@")[0], scan["rows"][args.key]["sites"]
        kinked = arrays(json.loads(Path(args.kinked_references).read_text())["kinked_references"][args.key])
        tag = args.tag or Path(scan["preset"]).stem.split("-")[-1]
        hand = "m" if args.key.endswith("~") else ""
        # with --sheared: each site's lowest sheet of examples/kink_sheet_shear.py (the cell free to shear along its chains), in place of the scan's upright one
        sheared = json.loads(Path(args.sheared).read_text())["rows"] if args.sheared else None
        with parameters_of(scan["preset"]).applied():
            polymer = DK.C.make_polymer(row, scan["rows"][args.key]["letters"])
            model, _, ch = SP.sheared_model(polymer, kinked, check=False) if sheared else DK.make_model(polymer, kinked, 1, 1)
            for site in [s for s in args.sites.split(",") if s]:
                state, formation, how, mark = sites[site]["state"], float(sites[site]["formation"]), "", ""
                if sheared:
                    best = min((t for t in sheared[f"{args.key} site {site} {scan['preset']}"]["starts"] if t.get("kink_kept") and "failed" not in t), key=lambda t: t["formation"])
                    state, formation, mark = best["state"], float(best["formation"]), "_sheared"
                    how = f", the cell free to shear: {best['shear_monomers'][0]:+.2f}, {best['shear_monomers'][1]:+.2f} monomers along a, b"
                a = atoms_of(model, np.asarray(state, dtype=float))
                name = f"sheet_{row}_{tag}_at{site}{hand}{mark}"
                a.info.update(family="kinkcell", frame=name, n_monomers=model.nc * ch.n_monomers, kmesh=args.kmesh, conformation=f"kink sheet at bond {site}", kinked_chains=model.nc,
                              source=f"polyfind {scan['preset']}, {args.key}, kink sheet at offset {site} (units {sites[site]['units']}){how}, pendants as built",
                              polyfind_kcal_per_kink_pair=formation)
                write(out / f"{name}.start.xyz", a, format="extxyz")
                print(f"{name}: {len(a)} atoms, cell {np.round(a.cell.cellpar(), 3).tolist()}, polyfind {formation:+.3f} kcal/mol per kink pair")
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


def kink_collect(args):
    """A kink's energy by DFT, by the potential that relaxed the cells and by polyfind: each ``--pairs NAME=ALL_TRANS,KINKED`` is
    two labelled frames (``*.dft.xyz``, paths) of one lattice, the difference per kinked chain."""
    from ase.io import read, write

    spec = importlib.util.spec_from_file_location("chain_split", Path(__file__).with_name("chain_split.py"))
    CS = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(CS)
    out = {"what": "a kink pair's energy (T G+ T G- T on every chain of the cell, in register), kcal/mol per kink pair: the kinked cell minus the all-trans cell, per "
                   "kinked chain. own_cells false: the two share a lattice (a kink at a held lattice); true: each was relaxed with its cell free (a kink sheet in its "
                   "own cell). DFT is a single point on each cell as the potential relaxed it (the force and the stress left are stored); polyfind is the preset's "
                   "own value for the structure in source. steps: the kinked cell relaxed on in the potential tilted to DFT's forces and stress "
                   "(examples/mlip_train.py relax --dft-gradient) and labelled again; dft_predicted is what that step expected of DFT",
           "recipe": "PBE-D3(BJ), 90/360 Ry, the anchor's (examples/mlip_frames.py label)", "lattices": {}}
    frames = []
    steps = {item.split("=")[0]: item.split("=")[1].split(",") for item in args.steps}
    for item in args.pairs:
        name, files = item.split("=")
        tt, kinked = (read(f) for f in files.split(","))
        later = []
        n, chains = int(tt.info["n_monomers"]), int(kinked.info["kinked_chains"])
        per = lambda key: (float(kinked.info[key]) * (n if "per_monomer" in key else 1) - float(tt.info[key]) * (n if "per_monomer" in key else 1)) * KCAL / chains  # noqa: E731
        torsions = [[float(v) for v in CS.backbone_torsions(CS.lone_chain(kinked, idx, P, T))["dihedrals"]] for idx, P, T in CS.chains_of(kinked)]
        rec = {"all_trans": str(tt.info["frame"]), "kinked": str(kinked.info["frame"]), "atoms": len(tt), "n_monomers": n, "kinked_chains": chains,
               "cell_A_deg": [float(v) for v in tt.cell.cellpar()], "volume_A3_per_monomer": float(tt.get_volume() / n),
               "dft": per("REF_energy"), "potential": per("potential_eV_per_monomer"),
               "largest_force_left_eV_per_A": float(max(np.abs(tt.arrays["REF_forces"]).max(), np.abs(kinked.arrays["REF_forces"]).max())),
               "potential_minus_dft_meV_per_monomer": {"all_trans": 1e3 * (float(tt.info["potential_eV_per_monomer"]) - float(tt.info["REF_energy"]) / n),
                                                       "kinked": 1e3 * (float(kinked.info["potential_eV_per_monomer"]) - float(kinked.info["REF_energy"]) / n)},
               "kink_dihedrals_deg": [t[:3] for t in torsions], "source": str(kinked.info.get("source", ""))}
        rec["own_cells"] = bool(np.abs(kinked.cell.array - tt.cell.array).max() > 1e-6)
        if rec["own_cells"]:
            rec["kinked_cell_A_deg"] = [float(v) for v in kinked.cell.cellpar()]
        if "REF_stress" in tt.info and "REF_stress" in kinked.info:
            rec["largest_stress_left_GPa"] = {"all_trans": float(np.abs(tt.info["REF_stress"]).max() * GPA), "kinked": float(np.abs(kinked.info["REF_stress"]).max() * GPA)}
        if "polyfind_kcal_per_kink_pair" in kinked.info and "carried to the lattice" not in rec["source"] and name not in args.moved:
            rec["polyfind"] = float(kinked.info["polyfind_kcal_per_kink_pair"])
        for k, f in enumerate(steps.get(name, []), start=1):  # the kinked cell after each Newton step toward DFT's own minimum
            b = read(f)
            moved = (b.get_scaled_positions(wrap=False) - kinked.get_scaled_positions(wrap=False)) @ b.cell.array
            turns = [[float(v) for v in CS.backbone_torsions(CS.lone_chain(b, idx, P, T))["dihedrals"]][:3] for idx, P, T in CS.chains_of(b)]
            rec.setdefault("steps", []).append({
                "frame": str(b.info["frame"]), "from": str(b.info.get("start_frame", "")), "cell_A_deg": [float(v) for v in b.cell.cellpar()],
                "dft": (float(b.info["REF_energy"]) - float(tt.info["REF_energy"])) * KCAL / chains,
                "dft_predicted": (float(b.info["dft_predicted_eV_per_monomer"]) * n - float(tt.info["REF_energy"])) * KCAL / chains,
                "largest_force_left_eV_per_A": float(np.abs(b.arrays["REF_forces"]).max()), "largest_stress_left_GPa": float(np.abs(b.info["REF_stress"]).max() * GPA),
                "largest_move_from_the_first_cell_A": float(np.linalg.norm(moved, axis=1).max()), "kink_dihedrals_deg": turns})
            later.append(b)
        out["lattices"][name] = rec
        frames += [tt, kinked, *later]
        print(f"{name}: cell {np.round(rec['cell_A_deg'], 3).tolist()}" + (f" (kinked {np.round(rec['kinked_cell_A_deg'], 3).tolist()})" if rec["own_cells"] else "")
              + f", {rec['volume_A3_per_monomer']:.2f} A3 a monomer; kink pair DFT {rec['dft']:+.3f}, the potential {rec['potential']:+.3f}"
              + (f", polyfind {rec['polyfind']:+.3f}" if "polyfind" in rec else "") + f"; force left {rec['largest_force_left_eV_per_A']:.2f} eV/A; kink at "
              + ", ".join(f"{v:.0f}" for v in rec["kink_dihedrals_deg"][0]))
        for k, st in enumerate(rec.get("steps", []), start=1):
            print(f"  step {k}: DFT {st['dft']:+.3f} (expected {st['dft_predicted']:+.3f}); force left {st['largest_force_left_eV_per_A']:.2f} eV/A, "
                  f"stress {st['largest_stress_left_GPa']:.2f} GPa; kink at " + ", ".join(f"{v:.0f}" for v in st["kink_dihedrals_deg"][0]))
    for item in args.also:  # numbers with no DFT of their own: NAME=polyfind kcal per kink pair, as chain1 at its own lattice
        name, value = item.split("=")
        out.setdefault("polyfind_only", {})[name] = float(value)
    Path(args.out).write_text(json.dumps(out, indent=1))
    if args.extxyz:
        seen, keep = set(), []
        for a in frames:
            if a.info["frame"] not in seen:
                seen.add(a.info["frame"])
                keep.append(a)
        write(args.extxyz, keep, format="extxyz")
    print(f"written {args.out}")


def shear_along_chains(atoms, monomer: float = 2.56) -> list:
    """How far the neighbours along the first two lattice vectors sit along the chain axis (the third), in monomers: the same
    number however the cell lies in space (examples/cell_tilt.py)."""
    cell = atoms.cell.array
    axis = cell[2] / float((cell[2] @ cell[2]) ** 0.5)
    return [float(cell[0] @ axis) / monomer, float(cell[1] @ axis) / monomer]


def sheets_potential(args):
    """Kink sheets relaxed in a potential and never labelled (``sheet_ROW_TAG_atSITE_in_NAME.start.xyz``), against each row's
    all-trans cell in the same potential (``examples/row_cells_frames.py potential``): kcal/mol per kink pair."""
    from ase.io import read

    spec = importlib.util.spec_from_file_location("chain_split", Path(__file__).with_name("chain_split.py"))
    CS = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(CS)
    all_trans = json.loads(Path(args.all_trans).read_text())["rows"]
    out = {"what": "the kink sheet (a kink pair on every chain, in register, at one site of the sequence) relaxed in the potential with the cell free, minus the row's all-trans "
                   "cell there, kcal/mol per kink pair; polyfind is the tables' value for the same site (examples/kink_site_scan.py). No DFT on these cells. "
                   "kink_kept: every chain still has two gauche bonds of opposite sense about a trans one", "potential": args.name, "all_trans": Path(args.all_trans).name, "rows": {}}
    # --pairs: each sheet against an all-trans cell named with it (both relaxed to the same stop), in place of the folder's sheets against the stored rows
    pairs = [(item.split("=")[0], *item.split("=")[1].split(",")) for item in args.pairs]
    if pairs:
        out["what"] = out["what"].replace("minus the row's all-trans cell there", "minus the all-trans cell named with it, relaxed there to the same stop")
        out["all_trans"] = "the first file of each pair"
    for key, flat, f in pairs or [(None, None, f) for f in sorted(Path(args.frames).glob(f"sheet_*_in_{args.name}.start.xyz"))]:
        a = read(f)
        row = str(a.info["frame"]).split("_")[1]
        n, chains = int(a.info["n_monomers"]), int(a.info["kinked_chains"])
        torsions = [[float(v) for v in CS.backbone_torsions(CS.lone_chain(a, idx, P, T))["dihedrals"]] for idx, P, T in CS.chains_of(a)]
        gauche = [[(k, t) for k, t in enumerate(ch) if abs(t) < 120.0] for ch in torsions]
        kept = all(len(g) == 2 and g[0][1] * g[1][1] < 0 and (g[1][0] - g[0][0]) in (2, len(ch) - 2) for g, ch in zip(gauche, torsions))
        e_flat = float(all_trans[row]["cells"]["TTTT"]["eV_per_monomer"]) if flat is None else float(read(flat).info["potential_eV_per_monomer"])
        rec = {"frame": str(a.info["frame"]), "atoms": len(a), "site": str(a.info["frame"]).split("_at")[1].split("_in_")[0].split("_")[0],
               "potential": (float(a.info["potential_eV_per_monomer"]) - e_flat) * n * KCAL / chains,
               "polyfind": float(a.info["polyfind_kcal_per_kink_pair"]), "kink_kept": bool(kept), "gauche_bonds_deg": [[round(t, 1) for _, t in g] for g in gauche],
               "cell_A_deg": [float(v) for v in a.cell.cellpar()], "shear_along_chains_monomers": shear_along_chains(a), "source": str(a.info.get("source", ""))}
        if flat is not None:
            b = read(flat)
            flat_torsions = [[float(v) for v in CS.backbone_torsions(CS.lone_chain(b, idx, P, T))["dihedrals"]] for idx, P, T in CS.chains_of(b)]
            rec.update(all_trans_frame=str(b.info["frame"]), all_trans_eV_per_monomer=e_flat, sheet_eV_per_monomer=float(a.info["potential_eV_per_monomer"]),
                       all_trans_cell_A_deg=[float(v) for v in b.cell.cellpar()], all_trans_is_all_trans=bool(all(abs(t) > 150.0 for ch in flat_torsions for t in ch)),
                       relaxed={"sheet": str(a.info.get("note", "")).split(", relaxed in ")[-1], "all_trans": str(b.info.get("note", "")).split(", relaxed in ")[-1]})
        out["rows"][key or (row if row not in out["rows"] else f"{row}@{rec['site']}")] = rec   # a row's second site is keyed with its site
        print(f"{row:26s} site {rec['site']:4s}: the potential {rec['potential']:+.3f}, polyfind {rec['polyfind']:+.3f} kcal/mol per kink pair; "
              f"{'kink kept' if kept else 'KINK NOT KEPT'}, gauche bonds {rec['gauche_bonds_deg']}")
    # --controls: an all-trans cell relaxed from another start (a tilt of the chain axis, say) against the one the sheets are read
    # against: whether the reference itself has a lower cell where a sheet went. Files named after the second are the legs it
    # passed through on the way, in order: where the cell was and what it cost there, when the start is not where it rests
    for key, flat, other, *legs in [(item.split("=")[0], *item.split("=")[1].split(",")) for item in args.controls]:
        a, b = read(flat), read(other)

        def chain_torsions(cell):
            return [[float(v) for v in CS.backbone_torsions(CS.lone_chain(cell, idx, P, T))["dihedrals"]] for idx, P, T in CS.chains_of(cell)]

        def is_all_trans(cell):
            return bool(all(abs(t) > 150.0 for ch in chain_torsions(cell) for t in ch))

        chains = len(chain_torsions(b))
        rec = {"frame": str(b.info["frame"]), "against": str(a.info["frame"]),
               "potential": (float(b.info["potential_eV_per_monomer"]) - float(a.info["potential_eV_per_monomer"])) * int(b.info["n_monomers"]) * KCAL / chains,
               "is_all_trans": is_all_trans(b), "cell_A_deg": [float(v) for v in b.cell.cellpar()],
               "against_cell_A_deg": [float(v) for v in a.cell.cellpar()], "shear_along_chains_monomers": shear_along_chains(b),
               "against_shear_along_chains_monomers": shear_along_chains(a), "source": str(b.info.get("source", "")),
               "relaxed": {"this": str(b.info.get("note", "")).split(", relaxed in ")[-1], "against": str(a.info.get("note", "")).split(", relaxed in ")[-1]}}
        if legs:
            rec["legs_before"] = [{"frame": str(c.info["frame"]), "potential": (float(c.info["potential_eV_per_monomer"]) - float(a.info["potential_eV_per_monomer"]))
                                   * int(c.info["n_monomers"]) * KCAL / chains, "is_all_trans": is_all_trans(c), "shear_along_chains_monomers": shear_along_chains(c),
                                   "relaxed": str(c.info.get("note", "")).split(", relaxed in ")[-1]} for c in map(read, legs)]
        out.setdefault("all_trans_controls", {"what": "an all-trans cell relaxed in the potential from another start, minus the all-trans cell the sheets "
                                                      "are read against, kcal/mol per chain (the unit of the sheets: one kink pair a chain); legs_before: "
                                                      "the cells it passed through on the way, in order"})[key] = rec
        print(f"control {key}: {rec['potential']:+.3f} kcal/mol per chain against {rec['against']}; {'all-trans' if rec['is_all_trans'] else 'NOT ALL-TRANS'}")
    Path(args.out).write_text(json.dumps(out, indent=1))
    print(f"written {args.out}")


def states_dft(args):
    """States of a row's cell labelled by DFT on cells a potential relaxed, each against the row's all-trans cell: DFT's energy, the
    potential's on the same cells, what DFT leaves on each cell, and DFT's own minimum as one step in the potential tilted to DFT's forces
    and stress predicts it (``examples/mlip_train.py relax --dft-gradient``).  kcal/mol per chain (a sheet: per kink pair) and per monomer."""
    from ase import units
    from ase.io import read

    labelled = {str(a.info["frame"]): a for f in args.frames for a in read(f, index=":")}
    stepped = {str(a.info["start_frame"]): a for f in args.stepped for a in read(f, index=":")}

    def cell(name):
        a = labelled[name]
        n, forces = int(a.info["n_monomers"]), a.arrays["REF_forces"]
        rec = {"frame": name, "atoms": len(a), "n_monomers": n, "dft_eV_per_monomer": float(a.info["REF_energy"]) / n,
               "potential_eV_per_monomer": float(a.info["potential_eV_per_monomer"]), "dft_largest_force_eV_per_A": float(np.abs(forces).max()),
               "dft_stress_xx_yy_zz_GPa": [float(v) for v in np.diag(np.asarray(a.info["REF_stress"]).reshape(3, 3)) / units.GPa],
               "dft_minutes": float(a.info["dft_seconds"]) / 60.0, "dft_kpoints": int(a.info["dft_kpoints"]), "cell_A_deg": [float(v) for v in a.cell.cellpar()],
               "shear_along_chains_monomers": shear_along_chains(a), "relaxed": str(a.info.get("note", "")).split(", relaxed in ")[-1]}
        if name in stepped:
            b = stepped[name]
            note = str(b.info.get("note", ""))
            rec["stepped"] = {"frame": str(b.info["frame"]), "dft_predicted_eV_per_monomer": float(b.info["dft_predicted_eV_per_monomer"]),
                              "came_to_rest": "NOT converged" not in note and "stopped when" not in note,
                              "steps": note.split("(")[-1].split(")")[0], "shear_along_chains_monomers": shear_along_chains(b)}
        return rec

    out = {"what": "states of a fourteen-monomer cell labelled by DFT (PBE-D3(BJ), the anchor's set-up) on cells the potential relaxed, against the row's all-trans cell: "
                   "kcal/mol per chain of the cell (for a kink sheet: per kink pair). dft: single points on the potential's cells. stepped: DFT's own minimum as one "
                   "relaxation in the potential tilted to DFT's forces and stress at the label predicts it, both cells stepped; came_to_rest false: the step was stopped "
                   "at the largest move allowed and the prediction is not settled. margin: all-trans minus the state, kcal/mol per monomer",
           "frames": [Path(f).name for f in args.frames], "stepped": [Path(f).name for f in args.stepped], "rows": {}}
    for item in args.rows:
        row, rest = item.split("=", 1)
        flat_name, *states = rest.split(",")
        flat = cell(flat_name)
        chains = 2  # the rows' cells hold two chains
        per_chain = flat["n_monomers"] / chains
        rec = {"all_trans": flat, "monomers_per_chain": per_chain, "states": {}}
        for state in states:
            label, name = state.split(":", 1)
            c = cell(name)
            entry = {"cell": c, "dft": (c["dft_eV_per_monomer"] - flat["dft_eV_per_monomer"]) * per_chain * KCAL,
                     "potential": (c["potential_eV_per_monomer"] - flat["potential_eV_per_monomer"]) * per_chain * KCAL}
            entry["potential_minus_dft"] = entry["potential"] - entry["dft"]
            entry["margin_dft_per_monomer"], entry["margin_potential_per_monomer"] = -entry["dft"] / per_chain, -entry["potential"] / per_chain
            if "stepped" in c and "stepped" in flat:
                entry["stepped"] = (c["stepped"]["dft_predicted_eV_per_monomer"] - flat["stepped"]["dft_predicted_eV_per_monomer"]) * per_chain * KCAL
                entry["stepped_settled"] = bool(c["stepped"]["came_to_rest"] and flat["stepped"]["came_to_rest"])
                entry["fall_on_stepping"] = {"state": (c["stepped"]["dft_predicted_eV_per_monomer"] - c["dft_eV_per_monomer"]) * per_chain * KCAL,
                                             "all_trans": (flat["stepped"]["dft_predicted_eV_per_monomer"] - flat["dft_eV_per_monomer"]) * per_chain * KCAL}
            rec["states"][label] = entry
            print(f"{row:26s} {label:22s} DFT {entry['dft']:+.3f}, the potential {entry['potential']:+.3f}"
                  + (f", stepped {entry['stepped']:+.3f}{'' if entry['stepped_settled'] else ' (not settled)'}" if "stepped" in entry else "")
                  + f"; DFT leaves {c['dft_largest_force_eV_per_A']:.2f} eV/A")
        helix = {k: v for k, v in rec["states"].items() if k.startswith("helix")}
        if helix:  # all-trans minus the lowest helix state, per monomer: what the tables call the margin
            low = min(helix, key=lambda k: helix[k]["dft"])
            rec.update(margin=helix[low]["margin_dft_per_monomer"], margin_potential=helix[low]["margin_potential_per_monomer"], margin_helix=low)
            if "stepped" in helix[low]:
                rec.update(margin_stepped=-helix[low]["stepped"] / per_chain, margin_stepped_settled=helix[low]["stepped_settled"])
        out["rows"][row] = rec
    Path(args.out).write_text(json.dumps(out, indent=1))
    print(f"written {args.out}")


def sheets_as_scan(args):
    """A kink-site scan (examples/kink_site_scan.py's file, which the rate model sums over) with each site's energy replaced by a
    potential's sheet at that site (``sheets-potential``), each row moved by a constant where one is given: the potential's map of the
    sites at another evaluator's level.  A site the potential has no sheet for is dropped; one whose kink it lost is marked so."""
    scan = json.loads(Path(args.scan).read_text())
    sheets = json.loads(Path(args.sheets).read_text())["rows"]
    shifts = {item.split("=")[0]: float(item.split("=")[1]) for item in args.shift}
    by_site = {(str(r["frame"]).split("_")[1], str(r["site"])): r for r in sheets.values()}
    rows, kept = {}, 0
    for key, row in scan["rows"].items():
        name, hand = key.split("@")[0], "m" if key.endswith("~") else ""
        sites = {}
        for site, rec in row["sites"].items():
            sheet = by_site.get((name, f"{site}{hand}"))
            if sheet is None:
                continue
            sites[site] = dict(rec, formation=float(sheet["potential"]) + shifts.get(name, 0.0), kink_kept=bool(sheet["kink_kept"]), converged=True,
                               potential=float(sheet["potential"]), tables_formation=float(rec["formation"]), shear_monomers=sheet["shear_along_chains_monomers"])
            kept += 1
        if sites:
            rows[key] = dict(row, sites=sites)
    out = {"what": f"the kink sites of {Path(args.scan).name}, each site's energy replaced by the sheet relaxed in the potential there ({Path(args.sheets).name}; "
                   "kcal/mol per kink pair against the row's all-trans cell in the potential)" + (", each row moved by shift_kcal_per_kink_pair" if shifts else "")
                   + "; potential: the sheet's own energy; tables_formation: the scan's",
           "energy": "kcal/mol per kink pair", "sheets": Path(args.sheets).name, "scan": Path(args.scan).name, "shift_kcal_per_kink_pair": shifts, "rows": rows}
    Path(args.out).write_text(json.dumps(out, indent=1))
    print(f"{kept} sites of {len(rows)} hands written to {args.out}")


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
    s.add_argument("--sheared", default="", help="with --scan: a result file of examples/kink_sheet_shear.py; each site's lowest sheet there (the cell sheared along its chains) "
                                                 "is written in place of the scan's upright one")
    s.add_argument("--sheet-of", default="", help="instead of --scan or --key: a key of --kinked-references (as vdf10@0); writes that kinked crystal and the row's all-trans cell "
                                                  "(the all-trans references are read from --dilute)")
    s.add_argument("--key", default="", help="the run's key there, as 'vdf4-trfe2-6@0 1x1 kinked 0,1'")
    s.add_argument("--stages", default="flexible_host", help="comma list of the stages whose kinked cell is written")
    s.add_argument("--kmesh", required=True)
    s.add_argument("--tag", default="")
    s.add_argument("--out", required=True)
    s = sub.add_parser("sheets-potential", help="kink sheets relaxed in a potential, against the rows' all-trans cells there")
    s.add_argument("--frames", default="", help="the folder of sheet_*_in_NAME.start.xyz")
    s.add_argument("--pairs", nargs="*", default=[], help="instead of --frames: KEY=ALL_TRANS.start.xyz,SHEET.start.xyz, both relaxed in the potential to the same stop")
    s.add_argument("--controls", nargs="*", default=[], help="KEY=ALL_TRANS.start.xyz,OTHER_ALL_TRANS.start.xyz[,LEG.start.xyz ...]: the second all-trans cell "
                                                           "(relaxed from another start, e.g. examples/cell_tilt.py like) against the first, kcal/mol per chain; "
                                                           "files after it are the legs it passed through, in order")
    s.add_argument("--name", default="all_s1")
    s.add_argument("--all-trans", default="docs/predictions/rows14_potential.json", help="the rows' cells in the same potential")
    s.add_argument("--out", required=True)
    s = sub.add_parser("kink-collect", help="a kink's energy from labelled frames")
    s.add_argument("--pairs", nargs="+", required=True, help="NAME=ALL_TRANS.dft.xyz,KINKED.dft.xyz: one lattice for both, or each cell its own")
    s.add_argument("--steps", nargs="*", default=[], help="NAME=KINKED_STEP1.dft.xyz[,STEP2...]: the pair's kinked cell relaxed on in the potential tilted to DFT's forces "
                                                           "(examples/mlip_train.py relax --dft-gradient) and labelled, in order")
    s.add_argument("--also", nargs="*", default=[], help="NAME=VALUE: a polyfind number with no DFT of its own, kept beside the others")
    s.add_argument("--moved", nargs="*", default=[], help="pairs whose kinked cell relaxed into another structure than the one polyfind's stored number is of: that number is left out")
    s.add_argument("--out", required=True)
    s.add_argument("--extxyz", default="", help="store the labelled frames here")
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
    s = sub.add_parser("states-dft", help="DFT labels on a row's states (all-trans, helix, kink sheets) against its all-trans cell, beside the potential")
    s.add_argument("--frames", nargs="+", required=True, help="files of labelled frames (each carries the potential's energy of the same cell)")
    s.add_argument("--stepped", nargs="*", default=[], help="files of the labelled cells stepped toward DFT's minimum (mlip_train.py relax --dft-gradient)")
    s.add_argument("--rows", nargs="+", required=True, help="ROW=ALL_TRANS_FRAME,LABEL:FRAME,LABEL:FRAME... (labels without spaces)")
    s.add_argument("--out", required=True)
    s = sub.add_parser("sheets-as-scan", help="a kink-site scan with a potential's sheets in place of its energies")
    s.add_argument("--scan", required=True, help="the scan whose sites, units and hands are kept (examples/kink_site_scan.py)")
    s.add_argument("--sheets", required=True, help="the potential's sheets at those sites (sheets-potential)")
    s.add_argument("--shift", nargs="*", default=[], help="ROW=DELTA, kcal/mol per kink pair added to every site of ROW (another evaluator minus the potential, measured at one site)")
    s.add_argument("--out", required=True)
    s = sub.add_parser("potential", help="cells relaxed in a potential and not labelled: its own energies and margins")
    s.add_argument("--frames", nargs="+", required=True, help="folders of NAME_in_POTENTIAL.start.xyz")
    s.add_argument("--rows", nargs="+", required=True)
    s.add_argument("--name", default="all_s1", help="the potential's name in the frames' names")
    s.add_argument("--out", required=True)
    args = ap.parse_args()
    {"write": write_cells, "collect": collect, "potential": potential, "kink": write_kink, "kink-collect": kink_collect, "sheets-potential": sheets_potential, "states-dft": states_dft, "sheets-as-scan": sheets_as_scan}[args.cmd](args)


if __name__ == "__main__":
    main()
