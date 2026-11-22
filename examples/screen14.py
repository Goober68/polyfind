"""A screen of fourteen-monomer compositions in the potential, each made of the film composition's own relaxed cells
(docs/MLIP_SECOND_OPINION.md 3.25).

3.24 left two things: a composition's margin has to be taken from its lowest all-trans cell, which is not the upright one
for the two compositions without a nitrile, and that low cell is not reached from upright chains.  So a candidate's cells
are not built afresh here.  Each is made of a cell of the film's composition (VVRVVRVCVVRVVR) that the potential has
already relaxed, by changing substituents where they stand (examples/registry_cells.py `morphed`): the upright all-trans
cell, the low all-trans cell, and the four placements of the helix along the chain.  A sequence is a word of unit letters,
one a monomer, read along the chain as polyfind builds it (every hand as its builder has it):

    V CH2-CF2   I CF2-CH2 (a VDF unit turned round)   R CHF-CF2   C CH2-CFCl   K CF2-CFCl   D CHCl-CF2
    E CH2-CH2   F CF2-CF2   X CHF-CHF   A CH2-CH(CN)

    python -I examples/screen14.py starts --frames FILES --source KIND=FRAME ... --row NAME=LETTERS ... [--as-is ROW:KIND=FRAME ...] --out STARTS.extxyz
    python -I examples/screen14.py table --starts STARTS.extxyz ... --relaxed DIR --out TABLE.json [--frames-out RELAXED.extxyz]
    python -I examples/screen14.py choose --table TABLE.json --window LO HI [--most N] [--skip ROW ...]
    python -I examples/screen14.py dft --table TABLE.json --labels FILES --out MARGINS.json
    python -I examples/screen14.py sheets --frames FILES --of ROW_OF_THE_SHEETS --row NAME=LETTERS ... --out SHEETS.extxyz
    python -I examples/screen14.py sheets-table --starts SHEETS.extxyz --relaxed DIR --out SHEETS.json [--frames-out RELAXED.extxyz]

`starts` writes a row's cells, row after row in the order given, one for each --source: a kind is tt_SOMETHING (an
all-trans cell) or hx_SOMETHING (a helix cell).  An all-trans chain's backbone bonds are turned a little off trans in a
pattern along the chain, the patterns are minima of their own a kcal/mol per chain apart, and a made cell keeps the
pattern of the cell it is made of (the matrix made of the film composition's upright cell rests 1.01 above the matrix's
own upright cell): so a row is given all-trans cells in several patterns, and its own where there is one (--as-is: a
cell of the row itself, taken as it is: polyfind's cell of it, or a labelled one).  `table` reads each start's relaxed
cell (examples/mlip_train.py relax: NAME_in_POTENTIAL.start.xyz) from every file of starts given and gives each row's
margin, all-trans less helix in kcal/mol per monomer, between its lowest all-trans cell and its lowest helix cell,
beside the margin from the film composition's upright pattern (tt_up) and every all-trans cell against the lowest.  `choose` prints the cells DFT is asked
about for the rows whose margin is inside a window: each one's lowest all-trans cell, its lowest helix cell and its other
all-trans cell, in two lines (the card's and the processors').  `dft` takes the margins again from the labelled cells.

`sheets` makes a row's kink sheets (a kink pair on every chain, in register) at every site of the sequence the same way:
out of the film composition's potential-relaxed sheets (docs/predictions/mlip/cells14_tight_in_all_s1.extxyz, 28 of
them: fourteen sites and the kink's two hands).  `sheets-table` gives each relaxed sheet against the row's lowest
all-trans cell in the same folder, kcal/mol per chain, with whether the kink is still one: where a row's kink sits and
how deep, which is what holds it (3.23: the film composition's sits 4.2 deep on its chlorine and its end moves over 7.70;
the leader's 0.7 deep on a run of VDF units and moves over 5.13).

The potential's own numbers: it is within 0.1 to 0.4 kcal/mol per chain of DFT on the low cells of the two compositions
without a nitrile and 1.3 to 1.8 off, with a sign wrong, on those with one (3.24).  Rows with an A are controls here, not
candidates.
"""

import argparse
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

KCAL = 23.0605
CHAINS = 2
# a unit's two backbone carbons in the order the chain is built: each one's substituents (in the alphabet) and their hand
UNITS = {"V": (("HH", 0), ("FF", 0)), "I": (("FF", 0), ("HH", 0)), "R": (("FH", 1), ("FF", 0)), "C": (("HH", 0), ("ClF", 1)), "K": (("FF", 0), ("ClF", 1)),
         "D": (("ClH", 1), ("FF", 0)), "E": (("HH", 0), ("HH", 0)), "F": (("FF", 0), ("FF", 0)), "X": (("FH", 1), ("FH", 1)), "A": (("HH", 0), ("CH", 1))}
KINDS = ("tt_low", "tt_up", "hx_TGTGp", "hx_GpTGT", "hx_GTGpT", "hx_TGpTG")


def _tool(name):
    if name not in sys.modules:
        spec = importlib.util.spec_from_file_location(name, Path(__file__).with_name(f"{name}.py"))
        module = importlib.util.module_from_spec(spec)
        sys.modules[name] = module
        spec.loader.exec_module(module)
    return sys.modules[name]


def cycle_of(letters):
    """A chain's cycle of substituents for a sequence of unit letters (as registry_cells.substituents reads one off a cell)."""
    return [pair for letter in letters for pair in UNITS[letter]]


def letters_of(cycle):
    """The sequence a cycle of substituents spells, read from its first carbon (None when it is not a whole number of known units)."""
    back = {pair: letter for letter, pair in UNITS.items()}
    pairs = [tuple(cycle[k:k + 2]) for k in range(0, len(cycle), 2)]
    return "".join(back[p] for p in pairs) if len(cycle) % 2 == 0 and all(p in back for p in pairs) else None


def made(source, letters):
    """The cell `source` with its chains made the sequence `letters`, and how many backbone carbons that changed."""
    reg = _tool("registry_cells")
    want = cycle_of(letters)
    if len(want) != len(reg.substituents(source)[0]):
        raise SystemExit(f"{letters}: {len(letters)} units, the cell's chains have {len(reg.substituents(source)[0]) // 2}")
    cell, changed = reg.morphed(source, want)
    for chain in reg.substituents(cell):
        if not reg.same_molecule(chain, want):
            raise SystemExit(f"{letters}: the cell made is not that chain")
    return cell, changed


def starts(args):
    from ase.io import read, write

    reg, path = _tool("registry_cells"), _tool("mlip_path")
    frames = {str(a.info["frame"]): a for f in args.frames for a in read(f, index=":")}
    sources = dict(item.split("=", 1) for item in args.source)
    if any(not k.startswith(("tt_", "hx_")) for k in sources):
        raise SystemExit("a kind is tt_SOMETHING or hx_SOMETHING")
    own = {}
    for item in args.as_is:                      # ROW:KIND=FRAME: a cell of the row itself, taken as it is
        key, frame = item.split("=", 1)
        own.setdefault(key.split(":")[0], []).append((key.split(":")[1], frame))
    out = []
    for item in args.row:
        name, letters = item.split("=", 1)
        first = len(out)
        for kind, frame in own.get(name, []):
            source = frames[frame]
            if not all(reg.same_molecule(chain, cycle_of(letters)) for chain in reg.substituents(source)):
                raise SystemExit(f"{frame} is not a cell of {letters}")
            cell = source.copy()
            for k in list(cell.arrays):
                if k not in ("numbers", "positions"):
                    del cell.arrays[k]
            cell.wrap()
            cell.info = {"n_monomers": int(source.info["n_monomers"]), "kmesh": str(source.info["kmesh"]).strip("[]"), "frame": f"scr_{name}_{kind}", "screen_row": name,
                         "screen_letters": letters, "screen_kind": kind, "screen_from": frame, "screen_changed": 0, "family": "relaxed",
                         "conformation": "all-trans" if kind.startswith("tt") else "helix", "source": f"examples/screen14.py: {frame} as it is, a start"}
            out.append(cell)
        for kind in sources:
            source = frames[sources[kind]]
            cell, changed = made(source, letters)
            before, after = reg.gauche(source, path), reg.gauche(cell, path)
            if before != after:
                raise SystemExit(f"{name} {kind}: {after} gauche torsions, the cell it was made of has {before}")
            cell.info = {"n_monomers": int(source.info["n_monomers"]), "kmesh": str(source.info["kmesh"]).strip("[]"), "frame": f"scr_{name}_{kind}", "screen_row": name,
                         "screen_letters": letters, "screen_kind": kind, "screen_from": sources[kind], "screen_changed": changed, "family": "relaxed",
                         "conformation": "all-trans" if kind.startswith("tt") else "helix",
                         "source": f"examples/screen14.py: {sources[kind]} with its chains made {letters} where they stand ({changed} carbons a cell), a start, not relaxed"}
            out.append(cell)
        print(f"{name:28s} {letters}  {out[-1].get_chemical_formula()}  " + " ".join(f"{a.info['screen_kind']}:{a.info['screen_changed']}" for a in out[first:]) + "  (kind: carbons changed)")
    write(args.out, out, format="extxyz")
    print(f"{len(out)} starts written to {args.out}")


def table(args):
    from ase.io import read, write

    reg, path = _tool("registry_cells"), _tool("mlip_path")
    rows, frames = {}, []
    every = [a for f in args.starts for a in read(f, index=":")]
    per_chain = int(every[0].info["n_monomers"]) / CHAINS * KCAL
    for start in every:
        name, kind = str(start.info["screen_row"]), str(start.info["screen_kind"])
        row = rows.setdefault(name, {"letters": str(start.info["screen_letters"]), "cells": {}})
        found = sorted(Path(args.relaxed).glob(f"{start.info['frame']}_in_*.start.xyz"))
        if not found:
            row["cells"][kind] = {"relaxed": False}
            continue
        a = read(found[0])
        note = str(a.info.get("note", ""))
        g, g0 = reg.gauche(a, path), reg.gauche(start, path)
        row["cells"][kind] = {"relaxed": True, "frame": str(a.info["frame"]), "converged": "NOT converged" not in note,
                              "steps": int(note.split("(")[-1].split(" steps")[0]) if " steps" in note else None,
                              "eV_per_monomer": float(a.info["potential_eV_per_monomer"]), "gauche_torsions": g, "gauche_torsions_start": g0,
                              "kept": g == 0 if kind.startswith("tt") else abs(g - g0) <= 2, "carbons_changed": int(start.info["screen_changed"]),
                              "shear_along_chains_monomers": path._shear(a), "cell": [float(v) for v in a.cell.cellpar()]}
        frames.append(a)
    for name, row in rows.items():
        ok = {k: c for k, c in row["cells"].items() if c.get("relaxed") and c["kept"]}
        flat, helix = {k: c for k, c in ok.items() if k.startswith("tt")}, {k: c for k, c in ok.items() if k.startswith("hx")}
        row["at_rest"] = sum(bool(c.get("converged")) for c in row["cells"].values())
        row["complete"] = row["at_rest"] == len(row["cells"]) and len(ok) == len(row["cells"])
        if not flat or not helix:
            print(f"{name:28s} {row['letters']}  no margin: {len(flat)} all-trans and {len(helix)} helix cells that kept their conformation")
            continue
        energy = lambda c: c["eV_per_monomer"]
        row["lowest_tt"], row["lowest_hx"] = min(flat, key=lambda k: energy(flat[k])), min(helix, key=lambda k: energy(helix[k]))
        floor = energy(helix[row["lowest_hx"]])
        row["margin"] = (energy(flat[row["lowest_tt"]]) - floor) * KCAL
        if "tt_up" in flat:
            row["margin_upright"] = (energy(flat["tt_up"]) - floor) * KCAL
        if "tt_low" in flat and "tt_up" in flat:
            row["low_less_upright_per_chain"] = (energy(flat["tt_low"]) - energy(flat["tt_up"])) * per_chain
        row["all_trans_per_chain"] = {k: (energy(c) - energy(flat[row["lowest_tt"]])) * per_chain for k, c in flat.items()}
        row["helix_placements_per_chain"] = {k: (energy(c) - floor) * per_chain for k, c in helix.items()}
        moving = [k for k, c in row["cells"].items() if c.get("relaxed") and not c["converged"]]
        print(f"{name:28s} {row['letters']}  margin {row['margin']:+.3f}" + (f" (upright {row['margin_upright']:+.3f})" if "margin_upright" in row else "")
              + f"; all-trans " + " ".join(f"{k[3:]} {v:+.2f}" for k, v in sorted(row["all_trans_per_chain"].items(), key=lambda kv: kv[1]))
              + f"; helix {row['lowest_hx'][3:]}; {row['at_rest']} of {len(row['cells'])} at rest" + (f", moving: {' '.join(moving)}" if moving else "")
              + ("" if len(ok) == sum(bool(c.get('relaxed')) for c in row['cells'].values()) else "; a cell left its conformation"))
    out = {"what": "fourteen-monomer compositions, each made of the film composition's potential-relaxed cells by changing substituents where they stand "
                   "(examples/screen14.py) and relaxed in the potential: margin is all-trans less helix in kcal/mol per monomer between the lowest all-trans cell "
                   "and the lowest helix cell, margin_upright from the upright all-trans cell; the potential's own energies with D3, no DFT",
           "starts": [str(f) for f in args.starts], "rows": rows}
    Path(args.out).write_text(json.dumps(out, indent=1) + "\n")
    if args.frames_out:
        write(args.frames_out, frames, format="extxyz")
    print(f"{len(frames)} cells of {len(rows)} rows; written {args.out}")


def chosen(rows, low, high, most=None, skip=()):
    """The rows whose margin is inside [low, high], complete ones first and then by how near the middle of the window, and for
    each the cells DFT is asked about: lowest all-trans, lowest helix, the next all-trans cell up."""
    inside = [(not r.get("complete", False), abs(r["margin"] - 0.5 * (low + high)), name) for name, r in rows.items() if "margin" in r and low <= r["margin"] <= high and name not in skip]
    picked = []
    for _, _, name in sorted(inside)[:most]:
        r = rows[name]
        flat = sorted((v, k) for k, v in r["all_trans_per_chain"].items() if k != r["lowest_tt"])     # the next all-trans cell up: a check on which is lowest
        picked.append((name, [r["cells"][k]["frame"] for k in [r["lowest_tt"], r["lowest_hx"]] + [k for _, k in flat[:1]]]))
    return picked


def choose(args):
    rows = json.loads(Path(args.table).read_text())["rows"]
    picked = chosen(rows, args.window[0], args.window[1], args.most, set(args.skip))
    cells = [frame for _, frames in picked for frame in frames]
    print(" ".join(cells[0::2]))
    print(" ".join(cells[1::2]))
    for name, frames in picked:
        print(f"# {name}: margin {rows[name]['margin']:+.3f} in the potential; {' '.join(frames)}", file=sys.stderr)


def dft(args):
    from ase.io import read

    stored = json.loads(Path(args.table).read_text())["rows"]
    labelled = {str(a.info["frame"]): a for f in args.labels for a in read(f, index=":")}
    out = {"what": "the screen's rows that have DFT labels (PBE-D3(BJ) single points on the cells the potential relaxed): margin, all-trans less helix in kcal/mol per "
                   "monomer, between the lowest labelled all-trans cell and the lowest labelled helix cell, in DFT and in the potential on the same two cells; "
                   "cells: eV per monomer", "rows": {}}
    for name, row in stored.items():
        have = {k: labelled[c["frame"]] for k, c in row["cells"].items() if c.get("relaxed") and c["frame"] in labelled}
        flat, helix = {k: a for k, a in have.items() if k.startswith("tt")}, {k: a for k, a in have.items() if k.startswith("hx")}
        if not flat or not helix:
            continue
        e = lambda a: float(a.info["REF_energy"]) / int(a.info["n_monomers"])
        low_tt, low_hx = min(flat, key=lambda k: e(flat[k])), min(helix, key=lambda k: e(helix[k]))
        per_chain = int(flat[low_tt].info["n_monomers"]) / CHAINS * KCAL
        rec = {"letters": row["letters"], "lowest_tt": low_tt, "lowest_hx": low_hx, "margin": (e(flat[low_tt]) - e(helix[low_hx])) * KCAL,
               "margin_potential_same_cells": (row["cells"][low_tt]["eV_per_monomer"] - row["cells"][low_hx]["eV_per_monomer"]) * KCAL, "margin_potential": row.get("margin"),
               "cells": {k: {"frame": str(a.info["frame"]), "dft_eV_per_monomer": e(a), "potential_eV_per_monomer": row["cells"][k]["eV_per_monomer"],
                             "dft_largest_force_eV_per_A": float(np.abs(a.arrays["REF_forces"]).max())} for k, a in have.items()}}
        if "tt_low" in flat and "tt_up" in flat:
            rec["low_less_upright_per_chain"] = (e(flat["tt_low"]) - e(flat["tt_up"])) * per_chain
            rec["low_less_upright_per_chain_potential"] = (row["cells"]["tt_low"]["eV_per_monomer"] - row["cells"]["tt_up"]["eV_per_monomer"]) * per_chain
        rec["all_trans_per_chain"] = {k: (e(a) - e(flat[low_tt])) * per_chain for k, a in flat.items()}
        rec["all_trans_per_chain_potential"] = {k: (row["cells"][k]["eV_per_monomer"] - row["cells"][low_tt]["eV_per_monomer"]) * per_chain for k in flat}
        out["rows"][name] = rec
        print(f"{name:28s} {row['letters']}  DFT margin {rec['margin']:+.3f} (the potential on the same cells {rec['margin_potential_same_cells']:+.3f}); all-trans "
              + " ".join(f"{k[3:]} {v:+.2f} ({rec['all_trans_per_chain_potential'][k]:+.2f})" for k, v in sorted(rec["all_trans_per_chain"].items(), key=lambda kv: kv[1])))
    Path(args.out).write_text(json.dumps(out, indent=1) + "\n")
    print(f"{len(out['rows'])} rows with a DFT margin; written {args.out}")


def sheets(args):
    import re

    from ase.io import read, write

    reg, path = _tool("registry_cells"), _tool("mlip_path")
    pattern = re.compile(rf"sheet_{re.escape(args.of)}_belt4_at(\d+m?)_sheared_in_all_s1")
    sources = {}
    for f in args.frames:
        for a in read(f, index=":"):
            m = pattern.fullmatch(str(a.info["frame"]))
            if m:
                sources[m.group(1)] = a
    if not sources:
        raise SystemExit(f"no sheet of {args.of} in {args.frames}")
    order = sorted(sources, key=lambda site: (int(site.rstrip("m")), site.endswith("m")))
    if args.sites:
        missing = [site for site in args.sites if site not in sources]
        if missing:
            raise SystemExit(f"no sheet of {args.of} at {' '.join(missing)}")
        order = list(args.sites)
    out = []
    for item in args.row:
        name, letters = item.split("=", 1)
        for site in order:
            source = sources[site]
            cell, changed = made(source, letters)
            if reg.gauche(cell, path) != reg.gauche(source, path):
                raise SystemExit(f"{name} at {site}: the kink did not come through")
            cell.info = {"n_monomers": int(source.info["n_monomers"]), "kmesh": str(source.info["kmesh"]).strip("[]"), "frame": f"scr_{name}_sheet{args.label}_at{site}", "screen_row": name,
                         "screen_letters": letters, "screen_kind": f"sheet_at{site}", "screen_site": site, "screen_from": str(source.info["frame"]), "screen_changed": changed,
                         "family": "relaxed", "conformation": "kink sheet",
                         "source": f"examples/screen14.py: {source.info['frame']} with its chains made {letters} where they stand ({changed} carbons a cell), a start, not relaxed"}
            out.append(cell)
        print(f"{name:28s} {letters}  {len(order)} sheets, sites {' '.join(order)}")
    write(args.out, out, format="extxyz")
    print(f"{len(out)} starts written to {args.out}")


def kink_of(atoms, path):
    """Each chain's gauche bonds as (place along the chain, torsion in degrees), and whether every chain holds one kink:
    two gauche bonds of opposite sense about a trans one."""
    quads = path.backbone_torsions(atoms)
    torsions = path.torsions_deg(atoms, quads)
    half = len(quads) // CHAINS
    out, kept = [], True
    for chain in range(CHAINS):
        t = torsions[chain * half:(chain + 1) * half]
        g = [int(k) for k in np.flatnonzero(np.abs(t) < 120.0)]
        out.append([(k, float(t[k])) for k in g])
        kept = kept and len(g) == 2 and (g[1] - g[0]) % half in (2, half - 2) and t[g[0]] * t[g[1]] < 0
    return out, bool(kept)


def sheets_table(args):
    from ase.io import read, write

    path = _tool("mlip_path")
    rows, frames = {}, []
    for start in read(args.starts, index=":"):
        name, site = str(start.info["screen_row"]), str(start.info["screen_site"])
        row = rows.setdefault(name, {"letters": str(start.info["screen_letters"]), "sites": {}})
        if "all_trans" not in row and not row.get("no_all_trans"):       # the row's lowest all-trans cell at rest in the same folder
            flat = {}
            for found in sorted(Path(args.relaxed).glob(f"scr_{name}_tt_*_in_*.start.xyz")):      # every pattern the row has been given (tt_low, tt_up, tt_reg, ...)
                a = read(found)
                if "NOT converged" not in str(a.info.get("note", "")) and not any(kink_of(a, path)[0]):
                    flat[str(a.info["screen_kind"])] = float(a.info["potential_eV_per_monomer"])
            if not flat:
                row["no_all_trans"] = True
                print(f"{name}: no all-trans cell at rest in {args.relaxed}: its sheets are left out")
            else:
                row["all_trans"] = min(flat, key=flat.get)
                row["all_trans_eV_per_monomer"] = flat[row["all_trans"]]
        if row.get("no_all_trans"):
            continue
        found = sorted(Path(args.relaxed).glob(f"{start.info['frame']}_in_*.start.xyz"))
        if not found:
            row["sites"][site] = {"relaxed": False}
            continue
        a = read(found[0])
        note = str(a.info.get("note", ""))
        gauche, kept = kink_of(a, path)
        per_chain = int(a.info["n_monomers"]) / CHAINS * KCAL
        row["sites"][site] = {"relaxed": True, "frame": str(a.info["frame"]), "converged": "NOT converged" not in note,
                              "energy": (float(a.info["potential_eV_per_monomer"]) - row["all_trans_eV_per_monomer"]) * per_chain, "kink_kept": kept,
                              "gauche_bonds": gauche, "units": str(start.info["screen_letters"])[(int(site.rstrip("m")) // 2) % len(str(start.info["screen_letters"]))]
                              + str(start.info["screen_letters"])[(int(site.rstrip("m")) // 2 + 1) % len(str(start.info["screen_letters"]))],
                              "shear_along_chains_monomers": path._shear(a)}
        frames.append(a)
    for name, row in rows.items():
        held = {k: v for k, v in row["sites"].items() if v.get("relaxed") and v["kink_kept"]}
        if held:
            row["cheapest"] = min(held, key=lambda k: held[k]["energy"])
            row["below_all_trans"] = sum(v["energy"] < 0 for v in held.values())
            c = held[row["cheapest"]]
            print(f"{name:28s} {row['letters']}  cheapest sheet {c['energy']:+.2f} per chain at {row['cheapest']} ({c['units']}), {row['below_all_trans']} of {len(held)} kept sheets "
                  f"below all-trans ({row['all_trans']}); {sum(bool(v.get('converged')) for v in row['sites'].values())} of {len(row['sites'])} at rest")
            for site, v in sorted(held.items(), key=lambda kv: kv[1]["energy"])[:6]:
                print(f"    at {site:4s} {v['units']}  {v['energy']:+6.2f}{'' if v['converged'] else '  NOT at rest'}")
    out = {"what": "a row's kink sheets (a kink pair on every chain, in register, at one site of the sequence), each made of the film composition's potential-relaxed sheet "
                   "at that site by changing substituents where they stand (examples/screen14.py sheets) and relaxed in the potential with the cell free: energy is "
                   "kcal/mol per chain against the row's lowest all-trans cell in the potential; kink_kept: every chain still has two gauche bonds of opposite sense about "
                   "a trans one; units: the two units the site lies on, by the film composition's numbering of sites. No DFT", "starts": str(args.starts), "rows": rows}
    Path(args.out).write_text(json.dumps(out, indent=1) + "\n")
    if args.frames_out:
        write(args.frames_out, frames, format="extxyz")
    print(f"{len(frames)} sheets of {len(rows)} rows; written {args.out}")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="step", required=True)
    s = sub.add_parser("starts")
    s.add_argument("--frames", nargs="+", required=True, help="files that hold the film composition's relaxed cells")
    s.add_argument("--source", nargs="*", default=[], metavar="KIND=FRAME", help="the relaxed cell each kind is made of (tt_low, tt_up, hx_TGTGp, ...)")
    s.add_argument("--as-is", nargs="*", default=[], metavar="ROW:KIND=FRAME", help="a cell of the row itself, taken as it is")
    s.add_argument("--row", nargs="+", required=True, metavar="NAME=LETTERS")
    s.add_argument("--out", required=True)
    s = sub.add_parser("table")
    s.add_argument("--starts", nargs="+", required=True)
    s.add_argument("--relaxed", required=True, help="the folder of NAME_in_POTENTIAL.start.xyz")
    s.add_argument("--out", required=True)
    s.add_argument("--frames-out", default="")
    s = sub.add_parser("choose")
    s.add_argument("--table", required=True)
    s.add_argument("--window", nargs=2, type=float, required=True, metavar=("LOW", "HIGH"), help="kcal/mol per monomer: the margins DFT is asked about")
    s.add_argument("--most", type=int, default=None, help="at most this many rows")
    s.add_argument("--skip", nargs="*", default=[], help="rows not to choose (they have DFT already)")
    s = sub.add_parser("dft")
    s.add_argument("--table", required=True)
    s.add_argument("--labels", nargs="+", required=True)
    s.add_argument("--out", required=True)
    s = sub.add_parser("sheets")
    s.add_argument("--frames", nargs="+", required=True, help="files that hold the relaxed sheets the rows' are made of")
    s.add_argument("--of", default="vdf9-trfe4-cfe1-14", help="the row whose sheets those are (the film's composition)")
    s.add_argument("--row", nargs="+", required=True, metavar="NAME=LETTERS")
    s.add_argument("--sites", nargs="*", default=[], help="only these sites of the row's sheets, in this order (18m 18 16m ...); all of them if not given")
    s.add_argument("--label", default="", help="put after `sheet` in the frames' names (scr_ROW_sheetLABEL_atSITE), to keep a row's sheets made of two rows' apart")
    s.add_argument("--out", required=True)
    s = sub.add_parser("sheets-table")
    s.add_argument("--starts", required=True)
    s.add_argument("--relaxed", required=True, help="the folder of NAME_in_POTENTIAL.start.xyz, the rows' all-trans cells among them")
    s.add_argument("--out", required=True)
    s.add_argument("--frames-out", default="")
    args = p.parse_args()
    {"starts": starts, "table": table, "choose": choose, "dft": dft, "sheets": sheets, "sheets-table": sheets_table}[args.step](args)


if __name__ == "__main__":
    main()
