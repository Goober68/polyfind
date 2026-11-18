"""The energy between rigid chains as the cell is strained and sheared: polyfind against DFT on the same atoms.

examples/shear_registry.py measured one thing, what it costs to move neighbouring chains along
each other at the cell polyfind holds, and found polyfind three times DFT.  The variables that set
that cost also hold the lattice (docs/MLIP_SECOND_OPINION.md 3.9), so a wall that is to give both
has to be fitted to both, on one footing.  This script makes that footing: plain PVDF in three
packings (polyfind's beta cell; the all-trans packing of the reversed-unit candidate with plain
chains in it; and the TGTG' helix cell, two monomers per chain), each with its chains held rigid
and the cell changed around them:

* the lateral cell edges scaled, one at a time and together (-6% to +6%), the chains' axes at the
  same fractional places, their setting angles and heights unchanged;
* the cell sheared along its chains, at the unstrained cell and at the cells 4% narrower and wider.

Everything inside a chain is the same in every frame of a packing, so an energy difference between
two frames is energy between chains and nothing else: what a non-bonded wall has to reproduce, with
no valence term in it.

    python examples/interchain_frames.py frames --out /opt/dft/interchain
    python examples/mlip_frames.py label --frames /opt/dft/interchain --pw .../pw.x --pseudo .../pseudo --minutes 100
    python examples/interchain_frames.py collect --frames /opt/dft/interchain
"""
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np

from polyfind.fitting import PRESET_PARAMETERS


def _load(name):
    spec = importlib.util.spec_from_file_location(name, Path(__file__).with_name(f"{name}.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


SR = _load("shear_registry")
SP = _load("shear_polish")
DK, C = SR.DK, SR.C
PRED = Path("docs/predictions")
PRESET = SR.PRESET
LATERAL = (0.94, 0.96, 0.98, 0.99, 1.01, 1.02, 1.04, 1.06)
BIAXIAL = (0.96, 0.98, 1.02, 1.04)
MESH = {"beta": "4 2 8", "candidate": "4 2 8", "alpha": "4 2 4", "hh": "4 2 4", "cand6": "4 2 2"}
# The packings built here by the packing search and its polish: (row name, unit letters, conformation).
# "hh" is the chain in which every second unit is reversed, -CH2-CF2-CF2-CH2-: head to head and tail to
# tail throughout.  Plain PVDF never puts fluorine against fluorine between chains or hydrogen against
# hydrogen; this chain does, in every cell, and the reversed-unit candidate's all-trans cell is decided
# by exactly those contacts (docs/MLIP_SECOND_OPINION.md 3.10).
BUILT = {"alpha": ("vdf2", "VV", "TGTG'"), "hh": ("vdf1-rvdf1-2", "VI", "TTTT")}
# A candidate's own polished cell, its own chains held rigid: (row, conformation).  "cand6" is the
# reversed-unit candidate's all-trans cell, sheared along the line DFT's relaxation took it (one of
# ``a`` to minus three of ``b``, to (0.30, -0.90) monomers): what DFT gives for that shear with
# nothing inside a chain allowed to move is the part of its 0.35 kcal/mol that is between chains.
STORED = {"cand6": ("vdf5-rvdf1-6", "TTTT")}


def variants(pack: str) -> list:
    """``[(scale of A, scale of B, shear of a, shear of b in monomers)]`` for a packing, the unstrained cell first."""
    out = [(1.0, 1.0, 0.0, 0.0)]
    out += [(s, 1.0, 0.0, 0.0) for s in LATERAL] + [(1.0, s, 0.0, 0.0) for s in LATERAL] + [(s, s, 0.0, 0.0) for s in BIAXIAL]
    if pack in STORED:  # seventy-two atoms a cell: the unsheared cell and DFT's line, nothing else
        return [(1.0, 1.0, 0.0, 0.0), (1.0, 1.0, 0.15, -0.45), (1.0, 1.0, 0.30, -0.90), (1.0, 1.0, 0.075, -0.225)]
    if pack == "alpha":  # the helix repeats every two monomers: a quarter, a half and a whole monomer of shear
        out += [(1.0, 1.0, f, 0.0) for f in (0.25, 0.5, 1.0)] + [(1.0, 1.0, 0.0, f) for f in (0.25, 0.5, 1.0)]
    elif pack == "hh":  # two monomers to the repeat and no mirror to halve the range: every quarter along each vector, and three mixed
        out += [(1.0, 1.0, f, 0.0) for f in (0.25, 0.5, 0.75, 1.0, 1.5)] + [(1.0, 1.0, 0.0, f) for f in (0.25, 0.5, 0.75, 1.0, 1.5)]
        out += [(1.0, 1.0, 0.5, 0.5), (1.0, 1.0, 0.5, -0.5), (1.0, 1.0, 1.0, 1.0), (0.96, 0.96, 0.5, 0.0), (0.96, 0.96, 0.0, 0.5), (0.96, 0.96, 1.0, 0.0)]
    else:  # the unstrained shears are examples/shear_registry.py's; here the same at the narrower and the wider cell
        out += [(s, s, fa, fb) for s in (0.96, 1.04) for fa, fb in ((0.25, 0.0), (0.5, 0.0), (0.0, 0.5))]
    return out


def built_cell(pack: str, stored: dict | None = None) -> dict:
    """The packing search's cell of a :data:`BUILT` packing, polished in the mixed cell
    (``examples/dilute_kink.py refined_cell``); ``stored`` is that dict from an earlier run."""
    row, letters, conf = BUILT[pack]
    polymer = C.make_polymer(row, letters)
    return DK.refined_cell(polymer, row, C.conformation(conf, polymer.bonds_per_repeat), f"{row} {conf}", stored)


def alpha_cell(stored: dict | None = None) -> dict:
    return built_cell("alpha", stored)


def model_of(pack: str, alpha: dict | None = None, check: bool = True):
    """``(model, x, chain, monomers per cell)`` of a packing at no strain and no shear; ``alpha`` is
    the stored cell of a :data:`BUILT` packing, and ``check`` false builds it under parameters
    other than the ones it was polished with."""
    if pack in STORED:
        row, conf = STORED[pack]
        polymer = C.make_polymer(row, DK.sequence(row)[0])
        d = SP.stored_cell(row, conf, json.load(open(SP.POLISHED)), json.load(open(SP.REFERENCES)))
        model, x, ch = SP.sheared_model(polymer, d, check)
        return model, x, ch, 2 * ch.n_monomers
    if pack in BUILT:
        polymer = C.make_polymer(*BUILT[pack][:2])
        d = {k: (np.asarray(v, dtype=float) if isinstance(v, list) else v) for k, v in alpha.items()}
        model, x, ch = SP.sheared_model(polymer, d, check)
        return model, x, ch, 2 * ch.n_monomers
    model, x, ch = SR.plain_model(pack)
    return model, x, ch, 2


def moved(model, x0, ch, variant, monomer: float):
    sA, sB, fa, fb = variant
    x = np.array(x0, dtype=float)
    x[0] *= sA
    x[1] *= sB
    x[model.index("shear")] = [fa * monomer, fb * monomer]
    return x


def frames(args):
    from ase import Atoms
    from ase.io import write

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    table_path = out / "polyfind.json"
    table = json.loads(table_path.read_text()) if table_path.exists() else {}
    with PRESET_PARAMETERS[PRESET].applied():
        for pack in [p for p in args.packings.split(",") if p]:
            cell = built_cell(pack, table.get(pack, {}).get("cell")) if pack in BUILT else None
            model, x0, ch, n_mon = model_of(pack, cell)
            st = model.state(x0)
            monomer = float(st.lattice[3]) / (n_mon // 2)
            e0, es0 = model.physical_energy(x0), SR.electrostatics(model, x0)
            el, keep = SR.real_atoms(model)
            block = {"cell_A_deg": st.lattice.tolist(), "setting_angles_deg": st.place[:, 2].tolist(), "monomers_per_cell": n_mon,
                     "monomer_A": monomer, "polyfind_kcal_per_monomer": e0 / n_mon, "frames": []}
            if cell is not None:
                block["cell"] = DK.as_json(cell)
            for i, v in enumerate(variants(pack)):
                x = moved(model, x0, ch, v, monomer)
                e, es = model.physical_energy(x), SR.electrostatics(model, x)
                latn, _, _, Pn, _ = model.build(x)
                a = Atoms(symbols=list(el[keep]), positions=Pn[keep], cell=latn, pbc=True)
                a.wrap()
                name = f"inter_{pack}_{i:02d}"
                a.info.update(family="interchain", frame=name, n_monomers=n_mon, kmesh=MESH[pack],
                              source=(f"polyfind {PRESET}, {STORED[pack][0]} {STORED[pack][1]}, its own polished cell, rigid chains" if pack in STORED
                                      else f"polyfind {PRESET}, plain PVDF, {pack} packing, rigid chains"),
                              scale_A=v[0], scale_B=v[1], shear_monomers_a=v[2], shear_monomers_b=v[3])
                write(out / f"{name}.start.xyz", a, format="extxyz")
                block["frames"].append({"frame": name, "scale": [v[0], v[1]], "shear_monomers": [v[2], v[3]], "polyfind": (e - e0) / n_mon,
                                        "polyfind_electrostatics": (es - es0) / n_mon, "polyfind_repulsion_dispersion": ((e - es) - (e0 - es0)) / n_mon})
            table[pack] = block
            f = block["frames"]
            print(f"{pack}: cell {np.round(st.lattice, 4).tolist()}, setting angles {np.round(st.place[:, 2], 1).tolist()}, {n_mon} monomers, {len(f)} frames; polyfind "
                  + (f"at -4% and +4% biaxial {f[17]['polyfind']:+.3f}, {f[20]['polyfind']:+.3f}" if len(f) > 20 else
                     ", ".join(f"{r['shear_monomers']} {r['polyfind']:+.3f}" for r in f[1:])) + " kcal/mol per monomer")
    table_path.write_text(json.dumps(table, indent=1))
    print(f"written {out}; label with examples/mlip_frames.py label --frames {out}")


def collect(args):
    from ase.io import read, write

    src, out_path = Path(args.frames), Path(args.out)
    table = json.loads((src / "polyfind.json").read_text())
    kept = []
    done = {}
    for pack, block in table.items():
        zero = None
        rows = []
        for row in block["frames"]:
            f = src / f"{row['frame']}.dft.xyz"
            if not f.exists():
                continue
            a = read(f)
            kept.append(a)
            e = float(a.info["REF_energy"]) / a.info["n_monomers"]
            zero = e if zero is None else zero
            rows.append({**row, "dft": (e - zero) * SR.EV_TO_KCAL})
        if not rows:
            continue
        done[pack] = {**{k: v for k, v in block.items() if k != "frames"}, "frames": rows}
        print(f"{pack}: {len(rows)} of {len(block['frames'])} frames labelled")
        for r in rows:
            print(f"   A x{r['scale'][0]:.2f} B x{r['scale'][1]:.2f} shear ({r['shear_monomers'][0]:+.2f}, {r['shear_monomers'][1]:+.2f}): "
                  f"polyfind {r['polyfind']:+.4f}, DFT {r['dft']:+.4f}")
    out_path.write_text(json.dumps({
        "what": "the energy between rigid chains as the cell is strained and sheared around them, polyfind against PBE-D3(BJ) single points on the same atoms: "
                "plain PVDF in three packings (beta, candidate, alpha), the chain with every second unit reversed (hh), and the reversed-unit candidate's own "
                "polished all-trans cell with its own chains (cand6)",
        "energy": "kcal/mol per monomer against the packing's own unstrained, unsheared cell", "preset": PRESET,
        "dft": "Quantum ESPRESSO 7.5, PBE-D3(BJ), SSSP 1.3.0 precision, 90/360 Ry; 4x2x8 on the one-monomer cells, 4x2x4 on the two-monomer cells, 4x2x2 on the six-monomer cell",
        "packings": done}, indent=1))
    write(args.extxyz, kept, format="extxyz")
    print(f"written {out_path} and {args.extxyz} ({len(kept)} frames)")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="step", required=True)
    s = sub.add_parser("frames")
    s.add_argument("--out", required=True, help="working directory for the DFT, outside the repository")
    s.add_argument("--packings", default="beta,candidate,alpha,hh,cand6", help="which packings to (re)build; the others' frames and entries are kept")
    s = sub.add_parser("collect")
    s.add_argument("--frames", required=True)
    s.add_argument("--out", default=str(PRED / "interchain_belt4.json"))
    s.add_argument("--extxyz", default=str(PRED / "mlip" / "interchain_frames_pbe_d3bj.extxyz"))
    args = ap.parse_args()
    {"frames": frames, "collect": collect}[args.step](args)


if __name__ == "__main__":
    main()
