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
MESH = {"beta": "4 2 8", "candidate": "4 2 8", "alpha": "4 2 4"}


def variants(pack: str) -> list:
    """``[(scale of A, scale of B, shear of a, shear of b in monomers)]`` for a packing, the unstrained cell first."""
    out = [(1.0, 1.0, 0.0, 0.0)]
    out += [(s, 1.0, 0.0, 0.0) for s in LATERAL] + [(1.0, s, 0.0, 0.0) for s in LATERAL] + [(s, s, 0.0, 0.0) for s in BIAXIAL]
    if pack == "alpha":  # the helix repeats every two monomers: a quarter, a half and a whole monomer of shear
        out += [(1.0, 1.0, f, 0.0) for f in (0.25, 0.5, 1.0)] + [(1.0, 1.0, 0.0, f) for f in (0.25, 0.5, 1.0)]
    else:  # the unstrained shears are examples/shear_registry.py's; here the same at the narrower and the wider cell
        out += [(s, s, fa, fb) for s in (0.96, 1.04) for fa, fb in ((0.25, 0.0), (0.5, 0.0), (0.0, 0.5))]
    return out


def alpha_cell(stored: dict | None = None) -> dict:
    """The packing search's TGTG' cell of plain PVDF, two monomers per chain, polished in the mixed
    cell (``examples/dilute_kink.py refined_cell``); ``stored`` is that dict from an earlier run."""
    polymer = C.make_polymer("vdf2", "VV")
    return DK.refined_cell(polymer, "vdf2", C.conformation("TGTG'", polymer.bonds_per_repeat), "vdf2 TGTG'", stored)


def model_of(pack: str, alpha: dict | None = None):
    """``(model, x, chain, monomers per cell)`` of a packing at no strain and no shear."""
    if pack == "alpha":
        polymer = C.make_polymer("vdf2", "VV")
        d = {k: (np.asarray(v, dtype=float) if isinstance(v, list) else v) for k, v in alpha.items()}
        model, x, ch = SP.sheared_model(polymer, d)
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
        alpha = alpha_cell(table.get("alpha", {}).get("cell"))
        for pack in ("beta", "candidate", "alpha"):
            model, x0, ch, n_mon = model_of(pack, alpha)
            st = model.state(x0)
            monomer = float(st.lattice[3]) / (n_mon // 2)
            e0, es0 = model.physical_energy(x0), SR.electrostatics(model, x0)
            el, keep = SR.real_atoms(model)
            block = {"cell_A_deg": st.lattice.tolist(), "setting_angles_deg": st.place[:, 2].tolist(), "monomers_per_cell": n_mon,
                     "monomer_A": monomer, "polyfind_kcal_per_monomer": e0 / n_mon, "frames": []}
            if pack == "alpha":
                block["cell"] = DK.as_json(alpha)
            for i, v in enumerate(variants(pack)):
                x = moved(model, x0, ch, v, monomer)
                e, es = model.physical_energy(x), SR.electrostatics(model, x)
                latn, _, _, Pn, _ = model.build(x)
                a = Atoms(symbols=list(el[keep]), positions=Pn[keep], cell=latn, pbc=True)
                a.wrap()
                name = f"inter_{pack}_{i:02d}"
                a.info.update(family="interchain", frame=name, n_monomers=n_mon, kmesh=MESH[pack],
                              source=f"polyfind {PRESET}, plain PVDF, {pack} packing, rigid chains",
                              scale_A=v[0], scale_B=v[1], shear_monomers_a=v[2], shear_monomers_b=v[3])
                write(out / f"{name}.start.xyz", a, format="extxyz")
                block["frames"].append({"frame": name, "scale": [v[0], v[1]], "shear_monomers": [v[2], v[3]], "polyfind": (e - e0) / n_mon,
                                        "polyfind_electrostatics": (es - es0) / n_mon, "polyfind_repulsion_dispersion": ((e - es) - (e0 - es0)) / n_mon})
            table[pack] = block
            f = block["frames"]
            print(f"{pack}: cell {np.round(st.lattice, 4).tolist()}, {n_mon} monomers, {len(f)} frames; polyfind at -4% and +4% biaxial "
                  f"{f[17]['polyfind']:+.3f}, {f[20]['polyfind']:+.3f} kcal/mol per monomer")
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
        "what": "the energy between rigid chains of plain PVDF as the cell is strained and sheared around them, polyfind against PBE-D3(BJ) single points on the same atoms",
        "energy": "kcal/mol per monomer against the packing's own unstrained, unsheared cell", "preset": PRESET,
        "dft": "Quantum ESPRESSO 7.5, PBE-D3(BJ), SSSP 1.3.0 precision, 90/360 Ry; 4x2x8 on the one-monomer cells, 4x2x4 on the helix cell",
        "packings": done}, indent=1))
    write(args.extxyz, kept, format="extxyz")
    print(f"written {out_path} and {args.extxyz} ({len(kept)} frames)")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="step", required=True)
    s = sub.add_parser("frames")
    s.add_argument("--out", required=True, help="working directory for the DFT, outside the repository")
    s = sub.add_parser("collect")
    s.add_argument("--frames", required=True)
    s.add_argument("--out", default=str(PRED / "interchain_belt4.json"))
    s.add_argument("--extxyz", default=str(PRED / "mlip" / "interchain_frames_pbe_d3bj.extxyz"))
    args = ap.parse_args()
    {"frames": frames, "collect": collect}[args.step](args)


if __name__ == "__main__":
    main()
