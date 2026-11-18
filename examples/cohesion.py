"""What holds a packing together: the cell with its chains rigid minus the same chains alone, DFT against polyfind.

A by-product of two sets of single points that already exist: the unstrained rigid cells of
examples/interchain_frames.py and the single chains of examples/chain_potential.py (the same
chain, alone in a 16 A box).  Their difference per monomer is the energy between the chains of a
packing, absolute, with everything inside a chain cancelled because the atoms of the chain are
the same in both.  No fit has seen it.

Polyfind's number is its lattice energy of the cell minus that of the same chains far apart, and
minus that of the chains on a 16 x 32 A lattice (each chain's nearest images 16 A away, as DFT's
box has them): the two bracket what the box does to a chain with a dipole.

    python examples/cohesion.py --chains /opt/dft/chainpot --cells /opt/dft/interchain --fit name=file.json [--fit ...]
"""
import argparse
import importlib.util
import json
import re
from pathlib import Path

import numpy as np

from polyfind import crystalfit as C
from polyfind import fitting as F

SPEC = importlib.util.spec_from_file_location("chain_potential", Path(__file__).with_name("chain_potential.py"))
CP = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CP)
IF = CP.IF
RY, EV = 313.7547, 23.060548
PAIRS = {"beta": ("ht", 1), "alpha": ("tgtg", 2), "hh": ("hh", 2)}  # packing: (its chain, monomers in the chain's cell)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--chains", required=True, help="working directory of examples/chain_potential.py")
    ap.add_argument("--cells", required=True, help="working directory of examples/interchain_frames.py")
    ap.add_argument("--fit", action="append", default=[], help="name=crystal-stage result file; belt4 is always there")
    ap.add_argument("--out", default="docs/predictions/cohesion_rigid_cells.json")
    args = ap.parse_args()
    from ase.io import read

    def chain(name, monomers):
        text = (Path(args.chains) / "work" / f"chain_{name}" / f"chain_{name}.out").read_text()
        return float(re.findall(r"^!\s+total energy\s+=\s+(-?[\d.]+)", text, re.M)[-1]) * RY / monomers

    def cell(pack):
        a = read(Path(args.cells) / f"inter_{pack}_00.dft.xyz")
        return a.info["REF_energy"] * EV / a.info["n_monomers"]

    dft = {p: {"cell": cell(p), "chain": chain(*c)} for p, c in PAIRS.items()}
    out = {"what": "the energy between the rigid chains of a packing: the unstrained cell minus the same chains alone, kcal/mol per monomer",
           "dft": {p: v["cell"] - v["chain"] for p, v in dft.items()},
           "dft_note": "PBE-D3(BJ) single points; the chains alone are in a 16 A box (examples/chain_potential.py)", "polyfind": {}}
    print("DFT: " + ", ".join(f"{p} {v:+.3f}" for p, v in out["dft"].items()) + f"; alpha minus beta {out['dft']['alpha'] - out['dft']['beta']:+.3f}")
    table = json.loads((Path(args.cells) / "polyfind.json").read_text())
    sets = {"belt4": F.PRESET_PARAMETERS[IF.PRESET]}
    for item in args.fit:
        name, path = item.split("=", 1)
        sets[name] = C.parameters(np.array(json.loads(Path(path).read_text())["x"], dtype=float))
    for name, prm in sets.items():
        rows = {}
        with prm.applied():
            for pack in PAIRS:
                model, x0, ch, n_mon = IF.model_of(pack, table[pack].get("cell"), check=False)
                here, far, box = {}, {}, {}
                model.value_and_grad(x0, here)
                x = np.array(x0)
                x[0] *= 8.0
                x[1] *= 8.0
                model.value_and_grad(x, far)
                x = np.array(x0)
                x[0], x[1], x[2] = 16.0, 32.0, 90.0
                model.value_and_grad(x, box)
                rows[pack] = {"against_chains_far_apart": (here["lattice"] - far["lattice"]) / n_mon, "against_a_16_A_lattice": (here["lattice"] - box["lattice"]) / n_mon}
        out["polyfind"][name] = rows
        print(f"{name:28s} " + ", ".join(f"{p} {v['against_chains_far_apart']:+.3f} | {v['against_a_16_A_lattice']:+.3f}" for p, v in rows.items())
              + f"; alpha minus beta {rows['alpha']['against_a_16_A_lattice'] - rows['beta']['against_a_16_A_lattice']:+.3f}")
    Path(args.out).write_text(json.dumps(out, indent=1))
    print(f"written {args.out}")


if __name__ == "__main__":
    main()
