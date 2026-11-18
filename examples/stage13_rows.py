"""The rows stage 13 of the crystal fit takes from outside ``polyfind.crystalfit``.

``Rows()(params)`` returns ``{name: (polyfind, reference, z)}`` for a parameter set, for
``crystalfit.evaluate(extra_rows=...)`` (``examples/crystal_stage.py --stage 13 --rows``).  Three
kinds, every one read from files in the repository:

* **the rigid cells** (``docs/predictions/interchain_belt4.json``, examples/interchain_frames.py):
  the energy between rigid chains of 122 strained and sheared cells against DFT on the same
  atoms, evaluated in the kernels under the parameter set, chains built with its bond sites.  A
  row is ``(polyfind - DFT) / hypot(0.03, 0.15 |DFT|)`` kcal/mol per monomer, and each packing's
  rows together count as :data:`COUNT` rows however many cells it has, so that a packing of four
  72-atom cells is not outvoted by one of thirty-seven small ones.
* **the potential around single chains** (``docs/predictions/chain_potential_pbe_samples_sites.npz``,
  examples/chain_potential.py): the rms error of polyfind's charges and the dipoles they induce
  against DFT's electrostatic potential, one row a chain, at a tolerance of 0.15 kcal/mol per e.
  The potential is linear in the two bond increments, the site's charge and the backbone pull, so
  it costs nothing per evaluation; the induced dipoles are those of the polarizabilities and the
  damping the chains were built with, which the stage does not vary.
* **the conformers** (the reference set named by ``$POLYFIND_TRAINSET``): the energy and force
  rms of the training systems' DFT frames under the parameter set's non-bonded terms and its
  torsion triple (``ValenceDesign.lattice``), against what the valence fit itself has, at
  tolerances of 0.15 kcal/mol and 0.15 kcal/mol/A.  The triple is the CF2-CH2 type's; the other
  three types' are held at :data:`OTHER_TORSIONS` (refitted under the potential's charges).
"""
import importlib.util
import json
import os
import time
from pathlib import Path

import numpy as np

from polyfind import fitting as F
from polyfind.fitting import PRESET_PARAMETERS

ROOT = Path(__file__).resolve().parents[1]
PRED = ROOT / "docs" / "predictions"


def _load(name):
    spec = importlib.util.spec_from_file_location(name, Path(__file__).with_name(f"{name}.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# how many rows each packing's rigid cells count as
COUNT = {"beta": 6.0, "candidate": 6.0, "alpha": 6.0, "hh": 12.0, "cand6": 9.0}
RIGID_TOL = (0.03, 0.15)        # kcal/mol per monomer, and a fraction of DFT's value
POTENTIAL_TOL = 0.15            # kcal/mol per e
CONFORMER_TOL = (0.15, 0.15)    # kcal/mol, kcal/mol/A
REFERENCE_PRESET = "pvdf-dft-valence-belt4"  # the charges the stored potential samples' patterns were built with


class Rows:
    def __init__(self, rigid=PRED / "interchain_belt4.json", potential=PRED / "chain_potential_pbe_samples_sites.npz",
                 torsions=PRED / "valence_torsions_stage13.json", trainset: str | None = None, log=None):
        self.log = log or (lambda *a, **k: None)
        t0 = time.time()
        self.CP = _load("chain_potential")
        self.IF, self.SR = self.CP.IF, self.CP.SR
        self.table = json.loads(Path(rigid).read_text())["packings"]
        # --- the potential: every pattern's potential at the sampled points, once
        data = self.CP.stored_blocks(Path(potential), list(self.CP.CHAINS), 2.2)
        if any(len(d["patterns"]) != 4 for d in data.values()):
            raise ValueError("the potential samples were not collected from chains built with bond sites (chain_potential.py collect --sites)")
        self.blocks = self.CP.potential_blocks(data)
        ref = {(a, b): float(d) for a, b, d in PRESET_PARAMETERS[REFERENCE_PRESET].charge_increments}
        self.ref_increments = (ref[("C", "H")], ref[("C", "F")])
        # --- the conformers
        self.design = None
        path = trainset or os.environ.get(F.REFERENCE_ENV)
        if path:
            frames = list(F.load_reference(path))
            train, _ = F.split_systems(frames)
            self.design = F.ValenceDesign([f for f in frames if f.system in train])
            stored = json.loads(Path(torsions).read_text())
            self.x_valence = np.array(stored["x"], dtype=float)
            self.reference = stored["valence_fit"]
        self.log(f"stage 13's outside rows: {sum(len(b['frames']) - 1 for b in self.table.values())} rigid cells in {len(self.table)} packings, "
                 f"{len(self.blocks)} chains' potentials, {'no conformers (set $' + F.REFERENCE_ENV + ')' if self.design is None else str(len(self.design.frames)) + ' conformer frames'}; "
                 f"{time.time() - t0:.0f} s")

    # ------------------------------------------------------------------ the three kinds
    def rigid(self, params) -> dict:
        out = {}
        previous = self.SR.C.SITES_FOR
        self.SR.C.SITES_FOR = params
        try:
            with params.applied():
                for pack, block in self.table.items():
                    model, x0, ch, n_mon = self.IF.model_of(pack, block.get("cell"), check=False)
                    monomer = float(model.state(x0).lattice[3]) / (n_mon // 2)
                    e0 = model.physical_energy(x0)
                    rows = block["frames"][1:]
                    weight = np.sqrt(COUNT[pack] / len(rows))
                    for r in rows:
                        v = (r["scale"][0], r["scale"][1], r["shear_monomers"][0], r["shear_monomers"][1])
                        e = (model.physical_energy(self.IF.moved(model, x0, ch, v, monomer)) - e0) / n_mon
                        out[f"rigid {r['frame'][6:]}"] = (e, r["dft"], weight * (e - r["dft"]) / np.hypot(RIGID_TOL[0], RIGID_TOL[1] * r["dft"]))
        finally:
            self.SR.C.SITES_FOR = previous
        return out

    def potential(self, params) -> dict:
        inc = {(a, b): float(d) for a, b, d in params.charge_increments}
        w = {"s_CH": inc[("C", "H")] / self.ref_increments[0], "s_CF": inc[("C", "F")] / self.ref_increments[1],
             "d_CC": -2.0 * dict(params.backbone_pull).get("F", 0.0), "M_H": -dict(params.bond_sites).get("H", 0.0)}
        out = {}
        for name, b in self.blocks.items():
            rms = self.CP.rms_error(b, w)
            out[f"potential around {name}"] = (rms, 0.0, rms / POTENTIAL_TOL)
        return out

    def conformers(self, params) -> dict:
        if self.design is None:
            return {}
        x = self.x_valence.copy()
        x[:3] = params.torsion
        self.design.lattice = params
        e = self.design.errors(x)
        return {"conformers, energy rms": (float(e["energy_rms"]), self.reference["energy_rms"], (float(e["energy_rms"]) - self.reference["energy_rms"]) / CONFORMER_TOL[0]),
                "conformers, force rms": (float(e["force_rms"]), self.reference["force_rms"], (float(e["force_rms"]) - self.reference["force_rms"]) / CONFORMER_TOL[1])}

    def __call__(self, params) -> dict:
        return {**self.rigid(params), **self.potential(params), **self.conformers(params)}


def summary(rows: dict) -> str:
    """The rigid rows by packing (rms error in kcal/mol per monomer, the sum of z squared), then the rest."""
    lines = []
    for pack in COUNT:
        mine = {k: v for k, v in rows.items() if k.startswith(f"rigid {pack}_")}
        if mine:
            err = np.array([v[0] - v[1] for v in mine.values()])
            lines.append(f"rigid {pack:10s} {len(mine):3d} cells: rms error {np.sqrt((err ** 2).mean()):.4f}, sum of z^2 {sum(v[2] ** 2 for v in mine.values()):.2f}")
    lines += [f"{k:28s} {v[0]:9.4f}  reference {v[1]:9.4f}  z {v[2]:+.2f}" for k, v in rows.items() if not k.startswith("rigid ")]
    return "\n".join(lines)


if __name__ == "__main__":
    import sys

    from polyfind import crystalfit as C

    fit = json.loads(Path(sys.argv[1]).read_text())
    C.BELT_GEOMETRY = "pbe"
    rows = Rows(log=print)
    t0 = time.time()
    out = rows(C.parameters(np.array(fit["x"])))
    print(summary(out))
    print(f"sum of z^2 {sum(v[2] ** 2 for v in out.values()):.2f}; {time.time() - t0:.1f} s")
