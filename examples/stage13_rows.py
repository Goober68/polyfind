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
* **the lone chain's line** (``docs/predictions/chain_stretch_pbe_d3bj.json``, examples/chain_stretch.py):
  plain PVDF's chain alone, held at each repeat DFT was asked at, its torsions, backbone angles
  and pendant groups relaxed: the helix from 6% shorter to 6% longer than it likes to be
  (:data:`LINE_STRETCH_MAX`), and the all-trans chain about its own repeat.  The same for the
  polymer with every second unit reversed (``chain_stretch_hh_pbe_d3bj.json``: its all-trans chain
  and its two folded ones, the CF2-CF2 and CH2-CH2 bonds gauche or the CF2-CH2 bonds), which is
  what those two bonds' own torsion triples are fitted to.  A row a point,
  the energy per monomer against the point where DFT is lowest, at ``hypot(0.05, 0.05 |DFT|)``
  kcal/mol; together they count as :data:`LINE_COUNT` rows.  Nothing between chains is in them:
  they are what holds the torsion terms and the backbone angle's stiffness to the chain itself,
  where the gap between the two crystals alone lets a torsion term pay for an error between chains
  (docs/CRYSTAL_STAGE.md 5.29).
* **the reversed unit's bonds** (``docs/predictions/chain_conformers_pbe_d3bj.json``,
  examples/chain_conformers.py): lone chains of the polymer with every second unit reversed, each
  conformer held at the torsions DFT settled it at, its energy per monomer against the all-trans
  chain's at ``hypot(0.05, 0.05 |DFT|)``; together :data:`REVERSED_COUNT` rows.  They are what a
  CF2-CF2 and a CH2-CH2 bond's torsion triples are fitted to.
* **the reversed-unit candidate's cells** (``docs/predictions/chain_split_pbe_d3bj.json``,
  examples/chain_split.py; only with ``frozen``): its two helix placements and its all-trans cell,
  each against the first placement, as a cell (the margin's tolerance, 0.05) and inside a chain
  (0.10 and 0.15: DFT's lone chain keeps the bond lengths the crystal gave it, which a chain of
  fixed bonds has no counterpart of).  Polyfind's cells are not relaxed in a row: they are the
  cells of ``frozen`` (``chain_split.py polyfind --freeze``, relaxed under an earlier parameter
  set, pendant tilts included), which is right to first order in how far the parameters have
  moved; the fit is run again from cells frozen at its result.
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
LINE_TOL = (0.05, 0.05)         # kcal/mol per monomer, and a fraction of DFT's value
LINE_COUNT = 8.0                # how many rows the lone chain's line counts as, a polymer
REVERSED_COUNT = 4.0            # and the reversed unit's lone conformers
CANDIDATE = ("vdf5-rvdf1-6", "TG'TG", {"TTTT": "cand_vdf5-rvdf1-6_TTTT_in_pvdfb_s2_upright", "GTG'T": "rvdf6_GTGpT_in_pvdfb_s2"}, "cand_vdf5-rvdf1-6_TGTGp_in_pvdfb_s2")
CANDIDATE_TOL = {"cell": 0.05, "inside": {"TTTT": 0.15, "GTG'T": 0.10}}  # kcal/mol per monomer
LINE_STRETCH_MAX = 1.065        # a branch is a row up to this many times the repeat where DFT has it lowest (the helix: 4.95 A; DFT's C-C
                                # bonds are 1 to 3% longer there, which a chain of fixed bond lengths cannot follow, and further out
                                # polyfind's helix has let go)
REFERENCE_PRESET = "pvdf-dft-valence-belt4"  # the charges the stored potential samples' patterns were built with


class Rows:
    def __init__(self, rigid=PRED / "interchain_belt4.json", potential=PRED / "chain_potential_pbe_samples_sites.npz",
                 torsions=PRED / "valence_torsions_stage13.json", trainset: str | None = None, log=None, chains=None, conformers: bool = True,
                 line=(PRED / "chain_stretch_pbe_d3bj.json", PRED / "chain_stretch_hh_pbe_d3bj.json"), reversed_chains=PRED / "chain_conformers_pbe_d3bj.json",
                 split=PRED / "chain_split_pbe_d3bj.json", frozen=None):
        self.log = log or (lambda *a, **k: None)
        t0 = time.time()
        self.CP = _load("chain_potential")
        self.IF, self.SR = self.CP.IF, self.CP.SR
        self.table = json.loads(Path(rigid).read_text())["packings"]
        # --- the potential: every pattern's potential at the sampled points, once
        data = self.CP.stored_blocks(Path(potential), list(self.CP.CHAINS) if chains is None else list(chains), 2.2)
        if any(len(d["patterns"]) != 4 for d in data.values()):
            raise ValueError("the potential samples were not collected from chains built with bond sites (chain_potential.py collect --sites)")
        self.blocks = self.CP.potential_blocks(data, extras=False)
        ref = {(a, b): float(d) for a, b, d in PRESET_PARAMETERS[REFERENCE_PRESET].charge_increments}
        self.ref_increments = (ref[("C", "H")], ref[("C", "F")])
        # --- the lone chain's line
        self.line = []  # one set of branches a polymer, each on its own scale
        for f in ([] if line is None else [line] if isinstance(line, (str, Path)) else list(line)):
            if Path(f).exists():
                self.ST = _load("chain_stretch")
                branches = {}
                for name, br in json.loads(Path(f).read_text())["branches"].items():
                    longest = LINE_STRETCH_MAX * min(br["points"], key=lambda q: q["energy"])["repeat"]
                    branches[name] = {**br, "points": [q for q in br["points"] if q["repeat"] <= longest]}
                self.line.append(branches)
        # --- the reversed unit: lone conformers, and the candidate's cells on frozen geometry
        self.reversed = None
        if reversed_chains is not None and Path(reversed_chains).exists():
            self.CC = _load("chain_conformers")
            self.reversed = json.loads(Path(reversed_chains).read_text())["conformers"]
        self.frozen = None
        if frozen is not None:
            self.CS = _load("chain_split")
            self.frozen = json.loads(Path(frozen).read_text())[CANDIDATE[0]]["cells"]
            self.split = json.loads(Path(split).read_text())["pairs"]
        # --- the conformers
        self.design = None
        path = (trainset or os.environ.get(F.REFERENCE_ENV)) if conformers else None
        if path:
            frames = list(F.load_reference(path))
            train, _ = F.split_systems(frames)
            self.design = F.ValenceDesign([f for f in frames if f.system in train])
            stored = json.loads(Path(torsions).read_text())
            self.x_valence = np.array(stored["x"], dtype=float)
            self.reference = stored["valence_fit"]
        self.log(f"stage 13's outside rows: {sum(len(b['frames']) - 1 for b in self.table.values())} rigid cells in {len(self.table)} packings, "
                 f"{len(self.blocks)} chains' potentials, {'no line of a lone chain' if not self.line else str(sum(len(b['points']) for fam in self.line for b in fam.values())) + ' points of lone chains'}, "
                 f"{'no conformers (set $' + F.REFERENCE_ENV + ')' if self.design is None else str(len(self.design.frames)) + ' conformer frames'}; "
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

    def lone(self, params) -> dict:
        out = {}
        previous = self.SR.C.SITES_FOR
        self.SR.C.SITES_FOR = params
        try:
            with params.applied():
                lines = [self.ST.polyfind_branches(self.IF.SP, self.IF, self.table, fam) for fam in self.line]
        finally:
            self.SR.C.SITES_FOR = previous
        for fam, mine in zip(self.line, lines):
            flat = [(name, q, r) for name, br in fam.items() for q, r in zip(br["points"], mine[name])]
            ref = int(np.argmin([q["energy"] for _, q, _ in flat]))
            weight = np.sqrt(LINE_COUNT / (len(flat) - 1))
            for k, (name, q, r) in enumerate(flat):
                if k != ref:
                    e, dft = r["energy"] - flat[ref][2]["energy"], q["energy"]
                    out[f"lone chain, {name} at {q['repeat']:.3f} A"] = (e, dft, weight * (e - dft) / np.hypot(LINE_TOL[0], LINE_TOL[1] * dft))
        return out

    def reversed_unit(self, params) -> dict:
        if self.reversed is None:
            return {}
        SP = self.IF.SP
        with params.applied():
            mine = self.CC.polyfind_conformers(SP, SP.C, self.reversed)
        zero = mine[self.CC.REFERENCE]["energy"]
        weight = np.sqrt(REVERSED_COUNT / (len(self.reversed) - 1))
        out = {}
        for name, r in self.reversed.items():
            if name != self.CC.REFERENCE:
                e, dft = mine[name]["energy"] - zero, r["energy"]
                out[f"lone reversed chain, {name}"] = (e, dft, weight * (e - dft) / np.hypot(LINE_TOL[0], LINE_TOL[1] * dft))
        return out

    def candidate(self, params) -> dict:
        if self.frozen is None:
            return {}
        row, base, others, base_frame = CANDIDATE
        SP = self.IF.SP
        pol, ref = json.load(open(SP.POLISHED)), json.load(open(SP.REFERENCES))
        e = {}
        with params.applied():
            polymer = SP.C.make_polymer(row, SP.DK.sequence(row)[0])
            for conf in (base, *others):
                model, _, ch = SP.sheared_model(polymer, SP.stored_cell(row, conf, pol, ref), check=False)
                f = self.frozen[conf]
                s = self.CS.polyfind_split(SP, model, np.array(f["x"], dtype=float), ch, tilts=np.array(f["tilts_deg"], dtype=float))
                e[conf] = (s["cell_tilted"] / (2 * ch.n_monomers), s["chain_tilted"] / (2 * ch.n_monomers))
        out = {}
        for conf, frame in others.items():
            dft = self.split[f"{frame} minus {base_frame}"]
            cell, inside = e[conf][0] - e[base][0], e[conf][1] - e[base][1]
            out[f"candidate {conf} minus {base}"] = (cell, dft["cell"], (cell - dft["cell"]) / CANDIDATE_TOL["cell"])
            out[f"candidate {conf} minus {base}, inside a chain"] = (inside, dft["inside_a_chain"], (inside - dft["inside_a_chain"]) / CANDIDATE_TOL["inside"][conf])
        return out

    def conformers(self, params) -> dict:
        if self.design is None:
            return {}
        x = self.x_valence.copy()
        x[:3] = params.torsion
        for name, *triple in params.torsion_types:  # a bond type the set has its own triple for (the conformer set has no CH2-CH2 bond)
            if name in F.REF_TORSION_TYPES:
                k = F.REF_TORSION_TYPES.index(name)
                x[3 * k:3 * k + 3] = triple
        self.design.lattice = params
        e = self.design.errors(x)
        return {"conformers, energy rms": (float(e["energy_rms"]), self.reference["energy_rms"], (float(e["energy_rms"]) - self.reference["energy_rms"]) / CONFORMER_TOL[0]),
                "conformers, force rms": (float(e["force_rms"]), self.reference["force_rms"], (float(e["force_rms"]) - self.reference["force_rms"]) / CONFORMER_TOL[1])}

    def __call__(self, params) -> dict:
        return {**self.rigid(params), **self.potential(params), **self.lone(params), **self.reversed_unit(params), **self.candidate(params), **self.conformers(params)}


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
