"""Lone chains in named conformations (examples/chain_conformers.py) and the torsion triples they are for."""

import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

from polyfind import crystalfit as C
from polyfind.fitting import PRESET_PARAMETERS

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("chain_conformers", ROOT / "examples" / "chain_conformers.py")
CC = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CC)


def test_a_pattern_gives_each_bond_type_its_angles_in_turn():
    types = ["CF2-CF2", "CF2-CH2", "other", "CF2-CH2"]
    assert CC.torsions_of(types, {}).tolist() == [180.0] * 4
    assert CC.torsions_of(types, {"CF2-CF2": (60.0,), "other": (-60.0,)}).tolist() == [60.0, 180.0, -60.0, 180.0]
    assert CC.torsions_of(types, {"CF2-CH2": (60.0, -60.0)}).tolist() == [180.0, 60.0, 180.0, -60.0]
    assert CC.torsions_of(types * 2, {"CF2-CF2": (60.0, -60.0)}).tolist() == [60.0, 180.0, 180.0, 180.0, -60.0, 180.0, 180.0, 180.0]


def test_dfts_dihedrals_go_onto_a_models_bonds_wherever_the_walk_began_and_either_way_round():
    types = ["CF2-CH2", "CF2-CH2", "CF2-CF2", "CF2-CH2", "other", "CF2-CH2"]
    # walking the same chain from two bonds further on
    bonds = ["CF2-CF2", "CF2-CH2", "CH2-CH2", "CH2-CF2", "CF2-CH2", "CH2-CF2"]
    assert CC.model_order([1, 2, 3, 4, 5, 6], bonds, types).tolist() == [5, 6, 1, 2, 3, 4]
    # and walking it the other way: the reversed unit's two bonds change places along the walk
    back = bonds[::-1]
    assert CC.model_order([6, 5, 4, 3, 2, 1], back, types).tolist() == [5, 6, 1, 2, 3, 4]
    with pytest.raises(ValueError, match="no turn"):
        CC.model_order([1, 2, 3, 4, 5, 6], ["CF2-CH2"] * 6, types)
    with pytest.raises(ValueError, match="dihedrals for a chain"):
        CC.model_order([1, 2, 3], bonds[:3], types)


def test_the_twenty_seven_variable_form_gives_the_reversed_units_bonds_their_own_triples():
    x21 = np.array(json.loads((ROOT / "docs" / "predictions" / "crystal_stage_v13_line1.json").read_text())["x"])
    x = np.concatenate([x21, [0.5, 0.3, 1.4], [0.0, 0.2, 1.3]])
    names, x0, bounds, sigma = C._config(x)
    assert names == C.VARIABLES13_TYPED and names[-6:] == ("V1_FF", "V2_FF", "V3_FF", "V1_HH", "V2_HH", "V3_HH") and len(x0) == len(bounds) == len(sigma) == 27
    prm = C.parameters(x)
    assert prm.torsion == pytest.approx(tuple(x21[16:19])) and prm.torsion_types == (("CF2-CF2", 0.5, 0.3, 1.4), ("other", 0.0, 0.2, 1.3))
    # everything else is the 21-variable set's
    base = C.parameters(x21)
    assert prm.lj == base.lj and prm.lj_pairs == base.lj_pairs and prm.bond_sites == base.bond_sites and prm.wall == base.wall
    # the ridge holds the new triples to the CF2-CH2 bond's of the same vector, weakly
    assert C.ridge_parts(x)[-6:] == pytest.approx((x[21:27] - np.tile(x21[16:19], 2)) / 2.0)
    same = np.concatenate([x21, x21[16:19], x21[16:19]])
    assert C.ridge(same) == pytest.approx(C.ridge(x21)) and all(set(v) == set() for v in [C.parameters(x21).torsion_types])
    assert all(n in C.FD_STEPS for n in names)


def test_a_chain_held_at_given_torsions_closes_and_the_reversed_polymers_fold_costs_under_belt4():
    """Polyfind's side: the chain with every second unit reversed at its nominal conformers.  Held at
    the torsions, the backbone angles and the repeat do the closing; folding it at the CF2-CF2 and
    CH2-CH2 bonds costs more than folding it at the CF2-CH2 bonds, as in DFT."""
    SP = CC._load("shear_polish")
    with PRESET_PARAMETERS["pvdf-dft-valence-belt4"].applied():
        polymer = SP.C.make_polymer("hh2", "VI")
        types = CC.chain_types(polymer, 4)
        assert sorted(types) == ["CF2-CF2", "CF2-CH2", "CF2-CH2", "other"]
        e = {}
        for conf, pattern in CC.CHAINS["hh2"][1].items():
            q = CC.polyfind_at_torsions(SP, polymer, CC.torsions_of(types, pattern), tilts=False)
            assert q["closure"] < 1e-6 and q["n_monomers"] == 2
            e[conf] = q["energy"] / 2
    assert e["special"] - e["tt"] > e["plain"] - e["tt"] > 1.0
