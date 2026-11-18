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
    assert prm.torsion == pytest.approx(tuple(x21[16:19])) and prm.torsion_types == (("CF2-CF2", 0.5, 0.3, 1.4), ("CH2-CH2", 0.0, 0.2, 1.3))
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


def test_the_chain1_preset_is_lattice_side_only_and_types_only_the_bonds_it_was_fitted_to():
    from polyfind import fitting as F
    from polyfind.forcefield import PRESETS

    name = "pvdf-dft-chain1"
    p = F.PRESET_PARAMETERS[name]
    assert name not in PRESETS and name in F.LATTICE_ONLY_PRESETS  # a SimpleFF would not carry its sites
    stored = json.loads((ROOT / "docs" / "predictions" / "crystal_stage_v13_typed1.json").read_text())
    assert np.array(stored["x"]) == pytest.approx(F.CRYSTAL13_CHAIN_FITTED_X) and tuple(stored["variables"]) == C.VARIABLES13_TYPED
    line = json.loads((ROOT / "docs" / "predictions" / "crystal_stage_v13_line1.json").read_text())
    assert np.array(line["x"]) == pytest.approx(F.CRYSTAL13_CHAIN_FITTED_X[:21])  # the 21 the joint fit set
    assert p.wall[0] == "b14-7" and p.bond_sites == (("H", pytest.approx(-0.3292, abs=1e-4)),) and [t[0] for t in p.torsion_types] == ["CF2-CF2", "CH2-CH2"]
    SP = CC._load("shear_polish")
    with p.applied():
        cand = SP.C.make_polymer("vdf5-rvdf1-6", "VVVVIV")
        cfe = SP.C.make_polymer("vdf7-trfe2-cfe1", SP.DK.sequence("vdf7-trfe2-cfe1")[0])
        types = F.polymer_lattice_torsion_types(cand, 12)
        assert types.count("CF2-CF2") == 1 and types.count("CH2-CH2") == 1 and types.count("CF2-CH2") == 10
        t = p.torsion_for(cand, 12)
        assert t[types.index("CF2-CF2")] == pytest.approx(p.torsion_types[0][1:]) and t[types.index("CH2-CH2")] == pytest.approx(p.torsion_types[1][1:])
        assert t[0] == pytest.approx(p.torsion)
        # a chlorinated or a CHF carbon's bonds are neither: they keep the CF2-CH2 bond's triple, as under every earlier preset
        n = 2 * len(SP.DK.sequence("vdf7-trfe2-cfe1")[0])
        other = [k for k, name in enumerate(F.polymer_lattice_torsion_types(cfe, n)) if name in ("other", "CF2-CHF")]
        assert other and all(p.torsion_for(cfe, n)[k] == pytest.approx(p.torsion) for k in other)
    # a triple given to the catch-all still reaches a CH2-CH2 bond that has none of its own
    from dataclasses import replace
    q = replace(p, torsion_types=(("other", 1.0, 2.0, 3.0),))
    with q.applied():
        assert q.torsion_for(SP.C.make_polymer("vdf5-rvdf1-6", "VVVVIV"), 12)[types.index("CH2-CH2")] == pytest.approx((1.0, 2.0, 3.0))


def test_the_acceptance_under_chain1_passes_both_halves():
    """docs/CRYSTAL_STAGE.md 5.29: the candidate's shear, its level against its helix and the order of its two helix placements."""
    a = json.loads((ROOT / "docs" / "predictions" / "candidate_acceptance_chain1.json").read_text())
    h = a["against_the_helix_dft_relaxed"]
    dft = h["dft"]
    assert a["all_trans_shear_gain"] == pytest.approx(-0.275, abs=0.003) and abs(h["helix_shear_gain"]) < 0.05
    assert abs(a["all_trans_shear_gain"] - (-0.346)) < 0.08
    for key, mine in (("margin_upright", 0.042), ("margin_sheared", -0.199), ("other_helix_minus_it", 0.242)):
        assert h[key] == pytest.approx(mine, abs=0.003) and abs(h[key] - dft[key]) < 0.07 and h[key] * dft[key] > 0
    s = json.loads((ROOT / "docs" / "predictions" / "chain_split_polyfind_chain1.json").read_text())["sets"]["pvdf-dft-chain1"]
    ref = json.loads((ROOT / "docs" / "predictions" / "chain_split_pbe_d3bj.json").read_text())["pairs"]
    level = s["vdf5-rvdf1-6"]["pairs"]["TTTT minus TG'TG"]
    want = ref["cand_vdf5-rvdf1-6_TTTT_in_pvdfb_s2_upright minus cand_vdf5-rvdf1-6_TGTGp_in_pvdfb_s2"]
    assert abs(level["inside_a_chain"] - want["inside_a_chain"]) < 0.05 and abs(level["between_chains"] - want["between_chains"]) < 0.08
    # what it does not have: the order of the placements is right as a sum, not in its parts
    order = s["vdf5-rvdf1-6"]["pairs"]["GTG'T minus TG'TG"]
    assert order["between_chains"] > 0 > ref["rvdf6_GTGpT_in_pvdfb_s2 minus cand_vdf5-rvdf1-6_TGTGp_in_pvdfb_s2"]["between_chains"]
