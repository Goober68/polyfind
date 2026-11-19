"""Rows' cells against DFT: the stored results (examples/row_cells_frames.py collect) and the tool's bookkeeping."""

import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
PRED = ROOT / "docs" / "predictions"
NEEDED_KCAL = 0.1


def _stored(row):
    return json.loads((PRED / f"row_cells_dft_{row}.json").read_text())


def test_the_trfe_cells_margin_is_met_by_both_presets():
    """VVRV, 25 mol% TrFE, 48 atoms: DFT on polyfind's cells relaxed in our potential, 0.08 eV/A left at most."""
    r = _stored("vdf3-trfe1-4")
    free, upright, sampler = (r["margins"][k] for k in ("free", "upright", "sampler"))
    assert free["dft"] == pytest.approx(0.678, abs=0.002) and free["largest_force_left_eV_per_A"] < 0.09
    assert upright["dft"] == pytest.approx(0.505, abs=0.002)
    # at the cells as polyfind has them the forces are 1.3 eV/A and the gap is not the margin
    assert free["dft_at_the_starts"] == pytest.approx(0.438, abs=0.002)
    for name, off in (("chain1", -0.007), ("belt4", -0.037)):
        assert r["polyfind"][name]["margin"] - free["dft"] == pytest.approx(off, abs=0.003)
        assert abs(r["polyfind"][name]["margin"] - free["dft"]) < NEEDED_KCAL
    # the two presets split the same margin differently: the helix's gain from shearing is -0.17 in DFT
    gain = free["dft"] - upright["dft"]
    assert gain == pytest.approx(0.174, abs=0.003)
    assert -r["polyfind"]["chain1"]["helix_shear_gain"] == pytest.approx(0.250, abs=0.003)
    assert -r["polyfind"]["belt4"]["helix_shear_gain"] == pytest.approx(0.064, abs=0.003)
    # polyfind's arrangement of the chains is the lower one in DFT, in both conformations
    last = {k: v["steps"][-1]["eV_per_monomer"] for k, v in r["lines"].items()}
    assert last["row_vdf3-trfe1-4_TTTT_chain1"] < last["units_24_0"] and last["row_vdf3-trfe1-4_TGpTG_chain1"] < last["units_25_0"]
    assert sampler["dft"] == pytest.approx(0.793, abs=0.002)
    # our potential found the geometries (0.08 eV/A left); its own margins are 0.07 to 0.10 below DFT's, at the edge of the needed
    for m in (free, upright, sampler):
        assert -0.10 < m["potential"] - m["dft"] < -0.065


def test_neither_preset_has_the_chlorine_cell():
    """VRVC, one CFE and one TrFE in four, 48 atoms: both presets' cells relax to the same two cells in the potential."""
    r = _stored("vdf2-trfe1-cfe1-4")
    m = r["margins"]["free"]
    assert m["dft"] == pytest.approx(0.682, abs=0.002) and m["largest_force_left_eV_per_A"] < 0.07
    assert m["dft_at_the_starts"] is None  # polyfind's own cells were not labelled here, only the relaxed ones
    assert r["polyfind"]["chain1"]["margin"] - m["dft"] == pytest.approx(-0.301, abs=0.003)
    assert r["polyfind"]["belt4"]["margin"] - m["dft"] == pytest.approx(+0.375, abs=0.003)
    assert m["potential"] - m["dft"] == pytest.approx(+0.086, abs=0.003)
    # the unit's effect in this cell, against the TrFE cell with a VDF in its place: nothing in DFT, large and opposite in the presets
    t = _stored("vdf3-trfe1-4")
    assert m["dft"] - t["margins"]["free"]["dft"] == pytest.approx(0.004, abs=0.003)
    assert r["polyfind"]["belt4"]["margin"] - t["polyfind"]["belt4"]["margin"] == pytest.approx(+0.416, abs=0.003)
    assert r["polyfind"]["chain1"]["margin"] - t["polyfind"]["chain1"]["margin"] == pytest.approx(-0.289, abs=0.003)


def test_a_relaxation_line_is_read_in_the_order_it_was_walked():
    pytest.importorskip("ase")
    spec = importlib.util.spec_from_file_location("row_cells_frames", ROOT / "examples" / "row_cells_frames.py")
    rcf = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(rcf)
    from ase import Atoms

    def frame(name):
        a = Atoms("H")
        a.info["frame"] = name
        return a

    lines = rcf.lines_of([frame(n) for n in ("a_in_p_in_p", "b", "a_in_p", "a", "b_in_p")])
    assert [x.info["frame"] for x in lines["a"]] == ["a", "a_in_p", "a_in_p_in_p"]
    assert [x.info["frame"] for x in lines["b"]] == ["b", "b_in_p"]
    assert rcf.safe("TG'TG") == "TGpTG"


def _splits():
    dft = json.loads((PRED / "chain_split_dft_rows_pbe_d3bj.json").read_text())["pairs"]
    poly = json.loads((PRED / "chain_split_polyfind_dft_rows.json").read_text())["sets"]
    t = "row_vdf3-trfe1-4_TTTT_chain1_in_all_s1 minus row_vdf3-trfe1-4_TGpTG_chain1_upright_in_all_s1"
    c = "row_vdf2-trfe1-cfe1-4_TTTT_chain1_in_all_s1 minus row_vdf2-trfe1-cfe1-4_GpTGT_chain1_in_all_s1"
    own = {name: (s["vdf3-trfe1-4"]["pairs"]["TTTT minus TG'TG"], s["vdf2-trfe1-cfe1-4"]["pairs"]["TTTT minus G'TGT"]) for name, s in poly.items()
           if "vdf3-trfe1-4" in s and "vdf2-trfe1-cfe1-4" in s}  # a set made for another unit holds that unit's rows only
    return dft[t], dft[c], own


def test_the_two_cells_split_alike_in_dft():
    """A lone chain cut from each relaxed cell: the margins are three inside a chain less two and a third between chains."""
    trfe, cfe, _ = _splits()
    assert trfe["inside_a_chain"] == pytest.approx(2.935, abs=0.003) and trfe["between_chains"] == pytest.approx(-2.430, abs=0.003)
    assert cfe["inside_a_chain"] == pytest.approx(2.979, abs=0.003) and cfe["between_chains"] == pytest.approx(-2.298, abs=0.003)
    # what the CFE unit does in DFT, against a VDF in its place: little inside the chain, little between chains
    assert cfe["inside_a_chain"] - trfe["inside_a_chain"] == pytest.approx(+0.044, abs=0.004)
    assert cfe["between_chains"] - trfe["between_chains"] == pytest.approx(+0.133, abs=0.004)


def test_the_chlorine_error_of_both_presets_is_inside_the_chain():
    trfe, cfe, own = _splits()
    unit = {k: cfe[k] - trfe[k] for k in ("inside_a_chain", "between_chains")}
    errors = {}
    for name, (t, c) in own.items():
        errors[name] = {k: (c[k] - t[k]) - unit[k] for k in unit}
    # between chains the unit is right in both (0.06 and 0.02); inside the chain chain1 has it 0.25 too low and belt4 0.28 too high
    assert errors["chain1"]["between_chains"] == pytest.approx(+0.059, abs=0.005) and errors["belt4"]["between_chains"] == pytest.approx(+0.022, abs=0.005)
    assert errors["chain1"]["inside_a_chain"] == pytest.approx(-0.249, abs=0.005) and errors["belt4"]["inside_a_chain"] == pytest.approx(+0.284, abs=0.005)
    # belt4 on the TrFE cell alone: two errors of 0.7 that cancel, as on the reversed-unit candidate; chain1's are 0.13 and 0.22
    t4, tc = own["belt4"][0], own["chain1"][0]
    assert t4["inside_a_chain"] - trfe["inside_a_chain"] == pytest.approx(+0.729, abs=0.005) and t4["between_chains"] - trfe["between_chains"] == pytest.approx(-0.657, abs=0.005)
    assert tc["inside_a_chain"] - trfe["inside_a_chain"] == pytest.approx(+0.134, abs=0.005) and tc["between_chains"] - trfe["between_chains"] == pytest.approx(-0.218, abs=0.005)


def test_the_ctfe_cell_is_outside_both_presets_too():
    """VRVK, one CTFE and one TrFE in four: the unit's CFCl-CH2 bond is the CFE unit's, its CF2-CF2 bond the reversed unit's."""
    r = _stored("vdf2-trfe1-ctfe1-4")
    m = r["margins"]["free"]
    assert m["dft"] == pytest.approx(0.315, abs=0.002) and m["largest_force_left_eV_per_A"] < 0.12
    assert r["polyfind"]["chain1"]["margin"] - m["dft"] == pytest.approx(-0.181, abs=0.003)
    assert r["polyfind"]["belt4"]["margin"] - m["dft"] == pytest.approx(+0.419, abs=0.003)
    split = next(iter(json.loads((PRED / "chain_split_dft_vdf2-trfe1-ctfe1-4.json").read_text())["pairs"].values()))
    # a third of the CFE cell's margin inside a chain, a third of it between chains: the unit folds the chain and parts the chains far less
    assert split["inside_a_chain"] == pytest.approx(1.091, abs=0.003) and split["between_chains"] == pytest.approx(-0.776, abs=0.003)
    assert split["cell"] == pytest.approx(m["dft"], abs=1e-4)


def test_the_cfcl_bonds_triple_is_fitted_to_lone_chains():
    """The alternating VDF-CFE chain in DFT: all-trans, and the helix in its two placements (chlorine on either side of the fold)."""
    from polyfind import fitting as F

    conf = json.loads((PRED / "chain_conformers_cfe_pbe_d3bj.json").read_text())["conformers"]
    assert sorted(conf) == ["cfe2_helix", "cfe2_helix2", "cfe2_tt"]
    assert all(c["largest_force"] < 0.01 and c["bonds"] == ["CF2-CH2", "CH2-CFCl", "CFCl-CH2", "CH2-CF2"] for c in conf.values())
    assert conf["cfe2_helix"]["energy"] == pytest.approx(-2.475, abs=0.002) and conf["cfe2_helix2"]["energy"] == pytest.approx(-3.189, abs=0.002)
    # the CFCl-CH2 bond is the gauche one in both, one way and the other; the CH2-CFCl bond stays trans
    assert conf["cfe2_helix"]["dihedrals"][2] == pytest.approx(-60, abs=2) and conf["cfe2_helix2"]["dihedrals"][2] == pytest.approx(+50, abs=2)
    assert all(abs(abs(c["dihedrals"][1]) - 180) < 5 for k, c in conf.items() if k != "cfe2_tt")
    fit = json.loads((PRED / "unit_torsion_cfe_chain1.json").read_text())
    assert fit["base"] == "pvdf-dft-chain1" and fit["fitted"] == ["cfe2_helix", "cfe2_helix2"]
    rows = fit["rows"]
    # chain1 has both helices too high, by 0.69 and 0.35 a monomer; the triple brings both inside 0.11
    assert rows["cfe2_helix"]["at_the_base"] - rows["cfe2_helix"]["dft"] == pytest.approx(+0.692, abs=0.003)
    assert rows["cfe2_helix2"]["at_the_base"] - rows["cfe2_helix2"]["dft"] == pytest.approx(+0.355, abs=0.003)
    assert rows["cfe2_helix"]["polyfind"] - rows["cfe2_helix"]["dft"] == pytest.approx(+0.073, abs=0.003)
    assert rows["cfe2_helix2"]["polyfind"] - rows["cfe2_helix2"]["dft"] == pytest.approx(-0.104, abs=0.003)
    # the preset is chain1 with that triple on the CH2-CFCl bonds, nothing else changed
    preset, base = F.PRESET_PARAMETERS["pvdf-dft-chain1-cfe"], F.PRESET_PARAMETERS["pvdf-dft-chain1"]
    assert preset == F.parameters_of(PRED / "unit_torsion_cfe_chain1.json")
    assert preset == F.typed_on(base, {"CH2-CFCl": F.UNIT_TORSION_CFCL}) and preset.torsion == base.torsion
    assert dict((t, v) for t, *v in preset.torsion_types)["CH2-CFCl"] == pytest.approx(fit["torsion_types"]["CH2-CFCl"])
    assert [t for t, *_ in preset.torsion_types] == [t for t, *_ in base.torsion_types] + ["CH2-CFCl"]


def test_the_cfcl_triple_meets_dft_on_cells_it_was_not_fitted_to():
    under = json.loads((PRED / "margins_chain1-cfe_dft_rows.json").read_text())
    assert under["preset"] == "pvdf-dft-chain1-cfe"
    rows = under["rows"]
    cfe, ctfe, trfe = (_stored(k)["margins"]["free"]["dft"] for k in ("vdf2-trfe1-cfe1-4", "vdf2-trfe1-ctfe1-4", "vdf3-trfe1-4"))
    # the CFE cell: +0.690 where DFT has +0.682 and chain1 +0.381, in DFT's helix placement
    assert rows["vdf2-trfe1-cfe1-4"]["margin"] - cfe == pytest.approx(+0.008, abs=0.003)
    assert rows["vdf2-trfe1-cfe1-4"]["helix_sheared"] == "G'TGT"
    # the CTFE cell, whose CFCl-CH2 bond takes the same triple: +0.400 against +0.315, inside the needed 0.1 and not by much
    assert rows["vdf2-trfe1-ctfe1-4"]["margin"] - ctfe == pytest.approx(+0.085, abs=0.003)
    assert rows["vdf2-trfe1-ctfe1-4"]["helix_sheared"] == "G'TGT"
    assert all(abs(rows[k]["margin"] - d) < NEEDED_KCAL for k, d in (("vdf2-trfe1-cfe1-4", cfe), ("vdf2-trfe1-ctfe1-4", ctfe), ("vdf3-trfe1-4", trfe)))
    # a cell with no CFCl carbon is chain1's to the last digit
    plain = json.loads((PRED / "margins_chain1_dft_rows.json").read_text())["rows"]["vdf3-trfe1-4"]
    assert rows["vdf3-trfe1-4"]["margin"] == pytest.approx(plain["margin"], abs=1e-9)
    # the split of the CFE cell: what the unit adds to the TrFE cell is now right inside the chain (0.05 off, from 0.25) and between chains (0.07)
    trfe_s, cfe_s, own = _splits()
    t, c = own["chain1-cfe"]
    assert (c["inside_a_chain"] - t["inside_a_chain"]) - (cfe_s["inside_a_chain"] - trfe_s["inside_a_chain"]) == pytest.approx(+0.048, abs=0.005)
    assert (c["between_chains"] - t["between_chains"]) - (cfe_s["between_chains"] - trfe_s["between_chains"]) == pytest.approx(+0.070, abs=0.005)


def test_the_nitrile_cell_has_all_trans_below_every_helix_and_the_presets_miss_it_between_chains():
    """VRVA, one acrylonitrile and one TrFE in four: DFT on the all-trans cell and three helix placements."""
    r = _stored("vdf2-trfe1-an1-4")
    m = r["margins"]
    got = {r["lines"][v["helix"]]["conformation"]: v["dft"] for v in m.values()}
    assert got == pytest.approx({"TGTG'": -0.386, "TG'TG": -0.404, "G'TGT": -0.678}, abs=0.002)
    assert max(v["largest_force_left_eV_per_A"] for v in m.values()) < 0.16
    # both presets have a helix cell below the all-trans one: the sign is wrong, by 0.64 and 0.74
    assert r["polyfind"]["chain1"]["margin"] - m["free"]["dft"] == pytest.approx(+0.644, abs=0.003)
    assert r["polyfind"]["belt4"]["margin"] - m["free"]["dft"] == pytest.approx(+0.739, abs=0.003)
    # the potential, which never saw this cell, has the sign and is within 0.08 of DFT on the lowest placement
    assert m["free"]["potential"] - m["free"]["dft"] == pytest.approx(-0.080, abs=0.003)
    split = next(iter(json.loads((PRED / "chain_split_dft_vdf2-trfe1-an1-4.json").read_text())["pairs"].values()))
    assert split["inside_a_chain"] == pytest.approx(2.128, abs=0.003) and split["between_chains"] == pytest.approx(-2.514, abs=0.003)
    poly = json.loads((PRED / "chain_split_polyfind_dft_rows.json").read_text())["sets"]
    off = {n: {k: poly[n]["vdf2-trfe1-an1-4"]["pairs"]["TTTT minus TGTG'"][k] - split[k] for k in ("inside_a_chain", "between_chains")} for n in ("chain1", "belt4")}
    # the error is between chains in both (0.59 and 0.37), where the CFE unit's was inside the chain: a torsion term does not reach it
    assert off["chain1"]["between_chains"] == pytest.approx(+0.590, abs=0.005) and off["chain1"]["inside_a_chain"] == pytest.approx(-0.144, abs=0.005)
    assert off["belt4"]["between_chains"] == pytest.approx(+0.365, abs=0.005) and off["belt4"]["inside_a_chain"] == pytest.approx(-0.277, abs=0.005)
