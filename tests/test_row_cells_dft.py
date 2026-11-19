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


def test_the_second_nitrile_cell_and_the_nitrile_fits_that_do_not_pass_it():
    """VVVA, one nitrile in four and no TrFE: the cell no nitrile fit has seen.  Two fits made on VRVA are stored and neither is a preset."""
    from polyfind import fitting as F

    r = _stored("vdf3-an1-4")
    m = r["margins"]["free"]
    assert m["dft"] == pytest.approx(-0.091, abs=0.002) and m["largest_force_left_eV_per_A"] < 0.2
    assert m["potential"] - m["dft"] == pytest.approx(+0.083, abs=0.003)
    # chain1 is 0.20 low and on another helix placement than DFT's (on DFT's it is 0.66 low); belt4 is 0.15 high
    assert r["polyfind"]["chain1"]["margin"] - m["dft"] == pytest.approx(-0.201, abs=0.003) and r["polyfind"]["chain1"]["helix_sheared"] != "TG'TG"
    assert r["polyfind"]["belt4"]["margin"] - m["dft"] == pytest.approx(+0.149, abs=0.003)
    split = next(iter(json.loads((PRED / "chain_split_dft_vdf3-an1-4.json").read_text())["pairs"].values()))
    assert split["inside_a_chain"] == pytest.approx(0.878, abs=0.003) and split["between_chains"] == pytest.approx(-0.969, abs=0.003)
    first = _stored("vdf2-trfe1-an1-4")
    dft = {first["lines"][v["helix"]]["conformation"]: v["dft"] for v in first["margins"].values()}
    # 1. the charge and nitrogen's wall fitted to the three energies: within 0.05, by two errors of 0.5 that cancel
    a = json.loads((PRED / "unit_nonbonded_an_margins_chain1.json").read_text())
    assert max(abs(v["polyfind"] - v["dft"]) for v in a["rows"].values() if v["fitted"]) < 0.05
    q = json.loads((PRED / "chain_split_polyfind_dft_rows.json").read_text())["sets"]["an-fitted-to-margins"]["vdf2-trfe1-an1-4"]["pairs"]["TTTT minus TGTG'"]
    s1 = next(iter(json.loads((PRED / "chain_split_dft_vdf2-trfe1-an1-4.json").read_text())["pairs"].values()))
    assert q["inside_a_chain"] - s1["inside_a_chain"] == pytest.approx(+0.500, abs=0.005) and q["between_chains"] - s1["between_chains"] == pytest.approx(-0.593, abs=0.005)
    # 2. fitted to the between-chains parts: 0.63 e, a smaller nitrogen; the parts stay 0.30 and 0.17 off (the cells change basin as the charge moves)
    b = json.loads((PRED / "unit_nonbonded_an_between_chain1.json").read_text())
    assert b["fitted_to"] == "between_chains" and b["charge_increments"]["C-N"] == pytest.approx(0.632, abs=0.001) and b["lj"]["N"] == pytest.approx([3.445, 0.0519], abs=0.001)
    assert {c: round(v["polyfind"] - v["dft"], 2) for c, v in b["rows"].items()} == {"TGTG'": -0.30, "TG'TG": -0.04, "G'TGT": 0.17}
    assert {c: round(v["polyfind"] - v["dft"], 2) for c, v in b["energies_shear_free"].items() if v["dft"] is not None} == {"TGTG'": 0.20, "TG'TG": -0.02, "G'TGT": 0.08}
    # 3. with a torsion triple for the nitrile's bonds on top, fitted to the lone chains (to 0.03, with coefficients of 3): the cells' energies are the test
    t = json.loads((PRED / "unit_torsion_an_chain1.json").read_text())
    assert max(abs(v["polyfind"] - v["dft"]) for v in t["rows"].values()) < 0.03 and max(abs(v) for v in t["torsion_types"]["CH2-CHCN"]) > 3.0
    both = json.loads((PRED / "margins_an_between_torsion_dft_rows.json").read_text())["rows"]
    assert both["vdf2-trfe1-an1-4"]["margin"] - min(dft.values(), key=abs) == pytest.approx(+0.432, abs=0.004)
    assert both["vdf3-an1-4"]["margin"] - m["dft"] == pytest.approx(+0.255, abs=0.004)
    only = json.loads((PRED / "margins_an_between_dft_rows.json").read_text())["rows"]
    assert only["vdf3-an1-4"]["margin"] - m["dft"] == pytest.approx(-0.458, abs=0.004)
    # so neither is accepted: no preset has a term for the nitrile
    for name in ("pvdf-dft-chain1", "pvdf-dft-chain1-cfe"):
        prm = F.PRESET_PARAMETERS[name]
        assert "N" not in {e for e, *_ in prm.lj} and "CH2-CHCN" not in {t for t, *_ in prm.torsion_types}
        assert dict(((a, b), v) for a, b, v in prm.charge_increments)[("C", "N")] == pytest.approx(0.370, abs=0.001)


def test_chain1_has_the_cdfe_cells_margin_and_not_its_parts():
    """VRVD, one CDFE and one TrFE in four: the sum is met to 0.01 with the two parts 0.7 off each way, which is not a pass."""
    r = _stored("vdf2-trfe1-cdfe1-4")
    m = r["margins"]["free"]
    assert m["dft"] == pytest.approx(0.968, abs=0.002) and m["largest_force_left_eV_per_A"] < 0.23
    assert r["polyfind"]["chain1"]["margin"] - m["dft"] == pytest.approx(+0.007, abs=0.003)
    assert r["polyfind"]["belt4"]["margin"] - m["dft"] == pytest.approx(+0.752, abs=0.003)
    split = next(iter(json.loads((PRED / "chain_split_dft_vdf2-trfe1-cdfe1-4.json").read_text())["pairs"].values()))
    assert split["inside_a_chain"] == pytest.approx(2.394, abs=0.003) and split["between_chains"] == pytest.approx(-1.426, abs=0.003)
    own = json.loads((PRED / "chain_split_polyfind_dft_rows.json").read_text())["sets"]["chain1"]["vdf2-trfe1-cdfe1-4"]["pairs"]["TTTT minus TGTG'"]
    assert own["inside_a_chain"] - split["inside_a_chain"] == pytest.approx(+0.717, abs=0.005)
    assert own["between_chains"] - split["between_chains"] == pytest.approx(-0.791, abs=0.005)
    # one helix placement of this small cell does not build at its ideal torsions, and the polish says so
    polished = json.loads((PRED / "polished_margins_belt4_dft_rows.json").read_text())["rows"]["vdf2-trfe1-cdfe1-4"]
    assert list(polished["not_built"]) == ["TG'TG"] and "TG'TG" not in polished["cells"]


def test_the_cdfe_lone_chains_say_the_cells_margin_is_met_by_cancellation():
    conf = json.loads((PRED / "chain_conformers_cdfe_pbe_d3bj.json").read_text())["conformers"]
    assert all(c["largest_force"] < 0.01 for c in conf.values())
    # the helix one way is 2.8 below the all-trans chain and the other way 3.8 above it
    assert conf["cdfe2_h2"]["energy"] == pytest.approx(-2.843, abs=0.003) and conf["cdfe2_h1"]["energy"] == pytest.approx(+3.750, abs=0.003)
    assert conf["cdfe2_h2"]["bonds"][0] == "CF2-CHCl" and conf["cdfe2_h2"]["dihedrals"][0] == pytest.approx(-62, abs=2)


def test_the_rows_under_chain1_cfe_are_the_25_without_a_nitrile_or_an_unchecked_unit():
    table = json.loads((PRED / "margins_chain1-cfe.json").read_text())
    assert table["preset"] == "pvdf-dft-chain1-cfe"
    rows = table["rows"]
    assert len(rows) == 25 and all("margin" in e and "TTTT" in e["cells"] for e in rows.values())
    # all four helix placements, but for the CDFE rows, whose polished table has two or three
    assert all(len(e["cells"]) == 5 for e in rows.values() if "D" not in e["letters"])
    assert sorted(len(e["cells"]) for e in rows.values() if "D" in e["letters"]) == [3, 3, 3, 4]
    assert not any(set(e["letters"]) & set("AEYX") for e in rows.values())
    assert [k for k, e in rows.items() if "F" in e["letters"]] == ["vdf9-trfe4-tfe1-14"]  # the TFE unit passed its check; its row with a nitrile is not here
    # carried from chain1's own run only where the two presets are the same potential: no CFE, no CTFE
    carried = [k for k, e in rows.items() if e.get("carried_from")]
    assert len(carried) == 8 and not any(set(rows[k]["letters"]) & set("CK") for k in carried)
    m = {k: e["margin"] for k, e in rows.items()}
    assert m["vdf10-trfe4-14"] == pytest.approx(0.222, abs=0.002) and m["vdf9-trfe4-cfe1-14"] == pytest.approx(0.359, abs=0.002)
    assert m["vdf13-trfe6-cfe1-20"] == pytest.approx(0.307, abs=0.002)
    # in the window of 55 V/um on this one axis: four rows, two of them of chemistry that has passed its check cleanly
    assert sorted(k for k, v in m.items() if 0.0 < v < 0.15) == ["vdf6-trfe4", "vdf8-trfe4-ctfe1-rvdf1-14", "vdf9-trfe4-cdfe1-14", "vdf9-trfe4-tfe1-14"]
    assert m["vdf9-trfe4-tfe1-14"] == pytest.approx(0.147, abs=0.002)
    assert sorted(k for k, v in m.items() if v < 0.0) == ["vdf12-rvdf2-14", "vdf5-rvdf1-6", "vdf7-trfe2-cdfe1", "vdf8-trfe4-cfe1-rvdf1-14", "vdf9-trfe4-rvdf1-14"]
    # a row's strain along the chain is stored beside its margin; the other three axes are not here
    assert all(0.08 < e["strain_along_chain"] < 0.13 for e in rows.values())


def test_the_tfe_cell_is_met_by_chain1_and_missed_by_the_potential():
    """VRVF, one TFE and one TrFE in four: all-trans the lower by a little.  The one cell so far on which the potential's own margin is far out."""
    r = _stored("vdf2-trfe1-tfe1-4")
    m = r["margins"]["free"]
    assert m["dft"] == pytest.approx(-0.076, abs=0.002) and m["largest_force_left_eV_per_A"] < 0.08
    assert r["polyfind"]["chain1"]["margin"] - m["dft"] == pytest.approx(+0.016, abs=0.003)
    assert r["polyfind"]["belt4"]["margin"] - m["dft"] == pytest.approx(+0.300, abs=0.003)
    assert m["potential"] - m["dft"] == pytest.approx(-0.358, abs=0.003)
    split = next(iter(json.loads((PRED / "chain_split_dft_vdf2-trfe1-tfe1-4.json").read_text())["pairs"].values()))
    assert split["inside_a_chain"] == pytest.approx(0.895, abs=0.003) and split["between_chains"] == pytest.approx(-0.971, abs=0.003)
    own = json.loads((PRED / "chain_split_polyfind_dft_rows.json").read_text())["sets"]["chain1"]
    q, t = own["vdf2-trfe1-tfe1-4"]["pairs"]["TTTT minus TGTG'"], own["vdf3-trfe1-4"]["pairs"]["TTTT minus TG'TG"]
    # chain1's two parts are 0.26 and 0.29 off, opposite ways: 0.12 and 0.08 more than its own errors on the TrFE cell (0.13, 0.22)
    assert q["inside_a_chain"] - split["inside_a_chain"] == pytest.approx(+0.255, abs=0.005) and q["between_chains"] - split["between_chains"] == pytest.approx(-0.293, abs=0.005)
    trfe = json.loads((PRED / "chain_split_dft_rows_pbe_d3bj.json").read_text())["pairs"]["row_vdf3-trfe1-4_TTTT_chain1_in_all_s1 minus row_vdf3-trfe1-4_TGpTG_chain1_upright_in_all_s1"]
    assert (q["inside_a_chain"] - split["inside_a_chain"]) - (t["inside_a_chain"] - trfe["inside_a_chain"]) == pytest.approx(+0.121, abs=0.006)


def test_the_potentials_own_margin_on_the_seven_small_cells():
    """Within 0.11 kcal/mol per monomer of DFT on six of them and 0.36 off on the TFE cell: an indication, not a referee."""
    off = {}
    for row in ("vdf3-trfe1-4", "vdf2-trfe1-cfe1-4", "vdf2-trfe1-ctfe1-4", "vdf2-trfe1-an1-4", "vdf3-an1-4", "vdf2-trfe1-cdfe1-4", "vdf2-trfe1-tfe1-4"):
        m = _stored(row)["margins"]["free"]
        off[row] = m["potential"] - m["dft"]
    assert max(abs(v) for k, v in off.items() if k != "vdf2-trfe1-tfe1-4") < 0.115
    assert off["vdf2-trfe1-tfe1-4"] == pytest.approx(-0.358, abs=0.003)
    # and it has the sign of every margin
    for row in off:
        m = _stored(row)["margins"]["free"]
        assert m["potential"] * m["dft"] > 0
