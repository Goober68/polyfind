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
    own = {name: (s["vdf3-trfe1-4"]["pairs"]["TTTT minus TG'TG"], s["vdf2-trfe1-cfe1-4"]["pairs"]["TTTT minus G'TGT"]) for name, s in poly.items()}
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
