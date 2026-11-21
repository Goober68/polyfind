"""DFT at fourteen monomers and every kink site of the four cells in the potential (docs/MLIP_SECOND_OPINION.md 3.21 and 3.22): the stored results."""

import importlib.util
import json
import re
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
PRED = ROOT / "docs" / "predictions"
MLIP = PRED / "mlip"
CFE, LEADER, AN, MATRIX = "vdf9-trfe4-cfe1-14", "vdf9-trfe3-an1-rvdf1-14", "vdf9-trfe4-an1-14", "vdf10-trfe4-14"
KCAL = 23.0605
ROWS_OF_THE_STATES = [   # single points on the cells the potential relaxed
    f"{CFE}=row_{CFE}_TTTT_chain1_in_all_s1,helix_TGTGp:row_{CFE}_TGTGp_chain1_in_all_s1,sheet_at_12:sheet_{CFE}_belt4_at12_sheared_r2,"
    f"sheet_at_14m:sheet_{CFE}_belt4_at14m_sheared_in_all_s1",
    f"{LEADER}=tt_leader_from_rows14an_r1,helix_TGTGp:helix_leader_TGTGp_r1,sheet_at_18m:sheet_{LEADER}_belt4_at18m_sheared_r9,"
    f"sheet_at_18:sheet_{LEADER}_belt4_at18_sheared_in_all_s1,sheet_at_20m:sheet_{LEADER}_belt4_at20m_sheared_in_all_s1,"
    f"sheet_at_16m:sheet_{LEADER}_belt4_at16m_sheared_in_all_s1",
    f"{AN}=tt_nitrile_from_rows14an_r1,helix_GpTGT:row_{AN}_GpTGT_cfe_in_all_s1,sheet_at_4:sheet_{AN}_belt4_at4_sheared_in_all_s1",
    f"{MATRIX}=row_{MATRIX}_TTTT_chain1_in_all_s1,helix_TGTGp:row_{MATRIX}_TGTGp_chain1_in_all_s1,sheet_at_18m:sheet_{MATRIX}_belt4_at18m_sheared_in_all_s1",
]
ROWS_OF_THE_STEPPED = [   # a second DFT calculation on cells stepped once toward DFT's minimum, against the stepped all-trans cells
    f"{CFE}=row_{CFE}_TTTT_chain1_in_all_s1_in_all_s1,helix_TGTGp:row_{CFE}_TGTGp_chain1_in_all_s1_in_all_s1,"
    f"sheet_at_12:sheet_{CFE}_belt4_at12_sheared_r2_in_all_s1,sheet_at_14m:sheet_{CFE}_belt4_at14m_sheared_in_all_s1_in_all_s1",
    f"{LEADER}=tt_leader_from_rows14an_r1_in_all_s1,helix_TGTGp:helix_leader_TGTGp_r1_in_all_s1,sheet_at_18m:sheet_{LEADER}_belt4_at18m_sheared_r9_in_all_s1",
]
DFT_AT_THREE = [f"{CFE}@12=-1.930", f"{CFE}@14m=-5.034", f"{LEADER}@18m=0.338"]
DFT_AT_EIGHT = DFT_AT_THREE + [f"{LEADER}@18=-0.715", f"{LEADER}@20m=-0.400", f"{LEADER}@16m=0.353", f"{AN}@4=-0.431", f"{MATRIX}@18m=0.899"]
# the potential's mean error on a cell's sheets that DFT has, taken off its other sites of that cell: a second reading of the map
MOVED = [f"{CFE}=0.947", f"{LEADER}=1.082", f"{AN}=-0.630", f"{MATRIX}=-0.301"]


def _tool():
    spec = importlib.util.spec_from_file_location("row_cells_frames", ROOT / "examples" / "row_cells_frames.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _states():
    return json.loads((PRED / "cells14_states_dft.json").read_text())["rows"]


def _stepped():
    return json.loads((PRED / "cells14_states_dft_stepped.json").read_text())["rows"]


def _sheets():
    by_row = {}
    for rec in json.loads((PRED / "sheets14_potential_card.json").read_text())["rows"].values():
        by_row.setdefault(rec["frame"].split("_")[1], {})[rec["site"]] = rec
    return by_row


def test_dft_has_the_film_cells_kinked_state_five_kcal_below_all_trans_and_three_below_the_rest():
    film = _states()[CFE]["states"]
    on_the_chlorine, other, helix = film["sheet_at_14m"], film["sheet_at_12"], film["helix_TGTGp"]
    assert [s["dft"] for s in (on_the_chlorine, other, helix)] == pytest.approx([-5.03, -1.93, -1.72], abs=0.006)
    assert [s["potential"] for s in (on_the_chlorine, other, helix)] == pytest.approx([-5.26, -3.60, -2.28], abs=0.006)
    assert on_the_chlorine["dft"] < other["dft"] - 3.0 and on_the_chlorine["dft"] < helix["dft"] - 3.0
    # stepped toward DFT's own minimum it is 5.24 and its step came to rest: the number is firm
    assert on_the_chlorine["stepped"] == pytest.approx(-5.24, abs=0.006) and on_the_chlorine["stepped_settled"]
    # the other sheet and the helix are not told apart: 0.2 on single points, 0.07 at the stepped estimates, and the helix's step did not come to rest
    assert abs(other["dft"] - helix["dft"]) < 0.25 and abs(other["stepped"] - helix["stepped"]) < 0.1
    assert other["stepped_settled"] and not helix["stepped_settled"]
    # the margin as the tables define it (all-trans minus the lowest helix) is a third of the gap to the lowest state
    assert _states()[CFE]["margin"] == pytest.approx(0.123, abs=0.001) and on_the_chlorine["margin_dft_per_monomer"] == pytest.approx(0.360, abs=0.001)


def test_dft_has_the_leaders_first_sheet_level_with_all_trans_and_its_margin_below_zero():
    # the first round's one sheet; the second round found two of its four below all-trans (the test on the nitrile cells below)
    leader = _states()[LEADER]
    helix, sheet = leader["states"]["helix_TGTGp"], leader["states"]["sheet_at_18m"]
    assert helix["dft"] == pytest.approx(0.65, abs=0.006) and helix["stepped"] == pytest.approx(0.74, abs=0.006) and helix["stepped_settled"]
    # the sheet the potential has 1.37 below all-trans: +0.34 on single points, -0.34 stepped; the correction is as large as the number
    assert sheet["potential"] == pytest.approx(-1.37, abs=0.006) and sheet["dft"] == pytest.approx(0.34, abs=0.006)
    assert sheet["stepped"] == pytest.approx(-0.34, abs=0.01) and sheet["stepped_settled"]
    assert abs(sheet["fall_on_stepping"]["state"] - sheet["fall_on_stepping"]["all_trans"]) > abs(sheet["dft"])
    assert min(sheet["dft"], sheet["stepped"]) > sheet["potential"] + 1.0
    assert leader["margin"] == pytest.approx(-0.047, abs=0.001) and leader["margin_stepped"] == pytest.approx(-0.053, abs=0.001) and leader["margin_stepped_settled"]
    assert leader["margin_potential"] == pytest.approx(-0.070, abs=0.001)


def test_the_potential_is_inside_a_tenth_on_the_margins_and_its_error_on_a_sheet_is_not_a_constant():
    states = _states()
    for row in (CFE, LEADER):  # margins: 0.04 and 0.02 kcal/mol per monomer off
        assert abs(states[row]["margin_potential"] - states[row]["margin"]) < 0.05
    for row in (AN, MATRIX):   # the second round's two: 0.03 and 0.05
        assert abs(states[row]["margin_potential"] - states[row]["margin"]) < 0.051
    off = [states[CFE]["states"]["sheet_at_12"]["potential_minus_dft"], states[LEADER]["states"]["sheet_at_18m"]["potential_minus_dft"],
           states[CFE]["states"]["sheet_at_14m"]["potential_minus_dft"]]
    assert off == pytest.approx([-1.67, -1.70, -0.22], abs=0.006)
    assert max(off) - min(off) > 1.4   # two sites agreeing at 1.7 was not an offset
    # per monomer, a sheet's error reaches 0.12: outside the 0.1 needed
    assert max(abs(v) for v in off) / 14 > 0.1


FIRST_ROUND = [f"row_{CFE}_TTTT_chain1_in_all_s1", f"sheet_{CFE}_belt4_at12_sheared_r2", f"row_{CFE}_TGTGp_chain1_in_all_s1", "tt_leader_from_rows14an_r1",
               f"sheet_{LEADER}_belt4_at18m_sheared_r9", "helix_leader_TGTGp_r1", f"sheet_{CFE}_belt4_at14m_sheared_in_all_s1"]


def test_the_labels_are_single_points_in_the_anchors_set_up_and_the_first_seven_carry_the_table():
    from ase.io import read

    frames = {str(a.info["frame"]): a for a in read(MLIP / "cells14_frames_pbe_d3bj.extxyz", index=":")}
    assert list(frames)[:7] == FIRST_ROUND and sorted({len(a) for a in frames.values()}) == [168, 170]
    for a in frames.values():
        assert int(a.info["dft_kpoints"]) == 6 and float(a.info["dft_conv_thr_ry"]) == 1e-9 and "REF_forces" in a.arrays
        # each is a cell the potential brought to rest: 0.004 eV/A, the film cell's all-trans cell 0.005
        assert ("to 0.004 eV/A" in a.info["note"] or "to 0.005 eV/A" in a.info["note"]) and "NOT converged" not in a.info["note"]
        assert 40 < a.info["dft_seconds"] / 60 < 100
    first = [frames[name] for name in FIRST_ROUND]
    assert all(55 < a.info["dft_seconds"] / 60 < 85 for a in first)
    flat, sheet = frames[f"row_{CFE}_TTTT_chain1_in_all_s1"], frames[f"sheet_{CFE}_belt4_at14m_sheared_in_all_s1"]
    per_kink_pair = (sheet.info["REF_energy"] - flat.info["REF_energy"]) / 28 * 14 * KCAL
    assert per_kink_pair == pytest.approx(_states()[CFE]["states"]["sheet_at_14m"]["dft"], abs=1e-6)
    left = [float(np.abs(a.arrays["REF_forces"]).max()) for a in first]
    assert 0.07 < min(left) and max(left) < 0.30   # DFT's forces on the potential's cells


def test_the_stored_tables_are_what_the_tool_writes_from_the_stored_frames(tmp_path):
    tool = _tool()
    for name, stepped, rows in (("cells14_states_dft.json", "cells14_stepped_to_dft.extxyz", ROWS_OF_THE_STATES),
                                ("cells14_states_dft_stepped.json", "cells14_stepped_twice_to_dft.extxyz", ROWS_OF_THE_STEPPED)):
        out = tmp_path / name
        tool.states_dft(SimpleNamespace(frames=[str(MLIP / "cells14_frames_pbe_d3bj.extxyz")], stepped=[str(MLIP / stepped)], rows=rows, out=str(out)))
        assert json.loads(out.read_text()) == json.loads((PRED / name).read_text())


def test_all_four_cells_have_a_margin_in_dft_and_the_potential_has_their_order_and_signs():
    states = _states()
    rows = (CFE, MATRIX, LEADER, AN)
    dft, potential = [states[r]["margin"] for r in rows], [states[r]["margin_potential"] for r in rows]
    assert dft == pytest.approx([0.123, 0.068, -0.047, -0.083], abs=0.001) and potential == pytest.approx([0.163, 0.018, -0.070, -0.109], abs=0.001)
    assert sorted(dft, reverse=True) == dft and sorted(potential, reverse=True) == potential
    assert all(abs(p - d) < 0.051 and (p > 0) == (d > 0) for p, d in zip(potential, dft))
    assert [states[r]["margin_helix"] for r in rows] == ["helix_TGTGp", "helix_TGTGp", "helix_TGTGp", "helix_GpTGT"]
    belt4 = json.loads((PRED / "polished_margins_belt4.json").read_text())["rows"]
    theirs = [belt4[r]["margin"] for r in rows]
    assert theirs == pytest.approx([0.317, 0.056, 0.037, 0.197], abs=0.001)
    assert theirs[3] - dft[3] > 0.27 and theirs[3] > theirs[1] > theirs[2]   # belt4 has the nitrile cell second of the four and DFT has it last
    assert sum(abs(m) < 0.1 for m in dft) == 3                                # three of the four within the 0.1 asked of an evaluator: their signs are not firm
    # one unit in fourteen, in cells that differ from the matrix's by that unit: the chlorine raises the margin, the nitrile lowers it
    assert dft[0] - dft[1] == pytest.approx(0.055, abs=0.001) and dft[3] - dft[1] == pytest.approx(-0.151, abs=0.001)


def test_both_nitrile_cells_have_shallow_kinked_states_below_all_trans_and_the_matrix_none_at_the_site_asked():
    states = _states()
    leader, nitrile, matrix = states[LEADER]["states"], states[AN]["states"], states[MATRIX]["states"]
    sites = ("18", "20m", "18m", "16m")   # all four on the leader's run of five VDF units
    assert [leader[f"sheet_at_{s}"]["dft"] for s in sites] == pytest.approx([-0.72, -0.40, 0.34, 0.35], abs=0.006)
    assert [leader[f"sheet_at_{s}"]["potential"] for s in sites] == pytest.approx([-1.25, -0.93, -1.37, -1.20], abs=0.006)
    # one site on the two hands of the kink: 1.05 apart in DFT, 0.12 in the potential
    assert leader["sheet_at_18m"]["dft"] - leader["sheet_at_18"]["dft"] == pytest.approx(1.05, abs=0.006)
    assert abs(leader["sheet_at_18m"]["potential"] - leader["sheet_at_18"]["potential"]) < 0.13
    assert nitrile["sheet_at_4"]["dft"] == pytest.approx(-0.43, abs=0.006) and nitrile["sheet_at_4"]["potential"] == pytest.approx(0.20, abs=0.006)
    assert nitrile["helix_GpTGT"]["dft"] == pytest.approx(1.16, abs=0.006) and leader["helix_TGTGp"]["dft"] == pytest.approx(0.65, abs=0.006)
    assert matrix["helix_TGTGp"]["dft"] == pytest.approx(-0.95, abs=0.006) and matrix["sheet_at_18m"]["dft"] == pytest.approx(0.90, abs=0.006)
    # the lowest state found, kcal/mol per monomer against all-trans: the film cell's is seven times the leader's
    deepest = [min(s["dft"] for s in states[row]["states"].values()) / 14 for row in (CFE, LEADER, AN, MATRIX)]
    assert deepest == pytest.approx([-0.360, -0.051, -0.031, -0.068], abs=0.001)
    assert min(states[MATRIX]["states"], key=lambda k: matrix[k]["dft"]) == "helix_TGTGp"   # the matrix's lowest state is its helix
    # the potential on the eight sheets DFT has: -1.71 to +0.63 per kink pair, rms 1.07, the sign against all-trans right at five
    sheets = [s for row in states.values() for label, s in row["states"].items() if label.startswith("sheet_at_")]
    off = [s["potential_minus_dft"] for s in sheets]
    assert len(off) == 8 and min(off) == pytest.approx(-1.71, abs=0.006) and max(off) == pytest.approx(0.63, abs=0.006)
    assert float(np.mean(off)) == pytest.approx(-0.66, abs=0.006) and float(np.sqrt(np.mean(np.square(off)))) == pytest.approx(1.07, abs=0.006)
    assert sum((s["potential"] < 0) == (s["dft"] < 0) for s in sheets) == 5
    # the helix cells: every one within 0.05 per monomer
    helix = [s["potential_minus_dft"] / 14 for row in states.values() for label, s in row["states"].items() if label.startswith("helix")]
    assert len(helix) == 4 and max(abs(v) for v in helix) < 0.051


def test_dft_on_the_stepped_cells_holds_the_film_cells_states_and_both_margins_and_levels_the_leaders_sheet():
    from ase.io import read

    stepped, single = _stepped(), _states()
    film, leader = stepped[CFE]["states"], stepped[LEADER]["states"]
    assert [film[k]["dft"] for k in ("sheet_at_14m", "sheet_at_12", "helix_TGTGp")] == pytest.approx([-5.27, -2.04, -1.70], abs=0.006)
    assert film["sheet_at_12"]["dft"] < film["helix_TGTGp"]["dft"] - 0.3           # told apart, in the single points' order (0.21 there)
    assert leader["sheet_at_18m"]["dft"] == pytest.approx(0.06, abs=0.006) and leader["helix_TGTGp"]["dft"] == pytest.approx(0.68, abs=0.006)
    assert stepped[CFE]["margin"] == pytest.approx(0.121, abs=0.001) and stepped[LEADER]["margin"] == pytest.approx(-0.048, abs=0.001)
    for row in (CFE, LEADER):   # a margin moves by 0.002 at most when DFT is let move its two cells
        assert abs(stepped[row]["margin"] - single[row]["margin"]) < 0.003
    # a sheet falls 0.1 to 0.3 further than its all-trans cell
    moved = [film["sheet_at_14m"]["dft"] - single[CFE]["states"]["sheet_at_14m"]["dft"], film["sheet_at_12"]["dft"] - single[CFE]["states"]["sheet_at_12"]["dft"],
             leader["sheet_at_18m"]["dft"] - single[LEADER]["states"]["sheet_at_18m"]["dft"]]
    assert moved == pytest.approx([-0.23, -0.11, -0.27], abs=0.006)
    cells = [rec["all_trans"] for rec in stepped.values()] + [s["cell"] for rec in stepped.values() for s in rec["states"].values()]
    assert len(cells) == 7 and all(c["stepped_from"] for c in cells) and max(c["dft_largest_force_eV_per_A"] for c in cells) < 0.10
    assert max(abs(v) for c in cells for v in c["dft_stress_xx_yy_zz_GPa"]) < 0.06
    # the step's predicted fall against the second calculation's: within 0.04 on four cells of seven, 0.10 to 0.42 too large on three
    frames = {str(a.info["frame"]): a for a in read(MLIP / "cells14_frames_pbe_d3bj.extxyz", index=":")}
    off = {}
    for name, a in frames.items():
        if "start_frame" in a.info:
            before = frames[str(a.info["start_frame"])].info["REF_energy"] / 28
            off[name] = (float(a.info["dft_predicted_eV_per_monomer"]) - a.info["REF_energy"] / 28) * 14 * KCAL
            assert a.info["REF_energy"] / 28 < before
    assert len(off) == 7 and sorted(abs(v) for v in off.values())[:4] == pytest.approx([0.001, 0.024, 0.025, 0.038], abs=0.002)
    too_large = sorted(v for v in off.values() if abs(v) > 0.05)
    assert too_large == pytest.approx([-0.423, -0.376, -0.099], abs=0.002)
    assert off[f"sheet_{CFE}_belt4_at14m_sheared_in_all_s1_in_all_s1"] == pytest.approx(0.0, abs=0.002)   # exact on the round's main number


def test_the_card_labels_what_the_processors_label():
    check = json.loads((PRED / "qe_gpu_check_a40.json").read_text())
    assert check["agrees"] and check["atoms"] == 168 and check["anchor"]["reproduced"]
    assert abs(check["energy_difference_ry_per_cell"]) < 1e-7 and check["largest_force_difference_eV_per_A"] < 1e-5
    assert check["seconds_on_the_card"] / 60 == pytest.approx(45, abs=1) and check["seconds_on_the_processors"] / 60 == pytest.approx(69, abs=1)
    from ase.io import read

    second = read(MLIP / "cells14_frames_pbe_d3bj.extxyz", index="7:")
    minutes = {run: sorted(a.info["dft_seconds"] / 60 for a in second if a.info["dft_run"] == run) for run in ("processors", "card")}
    assert len(minutes["processors"]) == 7 and len(minutes["card"]) == 9
    assert 60 < minutes["processors"][0] and minutes["processors"][-1] < 75 and 42 < minutes["card"][0] and minutes["card"][-1] < 52


def test_every_kink_site_of_the_four_cells_in_the_potential():
    sheets = _sheets()
    assert {row: len(sites) for row, sites in sheets.items()} == {CFE: 28, LEADER: 26, AN: 26, MATRIX: 28}
    kept = {row: {s: r for s, r in sites.items() if r["kink_kept"]} for row, sites in sheets.items()}
    assert {row: len(sites) for row, sites in kept.items()} == {CFE: 27, LEADER: 24, AN: 25, MATRIX: 28}
    below = {row: sorted(s for s, r in sites.items() if r["potential"] < 0) for row, sites in kept.items()}
    assert below[CFE] == ["10", "12", "14", "14m"] and len(below[LEADER]) == 7 and below[AN] == [] and below[MATRIX] == []
    lowest = {row: min(sites, key=lambda s: sites[s]["potential"]) for row, sites in kept.items()}
    assert lowest == {CFE: "14m", LEADER: "18m", AN: "4", MATRIX: "18m"}
    assert [kept[row][lowest[row]]["potential"] for row in (CFE, LEADER, AN, MATRIX)] == pytest.approx([-5.26, -1.32, 0.20, 1.20], abs=0.006)
    # belt4 names another site its cheapest in every one of the four
    for row, sites in kept.items():
        assert min(sites, key=lambda s: sites[s]["polyfind"]) != lowest[row]
    assert kept[CFE]["14m"]["polyfind"] == pytest.approx(-1.17, abs=0.006) and kept[LEADER]["18m"]["polyfind"] == pytest.approx(1.99, abs=0.006)


def test_the_card_finds_the_sheets_this_machines_processors_found():
    sheets = _sheets()
    before = json.loads((PRED / "sheets14_potential_sheared.json").read_text())["rows"]
    # from belt4's own starts in one run each: the leader's two known sheets and the film cell's, within the flat bottom's 0.05
    assert sheets[LEADER]["18m"]["potential"] == pytest.approx(before[f"{LEADER}@18m"]["potential"], abs=0.05)
    assert sheets[LEADER]["10"]["potential"] == pytest.approx(before[f"{LEADER}@10"]["potential"], abs=0.01)
    assert sheets[CFE]["12"]["potential"] == pytest.approx(before[f"{CFE}@12"]["potential"], abs=0.05)
    assert sheets[LEADER]["18m"]["shear_along_chains_monomers"] == pytest.approx(before[f"{LEADER}@18m"]["shear_along_chains_monomers"], abs=0.03)
    check = json.loads((PRED / "card_check_a40.json").read_text())
    assert check["device"] == "cuda" and len(check["cells"]) == 18 and check["largest_difference_meV_per_monomer"] < 1e-6


def test_the_124_cells_all_came_to_rest_at_the_tight_stop_and_the_queue_covered_each_start_once():
    from ase.io import read

    relaxed = [a.info for a in read(MLIP / "cells14_tight_in_all_s1.extxyz", index=":")]
    assert len(relaxed) == 124
    for info in relaxed:
        assert "to 0.004 eV/A" in info["note"] and "NOT converged" not in info["note"]
    starts = [str(a.info["frame"]) for a in read(MLIP / "cells14_starts_for_all_s1.extxyz", index=":")]
    assert sorted(f"{s}_in_all_s1" for s in starts) == sorted(str(i["frame"]) for i in relaxed)
    hits = {s: 0 for s in starts}
    for line in (ROOT / "examples" / "gpu_jobs" / "queue_20261120_dft14.txt").read_text().splitlines():
        if line.startswith("relax_card.sh"):
            for s in starts:
                hits[s] += bool(re.fullmatch(line.split()[2], s))
    assert set(hits.values()) == {1}


def test_the_potentials_margins_from_all_four_helix_placements():
    card = json.loads((PRED / "rows14_potential_card.json").read_text())["rows"]
    tight = json.loads((PRED / "rows14_potential_tight.json").read_text())["rows"]
    assert {row: rec["helix_placements"] for row, rec in card.items()} == {MATRIX: 4, CFE: 4, AN: 4, LEADER: 4}
    assert [card[row]["margin"] for row in (MATRIX, CFE, AN, LEADER)] == pytest.approx([0.018, 0.163, -0.109, -0.070], abs=0.001)
    # only the film cell's moves: its lowest placement at the tight stop is not the one the loose stop had lowest
    assert card[CFE]["helix"] in ("TGTG'", "GTG'T") and tight[CFE]["helix"] == "G'TGT"
    assert card[CFE]["all_trans_minus_helix"]["G'TGT"] == pytest.approx(tight[CFE]["margin"], abs=0.002)
    for row in (MATRIX, AN, LEADER):
        assert card[row]["margin"] == pytest.approx(tight[row]["margin"], abs=0.003)


def test_the_four_axes_at_dfts_margins_give_the_film_cell_a_row_and_the_leader_none():
    as_it_is = json.loads((PRED / "four_axes_dft_margins_card_sheets.json").read_text())
    lifted = json.loads((PRED / "four_axes_dft_margins_card_sheets_lifted.json").read_text())
    states = _states()
    for stored in (as_it_is, lifted):
        assert stored["replaced"] == pytest.approx({CFE: states[CFE]["margin"], LEADER: states[LEADER]["margin"]})
        film = stored["shifts"]["+0.00"]["rows"][CFE]["by_the_sum_over_kink_sites"]
        assert film["f_relax_hz"] == pytest.approx(364e3, rel=0.001) and film["tan_delta_34kHz"] == pytest.approx(0.065, abs=0.001)
        assert film["strain_ac_at_the_heat_budget_34kHz"] == pytest.approx(0.00073, abs=1e-5) and film["blocked_pn_nm2_through_the_crystal"] == pytest.approx(6.4, abs=0.05)
        assert stored["shifts"]["+0.00"]["rows"][CFE]["kink_sites"]["cheapest_kcal_per_kink_pair"] == pytest.approx(-5.034)   # DFT's own number at that site
        leader = stored["shifts"]["+0.00"]["rows"][LEADER]
        assert leader["all_trans_is_lower"] and leader["by_the_sum_over_kink_sites"]["margin_the_model_used"] == 0.0
    # the leader's readings span three orders of magnitude in frequency: the loss model outside its domain, not a row
    readings = [stored["shifts"]["+0.00"]["rows"][LEADER][model] for stored in (as_it_is, lifted) for model in ("by_the_sum_over_kink_sites", "by_the_occupancy_of_kink_sites")]
    hz = [r["f_relax_hz"] for r in readings]
    assert min(hz) < 300 and max(hz) == pytest.approx(415e3, rel=0.01) and max(hz) / min(hz) > 1000
    film = as_it_is["shifts"]["+0.00"]["rows"][CFE]["by_the_sum_over_kink_sites"]
    for r in readings:  # in no reading has the leader the film cell's strain or its blocked force
        assert r["strain_ac_at_the_heat_budget_34kHz"] < film["strain_ac_at_the_heat_budget_34kHz"]
        assert r["blocked_pn_nm2_through_the_crystal"] < 0.6 * film["blocked_pn_nm2_through_the_crystal"]


def test_at_dfts_four_margins_only_the_film_cells_row_holds_under_both_readings_of_the_kink_sites():
    states = _states()
    tables = [json.loads((PRED / name).read_text()) for name in ("four_axes_dft_margins_four_cells.json", "four_axes_dft_margins_four_cells_moved.json")]
    readings = ("by_the_sum_over_kink_sites", "by_the_occupancy_of_kink_sites")
    strain, force = "strain_ac_at_the_heat_budget_34kHz", "blocked_pn_nm2_through_the_crystal"
    for stored in tables:
        assert stored["replaced"] == pytest.approx({row: states[row]["margin"] for row in (CFE, LEADER, AN, MATRIX)})
        rows = stored["shifts"]["+0.00"]["rows"]
        assert [rows[row]["all_trans_is_lower"] for row in (CFE, MATRIX, LEADER, AN)] == [False, False, True, True]
        film = rows[CFE][readings[0]]
        assert rows[CFE][readings[1]] == film                                  # the readings are calibrated on this cell
        assert film["f_relax_hz"] == pytest.approx(364e3, rel=0.001) and film["tan_delta_34kHz"] == pytest.approx(0.065, abs=0.001)
        assert film[strain] == pytest.approx(0.00073, abs=1e-5) and film[force] == pytest.approx(6.4, abs=0.05)
        assert stored["shifts"]["+0.00"]["film_channel_as_calibrated"][strain] == pytest.approx(0.00066, abs=1e-5)
        for row in (LEADER, AN, MATRIX):   # one map, two readings: 260 to 1,100 times apart in frequency
            by_sum, by_share = rows[row][readings[0]], rows[row][readings[1]]
            assert by_sum["f_relax_hz"] < 1500 and 250 < by_share["f_relax_hz"] / by_sum["f_relax_hz"] < 1200
    # one corner of the twelve is ahead of the film cell on three axes: the nitrile cell by the sites' share of time, on the map moved by
    # the potential's error at its one site in DFT.  Its margin is below zero: the model outside its domain.
    film = tables[1]["shifts"]["+0.00"]["rows"][CFE][readings[0]]
    corner = tables[1]["shifts"]["+0.00"]["rows"][AN][readings[1]]
    assert corner["f_relax_hz"] == pytest.approx(534e3, rel=0.01) and corner["tan_delta_34kHz"] < film["tan_delta_34kHz"]
    assert corner[strain] > film[strain] and corner[force] < film[force]
    for i, stored in enumerate(tables):    # in the other eleven: slower, less strain and well under the film cell's blocked force
        for row in (LEADER, AN, MATRIX):
            for reading in readings:
                if (i, row, reading) != (1, AN, readings[1]):
                    r = stored["shifts"]["+0.00"]["rows"][row][reading]
                    assert r["f_relax_hz"] < film["f_relax_hz"] and r[strain] < film[strain] and r[force] < 0.6 * film[force]


def test_a_scan_with_the_potentials_sheets_keeps_the_scans_sites_and_takes_the_energies_given(tmp_path):
    tool = _tool()
    scan = {"rows": {"rowA@0": {"letters": "VV", "sites": {"0": {"formation": 1.0, "units": "VV"}, "2": {"formation": 2.0, "units": "VR"}, "4": {"formation": 3.0, "units": "RV"}}},
                     "rowA@0~": {"letters": "VV", "sites": {"2": {"formation": 5.0, "units": "VR"}}},
                     "rowB@0": {"letters": "VC", "sites": {"0": {"formation": 0.5, "units": "VC"}}}}}
    sheets = {"rows": {"rowA": {"frame": "sheet_rowA_belt4_at0_sheared_in_x", "site": "0", "potential": -1.0, "kink_kept": True, "shear_along_chains_monomers": [0.0, 1.0]},
                       "rowA@2": {"frame": "sheet_rowA_belt4_at2_sheared_in_x", "site": "2", "potential": 0.4, "kink_kept": False, "shear_along_chains_monomers": [0.1, 0.2]},
                       "rowA@2m": {"frame": "sheet_rowA_belt4_at2m_sheared_in_x", "site": "2m", "potential": -0.3, "kink_kept": True, "shear_along_chains_monomers": [0.5, 0.5]}}}
    (tmp_path / "scan.json").write_text(json.dumps(scan))
    (tmp_path / "sheets.json").write_text(json.dumps(sheets))
    tool.sheets_as_scan(SimpleNamespace(scan=str(tmp_path / "scan.json"), sheets=str(tmp_path / "sheets.json"), shift=["rowA=1.5"], set=["rowA@2m=0.25"], out=str(tmp_path / "out.json")))
    out = json.loads((tmp_path / "out.json").read_text())
    assert set(out["rows"]) == {"rowA@0", "rowA@0~"}                       # a row with no sheet is dropped, and so is a site with none
    assert set(out["rows"]["rowA@0"]["sites"]) == {"0", "2"}
    first, lost, other_hand = out["rows"]["rowA@0"]["sites"]["0"], out["rows"]["rowA@0"]["sites"]["2"], out["rows"]["rowA@0~"]["sites"]["2"]
    assert first["formation"] == pytest.approx(0.5) and first["potential"] == -1.0 and first["tables_formation"] == 1.0 and first["units"] == "VV"
    assert lost["kink_kept"] is False and lost["formation"] == pytest.approx(1.9)   # kept in the file, marked, and left out by the rate model
    assert other_hand["formation"] == 0.25 and other_hand["set_from_outside"] is True   # the other hand is the sheets' site 2m; its own number is not shifted
    assert out["shift_kcal_per_kink_pair"] == {"rowA": 1.5} and out["set_kcal_per_kink_pair"] == {"rowA@2m": 0.25}


def test_the_stored_scans_are_what_the_tool_writes(tmp_path):
    tool = _tool()
    for name, shift, known in (("kink_site_scan_potential_card_dft_at_three.json", [], DFT_AT_THREE),
                               ("kink_site_scan_potential_card_dft_at_three_lifted.json", [f"{row}=1.7" for row in (CFE, LEADER, AN, MATRIX)], DFT_AT_THREE),
                               ("kink_site_scan_potential_card_dft_at_eight.json", [], DFT_AT_EIGHT),
                               ("kink_site_scan_potential_card_dft_at_eight_moved.json", MOVED, DFT_AT_EIGHT)):
        out = tmp_path / name
        tool.sheets_as_scan(SimpleNamespace(scan=str(PRED / "kink_site_scan_belt4_shear.json"), sheets=str(PRED / "sheets14_potential_card.json"),
                                            shift=shift, set=known, out=str(out)))
        assert json.loads(out.read_text()) == json.loads((PRED / name).read_text())
        sites = json.loads(out.read_text())["rows"]
        assert sites[f"{CFE}@0~"]["sites"]["14"]["formation"] == -5.034 and sites[f"{LEADER}@0~"]["sites"]["18"]["formation"] == 0.338
    assert sum(len(row["sites"]) for row in sites.values()) == 108
    assert sum(bool(s.get("set_from_outside")) for row in sites.values() for s in row["sites"].values()) == 8


def test_the_sites_given_to_the_scans_are_the_stored_tables_and_the_shifts_are_the_potentials_mean_error():
    states = _states()
    given = {item.split("=")[0]: float(item.split("=")[1]) for item in DFT_AT_EIGHT}
    errors = {}
    for row, rec in states.items():
        for label, state in rec["states"].items():
            if label.startswith("sheet_at_"):
                assert given[f"{row}@{label[len('sheet_at_'):]}"] == pytest.approx(state["dft"], abs=5e-4)
                errors.setdefault(row, []).append(state["potential_minus_dft"])
    assert sum(len(v) for v in errors.values()) == len(given) == 8
    for item in MOVED:
        row, shift = item.split("=")
        assert float(shift) == pytest.approx(-np.mean(errors[row]), abs=2e-3)


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
