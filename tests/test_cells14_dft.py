"""DFT at fourteen monomers and every kink site of the four cells in the potential (docs/MLIP_SECOND_OPINION.md 3.21): the stored results."""

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
ROWS_OF_THE_STATES = [
    f"{CFE}=row_{CFE}_TTTT_chain1_in_all_s1,helix_TGTGp:row_{CFE}_TGTGp_chain1_in_all_s1,sheet_at_12:sheet_{CFE}_belt4_at12_sheared_r2,"
    f"sheet_at_14m:sheet_{CFE}_belt4_at14m_sheared_in_all_s1",
    f"{LEADER}=tt_leader_from_rows14an_r1,helix_TGTGp:helix_leader_TGTGp_r1,sheet_at_18m:sheet_{LEADER}_belt4_at18m_sheared_r9",
]
DFT_AT_THREE = [f"{CFE}@12=-1.930", f"{CFE}@14m=-5.034", f"{LEADER}@18m=0.338"]


def _tool():
    spec = importlib.util.spec_from_file_location("row_cells_frames", ROOT / "examples" / "row_cells_frames.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _states():
    return json.loads((PRED / "cells14_states_dft.json").read_text())["rows"]


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
    # the other sheet and the helix are not told apart: 0.2 on single points, 0.07 at the stepped estimates, and the helix's step did not come to rest
    assert abs(other["dft"] - helix["dft"]) < 0.25 and abs(other["stepped"] - helix["stepped"]) < 0.1
    assert other["stepped_settled"] and not helix["stepped_settled"]
    # the margin as the tables define it (all-trans minus the lowest helix) is a third of the gap to the lowest state
    assert _states()[CFE]["margin"] == pytest.approx(0.123, abs=0.001) and on_the_chlorine["margin_dft_per_monomer"] == pytest.approx(0.360, abs=0.001)


def test_dft_has_the_leaders_all_trans_cell_lowest_or_level_and_its_margin_below_zero():
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
    off = [states[CFE]["states"]["sheet_at_12"]["potential_minus_dft"], states[LEADER]["states"]["sheet_at_18m"]["potential_minus_dft"],
           states[CFE]["states"]["sheet_at_14m"]["potential_minus_dft"]]
    assert off == pytest.approx([-1.67, -1.70, -0.22], abs=0.006)
    assert max(off) - min(off) > 1.4   # two sites agreeing at 1.7 was not an offset
    # per monomer, a sheet's error reaches 0.12: outside the 0.1 needed
    assert max(abs(v) for v in off) / 14 > 0.1


def test_the_seven_labels_are_single_points_in_the_anchors_set_up_and_carry_the_table():
    from ase.io import read

    frames = {str(a.info["frame"]): a for a in read(MLIP / "cells14_frames_pbe_d3bj.extxyz", index=":")}
    assert len(frames) == 7 and sorted({len(a) for a in frames.values()}) == [168, 170]
    for a in frames.values():
        assert int(a.info["dft_kpoints"]) == 6 and float(a.info["dft_conv_thr_ry"]) == 1e-9 and "REF_forces" in a.arrays
        # each is a cell the potential brought to rest: 0.004 eV/A, the film cell's all-trans cell 0.005
        assert ("to 0.004 eV/A" in a.info["note"] or "to 0.005 eV/A" in a.info["note"]) and "NOT converged" not in a.info["note"]
        assert 55 < a.info["dft_seconds"] / 60 < 85
    flat, sheet = frames[f"row_{CFE}_TTTT_chain1_in_all_s1"], frames[f"sheet_{CFE}_belt4_at14m_sheared_in_all_s1"]
    per_kink_pair = (sheet.info["REF_energy"] - flat.info["REF_energy"]) / 28 * 14 * KCAL
    assert per_kink_pair == pytest.approx(_states()[CFE]["states"]["sheet_at_14m"]["dft"], abs=1e-6)
    left = [float(np.abs(a.arrays["REF_forces"]).max()) for a in frames.values()]
    assert 0.07 < min(left) and max(left) < 0.30   # DFT's forces on the potential's cells


def test_the_stored_table_is_what_the_tool_writes_from_the_stored_frames(tmp_path):
    tool = _tool()
    out = tmp_path / "states.json"
    tool.states_dft(SimpleNamespace(frames=[str(MLIP / "cells14_frames_pbe_d3bj.extxyz")], stepped=[str(MLIP / "cells14_stepped_to_dft.extxyz")],
                                    rows=ROWS_OF_THE_STATES, out=str(out)))
    assert json.loads(out.read_text()) == json.loads((PRED / "cells14_states_dft.json").read_text())


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
    for name, shift in (("kink_site_scan_potential_card_dft_at_three.json", []),
                        ("kink_site_scan_potential_card_dft_at_three_lifted.json", [f"{row}=1.7" for row in (CFE, LEADER, AN, MATRIX)])):
        out = tmp_path / name
        tool.sheets_as_scan(SimpleNamespace(scan=str(PRED / "kink_site_scan_belt4_shear.json"), sheets=str(PRED / "sheets14_potential_card.json"),
                                            shift=shift, set=DFT_AT_THREE, out=str(out)))
        assert json.loads(out.read_text()) == json.loads((PRED / name).read_text())
        sites = json.loads(out.read_text())["rows"]
        assert sites[f"{CFE}@0~"]["sites"]["14"]["formation"] == -5.034 and sites[f"{LEADER}@0~"]["sites"]["18"]["formation"] == 0.338
    assert sum(len(row["sites"]) for row in sites.values()) == 108


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
