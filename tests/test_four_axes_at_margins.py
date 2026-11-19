"""The four axes at another evaluator's margins (examples/four_axes_at_margins.py): the stored table for the potential's."""

import json
from pathlib import Path

import pytest

PRED = Path(__file__).resolve().parents[1] / "docs" / "predictions"
FILM, CFE, LEADER, AN, MATRIX = "film_channel_as_calibrated", "vdf9-trfe4-cfe1-14", "vdf9-trfe3-an1-rvdf1-14", "vdf9-trfe4-an1-14", "vdf10-trfe4-14"
SUM, MOVE = "by_the_sum_over_kink_sites", "by_one_chains_move"


def _stored():
    return json.loads((PRED / "four_axes_potential_margins.json").read_text())


def test_the_margins_are_the_potentials_and_the_rest_is_belt4s():
    stored = _stored()
    potential = json.loads((PRED / "rows14_potential.json").read_text())["rows"]
    assert stored["replaced"] == pytest.approx({k: v["margin"] for k, v in potential.items()})
    assert set(stored["shifts"]) == {"-0.10", "+0.00", "+0.05", "+0.10"}
    belt4 = json.loads((PRED / "pin_design_belt4.json").read_text())
    for entry in stored["shifts"].values():  # the film's calibrated channel does not move with a cell's margin
        film = entry[FILM]
        assert film["strain_ac_at_the_heat_budget_34kHz"] == pytest.approx(belt4["three_axes"]["rows"]["calibrated_film"]["at_the_heat_budget"]["strain_ac"][0])
        for row in (CFE, LEADER, AN, MATRIX):  # a row's frequency comes from the kink energies, not from its margin
            assert entry["rows"][row][SUM]["f_relax_hz"] == pytest.approx(belt4["sites"]["rows"][row]["three_axes"]["f_relax_hz"])


def test_the_cell_of_the_films_composition_lands_on_the_measured_film_over_the_whole_error_bar():
    stored = _stored()
    film = stored["shifts"]["+0.00"][FILM]
    assert (film["tan_delta_34kHz"], 100 * film["strain_ac_at_the_heat_budget_34kHz"]) == pytest.approx((0.060, 0.066), abs=0.001)
    at = stored["shifts"]["+0.00"]["rows"][CFE]
    assert at["margin"] == pytest.approx(+0.078, abs=0.002) and stored["belt4_margins"][CFE] == pytest.approx(+0.317, abs=0.002)
    assert 100 * at[SUM]["strain_ac_at_the_heat_budget_34kHz"] == pytest.approx(0.081, abs=0.001) and at[SUM]["tan_delta_34kHz"] == pytest.approx(0.088, abs=0.001)
    for entry in stored["shifts"].values():
        a = entry["rows"][CFE][SUM]
        assert 0.05 < 100 * a["strain_ac_at_the_heat_budget_34kHz"] < 0.09 and 0.02 < a["tan_delta_34kHz"] < 0.09
        assert a["f_relax_hz"] == pytest.approx(film["f_relax_hz"])   # the anchor: not a test
    # under belt4's margin the same cell had a twentieth of the film's strain
    belt4 = json.loads((PRED / "pin_design_belt4.json").read_text())["sites"]["rows"][CFE]["three_axes"]
    assert 100 * belt4["at_the_heat_budget"]["strain_ac"][0] < 0.005


def test_the_nitrile_cells_are_below_zero_in_the_potential_and_the_model_takes_that_as_zero():
    stored = _stored()
    at = stored["shifts"]["+0.00"]["rows"]
    assert at[LEADER]["all_trans_is_lower"] and at[AN]["all_trans_is_lower"] and at[MATRIX]["all_trans_is_lower"] and not at[CFE]["all_trans_is_lower"]
    assert all(at[r][m]["margin_the_model_used"] == 0.0 for r in (LEADER, AN, MATRIX) for m in (SUM, MOVE))
    # the leader crosses zero inside the bar; the nitrile alone does not
    assert not stored["shifts"]["+0.10"]["rows"][LEADER]["all_trans_is_lower"] and stored["shifts"]["+0.10"]["rows"][AN]["all_trans_is_lower"]
    # if it is not below zero the leader leads on every axis: strain, frequency, loss, and the crystal's reading of the blocked force
    film = stored["shifts"]["+0.00"][FILM]
    for entry in stored["shifts"].values():
        for model in (SUM, MOVE):
            a = entry["rows"][LEADER][model]
            assert 1.9 < a["strain_ac_at_the_heat_budget_34kHz"] / film["strain_ac_at_the_heat_budget_34kHz"] < 2.6
            assert 7.0 < a["f_relax_hz"] / film["f_relax_hz"] < 19.0 and a["tan_delta_34kHz"] < film["tan_delta_34kHz"] / 4.9
            assert a["blocked_pn_nm2_through_the_crystal"] > film["blocked_pn_nm2_through_the_crystal"]
