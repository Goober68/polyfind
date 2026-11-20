"""The four axes at another evaluator's margins (examples/four_axes_at_margins.py): the stored table for the potential's."""

import json
import math
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


def test_the_kink_sheets_in_the_potential_were_stopped_loose_and_the_lowest_is_the_one_cell_that_sheared():
    """The eleven sheets as stored (stop 0.02 eV/A).  Their ranking is withdrawn (docs/MLIP_SECOND_OPINION.md 3.19): the numbers below are what was found, not minima."""
    batches = [json.loads((PRED / f"{name}.json").read_text())["rows"] for name in ("sheets14_potential", "sheets14_potential_more_sites", "sheets14_potential_third_batch")]
    assert set(batches[0]) == {CFE, LEADER, AN, MATRIX} and [len(b) for b in batches] == [4, 3, 4]
    sheets = [dict(r, row=k.split("@")[0]) for b in batches for k, r in b.items()]
    lost = [r for r in sheets if not r["kink_kept"]]
    assert len(lost) == 1 and lost[0]["row"] == AN and lost[0]["site"] == "10m"   # one of belt4's sites holds no kink in the potential
    kept = [r for r in sheets if r["kink_kept"]]
    # site by site the potential and belt4 differ by -1.3 to +3.0: belt4's profile along the sequence is not confirmed
    diff = sorted(r["potential"] - r["polyfind"] for r in kept)
    assert diff[0] == pytest.approx(-1.28, abs=0.01) and diff[-1] == pytest.approx(+3.03, abs=0.01)
    # the lowest sheet found for each cell: leader < nitrile alone < matrix < the film's cell
    low = {row: min(r["potential"] for r in kept if r["row"] == row) for row in (LEADER, AN, MATRIX, CFE)}
    assert low == pytest.approx({LEADER: 0.655, AN: 2.951, MATRIX: 3.592, CFE: 4.162}, abs=0.005)
    # the leader's three sites: two near 2.4 to 2.5 and one for about kT
    assert sorted(r["potential"] for r in kept if r["row"] == LEADER) == pytest.approx([0.655, 2.411, 2.487], abs=0.003)
    # that one is the only cell whose third vector had tilted along the chains' neighbours (over 3 A along a and b); the other ten stopped within 1.2 A of upright
    tilt = lambda r: max(abs(r["cell_A_deg"][2] * math.cos(math.radians(r["cell_A_deg"][k]))) for k in (3, 4))  # noqa: E731
    far = [r for r in sheets if tilt(r) > 1.2]
    assert len(far) == 1 and far[0]["row"] == LEADER and far[0]["potential"] == pytest.approx(0.655, abs=0.003) and tilt(far[0]) == pytest.approx(4.0, abs=0.1)
    belt4_low = {LEADER: 1.395, AN: 2.148, MATRIX: 4.659, CFE: 3.418}   # the lowest site of each cell's scan, both hands
    scans = json.loads((PRED / "kink_site_scan_belt4.json").read_text())["rows"]
    for row, e in belt4_low.items():
        sites = [v["formation"] for key in (f"{row}@0", f"{row}@0~") for v in scans[key]["sites"].values() if v.get("converged", True) and v.get("kink_kept", True)]
        assert min(sites) == pytest.approx(e, abs=0.002)
    # the CFE's cut of the matrix's sheet (1.24 in belt4) is not there in the potential
    assert belt4_low[MATRIX] - belt4_low[CFE] == pytest.approx(+1.24, abs=0.01) and low[MATRIX] - low[CFE] == pytest.approx(-0.57, abs=0.01)
    stored = json.loads((PRED / "four_axes_potential_margins_sheets.json").read_text())
    assert stored["sheet_shifts_kcal_per_kink_pair"] == pytest.approx({k: low[k] - belt4_low[k] for k in low}, abs=0.002)
    at = stored["shifts"]["+0.00"]["rows"]
    khz = {k: at[k][SUM]["f_relax_hz"] / 1e3 for k in at}
    assert khz == pytest.approx({MATRIX: 697, CFE: 364, AN: 1590, LEADER: 82945}, rel=0.01)   # the leader's is the formula outside its range (a sheet for kT)
    assert all(MOVE not in r for r in at.values())   # the other rate model does not see the sites
    before = _stored()["shifts"]["+0.00"]["rows"]
    assert before[MATRIX][SUM]["f_relax_hz"] / 1e3 == pytest.approx(33, abs=1) and before[LEADER][SUM]["f_relax_hz"] / 1e3 == pytest.approx(6877, rel=0.01)
    assert at[LEADER][SUM]["tan_delta_34kHz"] < at[AN][SUM]["tan_delta_34kHz"] < at[MATRIX][SUM]["tan_delta_34kHz"] < at[CFE][SUM]["tan_delta_34kHz"]


def test_with_the_shear_free_the_films_sheets_fall_below_all_trans_and_the_leaders_stay_above():
    """Every kink site of the four cells with the cell free to shear along its chains (examples/kink_sheet_shear.py, belt4)."""
    scan = json.loads((PRED / "kink_site_scan_belt4_shear.json").read_text())
    records = json.loads((PRED / "kink_sheet_shear_belt4_sites.json").read_text())["rows"]
    assert scan["sites_left_upright"] == [] and len(records) == 108

    def sites(row):
        return [(v["formation"], v["upright_formation"], v["units"], v["shear_monomers"], f"{row}{h} site {o}") for h in ("@0", "@0~")
                for o, v in scan["rows"][row + h]["sites"].items() if v.get("converged", True) and v.get("kink_kept", True)]

    for row in (CFE, AN, LEADER, MATRIX):  # the scan's energy of a site is the lowest sheet of its record that kept its kink
        for formation, _, _, _, key in sites(row):
            rec = records[f"{key} pvdf-dft-valence-belt4"]
            assert formation == pytest.approx(min(s["formation"] for s in rec["starts"] if s.get("kink_kept") and "failed" not in s)) == pytest.approx(rec["lowest"])
    low = {row: min(sites(row)) for row in (CFE, AN, LEADER, MATRIX)}
    assert [len(sites(row)) for row in (CFE, AN, LEADER, MATRIX)] == [28, 26, 26, 28]
    # upright, the leader had the cheapest sheet of the four; with the shear free it has the second dearest
    assert {row: min(s[1] for s in sites(row)) for row in low} == pytest.approx({CFE: 3.42, AN: 2.15, LEADER: 1.40, MATRIX: 4.66}, abs=0.006)
    assert {row: low[row][0] for row in low} == pytest.approx({CFE: -2.27, AN: -2.08, LEADER: 0.65, MATRIX: 0.94}, abs=0.006)
    assert low[CFE][2] == "VC" and low[CFE][3] == pytest.approx([0.0, -1.78], abs=0.03)  # at the chlorine, the cell sheared by nearly two monomers along b
    assert {row: sum(1 for s in sites(row) if s[0] < 0.0) for row in low} == {CFE: 9, AN: 6, LEADER: 0, MATRIX: 0}
    # how many sites end in a cell sheared by half a monomer or more, and what the sites fall by on average: the leader's sheets mostly stay upright
    inclined = {row: sum(1 for s in sites(row) if max(abs(s[3][0]), abs(s[3][1])) > 0.5) for row in low}
    assert inclined == {CFE: 28, AN: 24, LEADER: 6, MATRIX: 24}
    fall = {row: sum(s[1] - s[0] for s in sites(row)) / len(sites(row)) for row in low}
    assert fall[CFE] == pytest.approx(4.34, abs=0.02) and fall[LEADER] == pytest.approx(0.78, abs=0.02) and fall[MATRIX] == pytest.approx(3.19, abs=0.02)


def test_on_the_shear_free_sheets_the_leader_is_with_the_matrix_on_frequency_and_loss():
    """The four axes at the potential's margins with the rate by those sites: stored, belt4's sheets, one model of the rate in two readings."""
    stored = json.loads((PRED / "four_axes_potential_margins_shear.json").read_text())
    assert stored["sites_from"] == "kink_site_scan_belt4_shear.json"
    occupancy = "by_the_occupancy_of_kink_sites"
    for entry in stored["shifts"].values():
        rows, film = entry["rows"], entry[FILM]
        assert rows[CFE][SUM]["f_relax_hz"] == pytest.approx(film["f_relax_hz"]) == pytest.approx(rows[CFE][occupancy]["f_relax_hz"])  # the anchor
        assert rows[CFE]["kink_sites"]["sites_below_all_trans"] == 9 and rows[LEADER]["kink_sites"]["sites_below_all_trans"] == 0
        for row in (LEADER, MATRIX):  # two decades below the film by the sum, one by the occupancy: at or below the 34 kHz drive either way
            assert rows[row]["kink_sites"]["sum_against_the_film_cell"] == pytest.approx(0.012, abs=0.002)
            assert rows[row]["kink_sites"]["occupancy_against_the_film_cell"] == pytest.approx(0.09, abs=0.012)
            assert rows[row][SUM]["f_relax_hz"] < 6e3 and 28e3 < rows[row][occupancy]["f_relax_hz"] < 36e3
            for model in (SUM, occupancy):
                assert rows[row][model]["tan_delta_34kHz"] > 0.5 and rows[row][model]["strain_ac_at_the_heat_budget_34kHz"] < 0.4 * film["strain_ac_at_the_heat_budget_34kHz"]
        for model in (SUM, occupancy):  # the nitrile alone stays with the film
            assert 0.8 < rows[AN][model]["f_relax_hz"] / film["f_relax_hz"] < 1.0 and rows[AN][model]["tan_delta_34kHz"] < 0.07
    zero = stored["shifts"]["+0.00"]["rows"]
    assert zero[LEADER][SUM]["f_relax_hz"] == pytest.approx(4.1e3, rel=0.03) and zero[LEADER][occupancy]["f_relax_hz"] == pytest.approx(30.4e3, rel=0.03)
    assert zero[AN][SUM]["f_relax_hz"] == pytest.approx(306e3, rel=0.02) and zero[AN][occupancy]["f_relax_hz"] == pytest.approx(335e3, rel=0.02)


def test_three_evaluators_put_the_leaders_cheapest_sheet_three_kcal_above_the_film_cells():
    """belt4's cheapest sheet of the film cell and of the leader, relaxed again under chain1-cfe and in the potential (tight stop)."""
    belt4 = json.loads((PRED / "kink_sheet_shear_belt4_sites.json").read_text())["rows"]
    chain1 = json.loads((PRED / "kink_sheet_shear_chain1cfe_sites.json").read_text())["rows"]
    potential = json.loads((PRED / "sheets14_potential_sheared.json").read_text())["rows"]
    film, leader = f"{CFE}@0 site 12", f"{LEADER}@0 site 10"
    b = {k: belt4[f"{k} pvdf-dft-valence-belt4"]["lowest"] for k in (film, leader)}
    c = {k: chain1[f"{k} pvdf-dft-chain1-cfe"]["lowest"] for k in (film, leader)}
    p = {film: potential[f"{CFE}@12"], leader: potential[f"{LEADER}@10"]}
    assert b == pytest.approx({film: -2.27, leader: 0.65}, abs=0.006) and c == pytest.approx({film: -1.65, leader: 1.66}, abs=0.006)
    assert p[film]["potential"] == pytest.approx(-3.60, abs=0.006) and p[leader]["potential"] == pytest.approx(0.05, abs=0.006)
    for rec in p.values():  # the kink is still there, the reference is all-trans, and both were relaxed to 0.005 eV/A or tighter
        assert rec["kink_kept"] and rec["all_trans_is_all_trans"] and all(float(v.split(" to ")[1].split(" eV")[0]) <= 0.005 for v in rec["relaxed"].values())
    gaps = [b[leader] - b[film], c[leader] - c[film], p[leader]["potential"] - p[film]["potential"]]
    assert gaps == pytest.approx([2.92, 3.31, 3.65], abs=0.01)
    # the nitrile alone, level with the film cell in belt4 (-2.08), is above all-trans in the potential: belt4 has the nitrile wrong (3.12)
    nitrile = potential[f"{AN}@10"]
    assert nitrile["potential"] == pytest.approx(0.36, abs=0.006) and nitrile["polyfind"] == pytest.approx(-2.08, abs=0.006) and nitrile["kink_kept"]
    # the film cell's sheets are below all-trans under chain1-cfe too: eight of the nine sites run
    sites = [r["lowest"] for k, r in chain1.items() if k.startswith(CFE)]
    assert len(sites) == 9 and sum(1 for v in sites if v < 0.0) == 8 and min(sites) == pytest.approx(-4.50, abs=0.006)


def test_at_a_tight_stop_the_potentials_margins_move_and_the_film_cells_sheet_is_below_its_helix():
    """The four cells' all-trans and lowest helix cells in the potential at 0.005 eV/A or tighter (3.14 stopped the helix cells at 0.02)."""
    tight = json.loads((PRED / "rows14_potential_tight.json").read_text())["rows"]
    loose = json.loads((PRED / "rows14_potential.json").read_text())["rows"]
    sheets = json.loads((PRED / "sheets14_potential_sheared.json").read_text())["rows"]
    assert {r: tight[r]["margin"] for r in tight} == pytest.approx({MATRIX: 0.016, CFE: 0.129, AN: -0.109, LEADER: -0.070}, abs=0.001)
    assert {r: loose[r]["margin"] for r in tight} == pytest.approx({MATRIX: -0.024, CFE: 0.078, AN: -0.118, LEADER: -0.055}, abs=0.001)
    for row in tight:  # every cell fell or stayed, by up to 0.05 kcal/mol per monomer: the size of the margins themselves
        for conf, cell in tight[row]["cells"].items():
            fall = (loose[row]["cells"][conf]["eV_per_monomer"] - cell["eV_per_monomer"]) * 23.0605
            assert 0.0 <= fall < 0.052
    assert tight[MATRIX]["margin"] > 0.0 > loose[MATRIX]["margin"]  # the matrix changes sign: its helix is the lower by a little
    # per chain of fourteen against all-trans, the lowest helix and the sheet at belt4's cheapest site: of these four sheets only the film cell's is below all-trans,
    # and it is below its helix (the leader's second site, in the test below, is below all-trans too)
    helix = {r: -14 * tight[r]["margin"] for r in tight}
    sheet = {CFE: sheets[f"{CFE}@12"]["potential"], MATRIX: sheets[f"{MATRIX}@24m"]["potential"], AN: sheets[f"{AN}@10"]["potential"], LEADER: sheets[f"{LEADER}@10"]["potential"]}
    assert helix == pytest.approx({MATRIX: -0.23, CFE: -1.81, AN: 1.53, LEADER: 0.98}, abs=0.01)
    assert sheet == pytest.approx({MATRIX: 3.13, CFE: -3.60, AN: 0.36, LEADER: 0.05}, abs=0.006)
    assert sheet[CFE] < helix[CFE] - 1.7 and all(rec["kink_kept"] and rec["all_trans_is_all_trans"] for rec in sheets.values())
    # so at belt4's cheapest sites the potential has the leader's sheet three below the matrix's and three and a half above the film cell's: belt4 had it level with the matrix's
    assert sheet[MATRIX] - sheet[LEADER] == pytest.approx(3.08, abs=0.01) and sheet[LEADER] - sheet[CFE] == pytest.approx(3.65, abs=0.01)
    stored = json.loads((PRED / "four_axes_potential_tight_margins_shear.json").read_text())
    film = stored["shifts"]["+0.00"]["rows"][CFE][SUM]  # at its tight margin the film cell's row is on the measured film's calibrated channel
    assert film["tan_delta_34kHz"] == pytest.approx(0.061, abs=0.001) and film["strain_ac_at_the_heat_budget_34kHz"] == pytest.approx(0.00071, abs=0.00001)


def test_in_the_potential_the_leader_has_a_sheet_below_all_trans_at_a_site_belt4_ranks_low():
    """The leader's second site in the potential: started from belt4's upright sheet it slides (the chain axis tilts) and rests below all-trans."""
    stored = json.loads((PRED / "sheets14_potential_sheared.json").read_text())
    sheets, control = stored["rows"], stored["all_trans_controls"][LEADER]
    first, second = sheets[f"{LEADER}@10"], sheets[f"{LEADER}@18m"]
    # belt4 has this site 1.3 above the leader's cheapest; the potential has it 1.4 below the site belt4 ranks cheapest
    assert second["potential"] == pytest.approx(-1.34, abs=0.006) and second["polyfind"] == pytest.approx(1.99, abs=0.006)
    assert second["potential"] < first["potential"] - 1.3 and first["polyfind"] < second["polyfind"] - 1.3
    assert second["kink_kept"] and second["all_trans_is_all_trans"] and "0.004 eV/A" in second["relaxed"]["sheet"] and "NOT" not in second["relaxed"]["sheet"]
    # the reference does not follow it: the leader's all-trans cell started at the sheet's tilt rests, still all-trans, within 0.05 of the upright one
    assert control["is_all_trans"] and abs(control["potential"]) < 0.05 and "NOT" not in control["relaxed"]["this"]
    assert control["against"] == second["all_trans_frame"] and control["cell_A_deg"][3] > control["against_cell_A_deg"][3] + 1.5  # it stayed tilted
    # so the leader too has a kinked state below its all-trans cell and its helix, and the film cell's sheet is 2.3 below it, not 3.65
    tight = json.loads((PRED / "rows14_potential_tight.json").read_text())["rows"]
    assert second["potential"] < 0.0 < -14 * tight[LEADER]["margin"]
    assert second["potential"] - sheets[f"{CFE}@12"]["potential"] == pytest.approx(2.26, abs=0.01)
