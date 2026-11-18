"""What it costs to move neighbouring chains along each other: polyfind against DFT, stored results.

Reads the stored files only (`examples/shear_registry.py` and `examples/shear_polish.py` write
them; the DFT is `examples/mlip_frames.py label`).  The shear itself is held to known answers in
`tests/test_mixedcell.py`.
"""

import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
PRED = ROOT / "docs" / "predictions"
EV_TO_KCAL = 23.060548


def _table():
    return json.loads((PRED / "shear_registry_belt4.json").read_text())


def _frames():
    """{frame: (atoms, family, DFT energy per monomer in eV)} of the stored labelled frames, without ase."""
    lines = (PRED / "mlip" / "shear_registry_frames_pbe_d3bj.extxyz").read_text().splitlines()
    out, i = {}, 0
    while i < len(lines):
        n = int(lines[i])
        head = lines[i + 1]
        value = lambda key: head.split(f"{key}=")[1].split()[0]  # noqa: E731
        out[value("frame")] = (n, value("family"), float(value("REF_energy")) / int(value("n_monomers")))
        i += n + 2
    return out


def test_the_stored_frames_are_the_cells_of_the_table():
    frames, t = _frames(), _table()
    assert len(frames) == 35
    families = [f[1] for f in frames.values()]
    assert families.count("shear") == 28 and families.count("e5") == 7
    assert {f[0] for f in frames.values() if f[1] == "shear"} == {12}
    for pack, block in t["packings"].items():
        rows = block["shears"]
        assert len(rows) == 14 and rows[0]["shear_monomers"] == [0.0, 0.0] and rows[0]["polyfind"] == 0.0 and rows[0]["dft"] == 0.0
        zero = frames[rows[0]["frame"]][2]
        for row in rows:  # the table's DFT column is the frames' energies, nothing else
            assert (frames[row["frame"]][2] - zero) * EV_TO_KCAL == pytest.approx(row["dft"], abs=1e-6)
    ref = frames["e5_beta_2a_reference"][2]
    for name, row in t["e5_like_for_like"].items():
        assert (frames[name][2] - ref) * EV_TO_KCAL == pytest.approx(row["dft"], abs=1e-6)


def test_polyfind_charges_about_three_times_dft_for_moving_plain_chains_along_each_other():
    """Every sheared cell of plain PVDF, in both packings: the sign is DFT's and the size is not."""
    t = _table()
    for pack, lo, hi in (("beta", 2.6, 3.6), ("candidate", 2.4, 3.1)):
        rows = t["packings"][pack]["shears"][1:]
        assert all(r["polyfind"] > r["dft"] > 0.0 for r in rows)
        assert lo < min(r["ratio"] for r in rows) and max(r["ratio"] for r in rows) < hi
    beta = {tuple(r["shear_monomers"]): r for r in t["packings"]["beta"]["shears"]}
    assert beta[(0.5, 0.0)]["polyfind"] == pytest.approx(2.343, abs=0.002) and beta[(0.5, 0.0)]["dft"] == pytest.approx(0.692, abs=0.002)
    assert beta[(0.0, 0.5)]["polyfind"] == pytest.approx(0.732, abs=0.002) and beta[(0.0, 0.5)]["dft"] == pytest.approx(0.240, abs=0.002)
    # the ratio grows with the displacement: the wall is too steep, not only too high
    along_a = [beta[(k / 8, 0.0)]["ratio"] for k in (1, 2, 3, 4)]
    assert along_a == sorted(along_a) and along_a[-1] - along_a[0] > 0.5
    assert t["cells"] == 32 and t["ratio_range"][0] == pytest.approx(2.40, abs=0.01) and t["ratio_range"][1] == pytest.approx(3.53, abs=0.01)
    assert t["ratio_median"] == pytest.approx(2.87, abs=0.02)


def test_the_registered_e5_predictions_miss_their_acceptance_by_size_alone():
    """docs/E1_E2_PREDICTIONS.md, E5: sign, order along the slide, one chain of four above one
    sublattice per moved monomer, magnitudes within a factor 1.5.  All but the last hold."""
    e5 = _table()["e5_like_for_like"]
    assert len(e5) == 6 and all(r["polyfind_registered"] > 0.0 and r["dft"] > 0.0 for r in e5.values())
    a = [e5[f"e5_beta_a_slide_{s}"] for s in ("c8", "c4", "3c8", "c2")]
    assert [r["dft"] for r in a] == sorted(r["dft"] for r in a) and [r["polyfind_registered"] for r in a] == sorted(r["polyfind_registered"] for r in a)
    assert [round(r["dft"], 3) for r in a] == [0.066, 0.230, 0.396, 0.467]
    assert [round(r["polyfind_registered"], 3) for r in a] == [0.172, 0.656, 1.251, 1.533]
    half_a, half_c = e5["e5_beta_a_slide_c2"], e5["e5_beta_2a_slide_c2"]
    # per moved monomer: one sublattice is half the cell's monomers, one chain of four a quarter
    assert 4 * half_c["dft"] / (2 * half_a["dft"]) == pytest.approx(2.01, abs=0.02)
    assert 4 * half_c["polyfind_registered"] / (2 * half_a["polyfind_registered"]) == pytest.approx(2.07, abs=0.02)
    assert all(2.6 < r["ratio"] < 3.4 for r in e5.values())
    assert min(r["ratio"] for r in e5.values()) > 1.5  # the acceptance's factor


def test_with_the_shear_free_polyfind_keeps_the_candidates_cell_nearly_upright_where_dft_does_not():
    c = _table()["candidate"]
    own = c["polyfind_shear_free_from_upright"]
    assert own["gain"] == pytest.approx(-0.029, abs=0.001)
    assert own["shear_monomers"][0] == pytest.approx(0.026, abs=0.001) and own["shear_monomers"][1] == pytest.approx(-0.078, abs=0.001)
    dft = c["dft_relaxed_cell"]
    assert dft["shear_monomers"][0] == pytest.approx(0.300, abs=0.002) and dft["shear_monomers"][1] == pytest.approx(-0.901, abs=0.002)
    # the same direction, one of a to minus three of b, a tenth of the way
    assert own["shear_A"][1] / own["shear_A"][0] == pytest.approx(-3.0, abs=0.01) and dft["shear_A"][1] / dft["shear_A"][0] == pytest.approx(-3.0, abs=0.01)
    assert dft["fit_rms_A"] < 0.02 and dft["fit_worst_A"] < 0.04  # DFT's cell is in the mixed cell's variables to the atom
    assert max(abs(x - y) for x, y in zip(dft["setting_angles_deg"], dft["polyfind_setting_angles_deg"])) < 0.01
    at = c["polyfind_at_dft_cell"]
    assert at["shear_held_rest_relaxed"] == pytest.approx(0.426, abs=0.002) and at["as_fitted"] > 1.0
    assert c["dft"]["sheared_below_upright"] == pytest.approx(-0.346, abs=0.001)
    # a start at DFT's shear with polyfind's own chains runs back to polyfind's minimum
    back = c["polyfind_from_dft_shear_with_its_own_chains"]
    assert back["relaxed"] == pytest.approx(own["gain"], abs=1e-6) and back["rigid_at_start"] == pytest.approx(1.576, abs=0.002)
    # and the size of the disagreement is the plain chains': 1.05 rigid there, 0.77 on the candidate relaxed
    plain = {tuple(r["shear_monomers"]): r for r in _table()["packings"]["candidate"]["shears"]}[(0.303, -0.909)]
    assert plain["polyfind"] - plain["dft"] == pytest.approx(1.05, abs=0.01)
    assert at["shear_held_rest_relaxed"] - c["dft"]["sheared_below_upright"] == pytest.approx(0.77, abs=0.01)


def test_the_excess_is_polyfinds_repulsion_and_dispersion_and_none_of_it_electrostatic():
    """Polyfind's cost of a shear split into its terms: the charges and the induced dipoles favour
    the sheared cell slightly, in every cell; the pair terms alone exceed DFT's whole cost."""
    for block in _table()["packings"].values():
        for r in block["shears"][1:]:
            es, rd = r["polyfind_electrostatics"], r["polyfind_repulsion_dispersion"]
            assert es + rd == pytest.approx(r["polyfind"], abs=1e-9)
            assert es < 0.0 < r["dft"] < rd and abs(es) < 0.13 * rd


def test_the_cost_is_set_by_the_variables_that_place_hydrogens_contacts():
    """One fitted variable stepped at a time: four move the cost of half a monomer of shear by
    more than 3%, all of them diameters or the belt wall, and each pushes the cell apart as well;
    the well depths, the carbon diameter and the charges do not move it."""
    h = _table()["handles"]
    steps, fit = h["steps"], h["fitted"]
    assert fit["half_a"] == pytest.approx(2.343, abs=0.002) and fit["half_b"] == pytest.approx(0.732, abs=0.002)
    moved = {k for k, r in steps.items() if abs(r["half_a_change"]) > 0.03}
    assert moved == {"x_H", "x_HF", "log10_D_X", "x_F"}
    assert steps["x_H"]["half_a_change"] == pytest.approx(0.135, abs=0.003) and steps["x_H"]["half_b_change"] == pytest.approx(0.143, abs=0.003)
    assert all(abs(steps[k]["half_a_change"]) < 0.011 and abs(steps[k]["half_b_change"]) < 0.011 for k in set(steps) - moved)
    for k in moved:  # the lateral force on the cell, kcal/mol/A: more negative is a cell that wants to be wider
        assert steps[k]["dE_dA"] < fit["dE_dA"] - 0.2 and steps[k]["dE_dB"] < fit["dE_dB"] - 0.15


def test_the_four_handles_move_the_cost_and_the_lattice_together():
    """The linear estimate the write-up draws from the stored steps: four targets (two costs
    down by 2.9, two lateral forces unchanged) on the four variables is an ill-conditioned
    problem, and the hydrogen diameter alone buys the factor with forces that close the cell."""
    import numpy as np

    h = _table()["handles"]
    fit, st = h["fitted"], h["steps"]
    names = ["x_H", "x_HF", "log10_D_X", "x_F"]
    M = np.array([[np.log1p(st[n]["half_a_change"]) / st[n]["step"] for n in names],
                  [np.log1p(st[n]["half_b_change"]) / st[n]["step"] for n in names],
                  [(st[n]["dE_dA"] - fit["dE_dA"]) / st[n]["step"] for n in names],
                  [(st[n]["dE_dB"] - fit["dE_dB"]) / st[n]["step"] for n in names]])
    assert np.linalg.cond(M) > 5000
    exact = np.linalg.solve(M, [-np.log(2.9), -np.log(2.9), 0.0, 0.0])
    assert abs(exact[0]) > 1.0  # a hydrogen diameter several times over: no such potential
    cheapest = np.linalg.lstsq(M[:2], np.array([-np.log(2.9)] * 2), rcond=None)[0]
    assert cheapest[0] == pytest.approx(-0.17, abs=0.01) and abs(cheapest[1]) < 0.03 and abs(cheapest[3]) < 0.03
    closing = M[2:] @ cheapest
    assert closing[0] == pytest.approx(5.1, abs=0.2) and closing[1] == pytest.approx(3.0, abs=0.2)


def _sheared():
    return json.loads((PRED / "sheared_margins_belt4.json").read_text())


def test_polyfinds_own_re_polish_moves_helix_cells_and_leaves_the_all_trans_ones_upright():
    """Every polished cell relaxed again with the shear free (examples/shear_polish.py): each
    starts from its stored upright cell, ends stationary and no higher, and what moves is the
    helix side."""
    rows = _sheared()["rows"]
    polished = json.loads((PRED / "polished_margins_belt4.json").read_text())["rows"]
    assert set(rows) == set(polished) and len(rows) == 35
    assert all(set(v["cells"]) == set(polished[r]["cells"]) for r, v in rows.items())
    cells = [(r, c, x) for r, v in rows.items() for c, x in v["cells"].items()]
    assert len(cells) == 168
    for r, c, x in cells:
        assert x["upright_per_monomer"] == pytest.approx(polished[r]["cells"][c]["rigid_per_monomer"], abs=1e-6)
        assert x["gain_per_monomer"] < 1e-6 and x["gradient_max"] < 1e-4 and not x.get("stalled")
        assert x["sheared_per_monomer"] - x["upright_per_monomer"] == pytest.approx(x["gain_per_monomer"], abs=1e-9)
    tt = [x["gain_per_monomer"] for _, c, x in cells if c == "TTTT"]
    tg = [x["gain_per_monomer"] for _, c, x in cells if c != "TTTT"]
    assert min(tt) == pytest.approx(-0.029, abs=0.001) and sum(g < -0.01 for g in tt) == 4
    assert min(tt) == rows["vdf5-rvdf1-6"]["cells"]["TTTT"]["gain_per_monomer"]  # the reversed-unit candidate's, DFT's -0.35
    assert min(tg) == pytest.approx(-0.294, abs=0.001) and sum(g < -0.05 for g in tg) == 38
    assert sum(tt) / len(tt) == pytest.approx(-0.004, abs=0.001) and sum(tg) / len(tg) == pytest.approx(-0.038, abs=0.001)


def test_the_rigid_margins_move_by_as_much_as_the_window_is_wide():
    rows = _sheared()["rows"]
    for v in rows.values():
        c = v["cells"]
        tgs = [k for k in c if k != "TTTT"]
        assert v["margin_rigid_upright"] == pytest.approx(c["TTTT"]["upright_per_monomer"] - min(c[k]["upright_per_monomer"] for k in tgs), abs=1e-9)
        assert v["margin_rigid_sheared"] == pytest.approx(c["TTTT"]["sheared_per_monomer"] - min(c[k]["sheared_per_monomer"] for k in tgs), abs=1e-9)
    change = {r: v["margin_change"] for r, v in rows.items()}
    assert min(change.values()) == pytest.approx(-0.028, abs=0.001) and max(change.values()) == pytest.approx(0.265, abs=0.001)
    assert min(change, key=change.get) == "vdf5-rvdf1-6" and max(change, key=change.get) == "vdf9-trfe3-an1-rvdf1-14"
    assert sum(abs(v) > 0.06 for v in change.values()) == 9 and sum(abs(v) > 0.12 for v in change.values()) == 4  # the window is 0.12 wide
    assert sum(change.values()) / len(change) == pytest.approx(0.044, abs=0.001)
    assert change["vdf9-trfe4-cfe1-14"] == pytest.approx(0.111, abs=0.001)  # the terpolymer the candidates are set against
    assert sum(v["tg_conformer_upright"] != v["tg_conformer_sheared"] for v in rows.values()) == 16
    # on the one row DFT has, polyfind's change is a twelfth of DFT's: -0.028 against -0.33
    assert change["vdf5-rvdf1-6"] / (-0.276 - 0.070) == pytest.approx(0.082, abs=0.005)
