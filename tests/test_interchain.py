"""The energy between rigid chains as the cell is strained and sheared: polyfind against DFT, stored results.

Reads the stored files only (`examples/interchain_frames.py` writes them; the DFT is
`examples/mlip_frames.py label`).  Everything inside a chain is the same in every frame of a
packing, so the numbers here are energy between chains and nothing else.
"""

import json
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
PRED = ROOT / "docs" / "predictions"
EV_TO_KCAL = 23.060548


def _table():
    return json.loads((PRED / "interchain_belt4.json").read_text())["packings"]


def _frames():
    """{frame: (atoms, family, DFT energy per monomer in eV)} of the stored labelled frames, without ase."""
    lines = (PRED / "mlip" / "interchain_frames_pbe_d3bj.extxyz").read_text().splitlines()
    out, i = {}, 0
    while i < len(lines):
        n = int(lines[i])
        head = lines[i + 1]
        value = lambda key: head.split(f"{key}=")[1].split()[0]  # noqa: E731
        out[value("frame")] = (n, value("family"), float(value("REF_energy")) / int(value("n_monomers")))
        i += n + 2
    return out


def _row(block, scale, shear=(0.0, 0.0)):
    return next(r for r in block["frames"] if r["scale"] == list(scale) and r["shear_monomers"] == list(shear))


def test_the_stored_frames_are_the_cells_of_the_table():
    frames, t = _frames(), _table()
    assert {f[1] for f in frames.values()} == {"interchain"}
    rows = {r["frame"]: (pack, r) for pack, b in t.items() for r in b["frames"]}
    assert set(rows) == set(frames)
    assert sum(pack == "beta" for pack, _ in rows.values()) == 27
    assert sum(pack == "alpha" for pack, _ in rows.values()) == 27
    assert sum(pack == "hh" for pack, _ in rows.values()) == 37
    assert sum(pack == "candidate" for pack, _ in rows.values()) == 27 and sum(pack == "cand6" for pack, _ in rows.values()) == 4
    atoms = {"beta": 12, "candidate": 12, "alpha": 24, "hh": 24, "cand6": 72}
    for name, (pack, r) in rows.items():
        assert frames[name][0] == atoms[pack] == 6 * t[pack]["monomers_per_cell"]
        zero = frames[f"inter_{pack}_00"][2]
        assert r["dft"] == pytest.approx((frames[name][2] - zero) * EV_TO_KCAL, abs=2e-6)
        # polyfind's energy is its electrostatics and its repulsion and dispersion, and nothing else moves
        assert r["polyfind"] == pytest.approx(r["polyfind_electrostatics"] + r["polyfind_repulsion_dispersion"], abs=1e-9)


def test_plain_pvdf_is_right_stretched_and_too_stiff_squeezed():
    """The beta cell under belt4: within a tenth of DFT wherever a lateral edge is 2% or more
    longer, 1.1 to 1.45 times DFT where one is 4% or more shorter."""
    beta = _table()["beta"]
    for scale in ((1.02, 1.0), (1.04, 1.0), (1.06, 1.0), (1.0, 1.04), (1.0, 1.06), (1.02, 1.02), (1.04, 1.04)):
        r = _row(beta, scale)
        assert 0.90 < r["polyfind"] / r["dft"] < 1.05, scale
    for scale in ((0.94, 1.0), (0.96, 1.0), (1.0, 0.94), (1.0, 0.96), (0.96, 0.96)):
        r = _row(beta, scale)
        assert 1.10 < r["polyfind"] / r["dft"] < 1.45, scale


def test_the_registry_cost_is_overstated_more_the_closer_the_chains():
    """Half a monomer along ``a`` in the beta cell, each cell against its own unsheared state: 2.4
    times DFT in a cell 4% wider, 5.4 times in one 4% narrower.  The wall is too steep close in,
    which no 12-6 coefficients repair."""
    beta = _table()["beta"]
    wide, narrow = _row(beta, (1.04, 1.04), (0.5, 0.0)), _row(beta, (0.96, 0.96), (0.5, 0.0))
    # each against its own unsheared cell
    cost = lambda r, s: (r["polyfind"] - _row(beta, s)["polyfind"]) / (r["dft"] - _row(beta, s)["dft"])  # noqa: E731
    assert 1.8 < cost(wide, (1.04, 1.04)) < 2.6
    assert 4.5 < cost(narrow, (0.96, 0.96)) < 6.0
    assert all(_row(beta, (s, s), sh)["polyfind"] > _row(beta, (s, s), sh)["dft"] for s in (0.96, 1.04) for sh in ((0.25, 0.0), (0.5, 0.0), (0.0, 0.5)))


def test_the_reversed_chain_is_missed_the_other_way():
    """Every second unit reversed: fluorine against fluorine and hydrogen against hydrogen between
    chains.  DFT's cell wants to be wider along ``a`` than polyfind's (its energy falls on
    stretching, polyfind's rises) and resists squeezing four times as hard; and half a monomer
    of shear one way along ``a`` costs three times what it costs the other way, which polyfind
    halves."""
    hh = _table()["hh"]
    assert _row(hh, (1.04, 1.0))["dft"] < -0.15 and _row(hh, (1.04, 1.0))["polyfind"] > 0.04
    assert 3.5 < _row(hh, (0.94, 1.0))["dft"] / _row(hh, (0.94, 1.0))["polyfind"] < 4.6
    there, back = _row(hh, (1.0, 1.0), (0.5, 0.0)), _row(hh, (1.0, 1.0), (1.5, 0.0))  # 1.5 is minus a half: the chain repeats every two monomers
    assert 3.0 < back["dft"] / there["dft"] < 3.5 and 1.9 < back["polyfind"] / there["polyfind"] < 2.2
    err = np.array([r["polyfind"] - r["dft"] for r in hh["frames"][1:]])
    dft = np.array([r["dft"] for r in hh["frames"][1:]])
    assert 0.33 < np.sqrt((err ** 2).mean()) < 0.40 and 0.80 < np.sqrt((dft ** 2).mean()) < 0.90


def test_the_candidates_own_cell_gains_as_much_rigid_as_relaxed():
    """`vdf5-rvdf1-6` all-trans, its own chains rigid, sheared along the line DFT's relaxation took
    it: DFT gains 0.34 kcal/mol per monomer, where its full relaxation gains 0.35 and belt4 loses
    1.58.  The gain is between chains."""
    t = _table()["cand6"]
    assert t["monomers_per_cell"] == 12
    end = _row(t, (1.0, 1.0), (0.30, -0.90))
    assert end["dft"] == pytest.approx(-0.338, abs=0.002)
    assert end["polyfind"] == pytest.approx(1.576, abs=0.002)
    half = _row(t, (1.0, 1.0), (0.15, -0.45))   # the rigid cell's own minimum is part way along the line
    assert half["dft"] == pytest.approx(-0.412, abs=0.002) and half["polyfind"] == pytest.approx(0.695, abs=0.002)
    quarter = _row(t, (1.0, 1.0), (0.075, -0.225))
    assert quarter["dft"] == pytest.approx(-0.359, abs=0.002) and quarter["polyfind"] == pytest.approx(0.103, abs=0.002)
    assert half["dft"] < quarter["dft"] < end["dft"] < 0.0
    relaxed = json.loads((PRED / "shear_registry_belt4.json").read_text())["candidate"]["dft"]["sheared_below_upright"]
    assert relaxed == pytest.approx(-0.346, abs=0.002) and abs(end["dft"] - relaxed) < 0.02


def test_the_predictions_registered_before_the_candidates_dft():
    """Four variants, each closer than the one before; the best is still 0.30 short."""
    p = json.loads((PRED / "interchain_cand6_predictions.json").read_text())
    dft = _row(_table()["cand6"], (1.0, 1.0), (0.30, -0.90))["dft"]
    names = ["belt4 as it is", "buffered wall, present charges", "buffered wall, bond site", "buffered wall, bond site and backbone pull"]
    assert list(p["predictions"]) == names
    at_end = [p["predictions"][n]["cells"]["inter_cand6_02"]["predicted"] for n in names]
    assert at_end == pytest.approx([1.576, 0.153, 0.029, -0.039], abs=0.002)
    miss = [v - dft for v in at_end]
    assert all(a > b > 0 for a, b in zip(miss, miss[1:])) and 0.28 < miss[-1] < 0.32
    # half way the same order holds, the best 0.18 short; and every wall variant has the cell lower half way than at the end, as DFT does
    half_dft = _row(_table()["cand6"], (1.0, 1.0), (0.15, -0.45))["dft"]
    at_half = [p["predictions"][n]["cells"]["inter_cand6_01"]["predicted"] for n in names]
    miss_half = [v - half_dft for v in at_half]
    assert all(a > b > 0 for a, b in zip(miss_half, miss_half[1:])) and 0.16 < miss_half[-1] < 0.20
    assert half_dft < dft and all(h < e for h, e in zip(at_half[1:], at_end[1:])) and at_half[0] < at_end[0]
    quarter_dft = _row(_table()["cand6"], (1.0, 1.0), (0.075, -0.225))["dft"]
    at_quarter = [p["predictions"][n]["cells"]["inter_cand6_03"]["predicted"] for n in names]
    assert 0.09 < at_quarter[-1] - quarter_dft < 0.12 and at_quarter[0] - quarter_dft > 0.45
    assert p["predictions"][names[0]]["cells"]["inter_cand6_02"]["predicted"] == pytest.approx(_row(_table()["cand6"], (1.0, 1.0), (0.30, -0.90))["polyfind"], abs=1e-6)
    # the charges of the last variant are the potential's, not fitted to any energy
    v = p["predictions"][names[-1]]["variables"]
    assert (v["s_CH"], v["s_CF"], v["d_CC"], v["M_H"]) == (0.300, 0.701, -0.031, 0.318)


def _trials():
    return json.loads((PRED / "interchain_wall_trials.json").read_text())["trials"]


def _rms(trial, group):
    return trial["after"][group]["rms_error"]


def test_the_stored_trials_are_fits_of_all_the_cells():
    """Every trial of `examples/interchain_wall.py` that docs/CRYSTAL_STAGE.md 5.26 quotes, on the
    150 rigid cells (143 of them rows: each group against its own first cell)."""
    t = _trials()
    assert len(t) == 11 and all(tr["after"]["all"]["cells"] == 143 for tr in t.values())
    groups = {"beta strain": 26, "beta shear": 13, "candidate strain": 26, "candidate shear": 13, "alpha strain": 26, "hh strain": 36, "cand6 strain": 3}
    for tr in t.values():
        assert {k: v["cells"] for k, v in tr["after"].items() if k != "all"} == groups
    table, frames = _table(), {}
    for pack, b in table.items():
        frames.update({r["frame"]: r["dft"] for r in b["frames"]})
    # the DFT the trials were fitted to is the stored table's, row for row
    first = next(iter(t.values()))
    for r in first["rows"]:
        if r["frame"] in frames:
            assert r["dft"] == pytest.approx(frames[r["frame"]], abs=2e-6)
    # and the first trial is belt4 itself
    belt4 = t["belt4 as it is"]
    for r in belt4["rows"]:
        if r["frame"].startswith("inter_cand6"):
            assert r["model"] == pytest.approx(next(x["polyfind"] for x in table["cand6"]["frames"] if x["frame"] == r["frame"]), abs=1e-6)


def test_what_the_trials_say():
    t = _trials()
    plain = ("beta strain", "beta shear", "candidate strain", "candidate shear", "alpha strain")
    wall = t["buffered wall, polyfind's charges"]
    # 12-6 cannot be tuned to the plain cells; the buffered wall can
    assert max(_rms(t["12-6, diameters and belt wall refitted"], g) for g in plain) > 0.3 and max(_rms(wall, g) for g in plain) < 0.09
    # the reversed chain is not a matter of the wall: its own F-F and H-H terms do not help, and end on bounds
    like = t["buffered wall, polyfind's charges, F-F and H-H pair terms"]
    assert 0.16 < _rms(wall, "hh strain") < 0.18 and _rms(like, "hh strain") > _rms(wall, "hh strain") - 0.01 and len(like["at_bound"]) == 3
    # the two charge terms, at the potential's strengths, nothing of them fitted to an energy
    site = t["buffered wall, bond site, charges from the potential"]
    both = t["buffered wall, bond site and backbone pull, charges from the potential"]
    assert _rms(site, "hh strain") == pytest.approx(0.142, abs=0.003) and _rms(both, "hh strain") == pytest.approx(0.091, abs=0.003)
    assert all(k not in both["free"] for k in ("s_CH", "s_CF", "d_CC", "M_H"))
    assert (both["variables"]["s_CH"], both["variables"]["s_CF"], both["variables"]["d_CC"], both["variables"]["M_H"]) == (0.300, 0.701, -0.031, 0.318)
    # with an F-F term the candidate's own line comes within 0.09, and weighting it within 0.04 at a cost to the reversed chain
    ff = t["the same with the F-F pair's own diameter and depth"]
    heavy = t["the same, the candidate's own cells counted four times"]
    assert _rms(ff, "cand6 strain") < 0.09 and _rms(ff, "hh strain") < 0.11 and max(_rms(ff, g) for g in plain) < 0.03
    assert _rms(heavy, "cand6 strain") < 0.04 and 0.12 < _rms(heavy, "hh strain") < 0.14 and max(_rms(heavy, g) for g in plain) < 0.06
    # the charge terms are what let one wall hold both: the same weighting without them
    without = t["polyfind's charges, F-F term, the candidate's own cells counted four times"]
    assert _rms(without, "hh strain") > 0.24
    # and energies cannot set the charges: fitted to them, the atoms' own increments go to zero
    free = t["charges fitted to the energies too (they leave the potential's values)"]
    assert free["variables"]["s_CH"] < 1e-6 and free["variables"]["s_CF"] < 1e-6 and free["variables"]["M_H"] > 0.7
