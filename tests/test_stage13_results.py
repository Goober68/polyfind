"""Stage 13's first run and the checks made on it, stored results (docs/CRYSTAL_STAGE.md 5.28)."""

import json
from pathlib import Path

import numpy as np
import pytest

from polyfind import crystalfit as C

ROOT = Path(__file__).resolve().parents[1]
PRED = ROOT / "docs" / "predictions"


def _fit():
    return json.loads((PRED / "crystal_stage_v13_run1.json").read_text())


def test_the_first_run_is_a_stage13_vector_and_says_it_is_not_finished():
    d = _fit()
    assert tuple(d["variables"]) == C.VARIABLES13 and len(d["x"]) == 21 and "Not a preset" in d["stopped"]["why"]
    x = dict(zip(d["variables"], d["x"]))
    prm = C.parameters(np.array(d["x"]))
    assert prm.bond_sites == (("H", pytest.approx(x["q_M"])),) and prm.backbone_pull == (("F", pytest.approx(x["pull_F"])),)
    assert prm.torsion == pytest.approx((x["V1"], x["V2"], x["V3"])) and prm.wall[0] == "b14-7"
    # the charges stayed where the potential put them
    assert x["q_M"] == pytest.approx(-0.318, abs=0.03) and x["pull_F"] == pytest.approx(0.0155, abs=0.002)
    r = d["residuals"]
    assert d["objective"] == pytest.approx(sum(v["z"] ** 2 for v in r.values()) + d["ridge"], abs=1e-6)
    # the lattice rows within a tolerance and a quarter, one stiffness constant apart
    lattice = ["beta a (polar)", "beta b (lateral)", "beta c (chain)", "alpha a", "alpha b", "alpha c", "beta |P|", "E(alpha)-E(beta)"]
    assert max(abs(r[k]["z"]) for k in lattice) < 1.0 and r["beta C_12"]["z"] > 2.0
    assert r["E(alpha)-E(beta)"]["polyfind"] == pytest.approx(-0.772, abs=0.002)
    # the rows from outside the module
    for pack, limit in (("beta", 0.07), ("candidate", 0.06), ("alpha", 0.06), ("hh", 0.13), ("cand6", 0.08)):
        err = np.array([v["polyfind"] - v["dft"] for k, v in r.items() if k.startswith(f"rigid {pack}_")])
        assert len(err) and np.sqrt((err ** 2).mean()) < limit, pack
    assert all(r[f"potential around {c}"]["polyfind"] < 0.23 for c in ("ht", "hh", "tgtg"))
    assert r["conformers, energy rms"]["polyfind"] == pytest.approx(1.317, abs=0.005)


def test_the_second_run_gives_nothing_back():
    one, two = _fit(), json.loads((PRED / "crystal_stage_v13_run2.json").read_text())
    assert tuple(two["variables"]) == C.VARIABLES13 and "Not a preset" in two["stopped"]["why"]
    assert two["objective"] < one["objective"] - 15.0
    r1, r2 = one["residuals"], two["residuals"]
    assert r2["E(alpha)-E(beta)"]["polyfind"] == pytest.approx(-0.730, abs=0.002)
    rms = lambda r, pack: float(np.sqrt(np.mean([(v["polyfind"] - v["dft"]) ** 2 for k, v in r.items() if k.startswith(f"rigid {pack}_")])))  # noqa: E731
    for pack in ("beta", "candidate", "alpha", "cand6"):
        assert rms(r2, pack) < rms(r1, pack)
    # the reversed chain's cells are where the terms the fit has stop
    assert rms(r2, "hh") == pytest.approx(rms(r1, "hh"), abs=0.015) and rms(r2, "hh") > 0.11
    registry = [k for k in r2 if k.startswith("registry ")]
    assert len(registry) == 8 and all(abs(r2[k]["z"]) < abs(r1[k]["z"]) for k in registry)
    assert all(0.75 < r2[k]["polyfind"] / r2[k]["dft"] < 0.87 for k in registry)
    # what it paid: beta's polar axis, the one cell row past its tolerance
    assert r2["beta a (polar)"]["z"] == pytest.approx(-1.01, abs=0.02)


def test_what_holds_a_packing_together():
    """The cell minus its own chains: belt4 binds too weakly by a quarter to a half, the buffered
    wall takes most of that away unasked, and the difference between packings is right throughout."""
    d = json.loads((PRED / "cohesion_rigid_cells.json").read_text())
    dft = d["dft"]
    assert dft["beta"] == pytest.approx(-7.815, abs=0.005) and dft["alpha"] == pytest.approx(-5.074, abs=0.005) and dft["hh"] == pytest.approx(-3.110, abs=0.005)
    box = lambda name, pack: d["polyfind"][name][pack]["against_a_16_A_lattice"]  # noqa: E731
    for pack in dft:
        assert 0.50 < box("belt4", pack) / dft[pack] < 0.75
        assert 0.70 < box("stage 13, first run", pack) / dft[pack] < 0.90
        assert abs(box("stage 13, first run", pack) - dft[pack]) < abs(box("belt4", pack) - dft[pack]) - 0.5
    for name in d["polyfind"]:
        assert box(name, "alpha") - box(name, "beta") == pytest.approx(dft["alpha"] - dft["beta"], abs=0.12)


def test_the_acceptance_on_the_candidate_has_two_outcomes():
    one = json.loads((PRED / "candidate_acceptance_stage13_run1.json").read_text())
    typed = json.loads((PRED / "candidate_acceptance_stage13_run1_typed.json").read_text())
    # the half the refit was for: the all-trans cell shears and gains what DFT gains, the helix cell does not
    assert one["all_trans_shear_gain"] == pytest.approx(-0.309, abs=0.003) and one["dft"]["all_trans_shear_gain"] == -0.346
    assert all(abs(one["cells"][c]["gain"]) < 0.01 for c in ("TG'TG", "GTG'T"))
    assert abs(typed["all_trans_shear_gain"] - one["all_trans_shear_gain"]) < 0.005
    # the half belt4 had: the levels against DFT's own helix cell, and the order of the two helix placements
    a = one["against_the_helix_dft_relaxed"]
    assert a["helix"] == "TG'TG" and a["dft"] == {"margin_upright": 0.070, "margin_sheared": -0.261, "other_helix_minus_it": 0.197}
    assert a["margin_upright"] == pytest.approx(-0.135, abs=0.003) and a["margin_sheared"] == pytest.approx(-0.444, abs=0.003)
    assert a["other_helix_minus_it"] == pytest.approx(-0.158, abs=0.003)
    # the shear's part of the margin is right; the level is what is off
    assert (a["margin_sheared"] - a["margin_upright"]) == pytest.approx(a["dft"]["margin_sheared"] - a["dft"]["margin_upright"], abs=0.03)
    # typed torsions from the conformer refit make the order of the placements worse
    assert typed["against_the_helix_dft_relaxed"]["other_helix_minus_it"] < a["other_helix_minus_it"] - 0.25


def test_with_the_lone_chains_line_in_the_rows_the_chain_is_dfts_and_the_candidates_level_is_back():
    """docs/CRYSTAL_STAGE.md 5.29: the joint fit with the line, and the acceptance no row of it had seen."""
    two, line = json.loads((PRED / "crystal_stage_v13_run2.json").read_text()), json.loads((PRED / "crystal_stage_v13_line1.json").read_text())
    assert tuple(line["variables"]) == C.VARIABLES13 and "Not a preset" in line["stopped"]["why"]
    r = line["residuals"]
    helix = [k for k in r if k.startswith("lone chain, helix")]
    assert len(helix) == 6 and max(abs(r[k]["polyfind"] - r[k]["dft"]) for k in helix) < 0.05
    x1, x2 = dict(zip(line["variables"], line["x"])), dict(zip(two["variables"], two["x"]))
    # the torsion no longer leans toward cis, and the backbone angle is softer
    assert x2["V1"] < -1.5 and 0.2 < x1["V1"] < 0.7 and x1["log10_k_angle"] < x2["log10_k_angle"] - 0.15
    # what is between chains is kept
    rms = lambda q, pack: float(np.sqrt(np.mean([(v["polyfind"] - v["dft"]) ** 2 for k, v in q.items() if k.startswith(f"rigid {pack}_")])))  # noqa: E731
    for pack in ("beta", "candidate", "alpha", "hh", "cand6"):
        assert rms(r, pack) < rms(two["residuals"], pack) + 0.005
    # what it paid shows in alpha's cell and in the gap
    assert r["alpha b"]["z"] < -2.0 and r["E(alpha)-E(beta)"]["polyfind"] == pytest.approx(-0.651, abs=0.003)
    a = json.loads((PRED / "candidate_acceptance_stage13_line1.json").read_text())
    h = a["against_the_helix_dft_relaxed"]
    assert h["margin_upright"] == pytest.approx(0.041, abs=0.003) and abs(h["margin_upright"] - h["dft"]["margin_upright"]) < 0.05
    assert h["margin_sheared"] == pytest.approx(-0.197, abs=0.003) and a["all_trans_shear_gain"] == pytest.approx(-0.272, abs=0.003)
    assert h["other_helix_minus_it"] < 0.0 < h["dft"]["other_helix_minus_it"]  # the one row that still fails
    s = json.loads((PRED / "chain_split_polyfind_line1.json").read_text())["sets"]["stage 13 with the line"]
    dft = json.loads((PRED / "chain_split_pbe_d3bj.json").read_text())["pairs"]
    mine = s["vdf5-rvdf1-6"]["pairs"]["TTTT minus TG'TG"]
    ref = dft["cand_vdf5-rvdf1-6_TTTT_in_pvdfb_s2_upright minus cand_vdf5-rvdf1-6_TGTGp_in_pvdfb_s2"]
    # both parts of the level, not only their sum
    assert abs(mine["inside_a_chain"] - ref["inside_a_chain"]) < 0.05 and abs(mine["between_chains"] - ref["between_chains"]) < 0.08
    plain = s["pvdf"]["pairs"]["alpha minus beta"]
    assert plain["lone_chain"] == pytest.approx(dft["pvdf_000 minus pvdf_001"]["lone_chain"], abs=0.05)


def test_dfts_split_of_the_cells():
    d = json.loads((PRED / "chain_split_pbe_d3bj.json").read_text())
    for pair in d["pairs"].values():
        assert pair["cell"] == pytest.approx(pair["inside_a_chain"] + pair["between_chains"], abs=1e-9)
    plain = d["pairs"]["pvdf_000 minus pvdf_001"]
    assert plain["cell"] == pytest.approx(-0.726, abs=0.002) and plain["inside_a_chain"] == pytest.approx(-3.575, abs=0.003)
    assert plain["lone_chain"] == pytest.approx(-3.449, abs=0.003)
    cells = d["cells"]
    assert cells["pvdf_001"]["between"] == pytest.approx(-7.818, abs=0.003) and cells["pvdf_000"]["between"] == pytest.approx(-4.969, abs=0.003)
    # the lone chain at the cell's geometry is pulled on, mostly inside its CF2 groups; relaxed alone it gains a few tenths
    assert all(0.5 < c["largest_force_on_the_lone_chain"] < 1.1 for c in cells.values())
    assert cells["pvdf_001"]["relaxation"] == pytest.approx(-0.315, abs=0.005) and abs(cells["pvdf_001"]["lone_predicted"]) < 0.02
    # DFT's helix cells of the candidate: every gauche bond between 60 and 77 degrees, the reversed unit's two bonds trans in one and gauche in the other
    one, other = cells["cand_vdf5-rvdf1-6_TGTGp_in_pvdfb_s2"], cells["rvdf6_GTGpT_in_pvdfb_s2"]
    special = lambda c: [abs(v) for b, v in zip(c["bonds"], c["dihedrals"]) if b in ("CF2-CF2", "CH2-CH2")]  # noqa: E731
    assert min(special(one)) > 170 and max(special(other)) < 70
