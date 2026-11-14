"""examples/kink_migration.py and the step tables of examples/kink_summary.py: one chain's kink pair one
monomer ahead of its sheet.  The pieces that can be checked without a relaxation are checked directly;
the driver's bookkeeping (one piece of a row per process, each writing its own fields; seeded starts;
checkpoints; a band started from another row's) is run once on the two-chain cell with rigid neighbours."""
import importlib.util
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from polyfind.fitting import PRESET_PARAMETERS

ROOT = Path(__file__).resolve().parents[1]


def _load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "examples" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


KM = _load("kink_migration")
KS = _load("kink_summary")
MIG = ROOT / "docs" / "predictions" / "kink_migration_belt4.json"
STACK = "vdf10@0 1x1 mover 0 ahead host flexible"  # the stored two-chain row, bands and all


def test_moving_a_kink_pair_carries_its_two_gauche_values_and_nothing_else():
    t = np.full(20, 180.0)
    t[[0, 2]] = 62.0, -65.0
    t[7] = 171.0  # a deflection that belongs to the sequence stays where it is
    assert KM.letters_of(t).split()[:5] == ["G+", "T", "G-", "T", "T"]
    ahead = KM.move_kink(t, 0, 2)
    assert KM.letters_of(ahead).split()[:5] == ["T", "T", "G+", "T", "G-"]
    assert (ahead[2], ahead[4], ahead[7]) == (62.0, -65.0, 171.0) and ahead[0] == ahead[1] == ahead[3] == 180.0
    back = KM.move_kink(ahead, 2, -2)
    np.testing.assert_array_equal(back, t)
    round_the_repeat = KM.move_kink(t, 0, -2)  # the repeat is periodic
    assert (round_the_repeat[18], round_the_repeat[0]) == (62.0, -65.0)


def test_step_response_reads_the_neighbours_displacements_from_the_stored_states():
    """A synthetic four-chain row: the moving chain shifted, one neighbour turned and another's torsion
    changed, a uniform axial shift of every chain that does not move (the held chain's share) taken out."""
    nc, mon = 4, 2
    P, R = 2 * mon, 2 * mon
    A, B, gam, c = 10.0, 12.0, 90.0, 5.0
    place = np.array([[0.0, 0.0, 10.0, 0.0], [0.25, 0.5, 100.0, 0.0], [0.5, 0.0, 10.0, 0.0], [0.75, 0.5, 100.0, 0.0]])
    inner = np.tile(np.concatenate([[60.0, 180.0, -60.0, 180.0], np.full(R, 112.0)]), (nc, 1))
    xA = np.concatenate([[A, B, gam, c], place.ravel(), inner.ravel()])
    p2, i2 = place.copy(), inner.copy()
    p2[0, 0] += 0.03          # the mover 0.3 A along a
    p2[0, 3] -= 0.40          # and 0.4 A back along its axis
    p2[[1, 2, 3], 3] += 0.10  # every other chain 0.1 A along the axis: the frame's, not theirs
    p2[2, 2] += 4.0           # chain 2 turned
    i2[1, 1] = 170.0          # chain 1: one torsion 10 degrees off
    i2[0, :P] = [180.0, 180.0, 70.0, 180.0]  # (a two-monomer toy: only the gauche list is read)
    row = {"chains": nc, "monomers_per_chain": mon, "mover": "0", "kinked_crystal": {"state": xA.tolist()},
           "B": {"state": np.concatenate([[A, B, gam, c], p2.ravel(), i2.ravel()]).tolist()}}
    r = KS.step_response(row, "B")
    assert r["mover"]["move_A"] == pytest.approx(0.30) and r["mover"]["axial_A"] == pytest.approx(-0.50)
    assert r["mover_gauche_deg"] == [70.0]
    shells = {round(s["distance_A"], 2): s for s in r["shells"]}
    assert sorted(shells) == [5.0, 6.5] and shells[6.5]["chains"] == 2 and shells[5.0]["chains"] == 1
    assert shells[5.0]["turn_deg"] == pytest.approx(4.0) and shells[5.0]["torsion_deg"] == pytest.approx(0.0)
    assert shells[6.5]["torsion_deg"] == pytest.approx(10.0) and shells[6.5]["axial_A"] == pytest.approx(0.0)
    assert shells[5.0]["move_A"] == pytest.approx(0.0, abs=1e-12)


def test_the_stored_sixteen_chain_end_states_are_the_states_they_are_named_for():
    """The single-chain step in sixteen chains (neighbours flexible): in I and in B the moving chain has
    its two gauche bonds where the step puts them, every other chain keeps its kink pair in the sheet,
    and the supercell check left the tiled sheet as the reference."""
    rows = json.load(open(MIG))["rows"]
    row = rows.get("vdf10@0 4x2 mover 0 ahead host flexible")
    if row is None or "I" not in row or "B" not in row:
        pytest.skip("the sixteen-chain row is not stored")
    nc, P = row["chains"], 2 * row["monomers_per_chain"]
    assert (nc, P) == (16, 20)
    check = row["kinked_crystal"]["supercell_check"]
    assert check["converged"] and abs(check["relaxed_above_tiled_per_chain"]) < 1e-6 and check["largest_torsion_change_deg"] < 0.01
    for name, bonds in (("I", (0, 4)), ("B", (2, 4))):
        rec = row[name]
        tors = np.asarray(rec["state"])[4 + 4 * nc:].reshape(nc, -1)[:, :P]
        letters = [KM.letters_of(t).split() for t in tors]
        assert [k for k, v in enumerate(letters[0]) if v != "T"] == list(bonds) and (letters[0][bonds[0]], letters[0][bonds[1]]) == ("G+", "G-")
        for other in letters[1:]:
            assert [k for k, v in enumerate(other) if v != "T"] == [0, 2]
        assert rec["converged"] and rec["gradient_max"] < 0.01
    assert row["B"]["energy"] < rows["vdf10@0 2x1 mover 0 ahead host flexible"]["B"]["energy"] - 1.0  # the cramped four-chain cell overstates B
    assert abs(row["I"]["energy"] - rows["vdf10@0 2x1 mover 0 ahead host flexible"]["I"]["energy"]) < 0.3  # and not I


def test_the_stored_sixteen_chain_bands_join_their_end_states_over_one_top_each():
    """Each band of the sixteen-chain row starts and ends on the stored end states, rises to one top
    and comes down again with no dip between, its top is where the moving chain's torsions say the
    gauche bond is in flight, and the barriers quoted for the row are the higher top against A and
    against B.  The completed step lies lower than in four chains and so does its barrier."""
    rows = json.load(open(MIG))["rows"]
    row, four = rows.get("vdf10@0 4x2 mover 0 ahead host flexible"), rows["vdf10@0 2x1 mover 0 ahead host flexible"]
    if row is None or "A-I" not in row or "I-B" not in row:
        pytest.skip("the sixteen-chain bands are not stored")
    nc = row["chains"]
    ends = {"A": (np.asarray(row["kinked_crystal"]["state"]), 0.0), "I": (np.asarray(row["I"]["state"]), row["I"]["energy"]), "B": (np.asarray(row["B"]["state"]), row["B"]["energy"])}
    periodic = np.zeros(4 + nc * 44, dtype=bool)
    periodic[6:4 + 4 * nc:4] = True  # the setting angles
    for s in range(nc):
        periodic[4 + 4 * nc + 40 * s:4 + 4 * nc + 40 * s + 20] = True  # the torsions: twenty of each chain's forty internal coordinates
    for seg, (a, b) in (("A-I", ("A", "I")), ("I-B", ("I", "B"))):
        e = row[seg]
        imgs, E = np.asarray(e["images"]), np.asarray(e["energies"])
        assert imgs.shape == (9, 4 + nc * 44) and e["force_max"] < 0.1 and not e["interior_minima"]  # 0.1: the criterion of every band of section 8.3
        for img, (x, energy), got in ((imgs[0], ends[a], E[0]), (imgs[-1], ends[b], E[-1])):
            diff = img - x
            diff[periodic] = (diff[periodic] + 180.0) % 360.0 - 180.0
            assert np.abs(diff).max() < 1e-8 and got == pytest.approx(energy, abs=1e-6)
        top = int(np.argmax(E))
        assert top == e["saddle_image"] and 0 < top < 8 and E[top] == pytest.approx(e["saddle"])
        assert np.all(np.diff(E[:top + 1]) > 0) and np.all(np.diff(E[top:]) < 0)
        np.testing.assert_array_equal(imgs[:, :4], np.tile(imgs[0, :4], (9, 1)))  # the cell is fixed along the band
    # in flight: at the first top the leading gauche bond is between bonds 2 and 4, at the second the trailing one has left bond 0
    t1, t2 = np.asarray(row["A-I"]["saddle_torsions_deg"]), np.asarray(row["I-B"]["saddle_torsions_deg"])
    assert 60 < t1[0] < 120 and -170 < t1[2] < -120 and -130 < t1[4] < -90
    assert 120 < t2[0] < 150 and abs(t2[2]) > 150 and -120 < t2[4] < -60
    top = max(row["A-I"]["saddle"], row["I-B"]["saddle"])
    assert row["barrier_forward"] == pytest.approx(top) and row["barrier_back"] == pytest.approx(top - row["B"]["energy"])
    assert row["I-B"]["saddle"] > row["A-I"]["saddle"] > row["I"]["energy"] > row["B"]["energy"] > 5.0
    assert row["barrier_forward"] < four["barrier_forward"] - 0.4 and row["barrier_back"] > four["barrier_back"] + 0.5


def test_the_stored_saddles_are_stationary_points_near_their_bands_tops():
    """Every settled saddle in the result file: no force left along the band's tangent nor across
    it, the curvature along the path negative where it was measured, above both ends of its band,
    and within 0.25 kcal/mol of the band's top unless the band had stopped far from convergence
    (a force over 0.5).  One saddle is 0.7 below a band that had met its criterion (the second
    chain's single-pass step: the band sat on the sequential side of a concerted pass); that one
    is held to more: descents from it either way end in the two states.  A row's barriers on the
    settled saddles are the higher of its two."""
    rows = json.load(open(MIG))["rows"]
    seen = 0
    for key, row in rows.items():
        ends = {"A": 0.0, "I": row.get("I", {}).get("energy"), "B": row.get("B", {}).get("energy")}
        for seg, (a, b) in (("A-I", ("A", "I")), ("I-B", ("I", "B"))):
            pol = row.get(seg, {}).get("polished")
            if not pol or not pol["converged"]:
                continue
            seen += 1
            assert abs(pol["force_along"]) < 0.01 and pol["gradient_across"] < 0.02, (key, seg)
            assert pol["curvature_along"] is None or pol["curvature_along"] < 0.0, (key, seg)
            assert pol["saddle"] > max(ends[a], ends[b]), (key, seg)
            assert pol["band_top"] == pytest.approx(row[seg]["saddle"]) and pol["band_force"] == pytest.approx(row[seg]["force_max"])
            if "connects" in pol:  # a saddle far from its band's top: shown by descents either way to join the two states it is said to
                ways = [pol["connects"][k]["largest_torsion_difference_from_the_stored_states_deg"] for k in ("towards_the_start", "towards_both_ahead")]
                assert ways[0]["A"] < 1.0 and ways[1]["B"] < 1.0 and ways[0]["B"] > 90.0 and ways[1]["A"] > 90.0, (key, seg)
                assert all(pol["connects"][k]["highest_energy_on_the_way"] < pol["saddle"] + 0.01 for k in ("towards_the_start", "towards_both_ahead")), (key, seg)
            elif row[seg]["force_max"] < 0.5:
                assert abs(pol["saddle"] - row[seg]["saddle"]) < 0.25, (key, seg)
        if all(row.get(seg, {}).get("polished", {}).get("converged") for seg in ("A-I", "I-B")):
            top = max(row["A-I"]["polished"]["saddle"], row["I-B"]["polished"]["saddle"])
            assert row["barrier_forward_polished"] == pytest.approx(top) and row["barrier_back_polished"] == pytest.approx(top - row["B"]["energy"])
    assert seen >= 9


def test_b_was_reached_twice_and_a_second_chain_beside_the_first_costs_less():
    """The sixteen-chain B from the sheet with the mover's torsions reset is the B relaxed from the
    rigid-chain-neighbour state; and two chains one monomer ahead side by side along the stack lie
    lower than two apart, both with their kink pairs where the step puts them."""
    rows = json.load(open(MIG))["rows"]
    one, two = rows.get("vdf10@0 4x2 mover 0 ahead host flexible"), rows.get("vdf10@0 4x2 mover 0,4 ahead host flexible")
    if one is None or two is None or "independent_start" not in one.get("B", {}):
        pytest.skip("the sixteen-chain pair or the second route to B is not stored")
    again = one["B"]["independent_start"]
    assert again["energy"] == pytest.approx(one["B"]["energy"], abs=1e-3) and again["start_energy"] > 100.0
    far = again["largest_difference_from_the_stored_B"]
    assert far["lateral_A"] < 0.01 and far["axial_A"] < 0.01 and far["torsion_deg"] < 0.5 and far["setting_angle_deg"] < 0.05
    assert two["B"]["converged"] and two["B"]["gradient_max"] < 0.01
    assert set(two["B"]["torsions_movers"]) == {"0", "4"}
    for letters in two["B"]["torsions_movers"].values():
        assert letters.split()[:5] == ["T", "T", "G+", "T", "G-"]
    tors = np.asarray(two["B"]["state"])[4 + 4 * 16:].reshape(16, 40)[:, :20]
    assert [s for s in range(16) if [k for k, v in enumerate(KM.letters_of(tors[s]).split()) if v != "T"] == [2, 4]] == [0, 4]
    assert 2 * two["B"]["energy"] < 2 * one["B"]["energy"] - 3.0  # side by side: more than 3 kcal/mol below two apart
    assert 2 * two["B"]["energy"] > one["B"]["energy"] + 1.0     # and the second chain is not free
    diag = rows.get("vdf10@0 4x2 mover 0,1 ahead host flexible")
    if diag is not None:  # a diagonal neighbour instead: cheaper than two apart, dearer than along the stack
        assert diag["B"]["converged"] and set(diag["B"]["torsions_movers"]) == {"0", "1"}
        assert all(v.split()[:5] == ["T", "T", "G+", "T", "G-"] for v in diag["B"]["torsions_movers"].values())
        assert two["B"]["energy"] < diag["B"]["energy"] < one["B"]["energy"] - 1.0


def test_the_completed_step_falls_as_the_cell_grows_and_the_way_back_does_not_change():
    """The stored rows of the single-chain step in cells of up to thirty-two chains: B falls with
    every doubling of the cell, by less each time along the stack; the wide kink does not move; in
    every B only the moving chain's kink pair has left the sheet.  And where the saddles were
    settled in the thirty-two-chain cell, the second has come down with B, so the barrier back is
    the sixteen-chain one to a tenth."""
    rows = json.load(open(MIG))["rows"]
    get = lambda cell: rows.get(f"vdf10@0 {cell} mover 0 ahead host flexible")  # noqa: E731
    cells = {c: get(c) for c in ("2x1", "4x1", "8x1", "2x2", "4x2", "8x2", "4x4")}
    if any(v is None for v in cells.values()):
        pytest.skip("the wider cells of the single-chain step are not stored")
    b = {c: v["B"]["energy"] for c, v in cells.items()}
    assert b["2x1"] > b["4x1"] > b["8x1"] and b["2x2"] > b["4x2"] > b["8x2"] and b["4x1"] > b["4x2"] > b["4x4"]
    assert b["2x1"] - b["4x1"] > 3.0 * (b["4x1"] - b["8x1"]) and b["2x2"] - b["4x2"] > 3.0 * (b["4x2"] - b["8x2"])  # closing along the stack
    assert 0.1 < b["4x2"] - b["8x2"] < 0.35 and 0.1 < b["4x2"] - b["4x4"] < 0.35
    wide = [v["I"]["energy"] for v in cells.values() if "I" in v]
    assert len(wide) >= 6 and max(wide) - min(wide) < 0.25
    for c in ("8x1", "8x2", "4x4"):
        row = cells[c]
        assert row["B"]["converged"] and row["B"]["gradient_max"] < 1e-3 and row["B"]["method"].startswith("on the closure surface")
        nc = row["chains"]
        tors = np.asarray(row["B"]["state"])[4 + 4 * nc:].reshape(nc, 40)[:, :20]
        where = [[k for k, v in enumerate(KM.letters_of(t).split()) if v != "T"] for t in tors]
        assert where[0] == [2, 4] and where[1:] == [[0, 2]] * (nc - 1)
    sd = cells["8x2"].get("saddles", {})
    if "I-B" in sd:
        small = cells["4x2"]
        assert sd["I-B"]["converged"] and abs(sd["I-B"]["force_along"]) < 0.01 and sd["I-B"]["gradient_across"] < 1e-3 and sd["I-B"]["curvature_along"] < 0.0
        fall = small["I-B"]["polished"]["saddle"] - sd["I-B"]["saddle"]
        assert 0.05 < fall < b["4x2"] - b["8x2"] + 0.05
        back16, back32 = small["I-B"]["polished"]["saddle"] - b["4x2"], sd["I-B"]["saddle"] - b["8x2"]
        assert abs(back32 - back16) < 0.1
    if "A-I" in sd:
        assert sd["A-I"]["converged"] and sd["A-I"]["saddle"] == pytest.approx(cells["4x2"]["A-I"]["polished"]["saddle"], abs=0.05)
        assert sd["A-I"]["saddle"] < sd.get("I-B", {}).get("saddle", np.inf)


def test_the_step_in_sixty_four_chains_and_the_thirty_two_chain_row_computed_twice():
    """The rows evaluated on a card in the batched form (docs/PERFORMANCE_REVIEW.md 17.9).  The
    thirty-two-chain row computed again from the same seeds is the host's to 1e-4 in all four
    numbers.  In sixty-four chains (8a x 4b) the completed step is where the two thirty-two-chain
    doublings put it if they add; the wide kink and the first saddle come down across the stack
    and not along it; the second saddle has NOT followed B down (reached by two routes, the same
    to 1e-5), is not additive in the two doublings, and stays within 0.2 over the four cells of
    sixteen to sixty-four chains, so the barrier ahead does not move and the way back grows.  In
    a hundred and twenty-eight chains B has come down again by a third to a quarter of the
    doubling before, each way."""
    rows = json.load(open(MIG))["rows"]
    get = lambda cell: rows.get(f"vdf10@0 {cell} mover 0 ahead host flexible")  # noqa: E731
    r16, r32, r32b, r64 = get("4x2"), get("8x2"), get("4x4"), get("8x4")
    if r64 is None or "again_on_device" not in (r32 or {}):
        pytest.skip("the rows evaluated on a card are not stored")
    again = r32["again_on_device"]
    for name in ("I", "B"):
        assert again[name]["converged"] and again[name]["energy"] == pytest.approx(r32[name]["energy"], abs=1e-6)
    for seg in ("A-I", "I-B"):
        assert again["saddles"][seg]["converged"] and again["saddles"][seg]["saddle"] == pytest.approx(r32["saddles"][seg]["saddle"], abs=1e-4)
    assert again["I"]["seconds"] < 0.05 * r32["I"]["seconds"] and "state" not in again["I"]

    nc = r64["chains"]
    assert nc == 64 and r64["supercell"] == [8, 4] and "batched" in r64["evaluated"]
    for name, gauche in (("I", [0, 4]), ("B", [2, 4])):
        st = r64[name]
        assert st["converged"] and st["gradient_max"] < 1e-3 and st["method"].startswith("on the closure surface")
        tors = np.asarray(st["state"])[4 + 4 * nc:].reshape(nc, 40)[:, :20]
        where = [[k for k, v in enumerate(KM.letters_of(t).split()) if v != "T"] for t in tors]
        assert where[0] == gauche and where[1:] == [[0, 2]] * (nc - 1)
    b = {c: v["B"]["energy"] for c, v in (("4x2", r16), ("8x2", r32), ("4x4", r32b), ("8x4", r64))}
    assert b["8x4"] < min(b["8x2"], b["4x4"])
    assert b["8x4"] == pytest.approx(b["8x2"] + b["4x4"] - b["4x2"], abs=0.01)  # the two doublings add
    assert 0.03 < r32["I"]["energy"] - r64["I"]["energy"] < 0.15
    sd = r64["saddles"]
    for seg in ("A-I", "I-B"):
        assert sd[seg]["converged"] and abs(sd[seg]["force_along"]) < 0.01 and sd[seg]["gradient_across"] < 1e-3 and sd[seg]["curvature_along"] < 0.0
        assert sd[seg]["from_cell"] == "8x2"
    assert 0.0 < r32["saddles"]["A-I"]["saddle"] - sd["A-I"]["saddle"] < 0.1
    tops = [r16["I-B"]["polished"]["saddle"], r32["saddles"]["I-B"]["saddle"], sd["I-B"]["saddle"]]
    assert sd["I-B"]["saddle"] > r32["saddles"]["I-B"]["saddle"] and max(tops) - min(tops) < 0.2
    assert r64["barrier_forward_polished"] == pytest.approx(sd["I-B"]["saddle"])
    assert r64["barrier_back_polished"] == pytest.approx(sd["I-B"]["saddle"] - b["8x4"])
    assert r64["barrier_back_polished"] > r32["barrier_back_polished"] + 0.15

    # the sixteen-chain cell doubled across the stack instead (4a x 4b): B is the host's, the rest is new
    assert "batched" in r32b["evaluated"] and r32b["on_host"]["B"]["energy"] == pytest.approx(b["4x4"], abs=1e-6)
    wide = {c: v["I"]["energy"] for c, v in (("4x2", r16), ("8x2", r32), ("4x4", r32b), ("8x4", r64))}
    assert 0.03 < wide["4x2"] - wide["4x4"] < 0.15 and 0.03 < wide["8x2"] - wide["8x4"] < 0.15     # across the stack
    assert abs(wide["4x2"] - wide["8x2"]) < 0.03 and abs(wide["4x4"] - wide["8x4"]) < 0.03         # along it
    for seg in ("A-I", "I-B"):
        assert r32b["saddles"][seg]["converged"] and abs(r32b["saddles"][seg]["force_along"]) < 0.01 and r32b["saddles"][seg]["curvature_along"] < 0.0
    other = r64["second_route"]["saddles"]["I-B"]
    assert other["converged"] and other["from_cell"] == "4x4" and "state" not in other
    assert other["saddle"] == pytest.approx(sd["I-B"]["saddle"], abs=1e-5)
    top = {"4x2": tops[0], "8x2": tops[1], "4x4": r32b["saddles"]["I-B"]["saddle"], "8x4": tops[2]}
    assert max(top.values()) - min(top.values()) < 0.2
    assert abs(top["8x4"] - (top["8x2"] + top["4x4"] - top["4x2"])) > 0.05  # where B's doublings add to 0.002

    r128a, r128b = get("16x4"), get("8x8")
    if r128a is None or r128b is None:
        return
    for row, small, before in ((r128a, b["4x4"], "along a"), (r128b, b["8x2"], "along b")):
        st, n = row["B"], row["chains"]
        assert n == 128 and st["converged"] and st["gradient_max"] < 1e-3
        tors = np.asarray(st["state"])[4 + 4 * n:].reshape(n, 40)[:, :20]
        where = [[k for k, v in enumerate(KM.letters_of(t).split()) if v != "T"] for t in tors]
        assert where[0] == [2, 4] and where[1:] == [[0, 2]] * (n - 1), before
        gain, gain_before = b["8x4"] - st["energy"], small - b["8x4"]
        assert 0.04 < gain < 0.08 and 2.5 < gain_before / gain < 5.0, before


def test_a_stepped_chain_is_recalled_not_joined():
    """The stored rows of the second chain's step (sixteen chains, energies against the state with
    chain 0 one monomer on).  Both stack neighbours end 2.44 above that start, the stored pair
    state; the one that has the stepped chain at -a goes through a wide kink (found from two
    starts) and two saddles, the second the higher; the one that has it at +a has no wide kink: one
    saddle, and from the wide-kink geometry a descent in short steps that never rises and ends in
    the completed step.  Either way the neighbour's barrier to join is far above the first chain's
    barrier to come back, and a pair parts more slowly than one chain returns."""
    rows = json.load(open(MIG))["rows"]
    one, pair = rows.get("vdf10@0 4x2 mover 0 ahead host flexible"), rows.get("vdf10@0 4x2 mover 0,4 ahead host flexible")
    hard, easy = rows.get("vdf10@0 4x2 mover 4 ahead host flexible beside 0"), rows.get("vdf10@0 4x2 mover 12 ahead host flexible beside 0")
    if any(v is None for v in (one, pair, hard, easy)):
        pytest.skip("the second chain's step is not stored")
    back = one["barrier_back_polished"]
    for row, mover in ((hard, 4), (easy, 12)):
        assert row["beside"]["chains"] == [0] and row["beside"]["above_the_sheet"] == pytest.approx(one["B"]["energy"], abs=1e-9)
        np.testing.assert_allclose(row["beside"]["state"], one["B"]["state"], atol=1e-12)
        assert row["B"]["converged"] and row["B"]["energy"] + row["beside"]["above_the_sheet"] == pytest.approx(2 * pair["B"]["energy"], abs=1e-3)
        nc = row["chains"]
        tors = np.asarray(row["B"]["state"])[4 + 4 * nc:].reshape(nc, 40)[:, :20]
        where = [[k for k, v in enumerate(KM.letters_of(t).split()) if v != "T"] for t in tors]
        assert [s for s, w in enumerate(where) if w == [2, 4]] == sorted([0, mover]) and sum(1 for w in where if w == [0, 2]) == nc - 2
        pol = row["A-I"]["polished"]
        assert pol["converged"] and abs(pol["force_along"]) < 0.01 and pol["gradient_across"] < 2e-3 and pol["curvature_along"] < 0.0
        assert pol["method"].startswith("on the closure surface") and pol["saddle"] < row["A-I"]["saddle"] + 1e-6
        assert row["barrier_back_polished"] == pytest.approx(row["barrier_forward_polished"] - row["B"]["energy"])
        assert row["barrier_forward_polished"] > back + 5.0  # joined far more rarely than recalled
        assert row["barrier_back_polished"] > back + 2.0     # and a pair holds: it parts more slowly than one chain comes back
    assert hard["B"]["energy"] == pytest.approx(easy["B"]["energy"], abs=1e-4)  # the two pairs are translates of each other
    # the neighbour that has the stepped chain at -a: a wide kink, by two starts, and two saddles
    assert hard["I"]["converged"] and hard["I"]["torsions"].split()[:5] == ["G+", "T", "T", "T", "G-"] and "not_a_minimum" not in hard["I"]
    again = hard["I"]["independent_start"]
    assert again["energy"] == pytest.approx(hard["I"]["energy"], abs=1e-4) and again["start_energy"] > 100.0 and again["largest_difference_from_the_stored_I"]["torsion_deg"] < 0.05
    second = hard["I-B"]["polished"]
    assert second["converged"] and abs(second["force_along"]) < 0.01 and second["gradient_across"] < 1e-3 and second["curvature_along"] < 0.0
    assert second["saddle"] > hard["A-I"]["polished"]["saddle"] > hard["I"]["energy"] > hard["B"]["energy"]
    assert hard["barrier_forward_polished"] == pytest.approx(second["saddle"]) and second["saddle"] > one["I-B"]["polished"]["saddle"] + 1.0
    # the neighbour that has it at +a: no wide kink, one saddle, then downhill all the way
    assert easy["I"]["not_a_minimum"] and easy["single_pass"] and "I-B" not in easy
    assert easy["I_relaxed_from_the_wide_kink"]["torsions"].split()[:5] == ["T", "T", "G+", "T", "G-"]  # relaxed from the wide kink it ends in the completed step
    assert easy["I_relaxed_from_the_wide_kink"]["energy"] == pytest.approx(easy["B"]["energy"], abs=1e-3)
    slide = easy["slide"]
    e = slide["energies_every_fifth_evaluation"]
    assert all(b <= a + 1e-9 for a, b in zip(e, e[1:])) and e[0] == pytest.approx(slide["start_energy"]) and e[0] > 9.0
    assert e[-1] == pytest.approx(easy["B"]["energy"], abs=1e-3) and abs(slide["end_against_the_stored_B"]["energy_difference"]) < 1e-3
    assert slide["end_against_the_stored_B"]["largest_torsion_difference_deg"] < 0.05
    turn = slide["trailing_gauche_bond_deg_every_fifth_evaluation"]
    assert turn[0] < 100.0 and turn[-1] > 170.0  # the trailing gauche bond goes from where a wide kink has it to trans
    assert easy["barrier_forward_polished"] == pytest.approx(easy["A-I"]["polished"]["saddle"]) and easy["I"]["energy"] < easy["A-I"]["polished"]["saddle"] - 2.0
    assert easy["barrier_forward_polished"] < hard["barrier_forward_polished"] - 1.5
    t = np.asarray(easy["A-I"]["polished"]["torsions_deg"])
    assert 120.0 < t[0] < 160.0 and (t[2] < -150.0 or t[2] > 150.0) and -140.0 < t[4] < -100.0  # at the one saddle both gauche bonds are on their way


@pytest.mark.parametrize("how", ["slsqp", "box", "surface"])
def test_saddle_polish_settles_a_band_top_into_a_stationary_point(how):
    """The first top of the stored two-chain band, polished: a point on the closure surface with no
    force left along the band's tangent nor across it (what a descent would feel there is under the
    tolerance), a little below the band's top, and above both of the band's end states.  With each
    minimisation confined to a box about its start the saddle is the same; and so it is with the
    minimisations done on the closure surface in the per-chain metric, in a fraction of the
    evaluations and with less force left across the path."""
    row = json.load(open(MIG))["rows"][STACK]
    with PRESET_PARAMETERS["pvdf-dft-valence-belt4"].applied():
        model = _two_chain_model()
        free, closed = KM.free_sets(model, [0], "flexible")
        seg = row["A-I"]
        imgs, i = np.asarray(seg["images"]), seg["saddle_image"]
        e_a = model.physical_energy(imgs[0])
        kw = {"slsqp": {}, "box": {"trust": 0.5}}.get(how)
        if kw is None:
            blocks = KM.MC.chain_blocks(model, imgs[0], [0, 1])
            kw = {"on_surface": True, "metric": KM.MC.chain_metric(model, imgs[0], free, closed, blocks)}
        x, info = KM.MC.saddle_polish(model, imgs[i], imgs[i + 1] - imgs[i - 1], free, closed, tol=0.01, patience=30, patience_tol=1e-6, **kw)
        assert info["converged"] and abs(info["force"]) < 0.01 and info["gradient_max"] < (1e-3 if how == "surface" else 0.02)
        e = info["energy"] - e_a
        assert e == pytest.approx(row["A-I"]["polished"]["saddle"], abs=2e-3)  # the stored one, found without a box
        assert seg["saddle"] - 0.05 < e <= seg["saddle"] + 1e-6 and e > row["I"]["energy"] + 0.3
        assert max(model.closure_residual(x)) < 1e-6
        if how == "surface":
            assert info["evaluations"] < 40 and all("success" in r["message"].lower() for r in info["rounds"])
        _, gmax, n = model.descend(x, free, closed, steps=1, gtol=1e9)
        assert n == 1 and gmax < 0.02  # the whole gradient along the closure surface, tangent included
        t = np.asarray(model.state(x).tors[0])
        assert KM.letters_of(t).split()[0] == "G+" and -170 < KM.wrap(t)[2] < -120 and -130 < KM.wrap(t)[4] < -90  # the leading gauche bond in flight


def test_the_sixteen_chain_step_was_settled_on_the_closure_surface():
    """The sixteen-chain row's two saddles, as stored: each found with every hyperplane minimised on
    the closure surface to its own test (no stall), a force across the path a hundred times smaller
    than SLSQP left, and not above what SLSQP gave (kept beside it as ``earlier``); and its two end
    states relaxed further the same way fall by less than a thousandth, so the energies the bands
    start and end on stand."""
    row = json.load(open(MIG))["rows"].get("vdf10@0 4x2 mover 0 ahead host flexible")
    if row is None or not all("method" in row.get(seg, {}).get("polished", {}) for seg in ("A-I", "I-B")):
        pytest.skip("the sixteen-chain saddles settled on the closure surface are not stored")
    for seg in ("A-I", "I-B"):
        pol = row[seg]["polished"]
        assert pol["method"].startswith("on the closure surface") and pol["converged"]
        assert abs(pol["force_along"]) < 0.01 and pol["gradient_across"] < 1e-3 and pol["curvature_along"] < 0.0
        assert all("success" in r["message"].lower() for r in pol["rounds"])
        before = pol["earlier"][-1]
        assert before["method"].startswith("SLSQP") if "method" in before else True
        assert pol["saddle"] < before["saddle"] + 2e-3 and pol["evaluations"] < before["evaluations"]
    assert row["I-B"]["polished"]["saddle"] > row["A-I"]["polished"]["saddle"] > row["I"]["energy"]
    again = row["A-I"]["polished"].get("found_again")
    if again is not None:  # a search started from a concerted geometry ends on the same first saddle: the single chain steps in two passes
        assert again["energy"] == pytest.approx(row["A-I"]["polished"]["saddle"], abs=1e-3) and again["start_energy"] > 20.0
        assert again["largest_torsion_difference_from_the_stored_saddle_deg"] < 1.0 and 80.0 < again["torsions_deg"][0] < 100.0
    for name in ("I", "B"):
        check = row[name]["relaxed_further"]
        assert 0.0 <= check["fall"] < 1e-3 and check["gradient_max"] < 2e-4
        assert check["largest_torsion_change_deg"] < 1.0  # soft: B's torsions move by up to 0.7 degrees for a fall of 0.0001
        assert "success" in check["message"].lower()


def test_the_pieces_of_a_row_run_one_per_process_and_merge(tmp_path):
    """The driver on the two-chain cell with the other chain fixed: I and B from separate processes into
    one row (each keeps the other's fields), I seeded from another row's and B checkpointed; then one
    band, started from the stored band of that row, a few iterations; the checkpoints are written and a
    second call of a finished piece reads it back instead of recomputing."""
    out = tmp_path / "step.json"
    stored = json.load(open(MIG))
    stored["rows"] = {STACK: stored["rows"][STACK]}
    out.write_text(json.dumps(stored))
    ck = tmp_path / "ck"
    env = dict(os.environ, OPENBLAS_NUM_THREADS="1")
    base = [sys.executable, str(ROOT / "examples" / "kink_migration.py"), "--row", "vdf10", "--nx", "1", "--ny", "1", "--mover", "0", "--host", "rigid",
            "--references", str(ROOT / "docs" / "predictions" / "dilute_kink_belt4.json"), "--out", str(out), "--checkpoints", str(ck)]

    def run(*extra):
        res = subprocess.run(base + list(extra), capture_output=True, text=True, env=env, cwd=ROOT, timeout=900)
        assert res.returncode == 0, res.stdout[-2000:] + res.stderr[-2000:]
        return res.stdout

    key = "vdf10@0 1x1 mover 0 ahead host rigid"
    log_i = run("--only", "I", "--seed-from", STACK)
    row = json.load(open(out))["rows"][key]
    assert "I" in row and "B" not in row and row["I"]["start"].startswith("I of the row") and row["kinked_crystal"]["message"].startswith("the polished")
    assert "end points only" in log_i
    run("--only", "B", "--patience", "25")
    row = json.load(open(out))["rows"][key]
    assert {"I", "B", "kinked_crystal", "chains", "host"} <= set(row) and row["host"] == "rigid" and row["chains"] == 2
    assert row["I"]["converged"] and row["B"]["converged"] and row["B"]["gradient_max"] < 0.01
    assert row["I"]["torsions"].split()[:5] == ["G+", "T", "T", "T", "G-"] and row["B"]["torsions"].split()[:5] == ["T", "T", "G+", "T", "G-"]
    # the other chain cannot answer: both states lie above the ones with every chain free
    flexible = stored["rows"][STACK]
    assert row["I"]["energy"] > flexible["I"]["energy"] - 1e-6 and row["B"]["energy"] > flexible["B"]["energy"]
    assert (ck / "vdf10_0_1x1_mover_0_ahead_host_rigid.B.npy").exists()
    log_band = run("--only", "A-I", "--band-from", STACK, "--maxiter", "3")
    row = json.load(open(out))["rows"][key]
    seg = row["A-I"]
    assert "from the stored run" in log_band and "start: the straight band bent as the band of" in log_band
    assert np.asarray(seg["images"]).shape == (9, 4 + 2 * 44) and not seg["converged"] and seg["iterations"] == 3
    np.testing.assert_allclose(seg["images"][0], row["kinked_crystal"]["state"], atol=1e-12)
    assert seg["energies"][0] == 0.0 and seg["energies"][-1] == pytest.approx(row["I"]["energy"], abs=1e-6)
    assert "barrier_forward" not in row and "I-B" not in row  # one band of two: no barrier yet
    assert not (ck / "vdf10_0_1x1_mover_0_ahead_host_rigid.A-I.npy").exists()  # a band's checkpoint goes when the band is stored: one that is there is an unfinished run's
    assert {"I", "B"} <= set(row) and stored["rows"][STACK] == json.load(open(out))["rows"][STACK]  # nothing else was touched
    # --beside: the other chain's step from the state in which this one has stepped.  In the two-chain cell that completes the
    # sheet's own move by one monomer, so it runs downhill from its start; not all the way to the sheet's energy, because the
    # first chain is held as it was when it stood alone.
    log_2 = subprocess.run(base[:base.index("--mover")] + ["--mover", "1", "--beside", "0"] + base[base.index("--mover") + 2:] + ["--only", "B", "--patience", "25"],
                           capture_output=True, text=True, env=env, cwd=ROOT, timeout=900)
    assert log_2.returncode == 0, log_2.stdout[-2000:] + log_2.stderr[-2000:]
    rows = json.load(open(out))["rows"]
    second = rows["vdf10@0 1x1 mover 1 ahead host rigid beside 0"]
    assert second["beside"]["chains"] == [0] and second["beside"]["row"] == key and second["beside"]["above_the_sheet"] == pytest.approx(row["B"]["energy"], abs=1e-9)
    np.testing.assert_allclose(second["beside"]["state"], row["B"]["state"], atol=1e-12)
    assert "energies below are against it" in log_2.stdout and second["B"]["converged"] and second["B"]["torsions"].split()[:5] == ["T", "T", "G+", "T", "G-"]
    assert second["B"]["energy"] < -1.0 and -1e-6 < second["beside"]["above_the_sheet"] + second["B"]["energy"] < second["beside"]["above_the_sheet"]
    assert rows[key] == row  # the first chain's row is read, not written
    refused = subprocess.run(base[:base.index("--mover")] + ["--mover", "0", "--beside", "0"] + base[base.index("--mover") + 2:] + ["--only", "B"],
                             capture_output=True, text=True, env=env, cwd=ROOT, timeout=900)
    assert refused.returncode != 0 and "cannot already have stepped" in refused.stderr


def test_free_sets_name_what_each_host_lets_move():
    model = _two_chain_model()
    n_int = model.P + model.R
    free, closed = KM.free_sets(model, [0], "rigid")
    assert closed == [0] and len(free) == n_int + 4  # the mover's internal coordinates and its placement; the other chain is the anchor
    free_b, closed_b = KM.free_sets(model, [0], "bodies")
    assert closed_b == [0] and set(free) < set(free_b) and len(free_b) == len(free) + 1  # the anchor may turn, and that is all a second chain adds
    free_f, closed_f = KM.free_sets(model, [0], "flexible")
    assert closed_f == [0, 1] and len(free_f) == len(free_b) + n_int
    assert not set(free_f) & set(model.index("lattice").tolist() + model.index("c").tolist())  # the cell is fixed in every case


def _two_chain_model():
    stored = json.load(open(MIG))["kinked_references"]["vdf10@0"]
    with PRESET_PARAMETERS["pvdf-dft-valence-belt4"].applied():
        letters, _ = KM.DK.sequence("vdf10")
        polymer = KM.DK.C.make_polymer("vdf10", letters)
        seq = [1, 0, 2] + [0] * (polymer.bonds_per_repeat - 3)
        d = KM.DK.refined_cell(polymer, "vdf10", seq, "vdf10 kink@0", stored)
        return KM.DK.make_model(polymer, d, 1, 1)[0]


def test_a_state_let_into_a_wider_cell_keeps_the_moving_chains_neighbours_and_lets_the_sheet_in():
    """examples/kink_wide_cell.py: the doubling behind every row of the cell-size series.  The sheet
    let into itself is the stored wider sheet; a completed step let in keeps chain 0's own layer and
    the layers beside it on both sides, and is the sheet where the narrower cell is farthest away."""
    KW = _load("kink_wide_cell")
    rows = json.loads(MIG.read_text())["rows"]
    key = lambda cell: f"vdf10@0 {cell} mover 0 ahead host flexible"
    sheet = np.asarray(rows[key("4x2")]["kinked_crystal"]["state"])
    stepped = np.asarray(rows[key("4x2")]["B"]["state"])
    for cell, along in (("8x2", "a"), ("4x4", "b")):
        assert np.array_equal(KW.let_in(sheet, sheet, 4, 2, along), np.asarray(rows[key(cell)]["kinked_crystal"]["state"]))
    wide_sheet = KW.let_in(sheet, sheet, 4, 2, "a")
    wide = KW.let_in(stepped, sheet, 4, 2, "a")
    nc, width = 16, (sheet.size - 4 - 4 * 16) // 16
    inner = lambda x, n: x[4 + 4 * n:].reshape(n, width)
    per = 4  # chains per layer along a in the 4x2 cell
    for layer2, layer in ((0, 0), (1, 1), (6, 2), (7, 3)):  # chain 0's side, and its far side through the boundary
        assert np.array_equal(inner(wide, 32)[layer2 * per:(layer2 + 1) * per], inner(stepped, nc)[layer * per:(layer + 1) * per])
    for layer2 in (2, 3, 4, 5):  # the sheet, let in between
        assert np.array_equal(inner(wide, 32)[layer2 * per:(layer2 + 1) * per], inner(wide_sheet, 32)[layer2 * per:(layer2 + 1) * per])
    assert wide[0] == 2.0 * sheet[0] and np.array_equal(wide[1:4], sheet[1:4])
    both, both_sheet, nx, ny = KW.through(stepped, sheet, 4, 2, ["a", "b"])
    assert (nx, ny) == (8, 4) and both.size == 4 + 64 * (4 + width)
    assert np.array_equal(both_sheet, KW.let_in(wide_sheet, wide_sheet, 8, 2, "b"))
    with pytest.raises(ValueError, match="no far side"):
        KW.let_in(np.zeros(4 + 4 * (4 + width)), np.zeros(4 + 4 * (4 + width)), 2, 1, "b")


def test_a_job_table_starts_from_the_stored_kinked_crystals_and_a_held_state_is_found_by_its_path(tmp_path):
    """examples/kink_wide_cell.py: ``table`` gives a job a table with the stored kinked crystals and
    none of the stored rows but those named (so a row that is stored is computed again), and
    ``held_state`` reads a state out of a row by a dotted path or out of the settled records."""
    import types

    KW = _load("kink_wide_cell")
    stored = json.loads(MIG.read_text())
    key = "vdf10@0 2x1 mover 0 ahead host flexible"
    out = tmp_path / "steps.json"
    KW.table(types.SimpleNamespace(source=str(MIG), out=str(out), row=None))
    first = json.loads(out.read_text())
    assert first["rows"] == {} and first["kinked_references"] == stored["kinked_references"] and first["preset"] == stored["preset"]
    KW.table(types.SimpleNamespace(source=str(MIG), out=str(out), row=[key]))
    KW.table(types.SimpleNamespace(source=str(MIG), out=str(out), row=[key]))  # again: nothing changes
    second = json.loads(out.read_text())
    assert list(second["rows"]) == [key] and second["rows"][key] == stored["rows"][key]
    row = stored["rows"][key]
    assert np.array_equal(KW.held_state(row, key, "I.state"), np.asarray(row["I"]["state"]))
    assert np.array_equal(KW.held_state(row, key, "I-B.images.3"), np.asarray(row["I-B"]["images"][3]))
    assert np.array_equal(KW.held_state(row, key, "A-I.polished.state"), np.asarray(row["A-I"]["polished"]["state"]))
    settled = tmp_path / "states.json"
    settled.write_text(json.dumps({"settled": {key: {"M": {"state": [1.0, 2.0]}}}}))
    assert KW.held_state(row, key, "settled:M", settled).tolist() == [1.0, 2.0]


def test_one_copolymer_chain_alone_pays_three_more_than_a_sublattice_and_the_chemistries_keep_their_distances():
    """The single-chain step of the terpolymer (64/29/7) and of its matrix (71/29) in four and in
    sixteen chains, computed on a card (docs/CANDIDATE_MATERIALS.md 8.3, 2026-11-14), with PVDF's
    sixteen-chain row computed again the same way as the control.  The end states are the states
    they are named for; from the two-chain protocol to one chain in sixteen every barrier rises
    by about three, so the differences between the materials stay; the terpolymer's sixteen-chain
    barrier is the route through the intermediate state M, its direct band's top not being a
    settled saddle; and the states M and N are what the stored file says they are."""
    rows = json.load(open(MIG))["rows"]
    get = lambda row, cell: rows.get(f"{row}@0 {cell} mover 0 ahead host flexible")  # noqa: E731
    ter, mat, pvdf = "vdf9-trfe4-cfe1-14", "vdf10-trfe4-14", "vdf10"
    states_path = MIG.with_name("kink_states_belt4.json")
    if get(ter, "4x2") is None or get(mat, "4x2") is None or not states_path.exists():
        pytest.skip("the copolymers' single-chain rows are not stored")
    states = json.loads(states_path.read_text())
    key = lambda row, cell: f"{row}@0 {cell} mover 0 ahead host flexible"  # noqa: E731

    again, host = get(pvdf, "4x2")["again_on_device"], get(pvdf, "4x2")  # the control: the protocol from nothing, on the card
    for name in ("I", "B"):
        assert again[name]["converged"] and again[name]["energy"] == pytest.approx(host[name]["energy"], abs=2e-4) and "state" not in again[name]
    for seg in ("A-I", "I-B"):
        assert again[seg]["converged"] and "images" not in again[seg] and again[seg]["polished"]["converged"]
        assert again[seg]["polished"]["saddle"] == pytest.approx(host[seg]["polished"]["saddle"], abs=2e-4)
        assert again[seg]["seconds"] < 0.05 * host[seg]["seconds"]

    for row in (ter, mat):
        for cell in ("2x1", "4x2"):
            r = get(row, cell)
            nc = r["chains"]
            assert nc == {"2x1": 4, "4x2": 16}[cell] and r["monomers_per_chain"] == 14 and "batched" in r["evaluated"]
            for name, gauche in (("I", [0, 4]), ("B", [2, 4])):
                st = r[name]
                assert st["converged"] and st["gradient_max"] < 1e-3 and st["start"].startswith("the kinked crystal")
                tors = np.asarray(st["state"])[4 + 4 * nc:].reshape(nc, 56)[:, :28]
                where = [[k for k, v in enumerate(KM.letters_of(t).split()) if v != "T"] for t in tors]
                assert where[0] == gauche and where[1:] == [[0, 2]] * (nc - 1)
            first = r["A-I"]["polished"]
            assert first["converged"] and 0.0 < first["saddle"] - r["I"]["energy"] < 0.25  # the wide kink is barely a minimum
            assert r["barrier_back_polished"] == pytest.approx(r["barrier_forward_polished"] - r["B"]["energy"])
    check = get(ter, "2x1")["kinked_crystal"]["supercell_check"]
    assert check["converged"] and abs(check["relaxed_above_tiled_per_chain"]) < 1e-6
    assert abs(get(mat, "2x1")["kinked_crystal"]["supercell_check"]["relaxed_above_tiled_per_chain"]) < 1e-3

    bar = {(row, cell): get(row, cell)["barrier_forward_polished"] for row in (ter, mat, pvdf) for cell in ("1x1", "2x1", "4x2")}
    for row in (ter, mat, pvdf):
        assert 2.5 < bar[row, "4x2"] - bar[row, "1x1"] < 3.3
    assert (bar[mat, "4x2"] - bar[ter, "4x2"]) == pytest.approx(bar[mat, "1x1"] - bar[ter, "1x1"], abs=0.1)
    assert all(2.4 < bar[ter, cell] - bar[pvdf, cell] < 3.0 for cell in ("1x1", "4x2"))
    assert all(get(row, "4x2")["B"]["energy"] > 1.7 * get(pvdf, "4x2")["B"]["energy"] for row in (ter, mat))
    assert all(bar[row, "4x2"] > bar[row, "2x1"] and get(row, "4x2")["B"]["energy"] < get(row, "2x1")["B"]["energy"] for row in (ter, mat))

    r = get(ter, "4x2")  # its barrier is the route's: the direct band's top is not a settled saddle
    assert "polished" not in r["I-B"] and not r["I-B"]["search_left_the_path"]["converged"] and r["I-B"]["converged"]
    way = r["by_way_of_M"]
    assert r["barrier_forward_polished"] == pytest.approx(max(way["saddles"].values())) and way["barrier_forward"] < r["I-B"]["saddle"] - 0.5
    assert all(v < 0.0 for v in way["curvature_along"].values()) and "by way of" in r["barrier_polished_note"]
    for row, cell in ((ter, "2x1"), (mat, "4x2")):  # where the direct band settled, the route through M ends on the same saddle
        assert get(row, cell)["by_way_of_M"]["barrier_forward"] == pytest.approx(get(row, cell)["I-B"]["polished"]["saddle"], abs=2e-3)

    for row, cell in ((ter, "2x1"), (ter, "4x2"), (mat, "2x1"), (mat, "4x2")):
        m, r = states["settled"][key(row, cell)]["M"], get(row, cell)
        t = m["torsions_deg"]
        assert m["converged"] and -110.0 < t[4] < -90.0 and 110.0 < t[2] < 150.0 and 150.0 < abs(t[0]) and r["I"]["energy"] < m["energy"]
    for row in (ter, mat, pvdf):
        n, r = states["settled"][key(row, "4x2")]["N"], get(row, "4x2")
        assert n["converged"] and n["gauche_bonds"] == [4, 6] and n["torsions_deg"][4] < 0.0 < n["torsions_deg"][6] and n["chains_off_the_sheet"].keys() == {"0"}
        assert (n["energy"] < r["I"]["energy"]) and ((n["energy"] < r["B"]["energy"]) == (row != pvdf))
        band = states["bands"][key(row, "4x2")]["I-N"]  # not converged: a bound, not a saddle
        assert not band["converged"] and not band["polished"]["a_saddle_between_the_two"] and band["saddle"] > max(r["I"]["energy"], n["energy"]) + 3.0
    legs = states["bands"][key(ter, "4x2")]
    assert legs["I-M"]["polished"]["a_saddle_between_the_two"] and legs["M-B"]["polished"]["a_saddle_between_the_two"]
