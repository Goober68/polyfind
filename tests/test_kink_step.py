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
    (a force over 0.5).  A row's barriers on the settled saddles are the higher of its two."""
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
            if row[seg]["force_max"] < 0.5:
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


def test_saddle_polish_settles_a_band_top_into_a_stationary_point():
    """The first top of the stored two-chain band, polished: a point on the closure surface with no
    force left along the band's tangent nor across it (what a descent would feel there is under the
    tolerance), a little below the band's top, and above both of the band's end states."""
    row = json.load(open(MIG))["rows"][STACK]
    with PRESET_PARAMETERS["pvdf-dft-valence-belt4"].applied():
        model = _two_chain_model()
        free, closed = KM.free_sets(model, [0], "flexible")
        seg = row["A-I"]
        imgs, i = np.asarray(seg["images"]), seg["saddle_image"]
        e_a = model.physical_energy(imgs[0])
        x, info = KM.MC.saddle_polish(model, imgs[i], imgs[i + 1] - imgs[i - 1], free, closed, tol=0.01, patience=30, patience_tol=1e-6)
        assert info["converged"] and abs(info["force"]) < 0.01 and info["gradient_max"] < 0.02
        e = info["energy"] - e_a
        assert seg["saddle"] - 0.05 < e <= seg["saddle"] + 1e-6 and e > row["I"]["energy"] + 0.3
        assert max(model.closure_residual(x)) < 1e-6
        _, gmax, n = model.descend(x, free, closed, steps=1, gtol=1e9)
        assert n == 1 and gmax < 0.02  # the whole gradient along the closure surface, tangent included
        t = np.asarray(model.state(x).tors[0])
        assert KM.letters_of(t).split()[0] == "G+" and -170 < KM.wrap(t)[2] < -120 and -130 < KM.wrap(t)[4] < -90  # the leading gauche bond in flight


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
    assert "from the stored run" in log_band and "start: the movers' internal coordinates from the band of" in log_band
    assert np.asarray(seg["images"]).shape == (9, 4 + 2 * 44) and not seg["converged"] and seg["iterations"] == 3
    np.testing.assert_allclose(seg["images"][0], row["kinked_crystal"]["state"], atol=1e-12)
    assert seg["energies"][0] == 0.0 and seg["energies"][-1] == pytest.approx(row["I"]["energy"], abs=1e-6)
    assert "barrier_forward" not in row and "I-B" not in row  # one band of two: no barrier yet
    assert not (ck / "vdf10_0_1x1_mover_0_ahead_host_rigid.A-I.npy").exists()  # a band's checkpoint goes when the band is stored: one that is there is an unfinished run's
    assert {"I", "B"} <= set(row) and stored["rows"][STACK] == json.load(open(out))["rows"][STACK]  # nothing else was touched


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
