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
    log_b = run("--only", "B", "--patience", "25")
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
    assert (ck / "vdf10_0_1x1_mover_0_ahead_host_rigid.A-I.npy").exists()
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
