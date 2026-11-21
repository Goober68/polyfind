"""The potential fitted on with the fourteen-monomer labels (docs/MLIP_SECOND_OPINION.md 3.23): two fits scored on the 23 labelled cells beside all_s1."""

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
PRED = ROOT / "docs" / "predictions"
KCAL = 23.0605
FITS = ("all_s1", "refit_s1", "refitb_s1")


def _rows():
    spec = importlib.util.spec_from_file_location("cells14", ROOT / "tests" / "test_cells14_dft.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    out = []
    for stepped, specs in ((False, module.ROWS_OF_THE_STATES), (True, module.ROWS_OF_THE_STEPPED)):
        for item in specs:
            row, rest = item.split("=", 1)
            flat, *states = rest.split(",")
            out += [(row, label, stepped, flat, name) for label, name in (s.split(":", 1) for s in states)]
    return out


def _scores():
    return {fit: json.loads((PRED / "mlip" / "refit" / f"cells14_in_{fit}.json").read_text()) for fit in FITS}


def _errors(fit, scores):
    """Each state against its row's all-trans cell, kcal/mol per chain: the fit's less DFT's, the states held out apart."""
    frames, held = scores[fit]["frames"], set(scores[fit]["held_out"])
    out = {"in": {}, "out": {}}
    for row, label, stepped, flat, name in _rows():
        ours = (frames[name]["potential_eV_per_monomer"] - frames[flat]["potential_eV_per_monomer"]) * 14 * KCAL
        dft = (frames[name]["dft_eV_per_monomer"] - frames[flat]["dft_eV_per_monomer"]) * 14 * KCAL
        out["out" if name in held else "in"][(row, label, stepped)] = (ours, dft)
    return out


def _rms(pairs):
    return float(np.sqrt(np.mean([(a - b) ** 2 for a, b in pairs.values()])))


def test_five_frames_were_held_out_of_both_fits_and_eighteen_heavy_ones_were_in():
    scores = _scores()
    prepare = json.loads((PRED / "mlip" / "refit" / "prepare.json").read_text())
    assert prepare["frames"] == 621 and prepare["train"] == 539 and prepare["heavy"] == 23 and prepare["heavy_in_the_fit"] == 18
    assert prepare["held_out_whole"] == {"VAVF|VRVV": 8}                       # the composition all_s1 was not shown either
    for fit in FITS:
        assert len(scores[fit]["frames"]) == 23 and sorted(scores[fit]["held_out"]) == sorted(prepare["held_out_by_name"])
        assert sum(not f["in_the_fit"] for f in scores[fit]["frames"].values()) == 5
    for fit in FITS:   # DFT's energies are the stored labels' in each file
        assert [f["dft_eV_per_monomer"] for f in scores[fit]["frames"].values()] == [f["dft_eV_per_monomer"] for f in scores["all_s1"]["frames"].values()]


def test_neither_fit_is_better_than_all_s1_on_the_states_held_out():
    scores = _scores()
    errors = {fit: _errors(fit, scores) for fit in FITS}
    assert [len(errors[fit]["out"]) for fit in FITS] == [5, 5, 5] and [len(errors[fit]["in"]) for fit in FITS] == [12, 12, 12]
    held = [_rms(errors[fit]["out"]) for fit in FITS]
    shown = [_rms(errors[fit]["in"]) for fit in FITS]
    assert held == pytest.approx([1.06, 1.14, 1.18], abs=0.006) and shown == pytest.approx([0.74, 0.61, 0.62], abs=0.006)
    assert held[0] < min(held[1:])                                              # the rule set before the scores: the paths are scanned in all_s1
    # a state the fits were shown, weighted tenfold, and did not learn: the leader's sheet at 18m
    key = ("vdf9-trfe3-an1-rvdf1-14", "sheet_at_18m", False)
    assert [errors[fit]["in"][key][0] for fit in FITS] == pytest.approx([-1.37, -0.95, -0.91], abs=0.006) and errors["all_s1"]["in"][key][1] == pytest.approx(0.34, abs=0.006)
    # held out: the film cell's sheet at 12 stays 1.7 to 1.9 too low and the nitrile cell's sheet keeps the wrong sign
    key = ("vdf9-trfe4-cfe1-14", "sheet_at_12", False)
    assert [errors[fit]["out"][key][0] for fit in FITS] == pytest.approx([-3.60, -3.82, -3.73], abs=0.006) and errors["all_s1"]["out"][key][1] == pytest.approx(-1.93, abs=0.006)
    key = ("vdf9-trfe4-an1-14", "sheet_at_4", False)
    assert all(errors[fit]["out"][key][0] > 0.19 for fit in FITS) and errors["all_s1"]["out"][key][1] == pytest.approx(-0.43, abs=0.006)
    # the paths' own ends are in the fit, and all_s1 has the film cell's as near as either fit
    key = ("vdf9-trfe4-cfe1-14", "sheet_at_14m", True)
    assert [errors[fit]["in"][key][0] for fit in FITS] == pytest.approx([-5.03, -4.84, -4.43], abs=0.006) and errors["all_s1"]["in"][key][1] == pytest.approx(-5.27, abs=0.006)
    # absolute energies did come nearer DFT's: 3.4 meV per monomer rms to 2.5 and 2.6; a kcal/mol per chain is 3.1
    absolute = [1e3 * float(np.sqrt(np.mean([(f["potential_eV_per_monomer"] - f["dft_eV_per_monomer"]) ** 2 for f in scores[fit]["frames"].values()]))) for fit in FITS]
    assert absolute == pytest.approx([3.40, 2.47, 2.62], abs=0.01)


def test_both_fits_keep_the_anchor_and_the_frames_all_s1_held_back():
    own = {fit: json.loads((PRED / f"mlip_own_{fit}.json").read_text()) for fit in FITS}
    margin = [own[fit]["margins"]["on_the_dft_geometries"]["margin_kcal"] for fit in FITS]
    assert margin == pytest.approx([0.724, 0.731, 0.677], abs=0.001) and own["all_s1"]["margins"]["dft"]["margin_kcal"] == 0.726
    assert all(abs(m - 0.726) < 0.05 for m in margin)
    forces = [own[fit]["held_out"]["force_rmse_meV_per_A"] for fit in FITS]
    assert forces == pytest.approx([28.8, 27.5, 25.8], abs=0.1)
    assert [own[fit]["prepare"]["frames"] for fit in FITS] == [558, 621, 621]
    for fit in ("refit_s1", "refitb_s1"):   # 40 epochs each, side by side, the record's last lines the error tables
        text = (PRED / "mlip" / "refit" / f"{fit}.epochs.txt").read_text()
        assert "Epoch 35" in text and "Error-table on TEST" in text and "heavy" in text


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
