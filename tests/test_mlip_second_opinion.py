"""The pretrained machine-learned potentials against the one margin held from DFT.

Reads the stored results only: torch, mace and ase are not dependencies of the package.
`examples/mlip_second_opinion.py` regenerates the file.
"""

import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "docs" / "predictions" / "mlip_second_opinion.json"
ANCHOR = ROOT / "deliverables" / "polar_margin_dft" / "polar_margin_dft.json"
KCAL_PER_MEV = 0.023060548


def _script():
    spec = importlib.util.spec_from_file_location("mlip_second_opinion", ROOT / "examples" / "mlip_second_opinion.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_the_test_quotes_the_anchor_as_the_dft_deliverable_has_it():
    quoted = _script().DFT
    gap = json.loads(ANCHOR.read_text())["margins"]["alpha - polar (beta)"]
    assert quoted["alpha_minus_beta_meV"] == pytest.approx(gap["total"]["scf_k1"]["meV"], abs=0.005)
    assert quoted["alpha_minus_beta_d3_part_meV"] == pytest.approx(gap["d3_part"]["scf_k1"]["meV"], abs=0.005)
    assert quoted["alpha_minus_beta_pbe_part_meV"] == pytest.approx(
        gap["total"]["scf_k1"]["meV"] - gap["d3_part"]["scf_k1"]["meV"], abs=0.005)
    assert quoted["margin_kcal"] == pytest.approx(-gap["total"]["scf_k1"]["meV"] * KCAL_PER_MEV, abs=0.001)
    anti = json.loads(ANCHOR.read_text())["margins"]["anti_free - polar"]["total"]["scf_k1"]["meV"]
    assert quoted["anti_minus_beta_meV"] == pytest.approx(anti, abs=0.005)


def test_the_dispersion_add_on_alone_is_the_dft_margins_dispersion_part():
    d3 = json.loads(RESULT.read_text())["d3_bj_alone"]
    assert d3["alpha_minus_beta_meV"] == pytest.approx(d3["dft_d3_part_meV"], abs=0.1)


def test_no_pretrained_potential_resolves_the_margin():
    """Each is off by more than ten times the 0.1 kcal/mol needed, on both readings."""
    result = json.loads(RESULT.read_text())
    needed = result["needed_accuracy_kcal"]
    assert needed == 0.1
    models = result["models"]
    assert {"medium-mpa-0", "medium-0b3", "medium", "egret-1", "egret-1t"} <= set(models)
    for name, m in models.items():
        for reading in ("on_the_dft_geometries", "relaxed_in_the_potential"):
            error = abs(m[reading]["margin_error_kcal"])
            assert error > 10 * needed, (name, reading, error)
            # and outside the 0.62 to 1.55 kcal/mol that published functionals span
            assert not 0.62 <= m[reading]["margin_kcal"] <= 1.55, (name, reading)
    usable = {k: m for k, m in models.items() if k != "egret-1e"}
    errors = [abs(m[r]["margin_error_kcal"]) for m in usable.values()
              for r in ("on_the_dft_geometries", "relaxed_in_the_potential")]
    assert 1.2 < min(errors) < 1.3 and 3.9 < max(errors) < 4.0
    # they do not agree on the sign: one puts all-trans below the helix
    signs = {m["on_the_dft_geometries"]["margin_kcal"] > 0 for m in usable.values()}
    assert signs == {True, False}
    # the antipolar/polar pair, 1.1 meV apart in DFT, is 10 to 220 meV apart on the DFT geometries
    pair = [abs(m["on_the_dft_geometries"]["anti_minus_beta_meV"]) for m in usable.values()]
    assert 10 < min(pair) and max(pair) < 220


def test_relaxing_in_the_potential_moves_the_volume_by_up_to_a_fifth():
    models = json.loads(RESULT.read_text())["models"]
    moved = 0
    for name in ("medium-mpa-0", "medium-0b3", "medium", "egret-1", "egret-1t"):
        cells = models[name]["relaxed_in_the_potential"]["cells"]
        assert all(cells[c]["converged"] for c in cells), name
        change = max(abs(cells[c]["volume_change_percent"]) for c in ("alpha", "beta"))
        assert change < 21
        moved += change > 7.5
    assert moved == 3
