"""Our own potential, fitted to the PVDF frames, against the DFT anchor.

Reads the stored results only: torch, mace and ase are not dependencies of the package, and
the model files are not in the repository. `examples/mlip_train.py` regenerates the results
from `docs/predictions/mlip/pvdf_frames_pbe_d3bj.extxyz`.
"""

import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
PRED = ROOT / "docs" / "predictions"
FRAMES = PRED / "mlip" / "pvdf_frames_pbe_d3bj.extxyz"
NEEDED_KCAL = 0.1


def _result(name):
    return json.loads((PRED / f"mlip_own_{name}.json").read_text())


def _frames():
    """(atoms, family, DFT energy per monomer in eV) of each stored frame, without ase."""
    lines = FRAMES.read_text().splitlines()
    out, i = [], 0
    while i < len(lines):
        n = int(lines[i])
        head = lines[i + 1]
        family = head.split("family=")[1].split()[0]
        energy = float(head.split("REF_energy=")[1].split()[0])
        monomers = int(head.split("n_monomers=")[1].split()[0])
        out.append((n, family, energy / monomers))
        i += n + 2
    return out


def test_the_frames_are_the_set_the_write_up_describes():
    frames = _frames()
    families = [f[1] for f in frames]
    assert len(frames) == 76
    assert {k: families.count(k) for k in set(families)} == {"end": 4, "strain": 24, "rattle": 18, "md": 30}
    assert {f[0] for f in frames} == {12, 24}
    # the relaxed cells carry the anchor's own gap: alpha below polar beta by 31.46 meV per monomer
    ends = [f[2] for f in frames if f[1] == "end"]
    assert 1e3 * (ends[0] - ends[1]) == pytest.approx(-31.46, abs=0.02)
    assert 1e3 * (ends[2] - ends[1]) == pytest.approx(-1.11, abs=0.02)


def test_the_first_fit_holds_the_margin_but_not_the_antipolar_pair_or_the_lateral_stiffness():
    r = _result("pvdf_ft_0b3")
    assert r["prepare"]["train"] == 64 and r["prepare"]["test"] == 12 and r["prepare"]["test_ends"] == 4
    held = r["held_out"]
    assert held["energy_rmse_meV_per_atom"] == pytest.approx(0.39, abs=0.01)
    assert held["force_rmse_meV_per_A"] / held["force_rms_of_the_labels_meV_per_A"] == pytest.approx(0.027, abs=0.001)
    assert held["stress_rmse_GPa"] == pytest.approx(0.40, abs=0.01)
    m = r["margins"]
    # the margin, inside the 0.1 needed on both readings
    assert m["on_the_dft_geometries"]["margin_kcal"] == pytest.approx(0.671, abs=0.001)
    assert m["relaxed_in_the_potential"]["margin_kcal"] == pytest.approx(0.700, abs=0.001)
    for reading in ("on_the_dft_geometries", "relaxed_in_the_potential"):
        assert abs(m[reading]["margin_error_kcal"]) < NEEDED_KCAL
    # the antipolar pair, 1.1 meV apart in DFT with the antipolar cell lower: the sign is wrong
    assert m["on_the_dft_geometries"]["anti_minus_beta_meV"] == pytest.approx(2.7, abs=0.1)
    assert m["relaxed_in_the_potential"]["anti_minus_beta_meV"] == pytest.approx(5.1, abs=0.1)
    # its own cells are too small across the chains
    volumes = {k: c["volume_change_percent"] for k, c in m["relaxed_in_the_potential"]["cells"].items()}
    assert volumes["alpha"] == pytest.approx(-4.6, abs=0.1) and volumes["beta"] == pytest.approx(-2.1, abs=0.1)
    e = r["beta_elastic_GPa"]
    assert e["ours"]["C_11"] / e["dft"]["C_11"] == pytest.approx(1.73, abs=0.01)
    assert e["ours"]["C_22"] / e["dft"]["C_22"] == pytest.approx(1.05, abs=0.01)
    assert e["ours"]["C_33"] / e["dft"]["C_33"] == pytest.approx(0.93, abs=0.01)


def test_dft_on_the_first_fits_own_cells():
    """DFT single points on the cells the potential relaxed, and on a packing it was not shown."""
    c = _result("pvdf_ft_0b3")["dft_on_the_potentials_cells"]
    cells, gaps = c["cells"], c["gaps"]
    assert set(cells) == {"alpha", "beta", "anti", "turned"}
    # the atoms are where DFT wants them: little force is left
    assert max(v["dft_largest_force_eV_per_A"] for v in cells.values()) < 0.08
    # the cells are not: DFT puts them 2.5 to 5.5 meV per monomer above its own minima
    above = [cells[k]["dft_above_its_own_minimum_meV_per_monomer"] for k in ("alpha", "beta", "anti")]
    assert all(v > 0 for v in above) and min(above) == pytest.approx(2.5, abs=0.1) and max(above) == pytest.approx(5.5, abs=0.1)
    # helix against all-trans on the same geometries: right to 0.02 kcal/mol, also for the unseen packing
    assert abs(gaps["alpha_minus_beta_meV"]["difference_kcal"]) < 0.02
    assert abs(gaps["turned_minus_beta_meV"]["difference_kcal"]) < 0.02
    # two packings of all-trans chains: off by more than the 0.1 needed
    assert gaps["anti_minus_beta_meV"]["difference_kcal"] == pytest.approx(0.17, abs=0.005)
    # DFT puts the turned-chain helix cell 3.5 meV per monomer below the alpha cell beside it
    assert gaps["turned_minus_beta_meV"]["dft"] - gaps["alpha_minus_beta_meV"]["dft"] == pytest.approx(-3.5, abs=0.1)
