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


def test_the_second_fit_is_no_better_and_the_two_differ_by_the_resolution_needed():
    """Thirty more epochs with energy and stress weighted up: stress better, energies worse."""
    a, b = _result("pvdf_ft_0b3"), _result("pvdf_ft_0b3_b")
    assert b["held_out"]["stress_rmse_GPa"] == pytest.approx(0.27, abs=0.01)
    assert b["held_out"]["stress_rmse_GPa"] < a["held_out"]["stress_rmse_GPa"]
    assert b["held_out"]["energy_rmse_meV_per_atom"] == pytest.approx(0.55, abs=0.01)
    assert b["held_out"]["energy_rmse_meV_per_atom"] > a["held_out"]["energy_rmse_meV_per_atom"]
    m = b["margins"]
    assert m["on_the_dft_geometries"]["margin_kcal"] == pytest.approx(0.617, abs=0.001)
    assert m["relaxed_in_the_potential"]["margin_kcal"] == pytest.approx(0.566, abs=0.001)
    assert all(abs(m[r]["margin_error_kcal"]) > NEEDED_KCAL for r in ("on_the_dft_geometries", "relaxed_in_the_potential"))
    e = b["beta_elastic_GPa"]
    assert e["ours"]["C_11"] / e["dft"]["C_11"] == pytest.approx(1.37, abs=0.01)
    assert e["ours"]["C_22"] / e["dft"]["C_22"] == pytest.approx(0.78, abs=0.01)
    gaps_a = a["dft_on_the_potentials_cells"]["gaps"]
    gaps_b = b["dft_on_the_potentials_cells"]["gaps"]
    assert gaps_b["alpha_minus_beta_meV"]["difference_kcal"] == pytest.approx(0.229, abs=0.002)
    assert gaps_b["anti_minus_beta_meV"]["difference_kcal"] == pytest.approx(-0.246, abs=0.002)
    assert gaps_b["turned_minus_beta_meV"]["difference_kcal"] == pytest.approx(0.105, abs=0.002)
    # fit against fit, on the same test: as large as the 0.1 the second opinion has to resolve
    spread = [abs(gaps_a[k]["difference_kcal"] - gaps_b[k]["difference_kcal"]) for k in gaps_a]
    assert min(spread) > NEEDED_KCAL
    # every gap of either fit is within 0.25 of DFT: the pretrained potentials were off by 1.2 to 3.9
    errors = [abs(g[k]["difference_kcal"]) for g in (gaps_a, gaps_b) for k in g]
    errors += [abs(r["margins"][x]["margin_error_kcal"]) for r in (a, b) for x in ("on_the_dft_geometries", "relaxed_in_the_potential")]
    assert max(errors) < 0.25 and min(errors) > 0.01


def test_dft_puts_the_turned_chain_helix_cell_level_with_alpha():
    """A by-product: the helix packing with the chains' dipoles adding is not above alpha in DFT."""
    alpha_minimum = [f[2] for f in _frames() if f[1] == "end"][0]
    turned = {name: _result(name)["dft_on_the_potentials_cells"]["cells"]["turned"]["dft_eV_per_monomer"]
              for name in ("pvdf_ft_0b3", "pvdf_ft_0b3_b")}
    assert 1e3 * (turned["pvdf_ft_0b3"] - alpha_minimum) == pytest.approx(2.0, abs=0.1)
    assert 1e3 * (turned["pvdf_ft_0b3_b"] - alpha_minimum) == pytest.approx(-0.4, abs=0.1)
    probes = (PRED / "mlip" / "pvdf_probe_frames_pbe_d3bj.extxyz").read_text()
    assert probes.count("family=probe") == 32  # two short fits, four long ones, two of those relaxed strictly


def test_the_anchors_cutoffs_hold_for_chlorine_and_the_nitrile():
    """A CFE pair and an acrylonitrile pair again at 100/800 Ry: nothing a margin would see."""
    check = json.loads((PRED / "mlip_cutoff_check.json").read_text())
    assert check["anchor_cutoffs_ry"] == [90, 360] and check["higher_cutoffs_ry"] == [100, 800]
    sources = [f["source"] for f in check["frames"].values()]
    assert sources == ["TT2 VC|VV", "TG2 VC|VV", "TT2 VA|VV", "TG2 VA|VV"]
    for frame in check["frames"].values():
        assert frame["force_change_rms_meV_per_A"] < 2.0
        assert frame["force_change_largest_meV_per_A"] < 10.0
        assert frame["stress_change_largest_GPa"] < 0.02
    changes = [abs(pair["change_kcal_per_monomer"]) for pair in check["pairs"].values()]
    assert len(changes) == 2 and max(changes) < 0.002 < NEEDED_KCAL / 10


LONG_FITS = ("pvdf_s1", "pvdf_s2", "pvdfb_s1", "pvdfb_s2")


def test_the_long_fits_are_the_dft_where_frames_are_near():
    """Four fits of 400 epochs on 130 frames (a rented card): the tests on the trained packings."""
    for name in LONG_FITS:
        r = _result(name)
        assert r["prepare"]["frames"] == 130 and r["prepare"]["train"] == 112 and r["prepare"]["test"] == 18
        held = r["held_out"]
        assert held["energy_rmse_meV_per_atom"] < 0.35 and held["stress_rmse_GPa"] < 0.11
        assert held["force_rmse_meV_per_A"] / held["force_rms_of_the_labels_meV_per_A"] < 0.02
        m = r["margins"]["on_the_dft_geometries"]
        assert abs(m["margin_error_kcal"]) < 0.007
        assert -1.3 < m["anti_minus_beta_meV"] < -0.5  # DFT -1.11: the sign the short fits had wrong
        e = r["beta_elastic_GPa"]
        for key, tolerance in (("C_11", 0.07), ("C_22", 0.07), ("C_33", 0.07)):
            assert e["ours"][key] / e["dft"][key] == pytest.approx(1.0, abs=tolerance), (name, key)
        assert 1.35 < e["ours"]["C_12"] / e["dft"]["C_12"] < 1.46  # the one constant that is not right
        gaps = r["dft_on_the_potentials_cells"]["gaps"]
        assert abs(gaps["anti_minus_beta_meV"]["difference_kcal"]) < 0.015
        cells = r["dft_on_the_potentials_cells"]["cells"]
        assert max(cells[k]["dft_largest_force_eV_per_A"] for k in ("alpha", "beta", "anti")) < 0.015
    alpha = [_result(n)["dft_on_the_potentials_cells"]["gaps"]["alpha_minus_beta_meV"]["difference_kcal"] for n in LONG_FITS]
    assert sorted(abs(v) for v in alpha)[:3] == pytest.approx([0.005, 0.006, 0.007], abs=0.001)
    assert min(alpha) == pytest.approx(-0.083, abs=0.002)  # the fit that let its helix cell shrink by 1%


def test_the_long_fits_are_not_the_dft_on_a_packing_they_were_not_shown():
    """The turned-chain helix cell: every fit too low, and DFT on their cells agrees with itself."""
    gaps = [_result(n)["dft_on_the_potentials_cells"]["gaps"]["turned_minus_beta_meV"] for n in LONG_FITS]
    errors = [g["difference_kcal"] for g in gaps]
    assert all(e < 0 for e in errors)
    assert min(errors) == pytest.approx(-0.210, abs=0.002) and max(errors) == pytest.approx(-0.068, abs=0.002)
    assert sum(abs(e) > NEEDED_KCAL for e in errors) == 3
    by_potential = [g["potential"] for g in gaps]
    by_dft = [g["dft"] for g in gaps]
    # the fits' spread flags the cell, but their mean is as far off as the spread is wide
    assert 0.0230605 * (max(by_potential) - min(by_potential)) == pytest.approx(0.143, abs=0.003)
    assert sum(errors) / 4 == pytest.approx(-0.132, abs=0.003)
    # DFT single points on four different relaxed geometries: the same gap to 1 meV per monomer
    assert max(by_dft) - min(by_dft) < 1.0
    forces = [_result(n)["dft_on_the_potentials_cells"]["cells"]["turned"]["dft_largest_force_eV_per_A"] for n in LONG_FITS]
    assert max(forces) < 0.04
    # and that packing is level with alpha in DFT, to within 1 meV per monomer
    alpha_minimum = [f[2] for f in _frames() if f[1] == "end"][0]
    turned = [1e3 * (_result(n)["dft_on_the_potentials_cells"]["cells"]["turned"]["dft_eV_per_monomer"] - alpha_minimum)
              for n in LONG_FITS]
    assert -0.2 < min(turned) and max(turned) < 1.0


def test_relaxed_strictly_the_potentials_cells_are_dfts():
    for name in ("pvdfb_s1_tight", "pvdfb_s2_tight"):
        r = _result(name)
        assert r["relaxation_fmax_eV_per_A"] == 0.002
        volumes = [abs(c["volume_change_percent"]) for c in r["margins"]["relaxed_in_the_potential"]["cells"].values()]
        assert max(volumes) < 0.8
        closure = r["dft_on_the_potentials_cells"]
        above = [closure["cells"][k]["dft_above_its_own_minimum_meV_per_monomer"] for k in ("alpha", "beta", "anti")]
        assert -0.2 < min(above) and max(above) < 0.5
        assert abs(closure["gaps"]["alpha_minus_beta_meV"]["difference_kcal"]) < 0.01
        assert abs(closure["gaps"]["anti_minus_beta_meV"]["difference_kcal"]) < 0.025
