"""examples/dft_cell_geometry.py: backbone angles, torsions, bond lengths and pendant-group tilts of a periodic cell by distance over images."""
import importlib.util
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("dft_cell_geometry", ROOT / "examples" / "dft_cell_geometry.py")
DG = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(DG)


def test_a_polyfind_cell_reads_back_its_built_geometry():
    """The rotated-setting beta cell of 5.17 was built at Dipole's PBE-D3 bond lengths with the
    backbone angle refined to 116.0 deg: the reader finds exactly four C-C, four C-H and four C-F
    bonds (one repeat along c joins to its own image), both backbone angles, and all-trans torsions."""
    cells = DG.read_extxyz(str(ROOT / "docs" / "predictions" / "basins_pvdf-dft-valence-belt2_beta_basin1.xyz"))
    (name, c), = cells.items()
    r = DG.analyse(c["cell_A"], c["positions_A"], c["symbols"])
    assert r["n_atoms"] == 12
    assert {k: v["n"] for k, v in r["bond_lengths_A"].items()} == {"C-C": 4, "C-H": 4, "C-F": 4}
    assert abs(r["bond_lengths_A"]["C-C"]["mean"] - 1.5331) < 2e-4 and abs(r["bond_lengths_A"]["C-F"]["mean"] - 1.3863) < 2e-4
    assert len(r["backbone_angles"]) == 4 and set(r["backbone_angle_by_substituents_deg"]) == {"FF", "HH"}
    assert abs(r["backbone_angle_mean_deg"] - 115.97) < 0.02
    assert r["torsion_abs_values_deg"] == [180.0]
    assert abs(r["cell"]["gamma_deg"] - 88.47) < 0.01


def test_a_bond_is_counted_once_across_images():
    """A one-atom chain along c bonded to its own image: one bond, not two."""
    cell = [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 1.5]]
    bl = DG.bonds(cell, [[0.0, 0.0, 0.0]], ["C"])
    assert len(bl) == 1 and abs(bl[0][3] - 1.5) < 1e-12


def _group(rock=0.0, wag=0.0, twist=0.0):
    """One CF2 carbon at the origin with its two backbone neighbours in the xy plane (bisector +y, chain
    along x) and two fluorines placed as the chain builder would, then rotated as a rigid pair about the
    frame's normal (rock), the chain direction (wag) or the bisector (twist), in a 12 A box."""
    theta, a = np.radians(112.0), np.radians(108.0)
    c1 = 1.53 * np.array([-np.sin(theta / 2), -np.cos(theta / 2), 0.0])
    c2 = 1.53 * np.array([np.sin(theta / 2), -np.cos(theta / 2), 0.0])
    u, v = np.array([0.0, 1.0, 0.0]), np.array([1.0, 0.0, 0.0])
    w = np.cross(u, v)
    f = [1.38 * (np.cos(a / 2) * u + s * np.sin(a / 2) * w) for s in (1.0, -1.0)]

    def rot(axis, deg):
        k, th = axis / np.linalg.norm(axis), np.radians(deg)
        K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
        return np.eye(3) + np.sin(th) * K + (1 - np.cos(th)) * (K @ K)

    R = rot(u, twist) @ rot(v, wag) @ rot(w, rock)
    f = [R @ x for x in f]
    pos = np.array([[0.0, 0.0, 0.0], c1, c2, f[0], f[1]]) + 6.0
    cell = np.diag([12.0, 12.0, 12.0])
    return DG.analyse(cell, pos, ["C", "C", "C", "F", "F"])


def test_the_builders_symmetric_pendant_placement_reads_as_zero_tilt():
    r = _group()
    (g,) = r["pendant_groups"]
    assert g["substituents"] == "FF" and abs(g["xcx_deg"] - 108.0) < 1e-9
    assert max(g["rock_deg"], g["wag_deg"], g["twist_deg"]) < 1e-9
    assert r["pendant_ccx_asymmetry_max_deg"] < 1e-9 and abs(r["backbone_angle_mean_deg"] - 112.0) < 1e-9


def test_each_pendant_tilt_is_read_back_and_shows_in_the_right_c_c_x_asymmetry():
    """Rock moves one X against the two neighbours and leaves the two X alike against either; wag does the
    opposite; twist shows in both.  The X-C-X angle is a rigid-pair invariant and stays 108 deg throughout."""
    for kw, expect in (({"rock": 5.0}, (5.0, 0.0, 0.0)), ({"wag": 4.0}, (0.0, 4.0, 0.0)), ({"twist": 3.0}, (0.0, 0.0, 3.0))):
        (g,) = _group(**kw)["pendant_groups"]
        got = (g["rock_deg"], g["wag_deg"], g["twist_deg"])
        assert np.allclose(got, expect, atol=1e-7), (kw, got)
        assert abs(g["xcx_deg"] - 108.0) < 1e-9
        if "rock" in kw:
            assert g["ccx_same_x_asymmetry_deg"] > 4.0 and g["ccx_same_neighbour_asymmetry_deg"] < 1e-9
        if "wag" in kw:
            assert g["ccx_same_neighbour_asymmetry_deg"] > 3.0 and g["ccx_same_x_asymmetry_deg"] < 1e-9
        if "twist" in kw:
            assert g["ccx_same_x_asymmetry_deg"] > 1.0 and g["ccx_same_neighbour_asymmetry_deg"] > 1.0
    (g,) = _group(rock=5.0, wag=4.0, twist=3.0)["pendant_groups"]
    assert np.allclose((g["rock_deg"], g["wag_deg"], g["twist_deg"]), (5.0, 4.0, 3.0), atol=0.3)


def test_a_polyfind_cell_has_no_pendant_tilt():
    """The chain builder's equal-split, in-plane placement is what the frame angles are measured from."""
    cells = DG.read_extxyz(str(ROOT / "docs" / "predictions" / "basins_pvdf-dft-valence-belt2_beta_basin1.xyz"))
    (name, c), = cells.items()
    r = DG.analyse(c["cell_A"], c["positions_A"], c["symbols"])
    assert len(r["pendant_groups"]) == 4
    assert max(r["pendant_rock_max_deg"], r["pendant_wag_max_deg"], r["pendant_twist_max_deg"]) < 1e-4  # the file holds 6 decimals
    assert set(r["pendant_xcx_by_substituents_deg"]) == {"FF", "HH"}
