"""examples/dft_cell_geometry.py: backbone angles, torsions and bond lengths of a periodic cell by distance over images."""
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
