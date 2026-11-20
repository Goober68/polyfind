"""examples/cell_tilt.py: starts for a cell whose chain axis tilts as it relaxes."""
import importlib.util
from pathlib import Path

import numpy as np
import pytest

ase = pytest.importorskip("ase")
ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("cell_tilt", ROOT / "examples" / "cell_tilt.py")
CT = importlib.util.module_from_spec(spec)
spec.loader.exec_module(CT)


def cell_of(third):
    from ase import Atoms

    a = Atoms("C2", positions=[[0.5, 0.5, 1.0], [1.0, 2.0, 20.0]], cell=[[5.0, 0.0, 0.0], [0.5, 9.0, 0.0], third], pbc=True)
    a.info.update(frame="start", source="built", note="relaxed", potential_eV_per_monomer=-1.0)
    return a


def test_the_shear_along_the_chains_does_not_depend_on_how_the_cell_lies_in_space():
    upright = cell_of([0.0, 0.0, 35.0]).cell.array
    assert CT.shear_along_chains(upright) == pytest.approx([0.0, 0.0])
    tilted = cell_of([-2.5, -3.6, 35.0]).cell.array
    shear = CT.shear_along_chains(tilted)
    assert all(shear < 0.0)  # a chain axis leaning back along a and b puts the neighbours along a and b behind along the chain
    turn = np.array([[np.cos(0.4), -np.sin(0.4), 0.0], [np.sin(0.4), np.cos(0.4), 0.0], [0.0, 0.0, 1.0]]) @ np.array([[1.0, 0.0, 0.0], [0.0, np.cos(0.7), -np.sin(0.7)], [0.0, np.sin(0.7), np.cos(0.7)]])
    assert CT.shear_along_chains(tilted @ turn.T) == pytest.approx(shear)  # the same cell turned in space


def test_the_axis_found_for_a_shear_gives_that_shear():
    cell = cell_of([0.3, -0.2, 35.0]).cell.array
    for shear in ([0.0, 0.0], [-0.5, -1.1], [0.4, -0.9]):
        axis = CT.axis_for(cell, np.array(shear))
        assert np.linalg.norm(axis) == pytest.approx(1.0) and axis @ cell[2] > 0.0
        assert CT.shear_along_chains(np.array([cell[0], cell[1], 35.0 * axis])) == pytest.approx(shear)
    with pytest.raises(ValueError):
        CT.axis_for(cell, np.array([3.0, 0.0]))  # 7.7 A along the chain for a neighbour 5 A away


def test_a_turned_axis_keeps_the_repeat_and_the_fractional_coordinates_and_drops_the_old_energy():
    a = cell_of([0.2, -0.1, 35.0])
    b = CT.with_shear(a, np.array([-0.5, -1.1]), "turned", "the chain axis turned")
    assert CT.shear_along_chains(b.cell.array) == pytest.approx([-0.5, -1.1])
    assert np.linalg.norm(b.cell.array[2]) == pytest.approx(np.linalg.norm(a.cell.array[2]))
    assert b.cell.array[:2] == pytest.approx(a.cell.array[:2]) and b.get_scaled_positions() == pytest.approx(a.get_scaled_positions())
    assert b.info["frame"] == "turned" and b.info["source"] == "built; the chain axis turned" and "potential_eV_per_monomer" not in b.info and "note" not in b.info
    assert a.info["frame"] == "start" and "potential_eV_per_monomer" in a.info  # the cell it was made from is untouched
