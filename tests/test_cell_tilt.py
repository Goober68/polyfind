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


def test_the_tilt_is_the_third_vector_in_units_of_the_first_two():
    c = cell_of([0.0, 0.0, 35.0]).cell.array
    assert CT.tilt(c) == pytest.approx([0.0, 0.0]) and CT.shear_along_chains(c) == pytest.approx([0.0, 0.0])
    c = cell_of(-0.5 * np.array([5.0, 0.0, 0.0]) - 0.4 * np.array([0.5, 9.0, 0.0]) + np.array([0.0, 0.0, 35.0])).cell.array
    assert CT.tilt(c) == pytest.approx([-0.5, -0.4])
    # a tilted axis puts the neighbours along the first two vectors off along the chain: that is the shear along the chains
    axis = c[2] / np.linalg.norm(c[2])
    assert CT.shear_along_chains(c) == pytest.approx([c[0] @ axis / 2.56, c[1] @ axis / 2.56]) and all(CT.shear_along_chains(c) < 0.0)


def test_a_turned_axis_keeps_the_repeat_and_the_fractional_coordinates_and_drops_the_old_energy():
    a = cell_of([0.2, -0.1, 35.0])
    b = CT.with_axis(a, a.cell.array[2] + np.array([-2.0, -3.0, 0.0]), "turned", "the chain axis turned")
    assert np.linalg.norm(b.cell.array[2]) == pytest.approx(np.linalg.norm(a.cell.array[2]))
    assert b.cell.array[:2] == pytest.approx(a.cell.array[:2]) and b.get_scaled_positions() == pytest.approx(a.get_scaled_positions())
    assert b.info["frame"] == "turned" and b.info["source"] == "built; the chain axis turned" and "potential_eV_per_monomer" not in b.info and "note" not in b.info
    assert a.info["frame"] == "start" and "potential_eV_per_monomer" in a.info  # the cell it was made from is untouched
