"""The crystal stage: its parameter map, the registered preset and the saved result agree."""
import json

import numpy as np

from polyfind import crystalfit as C
from polyfind import fitting as F
from polyfind.forcefield import PRESETS, SimpleFF


def test_identity_variables_reproduce_the_valence_preset():
    p = C.parameters(C.X0)
    assert p.lj == F.FITTED_VALENCE.lj
    assert p.charge_increments == F.FITTED_VALENCE.charge_increments
    assert p.torsion == F.FITTED_VALENCE.torsion and p.bond_terms == F.FITTED_VALENCE.bond_terms


def test_variables_move_only_what_they_name():
    x = np.array([0.95, 1.05, 1.0, 0.5, 2.0, 1.0, -0.02, 0.03])
    p = C.parameters(x)
    base = dict((e, (xi, di)) for e, xi, di in F.FITTED_VALENCE.lj)
    got = dict((e, (xi, di)) for e, xi, di in p.lj)
    assert np.isclose(got["C"][0], 0.95 * base["C"][0]) and np.isclose(got["C"][1], 0.5 * base["C"][1])
    assert np.isclose(got["H"][0], 1.05 * base["H"][0]) and np.isclose(got["H"][1], 2.0 * base["H"][1])
    assert got["F"] == base["F"]
    assert got["Cl"] == base["Cl"]
    inc = {(a, b): d for a, b, d in p.charge_increments}
    base_inc = {(a, b): d for a, b, d in F.FITTED_VALENCE.charge_increments}
    assert np.isclose(inc[("C", "H")], base_inc[("C", "H")] - 0.02)
    assert np.isclose(inc[("C", "F")], base_inc[("C", "F")] + 0.03)
    assert inc[("C", "Cl")] == base_inc[("C", "Cl")]
    assert p.torsion == F.FITTED_VALENCE.torsion


def test_registered_preset_matches_the_saved_fit():
    with open("docs/predictions/crystal_stage_v1.json") as fh:
        saved = json.load(fh)
    assert saved["variables"] == list(C.VARIABLES) == list(F.CRYSTAL_STAGE_VARIABLES)
    np.testing.assert_allclose(F.CRYSTAL_FITTED_X, saved["x"], rtol=0, atol=1e-12)
    p = F.FITTED_VALENCE_CRYSTAL
    assert [list(t) for t in p.lj] == [list(t) for t in saved["parameters"]["lj"]]
    assert [list(t) for t in p.charge_increments] == [list(t) for t in saved["parameters"]["charge_increments"]]
    assert "pvdf-dft-valence-crystal" in PRESETS
    assert SimpleFF.from_preset("pvdf-dft-valence-crystal") is not None
    assert F.PRESET_PARAMETERS["pvdf-dft-valence-crystal"] is p
