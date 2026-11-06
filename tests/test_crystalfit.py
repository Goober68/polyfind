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


def test_pair_term_reaches_the_packer_and_the_force_field():
    """A pair-specific H-F term overrides the geometric mean in both the lattice kernel and SimpleFF."""
    from polyfind import pack as P
    from polyfind.polymers import PVDF
    x = np.concatenate([C.X0, [1.1, 2.0]])
    params = C.parameters(x)
    assert params.lj_pairs and params.lj_pairs[0][:2] == ("F", "H")
    table = {e: (xi, di) for e, xi, di in params.lj}
    x_hf = 1.1 * np.sqrt(table["H"][0] * table["F"][0])
    d_hf = 2.0 * np.sqrt(table["H"][1] * table["F"][1])
    assert np.isclose(params.lj_pairs[0][2], x_hf) and np.isclose(params.lj_pairs[0][3], d_hf)
    from polyfind.polymers import THREE_STATE
    chain = P.periodic_chain(PVDF, [0, 0], THREE_STATE)
    el = list(chain.elements)
    i_h, i_f, i_c = el.index("H"), el.index("F"), el.index("C")
    with params.applied():
        pk = P.CrystalPacker(chain, n_chains=2, coulomb="ewald")
        A = np.asarray(pk._A_nn); B = np.asarray(pk._B_nn)
        # v = A r^-12 - B r^-6 with A = D x^12, B = 2 D x^6  ->  x^6 = 2A/B, D = B^2/(4A)
        x6 = 2.0 * A / B; D = B * B / (4.0 * A)
        assert np.isclose(x6[i_h, i_f] ** (1 / 6), x_hf) and np.isclose(D[i_h, i_f], d_hf)
        assert np.isclose(x6[i_c, i_f] ** (1 / 6), np.sqrt(table["C"][0] * table["F"][0]))  # untouched pair
    with C.parameters(C.X0).applied():
        pk0 = P.CrystalPacker(chain, n_chains=2, coulomb="ewald")
        # without the pair term the H-F entry is the geometric-mean mixing of the valence values
        assert np.isclose(np.asarray(pk0._A_nn)[i_h, i_f], np.sqrt(table["H"][1] * table["F"][1]) * (table["H"][0] * table["F"][0]) ** 6)
    # SimpleFF's own pair terms
    ff = params.simple_ff()
    assert ff.lj_pair_table() == {("F", "H"): (float(x_hf), float(d_hf))}
    iu, ju, lj_x, lj_d, qq = ff._pair_terms(el, np.full((len(el), len(el)), 9), np.zeros(len(el)))
    m = ((np.array(el)[iu] == "H") & (np.array(el)[ju] == "F")) | ((np.array(el)[iu] == "F") & (np.array(el)[ju] == "H"))
    assert m.any() and np.allclose(lj_x[m], x_hf) and np.allclose(lj_d[m], d_hf)
    assert ff._param_key() != C.parameters(C.X0).simple_ff()._param_key()


def test_identity_pair_multipliers_change_nothing():
    a = C.parameters(np.concatenate([C.X0, [1.0, 1.0]]))
    b = C.parameters(C.X0)
    table = {e: (xi, di) for e, xi, di in b.lj}
    assert np.isclose(a.lj_pairs[0][2], np.sqrt(table["H"][0] * table["F"][0]))
    assert np.isclose(a.lj_pairs[0][3], np.sqrt(table["H"][1] * table["F"][1]))
    assert C.ridge(np.concatenate([C.X0, [1.0, 1.0]])) == 0.0 == C.ridge(C.X0)


def test_predict_with_induced_dipoles_runs_on_beta():
    from polyfind.ewald import EwaldSpec
    from polyfind.polarizability import Polarizable
    pred = F.predict(F.FITTED_CASES[1], F.FITTED_VALENCE, coulomb="ewald", ewald=EwaldSpec(),
                     polarizable=Polarizable(), lj_cutoff="force", refine=False, n_refine=1)
    assert 4.0 < pred.a < 5.5 and 0.05 < pred.polarization < 0.4
    assert pred.rigid.polarization_magnitude == pred.polarization


def test_registered_stage2_preset_matches_the_saved_fit():
    with open("docs/predictions/crystal_stage_v2.json") as fh:
        saved = json.load(fh)
    assert saved["variables"] == list(C.VARIABLES2) == list(F.CRYSTAL_STAGE2_VARIABLES)
    np.testing.assert_allclose(F.CRYSTAL2_FITTED_X, saved["x"], rtol=0, atol=1e-12)
    p = F.FITTED_VALENCE_CRYSTAL2
    assert [list(t) for t in p.lj] == [list(t) for t in saved["parameters"]["lj"]]
    assert [list(t) for t in p.lj_pairs] == [list(t) for t in saved["parameters"]["lj_pairs"]]
    assert [list(t) for t in p.charge_increments] == [list(t) for t in saved["parameters"]["charge_increments"]]
    assert "pvdf-dft-valence-crystal2" in PRESETS and F.PRESET_PARAMETERS["pvdf-dft-valence-crystal2"] is p
    ff = SimpleFF.from_preset("pvdf-dft-valence-crystal2")
    assert ff.lj_pair_table() == {("F", "H"): tuple(saved["parameters"]["lj_pairs"][0][2:])}


def test_registered_stage4_preset_matches_the_saved_fit():
    with open("docs/predictions/crystal_stage_v4.json") as fh:
        saved = json.load(fh)
    assert saved["variables"] == list(C.VARIABLES4) == list(F.CRYSTAL_STAGE4_VARIABLES)
    np.testing.assert_allclose(F.CRYSTAL4_FITTED_X, saved["x"], rtol=0, atol=1e-12)
    p = F.FITTED_VALENCE_BELT
    assert [list(t) for t in p.lj] == [list(t) for t in saved["parameters"]["lj"]]
    assert [list(t) for t in p.lj_pairs] == [list(t) for t in saved["parameters"]["lj_pairs"]]
    assert [list(t) for t in p.charge_increments] == [list(t) for t in saved["parameters"]["charge_increments"]]
    assert p.charge_damping == saved["parameters"]["charge_damping"]
    assert "pvdf-dft-valence-belt" in PRESETS and F.PRESET_PARAMETERS["pvdf-dft-valence-belt"] is p
    ff = SimpleFF.from_preset("pvdf-dft-valence-belt")
    table = ff.lj_pair_table()
    assert table[("F", "H")] == tuple(saved["parameters"]["lj_pairs"][0][2:])
    # the wall is X-H only: the site sees no carbon, fluorine or other site
    x_depth = dict((e, d) for e, _, d in p.lj)["X"]
    assert x_depth == 10 ** saved["x"][10] and x_depth > 0
    for pair in (("C", "X"), ("F", "X"), ("X", "X")):
        assert table[pair][1] == 0.0


def test_polar_slope_removes_the_cell_orientation():
    from types import SimpleNamespace
    step = 2e-3
    up = {(0, 1): SimpleNamespace(polarization=np.array([0.1890, 0.0, 0.0])),
          (0, -1): SimpleNamespace(polarization=np.array([0.1960, 0.0, 0.0]))}
    down = {k: SimpleNamespace(polarization=-v.polarization) for k, v in up.items()}
    s_up, s_down = C.polar_slope(up, step), C.polar_slope(down, step)
    assert s_up == s_down < 0
    np.testing.assert_allclose(s_up, (0.1890 - 0.1960) / (2 * step))


def test_slope_anchor_is_dipoles_three_point_berry_pair():
    with open("docs/predictions/e1_dipole_existing.json") as fh:
        rec = json.load(fh)["systems"]["beta_qe_berry_axial_v1"]
    cells = rec["strained"]["xx (polar)"]
    centered = (cells["+0.010"]["P_C_m2"][0] - cells["-0.010"]["P_C_m2"][0]) / 0.02
    np.testing.assert_allclose(rec["dP_d_eta_polar_C_m2"]["centered"], centered, atol=1e-6)
    np.testing.assert_allclose(C.DFT_ANCHORS["beta_dP_deta_polar"], centered, atol=1e-4)
    assert C.STAGE6["slope"] and C.STAGE6["stiffness"] and C.STAGE6["polarizable"]
