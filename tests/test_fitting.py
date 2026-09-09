"""The parameter fit: what it is allowed to change, what it is scored on, what it holds out."""
import numpy as np
import pytest

from polyfind import fitting as F
from polyfind import pack as pack_mod
from polyfind import refine as refine_mod
from polyfind.forcefield import DEFAULT_TORSION, PRESETS, SimpleFF
from polyfind.pack import periodic_chain
from polyfind.polymers import PVDF, THREE_STATE, UFF_LJ

T, GP, GM = 0, 1, 2


def _alpha_chain():
    return periodic_chain(PVDF, [T, GP, T, GM], THREE_STATE)


# ------------------------------------------------------------------- parameters
def test_vector_round_trip_and_illustrative_start():
    x = np.array([2.5, 0.8, 1.05, 0.95, 0.4])
    assert np.allclose(F.vector_from_parameters(F.parameters_from_vector(x)), x)
    # X0 is exactly the illustrative potential: same SimpleFF energies as a bare one
    p0 = F.parameters_from_vector(F.X0)
    assert p0.eps_r == 1.0 and p0.charge_scale == 1.0 and p0.torsion == DEFAULT_TORSION
    assert p0.lj_dict == {"H": (2.886, 0.044), "F": (3.364, 0.050)}  # the UFF values themselves
    from polyfind.chain import build_chain

    s = build_chain(PVDF, np.full(8, 180.0))
    assert p0.simple_ff().energy(s) == SimpleFF().energy(s)


def test_applied_patches_pack_and_refine_and_puts_them_back():
    before = (pack_mod.lj_params, pack_mod.CrystalPacker, refine_mod.CrystalPacker)
    params = F.FFParameters(torsion=(0.5, 0.0, 1.0), eps_r=3.0, charge_scale=0.5,
                            lj=(("F", 4.0, 0.05),))
    with params.applied():
        assert (pack_mod.lj_params, pack_mod.CrystalPacker) != before[:2]
        packer = pack_mod.CrystalPacker(_alpha_chain(), n_chains=2)
        assert packer.eps_r == 3.0 and packer.torsion == (0.5, 0.0, 1.0)
        # charges scaled on the packer, not on the chain the caller passed in
        chain = _alpha_chain()
        assert np.allclose(packer._q_cell[: chain.n_atoms], chain.charges * 0.5)
        x, _ = pack_mod.lj_params(["C", "F", "H"])
        assert x[1] == 4.0 and x[0] == UFF_LJ["C"][0]
    assert (pack_mod.lj_params, pack_mod.CrystalPacker, refine_mod.CrystalPacker) == before


def test_applied_restores_even_when_the_body_raises():
    before = (pack_mod.lj_params, pack_mod.CrystalPacker, refine_mod.CrystalPacker)
    with pytest.raises(RuntimeError):
        with F.FFParameters(eps_r=9.0).applied():
            raise RuntimeError("boom")
    assert (pack_mod.lj_params, pack_mod.CrystalPacker, refine_mod.CrystalPacker) == before


def test_an_unpatched_packer_is_unchanged_by_a_fit():
    """The fit must not leak into anyone else's potential: same chain, same numbers."""
    chain = _alpha_chain()
    p = np.array([[5.09, 9.14, 90.0, 90.0, 270.0, 2.24, 1]])
    e0 = pack_mod.CrystalPacker(chain, n_chains=2).energy(p)[0]
    with F.FFParameters(eps_r=4.0, charge_scale=0.5, lj=(("F", 4.2, 0.05),)).applied():
        e1 = pack_mod.CrystalPacker(chain, n_chains=2).energy(p)[0]
    e2 = pack_mod.CrystalPacker(chain, n_chains=2).energy(p)[0]
    assert e1 != e0
    assert e2 == e0
    assert UFF_LJ["F"] == (3.364, 0.050)


# ------------------------------------------------------------ the antipolar cell
@pytest.mark.parametrize("phi,dz", [(0.0, 0.0), (37.0, 1.3), (90.0, 2.28), (211.0, 4.1)])
def test_flip_with_equal_setting_angles_is_exactly_unpolarized(phi, dz):
    """The premise of the alpha discriminator, measured rather than assumed.

    Flipping chain 2 reverses its axial dipole and equal setting angles put the two
    transverse dipoles head to head, so the cell dipole vanishes identically -- for
    every cell, not only at a minimum.
    """
    packer = pack_mod.CrystalPacker(_alpha_chain(), n_chains=2)
    p = np.array([[5.1, 9.2, 90.0, phi, phi, dz, 1]])
    assert np.abs(packer.polarization(p)).max() < 1e-12
    # and the polar comparison is a real one: unequal setting angles are not zero
    q = np.array([[5.1, 9.2, 90.0, phi, phi + 180.0, dz, 1]])
    assert np.linalg.norm(packer.polarization(q)) > 1e-3


def test_antipolar_cell_returns_a_zero_polarization_minimum():
    chain = _alpha_chain()
    packer = pack_mod.CrystalPacker(chain, n_chains=2)
    start = packer.result(np.array([5.09, 9.14, 90.0, 90.0, 270.0, 2.24, 1]))
    params, e_per_mon = F.antipolar_cell(packer, start, maxfev=120)
    assert params[3] == params[4] and params[6] == 1  # equal setting angles, flipped
    assert np.abs(packer.polarization(params[None])).max() < 1e-12
    assert e_per_mon <= packer.energy_per_monomer(params[None])[0] + 1e-9
    assert np.isfinite(e_per_mon)


# -------------------------------------------------------------------- objective
def _prediction(key, a, b, c, rho, e, p, gap=None):
    return F.Prediction(key, a, b, c, rho, e, p, gap)


def _exact_predictions():
    """Predictions that hit every fitted target exactly."""
    out = {}
    for case in F.FITTED_CASES:
        ea, eb, ec, erho = case.cell
        ea, eb = sorted((ea, eb))
        out[case.key] = _prediction(case.key, ea, eb, ec, erho, 0.0,
                                    0.0 if case.antipolar else (case.polarization or 0.0),
                                    0.0 if case.antipolar else None)
    out["beta"].polarization = F.P_BETA_EXPERIMENT
    out["alpha"].energy_per_monomer = F.ENERGY_GAP_TARGET  # beta sits at 0.0
    return out


def test_a_perfect_prediction_scores_zero():
    res = F.all_residuals(_exact_predictions())
    assert sum(r.contribution for r in res) == pytest.approx(0.0, abs=1e-12)
    assert {r.name for r in res} >= {"pe a", "beta |P|", "alpha antipolar-polar",
                                     "alpha-beta energy", "alpha |P| / beta |P|"}


def test_weights_are_the_documented_ones():
    res = {r.name: r for r in F.all_residuals(_exact_predictions())}
    assert res["pe a"].weight == 1.0 and res["pe a"].tol == pytest.approx(0.02 * 4.95)
    assert res["pe rho"].weight == F.DENSITY_WEIGHT == 0.25  # not an independent observation
    assert res["alpha antipolar-polar"].weight == 2.0
    assert res["alpha-beta energy"].weight == 2.0 and res["alpha-beta energy"].tol == 0.5
    assert res["beta |P|"].tol == 0.02


def test_the_antipolar_residual_is_one_sided():
    """Once the antipolar cell wins, winning by more is not worth anything."""
    preds = _exact_predictions()
    preds["alpha"].polar_gap = -2.0
    assert {r.name: r.contribution for r in F.all_residuals(preds)}["alpha antipolar-polar"] == 0.0
    preds["alpha"].polar_gap = +0.2
    worse = {r.name: r.contribution for r in F.all_residuals(preds)}["alpha antipolar-polar"]
    assert worse == pytest.approx(2.0 * (0.2 / 0.1) ** 2)


def test_the_alpha_polarization_residual_cannot_be_gamed_by_shrinking_the_charges():
    """It is alpha's polarization *relative to beta's*, so a uniform rescale is free."""
    preds = _exact_predictions()
    preds["alpha"].polarization, preds["beta"].polarization = 0.078, 0.140
    r1 = {r.name: r.value for r in F.all_residuals(preds)}["alpha |P| / beta |P|"]
    preds["alpha"].polarization, preds["beta"].polarization = 0.039, 0.070  # every charge halved
    r2 = {r.name: r.value for r in F.all_residuals(preds)}["alpha |P| / beta |P|"]
    assert r1 == pytest.approx(r2)
    assert r1 == pytest.approx(0.078 / 0.140)


def test_gamma_is_held_out_and_pvdc_is_not_a_target():
    keys = {c.key for c in F.FITTED_CASES}
    assert keys == {"pe", "beta", "alpha"}
    assert {c.key for c in F.HELD_OUT_CASES} == {"gamma"}
    assert not keys & {c.key for c in F.HELD_OUT_CASES}
    assert "pvdc" not in {c.polymer.name for c in F.ALL_CASES}  # see the module docstring
    preds = _exact_predictions()
    preds["gamma"] = _prediction("gamma", 1.0, 1.0, 1.0, 1.0, 0.0, 0.0)  # nonsense on purpose
    assert sum(r.contribution for r in F.all_residuals(preds)) == pytest.approx(0.0, abs=1e-12)
    assert all(not r.name.startswith("gamma") for r in F.all_residuals(preds))


def test_experimental_targets_come_from_the_repository():
    from polyfind.pipeline import EXPERIMENTAL_CELLS

    for case in F.ALL_CASES:
        assert case.cell == EXPERIMENTAL_CELLS[case.polymer.name][case.label]


# ---------------------------------------------------------------------- presets
def test_the_fitted_preset_is_available_and_changes_nothing_by_default():
    assert "pvdf-crystal-fit" in PRESETS
    fitted = SimpleFF.from_preset("pvdf-crystal-fit")
    default = SimpleFF()
    assert (default.torsion, default.eps_r, default.charge_scale, default.lj) == \
           (DEFAULT_TORSION, 1.0, 1.0, None)
    assert (fitted.torsion, fitted.eps_r, fitted.charge_scale) != \
           (default.torsion, default.eps_r, default.charge_scale) or fitted.lj != default.lj
    assert F.FITTED_PVDF.preset_kwargs() == PRESETS["pvdf-crystal-fit"]
    assert UFF_LJ["F"] == (3.364, 0.050)  # selecting a preset never edits the global table
