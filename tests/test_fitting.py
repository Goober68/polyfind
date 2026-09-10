"""The parameter fit: what it is allowed to change, what it is scored on, what it holds out."""
import os

import numpy as np
import pytest

from polyfind import fitting as F
from polyfind import pack as pack_mod
from polyfind import refine as refine_mod
from polyfind.forcefield import DEFAULT_TORSION, PRESETS, UFF_LJ_EXTRA, SimpleFF
from polyfind.pack import periodic_chain
from polyfind.polymers import PE, PVDF, THREE_STATE, UFF_LJ, get_polymer

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


# ============================================ the fit to first-principles data
def _synthetic_frames():
    """Frames built from the package's own chains, so no external data file is needed.

    Two "systems" of a few conformers each, with the reference energy set to what a
    *known* parameter vector says, so a fit has a target it can reach exactly and the
    machinery can be tested without asserting anything about the real data.
    """
    from polyfind.chain import build_chain
    from polyfind.forcefield import Frame

    rng = np.random.default_rng(11)
    conformers = [  # rigid-geometry chains fold onto themselves at arbitrary torsions;
        [178.0, -176.0, 179.0, 174.0, -178.0, 177.0, 180.0],   # see tests/test_forcefield.py
        [172.0, 62.0, 176.0, -58.0, 168.0, 65.0, -175.0],
        [175.0, 171.0, 178.0, 58.0, -174.0, 176.0, 172.0],
        [-172.0, 178.0, 64.0, 170.0, 176.0, -61.0, 174.0],
    ]
    out = []
    # PVDF and PE, not PVDC: with frozen backbone angles the chlorines of an all-trans
    # PVDC chain are inside each other's covalent radii (that same crowding is the 313
    # kcal/mol of strain polymers.py's PVDC comment records), so bond perception on a
    # *built* PVDC chain finds rings that the real, relaxed molecule does not have.
    for name, polymer in (("pvdf", PVDF), ("pe", PE)):
        for dih in conformers:
            s = build_chain(polymer, dih)
            f = Frame.from_geometry(s.elements, s.coords, system=name, source="conf")
            f.forces = rng.normal(0, 3.0, (f.n_atoms, 3))
            f.energy = 0.0
            out.append(f)
    return out


def test_reference_design_reproduces_the_reference_implementation():
    """The precomputed design and ``SimpleFF`` must be the same function, to rounding.

    The design exists only for speed.  If it drifted from ``energy_frame`` the fit would
    be optimising something the package does not ship, which is the one failure mode that
    would make every number downstream meaningless.
    """
    frames = _synthetic_frames()
    design = F.ReferenceDesign(frames)
    rng = np.random.default_rng(5)
    x = np.clip(F.REF_X0 + 0.4 * rng.normal(0, F.REF_SIGMA), F.REF_BOUNDS[:, 0], F.REF_BOUNDS[:, 1])
    energy, torque = design.evaluate(F.reference_unpack(x))
    p = F.reference_unpack(x)
    names = sorted({**UFF_LJ_EXTRA, **UFF_LJ})
    ff = SimpleFF(
        eps_r=1.0, scale14=float(p["scale14"]),
        lj={e: (float(p["lj_x"][names.index(e)]), float(p["lj_d"][names.index(e)])) for e in names},
        charge_increments=tuple((a, b, float(v)) for (a, b), v in zip(F.REF_INCREMENTS, p["increments"])))
    for i, f in enumerate(frames):
        V = p["torsion"][F.torsion_types(f)]
        assert energy[i] == pytest.approx(ff.energy_frame(f, torsion=V), abs=1e-8)
        sel = design.dih_frame == i
        assert np.allclose(torque[sel], ff.frame_torques(f, torsion=V), atol=1e-5)


def test_permittivity_and_charge_scale_stay_degenerate_on_relaxed_geometries():
    """DESIGN.md 5.7 hoped off-ideal geometries would separate them.  They do not, and
    they cannot: the degeneracy is exact in the functional form, not an accident of the
    ideal-angle chains that fit was run on."""
    design = F.ReferenceDesign(_synthetic_frames())
    d = F.permittivity_degeneracy(design)
    assert d["max_energy_change"] < 1e-9
    assert d["max_torque_change"] < 1e-9


def test_the_split_holds_out_whole_systems_and_keeps_the_anchor():
    frames = _synthetic_frames() + []
    train, test = F.split_systems(frames, anchor="pvdf", stride=2, offset=1)
    assert "pvdf" in train
    assert set(train) & set(test) == set()
    assert set(train) | set(test) == {f.system for f in frames}


def test_pvdf_torsions_are_all_one_type():
    """Every backbone bond of PVDF is a CF2-CH2 torsion, which is why one triple of the
    fitted table is the whole of the shipped preset's ``torsion``."""
    from polyfind.chain import build_chain
    from polyfind.forcefield import Frame

    s = build_chain(PVDF, np.full(7, 175.0))
    f = Frame.from_geometry(s.elements, s.coords)
    assert F.REF_TORSION_TYPES[0] == "CF2-CH2"
    assert set(F.torsion_types(f)) == {0}


def test_reference_x0_is_the_illustrative_potential():
    """The fit starts where the package already is, so "did fitting help" is answerable."""
    p = F.reference_ff_parameters(F.REF_X0)
    assert p.torsion == DEFAULT_TORSION
    assert p.scale14 == 0.5 and p.eps_r == 1.0
    assert dict(p.lj_dict)["F"] == pytest.approx(UFF_LJ["F"])
    inc = {(a, b): v for a, b, v in p.charge_increments}
    assert inc[("C", "H")] == pytest.approx(-0.10) and inc[("C", "F")] == pytest.approx(0.20)
    # and those increments are exactly the charges polymers.py gives PVDF
    from polyfind.chain import build_chain
    from polyfind.forcefield import bci_charges

    s = build_chain(PVDF, np.full(7, 180.0))
    q = bci_charges(s.elements, s.bonds, {("C", "H"): -0.10, ("C", "F"): 0.20})
    assert np.allclose(q[3:27], s.charges[3:27])


def test_charge_increments_reach_the_lattice_unchanged():
    """A fitted charge model has to arrive in ``pack``/``refine``, not only in ``SimpleFF``.

    Written as increments, the illustrative charges must give the packer exactly the
    charges the polymer table does -- otherwise the preset would mean one thing to the
    oligomer side and another to the crystal side.
    """
    chain = _alpha_chain()
    params = F.FFParameters(charge_increments=(("C", "H", -0.10), ("C", "F", 0.20)))
    with params.applied():
        packer = pack_mod.CrystalPacker(chain, n_chains=2)
        assert np.allclose(packer.chain.charges, chain.charges)


def test_charge_increment_bounds_follow_electronegativity():
    """Not a fact about any target here: carbon is more electronegative than hydrogen and
    less than F, Cl, N and O, so H must come out positive and those four negative."""
    lo = dict(zip(F.REF_VARIABLES, F.REF_BOUNDS[:, 0]))
    hi = dict(zip(F.REF_VARIABLES, F.REF_BOUNDS[:, 1]))
    assert hi["q_C-H"] < 0  # C gains charge from H: H positive
    for e in ("F", "Cl", "N", "O"):
        assert lo[f"q_C-{e}"] > 0
    assert lo["q_C-S"] < 0 < hi["q_C-S"]  # C and S are within 0.03 on the Pauling scale


def test_the_dft_preset_is_registered_and_the_defaults_are_untouched():
    assert "pvdf-dft-fit" in PRESETS
    assert F.FITTED_DFT.preset_kwargs() == PRESETS["pvdf-dft-fit"]
    default = SimpleFF()
    assert (default.torsion, default.eps_r, default.charge_scale, default.lj,
            default.charge_increments) == (DEFAULT_TORSION, 1.0, 1.0, None, None)
    fitted = SimpleFF.from_preset("pvdf-dft-fit")
    assert fitted.charge_increments is not None and fitted.torsion != DEFAULT_TORSION
    assert UFF_LJ["F"] == (3.364, 0.050)


@pytest.mark.skipif(not os.environ.get(F.REFERENCE_ENV),
                    reason=f"set ${F.REFERENCE_ENV} to the reference set to run this")
def test_the_reference_set_matches_the_package_chain_for_pvdf():
    """The one check the loader owes: on pure PVDF frames the inferred topology is
    bond-for-bond what ``build_chain`` builds for a five-monomer oligomer."""
    frames = [f for f in F.load_reference() if f.system == "pvdf"]
    assert frames
    for f in frames:
        s = f.to_structure(PVDF)  # raises if the elements or the bonds differ
        assert len(s.elements) == 32 and f.n_dihedrals == 7
