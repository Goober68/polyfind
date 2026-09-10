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


# ================================ valence terms, Cartesian forces, off-site charges
def _valence_x(seed=5, offset=0.16):
    """A parameter vector off the starting point, inside the bounds, with a real offset."""
    rng = np.random.default_rng(seed)
    x = np.clip(F.VAL_X0 + 0.3 * rng.normal(0, F.VAL_SIGMA), F.VAL_BOUNDS[:, 0], F.VAL_BOUNDS[:, 1])
    x[-1] = offset
    return x


def test_the_valence_vector_extends_the_previous_one_exactly():
    """The first 27 parameters are the previous fit's, unchanged and in the same order, so
    "what did the new terms buy" is a question about this vector and not about two."""
    assert F.VAL_VARIABLES[:len(F.REF_VARIABLES)] == F.REF_VARIABLES
    assert np.array_equal(F.VAL_X0[:len(F.REF_X0)], F.REF_X0)
    assert np.array_equal(F.VAL_BOUNDS[:len(F.REF_X0)], F.REF_BOUNDS)
    # and the valence block starts at textbook values, not at averages of the reference data
    p = F.valence_unpack(F.VAL_X0)
    assert p["offset"] == 0.0
    assert dict(zip(F.VAL_BOND_TYPES, p["bond_r0"]))["C-H"] == pytest.approx(1.090)
    assert F.valence_ff_parameters(F.VAL_X0).charge_offsets == ()  # zero offset ships as none


@pytest.mark.parametrize("offset", [0.0, 0.17, -0.21])
def test_valence_design_reproduces_the_reference_implementation(offset):
    """Energies, *Cartesian forces* and torques, against ``SimpleFF``'s own finite differences.

    Same contract as the rigid design: the fast design exists only for speed, and if it
    drifted from the shipped potential the fit would be optimising something else.  The
    forces are the new part -- they are what the objective now uses instead of the torsional
    projection -- so they are checked component by component.
    """
    frames = _synthetic_frames()
    design = F.ValenceDesign(frames)
    x = _valence_x(offset=offset)
    p = F.valence_unpack(x)
    energy, forces, torque = design.evaluate(p)
    ff = F.valence_ff_parameters(x).simple_ff()
    for i, f in enumerate(frames):
        V = p["torsion"][F.torsion_types(f)]
        assert energy[i] == pytest.approx(ff.energy_frame(f, torsion=V), abs=1e-8)
        assert np.allclose(forces[design.atom_frame == i], ff.forces_frame(f, torsion=V), atol=2e-5)
        assert np.allclose(torque[design.dih_frame == i], ff.frame_torques(f, torsion=V), atol=1e-4)


def test_the_valence_design_carries_the_full_force_signal():
    """The reason the objective changed: with valence terms the reference forces are no
    longer 93% unrepresentable, so the residual block is Cartesian and not a projection."""
    design = F.ValenceDesign(_synthetic_frames())
    x = _valence_x()
    r_no_force = design.residuals(x, force_weight=0.0, ridge=0.0)
    r_force = design.residuals(x, force_weight=0.03, ridge=0.0)
    assert r_no_force.size == design.n_frames
    assert r_force.size == design.n_frames + 3 * design.n_atoms
    assert design.bond_r.size > 0 and design.angle_theta.size > 0


def test_the_analytic_angle_and_dihedral_gradients_match_finite_differences():
    """A sign error in either would quietly corrupt every fitted valence and torsion term."""
    from polyfind.chain import angle as chain_angle
    from polyfind.forcefield import dihedral_angles

    rng = np.random.default_rng(2)
    c = rng.normal(size=(4, 3))
    theta, g = F._angle_gradient(c, np.array([[0, 1, 2]]))
    phi, gd = F._dihedral_gradient(c, np.array([[0, 1, 2, 3]]))
    assert np.degrees(theta[0]) == pytest.approx(chain_angle(c, 0, 1, 2))
    assert np.degrees(phi[0]) == pytest.approx(dihedral_angles(c, [[0, 1, 2, 3]])[0])
    h = 1e-6
    for i in range(4):
        for a in range(3):
            cp, cm = c.copy(), c.copy()
            cp[i, a] += h
            cm[i, a] -= h
            d = np.radians(dihedral_angles(cp, [[0, 1, 2, 3]])[0] - dihedral_angles(cm, [[0, 1, 2, 3]])[0])
            assert gd[0, i, a] == pytest.approx(d / (2 * h), abs=1e-6)
            if i < 3:
                da = np.radians(chain_angle(cp, 0, 1, 2) - chain_angle(cm, 0, 1, 2))
                assert g[0, i, a] == pytest.approx(da / (2 * h), abs=1e-6)


def test_the_charge_permittivity_degeneracy_survives_the_off_site_charge():
    """Measured, not assumed.  Moving a charge site does not touch the argument that only
    ``q_i q_j / eps_r`` enters an energy, so the degeneracy is expected to be exactly as
    exact as it was -- and the offset is expected to be separable from the charge, because
    only their *product* is fixed at long range and these contacts are not long range."""
    design = F.ValenceDesign(_synthetic_frames())
    d = F.valence_degeneracies(design, _valence_x())
    assert d["charge_permittivity_max_energy_change"] < 1e-9  # still exact
    assert d["charge_offset_max_energy_change"] > 0.1  # and this one is not a degeneracy
    assert d["max_bent_angle_deg"] < 170.0  # no harmonic is being asked to work near linear


def test_valence_terms_do_not_reach_the_packer_and_the_offset_does():
    """Inside a packing the valence energy is a constant, so ``applied`` deliberately does
    not forward it; the off-site charges do reach the lattice, as their dipole-equivalent
    projection, because a charge model that stopped at the oligomer would make the crystal
    numbers mean something different from the fitted ones."""
    chain = _alpha_chain()
    p = np.array([[5.09, 9.14, 90.0, 90.0, 270.0, 2.24, 1]])
    from dataclasses import replace

    base = F.FFParameters(charge_increments=(("C", "H", -0.10), ("C", "F", 0.20)))
    with_valence = replace(base, bond_terms=(("C-C", 620.0, 1.526), ("*", 500.0, 1.5)),
                           angle_terms=(("*", 120.0, 109.5),))
    with_offset = replace(base, charge_offsets=(("F", 0.18),))
    with base.applied():
        e0 = pack_mod.CrystalPacker(chain, n_chains=2).energy(p)[0]
        q0 = pack_mod.CrystalPacker(chain, n_chains=2).chain.charges.copy()
    with with_valence.applied():
        assert pack_mod.CrystalPacker(chain, n_chains=2).energy(p)[0] == e0
    with with_offset.applied():
        packer = pack_mod.CrystalPacker(chain, n_chains=2)
        assert packer.energy(p)[0] != e0
        assert packer.chain.charges.sum() == pytest.approx(q0.sum(), abs=1e-12)
        assert not np.allclose(packer.chain.charges, q0)


def test_the_valence_preset_is_registered_and_the_defaults_are_untouched():
    assert "pvdf-dft-valence" in PRESETS
    assert F.FITTED_VALENCE.preset_kwargs() == PRESETS["pvdf-dft-valence"]
    default = SimpleFF()
    assert (default.torsion, default.eps_r, default.charge_scale, default.lj,
            default.charge_increments, default.bond_terms, default.angle_terms,
            default.charge_offsets) == (DEFAULT_TORSION, 1.0, 1.0, None, None, None, None, None)
    fitted = SimpleFF.from_preset("pvdf-dft-valence")
    assert fitted.bond_terms and fitted.angle_terms and fitted.charge_offsets
    assert {t for t, _, _ in fitted.bond_terms} == set(F.VAL_BOND_TYPES)
    assert {t for t, _, _ in fitted.angle_terms} == set(F.VAL_ANGLE_TYPES)
    assert UFF_LJ["F"] == (3.364, 0.050)  # selecting a preset never edits the global table


def test_the_fitted_fluorine_charge_is_physical_and_off_its_bound():
    """The finding the valence terms were supposed to produce, pinned as a number.

    DESIGN.md 5.9 records the previous fit sitting on the 0.02 e floor for ``q_C-F`` and
    turning fluorine *positive* once released.  This asserts the shipped vector is on the
    right side of that and is not being held there by a bound.
    """
    i = F.VAL_VARIABLES.index("q_C-F")
    q = F.VAL_FITTED_X[i]
    lo, hi = F.VAL_BOUNDS[i]
    assert q > 0.10  # fluorine negative, carbon positive, and not marginally
    assert min(abs(q - lo), abs(q - hi)) > 0.05  # nowhere near either bound
    inc = {(a, b): v for a, b, v in F.FITTED_VALENCE.charge_increments}
    assert inc[("C", "H")] < 0 and inc[("C", "F")] > 0  # H positive, F negative
    # and the rigid vector really is the old model: no valence, charge on the nucleus
    assert F.valence_ff_parameters(F.rigid_vector()).charge_offsets == ()
    assert all(k == 0.0 for _, k, _ in F.valence_ff_parameters(F.rigid_vector()).bond_terms)


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
