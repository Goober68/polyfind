"""Tests for :mod:`polyfind.mechanics`.

Four groups: that adding this module changed no energy at all, that the strain machinery
is what it claims to be (the reachable/unreachable split, the map and its analytic
derivative), that the physics comes out on the rigid path -- the two piezoelectric routes
agree, a crystal with no dipoles has exactly no piezoelectric response, and the axial
constant there is a report of an invented bend constant rather than of the potential -- and
that the deformable path removes exactly that contamination: with fitted valence terms in
the kernel, ``C_33`` stops moving when the invented constant is swept, the reference becomes
axially stress-free, and the diagonal piezoelectric columns appear for a helix and provably
cannot for a planar zigzag.
"""
import numpy as np
import pytest

from polyfind import mechanics as M
from polyfind.fitting import FITTED_VALENCE
from polyfind.forcefield import SimpleFF
from polyfind.pack import CrystalPacker, periodic_chain
from polyfind.polymers import PE, PVDF, THREE_STATE

T, GP, GM = 0, 1, 2

# One fixed configuration and the energies a default CrystalPacker gives for it, without a
# field and with one.  Compared with ``==``, not a tolerance: mechanics.py must not have
# moved a single bit of the kernel.
#
# Re-recorded once, and the guard worked as intended when it happened.  These were first
# taken while PVDF's backbone C-C bond length was 1.54 A; correcting it to the measured
# 1.528 A moved every PVDF geometry and so every PVDF energy here, and these three
# assertions fired.  PE is untouched by that correction and its value is unchanged from
# the original recording, which is the useful control: it shows the shift came from the
# polymer definition and not from anything in this module.
FIXED_PARAMS = np.array([5.0, 9.6, 92.0, 20.0, 50.0, 1.0, 0.0])
FIXED_FIELD = (0.01, -0.02, 0.005)
RECORDED_ENERGY = {
    "PE": (-4.143020937018386, -4.143020937018386),
    "beta": (-8.784383912060322, -8.718173227934876),
    "alpha": (27.715227596863407, 27.29006207153844),
    "gamma": (27.70619594115115, 26.855864890501213),
}
CHAINS = {
    "PE": (PE, [T]),
    "beta": (PVDF, [T, T]),
    "alpha": (PVDF, [T, GP, T, GM]),
    "gamma": (PVDF, [T, T, T, GP, T, T, T, GM]),
}


@pytest.mark.parametrize("name", sorted(RECORDED_ENERGY))
def test_default_packer_energies_are_bit_for_bit_unchanged(name):
    polymer, seq = CHAINS[name]
    chain = periodic_chain(polymer, seq, THREE_STATE)
    e0, e1 = RECORDED_ENERGY[name]
    assert float(CrystalPacker(chain).energy(FIXED_PARAMS[None])[0]) == e0
    assert float(CrystalPacker(chain, field=FIXED_FIELD).energy(FIXED_PARAMS[None])[0]) == e1
    # and the analytic-gradient path still returns the same value from the same expressions
    assert CrystalPacker(chain).energy_and_grad(FIXED_PARAMS)[0] == e0


# --- what a strain can be -------------------------------------------------------------
@pytest.mark.parametrize("K", M.UNREACHABLE)
def test_chain_axis_shears_are_refused_by_name(K):
    eps = np.zeros(6)
    eps[K] = 1e-3
    with pytest.raises(ValueError, match=r"cannot be applied"):
        M.strained_cell([4.6, 8.6, 90.0, 0.0, 0.0, 0.0, 0.0], 2.63, eps)


def test_reachable_and_unreachable_partition_the_voigt_indices():
    assert sorted([*M.IN_PLANE, M.AXIAL, *M.UNREACHABLE]) == list(range(6))


@pytest.mark.parametrize("K", M.IN_PLANE)
@pytest.mark.parametrize("h", [3e-3, -3e-3])
def test_strain_map_round_trips(K, h):
    p0 = np.array([4.6, 8.59, 93.0, 10.0, 200.0, 1.0, 1.0])
    eps = np.zeros(6)
    eps[K] = h
    sc = M.strained_cell(p0, 2.63, eps)
    back = M.strain_from_cell(p0[0], p0[1], p0[2], sc.params[0], sc.params[1], sc.params[2])
    assert back == pytest.approx(eps, abs=1e-9)
    # a shear leaves a_vec off x, and the canonicalising rotation is folded into both phi
    assert sc.theta == pytest.approx(0.0 if K != 5 else np.degrees(-h / 2), abs=1e-5)
    assert sc.params[3] - p0[3] == pytest.approx(sc.theta)
    assert sc.params[4] - p0[4] == pytest.approx(sc.theta)


def test_axial_strain_is_the_only_thing_that_moves_c():
    p0 = np.array([4.6, 8.59, 90.0, 0.0, 0.0, 0.0, 0.0])
    for K in M.IN_PLANE:
        eps = np.zeros(6)
        eps[K] = 5e-3
        assert M.strained_cell(p0, 2.63, eps).c == pytest.approx(2.63, abs=1e-12)
    eps = np.zeros(6)
    eps[M.AXIAL] = 5e-3
    assert M.strained_cell(p0, 2.63, eps).c == pytest.approx(2.63 * 1.005)


# --- references, built once ------------------------------------------------------------
@pytest.fixture(scope="module")
def beta():
    ref, _ = M.refined_reference(PVDF, [T, T], label="beta")
    return M.relax_reference(ref)


@pytest.fixture(scope="module")
def pe():
    ref, _ = M.refined_reference(PE, [T], label="PE")
    return M.relax_reference(ref)


def test_analytic_stress_matches_a_finite_difference_of_the_energy(beta):
    _, sigma = M.stress(beta)
    V0 = beta.volume
    for K in (*M.IN_PLANE, M.AXIAL):
        h = 1e-5
        eps = np.zeros(6)
        eps[K] = h

        def energy(e):
            sc = M.strained_cell(beta.params, beta.c, e)
            same_c = abs(sc.c - beta.packer.chain.c) < 1e-12
            return float(beta.packer.energy(sc.params[None], coords=None if same_c else beta.packer.chain.coords,
                                            c=None if same_c else np.array([sc.c]))[0])

        fd = (energy(eps) - energy(-eps)) / (2 * h) / V0 * M.KCAL_MOL_A3_TO_GPA
        assert sigma[K] == pytest.approx(fd, abs=1e-4, rel=1e-4)
    assert np.isnan(sigma[list(M.UNREACHABLE)]).all()


def test_relaxed_reference_has_no_in_plane_stress(beta):
    _, sigma = M.stress(beta)
    assert np.abs(sigma[list(M.IN_PLANE)]).max() < 1e-3  # GPa


def test_elastic_block_is_symmetric_and_stable(beta):
    el = M.elastic_constants(beta)
    assert el.asymmetry < 0.2  # GPa, before symmetrising: a check on the finite difference
    assert np.allclose(el.block, el.block.T)
    assert np.linalg.eigvalsh(el.block).min() > 0.0  # the in-plane block is positive definite
    assert np.isnan(el.C[2, 2])  # C_33 is not reported
    for K in M.UNREACHABLE:
        assert np.isnan(el.C[K]).all() and np.isnan(el.C[:, K]).all()


def test_step_dependence_of_C22_is_the_lennard_jones_cutoff(beta):
    """The residual step dependence at the default cutoff is diagnosed, not tolerated.

    The LJ term is energy-shifted but not force-shifted, so its force jumps at ``r = rc``;
    beta's ``b = 8.59`` sits on the default 8 A cutoff, and straining ``b`` walks pairs
    across the jump.  Lengthening the cutoff should make the finite difference behave.
    """
    def spread(ref):
        return abs(M.elastic_constants(ref, step=1e-3).C[1, 1] - M.elastic_constants(ref, step=4e-3).C[1, 1])

    long_cut = M.relax_reference(M.reference_from_chain(beta.packer.chain, beta.params, n_chains=2, cutoff=12.0))
    assert spread(beta) > 1.0  # GPa, at the default cutoff
    assert spread(long_cut) < 0.3  # and gone at 12 A


def test_proper_and_improper_e_differ_by_exactly_the_polarization(beta):
    el = M.elastic_constants(beta)
    pz = M.piezoelectric(beta, el)
    P = beta.polarization()
    diff = pz.e_improper - pz.e
    # -P_i on the two diagonal columns (xx, yy) and nothing at all on the shear column
    assert diff[:, 0] == pytest.approx(-P, abs=2e-4)
    assert diff[:, 1] == pytest.approx(-P, abs=2e-4)
    assert diff[:, 2] == pytest.approx(np.zeros(3), abs=2e-4)


def test_direct_and_converse_piezoelectric_routes_agree(beta):
    el = M.elastic_constants(beta)
    pz = M.piezoelectric(beta, el)
    assert np.abs(pz.d_from_e).max() > 1.0  # there is something to check
    assert pz.relative_difference < 0.03
    assert pz.d_direct == pytest.approx(pz.d_from_e, abs=0.05 * np.abs(pz.d_from_e).max() + 1e-6)


def test_polyethylene_has_exactly_no_piezoelectric_response(pe):
    """The symmetry check, and not a trivial one.

    PE's atoms do carry charges (C -0.12, H +0.06); what vanishes is the *cell dipole*, and
    it vanishes for every configuration and not only at the minimum, because each CH2 group
    is neutral and the all-trans repeat is centrosymmetric.  So the piezoelectric response
    has to be zero to machine precision through a cancellation, which is a real test of the
    whole chain -- dipole, frame rotation, both derivative routes -- rather than of a
    multiplication by zero.
    """
    assert np.abs(pe.packer._q_cell).max() > 0.0
    assert pe.packer.cell_charge == 0.0
    rng = np.random.default_rng(1)
    for _ in range(4):
        p = np.array([rng.uniform(4, 9), rng.uniform(4, 9), rng.uniform(70, 110), rng.uniform(0, 360),
                      rng.uniform(0, 360), rng.uniform(0, pe.c), float(rng.integers(0, 2))])
        assert np.abs(pe.packer.dipole(p[None])[0]).max() < 1e-14
    assert np.abs(pe.polarization()).max() < 1e-14
    resp = M.electromechanical_response(pe, relax_first=False)
    assert np.abs(resp.piezo.e).max() < 1e-12
    assert np.abs(resp.piezo.e_improper).max() < 1e-12
    assert np.abs(resp.piezo.d_from_e).max() < 1e-9
    assert np.abs(resp.piezo.d_direct).max() < 1e-9
    assert np.abs(resp.actuator.free_strain).max() < 1e-12
    assert np.abs(resp.actuator.blocking_stress).max() < 1e-12
    assert resp.actuator.work_density < 1e-12
    # the elastic constants, which owe nothing to the charges, are perfectly ordinary
    assert np.linalg.eigvalsh(resp.elastic.block).min() > 0.0


def test_actuator_figures_follow_from_C_and_d(beta):
    el = M.elastic_constants(beta)
    pz = M.piezoelectric(beta, el)
    ac = M.actuator(el, pz, [0.0, 1.0, 0.0], 0.01)
    assert np.linalg.norm(ac.direction) == pytest.approx(1.0)
    assert ac.free_strain == pytest.approx(pz.d_from_e.T @ (ac.direction * ac.field_v_per_m) * 1e-12)
    assert ac.blocking_stress == pytest.approx(el.block @ ac.free_strain)
    w = float(ac.blocking_stress @ ac.free_strain)
    assert ac.work_density == pytest.approx(0.25 * w * M.GPA_TO_KJ_PER_M3)
    assert ac.elastic_energy == pytest.approx(2.0 * ac.work_density)
    # the polarization direction is a *bad* poling direction here: a rigid dipole already
    # aligned with the field feels no torque, and reorientation is the only channel there is
    along_P = M.actuator(el, pz, beta.polarization(), 0.01)
    assert along_P.work_density < 1e-6 * ac.work_density


def test_axial_c33_is_a_report_of_the_invented_bend_constant(beta):
    ax = M.axial_report(PVDF, beta, n_points=5)
    assert ax.geometrically_reachable and not ax.computable
    assert ax.strain_range[0] < 0.0 < ax.strain_range[1]  # c does move, through the angles
    assert ax.sigma_zz_rigid < -1.0  # GPa: the reference is nowhere near axially stress free
    ks = sorted(ax.c33)
    values = np.array([ax.c33[k] for k in ks])
    # C_33 is affine in the bend constant, to the precision of the parabola fit
    slope = (values[-1] - values[0]) / (ks[-1] - ks[0])
    assert values == pytest.approx(values[0] + slope * (np.array(ks) - ks[0]), rel=1e-3)
    assert ax.bend_fraction > 0.5  # most of it is the constant, not the potential


# --- the deformable chain: what valence terms in the kernel buy -----------------------------
def _valence_ff():
    import polyfind.fitting  # noqa: F401 - registers the fitted presets

    return SimpleFF.from_preset("pvdf-dft-valence")


def _deformable(polymer, seq, label, angle_stiffness=None):
    """Refine with valence terms in the kernel, then relax cell *and* chain to zero stress."""
    ff = _valence_ff()
    refine_kw = {"valence": ff}
    if angle_stiffness is not None:
        refine_kw["angle_stiffness"] = angle_stiffness
    ref, rr = M.refined_reference(polymer, seq, label=label, valence=ff, refine_kw=refine_kw)
    shape = M.shape_of(polymer, ref, angles=rr.angles)
    return M.relax_reference_deformable(ref, shape)


@pytest.fixture(scope="module")
def beta_deformable():
    with FITTED_VALENCE.applied():
        yield _deformable(PVDF, [T, T], "beta")


@pytest.fixture(scope="module")
def alpha_deformable():
    with FITTED_VALENCE.applied():
        yield _deformable(PVDF, [T, GP, T, GM], "alpha")


def test_a_shape_relaxation_without_valence_terms_is_refused(beta):
    """The rigid kernel has no restoring force for an angle, and says so rather than answering."""
    shape = M.shape_of(PVDF, beta, angles=[112.0, 112.0])
    with pytest.raises(ValueError, match="needs a packer with valence terms"):
        M.relax_deformable(beta, shape, beta.params, shape.x0)


def test_the_deformable_reference_is_stress_free_along_the_chain_too(beta_deformable):
    """What the rigid path could not have: a reference at zero *axial* stress.

    The multiplier of the constraint that holds ``c`` is ``dE/dc`` of the relaxed crystal,
    so a reference relaxed over ``c`` has it at zero -- against the rigid path's residual
    ``sigma_zz`` of several GPa, which was the sign that nothing resisted the deformation.
    """
    ref, shape = beta_deformable
    with FITTED_VALENCE.applied():
        el = M.elastic_constants(ref, shape=shape)
    assert np.abs(el.residual_stress[[0, 1, 2, 5]]).max() < 1e-4  # GPa, all four reachable
    assert el.reachable == M.WITH_AXIAL
    assert np.allclose(el.block, el.block.T)
    assert np.linalg.eigvalsh(el.block).min() > 0.0
    assert el.asymmetry < 0.05  # GPa: the multiplier route against the analytic-stress route
    assert el.C[2, 2] == pytest.approx(el.c33_from_energy, rel=2e-3)  # and against the curvature
    assert 100.0 < el.C[2, 2] < 1000.0


def test_c33_is_no_longer_affine_in_the_invented_bend_constant():
    """The test that the number has become physical.

    On the rigid path ``C_33`` was exactly affine in ``refine_crystal``'s invented
    ``angle_stiffness``, four fifths of beta's value coming from it.  Here the same constant
    is swept over the same range: it still moves the structure the refinement hands over --
    the assertion on the refined ``c`` checks that it does -- but the relaxation against the
    *fitted* valence terms puts it back, and ``C_33`` does not move at all.
    """
    with FITTED_VALENCE.applied():
        out = []
        for k in (0.0, 210.0):
            ref, shape = _deformable(PVDF, [T, T], "beta", angle_stiffness=k)
            out.append((ref.c, M.elastic_constants(ref, shape=shape).C[2, 2]))
    (c_soft, c33_soft), (c_stiff, c33_stiff) = out
    assert c33_soft == pytest.approx(c33_stiff, rel=1e-4)
    assert c_soft == pytest.approx(c_stiff, rel=1e-6)  # the relaxed repeat forgets it as well


def test_a_planar_zigzag_cannot_change_its_dipole_by_bending(beta_deformable):
    """Why beta's diagonal piezoelectric columns are still exactly zero.

    Not "the chain is rigid" any more -- it deforms.  With bond-charge-increment charges the
    cell dipole is a sum over bonds, every bond length is fixed, and the mirror symmetry of a
    planar all-trans zigzag pins the bisector of each pendant pair perpendicular to the chain
    axis whatever the backbone angle is.  So the one internal coordinate beta has moves ``c``
    and leaves ``mu`` exactly where it was, and no dilation can be piezoelectric.
    """
    ref, shape = beta_deformable
    with FITTED_VALENCE.applied():
        mus = []
        for dx in (-2.0, 0.0, 2.0):
            chain = shape.chains(np.full((1, shape.n), dx))[0]
            ref.packer.update_chain(chain)
            mus.append((chain.c, ref.packer.dipole(ref.params[None])[0].copy()))
    cs = [c for c, _ in mus]
    assert max(cs) - min(cs) > 0.05  # the chain really did change length
    assert np.abs(mus[0][1] - mus[2][1]).max() < 1e-10  # and the dipole did not move at all


def test_beta_has_no_diagonal_piezoelectric_column_and_alpha_does(beta_deformable, alpha_deformable):
    """The informative negative and the positive beside it.

    Beta's chain has one shape parameter and the fixed-``eps_zz`` constraint uses it up, and
    even the axial column cannot help (the test above).  Alpha's helix has three, so a
    dilation does relax the conformation and does change the dipole -- the first non-zero
    diagonal column this package has produced.
    """
    with FITTED_VALENCE.applied():
        out = {}
        for name, (ref, shape) in (("beta", beta_deformable), ("alpha", alpha_deformable)):
            el = M.elastic_constants(ref, shape=shape)
            out[name] = (el, M.piezoelectric(ref, el, shape=shape, converse=(name == "alpha")))
    diag = [M.WITH_AXIAL.index(K) for K in (0, 1, 2)]
    el_b, pz_b = out["beta"]
    el_a, pz_a = out["alpha"]
    assert np.abs(pz_b.e[:, diag]).max() < 1e-6  # C/m^2: exactly zero, not merely small
    assert np.abs(pz_a.e[:, diag]).max() > 1e-2
    # and the two routes still agree where there is something to agree about
    assert pz_a.relative_difference < 0.02
    # the improper column of a rigid dipole array is still exactly -P, on both
    for pz, ref_shape in ((pz_b, beta_deformable), (pz_a, alpha_deformable)):
        with FITTED_VALENCE.applied():
            P = ref_shape[0].polarization()
        for n in diag:
            assert pz.e_improper[:, n] - pz.e[:, n] == pytest.approx(-P, abs=2e-3)


def test_polyethylene_has_no_piezoelectric_response_when_the_chain_deforms_either():
    """The null control, on the harder path: the cancellation has to survive a relaxation."""
    with FITTED_VALENCE.applied():
        ref, shape = _deformable(PE, [T], "PE")
        resp = M.electromechanical_response(ref, shape=shape, relax_first=False)
    assert resp.elastic.reachable == M.WITH_AXIAL
    assert np.abs(resp.piezo.e).max() < 1e-12
    assert np.abs(resp.piezo.d_from_e).max() < 1e-9
    assert np.abs(resp.piezo.d_improper).max() < 1e-9
    # The converse route survives too, at 1e-12 pC/N -- but only because its Newton stops at
    # the reference's own residual stress rather than at a fixed 1e-8 GPa.  Chasing a
    # tolerance below the constrained relaxation's noise floor ran to the iteration cap and
    # came back at 1e-6 instead, which is how the floor came to be measured.
    assert np.abs(resp.piezo.d_direct).max() < 1e-9
    assert resp.actuator.work_density < 1e-12
    assert np.linalg.eigvalsh(resp.elastic.block).min() > 0.0  # the elastic block is ordinary
    assert 100.0 < resp.elastic.C[2, 2] < 1000.0


def test_d_improper_is_the_dimensional_term_and_nothing_else(beta_deformable):
    """``d_improper`` on a diagonal column is ``-P`` contracted with the compliance.

    That is the thickness effect of the Broadhurst-Davis model and not a piezoelectric
    constant; the test pins the identity so that nobody reads it as one.
    """
    ref, shape = beta_deformable
    with FITTED_VALENCE.applied():
        el = M.elastic_constants(ref, shape=shape)
        pz = M.piezoelectric(ref, el, shape=shape, converse=False)
        P = ref.polarization()
    diag = [M.WITH_AXIAL.index(K) for K in (0, 1, 2)]
    for n in diag:
        expect = -float(P[0]) * el.S[diag, n].sum() * M.C_PER_M2_PER_GPA_TO_PC_PER_N
        assert pz.d_improper[0, n] == pytest.approx(expect, rel=1e-3, abs=1e-6)
    # A fixed dipole array always loses P on dilation -- but only once the poling axis is
    # taken along +P.  Which way P points is the domain the packing happened to land in, so
    # the sign is a statement about the *projection* and not about a column of the array.
    n_hat = P / np.linalg.norm(P)
    assert (n_hat @ pz.d_improper[:, diag]).max() < 0.0


# --- the Lennard-Jones cutoff ---------------------------------------------------------------
def test_force_shifting_the_cutoff_removes_the_step_dependence(beta):
    """The diagnosis of :func:`test_step_dependence_of_C22_is_the_lennard_jones_cutoff`, acted on.

    Force-shifting makes the pair force continuous at ``rc``, and the finite difference stops
    seeing a jump.  It does *not* make the value converge faster -- it is further from the
    long-cutoff answer at the same ``rc``, because it subtracts a tail as well -- so it is a
    fix for the derivative, not for the truncation.
    """
    def spread(ref):
        return abs(M.elastic_constants(ref, step=1e-3).C[1, 1] - M.elastic_constants(ref, step=4e-3).C[1, 1])

    shifted = M.relax_reference(M.reference_from_chain(beta.packer.chain, beta.params, n_chains=2,
                                                       lj_cutoff="force"))
    assert spread(beta) > 1.0  # GPa, energy-shifted at the default cutoff
    assert spread(shifted) < 0.5  # and much smaller once the force is continuous too
    assert M.elastic_constants(shifted).C[1, 1] < M.elastic_constants(beta).C[1, 1] + 1.0
