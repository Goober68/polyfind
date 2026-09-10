"""Tests for :mod:`polyfind.mechanics`.

Three groups: that adding this module changed no energy at all, that the strain machinery
is what it claims to be (the reachable/unreachable split, the map and its analytic
derivative), and that the physics comes out -- the two piezoelectric routes agree, a
crystal with no dipoles has exactly no piezoelectric response, and the axial constant is a
report of an invented bend constant rather than of the potential.
"""
import numpy as np
import pytest

from polyfind import mechanics as M
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
