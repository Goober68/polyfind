import numpy as np
import pytest
from scipy.optimize import approx_fprime

from polyfind.pack import (
    CrystalPacker,
    FD_STEPS,
    cell_value_and_grad,
    periodic_chain,
    periodic_chain_from_torsions,
    repeat_chains_from_torsions,
)
from polyfind.polymers import PE, PVDF, THREE_STATE
from polyfind.refine import refine_crystal

T, GP, GM = 0, 1, 2


def _cells(chain, M=8, seed=0):
    rng = np.random.default_rng(seed)
    return np.column_stack([
        rng.uniform(4.6, 9.0, M), rng.uniform(4.6, 9.0, M), np.full(M, 90.0),
        rng.uniform(0, 360, M), rng.uniform(0, 360, M), rng.uniform(0, chain.c, M),
        rng.integers(0, 2, M).astype(float),
    ])


def test_update_chain_matches_a_fresh_packer():
    ideal = periodic_chain(PVDF, [T, GP, T, GM], THREE_STATE)
    packer = CrystalPacker(ideal)
    params = _cells(ideal)
    for delta in ([0.0, 0.0, 0.0, 0.0], [-3.0, 7.0, 2.0, -5.0], [10.0, -12.0, -4.0, 6.0]):
        chain = periodic_chain_from_torsions(PVDF, "TG+TG-", ideal.dihedrals + np.array(delta))
        packer.update_chain(chain)
        fresh = CrystalPacker(chain)
        assert packer.K == fresh.K and packer.e_torsion == pytest.approx(fresh.e_torsion, abs=1e-12)
        assert packer.e_intra == pytest.approx(fresh.e_intra, abs=1e-12)
        np.testing.assert_allclose(packer.energy(params), fresh.energy(params), rtol=1e-9, atol=1e-9)
    # the topology-dependent tables are shared, not rebuilt
    assert packer.same_site.keys() == CrystalPacker(ideal).same_site.keys()
    with pytest.raises(ValueError):
        packer.update_chain(periodic_chain(PE, [T], THREE_STATE))


def test_per_row_chain_batch_matches_individual_packers():
    ideal = periodic_chain(PVDF, [T, GP, T, GM], THREE_STATE)
    rng = np.random.default_rng(1)
    chains = repeat_chains_from_torsions(
        PVDF, "TG+TG-", ideal.dihedrals + rng.uniform(-12, 12, (6, 4)), align_to=ideal
    )
    params = _cells(ideal, M=6, seed=2)
    packer = CrystalPacker(ideal)
    batched = packer.energy(params, **packer.chain_batch(chains))
    one_by_one = [float(CrystalPacker(c).energy(params[i : i + 1])[0]) for i, c in enumerate(chains)]
    np.testing.assert_allclose(batched, one_by_one, rtol=1e-9, atol=1e-9)


def test_batched_builder_matches_the_single_chain_builder():
    ideal = periodic_chain(PVDF, [T, GP, T, GM], THREE_STATE)
    tors = ideal.dihedrals + np.array([[0.0, 0, 0, 0], [4.0, -3.0, 1.0, 2.0]])
    batch = repeat_chains_from_torsions(PVDF, "TG+TG-", tors)
    for i, ch in enumerate(batch):
        one = periodic_chain_from_torsions(PVDF, "TG+TG-", tors[i])
        np.testing.assert_allclose(ch.coords, one.coords, atol=1e-12)
        assert ch.c == pytest.approx(one.c, abs=1e-12)
        assert ch.rotation_error == pytest.approx(one.rotation_error, abs=1e-12)


@pytest.mark.parametrize("poly,seq,p", [
    (PE, [T], [4.9, 7.4, 90.0, 45.0, 135.0, 0.6, 0.0]),
    (PVDF, [T, GP, T, GM], [5.2, 9.2, 90.0, 210.0, 210.0, 0.5, 0.0]),
])
def test_batched_cell_gradient_matches_finite_differences(poly, seq, p):
    packer = CrystalPacker(periodic_chain(poly, seq, THREE_STATE))
    x = np.array(p)
    free = [0, 1, 3, 4, 5]
    f, grad = cell_value_and_grad(packer, x, free)
    assert f == pytest.approx(float(packer.energy(x[None])[0]), abs=1e-12)

    def fun(xf):
        y = x.copy()
        y[free] = xf
        return float(packer.energy(y[None])[0])

    ref = approx_fprime(x[free], fun, [1e-5, 1e-5, 1e-4, 1e-4, 1e-5])
    scale = max(np.abs(ref).max(), 1e-6)
    np.testing.assert_allclose(grad, ref, atol=0.02 * scale)


def test_batched_torsion_gradient_matches_finite_differences():
    ideal = periodic_chain(PVDF, [T, GP, T, GM], THREE_STATE)
    ref_chain = periodic_chain_from_torsions(PVDF, "TG+TG-", ideal.dihedrals)
    packer = CrystalPacker(ref_chain)
    p = np.array([5.2, 9.2, 90.0, 210.0, 210.0, 0.5, 0.0])
    t0 = ideal.dihedrals + np.array([-2.0, 5.0, 1.0, -4.0])

    def energy_of(t):
        chain = repeat_chains_from_torsions(PVDF, "TG+TG-", np.atleast_2d(t), align_to=ref_chain)[0]
        packer.update_chain(chain)
        return float(packer.energy(p[None])[0])

    h = 0.05
    eye = np.eye(len(t0))
    batch = [t0] + [t0 + h * e for e in eye] + [t0 - h * e for e in eye]
    chains = repeat_chains_from_torsions(PVDF, "TG+TG-", np.array(batch), align_to=ref_chain)
    E = packer.energy(np.repeat(p[None], len(batch), axis=0), **packer.chain_batch(chains))
    grad = np.array([(E[1 + j] - E[1 + len(t0) + j]) / (2 * h) for j in range(len(t0))])
    ref = approx_fprime(t0, energy_of, 1e-3)
    np.testing.assert_allclose(grad, ref, atol=0.05 * max(np.abs(ref).max(), 1e-6))
    # the batched centre row reproduces a plain single-chain evaluation
    assert E[0] == pytest.approx(energy_of(t0), abs=1e-9)


def test_alpha_refinement_lowers_the_energy_and_stays_commensurate():
    chain = periodic_chain(PVDF, [T, GP, T, GM], THREE_STATE)
    packer = CrystalPacker(chain)
    start = packer.result(np.array([5.29, 8.99, 90.0, 209.5, 209.5, 0.0, 0.0]))
    res = refine_crystal(PVDF, start, maxfev=25, maxiter=25)
    assert res.result.energy_per_monomer <= start.energy_per_monomer + 1e-9
    assert res.rotation_error < 0.5
    assert res.n_evaluations <= 35  # maxiter iterations plus a few line-search evaluations
    assert len(res.torsions) == 4
    assert np.abs(res.torsions - chain.dihedrals).max() <= 40.0 + 1e-6
    assert "commensurability" in res.summary()


def test_beta_chain_stays_all_trans():
    chain = periodic_chain(PVDF, [T, T], THREE_STATE)
    packer = CrystalPacker(chain)
    start = packer.result(np.array([4.65, 8.61, 90.0, 0.0, 0.0, 2.58, 0.0]))
    res = refine_crystal(PVDF, start, maxfev=30, maxiter=30)
    np.testing.assert_allclose(res.torsions, [180.0, 180.0], atol=0.5)
    assert res.rotation_error < 1e-3
    assert res.result.energy_per_monomer <= start.energy_per_monomer + 1e-9


def test_nelder_mead_refinement_still_available():
    chain = periodic_chain(PVDF, [T, GP, T, GM], THREE_STATE)
    packer = CrystalPacker(chain)
    start = packer.result(np.array([5.29, 8.99, 90.0, 209.5, 209.5, 0.0, 0.0]))
    res = refine_crystal(PVDF, start, maxfev=60, method="nelder-mead")
    assert 0 < res.n_evaluations <= 60
    assert res.result.energy_per_monomer <= start.energy_per_monomer + 1e-9
    with pytest.raises(ValueError):
        refine_crystal(PVDF, start, method="powell")


def test_refinement_without_torsions_only_moves_the_cell():
    chain = periodic_chain(PVDF, [T, T], THREE_STATE)
    packer = CrystalPacker(chain)
    start = packer.result(np.array([4.8, 8.3, 90.0, 10.0, 10.0, 1.0, 0.0]))
    res = refine_crystal(PVDF, start, refine_torsions=False, maxfev=20, maxiter=20)
    np.testing.assert_allclose(res.torsions, start.dihedrals, atol=1e-12)
    assert res.result.energy_per_monomer <= start.energy_per_monomer + 1e-9


# ------------------------------------------------- line-group parametrisation (item 5)
def test_linegroup_is_the_default_and_needs_no_penalty():
    chain = periodic_chain(PVDF, [T, GP, T, GM], THREE_STATE)
    packer = CrystalPacker(chain)
    start = packer.result(np.array([5.29, 8.99, 90.0, 209.5, 209.5, 0.0, 0.0]))
    lgr = refine_crystal(PVDF, start, maxfev=40, maxiter=40)
    free = refine_crystal(PVDF, start, maxfev=40, maxiter=40, parametrisation="free", refine_angles=False)
    assert lgr.parametrisation == "linegroup" and free.parametrisation == "free"
    # exactly periodic by construction, with fewer variables than the free path
    assert lgr.rotation_error < 1e-5
    assert lgr.n_variables < free.n_variables
    # the penalty term is inactive: it contributes less than a micro-kcal/mol
    assert 5.0 * lgr.rotation_error ** 2 < 1e-6
    # ... and a run with the penalty switched off entirely lands in the same minimum
    off = refine_crystal(PVDF, start, maxfev=40, maxiter=40, penalty=0.0)
    assert off.result.energy_per_monomer == pytest.approx(lgr.result.energy_per_monomer, abs=1e-3)
    assert off.rotation_error < 1e-5
    # the objective (lattice energy + bond-angle strain) is no worse than the free path's
    assert lgr.result.energy_per_monomer + lgr.angle_energy <= free.result.energy_per_monomer + 1e-6


def test_linegroup_refines_the_bond_angles_within_their_bounds():
    chain = periodic_chain(PVDF, [T, GP, T, GM], THREE_STATE)
    packer = CrystalPacker(chain)
    start = packer.result(np.array([5.29, 8.99, 90.0, 209.5, 209.5, 0.0, 0.0]))
    res = refine_crystal(PVDF, start, maxfev=60, maxiter=60, max_angle_change=5.0)
    assert len(res.angles) == PVDF.bonds_per_repeat
    assert np.abs(res.angles - 114.0).max() <= 5.0 + 1e-6
    assert np.abs(res.angles - 114.0).max() > 0.1  # they actually moved
    assert res.angle_energy > 0.0
    assert res.result.c != pytest.approx(chain.c, abs=1e-3)  # and the repeat followed them
    frozen = refine_crystal(PVDF, start, maxfev=60, maxiter=60, refine_angles=False, parametrisation="free")
    np.testing.assert_allclose(frozen.angles, 114.0, atol=1e-12)
    assert frozen.angle_energy == 0.0


def test_beta_linegroup_forces_equal_angles_and_exact_trans():
    """One monomer per repeat: a straight all-trans chain must have equal angles."""
    chain = periodic_chain(PVDF, [T, T], THREE_STATE)
    packer = CrystalPacker(chain)
    start = packer.result(np.array([4.65, 8.61, 90.0, 0.0, 0.0, 2.58, 0.0]))
    res = refine_crystal(PVDF, start, maxfev=40, maxiter=40)
    assert res.n_variables == 6  # 5 cell + a single free parameter, the common angle
    np.testing.assert_allclose(res.torsions, 180.0, atol=1e-5)
    assert res.angles[0] == pytest.approx(res.angles[1], abs=1e-5)
    assert res.rotation_error < 1e-5


def test_linegroup_falls_back_to_the_penalty_method_for_an_unpatterned_sequence():
    chain = periodic_chain(PVDF, [T, T, GP, T, GM, T], THREE_STATE)
    packer = CrystalPacker(chain)
    start = packer.result(np.array([5.5, 9.0, 90.0, 100.0, 100.0, 1.0, 0.0]))
    res = refine_crystal(PVDF, start, maxfev=15, maxiter=15)
    assert res.parametrisation == "free"
    assert res.n_variables == 5 + len(start.dihedrals) + PVDF.bonds_per_repeat
    assert res.result.energy_per_monomer <= start.energy_per_monomer + 1e-9


def test_unknown_parametrisation_rejected():
    chain = periodic_chain(PVDF, [T, T], THREE_STATE)
    packer = CrystalPacker(chain)
    start = packer.result(np.array([4.65, 8.61, 90.0, 0.0, 0.0, 2.58, 0.0]))
    with pytest.raises(ValueError):
        refine_crystal(PVDF, start, parametrisation="symmetry")


def test_nelder_mead_works_with_the_line_group_too():
    chain = periodic_chain(PVDF, [T, GP, T, GM], THREE_STATE)
    packer = CrystalPacker(chain)
    start = packer.result(np.array([5.29, 8.99, 90.0, 209.5, 209.5, 0.0, 0.0]))
    res = refine_crystal(PVDF, start, maxfev=80, method="nelder-mead")
    assert res.parametrisation == "linegroup"
    assert res.rotation_error < 1e-5
    assert res.result.energy_per_monomer <= start.energy_per_monomer + 1e-9
