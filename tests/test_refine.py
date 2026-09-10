import numpy as np
import pytest
from scipy.optimize import approx_fprime

from polyfind.pack import (
    CrystalPacker,
    EV_TO_KCAL,
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


# --------------------------------------------- analytic gradients (docs section 13)
@pytest.mark.parametrize("poly,seq,cell,par", [
    (PE, [T], [4.93, 7.40, 90.0, 45.0, 135.0, 0.60, 0.0], "linegroup"),
    (PE, [T], [4.93, 7.40, 90.0, 45.0, 135.0, 0.60, 0.0], "free"),
    (PVDF, [T, T], [4.65, 8.61, 90.0, 8.0, 26.0, 1.10, 0.0], "linegroup"),
    (PVDF, [T, T], [4.65, 8.61, 90.0, 8.0, 26.0, 1.10, 0.0], "free"),
    (PVDF, [T, GP, T, GM], [5.29, 8.99, 90.0, 209.5, 150.0, 0.40, 1.0], "linegroup"),
    (PVDF, [T, GP, T, GM], [5.29, 8.99, 90.0, 209.5, 150.0, 0.40, 1.0], "free"),
    (PVDF, [T, T, T, GP, T, T, T, GM], [4.96, 9.67, 90.0, 30.0, 200.0, 2.50, 0.0], "linegroup"),
    (PVDF, [T, T, T, GP, T, T, T, GM], [4.96, 9.67, 90.0, 30.0, 200.0, 2.50, 0.0], "free"),
    (PVDF, [T, T, GP, T, GM, T], [5.5, 9.0, 90.0, 100.0, 100.0, 1.0, 0.0], "free"),
])
@pytest.mark.parametrize("field", [None, (0.12, -0.07, 0.2)])
def test_analytic_refinement_gradient_matches_central_differences(poly, seq, cell, par, field):
    """Every component of the analytic gradient, against a central difference of the
    objective L-BFGS-B actually minimises (energy + restraints), at two points.

    The cell variables are analytic outright and agree to ~1e-7; the conformational ones
    still contract the exact ``dE/d(coords)`` with a finite-difference Jacobian of the
    chain build, so they agree to ~1e-5.  Steps are swept and the closest taken: a coarse
    step straddles the Lennard-Jones cutoff, where the *difference* has a jump the
    derivative correctly does not.
    """
    chain = periodic_chain(poly, seq, THREE_STATE)
    start = CrystalPacker(chain, field=field).result(np.array(cell))
    probe = {}
    refine_crystal(poly, start, maxfev=1, maxiter=1, parametrisation=par, field=field,
                   gamma_free=True, probe=probe)
    f, van, vfd = probe["value"], probe["value_and_grad_analytic"], probe["value_and_grad_fd"]
    n_cell, n_shape = probe["n_cell"], probe["n_shape"]
    rng = np.random.default_rng(5)
    for x in (probe["x0"],
              probe["x0"] + np.concatenate([rng.uniform(-0.15, 0.15, n_cell), rng.uniform(-3, 3, n_shape)])):
        E, g = van(x)
        E_fd, g_fd = vfd(x)
        assert E == f(x) == E_fd  # the value is the same function, bit for bit
        scale = max(np.abs(g).max(), 1e-6)
        for i in range(len(x)):
            steps = (3e-2, 1e-2, 3e-3, 1e-3) if (i >= n_cell or probe["labels"][i] not in ("a", "b", "dz")) \
                else (3e-3, 1e-3, 3e-4, 1e-4)
            best = None
            for h in steps:
                xp_, xm = x.copy(), x.copy()
                xp_[i] += h
                xm[i] -= h
                d = (f(xp_) - f(xm)) / (2 * h)
                if best is None or abs(d - g[i]) < abs(best - g[i]):
                    best = d
            tol = (1e-5 if i < n_cell else 1e-3) * abs(g[i]) + 1e-5 * scale
            assert abs(best - g[i]) <= tol, (probe["labels"][i], g[i], best)
        # ... and the two routes agree with each other on every component
        np.testing.assert_allclose(g, g_fd, atol=2e-3 * scale)


def test_refinement_reaches_the_same_minimum_on_either_gradient_route():
    chain = periodic_chain(PVDF, [T, GP, T, GM], THREE_STATE)
    packer = CrystalPacker(chain)
    start = packer.result(np.array([5.29, 8.99, 90.0, 209.5, 209.5, 0.0, 0.0]))
    an = refine_crystal(PVDF, start, maxfev=60, maxiter=60)
    fd = refine_crystal(PVDF, start, maxfev=60, maxiter=60, gradient="fd")
    assert an.result.energy_per_monomer == pytest.approx(fd.result.energy_per_monomer, abs=1e-3)
    assert an.result.a == pytest.approx(fd.result.a, abs=5e-3)
    assert an.result.b == pytest.approx(fd.result.b, abs=5e-3)
    assert an.result.c == pytest.approx(fd.result.c, abs=5e-3)
    np.testing.assert_allclose(an.torsions, fd.torsions, atol=0.5)
    with pytest.raises(ValueError):
        refine_crystal(PVDF, start, gradient="autodiff")


def test_analytic_gradient_costs_one_kernel_row_per_iteration():
    """The point of the change: rows per iteration drop from 1 + 2 n_vars to 1."""
    chain = periodic_chain(PVDF, [T, T, T, GP, T, T, T, GM], THREE_STATE)
    start = CrystalPacker(chain).result(np.array([4.96, 9.67, 90.0, 30.0, 200.0, 2.50, 0.0]))
    rows, evals, nvars = {}, {}, {}
    for g in ("fd", "analytic"):
        probe = {}
        res = refine_crystal(PVDF, start, maxfev=20, maxiter=20, gradient=g, probe=probe)
        rows[g] = probe["packer"].n_energy_rows
        evals[g] = res.n_evaluations
        nvars[g] = res.n_variables
    assert nvars["fd"] == nvars["analytic"]
    assert rows["fd"] >= (1 + 2 * nvars["fd"]) * evals["fd"]
    assert rows["analytic"] <= evals["analytic"] + 2  # the final result() evaluation
    assert rows["analytic"] * 10 < rows["fd"]


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


# ------------------------------------------------------------------ applied field
def test_refinement_honours_and_inherits_an_applied_field():
    """The field reaches the refinement, and is not silently dropped by a plain call."""
    chain = periodic_chain(PVDF, [T, T], THREE_STATE)
    p = np.array([4.65, 8.61, 90.0, 0.0, 0.0, 1.29, 0.0])
    plain = CrystalPacker(chain)
    mu = plain.dipole(p[None])[0]
    E = 0.3 * mu / np.linalg.norm(mu)
    packer = CrystalPacker(chain, field=E)
    start = packer.result(p)
    assert start.field == pytest.approx(tuple(E)) and start.field_energy < 0.0

    res = refine_crystal(PVDF, start, maxfev=20, maxiter=20)  # field inherited from start
    assert res.result.field == pytest.approx(tuple(E))
    assert res.result.energy_per_monomer <= start.energy_per_monomer + 1e-9
    # the refined result carries the coupling of its own (relaxed) dipole
    assert res.result.field_energy == pytest.approx(-float(res.result.dipole @ E) * EV_TO_KCAL, rel=1e-12)
    assert res.result.field_energy < 0.0

    # an explicit zero field overrides the inheritance
    off = refine_crystal(PVDF, start, field=(0.0, 0.0, 0.0), maxfev=20, maxiter=20)
    assert off.result.field == (0.0, 0.0, 0.0) and off.result.field_energy == 0.0
    # the field really changed the answer: it beats the fieldless minimum by roughly -mu.E
    assert res.result.energy_per_cell < off.result.energy_per_cell + res.result.field_energy + 1e-6

    # a fieldless start refines exactly as before: no field anywhere
    plain_start = plain.result(p)
    base = refine_crystal(PVDF, plain_start, maxfev=20, maxiter=20)
    assert base.result.field == (0.0, 0.0, 0.0) and base.result.field_energy == 0.0
    # ... and the opposed field costs the same as the aligned one gains, at the start cell
    against = CrystalPacker(chain, field=-E).result(p)
    assert against.field_energy == pytest.approx(-start.field_energy, rel=1e-12)


def test_nelder_mead_works_with_the_line_group_too():
    chain = periodic_chain(PVDF, [T, GP, T, GM], THREE_STATE)
    packer = CrystalPacker(chain)
    start = packer.result(np.array([5.29, 8.99, 90.0, 209.5, 209.5, 0.0, 0.0]))
    res = refine_crystal(PVDF, start, maxfev=80, method="nelder-mead")
    assert res.parametrisation == "linegroup"
    assert res.rotation_error < 1e-5
    assert res.result.energy_per_monomer <= start.energy_per_monomer + 1e-9
