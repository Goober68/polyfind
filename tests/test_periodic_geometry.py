"""Full affine image geometry: derivatives, reversal, and frame covariance."""
import importlib.util
from pathlib import Path

import numpy as np
import pytest

from polyfind import mechanics as M, pack as pack_mod
from polyfind.backend import to_numpy
from polyfind.ewald import charge_dipole_exclusion, exclusion_correction
from polyfind.forcefield import SimpleFF
from polyfind.pack import CrystalPacker, periodic_chain
from polyfind.periodic_geometry import PlacedCell, repeat_rows, repeat_vector
from polyfind.polarizability import Polarizable
from polyfind.polymers import PVDF, THREE_STATE
from polyfind.strain_response import frame_weights, proper_columns, rotate_columns


def fd(x, evaluate, h=1e-6):
    x = np.array(x, dtype=float)
    result = np.empty(np.asarray(evaluate(x)).shape + x.shape)
    for idx in np.ndindex(x.shape):
        plus, minus = x.copy(), x.copy()
        plus[idx] += h
        minus[idx] -= h
        result[(Ellipsis,) + idx] = (np.asarray(evaluate(plus)) - evaluate(minus)) / (2 * h)
    return result


def rotation():
    Q, _ = np.linalg.qr(np.random.default_rng(41).normal(size=(3, 3)))
    if np.linalg.det(Q) < 0:
        Q[:, 0] *= -1
    return Q


@pytest.fixture
def flux_packer():
    from polyfind.fitting import FITTED_VALENCE

    base = SimpleFF.from_preset("pvdf-dft-valence")
    flux = SimpleFF(**{**{k: v for k, v in base.__dict__.items() if k != "_cache"},
                       "charge_flux": (("C", "H", -0.7, 0.31), ("C", "F", 0.45, -0.22))})
    with FITTED_VALENCE.applied():
        chain = periodic_chain(PVDF, (0, 0), THREE_STATE)
        return pack_mod.CrystalPacker(chain, n_chains=2, valence=base, charge_flux=flux,
                                      coulomb="ewald", polarizable=Polarizable())


def placed(packer, flip=0.0):
    params = np.array([4.6, 8.6, 90.0, 12.0, 43.0, 0.4, flip])
    P, H = packer._place(params[None])
    return params, PlacedCell(to_numpy(P)[0], to_numpy(H)[0])


def test_affine_cell_preserves_fractional_nuclei_for_all_six_strains():
    X = np.random.default_rng(4).normal(size=(7, 3))
    H = np.array([[4., 0., 0.], [1., 6., 0.], [.3, -.7, 3.]])
    baseline = PlacedCell(X, H)
    for k in range(6):
        eps = np.zeros(6)
        eps[k] = .02
        cell = baseline.deformed(np.eye(3) + M.strain_tensor(eps))
        np.testing.assert_allclose(cell.coords @ np.linalg.inv(cell.lattice), X @ np.linalg.inv(H), atol=1e-15)
    X[:] = 0
    assert np.any(baseline.coords)
    with pytest.raises(ValueError):
        baseline.coords[0, 0] = 3


def test_repeat_batches_do_not_confuse_three_scalar_rows_with_one_vector():
    np.testing.assert_array_equal(repeat_rows([2., 3., 4.], 3)[:, 2], [2., 3., 4.])
    np.testing.assert_array_equal(repeat_rows(0., 3, repeat=[2., 3., 4.]), [[2., 3., 4.]] * 3)
    for bad in (True, [1., 2.], [1., 2., np.nan], [1j, 0, 2]):
        with pytest.raises(ValueError):
            repeat_vector(bad)
    for bad in (np.zeros((3, 3)), np.diag([1., 1., -1.]), np.full((3, 3), np.nan)):
        with pytest.raises(ValueError):
            PlacedCell(np.ones((2, 3)), bad)


def test_flux_repeat_and_coordinate_derivatives_include_boundary_images(flux_packer):
    flux = flux_packer._flux
    X = np.asarray(flux_packer.chain.coords, dtype=float)
    repeat = np.array([.2, -.15, flux_packer.chain.c])
    q, dq, dr = flux.charges_and_repeat_grad(X, repeat)
    np.testing.assert_allclose(q, flux.charges(X, repeat=repeat)[0], atol=1e-15)
    np.testing.assert_allclose(fd(X, lambda v: flux.charges(v, repeat=repeat)[0]), dq, atol=3e-9)
    np.testing.assert_allclose(fd(repeat, lambda v: flux.charges(X, repeat=v)[0]), dr, atol=3e-9)
    assert np.max(np.abs(dr[:, :2])) > 1e-4
    q0, d0, c0 = flux.charges_and_grad(X, flux_packer.chain.c)
    qv, dv, rv = flux.charges_and_repeat_grad(X, [0., 0., flux_packer.chain.c])
    np.testing.assert_array_equal(q0, qv)
    np.testing.assert_array_equal(d0, dv)
    np.testing.assert_array_equal(c0, rv[:, 2])
    Q = rotation()
    np.testing.assert_allclose(flux.charges(X @ Q, repeat=repeat @ Q)[0], q, atol=1e-14)


@pytest.mark.parametrize("dipolar", [False, True])
def test_exclusion_repeat_gradients_and_frame_covariance(dipolar):
    rng = np.random.default_rng(11)
    X, q, p = rng.normal(size=(5, 3)), rng.normal(size=5), rng.normal(size=(5, 3)) * .2
    S = rng.uniform(0., 1., (5, 5, 5))
    S = .5 * (S + S[::-1].transpose(0, 2, 1))
    repeat = np.array([.6, -.3, 2.6])
    def evaluate(coords=X, charges=q, vector=repeat, dipoles=p):
        if dipolar:
            return charge_dipole_exclusion(coords, charges, vector, S, dipoles=dipoles, grad=True, charge_grad=True)
        return exclusion_correction(coords, charges, vector, S, grad=True, charge_grad=True)
    result = evaluate()
    np.testing.assert_allclose(fd(X, lambda v: evaluate(coords=v)[0]), result[1], atol=2e-6)
    np.testing.assert_allclose(fd(q, lambda v: evaluate(charges=v)[0]), result[3], atol=2e-6)
    np.testing.assert_allclose(fd(repeat, lambda v: evaluate(vector=v)[0]), result[2], atol=2e-6)
    Q = rotation()
    rotated = evaluate(coords=X @ Q, vector=repeat @ Q, dipoles=p @ Q)
    assert rotated[0] == pytest.approx(result[0], abs=1e-11)
    np.testing.assert_allclose(rotated[1], result[1] @ Q, atol=1e-11)
    np.testing.assert_allclose(rotated[2], result[2] @ Q, atol=1e-11)
    if dipolar:
        np.testing.assert_allclose(fd(p, lambda v: evaluate(dipoles=v)[0]), -result[4], atol=2e-6)
        np.testing.assert_allclose(rotated[4], result[4] @ Q, atol=1e-11)


@pytest.mark.parametrize("flip", [0., 1.])
def test_placed_dipole_reproduces_canonical_and_rotates_as_vector(flux_packer, flip):
    params, cell = placed(flux_packer, flip)
    np.testing.assert_allclose(flux_packer.placed_charges(cell.coords, cell.lattice, flip), flux_packer._q_cell, atol=1e-14)
    moment = sum(flux_packer.placed_dipole(cell.coords, cell.lattice, flip))
    np.testing.assert_allclose(moment, flux_packer.dipole(params[None])[0], atol=1e-12)
    eps = np.array([.002, -.001, .003, .015, -.017, .006])
    cell = cell.deformed(np.eye(3) + M.strain_tensor(eps))
    moment = sum(flux_packer.placed_dipole(cell.coords, cell.lattice, flip))
    Q = rotation()
    transformed = sum(flux_packer.placed_dipole(cell.coords @ Q, cell.lattice @ Q, flip))
    np.testing.assert_allclose(transformed, moment @ Q, atol=1e-10)
    shift = np.array([1.3, -.7, .4])
    translated = sum(flux_packer.placed_dipole(cell.coords + shift, cell.lattice, flip))
    np.testing.assert_allclose(translated, moment, atol=1e-11)


def test_polarization_full_lattice_gradient_includes_tilted_exclusion_images(flux_packer):
    _, cell = placed(flux_packer, 1.)
    cell = cell.deformed(np.array([[1., .01, .02], [.01, 1., -.015], [.02, -.015, 1.]]))
    q = flux_packer.placed_charges(cell.coords, cell.lattice, 1.)
    result = flux_packer._polarize(cell.coords, q, cell.lattice, 1., grad=True, charge_grad=True)
    gP, gH, gq = result[3]
    np.testing.assert_allclose(fd(cell.lattice, lambda H: flux_packer._polarize(cell.coords, q, H, 1.)[0]), gH, atol=2e-6)
    np.testing.assert_allclose(fd(cell.coords, lambda X: flux_packer._polarize(X, q, cell.lattice, 1.)[0]), gP, atol=2e-6)
    np.testing.assert_allclose(fd(q, lambda v: flux_packer._polarize(cell.coords, v, cell.lattice, 1.)[0]), gq, atol=2e-6)


def test_complete_clamped_tensor_fixed_charge_geometric_correction_is_zero():
    path = Path(__file__).resolve().parents[1] / "examples" / "clamped_ion_columns.py"
    spec = importlib.util.spec_from_file_location("clamped_columns_control", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    packer = CrystalPacker(periodic_chain(PVDF, (0, 0), THREE_STATE), n_chains=2)
    params, cell = placed(packer)
    ref = M.Reference(packer, params, packer.chain.c)
    P = sum(packer.placed_dipole(cell.coords, cell.lattice)) / ref.volume * module.ISJ.E_PER_A2_TO_C_PER_M2
    for h in (.0025, .005):
        columns = module.clamped_columns(ref, h)
        assert set(columns) == set(range(6))
        record = module.to_provider_columns(columns, np.eye(3), P)
        for k, label in enumerate(module.PROVIDER_VOIGT):
            np.testing.assert_allclose(record[label]["vanderbilt_proper_C_m2"], np.zeros(3), atol=2e-13)
            np.testing.assert_allclose(columns[k], M.strain_tensor(np.eye(6)[k]) @ P, atol=2e-13)


def test_complete_response_rotates_all_three_indices_not_just_the_vector():
    rng = np.random.default_rng(23)
    response = rng.normal(size=(3, 6))
    L = rotation()
    transformed = rotate_columns(response, L)
    for eps in rng.normal(size=(8, 6)):
        old_tensor = L.T @ M.strain_tensor(eps) @ L
        old_eps = np.array([old_tensor[0, 0], old_tensor[1, 1], old_tensor[2, 2],
                            2*old_tensor[1, 2], 2*old_tensor[0, 2], 2*old_tensor[0, 1]])
        np.testing.assert_allclose(transformed @ eps, L @ response @ old_eps, atol=3e-15)
    P = rng.normal(size=3)
    np.testing.assert_allclose(proper_columns(transformed, L @ P),
                               rotate_columns(proper_columns(response, P), L), atol=2e-15)
    np.testing.assert_allclose(rotate_columns(transformed, L.T), response, atol=2e-15)
    with pytest.raises(ValueError):
        rotate_columns(response[:, :4], L)
    with pytest.raises(ValueError):
        frame_weights(L * 1.001)
