"""Geometry-derived torsion forces/stiffness, all atoms and repeat components."""
import dataclasses
import numpy as np
import pytest

from polyfind.forcefield import SimpleFF, dihedral_angles
from polyfind.mechanics import strain_tensor
from polyfind.pack import periodic_chain, chain_torsion, CrystalPacker
from polyfind.periodic_torsion import ChainTorsion
from polyfind.polymers import PVDF, PE, THREE_STATE
from polyfind.torsion_geometry import dihedral_value_gradient, fourier_terms, fourier_derivative


def fd(value, evaluate, h=1e-6):
    value = np.array(value, dtype=float)
    result = np.zeros_like(value)
    for idx in np.ndindex(value.shape):
        plus, minus = value.copy(), value.copy()
        plus[idx] += h
        minus[idx] -= h
        result[idx] = (evaluate(plus)-evaluate(minus))/(2*h)
    return result


@pytest.mark.parametrize("polymer,sequence", [(PVDF, (0, 0)), (PVDF, (0, 1, 0, 2)), (PE, (0,))])
def test_periodic_topology_reproduces_builder_angles_and_coefficients(polymer, sequence):
    chain = periodic_chain(polymer, sequence, THREE_STATE)
    ff = SimpleFF(torsion_by_bond=((1., 2., 3.), (2., 4., 1.)) if polymer.bonds_per_repeat == 2 else None)
    V = ff.torsion_coefficients(polymer, len(chain.dihedrals))
    term = chain_torsion(chain, V)
    phi = term.angles(chain.coords, chain.c)[0]
    np.testing.assert_allclose(np.exp(1j*phi), np.exp(1j*np.deg2rad(chain.dihedrals)), atol=3e-14)
    expected = float(fourier_terms(np.deg2rad(chain.dihedrals), V).sum())
    assert term.energy(chain.coords, chain.c)[0] == pytest.approx(expected, abs=1e-13)
    assert len(term.atoms) == len(chain.dihedrals)
    assert np.any(term.images)
    # Only topology/period length counts matter, not metadata angle values.
    other = dataclasses.replace(chain, dihedrals=np.full_like(chain.dihedrals, 20.))
    other_term = chain_torsion(other, V)
    np.testing.assert_array_equal(term.atoms, other_term.atoms)
    np.testing.assert_array_equal(term.images, other_term.images)
    assert term.energy(chain.coords, chain.c)[0] == other_term.energy(chain.coords, chain.c)[0]


@pytest.mark.parametrize("sequence", [(0, 0, 0, 0), (0, 1, 0, 2)])
def test_periodic_torsion_full_coordinate_repeat_and_strain_derivatives(sequence):
    chain = periodic_chain(PVDF, sequence, THREE_STATE)
    term = chain_torsion(chain, (1.3, -.05, 2.5))
    X = chain.coords + np.random.default_rng(6).normal(size=chain.coords.shape)*.03
    repeat = np.array([.03, -.02, chain.c])
    E, gX, gr = term.energy_and_repeat_grad(X, repeat)
    assert E == pytest.approx(term.energy(X, repeat=repeat)[0], abs=1e-14)
    np.testing.assert_allclose(fd(X, lambda v: term.energy(v, repeat=repeat)[0]), gX, atol=2e-8)
    np.testing.assert_allclose(fd(repeat, lambda v: term.energy(X, repeat=v)[0]), gr, atol=2e-8)
    np.testing.assert_allclose(gX.sum(axis=0), np.zeros(3), atol=2e-14)
    np.testing.assert_allclose(np.cross(X,gX).sum(axis=0)+np.cross(repeat,gr), np.zeros(3), atol=2e-13)
    analytical = [np.sum(gX*(X@strain_tensor(k)))+gr@(repeat@strain_tensor(k)) for k in np.eye(6)]
    for h in (1e-4, 5e-5):
        numerical = []
        for k in np.eye(6):
            energies = []
            for sign in (-1, 1):
                F = np.eye(3)+sign*h*strain_tensor(k)
                energies.append(term.energy(X@F, repeat=repeat@F)[0])
            numerical.append((energies[1]-energies[0])/(2*h))
        np.testing.assert_allclose(numerical, analytical, atol=2e-6)
    assert np.linalg.norm(gX) > .01
    assert np.linalg.norm(gr) > .01


def test_periodic_torsion_covariance_reversal_and_batches():
    chain = periodic_chain(PVDF, (0, 1, 0, 2), THREE_STATE)
    term = chain_torsion(chain, (1.3, -.05, 2.5))
    X = chain.coords + np.random.default_rng(6).normal(size=chain.coords.shape)*.03
    repeat = np.array([.03, -.02, chain.c])
    E, gx, gr = term.energy_and_repeat_grad(X, repeat)
    Q, _ = np.linalg.qr(np.random.default_rng(41).normal(size=(3,3)))
    for frame in (Q, np.diag([1., -1., -1.])):
        Er, gxr, grr = term.energy_and_repeat_grad(X@frame, repeat@frame)
        assert Er == pytest.approx(E, abs=2e-14)
        np.testing.assert_allclose(gxr, gx@frame, atol=3e-14)
        np.testing.assert_allclose(grr, gr@frame, atol=3e-14)
    Xs = np.stack([X, X@Q, X])
    repeats = np.stack([repeat, repeat@Q, repeat])
    np.testing.assert_allclose(term.energy(Xs, repeat=repeats), [E]*3, atol=2e-14)
    scalar = term.energy_and_grad(X, chain.c)
    vector = term.energy_and_repeat_grad(X, [0., 0., chain.c])
    assert scalar[0] == vector[0] and scalar[2] == vector[2][2]
    np.testing.assert_array_equal(scalar[1], vector[1])


def test_trans_state_has_nonzero_cartesian_torsional_stiffness():
    chain = periodic_chain(PVDF, (0, 0, 0, 0), THREE_STATE)
    term = chain_torsion(chain, (1.3, -.05, 2.5))
    X = chain.coords.copy()
    repeat = np.array([0., 0., chain.c])
    E, gradient, _ = term.energy_and_repeat_grad(X, repeat)
    assert E == pytest.approx(0., abs=1e-14)
    assert np.linalg.norm(gradient) < 1e-12
    displacement = np.zeros_like(X)
    displacement[chain.backbone[0], 1] = 1.
    points = X[term.atoms]+term.images[..., None]*repeat
    _, dphi = dihedral_value_gradient(points)
    phase_velocity = (dphi*displacement[term.atoms]).sum(axis=(1, 2))
    V1, V2, V3 = term.coefficients.T
    expected = float((.5*(V1+4*V2+9*V3)*phase_velocity**2).sum())
    stiffness = []
    for h in (1e-4, 5e-5):
        plus = term.energy_and_repeat_grad(X+h*displacement, repeat)[1]
        minus = term.energy_and_repeat_grad(X-h*displacement, repeat)[1]
        stiffness.append(float(np.sum(displacement*(plus-minus))/(2*h)))
    assert stiffness[0] > 1.
    assert stiffness[0] == pytest.approx(expected, rel=1e-7)
    assert stiffness[0] == pytest.approx(stiffness[1], rel=1e-7)


def test_primitive_trans_repeat_cannot_represent_chain_twisting():
    chain = periodic_chain(PVDF, (0, 0), THREE_STATE)
    term = chain_torsion(chain, (1.3, -.05, 2.5))
    assert len(chain.backbone) == 2
    X = chain.coords + np.random.default_rng(6).normal(size=chain.coords.shape)*.03
    repeat = np.array([.03, -.02, chain.c])
    # Two independent backbone sites and one translation always span a plane.
    phi = term.angles(X, repeat=repeat)[0]
    np.testing.assert_allclose(np.cos(phi), -np.ones(2), atol=1e-14)
    E, gx, gr = term.energy_and_repeat_grad(X, repeat)
    assert E == pytest.approx(0., abs=1e-14)
    np.testing.assert_allclose(gx, 0., atol=1e-13)
    np.testing.assert_allclose(gr, 0., atol=1e-13)
    for k in np.eye(6):
        F = np.eye(3)+.01*strain_tensor(k)
        assert term.energy(X@F, repeat=repeat@F)[0] == pytest.approx(0., abs=1e-14)


def test_shared_dihedral_gradient_and_fourier_derivative():
    points = np.random.default_rng(2).normal(size=(4,3))
    phi, gradient = dihedral_value_gradient(points)
    assert np.degrees(phi) == pytest.approx(dihedral_angles(points, [[0,1,2,3]])[0])
    V = np.array([1.3, -.05, 2.5])
    np.testing.assert_allclose(fd(points, lambda p: fourier_terms(dihedral_value_gradient(p)[0], V)),
                               fourier_derivative(phi,V)*gradient, atol=2e-8)
    packer = CrystalPacker(periodic_chain(PVDF, (0,0), THREE_STATE))
    angles = np.array([180., 60., -60., 17.])
    assert packer.torsion_energy(angles) == float(fourier_terms(np.deg2rad(angles), packer.torsion).sum())


def test_torsion_image_representative_and_atom_wrapping_are_gauge_invariant():
    chain = periodic_chain(PVDF, (0,1,0,2), THREE_STATE)
    term = chain_torsion(chain, (1.3, -.05, 2.5))
    X = chain.coords+np.random.default_rng(6).normal(size=chain.coords.shape)*.03
    repeat = np.array([.03, -.02, chain.c])
    E, gx, gr = term.energy_and_repeat_grad(X, repeat)
    shifted = ChainTorsion(term.n_atoms, term.atoms, term.images+3, term.coefficients)
    Es, gxs, grs = shifted.energy_and_repeat_grad(X, repeat)
    assert Es == pytest.approx(E, abs=2e-13)
    np.testing.assert_allclose(gxs, gx, atol=2e-13)
    np.testing.assert_allclose(grs, gr, atol=2e-13)
    wraps = np.random.default_rng(41).integers(-2,3,term.n_atoms)
    wrapped = ChainTorsion(term.n_atoms, term.atoms, term.images-wraps[term.atoms], term.coefficients)
    Ew, gxw, grw = wrapped.energy_and_repeat_grad(X+wraps[:,None]*repeat, repeat)
    assert Ew == pytest.approx(E, abs=2e-13)
    np.testing.assert_allclose(gxw, gx, atol=2e-13)
    # The repeat derivative is at fixed stored coordinates; wrapping changes it.
    np.testing.assert_allclose(grw+(gxw*wraps[:,None]).sum(axis=0), gr, atol=2e-13)


def test_torsion_ownership_refuses_bad_inputs_and_is_immutable():
    chain = periodic_chain(PVDF, (0,0), THREE_STATE)
    term = chain_torsion(chain, (1.3, -.05, 2.5))
    with pytest.raises(ValueError):
        term.atoms[0,0] = 3
    with pytest.raises(ValueError):
        term.energy(chain.coords, chain.c, repeat=[0.,0.,chain.c])
    with pytest.raises(ValueError):
        term.energy_and_repeat_grad(np.stack([chain.coords]*2), [0.,0.,chain.c])
    with pytest.raises(ValueError):
        ChainTorsion(chain.n_atoms, term.atoms.astype(float), term.images, term.coefficients)
    with pytest.raises(ValueError):
        ChainTorsion(chain.n_atoms, term.atoms, term.images, [True,True,False])
    with pytest.raises(ValueError):
        dihedral_value_gradient(np.zeros((4,3)))
    with pytest.raises(ValueError):
        term.energy(np.full_like(chain.coords, np.nan), chain.c)
