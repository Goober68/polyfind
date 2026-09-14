"""Periodic bonds/angles retain vector images and full six-strain derivatives."""
import numpy as np
import pytest

from polyfind.forcefield import SimpleFF
from polyfind.mechanics import strain_tensor
from polyfind.pack import chain_valence, periodic_chain
from polyfind.polymers import PVDF, THREE_STATE


def valence_case(sequence):
    import polyfind.fitting  # registers the fitted valence preset

    chain = periodic_chain(PVDF, sequence, THREE_STATE)
    ff = SimpleFF.from_preset("pvdf-dft-valence")
    return chain, chain_valence(chain, ff.bond_table(), ff.angle_table())


@pytest.mark.parametrize("sequence",[(0,0),(0,1,0,2)])
@pytest.mark.parametrize("component",["bond","angle"])
def test_component_gradients_share_complete_valence_and_full_repeat_geometry(sequence,component):
    chain,valence = valence_case(sequence)
    F = np.array([[1.01,.006,-.007],[.006,.99,.008],[-.007,.008,1.003]])
    X = np.asarray(chain.coords)@F
    repeat = np.array([0.,0.,chain.c])@F
    evaluate = getattr(valence,component+"_energy_and_repeat_grad")
    E,gX,gr = evaluate(X,repeat)
    for h in (1e-6,5e-7):
        np.testing.assert_allclose(fd(X,lambda v:evaluate(v,repeat)[0],h=h),gX,atol=4e-7)
        np.testing.assert_allclose(fd(repeat,lambda v:evaluate(X,v)[0],h=h),gr,atol=4e-7)
    # Component sums are independently scattered; the COMPLETE path retains
    # its original bond-then-angle scatter order, so gradients agree to roundoff.
    bE,bG,bR = valence.bond_energy_and_repeat_grad(X,repeat)
    aE,aG,aR = valence.angle_energy_and_repeat_grad(X,repeat)
    total,G,R = valence.energy_and_repeat_grad(X,repeat)
    assert total == bE+aE
    np.testing.assert_allclose(G,bG+aG,atol=2e-14,rtol=0.)
    np.testing.assert_allclose(R,bR+aR,atol=2e-14,rtol=0.)
    np.testing.assert_allclose(np.cross(X,gX).sum(axis=0)+np.cross(repeat,gr),0.,atol=3e-13)


def test_terminal_fluorine_bend_force_is_transverse_to_its_bond():
    chain,valence = valence_case((0,0,0,0))
    X = chain.coords+np.random.default_rng(7).normal(size=chain.coords.shape)*.003
    repeat = np.array([.01,-.02,chain.c])
    _,gA,_ = valence.angle_energy_and_repeat_grad(X,repeat)
    lengths = valence.bond_lengths(X,repeat=repeat)[0]
    _,gB,_ = valence.bond_energy_and_repeat_grad(X,repeat)
    for i,j,s,k,r0,r in zip(valence.bond_i,valence.bond_j,valence.bond_s,valence.bond_k,valence.bond_r0,lengths):
        if chain.elements[j] != "F":
            continue
        direction = (X[j]+s*repeat-X[i])/r
        assert gA[j]@direction == pytest.approx(0.,abs=3e-13)
        assert gB[j]@direction == pytest.approx(k*(r-r0),abs=3e-13)


def fd(value, evaluate, h=1e-6):
    value = np.array(value, dtype=float)
    result = np.zeros_like(value)
    for idx in np.ndindex(value.shape):
        plus, minus = value.copy(), value.copy()
        plus[idx] += h
        minus[idx] -= h
        result[idx] = (evaluate(plus) - evaluate(minus)) / (2*h)
    return result


@pytest.mark.parametrize("sequence", [(0, 0), (0, 1, 0, 2)])
def test_valence_full_repeat_and_atom_derivatives(sequence):
    chain, valence = valence_case(sequence)
    assert valence.n_bonds and valence.n_angles
    assert np.any(valence.bond_s) and (np.any(valence.ang_si) or np.any(valence.ang_sk))
    F = np.array([[1.01, .006, -.007], [.006, .99, .008], [-.007, .008, 1.003]])
    X = np.asarray(chain.coords) @ F
    repeat = np.array([0., 0., chain.c]) @ F
    E, gX, gr = valence.energy_and_repeat_grad(X, repeat)
    lengths = valence.bond_lengths(X, repeat=repeat)[0]
    angles = valence.angle_values(X, repeat=repeat)[0]
    metric_energy = .5*np.sum(valence.bond_k*(lengths-valence.bond_r0)**2) + .5*np.sum(valence.ang_kk*(angles-valence.ang_t0)**2)
    assert metric_energy == pytest.approx(E, abs=2e-14)
    assert E == pytest.approx(valence.energy(X, repeat=repeat)[0], abs=2e-14)
    np.testing.assert_allclose(fd(X, lambda v: valence.energy(v, repeat=repeat)[0]), gX, atol=3e-7)
    np.testing.assert_allclose(fd(repeat, lambda v: valence.energy(X, repeat=v)[0]), gr, atol=3e-7)
    assert np.linalg.norm(gr[:2]) > 1e-3
    np.testing.assert_allclose(gX.sum(axis=0), np.zeros(3), atol=2e-14)
    # Rotational invariance includes the periodic-image torque.
    np.testing.assert_allclose(np.cross(X, gX).sum(axis=0) + np.cross(repeat, gr), np.zeros(3), atol=3e-13)


@pytest.mark.parametrize("sequence", [(0, 0), (0, 1, 0, 2)])
def test_all_six_affine_strain_derivatives_include_repeat_motion(sequence):
    chain, valence = valence_case(sequence)
    F0 = np.array([[1.01, .006, -.007], [.006, .99, .008], [-.007, .008, 1.003]])
    X = np.asarray(chain.coords) @ F0
    repeat = np.array([0., 0., chain.c]) @ F0
    _, gX, gr = valence.energy_and_repeat_grad(X, repeat)
    analytic = []
    for component in np.eye(6):
        basis = strain_tensor(component)
        analytic.append(float(np.sum(gX * (X @ basis)) + gr @ (repeat @ basis)))
    for h in (1e-4, 5e-5):
        numerical = []
        for component in np.eye(6):
            energies = []
            for sign in (-1, 1):
                F = np.eye(3) + sign * h * strain_tensor(component)
                energies.append(valence.energy(X @ F, repeat=repeat @ F)[0])
            numerical.append((energies[1] - energies[0]) / (2*h))
        np.testing.assert_allclose(numerical, analytic, atol=2e-5, rtol=2e-5)
    assert np.linalg.norm(np.asarray(analytic)[3:5]) > 1e-2


def test_repeat_batch_semantics_and_axial_gradient_pullback():
    chain, valence = valence_case((0, 0))
    X = np.asarray(chain.coords)
    repeats = np.array([[.1, .2, chain.c], [.2, -.1, chain.c*1.01], [-.3, .1, chain.c*.99]])
    batch = valence.energy(np.stack([X]*3), repeat=repeats)
    scalar_rows = valence.energy(np.stack([X]*3), repeats[:, 2])
    for k in range(3):
        assert batch[k] == pytest.approx(valence.energy(X, repeat=repeats[k])[0], abs=1e-14)
        assert scalar_rows[k] == pytest.approx(valence.energy(X, repeats[k, 2])[0], abs=1e-14)
    E, gx, gc = valence.energy_and_grad(X, chain.c)
    Ev, gv, gr = valence.energy_and_repeat_grad(X, [0., 0., chain.c])
    assert E == Ev and gc == gr[2]
    np.testing.assert_array_equal(gx, gv)
    with pytest.raises(ValueError):
        valence.energy(X, repeat=[1., 2.])
    with pytest.raises(ValueError):
        valence.energy(X, chain.c, repeat=[0., 0., chain.c])


def test_valence_frame_covariance_and_reversed_chain_images():
    chain, valence = valence_case((0, 1, 0, 2))
    X = np.asarray(chain.coords) + np.random.default_rng(3).normal(size=chain.coords.shape)*.003
    repeat = np.array([.02, -.03, chain.c])
    E, gX, gr = valence.energy_and_repeat_grad(X, repeat)
    Q, _ = np.linalg.qr(np.random.default_rng(41).normal(size=(3, 3)))
    for frame in (Q, np.diag([1., -1., -1.])):
        Er, gx_r, gr_r = valence.energy_and_repeat_grad(X @ frame, repeat @ frame)
        assert Er == pytest.approx(E, abs=3e-13)
        np.testing.assert_allclose(gx_r, gX @ frame, atol=6e-13)
        np.testing.assert_allclose(gr_r, gr @ frame, atol=6e-13)
