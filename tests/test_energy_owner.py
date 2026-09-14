"""Consumer delegation and real nonprimitive torsional Hessian contribution."""
import numpy as np
from dataclasses import replace

from polyfind import phonon as PH
from polyfind.pack import CrystalPacker, periodic_chain
from polyfind.polymers import PVDF, THREE_STATE
from polyfind.torsion_geometry import dihedral_value_gradient
from polyfind.mechanics import _shape_scalar


def test_canonical_and_phonon_consumers_delegate_to_complete_owner(monkeypatch):
    pk = CrystalPacker(periodic_chain(PVDF,(0,0),THREE_STATE))
    p = np.array([4.7,8.5,88.,25.,-40.,.7,1.])
    calls = []
    def observed(P,H,flip,**kw):
        calls.append((P.copy(),H.copy(),flip))
        return 123.4,np.ones_like(P),np.zeros_like(H)
    monkeypatch.setattr(pk,"placed_energy_and_grad",observed)
    assert pk.energy(p[None])[0] == 123.4
    assert pk.energy_and_grad(p)[0] == 123.4
    P,H = PH.placed_coordinates(pk,p)
    assert PH.cell_energy(pk,p,P,H) == 123.4
    assert len(calls) == 3
    for actual in calls:
        np.testing.assert_array_equal(actual[0],P)
        np.testing.assert_array_equal(actual[1],H)
        assert actual[2] == 1.


def test_mechanics_shape_pullback_does_not_add_metadata_torsion_again():
    chain = periodic_chain(PVDF,(0,1,0,2),THREE_STATE)
    other = replace(chain,dihedrals=np.full_like(chain.dihedrals,20.))
    gX = np.random.default_rng(41).normal(size=chain.coords.shape)
    expected = float((gX*chain.coords).sum())+.7*chain.c+3.*chain.rotation_error**2
    assert _shape_scalar(chain,gX,.7,3.) == expected
    assert _shape_scalar(other,gX,.7,3.) == expected


def test_phonon_hessian_contains_stationary_nonprimitive_torsion_curvature():
    chain = periodic_chain(PVDF,(0,0,0,0),THREE_STATE)
    pk = CrystalPacker(chain)
    zero = CrystalPacker(chain,torsion=(0.,0.,0.))
    p = np.array([4.7,8.5,88.,25.,-40.,.7,1.])
    P,H = PH.placed_coordinates(pk,p)
    direction = np.zeros_like(P)
    direction[chain.backbone[0],1] = 1.
    curvatures = []
    for h in (1e-4,5e-5):
        derivative = []
        for owner in (pk,zero):
            plus = PH.cell_energy_and_grad(owner,p,P+h*direction,H)[1]
            minus = PH.cell_energy_and_grad(owner,p,P-h*direction,H)[1]
            derivative.append(np.sum(direction*(plus-minus))/(2*h))
        curvatures.append(derivative[0]-derivative[1])
    # Both Hamiltonians have identical nonbonded terms; only torsion differs.
    expected = []
    for h in (1e-4,5e-5):
        plus = pk._torsion.energy_and_repeat_grad(P[:pk.n]+h*direction[:pk.n],H[2])[1]
        minus = pk._torsion.energy_and_repeat_grad(P[:pk.n]-h*direction[:pk.n],H[2])[1]
        expected.append(np.sum(direction[:pk.n]*(plus-minus))/(2*h))
    np.testing.assert_allclose(curvatures,expected,rtol=1e-9,atol=1e-9)
    assert curvatures[0] > 1.
    np.testing.assert_allclose(curvatures[0],curvatures[1],rtol=1e-7)


def test_complete_nonprimitive_gamma_hessian_torsion_block_matches_analytic_matrix():
    chain = periodic_chain(PVDF,(0,0,0,0),THREE_STATE)
    pk = CrystalPacker(chain)
    zero = CrystalPacker(chain,torsion=(0.,0.,0.))
    p = np.array([4.7,8.5,88.,25.,-40.,.7,1.])
    P,H = PH.placed_coordinates(pk,p)
    jacobian = np.zeros((2*len(pk._torsion.atoms),3*pk.N))
    for s in range(2):
        X = P[s*pk.n:(s+1)*pk.n]
        repeat = H[2]*(-1. if s == 1 else 1.)
        points = X[pk._torsion.atoms]+pk._torsion.images[...,None]*repeat
        phi,dphi = dihedral_value_gradient(points)
        np.testing.assert_allclose(np.cos(phi),-1.,atol=1e-14)
        for k,quad in enumerate(pk._torsion.atoms):
            row = jacobian[s*len(pk._torsion.atoms)+k].reshape(pk.N,3)
            np.add.at(row,quad+s*pk.n,dphi[k])
    V1,V2,V3 = np.tile(pk._torsion.coefficients,(2,1)).T
    expected = jacobian.T@(.5*(V1+4*V2+9*V3)[:,None]*jacobian)
    numerical = []
    for h in (1e-4,5e-5):
        full = PH.gamma_hessian(pk,p,h=h,Pn=P,latn=H).H
        reference = PH.gamma_hessian(zero,p,h=h,Pn=P,latn=H).H
        numerical.append(full-reference)
        np.testing.assert_allclose(numerical[-1],expected,atol=3e-6)
    np.testing.assert_allclose(numerical[0],numerical[1],atol=3e-6)
    assert np.linalg.eigvalsh(expected).min() > -1e-12
    assert np.linalg.eigvalsh(expected).max() > 1.
