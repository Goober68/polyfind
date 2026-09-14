"""Declared chain/local maps, arbitrary independent-chain energy and dipoles."""
from dataclasses import replace
import numpy as np
import pytest

from polyfind import pack as pack_mod
from polyfind.ewald import EwaldSpec
from polyfind.forcefield import SimpleFF
from polyfind.mechanics import strain_tensor
from polyfind.pack import CrystalPacker,periodic_chain
from polyfind.periodic_geometry import ChainLayout
from polyfind.polarizability import Polarizable
from polyfind.polymers import PVDF,THREE_STATE
from polyfind.topology import build_cell


def model(full=False):
    if full:
        from polyfind.fitting import FITTED_VALENCE
        with FITTED_VALENCE.applied():
            ch = periodic_chain(PVDF,(0,0),THREE_STATE)
            return pack_mod.CrystalPacker(ch,coulomb="ewald",ewald=EwaldSpec(),
                                         valence=SimpleFF.from_preset("pvdf-dft-valence"),
                                         charge_flux=SimpleFF.from_preset("pvdf-dft-valence-flux-born"),
                                         polarizable=Polarizable(),lj_cutoff="force",field=(.02,-.03,.01))
    return CrystalPacker(periodic_chain(PVDF,(0,0,0,0),THREE_STATE),lj_cutoff="force")


def permuted(cell):
    order = np.random.default_rng(41).permutation(cell.n_atoms)
    return replace(cell,elements=list(np.asarray(cell.elements)[order]),coords=cell.coords[order],
                   chain_of=cell.chain_of[order],local_of=cell.local_of[order]),order


@pytest.mark.parametrize("full",[False,True])
def test_declared_noncontiguous_chain_maps_preserve_energy_forces_and_dipole(full):
    pk = model(full)
    p = np.array([4.7,8.5,88.,25.,-40.,.7,1.])
    cell = build_cell(pk.chain,p,nx=2,ny=2,packer=pk)
    # Independently shift/perturb each chain, not just a repeated two-chain motif.
    cell.coords = cell.coords+np.random.default_rng(6).normal(size=cell.coords.shape)*.004
    layout = ChainLayout.from_cell(cell)
    E,gP,gH = pk.chain_energy_and_grad(cell.coords,cell.lattice,layout)
    mu = pk.chain_dipole(cell.coords,cell.lattice,layout)
    other,order = permuted(cell)
    Eo,gPo,gHo = pk.chain_energy_and_grad(other.coords,other.lattice,ChainLayout.from_cell(other))
    assert Eo == E
    np.testing.assert_array_equal(gPo,gP[order])
    np.testing.assert_array_equal(gHo,gH)
    for actual,expected in zip(pk.chain_dipole(other.coords,other.lattice,ChainLayout.from_cell(other)),mu):
        np.testing.assert_array_equal(actual,expected)


@pytest.mark.parametrize("flip",[0.,1.])
def test_replication_preserves_energy_and_dipole_per_repeat(flip):
    pk = model()
    p = np.array([4.7,8.5,88.,25.,-40.,.7,flip])
    baseline = build_cell(pk.chain,p,packer=pk)
    base_layout = ChainLayout.from_cell(baseline)
    E = pk.chain_energy_and_grad(baseline.coords,baseline.lattice,base_layout)[0]
    dipole = pk.chain_dipole(baseline.coords,baseline.lattice,base_layout)[0]
    for nx,ny in ((2,1),(1,2),(2,2)):
        cell = build_cell(pk.chain,p,nx=nx,ny=ny,packer=pk)
        layout = ChainLayout.from_cell(cell)
        assert pk.chain_energy_and_grad(cell.coords,cell.lattice,layout)[0]/(nx*ny) == pytest.approx(E,abs=2e-11)
        np.testing.assert_allclose(pk.chain_dipole(cell.coords,cell.lattice,layout)[0]/(nx*ny),dipole,atol=2e-14)


def test_eight_chain_full_field_model_coordinate_lattice_and_strain_derivatives():
    pk = model(True)
    p = np.array([4.7,8.5,88.,25.,-40.,.7,1.])
    cell = build_cell(pk.chain,p,nx=2,ny=2,packer=pk)
    deformation = np.array([[1.01,.006,-.009],[.006,.99,.011],[-.009,.011,1.02]])
    P = cell.coords@deformation+np.random.default_rng(6).normal(size=cell.coords.shape)*.003
    H = cell.lattice@deformation
    layout = ChainLayout.from_cell(cell)
    E,gP,gH = pk.chain_energy_and_grad(P,H,layout)
    rng = np.random.default_rng(41)
    for i,d in [(int(rng.integers(len(P))),int(rng.integers(3))) for _ in range(12)]:
        plus,minus = P.copy(),P.copy()
        plus[i,d] += 1e-6
        minus[i,d] -= 1e-6
        difference = (pk.chain_energy_and_grad(plus,H,layout)[0]-pk.chain_energy_and_grad(minus,H,layout)[0])/2e-6
        assert difference == pytest.approx(gP[i,d],abs=3e-5)
    for i,d in np.ndindex(3,3):
        plus,minus = H.copy(),H.copy()
        plus[i,d] += 1e-6
        minus[i,d] -= 1e-6
        difference = (pk.chain_energy_and_grad(P,plus,layout)[0]-pk.chain_energy_and_grad(P,minus,layout)[0])/2e-6
        assert difference == pytest.approx(gH[i,d],abs=3e-5)
    for h in (1e-4,5e-5):
        for k in np.eye(6):
            S = strain_tensor(k)
            plus,minus = np.eye(3)+h*S,np.eye(3)-h*S
            derivative = (pk.chain_energy_and_grad(P@plus,H@plus,layout)[0]-pk.chain_energy_and_grad(P@minus,H@minus,layout)[0])/(2*h)
            assert derivative == pytest.approx(np.sum(gP*(P@S))+np.sum(gH*(H@S)),abs=1e-3)


def test_arbitrary_per_chain_reversal_is_a_declared_fact():
    pk = model(True)
    p = np.array([4.7,8.5,88.,25.,-40.,.7,0.])
    cell = build_cell(pk.chain,p,nx=2,ny=2,packer=pk)
    reversal = np.array([True,True,False,True,False,False,True,False])
    for s in np.flatnonzero(reversal):
        sl = slice(s*pk.n,(s+1)*pk.n)
        center = cell.coords[sl].mean(axis=0)
        cell.coords[sl] = (cell.coords[sl]-center)@np.diag([1.,-1.,-1.])+center
    cell.reversed_of = reversal
    layout = ChainLayout.from_cell(cell)
    E,gP,gH = pk.chain_energy_and_grad(cell.coords,cell.lattice,layout)
    assert np.isfinite(E) and np.isfinite(gP).all() and np.isfinite(gH).all()
    for i,d in [(0,1),(pk.n+1,2),(5*pk.n,0)]:
        plus,minus = cell.coords.copy(),cell.coords.copy()
        plus[i,d] += 1e-6
        minus[i,d] -= 1e-6
        derivative = (pk.chain_energy_and_grad(plus,cell.lattice,layout)[0]-pk.chain_energy_and_grad(minus,cell.lattice,layout)[0])/2e-6
        assert derivative == pytest.approx(gP[i,d],abs=3e-5)
    other,order = permuted(cell)
    Eo,gPo,gHo = pk.chain_energy_and_grad(other.coords,other.lattice,ChainLayout.from_cell(other))
    assert Eo == E
    np.testing.assert_array_equal(gPo,gP[order])
    np.testing.assert_array_equal(gHo,gH)


def test_layout_owns_complete_order_and_refuses_corrupted_metadata():
    pk = model()
    p = np.array([4.7,8.5,88.,25.,-40.,.7,0.])
    cell = build_cell(pk.chain,p,packer=pk)
    layout = ChainLayout.from_cell(cell)
    with pytest.raises(ValueError):
        layout.atoms[0,0] = 2
    bad = cell.local_of.copy()
    bad[0] = bad[1]
    with pytest.raises(ValueError):
        ChainLayout.from_cell(replace(cell,local_of=bad))
    with pytest.raises(ValueError):
        ChainLayout.from_cell(replace(cell,chain_of=cell.chain_of.astype(float)))
    with pytest.raises(ValueError):
        ChainLayout(layout.atoms,layout.reversed_of.astype(int),layout.elements)
    with pytest.raises(ValueError):
        ChainLayout(layout.atoms[:,::-1],layout.reversed_of,layout.elements[:1])
    wrong = replace(layout,elements=tuple(["Cl"]*cell.n_atoms))
    with pytest.raises(ValueError):
        pk.chain_energy_and_grad(cell.coords,cell.lattice,wrong)
