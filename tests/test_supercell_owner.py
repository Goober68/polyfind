"""Supercell model ownership, actual geometry, field and image diagnostics."""
from dataclasses import FrozenInstanceError, replace

import numpy as np
import pytest

from polyfind import backend as bk, pack as pack_mod
from polyfind.ewald import EwaldSpec
from polyfind.forcefield import SimpleFF
from polyfind.pack import EV_TO_KCAL, CrystalPacker, periodic_chain
from polyfind.periodic_geometry import ChainLayout
from polyfind.polarizability import Polarizable
from polyfind.polymers import PVDF, THREE_STATE
from polyfind.supercell import SupercellEnergy, cell_dipole
from polyfind.topology import build_cell


def make_model(full):
    if not full:
        return CrystalPacker(periodic_chain(PVDF,(0,0,0,0),THREE_STATE),
                             valence=SimpleFF(bond_terms=(("*",300.,1.5),),
                                              angle_terms=(("*",60.,110.),)),
                             lj_cutoff="force",field=(.02,-.03,.01))
    from polyfind.fitting import FITTED_VALENCE
    with FITTED_VALENCE.applied():
        chain = periodic_chain(PVDF,(0,0),THREE_STATE)
        return pack_mod.CrystalPacker(chain,coulomb="ewald",ewald=EwaldSpec(),
                                     valence=SimpleFF.from_preset("pvdf-dft-valence"),
                                     charge_flux=SimpleFF.from_preset("pvdf-dft-valence-flux-born"),
                                     polarizable=Polarizable(),lj_cutoff="force",field=(.02,-.03,.01))


def make_cell(pk):
    cell = build_cell(pk.chain,[4.7,8.5,88.,25.,-40.,.7,1.],nx=2,ny=1,packer=pk)
    F = np.array([[1.01,.006,-.009],[.006,.99,.011],[-.009,.011,1.02]])
    return replace(cell,coords=cell.coords@F+np.random.default_rng(6).normal(size=cell.coords.shape)*.004,
                   lattice=cell.lattice@F)


@pytest.mark.parametrize("full",[False,True])
def test_complete_terms_force_scatter_and_dipoles_on_independent_chains(full):
    pk = make_model(full)
    cell = make_cell(pk)
    order = np.random.default_rng(41).permutation(cell.n_atoms)
    shuffled = replace(cell,elements=list(np.asarray(cell.elements)[order]),coords=cell.coords[order],
                       chain_of=cell.chain_of[order],local_of=cell.local_of[order])
    evaluator = SupercellEnergy(pk)
    expected = pk.evaluate_chain_cell(cell.coords,cell.lattice,ChainLayout.from_cell(cell))
    actual = evaluator.evaluate(shuffled)
    assert actual.terms == expected.terms
    np.testing.assert_array_equal(actual.grad_coords,expected.grad_coords[order])
    np.testing.assert_array_equal(actual.grad_lattice,expected.grad_lattice)
    E,gP,gH = evaluator.energy_and_grad(shuffled)
    assert E == actual.terms.total == evaluator.energy(shuffled) == evaluator.terms(shuffled).total
    assert evaluator.energy_per_monomer(shuffled) == E/(cell.n_chains*pk.chain.n_monomers)
    np.testing.assert_array_equal(gP,actual.grad_coords)
    np.testing.assert_array_equal(gH,actual.grad_lattice)
    for actual_mu,expected_mu in zip(evaluator.dipole(shuffled),pk.chain_dipole(cell.coords,cell.lattice,ChainLayout.from_cell(cell))):
        np.testing.assert_array_equal(actual_mu,expected_mu)
    np.testing.assert_allclose(cell_dipole(shuffled,pk.chain.charges),cell_dipole(cell,pk.chain.charges),atol=2e-14)
    t = actual.terms
    assert t.total == pytest.approx(t.pair+t.torsion+t.ewald+t.exclusion+t.valence+t.field,abs=2e-11)
    assert t.valence > 0.
    assert t.field != 0.
    # Independent per-chain term reconstruction at the actual nuclei/repeat.
    ordered,q,states = pk._chain_charge_state(cell.coords,cell.lattice,ChainLayout.from_cell(cell))
    torsion = valence = 0.
    for sl,sign,_ in states:
        torsion += float(pk._torsion.energy_and_repeat_grad(ordered.coords[sl],sign*ordered.lattice[2])[0])
        valence += float(pk._valence.energy_and_repeat_grad(ordered.coords[sl],sign*ordered.lattice[2],n_atoms=pk.n)[0])
    assert t.torsion == torsion
    assert t.valence == valence
    assert t.field == -EV_TO_KCAL*float((q@ordered.coords)@pk.field)
    if not full:
        assert t.torsion > 0.  # Four-site repeat really twists; metadata stays TT.
        assert pk.e_torsion == pytest.approx(0.,abs=1e-12)
    else:
        assert t.ewald != 0. and t.exclusion != 0.
    with pytest.raises(FrozenInstanceError):
        actual.terms.total = 0.
    with pytest.raises(ValueError):
        actual.grad_coords[0,0] = 0.


@pytest.mark.parametrize("full",[False,True])
def test_supercell_field_model_atom_and_full_lattice_derivatives(full):
    evaluator = SupercellEnergy(make_model(full))
    cell = make_cell(evaluator.pk)
    _,gP,gH = evaluator.energy_and_grad(cell)
    for h in (1e-5,5e-6):
        for i,d in ((0,0),(cell.n_per_chain+1,2),(cell.n_atoms-1,1)):
            delta = np.zeros_like(cell.coords)
            delta[i,d] = h
            fd = (evaluator.energy(replace(cell,coords=cell.coords+delta))
                  - evaluator.energy(replace(cell,coords=cell.coords-delta)))/(2*h)
            assert fd == pytest.approx(gP[i,d],abs=3e-5)
        for i,d in np.ndindex(3,3):
            delta = np.zeros((3,3))
            delta[i,d] = h
            fd = (evaluator.energy(replace(cell,lattice=cell.lattice+delta))
                  - evaluator.energy(replace(cell,lattice=cell.lattice-delta)))/(2*h)
            assert fd == pytest.approx(gH[i,d],abs=3e-5)


def test_adapter_routes_every_energy_observable_to_complete_owner(monkeypatch):
    pk = make_model(False)
    cell = make_cell(pk)
    evaluator = SupercellEnergy(pk)
    owner = pk.evaluate_chain_cell
    calls = []
    def traced(*args,**kwargs):
        calls.append(kwargs.get("images"))
        return owner(*args,**kwargs)
    monkeypatch.setattr(pk,"evaluate_chain_cell",traced)
    for action in (evaluator.evaluate,evaluator.energy,evaluator.terms,evaluator.energy_per_monomer,
                   evaluator.energy_and_grad,evaluator.pair_energy,evaluator.column_correction):
        before = len(calls)
        action(cell)
        assert len(calls) == before+1
    evaluator.energy_is_converged(cell,extra=1)
    assert len(calls) == 9
    assert calls[-1].shape[1] == calls[-2].shape[1] == 3


def test_pair_diagnostic_uses_actual_geometry_flux_and_effective_images():
    pk = make_model(True)
    cell = make_cell(pk)
    evaluator = SupercellEnergy(pk)
    ordered,q,_ = pk._chain_charge_state(cell.coords,cell.lattice,ChainLayout.from_cell(cell))
    tables = pk._pair_charge_tables(q,q)
    A = np.tile(bk.to_numpy(pk._A_nn),(cell.n_chains,cell.n_chains))
    B = np.tile(bk.to_numpy(pk._B_nn),(cell.n_chains,cell.n_chains))
    raw = 0.
    for D,_ in ordered.pair_image_chunks(pk.rc):
        r2 = (D*D).sum(axis=-1)
        v = pk._pair_energy_and_dv(r2,A,B,*tables)[0]
        raw += .5*float(v.sum())
    result = evaluator.evaluate(cell)
    assert raw+.5*result.terms.pair_correction == pytest.approx(result.terms.pair,abs=2e-9)
    assert result.terms.pair_correction < -1e4
    assert evaluator.pair_energy(cell,evaluator.images(cell)) == result.terms.pair
    assert evaluator.energy_is_converged(cell,extra=2) < 2e-11


def test_image_grid_validation_and_topology_ownership():
    pk = make_model(False)
    cell = make_cell(pk)
    evaluator = SupercellEnergy(pk)
    for extra in (-1,True,.5):
        with pytest.raises(ValueError,match="shells"):
            evaluator.images(cell,extra)
    for images in ([[0,0,.5]],[[0,0,0],[0,0,0]],[],[[0,0,np.nan]],
                   np.array([[0,0,np.iinfo(np.int64).min]],dtype=np.int64)):
        with pytest.raises(ValueError,match="image representatives"):
            evaluator.terms(cell,images)
    other = periodic_chain(PVDF,(0,0,0,0),THREE_STATE)
    with pytest.raises(ValueError,match="topology belongs"):
        SupercellEnergy(pk,other)
    assert SupercellEnergy(pk,pk.chain).chain is pk.chain
    pk.update_chain(other)
    assert evaluator.chain is other
    invalid = replace(cell,local_of=np.zeros(cell.n_atoms,dtype=int))
    with pytest.raises(ValueError,match="exactly once"):
        evaluator.energy(invalid)


@pytest.mark.parametrize("full",[False,True])
def test_field_energy_derivative_is_actual_total_dipole(full):
    pk = make_model(full)
    cell = make_cell(pk)
    evaluator = SupercellEnergy(pk)
    permanent,induced = evaluator.dipole(cell)
    field = np.array(pk.field,copy=True)
    try:
        for h in (1e-5,5e-6):
            for d in range(3):
                delta = np.zeros(3)
                delta[d] = h
                pk.set_field(field+delta)
                plus = evaluator.energy(cell)
                pk.set_field(field-delta)
                minus = evaluator.energy(cell)
                assert (plus-minus)/(2*h) == pytest.approx(-EV_TO_KCAL*(permanent+induced)[d],abs=3e-6)
    finally:
        pk.set_field(field)
