"""Full independent-atom/six-strain chart and physically guarded relaxation."""
from dataclasses import replace

import numpy as np
import pytest

from polyfind.cartesian_mechanics import (
    CartesianChart, CartesianRelaxer, CartesianTolerance, RelaxationPhase, TranslationGauge,
)
from polyfind.mechanics import strain_tensor
from polyfind.pack import CellEnergy, CellEnergyTerms, CrystalPacker, periodic_chain
from polyfind.periodic_geometry import ChainLayout
from polyfind.polymers import PVDF, THREE_STATE
from polyfind.topology import Cell, build_cell, check_topology, intended_edges


def model_cell():
    from polyfind import pack as pack_mod
    from polyfind.ewald import EwaldSpec
    from polyfind.fitting import FITTED_VALENCE
    from polyfind.forcefield import SimpleFF
    from polyfind.polarizability import Polarizable
    with FITTED_VALENCE.applied():
        ch = periodic_chain(PVDF,(0,0,0,0),THREE_STATE)
        pk = pack_mod.CrystalPacker(ch,coulomb="ewald",ewald=EwaldSpec(),
                                    valence=SimpleFF.from_preset("pvdf-dft-valence"),
                                    charge_flux=SimpleFF.from_preset("pvdf-dft-valence-flux-born"),
                                    polarizable=Polarizable(),lj_cutoff="force",field=(.02,-.03,.01))
    return pk,build_cell(ch,[4.7,8.5,88.,25.,-40.,.7,1.],nx=2,ny=1,packer=pk)


class ElasticOwner:
    """Independent analytic atomic spring/full metric elasticity, not chemistry."""
    def evaluate_chain_cell(self,P,H,layout):
        delta = P[0]-P[1]-np.array([1.,0.,0.])
        atom = .5*20.*(delta@delta)
        gP = np.array([20.*delta,-20.*delta])
        metric = H.T@H/100.
        mismatch = metric-np.eye(3)
        lattice = .5*1000.*np.sum(mismatch*mismatch)
        gH = 20.*H@mismatch
        terms = CellEnergyTerms(atom+lattice,0.,0.,0.,0.,atom+lattice,0.,0.,1,1)
        return CellEnergy(terms,gP,gH)


def elastic_cell():
    return Cell(["C","C"],np.array([[1.2,.3,.1],[0.,0.,0.]]),
                np.array([[10.1,.2,.1],[.2,9.9,.15],[.1,.15,10.2]]),
                np.array([0,0]),np.array([0,1]),2,np.array([False]))


def test_translation_gauge_is_complete_orthonormal_and_linear_memory():
    masses = np.array([12.,1.,19.,35.])
    gauge = TranslationGauge(masses)
    rng = np.random.default_rng(9)
    x = rng.normal(size=(3,3))
    U = gauge.displacements(x)
    np.testing.assert_allclose(masses@U,0.,atol=3e-14)
    np.testing.assert_allclose(gauge.pullback(U),x,atol=5e-16)
    assert np.sum(U*U) == pytest.approx(np.sum(x*x),abs=2e-15)
    G = rng.normal(size=(4,3))
    assert np.sum(G*U) == pytest.approx(np.sum(gauge.pullback(G)*x),abs=2e-15)
    assert sum(v.size for v in gauge.__dict__.values()) == 2*len(masses)
    extreme = TranslationGauge(masses*1e300)
    np.testing.assert_allclose(extreme.displacements(x),U,atol=1e-15)


def test_full_model_chart_atom_and_six_log_strain_pullbacks_at_two_steps():
    pk,cell = model_cell()
    chart = CartesianChart(cell,strain_scale=3.)
    x = chart.x0+np.random.default_rng(5).normal(size=chart.n_dof)*.003
    x[-6:] = np.array([.01,-.02,.03,.014,-.019,.015])*chart.strain_scale
    geometry,Q,F,eta = chart.decode(x)
    evaluation = pk.evaluate_chain_cell(geometry.coords,geometry.lattice,chart.layout)
    gradient,stress,conjugate = chart.pullback(x,evaluation)
    np.testing.assert_allclose(chart.gauge.masses@(Q-chart.reference.coords),0.,atol=1e-13)
    assert np.linalg.det(geometry.lattice) > 0.
    for h in (1e-5,5e-6):
        for k in [0,4,chart.atomic_dof-1,*range(chart.atomic_dof,chart.n_dof)]:
            delta = np.zeros(chart.n_dof)
            delta[k] = h
            energies = []
            for value in (x+delta,x-delta):
                placed = chart.decode(value)[0]
                energies.append(pk.evaluate_chain_cell(placed.coords,placed.lattice,chart.layout).terms.total)
            assert (energies[0]-energies[1])/(2*h) == pytest.approx(gradient[k],abs=4e-5)
    # Physical stress follows CURRENT affine nuclei and all lattice rows,
    # not the log-coordinate force (the two differ at finite strain).
    volume = np.linalg.det(geometry.lattice)
    from polyfind.cartesian_mechanics import KCAL_MOL_A3_TO_MPA
    for h in (1e-5,5e-6):
        for k in np.eye(6):
            S = strain_tensor(k)
            values = []
            for stretch in (np.eye(3)+h*S,np.eye(3)-h*S):
                values.append(pk.evaluate_chain_cell(geometry.coords@stretch,geometry.lattice@stretch,chart.layout).terms.total)
            assert (values[0]-values[1])/(2*h)/volume*KCAL_MOL_A3_TO_MPA == pytest.approx(stress@k,abs=.003)
    assert np.max(abs(stress-conjugate)) > .01


def test_constraints_are_exact_log_components_not_inferred_lengths():
    cell = elastic_cell()
    fixed = np.array([.01,.02,.03,.04,.05,.06])
    mask = np.array([True,True,False,True,False,True])
    chart = CartesianChart(cell,mask,fixed)
    x = chart.x0.copy()
    x[-4:] += .2
    geometry,Q,F,eta = chart.decode(x)
    np.testing.assert_array_equal(eta[~mask],fixed[~mask])
    np.testing.assert_allclose(geometry.lattice,cell.lattice@F,atol=1e-14)
    assert chart.n_dof == 3*(cell.n_atoms-1)+mask.sum()
    # Large compressive log strain remains nonsingular without a physics bound.
    chart = CartesianChart(cell,np.zeros(6,dtype=bool),[-10.,0.,0.,0.,0.,0.])
    assert np.linalg.det(chart.decode(chart.x0)[0].lattice) > 0.


def test_small_optimizer_gradients_do_not_hide_atom_or_shear_stress():
    cell = elastic_cell()
    chart = CartesianChart(cell,strain_scale=1e9)
    placed = chart.decode(chart.x0)[0]
    evaluation = ElasticOwner().evaluate_chain_cell(placed.coords,placed.lattice,chart.layout)
    report = CartesianTolerance().assess(chart,chart.x0,evaluation)
    assert not report.converged
    assert report.max_atomic_force_kcal_mol_A > 1.
    assert abs(report.stress_MPa[3:]).max() > 100.
    assert abs(chart.pullback(chart.x0,evaluation)[0][-6:]).max() < 1e-6
    clamped = CartesianChart(cell,np.zeros(6,dtype=bool))
    report = CartesianTolerance(atomic_force_kcal_mol_A=100.).assess(clamped,clamped.x0,evaluation)
    assert report.converged
    assert abs(report.stress_MPa[3:]).max() > 100.
    assert report.max_selected_physical_stress_MPa == report.max_free_log_stress_MPa == 0.


def test_full_cell_relaxes_all_atoms_normals_and_shears_with_terminal_guards():
    chart = CartesianChart(elastic_cell())
    observed = []
    relaxer = CartesianRelaxer(ElasticOwner(),chart,CartesianTolerance(1e-5,.01),observed.append)
    result = relaxer.run(maxiter=300)
    assert result.checkpoint.phase == RelaxationPhase.CONVERGED
    assert result.checkpoint.stationarity.converged
    assert result.checkpoint.stationarity.max_atomic_force_kcal_mol_A <= 1e-5
    assert result.checkpoint.stationarity.max_selected_physical_stress_MPa <= .01
    np.testing.assert_allclose(result.geometry.lattice,10.*np.eye(3),atol=2e-7)
    np.testing.assert_allclose(result.geometry.coords[0]-result.geometry.coords[1],[1.,0.,0.],atol=5e-7)
    assert observed[0].phase == RelaxationPhase.RUNNING
    assert observed[-1].phase == RelaxationPhase.CONVERGED
    assert observed[-1].accepted_steps > 0
    assert result.evaluations >= result.checkpoint.accepted_steps
    with pytest.raises(RuntimeError,match="illegal"):
        relaxer.run()


def test_iteration_limit_and_same_chart_recovery_do_not_fabricate_convergence():
    chart = CartesianChart(elastic_cell())
    result = CartesianRelaxer(ElasticOwner(),chart).run(maxiter=0)
    assert result.checkpoint.phase == RelaxationPhase.UNCONVERGED
    assert not result.checkpoint.stationarity.converged
    # Explicit accepted checkpoint with the SAME chart/model, not a new helix.
    recovered = CartesianRelaxer(ElasticOwner(),chart).run(maxiter=300,initial=result.checkpoint.parameters)
    assert recovered.checkpoint.phase == RelaxationPhase.CONVERGED


def test_failed_publication_is_failed_not_terminal_success():
    chart = CartesianChart(elastic_cell())
    def failing(record):
        if record.phase != RelaxationPhase.RUNNING:
            raise OSError("checkpoint publication unavailable")
    driver = CartesianRelaxer(ElasticOwner(),chart,checkpoint_sink=failing)
    with pytest.raises(OSError,match="publication"):
        driver.run(maxiter=0)
    assert driver.phase == driver.checkpoint.phase == RelaxationPhase.FAILED
    assert "publication" in driver.failure_diagnosis
    with pytest.raises(RuntimeError,match="illegal"):
        driver.run()


def test_topology_expectations_follow_noncontiguous_source_maps():
    chain = periodic_chain(PVDF,(0,0,0,0),THREE_STATE)
    cell = build_cell(chain,[4.7,8.5,88.,25.,-40.,.7,1.],nx=2,ny=1)
    order = np.random.default_rng(41).permutation(cell.n_atoms)
    other = replace(cell,elements=list(np.asarray(cell.elements)[order]),coords=cell.coords[order],
                    chain_of=cell.chain_of[order],local_of=cell.local_of[order])
    source_to_other = np.argsort(order)
    expected = set()
    from polyfind.topology import _key
    for i,j,image in intended_edges(cell,chain):
        expected.add(_key(int(source_to_other[i]),int(source_to_other[j]),image))
    assert intended_edges(other,chain) == expected
    first,second = check_topology(cell,chain),check_topology(other,chain)
    assert first.ok == second.ok
    assert first.components == second.components and first.missing == second.missing and first.extra == second.extra
    assert first.safe_scale_window == second.safe_scale_window
    # The cell's count alone does not declare a doubled chemical graph.
    with pytest.raises(ValueError,match="atoms per repeat"):
        intended_edges(cell,PVDF)
    with pytest.raises(ValueError,match="complete chemical"):
        intended_edges(cell,replace(chain,backbone=chain.backbone[:-1]))
    with pytest.raises(ValueError,match="local atom order"):
        intended_edges(cell,replace(chain,backbone=chain.backbone[::-1]))


def test_nonfinite_invalid_metadata_and_controls_refused():
    cell = elastic_cell()
    for mask in ([1]*6,[True]*5):
        with pytest.raises(ValueError,match="mask"):
            CartesianChart(cell,mask)
    for value in (True,0.,-1.,np.inf,1j,"bad"):
        with pytest.raises(ValueError,match="scale"):
            CartesianChart(cell,strain_scale=value)
        with pytest.raises(ValueError,match="tolerances"):
            CartesianTolerance(stress_MPa=value)
    chart = CartesianChart(cell)
    with pytest.raises(ValueError,match="parameters"):
        chart.decode(chart.x0*np.nan)
    with pytest.raises(ValueError,match="iteration"):
        CartesianRelaxer(ElasticOwner(),chart).run(maxiter=-1)


def test_pilot_codec_owns_whole_schema_and_integrity():
    import importlib.util
    import json
    from pathlib import Path
    path = Path(__file__).resolve().parents[1]/"examples/cartesian_relaxation.py"
    spec = importlib.util.spec_from_file_location("cartesian_pilot_codec",path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    codec = module.PilotRecord
    payload = dict(schema=codec.SCHEMA,phase="unconverged",quantitatively_valid=False)
    raw = codec.encode(payload)
    assert codec.decode(raw) == payload
    mutated = json.loads(raw)
    mutated["payload"]["phase"] = "converged"
    with pytest.raises(ValueError,match="integrity"):
        codec.decode(codec.canonical(mutated))
    with pytest.raises(ValueError,match="schema"):
        codec.encode(dict(payload,schema="wrong"))
    with pytest.raises(ValueError,match="whole"):
        codec.decode(codec.canonical(dict(json.loads(raw),unknown=1)))
