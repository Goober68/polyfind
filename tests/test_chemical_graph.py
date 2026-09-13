import importlib.util
from dataclasses import FrozenInstanceError
from pathlib import Path

import pytest

from polyfind.chemical_graph import ChemicalBondGraph
from polyfind.chain import build_chain,build_chain_batch
from polyfind.polymers import AN,VDCN,PVDF,methoxy,nitrile,Pendant


def test_graph_owns_complete_immutable_ordered_bonds():
    edges=[[1,0,1],[1,2,3]]
    graph=ChemicalBondGraph(3,edges)
    edges[1][2]=1
    assert graph.edges==((0,1,1),(1,2,3))
    exposed=graph.bonds;exposed.clear()
    identity=graph.identity();identity['bond_orders_0idx'].clear()
    assert graph.bonds==[(0,1),(1,2)]
    with pytest.raises(FrozenInstanceError):graph.edges=()


@pytest.mark.parametrize('count,edges',[(0,()),(True,()),(2,((0,1),)),
    (2,((0,1,1),(1,0,3))),(2,((0,2,1),)),(2,((0,0,1),)),
    (2,((False,1,1),)),(2,((0,1,True),)),(2,((0,1,float('nan')),)),
    (2,((0,1,1.5),))])
def test_invalid_or_incomplete_chemical_facts_reject(count,edges):
    with pytest.raises(ValueError):ChemicalBondGraph(count,edges)


def test_pendant_owns_single_defaults_and_declares_nitrile_triple():
    cn=nitrile();assert cn.chemical_graph.edges==((0,1,3),)
    ether=methoxy();assert all(order==1 for _,_,order in ether.chemical_graph.edges)
    with pytest.raises(ValueError,match='every pendant bond'):
        Pendant(cn.atoms,cn.bonds,bond_orders=())


def test_pendant_snapshots_mutable_constructor_collections():
    cn=nitrile();atoms=list(cn.atoms);bonds=[[0,1]];orders=[3]
    pendant=Pendant(atoms,bonds,bond_orders=orders)
    atoms.clear();bonds[0].clear();bonds.append([0,0]);orders.clear()
    assert pendant.chemical_graph==cn.chemical_graph
    assert len(pendant.atoms)==2 and pendant.bonds==((0,1),)


@pytest.mark.parametrize('polymer,triples',[(PVDF,0),(AN,5),(VDCN,10)])
def test_one_builder_propagates_chemical_orders_to_scalar_and_batch(polymer,triples):
    source=build_chain(polymer,[180.]*7)
    assert sum(order==3 for _,_,order in source.chemical_graph.edges)==triples
    for i,j,order in source.chemical_graph.edges:
        if order==3:assert {source.elements[i],source.elements[j]}=={'C','N'}
    assert all(order==1 for i,j,order in source.chemical_graph.edges
               if i in source.backbone and j in source.backbone)
    template,_=build_chain_batch(polymer,[[180.]*7,[60.]*7])
    assert template.chemical_graph==source.chemical_graph


def test_new_boundary_metadata_leaves_all_archived_xyz_bytes_unchanged(tmp_path):
    root=Path(__file__).resolve().parents[1]
    spec=importlib.util.spec_from_file_location('graph_starts',root/'examples/boundary_sensitivity_starts.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    manifest=module.materialize(tmp_path);assert manifest['schema_version']==2
    for record in manifest['records']:
        old=root/'deliverables/boundary_sensitivity'/record['path']
        assert (tmp_path/record['path']).read_bytes().replace(b'\r\n',b'\n')==old.read_bytes().replace(b'\r\n',b'\n')
        graph=record['chemical_graph']
        assert graph['atom_count']==record['atom_count']
        assert [edge[:2] for edge in graph['bond_orders_0idx']]==record['bonds_0idx']
        expected={'pvdf':0,'an':1,'vdcn':2}[record['chemistry']]
        assert sum(edge[2]==3 for edge in graph['bond_orders_0idx'])==expected
