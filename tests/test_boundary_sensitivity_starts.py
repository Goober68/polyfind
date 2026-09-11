from __future__ import annotations

import importlib.util
import json
import tempfile
from pathlib import Path

import numpy as np


SCRIPT = Path(__file__).resolve().parents[1] / "examples" / "boundary_sensitivity_starts.py"
SPEC = importlib.util.spec_from_file_location("boundary_sensitivity_starts", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
starts = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(starts)


def test_length_ladder_holds_defect_count_and_orientation() -> None:
    atom_counts = []
    for monomers in starts.CHAIN_LENGTHS:
        structure, metadata = starts.build_start("vdcn", monomers)
        atom_counts.append(structure.n_atoms)
        assert metadata["defect_count"] == 1
        assert metadata["defect_monomer_index_0idx"] == monomers // 2
        assert metadata["sequence"].count("vdcn") == 1
        assert len(metadata["backbone_atom_ids_0idx"]) == 2 * monomers
        assert metadata["connected_component_sizes"] == [structure.n_atoms]
        endpoints = structure.coords[structure.backbone[-1]] - structure.coords[structure.backbone[0]]
        assert abs(endpoints[0]) < 1.0e-10
        assert abs(endpoints[1]) < 1.0e-10
        assert endpoints[2] > 0.0
    assert atom_counts == sorted(atom_counts)
    assert len(set(atom_counts)) == len(atom_counts)


def test_materialized_inputs_are_deterministic_and_self_describing() -> None:
    with tempfile.TemporaryDirectory() as directory:
        output = Path(directory)
        first = starts.materialize(output)
        first_payloads = {
            record["path"]: (output / record["path"]).read_bytes()
            for record in first["records"]
        }
        second = starts.materialize(output)
        assert first == second
        assert len(first["records"]) == 9
        assert json.loads((output / "manifest.json").read_text()) == first
        for record in second["records"]:
            payload = (output / record["path"]).read_bytes()
            assert payload == first_payloads[record["path"]]
            assert starts.digest_bytes(payload) == record["sha256"]
            assert record["connected_component_sizes"] == [record["atom_count"]]


def test_pvdf_has_no_defect_and_an_keeps_registered_stereochemistry() -> None:
    _, pvdf = starts.build_start("pvdf", 5)
    _, an = starts.build_start("an", 5)
    assert pvdf["defect_count"] == 0
    assert pvdf["defect_monomer_index_0idx"] is None
    assert an["defect_count"] == 1
    assert "isotactic enantiomer" in an["stereochemistry"]
    assert np.isclose(an["net_illustrative_point_charge_e"], 0.0, atol=1.0e-12)
