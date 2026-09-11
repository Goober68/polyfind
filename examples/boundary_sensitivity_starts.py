"""Build finite all-trans size-ladder inputs for Sarco boundary tests.

The chemistry owner emits one neutral, H-capped chain for each requested size.
AN and VDCN keep exactly one central defect while VDF host units are added
symmetrically, so chain length and defect count are not conflated.

Run with ``PYTHONPATH=src python examples/boundary_sensitivity_starts.py``.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

from polyfind.chain import Structure, build_chain
from polyfind.polymers import AN, PVDF, VDCN, VDF_UNIT, copolymer, monomer_of


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "deliverables" / "boundary_sensitivity"
CHAIN_LENGTHS = (5, 7, 9)
CHEMISTRIES = {"pvdf": PVDF, "an": AN, "vdcn": VDCN}
AXES = {
    "x": "transverse, completing the right-handed frame",
    "y": "transverse, toward the central substituted unit's non-H pendants",
    "z": "chain end-to-end",
}


def digest_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def chain_polymer(chemistry: str, monomers: int):
    if chemistry not in CHEMISTRIES:
        raise KeyError(f"unsupported chemistry {chemistry!r}")
    if monomers < 3 or monomers % 2 == 0:
        raise ValueError("chain length must be odd and at least three monomers")
    if chemistry == "pvdf":
        units = [VDF_UNIT] * monomers
        defect_index = None
    else:
        defect_index = monomers // 2
        defect = monomer_of(CHEMISTRIES[chemistry])
        units = [VDF_UNIT] * monomers
        units[defect_index] = defect
    polymer = copolymer(
        name=f"vdf-{chemistry}-finite-{monomers}",
        units=units,
        bond_length=PVDF.bond_length,
    )
    return polymer, defect_index


def orient_chain(structure: Structure, central_monomer: int) -> Structure:
    coords = np.asarray(structure.coords, dtype=float)
    backbone = np.asarray(structure.backbone, dtype=int)
    z = coords[backbone[-1]] - coords[backbone[0]]
    z /= np.linalg.norm(z)
    substituted = int(backbone[2 * central_monomer + 1])
    non_h = [
        atom
        for atom in structure.subs_of[substituted]
        if structure.elements[atom] != "H"
    ]
    if not non_h:
        raise RuntimeError("central substituted unit has no non-H pendant direction")
    y = coords[non_h].mean(axis=0) - coords[substituted]
    y -= z * np.dot(y, z)
    y /= np.linalg.norm(y)
    x = np.cross(y, z)
    x /= np.linalg.norm(x)
    basis = np.column_stack((x, y, z))
    structure.coords = (coords - coords.mean(axis=0)) @ basis
    return structure


def connected_components(atom_count: int, bonds: list[tuple[int, int]]) -> list[list[int]]:
    adjacency = [[] for _ in range(atom_count)]
    for first, second in bonds:
        adjacency[first].append(second)
        adjacency[second].append(first)
    unseen = set(range(atom_count))
    components = []
    while unseen:
        stack = [unseen.pop()]
        component = []
        while stack:
            atom = stack.pop()
            component.append(atom)
            for neighbor in adjacency[atom]:
                if neighbor in unseen:
                    unseen.remove(neighbor)
                    stack.append(neighbor)
        components.append(sorted(component))
    return sorted(components, key=lambda values: (len(values), values))


def build_start(chemistry: str, monomers: int) -> tuple[Structure, dict]:
    polymer, defect_index = chain_polymer(chemistry, monomers)
    dihedral_count = 2 * monomers - 3
    structure = build_chain(polymer, [180.0] * dihedral_count, cap=True)
    structure = orient_chain(structure, monomers // 2)
    bonds = [(int(first), int(second)) for first, second in structure.bonds]
    components = connected_components(structure.n_atoms, bonds)
    if len(components) != 1:
        raise RuntimeError(f"generated {chemistry}/{monomers} chain is disconnected")
    metadata = {
        "chemistry": chemistry,
        "host": "vdf",
        "monomer_count": monomers,
        "defect_count": 0 if defect_index is None else 1,
        "defect_monomer_index_0idx": defect_index,
        "sequence": [unit.name for unit in polymer.sequence],
        "stereochemistry": (
            "Polyfind polymer definition; AN pendant order selects its registered "
            "isotactic enantiomer; PVDF and VDCN are achiral"
        ),
        "conformation": "all-trans, every constructed backbone dihedral 180 degrees",
        "boundary": "finite neutral H-capped single chain in vacuum",
        "backbone_bond_length_A": polymer.bond_length,
        "atom_count": structure.n_atoms,
        "backbone_atom_ids_0idx": [int(value) for value in structure.backbone],
        "bonds_0idx": [list(pair) for pair in bonds],
        "connected_component_sizes": [len(component) for component in components],
        "net_illustrative_point_charge_e": float(np.sum(structure.charges)),
        "axes": AXES,
    }
    return structure, metadata


def materialize(output: Path = OUT) -> dict:
    output.mkdir(parents=True, exist_ok=True)
    records = []
    for chemistry in CHEMISTRIES:
        for monomers in CHAIN_LENGTHS:
            structure, metadata = build_start(chemistry, monomers)
            filename = f"{chemistry}_vdfhost_{monomers}mer_alltrans.xyz"
            comment = json.dumps(
                {
                    "chemistry": chemistry,
                    "monomers": monomers,
                    "defect_count": metadata["defect_count"],
                    "axes": "x=transverse,y=central-pendant,z=chain",
                },
                separators=(",", ":"),
            )
            payload = structure.to_xyz(comment=comment).encode("ascii")
            (output / filename).write_bytes(payload)
            records.append(
                {
                    "path": filename,
                    "sha256": digest_bytes(payload),
                    **metadata,
                }
            )
    manifest = {
        "schema_version": 1,
        "purpose": "Sarco stage-2 finite-chain boundary and size sensitivity inputs",
        "generator": "examples/boundary_sensitivity_starts.py",
        "length_ladder_monomers": list(CHAIN_LENGTHS),
        "length_invariant": (
            "AN/VDCN defect count remains one while symmetric VDF host units are added"
        ),
        "scope_limit": (
            "CNEPO is excluded because its bridging epoxide changes backbone topology "
            "and is not representable by Polyfind's pendant-fragment model"
        ),
        "records": records,
    }
    manifest_path = output / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return manifest


if __name__ == "__main__":
    result = materialize()
    print(json.dumps({"output": str(OUT), "structures": len(result["records"])}, indent=2))
