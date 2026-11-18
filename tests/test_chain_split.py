"""A cell's chains taken out alone (examples/chain_split.py)."""

import importlib.util
from pathlib import Path

import numpy as np
import pytest

ase = pytest.importorskip("ase")
SPEC = importlib.util.spec_from_file_location("chain_split", Path(__file__).resolve().parents[1] / "examples" / "chain_split.py")
CS = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CS)


def _cell(shift=0.0):
    """Two zigzag -CH2-CF2- chains in a cell whose chain axis leans; the second chain cut at another place."""
    from ase import Atoms

    c = 2.56
    unit = [("C", (0.0, 0.42, 0.0)), ("H", (0.88, 1.05, 0.0)), ("H", (-0.88, 1.05, 0.0)),
            ("C", (0.0, -0.42, c / 2)), ("F", (1.09, -1.20, c / 2)), ("F", (-1.09, -1.20, c / 2))]
    cell = np.array([[9.0, 0.0, 0.0], [0.5, 8.0, 0.0], [0.7, -0.4, c]])
    lean = cell[2] / c
    sym, pos = [], []
    for origin, dz in (((2.0, 2.0, 0.0), 0.0), ((6.5, 6.0, 0.0), shift)):
        for s, p in unit:
            z = p[2] + dz
            sym.append(s)
            pos.append(np.array(origin) + np.array([p[0], p[1], 0.0]) + z * lean)
    a = Atoms(symbols=sym, positions=pos, cell=cell, pbc=True)
    a.wrap()
    return a


def test_the_chains_are_found_whole_and_stood_along_z():
    a = _cell(shift=0.9)
    found = CS.chains_of(a)
    assert len(found) == 2 and all(len(idx) == 6 for idx, _, _ in found)
    for idx, P, T in found:
        assert T == pytest.approx(a.cell.array[2])
        lone = CS.lone_chain(a, idx, P, T)
        assert lone.cell.array == pytest.approx(np.diag([CS.BOX, CS.BOX, 2.56 * np.linalg.norm(a.cell.array[2] / 2.56)]))
        assert lone.positions[:, :2].mean(axis=0) == pytest.approx([CS.BOX / 2, CS.BOX / 2])
        # nothing inside the chain moved: every distance is the cell's own
        d0 = np.linalg.norm(P[:, None] - P[None], axis=-1)
        assert np.linalg.norm(lone.positions[:, None] - lone.positions[None], axis=-1) == pytest.approx(d0, abs=1e-10)
        # and the chain closes on its image at a bond's length
        cc = np.linalg.norm(lone.positions[0] - lone.positions[3] + [[0, 0, k * lone.cell[2, 2]] for k in (-1, 0, 1)], axis=1)
        assert np.sort(cc)[:2] == pytest.approx([np.hypot(0.84, 1.28 * np.linalg.norm(a.cell.array[2]) / 2.56)] * 2, rel=0.2)


def test_equal_chains_have_equal_signatures_wherever_they_were_cut():
    a = _cell(shift=0.9)
    lone = [CS.lone_chain(a, *c) for c in CS.chains_of(a)]
    assert np.abs(CS.signature(lone[0]) - CS.signature(lone[1])).max() < 1e-9
    b = a.copy()
    b.positions[7] += [0.0, 0.0, 0.2]  # a hydrogen of the second chain
    lone = [CS.lone_chain(b, *c) for c in CS.chains_of(b)]
    assert np.abs(CS.signature(lone[0]) - CS.signature(lone[1])).max() > 0.05


def test_polyfinds_lone_chain_is_the_two_chain_cell_with_equal_chains_on_the_same_lattice():
    """The split's chain alone on a 16 A square lattice, from a one-chain packer, against the
    cell's own model holding two copies of that chain on a 16 x 32 lattice: the same sum."""
    import json

    from polyfind import mixedcell as MC
    from polyfind.fitting import PRESET_PARAMETERS

    CA = CS._load("candidate_acceptance")
    SP = CA.SP
    root = Path(__file__).resolve().parents[1]
    pol, ref = json.loads((root / SP.POLISHED).read_text()), json.loads((root / SP.REFERENCES).read_text())
    row = "vdf5-rvdf1-6"
    with PRESET_PARAMETERS["pvdf-dft-valence-belt4"].applied():
        letters, _ = SP.DK.sequence(row)
        polymer = SP.C.make_polymer(row, letters)
        d = SP.stored_cell(row, "TTTT", pol, ref)
        model, x0, ch = SP.sheared_model(polymer, d, check=False)
        assert model.flips[0] == model.flips[1]
        s = CS.polyfind_split(SP, model, x0, ch, tilts=False)
        st = model.state(x0)
        x = np.array(x0)
        x[0], x[1], x[2] = CS.BOX, 2 * CS.BOX, 90.0
        x[model.index("place", 1)] = [st.place[0, 0], st.place[0, 1] + 0.5, st.place[0, 2], st.place[0, 3]]
        x[model.index("internal", 1)] = x[model.index("internal", 0)]
        two = model.physical_energy(x)
        pk1 = SP.pack_mod.CrystalPacker(ch, n_chains=1, **SP.DK.KW)
        one = MC.MixedModel(polymer, "ref", pk1, model.flips[:1], [ch], angle_stiffness=model.angle_stiffness)
        e1 = one.physical_energy(np.concatenate([[CS.BOX, CS.BOX, 90.0, st.lattice[3]], [0.0, 0.0, st.place[0, 2], 0.0], st.tors[0], st.angs[0]]))
    assert two == pytest.approx(2 * e1, abs=1e-8)
    assert s["cell"] == pytest.approx(s["chain"] + s["between"], abs=1e-9) and s["cell"] == pytest.approx(float(d["energy"]), abs=1e-6)
    # the polar all-trans cell holds its chains by about four kcal/mol per monomer under belt4
    assert s["between"] / (2 * ch.n_monomers) == pytest.approx(-4.0, abs=0.1)
