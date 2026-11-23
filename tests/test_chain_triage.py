"""examples/chain_triage.py: the parts that need no DFT run."""
import importlib.util
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]


def _load():
    spec = importlib.util.spec_from_file_location("chain_triage", ROOT / "examples" / "chain_triage.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


CT = _load()


def test_every_chain_names_its_backbone_ends_and_a_state():
    for name, chain in CT.CHAINS.items():
        assert 0 < chain["tail"] < chain["m"], name
        assert chain["states"], name
        for torsions, period in chain["states"].values():
            assert period >= 1 and len(torsions) >= 1


def test_pvdf_periods_are_the_known_ones():
    """All-trans: one monomer, about 2.5 A; TGTG': two monomers, about 2.2 to 2.3 A a monomer."""
    pytest.importorskip("rdkit")
    chain = CT.CHAINS["pvdf"]
    tt = CT.periodic(chain, "tt")
    assert len(tt) == 6 and sorted(tt.get_chemical_symbols()) == ["C", "C", "F", "F", "H", "H"]
    assert 2.4 < tt.cell.lengths()[2] < 2.7
    helix = CT.periodic(chain, "tgtg")
    assert len(helix) == 12
    assert 2.1 < helix.cell.lengths()[2] / 2 < 2.4
    assert helix.cell.lengths()[2] / 2 < tt.cell.lengths()[2]  # the helix is the shorter chain


def test_an_amide_chain_is_cut_whole():
    """Nylon-3's period is one repeat, C3H5NO, with every hydrogen of its heavy atoms."""
    pytest.importorskip("rdkit")
    a = CT.periodic(CT.CHAINS["nylon3"], "tt")
    symbols = a.get_chemical_symbols()
    assert (symbols.count("C"), symbols.count("N"), symbols.count("O"), symbols.count("H")) == (3, 1, 1, 5)
    assert 4.5 < a.cell.lengths()[2] < 5.3
    assert np.allclose(a.positions[:, :2].mean(axis=0), 7.0)


def test_moment_of_two_point_charges(tmp_path):
    """A lump of two electrons 1 A off an ion of charge 2 across the box: a dipole of 2 e A."""
    from ase import Atoms
    from ase.io.cube import write_cube

    cell = np.diag([8.0, 8.0, 3.0])
    n = (40, 40, 6)
    atoms = Atoms("He", positions=[[4.0, 4.0, 1.5]], cell=cell, pbc=True)
    rho = np.zeros(n)
    voxel = np.linalg.det(cell) / np.prod(n) / 0.529177210903**3  # bohr3
    rho[25, 20, :] = 2.0 / n[2] / voxel  # x = 25/40 * 8 = 5.0 A, y = 4.0 A
    cube = tmp_path / "rho.cube"
    with open(cube, "w") as fh:
        write_cube(fh, atoms, data=rho)
    saved = CT.PSEUDO.get("He")
    CT.PSEUDO["He"] = ("none", 4.0, 2)
    try:
        mu, electrons, ions = CT.moment(cube, atoms)
    finally:
        if saved is None:
            del CT.PSEUDO["He"]
    assert electrons == pytest.approx(2.0, abs=1e-6) and ions == 2.0
    assert mu[0] == pytest.approx(-2.0, abs=1e-6)  # electrons 1 A to +x of the ion: the dipole points to -x
    assert mu[1] == pytest.approx(0.0, abs=1e-6)


def test_read_berry_takes_the_last_value_and_its_direction(tmp_path):
    out = tmp_path / "berry.out"
    out.write_text(
        "         P =   0.5385472  (mod   4.8333333)  (e/Omega).bohr\n"
        "         P =   0.0027133  (mod   0.0243513)  e/bohr^2\n"
        "         P =   0.1551305  (mod   1.3922482)  C/m^2\n"
        "     The polarization direction is:  ( 0.00000 , 0.60000 , 0.80000 )\n"
        "     JOB DONE.\n")
    value, modulus, direction = CT.read_berry(out)
    assert value == pytest.approx(0.1551305) and modulus == pytest.approx(1.3922482)
    assert np.allclose(direction, [0.0, 0.6, 0.8])


def test_a_lattice_start_keeps_the_chain_and_its_period():
    from ase import Atoms

    chain = Atoms("CH2", positions=[[7.0, 7.0, 0.0], [7.9, 7.5, 0.0], [6.1, 7.5, 0.0]], cell=[14.0, 14.0, 2.5], pbc=True)
    chain.info.update(triage_chain="x", triage_state="tt", triage_repeats=1, triage_kind="chain")
    a = CT.lattice_start(chain, np.random.default_rng(3))
    assert a.info["triage_kind"] == "crystal" and a.info["triage_repeats"] == 1
    assert np.allclose(a.cell.array[2], [0.0, 0.0, 2.5])
    assert a.get_volume() < chain.get_volume()
    assert CT.closest(a) > 0.0
    assert CT.kmesh(a).split()[2] == "8" and CT.kmesh(chain) == "1 1 8"


def test_a_turned_cell_cancels_what_the_chain_carries_across_itself():
    """Doubled along an edge with the second chain half round: twice the atoms, the same chain axis,
    and a set of point charges whose moment across the chain is zero while the one along it doubles."""
    from ase import Atoms

    a = Atoms("CHF", positions=[[0.0, 0.0, 0.0], [0.9, 0.3, 0.4], [-0.8, 0.6, 1.1]], cell=[[4.0, 0.0, 0.3], [1.0, 4.5, -0.2], [0.0, 0.0, 2.5]], pbc=True)
    a.info.update(triage_chain="x", triage_state="tt", triage_repeats=1, triage_kind="crystal", kmesh="4 5 8")
    q = np.array([0.1, 0.4, -0.5])
    for edge, mesh in ((0, "2 5 8"), (1, "4 3 8")):
        b = CT.turned(a, edge)
        assert len(b) == 6 and b.info["triage_repeats"] == 2 and b.info["kmesh"] == mesh
        assert np.allclose(b.cell.array[edge], 2 * a.cell.array[edge]) and np.allclose(b.cell.array[2], a.cell.array[2])
        first = (q[:, None] * (b.positions[:3] - b.positions[:3].mean(axis=0))).sum(axis=0)
        second = (q[:, None] * (b.positions[3:] - b.positions[3:].mean(axis=0))).sum(axis=0)
        assert np.allclose((first + second)[:2], 0.0) and np.isclose((first + second)[2], 2 * first[2])
        assert np.allclose(b.positions[:3], a.positions)


def test_turning_by_nothing_and_by_half_a_turn_are_the_two_ends():
    from ase import Atoms

    a = Atoms("CHF", positions=[[0.0, 0.0, 0.0], [0.9, 0.3, 0.4], [-0.8, 0.6, 1.1]], cell=[[4.0, 0.0, 0.3], [1.0, 4.5, -0.2], [0.0, 0.0, 2.5]], pbc=True)
    a.info.update(triage_chain="x", triage_state="tt", triage_repeats=1, triage_kind="crystal", kmesh="4 5 8")
    assert np.allclose(CT.turned_by(a, 1, 0.0).positions, a.repeat([1, 2, 1]).positions)
    assert np.allclose(CT.turned_by(a, 1, np.pi).positions, CT.turned(a, 1).positions)
    quarter = CT.turned_by(a, 0, np.pi / 2)
    d0 = a.positions - a.positions.mean(axis=0)
    d1 = quarter.positions[3:] - quarter.positions[3:].mean(axis=0)
    assert np.allclose(np.linalg.norm(d0, axis=1), np.linalg.norm(d1, axis=1))  # a rigid turn
    assert np.allclose(d0[:, 2], d1[:, 2])  # about the chain's own axis
    assert np.allclose((d0[:, :2] * d1[:, :2]).sum(axis=1), 0.0)  # by a quarter


def test_rebased_puts_the_hydrogen_bond_along_the_first_edge_and_keeps_the_crystal():
    """A chain of one N-H and one O whose bond runs along the cell's diagonal: afterwards it runs
    along the first edge, with the same volume and the same atoms in the same places."""
    from ase import Atoms

    cell = np.array([[4.0, 0.0, 0.0], [0.0, 4.0, 0.0], [0.0, 0.0, 3.0]])
    a = Atoms("NHO", positions=[[1.0, 1.0, 0.0], [1.7, 0.3, 0.0], [-1.03, 3.03, 0.0]], cell=cell, pbc=True)
    a.info.update(triage_chain="x", triage_state="tt", triage_repeats=1, triage_kind="crystal", kmesh="5 5 7")
    # the H at (1.7, 0.3) is 1.8 A from the O of the image one cell along x and one back along y, at (2.97, -0.97)
    assert CT.bond_step(a) == (1, -1)
    b = CT.rebased(a)
    assert CT.bond_step(b) == (1, 0)
    assert np.isclose(abs(np.linalg.det(b.cell.array)), abs(np.linalg.det(cell)))
    assert np.allclose(b.positions, a.positions) and np.allclose(b.cell.array[0], cell[0] - cell[1])
    plain = Atoms("CH", positions=[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], cell=cell, pbc=True)
    assert CT.bond_step(plain) is None and CT.rebased(plain) is plain


def test_spun_turns_one_chain_of_a_row_and_leaves_the_rest():
    from ase import Atoms

    row = Atoms("CHCH", positions=[[0.0, 0.0, 0.0], [1.0, 0.2, 0.3], [4.0, 0.0, 0.0], [5.0, 0.2, 0.3]], cell=[8.0, 5.0, 2.5], pbc=True)
    before = row.positions.copy()
    CT.spun(row, 2, 1, np.pi, np.array([0.0, 0.0, 1.0]))
    assert np.allclose(row.positions[:2], before[:2])
    assert np.allclose(row.positions[2:].mean(axis=0), before[2:].mean(axis=0))
    assert np.allclose(row.positions[3] - row.positions[2], [-1.0, -0.2, 0.3])


def test_methylated_swaps_one_n_h_for_n_ch3_and_keeps_the_chains_in_place():
    from ase import Atoms

    row = Atoms("NHCNHC", positions=[[0.0, 0.0, 0.0], [1.02, 0.0, 0.0], [-1.3, 0.0, 0.6], [5.0, 0.0, 0.0], [6.02, 0.0, 0.0], [3.7, 0.0, 0.6]],
                cell=[10.0, 6.0, 4.0], pbc=True)
    out = CT.methylated(row, 3, 1)
    assert len(out) == 9 and out.get_chemical_symbols() == ["N", "H", "C", "N", "C", "C", "H", "H", "H"]
    assert np.allclose(out.positions[:4], row.positions[:4]) and np.allclose(out.positions[5], row.positions[5])
    carbon = out.positions[4]
    assert np.isclose(np.linalg.norm(carbon - out.positions[3]), 1.46)
    for h in out.positions[6:]:
        assert np.isclose(np.linalg.norm(h - carbon), 1.09)
        to_h, to_n = h - carbon, out.positions[3] - carbon
        assert np.isclose(np.degrees(np.arccos(to_h @ to_n / 1.09 / 1.46)), 109.47, atol=0.05)
