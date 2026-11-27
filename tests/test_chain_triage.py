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


def test_a_stretch_is_one_inside_none_outside_and_goes_round_the_ring():
    s = CT.stretch_of([0.0, 0.2, 0.3, 0.5, 0.7, 0.8, 0.95], 0.2, 0.7, 0.2)
    assert np.allclose(s, [0.0, 0.5, 1.0, 1.0, 0.5, 0.0, 0.0])
    round_the_end = CT.stretch_of([0.0, 0.05, 0.3, 0.5, 0.8, 0.9, 0.99], 0.8, 1.1, 0.1)
    assert np.allclose(round_the_end, [1.0, 1.0, 0.0, 0.0, 0.5, 1.0, 1.0])
    ramp = CT.stretch_of(np.linspace(0.1, 0.3, 21), 0.2, 0.7, 0.2)
    assert np.all(np.diff(ramp) >= 0) and ramp[0] == 0.0 and ramp[-1] == 1.0


def _zigzag():
    """A chain of N-H, C=O and CH2 along z, one to a cell 6 A square across: a repeat of 4 A.

    Its methylene is 1.4 A from the next nitrogen (a single bond) and its carbonyl 1.2 A long."""
    from ase import Atoms

    unit = [("N", [0.4, 0.0, 0.0]), ("H", [1.4, 0.0, 0.0]), ("C", [-0.4, 0.0, 1.25]), ("O", [-1.6, 0.0, 1.25]),
            ("C", [0.4, 0.0, 2.6]), ("H", [1.0, 0.9, 2.6]), ("H", [1.0, -0.9, 2.6])]
    a = Atoms([s for s, _ in unit], positions=[p for _, p in unit], cell=[6.0, 6.0, 4.0], pbc=True)
    a.info.update(triage_chain="x", triage_state="tt", triage_repeats=1, triage_kind="crystal", kmesh="4 4 6")
    return a


def test_a_twisted_stretch_is_half_round_and_the_rest_is_where_it_was():
    a = _zigzag()
    row = a.repeat([1, 1, 4])
    before = row.positions.copy()
    along = np.array([0.0, 0.0, 1.0])
    # from the methylene of the first repeat to that of the third, each twist half a repeat wide
    h = CT.hinge(a, along)
    assert h == pytest.approx(2.6 / 4.0)
    CT.twisted(row, slice(0, len(row)), along, h / 4, (2 + h) / 4, 0.5 / 4)
    centre = before.mean(axis=0)
    z = before[:, 2]
    inside = (z > 2.6 + 1.01) & (z < 10.6 - 1.01)  # each twist is half a repeat wide, 1 A each side of its methylene
    outside = (z < 2.6 - 1.01) | (z > 10.6 + 1.01)
    assert inside.sum() >= 8 and outside.sum() >= 8
    assert np.allclose(row.positions[outside], before[outside])
    assert np.allclose(row.positions[inside][:, :2] - centre[:2], -(before[inside][:, :2] - centre[:2]))
    assert np.allclose(row.positions[:, 2], before[:, 2])  # nothing moves along the chain
    assert np.allclose(np.linalg.norm(row.positions[:, :2] - centre[:2], axis=1), np.linalg.norm(before[:, :2] - centre[:2], axis=1))
    assert CT.facing(row, slice(0, len(row)), along) == pytest.approx([0.0, 180.0, 180.0, 0.0], abs=1e-6)


def test_a_kink_row_has_its_walls_and_its_part_turned_chain():
    a = _zigzag()
    n, long = len(a), 4
    along = np.array([0.0, 0.0, 1.0])
    reference = np.array([-1.0, 0.0, 0.0])  # the C=O of the chain as made

    def faces(row):
        return [[abs(round(v)) for v in CT.facing(row, np.where(row.arrays["chain"] == c)[0], along, reference)] for c in range(4)]

    rest = CT.kink_row(a, 4, long, long)
    assert len(rest) == 4 * long * n and rest.info["triage_repeats"] == 16 and rest.info["kmesh"] == "1 4 2"
    assert faces(rest) == [[0] * 4, [0] * 4, [180] * 4, [180] * 4]
    assert faces(CT.kink_row(a, 4, long, 0)) == [[0] * 4, [0] * 4, [0] * 4, [180] * 4]
    part = CT.kink_row(a, 4, long, 2, 0.5)
    assert faces(part) == [[0] * 4, [0] * 4, [0, 180, 180, 0], [180] * 4]
    assert np.allclose(part.positions[:2 * n * long], rest.positions[:2 * n * long])  # the chains that do not move
    assert np.allclose(part.positions[3 * n * long:], rest.positions[3 * n * long:])
    with pytest.raises(ValueError):
        CT.kink_row(a, 4, long, 1, 1.5)
    # the other wall: the last chain of the turned half is the one that moves, and the row at rest is the same cell
    assert np.allclose(CT.kink_row(a, 4, long, long, side=2).positions, rest.positions)
    assert faces(CT.kink_row(a, 4, long, 0, side=2)) == [[0] * 4, [0] * 4, [180] * 4, [0] * 4]
    assert faces(CT.kink_row(a, 4, long, 2, 0.5, side=2)) == [[0] * 4, [0] * 4, [180] * 4, [0, 180, 180, 0]]


def test_a_row_two_sheets_thick_moves_the_chain_of_its_first_sheet_only():
    """Stacked twice, the row at rest is the one-sheet row and its copy a sheet up; a part-turned
    state has the part turned in the first sheet and the chain above it as it is at rest."""
    a = _zigzag()
    n, long = len(a), 4
    along = np.array([0.0, 0.0, 1.0])
    reference = np.array([-1.0, 0.0, 0.0])
    one = CT.kink_row(a, 4, long, long)
    rest = CT.kink_row(a, 4, long, long, sheets=2)
    assert len(rest) == 2 * len(one) and rest.info["triage_repeats"] == 32 and rest.info["kmesh"] == "1 2 2"
    assert np.allclose(rest.cell.array[1], 2 * one.cell.array[1]) and np.allclose(rest.cell.array[[0, 2]], one.cell.array[[0, 2]])

    def sheet(row, k):
        return row.positions[row.arrays["sheet"] == k], row.arrays["chain"][row.arrays["sheet"] == k]

    for k in (0, 1):
        pos, chains = sheet(rest, k)
        assert np.allclose(pos, one.positions + k * one.cell.array[1]) and np.array_equal(chains, one.arrays["chain"])
    part = CT.kink_row(a, 4, long, 2, 0.5, sheets=2)
    mine = CT.moving_chain(part, 4)
    assert len(mine) == n * long and set(part.arrays["sheet"][mine]) == {0} and set(part.arrays["chain"][mine]) == {2}
    assert [abs(round(v)) for v in CT.facing(part, mine, along, reference)] == [0, 180, 180, 0]
    above = np.where((part.arrays["chain"] == 2) & (part.arrays["sheet"] == 1))[0]
    assert [abs(round(v)) for v in CT.facing(part, above, along, reference)] == [180] * 4
    still = np.ones(len(part), dtype=bool)
    still[mine] = False
    assert np.allclose(part.positions[still], rest.positions[still])
    # one sheet thick, the moving chain is the whole of its chain, array or none
    assert np.array_equal(CT.moving_chain(one, 4), np.where(one.arrays["chain"] == 2)[0])
    bare = one.copy()
    del bare.arrays["chain"]
    assert np.array_equal(CT.moving_chain(bare, 4), np.where(one.arrays["chain"] == 2)[0])
    assert np.array_equal(CT.moving_chain(one, 4, side=2), np.where(one.arrays["chain"] == 3)[0])


def test_an_etherified_chain_loses_a_hydrogen_and_keeps_its_carbonyls():
    """The N-H of the second repeat of the moving chain made O: an atom fewer, the same four C=O, and
    the chains after it still found by the cell's own array."""
    a = _zigzag()
    n, long = len(a), 4
    along = np.array([0.0, 0.0, 1.0])
    plain = CT.kink_row(a, 4, long, 2, 0.5)
    made = CT.kink_row(a, 4, long, 2, 0.5, ether=[(1, 0)])
    assert len(made) == len(plain) - 1
    assert made.get_chemical_symbols().count("O") == plain.get_chemical_symbols().count("O") + 1
    assert made.get_chemical_symbols().count("N") == plain.get_chemical_symbols().count("N") - 1
    assert [int(v) for v in np.bincount(made.arrays["chain"])] == [n * long, n * long, n * long - 1, n * long]
    mover = np.where(made.arrays["chain"] == 2)[0]
    reference = np.array([-1.0, 0.0, 0.0])
    assert [abs(round(v)) for v in CT.facing(made, mover, along, reference)] == [0, 180, 180, 0]
    last = np.where(made.arrays["chain"] == 3)[0]
    assert np.allclose(made.positions[last], plain.positions[3 * n * long:])
    with pytest.raises(ValueError):
        CT.kink_row(a, 4, long, 2, 0.5, ether=[(1, 1)])  # this toy chain has one nitrogen a repeat


def test_nearest_reads_a_collapse_across_the_cell_face():
    from ase import Atoms

    a = Atoms("HH", positions=[[0.1, 0.0, 0.0], [3.9, 0.0, 0.0]], cell=[4.0, 5.0, 5.0], pbc=True)
    assert CT.nearest(a) == pytest.approx(0.2)
    b = Atoms("HH", positions=[[1.0, 0.0, 0.0], [3.0, 0.0, 0.0]], cell=[4.0, 5.0, 5.0], pbc=True)
    assert CT.nearest(b) > 1.5


def test_a_held_carbonyl_is_pulled_round_the_axis_to_the_facing_asked_for():
    """The spring's energy is half the stiffness times the angle off, squared; its forces are that energy's slope."""
    from ase import Atoms
    from ase.calculators.calculator import Calculator, all_changes
    from ase.optimize import BFGS

    class Flat(Calculator):
        implemented_properties = ["energy", "forces"]

        def calculate(self, atoms=None, properties=("energy", "forces"), system_changes=all_changes):
            super().calculate(atoms, properties, system_changes)
            self.results = {"energy": 0.0, "forces": np.zeros((len(self.atoms), 3))}

    along, reference = np.array([0.0, 0.0, 1.0]), np.array([1.0, 0.0, 0.3])  # what the reference has along the chain does not count
    at = np.radians(50.0)
    pair = Atoms("CO", positions=[[5.0, 5.0, 5.0], [5.0 + 1.2 * np.cos(at), 5.0 + 1.2 * np.sin(at), 5.4]], cell=[10.0] * 3, pbc=True)
    pair.calc = CT.held(Flat(), 1, 0, along, reference, 30.0, stiff=40.0)
    assert pair.get_potential_energy() == pytest.approx(0.5 * 40.0 * np.radians(20.0) ** 2)
    forces, slope = pair.get_forces(), np.zeros((2, 3))
    for i in range(2):
        for k in range(3):
            e = []
            for d in (1e-5, -1e-5):
                moved = pair.copy()
                moved.positions[i, k] += d
                moved.calc = CT.held(Flat(), 1, 0, along, reference, 30.0, stiff=40.0)
                e.append(moved.get_potential_energy())
            slope[i, k] = -(e[0] - e[1]) / 2e-5
    assert slope == pytest.approx(forces, abs=1e-6)
    assert pair.get_forces().sum(axis=0) == pytest.approx(np.zeros(3), abs=1e-12)
    # a target the other side of half round is reached the short way
    pair.calc = CT.held(Flat(), 1, 0, along, reference, -170.0, stiff=40.0)
    pair.positions[1] = [5.0 + 1.2 * np.cos(np.radians(175.0)), 5.0 + 1.2 * np.sin(np.radians(175.0)), 5.4]
    assert pair.get_potential_energy() == pytest.approx(0.5 * 40.0 * np.radians(15.0) ** 2)
    BFGS(pair, logfile=None).run(fmax=1e-4, steps=200)
    assert CT.facing(pair, slice(0, 2), along, reference)[0] == pytest.approx(-170.0, abs=0.05)


def test_the_atoms_that_turn_with_a_urea_carbonyl_are_its_six():
    """Carbon, oxygen, two nitrogens and their two hydrogens; the methylene beside it stays."""
    from ase import Atoms

    #        O
    #        |
    #  C(H2)-N(H)-C-N(H)-C(H2), along z
    z = np.array([0.0, 0.0, 1.0])
    pos = {"Ca": [0, 0, 0.0], "Ha1": [0.9, 0.5, 0.0], "Ha2": [-0.9, 0.5, 0.0], "N1": [0, 0, 1.45], "H1": [-1.0, 0, 1.45], "C": [0, 0, 2.8],
           "O": [1.23, 0, 2.8], "N2": [0, 0, 4.15], "H2": [-1.0, 0, 4.15], "Cb": [0, 0, 5.6], "Hb1": [0.9, 0.5, 5.6], "Hb2": [-0.9, 0.5, 5.6]}
    chain = Atoms("CHHNHCONHCHH", positions=[np.array(v, dtype=float) + 5.0 for v in pos.values()], cell=[12.0] * 3, pbc=True)
    assert CT.carbonyls(chain, slice(0, len(chain))) == [(6, 5)]
    assert sorted(CT.urea_of(chain, slice(0, len(chain)), 5)) == [3, 4, 5, 6, 7, 8]
    assert CT.facing(chain, slice(0, len(chain)), z, np.array([0.0, 1.0, 0.0])) == pytest.approx([-90.0])


def test_a_methyl_is_found_with_the_carbon_it_is_on_and_what_turns_with_it():
    """The methyl's carbon and the backbone carbon; with them turn the hydrogens on both, and nothing of the nitrogens."""
    from ase import Atoms

    #  N - C(H)(CH3) - N(H) along z; the methyl straight up (y), the carbon's own hydrogen along x
    z = np.array([0.0, 0.0, 1.0])
    pos = {"N1": [0, 0, -1.45], "Ca": [0, 0, 0], "Ha": [1.09, 0, 0], "Cm": [0, 1.53, 0], "Hm1": [1.03, 1.89, 0], "Hm2": [-0.51, 1.89, 0.89],
           "Hm3": [-0.51, 1.89, -0.89], "N2": [0, 0, 1.45], "H2": [-1.0, 0, 1.45]}
    chain = Atoms("NCHCHHHNH", positions=[np.array(v, dtype=float) + 5.0 for v in pos.values()], cell=[12.0] * 3, pbc=True)
    part = slice(0, len(chain))
    assert CT.methyls(chain, part) == [(3, 1)]
    assert sorted(CT.side_of(chain, part, 1)) == [1, 2, 3, 4, 5, 6]
    assert CT.about(chain, CT.methyls(chain, part), z, np.array([1.0, 0.0, 0.0])) == pytest.approx([90.0])
    assert CT.about(chain, [(2, 1)], z, np.array([1.0, 0.0, 0.4])) == pytest.approx([0.0])  # the reference's part along the chain does not count
    assert CT.carbonyls(chain, part) == [] and CT.facing(chain, part, z) == []


def test_a_spread_cell_has_its_sheets_further_apart_and_its_atoms_where_they_were():
    from ase import Atoms

    cell = np.array([[4.6, 0.0, 0.0], [1.7, 3.3, 0.2], [0.1, 0.0, 4.7]])
    a = Atoms("CO", positions=[[1.0, 1.0, 1.0], [2.2, 1.0, 1.0]], cell=cell, pbc=True)
    a.info["triage_repeats"] = 1
    b = CT.spread_cell(a, 1.2)
    assert CT.sheet_gap(b) == pytest.approx(1.2 * CT.sheet_gap(a))
    assert b.get_volume() == pytest.approx(1.2 * a.get_volume())
    assert b.positions == pytest.approx(a.positions)
    assert b.cell.array[0] == pytest.approx(cell[0]) and b.cell.array[2] == pytest.approx(cell[2])
    # what the second edge has along the sheet is kept: the sheets move apart and do not slide
    normal = np.cross(cell[0], cell[2]) / np.linalg.norm(np.cross(cell[0], cell[2]))
    assert (b.cell.array[1] - cell[1]) == pytest.approx(0.2 * (cell[1] @ normal) * normal)


def test_a_slid_cell_has_its_sheets_as_far_apart_and_moved_along_the_chain():
    from ase import Atoms

    cell = np.array([[4.6, 0.0, 0.0], [1.7, 3.3, 0.2], [0.1, 0.0, 4.7]])
    a = Atoms("CO", positions=[[1.0, 1.0, 1.0], [2.2, 1.0, 1.0]], cell=cell, pbc=True)
    a.info["triage_repeats"] = 1
    b = CT.slid_cell(a, 0.5)
    assert CT.sheet_gap(b) == pytest.approx(CT.sheet_gap(a)) and b.get_volume() == pytest.approx(a.get_volume())
    assert b.positions == pytest.approx(a.positions)
    assert b.cell.array[0] == pytest.approx(cell[0]) and b.cell.array[2] == pytest.approx(cell[2])
    assert (b.cell.array[1] - cell[1]) == pytest.approx(0.5 * cell[2])
    # a whole repeat along is the same crystal
    c = CT.slid_cell(a, 1.0)
    assert np.linalg.solve(cell.T, c.cell.array[1]) == pytest.approx([0.0, 1.0, 1.0])


def test_a_cell_told_with_another_second_edge_is_the_same_crystal():
    from ase import Atoms

    cell = np.array([[4.6, 0.0, 0.0], [1.7, 3.3, 0.2], [0.1, 0.0, 4.7]])
    a = Atoms("CO", positions=[[1.0, 1.0, 1.0], [2.2, 1.0, 1.0]], cell=cell, pbc=True)
    a.info["triage_repeats"] = 1
    b = CT.edged_cell(a, -1)
    assert b.cell.array[1] == pytest.approx(cell[1] - cell[0]) and b.cell.array[0] == pytest.approx(cell[0])
    assert b.positions == pytest.approx(a.positions) and b.get_volume() == pytest.approx(a.get_volume())
    assert CT.sheet_gap(b) == pytest.approx(CT.sheet_gap(a))
    assert CT.edged_cell(a, 0).cell.array == pytest.approx(cell)


def test_a_row_of_domains_turns_each_chain_whole_by_its_own_angle():
    from ase import Atoms

    cell = np.array([[4.6, 0.0, 0.0], [-2.9, 3.6, 0.0], [0.0, 0.0, 4.7]])
    # a stand-in chain: C=O across the chain (along the row) with an N-H behind it
    a = Atoms("COHN", positions=[[1.0, 1.0, 1.0], [2.2, 1.0, 1.0], [-0.8, 1.0, 1.0], [0.2, 1.0, 1.0]], cell=cell, pbc=True)
    a.info.update(triage_repeats=1, kmesh="6 6 6", triage_chain="x", triage_state="tt")
    row = CT.domain_row(a, 3, 2, [0, 0, 90])
    assert len(row) == 24 and list(row.arrays["chain"]) == [0] * 8 + [1] * 8 + [2] * 8
    assert row.info["triage_repeats"] == 6 and row.info["kmesh"] == "2 6 3" and row.info["triage_kink_mover"] == 1
    assert row.cell.array == pytest.approx(cell * np.array([[3.0], [1.0], [2.0]]))
    z = np.array([0.0, 0.0, 1.0])
    for k, angle in enumerate((0.0, 0.0, 90.0)):
        part = np.where(row.arrays["chain"] == k)[0]
        assert CT.facing(row, part, z, np.array([1.0, 0.0, 0.0])) == pytest.approx([angle, angle])
        # a chain turns about its own axis: its middle stays where the untouched row has it
        assert row.positions[part].mean(axis=0) == pytest.approx(a.positions.mean(axis=0) + k * cell[0] + 0.5 * cell[2])
    with pytest.raises(ValueError):
        CT.domain_row(a, 3, 2, [0, 0])
    # the other hand: a hydrogen off the plane of the carbonyl and the axis goes to the other face, and the carbonyl faces as it did
    b = a.copy()
    b.info.update(a.info)
    b.positions[2] += [0.0, 0.5, 0.0]
    plain, other = CT.domain_row(b, 2, 1, ["0", "0"]), CT.domain_row(b, 2, 1, ["0", "m0"])
    assert other.positions[:4] == pytest.approx(plain.positions[:4])
    off = lambda q: q.positions[6, 1] - q.positions[4, 1]
    assert off(plain) == pytest.approx(0.5) and off(other) == pytest.approx(-0.5)
    assert CT.facing(other, np.arange(4, 8), z, np.array([1.0, 0.0, 0.0])) == pytest.approx([0.0])
    turned = CT.domain_row(b, 2, 1, ["0", "m90"])
    assert CT.facing(turned, np.arange(4, 8), z, np.array([1.0, 0.0, 0.0])) == pytest.approx([90.0])
    # the first two chains are bonded along the row, O of one to H of the next; the turned one is bonded to neither
    assert CT.joined(row, within=2.0) == {(1, 0): 2}


def test_a_mirrored_chain_faces_the_other_way_about_the_edge_and_has_the_other_hand():
    from ase import Atoms

    z = np.array([0.0, 0.0, 1.0])
    edge = np.array([1.0, 1.0, 0.3])  # 45 degrees across the chain; what it has along the chain does not count
    # C=O along x, and a hydrogen off the plane of the carbonyl and the axis, which gives the chain a hand
    a = Atoms("COH", positions=[[0.0, 0.0, 0.0], [1.2, 0.0, 0.0], [0.3, 0.9, 0.5]], cell=np.diag([6.0, 6.0, 4.7]), pbc=True)
    b = CT.mirrored_chain(a, edge, z)
    co = b.positions[1] - b.positions[0]
    assert co == pytest.approx([0.0, 1.2, 0.0], abs=1e-12)  # x and y change places across the line at 45 degrees
    assert b.positions.mean(axis=0) == pytest.approx(a.positions.mean(axis=0))
    assert b.positions[:, 2] == pytest.approx(a.positions[:, 2])  # nothing moves along the chain
    hand = lambda q: float(np.cross(q.positions[1] - q.positions[0], q.positions[2] - q.positions[0]) @ z)
    assert hand(b) == pytest.approx(-hand(a))
    # twice is the chain itself
    assert CT.mirrored_chain(b, edge, z).positions == pytest.approx(a.positions)
    # a methylene's two hydrogens change names with the reflection: each stays near where it was
    m = Atoms("COHH", positions=[[0.0, 0.0, 0.0], [5.0, 0.0, 0.0], [0.3, -1.0, 0.0], [0.9, -0.2, 0.3]], cell=np.diag([12.0, 12.0, 4.7]), pbc=True)
    assert CT.methylene_pairs(m) == [(2, 3)]
    n = CT.mirrored_chain(m, np.array([1.0, 0.0, 0.0]), z)  # through the plane y = the middle's
    shift = n.positions[0] - m.positions[0] * np.array([1.0, -1.0, 1.0])
    assert n.positions[2] - shift == pytest.approx([0.9, 0.2, 0.3]) and n.positions[3] - shift == pytest.approx([0.3, 1.0, 0.0])


def test_two_states_are_numbered_alike_by_what_moves_their_methylene_hydrogens_least():
    from ase import Atoms

    first = Atoms("CHH", positions=[[0.0, 0.0, 0.0], [0.6, -0.9, 0.0], [0.6, 0.9, 0.0]], cell=np.diag([8.0, 8.0, 8.0]), pbc=True)
    last = first.copy()
    last.positions[[1, 2]] = last.positions[[2, 1]] + [0.05, 0.0, 0.0]  # the same shape with the two exchanged, a little moved
    CT.numbered_alike(first, last)
    assert last.positions[1] == pytest.approx([0.65, -0.9, 0.0]) and last.positions[2] == pytest.approx([0.65, 0.9, 0.0])
    kept = last.positions.copy()
    CT.numbered_alike(first, last)  # already alike: nothing changes
    assert last.positions == pytest.approx(kept)
    # a cell thinner than the move: the carbon has gone 1.9 A across a cell 4 A thick with its hydrogens,
    # and the nearest copy of each hydrogen taken by itself would be the other one's
    thin = Atoms("CHH", positions=[[0.0, 0.0, 0.0], [0.6, -0.9, 0.0], [0.6, 0.9, 0.0]], cell=np.diag([8.0, 4.0, 8.0]), pbc=True)
    moved = thin.copy()
    moved.positions += [0.0, 1.9, 0.0]
    kept = moved.positions.copy()
    CT.numbered_alike(thin, moved)
    assert moved.positions == pytest.approx(kept)


def test_carried_sends_the_atoms_of_a_molecule_the_same_way_in_a_thin_cell():
    from ase import Atoms
    from ase.geometry import find_mic

    # two molecules in a cell 4 A thick; the first goes 2.1 A across it and turns a little, the second stays
    first = Atoms("CHHOH", positions=[[1.0, 0.5, 1.0], [1.6, -0.4, 1.0], [1.6, 1.4, 1.0], [5.0, 1.0, 1.0], [5.9, 1.0, 1.0]],
                  cell=np.diag([8.0, 4.0, 8.0]), pbc=True)
    way = np.array([[0.0, 2.1, 0.0], [0.1, 2.3, 0.0], [-0.1, 1.9, 0.0], [0.0, 0.0, 0.0], [0.0, 0.1, 0.0]])
    last = first.copy()
    last.positions += way
    last.wrap()
    alone = find_mic(last.positions - first.positions, first.cell.array)[0]
    assert alone[1, 1] == pytest.approx(-1.7) and alone[2, 1] == pytest.approx(1.9)  # one hydrogen each way: the bond would be pulled apart
    got = CT.carried(first, last)
    # the molecule goes together, the way its atom that moved least by itself went: the hydrogen 1.7 A from its nearest copy
    assert got[:3] == pytest.approx(way[:3] - [0.0, 4.0, 0.0]) and got[3:] == pytest.approx(way[3:])
    # filled in half way the bonds are still bonds; atom by atom to the nearest copy, one hydrogen is 2 A from its carbon
    mid, apart = first.copy(), first.copy()
    mid.positions += 0.5 * got
    apart.positions += 0.5 * alone
    assert max(mid.get_distance(0, 1), mid.get_distance(0, 2)) < 1.15
    assert apart.get_distance(0, 2) > 1.9


def test_a_twin_has_the_lattice_points_of_the_wall_and_a_row_of_it_is_the_crystal_sheared_along_the_wall():
    from ase import Atoms

    cell = np.array([[4.6, 0.2, 0.3], [1.6, 3.3, -0.1], [-0.3, 0.1, 4.7]])
    turn = CT.half_turn_across(cell)
    assert turn @ cell[1] == pytest.approx(-cell[1]) and turn @ cell[2] == pytest.approx(-cell[2])
    assert turn @ turn == pytest.approx(np.eye(3)) and np.linalg.det(turn) == pytest.approx(1.0)
    step = CT.twin_step(cell)
    # the step is the first edge's image and a whole number of second edges, and the change of step lies in the wall
    assert np.linalg.solve(cell.T, turn @ cell[0] - step) == pytest.approx([0.0, -1.0, 0.0], abs=1e-9)
    normal = np.cross(cell[1], cell[2])
    assert (cell[0] - step) @ normal == pytest.approx(0.0, abs=1e-9)
    assert np.linalg.norm(step) == pytest.approx(np.linalg.norm(cell[0] - cell[1]))  # the crystal's diagonal is the twin's row

    # a stand-in chain: C=O across the chain with an N-H behind it and a hydrogen off their plane
    a = Atoms("COHNH", positions=[[1.0, 1.0, 1.0], [2.2, 1.0, 1.0], [-0.8, 1.0, 1.0], [0.2, 1.0, 1.0], [1.0, 1.5, 3.0]], cell=cell, pbc=True)
    a.info.update(triage_repeats=1, kmesh="6 6 6", triage_chain="x", triage_state="tt")
    along = cell[2] / np.linalg.norm(cell[2])
    for twin in (0, 2, 3):
        row = CT.twin_row(a, 5, 1, twin)
        assert row.cell.array[0] == pytest.approx((5 - twin) * cell[0] + twin * step)
        assert row.cell.array[1:] == pytest.approx(cell[1:])
        assert row.get_volume() == pytest.approx(5 * a.get_volume())  # a shear along the wall keeps the volume
        assert row.info["triage_twin_chains"] == twin and row.info["kmesh"] == "1 6 6"
        v = a.positions[1] - a.positions[0]
        v = v - (v @ along) * along
        for k in range(5):
            part = np.where(row.arrays["chain"] == k)[0]
            co = row.positions[part[1]] - row.positions[part[0]]
            co = co - (co @ along) * along
            # a twin chain's carbonyl is where the half turn has it; the others are as the crystal has them
            assert co == pytest.approx(turn @ v if k >= 5 - twin else v, abs=1e-6)
            # the chains follow one another by the crystal's step, then by the twin's
            site = min(k, 5 - twin) * cell[0] + max(0, k - (5 - twin)) * step
            assert row.positions[part].mean(axis=0) - (a.positions.mean(axis=0) + site) == pytest.approx(np.zeros(3), abs=0.35)
    # the cell told for half a chain more in the twin: the first edge between the two, the atoms carried with it
    mid2, mid3 = CT.twin_row(a, 5, 1, 2, 2.5), CT.twin_row(a, 5, 1, 3, 2.5)
    assert mid2.cell.array == pytest.approx(mid3.cell.array)
    assert mid2.cell.array[0] == pytest.approx(2.5 * cell[0] + 2.5 * step)
    assert mid2.get_scaled_positions() == pytest.approx(CT.twin_row(a, 5, 1, 2).get_scaled_positions(), abs=1e-9)
    with pytest.raises(ValueError):
        CT.twin_row(a, 5, 1, 6)


def test_part_of_a_chain_set_over_moves_those_repeats_whole_and_nothing_else():
    from ase import Atoms

    # two chains of bonded atoms three repeats long (two atoms a repeat); in the second state the second chain is turned
    cell = np.diag([8.0, 4.0, 7.5])
    pos = []
    for k in range(2):
        for r in range(3):
            pos += [[4.0 * k + 0.5, 1.0, 2.5 * r], [4.0 * k + 0.5, 1.9, 2.5 * r + 1.25]]
    first = Atoms("CN" * 6, positions=pos, cell=cell, pbc=True)
    first.set_array("chain", np.repeat([0, 1], 6))
    last = first.copy()
    centre = np.array([4.5, 1.45, 0.0])
    t = np.radians(50.0)
    spin = np.array([[np.cos(t), -np.sin(t), 0.0], [np.sin(t), np.cos(t), 0.0], [0.0, 0.0, 1.0]])
    d = last.positions[6:] - centre
    last.positions[6:] = centre + np.column_stack([(d[:, :2] @ spin[:2, :2].T), d[:, 2]])
    out = CT.part_over(first, last, 1, [1], 2, 6)
    assert out.positions[8:10] == pytest.approx(last.positions[8:10])  # the middle repeat of the second chain, as the second state has it
    keep = [i for i in range(12) if i not in (8, 9)]
    assert out.positions[keep] == pytest.approx(first.positions[keep])
    assert CT.part_over(first, last, 1, [0, 1, 2], 2, 6).positions == pytest.approx(last.positions)
    with pytest.raises(ValueError):
        CT.part_over(first, last, 1, [3], 2, 6)


def test_a_methyl_on_every_methylene_of_one_chain_keeps_the_ureas_and_the_other_chains():
    from ase import Atoms

    # two stand-in chains two repeats long: N-H, C=O and a methylene to a repeat
    cell = np.diag([9.0, 6.0, 10.0])
    one = [("N", [0.0, 0.0, 0.0]), ("H", [0.0, 1.0, 0.0]), ("C", [0.0, -0.6, 1.2]), ("O", [0.0, -1.8, 1.2]),
           ("C", [0.0, 0.3, 1.9]), ("H", [0.9, 0.9, 1.9]), ("H", [-0.9, 0.9, 1.9])]
    symbols, pos = [], []
    for k in range(2):
        for r in range(2):
            for sym, xyz in one:
                symbols.append(sym)
                pos.append(np.array(xyz) + [4.5 * k, 0.0, 5.0 * r])
    row = Atoms(symbols, positions=pos, cell=cell, pbc=True)
    row.set_array("chain", np.repeat([0, 1], 14))
    out = CT.c_methylated(row, 14, 1, 0)
    assert len(out) == len(row) + 6 and out.get_chemical_symbols().count("C") == row.get_chemical_symbols().count("C") + 2
    assert list(out.arrays["chain"][len(row):]) == [1] * 6 and list(out.arrays["chain"][:len(row)]) == list(row.arrays["chain"])
    # the first chain is as it was; of the second, each methylene keeps one hydrogen and gains a carbon 1.53 A off
    assert out.positions[:14] == pytest.approx(row.positions[:14])
    for r in range(2):
        c_at, h0, h1 = 14 + 7 * r + 4, 14 + 7 * r + 5, 14 + 7 * r + 6
        assert out.symbols[h0] == "C" and out.symbols[h1] == "H"
        assert np.linalg.norm(out.positions[h0] - out.positions[c_at]) == pytest.approx(1.53)
        mine = out.positions[len(row) + 3 * r:len(row) + 3 * r + 3]
        assert np.linalg.norm(mine - out.positions[h0], axis=1) == pytest.approx([1.09] * 3)
    # the other face
    other = CT.c_methylated(row, 14, 1, 1)
    assert other.symbols[14 + 6] == "C" and other.symbols[14 + 5] == "H"
    # the N-H and the carbonyl are untouched
    assert CT.carbonyls(out, np.where(out.arrays["chain"] == 1)[0]) == CT.carbonyls(row, slice(14, 28))
    with pytest.raises(ValueError):
        CT.c_methylated(Atoms("CO", positions=[[0, 0, 0], [0, 1.2, 0]], cell=cell, pbc=True), 2, 0)


def test_a_twin_of_a_chain_without_a_carbonyl_is_the_chain_turned_and_not_reflected():
    from ase import Atoms

    # a stand-in for PVDF all-trans in a cell a little off hexagonal: fluorines one side, hydrogens the other,
    # one fluorine written a cell away as a packing search may leave it
    cell = np.array([[4.88, 0.0, 0.0], [-2.56, 4.16, 0.0], [0.0, 0.0, 2.58]])
    pos = np.array([[0.0, 0.0, 0.0], [0.0, 0.85, 1.29], [1.1, -0.75, 0.0], [-1.1, -0.75, 0.0], [0.9, 1.45, 1.29], [-0.9, 1.45, 1.29]])
    turn = np.radians(-20.0)  # the chain's dipole not along an edge
    spin = np.array([[np.cos(turn), -np.sin(turn), 0.0], [np.sin(turn), np.cos(turn), 0.0], [0.0, 0.0, 1.0]])
    a = Atoms("CCFFHH", positions=pos @ spin.T, cell=cell, pbc=True)
    a.positions[3] += cell[0]
    a.info.update(triage_repeats=1, kmesh="4 4 8", triage_chain="pvdf", triage_state="tt")
    along = np.array([0.0, 0.0, 1.0])
    v = CT.way_across(a)
    v = v - (v @ along) * along  # the fluorines and the hydrogens are half a repeat apart along the chain
    w = CT.half_turn_across(cell) @ v
    row = CT.twin_row(a, 4, 1, 2)
    assert row.info["triage_twin_chains"] == 2 and row.get_volume() == pytest.approx(4 * a.get_volume())
    whole = CT.whole_in_place(a)
    for k in range(4):
        part = np.where(row.arrays["chain"] == k)[0]
        got = CT.chain_facing(row, part, along, v)
        want = np.degrees(np.arctan2(np.cross(v, w) @ along, v @ w)) if k >= 2 else 0.0
        assert got == pytest.approx([want], abs=1e-3)  # the turn is handed on to four places
        # turned about its own axis and not reflected: every bond as long as it was
        d = row.positions[part] - row.positions[part[0]]
        assert np.linalg.norm(d, axis=1) == pytest.approx(np.linalg.norm(whole.positions - whole.positions[0], axis=1), abs=1e-9)
    # with a carbonyl the facing is the carbonyl's, as before
    b = Atoms("COH", positions=[[0.0, 0.0, 0.0], [0.0, 1.2, 0.0], [0.5, -0.9, 0.0]], cell=np.diag([5.0, 5.0, 2.6]), pbc=True)
    assert CT.chain_facing(b, slice(None), along, np.array([1.0, 0.0, 0.0])) == pytest.approx([90.0])


def test_a_state_let_go_keeps_its_mesh_and_is_told_against_the_crystal_it_was_made_of():
    from ase import Atoms

    cell = np.array([[4.6, 0.2, 0.3], [1.6, 3.3, -0.1], [-0.3, 0.1, 4.7]])
    parent = Atoms("CO", positions=[[1.0, 1.0, 1.0], [2.2, 1.0, 1.0]], cell=cell, pbc=True)
    # two chains across and three repeats along, as a states file DFT relaxed has it (the mesh read back as numbers)
    state = parent.repeat((2, 1, 3))
    state.info.update(triage_chain="x", triage_state="tt", triage_kind="crystal", triage_repeats=6, kmesh=np.array([2, 5, 1]),
                      triage_domain_faces="twin 1", triage_turn_of="x_w1", triage_turn_ev=0.0, triage_path_unit="for the cell",
                      triage_state_stem="t1", triage_dft_ev=-1.0)
    c = CT.let_go_start(state, parent, "x_tt_x1")
    assert c.info["kmesh"] == "2 5 1" and c.info["triage_let_go_of"] == "x_tt_x1" and c.info["triage_repeats_along"] == 3
    assert np.asarray(c.info["triage_let_go_cell"]).reshape(3, 3) == pytest.approx(state.cell.array)
    assert c.info["triage_state"] == "tt" and c.info["triage_repeats"] == 6 and c.info["triage_domain_faces"] == "twin 1"
    assert not [k for k in c.info if k.startswith(("triage_turn", "triage_path", "triage_dft", "triage_state_stem"))]
    assert c.positions == pytest.approx(state.positions) and c.calc is None
    # what dft writes with it is a relaxation of every edge
    assert c.info["triage_cell_free"] == "all"
    assert "calculation = 'vc-relax'" in CT.pw_input(c, "/p", "/o", "vc-relax", "all") and "cell_dofree = 'all'" in CT.pw_input(c, "/p", "/o", "vc-relax", "all")
    # or of the first edge only, the wall's plane held
    w = CT.let_go_start(state, parent, "x_tt_x1", "epitaxial_bc")
    assert w.info["triage_cell_free"] == "epitaxial_bc"
    assert "cell_dofree = 'epitaxial_bc'" in CT.pw_input(w, "/p", "/o", "vc-relax", w.info["triage_cell_free"])


def test_a_chain_without_a_carbonyl_faces_from_its_hydrogens_to_its_fluorines():
    from ase import Atoms

    z = np.array([0.0, 0.0, 1.0])
    cell = np.diag([5.0, 5.0, 2.6])
    # fluorines on +x, hydrogens on -x; one fluorine written a cell away, as a packing search may leave it
    a = Atoms("CCFFHH", positions=[[0.0, 0.0, 0.0], [0.0, 0.0, 1.2], [1.1, 0.8, 0.0], [1.1 - 5.0, -0.8, 0.0], [-0.6, 0.7, 1.2], [-0.6, -0.7, 1.2]], cell=cell, pbc=True)
    way = CT.way_across(a)
    assert way[:2] == pytest.approx([1.7, 0.0])
    assert CT.faces_at(a, z, np.array([1.0, 0.0, 0.0])) == pytest.approx(0.0)
    w = CT.whole_in_place(a)
    assert w.positions[3] == pytest.approx([1.1, -0.8, 0.0]) and w.positions[[0, 1, 2, 4, 5]] == pytest.approx(a.positions[[0, 1, 2, 4, 5]])
    CT.spun(w, len(w), 0, np.radians(60.0), z)
    assert CT.faces_at(w, z, np.array([1.0, 0.0, 0.0])) == pytest.approx(60.0)
    # with a carbonyl it is the carbonyl
    b = Atoms("COH", positions=[[0.0, 0.0, 0.0], [0.0, 1.2, 0.0], [0.5, -0.9, 0.0]], cell=cell, pbc=True)
    assert CT.way_across(b) == pytest.approx([0.0, 1.2, 0.0]) and CT.faces_at(b, z, np.array([1.0, 0.0, 0.0])) == pytest.approx(90.0)


def test_the_row_frame_is_square_with_the_chain_and_the_row_for_axes():
    cell = np.array([[4.6, 0.0, 0.3], [1.7, 3.3, 0.2], [0.0, 0.0, 4.7]])
    frame = CT.row_frame(cell)
    assert frame @ frame.T == pytest.approx(np.eye(3))
    assert frame[2] == pytest.approx([0.0, 0.0, 1.0])  # the chain
    assert frame[0] == pytest.approx([1.0, 0.0, 0.0])  # the row, without what it has along the chain
    assert np.linalg.det(frame) == pytest.approx(1.0)


def test_principal_values_come_with_their_axes_from_the_row():
    (lo, lo_at), (hi, hi_at) = CT.principal([[0.0, 2.0, 9.0], [2.0, 0.0, 9.0], [9.0, 9.0, 9.0]])  # a shear across the chain; the rest is not looked at
    assert (lo, hi) == pytest.approx((-2.0, 2.0))
    assert (lo_at, hi_at) == pytest.approx((135.0, 45.0))
    (lo, lo_at), (hi, hi_at) = CT.principal([[1.0, 0.0], [0.0, 3.0]])
    assert (lo, hi) == pytest.approx((1.0, 3.0)) and hi_at == pytest.approx(90.0) and lo_at % 180.0 == pytest.approx(0.0)


def test_stretches_take_one_cell_into_another_edge_for_edge_and_leave_a_turn_out():
    cell = np.array([[4.6, 0.0, 0.0], [1.7, 3.3, 0.0], [0.0, 0.0, 4.7]])
    # the same crystal turned about the chain: no stretch
    t = np.radians(40.0)
    turn = np.array([[np.cos(t), -np.sin(t), 0.0], [np.sin(t), np.cos(t), 0.0], [0.0, 0.0, 1.0]])
    assert [v for v, _ in CT.stretches(cell, cell @ turn.T)] == pytest.approx([0.0, 0.0, 0.0], abs=1e-12)
    # 10% longer across the row, 5% shorter along the chain, then turned: the stretches and their axes in the row's frame
    grown = cell * np.array([1.0, 1.1, 0.95])
    got = CT.stretches(cell, grown @ turn.T)
    assert [v for v, _ in got] == pytest.approx([-0.05, 0.0, 0.1])
    assert np.abs(got[0][1]) == pytest.approx([0.0, 0.0, 1.0], abs=1e-9) and np.abs(got[2][1]) == pytest.approx([0.0, 1.0, 0.0], abs=1e-9)
    # the sheets slid along the row: a shear, one way longer and the other shorter by nearly as much
    slid = cell.copy()
    slid[1, 0] -= 0.33
    lo, mid, hi = (v for v, _ in CT.stretches(cell, slid))
    assert lo == pytest.approx(-hi / (1 + hi)) and mid == pytest.approx(0.0, abs=1e-12) and hi == pytest.approx(0.0513, abs=2e-4)


def test_a_hydrogen_bond_is_told_with_the_edge_it_lies_along():
    from ase import Atoms

    cell = np.array([[4.6, 0.0, 0.0], [1.7, 3.9, 0.0], [0.0, 0.0, 4.7]])
    # N-H pointing along the row at the oxygen of the next cell over, 1.9 A from it; its own nitrogen is no bond of this kind
    a = Atoms("ONH", positions=[[0.5, 1.0, 1.0], [2.2, 1.0, 1.0], [3.2, 1.0, 1.0]], cell=cell, pbc=True)
    assert CT.hydrogen_bonds(a) == [(1.9, [1, 0])]
    # the same hydrogen turned to face up the stack
    b = a.copy()
    b.positions[2] = b.positions[0] + cell[1] - 1.8 * cell[1] / np.linalg.norm(cell[1])
    assert CT.hydrogen_bonds(b, within=2.2) == [(1.8, [0, 1])]
    # further out the oxygen of its own cell is in reach as well, the shorter first
    assert CT.hydrogen_bonds(b) == [(1.8, [0, 1]), (2.45, [0, 0])]


def test_resting_places_gather_the_starts_that_end_alike_and_close_the_circle():
    def state(start, rest, energy, rigid):
        return {"start_degrees": start, "rest_facing_degrees": rest, "rest_ev": energy, "rigid_ev": rigid, "hydrogen_bonds": [], "stress_gpa": np.zeros((3, 3)).tolist()}

    states = [state(0.0, -2.0, -10.00, -10.00), state(15.0, -1.9, -10.01, -9.9), state(60.0, 55.0, -9.9, -9.8), state(170.0, 178.0, -10.015, -9.95),
              state(185.0, -178.0, -10.02, -9.9)]
    places = CT.resting_places(states)
    assert [round(p["rest_facing_degrees"]) for p in places] == [-178, -2, 55]
    assert [sorted(p["starts"]) for p in places] == [[170.0, 185.0], [0.0, 15.0], [60.0]]
    assert [p["rest_ev"] for p in places] == [-10.02, -10.01, -9.9]  # each place by its lowest state
    # two that face alike and are not alike in energy are two places: a chain and its other hand
    more = CT.resting_places(states + [state(139.0, 57.0, -9.85, -9.7)])
    assert [round(p["rest_facing_degrees"]) for p in more] == [-178, -2, 55, 57] and [p["starts"] for p in more][2:] == [[60.0], [139.0]]


def test_a_cell_is_whole_only_while_it_is_bonded_as_the_repeat_is():
    """A chain whose bond across the period has been pulled apart is no longer the chain, and is told from it."""
    pytest.importorskip("rdkit")
    chain = CT.CHAINS["pvdf"]
    a = CT.periodic(chain, "tt")
    assert CT.repeat_bonding(chain) == ["C:CCFF", "C:CCHH", "F:C", "F:C", "H:C", "H:C"]
    assert CT.whole(a, chain)
    cell = a.cell.array.copy()
    cell[2, 2] += 0.9  # the atoms stay: the carbon at the top no longer reaches the next period's
    a.set_cell(cell, scale_atoms=False)
    assert not CT.whole(a, chain)
    # two repeats are counted as two
    assert CT.repeat_bonding(chain, 2) == sorted(CT.repeat_bonding(chain) * 2)


def test_a_methyl_takes_the_place_of_the_lone_hydrogen_on_a_carbon():
    """The two-methyl urea chain is made from the one-methyl chain's period, heavy atoms in the repeat's order."""
    pytest.importorskip("rdkit")
    a = CT.periodic(CT.CHAINS["urea1m"], "tt")
    assert CT.CHAINS["urea1mm"]["methyl_on"] == "urea1m"
    b = CT.methyl_for_hydrogen(a)
    assert len(b) == len(a) + 3
    assert b.get_chemical_symbols() == list("NCONCCC") + ["H"] * 8
    assert np.allclose(b.positions[:6], a.positions[:6]) and np.allclose(b.cell.array, a.cell.array)
    assert b.get_distance(4, 6) == pytest.approx(1.53, abs=1e-6)  # on the carbon that carried the hydrogen
    assert sorted(round(float(b.get_distance(6, k)), 2) for k in (12, 13, 14)) == [1.09, 1.09, 1.09]
    assert CT.cell_bonding(b).count("C:CHHH") == 2
    # a period with no such carbon, or with two, is refused
    with pytest.raises(ValueError):
        CT.methyl_for_hydrogen(b)

