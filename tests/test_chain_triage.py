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
