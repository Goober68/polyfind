"""Cells of many independent chains: does the stagger survive, and is it the same potential?

The load-bearing claim of ``examples/copolymer_staggered.py`` is that a staggered eight-chain
cell and the aligned cell it came from are scored by *one* function, so that a difference
between them is a difference between the structures.  That is asserted here against
:class:`polyfind.pack.CrystalPacker` itself rather than against a stored number: the energy
per monomer of a pair potential is invariant under tiling, so the multi-chain sum has no
freedom to disagree.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pytest

from polyfind.ewald import EwaldSpec
from polyfind.forcefield import SimpleFF
from polyfind.pack import CrystalPacker, periodic_chain
from polyfind.polymers import PVDF, THREE_STATE, VDF_VDCN_11_1 as CP
from polyfind.supercell import (
    STAGGER_PATTERNS, SupercellEnergy, cell_density, cell_dipole, cell_polarization,
    cell_to_cif, column_layout, stagger_pattern, staggered_cell,
)
from polyfind.topology import Cell, build_cell, check_topology

PARAMS = np.array([5.2, 9.8, 90.0, 35.0, 215.0, 1.7, 0.0])


@pytest.fixture(scope="module")
def chain():
    return periodic_chain(CP, [0] * CP.bonds_per_repeat, THREE_STATE)


@pytest.fixture(scope="module")
def packer(chain):
    return CrystalPacker(chain, n_chains=2)


@pytest.fixture(scope="module")
def se(packer):
    return SupercellEnergy(packer)


# ------------------------------------------------- 1. it is the packer's own potential
@pytest.mark.parametrize("nx,ny", [(1, 1), (2, 1), (1, 2), (2, 2)])
def test_reproduces_the_packer_under_any_tiling(se, nx, ny):
    assert se.check_against_packer(PARAMS, nx, ny) < 1e-8


@pytest.mark.parametrize("flip", [0.0, 1.0])
def test_reproduces_the_packer_for_an_antiparallel_pair(se, flip):
    """A flipped chain's repeat closes towards -c, so its bonded-scale stack is read at -k."""
    p = PARAMS.copy()
    p[6] = flip
    assert se.check_against_packer(p, 2, 1) < 1e-8


def test_reproduces_a_one_chain_packer(chain):
    pk = CrystalPacker(chain, n_chains=1)
    assert SupercellEnergy(pk).check_against_packer(PARAMS, 2, 2) < 1e-8


def test_reproduces_the_packer_for_a_homopolymer():
    ch = periodic_chain(PVDF, [0, 0], THREE_STATE)
    pk = CrystalPacker(ch, n_chains=2)
    assert SupercellEnergy(pk).check_against_packer(
        np.array([4.6, 8.6, 90.0, 90.0, 90.0, 1.28, 0.0]), 2, 2) < 1e-8


def test_reproduces_an_ewald_packer(chain):
    pk = CrystalPacker(chain, n_chains=2, coulomb="ewald", ewald=EwaldSpec())
    assert SupercellEnergy(pk).check_against_packer(PARAMS, 2, 1) < 1e-8


def test_the_image_bound_is_tight_and_sufficient(se, chain):
    """Pruned images must be empty ones: adding shells back may not move the energy."""
    cell = build_cell(chain, PARAMS, nx=2, ny=2)
    assert se.energy_is_converged(cell, extra=1) == 0.0
    assert se.energy_is_converged(cell, extra=2) == 0.0


def test_the_bonded_exclusions_are_actually_applied(se, chain):
    """Without them the intra-chain 1-2 and 1-3 pairs dominate, so this is a huge number."""
    cell = build_cell(chain, PARAMS)
    assert se.column_correction(cell) < -1e4


@pytest.mark.parametrize("bad,match", [
    ({"valence": SimpleFF()}, "valence"),
    ({"field": (0.01, 0.0, 0.0)}, "field"),
])
def test_refuses_a_potential_it_does_not_carry(chain, bad, match):
    with pytest.raises(ValueError, match=match):
        SupercellEnergy(CrystalPacker(chain, n_chains=2, **bad))


# --------------------------------------------------------------- 2. the stagger itself
def test_zero_stagger_is_the_tiling(chain, packer):
    a = build_cell(chain, PARAMS, nx=2, ny=2, packer=packer)
    b = staggered_cell(chain, PARAMS, np.zeros(8, int), nx=2, ny=2, packer=packer)
    assert np.array_equal(a.coords, b.coords)
    assert np.array_equal(a.lattice, b.lattice)
    assert a.formula() == b.formula()


def test_stagger_refuses_the_wrong_number_of_chains(chain, packer):
    with pytest.raises(ValueError, match="4 entries but the cell has 8"):
        staggered_cell(chain, PARAMS, np.zeros(4, int), nx=2, ny=2, packer=packer)


def test_a_whole_monomer_slide_moves_only_the_comonomer_site(chain):
    """The premise of an integer stagger: the skeleton is monomer-periodic, measured.

    68 of the 74 atoms land exactly on an identical atom; the six that do not are the two
    fluorines of the slot the comonomer occupies and the comonomer's own four atoms.
    """
    X, c = chain.coords, chain.c
    els = list(chain.elements)
    shifted = X + np.array([0.0, 0.0, c / chain.n_monomers])
    moved = []
    for i in range(chain.n_atoms):
        d = min((float(np.linalg.norm(X - (shifted[i] + np.array([0.0, 0.0, k * c])), axis=1).min()),
                 int(np.argmin(np.linalg.norm(X - (shifted[i] + np.array([0.0, 0.0, k * c])), axis=1))))
                for k in (-1, 0, 1))
        if d[0] > 1e-8 or els[d[1]] != els[i]:
            moved.append(els[i])
    assert sorted(moved) == ["C", "C", "F", "F", "N", "N"]


def test_the_stagger_changes_the_energy(se, chain, packer):
    """If the evaluator could not see a stagger the whole exercise would be vacuous."""
    flat = se.energy(build_cell(chain, PARAMS, nx=2, ny=2, packer=packer))
    lad = se.energy(staggered_cell(chain, PARAMS, stagger_pattern("ladder"),
                                   nx=2, ny=2, packer=packer))
    assert abs(lad - flat) > 1.0


def test_the_stagger_does_not_change_the_dipole(se, chain, packer):
    """A rigid slide of a neutral chain cannot move the cell dipole, so polarity is preserved.

    This is what lets the stagger and the polarity be chosen independently: a staggered
    antipolar cell is exactly as antipolar as the aligned one it came from.
    """
    flat = build_cell(chain, PARAMS, nx=2, ny=2, packer=packer)
    for name in STAGGER_PATTERNS:
        cell = staggered_cell(chain, PARAMS, stagger_pattern(name), nx=2, ny=2, packer=packer)
        assert cell_dipole(cell, chain.charges) == pytest.approx(
            cell_dipole(flat, chain.charges), abs=1e-9)


def test_an_antipolar_setting_angle_gives_a_zero_dipole_staggered(chain, packer):
    """phi2 = phi1 + 180 on a chain whose moment is transverse: four up, four down, exactly."""
    p = np.array([10.46, 4.80, 90.0, 90.0, 270.0, 3.81, 0.0])
    for name in STAGGER_PATTERNS:
        cell = staggered_cell(chain, p, stagger_pattern(name), nx=2, ny=2, packer=packer)
        assert np.max(np.abs(cell_polarization(cell, chain.charges))) < 1e-12


def test_the_stagger_keeps_the_bond_graph(chain, packer):
    """The point of the deliverable: a staggered cell must still read as eight 74-atom chains."""
    p = np.array([10.46, 4.80, 90.0, 90.0, 90.0, 2.53, 0.0])
    cell = staggered_cell(chain, p, stagger_pattern("spread", nx=2, ny=1), nx=2, ny=1, packer=packer)
    rep = check_topology(cell, CP, scales=(1.05, 1.15, 1.3))
    assert rep.ok and rep.components[1.15] == [74] * 4
    assert rep.min_interchain_distance > 2.0


def test_density_and_layout(chain, packer):
    cell = build_cell(chain, PARAMS, nx=2, ny=2, packer=packer)
    assert cell_density(cell, chain.mass) == pytest.approx(
        packer.density(PARAMS[0], PARAMS[1], PARAMS[2]))
    assert column_layout() == [(0, 0, 0), (0, 0, 1), (0, 1, 0), (0, 1, 1),
                              (1, 0, 0), (1, 0, 1), (1, 1, 0), (1, 1, 1)]


@pytest.mark.parametrize("name,want", [
    ("aligned", [0, 0, 0, 0, 0, 0, 0, 0]),
    ("ladder", [0, 3, 6, 9, 0, 3, 6, 9]),
    ("spread", [0, 3, 6, 9, 1, 4, 7, 10]),
    ("antiphase", [0, 6, 0, 6, 0, 6, 0, 6]),
    ("sheet", [0, 0, 0, 0, 6, 6, 6, 6]),
])
def test_stagger_patterns_are_what_they_say(name, want):
    assert stagger_pattern(name).tolist() == want


def test_antiphase_is_the_packers_own_dz(chain, packer, se):
    """``antiphase`` is not an independent structure, and the evaluator has to say so.

    ``dz`` already slides chain 2 against chain 1, so offsetting every ``k = 1`` chain by half
    the repeat is the same cell at ``dz + c/2``.  If it came out as anything else, either the
    stagger or ``_place`` would be wrong.
    """
    p = np.array([10.46, 4.80, 90.0, 90.0, 90.0, 2.53, 0.0])
    q = p.copy()
    q[5] = (p[5] + chain.c / 2) % chain.c
    a = se.energy(staggered_cell(chain, p, stagger_pattern("antiphase"), nx=2, ny=2, packer=packer))
    b = se.energy(build_cell(chain, q, nx=2, ny=2, packer=packer))
    assert a == pytest.approx(b, rel=1e-9)


def test_unknown_stagger_pattern_is_refused():
    with pytest.raises(ValueError, match="unknown stagger pattern"):
        stagger_pattern("diagonal")


def test_cell_to_cif_carries_the_supercell_lattice(chain, packer):
    cell = build_cell(chain, PARAMS, nx=2, ny=2, packer=packer)
    text = cell_to_cif(cell, title="t")
    assert "data_t" in text and "'P 1'" in text
    assert f"_cell_length_a {2 * PARAMS[0]:.4f}" in text
    assert f"_cell_length_c {chain.c:.4f}" in text
    assert len(text.strip().splitlines()) == 14 + cell.n_atoms


# ------------------------------------------- 3. the shipped files say what they contain
DELIVERABLES = Path(__file__).resolve().parents[1] / "deliverables"
SUMMARY = DELIVERABLES / "staggered_summary.json"


@pytest.fixture(scope="module")
def shipped():
    if not SUMMARY.exists():
        pytest.skip("run examples/copolymer_staggered.py first")
    with open(SUMMARY) as f:
        return json.load(f)


def read_extxyz(path: Path, n_per_chain: int) -> tuple[Cell, str]:
    """Parse an extended-XYZ cell back out of its own text, with no help from the generator."""
    lines = path.read_text().splitlines()
    n = int(lines[0])
    head = lines[1]
    lat = np.array([float(v) for v in re.search(r'Lattice="([^"]+)"', head).group(1).split()])
    elements, coords = [], []
    for line in lines[2 : 2 + n]:
        parts = line.split()
        elements.append(parts[0])
        coords.append([float(v) for v in parts[1:4]])
    n_chains = n // n_per_chain
    return Cell(elements=elements, coords=np.array(coords), lattice=lat.reshape(3, 3),
                chain_of=np.repeat(np.arange(n_chains), n_per_chain),
                local_of=np.tile(np.arange(n_per_chain), n_chains),
                n_per_chain=n_per_chain), head


@pytest.mark.parametrize("polarity", ["polar", "antipolar"])
@pytest.mark.parametrize("pattern", ["ladder", "spread"])
def test_the_shipped_staggered_file_re_derives_its_own_numbers(shipped, chain, se, polarity, pattern):
    """The deliverable's guarantee: the file, read as text, reproduces what we said about it.

    Re-reading rather than re-building is the point.  A consumer gets the text, not our
    objects, so the energy, the density, the topology verdict *and* the stagger pattern are all
    re-derived from what a parser gets back -- and the header's own key-value pairs are checked
    against the file's geometry rather than against the generator's memory of it.
    """
    ref = shipped["cells"][f"{polarity}_{pattern}"]
    if "files" not in ref:
        pytest.skip(f"{polarity}_{pattern} was not exported")
    cell, head = read_extxyz(DELIVERABLES / ref["files"]["xyz"], chain.n_atoms)
    assert cell.n_atoms == 592 and cell.n_chains == 8
    assert cell.formula() == "C208H192F176N16"

    # Eight decimals of coordinate is enough for the energy only because the bonded-exclusion
    # correction is taken from the cell's own atoms; see SupercellEnergy.column_correction.
    assert se.energy_per_monomer(cell) == pytest.approx(
        ref["energy_per_monomer_kcal_mol"], abs=1e-6)
    assert cell_density(cell, chain.mass) == pytest.approx(ref["density_g_cm3"], abs=1e-6)

    rep = check_topology(cell, CP, scales=tuple(shipped["meta"]["scales"]))
    assert rep.ok
    assert rep.components[1.2] == [74] * 8
    assert all(v == 0 for v in rep.missing.values())
    assert all(v == 0 for v in rep.extra.values())
    assert all(v == 0 for v in rep.extra_interchain.values())
    assert rep.safe_scale_window == pytest.approx(ref["topology"]["safe_scale_window"], abs=1e-6)
    assert rep.min_interchain_distance == pytest.approx(
        ref["topology"]["min_interchain_distance_A"], abs=1e-4)

    # the stagger, recovered from the coordinates alone.  Chains with k = 1 also carry the
    # packer's dz, so each sublattice is measured against its own first chain.
    per = chain.c / chain.n_monomers
    z = [cell.coords[t * chain.n_atoms, 2] for t in range(8)]
    declared = ref["stagger_monomers"]
    for k in (0, 1):
        grp = [t for t in range(8) if t % 2 == k]
        got = [int(round((z[t] - z[grp[0]]) / per)) % chain.n_monomers for t in grp]
        want = [(declared[t] - declared[grp[0]]) % chain.n_monomers for t in grp]
        assert got == want, (k, got, want)

    kv = dict(re.findall(r"(\w+)=(-?[\d.]+)", head))
    assert float(kv["energy_per_monomer_kcal_mol"]) == pytest.approx(
        se.energy_per_monomer(cell), abs=1e-5)
    assert float(kv["density_g_cm3"]) == pytest.approx(cell_density(cell, chain.mass), abs=1e-5)
    assert float(kv["min_interchain_distance_A"]) == pytest.approx(
        rep.min_interchain_distance, abs=1e-3)
    assert re.search(r'stagger_monomers="([\d ]+)"', head).group(1).split() == [
        str(v) for v in declared]


def test_the_shipped_antipolar_staggered_files_are_antipolar(shipped, chain):
    """Zero polarization, measured off the file rather than off the construction.

    Eight decimals of coordinate is not enough to cancel four chains up against four down to
    floating point, so the tolerance here is the file's precision, not the construction's --
    which is 3e-15, and is asserted separately on the relaxed cells.
    """
    for pattern in ("ladder", "spread"):
        ref = shipped["cells"][f"antipolar_{pattern}"]
        if "files" not in ref:
            continue
        cell, _ = read_extxyz(DELIVERABLES / ref["files"]["xyz"], chain.n_atoms)
        assert ref["polarization_magnitude_C_m2"] < 1e-12
        assert np.max(np.abs(cell_polarization(cell, chain.charges))) < 1e-8


def test_the_shipped_staggered_cells_did_not_beat_the_aligned_ones(shipped):
    """The verdict is a fact about these files and is asserted, not just narrated.

    If a future change made a staggered cell come out lower or denser, the README's verdict
    section would be wrong and this is what would say so.
    """
    for p in ("polar", "antipolar"):
        v = shipped["verdict"][p]
        assert v["density_recovered"] is False
        assert v["energy_lowered"] is False
        assert v["d_energy_per_monomer"] > 0.0
        assert "antiphase" in v["patterns_degenerate_with_aligned"]
        assert set(v["genuine_stagger_patterns"]) >= {"ladder", "spread"}
        base = shipped["cells"][f"{p}_aligned"]["energy_per_monomer_kcal_mol"]
        assert shipped["cells"][f"{p}_antiphase"]["energy_per_monomer_kcal_mol"] == pytest.approx(
            base, abs=1e-6)


def test_dipole_refuses_a_charged_cell(chain, packer):
    cell = build_cell(chain, PARAMS, packer=packer)
    q = np.asarray(chain.charges, dtype=float).copy()
    q[0] += 0.5
    with pytest.raises(ValueError, match="not a dipole moment"):
        cell_dipole(cell, q)
