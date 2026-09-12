"""Copolymer composition: an explicit multi-monomer repeat, and the topology check.

``docs/CHEMISTRY_EXTENSION.md`` section 3's first bullet scopes this: make
``Polymer.backbone`` hold a full explicit sequence rather than one cycled monomer, reuse
the comonomers' fitted energies for the background bonds, and fit only the junction terms.
The claim that it "mostly falls out of the existing bond-type indexing" is what most of
this file checks -- that a 24-bond repeat needs no new machinery anywhere, and that every
existing polymer is the degenerate one-monomer case of it, unchanged.

The rest checks :mod:`polyfind.topology`, which exists because the consumer of these
structures (``docs/NOTE_VDCN_CONVERGENCE.md``) rejected a seed whose chains changed bond
topology during relaxation.  A start that passes our own energy but trips a covalent-cutoff
detector is a failed deliverable, so the detector is tested in both directions: a sane
packing must pass over a wide band of cutoff scales, and a deliberately overlapped one must
fail in exactly the way theirs did.
"""
import warnings

import numpy as np
import pytest

from polyfind.chain import build_chain
from polyfind.forcefield import SimpleFF, fit_ris
from polyfind.pack import CrystalPacker, periodic_chain
from polyfind.polymers import (
    AN,
    CDFE,
    CFE,
    FANOME,
    PE,
    PVDC,
    PVDF,
    THREE_STATE,
    VDCN,
    VDCN_UNIT,
    VDF_UNIT,
    VDF_VDCN_11_1,
    Monomer,
    Polymer,
    copolymer,
    get_polymer,
    monomer_of,
)
from polyfind.ris import transfer_ris, transfer_summary
from polyfind.topology import build_cell, check_topology, intended_edges, repeat_bond_graph

HOMO = [PE, PVDF, PVDC, CFE, CDFE, AN, VDCN, FANOME]
CP = VDF_VDCN_11_1


# ----------------------------------------------------- 1. homopolymers are unchanged
@pytest.mark.parametrize("polymer", HOMO)
def test_homopolymer_is_the_degenerate_one_monomer_case(polymer):
    assert polymer.monomers_per_repeat == 1
    assert not polymer.is_copolymer
    assert polymer.sequence[0].backbone == polymer.backbone
    assert polymer.sequence[0].name == polymer.name
    assert polymer.composition() == {polymer.name: 1}
    # every count derived from the sequence is what the inline expression used to give
    for n in (polymer.bonds_per_repeat, 6 * polymer.bonds_per_repeat):
        assert polymer.monomer_count(n) == n // polymer.bonds_per_repeat
    assert polymer.monomer_of_atom == tuple([0] * polymer.bonds_per_repeat)


def test_monomer_count_refuses_a_partial_repeat():
    with pytest.raises(ValueError, match="whole number"):
        PVDF.monomer_count(3)


def test_periodic_chain_monomer_count_is_unchanged_for_a_homopolymer():
    ch = periodic_chain(PVDF, [0, 0], THREE_STATE)
    assert ch.n_monomers == 1
    ch4 = periodic_chain(PVDF, [0] * 4, THREE_STATE)
    assert ch4.n_monomers == 2


# --------------------------------------------------------- 2. the explicit 11:1 repeat
def test_composition():
    assert CP.bonds_per_repeat == 24
    assert CP.monomers_per_repeat == 12
    assert CP.is_copolymer
    assert CP.composition() == {"vdf": 11, "vdcn": 1}
    assert CP.mol_percent()["vdcn"] == pytest.approx(100 / 12)
    assert CP.atoms_per_repeat == 74
    assert not CP.is_chiral  # two identical nitriles are not a stereocentre
    assert get_polymer("vdf11-vdcn1") is CP
    assert CP.monomer_count(24) == 12 and CP.monomer_count(96) == 48


def test_only_one_backbone_atom_differs_from_pvdf():
    """The chemical content of the substitution, stated as a test.

    VDCN's own CH2 entry is field-for-field PVDF's, so in this model the 11:1 copolymer
    differs from the homopolymer at exactly ONE of the 24 backbone atoms -- the cyano
    carbon -- and the junction it makes is two backbone bonds wide.  Every "which bonds got
    which parameters" statement elsewhere rests on this.
    """
    differ = [k for k in range(24) if CP.backbone[k] != PVDF.backbone[k % 2]]
    assert differ == [13]
    assert CP.monomer_of_atom[13] == 6  # the VDCN unit sits at monomer index 6
    assert CP.sequence[6].name == "vdcn"
    assert CP.backbone[13].substituents == ("C", "C")  # the two nitrile carbons


def test_repeat_is_neutral_and_has_the_requested_formula():
    st = build_chain(CP, np.full(24, 180.0), cap=False)
    ch = periodic_chain(CP, [0] * 24, THREE_STATE)
    assert ch.n_atoms == 74
    assert abs(float(ch.charges.sum())) < 1e-12
    counts = {e: ch.elements.count(e) for e in set(ch.elements)}
    assert counts == {"C": 26, "H": 24, "F": 22, "N": 2}
    # eight such chains are exactly the consumer's C208H192F176N16 at 592 atoms
    assert {e: 8 * n for e, n in counts.items()} == {"C": 208, "H": 192, "F": 176, "N": 16}
    assert 8 * ch.n_atoms == 592
    assert st.n_atoms > 0


def test_chain_geometry_matches_the_homopolymer_backbone():
    """The copolymer's backbone is geometrically PVDF's: same bond, same angles.

    That is why the homopolymer packing is a legitimate seed, and why the repeat is exactly
    twelve times the homopolymer's c.
    """
    cp = periodic_chain(CP, [0] * 24, THREE_STATE)
    hp = periodic_chain(PVDF, [0, 0], THREE_STATE)
    assert cp.c == pytest.approx(12 * hp.c, abs=1e-9)
    assert cp.helix.label == "planar zigzag"


def test_polymer_refuses_a_sequence_that_does_not_match_its_backbone():
    with pytest.raises(ValueError, match="monomer sequence"):
        Polymer(name="bad", backbone=PVDF.backbone, bond_length=1.528,
                sequence=(Monomer("a", PVDF.backbone), Monomer("b", PVDF.backbone)))


def test_copolymer_refuses_a_silent_bond_length_choice():
    with pytest.raises(ValueError, match="backbone bond length"):
        copolymer("x", [PVDF, VDCN])  # 1.528 vs 1.54
    ok = copolymer("x", [PVDF, VDCN], bond_length=1.528)
    assert ok.bond_length == 1.528 and ok.monomers_per_repeat == 2


def test_monomer_of_round_trips():
    m = monomer_of(PVDF)
    assert m.name == "pvdf" and m.source == "pvdf" and m.backbone == PVDF.backbone
    assert VDF_UNIT.source == "pvdf" and VDF_UNIT.name == "vdf"
    assert VDCN_UNIT.source == "vdcn"
    with pytest.raises(ValueError, match="not a single-monomer repeat"):
        monomer_of(CP)
    with pytest.raises(ValueError, match="no backbone atoms"):
        Monomer("empty", ())


# ------------------------------------------------------- 3. the parameter provenance
@pytest.fixture(scope="module")
def sources():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        vdf = fit_ris(PVDF, SimpleFF(), step=20.0, n_monomers=4, third_order=True, adapt_angles=False)
        # angles="rigid": the transfer tests compare the copolymer's arrays with their source
        # arrays, whichever scan produced them, and VDCN's default (relaxed) scan is slow
        cn = fit_ris(VDCN, SimpleFF(), step=20.0, n_monomers=4, third_order=True, adapt_angles=False, angles="rigid")
    return [(PVDF, vdf.model), (VDCN, cn.model)]


def test_transfer_assigns_vdf_to_the_background_and_flags_the_junction(sources):
    model, report = transfer_ris(CP, sources)
    assert model.B == 24 and model.order == 3
    first = {t.bond: t for t in report if t.order == 1}
    # 20 of 24 dihedrals sit entirely inside a VDF stretch and carry PVDF's fitted values
    exact = sorted(b for b, t in first.items() if t.exact)
    assert exact == [b for b in range(24) if b not in (10, 11, 12, 13)]
    for b in exact:
        assert first[b].source == "pvdf" and first[b].source_bond == b % 2
    # the two bonds whose rotating bond touches the cyano carbon take VDCN's own values;
    # the two either side of them take VDF's.  None of the four is a fitted value for its
    # own environment -- each matches its source in three backbone atoms of four.
    assert (first[11].source, first[11].source_bond) == ("vdcn", 1)
    assert (first[12].source, first[12].source_bond) == ("vdcn", 0)
    assert (first[10].source, first[10].source_bond) == ("pvdf", 0)
    assert (first[13].source, first[13].source_bond) == ("pvdf", 1)
    for b in (10, 11, 12, 13):
        assert not first[b].exact and first[b].window_matches == 3 and first[b].core_matches == 2


def test_transferred_background_arrays_are_literally_the_source_arrays(sources):
    (_, vdf_model), (_, cn_model) = sources
    model, report = transfer_ris(CP, sources)
    for b in range(24):
        if b not in (10, 11, 12, 13):
            assert np.array_equal(model.first_order[b], vdf_model.first_order[b % 2])
    assert np.array_equal(model.first_order[11], cn_model.first_order[1])
    assert np.array_equal(model.first_order[12], cn_model.first_order[0])
    assert "TRANSFERRED" in transfer_summary(report)


def test_transferred_model_drives_the_existing_solvers(sources):
    """The point of the bullet: B = 24 needs no new mathematics."""
    model, _ = transfer_ris(CP, sources)
    e = model.energy([0] * 24, periodic=True)
    assert np.isfinite(e)
    best = model.cyclic_k_best(24, 2)
    assert len(best) >= 1 and len(best[0][1]) == 24
    e_min, seq = model.minimum(48)
    assert len(seq) == 48 and np.isfinite(e_min)


def test_transfer_refuses_mismatched_state_sets(sources):
    from polyfind.polymers import RISStates
    from polyfind.ris import RISModel

    two = RISStates(names=("T", "G"), angles=(180.0, 60.0), mirror=(0, 1))
    bad = RISModel(two, 2, np.zeros((2, 2)), np.zeros((2, 2, 2)), name="bad")
    with pytest.raises(ValueError, match="different state set"):
        transfer_ris(CP, [(PVDF, bad)])


def test_transfer_warns_when_the_source_angles_were_adapted():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        vdf = fit_ris(PVDF, SimpleFF(), step=20.0, n_monomers=4, adapt_angles=True)
    if vdf.model.states.angles != THREE_STATE.angles:
        with pytest.warns(UserWarning, match="adapted state angles"):
            transfer_ris(CP, [(PVDF, vdf.model)])


# --------------------------------------------------------------- 4. topology checking
def test_repeat_bond_graph_matches_the_built_chain():
    starts, bonds, n_rep = repeat_bond_graph(CP)
    ch = periodic_chain(CP, [0] * 24, THREE_STATE)
    assert n_rep == ch.n_atoms
    assert np.array_equal(starts, np.asarray(ch.backbone))
    # every intended bond is a real bond length, including the one closing onto +c
    for i, j, dk in bonds:
        d = float(np.linalg.norm(ch.coords[i] - (ch.coords[j] + np.array([0.0, 0.0, dk * ch.c]))))
        assert 0.9 < d < 1.9, (i, j, dk, d)
    assert sum(1 for _, _, dk in bonds if dk) == 1  # exactly one bond crosses the boundary
    # 24 backbone bonds + one bond per pendant + the nitrile C-N bonds
    assert len(bonds) == 24 + 48 + 2


@pytest.fixture(scope="module")
def chain():
    return periodic_chain(CP, [0] * 24, THREE_STATE)


def test_tiling_reaches_the_requested_size(chain):
    params = np.array([5.0, 9.5, 90.0, 0.0, 0.0, 0.0, 1.0])
    prim = build_cell(chain, params)
    assert prim.n_atoms == 148 and prim.n_chains == 2
    big = build_cell(chain, params, nx=2, ny=2)
    assert big.n_atoms == 592 and big.n_chains == 8
    assert big.formula() == "C208H192F176N16"
    assert big.lattice[0, 0] == pytest.approx(2 * prim.lattice[0, 0])
    assert big.lattice[2, 2] == pytest.approx(prim.lattice[2, 2])
    assert "Lattice=" in big.to_extxyz() and big.to_extxyz().splitlines()[0] == "592"


def test_reversed_chains_keep_their_closure_bond(chain):
    """A flipped chain's repeat closes towards -c, and the intended graph has to know.

    Without this the closure bond of every antiparallel chain reads as a missing bond and a
    1.53 A "non-bonded contact" -- a check that would fail a perfectly good structure.
    """
    for flip in (0.0, 1.0):
        params = np.array([5.0, 9.5, 90.0, 0.0, 0.0, 0.0, flip])
        cell = build_cell(chain, params)
        want = intended_edges(cell, CP)
        assert len(want) == 2 * (24 + 48 + 2)
        lat = cell.lattice
        for i, j, img in want:
            shift = img[0] * lat[0] + img[1] * lat[1] + img[2] * lat[2]
            d = float(np.linalg.norm(cell.coords[i] - (cell.coords[j] + shift)))
            assert 0.9 < d < 1.9, (i, j, img, d)


def test_a_sane_packing_passes_over_a_wide_band_of_cutoffs(chain):
    params = np.array([5.2, 9.8, 90.0, 0.0, 0.0, 0.0, 1.0])
    rep = check_topology(build_cell(chain, params), CP)
    assert rep.ok
    assert rep.components[1.2] == [74, 74]
    lo, hi = rep.safe_scale_window
    assert lo < 1.05 and hi > 1.4, (lo, hi)
    assert rep.min_interchain_distance > 1.5
    assert set(rep.min_interchain_by_element) >= {("C", "C"), ("F", "H"), ("C", "N")}


def test_an_overlapped_packing_fails_the_way_the_consumer_s_seed_did(chain):
    """The negative control: chains pushed together must be caught, not passed.

    The consumer's rejected seed showed eight 74-atom chains becoming components of 72, 224
    and 296 atoms with new interchain bonds.  A check that cannot reproduce that on a
    deliberately bad cell is not checking anything.
    """
    params = np.array([3.2, 6.0, 90.0, 0.0, 0.0, 0.0, 1.0])
    rep = check_topology(build_cell(chain, params), CP)
    assert not rep.ok
    assert rep.extra_interchain[1.2] > 0
    assert rep.components[1.2] != [74, 74]
    assert rep.min_interchain_distance < 1.4


def test_topology_check_refuses_the_wrong_polymer(chain):
    cell = build_cell(chain, np.array([5.2, 9.8, 90.0, 0.0, 0.0, 0.0, 1.0]))
    with pytest.raises(ValueError, match="atoms per repeat"):
        check_topology(cell, PVDF)


# ------------------------------------------------------ 5. the packer counts monomers
def test_energy_per_monomer_is_per_monomer(chain):
    assert chain.n_monomers == 12
    pk = CrystalPacker(chain, n_chains=2)
    params = np.array([5.2, 9.8, 90.0, 0.0, 0.0, 0.0, 1.0])
    e = float(pk.energy(params[None])[0])
    assert pk.energy_per_monomer(params[None])[0] == pytest.approx(e / 24)
    res = pk.result(params)
    assert res.energy_per_monomer == pytest.approx(e / 24)
    assert res.dipole is not None and res.polarization is not None
