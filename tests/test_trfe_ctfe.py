"""The TrFE (CHF-CF2) and CTFE (CF2-CFCl) units and the 70/30 VDF/TrFE model copolymer.

Added for the electrostriction plan's tier 2 (``docs/ELECTROSTRICTION_RESEARCH.md``).  As
with the CFE/CDFE tests, nothing here asserts a fitted energy; the machinery (registration,
chirality, neutral repeats, periodic chains of both conformations, the transition strain
and the chiral energy relation) is what is checked.
"""
import warnings

import numpy as np
import pytest

from polyfind import electrostriction as ES
from polyfind.chain import build_chain
from polyfind.forcefield import SimpleFF
from polyfind.helix import helix_parameters
from polyfind.pack import periodic_chain
from polyfind.polymers import CTFE, PVDF, THREE_STATE, TRFE, VDF_TRFE_7_3, get_polymer

T, GP, GM = 0, 1, 2


def test_registered_chiral_and_composed():
    assert get_polymer("trfe") is TRFE and get_polymer("ctfe") is CTFE
    assert get_polymer("vdf7-trfe3") is VDF_TRFE_7_3
    assert TRFE.is_chiral and CTFE.is_chiral and VDF_TRFE_7_3.is_chiral
    assert TRFE.bonds_per_repeat == CTFE.bonds_per_repeat == 2
    assert VDF_TRFE_7_3.bonds_per_repeat == 20 and VDF_TRFE_7_3.monomers_per_repeat == 10
    assert VDF_TRFE_7_3.mol_percent() == {"vdf": 70.0, "trfe": 30.0}
    assert VDF_TRFE_7_3.is_copolymer and not TRFE.is_copolymer
    assert VDF_TRFE_7_3.bond_length == PVDF.bond_length


@pytest.mark.parametrize("polymer,seq", [(TRFE, [T, T]), (CTFE, [T, T]), (VDF_TRFE_7_3, [T] * 20)])
def test_all_trans_periodic_chains_build_neutral_with_the_right_atoms(polymer, seq):
    ch = periodic_chain(polymer, seq, THREE_STATE)
    assert ch.n_atoms == polymer.atoms_per_repeat
    assert abs(float(np.sum(ch.charges))) < 1e-12
    els = "".join(ch.elements)
    if polymer is TRFE:
        assert els == "CHFCFF"
    elif polymer is CTFE:
        assert els == "CFFCFCl"
    else:
        assert els.count("Cl") == 0 and els.count("H") == 2 * 7 + 3 and els.count("F") == 2 * 7 + 3 * 3


def test_the_fluorinated_units_share_pvdfs_all_trans_geometry():
    """Same backbone angles, so the same rise per monomer per unit of bond length as
    beta-PVDF (the units carry the textbook 1.54 A against PVDF's measured 1.528)."""
    ref = helix_parameters(PVDF, np.array([T, T]), THREE_STATE)
    scaled = ref.rise_per_period / ref.monomers_per_period / PVDF.bond_length
    for polymer in (TRFE, CTFE):
        h = helix_parameters(polymer, np.array([T, T]), THREE_STATE)
        assert h.rise_per_period / h.monomers_per_period / polymer.bond_length == pytest.approx(scaled, rel=1e-9)


@pytest.mark.parametrize("polymer", [TRFE, CTFE])
def test_the_tgtg_to_tttt_transition_strain_is_positive(polymer):
    ts = ES.transition_strain(polymer, [T, GP, T, GM], [T, T], THREE_STATE)
    assert 0.07 < ts["axial"] < 0.16
    ch = periodic_chain(polymer, [T, GP, T, GM], THREE_STATE)  # the alpha-like chain builds too
    assert ch.n_atoms == 2 * polymer.atoms_per_repeat


@pytest.mark.parametrize("polymer", [TRFE, CTFE])
def test_the_potential_sees_the_stereocentre(polymer):
    """Reflection composed with chain reversal is a symmetry; reflection alone is not."""
    dih = np.array([180.0, 60.0, 180.0, -60.0, 180.0, 180.0, 60.0, 180.0])
    ff = SimpleFF()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        e = ff.energy(build_chain(polymer, dih))
        rev = ff.energy(build_chain(polymer, -dih[::-1]))
        plain = ff.energy(build_chain(polymer, -dih))
    assert np.isfinite(e)
    assert abs(e - rev) < 1e-6 * max(1.0, abs(e))
    assert abs(e - plain) > 1e-3


def test_the_copolymer_cell_carries_a_smaller_dipole_per_monomer_than_pvdf():
    """A TrFE unit replaces one +0.10 hydrogen by a -0.20 fluorine on the CH2 side, so the
    all-trans cell's transverse dipole per monomer falls below PVDF's; the model's charges
    are illustrative, but the direction of the change is not."""
    def dipole_per_monomer(polymer, seq):
        ch = periodic_chain(polymer, seq, THREE_STATE)
        mu = ch.charges @ ch.coords  # e.A, chain axis z
        return np.linalg.norm(mu[:2]) / polymer.monomer_count(len(seq))
    pvdf = dipole_per_monomer(PVDF, [T, T])
    cop = dipole_per_monomer(VDF_TRFE_7_3, [T] * 20)
    assert pvdf > 0.3
    assert cop < pvdf
