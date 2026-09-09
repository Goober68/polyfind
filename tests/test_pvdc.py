"""Phase 1 of the chemistry extension: PVDC as a drop-in for the existing model.

PVDC is the one other target chemistry the current monomer model can express
exactly, since it has two identical single-atom substituents per backbone atom.
These tests assert that its geometry is built correctly and that the whole
funnel runs on it. They deliberately assert nothing about its *energies*: with
the illustrative potential and rigid bond angles the all-trans chain is heavily
strained (see the module docstring note in ``polymers.py``), so its RIS numbers
are not meaningful. That is a finding about the potential, not a reason to
withhold the chemistry.
"""
import numpy as np
import pytest

from polyfind.chain import angle, build_chain, distance
from polyfind.enumerate import enumerate_periodic
from polyfind.forcefield import SimpleFF, fit_ris
from polyfind.helix import helix_parameters
from polyfind.pack import CrystalPacker, periodic_chain
from polyfind.polymers import PVDC, THREE_STATE, get_polymer


def test_pvdc_is_registered():
    assert get_polymer("pvdc") is PVDC
    assert PVDC.bonds_per_repeat == 2 and PVDC.atoms_per_repeat == 6


def test_pvdc_geometry_is_built_as_specified():
    st = build_chain(PVDC, np.full(6, 180.0))
    assert st.n_atoms == 9 + 18 + 2
    assert sum(1 for e in st.elements if e == "Cl") == 8  # four CCl2 groups
    assert abs(float(st.charges.sum())) < 1e-12
    for k in (1, 3, 5, 7):  # the CCl2 backbone atoms
        idx = st.backbone[k]
        s1, s2 = st.subs_of[idx]
        assert distance(st.coords, idx, s1) == pytest.approx(1.77, abs=1e-9)
        assert angle(st.coords, s1, idx, s2) == pytest.approx(110.0, abs=1e-6)
    # no two atoms occupy the same place
    d = np.linalg.norm(st.coords[:, None] - st.coords[None], axis=-1) + np.eye(st.n_atoms) * 9
    assert d.min() > 1.0


def test_pvdc_chain_repeats_are_sane():
    for seq, lo, hi in ((THREE_STATE and [0, 0], 2.4, 2.8), ([0, 1, 0, 2], 4.0, 5.0)):
        h = helix_parameters(PVDC, np.array(seq), THREE_STATE)
        assert h.c is not None and lo < h.c < hi
        assert h.radius_all > h.radius_backbone > 0


def test_pvdc_runs_the_whole_funnel():
    """Fit, enumerate, build a periodic chain and score a cell: the harness works."""
    rep = fit_ris(PVDC, SimpleFF(), step=30.0, n_monomers=4, third_order=False)
    m = rep.model
    assert m.B == 2 and m.second_order.shape == (2, 3, 3)
    mir = np.array(THREE_STATE.mirror)
    assert np.allclose(m.second_order, m.second_order[:, mir][:, :, mir])  # achiral

    cands = enumerate_periodic(PVDC, m, max_period=4, k_per_period=20)
    assert cands and all(c.helix.rise_per_bond > 0.3 for c in cands)

    ch = periodic_chain(PVDC, m.parse("TT"), m.states)
    pk = CrystalPacker(ch, n_chains=2)
    e = pk.energy(np.array([[6.0, 10.0, 90.0, 30.0, 200.0, 1.0, 0]]))
    assert np.isfinite(e).all()


def test_pvdc_flip_is_free_beyond_the_cutoff():
    """The bonded-exclusion fix must hold for a new chemistry too."""
    ch = periodic_chain(PVDC, [0, 0], THREE_STATE)
    pk = CrystalPacker(ch, n_chains=2)
    far = np.array([[60.0, 60.0, 90.0, 17.0, 143.0, 0.9, 0.0]])
    anti = far.copy()
    anti[0, 6] = 1.0
    assert float(pk.energy(anti)[0]) == pytest.approx(float(pk.energy(far)[0]), abs=1e-9)


def test_all_trans_pvdc_is_strained_relative_to_pvdf():
    """Documents the Phase 1 finding rather than hiding it.

    Chlorine is bulky (UFF x_i 3.95 A) and its C-Cl bond long (1.77 A), so with
    bond angles frozen the planar zigzag cannot relieve Cl...Cl contact. The
    all-trans chain is therefore not a sensible RIS reference for PVDC, which is
    consistent with PVDC not adopting the planar zigzag that PVDF's beta phase
    does. Variable bond angles are the principled fix.
    """
    ff = SimpleFF()
    from polyfind.polymers import PVDF

    lj = {p.name: ff.components(build_chain(p, np.full(10, 180.0)))["lj"] for p in (PVDF, PVDC)}
    assert lj["pvdf"] < 50.0
    assert lj["pvdc"] > 200.0  # heavily strained: rotating any bond lowers the energy
