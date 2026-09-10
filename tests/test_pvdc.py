"""Phase 1 of the chemistry extension: PVDC as a drop-in for the existing model.

PVDC is the one other target chemistry the current monomer model can express
exactly, since it has two identical single-atom substituents per backbone atom.
These tests assert that its geometry is built correctly and that the whole
funnel runs on it. They deliberately assert nothing about its *energies*: with
the illustrative potential and rigid bond angles the all-trans chain is heavily
strained (see the module docstring note in ``polymers.py``), so its RIS numbers
are not meaningful. That is a finding about the potential, not a reason to
withhold the chemistry.

One thing changed when PVDC's backbone angles were corrected to their measured
123 deg (CH2) and 114 deg (CCl2): with unequal angles no PVDC chain closes into
a crystallographic repeat at *ideal* torsions, because each two-bond repeat
carries a net 9 deg of curl. The chain that does close is the one with the small
torsion deflections the real polymer has, and it lands on the published
structure -- see ``test_measured_angles_recover_the_published_pvdc_glide_chain``.
So the packing tests here run on that chain rather than on all-trans.
"""
import numpy as np
import pytest

from polyfind.chain import angle, build_chain, distance
from polyfind.enumerate import enumerate_periodic
from polyfind.forcefield import SimpleFF, fit_ris
from polyfind.helix import helix_parameters
from polyfind.pack import CrystalPacker, periodic_chain, periodic_chain_from_torsions
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


def test_unequal_backbone_angles_make_the_planar_pvdc_chain_curl():
    """A consequence of the measured 123/114 deg angles, and it is not a small one.

    An all-trans chain turns by (180 - a) at each backbone atom, alternately left and
    right.  With equal angles those turns cancel and the chain is straight; with PVDC's
    measured C-CH2-C = 123 and C-CCl2-C = 114 they do not, and the chain acquires a net
    123 - 114 = 9 deg of rotation per two-bond repeat.  So at ideal torsions the planar
    PVDC chain is a **circular arc**, not a linear stem: 40 repeats close a full turn, its
    rise per repeat is zero, and it has no crystallographic repeat at all.  The same 9 deg
    defeats the ideal TGTG' glide the real crystal has.

    This is a real limit of the rigid-geometry model rather than a bug.  The published
    chain resolves the mismatch with torsions of 175 and 49 deg rather than the ideal 180
    and 60, which is a degree of freedom the ideal-angle builder does not have (the
    refinement stage and ``periodic_chain_from_torsions`` do).  It is also consistent with
    the crystallography: PVDC does not form a planar-zigzag stem.
    """
    a1, a2 = (b.backbone_angle for b in PVDC.backbone)
    assert (a1, a2) == (123.0, 114.0)
    h = helix_parameters(PVDC, np.array([0, 0]), THREE_STATE)  # all-trans
    assert h.c is None  # no commensurate repeat
    assert h.rotation_per_period == pytest.approx(a1 - a2, abs=1e-9)  # exactly the mismatch
    g = helix_parameters(PVDC, np.array([0, 1, 0, 2]), THREE_STATE)  # ideal TG+TG-
    assert g.c is None and abs(abs(g.rotation_per_period) - (a1 - a2)) < 0.02
    with pytest.raises(ValueError, match="not commensurate"):
        periodic_chain(PVDC, [0, 0], THREE_STATE)


# The glide TG+TG- chain, closed.  Torsions from a closure search over the four backbone
# dihedrals of the two-monomer repeat under the measured backbone angles; see
# ``test_measured_angles_recover_the_published_pvdc_glide_chain`` for what they mean.
GLIDE_TORSIONS = np.array([175.3, 49.4, 184.7, -49.4])


def _glide_chain():
    return periodic_chain_from_torsions(PVDC, "TG+TG-", GLIDE_TORSIONS)


def test_measured_angles_recover_the_published_pvdc_glide_chain():
    """What the corrected backbone angles buy, and it is more than they cost.

    Takahagi et al. (1988) report PVDC as a glide TGTG' chain with internal rotation
    angles of 175 deg (T) and 49 deg (G'), in a cell whose fibre axis is 4.68 A.  Ask this
    model for the torsions that make a TG+TG- repeat close (rotation error zero), holding
    the backbone angles fixed, and with the measured 123/114 it answers **175.3 and 49.4
    deg with a repeat of 4.677 A** -- the published torsion pair to within half a degree
    and the fibre repeat to 0.06%.

    With the equal 114/114 angles this package used to carry, the same closure search
    could only answer T = exactly 180 deg (no deflection at all) and G = 49.6, giving
    4.491 A, 4.0% short; and the published 175/49 did not close at all under those angles.
    So the wide CH2 angle is what *creates* the trans deflection: the 9 deg of curl per
    repeat that the unequal angles introduce is exactly what the 5 deg deflection of each
    trans bond cancels.  That is the mechanism the crystallography describes, reproduced
    here from geometry alone.
    """
    ch = _glide_chain()
    assert ch.rotation_error < 0.01  # it closes
    assert ch.c == pytest.approx(4.68, abs=0.02)  # the measured fibre repeat
    assert 174.0 < GLIDE_TORSIONS[0] < 177.0 and 47.0 < GLIDE_TORSIONS[1] < 52.0


def test_pvdc_chain_repeats_are_sane():
    """The ideal-angle conformations that *do* close still give sane helices."""
    for seq, lo, hi in (([0, 1], 6.5, 7.1), ([0, 1, 0, 1], 13.0, 14.2)):
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

    # The chain that gets packed is the closed glide, NOT an ideal-angle conformation,
    # and that is a real limitation rather than a convenience.  With the measured unequal
    # backbone angles no PVDC sequence closes at ideal torsions: every one carries the
    # 9 deg per-repeat curl, so ``periodic_chain`` refuses them all (the few the screw
    # decomposition calls commensurate are only approximately so, and fail its 1e-3 A
    # consistency check by 0.03 to 0.13 A).  What does close is the chain with the small
    # torsion deflections the real polymer has, which is what the crystal structure
    # reports and what ``periodic_chain_from_torsions`` is for.
    assert all(_ideal_chain_fails(PVDC, c.seq, m.states) for c in cands)
    ch = _glide_chain()
    pk = CrystalPacker(ch, n_chains=2)
    e = pk.energy(np.array([[6.0, 10.0, 90.0, 30.0, 200.0, 1.0, 0]]))
    assert np.isfinite(e).all()


def _ideal_chain_fails(polymer, seq, states) -> bool:
    """True when ``seq`` has no usable ideal-angle crystallographic repeat."""
    try:
        periodic_chain(polymer, seq, states)
    except (ValueError, RuntimeError):
        return True
    return False


def test_pvdc_flip_is_free_beyond_the_cutoff():
    """The bonded-exclusion fix must hold for a new chemistry too."""
    ch = _glide_chain()
    pk = CrystalPacker(ch, n_chains=2)
    far = np.array([[60.0, 60.0, 90.0, 17.0, 143.0, 0.9, 0.0]])
    anti = far.copy()
    anti[0, 6] = 1.0
    assert float(pk.energy(anti)[0]) == pytest.approx(float(pk.energy(far)[0]), abs=1e-9)


def test_all_trans_pvdc_is_strained_relative_to_pvdf():
    """Documents the Phase 1 finding, and what correcting the backbone angles did to it.

    Chlorine is bulky (UFF x_i 3.95 A) and its C-Cl bond long (1.77 A), so with bond
    angles frozen the planar zigzag cannot relieve Cl...Cl contact.  With both angles
    frozen at 114 deg -- PVDF's value, which an earlier version of ``polymers.py``
    borrowed without justification -- that cost 313 kcal/mol on a ten-bond oligomer
    against PVDF's 11.

    The audit predicted that this was mostly an artifact of the wrong angles, since the
    published structure opens C-CH2-C to 123 deg precisely to relieve that contact.  It
    was: giving the model the measured 123/114 takes the strain to **74 kcal/mol**, a
    factor of 4.2, and takes the drop available by rotating every bond away from trans
    from 79 kcal/mol to 4.0, so all-trans goes from grossly unphysical to nearly
    metastable.  Both halves are asserted here -- the counterfactual as well as the
    current value -- because the comparison is the finding.

    What it does *not* do is make all-trans a good RIS reference: PVDC's crystal
    conformation is a glide TGTG', not a planar zigzag.  Variable bond angles remain the
    principled fix.
    """
    from polyfind.polymers import PVDF, BackboneAtom, Polymer

    ff = SimpleFF()
    equal = Polymer(  # the counterfactual: PVDF's angles borrowed, as this used to be
        name="pvdc-equal-angles", formula="-(CH2-CCl2)n-", bond_length=PVDC.bond_length,
        backbone=(BackboneAtom("C", "H", 1.09, 114.0, 108.0, -0.20, +0.10),
                  BackboneAtom("C", "Cl", 1.77, 114.0, 110.0, +0.20, -0.10)),
    )
    trans = np.full(10, 180.0)
    lj = {p.name: ff.components(build_chain(p, trans))["lj"] for p in (PVDF, PVDC, equal)}
    assert lj["pvdf"] < 50.0
    assert lj["pvdc-equal-angles"] > 200.0  # the old, borrowed-angle number
    assert 50.0 < lj["pvdc"] < 150.0  # the measured angles remove three quarters of it
    assert lj["pvdc"] < 0.3 * lj["pvdc-equal-angles"]
    # and all-trans is far closer to a minimum than it was: the best rigid rotation away
    # from it now gains 4 kcal/mol rather than 79
    e_t = float(ff.energy(build_chain(PVDC, trans)))
    best = min(float(ff.energy(build_chain(PVDC, np.full(10, float(p))))) for p in range(-180, 181, 30))
    assert -10.0 < best - e_t < 0.0
