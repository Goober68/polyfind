"""Phase 3 of the chemistry extension: pendants that are small rigid fragments.

AN (CH2-CH(CN)), VDCN (CH2-C(CN)2) and FANOME (CH2-C(CN)(OCH3)) need a pendant to be a
group of atoms rather than one atom.  Covered here:

1. That the generalisation is invisible to everything that already worked -- PE, PVDF,
   PVDC, CFE and CDFE must build byte-identical coordinates and charges (digests taken
   from the code immediately before the phase-3 change), exactly as phase 2 asserted for
   the three that predated it.
2. That the fragments are built as specified: every atom at its stated distance and
   angle, nitriles actually linear, the methoxy at its stated frozen rotamer, a neutral
   repeat and no overlapping atoms inside a monomer -- through the single builder, the
   batched one and the backbone-angle override.
3. That ``Polymer.is_chiral`` is right for each: AN and FANOME carry two *different*
   pendants on one backbone carbon and are chiral, VDCN carries two identical ones and
   is not.  A multi-atom pendant does not by itself make a stereocentre.
4. That the funnel runs end to end on VDCN and AN.

As in ``test_pvdc.py`` and ``test_cfe_cdfe.py``, nothing is asserted about the *values*
of the fitted energies.  Quite the opposite: ``test_all_trans_is_a_bad_reference_state``
measures how badly strained the all-trans chain of each of these is and pins the finding
down, because with frozen bond angles a pendant that reaches 2.6 A cannot get out of its
own 1-3 neighbour's way.  That is a property of the rigid-geometry model, and it is
recorded rather than tuned away.
"""
import hashlib
import warnings

import numpy as np
import pytest

from polyfind.chain import (
    angle,
    build_chain,
    build_chain_batch,
    dihedral,
    distance,
    pendant_positions,
    set_dihedrals,
    substituent_positions,
)
from polyfind.enumerate import enumerate_periodic
from polyfind.forcefield import SimpleFF, fit_ris
from polyfind.pack import CrystalPacker, periodic_chain
from polyfind.polymers import (
    AN,
    CDFE,
    CFE,
    FANOME,
    METHOXY,
    NITRILE,
    PE,
    PVDC,
    PVDF,
    THREE_STATE,
    VDCN,
    BackboneAtom,
    Pendant,
    PendantAtom,
    get_polymer,
    lj_params,
    methoxy,
    nitrile,
    single_atom,
)

NEW = [AN, VDCN, FANOME]


# --------------------------------------------------------------- 1. no regression
def _digest(*arrays) -> str:
    h = hashlib.sha256()
    for a in arrays:
        a = np.ascontiguousarray(a)
        h.update(str(a.shape).encode())
        h.update(a.tobytes())
    return h.hexdigest()


# coordinates+charges of each polymer, taken before a pendant could be a fragment.
# The PE/PVDF/PVDC entries are the same strings as in ``test_cfe_cdfe.py``, i.e. they
# have now survived two changes to the pendant model unchanged.
BASELINE = {
    ("pe", True): "9fa3fd6b6d6b909c7570164f93b72b11758beb29d96a0f69b385d3705690ecdd",
    ("pe", False): "539a09e1cf75594f19c3321883bea185ff5b908f8da9ae9278873d3c238e01dd",
    ("pe", "batch"): "15c45ec4a2971fe5337ba834cd750c875d7122bd92c26fffa18448d21440d8a7",
    ("pe", "trans"): "8ae1ecce71f26cecfc20b71ac379ad11a3fa38c4eeaa67d8c64488d9d4cbb707",
    ("pvdf", True): "9650a4c802ac61a875e46a07f61aa6389c2c8baabc64d27dd2a82db05a32a0f6",
    ("pvdf", False): "ad46e666db004c26a2cdd0dce952d31dbba5a708ef1bf1ba6e8b234c076214a0",
    ("pvdf", "batch"): "796284f466d358b2d41373e776988d0799daea36ce666bc0a265d26cf55e16a5",
    ("pvdf", "trans"): "b51d418a7a926b5ab8e149ba8d262ae9446bb50c9a7027a337aa748305927e50",
    ("pvdc", True): "81b85857f57c855c014f13e75401424e2720aa9824870304113eb05778ac5f8e",
    ("pvdc", False): "cb6ed4d7838f9445fd81b3807240700b57f3cddb60ad6e060c2f5e3c1d8b4f82",
    ("pvdc", "batch"): "58f5b04d2d04af12c242dcf388ca06a6066ae4c005d2e9bc656243540cfa1f37",
    ("pvdc", "trans"): "331b546c3025fb8c9c23b099fe6068f62580885f523c574af0c6b52ca8888eec",
    ("cfe", True): "d03823aee5aef827e56aa74bbca230c8c8586ac47ec044864980ce1b63acefaf",
    ("cfe", False): "52d962b948dc69081b36e76f0d0eeb49be010a1f105221c06fb55953d0c91dc2",
    ("cfe", "batch"): "170c07f199cf420753aaa8f96208ef050dad0bca5dde349e3dc328e25cb7245b",
    ("cfe", "trans"): "44efc3fcb66e1095cffd6f229a6fe222c75b6f135779db24d4f75ed841c49d45",
    ("cdfe", True): "eee63607af900f92c2b5d8bdf2fc8581ecde1ba95d86701bc460366fb87b8356",
    ("cdfe", False): "79e794cdfc46a57f488d16aaf3a5a49183c0ef877a41b8e6ea01b71de337d98d",
    ("cdfe", "batch"): "6b714c74c95bad8d04c2c577f2a4a5ffddd77b5be82066659854c0156b06f9b7",
    ("cdfe", "trans"): "ffa5f5db586ba04af064dcab53949b190ad332574d71eb04e61dea8580b516c8",
}
ELEMENTS = {
    "pe": "CHHCHHCHHCHHCHHCHHCHHCHHCHHCHHCHHCHHCHHHH",
    "pvdf": "CHHCFFCHHCFFCHHCFFCHHCFFCHHCFFCHHCFFCHHHH",
    "pvdc": "CHHCClClCHHCClClCHHCClClCHHCClClCHHCClClCHHCClClCHHHH",
    "cfe": "CHHCFClCHHCFClCHHCFClCHHCFClCHHCFClCHHCFClCHHHH",
    "cdfe": "CHClCFFCHClCFFCHClCFFCHClCFFCHClCFFCHClCFFCHClHH",
}


@pytest.mark.parametrize("polymer", [PE, PVDF, PVDC, CFE, CDFE])
def test_existing_polymers_are_byte_identical(polymer):
    """Bit-for-bit: a single-atom pendant must still take the single-atom code path."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        dih = np.random.default_rng(1234).uniform(-180.0, 180.0, size=12)
        for cap in (True, False):
            s = build_chain(polymer, dih, cap=cap)
            assert _digest(s.coords, s.charges) == BASELINE[(polymer.name, cap)]
        dihs = np.random.default_rng(99).uniform(-180.0, 180.0, size=(6, 10))
        tpl, coords = build_chain_batch(polymer, dihs, cap=True)
        assert _digest(coords, tpl.charges) == BASELINE[(polymer.name, "batch")]
        s = build_chain(polymer, np.full(10, 180.0), cap=True)
        assert _digest(s.coords, s.charges) == BASELINE[(polymer.name, "trans")]
    assert "".join(s.elements) == ELEMENTS[polymer.name]
    assert polymer.atoms_per_repeat == 3 * polymer.bonds_per_repeat  # all single atoms


def test_a_single_atom_pendant_is_the_degenerate_fragment():
    """The old spelling still means the old thing, and now means a one-atom fragment.

    This is the backwards-compatibility contract in one place: a plain element symbol
    plus a bond length and a charge becomes a :class:`Pendant` of one atom sitting on
    the pendant bond axis, which is exactly where ``substituent_positions`` puts it.
    """
    spec = BackboneAtom("C", ("F", "Cl"), (1.35, 1.77), 114.0, 108.0, +0.30, (-0.20, -0.10))
    assert spec.pendants == (single_atom("F", 1.35, -0.20), single_atom("Cl", 1.77, -0.10))
    assert [len(p) for p in spec.pendants] == [1, 1]
    assert spec.substituents == ("F", "Cl")
    assert spec.sub_bonds == (1.35, 1.77) and spec.sub_charges == (-0.20, -0.10)
    assert spec.n_pendant_atoms == 2 and spec.is_stereocentre

    p, x, n = np.array([-1.2, 0.9, 0.0]), np.zeros(3), np.array([1.3, 0.8, 0.2])
    s1, s2 = substituent_positions(p, x, n, spec.sub_bond, spec.sub_angle)
    (g1,), (g2,) = pendant_positions(p, x, n, spec)
    assert g1.tobytes() == s1.tobytes() and g2.tobytes() == s2.tobytes()


def test_a_fragment_pendant_still_reports_a_numeric_bond_and_charge():
    """``pack`` and ``linegroup`` forward ``spec.sub_bond`` straight into
    ``substituent_positions``, so it has to stay a number (or a pair of numbers) even
    when the pendant carries its own geometry.  The fragment fills it in."""
    spec = FANOME.backbone[1]
    assert spec.sub_bond == (NITRILE.bond, METHOXY.bond) == (1.47, 1.41)
    assert spec.sub_charge == pytest.approx((-0.15, -0.05))
    assert spec.substituents == ("C", "O")  # the atom bonded to the backbone
    assert VDCN.backbone[1].sub_bond == 1.47  # both pendants alike: still a scalar

    # and the first atom of each pendant is exactly where substituent_positions puts it
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        st = build_chain(FANOME, np.full(6, 180.0), cap=False)
    k = 1
    idx, prev, nxt = (int(st.backbone[j]) for j in (k, k - 1, k + 1))
    s1, s2 = substituent_positions(st.coords[prev], st.coords[idx], st.coords[nxt],
                                   spec.sub_bond, spec.sub_angle)
    first_of_each = [st.subs_of[idx][0], st.subs_of[idx][len(spec.pendants[0])]]
    np.testing.assert_allclose(st.coords[first_of_each[0]], s1, atol=1e-12)
    np.testing.assert_allclose(st.coords[first_of_each[1]], s2, atol=1e-12)


def test_a_bad_pendant_spec_is_rejected():
    with pytest.raises(ValueError, match="bond axis"):
        Pendant(atoms=(PendantAtom("C", 0.0, (1.4, 0.1, 0.0)),))
    with pytest.raises(ValueError, match="at least one atom"):
        Pendant(atoms=())
    with pytest.raises(ValueError, match="intra-fragment bond"):
        Pendant(atoms=(PendantAtom("C", 0.0, (1.4, 0.0, 0.0)),), bonds=((0, 3),))
    with pytest.raises(ValueError, match="contradicts pendant"):
        BackboneAtom("C", NITRILE, 1.09, 114.0, 108.0, 0.0, None)
    with pytest.raises(ValueError, match="needs its own"):
        BackboneAtom("C", "H", None, 114.0, 108.0, 0.0, +0.10)


# ------------------------------------------------------------- 2. the new chemistries
def test_new_polymers_are_registered():
    assert get_polymer("an") is AN and get_polymer("vdcn") is VDCN
    assert get_polymer("fanome") is FANOME
    # backbone atom + pendant atoms, no longer 3 per backbone atom
    assert (AN.atoms_per_repeat, VDCN.atoms_per_repeat, FANOME.atoms_per_repeat) == (7, 8, 11)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for p in NEW:
            assert p.bonds_per_repeat == 2
            st = build_chain(p, np.full(8, 180.0), cap=False)  # 11 backbone atoms
            assert st.n_atoms == 5 * p.atoms_per_repeat + (1 + p.backbone[0].n_pendant_atoms)


def test_nitrogen_and_oxygen_are_parameterised():
    """UFF (Rappe et al. 1992), the same source as the C/H/F/Cl entries already there."""
    x, d = lj_params(["N", "O"])
    np.testing.assert_allclose(x, [3.660, 3.500])
    np.testing.assert_allclose(d, [0.069, 0.060])
    st = build_chain(FANOME, np.full(4, 180.0))
    lj_params(st.elements)  # every element of the new chemistries is scorable


@pytest.mark.parametrize("polymer", NEW)
def test_the_pendant_bonds_and_the_angle_between_them(polymer):
    """Each pendant's first atom at its own bond length, at the specified angle apart."""
    st = build_chain(polymer, np.full(8, 180.0), cap=False)
    B = polymer.bonds_per_repeat
    for k, idx in enumerate(st.backbone):
        spec = polymer.backbone[k % B]
        idx = int(idx)
        n1 = len(spec.pendants[0])
        i1, i2 = st.subs_of[idx][0], st.subs_of[idx][n1]
        b1, b2 = spec.sub_bonds
        assert distance(st.coords, idx, i1) == pytest.approx(b1, abs=1e-9)
        assert distance(st.coords, idx, i2) == pytest.approx(b2, abs=1e-9)
        assert angle(st.coords, i1, idx, i2) == pytest.approx(spec.sub_angle, abs=1e-6)
        assert (st.elements[i1], st.elements[i2]) == spec.substituents


@pytest.mark.parametrize("polymer", NEW)
def test_the_nitrile_is_linear_and_at_its_stated_lengths(polymer):
    """The easy fragment: two atoms collinear with the pendant bond, no new freedom.

    Linear means linear -- the nitrogen is on the pendant bond axis in every
    conformation, because both fragment atoms are placed along it rather than by an
    internal angle of their own.  Tested as exact collinearity (a vanishing cross
    product, and the two bond lengths summing to the C...N distance) rather than as an
    angle, since ``arccos`` near 180 deg is only good to about 1e-6 deg.
    """
    dih = np.random.default_rng(7).uniform(-180, 180, size=8)
    st = build_chain(polymer, dih, cap=False)
    B = polymer.bonds_per_repeat
    seen = 0
    for k, idx in enumerate(st.backbone):
        idx = int(idx)
        offset = 0
        for pendant in polymer.backbone[k % B].pendants:
            if pendant.label == "CN":
                c, n = st.subs_of[idx][offset], st.subs_of[idx][offset + 1]
                assert (st.elements[c], st.elements[n]) == ("C", "N")
                assert distance(st.coords, idx, c) == pytest.approx(1.47, abs=1e-9)
                assert distance(st.coords, c, n) == pytest.approx(1.16, abs=1e-9)
                assert distance(st.coords, idx, n) == pytest.approx(1.47 + 1.16, abs=1e-9)
                v1 = st.coords[c] - st.coords[idx]
                v2 = st.coords[n] - st.coords[c]
                assert np.linalg.norm(np.cross(v1, v2)) < 1e-12
                assert angle(st.coords, idx, c, n) == pytest.approx(180.0, abs=1e-4)
                seen += 1
            offset += len(pendant)
    assert seen >= 4


def test_the_methoxy_geometry_and_its_frozen_rotamer():
    """The hard fragment: right geometry, and the frozen rotamer is the stated one.

    The C-O torsion is a real degree of freedom that this model cannot carry, so it is
    frozen with the methyl anti to the other pendant, and the methyl itself is frozen
    staggered.  Both are asserted here so that a later change to the rotamer is a test
    failure and not a silent change of chemistry.

    The methyl dihedrals are checked up to an overall sign because pendant 2's frame is
    the mirror image of pendant 1's (see ``chain.pendant_frames``): a fragment placed as
    the second pendant is reflected, which for the methoxy -- mirror-symmetric about its
    own plane -- only swaps the two out-of-plane hydrogens.
    """
    dih = np.random.default_rng(8).uniform(-180, 180, size=8)
    st = build_chain(FANOME, dih, cap=False)
    for k, idx in enumerate(st.backbone):
        if k % 2 == 0:
            continue
        idx = int(idx)
        nit_c, _n, o, cme, h1, h2, h3 = (int(i) for i in st.subs_of[idx])
        assert [st.elements[i] for i in (o, cme, h1, h2, h3)] == ["O", "C", "H", "H", "H"]
        assert distance(st.coords, idx, o) == pytest.approx(1.41, abs=1e-9)
        assert distance(st.coords, o, cme) == pytest.approx(1.43, abs=1e-9)
        assert angle(st.coords, idx, o, cme) == pytest.approx(111.5, abs=1e-6)
        # frozen C-O torsion: the methyl carbon is anti to the other pendant
        assert abs(dihedral(st.coords, nit_c, idx, o, cme)) == pytest.approx(180.0, abs=1e-6)
        # frozen, staggered methyl
        psis = []
        for h in (h1, h2, h3):
            assert distance(st.coords, cme, h) == pytest.approx(1.09, abs=1e-9)
            assert angle(st.coords, o, cme, h) == pytest.approx(109.5, abs=1e-6)
            psis.append(dihedral(st.coords, idx, o, cme, h))
        assert sorted(abs(p) for p in psis) == pytest.approx([60.0, 60.0, 180.0], abs=1e-6)
        assert psis[1] * psis[2] < 0  # one gauche each way


@pytest.mark.parametrize("polymer", NEW)
def test_the_repeat_is_neutral_and_a_monomer_does_not_overlap_itself(polymer):
    st = build_chain(polymer, np.full(8, 180.0), cap=True)
    assert abs(float(st.charges.sum())) < 1e-12
    for spec in polymer.backbone:  # each backbone atom's own group is neutral too
        assert abs(spec.charge + sum(spec.sub_charges)) < 1e-12
        assert all(abs(p.charge - sum(a.charge for a in p.atoms)) < 1e-12 for p in spec.pendants)
    d = np.linalg.norm(st.coords[:, None] - st.coords[None], axis=-1)
    for idx in st.backbone:
        group = [int(idx)] + [int(i) for i in st.subs_of[int(idx)]]
        sub = d[np.ix_(group, group)] + np.eye(len(group)) * 9
        assert sub.min() > 1.0  # nothing in a monomer sits on top of anything else


@pytest.mark.parametrize("polymer", [AN, VDCN])
def test_nothing_overlaps_along_an_all_trans_chain(polymer):
    """True for AN and VDCN.  It is *not* true for FANOME; see the strain test."""
    st = build_chain(polymer, np.full(8, 180.0), cap=True)
    d = np.linalg.norm(st.coords[:, None] - st.coords[None], axis=-1) + np.eye(st.n_atoms) * 9
    assert d.min() > 1.0


@pytest.mark.parametrize("polymer", NEW)
@pytest.mark.parametrize("cap", [True, False])
def test_batched_builder_matches_the_single_one(polymer, cap):
    dih = np.random.default_rng(3).uniform(-180.0, 180.0, size=(5, 8))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        template, coords = build_chain_batch(polymer, dih, cap=cap)
        for m in range(5):
            single = build_chain(polymer, dih[m], cap=cap)
            assert single.elements == template.elements and single.bonds == template.bonds
            assert single.subs_of == template.subs_of
            np.testing.assert_allclose(coords[m], single.coords, atol=1e-9)


@pytest.mark.parametrize("polymer", NEW)
def test_bond_angle_override_still_works(polymer):
    """The per-repeat backbone-angle override must survive the fragment change too."""
    dih = np.random.default_rng(5).uniform(-180.0, 180.0, size=(3, 8))
    angles = np.array([118.0, 110.0])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        template, coords = build_chain_batch(polymer, dih, cap=True, bond_angles=angles)
        for m in range(3):
            single = build_chain(polymer, dih[m], cap=True, bond_angles=angles)
            np.testing.assert_allclose(coords[m], single.coords, atol=1e-9)
        st = build_chain(polymer, dih[0], cap=True, bond_angles=angles)
        frozen = build_chain(polymer, dih[0], cap=True)
    assert not np.allclose(st.coords, frozen.coords)
    # the whole fragment followed the opened-up backbone rigidly: internal geometry intact
    for k, idx in enumerate(st.backbone):
        spec = polymer.backbone[k % polymer.bonds_per_repeat]
        idx, offset = int(idx), 0
        for pendant in spec.pendants:
            atoms = st.subs_of[idx][offset:offset + len(pendant)]
            assert distance(st.coords, idx, atoms[0]) == pytest.approx(pendant.bond, abs=1e-9)
            for i, j in pendant.bonds:
                d0 = np.linalg.norm(np.subtract(pendant.atoms[i].offset, pendant.atoms[j].offset))
                assert distance(st.coords, atoms[i], atoms[j]) == pytest.approx(d0, abs=1e-9)
            offset += len(pendant)


@pytest.mark.parametrize("polymer", NEW)
def test_the_bond_graph_carries_the_fragment(polymer):
    """The nonbonded exclusions come from the bond graph, so the fragment's own bonds
    have to be in it -- otherwise a nitrile nitrogen would be scored against the carbon
    it is triple-bonded to."""
    st = build_chain(polymer, np.full(6, 180.0), cap=False)
    bonds = {tuple(sorted(b)) for b in st.bonds}
    B = polymer.bonds_per_repeat
    for k, idx in enumerate(st.backbone):
        idx, offset = int(idx), 0
        for pendant in polymer.backbone[k % B].pendants:
            atoms = st.subs_of[idx][offset:offset + len(pendant)]
            assert tuple(sorted((idx, atoms[0]))) in bonds
            for i, j in pendant.bonds:
                assert tuple(sorted((atoms[i], atoms[j]))) in bonds
            offset += len(pendant)
    assert len(bonds) == len(st.bonds)  # no duplicates


@pytest.mark.parametrize("polymer", NEW)
def test_set_dihedrals_still_infers_the_cap(polymer):
    """It used to count atoms as ``3 * (N + 3)``, which is wrong for a fragment."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for cap in (True, False):
            st = build_chain(polymer, np.full(6, 180.0), cap=cap)
            again = set_dihedrals(st, np.full(6, 60.0))
            assert again.n_atoms == st.n_atoms and again.elements == st.elements


# ------------------------------------------------------------------- 3. chirality
def test_is_chiral_is_right_for_every_polymer():
    """Two different pendants make a stereocentre; two identical ones do not, and a
    multi-atom pendant does not by itself make one.  AN (H and a nitrile) and FANOME (a
    nitrile and a methoxy) are chiral; VDCN (two nitriles) is not, for exactly the same
    reason PVDF is not."""
    assert AN.is_chiral and FANOME.is_chiral
    assert not VDCN.is_chiral
    assert [p.is_chiral for p in (PE, PVDF, PVDC, CFE, CDFE)] == [False, False, False, True, True]
    assert VDCN.backbone[1].pendants[0] == VDCN.backbone[1].pendants[1] == NITRILE
    assert AN.backbone[1].pendants != (NITRILE, NITRILE)
    assert nitrile() == NITRILE and methoxy() == METHOXY  # value equality, not identity


@pytest.mark.parametrize("polymer,chiral", [(VDCN, False), (AN, True), (FANOME, True)])
def test_the_potential_agrees_with_is_chiral(polymer, chiral):
    """Measured, not asserted from the definition.

    Reflecting an achiral chain gives the same molecule, so ``E(-phi) == E(phi)``.
    Reflecting a chiral one gives its enantiomer and the energies genuinely differ; what
    survives is reflection composed with chain reversal, which holds for both.
    """
    dih = np.array([180.0, 60.0, 180.0, -60.0, 180.0, 180.0, 60.0, 180.0])  # not a palindrome
    ff = SimpleFF()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        e = ff.energy(build_chain(polymer, dih))
        rev = ff.energy(build_chain(polymer, -dih[::-1]))
        plain = ff.energy(build_chain(polymer, -dih))
    assert abs(e - rev) < 1e-6 * max(1.0, abs(e))
    assert (abs(e - plain) > 1.0) if chiral else (abs(e - plain) < 1e-8)


def test_the_fit_uses_the_chirality_it_is_told():
    """``symmetrize='auto'`` mirrors VDCN, which is exact, and must not mirror AN."""
    ff = SimpleFF()
    mir = np.array(THREE_STATE.mirror)

    def residual(model):
        e2 = model.second_order
        return float(np.abs(e2 - e2[:, mir][:, :, mir]).max())

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        assert residual(fit_ris(VDCN, ff, step=45.0, n_monomers=3, third_order=False).model) < 1e-9
        assert residual(fit_ris(AN, ff, step=45.0, n_monomers=3, third_order=False).model) > 1.0


# ------------------------------------------------------- 4. the finding, and the funnel
def test_all_trans_is_a_bad_reference_state_for_all_three():
    """The phase-1 PVDC finding, made worse by bulkier pendants.  Recorded, not tuned.

    With bond angles frozen, a planar zigzag cannot relieve contact between substituents
    on 1-3 backbone atoms.  PVDF's all-trans chain carries about 9 kcal/mol of
    Lennard-Jones strain and PVDC's about 313; a nitrile reaches 2.63 A from the backbone
    and a methoxy 2.35 A, against 1.77 A for chlorine, so:

        AN      ~90 kcal/mol   (one nitrile per repeat)
        VDCN   ~176 kcal/mol   (two)
        FANOME  ~1e6           (a nitrile and a methoxy; two methyl hydrogens on
                                consecutive substituted carbons end up 0.80 A apart)

    AN and VDCN are strained but finite, and comparable with the CFE/CDFE numbers the
    pipeline is already run on; FANOME's all-trans chain is not a physical structure at
    all.  For none of the three is all-trans a sensible RIS reference: rotating away from
    it lowers the energy in every case.  FANOME's number is robust to the frozen methoxy
    rotamer -- every azimuth of the C-O torsion gives between 7e3 and 2e8 kcal/mol -- so
    it is the rigid backbone angles, not the rotamer choice, that produces it.  Variable
    backbone angles (and probably a reference state other than all-trans) are the fix.
    """
    ff = SimpleFF()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        lj = {p.name: ff.components(build_chain(p, np.full(10, 180.0)))["lj"]
              for p in (PVDF, PVDC, AN, VDCN, FANOME)}
        # every one of them is lowered by rotating away from all-trans
        best = {}
        for p in NEW:
            e_trans = ff.energy(build_chain(p, np.full(10, 180.0)))
            best[p.name] = min(ff.energy(build_chain(p, np.full(10, float(phi))))
                               for phi in range(-180, 181, 30)) - e_trans
    assert lj["pvdf"] < 50.0 and lj["pvdc"] > 200.0  # the phase-1 baseline, unchanged
    assert 50.0 < lj["an"] < 200.0
    assert 100.0 < lj["vdcn"] < 400.0
    assert lj["fanome"] > 1e5  # an outright atom-atom overlap, not merely close contact
    assert all(v < -10.0 for v in best.values()), best  # all-trans is not even a minimum


@pytest.mark.parametrize("polymer", [VDCN, AN])
def test_the_whole_funnel_runs(polymer):
    """Cheap fit, enumeration, a periodic chain and one packing energy.

    FANOME is left out on purpose: its all-trans reference is an overlapping structure
    (see above), so running the fit on it would produce numbers with no meaning at all.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        rep = fit_ris(polymer, SimpleFF(), step=45.0, n_monomers=3, third_order=False)
    m = rep.model
    assert m.B == 2 and m.second_order.shape == (2, 3, 3)
    assert np.isfinite(m.first_order).all() and np.isfinite(m.second_order).all()

    cands = enumerate_periodic(polymer, m, max_period=4, k_per_period=10)
    assert cands and all(c.helix.rise_per_bond > 0.3 for c in cands)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ch = periodic_chain(polymer, m.parse("TT"), m.states)
    assert len(ch.elements) == polymer.atoms_per_repeat * ch.n_monomers
    assert set(ch.elements) >= {"C", "H", "N"}  # the nitrogens made it into the cell
    pk = CrystalPacker(ch, n_chains=2)
    e = pk.energy(np.array([[6.5, 10.0, 90.0, 30.0, 200.0, 1.0, 0]]))
    assert np.isfinite(e).all()


def test_packing_and_refinement_paths_accept_multi_atom_pendants():
    """The block builders used to index atoms as 3k+{0,1,2}, so they only ever
    handled single-atom pendants. Densities also needed masses for N and O."""
    import warnings as _w

    import numpy as np

    from polyfind.pack import CrystalPacker, periodic_chain, repeat_chains_from_torsions
    from polyfind.polymers import THREE_STATE, get_polymer

    with _w.catch_warnings():
        _w.simplefilter("ignore")
        for name, n_atoms in (("vdcn", 8), ("an", 7)):
            poly = get_polymer(name)
            ch = periodic_chain(poly, [0, 0], THREE_STATE)
            assert ch.n_atoms == n_atoms
            assert ch.mass > 0  # needs N in the mass table
            pk = CrystalPacker(ch, n_chains=2)
            assert pk.density(7.0, 11.0, 90.0) > 0

            # the batched torsion path feeds refinement's gradients
            tors = np.tile(ch.dihedrals, (3, 1)) + np.array([[0.0], [2.0], [-2.0]])
            chains = repeat_chains_from_torsions(poly, ch.name, tors)
            assert len(chains) == 3
            assert all(c.n_atoms == n_atoms for c in chains)
            # the unperturbed row must reproduce the original chain exactly
            assert np.allclose(chains[0].coords, ch.coords, atol=1e-9)

            p = np.array([[7.0, 11.0, 90.0, 20.0, 200.0, 1.0, 0]])
            e = pk.energy(np.repeat(p, 3, axis=0), coords=np.stack([c.coords for c in chains]),
                          c=np.array([c.c for c in chains]))
            assert np.isfinite(e).all()
