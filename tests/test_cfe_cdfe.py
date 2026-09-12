"""Phase 2 of the chemistry extension: backbone atoms with two *different* pendants.

Covers three things:

1. That the change is invisible to the chemistries that already worked -- PE, PVDF and
   PVDC must build byte-identical coordinates and charges (digests below were taken
   from the code immediately before the phase-2 change).
2. That CFE (CH2-CFCl) and CDFE (CHCl-CF2) are built as specified: each pendant at its
   own bond length, the specified angle between them, a neutral repeat, no overlaps --
   through every construction path, including the batched builder and the coordinate
   builders in ``pack``/``linegroup`` that forward ``spec.sub_bond`` themselves.
3. That the chain these produce is the *isotactic* one, and therefore that ``fit_ris``'s
   ``symmetrize=True`` (G+/G- averaging) is no longer valid for them.  The tests assert
   the inequivalence directly rather than trusting the argument.

As with ``test_pvdc.py`` these assert nothing about the *values* of the fitted energies:
chlorine plus frozen bond angles makes the all-trans reference badly strained (the phase
1 finding), so the numbers are not meaningful even though the machinery runs.
"""
import hashlib
import warnings

import numpy as np
import pytest

import polyfind.chain as chain_mod
from polyfind.chain import build_chain, build_chain_batch, distance
from polyfind.enumerate import enumerate_periodic
from polyfind.forcefield import SimpleFF, fit_ris
from polyfind.pack import CrystalPacker, periodic_chain, repeat_chains_from_torsions
from polyfind.polymers import CDFE, CFE, PE, PVDC, PVDF, THREE_STATE, get_polymer


# --------------------------------------------------------------- 1. no regression
def _digest(*arrays) -> str:
    h = hashlib.sha256()
    for a in arrays:
        a = np.ascontiguousarray(a)
        h.update(str(a.shape).encode())
        h.update(a.tobytes())
    return h.hexdigest()


# coordinates+charges of each polymer, taken before BackboneAtom grew a second pendant.
# The PVDF and PVDC entries were RE-RECORDED when the batched geometry corrections landed
# (PVDF C-C 1.54 -> 1.528 A, PVDC backbone angles 114/114 -> 123/114; docs/REFERENCES.md).
# The PE entries are still the originals, and they are the ones that carry this test's
# actual claim: the single-atom pendant path is unchanged.  PVDF's and PVDC's digests
# moved because their declared geometry moved, not because the builder did.
BASELINE = {
    ("pe", True): "9fa3fd6b6d6b909c7570164f93b72b11758beb29d96a0f69b385d3705690ecdd",
    ("pe", False): "539a09e1cf75594f19c3321883bea185ff5b908f8da9ae9278873d3c238e01dd",
    ("pe", "batch"): "15c45ec4a2971fe5337ba834cd750c875d7122bd92c26fffa18448d21440d8a7",
    ("pe", "trans"): "8ae1ecce71f26cecfc20b71ac379ad11a3fa38c4eeaa67d8c64488d9d4cbb707",
    ("pvdf", True): "e251321676b51f72cf55cf87cfeec1fee84e59bdb259ef99196cf984f7f1d0ff",
    ("pvdf", False): "4671ab3ef1ff53a94b64aa483afb845751472d63b1bab417c418c9b2cd799003",
    ("pvdf", "batch"): "8ef0a6cf915c7e274169ed3ea14322001bce5c12dab50596ef4fb01ea5b3d6dd",
    ("pvdf", "trans"): "a1455e5c19761c6d3e6bd389d8756bbef838b42a176e6539d5252ebb6c889585",
    ("pvdc", True): "b1e01e4f429ac8e679454cc684951233908b325fa64dde8134d9b44923363a23",
    ("pvdc", False): "11720a0fef5c09503549c6f05c3b7274adf87198dc4cf82c2c4b890a4dceeb95",
    ("pvdc", "batch"): "b3dbf43408b3ea6fcfee4aca1397dee0e36510edaabf5eeea6fa93457551eda4",
    ("pvdc", "trans"): "7026baa2e03ff0c4fc488c01045aa8ecc41643df7fdfe89d1f37a8f6401d0f84",
}
ELEMENTS = {
    "pe": "CHHCHHCHHCHHCHHCHHCHHCHHCHHCHHCHHCHHCHHHH",
    "pvdf": "CHHCFFCHHCFFCHHCFFCHHCFFCHHCFFCHHCFFCHHHH",
    "pvdc": "CHHCClClCHHCClClCHHCClClCHHCClClCHHCClClCHHCClClCHHHH",
}


@pytest.mark.parametrize("polymer", [PE, PVDF, PVDC])
def test_existing_polymers_are_byte_identical(polymer):
    """Bit-for-bit, not approximately: the symmetric case must take the same code path."""
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


# ------------------------------------------------------------- 2. the new chemistries
def test_new_polymers_are_registered():
    assert get_polymer("cfe") is CFE and get_polymer("cdfe") is CDFE
    for p in (CFE, CDFE):
        assert p.bonds_per_repeat == 2 and p.atoms_per_repeat == 6
        assert p.is_chiral
    assert not any(p.is_chiral for p in (PE, PVDF, PVDC))


@pytest.mark.parametrize(
    "polymer,expected",
    [(CFE, {"H": 2, "F": 1, "Cl": 1}), (CDFE, {"H": 1, "F": 2, "Cl": 1})],
)
def test_repeat_unit_composition(polymer, expected):
    st = build_chain(polymer, np.full(8, 180.0), cap=False)  # 11 backbone atoms
    counts = {e: sum(1 for x in st.elements if x == e) for e in ("H", "F", "Cl")}
    # 11 backbone atoms = 5 whole repeats plus one extra atom of the first type
    extra = {}
    for e in polymer.backbone[0].substituents:
        extra[e] = extra.get(e, 0) + 1
    for e in ("H", "F", "Cl"):
        assert counts[e] == 5 * expected.get(e, 0) + extra.get(e, 0)
    assert sum(1 for x in st.elements if x == "C") == 11


@pytest.mark.parametrize("polymer", [CFE, CDFE])
def test_each_pendant_sits_at_its_own_bond_length_and_angle(polymer):
    st = build_chain(polymer, np.full(8, 180.0), cap=False)
    B = polymer.bonds_per_repeat
    seen_asymmetric = False
    for k, idx in enumerate(st.backbone):
        spec = polymer.backbone[k % B]
        i1, i2 = st.subs_of[int(idx)]
        assert (st.elements[i1], st.elements[i2]) == spec.substituents
        assert st.charges[i1] == spec.sub_charges[0] and st.charges[i2] == spec.sub_charges[1]
        b1, b2 = spec.sub_bonds
        assert distance(st.coords, int(idx), i1) == pytest.approx(b1, abs=1e-9)
        assert distance(st.coords, int(idx), i2) == pytest.approx(b2, abs=1e-9)
        # the angle between the two pendants is the specified one
        v1 = st.coords[i1] - st.coords[int(idx)]
        v2 = st.coords[i2] - st.coords[int(idx)]
        cos = v1 @ v2 / np.linalg.norm(v1) / np.linalg.norm(v2)
        assert np.degrees(np.arccos(cos)) == pytest.approx(spec.sub_angle, abs=1e-6)
        seen_asymmetric |= b1 != b2
    assert seen_asymmetric  # the test would be vacuous otherwise


@pytest.mark.parametrize("polymer", [CFE, CDFE])
def test_repeat_is_neutral_and_nothing_overlaps(polymer):
    st = build_chain(polymer, np.full(8, 180.0), cap=True)
    assert abs(float(st.charges.sum())) < 1e-12
    for spec in polymer.backbone:  # each backbone atom's own group is neutral too
        assert abs(spec.charge + sum(spec.sub_charges)) < 1e-12
    d = np.linalg.norm(st.coords[:, None] - st.coords[None], axis=-1) + np.eye(st.n_atoms) * 9
    assert d.min() > 0.9


@pytest.mark.parametrize("polymer", [CFE, CDFE])
@pytest.mark.parametrize("cap", [True, False])
def test_batched_builder_matches_the_single_one(polymer, cap):
    dih = np.random.default_rng(3).uniform(-180.0, 180.0, size=(5, 8))
    template, coords = build_chain_batch(polymer, dih, cap=cap)
    for m in range(5):
        single = build_chain(polymer, dih[m], cap=cap)
        assert single.elements == template.elements and single.bonds == template.bonds
        np.testing.assert_allclose(coords[m], single.coords, atol=1e-9)


@pytest.mark.parametrize("polymer", [CFE, CDFE])
def test_bond_angle_override_still_works(polymer):
    """The recently added per-repeat backbone-angle override must survive the change."""
    dih = np.random.default_rng(5).uniform(-180.0, 180.0, size=(3, 8))
    angles = np.array([118.0, 110.0])
    template, coords = build_chain_batch(polymer, dih, cap=True, bond_angles=angles)
    for m in range(3):
        single = build_chain(polymer, dih[m], cap=True, bond_angles=angles)
        np.testing.assert_allclose(coords[m], single.coords, atol=1e-9)
    # the substituents followed the opened-up backbone: bond lengths are still per pendant
    st = build_chain(polymer, dih[0], cap=True, bond_angles=angles)
    frozen = build_chain(polymer, dih[0], cap=True)
    assert not np.allclose(st.coords, frozen.coords)
    for k, idx in enumerate(st.backbone):
        spec = polymer.backbone[k % polymer.bonds_per_repeat]
        for i, b in zip(st.subs_of[int(idx)], spec.sub_bonds):
            assert distance(st.coords, int(idx), i) == pytest.approx(b, abs=1e-9)


@pytest.mark.parametrize("polymer", [CFE, CDFE])
def test_the_pack_coordinate_builder_sees_the_two_bond_lengths(polymer):
    """``pack``/``linegroup`` build coordinates themselves, forwarding ``spec.sub_bond``.

    They were not touched by phase 2, so this asserts that forwarding still gives the
    asymmetric geometry rather than silently placing both pendants at the first length.
    """
    tors = np.array([[180.0, 60.0]])
    ch = repeat_chains_from_torsions(polymer, "TG", tors)[0]
    bb = ch.backbone
    for local, idx in enumerate(bb):
        spec = polymer.backbone[local % polymer.bonds_per_repeat]
        # the two atoms following each backbone atom in the block are its pendants
        for off, b in zip((1, 2), spec.sub_bonds):
            assert np.linalg.norm(ch.coords[int(idx) + off] - ch.coords[int(idx)]) == pytest.approx(b, abs=1e-9)


# ------------------------------------------------------------------- 3. chirality
@pytest.mark.parametrize("polymer", [CFE, CDFE])
def test_the_chain_built_is_isotactic(polymer):
    """Every stereocentre gets the same configuration relative to the chain direction.

    The invariant is the sign of ``(s1 - C) . [(prev - C) x (next - C)]``: it is the
    local relative configuration, and it is conformation-independent.  It comes out the
    same at every backbone atom, which is what "isotactic" means; a syndiotactic chain
    would alternate.  Checked again the classical way on the planar zigzag, where the
    stereocentres all share a parity and so must put the same pendant on the same side
    of the backbone plane.
    """
    st = build_chain(polymer, np.random.default_rng(11).uniform(-180, 180, size=8), cap=False)
    bb = [int(i) for i in st.backbone]
    signs = []
    for k in range(1, len(bb) - 1):
        c = st.coords[bb[k]]
        s1 = st.coords[st.subs_of[bb[k]][0]]
        signs.append(np.sign((s1 - c) @ np.cross(st.coords[bb[k - 1]] - c, st.coords[bb[k + 1]] - c)))
    assert len(set(signs)) == 1 and signs[0] > 0

    flat = build_chain(polymer, np.full(8, 180.0), cap=False)
    fbb = [int(i) for i in flat.backbone]
    normal = np.cross(flat.coords[fbb[1]] - flat.coords[fbb[0]], flat.coords[fbb[2]] - flat.coords[fbb[0]])
    normal /= np.linalg.norm(normal)
    stereo = [k for k in range(1, len(fbb) - 1) if polymer.backbone[k % polymer.bonds_per_repeat].is_stereocentre]
    assert len(stereo) >= 2
    sides = [float((flat.coords[flat.subs_of[fbb[k]][0]] - flat.coords[fbb[k]]) @ normal) for k in stereo]
    assert all(abs(s - sides[0]) < 1e-9 for s in sides)  # same pendant, same side: isotactic


def test_a_user_is_warned_that_the_chain_has_a_tacticity():
    chain_mod._CHIRAL_WARNED.discard(CFE.name)
    with pytest.warns(UserWarning, match="ISOTACTIC"):
        build_chain(CFE, np.full(6, 180.0))
    with warnings.catch_warnings():  # once per polymer, not once per conformer
        warnings.simplefilter("error")
        build_chain(CFE, np.full(6, 180.0))


@pytest.mark.parametrize("polymer", [PVDF, CFE, CDFE])
def test_the_surviving_symmetry_is_mirror_composed_with_chain_reversal(polymer):
    """What replaces plain mirror symmetry for a chiral chain, and it is exact.

    Reflecting a conformer gives the *enantiomeric* chain in conformation ``-phi``.
    Reading that chain from the other end gives our chain back, because reversing the
    chain direction swaps ``prev`` and ``next`` and so flips the sign of the local
    frame that fixes each stereocentre's configuration.  Composing the two gives an
    identity that holds for chiral and achiral chains alike:

        E(phi_1 ... phi_N) == E(-phi_N ... -phi_1)

    For an achiral chain reversal is a symmetry on its own, so plain mirroring is
    exact too and ``fit_ris(symmetrize=True)`` is justified.  For a chiral chain only
    the composition survives -- which is why ``symmetrize`` would have to *transpose*
    the pair energies and swap the reversed bond types, not merely mirror the states.
    """
    dih = np.array([180.0, 60.0, 180.0, -60.0, 180.0, 180.0, 60.0, 180.0])  # not a palindrome
    ff = SimpleFF()
    e = ff.energy(build_chain(polymer, dih))
    assert abs(e - ff.energy(build_chain(polymer, -dih[::-1]))) < 1e-8 * max(1.0, abs(e))
    plain = abs(e - ff.energy(build_chain(polymer, -dih)))
    assert (plain > 1.0) if polymer.is_chiral else (plain < 1e-8)


@pytest.mark.parametrize("polymer", [CFE, CDFE])
def test_the_fitted_model_shows_the_same_thing(polymer):
    """The fitted RIS energies obey mirror+reversal and violate plain mirror.

    ``symmetrize=False`` therefore leaves a model that genuinely distinguishes G+ from
    G-, and the residual under the *correct* operation is only grid noise, which is the
    measurement that says the violation is physics and not a fitting artefact.
    ``adapt_angles=False`` so both models are read at the same state angles.
    ``angles="rigid"``: the relaxed scan CFE and CDFE now default to obeys the same
    relation (residual under 0.5 there too), but most of the rigid fit's 24 kcal/mol
    G+/G- asymmetry was frozen-angle strain -- relaxed, CFE's is 3.0 -- so the "gross
    violation" floor below is a property of the rigid scan this test was written on.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        rep = fit_ris(polymer, SimpleFF(), step=45.0, n_monomers=3, third_order=False,
                      symmetrize=False, adapt_angles=False, angles="rigid")
    mir = np.array(THREE_STATE.mirror)
    e1, e2 = rep.model.first_order, rep.model.second_order
    # mirror + reversal: for B = 2 reversal swaps the two bond types of the 1-D scan and
    # transposes the pair, so e1[b, s] = e1[1 - b, m(s)] and e2[b, s, s'] = e2[b, m(s'), m(s)]
    assert np.abs(e1 - e1[[1, 0]][:, mir]).max() < 0.5
    assert np.abs(e2 - np.transpose(e2[:, mir][:, :, mir], (0, 2, 1))).max() < 0.5
    # plain mirror, which is what symmetrize=True averages over, is violated grossly
    assert np.abs(e1 - e1[:, mir]).max() > 5.0
    assert np.abs(e2 - e2[:, mir][:, :, mir]).max() > 5.0


# --------------------------------------------------------------------- the funnel
@pytest.mark.parametrize("polymer", [CFE, CDFE])
def test_the_whole_funnel_runs(polymer):
    """Cheap fit, enumeration, a periodic chain and one packing energy: the harness works."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        rep = fit_ris(polymer, SimpleFF(), step=45.0, n_monomers=3, third_order=False, symmetrize=False)
    m = rep.model
    assert m.B == 2 and m.second_order.shape == (2, 3, 3)
    assert np.isfinite(m.first_order).all() and np.isfinite(m.second_order).all()

    cands = enumerate_periodic(polymer, m, max_period=4, k_per_period=10)
    assert cands and all(c.helix.rise_per_bond > 0.3 for c in cands)

    ch = periodic_chain(polymer, m.parse("TT"), m.states)
    assert len(ch.elements) == polymer.atoms_per_repeat * ch.n_monomers
    pk = CrystalPacker(ch, n_chains=2)
    e = pk.energy(np.array([[6.5, 10.0, 90.0, 30.0, 200.0, 1.0, 0]]))
    assert np.isfinite(e).all()


def test_symmetrize_auto_protects_chiral_polymers():
    """Mirror averaging is exact for an achiral chain and destructive for a chiral one.

    Reflecting a chiral chain gives its enantiomer, so G+ and G- genuinely differ;
    averaging them silently deletes that. The default must mirror PVDF and must not
    mirror CFE or CDFE.
    """
    import warnings as _w

    import numpy as np

    from polyfind.forcefield import SimpleFF, fit_ris
    from polyfind.polymers import get_polymer

    ff = SimpleFF()

    def mirror_residual(model):
        m = np.array(model.states.mirror)
        return float(np.abs(model.second_order - model.second_order[:, m][:, :, m]).max())

    # achiral: the default mirrors, and it was exact anyway
    pvdf = fit_ris(get_polymer("pvdf"), ff, step=30.0, n_monomers=4, third_order=False).model
    assert mirror_residual(pvdf) < 1e-9

    # chiral: the default must PRESERVE the asymmetry rather than average it away
    # (angles="rigid": this measures the symmetrisation, which does not depend on how the
    # backbone angles are treated, and the relaxed scan CFE and CDFE default to is slow)
    for name, floor in (("cfe", 1.0), ("cdfe", 1.0)):
        with _w.catch_warnings():
            _w.simplefilter("ignore")
            auto = fit_ris(get_polymer(name), ff, step=30.0, n_monomers=4, third_order=False, angles="rigid").model
            forced = fit_ris(get_polymer(name), ff, step=30.0, n_monomers=4, third_order=False, symmetrize=True,
                             angles="rigid").model
        assert mirror_residual(auto) > floor, f"{name}: real G+/G- asymmetry was averaged away"
        assert mirror_residual(forced) < 1e-9  # forcing it still works, but is wrong

    # and forcing it on a chiral polymer warns
    with pytest.warns(UserWarning, match="enantiomer"):
        fit_ris(get_polymer("cfe"), ff, step=60.0, n_monomers=3, third_order=False, symmetrize=True)


def test_chiral_fits_are_symmetrised_over_reflection_with_reversal():
    """The symmetry a chiral chain does have is reflection composed with reversal.

    Reversing the chain swaps each stereocentre's neighbours and undoes the
    reflection, so E(phi_1..phi_N) == E(-phi_N..-phi_1) holds for any linear repeat.
    The default uses it for chiral polymers, which restores noise averaging without
    destroying the real G+/G- difference.
    """
    import warnings as _w

    import numpy as np

    from polyfind.forcefield import SimpleFF, fit_ris
    from polyfind.polymers import get_polymer

    ff = SimpleFF()
    for name in ("cfe", "cdfe"):
        with _w.catch_warnings():
            _w.simplefilter("ignore")
            # rigid: the symmetrisation is what is measured, and the relaxed scan is slow
            m = fit_ris(get_polymer(name), ff, step=20.0, n_monomers=4, third_order=False, angles="rigid").model
        mir = np.array(m.states.mirror)
        e2 = m.second_order
        reversal = np.abs(e2 - np.transpose(e2[:, mir][:, :, mir], (0, 2, 1))).max()
        plain = np.abs(e2 - e2[:, mir][:, :, mir]).max()
        assert reversal < 1e-9, "reversal symmetry should be enforced exactly"
        assert plain > 1.0, "the real mirror asymmetry must survive"
