import numpy as np
import pytest

from polyfind.chain import (
    build_backbone,
    build_chain,
    build_chain_batch,
    distance,
    angle,
    dihedral,
    nerf,
    substituent_positions,
)
from polyfind.polymers import BackboneAtom, PVDF, PE


def test_backbone_geometry_reproduced():
    rng = np.random.default_rng(0)
    dih = rng.uniform(-180, 180, size=12)
    bb = build_backbone(PVDF, dih, xp=np)
    assert bb.shape == (15, 3)
    for k in range(14):
        # read the spec rather than repeating it: this asserts that the builder
        # reproduces whatever geometry the Polymer declares, which is the property
        # under test.  A literal here only re-asserted the value of the day, and broke
        # when PVDF's C-C was corrected from 1.54 to the DFT Form I 1.528 A.
        assert distance(bb, k, k + 1) == pytest.approx(PVDF.bond_length, abs=1e-9)
    for k in range(1, 14):
        expected = PVDF.backbone[k % 2].backbone_angle
        assert angle(bb, k - 1, k, k + 1) == pytest.approx(expected, abs=1e-6)
    for j in range(12):
        assert dihedral(bb, j, j + 1, j + 2, j + 3) == pytest.approx(dih[j], abs=1e-6)


def test_batched_build_equals_single():
    rng = np.random.default_rng(1)
    dih = rng.uniform(-180, 180, size=(5, 9))
    bb = build_backbone(PE, dih, xp=np)
    assert bb.shape == (5, 12, 3)
    for m in range(5):
        assert np.allclose(bb[m], build_backbone(PE, dih[m], xp=np))


def test_full_chain_substituents_and_caps():
    dih = np.full(6, 180.0)
    s = build_chain(PVDF, dih, cap=True)
    # 9 backbone atoms, 18 substituents, 2 caps
    assert s.n_atoms == 9 + 18 + 2
    assert sum(1 for e in s.elements if e == "F") == 8  # 4 CF2 groups: atoms 1,3,5,7
    for bb_idx in s.backbone:
        spec = PVDF.backbone[list(s.backbone).index(bb_idx) % 2]
        s1, s2 = s.subs_of[bb_idx]
        assert distance(s.coords, bb_idx, s1) == pytest.approx(spec.sub_bond, abs=1e-9)
        assert angle(s.coords, s1, bb_idx, s2) == pytest.approx(spec.sub_angle, abs=1e-6)
    # dihedrals unchanged by the virtual-atom extension
    for j in range(6):
        a, b, c, d = s.dihedral_atoms(j)
        assert dihedral(s.coords, a, b, c, d) == pytest.approx(180.0, abs=1e-6) or dihedral(s.coords, a, b, c, d) == pytest.approx(-180.0, abs=1e-6)
    # charges: neutral overall
    assert abs(s.charges.sum()) < 1e-9
    # no two atoms overlap
    d = np.linalg.norm(s.coords[:, None] - s.coords[None], axis=-1) + np.eye(s.n_atoms) * 10
    assert d.min() > 1.0


@pytest.mark.parametrize("polymer,n_dih", [(PE, 8), (PVDF, 8)])
@pytest.mark.parametrize("cap", [True, False])
def test_build_chain_batch_matches_build_chain(polymer, n_dih, cap):
    rng = np.random.default_rng(42)
    M = 7
    dih = rng.uniform(-180, 180, size=(M, n_dih))
    template, coords = build_chain_batch(polymer, dih, cap=cap)
    assert coords.shape == (M, template.n_atoms, 3)
    for m in range(M):
        single = build_chain(polymer, dih[m], cap=cap)
        assert single.n_atoms == template.n_atoms
        assert single.elements == template.elements
        assert single.bonds == template.bonds
        np.testing.assert_allclose(coords[m], single.coords, atol=1e-9)
    # row 0 is exactly the template used to build it
    np.testing.assert_allclose(coords[0], template.coords, atol=1e-9)


def test_substituent_positions_accepts_one_bond_length_or_two():
    """A pair of bond lengths puts each pendant at its own distance and keeps the angle.

    A single length must give exactly the old answer, bit for bit -- that is what lets
    the callers in ``pack``/``linegroup`` go on forwarding ``spec.sub_bond`` unchanged.
    """
    p, x, n = np.array([-1.2, 0.9, 0.0]), np.zeros(3), np.array([1.3, 0.8, 0.2])
    s1, s2 = substituent_positions(p, x, n, 1.35, 108.0)
    t1, t2 = substituent_positions(p, x, n, (1.35, 1.35), 108.0)
    assert s1.tobytes() == t1.tobytes() and s2.tobytes() == t2.tobytes()

    a1, a2 = substituent_positions(p, x, n, (1.35, 1.77), 108.0)
    assert np.linalg.norm(a1 - x) == pytest.approx(1.35, abs=1e-12)
    assert np.linalg.norm(a2 - x) == pytest.approx(1.77, abs=1e-12)
    # the inter-substituent angle is the specified one whatever the bond lengths...
    cos = (a1 - x) @ (a2 - x) / np.linalg.norm(a1 - x) / np.linalg.norm(a2 - x)
    assert np.degrees(np.arccos(cos)) == pytest.approx(108.0, abs=1e-9)
    # ...and the equal-angle split means only the length changed: same direction as before
    assert np.allclose(_unit_vec(a1 - x), _unit_vec(s1 - x), atol=1e-12)
    assert np.allclose(_unit_vec(a2 - x), _unit_vec(s2 - x), atol=1e-12)


def _unit_vec(v):
    return v / np.linalg.norm(v)


def test_substituent_positions_rejects_a_bad_pair():
    p, x, n = np.array([-1.2, 0.9, 0.0]), np.zeros(3), np.array([1.3, 0.8, 0.2])
    with pytest.raises(ValueError):
        substituent_positions(p, x, n, (1.0, 1.2, 1.4), 108.0)


def test_backbone_atom_scalar_and_pair_specs():
    sym = BackboneAtom("C", "Cl", 1.77, 114.0, 110.0, +0.20, -0.10)
    assert sym.substituents == ("Cl", "Cl")  # a two-letter symbol is one element, not two
    assert sym.sub_bonds == (1.77, 1.77) and sym.sub_charges == (-0.10, -0.10)
    assert not sym.is_stereocentre

    asym = BackboneAtom("C", ("F", "Cl"), (1.35, 1.77), 114.0, 108.0, +0.30, (-0.20, -0.10))
    assert asym.substituents == ("F", "Cl")
    assert asym.sub_bonds == (1.35, 1.77) and asym.sub_charges == (-0.20, -0.10)
    assert asym.is_stereocentre


def test_all_trans_pvdf_rise():
    """The closed form, and the experimental repeat it is supposed to reproduce.

    The rise per monomer of an all-trans chain is l (sin(a1/2) + sin(a2/2)), and for
    beta-PVDF it *is* the crystallographic c: 2.56 A (Hasegawa; docs/REFERENCES.md 1.1).
    Both halves are asserted, the second because it is the physical claim: with the DFT
    Form I C-C distance of 1.528 A the rise is 2.5630 A, +0.12% against experiment,
    where the textbook 1.54 A this used to carry gave 2.5831 A, +0.90%.
    """
    a1, a2 = (b.backbone_angle for b in PVDF.backbone)
    bb = build_backbone(PVDF, np.full(10, 180.0), xp=np)
    rise = np.linalg.norm(bb[12] - bb[2]) / 5
    expected = PVDF.bond_length * (np.sin(np.deg2rad(a1 / 2)) + np.sin(np.deg2rad(a2 / 2)))
    assert rise == pytest.approx(expected, abs=1e-6)
    assert rise == pytest.approx(2.56, abs=0.01)  # the experimental beta repeat
