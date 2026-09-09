import numpy as np
import pytest

from polyfind.chain import build_backbone, build_chain, build_chain_batch, distance, angle, dihedral, nerf
from polyfind.polymers import PVDF, PE


def test_backbone_geometry_reproduced():
    rng = np.random.default_rng(0)
    dih = rng.uniform(-180, 180, size=12)
    bb = build_backbone(PVDF, dih, xp=np)
    assert bb.shape == (15, 3)
    for k in range(14):
        assert distance(bb, k, k + 1) == pytest.approx(1.54, abs=1e-9)
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


def test_all_trans_pvdf_rise():
    # all-trans rise per monomer = l (sin(a1/2) + sin(a2/2))
    bb = build_backbone(PVDF, np.full(10, 180.0), xp=np)
    rise = np.linalg.norm(bb[12] - bb[2]) / 5
    expected = 1.54 * (np.sin(np.deg2rad(114 / 2)) + np.sin(np.deg2rad(114 / 2)))
    assert rise == pytest.approx(expected, abs=1e-6)
