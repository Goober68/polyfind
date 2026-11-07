"""polyfind.pendants: rigid tilts of the CH2/CF2 groups and the energy the fitted potential releases through them."""
import numpy as np
import pytest

from polyfind import pendants as PD
from polyfind.fitting import PRESET_PARAMETERS
from polyfind.mechanics import refined_reference
from polyfind.polymers import PVDF, belt_variant

PRESET = "pvdf-dft-valence-belt3"


def _reference(seq, label):
    polymer = belt_variant(PVDF, geometry="pbe")
    with PRESET_PARAMETERS[PRESET].applied():
        ref, rr = refined_reference(polymer, seq, label=label, n_chains=2)
    return ref, polymer


def test_a_group_turns_rigidly():
    """The displacement map is a proper rotation of each group about its carbon: bond lengths and the
    X-C-X angle are invariant, and a rock/wag/twist about the frame's own axes is read back by the
    same construction examples/dft_cell_geometry.py uses."""
    theta, a = np.radians(112.0), np.radians(108.0)
    c1 = 1.53 * np.array([-np.sin(theta / 2), -np.cos(theta / 2), 0.0])
    c2 = 1.53 * np.array([np.sin(theta / 2), -np.cos(theta / 2), 0.0])
    u, v = np.array([0.0, 1.0, 0.0]), np.array([1.0, 0.0, 0.0])
    w = np.cross(u, v)
    f = [1.38 * (np.cos(a / 2) * u + s * np.sin(a / 2) * w) for s in (1.0, -1.0)]
    P = np.array([[0.0, 0.0, 0.0], c1, c2, f[0], f[1]]) + 6.0
    lattice = np.diag([12.0, 12.0, 12.0])
    groups = PD.pendant_groups(lattice, P, ["C", "C", "C", "F", "F"])
    assert len(groups) == 1 and groups[0].substituents == "FF" and groups[0].carbon == 0
    g = groups[0]
    assert np.allclose(g.frame @ g.frame.T, np.eye(3))
    tilt = np.radians([5.0, -4.0, 3.0])
    D = PD.group_displacements(groups, tilt, len(P)).reshape(-1, 3)
    assert np.allclose(D[:3], 0.0)
    Q = P + D
    for j in g.x:
        assert abs(np.linalg.norm(Q[j] - Q[0]) - 1.38) < 1e-12
    d1, d2 = Q[g.x[0]] - Q[0], Q[g.x[1]] - Q[0]
    assert abs(np.degrees(np.arccos(d1 @ d2 / np.linalg.norm(d1) / np.linalg.norm(d2))) - 108.0) < 1e-9
    # rock alone: the bisector tilts toward a neighbour by exactly that angle
    D = PD.group_displacements(groups, np.radians([5.0, 0.0, 0.0]), len(P)).reshape(-1, 3)
    m = (P[g.x[0]] + D[g.x[0]] - P[0]) + (P[g.x[1]] + D[g.x[1]] - P[0])
    uu, vv, ww = g.frame
    assert abs(abs(np.degrees(np.arctan2(m @ vv, m @ uu))) - 5.0) < 1e-9 and abs(m @ ww) < 1e-12


def test_beta_releases_nothing_and_alpha_releases_the_documented_amount():
    """Beta's planar zigzag has a mirror through every group, so the symmetric placement is a
    stationary point and the relaxation leaves it; alpha's TGTG' does not, and under the belt3
    preset its groups turn (CF2 wag the largest) for -0.226 kcal/mol per monomer
    (docs/CRYSTAL_STAGE.md 5.22)."""
    ref_b, polymer = _reference([0, 0], "beta")
    with PRESET_PARAMETERS[PRESET].applied():
        tb = PD.tilt_relaxation(ref_b, polymer)
    assert tb.n_monomers == 2 and len(tb.groups) == 4
    assert abs(tb.dE_per_monomer) < 1e-6 and np.abs(tb.tilts_deg).max() < 1e-3
    ref_a, polymer = _reference([0, 1, 0, 2], "alpha")
    with PRESET_PARAMETERS[PRESET].applied():
        ta = PD.tilt_relaxation(ref_a, polymer)
    assert ta.n_monomers == 4 and len(ta.groups) == 8
    assert -0.26 < ta.dE_per_monomer < -0.19, ta.dE_per_monomer
    t = ta.tilts_by_substituents
    assert t["FF"]["wag"] > t["FF"]["rock"] > t["FF"]["twist"] > 0.3
    assert 2.0 < t["FF"]["wag"] < 3.5 and 2.0 < t["HH"]["wag"] < 3.2
    # the belt sites follow their CF2 groups: every site is still in its carbon's CF2 plane
    el = ta.elements
    assert (el == "X").sum() > 0
    J, atoms = PD.belt_site_constraint(ref_a.packer, ta.P1)  # raises if a site were off the plane
    assert J.shape == (3 * len(el), 3 * len(atoms))
    # the glide of TGTG' maps the groups onto one another: equal magnitudes within a substituent type
    for sub in ("FF", "HH"):
        mags = np.array([np.abs(x) for g, x in zip(ta.groups, ta.tilts_deg) if g.substituents == sub])
        assert np.ptp(mags, axis=0).max() < 0.05


def test_the_stretch_constant_floor_only_raises_the_soft_bonds():
    """The valence fit's C-C stretch constant sits at its bound (100 kcal/mol/A^2); a floor of 540
    raises it and leaves the stiffer C-H and C-F constants alone, and the pinned reference stays a
    stationary point (the equilibrium lengths are the reference's own)."""
    ref, polymer = _reference([0, 0], "beta")
    with PRESET_PARAMETERS[PRESET].applied():
        pk0, _, _ = PD.pinned_valence_packer(ref, polymer)
        pk1, _, _ = PD.pinned_valence_packer(ref, polymer, k_bond_min=540.0)
    k0, k1 = np.asarray(pk0._valence.bond_k), np.asarray(pk1._valence.bond_k)
    assert k0.min() < 150.0 and k1.min() >= 540.0
    assert np.all(k1 >= k0) and np.all(np.where(k0 >= 540.0, k1 == k0, True))
    np.testing.assert_allclose(pk1._valence.bond_r0, pk0._valence.bond_r0)
