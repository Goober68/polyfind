import numpy as np
import pytest

from polyfind.chain import angle as measure_angle
from polyfind.chain import build_backbone, build_chain, build_chain_batch, distance
from polyfind.linegroup import (
    LineGroupError,
    line_group,
    repeat_chains,
    repeat_rotvec,
    state_sequence,
    torsion_pattern,
    wrap180,
)
from polyfind.pack import periodic_chain, repeat_chains_from_torsions
from polyfind.polymers import PE, PVDF, THREE_STATE

T, GP, GM = 0, 1, 2

# The three PVDF polymorph sequences plus a screw and a plain zigzag.
SEQUENCES = {
    "TT": [T, T],
    "TG+TG-": [T, GP, T, GM],
    "TTTG+TTTG-": [T, T, T, GP, T, T, T, GM],
    "TG+": [T, GP],
}


# ------------------------------------------------------------------ bond angles (chain.py)
def test_bond_angle_override_is_used_and_defaults_to_the_polymer():
    dih = np.full(9, 180.0)
    plain = build_backbone(PVDF, dih, xp=np)
    same = build_backbone(PVDF, dih, xp=np, bond_angles=[114.0, 114.0])
    np.testing.assert_allclose(same, plain, atol=1e-12)
    bent = build_backbone(PVDF, dih, xp=np, bond_angles=[110.0, 119.0])
    for k in range(1, 11):
        assert measure_angle(bent, k - 1, k, k + 1) == pytest.approx([110.0, 119.0][k % 2], abs=1e-6)
        assert distance(bent, k, k + 1) == pytest.approx(PVDF.bond_length, abs=1e-12)


@pytest.mark.parametrize("cap", [True, False])
def test_build_chain_with_angles_matches_the_batched_builder_and_moves_substituents(cap):
    rng = np.random.default_rng(7)
    dih = rng.uniform(-180, 180, size=(4, 8))
    ba = np.array([[110.0, 119.0], [114.0, 114.0], [117.0, 108.0], [112.0, 116.0]])
    template, coords = build_chain_batch(PVDF, dih, cap=cap, bond_angles=ba)
    for m in range(4):
        single = build_chain(PVDF, dih[m], cap=cap, bond_angles=ba[m])
        np.testing.assert_allclose(coords[m], single.coords, atol=1e-9)
    # substituents follow the changed backbone: bond lengths and H-C-H angles are kept,
    # but the substituent positions differ from the frozen-angle chain
    frozen = build_chain(PVDF, dih[0], cap=cap)
    moved = build_chain(PVDF, dih[0], cap=cap, bond_angles=ba[0])
    assert np.abs(frozen.coords - moved.coords).max() > 0.05
    for k, idx in enumerate(moved.backbone):
        spec = PVDF.backbone[k % 2]
        s1, s2 = moved.subs_of[int(idx)]
        assert distance(moved.coords, int(idx), s1) == pytest.approx(spec.sub_bond, abs=1e-9)
        assert measure_angle(moved.coords, s1, int(idx), s2) == pytest.approx(spec.sub_angle, abs=1e-6)


def test_one_angle_set_shared_by_a_whole_batch():
    rng = np.random.default_rng(11)
    dih = rng.uniform(-180, 180, size=(3, 7))
    _, coords = build_chain_batch(PVDF, dih, cap=True, bond_angles=[111.0, 117.0])
    for m in range(3):
        single = build_chain(PVDF, dih[m], cap=True, bond_angles=[111.0, 117.0])
        np.testing.assert_allclose(coords[m], single.coords, atol=1e-12)


def test_polymer_definition_is_not_mutated():
    build_chain(PVDF, np.full(6, 180.0), bond_angles=[100.0, 130.0])
    assert [b.backbone_angle for b in PVDF.backbone] == [114.0, 114.0]


# ------------------------------------------------------------------ the pattern
def test_pattern_is_read_off_the_state_sequence():
    p = torsion_pattern(SEQUENCES["TT"], THREE_STATE, 2)
    assert (p.kind, p.period, p.n_params) == ("glide", 1, 1)
    assert p.signs == (1.0, -1.0)  # (180 - d, 180 + d)

    p = torsion_pattern(SEQUENCES["TG+TG-"], THREE_STATE, 2)
    assert (p.kind, p.period, p.n_params) == ("glide", 2, 2)
    assert p.group == (0, 1, 0, 1) and p.signs == (1.0, 1.0, -1.0, -1.0)  # (t, g, -t, -g)

    p = torsion_pattern(SEQUENCES["TTTG+TTTG-"], THREE_STATE, 2)
    assert (p.kind, p.period, p.n_params) == ("glide", 4, 4)
    assert p.signs == (1.0,) * 4 + (-1.0,) * 4

    p = torsion_pattern([T, GP] * 3, THREE_STATE, 2)  # a 3/1 screw
    assert (p.kind, p.period, p.n_params) == ("screw", 2, 2)
    assert p.signs == (1.0,) * 6 and p.group == (0, 1, 0, 1, 0, 1)

    p = torsion_pattern([T, T, GP, T, GM, T], THREE_STATE, 2)  # no symmetry
    assert p.kind == "free"


def test_pattern_reproduces_the_ideal_torsions():
    for name, seq in SEQUENCES.items():
        chain = periodic_chain(PVDF, seq, THREE_STATE)
        pattern = torsion_pattern(state_sequence(chain.name, THREE_STATE, len(chain.dihedrals)), THREE_STATE, 2)
        params, err = pattern.params_from_torsions(chain.dihedrals)
        assert err < 1e-9, name
        np.testing.assert_allclose(wrap180(pattern.torsions(params) - chain.dihedrals), 0.0, atol=1e-9)


def test_sequence_name_parsing():
    assert state_sequence("TG+TG-", THREE_STATE) == [T, GP, T, GM]
    assert state_sequence("TT", THREE_STATE, length=4) == [T, T, T, T]
    with pytest.raises(LineGroupError):
        state_sequence("TX", THREE_STATE)
    with pytest.raises(LineGroupError):
        state_sequence("TG+TG-", THREE_STATE, length=6)


# ------------------------------------------------------------------ the point of it all
@pytest.mark.parametrize("name", list(SEQUENCES))
def test_generated_chains_are_periodic_for_any_parameters(name):
    """The whole point: any free parameters in a sensible range give a periodic chain.

    The range differs by family: a glide's closure has a solution for wide torsion
    excursions, while a screw's closure *is* the commensurability of the helix, whose
    solutions run out once the free torsion has moved ~10 deg (the dependent torsion is
    by then 40 deg away from trans).  Outside that set the closure solve leaves its
    residual behind and the refinement's penalty term takes over.
    """
    chain = periodic_chain(PVDF, SEQUENCES[name], THREE_STATE)
    lg = line_group(PVDF, chain.name, chain.dihedrals)
    assert lg.n_params < len(chain.dihedrals) + PVDF.bonds_per_repeat  # a genuine reduction
    rng = np.random.default_rng(0)
    lim = np.where(lg.is_angle(), 8.0, 25.0 if lg.kind == "glide" else 10.0)
    X = np.vstack([np.zeros(lg.n_params), rng.uniform(-lim, lim, (40, lg.n_params))])
    tors, angles = lg.expand(X)
    for i in range(len(X)):
        # the exact residual: the rotation vector of the repeat transform
        assert np.linalg.norm(repeat_rotvec(PVDF, tors[i : i + 1], angles[i])[0]) < 1e-6
    # and as pack measures it, on the assembled chains (arccos near 1 floors at ~2e-6 deg)
    for ch in lg.chains(X[:8]):
        assert ch.rotation_error < 1e-5
        assert ch.c > 1.0


def test_zero_parameters_reproduce_the_ideal_chain():
    for name, seq in SEQUENCES.items():
        chain = periodic_chain(PVDF, seq, THREE_STATE)
        lg = line_group(PVDF, chain.name, chain.dihedrals)
        tors, angles = lg.expand(lg.x0)
        np.testing.assert_allclose(wrap180(tors[0] - chain.dihedrals), 0.0, atol=1e-8)
        np.testing.assert_allclose(angles[0], [114.0, 114.0], atol=1e-8)
        assert lg.chains(lg.x0)[0].c == pytest.approx(chain.c, abs=1e-6)


def test_all_trans_pvdf_closure_forces_equal_angles():
    """A one-monomer repeat cannot be straight with unequal backbone angles.

    The repeat transform of a planar zigzag is a rotation by theta1 - theta2 per
    monomer, so the closure condition removes both the deflection and the angle
    difference: only the mean angle is left free.
    """
    chain = periodic_chain(PVDF, [T, T], THREE_STATE)
    lg = line_group(PVDF, chain.name, chain.dihedrals)
    assert lg.n_params == 1
    tors, angles = lg.expand(np.array([[-5.0], [0.0], [6.0]]))
    np.testing.assert_allclose(tors, 180.0, atol=1e-7)
    np.testing.assert_allclose(angles[:, 0], angles[:, 1], atol=1e-7)
    # ... and the rotation error of an unequal-angle all-trans chain is exactly the
    # angle difference, which is why nothing else can work
    err = np.linalg.norm(repeat_rotvec(PVDF, np.array([[180.0, 180.0]]), [111.0, 117.0])[0])
    assert err == pytest.approx(6.0, abs=1e-6)


def test_no_pattern_raises_so_the_caller_can_fall_back():
    with pytest.raises(LineGroupError):  # TTG+TG-T has neither a glide nor a shift
        line_group(PVDF, "TTG+TG-T", [180.0, 180.0, 60.0, 180.0, -60.0, 180.0])
    with pytest.raises(LineGroupError):  # torsions that are not of the sequence's pattern
        line_group(PVDF, "TG+TG-", [180.0, 60.0, 150.0, -60.0])


def test_polyethylene_all_trans():
    chain = periodic_chain(PE, [T], THREE_STATE)
    lg = line_group(PE, chain.name, chain.dihedrals)
    assert lg.kind == "glide" and lg.pattern.period == 1
    tors, angles = lg.expand(np.array([[4.0]]))
    assert np.linalg.norm(repeat_rotvec(PE, tors, angles[0])[0]) < 1e-6


# ------------------------------------------------------------------ the chain builder
def test_repeat_chains_matches_pack_when_the_angles_are_the_polymer_s_own():
    ideal = periodic_chain(PVDF, [T, GP, T, GM], THREE_STATE)
    rng = np.random.default_rng(3)
    tors = ideal.dihedrals + rng.uniform(-12, 12, (5, 4))
    mine = repeat_chains(PVDF, "TG+TG-", tors, None, align_to=ideal)
    theirs = repeat_chains_from_torsions(PVDF, "TG+TG-", tors, align_to=ideal)
    explicit = repeat_chains(PVDF, "TG+TG-", tors, [114.0, 114.0], align_to=ideal)
    for a, b, c in zip(mine, theirs, explicit):
        np.testing.assert_allclose(a.coords, b.coords, atol=1e-12)
        np.testing.assert_allclose(c.coords, b.coords, atol=1e-12)
        assert a.c == pytest.approx(b.c, abs=1e-12)
        assert a.rotation_error == pytest.approx(b.rotation_error, abs=1e-12)


def test_repeat_chains_carries_per_row_angles():
    ideal = periodic_chain(PVDF, [T, GP, T, GM], THREE_STATE)
    tors = np.repeat(ideal.dihedrals[None], 2, axis=0)
    chains = repeat_chains(PVDF, "TG+TG-", tors, np.array([[112.0, 118.0], [114.0, 114.0]]))
    assert chains[0].c != pytest.approx(chains[1].c, abs=1e-6)
    assert chains[1].c == pytest.approx(ideal.c, abs=1e-6)
