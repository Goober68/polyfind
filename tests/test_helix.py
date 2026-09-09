import numpy as np
import pytest

from polyfind.helix import kabsch, screw_decompose, helix_parameters, canonical_sequence, sequence_images, is_primitive, rotation_to_z
from polyfind.polymers import PVDF, PE, THREE_STATE

T, GP, GM = 0, 1, 2


def test_kabsch_recovers_random_rigid_transform():
    rng = np.random.default_rng(0)
    P = rng.normal(size=(7, 3))
    q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    if np.linalg.det(q) < 0:
        q[:, 0] *= -1
    t = rng.normal(size=3)
    Q = P @ q.T + t
    R, tt = kabsch(P, Q)
    assert np.allclose(R, q) and np.allclose(tt, t)


def test_screw_decompose_known_screw():
    th = np.deg2rad(50.0)
    R = np.array([[np.cos(th), -np.sin(th), 0], [np.sin(th), np.cos(th), 0], [0, 0, 1.0]])
    p0 = np.array([1.0, 2.0, 0.0])
    d = 1.7
    t = p0 - R @ p0 + np.array([0, 0, d])  # screw about the line through p0 parallel to z
    n, theta, rise, p = screw_decompose(R, t)
    assert np.allclose(n, [0, 0, 1]) and theta == pytest.approx(50.0) and rise == pytest.approx(d)
    assert np.allclose(p[:2], p0[:2], atol=1e-9)


def test_pe_all_trans_is_21_helix():
    h = helix_parameters(PE, [T], THREE_STATE)
    assert abs(abs(h.rotation_per_period) - 180.0) < 1e-6
    assert h.periods_per_repeat == 2 and h.turns_per_repeat == 1
    assert h.c == pytest.approx(2 * 1.54 * np.sin(np.deg2rad(56)), abs=1e-6)  # 2.553 A
    assert h.label == "planar zigzag"


def test_pe_gauche_helix_handedness():
    hp = helix_parameters(PE, [GP], THREE_STATE)
    hm = helix_parameters(PE, [GM], THREE_STATE)
    assert hp.rotation_per_period == pytest.approx(-hm.rotation_per_period)
    assert hp.rise_per_period == pytest.approx(hm.rise_per_period)
    assert hp.rise_per_period < 1.27  # more compact than trans
    assert 60 < abs(hp.rotation_per_period) < 120


def test_pe_TG_is_31_helix():
    # alternating TG (the isotactic-polypropylene-type chain) is the classic 3/1 helix
    h = helix_parameters(PE, [T, GP], THREE_STATE)
    assert abs(abs(h.rotation_per_period) - 120.0) < 3.0
    assert h.periods_per_repeat == 3 and h.turns_per_repeat == 1


def test_pvdf_polymorph_chain_repeats():
    beta = helix_parameters(PVDF, [T, T], THREE_STATE)
    alpha = helix_parameters(PVDF, [T, GP, T, GM], THREE_STATE)
    gamma = helix_parameters(PVDF, [T, T, T, GP, T, T, T, GM], THREE_STATE)
    assert abs(beta.rotation_per_period) < 1e-6 and beta.c == pytest.approx(2.58, abs=0.02)
    assert abs(alpha.rotation_per_period) < 1e-6 and 4.3 < alpha.c < 4.9  # expt 4.62
    assert abs(gamma.rotation_per_period) < 1e-6 and 8.8 < gamma.c < 9.5  # expt 9.20
    assert beta.radius_all > beta.radius_backbone > 0
    assert alpha.radius_backbone > beta.radius_backbone


def test_pvdf_TG_is_a_true_helix():
    h = helix_parameters(PVDF, [T, GP], THREE_STATE)
    assert abs(h.rotation_per_period) > 10
    assert h.periods_per_repeat is None or h.periods_per_repeat >= 2


def test_canonical_form_symmetries():
    st = THREE_STATE
    a = canonical_sequence([T, GP, T, GM], st, 2)
    assert a == canonical_sequence([T, GM, T, GP], st, 2)  # shift by one monomer / mirror
    assert a == canonical_sequence([GM, T, GP, T], st, 2)  # reversal
    b = canonical_sequence([T, T, GP, GP], st, 2)
    assert b != a
    imgs = sequence_images([T, GP, T, T], st, 2)
    assert (GP, T, T, T) in imgs or (T, GP, T, T) in imgs
    assert len(imgs) <= 8


def test_is_primitive():
    assert is_primitive([T, GP, T, GM], 2)
    assert not is_primitive([T, T, T, T], 2)
    assert not is_primitive([T, GP, T, GP], 2)
    assert is_primitive([T, GP], 2)


def test_rotation_to_z():
    rng = np.random.default_rng(3)
    for _ in range(5):
        n = rng.normal(size=3)
        n /= np.linalg.norm(n)
        R = rotation_to_z(n)
        assert np.allclose(R @ n, [0, 0, 1]) and np.allclose(np.linalg.det(R), 1)
