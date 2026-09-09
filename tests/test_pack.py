import numpy as np
import pytest

from polyfind.pack import periodic_chain, periodic_chain_from_torsions, CrystalPacker, pack, to_cif, default_bounds
from polyfind.polymers import PE, PVDF, THREE_STATE

T, GP, GM = 0, 1, 2


@pytest.mark.parametrize("poly,seq,c_expt", [(PE, [T], 2.55), (PVDF, [T, T], 2.56), (PVDF, [T, GP, T, GM], 4.62), (PVDF, [T, T, T, GP, T, T, T, GM], 9.20)])
def test_periodic_chain_repeat_lengths(poly, seq, c_expt):
    ch = periodic_chain(poly, seq, THREE_STATE)
    assert abs(ch.c - c_expt) / c_expt < 0.03
    assert ch.coords[:, 2].min() == pytest.approx(0.0)
    assert ch.coords[:, 2].max() < ch.c + 2.0  # pendant atoms may overhang the repeat
    # no atom pair within a block is unphysically close
    d = np.linalg.norm(ch.coords[:, None] - ch.coords[None], axis=-1) + np.eye(ch.n_atoms) * 9
    assert d.min() > 1.0


def test_periodic_chain_from_torsions_matches_ideal():
    ch_ideal = periodic_chain(PVDF, [T, GP, T, GM], THREE_STATE)
    ch_t = periodic_chain_from_torsions(PVDF, "TG+TG-", ch_ideal.dihedrals)
    assert ch_t.rotation_error < 1e-6
    assert ch_t.c == pytest.approx(ch_ideal.c, abs=1e-6)
    pk1, pk2 = CrystalPacker(ch_ideal), CrystalPacker(ch_t)
    p = np.array([[5.0, 9.6, 90.0, 20.0, 200.0, 1.0, 1]])
    assert pk1.energy(p)[0] == pytest.approx(pk2.energy(p)[0], abs=1e-6)
    # a deflected all-trans PVDF chain is not commensurate at ideal angles
    ch_d = periodic_chain_from_torsions(PVDF, "TT", np.array([170.0, 170.0]))
    assert ch_d.rotation_error > 1.0
    ch_d2 = periodic_chain_from_torsions(PVDF, "TT", np.array([180.0, 180.0]))
    assert ch_d2.rotation_error < 1e-6


def test_kernel_batch_equals_single_and_is_invariant():
    ch = periodic_chain(PE, [T], THREE_STATE)
    pk = CrystalPacker(ch, n_chains=2)
    rng = np.random.default_rng(0)
    M = 12
    params = np.column_stack([rng.uniform(4.5, 9, M), rng.uniform(4.5, 9, M), np.full(M, 90.0), rng.uniform(0, 360, M), rng.uniform(0, 360, M), rng.uniform(0, ch.c, M), rng.integers(0, 2, M)])
    E = pk.energy(params)
    for i in range(M):
        assert pk.energy(params[i : i + 1])[0] == pytest.approx(E[i], rel=1e-9)
    # rotating both chains by the same angle and swapping a<->b (with 90 deg rotation) leaves E unchanged
    p = params[0].copy()
    q = p.copy()
    q[3] += 37.0
    q[4] += 37.0
    assert pk.energy(p[None])[0] != pytest.approx(pk.energy(q[None])[0], abs=1e-6)  # setting angle matters
    r = p.copy()
    r[0], r[1] = p[1], p[0]
    r[3], r[4] = p[3] + 90, p[4] + 90
    assert pk.energy(r[None])[0] == pytest.approx(pk.energy(p[None])[0], rel=1e-6)
    # z-translation of chain 2 by a full repeat is a symmetry
    s = p.copy()
    s[5] = (p[5] + ch.c)
    assert pk.energy(s[None])[0] == pytest.approx(pk.energy(p[None])[0], rel=1e-6)


def test_pe_packing_recovers_orthorhombic_cell():
    ch = periodic_chain(PE, [T], THREE_STATE)
    res = pack(ch, n_chains=2, n_random=600, n_refine=3, rng=np.random.default_rng(0), maxfev=400)
    best = res[0]
    dims = sorted([best.a, best.b])
    # experimental 4.95 x 7.42 A; the illustrative potential should land within ~10%
    assert abs(dims[0] - 4.95) / 4.95 < 0.12 and abs(dims[1] - 7.42) / 7.42 < 0.12
    assert 0.85 < best.density < 1.15
    cif = to_cif(best)
    assert "_cell_length_a" in cif and cif.count("\nC") == 4


def test_default_bounds_cover_known_cells():
    b = default_bounds(periodic_chain(PVDF, [T, GP, T, GM], THREE_STATE))
    assert b["a"][0] < 4.96 < b["a"][1] and b["b"][0] < 9.64 < b["b"][1]
