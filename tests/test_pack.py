import numpy as np
import pytest

from polyfind.pack import periodic_chain, periodic_chain_from_torsions, CrystalPacker, pack, polish, to_cif, default_bounds
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


@pytest.mark.parametrize("poly,seq", [(PE, [T]), (PVDF, [T, T]), (PVDF, [T, GP, T, GM]), (PVDF, [T, T, T, GP, T, T, T, GM])])
def test_flipping_a_chain_beyond_the_cutoff_costs_nothing(poly, seq):
    """A flip is an isometry of the chain, so at infinite separation it is free.

    The chain-2 coordinates are mirrored to (x, -y, -z), which maps its z-image k onto
    image -k; scoring it with the unmirrored exclusion matrices left its bonded pairs
    unexcluded and put thousands of kcal/mol on every antiparallel cell, so the
    antiparallel half of the search space was unreachable.
    """
    ch = periodic_chain(poly, seq, THREE_STATE)
    pk = CrystalPacker(ch, n_chains=2)
    far = 60.0  # >> cutoff + 2 * chain radius: the two chains cannot interact
    para = np.array([[far, far, 90.0, 37.0, 111.0, 0.7, 0.0]])
    anti = para.copy()
    anti[0, 6] = 1.0
    e_para, e_anti = float(pk.energy(para)[0]), float(pk.energy(anti)[0])
    assert e_anti == pytest.approx(e_para, abs=1e-9)
    # and both are exactly the two isolated chains: intra + torsion, no interaction
    assert e_para == pytest.approx(pk.e_intra + pk.e_torsion, abs=1e-9)
    single = CrystalPacker(ch, n_chains=1)
    assert e_para == pytest.approx(2.0 * float(single.energy(np.array([[far, far, 90.0, 0.0, 0.0, 0.0, 0.0]]))[0]), abs=1e-9)


def test_antiparallel_packing_is_competitive_at_contact():
    """At a realistic cell the two orientations differ by ordinary packing energies.

    The exclusion bug put ~3200 (PE) to ~6290 (PVDF) kcal/mol per cell on flip=1, so
    antiparallel could never win; the honest difference is a few kcal/mol per monomer,
    and for alpha-PVDF antiparallel is the lower of the two once relaxed.
    """
    keys = ["a", "b", "gamma", "phi1", "phi2", "dz"]
    for poly, seq, cell in ((PE, [T], [4.95, 7.42]), (PVDF, [T, GP, T, GM], [4.96, 9.64])):
        ch = periodic_chain(poly, seq, THREE_STATE)
        pk = CrystalPacker(ch)
        bd = default_bounds(ch)
        lo = np.array([bd[k][0] for k in keys])
        hi = np.array([bd[k][1] for k in keys])
        free = [i for i in range(6) if hi[i] > lo[i]]
        e = {}
        for flip in (0.0, 1.0):
            p = np.array([cell[0], cell[1], 90.0, 45.0, 135.0, 0.5 * ch.c, flip])
            assert abs(float(pk.energy(p[None])[0])) < 500.0  # not the bug's ~3000-6300
            x = polish(pk, p, lo, hi, free)
            e[flip] = float(pk.energy(x[None])[0]) / (pk.n_chains * ch.n_monomers)
        assert abs(e[1.0] - e[0.0]) < 5.0  # kcal/mol per monomer


def test_polish_lbfgs_reaches_the_nelder_mead_minimum_with_fewer_calls():
    ch = periodic_chain(PE, [T], THREE_STATE)
    packer = CrystalPacker(ch, n_chains=2)
    b = default_bounds(ch)
    keys = ["a", "b", "gamma", "phi1", "phi2", "dz"]
    lo = np.array([b[k][0] for k in keys])
    hi = np.array([b[k][1] for k in keys])
    free = [i for i in range(6) if hi[i] > lo[i]]
    rng = np.random.default_rng(0)
    cont = lo + (hi - lo) * rng.random((300, 6))
    params = np.concatenate([cont, np.zeros((300, 1))], axis=1)
    starts = params[np.argsort(packer.energy(params))[:2]]
    for start in starts:
        packer.n_energy_calls = 0
        x_nm = polish(packer, start, lo, hi, free, 1500, method="nelder-mead")
        calls_nm = packer.n_energy_calls
        packer.n_energy_calls = 0
        x_lb = polish(packer, start, lo, hi, free, 1500, method="lbfgs")
        calls_lb = packer.n_energy_calls
        e_nm = float(packer.energy(x_nm[None])[0])
        e_lb = float(packer.energy(x_lb[None])[0])
        assert e_lb <= e_nm + 1e-3
        assert calls_lb * 5 <= calls_nm  # at least a 5x reduction in kernel calls
        assert np.abs(x_lb[:2] - x_nm[:2]).max() < 0.02  # a, b
        for i in (3, 4):  # setting angles, modulo 360
            assert abs(((x_lb[i] - x_nm[i] + 180.0) % 360.0) - 180.0) < 0.5
        assert lo[0] <= x_lb[0] <= hi[0] and 0.0 <= x_lb[3] < 360.0 and 0.0 <= x_lb[5] < ch.c
    with pytest.raises(ValueError):
        polish(packer, starts[0], lo, hi, free, method="powell")


def test_pack_accepts_either_polish_method():
    ch = periodic_chain(PE, [T], THREE_STATE)
    kw = dict(n_chains=2, n_random=300, n_refine=2, maxfev=400, screen="random")
    a = pack(ch, rng=np.random.default_rng(0), method="lbfgs", **kw)
    b = pack(ch, rng=np.random.default_rng(0), method="nelder-mead", **kw)
    assert a[0].energy_per_cell <= b[0].energy_per_cell + 1e-3


def test_default_bounds_cover_known_cells():
    b = default_bounds(periodic_chain(PVDF, [T, GP, T, GM], THREE_STATE))
    assert b["a"][0] < 4.96 < b["a"][1] and b["b"][0] < 9.64 < b["b"][1]


# ------------------------------------------------------------------ the table screen
@pytest.mark.parametrize("poly,seq", [(PE, [T]), (PVDF, [T, T])])
def test_table_screen_is_no_worse_than_the_random_screen(poly, seq):
    """The exhaustive screen must not lose minima the 6,000-cell random sample finds.

    Both paths polish their starts with the exact kernel, so the comparison is between
    final, exact energies; the table only decides which cells are polished.
    """
    ch = periodic_chain(poly, seq, THREE_STATE)
    rnd = pack(ch, screen="random", n_random=1200, n_refine=4, rng=np.random.default_rng(0), maxfev=600)
    tab = pack(ch, screen="table", n_refine=4, maxfev=600, table_cache_dir=None)
    assert tab[0].energy_per_monomer <= rnd[0].energy_per_monomer + 1e-3
    # the table path returns the same kind of answer: sorted, deduplicated, exact
    assert [r.energy_per_cell for r in tab] == sorted(r.energy_per_cell for r in tab)
    assert all(abs(u.energy_per_cell - v.energy_per_cell) > 1e-3 for i, u in enumerate(tab) for v in tab[i + 1:])
    pk = CrystalPacker(ch, n_chains=2)
    assert pk.energy(tab[0].params[None])[0] == pytest.approx(tab[0].energy_per_cell, abs=1e-9)


def test_pair_table_cache_reuses_the_table_in_memory_and_on_disk(tmp_path):
    from polyfind import lattice_table as lt

    ch = periodic_chain(PE, [T], THREE_STATE)
    kw = dict(n_angle=24, n_z=4, dr=0.5)
    lt.clear_pair_table_cache()
    n0 = lt.pair_table_cache_info()["builds"]
    first = lt.pair_table(ch, cache_dir=str(tmp_path), **kw)
    assert lt.pair_table_cache_info()["builds"] == n0 + 1
    second = lt.pair_table(ch, cache_dir=str(tmp_path), **kw)
    assert second is first  # identical table, not rebuilt
    assert lt.pair_table_cache_info()["builds"] == n0 + 1
    # a fresh process (here: an empty in-process cache) picks the table up from disk
    files = list(tmp_path.glob("*.npz"))
    assert len(files) == 1 and lt.table_key(ch, **kw) in files[0].name
    lt.clear_pair_table_cache()
    third = lt.pair_table(ch, cache_dir=str(tmp_path), **kw)
    assert lt.pair_table_cache_info()["builds"] == n0 + 1
    assert np.array_equal(third.W, first.W) and third.W.dtype == first.W.dtype
    assert np.array_equal(third.e_intra, first.e_intra)
    for attr in ("r_min", "dr", "n_angle", "n_z", "c", "rc", "cap"):
        assert getattr(third, attr) == getattr(first, attr)
    # the key covers the potential and grid, not the cell: a different cutoff is a miss
    assert lt.table_key(ch, **kw) != lt.table_key(ch, cutoff=7.0, **kw)
    assert lt.table_key(ch, **kw) == lt.table_key(ch, symmetry=False, n_threads=1, **kw)


def test_pack_api_is_unchanged():
    """Every pre-existing call signature keeps working and keeps its old behaviour."""
    ch = periodic_chain(PE, [T], THREE_STATE)
    kw = dict(n_chains=2, n_random=200, n_refine=2, maxfev=300)
    a = pack(ch, screen="random", rng=np.random.default_rng(0), **kw)
    b = pack(ch, screen="random", rng=np.random.default_rng(0), **kw)
    assert [r.energy_per_cell for r in a] == [r.energy_per_cell for r in b]  # deterministic
    assert a[0].chain == "T" and a[0].n_chains == 2 and a[0].cell_coords.shape == (12, 3)
    # the positional form that predates the screen argument (screen is keyword-only in effect)
    assert pack(ch, 2, 200, 1, None, False, (0, 1), np.random.default_rng(1), 8.0, 1.0, 300, False, "lbfgs",
                screen="random")
    # a one-chain cell has no chain-pair table; it falls back to the random screen
    one = pack(ch, n_chains=1, n_random=200, n_refine=1, maxfev=300, rng=np.random.default_rng(0))
    assert one[0].n_chains == 1 and one[0].phi2 == 0.0 and one[0].dz == 0.0
    with pytest.raises(ValueError):
        pack(ch, screen="exhaustive", **kw)
