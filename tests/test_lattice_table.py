import numpy as np
import pytest

from polyfind.lattice_table import (
    PairTable,
    chain_pair_energy,
    direct_lattice_energy,
    fft_screen,
    intra_constants,
    intra_energy,
    lattice_sites,
    screen,
    table_energy,
)
from polyfind.pack import CrystalPacker, default_bounds, periodic_chain
from polyfind.polymers import PE, PVDF, THREE_STATE

T, GP, GM = 0, 1, 2
PE_T = (PE, [T])
PE_TG = (PE, [T, GP])  # a 3/1 helix: chiral, so the flipped self term is not symmetric
PVDF_TT = (PVDF, [T, T])
PVDF_A = (PVDF, [T, GP, T, GM])
PVDF_G = (PVDF, [T, T, T, GP, T, T, T, GM])


_TABLES: dict = {}


def _chain(case):
    poly, seq = case
    return periodic_chain(poly, seq, THREE_STATE)


def _table(case, **kw):
    """Tables are expensive; build each configuration once per session."""
    key = (case[0].name, tuple(case[1]), tuple(sorted(kw.items())))
    if key not in _TABLES:
        _TABLES[key] = PairTable.build(_chain(case), **kw)
    return _TABLES[key]


def _random_cells(chain, n, rng, gammas=(90.0,)):
    b = default_bounds(chain)
    g = rng.choice(np.asarray(gammas, dtype=float), n)
    return np.column_stack([
        rng.uniform(*b["a"], n), rng.uniform(*b["b"], n), g,
        rng.uniform(0, 360, n), rng.uniform(0, 360, n), rng.uniform(0, chain.c, n),
        rng.integers(0, 2, n).astype(float),
    ])


# ------------------------------------------------------------------ (a) decomposition
@pytest.mark.parametrize("case", [PE_T, PE_TG, PVDF_A, PVDF_G])
def test_lattice_decomposition_is_exact(case):
    """E(cell) = E_intra + 1/2 sum_self W + sum_centred W, with W evaluated directly."""
    chain = _chain(case)
    pk2 = CrystalPacker(chain, n_chains=2)
    pk1 = CrystalPacker(chain, n_chains=1)
    rng = np.random.default_rng(4)
    params = _random_cells(chain, 6, rng, gammas=(90.0, 78.0, 104.0))
    params[:2, 6] = [0.0, 1.0]  # make sure both flips appear
    assert np.abs(direct_lattice_energy(pk1, params) - pk2.energy(params)).max() < 1e-6
    # the parallel cell's constant is twice one chain's intra-chain energy
    assert intra_constants(pk2)[0] == pytest.approx(2 * intra_energy(pk1), rel=1e-9)


def test_flipped_self_term_sign_matters():
    """The flipped chain's self term uses the *negated* setting angle.  For a mirror-
    symmetric chain (PE all-trans, alpha, gamma) both signs agree; for a chiral one
    (the TG+ 3/1 helix) the naive sign is wrong, so the rule has teeth."""
    chain = _chain(PE_TG)
    pk1 = CrystalPacker(chain, n_chains=1)
    pk2 = CrystalPacker(chain, n_chains=2)
    r0 = chain.radius + 1.5
    p = np.array([1.7 * r0, 2.6 * r0, 90.0, 37.0, 211.0, 1.3, 1.0])
    a, b, gam, phi1, phi2, dz, _ = p
    rs, ts = lattice_sites(a, b, gam, pk1.rc + 2 * chain.radius)
    rho = np.stack([rs, np.zeros_like(rs)], axis=1)
    e_intra = intra_constants(pk2)
    wrong = float(e_intra[1])
    for phi, sgn in ((phi1, 1.0), (phi2, +1.0)):  # +1 for chain 2 is the wrong sign
        al = sgn * (phi - ts)
        wrong += 0.5 * chain_pair_energy(pk1, rho, 0.0, al, al, 0.0).sum()
    rc_, tc = lattice_sites(a, b, gam, pk1.rc + 2 * chain.radius, centred=True)
    wrong += chain_pair_energy(pk1, np.stack([rc_, np.zeros_like(rc_)], 1), dz, phi1 - tc, phi2 - tc, 1.0).sum()
    ref = pk2.energy(p[None])[0]
    assert abs(direct_lattice_energy(pk1, p[None])[0] - ref) < 1e-6
    assert abs(wrong - ref) > 1e-3


@pytest.mark.parametrize("case", [PE_T, PVDF_A])
def test_chain_pair_exchange_symmetry(case):
    """Swapping the two chains: rho -> -rho, alpha1 <-> alpha2, dz -> -dz (flip 0) and
    the in-plane mirror for flip 1.  These are what the builder uses to halve its work."""
    chain = _chain(case)
    pk1 = CrystalPacker(chain, n_chains=1)
    rng = np.random.default_rng(5)
    n = 8
    r = rng.uniform(4.0, 9.0, n)
    a1, a2 = rng.uniform(0, 360, n), rng.uniform(0, 360, n)
    dz = rng.uniform(0, chain.c, n)
    rho = np.stack([r, np.zeros_like(r)], axis=1)
    w0 = chain_pair_energy(pk1, rho, dz, a1, a2, 0.0)
    assert np.abs(w0 - chain_pair_energy(pk1, rho, -dz, a2 + 180, a1 + 180, 0.0)).max() < 1e-8
    w1 = chain_pair_energy(pk1, rho, dz, a1, a2, 1.0)
    assert np.abs(w1 - chain_pair_energy(pk1, rho, dz, 180 - a2, 180 - a1, 1.0)).max() < 1e-8


# ------------------------------------------------------------------ (b) the table
def test_table_grid_points_match_direct_kernel():
    chain = _chain(PE_T)
    pk1 = CrystalPacker(chain, n_chains=1)
    tab = _table(PE_T, n_angle=24, n_z=4, dr=0.4)
    rng = np.random.default_rng(6)
    n = 40
    ir = rng.integers(0, tab.n_r, n)
    i1, i2 = rng.integers(0, tab.n_angle, n), rng.integers(0, tab.n_angle, n)
    iz, fl = rng.integers(0, tab.n_z, n), rng.integers(0, 2, n)
    r = tab.r_min + tab.dr * ir
    ref = chain_pair_energy(pk1, np.stack([r, np.zeros_like(r)], 1), chain.c * iz / tab.n_z,
                            360.0 * i1 / tab.n_angle, 360.0 * i2 / tab.n_angle, fl)
    got = tab.W[fl, ir, i1, i2, iz]
    uncapped = ref < tab.cap - 1e-3
    assert uncapped.sum() > 10
    assert np.abs(got[uncapped] - ref[uncapped]).max() < 1e-3  # float32 accumulation
    assert (got <= tab.cap + 1e-4).all()
    # interpolation reproduces the stored values at grid points, and vanishes beyond r_max
    interp = tab.interpolate(r, 360.0 * i1 / tab.n_angle, 360.0 * i2 / tab.n_angle, chain.c * iz / tab.n_z, fl)
    assert np.abs(interp - got).max() < 1e-6
    assert tab.interpolate(tab.r_max + 0.5, 30.0, 60.0, 0.3, 0) == 0.0
    # angles and dz are periodic
    assert tab.interpolate(6.0, 30.0, 60.0, 0.3, 0) == pytest.approx(
        tab.interpolate(6.0, 390.0, -300.0, 0.3 + chain.c, 0), abs=1e-9)


def test_table_symmetry_shortcut_agrees_with_full_build():
    chain = _chain(PE_T)
    kw = dict(n_angle=24, n_z=4, dr=0.5)
    full = PairTable.build(chain, symmetry=False, **kw)
    fast = PairTable.build(chain, symmetry=True, **kw)
    assert np.abs(full.W - fast.W).max() < 1e-3


@pytest.mark.parametrize("case", [PE_T, PVDF_TT])
def test_table_energy_matches_direct_kernel(case):
    """Interpolated lattice energies against the exact kernel over random cells.

    The tolerances are for the *coarse* table used here (10 deg angles, c/8 in dz,
    0.1 A in r); the default table (5 deg, c/16, 0.05 A) is an order of magnitude
    better: median 0.005, p90 0.36 kcal/mol per cell over the same sample.
    """
    chain = _chain(case)
    pk2 = CrystalPacker(chain, n_chains=2)
    tab = _table(case, n_angle=36, n_z=8, dr=0.1)
    rng = np.random.default_rng(8)
    params = _random_cells(chain, 60, rng)
    Eref = pk2.energy(params)
    d = np.abs(table_energy(tab, params) - Eref)
    near = Eref < 50.0  # away from the repulsive wall, where the table is capped
    bound = Eref < 0.0
    assert bound.sum() > 10
    med, p90 = np.median(d[near]), np.percentile(d[near], 90)
    assert med < 0.35 and p90 < 3.0, f"E<50: median {med:.3f} p90 {p90:.3f}"
    assert np.median(d[bound]) < 0.15, f"bound cells: median {np.median(d[bound]):.3f}"


# ------------------------------------------------------------------ (c) the FFT screen
def test_fft_grid_equals_table_lookups():
    """The FFT landscape is the same function as table_energy, to float32."""
    chain = _chain(PE_T)  # for its repeat length c
    tab = _table(PE_T, n_angle=24, n_z=8, dr=0.2)
    av, bv = np.arange(4.4, 5.5, 0.25), np.arange(6.8, 7.9, 0.25)
    res = fft_screen(tab, av, bv, flips=(0, 1), n_top=3)
    rows, ref = [], []
    for fi, f in enumerate(res.flips):
        for ia in range(len(av)):
            for ib in range(len(bv)):
                rows.append([av[ia], bv[ib], 90.0, res.phi1[ia, ib, fi], res.phi2[ia, ib, fi], res.dz[ia, ib, fi], f])
                ref.append(res.energy[ia, ib, fi])
    d = np.abs(table_energy(tab, np.array(rows)) - np.array(ref))
    assert d.max() < 1e-3, f"max deviation {d.max():.2e}"
    # the reported grid minimum is the lowest of its own (phi1, phi2, dz) grid
    p = res.top[0].copy()
    for k, step in ((3, 360.0 / tab.n_angle), (4, 360.0 / tab.n_angle), (5, chain.c / tab.n_z)):
        for s in (step, -step):
            q = p.copy()
            q[k] += s
            assert table_energy(tab, q[None])[0] > res.top_energy[0] - 1e-3


# ------------------------------------------------------------------ (d), (e) known cells
def _rank_of_cell(top, a_expt, b_expt, tol=0.3):
    for i, p in enumerate(top):
        a, b = sorted(p[:2])
        if abs(a - min(a_expt, b_expt)) <= tol and abs(b - max(a_expt, b_expt)) <= tol:
            return i
    return -1


def test_pe_screen_finds_orthorhombic_cell_with_both_settings():
    """PE all-trans: the 4.6 x 7.35 cell, with a parallel and a herringbone setting."""
    tab = _table(PE_T, n_angle=36, n_z=8, dr=0.1)
    # flip 1 is dominated by the packer's intra-chain constant for antiparallel cells,
    # so the parallel branch alone decides the ranking here (and is half the work)
    res = fft_screen(tab, np.arange(4.0, 8.81, 0.2), np.arange(4.0, 8.81, 0.2), flips=(0,), n_top=10,
                     min_separation=0.35)
    rank = _rank_of_cell(res.top, 4.6, 7.35, tol=0.3)
    assert rank >= 0, f"known PE cell not in the top 10: {np.round(res.top[:, :2], 2)}"
    dphi = np.abs((res.top[:, 3] - res.top[:, 4] + 90.0) % 180.0 - 90.0)  # setting-angle difference
    assert (dphi[:10] < 25.0).any(), "no parallel setting in the top 10"
    assert ((dphi[:10] > 55.0) & (dphi[:10] < 125.0)).any(), "no herringbone setting in the top 10"


def test_beta_pvdf_screen_finds_known_cell():
    tab = _table(PVDF_TT, n_angle=36, n_z=8, dr=0.1)
    res = fft_screen(tab, np.arange(4.2, 9.41, 0.2), np.arange(4.2, 9.41, 0.2), flips=(0,), n_top=10,
                     min_separation=0.35)
    assert _rank_of_cell(res.top, 4.65, 8.6, tol=0.35) >= 0, f"top cells {np.round(res.top[:, :2], 2)}"


def test_interpolation_survives_the_angle_wrap():
    """A tiny negative angle wraps to exactly 360.0, whose grid index is 0 but whose
    interpolation weight must be 0 and not n_angle.  Taking the weight from the wrapped
    index put -47/+48 on two table entries; oblique lattices hit it routinely, because
    ``arctan2`` puts a site angle a rounding error above a grid angle."""
    tab = _table(PE_T, n_angle=24, n_z=4, dr=0.4)
    assert (-1e-16) % 360.0 == 360.0 and (-1e-17) % tab.c == tab.c  # the trap is real
    for bad, good in ((-1e-16, 0.0), (-1e-17, 0.0)):
        assert tab.interpolate(6.0, bad, 40.0, 0.5, 0) == pytest.approx(
            tab.interpolate(6.0, good, 40.0, 0.5, 0), abs=1e-9)
        assert tab.interpolate(6.0, 40.0, bad, 0.5, 0) == pytest.approx(
            tab.interpolate(6.0, 40.0, good, 0.5, 0), abs=1e-9)
    assert tab.interpolate(6.0, 30.0, 60.0, -1e-17, 0) == pytest.approx(
        tab.interpolate(6.0, 30.0, 60.0, 0.0, 0), abs=1e-9)
    # every interpolated value stays inside the table's own range
    rng = np.random.default_rng(11)
    n = 4000
    v = tab.interpolate(rng.uniform(3.0, tab.r_max, n), rng.uniform(-720, 720, n),
                        rng.uniform(-720, 720, n), rng.uniform(-3 * tab.c, 3 * tab.c, n),
                        rng.integers(0, 2, n))
    assert v.min() >= tab.W.min() - 1e-4 and v.max() <= tab.cap + 1e-4


@pytest.mark.parametrize("gamma", [60.0, 75.0, 105.0, 120.0])
def test_oblique_cells_score_correctly(gamma):
    """Oblique lattices put site angles on the grid to within a rounding error; the screen
    must still agree with the exact kernel there (it used to be hundreds of kcal/mol low)."""
    chain = _chain(PE_T)
    pk2 = CrystalPacker(chain, n_chains=2)
    tab = _table(PE_T, n_angle=36, n_z=8, dr=0.1)
    av = np.arange(4.4, 9.61, 0.4)
    res = fft_screen(tab, av, av, gamma=gamma, flips=(0, 1), n_top=4)
    assert np.isfinite(res.top_energy).all()
    d = np.abs(pk2.energy(res.top) - res.top_energy)
    assert d.max() < 1.0, f"screen vs exact: {np.round(d, 3)}"
    assert np.abs(table_energy(tab, res.top) - res.top_energy).max() < 1e-3


def test_ab_symmetry_halves_the_grid_without_changing_the_answer():
    """(a, b) and (b, a) are the same lattice rotated by 90 deg, so half the FFTs are free."""
    tab = _table(PE_T, n_angle=24, n_z=4, dr=0.4)
    av = np.arange(4.4, 8.01, 0.4)
    full = fft_screen(tab, av, av, flips=(0, 1), n_top=5)
    half = fft_screen(tab, av, av, flips=(0, 1), n_top=5, ab_symmetry=True)
    assert np.abs(half.top_energy - full.top_energy).max() < 1e-4
    assert np.isinf(half.energy[-1, 0, 0]) and np.isfinite(half.energy[0, -1, 0])
    for p in half.top:
        assert p[1] >= p[0] - 1e-9
    # every screened cell agrees with its transpose in the full screen
    fin = np.isfinite(half.energy)
    assert np.abs(half.energy[fin] - np.swapaxes(full.energy, 0, 1)[fin]).max() < 1e-3
    with pytest.raises(ValueError):
        fft_screen(tab, av, av, gamma=100.0, ab_symmetry=True)


def test_screen_wrapper_returns_rows_with_energies():
    chain = _chain(PE_T)
    tab = _table(PE_T, n_angle=24, n_z=4, dr=0.4)
    params, energies, res = screen(chain, n_top=5, a_range=(4.5, 5.5), b_range=(6.5, 7.5), da=0.25, db=0.25,
                                   flips=(0,), table=tab)
    assert params.shape == (5, 7) and energies.shape == (5,)
    assert (np.diff(energies) >= -1e-9).all()
    assert np.abs(table_energy(tab, params) - energies).max() < 1e-3
    assert res.energy.shape == (5, 5, 1)
