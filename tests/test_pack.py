import dataclasses

import numpy as np
import pytest

from polyfind.forcefield import SimpleFF
from polyfind.pack import (
    EV_TO_KCAL,
    CrystalPacker,
    chain_valence,
    default_bounds,
    pack,
    periodic_chain,
    periodic_chain_from_torsions,
    polish,
    repeat_chains_from_torsions,
    to_cif,
)
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


def test_polish_analytic_and_fd_gradients_reach_the_same_minimum():
    """The analytic polish lands where the finite-difference one does, in far fewer rows."""
    ch = periodic_chain(PVDF, [T, GP, T, GM], THREE_STATE)
    packer = CrystalPacker(ch, n_chains=2)
    b = default_bounds(ch)
    keys = ["a", "b", "gamma", "phi1", "phi2", "dz"]
    lo = np.array([b[k][0] for k in keys])
    hi = np.array([b[k][1] for k in keys])
    free = [i for i in range(6) if hi[i] > lo[i]]
    for start in ([5.3, 9.0, 90.0, 209.0, 209.0, 0.1, 0.0], [5.0, 9.8, 90.0, 40.0, 220.0, 1.5, 1.0]):
        rows = {}
        x = {}
        for g in ("fd", "analytic"):
            packer.n_energy_rows = 0
            x[g] = polish(packer, np.array(start), lo, hi, free, 1500, gradient=g)
            rows[g] = packer.n_energy_rows
        e = {g: float(packer.energy(x[g][None])[0]) for g in x}
        assert e["analytic"] <= e["fd"] + 1e-6
        assert np.abs(x["analytic"][:2] - x["fd"][:2]).max() < 0.01
        for i in (3, 4):
            assert abs(((x["analytic"][i] - x["fd"][i] + 180.0) % 360.0) - 180.0) < 0.2
        # 1 row per evaluation instead of 1 + 2 * 5
        assert rows["analytic"] * 5 <= rows["fd"]
    with pytest.raises(ValueError):
        polish(packer, np.array(start), lo, hi, free, gradient="autodiff")


@pytest.mark.parametrize("poly,seq", [
    (PE, [T]),
    (PVDF, [T, T]),
    (PVDF, [T, GP, T, GM]),
    (PVDF, [T, T, T, GP, T, T, T, GM]),
])
def test_energy_and_grad_value_is_bit_identical_to_energy(poly, seq):
    ch = periodic_chain(poly, seq, THREE_STATE)
    rng = np.random.default_rng(4)
    tors = ch.dihedrals + rng.uniform(-5, 5, len(ch.dihedrals))
    other = repeat_chains_from_torsions(poly, ch.name, tors[None], align_to=ch)[0]
    for field in (None, (0.14, -0.06, 0.19)):
        for n_chains in (1, 2):
            pk = CrystalPacker(ch, n_chains=n_chains, field=field)
            for flip in (0.0, 1.0) if n_chains == 2 else (0.0,):
                p = np.array([5.1, 9.3, 97.0, 37.0, 212.0, 0.4 * ch.c, flip])
                E, _, _, _ = pk.energy_and_grad(p)
                assert E == float(pk.energy(p[None])[0])  # exactly, not to a tolerance
                kw = dict(coords=other.coords, c=np.array([other.c]))
                E2, _, _, _ = pk.energy_and_grad(p, **kw)
                assert E2 == float(pk.energy(p[None], **kw)[0])
    with pytest.raises(ValueError):
        CrystalPacker(ch).energy_and_grad(np.zeros((2, 7)))


@pytest.mark.parametrize("poly,seq", [
    (PE, [T]),
    (PVDF, [T, T]),
    (PVDF, [T, GP, T, GM]),
    (PVDF, [T, T, T, GP, T, T, T, GM]),
])
def test_energy_and_grad_matches_central_differences(poly, seq):
    """Every component of the analytic gradient, against a central difference.

    The steps are swept and the closest agreement is taken: the Lennard-Jones term is
    energy-shifted but not force-shifted, so a step that straddles the cutoff sees a jump
    in the *difference* that the derivative correctly does not have.  Away from that, the
    agreement is at the finite-difference floor, ~1e-6 relative.
    """
    ch = periodic_chain(poly, seq, THREE_STATE)
    rng = np.random.default_rng(7)
    tors = ch.dihedrals + rng.uniform(-5, 5, len(ch.dihedrals))
    other = repeat_chains_from_torsions(poly, ch.name, tors[None], align_to=ch)[0]
    coords, cz = other.coords, np.array([other.c])
    for field in (None, (0.14, -0.06, 0.19)):
        pk = CrystalPacker(ch, field=field)
        for flip in (0.0, 1.0):
            p = np.array([5.1, 9.3, 97.0, 37.0, 212.0, 0.4 * ch.c, flip])
            _, g_cell, g_coords, g_c = pk.energy_and_grad(p, coords=coords, c=cz)
            scale = max(np.abs(g_cell).max(), np.abs(g_coords).max(), 1e-6)

            def fd(bump, steps):
                best = None
                for h in steps:
                    d = (bump(h) - bump(-h)) / (2 * h)
                    if best is None or abs(d - ref) < abs(best - ref):
                        best = d
                return best

            def energy_of(pp=None, cc=None, zz=None):
                return float(pk.energy((p if pp is None else pp)[None],
                                       coords=coords if cc is None else cc,
                                       c=cz if zz is None else zz)[0])

            # 1e-4 of the component, plus 1e-6 of the gradient's own scale so that a
            # component near zero is not held to an absolute impossibility
            def close(got, ref):
                assert abs(got - ref) <= 1e-4 * abs(ref) + 1e-6 * scale, (got, ref)

            for i in range(6):
                steps = (3e-2, 1e-2, 3e-3, 1e-3) if i in (2, 3, 4) else (3e-3, 1e-3, 3e-4, 1e-4)
                ref = g_cell[i]
                close(fd(lambda h, i=i: energy_of(pp=p + h * np.eye(7)[i]), steps), ref)
            for i in rng.choice(coords.shape[0], min(4, coords.shape[0]), replace=False):
                for d in range(3):
                    ref = g_coords[i, d]
                    close(fd(lambda h, i=int(i), d=d: energy_of(cc=coords + h * np.eye(coords.size)[
                        3 * i + d].reshape(coords.shape)), (3e-3, 1e-3, 3e-4, 1e-4)), ref)
            ref = g_c
            close(fd(lambda h: energy_of(zz=cz + h), (3e-3, 1e-3, 3e-4, 1e-4)), ref)


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


# ------------------------------------------------------------------ dipole and field coupling
# Lattice energies measured with the fieldless code (commit 12a817d) at fixed cells, as a
# regression guard that adding the field term left the zero-field kernel alone.
#
# The three PVDF rows were RE-RECORDED when PVDF's C-C bond length was corrected from
# 1.54 to the DFT Form I 1.528 A (docs/REFERENCES.md): that moves every atom of the
# chain, so the energies at these fixed cells necessarily move with it.  PE's row is
# untouched -- PE's geometry did not change -- and it is what still holds the original
# recording, and hence the kernel itself, in place.
FIELDLESS_ENERGIES = {
    "pe": (-5.610711164309272, -4.205273674665557),
    "beta": (-0.5491161941046234, -9.17687306724556),
    "alpha": (96.67441893488612, 16.146841685705),
    "gamma": (87.07041823166938, 159.77466693085987),
}


def _fixed_cells(c):
    return np.array([[4.9, 8.7, 90.0, 30.0, 210.0, 0.4 * c, 0.0],
                     [5.3, 9.1, 97.0, 145.0, 25.0, 0.8 * c, 1.0]])


@pytest.mark.parametrize("tag,poly,seq", [
    ("pe", PE, [T]), ("beta", PVDF, [T, T]),
    ("alpha", PVDF, [T, GP, T, GM]), ("gamma", PVDF, [T, T, T, GP, T, T, T, GM]),
])
def test_zero_field_energies_are_exactly_the_fieldless_ones(tag, poly, seq):
    """No field, no change: the kernel must not even see the coupling term."""
    ch = periodic_chain(poly, seq, THREE_STATE)
    p = _fixed_cells(ch.c)
    e = CrystalPacker(ch).energy(p)
    np.testing.assert_allclose(e, FIELDLESS_ENERGIES[tag], rtol=0, atol=1e-9)
    # an explicitly zero field is the same code path and gives bit-identical numbers
    for field in (None, (0.0, 0.0, 0.0), np.zeros(3)):
        assert np.array_equal(CrystalPacker(ch, field=field).energy(p), e)
    assert np.array_equal(CrystalPacker(ch, field=(0.0, 0.0, 0.0)).field_energy(p), np.zeros(2))
    assert CrystalPacker(ch).result(p[0]).field_energy == 0.0


def test_cell_dipole_is_neutral_origin_and_shift_independent():
    ch = periodic_chain(PVDF, [T, GP, T, GM], THREE_STATE)
    pk = CrystalPacker(ch)
    assert pk.cell_charge == pytest.approx(0.0, abs=1e-12) and pk.is_neutral
    p = _fixed_cells(ch.c)
    mu = pk.dipole(p)
    # dz translates a neutral chain, which cannot change the cell dipole
    q = p.copy()
    q[:, 5] += 1.234
    np.testing.assert_allclose(pk.dipole(q), mu, atol=1e-12)
    # rotating both chains together rotates mu rigidly about z
    d = 37.0
    r = p.copy()
    r[:, 3] += d
    r[:, 4] += d
    ca, sa = np.cos(np.deg2rad(d)), np.sin(np.deg2rad(d))
    expect = np.column_stack([ca * mu[:, 0] - sa * mu[:, 1], sa * mu[:, 0] + ca * mu[:, 1], mu[:, 2]])
    np.testing.assert_allclose(pk.dipole(r), expect, atol=1e-10)
    # P = mu / V with 1 e/A^2 = 16.0218 C/m^2
    V = p[:, 0] * p[:, 1] * np.sin(np.deg2rad(p[:, 2])) * ch.c
    np.testing.assert_allclose(pk.polarization(p), mu / V[:, None] * 16.0218, rtol=1e-12)


def test_a_charged_cell_has_no_dipole_and_no_field():
    """sum q_i r_i moves with the origin unless the cell is neutral, so it is refused."""
    ch = periodic_chain(PVDF, [T, T], THREE_STATE)
    charged = dataclasses.replace(ch, charges=ch.charges + 0.05)
    pk = CrystalPacker(charged)
    assert not pk.is_neutral and pk.cell_charge == pytest.approx(2 * 0.05 * ch.n_atoms)
    with pytest.raises(ValueError, match="origin"):
        pk.dipole(_fixed_cells(ch.c))
    with pytest.raises(ValueError, match="origin"):
        CrystalPacker(charged, field=(0.1, 0.0, 0.0))
    assert pk.result(_fixed_cells(ch.c)[0]).polarization is None


def test_pe_polarization_is_exactly_zero_by_symmetry():
    """All-trans PE has a CH2 dipole at every backbone atom and they cancel in pairs.

    The repeat holds two CH2 groups whose bisectors point opposite ways, so the sum is
    zero identically -- for every setting angle, shift and flip, not just at the minimum.
    """
    ch = periodic_chain(PE, [T], THREE_STATE)
    for n_chains in (1, 2):
        pk = CrystalPacker(ch, n_chains=n_chains)
        rng = np.random.default_rng(3)
        M = 10
        p = np.column_stack([rng.uniform(4.5, 9, M), rng.uniform(4.5, 9, M), rng.uniform(80, 100, M),
                             rng.uniform(0, 360, M), rng.uniform(0, 360, M), rng.uniform(0, ch.c, M),
                             rng.integers(0, 2, M).astype(float)])
        assert np.abs(pk.dipole(p)).max() < 1e-12
        assert np.abs(pk.polarization(p)).max() < 1e-12
        # ... so no field can do work on a PE cell either
        assert np.abs(CrystalPacker(ch, n_chains=n_chains, field=(0.4, -0.3, 0.2)).energy(p) - pk.energy(p)).max() < 1e-9
    assert pk.result(p[0]).polarization_magnitude < 1e-12


def test_a_field_along_the_dipole_lowers_the_energy_and_against_it_raises_it_equally():
    ch = periodic_chain(PVDF, [T, T], THREE_STATE)
    p = np.array([[4.65, 8.61, 90.0, 0.0, 0.0, 1.29, 0.0]])
    pk = CrystalPacker(ch)
    mu = pk.dipole(p)[0]
    assert np.linalg.norm(mu) > 0.1  # beta really is polar
    u = mu / np.linalg.norm(mu)
    e0 = float(pk.energy(p)[0])
    s = 0.25
    down = float(CrystalPacker(ch, field=s * u).energy(p)[0])
    up = float(CrystalPacker(ch, field=-s * u).energy(p)[0])
    expect = s * np.linalg.norm(mu) * EV_TO_KCAL
    assert down < e0 < up
    assert e0 - down == pytest.approx(expect, rel=1e-9)
    assert up - e0 == pytest.approx(expect, rel=1e-9)
    # reversing the field reverses the coupling; doubling it doubles it; and it is
    # exactly the difference the kernel makes
    for E in (np.array([0.3, -0.1, 0.05]), s * u):
        pk_f = CrystalPacker(ch, field=E)
        pk_r = CrystalPacker(ch, field=-E)
        pk_2 = CrystalPacker(ch, field=2 * E)
        w = float(pk_f.field_energy(p)[0])
        assert w != pytest.approx(0.0, abs=1e-6)
        assert float(pk_r.field_energy(p)[0]) == pytest.approx(-w, rel=1e-12)
        assert float(pk_2.field_energy(p)[0]) == pytest.approx(2 * w, rel=1e-12)
        assert float(pk_f.energy(p)[0]) - e0 == pytest.approx(w, abs=1e-9)
        assert pk_f.result(p[0]).field_energy == pytest.approx(w, abs=1e-9)


def test_the_field_term_follows_the_configuration_not_the_chain_alone():
    """The dipole must be recomputed per row: setting angle and flip move the atoms."""
    ch = periodic_chain(PVDF, [T, GP, T, GM], THREE_STATE)
    pk = CrystalPacker(ch, field=(0.2, 0.1, -0.15))
    rng = np.random.default_rng(5)
    M = 9
    p = np.column_stack([rng.uniform(4.8, 9, M), rng.uniform(4.8, 9, M), np.full(M, 90.0),
                         rng.uniform(0, 360, M), rng.uniform(0, 360, M), rng.uniform(0, ch.c, M),
                         rng.integers(0, 2, M).astype(float)])
    w = pk.field_energy(p)
    assert w.std() > 0.1  # it genuinely varies over the batch
    # one batched call agrees with row-by-row calls, and with energy() minus the fieldless energy
    for i in range(M):
        assert float(pk.field_energy(p[i : i + 1])[0]) == pytest.approx(w[i], rel=1e-12)
    np.testing.assert_allclose(pk.energy(p) - CrystalPacker(ch).energy(p), w, atol=1e-9)
    # per-row chain geometries carry their own dipoles too
    chains = repeat_chains_from_torsions(PVDF, "TG+TG-", ch.dihedrals + rng.uniform(-10, 10, (M, 4)), align_to=ch)
    kw = pk.chain_batch(chains)
    got = pk.energy(p, **kw) - CrystalPacker(ch).energy(p, **kw)
    np.testing.assert_allclose(got, pk.field_energy(p, coords=kw["coords"], c=kw["c"]), atol=1e-9)


def test_pack_honours_the_field_and_reports_the_polarization():
    ch = periodic_chain(PVDF, [T, T], THREE_STATE)
    kw = dict(n_chains=2, n_random=500, n_refine=2, maxfev=400, screen="random")
    zero = pack(ch, rng=np.random.default_rng(0), **kw)[0]
    assert zero.polarization is not None and zero.field == (0.0, 0.0, 0.0)
    assert f"|P|={zero.polarization_magnitude:5.3f}" in zero.row()
    u = zero.polarization / np.linalg.norm(zero.polarization)
    with_field = pack(ch, rng=np.random.default_rng(0), field=0.3 * u, **kw)[0]
    assert with_field.field == pytest.approx(tuple(0.3 * u))
    assert with_field.field_energy < 0.0  # the cell it picks is aligned with the field
    # the reported energy is the field-aware kernel's own value
    pk = CrystalPacker(ch, field=0.3 * u)
    assert float(pk.energy(with_field.params[None])[0]) == pytest.approx(with_field.energy_per_cell, abs=1e-9)
    assert with_field.energy_per_cell < zero.energy_per_cell


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


# --- valence terms on the lattice side (opt-in) -----------------------------------------
def _valence_ff():
    import polyfind.fitting  # noqa: F401 - registers the fitted presets

    return SimpleFF.from_preset("pvdf-dft-valence")


@pytest.mark.parametrize("poly,seq,n_bonds,n_angles", [
    # one repeat of PE is two CH2 groups: 2 C-C bonds + 4 C-H, and 6 angles at each carbon
    (PE, [T], 6, 12),
    (PVDF, [T, T], 6, 12),
    (PVDF, [T, GP, T, GM], 12, 24),
])
def test_periodic_valence_counts_one_representative_of_every_term(poly, seq, n_bonds, n_angles):
    """The bonds and angles of the *infinite* chain, counted once each per repeat.

    A backbone bond joins the last atom of the repeat to the first atom of the next one, so
    the terms do not stop at the block's edges and the count is a real check on the image
    bookkeeping rather than a formality.
    """
    ch = periodic_chain(poly, seq, THREE_STATE)
    ff = _valence_ff()
    cv = chain_valence(ch, ff.bond_table(), ff.angle_table())
    assert (cv.n_bonds, cv.n_angles) == (n_bonds, n_angles)
    assert set(cv.bond_s) <= {-1.0, 0.0, 1.0} and set(cv.ang_si) <= {-1.0, 0.0, 1.0}
    assert (cv.ang_j < ch.n_atoms).all() and (cv.bond_i < ch.n_atoms).all()
    # translating the repeat cannot change it, and neither can a whole-cell rotation about z
    e0 = cv.energy(ch.coords[None], np.array([ch.c]))[0]
    shifted = ch.coords + np.array([0.31, -0.22, 0.53])
    assert cv.energy(shifted[None], np.array([ch.c]))[0] == pytest.approx(e0, abs=1e-12)
    assert chain_valence(ch, {}, {}) is None


@pytest.mark.parametrize("poly,seq", [(PE, [T]), (PVDF, [T, T]), (PVDF, [T, GP, T, GM])])
def test_periodic_valence_gradient_matches_central_differences(poly, seq):
    ch = periodic_chain(poly, seq, THREE_STATE)
    ff = _valence_ff()
    cv = chain_valence(ch, ff.bond_table(), ff.angle_table())
    E, gX, gc = cv.energy_and_grad(ch.coords, ch.c)
    assert E == pytest.approx(float(cv.energy(ch.coords[None], np.array([ch.c]))[0]), abs=1e-12)
    h = 1e-6
    for i in range(ch.n_atoms):
        for d in range(3):
            X = ch.coords.copy()
            X[i, d] += h
            ep = cv.energy(X[None], np.array([ch.c]))[0]
            X[i, d] -= 2 * h
            em = cv.energy(X[None], np.array([ch.c]))[0]
            assert gX[i, d] == pytest.approx((ep - em) / (2 * h), abs=1e-6)
    fd = (cv.energy(ch.coords[None], np.array([ch.c + h]))[0]
          - cv.energy(ch.coords[None], np.array([ch.c - h]))[0]) / (2 * h)
    assert gc == pytest.approx(fd, abs=1e-6)
    assert abs(gc) > 1e-3  # c really does have a restoring force now


def test_valence_terms_change_the_energy_only_when_asked_for():
    """The default packer is untouched, and an empty valence table changes nothing either."""
    ch = periodic_chain(PVDF, [T, T], THREE_STATE)
    p = np.array([[5.0, 9.6, 92.0, 20.0, 50.0, 1.0, 0.0]])
    plain = float(CrystalPacker(ch).energy(p)[0])
    assert float(CrystalPacker(ch, valence=SimpleFF()).energy(p)[0]) == plain  # no terms: no change
    pk = CrystalPacker(ch, valence=_valence_ff())
    assert pk.e_valence > 0.0
    assert float(pk.energy(p)[0]) == pytest.approx(plain + pk.e_valence, abs=1e-9)
    # and the gradient path returns the same value bit for bit, as it does without them
    assert pk.energy_and_grad(p[0])[0] == float(pk.energy(p)[0])


@pytest.mark.parametrize("poly,seq", [(PE, [T]), (PVDF, [T, GP, T, GM])])
def test_valence_gradient_reaches_the_kernel_gradient(poly, seq):
    """``energy_and_grad`` carries the valence derivative in ``g_coords`` and ``g_c``.

    That is the whole of the plumbing: a refinement's chain rule and the axial relaxation
    read those two arrays and nothing else, so if the valence term is in them it is in
    every consumer.
    """
    ch = periodic_chain(poly, seq, THREE_STATE)
    ff = _valence_ff()
    pk_plain, pk_val = CrystalPacker(ch), CrystalPacker(ch, valence=ff)
    cv = chain_valence(ch, ff.bond_table(), ff.angle_table())
    p = np.array([5.1, 9.3, 97.0, 37.0, 212.0, 0.4 * ch.c, 1.0])
    kw = dict(coords=ch.coords, c=np.array([ch.c]))
    _, gc0, gX0, gz0 = pk_plain.energy_and_grad(p, **kw)
    _, gc1, gX1, gz1 = pk_val.energy_and_grad(p, **kw)
    _, vX, vz = cv.energy_and_grad(ch.coords, ch.c)
    assert gc1 == pytest.approx(gc0, abs=1e-9)  # the cell gradient cannot see an internal term
    assert gX1 - gX0 == pytest.approx(pk_val.n_chains * vX, abs=1e-9)
    assert gz1 - gz0 == pytest.approx(pk_val.n_chains * vz, abs=1e-9)


# --- the Lennard-Jones cutoff -------------------------------------------------------------
def test_force_shifted_lennard_jones_reaches_the_cutoff_with_zero_force():
    """Value *and* force vanish at ``rc``; with the default form only the value does."""
    ch = periodic_chain(PVDF, [T, T], THREE_STATE)
    r2 = np.array([[(8.0 - 1e-6) ** 2]])
    for shift, force_at_cutoff in (("energy", False), ("force", True)):
        pk = CrystalPacker(ch, lj_cutoff=shift)
        args = (pk._A_nn[:1, :1], pk._B_nn[:1, :1], pk._qq_nn[:1, :1], pk._qqf_nn[:1, :1],
                pk._const_nn[:1, :1])
        v, dv = pk._pair_energy_and_dv(r2, *args)
        assert abs(float(v[0, 0])) < 1e-7  # both forms are energy-continuous
        # dV/dr = 2 r dV/d(r^2); the LJ well depth sets the scale of what "not zero" means
        dvdr = abs(2.0 * 8.0 * float(dv[0, 0]))
        assert (dvdr < 1e-6) == force_at_cutoff
    with pytest.raises(ValueError, match="unknown lj_cutoff"):
        CrystalPacker(ch, lj_cutoff="switch")


def test_force_shift_changes_the_energy_and_only_when_asked_for():
    ch = periodic_chain(PVDF, [T, T], THREE_STATE)
    p = np.array([[4.6, 8.6, 90.0, 0.0, 0.0, 1.0, 0.0]])
    e_shift = float(CrystalPacker(ch).energy(p)[0])
    assert float(CrystalPacker(ch, lj_cutoff="energy").energy(p)[0]) == e_shift
    assert abs(float(CrystalPacker(ch, lj_cutoff="force").energy(p)[0]) - e_shift) > 0.1


# --- the rigid table and the deformable kernel stay apart ----------------------------------
def test_the_table_screen_refuses_a_packer_it_cannot_describe():
    """The tabulated interaction is rigid and energy-shifted; anything else is refused.

    ``chain_pair_energy`` rebuilds the pair potential from ``lj_shift`` and the two DSF
    constants and knows nothing about a valence term or a force-shifted tail, so a table
    screen of such a packer would be screening a different potential from the one polished.
    """
    ch = periodic_chain(PE, [T], THREE_STATE)
    kw = dict(n_chains=2, n_refine=1, maxfev=200, screen="table")
    from polyfind.pack import CrystalPacker as _CP, _table_starts

    for bad in (dict(valence=_valence_ff()), dict(lj_cutoff="force")):
        pk = _CP(ch, n_chains=2, **bad)
        with pytest.raises(ValueError, match="screen='table' cannot be used"):
            _table_starts(pk, ch, np.array([4.0, 6.0, 90.0, 0.0, 0.0, 0.0]),
                          np.array([5.0, 8.0, 90.0, 360.0, 360.0, ch.c]), 1, [0, 1], 0.5, [90.0],
                          None, 8.0, 0.2, 1.0, None, None, False)
    assert pack(ch, **kw)  # and the rigid, energy-shifted packer is still screened by table


def test_the_table_screen_refuses_a_table_built_for_another_conformation():
    from polyfind.lattice_table import pair_table

    ch = periodic_chain(PVDF, [T, T], THREE_STATE)
    bent = repeat_chains_from_torsions(PVDF, ch.name, (ch.dihedrals - 6.0)[None], align_to=ch)[0]
    table = pair_table(ch, cutoff=8.0, cache_dir=None, n_angle=16, n_z=4, dr=0.4)
    assert pack(ch, n_refine=1, maxfev=200, screen="table", table=table)
    with pytest.raises(ValueError, match="different chain coordinates"):
        pack(bent, n_refine=1, maxfev=200, screen="table", table=table)
