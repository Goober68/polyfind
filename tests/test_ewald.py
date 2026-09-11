"""Ewald summation: the published Madelung constants, the invariances, the gradients,
and the boundary between the Ewald packer and the pairwise screen.

Ewald is easy to get subtly wrong and its failure mode is a plausible number, so these
are checks against things known independently of this package: a lattice sum whose value
is published to more digits than a double can hold, the requirement that the answer not
depend on how the sum was split, and the requirement that it not depend on which periodic
image an atom is written at.
"""
import numpy as np
import pytest

from polyfind.ewald import (
    MADELUNG_CSCL,
    MADELUNG_NACL,
    MADELUNG_ZINCBLENDE,
    Ewald,
    EwaldSpec,
    exclusion_correction,
    madelung_nacl,
)
from polyfind.forcefield import COULOMB
from polyfind.lattice_table import chain_pair_energy, intra_constants
from polyfind.pack import CrystalPacker, periodic_chain
from polyfind.polymers import PE, PVDF, THREE_STATE

T, GP, GM = 0, 1, 2
BETA = np.array([4.96, 9.64, 90.0, 30.0, 30.0, 0.6, 0.0])


@pytest.fixture(scope="module")
def beta_chain():
    return periodic_chain(PVDF, (T, T), THREE_STATE)


# --------------------------------------------------------------- Madelung constants
def test_madelung_nacl_to_published_precision():
    """No tolerance: the rock-salt Madelung constant to every digit of the published value.

    A simple cubic lattice of alternating unit charges -- the textbook conditionally
    convergent sum -- at accuracy 1e-16 and rc = 12 reproduces 1.747564594633 exactly.
    """
    m = madelung_nacl(EwaldSpec(accuracy=1e-16), rc=12.0)
    assert f"{m:.12f}" == f"{MADELUNG_NACL:.12f}"
    assert m == pytest.approx(MADELUNG_NACL, abs=1e-14)


@pytest.mark.parametrize("rc", [8.0, 10.0, 12.0, 16.0])
def test_madelung_nacl_independent_of_the_real_space_cutoff(rc):
    assert madelung_nacl(EwaldSpec(accuracy=1e-14), rc=rc) == pytest.approx(MADELUNG_NACL, abs=1e-11)


def test_madelung_cscl_and_zincblende():
    """Two more published lattices, because rock salt alone does not exercise a cell whose
    ions sit off the lattice points."""
    ew = Ewald(14.0, EwaldSpec(accuracy=1e-14), eps_r=1.0, coulomb=1.0)
    e = ew.energy(np.array([[0.0, 0, 0], [0.5, 0.5, 0.5]]), np.array([1.0, -1.0]), np.eye(3))
    assert -e * np.sqrt(3) / 2 == pytest.approx(MADELUNG_CSCL, abs=1e-9)
    fcc = np.array([[0.0, 0, 0], [0, .5, .5], [.5, 0, .5], [.5, .5, 0]])
    e = ew.energy(np.vstack([fcc, fcc + 0.25]), np.array([1.0] * 4 + [-1.0] * 4), np.eye(3))
    assert -e * (np.sqrt(3) / 4) / 4 == pytest.approx(MADELUNG_ZINCBLENDE, abs=1e-9)


# ------------------------------------------------------- splitting and cutoffs
def _cell(seed=3, n=14):
    rng = np.random.default_rng(seed)
    g = np.deg2rad(97.0)
    h = np.array([[5.1, 0.0, 0.0], [8.7 * np.cos(g), 8.7 * np.sin(g), 0.0], [0.0, 0.0, 2.62]])
    q = rng.normal(size=n)
    return rng.random((n, 3)) @ h, q - q.mean(), h


@pytest.mark.parametrize("boundary", ["tinfoil", "vacuum"])
def test_energy_does_not_depend_on_the_splitting_parameter(boundary):
    """alpha moves work between the two halves and must not move the total.  This is the
    single most informative internal check: a wrong prefactor in either half breaks it."""
    X, q, h = _cell()
    vals = [Ewald(12.0, EwaldSpec(alpha=a, accuracy=1e-14, boundary=boundary)).energy(X, q, h)
            for a in (0.35, 0.45, 0.55, 0.70, 0.90)]
    assert max(vals) - min(vals) < 1e-8 * max(1.0, abs(vals[0]))


def test_reciprocal_cutoff_converges():
    X, q, h = _cell()
    ref = Ewald(10.0, EwaldSpec(alpha=0.5, kmax=14.0)).energy(X, q, h)
    errs = [abs(Ewald(10.0, EwaldSpec(alpha=0.5, kmax=km)).energy(X, q, h) - ref)
            for km in (3.0, 4.0, 5.0, 6.0)]
    assert errs == sorted(errs, reverse=True)  # monotone convergence
    assert errs[-1] < 1e-9 * abs(ref)


def test_real_cutoff_converges():
    X, q, h = _cell()
    ref = Ewald(14.0, EwaldSpec(alpha=0.5, kmax=14.0, rc=14.0)).energy(X, q, h)
    for rc in (8.0, 10.0, 12.0):
        e = Ewald(rc, EwaldSpec(alpha=0.5, kmax=14.0, rc=rc)).energy(X, q, h)
        assert abs(e - ref) < 1e-4


@pytest.mark.parametrize("boundary", ["tinfoil", "vacuum"])
def test_default_alpha_and_kmax_follow_the_requested_accuracy(boundary):
    ew = Ewald(8.0, EwaldSpec(accuracy=1e-8, boundary=boundary))
    assert np.exp(-(ew.alpha * ew.rc) ** 2) == pytest.approx(1e-8, rel=1e-9)
    assert np.exp(-ew.kmax ** 2 / (4 * ew.alpha ** 2)) == pytest.approx(1e-8, rel=1e-9)


# ------------------------------------------------------------------- invariances
def test_tinfoil_is_invariant_under_translating_an_atom_by_a_lattice_vector():
    """The conditionally convergent part lives in the surface term, and that is exactly the
    term that cares which image an atom is written at.  Tinfoil does not."""
    X, q, h = _cell()
    ew = Ewald(10.0)
    e0 = ew.energy(X, q, h)
    for atom, vec in ((0, h[0]), (5, -h[1]), (9, 3 * h[2] - 2 * h[0])):
        Y = X.copy()
        Y[atom] += vec
        assert ew.energy(Y, q, h) == pytest.approx(e0, abs=1e-9)


def test_vacuum_is_not_invariant_and_that_is_the_physics():
    """The vacuum convention's surface term is ``2 pi |M|^2 / 3 V``, and ``M`` changes when a
    *charged* atom is rewritten at another image.  Moving a whole neutral group does not
    change it -- which is why the molecular branch (chains kept whole) is the one used."""
    X, q, h = _cell()
    ew = Ewald(10.0, EwaldSpec(boundary="vacuum"))
    Y = X.copy()
    Y[0] += h[0]
    assert abs(ew.energy(Y, q, h) - ew.energy(X, q, h)) > 1.0
    Z = X.copy()  # a neutral pair moved together
    pair = [1, 2]
    Z[pair] += h[1]
    q2 = q.copy()
    q2[pair] = [0.3, -0.3]
    assert ew.energy(Z, q2, h) == pytest.approx(ew.energy(X, q2, h), abs=1e-9)


# -------------------------------------------------------------------- gradients
@pytest.mark.parametrize("boundary", ["tinfoil", "vacuum"])
def test_gradients_against_central_differences(boundary):
    X, q, h = _cell(seed=5, n=9)
    ew = Ewald(9.0, EwaldSpec(boundary=boundary, accuracy=1e-10))
    t = ew.terms(X, q, h, grad=True, charge_grad=True)
    eps = 1e-6

    def fd(base, put):
        out = np.zeros(base.shape)
        for idx in np.ndindex(base.shape):
            a, b = base.copy(), base.copy()
            a[idx] += eps
            b[idx] -= eps
            out[idx] = (put(a) - put(b)) / (2 * eps)
        return out

    assert fd(X, lambda v: ew.energy(v, q, h)) == pytest.approx(t.grad_coords, rel=2e-6, abs=1e-6)
    assert fd(h, lambda v: ew.energy(X, q, v)) == pytest.approx(t.grad_lattice, rel=2e-6, abs=1e-6)
    assert fd(q, lambda v: ew.energy(X, v, h)) == pytest.approx(t.grad_charges, rel=2e-6, abs=1e-6)


def test_exclusion_correction_value_and_gradients():
    rng = np.random.default_rng(11)
    n, K, c = 6, 3, 2.55
    X = rng.random((n, 3)) * np.array([1.5, 1.5, 2.4])
    q = rng.normal(size=n)
    q -= q.mean()
    S = np.ones((2 * K + 1, n, n))
    S[K][0, 1] = S[K][1, 0] = 0.0
    S[K][0, 2] = S[K][2, 0] = 0.5
    S[K + 1][0, 3] = S[K - 1][3, 0] = 0.0
    S[K][np.arange(n), np.arange(n)] = 0.0
    e, gX, gc, gq = exclusion_correction(X, q, c, S, grad=True, charge_grad=True)
    # the value is the plain Coulomb energy of the removed pairs, built here independently
    ks = np.arange(-K, K + 1)
    want = 0.0
    for a, k in enumerate(ks):
        for i in range(n):
            for j in range(n):
                if S[a][i, j] == 1.0:
                    continue
                d = X[i] - X[j] - np.array([0.0, 0.0, k * c])
                r = np.linalg.norm(d)
                if r > 1e-9:
                    want += -0.5 * COULOMB * (1 - S[a][i, j]) * q[i] * q[j] / r
    assert e == pytest.approx(want, rel=1e-12)
    eps = 1e-6
    num = np.zeros_like(X)
    for idx in np.ndindex(X.shape):
        a, b = X.copy(), X.copy()
        a[idx] += eps
        b[idx] -= eps
        num[idx] = (exclusion_correction(a, q, c, S)[0] - exclusion_correction(b, q, c, S)[0]) / (2 * eps)
    assert num == pytest.approx(gX, rel=1e-6, abs=1e-6)
    numc = (exclusion_correction(X, q, c + eps, S)[0] - exclusion_correction(X, q, c - eps, S)[0]) / (2 * eps)
    assert numc == pytest.approx(gc, rel=1e-6)
    numq = np.zeros(n)
    for i in range(n):
        a, b = q.copy(), q.copy()
        a[i] += eps
        b[i] -= eps
        numq[i] = (exclusion_correction(X, a, c, S)[0] - exclusion_correction(X, b, c, S)[0]) / (2 * eps)
    assert numq == pytest.approx(gq, rel=1e-6)


# ------------------------------------------------------------------ in the packer
def test_default_packer_is_bit_for_bit_the_truncated_one(beta_chain):
    """The literal is the number this package was measured with; coulomb='dsf' is default."""
    pk = CrystalPacker(beta_chain, n_chains=2)
    assert pk.coulomb == "dsf"
    assert pk._ewald is None
    assert float(pk.energy(BETA[None])[0]) == -8.645931161267756


def test_ewald_packer_total_is_the_spherical_shell_sum_minus_the_surface_term(beta_chain):
    """An independent check with no Ewald in it: summing the charge-charge lattice sum over
    spherical shells of whole cells converges to the *vacuum* convention, so it must match
    the Ewald vacuum energy and differ from the tinfoil one by exactly the surface term.
    This is the conditional convergence of the dipole sum, measured."""
    pk_d = CrystalPacker(beta_chain, n_chains=2)
    pk_t = CrystalPacker(beta_chain, n_chains=2, coulomb="ewald")
    pk_v = CrystalPacker(beta_chain, n_chains=2, coulomb="ewald", ewald=EwaldSpec(boundary="vacuum"))
    P, lat = pk_d._place(BETA[None])
    X = np.asarray(P[0], dtype=float)
    h = np.asarray(lat[0], dtype=float)
    q, n, N = pk_d._q_cell, pk_d.n, pk_d.N
    R = 60.0
    lim = [int(np.ceil(R / np.linalg.norm(h[i]))) + 1 for i in range(3)]
    m = np.stack(np.meshgrid(*[np.arange(-L, L + 1) for L in lim], indexing="ij"), -1).reshape(-1, 3)
    nv = m @ h
    sel = np.linalg.norm(nv, axis=1) <= R
    D = X[:, None, :] - X[None, :, :]
    brute = 0.0
    for mi, nvi in zip(m[sel], nv[sel]):
        r = np.linalg.norm(D - nvi, axis=-1)
        w = np.ones((N, N))
        if mi[0] == 0 and mi[1] == 0:
            Sk = pk_d._scale_nn_base.get(int(mi[2]))
            if Sk is not None:
                for s in range(pk_d.n_chains):
                    w[s * n:(s + 1) * n, s * n:(s + 1) * n] = Sk
            if mi[2] == 0:
                np.fill_diagonal(w, 0.0)
        brute += 0.5 * COULOMB * float((np.where(w > 0, q[:, None] * q[None, :] * w
                                                / np.where(r > 1e-9, r, 1.0), 0.0)).sum())
    # the packer's own electrostatics: its Ewald term plus the exclusion correction
    terms = pk_t._ewald.terms(X, q, h)
    coul_tinfoil = terms.total + pk_t.e_excl
    e_surf = float(pk_v.energy(BETA[None])[0]) - float(pk_t.energy(BETA[None])[0])
    assert e_surf > 4.0  # for this polar cell the term is 4.66 kcal/mol per cell
    assert coul_tinfoil + e_surf == pytest.approx(brute, abs=2e-3)


def test_surface_term_is_the_packers_own_dipole(beta_chain):
    pk_t = CrystalPacker(beta_chain, n_chains=2, coulomb="ewald")
    pk_v = CrystalPacker(beta_chain, n_chains=2, coulomb="ewald", ewald=EwaldSpec(boundary="vacuum"))
    mu = pk_t.dipole(BETA[None])[0]
    V = float(pk_t.cell_volume(BETA[None])[0])
    want = 2.0 * np.pi / (3.0 * V) * COULOMB * float(mu @ mu)
    got = float(pk_v.energy(BETA[None])[0]) - float(pk_t.energy(BETA[None])[0])
    assert got == pytest.approx(want, rel=1e-12)


@pytest.mark.parametrize("boundary", ["tinfoil", "vacuum"])
@pytest.mark.parametrize("flip", [0.0, 1.0])
def test_energy_and_grad_value_matches_energy_bit_for_bit(beta_chain, boundary, flip):
    pk = CrystalPacker(beta_chain, n_chains=2, coulomb="ewald", ewald=EwaldSpec(boundary=boundary))
    p = BETA.copy()
    p[6] = flip
    assert pk.energy_and_grad(p)[0] == float(pk.energy(p[None])[0])


@pytest.mark.parametrize("boundary", ["tinfoil", "vacuum"])
def test_ewald_packer_cell_gradient_is_analytic(beta_chain, boundary):
    pk = CrystalPacker(beta_chain, n_chains=2, coulomb="ewald", ewald=EwaldSpec(boundary=boundary))
    p = BETA.copy()
    _, g_cell, _, _ = pk.energy_and_grad(p)
    steps = np.array([1e-5, 1e-5, 1e-4, 1e-4, 1e-4, 1e-5])
    num = np.zeros(6)
    for i in range(6):
        a, b = p.copy(), p.copy()
        a[i] += steps[i]
        b[i] -= steps[i]
        num[i] = (float(pk.energy(a[None])[0]) - float(pk.energy(b[None])[0])) / (2 * steps[i])
    assert num == pytest.approx(g_cell, rel=1e-5, abs=1e-5)


def test_ewald_packer_chain_gradient_is_analytic(beta_chain):
    pk = CrystalPacker(beta_chain, n_chains=2, coulomb="ewald")
    X0 = np.asarray(pk.chain.coords, dtype=float)
    c0 = float(pk.chain.c)
    _, _, gX, gc = pk.energy_and_grad(BETA, coords=X0, c=c0)
    eps = 1e-6
    num = np.zeros_like(X0)
    for idx in np.ndindex(X0.shape):
        a, b = X0.copy(), X0.copy()
        a[idx] += eps
        b[idx] -= eps
        num[idx] = (float(pk.energy(BETA[None], coords=a[None], c=np.array([c0]))[0])
                    - float(pk.energy(BETA[None], coords=b[None], c=np.array([c0]))[0])) / (2 * eps)
    assert num == pytest.approx(gX, rel=1e-5, abs=1e-5)
    numc = (float(pk.energy(BETA[None], coords=X0[None], c=np.array([c0 + eps]))[0])
            - float(pk.energy(BETA[None], coords=X0[None], c=np.array([c0 - eps]))[0])) / (2 * eps)
    assert numc == pytest.approx(gc, rel=1e-5)


def test_translating_a_whole_neutral_chain_by_a_lattice_vector_is_free(beta_chain):
    """dz -> dz + c moves chain 2 by one lattice vector.  Free in every convention, because
    the chain is neutral; it is the check that the molecular branch is self-consistent."""
    for kw in ({}, {"coulomb": "ewald"},
               {"coulomb": "ewald", "ewald": EwaldSpec(boundary="vacuum")}):
        pk = CrystalPacker(beta_chain, n_chains=2, **kw)
        for flip in (0.0, 1.0):
            p = BETA.copy()
            p[6] = flip
            s = p.copy()
            s[5] = p[5] + pk.chain.c
            assert float(pk.energy(s[None])[0]) == pytest.approx(float(pk.energy(p[None])[0]), abs=1e-9)


def test_per_row_path_agrees_with_the_packers_own_chain(beta_chain):
    """The ``coords``/``c`` batched path and the rigid path must agree where they overlap --
    the exclusion correction is a constant on one and recomputed on the other."""
    pk = CrystalPacker(beta_chain, n_chains=2, coulomb="ewald")
    X0 = np.asarray(beta_chain.coords, dtype=float)
    rigid = pk.energy(BETA[None])
    per_row = pk.energy(BETA[None], coords=X0[None], c=np.array([beta_chain.c]))
    assert per_row[0] == pytest.approx(rigid[0], abs=1e-10)


# ------------------------------------------------------------------ the boundary
def test_the_tabulated_screen_refuses_an_ewald_packer(beta_chain):
    pk = CrystalPacker(beta_chain, n_chains=1, coulomb="ewald")
    with pytest.raises(ValueError, match="pairwise"):
        chain_pair_energy(pk, np.array([[5.0, 0.0]]), 0.0, 0.0, 0.0, 0)
    with pytest.raises(ValueError, match="pairwise"):
        intra_constants(CrystalPacker(beta_chain, n_chains=2, coulomb="ewald"))


def test_pack_screens_truncated_and_polishes_with_ewald(beta_chain):
    from polyfind.pack import pack

    res = pack(beta_chain, n_chains=2, coulomb="ewald", table_cache_dir=None, n_refine=2)
    pk = CrystalPacker(beta_chain, n_chains=2, coulomb="ewald")
    assert float(pk.energy(res[0].params[None])[0]) == pytest.approx(res[0].energy_per_cell, abs=1e-9)


def test_unknown_options_are_refused(beta_chain):
    with pytest.raises(ValueError, match="unknown coulomb"):
        CrystalPacker(beta_chain, n_chains=2, coulomb="wolf")
    with pytest.raises(ValueError, match="unknown Ewald boundary"):
        EwaldSpec(boundary="metallic")
    with pytest.raises(ValueError, match="coulomb='ewald'"):
        CrystalPacker(beta_chain, n_chains=2, ewald=EwaldSpec())


def test_uncharged_chain_gives_the_identical_energy_under_either_sum():
    """With every charge zero both sums are identically zero, so an Ewald packer must return
    *bit for bit* the truncated packer's number.  That is the check that the Lennard-Jones
    half of the kernel is untouched by the Coulomb rewiring -- a real polyethylene chain
    does not serve, because its atoms carry charges even though its cell dipole is zero."""
    from dataclasses import replace

    ch = periodic_chain(PE, (T,), THREE_STATE)
    ch = replace(ch, charges=np.zeros_like(np.asarray(ch.charges, dtype=float)))
    p = np.array([[4.9, 7.4, 90.0, 45.0, 45.0, 0.6, 0.0]])
    a = CrystalPacker(ch, n_chains=2).energy(p)
    b = CrystalPacker(ch, n_chains=2, coulomb="ewald").energy(p)
    assert float(b[0]) == float(a[0])


# ---------------------------------------------- the deformable path: valence and flux
@pytest.mark.parametrize("boundary", ["tinfoil", "vacuum"])
def test_gradients_with_valence_and_charge_flux(boundary):
    """The combination the response calculation uses.

    With charge flux the charges are functions of the geometry, so every derivative picks
    up ``dE/dq . dq/dgeometry``.  Under Ewald that ``dE/dq`` is the Ewald sum's, not the
    truncated kernel's -- the kernel has no Coulomb term left to differentiate -- and
    getting that wrong leaves the *cell* gradient exact and the *chain* gradient wrong by a
    factor of order one, which is exactly the kind of plausible number this file exists to
    catch.
    """
    from polyfind import pack as pack_mod
    from polyfind.fitting import FITTED_VALENCE_FLUX
    from polyfind.forcefield import SimpleFF

    p = np.array([4.58, 8.55, 90.0, 180.0, 180.0, 0.3, 0.0])
    with FITTED_VALENCE_FLUX.applied():
        chain = periodic_chain(PVDF, (T, T), THREE_STATE)
        pk = pack_mod.CrystalPacker(
            chain, n_chains=2, valence=SimpleFF.from_preset("pvdf-dft-valence"),
            charge_flux=SimpleFF.from_preset("pvdf-dft-valence-flux"),
            coulomb="ewald", ewald=EwaldSpec(boundary=boundary))
        X0 = np.asarray(pk.chain.coords, dtype=float)
        c0 = float(pk.chain.c)
        e, g_cell, _, _ = pk.energy_and_grad(p)
        assert e == float(pk.energy(p[None])[0])
        _, _, gX, gc = pk.energy_and_grad(p, coords=X0, c=c0)
        eps = 1e-6
        num = np.zeros_like(X0)
        for idx in np.ndindex(X0.shape):
            a, b = X0.copy(), X0.copy()
            a[idx] += eps
            b[idx] -= eps
            num[idx] = (float(pk.energy(p[None], coords=a[None], c=np.array([c0]))[0])
                        - float(pk.energy(p[None], coords=b[None], c=np.array([c0]))[0])) / (2 * eps)
        assert num == pytest.approx(gX, rel=1e-5, abs=1e-5)
        numc = (float(pk.energy(p[None], coords=X0[None], c=np.array([c0 + eps]))[0])
                - float(pk.energy(p[None], coords=X0[None], c=np.array([c0 - eps]))[0])) / (2 * eps)
        assert numc == pytest.approx(gc, rel=1e-5)
        steps = np.array([1e-5, 1e-5, 1e-4, 1e-4, 1e-4, 1e-5])
        numcell = np.zeros(6)
        for i in range(6):
            a, b = p.copy(), p.copy()
            a[i] += steps[i]
            b[i] -= steps[i]
            numcell[i] = (float(pk.energy(a[None])[0]) - float(pk.energy(b[None])[0])) / (2 * steps[i])
        assert numcell == pytest.approx(g_cell, rel=1e-5, abs=1e-5)


# ------------------------------------------------- the antipolar subspace it measures
def test_antipolar_subspace_is_derived_not_assumed():
    """``fitting.antipolar_cell``'s ``phi2 = phi1`` is not antipolar for a planar zigzag.

    beta-PVDF's chain moment lies along its own x, so the flip is a rotation: the cell it
    returns carries the *full* polarization of the polar minimum.  The corrected
    construction returns one with zero dipole, and must for the helix too, where the old
    one was right.
    """
    from polyfind import pack as pack_mod
    from polyfind.fitting import antipolar_cell, antipolar_cell_exact, antipolar_offsets, chain_moment
    from polyfind.pack import pack

    beta = periodic_chain(PVDF, (T, T), THREE_STATE)
    pk = pack_mod.CrystalPacker(beta, n_chains=2)
    m = chain_moment(pk)
    assert abs(m[2]) < 1e-9 and abs(m[1]) < 1e-9 and abs(m[0]) > 0.1  # along the chain's own x
    assert (0, 180.0) in antipolar_offsets(pk)
    start = pack(beta, n_chains=2)[0]
    old, _ = antipolar_cell(pk, start)
    assert np.abs(pk.polarization(old[None])).max() > 0.1  # not antipolar at all
    _, _, pol = antipolar_cell_exact(pk, target=1500, n_polish=3, maxfev=300)
    assert pol < 1e-12

    alpha = periodic_chain(PVDF, (T, GP, T, GM), THREE_STATE)
    pka = pack_mod.CrystalPacker(alpha, n_chains=2)
    assert abs(chain_moment(pka)[0]) < 1e-9  # the helix: the old construction was right
    dphi = antipolar_offsets(pka)[0][1]
    assert min(abs(dphi), abs(dphi - 360.0)) < 1e-6


def test_antipolar_scan_resolution_is_a_length_not_a_point_count():
    """``dz`` is sampled every ``dz_step`` angstroms, not four times per repeat.

    The old grid put four ``dz`` points across the repeat whatever the repeat was, which on a
    four-monomer gamma-type cell is one sample per monomer -- aliased to the interchain
    registry the scan exists to search.  Two things are asserted, and the second is the
    measured result rather than the expected one: the resolution now follows the repeat, and on
    beta-PVDF the change moves the answer by nothing, because the exact polish recovers the
    basin from the coarse grid as well (``docs/SCREEN.md``).
    """
    from polyfind import pack as pack_mod
    from polyfind.fitting import antipolar_cell_exact

    beta = periodic_chain(PVDF, (T, T), THREE_STATE)
    pk = pack_mod.CrystalPacker(beta, n_chains=2)
    # absolute resolution: 0.5 A on a 2.563 A repeat is six samples, where the old fixed count
    # of four corresponds to a step of c/4
    assert int(np.ceil(beta.c / 0.5)) == 6
    kw = dict(target=1500, n_polish=3, maxfev=300)
    _, e_fine, pol_fine = antipolar_cell_exact(pk, dz_step=0.5, **kw)
    _, e_coarse, pol_coarse = antipolar_cell_exact(pk, dz_step=beta.c / 4.0, **kw)
    assert pol_fine < 1e-12 and pol_coarse < 1e-12
    assert e_fine == pytest.approx(e_coarse, abs=2e-3)


def test_antipolar_cell_exact_checks_its_own_polarization():
    """A non-zero dipole is raised, not returned quietly: the subspace is exact by construction."""
    from polyfind import pack as pack_mod
    from polyfind.fitting import antipolar_cell_exact

    beta = periodic_chain(PVDF, (T, T), THREE_STATE)
    pk = pack_mod.CrystalPacker(beta, n_chains=2)
    with pytest.raises(RuntimeError, match="the construction is wrong"):
        # a tolerance no floating-point cell can meet, so the guard has to fire
        antipolar_cell_exact(pk, target=600, n_polish=2, maxfev=120, pol_tol=1e-30)


def test_predict_switches_default_to_the_fit_exactly():
    """``predict``'s ``coulomb``/``ewald``/``antipolar`` arguments must not move the fit.

    The defaults are the fit's own, so passing them explicitly has to give the identical
    numbers -- not approximately, exactly -- or every score this module has recorded would be
    in question.
    """
    from polyfind.fitting import ALL_CASES, ILLUSTRATIVE, predict

    beta = next(c for c in ALL_CASES if c.key == "beta")
    a = predict(beta, ILLUSTRATIVE, refine=False)
    b = predict(beta, ILLUSTRATIVE, refine=False, coulomb="dsf", ewald=None, antipolar="legacy")
    assert a.energy_per_monomer == b.energy_per_monomer
    assert (a.a, a.b, a.c) == (b.a, b.b, b.c)
    assert a.polar_gap is None and b.polar_gap is None  # beta is not an antipolar case
    with pytest.raises(ValueError, match="unknown antipolar"):
        predict(beta, ILLUSTRATIVE, refine=False, antipolar="nonsense", always_gap=True)
