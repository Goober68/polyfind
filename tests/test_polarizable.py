"""Induced point dipoles: the Ewald dipole sum, the self-consistent solve, and the packer.

Every check here is against something known independently of this package, because the
failure mode of a polarizable model is a plausible number: the Lorentz field of a simple
cubic dipole lattice, the Clausius-Mossotti dielectric constant of a cubic lattice of
polarizable points, independence of the Ewald splitting parameter, gradients against
central differences, the stationarity of the energy in the induced dipoles (which is what
makes the Hellmann-Feynman gradient exact), a centrosymmetric arrangement carrying no
induced moment, polyethylene staying non-piezoelectric, and the default kernel staying
bit-for-bit what it was.
"""
import numpy as np
import pytest

from polyfind import mechanics as M
from polyfind.ewald import Ewald, EwaldSpec, Thole, charge_dipole_exclusion
from polyfind.forcefield import COULOMB, SimpleFF
from polyfind.pack import EV_TO_KCAL, CrystalPacker, periodic_chain
from polyfind.polarizability import VDS98, Polarizable, clausius_mossotti
from polyfind.polymers import PE, PVDF, THREE_STATE

T = 0
BETA = np.array([4.58, 8.55, 90.0, 180.0, 180.0, 0.3, 0.0])


def _random_cell(seed=3, n=7):
    rng = np.random.default_rng(seed)
    h = np.array([[5.1, 0.0, 0.0], [0.9, 6.3, 0.0], [0.2, -0.4, 4.7]])
    X = rng.uniform(0, 1, (n, 3)) @ h
    q = rng.normal(size=n)
    q -= q.mean()
    p = 0.3 * rng.normal(size=(n, 3))
    alpha = rng.uniform(0.4, 1.4, n)
    s = (alpha[:, None] * alpha[None, :]) ** (1.0 / 6.0)
    return X, q, p, h, s


def _fd(base, put, eps=1e-6):
    out = np.zeros(base.shape)
    for idx in np.ndindex(base.shape):
        a, b = base.copy(), base.copy()
        a[idx] += eps
        b[idx] -= eps
        out[idx] = (put(a) - put(b)) / (2 * eps)
    return out


# ------------------------------------------------------------------ the dipole Ewald sum
def test_lorentz_field_of_a_simple_cubic_dipole_lattice():
    """Tinfoil: the field at a site from all its images is 4 pi p / 3V; vacuum: zero.

    The textbook result -- the Lorentz local field of a cubic lattice -- and the sharpest
    check that the reciprocal, real and self parts of the dipole sum add up to the right
    thing, since each is of order one and they nearly cancel.
    """
    a0 = 3.0
    h = np.eye(3) * a0
    for boundary, want in (("tinfoil", 4 * np.pi / (3 * a0 ** 3)), ("vacuum", 0.0)):
        ew = Ewald(12.0, EwaldSpec(boundary=boundary, accuracy=1e-12), coulomb=1.0)
        Tm = ew.terms(np.zeros((1, 3)), np.zeros(1), h, tensor=True).tensor[0, 0]
        assert Tm == pytest.approx(want * np.eye(3), abs=1e-10)


def test_clausius_mossotti_from_the_field_response():
    """A cubic lattice of one polarizable point, field applied under tinfoil: (eps-1)/(eps+2) = 4 pi alpha / 3V."""
    a0, alpha = 3.0, 2.0
    h = np.eye(3) * a0
    ew = Ewald(12.0, EwaldSpec(accuracy=1e-12), coulomb=1.0)
    Tm = ew.terms(np.zeros((1, 3)), np.zeros(1), h, tensor=True).tensor[0, 0]
    E0 = np.array([1e-3, 0.0, 0.0])
    p = np.linalg.solve(np.eye(3) / alpha - Tm, E0)
    eps = 1.0 + 4 * np.pi * p[0] / (a0 ** 3 * E0[0])
    assert eps == pytest.approx(clausius_mossotti(alpha, a0 ** 3), rel=1e-10)


@pytest.mark.parametrize("boundary", ["tinfoil", "vacuum"])
@pytest.mark.parametrize("thole", [None, Thole()])
def test_charge_dipole_energy_does_not_depend_on_the_splitting_parameter(boundary, thole):
    X, q, p, h, s = _random_cell()
    kw = dict(thole=thole, thole_scale=None if thole is None else s)
    vals = [Ewald(9.0, EwaldSpec(boundary=boundary, accuracy=1e-10, alpha=a)).terms(X, q, h, dipoles=p, **kw).total
            for a in (0.6, 0.75, 0.9)]
    assert max(vals) - min(vals) < 1e-6 * abs(vals[0])


@pytest.mark.parametrize("boundary", ["tinfoil", "vacuum"])
def test_charge_dipole_gradients_against_central_differences(boundary):
    X, q, p, h, s = _random_cell()
    ew = Ewald(9.0, EwaldSpec(boundary=boundary, accuracy=1e-10))
    kw = dict(thole=Thole(), thole_scale=s)
    t = ew.terms(X, q, h, dipoles=p, grad=True, charge_grad=True, dipole_grad=True, field=True, tensor=True, **kw)

    def E(XX=X, qq=q, hh=h, pp=p):
        return ew.terms(XX, qq, hh, dipoles=pp, **kw).total

    assert _fd(X, lambda v: E(XX=v)) == pytest.approx(t.grad_coords, rel=1e-6, abs=1e-6)
    assert _fd(h, lambda v: E(hh=v)) == pytest.approx(t.grad_lattice, rel=1e-6, abs=1e-6)
    assert _fd(q, lambda v: E(qq=v)) == pytest.approx(t.grad_charges, rel=1e-6, abs=1e-6)
    assert _fd(p, lambda v: E(pp=v)) == pytest.approx(t.grad_dipoles, rel=1e-6, abs=1e-6)
    # dE/dp is minus the total field, the field and the tensor being what the solve uses
    Tp = np.einsum("ijab,jb->ia", t.tensor, p)
    assert t.grad_dipoles == pytest.approx(-(t.field + Tp), abs=1e-10)
    assert t.tensor == pytest.approx(t.tensor.transpose(1, 0, 3, 2), abs=1e-10)


def test_charge_only_path_is_untouched():
    """Asking for nothing dipole-related gives the same five terms, expression for expression."""
    X, q, _, h, _ = _random_cell()
    ew = Ewald(9.0, EwaldSpec(accuracy=1e-10))
    a = ew.terms(X, q, h, grad=True, charge_grad=True)
    b = ew.terms(X, q, h, grad=True, charge_grad=True, dipoles=np.zeros_like(X))
    assert (a.real, a.recip, a.self_, a.surface) == (b.real, b.recip, b.self_, b.surface)
    assert a.grad_coords == pytest.approx(b.grad_coords, abs=1e-12)
    assert a.grad_lattice == pytest.approx(b.grad_lattice, abs=1e-12)
    assert a.grad_charges == pytest.approx(b.grad_charges, abs=1e-12)


def test_charge_dipole_exclusion_gradients():
    rng = np.random.default_rng(11)
    n, K, c = 5, 2, 2.6
    X = rng.normal(size=(n, 3)) * 1.2
    S = rng.uniform(0, 1, (2 * K + 1, n, n))
    S[S > 0.6] = 1.0
    q, p = rng.normal(size=n), 0.2 * rng.normal(size=(n, 3))
    e, gX, gc, gq, F = charge_dipole_exclusion(X, q, c, S, dipoles=p, grad=True, charge_grad=True)

    def E(XX=X, qq=q, cc=c, pp=p):
        return charge_dipole_exclusion(XX, qq, cc, S, dipoles=pp)[0]

    assert _fd(X, lambda v: E(XX=v)) == pytest.approx(gX, rel=1e-6, abs=1e-7)
    assert _fd(q, lambda v: E(qq=v)) == pytest.approx(gq, rel=1e-6, abs=1e-7)
    assert (E(cc=c + 1e-6) - E(cc=c - 1e-6)) / 2e-6 == pytest.approx(gc, rel=1e-6)
    assert _fd(p, lambda v: E(pp=v)) == pytest.approx(-F, rel=1e-6, abs=1e-7)
    # the field correction does not depend on the dipoles
    assert charge_dipole_exclusion(X, q, c, S)[4] == pytest.approx(F, abs=1e-14)


def test_thole_damping_forms_go_to_one_and_zero():
    for form in ("exp", "cubic"):
        th = Thole(2.1304 if form == "exp" else 0.39, form)
        f3, f5, d3, d5 = th.damping(np.array([1e-6, 20.0]))
        assert f3[0] == pytest.approx(1.0, abs=1e-5) and f5[0] == pytest.approx(1.0, abs=1e-9)
        assert abs(f3[1]) < 1e-6 and abs(f5[1]) < 1e-6
        u = 1.3
        h = 1e-6
        f3p, f5p, _, _ = th.damping(np.array([u + h]))
        f3m, f5m, _, _ = th.damping(np.array([u - h]))
        _, _, d3, d5 = th.damping(np.array([u]))
        assert (f3p - f3m) / (2 * h) == pytest.approx(d3, rel=1e-6)
        assert (f5p - f5m) / (2 * h) == pytest.approx(d5, rel=1e-6)


# ------------------------------------------------------------------------ the packer
@pytest.fixture(scope="module")
def beta_chain():
    return periodic_chain(PVDF, (T, T), THREE_STATE)


def test_polarizable_needs_ewald(beta_chain):
    with pytest.raises(ValueError, match="ewald"):
        CrystalPacker(beta_chain, polarizable=Polarizable())


def test_default_packer_is_bit_for_bit_unchanged(beta_chain):
    """``polarizable=None`` (the default) must not have moved a single bit."""
    p = BETA.copy()
    a = CrystalPacker(beta_chain)
    b = CrystalPacker(beta_chain, coulomb="ewald", ewald=EwaldSpec())
    assert float(a.energy(p[None])[0]) == a.energy_and_grad(p)[0]
    assert float(b.energy(p[None])[0]) == b.energy_and_grad(p)[0]
    assert a.polarizable is None and b.polarizable is None


@pytest.mark.parametrize("flip", [0.0, 1.0])
@pytest.mark.parametrize("field", [None, (0.01, -0.02, 0.005)])
def test_polarizable_value_and_gradients(beta_chain, flip, field):
    """``energy`` and ``energy_and_grad`` agree bit for bit; the gradient is the central difference.

    Valence terms, charge flux, Ewald and induced dipoles all on: the combination the response
    calculation uses.  The gradient is taken at fixed induced dipoles and is only the total
    derivative because the solve is stationary, so this is the test of that argument.
    """
    from polyfind import pack as pack_mod
    from polyfind.fitting import FITTED_VALENCE_FLUX

    p = BETA.copy()
    p[6] = flip
    with FITTED_VALENCE_FLUX.applied():
        chain = periodic_chain(PVDF, (T, T), THREE_STATE)
        pk = pack_mod.CrystalPacker(
            chain, n_chains=2, valence=SimpleFF.from_preset("pvdf-dft-valence"),
            charge_flux=SimpleFF.from_preset("pvdf-dft-valence-flux"),
            coulomb="ewald", ewald=EwaldSpec(), polarizable=Polarizable(), field=field)
        e, g_cell, _, _ = pk.energy_and_grad(p)
        assert e == float(pk.energy(p[None])[0])
        X0 = np.asarray(pk.chain.coords, dtype=float)
        c0 = float(pk.chain.c)
        _, _, gX, gc = pk.energy_and_grad(p, coords=X0, c=c0)
        eps = 1e-5

        def E(XX=X0, cc=c0):
            return float(pk.energy(p[None], coords=XX[None], c=np.array([cc]))[0])

        assert _fd(X0, lambda v: E(XX=v), eps) == pytest.approx(gX, rel=1e-5, abs=1e-5)
        assert (E(cc=c0 + eps) - E(cc=c0 - eps)) / (2 * eps) == pytest.approx(gc, rel=1e-5)
        steps = np.array([1e-5, 1e-5, 1e-4, 1e-4, 1e-4, 1e-5])
        num = np.zeros(6)
        for i in range(6):
            a, b = p.copy(), p.copy()
            a[i] += steps[i]
            b[i] -= steps[i]
            num[i] = (float(pk.energy(a[None])[0]) - float(pk.energy(b[None])[0])) / (2 * steps[i])
        assert num == pytest.approx(g_cell, rel=1e-5, abs=1e-5)


@pytest.mark.parametrize("flip", [0.0, 1.0])
def test_flipped_chain_uses_the_reversed_exclusion_stack(flip):
    """A flipped chain 2 is mirrored z -> -z, so its image k is the chain frame's image -k.

    On a helix (no two-fold axis along x, so the flip is not a symmetry) the excluded-charge
    field of the placed chain must equal the chain frame's, transformed by the same isometry;
    with the unreversed stack it is off by half its size.
    """
    T, GP, GM = 0, 1, 2
    chain = periodic_chain(PVDF, (T, GP, T, GM), THREE_STATE)
    pk = CrystalPacker(chain, n_chains=2, coulomb="ewald", ewald=EwaldSpec(), polarizable=Polarizable())
    n = pk.n
    p = np.array([5.0, 9.6, 90.0, 20.0, 130.0, 1.1, flip])
    P, lat = pk._place(p[None])
    Pn, latn = np.asarray(P[0], dtype=float), np.asarray(lat[0], dtype=float)
    cz = float(latn[2, 2])
    Sk = pk._scale_column_np(pk.K)
    right = Sk[::-1] if flip > 0.5 else Sk
    Fa = charge_dipole_exclusion(Pn[n:], pk._q_cell[n:], cz, right, pref=1.0)[4]
    Fc = charge_dipole_exclusion(np.asarray(chain.coords, dtype=float), pk._q_cell[:n], cz, Sk, pref=1.0)[4]
    if flip > 0.5:
        Fc = Fc * np.array([1.0, -1.0, -1.0])
    ang = np.deg2rad(p[4])
    R = np.array([[np.cos(ang), -np.sin(ang), 0.0], [np.sin(ang), np.cos(ang), 0.0], [0.0, 0.0, 1.0]])
    Fb = Fc @ R.T
    assert Fa == pytest.approx(Fb, abs=1e-12)
    wrong = charge_dipole_exclusion(Pn[n:], pk._q_cell[n:], cz, Sk if flip > 0.5 else Sk[::-1], pref=1.0)[4]
    assert np.abs(wrong - Fb).max() > 0.1 * np.abs(Fb).max()
    # and translating chain 2 by one repeat is free, polarization energy included
    p2 = p.copy()
    p2[5] += cz
    assert float(pk.energy(p2[None])[0]) == pytest.approx(float(pk.energy(p[None])[0]), abs=1e-10)


def test_energy_is_stationary_in_the_induced_dipoles(beta_chain):
    """U(p) = sum p^2/2a - p.E0 - (1/2) p.T.p is minimised by the solve, and U* = -(1/2) p.E0.

    The first-order change along a random direction must vanish, the second-order one must
    be positive (it is a minimum), and the closed form must equal the functional.
    """
    pk = CrystalPacker(beta_chain, n_chains=2, coulomb="ewald", ewald=EwaldSpec(),
                       polarizable=Polarizable(), field=(0.01, 0.0, 0.0))
    P, lat = pk._place(BETA[None])
    Pn, latn = np.asarray(P[0], dtype=float), np.asarray(lat[0], dtype=float)
    e_es, p, E0, _ = pk._polarize(Pn, pk._q_cell, latn, 0.0)
    n = pk.n

    def U(pv):
        t = pk._ewald.terms(Pn, pk._q_cell, latn, dipoles=pv, thole=pk._thole, thole_scale=pk._thole_scale)
        ex = sum(charge_dipole_exclusion(Pn[s * n:(s + 1) * n], pk._q_cell[s * n:(s + 1) * n], latn[2, 2],
                                         pk._scale_column_np(pk.K), dipoles=pv[s * n:(s + 1) * n])[0] for s in range(2))
        self_e = float((COULOMB / pk._alpha_cell * (pv * pv).sum(axis=1) / 2).sum())
        return t.total + ex + self_e - EV_TO_KCAL * float((pv * np.asarray(pk.field)).sum())

    U0 = U(p)
    assert U0 == pytest.approx(e_es, abs=1e-10)
    assert U0 == pytest.approx(pk._ewald.terms(Pn, pk._q_cell, latn).total - 0.5 * float((p * E0).sum()), abs=1e-10)
    rng = np.random.default_rng(1)
    d = rng.normal(size=p.shape)
    d /= np.linalg.norm(d)
    for h in (1e-2, 1e-3):
        first = (U(p + h * d) - U(p - h * d)) / (2 * h)
        second = (U(p + h * d) + U(p - h * d) - 2 * U0) / h ** 2
        assert abs(first) < 1e-8 * second
        assert second > 0


def test_centrosymmetric_sites_carry_no_induced_dipole():
    """Every site of rock salt is an inversion centre: the permanent field vanishes there, so p = 0."""
    from polyfind.ewald import Ewald

    a = 2.8
    idx = np.array([[i, j, k] for i in (0, 1) for j in (0, 1) for k in (0, 1)], dtype=float)
    X = idx * a
    q = np.array([1.0 if (i + j + k) % 2 == 0 else -1.0 for i, j, k in idx.astype(int)])
    ew = Ewald(12.0, EwaldSpec(accuracy=1e-12))
    t = ew.terms(X, q, np.eye(3) * 2 * a, field=True, tensor=True)
    assert np.abs(t.field).max() < 1e-10
    alpha = np.where(q > 0, 0.3, 1.5)
    A = -t.tensor.transpose(0, 2, 1, 3).reshape(24, 24)
    A[np.arange(24), np.arange(24)] += np.repeat(COULOMB / alpha, 3)
    p = np.linalg.solve(A, t.field.ravel())
    assert np.abs(p).max() < 1e-12


def test_polyethylene_induced_moment_vanishes_and_it_stays_non_piezoelectric():
    """PE's cell is centrosymmetric: the induced dipoles sum to zero at zero field and e, d stay at machine zero."""
    chain = periodic_chain(PE, (T,), THREE_STATE)
    pk = CrystalPacker(chain, n_chains=2, coulomb="ewald", ewald=EwaldSpec(), polarizable=Polarizable())
    for params in (np.array([7.4, 4.95, 90.0, 48.0, -48.0, 1.27, 0.0]), np.array([7.4, 4.95, 90.0, 48.0, 132.0, 1.27, 1.0])):
        p = pk.induced_dipoles(params[None])[0]
        assert np.abs(p).max() > 1e-4  # the atoms are polarised ...
        assert np.linalg.norm(p.sum(axis=0)) < 1e-12  # ... and the cell is not
        assert np.linalg.norm(pk.dipole(params[None])[0]) < 1e-12
    ref, _ = M.refined_reference(PE, [T], label="PE", coulomb="ewald", ewald=EwaldSpec(), polarizable=Polarizable())
    ref = M.relax_reference(ref)
    el = M.elastic_constants(ref)
    pz = M.piezoelectric(ref, el, converse=True)
    assert np.abs(pz.e).max() < 1e-9
    assert np.abs(pz.d_from_e).max() < 1e-6 and np.abs(pz.d_direct).max() < 1e-6
    # ... while its dielectric constant is no longer 1
    de = M.dielectric_tensor(ref, relax=False)
    assert np.all(np.diag(de.clamped) > 1.5) and np.all(np.diag(de.clamped) < 3.5)


def test_dielectric_tensor_is_one_without_polarizability_and_physical_with_it(beta_chain):
    ref = M.reference_from_chain(beta_chain, BETA, coulomb="ewald", ewald=EwaldSpec())
    assert M.dielectric_tensor(ref, relax=False).clamped == pytest.approx(np.eye(3), abs=1e-9)
    ref = M.reference_from_chain(beta_chain, BETA, coulomb="ewald", ewald=EwaldSpec(), polarizable=Polarizable())
    de = M.dielectric_tensor(ref, relax=False)
    d = np.diag(de.clamped)
    # beta-PVDF's electronic dielectric tensor from periodic DFPT is diag(2.24, 2.25, 2.60) with the
    # chain axis highest; the literature polarizabilities land within 0.15 of each component and
    # reproduce the ordering (docs/ELECTROMECHANICS.md section 5.8)
    assert np.abs(de.clamped - np.diag(d)).max() < 1e-6
    assert 1.9 < d[0] < 2.5 and 1.9 < d[1] < 2.5 and 2.3 < d[2] < 2.9
    assert d[2] > d[0] and d[2] > d[1]
    # linear: the same tensor at half the field
    assert M.dielectric_tensor(ref, field=1e-3, relax=False).clamped == pytest.approx(de.clamped, abs=1e-6)


def test_table_is_van_duijnen_swart_in_cubic_angstrom():
    assert VDS98["H"] == pytest.approx(0.4138, abs=1e-4)
    assert VDS98["C"] == pytest.approx(1.2886, abs=1e-4)
    assert VDS98["F"] == pytest.approx(0.3879, abs=1e-4)
    assert Polarizable().per_atom(["C", "H", "F"]) == pytest.approx([VDS98["C"], VDS98["H"], VDS98["F"]])
    with pytest.raises(KeyError):
        Polarizable().per_atom(["Xx"])
    with pytest.raises(ValueError):
        Polarizable(table="nope")
