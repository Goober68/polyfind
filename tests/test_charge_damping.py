"""Charge penetration in the Ewald sum: off is the point-charge sum bit for bit; on, the energy is
the point-charge sum minus the direct damped-term sum, the exclusions stay excluded, and the
gradients are the derivatives of the energy."""
import numpy as np

from polyfind import pack as P
from polyfind.ewald import COULOMB, Ewald, EwaldSpec, exclusion_correction
from polyfind.polymers import PVDF, THREE_STATE


def _packers(b):
    chain = P.periodic_chain(PVDF, [0, 0], THREE_STATE)
    pk0 = P.CrystalPacker(chain, n_chains=2, coulomb="ewald", lj_cutoff="force")
    pk1 = P.CrystalPacker(chain, n_chains=2, coulomb="ewald", lj_cutoff="force", ewald=EwaldSpec(charge_damping=b))
    return chain, pk0, pk1


def _g(r, b):
    return (1.0 + 0.5 * b * r) * np.exp(-b * r) / r


def test_off_is_bit_identical():
    chain, pk0, _ = _packers(2.0)
    pk_none = P.CrystalPacker(chain, n_chains=2, coulomb="ewald", lj_cutoff="force", ewald=EwaldSpec())
    p = np.array([4.6, 8.6, 90.0, 20.0, 200.0, 0.3, 0.0])
    assert pk_none._ewald.charge_damping is None
    e0, g0 = pk0.energy_and_grad(p[None])[:2]
    e1, g1 = pk_none.energy_and_grad(p[None])[:2]
    assert np.array_equal(np.asarray(e0), np.asarray(e1)) and np.array_equal(np.asarray(g0), np.asarray(g1))


def test_damped_energy_is_point_energy_minus_the_direct_term():
    """On a placed two-chain cell, E(damped) - E(point) = -1/2 pref sum q_i q_j g(r) over every
    non-excluded pair-image inside the cutoff (1-4 pairs at their scale)."""
    b = 2.0
    chain, pk0, pk1 = _packers(b)
    p = np.array([4.6, 8.6, 90.0, 20.0, 200.0, 0.3, 0.0])
    e0 = float(pk0.energy(p[None])[0])
    e1 = float(pk1.energy(p[None])[0])
    # brute force on the placed cell
    P_, L_ = pk0._place(p[None], None, np.array([chain.c]))
    Pn = np.asarray(P_[0], dtype=float); latn = np.asarray(L_[0], dtype=float)
    q = pk0._q_cell
    N = len(q); n = chain.n_atoms
    rc = pk0.rc
    # bonded-exclusion scales of one chain: (2K+1, n, n), k = -K..K; across chains everything is full
    direct = 0.0
    rng = range(-6, 7)  # lateral images: 6 a > rc
    for i in range(N):
        for j in range(N):
            for ka in rng:
                for kb in rng:
                    for kc in range(-4, 5):  # 4 c > rc
                        if i == j and ka == kb == kc == 0:
                            continue
                        d = Pn[i] - Pn[j] - ka * latn[0] - kb * latn[1] - kc * latn[2]
                        r = float(np.linalg.norm(d))
                        if r >= rc:
                            continue
                        s = 1.0
                        if (i < n) == (j < n) and ka == 0 and kb == 0:
                            # same chain, k-th z image: the packer's own bonded scale
                            S = pk0._scale_column_np(abs(kc) + 1)
                            s = float(S[kc + abs(kc) + 1, i % n, j % n])
                        direct += 0.5 * s * q[i] * q[j] * _g(r, b)
    assert np.isclose(e1 - e0, -COULOMB * direct, rtol=1e-6, atol=1e-8), (e1 - e0, -COULOMB * direct)


def test_damped_gradient_matches_finite_differences():
    chain, pk0, pk1 = _packers(1.5)
    p = np.array([4.6, 8.6, 90.0, 20.0, 200.0, 0.3, 0.0])
    e, g = pk1.energy_and_grad(p[None])[:2]
    g = np.asarray(g, dtype=float).reshape(-1)
    for k, h in ((0, 1e-5), (1, 1e-5), (2, 1e-4), (3, 1e-4), (4, 1e-4), (5, 1e-5)):
        pp, pm = p.copy(), p.copy()
        pp[k] += h; pm[k] -= h
        fd = (float(pk1.energy(pp[None])[0]) - float(pk1.energy(pm[None])[0])) / (2 * h)
        assert np.isclose(g[k], fd, rtol=2e-5, atol=2e-6), (k, g[k], fd)


def test_exclusion_correction_removes_the_damped_interaction():
    """With damping the correction subtracts (1 - s) q_i q_j [1/r - g(r)] for the excluded pairs."""
    b = 2.0
    chain, pk0, _ = _packers(b)
    X = np.asarray(chain.coords, dtype=float); q = np.asarray(chain.charges, dtype=float)
    K = 3
    S = pk0._scale_column_np(K)
    e0, gX0, gc0, gq0 = exclusion_correction(X, q, chain.c, S, grad=True, charge_grad=True)
    e1, gX1, gc1, gq1 = exclusion_correction(X, q, chain.c, S, grad=True, charge_grad=True, charge_damping=b)
    n = X.shape[0]
    direct = 0.0
    for k in range(-K, K + 1):
        for i in range(n):
            for j in range(n):
                s = S[k + K, i, j]
                if abs(1.0 - s) == 0.0:
                    continue
                d = X[i] - X[j] - np.array([0.0, 0.0, k * chain.c])
                r = float(np.linalg.norm(d))
                if r < 1e-6:
                    continue
                direct += 0.5 * (1.0 - s) * q[i] * q[j] * _g(r, b)
    # e = -1/2 pref sum (1-s) qq h with h = 1/r - g  ->  e1 - e0 = +1/2 pref sum (1-s) qq g
    assert np.isclose(e1 - e0, COULOMB * direct, rtol=1e-9)
    # gradient of the damped correction against finite differences of its own energy
    h = 1e-6
    for a in range(min(n, 3)):
        for c in range(3):
            Xp, Xm = X.copy(), X.copy(); Xp[a, c] += h; Xm[a, c] -= h
            fd = (exclusion_correction(Xp, q, chain.c, S, charge_damping=b)[0] - exclusion_correction(Xm, q, chain.c, S, charge_damping=b)[0]) / (2 * h)
            assert np.isclose(gX1[a, c], fd, rtol=1e-5, atol=1e-7), (a, c, gX1[a, c], fd)
    fd_c = (exclusion_correction(X, q, chain.c + h, S, charge_damping=b)[0] - exclusion_correction(X, q, chain.c - h, S, charge_damping=b)[0]) / (2 * h)
    assert np.isclose(gc1, fd_c, rtol=1e-5, atol=1e-7)
