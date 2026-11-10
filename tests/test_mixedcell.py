"""polyfind.mixedcell: a cell whose chains differ, held to the packer by known answers."""
import numpy as np
import pytest

from polyfind import mixedcell as MC
from polyfind import pack as pack_mod
from polyfind.ewald import EwaldSpec
from polyfind.fitting import PRESET_PARAMETERS
from polyfind.linegroup import repeat_chains, repeat_chains_and_jacobian
from polyfind.polarizability import Polarizable
from polyfind.polymers import PVDF, THREE_STATE, pvdf_belt

T, GP, GM = 0, 1, 2
BELT4 = "pvdf-dft-valence-belt4"
EWALD = dict(lj_cutoff="force", coulomb="ewald", ewald=EwaldSpec(), polarizable=Polarizable())


def _chain(polymer, seq):
    tors = np.array([THREE_STATE.angles[s] for s in seq], dtype=float)
    ang0 = np.array([polymer.backbone[k].backbone_angle for k in range(polymer.bonds_per_repeat)], dtype=float)
    name = "".join(THREE_STATE.names[s] for s in seq)
    return repeat_chains(polymer, name, tors[None], ang0)[0], name, tors, ang0


def _model(polymer, seq, params, nx=1, ny=1, pen_rot=0.0, k_c=0.0, **kw):
    ch, name, tors, ang0 = _chain(polymer, seq)
    nc = 2 * nx * ny
    pk = pack_mod.CrystalPacker(ch, n_chains=nc, **kw)
    st, flips = MC.supercell_state(params, ch.c, tors, ang0, nx, ny)
    return MC.MixedModel(polymer, name, pk, flips, [ch] * nc, angle_stiffness=105.0, pen_rot=pen_rot, k_c=k_c), st.vector(), ch


CASES = [
    (PVDF, [T, GP, T, GM], [5.0, 9.6, 90.0, 30.0, 200.0, 1.1, 1.0], None),
    (PVDF, [T, T], [8.6, 4.9, 93.0, 40.0, 40.0, 0.3, 0.0], None),
    ("belt", [T, T], [4.75, 8.48, 90.0, 180.0, 180.0, 1.3, 1.0], BELT4),
    ("belt", [T, GP, T, GM], [5.0, 9.5, 91.0, 30.0, 215.0, 2.2, 1.0], BELT4),
]


def _run(case, body):
    polymer, seq, params, preset = case
    if preset is None:
        return body(polymer, seq, np.asarray(params, dtype=float), {})
    with PRESET_PARAMETERS[preset].applied():
        return body(pvdf_belt(0.7, "pbe"), seq, np.asarray(params, dtype=float), EWALD)


@pytest.mark.parametrize("case", CASES)
def test_two_chain_mixed_cell_is_the_packer(case):
    """Energy and every gradient of the packer's own cell, to rounding."""
    def body(polymer, seq, params, kw):
        m, x, ch = _model(polymer, seq, params, **kw)
        E, g_cell, gX, g_c = m.pk.energy_and_grad(params)
        parts: dict = {}
        _, g = m.value_and_grad(x, parts)
        assert parts["physical"] == pytest.approx(E, rel=1e-11, abs=1e-9)
        assert parts["bend"] == 0.0
        np.testing.assert_allclose(g[:3], g_cell[:3], atol=1e-9)
        mine = [g[m.index("phi", 0)[0]], g[m.index("phi", 1)[0]], g[m.index("z", 1)[0]]]
        np.testing.assert_allclose(mine, g_cell[3:6], atol=1e-9)
        assert g[3] == pytest.approx(g_c, abs=1e-8)
        _, dX, dc, _, _ = repeat_chains_and_jacobian(polymer, m.name, ch.dihedrals[None], m.angles0, ch)
        folded = np.einsum("ni,vni->v", gX, dX[0]) + g_c * dc[0]
        folded[: m.P] += m.pk.torsion_energy_grad(ch.dihedrals) * 2
        np.testing.assert_allclose(g[m.index("internal", 0)] + g[m.index("internal", 1)] + g[3] * dc[0], folded, atol=1e-8)
    _run(case, body)


@pytest.mark.parametrize("case", CASES)
@pytest.mark.parametrize("tiling", [(2, 1), (1, 2), (2, 2)])
def test_a_tiling_has_the_same_energy_per_chain(case, tiling):
    def body(polymer, seq, params, kw):
        m2, x2, _ = _model(polymer, seq, params, **kw)
        e2 = m2.physical_energy(x2)
        m, x, _ = _model(polymer, seq, params, nx=tiling[0], ny=tiling[1], **kw)
        assert m.physical_energy(x) * 2 / m.nc == pytest.approx(e2, rel=1e-10, abs=1e-8)
    _run(case, body)


@pytest.mark.parametrize("case,tiling", [(CASES[0], (2, 1)), (CASES[3], (2, 1)), (CASES[2], (2, 2))])
def test_gradient_of_a_cell_of_different_chains_matches_central_differences(case, tiling):
    def body(polymer, seq, params, kw):
        m, x, _ = _model(polymer, seq, params, nx=tiling[0], ny=tiling[1], pen_rot=5.0, k_c=50.0, **kw)
        rng = np.random.default_rng(1)
        x[m.index("tors")] += rng.uniform(-6, 6, m.nc * m.P)
        x[m.index("angs")] += rng.uniform(-1.5, 1.5, m.nc * m.R)
        x[m.index("xy")] += rng.uniform(-0.02, 0.02, 2 * m.nc)
        x[m.index("phi")] += rng.uniform(-8, 8, m.nc)
        x[m.index("z")] += rng.uniform(-0.3, 0.3, m.nc)
        x[:4] += np.array([0.07, -0.05, 1.3, 0.04])
        parts: dict = {}
        _, g = m.value_and_grad(x, parts)
        assert parts["penalty_rotation"] > 1.0 and parts["penalty_c"] > 0.0  # both penalties are exercised
        for i in np.concatenate([np.arange(4), rng.choice(np.arange(4, m.n_var), size=12, replace=False)]):
            h = 1e-5 * max(1.0, abs(x[i])) if i >= 4 else 1e-5
            xp_, xm_ = x.copy(), x.copy()
            xp_[i] += h
            xm_[i] -= h
            fd = (m.value_and_grad(xp_)[0] - m.value_and_grad(xm_)[0]) / (2 * h)
            assert fd == pytest.approx(g[i], rel=2e-5, abs=2e-5)
    _run(case, body)


def test_state_and_index_roundtrip():
    m, x, _ = _model(PVDF, [T, GP, T, GM], [5.0, 9.6, 90.0, 30.0, 200.0, 1.1, 1.0], nx=2, ny=1)
    st = m.state(x)
    np.testing.assert_array_equal(st.vector(), x)
    assert st.lattice.tolist() == [10.0, 9.6, 90.0, pytest.approx(x[3])]
    assert m.n_var == 4 + 4 * (4 + 4 + 2)
    every = np.concatenate([m.index("lattice"), m.index("c"), m.index("place"), m.index("internal")])
    assert sorted(every.tolist()) == list(range(m.n_var))
    np.testing.assert_array_equal(m.index("tors", 2), 4 + 16 + 2 * 6 + np.arange(4))
    assert m.flips.tolist() == [False, True, False, True]


def test_closure_residual_jacobian_and_projection():
    m, x, _ = _model(PVDF, [T, GP, T, GM], [5.0, 9.6, 90.0, 30.0, 200.0, 1.1, 1.0], nx=2, ny=1)
    assert max(m.closure_residual(x)) < 1e-9  # the built chain is periodic
    rng = np.random.default_rng(2)
    x[m.index("tors")] += rng.uniform(-6, 6, m.nc * m.P)
    x[m.index("angs")] += rng.uniform(-1.5, 1.5, m.nc * m.R)
    x[3] += 0.05
    h, J = m.closure(x)
    assert h.shape == (4 * m.nc,) and J.shape == (4 * m.nc, m.n_var) and np.abs(h).max() > 0.1
    for i in np.concatenate([[3], rng.choice(m.index("internal"), size=8, replace=False)]):
        xp_, xm_ = x.copy(), x.copy()
        xp_[i] += 1e-5
        xm_[i] -= 1e-5
        np.testing.assert_allclose((m.closure(xp_)[0] - m.closure(xm_)[0]) / 2e-5, J[:, i], atol=2e-6)
    h1, J1 = m.closure(x, [2])  # one chain's four residuals are that chain's rows
    np.testing.assert_array_equal(h1, h[8:12])
    np.testing.assert_array_equal(J1, J[8:12])
    free = m.index("internal")
    y = m.project(x, free)
    assert max(m.closure_residual(y)) < 1e-8
    np.testing.assert_array_equal(y[:4 + 4 * m.nc], x[:4 + 4 * m.nc])  # only the free variables moved
    assert np.abs(y - x).max() < 6.0


def test_constrained_relaxation_lowers_the_energy_on_the_closure_surface():
    m, x, _ = _model(PVDF, [T, GP, T, GM], [5.0, 9.6, 90.0, 30.0, 200.0, 1.1, 1.0])
    e0 = m.physical_energy(x)
    x1, info = m.relax(x, maxiter=60)
    assert info["closed_chains"] == [0, 1]
    assert m.physical_energy(x1) < e0 - 1.0
    assert info["rotation_error_max"] < 1e-6 and info["c_mismatch_max"] < 1e-7  # periodic, exactly
    held, info = m.relax(x, free=m.index("place"), maxiter=10)  # only the placements move: nothing to close
    assert info["closed_chains"] == []
    np.testing.assert_array_equal(held[m.index("internal")], x[m.index("internal")])
    np.testing.assert_array_equal(held[:4], x[:4])
    assert m.physical_energy(held) < e0


def test_penalty_terms_are_off_by_default_and_available():
    m, x, _ = _model(PVDF, [T, GP, T, GM], [5.0, 9.6, 90.0, 30.0, 200.0, 1.1, 1.0])
    x[m.index("tors")] += 3.0
    parts: dict = {}
    total, _ = m.value_and_grad(x, parts)
    assert parts["penalty_rotation"] == 0.0 and parts["penalty_c"] == 0.0 and total == parts["physical"]
    assert parts["rotation_error"].max() > 0.1
    m.pen_rot, m.k_c = 5.0, 50.0
    parts2: dict = {}
    total2, _ = m.value_and_grad(x, parts2)
    assert total2 > parts2["physical"] == parts["physical"]


def test_tilt_relaxation_matches_the_packers_own():
    from polyfind import mechanics as M
    from polyfind import pendants as PD
    from polyfind.refine import resolve_angle_stiffness

    with PRESET_PARAMETERS[BELT4].applied():
        pb = pvdf_belt(0.7, "pbe")
        kw = {"pack_kw": {"table_cache_dir": None}, **EWALD}
        ref, rr = M.refined_reference(pb, [T, GP, T, GM], label="alpha", refine_kw=dict(EWALD), **kw)
        ref = M.relax_reference(ref)
        ch = ref.packer.chain
        st, flips = MC.supercell_state(ref.params, ref.c, rr.torsions, rr.angles, 1, 1)
        m = MC.MixedModel(pb, "alpha", ref.packer, flips, [ch, ch], angle_stiffness=resolve_angle_stiffness(None, None))
        x = st.vector()
        e_ref = float(ref.packer.energy(np.asarray(ref.params)[None])[0]) + rr.angle_energy * 2 * ch.n_monomers
        assert m.physical_energy(x) == pytest.approx(e_ref, abs=1e-8)  # the refinement's own number
        dE, tilts, groups = MC.tilt_relaxation(m, x)
        theirs = PD.tilt_relaxation(ref, pb)
        assert dE == pytest.approx(theirs.dE_per_cell, abs=1e-7) and dE < -0.5
        assert len(groups) == len(theirs.groups) and tilts.shape == (len(groups), 3)


def _muller_brown(X):
    A = np.array([-200.0, -100.0, -170.0, 15.0])
    a, b, c = np.array([-1.0, -1.0, -6.5, 0.7]), np.array([0.0, 0.0, 11.0, 0.6]), np.array([-10.0, -10.0, -6.5, 0.7])
    x0, y0 = np.array([1.0, 0.0, -0.5, -1.0]), np.array([0.0, 0.5, 1.5, 1.0])
    X = np.atleast_2d(X)
    dx, dy = X[:, :1] - x0[None], X[:, 1:2] - y0[None]
    t = A[None] * np.exp(a[None] * dx * dx + b[None] * dx * dy + c[None] * dy * dy)
    return t.sum(axis=1), np.stack([(t * (2 * a[None] * dx + b[None] * dy)).sum(axis=1), (t * (b[None] * dx + 2 * c[None] * dy)).sum(axis=1)], axis=1)


def test_neb_finds_the_muller_brown_saddle():
    a, b = np.array([-0.5582236, 1.4417258]), np.array([-0.0500108, 0.4666941])  # two of its minima
    images = np.linspace(a, b, 11)
    res = MC.neb(_muller_brown, images, k_spring=200.0, climb_after=40, fmax=0.05, maxiter=2000, max_step=0.02)
    assert res.converged and res.climbing > 0
    assert res.energies.max() == pytest.approx(-40.6648, abs=0.01)  # the saddle between them
    np.testing.assert_allclose(res.images[res.climbing], [-0.8220, 0.6243], atol=2e-3)
    np.testing.assert_array_equal(res.images[[0, -1]], images[[0, -1]])


def test_neb_on_a_constraint_surface_is_the_neb_of_the_reduced_problem():
    """Three variables, the third tied to the others by z = 0.3 x y, and an energy that depends on
    all three: the band on that surface has the saddle of the two-variable function it reduces to."""
    from scipy.optimize import minimize

    slope = 25.0

    def reduced(X):
        X = np.atleast_2d(X)
        e, g = _muller_brown(X)
        return e + slope * 0.3 * X[:, 0] * X[:, 1], g + slope * 0.3 * X[:, ::-1]

    def embedded(X):
        X = np.atleast_2d(X)
        e, g = _muller_brown(X[:, :2])
        return e + slope * X[:, 2], np.concatenate([g, np.full((len(X), 1), slope)], axis=1)

    def surface(x):
        return np.array([x[2] - 0.3 * x[0] * x[1]]), np.array([[-0.3 * x[1], -0.3 * x[0], 1.0]])

    ends = [minimize(lambda p: float(reduced(p)[0][0]), p0, jac=lambda p: reduced(p)[1][0], method="BFGS", options={"gtol": 1e-10}).x
            for p0 in ([-0.5582236, 1.4417258], [-0.0500108, 0.4666941])]
    flat = np.linspace(ends[0], ends[1], 11)
    ref = MC.neb(reduced, flat, k_spring=200.0, climb_after=40, fmax=0.05, maxiter=3000, max_step=0.02)
    assert ref.converged and ref.climbing > 0
    lifted = np.concatenate([flat, 0.3 * flat[:, :1] * flat[:, 1:2]], axis=1)
    res = MC.neb(embedded, lifted, k_spring=200.0, climb_after=40, fmax=0.05, maxiter=3000, max_step=0.02, constraint=surface)
    assert res.converged and res.climbing > 0
    z = res.images[:, 2] - 0.3 * res.images[:, 0] * res.images[:, 1]
    assert np.abs(z).max() < 1e-7  # every image is on the surface
    assert res.energies.max() == pytest.approx(ref.energies.max(), abs=0.02)
    np.testing.assert_allclose(res.images[res.climbing, :2], ref.images[ref.climbing], atol=5e-3)
    assert abs(ref.energies.max() - (-40.6648)) > 1.0  # and it is not the unconstrained surface's saddle


def test_tiling_a_mixed_state_keeps_the_energy_per_chain():
    """A two-chain cell whose chains differ, tiled: the same energy per chain, and the state
    vector round trip."""
    def body(polymer, seq, params, kw):
        m, x, ch = _model(polymer, seq, params, **kw)
        rng = np.random.default_rng(5)
        x = x.copy()
        x[m.index("tors")] += rng.uniform(-4, 4, m.nc * m.P)
        x[m.index("angs")] += rng.uniform(-1, 1, m.nc * m.R)
        x[m.index("place")] += rng.uniform(-0.02, 0.02, 4 * m.nc)
        x = m.project(x, m.index("internal"))
        st = MC.MixedState.from_vector(x, m.nc, m.P)
        np.testing.assert_array_equal(st.vector(), x)
        e2 = m.physical_energy(x)
        for nx, ny in ((2, 1), (1, 2)):
            big, flips = MC.tile_state(st, m.flips, nx, ny)
            pk = pack_mod.CrystalPacker(ch, n_chains=2 * nx * ny, **kw)
            mm = MC.MixedModel(polymer, m.name, pk, flips, [ch] * (2 * nx * ny), angle_stiffness=105.0)
            assert mm.physical_energy(big.vector()) == pytest.approx(e2 * nx * ny, rel=1e-11, abs=1e-8)
    _run(CASES[3], body)


def test_descent_stays_on_the_closure_surface_and_in_step():
    """The damped pre-relaxation: energy down, every chain still periodic, no variable moved
    further than the step limit times the number of steps."""
    m, x, _ = _model(PVDF, [T, GP, T, GM], [5.0, 9.6, 90.0, 30.0, 200.0, 1.1, 1.0])
    rng = np.random.default_rng(3)
    x = x.copy()
    x[m.index("tors")] += rng.uniform(-8, 8, m.nc * m.P)
    x = m.project(x, m.index("internal"))
    free = np.setdiff1d(np.arange(m.n_var), np.concatenate([m.index("xy", 1), m.index("z", 1)]))
    e0 = m.physical_energy(x)
    y, gmax, n = m.descend(x, free, [0, 1], steps=25, max_step=0.05, gtol=1e-3)
    assert n == 25 and m.physical_energy(y) < e0
    assert max(m.closure_residual(y)) < 1e-7
    assert np.abs((y - x) / m.scale(x)).max() <= 25 * 0.05 + 1e-6
    np.testing.assert_array_equal(y[m.index("xy", 1)], x[m.index("xy", 1)])  # what is not free does not move


def test_the_climbing_image_waits_for_a_relaxed_band():
    """With ``climb_force`` set the climbing image is not switched on while the band's largest
    force is above it, whatever ``climb_after`` says; once the band is under it, it is."""
    a, b = np.array([-0.5582236, 1.4417258]), np.array([-0.0500108, 0.4666941])
    images = np.linspace(a, b, 11)
    held = MC.neb(_muller_brown, images, k_spring=200.0, climb_after=0, climb_force=1e-12, fmax=0.05, maxiter=30, max_step=0.02)
    assert held.climbing == -1 and not held.converged
    res = MC.neb(_muller_brown, images, k_spring=200.0, climb_after=0, climb_force=20.0, fmax=0.05, maxiter=3000, max_step=0.02)
    assert res.converged and res.climbing > 0
    assert res.energies.max() == pytest.approx(-40.6648, abs=0.01)


def test_a_trust_box_bounds_one_relaxation_and_its_rounds_find_the_basins_minimum():
    """``relax(trust=t)`` moves no free variable further than ``t`` scaled units from its start and
    says when it stopped on the box.  ``relax_in_basin`` repeats such rounds from each result: the
    energy never rises between rounds, the last round is not on its box, and from a start inside
    a basin the rounds end at that basin's minimum (here the one the unbounded relaxation of a
    nearby start reaches).  The quasi-Newton iteration without a box can step over the rim of a
    shallow basin (examples/dilute_kink.py met it); with one, it cannot move past the box."""
    m, x, _ = _model(PVDF, [T, GP, T, GM], [5.0, 9.6, 90.0, 30.0, 200.0, 1.1, 1.0])
    free = np.setdiff1d(np.arange(m.n_var), np.concatenate([m.index("xy", 1), m.index("z", 1)]))
    sc = m.scale(x)
    x1, info = m.relax(x, free=free, trust=0.2, presteps=0, first_step=None, maxiter=80)
    assert info["at_trust_bound"] and m.physical_energy(x1) < m.physical_energy(x)
    assert np.abs((x1 - x) / sc)[free].max() <= 0.2 + 0.02  # the box, and the final projection onto the closure
    np.testing.assert_array_equal(x1[m.index("xy", 1)], x[m.index("xy", 1)])
    # a minimum, a start inside its basin, and the rounds back to it
    xm, im = m.relax(x, free=free, maxiter=600)
    xm, im = m.relax(xm, free=free, maxiter=600)
    assert "success" in im["message"].lower()
    start = xm.copy()
    start[free] += np.random.default_rng(5).uniform(-0.25, 0.25, len(free)) * sc[free]
    start = m.project(start, m.index("internal"))
    energies = [m.physical_energy(start)]
    xb, ib = m.relax_in_basin(start, free=free, trust=0.1, rounds=60, maxiter=200, callback=lambda k, y, i: energies.append(m.physical_energy(y)))
    assert ib["converged"] and not ib["at_trust_bound"] and ib["rounds"] >= 3
    assert all(b <= a + 1e-7 for a, b in zip(energies, energies[1:]))
    assert abs(m.physical_energy(xb) - m.physical_energy(xm)) < 2e-3
    assert max(m.closure_residual(xb)) < 1e-6
