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


def test_relaxation_lowers_the_energy_and_closes_the_chains():
    m, x, _ = _model(PVDF, [T, GP, T, GM], [5.0, 9.6, 90.0, 30.0, 200.0, 1.1, 1.0], pen_rot=5.0, k_c=2000.0)
    e0 = m.terms(x)
    x1, info = m.relax(x, maxiter=150)
    e1 = m.terms(x1)
    assert e1["total"] < e0["total"] - 1.0
    assert np.abs(e1["c_mismatch"]).max() < 0.02 and e1["rotation_error"].max() < 0.5
    held, _ = m.relax(x, free=m.index("place"), maxiter=10)  # only the placements move
    np.testing.assert_array_equal(held[m.index("internal")], x[m.index("internal")])
    np.testing.assert_array_equal(held[:4], x[:4])


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
