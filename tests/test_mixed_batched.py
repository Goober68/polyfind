"""The batched evaluation of a mixed cell, held to the loop form it replaces on a device.

``polyfind.mixedcell`` evaluates a cell either chain by chain and chain pair by chain pair -- the
form every stored number was computed with -- or as a fixed number of array operations whatever
the number of chains, with the pair sums on an array module.  Everything here is the second
against the first: the placed cell's energy and both gradients, the model's value and gradient
with every term switched on, and each new piece against the routine it stands in for.  The array
module is NumPy or ``device_double`` (NumPy underneath, as strict as CuPy on top); what only a
real device can show is ``tests/test_gpu.py``.
"""
import device_double as dd
import numpy as np
import pytest
from test_mixedcell import CASES, _model, _run

from polyfind import chainjac as CJ
from polyfind import ewald as EW
from polyfind import mixedcell as MC
from polyfind import pack as pack_mod
from polyfind.ewald import EwaldSpec
from polyfind.fitting import PRESET_PARAMETERS
from polyfind.forcefield import COULOMB, SimpleFF
from polyfind.linegroup import _windows, repeat_chains, repeat_chains_and_jacobian
from polyfind.polarizability import Polarizable
from polyfind.polymers import PVDF, THREE_STATE, pvdf_belt

T, GP, GM = 0, 1, 2
BELT4 = "pvdf-dft-valence-belt4"
TILINGS = [(1, 1), (2, 1), (2, 2)]


def _moved(m, x, seed=0, size=0.02):
    return x + size * np.random.default_rng(seed).standard_normal(x.size) * m.scale(x)


def _close(a, b, rel=1e-11):
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    np.testing.assert_allclose(a, b, rtol=0, atol=rel * max(1.0, float(np.abs(b).max())))


# ----------------------------------------------------------------- the placed cell
@pytest.mark.parametrize("case", range(len(CASES)))
@pytest.mark.parametrize("tiling", TILINGS)
@pytest.mark.parametrize("module", ["numpy", "double"])
def test_batched_placed_energy_is_the_loop_form(case, tiling, module):
    """Energy, Cartesian gradient and all nine lattice entries: truncated Coulomb and Ewald with
    induced dipoles, flipped chains, repeats short enough that several images along the chain are
    inside the cutoff."""
    xp = np if module == "numpy" else dd

    def body(polymer, seq, params, kw):
        m, x, _ = _model(polymer, seq, params, nx=tiling[0], ny=tiling[1], **kw)
        latn, frames, _, Pn, _ = m.build(_moved(m, x, case))
        E0, g0, l0 = MC.placed_energy_and_grad(m.pk, frames, Pn, latn)
        timings: dict = {}
        E1, g1, l1 = MC.placed_energy_and_grad(m.pk, frames, Pn, latn, xp=xp, timings=timings)
        assert E1 == pytest.approx(E0, rel=1e-11, abs=1e-9)
        _close(g1, g0)
        _close(l1, l0)
        assert isinstance(g1, np.ndarray) and isinstance(l1, np.ndarray) and "pair sums" in timings

    _run(CASES[case], body)


def test_pair_blocks_in_small_stacks_give_the_same_sums():
    """The stacks of blocks are a matter of memory, not of the result."""
    def body(polymer, seq, params, kw):
        m, x, _ = _model(polymer, seq, params, nx=2, ny=2, **kw)
        latn, _, _, Pn, _ = m.build(_moved(m, x, 4))
        P3 = Pn.reshape(m.nc, m.n, 3)
        K = int(np.ceil(m.pk.rc / latn[2, 2])) + 1
        whole = MC._pair_sums(m.pk, P3, latn, m.flips, K, np, block_elems=10 ** 9)
        for elems in (1, 3 * m.n * m.n):
            part = MC._pair_sums(m.pk, P3, latn, m.flips, K, np, block_elems=elems)
            assert part[0] == pytest.approx(whole[0], rel=1e-12)
            _close(part[1], whole[1], rel=1e-12)
            _close(part[2], whole[2], rel=1e-12)

    _run(CASES[3], body)


def test_chain_pair_blocks_hold_every_pair_inside_the_cutoff():
    """No pair of sites inside the cutoff lies outside the blocks, and each block is listed once."""
    def body(polymer, seq, params, kw):
        m, x, _ = _model(polymer, seq, params, nx=2, ny=1, **kw)
        latn, _, _, Pn, _ = m.build(_moved(m, x, 2))
        nc, n, rc = m.nc, m.n, float(m.pk.rc)
        P3 = Pn.reshape(nc, n, 3)
        bs, bt, bi, bj, sh = MC._chain_pair_blocks(P3, latn, rc)
        listed = set(zip(bs.tolist(), bt.tolist(), bi.tolist(), bj.tolist()))
        assert len(listed) == bs.size
        np.testing.assert_allclose(sh, bi[:, None] * latn[0, :2] + bj[:, None] * latn[1, :2], atol=1e-12)
        for s in range(nc):
            for t in range(s, nc):
                for i in range(-6, 7):
                    for j in range(-6, 7):
                        if s == t and (i < 0 or (i == 0 and j <= 0)):
                            assert (s, t, i, j) not in listed
                            continue
                        d = P3[s][:, None, :2] - P3[t][None, :, :2] - (i * latn[0, :2] + j * latn[1, :2])
                        if float((d * d).sum(axis=-1).min()) < rc * rc:
                            assert (s, t, i, j) in listed

    _run(CASES[2], body)


# ------------------------------------------------------------------------ the model
@pytest.mark.parametrize("case", range(len(CASES)))
@pytest.mark.parametrize("module", ["numpy", "double"])
def test_batched_model_is_the_loop_model(case, module):
    """Value, every component of the gradient and every reported part, with both closure
    penalties and their multipliers on."""
    xp = np if module == "numpy" else dd

    def body(polymer, seq, params, kw):
        m, x, _ = _model(polymer, seq, params, nx=2, ny=1, pen_rot=0.7, k_c=3.0, **kw)
        rng = np.random.default_rng(case)
        m.set_multipliers((0.3 * rng.standard_normal((m.nc, 3)), 0.2 * rng.standard_normal(m.nc)))
        x = _moved(m, x, case)
        assert m.batched_module() is None
        p0: dict = {}
        e0, g0 = m.value_and_grad(x, p0)
        m.xp = xp
        m.timings = {}
        p1: dict = {}
        e1, g1 = m.value_and_grad(x, p1)
        assert e1 == pytest.approx(e0, rel=1e-11, abs=1e-9)
        _close(g1, g0)
        assert set(p1) == set(p0)
        for key, value in p0.items():
            _close(p1[key], value, rel=1e-10)
        assert {"build", "pair sums", "assembly"} <= set(m.timings)

    _run(CASES[case], body)


def test_batched_gradient_is_the_derivative_of_the_batched_energy():
    def body(polymer, seq, params, kw):
        m, x, _ = _model(polymer, seq, params, nx=2, ny=1, pen_rot=0.4, k_c=2.0, **kw)
        m.xp = np
        x = _moved(m, x, 7)
        _, g = m.value_and_grad(x)
        d = np.random.default_rng(1).standard_normal(x.size) * m.scale(x)
        h = 1e-4
        e = [m.value_and_grad(x + k * h * d)[0] for k in (-2, -1, 1, 2)]
        slope = (e[0] - 8.0 * e[1] + 8.0 * e[2] - e[3]) / (12.0 * h)
        assert slope == pytest.approx(float(g @ d), rel=2e-6)

    _run(CASES[3], body)


def test_model_follows_the_module_of_its_ewald_sums(monkeypatch):
    """The loops by default; the batched form where the Ewald sums are on a device (so the device
    is not left waiting for the host), or on the host when asked for by the switch."""
    def body(polymer, seq, params, kw):
        m, _, _ = _model(polymer, seq, params, **kw)
        assert m.batched_module() is None
        ew = m.pk._ewald
        try:
            ew.xp = dd
            assert m.batched_module() is None  # too few sites for the device to be asked
            ew.device_min_sites = 0
            assert m.batched_module() is dd
        finally:
            ew.xp, ew.device_min_sites = np, 256
        monkeypatch.setenv("POLYFIND_MIXED_BATCHED", "1")
        assert m.batched_module() is np
        m.xp = "loops"  # the loops, whatever the switch says
        assert m.batched_module() is None
        monkeypatch.delenv("POLYFIND_MIXED_BATCHED")
        m.xp = dd
        assert m.batched_module() is dd
        m.xp = "cuda"
        with pytest.raises(ValueError, match="array module"):
            m.batched_module()

    _run(CASES[2], body)


def test_an_evaluation_on_a_device_carries_per_site_arrays_only():
    """After the first evaluation has sent the cell's tables, an evaluation sends and brings back
    nothing larger than a per-site array, in a number of transfers that does not grow with the
    number of chains."""
    counts = {}
    with PRESET_PARAMETERS[BELT4].applied():
        polymer = pvdf_belt(0.7, "pbe")
        for nx in (1, 2):
            m, x, _ = _model(polymer, [T, GP, T, GM], np.array([5.0, 9.5, 91.0, 30.0, 215.0, 2.2, 1.0]), nx=nx, ny=1,
                             lj_cutoff="force", coulomb="ewald", ewald=EwaldSpec(), polarizable=Polarizable())
            x = _moved(m, x, 3)
            ew = m.pk._ewald
            ew.xp, ew.device_min_sites = dd, 0
            try:
                m.value_and_grad(x)
                dd.reset_traffic()
                e, g = m.value_and_grad(x + 1e-3)
                seen = dict(dd.traffic)
            finally:
                ew.xp, ew.device_min_sites = np, 256
                ew.release()
            N = m.pk.N
            assert isinstance(g, np.ndarray) and np.isfinite(e)
            assert seen["largest_to_host"] <= 8 * 3 * N
            counts[nx] = (seen["to_device"], seen["to_host"])
    assert counts[1] == counts[2]


# ------------------------------------------------------------------- the valence terms
def _valence_packer(n_chains=4):
    import polyfind.fitting  # noqa: F401 - registers the fitted presets

    ff = SimpleFF.from_preset("pvdf-dft-valence")
    seq = [T, GP, T, GM]
    tors = np.array([THREE_STATE.angles[s] for s in seq], dtype=float)
    ang0 = np.array([PVDF.backbone[k].backbone_angle for k in range(PVDF.bonds_per_repeat)], dtype=float)
    ch = repeat_chains(PVDF, "TGTG'", tors[None], ang0)[0]
    return ch, tors, ang0, pack_mod.CrystalPacker(ch, n_chains=n_chains, valence=ff, lj_cutoff="force")


def test_valence_terms_of_many_geometries_at_once():
    """One call for all chains is the sum of one call per chain, in each chain's own frame and --
    with the image sign -- in the placed frame of a flipped chain."""
    ch, _, _, pk = _valence_packer()
    val = pk._valence
    assert val is not None and val.n_bonds and val.n_angles
    rng = np.random.default_rng(0)
    X = ch.coords[None] + 0.05 * rng.standard_normal((5, ch.n_atoms, 3))
    c = ch.c * 1.01
    single = [val.energy_and_grad(X[k], c) for k in range(5)]
    E, g, gc = val.energy_and_grad_batch(X, c)
    assert E == pytest.approx(sum(s[0] for s in single), rel=1e-13)
    _close(g, np.stack([s[1] for s in single]), rel=1e-13)
    assert gc == pytest.approx(sum(s[2] for s in single), rel=1e-12)
    # the same chains placed: turned about z, and every second one flipped
    flips = np.array([False, True, False, True, True])
    frames = np.stack([MC.placing_frame(40.0 * k + 10.0, bool(flips[k])) for k in range(5)])
    placed = np.einsum("sni,sji->snj", X, frames) + rng.standard_normal((5, 1, 3))
    Ep, gp, gcp = val.energy_and_grad_batch(placed, c, image_sign=np.where(flips, -1.0, 1.0))
    assert Ep == pytest.approx(E, rel=1e-12) and gcp == pytest.approx(gc, rel=1e-10)
    _close(np.einsum("sni,sij->snj", gp, frames), g, rel=1e-11)


@pytest.mark.parametrize("per_chain", [False, True])
def test_batched_cell_with_valence_terms_is_the_loop_form(per_chain):
    ch, tors, ang0, pk = _valence_packer()
    st, flips = MC.supercell_state([5.0, 9.6, 90.0, 30.0, 200.0, 1.1, 1.0], ch.c, tors, ang0, 2, 1)
    valences = [pk._valence, None, pk._valence, pk._valence] if per_chain else None
    m = MC.MixedModel(PVDF, "TGTG'", pk, flips, [ch] * 4, angle_stiffness=105.0, valences=valences)
    assert flips.any()
    latn, frames, _, Pn, _ = m.build(_moved(m, st.vector(), 5))
    E0, g0, l0 = MC.placed_energy_and_grad(pk, frames, Pn, latn, valences=valences)
    for xp in (np, dd):
        timings: dict = {}
        E1, g1, l1 = MC.placed_energy_and_grad(pk, frames, Pn, latn, valences=valences, xp=xp, timings=timings)
        assert E1 == pytest.approx(E0, rel=1e-11)
        _close(g1, g0)
        _close(l1, l0)
        assert timings["valence"] > 0.0


# ------------------------------------------------------------------ the exclusion lists
@pytest.mark.parametrize("damping", [None, 3.5])
def test_exclusion_lists_are_the_dense_corrections(damping):
    """Both bonded-exclusion corrections of all chains from the flat lists, against the dense
    routines chain by chain, with flipped chains and a charge penetration term."""
    with PRESET_PARAMETERS[BELT4].applied():
        polymer = pvdf_belt(0.7, "pbe")
        m, x, _ = _model(polymer, [T, GP, T, GM], np.array([5.0, 9.5, 91.0, 30.0, 215.0, 2.2, 1.0]), nx=2, ny=1,
                         lj_cutoff="force", coulomb="ewald", ewald=EwaldSpec(charge_damping=damping), polarizable=Polarizable())
        pk = m.pk
        latn, frames, _, Pn, _ = m.build(_moved(m, x, 9))
        n, N, nc = pk.n, pk.N, pk.n_chains
        cz = float(latn[2, 2])
        K = int(np.ceil(pk.rc / cz)) + 1
        flips = frames[:, 2, 2] < 0.0
        assert flips.any() and not flips.all() and pk._es_idx is not None
        q_cell = np.asarray(pk._q_cell, dtype=float)
        pref = COULOMB / pk.eps_r
        Sk = pk._scale_column_np(K)
        active = np.zeros(N, dtype=bool)
        active[pk._es_idx] = True
        pairs = pk._exclusion_pairs(K, flips)
        assert pk._exclusion_pairs(K, flips) is pairs and len(pairs) > 0
        assert active[pairs.i].all() and active[pairs.j].all()
        p = np.zeros((N, 3))
        p[active] = 0.05 * np.random.default_rng(5).standard_normal((int(active.sum()), 3))
        e1 = c1 = e2 = c2 = 0.0
        g1, g2, F2 = np.zeros((N, 3)), np.zeros((N, 3)), np.zeros((N, 3))
        for s in range(nc):
            sl = slice(s * n, (s + 1) * n)
            e, gX, gc, _ = pk._exclusion(Pn[sl] @ frames[s], q_cell[:n], cz, K, grad=True)
            e1, c1 = e1 + e, c1 + gc
            g1[sl] = gX @ frames[s].T
            e, gX, gc, _, F = EW.charge_dipole_exclusion(Pn[sl], q_cell[sl], cz, Sk[::-1] if flips[s] else Sk,
                                                         dipoles=p[sl], pref=pref, grad=True)
            e2, c2 = e2 + e, c2 + gc
            g2[sl], F2[sl] = gX, F
        E1, G1, C1 = EW.exclusion_correction_pairs(Pn, q_cell, cz, pairs, pref=pref, charge_damping=pk._ewald.charge_damping)
        E2, G2, C2, FF = EW.charge_dipole_exclusion_pairs(Pn, q_cell, cz, pairs, dipoles=p, pref=pref, grad=True)
    assert E1 == pytest.approx(e1, rel=1e-12) and C1 == pytest.approx(c1, rel=1e-11)
    _close(G1, g1, rel=1e-12)
    assert E2 == pytest.approx(e2, rel=1e-12) and C2 == pytest.approx(c2, rel=1e-11)
    _close(G2, g2, rel=1e-12)
    _close(FF[active], F2[active], rel=1e-12)
    assert not FF[~active].any()  # no dipole there for a field to act on
    only_field = EW.charge_dipole_exclusion_pairs(Pn, q_cell, cz, pairs, pref=pref)
    assert only_field[0] == 0.0 and only_field[1] is None and np.array_equal(only_field[3], FF)


def test_induced_dipoles_with_the_exclusion_lists_are_the_dense_ones():
    """The packer's own dipole routine: the lists against the per-chain stacks."""
    with PRESET_PARAMETERS[BELT4].applied():
        polymer = pvdf_belt(0.7, "pbe")
        m, x, _ = _model(polymer, [T, T], np.array([4.75, 8.48, 90.0, 180.0, 180.0, 1.3, 1.0]), nx=2, ny=1,
                         lj_cutoff="force", coulomb="ewald", ewald=EwaldSpec(), polarizable=Polarizable())
        pk = m.pk
        latn, frames, _, Pn, _ = m.build(_moved(m, x, 2))
        flips = frames[:, 2, 2] < 0.0
        K = int(np.ceil(pk.rc / latn[2, 2])) + 1
        dense = pk._polarize(Pn, pk._q_cell, latn, flips.astype(float), grad=True)
        lists = pk._polarize(Pn, pk._q_cell, latn, flips.astype(float), grad=True, pairs=pk._exclusion_pairs(K, flips))
        with pytest.raises(ValueError, match="charge derivative"):
            pk._polarize(Pn, pk._q_cell, latn, flips.astype(float), grad=True, charge_grad=True, pairs=pk._exclusion_pairs(K, flips))
    assert lists[0] == pytest.approx(dense[0], rel=1e-12)
    _close(lists[1], dense[1], rel=1e-11)
    _close(lists[2][pk._es_idx], dense[2][pk._es_idx], rel=1e-12)
    _close(lists[3][0], dense[3][0], rel=1e-11)
    _close(lists[3][1], dense[3][1], rel=1e-11)
    assert lists[3][3] == pytest.approx(dense[3][3], rel=1e-10)


# ------------------------------------------------------------ the chain build, read backwards
def _chains_and_refs(nc=5, seed=3):
    seq = [T, GP, T, GM, T, T]
    polymer = pvdf_belt(0.7, "pbe")
    tors0 = np.array([THREE_STATE.angles[s] for s in seq], dtype=float)
    ang0 = np.array([polymer.backbone[k].backbone_angle for k in range(polymer.bonds_per_repeat)], dtype=float)
    ref = repeat_chains(polymer, "ref", tors0[None], ang0)[0]
    rng = np.random.default_rng(seed)
    tors = tors0[None] + 3.0 * rng.standard_normal((nc, tors0.size))
    angs = ang0[None] + 1.0 * rng.standard_normal((nc, ang0.size))
    refs = []
    for _ in range(nc):
        a = rng.uniform(0.0, 2.0 * np.pi)
        Rz = np.array([[np.cos(a), -np.sin(a), 0.0], [np.sin(a), np.cos(a), 0.0], [0.0, 0.0, 1.0]])
        refs.append(ref.coords @ Rz.T + np.array([0.0, 0.0, rng.uniform(-1.0, 1.0)]))
    return polymer, tors, angs, refs


def test_torque_sums_are_the_contraction_with_the_jacobian():
    """A gradient carried back to the torsions and angles through the chain build, against the
    product with the full Jacobian of the coordinates -- each chain aligned to its own reference."""
    polymer, tors, angs, refs = _chains_and_refs()
    nc, P = tors.shape
    R = angs.shape[1]
    forward = [repeat_chains_and_jacobian(polymer, "x", tors[s][None], angs[s], refs[s]) for s in range(nc)]
    coords, bb, owner, is_pendant = CJ.oligomer(polymer, tors, angs)
    built = pack_mod._finish_repeat_chains(polymer, "x", tors, coords, 0.5, np.stack(refs))
    tpl = pack_mod._template(polymer, P)
    atoms0 = np.asarray(pack_mod._block_atoms(tpl, P, 0)[1])
    sx, sy = _windows(P)
    back = np.asarray(tpl.backbone)
    window = np.concatenate([back[sx], back[sy]])
    nw = back[sx].size
    jw = CJ.oligomer_jacobian(coords, bb, owner, is_pendant, P, R, atoms=window)
    tr = CJ.window_transform_jacobian(bb[:, sx], bb[:, sy], jw[:, :, :nw], jw[:, :, nw:])
    G = np.random.default_rng(1).standard_normal((nc, atoms0.size, 3))
    gv, Gp = CJ.aligned_block_vjp(coords[:, atoms0], tr["axis"], tr["daxis"], np.stack(refs), G)
    gv = gv + CJ.oligomer_vjp(coords, bb, owner, is_pendant, P, R, atoms0, Gp)
    for s in range(nc):
        chains, dX, dc, w, dw = forward[s]
        # the chains and the repeat transform are the per-chain build's own, bit for bit
        assert np.array_equal(built[s].coords, chains[0].coords) and built[s].c == chains[0].c
        assert np.array_equal(tr["dc"][s], dc[0])
        assert np.array_equal(np.degrees(tr["w"][s]), w[0]) and np.array_equal(np.degrees(tr["dw"][s]), dw[0])
        want = np.einsum("ni,vni->v", G[s], dX[0])
        np.testing.assert_allclose(gv[s], want, rtol=0, atol=1e-12 * float(np.abs(want).max()))


def test_jacobian_of_some_atoms_is_those_columns_of_the_whole():
    polymer, tors, angs, _ = _chains_and_refs(nc=3)
    P, R = tors.shape[1], angs.shape[1]
    coords, bb, owner, is_pendant = CJ.oligomer(polymer, tors, angs)
    full = CJ.oligomer_jacobian(coords, bb, owner, is_pendant, P, R)
    whole = CJ.oligomer_and_jacobian(polymer, tors, angs)
    assert np.array_equal(whole[2], full) and np.array_equal(whole[0], coords) and np.array_equal(whole[3], owner)
    atoms = np.array([0, 7, 33, coords.shape[1] // 2, coords.shape[1] - 1])
    assert np.array_equal(CJ.oligomer_jacobian(coords, bb, owner, is_pendant, P, R, atoms=atoms), full[:, :, atoms])


def test_closure_is_the_per_chain_closure_bit_for_bit():
    """The closure residuals and their Jacobian come from one light build of all chains; the
    numbers are those of the chain-by-chain build with the whole Jacobian."""
    def body(polymer, seq, params, kw):
        m, x, _ = _model(polymer, seq, params, nx=2, ny=1, **kw)
        x = _moved(m, x, 6)
        st = m.state(x)
        P, R, base = m.P, m.R, 4 + 4 * m.nc
        for chains in (None, [3, 0], [2]):
            ch = range(m.nc) if chains is None else chains
            h0, J0 = np.zeros(4 * len(ch)), np.zeros((4 * len(ch), m.n_var))
            for k, s in enumerate(ch):
                built, _, dc, w, dw = repeat_chains_and_jacobian(m.polymer, m.name, st.tors[s][None], st.angs[s], m.refs[s])
                cols = slice(base + s * (P + R), base + (s + 1) * (P + R))
                h0[4 * k:4 * k + 3], J0[4 * k:4 * k + 3, cols] = w[0], dw[0].T
                h0[4 * k + 3] = MC.CLOSURE_C * (float(built[0].c) - st.lattice[3])
                J0[4 * k + 3, cols], J0[4 * k + 3, 3] = MC.CLOSURE_C * dc[0], -MC.CLOSURE_C
            h1, J1 = m.closure(x, chains)
            assert np.array_equal(h1, h0) and np.array_equal(J1, J0)
            # beside the batched evaluation the derivative is taken from running sums along the
            # backbone: the same residuals, the same Jacobian to rounding
            m.xp = np
            h2, J2 = m.closure(x, chains)
            m.xp = None
            assert np.array_equal(h2, h0)
            np.testing.assert_allclose(J2, J0, rtol=0, atol=1e-9 * float(np.abs(J0).max()))
        assert m.closure(x, [])[0].size == 0

    _run(CASES[3], body)


@pytest.mark.parametrize("n_bonds", [1, 2, 4, 6])
def test_transform_derivative_from_running_sums_is_the_one_from_the_window_jacobian(n_bonds):
    """Every size of repeat, the two-block windows of the short ones included."""
    polymer = PVDF
    rng = np.random.default_rng(n_bonds)
    R = polymer.bonds_per_repeat
    tors = rng.uniform(-180.0, 180.0, (4, n_bonds))
    angs = np.array([polymer.backbone[k].backbone_angle for k in range(R)])[None] + rng.standard_normal((4, R))
    coords, bb, owner, is_pendant = CJ.oligomer(polymer, tors, angs)
    sx, sy = _windows(n_bonds)
    back = np.asarray(pack_mod._template(polymer, n_bonds).backbone)
    nw = back[sx].size
    jw = CJ.oligomer_jacobian(coords, bb, owner, is_pendant, n_bonds, R, atoms=np.concatenate([back[sx], back[sy]]))
    a = CJ.window_transform_jacobian(bb[:, sx], bb[:, sy], jw[:, :, :nw], jw[:, :, nw:])
    b = CJ.backbone_transform_jacobian(bb, sx, sy, n_bonds, R)
    for key in ("R", "t", "c", "axis", "w"):
        assert np.array_equal(a[key], b[key])
    for key in ("delta", "dt", "dc", "daxis", "dw"):
        np.testing.assert_allclose(b[key], a[key], rtol=0, atol=1e-10 * max(1e-3, float(np.abs(a[key]).max())))


# ------------------------------------------------------------------ a copolymer cell
@pytest.mark.parametrize("row", ["vdf9-trfe4-cfe1-14", "vdf10-trfe4-14"])
@pytest.mark.parametrize("module", ["numpy", "double"])
def test_batched_form_on_a_copolymer_cell_is_the_loop_form(row, module):
    """Everything above is PVDF.  The stored kinked cells of the terpolymer and of its matrix,
    tiled to four chains and displaced in every variable (examples/mixed_forms_check.py): sites
    that differ along the chain, a 28-bond repeat, a chiral chain, belt sites on some atoms only.
    Energy and gradient to rounding, the closure bit for bit."""
    import importlib.util
    import json
    from pathlib import Path

    root = Path(__file__).resolve().parents[1]
    stored = root / "docs" / "predictions" / "kink_migration_belt4.json"
    if not stored.exists() or f"{row}@0" not in json.loads(stored.read_text()).get("kinked_references", {}):
        pytest.skip("no stored kinked cell of this row")
    spec = importlib.util.spec_from_file_location("mixed_forms_check", root / "examples" / "mixed_forms_check.py")
    FC = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(FC)
    table = json.loads(stored.read_text())
    with PRESET_PARAMETERS[table["preset"]].applied():
        r = FC.forms(row, 2, 1, table=table, module=np if module == "numpy" else dd)
    assert r["chains"] == 4 and r["variables"] == 4 + 4 * (4 + 56) and r["module"] == ("numpy" if module == "numpy" else dd.__name__)
    assert r["above_the_sheet"] > 100.0  # a displaced geometry, far from the minimum the cell was stored at
    assert FC.agree(r), r
    assert r["closure_difference"] == 0.0
