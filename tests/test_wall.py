"""The buffered 14-7 pair term (polyfind.pack.wall_b147) as a second form beside 12-6.

What is held here: the term itself, that a packer with it is a consistent energy with exact
gradients in every form the sums are written in (the packer's kernels, the mixed cell's loops and
its batched kernel, the supercell sum), that nothing changes when it is not asked for, and that a
parameter set carries it to every packer built under it.  That it is the better wall is a matter
of DFT, in tests/test_interchain.py.
"""
import numpy as np
import pytest

from polyfind import crystalfit as C
from polyfind import fitting as F
from polyfind import mixedcell as MC
from polyfind import pack as pack_mod
from polyfind.ewald import EwaldSpec
from polyfind.linegroup import repeat_chains
from polyfind.pack import WALL_GAMMA, wall_b147
from polyfind.polarizability import Polarizable
from polyfind.polymers import PVDF, THREE_STATE, pvdf_belt

T, GP, GM = 0, 1, 2
BELT4 = "pvdf-dft-valence-belt4"
EWALD = dict(lj_cutoff="force", coulomb="ewald", ewald=EwaldSpec(), polarizable=Polarizable())
WALL = ("b14-7", 0.2)


def _chain(polymer, seq):
    tors = np.array([THREE_STATE.angles[s] for s in seq], dtype=float)
    ang0 = np.array([polymer.backbone[k].backbone_angle for k in range(polymer.bonds_per_repeat)], dtype=float)
    name = "".join(THREE_STATE.names[s] for s in seq)
    return repeat_chains(polymer, name, tors[None], ang0)[0], name, tors, ang0


def _model(polymer, seq, params, nx=1, ny=1, xp=None, shear=False, **kw):
    ch, name, tors, ang0 = _chain(polymer, seq)
    nc = 2 * nx * ny
    pk = pack_mod.CrystalPacker(ch, n_chains=nc, **kw)
    st, flips = MC.supercell_state(params, ch.c, tors, ang0, nx, ny)
    if shear:
        st.shear = np.zeros(2)
    return MC.MixedModel(polymer, name, pk, flips, [ch] * nc, angle_stiffness=105.0, xp=xp, shear=shear), st.vector()


def test_the_term_has_the_wells_depth_at_the_diameter_and_a_softer_wall_inside_it():
    x, D = 3.0, 0.05
    r = np.linspace(1.2, 7.5, 600)
    for delta in (0.07, 0.2, 0.35):
        v, dv = wall_b147(np, r, D, 1.0 / x, delta, grad=True)
        assert wall_b147(np, x, D, 1.0 / x, delta) == pytest.approx(-D, rel=1e-14)
        i = int(np.argmin(v))
        assert abs(r[i] / x - 1.0) < 0.02 and v[i] == pytest.approx(-D, rel=1e-2)  # the minimum is at the diameter to 2%, less for small buffers
        h = 1e-6
        fd = (wall_b147(np, r + h, D, 1.0 / x, delta) - wall_b147(np, r - h, D, 1.0 / x, delta)) / (2 * h)
        np.testing.assert_allclose(dv, fd, rtol=1e-6, atol=1e-8)
        assert np.isfinite(wall_b147(np, 0.0, D, 1.0 / x, delta)) and wall_b147(np, 0.0, D, 1.0 / x, delta) > 1e3  # finite and repulsive at r = 0
    lj = lambda q: D * ((x / q) ** 12 - 2 * (x / q) ** 6)  # noqa: E731
    inside = np.array([0.6, 0.7, 0.8]) * x  # well inside the minimum; at 0.9 x the wells themselves differ
    soft = [wall_b147(np, inside, D, 1.0 / x, d) for d in (0.07, 0.2, 0.35)]
    assert np.all(soft[0] < lj(inside)) and np.all(soft[1] < soft[0]) and np.all(soft[2] < soft[1])  # softer with the buffer
    # the attraction is the form's own: weaker than 12-6's outside the well, and falling as r^-7, not r^-6
    far = wall_b147(np, np.array([3.0, 6.0]) * x, D, 1.0 / x, 0.07)
    assert far[0] / far[1] == pytest.approx(2.0 ** 7, rel=0.1) and 0.7 < wall_b147(np, 1.4 * x, D, 1.0 / x, 0.07) / lj(1.4 * x) < 0.9
    assert wall_b147(np, r, 0.0, 1.0 / x, 0.2).tolist() == [0.0] * r.size


def test_a_negative_depth_is_a_wall_alone():
    """How a repulsion-only site reaches a kernel with two tables: no ``- 2``, positive everywhere."""
    x, D, delta = 2.6, 0.02, 0.2
    r = np.linspace(1.0, 7.0, 300)
    w, dw = wall_b147(np, r, -D, 1.0 / x, delta, grad=True)
    v = wall_b147(np, r, D, 1.0 / x, delta)
    a7 = ((1.0 + delta) / (r / x + delta)) ** 7
    np.testing.assert_allclose(w - v, 2.0 * D * a7, rtol=1e-12)  # the two differ by the dispersion part alone
    assert np.all(w > 0.0) and np.all(dw < 0.0) and wall_b147(np, x, -D, 1.0 / x, delta) == pytest.approx(D, rel=1e-14)
    h = 1e-6
    np.testing.assert_allclose(dw, (wall_b147(np, r + h, -D, 1.0 / x, delta) - wall_b147(np, r - h, -D, 1.0 / x, delta)) / (2 * h), rtol=1e-6, atol=1e-9)
    mixed = wall_b147(np, np.array([2.0, 2.0]), np.array([D, -D]), 1.0 / x, delta)  # per pair, in one array
    assert mixed[0] == wall_b147(np, 2.0, D, 1.0 / x, delta) and mixed[1] == wall_b147(np, 2.0, -D, 1.0 / x, delta)


def test_a_packer_without_a_wall_is_the_packer_it_always_was():
    ch = _chain(PVDF, [T, GP, T, GM])[0]
    params = np.array([[5.0, 9.6, 90.0, 30.0, 200.0, 1.1, 1.0]])
    a = pack_mod.CrystalPacker(ch, n_chains=2)
    b = pack_mod.CrystalPacker(ch, n_chains=2, wall=None)
    c = pack_mod.CrystalPacker(ch, n_chains=2, wall=())
    assert a.wall is None and b.wall is None and c.wall is None
    assert a.energy(params)[0] == b.energy(params)[0] == c.energy(params)[0]
    w = pack_mod.CrystalPacker(ch, n_chains=2, wall=WALL)
    assert w.wall == ("b14-7", 0.2, WALL_GAMMA) and abs(w.energy(params)[0] - a.energy(params)[0]) > 0.1
    with pytest.raises(ValueError, match="unknown wall"):
        pack_mod.CrystalPacker(ch, n_chains=2, wall=("exp-6", 12.0))


@pytest.mark.parametrize("setting", ["dsf", "ewald"])
def test_the_packers_gradient_with_the_wall_matches_central_differences(setting):
    """The packer's own cell: energy_and_grad against the energy, every cell parameter."""
    def body(polymer, kw):
        ch = _chain(polymer, [T, GP, T, GM])[0]
        pk = pack_mod.CrystalPacker(ch, n_chains=2, wall=WALL, **kw)
        p = np.array([5.05, 9.55, 91.0, 30.0, 215.0, 2.2, 1.0])
        e, g, _, _ = pk.energy_and_grad(p)
        assert e == pytest.approx(pk.energy(p[None])[0], rel=1e-10, abs=1e-8)
        for i in range(6):
            h = 1e-5
            pp, pm = p.copy(), p.copy()
            pp[i] += h
            pm[i] -= h
            assert (pk.energy(pp[None])[0] - pk.energy(pm[None])[0]) / (2 * h) == pytest.approx(g[i], rel=2e-5, abs=2e-5)
    if setting == "dsf":
        body(PVDF, dict(lj_cutoff="force"))
    else:
        with F.PRESET_PARAMETERS[BELT4].applied():
            body(pvdf_belt(0.7, "pbe"), EWALD)


@pytest.mark.parametrize("batched", [False, True])
def test_the_mixed_cell_with_the_wall_has_exact_gradients_in_both_forms(batched):
    """Belt sites and all (a wall alone against hydrogen), sheared, every kind of variable."""
    with F.PRESET_PARAMETERS[BELT4].applied():
        m, x = _model(pvdf_belt(0.7, "pbe"), [T, GP, T, GM], [5.0, 9.5, 91.0, 30.0, 215.0, 2.2, 1.0], xp=np if batched else None, shear=True, wall=WALL, **EWALD)
        assert m.pk.wall == ("b14-7", 0.2, WALL_GAMMA) and "X" in m.pk.elements
        rng = np.random.default_rng(3)
        x[m.index("tors")] += rng.uniform(-6, 6, m.nc * m.P)
        x[m.index("angs")] += rng.uniform(-1.5, 1.5, m.nc * m.R)
        x[m.index("phi")] += rng.uniform(-8, 8, m.nc)
        x[m.index("z")] += rng.uniform(-0.3, 0.3, m.nc)
        x[:4] += np.array([0.07, -0.05, 1.3, 0.04])
        x[m.index("shear")] = [0.6, -1.1]
        _, g = m.value_and_grad(x)
        rest = np.setdiff1d(np.arange(4, m.n_var), m.index("shear"))
        for i in np.concatenate([np.arange(4), m.index("shear"), rng.choice(rest, size=10, replace=False)]):
            xp_, xm_ = x.copy(), x.copy()
            xp_[i] += 1e-5
            xm_[i] -= 1e-5
            assert (m.value_and_grad(xp_)[0] - m.value_and_grad(xm_)[0]) / 2e-5 == pytest.approx(g[i], rel=2e-5, abs=2e-5)


def test_every_form_of_the_sums_gives_one_energy_with_the_wall():
    """The packer's cell in the mixed cell (loops and batched), its tiling, and a sheared cell
    against the upright cell twice as wide: the known answers of the 12-6 path, with the wall."""
    with F.PRESET_PARAMETERS[BELT4].applied():
        polymer, seq, params = pvdf_belt(0.7, "pbe"), [T, GP, T, GM], np.array([5.0, 9.5, 91.0, 30.0, 215.0, 2.2, 1.0])
        m, x = _model(polymer, seq, params, wall=WALL, **EWALD)
        e_packer = m.pk.energy_and_grad(params)[0]
        parts: dict = {}
        m.value_and_grad(x, parts)
        assert parts["physical"] == pytest.approx(e_packer, rel=1e-11, abs=1e-9)
        mb, _ = _model(polymer, seq, params, xp=np, wall=WALL, **EWALD)
        assert mb.physical_energy(x) == pytest.approx(parts["physical"], rel=1e-12, abs=1e-9)
        m2, x2 = _model(polymer, seq, params, nx=2, wall=WALL, **EWALD)
        assert m2.physical_energy(x2) / 2 == pytest.approx(parts["physical"], rel=1e-10, abs=1e-8)
        ms, xs = _model(polymer, seq, params, shear=True, wall=WALL, **EWALD)
        st = ms.state(xs)
        st.shear = np.array([0.5, 0.0]) * float(st.lattice[3])
        e_sheared = ms.physical_energy(st.vector())
        big, fl = MC.tile_state(st, ms.flips, 2, 1)
        up = MC.reduce_shear(big)
        assert up.shear is None and abs(e_sheared - parts["physical"]) > 0.05
        assert m2.physical_energy(up.vector()) / 2 == pytest.approx(e_sheared, rel=1e-10, abs=1e-8)


def test_the_supercell_sum_follows_the_packers_wall():
    from polyfind import supercell as SC

    ch = _chain(PVDF, [T, T])[0]
    pk = pack_mod.CrystalPacker(ch, n_chains=2, wall=WALL)
    sums = SC.SupercellEnergy(pk)
    assert sums._wall == pk._wall
    params = np.array([8.6, 4.9, 93.0, 40.0, 40.0, 0.3, 0.0])
    assert abs(sums.check_against_packer(params)) < 1e-8  # the supercell sum of the packer's own cell, per monomer
    assert abs(sums.check_against_packer(params, nx=2, ny=1)) < 1e-8


def test_a_parameter_set_carries_the_wall_to_every_packer_built_under_it():
    x14 = np.concatenate([F.CRYSTAL8_TILT_FITTED_X, [0.2]])
    prm = C.parameters(x14)
    assert prm.wall == ("b14-7", 0.2) and "b14-7" in prm.describe()
    assert C.parameters(np.array(F.CRYSTAL8_TILT_FITTED_X)).wall == () and F.PRESET_PARAMETERS[BELT4].wall == ()
    ch = _chain(pvdf_belt(0.7, "pbe"), [T, T])[0]
    with prm.applied():
        inside = pack_mod.CrystalPacker(ch, n_chains=2, **EWALD)
    with F.PRESET_PARAMETERS[BELT4].applied():
        before = pack_mod.CrystalPacker(ch, n_chains=2, **EWALD)
    outside = pack_mod.CrystalPacker(ch, n_chains=2, **EWALD)
    assert inside.wall == ("b14-7", 0.2, WALL_GAMMA) and before.wall is None and outside.wall is None
    # the same diameters, depths and charges as the 12-6 preset it extends: only the form differs
    np.testing.assert_array_equal(inside._q_cell, before._q_cell)
    np.testing.assert_array_equal(np.asarray(inside.lj_x), np.asarray(before.lj_x))
    # the stage's bookkeeping: fourteen variables, the buffer the last, its ridge from Halgren's 0.07
    names, x0, bounds, sigma = C._config(x14)
    assert names[-1] == "wall_delta" and len(names) == 14 and x0[-1] == 0.07 and bounds[-1].tolist() == [0.0, 0.5]
    assert C.ridge_parts(x14)[-1] == pytest.approx((0.2 - 0.07) / sigma[-1])
    assert C.STAGE12["registry"] and C.STAGE12["angles"]


def test_the_conformer_design_can_see_the_lattice_sides_non_bonded_terms():
    """``ValenceDesign.lattice``: a crystal-stage parameter set in place of the vector's own
    Lennard-Jones terms, increments, pair terms and wall, so that the torsions can be scored and
    fitted again under the potential the lattice has.  ``None`` is the vector as it always was."""
    import importlib.util
    from pathlib import Path

    spec = importlib.util.spec_from_file_location("tf", Path(__file__).with_name("test_fitting.py"))
    TF = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(TF)
    design = F.ValenceDesign(TF._synthetic_frames())
    x = TF._valence_x()
    assert design.lattice is None
    p0 = design.unpack(x)
    assert p0.get("wall") is None and not p0.get("lj_pairs")
    e0, f0, _ = design.evaluate(p0)
    names = sorted(F._UFF_ALL)
    belt4 = F.PRESET_PARAMETERS[BELT4]
    design.lattice = belt4
    p1 = design.unpack(x)
    lj = belt4.lj_dict
    for e in ("C", "H", "F"):
        assert p1["lj_x"][names.index(e)] == lj[e][0] and p1["lj_d"][names.index(e)] == lj[e][1]
    assert p1["wall"] is None and list(p1["lj_pairs"]) == [(names.index("F"), names.index("H"))]  # the belt sites are not atoms
    assert p1["increments"][F.REF_INCREMENTS.index(("C", "F"))] == dict(((a, b), v) for a, b, v in belt4.charge_increments)[("C", "F")]
    np.testing.assert_array_equal(p1["torsion"], p0["torsion"])  # the torsions stay the vector's
    design.lattice = C.parameters(np.concatenate([F.CRYSTAL8_TILT_FITTED_X, [0.2]]))
    p2 = design.unpack(x)
    assert p2["wall"] == (0.2, WALL_GAMMA)
    e2, f2, _ = design.evaluate(p2)
    e1, f1, _ = design.evaluate(p1)
    assert np.abs(e2 - e1).max() > 0.05 and np.abs(e1 - e0).max() > 0.01  # each change is felt
    # the forces are the energy's: one atom moved in one frame, with the wall, against the design's own force
    frames = TF._synthetic_frames()
    from polyfind.forcefield import Frame
    f = frames[1]
    atom, h = 5, 1e-5
    def energy_of(dx):
        c = np.array(f.coords, dtype=float)
        c[atom, 0] += dx
        moved = F.ValenceDesign([frames[0], Frame(**{**f.__dict__, "coords": c})] + frames[2:])
        moved.lattice = design.lattice
        return moved.evaluate(moved.unpack(x))[0][1]
    fd = -(energy_of(h) - energy_of(-h)) / (2 * h)
    first = int(np.flatnonzero(design.atom_frame == 1)[0])
    assert fd == pytest.approx(f2[first + atom, 0], rel=1e-5, abs=1e-5)


def test_the_isolated_chain_calculator_uses_the_parameter_sets_wall():
    """One pair term on both sides: a parameter set with a wall hands it to its SimpleFF, whose
    non-bonded energy is then wall_b147 over the same pairs, diameters and (1-4 scaled) depths."""
    from polyfind.chain import build_chain
    from polyfind.forcefield import SimpleFF

    prm = C.parameters(np.concatenate([F.CRYSTAL8_TILT_FITTED_X, [0.2]]))
    ff, ff12 = prm.simple_ff(), F.PRESET_PARAMETERS[BELT4].simple_ff()
    assert ff.wall == ("b14-7", 0.2) and ff12.wall is None and prm.preset_kwargs()["wall"] == ("b14-7", 0.2)
    chain = build_chain(PVDF, np.radians([172.0, 62.0, 176.0, -58.0, 168.0, 65.0, -175.0]))  # a conformer that does not fold onto itself
    a, b = ff.components(chain), ff12.components(chain)
    top = ff._topology(chain)
    r = np.linalg.norm(chain.coords[top.pairs_i] - chain.coords[top.pairs_j], axis=1)
    assert a["lj"] == pytest.approx(float(wall_b147(np, r, top.lj_d, 1.0 / top.lj_x, 0.2).sum()), rel=1e-12)
    assert abs(a["lj"] - b["lj"]) > 0.05 and a["coulomb"] == pytest.approx(b["coulomb"], rel=1e-12)
    assert ff.energy(chain) - ff12.energy(chain) == pytest.approx(a["lj"] - b["lj"], rel=1e-9, abs=1e-9)
    with pytest.raises(ValueError, match="unknown wall"):
        SimpleFF(wall=("exp-6", 12.0)).components(chain)
