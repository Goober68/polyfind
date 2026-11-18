"""A charge site on a bond (`polyfind.polymers.bond_site_atom`), held to known answers in the kernels.

The site is a point charge on a C-H bond and nothing else: no Lennard-Jones term, no
polarizability, no mass, no bond.  What has to be true of it: a chain that carries uncharged sites
is the chain without them, exactly, in every energy, gradient and placement; a charged site enters
the lattice sums as the point charge it is; the atom it belongs to decides its exclusions, its
damping and its motion; and the forms the cell's energy takes (loops, batched, the pendant
relaxation) agree with one another.
"""

import importlib.util
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from polyfind import mixedcell as MC
from polyfind import polymers as PM
from polyfind.fitting import PRESET_PARAMETERS

ROOT = Path(__file__).resolve().parents[1]
BELT4 = "pvdf-dft-valence-belt4"
SITE = -0.296  # e on each site: what DFT's potential around three chains asks for (docs/CRYSTAL_STAGE.md 5.26)


def _load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "examples" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def tools():
    mp = pytest.MonkeyPatch()
    mp.chdir(ROOT)  # the examples read docs/predictions by relative path
    IW = _load("interchain_wall")
    yield IW
    mp.undo()


def _sited(IW):
    return PM.with_bond_sites(IW.SR.C.make_polymer("vdf1", "V"), 0.65)


def _moved(m, x, shear=(0.3, -0.2), scale=0.98):
    x = np.array(x, dtype=float)
    x[m.index("shear")] = shear
    x[0] *= scale
    return x


def test_the_polymer_with_sites_is_the_polymer_with_a_point_on_each_c_h_bond():
    plain = PM.pvdf_belt(0.7, "pbe")
    sited = PM.with_bond_sites(plain, 0.65)
    assert sited.name == plain.name + "-bs0.65" and sited.bond_length == plain.bond_length
    ch2, cf2 = sited.backbone
    assert [a.element for p in ch2.pendants for a in p.atoms] == ["H", "M", "H", "M"]
    assert cf2.pendants == plain.backbone[1].pendants          # only the named element gets one
    h, m = ch2.pendants[0].atoms
    assert m.charge == 0.0 and h.charge == plain.backbone[0].pendants[0].atoms[0].charge
    assert m.offset == pytest.approx((0.65 * h.offset[0], 0.0, 0.0))
    assert ch2.charge == plain.backbone[0].charge and ch2.sub_charge == plain.backbone[0].sub_charge
    with pytest.raises(ValueError, match="on the bond"):
        PM.bond_site_atom("H", 1.09, 0.1, 1.2)
    # a copolymer keeps its sequence, every unit sited
    co = PM.with_bond_sites(PM.copolymer(name="vv", units=[PM.monomer_of(plain, "vdf")] * 2, bond_length=plain.bond_length))
    assert len(co.sequence) == 2 and all("M" in [a.element for p in b.pendants for a in p.atoms] for mon in co.sequence for b in mon.backbone[:1])


def test_uncharged_sites_change_nothing(tools):
    """Energy, every component of the gradient, the atoms' places at changed torsions and angles,
    and the pendant relaxation: the chain with sites is the chain without them."""
    SR = tools.SR
    with PRESET_PARAMETERS[BELT4].applied():
        m0, x0, _ = SR.plain_model("beta")
        m1, x1, _ = SR.plain_model("beta", _sited(tools))
        assert m1.n == m0.n + 2 and np.array_equal(x0, x1)
        atoms = np.array([e != "M" for e in list(m1.pk.elements) * 2])
        rng = np.random.default_rng(1)
        for trial in range(3):
            x = _moved(m0, x0) if trial else np.array(x0)
            if trial == 2:  # each chain's own torsions and angles moved: the frame must not see the sites
                x[m0.index("tors", 0)] += rng.normal(0.0, 2.0, m0.index("tors", 0).size)
                x[m0.index("angs", 1)] += rng.normal(0.0, 1.0, m0.index("angs", 1).size)
            e0, g0 = m0.value_and_grad(x)
            e1, g1 = m1.value_and_grad(x)
            assert e1 == pytest.approx(e0, abs=1e-10)
            assert np.abs(g1 - g0).max() < 1e-9
            assert np.abs(m1.build(x)[3][atoms] - m0.build(x)[3]).max() < 1e-12
        t0, t1 = MC.tilt_relaxation(m0, x0), MC.tilt_relaxation(m1, x1)
        assert t1[0] == pytest.approx(t0[0], abs=1e-9) and np.abs(np.asarray(t1[1]) - np.asarray(t0[1])).max() < 1e-6


def test_a_stored_frame_without_sites_places_the_chain_with_them(tools):
    """The cells polished before the sites existed carry their chains' frames on the old sites:
    the chain with sites aligns to the same frame, its atoms where they were."""
    DK = tools.IF.DK
    import json

    d = json.loads((ROOT / "docs" / "predictions" / "polished_margins_belt4.json").read_text())["references"]
    key = next(k for k in d if k.startswith("vdf5-rvdf1-6 "))
    cell = {k: (np.asarray(v, dtype=float) if isinstance(v, list) else v) for k, v in d[key].items()}
    assert cell.get("frame") is not None
    row, letters = "vdf5-rvdf1-6", DK.sequence("vdf5-rvdf1-6")[0]
    plain = DK.C.make_polymer(row, letters)
    ch0 = DK.reference_chain(plain, cell)
    ch1 = DK.reference_chain(PM.with_bond_sites(plain), cell)
    atoms = np.array([e != "M" for e in ch1.elements])
    assert atoms.sum() == len(ch0.elements) and np.abs(ch1.coords[atoms] - ch0.coords).max() < 1e-10
    # and each site is on its bond, 0.65 of the way out
    el = list(ch1.elements)
    for i in [k for k, e in enumerate(el) if e == "M"]:
        c = min((k for k, e in enumerate(el) if e == "C"), key=lambda k: np.linalg.norm(ch1.coords[k] - ch1.coords[i - 1]))
        assert np.allclose(ch1.coords[i], ch1.coords[c] + 0.65 * (ch1.coords[i - 1] - ch1.coords[c]), atol=1e-10)


def test_the_parameter_set_charges_the_sites_from_their_atoms(tools):
    SR = tools.SR
    base = PRESET_PARAMETERS[BELT4]
    params = replace(base, bond_sites=(("H", SITE),))
    assert "site on each C-H bond" in params.describe() and "site" not in base.describe()
    assert params.polymer(SR.C.make_polymer("vdf1", "V")).name.endswith("-bs0.65")
    assert base.polymer(PM.PVDF) is PM.PVDF
    plain = SR.C.make_polymer("vdf1", "V")          # built outside any block: no sites
    with base.applied():
        q0 = np.asarray(SR.plain_model("beta")[0].pk._q_cell)[:10]
    with params.applied():
        m, x, _ = SR.plain_model("beta", _sited(tools))
        q = np.asarray(m.pk._q_cell)[:m.n]
        el = list(m.pk.elements)
        assert abs(q.sum()) < 1e-12
        assert all(q[i] == pytest.approx(SITE) for i, e in enumerate(el) if e == "M")
        # each hydrogen with its site is the hydrogen there was; nothing else moved
        h = [i for i, e in enumerate(el) if e == "H"]
        assert all(q[i] + q[i + 1] == pytest.approx(q0[1]) for i in h)
        assert np.allclose(q[[i for i, e in enumerate(el) if e in ("C", "F", "X")]], q0[[0, 3, 4, 5, 6, 7, 8, 9]])
        # a chain without sites is refused, not given the charges without the term
        with pytest.raises(ValueError, match="with_bond_sites"):
            SR.plain_model("beta", plain)
        # while a chain the screen builds inside the block carries them without being asked
        assert SR.plain_model("beta")[0].n == m.n


def test_a_site_is_excluded_and_scaled_as_its_atom_is(tools):
    SR = tools.SR
    with replace(PRESET_PARAMETERS[BELT4], bond_sites=(("H", SITE),)).applied():
        m, _, _ = SR.plain_model("beta", _sited(tools))
        el = list(m.pk.elements)
        S = m.pk._scale_column_np(2)          # one (n, n) table per image along the chain, -2 .. 2
        for k in range(S.shape[0]):
            for i in [j for j, e in enumerate(el) if e == "M"]:
                rest = [j for j in range(len(el)) if j not in (i, i - 1)]
                assert np.array_equal(S[k][i, rest], S[k][i - 1, rest]) and np.array_equal(S[k][rest, i], S[k][rest, i - 1])
        own = S[2]
        for i in [j for j, e in enumerate(el) if e == "M"]:
            assert own[i, i - 1] == 0.0 and own[i - 1, i] == 0.0       # the two never see each other
        # no Lennard-Jones with anything, and its Thole length is its hydrogen's
        A, B = np.asarray(m.pk._A_nn), np.asarray(m.pk._B_nn)
        sites = np.array([e == "M" for e in el])
        assert not A[sites].any() and not A[:, sites].any() and not B[sites].any()
        i = el.index("M")
        assert np.array_equal(m.pk._thole_scale[i], m.pk._thole_scale[i - 1])
        assert m.pk._alpha_cell[i] == 0.0


def test_charged_sites_are_point_charges_in_the_lattice_sum(tools):
    """The fixed-charge Ewald sum of the cell against a plain Ewald sum written apart from it
    (`examples/interchain_wall.py`), as differences between rigid cells; point charges, so the
    penetration damping is set aside."""
    SR, IW = tools.SR, tools
    params = replace(PRESET_PARAMETERS[BELT4], bond_sites=(("H", SITE),), charge_damping=0.0)
    with params.applied():
        m, x0, _ = SR.plain_model("beta", _sited(tools))
        pk = m.pk
        idx, q = pk._es_idx, np.asarray(pk._q_cell, dtype=float)
        assert len(idx) == 2 * 8      # the belt sites stay out; the bond sites are in

        def both(x):
            latn, _, _, Pn, _ = m.build(x)
            return float(pk._ewald.terms(Pn[idx], q[idx], latn).total), 0.5 * IW.COULOMB * q[idx] @ IW.ewald_green(Pn[idx], latn) @ q[idx]

        a = both(x0)
        for shear, scale in (((0.6, 0.0), 1.0), ((0.5, -0.4), 0.96), ((0.0, 0.0), 1.05)):
            b = both(_moved(m, x0, shear, scale))
            assert b[0] - a[0] == pytest.approx(b[1] - a[1], abs=2e-6)
            assert abs(b[0] - a[0]) > 0.05


def test_the_cell_with_charged_sites_has_the_gradient_of_its_energy_in_both_forms(tools):
    """Induced dipoles, damping and the buffered wall all on: finite differences of the energy,
    the batched evaluation against the loops, and no dipole on a site."""
    SR = tools.SR
    params = replace(PRESET_PARAMETERS[BELT4], bond_sites=(("H", SITE),), wall=("b14-7", 0.2))
    with params.applied():
        m, x0, _ = SR.plain_model("beta", _sited(tools))
        x = _moved(m, x0)
        e, g = m.value_and_grad(x)
        for k in range(len(x)):
            h = 1e-5
            xp, xm = x.copy(), x.copy()
            xp[k] += h
            xm[k] -= h
            fd = (m.value_and_grad(xp)[0] - m.value_and_grad(xm)[0]) / (2 * h)
            assert fd == pytest.approx(g[k], abs=2e-6 * max(1.0, abs(g[k])))
        latn, frames, _, Pn, _ = m.build(x)
        E0, g0, l0 = MC.placed_energy_and_grad(m.pk, frames, Pn, latn)
        E1, g1, l1 = MC.placed_energy_and_grad(m.pk, frames, Pn, latn, xp=np)
        assert E1 == pytest.approx(E0, abs=1e-9) and np.abs(g1 - g0).max() < 1e-9 and np.abs(l1 - l0).max() < 1e-9
        _, p, _, _ = m.pk._polarize(Pn, np.asarray(m.pk._q_cell, dtype=float), latn, m.flips.astype(float))
        el = np.array(list(m.pk.elements) * 2)
        assert not p[el == "M"].any() and np.abs(p[el == "C"]).max() > 1e-3
        # the pendant relaxation runs with the sites following their bonds
        dE, tilts, _ = MC.tilt_relaxation(m, x0)
        assert dE <= 1e-9 and np.isfinite(np.asarray(tilts)).all()


def test_charge_drawn_along_the_backbone_tells_a_cf2_by_its_neighbours(tools):
    """`FFParameters.backbone_pull`: each fluorine draws charge onto its carbon from that carbon's
    backbone neighbours.  In plain PVDF every CF2 has two CH2 neighbours and gains twice what one
    bond carries; where two CF2 meet they draw equally on each other and that bond carries
    nothing; the chain stays neutral and nothing but backbone carbons changes."""
    import json

    IF = tools.IF
    base = PRESET_PARAMETERS[BELT4]
    pull = 0.0155
    params = replace(base, backbone_pull=(("F", pull),))
    assert "draws 0.0155 e" in params.describe()
    table = json.loads((ROOT / "docs" / "predictions" / "dilute_kink_belt4.json").read_text())
    assert "all_trans_references" in table

    def charges(prm, row, letters):
        with prm.applied():
            polymer = IF.C.make_polymer(row, letters)
            ch = IF.DK.repeat_chains(polymer, "TT", np.full((1, polymer.bonds_per_repeat), 180.0),
                                     np.array([polymer.backbone[k].backbone_angle for k in range(polymer.bonds_per_repeat)], dtype=float))[0]
            from polyfind import pack as pack_mod
            pk = pack_mod.CrystalPacker(ch, n_chains=2, **IF.DK.KW)
            return list(pk.elements), np.asarray(pk._q_cell)[:pk.n].copy()

    for row, letters, expect in (("vdf1", "V", [+2, -2]), ("vdf1-rvdf1-2", "VI", [+1, -1, -1, +1]), ("vdf5-rvdf1-6", "VVVVIV", None)):
        el, q0 = charges(base, row, letters)
        _, q = charges(params, row, letters)
        carbons = [i for i, e in enumerate(el) if e == "C"]
        rest = [i for i in range(len(el)) if i not in carbons]
        assert abs(q.sum()) < 1e-12 and np.array_equal(q[rest], q0[rest])
        moved = (q - q0)[carbons] / (2 * pull)   # in units of what one CH2-CF2 bond carries, onto the CH2 carbon
        if expect is not None:
            assert np.allclose(moved, expect)
        else:  # one reversed unit in six: the CF2-CF2 and CH2-CH2 bonds carry nothing
            assert sorted(np.round(moved, 6)) == sorted([2.0] * 4 + [1.0] * 2 + [-2.0] * 4 + [-1.0] * 2)
        # the same thing as the trial's pattern (examples/interchain_wall.py charge_patterns), at d_CC = -2 pull
        pat = tools.charge_patterns(el, q0, 1)
        assert np.allclose(q - q0, -2 * pull * pat[2])


def test_the_conformer_design_carries_the_two_charge_terms():
    """`ValenceDesign.lattice` with bond sites and a backbone pull: the conformers see the
    lattice's charges, so the torsions can be fitted again under them.  The Coulomb energy is a
    plain sum over sites written out here; the forces are the energy's; and without the two terms
    nothing changes."""
    spec = importlib.util.spec_from_file_location("tf", ROOT / "tests" / "test_fitting.py")
    TF = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(TF)
    from polyfind import fitting as F
    from polyfind.forcefield import COULOMB, Frame

    frames = TF._synthetic_frames()
    x = TF._valence_x()
    belt4 = PRESET_PARAMETERS[BELT4]
    pull, site, frac = 0.0155, SITE, 0.65
    prm = replace(belt4, bond_sites=(("H", site),), backbone_pull=(("F", pull),))
    design = F.ValenceDesign(frames)
    design.lattice = belt4
    e_plain = design.evaluate(design.unpack(x))[0]
    design.lattice = prm
    p = design.unpack(x)
    assert (p["offset"], p["site_q"], p["pull"]) == (0.0, site, pull) and p["project"] == dict(belt4.charge_offsets)["F"]
    e, forces, _ = design.evaluate(p)
    assert np.abs(e - e_plain).max() > 0.1

    # the Coulomb energy of one frame, as a plain sum over the sites
    f = frames[1]                                   # PVDF with gauche bonds
    el, c, n = list(f.elements), np.asarray(f.coords, dtype=float), f.n_atoms
    nb = [[] for _ in range(n)]
    for a, b in f.bonds:
        nb[a].append(b)
        nb[b].append(a)
    inc = {(a, b): v for a, b, v in belt4.charge_increments}
    q = np.zeros(n)
    for a, b in f.bonds:
        for i, j in ((a, b), (b, a)):               # the first of an increment's pair gains it, the second loses it
            q[i] += inc.get((el[i], el[j]), 0.0) - inc.get((el[j], el[i]), 0.0)
    for a in range(n):                              # each fluorine draws onto its carbon from the carbon's carbon neighbours
        if el[a] == "C":
            n_f = sum(el[b] == "F" for b in nb[a])
            for b in nb[a]:
                if el[b] == "C":
                    q[a] -= pull * n_f
                    q[b] += pull * n_f
    t = dict(belt4.charge_offsets)["F"]
    for a in range(n):                              # the off-nucleus charge on the nucleus, the bond's dipole kept
        if el[a] == "F":
            moved = q[a] * t / F.PROJECT_BOND
            q[a] += moved
            q[nb[a][0]] -= moved
    assert abs(q.sum()) < 1e-12
    pos, charge, host = [c[i] for i in range(n)], list(q), list(range(n))
    for h in range(n):                              # a site on every C-H bond, charged from its hydrogen
        if el[h] == "H" and el[nb[h][0]] == "C":
            charge[h] -= site
            pos.append(c[nb[h][0]] + frac * (c[h] - c[nb[h][0]]))
            charge.append(site)
            host.append(h)
    dist = f.graph_distances()
    coulomb = 0.0
    for i in range(len(pos)):
        for j in range(i + 1, len(pos)):
            gd = dist[host[i], host[j]]
            if host[i] == host[j] or gd < 3:        # a site is excluded exactly where its hydrogen is
                continue
            coulomb += (p["scale14"] if gd == 3 else 1.0) * COULOMB * charge[i] * charge[j] / np.linalg.norm(pos[i] - pos[j])  # the 1-4 scale is the vector's
    # the design's Coulomb energy alone: the same parameters with every other term switched off
    off = dict(p, lj_d=np.zeros_like(p["lj_d"]), lj_pairs={}, torsion=np.zeros_like(p["torsion"]),
               bond_k=np.zeros_like(p["bond_k"]), angle_k=np.zeros_like(p["angle_k"]))
    assert design.evaluate(off)[0][1] == pytest.approx(coulomb, abs=1e-9)

    # the forces are the energy's: atoms of both kinds moved, a sited hydrogen and its carbon among them
    first = int(np.flatnonzero(design.atom_frame == 1)[0])
    h_atom = next(i for i in range(n) if el[i] == "H")
    for atom, k in ((h_atom, 0), (nb[h_atom][0], 1), (next(i for i in range(n) if el[i] == "F"), 2)):
        def energy_of(dx):
            cc = c.copy()
            cc[atom, k] += dx
            moved_design = F.ValenceDesign([frames[0], Frame(**{**f.__dict__, "coords": cc})] + frames[2:])
            moved_design.lattice = prm
            return moved_design.evaluate(moved_design.unpack(x))[0][1]
        fd = -(energy_of(1e-5) - energy_of(-1e-5)) / 2e-5
        assert fd == pytest.approx(forces[first + atom, k], rel=1e-5, abs=1e-5)
    # a design built for sites elsewhere on the bond refuses these parameters
    other = F.ValenceDesign(frames[:2], bond_site_fraction=0.5)
    other.lattice = prm
    with pytest.raises(ValueError, match="of the way out"):
        other.unpack(x)


def test_a_polymer_built_inside_a_sets_block_is_built_as_the_set_needs_it(tools):
    """`fitting.active_parameters` and the screen's `make_polymer`: inside the block of a set with
    bond sites the row's chain carries them; outside, or inside a set without, it does not."""
    from polyfind.fitting import active_parameters

    C = tools.SR.C
    base = PRESET_PARAMETERS[BELT4]
    sited = replace(base, bond_sites=(("H", SITE),))
    assert active_parameters() is None and C.SITES_FOR is None
    assert not C.make_polymer("vdf1", "V").name.endswith("-bs0.65")
    with base.applied():
        assert active_parameters() is base and not C.make_polymer("vdf1", "V").name.endswith("-bs0.65")
    with sited.applied():
        assert active_parameters() is sited
        assert C.make_polymer("vdf1", "V").name.endswith("-bs0.65")
        with base.applied():          # blocks do not nest: the outer set's charges are the ones a packer gets
            assert active_parameters() is sited
    assert active_parameters() is None
