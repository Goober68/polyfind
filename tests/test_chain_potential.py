"""The electrostatic potential around one chain, polyfind's sources against DFT's: the tools and the stored result.

`examples/chain_potential.py` samples DFT's potential around a chain alone in a box and fits the
strengths of point sources to it; `examples/interchain_wall.py` carries the same sources into the
energies of rigid cells.  The sums are held to known answers here, and the stored fit
(`docs/predictions/chain_potential_belt4.json`) to the stored samples it was made from.
"""

import importlib.util
import json
from itertools import product
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
PRED = ROOT / "docs" / "predictions"


def _load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "examples" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def CP(monkeypatch_module):
    monkeypatch_module.chdir(ROOT)  # the examples read docs/predictions by relative path
    return _load("chain_potential")


@pytest.fixture(scope="module")
def monkeypatch_module():
    mp = pytest.MonkeyPatch()
    yield mp
    mp.undo()


def _quadrupoles(rng, n=4):
    """Neutral sources with no dipole: +1, -2, +1 on a line, at random places and directions."""
    pos, q = [], []
    for _ in range(n):
        c, u = rng.uniform(2.0, 5.0, 3), rng.normal(size=3)
        u *= 0.4 / np.linalg.norm(u)
        pos += [c + u, c, c - u]
        q += [1.0, -2.0, 1.0]
    return np.array(pos), np.array(q)


def test_the_periodic_potential_does_not_depend_on_the_splitting_and_is_the_image_sum(CP):
    rng = np.random.default_rng(3)
    cell = np.array([7.0, 8.0, 6.0])
    pos, q = _quadrupoles(rng)
    pts = rng.uniform(0.0, 1.0, (40, 3)) * cell
    a = CP.ewald_potential(pts, pos, q, cell, alpha=0.45)
    b = CP.ewald_potential(pts, pos, q, cell, alpha=0.6)
    assert np.abs(a - b).max() < 1e-8
    # quadrupoles: the plain sum over images converges absolutely, to the same function up to one constant
    direct = np.zeros(len(pts))
    R = 9
    for i, j, k in product(range(-R, R + 1), repeat=3):
        d = pts[:, None, :] - (pos + np.array([i, j, k]) * cell)[None, :, :]
        direct += (q[None, :] / np.sqrt((d * d).sum(axis=-1))).sum(axis=1)
    far = np.array([np.linalg.norm((p - pos) - cell * np.round((p - pos) / cell), axis=1).min() > 0.8 for p in pts])
    diff = (a - direct)[far]
    assert diff.std() < 2e-4 * a[far].std()


def test_the_pair_matrix_gives_one_energy_however_the_lattice_is_written(CP):
    IW = CP.IW
    rng = np.random.default_rng(5)
    pos, q = _quadrupoles(rng)
    lat = np.array([[7.0, 0.0, 0.0], [1.0, 8.0, 0.0], [0.5, -0.7, 6.0]])
    e = 0.5 * q @ IW.ewald_green(pos, lat) @ q
    assert e == pytest.approx(0.5 * q @ IW.ewald_green(pos, lat, alpha=0.6) @ q, abs=1e-8)
    # the same lattice on another basis, and the same sources moved by lattice vectors
    other = np.array([lat[0] + lat[2], lat[1] - lat[0], lat[2]])
    assert e == pytest.approx(0.5 * q @ IW.ewald_green(pos, other) @ q, abs=1e-8)
    moved = pos + np.outer(rng.integers(-2, 3, len(pos)), lat[0]) + np.outer(rng.integers(-2, 3, len(pos)), lat[2])
    assert e == pytest.approx(0.5 * q @ IW.ewald_green(moved, lat) @ q, abs=1e-8)
    # and its rows are the potential: at one source, that of the others
    cell = np.array([7.0, 8.0, 6.0])
    G = IW.ewald_green(pos, np.diag(cell))
    for i in (0, 4):
        others = np.delete(np.arange(len(q)), i)
        neutral = q[others] - q[others].sum() / len(others)  # the potential is defined for a neutral set
        probe = CP.ewald_potential(pos[i:i + 1], pos[others], neutral, cell)[0]
        assert probe == pytest.approx(G[i, others] @ neutral, abs=1e-7)


def test_the_charge_patterns_split_the_charges_and_count_unlike_backbone_bonds(CP):
    IW = CP.IW
    plain = ["C", "H", "H", "C", "F", "X", "X", "F", "X", "X"]
    q = np.array([-0.3, 0.15, 0.15, 0.4, -0.2, 0.0, 0.0, -0.2, 0.0, 0.0])
    pat = IW.charge_patterns(plain, q, 2)
    assert pat.shape == (3, 20)
    assert np.allclose(pat[0] + pat[1], np.tile(q, 2))
    assert np.allclose(pat[0, :10], [-0.3, 0.15, 0.15, 0, 0, 0, 0, 0, 0, 0])
    # head to tail: every carbon has two unlike neighbours (the chain closes on its image)
    assert np.allclose(pat[2, :10], [-2, 0, 0, 2, 0, 0, 0, 0, 0, 0])
    # every second unit reversed, CH2-CF2-CF2-CH2: one unlike neighbour each
    reversed_ = plain + ["C", "F", "X", "X", "F", "X", "X", "C", "H", "H"]
    q2 = np.concatenate([q, [0.4, -0.2, 0, 0, -0.2, 0, 0, -0.3, 0.15, 0.15]])
    pat2 = IW.charge_patterns(reversed_, q2, 1)
    assert np.allclose(pat2[2][[0, 3, 10, 17]], [-1, 1, 1, -1])
    assert abs(pat2[2].sum()) < 1e-12 and np.allclose(pat2[0] + pat2[1], q2)


def test_the_extra_sources_are_what_their_names_say(CP):
    IW = CP.IW
    el = ["C", "H", "H", "C", "F", "X", "X", "F", "X", "X"]
    sites = np.array([[0.0, 0.0, 0.0], [0.9, 0.6, 0.0], [-0.9, 0.6, 0.0], [0.0, -0.8, 1.3], [1.1, -1.6, 1.3], [1.6, -1.1, 1.3], [0.6, -2.1, 1.3],
                      [-1.1, -1.6, 1.3], [-0.6, -2.1, 1.3], [-1.6, -1.1, 1.3]])
    src = IW.extra_sources(el, sites)
    assert set(src) == set(IW.EXTRA)
    for name, (pos, q) in src.items():
        assert abs(q.sum()) < 1e-12, name
    dipole = lambda name: (src[name][1][:, None] * src[name][0]).sum(axis=0)  # noqa: E731
    u_h = [(sites[i] - sites[0]) / np.linalg.norm(sites[i] - sites[0]) for i in (1, 2)]
    u_f = [(sites[i] - sites[3]) / np.linalg.norm(sites[i] - sites[3]) for i in (4, 7)]
    assert np.allclose(dipole("mu_H"), u_h[0] + u_h[1])          # 1 e A on each hydrogen, away from its carbon
    assert np.allclose(dipole("mu_F"), u_f[0] + u_f[1])
    for name in ("Q_H", "Q_F", "X_belt"):
        assert np.allclose(dipole(name), 0.0, atol=1e-12), name
    # the bond's own site: -1 on it and +1 on the hydrogen, the site BOND_SITE of the way out
    bond = np.linalg.norm(sites[1] - sites[0])
    assert np.allclose(dipole("M_H"), (1.0 - IW.BOND_SITE) * bond * (u_h[0] + u_h[1]))
    assert len(src["M_H"][1]) == 4 and len(src["X_belt"][1]) == 6


def _stored():
    return json.loads((PRED / "chain_potential_belt4.json").read_text())


def test_what_the_stored_fit_says():
    t = _stored()
    c = t["chains"]
    assert set(c) == {"ht", "hh", "tgtg"}
    # polyfind's potential has about the same absolute error around all three chains ...
    assert all(0.35 < c[n]["rms_error"] < 0.45 for n in c)
    # ... which is a few per cent of the potential of a chain with a dipole and half of it for the chain without one
    assert c["ht"]["rms_error"] / c["ht"]["rms_dft"] < 0.04 and c["tgtg"]["rms_error"] / c["tgtg"]["rms_dft"] < 0.06
    assert 0.4 < c["hh"]["rms_error"] / c["hh"]["rms_dft"] < 0.5
    # the induced dipoles are part of what makes the plain chain right
    assert c["ht"]["charges_alone_rms_error"] > 8 * c["ht"]["rms_error"]
    fits = t["fits"]
    charges = fits["together, the three charge patterns"]["rms_error"]
    assert charges["hh"] > 0.39                       # no choice of atom charges repairs the reversed chain
    for term in ("X_belt", "Q_F", "mu_F"):            # nor does any one term on fluorine
        assert fits[f"together, the three charge patterns and {term}"]["rms_error"]["hh"] > 0.39, term
    for term in ("Q_H", "mu_H", "M_H"):               # a term on hydrogen halves the error on every chain
        e = fits[f"together, the three charge patterns and {term}"]["rms_error"]
        assert e["ht"] < 0.19 and e["hh"] < 0.23 and e["tgtg"] < 0.15, term
    m = fits["together, the three charge patterns and M_H"]["multipliers"]
    assert m["M_H"] == pytest.approx(0.296, abs=2e-3) and m["s_CH"] == pytest.approx(0.267, abs=2e-3) and m["s_CF"] == pytest.approx(0.771, abs=2e-3)


def test_the_stored_fit_is_the_stored_samples(CP):
    """One chain recomputed from the samples: polyfind as it is, and with the fitted bond site."""
    t = _stored()
    data = CP.stored_blocks(PRED / "chain_potential_pbe_samples.npz", ["hh"], 2.2)
    assert len(data["hh"]["v"]) == t["chains"]["hh"]["points"]
    blocks = CP.potential_blocks(data)
    assert CP.rms_error(blocks["hh"], CP.START) == pytest.approx(t["chains"]["hh"]["rms_error"], abs=2e-4)
    fit = t["fits"]["together, the three charge patterns and M_H"]
    assert CP.rms_error(blocks["hh"], fit["multipliers"]) == pytest.approx(fit["rms_error"]["hh"], abs=2e-4)
    # the samples keep to the shell they were asked for
    assert data["hh"]["near"].min() >= 2.2 - 1e-5 and data["hh"]["near"].max() <= 5.0 + 1e-5
