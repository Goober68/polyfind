"""The energy between rigid chains as the cell is strained and sheared: polyfind against DFT, stored results.

Reads the stored files only (`examples/interchain_frames.py` writes them; the DFT is
`examples/mlip_frames.py label`).  Everything inside a chain is the same in every frame of a
packing, so the numbers here are energy between chains and nothing else.
"""

import json
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
PRED = ROOT / "docs" / "predictions"
EV_TO_KCAL = 23.060548


def _table():
    return json.loads((PRED / "interchain_belt4.json").read_text())["packings"]


def _frames():
    """{frame: (atoms, family, DFT energy per monomer in eV)} of the stored labelled frames, without ase."""
    lines = (PRED / "mlip" / "interchain_frames_pbe_d3bj.extxyz").read_text().splitlines()
    out, i = {}, 0
    while i < len(lines):
        n = int(lines[i])
        head = lines[i + 1]
        value = lambda key: head.split(f"{key}=")[1].split()[0]  # noqa: E731
        out[value("frame")] = (n, value("family"), float(value("REF_energy")) / int(value("n_monomers")))
        i += n + 2
    return out


def _row(block, scale, shear=(0.0, 0.0)):
    return next(r for r in block["frames"] if r["scale"] == list(scale) and r["shear_monomers"] == list(shear))


def test_the_stored_frames_are_the_cells_of_the_table():
    frames, t = _frames(), _table()
    assert {f[1] for f in frames.values()} == {"interchain"}
    rows = {r["frame"]: (pack, r) for pack, b in t.items() for r in b["frames"]}
    assert set(rows) == set(frames)
    assert sum(pack == "beta" for pack, _ in rows.values()) == 27
    assert sum(pack == "alpha" for pack, _ in rows.values()) == 27
    assert sum(pack == "hh" for pack, _ in rows.values()) == 37
    atoms = {"beta": 12, "candidate": 12, "alpha": 24, "hh": 24, "cand6": 72}
    for name, (pack, r) in rows.items():
        assert frames[name][0] == atoms[pack] == 6 * t[pack]["monomers_per_cell"]
        zero = frames[f"inter_{pack}_00"][2]
        assert r["dft"] == pytest.approx((frames[name][2] - zero) * EV_TO_KCAL, abs=2e-6)
        # polyfind's energy is its electrostatics and its repulsion and dispersion, and nothing else moves
        assert r["polyfind"] == pytest.approx(r["polyfind_electrostatics"] + r["polyfind_repulsion_dispersion"], abs=1e-9)


def test_plain_pvdf_is_right_stretched_and_too_stiff_squeezed():
    """The beta cell under belt4: within a tenth of DFT wherever a lateral edge is 2% or more
    longer, 1.1 to 1.45 times DFT where one is 4% or more shorter."""
    beta = _table()["beta"]
    for scale in ((1.02, 1.0), (1.04, 1.0), (1.06, 1.0), (1.0, 1.04), (1.0, 1.06), (1.02, 1.02), (1.04, 1.04)):
        r = _row(beta, scale)
        assert 0.90 < r["polyfind"] / r["dft"] < 1.05, scale
    for scale in ((0.94, 1.0), (0.96, 1.0), (1.0, 0.94), (1.0, 0.96), (0.96, 0.96)):
        r = _row(beta, scale)
        assert 1.10 < r["polyfind"] / r["dft"] < 1.45, scale


def test_the_registry_cost_is_overstated_more_the_closer_the_chains():
    """Half a monomer along ``a`` in the beta cell, each cell against its own unsheared state: 2.4
    times DFT in a cell 4% wider, 5.4 times in one 4% narrower.  The wall is too steep close in,
    which no 12-6 coefficients repair."""
    beta = _table()["beta"]
    wide, narrow = _row(beta, (1.04, 1.04), (0.5, 0.0)), _row(beta, (0.96, 0.96), (0.5, 0.0))
    # each against its own unsheared cell
    cost = lambda r, s: (r["polyfind"] - _row(beta, s)["polyfind"]) / (r["dft"] - _row(beta, s)["dft"])  # noqa: E731
    assert 1.8 < cost(wide, (1.04, 1.04)) < 2.6
    assert 4.5 < cost(narrow, (0.96, 0.96)) < 6.0
    assert all(_row(beta, (s, s), sh)["polyfind"] > _row(beta, (s, s), sh)["dft"] for s in (0.96, 1.04) for sh in ((0.25, 0.0), (0.5, 0.0), (0.0, 0.5)))


def test_the_reversed_chain_is_missed_the_other_way():
    """Every second unit reversed: fluorine against fluorine and hydrogen against hydrogen between
    chains.  DFT's cell wants to be wider along ``a`` than polyfind's (its energy falls on
    stretching, polyfind's rises) and resists squeezing four times as hard; and half a monomer
    of shear one way along ``a`` costs three times what it costs the other way, which polyfind
    halves."""
    hh = _table()["hh"]
    assert _row(hh, (1.04, 1.0))["dft"] < -0.15 and _row(hh, (1.04, 1.0))["polyfind"] > 0.04
    assert 3.5 < _row(hh, (0.94, 1.0))["dft"] / _row(hh, (0.94, 1.0))["polyfind"] < 4.6
    there, back = _row(hh, (1.0, 1.0), (0.5, 0.0)), _row(hh, (1.0, 1.0), (1.5, 0.0))  # 1.5 is minus a half: the chain repeats every two monomers
    assert 3.0 < back["dft"] / there["dft"] < 3.5 and 1.9 < back["polyfind"] / there["polyfind"] < 2.2
    err = np.array([r["polyfind"] - r["dft"] for r in hh["frames"][1:]])
    dft = np.array([r["dft"] for r in hh["frames"][1:]])
    assert 0.33 < np.sqrt((err ** 2).mean()) < 0.40 and 0.80 < np.sqrt((dft ** 2).mean()) < 0.90
