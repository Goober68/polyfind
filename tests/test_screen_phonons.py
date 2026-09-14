"""examples/screen_phonons.py: the tied antipolar polish never leaves the antipolar subspace, the
gamma-free search only ever lowers the gamma = 90 answer, and the summary reads the record it is
given."""
import importlib.util
from pathlib import Path

import numpy as np
import pytest

from polyfind.fitting import antipolar_offsets
from polyfind.pack import CrystalPacker, default_bounds, periodic_chain
from polyfind.polymers import PVDF, THREE_STATE

SCRIPT = Path(__file__).resolve().parents[1] / "examples" / "screen_phonons.py"
SPEC = importlib.util.spec_from_file_location("screen_phonons", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
SP = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(SP)

T = 0


@pytest.fixture(scope="module")
def packer():
    return CrystalPacker(periodic_chain(PVDF, [T, T], THREE_STATE), n_chains=2)


def _lohi(packer, gamma_free):
    b = default_bounds(packer.chain, gamma_free=gamma_free)
    return (np.array([b[k][0] for k in SP.CELL_KEYS]), np.array([b[k][1] for k in SP.CELL_KEYS]))


def _antipolar_start(packer):
    flip, dphi = antipolar_offsets(packer)[0]
    return np.array([4.6, 9.9, 90.0, 20.0, 20.0 + dphi, 1.0, float(flip)])


def test_tied_polish_stays_exactly_antipolar_and_lowers_the_energy(packer):
    p0 = _antipolar_start(packer)
    e0 = float(packer.energy(p0[None])[0])
    assert packer.result(p0).polarization_magnitude < 1e-9, "the start itself must be antipolar"
    for gamma_free in (False, True):
        lo, hi = _lohi(packer, gamma_free)
        p, e = SP.polish_antipolar(packer, p0, lo, hi, gamma_free=gamma_free)
        assert e <= e0 + 1e-9
        assert e == pytest.approx(float(packer.energy(p[None])[0]))
        assert packer.result(p).polarization_magnitude < 1e-9
        assert (p[4] - p[3]) % 360.0 == pytest.approx((p0[4] - p0[3]) % 360.0, abs=1e-9)
        assert p[6] == p0[6]
        assert lo[0] - 1e-9 <= p[0] <= hi[0] + 1e-9 and lo[1] - 1e-9 <= p[1] <= hi[1] + 1e-9
        if gamma_free:
            assert lo[2] - 1e-9 <= p[2] <= hi[2] + 1e-9
        else:
            assert p[2] == 90.0


def test_a_start_outside_the_bounds_is_not_projected_uphill(packer):
    """antipolar_cell_exact may return a cell just outside default_bounds; the polish must not
    start by clipping it (which raised CDFE's antipolar energy by 0.15 kcal/mol per monomer)."""
    lo, hi = _lohi(packer, False)
    p0 = _antipolar_start(packer)
    p0[0] = lo[0] - 0.15  # below the screen's lower bound on a
    e0 = float(packer.energy(p0[None])[0])
    p, e = SP.polish_antipolar(packer, p0, lo, hi, gamma_free=False)
    assert e < e0 - 1e-6, "a start outside the bounds with a downhill gradient must move, not sit on the bound"
    assert packer.result(p).polarization_magnitude < 1e-9


def test_gamma_free_search_never_returns_above_the_gamma_90_answer(packer):
    lo90, hi90 = _lohi(packer, False)
    lof, hif = _lohi(packer, True)
    p90, e90 = SP.polish_antipolar(packer, _antipolar_start(packer), lo90, hi90, gamma_free=False)
    n_mon = packer.chain.n_monomers * packer.n_chains
    p, e_mon, seeds = SP.antipolar_gamma_free(packer, p90, lof, hif)
    assert e_mon <= e90 / n_mon + 1e-9
    assert packer.result(p).polarization_magnitude < 1e-9
    assert len(seeds) == 2 * (1 + len(SP.GAMMA_SEEDS))
    assert min(s["e_per_monomer"] for s in seeds) == pytest.approx(e_mon, abs=1e-4)


def test_closest_contact_is_below_the_chain_separation_and_symmetric(packer):
    p = np.array([4.6, 8.6, 90.0, 0.0, 0.0, 0.4, 0.0])
    d = SP.closest_contact(packer, p)
    assert 0.5 < d < np.hypot(4.6 / 2, 8.6 / 2)  # atoms sit off-axis, so closer than the axes are
    wide = p.copy()
    wide[:2] = 40.0
    assert SP.closest_contact(packer, wide) > d  # pulled apart, the contact opens


def test_tied_polish_evaluates_the_offset_inside_the_repeat(packer):
    """A start with dz far outside [0, c) must be judged on the wrapped offset, where the
    packer's energy is valid, and come back inside it."""
    lo, hi = _lohi(packer, True)
    p0 = _antipolar_start(packer)
    c = float(packer.chain.c)
    far = p0.copy()
    far[5] = p0[5] + 9 * c
    p_in, e_in = SP.polish_antipolar(packer, p0, lo, hi, gamma_free=True)
    p_far, e_far = SP.polish_antipolar(packer, far, lo, hi, gamma_free=True)
    assert 0.0 <= p_far[5] < c
    assert e_far == pytest.approx(e_in, abs=1e-6)
    assert e_far == pytest.approx(float(packer.energy(p_far[None])[0]), abs=1e-9)


def test_summary_reads_the_record():
    cell = lambda e, low, slid: {"atoms": 12, "relaxed": {"energy": e, "lowest_optical": low, "max_displacement": slid}}  # noqa: E731
    rec = {"polymer": "x", "screened": True, "n_monomers_per_cell": 2,
           "screen": {"gamma90": {"antipolar_gap": 1.0, "arrangement": "polar"},
                      "gamma_free": {"antipolar_gap": 0.1, "arrangement": "not resolved"}},
           "cells": {"polar90": cell(-10.0, [30.0, 40.0, 50.0], 0.05), "polar_free": cell(-10.0, [30.0, 40.0, 50.0], 0.05),
                     "anti90": cell(-8.0, [10.0, 20.0, 30.0], 1.0), "anti_free": cell(-9.0, [35.0, 45.0, 55.0], 0.1)}}
    assert SP.relaxed_gap(rec, "polar_free", "anti_free") == pytest.approx(0.5)
    assert SP.relaxed_gap({"cells": {}}, "polar_free", "anti_free") is None
    table = SP.summary_table([rec, {"polymer": "y", "screened": False, "reason": "no repeat"}])
    lines = table.splitlines()
    assert len(lines) == 3
    assert "+1.000" in lines[1] and "+0.100" in lines[1] and "+0.500" in lines[1] and "35.0 45.0 55.0" in lines[1]
    assert lines[2].startswith("y") and "no repeat" in lines[2]
