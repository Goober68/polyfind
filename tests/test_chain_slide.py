"""examples/chain_slide.py: one chain of beta-PVDF held at a slip along its axis, the lattice around it relaxed."""
import importlib.util
import json
from pathlib import Path

import numpy as np

from polyfind.fitting import PRESET_PARAMETERS

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("chain_slide", ROOT / "examples" / "chain_slide.py")
CS = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CS)
DATA = json.load(open(ROOT / "docs" / "predictions" / "chain_slide_belt4.json"))
E5 = json.load(open(ROOT / "docs" / "predictions" / "e5_polyfind_belt4.json"))
TEN = json.load(open(ROOT / "docs" / "predictions" / "dilute_kink_belt4.json"))["all_trans_references"]["vdf10"]
PRESET = PRESET_PARAMETERS["pvdf-dft-valence-belt4"]


def _cell(repeat: int, nx: int, ny: int):
    """``(model, x, monomer period, slip index)`` of a tiling of the stored short-repeat cell (call inside the preset)."""
    polymer, _, d = CS.reference("V" * repeat, DATA["references"][f"vdf{repeat}"])
    model, x, ch = CS.DK.make_model(polymer, d, nx, ny)
    return model, x, float(d["state"][3]) / repeat, model.index("z", 0)[0]


def _slid(model, x, iz, u):
    y = x.copy()
    y[iz] += u
    return y


def test_a_short_repeat_is_the_ten_monomer_crystal():
    """A one- or two-monomer repeat gives the energy per monomer of the polished ten-monomer cell and the
    rigid numbers registered for E5: half a monomer of one sublattice (two-chain cell) and of one chain
    in four.  So the large lateral supercells can be built from one monomer per chain."""
    with PRESET.applied():
        for repeat in (1, 2):
            model, x, cm, iz = _cell(repeat, 1, 1)
            assert abs(model.physical_energy(x) / (2 * repeat) - float(TEN["energy"]) / 20) < 1e-5
            sub = (model.physical_energy(_slid(model, x, iz, 0.5 * cm)) - model.physical_energy(x)) / repeat
            assert abs(sub - E5["own_cell"]["a"]["slides"]["c/2"]["per_moved_monomer_kcal"]) < 1e-3
            model, x, cm, iz = _cell(repeat, 2, 1)
            one = (model.physical_energy(_slid(model, x, iz, 0.5 * cm)) - model.physical_energy(x)) / repeat
            assert abs(one - E5["own_cell"]["c"]["slides"]["c/2"]["per_moved_monomer_kcal"]) < 1e-3
            assert abs(one - DATA["rows"][f"vdf{repeat} 2x1"]["rigid"]) < 1e-4


def test_the_half_slip_is_stationary_along_the_slide_and_nowhere_else():
    """The all-trans chain has a mirror plane normal to its axis: at half a monomer the axial force on
    every chain vanishes, which is why the relaxed half-slip is a saddle and needs no path; at a
    quarter the slid chain is pushed back and its neighbours forward by as much in total."""
    with PRESET.applied():
        model, x, cm, iz = _cell(1, 2, 1)
        zs = model.index("z")
        _, g = model.value_and_grad(_slid(model, x, iz, 0.5 * cm))
        assert np.abs(g[zs]).max() < 1e-3
        _, g = model.value_and_grad(_slid(model, x, iz, 0.25 * cm))
        assert g[zs][0] > 4.0 and abs(g[zs].sum()) < 1e-6


def test_the_variables_of_each_level():
    """The slip and the repeat are never free; the anchor chain keeps its lateral position; drag frees
    the other chains' axial shifts but the anchor's; the cell option frees a, b and gamma only."""
    with PRESET.applied():
        model, x, cm, iz = _cell(1, 4, 2)
        dist = CS.distances_from(model, x, 0)
        anchor = int(np.argmax(dist))
        near = np.sort(dist[1:])
        assert np.allclose(near[:2], x[0] / 4, atol=1e-6) and np.allclose(near[2:6], np.hypot(x[0] / 8, x[1] / 4), atol=1e-6)
        turn, shift = CS.free_set(model, "turn", anchor), CS.free_set(model, "shift", anchor)
        assert set(turn) == set(model.index("phi")) and set(turn) < set(shift)
        assert not set(shift) & set(model.index("xy", anchor)) and not set(shift) & set(model.index("z")) and not set(shift) & {0, 1, 2, 3}
        drag = CS.free_set(model, "shift", anchor, drag=True)
        assert set(drag) - set(shift) == set(model.index("z")) - {iz, model.index("z", anchor)[0]}
        assert set(CS.free_set(model, "shift", anchor, cell=True)) - set(shift) == {0, 1, 2}


def test_a_start_on_the_symmetric_lattice_stays_on_it():
    """In the four-chain cell the relaxation started exactly on the half-slipped lattice stops at a
    stationary point held by symmetry (6.63 kcal/mol per monomer); a perturbed start finds the state
    below it (6.07), where the chains have turned by up to 11.5 degrees.  Turning alone gains nothing
    in this cell.  The stored row is reproduced."""
    row = DATA["rows"]["vdf1 2x1"]
    with PRESET.applied():
        model, x, cm, iz = _cell(1, 2, 1)
        anchor = int(np.argmax(CS.distances_from(model, x, 0)))
        xs = _slid(model, x, iz, 0.5 * cm)
        _, turn = CS.held_slip(model, x, xs, "turn", anchor, [0], 300, 1, lambda s: None)
        xb, shift = CS.held_slip(model, x, xs, "shift", anchor, [0], 300, 1, lambda s: None)
    assert abs(turn["barrier"] - row["rigid"]) < 1e-3
    assert abs(shift["starts"]["as started"] - 6.6306) < 2e-3 and abs(shift["barrier"] - 6.0737) < 2e-3
    assert abs(shift["barrier"] - row["shift"]["barrier"]) < 2e-3
    assert abs(shift["largest_turn_deg"] - 11.5) < 0.2 and abs(shift["axial_force_on_slid_chain"]) < 1e-3


def test_the_lattice_giving_way_halves_the_barrier():
    """Sixteen chains: with the neighbours free to turn and to move aside the half-slip costs 3.70
    kcal/mol per monomer against 6.80 in a rigid lattice.  The two neighbours along the short axis move
    by 0.17 A, the four diagonal ones by 0.14 A and turn by 5 degrees; no chain feels an axial force, so
    the state is the saddle of the slide in this cell."""
    row = DATA["rows"]["vdf1 4x2"]
    with PRESET.applied():
        model, x, cm, iz = _cell(1, 4, 2)
        anchor = int(np.argmax(CS.distances_from(model, x, 0)))
        xb, r = CS.held_slip(model, x, _slid(model, x, iz, 0.5 * cm), "shift", anchor, [], 300, 1, lambda s: None)
    assert abs(r["barrier"] - row["shift"]["barrier"]) < 2e-3 and r["barrier"] < 0.56 * row["rigid"]
    shells = list(r["shells_A"].values())
    assert shells[0]["chains"] == 2 and abs(shells[0]["largest_move_A"] - 0.166) < 0.01
    assert shells[1]["chains"] == 4 and abs(shells[1]["largest_move_A"] - 0.137) < 0.01 and abs(shells[1]["largest_turn_deg"] - 5.0) < 0.3
    assert abs(r["axial_force_on_slid_chain"]) < 1e-3 and r["largest_axial_force_on_another"] < 1e-3


def test_the_two_bounds_meet():
    """The stored convergence: at fixed cell the barrier falls with the number of chains, with the
    lateral cell free it rises, and the two series extrapolate (linearly in one over the number of
    chains, from their two largest cells) to the same dilute limit within 0.03 kcal/mol per monomer."""
    rows = {int(r["chains"]): r for k, r in DATA["rows"].items() if k.startswith("vdf1 ") and k.split()[1] in ("2x1", "4x2", "6x3", "8x4")}
    fixed = {n: r["shift"]["barrier"] for n, r in rows.items() if "shift" in r}
    free = {n: r["shift_cell_free"]["barrier"] for n, r in rows.items() if "shift_cell_free" in r}
    assert all(fixed[a] > fixed[b] for a, b in zip(sorted(fixed), sorted(fixed)[1:]))
    assert all(free[a] < free[b] for a, b in zip(sorted(free), sorted(free)[1:]))

    def limit(series):
        n1, n2 = sorted(series)[-2:]
        slope = (series[n1] - series[n2]) / (1 / n1 - 1 / n2)
        return series[n2] - slope / n2

    assert max(free.values()) < limit(fixed) + 0.03 and min(fixed.values()) > limit(free) - 0.03
    assert abs(limit(fixed) - limit(free)) < 0.03
    assert 0.44 < limit(fixed) / rows[64]["rigid"] < 0.48
