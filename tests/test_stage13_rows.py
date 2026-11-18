"""The rows stage 13 of the crystal fit takes from outside the module (`examples/stage13_rows.py`), on the stored files."""

import importlib.util
import json
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from polyfind.fitting import PRESET_PARAMETERS

ROOT = Path(__file__).resolve().parents[1]
PRED = ROOT / "docs" / "predictions"
BELT4 = "pvdf-dft-valence-belt4"


@pytest.fixture(scope="module")
def rows():
    mp = pytest.MonkeyPatch()
    mp.chdir(ROOT)  # the examples read docs/predictions by relative path
    spec = importlib.util.spec_from_file_location("stage13_rows", ROOT / "examples" / "stage13_rows.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    yield mod, mod.Rows(chains=("hh",), conformers=False)
    mp.undo()


def test_the_rigid_rows_under_belt4_are_the_stored_table(rows):
    """Every rigid cell evaluated in the kernels under the preset the table was made with gives the
    table's own number, and DFT's beside it; a packing's rows count as a fixed number of rows."""
    mod, r = rows
    table = json.loads((PRED / "interchain_belt4.json").read_text())["packings"]
    out = r.rigid(PRESET_PARAMETERS[BELT4])
    stored = {f["frame"][6:]: (pack, f) for pack, b in table.items() for f in b["frames"][1:]}
    assert set(out) == {f"rigid {k}" for k in stored} and len(out) == 117
    for name, (pack, f) in stored.items():
        e, dft, z = out[f"rigid {name}"]
        assert e == pytest.approx(f["polyfind"], abs=2e-7) and dft == f["dft"]
        n = len(table[pack]["frames"]) - 1
        assert z == pytest.approx(np.sqrt(mod.COUNT[pack] / n) * (e - dft) / np.hypot(0.03, 0.15 * dft), abs=1e-9)
    # a set with bond sites is evaluated on chains that carry them: uncharged, the same numbers
    sited = r.rigid(replace(PRESET_PARAMETERS[BELT4], bond_sites=(("H", 0.0),)))
    assert max(abs(sited[k][0] - out[k][0]) for k in out) < 1e-8
    charged = r.rigid(replace(PRESET_PARAMETERS[BELT4], bond_sites=(("H", -0.318),)))
    assert max(abs(charged[k][0] - out[k][0]) for k in out) > 0.05
    assert "rigid hh" in mod.summary(out) and "rigid cand6" in mod.summary(out)


def test_the_potential_row_is_linear_in_the_charges(rows):
    """At belt4's own charges the row is the stored comparison's error; the site and the pull enter
    as the potential's fourth and third patterns."""
    mod, r = rows
    base = PRESET_PARAMETERS[BELT4]
    stored = json.loads((PRED / "chain_potential_belt4.json").read_text())["chains"]["hh"]
    row = r.potential(base)["potential around hh"]
    assert row[0] == pytest.approx(stored["rms_error"], abs=2e-4) and row[2] == pytest.approx(row[0] / mod.POTENTIAL_TOL)
    inc = {(a, b): d for a, b, d in base.charge_increments}
    scaled = tuple((a, b, d * {("C", "H"): 0.300, ("C", "F"): 0.701}.get((a, b), 1.0)) for a, b, d in base.charge_increments)
    better = r.potential(replace(base, charge_increments=scaled, bond_sites=(("H", -0.318),), backbone_pull=(("F", 0.0155),)))["potential around hh"]
    assert better[0] == pytest.approx(0.210, abs=0.003) and inc[("C", "H")] < 0
    assert r.conformers(base) == {}   # built without the reference set


def test_the_lone_chains_line_is_a_row_a_point_on_one_scale(rows):
    """The helix up to 6% stretched and the all-trans chain about its repeat, each against the point
    where DFT is lowest; under belt4 the helix is two to three times too stiff when pulled and
    falls when pushed, which is what the rows are there to say."""
    mod, r = rows
    line = json.loads((PRED / "chain_stretch_pbe_d3bj.json").read_text())["branches"]
    out = r.lone(PRESET_PARAMETERS[BELT4])
    helix = [q for q in line["helix"]["points"] if q["repeat"] <= mod.LINE_HELIX_MAX and q["energy"] > 0]
    assert set(out) == {f"lone chain, helix at {q['repeat']:.3f} A" for q in helix} | {f"lone chain, trans at {q['repeat']:.3f} A" for q in line["trans"]["points"]}
    assert len(out) == 11
    for q in helix:
        e, dft, z = out[f"lone chain, helix at {q['repeat']:.3f} A"]
        assert dft == q["energy"] and z == pytest.approx(np.sqrt(mod.LINE_COUNT / 11) * (e - dft) / np.hypot(0.05, 0.05 * dft), abs=1e-9)
    pulled, pushed = out["lone chain, helix at 4.850 A"], out["lone chain, helix at 4.550 A"]
    assert 1.6 < pulled[0] / pulled[1] < 2.4 and pushed[0] < 0.0 < pushed[1]
    # the all-trans chain alone against the helix alone: belt4 has it 0.3 too high
    low = out["lone chain, trans at 2.591 A"]
    assert low[1] == pytest.approx(3.449, abs=0.005) and low[0] == pytest.approx(3.78, abs=0.05)


def test_the_candidates_rows_on_frozen_cells_are_the_relaxed_split_at_the_set_that_froze_them(rows):
    """Cells relaxed under belt4 and frozen, evaluated under belt4: the rows are the split of
    examples/chain_split.py itself (cell +0.064 and +0.120, inside a chain +1.71 and +0.27), DFT's
    beside them; under another set the same cells give other numbers without being relaxed."""
    mod, _ = rows
    r = mod.Rows(chains=("hh",), conformers=False, line=None, reversed_chains=None, frozen=PRED / "chain_split_frozen_belt4.json")
    out = r.candidate(PRESET_PARAMETERS[BELT4])
    want = {"candidate TTTT minus TG'TG": (0.064, 0.070), "candidate TTTT minus TG'TG, inside a chain": (1.708, 1.060),
            "candidate GTG'T minus TG'TG": (0.120, 0.197), "candidate GTG'T minus TG'TG, inside a chain": (0.270, 0.415)}
    assert set(out) == set(want)
    for k, (mine, dft) in want.items():
        assert out[k][0] == pytest.approx(mine, abs=0.002) and out[k][1] == pytest.approx(dft, abs=0.001)
    assert out["candidate TTTT minus TG'TG"][2] == pytest.approx((out["candidate TTTT minus TG'TG"][0] - out["candidate TTTT minus TG'TG"][1]) / 0.05)
    assert mod.Rows(chains=("hh",), conformers=False, line=None, reversed_chains=None).candidate(PRESET_PARAMETERS[BELT4]) == {}
