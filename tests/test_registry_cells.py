"""Other registries of a two-chain cell (examples/registry_cells.py; docs/MLIP_SECOND_OPINION.md 3.24): the builder, the stored starts and the choice."""

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
MLIP = ROOT / "docs" / "predictions" / "mlip"
FILM, MATRIX, LEADER, NITRILE = "vdf9-trfe4-cfe1-14", "vdf10-trfe4-14", "vdf9-trfe3-an1-rvdf1-14", "vdf9-trfe4-an1-14"
CELLS = {"film": (f"row_{FILM}_TTTT_chain1_in_all_s1", f"row_{FILM}_TGTGp_chain1_in_all_s1"), "matrix": (f"row_{MATRIX}_TTTT_chain1_in_all_s1", f"row_{MATRIX}_TGTGp_chain1_in_all_s1"),
         "leader": ("tt_leader_from_rows14an_r1", "helix_leader_TGTGp_r1"), "nitrile": ("tt_nitrile_from_rows14an_r1", f"row_{NITRILE}_GpTGT_cfe_in_all_s1")}


def _load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "examples" / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def labelled():
    from ase.io import read

    return {str(a.info["frame"]): a for a in read(MLIP / "cells14_frames_pbe_d3bj.extxyz", index=":")}


def _bonds(atoms):
    """Every bond's length, sorted (the pairs examples/registry_cells.py takes for bonds)."""
    from ase.neighborlist import neighbor_list

    return np.sort(neighbor_list("d", atoms, _load("registry_cells").BONDS))


def _same_places(a, b, cell):
    """Two sets of positions that differ by whole cell vectors only."""
    d = np.linalg.solve(cell.T, (b - a).T).T
    return np.abs(d - np.round(d)).max() < 1e-8


def test_a_sheared_cell_has_its_chains_whole_where_they_were_and_its_neighbours_along(labelled):
    tool, path = _load("registry_cells"), _load("mlip_path")
    for flat, helix in CELLS.values():
        for name, (sa, sb) in ((flat, (1, -1)), (helix, (-2, 2))):
            a = labelled[name]
            b = tool.sheared(a, sa, sb)
            assert np.array(path._shear(b)) - np.array(path._shear(a)) == pytest.approx((sa, sb), abs=1e-9)
            assert b.get_volume() == pytest.approx(a.get_volume(), rel=1e-12) and b.cell.array[2] == pytest.approx(a.cell.array[2])
            assert _bonds(b) == pytest.approx(_bonds(a), abs=1e-9) and len(tool.chains(b)[0]) == 2          # no bond cut or made
            # every atom is where it was, up to whole cell vectors of the old cell and whole monomers along the chains for those carried across a face
            monomer = a.cell.array[2] / 14
            d = np.linalg.solve(np.vstack([a.cell.array[0], a.cell.array[1], monomer]).T, (b.positions - a.positions).T).T
            assert np.abs(d - np.round(d)).max() < 1e-8
    a = labelled[CELLS["film"][0]]
    assert _same_places(a.positions, tool.sheared(a, 0, 0).positions, a.cell.array)


def test_a_turned_cell_has_every_chain_turned_about_its_own_axis(labelled):
    tool = _load("registry_cells")
    for flat, _ in CELLS.values():
        a = labelled[flat]
        before = tool.dipole_from_a(a)
        assert 70 < before < 90                                                    # the upright cells: dipoles along the long edge
        for angle in (117, -117):
            b = tool.turned(a, angle)
            turn = (tool.dipole_from_a(b) - before - angle + 180) % 360 - 180
            assert abs(turn) < 0.5 and _bonds(b) == pytest.approx(_bonds(a), abs=1e-9)
            assert b.cell.array == pytest.approx(a.cell.array)
        assert _same_places(a.positions, tool.turned(tool.turned(a, 117), -117).positions, a.cell.array)   # and back
    # the film composition's low cell of 3.23 has its dipoles 21 degrees from the short edge and its neighbours a monomer along
    from ase.io import read

    low = read(MLIP / "paths" / "all_trans_tight.extxyz", index=0)
    assert abs(abs(tool.dipole_from_a(low)) - 180) == pytest.approx(21, abs=1.0)


def test_a_start_whose_chains_touch_is_given_room(labelled):
    tool, path = _load("registry_cells"), _load("mlip_path")
    a = labelled[CELLS["leader"][0]]
    tight = tool.sheared(tool.turned(a, 117), 1, -1)
    assert tool.closest_between_chains(tight) < 1.3                                # the nitrile groups meet
    wide = tool.spread(tight, 1.1)
    assert tool.closest_between_chains(wide) > tool.closest_between_chains(tight) + 0.1
    assert _bonds(wide) == pytest.approx(_bonds(tight), abs=1e-9) and wide.get_volume() == pytest.approx(1.21 * tight.get_volume(), rel=1e-9)
    along = lambda x: [float(x.cell.array[i] @ x.cell.array[2]) for i in (0, 1)]
    assert along(wide) == pytest.approx(along(tight), abs=1e-9)                     # what the edges carry along the chains is kept


def test_the_stored_starts_are_the_scripts_and_none_has_chains_that_touch(labelled):
    from ase.io import read

    tool, path = _load("registry_cells"), _load("mlip_path")
    starts = read(MLIP / "registry_starts.extxyz", index=":")
    names = [str(a.info["frame"]) for a in starts]
    assert len(names) == 34 and len(set(names)) == 34 and names[:3] == ["reg_matrix_tt_p1m1_tp117", "reg_leader_tt_p1m1_tm117", "reg_nitrile_tt_p1m1_tm117"]
    assert sum("_tt_" in n for n in names) == 18 and sum("_hx_" in n for n in names) == 16
    assert "reg_film_tt_p1m1_tp117" not in names                                   # 3.23's cell: labelled already
    for a in starts:
        label, conformation = str(a.info["registry_label"]), str(a.info["registry_conformation"])
        ref = labelled[str(a.info["registry_of"])]
        assert str(a.info["registry_of"]) == CELLS[label][0 if conformation == "tt" else 1]
        sa, sb = tool._ints(a.info["registry_shear"])
        turn = int(a.info["registry_turn"])
        assert str(a.info["frame"]) == tool.name_of(label, conformation, sa, sb, turn)
        assert turn == (0 if conformation == "hx" else turn) and turn in (0, 117, -117) and (turn != 117 or label in ("film", "matrix")) and (turn != -117 or label in ("leader", "nitrile"))
        assert np.array(path._shear(a)) - np.array(path._shear(ref)) == pytest.approx((sa, sb), abs=1e-6)
        assert len(a) == len(ref) and int(a.info["n_monomers"]) == 28 and list(np.atleast_1d(a.info["kmesh"])) == [4, 2, 1]
        assert tool.gauche(a, path) == tool.gauche(ref, path) and (tool.gauche(a, path) == 0 if conformation == "tt" else tool.gauche(a, path) >= 24)
        close = tool.closest_between_chains(a)
        assert close is None or close > 1.9
        assert 1.0 <= float(a.info["registry_widened"]) <= 1.08
    queue = (ROOT / "examples" / "gpu_jobs" / "queue_20261122_registry.txt").read_text().splitlines()
    first = next(line for line in queue if line and not line.startswith("#")).split()
    assert first[0] == "registry_card.sh" and (ROOT / first[1]).samefile(MLIP / "registry_starts.extxyz") and [int(v) for v in first[2:5]] == [4, 400, 2700]
    assert first[6] in {str(a.info["frame"]) for a in read(ROOT / first[5], index=":")}
    for line in queue:
        if line and not line.startswith("#"):
            assert (ROOT / "examples" / "gpu_jobs" / line.split()[0]).exists() or line.split()[0] == "stop"


def test_the_choice_is_the_lowest_cell_at_rest_of_each_composition_and_conformation():
    tool = _load("registry_cells")
    row = lambda label, conformation, energy, converged=True, kept=True: {"relaxed": True, "frame": f"{label}_{conformation}_{energy}", "label": label, "conformation": conformation,
                                                                           "energy": energy, "converged": converged, "conformation_kept": kept}
    rows = {"a": row("matrix", "tt", -0.4), "b": row("matrix", "tt", -0.9), "c": row("matrix", "tt", -2.0, converged=False), "d": row("matrix", "hx", 0.2),
            "e": row("matrix", "hx", -3.0, kept=False), "f": row("leader", "tt", 0.5), "g": row("film", "tt", -1.1), "h": row("film", "tt", -0.6), "i": {"relaxed": False},
            "j": row("nitrile", "hx", -0.1)}
    picked = tool.chosen(rows, 0.3, {"film:tt": -1.0})
    assert picked == [("film_tt_-1.1", -1.1), ("matrix_tt_-0.9", -0.9), ("nitrile_hx_-0.1", -0.1), ("matrix_hx_0.2", 0.2)]
    assert tool.chosen(rows, 0.3, {"film:tt": -1.2})[0] == ("matrix_tt_-0.9", -0.9)
    assert tool.name_of("film", "tt", 1, -1, 117) == "reg_film_tt_p1m1_tp117" and tool.name_of("leader", "hx", -2, 2, 0) == "reg_leader_hx_m2p2_tz0"
    assert tool._ints("1 -1") == [1, -1] and tool._ints(np.array([2, -2])) == [2, -2]


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
