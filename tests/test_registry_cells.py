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
    assert first[0] == "registry_card.sh" and (ROOT / first[1]).samefile(MLIP / "registry_starts.extxyz") and [int(v) for v in first[2:5]] == [4, 600, 2700]
    assert first[6] in {str(a.info["frame"]) for a in read(ROOT / first[5], index=":")}
    more = [line.split() for line in queue if line.startswith("registry_more.sh ")]
    assert (ROOT / more[0][1]).samefile(MLIP / "registry_starts.extxyz") and more[0][2:] == ["4", "900", "2400", "auto:4.0"]
    assert (ROOT / more[1][1]).samefile(MLIP / "registry_starts_like.extxyz") and more[1][2:] == ["4", "900", "1200", "auto:6.0"]
    jobs = [line.split()[0] for line in queue if line and not line.startswith("#")]
    assert jobs[:4] == ["registry_card.sh", "registry_more.sh", "registry_card.sh", "registry_more.sh"]
    assert (ROOT / "examples" / "gpu_jobs" / "registry_more.cap").read_text().split() == ["02", "600", "04", "60"]     # the second job cut to ten minutes, the fourth to one
    for line in queue:
        if line and not line.startswith("#"):
            assert (ROOT / "examples" / "gpu_jobs" / line.split()[0]).exists() or line.split()[0] == "stop"


def test_a_cell_made_like_the_film_compositions_low_cell_has_its_lattice_turn_and_registry(labelled):
    from ase.io import read

    tool, path = _load("registry_cells"), _load("mlip_path")
    low = read(MLIP / "paths" / "all_trans_tight.extxyz", index=0)
    assert str(low.info["frame"]) == "film_r_released_tight"
    upright = labelled[CELLS["film"][0]]
    shape = lambda a: (*np.linalg.norm(tool.lateral(a)[0], axis=1), float(np.degrees(np.arccos(tool.lateral(a)[0][0] @ tool.lateral(a)[0][1] / np.prod(np.linalg.norm(tool.lateral(a)[0], axis=1))))))
    mine, turn = tool.like(upright, upright, low)
    assert shape(mine) == pytest.approx(shape(low), abs=0.06) and shape(low) == pytest.approx((5.03, 9.20, 93.0), abs=0.03)      # the film composition's own: the low cell's lattice
    assert path._shear(mine) == pytest.approx(path._shear(low), abs=1e-6) and (turn + 180) % 360 - 180 == pytest.approx(117, abs=1.0)
    assert abs((tool.dipole_from_a(mine) - tool.dipole_from_a(low) + 180) % 360 - 180) < 1.0 and _bonds(mine) == pytest.approx(_bonds(upright), abs=1e-9)
    other, _ = tool.like(labelled[CELLS["matrix"][0]], upright, low)
    ratio = np.array(shape(other)[:2]) / np.array(shape(labelled[CELLS["matrix"][0]])[:2])
    assert ratio == pytest.approx(np.array(shape(low)[:2]) / np.array(shape(upright)[:2]), abs=0.005)                            # the same strain, of its own lattice
    starts = read(MLIP / "registry_starts_like.extxyz", index=":")
    assert [str(a.info["frame"]) for a in starts] == [f"reg_{c}_tt_like_film" for c in ("matrix", "film", "nitrile", "leader")]
    for a in starts:
        label = str(a.info["registry_label"])
        assert str(a.info["registry_of"]) == CELLS[label][0] and tool._ints(a.info["registry_shear"]) == [1, -1] and int(a.info["registry_turn"]) == 118
        assert tool.gauche(a, path) == 0 and float(a.info["registry_widened"]) == (1.2 if label in ("nitrile", "leader") else 1.0)   # the nitrile groups meet in this setting
    queue = [line.split() for line in (ROOT / "examples" / "gpu_jobs" / "queue_20261122_registry.txt").read_text().splitlines() if line.startswith("registry_card.sh ")]
    assert len(queue) >= 2 and (ROOT / queue[1][1]).samefile(MLIP / "registry_starts_like.extxyz") and queue[1][2:] == ["4", "600", "900"]


def test_the_matrix_is_the_film_compositions_low_cell_with_its_chlorine_made_fluorine(labelled):
    from ase.io import read

    tool, path = _load("registry_cells"), _load("mlip_path")
    low = read(MLIP / "paths" / "all_trans_tight.extxyz", index=0)
    film, matrix = labelled[CELLS["film"][0]], labelled[CELLS["matrix"][0]]
    # a chain as the cycle of its backbone carbons' substituents, each with a hand: the low cell's chains are the film composition's own
    mine, theirs = tool.substituents(low), tool.substituents(film)
    assert all(len(c) == 28 for c in mine) and sum(k == "ClF" for k, _ in mine[0]) == 1 and sum(h != 0 for _, h in mine[0]) == 5     # one CFCl and four CHF carbons a chain
    assert all(any(tool.same_molecule(m, t) for t in theirs) for m in mine)
    assert not tool.same_molecule(mine[0], [(k, -h) for k, h in mine[0]])                                                    # its mirror image is another chain
    swapped = tool.swapped(low, "Cl", "F", 1.36)
    assert sorted(swapped.get_chemical_symbols()) == sorted(matrix.get_chemical_symbols()) and "Cl" not in swapped.get_chemical_symbols()
    assert all(any(tool.same_molecule(m, t) for t in tool.substituents(matrix)) for m in tool.substituents(swapped))             # the matrix's own chains, sequence and hands
    moved = np.linalg.norm(swapped.positions - low.positions, axis=1)
    d = np.linalg.solve(low.cell.array.T, (swapped.positions - low.positions).T).T
    assert sum(np.abs(row - np.round(row)).max() > 1e-8 for row in d) == 2                                                     # two atoms moved: the two chlorines
    starts = read(MLIP / "registry_starts_as_film.extxyz", index=":")
    names = [str(a.info["frame"]) for a in starts]
    assert names == ["reg_matrix_tt_as_film_low", "reg_film_tt_low_registry_back", "reg_film_tt_low_chains_in_upright", "reg_matrix_tt_as_film_low_registry_back"]
    expected = {names[0]: ("matrix", [1, -1]), names[1]: ("film", [0, 0]), names[2]: ("film", [0, 0]), names[3]: ("matrix", [0, 0])}
    for a in starts:
        label, shear = expected[str(a.info["frame"])]
        ref = labelled[CELLS[label][0]]
        assert str(a.info["registry_of"]) == CELLS[label][0] and tool._ints(a.info["registry_shear"]) == shear and tool.gauche(a, path) == 0
        assert all(any(tool.same_molecule(m, t) for t in tool.substituents(ref)) for m in tool.substituents(a))
        assert np.array(path._shear(a)) - np.array(path._shear(ref)) == pytest.approx(shear, abs=0.07)
    # the chains of the low cell are kept in all four: the pattern of the twists off trans is the low cell's, not the upright cell's
    twists = lambda a: np.where(path.torsions_deg(a, path.backbone_torsions(a)) > 0, -180.0, 180.0) + path.torsions_deg(a, path.backbone_torsions(a))
    assert np.abs(twists(low)).max() < 25 and np.abs(twists(film)).max() < 25
    for a in starts:
        assert np.abs(twists(a) - twists(low)).max() < 1.0 and np.abs(twists(a) - twists(film)).max() > 15
    lines = [l.split() for l in (ROOT / "examples" / "gpu_jobs" / "queue_20261122_registry.txt").read_text().splitlines() if "registry_starts_as_film" in l and not l.startswith("#")]
    assert lines[0][0] == "registry_card.sh" and lines[0][2:] == ["4", "1000", "1500"] and lines[1][1].endswith("registry_starts_as_film_nitrile.extxyz") and lines[1][2:] == ["2", "1000", "1200"]
    # the nitrile cell likewise: the CFE unit made an acrylonitrile unit, the nitrile where the chlorine was; the other hand is another molecule
    nitrile = labelled[CELLS["nitrile"][0]]
    made = tool.nitrile_for(low, "Cl", "F")
    assert sorted(made.get_chemical_symbols()) == sorted(nitrile.get_chemical_symbols()) and len(made) == 170
    assert all(any(tool.same_molecule(m, t) for t in tool.substituents(nitrile)) for m in tool.substituents(made))
    assert not all(any(tool.same_molecule(m, t) for t in tool.substituents(nitrile)) for m in tool.substituents(tool.nitrile_for(low, "F", "Cl")))
    from ase.neighborlist import neighbor_list
    i, j, d = neighbor_list("ijd", made, 1.3)
    cn = sorted(float(x) for a, b, x in zip(i, j, d) if made.symbols[a] == "N" and made.symbols[b] == "C")
    assert cn == pytest.approx([1.16, 1.16], abs=1e-6)
    two = read(MLIP / "registry_starts_as_film_nitrile.extxyz", index=":")
    assert [str(a.info["frame"]) for a in two] == ["reg_nitrile_tt_as_film_low", "reg_nitrile_tt_as_film_low_registry_back"]
    assert [tool._ints(a.info["registry_shear"]) for a in two] == [[1, -1], [0, 0]] and all(str(a.info["registry_of"]) == CELLS["nitrile"][0] and tool.gauche(a, path) == 0 for a in two)


def test_chains_with_the_upright_cells_twists_rest_higher_in_every_other_registry_tried(labelled):
    """The first 34 starts on the fourth pod's card: 21 relaxed, 18 at rest, none below the cell it was made from."""
    from ase.io import read

    table = json.loads((ROOT / "docs" / "predictions" / "registry_shifted.json").read_text())["rows"]
    names = [str(a.info["frame"]) for a in read(MLIP / "registry_starts.extxyz", index=":")]
    assert list(table) == names
    done = {k: v for k, v in table.items() if v["relaxed"]}
    rest = {k: v for k, v in done.items() if v["converged"]}
    assert (len(done), len(rest)) == (21, 18) and all(v["conformation_kept"] for v in done.values())
    flat = sorted(v["energy"] for v in rest.values() if v["conformation"] == "tt")
    helix = sorted(v["energy"] for v in rest.values() if v["conformation"] == "hx")
    assert (len(flat), len(helix)) == (11, 7) and (flat[0], flat[-1]) == pytest.approx((0.45, 2.30), abs=0.006) and (helix[0], helix[-1]) == pytest.approx((1.40, 5.75), abs=0.006)
    assert min(v["energy"] for v in done.values()) > 0.44                                   # the three not at rest are above too, so far
    low = {c: min(v["energy"] for v in rest.values() if v["label"] == c and v["conformation"] == "tt") for c in ("film", "matrix", "leader", "nitrile")}
    assert low == pytest.approx({"film": 0.62, "matrix": 1.23, "leader": 2.30, "nitrile": 0.45}, abs=0.006)
    # the film composition's own low registry made from upright chains (the other hand of it was run here) rests 1.86 above, where 3.23's cell is 1.27 below
    assert rest["reg_film_tt_m1p1_tp117"]["energy"] == pytest.approx(1.86, abs=0.006) and rest["reg_film_tt_m1p1_tp117"]["shear_along_chains_monomers"] == pytest.approx([-0.98, 1.03], abs=0.006)
    # each cell kept the registry it was given, within a third of a monomer, and the stored frames are the table's
    for v in done.values():
        got = np.array(v["shear_along_chains_monomers"]) - np.array(v["shear_along_chains_monomers_labelled"])
        assert np.abs(got - np.array(v["shear_asked"])).max() < 0.75
    frames = read(MLIP / "registry_shifted_relaxed.extxyz", index=":")
    assert [str(a.info["frame"]) for a in frames] == [v["frame"] for v in done.values()]
    for a in frames:
        ref = labelled[str(a.info["registry_of"])]
        assert (float(a.info["potential_eV_per_monomer"]) - float(ref.info["potential_eV_per_monomer"])) * 14 * 23.0605 == pytest.approx(table[str(a.info["frame"])[:-len("_in_all_s1")]]["energy"], abs=1e-6)


def test_upright_chains_in_the_low_cells_lattice_and_registry_do_not_find_the_low_cell():
    """The control: the film composition's cell made like its own low cell, from upright chains, rests 0.57 above upright."""
    table = json.loads((ROOT / "docs" / "predictions" / "registry_like_film.json").read_text())["rows"]
    control = table["reg_film_tt_like_film"]
    assert control["converged"] and control["energy"] == pytest.approx(0.57, abs=0.006) and control["shear_along_chains_monomers"] == pytest.approx([1.03, -1.04], abs=0.006)
    assert abs(abs(control["dipole_from_a_deg"]) - 159) < 3                                   # the low cell's registry and turn, kept
    low = json.loads((ROOT / "docs" / "predictions" / "kink_paths_film.json").read_text())["points"]
    below = low["all-trans, the sheet's registry, relaxed on"]["potential"] - low["all-trans, upright"]["potential"]
    assert below == pytest.approx(-1.27, abs=0.006) and control["energy"] - below == pytest.approx(1.84, abs=0.01)   # 1.84 above the cell it was made like
    assert table["reg_matrix_tt_like_film"]["converged"] and table["reg_matrix_tt_like_film"]["energy"] == pytest.approx(1.40, abs=0.006)
    assert not table["reg_nitrile_tt_like_film"]["converged"] and not table["reg_leader_tt_like_film"]["converged"]    # the nitrile groups meet in this setting: cut, not at rest


def test_one_compositions_chains_made_anothers_carbon_by_carbon(labelled):
    from ase.io import read

    tool, path = _load("registry_cells"), _load("mlip_path")
    low = read(MLIP / "paths" / "all_trans_tight.extxyz", index=0)
    counts = {}
    for label in ("matrix", "nitrile", "leader"):
        target = labelled[CELLS[label][0]]
        made, counts[label] = tool.morphed(low, target)
        assert sorted(made.get_chemical_symbols()) == sorted(target.get_chemical_symbols())
        assert all(any(tool.same_molecule(m, t) for t in tool.substituents(target)) for m in tool.substituents(made)) and tool.gauche(made, path) == 0
    assert counts == {"matrix": 2, "nitrile": 2, "leader": 8}          # a chain: the CFE carbon; and for the leader the reversed unit's two carbons and a TrFE's as well
    # the same cells as the swaps made: the nitrile's to the atom, the matrix's but for the other fluorine's bond length on the changed carbon
    nitrile = read(MLIP / "registry_starts_as_film_nitrile.extxyz", index=0)
    assert _same_places(nitrile.positions, tool.morphed(low, labelled[CELLS["nitrile"][0]])[0].positions, low.cell.array)
    matrix = read(MLIP / "registry_starts_as_film.extxyz", index=0)
    assert np.abs(matrix.positions - tool.morphed(low, labelled[CELLS["matrix"][0]])[0].positions).max() < 0.03
    leader = read(MLIP / "registry_starts_as_film_leader.extxyz", index=":")
    assert [str(a.info["frame"]) for a in leader] == ["reg_leader_tt_as_film_low", "reg_leader_tt_as_film_low_registry_back"]
    assert tool._ints(leader[0].info["registry_shear"]) == [1, -1] and tool._ints(leader[1].info["registry_shear"]) == [0, 0]
    assert str(leader[0].info["registry_of"]) == CELLS["leader"][0] and leader[0].get_chemical_formula() == labelled[CELLS["leader"][0]].get_chemical_formula()
    # the leader's chain against the film composition's: a reversed VDF first, the nitrile for the chlorine, a VDF for the fourth TrFE
    word = lambda a: " ".join(k + ("+" if h > 0 else "-" if h < 0 else "") for k, h in tool.substituents(a)[0])
    assert word(labelled[CELLS["film"][0]]) == "HH FF HH FF FH+ FF HH FF HH FF FH+ FF HH FF HH ClF+ HH FF HH FF FH+ FF HH FF HH FF FH+ FF"
    assert word(labelled[CELLS["leader"][0]]) == "FF HH HH FF FH+ FF HH FF HH FF FH+ FF HH FF HH CH+ HH FF HH FF HH FF HH FF HH FF FH+ FF"


def test_the_low_cell_is_its_chains_and_its_registry_together_and_the_matrix_has_it_too(labelled):
    """The four cells made of the film composition's low cell, at rest in the potential on the fourth pod's card."""
    from ase.io import read

    tool, path = _load("registry_cells"), _load("mlip_path")
    table = json.loads((ROOT / "docs" / "predictions" / "registry_as_film.json").read_text())["rows"]
    assert all(v["converged"] and v["conformation_kept"] for v in table.values()) and len(table) == 4
    energy = {k: v["energy"] for k, v in table.items()}
    assert energy == pytest.approx({"reg_matrix_tt_as_film_low": -1.50, "reg_film_tt_low_registry_back": 0.18, "reg_film_tt_low_chains_in_upright": 0.54,
                                    "reg_matrix_tt_as_film_low_registry_back": 0.83}, abs=0.006)
    # each kept the registry it was given
    assert table["reg_matrix_tt_as_film_low"]["shear_along_chains_monomers"] == pytest.approx([1.06, -1.08], abs=0.006)
    for k in ("reg_film_tt_low_registry_back", "reg_film_tt_low_chains_in_upright", "reg_matrix_tt_as_film_low_registry_back"):
        assert np.abs(table[k]["shear_along_chains_monomers"]).max() < 0.05
    # beside them: the film composition's own low cell 1.27 below upright, and upright chains in its lattice and registry 0.57 above
    film = json.loads((ROOT / "docs" / "predictions" / "kink_paths_film.json").read_text())["points"]
    assert film["all-trans, the sheet's registry, relaxed on"]["potential"] - film["all-trans, upright"]["potential"] == pytest.approx(-1.27, abs=0.006)
    assert json.loads((ROOT / "docs" / "predictions" / "registry_like_film.json").read_text())["rows"]["reg_film_tt_like_film"]["energy"] == pytest.approx(0.57, abs=0.006)
    # the relaxed matrix cell is still the matrix's own chains, and carries the low cell's twists, not the matrix's upright cell's
    frames = {str(a.info["frame"]): a for a in read(MLIP / "registry_as_film_relaxed.extxyz", index=":")}
    matrix_low, matrix_up = frames["reg_matrix_tt_as_film_low_in_all_s1"], labelled[CELLS["matrix"][0]]
    assert all(any(tool.same_molecule(m, t) for t in tool.substituents(matrix_up)) for m in tool.substituents(matrix_low))
    low = read(MLIP / "paths" / "all_trans_tight.extxyz", index=0)
    twists = lambda a: np.where(path.torsions_deg(a, path.backbone_torsions(a)) > 0, -180.0, 180.0) + path.torsions_deg(a, path.backbone_torsions(a))
    assert np.abs(twists(matrix_low) - twists(low)).max() < 15 and np.abs(twists(matrix_low) - twists(matrix_up)).max() > 20


def test_the_nitrile_composition_has_the_low_cell_too(labelled):
    """The nitrile cell made of the film composition's low cell (the chlorine a nitrile), at rest in the potential on the fourth pod's card."""
    from ase.io import read

    tool, path = _load("registry_cells"), _load("mlip_path")
    table = json.loads((ROOT / "docs" / "predictions" / "registry_as_film_nitrile.json").read_text())["rows"]
    assert all(v["converged"] and v["conformation_kept"] for v in table.values()) and len(table) == 2
    assert {k: v["energy"] for k, v in table.items()} == pytest.approx({"reg_nitrile_tt_as_film_low": -1.51, "reg_nitrile_tt_as_film_low_registry_back": -0.05}, abs=0.006)
    assert table["reg_nitrile_tt_as_film_low"]["shear_along_chains_monomers"] == pytest.approx([1.12, -1.24], abs=0.006)
    assert np.abs(table["reg_nitrile_tt_as_film_low_registry_back"]["shear_along_chains_monomers"]).max() < 0.05
    frames = {str(a.info["frame"]): a for a in read(MLIP / "registry_as_film_nitrile_relaxed.extxyz", index=":")}
    up = labelled[CELLS["nitrile"][0]]
    for a in frames.values():
        assert all(any(tool.same_molecule(m, t) for t in tool.substituents(up)) for m in tool.substituents(a)) and tool.gauche(a, path) == 0
    # with its neighbours moved back it is level with the upright cell, its dipoles along the short edge and not the long one
    assert abs(abs(table["reg_nitrile_tt_as_film_low_registry_back"]["dipole_from_a_deg"]) - 180) < 30 and abs(abs(table["reg_nitrile_tt_as_film_low_registry_back"]["dipole_from_a_deg_labelled"]) - 90) < 15


def test_the_helix_cells_shaken_at_350_kelvin_come_back_to_where_they_started(labelled):
    """The fourth pod's eighth job: 1.5 ps at 350 K in the potential in each helix cell's own cell, then relaxed with the cell free."""
    from ase.io import read

    tool, path = _load("registry_cells"), _load("mlip_path")
    stored = json.loads((ROOT / "docs" / "predictions" / "registry_helix_annealed.json").read_text())
    table, shaken = stored["rows"], stored["shaken"]
    assert set(table) == set(shaken) == {f"reg_{c}_hx_annealed" for c in CELLS}
    for name, row in table.items():
        assert row["converged"] and row["conformation_kept"] and row["gauche_torsions"] == row["gauche_torsions_labelled"]
        assert abs(row["energy"]) < 0.015                                                    # kcal/mol per chain against the labelled helix cell: level
        assert row["shear_along_chains_monomers"] == pytest.approx(row["shear_along_chains_monomers_labelled"], abs=0.02)
        heat = shaken[name]
        assert (heat["md_steps"], heat["kelvin"]) == (1500, 350) and 380 < heat["hottest_kelvin"] < 420 and heat["converged"]
        assert abs(heat["meV_per_monomer_from_the_start"]) * 1e-3 * 14 * 23.0605 == pytest.approx(abs(row["energy"]), abs=0.006)
    frames = {str(a.info["frame"]): a for a in read(MLIP / "registry_helix_annealed_relaxed.extxyz", index=":")}
    for c, (_, helix) in CELLS.items():
        a, was = frames[f"reg_{c}_hx_annealed_in_all_s1"], labelled[helix]
        off = np.abs(a.cell.cellpar() - was.cell.cellpar())
        assert tool.gauche(a, path) == tool.gauche(was, path) and off[:3].max() < 0.03 and off[3:].max() < 0.3     # the same cell again, to 0.03 A and 0.3 degrees


def test_the_leader_does_not_have_the_low_cell_in_the_potential(labelled):
    """The leader made of the film composition's low cell, at rest in the potential (relaxed on this machine's processors, 519 steps)."""
    from ase.io import read

    tool, path = _load("registry_cells"), _load("mlip_path")
    table = json.loads((ROOT / "docs" / "predictions" / "registry_as_film_leader.json").read_text())["rows"]
    row = table["reg_leader_tt_as_film_low"]
    assert len(table) == 1 and row["converged"] and row["conformation_kept"] and row["gauche_torsions"] == 0
    assert row["energy"] == pytest.approx(1.02, abs=0.006) and row["shear_along_chains_monomers"] == pytest.approx([0.98, -0.96], abs=0.006)     # above its upright cell, in the low cell's registry
    a, up = read(MLIP / "registry_as_film_leader_relaxed.extxyz", index=0), labelled[CELLS["leader"][0]]
    assert str(a.info["frame"]) == "reg_leader_tt_as_film_low_in_all_s1" and all(any(tool.same_molecule(m, t) for t in tool.substituents(up)) for m in tool.substituents(a))
    # the other three compositions rest below their upright cells in the same cell
    others = {"film": json.loads((ROOT / "docs" / "predictions" / "kink_paths_film.json").read_text())["points"],
              "matrix": json.loads((ROOT / "docs" / "predictions" / "registry_as_film.json").read_text())["rows"]["reg_matrix_tt_as_film_low"]["energy"],
              "nitrile": json.loads((ROOT / "docs" / "predictions" / "registry_as_film_nitrile.json").read_text())["rows"]["reg_nitrile_tt_as_film_low"]["energy"]}
    film = others["film"]["all-trans, the sheet's registry, relaxed on"]["potential"] - others["film"]["all-trans, upright"]["potential"]
    assert (film, others["matrix"], others["nitrile"]) == pytest.approx((-1.27, -1.50, -1.51), abs=0.006)
    # one backbone bond a chain of the leader's cell is 44 and 46 degrees off trans there (the others' largest: 23 to 29)
    off = lambda x: np.abs(np.where(path.torsions_deg(x, path.backbone_torsions(x)) > 0, -180.0, 180.0) + path.torsions_deg(x, path.backbone_torsions(x))).max()
    assert off(a) == pytest.approx(46, abs=1.5) and off(read(MLIP / "registry_as_film_nitrile_relaxed.extxyz", index=0)) < 30 and off(read(MLIP / "paths" / "all_trans_tight.extxyz", index=0)) < 24


def test_the_helix_starts_one_monomer_along_are_the_last_eight_of_the_first_file():
    from ase.io import read

    first = {str(a.info["frame"]): a for a in read(MLIP / "registry_starts.extxyz", index=":")}
    own = read(MLIP / "registry_starts_helix_one.extxyz", index=":")
    assert [str(a.info["frame"]) for a in own] == list(first)[-8:] and all("_hx_" in str(a.info["frame"]) and str(a.info["frame"]).endswith("1_tz0") for a in own)
    for a in own:
        b = first[str(a.info["frame"])]
        assert np.allclose(a.positions, b.positions) and np.allclose(a.cell.array, b.cell.array) and a.get_chemical_symbols() == b.get_chemical_symbols()


def test_the_helix_cells_to_be_shaken_are_the_labelled_ones(labelled):
    from ase.io import read

    tool, path = _load("registry_cells"), _load("mlip_path")
    starts = read(MLIP / "registry_starts_helix_annealed.extxyz", index=":")
    assert [str(a.info["frame"]) for a in starts] == [f"reg_{c}_hx_annealed" for c in ("film", "matrix", "leader", "nitrile")]
    for a in starts:
        ref = labelled[CELLS[str(a.info["registry_label"])][1]]
        assert str(a.info["registry_of"]) == str(ref.info["frame"]) and tool._ints(a.info["registry_shear"]) == [0, 0] and int(a.info["registry_turn"]) == 0
        assert _same_places(ref.positions, a.positions, ref.cell.array) and a.cell.array == pytest.approx(ref.cell.array)
    line = next(l.split() for l in (ROOT / "examples" / "gpu_jobs" / "queue_20261122_registry.txt").read_text().splitlines() if l.startswith("anneal_card.sh "))
    assert (ROOT / line[1]).samefile(MLIP / "registry_starts_helix_annealed.extxyz") and line[2:] == ["4", "1500", "350", "700", "3000"]


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
    # a cell that is not at rest counts only when it is already well below its labelled cell; a cell in a DFT line already is not chosen again
    assert tool.chosen(rows, 0.3, {"film:tt": -1.0}, unsettled_below=-0.3)[0] == ("matrix_tt_-2.0", -2.0)
    assert tool.chosen(rows, 0.3, {"film:tt": -1.0}, unsettled_below=-2.5)[1] == ("matrix_tt_-0.9", -0.9)
    assert [n for n, _ in tool.chosen(rows, 0.3, {"film:tt": -1.0}, skip={"matrix_tt_-0.9", "film_tt_-1.1"})] == ["matrix_tt_-0.4", "nitrile_hx_-0.1", "matrix_hx_0.2"]
    assert tool.name_of("film", "tt", 1, -1, 117) == "reg_film_tt_p1m1_tp117" and tool.name_of("leader", "hx", -2, 2, 0) == "reg_leader_hx_m2p2_tz0"
    assert tool._ints("1 -1") == [1, -1] and tool._ints(np.array([2, -2])) == [2, -2]


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
