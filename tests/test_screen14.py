"""The screen of fourteen-monomer compositions made of the film composition's relaxed cells (examples/screen14.py;
docs/MLIP_SECOND_OPINION.md 3.25)."""
import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("ase")

ROOT = Path(__file__).resolve().parents[1]
MLIP = ROOT / "docs" / "predictions" / "mlip"
FILM = "vdf9-trfe4-cfe1-14"
KCAL = 23.0605
ROWS_A = {"vdf9-trfe4-cfe1-14": "VVRVVRVCVVRVVR", "vdf8-trfe4-cfe1-rvdf1-14": "IVRVVRVCVVRVVR", "vdf9-trfe3-cfe1-rvdf1-14": "IVRVVRVCVVVVVR",
          "vdf10-trfe3-cfe1-14": "VVRVVRVCVVVVVR", "vdf10-trfe4-14": "VVRVVRVVVVRVVR", "vdf8-trfe4-ctfe1-rvdf1-14": "IVRVVRVKVVRVVR",
          "vdf9-trfe3-ctfe1-rvdf1-14": "IVRVVRVKVVVVVR", "vdf9-trfe4-ctfe1-14": "VVRVVRVKVVRVVR"}
ROWS_B = {"vdf9-trfe4-rvdf1-first-14": "IVRVVRVVVVRVVR", "vdf10-trfe3-rvdf1-14": "IVRVVRVVVVVVVR", "vdf11-trfe3-14": "VVRVVRVVVVVVVR", "vdf9-trfe5-14": "VVRVVRVRVVRVVR",
          "vdf9-trfe4-cdfe1-14": "VVRVVRVDVVRVVR", "vdf8-trfe4-cfe1-cdfe1-14": "DVRVVRVCVVRVVR", "vdf9-trfe4-rvdf1-14": "VVRVVRVIVVRVVR", "vdf9-trfe4-e1-14": "VVRVVRVEVVRVVR",
          "vdf9-trfe3-an1-rvdf1-14": "IVRVVRVAVVVVVR"}


def _load(name):
    if name not in sys.modules:
        spec = importlib.util.spec_from_file_location(name, ROOT / "examples" / f"{name}.py")
        module = importlib.util.module_from_spec(spec)
        sys.modules[name] = module
        spec.loader.exec_module(module)
    return sys.modules[name]


@pytest.fixture(scope="module")
def sources():
    from ase.io import read

    frames = {}
    for f in ("cells14_tight_in_all_s1.extxyz", "paths/all_trans_tight.extxyz", "cells14_frames_pbe_d3bj.extxyz"):
        for a in read(MLIP / f, index=":"):
            frames[str(a.info["frame"])] = a
    return {"tt_low": frames["film_r_released_tight"], "tt_up": frames[f"row_{FILM}_TTTT_chain1_in_all_s1"], "hx_TGTGp": frames[f"row_{FILM}_TGTGp_chain1_in_all_s1"],
            "hx_GpTGT": frames[f"row_{FILM}_GpTGT_chain1_in_all_s1"], "hx_GTGpT": frames[f"row_{FILM}_GTGpT_chain1_in_all_s1"], "hx_TGpTG": frames[f"row_{FILM}_TGpTG_chain1_in_all_s1"],
            "matrix_tt": frames["row_vdf10-trfe4-14_TTTT_chain1_in_all_s1"], "leader_tt": frames["tt_leader_from_rows14an_r1"], "nitrile_tt": frames["tt_nitrile_from_rows14an_r1"]}


@pytest.fixture(scope="module")
def starts():
    from ase.io import read

    return read(MLIP / "screen14_starts_a.extxyz", index=":"), read(MLIP / "screen14_starts_b.extxyz", index=":")


# the rows polyfind has a cell of under chain1-cfe, and the leader, whose own upright cell is labelled
OWN = ("vdf8-trfe4-cfe1-rvdf1-14", "vdf8-trfe4-ctfe1-rvdf1-14", "vdf9-trfe4-ctfe1-14", "vdf9-trfe5-14", "vdf9-trfe4-cdfe1-14", "vdf8-trfe4-cfe1-cdfe1-14", "vdf9-trfe4-rvdf1-14",
       "vdf9-trfe4-e1-14", "vdf9-trfe3-an1-rvdf1-14")


def test_the_third_file_gives_every_row_all_trans_cells_in_more_patterns(sources):
    """tt_reg: made of the matrix's own upright cell (the regular pattern of twists); tt_own: the row's own upright cell where there is one."""
    from ase.geometry import find_mic
    from ase.io import read

    scr, reg, path = _load("screen14"), _load("registry_cells"), _load("mlip_path")
    cells = read(MLIP / "screen14_starts_c.extxyz", index=":")
    rows = dict(ROWS_A, **ROWS_B)
    order = [r for r in list(ROWS_A) if r not in (FILM, "vdf10-trfe4-14")] + [FILM] + list(ROWS_B)       # the matrix made of its own cell is its own cell: not here
    assert [str(a.info["frame"]) for a in cells] == [f"scr_{r}_{k}" for r in order for k in (("tt_own", "tt_reg") if r in OWN else ("tt_reg",))]
    for a in cells:
        row, kind = str(a.info["screen_row"]), str(a.info["screen_kind"])
        assert all(reg.same_molecule(chain, scr.cycle_of(rows[row])) for chain in reg.substituents(a)) and reg.gauche(a, path) == 0
        if kind == "tt_reg":
            source = sources["matrix_tt"]
            assert str(a.info["screen_from"]) == "row_vdf10-trfe4-14_TTTT_chain1_in_all_s1" and np.abs(a.cell.array - source.cell.array).max() < 1e-12
            backbone = sorted({q[0] for q in path.backbone_torsions(a)})
            assert np.abs(find_mic(a.positions[backbone] - source.positions[backbone], a.cell, pbc=True)[0]).max() < 1e-8
            assert int(a.info["screen_changed"]) == 2 * sum(x != y for x, y in zip(scr.cycle_of(rows[row]), scr.cycle_of("VVRVVRVVVVRVVR")))
        else:
            assert int(a.info["screen_changed"]) == 0 and str(a.info["screen_from"]) == ("tt_leader_from_rows14an_r1" if "an1" in row else f"row_{row}_TTTT_cfe")
    leader = {str(a.info["frame"]): a for a in cells}["scr_vdf9-trfe3-an1-rvdf1-14_tt_own"]
    assert np.abs(find_mic(leader.positions - sources["leader_tt"].positions, leader.cell, pbc=True)[0]).max() < 1e-8
    # the two patterns the made cells start in are not one: the film composition's upright cell and the matrix's own differ by 16 degrees rms in their backbone twists
    twists = lambda x: np.where(path.torsions_deg(x, path.backbone_torsions(x)) > 0, -180.0, 180.0) + path.torsions_deg(x, path.backbone_torsions(x))
    assert np.sqrt(((twists(sources["tt_up"]) - twists(sources["matrix_tt"])) ** 2).mean()) == pytest.approx(15.8, abs=0.3)


def test_a_sequences_letters_spell_the_chain_the_labelled_cells_have(sources):
    scr, reg = _load("screen14"), _load("registry_cells")
    for letters in list(ROWS_A.values()) + list(ROWS_B.values()) + ["XFEDKCARIV"]:
        assert scr.letters_of(scr.cycle_of(letters)) == letters
    # the labelled cells of the four compositions read as their sequences, hands and all (the builder's hand is + on every chiral carbon)
    for cell, letters in (("tt_up", "VVRVVRVCVVRVVR"), ("tt_low", "VVRVVRVCVVRVVR"), ("hx_GpTGT", "VVRVVRVCVVRVVR"), ("matrix_tt", "VVRVVRVVVVRVVR"),
                          ("leader_tt", "IVRVVRVAVVVVVR"), ("nitrile_tt", "VVRVVRVAVVRVVR")):
        for chain in reg.substituents(sources[cell]):
            assert reg.same_molecule(chain, scr.cycle_of(letters)), cell
    assert scr.letters_of(reg.substituents(sources["tt_up"])[0]) == "VVRVVRVCVVRVVR"
    assert scr.letters_of(scr.cycle_of("VVR")[1:] + scr.cycle_of("VVR")[:1]) is None       # read from the middle of a unit it is no sequence


def test_the_starts_are_the_film_compositions_cells_with_other_substituents(starts, sources):
    from ase.geometry import find_mic
    from ase.io import read

    scr, reg, path = _load("screen14"), _load("registry_cells"), _load("mlip_path")
    for cells, rows in zip(starts, (ROWS_A, ROWS_B)):
        assert [str(a.info["frame"]) for a in cells] == [f"scr_{row}_{kind}" for row in rows for kind in scr.KINDS]
        for a in cells:
            row, kind = str(a.info["screen_row"]), str(a.info["screen_kind"])
            source = sources[kind]
            assert str(a.info["screen_letters"]) == rows[row] and all(reg.same_molecule(chain, scr.cycle_of(rows[row])) for chain in reg.substituents(a))
            assert reg.gauche(a, path) == reg.gauche(source, path) == (0 if kind.startswith("tt") else 28)
            assert np.abs(a.cell.array - source.cell.array).max() < 1e-12 and int(a.info["n_monomers"]) == 28
            changed = sum(x != y for x, y in zip(scr.cycle_of(rows[row]), scr.cycle_of(ROWS_A[FILM])))
            assert int(a.info["screen_changed"]) == 2 * changed                                # two chains a cell
            backbone = sorted({q[0] for q in path.backbone_torsions(a)})
            moved, _ = find_mic(a.positions[backbone] - source.positions[backbone], a.cell, pbc=True)
            assert np.abs(moved).max() < 1e-8                                                   # no backbone carbon moves
            assert reg.closest_between_chains(a) is None or reg.closest_between_chains(a) > 2.0
    a, b = starts
    # the control: the film's composition made of itself is its six sources, wherever the cell's edge cuts them
    for cell in a[:6]:
        moved, _ = find_mic(cell.positions - sources[str(cell.info["screen_kind"])].positions, cell.cell, pbc=True)
        assert np.abs(moved).max() < 1e-8 and int(cell.info["screen_changed"]) == 0
    by = {str(x.info["frame"]): x for x in a + b}
    assert by["scr_vdf10-trfe4-14_tt_low"].get_chemical_formula() == "C56H48F64" and by["scr_vdf9-trfe3-ctfe1-rvdf1-14_tt_up"].get_chemical_formula() == "C56H46Cl2F64"
    assert by["scr_vdf9-trfe3-an1-rvdf1-14_hx_TGTGp"].get_chemical_formula() == "C58H52F58N2"
    # 3.24's cells again: the matrix in the low cell to a bond length's difference on the changed carbon, the leader in it to the atom
    matrix = read(MLIP / "registry_starts_as_film.extxyz", index=0)
    assert np.abs(find_mic(matrix.positions - by["scr_vdf10-trfe4-14_tt_low"].positions, matrix.cell, pbc=True)[0]).max() < 0.03
    leader = read(MLIP / "registry_starts_as_film_leader.extxyz", index=0)
    assert np.abs(find_mic(leader.positions - by["scr_vdf9-trfe3-an1-rvdf1-14_tt_low"].positions, leader.cell, pbc=True)[0]).max() < 1e-8


def _relaxed(folder, cells, energies, moving=()):
    """Cells written as the relaxer leaves them: NAME_in_all_s1.start.xyz with the potential's energy and a note."""
    from ase.io import write

    folder.mkdir(parents=True, exist_ok=True)
    for a in cells:
        b = a.copy()
        name = str(a.info["frame"])
        b.info = dict(a.info, frame=f"{name}_in_all_s1", potential_eV_per_monomer=energies[name],
                      note=f"{name}, relaxed in all_s1 with D3 (37 steps{', NOT converged' if name in moving else ''})")
        write(folder / f"{name}_in_all_s1.start.xyz", b, format="extxyz")


def test_the_table_gives_the_film_compositions_known_margins_from_its_own_cells(starts, sources, tmp_path):
    """The control row with the energies its six cells have in the potential: +0.072 from the low cell, +0.163 from the upright one, the low cell 1.27 per chain below."""
    scr = _load("screen14")
    a, _ = starts
    film = a[:6]
    energies = {str(c.info["frame"]): float(sources[str(c.info["screen_kind"])].info["potential_eV_per_monomer"]) for c in film}
    _relaxed(tmp_path / "relaxed", film, energies)
    subprocess.run([sys.executable, "-I", str(ROOT / "examples" / "screen14.py"), "table", "--starts", str(MLIP / "screen14_starts_a.extxyz"), "--relaxed", str(tmp_path / "relaxed"),
                    "--out", str(tmp_path / "table.json"), "--frames-out", str(tmp_path / "relaxed.extxyz")], check=True, capture_output=True, cwd=ROOT)
    rows = json.loads((tmp_path / "table.json").read_text())["rows"]
    assert list(rows) == list(ROWS_A) and all("margin" not in rows[r] and not rows[r]["cells"]["tt_low"]["relaxed"] for r in list(ROWS_A)[1:])
    row = rows[FILM]
    assert row["complete"] and row["at_rest"] == 6 and (row["lowest_tt"], row["lowest_hx"]) == ("tt_low", "hx_TGTGp")
    assert (row["margin"], row["margin_upright"], row["low_less_upright_per_chain"]) == pytest.approx((0.072, 0.163, -1.27), abs=0.006)
    # the four placements of the helix: two pairs, 0.46 per chain apart
    assert sorted(round(v, 2) for v in row["helix_placements_per_chain"].values()) == [0.0, 0.0, 0.46, 0.46]
    assert {k: round(v, 2) for k, v in row["all_trans_per_chain"].items()} == {"tt_low": 0.0, "tt_up": 1.27}
    # inside the window it is chosen with its lowest all-trans cell, its lowest helix cell and the next all-trans cell up; a control is left out by name
    assert scr.chosen(rows, -0.06, 0.14) == [(FILM, [f"scr_{FILM}_tt_low_in_all_s1", f"scr_{FILM}_hx_TGTGp_in_all_s1", f"scr_{FILM}_tt_up_in_all_s1"])]
    assert scr.chosen(rows, -0.06, 0.14, skip={FILM}) == [] and scr.chosen(rows, 0.10, 0.14) == []
    # a cell that leaves its conformation does not set a margin, and one not at rest makes the row incomplete
    energies[f"scr_{FILM}_hx_TGTGp"] += 0.01
    _relaxed(tmp_path / "relaxed", film, energies, moving={f"scr_{FILM}_hx_GTGpT"})
    subprocess.run([sys.executable, "-I", str(ROOT / "examples" / "screen14.py"), "table", "--starts", str(MLIP / "screen14_starts_a.extxyz"), "--relaxed", str(tmp_path / "relaxed"),
                    "--out", str(tmp_path / "table.json")], check=True, capture_output=True, cwd=ROOT)
    row = json.loads((tmp_path / "table.json").read_text())["rows"][FILM]
    assert row["lowest_hx"] == "hx_GTGpT" and not row["complete"] and row["at_rest"] == 5


def test_the_margins_taken_again_from_labelled_cells_are_dfts(starts, sources, tmp_path):
    """The control row with the DFT labels its cells have (3.22, 3.23): +0.063 from the low cell, which DFT has 0.84 per chain below the upright one."""
    from ase.io import read, write

    a, _ = starts
    film = a[:6]
    energies = {str(c.info["frame"]): float(sources[str(c.info["screen_kind"])].info["potential_eV_per_monomer"]) for c in film}
    _relaxed(tmp_path / "relaxed", film, energies)
    run = [sys.executable, "-I", str(ROOT / "examples" / "screen14.py")]
    subprocess.run(run + ["table", "--starts", str(MLIP / "screen14_starts_a.extxyz"), "--relaxed", str(tmp_path / "relaxed"), "--out", str(tmp_path / "table.json")],
                   check=True, capture_output=True, cwd=ROOT)
    low = {str(x.info["frame"]): x for x in read(MLIP / "cells14_third_pod_frames_pbe_d3bj.extxyz", index=":")}["film_r_released_tight"]
    labels = []
    for kind, cell in (("tt_low", low), ("tt_up", sources["tt_up"]), ("hx_TGTGp", sources["hx_TGTGp"])):
        b = cell.copy()
        b.info = dict(cell.info, frame=f"scr_{FILM}_{kind}_in_all_s1")
        labels.append(b)
    write(tmp_path / "labels.extxyz", labels, format="extxyz")
    subprocess.run(run + ["dft", "--table", str(tmp_path / "table.json"), "--labels", str(tmp_path / "labels.extxyz"), "--out", str(tmp_path / "dft.json")],
                   check=True, capture_output=True, cwd=ROOT)
    rows = json.loads((tmp_path / "dft.json").read_text())["rows"]
    assert list(rows) == [FILM]
    r = rows[FILM]
    assert (r["lowest_tt"], r["lowest_hx"]) == ("tt_low", "hx_TGTGp")
    assert (r["margin"], r["margin_potential_same_cells"], r["margin_potential"]) == pytest.approx((0.063, 0.072, 0.072), abs=0.0015)
    assert (r["low_less_upright_per_chain"], r["low_less_upright_per_chain_potential"]) == pytest.approx((-0.84, -1.27), abs=0.006)


def test_a_rows_kink_sheets_made_of_the_film_compositions_and_the_control_gives_its_known_sites(starts, sources, tmp_path):
    """The film composition's 28 relaxed sheets made of themselves, with the energies they have in the potential: each site as the card's table has it,
    1.27 per chain higher because it is measured from the low all-trans cell and not the upright one; the cheapest 3.99 deep on the chlorine."""
    from ase.geometry import find_mic
    from ase.io import read

    scr, reg, path = _load("screen14"), _load("registry_cells"), _load("mlip_path")
    run = [sys.executable, "-I", str(ROOT / "examples" / "screen14.py")]
    subprocess.run(run + ["sheets", "--frames", str(MLIP / "cells14_tight_in_all_s1.extxyz"), "--row", f"{FILM}=VVRVVRVCVVRVVR", "vdf9-trfe3-cfe1-rvdf1-14=IVRVVRVCVVVVVR",
                          "--out", str(tmp_path / "sheets.extxyz")], check=True, capture_output=True, cwd=ROOT)
    made = read(tmp_path / "sheets.extxyz", index=":")
    stored = {str(a.info["frame"]): a for a in read(MLIP / "cells14_tight_in_all_s1.extxyz", index=":")}
    sites = [f"{k}{m}" for k in range(0, 28, 2) for m in ("", "m")]
    assert [str(a.info["frame"]) for a in made] == [f"scr_{row}_sheet_at{site}" for row in (FILM, "vdf9-trfe3-cfe1-rvdf1-14") for site in sites]
    for a in made:
        source = stored[str(a.info["screen_from"])]
        assert str(a.info["screen_from"]) == f"sheet_{FILM}_belt4_at{a.info['screen_site']}_sheared_in_all_s1" and reg.gauche(a, path) == reg.gauche(source, path)
        assert all(reg.same_molecule(chain, scr.cycle_of(str(a.info["screen_letters"]))) for chain in reg.substituents(a))
    for a in made[:28]:
        assert np.abs(find_mic(a.positions - stored[str(a.info["screen_from"])].positions, a.cell, pbc=True)[0]).max() < 1e-8
    # the control's table: the film composition's sheets and its two all-trans cells as the relaxer would leave them
    energies = {str(a.info["frame"]): float(stored[str(a.info["screen_from"])].info["potential_eV_per_monomer"]) for a in made[:28]}
    _relaxed(tmp_path / "relaxed", made[:28], energies)
    flat = starts[0][:2]
    _relaxed(tmp_path / "relaxed", flat, {str(c.info["frame"]): float(sources[str(c.info["screen_kind"])].info["potential_eV_per_monomer"]) for c in flat})
    subprocess.run(run + ["sheets-table", "--starts", str(tmp_path / "sheets.extxyz"), "--relaxed", str(tmp_path / "relaxed"), "--out", str(tmp_path / "sheets.json")],
                   check=True, capture_output=True, cwd=ROOT)
    rows = json.loads((tmp_path / "sheets.json").read_text())["rows"]
    assert rows["vdf9-trfe3-cfe1-rvdf1-14"].get("no_all_trans") and "cheapest" not in rows["vdf9-trfe3-cfe1-rvdf1-14"]     # nothing of that row is relaxed here
    film = rows[FILM]
    assert (film["all_trans"], film["cheapest"], film["below_all_trans"]) == ("tt_low", "14m", 3)
    assert (film["sites"]["14m"]["energy"], film["sites"]["14m"]["units"], film["sites"]["12"]["units"]) == (pytest.approx(-3.99, abs=0.006), "CV", "VC")
    card = json.loads((ROOT / "docs" / "predictions" / "sheets14_potential_card.json").read_text())["rows"]
    for site, v in film["sites"].items():
        was = card[FILM + ("" if site == "0" else "@" + site)]
        assert v["energy"] - was["potential"] == pytest.approx(1.27, abs=0.006) and v["kink_kept"] == was["kink_kept"]


def test_the_first_files_margins_in_the_potential_and_its_two_controls(sources):
    """The fifth pod's card, file A, at rest at 0.004 eV/A: all-trans less helix per monomer from the lower of the two all-trans cells that file has."""
    from ase.io import read

    stored = json.loads((ROOT / "docs" / "predictions" / "screen14_potential_a.json").read_text())["rows"]
    assert list(stored) == list(ROWS_A) and all(r["complete"] and r["at_rest"] == 6 for r in stored.values())
    margins = {k: round(r["margin"], 3) for k, r in stored.items()}
    assert margins == {"vdf9-trfe4-cfe1-14": 0.072, "vdf8-trfe4-cfe1-rvdf1-14": -0.033, "vdf9-trfe3-cfe1-rvdf1-14": 0.088, "vdf10-trfe3-cfe1-14": 0.410,
                       "vdf10-trfe4-14": -0.087, "vdf8-trfe4-ctfe1-rvdf1-14": -0.491, "vdf9-trfe3-ctfe1-rvdf1-14": -0.246, "vdf9-trfe4-ctfe1-14": -0.305}
    # the controls: the film's composition made of itself comes back as it is (+0.072 and +0.163 are its known margins; its upright cell
    # relaxed sixty steps further here), and the matrix made of the film composition's cells has the matrix's known margin from its low cell (-0.089)
    film, matrix = stored[FILM], stored["vdf10-trfe4-14"]
    for kind, cell in film["cells"].items():
        assert cell["eV_per_monomer"] == pytest.approx(float(sources[kind].info["potential_eV_per_monomer"]), abs=2e-4) and cell["steps"] <= 60
    assert (film["lowest_tt"], film["lowest_hx"], matrix["lowest_tt"]) == ("tt_low", "hx_TGTGp", "tt_low")
    e = lambda a: float(a.info["potential_eV_per_monomer"])
    own = {str(a.info["frame"]): a for a in read(MLIP / "cells14_frames_pbe_d3bj.extxyz", index=":")}
    own_tt, own_hx = e(own["row_vdf10-trfe4-14_TTTT_chain1_in_all_s1"]), e(own["row_vdf10-trfe4-14_TGTGp_chain1_in_all_s1"])
    per_chain = 14 * KCAL
    made = {k: (c["eV_per_monomer"] - (own_tt if k.startswith("tt") else own_hx)) * per_chain for k, c in matrix["cells"].items()}
    # against the matrix's own labelled cells: its low cell as 3.24 has it, its helix to 0.03, and the cell in the film composition's upright pattern 1.01 above
    assert (made["tt_low"], made["tt_up"]) == pytest.approx((-1.50, 1.01), abs=0.006)
    assert sorted(round(v, 2) for k, v in made.items() if k.startswith("hx")) == [-0.03, -0.03, 0.08, 0.08]      # the matrix's own four: 0.00, 0.00, 0.14, 0.14
    # the row in the window: its lowest all-trans cell is the one in the upright pattern, the low cell 0.18 above it
    row = stored["vdf9-trfe3-cfe1-rvdf1-14"]
    assert (row["lowest_tt"], row["lowest_hx"]) == ("tt_up", "hx_GTGpT") and row["low_less_upright_per_chain"] == pytest.approx(0.18, abs=0.006)
    # the CTFE row's low cell is 6.35 below its upright one in the potential, in another registry (2, -1) and with a bond 43 degrees off trans: not believed without DFT
    ctfe = stored["vdf9-trfe4-ctfe1-14"]
    assert ctfe["low_less_upright_per_chain"] == pytest.approx(-6.35, abs=0.006) and [round(x) for x in ctfe["cells"]["tt_low"]["shear_along_chains_monomers"]] == [2, -1]
    frames = read(MLIP / "screen14_relaxed_a.extxyz", index=":")
    assert len(frames) == 48 and {str(a.info["frame"]) for a in frames} == {c["frame"] for r in stored.values() for c in r["cells"].values()}


def test_the_fourth_file_gives_every_row_cells_in_the_leaders_packings():
    """hx_ldr and tt_ldr: made of the leader's own lowest helix cell and of its own upright all-trans cell; hx_an (six rows): of the nitrile composition's
    lowest helix cell.  A nitrile made something else loses its nitrogen.  The rows in the order they matter: the pod's job stops at its budget."""
    from ase.geometry import find_mic
    from ase.io import read

    scr, reg, path = _load("screen14"), _load("registry_cells"), _load("mlip_path")
    cells = read(MLIP / "screen14_starts_d.extxyz", index=":")
    rows = dict(ROWS_A, **ROWS_B)
    leader = "vdf9-trfe3-an1-rvdf1-14"
    three = ["vdf9-trfe3-cfe1-rvdf1-14", "vdf8-trfe4-cfe1-rvdf1-14", FILM, "vdf10-trfe4-14", "vdf10-trfe3-rvdf1-14"]
    names = [str(a.info["frame"]) for a in cells]
    assert names[:8] == [f"scr_{r}_{k}" for r in three[:2] for k in ("hx_ldr", "tt_ldr", "hx_an")] + [f"scr_{leader}_hx_ldr", f"scr_{leader}_hx_an"]
    assert names[8:17] == [f"scr_{r}_{k}" for r in three[2:] for k in ("hx_ldr", "tt_ldr", "hx_an")] and len(cells) == 39
    rest = [r for r in rows if r not in three and r != leader]
    assert sorted(names[17:]) == sorted(f"scr_{r}_{k}" for r in rest for k in ("hx_ldr", "tt_ldr"))
    assert names[17:23] == [f"scr_{r}_{k}" for r in ("vdf9-trfe4-rvdf1-first-14", "vdf9-trfe3-ctfe1-rvdf1-14", "vdf8-trfe4-ctfe1-rvdf1-14") for k in ("hx_ldr", "tt_ldr")]
    own = {str(a.info["frame"]): a for f in ("cells14_tight_in_all_s1.extxyz", "cells14_frames_pbe_d3bj.extxyz") for a in read(MLIP / f, index=":")}
    froms = {"hx_ldr": "helix_leader_TGTGp_r1_in_all_s1", "tt_ldr": "tt_leader_from_rows14an_r1", "hx_an": "row_vdf9-trfe4-an1-14_GpTGT_cfe_in_all_s1"}
    for a in cells:
        row, kind = str(a.info["screen_row"]), str(a.info["screen_kind"])
        source = own[froms[kind]]
        assert str(a.info["screen_from"]) == froms[kind] and np.abs(a.cell.array - source.cell.array).max() < 1e-12
        assert all(reg.same_molecule(chain, scr.cycle_of(rows[row])) for chain in reg.substituents(a))
        assert reg.gauche(a, path) == reg.gauche(source, path) and (reg.gauche(a, path) == 0) == kind.startswith("tt")     # the cell it was made of, bond for bond
        assert int((a.numbers == 7).sum()) == (2 if "A" in rows[row] else 0) and len(a) == (170 if "A" in rows[row] else 168)
        backbone, was = sorted({q[0] for q in path.backbone_torsions(a)}), sorted({q[0] for q in path.backbone_torsions(source)})
        assert np.abs(find_mic(a.positions[backbone] - source.positions[was], a.cell, pbc=True)[0]).max() < 1e-8
    same = cells[names.index(f"scr_{leader}_hx_ldr")]
    assert int(same.info["screen_changed"]) == 0 and np.abs(find_mic(same.positions - own[froms["hx_ldr"]].positions, same.cell, pbc=True)[0]).max() < 1e-8


def test_the_second_files_margins_and_what_the_leader_as_a_control_says_of_them():
    """File B at rest (the table of files A and B): margins from the film composition's two all-trans patterns and four helix cells.  The leader, whose own
    cells the potential has, comes back 0.27 a monomer too low: 0.61 per chain on its all-trans side and 4.42 on its helix side."""
    from ase.io import read

    stored = json.loads((ROOT / "docs" / "predictions" / "screen14_potential_b.json").read_text())["rows"]
    assert list(stored) == list(ROWS_A) + list(ROWS_B) and all(r["at_rest"] == 6 for r in stored.values())
    margins = {k: round(r["margin"], 3) for k, r in stored.items() if k in ROWS_B}
    assert margins == {"vdf9-trfe4-rvdf1-first-14": -0.178, "vdf10-trfe3-rvdf1-14": -0.091, "vdf11-trfe3-14": 0.243, "vdf9-trfe5-14": -0.102, "vdf9-trfe4-cdfe1-14": -0.110,
                       "vdf8-trfe4-cfe1-cdfe1-14": 0.313, "vdf9-trfe4-rvdf1-14": -0.491, "vdf9-trfe4-e1-14": -0.249, "vdf9-trfe3-an1-rvdf1-14": -0.342}
    # inside the window the pod asks DFT about (-0.06 to +0.14): two rows of file A and none of file B
    assert [k for k, r in stored.items() if -0.06 <= r["margin"] <= 0.14 and k not in (FILM, "vdf10-trfe4-14")] == ["vdf8-trfe4-cfe1-rvdf1-14", "vdf9-trfe3-cfe1-rvdf1-14"]
    own = {str(a.info["frame"]): a for a in read(MLIP / "cells14_frames_pbe_d3bj.extxyz", index=":")}
    e = lambda a: float(a.info["potential_eV_per_monomer"])
    tt, hx = e(own["tt_leader_from_rows14an_r1"]), e(own["helix_leader_TGTGp_r1"])
    assert (tt - hx) * KCAL == pytest.approx(-0.070, abs=0.002)
    leader, per_chain = stored["vdf9-trfe3-an1-rvdf1-14"], 14 * KCAL
    made = {k: (c["eV_per_monomer"] - (tt if k.startswith("tt") else hx)) * per_chain for k, c in leader["cells"].items() if c.get("relaxed")}       # file C's were not yet
    assert (made["tt_up"], made["tt_low"]) == pytest.approx((0.61, 1.02), abs=0.006)
    assert [round(made[k], 2) for k in ("hx_GpTGT", "hx_GTGpT", "hx_TGTGp", "hx_TGpTG")] == [4.42, 4.87, 6.48, 9.32]
    assert leader["margin"] == pytest.approx((tt - hx) * KCAL + (made["tt_up"] - made["hx_GpTGT"]) / 14, abs=1e-6)
    frames = read(MLIP / "screen14_relaxed_b.extxyz", index=":")
    assert len(frames) == 54 and {str(a.info["frame"]) for a in frames} == {c["frame"] for k, r in stored.items() if k in ROWS_B for c in r["cells"].values() if c.get("relaxed")}


def test_the_film_compositions_chains_turned_without_their_registry_are_level_with_upright_in_dft():
    """The fifth pod's first label: the film composition's all-trans cell with every chain turned 117 degrees and its neighbours where the upright cell has
    them.  DFT has it level with the upright cell (the potential: -0.02), so the 0.84 of the low cell is not the turn's."""
    from ase.io import read

    cell = {str(a.info["frame"]): a for a in read(MLIP / "screen14_frames_pbe_d3bj.extxyz", index=":")}["film_all_trans_turned_117_tight"]
    up = {str(a.info["frame"]): a for a in read(MLIP / "cells14_frames_pbe_d3bj.extxyz", index=":")}[f"row_{FILM}_TTTT_chain1_in_all_s1"]
    dft = lambda a: float(a.info["REF_energy"]) / int(a.info["n_monomers"])
    pot = lambda a: float(a.info["potential_eV_per_monomer"])
    assert (dft(cell) - dft(up)) * 14 * KCAL == pytest.approx(0.00, abs=0.006) and (pot(cell) - pot(up)) * 14 * KCAL == pytest.approx(-0.02, abs=0.006)
    assert np.abs(cell.arrays["REF_forces"]).max() < 0.1


def test_the_second_sheets_file_is_two_rows_made_of_the_leaders_own_sheets():
    """IVRVVRVCVVVVVR and IVRVVRVVVVVVVR at the leader's cheapest kink sites (and, the first, the two beside the leader's nitrile), each made of the
    leader's relaxed sheet there: the nitrile a chlorine unit's or a VDF's substituents, its nitrogen gone, the kink bond for bond."""
    from ase.geometry import find_mic
    from ase.io import read

    scr, reg, path = _load("screen14"), _load("registry_cells"), _load("mlip_path")
    cells = read(MLIP / "screen14_sheets_2.extxyz", index=":")
    rows = {"vdf9-trfe3-cfe1-rvdf1-14": ("IVRVVRVCVVVVVR", ["18m", "18", "16m", "20m", "12", "14", "16", "20", "4"]),
            "vdf10-trfe3-rvdf1-14": ("IVRVVRVVVVVVVR", ["18m", "18", "16m", "20m", "16", "20", "4"])}
    names = [str(a.info["frame"]) for a in cells]
    assert len(cells) == 16 and all([n for n in names if r in n] == [f"scr_{r}_sheetldr_at{s}" for s in sites] for r, (_, sites) in rows.items())
    assert names[:4] == ["scr_vdf9-trfe3-cfe1-rvdf1-14_sheetldr_at18m", "scr_vdf10-trfe3-rvdf1-14_sheetldr_at18m", "scr_vdf9-trfe3-cfe1-rvdf1-14_sheetldr_at18",
                         "scr_vdf10-trfe3-rvdf1-14_sheetldr_at18"]                                    # the two rows in turn
    tight = {str(a.info["frame"]): a for a in read(MLIP / "cells14_tight_in_all_s1.extxyz", index=":")}
    for a in cells:
        row, site = str(a.info["screen_row"]), str(a.info["screen_site"])
        source = tight[f"sheet_vdf9-trfe3-an1-rvdf1-14_belt4_at{site}_sheared_in_all_s1"]
        assert str(a.info["screen_from"]) == str(source.info["frame"]) and int(a.info["screen_changed"]) == 2 and len(a) == 168 and int((a.numbers == 7).sum()) == 0
        assert all(reg.same_molecule(chain, scr.cycle_of(rows[row][0])) for chain in reg.substituents(a)) and reg.gauche(a, path) == reg.gauche(source, path) > 0
        backbone, was = sorted({q[0] for q in path.backbone_torsions(a)}), sorted({q[0] for q in path.backbone_torsions(source)})
        assert np.abs(find_mic(a.positions[backbone] - source.positions[was], a.cell, pbc=True)[0]).max() < 1e-8


def test_the_fourth_files_table_the_leader_comes_back_and_three_rows_move():
    """Every row given the leader's own packings too (file D, 31 of its 39 cells relaxed before the job's 3000 s).  The leader as a control now has its
    known margin; the rows with its reversed unit first and its VDF eleventh take its helix packing, the others do not."""
    from ase.io import read

    stored = json.loads((ROOT / "docs" / "predictions" / "screen14_potential_d.json").read_text())["rows"]
    margins = {k: round(r["margin"], 3) for k, r in stored.items()}
    assert margins == {"vdf9-trfe4-cfe1-14": 0.072, "vdf8-trfe4-cfe1-rvdf1-14": -0.033, "vdf9-trfe3-cfe1-rvdf1-14": 0.201, "vdf10-trfe3-cfe1-14": 0.410, "vdf10-trfe4-14": -0.087,
                       "vdf8-trfe4-ctfe1-rvdf1-14": -0.459, "vdf9-trfe3-ctfe1-rvdf1-14": -0.053, "vdf9-trfe4-ctfe1-14": -0.305, "vdf9-trfe4-rvdf1-first-14": -0.178,
                       "vdf10-trfe3-rvdf1-14": 0.071, "vdf11-trfe3-14": 0.243, "vdf9-trfe5-14": -0.102, "vdf9-trfe4-cdfe1-14": -0.110, "vdf8-trfe4-cfe1-cdfe1-14": 0.313,
                       "vdf9-trfe4-rvdf1-14": -0.491, "vdf9-trfe4-e1-14": -0.249, "vdf9-trfe3-an1-rvdf1-14": -0.069}
    before = {k: round(r["margin"], 3) for k, r in json.loads((ROOT / "docs" / "predictions" / "screen14_potential_c.json").read_text())["rows"].items()}
    assert {k: (before[k], v) for k, v in margins.items() if abs(v - before[k]) > 0.02} == {
        "vdf9-trfe3-cfe1-rvdf1-14": (0.088, 0.201), "vdf8-trfe4-ctfe1-rvdf1-14": (-0.513, -0.459), "vdf9-trfe3-ctfe1-rvdf1-14": (-0.246, -0.053),
        "vdf10-trfe3-rvdf1-14": (-0.091, 0.071), "vdf9-trfe3-an1-rvdf1-14": (-0.386, -0.069)}
    # the leader: its own cells in this potential give -0.070
    own = {str(a.info["frame"]): a for a in read(MLIP / "cells14_frames_pbe_d3bj.extxyz", index=":")}
    e = lambda a: float(a.info["potential_eV_per_monomer"])
    assert (e(own["tt_leader_from_rows14an_r1"]) - e(own["helix_leader_TGTGp_r1"])) * KCAL == pytest.approx(-0.070, abs=0.002)
    leader = stored["vdf9-trfe3-an1-rvdf1-14"]
    assert (leader["lowest_tt"], leader["lowest_hx"]) == ("tt_own", "hx_ldr")
    # the leader's helix packing against each row's lowest helix cell of the film composition's four, kcal/mol per chain
    def packing(row):
        hx = stored[row]["helix_placements_per_chain"]
        return round(hx["hx_ldr"] - min(v for k, v in hx.items() if k not in ("hx_ldr", "hx_an")), 2)
    assert {r: packing(r) for r in ("vdf9-trfe3-an1-rvdf1-14", "vdf9-trfe3-cfe1-rvdf1-14", "vdf10-trfe3-rvdf1-14", "vdf9-trfe3-ctfe1-rvdf1-14", "vdf8-trfe4-ctfe1-rvdf1-14")} == {
        "vdf9-trfe3-an1-rvdf1-14": -4.43, "vdf9-trfe3-cfe1-rvdf1-14": -1.62, "vdf10-trfe3-rvdf1-14": -2.26, "vdf9-trfe3-ctfe1-rvdf1-14": -2.70, "vdf8-trfe4-ctfe1-rvdf1-14": -0.75}
    assert {r: packing(r) for r in (FILM, "vdf10-trfe4-14", "vdf8-trfe4-cfe1-rvdf1-14", "vdf9-trfe4-rvdf1-first-14")} == {
        FILM: 1.84, "vdf10-trfe4-14": 1.66, "vdf8-trfe4-cfe1-rvdf1-14": 1.85, "vdf9-trfe4-rvdf1-first-14": 1.90}
    # the rows file D did not reach with the leader's helix packing, and the cells it left moving
    waiting = {k: sorted(c for c, v in r["cells"].items() if not v.get("relaxed")) for k, r in stored.items()}
    assert {k: v for k, v in waiting.items() if "hx_ldr" in v} == {"vdf9-trfe4-ctfe1-14": ["hx_ldr"], "vdf8-trfe4-cfe1-cdfe1-14": ["hx_ldr", "tt_ldr"], "vdf9-trfe4-rvdf1-14": ["hx_ldr", "tt_ldr"]}
    frames = read(MLIP / "screen14_relaxed_d.extxyz", index=":")
    assert len(frames) == 31 and sum(len(v) for v in waiting.values()) == 8


def test_where_the_kink_sits_in_the_two_rows_with_the_leaders_other_features():
    """The kink sheets relaxed on the card, kcal/mol per chain against each row's lowest all-trans cell: with a chlorine unit the kink is cheapest on the
    chlorine, as in the film's composition; without a pin it is cheapest on the run of VDF units, as in the leader."""
    one = json.loads((ROOT / "docs" / "predictions" / "screen14_sheets_1.json").read_text())["rows"]["vdf9-trfe3-cfe1-rvdf1-14"]
    two = json.loads((ROOT / "docs" / "predictions" / "screen14_sheets_2.json").read_text())["rows"]
    kept = lambda row: {k: round(v["energy"], 2) for k, v in row["sites"].items() if v.get("relaxed") and v["kink_kept"] and v["converged"]}
    assert one["all_trans"] == "tt_ldr" and one["cheapest"] == "12" and sum(bool(v.get("relaxed")) for v in one["sites"].values()) == 25
    k1 = kept(one)
    assert (k1["12"], k1["14m"], k1["16"], k1["18"]) == (-4.84, -4.21, -1.84, -1.55) and one["sites"]["12"]["units"] == "VC" and one["sites"]["14m"]["units"] == "CV"
    cfe, nopin = two["vdf9-trfe3-cfe1-rvdf1-14"], two["vdf10-trfe3-rvdf1-14"]
    k2, k3 = kept(cfe), kept(nopin)
    assert (k2["12"], k2["16m"], k2["18m"], k2["20m"]) == (-2.92, -2.21, -2.04, -1.78)                  # the same row in the leader's packing: the chlorine site still the cheapest
    assert nopin["all_trans"] == "tt_low" and nopin["cheapest"] == "18m" and (k3["18m"], k3["20m"], k3["16m"]) == (-1.87, -1.72, -1.60)
    assert all(nopin["sites"][s]["units"] == "VV" for s in ("18m", "20m", "16m"))


def test_the_candidates_margin_in_dft_is_the_potentials_on_its_lowest_cells_and_not_on_the_others():
    """IVRVVRVCVVVVVR, four cells labelled on the fifth pod (PBE-D3(BJ) single points on the potential's cells).  From its lowest labelled cells DFT has
    +0.196 and the potential +0.204; on its upright all-trans cell and its helix in the film composition's packing DFT has -0.035 and the potential +0.088."""
    from ase.io import read

    row = json.loads((ROOT / "docs" / "predictions" / "screen14_dft.json").read_text())["rows"]["vdf9-trfe3-cfe1-rvdf1-14"]
    assert (row["lowest_tt"], row["lowest_hx"]) == ("tt_up", "hx_ldr") and sorted(row["cells"]) == ["hx_GTGpT", "hx_ldr", "tt_low", "tt_up"]
    assert (row["margin"], row["margin_potential_same_cells"], row["margin_potential"]) == pytest.approx((0.196, 0.204, 0.201), abs=0.0006)
    dft = {k: c["dft_eV_per_monomer"] for k, c in row["cells"].items()}
    pot = {k: c["potential_eV_per_monomer"] for k, c in row["cells"].items()}
    per_chain = lambda e, k: (e[k] - e["tt_up"]) * 14 * KCAL
    assert [round(per_chain(dft, k), 2) for k in ("tt_low", "hx_GTGpT", "hx_ldr")] == [1.51, 0.49, -2.75]
    assert [round(per_chain(pot, k), 2) for k in ("tt_low", "hx_GTGpT", "hx_ldr")] == [0.18, -1.23, -2.85]
    assert (dft["tt_up"] - dft["hx_GTGpT"]) * KCAL == pytest.approx(-0.035, abs=0.0006) and (pot["tt_up"] - pot["hx_GTGpT"]) * KCAL == pytest.approx(0.088, abs=0.0006)
    assert max(c["dft_largest_force_eV_per_A"] for c in row["cells"].values()) < 0.23
    frames = read(MLIP / "screen14_frames_pbe_d3bj.extxyz", index=":")
    assert [str(a.info["dft_run"]) for a in frames] == ["processors"] * 4 + ["card"] and str(frames[-1].info["frame"]) == "scr_vdf9-trfe3-cfe1-rvdf1-14_hx_ldr_in_all_s1"


def test_what_a_cells_chains_are_from_their_gauche_bonds():
    scr = _load("screen14")
    assert scr._state([[], []], 28) == "all-trans" and scr._state([[(18, -76.0)], [(18, -76.0)]], 28) == "one gauche bond"
    assert scr._state([[(16, -55.0), (18, 60.0)]] * 2, 28) == "kink" and scr._state([[(16, -55.0), (20, 57.0)]] * 2, 28) == "stretched kink"
    assert scr._state([[(0, 60.0), (26, -60.0)]] * 2, 28) == "kink" and scr._state([[(0, 60.0), (24, -60.0)]] * 2, 28) == "stretched kink"      # round the end of the chain
    assert scr._state([[(16, 55.0), (20, 57.0)]] * 2, 28) == "other" and scr._state([[], [(18, 60.0)]], 28) == "mixed"


def test_neither_rows_cheapest_kink_moves_as_the_leaders_does():
    """Step scans on the fifth pod's card, kcal/mol per chain against each relaxed sheet.  A hop is a first leg that arrives at a stretched kink.  The stored
    scans of the film's composition and the leader (3.23) go through the same summary: 7.46 from its scan's start (7.81 from its relaxed sheet) and 5.67."""
    paths = json.loads((ROOT / "docs" / "predictions" / "screen14_hops.json").read_text())["paths"]
    hop = {k: round(v["hop"], 2) for k, v in paths.items()}
    assert hop == {"nopin": 7.55, "nopin18": 5.58, "lcfe": 6.45, "lcfev": 4.87, "film": 7.46, "leader": 5.67}
    arrives = {k: {o: (round(v["first_leg"]["top"], 2), v["first_leg"]["arrives_as"], round(v["first_leg"]["arrived"], 2)) for o, v in p["orders"].items()} for k, p in paths.items()}
    assert arrives["nopin"] == {"a": (7.55, "stretched kink", 0.68), "b": (7.62, "one gauche bond", 6.22)}
    assert arrives["nopin18"] == {"a": (3.41, "all-trans", 1.09), "b": (5.58, "stretched kink", 0.79)}          # at the leader's site and hand: the leader's hop, and a kink that is barely held
    assert arrives["lcfe"] == {"a": (6.45, "stretched kink", 0.26), "b": (7.28, "stretched kink", 1.66)}
    assert arrives["lcfev"] == {"a": (4.60, "one gauche bond", 2.87), "b": (4.87, "stretched kink", -0.34)}
    assert arrives["leader"] == {"a": (7.01, "stretched kink", 0.87), "b": (5.67, "stretched kink", 0.85)}
    assert {k: round(v["removed_over"], 2) for k, v in paths.items()} == {"nopin": 7.69, "nopin18": 3.41, "lcfe": 10.39, "lcfev": 4.60, "film": 8.42, "leader": 7.25}
    assert paths["nopin"]["sheet"] == "scr_vdf10-trfe3-rvdf1-14_sheetldr_at18m_in_all_s1" and paths["lcfe"]["sheet"] == "scr_vdf9-trfe3-cfe1-rvdf1-14_sheet_at12_in_all_s1"
    assert paths["nopin18"]["orders"]["a"]["first_leg"]["not_converged"] == ["i12"] and all(not o["first_leg"]["not_converged"] for k in ("nopin", "lcfe", "lcfev") for o in paths[k]["orders"].values())
    # where each scanned sheet stands against the row's cheapest: the no-pin row's other hand 2.23 above, the CFE row's VDF run 2.63 above its chlorine
    two = json.loads((ROOT / "docs" / "predictions" / "screen14_sheets_2.json").read_text())["rows"]
    one = json.loads((ROOT / "docs" / "predictions" / "screen14_sheets_1.json").read_text())["rows"]["vdf9-trfe3-cfe1-rvdf1-14"]
    nopin = two["vdf10-trfe3-rvdf1-14"]["sites"]
    assert nopin["18"]["energy"] - nopin["18m"]["energy"] == pytest.approx(2.23, abs=0.006)
    assert two["vdf9-trfe3-cfe1-rvdf1-14"]["sites"]["16m"]["energy"] - one["sites"]["12"]["energy"] == pytest.approx(2.63, abs=0.006)


def test_the_two_rows_on_the_four_axes():
    """Each row's margin and hop on the inputs of the two rows nearest it (examples/screen14_axes.sh): the CFE row is ten times the film composition's
    rate at a seventieth of its loss and seven to nine tenths of its strain, and its strain moves eightfold over 0.05 of margin either way; the no-pin row is
    the film's composition again."""
    out = {}
    for which in "ab":
        d = json.loads((ROOT / "docs" / "predictions" / f"four_axes_screen14_{which}.json").read_text())
        assert d["barriers_kcal_per_chain"]["vdf9-trfe4-cfe1-14"] == 7.81
        for shift, table in d["shifts"].items():
            for name, r in table["rows"].items():
                b = r["by_the_barrier"]
                through = b["blocked_pn_nm2_through_the_crystal"]
                out[(d["stands_in"].get(name, "film"), shift)] = (round(r["margin"], 3), round(b["f_relax_hz"] / 1e3), round(b["tan_delta_34kHz"], 4),
                                                                  round(b["strain_ac_at_the_heat_budget_34kHz"] * 100, 3), None if through is None else round(through, 1))
    assert out[("film", "+0.00")] == (0.063, 364, 0.0877, 0.081, 7.6)
    assert out[("IVRVVRVCVVVVVR on IVRVVRVCVVRVVR's inputs", "+0.00")] == (0.196, 3569, 0.0013, 0.073, None)
    assert out[("IVRVVRVCVVVVVR on the leader's inputs", "+0.00")] == (0.196, 3569, 0.0012, 0.057, 2.7)
    assert out[("IVRVVRVVVVVVVR on the leader's inputs", "+0.00")] == (0.071, 563, 0.0542, 0.087, 4.4)
    assert out[("IVRVVRVVVVVVVR on the matrix's inputs", "+0.00")] == (0.071, 563, 0.0571, 0.094, 9.8)
    strain = lambda who, shift: out[(who, shift)][3]
    assert (strain("IVRVVRVCVVVVVR on the leader's inputs", "-0.05"), strain("IVRVVRVCVVVVVR on the leader's inputs", "+0.05")) == (0.131, 0.014)
    assert (strain("IVRVVRVCVVVVVR on IVRVVRVCVVRVVR's inputs", "-0.05"), strain("IVRVVRVCVVVVVR on IVRVVRVCVVRVVR's inputs", "+0.05")) == (0.143, 0.018)
    assert (strain("IVRVVRVVVVVVVR on the leader's inputs", "-0.05"), strain("IVRVVRVVVVVVVR on the leader's inputs", "+0.05")) == (0.083, 0.075)


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
