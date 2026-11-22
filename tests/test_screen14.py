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


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
