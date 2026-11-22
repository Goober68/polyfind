"""A path between two labelled cells in a potential (examples/mlip_path.py; docs/MLIP_SECOND_OPINION.md 3.23): the parts that need no potential."""

import importlib.util
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
MLIP = ROOT / "docs" / "predictions" / "mlip"
CFE, LEADER = "vdf9-trfe4-cfe1-14", "vdf9-trfe3-an1-rvdf1-14"
FILM = (f"row_{CFE}_TTTT_chain1_in_all_s1_in_all_s1", f"sheet_{CFE}_belt4_at14m_sheared_in_all_s1_in_all_s1")
LEAD = ("tt_leader_from_rows14an_r1_in_all_s1", f"sheet_{LEADER}_belt4_at18_sheared_in_all_s1")


def _tool():
    spec = importlib.util.spec_from_file_location("mlip_path", ROOT / "examples" / "mlip_path.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def frames():
    from ase.io import read

    return {str(a.info["frame"]): a for a in read(MLIP / "cells14_frames_pbe_d3bj.extxyz", index=":")}


def test_the_backbone_is_two_chains_of_28_carbons_and_a_kink_pair_is_two_torsions_a_chain(frames):
    tool = _tool()
    for (flat, sheet), signs in ((FILM, (+1, -1)), (LEAD, (-1, +1))):
        quads = tool.backbone_torsions(frames[flat])
        assert len(quads) == 56 and tool.backbone_torsions(frames[sheet]) == quads   # the same atoms in the same order, the nitrile's carbon left off
        assert len({i for q in quads for i in q}) == 56
        a0, a1 = tool.torsions_deg(frames[flat], quads), tool.torsions_deg(frames[sheet], quads)
        assert np.abs(a0).min() > 155                                               # all-trans
        turn = tool.arc(a0, a1)
        led = [q for q in range(56) if abs(turn[q]) > 60]
        assert len(led) == 4 and [led[1] - led[0], led[3] - led[2]] == [2, 2]       # G T G' on each chain, the two chains at the same site
        assert led[2] - led[0] == 28
        assert np.abs(a1[led]) == pytest.approx([62, 62, 62, 62], abs=13)           # gauche
        assert [int(np.sign(turn[q])) for q in led] == [signs[0], signs[1], signs[0], signs[1]]
        assert 95 < min(abs(turn[q]) for q in led) and max(abs(turn[q]) for q in led) < 140


def test_the_short_way_between_two_angles():
    tool = _tool()
    assert tool.arc(-165.0, 60.8) == pytest.approx(-134.2) and tool.arc(170.0, -170.0) == pytest.approx(20.0)
    assert tool.arc(np.array([10.0, -170.0]), np.array([20.0, 170.0])) == pytest.approx([10.0, -20.0])


def test_the_restraints_forces_and_stress_are_its_energys_derivatives(frames):
    tool = _tool()
    sheet = frames[FILM[1]]
    quads = tool.backbone_torsions(sheet)
    a1 = tool.torsions_deg(sheet, quads)
    led = [14, 16, 42, 44]
    calc = tool.restraint_calculator([quads[q] for q in led], [a1[q] - 25.0 + 7.0 * i for i, q in enumerate(led)], 10.0)
    a = sheet.copy()
    a.calc = calc
    energy, forces, stress = a.get_potential_energy(), a.get_forces(), a.get_stress(voigt=False)
    assert energy > 0.5 and np.abs(forces.sum(axis=0)).max() < 1e-12
    assert set(np.flatnonzero(np.abs(forces).max(axis=1) > 0)) <= {i for q in led for i in quads[q]}
    h = 1e-5
    for atom, x in ((quads[14][0], 0), (quads[14][1], 2), (quads[44][3], 1)):
        b = a.copy()
        b.calc = calc
        b.positions[atom, x] += h
        up = b.get_potential_energy()
        b.positions[atom, x] -= 2 * h
        assert -(up - b.get_potential_energy()) / (2 * h) == pytest.approx(forces[atom, x], abs=1e-7)
    for i, j in ((0, 0), (1, 2), (0, 2)):   # a torsion changes when the cell is sheared: the restraint has a stress
        eps = np.zeros((3, 3))
        eps[i, j] += h / 2
        eps[j, i] += h / 2
        es = []
        for sign in (1, -1):
            b = a.copy()
            b.calc = calc
            b.set_cell(a.cell.array @ (np.eye(3) + sign * eps), scale_atoms=True)
            es.append(b.get_potential_energy())
        assert (es[0] - es[1]) / (2 * h) / a.get_volume() == pytest.approx(stress[i, j], abs=1e-9)
    assert abs(stress[1, 2]) > 1e-4
    # at its targets the restraint is nothing
    b = sheet.copy()
    b.calc = tool.restraint_calculator([quads[q] for q in led], [a1[q] for q in led], 10.0)
    assert b.get_potential_energy() == pytest.approx(0.0, abs=1e-12) and np.abs(b.get_forces()).max() < 1e-6


def _scan(name, start_torsions, top, top_at, arrived, end_itself, off_the_end=5.0, points=6):
    images = [{"image": f"i{r:02d}", "energy": (top if r == top_at else 0.1 * r)} for r in range(points + 1)] + [{"image": "released", "energy": arrived}]
    return {"name": name, "start_torsions_deg": start_torsions, "top": {"image": f"i{top_at:02d}", "energy": top, "index": top_at}, "images": images,
            "arrived": {"energy": arrived, "end_in_the_potential": end_itself, "torsions_from_the_end_deg": [off_the_end, -off_the_end]}}


def test_the_scan_chosen_arrived_where_it_was_led_and_has_the_lower_top_above_all_trans(tmp_path, capsys):
    tool = _tool()
    flat, kinked = [-165.0, 163.0], [-64.0, 50.0]
    # from all-trans: top 6.0; from the sheet (5.0 below all-trans): top 10.2 above the sheet, 5.1 above the all-trans cell it arrives at
    cases = {"both arrive, the one from the sheet is lower": (_scan("p_f", flat, 6.0, 3, -5.2, -5.0), _scan("p_r", kinked, 10.2, 4, 5.1, 5.0), "p_r", ["p_r_i04", "p_r_i03", "p_r_i05"]),
             "the lower one did not arrive": (_scan("p_f", flat, 6.0, 3, -5.2, -5.0), _scan("p_r", kinked, 10.2, 4, 5.1, 5.0, off_the_end=80.0), "p_f", ["p_f_i03", "p_f_i02", "p_f_i04"]),
             "the top is the first led image": (_scan("p_f", flat, 6.0, 1, -5.2, -5.0), _scan("p_r", kinked, 14.0, 4, 5.1, 5.0), "p_f", ["p_f_i01", "p_f_i02"]),
             "the top is the last image": (_scan("p_f", flat, 6.0, 6, -5.2, -5.0), _scan("p_r", kinked, 14.0, 4, 5.1, 5.0), "p_f", ["p_f_i06", "p_f_i05"])}
    for what, (forward, back, chosen, names) in cases.items():
        for t in (forward, back):
            (tmp_path / f"{t['name']}.json").write_text(json.dumps(t))
        tool.choose(SimpleNamespace(scans=[str(tmp_path / "p_f.json"), str(tmp_path / "p_r.json")]))
        words = capsys.readouterr().out.split()
        assert words[0] == str(tmp_path / f"{chosen}.json") and words[1:] == names, what
    assert tool._above_all_trans(forward) == (6.0, True) and tool._above_all_trans(back)[0] == pytest.approx(14.0 - 5.1)


def test_the_queues_fit_line_names_files_and_frames_that_are_there(frames):
    lines = (ROOT / "examples" / "gpu_jobs" / "queue_20261121_path.txt").read_text().splitlines()
    fit = next(line for line in lines if line.startswith("refit.sh ")).split()
    at = {w: fit.index(w) for w in ("frames", "prepare", "train", "also")}
    files = fit[at["frames"] + 1:at["prepare"]]
    assert len(files) == 24 and all((ROOT / f).exists() for f in files) and files[-1].endswith("cells14_frames_pbe_d3bj.extxyz")
    ten = (ROOT / "examples" / "gpu_jobs" / "queue_20261118_mlip.txt").read_text().splitlines()[-1].split()
    assert files[:10] == ten[ten.index("frames") + 1:ten.index("train")]   # all_s1's files first, in its order: its held-back frames are held back again
    assert not any("row_cells_lone_chain" in f for f in files)
    prepare = fit[at["prepare"] + 1:at["train"]]
    import re

    heavy = prepare[prepare.index("--heavy") + 1]
    assert sorted(n for n in frames if re.fullmatch(heavy, n)) == sorted(frames)                     # the 23 cells, and
    from ase.io import read

    others = [str(a.info["frame"]) for f in files[:-1] for a in read(ROOT / f, index=":")]
    assert len(others) == 598 and not any(re.fullmatch(heavy, n) for n in others)                    # nothing else
    held = prepare[prepare.index("--hold-out-frames") + 1:]
    assert len(held) == 5 and set(held) <= set(frames)
    for line in lines:
        if line.startswith("refit_score.sh "):
            assert line.split()[2:] == held
    for section in ("train", "also"):
        words = fit[at[section] + 1:at["also"]] if section == "train" else fit[at["also"] + 1:]
        assert words[:2] == ["--start", "docs/predictions/mlip/models/all_s1.model"] and "--heavy-weight" in words
    for line in lines:
        if line.startswith("path_card.sh "):
            for spec in line.split()[3:]:
                _, start, end = spec.split(":")
                assert start in frames and end in frames


SCANS = ("f", "r", "a1", "a2", "b1", "b2")
# each scan: its top (image, kcal/mol per chain above its first image), where its released last image rests, and how far its led torsions rest from the end (degrees)
STORED = {"film_f": ("i06", 10.31, 2.56, 9.4), "film_r": ("i07", 8.37, 5.13, 4.4), "film_a1": ("i07", 7.70, 3.86, 10.0), "film_a2": ("i11", 4.55, 3.66, 16.3),
          "film_b1": ("i07", 7.46, 1.31, 15.5), "film_b2": ("i07", 9.00, 4.39, 0.2), "leader_f": ("i07", 5.62, 1.50, 10.6), "leader_r": ("i04", 8.68, 3.08, 7.5),
          "leader_a1": ("i06", 7.01, 0.87, 4.9), "leader_a2": ("i06", 6.39, 3.53, 1.6), "leader_b1": ("i08", 5.67, 0.85, 10.9), "leader_b2": ("i07", 6.41, 2.41, 1.8)}
ARRIVE_ALL_TRANS = ("film_r", "film_a2", "film_b2", "leader_r", "leader_a2", "leader_b2")   # every backbone torsion trans at rest


def _stored(name):
    return json.loads((MLIP / "paths" / f"{name}.json").read_text())


def test_the_stored_scans_tops_and_where_each_arrived():
    from ase.io import read

    for name, (top, energy, rest, off) in STORED.items():
        t = _stored(name)
        assert t["potential"] == "all_s1" and t["restraint_eV"] == 10.0 and t["fmax_eV_per_A"] == 0.02
        rows = t["images"]
        assert [r["image"] for r in rows] == [f"i{i:02d}" for i in range(13)] + ["released"] and all(r["converged"] for r in rows)
        assert rows[0]["energy"] == pytest.approx(0.0, abs=1e-9)                        # every energy is against the scan's first image, as labelled
        assert (t["top"]["image"], t["top"]["energy"]) == (top, pytest.approx(energy, abs=0.006)) and energy == pytest.approx(max(r["energy"] for r in rows[:13]), abs=0.006)
        assert t["arrived"]["energy"] == pytest.approx(rest, abs=0.006) and t["arrived"]["converged"]
        assert np.abs(t["arrived"]["torsions_from_the_end_deg"]).max() == pytest.approx(off, abs=0.06)
        frames = read(MLIP / "paths" / f"{name}.extxyz", index=":")
        tool = _tool()
        at_rest = tool.torsions_deg(frames[-1], tool.backbone_torsions(frames[-1]))
        assert (np.abs(at_rest).min() > 150) == (name in ARRIVE_ALL_TRANS), name
        assert [str(a.info["frame"]) for a in frames] == [f"{name}_{r['image']}" for r in rows]
        e0 = float(frames[0].info["potential_eV_per_monomer"])
        assert [(float(a.info["potential_eV_per_monomer"]) - e0) * 14 * 23.0605 for a in frames] == pytest.approx([r["energy"] for r in rows], abs=2e-3)
    # from upright all-trans the kink's torsions arrive at a kinked cell that is still upright; from the sheet at an all-trans cell still sheared
    for cell, back in (("film", (0.96, -0.96)), ("leader", (-0.47, -1.07))):
        assert np.abs(_stored(f"{cell}_f")["arrived"]["shear_along_chains_monomers"]).max() < 0.21
        assert _stored(f"{cell}_r")["arrived"]["shear_along_chains_monomers"] == pytest.approx(back, abs=0.006)
    # one torsion of the pair led alone: the other gauche bond moves a monomer along and the kink is stretched, G T T T G'
    for cell in ("film", "leader"):
        one, rest = _stored(f"{cell}_b1"), _stored(f"{cell}_b2")
        assert one["only"] == [1, 3] and len(one["led_torsions"]) == 2
        first = [q[0] for q in rest["led_torsions"]]
        a = np.array(rest["start_torsions_deg"])
        gauche = [i for i in range(len(first) // 2) if abs(a[i]) < 75]
        assert len(gauche) == 2 and np.sign(a[gauche[0]]) != np.sign(a[gauche[1]])
        assert first[gauche[1]] - first[gauche[0]] == 12                                 # four backbone bonds apart (three atoms a backbone carbon here), where the kink's are two
        kink = _stored(f"{cell}_r")["led_torsions"]
        assert kink[1][0] - kink[0][0] == 6
    # the fine scans start at the image before the top and find it higher: 8.55 and 9.18
    assert [(_stored(f"{c}_r_fine")["top"]["energy"]) for c in ("film", "leader")] == pytest.approx([8.55, 9.18], abs=0.006)
    # released with no restraint, the images beside the film composition's top fall to its two ends; the leader's image after its top falls to neither
    assert [v["fell_to"] for v in _stored("film_r_free")["images"].values()] == ["start", "end"]
    assert [v["fell_to"] for v in _stored("leader_r_free")["images"].values()] == ["start", "neither"]


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-q"]))
