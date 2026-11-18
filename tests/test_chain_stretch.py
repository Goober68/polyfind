"""The lone chain's line (examples/chain_stretch.py): the measurement of a frame's backbone, and polyfind's side."""

import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

ase = pytest.importorskip("ase")
ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("chain_stretch", ROOT / "examples" / "chain_stretch.py")
ST = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(ST)


def _helix(gauche):
    """Polyfind's TGTG' chain of plain PVDF with this gauche angle, as a frame in a box."""
    from ase import Atoms

    from polyfind.pack import periodic_chain_from_torsions
    from polyfind.polymers import PVDF

    ch = periodic_chain_from_torsions(PVDF, "helix", [180.0, gauche, 180.0, -gauche])
    keep = [k for k, e in enumerate(ch.elements) if e in ("C", "H", "F")]
    P = ch.coords[keep] + [8.0, 8.0, 0.0]
    return Atoms(symbols=[ch.elements[k] for k in keep], positions=P, cell=[16.0, 16.0, ch.c], pbc=True), ch


@pytest.mark.parametrize("gauche", [45.0, 60.0, 85.0])
def test_a_frames_backbone_dihedrals_are_the_torsions_the_chain_was_built_with(gauche):
    a, ch = _helix(gauche)
    d = np.abs(ST.dihedrals(a)[:4])
    assert sorted(d) == pytest.approx([gauche, gauche, 180.0, 180.0], abs=0.5)
    a.wrap()  # wherever the atoms sit in the box
    assert sorted(np.abs(ST.dihedrals(a)[:4])) == pytest.approx([gauche, gauche, 180.0, 180.0], abs=0.5)


def test_the_stored_line_is_a_well_about_the_lone_helix_and_dft_was_where_the_tilted_step_said():
    d = json.loads((ROOT / "docs" / "predictions" / "chain_stretch_pbe_d3bj.json").read_text())
    helix, trans = d["branches"]["helix"], d["branches"]["trans"]
    assert helix["monomers"] == 2 and trans["monomers"] == 1
    pts = {round(q["repeat"], 2): q for q in helix["points"]}
    assert pts[4.65]["energy"] == 0.0 and pts[4.55]["energy"] == pytest.approx(pts[4.75]["energy"], abs=0.02) and 0.2 < pts[4.75]["energy"] < 0.3
    # the helix keeps its shape: the gauche angle opens by 5 degrees for every 2% of length
    gauche = [float(np.abs(pts[c]["dihedrals"]).min()) for c in (4.55, 4.65, 4.75, 4.85, 4.95)]
    assert gauche == pytest.approx([50, 56, 61, 65, 70], abs=1.5)
    # it never lets go: at the all-trans chain's length it is still a helix, two kcal/mol above that chain
    assert float(np.abs(pts[5.18]["dihedrals"]).min()) < 95 and pts[5.18]["energy"] > 5.0
    low = min(trans["points"], key=lambda q: q["energy"])
    assert low["repeat"] == pytest.approx(2.591, abs=0.002) and low["energy"] == pytest.approx(3.449, abs=0.005)
    # the Newton step's prediction of DFT's energy, and what DFT's forces had left
    every = helix["points"] + trans["points"]
    assert max(abs(q["predicted_minus_dft"]) for q in every) < 0.04 and max(q["largest_force"] for q in every) < 0.06


def test_polyfinds_lone_chain_at_a_repeat_is_closed_at_that_repeat_and_walks_from_the_cells():
    """The line's polyfind side under belt4: the chain held at each repeat really has it, the
    gauche angle opens as the chain is pulled, and its energy is lowest near the cell's repeat."""
    from polyfind.fitting import PRESET_PARAMETERS

    IF = ST._load("interchain_frames")
    table = json.loads((ROOT / "docs" / "predictions" / "interchain_belt4.json").read_text())["packings"]
    with PRESET_PARAMETERS["pvdf-dft-valence-belt4"].applied():
        model, x0, ch, _ = IF.model_of("alpha", table["alpha"].get("cell"), check=False)
        line = ST.polyfind_line(IF.SP, model, x0, ch, [4.55, 4.65, 4.75, 4.85], tilts=False)
    assert max(q["c_mismatch"] for q in line) < 1e-4
    gauche = [float(np.abs(np.mod(np.array(q["torsions"]) + 180.0, 360.0) - 180.0).min()) for q in line]
    assert gauche == sorted(gauche) and 35 < gauche[0] < gauche[-1] < 85
    e = np.array([q["energy"] for q in line]) / ch.n_monomers
    assert e[3] - e[1] > 1.0 and abs(e[0] - e[1]) < 0.2
