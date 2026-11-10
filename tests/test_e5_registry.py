"""E5, the axial registry of beta-PVDF: the cells built from Dipole's accepted record (examples/axial_registry.py)
and the ingester that recognises Dipole's records by their geometry (examples/crosscheck_e5.py)."""
import importlib.util
import json
from pathlib import Path

import numpy as np

from polyfind.fitting import PRESET_PARAMETERS

ROOT = Path(__file__).resolve().parents[1]


def _load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "examples" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


AR = _load("axial_registry")
CX = _load("crosscheck_e5")
IDX = json.load(open(ROOT / "docs" / "predictions" / "e5_structures" / "e5_structures.json"))
PRED = json.load(open(ROOT / "docs" / "predictions" / "e5_polyfind_belt4.json"))
EXPECTED = {"beta_a_slide_c8": ("a", "c/8"), "beta_a_slide_c4": ("a", "c/4"), "beta_a_slide_3c8": ("a", "3c/8"), "beta_a_slide_c2": ("a", "c/2"),
            "beta_2a_reference": ("c0", None), "beta_2a_slide_c4": ("c", "c/4"), "beta_2a_slide_c2": ("c", "c/2")}


def _reference():
    r = IDX["reference"]
    return {"case": IDX["reference_case"], "cell_A": np.array(r["cell_A"]), "positions_A": np.array(r["positions_A"]), "symbols": list(r["symbols"]),
            "energy_Eh": IDX["reference_energy_Eh"], "result_sha256": IDX["reference_result_sha256"]}


def test_the_slid_cells_differ_from_the_accepted_cell_only_in_the_moved_chain():
    """Part a moves the six atoms of the chain at (1/2, 1/2) along z by the stated fraction of c and
    nothing else; part c doubles the cell along its short axis and moves the first chain alone.  The
    committed index holds exactly these cells."""
    ref = _reference()
    cells = AR.slid_structures(ref)
    assert set(cells) == set(EXPECTED) == set(IDX["structures"])
    c = ref["cell_A"][2, 2]
    for name, f in (("beta_a_slide_c8", 0.125), ("beta_a_slide_c4", 0.25), ("beta_a_slide_3c8", 0.375), ("beta_a_slide_c2", 0.5)):
        d = cells[name]["positions_A"] - ref["positions_A"]
        assert np.abs(d[:6]).max() == 0.0 and np.abs(d[6:, :2]).max() == 0.0
        assert np.allclose(d[6:, 2], f * c, atol=1e-12)
        assert np.array_equal(cells[name]["cell_A"], ref["cell_A"])
    doubled = np.concatenate([ref["positions_A"], ref["positions_A"] + ref["cell_A"][0]])
    assert np.array_equal(cells["beta_2a_reference"]["positions_A"], doubled)
    for name, f in (("beta_2a_slide_c4", 0.25), ("beta_2a_slide_c2", 0.5)):
        s = cells[name]
        d = s["positions_A"] - doubled
        assert len(s["symbols"]) == 24 and np.allclose(s["cell_A"][0], 2 * ref["cell_A"][0]) and np.array_equal(s["cell_A"][1:], ref["cell_A"][1:])
        assert np.allclose(d[:6, 2], f * c, atol=1e-12) and np.abs(d[6:]).max() == 0.0 and np.abs(d[:6, :2]).max() == 0.0
    for name, s in cells.items():
        assert np.allclose(IDX["structures"][name]["initial"]["positions_A"], s["positions_A"], atol=1e-12)
        assert IDX["structures"][name]["moved_atoms"] == s["moved_atoms"]


def test_every_cell_is_recognised_by_its_geometry_alone():
    """The ingester takes no names: each built cell is classified from its bonds and the heights of its
    chains, also with its atoms shuffled and wrapped into the cell, and the c/2 cell in another lattice
    (b opened, a closed, fractional coordinates kept) is part b."""
    ref = _reference()
    assert CX.classify({"cell_A": ref["cell_A"], "positions_A": ref["positions_A"], "symbols": ref["symbols"]}, ref["cell_A"]) == ("zero", None)
    rng = np.random.default_rng(0)
    for name, want in EXPECTED.items():
        init = IDX["structures"][name]["initial"]
        assert CX.classify(init, ref["cell_A"]) == want, name
        cell, pos, sym = np.array(init["cell_A"]), np.array(init["positions_A"]), list(init["symbols"])
        order = rng.permutation(len(sym))
        frac = pos[order] @ np.linalg.inv(cell)
        shuffled = {"cell_A": cell, "positions_A": (frac - np.floor(frac)) @ cell, "symbols": [sym[i] for i in order]}
        assert CX.classify(shuffled, ref["cell_A"]) == want, name
    init = IDX["structures"]["beta_a_slide_c2"]["initial"]
    cell = np.array(init["cell_A"])
    relaxed = cell * np.array([[0.995], [1.07], [1.0]])
    frac = np.array(init["positions_A"]) @ np.linalg.inv(cell)
    # an affine change of the cell stretches the chains across b; the bonds along the chain survive, which is what the reader uses
    assert CX.classify({"cell_A": relaxed, "positions_A": frac @ relaxed, "symbols": init["symbols"]}, ref["cell_A"]) == ("b", "c/2")


def test_the_ingester_returns_the_energies_it_is_given():
    """Records carrying the accepted energy plus polyfind's own registered differences come back as those
    differences, per monomer of the cell, with the zero of part c taken from the unslid doubled record."""
    ref = _reference()
    own = PRED["own_cell"]
    e_ref = ref["energy_Eh"]
    records = [{"name": "zero", "energy_Eh": e_ref, "cell_A": ref["cell_A"], "positions_A": ref["positions_A"], "symbols": ref["symbols"]}]
    for name, (part, key) in EXPECTED.items():
        init = IDX["structures"][name]["initial"]
        if part == "a":
            e = e_ref + own["a"]["slides"][key]["dE_cell_mEh"] * 1e-3
        elif part == "c0":
            e = 2 * e_ref + 3e-9  # a doubled cell's SCF differs from twice the single cell's by noise; the zero of c absorbs it
        else:
            e = 2 * e_ref + 3e-9 + own["c"]["slides"][key]["dE_cell_mEh"] * 1e-3
        records.append({"name": name, "energy_Eh": e, **init})
    out = CX.crosscheck(PRED, IDX["structures"], records, log=lambda s: None)
    for key in ("c/8", "c/4", "3c/8", "c/2"):
        assert abs(out["a"][key]["dipole"] - own["a"]["slides"][key]["per_cell_monomer_kcal"]) < 1e-6 and out["a"][key]["rising"]
    for key in ("c/4", "c/2"):
        assert abs(out["c"][key]["dipole"] - own["c"]["slides"][key]["per_cell_monomer_kcal"]) < 1e-6
    assert abs(out["anisotropy"]["dipole"] - own["anisotropy_c_over_a_per_moved_monomer_at_half"]) < 1e-6
    # the plain table: names instead of geometries, part b with its cell
    table = {"structures": {"beta_a_slide_c2": {"energy_Eh": e_ref + 2e-3}, "beta_b_relaxed": {"energy_Eh": e_ref + 1e-3, "cell_A": (ref["cell_A"] * np.array([[0.99], [1.05], [1.0]])).tolist()}}}
    path = ROOT / "tests" / "data" / "_e5_table_tmp.json"
    path.write_text(json.dumps(table))
    try:
        out = CX.crosscheck(PRED, IDX["structures"], CX.read_records([str(path)], IDX["structures"]), log=lambda s: None)
    finally:
        path.unlink()
    assert abs(out["a"]["c/2"]["dipole"] - 2e-3 * CX.KCAL_PER_EH / 2) < 1e-9 and abs(out["b"]["dipole"] - 1e-3 * CX.KCAL_PER_EH / 2) < 1e-9
    assert abs(out["b"]["cell_change"]["b"] - 0.05) < 1e-9 and abs(out["b"]["cell_change"]["a"] + 0.01) < 1e-9


def test_the_registry_is_even_and_periodic_in_the_slide():
    """E5 asks for half a period because the all-trans chain has a mirror plane normal to its axis: sliding
    one sublattice by +u, by -u and by u plus one monomer period costs the same.  The registered numbers
    rise along the slide and one chain of four costs more per moved monomer than one sublattice."""
    letters, _ = AR.DK.sequence(AR.ROW)
    with PRESET_PARAMETERS["pvdf-dft-valence-belt4"].applied():
        polymer = AR.DK.C.make_polymer(AR.ROW, letters)
        stored = json.load(open(ROOT / "docs" / "predictions" / "dilute_kink_belt4.json"))["all_trans_references"][AR.ROW]
        d = AR.DK.all_trans_reference(polymer, AR.ROW, stored)
        model, x, ch = AR.DK.make_model(polymer, d, 1, 1)
        cm = float(model.state(x).lattice[3]) / ch.n_monomers
        assert abs(AR.slip(model, x, cm, 0, 1)) < 1e-6
        e = AR.rigid_slides(model, x, 0, cm, (0.25, -0.25, 1.25), ch.n_monomers, 2)
        y = x.copy()
        y[model.index("z", 0)[0]] += 0.25 * cm
        assert abs(AR.slip(model, y, cm, 0, 1) - 0.25 * cm) < 1e-6
    v = [e[k]["per_cell_monomer_kcal"] for k in ("0.2500", "-0.2500", "1.2500")]
    assert abs(v[0] - v[1]) < 1e-6 and abs(v[0] - v[2]) < 1e-6
    assert abs(v[0] - PRED["own_cell"]["a"]["slides"]["c/4"]["per_cell_monomer_kcal"]) < 1e-5
    for block in (PRED["own_cell"], PRED["dipole_cell"]):
        a = [block["a"]["slides"][k]["per_cell_monomer_kcal"] for k in ("c/8", "c/4", "3c/8", "c/2")]
        assert all(q > p > 0 for p, q in zip([1e-9] + a[:-1], a))
        assert block["anisotropy_c_over_a_per_moved_monomer_at_half"] > 1.5
