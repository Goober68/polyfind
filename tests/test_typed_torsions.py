"""One torsion triple a bond type on the lattice side (`FFParameters.torsion_types`).

Every earlier preset gives every backbone bond of every chain the CF2-CH2 type's triple.  A
reversed unit has a CF2-CF2 and a CH2-CH2 bond and TrFE a CF2-CHF one, and the valence fit has a
triple for each.  Held here: which bond a dihedral turns, the typing of a polymer's repeat, the
packer's energy and its gradient, both forms of the mixed cell, and that nothing changes when no
type has its own triple.
"""

import importlib.util
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from polyfind import fitting as F
from polyfind import mixedcell as MC
from polyfind import pack as pack_mod
from polyfind.chain import build_chain
from polyfind.fitting import PRESET_PARAMETERS
from polyfind.linegroup import repeat_chains

ROOT = Path(__file__).resolve().parents[1]
BELT4 = "pvdf-dft-valence-belt4"
TYPES = (("CF2-CHF", -0.3, -0.5, 0.8), ("CF2-CF2", -1.7, 0.1, 1.8), ("other", -0.1, -0.6, 1.5))


@pytest.fixture(scope="module")
def DK():
    mp = pytest.MonkeyPatch()
    mp.chdir(ROOT)
    spec = importlib.util.spec_from_file_location("dilute_kink", ROOT / "examples" / "dilute_kink.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    yield mod
    mp.undo()


def _dihedral(p0, p1, p2, p3):
    b0, b1, b2 = p1 - p0, p2 - p1, p3 - p2
    n1, n2 = np.cross(b0, b1), np.cross(b1, b2)
    return np.degrees(np.arctan2(np.linalg.norm(b1) * (b0 @ n2), n1 @ n2))


def test_a_dihedral_turns_the_bond_one_on_from_its_index(DK):
    """Dihedral ``j`` of a built chain is the torsion about the bond between backbone atoms
    ``j + 1`` and ``j + 2`` of the cycled repeat, read off built coordinates; and the typing of
    the repeat follows it."""
    polymer = DK.C.make_polymer("vdf1-rvdf1-2", "VI")      # CH2, CF2, CF2, CH2
    B = polymer.bonds_per_repeat
    assert ["".join(sorted(b.substituents)) for b in polymer.backbone] == ["HH", "FF", "FF", "HH"]
    for j in range(B):
        tors = np.full(B, 180.0)
        tors[j] = 65.0
        s = build_chain(polymer, np.tile(tors, 3), cap=False)
        X = s.coords[np.asarray(s.backbone)]
        gauche = [k for k in range(1, len(X) - 2) if abs(abs(_dihedral(X[k - 1], X[k], X[k + 1], X[k + 2])) - 65.0) < 1.0]
        assert gauche and {(k % B, (k + 1) % B) for k in gauche} == {((j + 1) % B, (j + 2) % B)}
    names = [F.REF_TORSION_TYPES[k] for k in F.polymer_torsion_types(polymer, B)]
    assert names == ["CF2-CF2", "CF2-CH2", "other", "CF2-CH2"]
    assert list(F.polymer_torsion_types(DK.C.make_polymer("vdf1", "V"), 6)) == [0] * 6
    trfe = [F.REF_TORSION_TYPES[k] for k in F.polymer_torsion_types(DK.C.make_polymer("vdf1-trfe1-2", "VR"), 4)]
    # V then R is CH2, CF2, CHF, CF2: its bonds are CH2-CF2, CF2-CHF, CHF-CF2 and CF2-CH2
    assert sorted(trfe) == ["CF2-CH2", "CF2-CH2", "CF2-CHF", "CF2-CHF"]


def test_the_parameter_set_gives_each_dihedral_its_bonds_triple(DK):
    base = PRESET_PARAMETERS[BELT4]
    polymer = DK.C.make_polymer("vdf5-rvdf1-6", "VVVVIV")
    P = polymer.bonds_per_repeat
    assert base.torsion_for(polymer, P) == tuple(base.torsion)          # no type has its own: the one triple, as ever
    typed = replace(base, torsion_types=TYPES)
    V = typed.torsion_for(polymer, P)
    assert V.shape == (P, 3) and "CF2-CF2 torsion" in typed.describe()
    own = {t: (a, b, c) for t, a, b, c in TYPES}
    for j, k in enumerate(F.polymer_torsion_types(polymer, P)):
        name = F.REF_TORSION_TYPES[k]
        assert tuple(V[j]) == (tuple(base.torsion) if name == "CF2-CH2" else own[name])
    assert sorted(F.REF_TORSION_TYPES[k] for k in F.polymer_torsion_types(polymer, P)).count("CF2-CF2") == 1
    with pytest.raises(ValueError, match="unknown torsion types"):
        replace(base, torsion_types=(("CF2-CCl2", 0.0, 0.0, 0.0),)).torsion_for(polymer, P)


def _cell(DK, params, row="vdf1-rvdf1-2", letters="VI", tors=(64.0, 178.0, -66.0, 171.0)):
    with params.applied():
        polymer = DK.C.make_polymer(row, letters)
        P = polymer.bonds_per_repeat
        angs = np.array([polymer.backbone[k].backbone_angle for k in range(P)], dtype=float)
        ch = repeat_chains(polymer, "t", np.array(tors, dtype=float)[None], angs)[0]
        pk = pack_mod.CrystalPacker(ch, n_chains=2, **DK.KW)
        st = MC.MixedState(np.array([5.2, 9.4, 90.0, ch.c]), np.array([[0.0, 0.0, 10.0, 0.0], [0.5, 0.5, 170.0, 0.4]]),
                           np.tile(np.array(tors, dtype=float), (2, 1)), np.tile(angs, (2, 1)), np.zeros(2))
        model = MC.MixedModel(polymer, "t", pk, np.array([False, False]), [ch, ch], angle_stiffness=105.0, shear=True)
        return polymer, pk, model, st.vector()


def test_the_packer_and_the_mixed_cell_use_the_typed_terms(DK):
    base = PRESET_PARAMETERS[BELT4]
    typed = replace(base, torsion_types=TYPES)
    tors = np.array([64.0, 178.0, -66.0, 171.0])  # gauche about the CF2-CF2 bond and about the CH2-CH2 one
    polymer, pk0, m0, x = _cell(DK, base)
    _, pk1, m1, _ = _cell(DK, typed)
    assert isinstance(pk0.torsion, tuple) and pk1.torsion.shape == (4, 3)
    # the energy is each dihedral's own Fourier series
    V = typed.torsion_for(polymer, 4)
    ph = np.deg2rad(tors)
    by_hand = float((0.5 * (V[:, 0] * (1 + np.cos(ph)) + V[:, 1] * (1 - np.cos(2 * ph)) + V[:, 2] * (1 + np.cos(3 * ph)))).sum())
    assert pk1.torsion_energy(tors) == pytest.approx(by_hand, abs=1e-12)
    assert abs(pk1.torsion_energy(tors) - pk0.torsion_energy(tors)) > 0.05
    g = pk1.torsion_energy_grad(tors)
    for j in range(4):
        up, dn = tors.copy(), tors.copy()
        up[j] += 1e-4
        dn[j] -= 1e-4
        assert (pk1.torsion_energy(up) - pk1.torsion_energy(dn)) / 2e-4 == pytest.approx(g[j], abs=1e-8)
    with pytest.raises(ValueError, match="one torsion triple for each of 4"):
        pk1.torsion_energy(np.tile(tors, 2))
    # the mixed cell: the two parameter sets differ by exactly the torsion term, in both of its forms, and the gradient is the energy's
    for xp in (None, np):
        m0.xp, m1.xp = xp, xp
        p0, p1 = {}, {}
        e0, _ = m0.value_and_grad(x, p0)
        e1, g1 = m1.value_and_grad(x, p1)
        assert p1["torsion"] == pytest.approx(2 * by_hand, abs=1e-10) and p0["torsion"] == pytest.approx(2 * pk0.torsion_energy(tors), abs=1e-10)
        assert e1 - e0 == pytest.approx(p1["torsion"] - p0["torsion"], abs=1e-8)
        for k in m1.index("tors", 0):
            up, dn = x.copy(), x.copy()
            up[k] += 1e-5
            dn[k] -= 1e-5
            assert (m1.value_and_grad(up)[0] - m1.value_and_grad(dn)[0]) / 2e-5 == pytest.approx(g1[k], abs=2e-5 * max(1.0, abs(g1[k])))
    # every type given the one triple is the untyped set
    same = replace(base, torsion_types=tuple((t,) + tuple(base.torsion) for t in ("CF2-CHF", "CF2-CF2", "other")))
    _, _, m2, _ = _cell(DK, same)
    assert m2.value_and_grad(x)[0] == pytest.approx(m0.value_and_grad(x)[0], abs=1e-10)


def test_a_units_own_bonds_have_names_and_take_a_triple_put_on_a_preset():
    """CH2-CFCl and the like are named apart from the catch-all; a triple put on a preset lands on those bonds only."""
    import importlib.util
    from pathlib import Path

    import numpy as np

    from polyfind import fitting as F

    spec = importlib.util.spec_from_file_location("candidate_screen", Path(__file__).resolve().parents[1] / "examples" / "candidate_screen.py")
    C = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(C)
    base = F.PRESET_PARAMETERS["pvdf-dft-chain1"]
    with base.applied():
        names = {letters: F.polymer_lattice_torsion_types(C.make_polymer("t", letters), 4) for letters in ("VC", "VK", "VD", "VA", "VY", "VX", "VE", "VR")}
        cfe = C.make_polymer("t", "VC")
    assert names["VC"] == ["CF2-CH2", "CH2-CFCl", "CH2-CFCl", "CF2-CH2"]
    assert names["VK"] == ["CF2-CF2", "CF2-CFCl", "CH2-CFCl", "CF2-CH2"]
    assert names["VD"] == ["CF2-CHCl", "CF2-CHCl", "CF2-CH2", "CF2-CH2"]
    assert names["VA"][1:3] == ["CH2-CHCN", "CH2-CHCN"] and names["VY"][1:3] == ["CH2-CHF", "CH2-CHF"] and names["VX"][1:3] == ["CHF-CHF", "CH2-CHF"]
    assert names["VE"][1:3] == ["CH2-CH2", "CH2-CH2"] and names["VR"] == ["CF2-CHF", "CF2-CHF", "CF2-CH2", "CF2-CH2"]
    assert set(sum(names.values(), [])) <= set(F.LATTICE_TORSION_TYPES)
    # chain1 has no triple for a CH2-CFCl bond: it takes the CF2-CH2 bond's
    plain = np.array(base.torsion)
    assert np.allclose(base.torsion_for(cfe, 4), [plain] * 4)
    typed = F.typed_on(base, {"CH2-CFCl": (0.1, -0.2, 0.3)})
    assert np.allclose(typed.torsion_for(cfe, 4), [plain, (0.1, -0.2, 0.3), (0.1, -0.2, 0.3), plain])
    assert dict((t, v) for t, *v in typed.torsion_types)["CF2-CF2"] == list(dict((t, v) for t, *v in base.torsion_types)["CF2-CF2"])
    with pytest.raises(ValueError):
        F.typed_on(base, {"CH2-CBr2": (0.0, 0.0, 0.0)})


def test_a_parameter_set_is_read_from_a_name_or_a_file(tmp_path):
    import json

    from polyfind import fitting as F

    assert F.parameters_of("pvdf-dft-chain1") is F.PRESET_PARAMETERS["pvdf-dft-chain1"]
    f = tmp_path / "typed.json"
    f.write_text(json.dumps({"base": "pvdf-dft-chain1", "torsion_types": {"CH2-CFCl": [0.1, -0.2, 0.3]}}))
    prm = F.parameters_of(f)
    assert ("CH2-CFCl", 0.1, -0.2, 0.3) in prm.torsion_types and prm.lj == F.PRESET_PARAMETERS["pvdf-dft-chain1"].lj
    g = tmp_path / "stage.json"
    g.write_text(json.dumps({"x": list(F.CRYSTAL13_CHAIN_FITTED_X)}))
    assert F.parameters_of(g).torsion_types == F.PRESET_PARAMETERS["pvdf-dft-chain1"].torsion_types


def test_a_units_own_wall_and_bond_charge_go_on_a_preset(tmp_path):
    import json

    from polyfind import fitting as F

    base = F.PRESET_PARAMETERS["pvdf-dft-chain1-cfe"]
    assert "N" not in {e for e, *_ in base.lj}  # nitrogen is at its UFF values in the preset
    prm = F.with_terms(base, lj={"N": (3.5, 0.08)}, charge_increments={"C-N": 0.55})
    assert ("N", 3.5, 0.08) in prm.lj and ("C", "N", 0.55) in prm.charge_increments
    assert [e for e, *_ in prm.lj][:-1] == [e for e, *_ in base.lj] and len(prm.charge_increments) == len(base.charge_increments)
    assert prm.torsion_types == base.torsion_types and prm.wall == base.wall and F.with_terms(base) == base
    # a bond named the other way round moves the opposite charge
    assert F.with_terms(base, charge_increments={"N-C": -0.55}).charge_increments == prm.charge_increments
    assert F.with_terms(prm, lj={"N": (3.7, 0.05)}).lj[-1] == ("N", 3.7, 0.05)
    f = tmp_path / "unit.json"
    f.write_text(json.dumps({"base": "pvdf-dft-chain1-cfe", "lj": {"N": [3.5, 0.08]}, "charge_increments": {"C-N": 0.55}, "torsion_types": {}}))
    assert F.parameters_of(f) == prm
    g = tmp_path / "both.json"
    g.write_text(json.dumps({"base": str(f), "torsion_types": {"CH2-CHCN": [0.1, 0.2, 0.3]}}))
    assert F.parameters_of(g) == F.typed_on(prm, {"CH2-CHCN": (0.1, 0.2, 0.3)})
