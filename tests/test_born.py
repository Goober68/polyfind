"""Born effective charges: the dipole derivative the model actually has, checked against what it must be.

A fixed-charge cell's Born tensor is the charge times the identity; with charge flux it is the
analytic charge gradient contracted with the positions; with two chains it sums to the
derivative of the packer's own dipole with respect to the repeat coordinate, and the two
chains' tensors are related by the operation that maps one chain onto the other; the
acoustic sum vanishes identically; and the stretch channel's zero on the lattice path is the
built bond length, so a ``k_bond`` changes the Born tensor and nothing at rest.
"""
import json

import numpy as np
import pytest

from polyfind import pack as pack_mod
from polyfind.born import atom_labels, born_charges, canonical, load_reference
from polyfind.ewald import EwaldSpec
from polyfind.forcefield import SimpleFF
from polyfind.pack import CrystalPacker, periodic_chain
from polyfind.polarizability import Polarizable
from polyfind.polymers import PVDF, THREE_STATE

T = 0
P0 = np.array([4.6, 8.6, 90.0, 0.0, 0.0, 0.4, 0.0])


def _flux_ff(flux):
    import polyfind.fitting  # noqa: F401 - registers the presets

    base = SimpleFF.from_preset("pvdf-dft-valence")
    return SimpleFF(**{**{k: v for k, v in base.__dict__.items() if k != "_cache"}, "charge_flux": flux})


def _fitted_packer(flux, n_chains: int = 2, **kw):
    """A beta-PVDF chain carrying ``pvdf-dft-valence``'s charges, and a packer with the given flux.

    ``FFParameters.applied`` rebinds ``pack.CrystalPacker`` to a class that puts the preset's
    charges on the chain, so the packer has to be looked up through the module inside the block.
    """
    from polyfind.fitting import FITTED_VALENCE

    with FITTED_VALENCE.applied():
        ch = periodic_chain(PVDF, (T, T), THREE_STATE)
        pk = pack_mod.CrystalPacker(ch, n_chains=n_chains, valence=SimpleFF.from_preset("pvdf-dft-valence"),
                                    charge_flux=None if flux is None else _flux_ff(flux), **kw)
    return pk.chain, pk


def test_fixed_charges_give_the_charge_times_the_identity():
    pk = CrystalPacker(periodic_chain(PVDF, [T, T], THREE_STATE), n_chains=2)
    b = born_charges(pk, P0)
    for i in range(pk.N):
        assert b.Z[i] == pytest.approx(pk._q_cell[i] * np.eye(3), abs=1e-9)
    assert np.abs(b.asr).max() < 1e-9
    assert b.labels[:6].count("F") == 2 and "C(F2)" in b.labels and "C(H2)" in b.labels


def test_flux_born_tensor_is_the_analytic_charge_gradient():
    """One chain, no rotation: ``Z_i = q_i I + sum_a r_a (dq_a/du_i)^T`` with the closed-form ``dq``."""
    ch, pk = _fitted_packer((("C", "H", -0.7, 0.31), ("C", "F", 0.45, -0.22)), n_chains=1)
    p = np.array([4.6, 8.6, 90.0, 0.0, 0.0, 0.0, 0.0])
    b = born_charges(pk, p)
    q, dq, _ = pk._flux.charges_and_grad(ch.coords, ch.c)
    X = np.asarray(ch.coords, dtype=float)
    for i in range(pk.n):
        assert b.Z[i] == pytest.approx(q[i] * np.eye(3) + X.T @ dq[:, i, :], abs=1e-7)
    assert np.abs(b.asr).max() < 1e-8


def test_two_chains_sum_to_the_repeat_dipole_derivative_with_flux_and_dipoles():
    """Flux, Ewald and induced dipoles on: moving the repeat coordinate moves the atom in both chains."""
    ch, pk = _fitted_packer((("C", "H", -1.2, 0.3), ("C", "F", 0.1, 0.6)),
                            coulomb="ewald", ewald=EwaldSpec(), polarizable=Polarizable())
    b = born_charges(pk, P0)
    n, h = pk.n, 1e-4
    X0 = np.asarray(ch.coords, dtype=float)
    cz = np.array([ch.c])
    for i in range(n):
        for d in range(3):
            Xp, Xm = X0.copy(), X0.copy()
            Xp[i, d] += h
            Xm[i, d] -= h
            num = (pk.dipole(P0[None], coords=Xp[None], c=cz)[0] - pk.dipole(P0[None], coords=Xm[None], c=cz)[0]) / (2 * h)
            assert num == pytest.approx(b.Z[i, :, d] + b.Z[n + i, :, d], abs=1e-6)
    assert np.abs(b.asr).max() < 1e-6


@pytest.mark.parametrize("flip, phi", [(0.0, 30.0), (1.0, 0.0), (1.0, 90.0)])
def test_chain_two_is_the_image_of_chain_one(flip, phi):
    """``Z_2 = M Z_1 M^T`` for ``M`` the operation placing chain 2, where that operation is a symmetry.

    A translation by ``(a + b) / 2`` always is (``dz = 0``); the flip's mirror is one only when
    its line lies along a cell axis, so the flipped cases use ``phi`` = 0 and 90.
    """
    ch, pk = _fitted_packer((("C", "H", -1.2, 0.3), ("C", "F", 0.1, 0.6)),
                            coulomb="ewald", ewald=EwaldSpec(), polarizable=Polarizable())
    p = np.array([4.6, 8.6, 90.0, phi, phi, 0.0, flip])
    b = born_charges(pk, p)
    c, s = np.cos(np.deg2rad(phi)), np.sin(np.deg2rad(phi))
    R = np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
    Mf = np.diag([1.0, -1.0, -1.0]) if flip > 0.5 else np.eye(3)
    M = R @ Mf @ R.T
    n = pk.n
    for i in range(n):
        assert b.Z[n + i] == pytest.approx(M @ b.Z[i] @ M.T, abs=1e-7)
    can = canonical(b, cz=float(ch.c))
    assert set(can.mean) == {"C(F2)", "F", "C(H2)", "H"}
    if phi == 0.0:  # the canonical reflections assume the chain's mirror lies along the cell's y
        assert max(can.spread.values()) < 1e-7


def test_the_stretch_channels_zero_is_the_built_bond_length():
    """``k_bond`` changes the Born tensor and nothing at rest on the lattice path."""
    ch, a = _fitted_packer((("C", "F", 0.1, 0.0),))
    _, b = _fitted_packer((("C", "F", 0.1, 0.6),))
    assert np.abs(a._q_cell - b._q_cell).max() < 1e-14  # r - r0 is the rounding of one bond length
    assert float(a.energy(P0[None])[0]) == pytest.approx(float(b.energy(P0[None])[0]), rel=1e-13)
    Za, Zb = born_charges(a, P0).Z, born_charges(b, P0).Z
    F = [i for i, e in enumerate(b.elements) if e == "F"]
    assert np.abs(Za[F] - Zb[F]).max() > 0.3
    # the difference on a fluorine is -k_bond (r + d) u u^T: rank one, along its own bond, with the
    # transferred charge sitting the off-site distance d beyond the nucleus
    ff = SimpleFF.from_preset("pvdf-dft-valence")
    i = F[0]
    D = Zb[i] - Za[i]
    w, V = np.linalg.eigh(0.5 * (D + D.T))
    assert w[0] == pytest.approx(-0.6 * (1.35 + ff.offset_table()["F"]) * ff.charge_scale, rel=1e-6)
    assert abs(w[1]) < 1e-7 and abs(w[2]) < 1e-7  # the central difference's O(h^2) floor
    carbon = min((j for j, e in enumerate(b.elements) if e == "C"),
                 key=lambda j: np.linalg.norm(born_charges(b, P0).positions[i] - born_charges(b, P0).positions[j]))
    u = born_charges(b, P0).positions[i] - born_charges(b, P0).positions[carbon]
    assert abs(abs(V[:, 0] @ (u / np.linalg.norm(u))) - 1.0) < 1e-6
    # the molecule path keeps the valence fit's r0 and therefore shifts the charges
    from polyfind.forcefield import infer_bonds

    el, X = list(ch.elements), np.asarray(ch.coords, dtype=float)
    bonds = infer_bonds(el, X)
    qa = _flux_ff((("C", "F", 0.1, 0.0),)).charges_at(el, bonds, X)
    qb = _flux_ff((("C", "F", 0.1, 0.6),)).charges_at(el, bonds, X)
    assert np.abs(qa - qb).max() > 0.01


def test_built_length_reference_refuses_a_type_with_two_lengths():
    from polyfind.chain import build_chain
    from polyfind.pack import built_bond_lengths

    tpl = build_chain(PVDF, np.zeros(10), cap=False)
    ff = _flux_ff((("C", "F", 0.0, 0.6),))
    assert built_bond_lengths(tpl, ff) == {"C-F": (0.0, 1.35)}
    assert built_bond_lengths(tpl, _flux_ff((("C", "F", 0.3, 0.0),))) == {}
    import dataclasses

    coords = np.asarray(tpl.coords, dtype=float).copy()
    Fs = [i for i, e in enumerate(tpl.elements) if e == "F"]
    carbon = [j for j in range(len(tpl.elements)) if tpl.elements[j] == "C" and
              np.linalg.norm(coords[j] - coords[Fs[0]]) < 1.5][0]
    coords[Fs[0]] += 0.05 * (coords[Fs[0]] - coords[carbon])  # one C-F bond 5% longer than the rest
    with pytest.raises(ValueError, match="no single stretch reference"):
        built_bond_lengths(dataclasses.replace(tpl, coords=coords), ff)


def test_load_reference_reflects_pendants_onto_the_plus_a_representative(tmp_path):
    """The provider's file: labels from the geometry, the ``+a`` pendant is the representative."""
    el = ["C", "H", "H", "C", "F", "F"]
    X = np.array([[0.0, -0.05, 0.0], [0.9, -0.7, 0.0], [-0.9, -0.7, 0.0],
                  [0.0, 0.78, 1.29], [1.1, 1.63, 1.29], [-1.1, 1.63, 1.29]])
    R = np.diag([-1.0, 1.0, 1.0])
    ZF = np.array([[-1.0, 0.5, 0.0], [0.4, -0.8, 0.0], [0.0, 0.0, -0.4]])
    ZH = np.array([[0.14, 0.03, 0.0], [0.02, 0.15, 0.0], [0.0, 0.0, 0.13]])
    ZCF, ZCH = np.diag([1.8, 1.5, 1.2]), np.diag([-0.1, -0.2, -0.7])
    tensors = [ZCH, ZH, R @ ZH @ R, ZCF, ZF, R @ ZF @ R]
    with open(tmp_path / "g.xyz", "w") as fh:
        fh.write(f"{len(el)}\nLattice\n" + "\n".join(f"{e} {x} {y} {z}" for e, (x, y, z) in zip(el, X)) + "\n")
    d = {"geometry": "g.xyz", "protocol": "synthetic",
         "atoms": [{"index": i + 1, "element": e, "born_effective_charge_e": t.tolist(),
                    "born_effective_charge_asr_corrected_e": (t + 0.01 * np.eye(3)).tolist()}
                   for i, (e, t) in enumerate(zip(el, tensors))],
         "acoustic_sum_tensor_e": np.zeros((3, 3)).tolist(),
         "acoustic_sum_tensor_asr_corrected_e": np.zeros((3, 3)).tolist(),
         "dielectric_tensor_electronic": np.eye(3).tolist()}
    with open(tmp_path / "born_results.json", "w") as fh:
        json.dump(d, fh)
    ref = load_reference(str(tmp_path / "born_results.json"))
    assert set(ref.raw) == {"C(F2)", "F", "C(H2)", "H"} and ref.count["F"] == 2
    assert ref.raw["F"] == pytest.approx(ZF)  # the -a fluorine reflected back onto the +a one
    assert ref.raw["H"] == pytest.approx(ZH)
    assert ref.corrected["C(F2)"] == pytest.approx(ZCF + 0.01 * np.eye(3))
    assert ref.component("F", "ab") == pytest.approx(0.5) and ref.component("F", "ba", corrected=False) == pytest.approx(0.4)
    assert ref.component("F", "aa") == pytest.approx(-0.99) and ref.component("F", "aa", corrected=False) == pytest.approx(-1.0)


def test_labels_name_carbons_by_their_pendants():
    ch = periodic_chain(PVDF, [T, T], THREE_STATE)
    labels = atom_labels(ch)
    assert sorted(labels) == ["C(F2)", "C(H2)", "F", "F", "H", "H"]
