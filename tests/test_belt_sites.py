"""Fluorine belt sites (polymers.belt_fluorine, pseudo-element "X"): placed perpendicular to the C-F bond in
the CF2 plane, bonded to nothing, charged with nothing, inert at zero well depth, and a wall on the X-H pair
that acts only between chains."""
from dataclasses import replace

import numpy as np

from polyfind import pack as P
from polyfind import fitting as F
from polyfind.forcefield import infer_bonds, offset_sites
from polyfind.polymers import PVDF, THREE_STATE, pvdf_belt
from polyfind.ewald import EwaldSpec
from polyfind.polarizability import Polarizable

KW = {"lj_cutoff": "force", "coulomb": "ewald", "ewald": EwaldSpec(), "polarizable": Polarizable()}


def _chains(d=0.7):
    return P.periodic_chain(PVDF, [0, 0], THREE_STATE), P.periodic_chain(pvdf_belt(d), [0, 0], THREE_STATE)


def test_sites_sit_in_the_belt_and_are_bonded_to_nothing():
    ch0, ch1 = _chains(0.7)
    el = list(ch1.elements)
    assert ch1.n_atoms == ch0.n_atoms + 4 and el.count("X") == 4
    X = np.asarray(ch1.coords)
    iC = [i for i, e in enumerate(el) if e == "C"]; iF = [i for i, e in enumerate(el) if e == "F"]
    for x in (i for i, e in enumerate(el) if e == "X"):
        f = min(iF, key=lambda j: np.linalg.norm(X[x] - X[j]))
        c = min(iC, key=lambda j: np.linalg.norm(X[f] - X[j]))
        u, v = X[f] - X[c], X[x] - X[f]
        assert np.isclose(np.linalg.norm(v), 0.7, atol=1e-9)
        assert abs(u @ v) < 1e-9  # perpendicular to the C-F bond
        assert abs(v[2]) < 1e-9  # in the CF2 plane (perpendicular to the chain axis)
    bonds = infer_bonds(el, X)
    assert all("X" not in (el[a], el[b]) for a, b in bonds)
    # the fluorine stays terminal for the off-site charge even with sites beside it
    sites = offset_sites(el, bonds, {"F": 0.18})
    assert len(sites) == 2 and all(el[a] == "F" for a, _, _ in sites)


def test_sites_carry_no_charge_and_no_lattice_interaction_by_default():
    ch0, ch1 = _chains()
    with F.FITTED_VALENCE.applied():
        pk0 = P.CrystalPacker(ch0, n_chains=2, **KW)
        pk1 = P.CrystalPacker(ch1, n_chains=2, **KW)
        n = ch1.n_atoms
        el = np.asarray(list(ch1.elements)); xm = (el[:, None] == "X") | (el[None, :] == "X")
        assert np.all(pk1._q_cell[:n][el == "X"] == 0.0)
        for name in ("_A_nn", "_B_nn", "_const_nn", "_qqf_nn"):
            assert np.all(np.asarray(getattr(pk1, name))[xm] == 0.0), name
        assert np.all(pk1._scale_column_np(4)[:, xm] == 0.0)  # nothing within its own chain
        assert pk1.e_intra == pk0.e_intra and np.isclose(pk1.e_excl, pk0.e_excl, rtol=1e-12, atol=0)  # the zero charges only reorder a sum
        p = np.array([4.52, 8.52, 90.0, 0.0, 0.0, 0.0, 0.0])
        e0, e1 = float(pk0.energy(p[None])[0]), float(pk1.energy(p[None])[0])
        assert abs(e1 - e0) < 1e-4, (e0, e1)  # the massless sites move the chain frame by a hair, nothing else
    ff = F.FITTED_VALENCE.simple_ff()
    iu, ju, lj_x, lj_d, qq = ff._pair_terms(list(ch1.elements), np.full((n, n), 9), np.zeros(n))
    assert not np.any((el[iu] == "X") | (el[ju] == "X"))  # the force field's own pairs skip the sites


def test_a_wall_on_the_x_h_pair_acts_between_chains_only():
    ch0, ch1 = _chains()
    base = F.FITTED_VALENCE
    wall = replace(base, lj=tuple(base.lj) + (("X", 2.6, 0.05),),
                   lj_pairs=(("C", "X", 2.6, 0.0), ("F", "X", 2.6, 0.0), ("X", "X", 2.6, 0.0)))
    with wall.applied():
        pk = P.CrystalPacker(ch1, n_chains=2, **KW)
        el = np.asarray(list(ch1.elements))
        hx = (el[:, None] == "X") & (el[None, :] == "H") | (el[:, None] == "H") & (el[None, :] == "X")
        A = np.asarray(pk._A_nn); B = np.asarray(pk._B_nn)
        assert np.all(A[hx] > 0.0) and np.all(B[hx] == 0.0)  # repulsion only
        cx = ((el[:, None] == "X") & (el[None, :] != "H") & (el[None, :] != "X"))
        assert np.all(A[cx] == 0.0)  # and only against hydrogen
        assert pk.e_intra == P.CrystalPacker(ch0, n_chains=2, **KW).e_intra  # its own chain unchanged
        # the wall raises the energy of a close-packed cell and vanishes for a wide one
        tight = np.array([4.3, 8.3, 90.0, 0.0, 0.0, 0.0, 0.0]); wide = np.array([7.0, 12.0, 90.0, 0.0, 0.0, 0.0, 0.0])
    with base.applied():
        pk_plain = P.CrystalPacker(ch0, n_chains=2, **KW)
    e_tight = float(pk.energy(tight[None])[0]) - float(pk_plain.energy(tight[None])[0])
    e_wide = float(pk.energy(wide[None])[0]) - float(pk_plain.energy(wide[None])[0])
    assert e_tight > 1.0 and abs(e_wide) < 1e-2, (e_tight, e_wide)  # the force-shifted tail is all that is left at 7 x 12 A


def test_inert_sites_are_left_out_of_the_electrostatics():
    """Belt sites have no charge and no polarizability, so the Ewald sums and the dipole solve run
    on the atoms alone (``CrystalPacker._es_idx``): the sites' induced dipoles are zero, the charge
    sum over every site is the sum over the active ones (a zero charge adds nothing), and the
    analytic gradient still matches finite differences."""
    _, ch1 = _chains()
    with F.FITTED_VALENCE_BELT.applied():
        pk = P.CrystalPacker(ch1, n_chains=2, **KW)
        el = np.tile(np.asarray(list(ch1.elements)), 2)
        assert pk._es_idx is not None and np.array_equal(pk._es_idx, np.flatnonzero(el != "X"))
        assert pk._alpha_cell[el == "X"].max() == 0.0 and pk._alpha_es.min() > 0.0
        params = np.array([4.73, 8.52, 90.0, 0.0, 0.0, 0.0, 0.0])
        P_, lat = pk._place(params[None])
        Pn, latn = np.asarray(P_[0], dtype=float), np.asarray(lat[0], dtype=float)
        e_es, p, E0, _ = pk._polarize(Pn, pk._q_cell, latn, 0.0)
        assert np.all(p[el == "X"] == 0.0) and np.abs(p[el != "X"]).max() > 1e-3
        full = pk._ewald.terms(Pn, pk._q_cell, latn).total
        act = pk._ewald.terms(Pn[pk._es_idx], pk._q_cell[pk._es_idx], latn).total
        assert abs(full - act) < 1e-9
        E, g_cell, _, _ = pk.energy_and_grad(params)
        assert float(pk.energy(params[None])[0]) == E
        for i, h in ((0, 1e-4), (1, 1e-4), (3, 1e-3)):
            pp, pm = params.copy(), params.copy()
            pp[i] += h
            pm[i] -= h
            fd = (float(pk.energy(pp[None])[0]) - float(pk.energy(pm[None])[0])) / (2 * h)
            assert abs(fd - g_cell[i]) < 2e-5 * max(1.0, abs(fd)), (i, fd, g_cell[i])
