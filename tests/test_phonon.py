"""Gamma-point lattice dynamics: the assembled all-atom energy is the packer's, its gradient is exact,
and the Hessian built from it is symmetric, translation-invariant and, at a minimum, positive.

The known-answer identity ``cell_energy == packer.energy`` is checked on a plain LJ/DSF packer,
an Ewald one, a flipped cell, and the fitted flux + polarizable + valence beta-PVDF packer.
"""
import numpy as np
import pytest

from polyfind import pack as pack_mod
from polyfind import phonon as PH
from polyfind.ewald import EwaldSpec
from polyfind.forcefield import SimpleFF
from polyfind.pack import CrystalPacker, periodic_chain
from polyfind.polarizability import Polarizable
from polyfind.polymers import PVDF, THREE_STATE

T, GP, GM = 0, 1, 2
P0 = np.array([4.6, 8.6, 90.0, 0.0, 0.0, 0.4, 0.0])
P1 = np.array([4.7, 8.5, 88.0, 25.0, -40.0, 0.7, 1.0])  # setting angles, flip, oblique cell


def _fitted(flux=None, **kw):
    """beta-PVDF with ``pvdf-dft-valence``'s charges and valence terms; see tests/test_born.py."""
    import polyfind.fitting  # noqa: F401 - registers the presets
    from polyfind.fitting import FITTED_VALENCE

    with FITTED_VALENCE.applied():
        ch = periodic_chain(PVDF, (T, T), THREE_STATE)
        pk = pack_mod.CrystalPacker(ch, n_chains=2, valence=SimpleFF.from_preset("pvdf-dft-valence"),
                                    charge_flux=None if flux is None else SimpleFF.from_preset(flux), **kw)
    return pk


@pytest.fixture(scope="module")
def beta_chain():
    return periodic_chain(PVDF, [T, T], THREE_STATE)


@pytest.fixture(scope="module")
def full_packer():
    """Flux + induced dipoles + Ewald + valence: every term the module has to assemble."""
    return _fitted("pvdf-dft-valence-flux-born", coulomb="ewald", ewald=EwaldSpec(), polarizable=Polarizable())


def _fd_all_atom(pk, params, Pn, latn, gP, rng, n=6, h=1e-5):
    worst = 0.0
    for _ in range(n):
        i, d = int(rng.integers(pk.N)), int(rng.integers(3))
        Pp, Pm = Pn.copy(), Pn.copy()
        Pp[i, d] += h
        Pm[i, d] -= h
        fd = (PH.cell_energy(pk, params, Pp, latn) - PH.cell_energy(pk, params, Pm, latn)) / (2 * h)
        worst = max(worst, abs(fd - gP[i, d]))
    return worst


def _fd_symmetric(pk, params, gX, rng, n=4, h=1e-5):
    """A repeat-coordinate displacement, placed into every chain, against ``g_coords``."""
    X0 = np.asarray(pk.chain.coords, dtype=float)
    worst = 0.0
    for _ in range(n):
        i, d = int(rng.integers(pk.n)), int(rng.integers(3))
        es = []
        for sgn in (1.0, -1.0):
            X = X0.copy()
            X[i, d] += sgn * h
            P, lat = pk._place(params[None], coords=X[None], c=np.array([pk.chain.c]))
            es.append(PH.cell_energy(pk, params, np.asarray(P)[0], np.asarray(lat)[0]))
        worst = max(worst, abs((es[0] - es[1]) / (2 * h) - gX[i, d]))
    return worst


@pytest.mark.parametrize("params", [P0, P1])
@pytest.mark.parametrize("kw", [{}, dict(coulomb="ewald", ewald=EwaldSpec()), dict(lj_cutoff="force", field=(0.1, -0.05, 0.02))])
def test_known_answer_and_gradient_on_plain_packers(beta_chain, params, kw):
    """``cell_energy`` is ``packer.energy`` to rounding; its gradient is exact and folds to ``g_coords``."""
    pk = CrystalPacker(beta_chain, n_chains=2, **kw)
    e_cell, e_pack, rel = PH.check_known_answer(pk, params)
    assert rel < 1e-12
    Pn, latn = PH.placed_coordinates(pk, params)
    E, gP = PH.cell_energy_and_grad(pk, params, Pn, latn)
    assert E == pytest.approx(e_pack, rel=1e-12)
    rng = np.random.default_rng(1)
    assert _fd_all_atom(pk, params, Pn, latn, gP, rng) < 1e-7
    _, _, gX, _ = pk.energy_and_grad(params)
    assert np.abs(PH.fold_gradient(pk, params, gP) - gX).max() < 1e-10
    assert _fd_symmetric(pk, params, gX, rng) < 1e-6


def test_known_answer_one_chain_and_alpha(beta_chain):
    pk = CrystalPacker(beta_chain, n_chains=1)
    assert PH.check_known_answer(pk, np.array([4.6, 8.6, 90.0, 20.0, 0.0, 0.0, 0.0]))[2] < 1e-12
    alpha = periodic_chain(PVDF, [T, GP, T, GM], THREE_STATE)
    pk = CrystalPacker(alpha, n_chains=2)
    p = np.array([5.0, 9.6, 90.0, 0.0, 0.0, 0.0, 1.0])
    assert PH.check_known_answer(pk, p)[2] < 1e-12
    Pn, latn = PH.placed_coordinates(pk, p)
    _, gP = PH.cell_energy_and_grad(pk, p, Pn, latn)
    _, _, gX, _ = pk.energy_and_grad(p)
    assert np.abs(PH.fold_gradient(pk, p, gP) - gX).max() < 1e-10


@pytest.mark.parametrize("params", [P0, P1])
def test_known_answer_and_gradient_with_flux_dipoles_ewald_and_valence(full_packer, params):
    """The mandatory check on beta-PVDF with every term on, and the symmetric-displacement check."""
    pk = full_packer
    assert pk._flux is not None and pk.polarizable is not None and pk._valence is not None
    e_cell, e_pack, rel = PH.check_known_answer(pk, params)
    assert rel < 1e-12
    Pn, latn = PH.placed_coordinates(pk, params)
    E, gP = PH.cell_energy_and_grad(pk, params, Pn, latn)
    rng = np.random.default_rng(2)
    assert _fd_all_atom(pk, params, Pn, latn, gP, rng) < 1e-6
    _, _, gX, _ = pk.energy_and_grad(params)
    assert np.abs(PH.fold_gradient(pk, params, gP) - gX).max() < 1e-9
    assert _fd_symmetric(pk, params, gX, rng) < 1e-6


def test_gradient_is_exact_on_an_asymmetrically_displaced_cell(full_packer):
    """Off the symmetric geometry the two chains carry different fluxed charges and different induced
    dipoles; the gradient re-solves both per chain, and a central difference confirms it."""
    pk = full_packer
    Pn, latn = PH.placed_coordinates(pk, P1)
    rng = np.random.default_rng(3)
    Pn = Pn + 0.03 * rng.standard_normal(Pn.shape)
    _, gP = PH.cell_energy_and_grad(pk, P1, Pn, latn)
    assert _fd_all_atom(pk, P1, Pn, latn, gP, rng, n=8) < 1e-6
    # the two chains' charges are now their own, not one repeat's tiled
    n, cz = pk.n, float(latn[2, 2])
    from polyfind.born import chain_frame_coords
    q1 = pk._flux.charges(chain_frame_coords(P1, latn, Pn[:n], 0), cz)[0]
    q2 = pk._flux.charges(chain_frame_coords(P1, latn, Pn[n:], 1), cz)[0]
    assert np.abs(q1 - q2).max() > 1e-4


def test_built_bond_length_packer_is_consistent_and_leaves_the_original_alone():
    pk = _fitted()
    e0 = float(pk.energy(P0[None])[0])
    pk3 = PH.packer_with_built_bond_lengths(pk)
    assert float(pk.energy(P0[None])[0]) == e0 and pk._valence is not pk3._valence
    assert PH.check_known_answer(pk3, P0)[2] < 1e-12
    ft0, ft3 = PH.force_terms(pk, P0, *PH.placed_coordinates(pk, P0)), PH.force_terms(pk3, P0, *PH.placed_coordinates(pk3, P0))
    assert ft0["bond_strain"] > 0.05 and ft3["bond_strain"] < 1e-12 and ft3["bonds"] < 1e-9
    assert ft3["angles"] == pytest.approx(ft0["angles"])


def test_hessian_is_symmetric_translation_invariant_and_matches_energy_differences():
    pk = _fitted()  # valence, DSF: the cheap packer with internal stiffness
    hs = PH.gamma_hessian(pk, P0)
    H, Pn, latn = hs  # unpacks as the three-tuple the design asks for
    assert H.shape == (3 * pk.N, 3 * pk.N) and np.allclose(H, H.T)
    assert hs.asymmetry < 1e-5
    scale = np.abs(H).max()
    asr = H.reshape(pk.N, 3, pk.N, 3).sum(axis=2)
    assert np.abs(asr).max() < 1e-6 * scale
    pairs = [(0, 0), (3, 4), (5, 20), (17, 17), (2, 35)]
    He = PH.hessian_by_energy(pk, P0, pairs=pairs)
    for k, l in pairs:
        assert H[k, l] == pytest.approx(He[k, l], rel=2e-6, abs=1e-6)
    ph = PH.phonons_from_hessian(pk, hs)
    assert ph.n_zero == 3 and len(ph.acoustic) == 3
    assert all(ph.mode(i).translation > 0.99 for i in ph.acoustic)
    assert ph.freq_cm1.shape == (3 * pk.N,) and np.all(np.diff(ph.freq_cm1) >= 0)
    assert ph.freq_thz == pytest.approx(ph.freq_cm1 * (PH.SQRT_KCAL_A2_AMU_TO_THZ / PH.SQRT_KCAL_A2_AMU_TO_CM1))


def test_relaxed_cell_has_three_zero_modes_and_no_imaginary_ones():
    """At a minimum over every atom the Hessian is positive semi-definite with exactly the three
    translations null; projecting the acoustic sum rule then moves no optical mode."""
    pk = _fitted()
    Pn, E, gmax = PH.relax_all_atom(pk, P0)
    assert gmax < 1e-4
    latn = PH.placed_coordinates(pk, P0)[1]
    raw = PH.phonons_gamma(pk, P0, Pn=Pn, latn=latn)
    assert raw.n_imaginary == 0 and raw.n_zero == 3
    assert raw.optical.min() > 5.0
    proj = PH.phonons_gamma(pk, P0, Pn=Pn, latn=latn, asr="project")
    assert proj.raw is not None and np.abs(proj.freq_cm1[proj.acoustic]).max() < 1e-3  # eigh's rounding on a ~1e3 matrix, in cm^-1
    assert np.abs(proj.optical - raw.optical).max() < 1e-2
    m = proj.lowest_optical[0]
    assert m.transverse + m.axial == pytest.approx(1.0, abs=1e-9)
    assert sum(m.by_type.values()) == pytest.approx(1.0, abs=1e-9)
    assert set(m.by_type) == {"C(F2)", "F", "C(H2)", "H"}
    assert "cm^-1" in proj.report()


def test_unit_conversion_is_derived_from_the_stated_constants():
    kcal_J = 4184.0 / 6.02214076e23
    omega2 = kcal_J / (1.66053906660e-27 * 1e-20)
    assert PH.SQRT_KCAL_A2_AMU_TO_CM1 == pytest.approx(np.sqrt(omega2) / (2 * np.pi * 2.99792458e10), rel=1e-12)
    assert PH.SQRT_KCAL_A2_AMU_TO_CM1 == pytest.approx(108.59, abs=0.01)
    assert PH.SQRT_KCAL_A2_AMU_TO_THZ == pytest.approx(3.2556, abs=1e-3)


def test_asr_must_be_none_or_project():
    pk = CrystalPacker(periodic_chain(PVDF, [T, T], THREE_STATE), n_chains=2)
    with pytest.raises(ValueError):
        PH.phonons_gamma(pk, P0, asr="fix")
