"""Infrared intensities: the contraction is right on a case with a closed form, the Born tensors from the
field derivative are the dipole derivative the model has, and the symmetry bookkeeping enforces what it must.

A diatomic's optical mode carries ``|d|^2 = q^2 / mu`` and its translation nothing; the absolute scale
is rederived from the stated constants; a fixed-charge cell's field-derivative Born tensor is the charge
times the identity; with flux, Ewald and induced dipoles it equals :func:`polyfind.born.born_charges`
at the placed geometry and the dipole-route difference at a displaced one; the packer's own field is put
back; the acoustic modes of a relaxed cell carry no intensity; and in a centrosymmetric cell every
inversion-even mode is dark (mutual exclusion), with the group found from the coordinates alone.
"""
from types import SimpleNamespace

import numpy as np
import pytest

from polyfind import infrared as IR
from polyfind import pack as pack_mod
from polyfind import phonon as PH
from polyfind.born import born_charges
from polyfind.ewald import EwaldSpec
from polyfind.forcefield import SimpleFF
from polyfind.pack import CrystalPacker, periodic_chain
from polyfind.polarizability import Polarizable
from polyfind.polymers import PVDF, THREE_STATE

T, GP, GM = 0, 1, 2
P0 = np.array([4.6, 8.6, 90.0, 0.0, 0.0, 0.4, 0.0])
ALPHA_C2H = np.array([5.0, 9.6, 90.0, 0.0, 0.0, 0.0, 1.0])  # the alpha placement test_phonon checks; C2h


def _fitted(flux=None, **kw):
    """beta-PVDF with ``pvdf-dft-valence``'s charges and valence terms, as tests/test_phonon.py builds it."""
    import polyfind.fitting  # noqa: F401 - registers the presets
    from polyfind.fitting import FITTED_VALENCE

    with FITTED_VALENCE.applied():
        ch = periodic_chain(PVDF, (T, T), THREE_STATE)
        pk = pack_mod.CrystalPacker(ch, n_chains=2, valence=SimpleFF.from_preset("pvdf-dft-valence"),
                                    charge_flux=None if flux is None else SimpleFF.from_preset(flux), **kw)
    return pk


def _fitted_alpha(flux=None, **kw):
    """The same potential on alpha's TGTG' chain."""
    import polyfind.fitting  # noqa: F401 - registers the presets
    from polyfind.fitting import FITTED_VALENCE

    with FITTED_VALENCE.applied():
        ch = periodic_chain(PVDF, (T, GP, T, GM), THREE_STATE)
        pk = pack_mod.CrystalPacker(ch, n_chains=2, valence=SimpleFF.from_preset("pvdf-dft-valence"),
                                    charge_flux=None if flux is None else SimpleFF.from_preset(flux), **kw)
    return pk


def test_absolute_scale_is_derived_from_the_stated_constants():
    """N_A / (12 eps0 c^2) (e^2/amu) in km/mol; the familiar 42.2561 km/mol per (D/A)^2/amu times (4.80320 D/eA)^2."""
    e, NA, eps0, c, amu = 1.602176634e-19, 6.02214076e23, 8.8541878128e-12, 299792458.0, 1.66053906660e-27
    assert IR.KM_PER_MOL_PER_E2_AMU == pytest.approx(NA / (12 * eps0 * c * c) * e * e / amu / 1e3, rel=1e-12)
    assert IR.KM_PER_MOL_PER_E2_AMU == pytest.approx(42.2561 * 4.80320 ** 2, rel=1e-4)
    assert IR.INV_EPS0_V_A == pytest.approx(180.95, abs=0.01)


def test_diatomic_known_answer():
    """Charges +q, -q on masses m1, m2 held by a spring of a different stiffness along each axis: each of
    the three optical modes carries d = q / sqrt(mu) along its own axis, the three translations nothing."""
    m1, m2, q = 12.011, 18.998, 0.37
    ks = (300.0, 50.0, 80.0)
    H = np.zeros((6, 6))
    for a in range(3):
        H[np.ix_([a, 3 + a], [a, 3 + a])] = ks[a] * np.array([[1.0, -1.0], [-1.0, 1.0]])
    m = np.array([m1, m2])
    W = np.repeat(1.0 / np.sqrt(m), 3)
    lam, V = np.linalg.eigh(W[:, None] * H * W[None, :])
    Z = np.array([q * np.eye(3), -q * np.eye(3)])
    d = IR.mode_dipoles(Z, V, m)
    mu = m1 * m2 / (m1 + m2)
    for k in range(3, 6):
        axis = int(np.argmax(np.abs(V[:, k]).reshape(2, 3).sum(axis=0)))
        assert float(d[k] @ d[k]) == pytest.approx(q * q / mu, rel=1e-12)
        assert abs(d[k, axis]) == pytest.approx(q / np.sqrt(mu), rel=1e-12)
        assert lam[k] == pytest.approx(ks[axis] / mu, rel=1e-12)
    assert np.abs(d[:3]).max() < 1e-14  # the translations
    ph = SimpleNamespace(modes=V, masses=m, freq_cm1=np.sign(lam) * np.sqrt(np.abs(lam)) * PH.SQRT_KCAL_A2_AMU_TO_CM1)
    ir = IR.infrared(ph, Z)
    top = ir.modes[5]
    assert top.km_per_mol == pytest.approx(q * q / mu * IR.KM_PER_MOL_PER_E2_AMU, rel=1e-12)
    assert top.shares == pytest.approx([1.0, 0.0, 0.0], abs=1e-12)  # the stiffest spring is along x
    assert np.abs(ir.asr).max() < 1e-15 and top.parity is None and top.allowed is None


def test_fixed_charges_give_the_charge_times_the_identity_and_the_field_is_restored():
    pk = CrystalPacker(periodic_chain(PVDF, [T, T], THREE_STATE), n_chains=2, field=(0.02, 0.0, -0.01))
    Pn, latn = PH.placed_coordinates(pk, P0)
    Z = IR.born_charges_by_field(pk, P0, Pn, latn)
    q = np.asarray(pk._q_cell, dtype=float)
    assert np.abs(Z - q[:, None, None] * np.eye(3)[None]).max() < 1e-9
    assert np.allclose(pk.field, [0.02, 0.0, -0.01]) and pk._field_on
    pk.set_field(None)
    IR.born_charges_by_field(pk, P0, Pn, latn)
    assert pk.field is None and not pk._field_on


def test_field_route_is_the_models_born_tensor_with_flux_ewald_and_dipoles():
    """At the placed geometry ``born.born_charges`` is valid on every tree and the field route must equal it;
    off it, the field route must equal central differences of the cell dipole with per-chain charges."""
    pk = _fitted("pvdf-dft-valence-flux-born", coulomb="ewald", ewald=EwaldSpec(), polarizable=Polarizable())
    Pn, latn = PH.placed_coordinates(pk, P0)
    Z = IR.born_charges_by_field(pk, P0, Pn, latn)
    assert np.abs(Z - born_charges(pk, P0).Z).max() < 1e-6
    assert np.abs(Z.sum(axis=0)).max() < 1e-9
    q = np.asarray(pk._q_cell, dtype=float)
    assert np.abs(Z - q[:, None, None] * np.eye(3)[None]).max() > 0.05  # the flux and the dipoles act
    # displaced: the two chains now carry different fluxed charges and different induced dipoles
    rng = np.random.default_rng(5)
    Pd = Pn + 0.02 * rng.standard_normal(Pn.shape)
    n, cz = pk.n, float(latn[2, 2])

    def dipole(P):
        if hasattr(pk, "placed_dipole"):
            return np.sum(pk.placed_dipole(P, latn, float(P0[6])), axis=0)
        from polyfind.born import chain_frame_coords
        q = np.concatenate([pk._flux.charges(chain_frame_coords(P0, latn, P[s * n:(s + 1) * n], s), cz)[0] for s in range(2)])
        return q @ P + pk._polarize(P, q, latn, float(P0[6]))[1].sum(axis=0)

    Zf = IR.born_charges_by_field(pk, P0, Pd, latn)
    Zd = IR.born_charges_by_dipole(dipole, Pd)
    assert np.abs(Zf - Zd).max() < 1e-6
    # the clamped-ion dielectric tensor is the induced dipoles' and nothing else
    eps = IR.electronic_dielectric(pk, P0, Pn, latn)
    assert np.allclose(eps, eps.T, atol=1e-8) and np.all(np.diag(eps) > 1.5) and np.all(np.diag(eps) < 3.5)
    assert np.allclose(IR.electronic_dielectric(_fitted(), P0, Pn, latn), np.eye(3), atol=1e-8)


def test_acoustic_modes_of_a_relaxed_cell_are_dark_and_the_polar_cell_is_c2v():
    """Beta-PVDF (polar, no inversion): the translations carry no intensity because the sum rule holds,
    the group from the coordinates is C2v with its two-fold along the polar axis, and every mode's dipole
    lies along the axes its characters allow."""
    pk = _fitted("pvdf-dft-valence-flux-born")
    rel = PH.relax_all_atom(pk, P0)
    assert rel.converged
    latn = PH.placed_coordinates(pk, P0)[1]
    ph = PH.phonons_gamma(pk, P0, Pn=rel.Pn, latn=latn, asr="project")
    Z = IR.born_charges_by_field(pk, P0, rel.Pn, latn)
    ops = IR.space_group(rel.Pn, latn, list(pk.elements) * 2)
    assert IR.point_group_name(ops) == "C2v"
    two = [o for o in ops if o.kind == "2-fold"]
    assert two and np.allclose(np.abs(two[0].axis), [1.0, 0.0, 0.0])
    ir = IR.infrared(ph, Z, ops)
    top = max(m.intensity for m in ir.modes)
    assert max(ir.modes[i].intensity for i in ph.acoustic) < 1e-20 * top
    assert ir.born_symmetry < 1e-6 and not ir.centrosymmetric
    assert ir.forbidden_leak(ph.acoustic) < 1e-12
    optical = [k for k in range(ir.chi.shape[1]) if k not in ph.acoustic]  # the translations are degenerate
    assert np.abs(np.abs(ir.chi[:, optical]) - 1.0).max() < 1e-4
    assert any(m.allowed == "" for m in ir.modes) and any(m.intensity > 1e-3 * top for m in ir.modes if m.allowed)


def test_mutual_exclusion_in_a_centrosymmetric_alpha_cell():
    """The alpha placement at ALPHA_C2H is C2h with an inversion.  The Hessian there need not be at a
    minimum -- it is a function of the coordinates alone, so it commutes with the group -- and with the
    fluxing charges' Born tensors every inversion-even mode is dark, Au along the two-fold axis, Bu across."""
    pk = _fitted_alpha("pvdf-dft-valence-flux-born")
    Pn, latn = PH.placed_coordinates(pk, ALPHA_C2H)
    ops = IR.space_group(Pn, latn, list(pk.elements) * 2)
    assert IR.point_group_name(ops) == "C2h"
    inv = [o for o in ops if o.kind == "inversion"][0]
    assert inv.error < 1e-9
    ph = PH.phonons_gamma(pk, ALPHA_C2H, h=1e-3, asr="project")
    Z = IR.born_charges_by_field(pk, ALPHA_C2H, Pn, latn)
    ir = IR.infrared(ph, Z, ops)
    assert ir.centrosymmetric and ir.born_symmetry < 1e-9
    clean = [m for m in ir.modes if abs(abs(m.parity) - 1.0) < 1e-4 and m.index not in ph.acoustic]
    assert len(clean) > 0.9 * (len(ir.modes) - 3)
    top = max(m.intensity for m in ir.modes)
    assert max(m.intensity for m in clean if m.parity > 0) < 1e-12 * top
    assert max(m.intensity for m in clean if m.parity < 0) > 1e-3 * top
    labels = {m.irrep for m in clean}
    assert labels <= {"Ag", "Bg", "Au", "Bu"} and {"Au", "Bu"} <= labels
    # Au along the two-fold axis, Bu perpendicular to it
    axis = np.abs([o for o in ops if o.kind == "2-fold"][0].axis)
    for m in clean:
        if m.irrep in ("Au", "Bu") and m.intensity > 1e-6 * top:
            along = float((m.dipole @ axis) ** 2 / m.intensity)
            assert along == pytest.approx(1.0 if m.irrep == "Au" else 0.0, abs=1e-6)


def test_space_group_rejects_a_broken_inversion():
    alpha = periodic_chain(PVDF, [T, GP, T, GM], THREE_STATE)
    pk = CrystalPacker(alpha, n_chains=2)
    Pn, latn = PH.placed_coordinates(pk, ALPHA_C2H)
    Pb = Pn.copy()
    Pb[4, 0] += 0.01  # one fluorine, 0.01 A
    assert pk.elements[4] == "F"
    assert IR.point_group_name(IR.space_group(Pb, latn, list(pk.elements) * 2, tol=1e-3)) == "C1"
    assert IR.point_group_name(IR.space_group(Pb, latn, list(pk.elements) * 2, tol=2e-2)) == "C2h"


def test_lo_shift_is_zero_for_a_dark_mode_and_grows_with_the_dipole():
    eps = np.diag([2.0, 2.0, 2.5])
    assert IR.lo_shift(60.0, np.zeros(3), eps, 200.0) == 60.0
    lo1 = IR.lo_shift(60.0, np.array([0.05, 0.0, 0.0]), eps, 200.0)
    lo2 = IR.lo_shift(60.0, np.array([0.10, 0.0, 0.0]), eps, 200.0)
    assert 60.0 < lo1 < lo2
    lam = (60.0 / PH.SQRT_KCAL_A2_AMU_TO_CM1) ** 2 + 4 * np.pi * 332.0637 * 0.01 / (200.0 * 2.0)
    assert lo2 == pytest.approx(np.sqrt(lam) * PH.SQRT_KCAL_A2_AMU_TO_CM1, rel=1e-12)
