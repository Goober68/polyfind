"""The torsion gate at a chain end and the registry corrugation of a crystal."""
import numpy as np
import pytest

import polyfind.fitting  # noqa: F401  (registers the fitted SimpleFF presets)
from polyfind import chainend as CE
from polyfind import mechanics as M
from polyfind.forcefield import SimpleFF
from polyfind.polymers import PE, PVDF, TRFE

T = 0


@pytest.fixture(scope="module")
def fitted():
    return SimpleFF.from_preset("pvdf-dft-fit")


def test_vdf_end_gate_is_writable_below_two_volts_per_nanometre(fitted):
    g = CE.torsion_gate(PVDF, fitted)
    assert g.bond == 15 and g.gauche in (60.0, -60.0)
    assert np.isfinite(g.delta_e)
    assert 0.3 < g.dmu_perp < 1.0
    assert abs(np.linalg.norm(g.field_axis) - 1.0) < 1e-12
    # trans is favoured at the end at zero field, so holding trans costs less than writing gauche
    assert g.delta_e > 0.0
    assert 0.0 <= g.field_to_trans < g.field_to_gauche < 0.3  # V/A; the design's limit is 0.2 (2 V/nm)
    # the margin asked is 2 kT at 298 K
    assert g.margin_kcal == pytest.approx(2 * 8.617333262e-5 * 298.15 * 23.060548)
    # the field scales with the margin: a one-kT bias needs less
    g1 = CE.torsion_gate(PVDF, fitted, margin_kt=1.0)
    assert g1.field_to_gauche < g.field_to_gauche


def test_an_interior_bond_and_a_chiral_unit_are_handled(fitted):
    inner = CE.torsion_gate(PVDF, fitted, bond=14)
    assert inner.bond == 14 and np.isfinite(inner.delta_e) and inner.dmu_perp > 0.3
    t = CE.torsion_gate(TRFE, fitted)
    assert np.isfinite(t.delta_e) and t.dmu_perp > 0.3
    assert (t.field_to_gauche == 0.0) != (t.field_to_trans == 0.0) or (t.field_to_gauche > 0 and t.field_to_trans > 0)


def test_polyethylene_has_no_dipole_to_gate(fitted):
    g = CE.torsion_gate(PE, fitted)
    assert g.dmu_perp < 1e-9  # the CH2 dipoles cancel in every conformation
    assert np.isinf(g.field_to_gauche) and np.isinf(g.field_to_trans)


def test_registry_corrugation_of_beta_is_periodic_symmetric_and_minimal_at_the_refined_offset():
    ref, _ = M.refined_reference(PVDF, [T, T], label="beta", lj_cutoff="force", pack_kw={"table_cache_dir": None})
    ref = M.relax_reference(ref)
    cor = CE.registry_corrugation(ref, monomers_per_cell=2, n=41)
    assert cor.amplitude > 0.1
    assert cor.energy[0] == pytest.approx(0.0, abs=1e-6)  # the refined dz is the minimum
    assert cor.energy[-1] == pytest.approx(0.0, abs=1e-6)  # one full repeat later
    assert np.allclose(cor.energy, cor.energy[::-1], atol=1e-6)  # E(dz) = E(c - dz)
    assert cor.half_slip == pytest.approx(cor.amplitude, rel=1e-6)  # the half-slipped cell is the top
    assert cor.peak_force_pn == pytest.approx(cor.peak_slope * 69.478, rel=1e-3)
    assert cor.c == pytest.approx(ref.c)
