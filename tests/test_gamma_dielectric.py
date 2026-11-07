"""examples/gamma_dielectric.py: the relaxed-ion dielectric tensor by the Gamma route (docs/CRYSTAL_STAGE.md 5.16)."""
import importlib.util
from pathlib import Path

import numpy as np

from polyfind import fitting as F
from polyfind.polymers import PVDF, belt_variant

SCRIPT = Path(__file__).resolve().parents[1] / "examples" / "gamma_dielectric.py"
SPEC = importlib.util.spec_from_file_location("gamma_dielectric", SCRIPT)
GD = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(GD)


def test_belt_sites_follow_their_cf2_group_through_any_rigid_motion():
    """The constraint places each site as a fixed affine combination of its carbon and the two
    fluorines on it, so a finite rotation plus translation of a whole chain carries the sites to
    exactly where the chain builder would put them."""
    from polyfind import phonon as PH
    from polyfind.refine import resolve_angle_stiffness
    polymer = belt_variant(PVDF, geometry="pbe")
    with F.PRESET_PARAMETERS["pvdf-dft-valence-belt2"].applied():
        ref, _ = GD.build(polymer, [0, 0], "beta-belt")
        pk, _, _ = GD.pinned_valence_packer(ref, polymer, resolve_angle_stiffness(None, None))
    P0, _ = PH.placed_coordinates(pk, np.asarray(ref.params, dtype=float))
    P0 = np.asarray(P0, dtype=float)
    J, atoms = GD.site_constraint(pk, P0)
    el = np.asarray(list(pk.elements) * pk.n_chains)
    assert J.shape == (3 * pk.N, 3 * len(atoms)) and len(atoms) == 12 and pk.N == 20
    # rigid motion of chain 0 about its own axis plus a translation
    th = np.deg2rad(7.0)
    R = np.array([[np.cos(th), -np.sin(th), 0.0], [np.sin(th), np.cos(th), 0.0], [0.0, 0.0, 1.0]])
    axis = P0[:pk.n][el[:pk.n] == "C"].mean(axis=0) * np.array([1.0, 1.0, 0.0])
    P1 = P0.copy()
    P1[:pk.n] = (P0[:pk.n] - axis) @ R.T + axis + np.array([0.13, -0.07, 0.21])
    u_atoms = (P1 - P0)[atoms].ravel()
    P_from_J = P0 + (J @ u_atoms).reshape(-1, 3)
    assert np.abs(P_from_J - P1).max() < 1e-9  # the sites included


def test_gamma_route_agrees_with_the_rigid_chain_route_on_the_plain_chain():
    """On the plain chain the rigid-chain route (chains reorienting under the field) and the full
    Gamma route agree on the lateral axis to a few hundredths: the libration carries it in both, the
    extra intra-chain modes add little, and no mode is imaginary at the relaxed reference."""
    from polyfind.refine import resolve_angle_stiffness
    with F.PRESET_PARAMETERS["pvdf-dft-valence"].applied():
        ref, _ = GD.build(PVDF, [0, 0], "beta")
        out = GD.gamma_dielectric(ref, PVDF, resolve_angle_stiffness(None, None))
    rigid, gamma, ion = (np.diag(np.asarray(out[k])) for k in ("eps_0_rigid_chain", "eps_0_gamma", "eps_ion_gamma"))
    assert out["n_negative_modes"] == 0 and out["relaxation"]["max_force_after"] < 1e-6
    assert out["relaxation"]["max_displacement_A"] < 0.05
    assert abs(gamma[1] - rigid[1]) < 0.1 and gamma[1] > rigid[1]  # the full relaxation is at least as soft
    assert ion[1] > 0.5 and ion[0] < 0.15 and ion[2] < 0.3
    top = out["modes_by_contribution"][0]
    assert 40 < top["freq_cm1"] < 90 and top["d_eps_diag"][1] > 0.9 * ion[1]  # the libration carries the lateral axis
    assert np.allclose(np.asarray(out["eps_0_gamma"]), np.asarray(out["eps_0_gamma"]).T, atol=1e-6)
