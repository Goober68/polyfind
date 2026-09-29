import numpy as np
import pytest

from polyfind import film as FI

# A crystal of the shape the deformable path returns: orthorhombic, stiff chain axis z, polar x.
S_CRYSTAL = np.array([[0.04221, -0.00566, -0.00177, 0.0],
                      [-0.00566, 0.04942, -0.00022, 0.0],
                      [-0.00177, -0.00022, 0.00305, 0.0],
                      [0.0, 0.0, 0.0, 0.2307]])
D_CRYSTAL = np.array([[-9.864, -0.587, 0.422, 0.0],
                      [0.0, 0.0, 0.0, 4.1],
                      [0.0, 0.0, 0.0, 0.0]])


def test_all_crystal_is_the_crystal_and_all_amorphous_is_the_amorphous():
    f = FI.laminate(S_CRYSTAL, D_CRYSTAL, 1.0, 1.0, 0.45)
    assert f.S == pytest.approx(S_CRYSTAL, rel=1e-10, abs=1e-14)
    assert f.d[:2] == pytest.approx(D_CRYSTAL[:2], rel=1e-10, abs=1e-12)
    assert np.isnan(f.d[2]).all()  # no field along the layer normal
    g = FI.laminate(S_CRYSTAL, D_CRYSTAL, 0.0, 1.3, 0.3)
    assert g.S == pytest.approx(FI.isotropic_compliance(1.3, 0.3), rel=1e-10)
    assert np.abs(g.d[:2]).max() < 1e-12


def test_identical_layers_homogenise_to_themselves():
    """Splitting a crystal into two laminae of itself changes nothing, at any fraction --
    with a crystal standing in for the amorphous phase via its isotropic twin."""
    S_iso = FI.isotropic_compliance(2.0, 0.35)
    for phi in (0.2, 0.5, 0.8):
        f = FI.laminate(S_iso, np.zeros((3, 4)), phi, 2.0, 0.35)
        assert f.S == pytest.approx(S_iso, rel=1e-12)


def test_equal_poisson_isotropic_laminate_is_voigt_in_plane():
    """Two isotropic layers with the same Poisson ratio: under in-plane stress both are in
    plane stress with a common strain, so the in-plane Young's modulus is exactly the volume
    average and the in-plane Poisson ratio is unchanged."""
    E1, E2, v, phi = 5.0, 0.5, 0.3, 0.4
    f = FI.laminate(FI.isotropic_compliance(E1, v), np.zeros((3, 4)), phi, E2, v)
    assert 1.0 / f.S[0, 0] == pytest.approx(phi * E1 + (1 - phi) * E2, rel=1e-12)
    assert -f.S[0, 1] / f.S[0, 0] == pytest.approx(v, rel=1e-12)
    # and along the normal the layers are in series, with the Poisson correction of plane strain
    assert 1.0 / f.S[2, 2] < phi * E1 + (1 - phi) * E2


def test_film_compliance_is_symmetric_and_positive_definite():
    for phi in (0.3, 0.5, 0.7):
        for v in (0.3, 0.45, 0.49):
            f = FI.laminate(S_CRYSTAL, D_CRYSTAL, phi, 0.5, v)
            assert np.abs(f.S - f.S.T).max() < 1e-14
            assert np.linalg.eigvalsh(f.S).min() > 0.0


def test_an_incompressible_amorphous_layer_turns_lateral_contraction_into_draw_extension():
    """The mechanism the module exists for.  A field along P contracts the crystal in-plane
    (d_x,xx + d_x,yy < 0); the soft amorphous layers are bonded to it in-plane and nearly
    incompressible, so they extend along the draw axis.  The film's d_31 goes positive and
    grows with the amorphous Poisson ratio, far beyond the crystal's own +0.42."""
    fc = {}
    for v in (0.3, 0.45, 0.49):
        fc[v] = FI.laminate(S_CRYSTAL, D_CRYSTAL, 0.5, 0.5, v).film_coefficients(polar=0)
    assert fc[0.3]["d31"] < fc[0.45]["d31"] < fc[0.49]["d31"]
    assert fc[0.49]["d31"] > 5.0 * D_CRYSTAL[0, 2]
    # with a soft amorphous phase the crystal carries the in-plane stiffness, so the film's
    # in-plane coefficients stay near the crystal's
    assert fc[0.49]["d33"] == pytest.approx(D_CRYSTAL[0, 0], rel=0.15)
    # a passive crystal gives a passive film
    z = FI.laminate(S_CRYSTAL, np.zeros((3, 4)), 0.5, 0.5, 0.49)
    assert np.abs(z.d[:2]).max() < 1e-12


def test_inputs_are_checked():
    with pytest.raises(ValueError):
        FI.laminate(S_CRYSTAL, D_CRYSTAL, 1.2, 0.5, 0.4)
    with pytest.raises(ValueError):
        FI.laminate(S_CRYSTAL, D_CRYSTAL, 0.5, 0.5, 0.5)
    with pytest.raises(ValueError):
        FI.laminate(S_CRYSTAL[:3, :3], D_CRYSTAL, 0.5, 0.5, 0.4)
    with pytest.raises(ValueError):
        FI.laminate(S_CRYSTAL, D_CRYSTAL, 0.5, 0.5, 0.4).film_coefficients(polar=2)
