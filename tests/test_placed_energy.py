"""Complete placed-cell forces and lattice derivatives, not material validation."""
import numpy as np
import pytest

from polyfind.ewald import EwaldSpec
from polyfind.forcefield import SimpleFF
from polyfind.mechanics import strain_tensor
from polyfind.pack import CrystalPacker, periodic_chain
from polyfind.periodic_geometry import PlacedCell
from polyfind.polarizability import Polarizable
from polyfind.polymers import PVDF, THREE_STATE


def packer(mode, sequence=(0,0)):
    if mode == "full":
        from polyfind.fitting import FITTED_VALENCE
        from polyfind import pack as pack_mod
        with FITTED_VALENCE.applied():
            chain = periodic_chain(PVDF,sequence,THREE_STATE)
            return pack_mod.CrystalPacker(chain,valence=SimpleFF.from_preset("pvdf-dft-valence"),
                                 charge_flux=SimpleFF.from_preset("pvdf-dft-valence-flux-born"),
                                 coulomb="ewald",ewald=EwaldSpec(),polarizable=Polarizable(),
                                 lj_cutoff="force",field=(.02,-.03,.01))
    return CrystalPacker(periodic_chain(PVDF,sequence,THREE_STATE),lj_cutoff="force",
                         **({"coulomb":"ewald","ewald":EwaldSpec()} if mode == "ewald" else {}))


def fd(value, evaluate, h=1e-6):
    result = np.zeros_like(value)
    for idx in np.ndindex(value.shape):
        plus,minus = value.copy(),value.copy()
        plus[idx] += h
        minus[idx] -= h
        result[idx] = (evaluate(plus)-evaluate(minus))/(2*h)
    return result


@pytest.mark.parametrize("mode",["plain","ewald","full"])
@pytest.mark.parametrize("flip",[0.,1.])
def test_full_placed_energy_coordinate_lattice_and_strain_derivatives(mode,flip):
    pk = packer(mode)
    params = np.array([4.7,8.5,88.,25.,-40.,.7,flip])
    P,H = (v[0] for v in pk._place(params[None]))
    baseline = pk.placed_energy_and_grad(P,H,flip)[0]
    assert baseline == pytest.approx(pk.energy(params[None])[0],abs=3e-11)
    F = np.array([[1.01,.006,-.009],[.006,.99,.011],[-.009,.011,1.02]])
    P = P@F+np.random.default_rng(4).normal(size=P.shape)*.003
    H = H@F
    E,gP,gH = pk.placed_energy_and_grad(P,H,flip)
    np.testing.assert_allclose(fd(P,lambda X:pk.placed_energy_and_grad(X,H,flip)[0]),gP,atol=2e-5)
    np.testing.assert_allclose(fd(H,lambda L:pk.placed_energy_and_grad(P,L,flip)[0]),gH,atol=2e-5)
    analytical = [np.sum(gP*(P@strain_tensor(k)))+np.sum(gH*(H@strain_tensor(k))) for k in np.eye(6)]
    for h in (1e-4,5e-5):
        numerical = []
        for k in np.eye(6):
            energies = []
            for sign in (-1,1):
                deformation = np.eye(3)+sign*h*strain_tensor(k)
                energies.append(pk.placed_energy_and_grad(P@deformation,H@deformation,flip)[0])
            numerical.append((energies[1]-energies[0])/(2*h))
        np.testing.assert_allclose(numerical,analytical,atol=2e-4)


def test_nonprimitive_cartesian_torsion_is_in_complete_energy():
    pk = packer("plain",(0,0,0,0))
    params = np.array([4.7,8.5,88.,25.,-40.,.7,1.])
    P,H = (v[0] for v in pk._place(params[None]))
    P = P.copy()
    P[pk.chain.backbone[0],1] += .07
    E,gP,gH = pk.placed_energy_and_grad(P,H,1.)
    np.testing.assert_allclose(fd(P,lambda X:pk.placed_energy_and_grad(X,H,1.)[0]),gP,atol=2e-6)
    np.testing.assert_allclose(fd(H,lambda L:pk.placed_energy_and_grad(P,L,1.)[0]),gH,atol=2e-6)
    Et,gt,_ = pk._torsion.energy_and_repeat_grad(P[:pk.n],H[2])
    assert Et > .01 and np.linalg.norm(gt) > .1


@pytest.mark.parametrize("mode",["plain","ewald","full"])
def test_placed_energy_rotation_translation_and_torque(mode):
    pk = packer(mode,(0,1,0,2))
    params = np.array([5.,9.6,88.,25.,-40.,.7,1.])
    P,H = (v[0] for v in pk._place(params[None]))
    E,gP,gH = pk.placed_energy_and_grad(P,H,1.)
    field = np.zeros(3) if pk.field is None else np.array(pk.field)
    Q,_ = np.linalg.qr(np.random.default_rng(41).normal(size=(3,3)))
    pk.set_field(field@Q)
    Er,gPr,gHr = pk.placed_energy_and_grad(P@Q,H@Q,1.)
    assert Er == pytest.approx(E,abs=3e-11)
    np.testing.assert_allclose(gPr,gP@Q,atol=3e-11)
    np.testing.assert_allclose(gHr,gH@Q,atol=3e-11)
    pk.set_field(field)
    Et,gPt,gHt = pk.placed_energy_and_grad(P+[.31,-.27,.44],H,1.)
    assert Et == pytest.approx(E,abs=3e-11)
    np.testing.assert_allclose(gPt,gP,atol=3e-11)
    np.testing.assert_allclose(gHt,gH,atol=3e-11)
    np.testing.assert_allclose(gP.sum(axis=0),0.,atol=3e-11)
    if not pk._field_on:
        np.testing.assert_allclose(np.cross(P,gP).sum(axis=0)+np.cross(H,gH).sum(axis=0),0.,atol=3e-11)


def test_triclinic_pair_images_cover_bruteforce_cutoff_and_wrapping():
    P = np.random.default_rng(6).normal(size=(5,3))
    H = np.array([[2.5,0.,0.],[1.2,2.8,0.],[.8,.7,2.1]])
    cutoff = 2.2
    def observed(X):
        result = []
        for D,integers in PlacedCell(X,H).pair_image_chunks(cutoff,31):
            np.testing.assert_allclose(D,X[None,:,None,:]-X[None,None,:,:]+integers@H,atol=3e-15)
            result.extend(np.round(D[(D*D).sum(-1)<cutoff**2],10).tolist())
        return sorted(map(tuple,result))
    grid = np.stack(np.meshgrid(*[np.arange(-5,6)]*3,indexing="ij"),axis=-1).reshape(-1,3)
    D = P[None,:,None,:]-P[None,None,:,:]+(grid@H)[:,None,None,:]
    expected = sorted(map(tuple,np.round(D[(D*D).sum(-1)<cutoff**2],10).tolist()))
    assert observed(P) == expected
    assert observed(P+np.random.default_rng(4).integers(-5,6,(5,3))@H) == expected


@pytest.mark.parametrize("mode",["plain","ewald"])
def test_whole_chain_lattice_wrapping_preserves_energy_and_affine_derivative(mode):
    pk = packer(mode,(0,0,0,0))
    params = np.array([4.7,8.5,88.,25.,-40.,.7,1.])
    P,H = (v[0] for v in pk._place(params[None]))
    E,gP,gH = pk.placed_energy_and_grad(P,H,1.)
    wraps = np.zeros((pk.N,3))
    wraps[pk.n:] = [2,-3,4]
    Ew,gPw,gHw = pk.placed_energy_and_grad(P+wraps@H,H,1.)
    assert Ew == pytest.approx(E,abs=3e-11)
    np.testing.assert_allclose(gPw,gP,atol=3e-11)
    np.testing.assert_allclose(gHw+wraps.T@gPw,gH,atol=3e-11)


def test_placed_energy_owner_refuses_invalid_geometry_and_reversal():
    pk = packer("plain")
    params = np.array([4.7,8.5,88.,25.,-40.,.7,1.])
    P,H = (v[0] for v in pk._place(params[None]))
    for coords,lattice,flip in [(P[:-1],H,1.),(P,H,.5),(P,H,np.nan),(P,-H,1.),(P,H*0.,1.)]:
        with pytest.raises(ValueError):
            pk.placed_energy_and_grad(coords,lattice,flip)
    for cutoff,chunk in [(True,10),(-1.,10),(np.nan,10),(2.,False),(2.,0)]:
        with pytest.raises(ValueError):
            list(PlacedCell(P,H).pair_image_chunks(cutoff,chunk))
