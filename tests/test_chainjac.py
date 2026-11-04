"""The closed-form Jacobian of the chain build (polyfind.chainjac) against finite differences."""
import numpy as np
import pytest

from polyfind.chainjac import aligned_block_jacobian, oligomer_and_jacobian, repeat_transform_jacobian
from polyfind.linegroup import _block_coords, _windows, line_group, repeat_chains, repeat_chains_and_jacobian, repeat_rotvec
from polyfind.pack import CrystalPacker, _block_atoms, _template, periodic_chain
from polyfind.polymers import PE, PVDF, THREE_STATE
from polyfind.refine import refine_crystal

T, GP, GM = 0, 1, 2
CASES = [
    (PE, np.array([180.0, 172.0]), np.array([112.0, 114.0])),
    (PVDF, np.array([175.0, -176.0]), np.array([113.0, 116.0])),
    (PVDF, np.array([178.0, 52.0, -177.0, -49.0]), np.array([114.5, 112.0])),
    (PVDF, np.array([175.0, 178.0, -179.0, 60.0, -176.0, -177.0, 179.0, -58.0]), np.array([115.0, 113.0])),
    (PVDF, np.array([179.0, 176.0, 55.0, 178.0, -62.0, 175.0]), np.array([114.0, 116.0])),
]


def _displaced(tors, angs, v, h):
    tp, tm, ap, am = tors.copy(), tors.copy(), angs.copy(), angs.copy()
    if v < len(tors):
        tp[v] += h
        tm[v] -= h
    else:
        ap[v - len(tors)] += h
        am[v - len(tors)] -= h
    return tp, ap, tm, am


@pytest.mark.parametrize("poly,tors,angs", CASES)
def test_oligomer_coordinates_are_bit_identical_and_their_jacobian_matches_finite_differences(poly, tors, angs):
    C, bb, J, owner = oligomer_and_jacobian(poly, tors[None], angs[None])
    assert np.array_equal(C, _block_coords(poly, tors[None], angs[None]))
    N = 5 * len(tors)
    keep = ~((owner == 0) | (owner == N + 2))  # the end atoms' pendants are documented as not covered
    h = 1e-4
    for v in range(len(tors) + len(angs)):
        tp, ap, tm, am = _displaced(tors, angs, v, h)
        fd = (_block_coords(poly, tp[None], ap[None])[0] - _block_coords(poly, tm[None], am[None])[0]) / (2 * h)
        np.testing.assert_allclose(J[0, v][keep], fd[keep], atol=1e-7)


@pytest.mark.parametrize("poly,tors,angs", CASES)
def test_repeat_transform_and_aligned_coordinates_match_finite_differences(poly, tors, angs):
    P = len(tors)
    C, bb, J, _ = oligomer_and_jacobian(poly, tors[None], angs[None])
    tpl = _template(poly, P)
    sx, sy = _windows(P)
    tr = repeat_transform_jacobian(bb, J[:, :, tpl.backbone], sx, sy)
    # the rotation vector is the one repeat_rotvec computes, bit for bit
    assert np.array_equal(np.degrees(tr["w"]), repeat_rotvec(poly, tors[None], angs[None]))
    ref = repeat_chains(poly, "x", tors[None], angs)[0]
    chains, dX, dc, w, dw = repeat_chains_and_jacobian(poly, "x", tors[None], angs, ref)
    plain = repeat_chains(poly, "x", tors[None], angs, align_to=ref)[0]
    assert np.array_equal(chains[0].coords, plain.coords) and chains[0].c == plain.c
    assert np.linalg.norm(w[0]) == pytest.approx(plain.rotation_error, abs=1e-9)
    h = 1e-4
    for v in range(P + len(angs)):
        tp, ap, tm, am = _displaced(tors, angs, v, h)
        fdw = (repeat_rotvec(poly, tp[None], ap[None])[0] - repeat_rotvec(poly, tm[None], am[None])[0]) / (2 * h)
        np.testing.assert_allclose(dw[0, v], fdw, atol=1e-7)
        cp = repeat_chains(poly, "x", tp[None], ap, align_to=ref)[0]
        cm = repeat_chains(poly, "x", tm[None], am, align_to=ref)[0]
        assert dc[0, v] == pytest.approx((cp.c - cm.c) / (2 * h), abs=1e-8)
        np.testing.assert_allclose(dX[0, v], (cp.coords - cm.coords) / (2 * h), atol=1e-7)


@pytest.mark.parametrize("seq", [[T, T], [T, GP, T, GM], [T, T, T, GP, T, T, T, GM], [T, GP] * 3])
def test_closure_jacobian_is_the_finite_difference_one_and_the_solve_needs_no_extra_builds(seq):
    chain = periodic_chain(PVDF, seq, THREE_STATE)
    lg = line_group(PVDF, chain.name, chain.dihedrals, states=PVDF.states, bond_angles=np.array([114.0, 114.0]))
    rng = np.random.default_rng(1)
    Y = np.repeat(lg.y0[None], 3, axis=0)
    Y[:, lg.free] += rng.uniform(-10, 10, (3, lg.n_params))
    r, J = lg._residual_and_jac(Y)
    assert np.array_equal(r, lg._residual(Y))
    np.testing.assert_allclose(J[:, :, lg.dep], lg._jacobian(Y), atol=1e-6, rtol=1e-6)
    lg.n_builds = 0
    Yc, it, Jc = lg._close(Y)
    tors, angs = lg._split(Yc)
    assert np.linalg.norm(repeat_rotvec(PVDF, tors, angs), axis=1).max() < 1e-9
    assert it < lg.max_newton
    assert lg.n_builds <= it + 1 + 3  # one build per residual; at most a few damping trials


def test_analytic_route_reproduces_the_geometry_difference_route_and_costs_one_build():
    chain = periodic_chain(PVDF, [T, GP, T, GM], THREE_STATE)
    start = CrystalPacker(chain).result(np.array([5.29, 8.99, 90.0, 209.5, 150.0, 0.40, 1.0]))
    for par in ("linegroup", "free"):
        probe = {}
        refine_crystal(PVDF, start, maxfev=1, maxiter=1, parametrisation=par, gamma_free=True, probe=probe)
        x = probe["x0"] + 0.5
        E, g = probe["value_and_grad_analytic"](x)
        Eg, gg = probe["value_and_grad_fdgeom"](x)
        assert E == Eg == probe["value"](x)
        np.testing.assert_allclose(g, gg, atol=2e-5 * max(np.abs(g).max(), 1e-6))
    with pytest.raises(ValueError):
        refine_crystal(PVDF, start, gradient="fd-geom")
