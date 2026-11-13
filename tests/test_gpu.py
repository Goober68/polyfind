"""GPU backend checks. Skipped unless CuPy and a CUDA device are available.

Run on a GPU machine with:  POLYFIND_DEVICE=cuda pytest tests/test_gpu.py
These compare every batched kernel against the NumPy result.
"""
import numpy as np
import pytest

from polyfind import backend as bk

try:
    import cupy  # noqa: F401

    _has_cupy = True
    try:
        _has_device = cupy.cuda.runtime.getDeviceCount() > 0
    except Exception:
        _has_device = False
except Exception:
    _has_cupy = _has_device = False

pytestmark = pytest.mark.skipif(not (_has_cupy and _has_device), reason="CuPy with a CUDA device is required")


@pytest.fixture
def cuda_backend():
    bk.set_backend("cuda")
    yield
    bk.set_backend("cpu")


def test_backend_selected(cuda_backend):
    assert bk.device_name() == "cuda"
    assert bk.get_backend().__name__ == "cupy"


def test_packing_kernel_matches_cpu(cuda_backend):
    from polyfind.pack import CrystalPacker, periodic_chain
    from polyfind.polymers import PVDF, THREE_STATE

    ch = periodic_chain(PVDF, [0, 1, 0, 2], THREE_STATE)
    rng = np.random.default_rng(0)
    M = 64
    params = np.column_stack([rng.uniform(4.5, 10, M), rng.uniform(4.5, 10, M), np.full(M, 90.0), rng.uniform(0, 360, M), rng.uniform(0, 360, M), rng.uniform(0, ch.c, M), rng.integers(0, 2, M)])
    e_gpu = CrystalPacker(ch).energy(params)
    bk.set_backend("cpu")
    e_cpu = CrystalPacker(ch).energy(params)
    assert np.allclose(e_gpu, e_cpu, rtol=1e-3, atol=1e-2)  # float32 on GPU


def test_sampling_and_build_on_gpu(cuda_backend):
    from polyfind.amorphous import ensemble_stats, sample_ensemble
    from polyfind.polymers import PE
    from polyfind.ris import polyethylene_like_model

    m = polyethylene_like_model()
    states, coords = sample_ensemble(PE, m, 100, 5000, 400.0, rng=bk.default_rng(0))
    assert type(states).__module__.startswith("cupy")
    st = ensemble_stats(PE, m, states, coords, 400.0)
    p = m.marginals(100, 400.0).mean(axis=0)
    assert abs(st.state_fractions["T"] - p[0]) < 0.02


def test_forcefield_batch_on_gpu(cuda_backend):
    from polyfind.chain import build_chain
    from polyfind.forcefield import SimpleFF
    from polyfind.polymers import PVDF

    ff = SimpleFF()
    rng = np.random.default_rng(0)
    structs = [build_chain(PVDF, rng.uniform(-180, 180, 8)) for _ in range(16)]
    e = ff.energy_batch(structs)
    for s, ei in zip(structs, e):
        assert abs(ff.energy(s) - ei) < 1e-3 * max(1.0, abs(ei))


# --------------------------------------------------------------------------------------------
# The Ewald sum and the induced-dipole solve on the device (float64; asked for per Ewald object
# or with POLYFIND_EWALD_DEVICE=cuda, independently of the backend of the batched kernels above).
# ``tests/test_ewald_device.py`` runs the same comparisons everywhere on a NumPy-backed stand-in;
# these are the ones only a real device can make.
@pytest.fixture
def host_kernels():
    """The batched kernels on the host in float64, as a mixed cell needs them, whatever the machine
    would choose by itself; put back afterwards."""
    import os

    before = os.environ.get("POLYFIND_DEVICE")
    bk.set_backend("cpu")
    yield
    if before is None:
        os.environ.pop("POLYFIND_DEVICE", None)
        bk._backend = bk._backend_name = None
    else:
        bk.set_backend(before)


def test_ewald_terms_on_the_device_are_the_hosts():
    import cupy
    from test_ewald_device import _cell, _on, _same

    from polyfind.ewald import EwaldSpec, Thole

    X, q, h, p, scale = _cell(60)
    for boundary in ("tinfoil", "vacuum"):
        spec = EwaldSpec(boundary=boundary, charge_damping=4.0)
        for th in (Thole(), Thole(a=0.39, form="cubic"), None):
            kw = {} if th is None else dict(thole=th, thole_scale=scale)
            for flags in (dict(dipoles=p, grad=True, charge_grad=True, dipole_grad=True, field=True, tensor=True),
                          dict(dipoles=p, grad=True), dict(field=True, tensor=True)):
                _same(_on(np, spec).terms(X, q, h, **flags, **kw), _on(cupy, spec).terms(X, q, h, **flags, **kw), rel=1e-10)
        for flags in (dict(), dict(grad=True, charge_grad=True)):
            _same(_on(np, spec).terms(X, q + 0.01, h, **flags), _on(cupy, spec).terms(X, q + 0.01, h, **flags), rel=1e-10)
    kept = _on(cupy).terms(X, q, h, thole=Thole(), thole_scale=cupy.asarray(scale), tensor=True, device_tensor=True)
    assert isinstance(kept.tensor, cupy.ndarray) and kept.tensor.dtype == np.float64


def test_ewald_phases_are_kept_on_the_device_between_the_two_calls_of_an_evaluation():
    import cupy
    from test_ewald_device import _on, _wide_cell

    X, q, h, p = _wide_cell()
    ew = _on(cupy)
    first = ew.terms(X, q, h)
    kept = ew._phase_cache[2]
    assert isinstance(kept, cupy.ndarray)
    second = ew.terms(X, q, h, dipoles=p, grad=True)
    assert ew._phase_cache[2] is kept
    host = _on(np).terms(X, q, h, dipoles=p, grad=True)
    assert second.total == pytest.approx(host.total, rel=1e-11)
    np.testing.assert_allclose(second.grad_coords, host.grad_coords, rtol=0, atol=1e-9 * float(np.abs(host.grad_coords).max()))
    assert first.total == pytest.approx(_on(np).terms(X, q, h).total, rel=1e-11)
    ew.release()
    assert ew._phase_cache is None


def test_ewald_induced_dipole_evaluation_on_the_device_is_the_hosts(host_kernels):
    """Energy and every gradient of a four-chain cell with the sums, the solve and the gradient on
    the device, against the same evaluation on the host."""
    import cupy
    from test_ewald_device import BELT4, _model

    from polyfind.fitting import PRESET_PARAMETERS

    with PRESET_PARAMETERS[BELT4].applied():
        m, x = _model(2, 1)
        x = x + 0.02 * np.random.default_rng(0).standard_normal(x.size) * m.scale(x)
        e_host, g_host = m.value_and_grad(x)
        ew = m.pk._ewald
        ew.xp, ew.device_min_sites = cupy, 0
        try:
            e_dev, g_dev = m.value_and_grad(x)
            assert isinstance(m.pk._thole_scale_dev[2], cupy.ndarray)
        finally:
            ew.xp, ew.device_min_sites = np, 256
            ew.release()
    assert e_dev == pytest.approx(e_host, rel=1e-11)
    np.testing.assert_allclose(g_dev, g_host, rtol=0, atol=1e-9 * max(1.0, float(np.abs(g_host).max())))


def test_ewald_device_is_asked_for_by_its_own_switch(monkeypatch):
    import cupy

    from polyfind.ewald import Ewald

    monkeypatch.setenv("POLYFIND_EWALD_DEVICE", "cuda")
    assert bk.ewald_backend() is cupy and Ewald(8.0).xp is cupy
    monkeypatch.delenv("POLYFIND_EWALD_DEVICE")
    assert Ewald(8.0).xp is np  # whatever the backend of the batched kernels is


# A mixed cell's whole evaluation on the device: the batched form (``tests/test_mixed_batched.py``
# holds it to the loop form on the host and on the stand-in), with CuPy's own functions.
@pytest.mark.parametrize("case", [0, 1, 2, 3])
def test_mixed_cell_on_the_device_is_the_loop_form(host_kernels, case):
    """Placed energy, both gradients, and the model's value and gradient with the pair sums on the
    card: truncated Coulomb and Ewald with induced dipoles, flipped chains, several images along
    short repeats."""
    import cupy
    from test_mixedcell import CASES, _model, _run

    from polyfind import mixedcell as MC

    def body(polymer, seq, params, kw):
        m, x, _ = _model(polymer, seq, params, nx=2, ny=1, pen_rot=0.7, k_c=3.0, **kw)
        x = x + 0.02 * np.random.default_rng(case).standard_normal(x.size) * m.scale(x)
        latn, frames, _, Pn, _ = m.build(x)
        E0, g0, l0 = MC.placed_energy_and_grad(m.pk, frames, Pn, latn)
        e0, v0 = m.value_and_grad(x)
        ew = m.pk._ewald
        try:
            if ew is not None:
                ew.xp, ew.device_min_sites = cupy, 0
            E1, g1, l1 = MC.placed_energy_and_grad(m.pk, frames, Pn, latn, xp=cupy)
            m.xp = cupy
            e1, v1 = m.value_and_grad(x)
            small = MC._pair_sums(m.pk, Pn.reshape(m.nc, m.n, 3), latn, m.flips, int(np.ceil(m.pk.rc / latn[2, 2])) + 1, cupy,
                                  block_elems=2 * m.n * m.n)
            whole = MC._pair_sums(m.pk, Pn.reshape(m.nc, m.n, 3), latn, m.flips, int(np.ceil(m.pk.rc / latn[2, 2])) + 1, cupy)
        finally:
            m.xp = None
            if ew is not None:
                ew.xp, ew.device_min_sites = np, 256
                ew.release()
        assert E1 == pytest.approx(E0, rel=1e-11, abs=1e-9) and e1 == pytest.approx(e0, rel=1e-11, abs=1e-9)
        for a, b in ((g1, g0), (l1, l0), (v1, v0), (small[1], whole[1]), (small[2], whole[2])):
            assert isinstance(a, np.ndarray)
            np.testing.assert_allclose(a, b, rtol=0, atol=1e-10 * max(1.0, float(np.abs(b).max())))
        assert small[0] == pytest.approx(whole[0], rel=1e-12)

    _run(CASES[case], body)


def test_device_cholesky_solve_is_the_general_solve_and_says_when_it_cannot():
    """``polyfind.backend.cholesky_solve``: the general solve's answer for a positive definite
    matrix, and ``None`` -- not numbers -- for one that is not."""
    import cupy

    rng = np.random.default_rng(0)
    M = rng.standard_normal((300, 300))
    a = cupy.asarray(M @ M.T + 300.0 * np.eye(300))
    b = cupy.asarray(rng.standard_normal(300))
    x = bk.cholesky_solve(cupy, a, b)
    assert isinstance(x, cupy.ndarray) and x.shape == (300,)
    np.testing.assert_allclose(x.get(), cupy.linalg.solve(a, b).get(), rtol=0, atol=1e-12)
    indefinite = a - 2000.0 * cupy.eye(300)
    assert bk.cholesky_solve(cupy, indefinite, b) is None
    assert bk.cholesky_solve(np, M @ M.T, rng.standard_normal(300)) is None  # the host has its own solve
    import cupyx

    assert cupyx.geterr()["linalg"] == "ignore"  # the status was asked for inside, and the setting put back
