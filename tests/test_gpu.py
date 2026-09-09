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
