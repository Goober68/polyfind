"""Array backend selection (NumPy on CPU, CuPy on CUDA GPUs).

Only the *batched* kernels in polyfind use this module: packing-energy landscapes
evaluated over thousands of candidate cells at once, Boltzmann sampling of many
chains at once, and batched chain construction. Everything sequential along a
single chain (Viterbi, forward-backward, screw decomposition) stays in NumPy
because the state space is tiny (3-5 states) and the recursion cannot be
parallelised along the chain; a GPU would be slower there.

Selection: ``POLYFIND_DEVICE=cuda`` forces CuPy (import error if unavailable),
``POLYFIND_DEVICE=cpu`` forces NumPy, unset = CuPy if importable else NumPy.

All kernels are written against the subset of the NumPy API that CuPy
implements identically, so the same source runs on both. Use :func:`to_numpy`
before handing results to SciPy or to code that expects host arrays.
"""
from __future__ import annotations

import os
from types import ModuleType

import numpy as _np

_backend: ModuleType | None = None
_backend_name: str | None = None


def _select() -> tuple[ModuleType, str]:
    want = os.environ.get("POLYFIND_DEVICE", "auto").lower()
    if want == "cpu":
        return _np, "cpu"
    try:
        import cupy as _cp  # type: ignore

        if want == "cuda":
            return _cp, "cuda"
        # auto: only use CuPy if a device is actually present
        try:
            if _cp.cuda.runtime.getDeviceCount() > 0:
                return _cp, "cuda"
        except Exception:
            pass
    except Exception:
        if want == "cuda":
            raise
    return _np, "cpu"


def get_backend() -> ModuleType:
    """Return the array module (``numpy`` or ``cupy``)."""
    global _backend, _backend_name
    if _backend is None:
        _backend, _backend_name = _select()
    return _backend


def device_name() -> str:
    get_backend()
    return _backend_name or "cpu"


def set_backend(name: str) -> None:
    """Force a backend at runtime ('cpu' or 'cuda'). Mainly for tests/benchmarks."""
    global _backend, _backend_name
    os.environ["POLYFIND_DEVICE"] = name
    _backend = None
    _backend_name = None
    get_backend()


def to_numpy(a):
    """Bring an array to host memory as a NumPy array."""
    if isinstance(a, _np.ndarray):
        return a
    if hasattr(a, "get"):  # CuPy
        return a.get()
    return _np.asarray(a)


def default_rng(seed=None):
    """Random generator for the active backend (both expose ``.random(size)``)."""
    xp = get_backend()
    return xp.random.default_rng(seed)


def float_dtype():
    """float32 on GPU (throughput), float64 on CPU (simplicity/accuracy)."""
    return _np.float32 if device_name() == "cuda" else _np.float64
