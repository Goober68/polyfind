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


def ewald_backend() -> ModuleType:
    """The array module of the Ewald sums and the induced-dipole solve (:mod:`polyfind.ewald`).

    NumPy unless ``POLYFIND_EWALD_DEVICE=cuda``, which asks for CuPy (and raises if it cannot
    be imported).  This is a separate switch from ``POLYFIND_DEVICE`` on purpose: that one
    turns itself on wherever a device is found and runs the batched screen in float32, while
    this path is float64, serves cells of a few dozen sites as well as of thousands, and
    gives the stored numbers of every large-cell result -- so it moves to a device only when
    asked, and then only for cells large enough to gain (``Ewald.device_min_sites``).
    """
    want = os.environ.get("POLYFIND_EWALD_DEVICE", "cpu").lower()
    if want in ("", "cpu"):
        return _np
    if want == "cuda":
        import cupy as _cp  # type: ignore

        return _cp
    raise ValueError(f"POLYFIND_EWALD_DEVICE must be 'cpu' or 'cuda', got {want!r}")


def erfc_of(xp):
    """The complementary error function of an array module.

    SciPy's for NumPy and ``cupyx.scipy.special``'s for CuPy; any other module (a test's
    stand-in for a device) carries its own as ``xp.erfc``.
    """
    if xp is _np:
        from scipy.special import erfc

        return erfc
    own = getattr(xp, "erfc", None)
    if own is not None:
        return own
    from cupyx.scipy.special import erfc  # type: ignore

    return erfc


def synchronize(xp) -> None:
    """Wait until the device has finished what it was given (nothing to wait for on the host).

    A device call returns when the work is queued, not when it is done, so a clock read
    without this measures the queueing.
    """
    cuda = getattr(xp, "cuda", None)
    if cuda is not None:
        cuda.get_current_stream().synchronize()


def cholesky_solve(xp, a, b):
    """``a x = b`` by Cholesky for a symmetric positive definite ``a`` on a CuPy device, or ``None``.

    ``None`` means "use the general solve": the module is not CuPy, or the factorization failed
    (``a`` is not positive definite).  CuPy's general solve is an LU factorization, measured on an
    A100 at 16.9 ms for 2,880 unknowns and 39.5 ms for 5,760 against 3.1 and 8.6 ms for this one,
    the two solutions agreeing to 2e-15 (docs/predictions/ewald_device_a100_20261113).  By default
    CuPy does not look at the factorization's status -- a failed one returns numbers without a
    word -- so the status is asked for here, which is what makes ``None`` reliable.
    """
    if getattr(xp, "__name__", "") != "cupy":
        return None
    try:
        import cupyx  # type: ignore
        from cupyx import lapack  # type: ignore
    except ImportError:
        return None
    try:
        with cupyx.errstate(linalg="raise"):
            return lapack.posv(a, b)
    except _np.linalg.LinAlgError:
        return None


def gram(xp, m):
    """``m @ m.T`` for a C-ordered ``(N, K)`` array ``m`` of the module ``xp``.

    NumPy's product of an array with its own transpose is the BLAS's symmetric rank-k update
    already.  CuPy's is a general product, twice the arithmetic; its ``cupy.cublas.syrk`` is the
    update, and given the transposed *view* (column-ordered, so nothing is copied) it fills the
    upper triangle of the result and leaves the lower at zero.  Measured on an A100 against the
    general product: 51 ms against 88 for 3,840 x 57,886 and 0.43 s against 0.85 for 7,680 x
    115,706, where six of these are most of an evaluation; 1.7 against 1.8 ms for 960 x 14,486
    (docs/PERFORMANCE_REVIEW.md 17.8).  Any other module takes the plain product.
    """
    if getattr(xp, "__name__", "") != "cupy":
        return m @ m.T
    from cupy import cublas  # type: ignore

    out = cublas.syrk("T", m.T)
    out = out + out.T
    d = xp.arange(out.shape[0])
    out[d, d] *= 0.5
    return out


def default_rng(seed=None):
    """Random generator for the active backend (both expose ``.random(size)``)."""
    xp = get_backend()
    return xp.random.default_rng(seed)


def float_dtype():
    """float32 on GPU (throughput), float64 on CPU (simplicity/accuracy)."""
    return _np.float32 if device_name() == "cuda" else _np.float64
