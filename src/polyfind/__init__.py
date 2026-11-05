"""polyfind: efficient search for stable atomic arrangements in semi-crystalline polymers."""


def _limit_blas_threads() -> None:
    """One BLAS thread per process, which is what this package's workload wants.

    polyfind's linear algebra is many small operations (360 x 360 polarization solves,
    (120, 2000) @ (2000, 120) Ewald products), and OpenBLAS's threads make those *slower*:
    measured 2026-11-05 on a four-core box, a 360 x 360 ``numpy.linalg.solve`` took 400-1000
    ms with the default threads and 2 ms with one, the spin-waiting threads contending with
    each other and with any other process on the machine.  Parallelism in this package comes
    from independent relaxations run in separate processes, not from BLAS threads.
    ``POLYFIND_BLAS_THREADS`` overrides the limit (``0`` leaves the library's default).
    """
    import os
    want = os.environ.get("POLYFIND_BLAS_THREADS", "1")
    try:
        n = int(want)
    except ValueError:
        return
    if n <= 0:
        return
    # the environment variables reach every BLAS library loaded after this point (scipy's
    # own OpenBLAS loads with scipy.linalg); threadpoolctl reaches the ones already loaded
    for var in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ.setdefault(var, str(n))
    try:
        import threadpoolctl
        threadpoolctl.threadpool_limits(limits=n, user_api="blas")
    except Exception:
        pass


_limit_blas_threads()

from .polymers import PVDF, PE, Monomer, Polymer, copolymer, get_polymer  # noqa: E402
from .ris import RISModel

__all__ = ["PVDF", "PE", "Monomer", "Polymer", "copolymer", "get_polymer", "RISModel"]
__version__ = "0.1.0"
