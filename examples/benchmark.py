"""Timing of the batched kernels on the active backend.

    python examples/benchmark.py            # auto: CuPy if a CUDA device is present, else NumPy
    POLYFIND_DEVICE=cpu python examples/benchmark.py
    POLYFIND_DEVICE=cuda python examples/benchmark.py

Prints the device, then per-kernel throughput, so CPU and GPU runs can be
compared line by line.  The sequential recursions (Viterbi, forward pass) are
timed too, to show that they are already negligible and gain nothing from a GPU.
"""
import time

import numpy as np

from polyfind import backend as bk
from polyfind.amorphous import sample_ensemble
from polyfind.chain import build_chain
from polyfind.forcefield import SimpleFF
from polyfind.lattice_table import PairTable, physical_cores, shutdown_table_pool
from polyfind.pack import SCREEN_TABLE, CrystalPacker, periodic_chain
from polyfind.polymers import PVDF, THREE_STATE, PE
from polyfind.ris import RISModel, polyethylene_like_model


def timeit(fn, repeat=3):
    fn()  # warm-up (JIT/kernel compilation on GPU, cache on CPU)
    best = float("inf")
    for _ in range(repeat):
        if bk.device_name() == "cuda":
            bk.get_backend().cuda.Stream.null.synchronize()
        t0 = time.perf_counter()
        fn()
        if bk.device_name() == "cuda":
            bk.get_backend().cuda.Stream.null.synchronize()
        best = min(best, time.perf_counter() - t0)
    return best


def main():
    print(f"backend: {bk.device_name()} ({bk.get_backend().__name__} {bk.get_backend().__version__}), float dtype {bk.float_dtype().__name__}")
    T, GP, GM = 0, 1, 2
    m = polyethylene_like_model()
    m_pvdf = RISModel(THREE_STATE, 2, np.zeros((2, 3)), np.zeros((2, 3, 3)))

    print("\nsequential recursions (always NumPy, single core):")
    for N in (1000, 10000):
        print(f"  Viterbi minimum, {N:6d} bonds: {1e3 * timeit(lambda: m.minimum(N)):7.1f} ms")
    print(f"  forward pass, 1000 bonds:      {1e3 * timeit(lambda: m.forward(1000, 400.0)):7.1f} ms")

    print("\nexact Boltzmann sampling (backward pass batched over chains):")
    for n in (10_000, 100_000):
        rng = bk.default_rng(0)
        print(f"  {n:7d} chains x 200 bonds: {timeit(lambda: m.sample(200, 400.0, n, rng=rng), 2):6.2f} s")

    print("\nbatched chain building (NeRF, vectorised over chains):")
    print(f"  10,000 PVDF backbones x 200 bonds: {timeit(lambda: sample_ensemble(PVDF, m_pvdf, 200, 10_000, 450.0, rng=bk.default_rng(0)), 2):6.2f} s (incl. sampling)")

    print("\nforce-field scan batch (fit_ris workload):")
    ff = SimpleFF()
    rng = np.random.default_rng(0)
    structs = [build_chain(PVDF, rng.uniform(-180, 180, 12)) for _ in range(2000)]
    print(f"  2,000 PVDF 12-bond oligomers: {1e3 * timeit(lambda: ff.energy_batch(structs)):7.1f} ms")

    print("\ncrystal packing kernel (lattice energy of M cells at once):")
    for poly, seq, lbl in [(PE, [T], "PE all-trans   (6 atoms/chain)"), (PVDF, [T, GP, T, GM], "PVDF alpha    (12 atoms/chain)"), (PVDF, [T, T, T, GP, T, T, T, GM], "PVDF gamma    (24 atoms/chain)")]:
        ch = periodic_chain(poly, seq, THREE_STATE)
        pk = CrystalPacker(ch)
        for M in (512, 4096):
            r = np.random.default_rng(0)
            params = np.column_stack([r.uniform(4.5, 10, M), r.uniform(4.5, 10, M), np.full(M, 90.0), r.uniform(0, 360, M), r.uniform(0, 360, M), r.uniform(0, ch.c, M), r.integers(0, 2, M)])
            dt = timeit(lambda: pk.energy(params), 2)
            print(f"  {lbl}: M={M:5d}: {dt:6.2f} s  = {1e3 * dt / M:6.2f} ms/cell")

    # The table build is the packing path's one real cost, and it is the stage process
    # parallelism helps: the radial axis splits cleanly and the result is bit-identical,
    # so the only question is how far the machine's memory bandwidth lets it scale.
    cores = physical_cores()
    counts = sorted({1, 2, max(2, cores // 2), cores, cores * 2})
    print(f"\nchain-pair table build, screen grid ({cores} physical cores; radial axis split across processes):")
    for poly, seq, lbl in [(PE, [T], "PE all-trans  "), (PVDF, [T, T, T, GP, T, T, T, GM], "PVDF gamma    ")]:
        ch = periodic_chain(poly, seq, THREE_STATE)
        ref, base, row = None, None, []
        for p in counts:
            best = float("inf")
            for _ in range(2 if p > 1 else 1):  # the second pass has the pool already up
                t0 = time.perf_counter()
                tab = PairTable.build(ch, n_procs=p, **SCREEN_TABLE)
                best = min(best, time.perf_counter() - t0)
            if ref is None:
                ref, base = tab.W, best
            elif not np.array_equal(ref, tab.W):
                row.append(f"{p}: NOT BIT-IDENTICAL")
                continue
            row.append(f"{p:2d} proc {best:6.2f} s ({base / best:4.2f}x)")
        print(f"  {lbl} {tab.nbytes / 1e6:3.0f} MB:  " + "  ".join(row))
    shutdown_table_pool()


if __name__ == "__main__":
    main()
