#!/usr/bin/env bash
# The large-cell electrostatics on a GPU: checks and timings on a fresh CUDA machine (for example
# a RunPod pod with one A100, started from any CUDA 12 image).  Run inside the pod:
#
#   bash <(curl -sL https://raw.githubusercontent.com/Goober68/polyfind/claude/polymeric-stable-arrangements-uh06b1/examples/runpod_ewald_bench.sh)
#
# or clone the repository and run  bash examples/runpod_ewald_bench.sh
#
# It installs the package and CuPy (into a virtual environment of its own where the image lets it
# make one), makes sure CuPy can reach the device and its libraries, runs the device tests (the
# Ewald sum, the dipole solve, the batched evaluation of a cell, each against the host), and then
# the benchmark (examples/ewald_device_benchmark.py) at sixteen, thirty-two, sixty-four and a
# hundred and twenty-eight chains on the device, held to the stored two-chain cell and to its own
# energy differences.  The whole run is a few minutes; apart from the installation and the
# building of each cell's tables it is the card's work.  Everything printed is also in
# ewald_bench_<time>/log.txt, and the JSON files beside it are the result.
#
#   EWALD_BENCH_OUT=/some/dir    where the log and the JSON files go (default ewald_bench_<time>/)
#   EWALD_BENCH_CHAINS="16 32"   stop at thirty-two chains
#   EWALD_BENCH_LARGE=0          skip the device-only 128-chain run
#   EWALD_BENCH_HOST_MAX=16      also time the host up to sixteen chains (CPU work, the card
#                                idle meanwhile: minutes per size beyond sixteen; default: not at all)
#   EWALD_BENCH_HOST_THREADS="1 8 0"   the BLAS thread counts the host is timed at (0: every core
#                                the process is allotted, which is the default)
#   EWALD_BENCH_LOOP_FORM=1      also run the loop form with the sums on the device (the organisation
#                                before the batched form: host loops, the card at 5 to 10%)
#   EWALD_BENCH_SETUP_ONLY=1     install, check that CuPy reaches the device, and stop
#   EWALD_BENCH_DRY=1            no GPU needed: the same steps on small cells with the tests'
#                                NumPy-backed stand-in for the device (checks this script itself)
set -euo pipefail

BRANCH="${POLYFIND_BRANCH:-claude/polymeric-stable-arrangements-uh06b1}"
DRY="${EWALD_BENCH_DRY:-0}"
if [ ! -f pyproject.toml ] || ! grep -q '^name = "polyfind"' pyproject.toml; then
  git clone --depth 1 --branch "$BRANCH" https://github.com/Goober68/polyfind.git polyfind
  cd polyfind
fi
OUT="${EWALD_BENCH_OUT:-ewald_bench_$(date +%Y%m%d_%H%M%S)}"
mkdir -p "$OUT"
exec > >(tee "$OUT/log.txt") 2>&1

if [ "$DRY" != "1" ]; then
  nvidia-smi || { echo "no NVIDIA GPU visible (nvidia-smi failed)"; exit 1; }
  CUDA_MAJOR=$(nvidia-smi | grep -oE "CUDA Version: [0-9]+" | grep -oE "[0-9]+" || echo 12)
fi

# A Python of our own, so that nothing the image carries (an older NumPy under its PyTorch, say)
# is replaced under it; where the image cannot make one, its own Python is used.
PY="$(command -v python3 || command -v python)"
if "$PY" -m venv "${EWALD_BENCH_VENV:-.venv_ewald_bench}" >/dev/null 2>&1; then
  PY="$(pwd)/${EWALD_BENCH_VENV:-.venv_ewald_bench}/bin/python"
  "$PY" -m pip install -q -U pip
else
  echo "== no virtual environment could be made here: installing into the image's Python"
  export PIP_BREAK_SYSTEM_PACKAGES=1
fi
"$PY" --version
# CuPy 14 wants NumPy 2.0 to 2.5; SciPy has to be one built for NumPy 2.
"$PY" -m pip install -q -U "numpy>=2.0,<2.6" "scipy>=1.14"
"$PY" -m pip install -q -e ".[dev]"

smoke() {
  "$PY" - <<'PY'
import cupy as cp, numpy as np
from cupyx.scipy.special import erfc
n = cp.cuda.runtime.getDeviceCount()
print("cupy", cp.__version__, "| numpy", np.__version__, "| devices", n, "|", cp.cuda.runtime.getDeviceProperties(0)["name"].decode())
a = cp.random.default_rng(0).standard_normal((300, 400), dtype=cp.float64)
s = a @ a.T + 400.0 * cp.eye(300)                       # cuBLAS
x = cp.linalg.solve(s, cp.ones(300))                    # cuSOLVER
assert float(abs(s @ x - 1.0).max()) < 1e-9
assert abs(float(erfc(cp.asarray(0.5))) - 0.4795001221869535) < 1e-14   # a compiled kernel
b = cp.bincount(cp.asarray(np.array([0, 2, 2])), weights=cp.asarray(np.array([1.0, 2.0, 3.0])), minlength=4)
assert b.dtype == np.float64 and b.get().tolist() == [1.0, 0.0, 5.0, 0.0]
print("CuPy reaches the device, cuBLAS, cuSOLVER and the kernel compiler")
PY
}

if [ "$DRY" != "1" ]; then
  # CuPy's wheels are per CUDA major version: 12 and 13 in the current release, 11 in the 13.x line.
  # The one to take is the one whose libraries the image carries, which need not be the newest the
  # driver would run (three pods of 2026-11-13: driver CUDA 13.2, image CUDA 12.1; the wheel for 13
  # found no libraries and the runtime wheels below had to be downloaded, most of a minute each time).
  TOOLKIT=""
  for v in 13 12 11; do
    # (grep without -q: it reads to the end, so the listing is never cut off under pipefail)
    if ldconfig -p 2>/dev/null | grep "libcublas\.so\.$v" >/dev/null || ls /usr/local/cuda*/lib64/libcublas.so.$v >/dev/null 2>&1; then TOOLKIT=$v; break; fi
  done
  if [ -n "$TOOLKIT" ]; then WHEEL="cupy-cuda${TOOLKIT}x"
  elif [ "$CUDA_MAJOR" -ge 13 ]; then WHEEL=cupy-cuda13x; elif [ "$CUDA_MAJOR" -ge 12 ]; then WHEEL=cupy-cuda12x; else WHEEL=cupy-cuda11x; fi
  echo "== CuPy wheel: $WHEEL (CUDA libraries in the image: ${TOOLKIT:-none found}; the driver runs up to CUDA $CUDA_MAJOR)"
  "$PY" -m pip install -q "$WHEEL"
  if ! smoke; then
    # no CUDA libraries where CuPy looks for them: take NVIDIA's own wheels of them
    echo "== CuPy could not use the device with the libraries it found; installing the CUDA runtime wheels"
    "$PY" -m pip install -q "${WHEEL}[ctk]"
    smoke
  fi
  DEVICE=auto
else
  DEVICE=double
fi

# The package's batched kernels (the screen) go to a GPU in float32 by themselves; the large cells
# are float64 work on the host with their Ewald sums on the device, which is asked for per object.
export POLYFIND_DEVICE=cpu

if [ "${EWALD_BENCH_SETUP_ONLY:-0}" = "1" ]; then
  # for a caller that wants the installation and nothing else (examples/runpod_jobs_pod.sh)
  echo "$PY" > .ewald_bench_python
  echo "== set up in $(pwd): $PY, device $DEVICE"
  exit 0
fi

echo "== device tests =="
"$PY" -m pytest tests/test_ewald_device.py tests/test_mixed_batched.py tests/test_gpu.py -q -k "ewald or double or device or dipole or batched or mixed" -p no:cacheprovider

echo "== benchmark =="
if [ "$DRY" = "1" ]; then CHAINS="${EWALD_BENCH_CHAINS:-2 4}"; else CHAINS="${EWALD_BENCH_CHAINS:-16 32 64}"; fi
# A host baseline is CPU work with the card idle, and so is the loop form the batched form
# replaces; a rented card is paid for by the minute and whoever watches its meter sees an idle
# card.  Both comparisons are on record (docs/PERFORMANCE_REVIEW.md 17.6 and 17.8), so by default
# neither is run: the device is held to the stored two-chain cell and to its own energy
# differences, and everything timed is the card's.  To have them again:
#   EWALD_BENCH_HOST_MAX=16 EWALD_BENCH_HOST_THREADS="1 8 0"     (0 stands for every allotted core)
#   EWALD_BENCH_LOOP_FORM=1
if [ "$DRY" = "1" ]; then HOSTMAX="${EWALD_BENCH_HOST_MAX:-4}"; else HOSTMAX="${EWALD_BENCH_HOST_MAX:-0}"; fi
THREADS="${EWALD_BENCH_HOST_THREADS:-0}"
LOOPS=""
if [ "${EWALD_BENCH_LOOP_FORM:-$DRY}" = "1" ]; then LOOPS="--loop-form"; fi
FIRST="${CHAINS%% *}"
"$PY" examples/ewald_device_benchmark.py --device "$DEVICE" --chains $FIRST --host-max-chains "$HOSTMAX" --host-threads $THREADS $LOOPS \
  --out "$OUT/ewald_device_benchmark_$FIRST.json"
REST="${CHAINS#"$FIRST"}"
if [ -n "${REST// /}" ]; then
  "$PY" examples/ewald_device_benchmark.py --device "$DEVICE" --chains $REST --host-max-chains "$HOSTMAX" --host-threads $THREADS $LOOPS \
    --all-devices --out "$OUT/ewald_device_benchmark.json"
fi
if [ "${EWALD_BENCH_LARGE:-1}" = "1" ]; then
  # no host run at this size (it would take the better part of an hour); the device is held to the
  # stored two-chain cell and to its own energy differences instead
  if [ "$DRY" = "1" ]; then LARGE=8; else LARGE=128; fi
  "$PY" examples/ewald_device_benchmark.py --device "$DEVICE" --chains $LARGE --host-max-chains 0 $LOOPS --out "$OUT/ewald_device_benchmark_$LARGE.json" || true
fi

echo
echo "== done.  The results are the JSON files in $(pwd)/$OUT ; log.txt beside them is this output."
ls -la "$OUT"
