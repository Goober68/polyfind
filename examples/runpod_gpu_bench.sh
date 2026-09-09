#!/usr/bin/env bash
# One-shot GPU validation of polyfind on a fresh CUDA machine (e.g. a RunPod pod
# started from the "runpod/pytorch" or any CUDA 12 image).  Run inside the pod:
#
#   bash <(curl -sL https://raw.githubusercontent.com/Goober68/polyfind/claude/polymeric-stable-arrangements-uh06b1/examples/runpod_gpu_bench.sh)
#
# or clone the repo and run  bash examples/runpod_gpu_bench.sh
#
# It installs the package and CuPy, runs the GPU tests (which compare every
# batched kernel against the NumPy result), then prints CPU and GPU timings
# side by side.
set -euo pipefail

BRANCH="${POLYFIND_BRANCH:-claude/polymeric-stable-arrangements-uh06b1}"
if [ ! -f pyproject.toml ] || ! grep -q '^name = "polyfind"' pyproject.toml; then
  git clone --depth 1 --branch "$BRANCH" https://github.com/Goober68/polyfind.git polyfind
  cd polyfind
fi

nvidia-smi || { echo "no NVIDIA GPU visible (nvidia-smi failed)"; exit 1; }
CUDA_MAJOR=$(nvidia-smi | grep -oE "CUDA Version: [0-9]+" | grep -oE "[0-9]+" || echo 12)
python -m pip install -q -e ".[dev]"
if [ "$CUDA_MAJOR" -ge 12 ]; then python -m pip install -q cupy-cuda12x; else python -m pip install -q cupy-cuda11x; fi
python -c "import cupy; print('cupy', cupy.__version__, 'devices', cupy.cuda.runtime.getDeviceCount()); print(cupy.cuda.runtime.getDeviceProperties(0)['name'].decode())"

echo "== GPU tests =="
POLYFIND_DEVICE=cuda python -m pytest tests/test_gpu.py -q
echo "== full test suite on the GPU backend =="
POLYFIND_DEVICE=cuda python -m pytest tests -q -x

echo "== benchmark: CPU =="
POLYFIND_DEVICE=cpu python examples/benchmark.py
echo "== benchmark: GPU =="
POLYFIND_DEVICE=cuda python examples/benchmark.py
