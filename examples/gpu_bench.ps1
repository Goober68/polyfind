# GPU validation of polyfind on Windows (PowerShell), from your local clone, e.g.
#   cd e:\source\gh\polyfind
#   git fetch origin claude/polymeric-stable-arrangements-uh06b1
#   git checkout claude/polymeric-stable-arrangements-uh06b1
#   .\examples\gpu_bench.ps1
#
# Installs the package and CuPy (CUDA 12 wheel by default; pass -Cuda 11 for CUDA 11),
# runs the GPU tests (every batched kernel compared with its NumPy result), the full
# test suite on the CuPy backend, then prints CPU and GPU timings side by side.
param(
    [int]$Cuda = 12,
    [switch]$SkipInstall
)
$ErrorActionPreference = "Stop"

nvidia-smi
if ($LASTEXITCODE -ne 0) { throw "nvidia-smi failed: no NVIDIA GPU / driver visible" }

if (-not $SkipInstall) {
    python -m pip install -q -e ".[dev]"
    if ($Cuda -ge 12) { python -m pip install -q cupy-cuda12x } else { python -m pip install -q cupy-cuda11x }
}
python -c "import cupy; print('cupy', cupy.__version__, 'devices', cupy.cuda.runtime.getDeviceCount()); print(cupy.cuda.runtime.getDeviceProperties(0)['name'].decode())"

Write-Host "== GPU tests =="
$env:POLYFIND_DEVICE = "cuda"
python -m pytest tests/test_gpu.py -q
if ($LASTEXITCODE -ne 0) { throw "GPU tests failed" }

Write-Host "== full test suite on the GPU backend =="
python -m pytest tests -q -x
if ($LASTEXITCODE -ne 0) { throw "test suite failed on the GPU backend" }

Write-Host "== benchmark: CPU =="
$env:POLYFIND_DEVICE = "cpu"
python examples/benchmark.py

Write-Host "== benchmark: GPU =="
$env:POLYFIND_DEVICE = "cuda"
python examples/benchmark.py
