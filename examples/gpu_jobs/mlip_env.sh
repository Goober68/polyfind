# Sourced by a job of examples/runpod_jobs_pod.sh that needs the potential: leaves PY as a Python
# with torch, mace, torch-dftd and ase, made once and kept beside the jobs (examples/gpu_jobs/mlip_fit.sh
# makes the same one the same way).  MLIP_ENV names another place for it.
PY="${PY:-python3}"
# mace-torch 0.3.17 asks for Python 3.11 or later and the pod image's is 3.10: an environment of
# its own, with a Python that uv fetches.
if ! "$PY" -c "import sys; sys.exit(0 if sys.version_info >= (3, 11) else 1)"; then
  MLIP_ENV_DIR="${MLIP_ENV:-$(dirname "$JOB_OUT")/shared/mlip_env}"
  if [ ! -x "$MLIP_ENV_DIR/bin/python" ]; then
    mkdir -p "$(dirname "$MLIP_ENV_DIR")"
    "$PY" -m pip install -q uv
    "$PY" -m uv venv -q --seed --python "${MLIP_PYTHON:-3.12}" "$MLIP_ENV_DIR"
  fi
  PY="$MLIP_ENV_DIR/bin/python"
fi
# It asks for torch 2.5 or later; the pod image is CUDA 12.1, and a CUDA 12.1 build of torch runs
# on any driver from 530 on, where the newest wheel may ask for a newer driver.
if ! "$PY" -c "import torch" 2>/dev/null; then
  "$PY" -m pip install -q "torch==${MLIP_TORCH:-2.5.1}" --index-url "https://download.pytorch.org/whl/${MLIP_CUDA:-cu121}"
fi
if ! "$PY" -c "import mace, torch_dftd, ase" 2>/dev/null; then
  "$PY" -m pip install -q "mace-torch==0.3.17" "torch-dftd==0.5.3" ase
fi
export PY
