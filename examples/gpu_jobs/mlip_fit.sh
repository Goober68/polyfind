#!/usr/bin/env bash
# Fit our own interatomic potential on the card, from one or more seeds at once, and test each
# fit against the DFT anchor (examples/mlip_train.py; docs/MLIP_SECOND_OPINION.md 3.3 and 3.5).
# A job of examples/runpod_jobs_pod.sh.  Arguments:
#
#   mlip_fit.sh NAME EPOCHS SEED [SEED...] [frames FILE...] [prepare ARGUMENT...] [train ARGUMENT...] [also ARGUMENT...]
#
# NAME_sSEED names each fit.  The frames default to the two PVDF files; what follows `prepare` is
# given to `mlip_train.py prepare` (e.g. `prepare --heavy REGEX --hold-out-frames NAME NAME`);
# what follows `train` is
# given to `mlip_train.py train` (e.g. `train --stress-weight 100`); what follows `also` is a
# second set of arguments, fitted from the same seeds as NAMEb_sSEED.  The fits are made
# side by side on the one card (a fit of cells this small does not fill it), then tested one
# after another on the host's processors: held-out errors, beta minus alpha, the antipolar
# pair, beta's elastic constants, and the cells relaxed in the potential written out for DFT
# single points (those are run elsewhere: the pod has no DFT program).
#
# Into $JOB_OUT: each fit's log, its model, its results (mlip_own_NAME_sSEED.json), its relaxed
# cells (probe/), and seconds.txt.  MLIP_DEVICE=cpu and a small EPOCHS try the script without a
# card; MLIP_SKIP_TEST=1 stops after the fits.
set -euo pipefail
NAME="$1"; EPOCHS="$2"; shift 2
SEEDS=(); FRAMES=(); PREPARE=(); TRAIN=(); ALSO=(); SECOND=0
while [ $# -gt 0 ] && [ "$1" != "frames" ] && [ "$1" != "prepare" ] && [ "$1" != "train" ] && [ "$1" != "also" ]; do SEEDS+=("$1"); shift; done
if [ "${1:-}" = "frames" ]; then shift; while [ $# -gt 0 ] && [ "$1" != "prepare" ] && [ "$1" != "train" ] && [ "$1" != "also" ]; do FRAMES+=("$1"); shift; done; fi
if [ "${1:-}" = "prepare" ]; then shift; while [ $# -gt 0 ] && [ "$1" != "train" ] && [ "$1" != "also" ]; do PREPARE+=("$1"); shift; done; fi
if [ "${1:-}" = "train" ]; then shift; while [ $# -gt 0 ] && [ "$1" != "also" ]; do TRAIN+=("$1"); shift; done; fi
if [ "${1:-}" = "also" ]; then shift; ALSO=("$@"); SECOND=1; fi
[ ${#FRAMES[@]} -gt 0 ] || FRAMES=(docs/predictions/mlip/pvdf_frames_pbe_d3bj.extxyz docs/predictions/mlip/pvdf_near_frames_pbe_d3bj.extxyz)
PY="${PY:-python3}"
DEVICE="${MLIP_DEVICE:-cuda}"
THREADS="${MLIP_THREADS:-8}"   # a pod shows far more processors than it is given
export OMP_NUM_THREADS="$THREADS" MKL_NUM_THREADS="$THREADS"
DATA="$JOB_OUT/data"

# mace-torch 0.3.17 asks for Python 3.11 or later and the pod image's is 3.10: an environment of
# its own, with a Python that uv fetches (kept beside the jobs, so a later job finds it made).
if ! "$PY" -c "import sys; sys.exit(0 if sys.version_info >= (3, 11) else 1)"; then
  ENV="${MLIP_ENV:-$(dirname "$JOB_OUT")/shared/mlip_env}"
  if [ ! -x "$ENV/bin/python" ]; then
    mkdir -p "$(dirname "$ENV")"
    "$PY" -m pip install -q uv
    "$PY" -m uv venv -q --seed --python "${MLIP_PYTHON:-3.12}" "$ENV"
  fi
  PY="$ENV/bin/python"
fi
# It asks for torch 2.5 or later; the pod image is CUDA 12.1, and a CUDA 12.1 build of torch runs
# on any driver from 530 on, where the newest wheel may ask for a newer driver.
if ! "$PY" -c "import torch" 2>/dev/null; then
  "$PY" -m pip install -q "torch==${MLIP_TORCH:-2.5.1}" --index-url "https://download.pytorch.org/whl/${MLIP_CUDA:-cu121}"
fi
if ! "$PY" -c "import mace, torch_dftd, ase" 2>/dev/null; then
  "$PY" -m pip install -q "mace-torch==0.3.17" "torch-dftd==0.5.3" ase
fi
"$PY" - "$DEVICE" <<'PY'
import sys
import torch, mace, ase
print("torch", torch.__version__, "mace", mace.__version__, "ase", ase.__version__)
print("card:", torch.cuda.get_device_name(0) if torch.cuda.is_available() else "none")
if sys.argv[1] == "cuda" and not torch.cuda.is_available():
    sys.exit("torch does not reach the card: nothing is fitted")
PY

# MLIP_HOLD_OUT: sources kept whole out of the fit, e.g. "VAVF|VRVV" (see mlip_train.py prepare)
read -r -a HOLD <<< "${MLIP_HOLD_OUT:-}"
"$PY" examples/mlip_train.py prepare --frames "${FRAMES[@]}" --data "$DATA" --hold-out "${HOLD[@]}" "${PREPARE[@]}"
t0=$(date +%s)
FITS=()
for seed in "${SEEDS[@]}"; do FITS+=("${NAME}_s${seed}"); done
if [ "$SECOND" = "1" ]; then for seed in "${SEEDS[@]}"; do FITS+=("${NAME}b_s${seed}"); done; fi
each=$(( THREADS / ${#FITS[@]} )); [ "$each" -ge 1 ] || each=1
for seed in "${SEEDS[@]}"; do
  OMP_NUM_THREADS="$each" MKL_NUM_THREADS="$each" "$PY" examples/mlip_train.py train --data "$DATA" \
    --name "${NAME}_s${seed}" --seed "$seed" --epochs "$EPOCHS" --device "$DEVICE" "${TRAIN[@]}" \
    > "$JOB_OUT/${NAME}_s${seed}.launch.txt" 2>&1 &
  if [ "$SECOND" = "1" ]; then
    OMP_NUM_THREADS="$each" MKL_NUM_THREADS="$each" "$PY" examples/mlip_train.py train --data "$DATA" \
      --name "${NAME}b_s${seed}" --seed "$seed" --epochs "$EPOCHS" --device "$DEVICE" "${ALSO[@]}" \
      > "$JOB_OUT/${NAME}b_s${seed}.launch.txt" 2>&1 &
  fi
done
wait || true
echo "fits: $(( $(date +%s) - t0 )) s for ${#FITS[@]} fits of $EPOCHS epochs on $DEVICE" | tee "$JOB_OUT/seconds.txt"
for fit in "${FITS[@]}"; do
  grep -E "Epoch|Error-table|Default" "$DATA/$fit.train.log" | tail -n 60 > "$JOB_OUT/$fit.epochs.txt" || true
done
[ "${MLIP_SKIP_TEST:-0}" = "1" ] && exit 0

mkdir -p "$JOB_OUT/probe"
for fit in "${FITS[@]}"; do
  cat "$JOB_OUT/$fit.launch.txt"
  model="$DATA/${fit}_stagetwo.model"
  [ -f "$model" ] || model="$DATA/$fit.model"
  if [ ! -f "$model" ]; then echo "$fit: NO MODEL (see $fit.train.log)"; cp "$DATA/$fit.train.log" "$JOB_OUT/" || true; continue; fi
  cp "$model" "$JOB_OUT/$fit.model"
  t1=$(date +%s)
  "$PY" examples/mlip_train.py evaluate --data "$DATA" --model "$model" --name "$fit" --probe "$JOB_OUT/probe" --threads "$THREADS"
  cp "docs/predictions/mlip_own_$fit.json" "$JOB_OUT/"
  echo "$fit tested in $(( $(date +%s) - t1 )) s" | tee -a "$JOB_OUT/seconds.txt"
done
