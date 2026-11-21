#!/usr/bin/env bash
# A fit's energies on the labelled fourteen-monomer cells, beside DFT's (examples/mlip_train.py
# score).  A job of examples/runpod_jobs_pod.sh, after refit.sh:
#
#   refit_score.sh FIT [HELD...]
#
# FIT is a fit's name (refit_s1: the model the last refit.sh left) or a model file of the
# repository; HELD are the frames the fit did not see.  Into $JOB_OUT: cells14_in_FIT.json and
# the frames with the fit's energy on each (cells14_in_FIT.extxyz).
set -euo pipefail
FIT="$1"; shift
source "$(dirname "$0")/mlip_env.sh"
export MLIP_DEVICE="${MLIP_DEVICE:-cuda}" OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
if [ -f "$FIT" ]; then MODEL="$FIT"; FIT="$(basename "$FIT" .model)"; else MODEL="$(ls "$(dirname "$JOB_OUT")"/*_refit/"$FIT.model" | tail -n 1)"; fi
"$PY" examples/mlip_train.py score --model "$MODEL" --name "$FIT" --frames docs/predictions/mlip/cells14_frames_pbe_d3bj.extxyz \
  --held-out "$@" --out "$JOB_OUT/cells14_in_$FIT.json" --threads 1
