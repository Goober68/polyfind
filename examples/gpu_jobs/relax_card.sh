#!/usr/bin/env bash
# Stored cells relaxed in the potential on the card, atoms and lattice, to a tight stop
# (examples/mlip_train.py relax).  A job of examples/runpod_jobs_pod.sh:
#
#   relax_card.sh FILE PATTERN [STEPS] [FMAX]
#
# FILE is a file of frames in the repository; the frames whose whole name matches the regular
# expression PATTERN are relaxed one after another in one process, in the file's order, to FMAX
# eV/A (default 0.004) or STEPS steps (default 3000).  Each cell is written to
# $JOB_OUT/relaxed/NAME_in_all_s1.start.xyz when it ends, and every 250 steps on the way (marked
# NOT converged), so a job that is cut leaves every cell where it had got to.  A cell already
# there and converged is left.  One processor thread: the rest of the pod's are DFT's.
# MLIP_DEVICE=cpu tries the script without a card.
set -euo pipefail
FILE="$1"; PATTERN="$2"; STEPS="${3:-3000}"; FMAX="${4:-0.004}"
source "$(dirname "$0")/mlip_env.sh"
export MLIP_DEVICE="${MLIP_DEVICE:-cuda}"
THREADS="${MLIP_THREADS:-1}"
export OMP_NUM_THREADS="$THREADS" MKL_NUM_THREADS="$THREADS"
t0=$(date +%s)
"$PY" examples/mlip_train.py relax --model docs/predictions/mlip/models/all_s1.model --name all_s1 \
  --cells "$FILE" --pick "$PATTERN" --out "$JOB_OUT/relaxed" --threads "$THREADS" \
  --fmax "$FMAX" --steps "$STEPS" --every 250 --skip-done
echo "relaxed in $(( $(date +%s) - t0 )) s on $MLIP_DEVICE" | tee "$JOB_OUT/seconds.txt"
