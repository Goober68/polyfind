#!/usr/bin/env bash
# Does the card give the potential's numbers?  Stored cells (each carries the energy the potential
# gave it on a processor) are evaluated on the card: examples/mlip_train.py device-check.
# A job of examples/runpod_jobs_pod.sh:
#
#   card_check.sh [FILE]      default: the eighteen fourteen-monomer cells relaxed tight
#
# Into $JOB_OUT: card_check.json (the differences in meV per monomer and the seconds an evaluation
# takes).  The job fails when a cell differs by more than 0.01 meV per monomer.  MLIP_DEVICE=cpu
# tries the script without a card.
set -euo pipefail
FILE="${1:-docs/predictions/mlip/sheets14_sheared_in_all_s1.extxyz}"
source "$(dirname "$0")/mlip_env.sh"
export MLIP_DEVICE="${MLIP_DEVICE:-cuda}"
nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader 2>/dev/null || echo "no card visible"
"$PY" examples/mlip_train.py device-check --model docs/predictions/mlip/models/all_s1.model --cells "$FILE" \
  --out "$JOB_OUT/card_check.json" --threads "${MLIP_THREADS:-1}"
