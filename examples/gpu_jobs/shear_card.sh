#!/usr/bin/env bash
# The cell's part of a path in a potential on the card (examples/mlip_path.py shear): from one
# frame's cell to another's in equal steps, the atoms relaxed in each.  A job of
# examples/runpod_jobs_pod.sh:
#
#   shear_card.sh MODEL POINTS NAME:START:END [NAME:START:END...]
#
# MODEL is a model file of the repository.  START and END are frames of the labelled
# fourteen-monomer cells or of the scans stored under docs/predictions/mlip/paths/.  All the
# scans run side by side.  Into $JOB_OUT: NAME.json/.extxyz/.log.
set -uo pipefail
MODEL="$1"; POINTS="$2"; shift 2
source "$(dirname "$0")/mlip_env.sh"
export MLIP_DEVICE="${MLIP_DEVICE:-cuda}" OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
POT="$(basename "$MODEL" .model)"
FILES=(docs/predictions/mlip/cells14_frames_pbe_d3bj.extxyz docs/predictions/mlip/paths/*.extxyz)
t0=$(date +%s)
for spec in "$@"; do
  IFS=: read -r name start end <<< "$spec"
  "$PY" examples/mlip_path.py shear --model "$MODEL" --potential "$POT" --frames "${FILES[@]}" --start "$start" --end "$end" \
    --name "$name" --out "$JOB_OUT" --points "$POINTS" --steps "${SHEAR_STEPS:-300}" > "$JOB_OUT/$name.log" 2>&1 &
done
wait
echo "shear scans: $(( $(date +%s) - t0 )) s on $MLIP_DEVICE" | tee "$JOB_OUT/seconds.txt"
tail -n 1 "$JOB_OUT"/*.log
exit 0
