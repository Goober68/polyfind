#!/usr/bin/env bash
# A kink pair's two torsions led one after the other, in a potential on the card
# (examples/mlip_path.py scan --only): a scan that leads both together is an upper bound when
# they do not go together.  A job of examples/runpod_jobs_pod.sh:
#
#   steps_card.sh MODEL POINTS NAME:START:END [NAME:START:END...]
#
# MODEL is a model file of the repository; START and END are frames of the labelled
# fourteen-monomer cells or of the scans stored under docs/predictions/mlip/paths/.  For each
# path and each order (a: the first torsion of each chain first; b: the second first) there are
# two legs: NAME_a1 leads that torsion all the way with the other free and is released;
# NAME_a2 goes on from where that arrived, leading whatever still differs from END.  The two
# orders of every path run side by side.  Into $JOB_OUT: each leg's .json/.extxyz/.log.
set -uo pipefail
MODEL="$1"; POINTS="$2"; shift 2
source "$(dirname "$0")/mlip_env.sh"
export MLIP_DEVICE="${MLIP_DEVICE:-cuda}" OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
POT="$(basename "$MODEL" .model)"
FILES=(docs/predictions/mlip/cells14_frames_pbe_d3bj.extxyz docs/predictions/mlip/paths/*.extxyz)
STEPS="${PATH_STEPS:-300}"
t0=$(date +%s)
one() {  # NAME START END ORDER FIRST...
  local name="$1" start="$2" end="$3" order="$4"; shift 4
  "$PY" examples/mlip_path.py scan --model "$MODEL" --potential "$POT" --frames "${FILES[@]}" --start "$start" --end "$end" \
    --name "${name}_${order}1" --out "$JOB_OUT" --points "$POINTS" --steps "$STEPS" --only "$@" > "$JOB_OUT/${name}_${order}1.log" 2>&1
  "$PY" examples/mlip_path.py scan --model "$MODEL" --potential "$POT" --frames "${FILES[@]}" "$JOB_OUT/${name}_${order}1.extxyz" \
    --start "${name}_${order}1_released" --end "$end" --differ 30 \
    --name "${name}_${order}2" --out "$JOB_OUT" --points "$POINTS" --steps "$STEPS" > "$JOB_OUT/${name}_${order}2.log" 2>&1
}
for spec in "$@"; do
  IFS=: read -r name start end <<< "$spec"
  one "$name" "$start" "$end" a 0 2 &
  one "$name" "$start" "$end" b 1 3 &
done
wait
echo "legs: $(( $(date +%s) - t0 )) s on $MLIP_DEVICE" | tee "$JOB_OUT/seconds.txt"
tail -n 1 "$JOB_OUT"/*.log
exit 0
