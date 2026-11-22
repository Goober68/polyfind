#!/usr/bin/env bash
# Cells of the registry search that had not come to rest, relaxed on from where they stood, on the card
# (docs/MLIP_SECOND_OPINION.md 3.24).  A job of examples/runpod_jobs_pod.sh:
#
#   registry_more.sh STARTS SIDE STEPS BUDGET WHICH [WHICH...]
#
# WHICH are names of starts (reg_matrix_tt_p1m1_tp117), or one word auto:KCAL: every cell of the last
# registry job that is not at rest, kept its conformation and is below KCAL kcal/mol per chain against
# its labelled cell, the lowest first.  They are dealt in turn to SIDE processes side by side, each
# relaxing its cells one after another for STEPS more steps at most (to 0.004 eV/A); after BUDGET
# seconds whatever is still running is stopped.  The cells are written where they were (the last
# registry_card.sh job's relaxed/), and examples/gpu_jobs/registry_finish.sh makes the table and the
# choice again and puts the chosen on the DFT lines.
set -uo pipefail
STARTS="$1"; SIDE="$2"; STEPS="$3"; BUDGET="$4"; shift 4
HERE="$(dirname "$0")"
# registry_more.cap (beside this script): lines "NN SECONDS".  A line of a pod's queue cannot be changed once it is
# there; a job whose number is NN gets at most SECONDS, whatever its line says.
CAP="$(awk -v n="$(basename "$JOB_OUT" | cut -c1-2)" '$1 == n {print $2}' "$HERE/registry_more.cap" 2>/dev/null | tail -n 1)"
if [ -n "$CAP" ] && [ "$CAP" -lt "$BUDGET" ]; then echo "this job's time is cut from $BUDGET to $CAP s (registry_more.cap)"; BUDGET="$CAP"; fi
source "$HERE/mlip_env.sh"
export MLIP_DEVICE="${MLIP_DEVICE:-cuda}" OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
MODEL=docs/predictions/mlip/models/all_s1.model
LABELS=docs/predictions/mlip/cells14_frames_pbe_d3bj.extxyz
DIR="$( (ls -d "$(dirname "$JOB_OUT")"/*_registry_card 2>/dev/null || true) | tail -n 1)"
[ -d "$DIR/relaxed" ] || { echo "no relaxed registry cells on this pod"; exit 1; }
t0=$(date +%s)
if [ "${1#auto:}" != "$1" ]; then
  "$PY" examples/registry_cells.py table --starts "$STARTS" --relaxed "$DIR/relaxed" --labels "$LABELS" --out "$JOB_OUT/registry_table_before.json" > "$JOB_OUT/table_before.log" || exit 1
  NAMES=($("$PY" examples/registry_cells.py unsettled --table "$JOB_OUT/registry_table_before.json" --below "${1#auto:}"))
else
  NAMES=("$@")
fi
echo "${#NAMES[@]} cells relaxed on: ${NAMES[*]}"
PIDS=()
for k in $(seq 0 $((SIDE - 1))); do
  : > "$JOB_OUT/more_$k.sh"
  for i in "${!NAMES[@]}"; do
    if [ $((i % SIDE)) -eq "$k" ]; then
      name="${NAMES[$i]}"
      echo "\"$PY\" examples/mlip_train.py relax --model \"$MODEL\" --name all_s1 --cells \"$DIR/relaxed/${name}_in_all_s1.start.xyz\" --frame \"${name}_in_all_s1\" --out \"$DIR/relaxed\" --threads 1 --fmax 0.004 --steps \"$STEPS\" --every 100" >> "$JOB_OUT/more_$k.sh"
    fi
  done
  timeout "$BUDGET" bash "$JOB_OUT/more_$k.sh" > "$JOB_OUT/relax_$k.log" 2>&1 &
  PIDS+=($!)
done
for pid in "${PIDS[@]}"; do wait "$pid" || true; done
echo "relaxed on for $(( $(date +%s) - t0 )) s on $MLIP_DEVICE, $SIDE side by side" | tee "$JOB_OUT/seconds.txt"
grep -h "steps, " "$JOB_OUT"/relax_*.log || true
bash "$HERE/registry_finish.sh" "$STARTS" "$DIR/relaxed"
