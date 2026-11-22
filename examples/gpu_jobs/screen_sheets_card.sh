#!/usr/bin/env bash
# A screened row's kink sheets at every site, relaxed in the potential on a pod's card (examples/screen14.py sheets,
# sheets-table; docs/MLIP_SECOND_OPINION.md 3.25).  A job of examples/runpod_jobs_pod.sh:
#
#   screen_sheets_card.sh STARTS SIDE STEPS BUDGET
#
# STARTS is a file of sheet starts in the repository (28 a row).  They are dealt in turn to SIDE processes side by side,
# each relaxing its cells one after another to 0.004 eV/A or STEPS steps; after BUDGET seconds whatever is still running
# is stopped.  The relaxed sheets go where the screen's other cells are (<served>/jobs/shared/screen/relaxed), so that
# each is measured against its row's lowest all-trans cell there.  Into $JOB_OUT: sheets_table.json and
# sheets_relaxed.extxyz.  Nothing is put to DFT: which sheet is asked about is chosen from the table.
set -uo pipefail
STARTS="$1"; SIDE="$2"; STEPS="$3"; BUDGET="$4"
HERE="$(dirname "$0")"
source "$HERE/mlip_env.sh"
export MLIP_DEVICE="${MLIP_DEVICE:-cuda}" OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
MODEL=docs/predictions/mlip/models/all_s1.model
RELAXED="$(dirname "$JOB_OUT")/shared/screen/relaxed"
mkdir -p "$RELAXED"
t0=$(date +%s)
mapfile -t PICKS < <("$PY" - "$STARTS" "$SIDE" <<'PY'
import sys
from ase.io import read
names = [str(a.info["frame"]) for a in read(sys.argv[1], index=":")]
side = int(sys.argv[2])
for k in range(side):
    if names[k::side]:
        print("(" + "|".join(names[k::side]) + ")")
PY
)
PIDS=()
for k in "${!PICKS[@]}"; do
  timeout "$BUDGET" "$PY" examples/mlip_train.py relax --model "$MODEL" --name all_s1 --cells "$STARTS" --pick "${PICKS[$k]}" \
    --out "$RELAXED" --threads 1 --fmax 0.004 --steps "$STEPS" --every 100 --skip-done > "$JOB_OUT/relax_$k.log" 2>&1 &
  PIDS+=($!)
done
for pid in "${PIDS[@]}"; do wait "$pid" || true; done
echo "relaxed for $(( $(date +%s) - t0 )) s on $MLIP_DEVICE, $SIDE side by side" | tee "$JOB_OUT/seconds.txt"
grep -h "steps, " "$JOB_OUT"/relax_*.log || true
"$PY" examples/screen14.py sheets-table --starts "$STARTS" --relaxed "$RELAXED" --out "$JOB_OUT/sheets_table.json" --frames-out "$JOB_OUT/sheets_relaxed.extxyz"
