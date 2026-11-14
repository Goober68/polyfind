#!/usr/bin/env bash
# The batched form of the evaluation, on the array module the cell chooses (the card), against the
# host's loops, on the rows' own supercells (examples/mixed_forms_check.py).  A job of
# examples/runpod_jobs_pod.sh.  Arguments: rows; the cells are 2x1 and 4x2 (FORMS_CELLS to change).
set -euo pipefail
ROWS=()
for r in "$@"; do ROWS+=(--row "$r"); done
status=0
for cell in ${FORMS_CELLS:-2x1 4x2}; do
  "$PY" examples/mixed_forms_check.py "${ROWS[@]}" --nx "${cell%x*}" --ny "${cell#*x}" --out "$JOB_OUT/forms_$cell.json" || status=1
done
exit $status
