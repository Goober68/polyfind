#!/usr/bin/env bash
# Relax a state held in a row of band_row.sh's table to the minimum it belongs to and say which
# state that is (examples/kink_wide_cell.py settle).  A job of examples/runpod_jobs_pod.sh.
#
#   settle.sh ROW NX NY STATE NAME [mover CELL NAME2 | torsions K:DEG,...]
#
# STATE is a dotted path in the row (I.state, I-B.polished.state, I-B.images.5) or settled:NAME.
# "mover CELL NAME2": first give the moving chain the coordinates it has in the settled record
# NAME2 of the same sequence's CELL (a state found in a smaller cell, looked for in this one).
# "torsions 0:180,6:100": first set those torsions of the moving chain (a state looked for by its
# gauche bonds).
# Records go to jobs/shared/states.json, which starts from docs/predictions/kink_states_belt4.json.
# SETTLE_EXTRA: further arguments for the relaxation (a dry run's --maxiter 2).
set -euo pipefail
ROW="$1"; NX="$2"; NY="$3"; STATE="$4"; NAME="$5"; shift 5
KINK="${STEP_KINK:-0}"
SHARED="$(dirname "$JOB_OUT")/shared"
mkdir -p "$SHARED/ck"
STATES="$SHARED/states.json"
[ -f "$STATES" ] || cp docs/predictions/kink_states_belt4.json "$STATES"
# a row no earlier job of this pod computed is taken from the stored table
if ! { [ -f "$SHARED/steps.json" ] && "$PY" -c "import json,sys; sys.exit(0 if sys.argv[2] in json.load(open(sys.argv[1]))['rows'] else 1)" "$SHARED/steps.json" "$ROW@$KINK ${NX}x${NY} mover 0 ahead host flexible"; }; then
  "$PY" examples/kink_wide_cell.py table --out "$SHARED/steps.json" --row "$ROW@$KINK ${NX}x${NY} mover 0 ahead host flexible"
fi
key() { echo "$ROW@$KINK $1 mover 0 ahead host flexible"; }
SLUG="$(echo "$ROW@$KINK" | sed 's/[^[:alnum:]]/_/g')"
MOVER=()
if [ "${1:-}" = "mover" ]; then MOVER=(--mover-from "$(key "$2")" "$3"); fi
if [ "${1:-}" = "torsions" ]; then MOVER=(--mover-torsions "$2"); fi
"$PY" examples/kink_wide_cell.py settle --source "$SHARED/steps.json" --row "$(key "${NX}x${NY}")" --state "$STATE" --name "$NAME" \
  "${MOVER[@]}" --blocks "$SHARED/ck/${SLUG}.chain_blocks.npz" --out "$STATES" ${SETTLE_EXTRA:-}
cp "$STATES" "$JOB_OUT/"
