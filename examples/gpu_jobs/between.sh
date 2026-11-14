#!/usr/bin/env bash
# A band between two states of a row of band_row.sh's table that are not its named end points,
# and the saddle settled from its top (examples/kink_wide_cell.py between).  A job of
# examples/runpod_jobs_pod.sh.
#
#   between.sh ROW NX NY START END NAME [plain]
#
# START and END are dotted paths in the row (I.state, B.state) or settled:NAME (a record of
# settle.sh, or of docs/predictions/kink_states_belt4.json).  The same line again continues a band
# that has not converged.  Records go to jobs/shared/states.json.  The band runs in the per-chain
# metric (polyfind.mixedcell.neb, metric) unless "plain" follows.
set -euo pipefail
ROW="$1"; NX="$2"; NY="$3"; START="$4"; END="$5"; NAME="$6"
PLAIN=(); [ "${7:-}" = "plain" ] && PLAIN=(--plain)   # the band in per-variable units, not in the per-chain metric
KINK="${STEP_KINK:-0}"
SHARED="$(dirname "$JOB_OUT")/shared"
mkdir -p "$SHARED/ck"
STATES="$SHARED/states.json"
[ -f "$STATES" ] || cp docs/predictions/kink_states_belt4.json "$STATES"
# a row no earlier job of this pod computed is taken from the stored table
if ! { [ -f "$SHARED/steps.json" ] && "$PY" -c "import json,sys; sys.exit(0 if sys.argv[2] in json.load(open(sys.argv[1]))['rows'] else 1)" "$SHARED/steps.json" "$ROW@$KINK ${NX}x${NY} mover 0 ahead host flexible"; }; then
  "$PY" examples/kink_wide_cell.py table --out "$SHARED/steps.json" --row "$ROW@$KINK ${NX}x${NY} mover 0 ahead host flexible"
fi
SLUG="$(echo "$ROW@$KINK" | sed 's/[^[:alnum:]]/_/g')"
"$PY" examples/kink_wide_cell.py between --source "$SHARED/steps.json" --row "$ROW@$KINK ${NX}x${NY} mover 0 ahead host flexible" \
  --start "$START" --end "$END" --name "$NAME" --maxiter "${STEP_BAND_MAXITER:-300}" --fmax "${STEP_FMAX:-0.1}" \
  --blocks "$SHARED/ck/${SLUG}.chain_blocks.npz" "${PLAIN[@]}" --out "$STATES"
cp "$STATES" "$JOB_OUT/"
