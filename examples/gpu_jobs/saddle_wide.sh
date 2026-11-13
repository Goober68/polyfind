#!/usr/bin/env bash
# A saddle of the single-chain step in a wider cell, from a narrower cell's saddle and band
# (examples/kink_wide_cell.py polish).  A job of examples/runpod_jobs_pod.sh.  Arguments:
#
#   saddle_wide.sh SEGMENT NX NY band CELL ALONG... start CELL ALONG...
#
# e.g.  saddle_wide.sh I-B 8 4 band 4x2 a b start 8x2 b : the tangent of the 4x2 band let in along
# a then b, the start the 8x2 cell's settled saddle let in along b.  The start row is read from the
# shared table of earlier jobs if its saddle is there, else from the stored one.
set -euo pipefail
ROW="${STEP_ROW:-vdf10}"
KINK="${STEP_KINK:-0}"
SEG="$1"; NX="$2"; NY="$3"; shift 3
[ "$1" = "band" ] || { echo "expected 'band'"; exit 2; }
BAND="$2"; shift 2
BPATH=()
while [ "$1" != "start" ]; do BPATH+=("$1"); shift; done
START="$2"; shift 2
SPATH=("$@")
SHARED="$(dirname "$JOB_OUT")/shared"
mkdir -p "$SHARED/ck"
SADDLES="$SHARED/saddles_${NX}x${NY}.json"
key() { echo "$ROW@$KINK $1 mover 0 ahead host flexible"; }
ARGS=()
PRIOR="$SHARED/start_${START}.json"   # an earlier job's saddles of the start cell, as a table with that row
if [ -f "$SHARED/saddles_${START}.json" ]; then
  "$PY" - "$SHARED/saddles_${START}.json" "$SHARED/wide_steps.json" "$(key "$START")" "$PRIOR" <<'PY'
import json, sys
saddles, table, key, out = sys.argv[1:5]
row = dict(json.load(open(table))["rows"][key])
row["saddles"] = json.load(open(saddles))["saddles"]
json.dump({"rows": {key: row}}, open(out, "w"))
PY
  ARGS+=(--start-file "$PRIOR")
fi
"$PY" examples/kink_wide_cell.py polish --segment "$SEG" --nx "$NX" --ny "$NY" --band-row "$(key "$BAND")" --band-path "${BPATH[@]}" \
  --start-row "$(key "$START")" --start-path "${SPATH[@]}" "${ARGS[@]}" --blocks "$SHARED/ck/${ROW}_${KINK}.chain_blocks.npz" \
  --checkpoint "$JOB_OUT/running_${SEG}.npy" --gtol "${STEP_GTOL:-1e-4}" --out "$SADDLES"
cp "$SADDLES" "$JOB_OUT/"
