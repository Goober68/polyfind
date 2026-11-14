#!/usr/bin/env bash
# A saddle of the single-chain step in a wider cell, from a narrower cell's saddle and band
# (examples/kink_wide_cell.py polish).  A job of examples/runpod_jobs_pod.sh.  Arguments:
#
#   saddle_wide.sh [row ROW] SEGMENT NX NY band CELL ALONG... start CELL ALONG... [as TAG]
#
# e.g.  saddle_wide.sh I-B 8 4 band 4x2 a b start 8x2 b : the tangent of the 4x2 band let in along
# a then b, the start the 8x2 cell's settled saddle let in along b.  The start row is read from the
# shared table of earlier jobs if its saddle is there, else from the stored one.  With "as TAG" the
# result goes to saddles_NXxNY_TAG.json, beside the cell's own file and not into it: a second route
# to a saddle that is there already.  ROW is the sequence row (default vdf10); a band or a start
# that band_row.sh computed (jobs/shared/steps.json) is taken from there, else from the stored table.
set -euo pipefail
ROW="${STEP_ROW:-vdf10}"
KINK="${STEP_KINK:-0}"
if [ "${1:-}" = "row" ]; then ROW="$2"; shift 2; fi
SEG="$1"; NX="$2"; NY="$3"; shift 3
[ "$1" = "band" ] || { echo "expected 'band'"; exit 2; }
BAND="$2"; shift 2
BPATH=()
while [ "$1" != "start" ]; do BPATH+=("$1"); shift; done
START="$2"; shift 2
SPATH=(); TAG=""
while [ $# -gt 0 ]; do
  if [ "$1" = "as" ]; then TAG="_$2"; shift 2; else SPATH+=("$1"); shift; fi
done
SHARED="$(dirname "$JOB_OUT")/shared"
mkdir -p "$SHARED/ck"
NAME="saddles"; [ "$ROW" = "vdf10" ] || NAME="saddles_${ROW}"   # PVDF's files keep the names they had
SADDLES="$SHARED/${NAME}_${NX}x${NY}${TAG}.json"
SLUG="$(echo "$ROW@$KINK" | sed 's/[^[:alnum:]]/_/g')"            # the name kink_migration.py gives the row's curvature blocks
key() { echo "$ROW@$KINK $1 mover 0 ahead host flexible"; }
ARGS=()
if [ -f "$SHARED/steps.json" ] && "$PY" -c "import json,sys; sys.exit(0 if 'A-I' in json.load(open(sys.argv[1]))['rows'].get(sys.argv[2], {}) else 1)" "$SHARED/steps.json" "$(key "$BAND")"; then
  ARGS+=(--source "$SHARED/steps.json")   # the band (and a start in the same cell) from band_row.sh's table
fi
PRIOR="$SHARED/start_${NAME}_${START}.json"   # an earlier job's saddles of the start cell, as a table with that row
if [ -f "$SHARED/${NAME}_${START}.json" ]; then
  "$PY" - "$SHARED/${NAME}_${START}.json" "$SHARED/wide_steps.json" "$(key "$START")" "$PRIOR" <<'PY'
import json, sys
saddles, table, key, out = sys.argv[1:5]
row = dict(json.load(open(table))["rows"][key])
row["saddles"] = json.load(open(saddles))["saddles"]
json.dump({"rows": {key: row}}, open(out, "w"))
PY
  ARGS+=(--start-file "$PRIOR")
fi
"$PY" examples/kink_wide_cell.py polish --segment "$SEG" --nx "$NX" --ny "$NY" --band-row "$(key "$BAND")" --band-path "${BPATH[@]}" \
  --start-row "$(key "$START")" --start-path "${SPATH[@]}" "${ARGS[@]}" --blocks "$SHARED/ck/${SLUG}.chain_blocks.npz" \
  --checkpoint "$JOB_OUT/running_${SEG}.npy" --gtol "${STEP_GTOL:-1e-4}" --out "$SADDLES"
cp "$SADDLES" "$JOB_OUT/"
