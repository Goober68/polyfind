#!/usr/bin/env bash
# The single-chain step's end states I and B in a wider cell, started from a narrower one
# (examples/kink_wide_cell.py seed, then examples/kink_migration.py --surface --seed-from).
# A job of examples/runpod_jobs_pod.sh.  Arguments:
#
#   step_wide.sh FROM ALONG... to NX NY [PIECES]
#
# FROM is the narrower cell ("4x2"), ALONG the doublings ("a", "b", or both in order), NX NY the
# cell they lead to, PIECES what to relax (default "I,B").  The narrower row is read from the
# shared table of earlier jobs if it is there, else from the stored one; the result goes to the
# shared table (jobs/shared/wide_steps.json), so a later job can start from it.
set -euo pipefail
ROW="${STEP_ROW:-vdf10}"
KINK="${STEP_KINK:-0}"
FROM="$1"; shift
ALONG=()
while [ "$1" != "to" ]; do ALONG+=("$1"); shift; done
shift
NX="$1"; NY="$2"; PIECES="${3:-I,B}"
SHARED="$(dirname "$JOB_OUT")/shared"
mkdir -p "$SHARED/ck"
TABLE="$SHARED/wide_steps.json"
KEY="$ROW@$KINK $FROM mover 0 ahead host flexible"
SOURCE="docs/predictions/kink_migration_belt4.json"
if [ -f "$TABLE" ] && "$PY" -c "import json,sys; sys.exit(0 if sys.argv[2] in json.load(open(sys.argv[1]))['rows'] and 'B' in json.load(open(sys.argv[1]))['rows'][sys.argv[2]] else 1)" "$TABLE" "$KEY"; then
  SOURCE="$TABLE"
fi
SEED="seed ${NX}x${NY} from $FROM"
echo "== ${NX}x${NY} from the $FROM row of $SOURCE, let in along ${ALONG[*]}; pieces $PIECES"
"$PY" examples/kink_wide_cell.py seed --source "$SOURCE" --row "$KEY" --along "${ALONG[@]}" --name "$SEED" --out "$TABLE"
"$PY" examples/kink_migration.py --row "$ROW" --kink-at "$KINK" --nx "$NX" --ny "$NY" --mover 0 --host flexible --surface --trust-tiled \
  --ends-only --only "$PIECES" --seed-from "$SEED" --patience "${STEP_PATIENCE:-40}" --surface-gtol "${STEP_GTOL:-3e-4}" \
  --checkpoints "$SHARED/ck" --out "$TABLE"
"$PY" - "$TABLE" "$ROW@$KINK ${NX}x${NY} mover 0 ahead host flexible" "$JOB_OUT/result.json" <<'PY'
import json, sys
table, key, out = sys.argv[1:4]
row = json.load(open(table))["rows"][key]
res = {"row": key, "chains": row["chains"]}
for name in ("I", "B"):
    if name in row:
        res[name] = {k: v for k, v in row[name].items() if k != "state"}
        print(f"   {name}: {row[name]['energy']:+.4f} kcal/mol in {row['chains']} chains; {row[name]['evaluations']} evaluations, {row[name]['seconds']:.0f} s; "
              f"gradient {row[name]['gradient_max']:.1e}; {row[name]['torsions']}")
json.dump(res, open(out, "w"), indent=1)
PY
