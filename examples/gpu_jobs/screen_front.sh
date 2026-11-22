#!/usr/bin/env bash
# The lowest cells of rows of the screen, as the last screen job's table has them, put at the FRONT of a DFT line.
# A job of examples/runpod_jobs_pod.sh, for when a later file of starts may have found a row a lower cell than the
# one in line for it:
#
#   screen_front.sh card|processors hx|tt|both ROW [ROW...]
#
# For each ROW (vdf9-trfe3-cfe1-rvdf1-14) its lowest helix cell at rest (hx), its lowest all-trans cell (tt), or
# both with the helix first.  A cell that a line has had already (waiting, running or done) is left where it is.
set -euo pipefail
LINE="$1"; WHICH="$2"; shift 2
HERE="$(dirname "$0")"
source "$HERE/mlip_env.sh"
SHARED="$(dirname "$JOB_OUT")/shared"
TABLE="$( (ls "$(dirname "$JOB_OUT")"/*_screen_card/screen_table.json 2>/dev/null || true) | sort | tail -n 1)"
[ -n "$TABLE" ] || { echo "no screen job on this pod has a table"; exit 1; }
INLINE="$(cat "$SHARED"/dft/line.txt "$SHARED"/dft_card/line.txt 2>/dev/null | awk '{print $3}' | tr '\n' ' ')"
# last named is first run: `next` puts each in front of what is there
NAMES="$("$PY" - "$TABLE" "$WHICH" "$@" <<'PY'
import json, sys
rows = json.load(open(sys.argv[1]))["rows"]
which = sys.argv[2]
out = []
for name in sys.argv[3:][::-1]:
    row = rows[name]
    for key in {"hx": ["lowest_hx"], "tt": ["lowest_tt"], "both": ["lowest_tt", "lowest_hx"]}[which]:
        kind = row.get(key)
        if kind:
            out.append(row["cells"][kind]["frame"])
            print(f"# {name}: {key} is {kind} ({row['cells'][kind]['frame']}), margin {row.get('margin', float('nan')):+.3f} in the potential", file=sys.stderr)
print(" ".join(out))
PY
)"
echo "table: $TABLE"
for name in $NAMES; do
  case " $INLINE " in
    *" $name "*) echo "$name: a line has had it already" ;;
    *) bash "$HERE/screen_dft.sh" "$LINE" next "$name" ;;
  esac
done
