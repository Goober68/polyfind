#!/usr/bin/env bash
# The cells of a row that DFT is asked about, taken from the last screen job's table and put on a DFT line.
# A job of examples/runpod_jobs_pod.sh:
#
#   screen_pick_dft.sh card|processors start|next ROW [MOST]
#
# examples/screen14.py pick: the row's lowest helix cell, its lowest all-trans cell, then the lowest helix cell and
# the lowest all-trans cell of another packing than those (at most MOST of the four, default all).  No person has
# to read the table first, so the line starts when the relaxations end.
set -uo pipefail
LINE="$1"; MODE="$2"; ROW="$3"; MOST="${4:-4}"
HERE="$(dirname "$0")"
source "$HERE/mlip_env.sh"
TABLE="$( (ls "$(dirname "$JOB_OUT")"/*_screen_card/screen_table.json 2>/dev/null || true) | sort | tail -n 1)"
[ -n "$TABLE" ] || { echo "no screen job on this pod has a table"; exit 1; }
echo "table: $TABLE"
NAMES="$("$PY" examples/screen14.py pick --table "$TABLE" --row "$ROW" | tr ' ' '\n' | head -n "$MOST" | tr '\n' ' ')"
[ -n "${NAMES// /}" ] || { echo "$ROW has no cell at rest to ask DFT about"; exit 1; }
echo "picked: $NAMES"
# shellcheck disable=SC2086
bash "$HERE/screen_dft.sh" "$LINE" "$MODE" $NAMES
