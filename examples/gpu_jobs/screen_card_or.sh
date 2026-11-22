#!/usr/bin/env bash
# A cell of the screen on the card's DFT line unless the processors' line has taken it; then other cells instead.
# A job of examples/runpod_jobs_pod.sh, for when a cell was put at the front of the processors' line by a job that
# may have come too late for the processors' next turn:
#
#   screen_card_or.sh NAME -- [OTHER...]
#
# If the processors' line has started NAME (or labelled it), the card's line is given OTHER (nothing, if none are
# named); if not, the card's line is given NAME.  Names as examples/gpu_jobs/screen_dft.sh takes them.
set -uo pipefail
HERE="$(dirname "$0")"
SHARED="$(dirname "$JOB_OUT")/shared"
NAME="$1"; shift
[ "${1:-}" = "--" ] && shift
if grep -q " $NAME: started" "$SHARED/dft/log.txt" 2>/dev/null || [ -f "$SHARED/dft/$NAME.dft.xyz" ]; then
  echo "the processors have $NAME: the card's line takes: ${*:-nothing}"
  [ $# -gt 0 ] && bash "$HERE/screen_dft.sh" card start "$@"
else
  echo "the processors have not started $NAME: the card's line takes it"
  bash "$HERE/screen_dft.sh" card start "$NAME"
fi
exit 0
