#!/usr/bin/env bash
# Cells of the screen relaxed on this pod (examples/gpu_jobs/screen_card.sh) put in line for a DFT single point each, by
# name: the relaxed cells are under the pod's served directory, which a queue line cannot name.  A job of
# examples/runpod_jobs_pod.sh:
#
#   screen_dft.sh card|processors start|next NAME [NAME...]
#
# NAME are frames of the last screen_card.sh job's screen_relaxed.extxyz that holds them all
# (scr_vdf8-trfe4-cfe1-rvdf1-14_tt_low_in_all_s1).  `card` is the line of the DFT program built for the card (one
# process), `processors` the other (eight processes in two pools); start and next are examples/gpu_jobs/dft_cells.sh's:
# the end of the line, or its front.
set -euo pipefail
LINE="$1"; MODE="$2"; shift 2
HERE="$(dirname "$0")"
source "$HERE/mlip_env.sh"
FILE=""
for f in $( (ls "$(dirname "$JOB_OUT")"/*_screen_card/screen_relaxed.extxyz 2>/dev/null || true) | sort -r); do
  if "$PY" - "$f" "$@" <<'PY'
import sys
from ase.io import read
names = {str(a.info["frame"]) for a in read(sys.argv[1], index=":")}
sys.exit(0 if all(n in names for n in sys.argv[2:]) else 1)
PY
  then FILE="$f"; break; fi
done
[ -n "$FILE" ] || { echo "no screen job on this pod has relaxed all of: $*"; exit 1; }
if [ "$LINE" = "card" ]; then
  bash "$HERE/dft_cells.sh" in=dft_card pw=/opt/qe_gpu/bin/pw.x "$MODE" 1 1 "$FILE" "$@"
else
  bash "$HERE/dft_cells.sh" "$MODE" 8 2 "$FILE" "$@"
fi
