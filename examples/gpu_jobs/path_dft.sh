#!/usr/bin/env bash
# Images of a path scanned on this pod (examples/gpu_jobs/path_card.sh, shear_card.sh) put in line for a DFT
# single point each, by the scan's name: the scan's file is under the pod's served directory,
# which a queue line cannot name.  A job of examples/runpod_jobs_pod.sh:
#
#   path_dft.sh card|processors start|next SCAN NAME [NAME...]
#
# SCAN is a scan's name (film_f: the last path_card.sh job's film_f.extxyz); NAME are frames of
# it (film_f_i06, film_f_released).  `card` is the line of the DFT program built for the card
# (one process), `processors` the other (eight processes in two pools); start and next are
# examples/gpu_jobs/dft_cells.sh's: the end of the line, or its front.
set -euo pipefail
LINE="$1"; MODE="$2"; SCAN="$3"; shift 3
HERE="$(dirname "$0")"
FILE="$( (ls "$(dirname "$JOB_OUT")"/*_path_card/"$SCAN.extxyz" "$(dirname "$JOB_OUT")"/*_shear_card/"$SCAN.extxyz" 2>/dev/null || true) | tail -n 1)"
[ -f "$FILE" ] || { echo "no scan $SCAN on this pod"; exit 1; }
if [ "$LINE" = "card" ]; then
  bash "$HERE/dft_cells.sh" in=dft_card pw=/opt/qe_gpu/bin/pw.x "$MODE" 1 1 "$FILE" "$@"
else
  bash "$HERE/dft_cells.sh" "$MODE" 8 2 "$FILE" "$@"
fi
