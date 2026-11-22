#!/usr/bin/env bash
# Cells of the registry search relaxed on this pod (examples/gpu_jobs/registry_card.sh) put in line
# for a DFT single point each, by name: the relaxed cells are under the pod's served directory,
# which a queue line cannot name.  A job of examples/runpod_jobs_pod.sh:
#
#   registry_dft.sh card|processors start|next NAME [NAME...]
#
# NAME are frames of the last registry_card.sh job's registry_relaxed.extxyz
# (reg_matrix_tt_p1m1_tp117_in_all_s1).  `card` is the line of the DFT program built for the card
# (one process), `processors` the other (eight processes in two pools); start and next are
# examples/gpu_jobs/dft_cells.sh's: the end of the line, or its front.
set -euo pipefail
LINE="$1"; MODE="$2"; shift 2
HERE="$(dirname "$0")"
FILE="$( (ls "$(dirname "$JOB_OUT")"/*_registry_card/registry_relaxed.extxyz 2>/dev/null || true) | tail -n 1)"
[ -f "$FILE" ] || { echo "no relaxed registry cells on this pod"; exit 1; }
if [ "$LINE" = "card" ]; then
  bash "$HERE/dft_cells.sh" in=dft_card pw=/opt/qe_gpu/bin/pw.x "$MODE" 1 1 "$FILE" "$@"
else
  bash "$HERE/dft_cells.sh" "$MODE" 8 2 "$FILE" "$@"
fi
