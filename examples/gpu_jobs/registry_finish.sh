#!/usr/bin/env bash
# The end of a registry job on a pod (examples/gpu_jobs/registry_card.sh, registry_more.sh): the table of
# the relaxed cells, the choice, and the chosen onto the two DFT lines.
#
#   registry_finish.sh STARTS RELAXED
#
# Into $JOB_OUT: registry_table.json, registry_relaxed.extxyz, chosen.txt (the card's line, then the
# processors').  A cell is chosen when it is the lowest of its composition and conformation, kept its
# conformation, and is at rest within REGISTRY_WITHIN (0.3) kcal/mol per chain of its labelled cell, or
# not at rest but already REGISTRY_UNSETTLED (-0.3) below it; the film composition's all-trans only below
# -1.6 (3.23's cell, 1.27 below, is labelled).  Cells in a DFT line already are not chosen again.
# REGISTRY_NO_DFT=1 puts nothing in line.
set -uo pipefail
STARTS="$1"; RELAXED="$2"
HERE="$(dirname "$0")"
LABELS=docs/predictions/mlip/cells14_frames_pbe_d3bj.extxyz
GPUPW=/opt/qe_gpu/bin/pw.x
SHARED="$(dirname "$JOB_OUT")/shared"
"$PY" examples/registry_cells.py table --starts "$STARTS" --relaxed "$RELAXED" --labels "$LABELS" \
  --out "$JOB_OUT/registry_table.json" --frames-out "$JOB_OUT/registry_relaxed.extxyz" || exit 1
INLINE="$(cat "$SHARED"/dft/line.txt "$SHARED"/dft_card/line.txt 2>/dev/null | awk '{print $3}' | tr '\n' ' ')"
"$PY" examples/registry_cells.py choose --table "$JOB_OUT/registry_table.json" --within "${REGISTRY_WITHIN:-0.3}" --below film:tt=-1.6 \
  --unsettled-below "${REGISTRY_UNSETTLED:--0.3}" --skip $INLINE > "$JOB_OUT/chosen.txt" || exit 1
echo "chosen (the card's line, then the processors'):"; cat "$JOB_OUT/chosen.txt"
[ "${REGISTRY_NO_DFT:-0}" = "1" ] && exit 0
CARD="$(sed -n 1p "$JOB_OUT/chosen.txt")"; REST="$(sed -n 2p "$JOB_OUT/chosen.txt")"
if [ -x "$GPUPW" ]; then
  [ -n "$CARD" ] && bash "$HERE/dft_cells.sh" in=dft_card pw="$GPUPW" start 1 1 "$JOB_OUT/registry_relaxed.extxyz" $CARD
  [ -n "$REST" ] && bash "$HERE/dft_cells.sh" start 8 2 "$JOB_OUT/registry_relaxed.extxyz" $REST
else
  [ -n "$CARD$REST" ] && bash "$HERE/dft_cells.sh" start 8 2 "$JOB_OUT/registry_relaxed.extxyz" $CARD $REST
fi
exit 0
