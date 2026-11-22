#!/usr/bin/env bash
# The screen of docs/MLIP_SECOND_OPINION.md 3.25 on a pod's card: the starts of examples/screen14_starts.sh relaxed in
# the potential, each row's margin from its lowest all-trans and helix cells, and the rows whose margin is inside a window
# put to DFT.  A job of examples/runpod_jobs_pod.sh:
#
#   screen_card.sh STARTS SIDE STEPS BUDGET [MOST]
#
# STARTS is a file of starts in the repository (six cells a row, row after row).  They are dealt in turn to SIDE processes
# side by side on the card, each relaxing its cells one after another to 0.004 eV/A or STEPS steps; after BUDGET seconds
# whatever is still running is stopped (a cell is written every 100 steps, marked NOT converged).  The relaxed cells of
# every screen job of the pod are kept in one folder (<served>/jobs/shared/screen/relaxed), so a later job with the same
# STARTS goes on where the earlier one stopped.  Meanwhile, on the first job, the processors build the DFT program for the
# card (examples/gpu_jobs/qe_gpu_build.sh).
#
# Then into $JOB_OUT: screen_table.json (examples/screen14.py table), screen_relaxed.extxyz, chosen.txt.  The rows whose
# margin in the potential is between SCREEN_LOW and SCREEN_HIGH (-0.06 and +0.14 kcal/mol per monomer: the potential has
# the film's composition 0.01 above DFT and the matrix 0.06 below, and the margins wanted are DFT's 0 to +0.12), MOST of
# them at most (default 2), the ones with every cell at rest first, have their lowest all-trans cell, lowest helix cell
# and other all-trans cell put on the two DFT lines (the card's line, then the processors').  Rows in SCREEN_SKIP (the
# controls: they have DFT) and cells in a DFT line already are left out.  SCREEN_NO_DFT=1 puts nothing in line.
set -uo pipefail
STARTS="$1"; SIDE="$2"; STEPS="$3"; BUDGET="$4"; MOST="${5:-2}"
HERE="$(dirname "$0")"
source "$HERE/mlip_env.sh"
export MLIP_DEVICE="${MLIP_DEVICE:-cuda}" OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
MODEL=docs/predictions/mlip/models/all_s1.model
GPUPW=/opt/qe_gpu/bin/pw.x
SHARED="$(dirname "$JOB_OUT")/shared"
RELAXED="$SHARED/screen/relaxed"
mkdir -p "$RELAXED"
t0=$(date +%s)

if [ "${SCREEN_NO_DFT:-0}" != "1" ] && [ ! -x "$GPUPW" ]; then
  ( bash "$HERE/qe_gpu_build.sh" "${REGISTRY_BUILD_JOBS:-4}" > "$JOB_OUT/qe_gpu_build.log" 2>&1 \
      && echo "the DFT program for the card is built" || echo "THE BUILD FOR THE CARD FAILED (qe_gpu_build.log): DFT is the processors' alone" ) > "$JOB_OUT/processors.log" 2>&1 &
  BUILD=$!
fi

mapfile -t PICKS < <("$PY" - "$STARTS" "$SIDE" <<'PY'
import sys
from ase.io import read
names = [str(a.info["frame"]) for a in read(sys.argv[1], index=":")]
side = int(sys.argv[2])
for k in range(side):
    if names[k::side]:
        print("(" + "|".join(names[k::side]) + ")")
PY
)
PIDS=()
for k in "${!PICKS[@]}"; do
  timeout "$BUDGET" "$PY" examples/mlip_train.py relax --model "$MODEL" --name all_s1 --cells "$STARTS" --pick "${PICKS[$k]}" \
    --out "$RELAXED" --threads 1 --fmax 0.004 --steps "$STEPS" --every 100 --skip-done > "$JOB_OUT/relax_$k.log" 2>&1 &
  PIDS+=($!)
done
for pid in "${PIDS[@]}"; do wait "$pid" || true; done
echo "relaxed for $(( $(date +%s) - t0 )) s on $MLIP_DEVICE, $SIDE side by side" | tee "$JOB_OUT/seconds.txt"
grep -h "steps, " "$JOB_OUT"/relax_*.log || true
[ -n "${BUILD:-}" ] && { wait "$BUILD" || true; tail -n 2 "$JOB_OUT/processors.log"; }

"$PY" examples/screen14.py table --starts "$STARTS" --relaxed "$RELAXED" --out "$JOB_OUT/screen_table.json" --frames-out "$JOB_OUT/screen_relaxed.extxyz" || exit 1
"$PY" examples/screen14.py choose --table "$JOB_OUT/screen_table.json" --window "${SCREEN_LOW:--0.06}" "${SCREEN_HIGH:-0.14}" --most "$MOST" \
  --skip ${SCREEN_SKIP:-vdf9-trfe4-cfe1-14 vdf10-trfe4-14 vdf9-trfe3-an1-rvdf1-14} > "$JOB_OUT/chosen_all.txt" || exit 1
# cells in a DFT line already are not put there again
INLINE="$(cat "$SHARED"/dft/line.txt "$SHARED"/dft_card/line.txt 2>/dev/null | awk '{print $3}' | tr '\n' ' ')"
for n in 1 2; do
  LINE=""
  for name in $(sed -n "${n}p" "$JOB_OUT/chosen_all.txt"); do
    case " $INLINE " in *" $name "*) ;; *) LINE="$LINE $name" ;; esac
  done
  echo "$LINE" | sed 's/^ //' >> "$JOB_OUT/chosen.txt"
done
echo "chosen (the card's line, then the processors'):"; cat "$JOB_OUT/chosen.txt"
[ "${SCREEN_NO_DFT:-0}" = "1" ] && exit 0
CARD="$(sed -n 1p "$JOB_OUT/chosen.txt")"; REST="$(sed -n 2p "$JOB_OUT/chosen.txt")"
if [ -x "$GPUPW" ]; then
  [ -n "$CARD" ] && bash "$HERE/dft_cells.sh" in=dft_card pw="$GPUPW" start 1 1 "$JOB_OUT/screen_relaxed.extxyz" $CARD
  [ -n "$REST" ] && bash "$HERE/dft_cells.sh" start 8 2 "$JOB_OUT/screen_relaxed.extxyz" $REST
else
  [ -n "$CARD$REST" ] && bash "$HERE/dft_cells.sh" start 8 2 "$JOB_OUT/screen_relaxed.extxyz" $CARD $REST
fi
exit 0
