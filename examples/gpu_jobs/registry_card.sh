#!/usr/bin/env bash
# Other registries of labelled cells relaxed in the potential on the card, and the lowest of each
# composition and conformation put to DFT (examples/registry_cells.py; docs/MLIP_SECOND_OPINION.md
# 3.24).  A job of examples/runpod_jobs_pod.sh:
#
#   registry_card.sh STARTS SIDE STEPS BUDGET [FILE NAME...]
#
# STARTS is a file of starts in the repository (examples/registry_starts.sh), in the order they
# matter.  They are dealt in turn to SIDE processes side by side on the card, each relaxing its
# cells one after another to 0.004 eV/A or STEPS steps; after BUDGET seconds whatever is still
# running is stopped (a cell is written every 100 steps, marked NOT converged).  Meanwhile the
# processors build the DFT program for the card (examples/gpu_jobs/qe_gpu_build.sh) and, when
# it is built, start labelling the named frames of FILE on their own line (cells known before
# this job, so that the processors do not wait for it).
#
# Then: registry_table.json (each cell against the labelled cell it was made from),
# registry_relaxed.extxyz (the relaxed cells), chosen.txt (examples/registry_cells.py choose: the
# card's line, then the processors'), and the chosen go to the two DFT lines.
#
# REGISTRY_NO_DFT=1 puts nothing in line and builds nothing; MLIP_DEVICE=cpu tries it without a card.
set -uo pipefail
STARTS="$1"; SIDE="$2"; STEPS="$3"; BUDGET="$4"; shift 4
HERE="$(dirname "$0")"
source "$HERE/mlip_env.sh"
export MLIP_DEVICE="${MLIP_DEVICE:-cuda}" OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
MODEL=docs/predictions/mlip/models/all_s1.model
LABELS=docs/predictions/mlip/cells14_frames_pbe_d3bj.extxyz
GPUPW=/opt/qe_gpu/bin/pw.x
t0=$(date +%s)

if [ "${REGISTRY_NO_DFT:-0}" != "1" ]; then
  (
    if [ ! -x "$GPUPW" ]; then
      bash "$HERE/qe_gpu_build.sh" "${REGISTRY_BUILD_JOBS:-4}" > "$JOB_OUT/qe_gpu_build.log" 2>&1 \
        && echo "the DFT program for the card is built" || echo "THE BUILD FOR THE CARD FAILED (qe_gpu_build.log): DFT is the processors' alone"
    fi
    if [ $# -ge 2 ]; then file="$1"; shift; bash "$HERE/dft_cells.sh" start 8 2 "$file" "$@"; fi
  ) > "$JOB_OUT/processors.log" 2>&1 &
  PROCESSORS=$!
fi

mapfile -t PICKS < <("$PY" - "$STARTS" "$SIDE" <<'PY'
import sys
from ase.io import read
names = [str(a.info["frame"]) for a in read(sys.argv[1], index=":")]
side = int(sys.argv[2])
for k in range(side):
    print("(" + "|".join(names[k::side]) + ")")
PY
)
PIDS=()
for k in "${!PICKS[@]}"; do
  timeout "$BUDGET" "$PY" examples/mlip_train.py relax --model "$MODEL" --name all_s1 --cells "$STARTS" --pick "${PICKS[$k]}" \
    --out "$JOB_OUT/relaxed" --threads 1 --fmax 0.004 --steps "$STEPS" --every 100 --skip-done > "$JOB_OUT/relax_$k.log" 2>&1 &
  PIDS+=($!)
done
for pid in "${PIDS[@]}"; do wait "$pid" || true; done
echo "relaxed for $(( $(date +%s) - t0 )) s on $MLIP_DEVICE, $SIDE side by side" | tee "$JOB_OUT/seconds.txt"
grep -h "steps, " "$JOB_OUT"/relax_*.log || true

"$PY" examples/registry_cells.py table --starts "$STARTS" --relaxed "$JOB_OUT/relaxed" --labels "$LABELS" \
  --out "$JOB_OUT/registry_table.json" --frames-out "$JOB_OUT/registry_relaxed.extxyz" || exit 1
# the film composition's low all-trans cell is labelled already (3.23: 1.27 below its upright cell here): only a lower one is news
"$PY" examples/registry_cells.py choose --table "$JOB_OUT/registry_table.json" --within "${REGISTRY_WITHIN:-0.3}" --below film:tt=-1.0 > "$JOB_OUT/chosen.txt" || exit 1
echo "chosen (the card's line, then the processors'):"; cat "$JOB_OUT/chosen.txt"
[ "${REGISTRY_NO_DFT:-0}" = "1" ] && exit 0

wait "$PROCESSORS" || true
cat "$JOB_OUT/processors.log" | tail -n 5
CARD="$(sed -n 1p "$JOB_OUT/chosen.txt")"; REST="$(sed -n 2p "$JOB_OUT/chosen.txt")"
if [ -x "$GPUPW" ]; then
  [ -n "$CARD" ] && bash "$HERE/dft_cells.sh" in=dft_card pw="$GPUPW" start 1 1 "$JOB_OUT/registry_relaxed.extxyz" $CARD
  [ -n "$REST" ] && bash "$HERE/dft_cells.sh" start 8 2 "$JOB_OUT/registry_relaxed.extxyz" $REST
else
  [ -n "$CARD$REST" ] && bash "$HERE/dft_cells.sh" start 8 2 "$JOB_OUT/registry_relaxed.extxyz" $CARD $REST
fi
exit 0
