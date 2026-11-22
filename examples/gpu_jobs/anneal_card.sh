#!/usr/bin/env bash
# Labelled cells shaken at a temperature in the potential on the card and relaxed again
# (examples/mlip_anneal.py; docs/MLIP_SECOND_OPINION.md 3.24).  A job of examples/runpod_jobs_pod.sh:
#
#   anneal_card.sh STARTS SIDE MD_STEPS KELVIN RELAX_STEPS BUDGET
#
# STARTS is a file of starts in the repository whose frames carry registry_* entries
# (examples/registry_cells.py derive).  They are dealt in turn to SIDE processes side by side; each cell
# is run MD_STEPS femtoseconds at KELVIN in its own cell and relaxed (0.004 eV/A, RELAX_STEPS steps at
# most); after BUDGET seconds whatever is still running is stopped.  Then examples/gpu_jobs/
# registry_finish.sh: the table against the cells they were made from, and a cell that came to rest more
# than 0.3 kcal/mol per chain BELOW the cell it started as goes to DFT.
set -uo pipefail
STARTS="$1"; SIDE="$2"; MD="$3"; KELVIN="$4"; STEPS="$5"; BUDGET="$6"
HERE="$(dirname "$0")"
source "$HERE/mlip_env.sh"
export MLIP_DEVICE="${MLIP_DEVICE:-cuda}" OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
MODEL=docs/predictions/mlip/models/all_s1.model
t0=$(date +%s)
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
  timeout "$BUDGET" "$PY" examples/mlip_anneal.py --model "$MODEL" --name all_s1 --cells "$STARTS" --pick "${PICKS[$k]}" --out "$JOB_OUT/relaxed" \
    --kelvin "$KELVIN" --md-steps "$MD" --steps "$STEPS" --seed "$((k + 1))" > "$JOB_OUT/anneal_$k.log" 2>&1 &
  PIDS+=($!)
done
for pid in "${PIDS[@]}"; do wait "$pid" || true; done
echo "annealed and relaxed for $(( $(date +%s) - t0 )) s on $MLIP_DEVICE, $SIDE side by side" | tee "$JOB_OUT/seconds.txt"
grep -h "steps at" "$JOB_OUT"/anneal_*.log || true
export REGISTRY_WITHIN=-0.3
bash "$HERE/registry_finish.sh" "$STARTS" "$JOB_OUT/relaxed"
