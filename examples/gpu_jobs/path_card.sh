#!/usr/bin/env bash
# Paths between labelled cells in a potential on the card (examples/mlip_path.py): a kink sheet
# forming from the all-trans cell.  A job of examples/runpod_jobs_pod.sh:
#
#   path_card.sh FIT POINTS NAME:START:END [NAME:START:END...]
#
# FIT is a fit's name (the model the last refit.sh left) or a model file of the repository.
# START and END are frames of docs/predictions/mlip/cells14_frames_pbe_d3bj.extxyz, the
# all-trans cell first.  Each path is scanned from both ends at once (NAME_f from START, NAME_r
# from END), all scans side by side on the card, POINTS steps each.  Then, for each path, the
# scan that arrived where it was led and has the lower top above all-trans is taken: its top
# image goes to DFT first (the card's line if the DFT program was built for it, else the
# processors'), the image before it to the processors' line and the image after it behind the
# tops; and the same scan is refined between the images beside its top, and those two images
# are released (does one fall back and the other go on?).
#
# Into $JOB_OUT: NAME_f.json/.extxyz/.log and NAME_r..., NAME_x_fine..., NAME_x_free..., and
# chosen.txt.  PATH_NO_DFT=1 puts nothing in line; MLIP_DEVICE=cpu tries it without a card.
set -uo pipefail
FIT="$1"; POINTS="$2"; shift 2
HERE="$(dirname "$0")"
source "$HERE/mlip_env.sh"
export MLIP_DEVICE="${MLIP_DEVICE:-cuda}" OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
if [ -f "$FIT" ]; then MODEL="$FIT"; FIT="$(basename "$FIT" .model)"; else MODEL="$(ls "$(dirname "$JOB_OUT")"/*_refit/"$FIT.model" | tail -n 1)"; fi
[ -f "$MODEL" ] || { echo "no model for $FIT"; exit 1; }
LABELS=docs/predictions/mlip/cells14_frames_pbe_d3bj.extxyz
STEPS="${PATH_STEPS:-300}"
t0=$(date +%s)
for spec in "$@"; do
  IFS=: read -r name start end <<< "$spec"
  "$PY" examples/mlip_path.py scan --model "$MODEL" --potential "$FIT" --frames "$LABELS" --start "$start" --end "$end" \
    --name "${name}_f" --out "$JOB_OUT" --points "$POINTS" --steps "$STEPS" > "$JOB_OUT/${name}_f.log" 2>&1 &
  "$PY" examples/mlip_path.py scan --model "$MODEL" --potential "$FIT" --frames "$LABELS" --start "$end" --end "$start" \
    --name "${name}_r" --out "$JOB_OUT" --points "$POINTS" --steps "$STEPS" > "$JOB_OUT/${name}_r.log" 2>&1 &
done
wait
echo "scans: $(( $(date +%s) - t0 )) s on $MLIP_DEVICE" | tee "$JOB_OUT/seconds.txt"
tail -n 2 "$JOB_OUT"/*_f.log "$JOB_OUT"/*_r.log

GPUPW=/opt/qe_gpu/bin/pw.x
: > "$JOB_OUT/chosen.txt"
AFTER=()
for spec in "$@"; do
  IFS=: read -r name start end <<< "$spec"
  [ -f "$JOB_OUT/${name}_f.json" ] && [ -f "$JOB_OUT/${name}_r.json" ] || { echo "$name: a scan left no table, nothing chosen"; continue; }
  read -r scan top beside <<< "$("$PY" examples/mlip_path.py choose --scans "$JOB_OUT/${name}_f.json" "$JOB_OUT/${name}_r.json")"
  echo "$name $scan $top $beside" | tee -a "$JOB_OUT/chosen.txt"
  [ "${PATH_NO_DFT:-0}" = "1" ] && continue
  read -r before after <<< "$beside"
  file="${scan%.json}.extxyz"
  if [ -x "$GPUPW" ]; then
    bash "$HERE/dft_cells.sh" in=dft_card pw="$GPUPW" start 1 1 "$file" "$top"
    [ -n "$before" ] && bash "$HERE/dft_cells.sh" start 8 2 "$file" "$before"
    [ -n "$after" ] && AFTER+=("$file" "$after")
  else
    bash "$HERE/dft_cells.sh" start 8 2 "$file" "$top"
    [ -n "$before" ] && AFTER+=("$file" "$before")
    [ -n "$after" ] && AFTER+=("$file" "$after")
  fi
done
i=0
while [ $i -lt ${#AFTER[@]} ]; do   # the images after the tops, behind every top
  if [ -x "$GPUPW" ]; then bash "$HERE/dft_cells.sh" in=dft_card pw="$GPUPW" start 1 1 "${AFTER[$i]}" "${AFTER[$((i+1))]}"
  else bash "$HERE/dft_cells.sh" start 8 2 "${AFTER[$i]}" "${AFTER[$((i+1))]}"; fi
  i=$((i+2))
done

# the chosen scans in finer steps about their tops, and the images beside the tops released
t1=$(date +%s)
while read -r name scan top beside; do
  base="$(basename "$scan" .json)"
  "$PY" examples/mlip_path.py scan --model "$MODEL" --potential "$FIT" --frames "$LABELS" --refine "$scan" \
    --name "${base}_fine" --out "$JOB_OUT" --points "${PATH_FINE:-8}" --steps "$STEPS" > "$JOB_OUT/${base}_fine.log" 2>&1 &
  "$PY" examples/mlip_path.py release --model "$MODEL" --scan "$scan" --out "$JOB_OUT" > "$JOB_OUT/${base}_free.log" 2>&1 &
done < "$JOB_OUT/chosen.txt"
wait
echo "refinements and releases: $(( $(date +%s) - t1 )) s" | tee -a "$JOB_OUT/seconds.txt"
tail -n 2 "$JOB_OUT"/*_fine.log "$JOB_OUT"/*_free.log 2>/dev/null
exit 0
