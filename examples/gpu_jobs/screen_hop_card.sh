#!/usr/bin/env bash
# The hop of a row's cheapest kink, without a person reading the sheet table first.  A job of
# examples/runpod_jobs_pod.sh:
#
#   screen_hop_card.sh ROW NAME [POINTS]
#
# From the last sheets job's table (examples/gpu_jobs/screen_sheets_card.sh) it takes ROW's cheapest sheet that is
# at rest with its kink kept, and from the cells the screen relaxed on this pod the row's lowest all-trans cell at
# rest that has the sheet's atoms in the sheet's order (a scan needs that; in a cell of two different chains it is
# the all-trans cell with the chains the same way round), and runs examples/gpu_jobs/steps_card.sh on the pair: the
# scans NAME_a1, NAME_a2, NAME_b1, NAME_b2 land in this job's folder.
set -uo pipefail
ROW="$1"; NAME="$2"; POINTS="${3:-12}"
HERE="$(dirname "$0")"
source "$HERE/mlip_env.sh"
SHEETS="$( (ls "$(dirname "$JOB_OUT")"/*_screen_sheets_card/sheets_table.json 2>/dev/null || true) | sort | tail -n 1)"
[ -n "$SHEETS" ] || { echo "no sheets job on this pod has a table"; exit 1; }
RELAXED="$(dirname "$JOB_OUT")/shared/screen/relaxed"
PAIR="$("$PY" - "$SHEETS" "$ROW" "$RELAXED" <<'PY'
import glob, json, sys
from ase.io import read
row = json.load(open(sys.argv[1]))["rows"].get(sys.argv[2])
if not row:
    sys.exit(0)
held = {k: v for k, v in row["sites"].items() if v.get("relaxed") and v["kink_kept"] and v["converged"]}
if not held:
    sys.exit(0)
best = min(held, key=lambda k: held[k]["energy"])
sheet = read(glob.glob(f"{sys.argv[3]}/{held[best]['frame']}.start.xyz")[0])
flat = []
for f in sorted(glob.glob(f"{sys.argv[3]}/scr_{sys.argv[2]}_tt_*_in_*.start.xyz")):
    a = read(f)
    if "NOT converged" not in str(a.info.get("note", "")) and len(a) == len(sheet) and (a.numbers == sheet.numbers).all():
        flat.append((float(a.info["potential_eV_per_monomer"]), str(a.info["frame"])))
if not flat:
    sys.exit(0)
print(f"# {sys.argv[2]}: cheapest sheet at rest {held[best]['energy']:+.2f} per chain at {best} ({held[best]['units']}); "
      f"{len(flat)} all-trans cells at rest with its atoms in its order, the lowest {min(flat)[1]}", file=sys.stderr)
print(held[best]["frame"], min(flat)[1])
PY
)"
[ -n "$PAIR" ] || { echo "$ROW has no sheet at rest with its kink kept in $SHEETS, or no all-trans cell at rest with that sheet's atoms in $RELAXED"; exit 1; }
read -r START END <<< "$PAIR"
echo "table: $SHEETS"
echo "the hop of $START toward $END"
bash "$HERE/steps_card.sh" docs/predictions/mlip/models/all_s1.model "$POINTS" "$NAME:$START:$END"
