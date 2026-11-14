#!/usr/bin/env bash
# One row of the single-chain step from nothing (examples/kink_migration.py): the end states I and
# B from the kinked crystal with the moving chain's torsions reset, the two bands, and the two
# saddles settled from the bands' tops; all minimisations on the closure surface.  A job of
# examples/runpod_jobs_pod.sh.  Arguments:
#
#   band_row.sh ROW NX NY [from CELL] [tiled]
#
# ROW is a sequence row with a stored kinked crystal ("vdf9-trfe4-cfe1-14"), NX NY the cell.
# "from CELL": start each band from the band of the same row in that cell (an earlier job's, in
# the shared table).  "tiled": take the tiled two-chain cell as the sheet without the check that
# relaxes the perturbed supercell with every chain and the cell free.  The row goes to the shared
# table (jobs/shared/steps.json), which starts from the stored table's kinked crystals and holds
# none of its rows: a row that is stored already is computed again.  A band gets 300 iterations
# (STEP_BAND_MAXITER) and need not converge: the saddle is settled from its top afterwards, and in
# the four-chain PVDF cell forty iterations (tops 1.1 and 3.8 kcal/mol too high, forces of 4 and
# 6) settle to the stored saddles in the fourth decimal.
set -euo pipefail
ROW="$1"; NX="$2"; NY="$3"; shift 3
FROM=""; TILED=()
while [ $# -gt 0 ]; do
  case "$1" in
    from) FROM="$2"; shift 2 ;;
    tiled) TILED=(--trust-tiled); shift ;;
    *) echo "band_row.sh: unknown argument $1"; exit 2 ;;
  esac
done
KINK="${STEP_KINK:-0}"
SHARED="$(dirname "$JOB_OUT")/shared"
mkdir -p "$SHARED/ck"
TABLE="$SHARED/steps.json"
key() { echo "$ROW@$KINK $1 mover 0 ahead host flexible"; }
"$PY" examples/kink_wide_cell.py table --out "$TABLE"
COMMON=(--row "$ROW" --kink-at "$KINK" --nx "$NX" --ny "$NY" --mover 0 --host flexible --surface --workers 1
        --patience "${STEP_PATIENCE:-40}" --surface-gtol "${STEP_GTOL:-3e-4}" --end-maxiter "${STEP_END_MAXITER:-600}"
        --checkpoints "$SHARED/ck" --out "$TABLE")
BAND=()
[ -n "$FROM" ] && BAND=(--band-from "$(key "$FROM")")
echo "== $(key "${NX}x${NY}"): end states and bands${FROM:+, each band started from that of the $FROM cell}"
"$PY" examples/kink_migration.py "${COMMON[@]}" "${TILED[@]}" --maxiter "${STEP_BAND_MAXITER:-300}" --fmax "${STEP_FMAX:-0.1}" "${BAND[@]}"
echo "== $(key "${NX}x${NY}"): the saddles, settled from the tops of the bands"
"$PY" examples/kink_migration.py "${COMMON[@]}" --trust-tiled --polish
"$PY" - "$TABLE" "$(key "${NX}x${NY}")" "$JOB_OUT/result.json" <<'PY'
import json, sys
table, key, out = sys.argv[1:4]
row = json.load(open(table))["rows"][key]
drop = ("state", "images", "rounds", "torsions_deg", "saddle_torsions_deg", "state_relaxed")
def light(v):
    if isinstance(v, dict):
        return {k: light(w) for k, w in v.items() if k not in drop}
    return v
json.dump(light(row), open(out, "w"), indent=1)
kc = row["kinked_crystal"]
check = kc.get("supercell_check")
print(f"   {key}: {row['chains']} chains")
if check:
    print(f"   sheet check: the perturbed supercell relaxed to {check['relaxed_above_tiled_per_chain']:+.2e} kcal/mol per chain against the tiled cell; converged {check['converged']}")
for name in ("I", "B"):
    st = row[name]
    print(f"   {name}: {st['energy']:+.4f}; {st['evaluations']} evaluations, {st['seconds']:.0f} s; gradient {st['gradient_max']:.1e}; converged {st['converged']}; {st['torsions']}")
for seg in ("A-I", "I-B"):
    b = row.get(seg)
    if not b:
        continue
    p = b.get("polished", {})
    print(f"   {seg}: band top {b['saddle']:+.4f} at force {b['force_max']:.3f} after {b['iterations']} iterations ({b['seconds']:.0f} s, converged {b['converged']})"
          + (f"; settled saddle {p['saddle']:+.4f}, along {p['force_along']:+.4f}, across {p['gradient_across']:.1e}, {p['evaluations']} evaluations, converged {p['converged']}" if p else ""))
if "barrier_forward_polished" in row:
    print(f"   barrier on the settled saddles: {row['barrier_forward_polished']:.4f} ahead, {row['barrier_back_polished']:.4f} back")
elif "barrier_forward" in row:
    print(f"   barrier on the bands' tops: {row['barrier_forward']:.4f} ahead, {row['barrier_back']:.4f} back")
PY
