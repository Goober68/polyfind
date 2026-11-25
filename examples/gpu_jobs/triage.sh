#!/usr/bin/env bash
# One subcommand of examples/chain_triage.py as a job of examples/mi100_runner.sh:
#
#   triage.sh SUBCOMMAND [its arguments but --frames]
#
# as, in a queue,
#
#   triage.sh kink --names urea1_tt_x1 --device cuda --long 6 --turned 2,3,4
#   triage.sh points --names urea1_tt_x1_kink4 --relax
#
# The frames folder is the machine's own (TRIAGE_FRAMES), first given the relaxed cells, records and
# paths the repository keeps (the `seed` subcommand: nothing there is overwritten).  The job's output
# is its log; the records and paths it made or changed (small files) are copied to $JOB_OUT, and the
# DFT program's own outputs stay in the frames folder but for their last lines.
set -u
# shellcheck disable=SC1091
. "${MI100_ROOT:-$HOME/polyfind_mi100}/env.sh"
[ $# -ge 1 ] || { echo "no subcommand"; exit 2; }
mkdir -p "$TRIAGE_FRAMES" "$JOB_OUT"
"$PY" examples/chain_triage.py seed --frames "$TRIAGE_FRAMES" 2>&1 | tail -n 1 | cut -c1-300
touch "$JOB_OUT/.started"
SUB="$1"; shift
"$PY" examples/chain_triage.py "$SUB" --frames "$TRIAGE_FRAMES" "$@" 2>&1 | grep -v -i "warn\|torch.load\|_Jd, _W3j\|cuequivariance"
code="${PIPESTATUS[0]}"
find "$TRIAGE_FRAMES" -maxdepth 1 -type f \( -name '*.json' -o -name '*.path.xyz' -o -name '*.relaxed.xyz' \) -newer "$JOB_OUT/.started" -size -5M -exec cp {} "$JOB_OUT/" \;
# where each DFT run of this job had got to: its energies and its last lines
find "$TRIAGE_FRAMES/work" -type f -name '*.out' -newer "$JOB_OUT/.started" 2>/dev/null | sort | while read -r out; do
  echo "-- ${out#"$TRIAGE_FRAMES"/work/}: $(grep -c '^!' "$out") energies, last $(grep '^!' "$out" | tail -n 1 | awk '{print $(NF-1)}') Ry, $(grep -E 'PWSCF .*WALL' "$out" | tail -n 1 | sed 's/^ *//')"
done
rm -f "$JOB_OUT/.started"
exit "$code"
