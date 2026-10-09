#!/bin/bash
# usage: worker.sh <name> <log-to-wait-on> <pattern>  -- starts after the pattern appears, then pulls jobs
name=$1
until grep -q "$3" "$2" 2>/dev/null; do sleep 60; done
while true; do
  job=$(flock ~/polar_margin/queue.lock bash -c "l=\$(head -1 ~/polar_margin/queue.txt); [ -n \"\$l\" ] && sed -i 1d ~/polar_margin/queue.txt; echo \"\$l\"")
  [ -z "$job" ] && break
  cell=${job%% *}
  echo "[$(date -Is)] $name takes: $job"
  if [ "$cell" != alpha ]; then
    until grep -q "pipeline $cell done: RELAX FINALCHECK" ~/polar_margin/drive_$cell.log 2>/dev/null; do sleep 60; done
  fi
  NP=3 NK=3 ~/polar_margin/bin/drive_cell.sh $job || echo "[$(date -Is)] $name: job FAILED: $job"
done
echo "[$(date -Is)] $name worker done"
