#!/bin/bash
# alpha is optional: relaxed at the 2x1x2 (k0) mesh only; its k1 relaxation would take 4-8 h on 3 ranks.
# Stop the automatic relax1 the moment it starts, then evaluate alpha at the k1 mesh on the relax0 geometry.
until grep -q "run alpha/relax1" ~/polar_margin/drive_alpha.log 2>/dev/null; do sleep 5; done
sleep 3
pkill -f "[p]rterun --bind-to none -np 3 pw.x -nk 3 -in relax1.in"
pkill -f "[p]w.x -nk 3 -in relax1.in"
until grep -q "W2 worker done" ~/polar_margin/worker_W2.log; do sleep 5; done
rm -rf ~/polar_margin/runs/alpha/relax1 ~/polar_margin/scratch/alpha_relax1
echo "[$(date -Is)] alpha relax1 stopped by design (reduced protocol); running scf_k1 scf_k1_3body on the relax0 geometry" | tee -a ~/polar_margin/drive_alpha.log
NP=3 NK=3 ~/polar_margin/bin/drive_cell.sh alpha scf_k1 scf_k1_3body
