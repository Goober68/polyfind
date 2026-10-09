#!/bin/bash
# usage: run_pw.sh <workdir> <input> [nk] [np]
# Runs pw.x niced, 1 OpenMP thread, at most 6 MPI ranks; writes <input%.in>.out next to the input.
set -u
dir=$1; inp=$2; nk=${3:-6}; np=${4:-6}
export PATH=$HOME/.local/share/mamba/envs/qe/bin:$PATH
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd "$dir" || exit 2
out=${inp%.in}.out
echo "[$(date -Is)] start $dir/$inp nk=$nk np=$np" >> $HOME/polar_margin/runlog.txt
nice -n 10 mpirun --bind-to none -np $np pw.x -nk $nk -in "$inp" > "$out" 2> "${inp%.in}.err"
rc=$?
echo "[$(date -Is)] end   $dir/$inp rc=$rc $(grep -c "JOB DONE" "$out")" >> $HOME/polar_margin/runlog.txt
exit $rc
