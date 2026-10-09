#!/bin/bash
# usage: NP=3 NK=3 drive_cell.sh <cell> <stage|RELAX|FINALCHECK>...
#   RELAX      = relax0, relax1, then repeats relax1_r1, relax1_r2 while the last one fails the criteria
#   FINALCHECK = scf_k2; if it fails the force/stress criteria, relax2 (k2 mesh) then scf_k2 again
set -u
cell=$1; shift
NP=${NP:-3}; NK=${NK:-3}
PY=$HOME/.local/share/mamba/envs/pmpy/bin/python
S=/mnt/e/source/gh/polyfind-wt-margin/examples/polar_margin_dft.py
ROOT=$HOME/polar_margin/runs
LOG=$HOME/polar_margin/drive_${cell}.log
say() { echo "[$(date -Is)] $*" | tee -a $LOG; }
done_ok() { [ -f $ROOT/$cell/$1/$1.out ] && grep -q "JOB DONE" $ROOT/$cell/$1/$1.out; }
run() {
  local st=$1
  if done_ok $st; then say "skip $cell/$st (already done)"; return 0; fi
  $PY $S inputs --cell $cell --stage $st --root $ROOT >> $LOG 2>&1 || { say "inputs failed $cell/$st"; return 1; }
  say "run $cell/$st np=$NP nk=$NK"
  $HOME/polar_margin/bin/run_pw.sh $ROOT/$cell/$st $st.in $NK $NP; rc=$?
  rm -rf $HOME/polar_margin/scratch/${cell}_${st}
  $PY $S check --cell $cell --stage $st --root $ROOT >> $LOG 2>&1
  say "finished $cell/$st rc=$rc"
  return $rc
}
chk() { $PY $S check --cell $cell --stage $1 --root $ROOT > /dev/null 2>&1; }
for st in "$@"; do
  case $st in
    RELAX)
      run relax0 || exit 1
      run relax1 || exit 1
      prev=relax1
      for r in relax1_r1 relax1_r2; do
        if chk $prev; then break; fi
        run $r || exit 1; prev=$r
      done ;;
    FINALCHECK)
      run scf_k2 || exit 1
      if ! chk scf_k2; then
        say "scf_k2 fails the criteria: relax2 at k2"
        mv $ROOT/$cell/scf_k2 $ROOT/$cell/scf_k2_prerelax2
        run relax2 || exit 1
        chk relax2 || run relax2_r1 || exit 1
        run scf_k2 || exit 1
      fi ;;
    *) run $st || exit 1 ;;
  esac
done
say "pipeline $cell done: $*"
