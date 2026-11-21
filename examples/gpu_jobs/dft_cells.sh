#!/usr/bin/env bash
# DFT single points (PBE-D3(BJ), the anchor's set-up: examples/mlip_frames.py label) on stored
# cells, on the pod's processors, while its card does the jobs that follow in the queue.
# A job of examples/runpod_jobs_pod.sh:
#
#   dft_cells.sh start NP POOLS FILE NAME [NAME...]
#       put the named frames of FILE (a file of frames in the repository) in line, start the
#       worker if it is not running, and return: the queue goes on
#   dft_cells.sh next NP POOLS FILE NAME [NAME...]
#       the same, but the named frames go to the front of the line: they are the next cells run
#   dft_cells.sh last NAME [NAME...]
#       move cells already in line to its end
#   dft_cells.sh usepw PW NP POOLS
#       every cell still waiting in line is run with the pw.x at PW on NP processes in POOLS pools
#       (a build for the card, examples/gpu_jobs/qe_gpu_build.sh: one process, one pool)
#   dft_cells.sh drop NAME [NAME...]
#       take cells that are still waiting out of the line (another line has them)
#   dft_cells.sh wait
#       return when the line is empty (so that the pod does not count itself idle and leave)
#
# Before the mode, `in=NAME` names another line (its own folder <served>/jobs/shared/NAME/, its
# own worker and log) and `pw=PATH` the pw.x that line's cells are run with: a second line for
# the card beside the processors' one, as
#   dft_cells.sh in=dft_card pw=/opt/qe_gpu/bin/pw.x start 1 1 FILE NAME [NAME...]
#
# The first start installs Quantum ESPRESSO 7.5 (conda-forge, as on the machine the frames were
# made on) and fetches the five SSSP 1.3.0 PBE precision files, which must have the checksums of
# the ones every label so far was computed with.  Before the first label the anchor's own
# beta-PVDF input is run and must return the anchor's energy (to DFT_CHECK_TOL Ry, default 1e-6:
# another processor moves the ninth decimal; the difference is printed).
#
# Everything is under <served>/jobs/shared/dft/: log.txt (one line a cell: energy, largest force,
# seconds), NAME.dft.xyz (the labelled frame), NAME.pw.out, and work/NAME/NAME.out while a cell
# runs.  The cells are run one after another, each on NP processes in POOLS k-point pools.
# DFT_QE and DFT_PSEUDO name an installation and the files already on the machine; DFT_DRY=1
# writes the worker and does not start it (to try this script where nothing may be left running).
set -euo pipefail
SUBDIR="dft"; OWN_PW=""
while [ $# -gt 0 ]; do
  case "$1" in
    in=*) SUBDIR="${1#in=}"; shift ;;
    pw=*) OWN_PW="${1#pw=}"; shift ;;
    *) break ;;
  esac
done
MODE="$1"; shift
SHARED="$(dirname "$JOB_OUT")/shared"
DFT="$SHARED/$SUBDIR"
QE="${DFT_QE:-/opt/qe}"
PSEUDO="${DFT_PSEUDO:-/opt/pseudo}"
mkdir -p "$DFT"

running() { [ -f "$DFT/worker.pid" ] && kill -0 "$(cat "$DFT/worker.pid")" 2>/dev/null; }

if [ "$MODE" = "wait" ]; then
  while running; do sleep 30; done
  cat "$DFT/log.txt"
  exit 0
fi
if [ "$MODE" = "last" ]; then
  for name in "$@"; do
    awk -v n="$name" '$3 != n' "$DFT/line.txt" > "$DFT/line.new"
    awk -v n="$name" '$3 == n' "$DFT/line.txt" >> "$DFT/line.new"
    mv "$DFT/line.new" "$DFT/line.txt"
  done
  echo "the line is now:"; cat "$DFT/line.txt"
  exit 0
fi
if [ "$MODE" = "drop" ]; then
  for name in "$@"; do
    if [ -f "$DFT/$name.dft.xyz" ] || grep -qx "$name" "$DFT/tried.txt" 2>/dev/null; then echo "$name: already started here, left in this line"
    else echo "$name" >> "$DFT/tried.txt"; echo "$name: out of this line"; fi
  done
  exit 0
fi
if [ "$MODE" = "usepw" ]; then
  PW="$1"; NP="$2"; POOLS="$3"
  [ -x "$PW" ] || { echo "no pw.x at $PW: the line is left as it is"; exit 1; }
  : > "$DFT/line.new"
  while read -r np pools name pw; do
    if [ -f "$DFT/$name.dft.xyz" ] || grep -qx "$name" "$DFT/tried.txt" 2>/dev/null; then echo "$np $pools $name $pw" >> "$DFT/line.new"
    else echo "$NP $POOLS $name $PW" >> "$DFT/line.new"; fi
  done < "$DFT/line.txt"
  mv "$DFT/line.new" "$DFT/line.txt"
  echo "the line is now:"; cat "$DFT/line.txt"
  exit 0
fi
[ "$MODE" = "start" ] || [ "$MODE" = "next" ] || { echo "dft_cells.sh start|next NP POOLS FILE NAME... | last NAME... | usepw PW NP POOLS | wait"; exit 2; }
NP="$1"; POOLS="$2"; FILE="$3"; shift 3
source "$(dirname "$0")/mlip_env.sh"

echo "processors shown $(nproc), allotted $(cat /sys/fs/cgroup/cpu.max 2>/dev/null || echo unknown); memory:"
free -g | head -2

if [ ! -x "$QE/bin/pw.x" ]; then
  mkdir -p /opt/mm
  curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xj -C /opt/mm bin/micromamba
  /opt/mm/bin/micromamba create -y -q -p "$QE" -c conda-forge "qe=7.5"
fi
[ -x "$QE/bin/pw.x" ] || { echo "no pw.x at $QE/bin: nothing is labelled"; exit 1; }
DEFAULT_PW="${OWN_PW:-$QE/bin/pw.x}"
[ -x "$DEFAULT_PW" ] || { echo "no pw.x at $DEFAULT_PW: nothing is labelled"; exit 1; }

# name and sha256 of the files in use since 2026-11-16 (the first three are also in deliverables/polar_margin_dft/raw_stages.json)
SUMS="ba6ca6af365222805e3ce85564616b50a2feb562106dfbf1d725170cc7afd7b8  F.oncvpsp.upf
69d631780c0178c5239ed3dcb4a3e22c13bd13cde9beb97e123d036949da5c09  H_ONCV_PBE-1.0.oncvpsp.upf
36f1fc9a24ff5c788a536d4c1a53f82adca7b1fd132827fc654c124d69f274ca  N.oncvpsp.upf
9900d1efd50b9848e31849f39094b33348486b400ee51e0f3922f716137cf3d7  C.pbe-n-kjpaw_psl.1.0.0.UPF
a673d3adf336817175b7ffcbce2f811dd5929176c94a2c829cd1fc578ca924b8  Cl.pbe-n-rrkjus_psl.1.0.0.UPF"
if ! (cd "$PSEUDO" 2>/dev/null && echo "$SUMS" | sha256sum -c --quiet 2>/dev/null); then
  mkdir -p "$PSEUDO" /opt/pseudo_fetch
  curl -sL "https://archive.materialscloud.org/api/records/rcyfm-68h65/files/SSSP_1.3.0_PBE_precision.tar.gz/content" -o /opt/pseudo_fetch/sssp.tar.gz
  tar -xzf /opt/pseudo_fetch/sssp.tar.gz -C /opt/pseudo_fetch
  for f in F.oncvpsp.upf H_ONCV_PBE-1.0.oncvpsp.upf N.oncvpsp.upf C.pbe-n-kjpaw_psl.1.0.0.UPF Cl.pbe-n-rrkjus_psl.1.0.0.UPF; do
    cp "$(find /opt/pseudo_fetch -name "$f" | head -n 1)" "$PSEUDO/"
  done
  (cd "$PSEUDO" && echo "$SUMS" | sha256sum -c) || { echo "the pseudopotential files are not the ones the labels were computed with: nothing is labelled"; exit 1; }
fi
echo "pseudopotentials: the five files, checksums as stored"

# the named frames as start files; a name the file does not have stops the job before anything runs
"$PY" - "$FILE" "$DFT" "$@" <<'PY'
import sys
from pathlib import Path
from ase.io import read, write

frames = {str(a.info["frame"]): a for a in read(sys.argv[1], index=":")}
missing = [n for n in sys.argv[3:] if n not in frames]
if missing:
    sys.exit(f"not in {sys.argv[1]}: {missing}")
for n in sys.argv[3:]:
    write(Path(sys.argv[2]) / f"{n}.start.xyz", frames[n], format="extxyz")
    print(f"{n}: {len(frames[n])} atoms, k-points {frames[n].info['kmesh']}")
PY
if [ "$MODE" = "next" ]; then   # in front of everything in line: the worker takes the first cell it has not tried
  : > "$DFT/line.new"
  for name in "$@"; do echo "$NP $POOLS $name" >> "$DFT/line.new"; done
  cat "$DFT/line.txt" >> "$DFT/line.new" 2>/dev/null || true
  mv "$DFT/line.new" "$DFT/line.txt"
else
  for name in "$@"; do echo "$NP $POOLS $name" >> "$DFT/line.txt"; done
fi

if running; then
  echo "the worker is running: $# cells added to its line"
  exit 0
fi
cat > "$DFT/worker.sh" <<WORKER
#!/usr/bin/env bash
# made by examples/gpu_jobs/dft_cells.sh: labels the cells of line.txt one after another
cd "$(pwd)"
# a container: no remote shell to look for, and no binding to processors the pod is not allotted
export OMPI_MCA_plm_rsh_agent="" PRTE_MCA_plm_ssh_agent="" OMPI_MCA_hwloc_base_binding_policy=none PRTE_MCA_hwloc_default_binding_policy=none
while true; do
  next=""
  while read -r np pools name pw; do
    if [ ! -f "$DFT/\$name.dft.xyz" ] && ! grep -qx "\$name" "$DFT/tried.txt" 2>/dev/null; then next="\$np \$pools \$name \$pw"; break; fi
  done < "$DFT/line.txt"
  [ -n "\$next" ] || break
  read -r np pools name pw <<< "\$next"
  pw="\${pw:-$DEFAULT_PW}"   # a fourth word in the line: another pw.x for that cell (dft_cells.sh usepw)
  echo "\$name" >> "$DFT/tried.txt"
  echo "\$(date -u +%FT%TZ) \$name: started on \$np processes in \$pools pools with \$pw" >> "$DFT/log.txt"
  "$PY" examples/mlip_frames.py label --frames "$DFT" --names "\$name" --pw "\$pw" --pseudo "$PSEUDO" \\
    --np "\$np" --pools "\$pools" --check-tol "${DFT_CHECK_TOL:-1e-6}" >> "$DFT/log.txt" 2>&1 \\
    || echo "\$(date -u +%FT%TZ) \$name: the labelling step failed" >> "$DFT/log.txt"
done
echo "\$(date -u +%FT%TZ) the line is empty" >> "$DFT/log.txt"
WORKER
if [ "${DFT_DRY:-0}" = "1" ]; then echo "DFT_DRY=1: the worker is written and not started ($DFT/worker.sh)"; exit 0; fi
setsid nohup bash "$DFT/worker.sh" > "$DFT/worker.out" 2>&1 < /dev/null &
echo $! > "$DFT/worker.pid"
echo "the worker is started (its log: jobs/shared/$SUBDIR/log.txt): $# cells in line"
