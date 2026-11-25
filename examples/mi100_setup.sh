#!/usr/bin/env bash
# Sets up a Linux machine with an AMD card (the PI's Instinct MI100) for this repository's chain
# jobs.  Nothing is installed outside one folder and nothing needs root:
#
#   bash examples/mi100_setup.sh            (examples/mi100_runner.sh runs it first)
#
#   MI100_ROOT         the folder (default ~/polyfind_mi100)
#   MI100_TORCH        PyTorch's version (default 2.9.1)
#   MI100_TORCH_INDEX  where its build for the card is fetched from
#                      (default https://download.pytorch.org/whl/rocm6.4)
#
# What it puts under MI100_ROOT, each once:
#
#   mm/qe       Quantum ESPRESSO 7.5 from conda-forge, the build every number of
#               examples/chain_triage.py was computed with.  Its card support is NVIDIA's only: on
#               this machine DFT is the processors' work.
#   mm/py       a Python 3.12 with PyTorch's ROCm build (ROCm is AMD's CUDA; the wheel carries its
#               own ROCm libraries and needs only the kernel's driver, /dev/kfd), mace-torch 0.3.17,
#               ase and rdkit.  The card answers to the device name "cuda" in it.
#   egret-public  the pretrained potential (rowansci/egret-public at the commit used so far)
#   pseudo, pseudo_nc  the pseudopotential files, by checksum: five from the SSSP 1.3.0 PBE
#               precision library, three from examples/pseudo/
#   env.sh      what a job sources: PY, TRIAGE_PW, TRIAGE_PSEUDO, TRIAGE_PSEUDO_NC, TRIAGE_MODEL,
#               TRIAGE_NP (the physical cores), TRIAGE_FRAMES
#
# Whether the card is seen is NOT decided here: examples/gpu_jobs/mi100_check.sh reports it, with
# the processors, the memory and timings, and a card that is not seen leaves the DFT jobs standing.
set -euo pipefail

ROOT="${MI100_ROOT:-$HOME/polyfind_mi100}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MM="$ROOT/mm"
mkdir -p "$ROOT"
say() { echo "$(date -u +%FT%TZ) $*"; }

[ "$(uname -s)" = "Linux" ] && [ "$(uname -m)" = "x86_64" ] || { say "this set-up is for Linux on x86_64 (this is $(uname -s) $(uname -m))"; exit 1; }
for tool in curl tar git; do command -v "$tool" >/dev/null || { say "no $tool on this machine"; exit 1; }; done

if [ ! -x "$MM/bin/micromamba" ]; then
  mkdir -p "$MM"
  curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xj -C "$MM" bin/micromamba
fi
export MAMBA_ROOT_PREFIX="$MM/root"
if [ ! -x "$MM/qe/bin/pw.x" ]; then
  "$MM/bin/micromamba" create -y -q -p "$MM/qe" -c conda-forge "qe=7.5"
fi
# (the program is not run for its version: started without an input it leaves a CRASH file where it stands)
say "DFT: $(ls "$MM/qe/conda-meta" | grep -o '^qe-[0-9.]*' | head -n 1) at $MM/qe/bin/pw.x"

if [ ! -x "$MM/py/bin/python" ]; then
  "$MM/bin/micromamba" create -y -q -p "$MM/py" -c conda-forge "python=3.12" pip
fi
PY="$MM/py/bin/python"
if ! "$PY" -c "import torch" 2>/dev/null; then
  "$PY" -m pip install -q "torch==${MI100_TORCH:-2.9.1}" --index-url "${MI100_TORCH_INDEX:-https://download.pytorch.org/whl/rocm6.4}"
fi
if ! "$PY" -c "import mace, ase, rdkit, pytest" 2>/dev/null; then
  "$PY" -m pip install -q "mace-torch==0.3.17" ase rdkit pytest
fi
say "Python: $("$PY" -c "import torch, mace, ase; print('torch', torch.__version__, 'mace', mace.__version__, 'ase', ase.__version__)" 2>/dev/null | tail -n 1)"

EGRET="$ROOT/egret-public"
EGRET_COMMIT="ed85f3c7a792b5426dd73925e7f346590b60b2ad"
if [ ! -f "$EGRET/compiled_models/EGRET_1.model" ]; then
  rm -rf "$EGRET.part"
  git clone -q https://github.com/rowansci/egret-public "$EGRET.part"
  git -C "$EGRET.part" checkout -q "$EGRET_COMMIT"
  mv "$EGRET.part" "$EGRET"
fi
say "potential: $EGRET/compiled_models/EGRET_1.model ($(git -C "$EGRET" rev-parse --short HEAD))"

# the five files every DFT number of the project was computed with (examples/gpu_jobs/dft_cells.sh)
PSEUDO="$ROOT/pseudo"
SUMS="ba6ca6af365222805e3ce85564616b50a2feb562106dfbf1d725170cc7afd7b8  F.oncvpsp.upf
69d631780c0178c5239ed3dcb4a3e22c13bd13cde9beb97e123d036949da5c09  H_ONCV_PBE-1.0.oncvpsp.upf
36f1fc9a24ff5c788a536d4c1a53f82adca7b1fd132827fc654c124d69f274ca  N.oncvpsp.upf
9900d1efd50b9848e31849f39094b33348486b400ee51e0f3922f716137cf3d7  C.pbe-n-kjpaw_psl.1.0.0.UPF
a673d3adf336817175b7ffcbce2f811dd5929176c94a2c829cd1fc578ca924b8  Cl.pbe-n-rrkjus_psl.1.0.0.UPF"
if ! (cd "$PSEUDO" 2>/dev/null && echo "$SUMS" | sha256sum -c --quiet 2>/dev/null); then
  mkdir -p "$PSEUDO" "$ROOT/pseudo_fetch"
  curl -sL "https://archive.materialscloud.org/api/records/rcyfm-68h65/files/SSSP_1.3.0_PBE_precision.tar.gz/content" -o "$ROOT/pseudo_fetch/sssp.tar.gz"
  tar -xzf "$ROOT/pseudo_fetch/sssp.tar.gz" -C "$ROOT/pseudo_fetch"
  for f in F.oncvpsp.upf H_ONCV_PBE-1.0.oncvpsp.upf N.oncvpsp.upf C.pbe-n-kjpaw_psl.1.0.0.UPF Cl.pbe-n-rrkjus_psl.1.0.0.UPF; do
    cp "$(find "$ROOT/pseudo_fetch" -name "$f" | head -n 1)" "$PSEUDO/"
  done
  (cd "$PSEUDO" && echo "$SUMS" | sha256sum -c --quiet) || { say "the fetched pseudopotential files are not the ones the numbers were computed with"; exit 1; }
  rm -rf "$ROOT/pseudo_fetch"
fi
# the three that library does not have (examples/pseudo/README.md), and the norm-conserving set
OWN="2b6ba2a535fad72e0dd8046da3f14876a13829f33e398447fc366012bd7b8999  O.pbe-kjpaw.UPF
a30eca0fbeec806e3015e87dd58d945146f7e10b3337380053d9bf4089e71504  C.pbe-mt_gipaw.UPF
c8d54f0c9425650b32bb4f1e9e56b009839ab50f524176a5122d8a952ddac3ba  O.pbe-tm-triage.UPF"
(cd "$HERE/examples/pseudo" && echo "$OWN" | sha256sum -c --quiet) || { say "the files of examples/pseudo are not the ones stored"; exit 1; }
PSEUDO_NC="$ROOT/pseudo_nc"
mkdir -p "$PSEUDO_NC"
cp "$HERE/examples/pseudo/O.pbe-kjpaw.UPF" "$PSEUDO/"
cp "$HERE/examples/pseudo/C.pbe-mt_gipaw.UPF" "$HERE/examples/pseudo/O.pbe-tm-triage.UPF" "$PSEUDO_NC/"
cp "$PSEUDO/H_ONCV_PBE-1.0.oncvpsp.upf" "$PSEUDO/F.oncvpsp.upf" "$PSEUDO/N.oncvpsp.upf" "$PSEUDO_NC/"
say "pseudopotentials: eight files, checksums as stored"

# physical cores: DFT's processes gain nothing from a core's second thread
CORES="$(lscpu -p=CORE,SOCKET 2>/dev/null | grep -v '^#' | sort -u | wc -l)"
[ "${CORES:-0}" -ge 1 ] || CORES="$(nproc)"
cat > "$ROOT/env.sh" <<EOF
export MI100_ROOT="$ROOT"
export PY="$PY"
export TRIAGE_PW="$MM/qe/bin/pw.x"
export TRIAGE_PSEUDO="$PSEUDO"
export TRIAGE_PSEUDO_NC="$PSEUDO_NC"
export TRIAGE_MODEL="$EGRET/compiled_models/EGRET_1.model"
export TRIAGE_NP="\${MI100_NP:-$CORES}"
export TRIAGE_FRAMES="$ROOT/triage"
EOF
say "set up under $ROOT ($CORES physical cores); a job sources $ROOT/env.sh"
