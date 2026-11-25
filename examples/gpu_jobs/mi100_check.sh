#!/usr/bin/env bash
# What the PI's machine is, and whether its card, the potential and DFT give on it the numbers they
# gave where the work was done.  A job of examples/mi100_runner.sh, the first of its queue:
#
#   mi100_check.sh
#
# Four parts, each standing if the one before it fails:
#   1. the machine: processors, memory, disk, system
#   2. the card: the kernel's device (/dev/kfd), the vendor's tools if there are any, and what
#      PyTorch's ROCm build sees; a product of two matrices on it against the processors'
#   3. the potential: the stored long-cell states of the urea row (relaxed on the other machine's
#      processors, energies in double precision) evaluated again on the processors and on the card,
#      with timings at 144 and 576 atoms
#   4. DFT: two stored single points of 36 atoms run again, against the energies stored for them
#
# MI100_CHECK_SKIP_DFT=1 leaves part 4 out (it takes minutes).
set -u
# shellcheck disable=SC1091
. "${MI100_ROOT:-$HOME/polyfind_mi100}/env.sh"

echo "== 1. the machine"
uname -srm
( . /etc/os-release 2>/dev/null && echo "${PRETTY_NAME:-}" )
lscpu 2>/dev/null | grep -E "^(Model name|Socket\(s\)|Core\(s\) per socket|Thread\(s\) per core|CPU\(s\)|NUMA node\(s\))" | sed 's/  */ /g'
free -g | sed -n 1,2p
df -h "$MI100_ROOT" | tail -n 1 | awk '{print "disk under the working folder: " $4 " free of " $2}'
echo "DFT will run on $TRIAGE_NP processes"

echo
echo "== 2. the card"
if [ -e /dev/kfd ]; then
  echo "/dev/kfd is there ($(stat -c '%A, group %G' /dev/kfd))"
  [ -r /dev/kfd ] && [ -w /dev/kfd ] || echo "THIS USER CANNOT OPEN /dev/kfd: add the user to that group (usually render and video) and log in again"
else
  echo "NO /dev/kfd: the kernel's driver for the card's computing side is not loaded (in WSL it never is)"
fi
[ -r /sys/module/amdgpu/version ] && echo "amdgpu driver $(cat /sys/module/amdgpu/version)"
(lspci 2>/dev/null | grep -i -E "vga|display|accelerat" | cut -d: -f3- | sed 's/^ *//') || true
for tool in rocm-smi amd-smi rocminfo; do
  if command -v "$tool" >/dev/null; then
    echo "-- $tool"
    case "$tool" in
      rocm-smi) rocm-smi --showproductname --showmeminfo vram 2>&1 | grep -v '^=*$' | head -n 20 ;;
      amd-smi) amd-smi version 2>&1 | head -n 5 ;;
      rocminfo) rocminfo 2>&1 | grep -E "Marketing Name|Name: *gfx|Compute Unit" | head -n 8 ;;
    esac
  fi
done
"$PY" - <<'EOF'
import time

import torch

print("PyTorch", torch.__version__, "| built for ROCm", getattr(torch.version, "hip", None), "| card seen:", torch.cuda.is_available())
if torch.cuda.is_available():
    props = torch.cuda.get_device_properties(0)
    print("card:", torch.cuda.get_device_name(0), f"{props.total_memory / 2**30:.1f} GB")
    try:
        a = torch.randn(3000, 3000, dtype=torch.float64, generator=torch.Generator().manual_seed(1))
        t0 = time.time()
        on_cpu = float((a @ a).abs().sum())
        t_cpu = time.time() - t0
        b = a.cuda()
        float((b @ b).abs().sum())  # the first call loads the card's routines
        t0 = time.time()
        on_card = float((b @ b).abs().sum())
        torch.cuda.synchronize()
        t_card = time.time() - t0
        print(f"a 3000 x 3000 product in double precision: processors {t_cpu:.3f} s, card {t_card:.3f} s, the two sums differ by {abs(on_card - on_cpu) / on_cpu:.1e} of the sum")
    except Exception as err:  # a card that is seen and cannot compute says how
        print("THE CARD IS SEEN AND A PRODUCT ON IT FAILED:", type(err).__name__, str(err)[:300])
else:
    print("THE CARD IS NOT SEEN BY THIS BUILD OF PYTORCH: the potential's jobs can still run on the processors (--device cpu)")
EOF

echo
echo "== 3. the potential"
"$PY" examples/chain_triage.py seed --frames "$TRIAGE_FRAMES" 2>&1 | tail -n 1 | cut -c1-300
"$PY" - <<'EOF' 2>&1 | grep -v -i "warn\|torch.load\|_Jd, _W3j\|cuequivariance"
import importlib.util
import os
import time
from pathlib import Path

import numpy as np
import torch
from ase.io import read

spec = importlib.util.spec_from_file_location("chain_triage", "examples/chain_triage.py")
ct = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ct)
frames = Path(os.environ["TRIAGE_FRAMES"])
devices = ["cpu"] + (["cuda"] if torch.cuda.is_available() else [])
stored = frames / "urea1_tt_x1_kink4.path.xyz"
calcs = {}
for device in devices:
    for dtype in ("float64", "float32"):
        try:
            calcs[device, dtype] = ct.model(ct.MODEL_DEFAULT, dtype=dtype, device=device)
        except Exception as err:
            print(f"THE POTENTIAL DID NOT LOAD on {device} in {dtype}:", type(err).__name__, str(err)[:300])
if stored.exists():
    for image in read(stored, index=":"):
        was = float(image.info["triage_kink_ev"])
        for (device, dtype), calc in calcs.items():
            if dtype != "float64":
                continue
            try:
                image.calc = calc
                e = float(image.get_potential_energy())
                f = float(np.linalg.norm(image.get_forces(), axis=1).max())
                print(f"stored state, {int(image.info['triage_kink_turned'])} repeats turned, {device}: {e - was:+.2e} eV from the stored energy of {len(image)} atoms, largest force {f:.3f} eV/A")
            except Exception as err:
                print(f"THE POTENTIAL FAILED on {device}:", type(err).__name__, str(err)[:300])
else:
    print("no stored long-cell states in the repository yet: timings only")
a = read(frames / "urea1_tt_x1.relaxed.xyz")
for rows, long in ((4, 4), (8, 8)):
    row = ct.kink_row(a, rows, long, long)
    reference = None
    for (device, dtype), calc in calcs.items():
        try:
            row.calc = calc
            t0 = time.time()
            e, f = float(row.get_potential_energy()), row.get_forces()
            first = time.time() - t0
            moved = row.copy()
            moved.rattle(0.01, seed=1)
            moved.calc = calc
            t0 = time.time()
            moved.get_potential_energy(), moved.get_forces()
            again = time.time() - t0
            if reference is None:
                reference = (e, f)
            print(f"{len(row)} atoms on {device} in {dtype}: {first:.2f} s the first time, {again:.3f} s after; against the first line {e - reference[0]:+.2e} eV, forces {np.abs(f - reference[1]).max():.1e} eV/A")
        except Exception as err:
            print(f"THE POTENTIAL FAILED on {device} in {dtype} at {len(row)} atoms:", type(err).__name__, str(err)[:300])
EOF

echo
echo "== 4. DFT"
if [ -n "${MI100_CHECK_SKIP_DFT:-}" ]; then
  echo "left out (MI100_CHECK_SKIP_DFT)"
else
  "$PY" examples/chain_triage.py points --frames "$TRIAGE_FRAMES" --names urea1_tt_x1_wall --ends 2>&1 | tail -n 1 | cut -c1-300
  grep -h -E "PWSCF .*WALL" "$TRIAGE_FRAMES"/work/urea1_tt_x1_wall_ends/image_*.out 2>/dev/null | sed 's/^ */time: /'
  "$PY" - <<'EOF'
import json
import os
from pathlib import Path

here = json.loads((Path(os.environ["TRIAGE_FRAMES"]) / "urea1_tt_x1_wall.path.json").read_text())["energy_ev"]
for q in json.loads(Path("docs/predictions/chain_triage.json").read_text())["paths"]:
    if q["name"] == "urea1_tt_x1_wall":
        was = [q["energy_ev"][0], q["energy_ev"][-1]]
        for label, e, w in zip(("first image", "last image"), (here[0], here[-1]), was):
            print(f"{label}: {e:.6f} eV here, {w:.6f} stored, {(e - w) / 13.605693123:+.1e} Ry apart (another processor moves the eighth or ninth decimal)")
EOF
fi
