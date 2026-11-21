#!/usr/bin/env bash
# Does the card's pw.x (examples/gpu_jobs/qe_gpu_build.sh) give the processors' DFT?
# A job of examples/runpod_jobs_pod.sh:
#
#   qe_gpu_check.sh FILE NAME LABELS [NP [POOLS]]
#
# First the anchor's own beta-PVDF input must return the anchor's energy (to 1e-5 Ry here; the
# difference is printed).  Then the frame NAME of FILE, which the processors' pw.x has already
# labelled (the same frame in LABELS), is labelled by the card's on NP processes (default 1) in
# POOLS pools (default 1), and the two are compared: energy, forces, stress and seconds.
#
# Into $JOB_OUT: gpu_check.json, the label and its pw.x output (jobs/shared/dft_gpu/ while it
# runs).  The job fails when the energies differ by more than 1e-5 Ry per cell or a force by more
# than 0.005 eV/A: then the card's pw.x labels nothing.
set -euo pipefail
FILE="$1"; NAME="$2"; LABELS="$3"; NP="${4:-1}"; POOLS="${5:-1}"
PW="${QE_GPU_PREFIX:-/opt/qe_gpu}/bin/pw.x"
PSEUDO="${DFT_PSEUDO:-/opt/pseudo}"
WORK="$(dirname "$JOB_OUT")/shared/dft_gpu"
[ -x "$PW" ] || { echo "no pw.x at $PW (qe_gpu_build.sh has not built it)"; exit 1; }
source "$(dirname "$0")/mlip_env.sh"
mkdir -p "$WORK"
export OMPI_MCA_plm_rsh_agent="" PRTE_MCA_plm_ssh_agent="" OMPI_MCA_hwloc_base_binding_policy=none PRTE_MCA_hwloc_default_binding_policy=none
"$PY" - "$FILE" "$WORK" "$NAME" <<'PY'
import sys
from pathlib import Path
from ase.io import read, write

frames = {str(a.info["frame"]): a for a in read(sys.argv[1], index=":")}
write(Path(sys.argv[2]) / f"{sys.argv[3]}.start.xyz", frames[sys.argv[3]], format="extxyz")
PY
"$PY" examples/mlip_frames.py label --frames "$WORK" --names "$NAME" --pw "$PW" --pseudo "$PSEUDO" --np "$NP" --pools "$POOLS" --check-tol 1e-5
grep -m 1 -E "GPU acceleration is ACTIVE|GPU-aware MPI|GPU acceleration" "$WORK/$NAME.pw.out" || echo "the output does not say the card was used"
cp "$WORK/$NAME.dft.xyz" "$WORK/$NAME.pw.out" "$WORK/check/check.json" "$JOB_OUT/" 2>/dev/null || true
"$PY" - "$WORK/$NAME.dft.xyz" "$LABELS" "$NAME" "$JOB_OUT/gpu_check.json" "$WORK/check/check.json" <<'PY'
import json
import sys

import numpy as np
from ase import units
from ase.io import read

card = read(sys.argv[1])
host = {str(a.info["frame"]): a for a in read(sys.argv[2], index=":")}[sys.argv[3]]
n = int(card.info["n_monomers"])
d_ry = float(card.info["dft_energy_ry"]) - float(host.info["dft_energy_ry"])
d_f = float(np.abs(card.arrays["REF_forces"] - host.arrays["REF_forces"]).max())
d_s = float(np.abs(np.asarray(card.info["REF_stress"]) - np.asarray(host.info["REF_stress"])).max() / units.GPa)
out = {"what": "one cell labelled by the pw.x built for the card against the same cell labelled by the processors' pw.x",
       "frame": sys.argv[3], "atoms": len(card), "energy_difference_ry_per_cell": d_ry, "energy_difference_meV_per_monomer": 1e3 * d_ry * units.Ry / n,
       "largest_force_difference_eV_per_A": d_f, "largest_stress_difference_GPa": d_s,
       "seconds_on_the_card": float(card.info["dft_seconds"]), "seconds_on_the_processors": float(host.info["dft_seconds"]),
       "anchor": json.load(open(sys.argv[5]))}
out["agrees"] = bool(abs(d_ry) < 1e-5 and d_f < 0.005)
json.dump(out, open(sys.argv[4], "w"), indent=1)
print(f"{sys.argv[3]}: the card's energy {d_ry:+.2e} Ry from the processors' ({out['energy_difference_meV_per_monomer']:+.4f} meV per monomer), forces within {d_f:.4f} eV/A, "
      f"stress within {d_s:.4f} GPa; {out['seconds_on_the_card']:.0f} s on the card, {out['seconds_on_the_processors']:.0f} s on the processors: "
      f"{'AGREES' if out['agrees'] else 'DOES NOT AGREE'}")
sys.exit(0 if out["agrees"] else 1)
PY
