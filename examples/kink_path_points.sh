#!/usr/bin/env bash
# The states and tops of the kink paths of two fourteen-monomer cells (docs/MLIP_SECOND_OPINION.md 3.23) against each
# cell's sheet as the potential relaxed it, in the potential and in DFT, the steps between them, and the four axes with
# the frequency read from each step's barrier:
#
#   bash examples/kink_path_points.sh [OUTDIR]        (default docs/predictions)
#
# writes kink_paths_film.json, kink_paths_leader.json and four_axes_by_barrier_{forming,hop,removal}.json.  Every point
# but the two marked "stepped" is a cell relaxed in the potential, and DFT is a single point on it: one footing.
set -euo pipefail
cd "$(dirname "$0")/.."
OUT="${1:-docs/predictions}"
P=docs/predictions/mlip/paths
L="docs/predictions/mlip/cells14_frames_pbe_d3bj.extxyz docs/predictions/mlip/cells14_third_pod_frames_pbe_d3bj.extxyz"
CFE=vdf9-trfe4-cfe1-14; LEADER=vdf9-trfe3-an1-rvdf1-14
STEPS=("the kink forms, both torsions together=all-trans, the sheet's registry>top: both torsions together>sheet"
  "the second gauche bond moves a monomer along=sheet>top: sheet to stretched kink>stretched kink"
  "the stretched kink forms=all-trans, the sheet's registry>top: stretched kink to all-trans>stretched kink"
  "upright, the kink forms=all-trans, upright>top: upright all-trans to upright kinked>kinked, upright")
python -I examples/mlip_path.py points --name "film's composition" --scans $P/film_*.extxyz $P/all_trans_tight.extxyz --labels $L \
  --zero "sheet_${CFE}_belt4_at14m_sheared_in_all_s1" --out "$OUT/kink_paths_film.json" --point \
  "sheet=sheet_${CFE}_belt4_at14m_sheared_in_all_s1" \
  "sheet, stepped toward DFT's minimum=film_r_i00|sheet_${CFE}_belt4_at14m_sheared_in_all_s1_in_all_s1" \
  "all-trans, upright=row_${CFE}_TTTT_chain1_in_all_s1" \
  "all-trans, upright, stepped toward DFT's minimum=film_f_i00|row_${CFE}_TTTT_chain1_in_all_s1_in_all_s1" \
  "all-trans, the sheet's registry=film_r_released" \
  "all-trans, the sheet's registry, relaxed on=film_r_released_tight" \
  "all-trans, upright, every chain turned=film_all_trans_turned_117_tight" \
  "kinked, upright=film_f_released" \
  "stretched kink=film_b1_released" \
  "top: upright all-trans to upright kinked=film_f_i06" \
  "top: both torsions together=film_r_i07" \
  "the image before that top=film_r_i06" \
  "the image after that top=film_r_i08" \
  "that top in finer steps=film_r_fine_i04" \
  "top: sheet to stretched kink=film_b1_i07" \
  "top: stretched kink to all-trans=film_b2_i07" \
  --between "${STEPS[@]}"
python -I examples/mlip_path.py points --name "leader" --scans $P/leader_*.extxyz $P/all_trans_tight.extxyz --labels $L \
  --zero "leader_r_i00|sheet_${LEADER}_belt4_at18_sheared_in_all_s1" --out "$OUT/kink_paths_leader.json" --point \
  "sheet=leader_r_i00|sheet_${LEADER}_belt4_at18_sheared_in_all_s1" \
  "all-trans, upright=tt_leader_from_rows14an_r1" \
  "all-trans, upright, stepped toward DFT's minimum=leader_f_i00|tt_leader_from_rows14an_r1_in_all_s1" \
  "all-trans, the sheet's registry=leader_r_released" \
  "all-trans, the sheet's registry, relaxed on=leader_r_released_tight" \
  "kinked, upright=leader_f_released" \
  "stretched kink=leader_b1_released" \
  "top: upright all-trans to upright kinked=leader_f_i07" \
  "top: both torsions together=leader_r_i04" \
  "the image before that top=leader_r_i03" \
  "the image after that top=leader_r_i05" \
  "that top in finer steps=leader_r_fine_i05" \
  "top: sheet to stretched kink=leader_b1_i08" \
  "top: stretched kink to all-trans=leader_b2_i07" \
  "the image after the top from the stretched kink=leader_b2_i08" \
  --between "${STEPS[@]}"
# the four axes by each event's barrier in DFT: the film composition's row is the anchor of each
python - "$OUT" <<'PY'
import json, subprocess, sys
out = sys.argv[1]
film, leader = (json.load(open(f"{out}/kink_paths_{c}.json")) for c in ("film", "leader"))
def dft(t, label):
    return t["points"][label]["dft"]
events = {"forming": (film["steps"]["the kink forms, both torsions together"]["dft"]["forward"],
                      dft(leader, "top: stretched kink to all-trans") - dft(leader, "all-trans, the sheet's registry")),
          "hop": (dft(film, "top: sheet to stretched kink"), dft(leader, "top: sheet to stretched kink")),
          "removal": (dft(film, "top: both torsions together"), dft(leader, "top: stretched kink to all-trans"))}
for name, (f, l) in events.items():
    print(f"== the kink's {name}: the film's composition {f:.2f}, the leader {l:.2f} kcal/mol per chain (DFT)")
    subprocess.run([sys.executable, "examples/four_axes_at_margins.py", "--margins-from", "docs/predictions/cells14_states_dft.json", "--shifts=0", "--barriers",
                    f"vdf9-trfe4-cfe1-14={f:.4f}", f"vdf9-trfe3-an1-rvdf1-14={l:.4f}", "--out", f"{out}/four_axes_by_barrier_{name}.json"], check=True)
PY
