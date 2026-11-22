#!/usr/bin/env bash
# The four fourteen-monomer cells' margins from the lowest labelled all-trans cell and the lowest labelled helix cell of
# any registry (docs/MLIP_SECOND_OPINION.md 3.24), and the four axes at them:
#
#   bash examples/registry_margins.sh [OUTDIR]        (default docs/predictions)
#
# writes cells14_margins_registries.json, four_axes_registries.json and four_axes_registries_by_barrier_{hop,removal}.json.
set -euo pipefail
cd "$(dirname "$0")/.."
OUT="${1:-docs/predictions}"
M=docs/predictions/mlip
LABELS="$M/cells14_frames_pbe_d3bj.extxyz $M/cells14_third_pod_frames_pbe_d3bj.extxyz"
[ -f "$M/cells14_fourth_pod_frames_pbe_d3bj.extxyz" ] && LABELS="$LABELS $M/cells14_fourth_pod_frames_pbe_d3bj.extxyz"
python -I examples/registry_cells.py margins --labels $LABELS --out "$OUT/cells14_margins_registries.json" --cell \
  film:tt=row_vdf9-trfe4-cfe1-14_TTTT_chain1_in_all_s1 film:hx=row_vdf9-trfe4-cfe1-14_TGTGp_chain1_in_all_s1 \
  matrix:tt=row_vdf10-trfe4-14_TTTT_chain1_in_all_s1 matrix:hx=row_vdf10-trfe4-14_TGTGp_chain1_in_all_s1 \
  leader:tt=tt_leader_from_rows14an_r1 leader:hx=helix_leader_TGTGp_r1 \
  nitrile:tt=tt_nitrile_from_rows14an_r1 nitrile:hx=row_vdf9-trfe4-an1-14_GpTGT_cfe_in_all_s1 \
  --also film:tt=film_r_released_tight film:tt=film_r_released leader:tt=leader_r_released \
  --row film=vdf9-trfe4-cfe1-14 matrix=vdf10-trfe4-14 leader=vdf9-trfe3-an1-rvdf1-14 nitrile=vdf9-trfe4-an1-14
# the four axes at those margins: the kink sites as in 3.22 (the potential's map with DFT's energy at the eight sites that have one), and
# the film's composition and the leader with the frequency read from the barriers of 3.23 (the kink's end moving, and its removal).
# A kink site's energy is against its cell's all-trans cell, the upright one: where a lower all-trans cell is labelled every site of that
# row stands that much higher against it (--sheet-shift, DFT's difference between the two all-trans cells, per chain).
SHIFTS="$(python -I - "$OUT/cells14_margins_registries.json" <<'PY'
import json, sys
for c in json.load(open(sys.argv[1]))["compositions"].values():
    lower = -c["cells"]["tt"][c["lowest_tt"]]["dft"]
    if lower > 0:
        print(f"{c['row']}={lower:.4f}")
PY
)"
python examples/four_axes_at_margins.py --margins-from "$OUT/cells14_margins_registries.json" --shifts=0 --sheet-shift $SHIFTS \
  --sites docs/predictions/kink_site_scan_potential_card_dft_at_eight.json --out "$OUT/four_axes_registries.json"
python - "$OUT" <<'PY'
import json, subprocess, sys
out = sys.argv[1]
film, leader = (json.load(open(f"docs/predictions/kink_paths_{c}.json"))["points"] for c in ("film", "leader"))
top = max(leader["top: stretched kink to all-trans"]["dft"], leader["the image after the top from the stretched kink"].get("dft", -1e9))   # DFT's own top (3.24)
events = {"hop": (film["top: sheet to stretched kink"]["dft"], leader["top: sheet to stretched kink"]["dft"]),
          "removal": (film["top: both torsions together"]["dft"], top)}
for name, (f, l) in events.items():
    subprocess.run([sys.executable, "examples/four_axes_at_margins.py", "--margins-from", f"{out}/cells14_margins_registries.json", "--shifts=0", "--barriers",
                    f"vdf9-trfe4-cfe1-14={f:.4f}", f"vdf9-trfe3-an1-rvdf1-14={l:.4f}", "--out", f"{out}/four_axes_registries_by_barrier_{name}.json"], check=True)
PY
