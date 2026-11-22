#!/usr/bin/env bash
# The film composition's low all-trans cell of 3.23 given to the other three compositions: each one's upright all-trans
# cell with its chains turned as the film composition's are, its lattice strained across the chains as the film
# composition's is, and its neighbours the same monomers along (docs/MLIP_SECOND_OPINION.md 3.24).  The film
# composition's own is the control: it must come back to the cell 3.23 labelled.
#
#   bash examples/registry_starts_like.sh [OUT]      (default docs/predictions/mlip/registry_starts_like.extxyz)
set -euo pipefail
cd "$(dirname "$0")/.."
OUT="${1:-docs/predictions/mlip/registry_starts_like.extxyz}"
python -I examples/registry_cells.py build --labels docs/predictions/mlip/cells14_frames_pbe_d3bj.extxyz docs/predictions/mlip/paths/all_trans_tight.extxyz \
  --out "$OUT" --like row_vdf9-trfe4-cfe1-14_TTTT_chain1_in_all_s1 film_r_released_tight --cell \
  film:tt=row_vdf9-trfe4-cfe1-14_TTTT_chain1_in_all_s1 matrix:tt=row_vdf10-trfe4-14_TTTT_chain1_in_all_s1 \
  leader:tt=tt_leader_from_rows14an_r1 nitrile:tt=tt_nitrile_from_rows14an_r1 \
  --variant matrix:tt:like film:tt:like nitrile:tt:like leader:tt:like
