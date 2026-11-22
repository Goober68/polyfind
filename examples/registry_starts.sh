#!/usr/bin/env bash
# The starts of the registry search of docs/MLIP_SECOND_OPINION.md 3.24, in the order they are relaxed in:
#
#   bash examples/registry_starts.sh [OUT]      (default docs/predictions/mlip/registry_starts.extxyz)
#
# From each of the four fourteen-monomer cells' labelled all-trans cell (tt) and lowest helix cell (hx):
# the registry 3.23 found by accident in the film's composition (neighbours +1 and -1 monomers along the
# two cell edges, every chain turned) for the other three; the helix cells with their neighbours
# two monomers along (a helix repeats in two); the mirror registry; the all-trans registries without the
# turn; the turn alone; the helix cells with their neighbours one monomer along.
set -euo pipefail
cd "$(dirname "$0")/.."
OUT="${1:-docs/predictions/mlip/registry_starts.extxyz}"
FOUR="film matrix leader nitrile"
# the turn: +117 degrees put the film composition's dipoles where its low cell has them; the two cells with a nitrile are
# turned the other way (-117, the dipoles as near the short edge from its other side), where the nitrile groups do not meet
turn() { case "$1" in leader|nitrile) echo -117 ;; *) echo 117 ;; esac; }
V=()
for c in matrix leader nitrile; do V+=("$c:tt:1,-1,$(turn $c)"); done
for c in $FOUR; do V+=("$c:hx:2,-2,0" "$c:hx:-2,2,0"); done
for c in matrix leader nitrile film; do V+=("$c:tt:-1,1,$(turn $c)"); done
for c in $FOUR; do V+=("$c:tt:1,-1,0" "$c:tt:-1,1,0"); done
for c in matrix leader nitrile; do V+=("$c:tt:0,0,$(turn $c)"); done
for c in $FOUR; do V+=("$c:hx:1,-1,0" "$c:hx:-1,1,0"); done
python -I examples/registry_cells.py build --labels docs/predictions/mlip/cells14_frames_pbe_d3bj.extxyz --out "$OUT" --cell \
  film:tt=row_vdf9-trfe4-cfe1-14_TTTT_chain1_in_all_s1 film:hx=row_vdf9-trfe4-cfe1-14_TGTGp_chain1_in_all_s1 \
  matrix:tt=row_vdf10-trfe4-14_TTTT_chain1_in_all_s1 matrix:hx=row_vdf10-trfe4-14_TGTGp_chain1_in_all_s1 \
  leader:tt=tt_leader_from_rows14an_r1 leader:hx=helix_leader_TGTGp_r1 \
  nitrile:tt=tt_nitrile_from_rows14an_r1 nitrile:hx=row_vdf9-trfe4-an1-14_GpTGT_cfe_in_all_s1 \
  --variant "${V[@]}"
