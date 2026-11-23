#!/usr/bin/env bash
# The starts of a cell of two different chains (docs/MLIP_SECOND_OPINION.md 3.26): one chain of IVRVVRVCVVVVVR and one
# of IVRVVRVCVVRVVR, the two rows of the screen that differ by one TrFE in fourteen and sit at +0.20 and -0.03 in the
# potential (3.25).  Seven TrFE in twenty-eight: 25 mol%, in a cell of the screen's size.
#
#   bash examples/screen14_mixed_starts.sh [OUT_E [OUT_SHEETS]]
#       (default docs/predictions/mlip/screen14_starts_e.extxyz and screen14_sheets_3.extxyz)
#
# Made as the screen's rows were (examples/screen14.py: substituents changed where they stand in a relaxed cell), in
# every packing the screen learned a row must be given (3.25): all-trans in the film composition's low pattern and its
# upright one, the matrix's and the leader's; helix in the film composition's three lower placements and in the
# leader's packing.  A cell of two chains has two ways of putting two different chains in it, and both are made: a
# kind ending in x has the chains the other way round.  Sixteen cells, the ones that matter first.
# The kink sheets: on the chlorine, where both parents' cheapest kink is (sites 12 and 14m of the film composition's
# sheets; 12 of the leader's), and on the run of VDF units (16m of the leader's), both ways round: eight sheets, their
# label f or l for whose sheet it was made of, x for the chains the other way round.
set -euo pipefail
cd "$(dirname "$0")/.."
OUT_E="${1:-docs/predictions/mlip/screen14_starts_e.extxyz}"
OUT_S="${2:-docs/predictions/mlip/screen14_sheets_3.extxyz}"
M=docs/predictions/mlip
FILM=vdf9-trfe4-cfe1-14
ROW=vdf17-trfe7-cfe2-rvdf2-28
A=IVRVVRVCVVVVVR; B=IVRVVRVCVVRVVR
FRAMES="$M/cells14_tight_in_all_s1.extxyz $M/paths/all_trans_tight.extxyz $M/cells14_frames_pbe_d3bj.extxyz"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
sources() {  # the kinds in the order they matter, each with the suffix given
  local x="$1"
  echo "hx_ldr$x=helix_leader_TGTGp_r1_in_all_s1 hx_GTGpT$x=row_${FILM}_GTGpT_chain1_in_all_s1 tt_low$x=film_r_released_tight tt_up$x=row_${FILM}_TTTT_chain1_in_all_s1" \
       "tt_ldr$x=tt_leader_from_rows14an_r1 hx_TGTGp$x=row_${FILM}_TGTGp_chain1_in_all_s1 hx_GpTGT$x=row_${FILM}_GpTGT_chain1_in_all_s1 tt_reg$x=row_vdf10-trfe4-14_TTTT_chain1_in_all_s1"
}
# shellcheck disable=SC2046
python -I examples/screen14.py starts --frames $FRAMES --source $(sources "") --out "$TMP/ab.extxyz" --row "$ROW=$A+$B"
# shellcheck disable=SC2046
python -I examples/screen14.py starts --frames $FRAMES --source $(sources x) --out "$TMP/ba.extxyz" --row "$ROW=$B+$A"
# the two ways round in turn, so that a job stopped at its budget has both of the cells that matter most
python -I - "$TMP/ab.extxyz" "$TMP/ba.extxyz" "$OUT_E" <<'PY'
import sys
from ase.io import read, write
a, b = read(sys.argv[1], index=":"), read(sys.argv[2], index=":")
out = [x for pair in zip(a, b) for x in pair]
write(sys.argv[3], out, format="extxyz")
print(f"{len(out)} starts written to {sys.argv[3]}")
PY
S="$M/cells14_tight_in_all_s1.extxyz"
python -I examples/screen14.py sheets --frames $S --of "$FILM" --label f --sites 12 14m --row "$ROW=$A+$B" --out "$TMP/f.extxyz"
python -I examples/screen14.py sheets --frames $S --of "$FILM" --label fx --sites 12 14m --row "$ROW=$B+$A" --out "$TMP/fx.extxyz"
python -I examples/screen14.py sheets --frames $S --of vdf9-trfe3-an1-rvdf1-14 --label l --sites 12 16m --row "$ROW=$A+$B" --out "$TMP/l.extxyz"
python -I examples/screen14.py sheets --frames $S --of vdf9-trfe3-an1-rvdf1-14 --label lx --sites 12 16m --row "$ROW=$B+$A" --out "$TMP/lx.extxyz"
python -I - "$TMP/f.extxyz" "$TMP/fx.extxyz" "$TMP/l.extxyz" "$TMP/lx.extxyz" "$OUT_S" <<'PY'
import sys
from ase.io import read, write
files = [read(f, index=":") for f in sys.argv[1:5]]
out = [x for group in zip(*files) for x in group]
write(sys.argv[5], out, format="extxyz")
print(f"{len(out)} sheets written to {sys.argv[5]}")
PY
