#!/usr/bin/env bash
# The starts of the screen of docs/MLIP_SECOND_OPINION.md 3.25: fourteen-monomer compositions, each made of the film
# composition's potential-relaxed cells (examples/screen14.py), six cells a row.
#
#   bash examples/screen14_starts.sh [OUT_A [OUT_B [OUT_C [OUT_D]]]]
#       (default docs/predictions/mlip/screen14_starts_a.extxyz, screen14_starts_b.extxyz, ..._c.extxyz and ..._d.extxyz)
#
# File A: two controls, then the rows that put a chlorine unit and a reversed VDF unit together and the steps between them
#   and the leader.  The film's composition made of itself (nothing changes: every cell must come back as it is), and the
#   matrix (its margins in the potential and in DFT are known from its own cells: 3.22, 3.24).  The leader is
#   IVRVVRVAVVVVVR: a reversed unit first, the nitrile eighth, a VDF eleventh where the film's composition has its
#   fourth TrFE.  Its nitrile is what the potential cannot do (3.24), so the candidates take the leader's other two
#   features onto a chlorine unit: CFE (C) and CTFE (K), with the reversed unit, with the eleventh unit a VDF, and with both.
# File B: the same features with no pin at all, more TrFE, the other units the rows have in place of the CFE, and the
#   leader itself as a control (the potential's own numbers for it are known; it is not a candidate here).
# File C: more all-trans cells for every row.  A made cell keeps the pattern of backbone twists of the cell it is made of,
#   and the patterns are minima of their own (the matrix made of the film composition's upright cell rests 1.01 kcal/mol
#   per chain above the matrix's own upright cell).  Files A and B give each row the film composition's two patterns (its
#   low cell's and its upright cell's).  Here: tt_reg, made of the matrix's own upright cell (the regular pattern of a chain
#   with no pin), and tt_own, the row's own upright cell where there is one: polyfind's cell of it under chain1-cfe
#   (examples/row_cells_frames.py write, the cell free to shear, as built), or the labelled one (the leader's).
# File D: more helix cells for every row.  The same holds on the helix side: the leader's own lowest helix cell is a
#   packing of its own (5.25 by 9.84 A with the angle between the chain and the long edge 94 degrees; the film
#   composition's four are 5.14 to 5.25 by 9.9 to 10.1 with the cell leaning 3 to 8 degrees two ways), and the leader made
#   of the film composition's four helix cells rests 4.4 kcal/mol per chain above it at the least.  So: hx_ldr, every row
#   made of the leader's lowest helix cell (the nitrile made a chlorine or what the row has there, its nitrogen removed),
#   and for six rows hx_an, the same of the nitrile composition's lowest helix cell (a third packing: 5.58 by 9.75, the
#   other angle 74).  The leader's hx_ldr is its own cell: nothing changes.  And the all-trans side once more: the leader's
#   own upright cell is 0.61 below the leader made of the film composition's, so tt_ldr, every other row made of the
#   leader's own upright all-trans cell.  The rows come in the order they matter (the two in the window first): the
#   pod's job stops at its budget.
set -euo pipefail
cd "$(dirname "$0")/.."
OUT_A="${1:-docs/predictions/mlip/screen14_starts_a.extxyz}"
OUT_B="${2:-docs/predictions/mlip/screen14_starts_b.extxyz}"
OUT_C="${3:-docs/predictions/mlip/screen14_starts_c.extxyz}"
OUT_D="${4:-docs/predictions/mlip/screen14_starts_d.extxyz}"
M=docs/predictions/mlip
FILM=vdf9-trfe4-cfe1-14
FRAMES="$M/cells14_tight_in_all_s1.extxyz $M/paths/all_trans_tight.extxyz $M/cells14_frames_pbe_d3bj.extxyz"
SOURCES=(tt_low=film_r_released_tight tt_up=row_${FILM}_TTTT_chain1_in_all_s1 hx_TGTGp=row_${FILM}_TGTGp_chain1_in_all_s1
  hx_GpTGT=row_${FILM}_GpTGT_chain1_in_all_s1 hx_GTGpT=row_${FILM}_GTGpT_chain1_in_all_s1 hx_TGpTG=row_${FILM}_TGpTG_chain1_in_all_s1)
python -I examples/screen14.py starts --frames $FRAMES --source "${SOURCES[@]}" --out "$OUT_A" --row \
  vdf9-trfe4-cfe1-14=VVRVVRVCVVRVVR \
  vdf8-trfe4-cfe1-rvdf1-14=IVRVVRVCVVRVVR \
  vdf9-trfe3-cfe1-rvdf1-14=IVRVVRVCVVVVVR \
  vdf10-trfe3-cfe1-14=VVRVVRVCVVVVVR \
  vdf10-trfe4-14=VVRVVRVVVVRVVR \
  vdf8-trfe4-ctfe1-rvdf1-14=IVRVVRVKVVRVVR \
  vdf9-trfe3-ctfe1-rvdf1-14=IVRVVRVKVVVVVR \
  vdf9-trfe4-ctfe1-14=VVRVVRVKVVRVVR
python -I examples/screen14.py starts --frames $FRAMES --source "${SOURCES[@]}" --out "$OUT_B" --row \
  vdf9-trfe4-rvdf1-first-14=IVRVVRVVVVRVVR \
  vdf10-trfe3-rvdf1-14=IVRVVRVVVVVVVR \
  vdf11-trfe3-14=VVRVVRVVVVVVVR \
  vdf9-trfe5-14=VVRVVRVRVVRVVR \
  vdf9-trfe4-cdfe1-14=VVRVVRVDVVRVVR \
  vdf8-trfe4-cfe1-cdfe1-14=DVRVVRVCVVRVVR \
  vdf9-trfe4-rvdf1-14=VVRVVRVIVVRVVR \
  vdf9-trfe4-e1-14=VVRVVRVEVVRVVR \
  vdf9-trfe3-an1-rvdf1-14=IVRVVRVAVVVVVR
# file C
OWN="$(mktemp -d)"
trap 'rm -rf "$OWN"' EXIT
HAVE="vdf8-trfe4-cfe1-rvdf1-14 vdf8-trfe4-ctfe1-rvdf1-14 vdf9-trfe4-ctfe1-14 vdf9-trfe5-14 vdf9-trfe4-cdfe1-14 vdf8-trfe4-cfe1-cdfe1-14 vdf9-trfe4-rvdf1-14 vdf9-trfe4-e1-14"
AS_IS=(vdf9-trfe3-an1-rvdf1-14:tt_own=tt_leader_from_rows14an_r1)
for row in $HAVE; do
  python examples/row_cells_frames.py write --margins docs/predictions/margins_chain1-cfe.json --row "$row" --confs TTTT --kmesh "4 2 1" --out "$OWN" > /dev/null
  AS_IS+=("$row:tt_own=row_${row}_TTTT_cfe")
done
python -I examples/screen14.py starts --frames $FRAMES "$OWN"/*.start.xyz --source tt_reg=row_vdf10-trfe4-14_TTTT_chain1_in_all_s1 --as-is "${AS_IS[@]}" --out "$OUT_C" --row \
  vdf8-trfe4-cfe1-rvdf1-14=IVRVVRVCVVRVVR \
  vdf9-trfe3-cfe1-rvdf1-14=IVRVVRVCVVVVVR \
  vdf10-trfe3-cfe1-14=VVRVVRVCVVVVVR \
  vdf8-trfe4-ctfe1-rvdf1-14=IVRVVRVKVVRVVR \
  vdf9-trfe3-ctfe1-rvdf1-14=IVRVVRVKVVVVVR \
  vdf9-trfe4-ctfe1-14=VVRVVRVKVVRVVR \
  vdf9-trfe4-cfe1-14=VVRVVRVCVVRVVR \
  vdf9-trfe4-rvdf1-first-14=IVRVVRVVVVRVVR \
  vdf10-trfe3-rvdf1-14=IVRVVRVVVVVVVR \
  vdf11-trfe3-14=VVRVVRVVVVVVVR \
  vdf9-trfe5-14=VVRVVRVRVVRVVR \
  vdf9-trfe4-cdfe1-14=VVRVVRVDVVRVVR \
  vdf8-trfe4-cfe1-cdfe1-14=DVRVVRVCVVRVVR \
  vdf9-trfe4-rvdf1-14=VVRVVRVIVVRVVR \
  vdf9-trfe4-e1-14=VVRVVRVEVVRVVR \
  vdf9-trfe3-an1-rvdf1-14=IVRVVRVAVVVVVR
# file D
LDR="hx_ldr=helix_leader_TGTGp_r1_in_all_s1 tt_ldr=tt_leader_from_rows14an_r1"
AN="hx_an=row_vdf9-trfe4-an1-14_GpTGT_cfe_in_all_s1"
python -I examples/screen14.py starts --frames $FRAMES --source $LDR $AN --out "$OWN/d1.extxyz" --row \
  vdf9-trfe3-cfe1-rvdf1-14=IVRVVRVCVVVVVR \
  vdf8-trfe4-cfe1-rvdf1-14=IVRVVRVCVVRVVR
python -I examples/screen14.py starts --frames $FRAMES --source hx_ldr=helix_leader_TGTGp_r1_in_all_s1 $AN --out "$OWN/d2.extxyz" --row \
  vdf9-trfe3-an1-rvdf1-14=IVRVVRVAVVVVVR
python -I examples/screen14.py starts --frames $FRAMES --source $LDR $AN --out "$OWN/d3.extxyz" --row \
  vdf9-trfe4-cfe1-14=VVRVVRVCVVRVVR \
  vdf10-trfe4-14=VVRVVRVVVVRVVR \
  vdf10-trfe3-rvdf1-14=IVRVVRVVVVVVVR
python -I examples/screen14.py starts --frames $FRAMES --source $LDR --out "$OWN/d4.extxyz" --row \
  vdf9-trfe4-rvdf1-first-14=IVRVVRVVVVRVVR \
  vdf9-trfe3-ctfe1-rvdf1-14=IVRVVRVKVVVVVR \
  vdf8-trfe4-ctfe1-rvdf1-14=IVRVVRVKVVRVVR \
  vdf9-trfe5-14=VVRVVRVRVVRVVR \
  vdf9-trfe4-cdfe1-14=VVRVVRVDVVRVVR \
  vdf10-trfe3-cfe1-14=VVRVVRVCVVVVVR \
  vdf9-trfe4-ctfe1-14=VVRVVRVKVVRVVR \
  vdf11-trfe3-14=VVRVVRVVVVVVVR \
  vdf8-trfe4-cfe1-cdfe1-14=DVRVVRVCVVRVVR \
  vdf9-trfe4-e1-14=VVRVVRVEVVRVVR \
  vdf9-trfe4-rvdf1-14=VVRVVRVIVVRVVR
cat "$OWN"/d[1234].extxyz > "$OUT_D"
echo "file D: $(grep -c '^1[67][0-9]$' "$OUT_D") cells in $OUT_D"
