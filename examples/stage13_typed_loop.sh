#!/bin/bash
# The reversed unit's two torsion triples: fit them on the candidate's frozen cells
# (examples/stage13_subfit.py), relax the cells again under the result
# (examples/chain_split.py polyfind --freeze), fit again.  docs/CRYSTAL_STAGE.md 5.29.
#
#   POLYFIND_TRAINSET=... examples/stage13_typed_loop.sh START.json FROZEN.json TAG OUTDIR [ROUNDS]
#
# START is a crystal-stage result (21 or 27 variables), FROZEN the cells frozen under it.
# Run from the repository's root.  Each round writes OUTDIR/TAG_k.json (the vector),
# OUTDIR/frozen_TAG_k.json and OUTDIR/chainsplit_TAG_k.json.
set -e
start=$1; frozen=$2; tag=$3; out=$4; rounds=${5:-3}
for k in $(seq 1 "$rounds"); do
  python -W ignore examples/stage13_subfit.py --x0 "$start" --vary V1_FF,V2_FF,V3_FF,V1_HH,V2_HH,V3_HH \
      --rows lone,candidate,conformers --only "hh_,candidate,conformers" --frozen "$frozen" --maxfev 10 \
      --out "$out/${tag}_$k.json" > "$out/${tag}_$k.log" 2>&1
  python -W ignore examples/chain_split.py polyfind --fit "$tag $k=$out/${tag}_$k.json" --freeze "$out/frozen_${tag}_$k.json" \
      --out "$out/chainsplit_${tag}_$k.json" > "$out/chainsplit_${tag}_$k.log" 2>&1
  start="$out/${tag}_$k.json"; frozen="$out/frozen_${tag}_$k.json"
  grep "^end:" "$out/${tag}_$k.log"; grep "minus" "$out/chainsplit_${tag}_$k.log"
done
