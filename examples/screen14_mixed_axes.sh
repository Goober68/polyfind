#!/usr/bin/env bash
# The four axes of the cell of two different chains (docs/MLIP_SECOND_OPINION.md 3.26): strain, tan delta, frequency
# and blocked force per area, its margin and its kink's hop put through the tables of docs/CANDIDATE_MATERIALS.md
# 5.6 (examples/four_axes_at_margins.py --barriers), as examples/screen14_axes.sh does for the screen's two rows.
#
#   bash examples/screen14_mixed_axes.sh [OUTDIR]        (default docs/predictions)
#
# The cell has no dipole, switch strain or stiffness of its own (belt4 has not built it), so it is put on the
# inputs of its parent IVRVVRVCVVRVVR, which belt4 has (no crystal stiffness), and on the leader's.  Its margin is
# DFT's from its lowest labelled cells (screen14_mixed_dft.json); its hop is the potential's
# (screen14_mixed_hops.json), against the film composition's 7.81 in the same potential.  The margin is also moved
# by 0.05 either way.  The film's composition is in the table at its DFT margin from its lowest labelled cells.
set -euo pipefail
cd "$(dirname "$0")/.."
OUT="${1:-docs/predictions}"
ROW=vdf17-trfe7-cfe2-rvdf2-28
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
python - "$OUT" "$TMP" "$ROW" <<'PY'
import json, sys
out, tmp, row = sys.argv[1:4]
margin = json.load(open(f"{out}/screen14_mixed_dft.json"))["rows"][row]["margin"]
hop = json.load(open(f"{out}/screen14_mixed_hops.json"))["paths"]["mix"]["hop"]
json.dump({"rows": {"vdf8-trfe4-cfe1-rvdf1-14": {"margin": margin}, "vdf9-trfe3-an1-rvdf1-14": {"margin": margin}, "vdf9-trfe4-cfe1-14": {"margin": 0.063}}}, open(f"{tmp}/m.json", "w"))
open(f"{tmp}/hop.txt", "w").write(f"{hop:.4f}\n")
PY
HOP="$(cat "$TMP/hop.txt")"
python examples/four_axes_at_margins.py --margins-from "$TMP/m.json" --shifts=-0.05,0,0.05 \
  --barriers vdf9-trfe4-cfe1-14=7.81 vdf8-trfe4-cfe1-rvdf1-14=$HOP vdf9-trfe3-an1-rvdf1-14=$HOP --out "$OUT/four_axes_screen14_mixed.json"
python - "$OUT" <<'PY'
import json, sys
p = f"{sys.argv[1]}/four_axes_screen14_mixed.json"
d = json.load(open(p))
d["stands_in"] = {"vdf8-trfe4-cfe1-rvdf1-14": "the mixed cell on IVRVVRVCVVRVVR's inputs", "vdf9-trfe3-an1-rvdf1-14": "the mixed cell on the leader's inputs"}
d["note"] = ("examples/screen14_mixed_axes.sh: the rows named in stands_in carry the margin and hop of the cell of two different chains "
             "(IVRVVRVCVVVVVR+IVRVVRVCVVRVVR, which belt4 has not built): their dipole, switch strain and stiffness are the named row's own")
json.dump(d, open(p, "w"), indent=1)
PY
