#!/usr/bin/env bash
# The four axes of the screen's two rows (docs/MLIP_SECOND_OPINION.md 3.25): strain, tan delta, frequency and blocked
# force per area, each row's margin and its kink's hop put through the tables of docs/CANDIDATE_MATERIALS.md 5.6
# (examples/four_axes_at_margins.py --barriers).
#
#   bash examples/screen14_axes.sh [OUTDIR]        (default docs/predictions)
#
# A new row has no dipole, switch strain or stiffness of its own (belt4 has not built it), so each is put on the
# inputs of the two rows nearest it and the table is given twice:
#   IVRVVRVCVVVVVR (margin: DFT, from its lowest labelled cells; screen14_dft.json): on IVRVVRVCVVRVVR's inputs (a TrFE
#     eleventh for its VDF) and on the leader's (a nitrile eighth for its CFE);
#   IVRVVRVVVVVVVR (margin: the potential's; screen14_potential_d.json; it has no DFT): on the leader's inputs (a
#     nitrile eighth for its VDF) and on the matrix's (no reversed unit, a TrFE eleventh).
# The frequency is the measured film's 364 kHz times exp(-(hop - the film composition's hop) / kT) with every hop in
# the potential (the film composition's 7.81, 3.23; DFT has 7.70 for it and has not been asked about the new rows').
# Each margin is also moved by 0.05 either way: what the potential's and DFT's margins differ by on the cells that
# have both.  The film's composition is in every table at its DFT margin from its lowest labelled cells (+0.063).
set -euo pipefail
cd "$(dirname "$0")/.."
OUT="${1:-docs/predictions}"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
python - "$OUT" "$TMP" <<'PY'
import json, sys
out, tmp = sys.argv[1:3]
dft = json.load(open(f"{out}/screen14_dft.json"))["rows"]["vdf9-trfe3-cfe1-rvdf1-14"]["margin"]
pot = json.load(open(f"{out}/screen14_potential_d.json"))["rows"]["vdf10-trfe3-rvdf1-14"]["margin"]
hops = json.load(open(f"{out}/screen14_hops.json"))["paths"]
film = 0.063
json.dump({"rows": {"vdf8-trfe4-cfe1-rvdf1-14": {"margin": dft}, "vdf9-trfe3-an1-rvdf1-14": {"margin": pot}, "vdf9-trfe4-cfe1-14": {"margin": film}}}, open(f"{tmp}/a.json", "w"))
json.dump({"rows": {"vdf9-trfe3-an1-rvdf1-14": {"margin": dft}, "vdf10-trfe4-14": {"margin": pot}, "vdf9-trfe4-cfe1-14": {"margin": film}}}, open(f"{tmp}/b.json", "w"))
open(f"{tmp}/hops.txt", "w").write(f"{hops['lcfe']['hop']:.4f} {hops['nopin']['hop']:.4f}\n")
PY
read -r CFE NOPIN < "$TMP/hops.txt"
FILM=7.81
python examples/four_axes_at_margins.py --margins-from "$TMP/a.json" --shifts=-0.05,0,0.05 \
  --barriers vdf9-trfe4-cfe1-14=$FILM vdf8-trfe4-cfe1-rvdf1-14=$CFE vdf9-trfe3-an1-rvdf1-14=$NOPIN --out "$OUT/four_axes_screen14_a.json"
python examples/four_axes_at_margins.py --margins-from "$TMP/b.json" --shifts=-0.05,0,0.05 \
  --barriers vdf9-trfe4-cfe1-14=$FILM vdf9-trfe3-an1-rvdf1-14=$CFE vdf10-trfe4-14=$NOPIN --out "$OUT/four_axes_screen14_b.json"
python - "$OUT" <<'PY'
import json, sys
out = sys.argv[1]
names = {"a": {"vdf8-trfe4-cfe1-rvdf1-14": "IVRVVRVCVVVVVR on IVRVVRVCVVRVVR's inputs", "vdf9-trfe3-an1-rvdf1-14": "IVRVVRVVVVVVVR on the leader's inputs"},
         "b": {"vdf9-trfe3-an1-rvdf1-14": "IVRVVRVCVVVVVR on the leader's inputs", "vdf10-trfe4-14": "IVRVVRVVVVVVVR on the matrix's inputs"}}
for which, stands in in names.items():
    p = f"{out}/four_axes_screen14_{which}.json"
    d = json.load(open(p))
    d["stands_in"] = stands_in
    d["note"] = ("examples/screen14_axes.sh: the rows named in stands_in carry ANOTHER row's margin and hop (the screen's new rows, which belt4 has not built): "
                 "their dipole, switch strain and stiffness are the named row's own")
    json.dump(d, open(p, "w"), indent=1)
PY
