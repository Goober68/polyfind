#!/usr/bin/env bash
# The film composition's low all-trans cell of 3.23 itself as the start (docs/MLIP_SECOND_OPINION.md 3.24): its
# chains do not have the upright cell's conformation (the small twists of an all-trans chain's bonds fall in another
# pattern), and a cell built from upright chains does not find it.  Four cells:
#   the low cell with its chlorine made fluorine: the matrix, in the film composition's low cell;
#   the low cell with its neighbours moved back to the upright registry (the chains and the turn kept): is it the registry?
#   the low cell given back the upright cell's turn, lattice and registry (the chains kept): is it the chains?
#   the matrix made of the low cell, moved back to the upright registry.
#
# and, in a second file, the low cell with its CFE unit made an acrylonitrile unit (the nitrile where the chlorine was, which
# is the hand the labelled nitrile cell has): the nitrile cell, in the film composition's low cell, and the same moved back
# to the upright registry.
#
# and, in a third, the low cell with every backbone carbon's substituents made the leader's (a reversed VDF unit, the
# nitrile, and a VDF unit where the film's composition has its fourth TrFE): the leader, in the film composition's low cell,
# and the same moved back to the upright registry.
#
#   bash examples/registry_starts_as_film.sh [OUT [OUT2 [OUT3]]]
#       (default docs/predictions/mlip/registry_starts_as_film.extxyz, ..._nitrile.extxyz and ..._leader.extxyz)
set -euo pipefail
cd "$(dirname "$0")/.."
OUT="${1:-docs/predictions/mlip/registry_starts_as_film.extxyz}"
OUT2="${2:-docs/predictions/mlip/registry_starts_as_film_nitrile.extxyz}"
OUT3="${3:-docs/predictions/mlip/registry_starts_as_film_leader.extxyz}"
UP=row_vdf9-trfe4-cfe1-14_TTTT_chain1_in_all_s1
python -I examples/registry_cells.py derive --labels docs/predictions/mlip/cells14_frames_pbe_d3bj.extxyz docs/predictions/mlip/paths/all_trans_tight.extxyz \
  --out "$OUT" --cell film:tt=$UP matrix:tt=row_vdf10-trfe4-14_TTTT_chain1_in_all_s1 --make \
  reg_matrix_tt_as_film_low=film_r_released_tight:matrix:tt:swap.Cl.F.1.36 \
  reg_film_tt_low_registry_back=film_r_released_tight:film:tt:shear.-1.1 \
  reg_film_tt_low_chains_in_upright=film_r_released_tight:film:tt:back.film_r_released_tight.$UP \
  reg_matrix_tt_as_film_low_registry_back=film_r_released_tight:matrix:tt:swap.Cl.F.1.36:shear.-1.1
python -I examples/registry_cells.py derive --labels docs/predictions/mlip/cells14_frames_pbe_d3bj.extxyz docs/predictions/mlip/paths/all_trans_tight.extxyz \
  --out "$OUT2" --cell nitrile:tt=tt_nitrile_from_rows14an_r1 --make \
  reg_nitrile_tt_as_film_low=film_r_released_tight:nitrile:tt:nitrile.Cl.F \
  reg_nitrile_tt_as_film_low_registry_back=film_r_released_tight:nitrile:tt:nitrile.Cl.F:shear.-1.1
python -I examples/registry_cells.py derive --labels docs/predictions/mlip/cells14_frames_pbe_d3bj.extxyz docs/predictions/mlip/paths/all_trans_tight.extxyz \
  --out "$OUT3" --cell leader:tt=tt_leader_from_rows14an_r1 --make \
  reg_leader_tt_as_film_low=film_r_released_tight:leader:tt:morph.tt_leader_from_rows14an_r1 \
  reg_leader_tt_as_film_low_registry_back=film_r_released_tight:leader:tt:morph.tt_leader_from_rows14an_r1:shear.-1.1
