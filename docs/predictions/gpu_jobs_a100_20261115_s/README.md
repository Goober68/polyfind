# One chain's moves in the unit screen's leaders, on an A100, 2026-11-15

One rented pod (Runpod, Secure Cloud, one NVIDIA A100-SXM4-80GB at $2.19 per hour) ran
`examples/runpod_jobs_pod.sh` on `examples/gpu_jobs/queue_20261114_s.txt`: 23:59 UTC on 2026-11-14
to 00:41 UTC on 2026-11-15, 42 minutes, about $1.53; deleted at the end and the deletion confirmed
by a read.  Sixteen chains, one moving, neighbours flexible, cell fixed; bands in the per-chain
metric.  Two lattices: the 64/29/7 cell with acrylonitrile in place of CFE
(`vdf9-trfe4-an1-14`) and the cell with one CTFE and one reversed VDF unit per fourteen
(`vdf8-trfe4-ctfe1-rvdf1-14`), both from the stored polished kinked cells at kink offset 0 (a
queue line cannot set another offset).

`STATUS.txt`: each job's start (with its commit) and exit code.  `NN_<script>_log.txt`: a job's
output (`/srv/bench/<token>/` is the pod's served directory, the token removed).
`NN_band_row_result.json`: the row's end states, without the bands.  The rows are in
`docs/predictions/kink_migration_belt4.json` (with their `to_N` and `to_L`), the states N and L and
the bands A-N and N-L in `docs/predictions/kink_states_belt4.json`, tabulated in
`docs/predictions/kink_summary_belt4.md`.

Energies in kcal/mol against the sheet.

| job | what | result |
|---|---|---|
| 1 | AN cell: the wide kink I, the completed step B, the bands | I +8.869, B +11.859; A-I top +8.925 (settled +8.928); I-B top +14.977, converged; the search for its settled saddle slid off the path (to +7.85, not converged), so the top is the number |
| 2, 3 | AN cell: N (gauche bonds 4 and 6, the other hand), and the band from the sheet to it | N +7.950; top +10.572, settled saddle **+10.562** (a shallow minimum on the way, image 3) |
| 4 | CTFE cell with a reversed unit: I, B, the bands | **I is not a state: it relaxes back to the sheet** (+0.000); B +11.029; the "I-B" band is then a straight path from the sheet to B across cis at bond 2: 900 iterations, not converged (tops +115.6, +158.8, +145.9), and the search from its top ended at +14.75 with positive curvature along the path. No path to B was found. About 17 of the pod's minutes |
| 5, 6 | the same cell: N, and the band from the sheet to it | N +7.248; top +9.824, settled saddle **+9.826** (a shallow minimum on the way, image 5) |
| 7, 8 | AN cell: L (gauche bonds 3 and 5, the sheet's hand), and the band from N to it | L +4.448; top +10.002, settled +10.001: 2.05 above N |
| 9, 10 | CTFE cell with a reversed unit: L, and the band from N to it | L +3.405; top +10.192, settled +10.185: 2.94 above N |

What it says is in `docs/CANDIDATE_MATERIALS.md` 5.6 ("The leaders by the other rate model") and
8.3 ("The same moves in the unit screen's leaders").
