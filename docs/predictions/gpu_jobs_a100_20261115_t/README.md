# One chain's moves in the nitrile cell with a reversed unit on 21 mol% TrFE, on an A100, 2026-11-15

One rented pod (Runpod, Secure Cloud, one NVIDIA A100-SXM4-80GB at $2.19 per hour) ran
`examples/runpod_jobs_pod.sh` on `examples/gpu_jobs/queue_20261115_t.txt`: 03:36 to 04:01 UTC, 25
minutes, about $0.91; deleted at the end and the deletion confirmed by a read.  The cell is
`vdf9-trfe3-an1-rvdf1-14` (`IVRVVRVAVVVVVR`), sixteen chains, one moving, neighbours flexible, cell
fixed, bands in the per-chain metric.  Two start sites, written `ROW@OFFSET` in the queue: the
sheet at bond offset 18 (a VDF-VDF pair in the five-VDF run, the cheapest site of the G+ T G-
hand, +1.86 against all-trans) and at offset 0 (the pair on the reversed unit, +2.60).

`STATUS.txt`: each job's start (with its commit) and exit code.  `NN_<script>_log.txt`: a job's
output (`/srv/bench/<token>/` is the pod's served directory, the token removed).
`NN_band_row_result.json`: a row's end states, and its bands' summaries where they ran.  The rows
are in `docs/predictions/kink_migration_belt4.json`, the states N and L and the bands A-N and N-L
in `docs/predictions/kink_states_belt4.json`, tabulated in `docs/predictions/kink_summary_belt4.md`.

Energies in kcal/mol against the sheet at the start site.

| job | what | result |
|---|---|---|
| 1 | offset 18: the wide kink I and the completed step B (end states only) | I +8.458, B +6.790 |
| 2, 3 | offset 18: N, and the band from the sheet to it | N +7.027; top +10.219, settled saddle **+10.218** |
| 4, 5 | offset 18: L, and the band from N to it | L +6.118; the band did not converge in 300 iterations (top +11.68, two minima on the way) and the search from its top is not a saddle: no path from N to L here |
| 6 | offset 0: I and B (end states only) | **I is not a state** (it relaxes to the sheet, as at the reversed unit's site of the CTFE cell); B +12.846 |
| 7, 8 | offset 0: N, and the band from the sheet to it | N +7.980; top and settled saddle **+10.660** |
| 9, 10 | offset 0: L, and the band from N to it | L +4.480; top and settled saddle +10.111, 2.13 above N |
| 11 | offset 18: the one-monomer step's two bands and saddles (queued after job 1 showed B as low as PVDF's) | A-I +9.561, settled +9.554; I-B +10.354, settled **+10.350**: 10.35 ahead, 3.56 back |

What it says is in `docs/CANDIDATE_MATERIALS.md` 5.6 ("The reversed unit beside the nitrile, and
beside CFE") and 8.3 ("The same moves in the unit screen's leaders").
