# What N does next, on an A100, 2026-11-14

One rented pod (Runpod, Secure Cloud, one NVIDIA A100-SXM4-80GB at $2.19 per hour) ran
`examples/runpod_jobs_pod.sh` on `examples/gpu_jobs/queue_20261114_n4.txt`: 05:41 to 06:07 UTC, 26 minutes, about $0.94; deleted at the end and the
deletion confirmed by a read.  Sixteen chains, bands in the per-chain metric; the rows and the
earlier states are the stored ones, copied in by the jobs.

`STATUS.txt`: each job's start (with its commit) and exit code.  `NN_<script>_log.txt`: a job's
output (`/srv/bench/<token>/` is the pod's served directory, the token removed).  The states and
bands are in `docs/predictions/kink_states_belt4.json`, the rows' summaries in their `to_N.then`
and `to_L` (`docs/predictions/kink_migration_belt4.json`), tabulated in
`docs/predictions/kink_summary_belt4.md`.

| job | what | result |
|---|---|---|
| 1 to 3 | L built from the sheet in 71/29, 64/29/7 and PVDF | +3.7440 (the dip's, to four decimals), +4.9692, +5.4227 |
| 4 to 6 | from the sheet straight to L | 71/29 top +15.58 (converged; its search slid), 64/29/7 saddle +15.290, PVDF saddle +11.821 |
| 7 to 9 | from N to L | saddles +10.557 (71/29), +10.493 (64/29/7), +13.657 (PVDF) |

What it says is in `docs/CANDIDATE_MATERIALS.md` 8.3 ("What N does next").
