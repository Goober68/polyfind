# Bands in the per-chain metric on an A100, 2026-11-14

One rented pod (Runpod, Secure Cloud, one NVIDIA A100-SXM4-80GB at $2.19 per hour) ran
`examples/runpod_jobs_pod.sh` on `examples/gpu_jobs/queue_20261114_n2.txt`: 03:57 to 04:09 UTC, 11 minutes, about $0.42; deleted at the end and the
deletion confirmed by a read.  Every cell is sixteen chains; the rows and the states M and N are
the stored ones, copied in by the jobs.

`STATUS.txt`: each job's start (with its commit) and exit code.  `NN_<script>_log.txt`: a job's
output (`/srv/bench/<token>/` is the pod's served directory, the token removed).  The states and
bands are in `docs/predictions/kink_states_belt4.json`, each row's way to N in its `to_N`
(`docs/predictions/kink_migration_belt4.json`), both tabulated in
`docs/predictions/kink_summary_belt4.md`.

| job | what | result |
|---|---|---|
| 1 | 64/29/7, from the sheet straight to N | converged in 148 iterations; saddle +11.3960, the same as the piece from Q to N |
| 2 | 71/29, from the sheet straight to N | top settled at +11.893 (the band's plateau behind it not converged in 300 iterations) |
| 3 | PVDF, from the sheet straight to N | converged in 128 iterations; saddle +12.387, above its step's +11.637 |

What it says is in `docs/CANDIDATE_MATERIALS.md` 8.3 and `docs/PERFORMANCE_REVIEW.md` 17.10.
