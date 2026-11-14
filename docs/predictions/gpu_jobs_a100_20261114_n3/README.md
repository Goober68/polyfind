# What N does next, on an A100, 2026-11-14

One rented pod (Runpod, Secure Cloud, one NVIDIA A100-SXM4-80GB at $2.19 per hour) ran
`examples/runpod_jobs_pod.sh` on `examples/gpu_jobs/queue_20261114_n3.txt`: 05:02 to 05:46 UTC, 43 minutes, about $1.57; deleted at the end and the
deletion confirmed by a read.  Sixteen chains, bands in the per-chain metric; the rows and the
earlier states are the stored ones, copied in by the jobs.

`STATUS.txt`: each job's start (with its commit) and exit code.  `NN_<script>_log.txt`: a job's
output (`/srv/bench/<token>/` is the pod's served directory, the token removed).  The states and
bands are in `docs/predictions/kink_states_belt4.json`, the rows' summaries in their `to_N.then`
and `to_L` (`docs/predictions/kink_migration_belt4.json`), tabulated in
`docs/predictions/kink_summary_belt4.md`.

| job | what | result |
|---|---|---|
| 1 to 4 | 64/29/7 from N: N2, B2, NN4, NN12 relaxed | +14.884, +12.398, +19.896, +19.991 |
| 5 to 7 | 64/29/7: bands from N to N2 and to the two pairs | tops +16.04, +21.34, +20.51 (all converged) |
| 8 to 11 | 71/29 from N: N2, B2, NN4, NN12 | +12.660, +12.591, +18.730, +18.730 |
| 12 to 14, 16 | 71/29: the same bands | N-N2 not converged in 600 iterations (top +14.59); N-NN4 not converged (+22.15); N-NN12 converged (+23.05) |
| 15 | 71/29: the dip under N on the band towards N2, relaxed | L: gauche bonds 3 and 5, +3.744 |

What it says is in `docs/CANDIDATE_MATERIALS.md` 8.3 ("What N does next").
