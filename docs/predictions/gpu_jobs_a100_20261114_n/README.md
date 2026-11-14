# Bands in the per-chain metric on an A100, 2026-11-14

One rented pod (Runpod, Secure Cloud, one NVIDIA A100-SXM4-80GB at $2.19 per hour) ran
`examples/runpod_jobs_pod.sh` on `examples/gpu_jobs/queue_20261114_n.txt`: 03:29 to 04:04 UTC, 35 minutes, about $1.29; deleted at the end and the
deletion confirmed by a read.  Every cell is sixteen chains; the rows and the states M and N are
the stored ones, copied in by the jobs.

`STATUS.txt`: each job's start (with its commit) and exit code.  `NN_<script>_log.txt`: a job's
output (`/srv/bench/<token>/` is the pod's served directory, the token removed).  The states and
bands are in `docs/predictions/kink_states_belt4.json`, each row's way to N in its `to_N`
(`docs/predictions/kink_migration_belt4.json`), both tabulated in
`docs/predictions/kink_summary_belt4.md`.

| job | what | result |
|---|---|---|
| 1 | 64/29/7, the band from I to N continued in the metric (it stood at +12.50 after 600 iterations in per-variable units) | restructured in 300 iterations; highest saddle settled at +11.519; two dips behind it |
| 2 | 64/29/7, the direct band from I to B from nothing in the metric | converged in 106 iterations (242 in per-variable units), top +15.307, settled +15.293: a saddle, and still the higher of the two ways (through M: +14.49) |
| 3 | 71/29, I to N continued | converged, top +13.81 on the image next to I; the search from there slid back |
| 4 | PVDF, I to N continued | runs through B: its highest saddle is the step's own, +11.637 |
| 5, 6 | 64/29/7, the two dips relaxed | P is the sheet itself (+0.0000); Q is an N-like state at +10.333 |
| 7 to 9 | 64/29/7 in pieces: I to the sheet, the sheet to Q, Q to N | saddles +9.109 (the first saddle), +10.626 and +11.396 |
| 10 to 12 | 71/29: the dip (P, N-like, +11.167), I to P, P to N | saddles +13.786 and a band top of +12.35 (its search not converged) |

What it says is in `docs/CANDIDATE_MATERIALS.md` 8.3 and `docs/PERFORMANCE_REVIEW.md` 17.10.
