# The copolymers' single-chain step on an A100, 2026-11-14

One rented pod (Runpod, Secure Cloud, one NVIDIA A100-SXM4-80GB at $2.19 per hour, image
`runpod/pytorch:2.2.0-py3.10-cuda12.1.1-devel-ubuntu22.04`, CuPy 14.2.0) ran
`examples/runpod_jobs_pod.sh` on `examples/gpu_jobs/queue_20261114.txt` from 01:29 to 03:07 UTC:
98 minutes, about $3.60; deleted at the end and the deletion confirmed by a read.  A first pod
was deleted 98 seconds in ($0.06): its second job was the four-chain PVDF row, 240 charged sites,
which the code evaluates on the host whatever the machine has, and a rented card would have sat
idle for 25 minutes.

| file | what it is |
|---|---|
| `STATUS.txt` | one line per event: each job's start (with its commit) and exit code |
| `NN_<script>_log.txt` | a job's output; `/srv/bench/<token>/` in a path is the pod's served directory, the token removed |
| `01_forms_check_forms_*.json` | the batched form on the card against the host's loops on the rows' own cells |
| `NN_band_row_result.json` | a row as a `band_row.sh` job left it, without state vectors and band images |

The rows are in `docs/predictions/kink_migration_belt4.json`, the states M and N and the bands
through them in `docs/predictions/kink_states_belt4.json`, both tabulated in
`docs/predictions/kink_summary_belt4.md`.

| job | what | result |
|---|---|---|
| 1 | the evaluation's two forms on the terpolymer's, the matrix's and PVDF's cells of four and sixteen chains | agree to rounding (energy 1e-13 kcal/mol at a physical displacement, 1e-14 of its size at an unphysical one) |
| 2, 7 | 64/29/7, one chain in four | I +9.045, B +12.999, saddles +9.084 and +13.307; the second band needed 1101 iterations to find the dip at M |
| 3, 8 | 64/29/7, one chain in sixteen | I +9.053, B +11.795, first saddle +9.109; the second band converged with its top at +15.41 and the saddle search from it left the path, twice |
| 4 | 71/29, one chain in four | I +9.525, B +12.311, saddles +9.583 and +13.675 |
| 5 | 71/29, one chain in sixteen | I +9.355, B +10.703, saddles +9.526 and +15.000 |
| 6 | PVDF, one chain in sixteen, from nothing: the control | I +9.9626, B +5.9110, saddles +10.6600 and +11.6365 (stored host row: +9.9626, +5.9111, +10.6600, +11.6366); 104 and 88 s for the bands the host took 4976 and 4803 s for |
| 9, 10 | 64/29/7 in four chains through M (+11.204) | saddles +12.999 and +13.307: the direct band's, by another route |
| 11 to 13 | 64/29/7 in sixteen chains through M (+13.089) | saddles +14.491 and +13.982: the row's barrier |
| 14 | where the sixteen-chain search had ended, relaxed | N, gauche bonds 4 and 6, +8.820 |
| 15 to 18 | 71/29: M in four chains (+11.680) and sixteen (+13.432), the bands through it | saddles +15.000 (the direct band's) and +13.567 |
| 19, 24 | 64/29/7 in sixteen chains, I to N | band not converged in 600 iterations: top +12.50 |
| 20 to 23 | N in 71/29 (+8.944) and PVDF (+7.722), and the bands to it | not converged in 300 iterations: tops +14.82 and +13.32 |

What it says is in `docs/CANDIDATE_MATERIALS.md` 8.3 ("One chain alone in the copolymer lattices").
