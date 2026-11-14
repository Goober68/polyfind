# A queue of jobs on an A100, 2026-11-13 to 14

One rented pod (Runpod, Secure Cloud, one NVIDIA A100-SXM4-80GB at $2.19 per hour, driver
595.91.07, image `runpod/pytorch:2.2.0-py3.10-cuda12.1.1-devel-ubuntu22.04`, CuPy 14.2.0 from the
`cupy-cuda12x` wheel, host AMD EPYC 7742 with 13 cores allotted of 128 visible) ran
`examples/runpod_jobs_pod.sh` from 23:25 to 00:28 UTC: it installed once (28 s) and then ran the
lines of `examples/gpu_jobs/queue_20261113.txt` one after another, each at the branch's commit
of the moment, so that a job could be decided on the result of the one before it by pushing a
line.  63 minutes, about $2.30; the pod was deleted at the end and the deletion confirmed by a read.

| file | what it is |
|---|---|
| `STATUS.txt` | one line per event: each job's start (with its commit) and exit code |
| `setup_log.txt` | the installation |
| `NN_<script>_log.txt` | a job's output; `/srv/bench/<token>/` in a path is the pod's served directory, the token removed |
| `NN_step_wide_result.json` | the end states a `step_wide.sh` job relaxed, without their state vectors |
| `08_device_check_*.json`, `08_device_check_summary.md` | the device check (73 device tests passed, then `examples/ewald_device_benchmark.py`): result files and the tables `examples/ewald_device_summary.py` makes from them |

The state vectors and the saddles are not here: they are rows of
`docs/predictions/kink_migration_belt4.json` (`examples/kink_wide_cell.py store`), tabulated in
`docs/predictions/kink_summary_belt4.md`.

| job | what | result |
|---|---|---|
| 1 | device check | failed two tests that assert the default array module is NumPy, because the pod exports `POLYFIND_EWALD_DEVICE=cuda` to its jobs; the script now clears it for the tests (job 8) |
| 2 to 4 | the thirty-two-chain row (8a x 2b) again, from the same seeds as the host's | I +9.9740, B +5.6921, saddles +10.6601 and +11.4726: the host's to 5e-9, 3e-8, 2e-8 and 5e-5 kcal/mol; 0.25 s an evaluation against 19.7 s |
| 5 to 7 | sixty-four chains (8a x 4b), from the thirty-two-chain row let in along b | I +9.8875, B +5.4950, saddles +10.6283 and +11.5330 |
| 8 | device check | one evaluation 52 ms, 0.13 s, 0.57 s and 3.67 s at 16, 32, 64 and 128 chains with the tensor's products as rank-k updates; one to four processes on the card 18.7, 23.0, 28.5, 28.4 evaluations per second at 16 chains and 7.4, 7.6, 8.6, 8.6 at 32 |
| 9, 10 | the completed step in a hundred and twenty-eight chains | B +5.4305 as 8a x 8b, +5.4356 as 16a x 4b; 7.5 s an evaluation |
| 11 to 13 | the other thirty-two-chain cell (4a x 4b), whose wide kink and saddles had not been computed | I +9.8998, B +5.7122 (the host's: +5.7122), saddles +10.6353 and +11.6200 |
| 14 | the sixty-four-chain second saddle by the other route, through 4a x 4b | +11.5330 again, 1e-7 from job 7's |

What it says is in `docs/PERFORMANCE_REVIEW.md` 17.9 (the device) and `docs/CANDIDATE_MATERIALS.md`
8.3 (the step against the size of the cell).
