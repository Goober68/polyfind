# The large-cell electrostatics on an A100, 2026-11-13

Two rented pods (Runpod, Secure Cloud, one NVIDIA A100-SXM4-80GB each at $2.19 per hour, driver
595.91.07, image `runpod/pytorch:2.2.0-py3.10-cuda12.1.1-devel-ubuntu22.04`, CuPy 14.2.0 from the
`cupy-cuda13x[ctk]` wheels, NumPy 2.2.6, host AMD EPYC 7742) ran `examples/runpod_ewald_pod.sh`,
which runs `examples/runpod_ewald_bench.sh`: the device tests, then `examples/ewald_device_benchmark.py`
on tilings of the stored two-chain kinked cell (`vdf10@0`, the belt4 preset).

| file | what it is |
|---|---|
| `pod1_16chains.json`, `pod1_log.txt` | first pod, commit 2fe1c9e, 21:22 to 21:25 UTC. Its "all cores" host run used 256 BLAS threads under a 32-core allotment and is not a valid host time (slower than one thread); the pod was removed for that and the benchmark changed (73b556f). Its one-thread host run and its device run stand. |
| `pod2_16chains.json`, `pod2_32chains.json`, `pod2_log.txt` | second pod, commit 73b556f, 21:27 to 21:55 UTC, 13 cores allotted. Complete for 16 and 32 chains. Stopped by the PI during the 64-chain host baselines, so 64 chains has two host lines in the log and no device row. |
| `summary_pod2.md` | `python examples/ewald_device_summary.py pod2_16chains.json pod2_32chains.json` |

What it measured, in one paragraph: the three ported steps (tensor and field, the dipole solve, the
gradient contraction) run 46 times (16 chains) and 65 times (32 chains) faster on the card than on
the 13 host cores, agreeing with the host to 8e-15 relative; but a whole evaluation is only 5.3 and
9.5 times faster, because 70% of the device evaluation is the part that was still on the host (0.36 s
of 0.51 s at 16 chains, 0.77 s of 1.14 s at 32).  That part is what the next port moves.

Notes on reading the files:

* The host times on these pods are noisy.  The container sees 128 or 256 CPUs and is held to its
  allotment by a quota; at one BLAS thread the 16-chain gradient step took 9 to 11 s on the second
  pod (0.46 s on the first, same CPU model) and went back to normal at 8 threads and at 32 chains.
  The comparison uses the fastest host run of each size.
* 64 chains, host only, from `pod2_log.txt`: 462.6 s at one thread and 59.4 s at eight for a whole
  evaluation (3,840 sites, 57,886 k vectors).
* The device tests passed on both pods (17 passed).
* Cost of the two pods together: about $1.14 (3 minutes and 28 minutes at $2.19 per hour).
