Device: NVIDIA A100-SXM4-80GB, 85.1 GB; cupy 14.2.0; host AMD EPYC 7742 64-Core Processor, 13 cores allotted (the container sees 128); NumPy 2.2.6.

One whole evaluation (energy and gradient):

| chains | sites | k vectors | device | of it outside the dipole step | device against fastest host | largest difference | card memory |
|---|---|---|---|---|---|---|---|
| 16 | 960 | 14,486 | 52.1 ms | 16.3 ms (31%) | - | no host run | 1.4 GB |
| 32 | 1,920 | 28,920 | 0.13 s | 25.3 ms (19%) | - | no host run | 5.5 GB |
| 64 | 3,840 | 57,886 | 0.57 s | 43.5 ms (8%) | - | no host run | 21.1 GB |
| 128 | 7,680 | 115,706 | 3.67 s | 85.6 ms (2%) | - | no host run | 82.5 GB |

The three ported steps, timed apart (fastest host run against the device):

| chains | tensor and field | solve | gradient | the three together | device against host |
|---|---|---|---|---|---|
| 16 | - / 19.9 ms | - / 3.5 ms | - / 9.1 ms | - / 32.5 ms | - |
| 32 | - / 70.9 ms | - / 9.9 ms | - / 19.2 ms | - / 100.0 ms | - |
| 64 | - / 0.41 s | - / 46.1 ms | - / 58.2 ms | - / 0.51 s | - |
| 128 | - / 3.06 s | - / 0.29 s | - / 0.21 s | - / 3.57 s | - |

One device evaluation part by part (the device waited for at the end of each part; build, valence, exclusions and assembly are the host's):

| chains | build | valence | pair sums | electrostatics | exclusions | assembly | the host's parts | of the evaluation |
|---|---|---|---|---|---|---|---|---|
| 16 | 8.0 ms | 0.0 ms | 4.7 ms | 33.7 ms | 1.0 ms | 2.6 ms | 11.5 ms | 23% |
| 32 | 12.3 ms | 0.0 ms | 6.1 ms | 0.11 s | 1.8 ms | 4.3 ms | 18.5 ms | 14% |
| 64 | 24.9 ms | 0.0 ms | 10.7 ms | 0.52 s | 3.6 ms | 7.8 ms | 36.3 ms | 6% |
| 128 | 41.0 ms | 0.0 ms | 22.9 ms | 3.59 s | 7.4 ms | 14.8 ms | 63.1 ms | 2% |

Evaluations back to back, each at a new geometry, and the card's own utilisation figure meanwhile:

| chains | form | one evaluation | evaluations per second | card utilisation, mean (lowest to highest sample) |
|---|---|---|---|---|
| 16 | batched: pair sums and electrostatics on the card | 53.4 ms | 18.7 | 45% (42 to 48%) |
| 32 | batched: pair sums and electrostatics on the card | 0.14 s | 7.4 | 71% (62 to 80%) |
| 64 | batched: pair sums and electrostatics on the card | 0.57 s | 1.8 | 89% (70 to 100%) |
| 128 | batched: pair sums and electrostatics on the card | 3.63 s | 0.3 | 89% (0 to 100%) |

Several processes on the one card, each with its own cell, evaluating back to back over a common window:

| chains | processes | evaluations per second, all together | one evaluation, per process | card utilisation |
|---|---|---|---|---|
| 16 | 1 | 18.7 | 53.6 ms | 46% |
| 16 | 2 | 23.0 | 86.9 ms | 79% |
| 16 | 3 | 28.5 | 0.11 s | 98% |
| 16 | 4 | 28.4 | 0.14 s | 100% |
| 32 | 1 | 7.4 | 0.14 s | 72% |
| 32 | 2 | 7.6 | 0.26 s | 89% |
| 32 | 3 | 8.6 | 0.35 s | 100% |
| 32 | 4 | 8.6 | 0.47 s | 100% |

Checks that need no host run, and the first call:

| chains | tiling against the stored two-chain cell | gradient against the device's own energy difference | first device evaluation | device memory in use |
|---|---|---|---|---|
| 16 | 2.3e-13 kcal/mol | 2.2e-05 | 0.55 s | 2.0 of 85.1 GB |
| 32 | 4.5e-13 kcal/mol | 2.4e-06 | 0.16 s | 6.0 of 85.1 GB |
| 64 | 6.8e-13 kcal/mol | 8.3e-06 | 0.65 s | 21.6 of 85.1 GB |
| 128 | 6.8e-13 kcal/mol | 2.7e-05 | 7.66 s | 83.0 of 85.1 GB |

Side measurements (not used by the port):

| chains | general product / rank-k update | double-precision rate of the product | LU solve / Cholesky solve | unknowns |
|---|---|---|---|---|
| 16 | 1.5 ms / 1.4 ms | 17.5 TFLOP/s | 12.5 ms / 2.6 ms | 2,880 |
| 32 | 11.2 ms / 7.7 ms | 19.1 TFLOP/s | 35.8 ms / 8.7 ms | 5,760 |
| 64 | 94.2 ms / 55.2 ms | 18.1 TFLOP/s | 0.13 s / 43.1 ms | 11,520 |
| 128 | 0.87 s / 0.46 s | 15.6 TFLOP/s | 0.67 s / 0.28 s | 23,040 |
