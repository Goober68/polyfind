Device: NVIDIA A100-SXM4-80GB, 85.1 GB; cupy 14.2.0; host Intel(R) Xeon(R) Platinum 8470, 20 cores allotted (the container sees 192); NumPy 2.2.6.

One whole evaluation (energy and gradient):

| chains | sites | k vectors | host, 20 threads | device | of it outside the dipole step | device against fastest host | largest difference | card memory |
|---|---|---|---|---|---|---|---|---|
| 16 | 960 | 14,486 | 2.88 s | 48.8 ms | 13.8 ms (28%) | 59.2x | 8.6e-15 | 1.4 GB |
| 32 | 1,920 | 28,920 | - | 0.14 s | 20.7 ms (14%) | - | no host run | 5.5 GB |
| 64 | 3,840 | 57,886 | - | 0.76 s | 36.4 ms (5%) | - | no host run | 21.1 GB |
| 128 | 7,680 | 115,706 | - | 6.00 s | 73.8 ms (1%) | - | no host run | 82.5 GB |

The three ported steps, timed apart (fastest host run against the device):

| chains | tensor and field | solve | gradient | the three together | device against host |
|---|---|---|---|---|---|
| 16 | 1.39 s / 19.4 ms | 0.40 s / 3.4 ms | 0.58 s / 8.9 ms | 2.38 s / 31.8 ms | 75x |
| 32 | - / 88.8 ms | - / 9.4 ms | - / 18.9 ms | - / 0.12 s | - |
| 64 | - / 0.61 s | - / 45.6 ms | - / 58.7 ms | - / 0.71 s | - |
| 128 | - / 5.40 s | - / 0.28 s | - / 0.21 s | - / 5.89 s | - |

One device evaluation part by part (the device waited for at the end of each part; build, valence, exclusions and assembly are the host's):

| chains | build | valence | pair sums | electrostatics | exclusions | assembly | the host's parts | of the evaluation |
|---|---|---|---|---|---|---|---|---|
| 16 | 6.0 ms | 0.0 ms | 4.4 ms | 34.7 ms | 0.9 ms | 2.3 ms | 9.2 ms | 19% |
| 32 | 9.4 ms | 0.0 ms | 5.9 ms | 0.12 s | 1.7 ms | 3.7 ms | 14.8 ms | 10% |
| 64 | 17.9 ms | 0.0 ms | 10.7 ms | 0.72 s | 3.6 ms | 6.3 ms | 27.8 ms | 4% |
| 128 | 34.9 ms | 0.0 ms | 23.3 ms | 5.90 s | 8.2 ms | 12.3 ms | 55.3 ms | 1% |

Evaluations back to back, each at a new geometry, and the card's own utilisation figure meanwhile:

| chains | form | one evaluation | evaluations per second | card utilisation, mean (lowest to highest sample) |
|---|---|---|---|---|
| 16 | batched: pair sums and electrostatics on the card | 49.5 ms | 20.2 | 49% (47 to 50%) |
| 16 | loops: electrostatics on the card, the rest in host loops | 0.50 s | 2.0 | 5% (0 to 13%) |
| 32 | batched: pair sums and electrostatics on the card | 0.15 s | 6.8 | 77% (67 to 84%) |
| 32 | loops: electrostatics on the card, the rest in host loops | 1.11 s | 0.9 | 10% (0 to 62%) |
| 64 | batched: pair sums and electrostatics on the card | 0.76 s | 1.3 | 91% (0 to 100%) |
| 128 | batched: pair sums and electrostatics on the card | 5.96 s | 0.2 | 92% (0 to 100%) |
| 128 | loops: electrostatics on the card, the rest in host loops | 10.9 s | 0.1 | 56% (0 to 100%) |

The batched form on the host (same BLAS threads as the fastest loop run):

| chains | loops | batched | difference in energy | in the gradient |
|---|---|---|---|---|
| 16 | 2.88 s | 2.41 s | 1.3e-15 | 6.0e-15 |

Checks that need no host run, and the first call:

| chains | tiling against the stored two-chain cell | gradient against the device's own energy difference | first device evaluation | device memory in use |
|---|---|---|---|---|
| 16 | 2.3e-13 kcal/mol | 2.2e-05 | 0.38 s | 1.9 of 85.1 GB |
| 32 | 2.3e-13 kcal/mol | 2.4e-06 | 0.39 s | 5.9 of 85.1 GB |
| 64 | 2.3e-13 kcal/mol | 8.3e-06 | 0.85 s | 21.6 of 85.1 GB |
| 128 | 6.8e-13 kcal/mol | 2.7e-05 | 8.00 s | 83.0 of 85.1 GB |

Side measurements (not used by the port):

| chains | general product / rank-k update | double-precision rate of the product | LU solve / Cholesky solve | unknowns |
|---|---|---|---|---|
| 16 | 1.8 ms / 1.7 ms | 14.8 TFLOP/s | 16.7 ms / 3.1 ms | 2,880 |
| 32 | 11.0 ms / 8.3 ms | 19.3 TFLOP/s | 38.9 ms / 8.6 ms | 5,760 |
| 64 | 88.1 ms / 51.1 ms | 19.4 TFLOP/s | 0.14 s / 42.7 ms | 11,520 |
| 128 | 0.85 s / 0.43 s | 16.1 TFLOP/s | 0.66 s / 0.27 s | 23,040 |
