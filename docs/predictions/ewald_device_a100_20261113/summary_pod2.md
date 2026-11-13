Device: NVIDIA A100-SXM4-80GB, 85.1 GB; cupy 14.2.0; host AMD EPYC 7742 64-Core Processor, 13 cores allotted (the container sees 128); NumPy 2.2.6.

One whole evaluation (energy and gradient):

| chains | sites | k vectors | host, 1 thread | host, 8 threads | host, 13 threads | device | of it on the host | device against fastest host | largest difference | card memory |
|---|---|---|---|---|---|---|---|---|---|---|
| 16 | 960 | 14,486 | 15.8 s | 2.80 s | 2.73 s | 0.51 s | 0.36 s (70%) | 5.3x | 7.3e-15 | 1.3 GB |
| 32 | 1,920 | 28,920 | 23.2 s | 11.7 s | 10.8 s | 1.14 s | 0.77 s (68%) | 9.5x | 8.3e-15 | 5.1 GB |

The three ported steps, timed apart (fastest host run against the device):

| chains | tensor and field | solve | gradient | the three together | device against host |
|---|---|---|---|---|---|
| 16 | 1.38 s / 20.7 ms | 0.44 s / 17.3 ms | 0.38 s / 9.6 ms | 2.20 s / 47.6 ms | 46x |
| 32 | 5.69 s / 92.3 ms | 3.04 s / 40.4 ms | 1.21 s / 19.6 ms | 9.94 s / 0.15 s | 65x |

Checks that need no host run, and the first call:

| chains | tiling against the stored two-chain cell | gradient against the device's own energy difference | first device evaluation | device memory in use |
|---|---|---|---|---|
| 16 | 3.4e-13 kcal/mol | 2.2e-05 | 0.89 s | 1.8 of 85.1 GB |
| 32 | 6.8e-13 kcal/mol | 2.4e-06 | 1.36 s | 5.6 of 85.1 GB |

Side measurements (not used by the port):

| chains | general product / rank-k update | double-precision rate of the product | LU solve / Cholesky solve | unknowns |
|---|---|---|---|---|
| 16 | 1.8 ms / 1.6 ms | 14.5 TFLOP/s | 16.9 ms / 3.1 ms | 2,880 |
| 32 | 11.1 ms / 8.5 ms | 19.3 TFLOP/s | 39.5 ms / 8.6 ms | 5,760 |
