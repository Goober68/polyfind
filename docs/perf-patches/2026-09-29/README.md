# Verified performance prototypes, 2026-09-29

Each file is the unified diff of one prototype from the second performance review
(`docs/PERFORMANCE_REVIEW.md` section 15), made in an isolated worktree against
branch head `bd5ef12` and measured back to back with the unchanged code.  They are
kept here as patches rather than applied, because other work on this branch depends
on the source staying as it is; every one of them applied cleanly to `bd5ef12` with
`git apply --check` on the day it was made.

| file | change | measured effect | numerics vs current code |
|---|---|---|---|
| `G1.diff` | NeRF builder on component arrays, vectorised pendant placement, batched per-row `repeat_chains` | refinement stage 15.15 -> 8.93 s (1.70x) | bit-identical (coordinates, energies, L-BFGS trace) |
| `G2.diff` | analytic torsion Jacobian of the chain build (`chainjac.py`), closure Newton without re-differencing | refinement stage 14.8 -> 11.9 s (1.24x) alone; chain builds per evaluation 98 -> 3 | gradient to 8.8e-7 of scale; free-parametrisation results move at 1e-5 kcal/mol per monomer |
| `G3.diff` | table-build kernel: buffered in-place passes, pair-type grouping, larger chunk; splined potential as an opt-in that does not pay | table stage 24.05 -> 10.55 s (2.28x, 4 processes) | tables differ by <= 8e-4 kcal/mol at wall-adjacent entries, same accuracy against the exact kernel; serial == parallel still bit-identical |
| `G4.diff` | `fft_screen` vectorised, thread pool removed | screen stage 12.76 -> 2.26 s (5.6x) | bit-identical |
| `G5.diff` | one chain build per candidate, batched helix analysis, vectorised k-best | enumeration 3.95 -> 0.06 s (60x) | exact |
| `G6.diff` | table build on the array backend (CuPy path), GPU tests | CPU path 1.0x, bit-identical; GPU projected only, with four defects to fix first (section 15.4) | see section 15.4 |

Apply one at a time on a fresh checkout of `bd5ef12` (or later) with
`git apply docs/perf-patches/2026-09-29/<id>.diff`; `G3`, `G4` and `G6` all edit
`lattice_table.py` and were made independently, so applying more than one of them
needs a manual merge.  `G1` and `G2` both touch `pack.py` and `linegroup.py` and were
also made independently; their combined effect on the refinement was not measured.
The benchmark scripts and logs behind every number are named in section 15.
