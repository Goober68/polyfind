# polyfind: algorithmic performance review

Scope: the algorithms, not the constant factors. The question for each stage is
whether the work being done is the *right* work, and what formulation would make
the same answer cheap. The numbers the review was *written* against are the
funnel run recorded in DESIGN.md section 5.3 (4-core cloud box): fit 7.6 s,
enumeration 3.6 s, packing of four conformations 81 s, refinement 88 s,
amorphous statistics 0.4 s; the expected gains in sections 1-9 are the reasoned
estimates made at that time, not measurements. **Every "today" number in
sections 1-9 is stale: section 10 carries a measured stage-by-stage profile of
the code as it was after that plan, and section 15 a second one of the current head
with the next round of verified prototypes.** What it says, in one
line: with analytic gradients in the refinement and the polish, no single stage
dominates any more -- refinement is 36% of the funnel and the chain-pair table
build 34%.

Hardware note: the target machine has an AMD Radeon RX 7700 XT. The current
GPU backend is CuPy, which needs CUDA (or ROCm on Linux) and cannot run on
this GPU under Windows. torch + torch-directml does run on it (verified: the
device is visible and executes kernels). Section 9 covers what that implies.

## 1. Packing: replace the direct lattice sum with a tabulated chain-pair interaction

**What the code does now.** `CrystalPacker.energy` evaluates, for every
candidate cell, every atom pair between the primary chains and every periodic
image within the cutoff: cost O(M x I x N^2) with I ~ 100-300 images and N =
12-48 atoms, i.e. 10^5-10^6 pair terms per cell. The search then samples 6,000
random cells and polishes a handful with Nelder-Mead, so the screen covers a
6-dimensional space with 6,000 points; the polyethylene herringbone minimum was
found only as the second-best start.

**The structural fact being ignored.** The chains are rigid and 1D-periodic,
and the energy is a sum of *pair-of-chains* interactions. The interaction
between chain A and an image of chain B depends only on their relative geometry:
lateral separation vector rho, z-shift dz, the two setting angles and the flip.
Rotating the pair as a whole about z leaves the energy unchanged, so the
interaction is a function of four continuous variables

    W(|rho|, phi1 - theta_rho, phi2 - theta_rho, dz)         (+ a flip flag)

where theta_rho is the direction of rho. A cell energy is then a sum of W over
the lattice sites within reach:

    E = 1/2 sum over sites s of L (s != 0) of [ W(rho_s; phi1-theta_s, phi1-theta_s, 0) + (same with phi2) ]
        + sum over sites s of L + (1/2, 1/2) of W(rho_s; phi1-theta_s, phi2-theta_s, dz, flip)
        + E_intra(chain)                                    (a constant per chain)

with 10-20 sites inside the cutoff for realistic cells.

**Proposed change.** Tabulate W once per chain conformation on a grid in
(|rho|, alpha1, alpha2, dz) using the existing pair kernel (exact at the grid
points), then evaluate cell energies as lattice sums of interpolated table
values. Sizes: 130 radii (0.05 A from 3.5 to 10 A) x 72 x 72 angles (5 deg) x
16 z-shifts x 2 flips = 2.2 x 10^7 float32 entries (86 MB), each an N^2 x (2K+1)
pair sum; ~10^10 pair terms for the gamma chain, which is seconds on a GPU and
about a minute on 12 cores, done once per conformation. Exchange symmetry
(swap the two chains) halves it. After that a cell costs ~20 interpolations:
microseconds instead of milliseconds.

**What it buys beyond speed.** With W tabulated, the search can be exhaustive
rather than random. For fixed (a, b, gamma) the energy as a function of
(phi1, phi2, dz) is a lattice sum of rotated copies of W; in Fourier space a
rotation by theta_s is a phase factor, so

    E_hat(m1, m2, n) = sum_s exp(-i (m1 + m2) theta_s) W_hat(|rho_s|; m1, m2, n)

and one inverse FFT (72 x 72 x 16 points) gives E on the whole
(phi1, phi2, dz) grid. A 60 x 60 grid in (a, b) then costs 3,600 small FFTs:
the *entire* 7 x 10^7-point search space in well under a second on the CPU,
with no basin missed. The best few grid minima are polished with the direct
kernel as now. Interpolation error only affects which points are selected for
polishing, never the final energies.

**Expected effect.** Screen: from ~20-60 s to < 1 s per conformation, after a
one-off table build. Coverage: global instead of 6,000 random points.
Effort: medium (table builder reusing `CrystalPacker`, interpolation, FFT
lattice sum, symmetry bookkeeping). Risk: interpolation near the repulsive
wall; mitigated by the fine radial grid and by polishing with the exact kernel.
GPU: the table build is embarrassingly parallel and the natural GPU job.

**Measured (screen grid of `pack.SCREEN_TABLE`, 6 worker processes, cold cache).**
PE all-trans 1.9 s, beta-PVDF 1.6 s, alpha-PVDF 3.7 s, gamma-PVDF 8.2 s; the FFT
screen itself 1.4-1.8 s per conformation for a 10^8-point landscape. On the finer
grid that `PairTable.build` defaults to, PE 18 s and gamma 87 s. See section 10
for the split across worker counts and for what the build was before.

## 2. Polishing and refinement: gradients, batched, instead of Nelder-Mead one cell at a time

**Now.** `polish` runs Nelder-Mead for up to 1,500 sequential single-cell
evaluations per start; `refine_crystal` runs Nelder-Mead in 6 + P dimensions
(up to 14 for the gamma chain) for 1,500-2,500 evaluations, and every
evaluation rebuilds a five-block oligomer twice, redoes the Kabsch fit and
orientation, re-tiles the parameter matrices and rebuilds the scale tensor
before evaluating one cell.

**Change.** The lattice energy is smooth in all variables, so use a
quasi-Newton method (L-BFGS-B or a trust-region method) with gradients. Two
routes:

* *Now, with NumPy:* finite-difference gradients evaluated as **one batched
  kernel call** per iteration (1 + 2 x n_vars configurations, 13 for the cell,
  up to 29 with torsions). Thirty iterations = thirty batched calls instead of
  1,500-2,500 sequential ones. The batched kernel already exists.
* *Later, with torch:* implement the kernel in torch and take exact gradients
  by automatic differentiation, including through the NeRF chain build. The
  same implementation runs on the AMD GPU via DirectML (section 9), so one
  kernel serves batching, gradients and the GPU.

Separately, hoist the invariants: topology, bonded-exclusion scales, LJ/charge
tiling and the image table depend on the chain *topology*, not on the torsion
values, and should be built once per refinement, with only coordinates updated
per evaluation.

**Expected effect.** Polish: 10-50x fewer sequential kernel calls, and a
gradient method converges to a tighter minimum than Nelder-Mead at xatol 1e-3.
Refinement: 20-50x fewer evaluations and roughly 5x cheaper evaluations, so
the 88 s stage becomes a few seconds. Effort: small for batched finite
differences, medium for torch autodiff. Risk: none to results; the minima are
the same or better.

**Measured, and the second route was not needed.** Both stages now run on
*analytic* gradients (`CrystalPacker.energy_and_grad`), and torch was not
involved: the derivatives have closed forms, one kernel row per iteration instead
of `1 + 2 n_vars`, 10.4x fewer rows in the polish and 25.6x in the refinement for
4.9x and 7.1x in time. The minima are the same or lower. Section 13 carries the
derivation, what stayed numerical and why, and the verification; section 10 the
numbers. The autodiff bullet above therefore survives only as the GPU argument of
section 9, not as the gradient argument.

## 3. Refinement variables: impose the chain's line-group symmetry, add bond angles, drop the penalty

**Now.** All P torsions of the repeat are free, and periodicity is enforced by
a quadratic penalty on the residual rotation of the repeat transform. The
penalty makes the objective stiff, the optimiser spends evaluations fighting
it, and bond angles are frozen, which is why the beta chain cannot deflect
(real beta-PVDF has CCC angles near 112 and 118 deg with dihedrals near
+/-172 deg; with equal rigid angles the only straight all-trans chain is
exactly 180 deg).

**Change.** Every RIS sequence has a line-group symmetry that can be read off
the state sequence: TT and T are 2/1 screws or mirror zigzags, TG+TG- has a
glide, TTTG+TTTG- has a glide, TG+TG+ is a screw. Parametrise the repeat with
that symmetry imposed (e.g. alpha: two free torsions (t, g) with the pattern
(t, g, t, -g); beta: one deflection and two bond angles), so that the penalty
has far less work to do.

**Correction, established during implementation.** The claim originally made
here, that imposing the symmetry makes the repeat transform a pure translation
*by construction*, is wrong. A glide's repeat transform is the square of an
improper isometry, which is a rotation by twice the glide's angle, not the
identity. What the symmetry actually buys is a reduction in the number of
closure conditions, from three scalar equations to one, and that one still has
to be solved. The implementation therefore measures the rank of the closure
Jacobian and solves the remaining conditions by a damped batched Newton step,
which reaches a rotation residual of 1e-14 to 1e-11 degrees, against the ~2e-6
floor of the arccos-based measure used elsewhere.
Add the backbone bond angles of the repeat as variables (two for PVDF). The
variable count drops from 6 + P to 6 + 2-4, and the refined geometry can
finally reproduce the deflected zigzag.

**Expected effect.** Fewer variables and a well-conditioned objective: a
further 2-4x on refinement on top of section 2, and physically better chains.
Effort: small-medium (a symmetry table keyed by the canonical sequence; the
generic path stays as fallback). Risk: low; the constrained space is a subset
of the current one plus the angles.

## 4. Run independent candidates in parallel

Candidates are independent: pack and refine each in its own process. This
machine has 12 cores; the pipeline currently uses one for packing and
refinement. A `concurrent.futures.ProcessPoolExecutor` over candidates (with
the main-module guard that Windows' spawn model needs, and the RIS model
passed by JSON) gives 4-8x on the pipeline's wall time immediately. On a GPU
the same independence becomes batching across candidates. Effort: small.
Risk: none.

**Measured, and twice wrong.** It gave 1.3x, because candidate costs are
unequal (section 12). And "12 cores" is 6 physical cores with SMT: the kernels
are bound by memory bandwidth, not by instruction issue, so `cpu_count() - 1`
workers oversubscribe. Measured on the PVDF funnel, same code and same cold
cache: 5 candidate workers 223 s, 3 workers 179 s. The pipeline's default
worker count should be the *physical* core count at most, and fewer when the
inner stages are themselves parallel: one worker per physical core already
saturates the memory system, and every level of nesting after that is pure
contention (5 candidate workers x 2 table workers ran past 600 s before it was
abandoned, against 223 s unnested).

## 5. RIS fitting: make the expensive potential's job small

The fit is the *only* stage that calls the expensive potential, so its
algorithm decides the cost of production runs with an MLIP or DFT, where a
single point costs 0.1-100 s rather than microseconds. Today it is a dense
36 x 36 scan per bond type plus 4 S^3 triplet points: 2,881 evaluations for
PVDF, of which the vast majority are far from any basin minimum and never used.

**Change.** Coarse-to-fine: a 30 deg grid (12 x 12 per bond type) to locate
the S x S basins, then a local 2-D minimisation per basin (Nelder-Mead or
Powell, ~25 evaluations, or a 3 x 3 stencil with quadratic interpolation),
then the triplet corrections. Roughly 2 x (144 + 9 x 25) + 216 = 950
evaluations, 3x fewer, and nearly all of them informative; with the basins'
locations known from the illustrative potential, the coarse grid can be
skipped for the expensive potential and the count drops to ~650. Independent
of that, build all conformers of a scan in one batched NeRF call (the
per-conformer Python `build_chain` is the whole cost of the fit with the cheap
potential) and hand the batch to `energy_batch`, which for an MLIP means one
batched forward pass.

Expected effect: 3-4x fewer expensive evaluations; 10x on the fit with the
built-in potential. Effort: small. Risk: a basin missed by a coarse grid;
keep the dense scan as an option.

## 6. Kernel-level algorithmic items (for whatever direct evaluations remain)

* **Intra-chain terms are a constant for a rigid chain.** Compute once per
  chain and add; remove the (0,0,k) column and its scale tensor from every
  call. Saves the same-site work and the bonded-exclusion bookkeeping.
* **Per-pair z-window instead of a global K.** A pair at lateral distance
  d_xy can only interact with z-images satisfying |z_p - z_q - k c| <
  sqrt(rc^2 - d_xy^2); on average far fewer than 2K + 1. With the table of
  section 1 this becomes moot for the screen; it still helps the polish.
* **Tabulated pair potentials.** The pair *types* are fixed, so V_type(r) =
  LJ + shifted-force Coulomb can be a cubic spline in r^2 per element pair
  (3-6 tables): one gather and a polynomial instead of pow, exp and the erfc
  approximation. Standard in MD codes; 2-3x on the kernel's arithmetic.
* **Image selection is a lattice-geometry computation**, independent of
  setting angles and dz; compute it per (a, b, gamma) triple rather than per
  configuration, and vectorise the per-configuration Python loop.
* **float32 for the screen** (results within 1e-3 kcal/mol); float64 only in
  the final polish. Halves memory traffic, which is what bounds the kernel.
  Memory traffic being the bound is now measured rather than asserted (section
  10), and the *working-set size* turned out to matter more than the element
  width: the table build's inner chunk went from 131,072 elements to 16,384 and
  gained 1.3x on one core and 4x on six, because 1 MB of live temporaries per
  worker fits in a shared 12 MB L3 and 8 MB does not. `pack.CPU_CHUNK_ELEMS`
  (60,000) has never been re-examined under that light and is the cheapest
  experiment left in this section.
* **Search-space symmetry.** Use the chain's own screw symmetry (2/1: phi ->
  phi + 180 with dz -> dz + c/2), the a <-> b swap with a 90 deg rotation, and
  the periodicity of dz to cut the search domain by 4-8x for any sampling
  strategy.

## 7. Amorphous statistics: exact where exact is cheap

Diad, tetrad and run-length probabilities of a Markov chain are computable
exactly from the forward-backward pair marginals in O(N S^2); run-length
distributions come from a small auxiliary DP (probability that a run of
length >= k starts at bond i). This replaces 10^5-chain sampling for the
summary statistics, removes sampling noise, and reduces the sampled ensemble to
what only samples can give (3-D chain statistics, interface geometry). The
per-row Python loop in `run_lengths` should in any case be a single vectorised
diff over the padded boolean matrix. Speed matters little here (0.4 s), but
the exactness is worth having.

## 8. Fewer candidates reach the expensive stages

* **Cross-section pre-screen.** A chain's 2-D projection bounds the density
  any packing can reach; candidates that cannot reach 85% of the best
  candidate's bound do not need packing.
* **Cluster by helix.** Sequences with the same c, radius and cross-section
  moments pack the same way; pack one representative per cluster and only
  refine the members.
* **Enumeration cost is already negligible**; `helix_parameters` builds a chain
  three times per candidate and could build once, but it is 3 s of a 3-minute
  run.

## 9. GPU: what it should be on this machine

* **Backend.** CuPy is CUDA-only here; the practical GPU path on this Windows
  machine is torch, either through torch-directml (works today: float32 only,
  higher per-launch overhead than CUDA, most elementwise / reduction / einsum /
  gather ops available) or through AMD's ROCm-on-Windows PyTorch builds for
  RDNA3 if the driver stack supports them (lower overhead, float64). A torch
  implementation of the packing kernel is worth having regardless of device
  because it also delivers autodiff gradients for section 2. Keep `xp` as the
  abstraction, but a torch backend needs a small shim: `concatenate -> cat`,
  `sum(axis=) -> sum(dim=)`, `asarray -> as_tensor(device=)`, index arrays as
  tensors on the device, no in-place assignment into views that DirectML
  cannot alias.
* **What runs on the GPU.** The W-table build (section 1), the polish and
  refinement gradient batches (section 2), batched packing across candidates
  (section 4), and the conformer batches of the fit (section 5). On the
  RX 7700 XT (about 576 GB/s), the element-wise kernels are 50-100x faster
  than 12 CPU cores once batches are large enough to hide launch overhead,
  which on DirectML means batches of tens of thousands of cells or table
  entries, not one cell at a time.
* **What stays on the CPU.** The Viterbi / forward-backward recursions (tiny
  state vectors, sequential along the chain), the NeRF loop for single chains,
  and the FFT lattice sums of section 1 (microseconds each). The backward
  sampling loop is batched over chains but sequential in bonds: 200 steps of
  a few launches each is fine on CUDA and marginal on DirectML; keep it on the
  CPU unless ensembles exceed ~10^5 chains.
* **Projected end-to-end effect.** This bullet predicted that the CPU pipeline
  would end up dominated by the table builds and the final polishes, and that a
  GPU would then matter most for the table builds. It was then briefly wrong --
  with the build split across processes it was 10% of the funnel against
  refinement's 69% -- and is now right again for a different reason: the analytic
  gradients of section 13 took refinement's kernel rows away, so the table build
  is 34% of the funnel and its largest single stage (section 10). A GPU should
  take the table build; there is no gradient *batch* left for it to take, because
  a gradient is now one row. Still projected: nothing has been run on the GPU
  beyond a device test.

## 10. Where the time actually goes (measured), and what to do next

All numbers below were measured on the target machine -- Windows 11, Python
3.12, **Intel i7-8700K, 6 physical cores / 12 threads**, `OMP_NUM_THREADS=1`,
with another agent's work running alongside -- as the minimum of repeats, and
each "before" was taken back to back with its "after". The workload is the
default `PipelineConfig` for PVDF with `SimpleFF`: a fit, enumeration to period
8, and five conformations packed and refined (TG+TG+TG-TG-, TTTG+TG-,
TTG+TG-G-TG+, TG+TG-, TT).

### The funnel, stage by stage

One process, no candidate parallelism -- so the columns are CPU work, not a wall
time with contention in it. The two gradient columns were measured in the *same*
process, from the same screen starts, with the chain-pair tables built once and
shared, so nothing but the gradient differs between them. The "central
differences" column is a fresh measurement of the same code that previously read
34.7 s and 165.0 s for the polish and the refinement; it comes out 10% higher
because the machine was busier. Compare within a row pair, never across runs.

| Stage | `bb2e4de` | central differences | **analytic gradients** | share now |
|---|---|---|---|---|
| RIS fit (`fit_ris`, step 10 deg, 6 monomers, third order) | 0.31 s | 0.33 s | 0.33 s | 0.5% |
| enumeration (period <= 8, k = 60) | 4.9 s | 4.6 s | 4.6 s | 6% |
| periodic chain construction (5) | 0.12 s | 0.13 s | 0.13 s | 0.2% |
| **chain-pair table builds (5, screen grid)** | 66.1 s | 24.3 s | **24.3 s** | **34%** |
| FFT screen (5 conformations, ~10^8 landscape points each) | 8.1 s | 8.3 s | 8.3 s | 12% |
| exact-kernel polish (4 starts per conformation) | 33.8 s | 38.6 s | **8.0 s** (4.9x) | 11% |
| **torsion + cell refinement (5)** | 155.7 s | 181.2 s | **25.7 s** (7.1x) | **36%** |
| amorphous + interface statistics | 0.38 s | 0.48 s | 0.48 s | 0.7% |
| total | 269.4 s | 258.0 s | **71.8 s** | |

In kernel rows -- one row being one configuration through the pair kernel, the
unit the two earlier sections argue about -- the collapse is larger than the time:

| Stage | central differences | analytic | rows saved |
|---|---|---|---|
| polish (20 L-BFGS runs) | 3,584 rows | 343 rows | 10.4x |
| refinement (5 candidates) | 13,115 rows | 512 rows | 25.6x |

The time gains less than the rows because an analytic row is not a plain row: it
carries the pair derivative and the reductions onto atoms and lattice vectors, and
measures about 2.2x a plain row (1.7-2.3x over polyethylene, alpha and gamma and
repeated runs; `examples/benchmark.py` prints it). A gradient therefore costs ~2.2
rows against 11 for a cell-only difference and 29 for the gamma chain's full
variable set, which is the 4.9x and 13x that benchmark measures directly and the
4.9x and 7.1x the two stages show once the geometry builds are paid too.

Per candidate, refinement before and after (time / rows / L-BFGS evaluations):

| Candidate | atoms | vars | central differences | analytic | speedup |
|---|---|---|---|---|---|
| TG+TG+TG-TG- | 24 | 10 (linegroup) | 24.4 s / 1,198 / 57 | 5.8 s / 69 / 68 | 4.2x |
| TTTG+TG- | 18 | 13 (free) | 30.2 s / 4,510 / 167 | 5.6 s / 180 / 179 | 5.4x |
| **TTG+TG-G-TG+** | 24 | 15 (free) | **124.2 s / 6,697 / 216** | **13.2 s / 218 / 217** | **9.4x** |
| TG+TG- (alpha) | 12 | 8 (linegroup) | 2.3 s / 579 / 34 | 1.0 s / 34 / 33 | 2.2x |
| TT (beta) | 6 | 6 (linegroup) | 0.24 s / 131 / 10 | 0.12 s / 11 / 10 | 2.0x |

The lopsidedness is gone with it: the worst candidate was 47% of the whole funnel
and is now 18% of a funnel 3.6x smaller.

### End to end

Same script, same conditions, run back to back; `workers=3` is the current
default (one per physical core, minus the main process):

| Run | central differences | analytic | speedup |
|---|---|---|---|
| PVDF funnel, `workers=3`, cold table cache | 188.8 s | **77.2 s** | 2.45x |
| same, warm on-disk table cache | 162.5 s | **25.4 s** | 6.40x |

The answer is unchanged: beta-PVDF TT wins at -8.153 kcal/mol per monomer,
a = 4.60, b = 8.59, c = 2.63 A, in every run, and every candidate's refined
energy agrees to better than 1e-3 kcal/mol per monomer with two exceptions worth
naming, because neither is a different minimum:

* **TT (beta) is reported "para" on the analytic path and "anti" on the
  finite-difference one, at the same energy (1.4e-6 kcal/mol per monomer apart)
  and the same cell.** For an all-trans chain the flip `(x, -y, -z)` *is* a
  translation along z -- measured residual 3.6e-15 A for beta-PVDF and 2.9e-15 A
  for polyethylene, against 0.8-1.5 A for alpha and gamma -- so `flip=0` and
  `flip=1` describe the identical crystal and the label is decided by a tie-break
  at the 1e-11 level. The polarity of an all-trans cell is set by `phi` and `dz`,
  not by the flip flag, and any reading of that flag for beta (including this
  document's earlier "antipolar") says nothing.
* **TTG+TG-G-TG+ ends 8.4e-3 kcal/mol per monomer *lower* on the analytic path at
  the pipeline's default budget.** Neither route has converged there: both stop at
  `maxiter=200` (216 and 217 evaluations), and the analytic gradient simply gets
  further down the same valley per iteration. Lifting the cap to 2,000 settles it
  -- both land on the same minimum, -5.4381 against -5.4377 kcal/mol per monomer,
  a within 1e-4 A, b within 1e-4, c within 2e-4, every torsion within 0.024 deg --
  and the analytic route takes 18.4 s to the finite difference's 208.5 s (11.3x).

### The table build against worker count

gamma-PVDF, screen grid, warm pool, **bit-identical `W` at every worker count**:

| workers | 1 (serial) | 4 threads | 2 proc | 3 proc | 4 proc | 6 proc | 8 proc | 12 proc |
|---|---|---|---|---|---|---|---|---|
| time | 35.2 s | 61.3 s | 15.1 s | 10.7 s | 7.6 s | 7.9 s | 7.4 s | 7.7 s |
| speedup | 1.00x | 0.57x | 2.34x | 3.31x | 4.62x | 4.48x | 4.75x | 4.60x |

Not linear, and it saturates at four workers on six cores: the build is bound by
memory bandwidth, not by cores. The default is one worker per *physical* core
(`lattice_table.default_table_procs`), 1 inside an outer process pool, and
`$POLYFIND_TABLE_PROCS` overrides both; `$POLYFIND_TABLE_CACHE` still names the
directory the finished tables are kept in, and has to be set for the warm-cache
row above to happen at all. Two things follow from the shape of the curve, and
both were surprises.

First, **the chunk size was the real bottleneck, not the GIL.** The inner kernel
held a 131,072-element working set, about 8 MB of live temporaries, so even one
worker streamed from DRAM. A probe of pure headroom -- N processes each running
the *whole* build concurrently -- showed six processes delivering only 1.8x one
process's throughput at that chunk, and 3.9x once the chunk was cut to 16,384
elements (~1 MB, L3-resident for six workers). Cutting it also made the serial
build faster (gamma 30.4 s -> 28.2 s). Threads were never going to get past
1.3x, and processes would not have either.

Second, **with the smaller chunk, threads are worse than useless**: each NumPy
call is eight times shorter, so the loop is pure GIL contention (gamma 35.2 s on
one thread against 61.3 s on four). The in-process fallback is now single-threaded.

Cold builds, screen grid, before (4 threads) vs after (6 processes): PE all-trans
5.14 -> 1.86 s, beta-PVDF 4.98 -> 1.59 s, alpha-PVDF 13.81 -> 3.65 s, gamma-PVDF
23.70 -> 8.22 s. On the finer default grid: PE 40.1 -> 18.3 s, gamma 157.4 ->
87.1 s (the 2x rather than 4x there is a first-call pool spawn of ~2 s plus a
72 x 72 angle grid that leaves only three atom pairs per kernel call).

### Next, in order

1. **Chain-pair table builds (34%).** Now the largest single stage, and the one
   whose scaling is already understood (section 13): bandwidth-bound, saturating
   at about four worker processes on six physical cores, and already 1 worker
   inside an outer candidate pool. Two levers are left -- the float32/spline
   arithmetic of section 6, and the GPU of section 9, which is where the table
   build becomes the natural job again now that the gradient batches have gone.
2. **Refinement (36%).** The gradient is no longer the cost; the cost is now the
   *geometry* -- each iteration still builds `1 + 2 n_shape` chains (one NeRF
   pass, plus a closure solve in the line-group parametrisation) to contract
   against the exact `dE/d(coords)`. Making the torsion and bond-angle Jacobian
   analytic as well (`dr_i/dtheta_j = u_j x (r_i - p_j)` through the NeRF build,
   the closure solve by the implicit function theorem, and the Kabsch/alignment
   stage by differentiating the SVD) would remove it; so would simply using the
   *line-group* parametrisation everywhere, since it has the fewest shape
   variables. Nothing here is a kernel problem any more.
3. **Sequential L-BFGS runs.** The polish is four independent L-BFGS runs per
   conformation and is parallel across starts; at 11% of the funnel that is now a
   small prize.
4. **GPU (section 9).** Back to the table build, not the gradients.

### Status of the earlier plan

| Step | Change | state |
|---|---|---|
| 1 | section 4: parallel candidates | done; 1.3x, not 4-8x, and oversubscribed at the default worker count |
| 2 | section 2: batched-gradient polish and refinement, hoisted invariants | done |
| 3 | section 6: intra constant, z-windows, float32 screen | done; tabulated pair potentials not done |
| 4 | section 1: tabulated W + FFT lattice sums, exhaustive screen | done; the screen is 3% of the funnel |
| 5 | section 3: symmetry-parametrised refinement with bond angles | done; 69% of the funnel then, 36% after step 9 |
| 6 | section 5: coarse-to-fine fit, batched conformers | done; the fit is 0.1% of the funnel |
| 7 | section 9: torch backend on the AMD GPU | not done |
| 8 | process-parallel table build (this section) | done; 2.8-3.8x on the build, 1.30x end to end |
| 9 | section 13: analytic gradients for the polish and the refinement | done; 4.9x on the polish, 7.1x on refinement, 2.4-6.4x end to end |

## 11. Open points

* The Fourier lattice sum assumes W is band-limited enough in the angles and
  dz for a 72 x 72 x 16 grid; the repulsive wall at contact may need 144
  angular points for chains with large substituents. Cheap to test on PE and
  beta-PVDF against the direct kernel.
* Line-group parametrisation of the refinement needs a table of symmetry
  patterns per canonical sequence; sequences without a recognised symmetry
  fall back to the penalty method.
* The torch-directml route has not been benchmarked here; the RX 7700 XT's
  ROCm-on-Windows support should be checked before committing to DirectML.

## 12. Findings from implementation

Implementing section 5's coarse-to-fine fit (`fit_ris(scan="adaptive")`) surfaced a
modelling assumption the RIS discretisation makes silently: that each state pair's
energy, taken as a basin minimum, sits at an interior critical point of a well. For
PVDF's CF2-centred (G+,G-) pair under `SimpleFF` this is false; the energy falls off
monotonically as the dihedral approaches the neighbouring trans basin, with no turning
point. The reported pair energy is then a property of the basin convention rather than
of the potential alone: 8.17 kcal/mol on a step=10 deg grid, versus 4.85 kcal/mol as the
true constrained infimum at the basin edge. A principled fix, if this matters
scientifically, is to define each RIS pair energy as -kT ln of the Boltzmann integral
over its basin rather than as a basin minimum; that is the classical RIS definition, and
it is insensitive to where (or whether) the minimum sits.

Implementing section 1's tabulated chain-pair interaction surfaced a correctness bug in
the packing kernel, described in DESIGN.md section 5.6: every antiparallel configuration
carried a spurious constant of thousands of kcal/mol, so the antipolar half of the search
space had never been reachable and every packed result in the documentation was parallel
by construction rather than by energetics. It was found precisely because the tabulated
formulation forces the intra-chain term to be isolated and named, where the direct kernel
folded it into a per-configuration sum and hid it. That is worth recording as an argument
for the reformulation independent of its speed: making a decomposition explicit exposes
terms that an aggregate never checks. The fix and its regression tests are in `pack.py`.

Section 4's process parallelism measured 1.3x, not the 4-8x estimated here. The estimate
assumed candidates of comparable cost; in practice one candidate dominates (the gamma
chain has 24 atoms per repeat against beta's 6), so wall time is set by the slowest
candidate and adding workers cannot help. Speeding that stage further has to come from
making the expensive candidate cheaper, which is what sections 1 and 2 do, or from
parallelising within a candidate across its polish starts.

Section 5's estimate was also wrong in an instructive direction. The projected win was a
3-4x reduction in expensive-potential calls from the coarse-to-fine scan, and that held
(2.84x). But it was not the bottleneck: building each scan conformer in Python cost far
more than evaluating it, and batching the builds alone made the fit 30x faster with
bit-identical energies. The adaptive scan therefore ships as an opt-in for genuinely
expensive calculators, where its evaluation count is what matters, rather than as the
default.

Section 3's implementation also settled a question the section raised, in the
negative. Adding bond angles as variables does **not** let the beta chain
reproduce its experimental geometry, and the reason is a constraint rather than
a limitation of the optimiser. With one monomer per repeat, an all-trans
chain's rotation residual is exactly the difference between its two backbone
angles, and no choice of torsions can cancel it, so exact periodicity forces
equal angles and exact trans. Allowing a two-monomer repeat does admit unequal
angles through the glide, but the required torsion deflection is steep: the
experimental 6-degree angle difference needs dihedrals near 153 and 207 degrees
rather than the observed +/-172, at a cost of about 2.2 kcal/mol per monomer.
With rigid, equal-length bonds the experimental combination of angles and
torsions is simply not realisable. Relaxing bond lengths, not just angles, is
what that would take.

Wiring the tabulated screen into `pack()` produced the strongest scientific
result of the exercise, and it depended on the flip fix above. With the
antiparallel branch reachable for the first time, alpha-PVDF now packs
antiparallel and wins outright, 0.62 against 0.80 kcal/mol per monomer for the
best parallel cell, at the experimental antipolar cell. That is the correct
answer for alpha-PVDF and the code could not previously express it. Both
screens find it, so it is a consequence of the correctness fix rather than of
the exhaustive search.

The exhaustive screen matched the random screen's global minimum on all four
chains while using five to nine times fewer exact-kernel evaluations, and it is
three to six times faster once a table is cached. Its advantage is not a lower
minimum but determinism: it cannot be unlucky. The two screens differ in which
*secondary* minima they surface, and interestingly the random screen sometimes
finds shallow distant basins the exhaustive one does not, because the latter's
best cells concentrate near the deepest region.

A second latent bug surfaced during that work, again through a reformulation
rather than a test: the table's angular interpolation took its weight from the
wrapped index, so an angle wrapping to exactly 360 degrees received weights of
-47 and +48. Oblique cells scored -670 kcal/mol where the exact kernel gives
-6.5. Right-angled cells were unaffected, so no previously reported packing
result changed, but every non-orthogonal cell would have been nonsense.

A note on how that last finding was nearly lost. The subagent that found the
chirality problem also proposed the correct replacement symmetry and gave its
index form. Checking it, I mapped the bond-type index on the pair term when the
relation does not, got a large residual on a polymer where it must vanish, and
recorded the relation as unresolved. The agent was right and the check was
wrong. Two lessons: a verification that fails on the case where the answer is
known is evidence about the verification, not only about the claim; and the
fix was to make the code validate the relation itself at fit time, so it is
applied only where it measurably holds, rather than to trust either of us.

The chirality work corrected my brief a second time, in the same direction. I
told that agent that "chain reversal stays a symmetry in both cases". It does
not. Reading an isotactic chain backwards swaps each stereocentre's neighbours,
so plain reversal returns the enantiomer just as reflection does; measured on
24-bond oligomers, reversal alone and reflection alone give bit-identical
energies, both differing from the original by up to 749 kcal/mol for CFE and
5103 for CDFE, while the two composed give zero to machine precision. All three
are zero for PE, PVDF and PVDC, which is why the achiral cases never exposed
it. The symmetry group of a chiral repeat is therefore shifts plus
reflection-with-reversal, and nothing else - exactly the operation the fitting
code had already been driven to.

Implementing the analytic gradients corrected a claim this document had been
repeating. Section 10 reported the winning beta-PVDF cell as "antipolar", read off
the `flip` flag of the packed result. For an all-trans chain that flag carries no
information: the flip `(x, -y, -z)` is *exactly* a translation along z (measured
residual 3.6e-15 A for beta-PVDF, 2.9e-15 A for polyethylene, against 0.8 A for
alpha and 1.5 A for gamma), so `flip=0` and `flip=1` describe the identical
crystal and which one a run reports is a tie-break between energies 1e-11 apart.
It surfaced because the analytic polish perturbs those last bits and the reported
label changed while the energy, a, b and c did not. The polarity of an all-trans
cell lives in `phi1`, `phi2` and `dz`; anywhere a conclusion rests on the flip
flag, it has to be checked against the chain's own symmetry first. The same
arithmetic says the flag *is* meaningful for alpha and gamma, which is where the
antipolar result of this document's earlier work was obtained, so that result
stands.

The same agent's own achiral control caught a false alarm before it reached me.
Its first angle check reported a 120-degree symmetry violation for PVDC, an
achiral polymer where the residual must vanish. Rather than report it, it looked
again and found the cause: that polymer's trans basin is a degenerate double
well whose two minima differ by 4e-13 kcal/mol, so the choice between them is an
argmin tie-break, not a physical asymmetry. The check was wrong, not the claim.
That is the second time on this project that a verification failing on a known
case was evidence about the verification, and the first time the discipline was
applied without me having to intervene.

Splitting the table build across processes produced the most instructive
correction of the exercise, and it corrected a *diagnosis* rather than an
estimate. The build's own note said it was "threaded but GIL-bound, so eleven of
twelve cores sit idle", and the obvious reading -- replace the threads with
processes and the cores fill up -- is what the work set out to do. It did not
work: the first correct process build was no faster than the threaded one, and at
twelve workers it was slower. The GIL was the second bottleneck, not the first.
The inner kernel carried a 131,072-element working set, about 8 MB of live
temporaries, so a single worker was already streaming from DRAM and there was
nothing for more workers to use. One probe settled it: N processes each running
the *whole* build concurrently, which at that chunk size gave six processes only
1.8x the throughput of one. Cutting the chunk to 16,384 elements -- about 1 MB,
so six workers' temporaries fit in the 12 MB L3 -- raised that to 3.9x and made
the serial build faster as well, and only then did splitting the radial axis give
4.6x. The same change then inverted the threading result: with each NumPy call
eight times shorter the GIL really does dominate, and four threads became 1.7x
slower than one, so the in-process fallback is now single-threaded. "Parallelise
it" was the right instruction and the wrong first step, and the experiment that
would have said so in a minute -- run the existing code N times at once and look
at the throughput, not at one worker's time -- is the measurement to take before
parallelising any numerical kernel, because it separates "no parallelism" from
"no headroom".

Two smaller findings came with it. Bit-identity was free rather than hard-won:
radii are independent in the build, including under both exchange symmetries and
the cap, so any partition gives the same bits and the test asserts exact equality
rather than a tolerance -- worth recording because the instinct with float32
accumulation is to reach for `allclose`, which would have hidden a real coupling
had one existed. And "12 cores" on this machine is 6 physical cores with SMT,
which matters for a bandwidth-bound kernel: the build saturates at four workers,
and the pipeline's default of `cpu_count() - 1` candidate workers is a 1.24x
*loss* against three (section 4). Nesting the two levels of parallelism is worse
again. The on-disk table cache, by contrast, needed nothing: its writes were
already atomic, its keys are per conformation so a multi-candidate run has no
duplicates to remove, and a warm cache takes the PVDF funnel from 179 s to 146 s.

## 13. The bandwidth lesson, a lead that did not transfer, and the rows that did

Parallelising the interaction-table build produced the most transferable finding
of this whole effort, and it was not the parallelism. Splitting the build across
processes bought *nothing* at first. The inner kernel held an 8 MB working set,
so a single worker was already streaming from main memory and additional workers
only contended for bandwidth it did not have. Shrinking the chunk to about 1 MB,
small enough to stay in cache, raised throughput by 3.9x and made the *serial*
build faster as well. Only after that did splitting pay, and the same change
inverted an earlier result: threads went from helping to being 1.7x worse than
serial. Scaling then saturated near four workers on six physical cores, which is
a bandwidth ceiling rather than a core count.

That single constant was worth more than the parallelism it enabled, which is
worth remembering before reaching for more workers anywhere in this package.

**The obvious follow-up does not transfer, and I measured it rather than
assuming.** The packing kernel has the same kind of constant, and since
refinement and polishing together are 84% of what remains, a similar win there
would have been large. Sweeping it from 4,000 to 500,000 elements across the
polyethylene, alpha and gamma chains, for both a gradient-sized batch of 29 rows
and a screen-sized batch of 512, moves the cost by at most 5%, which is inside
the noise; the present 60,000 is already near optimal, and only the extreme
500,000 is clearly bad. The difference from the table build is that the table
kernel was being run in parallel, where bandwidth contention dominates, while
this measurement is single-process.

**Where the time was, and the lever that moved it.** Before this section's work
the profile was refinement 69%, polishing 15%, table build 10%, FFT screen 3%,
everything else 2%, and 112 of refinement's 165 seconds sat in a single candidate
whose per-row kernel cost is 18 ms against 2.3 for alpha and 0.9 for
polyethylene - an N-squared effect from its 24-atom repeat that no constant will
fix.

The lever was not the kernel but the *number of rows*. Refinement got its
gradients by central differences, so every L-BFGS iteration cost one plus twice
the variable count - 29 rows for the gamma chain - to extract a single gradient.
**That is now done analytically, and the prediction in the paragraph this
replaces was about right in both directions:** 25.6x fewer rows in refinement
(13,115 to 512) and 10.4x fewer in the polish, at 1.7-2.2x the cost per row, for
7.1x on refinement and 4.9x on the polish. The funnel went from 258.0 s to
71.8 s of CPU work and the end-to-end run from 188.8 s to 77.2 s cold, 162.5 s to
25.4 s warm. Section 10 carries the per-candidate numbers.

**What is analytic and what is not.** The lattice energy is a sum of pair terms,
so `dE/dr_i` is a scalar times the separation vector, and
`CrystalPacker.energy_and_grad` returns, from one kernel row, the energy together
with `dE/d(a, b, gamma, phi1, phi2, dz)`, `dE/d(chain coordinates)` and `dE/dc`.
The cell variables are closed-form throughout: the lengths and the cell angle
through the lattice vectors and chain 2's offset, `dz` through that offset alone,
the setting angles through `dr/dphi = z_hat x r`. The *conformational* variables
are not: `d(coords, c)/d(shape)` runs through the NeRF build, the line group's
closure solve and the Kabsch alignment, and that Jacobian is still a central
difference - but of the **geometry only**, contracted against the exact
`dE/d(coords)`, so it costs no kernel row at all. That is why the row count
collapses by 25x while the gradient of a torsion is still numerical.

**Verification was the deliverable, not the speed.** A gradient wrong in one
component is worse than no gradient, because L-BFGS converges confidently to the
wrong place. Every component was checked against a central difference of the very
function the optimiser minimises - not a re-implementation of it; `refine_crystal`
grew a `probe` argument precisely so that the objective and both gradient routes
can be driven from a test. Over polyethylene, beta-, alpha- and gamma-PVDF and an
unpatterned sequence, with and without a field, in both parametrisations, at two
points each: worst relative error **2.3e-7** on the analytic cell variables and
**2.3e-5** on the geometry-differenced shape variables. The energy
`energy_and_grad` returns is bit-identical to `energy()`, asserted with `==`.

Two things that verification taught, both worth keeping:

* **Differentiate the approximation, not the function it approximates.** The
  Coulomb term uses the usual polynomial `erfc`, whose *value* is good to 1.5e-7 -
  but whose *derivative* differs from the true `-2/sqrt(pi) exp(-x^2)` by 1.2e-5
  relative over the range the kernel uses. Substituting the true derivative and
  re-measuring: the error on the `a` gradient goes from 4e-7 to 5.7e-6 and on
  `phi1` from 9e-10 to 4e-8, i.e. an order of magnitude of accuracy thrown away,
  and enough to sit above the finite-difference floor and look like a bug in the
  chain rule. The gradient of a kernel is the gradient of the kernel, not of the
  mathematics it stands in for.
* **The kernel has a real kink, and the finite difference is the thing that is
  wrong at it.** The Lennard-Jones term is energy-shifted at the cutoff but not
  force-shifted, so its force jumps at `r = rc`. A step that straddles that point
  makes the *difference* disagree with the derivative by ~1e-4 relative, and the
  disagreement vanishes at smaller steps rather than shrinking as `h^2`. Three
  components looked wrong until the step was swept; moving the cutoff by 0.04 A
  moved the anomaly with it. Force-shifting the Lennard-Jones term, or splining
  it to zero over the last half-angstrom, would remove a genuine non-smoothness
  from the objective that the old finite-difference refinement was walking over
  blind.

A smaller fix was applied alongside: the pipeline's worker count defaulted to
logical cores minus one, which oversubscribes a bandwidth-bound workload. On
this machine three workers beat five by 1.24x, so the default is now derived
from physical cores.

## 14. GPU status on this machine, measured, and the decision

Probed rather than assumed, so nobody repeats it:

| Check | Result |
|---|---|
| ROCm / HIP tooling on Windows | absent: no `rocminfo`, `hipcc`, `rocm-smi`, no ROCm install directory |
| PyTorch build present | 2.4.1+cpu; `torch.version.hip` is None, `cuda.is_available()` False |
| WSL2 Ubuntu 24.04 | present, exposes `/dev/dxg` but **not** `/dev/kfd` |
| torch-directml | works; sees the RX 7700 XT and executes kernels |

`/dev/kfd` is the ROCm compute device, so ROCm cannot run in this WSL instance
whatever is installed; what is exposed is DirectX paravirtualisation, which is
what DirectML uses. The card is gfx1101, which generally needs a version override
to pass as the officially supported gfx1100.

**Decision: no GPU work here.** A separate machine has HIP and ROCm working and
is committed to simulation work; large runs will go to rented GPUs when they are
needed. The local payoff also shrank sharply during this effort. A GPU looked
attractive while packing and refinement took minutes; after the analytic
gradients and the parallel table build the whole pipeline is 25 s warm and 77 s
cold, split between the table build at 34% and refinement at 36%. Refinement is
now latency-bound with few kernel rows, which is the shape a GPU serves worst, so
even a tenfold win on the table build saves about eight seconds.

The payoff returns with scale: many chemistries crossed with strain and field
grids multiplies the throughput-bound work and the table build dominates again.
At that point the torch route is the lower-risk one, since it already runs on
this card through DirectML, and it doubles as the automatic-differentiation path
section 2 wanted. CuPy remains unusable on an AMD card and should not be the
plan.

## 15. Second review (2026-09-29): what is left, measured on the current code

Scope: execution speed of the code at branch head `bd5ef12`, after sections 1-14
landed.  Every "before" below was measured on that head; every "after" is a
prototype made in an isolated worktree and run back to back with its baseline on the
same machine, a 4-core cloud box (Intel Xeon, NumPy 2.4, no GPU).  The machine is
slower and narrower than the i7-8700K of section 10, so compare shares and ratios,
not seconds, across the two sections.  The profile script is
`examples/profile_funnel.py`; the prototypes are kept as patches under
`docs/perf-patches/2026-09-29/` with an index of what each does and costs.

### 15.1 Where the time goes now

Default `PipelineConfig` for PVDF with `SimpleFF`, single process, no table cache:

| stage | seconds | share |
|---|---|---|
| RIS fit (2,881 single points) | 0.58 | 0.9% |
| enumeration (120 candidates, period <= 8, k = 60) | 3.95 | 5.9% |
| periodic chain construction (5) | 0.11 | 0.2% |
| **chain-pair table builds (5, screen grid, 4 worker processes)** | **25.55** | **38.5%** |
| **FFT screen (5)** | **12.04** | **18.1%** |
| exact-kernel polish (4 starts per conformation) | 2.57 | 3.9% |
| **torsion + cell refinement (5)** | **21.60** | **32.5%** |
| total | 66.40 | |

The same three stages as section 10, in a different order: this box has four cores
for the table build, and the thread pool inside `fft_screen` ran under the GIL.
The refinement row was taken under cProfile, which inflates Python-call-heavy code;
the same five refinements without the profiler, on a quiet machine, take 15.2 s
(finding G1 below), so the true share of refinement is nearer 25% and the true total
nearer 60 s.

Per candidate (seconds: table / screen / polish / refine; refinement evaluations,
variables, parametrisation):

| candidate | atoms | table | screen | polish | refine | evals / vars / param. |
|---|---|---|---|---|---|---|
| TG+TG+TG-TG- | 24 | 7.74 | 3.08 | 0.60 | 6.06 | 70 / 10 / linegroup |
| TTTG+TG- | 18 | 5.22 | 2.39 | 0.78 | 5.21 | 161 / 13 / free |
| TTG+TG-G-TG+ | 24 | 7.34 | 2.58 | 0.85 | 8.85 | 216 / 15 / free |
| TG+TG- (alpha) | 12 | 3.54 | 2.09 | 0.17 | 1.30 | 29 / 8 / linegroup |
| TT (beta) | 6 | 1.71 | 1.90 | 0.16 | 0.18 | 10 / 6 / linegroup |

**Refinement was mostly geometry, not kernel.**  cProfile of the five refinements,
cumulative seconds summed over candidates: `build_backbone` (the per-atom Python
NeRF loop, about fifteen small NumPy calls per atom) 7.1; pendant placement inside
`_block_coords` 4.1; the per-row Python loop of `repeat_chains` (Kabsch, `arccos`,
`_orient_block` with an `eigh`, `_align_about_z`, one `PeriodicChain` per row) about
3.6; the closure Newton `LineGroup._close`, which re-differences its Jacobian at every
iteration after the first, 4.2 (itself mostly `build_backbone`); the lattice energy
`CrystalPacker.energy_and_grad` 4.7.  Section 10 said "the cost is now the geometry";
this is that statement with numbers, and it says which geometry.

### 15.2 The earlier verification round, reconciled with the branch

A first round of this review (2026-09-09, against the code before sections 1-14)
verified three findings by prototype before its budget ran out: a compacted packing
kernel (12-17x on batched cells), a fast NeRF builder (10-55x on `build_backbone`) and
a candidate process pool (2.2-2.4x).  The branch has since taken the process pool
(section 4), replaced the random screen the kernel served with the table screen
(section 1) and made the polish a one-row gradient (section 13), so the kernel finding
is now worth about 3 s of 66 and was not pursued again.  The fast builder was never
taken; 15.1 shows it sitting under the largest remaining stage, and it is G1 below.

### 15.3 Findings, measured

Every prototype ran the full suite with `-k "not byte_identical"`; the ten failures
that remain are the same ten on the untouched head in this container (see 15.6) and
none touches the code a finding changes.

**G1. Fast chain geometry for the refinement** (`chain.py`, `pack.py`,
`linegroup.py`; `G1.diff`, 702 lines).  The NeRF recursion runs on component arrays
of shape `(N+3, M)` with an explicit cross product, the per-atom-type scalars hoisted
and the dihedral trig done once for all steps; a Python-float path serves `M == 1`;
pendant frames are built for every backbone atom in one call and each pendant atom is
placed with one gather and one scatter; the per-row tail of `repeat_chains` (Kabsch,
rotation to z, orientation, alignment) is batched.  Measured, quiet machine, back to
back: refinement stage 15.15 -> 8.93 s (1.70x); per candidate 1.5-1.9x; `build_backbone`
65x at `M = 1` and 3.3x at `M = 31`; `repeat_chains` on 31 rows 6.5-9.9x.  Bit-identical
outputs everywhere (155 dumped geometry arrays, seven polymers including the multi-atom
pendants of AN/VDCN/FANOME, every refined energy, and the full 216-evaluation L-BFGS-B
trace of the free-mode candidate).  That bit-identity cost something worth recording:
an intermediate version equal only to 1e-13 A moved one free-mode refinement onto a
neighbouring minimum, 0.01 kcal/mol per monomer away, through a one-ulp gradient
difference, so the final version reproduces NumPy's operation order exactly (the 1-D
norm as a BLAS dot, `s ** 2` through a Python float).  That order is a property of
this NumPy/BLAS build; on another the batched and per-row forms may differ by an ulp
again.  The refinement objective is chaotic at that level, which is a property of the
original, not of the change.  After G1 the refinement is kernel-bound: the analytic
row is 63-87% of each candidate, and the line-group closure is second (0.8 of 1.9 s
for the 8-bond glide chain).  A one-line fix found on the way, zeroing the closure
step for rows already inside tolerance, removes a batch-size dependence of the closure
that currently rests on rounding; it was reverted to keep outputs identical and should
be adopted on its own.

**G2. Fewer geometry builds per refinement iteration** (new `chainjac.py`,
`linegroup.py`, `pack.py`, `refine.py`, `tests/test_chainjac.py`; `G2.diff`, 733
lines).  The derivative of every oligomer coordinate with respect to every torsion is
closed-form (`dr_i/dtheta_j = u_j x (r_i - p_j)` for downstream atoms, carried through
the virtual end atoms and the rigid pendant frames), then through the repeat transform
(Kabsch derivative), the alignment and `c`; the closure solve reuses its Jacobian and
freezes converged rows.  Chain builds per L-BFGS evaluation: 98 -> 3 (TG+TG+TG-TG-),
17 -> 1 and 21 -> 1 (free mode), 59 -> 3 (alpha).  Gradient verified against a
Richardson-extrapolated central difference of the very objective the optimiser
minimises: 8.8e-7 of scale over 44 evaluations (PE, beta, alpha, gamma, the unpatterned
chain; both parametrisations; field on and off).  Measured alone: refinement stage
14.8 -> 11.9 s (1.24x).  The build count fell 30x and the time 1.24x because, on the
unchanged builder, one build costs the same for 1 row as for 31; G2 is worth much more
on top of G1, where the builds are cheap and the closure's remaining re-differencing
is the next cost, but the two were made independently and their combination was not
measured.  Linegroup results agree with the original to 2e-9 kcal/mol per monomer;
free-parametrisation results end in the same basin at a different flat-valley point,
1e-5 to 3e-5 kcal/mol per monomer in the objective (up to 6e-4 in the lattice energy
alone, traded against the bend restraint), so anything pinning refined numbers to many
digits will move.

**G3. Cheaper table-build kernel** (`lattice_table.py`; `G3.diff`, 475 lines).  The
premise was partly wrong and the measurement says how: the trig was about 1% of the
kernel, and a splined pair potential does not pay in NumPy (a gather costs about five
arithmetic passes, and the Horner on gathered rows is strided), so it ships as an
opt-in that is off.  What did pay: a buffered kernel with in-place passes and a single
mask, grouping pairs by element type, forming `r^2` from `cos(u - v)` (1.9-2.0x on the
exact path), and, enabled by the smaller working set, a chunk of 65,536 elements
instead of 16,384 (a further 1.2-1.3x; re-sweep on the i7-8700K with six workers, its
12 MB L3 may prefer 32,768).  Measured: serial build PE 5.20 -> 2.60 s, gamma 24.7 ->
11.4 s; with four processes PE 1.94 -> 0.67 s, gamma 7.51 -> 3.06 s; table stage of the
funnel 24.05 -> 10.55 s (2.28x).  Not bit-identical to the old tables: max 8e-4
kcal/mol at wall-adjacent entries, with the same accuracy against the float64 kernel
(max 1.5e-4 to 2.2e-4 at uncapped grid points, medians unchanged) and identical
polished energies on all twenty starts of the five candidates; serial and parallel
builds remain bit-identical.  On-disk caches built before the change stay valid but
will not match bit for bit.  Inside the new kernel the erfc polynomial with its `exp`
is about a third and the whole thing runs at about 7 ns per grid point against a
~5.5 ns floor for this many passes, so NumPy-level work is nearly exhausted there;
what is left is fewer grid points (20-40% of each chunk lies beyond the cutoff and is
still evaluated) or a compiled or fused kernel.

**G4. FFT screen without the thread pool** (`lattice_table.py`; `G4.diff`, 323
lines).  The premise was off by an order of magnitude: with `ab_symmetry` the screen
is 756 landscapes of 18,432 points per candidate, 14 million, not 10^8, and the 12 s
was GIL contention plus per-cell Python (two `lattice_sites` calls, a 16-term
interpolation, `np.unique`, a Python loop over sites with an outer product each,
`argmin` and `unravel_index`), not FFT work.  The prototype precomputes the sites and
shift multipliers for the whole `(a, b)` grid, interpolates all self terms in one call,
accumulates the Fourier sums for a batch of cells with array operations, and drops
the threads.  Measured: `fft_screen` 4.2-6.9x per candidate; screen stage 12.76 ->
2.26 s (5.6x).  Bit-identical `ScreenResult.energy`, angles, `dz`, top rows and
`_table_starts` output on all five candidates and nine edge cases (below `r_min`,
no sites, oblique gammas, single flips, 1x1 grids, all cells rejected).  What is left
in the 0.3-0.5 s per candidate: `irfftn` 0.12 s (scipy's with workers would halve it
but differs at the float32 ulp), the complex64 accumulation 0.12-0.15 s, and the cold
`table.fourier` of a fresh table, now 25-45% of the stage.

**G5. Enumeration** (`enumerate.py`, `helix.py`, `ris.py`; `G5.diff`, 448 lines).
One chain build per candidate and a batched helix analysis of all candidates of a
period (327 builds, 171 of them redundant, become four batched calls); `k_best` works
on one `(S, S*k)` candidate array per step with a vectorised backtrace, and
`cyclic_k_best` runs all start states in one pass.  Measured: `enumerate_periodic`
3.95 -> 0.06 s (60x); exact to 0.0 in energies, sequences, order, ranks and every
helix descriptor over 173 sequences, including chiral CFE.  The stage was 5.9% of the
funnel, which caps what this can ever deliver.

**G6. Table build on the array backend, for CUDA rentals** (`lattice_table.py`,
`tests/test_gpu.py`; `G6.diff`, 733 lines).  `_build_block` and `_pot_batch` used
`np` explicitly, so the table build could not run on a GPU even though `backend.py`
selects CuPy for the packing kernel.  The port takes the array module as an argument,
keeps the process-pool CPU path bit-identical (verified: SHA-256 of the tables equal
for PE, alpha and gamma, serial and four processes) and adds skipped-here GPU checks.
Nothing on a GPU was executed.  The agent's projection and an adversarial re-check of
it are both recorded, and the re-check wins: per gamma table, 6.9 s on four CPU
processes -> about 2-2.5 s unfused on an A10-class card (2.5-3.5x), 1-1.5 s if
`cupy.fuse` traces the kernel (5-7x, host-bound, and the first build pays a 1-3 s
NVRTC compile the agent had omitted); the term count in the brief was 7x too high
(1.3e9 grid elements, not 10^10), which shrinks the absolute payoff to a few seconds
per table.  The re-check also found four defects that must be fixed before the port
is used: `device="cuda"` with CuPy importable but no usable device returns the CUDA
module instead of the documented `ImportError` and then dies in `xp.zeros`; the
pipeline's `ProcessPoolExecutor` uses `fork` on Linux after the parent has touched the
CUDA runtime, so every worker's build raises and `_pack_and_refine_one` turns that
into "skip" for every candidate (needs a `spawn` context or `device="cpu"` inside
workers); `table_key` does not record the device, so a disk cache mixes CPU and GPU
tables that differ at 1e-4 kcal/mol; a new `cupy.fuse` object is created per block
call.  The torch-directml route for the RX 7700 XT remains a sketch (float32 only, no
fusion, op availability unverified).

### 15.4 GPU, again

The facts have not moved since section 14, and this round adds one measurement to
them: the user's card is AMD, CuPy cannot run on it, and every prototype above is a
CPU gain that the AMD machine collects in full.  On a rented CUDA machine the table
build is the one stage the port makes GPU-resident, and its projected effect on the
funnel is bounded by the stage: on the section-10 numbers (CPU work 71.8 s, tables
24.3 s) about 1.3-1.4x for a serial pipeline, nothing at all for the warm-cache run,
and less than G3 delivers on the CPU without a GPU.  The order therefore stands: GPU
work waits for the workload that multiplies table builds (many chemistries crossed
with strain and field grids, or the phonon supercells), and when it comes the four
defects above are the first day's work.

### 15.5 Order of implementation and expected effect

1. **G4** (bit-identical, one file): the screen stage 12.8 -> 2.3 s here.
2. **G5** (exact, three files, isolated): enumeration 4.0 -> 0.06 s.
3. **G1** (bit-identical on this build; re-check the ulp-level identity on the
   i7-8700K's NumPy before trusting refined digits): refinement 15.2 -> 8.9 s.
4. **G3** (changes table values below 1e-3 kcal/mol; clear `$POLYFIND_TABLE_CACHE`
   and re-sweep the chunk size on the target machine): table stage 24.1 -> 10.6 s.
5. **G2** on top of G1, with the closure-freeze fix: expected to take the refinement
   below 8 s here, unmeasured in combination; adopt once its 1e-5 movement of
   free-mode results is acceptable, and re-record any pinned refined values.
6. **G6** only when a CUDA machine is in play, after its four defects are fixed.

Projected from the measured stage ratios, single process on this box, with the
refinement row taken without the profiler:

| stage | now | after G1, G3, G4, G5 |
|---|---|---|
| fit + chains + polish + amorphous | 3.7 | 3.7 |
| enumeration | 3.95 | 0.06 |
| table builds | 25.55 | 10.55 |
| FFT screen | 12.04 | 2.26 |
| refinement | 15.15 | 8.93 |
| total | 60.4 | 25.5 (2.4x) |

On the user's machine (section 10, CPU work 71.8 s: tables 24.3, screen 8.3,
refinement 25.7, enumeration 4.6, the rest 9.0) the same ratios give about 36 s, a
little over 2x, and the warm-cache end-to-end run, which no longer pays the tables,
gains from the screen, the refinement and the enumeration alone.  After all five the
funnel here would be table builds 41%, refinement 35%, screen 9%, and the next levers
are the ones each finding names: a compiled or fused table kernel and fewer evaluated
grid points; the closure solve and then the kernel row itself in the refinement.

### 15.6 A finding about the tests

Ten tests fail on the untouched head in this container and pass on the machine the
branch is developed on: `test_ewald::test_default_packer_is_bit_for_bit_the_truncated_one`,
`test_forcefield::test_default_simpleff_is_bit_for_bit_the_illustrative_potential`
(four cases), `test_forcefield::test_third_order_terms_obey_reflection_with_reversal[pvdc]`,
`test_mechanics::test_default_packer_energies_are_bit_for_bit_unchanged` (four cases),
plus the eight `byte_identical` digest tests deselected above.  All compare `==`
against float64 literals or SHA digests recorded on one machine and differ here in
the last one or two digits (`-8.645931161267757` against `...756`).  They guard real
invariants, but as written they test the BLAS build rather than the code; recording
them to a tolerance of a few ulps, or per platform, would keep the guard and lose the
false alarms.  Every prototype above reproduced exactly these ten failures and no
others.  On the head this section was pushed to (`0586a35` plus this document), the
suite here gives 589 passed and those ten plus one more:
`test_mechanics::test_direct_and_converse_piezoelectric_routes_agree`, added the same
day in `dd86575`, asks the two routes to agree within 1% and they agree within 1.12%
on this box (relative difference 0.0112).  It is a tolerance test on a
finite-difference route, not a golden value, and no source changed here, so it is
recorded rather than diagnosed.

*Addressed 2026-10-08.* The literal comparisons are held to 1e-12 relative (jittering every
input by 4 ulps moves the recorded energies by at most 9e-14), and the eight digest tests
compare against the arrays the digests were taken of, stored in `tests/data/pendant_builds.npz`
and checked against the digests (charges exactly, coordinates to 1e-12 of their scale); the
comparisons between two computations in one process stay exact. Under an emulated libm whose
sin and cos are off by up to 2 ulps the old forms of these tests fail and the new ones pass.
One correction to the paragraph above: `test_third_order_terms_obey_reflection_with_reversal[pvdc]`
compares no literal. PVDC's trans basin is a symmetric double well whose minima at -120 and
+120 deg differ by 1e-13 kcal/mol, so `fit_ris`'s `argmin` picks between them by rounding; this
build resolves the two bond types opposite ways and `adapt_angles` puts T at 180, a build that
resolves them the same way puts T at +/-120 and the measured relation fails by the 50 kcal/mol
cap. No other decision in the test is that close (every other basin minimum leads by
0.02 kcal/mol or more), so this is the likely failure there; the test now measures PVDC at
the nominal state angles. The tie itself is in `fit_ris` and would move any rigid PVDC fit
the same way on such a build; it is recorded, not fixed.

### 15.7 Measured on the merged tree (2026-11-04)

G4, G5, G1, G3 and G2 are applied on the branch in that order (commits `31a2456`,
`bdbcb8f`, `e364e0b`, `1700515`); G3 and G2 were merged by hand onto the files G4 and
G1 had changed, and G6 stays a patch.  `examples/profile_funnel.py` was run on the
branch head and, from a second worktree with its own `src` on `PYTHONPATH`, on the
untouched base `0cb880e`: the box of 15.1 (four cores, NumPy 2.5.3), single process, no
table cache, the two runs back to back on an otherwise idle machine.  The refinement
row is under cProfile in both runs, as it was in 15.1, and is given again without it.

| stage | base `0cb880e` | merged | ratio | projected in 15.5 |
|---|---|---|---|---|
| RIS fit (2,881 single points) | 0.46 | 0.35 | | |
| enumeration (120 candidates, period <= 8, k = 60) | 3.00 | 0.05 | 60x | 0.06 |
| periodic chain construction (5) | 0.10 | 0.03 | | |
| chain-pair table builds (5, screen grid, 4 worker processes) | 23.36 | 10.92 | 2.1x | 10.55 |
| FFT screen (5) | 12.92 | 2.00 | 6.5x | 2.26 |
| exact-kernel polish (4 starts per conformation) | 3.32 | 3.34 | 1.0x | |
| torsion + cell refinement (5), under cProfile | 18.73 | 8.22 | 2.3x | |
| total, with that refinement row | 61.90 | 24.91 | 2.5x | |
| refinement (5) without the profiler | 14.9 | 9.4-10.9 | 1.4-1.6x | 8.93 |
| total, with that refinement row | 58.1 | 26.1-27.6 | 2.1-2.2x | 25.5 (2.4x) |

The profiler-free refinement row is the five refinements of the same polished starts,
timed on their own: 14.9 s on the base (15.2 s in 15.1); on the merged refinement code
9.4-9.9 s as the minimum of three runs in each of two orderings against G1 alone (which
took 9.7 s in both: G2 buys no wall time over G1 on this box), and 10.9 s in a single
run of the branch head.  Every stage landed where 15.5 projected it except the
refinement, where the projection was G1's own measurement and holds; G2's projected
"below 8 s" did not happen, for the reason given under *what G2 changed* below.

Per candidate, base -> merged (seconds; the refinement column under cProfile):

| candidate | atoms | table | screen | polish | refine | evals / vars / param. |
|---|---|---|---|---|---|---|
| TG+TG+TG-TG- | 24 | 6.69 -> 3.02 | 2.87 -> 0.54 | 0.90 -> 0.92 | 4.85 -> 1.56 | 70 -> 66 / 10 / linegroup |
| TTTG+TG- | 18 | 4.67 -> 2.31 | 2.66 -> 0.42 | 0.98 -> 0.92 | 3.89 -> 2.27 | 128 -> 143 / 13 / free |
| TTG+TG-G-TG+ | 24 | 7.14 -> 3.30 | 2.96 -> 0.37 | 1.12 -> 1.18 | 8.77 -> 3.86 | 211 -> 218 / 15 / free |
| TG+TG- (alpha) | 12 | 3.19 -> 1.51 | 2.27 -> 0.34 | 0.18 -> 0.22 | 1.04 -> 0.45 | 29 -> 31 / 8 / linegroup |
| TT (beta) | 6 | 1.67 -> 0.78 | 2.16 -> 0.33 | 0.14 -> 0.11 | 0.17 -> 0.08 | 10 -> 10 / 6 / linegroup |

**The refinement is now the kernel.**  cProfile cumulative seconds summed over the five
candidates, base -> merged: `CrystalPacker.energy_and_grad` 5.5 -> 5.8 (unchanged, it is
the same code); `build_backbone` 5.3 -> 0.07; `repeat_chains` 9.4 -> 0.36; the closure
Newton `LineGroup._close` 3.1 -> 0.42; `_jacobian` 1.0 -> 1.7 (it is now the closed-form
chain Jacobian of every candidate rather than the differenced closure Jacobian of the
line-group ones).  The geometry that was three quarters of the refinement in 15.1 is
now under a third, and the kernel row is about 70% of what is left; the next lever is
the one 15.5 named, the kernel row itself.

**Tests.**  The merged tree gives 630 passed, 5 skipped, 8 deselected (the
`byte_identical` digest tests of 15.6, which depend on the BLAS build) and no failures
in 235 s; the base gives 615 passed on the same box, and the fifteen new tests are
`tests/test_chainjac.py`.

**What moved numerically.**  G4, G5 and G1 are bit-identical or exact (15.3), and G3
changes table values by at most 8.7e-4 kcal/mol over the five screen tables (clear any
on-disk `$POLYFIND_TABLE_CACHE`; the build ran at its default 65,536-point chunk, and the
re-sweep on the target machine that 15.5 asks for is still to be done).
Refined from the same polished starts, the base and the branch with everything but G2
give the same end points to the ninth decimal, so every movement of a refined value is
G2's, and G2 moves them as follows (same starts, base -> merged, per monomer):

| candidate | param. | d objective | d lattice E | d bend | max d torsion (deg) | max d cell |
|---|---|---|---|---|---|---|
| TG+TG+TG-TG- | linegroup | -2.1e-9 | -2.9e-5 | +2.9e-5 | 1.0e-3 | 0 (phi2 reported as 0 and 360) |
| TG+TG- | linegroup | -1.9e-10 | -9.8e-6 | +9.8e-6 | 5.1e-5 | 1.2e-5 |
| TT | linegroup | +2.3e-13 | -2.0e-8 | +2.0e-8 | 2.8e-14 | 1.3e-9 |
| TTTG+TG- | free | +1.7e-5 | -9.6e-5 | +1.1e-4 | 7.6e-3 | 5.2e-3 |
| TTG+TG-G-TG+ | free | -1.9e-4 | +7.6e-3 | -7.7e-3 | 0.35 | 0.19 (setting angle, deg) |

The line-group end points are unchanged for practical purposes.  The two free-mode
candidates stop at different places along one flat valley: in TTG+TG-G-TG+ the lattice
energy and the bend restraint move by 8e-3 each in opposite directions and the full
objective (energy, commensurability penalty, bend, caps) is 5.9e-4 higher at the new
stop.  Continuing the optimisation from either stop with the analytic gradient reaches
the same point (objective -4.37772 against -4.37732 and -4.37673 at the two stops;
torsions agree to 2e-4 deg, objective to 2e-10), 0.36 and 0.58 deg of torsion further
down the valley: both are early stops of L-BFGS-B's relative-reduction test, not two
minima.  They differ because the previous route's finite-difference geometry gradient
is off by 13% of the gradient at the TTTG+TG- stop and by 2.7% at the TTG+TG-G-TG+ one
(where the gradients are 2.5e-4 and 1.4e-2, against 6e-6 and 2e-5 of gradients of 7 and
19 at the starts), while the analytic gradient agrees
with a finite difference of the whole objective to 5.7e-7 of the gradient scale over
380 components (18 cases of polymer, sequence, parametrisation and field) and with the
previous route to 3.2e-6 where the gradient is large.  `refine_crystal(gradient=
"fd-geometry")` keeps the previous route and reproduces the base end points (to 1e-8
kcal/mol for the line-group candidates, to the printed digits for the free ones).  In
the full funnel run above, whose starts come from the G3 tables, TTG+TG-G-TG+ refines
to -5.2941 kcal/mol per monomer against -5.3075 on the base, the same trade-off from a
slightly different start; a pinned refined value of a free-mode candidate should be
compared on the objective, or re-recorded, not held to the old lattice energy.

**What G2 changed, and why it bought no time.**  Per refinement evaluation the chain
build went from 7.3 backbone calls and 98 rows (line group) or 17-21 rows (free) to
3.2 calls and 3 rows or 1.1 calls and 1 row, and the oligomer build from 1.06 calls per
evaluation to 0.06; the finite differences G2 removes were already cheap once G1 had
made each build cheap, and what is left per evaluation is the kernel row, which G2 does
not touch.  G2 is kept for the gradient quality above and the 30x fewer builds, which
matter more on a slower chain build (a larger repeat, or the complete placed-cell
owner of the `vector-repeat-merge` branch) than they do here.

**GPU.**  Nothing in this section ran on a GPU; 15.4 stands and G6 stays a patch.

## 16. The large-cell response path (2026-11-05): the 120-atom cell an order of magnitude faster

Scope: the E1 prediction path of `examples/e1_predictions.py` on its largest system, the
all-trans 70/30 VDF/TrFE model cell (two chains of 60 atoms, 120 atoms and 20 monomers in the cell), on the fitted valence potential with
Ewald electrostatics and induced dipoles: build and relax the reference, then the
stiffness, the clamped and relaxed-ion dielectric tensors, the proper and converse
piezoelectric routes, the electrostriction routes and the true dielectric tensor at two
strained cells.  Same 4-core box as section 15, no GPU.  Before this section the full job
on that cell did not finish inside the two hours the background runner allows, and the
reduced job (without the electrostriction and converse routes) took 1,410 s.  The three
smaller systems of the E1 set (beta-PVDF, PE, alpha-PVDF) were never the problem: 43, 119
and 180 s each before any of this.

### 16.1 Where the time went

The stages of the path, serial, under cProfile (which inflates Python-call-heavy code
but not the two rows that matter here):

| stage | seconds |
|---|---|
| `refined_reference` (RIS fit, enumeration, pack) | 28.3 |
| `relax_reference` | 43.7 |
| `elastic_constants` (6 strained, relaxed states) | 90.5 |
| one tight strained state, no field | 15.3 |
| one tight strained state, field along x | 45.0 |
| `dielectric_tensor` (6 relaxed states) | 91.0 |
| total profiled | 313.8 |

Cumulative time, same run: `ewald.terms` 261.1 s of 313.8 (83%), of which the
three-operand `einsum` that assembled the reciprocal-space dipole tensor 152.9 s (49%
of everything) and the routine's own time 80.8 s plus 19.4 s in `damping`; `_polarize`
266.9 s; the Lennard-Jones and valence pair terms (`_pair_energy_and_dv`) 29.4 s.
`terms` was called 302 times at 0.865 s a call, 21 L-BFGS-B minimisations in all.  The
response path is one Ewald kernel evaluated a few hundred times, and the kernel was
dense: every real-space kernel was evaluated on every pair-image combination and masked
afterwards, and on this cell about 98% of those combinations lie outside the 8 A cutoff.

### 16.2 Four changes

**Neighbour-list real space and a BLAS reciprocal tensor** (`ewald.py`).  One pass
over the distances finds the pair-image list; every kernel (charge, dipole, Thole
damping, their gradients) runs on the list, sums over images and sites are
scatter-adds, and the two lattice-derivative moments accumulate as 3 x 3 matrices.  The
reciprocal dipole tensor, an einsum over `(N, N, K, 3, 3)`, becomes six `(N, K) @ (K, N)`
products.  Every output agrees with the previous routine to 1e-13 on the beta cell and
on this cell (energy, dipoles, field, tensor, coordinate, lattice and charge
gradients), and the 63 Ewald and polarizable tests pass.  On this cell a
field-and-tensor pass goes 1,289 -> 171 ms and a polarization solve with gradients
2,058 -> 947 ms, the remainder of which was one 360 x 360 linear solve.

**One BLAS thread per process** (`polyfind/__init__.py`).  That 360 x 360 solve took
400-1,000 ms with OpenBLAS's default thread count and 2 ms with one thread: the pool's
spin-waiting threads contend with each other and with every other process on the box,
and polyfind's linear algebra is many small operations, none of which a thread pool
helps.  The package now limits BLAS to one thread at import (environment variables for
libraries loaded later, `threadpoolctl` for those already loaded; `POLYFIND_BLAS_THREADS`
overrides).  With both changes, on this cell: a polarization solve with gradients 2,495 ->
168 ms, `relax_reference` 44 -> 7 s, `elastic_constants` 90 -> 14 s.

**Independent relaxations in parallel processes** (`mechanics.parallel_map`).  The
states behind the tensors do not depend on one another: 6 strained states for the
stiffness and the proper piezoelectric tensor, 6 field states for the relaxed-ion
dielectric tensor, 6 Newton chains for the converse route, 6 field states for each
strained true dielectric tensor, and in the electrostriction routes 55 (one at zero
stress and field, 12 at fields for the strain and the free susceptibility, 36 at
+/- stress and +/- field for the polarization route, 6 for the clamped susceptibility).
Each group now runs in forked worker processes, one BLAS thread each (the reason the
previous change had to come first), with what the jobs need inherited through the
fork rather than pickled; `POLYFIND_WORKERS` sets the count (default: the CPU count),
and the map is serial for one worker, one item, inside a worker (no nesting) or off
Linux.  The results equal the serial loops' (`tests/test_parallel.py`: stiffness to
1e-9 GPa, tensors to 1e-10, the electrostriction coefficients exactly), because the
jobs are the same function calls in another process.  What stays serial: the reference
build and relaxation, the zero-stress state, and the iterates inside each Newton
(`free_strain`, `state_at_stress`), which depend on the previous one.

**No Coulomb work in the pair kernel under Ewald** (`pack.py`, after the profile of 16.3).
The Lennard-Jones pair kernel also carries the damped-shifted-force Coulomb term of the
non-Ewald setting; under Ewald the charge product it multiplies is zeroed (every
electrostatic term comes from `ewald.terms`), but the `erfc` and its derivative were still
evaluated on every pair-image combination and multiplied by zero: 16.5 s of the 68.9 s
serial profile below. The kernel now skips them when the packer has an Ewald sum. Adding an
exact zero is bit-identical to leaving it out, and the A/B confirms it: energies and
gradients of the beta and copolymer packers agree with the previous module to the last bit
(fixed-charge and polarizable), and the 120-atom `energy_and_grad` goes 269 -> 217 ms. The
end-to-end numbers of 16.3 were measured before this change.

### 16.3 End to end

Full E1 job on this cell (`--systems vdftrfe`, electrostriction and converse routes
included), wall time, same box.  The crystal-stage fit of `docs/CRYSTAL_STAGE.md` was
running on another core during the first two measured rows and not during the last:

| version | wall |
|---|---|
| before this section, reduced job (no electrostriction, no converse) | 1,410 s |
| before this section, full job | > 7,200 s (killed at the runner's cap) |
| Ewald + BLAS changes, serial | 1,139 s (build, tensors and electrostriction routes 1,014 s; two strained true dielectric tensors 69 and 56 s) |
| all three changes, 4 workers, the fit holding about 1.5 cores | 512 s (444 s; 34 and 33 s); results identical to the serial run in every one of the 40 stored quantities |
| all three changes, 4 workers, idle box | 359 s (304 s; two strained tensors 30 s each); results identical to the serial run |

The same profile script after the three changes, serial (`POLYFIND_WORKERS=1`, so the
profiler sees every relaxation), same box:

| stage | before | after |
|---|---|---|
| `refined_reference` | 28.3 | 29.1 |
| `relax_reference` | 43.7 | 7.1 |
| `elastic_constants` | 90.5 | 13.1 |
| one tight strained state, no field | 15.3 | 2.4 |
| one tight strained state, field along x | 45.0 | 5.2 |
| `dielectric_tensor` | 91.0 | 12.1 |
| total profiled | 313.8 | 68.9 |

Cumulative, after: `ewald.terms` 20.1 s (was 261.1) in 290 calls of 69 ms (was 865 ms), the
einsum 3.5 s (was 152.9), `_polarize` 21.5 s (was 266.9), the 360 x 360 solves 0.33 s for 357
calls. The largest kernel is now the Lennard-Jones and valence pair routine
`_pair_energy_and_dv`, 29.4 s, unchanged in absolute terms (it was 29.4 s before too), of
which 16.5 s is the complementary error function and its derivative (`_erfc_pos`,
`_derfc_pos`); and the largest stage is the reference build, `refined_reference`, 29 s, of
which `pack` (the chain-pair tables with their worker pool, the FFT screen and the polish) is
20 s. Everything the response path repeats is 4.6x faster serially; the one-off build is
untouched and is now 42% of a serial profile.

### 16.4 What is left

- **The pair kernel.** With the Coulomb work gone, `_pair_energy_and_dv` still evaluates
  the Lennard-Jones terms on every pair-image combination, dense, exactly as the Ewald real
  space did before this section, and is the largest kernel row left (about 13 s of the
  serial profile). The same neighbour-list treatment would take most of it away; it is the
  next kernel row on this path.
- **The reference build.** 29 s per system, serial with respect to everything after it,
  and 42% of what is left on one process; on four workers it is a tenth of the full E1 job
  on this cell and most of a job that stops at the tensors. Its table build already uses a
  worker pool; the screen and the polish of the best starts run one after another, and the
  polish is L-BFGS-B relaxations of the same kind the rest of the path now runs in parallel.
- **Amdahl on the E1 path.** The serial parts of the end-to-end job are the build and
  reference relaxation (36 s), the zero-stress state, and the iterates inside each Newton
  (`free_strain`, `state_at_stress`), which depend on the previous one; the biggest group,
  the 55 states of the electrostriction routes, is 49 Newton chains of a few relaxations
  each plus 6 single relaxations, so it scales to four workers but not much beyond without
  splitting the chains.
- **GPU.** Nothing here ran on one; the kernels are now neighbour-list scatter-adds and
  small dense products, which is the shape section 14 said a GPU would want, and the
  decision there stands until a GPU is on the box.

### 16.5 The pair kernel on a list: measured, not kept

A neighbour-list form of the Lennard-Jones kernel of the gradient path (the pair potential
evaluated on the combinations inside the cutoff, found by one `nonzero` over the dense
distance array, the gradient assembled by `bincount` scatter-adds) was built and A/B'd:
energies equal to the dense kernel's, gradients to 4e-14, and 183 -> 169 ms per 120-atom
`energy_and_grad` (8%).  It is not kept.  The gain is small because the `(images, N, N)`
distance arrays are still formed densely, and that, not the pair potential, is most of
what the kernel costs once the Coulomb work is gone (16.2); and the different summation
order breaks the invariant `tests/test_pack.py` guards, that `energy_and_grad`'s value is
bit for bit the one `energy` returns for the same row, which the rest of the module leans
on.  Two per cent of the end-to-end time is not worth that invariant.

What a gradient call is made of now, from cProfile of ten `energy_and_grad` calls on the
copolymer packer (induced dipoles, Ewald), 1.55 s: `ewald.terms` 1.25 s cumulative (two calls
per energy, 0.90 s own time), the reciprocal-space `einsum`s 0.23 s, the packer's own time
0.16 s, the list kernel 0.004 s.  A gradient call is 80% Ewald, and the next lever on this
path is inside `terms`: the real-space list is rebuilt from a dense distance array at every
call (a Verlet skin would reuse it across the small steps of a relaxation), and the
reciprocal-space structure factors are recomputed for every k-vector.  Neither was taken
here; the end-to-end numbers of 16.3 stand as the state of the speed stage.

### 16.6 Inert sites out of the electrostatics (2026-11-07)

The belt sites of `polymers.belt_fluorine` (pseudo-element "X", `docs/CRYSTAL_STAGE.md` 5.7) carry
no charge and, until now, a 1e-6 A^3 placeholder polarizability so that the dipole solve's
diagonal stayed finite; they rode through the Ewald sums and the induced-dipole solve like atoms.
On the 70/30 copolymer cell that is 92 sites beside 120 atoms, and one `energy_and_grad` call
cost 0.64 s against 0.17 s for the plain chain, 56% of it in `ewald.terms` (cProfile of three
calls: `terms` 1.07 s cumulative, the Lennard-Jones kernel 0.43 s, the gradient assembly 0.27 s).

A zero charge adds nothing to any term of the charge sum, and a site without a polarizability has
no unknown in the solve, so the packer now leaves every site that is chargeless and unpolarizable
out of both (`CrystalPacker._es_idx`, set once from the topology; the Ewald sums and the
`(3N)^3` solve run on the active sites and the inert ones get the zero dipoles and gradients they
would have had; `_set_charge_tables` refuses a charge that later lands on one).  The placeholder
is now a true zero (`polarizability.VDS98["X"] = 0.0`; `Polarizable.per_atom` accepts zero as "no
induced dipole on this site").  Measured: the belt copolymer call 0.64 -> 0.42 s, `ewald.terms`
1.07 -> 0.38 s per three calls, i.e. back to the plain chain's 0.35; the stage-7 fit point
re-evaluated moves by at most 4.3e-6 relative (the placeholder dipoles' energy), and the gradient
matches finite differences as before (`tests/test_belt_sites.py::test_inert_sites_are_left_out_of_the_electrostatics`;
the Ewald total over all sites equals the total over the active ones to 1e-9).

What is left of the belt cost is the pair kernel: 0.46 s per three calls against 0.08 on the plain
chain, from the dense `(images, N, N)` arrays over 212 sites of which the 92 inert ones interact
with hydrogen alone.  That is the finding of 16.5 again, and the same answer: a pair list would
take it out at the price of the bit-for-bit `energy == energy_and_grad` invariant.  The belt cells
are 2.5x the plain chain's call now, from 3.8x.  End to end, the belt2 copolymer E1 cell
(`examples/e1_predictions.py --systems vdftrfe`, build, tensors and electrostriction) runs in 776 s
on a quiet machine, 601 s of it build and tensors, against 1,356 s for the belt set alone before
the change (1,719 s under load): 1.75x, with the cell, polarization and permittivities unchanged
and the stiffness block within 2e-4 relative, the relaxation's noise floor.


## 17. The many-chain cell on a device (2026-11-13): the port, and what it has not been measured on

Scope: one energy-and-gradient evaluation of a `polyfind.mixedcell` cell of sixteen and more
chains of the ten-monomer PVDF repeat (960 electrostatically active sites per sixteen chains),
the unit of work of every relaxation, band and saddle search of `docs/CANDIDATE_MATERIALS.md`
section 8.  The PI can rent an A100 and can buy four to six.  17.1 to 17.5 are the first port
(the Ewald sums and the dipole solve) and its measurement plan, written before any device time
existed; 17.6 is that port measured on a rented A100, which showed that it was not enough; 17.7
is the rest of the evaluation, moved after that; 17.8 is the whole measured again.  This
container has four cores and no GPU: every device number below is from a pod, and its files are
in `docs/predictions/`.

### 17.1 Where the time goes on the host

Measured here with `examples/ewald_device_benchmark.py --device none`, each evaluation on a
geometry the Ewald object had not seen.  The package holds the BLAS to one thread
(`polyfind/__init__.py`; `POLYFIND_BLAS_THREADS` lifts it), so both settings are given; the
sixteen-chain rows are the range of three runs each, taken in alternation:

| chains | sites | k vectors | BLAS threads | whole evaluation | tensor and field | solve | gradient |
|---|---|---|---|---|---|---|---|
| 16 | 960 | 14,486 | 1 | 5.3 to 6.5 s | 3.0 to 3.1 s | 0.4 to 0.5 s | 1.0 to 1.4 s |
| 16 | 960 | 14,486 | 4 | 4.4 to 4.7 s | 2.2 to 2.4 s | 0.2 s | 0.8 to 0.9 s |
| 32 | 1,920 | 28,920 | 1 | 34.5 s | 19.2 s | 3.1 s | 10.9 s |
| 32 | 1,920 | 28,920 | 4 | 24.7 s | 11.9 s | 1.4 s | 7.5 s |

(The three steps are timed apart on the evaluation's own arguments.  This container's speed is
not steady, as the ranges show.)

**A correction.**  RESUME's note of this morning, and the answer given to the PI on it, said
that four BLAS threads change nothing (26.7 s against 28.2 at thirty-two chains) because the
routine is elementwise.  That comparison did not vary the threads: the package sets the loaded
BLAS to one thread when it is imported, and `OPENBLAS_NUM_THREADS` in the environment does not
undo that in a script that has imported NumPy first, as every script here does.  Both runs were
at one.
Varied properly (`POLYFIND_BLAS_THREADS=4`), four threads make the whole evaluation 1.2 to 1.4
times faster: the solve halves, the six tensor products gain a quarter to a third, the
elementwise arithmetic gains nothing.  A large-cell job that has the machine to itself should
set it; jobs run side by side should not.

The dipole step is `CrystalPacker._polarize`: `Ewald.terms` for the permanent field and the
`(N, N, 3, 3)` dipole tensor, the dense solve for the induced dipoles, `Ewald.terms` again for
the gradient at the solved dipoles, and the bonded-exclusion corrections chain by chain.  Inside
`terms` the time is the reciprocal sum (line profile, sixteen chains, one thread, 3.2 s for the
tensor call: the six symmetric products that assemble the tensor 1.3 s, `cos` and `sin` of the
`(N, K)` phases 0.5 s, the field 0.3 s; 1.2 s for the gradient call, all of it elementwise
arithmetic on `(N, K)` arrays and products with the k vectors).  The real-space neighbour list
does not show.

What the port leaves on the host, from cProfile of one whole evaluation at one thread: 0.9 s
of 6.7 at sixteen chains (the short-range pair sums 0.39 s, the chains' geometry and its
Jacobians 0.26 s, the bonded-exclusion corrections 0.19 s) and 1.9 s of 33.7 at thirty-two.

The work grows as the cube of the cell: the number of k vectors `K` is proportional to the
volume, so the tensor is `N^2 K ~ N^3`, like the solve; the gradient call is `N K ~ N^2`.

### 17.2 What moved

`Ewald.terms` is written once against an array module `xp` (`ewald.py`).  With NumPy, the
default, it is the routine it was: every result bit for bit equal to the previous version's on
six reference cases (the sixteen-chain cell's tensor-and-field, gradient-with-dipoles and
charge-only calls, and a small skewed cell with every flag, charges alone, and the undamped
tensor).  With a device module the arrays that grow with the cell are built and kept there (the
pair separations, the neighbour list and its kernels, the `(N, K)` phases, the tensor); the 3 x 3
lattice algebra and the list of k vectors are made on the host and sent; the per-site results
come back, and the tensor comes back only if the caller does not ask to keep it
(`device_tensor=True`).  `_polarize` (`pack.py`) asks, builds the `3N x 3N` matrix on the
device, solves there, and brings back the dipoles; the `(N, N)` table of Thole lengths goes
over once per packer.  Per sixteen-chain evaluation 2.0 MB goes to the device in 26 pieces
(coordinates, k vectors and their weights) and 69 kB comes back in 13, the largest of them a
per-site array of 23 kB (counted by the stand-in of 17.4).

The switch is its own: `POLYFIND_EWALD_DEVICE=cuda`, or `Ewald(..., xp=cupy)`.  It is not
`POLYFIND_DEVICE`, which turns itself on wherever CuPy finds a device and runs the batched screen
in float32.  This path is float64, serves cells of a few dozen sites as well as of thousands,
and is the arithmetic of every stored large-cell number, so it moves only when asked and only
for cells of at least `Ewald.device_min_sites` (256) sites.  A device sums in another order,
and CuPy's weighted `bincount` is an atomic add whose order is not fixed, so device numbers
agree with the host's to rounding and are not reproducible bit for bit from run to run.

One trap found on the way and closed: on a machine with a GPU the packer's batched kernels go
to it in float32 unasked, and a mixed cell built on such a packer would have run its float64
host loops through them, block by block across the bus, at seven digits.  `MixedModel` now
refuses such a packer and says what to set (`POLYFIND_DEVICE=cpu`).

What did not move: the short-range pair sums (`mixedcell.placed_energy_and_grad`, a host loop
over chain pairs), the chains' geometry and its Jacobians, and the bonded-exclusion
corrections: the 0.9 s and 1.9 s of 17.1.  They grow linearly with the cell, and they are the
ceiling of this port: with the three ported steps at zero cost an evaluation still takes that
long, against 4.4 s and 24.7 s for the host on four threads.  Five times at sixteen chains and
thirteen at thirty-two is therefore the most a card can give over this host until those move
too.  The whole evaluation on the host is the same arithmetic as before the port: energy and
gradient of the two-, four- and sixteen-chain cells are bit for bit those of the previous
commit.

### 17.3 Read from CuPy's source, not assumed

Against the 14.2.0 release (the current one; wheels `cupy-cuda12x` and `cupy-cuda13x`, NumPy
2.0 to 2.5):

* `einsum` with batch axes and no cuTENSOR becomes `matmul` of stacks of matrices: the
  neighbour-distance step `einsum("mijc,mijc->mij", D, D)` would have been one 1 x 3 by 3 x 1
  product per pair-image.  On a device it is `(D * D).sum(-1)`.
* `a @ a.T` is a general `gemm`; nothing routes it to `syrk`.  The six products of the
  reciprocal tensor were arranged as symmetric rank-k updates because the host's BLAS forms
  those at half the cost (16.2); on a device they are general products, twice the arithmetic
  needed.  `cupy.cublas.syrk` exists and fills one triangle.  Not used yet: the benchmark
  times it beside the general product on the cell's own shapes.
* `linalg.solve` is LU (`gesv`) on copies of both operands, and by default does not report a
  singular matrix.  `cupyx.lapack.posv` (Cholesky) exists; the matrix is symmetric positive
  definite at any state short of a polarization catastrophe.  Also timed, not used.
* `bincount` with weights returns float64 from an `atomicAdd` kernel and synchronises twice.
* An array has no `setflags` and its flags are read-only: the kept phases cannot be marked
  read-only on a device as they are on the host.
* A reduction returns a zero-dimensional array; `float()` of it is a transfer; a NumPy array
  met in arithmetic, or converted implicitly, is a `TypeError`.

### 17.4 Tested without a device

`tests/device_double.py` is an array module that is NumPy underneath and has those last
rules on top: its arrays are not NumPy arrays, do not become one except through `get()` or
`float()`, do not mix with one in arithmetic, in a function or as an index, return
zero-dimensional arrays from reductions, and have no `setflags`; and it counts the bytes that
cross.  `tests/test_ewald_device.py` runs the device path on it: every term and gradient with
both boundary conditions, both Thole forms, charge penetration and a charged cell agrees with
the host to 1e-12 of the result's size (the only arithmetic that differs is the two places
above); the tensor stays only when asked; nothing of the size of `N^2` or `N K` crosses in a
call; a small cell stays on the host whatever the module; the phases are kept per module; an
`Ewald` object pickles without its module; and a whole four-chain evaluation, solve included,
is the host's to 1e-12 with the Thole table sent once.  The benchmark itself runs on the
double (`--device double`): at sixteen chains the largest relative difference from the host
over every result is 3e-15.

What this cannot show is that CuPy's own functions behave as read, and how long they take.
`tests/test_gpu.py` holds the same comparisons for a real device; they are skipped here and are
the first thing the pod script runs.

### 17.5 The measurement

`bash examples/runpod_ewald_bench.sh` on a fresh pod: installs, checks that CuPy reaches the
device, cuBLAS, cuSOLVER and the kernel compiler (taking NVIDIA's runtime wheels if the image
has no libraries CuPy can find), runs the device tests, then
`examples/ewald_device_benchmark.py` at 16, 32 and 64 chains with the host beside it, and at 128
on the device alone.  For a cell the host is not run on, the device is held to two things that
need no host run: the untouched tiling must have the stored two-chain cell's energy per chain
(1e-13 kcal/mol on the host at sixteen chains), and the gradient must agree with the device's
own fourth-order energy difference (which has a floor of about 2e-5 at sixteen chains from pairs
crossing the real-space cutoff, on the host as anywhere).

What each number it reports decides:

| number | decides |
|---|---|
| largest difference from the host, per size | whether the port is right on a real device |
| whole evaluation, device against the host at one thread and at every core | what one card buys a relaxation, which is sequential |
| the share left on the host | whether the pair sums are the next thing to move |
| memory pool at each size | which cells fit a 40 GB card and which need 80 |
| general product against `syrk`, LU against Cholesky | the next two changes on the device, or not |
| first call | the one-time kernel compilation a short job pays |
| `--all-devices` (more than one card) | whether cards add up for independent evaluations |

A farm of cards multiplies what is independent: the images of a band, the columns of a
difference Hessian, separate starts, the rows of a fit.  It does nothing for one relaxation,
whose evaluations follow one another; there the gain is one card's, capped by what 17.1 says
stays on the host.  Past about a hundred chains the cube is the issue and the remedy is not hardware: a
matrix-free solve of the dipoles on a particle mesh is `N log N` per iteration.  None of that is
designed here.

One limit of the present code on a device: a process that has used the device cannot fork
workers that use it (`examples/kink_migration.py --workers`); independent evaluations on
several cards want one process per card started fresh, which is how `--all-devices` runs.

### 17.6 Measured on an A100: the three steps alone (2026-11-13)

Two rented pods, one NVIDIA A100-SXM4-80GB each (driver 595.91.07, CuPy 14.2.0, NumPy 2.2.6, host
AMD EPYC 7742 with 13 cores allotted of the 128 the container sees), ran
`examples/runpod_ewald_pod.sh` at commits 2fe1c9e and 73b556f.  Files, logs and the tables below
(`examples/ewald_device_summary.py` on the JSON; nothing typed in):
`docs/predictions/ewald_device_a100_20261113/`.  The device tests passed on both (17 passed).

One whole evaluation, energy and gradient:

| chains | sites | k vectors | host, 1 thread | host, 8 threads | host, 13 threads | device | of it outside the three steps | device against fastest host | largest difference | card memory |
|---|---|---|---|---|---|---|---|---|---|---|
| 16 | 960 | 14,486 | 15.8 s | 2.80 s | 2.73 s | 0.51 s | 0.36 s (70%) | 5.3x | 7.3e-15 | 1.3 GB |
| 32 | 1,920 | 28,920 | 23.2 s | 11.7 s | 10.8 s | 1.14 s | 0.77 s (68%) | 9.5x | 8.3e-15 | 5.1 GB |

The three ported steps, timed apart (fastest host run against the device):

| chains | tensor and field | solve | gradient | the three together | device against host |
|---|---|---|---|---|---|
| 16 | 1.38 s / 20.7 ms | 0.44 s / 17.3 ms | 0.38 s / 9.6 ms | 2.20 s / 47.6 ms | 46x |
| 32 | 5.69 s / 92.3 ms | 3.04 s / 40.4 ms | 1.21 s / 19.6 ms | 9.94 s / 0.15 s | 65x |

Side measurements: the general product against the rank-k update 1.8 / 1.6 ms and 11.1 / 8.5 ms
(14.5 and 19.3 TFLOP/s in double precision); the LU solve against Cholesky 16.9 / 3.1 ms and
39.5 / 8.6 ms, the two solutions agreeing to 2e-15 and the matrix exactly symmetric.

What it says:

* **The port is right on a real device**: 8e-15 of each result's size against the host, the
  tiling check 3e-13 and 7e-13 kcal/mol, the gradient against the device's own energy difference
  2e-5 and 2e-6 (the cutoff floor of 17.5).
* **The three steps are no longer the evaluation.**  They run 46 and 65 times faster than
  thirteen host cores, and the evaluation 5.3 and 9.5 times, because 70% of the device
  evaluation is what was left on the host: 0.36 s of 0.51 s and 0.77 s of 1.14 s, the card idle
  meanwhile.  (More than 70%, in fact: the dipole step still held the bonded-exclusion
  corrections chain by chain on the host, 0.11 s of its 0.16 s at sixteen chains.  The three
  steps themselves are 48 ms.)  The PI saw it as 0% utilisation on the card's own meter.
  17.7 is the consequence.
* **Cholesky, yes; the rank-k update, no.**  The solve is a third of the three steps and
  Cholesky is a fifth of it; the rank-k update gains a tenth to a quarter of a step that is 2 to
  11 ms.
* **Memory**: 1.3 GB of pool at sixteen chains and 5.1 GB at thirty-two (1.8 and 5.6 GB in use
  on the card), growing as the square of the cell: sixty-four chains want about 20 GB and a
  hundred and twenty-eight about 80.
* **Host times on a rented pod are noisy** and are not the reference for anything but this
  comparison: the container is held to its allotment by a quota, and one run of the first pod
  on "all cores" (256 BLAS threads under a 32-core allotment) was slower than one thread, which
  is why the benchmark now reads the allotment.  Sixty-four chains, host only, from the second
  pod's log: 462.6 s at one thread and 59.4 s at eight.  The second pod was stopped there, before
  its sixty-four-chain device row, because the card was idle through those minutes.

### 17.7 The rest of the evaluation: the batched form (2026-11-13)

PI, on the first run: "we're using a GPU, let's actually use the GPU", then "do the port
properly".  What was left on the host was not one routine but everything organised as a loop
over chains or chain pairs.  cProfile of one host evaluation at sixteen chains and one thread
(3.45 s; `Ewald.terms` 2.61 s and the solve 0.30 s are the ported steps), the rest by function:

| what | seconds | how it was organised |
|---|---|---|
| the chains and their Jacobians (`MixedModel.build`) | 0.139 | one `repeat_chains_and_jacobian` per chain; 0.108 s of it the Jacobian of the five-block oligomer, of which one block is read |
| chain-pair sums (`_pair_and_dv` and the loop around it) | 0.172 | 144 blocks of `(5 images, 100, 100)` through the packer's kernel, one by one |
| same-site columns (`_intra_column_and_grad`) | 0.024 | one per chain |
| dipole exclusion corrections (`charge_dipole_exclusion`) | 0.092 | two dense `(5, 100, 100, 3)` evaluations per chain, of which a few hundred terms are not zero |
| charge exclusion correction (`exclusion_correction`) | 0.019 | one dense evaluation per chain |

All of it grows with the number of chains and none of it is the cube: on a card the cube takes
tens of milliseconds, and these took the rest.  The batched form
(`placed_energy_and_grad(..., xp=)`, `MixedModel._value_and_grad_batched`) is the same energy to
rounding as a fixed number of array operations whatever the number of chains.

**Pair sums: two stacks of blocks through one kernel, on the device.**  The blocks are the
loop's own (every chain pair and lateral image whose axes are within reach), chosen on the host
from the chains' centres with one further test that needs no pair distances (the nearest site
of one chain to the moved centre of the other): 144 of the 158 within reach at sixteen chains,
the same 144 the loop's exact test keeps.  Each site pair is evaluated at its nearest image
along the chain and as many further images either side as the cutoff reaches -- none once `c`
is at least twice the cutoff, which is every cell of this section -- where the loop evaluates a
range of images for the whole block (five at `c` = 25 A): a fifth of the elements.  The
same-site columns are one more stack, at half weight, the bonded scale of each element gathered
by its image (minus the image for a flipped chain).  The kernel is the packer's expression with
its mask, written against the array module, truncated Coulomb included; per-site gradients come
back through one product with a 0/1 chains-by-blocks matrix, the lattice gradient from the
block sums.  Sixteen chains are 1.4 million pair elements in the first stack, of which 0.26
million are inside the cutoff, and 0.16 million in the second.  A neighbour list would evaluate only
those; it would also need a list kept valid between evaluations, and a wrong list is a wrong
energy with no symptom.  The stack is rebuilt from the coordinates every time and cannot be
stale.

**Exclusion corrections: flat lists.**  The bonded pairs of a chain are a few hundred `(image,
i, j)` entries fixed by the topology (720 per chain for this repeat, among the electrostatically
active sites); `polyfind.ewald.exclusion_pairs` lists them for every chain of the cell, with
the image index in the placed frame, and `exclusion_correction_pairs` /
`charge_dipole_exclusion_pairs` evaluate each correction for all chains in one pass: 86 ms
chain by chain, 2.3 ms from the lists at sixteen chains, the same numbers to 5e-15.  They stay
on the host: 11,520 entries is below the size at which a device call is worth its launch.

**Valence terms**: `ChainValence.energy_and_grad_batch`, all chains that share a valence object
in one array, on the placed coordinates (a flipped chain's image indices change sign).

**The chains: one build, and torques in place of the Jacobian.**  The builder was already
batched over rows; batching alone gained nothing (117 ms chain by chain, 97 ms in one call),
because the time is the Jacobian of the coordinates -- forty variables times the 513 atoms of
five blocks, of which the evaluation contracts one block with a gradient.  That contraction does
not need the Jacobian: every variable is a rotation generator, the derivative of an energy along
it is the torque of the Cartesian gradient about its axis, `u . (T - a x F)`, and `F` and `T`
for all hundred bonds are suffix sums along the chain (`polyfind.chainjac.oligomer_vjp`).  The
alignment of the repeat to its reference is carried back the same way
(`aligned_block_vjp`: the chain rule of `aligned_block_jacobian` read backwards), and the repeat
transform's own derivative, which the closure terms need, comes from running sums over the two
backbone windows (`backbone_transform_jacobian`).  The chains of an evaluation, on this
container: 112 ms to 4.4 ms at sixteen chains and 432 ms to 13 ms at sixty-four; the gradient
agrees with the Jacobian's to 2e-14 of its size.  The sequential part -- the backbone
recursion, a hundred steps -- stays on the host, where a hundred small steps cost less than a
hundred launches.

**The closure constraints** rebuilt every chain's whole Jacobian to read two windows of it:
112 ms per call at sixteen chains and 432 ms at sixty-four, and a constrained relaxation calls it
as often as the energy.  `MixedModel.closure` now differentiates the window atoms alone, in one
build for all chains: bit for bit the same residuals and Jacobian (`tests/test_mixed_batched.py`
asserts it with `==`), 10 ms and 51 ms; beside the batched evaluation it uses the running sums,
4.6 ms and 14 ms, the same to 1e-10.  The first of these is a gain for every constrained
relaxation on the host, whatever form its energy takes.

**The solve is Cholesky on a device** (17.6): `polyfind.backend.cholesky_solve`.  CuPy does not
look at a factorization's status unless told to (`cupyx.errstate(linalg="raise")`; the default
is `ignore`, read in its source), so the status is asked for and a matrix that is not positive
definite goes to the general solve as before.

**Which form runs.**  `MixedModel(xp=None)`, the default, follows the packer's Ewald sums: on a
device when they are (`POLYFIND_EWALD_DEVICE=cuda` and a cell of at least 256 active sites),
otherwise the loops -- the form every stored number was computed with, unchanged and bit for
bit what it was (the closure included).  `POLYFIND_MIXED_BATCHED=1`, or `xp=numpy`, asks for the
batched form on the host; `xp="loops"` for the loops wherever the sums are.

**Agreement.**  Batched against loops, on the host: energy, every component of the gradient and
all nine lattice entries to 1e-15 of their size on twelve small cells (truncated Coulomb and
Ewald with induced dipoles, flipped chains, repeats short enough for seven images along the
chain) and on the sixteen-chain kinked cell, with both closure penalties and their multipliers
on.  `tests/test_mixed_batched.py` (50 tests) holds each piece to the routine it replaces, on
NumPy and on the device stand-in, where an evaluation after the first sends and brings back
nothing larger than a per-site array in a number of transfers that does not grow with the
cell; `tests/test_gpu.py` has the same comparison for a real device.

### 17.8 Measured on an A100: the whole evaluation in the batched form (2026-11-13)

A third pod (one A100-SXM4-80GB, CuPy 14.2.0, host Intel Xeon Platinum 8470 with 20 cores
allotted of 192 visible) ran commit 4349f40 for six minutes.  Files and the tables below
(`examples/ewald_device_summary.py`; nothing typed in): `pod3_*` in
`docs/predictions/ewald_device_a100_20261113/`.  The device tests passed (72, the batched form
against the loops on the card among them).

One whole evaluation, energy and gradient:

| chains | sites | k vectors | host loops, 20 threads | device | device against host | largest difference | card memory |
|---|---|---|---|---|---|---|---|
| 16 | 960 | 14,486 | 2.88 s | 48.8 ms | 59x | 8.6e-15 | 1.4 GB |
| 32 | 1,920 | 28,920 | - | 0.14 s | - | no host run | 5.5 GB |
| 64 | 3,840 | 57,886 | - | 0.76 s | - | no host run | 21.1 GB |
| 128 | 7,680 | 115,706 | - | 6.00 s | - | no host run | 82.5 GB |

Against 17.6 on the same card model: 511 ms to 49 ms at sixteen chains and 1.14 s to 0.14 s at
thirty-two.  Against the host times of 17.6 (another pod's cores, so to one figure): thirty-two
chains 10.8 s on thirteen cores, 75 times; sixty-four chains 59 s on eight threads, 78 times.

One device evaluation part by part (the device waited for at the end of each part; build,
exclusions and assembly are the host's):

| chains | build | pair sums | electrostatics | exclusions | assembly | the host's parts | of the evaluation |
|---|---|---|---|---|---|---|---|
| 16 | 6.0 ms | 4.4 ms | 34.7 ms | 0.9 ms | 2.3 ms | 9.2 ms | 19% |
| 32 | 9.4 ms | 5.9 ms | 0.12 s | 1.7 ms | 3.7 ms | 14.8 ms | 10% |
| 64 | 17.9 ms | 10.7 ms | 0.72 s | 3.6 ms | 6.3 ms | 27.8 ms | 4% |
| 128 | 34.9 ms | 23.3 ms | 5.90 s | 8.2 ms | 12.3 ms | 55.3 ms | 1% |

Evaluations back to back, each at a new geometry as a relaxation's are, with the card's own
utilisation figure (`nvidia-smi`) read five times a second meanwhile; "loops" is the
organisation of 17.6, the sums on the card and everything else in host loops:

| chains | form | one evaluation | per second | card utilisation, mean (lowest to highest sample) |
|---|---|---|---|---|
| 16 | batched | 49.5 ms | 20.2 | 49% (47 to 50%) |
| 16 | loops | 0.50 s | 2.0 | 5% (0 to 13%) |
| 32 | batched | 0.15 s | 6.8 | 77% (67 to 84%) |
| 32 | loops | 1.11 s | 0.9 | 10% (0 to 62%) |
| 64 | batched | 0.76 s | 1.3 | 91% (0 to 100%) |
| 128 | batched | 5.96 s | 0.2 | 92% (0 to 100%) |
| 128 | loops | 10.9 s | 0.1 | 56% (0 to 100%) |

Checks that need no host run: the untouched tiling against the stored two-chain cell 2.3e-13
kcal/mol at 16, 32 and 64 chains and 6.8e-13 at 128; the gradient against the device's own
energy difference 2.2e-5, 2.4e-6, 8.3e-6 and 2.7e-5 (the cutoff floor of 17.5).  The batched
form on the host at sixteen chains and twenty threads: 2.41 s against the loops' 2.88 s, the
energy to 1.3e-15 and the gradient to 6.0e-15.  The three steps apart: tensor and field 19.4 ms,
88.8 ms, 0.61 s and 5.40 s; the solve, now Cholesky, 3.4 ms, 9.4 ms, 45.6 ms and 0.28 s; the
gradient 8.9 ms, 18.9 ms, 58.7 ms and 0.21 s.  First evaluation of a process (kernels compiled
once per machine): 0.4 s to 0.9 s, and 8 s at 128 chains.

What it says:

* **The card now does the work.**  The host's share of an evaluation is a fifth at sixteen
  chains and a hundredth at a hundred and twenty-eight, and the card's meter reads 49%, 77%, 91%
  and 92% where the loop form read 5% and 10%.
* **Sixteen chains is a small job for this card.**  At 49 ms the evaluation is bounded by the
  issuing of work, not by arithmetic: the card is idle half the time between launches.  Two
  consequences.  A second process on the same card (another image of a band, another column of a
  Hessian) should cost little; that is not measured yet and is the first thing to measure
  before sizing a farm.  And further gains at this size are in the number of device calls the
  Ewald routine makes, not in what any of them computes.
* **From sixty-four chains the evaluation is the cube**: the tensor's reciprocal sum is 0.61 s of
  0.76 s, and 5.4 s of 6.0 s at a hundred and twenty-eight.  Inside it the six symmetric products
  are general matrix products, 88 ms each at sixty-four chains and 0.85 s at a hundred and
  twenty-eight, where the rank-k update takes 51 ms and 0.43 s (side measurement, same sizes).
  That reverses 17.6's verdict for the large cells: written as six rank-k updates (the mixed
  products through `k_a k_b = ((k_a + k_b)^2 - k_a^2 - k_b^2) / 2`, so that every weight is
  non-negative), the evaluation would lose about 0.2 s of 0.76 s and 2.5 s of 6.0 s.  Not done
  here: it can only be tested on a card.
* **Memory is the square, and a hundred and twenty-eight chains fill an 80 GB card**: 1.4, 5.5,
  21.1 and 82.5 GB of pool.  A 40 GB card holds sixty-four chains and not a hundred and
  twenty-eight; beyond that the tensor has to stop being stored, which is the same particle-mesh
  change that removes the cube (17.5).
* **What it buys the work of section 8 of the candidate document**: a sixteen-chain relaxation of
  five hundred evaluations is 25 s where it was 37 minutes on this container's four cores (4.4 s
  an evaluation); a nine-image band iteration half a second; the closure constraints beside each
  evaluation 5 ms.

**The farm question** (PI: four to six A100s for cheap), on these numbers.  One card already
changes what is feasible; more cards multiply only what is independent -- images of a band,
columns of a difference Hessian, separate starts, rows of a fit -- at one process per card, or
several per card at the sizes where one process leaves the card half idle.  A single
relaxation's evaluations follow one another and use one card.  Each process wants one fast host
core (the host's part is single-threaded: 9 to 55 ms an evaluation) and the host memory of its
cell's tables (the packer keeps about a dozen site-by-site tables: 1.3 GB each at a hundred and
twenty-eight chains).  Card memory decides the
largest cell: 80 GB cards for a hundred and twenty-eight chains, 40 GB for sixty-four.  So: two
80 GB cards cover everything in the present plan with bands at full speed; four to six pay off
only if bands or Hessians at sixty-four chains and more become routine, and before buying on
that expectation the two things to measure on a rented card are several processes on one card
and the rank-k update.

