# polyfind: algorithmic performance review

Scope: the algorithms, not the constant factors. The question for each stage is
whether the work being done is the *right* work, and what formulation would make
the same answer cheap. The numbers the review was *written* against are the
funnel run recorded in DESIGN.md section 5.3 (4-core cloud box): fit 7.6 s,
enumeration 3.6 s, packing of four conformations 81 s, refinement 88 s,
amorphous statistics 0.4 s; the expected gains in sections 1-9 are the reasoned
estimates made at that time, not measurements. **Every "today" number in
sections 1-9 is stale: section 10 now carries a measured stage-by-stage profile
of the current code and is the only profile to plan from.** What it says, in one
line: the table build is no longer the bottleneck and torsion refinement is
69% of the funnel.

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
  GPU would then matter most for the table builds. The first half is now wrong
  and the second follows it: with the build split across processes it is 10% of
  the funnel and refinement is 69% (section 10), so the gradient batches of
  section 2, not the table build, are what a GPU should take. Still projected:
  nothing has been run on the GPU beyond a device test.

## 10. Where the time actually goes (measured), and what to do next

All numbers below were measured on the target machine -- Windows 11, Python
3.12, **Intel i7-8700K, 6 physical cores / 12 threads**, `OMP_NUM_THREADS=1`,
with another agent's work running alongside -- as the minimum of repeats, and
each "before" was taken back to back with its "after". The "before" column is
commit `bb2e4de`. The workload is the default `PipelineConfig` for PVDF with
`SimpleFF`: a fit, enumeration to period 8, and five conformations packed and
refined (TG+TG+TG-TG-, TTTG+TG-, TTG+TG-G-TG+, TG+TG-, TT).

### The funnel, stage by stage

One process, no candidate parallelism, cold table cache -- so the columns are
CPU work, not a wall time with contention in it:

| Stage | before | after | share of the funnel now |
|---|---|---|---|
| RIS fit (`fit_ris`, step 10 deg, 6 monomers, third order) | 0.31 s | 0.30 s | 0.1% |
| enumeration (period <= 8, k = 60) | 4.9 s | 4.4 s | 1.9% |
| periodic chain construction (5) | 0.12 s | 0.14 s | 0.1% |
| **chain-pair table builds (5, screen grid)** | **66.1 s** | **24.8 s** | **10%** |
| FFT screen (5 conformations, ~10^8 landscape points each) | 8.1 s | 7.9 s | 3% |
| exact-kernel polish (4 starts per conformation) | 33.8 s | 34.7 s | 15% |
| **torsion + cell refinement (5)** | 155.7 s | **165.0 s** | **69%** |
| amorphous + interface statistics | 0.38 s | 0.38 s | 0.2% |
| total | 269.4 s | 237.5 s | |

Refinement is now the stage to attack, and it is lopsided: of its 165 s, 112 s
is one candidate (TTG+TG-G-TG+, an 8-bond 24-atom repeat), i.e. 47% of the whole
funnel sits in one L-BFGS run. The polish is the second target at 15%.

### End to end

| Run | before | after |
|---|---|---|
| `run_pipeline(pvdf)`, default workers (5), cold cache | 289.1 s | 222.7 s (**1.30x**) |
| same, `workers=3` (one per physical core, minus the main process) | - | 179.0 s |
| same, `workers=3`, warm on-disk table cache | - | 146.3 s |

The refined cell is unchanged: beta-PVDF TT wins at -8.153 kcal/mol per monomer,
a = 4.60, b = 8.59, c = 2.63 A, antipolar, in all four runs.

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

1. **Refinement (69%).** Its L-BFGS iterations are sequential, but each iteration
   is one batched kernel call of `1 + 2 n_vars` configurations with per-row
   coordinates -- embarrassingly parallel across rows, exactly like the radial
   axis of the table. Before splitting it, repeat the chunk-size experiment on
   `pack.CPU_CHUNK_ELEMS` (60,000, i.e. ~4 MB of temporaries): the table build
   gained 1.3x serially and 4x in parallel from that one constant, and the
   packing kernel has the same shape. Cutting the iteration count is the other
   half: 112 s in one candidate is an optimiser-conditioning problem.
2. **Polish (15%).** Four independent L-BFGS runs per conformation, each a
   sequence of small batched calls: parallel across starts, with the same
   caveat about the kernel's working set.
3. **Pipeline worker count.** `cpu_count() - 1` is 11 on a 6-core machine and
   oversubscribes it; 3 workers beat 5 by 1.24x (section 4).
4. **GPU (section 9).** The table build is no longer where a GPU would pay;
   the refinement gradient batches are.

### Status of the earlier plan

| Step | Change | state |
|---|---|---|
| 1 | section 4: parallel candidates | done; 1.3x, not 4-8x, and oversubscribed at the default worker count |
| 2 | section 2: batched-gradient polish and refinement, hoisted invariants | done |
| 3 | section 6: intra constant, z-windows, float32 screen | done; tabulated pair potentials not done |
| 4 | section 1: tabulated W + FFT lattice sums, exhaustive screen | done; the screen is 3% of the funnel |
| 5 | section 3: symmetry-parametrised refinement with bond angles | done; still 69% of the funnel |
| 6 | section 5: coarse-to-fine fit, batched conformers | done; the fit is 0.1% of the funnel |
| 7 | section 9: torch backend on the AMD GPU | not done |
| 8 | process-parallel table build (this section) | done; 2.8-3.8x on the build, 1.30x end to end |

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

## 13. The bandwidth lesson, and a lead that did not transfer

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

**Where the time now is.** The profile has moved enough that the old
recommendations are stale:

| Stage | Share of pipeline |
|---|---|
| refinement | 69% |
| polishing | 15% |
| table build | 10% |
| FFT screen | 3% |
| fit, enumeration, amorphous | about 2% together |

Refinement is the target, and 112 of its 165 seconds is a single candidate, the
gamma chain, whose per-row kernel cost is 18 ms against 2.3 for alpha and 0.9
for polyethylene - an N-squared effect from its 24-atom repeat that no constant
will fix.

The lever there is not the kernel but the *number of rows*. Refinement gets its
gradients by central differences, so every L-BFGS iteration costs one plus twice
the variable count, 29 rows for the gamma chain, to extract a single gradient.
Analytic or automatically-differentiated gradients would collapse that to one
evaluation per iteration at some increased cost per evaluation - plausibly five
to ten times overall, against roughly 30 times fewer rows. Section 2 of this
document already pointed at automatic differentiation as the route, for the
separate reason that the same implementation would run on a GPU; the profile now
says it is also where the remaining CPU time is. That is the next thing to do.

A smaller fix is already applied: the pipeline's worker count defaulted to
logical cores minus one, which oversubscribes a bandwidth-bound workload. On
this machine three workers beat five by 1.24x, so the default is now derived
from physical cores.
