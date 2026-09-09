# polyfind: algorithmic performance review

Scope: the algorithms, not the constant factors. The question for each stage is
whether the work being done is the *right* work, and what formulation would make
the same answer cheap. Numbers quoted for the current code are the ones from
the funnel run recorded in DESIGN.md section 5.3 (4-core cloud box): fit 7.6 s,
enumeration 3.6 s, packing of four conformations 81 s, refinement 88 s,
amorphous statistics 0.4 s. Packing and refinement are 95% of the wall time,
so the review starts there. Expected gains below are reasoned estimates, not
measurements.

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
(t, g, t, -g); beta: one deflection and two bond angles), so that the repeat
transform is a pure translation *by construction* and the penalty disappears.
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
* **Projected end-to-end effect.** With sections 1-4 implemented, the CPU
  pipeline is already dominated by the table builds and final polishes; the
  GPU then turns the table builds from tens of seconds into seconds and
  makes per-candidate refinement with autodiff interactive. All of this is
  projected: nothing has been run on the GPU beyond a device test.

## 10. Recommended order and expected cumulative effect

| Step | Change | pack() per conformation | refine per candidate | fit | full PVDF run |
|---|---|---|---|---|---|
| 0 | today | 20-60 s | 3-90 s | 7.6 s | ~3 min |
| 1 | section 4: parallel candidates | same | same | same | ~40 s |
| 2 | section 2: batched-gradient polish and refinement, hoisted invariants | 5-15 s | 1-5 s | same | ~15 s |
| 3 | section 6: intra constant, z-windows, tabulated pair potentials, float32 screen | 2-6 s | 0.5-2 s | same | ~8 s |
| 4 | section 1: tabulated W + FFT lattice sums, exhaustive screen | table 10-60 s once, then < 1 s | same | same | dominated by table builds |
| 5 | section 3: symmetry-parametrised refinement with bond angles | same | 0.2-1 s, better minima | same | |
| 6 | section 5: coarse-to-fine fit, batched conformers | | | 0.7 s; 3-4x fewer expensive calls | |
| 7 | section 9: torch backend on the AMD GPU | table builds in seconds; batched refinement | | batched MLIP calls | seconds |

Steps 1-3 are small changes to existing code; step 4 is the one real
reformulation and the largest algorithmic gain; step 7 is where the GPU
finally matters, and it should be built as a torch kernel rather than a CuPy
port.

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
the packing kernel, described in DESIGN.md section 5.5: every antiparallel configuration
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
