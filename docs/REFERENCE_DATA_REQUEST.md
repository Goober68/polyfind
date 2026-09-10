# Reference data request: bulk dipole response to strain

## Why this is the binding constraint

`polyfind` now reproduces both experimental piezoelectric **signs** for beta-PVDF
(d33 negative, d31 positive) where before d31 had the wrong sign, and it does so
out of sample, since the fit never saw a measured piezoelectric constant. The
**magnitudes** are 5x and 9x short and should not be trusted.

The limit is no longer the search, the speed, or the form of the potential. It is
the calibration data. The charge-flux coefficient that produces the response is
fitted against `sarco/materials/gpu_bundle/results/field_neighborhood_refined`,
which is exploratory GFN2-xTB on a **finite two-chain pair in vacuum with terminal
backbone atoms pinned**, at two strain values, with strain measured against a
fixed seed span rather than a stress-free lattice. Its own `interpretation` field
says as much.

Symptoms of that data being the limit, all measured (`docs/BENCHMARK.md`):

* the fit's R-squared is 0.39;
* a bond-stretch channel fits the same data better (0.86) and a plain charge
  rescale better still (0.75), and neither is a usable mechanism here;
* adding a bond channel alongside the angle channel changes the coefficient
  sixfold and **flips both signs back**;
* leave-one-chemistry-out spans the coefficient over a factor of eleven.

A bulk reference at the level the energy model is fitted to would settle all four.

## Provider status — 2026-09-10

The beta-PVDF PBE-D3(BJ) campaign is active in
`sarco/materials/gpu_bundle/periodic_reference/beta_pvdf`. The stress-free
periodic reference and all 13 unique relaxed-ion cells (zero plus +/-1% and
+/-2% on each of three axes) are complete. The accepted reference is
`a=8.3582753719`, `b=4.7314155341`, `c=2.5801545869` A with residual pressure
-0.554 MPa. Sarco axes are `x=a`, `y=b` (polar) and `z=c` (chain); Polyfind's
canonical packing frame maps its polar axis to `x` and chain axis to `z`.

The preliminary relaxed-ion normal stiffness diagonal is
`(111.934, 25.185, 315.851)` GPa in Sarco x/y/z order. The full raw normal
stiffness matrix is
`[[111.934, 8.547, 1.335], [10.054, 25.185, 2.817],
[1.249, 2.324, 315.851]]` GPa. Direct-fit R-squared is at least 0.999266;
the maximum reciprocity mismatch is 1.507 GPa and remains a reported numerical
diagnostic rather than being hidden by symmetrization.

The first 4x8x4 Berry sweep is rejected and cannot be used for the charge-flux
fit. Exact-checkpoint repeats isolated nondeterminism in parallel Quantum
ESPRESSO 7.6 Berry evaluation; restricting it to one MPI rank was insufficient
while multiple OpenMP threads remained. Two repeats on the full unsymmetrized
k-point grid at one MPI rank and one OpenMP thread both returned
`P_y=-0.4086832 C/m2` for the x/+1% cell. The runner now enforces that
execution shape, requires every Berry stage to load the saved SCF wavefunctions,
and keeps independent checkpoint copies for x/y/z. A clean 13-cell sweep is the
next running gate; polarization derivatives and the Polyfind refit remain
deliberately unpublished until branch unwrapping, linearity and curvature checks
pass. The zero-strain Born-effective-charge/phonon workflow is prepared behind
that gate, followed by polyethylene and alpha/gamma-PVDF controls.

## What is needed

**A periodic-crystal polarization response to strain.** Concretely, for
beta-PVDF:

1. **Periodic bulk cell**, not a finite cluster. The beta cell is about twelve
   atoms, so this is small.
2. **Polarization by Berry phase** (or Wannier centres). This is not optional: in
   a periodic system the dipole of a cell is not well defined, and a naive sum of
   charge times position depends on the choice of origin and cell. Take
   differences along a continuous strain path and watch for jumps by the
   polarization quantum `eR/V`; a jump misread as a response is the classic way
   this calculation goes wrong.
3. **Several strain points per component, not two.** Suggested
   `-0.02, -0.01, 0, +0.01, +0.02` for each of `eps_zz` (along the chain),
   `eps_xx` and `eps_yy`. Five points gives a linearity check and a curvature
   estimate; two gives a secant and no way to know it is one.
4. **Ions relaxed at each fixed strain**, with the cell held at the strained
   shape. The relaxed-ion response is the physical piezoelectric one. If the
   clamped-ion values come free, they are useful too, because the difference
   between them is exactly the internal-strain contribution that the charge-flux
   model is trying to represent.
5. **Report, per strain point**: the strain tensor, the Berry-phase polarization
   vector, the total energy, the relaxed atomic positions, and the residual
   forces and stress.
6. **Level of theory: PBE-D3**, to match the energy training set
   (`trainset_v1_566.xyz`). Matching matters more than absolute quality here: a
   flux fitted against one functional and an energy model fitted against another
   will not compose. If a hybrid is affordable, one point at PBE0-D3 as a spot
   check on the functional dependence would be valuable.

That is fifteen small periodic calculations, plus relaxations. On the hardware
described it should be hours, not days.

## Run it on CPU, not GPU

The cell is about twelve atoms, which is far below where plane-wave
density-functional theory starts to benefit from a GPU. At this size the
transforms are too small to saturate a device and the calculation is
latency-bound, so a GPU mostly adds transfer overhead; the crossover is in the
hundreds of atoms. GPU support in the common plane-wave codes is also
CUDA-oriented, with ROCm less mature.

The scheduling consequence is the useful part: **this work does not compete with
GPU-bound simulation work**. It wants CPU cores and can run alongside rather than
queueing behind whatever is occupying the accelerator.

## Worth more than the above, if it is cheap to add

**Born effective charge tensors** `Z*_ij` per atom, from density-functional
perturbation theory at zero strain. These are the exact physical content of
charge flux: the difference between `Z*` and the static charge is precisely what
a fixed-charge model is missing, so they calibrate the mechanism directly instead
of through a fitted response. Most periodic codes produce them alongside the
piezoelectric tensor from the same run.

If the piezoelectric tensor itself comes out of that run, it is the direct answer
and the fit becomes a validation rather than a calibration.

## Controls that make the result trustworthy

* **Polyethylene, same protocol.** Its response must be exactly zero by symmetry.
  If it is not, something is wrong with the polarization branch tracking rather
  than with the polymer.
* **alpha-PVDF**, whose antipolar structure should give a much smaller response
  than beta, and **gamma**, which is polar but weaker. Getting the ordering right
  is a cheap sanity check on the whole protocol.
* Confirm the reference structure is **stress-free** before straining it. The
  existing data's strain is defined against a pinned span, which is part of why
  it is hard to use.

## Cross-check that costs nothing

Published work already computes polarization and Born effective charges for
beta-PVDF, including a comparison of the planar and alternately-deflected
structures (`docs/REFERENCES.md`, section 5.8's sources). Comparing any new
calculation against those before trusting it is worth the ten minutes.

## What we will do with it

Refit the charge-flux coefficients against `dPolarization/dStrain` for a bulk
crystal at the same level as the energy model, then recompute d33 and d31. The
test is already written and is out of sample: the fit does not see a measured
piezoelectric constant, so agreement with the experimental values, in sign and
then in magnitude, is a genuine prediction. If the magnitudes remain short after
this, the deficiency is in the model rather than the data, which is a different
and equally useful thing to learn.

## Cheaper fallbacks, in order

1. **Same GFN2 protocol but periodic and stress-free**, with five strain points.
   Keeps the level, removes the finite-size and pinning problems, and would show
   how much of the current scatter those cause.
2. **PBE-D3 at three strain points on beta only.** Enough to fit one coefficient
   against a bulk reference, which is the single most valuable number.
3. **Born effective charges at zero strain only.** No strain sweep, but it
   calibrates the mechanism directly and is one calculation.

Any of these is a large improvement on what the coefficient currently rests on.

## Consumer response — elastic cross-check, 2026-09-10

Thank you, and the Berry nondeterminism catch is the important part of that
status: a silent parallelism-dependent result in the one quantity this whole
request exists to obtain is exactly the failure that would have propagated into
our fit unnoticed. Enforcing one rank and one thread, and requiring each Berry
stage to reload the saved wavefunctions, is the right response.

Your normal stiffnesses let us cross-check immediately, and the result is two
agreements and one outlier worth your attention.

**Axis mapping used.** Please confirm it, because reversing it inverts the
conclusion. Your cell is `a=8.3583, b=4.7314 (polar), c=2.5802`. Ours is
`a=4.542, b=8.546, c=2.600`, with the polar axis first. So we read your **b**
against our **a**, your **a** against our **b**, and chain against chain.

| direction | yours (GPa) | ours, illustrative | ours, fitted | published DFT |
|---|---|---|---|---|
| chain (your c, our c) | 315.9 | - | **328.4** | 287-341 |
| polar transverse (your b, our a) | 25.2 | 45.3 | **21.9** | about 20-25 |
| other transverse (your a, our b) | **111.9** | 38.4 | **17.1** | about 20-25 |

Cells agree well: our polar axis 4.542 against your 4.731 (4% low), our long axis
8.546 against your 8.358 (2% high), chain 2.600 against your 2.580 (0.8% high).

**The chain axis agreeing to 3.9% is the headline.** Two independent routes, a
fitted classical potential here and PBE-D3 there, on the quantity we had only
been able to compare against literature. That is mutual validation.

**The polar transverse agrees to 13%**, and both land in the published range.

**Your a-direction constant of 111.9 GPa is the outlier, and we would query it
before using it.** It is about five times our fitted value and about five times
the published DFT transverse constants for this phase, while your b-direction
constant matches both. A van der Waals-bonded direction in a polymer crystal
being stiffer than 100 GPa is surprising on its face.

The hypothesis we would test first is **under-relaxation of chain
reorientation**. Transverse compliance in these crystals comes mostly from the
chains rotating about their own axes and sliding, rather than from compressing
contacts. If the relaxation at fixed strain did not fully release the setting
angles and the relative axial shift, the result approaches a clamped constant,
which is stiff. Two cheap diagnostics: compare the relaxed-ion value against the
clamped-ion one for that axis specifically, since a large gap there and a small
gap on b would confirm it; and check how far the setting angles actually moved
between the zero-strain and plus-two-percent a cells.

It may equally be real and our transverse constants may be the wrong ones. But
the asymmetry between your two transverse directions is the thing we would want
explained before either of us builds on it.

**On our side**, our own transverse constants carry a caveat in the other
direction: the fitted potential roughly halves them against the illustrative one
and moves them toward the published range, which is consistent with the
illustrative potential being too stiff rather than with our fit being validated.

**No change to the request.** A clean 13-cell Berry sweep remains what we need,
and holding the polarization derivatives back until branch unwrapping, linearity
and curvature checks pass is the right call. We would rather wait than refit
against a rejected sweep.
