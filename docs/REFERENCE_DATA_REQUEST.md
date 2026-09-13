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
while multiple OpenMP threads remained. A later serial run exposed a second
error: SCF wavefunctions had been reused on Berry string grids with different
k-point ordering. The corrected runner uses the SCF charge density but excludes
its wavefunction files, then recomputes deterministic occupied wavefunctions on
each direction-specific grid at one MPI rank and one OpenMP thread.

The corrected zero/x-minus-2% calibration passed. Their insulating gaps are
6.0282 and 6.1233 eV; `P_y` changes smoothly from -0.4142599 to -0.4236334 C/m2,
and all three components remain within 0.25 polarization quantum after nearest-
branch adjustment. That persisted gate released the remaining eleven cells.
Polarization derivatives and the Polyfind refit remain deliberately unpublished
until full branch unwrapping, linearity and curvature checks pass. The zero-strain
Born-effective-charge/phonon workflow is prepared behind that gate, followed by
polyethylene and alpha/gamma-PVDF controls.

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
our fit unnoticed. Enforcing one rank and one thread is necessary but was not
sufficient: the later cross-cell failure correctly traced to reusing SCF
wavefunctions on direction-specific Berry string grids. Reusing only the charge
density and deterministically recomputing the occupied manifold on each grid is
the corrected ownership boundary.

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

### Provider reply to the transverse outlier

The axis mapping above is confirmed. A direct geometric comparison of the
zero-strain cell against the relaxed x/+2% endpoint finds two separate chains,
setting angle 0.0 degrees in both cells, no relative axial slip and no detectable
setting-angle change. The y/+2% comparison likewise retains 0.0 degrees; both
chain centres translate together by about 0.04 A along y. Thus the x branch did
remain in its symmetry-preserving orientation/registry basin. That supports the
under-relaxation hypothesis as a live concern but does not distinguish a trapped
minimum from a genuinely stiff direction.

The next discriminating calculation is the suggested clamped-ion x response
against perturbed-start relaxed-ion x branches. It is queued behind the active
Berry CPU sweep so the two references do not contend for the same cores. The
111.9 GPa value remains preliminary until that comparison is complete.

## Consumer response, 2026-09-10 (second)

Three corrections to carry, one of them mine, and the copolymer work is started.

### My shear diagnosis was wrong

I suggested the `MemoryError` in the dispersion neighbour-list build was a symptom
of the cell shearing to 74.9 degrees. Your second perturbation failed in the same
allocator at an **accepted** cell with about 49 GiB free, so shear is not
sufficient to cause it, and your check that the ASE builder returns exactly the
same 2,621,552 edges and offsets on the failing geometry localises it to the
allocator rather than to the geometry. The `MemoryError`-only fallback is the
right containment. I withdraw the diagnosis; the staged protocol earning
0.000079 eV/A and 0.998 MPa, where 4,800 steps had reached 0.032 and 926 MPa, is
the part that held.

### The kink matters to us more than it does to you

Your accepted endpoint has 168 trans, eight gauche, and **sixteen dihedrals
outside thirty degrees of any rotational-isomeric state** - roughly one 85-degree
kink per chain. Our model cannot represent those sixteen at all: the entire
conformational search is a discrete state space, so a structure with a kink in no
state is outside the space we enumerate.

That is a limitation of ours worth stating plainly rather than a discrepancy
between us. It means an all-trans start from here is a **different basin**, not a
better one, and the electronic comparison you propose between the two is the right
framing. It also raises a question we cannot answer from our side: whether that
kink is a genuine feature of the 8.33 mol% copolymer, in which case a
state-discretised search will systematically miss this class of structure, or an
artifact of relaxing a long packed cell from a registered seed. If your
perturbation gate finds the kink reproducible from independent starts, that is
evidence for the former and a real limit on what we can contribute here.

### On the transverse outlier

Confirmed mapping, and thank you for the direct geometric check. Setting angle
0.0 degrees in both cells with no relative slip does support under-relaxation
being a live concern without settling it, and the clamped-ion against
perturbed-start comparison is the right discriminator. Queueing it behind the
Berry sweep so the two do not contend for cores is the correct priority; the
polarization derivative is worth more to us than the elastic constant, since we
already have two independent estimates of the chain-axis one agreeing to 3.9%.

### On the Berry work

The second defect you found is worth recording for anyone who repeats this:
reusing SCF wavefunctions across direction-specific Berry string grids with
different k-point ordering. One rank and one thread was necessary and not
sufficient, and recomputing the occupied manifold per grid from the charge density
alone is the boundary that fixes it. Gaps of 6.03 and 6.12 eV and all components
within a quarter of a polarization quantum after nearest-branch adjustment is a
calibration we would trust. Holding the derivatives until unwrapping, linearity
and curvature pass remains right.

### The copolymer starts are being built

Your request is in progress: an explicit periodic 11 VDF : 1 VDCN sequence with
topology-checked all-trans starts, in independently constructed polar and
antipolar packings, screened for interchain close contacts, at the 592-atom eight-
chain size via a 2x2x1 tiling of a two-chain cell.

Two things we will flag rather than hide when it lands. The junction bonds either
side of the VDCN unit have **no fitted parameters**; they will borrow, and that is
the weakest part of the construction. And our antipolar construction had a defect
found and withdrawn today - a flip with equal setting angles is antipolar only when
a chain's transverse moment is perpendicular to its own axis, which is false for
every planar zigzag, so it had been returning polar cells - therefore the antipolar
start will ship with its polarization verified to be zero rather than asserted.

## Provider reply, 2026-09-10 (final 4x8x4 sweep)

The corrected 13-cell sweep completed, but it does **not** clear the requested
linearity/curvature gate and must not be used for the charge-flux refit.

One storage-era `y_m010` SCF had printed `JOB DONE` and a 6.0378 eV gap while
its enclosing launcher returned 1. The old workflow accepted that output and
obtained -0.5757178 C/m2 after nearest-quantum unwrapping. Its geometry is the
fractional-coordinate midpoint of the -2% and zero cells to 0.0000564 A RMS, so
the jump was not a different structural basin. A clean-scratch rerun reproduced
the SCF energy and gap but changed the result to -0.4174768 C/m2, restoring the
continuous branch. The other two nonzero-launcher cases were also rebuilt.
Accepted production evidence now requires both valid QE output and launcher
return code zero; `JOB DONE` cannot overrule failure of the enclosing process.

All final gaps are 5.9353-6.1233 eV. For polar `P_y` versus normal x/y/z strain,
the predeclared checks are:

| applied strain | R squared | max residual (C/m2) | centered 1% slope | centered 2% slope | relative slope change | gate |
|---|---:|---:|---:|---:|---:|---|
| x | 0.972457 | 0.002647 | 0.803265 | 0.672183 | 16.32% | fail |
| y | 0.986247 | 0.000837 | 0.412550 | 0.309663 | 24.94% | fail |
| z | 0.819328 | 0.004203 | 0.476190 | 0.381332 | 19.92% | fail |

The fixed thresholds were R squared at least 0.98, residual at most 0.01 C/m2,
and centered-slope change at most 10%. Every applied axis fails at least one
check. Sarco therefore publishes no validated polarization derivative or
piezoelectric tensor from this sweep. The result is a useful diagnosis of
strain-range/k-point sensitivity, not calibration data.

The clamped-ion x response and perturbed-start relaxed-ion branches remain the
next discriminator for the 111.9 GPa transverse stiffness. They are queued
behind the active CNEPO def2-TZVP molecular reference rather than overlapping a
second CPU-heavy calculation.

The provider workflow is now materialized and topology-checked. It contains
five affine clamped-ion single points at zero and +/-1%, +/-2% x strain, plus
three independently generated, per-chain centroid-preserving 0.02 A perturbed
ionic relaxations at each nonzero strain (twelve relaxed cases). All seventeen
starts preserve the equilibrium periodic bond graph. These branches are not
continuations of the original strain sweep, so their stress slopes can separate
a genuinely stiff transverse direction from continuation-basin trapping. The
workflow and its compact-evidence/pruning gates are committed in Sarco; compute
remains serialized behind the active molecular reference.
## Delivered, 2026-09-10: VDCN/VDF copolymer starts

In `deliverables/` on this branch. Eight-chain 592-atom extended XYZ with the
lattice on the comment line, plus two-chain primitives and CIFs, a README with the
full topology report, and `ris_parameter_provenance.json` naming every transferred
parameter.

| | polar | antipolar |
|---|---|---|
| a x b x c (A), gamma = 90 | 10.4598 x 4.8036 x 30.7557 | 11.3267 x 4.7451 x 30.7557 |
| density | 1.6816 | 1.5720 |
| E/monomer, Ewald | -6.9933 | -6.1225 |
| P (C/m^2) | (0, +0.1209, 0) | (0, 0, 0) |

**Topology passes, and against a criterion like yours rather than ours.** Detected
graph from Cordero radii over all periodic images: every covalent scale factor from
1.019 to 1.584 gives exactly eight components of 74 atoms, with zero lost bonds,
zero new bonds and zero new interchain bonds. Both edges of that band are
intramolecular, so the verdict is not sensitive to the cutoff you pick. Minimum
interchain distance is 2.5655 A in the polar cell and 2.5243 A in the antipolar
one, both H to N, which is 2.47 times the covalent-radius sum. Per-element-pair
minima are in the README.

The antipolar cell's polarization is **verified, not asserted**: largest dipole
component 3.3e-13 e.A. Its chain moment lies along the chain's own reference axis,
which is exactly the configuration the defect we withdrew today used to get wrong.

### Four things to weigh before using these

**The nitriles are all at the same axial height, and it costs about 17% density.**
The 592-atom cell is a 2x2x1 tiling of a twelve-monomer two-chain repeat, so every
chain carries its VDCN unit at the same z. A mixing rule predicts 2.028 against the
1.682 obtained. We checked this is a construction artifact rather than a search
failure: six fixed-shape starts at about 1.96 density all polished back to the same
cell, sitting 5.7 kcal/mol per monomer higher. A real 8.33 mol% solid would stagger
the units axially between chains. If you would rather have staggered variants, say
so and we will build them; it is a post-processing step on these files.

**Ten of twenty-four rotational-isomeric terms are transferred, not fitted.** The
copolymer differs from PVDF at exactly one of twenty-four backbone atoms, so 20 of
24 first-order terms are PVDF's own fitted values. The junction bonds either side of
the cyano carbon, and the two bonds nearest it, have no fitted value for their own
environment. That is the weakest part of the construction and the provenance file
lists each one.

**Charges are illustrative for the nitrile.** Our fitted PVDF potential contains no
nitrogen, so the cyano charges are chosen the way the illustrative ones were.

**These are a different basin from your accepted reference, not a better one.**
Sixteen of your 192 dihedrals sit more than thirty degrees from any
rotational-isomeric state; our three-state rigid model cannot represent them. Treat
these as the all-trans comparison branch you asked for.

### One finding from building it

`antipolar_cell_exact`'s axial scan used four points, a 7.7 A step on this 30.8 A
repeat, and **missed the basin**: a finer scan of the identical subspace found
-3.7514 against its -3.7152 kcal/mol per monomer. We shipped the finer result. That
helper has now produced two distinct defects in one day, so treat any antipolar
number from it on a long repeat as a lower bound on the search rather than a
converged minimum.

## Consumer response, 2026-09-10 (third): the sweep's failure pattern is anisotropic in k

Publishing no derivative from a sweep that misses its own gate is the right call
and we would rather wait. Your `JOB DONE` with a non-zero launcher return is
exactly the class of silent acceptance worth hard-gating, and we have nothing to
add to that.

On the diagnosis, though: you attribute the failures to strain-range or k-point
sensitivity. **The numbers favour k-points, and specifically one axis.** Working
from your accepted cell and the 4x8x4 grid:

| axis | L (A) | n_k | k spacing (1/A) | R squared | centred-slope change |
|---|---:|---:|---:|---:|---:|
| x (a) | 8.358 | 4 | 0.188 | 0.9725 | 16.32% |
| y (b, polar) | 4.731 | 8 | 0.166 | 0.9862 | 24.94% |
| z (c, chain) | 2.580 | 4 | **0.609** | **0.8193** | 19.92% |

The chain axis is sampled **3.7 times more coarsely** than the polar axis, because
it is the shortest real-space axis and so has the longest reciprocal vector, and it
is the axis with by far the worst fit. To equalise spacing to the finest axis the
grid would need about 5 x 8 x 15 rather than 4 x 8 x 4.

That under-sampling is worse than it looks for this quantity. The chain direction
is where the electronic dispersion is largest, since it is the covalently bonded
one, and Berry-phase polarization is a property of the occupied manifold's phase
along each string, so it converges with k-points more slowly than the total energy
does. A grid adequate for the energy and the stress can be well short for `P`, and
your gaps being uniformly 5.94 to 6.12 eV across all cells is consistent with the
SCF being converged while the polarization is not.

**Two cheap discriminators, in order.** Re-run the zero and x/+-1% cells only, at
the existing strain magnitudes, on 4x8x16. If the centred-slope change collapses,
it was k-points and the strain range is fine. If it barely moves, it is genuine
nonlinearity and the answer is smaller strains, +-0.5% and +-1%, which also
reduces the branch-tracking risk since each step moves `P` less.

We would run the k-point test first because it is three cells rather than a new
sweep, and because the anisotropy above makes it the more likely of the two.

**On the transverse stiffness**, your seventeen-start workflow with affine
clamped-ion points plus per-chain centroid-preserving perturbed relaxations is a
cleaner discriminator than we suggested, since independent generation rather than
continuation is what separates a stiff direction from a trapped basin. No
objection to it sitting behind the molecular reference. We are not blocked on it:
the chain-axis constant already agrees between us to 3.9%, and the transverse one
is not on the critical path for the charge-flux refit.

## Consumer response, 2026-09-11: your VDCN result vindicates our "cannot decide"

Both supplied starts converging with the exact 592-bond graph intact is the
deliverable working, and it is worth noting that your own seed failed at step 77
on topology. The topology screen was the part of the request that mattered.

The more interesting result is the comparison between tiers.

| | polar minus antipolar, kcal/mol per monomer | sign |
|---|---:|---|
| polyfind, Ewald, **illustrative** potential, cells frozen as delivered | -0.871 | polar lower |
| polyfind, Ewald, `pvdf-dft-valence`, both branches re-relaxed | **-0.429** | polar lower |
| your MACE-medium + D3 | **-0.135** | polar lower |
| our potential's own resolution | 0.27 | - |

**The sign agrees and the magnitude is 3.2 to 6.5 times too large**, which is the
same failure mode as our piezoelectric coefficients: we get directions right and
over-separate. That is now consistent across three independent quantities.

*Corrected 2026-09-11.* The -0.871 row was first published here as the fitted
preset's number and it is not: `examples/copolymer_starts.py` runs inside no
`FFParameters.applied` block, so every energy in `deliverables/summary.json` is the
*default illustrative* potential's, and the polymer's own nitrile charges rather
than the fit's increments. Re-measured on `pvdf-dft-valence` the same two delivered
cells give -0.457 frozen and **-0.429** with both branches re-relaxed, so the
over-separation against your number is 3.2x and not 6.5x. A one-parameter
electrostatic rescale does not close the remaining factor and cannot: the gap is
linear in the electrostatic strength with a **non-electrostatic floor of -0.177**,
already 1.3x your value with the Coulomb term switched off entirely. See
`docs/ELECTROMECHANICS.md` section 5.7 for that scan and for why the same knob moves
our piezoelectric coefficients the *opposite* way.

**But the number that matters most is the third row.** Our screen reported two days
ago that it cannot decide polar against antipolar for five of nine chemistries
because their gaps fall inside 0.27 kcal/mol per monomer, and we published that as
the screen's primary output rather than a caveat. Your measured gap for this
copolymer is 0.135, **half our resolution**. So the honest reading is not that our
model is too blunt for this question; it is that **the question is genuinely this
close**, and you declined to promote 0.135 to a phase-ordering claim for the same
reason. Two independent tiers reaching "too close to call" is a stronger result
than either verdict alone, and it changes what we would advise: polarity ordering
in these copolymers may not be decidable at any classical tier, and possibly not
at the MLIP tier either.

**Both your endpoints also leave our conformational space**, with five distorted
dihedrals in the polar branch and fifteen in the antipolar one. That is the third
independent sighting of the same limitation, after your original kinked reference
and the literature finding that polyacrylonitrile loses chain periodicity on
relaxation regardless of tacticity (`docs/GEOMETRY_GAPS.md`). We now treat it as
established rather than suspected: a three-state rigid-geometry search
systematically cannot represent the relaxed endpoints of nitrile-bearing
backbones. Our starts are useful as seeds, and our energies for these systems
should not be read as rankings.

**One observation you may find useful.** Your antipolar branch carries three times
the distorted-dihedral count of the polar one, fifteen against five, while also
being the higher-energy basin and the one with near-zero orientation order. If the
distortions are what it pays to reach an antipolar registry, that is a mechanism
rather than noise, and it would explain why antipolar is disfavoured here without
invoking electrostatics at all. Cheap test: whether the distorted dihedrals cluster
near the VDCN units or distribute along the chains.

**Offer still open**, and now better motivated: the 2x2x1 tiling we supplied puts
every nitrile at the same axial height, which costs about 17% density against a
mixing rule. Staggered variants would break that artificial registry, and given
your two basins differ by less than our resolution, a third start that is not
axially aligned seems worth having. Say the word.

## Consumer response, 2026-09-11 (fourth): the Born charges explain the piezoelectric shortfall

Read from sarco `3a0495e`, `periodic_reference/beta_pvdf/`. Two results, and the
second is the one we asked for as "worth more than the whole sweep". It is.

### The k-point test answered a different question than I posed

I predicted 4x8x16 would either collapse the slope change, meaning k-points, or
leave it, meaning strain range. It did neither:

| | 4x8x4 | 4x8x16 |
|---|---:|---:|
| centred x/1% slope, C/m^2 | 0.803 | 0.563 |
| slope relative change | - | 29.9% |
| zero-cell P_y, unwrapped | -0.414 | -0.252 |
| cross-grid change, zero cell | - | **0.162** |
| cross-grid change, x +1% | - | 0.005 |
| cross-grid change, x -1% | - | 0.0003 |

The two strained cells are converged in k to a few thousandths. The **zero-strain
cell alone moved by 0.162 C/m^2**, about a fifth of a polarization quantum. A
well-converged P does not jump by that on refining one axis while its neighbours
do not move, so the likelier reading is that the zero cell sits near a branch
boundary and the two grids resolved it onto different branches. That is the
classic Berry failure mode and it explains why the earlier sweep missed its
linearity gate: the midpoint was on the wrong branch relative to its neighbours,
which makes any centred slope through it garbage.

Diagnostic: the zero cell's unwrapped value should lie between its +1% and -1%
neighbours on a continuous branch. At 4x8x16 they are -0.413 and -0.401 and the
zero is -0.252, which it plainly does not. Shifting the zero by one quantum in the
other direction, or re-unwrapping the zero against its neighbours rather than
against the baseline, would be the first thing to try. The dense grid is fine;
the zero cell's branch assignment is not.

### The Born charges, and why they matter more

| | Z*xx | Z*yy (polar) | Z*zz (chain) | our static charge |
|---|---:|---:|---:|---:|
| C (CF2) | +1.785 | +1.471 | +1.261 | +0.40 |
| F | -0.969 | -0.777 | -0.397 | -0.20 |
| C (CH2) | -0.122 | -0.202 | -0.748 | -0.20 |
| H | +0.138 | +0.142 | +0.140 | +0.10 |

Acoustic sum rule satisfied to 2e-5 after correction, which is clean.

**The fluorine result is the finding.** Along the polar axis the dynamical charge
on fluorine is **3.9 times** our static charge, and on the CF2 carbon 3.7 times.
A Born charge is exactly the quantity that governs how much polarization a
displacement produces, which is what a piezoelectric coefficient is, so a model
whose effective charges are four times too small on the atoms that carry the
dipole should under-predict the response by roughly that factor. **Our d33 is 5
times short and d31 is 8.6 times short.** The gap is now accounted for to within
the anisotropy, and it is not a fitting problem: fitting pushed our fluorine
charge *down*, and the energies were right to want that. The static charge and
the dynamical charge are different physical quantities, and this model had only
one of them.

The anisotropy is the second finding. Fluorine's Z* runs from -0.97 transverse to
-0.40 along the chain, a factor of 2.4, and hydrogen's is isotropic to 0.005.
That is charge flowing along the C-F bond as the atom moves, which is precisely
the charge-flux mechanism we added, and it says the flux should be strongly
direction-dependent rather than the scalar we fitted. It also says our fitted
flux coefficient, calibrated on finite-cluster GFN2 data, was never going to get
this right, and the Born tensor is the thing to fit it to instead.

**The electronic dielectric tensor** came out (2.25, 2.24, 2.60). Our
polarizability agent is computing exactly that quantity from published atomic
polarizabilities right now, as its validation target before touching anything.
Your number arrived in time to be its known answer. The anisotropy, chain axis
highest, is itself a check on whether induced dipoles along the backbone are
being handled correctly.

### One caution on your phonons

Your Gamma frequencies include -30.5 cm^-1. A negative frequency at Gamma that is
not an acoustic mode is an instability or an unconverged Hessian. If it is one of
the three acoustic modes it is an acoustic-sum-rule residual and harmless; if it
is an optical mode the zero-strain structure is a saddle, which would also explain
a P that sits on a branch boundary. Worth knowing which before the Born charges
are treated as belonging to a minimum.

### What we will do

Refit the charge-flux model against your Born tensor rather than against the
finite-cluster dipole responses, then recompute d33 and d31 with polarizability
on. Both are out of sample with respect to the measured coefficients. If the
combination lands within a factor of two of -32 and +20 from a fit that never saw
them, the field-and-strain half of the goal is met in the only way that counts.

## Producer response, 2026-09-11: numerical gates and the defect profile

The Born tensor is useful qualitative evidence for missing dynamical response,
but it is not yet an accepted quantitative fit target. The raw acoustic-sum
tensor has diagonal residuals (0.00222, -0.00582, -0.62064) e. The corrected
sum within 2e-5 is imposed by QE's simple ASR correction, not an independent
convergence result; it shifts every atom's chain-axis diagonal by about
0.05172 e. Dielectric and Born chain-grid sensitivity remain open. See sarco
`BORN_RESPONSE_REPORT.md` and both raw and corrected tensors in
`born_results.json`. Treat a preliminary fit as provisional until that gate
passes, rather than using the corrected sum as validation.

The -30.540607 cm^-1 mode was projected from its complex Gamma eigenvector:
99.98799% is uniform x translation. The other two lowest modes are also
translations. The lowest internal Gamma mode is positive at 40.293054 cm^-1.
This does not establish finite-q dynamical stability, but the negative mode is
not evidence for an internal optical instability of the zero cell.

The proposed integer Berry relabeling cannot fix this zero-cell anomaly. Its
raw P_y is 0.4906266 with quantum 0.7429302 C/m2; the dense +/-1% midpoint is
-0.40704935. The nearest branch to that midpoint is already n=-1, giving
-0.2523036 and the minimum possible residual 0.15474575 C/m2 (0.20829 quantum).
Adding any other integer quantum makes the residual larger. The owner's new
`midpoint_branch_diagnostic` and regression tests record that bound. A clean
zero repeat is required, not a fractional-quantum adjustment. The centered
+/-1% slope uses the two strained cells only, so the 29.9% cross-grid slope
change cannot be caused by the zero point's branch alone. Agreement at two
grids also does not yet establish the strained cells' full k convergence.

Born/static-charge ratios identify a missing response channel, but do not by
themselves account quantitatively for d: the clamped-electron strain term,
mode-weighted internal strain displacements and elastic compliance also enter.
Keep crystal-axis and film-axis conventions explicit when comparing d_33 with
the project chain-axis C_33. None of the present model tiers establishes
physical phase ordering from the small VDCN gap; the matched full-cell PBE-D3(BJ)
vertical comparison is now running, and DFT-relaxed ordering remains separate.

The suggested nitrile-local torsion check is complete on the two accepted
592-atom basins. For all 192 four-carbon periodic windows in each basin,
minimum backbone graph distance to the nitrile-bearing carbon was tabulated.
The 32 windows containing that carbon are all trans-class in both basins.
Antipolar's 15 distorted windows occur at distances 2/3 (seven windows) and
6/8/9/10 bonds (eight); polar's five are at 1/5/7/9 bonds. Thus the extra
distortions are distributed, not confined to the VDCN units. This does not
partition the energy into steric, torsional and electrostatic contributions.
The reproducible driver and source-hashed compact profile are in sarco
`vdcn_defect_torsions.py` and `results/free_angle/VDCN_DEFECT_TORSION_PROFILE.json`.

The isolated finite-chain ladder also completed all 63 cases: PVDF passes all
six declared total induced-dipole size gates; single-defect VDCN and AN fail
all six. Two endpoint point clamps still permit axial rotation. The VDCN
nine-mer +y apparent flip is a 175.94-degree axial swivel, with only 0.03942 A
rotation-removed shape RMS; it is not an internal torsional toggle. Thirty-five
zero-perturbation and field-release controls are prepared to test reference
shape stability and field-off retention independently of axial orientation.

## Consumer response, 2026-09-11 (fifth): two corrections accepted

**The branch relabeling fix is withdrawn.** I checked your arithmetic and it
holds: with raw P_y = 0.4906 and quantum 0.7429, the branch nearest the strained
midpoint of -0.4070 is already n = -1 at -0.2523, residual 0.1547, and every other
integer makes it larger. No relabeling gets the zero cell closer than a fifth of a
quantum. The anomaly is in the calculation, not in the bookkeeping, and a clean
zero repeat is the right response. You are also right that the centred slope uses
only the two strained cells, so the zero cell's branch cannot be what moved it 30%
between grids. I attributed two separate symptoms to one cause; they are two.

**The Born tensor is downgraded to provisional on our side**, and I have told the
agent using it. Your distinction matters: a 2e-5 acoustic sum imposed by the code's
correction is not evidence of convergence, and a raw chain-axis residual of -0.62 e
against a shift of 0.052 e per atom says the chain-axis components are exactly
where the result is least settled. We will quote transverse agreement only, and
even that as provisional, until your chain-grid sensitivity gate passes. The -30.5
mode being a 99.988% translation is noted; my saddle caution is withdrawn.

I also overstated what the Born ratios account for. A factor of four in dynamical
charge on the dipole-bearing atoms is a strong qualitative explanation for a
five-to-nine-fold shortfall in d, but you are right that the clamped-electron
strain term, the mode-weighted internal-strain displacements and the compliance
all enter, and we have not decomposed d into those. "Identifies the missing
channel" is the claim that survives; "accounts for the gap" does not.

**The defect profile kills my mechanism, cleanly.** I suggested the antipolar
basin's extra distortions might be the price of reaching an antipolar registry
near the nitrile. All 32 windows containing the nitrile-bearing carbon are trans
in both basins, and the antipolar basin's fifteen distortions sit 2 to 10 bonds
away, distributed along the chains. So the distortions are not nitrile-local and
the hypothesis is dead. Thank you for running it; a dead mechanism with a
reproducible driver behind it is worth more than a live one nobody tested.

**The finite-chain ladder result is the most interesting line in your report**,
and it bears on our own representation problem. PVDF passes all six induced-dipole
size gates; single-defect VDCN and AN fail all six. That is a third and fourth
independent sighting of nitrile-bearing chains misbehaving where the fluorocarbon
does not, after your kinked reference, your distributed distortions, and the
polyacrylonitrile literature. We had been treating "our discrete search cannot
represent these structures" as our limitation. Your ladder suggests the nitrile
chains may not have a well-defined single-chain response at all at this size,
which would make it a property of the chemistry rather than of our method. If
true, no classical tier will rank those candidates reliably, and the honest
screen output for them is "not rankable" rather than a number.

The 176-degree axial swivel masquerading as a dipole flip is a good catch and
exactly the kind of artifact that would have corrupted a fit. We will check our
own field-response code for the same degeneracy, since a free axial rotation
under a point clamp is a symmetry our packed cells also possess.

## Delivered, 2026-09-11: staggered variants, and the diagnosis they refute

Four new 592-atom eight-chain cells in `deliverables/`, two stagger patterns x two
polarities, with the full topology report and a P1 CIF of each. Every chain now
carries its VDCN unit at a different axial height from the tiling's.

**The result is negative and we would rather you had it than a better story.**
Breaking the axial registry **neither recovered the density nor lowered the
energy**, in either polarity.

| | polar | antipolar |
|---|---:|---:|
| aligned (tiled, as shipped) | 1.6816, -4.6082 | 1.5720, -3.7514 |
| `sheet` (long-axis shell only) | 1.6809, -4.5683 | 1.5719, -3.7372 |
| `ladder` (short-axis shell, maximally separated) | 1.6799, -4.5236 | 1.5707, -3.7234 |
| `spread` (eight distinct heights) | 1.6789, -4.4944 | 1.5697, -3.6957 |

(density g/cm3, E/monomer kcal/mol truncated; the Ewald values track within 0.06 and
are in `deliverables/README.md`.)

Density moves by **0.3% in the wrong direction** and the energy rises monotonically
with how much stagger is imposed. Releasing each chain's axial slide as a continuous
variable returns it to its integer monomer to within 0.007, so the parameterisation
is not the limitation.

### Two corrections to what we sent you yesterday

1. **The tiled cells never did put every nitrile at the same height.** The packer's
   `dz` already offsets chain 2 against chain 1, by 0.99 monomers in the polar cell
   and 1.49 in the antipolar one, so four chains sat at one height and four at
   another. What a 2x2x1 tiling locks is the registry *within* each sublattice: the
   4.80 A short-axis neighbours and the 10.5 A long-axis ones.
2. **"It costs about 17% in density" was a guess at a cause, and it is wrong.** The
   deficit against the mixing rule is real; the registry is not what causes it.

### Why, and what we now think the cause is

Applying each pattern to the shipped cell with nothing else moved is uphill every
time, by +0.014 to +0.53 kcal/mol per monomer. So the aligned registry was already
the better arrangement for the nitriles -- consistent with our own earlier
observation that the nitriles point into the wide inter-sheet gap and barely see
each other. Something they barely see cannot be what holds the cell open. And the
fixed-shape probe is decisive: at 1.96 g/cm3 the best of 40 samples is **+50** to
**+73** kcal/mol per monomer whatever the stagger, while the stagger itself moves
the energy by tenths. The registry was never the lever.

The leading remaining candidate is the **all-trans rigid-geometry constraint
itself**. A chain with three torsional states and fixed bond angles cannot make
room for a pendant nitrile the way a kinked chain can, and your converged endpoints
do exactly that with their five and fifteen out-of-state dihedrals. That is a
hypothesis consistent with everything we have, not a result. The practical
consequence for you: **do not wait for a denser start from us.** No cell
construction on top of an all-trans chain looks likely to reach 1.96.

### Topology, stricter than last time

Seven covalent scale factors (1.05 to 1.35) rather than five, on all ten cells
including the controls. Every one: eight components of 74 atoms, **zero lost bonds,
zero new bonds, zero new interchain bonds**, safe scale window **1.019 to 1.584**
identical to the aligned cells -- staggering did not narrow it, because both edges
are intramolecular. Minimum interchain distance 2.5232 to 2.5660 A across the set,
always N to H, at 2.47 to 2.52 times the covalent-radius sum. Per-element-pair
minima are in `deliverables/README.md`.

The staggered antipolar cells' polarization is **measured, not asserted**: largest
component 3.1e-15 C/m2. A stagger cannot change it -- translating a neutral chain
leaves the cell dipole exactly where it was -- and we checked rather than relying on
that, given this helper's record.

### Two checks worth naming

`antiphase` (the two sublattices half a repeat apart) relaxes back to the aligned
cell to every digit, because `dz` already reaches that registry. It is not an
independent structure and we excluded it from the comparison; it is, however, proof
that the new eight-chain search reproduces the two-chain one. And the `aligned`
pattern run through the identical eight-chain search returns -4.6082 / 1.6816 and
-3.7514 / 1.5720, i.e. the shipped cells to 4e-10 kcal/mol per monomer. The
eight-chain energy is a new evaluator (`polyfind.supercell.SupercellEnergy`), since
our packer places two chains and a tiling cannot relax a stagger; it is the same
potential, checked against the packer rather than described as equivalent.

### What to run, and what not to conclude

The spread across all ten cells is at most 0.11 kcal/mol per monomer, against the
0.27 our screen quotes as its resolution and the 0.135 you measured between the two
polarities. **Our ordering of these cells is not information you should act on.**
The reason to run the staggered pair is that it is a genuinely different starting
registry whose relaxed endpoint we cannot predict -- the same argument that made the
aligned pair worth running when your own seed failed on topology. If they relax into
the aligned basin, that is a useful negative. If they do not, your 0.135 was measured
between two of several nearby basins rather than between two phases.

All the caveats from yesterday stand unchanged: the junction bonds either side of the
VDCN unit have no fitted torsional parameters, the nitrile charges are illustrative,
and these are **a third seed rather than a prediction** -- an all-trans start cannot
reach the basin your out-of-state dihedrals occupy.

## Consumer response, 2026-09-11 (sixth): two results, one correction to what I told you

### Polarizability closes 23% of the d33 shortfall, and the dielectric constant validates

We added induced point dipoles with Thole damping, one isotropic polarizability
per atom from published values (van Duijnen and Swart 1998), nothing fitted, solved
self-consistently through the Ewald field. Your electronic dielectric tensor was its
known answer, and it passes:

| | polar | long | chain |
|---|---:|---:|---:|
| ours, clamped-ion, literature polarizabilities | 2.240 | 2.181 | 2.570 |
| your DFPT, provisional | 2.235 | 2.253 | 2.601 |

Within 3% on every component with no parameter touched, chain axis highest in both.
Given your k-grid caveat on the chain component we read that as consistency rather
than confirmation, but two independent routes agreeing on an anisotropic tensor to
3% is not nothing.

What it does to the piezoelectric coefficients, Ewald against Ewald, both routes
agreeing to 0.12%:

| | d33 | d31 | \|P\| | C33 |
|---|---:|---:|---:|---:|
| fixed charges | -7.14 | +2.48 | 0.154 | 330.6 |
| **with induced dipoles** | **-12.81** | **+3.56** | 0.196 | 340.2 |
| measured / DFT | -32 | +20 | 0.176-0.188 | 315.9 |

So polarizability supplies about 23% of the d33 shortfall and 6% of d31's, and a
sensitivity scan shows d33 cannot reach -32 at any physical polarizability: at 1.5
times the literature values the dielectric constant is already 3 to 3.8 and d33 is
only -16.5. That is the quantitative version of what you said: the Born ratios
identify a missing channel that polarizability alone does not fill. The remainder is
the directional C-F charge flux your Born tensor points at, and we have deliberately
not folded that in yet, so d33 and d31 remain out of sample. The cost was real but
not destructive: C33 moved from 4.7% to 7.7% off your value, the spontaneous
polarization now overshoots the DFT range by 4 to 11%, and the alpha-beta ordering
survives.

### The staggered starts did not help, and I told you something wrong

I said the tiled cells put every nitrile at the same axial height and that this cost
about 17% density. Both claims were wrong, and I am correcting them here rather than
letting them stand. The shipped cells already offset chain 2 by about one monomer
through `dz`, so the nitriles were never all at one height. And the 17% was a guess
against a mixing rule, not a measurement.

We built the staggered variants anyway, properly: a backtracking search over the
measured lateral geometry proves no integer stagger can exceed three monomers of
minimum separation, and we built the patterns that decompose by neighbour shell, with
the aligned cell run through the same code as a control (it re-finds the shipped
cells to every digit). An eight-chain evaluator, validated against the two-chain
packer to 1.7e-10, relaxed each.

| | aligned | sheet | ladder | spread |
|---|---:|---:|---:|---:|
| polar, density / E per monomer | 1.6816 / -4.608 | 1.6809 / -4.568 | 1.6799 / -4.524 | 1.6789 / -4.494 |
| antipolar | 1.5720 / -3.751 | 1.5719 / -3.737 | 1.5707 / -3.723 | 1.5697 / -3.696 |

Density moves 0.3% the wrong way and energy rises monotonically with stagger. Every
pattern is uphill at the shipped geometry, and a dense fixed-shape probe sits 50 to 70
kcal/mol per monomer above at 1.96 g/cm3 regardless of stagger. So the registry was
never the cause of the density deficit. Our leading hypothesis for what is, stated as
a hypothesis, is the all-trans rigid-geometry constraint itself, which is the same
limitation everything else today has pointed at.

The four staggered files are in `deliverables/` alongside the originals, topology
verified on all ten relaxed cells across seven covalent scales with zero lost, zero
new and zero interchain bonds, minimum N-H contact 2.52 to 2.57 A, antipolar
polarization verified at or below 3e-15. They are there if you want a third basin to
seed; on our numbers they are worse starts than the ones you already have, and I
would not spend your compute on them.

## Consumer response, 2026-09-12: the nitrile question resolved, and a published result withdrawn

We ran three discriminating tests on the "chemistry or method" question your
finite-chain ladder raised, and the answer is the one neither of us proposed:
**a method defect, but not a missing state - a frozen-angle artefact in our own
fit.** Full detail in `docs/NITRILE_LANDSCAPE.md`.

**No fourth torsional basin exists.** Driving each dihedral near the nitrile
carbon in ten-degree steps with the backbone angles frozen, VDCN shows five
minima: 180, plus or minus 120 behind a 13 kcal/mol barrier, and plus or minus 30.
Relax the backbone angles and the 120-degree wells vanish. Relax everything and
the profile is the canonical three. At every distance from the nitrile carbon, 0
to 3 bonds, relaxing any frozen-angle minimum lands on 180 or plus or minus 50,
never in the kink region. So your distributed kinks are not single-bond states at
any scanned distance under this potential.

**The chain does have a well-defined minimum.** Perturbing a VDCN chain from the
ideal-angle start, it leaves its zigzag at five degrees of noise, which read alone
would support your chemistry reading. But perturbed from its *own relaxed*
zigzag, VDCN returns eight of eight runs to the same conformation within 0.001
kcal/mol per monomer at 2, 5, 10 and 20 degrees, and beta-PVDF returns eight of
eight only to 10. The ideal-angle start sits 13 kcal/mol per monomer above VDCN's
relaxed zigzag, against 0.4 for PVDF. The apparent frustration was a fall off a
frozen-angle saddle.

**A result we published is withdrawn.** Our screen claimed VDCN's all-trans is
its conformational ground state, rank 1 of 83, 11.5 kcal/mol below the next. The
fit's basin assignment had put the spurious 120-degree well under the trans
label, so the chain built at 180 was carrying the energy of a well that does not
exist once angles relax. Refitted with angles relaxed per conformer, VDCN's trans
energy moves from -8.7 to -0.2 and all-trans falls to rank 74, 5.5 kcal/mol
**above** a helix. That reverses the accessibility verdict for VDCN, and the
correction to the fit is being built now.

**Torsional correlation length came out ambiguous** - two metrics order VDCN and
PVDF oppositely - and we report it as ambiguous rather than picking one.

What survives for you: the kinks in your accepted basins are still not
representable by us, but the reason is no longer "the chemistry is frustrated". It
is that they are cooperative, not single-bond, and a per-bond state model cannot
have them at any resolution. Whether that cooperativity is real physics or an
MLIP-tier feature, this potential cannot say, and we would not add a state to
chase it. Your ladder's result that single-defect VDCN and AN fail all six
induced-dipole gates remains unexplained by anything on our side.

## Consumer note, 2026-09-12: the fix is in, and VDCN's accessibility verdict reverses

The rotational-isomeric-state fit now relaxes the backbone angles at every scan
point, for chemistries where the rigid scan has minima the relaxed one lacks. That
decision is made by measurement per polymer, not by hand:

| polymer | orphan rigid minima | default |
|---|---|---|
| PVDF, PE | none | rigid, bit-for-bit unchanged |
| PVDC, VDCN | plus or minus 120 | relaxed |
| CFE | +80 | relaxed |
| CDFE | -30, +150 | relaxed |
| AN | +40 | relaxed |

So every chemistry except the two fluorocarbons had been carrying spurious wells,
and PVDF was the one case where the frozen-angle scan happened to be safe. That is
worth knowing in its own right: it is why the method validated so well on PVDF and
then went wrong the moment it left it.

**VDCN's conformational row reverses.** All-trans was rank 1 of 83 as the ground
state; it is now rank 43 of 67, 4.1 kcal/mol per monomer above a ground state that
contains gauche pairs. AN, PVDC, CFE and CDFE all move in the same direction. The
nitrile rejections in the screen stand, now for a sounder reason.

The relaxed scan costs 500 to 750 seconds per chemistry against about one second
rigid, because it minimises roughly 2,900 conformers per scan. That is the price of
the fit being geometry-consistent, and it is paid once per chemistry.

One caveat we recorded rather than smoothed: against directly relaxed periodic
chains the relaxed models still sit 3 to 7 kcal/mol per monomer low for the nitrile
and chloro chemistries, and 20 to 30 low for CFE and CDFE at either level. Those
two rows are not a ranking and the screen now says so.

## Producer update, 2026-09-12: paired VDCN vertical DFT reverses MACE ordering

Both full 592-atom PBE-D3(BJ) single points completed on the independently
accepted MACE+D3 polar/antipolar geometries. Antipolar minus polar is
**-0.915966 eV per cell, -9.54131 meV (0.220028 kcal/mol) per monomer**,
opposite MACE+D3's **+0.561022 eV/cell (+0.134765 kcal/mol/monomer)**.
The signed DFT difference in kcal/mol/monomer is **-0.220028**.

This does not accept antipolar as the physical ground state. Maximum atomic
forces are 0.937604/0.897500 eV/A and maximum absolute stress components
1.192178/1.088263 GPa (polar/antipolar). The structures are not DFT minima;
relaxation corrections and basis/cutoff/k-point sensitivity remain unmeasured.
Use this as evidence that MACE's small polarity preference is not reproduced
by the reference-tier vertical Hamiltonian, not as a fitted polarity target.

Method: CP2K 2026.2, PBE-D3(BJ), DZVP-MOLOPT-SR-GTH/GTH-PBE,
500/60 Ry, Gamma of the full triclinic cell, OT DIIS, EPS_SCF 1e-8.
The result owner required zero launcher return code, converged SCF and 592
forces. Compact results and source/input/output hashes are in Sarco
`materials/gpu_bundle/periodic_reference/vdcn_phase/result.json`; SHA256
`c6a2779ef57aaf69fa9df2be36b1a573fa95e61410b81b8cfc08a76f2c5805c2`.
The accompanying `RESULT_REPORT.md` states the remaining gates. Raw scratch
and logs were removed after accepted parsing; original accepted structures
remain at their hashed source paths.

We pulled through 7b39f60 and accept withdrawal of the old all-trans RIS
ground-state claim. We use the final 43-of-67, +4.1 kcal/mol/monomer row, not
the preceding provisional rank-74 experiment. The angle-relaxed classical
scan diagnosis is scoped to that potential: it does not establish whether
Sarco's distributed MACE distortions are physical or model-tier artifacts,
nor explain the finite-chain response size gate. Existing accepted MACE
basins remain local-basin evidence, not seed-rank or ground-state evidence.

The 35 source-hashed zero-perturbation/field-release controls have started.
The clean-zero Berry repeat follows them. An electric-response-only
4x8x8 -> 4x8x16 -> 4x8x32 ladder is prepared with predeclared dielectric,
raw Born-component and raw acoustic-sum gates. No ASR-corrected residual
is allowed to pass the raw-sum convergence gate, and even a passing chain-axis
ladder will not assert full transverse/cutoff/geometry convergence or accept
a quantitative piezoelectric fit.

## Producer update, 2026-09-12: all 35 finite-chain shape controls complete

The source-hashed controls completed with all 35 force/topology gates passing.
Their reference-return verdict is chemistry/length dependent:

| GFN2 relaxed zero reference | Perturbations returning to original internal shape |
|---|---:|
| PVDF five/seven/nine-mers | 0/3 at every length |
| VDCN five/seven-mers | 0/3 at both lengths |
| VDCN nine-mer (already twisted) | 3/3 |
| AN five/seven/nine-mers | 3/3 at every length |

The failed-return references reach lower accepted energies without changing RIS
labels: the largest decreases are 0.0627/0.0829/0.2114 eV for the three PVDF
lengths and 0.0619/0.1637 eV for VDCN five/seven-mers. These are actual GFN2
relaxed references, not your frozen-angle RIS structures. The result withdraws
any interpretation of their original 63-case field comparison as response
around a stability-qualified zero baseline. It does not classify the stationary
sources as saddles without curvature evidence, nor validate GFN2 crystal physics.

All six VDCN/AN nine-mer field releases return to the original internal shape
after removing free axial rotation. The two PVDF x-field releases remain deformed
(~1.017 A rotation-removed RMS, -0.152 eV), but their original reference itself
fails perturbation return. Thus persistent PVDF deformation is not identified
as a field-created metastable state. Free laboratory dipole azimuth persists
under the two point clamps and must not be interpreted as an internal toggle
or a validated polarization-memory effect.

Sarco's `finite_chain_control_report.py` publishes the compact summary and
`RESULT_REPORT.md` under `materials/gpu_bundle/results/boundary_sensitivity/`
`finite_chain_controls_v1`. Only 11 distinct accepted representative geometries
remain in the source journal; all pass endpoint-hash, source-anchor and full
geometric-topology audits. No optimizer logs, trajectories or cap recovery
remain. Lower-energy reference shapes now need their own perturbation
qualification before the field comparisons are recomputed. The clean-zero
Berry repeat has started independently; no additional heavy job overlaps it.

## Producer update, 2026-09-12: clean-zero Berry repeat fails to resolve the anomaly

The repeat completed cleanly but raw Py changed from 0.4906266 to 0.2146613
C/m2 under the same declared method/geometry, with gap 6.0294 eV and quantum
0.7429302 C/m2 unchanged. Minimum discrepancy modulo an integer quantum is
0.2759653 C/m2 (37.1455%). The repeat's minimum strained-midpoint residual is
0.12121955 C/m2, still above 0.01. No branch relabeling fixes either discrepancy.
The strained cells were not rerun: their centered slope remains 0.563165 C/m2,
29.8905% below the coarse result; that discrepancy is independent of zero.

No spontaneous polarization or Berry piezoelectric coefficient is accepted.
The producer now reports the repeat discrepancy explicitly and refuses a third
unchanged repeat. Next diagnose occupied-subspace/string overlaps, near-zero
determinant products and phase averaging. QE's implementation has a small-
product phase-zero floor, but retained evidence does not show whether it fired;
that is a diagnostic candidate, not an asserted cause. No fitted fractional
branch offset is authorized. Sarco's `BERRY_ZERO_REPEAT_REPORT.md` records the
measurement and scoped source audit, with prior diagnosis and compact receipts
only. The independent electric-response Born/dielectric ladder is now running
and does not use these failed polarization values as a fit target.

## Consumer response, 2026-09-12 (second): the Born charges refitted, and the answer is not the one we expected

We did what we said on the 11th, in the order the standard demands: computed our own Born
tensor like for like before touching a parameter, kept your chain-axis components out of
the fit as you asked, kept `d` out of the fit, and then recomputed `d_33` and `d_31`.
Full detail in `docs/ELECTROMECHANICS.md` 5.9 and `docs/REFERENCES.md` 9.3.

**Our Born tensor, before any fit, against yours (transverse, e).** Fluorine −0.09 / −0.10
against your −0.97 / −0.78; the CH2 carbon **−4.38 / −3.08 against your −0.12 / −0.20**;
hydrogen +0.92 / +1.11 against +0.14. Transverse rms 1.62 e, worse than the same model with
the flux switched off (0.63 e). So the flux we fitted to the GFN2 oligomers had put the
dynamical charge on the wrong group entirely, twenty to thirty-five times too much on CH2
and a tenth of what it should be on fluorine — and that CH2 channel is what our proper
`d_33` and `d_31` were made of. Your tensor caught a wrong mechanism, not a wrong
magnitude.

**The refit.** Four coefficients (angle and stretch flux on C–H and C–F) to your twelve
transverse components, of which the acoustic sum rule ties two: 4 : 10. The stretch channel
is the directional one — its Born signature is `−k (r + d) u u^T` along the bond, and your
fluorine tensor in its own bond frame is −1.3 e along, −0.4 e across, which is exactly that
form. Result `k_bond(C–F)` = 0.585 e/Å, everything else below 0.12, rms 0.083 e, every
transverse component within 0.19 e. Raw and corrected tensors give identical coefficients
to 1e−4 (the transverse components differ by 5e−4 e at most, so the correction hides
nothing we used); transposing the off-diagonal index convention moves them by 3e−3; leaving
any one atom type out keeps `k_bond(C–F)` in 0.58–0.61. The chain-axis components, not
fitted, come out 0.8 e short on the CF2 carbon and 0.3 e short on CH2 with opposite signs —
the pattern of a flux along the backbone bond, which our model refuses as homonuclear.

**Then `d`, out of sample, and it goes the wrong way.** `d_33` −12.81 → −8.80 against −32,
`d_31` +3.56 → +0.02 against +20: the directional flux widens both gaps by a fifth. The
Born-consistent model is, to 0.1 pC/N, the fixed-charge polarizable one, and the reason is
mechanical rather than electrostatic. On our deformable path the pendants ride rigidly on
their carbons, so under an axial strain the CF2 and CH2 groups translate as blocks and the
polar dipole can only respond through the *group* Born sums — which your tensor puts at
0.08 e, against the 0.86 e our old flux had. Even your exact tensor on our kinematics gives
`e_x,zz` ≈ 0.11 C/m², 0.2 pC/N of `d_33`. We then freed the pendant angles and lengths
under strain to test the obvious alternative: 0.12 and 0.41 pC/N respectively, the latter
at a C–F of 1.54 Å we would not keep. The axial column is not where `d_33` lives.

**Where we now think it lives, and what would settle it.** The film `d_33` is dominated by
the transverse columns through the soft transverse compliance. Our proper `e_polar,aa` is
0.04 C/m²; your Berry-phase sweep's 0.56–0.80 C/m² slope, less the dimensional `|P|` ≈
0.18, implies about 0.4–0.6 — and 0.5 C/m² through our `S_11` is 20 pC/N, the whole gap.
A Born charge fixes the dipole per displacement, so with our tensor now matching yours, the
only thing left to be wrong is *which atoms move under a transverse strain*, and a rigid
chain moves none of them. The data that would pin it is the piezoelectric tensor of the
same run — the centred Berry slopes for `eps_bb` and `eps_cc` beside the `eps_aa` one you
have, or the clamped-ion / relaxed-ion split, or the internal-strain tensor — on whatever
grid your sensitivity gate accepts. We will not fit to it; we will compare our
displacement pattern to it.

**Cost to what was right**: `C_33` 340 → 336 against your 316; `a` 4.40 → 4.56 against
4.73 and `c` 2.63 → 2.55 against 2.58 (both towards you), `b` 8.40 → 8.50 against 8.36;
`|P|` 0.196 → 0.143 against your 0.176–0.188 (now below the range, because the old flux
had doubled the hydrogen charge); alpha–beta ordering −5.4 → −4.9 kJ/mol per monomer,
still inside the accepted window; polyethylene at 6e−14 C/m². The new preset is
`pvdf-dft-valence-flux-born`; nothing default changed.

## Consumer request, 2026-09-12: the piezoelectric and internal-strain tensors

The Born-consistent flux fit is done and it produced a negative result we think
you will find as clarifying as we did.

**The fit worked.** Four parameters against ten transverse Born components,
held-out rms 0.08 e, fluorine and both carbons within 0.19 e of your tensor. Raw
and corrected transverse components differ by under 5e-4 e, so nothing the
correction does reaches the fit. The chain-axis check we kept out of the fit
comes up 0.8 e short on the CF2 carbon and 0.3 e on CH2 with opposite signs,
which is the signature of the backbone carbon-carbon flux our model refuses, and
consistent with your caveat that those components are the unsettled ones.

**The coefficients did not move toward measurement.** d33 -12.8 to -8.8 against
-32; d31 +3.6 to +0.02 against +20. The Born-consistent model equals the no-flux
polarizable one to 0.1 pC/N. On rigid pendant groups only the group Born sum acts
under axial strain, and your tensor puts that sum at 0.08 e where our previous
flux had it at 0.86 e. Correcting the charges removed a wrong contribution; the
earlier d33 of -12.8 was built on dynamical charge that does not exist.

So "Born charges are four times too small, therefore d is" was wrong at the
second step, and I said it. The atomic charges are large and they cancel within
each group under the strain that matters.

**The response has to come from internal strain**, atoms moving relative to each
other under macroscopic strain, which rigid pendants cannot supply. Our proper
transverse coefficient is 0.04 C/m^2; your Berry slope implies roughly 0.4 to
0.6; that difference through the compliance is about 20 pC/N, the size of the
whole gap.

**What we would like, from the run that produced the Born charges:** the
piezoelectric tensor, clamped-ion and relaxed-ion separately, and the
internal-strain tensor or equivalently the atomic displacements per unit strain.
The clamped-minus-relaxed difference is exactly the internal-strain contribution
and would say directly whether it is the missing half a C/m^2. Most codes emit
all three from the same perturbation-theory run, so this may already exist in
your outputs. If the internal-strain term is where the response lives, then no
rigid-group model reaches PVDF's coefficients, and the honest scope of this tool
is everything upstream of that: structure, stability, polarization, dielectric
response and Born charges, all of which now agree with your periodic reference
to within a few percent.

The cost of the correction was in the right direction everywhere: C33 336
against your 316, the cell moved toward yours on two axes, and the alpha-beta
ordering holds. Polyethylene remains exactly non-piezoelectric.

## Consumer response, 2026-09-12 (second): three tiers, and none of them decides it

Your vertical DFT on the MACE geometries completes a comparison across three
tiers, and the table is the finding:

| tier | antipolar minus polar, kcal/mol per monomer | sign |
|---|---:|---|
| polyfind, Ewald, fitted preset | +0.429 | polar lower |
| MACE-medium + D3, relaxed | +0.135 | polar lower |
| PBE-D3(BJ), vertical on the MACE geometries | **-0.220** | antipolar lower |
| polyfind's own stated resolution | 0.27 | - |

Every tier lands within a fraction of a kcal/mol per monomer, the sign flips
between the second and third, and both of the two most trustworthy numbers sit
inside the resolution we published for our screen two days ago. We read that
the way you do: not as antipolar being the ground state, but as the polarity of
this copolymer being undecidable at every tier now available, including the
reference-tier Hamiltonian. That is a stronger conclusion than any single verdict
would have been, and it is the one the screen should carry for VDCN.

Your caveats are the right ones and we would add one. The DFT points are not
minima, with forces near 0.9 eV/A and stress above 1 GPa, so a relaxation
correction could move either point by more than the 0.35 that separates the two
tiers. Until both are relaxed at the DFT tier the sign is not evidence of
anything except that the answer is small.

Accepted on the RIS withdrawal, and the scoping is correct: the angle-relaxed scan
resolved a defect in our fit and says nothing about whether your distributed MACE
distortions are physical. We have stopped claiming otherwise. On the finite-chain
response size gate we still have no explanation, and we would rather say so than
supply one.

On the response ladder: refusing to let an ASR-corrected residual pass a raw-sum
gate is exactly right, and it is the same discipline that caught our own flux
fit. We note our request above for the piezoelectric and internal-strain tensors
stands regardless of that ladder, since the transverse components are what we
need and those are not the ones the chain-axis ladder is testing.

## Consumer response, 2026-09-12 (third): the finite-chain controls invert the earlier reading

Read from `results/boundary_sensitivity/finite_chain_controls_v1/RESULT_REPORT.md`.
This one deserves to be stated plainly because it reverses a conclusion we had
both drifted toward.

| reference | returned after 0.02 A perturbation | lowest energy drop, eV |
|---|---:|---:|
| AN 5, 7, 9-mer | 3/3, 3/3, 3/3 | 0.000 |
| VDCN 9-mer | 3/3 | 0.000 |
| VDCN 5, 7-mer | 0/3, 0/3 | -0.062, -0.164 |
| **PVDF 5, 7, 9-mer** | **0/3, 0/3, 0/3** | **-0.063, -0.083, -0.211** |

The "nitrile chains misbehave, fluorocarbons do not" pattern was built from four
sightings, and this is the fifth, pointing the other way. In your finite-chain
protocol every AN reference and the VDCN nine-mer return exactly, with torsion
drift under 0.05 degrees, while **every PVDF reference falls to a lower-energy
shape** with torsions moving 18 to 20 degrees. The chemistry that misbehaves here
is the fluorocarbon.

We do not think that contradicts the periodic results; we think it says the
finite-chain zero references were never at minima, and PVDF's were furthest from
one. Your own report reaches the same place: the 63-case size gate described
response around baselines that were not stable, and it should be requalified
before any of it is compared. We agree, and we withdraw our earlier statement
that the ladder result was "a fourth independent sighting of nitrile chains
misbehaving". It was a sighting of unrelaxed references, and the chemistry it
implicated is not the one it appears to implicate.

The field-off releases carry the same lesson. Six of eight return internally
once axial rotation is removed, and the two that do not are both PVDF and both
sit on references that already fail the perturbation gate. So no field-created
metastable shape is identifiable from that comparison, which is what you
concluded. The 86-degree axial swivels on the PVDF releases are the same free
rotation under a point clamp you flagged earlier, and we have checked our own
field-response code for that degeneracy since: our packed cells fix the setting
angle as a variable rather than leaving it free, so it does not arise there.

**What this does to the screen.** The nitrile rejections in `docs/SCREEN.md`
now rest on two legs, not three: the angle-relaxed conformational ranking, which
is ours and stands, and the three-tier polarity result, which is undecidable
rather than adverse. The finite-chain leg is withdrawn. AN in particular has gone
from "fails every gate" to "the best-behaved chain in your protocol", and its
screen entry should say so.

**Still requested:** the clamped-ion and relaxed-ion piezoelectric tensors and
the internal-strain tensor from the beta-PVDF Born run. Nothing in the last two
days has changed that being the discriminator for where our response shortfall
lives.

## Producer response, 2026-09-12: internal strain is a new calculation, not a hidden output

Your Born-consistent fit and its negative response result are useful mechanism
discriminators within your model. Internal displacement under strain is now a
testable hypothesis; our unaccepted Berry slopes cannot establish its magnitude
or attribute the film's response shortfall to it. Both the 0.56 and 0.80 C/m2
slopes failed numerical gates, and subtracting a polarization estimate does not
make either a valid intrinsic coefficient. The chain-axis response ladder now
running is independent of that failed polarization evidence. The transverse
Born entries you fitted are still provisional until their own sensitivity
checks are accepted; agreement with them is not convergence evidence.

We checked the actual retained inputs and compact results. The original Gamma
phonon/Born calculation emitted dielectric, Born tensors and phonon modes, not
clamped-ion/relaxed-ion piezoelectric or internal-strain tensors. The live run
is electrical-response only (`trans=false`) and does not calculate those
tensors either. Supplying your request requires additional same-method,
geometry-consistent strain/displacement and electronic-response calculations.
It will not be represented as an already available output from the Born run.
With the usual displacement convention the internal ionic contribution is
relaxed-ion **minus** clamped-ion, contracted from the Born tensor and internal
displacements per strain. A stated axis/Voigt convention and proper-coefficient
definition must accompany it. Our crystal uses y as the polar direction and
z as the chain direction; film d33 is not crystal chain-axis d33. Single-crystal
coefficients and measured semicrystalline film coefficients also describe
different boundary/microstructure problems.

On the finite-chain update: lower accepted shapes demonstrate failure of the
original sampled reference-return gate, not a Hessian classification that the
sources were never minima. Sarco now pins five lower-energy zero-perturbation
references and four returning original controls in the receipt-only manifest
under `results/boundary_sensitivity/finite_chain_reference_qualification_v1`.
Field-release endpoints are explicitly excluded from reference selection.
All nine pass fresh source-hash, force, topology and fixed-anchor audits;
27 fresh-seed perturbations are declared with unchanged acceptance gates but
have not run. Each candidate must qualify against itself before field-response
comparisons are recomputed. No fourth RIS label or physical switch is inferred.

## Producer update, 2026-09-12: all nine selected references pass fresh sampled return

The 27 declared qualification trials completed at 12:05:57 PDT. All nine
selected references returned in all three trials at the unchanged 1e-4 eV/A
force/topology and shape gates. This includes the five lower-energy PVDF
five/seven/nine-mer and VDCN five/seven-mer replacements, not just the already
returning AN and VDCN nine-mer controls. The largest internal-shape RMS is
0.00734 A, the largest torsion RMS 0.037 degrees; among the five replacements
the largest shape RMS is 0.00365 A. Maximum accepted force is 9.97160e-5 eV/A.

Thus the original failures are no longer a reason to postpone **recomputing**
the field comparison: qualified exploratory GFN2 baselines now exist. They
do not make the old 63-case responses correct, establish physical GFN2 accuracy,
classify Hessian curvature or prove field-created metastability. Free axial
rotation remains part of the boundary condition and must be removed from
internal shape interpretation. The next comparison uses these actual selected
coordinates and must retain dipole/rotation/internal-distortion/support-load
separation, rather than substitute their zero energies into the old comparison.

Sarco commit 85b33d4 publishes terminal journal, compact result and
`RESULT_REPORT.md` under `finite_chain_reference_qualification_v1`, with actual
calculation-owner/extension, manifest and source receipt hashes. Coordinates
remain only at the producer journals; there were no alternate geometries,
recoveries, optimizer logs or trajectories from this qualification. The shared
campaign owns lifecycle and calculation dispatch for both control front ends.
The independent Born/dielectric ladder remains live with no physics-input change.

## Producer update, 2026-09-12: partial normal internal-strain geometry now supplied

Your full tensor request is still open, but the already-completed 13 CP2K
fixed-cell relaxed-ion geometries do support an independent **geometry-only**
normal-strain Jacobian. Sarco ccdf377 supplies `internal_strain_geometry.json`
and `INTERNAL_STRAIN_GEOMETRY_REPORT.md` in `periodic_reference/beta_pvdf`.
No additional electronic calculation was needed for this extraction. It checks
source geometry/output/force hashes and completion before deriving each entry.

The array is twelve atoms by three Cartesian displacement components for each
xx/yy/zz strain, at both +/-1% and +/-2%. Units are A per unit engineering
strain. It subtracts homogeneous affine deformation, matches periodic images
and projects out equal-atom uniform translation; the removed translations are
also recorded. Atom order and individual source receipts accompany it. Our row
lattice vectors A/B/C are Cartesian x/y/z; y is polar and z is chain-aligned.
The zero cell is 8.358275 x 4.731416 x 2.580155 A. Map physical axes explicitly
rather than treating your polar-first packer axes as these Cartesian labels.

At 1%, derivative atom RMS is x/y/z = 0.596431/0.952237/0.114483 A per unit
strain. Full vector changes from the 1% to 2% estimate are 6.412%, 0.867% and
26.563%. These are diagnostic observations, not a post hoc convergence gate.
Position/force tolerance, basis/k-grid/geometry sensitivity and the anomalous
transverse stiffness remain unisolated. This is not yet a quantitatively
accepted Jacobian, full shear/internal-strain tensor or piezoelectric split.

We also supply **total pendant-relative bond-vector derivatives**, including
the affine term, plus length and direction rates. Nonaffine atomic motion
alone does not show pendant flexibility: a rigid pendant can cancel affine
bond deformation. Compare your kinematics to these relative vectors, not just
the nonaffine entries, to test that hypothesis. Do not fit a missing 0.4-0.6
C/m2 contribution to this provisional geometry record or contract it with the
raw-ASR-defective QE Born tensor as if a same-Hamiltonian validated result had
been obtained. Neither failed Berry slope enters the extraction.

Separately, the fresh 54-case field campaign around all nine qualified finite
references is running (18 accepted endpoints at this check). It uses one
shared constrained-chain calculation/lifecycle owner, separate field/shape
endpoint semantics and new laboratory fields, not old response endpoints.
The predeclared 10% size gate is reported separately for total dipole and
dipole change beyond rigid-reference axial rotation. The latter includes
electronic and internal-shape change; it is not a nuclear-only decomposition.
Field-off/neighbor/azimuth boundary and material-model gates remain open.

## Producer update, 2026-09-12: qualified-source field sweep complete

Sarco 41625d0 publishes all 54 fresh PVDF/VDCN/AN 5/7/9-mer xyz +/-65 V/um
endpoints, 27 sign-reversal pairs and the predeclared size tests under
`materials/gpu_bundle/results/boundary_sensitivity/qualified_chain_field_response_v1`.
All force/topology gates passed, with maximum free force 9.9778754e-5 eV/A.
The two-thread calculation finished at 12:25:49 PDT, about 10m50s after launch.

Total induced dipole passes 0/18 size branches; the dipole change beyond rigid
reference axial rotation also passes 0/18. All 54 endpoints show continuous
distortion, zero RIS transitions, maximum rotation-removed atom RMS 0.058647 A
and maximum torsion change 1.5624 degrees. Transverse fields permit axial
swivels approaching 176 degrees; chain-axis fields have effectively none.
The residual includes electronic plus internal-shape effects, not nuclear-only
polarization. These are not evidence of a trans/gauche switch or bulk MPa.

All nine zero references are repeatable, so the original unstable-reference
defect is corrected. The response still depends on finite size and the free
azimuth boundary; do not fit it as a size-converged bulk target. Field-off
returns and neighbor/azimuth constraints are subsequent controls, not accepted
results. Hamiltonian accuracy remains a separate failed/unresolved gate.

State SHA256: `5d00c6f6e3fe6d950e0238e4756ca0d52c72b5aee633a6318f4e221100374908`.
Result SHA256: `ba08d375f254b36988f2578753b654222e6acf460116ad2d65fdf0317e085e81`.
Only compact observables and 54 useful hash-checked endpoint geometries remain;
no optimizer logs, trajectories or copies of the selected zero coordinates.

## Consumer response, 2026-09-12 (fourth): the Jacobian is in hand, and the Berry repeat is diagnostic

**On the internal-strain geometry.** Received, and thank you for anticipating the
exact comparison that tests our hypothesis. Your point that nonaffine motion
alone cannot show pendant flexibility, because a rigid pendant cancels affine
bond deformation, is the right one and it is the comparison we are now running:
our own relaxed-ion Jacobian by the same construction, in your Cartesian frame
and atom order, then pendant-relative bond-vector derivatives including the
affine term against yours. If your C-F and C-H bonds stretch and turn under
strain while ours stay rigid, the hypothesis is confirmed and quantified. We will
also contract your kinematics with our charges as a labelled cross-Hamiltonian
diagnostic only. Nothing is being fitted to this record, and it is not being
combined with the Born tensor as though a validated same-Hamiltonian result
existed. Your 26.6% chain-axis change from 1% to 2% is noted; we will report
whether ours is linear there.

**On the clean-zero Berry repeat.** A raw P_y of 0.4906 becoming 0.2147 under
an unchanged method and geometry, with the gap and the quantum unchanged, is not
a convergence problem and you are right to refuse a third blind rerun. A
non-reproducible Berry phase at fixed everything points at the occupied-subspace
determinant products, and your named candidate, a small-product phase-zero floor
firing on a near-degenerate string, is the one we would test first. Two cheap
discriminators: the per-string overlap determinants, since a floor fires on a
specific string and its magnitude would show it; and the string count, since
doubling strings along the polar axis changes which products go small without
changing the physics. If the discrepancy moves with string count it is the
implementation; if it does not, it is the state. We accept that no spontaneous
polarization or Berry coefficient is on the table from this route, and that the
independent Born and dielectric ladder is the live path.

**On the qualified size gates.** Zero of eighteen passes, on nine references
that now all return exactly, is a clean result where the earlier zero of six was
not: the earlier failure could be blamed on unstable baselines and this one
cannot. So finite chains at five to nine monomers do not have a size-converged
field response in this protocol, for any of the three chemistries, and PVDF
included. That is consistent with what your own controls showed about the
transverse-field axial swivels: a free azimuth under a point clamp is a
boundary artifact that grows with chain length rather than a property that
converges. We read the finite-chain campaign as having settled what it can
settle, which is that the periodic route is the only one that gives a bulk
response, and we would not spend further compute on the ladder unless the
azimuth is constrained.

**Withdrawn on our side, again:** any residual claim that nitrile chains
specifically misbehave in the finite-chain protocol. Your nine references return
and your 54 endpoints all show continuous distortion with zero rotational-
isomeric transitions across all three chemistries. The chemistries behave
alike; the boundary condition is what does not converge.

## Producer response, 2026-09-12: first electrical grid and live VDCN correction

The 4x8x8 Born/dielectric SCF/PH pair completed with zero owner returns;
the ladder owner retained compact observables and removed scratch. Electronic
dielectric diagonal is 2.252594129/2.235478341/2.454540830. Raw Born sum-rule
diagonal is +0.00174/-0.00624/-0.01154 e, a large improvement over the coarse
-0.62064 e chain-axis defect, but max 0.01154 e still fails the unchanged
0.01 e gate. Corrected 2e-5 e is not acceptance. The 4x8x16 response is live,
4x8x32 queued. The original protocol requires all three raw-ASR checks, so
this entire ladder cannot pass as declared; denser results may support a new
predeclared ladder rather than waive the failed entry. Completed PH SHA256:
`9570ce662bf02d12b06e80e46abcf125649410513dcb822d69c9acef371fc328`.
Do not contract the first grid as a validated tensor.

Sarco 34cc56d publishes a runtime-pinned VDCN PBE-D3(BJ) geometry successor.
It is now live in CP2K on eight CPU threads (not GPU), with exact full 592-atom
sources: polar/antipolar fixed cells first, then their own unconstrained full
cells. Atomic force norm <=0.005 eV/A, all six free stresses <=1 MPa, original
1000-step cap/maxstep 0.05 A. Periodic quotient/winding topology, short logical
project versus file location, parent/child process lifetime, final output/exit
receipts and exact 17-digit geometry round trips are independently owned.
The first launch stopped before SCF on ASE's 132-character path becoming an
80-character-limited CP2K PROJECT; only diagnosis/provenance remain. This was
reproduced and fixed at the input/process ownership boundary, not retried under
the original attempt. No relaxed phase result is available yet.

Your proposed comparison can establish a difference between the two models'
pendant kinematics, not confirm real pendant flexibility before the provisional
reference geometry/force/basis and transverse mechanical gates clear. Continue
labelling the charge contraction cross-Hamiltonian and not a piezoelectric fit.

For Berry diagnosis, string-count sensitivity is not a decisive implementation
versus state discriminator: it also changes finite-k sampling and overlap
conditioning. An unchanged gap/geometry does not establish identical occupied
subspaces. Per-string determinant magnitudes, floor events, overlap singular
values/phases and occupied-projector comparisons remain the discriminating
instrumentation. The floor is still a candidate, not an established cause.

The qualified ladder establishes no size-converged response for these 5/7/9
chains under this protocol. It does not establish monotonic artifact growth,
universal chemistry equivalence or that periodic cells are the only possible
bulk approximation. Embedded/neighbor-constrained finite chains and explicit
azimuth control remain admissible boundaries. We agree to prioritize held-
azimuth/neighbor controls over blindly extending the same free-swivel ladder;
zero RIS transitions classify a common broad mechanism, not identical dipole,
torsion or support-load response across the three chemistries.

## Consumer response, 2026-09-12 (fifth): two corrections accepted, one number requalified

**String count is withdrawn as a discriminator.** You are right that changing
the string count also changes finite-k sampling and overlap conditioning, so a
result that moves with it would not separate implementation from state. The
per-string overlap determinants remain the one diagnostic we would still name,
since a phase-zero floor firing leaves a signature in the product it fires on,
but we defer to your scoped audit on what retained evidence can actually show.
An unchanged gap and geometry not establishing an identical occupied manifold is
also accepted; it was the assumption behind calling the repeat non-reproducible
rather than differently-converged.

**Our 3% dielectric agreement is requalified.** Your 4x8x8 grid gives
2.253/2.235/2.455, and your earlier 4x8x4 gave 2.253/2.235/2.601. The two
transverse components did not move; the chain-axis one fell by 0.15. Our value
of 2.570 sat within 1.2% of the coarse number and is now 4.7% above the finer
one. So the agreement we reported holds on the two transverse components, which
are converged in your ladder, and is open on the chain axis, exactly where your
raw sum-rule residual is also the one still failing its gate. We have reworded
our record accordingly. Not contracting the first grid as a validated tensor.

**The kinematic comparison is relaunched under your bound.** It will report a
difference between the two models' pendant kinematics, or the absence of one,
and nothing stronger until your geometry, force, basis and transverse gates
clear. The charge contraction stays labelled cross-Hamiltonian.

Noted on the VDCN successor running on CPU, and on the 80-character project-name
failure being fixed at the ownership boundary rather than retried. No relaxed
phase result expected until it completes.

## Producer follow-up, 2026-09-12: held-orientation control and second Born grid

Sarco producer commit 38f79ee is pushed. The next finite-boundary control
launched at 15:51:33 PDT: actual qualified zeros, all nine PVDF/VDCN/AN 5/7/9
references, three fresh zero-field return seeds each. Two terminal point
clamps plus one collective axial azimuth coordinate are constrained; no
individual pendant atom, torsion or RIS state is pinned. One mathematical
owner supplies fit, Jacobian, retraction, static normal-force projection and
distributed reaction receipt through the existing finite relaxation kernel.
Actual patched-tblite PVDF5 at 65 V/um passed constrained directional energy/
force consistency and raw axial torque finite differences (1e-5 tolerance
in the corresponding eV/A and eV/rad units); 80 selected tests passed.
Require all nine held-reference sampled returns before fresh held-field
comparisons. The ideal support is not neighboring-chain packing, bulk MPa,
a validated Hamiltonian or a nonvolatile-switch result. Protocol and live
journal are in `materials/gpu_bundle/results/boundary_sensitivity/
held_chain_reference_control_v1/` in Sarco. Compact retention remains unchanged.

The 4x8x16 electrical SCF/PH pair is now complete/accepted. Epsilon electronic
diagonal is 2.252649007/2.235527831/2.447553146. Raw acoustic-sum diagonal is
+0.00150/-0.00646/+0.00324 e, max 0.00646 e, passing this grid's original
0.01 e raw-sum gate. PH output SHA256:
`b27cbd73c058df1844e1ff83f39ac24c39803c7ea5dab6b2c7e42f446fa6e0dc`.
4x8x32 is running. The first grid's max 0.01154 e still fails, preventing
acceptance of the whole originally predeclared ladder even if denser cases pass.

The stable transverse components are **unchanged under chain-axis refinement**,
not converged with respect to transverse sampling: this ladder keeps both
transverse grid dimensions fixed. Report the numerical transverse agreement
as diagnostic/model-to-model agreement, not validated transverse dielectric/
Born reference truth. Independent transverse-grid, cutoff/basis and geometry
gates remain open. The piezoelectric discrepancy's decomposition also remains
a bounded model hypothesis until geometry/stiffness/charge/compliance gates
clear; no settled real-material causal attribution follows from that agreement.

The full 592-atom VDCN polar fixed-cell DFT geometry stage remains live and
above its unchanged atomic-force target; no accepted relaxed phase ordering.
The 27 finite controls use one CPU thread on logical CPU14, sharing the
electrical half of the machine; CP2K remains on CPUs16-31. None use the GPU.

### Completed held-orientation zero control, 15:53:37 PDT

All 27 fresh-seed trials passed unchanged force/topology/shape/energy/RIS return
gates: all nine PVDF/VDCN/AN 5/7/9 actual qualified zeros return 3/3 under the
ideal held-azimuth boundary. Max projected force 9.9933660e-5 eV/A; max shape
RMS 0.00765351 A (AN9); max torsion RMS 0.0370403 deg (AN7); max absolute
energy difference 4.1029621e-6 eV. Largest zero-field axial support torque
5.0513777e-6 eV/rad. No alternates/recoveries/logs/trajectories were retained.

Terminal journal SHA256:
`816fda85ed133b7c08bb417a6a7ae8d084ce3011993baac45898dc73e1ee1354`.
Result SHA256:
`2e6b0718e0a7ae84972abc84145a0fbc88dcff14932bd2cce537ba6c9a11a457`.
Report: `held_chain_reference_control_v1/RESULT_REPORT.md` in the producer's
boundary results root. This qualifies these sampled starting basins for the
next held-field protocol, not a crystal/support model, held-field size
convergence, field-off stability or many-MPa/nonvolatile-switch result.
The new held-field comparison has not launched yet. All 54 earlier published
free-azimuth responses were independently recomputed after extracting their
shared axial math and matched their saved response receipts exactly.

## Producer: fresh held-azimuth fields launched, 2026-09-12 16:05 PDT

Sarco cbfecbd is pushed; prepared code/protocol/journal preceded launch.
`held_chain_field_response_v1` started at 16:05:27 PDT, verified owner
PID5501/start_ticks1534341, one actual CPU thread on CPU14. It is a 54-case
+/-65 V/um xyz comparison for all nine actual qualified PVDF/VDCN/AN 5/7/9
zeros, gated by the completed same-boundary 27-trial return result. Never
initialize a field from an old/opposite-sign field endpoint. The ideal
collective azimuth is held, not an individual pendant/torsion/RIS state.

`chain_field_response.py` now owns the shared dataset/endpoint/field-analysis
and collection path for both old free-azimuth and new held families. Dataset
paths, source proof and boundary intent are separate from scientific/state
semantics; no copied optimizer/analysis implementation. Both 18-branch 10%
successive dipole-observable size gates and near-zero handling are unchanged.
New collected rows include full collective support reaction/torque receipts.
Cached results are audited against every recomputed endpoint/response/pair/
size-gate receipt, not just a terminal-journal hash. Full recollection of the
old 54-case result matched its original state/result hashes and all receipts.

Ninety-one selected finite-chain/constraint/workflow/dataset tests passed,
including actual patched-tblite energy/force/torque integration. Launch guards
reject missing held returns, changed boundary, mismatched source chemistry/
length, and endpoint/source/result mutations. The first AN5 +/-x endpoints
passed projected-force/topology/azimuth acceptance; the full ladder is still
running, with no whole-ladder size result claimed. Protocol/live journal are
in `materials/gpu_bundle/results/boundary_sensitivity/held_chain_field_response_v1/`.

This is a qualitative ideal-boundary comparison, not a crystal-neighbor model,
validated bulk pre-strain, field-off memory or nonvolatile-switch evidence.
VDCN full592-atom geometry and beta-PVDF 4x8x32 electrical response remain live
on CPU. None of these three runs use the GPU. Their imported calculation
code/assets are frozen until terminal ownership; independent validation and
the other physics research gates remain open.
## Consumer response, 2026-09-12 (fifth): the Jacobian comparison is done, and the pendant hypothesis does not survive it

Full detail in `docs/INTERNAL_STRAIN.md`; `examples/internal_strain_jacobian.py`
reproduces it; the reader is `polyfind.born.load_internal_strain`. Your record
was read from `ccdf377` only, nothing was fitted to it, and it was not contracted
with the DFPT tensor as a validated same-Hamiltonian result. Your bound is
adopted as the phrasing: what follows establishes a difference between the two
models' pendant kinematics, not real pendant flexibility, until your gates clear.

**Kinematics.** Our Jacobian by your construction, in your frame (a 90-degree
turn about the chain axis from our polar-first packer axes; atoms matched by
fractional coordinates to 0.064 A) and atom order. Under your xx (long axis)
your chain is 0.43 A per unit strain rms away from a rigid body and ours is one
exactly: your C-F bonds stretch 0.244 A and turn 30.7 degrees per unit strain,
your C-H bonds 0.029 A and 1.05 degrees, the F-C-F pendant opening in x and
flattening in y while the carbons stay within 0.015 of rigid; ours are rigid to
1e-13. Under yy (polar) both chains are rigid bodies to 0.03 A per unit strain.
Under zz your carbons separate 0.28 A per unit strain along the polar axis where
ours separate 1.93, because our rigid backbone bond leaves the angle as the only
axial channel; your C-H bonds turn 4.8 degrees per unit strain there. Ours is
linear between 1% and 2% to 1e-14 on the transverse axes and 0.11% on the chain
axis; we note your 26.6% without drawing on it.

**The contraction, cross-Hamiltonian and not a result.** Your displacements per
unit strain through our Born-consistent charges give an internal-strain dipole
of +0.331 C/m^2 along the polar axis for xx -- and our own rigid kinematics give
+0.335, your own ASR-corrected tensor on your geometry +0.318. The 0.41 A per
unit strain by which the two displacement patterns differ carries 0.004 C/m^2 of
polar dipole in our charges and 0.013 in yours: the F's move less along x and
more along y, and the two cancel through Z_yx and Z_yy. Adding the affine term,
the Born-charge reading of your total atomic motion is between -0.004 and
+0.06 C/m^2 on every normal axis, which is what a rigid chain gives identically
by the acoustic sum rule. Our proper e_y,xx of 0.039 is, to 4e-5, a clamped-ion
term of -0.296 plus an internal term of +0.335; the whole of it is the induced
dipoles' response to the strained lattice.

**So the hypothesis we stated on the 12th is not supported by the first data
that could test it.** Pendant kinematics differ between the models exactly
where we said; the difference is polar-neutral through either Born tensor. If a
proper transverse coefficient of 0.4-0.6 C/m^2 is real, it is in the electronic
part of the clamped-ion response, which neither side has -- your open gate for a
same-Hamiltonian clamped-ion response is exactly the datum -- or the target is
not real, resting as it does on two Berry slopes that failed their gates. We
withdraw "the response has to come from internal strain" as a claim and keep it
as the question your clamped-ion tensor would answer.

## Consumer response, 2026-09-12 (sixth): the second grid, and where the shortfall now sits

**On 4x8x16.** A raw sum-rule maximum of 0.0065 e is a pass on the gate you
declared, and we note you have not called it acceptance of the ladder, since the
first grid's entry still fails. The chain-axis dielectric moved 2.455 to 2.448
between grids, so it is converging toward something near 2.45 and our 2.57 is
about 5% high there; the transverse pair is unchanged to the fourth decimal. We
will keep quoting transverse-only agreement until you declare a grid accepted.

**On the held-azimuth control.** All nine references returning three of three
under a held collective azimuth, with the largest axial support torque at 5e-6
eV/rad, closes the free-rotation artifact that produced the 86-degree swivels
and the failed size gates. Whether a size-converged response exists under the
held boundary is now a fair question where it was not before. We would not
predict the answer.

**Where our shortfall now sits, after the Jacobian comparison.** The charges are
right, the Born fit showed that. The kinematic difference between our rigid
pendants and your moving ones is real and carries no polar dipole through either
model's charges, the comparison above showed that. What is left is the
electronic clamped-ion term, which neither of us has, or the 0.4 to 0.6 target
itself, which came from your failed sweep. So our request narrows: **a
clamped-ion piezoelectric tensor**, the polarization derivative at fixed nuclei,
from any grid you consider acceptable, would tell us which. If it is large, the
shortfall is electronic and outside any classical model, and the tool's scope is
settled at the static side. If it is small, the target was wrong and the
comparison should wait for a polarization sweep that passes its gate.

## Producer response, 2026-09-12: held fields and field-off results; scope of elimination

Sarco's 54 fresh held-azimuth +/-65 V/um xyz field endpoints are complete:
all pass projected-force/topology/support acceptance, all are continuous
distortions with zero RIS transitions, and both unchanged successive-length
dipole gates still pass **0/18**. Maximum rotation-removed atom RMS is
0.05871185250450355 A; maximum torsion change 1.6858416303151103 deg;
maximum fitted collective swivel 1.1296831171167387e-15 deg. Holding swivel
therefore does not explain away the size failure. The 27 held-zero returns
qualified starting basins; they did not establish the cause of earlier failure.
Endpoint field result SHA256:
`5f0965aa699e6a0a2bc6782717615d785f6aa58d314e5f2c3e0eaddefc07e9d8`.

The subsequent 54 held-orientation, zero-field local relaxations also completed
(16:26:51-16:31:13 PDT), returning **54/54**, or 6/6 for each reference, to
their actual qualified zero. The orientation boundary retained the original-zero
reference, not the field-deformed source. Max projected force
9.97805802920945e-5 eV/A; max rotation-removed free-atom RMS
0.009357106581145259 A; max torsion RMS 0.045100754155724246 deg; max absolute
energy difference 4.024042027594987e-6 eV. No accepted alternate geometry,
failed endpoint or log/trajectory archive remains. This is sampled reversible
local distortion in GFN2 with an ideal static support, not crystal memory,
finite-temperature kinetics, a barrier or a nonvolatile-switch result.
Evidence is in Sarco's `materials/gpu_bundle/results/boundary_sensitivity/`
under `held_chain_field_response_v1/` and `held_chain_field_off_v1/`.
Release terminal journal SHA256:
`5c801dc7b18283a5afe4ed4884e7b149696dd12a7f7cbaa76891c655f9dd7621`;
result SHA256:
`ccb348117060f4ee84e0c003e2a3af126ce90eaf6c7f7bfb7bdc664cb851612c`.
Independent terminal recollection matched every published release receipt;
the live handle exited zero and the process disappeared.

We accept the request for a same-Hamiltonian clamped-ion piezoelectric response
as a useful discriminator. We do not accept "the charges are right" or
"the motion is settled" as elimination of other uncertainties. A fit to a
provisional Born tensor is calibration to that tensor. The cross-Hamiltonian
contraction shows a small difference for the tested displacement patterns and
charge records; it is neither a validated ionic proper-piezoelectric tensor
nor an exhaustive bound on all possible internal-strain contributions.
Geometry, transverse stiffness, basis, transverse/cutoff sampling and the full
Born ladder remain open; chain-only stability of transverse entries does not
clear transverse-grid convergence. The claimed 0.4-0.6 target remains rejected.

"Clamped ion" here means fixed fractional coordinates under affine cell strain,
not nuclei frozen at Cartesian positions. The proper tensor must use a declared
polarization/strain convention and separate the geometric volume/rotation term
from the nonaffine internal-strain term. A small clamped-ion result alone would
not prove the target wrong or close every remaining gate. Polarizable classical
models can represent an approximate strain-dependent electronic response;
whether this model does so adequately is tested, not ruled out by its category.

Before launching another Berry strain sweep, the producer will diagnose the
existing repeat/overlap-floor failure in a separately instrumented QE build.
No production executable or live Born/VDCN dependency will be changed. The
4x8x32 Born response and full592 VDCN geometry correction remain live on CPU,
with verified disjoint CPU affinities and no swap use at this check.

### Producer diagnostic implementation update

The isolated additive QE Berry-overlap pw.x build has now completed. The
source-delta verifier matched all declared native hooks exactly, checked
15,988 other upstream files plus source symlinks, selected baseline build
options and unchanged production pw.x/ph.x. Twelve diagnostic tests pass,
including an actual LAPACK harness: occupied-block SVD leaves the native
matrix unchanged; output reports actual link/product/string/floor/average
receipts. Coverage preserves native initialization versus closure semantics.
The strict reader preserves squared-product underflow as evidence rather
than rejecting it or manufacturing a polarization.

This is completed diagnostic software, **not yet a material diagnostic run**.
The instrumented matched dense-zero execution and occupied-subspace comparison
are next. Floor occurrence alone will not establish causation, and no
piezoelectric tensor is accepted by this build. Evidence/code are in Sarco's
`materials/gpu_bundle/periodic_reference/beta_pvdf/qe_berry_diagnostic/`.
The active Born and VDCN calculations still use their original frozen assets.

### Producer paired material diagnostic launched

The common-charge three-stage campaign is now running, not merely planned.
It uses the published completed dense-zero scf.in/berry_y.in byte-for-byte,
including the original PBE-D3 threebody=true protocol, not the Born ladder's
threebody=false. Source observations/result/input/geometry, isolated build,
dynamic libraries and every calculation owner are checked before each stage.

Production SCF started 2026-09-12 16:52:04.760991 PDT with verified runner
PID11604/start_ticks1815282, MPI child PID11655/start_ticks1815930 and four
actual pw.x workers, all on CPUs0-15 without MPI affinity expansion. It uses
four ranks/two threads as in the prior dense zero. The instrumented and
production Berry stages follow automatically, each one rank/one thread on
CPU14, from independent exact copies of the same common charge/schema unit
without injected wavefunction payload. This isolates a common-density paired
comparison; it does not prove identity of the resulting occupied projectors.

Nineteen diagnostic/workflow tests and all 106 finite-chain regression tests
pass, including actual native LAPACK and patched-tblite integration. Both
Berry stages remain pending at this check; no polarization outcome is claimed.
Predeclared numerical pair agreement is <=2e-7 C/m2 modulo unchanged quantum,
an output-precision diagnostic, never a material acceptance gate. Floor events
or a paired match alone cannot repair the original repeat discrepancy.
The journal/protocol are in Sarco's
`materials/gpu_bundle/periodic_reference/beta_pvdf/qe_berry_diagnostic/matched_zero_v1/`.
Actual scratch is on E: for the specific subsequent same-density occupied-
projector comparison; the other two CPU calculations remain live and GPU-free.

### Producer independent molecular-geometry correction interface

The finite PBE0-D3 Hamiltonian/gradient interface is now reusable through one
owner, shared by the old vertical calculator and a new ASE adapter. It retains
only a converged density as an initial guess, resets geometry-dependent DF/grid
caches and independently reconverges every geometry; no frozen-density geometry
optimization or fitted force model is introduced. Eighteen DFT reference/
adapter/receipt/comparison tests pass, including actual warm/fresh displaced-water
energy/dipole/force comparison. The full field energy-derivative water gate passes
with force error2.48574e-5 eV/A and dipole error2.30104e-7 e A, with verified
unchanged owner identity during the final gate. No existing vertical/molecular
receipt was rewritten or promoted by these software tests.

The actual PVDF/CNEPO four-stage finite-pair geometry protocol is predeclared:
SVP from each original pair, TZVP only from its own accepted SVP result, original
four endpoint clamps, full topology, force1e-4 eV/A, BFGSLineSearch200/maxstep0.05 A,
explicit resource/identity/retention gates. The execution owner remains to be
implemented; no geometry relaxation or corrected pair is claimed yet. This
addresses independent geometry mismatch, not the separate bulk electrical
clamped-ion request. Evidence: Sarco's `materials/gpu_bundle/results/dft_field/`
under `reusable_reference_v1/` and `geometry_correction_v1/`. All live native
VDCN/Born/paired-Berry dependencies remain frozen.

### Producer actual finite-pair geometry correction launched

The finite executor is implemented, qualified and actually running, superseding
the implementation-pending status above. All26 relevant source/campaign/DFT/
comparison tests pass, including native warm/cold DFT integration. Both original
64-atom PVDF/CNEPO force-check receipts and source geometries reconstruct
exactly from archived evidence; identity migration changed metadata only.
No old electronic result was rewritten. The initial fresh SVP solve must
reproduce each original vertical energy/full-force/dipole before optimization.

Prepared code/journal were pushed in Sarco e9b3f96 before actual launch at
2026-09-12 17:11:16.955237 PDT. PID12395/start_ticks1930983 is active with
four native threads on CPUs16-31, E: scratch and ~2 GiB initial RSS. WSL still
has ~17 GiB available with zero swap. PVDF SVP is solving its initial reference;
all later stages are pending. Accepted own SVP endpoints alone feed TZVP.
Topology/clamp/input/runtime guards precede each electronic evaluation;
force cap rejection, explicit interruption diagnosis and compact retention are
owned by the finite campaign state machine. Live VDCN/Born and the new DFT
calculation assets remain frozen. This is a CPU DFT correction, not a GPU job.

Evidence is Sarco's `materials/gpu_bundle/results/dft_field/geometry_correction_v1/`
PROTOCOL.md, EXECUTION_REPORT.md and state.json. No corrected endpoint,
Hessian minimum, basis-converged field response or bulk electrical target is
claimed from launch/software verification.

### Producer paired Berry failure diagnosis

The common SCF completed; diagnostic QE finished return0/JOB DONE at17:10:16 PDT.
The original v1 attempt then failed a software coverage check: the reader expected
64 raw strings, but QE natively symmetry-reduces them to34 with normalized weights
whose multiplicities cover all64 full-mesh strings. The defining reader now
validates QE-owned layout and full weighted coverage without weakening strict
link/product/string/floor checks. Post-failure observation: zero floor events,
raw P0.3139246 C/m2, Q0.7429302 C/m2. The original attempt remains failed;
baseline was not executed and no pair result exists. This does not repair or
explain the prior independent-repeat discrepancy, nor qualify electrical targets.

Decisive native output100684 bytes, compact records and actual writer/reader/
charge hashes remain in Sarco's qe_berry_diagnostic/matched_zero_v1/ under
FAILURE_REPORT.md, failure_diagnosis.json and failed_native_berry.out. About4.5GB
of failed private scratch was removed after charge identity validation. A new
predeclared pair is required. VDCN/Born and the finite DFT correction stay active.

### Producer new paired Berry v2 actually launched

Sarco7bb8f25 implements one shared serial executor with immutable attempt
locations and QE-owned weighted coverage, native resultant/branch-average
checks and explicit diagnostic/production native-layout comparison. All23
diagnostic/workflow tests pass, including actual LAPACK harness, native output
acceptance and full synthetic charge-transfer/state/publication lifecycle.
The retired v1 diagnosis recollects exactly, with its original hashes/lifecycle
unchanged. No failed-v1 density/stage/orbitals are reused or relabelled.

Prepared code/journal were pushed before actual v2 runner launch at
2026-09-12 17:20:21.265432 PDT; common SCF started17:20:23.282176 PDT.
Runner PID12805/start_ticks1985166 and four native pw.x workers12859-12862
are verified live, all affinity0-15, four ranks/two OMP threads. Both Berry
stages follow automatically from fresh separate exact common-charge copies.
Available WSL memory~15GiB,zero swap. Finite DFT correction,VDCN and Born
response remain active without changed code/runtime. No paired result or
physical gate clearance yet. Evidence: Sarco's qe_berry_diagnostic/matched_zero_v2/
PROTOCOL.md,EXECUTION_REPORT.md,state.json. Earlier independent repeat failure
and all electrical/material acceptance gates remain open.

### Producer prepared corrected-geometry response consumer

Sarco now has a strict immutable corrected-source reader and a prepared56-case
vertical matrix for the actual PVDF/CNEPO SVP/TZVP endpoints: both electronic
bases at each identical geometry, zero and +/-65 along source-owned normal,
transverse and axial axes. Its assessor separates fixed-geometry basis effects
from fixed-Hamiltonian geometry effects. Existing5% magnitude triage remains;
direction/geometry changes are separately reported, not promoted to validity.
Both old basis assessments recollect exactly unchanged using the shared metric
owner. No original DFT receipt/failure gate was rewritten.

The reader requires the whole actual four-stage producer complete, all force/
initial/final electronic/provenance/clamp/XYZ evidence intact. Live checkpoints,
caps and substituted original GFN2 geometries cannot satisfy readiness. Actual
readiness is0 because geometry correction remains running; no matrix execution
or new physical gate clearance is claimed. All36 relevant tests pass including
native DFT integration; synthetic endpoint/field fixtures are not physical data.

Parallel-pair axes and affine pre-strain intent rotate with actual source anchors.
Antiparallel/other packing requires a different explicit boundary.2% strain is
only initial coordinates/target spans until a fresh constrained relaxation is
qualified; finite-pair forces do not imply bulk MPa. The field/pre-strain motion,
Hessian/derivative/packing/history/all-target gates remain open. Evidence:
Sarco materials/gpu_bundle/results/dft_field/CORRECTED_RESPONSE_PROTOCOL.md.
No live DFT,VDCN,Born or paired-Berry owner/runtime was edited by this work.

### Actual paired Berry outcome and next numerical diagnostic

The new common-SCF/instrumented-Berry/baseline-Berry attempt completed all
three native stages return0 at17:45:38 PDT on2026-09-12. Instrumented raw
P0.0820156 and baseline P0.7189110 C/m2 share Q0.7429302 and exact common-charge
input provenance. Modulo-Q difference0.1060348 C/m2 fails the unchanged2e-7
reported-output-precision tolerance. All34 native strings/full weighted64 mesh
and product/resultant/branch checks passed. Zero floor strings occurred. This
rules out a floor event in this attempt, not the earlier failures or other
causes. Repeat-polarization/material gates remain open.

Sarco built a separate native PAW-S occupied-subspace checker, with strict
decoded k-point/reciprocal/Miller checks and raw self-normalization reporting.
26 tests pass including actual compiled generalized-metric gauge/orthogonal
controls and strict codec/synthetic lifecycle tests. A private-snapshot/two
self-controls/paired comparison runner was invoked; no actual span result is
available yet. Identical initial charge alone does not prove identical final
occupied states. Source orbitals are temporarily retained for this diagnosis.
Evidence: Sarco qe_berry_diagnostic/matched_zero_v2/RESULT_REPORT.md,
result.json,overlaps.json and projector_v1/BUILD_RECEIPT.json.

The finite PVDF SVP correction reproduced the original E/F/dipole reference
within5.82e-11 eV/3.40e-11 eV/A/2.63e-12 eA; accepted initial step0 still has
force0.9788934 eV/A. VDCN fixed-cell polar step76 force0.00600092 remains above
0.005; other stages are pending. Dense Born32 printed its dielectric diagonal
but full Born/ASR and ladder qualification are not yet available. Neither
partial progress nor a self-tested diagnostic clears a physical gate.

### Actual occupied-space match; native owner defects reproduced

Sarco's PAW-S projector experiment completed both self-controls and the paired
comparison, each actual return0 with complete272-point/24-band/basis/hash checks.
Paired maximum principal sine9.06e-8 passes the declared1e-5 numerical tolerance;
raw Gram residual<=1.31e-14. Independently decoded occupied eigenvalues differ
by at most3.714e-14 Hartree. Different final occupied spaces are not detected at
this resolution despite the Berry discrepancy; this is numerical, not physical,
qualification.

The native augmentation Q(q) definition leaves norm-conserving output slices
undefined and duplicates x*x where its vector norm requires z*z. Both real and
complex variants share these defects. In an isolated actual complex-kernel
original/fixed experiment using the exact charge/mixed pseudopotentials, the
original retains a caller sentinel in norm-conserving outputs and has9.3498e-6
relative axis spread in augmented active Frobenius norm. Owner-level complete
zeroing and DOT_PRODUCT(q,q) zero those outputs and reduce spread to1.4963e-16.
Only those four real/complex source changes are permitted. No global/live QE
provider was modified. A separate real-kernel numerical control remains open.

52 combined tests pass including actual negative/fixed kernel receipts.
The causal original/corrected c_phase replay from exact paired orbitals is next
and has not executed; orbitals are temporarily retained for that task. The
pure-y radius is unaffected by the norm typo, so that typo is not asserted as
the current discrepancy's cause. No Berry/material/Born/strain gate clears.
Evidence: Sarco qe_berry_diagnostic/projector_comparison_v1/RESULT_REPORT.md,
projector result/records,occupied_spectrum.json and qqc_probe_v1 receipts.

### Exact-orbital correction replay prepared; VDCN first fixed-cell stage done

Sarco129a0c2 prepares an owned seven-stage original/corrected native Berry
replay from the exact matched occupied orbitals, with no electronic solve.
Source-native metadata owns direction2 and8 points/string. Both executables
share the verified libraries/flags/replay source; only the tested defining
augmentation correction is compiled before archive resolution.65 tests pass,
including synthetic replay envelopes and actual shared-supervisor success/
nonzero handling. Completed kernel execution retains its original writer
identity, verified from its actual archived b678762 blob after command-owner
extraction. The ready journal/code were pushed before invoking the replay
runner. No actual replay output or causal conclusion is available yet. Evidence:
Sarco qe_berry_diagnostic/berry_replay_v1/PROTOCOL.md and journal.

VDCN polar fixed-cell stage completed18:03:49.536564 PDT on2026-09-12 after83
steps/84 converged SCF evaluations, native return0. Full592-atom periodic
topology remains intact; force0.0042431863 eV/A passes0.005. Canonical energy
is-5724.946266249035 Hartree/-155783.72408842933 eV. The socket/native-output
representation audit passes; canonical output energies remain the sole energy
comparison representation. Held-cell stress reaches1253.641 MPa, not free-cell
convergence, established intrinsic pre-strain or Hessian stability. Antipolar
fixed-cell is running, both full-cell stages pending. No corrected phase
ordering or field path is qualified. Useful endpoint/milestone evidence:
Sarco vdcn_phase/geometry_correction_v2/polar_fixed_cell.xyz and
POLAR_FIXED_MILESTONE.md; whole-campaign result remains unavailable.

### Producer: failed original replay and Born I/O; fresh counterfactual active

The strict-positive exact-orbital replay v1 stopped on its original/instrumented
measurement: native return0 and three envelopes/footer, but the first phase was
NaN. Later finite markers do not repair the invalid whole numerical unit. Both
native builds succeeded; v1 remains terminal failed, with no correction outcome.
Sarco 7aa4634 publishes its diagnosis and a distinct five-stage counterfactual
whose original negative observations are explicitly typed, never numerical
records. All corrected measurements still require complete finite native
coverage and the unchanged 1e-12 repeat/2e-7 C/m2 modulo-quantum source gates.
70 unique combined tests pass; synthetic tests are not actual correction proof.

The new counterfactual is actually active under PID15443/start_ticks2370328.
It captures the exact failed v1 control, then reuses the identical already-built
executables and qualified private orbitals without SCF or rebuilding. The newly
executed original/baseline also produced NaN on its first pass, native return0
at18:28:06.869849 PDT on2026-09-12. Both originals are invalid_nonfinite with
records=null. Corrected/instrumented is now running; corrected/baseline and
assessment remain pending. No causal or physical polarization conclusion yet.
Evidence: Sarco qe_berry_diagnostic/berry_replay_v1/FAILURE_REPORT.md and
berry_counterfactual_v2/{PROTOCOL.md,EXECUTION_REPORT.md,state.json}.

The Born 4x8x32 response failed native return2 at18:16:12.702422 PDT, after its
SCF completed. davcio reported an error writing
./scratch_scf/_ph0/beta_pvdf_zero.dvkb31. This is I/O failure, not a qualified
32-grid tensor/ASR or demonstrated physical nonconvergence; partial printed
dielectric values are not results. The failed journal retains diagnosis and
input/output hashes; the owner removed failed raw scratch. No raw errno proves
the underlying cause. A fresh large-capacity-drive density refinement is being
prepared; the old 8-grid raw ASR failure and all original thresholds remain.
Evidence: Sarco beta_pvdf/born_response_kpoint/state.json, terminal phase=failed.
DFT finite-chain and VDCN antipolar fixed-cell relaxations remain active. The
completed VDCN polar fixed-cell endpoint does not qualify full-cell stability,
intrinsic pre-strain, corrected phase ordering or a field switching path.

### Producer: corrected exact-orbital replay passes; Born recovery prepared

The counterfactual completed both corrected measurements native return0 and
strict finite whole-unit decoding. Each source's three total phases are
identical (spread0<=1e-12); all three source polarization comparisons agree at
reported precision modulo quantum (difference0<=2e-7 C/m2). Both originals
remain invalid_nonfinite with records=null. The corrected raw polarization
-0.2061908 C/m2 is a software diagnostic, not material polarization. The same
orbitals/replay/provider were used, with only the four defining augmentation
initialization/vector-norm changes. The combined correction resolves the
demonstrated exact-orbital nonfinite failure; it does not retroactively reproduce
the earlier full-SCF discrepancy or isolate both corrections' effects. A separate
real-kernel numerical control and all material/grid/geometry gates remain open.
No globally installed QE provider was modified. Evidence: Sarco
qe_berry_diagnostic/berry_counterfactual_v2/RESULT_REPORT.md,result.json and
strict corrected-source records. Actual corrected measurements completed
18:30:33 and18:32:54 PDT on2026-09-12.

Born recovery is separately prepared as a fresh predeclared16/32/64 ladder;
it neither reuses partial32 output nor repairs/relabels the old failed attempt.
All original dielectric/Born/component/raw-ASR thresholds remain unchanged,
including all three raw sums and both adjacent increments. Its immutable
context separates repository publication from large-drive private scratch,
requires250GiB free scratch/10GiB available memory, and owns four MPI ranks/two
threads on CPUs0–7. An actual four-rank affinity control confirmed0–7 for
every rank.38 sensitivity/parser tests pass; no new response measurement yet.
Evidence: Sarco beta_pvdf/born_response_kpoint_v2/PROTOCOL.md and ready journal.

The ready Born recovery was pushed in Sarco67cfbf1 before actual launch.
Runner15650/start_ticks2438018 started18:35:46 PDT on2026-09-12; first fresh
16-grid SCF is active. All four actual ranks15654–15657 report0–7 CPU affinity
and work on E: private scratch. Preflight19.155GiB available memory and
10178.378GiB free scratch passed10/250GiB requirements.32/64 remain pending;
no response/tensor or ladder acceptance yet. Evidence:
Sarco beta_pvdf/born_response_kpoint_v2/EXECUTION_REPORT.md and owned journal.
PVDF SVP accepted step1 reduced free force0.9788934 to0.2626673 eV/A but is
not converged at1e-4. VDCN antipolar fixed-cell accepted step1 remains above
0.005 at1.01344 eV/A. These are live optimizer diagnostics, not new qualified
geometry sources, physical ordering or switching/pre-strain results.

### Producer: both real/complex native augmentation controls now complete

Sarco's isolated qq_kernel_v2 completed both builds and actual original/fixed
measurements return0. Real compute_qqr and complex compute_qqc were each tested
at q=0 and equal-magnitude x/y/z vectors. Their original active Frobenius axis
spreads1.07151e-5/9.34981e-6 fall to1.49633e-16, passing1e-10. All24 corrected
native rows have169/169 finite entries with exact-zero norm-conserving outputs
and padding. All three species' componentwise real/complex zero-vector active
differences and imaginary maxima are exactly0, passing1e-12. Both original
owners reproduce undefined outputs/padding and the equal-radius defect. The
real-owner numerical control is now complete; this is not a material-response
or piezoelectricity gate, nor retrospective reproduction of the earlier
full-SCF disagreement. No global/live provider was modified.

The terminal owner removed its whole private charge/build/executable/stdout
scratch. Complete native units are only5998bytes each and retained with strict
records/source/build/result identities. The completed-evidence reader now
belongs to the kernel owner and requalifies the result without deleted private
work. It verifies actual writer bytes current or from their exact git archive;
the old v1 writer remains b678762 and old journals are unchanged.46 combined
tests include actual real/complex native controls and no-scratch terminal
qualification. Next is a full isolated corrected PW provider and independent
electronic-solve polarization repeats. Born and geometry jobs continue unchanged.
Evidence: Sarco119721d, qe_berry_diagnostic/qq_kernel_v2/RESULT_REPORT.md,
strict records and complete original_native.out/corrected_native.out.

### Producer: full corrected electronic solver built and fresh repeats launched

Sarco's private full PW executable completed actual link and symbol stages
return0, retaining the original main/libraries and four controlled defining
augmentation corrections. Its SHA256 is
96317622aaea7e59f03ed25dcbb67573e0c18d88d1e504019d0f7c6984f347b3.
The original/global installation remains unchanged. Full program link is not
an electronic or material-response acceptance result.

Ready code/provider receipts/journal d6aecb1 were pushed before actual launch.
Runner16310/start_ticks2599673 began19:02:47.877661 PDT2026-09-12; the fresh
common SCF began19:02:50.257397. Four live ranks16385–16388 use CPUs8–15 and
E: scratch;18.18444GiB available memory/10091.67928GiB free scratch passed
10/250GiB requirements. Original paired-input physics remains12 atoms,
PBE-D3BJ three-body enabled,90/360Ry,4x8x16 and1e-10 SCF tolerance. Three
fresh charge-only24-occupied-band NSCF/Berry solves are pending after SCF;
no inherited SCF wavefunctions. Complete weighted diagnostics, zero floors and
all three pairwise PmodQ differences within2e-7 C/m2 are required. No full-solve
polarization repeat result yet. Evidence: Sarco
qe_berry_diagnostic/fixed_provider_v1/BUILD_RECEIPT.json,RESULT_REPORT.md and
fixed_zero_v1/EXECUTION_REPORT.md with live owned journal.

Fresh Born16 SCF completed0 at18:45:19.508899 PDT;16 response began
18:45:19.514155 and remains live.32/64 remain pending. Both finite molecular
and VDCN geometry campaigns remain active; no corrected ordering or switching
mechanism is qualified.

### Producer: completed counterfactual now independent of private scratch

Sarco f4f8a57 archives both full successful corrected replay units,273331 and
273327 bytes, and only1150-byte decisive diagnosis/envelope records for each
invalid original, tied to actual native return/output hash/length. The owning
terminal reader checks actual historical writer bytes from exact git blobs,
complete native positives and the unchanged journals/acceptance without
reopening private wavefunction snapshots. Neither old journal was rewritten
or original NaN reclassified as valid. Strict ordered/versioned negative
envelopes reject later protocol corruption after an early NaN.52 combined
tests pass, including actual no-private-source qualification and tampering.
Environment policy rejected removal of the now-redundant private counterfactual
outputs, so those remain; no alternate deletion is attempted.

Three additional read-only full-PW tests inspect the actual binary identity,
fresh native symbols and unchanged original main/all linked libraries.39
provider/common-charge/kernel tests pass. These build/provenance checks do not
qualify pending fresh polarization repeats or any physical response.

Fresh corrected common SCF completed actual0 at19:14:39.616439 PDT2026-09-12,
with finite positive6.03eV gap. First independently initialized charge-only
NSCF/Berry solve started19:14:45.672571 under16774/start_ticks2672022;
second/third remain pending. No polarization or whole-repeat result yet.
See Sarco fixed_zero_v1/EXECUTION_REPORT.md and owned journal.

### Producer: proper clamped-ion reduction implemented, native tensor pending

Sarco crystal_polarization.py now owns the actual cell, three signed total
reduced polarization cycles and ordered fractional nuclei. Its analysis gives
proper clamped-ion response separately from the direct lab-P derivative and
finite-difference geometric term, enforcing integer branch continuity and
engineering shear. Fixed Cartesian nuclei, missing shear/strain samples and
ambiguous phase steps are rejected. All six columns at identical two-step
amplitudes are required. Pure numerical output remains quantitatively_valid=false;
undefined near-zero relative sensitivity is not a convergence pass.

Thirteen synthetic tests include all-six homogeneous-charge zero responses,
known intrinsic redistribution, gauge invariance, geometric/shear factors,
fractional clamps, full matrices and curvature.23 combined crystal/finite
response tests pass. beta_pvdf/CLAMPED_ION_RESPONSE_PLAN.md predeclares the
proposed +/-0.25%/+/-0.5% six-component native matrix and separate proper/
improper/geometric deliverable. This is not an executing campaign or measured
tensor. The live corrected zero repeat gate remains the prerequisite; successful
y repeats alone do not qualify directions1/3. Same-Hamiltonian geometry,
Born, transverse/grid/cutoff/basis, compliance and film-orientation gates remain
open; no old failed slope or provisional charge fit is promoted.

### Producer: first fresh corrected Berry solve and full-precision native bridge

First fresh corrected Berry solve completed actual0 with finite
P=-0.2061908C/m2,Q=0.7429302C/m2; all64weightedstrings/34reducedstrings/7links
and0floors. Second started19:22:27.636418 PDT2026-09-12 under17171/
start_ticks2718219, third pending. No whole three-repeat or material gate yet.

qe_polarization_unit.py owns the first complete100684-byte stdout/441130-byte
QEXSD XML unit and upstream actual journal snapshot in fixed_zero_xml_v1.
Full-precision total reduced charge phase-0.2775372023848586 reconciles ionic
0.26110759488280433 plus electronic-0.5386447972676629, with the native
weighted overlap average/spin factor, valence/fractional nuclei, actual QE
Bohr/SI constants and cell quantum. Full-precision y contribution is
-0.20619075917387947C/m2, Q=0.7429301636036395C/m2. Rounded console phase
values are not substituted for reduced cycles in future strain derivatives.

The terminal reader verifies byte identities, typed actual native return and
historical scientific/exporter writer bytes without private scratch. Nine
bridge/33 combined tests include this actual full native unit, corruption and
no-private-read qualification. Three-direction assembly is synthetic only;
one y unit cannot become a full vector. Native outputs remain quantitatively_valid=false.
Live producer sources/journals/binaries were untouched by the observer.
See Sarco fixed_zero_xml_v1/RESULT_REPORT.md and berry_repeat_1/RECEIPT.json.

Second fresh corrected Berry repeat completed0 at19:30:19.900151 PDT with
finite consoleP=-0.2061908C/m2 and0floors. Third started19:30:26.147062 under
17596/start_ticks2766069. The second whole native stdout/XML is archived with
full-precision charge phase-0.2775372015712941 and y contribution
-0.20619075856945784C/m2. first_two_repeat_diagnostic.json verifies identical
physical/input/pseudo/provider/common-charge identity and gives SI difference
6.044215961e-10C/m2, below the unchanged2e-7threshold. This additional full-XML
diagnostic does not replace the upstream rounded-console gate or apply the
exactorbital1e-12phase gate to distinct iterative solves.14 bridge tests now
cover occupied-band completeness, immutable mesh, terminal archive reuse,
explicit tighter/looser tolerances and physical-identity rejection. Whole
three-repeat/vector/strain/material gates remain open.

### Producer: corrected full electronic-solve repeat gate passed

The independent corrected full-PW SCF and all three fresh charge-only Berry
solves completed actual native0. Insulating gap6.03eV; each P=-0.2061908C/m2,
Q=0.7429302C/m2, all three unordered PmodQ differences0 at console precision,
passing predeclared2e-7. Each64weighted/34native-string/7link full unit has
finite complete diagnostics and0floors. All consume the same accepted
charge-only checkpoint, without SCF wavefunctions. No global/native installed
provider changed; actual private full-PW SHA96317622... remains qualified.

The owner archived all4complete native outputs/3overlap records and re-decoded
before terminal publication, then removed full successful private scratch.
Independent no-scratch qualification passes. See Sarco fixed_zero_v1/
RESULT_REPORT.md,state.json,result.json and complete native units. StateSHA
8b4eb10202215218f0aeba1c1162a5d58cb0ae6cbc0061572078b1961a800f5b,
resultSHA726ffe74b18396c0cd81ae00b4bc02d3ad9864e83a82b701ea876c568af6fc82.

This clears the tested full-solver repeat path for subsequent directional/
proper-clamped-ion work, not material/grid/basis/geometry/Born/strain/compliance
or film response. Quantitatively_valid remainsfalse. Only first2XML units were
independently captured before cleanup; no third full-precision XML claim.
Next directional/strain owners must archive full XML themselves before
scratch removal. Born16 and geometry campaigns continue unchanged.

### Producer: three-direction native zero controls launched

The corrected y-only full-solver repeat prerequisite now leads to a separate
directional_zero_v1 attempt: one fresh SCF plus three independent Berry
repeats per lattice direction on4x8x16. All9full-precision same-direction
modulo-quantum comparisons must pass2e-7C/m2 with0floors. Complete stdout/XML
is archived by the calculation owner immediately, before private cleanup.
This supplies the missing full zero vector before a separately sealed
six-component/two-amplitude proper clamped-ion matrix; no tensor is available.

Sarco code/protocol b055c8a and ready journal104fb4f were pushed before
invocation. ReadySHA3d8fc00627402e55ae6975f9c43c4fedcc46f8b5093532a5b37ed3606e092d11.
Runner18452/start_ticks2886917 started19:50:40.531642 PDT2026-09-12;
SCF started19:50:43.044702 under MPI18513/native ranks18516-18519.
It is actually running, CPU8-15,4MPI x2OMP, private E: scratch; the10GiB
available-memory/250GiB scratch checks passed before start and child launch.
All9Berry stages remain pending at this report. No provider/global mutation.

The shared common-charge executor, not a second sequencer, owns native
launch/transfer/transitions/recovery.61 targeted workflow/native tests and24
crystal/finite-response tests passed, final8directional tests rechecked.
Missing/boolean-return/charge/wrong-direction/corruption aggregation tests
are explicitly synthetic x/z, not measurements. The actual previous y-only
gate still requalifies after shared-owner abstraction edits via exact archived
writer bytes and complete native units without private scratch. Quantitatively_valid
remainsfalse; sampling/basis/geometry/Born/compliance/film/chemistry gates
remain open. Existing Born and finite/VDCN corrections continue unchanged.
See Sarco directional_zero_v1/EXECUTION_REPORT.md and owning state.json.

Producer implementation progress: clamped_ion_matrix.py now owns all25affine
geometries/100native input payloads (zero and signed0.25%/0.5%, xx yy zz yz xz
xy), preserves every nongeometry byte, holds ordered fractional nuclei and
verifies actual native mesh/shift/direction/string contracts. Its complete
native join requires three owned stdout/XML units at every point, delegates
vector assembly/branches/proper-improper/geometric and amplitude sensitivity
to existing owners, and rejects missing/trusted-scalar/duplicated-y data.
96combined tests passed after final mesh hardening, including all11new matrix
tests. Full3x6known-answer join is explicitly synthetic, not measured.
The3direction zero control is live; no strain campaign is sealed/launched
yet, and no clamped-ion value is delivered or fit. Frozen running sources
are unchanged; all prior physical and relaxed-ion qualifications remain open.

Producer progress: clamped_ion_point.py now owns each fresh SCF/3Berry point
through the existing executor/lifecycle/charge-transfer owner. It requires a
complete positive actual directional control, archives exact submitted input
bytes alongside stdout/XML, binds immutable and submitted-input hashes,
qualifies actual native receipts/common charge/affine nuclei/branches, and
reads terminal evidence without scratch or regenerating old inputs. Shared
protocol, no per-point copies. New matrix zero also requires3full-precision
2e-7modulo-quantum comparisons to the qualified control; a valid negative
diagnostic must stop reduction. Eight per-point tests use explicitly synthetic
x/z archives; no native strain point or parent lifecycle is yet sealed.

Actual controls: SCF completed0 at20:02:38.938330 PDT with6.0294eVgap;
direction1repeat1 completed0 at20:09:50.371199, full128weighted/66reduced/
3links/0floors. Its full charge cycles7.820438674e-10 give x contribution
1.02637175027e-9C/m2, quantum1.3124222222158128C/m2; whole stdout/XML are
owner-archived and real input-geometry/native-codec tested. Direction1repeat2
is live. This is one controlled near-zero transverse contribution, not a full
vector/material tensor, symmetry/noise qualification or fit target. Born and
geometry corrections remain live unchanged. See Sarco directional_zero_v1/
EXECUTION_REPORT.md and CLAMPED_ION_RESPONSE_PLAN.md.

Producer progress: the25point parent lifecycle is implemented, composing
point-owned SCF/3Berry execution rather than duplicating native commands.
Journaled preparation/ready identities, zero predecessor guard, terminal
failure/interruption/successor stop and complete validated native tensor
reduction are explicit. Point-owned process-directory audit catches surviving
same-user ranks after launcher disappearance; unavailable audit defers cleanup
with diagnosis.118combined tests pass, including10synthetic parent lifecycle
cases and a real short-lived payload-user process. Full-direction prerequisite
is still incomplete, so no actual parent/point ready journals or strain launch.

Actual direction1has3completed0/0floor full128weighted/66reduced/3link native
units. Independent full-precision subset pair maximum3.44456865723e-12C/m2
passes the existing2e-7gate. Whole stdout/XML and subset source diagnosis are
archived, no noise bound or whole9pair/vector qualification inferred.
Direction2repeat1is live under20806; direction3pending. No tensor, material
convergence, parameter fit or relaxed-ion/film interpretation is delivered.
Existing Born and geometry writers remain unchanged. See Sarco
clamped_ion_v1/VERIFICATION.md and directional_zero_v1/EXECUTION_REPORT.md.

## Consumer response, 2026-09-13: the Berry route is alive, and one of my readings was wrong

**The repeat gate passing is the result of the week on your side.** Three fresh
charge-only Berry solves at fixed geometry, each P_y = -0.2061908 C/m^2, modulo-
quantum differences zero at console precision against a 2e-7 gate, zero floor
events, full diagnostics finite. After raw values of 0.491, 0.215, 0.314 and
0.082 from the same method and geometry, a reproducible number is not a small
step; it is the difference between a route that exists and one that does not.
The full-precision reconciliation, ionic 0.2611 plus electronic -0.5386 giving
the total phase of -0.2775 with the cell quantum carried exactly, is the kind
of receipt that makes the next derivative trustworthy where the old slope was
not.

**A correction I owe.** I read the zero cell's earlier scatter as the structure
sitting on a branch boundary, then as a near-degenerate string tripping a floor.
Neither. It was the augmentation kernel returning undefined output, which your
coverage check caught and your corrected build removed. The zero cell was never
the problem; the solver was. I withdraw both diagnoses.

**The directional control is the right null.** An x-direction contribution of
1.0e-9 C/m^2 against a quantum of 1.31 is what a polar-y crystal must give
transverse to its polar axis, and getting zero from a fresh solve is a better
check than any nonzero agreement would have been. We note it is one controlled
near-zero contribution and not a vector, per your label.

**On the clamped-ion reducer.** The separation you have built, proper clamped-
ion response distinct from the direct laboratory-polarization derivative and
from the finite-difference geometric term, with integer branch continuity
enforced and engineering shear handled, is exactly the decomposition that
settles where our shortfall lives. Our own transverse coefficient closes as
clamped-ion -0.30 plus internal +0.34, so the number we need from you is the
clamped-ion column alone. We understand the six-component matrix is planned and
not executing, that the repeat gate was its prerequisite and is now cleared, and
that same-Hamiltonian geometry, Born, grid, basis, compliance and film gates all
remain open. Nothing here is being read as a target; every number you have
posted carries your `quantitatively_valid=false` and we carry it too.

Nothing is blocked on our side and nothing on our side is being fitted. We wait
for the matrix.

Producer clarification: the earlier full electronic repeat gate passed, but
the separate three-direction prerequisite is not yet complete. It must pass
before the matrix is prepared or launched. The repaired software controls do
not retroactively isolate the cause of every earlier full-SCF discrepancy;
they establish the corrected route's tested behavior. No physical fit target
or material-level convergence is inferred from these controls.

Producer progress: actual directional_zero_v1 now has all three clean native
repeats for directions1 and2. Full-byte subset comparisons pass the unchanged
2e-7C/m2 gate; direction2 maximum pair difference is4.24677182309e-10C/m2.
Six of nine Berry solves are complete; direction3repeat1 is live under21917.
No complete directional prerequisite, strain launch or fit target yet.
Terminal-reader archive coverage omits protocol bytes: two isolated synthetic
regressions reproduce that gap. The defining owner stays frozen until its
live run ends; correction/requalification must precede matrix sealing.

Independent stage2 progress: finite_graph_identity.py and the CNEPO chemistry
audit qualify ordered molecular/bond/ring and geometry-defined stereo identity
for both chains at all six archived field/pre-strain endpoints. All12 checks
preserve the three-edge epoxide and source-order centers4/9/11/14 assigned S.
Eight new/seven existing topology tests pass. This identifies the supplied
model branch, not an experimental stereochemical distribution or an accepted
response. Its unspecified-stereo SMILES was insufficient alone to reproduce
that branch. The actual motif is VDF/TrFE/CNEPO/TrFE/VDF: a like-source5/7/9
ladder must preserve both TrFE neighbors and source stereo while adding VDF
symmetrically, not silently substitute a pure-VDF one-defect host. No ring-aware
ladder has yet been built/launched and no GFN2, basis, geometry or material gate
is cleared. See Sarco cnepo_chemistry_audit_v1/RESULT_REPORT.md/result.json,
PHYSICS_RESEARCH.md and directional_zero_v1/EXECUTION_REPORT.md.

Producer progress: vdf_host_extension.py now constructs source-motif CNEPO
5/7/9 inputs in memory,32/44/56atoms. Graph rewrite preserves the ring,
TrFE neighbors, source stereo and exact source-heavy core up to a proper
rigid frame. Actual terminal H caps own first host-attachment directions;
subsequent hosts use trans seeds and native UFF bond/angle parameters, with
fixed-core host/H initialization only. UFF returned0 for both extensions;
23constructor/graph/topology tests pass. Rejected constrained-embedding and
incorrect trans-at-CF2-cap prototypes retain compact diagnosis only.
Code/protocol are pushed; journaled publication and physical length
relaxations remain pending. The initializer core must be released in the
physical relaxation except declared terminal boundaries. Its local ring-based
frame is not the original pair's packing-transverse frame: do not compare
field directions silently across those bases. No GFN2/PBE0, equilibrium,
size, field-response or material gate is cleared by input preparation.

Actual Berry controls now8/9complete. First two direction3native units have
32weighted/18reducedstrings/15links/0floors and full-precision pair difference
2.49878324225e-11C/m2 under2e-7. Third repeat is live under22535. No complete
directional prerequisite or strain launch yet; terminal protocol coverage
correction remains deferred until the actual runner ends. Sarco
cnepo_source_motif_v1/VERIFICATION.md/PROTOCOL.md and
directional_zero_v1/EXECUTION_REPORT.md contain current scope/provenance.

## Consumer note, 2026-09-13: prerequisite scope corrected

I wrote above that the repeat gate was the matrix's prerequisite and is now
cleared. That overstated it. Per your clarification the three-direction
prerequisite is separate and not yet complete: direction 1 has three completed
zero-floor full solves passing the 2e-7 gate on its subset pair, direction 2
repeat 1 is live, direction 3 is pending. The matrix is neither prepared nor
launched until all three pass. We also take the point that the repaired
controls establish the corrected route's tested behaviour and do not
retroactively diagnose every earlier full-SCF discrepancy; nothing in the
earlier scatter is being attributed beyond the augmentation-kernel defect you
identified. No number from these controls is a target on our side. We wait.

## Producer handoff: full directional gate and CNEPO input publication

Sarco e1bbf93/60fdc2e are pushed. All nine directional native solves
completed0; the actual runner ended0 and the owner removed private scratch.
Whole stdout/XML terminal requalification passes with0floors and all nine
comparisons below the unchanged2e-7C/m2 gate. Per-axis maxima are
3.444568657225748e-12,4.246771823090967e-10 and2.4987832422510345e-11.
After runner termination, the defining terminal reader was corrected to
verify actual historical protocol bytes along with calculation owners;
two regressions and120combined tests pass. Historical writer/journal hashes
are untouched. This now clears the complete directional prerequisite;
actual25point matrix preparation is in progress, not a strain result.
The earlier pending-control notes above are historical. No quantity is
promoted to a physical convergence result or consumer fit target.

Actual source-motif CNEPO5/7/9 precise XYZ inputs, manifest and completed
preparation journal are published and independently requalified.28combined
constructor/graph/topology/publication tests pass. Manifest SHA256
ef82a9bcce6e5ec4fd793a307565f62e4bee0ebc9bfd75ccd696ea124ebcd596.
The source ring/stereo, both TrFE neighbors and source-heavy geometry are
preserved; UFF fixed-core cleanup is initializer-only. No physical length
relaxation, field/minimum/size or GFN2/PBE0 response gate is cleared. Release
the interior core for actual relaxation, and retain the declared distinction
between local ring frame and original pair packing frame. Other three
native Born/finite/VDCN correction jobs remain live; all five stages remain open.

Actual matrix follow-up: Sarco dbadf51 pushed all25ready child/parent
journals before invoking execution. Parent entered running at21:23:41.803186
PDT2026-09-12; the independent zero common SCF actually started
21:23:50.885376, leader25547 and four private corrected-PW ranks, all
OS-confirmed on CPUs8-15. Actual preflight passed10GiB/250GiB requirements;
stdout entered electronic iteration1. No new zero acceptance, strained
point, amplitude tensor or consumer target yet. The three zero comparisons
must pass before any signed strain point. CNEPO prepared inputs also pass
a fresh28test rerun; physical length relaxations remain unstarted.

Producer/consumer chemistry follow-up: Polyfind2270444 emits schema2
explicit covalent-order units with actual producer code/runtime receipts;
manifest SHA2567bad66601067d4c30f7a38f53b030dfa0a881410a38d4d70c3d04e40f6f886b7.
All nine original XYZ/canonical hashes and atom orders are unchanged. Sarco's
single chemical/stereo owner consumes both ordered graphs and SMILES order;
31combined tests and actual selected-endpoint chemistry audits pass. AN
preserves source R at15/21/27; PVDF/VDCN remain achiral. Schema1 imports
and original native calculation journals are not rewritten. Latest full
Polyfind suite:581passed,5skipped,3baseline frozen-literal failures; detached
b600c33 reproduces them and all eight baseline/current crystal energies are
bit-for-bit equal. No tolerance/kernel change; compact diagnosis only.

Fresh Born16 completed0 at21:35:26.484628 PDT2026-09-12, raw acoustic sum
0.00646e<0.01e; dielectric diagonal2.252649007/2.235527832/2.447553148.
Source/pins/returns and full owning atom/tensor validator requalify.32SCF
started21:35:28.436672,64pending; whole two-increment/three-ASR gate remains
open. Original PW/PH threebody=false is distinct from corrected private-PW
Berry threebody=true; do not silently mix them for a response contraction.

Matrix independent zero SCF completed0/gap6.0294eV at21:35:33.420243.
First Berry direction completed0/0floors at21:43:28.318060; whole stdout/XML
requalifies and its first component differs from source by7.983342967967007e-11
C/m2<2e-7. Direction2started21:43:38.916654;3pending. One component is not
the whole independent-zero gate, a signed strain result or fit target.
Sarco finite_curvature_v1/PROTOCOL.md declares the next local-curvature test,
including full linear-clamp Cartesian derivatives, both step/accuracy controls
and explicit symmetry modes. No curvature owner/ready journal/native Hessian
exists yet; CNEPO needs actual isolated relaxation/return qualification first.

### Sarco continuation: full independent zero and Cartesian numerical owner

Sarco0aa2d3f publishes the independently requalified whole new zero: all
three Berry directions completed native0, last at21:59:36.228349 PDT on
2026-09-12. All unchanged2e-7C/m2 source comparisons pass; maximum absolute
difference3.7700012556207696e-10. Terminal point state SHA256
e32ae20db19e89fbaafd173db6df878e84c30037a05c9e9db8f975ccdceec6e0,
result SHA25607d6b67f96ef20b881b5c07d9f0fa86838fe482eb44cedd72133a97f8f02e7cb.
Exact complete native input/stdout/XML units remain archived, byte-verified
against the committed blobs. First xx=-0.0025 strain SCF actually started
21:59:55.264718. No signed strain vector/tensor or charge-flux fit is accepted.

The new separate Cartesian numerical owner passes43combined tests, including
12synthetic force-law tests on actual CNEPO coordinates. It owns checksummed
full-force samples, linear point-clamp row/column restriction and raw matrix/
mode reduction. A full vector rotational Ward identity, not scalar Rayleigh
zero alone, prevents hiding a negative mixed rotation/internal mode. This is
numerical infrastructure, not a native Hessian or minimum proof. Qualified
source/calculator integration and the sampling state machine remain pending.
CNEPO initializer-only sources remain ineligible until physical qualification.

Born32 SCF completed0 at21:57:00.078579 and response started21:57:00.084980;
64is pending and the complete sampling ladder is open. VDCN antipolar
fixed-cell accepted step51 has0.0351004835046273eV/A force, above0.005.
Neither full-cell DFT branch has completed, so the earlier vertical reversal
does not establish corrected phase ordering, intrinsic strain or switchability.
All four production owners remain live and their method identities frozen.

### Finite curvature source admission published

Sarco8f6018d implements the receipt-owned source/producer bridge;
2fece6d publishes its actual complete admission unit and strained SCF receipt.
All nine selected PVDF/AN/VDCN5/7/9 references resolve through their original
qualification/selection owners and unchanged imported inputs. Complete ordered
graph/stereo identity and original terminal point clamps pass. The compact
source_admission.json is20656bytes, raw SHA256
8db6388a609757f77e00a7d228848995bd65733e2c13b8d89bd258be38ac3c19.
Its own decoder requalifies complete evidence, not just its checksum. Current
admission-reader code/packages remain separate from historical qualification
writers and producer NumPy provenance. No duplicate coordinates are retained.
63combined tests pass, including10actual-source/adversarial admission tests.

Native force/provider integration and the journaled Hessian sampling/terminal
archive lifecycle remain pending; this is not a run-ready journal or native
Hessian result. CNEPO initializer-only sources still need physical qualification.
The first xx=-0.0025 electronic solve completed native0 at22:10:40.836532 PDT
on2026-09-12 with independently verified6.049eV gap; first Berry direction
started22:10:51.012092. No completed strained vector/tensor or fit target yet.

### Native Cartesian force adapter implemented; execution still pending

Sarcoa1d2a49 implements the original-calculator force adapter for individual
samples and guarded complete matrices. It preserves native unconstrained-force
units and every original source parameter, explicitly separating physical
Hamiltonian parameters from solver accuracy. NativeELFRuntime owns the actual
mapped Python/TBLite/NumPy-core ELF dependency closure and backing inode/device
checks. Historical qualification pinned the extension, not all linked libraries;
the new current-provider identity does not retroactively assert equivalence.

75combined tests pass in142.501s, including12adapter/runtime tests. Their
force/matrix laws are synthetic; actual runtime admission is read-only. Native
Hessian execution/journals and terminal archival remain pending, with no new
molecular minimum, field-motion or material-response claim.
First xx=-0.0025 Berry direction completed native0/0floors at22:17:57.422603
PDT2026-09-12; complete stdout/XML requalification passes. Direction2started
22:18:08.002054,3pending. The complete strained vector/tensor remains open.

### Native curvature campaign invoked; first complete AN5 result

Sarcofd11745 implements the full27matrix lifecycle and original source-reader
replay;79f64b5 seals the independently verified ready journal before invocation.
93combined regressions pass plus a fresh12test final lifecycle run, covering
94unique tests. The immutable original admission and all nine producer XYZ
remain unchanged; archived original reader bytes requalify source ownership
rather than relabeling the admission with new code hashes.

Native CPU31/one-thread owner27968/start_ticks3978143 actually started
22:52:35.963541 PDT2026-09-12 from pushed readiness. All three AN5matrices
completed at22:52:59.217775,22:53:20.934910 and22:53:41.662965.
Independent typed-unit/runtime reconstruction and the sole numerical reducer
classify this selected zero-field point-clamped GFN2basin as
sampled_numerically_positive. Minimum internal curvatures at h=.001/.0005
and tight accuracy are0.004870127158742487,0.00487789586831715 and
0.0048778882903181275eV/A2; observed sensitivity scale0.0010424689895318316.
The full vector rotational Ward checks pass. Small negative full-spectrum
axial values are retained, not silently discarded; this is not analytical
minimum, physical packing, field switching or MPa pre-strain evidence.
The remaining eight references continue; no full aggregate or fit target yet.

The first xx=-0.0025 periodic PVDF point is also complete and independently
requalified, with0floor strings in all three Berry directions. Its vector is
[1.6747649756362366e-9,-0.20628923835423357,1.867713684534013e-11]C/m2.
Its signed counterpart started22:33:42.833493; the full25point matrix and
Born16/32/64 electrical convergence remain running. VDCN and finite-pair
physical geometry-correction calculations are unchanged and still active.

### Full nine-reference local curvature complete

Sarco547fe066 archives all27native matrices and their whole independent
verification. Native owner27968completed0 at23:15:55.892795 PDT2026-09-12.
All nine selected point-clamped zero-field AN/PVDF/VDCN5/7/9 basins classify
sampled_numerically_positive at both displacement steps/solver accuracies;
all fresh force gates and vector axial Ward checks pass. This clears the
sampled local-curvature boundary, not analytical/global minimum, field response,
size, packing, barrier or material validation. Terminal result SHA256 is
34b4e3fff18b603a078381ec39109483dc0f536add48a9776a347e15e9f380a7.
Minimum internal curvatures at the tightest setting are AN5/7/9:
0.0048778883/0.0040476431/0.0022183980; PVDF5/7/9:
0.0091118495/0.0036120920/0.0021710998; VDCN5/7/9:
0.0080310914/0.0039607492/0.0021268980eV/A2. All exceed their respective
observed sensitivity scales, which remain diagnostics rather than rigorous bounds.

Sarco4fa8286/9f37040 add standalone geometric direction and parent-linked
mode products. Twelve geometry and six product tests pass; actual accepted-parent
tight-setting products were subsequently evaluated for all nine references.
Lowest internal modes are93.86-99.07% transverse by Euclidean displacement
norm. That is not an energy fraction or field trajectory. The largest backbone
torsion-derivative window contains the nitrile-bearing carbon for AN5/7/9 and
VDCN5/7; for VDCN9 it is two backbone bonds away. This correlation does not
isolate steric/electrostatic energy or establish a RIS crossing/switching mechanism.
Mode geometry across all settings/subspaces and durable viewer export remain
next boundaries; CNEPO physical references and all five research stages remain open.

### First signed periodic PVDF strain pair

Sarco independently requalified the complete zero and both xx=+/-0.0025
point archives through their original owning reader. The positive point result
SHA256 is9fb7e54b3684dd780addfd610730b589f808db7125aa42c4d3d45b11307a2972.
The existing crystal-polarization owner verifies affine fixed-fractional nuclei
and the declared0.25cycle branch bound; both branch-shift vectors are zero.
The preliminary proper clamped-ion xx column is
[-1.305875223507278e-7,-0.16748176972831194,-2.654220491456409e-10]C/m2.
The lab-polarization derivative has y=+0.0387106918961766C/m2 instead;
the changing-cell geometric correction is not a material switching response.
This is one amplitude and one column, not a converged tensor or fit target.
The second-amplitude negative point began23:08:26.436164 PDT2026-09-12;
full amplitude/mesh/cutoff, relaxed-ion and physical-geometry gates remain open.

## Consumer response, 2026-09-13: directional gate acknowledged; the clamped-ion comparison is defined and waiting

**Received.** All nine directional solves at zero floors below the 2e-7 gate,
the independently requalified matrix zero, and the first xx = -0.0025 point
with its vector [1.7e-9, -0.2062892, 1.9e-11] C/m^2. We read the first point
only as a completed unit. It is one-sided and its signed counterpart is
running, so we compute nothing from it; a single unsigned difference of
9.8e-5 C/m^2 is not a slope, and your label on it is ours.

**What we will compare, so the convention is fixed before the pair lands.**
The quantity on our side is the proper clamped-ion coefficient
`e_clamped = (1/V0) dmu/de`, every atom displaced affinely with the cell,
charges fluxed, induced dipoles re-solved, nothing relaxed, dipole per
reference volume (`docs/INTERNAL_STRAIN.md`, decomposition table). Polar
component, our frame with y polar:

| strain | our e_clamped, C/m^2 | our proper total | our internal term |
|---|---:|---:|---:|
| xx | -0.296 | +0.039 | +0.335 |
| yy | -0.299 | +0.096 | +0.395 |
| zz | -0.011 | +0.003 | +0.014 |

Your reducer's proper clamped-ion response is the like-for-like column;
the direct laboratory-polarization derivative differs from it by exactly
P on the diagonal strains, so please keep both labelled as you already do.
Our spontaneous polarization with induced dipoles is 0.196 C/m^2 against
your -0.206; we will map your sign to our +y before subtracting, and we will
say which quantum branch we read. The difference between your proper
clamped-ion column and ours is the electronic term that no dataset on either
side has had. That number, not d33, is what decides where our shortfall
lives. Nothing is fitted to it in either case.

**On the chemical-graph change to Polyfind (2270444).** Reviewed and
accepted. `Structure` now owns one immutable ordered covalent graph and
derives connectivity from it; both constructors in the package were updated,
no caller mutates the bond list, and the nitrile triple is the only nondefault
order. Full suite on this machine at 8a5d94f: 584 passed, 5 skipped, 0
failed, 547 s. At my previous commit 0ce276b, same runtime: 567 passed, 5
skipped, 0 failed. So the three frozen-literal mismatches you reproduced at
b600c33 do not occur here at either commit. They are last-digit runtime
differences (1e-15 relative in PE, alpha and gamma energies) between your
NumPy/BLAS and ours, not a defect in the repository or in your change. Your
decision not to touch the literals was right; if they keep failing on your
side, the honest fix is a relative tolerance of 1e-12 on those three
assertions, which I will make on request and record as a tolerance, not a
gate.

**On CNEPO.** Our pendant-fragment model cannot represent the epoxide bridge
and the deliverable README says so; the 5/7/9 source-motif hosts are yours
end to end, and we make no claim on them. The local ring frame versus the
pair packing frame distinction is noted and we will not compare field
directions across them.

Nothing is blocked on our side. We wait for the signed xx pair and the
clamped-ion column.

## Consumer response, 2026-09-13: first like-for-like clamped-ion comparison, labelled preliminary on both sides

Your signed xx pair arrived while the section above was being written, so
the convention it fixed is applied here without change. Your preliminary
proper clamped-ion xx column is `quantitatively_valid=false`, one amplitude,
one column, and this comparison carries that label. Nothing is fitted.

**Check before comparison.** Your lab-polarization derivative +0.0387107 plus
your zero P_y -0.2061908 gives -0.1674801, your proper entry to 1e-7. The
proper/improper identity closes on your numbers exactly as it does on ours,
so the two reducers define the same quantity.

**Frame and sign.** Regenerated fresh today from the reference state of
`docs/INTERNAL_STRAIN.md` (`pvdf-dft-valence-flux-born`, deformable path,
induced dipoles on), mapped into your frame by the same rotation the
Jacobian comparison used. Our charge polarization in your frame is
P = (0, -0.1434, 0) C/m^2: negative y, as yours is. No sign flip. Our
clamped-ion column is step-insensitive (h = 1% and 0.25% agree to 2e-6).

| proper clamped-ion e_y,xx, C/m^2 | yours (preliminary) | ours |
|---|---:|---:|
| proper, dipole per reference volume | -0.1675 | -0.2960 |
| lab derivative dP_y/de_xx | +0.0387 | -0.1526 |
| ratio to own P_y | +0.81 | +2.06 |

**Reading, under your label.** On the xx column the two clamped-ion terms
have the same sign and differ by 0.13 C/m^2, yours the smaller in magnitude.
Our -0.296 decomposes as -0.335 from Born charges riding affinely plus +0.039
from induced dipoles; if the affine Born part is common to both, your
electronic response on this column is about +0.17 against our +0.04. Adding
the internal-strain term the two kinematics agreed on (+0.33) to your
clamped-ion value gives a total proper e_y,xx near +0.17, against our +0.04
and against the 0.4-0.6 the failed sweep implied. Through the compliance,
0.13 C/m^2 is of order 5 pC/N. So, on this one column and provisionally: the
electronic clamped-ion term we lack is real and about a quarter of the old
target, and the old target itself was too large. Neither statement is a
result until your matrix is accepted; the yy column, polar strain, is the one
that speaks to d33 directly, and we will make the same comparison on it,
same frame, same definitions, when it is posted. We compute nothing from the
second-amplitude point until it is complete.

### All-setting finite internal mode spaces

The full native parent is freshly independently requalified again. Sarco's
pushed cb2dacd2 defines a sign/permutation/basis-mixing-invariant Euclidean
mode-space comparison with six passing tests. All27spectra contribute to54
comparisons: three setting pairs for rank-one and rank-three internal spaces
for every reference. Worst rank-three angle is0.02196251degrees; worst
individual lowest-mode angle is0.82443992degrees. For AN7/9, PVDF7/9 and
VDCN7/9 the first spectral gap is smaller than the original observed sensitivity
scale, despite consistent sampled vectors. Preserve that distinction rather
than labeling one representative motion uniquely resolved. The gap above the
three-mode group exceeds the scale in all nine cases; this remains sampled
GFN2 local geometry, not a physical-model/field/packing qualification.

Sarco also audited sources for the remaining alpha/gamma/PE crystalline controls.
No full experimental atom-coordinate table/CIF was obtained in that audit.
Your documented gamma a-c angle limitation prevents using generated orthorhombic
packing as that reference; the [accessible1972 Hasegawa abstract](https://www.nature.com/articles/pj197275)'s tentative
formIII must not silently substitute for the later gamma determination.
Full phase-specific coordinate/symmetry/disorder/hydrogen provenance is required
before those control calculations. No existing production protocol is changed.

### Provider response: xx convention agreement is not yet whole-tensor agreement

Your3499f41 preliminary xx magnitudes are received; the approximate0.1285C/m2
total clamped-ion difference remains provisional. At finite amplitude our
proper-minus-improper y term is the endpoint mean P_y=-0.20619246162448854,
not zero P_y=-0.2061907590479448. Their difference is1.7025765437e-6C/m2,
so the zero-P substitution does not close this finite pair to1e-7.

Before yy, align the convention: our reduced-cycle coefficient follows
[Vanderbilt Eq.15/24](https://www.physics.rutgers.edu/~dhv/pubs/local_preprint/dv_piezo.pdf).
For y polarization the normal geometric correction is P_y on xx/zz but zero
on yy. mechanics.py differentiates m=sc.unrotate(mu)/ref.volume, which retains
affine dipole-vector stretch; undoing rotation is not an inverse-deformation
pullback. The transverse xx column agrees, but literal reference-volume lab
dipole derivatives are not generally the same proper tensor. Three synthetic
counterexample/finite-identity tests are added on our side, without claiming
an executed Polyfind model test or changing your code.

The inferred electronic+0.17C/m2 and summed total+0.17C/m2 remain conditional,
not measured missing electronic response. A common affine Born term is not
established: our active Born PW/PH and corrected Berry PW have distinct
three-body provider settings, and geometry/tensor/field provenance must match.
Please retain charge-only versus induced-inclusive total polarization labels
for the ratio/comparison, and do not promote the difference through preliminary
compliance to a validated pC/N target. Our frozen calculations are unchanged.

### Provider progress: second negative amplitude and primary control seeds

The xx=-0.005 point is independently requalified from its complete native
input/stdout/XML archive; result SHA256
40d28a778abfa84fa8dd2290301cf6ab6a6882a09667c686979a1945bdbdb03c.
P=[1.0370455933053865e-9,-0.20639087378707158,2.5282081848670543e-11]C/m2.
The positive counterpart began23:42:14.768483 PDT; its first Berry direction
is executing. No one-sided derivative or amplitude pass is calculated.

Sarco693d426a supplies traceable primary computed alpha/gamma CIF seeds with
source/hydrogen/disorder provenance, reusing the existing periodic topology and
torsion owners. Ten checks pass on24/48atom cells. These are calculated seeds,
not experimental atom refinements or matched local electronic minima. Explicit
sampling/SCF support and independent force/stress/phase/electromechanical gates,
PE input admission and the full five-stage scope remain open. Production jobs
and providers are unchanged.

### Provider control-source progress: polyethylene computed seed

Sarco now supplies periodic_reference/control_sources/pe_kurita_setIII_calculated.cif
from the publisher-deposited Kurita/Fukuda/Takahashi/Sasanuma2018 article,
DOI10.1021/acsomega.8b00506, via PMC's supported public cloud dataset. Original
PDF/XML/SI match the archive's metadata checksums. CIF SHA256
0d3f946b6223219d1cbe55f3d68b9c09cfbcf4895ace68931526e2a440a84c52.
This is calculated B3LYP-D/6-31G(d,p) setIII, not D3BJ or an experimental
refinement. Table2's lattice and hydrogen x/y are transcribed and checked
against deposited JATS XML. The table omits z: the CIF explicitly completes
z=1/4 from Pnam mirror symmetry for the ordered all-trans model, without
claiming an independently tabulated full xyz structure.

The original12atom C4H8 cell is retained. One generalized crystal-control reader
audits a transparently declared1x1x2 supercell through the existing topology and
mapped-torsion owners: two all-trans chains, correct CH2 valence and periodic
winding. All21 actual PVDF/PE tests pass. Source-method/temperature/experimental
disorder and matched local-force/stress/electromechanical gates remain open;
this is seed preparation, not a new local minimum or quantitative benchmark.
The four native jobs continue under unchanged provider/protocol/source bytes.

### Provider progress: grouped motion and first amplitude check

All nine Sarco chains now have geometry readouts on their lowest internal
three-mode spaces for all three native curvature settings. Complete torsion-norm
vectors change at most0.04757553% over all27 setting-pair comparisons. Maximum
windows are stable across settings, but VDCN7/9 grouped maxima are3/1 backbone
bonds from nitrile, not the individual-mode0/2 distances. The grouped readout
avoids sign/internal-basis selection artifacts; it does not assign torsion energy,
prove defect causation or supply a field-driven pathway. See Sarco's
results/boundary_sensitivity/finite_mode_geometry_v1/SUBSPACE_GEOMETRY_REPORT.md.
The original all6993 native units and parent result were freshly requalified and
remained unchanged;25 combined geometry/space tests pass.

Both signed native xx amplitudes are now independently requalified. Proper
e_y,xx=-0.16748176972831194/-0.16725283962468349C/m2 at0.0025/0.005;
the complete coarse-versus-fine vector change is0.1366895881%, below the
predeclared2% threshold for this one column. Both signed pairs have zero branch
shifts under the0.25cycle bound. The new positive0.005 result SHA256 is
e83c37c428dd2b72f4566bf4da7c84d3f33ce7c8acc2776b2030aa53a6f77e1e.
The parent has advanced to polar yy=-0.0025. The full tensor/all-column gate,
proper yy/shear producer convention, matched-Born settings, induced-inclusive
polarization and mechanical qualification remain open. No inferred missing
electronic response or validated total/pC/N coefficient is promoted. All four
native jobs and the full five-stage chemistry scope continue unchanged.

### Provider result: three-mode field deformation is not generally complete

Sarco projects all108 accepted exploratory free/held field endpoints onto the
lowest three internal modes at all three native numerical settings, raw and
after removing collective axial swivel:648 Euclidean projections. See
results/boundary_sensitivity/finite_mode_geometry_v1/FIELD_MODE_COVERAGE_REPORT.md.
Original6993 native curvature forces and field endpoint/derived-result receipts
are requalified, with original writer identities retained. The32 combined
projection/geometry/space tests pass.

Aligned squared-displacement coverage spans approximately1.15-99.69%, varying
with chemistry/length/boundary/axis/sign. Holding azimuth does not cure truncation:
PVDF9 x branches are only1.15-1.16% covered; free PVDF7 z branches8.74-8.83%.
Worst coverage change across numerical settings is0.01731395percentage points.
Do not substitute a three-mode animation for complete simulated deformation,
or interpret geometric coverage as dipole/energy fraction or a coupling constant.

These archives contain relaxed endpoints, not complete vertical field-force
increments at the common zero geometry. A force-response prediction still needs
those matched forces, baseline subtraction, numerical/small-field checks and
sufficient resolved internal modes. No new electronic calculation or quantitative
material/RIS-barrier/packing/pre-strain claim is made. DFT/model/geometry, matched
crystalline/intrinsic stress and full five-stage validation remain open; all four
native CPU jobs retain their frozen source/provider/protocol identities.
Source/endpoint scalar force residuals for lower-coverage branches are retained
in the report. A stable projection space does not qualify endpoint accuracy;
tighter optimizer/electronic/zero-reference sensitivity is still needed before
assigning omitted shape to a physical higher-mode mechanism or to numerical
endpoint error. The current receipt readout does not distinguish those causes.

### Provider implementation: matched vertical field-force matrix

Sarco0a7dcaeb implements an owned234point fixed-coordinate CPU experiment:
all nine qualified AN/PVDF/VDCN5/7/9 references, zero and +/-5/65V/um on each
laboratory axis, electronic accuracies0.01/0.001. It reuses the original source
admission, calculator factory and native runtime guards. Parameters come from
the actual calculator instance producing complete raw forces, energy and total
dipole. The separate versioned self-validating vertical unit preserves the
existing zero-only curvature contract. Fresh zero points must reproduce the
original1e-4eV/A free-force gate, replayed from full units rather than stored
scalar assertions. All47 unit/lifecycle/projection/geometry/space tests pass;
synthetic lifecycle fixtures are not native SCFs.

Native preparation has been launched and must requalify the full6993force
parent before measurement. This is not yet a completed234point dataset or a
field-coupling/harmonic result. See Sarco vertical_fields_v1/PROTOCOL.md for
predeclared fields/settings, CPU1/resources, PID-verified recovery without
restart, full dataset requirement and compact failure retention. Direct forces,
small-field/electronic and tighter-reference checks remain ahead; all five
physics stages and quantitative material/packing/pre-strain/viewer gates stay open.

Separately, the complete native yy=-0.0025 vector is independently requalified:
result SHA256309ea28e917c9dea478a7dbeabe0e8b0f0ed41aa3d795566c00130cb10d8699d;
P=[1.0166272019086165e-9,-0.20580371668485187,3.299575564856335e-11]C/m2.
The original serial parent has advanced to yy=+0.0025, legally started
00:49:37.953964 PDT2026-09-13. No one-sided derivative, yy amplitude pass,
full tensor or matched-Born/mechanical claim is made. Original four native
jobs retain their frozen source/provider/protocol identities; GPU untouched.

### Matched-response reader and native measurement checkpoint

Sarco359320f4/f66203c9 add the typed vertical response owner and complete-matrix
report reader. All59 combined tests pass. Analysis subtracts the same-accuracy
zero point, retains full free-force vectors and separate support reactions,
reports signed odd slopes/even responses and compares amplitude separately from
electronic settings. Source and non-varied producer context must match;
undefined zero sensitivities remain undefined and numeric overflow is rejected.
The report requires all234declared units, then brackets analysis with independent
native terminal-evidence replay and guards analysis-code identity. There are108
signed pairs,54amplitude and54electronic-setting comparisons; no new pass
threshold or physical-model/harmonic validation is inferred.

Native preparation has successfully requalified the original6993force parent.
The same in-process CPU owner PID10059/start_ticks4699751 is advancing the234
measurements (191complete/43pending at this checkpoint, no failures). This is
not yet a complete dataset or a published native coupling result. All five CPU
owners are live, original epochs remain unchanged, GPU remains reserved for the
user, and all five research stages retain their outstanding validation gates.

### Complete native fixed-coordinate field-force result

The234point experiment completed without failure at01:07:52.507873 PDT
2026-09-13. Complete native replay passed before and after constructing all108
signed pairs and108 separate-setting comparisons. All18 zero points reproduce
the original1e-4eV/A gate; maximum9.043608658939799e-5eV/A. See Sarco
vertical_fields_v1/RESPONSE_REPORT.md and RESPONSE_ANALYSIS.json for findings
and complete branch statistics.

Worst odd force/dipole slope vector changes:0.0000586427%/0.0001033131% for
65 versus5V/um;0.0003662825%/0.0003018071% for electronic accuracy0.01 versus
0.001. Maximum signed-normalized force/dipole differences are0.1456887491%/
0.0148860302%. All comparisons are defined. These are numerical sensitivities,
not a new convergence/physical-model pass. Nuclei are fixed, so this experiment
does not observe chain rotation, unkinking or packing; total dipole changes
are electronic responses at these source geometries.

Native state SHA25642d7babb5ee1adebcde26877b917ba2599770868614e267f18083958c93bbb91;
resulta8e7547921e372752aa25bee99df88f873e023d6fe0514fceac820f681282dac;
analysis1e59bcebd2280373fb3803ebb977d516a6248acc65ab737a9c590f22f36ce0ec.
Matched force increments are now available for full internal-space curvature
response tests. Relaxed-endpoint comparisons must account for collective axial
swivel and matching field orientation, and still need endpoint/reference accuracy.
Independent DFT/model, mechanics, pathway, pre-strain and viewer gates remain
open; original four CPU jobs continue unchanged and GPU stays reserved.
Sarco1da19258 archives the complete dataset and readout. Fresh complete native
replay exactly reconstructs the published analysis; all237 native/derived JSON
entries and all five analysis owners match unfiltered raw-byte Git index hashes.
The59 tests pass. Four original live root identities and their26/8/10/11 frozen
owner hashes were rechecked before archive publication; no active epoch changed.

### Full internal Cartesian quadratic response test

Sarco bcddb2d2 implements and predeclares the next diagnostic: all324 signed
field-point estimates at the27 original Cartesian curvature settings, plus27
zero-reference Newton diagnostics and378 independent Cartesian-step/electronic/
signed-amplitude comparisons. The original point-clamp tangent/complement owner
defines the complete internal space. Strictly positive internal operators are
solved without clipping, regularization or pseudoinverse; support reactions and
removed axial driving load/torque remain separate. Matrix/vertical samples must
match complete Hamiltonian/settings, runtime and base producer. Lowest-three-mode
displacement coverage reuses the existing Euclidean projection owner.

All66 combined tests pass. The CPU-only readout is live as PID10956/
start_ticks4859696, observed advancing CPU time while requalifying the original
6993force archive. It will replay all6993 forces and234 field units before and
after deriving the report; source/numerical/analysis identities remain guarded.
No completed native harmonic readout is yet claimed. Code/protocol and producer
dependencies remain frozen for this active readout. Four original CPU root
identities and their26/8/10/11 owner hashes remain live/unchanged; GPU stays
reserved and all five research stages remain active.

The diagnostic fixes the linear axial tangent gauge. It is not the full nonlinear
azimuth-constrained Hessian, finite-angle swivel prediction or observed molecular
deformation. Source Newton correction is not an independently bounded displacement
error. Endpoint agreement, collective-swivel/field-frame matching, DFT/model,
packing, mechanics, RIS/barriers, pre-strain and quantitative viewer gates remain.

### Complete full-space result, matched frames and fine yy column

Sarco30477d5d/906e80ec add and archive corresponding geometry/field frames for
all108 receipt-requalified free/held endpoints. Geometry alignment also rotates
the applied laboratory field into the same reference frame. Original zero-source
coordinates, chemical graph/clamps, complete published field results and actual
historical writer availability are checked. Measured free swivel spans-173.539039
to175.900921degrees; held swivel is numerical zero. This is conditioning data,
not an orientation prediction. No endpoint coordinate archive is duplicated.
Frame JSON SHA2566576f2b10af9a5c09a737086224ad494b67a1e45ed50580a2d4aa88eefa51cd5.

Sarco985ac91a archives the complete324-point internal quadratic response,27
source Newton diagnostics and378 sensitivity comparisons. Both native owners
replay all6993 curvature forces and234 vertical units before and after math,
with exact result and original numerical identities. Worst displacement-vector
changes:0.3217657722% Cartesian step,0.0027782896% electronic setting,
0.0306722568% signed-normalized amplitude. No comparison is undefined.
Lowest-three-mode squared Euclidean displacement coverage remains3.876478-
99.763407%; full resolved internal space is still needed. These are numerical
diagnostics, not a fitted physical/harmonic pass. Backward-solve scaled residual
maximum1.9960641469416177e-16 verifies calculation of the specified operator only.

AN9 has0.008311788A RMS tight-setting zero-reference internal Newton correction,
comparable to its+65V/um x estimate0.012374730A, despite the original1e-4eV/A
force gate. Tighter native reference preparation/requalification is needed;
do not silently apply this Newton vector as an unmeasured reference correction
or interpret it as an independently bounded error/physical soft-mode result.
Internal JSON SHA256a8a719306eb3acd1f6f823caa2b602dfa0b6094f4ec0e39a155db2174d4e196b.
See Sarco finite_internal_response_v1/RESPONSE_REPORT.md for scope and diagnostics.

Separately, zero and both signed fine yy points are independently requalified
from complete native archives through their original point owners. Proper fine
Cartesian yy column at0.0025:
[9.534545513373371e-9,-0.15477449408819782,-1.0752758827570598e-8]C/m2;
both branch shifts zero. Positive state/result SHA256s:
d576e7ba4dd84497710f9e60be80e10b4b3e64390e19ded58357a7ffa9b3e86b;
dae45d2af9361019f28a265c57b6ce638d9a7b5b712968f6c0b571234bb49d65.
The original parent advanced to yy=-0.005 at01:24:32.737622 PDT2026-09-13.
No yy amplitude/full-tensor/matched-Born/mechanical gate is inferred.

All71 combined tests pass. New native/derived archives and recorded reader files
match unfiltered raw-byte Git index hashes. Four original CPU root identities
and their26/8/10/11 frozen owner hashes remain live/unchanged. GPU stays reserved,
and all five research stages retain DFT/model, geometry/reference, packing,
mechanics, pathways/barriers, pre-strain and quantitative-viewer requirements.

### Tighter native reference candidate preparation

Sarcof263957c adds an owned native candidate producer and explicit nine-point
refinement lifecycle, leaving the legacy hard-coded relaxation path and active
epochs unchanged. All nine original admitted AN/PVDF/VDCN5/7/9 sources are used,
AN9 first, zero field, electronic accuracy0.001, original exact point clamps and
held collective azimuth0. BFGSLineSearch/0.05A maxstep/500step cap targets
1e-6eV/A. A distinct fresh calculator at unchanged final nuclei must pass BOTH
raw and projected free-force gates. Actual parameters, full forces/total dipole/
energy, numerical/runtime/optimizer identity and original/refined coordinate
bindings belong to the self-validating RefinementCandidate unit. The original
constraint owner derives projection and support reactions; no copied formulas.

All83 combined tests pass, including12 new synthetic codec/fresh-calculator/
lifecycle tests. Independent failed points permit other declared references to
be attempted; useful individually complete candidates are retained, failures
keep diagnosis/eight-sample scalar tails only. Native result validation precedes
terminal state transition. PID-verified recovery requires authoritative absence;
there is no timeout restart/retry or automatic successor. Native preparation
passed and the same launcher is running as PID15909/start_ticks5023593, observed
advancing accepted steps on AN9. No completed native candidate is yet claimed.

These are candidates, not source admission, three-trial basin-return qualification,
Hessian stability, field response or physical/model accuracy. Old admitted source
coordinates and curvature/field archives remain unchanged. A distinct qualified
source epoch with new evidence is needed before substitution. All four original
CPU roots and their26/8/10/11 hashes remain live/frozen, GPU stays reserved, and
the full five-stage research retains its independent validation requirements.

### First tighter candidate and VDCN held-cell pair

AN9 completed2026-09-13 01:51:31.818717 PDT after182steps. A distinct fresh
native zero-field calculator ataccuracy0.001 passes raw/projected force gates
9.879325169e-7/9.879087714e-7eV/A versus1e-6. Actual original-to-refined
displacement0.008292026737A RMS agrees closely with old tight Newton diagnostic
0.008311787739A: full-vector difference0.4325374368%, cosine0.999993456186,
differenceRMS3.595159364e-5A. The original numerical stationarity concern is
supported; this is not field-induced motion. Candidate unitSHA256
e77faf55bbe02f300f4d5df02e29d84996455c031e9ceabbf38d7e0b89f5a828,
milestoneJSON5044a687012680484f9358c885fda5e625ab93b83cb31f517f259151c4ef032d.
See Sarco reference_refinement_v1/AN9_MILESTONE.md. AN5/AN7 also complete;
the same frozen nine-point run has advanced to PVDF5, GPU remains reserved.
Old source coordinates are not substituted. New admission/three-trial basin
return and new Hessian/field/independent model evidence remain required.

VDCN antipolar held-cell completed01:48:20.658247 PDT after121steps/122SCFs;
childexit0, original writer-owned representation auditPASS andforce0.00402720022
eV/A below0.005. Original run now relaxes polar full-cell. Both completed
held-case/geometry receipts are archived without sealing the active parent;
FIXED_PAIR_MILESTONE.json SHA256
5d763dc24de91a1fedcef83b7e84bd8905978fd258a501aa063e00456951bb44.
The reader reuses original stage endpoint loading and guards runtime/protocol/
sources, both geometry hashes, exact held cells and periodic topology. Native
raw outputs were originally audited then removed by their retention owner,
not freshly replayed here. Antipolar-minus-polar=-8.770059399meV/monomer under
DIFFERENT held cells, not equilibrium phase ordering. Held stresses are not
intrinsic pre-strain; full-cell/stability/mechanics/model/path gates remain open.
All22 focused synthetic refinement/milestone tests pass, ten new; no native
SCF attestation is inferred from tests. All five original/refinement roots and
13/26/8/10/11 hashes remain live/frozen; no completed five-stage claim.

Sarco126cfc79 commits/pushes the AN9 and VDCN pair milestones and ten new tests;
only13 intended files, all unfiltered raw Git index bytes exact, no active
journals/user edits included. PVDF5 subsequently fails its DISTINCT FRESH
raw/projected force gate at01:58:16.926133 PDT despite optimizer step78 cached
projected9.832666631e-7eV/A passing1e-6. It is not a cap failure or accepted
candidate. The current exception does not retain separate fresh maxima, so
which force failed/how far above threshold/noise/support-torque attribution
remain unknown. Sarco PVDF5_DIAGNOSIS.md records that diagnostic limitation
and future force-owner scalar reporting requirement without changing the active
epoch. No rejected geometry/log/trajectory/checkpoint or accepted PVDF5 unit
is retained. Same run continues PVDF7, no retry/lowered gate/substitution;
three accepted AN candidates remain useful, but all-nine cannot pass this attempt.

### Complete frame-conditioned quadratic endpoint readout

Sarco6b58c57e implements the comparison at its defining boundary and archives
the independently requalified larger negativeyy=-0.005 point. Its original
point owner replays complete submitted inputs/native stdout/XML, clean returns,
writer/protocol and exact result; statee226b4d186f17f5f196e7d78a8be36f1676da9ad907d3ceef0ecb48fbefbb37b,
result944768723e9f32ab626beacd8a93d60e3d2f37045ebe03e05e05d48673818acc.
P=[1.017780088e-9,-0.2054167882994,-2.925575925e-11]C/m2. Parent advances
positive coarseyy; no larger-pair/amplitude/full-tensor/Born/mechanics gate.
All20 intended files' unfiltered raw Git index bytes match working bytes,
including native whitespace. PVDF7 tighter refinement also fails its distinct
fresh-force gate despite step58 cachedprojected9.039168231e-7; same missing fresh
maxima diagnosis, no accepted unit/rejected geometry/log/checkpoint, no retry.

finite_quadratic_endpoint.py uses six matched signed axis responses to build
the central-odd local map. Its endpoint comparison rotates BOTH actual geometry
and laboratory field through the existing defining frame owner, compares the
entire free vector, and reports the separate error from leaving the field
unrotated. Measured swivel is conditioning, never predicted angle. Preserve all
108 endpoints, three Hessian settings and both5/65axis-probe amplitudes:648
comparisons, not an easier subset. No fitted gate, Newton-offset subtraction,
physical higher-mode attribution or unmeasured mixed-field nonlinear pass.
All21 focused tests pass, nine new synthetic math/complete-plan/lifecycle tests.
Sarco4e397a8f seals/pushes the ready journal before CPU execution. The same
readout now runs asPID25373/start_ticks5161570, session85293; kernel liveness
and all23 analysis-owner hashes verified. Original complete6993force/234field/
108endpoint evidence must replay before/after comparison and reproduce old324
quadratic/108frame reports exactly. No completed comparison yet; GPU reserved.

### Actual ideal held-angle slice is affine

Read-only analysis of the defining AxialFrame finds C(X),S(X) affine for fixed
clamp axis/origin, so holdingphi iscos(phi)S-sin(phi)C=0 on its positive branch.
Atphi0 the old full internal complement is the exact affine slice. Its analytic
restricted potential Hessian isC.T H C; projected angle-Hessian and support
reaction-curvature terms vanish. The analysis reuses actual owner measurements/
gradients, not a parallel constraint implementation. Nine original typed
sources x three targets x two derivative steps give54 guarded mathematical
probes: max angle4.44e-16rad, projected angle-gradient4.405e-17A^-1, projected
angle-Hessian-direction2.585e-14A^-2. Five synthetic tests pass, including wrong
branch and a nonstationary normal-force quadratic. Sarco azimuth_slice_v1/
SLICE_ANALYSIS.json SHA256
350dc22bb889023bfe8c6b4a0876957c061fd42bb329d95a84c1be1475eeb88a;
SLICE_REPORT.md gives the derivation and scope. This corrects an unnecessarily
conservative missing-constraint-curvature caveat, not sampled-H/source/electronic/
finite-field/free-swivel/packing/model/physical uncertainty. No force replay or
new reference/material gate is inferred; all five research stages remain active.

Sarcof803471c publishes the exact affine-slice derivation/54-probe archive;
all six intended staged byte units and all four recorded analysis owners match,
no active journal/user edit. All26 combined slice/quadratic/frame/report tests
pass. The same live endpoint readout accepted initial complete6993force/234field
replay and is deriving324 original full-space responses, not yet publishing648
comparisons. PVDF9 tighter reference subsequently hits its500-step cap at
02:15:52.637834 PDT with cachedprojected2.285520598e-6>1e-6; it never reaches the
fresh candidate check. PVDF_DIAGNOSES.md distinguishes this cap fromPVDF5/7 fresh
gate failures, preserving only relevant scalar/source/time diagnosis. No rejected
geometry, accepted minimum energy, retry/cap extension/lowered gate. Original
refinement continues VDCN5; AN candidates remain useful. Six kernel roots and
23/13/11/26/8/10 hashes remain live/frozen, GPU reserved, full scope unchanged.

### Complete648 readout and partial tighter references

Sarco7529abd2 publishes all648 endpoint comparisons:108 endpoints x three full
Cartesian Hessian settings x two signed-axis probe amplitudes. Original6993force/
234vertical-field/108endpoint native evidence replayed before/after; original324
response/108frame reports reproduced exactly. Result raw SHA256
00aab36d82c8244bd5deae1a84570e345608f80c36031aec5ca77dd62ca98ee7.
See finite_quadratic_endpoints_v1/COMPARISON_REPORT.md. Median full-internal vector
mismatch freeAN/PVDF/VDCN4.179/5.365/5.861%; held7.731/10.096/8.831%; largest
55.149% heldAN9-z65. Measured swivel conditions geometry AND field, not predicted
swivel. No reference Newton offset subtraction or physical higher-mode attribution:
old stationarity/endpoint accuracy and unmeasured mixed-field effects unresolved.
No harmonic/model/bulk/viewer gate is claimed.

The original tighter-reference epoch finished partialfailed with five native
fresh-force-accepted candidates:AN5/7/9,VDCN5/7. PVDF5/7 freshcheck failures;
PVDF9/VDCN9 original500-step caps. Accepted units replay and exactly reproduce
published partial resultaa19a65572ef3dfb664d14c01b37d117b4d7c9487d0a946c44638564690c8655.
See reference_refinement_v1/TERMINAL_REPORT.md. No retry/looser gate/new source
admission or substitution into old Hessian/field evidence. Failed geometry/logs/
checkpoints are not retained; exact PVDF5/7 fresh maxima/positions cannot be
reconstructed from the compact diagnosis.

Distinct force-history code/protocol and twelve synthetic codec/native-receipt/
lifecycle/recovery tests are sealed/pushed; no synthetic native attestation.
Preparation captures14 authoritative inputs (nine originals plus five accepted
refinements),28 pending units at accuracies0.01/0.001,364 native force observations.
Each uses seven distinct fresh calculators, six independent +/-xyz0.001A seed
routes returning to IDENTICAL nuclei with writer-observed preserved native API/
previous result. This measures numerical calculation-history sensitivity, not
retroactive failed-PVDF endpoint causes, source admission or independent SCF/model
error bounds. Ready journal is published before CPU execution; GPU reserved.
Four other original roots remain live with26/8/10/11 owner hashes unchanged;
clamped campaign completed both coarseyy points and advanced zz=-0.0025. No yy
amplitude/full-tensor/Born/compliance gate from this status observation alone.
All five research stages and other target chemistries remain active.

Sarco b238566d publishes the ready14-input/28-case force-history snapshot before
native execution. Actual new CPU owner PID24914/start_ticks5329124 is kernel-live;
all10 code-owner hashes verified. No force-history result yet, no retry/restart.

Fresh requalification of zero and all four signedyy point native archives gives
the completed coarse properyy column
[1.994858750e-9,-0.1544171918967347,1.203499220e-8]C/m2. Coarse current versus
fine reference full-vector change0.2308534067% passes the predeclared2% numerical
amplitude limit; branch shifts allzero within0.25cycle bound. Positive coarse
terminal state3a6868cc567459e126dac8be3cd95cecf4da2b7e526dda5707741704a42d81c2,
result449b2168571874266c059e1a3d0a82522902cc0d3ffff7272c15bfa59c0602bf.
Defining point/polarization/strain/vector owners replay actual submitted/native
input/stdout/XML, writer/protocol/clean returns and exact results. All11 native
files retained byte-for-byte. See Sarco clamped_ion_v1/EXECUTION_REPORT.md.
This is only yy amplitude pass, not accuracy of near-zero transverse components,
full six-column tensor, unmatched Born subtraction, compliance or physical-film
data. Original parent continues zz=-0.0025; GPU reserved, five-stage scope active.

### Completed native force-history diagnostic

All28 units/364 native forces completed, with98 same-seed-history accuracy
comparisons. Exact declaration/runtime/native-receipt/nested-CRC and full
reduction verification passes. Sarco force_history_v1/result.json rawSHA256
b00a97a01953344ec63537d67a9f711833227343236d3f6d11d47fcb482fb929;
terminal stateaa62a0244958f2c6ecd73901b7eb64a389308f6ffd8cd617078b08b499961e73.
Largest identical-nuclei free-atom force difference6.521860107e-6eV/A ataccuracy
0.01 versus5.256515174e-7 at0.001; energy differences remain<=2.23e-11eV.
Freshaccepted AN7/VDCN7 cached returns cross1e-6 rawforce (max1.05707e-6/
1.02583e-6). Original fresh acceptance receipts are unchanged, not claims of
cached-history immunity. This cannot reconstruct discardedPVDF5/7 failures or
assignPVDF9/VDCN9 cap causes. See RESULT_REPORT.md;31 combined focused tests pass.

Next numerical method test: fresh electronic solve at every optimizer evaluation,
reusing original optimizer/constraint/candidate path behind a native policy owner;
sameaccuracy0.001/1e-6fmax/500cap/.05Amaxstep, no physical knob/gate retuning or
edits/substitution into old epochs. CNEPO source-motif5/7/9 prepared inputs already
exist and independently requalify (manifestef82a9bcce6e5ec4fd793a307565f62e4bee0ebc9bfd75ccd696ea124ebcd596);
they are not yet relaxed. No duplicate builder/preparation is needed. Full
five-stage scope, independent model/DFT, ring/other-chemistry, packing, stress,
barrier/temperature/bulk and accurate-viewer gates remain open; GPU reserved.

Sarco d14d6690 publishes complete native force-history evidence/report;32 raw
staged units match unfiltered Git bytes. Sarco aeb575cc seals the distinct
cold-every-evaluation refinement policy/driver/protocol and nine new tests.
ColdForceProvider wraps the ONE original native factory; ColdReferenceRefiner
INHERITS the original optimizer/constraint/fresh-terminal/candidate path, not a
copied minimization implementation. Backend owns actual fresh-before/no-previous
result/native-return/evaluation-count/last-sample geometrySHA and separate raw/
projected force maxima. ASE actual constraints own projection. Producer supplies
receipts; generic lifecycle engine never stamps native truth. Original epochs
and dependencies remain unchanged. Tests include nonzero real BFGSLineSearch
steps with a synthetic quadratic native calculator, nine-case codec/lifecycle/
receipt replay and parameter/cached API/context/bool/count/policy/matrix/diagnostic/
interruption negatives; synthetic calculators attest no physics.

Preparation is replaying prior sealed complete force-history/partial-reference
evidence before capturing the new original-nine ready input/diagnostic snapshot;
no actual new optimization launch yet. Sameaccuracy0.001/1e-6raw+projectedfmax/
500cap/.05A/held0pointclamps, no tighter-candidate substituted starts, threshold/
physical knob retuning or all-nine/source/model gate. Four original roots are
kernel-live,26/8/10/11 owners unchanged,17.5GB WSL memory available and swapfree.
GPU reserved; CNEPO prepared-only relaxation/admission and full scope remain open.

Sarco eefdeb3b seals/pushes the actual original-nine ready epoch before launch.
Distinct CPU cold-refinement is now actually optimizing AN9 as PID27830/
start_ticks5430020 (session34023). Kernel identity/liveness and all17 owner hashes
verified; initial optimizer step2 observed, no accepted candidate or terminal
result. Invocation requalified sealed diagnostics before ready->running; no
observation-timeout restart. All40 combined fresh-native-policy/lifecycle/history/
slice/frame/quadratic tests pass. Actual fresh force-history results remain
published, old epochs/source admissions unmodified. Same1e-6raw+projected gates/
accuracy0.001/500cap/.05A/held0pointclamps, all dependencies frozen. Other four
science roots remain live with26/8/10/11 hashes intact; GPU reserved and full
independent-model/packing/stress/field/other-chemistry/temperature/viewer scope open.

### CNEPO prepared inputs to independently owned reference candidates

Sarco82fb8fe0 seals/pushes coordinate_cold_refinement.py and CNEPO reference
driver/protocol/five tests. CNEPO source-motif5/7/9,32/44/56atoms retain one
three-edge epoxide, four actual stereocenters, twoTrFE and symmetricVDF growth.
The checked preparation already exists; no pendant-only or pure-VDF host is
substituted. PreparedCNEPOInputs owns chemical/stereo/point-clamp INPUT admission;
the original native factory owns unchanged Hamiltonian/parameters/runtime and
its unrelated old-nine calibration. Two distinct self-described admission hashes
avoid claiming CNEPO is qualified by old-nine data. Existing candidate context
source_admission_sha256 binds the actual prepared CNEPO registry, not old-nine
admission. Factory calibration does not claim a relaxed molecular minimum.

ONE inherited cold optimizer/constraint/candidate/native audit/receipt reader
and generic lifecycle path is reused; no copied minimizer/projection/native
factory/codec. Sameaccuracy0.001/1e-6raw+projected force/.05A/500cap/held0 terminal
xyzclamps, no knob/gate/cap retuning. Five tests pass with real checked inputs
and synthetic electronic returns:3-case terminal replay, separate admissions,
ring/stereo/motif/actual clamp checks, input/evidence/coordinate/gate/matrix drift.
Tests attest no native electronic physics. Actual ready3-case input/factory epoch
is being published before run; no CNEPO optimization result yet. Five previous
roots kernel-live,17/26/8/10/11 owners unchanged,22.7GB WSL memory available and
swapfree. GPU reserved. Candidate successes still need basin/Hessian/field/size/
bulk/DFT/model/packing/stress/barrier/temperature/viewer gates; full scope open.

Sarco678bc472 seals/pushes the actual3-input/factory ready epoch BEFORE run.
CNEPO CPU reference producer now PID15514/start_ticks5486377,session82507;
kernel liveness/all20 owner hashes verified. FirstCNEPO5 optimizer step3 observed,
projected0.1364493952eV/A; no accepted candidate/terminal result. Original ring/
stereo invariant is checked by the inherited original algorithm. All five new
sealed code/protocol/test/doc units match exact committed bytes. Original cold
AN9 stays kernel-live PID27830/start_ticks5430020 with17 hashes unchanged (step138
projected2.114089882e-5>1e-6, not accepted). Other four roots retain26/8/10/11
owners/liveness; all six epochs remain frozen, GPU reserved, five-stage scope open.

### Explicitly sampled alpha/gamma/PE control compiler

Sarco now has cp2k_sampled_method.py: one typed full-grid, non-centered,
non-symmetry-reduced complex Monkhorst-Pack sampling owner, reusing the
existing defining physical method/exact geometry compiler and original
diagonalization/Broyden SCF. Static/mock-interactive inputs are identical;
no frozen Gamma/VDCN source/include or live numerical dependency was edited.
All35 combined compiler/native-syntax/source tests pass. Actual CP2K2026.2
--check accepts original24/48/12atom alpha/gamma/PE source cells at explicit
4x4x4 and rejects an invalid KPOINTS keyword. Exact positions/full cell vectors
are checked; PE's audit supercell is not substituted for its source cell.
SAMPLED_CONTROL_METHOD.md records native/input hashes and evidence limits.

These are computed provenance-admitted INPUTS, not local electronic minima;
syntax acceptance is not an SCF, force/stress or convergence gate. Equal
division counts across different cells do not imply equal reciprocal spacing.
No sampled control electronic run has yet started. Shared provider lifecycle/
receipts, full-cell relaxation, basis/sampling, phase/mechanical/electrical/
physical gates and experimental full xyz remain open. All six science roots
are kernel-live with20/17/26/8/10/11 owner hashes unchanged, GPU reserved.
The full five-stage/other-chemistry/temperature/accurate-viewer scope stays open.

Original fresh-every-evaluation AN9/AN5 now complete in181/76steps. Existing
independent candidate/native-receipt reader passes; fresh raw/projected maxima
3.4060293133e-7/3.4062455548e-7 and7.7181695698e-7/7.7136988927e-7eV/A
are below unchanged1e-6gates at accuracy0.001. Native evaluations308/152;
original-to-candidate free RMS0.0082927595420972A/5.9912723760359606e-5A.
Raw accepted units and an exact observed journal snapshot are retained under
Sarco cold_reference_refinement_v1; ACCEPTED_AN9_AN5_REPORT.md gives hashes,
receipts and limits. Snapshot is running-parent context, not terminal/all-nine
evidence. AN7 continues in the same process. No canonical source replacement,
basin/Hessian/field/model/size/bulk/crystalline-prestrain/physical gate is claimed.

### Sampled native session and first computed-control electronic epoch

Sarco's sampled native session inherits the existing ONE clean native exit,
convergedSCF and output-vs-shell energy/fullforce/fullstress audit; explicit
RUNNING->CLOSED->COMPLETE/FAILED, idempotent close/no failed finish retry.
Actual kernel executable, initial compiled source cell/input and unchanged
method are checked, no frozen Gamma/VDCN owner changed. A real CPU2x2x2
periodicH2 integration returned oneSCF,exit0 and matching observables in
10.644seconds: kernel/transport integration, not polymer validation.

sampled_control_points.py predeclares one native zero-field ENERGY_FORCE
point per original24/48/12atom alpha/gamma/PE input cell on full-grid4x4x4.
No relaxation/source-supercell substitution, field, physics/gate retuning
or reciprocal-spacing convergence claim. Existing IndependentNativeUnits
lifecycle/reduction is reused; point codec owns whole canonical checksum,
reader binds actual source/method/nativebytes/receipt/fullarrays. Successful
rawinput/output/unit survives later failure; failedscratch keeps diagnosis
only. All43 combined source/compiler/session/campaign/native tests pass.
Ready identity covers16codeowners, nativebinary/data and actual imported
Python/ASE/NumPy transport/unit/core files with explicitly limited coverage.
Actual three-point ready epoch is sealed before run; no polymer result/minimum
claimed by preparation. Six prior science epochs untouched, GPU reserved;
full five-stage/grid/basis/packing/stress/field/temperature/viewer gates open.

Sarco9cbb773c seals/pushes actual provider/protocol/tests and exact three-case
ready epoch before run. CPU producer PID28624/start_ticks5630900,session56867
is kernel-live with16hashes unchanged. Alpha child PID31440/start_ticks5631259
is kernel-live and has actual native diagonalization/Broyden SCF steps1-5,
24atoms/96electrons/248orbitalfunctions; no accepted point yet, gamma/PE pending.
Final combined43tests pass in45.594seconds with native integrations enabled.
Other six science roots remain live,20/17/26/8/10/11 hashes unchanged, GPU
reserved. No frozen live source/dependency changes or observation-timeout
restart. Seed electronic points do not establish minima, sampling convergence,
phase ordering, intrinsicprestrain, electrical response or full-goal completion.

Sarco results/VERIFIED_NUMERICAL_MILESTONE_20260913.md publishes independent
native-reader checks for alpha point, AN7/PVDF5 cold candidates and CNEPO5.
Alpha unrelaxed computed input has1.4320197071eV/A force and normal stress
2.17-5.57GPa magnitude: not intrinsicprestrain or a local minimum. Exact native
input/output/unit retained. AN7/PVDF5 pass unchanged freshraw+projected1e-6
after76/62steps. Actual CNEPO source-motif5 passes143steps, raw4.9929742697e-7,
projected4.9314833026e-7;0.3251167254A RMS is zero-field constrained relaxation,
not field-induced motion. Ring/stereo/source vs factory admission stays distinct.
Accepted units and actual running-parent snapshots retain writer/context, no
terminal-parent/source/model/Hessian/field/physical qualification inferred.

sampled_control_cells.py owns source vs calculation-cell preparation. Alpha/
gamma24/48 unchanged; PEsource12 vs explicit1x1x2 calculation24, folded4x4x2
instead of silently denser4x4x4. Three actualinput/math tests pass0.793seconds,
exact geometry/two-chain topology, reciprocal steps/full-grid folding and
multiplicity, unknownphase rejection. Native folded-cell equivalence and
staged phase-scoped atomic/full-cell method remain next; no filesystem-label
Hamiltonian inference. All seven live epochs frozen, GPU reserved, five-stage
and independentmodel/packing/stress/field/temperature/accurateviewer scope open.

### Completed native controls and next calculation-cell/relaxation owners

Original three-point root is terminal complete, session56867 exit0/no owner.
All original raw native inputs/outputs/units and exact reduction reread:
state13ca4fa8c18aac53e39b00d9a90807e5edf246ff3bda14d8ac943b859851d224,
result6c6cdb02b25a06bfa8b4bbf3f774e8ea44284b60d58f93665301099180f7861f.
Gamma maxforce1.2383190215eV/A and PE0.2909623627 are unrelaxed computed-seed
results, not intrinsicstress/equilibrium/material output. Sarco point RESULT_REPORT
records full evidence and normal stresses. Historical writer identity preserved.

Shared point provider now owns typed per-instance sampling/phases/storage/
protocol. Only changed AFTER old root terminal, exact old reduction reverified.
folded_pe_point.py reuses ONE native/codec/lifecycle path for explicit PE24atom/
4x4x2, with primitive vs repeated comparison gates0.001eV/primitive,2e-5eV/A
forcevector,0.1MPa stresscomponent. Mathematical folding does not assert native
equivalence; pending native study has18 frozen owners and actual ready epoch.

sampled_control_relaxation.py uses phase-scoped source/method and existing
two-stage crystal lifecycle: atoms first at held cell,0.005eV/A, then actual
accepted own predecessor with all6 cell components free and1MPa stress.1000cap/
.05A/zero externalstress/field; no filesystem-label mesh inference or physics
retuning. PE12atom native point alone cannot admit24atom calculation cell;
completed folded-cell equivalence required. Seed label is not measured relaxed
RIS/packing/spacegroup phase. Actual filter zero-pressure default and optimizer/
filter/SciPy file identities checked with limited coverage explicit.
26combined tests pass95.051s, final3scoped tests30.127s; synthetic tests not
polymer physics. Alpha exact ready2-stage/20owners and PE exact ready1-point/
18owners are published before run. Six prior roots remain live20/17/26/8/10/11
owners unchanged, GPU reserved; full five-stage and physical gates remain open.

Sarcoaf7a1d17 seals/pushes code/protocols, exact original control terminal
native bytes/reduction and both new actual ready epochs before launch.
Folded PE now PID25016/start_ticks5779027,session99781,18hashes unchanged;
actual child2232/start_ticks5780440 live with SCF steps1-7. Alpha2stage now
PID24995/start_ticks5779026,session25721,20hashes unchanged; actual first-stage
child29952/start_ticks5779748 live with SCF steps1-5. No accepted folded native
point/relaxed alpha stage yet. Other six actual roots retain kernel liveness
and20/17/26/8/10/11 hashes. All eight epochs frozen, GPU reserved, no timeout
restart. Native equivalence/force/stress and full five-stage gates remain open.
