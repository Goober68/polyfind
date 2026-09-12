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
