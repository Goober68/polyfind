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
