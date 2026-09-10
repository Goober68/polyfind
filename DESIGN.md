# polyfind: efficient determination of stable atomic arrangements in semi-crystalline polymers

Target case: poly(vinylidene fluoride), PVDF, -(CH2-CF2)n-.  About half
crystalline; near-degenerate polymorphs whose chains are the torsion sequences
alpha = TG+TG-, beta = TTTT (all-trans), gamma = TTTG+TTTG- (delta and epsilon
are the polar/antipolar partners of alpha and gamma); strong F...F repulsion and
CF2 dipoles; an amorphous fraction whose local conformational statistics decide
which polymorph nucleates.  "Stable atomic arrangement" therefore means three
different things, and the method treats them as three levels of one problem:

1. the low-energy *chain conformations* (torsion sequences),
2. the *crystal packings* of periodic chains (polymorphs),
3. the *equilibrium ensemble* of the amorphous fraction and its interfaces with
   crystalline lamellae.

## 1. Why the usual approaches are slow

| Method | Cost driver | Failure mode for PVDF |
|---|---|---|
| Molecular dynamics / simulated annealing | 10^6-10^8 force calls; equilibration is impossible below Tg and slow above it | Trapped in local minima; polymorph ranking never converges; ns per chain |
| Cartesian random / evolutionary crystal structure prediction | Search dimension 3N per cell | Chain connectivity is a hard constraint the search does not know about; most trial structures are broken chains |
| Periodic DFT on guessed cells | Hours per structure | Only as good as the guess; cannot touch the amorphous fraction |
| Classical RIS theory (Flory) | Analytic | Gives statistics of isolated chains only; no route to 3D structures |

## 2. The approach

Solve the problem at the level where it is low-dimensional: the sequence of
backbone torsional states.  Every step that is normally stochastic becomes exact
and discrete, chain symmetry makes the 3D packing low-dimensional, and the
expensive potential is used where it changes the answer and nowhere else.

### 2.1 Conformation search is an exact dynamic program (`ris.py`)

Backbone dihedrals occupy discrete rotational isomeric states (T, G+, G-;
more states are allowed).  With first-, second- and optionally third-neighbour
interaction energies the chain is a Markov field, so all of the following are
exact and cost O(N S^2) via transfer-matrix / Viterbi recursions:

* minimum-energy conformation (Viterbi),
* the k best conformations (k-best Viterbi),
* partition function and conformational free energy; Flory's largest-eigenvalue
  free energy per repeat of the infinite chain,
* per-bond marginals (forward-backward),
* exact Boltzmann samples, batched, with arbitrary bonds clamped.

Heteropolymers get one pair matrix per bond type: for PVDF the pairs alternate
between CH2-centred and CF2-centred, which is where the F...F repulsion and the
pentane effect enter.  Third-neighbour terms are handled by *state
augmentation* (pairs of consecutive states become the states of an equivalent
second-order model), so the same solvers apply unchanged.  Those terms matter
for PVDF: with pair interactions only, TG+TG+ (a 3/1 helix) and TG+TG-
(alpha) are exactly degenerate; the triplet G+ T G- versus G+ T G+ is what
separates them.

### 2.2 Crystalline chains are periodic sequences (`enumerate.py`, `helix.py`)

Polymorph chains are cycles of length P in state space.  A cyclic k-best
dynamic program returns the lowest-energy cycles for each period directly, so
the 3^P sequences are never enumerated.  Cyclic shift by a repeat unit, chain
reversal and mirror (G+ <-> G-) are symmetries of the isolated chain's energy;
candidates are canonicalised under them, repetitions of shorter periods are
dropped, and sequences that close on themselves (zero rise) are discarded.

### 2.3 Helical symmetry collapses packing to a handful of variables (`helix.py`, `pack.py`)

Superimposing one period of the built chain onto the next (Kabsch) yields the
rigid screw motion that generates the infinite chain.  Decomposing it gives the
chain axis, the rotation and rise per period, hence the crystallographic repeat
c, the number of periods per turn and the chain radius, with no fitting.  The
chain becomes a rigid helix, and a two-chain cell (PE, alpha-, beta-, gamma-
PVDF all have chain 2 at (1/2, 1/2)) is fully described by

    a, b, gamma, phi1, phi2, dz, parallel/antiparallel

six continuous variables and one flag, instead of 3N atomic coordinates.  The
lattice energy (Lennard-Jones + damped-shifted-force Coulomb over all periodic
images, including the intra-chain terms with bonded exclusions across the
periodic boundary so that different chain conformations are comparable) is
evaluated by one batched kernel for thousands of candidate cells at once; a
coarse random screen is followed by Nelder-Mead polishing of the best few.

### 2.4 Continuous refinement inside the basin (`refine.py`)

The discrete stage uses ideal angles and rigid geometry; real chains deflect
from the ideal RIS angles.  (This sentence used to give two examples,
"beta-PVDF dihedrals near +/-172 deg, alpha gauche near +/-45 deg".  The first
is withdrawn - beta's torsions are exactly trans, section 5.8 - and the second is
unverified, sources putting alpha's gauche anywhere from 45 to 62 deg:
`docs/REFERENCES.md` sections 2.2 and 2.4.)  The
torsions of one crystallographic repeat and the cell are relaxed together
against the lattice energy, with the chain kept periodic by a quadratic penalty
on the residual rotation of the repeat transform (the linked-atom / helical-
constraint idea of fibre-diffraction refinement).  c is not a free variable: it
follows from the torsions.

### 2.5 Multi-fidelity funnel (`forcefield.py`, `pipeline.py`)

| Stage | Evaluations of the expensive potential |
|---|---|
| Fit RIS energies: 1D and 2D dihedral scans of a short oligomer, plus 4 S^3 triplet corrections per bond type | ~2.9k single points for PVDF (step 10 deg, third order) |
| Enumerate and rank periodic sequences (DP) | 0 |
| Pack the top candidates and the known polymorph chains (batched kernel) | 0 |
| Refine torsions + cell in each basin | 0 (cheap lattice energy) |
| Final re-scoring / relaxation of the top few | O(10) |

The scan and the re-scoring can be run with the built-in potential, with any ASE
calculator via `ASECalculator`, or with DFT; everything downstream inherits that
accuracy.  The scan is a batch of independent single points and goes through
`energy_batch`, so a GPU-resident potential evaluates it as one batch.

**Decision: no machine-learned-potential dependency, deliberately.**  An earlier
version of this section called a foundation MLIP "the intended production path".
That is withdrawn.  The reasoning, recorded here so it is not relitigated by
someone arriving without the context:

* The novelty of this package is the *search* - exact dynamic programming over
  torsion states, the screw decomposition, the tabulated chain-pair interaction
  with FFT lattice sums, the line-group parametrisation.  The potential is an
  input to that machinery, not part of it, and reaching for a heavier input is
  not progress on the method.
* A foundation MLIP is not obviously more trustworthy here than a classical
  potential.  Such models are trained overwhelmingly on inorganic crystals;
  semicrystalline fluoropolymers are far outside that distribution, so an
  unvalidated foundation model is a slower and more opaque unknown, not a
  better one.  It has to be validated against known structures before it has
  any authority, and if it is being validated against the same reference data a
  classical potential could be fitted to, fitting is the cheaper answer.
* What is actually wrong with `SimpleFF` is not that it is cheap, it is that it
  is *unfitted*.  The funnel asks very little of a potential: torsion profiles
  for a few bond types and nonbonded parameters for a handful of element pairs.
  That is a small, well-posed fitting problem - a few hundred well-chosen
  reference points, not thousands - and the result stays at microseconds per
  evaluation, which is what keeps the exhaustive search viable.

So `ASECalculator` stays, because supporting an external calculator costs
nothing and an MLIP is a reasonable *validator* at the end of the funnel.  But
no MLIP is a dependency of this package, none is installed, and the productive
direction is fitting the cheap potential to the chemistry class of interest
against reference data that already exists.

### 2.6 The amorphous fraction from the same transfer matrix (`amorphous.py`)

Exact backward sampling gives the equilibrium single-chain ensemble at any
temperature with no equilibration run.  From it: state and diad fractions,
trans-run statistics (the fraction of bonds in runs of >= k trans states is the
supply of beta-like nuclei), the TTTT tetrad fraction, end-to-end statistics
and the characteristic ratio, and the conformational free energy per monomer.
Clamping a segment to all-trans and sampling the rest models chains leaving a
crystalline stem: mean advance along the stem axis, turn-back probability
(loops), lateral excursion, and the trans excess in the first free bonds.

## 3. Complexity and measured performance

Measured on this machine (4 CPU cores, NumPy, no GPU), PVDF three-state model.

| Task | Conventional | polyfind |
|---|---|---|
| Minimum-energy conformation, 1,000-bond chain | MD/MC: 10^6-10^8 force calls, no guarantee | 5 ms, exact (3^1000 ~ 10^477 conformations); 51 ms for 10,000 bonds |
| 50 lowest conformations, 1,000 bonds | not available | 0.12 s |
| Periodic candidates up to period 12 | 3^12 = 531,441 sequences, each built in 3D | 42 ms for the 540 best cycles, then symmetry dedupe |
| Partition function / marginals, 1,000 bonds | thermodynamic integration | 0.10 s / 0.20 s, exact |
| Equilibrium ensemble, 10^5 chains x 200 bonds | ~us of MD per chain, not equilibrated below Tg | 5.5 s, exact (backward pass is the batched part) |
| Lattice energy of one packed cell | 3N-coordinate relaxation | 1 ms (PE), 3 ms (alpha-PVDF), 10 ms (gamma-PVDF) |
| Crystal packing of one chain conformation | CSP over 3N coordinates | 6,000 cells + polishing: 11 s (PE, beta), 24 s (alpha), 62 s (gamma) |
| High-fidelity potential calls | 10^6+ | ~3k for the fit + O(10) at the end |

All dynamic-programming routines are checked against brute-force enumeration in
the test suite (290 tests, 4 skipped).

## 4. GPU design

Where a GPU helps and where it does not was decided by the shape of each
computation, not by default.

**Not on the GPU.**  The Viterbi / forward-backward recursions run along the
chain with a 3-9 element state vector per step; they are sequential in N and
tiny in S, so a GPU would spend its time on launch latency.  They stay in NumPy
and are already milliseconds.  Screw decomposition, symmetry canonicalisation
and Nelder-Mead polishing are likewise scalar work.

**On the GPU (batched kernels, `backend.py`).**  Everything written against the
NumPy-subset API that CuPy implements identically, so the same source runs on
both; `POLYFIND_DEVICE=cuda|cpu|auto` selects the backend, `to_numpy` brings
results back for SciPy.

| Kernel | Batch dimension | Why it is GPU-shaped |
|---|---|---|
| `CrystalPacker.energy` | thousands of candidate cells | M x I x N x N element-wise tensor: distances, LJ, erfc (elementary-op approximation, no special-function dependency), masks, reductions.  Per configuration only lateral images whose chain axes can lie within the cutoff are selected; the (0,0,k) column carries the bonded-exclusion scales.  Chunk size 5x10^5 elements on CPU (cache-resident, measured 4x faster than large chunks) and 6x10^7 on GPU (launch-bound); float32 on GPU. |
| `RISModel.sample` backward pass | 10^4-10^6 chains | Per bond: one (M, S) softmax + inverse-CDF draw; the sequential part is the shared forward pass |
| `build_backbone` | chains | NeRF step vectorised over the batch |
| `SimpleFF.energy_batch` / `fit_ris` | ~3k oligomer conformers | pair distances (M, n_pairs) with cached topology; an MLIP evaluates the same batch natively |

Expected effect: these kernels are memory-bandwidth bound element-wise work,
so the packing screen that takes 10-60 s per conformation on this CPU becomes
sub-second on a data-centre GPU, at which point the coarse screen can afford
10^5-10^6 cells and a finer grid, and the polishing stage mostly disappears.
The 10^5-chain amorphous ensemble similarly drops from seconds to well under a
second.  These are projections.  The container this was built in has no GPU
(no CUDA driver or device), so the CuPy path has been reviewed but not
executed.  Validating it on a GPU machine (a RunPod pod, for example) is one
command:

    bash examples/runpod_gpu_bench.sh

which installs the package and CuPy, runs `tests/test_gpu.py` (every batched
kernel compared with its NumPy result), the full test suite on the CuPy
backend, and `examples/benchmark.py` on both backends for a side-by-side
timing table.  Candidates are independent of each other, so multi-GPU
parallelism is a loop over candidates (not implemented).

## 5. Validation

Every experimental number quoted in this section is sourced, claim by claim,
in `docs/REFERENCES.md`, together with a verdict on whether it survived the
check.  Read that file before building an argument on any of them.

### 5.1 Chain repeats from the screw decomposition (rigid textbook geometry)

| Chain | polyfind c (A) | experiment (A) |
|---|---|---|
| PE all-trans (2/1) | 2.55 | 2.55 |
| PVDF beta TT | 2.58 | 2.56 |
| PVDF alpha TG+TG- | 4.56 | 4.62 |
| PVDF gamma TTTG+TTTG- | 9.11 | 9.20 |
| PE TG (isotactic-polypropylene-type 3/1 helix) | 6.35 | 6.50 (iPP) |

One note on the "experiment" column: the last row compares a hypothetical
polyethylene helix with a *different* polymer, whose backbone angles are wider.
6.50 A is the confirmed chain-axis repeat of alpha-iPP, but the row tests the
screw decomposition rather than predicting a structure.  Everything else in the
column is sourced and confirmed, including the PE 2.55 A, which is the modern
room-temperature value (Bunn's older 2.534 A is the outlier, not the standard).
See `docs/REFERENCES.md` sections 1.4 and 8.

### 5.2 Packing the known chain conformations (built-in potential, ideal angles)

`examples/pvdf_polymorphs.py`, Part A.  Axes are listed as sorted pairs since
the search does not know which is a and which is b.

| Chain | predicted a x b x c (A), density | orientation | experiment (A), density |
|---|---|---|---|
| PE T | 4.61 x 7.35 x 2.55, 1.08 | degenerate | 4.95 x 7.42 x 2.55, 1.00 |
| PVDF beta TT | 4.65 x 8.61 x 2.58, 2.06 | degenerate | 4.91 x 8.58 x 2.56, 1.97 |
| PVDF alpha TG+TG- | 5.09 x 9.14 x 4.56, 2.01 | **antiparallel** | 4.96 x 9.64 x 4.62, 1.92 |
| PVDF gamma TTTG+TTTG- | 5.40 x 8.92 x 9.11, 1.94 | **parallel** | 4.96 x 9.67 x 9.20, 1.94 |

Cell edges are within about 7% and densities within about 5% with a potential
that was never fitted to any of this.  For PE the herringbone arrangement
(setting angles +/-48 deg from a, chain 2 offset by c/2) is found as the second
minimum, 0.01 kcal/mol per CH2 above a parallel arrangement, i.e. within the
potential's accuracy.  Two qualifications on that parenthesis.  The setting angle
has no single accepted value: measured from a - and naming the axis matters,
because much of the polyethylene literature measures from b instead, and near
45 deg the two conventions are indistinguishable - published determinations run
41 to 48.8 deg, so 48 deg is inside the range and near its room-temperature end.
The "chain 2 offset by c/2", on the other hand, is not a prediction: the accepted
Pnam structure puts both chains' carbons at the same two z levels, with no
stagger, and an all-trans chain's own 2_1 screw makes a c/2 shift identical to a
180-degree rotation about the chain axis, so the offset is degenerate with the
setting angle in this parametrisation rather than a feature of the crystal.  See
`docs/REFERENCES.md` section 1.5.

The orientation column is now meaningful, which it was not before the defect of
section 5.6 was fixed, and it is worth reading carefully because two of the four
entries are not a prediction at all.  For PE and for beta the two orientations
are *exactly* degenerate, to the last digit.  That is a symmetry, not a
coincidence: both chains are mirror-symmetric about their own axis, so flipping
one is the same as rotating it, and the setting angle already covers rotations.
The flag the search reports for them is arbitrary.  Note also that
crystallographic polarity is not the same question as this flag: two chains
pointing the same way can still oppose their transverse dipoles through their
setting angles.

For the two chains where the flip is a real degree of freedom, alpha prefers
antiparallel by 0.18 kcal/mol per monomer and gamma prefers parallel by 0.78.
Neither answer was reachable before the fix.

**A correction to what that means for alpha.**  An earlier version of this
paragraph read the antiparallel preference as reproducing the antipolar alpha
phase.  Measuring the cell dipole directly (section 5.5) shows it does not.
Flipping chain 2 reverses its *axial* dipole, and those do cancel exactly, but
the two chains' *transverse* dipoles still align, leaving a net polarization of
0.078 C/m^2.  The genuinely antipolar arrangement is reachable in this
parametrisation - flip with equal setting angles cancels the transverse part
too, and a constrained polish reaches zero - but it costs 1.76 kcal/mol per
monomer more, so the illustrative potential prefers the polar one.  The chain
orientation flag and crystallographic polarity are different questions, and
only the second is the one the alpha phase is named for.

Continuous refinement then does what it is meant to: alpha relaxes to c = 4.70
A with gauche angles at the potential's own minimum (+/-80 deg); gamma relaxes
to c = 9.27 A (experiment 9.20) with deflected trans angles (171-189 deg) and
reduced gauche (+/-64 deg).  An earlier version of this sentence added "that the
real gamma chain has"; that is withdrawn, because no published torsion set for
the gamma chain was found to support it - the one crystal-structure figure
available, a DFT internal rotation of 59.2 deg, does not
(`docs/REFERENCES.md` section 2.5).

### 5.3 One full funnel run (PVDF, built-in potential, third-order RIS)

`examples/pvdf_polymorphs.py`, Part B, on 4 CPU cores.  Stage wall times:
fit 7.6 s (2,881 single points), enumeration 3.6 s (120 distinct periodic
conformations up to period 8), packing of 4 conformations 81 s, refinement
88 s, amorphous statistics 0.4 s.  Packing and refinement are where a GPU
matters; everything else is already seconds.

The fitted model moves the gauche state to +/-80 deg, at which the TG+ chain
becomes an 18/5 helix and gamma is no longer commensurate within 8 periods, so
neither is packed; three candidates with 72-144 atoms per chain repeat are
listed but skipped as too large for the packing kernel.  Lattice energies per
monomer relative to beta: T3G+TG- +2.8, alpha +3.9, TG+TG+TG-TG- +4.4.  All
four came out as parallel packings, but see section 5.6: antiparallel was
unreachable when this run was made, so that is an artifact, not a result.  Refinement changes the
alpha energy by -0.12 kcal/mol and the 8-bond glide chain by -0.45 kcal/mol
per monomer, with commensurability residuals below 0.4 deg.

The melt-like ensemble at 450 K has 51% trans bonds, a mean trans run of 1.2
bonds and a TTTT tetrad fraction of 0.002; chains leaving an all-trans stem
advance 9.3 A on average along the stem axis, turn back with probability 0.32,
and show no trans excess in their first free bond.  These numbers are what the
funnel produces for *this* potential; with a properly fitted one the same run
costs the same three minutes.

### 5.4 What the illustrative potential gets wrong

`SimpleFF` (UFF Lennard-Jones, unscreened point charges, one Fourier torsion)
is there so the pipeline runs and can be tested; its energy *differences* are
not quantitative:

* it puts beta 3.9 kcal/mol per monomer below alpha in the crystal (the real
  ordering is nearly degenerate, alpha slightly favoured: across five exchange-
  correlation functionals and four independent studies beta sits 2.6 to 6.5
  kJ/mol per monomer *above* alpha, i.e. 0.6 to 1.6 kcal/mol, with delta
  essentially degenerate with alpha and gamma in between -
  `docs/REFERENCES.md` section 5), because unscreened
  dipole alignment in the polar beta cell is over-rewarded;
* its isolated-chain RIS ranking prefers the TG+ 3/1-type helix, which PVDF
  does not form;
* it puts the polar beta cell too low for the same reason.

**The screening diagnosis above is probably wrong, and a fit says so.**  Both
bullets blame unscreened electrostatics over-rewarding dipole alignment.  If
that were the mechanism, fitting the potential against the crystal data should
have pulled the electrostatics *down*.  It did the opposite: the fit chose a
stronger net charge-to-permittivity ratio, 1.42 against 1.0, and alpha stayed
polar at every point in the search box.  So whatever makes alpha come out polar
here is not simply a missing screening factor, and section 5.5's explanation
should be treated as unproven.  A likelier candidate is the absent
depolarisation energy noted in section 6, which no adjustment of the charges
can supply.

The first two are properties of the potential, not of the search, and the
design routes a better potential to exactly the two places that fix them: the
RIS fit and the final re-scoring.  A third claim once stood here, that the
search predicts the alpha packing as polar rather than antipolar; it has been
withdrawn, because the comparison it rested on was impossible to make at the
time (section 5.6).

### 5.5 Polarization, and a field term

The cell dipole is well defined here only because every repeat unit is built
neutral, so the sum of charge times position does not move with the origin;
`CrystalPacker` asserts that.  Spontaneous polarizations of the packed
reference chains, with the illustrative potential:

| Chain | \|P\| (C/m^2) | direction | expectation |
|---|---|---|---|
| PE all-trans | 0.0000 | - | exactly zero, and for every configuration, not only the minimum |
| PVDF beta TT | 0.1405 | perpendicular to c | 0.13 rigid-dipole, 0.176-0.188 from DFT, 0.050-0.100 measured |
| PVDF alpha TG+TG- | 0.0777 | perpendicular to c | should be near zero: **disagrees** |
| PVDF gamma T3GT3G' | 0.0910 | 34 deg out of ab | polar, weaker than beta (DFT 0.071) |

**A correction to the beta target.**  This table and sections 5.4 and 5.7 used
to call 0.13 C/m^2 the *experimental* polarization of beta-PVDF.  It is not a
measurement: it is the oldest and crudest calculated estimate, the sum of rigid
monomer dipoles over the cell volume.  Modern Berry-phase and Wannier DFT put a
perfect beta crystal at 0.176-0.188 C/m^2, some 35-45% higher, while *measured*
remanent polarizations of poled PVDF films run 0.050-0.100 C/m^2 (up to 0.140
for a biaxially oriented pure-beta film).  0.13 sits between the two and belongs
to neither.  See `docs/REFERENCES.md` section 3 for the full ladder of estimates,
the measurements, and their sources.  Nothing computed changes - 0.13 only ever
appeared as a comparison target in prose and as one fit target in section 5.7 -
but the agreement claimed below is with a calculation, and not a good one.

Three of the four are right in sign and rough magnitude, and beta landing near
the rigid-dipole estimate is better than this potential deserves.  Alpha is
wrong, and
the reason is already documented: unscreened Coulomb over-rewards dipole
alignment (section 5.4), and the damped-shifted-force sum carries no
depolarisation energy (section 6).  It is a property of the charges, not of the
search, since the antipolar cell is reachable and merely scores worse.

A uniform field enters as minus the dipole dotted with the field, computed per
configuration from the placed coordinates so it follows the setting angles, the
flip and the torsions, and `pack` and `refine_crystal` both optimise in it.
Beta responds linearly along its own polar axis and does not move, being already
saturated; a field opposing it finds the 180-degree-rotated cell at equal energy,
which is polarization reversal.  A transverse field does move the structure: at
0.5 V/A the setting angles rotate by 6.6 degrees and b opens from 8.61 to 8.85 A.

### 5.6 A defect that invalidated every polarity result

Chain 2 of a two-chain cell is made antiparallel by mirroring it to
(x, -y, -z).  That mirror maps the chain's z-image k onto image -k, so the
bonded-exclusion matrices, which are indexed by the unmirrored k, stopped
excluding the flipped chain's own 1-2 and 1-3 pairs.  Those pairs sit 1.1 to
1.5 A apart, so their Lennard-Jones repulsion landed in the total: a constant
of about 3,200 kcal/mol per cell for PE and 6,290 for the PVDF chains, added
to every antiparallel configuration and to nothing else.

The test that exposes it is one line of physics: two chains 60 A apart cannot
interact, so flipping one must cost nothing.  It cost 3,155 kcal/mol.

Consequence: the antipolar half of the search space was unreachable in every
run recorded above, which is why every packed result is parallel.  Nothing
else was affected, since parallel configurations never touched the faulty
term, and all previously reported energies for them are unchanged to nine
decimal places.

The fix follows from the same symmetry: a flip is an isometry of the chain, so
its intra-chain sum cannot depend on orientation.  It is now computed once per
chain and added as a constant, which also removes the exclusion bookkeeping
from the per-configuration kernel entirely.  Regression tests assert that a
flip beyond the cutoff is free and that an isolated pair of chains costs
exactly twice one chain.

With the antiparallel branch reachable, the answer changes: **alpha-PVDF packs
antiparallel**, and wins outright at 0.62 against 0.80 kcal/mol per monomer for
the best parallel cell, at the experimental antipolar cell.  That is the
correct result for alpha-PVDF, and it was unreachable before.  For gamma the
parallel packing still wins.  For PE and beta the two orientations are exactly
degenerate at the minimum, which is not a coincidence: both chains are
mirror-symmetric about their own axis, so for them a flip is a rotation.

### 5.7 Fitting the potential to crystal data, and why it did not work

The search is now fast enough to sit inside a parameter fit, which it was not
before: one pack-and-refine of a reference chain is seconds, so a few hundred
objective evaluations are affordable.  Five parameters were fitted - relative
permittivity, a charge scale, the hydrogen and fluorine Lennard-Jones radii,
and the leading torsion coefficient - against PE, beta and alpha cell edges,
densities, the beta polarization and the alpha/beta energy gap, with **gamma
held out entirely**.  Eighty-six evaluations, twenty-five minutes.

Two of the five parameters turned out to be structurally unidentifiable before
any data was involved.  At ideal torsion angles the third Fourier coefficient
vanishes identically and only the sum of the first two enters, so they cannot
be separated; and the permittivity and the charge scale enter every energy only
through the ratio of charge squared to permittivity, so they are exactly
degenerate there and separable only by the polarizations.

The fit halved its objective, from 100.8 to 56.5, and **it does not
generalise**.  Nearly all of the gain sits in the three targets that each have
their own private knob: PE's a axis is set by the hydrogen radius, alpha's b
axis by the fluorine radius, and the alpha-beta energy gap by the torsion
coefficient.  Targets without a dedicated knob got *worse*, including beta's
polarization, which was a fitted target and drifted from 0.140 to 0.160 against
the 0.13 target it was fitted to.  (That 0.13 is the rigid-dipole estimate, not
an experimental value; see section 5.5 and `docs/REFERENCES.md` section 3.
Against the DFT value of 0.176-0.188 the drift was towards the truth, which
makes this a weaker piece of evidence than it looked, not a stronger one.)

An ablation settles which parameters are real.  The two Lennard-Jones radii
alone deliver essentially the whole held-out improvement (gamma root-mean-square
error 3.73% against 3.48% for the full fit) while keeping beta's polarization at
0.128.  The other three alone cut the fitted objective to 89.0 while making the
held-out gamma *worse than not fitting at all*, 5.96% against 5.80%.  So two
parameters transfer and three fit noise, on twelve independent observations.

The conclusion is about the data, not the optimiser: **crystal structures alone
cannot determine this potential**.  Twelve numbers, several of them not
independent, cannot pin five parameters, and the quantities that would
discriminate - torsion profiles, and the relative energies of conformers - are
simply not in the training set.  The isolated-chain conformational ranking,
which the objective never constrained, degraded further during the fit.

What follows is a sharpening of section 2.5's decision rather than a reversal of
it.  That decision was not to make a machine-learned potential a *runtime
dependency*, and that stands.  Using a high-quality method once, offline, to
generate a few hundred reference torsion energies to fit a classical potential
against is a different activity: it is ordinary force-field development, the
expensive method appears as a data source rather than as the engine, and the
fitted result still evaluates in microseconds.  That is the missing ingredient,
and this fit is the evidence that crystal data by itself will not substitute for
it.

The fitted parameters ship as an opt-in preset carrying this warning.  The
defaults are unchanged.

### 5.8 A withdrawn claim about beta-PVDF's geometry

This document twice asserted that beta-PVDF's experimental geometry is
unreachable under rigid bonds, because unequal backbone angles near 112 and 118
degrees combined with torsions near +/-172 degrees cannot give a straight chain.
Checked against the literature, the premise was wrong and the claim is
withdrawn.

The accepted Form I structure has **equal** backbone angles at both carbons and
torsions of exactly 180 degrees.  A density-functional study of the crystal
gives a C-C-C angle of 114.4 degrees at both backbone carbons and an internal
rotation angle of 180 degrees, against 112 degrees and a 2.534 A repeat for
polyethylene, the difference attributed to repulsion between fluorines on
neighbouring carbons.  So the equal 114-degree angles and exact trans this
package uses are not a compromise at all; they are the structure.

The alternately-deflected zigzag, in which the CF2 carbons alternate out of the
plane to relieve crowding between fluorines 2.56 A apart against a van der Waals
contact of 2.70 A, is a proposal rather than a refinement result.  The
deflection would double the chain repeat to 5.12 A, and the intermediate layer
line that doubling requires is absent from the diffraction pattern; the original
fit used a statistically disordered structure at 2.56 A instead.  Treating it as
the experimental geometry, as this document did, was a misreading.

What survives is smaller and concrete.  The C-C distance in the DFT structure is
1.528 A rather than the 1.54 A used here, and adopting it moves the computed
chain repeat from 2.583 to 2.563 A against an experimental 2.56.  That change is
queued rather than applied: it shifts every PVDF energy and would invalidate the
tables above, and the unfitted potential contributes cell errors of 3 to 6%,
which dwarf this 0.9%.  It belongs with the next full re-measurement.

The general lesson is the same one section 5.7 records about the fit.  Both this
and the "unreachable geometry" argument were careful reasoning from a number
written in a code comment that nobody had checked.  Published structures for
these polymers exist and are not hard to consult.

That consultation has since been done for every load-bearing number in this
document, and the result is `docs/REFERENCES.md`: one entry per claim, with the
source and a verdict of confirmed, corrected, unverified or wrong.  Every
statement in this section 5.8 is confirmed there (sections 2.2 and 2.3).

### 5.9 Refitting against first-principles data, and where the wall moved to

Section 5.7's fit failed for a diagnosable reason: crystal data does not contain
the quantities that discriminate.  Those quantities existed in a sibling
project - 566 frames of PBE-D3 energies and forces from relaxed torsion scans
and conformer sampling on five-monomer oligomers across 37 systems, every one of
them a VDF backbone carrying a single substituted unit.  Refitting against them,
with ten chemistries held out **whole** rather than frames held out at random:

| | section 5.7 (crystal data) | this fit (DFT data) |
|---|---|---|
| observations | 12, several dependent | 467 frames, 31 systems |
| held out | one polymorph | ten entire chemistries |
| held-out error | worse than not fitting | 3.32 -> 1.76 kcal/mol |
| parameters that transfer | 2 of 5 | all 27 groups contribute |

So the data diagnosis was right, and the fit now generalises.  The headline
physical result is one the objective never saw: **alpha now falls below beta by
4.54 kJ/mol per monomer**, inside the 2.6 to 6.5 range that four independent
studies agree on, where the illustrative potential had beta below alpha by 7.38.
That is the polymorph ordering this package exists to get right.

Two other acceptance tests still fail, and the way they fail is the finding.
The isolated-chain ranking still prefers a 3/1 helix PVDF does not form, and
alpha's cell still packs polar rather than antipolar.  An ablation shows why
more fitting will not help: **different parameter groups carry different tests
and no setting carries both.**  Lennard-Jones parameters alone fix the chain
ranking and halve the polarity error, but only the full parameter vector fixes
the alpha/beta ordering, and doing that puts the wrong helix back on top.  The
three targets are in tension inside this functional form.

Two independent diagnostics point the same way.  Five of the 27 parameters sit
on their bounds, and released from those bounds the fit drives **fluorine's
partial charge positive** - physically absurd for the most electronegative
element, and a clear sign that atom-centred charges are being used to do work
the repulsion and dispersion terms should be doing.  Separately, the forces
contributed almost nothing, 0.3% of the held-out improvement, because 93% of the
force signal is bond stretching, which a rigid-geometry potential cannot
represent at all.

**The binding constraint has therefore moved from the data to the functional
form.**  That is progress, and it points somewhere specific.  A form that could
satisfy these targets together would need at least: explicit bond and angle
terms, which would make the force data usable instead of 93% wasted and would
also remove the strained-reference problem section 5.8 and the chemistry
extension both ran into; and electrostatics richer than fixed atom-centred point
charges, since the carbon-fluorine dipole is exactly what such charges
represent worst.  Both are ordinary force-field technology, and both are larger
changes than a refit.

The fitted parameters ship as the preset `pvdf-dft-fit`.  `SimpleFF()`'s
defaults are unchanged and still the default, so nothing silently moves.

## 6. Limitations and roadmap

* RIS uses rigid bond geometry and discrete states; the continuous refinement
  stage relaxes torsions but not bond angles.  Adding the two backbone angles
  per repeat to the refinement variables is straightforward.
* Third-order terms are fitted by inclusion-exclusion at the state angles, not
  from a full 3D scan.
* Thermodynamic ranking says nothing about which polymorph forms kinetically;
  the amorphous-ensemble statistics (trans-run supply) are the nucleation-side
  proxy provided here.
* Head-to-head / tail-to-tail defects (4 to 6% in typical PVDF, with 5 to 10%
  quoted across the literature) need an extra bond type in the model; the code
  supports arbitrary bond-type periods but no defect placement yet.  Whether
  they *stabilise* beta, as this bullet used to assert flatly, is disputed in
  the sources: see `docs/REFERENCES.md` section 7.1.
* The lattice energy uses a damped-shifted-force Coulomb sum, not Ewald; fine
  for ranking neutral chains, not for absolute lattice energies or the
  depolarisation energy of polar cells.
* One and two chains per cell with chain 2 at (1/2, 1/2) are implemented;
  general chain positions and more chains per cell are a small extension of
  the parameter vector.
* Free energies of crystals (phonons, thermal expansion) are not computed.
* The CuPy backend is untested on hardware.

## 7. Package layout

```
src/polyfind/
  polymers.py     monomer templates (PVDF, PE), RIS state sets, UFF nonbonded parameters
  ris.py          transfer matrices; Viterbi, k-best, cyclic k-best, partition function,
                  marginals, exact batched sampling; third order by state augmentation
  chain.py        torsion sequence -> Cartesian coordinates (NeRF), batched
  helix.py        screw decomposition, helix descriptors, symmetry canonicalisation
  enumerate.py    periodic candidate generation, dedupe, ranking, known-polymorph labels
  forcefield.py   SimpleFF, ASECalculator adapter, fit_ris (1D/2D scans + triplet corrections)
  pack.py         periodic chain construction, batched lattice-energy kernel, search, CIF
  refine.py       torsions + cell refinement with a commensurability penalty
  amorphous.py    Boltzmann ensembles, run statistics, lamella-interface sampling
  backend.py      NumPy / CuPy selection for the batched kernels
  pipeline.py     the funnel and the report
  cli.py          polyfind fit | enumerate | pack | sample | pipeline
tests/            290 tests: every DP routine vs brute force, geometry, helices, packing, pipeline
examples/         pvdf_polymorphs.py reproduces Section 5
```
