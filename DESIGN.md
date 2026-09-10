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
the test suite (362 tests, 5 skipped).

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

No potential enters this table; it is the declared bond lengths and angles put
through the screw decomposition.  Re-measured with PVDF's C-C distance at the
DFT Form I value of 1.528 A (it was the textbook 1.54 A until the batch of
corrections recorded in `docs/REFERENCES.md` was applied):

| Chain | polyfind c (A) | experiment (A) | error | error at C-C = 1.54 A |
|---|---|---|---|---|
| PE all-trans (2/1) | 2.553 | 2.55 | +0.13% | +0.13% (PE unchanged) |
| PVDF beta TT | 2.563 | 2.56 | **+0.12%** | +0.90% |
| PVDF alpha TG+TG- | 4.520 | 4.62 | **-2.17%** | -1.40% |
| PVDF gamma TTTG+TTTG- | 9.040 | 9.20 | **-1.74%** | -0.97% |
| PE TG (isotactic-polypropylene-type 3/1 helix) | 6.351 | 6.50 | -2.30% | -2.30% (PE unchanged) |

**The correction cuts beta's error by a factor of seven and makes alpha's and
gamma's worse, and that has to be said plainly.**  1.528 A is the C-C distance
of the *Form I* (beta) crystal, so beta is the row it was derived from, and
there it lands almost exactly on the experimental 2.56 A.  Every repeat is
proportional to the bond length, so shortening it shortens all of them; alpha
and gamma were already under-predicted, and they move further down.  The
root-mean-square error over the three PVDF rows therefore rises from 1.11% to
1.61% even though the one row with a matching reference structure improves.
Both facts are consequences of a single rigid bond length being asked to serve
three conformations whose measured C-C distances are not in fact identical.

One note on the "experiment" column: the last row compares a hypothetical
polyethylene helix with a *different* polymer, whose backbone angles are wider.
6.50 A is the confirmed chain-axis repeat of alpha-iPP, but the row tests the
screw decomposition rather than predicting a structure.  Everything else in the
column is sourced and confirmed, including the PE 2.55 A, which is the modern
room-temperature value (Bunn's older 2.534 A is the outlier, not the standard).
See `docs/REFERENCES.md` sections 1.4 and 8.

### 5.2 Packing the known chain conformations (built-in potential, ideal angles)

`examples/pvdf_polymorphs.py`, Part A, with the default `SimpleFF()`
(illustrative, unfitted).  Axes are listed as sorted pairs since the search does
not know which is a and which is b.  Re-measured after the geometry corrections
of `docs/REFERENCES.md`; the previous numbers are in the last column for
comparison.

| Chain | predicted a x b x c (A), density | orientation | experiment (A), density | before the corrections |
|---|---|---|---|---|
| PE T | 4.61 x 7.35 x 2.55, 1.076 | degenerate | 4.95 x 7.42 x 2.55, 1.00 | unchanged |
| PVDF beta TT | 4.64 x 8.62 x 2.56, 2.073 | degenerate | 4.91 x 8.58 x 2.56, 1.97 | 4.65 x 8.61 x 2.58, 2.058 |
| PVDF alpha TG+TG- | 5.09 x 9.11 x 4.52, 2.027 | **antiparallel** | 4.96 x 9.64 x 4.62, 1.92 | 5.09 x 9.14 x 4.56, 2.008 |
| PVDF gamma TTTG+TTTG- | 5.39 x 8.91 x 9.04, 1.960 | **parallel** | 4.96 x 9.67 x 9.20, 1.93 | 5.40 x 8.92 x 9.11, 1.939 |

The experimental gamma density reads 1.93, not the 1.94 an earlier version of
this table carried: 1.93 is what the tabulated gamma cell actually gives, while
1.94 came from a different determination (`docs/REFERENCES.md` section 1.3).
That correction makes the *apparent* agreement on gamma's density worse rather
than better - the predicted 1.960 is +1.5% against 1.93, where 1.939 against
1.94 looked like -0.1% - and the looking-good was an accident of comparing a
density from one structure with a cell from another.  The gamma reference cell is
also monoclinic, with beta = 93 deg, and **that angle is not expressible in this
packing model at all**: the chain axis is the cell's c by construction and a and
b are built perpendicular to it, so only the in-plane a-to-b angle is a search
variable (`pack(..., gamma_free=True)`).  The package therefore approximates
gamma as orthorhombic.  The cost is small and measurable: the monoclinic cell's
volume differs by sin(93 deg) = 0.14%, so the density the model could at best
reproduce differs by that much, against cell-edge errors of 3 to 9% here.  See
`pipeline.REFERENCE_CELL_ANGLES`, which records the angles and this limitation.

Cell edges are within about 9% and densities within about 5% with a potential
that was never fitted to any of this.  Every density moved *up* with the shorter
C-C bond, so beta's density error grows from +4.5% to +5.3% and alpha's from
+4.6% to +5.6% while their edges barely move; beta's c error falls from +0.90% to
+0.12%.  For PE the herringbone arrangement
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
are *exactly* degenerate, to the last digit, and the reason is sharper than an
earlier version of this paragraph claimed.

Measured: applying the flip to a planar-zigzag chain maps its atom set onto
itself, residual 3e-15 A, with **no** accompanying rotation or z-shift at all.
The flip is simply a symmetry operation of that chain.  It has to be: a planar
zigzag lies in a plane, so reflecting across that plane does nothing to it, and
the remaining reversal along the chain is absorbed by its own two-fold screw
axis.  For alpha and gamma the flip is not a symmetry - the best match over all
rotations and z-shifts leaves 0.8 A and 3.1 A - so there it is a real degree of
freedom.

This paragraph previously said the chains were mirror-symmetric about their axis
so that a flip was equivalent to a rotation absorbed by the setting angle.  The
intuition about the mirror plane was right and the mechanism was wrong; no
rotation is involved.  Note separately that crystallographic polarity is not the
same question as this flag: two chains pointing the same way can still oppose
their transverse dipoles through their setting angles.

For the two chains where the flip is a real degree of freedom, alpha prefers
antiparallel by 0.18 kcal/mol per monomer and gamma prefers parallel by 0.79.
Neither answer was reachable before the fix.  Both margins are unchanged by the
geometry corrections (they were 0.18 and 0.78), which is worth noting: the
orientation preferences are a property of the packing and not of the 0.8% of
bond length that moved.

**A correction to what that means for alpha.**  An earlier version of this
paragraph read the antiparallel preference as reproducing the antipolar alpha
phase.  Measuring the cell dipole directly (section 5.5) shows it does not.
Flipping chain 2 reverses its *axial* dipole, and those do cancel exactly, but
the two chains' *transverse* dipoles still align, leaving a net polarization of
0.079 C/m^2.  The genuinely antipolar arrangement is reachable in this
parametrisation - flip with equal setting angles cancels the transverse part
too, and a constrained polish (`fitting.antipolar_cell`) reaches exactly zero -
but it costs **0.20 kcal/mol per monomer** more after refinement, so the
illustrative potential prefers the polar one.  (An earlier version of this
sentence quoted 1.76 kcal/mol from a different and unrecorded measurement at the
rigid ideal-angle cell; 0.20 is the figure `fitting.acceptance_tests` reports and
is reproducible.  It is unmoved by the geometry corrections: +0.200 before and
after.)  The chain orientation flag and crystallographic polarity are different
questions, and only the second is the one the alpha phase is named for.

Continuous refinement then does what it is meant to: alpha relaxes to c = 4.80 A
with gauche angles at +/-65.7 deg; gamma relaxes to c = 9.51 A (experiment 9.20)
with deflected trans angles (173.8-186.2 deg) and gauche at +/-65.5 deg.  Two
notes on that sentence.  First, it previously read "c = 4.70 A ... (+/-80 deg)"
and "c = 9.27 A ... (171-189 deg) ... (+/-64 deg)", and those numbers had already
gone stale before these corrections - re-measuring the unchanged code gave 4.83
and 9.57 A - because the refinement stage acquired line-group constraints since
they were recorded.  The corrections then moved them down to 4.80 and 9.51, which
for gamma is an *improvement* against the experimental 9.20 (error +3.98% to
+3.38%) and for alpha likewise (+4.52% to +3.96%): refinement over-expands both
cells, and a shorter bond partly cancels that.  Second, an earlier version added
"that the real gamma chain has"; that is withdrawn, because no published torsion
set for the gamma chain was found to support it - the one crystal-structure
figure available, a DFT internal rotation of 59.2 deg, does not
(`docs/REFERENCES.md` section 2.5).

### 5.3 One full funnel run (PVDF, built-in potential, third-order RIS)

`examples/pvdf_polymorphs.py`, Part B, on 6 physical cores, re-measured after the
geometry corrections.  Stage wall times: fit 0.3 s (2,881 single points),
enumeration 4.5 s (120 distinct periodic conformations up to period 8), packing
and refinement of 4 conformations 49 s cold / 17 s with a warm interaction-table
cache, amorphous statistics 0.4 s.  Packing and refinement are where a GPU
matters; everything else is already seconds.

Three of those stage times had **already** gone stale before these corrections,
for reasons of their own, and re-measuring the unchanged code gives the same
numbers to within run-to-run noise (fit 0.4 s, enumeration 5.0 s, packing and
refinement 53 s cold / 14 s warm).  The fit stage is some twenty times faster
than the 7.6 s once recorded here; packing and refinement are fused per
candidate, so they no longer have separate wall times; and the exhaustive table
screen with its cache replaced random sampling.  Whatever produced those changes,
it was not these corrections, which move no stage time outside run-to-run
spread.

The fitted model moves the gauche state to +/-80 deg, at which the TG+ chain
becomes a helix with more than 8 periods per repeat and gamma is no longer
commensurate within 8 periods, so neither is packed; three candidates with 72 to
144 atoms per chain repeat are listed but skipped as too large for the packing
kernel.  What is packed is beta, alpha and two further glide chains.  Lattice
energies per monomer relative to beta, after refinement: alpha +1.13,
TTTG+TG- +1.41, TG+TG+TG-TG- +2.87 (they were +1.21, +1.46 and +2.92 before the
corrections).  An earlier version of this paragraph recorded +2.8, +3.9 and +4.4
from a run in which antiparallel packings were unreachable (section 5.6) and the
refinement had no line-group constraints; those numbers are superseded.  Alpha and
the 8-bond glide chain now come out antiparallel and beta and TTTG+TG- parallel,
which the earlier run could not have found at all; refinement lowers the alpha
energy by 3.67 kcal/mol per monomer and the 8-bond glide chain by 3.08, with
commensurability residuals below 0.2 deg.

The melt-like ensemble at 450 K has 52% trans bonds, a mean trans run of 1.2
bonds and a TTTT tetrad fraction of 0.002; chains leaving an all-trans stem
advance 8.8 A on average along the stem axis, turn back with probability 0.33,
and show a marked *deficit* of trans in their first free bond (18% against 51%
in bulk).  Before the corrections those read 51%, 1.2 bonds, 0.002, 9.0 A, 0.33
and 17% against 51% - so this paragraph's earlier "advance 9.3 A ... and show no
trans excess in their first free bond" was stale on the first number and too weak
on the last.  These numbers are what the funnel produces for *this* potential;
with a properly fitted one the same run costs the same minute.

### 5.4 What the illustrative potential gets wrong

`SimpleFF` (UFF Lennard-Jones, unscreened point charges, one Fourier torsion)
is there so the pipeline runs and can be tested; its energy *differences* are
not quantitative:

* it puts beta 1.8 kcal/mol per monomer below alpha in the crystal, measured as
  the refined alpha/beta gap of `fitting.acceptance_tests` (+7.34 kJ/mol per
  monomer; it was +7.38 before the geometry corrections, and an earlier version
  of this bullet quoted 3.9 kcal/mol from the superseded funnel run of section
  5.3).  The real ordering is nearly degenerate, alpha slightly favoured: across
  five exchange-correlation functionals and four independent studies beta sits
  2.6 to 6.5 kJ/mol per monomer *above* alpha, i.e. 0.6 to 1.6 kcal/mol, with
  delta essentially degenerate with alpha and gamma in between
  (`docs/REFERENCES.md` section 5).  The stated reason is that unscreened
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

| Chain | \|P\| (C/m^2) | before the corrections | direction | expectation |
|---|---|---|---|---|
| PE all-trans | 0.0000 | 0.0000 | - | exactly zero, and for every configuration, not only the minimum |
| PVDF beta TT | 0.1416 | 0.1405 | perpendicular to c | 0.13 rigid-dipole, 0.176-0.188 from DFT, 0.050-0.100 measured |
| PVDF alpha TG+TG- | 0.0785 | 0.0777 | perpendicular to c | should be near zero: **disagrees** |
| PVDF gamma T3GT3G' | 0.0920 | 0.0910 | 34 deg out of ab | polar, weaker than beta (DFT 0.071) |

Every PVDF polarization rose by about 0.8%, which is the shorter C-C bond packing
the same charges into a slightly smaller cell and nothing more.  Which way that
counts as agreement depends on the target, and the targets disagree by 40%: beta
moves a shade closer to the DFT 0.176-0.188 and a shade further from both the
rigid-dipole 0.13 and the measured 0.050-0.100, while gamma moves further from
the DFT 0.071.  None of the shifts is large enough to matter against that spread.
Alpha's disagreement is unaffected.

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
0.5 V/A the setting angles rotate by 6.7 degrees and b opens from 8.62 to 8.86 A
(6.6 degrees and 8.61 to 8.85 before the geometry corrections).  That is the
*local* response - the zero-field minimum polished in the field.  Re-screening
the whole cell in the same field instead finds a different basin altogether, with
the setting angles turned by 60 degrees, a opened to 5.01 and b closed to 8.03 A,
which is worth knowing before quoting the local number as "the" response.

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
antiparallel**, and wins outright at 2.05 against 2.23 kcal/mol per monomer for
the best parallel cell, at the experimental antipolar cell.  (Re-measured with
the corrected geometry; the same two numbers were 0.62 and 0.80 before it, so
what matters - the 0.18 kcal/mol margin and its sign - is unchanged, while the
absolute lattice energies of the rigid ideal-angle alpha cell moved by 1.4
kcal/mol per monomer because the shorter C-C bond tightens its internal
contacts.)  That is the
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
polarization, which was a fitted target and drifted from 0.141 to 0.162 against
the 0.13 target it was fitted to.  (That 0.13 is the rigid-dipole estimate, not
an experimental value; see section 5.5 and `docs/REFERENCES.md` section 3.
Against the DFT value of 0.176-0.188 the drift was towards the truth, which
makes this a weaker piece of evidence than it looked, not a stronger one.)

**What was and was not re-measured here after the geometry corrections.**  The
fit itself - 86 objective evaluations, twenty-five minutes - was run at the old
geometry and has not been repeated, so the objective values (100.8 to 56.5, and
89.0 in the ablation) and the held-out gamma percentages below are as they were.
The quantities that are a straight evaluation of the shipped preset were
re-measured and moved by about 1%: beta's polarization under `pvdf-crystal-fit`
is 0.1615 C/m^2 against 0.1604 before, and its alpha/beta gap +1.09 kJ/mol per
monomer against +1.28.  Nothing in the conclusion of this section turns on
differences that size.

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
1.528 A rather than the 1.54 A this package used, and adopting it moves the
computed beta chain repeat from 2.583 to 2.563 A against an experimental 2.56.
**That change has now been applied**, together with the rest of the batch in
`docs/REFERENCES.md`, and the tables above are the re-measurement.  Two things it
did that were not anticipated when it was queued: it made alpha's and gamma's
chain repeats *worse* (section 5.1), because 1.528 A is Form I's C-C and every
repeat scales with it; and it raised every PVDF density by about 1%, which the
unfitted potential already over-predicts.  What it improved is beta's repeat, by
a factor of seven, and - through the refinement stage, which over-expands - both
refined polymorph cells.

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
4.81 kJ/mol per monomer**, inside the 2.6 to 6.5 range that four independent
studies agree on, where the illustrative potential had beta below alpha by 7.34.
That is the polymorph ordering this package exists to get right.  (Both figures
re-measured after the geometry corrections; they were 4.54 and 7.38 before, so
the correction moved the fitted result a little further into the accepted range.
The held-out error table above is from the fit itself, which was run at the old
geometry and has not been repeated.)

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

### 5.10 Valence terms: the forces become usable, and a data problem surfaces

Section 5.9 concluded that the binding constraint had moved from the data to the
functional form, and named two gaps: no bond or angle terms, and electrostatics
too crude.  Both were addressed.  Harmonic stretching and bending were added with
fitted equilibria and stiffnesses rather than invented ones, and fluorine's charge
was allowed to sit off the nucleus while its Lennard-Jones centre stayed on it.
Held-out results, ten chemistries held out whole:

| | held-out energy | held-out force | force correlation |
|---|---|---|---|
| illustrative | 3.32 | 14.37 | +0.06 |
| DFT fit, rigid form | 1.76 | 14.06 | +0.10 |
| rigid form refit to full forces | 2.00 | 13.86 | +0.17 |
| with valence terms | **1.36** | **5.23** | **+0.93** |

The third row is the control that makes the point: refitting the *rigid* form
against the same full force data gets nowhere, 13.86 against the 14.05 that
predicting zero scores.  **The valence terms are what made the forces usable**,
and the force correlation going from +0.10 to +0.93 is the clearest single number
in this document.

Three further checks came out well.  Fluorine's partial charge is now physical
and off its bound at +0.114 e, and released from bounds entirely the new form
picks the electronegativity-correct sign for every charge increment unprompted,
where the rigid form got three of six backwards.  The strain that made all-trans
an unusable reference state largely dissolves once angles can relax: PVDC from
47 kcal/mol to 12, CFE 120 to 12, AN 27 to 2, VDCN 57 to minus 1.  And PVDC's
fitted equilibrium angle at the CH2 carbon comes out at 123.1 degrees against a
measured 123, recovered from DFT data alone with no crystallographic input - an
independent confirmation of both the valence terms and the strain diagnosis.  Its
angle at the CCl2 carbon is 102 against a measured 114, so the agreement is one
sided.

PVDC's first figure was 231 kcal/mol until its backbone angles were corrected from
114/114 to the measured 123/114 (`docs/REFERENCES.md` section 6.2).  That is worth
pausing on rather than editing quietly: the fitted valence model's own answer,
123.1 degrees, was already telling the rigid model what its input should be, and
giving the rigid model that input removed four fifths of the gap the valence terms
were being credited with closing.  The remaining 47-to-12 is what the angle terms
genuinely buy.  The relaxed geometry is unchanged either way, at 123.1 and 102.1
degrees, which is the check that this is a starting-point effect and not a
different minimum.

Three things did not work, and they matter as much.

**The off-site charge is a null result.**  Setting its displacement to zero gives
essentially the same fit.  It was the valence terms, not the electrostatics, that
fixed fluorine's charge, and the charge-permittivity degeneracy is still exact.

**The forces still do not help the energies**, 1.16 kcal/mol held out without
them against 1.36 with.  The diagnosis is a provenance problem in the reference
data rather than a modelling one: the fitted equilibrium bond lengths come out at
1.467 A for C-F against a measured 1.358, and 1.458 for C-C against 1.535.  Those
geometries are not at a PBE-D3 minimum - the reference force magnitude confirms
it, since a relaxed structure would show forces near zero and these average
14 kcal/mol per angstrom.  The dataset's geometries were optimised at one level
and labelled at another, so fitting bond parameters to its forces bends them
toward the wrong level.  Forces are kept at a low weight because they improve the
conditioning of the fit by a factor of forty and because refinement needs
gradients, not because they improve the energies.

**The acceptance tests still stand at one of three**, re-measured with the
corrected geometry:

| preset | E(alpha) - E(beta), kJ/mol per monomer | RIS top (margin, kcal/mol) | E(anti) - E(polar) | beta \|P\| |
|---|---|---|---|---|
| illustrative | +7.34 (was +7.38) | TG+ (1.58, was 1.64) | +0.200 (was +0.200) | 0.1408 |
| `pvdf-crystal-fit` | +1.09 (was +1.28) | TG+ (3.02, was 3.04) | +0.165 (was +0.166) | 0.1615 |
| `pvdf-dft-fit` | **-4.81** (was -4.54) | TG+ (0.49, was 0.34) | +0.092 (was +0.091) | 0.1148 |
| `pvdf-dft-valence` | **-3.13** (was -2.96) | TG+ (0.43, was 0.18) | +0.084 (was +0.085) | 0.1149 |

Alpha stays below beta, now at 3.13 kJ/mol per monomer, a little further inside
the accepted range than the 2.96 this fit gave at the old geometry.  The chain
ranking still prefers a helix PVDF does not form, and here the geometry
correction made things **worse**: the margin over the next candidate, which had
narrowed from 1.64 to 0.18 kcal/mol, widens back to 0.43.  Alpha's cell still
packs polar, unmoved.  So the functional form improved
substantially on every physical diagnostic while the three targets stayed in
tension, which suggests the remaining obstacle is not the terms in the potential
but something about how those two particular quantities are computed.

One thing found along the way is worth recording for its own sake.  Testing the
chain ranking exposed a spurious third-order term of about minus 50 kcal/mol,
produced by inclusion-exclusion across a steric overlap of order a million, which
would have faked a pass on that test.  The fitting code now refuses such terms.
An agent found it while checking its own acceptance test rather than reporting
the pass.

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
tests/            362 tests: every DP routine vs brute force, geometry, helices, packing, pipeline
examples/         pvdf_polymorphs.py reproduces Section 5
```
