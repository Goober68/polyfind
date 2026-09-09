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
(beta-PVDF dihedrals near +/-172 deg, alpha gauche near +/-45 deg).  The
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
calculator (a MACE or similar machine-learned potential is the intended
production path: `ASECalculator`), or with DFT; everything downstream inherits
that accuracy.  The scan is a batch of independent single points and goes
through `energy_batch`, so a GPU-resident potential evaluates it as one batch.

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
the test suite (72 tests).

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
second.  These are projections: no GPU was available in the environment where
this was built, and the CuPy path has not been executed; the kernels use only
functions CuPy provides, but the first GPU run should start with the tests
under `POLYFIND_DEVICE=cuda`.  Candidates are independent of each other, so
multi-GPU parallelism is a loop over candidates (not implemented).

## 5. Validation

### 5.1 Chain repeats from the screw decomposition (rigid textbook geometry)

| Chain | polyfind c (A) | experiment (A) |
|---|---|---|
| PE all-trans (2/1) | 2.55 | 2.55 |
| PVDF beta TT | 2.58 | 2.56 |
| PVDF alpha TG+TG- | 4.56 | 4.62 |
| PVDF gamma TTTG+TTTG- | 9.11 | 9.20 |
| PE TG (isotactic-polypropylene-type 3/1 helix) | 6.35 | 6.50 (iPP) |

### 5.2 Packing the known chain conformations (built-in potential, ideal angles)

`examples/pvdf_polymorphs.py`, Part A.  Axes are listed as sorted pairs since
the search does not know which is a and which is b.

| Chain | predicted a x b x c (A), density | experiment (A), density |
|---|---|---|
| PE T | 4.61 x 7.35 x 2.55, 1.08 | 4.95 x 7.42 x 2.55, 1.00 |
| PVDF beta TT | 4.65 x 8.61 x 2.58, 2.06 | 4.91 x 8.58 x 2.56, 1.97 |
| PVDF alpha TG+TG- | 5.08 x 9.32 x 4.56, 1.97 | 4.96 x 9.64 x 4.62, 1.92 |
| PVDF gamma TTTG+TTTG- | 5.40 x 8.92 x 9.11, 1.94 | 4.96 x 9.67 x 9.20, 1.94 |

Cell edges are within about 7% and densities within about 5% with a potential
that was never fitted to any of this.  The beta cell comes out with parallel
(polar) chains, as it should.  For PE the herringbone arrangement (setting
angles +/-48 deg from a, chain 2 offset by c/2) is found as the second minimum,
0.01 kcal/mol per CH2 above a parallel arrangement, i.e. within the potential's
accuracy.

Continuous refinement then does what it is meant to: alpha relaxes to c = 4.70
A with gauche angles at the potential's own minimum (+/-80 deg); gamma relaxes
to c = 9.27 A (experiment 9.20) with the deflected trans angles (171-189 deg)
and reduced gauche (+/-64 deg) that the real gamma chain has.

### 5.3 What the illustrative potential gets wrong

`SimpleFF` (UFF Lennard-Jones, unscreened point charges, one Fourier torsion)
is there so the pipeline runs and can be tested; its energy *differences* are
not quantitative:

* it puts beta 3.9 kcal/mol per monomer below alpha in the crystal (the real
  ordering is nearly degenerate, alpha slightly favoured), because unscreened
  dipole alignment in the polar beta cell is over-rewarded;
* its isolated-chain RIS ranking prefers the TG+ 3/1-type helix, which PVDF
  does not form;
* it predicts the alpha packing as polar rather than antipolar, the antipolar
  arrangement being 0.17 kcal/mol per monomer higher.

Every one of these is a property of the potential, not of the search, and the
design routes a better potential to exactly the two places that fix them: the
RIS fit and the final re-scoring.

## 6. Limitations and roadmap

* RIS uses rigid bond geometry and discrete states; the continuous refinement
  stage relaxes torsions but not bond angles.  Adding the two backbone angles
  per repeat to the refinement variables is straightforward.
* Third-order terms are fitted by inclusion-exclusion at the state angles, not
  from a full 3D scan.
* Thermodynamic ranking says nothing about which polymorph forms kinetically;
  the amorphous-ensemble statistics (trans-run supply) are the nucleation-side
  proxy provided here.
* Head-to-head / tail-to-tail defects (a few percent in PVDF, known to
  stabilise beta) need an extra bond type in the model; the code supports
  arbitrary bond-type periods but no defect placement yet.
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
tests/            72 tests: every DP routine vs brute force, geometry, helices, packing, pipeline
examples/         pvdf_polymorphs.py reproduces Section 5
```
