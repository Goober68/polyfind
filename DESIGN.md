# polyfind: efficient determination of stable atomic arrangements in semi-crystalline polymers

Target case: poly(vinylidene fluoride), PVDF, -(CH2-CF2)n-. About half crystalline, with
near-degenerate polymorphs (alpha = TGTG', beta = TTTT, gamma = TTTGTTTG', delta = polar alpha),
strong F...F repulsion and dipolar interactions, and an amorphous fraction whose local
conformational statistics decide which polymorph nucleates.

## Why the usual approaches are slow

| Method | Cost driver | Failure mode for PVDF |
|---|---|---|
| Molecular dynamics / annealing | 10^6-10^8 force calls; equilibration above Tg | Trapped in local minima; polymorph ranking never converges |
| Cartesian random / evolutionary CSP | Search dimension 3N | Chain connectivity is a constraint the search does not know about |
| Periodic DFT on guessed cells | Hours per structure | Cannot sample the amorphous fraction at all |

## The approach

Solve the problem at the level where it is low-dimensional: the sequence of backbone
torsional states. Every step that is normally stochastic becomes exact and discrete.

### 1. Conformation search is a dynamic program

Backbone dihedrals occupy discrete rotational isomeric states (T, G+, G-). With first- and
second-neighbour interaction energies the chain is a Markov field, so the following are all
O(N * S^2) via transfer-matrix / Viterbi recursions:

- minimum-energy conformation,
- the k best conformations,
- partition function and conformational free energy,
- exact Boltzmann samples and marginals.

For PVDF the pair matrices alternate between CH2-centred and CF2-centred pairs, which is where
the F...F repulsion and the pentane effect enter. Third-neighbour terms, if needed, are added by
treating pairs of states as the state (state augmentation), keeping the same algorithms.

### 2. Crystalline chains are periodic sequences

Polymorph chains are cycles of length P in state space (TGTG', TTTT, TTTGTTTG'). A k-best
cyclic dynamic program returns the lowest-energy cycles for each period directly. Cyclic shift
by the monomer repeat, chain reversal, and mirror (G+ <-> G-) are symmetries of the energy, so
candidates are canonicalised and deduplicated before any 3D work.

### 3. Helical symmetry collapses packing to a handful of variables

Superimposing one period of the built chain onto the next (Kabsch) yields a rigid screw motion.
Its axis, rotation angle and rise give the chain axis, the crystallographic c repeat and the
chain radius with no fitting. Chains are then rigid helices, and a two-chain cell has about six
degrees of freedom: a, b, cell angle, two setting angles, a z-offset, and a parallel /
antiparallel flag. Packing a polymorph is a cheap multistart local optimisation over these,
not over 3N coordinates.

### 4. Multi-fidelity funnel

| Stage | Evaluations of the expensive potential |
|---|---|
| Fit RIS energies from a 2D dihedral scan on a short oligomer | ~10^3 single points |
| Enumerate and rank periodic sequences (DP) | 0 |
| Pack top candidates with a fast intermolecular potential | 0 |
| Continuous refinement of torsions and cell within each basin | ~10^2 per candidate |
| Final relaxation / validation of the top few | ~10^1 relaxations |

The scan can be run with a classical force field, a machine-learned potential (MACE etc.) or
DFT; everything downstream inherits that accuracy.

### 5. The amorphous fraction from the same transfer matrix

Exact backward sampling gives the equilibrium single-chain ensemble at any temperature with no
equilibration run. Clamping a segment to all-trans and sampling the rest models chains leaving
a crystalline lamella (loop and tie-chain statistics). The frequency of long trans runs is a
direct measure of beta-nucleation propensity relative to alpha.

## What is and is not new

RIS theory, Viterbi, screw decomposition and crystal structure prediction all exist. The
contribution is the combination: exact discrete steps in place of stochastic ones, chain
symmetry to make packing low-dimensional, and a funnel that calls the costly potential
thousands of times rather than millions.

## Limitations

- RIS uses rigid bond geometry and discrete states; a continuous refinement stage is required.
  Real PVDF states are deflected (beta dihedrals near +/-172 deg, alpha gauche near +/-45 deg).
- Pairwise interactions only in the base model.
- Thermodynamic ranking says nothing about which polymorph forms kinetically.
- Head-to-head / tail-to-tail defects (common in PVDF) need an additional bond type.
- The prototype intermolecular potential uses a damped shifted-force Coulomb, not Ewald.

## Planned package layout

```
src/polyfind/
  polymers.py     monomer templates (PVDF, PE): geometry, charges, RIS states
  ris.py          transfer matrices, Viterbi, k-best, cyclic k-best, partition function, sampling
  chain.py        torsion sequence -> Cartesian coordinates (NeRF), substituent placement
  helix.py        screw decomposition, helix descriptors, symmetry canonicalisation
  enumerate.py    periodic candidate generation, dedupe, ranking, classification
  forcefield.py   simple intramolecular potential; RIS fitting from any calculator
  pack.py         rigid-helix crystal packing search
  amorphous.py    Boltzmann sampling, clamped sampling, conformational statistics
  calculators.py  pluggable energy backends (built-in, optional ASE/MACE adapter)
  pipeline.py     the funnel
  cli.py
tests/            DP vs brute force, geometry, helix parameters, PE and PVDF validation
```

Validation targets: polyethylene (all-trans, orthorhombic 7.42 x 4.95 x 2.55 A) and PVDF
beta (8.58 x 4.91 x 2.56 A), alpha (4.96 x 9.64 x 4.62 A), gamma (4.96 x 9.67 x 9.20 A).
