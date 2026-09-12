# polyfind

Efficient determination of stable atomic arrangements in semi-crystalline
polymers, with poly(vinylidene fluoride) (PVDF) as the target case.

The idea, in one paragraph: a polymer chain's conformation is a *sequence of
discrete torsional states*, so the search for stable arrangements is done there
with exact dynamic programming instead of stochastic search in Cartesian space.
Crystalline chains are periodic sequences (enumerated as cycles with symmetry
pruning); each becomes a rigid helix whose axis, repeat and radius follow from a
screw decomposition, so crystal packing needs only about six variables; the
amorphous fraction is sampled exactly from the same transfer matrix; and the
expensive potential (force field, machine-learned potential or DFT) is only
called to fit the torsional energies and to re-score the last few candidates.
See [DESIGN.md](DESIGN.md) for the full rationale, complexity analysis, GPU
design and validation results.

## Install

```bash
pip install -e .          # numpy, scipy
pip install -e .[dev]     # + pytest
pip install cupy-cuda12x  # optional: GPU backend for the batched kernels
```

## Quick start

```bash
polyfind pipeline --polymer pvdf            # whole funnel, prints a report
polyfind enumerate --polymer pvdf --top 15  # ranked periodic chain conformations
polyfind pack --polymer pvdf "TG+TG-" --refine --cif alpha.cif
polyfind sample --polymer pvdf --temperature 450 --n-chains 5000
polyfind fit --polymer pe --third-order --out pe_ris.json
```

Python:

```python
from polyfind.polymers import PVDF
from polyfind.forcefield import SimpleFF, fit_ris
from polyfind.enumerate import enumerate_periodic, table
from polyfind.pack import periodic_chain, pack
from polyfind.refine import refine_crystal
from polyfind.amorphous import sample_ensemble, ensemble_stats

model = fit_ris(PVDF, SimpleFF(), third_order=True).model      # ~3k single points
# angles="auto" (the default) keeps PVDF's backbone angles frozen and relaxes them per
# conformer for chemistries whose frozen-angle scan has spurious wells (VDCN, AN, PVDC,
# CFE, CDFE): forcefield.ANGLE_RELAXATION_DEFAULTS, docs/NITRILE_LANDSCAPE.md
cands = enumerate_periodic(PVDF, model, max_period=8)           # exact k-best cycles
print(table(cands, top=10))

chain = periodic_chain(PVDF, model.parse("TG+TG-"), model.states)
cells = pack(chain)                                             # batched search
best = refine_crystal(PVDF, cells[0])                           # torsions + cell
print(best.summary())

states, coords = sample_ensemble(PVDF, model, n_bonds=200, n_chains=5000, T=450.0)
print(ensemble_stats(PVDF, model, states, coords, 450.0).summary())
```

To use a machine-learned potential, wrap any ASE calculator:

```python
from polyfind.forcefield import ASECalculator
from mace.calculators import mace_mp
model = fit_ris(PVDF, ASECalculator(mace_mp()), third_order=True).model
```

## Package layout

| module | role |
|---|---|
| `polymers.py` | monomer templates (PVDF, PE), RIS state sets, UFF nonbonded parameters |
| `ris.py` | transfer-matrix model: Viterbi, k-best, cyclic k-best, partition function, marginals, exact batched sampling; third-order terms via state augmentation |
| `chain.py` | torsion sequence to coordinates (NeRF), batched |
| `helix.py` | screw decomposition, helix descriptors, symmetry canonicalisation |
| `enumerate.py` | periodic candidate generation, dedupe, ranking |
| `forcefield.py` | `SimpleFF` (illustrative), `ASECalculator` adapter, `fit_ris` |
| `pack.py` | rigid-helix crystal packing: batched lattice-energy kernel, search, CIF export |
| `ewald.py` | Ewald lattice electrostatics (`CrystalPacker(coulomb="ewald")`), tinfoil or vacuum boundary |
| `refine.py` | continuous refinement of torsions + cell with a commensurability penalty |
| `amorphous.py` | Boltzmann ensembles, run statistics, lamella-interface sampling |
| `backend.py` | NumPy / CuPy selection for the batched kernels |
| `pipeline.py`, `cli.py` | the funnel and the command line |

## Caveats

`SimpleFF` is a small illustrative potential (UFF Lennard-Jones, point charges,
one Fourier torsion) chosen so the pipeline runs without heavy dependencies. It
reproduces cell dimensions of PE and of alpha-, beta- and gamma-PVDF to within
about 10% and the right chain repeats, but its energy *differences* between
polymorphs are not quantitative. The design puts a real potential in exactly
one place (the RIS fit and final re-scoring), which is where it should go.

Tests: `pytest` (76 tests, including brute-force checks of every
dynamic-programming routine; `tests/test_gpu.py` runs only when CuPy and a CUDA
device are present).  `python examples/benchmark.py` prints per-kernel timings
for the active backend; `bash examples/runpod_gpu_bench.sh` (Linux, e.g. a RunPod
pod) or `.\examples\gpu_bench.ps1` (Windows PowerShell) does the whole GPU
validation: CuPy install, GPU tests, full suite on the CuPy backend, CPU vs GPU
timings.
