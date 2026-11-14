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
pip install cupy-cuda12x  # optional: GPU backend (cupy-cuda13x for CUDA 13)
```

The GPU backend serves two separate things. The batched kernels of the screen move to a
device by themselves when one is present (`POLYFIND_DEVICE=cuda|cpu`, float32). The Ewald
sums and the induced-dipole solve of the large cells are float64 and move only when asked
(`POLYFIND_EWALD_DEVICE=cuda`), and only for cells large enough to gain; a many-chain cell
(`mixedcell.py`) wants `POLYFIND_DEVICE=cpu` beside it and says so, and with its Ewald sums on
the device it evaluates everything else in its batched form, the pair sums on the device too
(the same energy to rounding; `POLYFIND_MIXED_BATCHED=1` asks for that form on the host).

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
| `ewald.py` | Ewald lattice electrostatics (`CrystalPacker(coulomb="ewald")`), tinfoil or vacuum boundary; NumPy, or CuPy for large cells |
| `refine.py` | continuous refinement of torsions + cell with a commensurability penalty |
| `amorphous.py` | Boltzmann ensembles, run statistics, lamella-interface sampling |
| `backend.py` | NumPy / CuPy selection: the batched kernels, and separately the Ewald path |
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
timings.  For the large cells, `bash examples/runpod_ewald_bench.sh` checks the Ewald
sums, the dipole solve and a whole evaluation on the device against the host, times
them at 16 to 128 chains and reads the card's utilisation while evaluations run
(`examples/ewald_device_benchmark.py`; `docs/PERFORMANCE_REVIEW.md` section 17);
`tests/test_ewald_device.py` and `tests/test_mixed_batched.py` run the same device code
path on any machine, on a NumPy-backed stand-in that is as strict about mixing host and
device arrays as CuPy is.  To put a rented card to work rather than time it,
`examples/runpod_jobs_pod.sh` is a pod's start command that follows a branch and runs the
next line of a queue file on it (`examples/gpu_jobs/`: each job a script, the next one
decided by pushing a line, `stop` to end); `docs/PERFORMANCE_REVIEW.md` 17.9 and 17.10 are
what the first two such queues did.  `python examples/mixed_forms_check.py --row <sequence
row>` holds the batched form of a mixed-cell evaluation to the loop form on that row's own
cell, which is the check to run before a new chemistry goes to a card.
