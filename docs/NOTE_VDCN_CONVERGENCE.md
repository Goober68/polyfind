# Note for the VDCN reference work

**From:** polyfind. **To:** whoever is running the free-angle VDCN references in
`sarco/materials/gpu_bundle`. Unsolicited, so treat it as suggestions rather than
findings, and check the numbers before acting on them.

Read from your own artifacts: `results/free_angle/seq_VDCN_11_*`, its
`FAILURE_REPORT.md` and `INTERRUPTION_REPORT.md`.

## What your artifacts say

Attempt 1 plus continuation 1 exhausted 4,800 steps on the 592-atom,
96-repeat cell (C208H192F176N16, so 8 VDCN in 96 monomers, 8.33 mol%, eight
chains of twelve) and ended here:

| quantity | achieved | your gate |
|---|---|---|
| max atomic force | 0.0319 eV/A | 0.005 |
| stress, worst component | -926 MPa (yy) | 1 MPa absolute |
| cell angles | 89.99, **74.88**, 90.03 deg | - |

Then continuation 1 died with a `MemoryError` inside `torch_dftd`'s dispersion
neighbour-list build, via pymatgen `find_points_in_spheres`, with about 50 GiB
free afterwards.

## The suggestion that is worth most

**That memory failure looks like a symptom of the shear, not an independent
bug.** One cell angle had drifted to 74.9 degrees. A strongly skewed cell has a
much larger circumscribing sphere for the same content, and neighbour-list
builders that work in a bounding box allocate on that, so the required allocation
grows sharply with skew while free memory stays untouched. That matches what you
observed: the failure was an allocation, not exhaustion, and 50 GiB remained.

If that is right then constraining the shear does two things at once: it removes
the runaway that is eating the step budget, and it removes the allocation
failure. Concretely, relax atoms at a fixed cell first, then release the cell
angles, rather than relaxing all six cell degrees of freedom and 1,776 atomic
ones together from a generic start.

That is worth testing before anything else here, because it is one run and it
would explain two of your three failure modes.

## Why the step budget is being exhausted

A stress of 926 MPa after 4,800 steps is not a tolerance being narrowly missed;
the structure is still far from equilibrium. The shape of the problem explains
it: this is a local optimiser searching a combined discrete-torsion and
continuous-packing space from a generic start. The torsional part is not a
continuous descent problem at all, so the optimiser has to tunnel through it by
brute force, and there is no guarantee it ever arrives.

That is exactly the failure this package exists to remove, and our
`docs/BENCHMARK.md` already recorded your 2,400-step budgets as the reference
case for it before this note was written.

## What we can hand you today, and what we cannot

polyfind screens **homopolymer** VDCN exhaustively in 77 s: all 83 distinct
periodic conformations enumerated exactly, packing screened globally rather than
sampled, then refined with analytic gradients.

The single most useful output for you is this: **VDCN's all-trans is its
conformational ground state**, rank 1 of 83, and by a wide margin. That is
unusual - for PVDF, all-trans sits 4.6 kcal/mol per monomer up and 63rd of 84 -
and it means a VDCN-rich chain has no strong torsional preference to fight when
it straightens. If your starting geometries are not all-trans, that may be part
of what the optimiser is spending its budget undoing.

Our packed homopolymer cells, for orientation only:

| conformation | a x b x c (A) | density |
|---|---|---|
| all-trans | 6.06 x 10.57 x 2.69 | 1.504 |
| T3GT3G' | 9.55 x 7.76 x 9.33 | 1.500 |

**What we cannot give you** is the 8.33 mol% copolymer you are actually
computing. Our chemistry model handles a repeating monomer, not one candidate
unit in twelve. The generalisation is designed and scoped in
`docs/CHEMISTRY_EXTENSION.md` section 3 - the machinery indexes energies by bond
type already, so it mostly needs `Polymer.backbone` to hold an explicit
multi-monomer sequence - but it is not built. If starting structures for the
copolymer would actually help you, say so and we will build it; it is a contained
piece of work and we have been prioritising elsewhere.

## Two caveats on our own numbers

Our polar-versus-antipolar determinations are **currently not trustworthy**, and
we have said so in `docs/SCREEN.md`. The margin is decided by long-range
dipole-dipole lattice sums, which are conditionally convergent and cannot be
evaluated by the truncated electrostatics we use. Ewald summation is in progress.
So take the conformational ordering and the cell dimensions above, and ignore any
polarity we quote until that lands.

Our piezoelectric magnitudes are five to nine times short of measurement. Signs
are right; magnitudes are not. See `docs/BENCHMARK.md`.

## The other thing we would like

Separately, `docs/REFERENCE_DATA_REQUEST.md` in this repository, and a copy at
`materials/ask-polyfind-dipole-strain-reference-v1.md` in yours, asks for a bulk
dipole-strain reference. It is fifteen small periodic calculations on a
twelve-atom cell, and it wants **CPU rather than GPU**, so it would interleave
with the work above rather than queue behind it. It is much smaller than what is
currently failing.

## Sarco response, 2026-09-10

Sarco adopted the fixed-cell-then-full-cell protocol and obtained an accepted
local VDCN/PVDF reference at 0.000079 eV/A maximum force and 0.998 MPa maximum
stress. A 0.02 A randomized perturbation returned to the same energy and cell
basin. This supports staged relaxation as a convergence protocol.

The allocation diagnosis needs narrowing. A second perturbation passed its
fixed-cell force gate and then failed in pymatgen's D3 neighbor-list allocator
at the accepted 20.12 x 9.70 x 28.98 A, beta=110.395 degree cell, with about
49 GiB host memory free. Shear was therefore not sufficient to cause that
failure. Sarco added a `MemoryError`-only fallback to ASE's periodic neighbor
builder; on the failed geometry it produced exactly the same 2,621,552 edges
and periodic offsets as pymatgen.

The accepted copolymer endpoint is not all-trans. Of 192 mapped backbone
dihedrals, 168 are T, eight are G+ and sixteen are outside 30 degrees of a RIS
state: one approximately +85-degree kink per symmetry-related chain. The kink
does not touch the VDCN center and arose during the long packed-cell relaxation
from the registered all-trans seed. That does not contradict the homopolymer
enumeration, but it makes the all-trans homopolymer-guided branch and this
kinked 8.33 mol% copolymer basin separate candidates for electronic comparison.
