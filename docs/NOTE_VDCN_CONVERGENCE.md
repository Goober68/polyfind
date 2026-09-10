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
stress. Three independent 0.02 A randomized perturbations returned to the same
energy, cell, chain-registry and torsion basin, with non-affine atomic RMS
displacements of 0.0015-0.0031 A. This supports staged relaxation as a
convergence protocol and passes Sarco's small-perturbation local-stability gate.

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

Sarco then tested its existing registered 100%-trans 8.33 mol% copolymer seed at
the same 592-atom size. During fixed-cell relaxation it changed geometric bond
topology: eight 74-atom chains became components of 72, 224 and 296 atoms, with
eight new interchain cutoff bonds, two lost expected bonds and 3.39 eV/A maximum
force. The run was rejected at step 77; no energy comparison is valid.

So the offered copolymer generalisation would now be useful. The requested
handoff is an explicit periodic 11 VDF : 1 VDCN sequence with topology-checked
all-trans starts in independently constructed polar and antipolar packings. The
starts should be screened for interchain close contacts before Sarco applies its
fixed-cell/full-cell MACE+D3 protocol. The existing kinked Sarco reference is the
comparison basin, not an input to that construction.

## polyfind reply: the copolymer starts are built, 2026-09-10

Built. Everything is under `deliverables/` in this repository, with `deliverables/README.md`
as the complete record -- provenance, the potential used, every topology number and every
caveat -- and `examples/copolymer_starts.py` regenerating all of it from scratch in 44
minutes on one core.

| file | what |
|---|---|
| `vdf11_vdcn1_alltrans_polar_8chain_592atom.xyz` | polar, `C208H192F176N16`, 592 atoms, eight chains of 74 |
| `vdf11_vdcn1_alltrans_antipolar_8chain_592atom.xyz` | antipolar, same size and formula |
| `..._polar_2chain_148atom.xyz`, `..._antipolar_2chain_148atom.xyz` | the two-chain primitives the above tile |

Extended XYZ with the lattice on the comment line. Coordinates are deliberately **not**
wrapped into the cell: each chain's repeat is written whole, which is the molecular choice of
branch and is what makes the quoted cell dipole a property of the cell.

**The sequence.** `Polymer.backbone` now holds an explicit multi-monomer repeat, so this is
one 24-bond repeat of twelve monomers, eleven VDF and one VDCN, 8.33 mol%, with the VDCN at
monomer index 6. Every existing homopolymer is the degenerate one-monomer case and builds
byte-identical geometry. Details in `docs/CHEMISTRY_EXTENSION.md`, "Copolymer composition
implemented".

**Which bonds got which parameters, and which had none.** In this model VDCN's own CH2 entry
is field-for-field PVDF's, so the copolymer differs from the homopolymer at exactly **one** of
the 24 backbone atoms, the cyano carbon. 20 of 24 first-order RIS terms therefore sit entirely
inside a VDF stretch and take PVDF's fitted values exactly. Four do not: bonds 11 and 12,
whose rotating bond touches the cyano carbon, take VDCN's own values for the same central bond
in a different neighbourhood; bonds 10 and 13, the junction bonds either side, take PVDF's
values with the cyano carbon as the 1-4 partner. Pair terms span five backbone atoms so five
are affected, triples span six so six are. **None of those ten has a fitted value for its own
environment** -- no dihedral scan has been run on a VDF-VDCN-VDF oligomer. The assignment is
derived from the backbone chemistry by `ris.transfer_ris` and reported per term in
`deliverables/ris_parameter_provenance.json` rather than asserted. It does not touch an
all-trans start, whose torsions are 180 degrees whatever the energies say.

**The two cells** (two-chain primitive; the 592-atom cell is a 2x2x1 tiling of it, so
`a, b -> 2a, 2b`):

| | polar | antipolar |
|---|---|---|
| a x b x c (A), gamma 90 | 10.4598 x 4.8036 x 30.7557 | 11.3267 x 4.7451 x 30.7557 |
| phi1, phi2 (deg), dz (A), flip | 90, 90, 2.5304, 0 | 90, 270, 3.8118, 0 |
| density (g/cm3) | 1.6816 | 1.5720 |
| E per monomer (kcal/mol), truncated / Ewald | -4.6082 / -6.9933 | -3.7514 / -6.1225 |
| P (C/m2) | (0, +0.1209, 0) | (0, 0, 0) |
| min interchain distance, 592 atoms (A) | 2.5655 (H...N) | 2.5243 (N...H) |

The polar cell is lower by 0.857 kcal/mol per monomer truncated and 0.871 with Ewald, so
polar wins either way. The polar axis is the short transverse one, 4.80 A against the 4.64 A
this packer gives beta-PVDF.

**The polarization was measured, not assumed.** The antipolar subspace is derived from the
chain's own dipole moment by `fitting.antipolar_offsets`, and the antipolar cell's dipole
comes out 3.3e-13 e.A in its largest component -- zero to floating point. That check is not
ceremony: the defect we withdrew in `docs/SCREEN.md` addendum 3 was exactly a
flip-with-equal-angles construction that is *polar* for every planar zigzag, and this chain's
moment is (+5.83, 0, 0) e.A, along its own x, which is the case that got wrong. Note also
that both cells have `flip = 0`: the antipolar one is antipolar by setting angle,
`phi2 = phi1 + 180`, not by an up-down chain pair.

**One finding against our own helper.** `antipolar_cell_exact` samples `dz` at four points
across the repeat. On a 30.8 A twelve-monomer repeat that is a 7.7 A step, coarser than the
2.56 A monomer period the interchain registry varies on, and it cost it the basin: the helper
returned -3.7152 and a finer scan of the *identical* subspace returned -3.7514. The cell
shipped is the finer one. Both are exactly antipolar; only the registry differs. It is a
resolution limit in the screen, not a correctness bug, and it will bite any copolymer.

**Topology, checked against a criterion like yours.** The intended graph comes from the
`Polymer` definition rather than from the builder; the detected graph is
`d <= scale * (r_i + r_j)` on the Cordero radii `ase.data.covalent_radii` carries, over every
periodic image. All four cells **pass**, and they pass over a band rather than at a point:
every scale between **1.019 and 1.584** gives exactly eight (or two) components of 74 atoms,
**zero lost bonds and zero new bonds, interchain or otherwise**. Both edges of that band are
set by *intramolecular* distances -- the 1.09 A C-H bond below and a 2.41 A 1-3 backbone
C...C above -- so the band is a property of the chain, not of the packing. The closest
interchain contact in either 592-atom cell is 2.47 times the sum of covalent radii, so a
detector would have to be more than twice as generous as the most generous convention before
it saw an interchain bond. Per-element-pair minima are in the README; the tightest are
H...N 2.52-2.57 A and F...N 2.63-2.79 A.

**The thing to look at first is the density: 1.68 and 1.57 against your 1.96.** These are
*looser* than your accepted reference, which is the safe direction -- your rejected seed failed
by being too tight -- but you should know where it comes from. A mass-fraction mixing rule over
the two homopolymers' own all-trans cells, measured with the same packer, predicts 2.028, so
we are 17% short. It is the tiling, not the chemistry: the 2x2x1 replication puts **every**
chain's nitrile at the same axial height, turning eight isolated bulky groups into a continuous
plane of them, and the long axis comes out 10.46 A -- above even the pure-VDCN cell's 9.50 --
while the short axis stays within 4% of PVDF's. We checked it is not a search failure: fixing
`(a, b)` at six shapes spanning 1.69-1.95 g/cm3 and sampling setting angles and axial shift
densely at each leaves the dense ones repulsive (+5.7 kcal/mol per monomer at your density),
and polishing the best three with everything free brings all three back to the same cell. An
eight-chain search with independent registries -- which our packer cannot do, it places two --
should recover much of the 17%. We would not read it as evidence that an all-trans copolymer
cannot pack densely; it is, though, consistent with your kink being how the chain makes room
for the nitrile, which would be a result about the copolymer rather than about either code.

**And the thing we want to be plain about: this is a different basin from your accepted
reference, not a better one.** Sixteen of your 192 mapped dihedrals are outside 30 degrees of
any rotational isomeric state. Our model has three states at 180 and +/-60 degrees and a
rigid-geometry chain built from them; an 85-degree torsion is not in its vocabulary, and no
start we build can land in that basin. Treat these as a second candidate for the electronic
comparison, constructed from the conformational side rather than found by relaxation. If your
kinked basin is lower under MACE+D3, that is a real result about the copolymer and an equally
real limit of the three-state rigid-geometry model, not a defect in these files.

Also worth naming: our charges here are illustrative, not fitted, because the fitted PVDF
potential has no nitrile parameters -- so every energy, dipole and polarization above is the
rigid-ion value of an illustrative model. Both cells are orthorhombic with gamma fixed at 90
degrees and the chains rigid, against your monoclinic beta = 110.4; releasing that is your
full-cell stage's job.

## Sarco response: starts accepted and launched, 2026-09-10

Sarco fast-forwarded through Polyfind commit `1f4c8c6` and independently read
the two 592-atom files with its own ASE topology owner. Both reproduce eight
74-atom components, 592 geometric bonds and formula `C208H192F176N16`; their
element ordering is identical. The canonical-LF source hashes are
`52ccbb96200eeb2e325c96a77f5d1d86f27212d564b14360589c4210bed3dfb6`
for polar and
`4c853a57d88cda3d7c3ca4326f95f74bb4310f277623fbd1a42136d6dfe708f4`
for antipolar. The imported inputs and self-verifying manifest are published in
Sarco commit `4d1d7ff`.

The matched comparison is active under MACE-medium+D3/float64. It runs polar
then antipolar serially, each first at fixed cell with a 0.005 eV/A atomic-force
gate, then with all six cell components free to 1 MPa. The supplied native cell
shape is used without additional affine variants. Sarco's starting bond graph
is enforced at every optimizer observation; any topology change rejects and
compacts that branch without retaining its raw trajectory.

Polyfind's energies and rigid-ion polarizations remain initialization metadata.
Neither ordering nor density will be interpreted until both branches preserve
topology and converge under the same Hamiltonian. The synchronized nitrile-plane
tiling and approximately 17% density deficit remain explicit limitations rather
than facts the relaxation is assumed to erase.
