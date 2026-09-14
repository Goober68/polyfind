# Independent atoms/full-six-strain relaxation: stationarity passes, geometry gate fails

2026-09-13 PDT, isolated `physics-vector-repeat`, based on `a89b8ff`.
Previous goal turn made implementation progress. This turn implements and
actually executes the coupled model relaxation; it does not close any of
the five physics research stages. C: admitted native dependencies remain
unchanged, native sessions21123/8143 re-polled live, GPU reserved/CPU-only.

## Shape and boundary conditions

`CartesianChart` owns independent atomic displacements and six logarithmic
strain coordinates in xx,yy,zz,2yz,2xz,2xy order. Reference nuclei and lattice
rows transform by the symmetric positive stretch `exp(strain_tensor(eta))`.
Macroscopic polar rotation is clamped relative to the reference; individual
chains can move/twist independently. No helical symmetry is imposed.
An O(N) orthonormal Householder gauge removes only mass-centre translation,
leaving3N-3 atomic variables. There is no dense N-by-N null-space matrix.

The exact strain pullback uses the matrix-exponential Frechet derivative.
Scaling optimizer strain coordinates changes conditioning only. Fixed log
components are explicit constraints, not automatically fixed lattice lengths
when free shears couple them. The physical current-affine Cauchy stress and
reference-volume log-conjugate stress are different finite-strain observables;
both full six-vectors are reported. No hidden xz/yz omissions.

`CartesianTolerance` independently requires max atom FORCE VECTOR norm,
selected physical-stress residual and free log-coordinate stress residual.
The latter prevents declaring an optimizer minimum only from a misleading
scaled gradient. Physical selected-component acceptance is an additional
gate for partially clamped log-strain charts, not a claim that log constraints
are equivalent to prescribed finite Cauchy loads. Reactions remain visible.
Default model criteria:1e-4kcal/(mol A) atom force and1MPa stress. Native
DFT criteria, chemistry/topology gates and calibration tolerances unchanged.

`CartesianRelaxer` represents CREATED -> RUNNING -> CONVERGED/UNCONVERGED/
FAILED. It publishes accepted checkpoints, not line-search trial geometries.
Iteration limits and solver stop flags never attest convergence. Terminal
publication failure is FAILED, not success. Recovery explicitly starts a
new owner from accepted parameters with the SAME chart/model; no canonical
chain reconstruction. Checkpoint persistence is supplied by its owner/sink;
this is not yet a source-authenticated durable multi-target production
campaign. `phonon.relax_all_atom` now uses this same fixed-cell relaxation
path rather than a second optimizer plus post-hoc recentering. Its legacy
tuple remains an observation API, not a stability claim.

The topology owner now gathers declared noncontiguous source atom IDs.
A PeriodicChain can explicitly declare a period spanning multiple chemical
repeats via its backbone sites; the same graph owner closes the actual full
period. A primitive Polymer plus a larger atom count still cannot silently
infer that declaration. Invalid backbone/local order is refused. No physical
distance gate or original covalent-radius multiplier changed.

## Actual exploratory polymer calculation

`examples/cartesian_relaxation.py --maxiter 500`, CPU NumPy2.4.6/SciPy1.18.1/
Python3.12.10, actually solves two independently perturbed PVDF chains, each
with four independent backbone sites/two chemical periods (24 atoms).
Initial motif[4.6,8.6,90,0,0,.4,0], independent0.003A perturbations/RNG7,
zero field, all six strain components free, fixed macroscopic polar rotation.
The Hamiltonian carries fitted valence/flux, Ewald accuracy1e-8, stationary
induced dipoles, actual periodic torsion and force-shifted LJ.

Both executions stop at accepted step173 after360 actual owner evaluations:
E=-28.39082734024751kcal/mol/cell, atom-force norm6.053288113469552e-5
kcal/(mol A), maximum selected physical stress0.01170084031469908MPa,
maximum free log stress0.011675476885412355MPa. SciPy status99 is the deliberate
physics-assessor StopIteration, not a failed solve or a convergence assertion
from `success`. Numerical stationarity passes; chemistry admission FAILS.

Worst declared bond: C/F source atoms3/5, image(0,0,0), length1.5366A,
covalent-radius ratio1.1553615883502308. Eight C/F edges are missing at
each existing scale1.05/1.10/1.15; there are no extra edges. At those scales
the components are eight isolated atoms plus two8-atom components;1.20/1.30
recover the two12-atom chains. The nonbond ratio is1.6950075701993217.
The fact a looser multiplier would recover the graph is NOT a pass: the
existing gate remains unchanged. No field/strain branch is seeded from this
unadmitted geometry, and no raw rejected geometry/trajectory/checkpoint is
saved. Only this reproduction protocol and relevant diagnosis are retained.

The actual fitted harmonic C/F parameters are k=259.5906378151254
kcal/(mol A^2), r0=1.4673820422843398A. Harmonic r0 is one Hamiltonian term's
parameter, not its full-model equilibrium bond length. Do not alter r0,
charges, scales or acceptance criteria to make this result pass. Diagnose
the fit's geometry/force balance against independently qualified references.

Diagnostic execution identities:

* original output payload SHA256:
  `7da2ba995e44c03be23db6c8a63d8099d2824d43983833863405f0c4a56f6fb5`;
* detailed-diagnosis output payload SHA256:
  `53a1e81313d0f738d74791f8360a484c01562c5922346070e33557c03acbf20a`;
* actual Cartesian owner source SHA256:
  `9ccba2da1c75f327887169b46e9977b3661f3633517a939a3b49cba4728d305d`;
* detailed topology source SHA256:
  `b5be3d98846098dbfce947abb48363a321d6dfae45b9c2a6f1710f3f73875dcb`;
* detailed writer source SHA256:
  `6baa598b6e4f6ef1525b9b2107d23534d798364ee0d6dbb1993b15297f7adeb7`;
* initial full-geometry identity:
  `bc53f5cc59148a1e0f009ef115049f7d2313d7e3a528e05ac677186131e4d9a8`.

The second writer adds diagnosis; it does not relabel the first writer.
Its self-validating JSON codec owns schema/canonical bytes/whole payload
integrity and refuses mutated fields/unknown envelopes. Failed work emits
no geometry. Model execution is not independently calibrated native science.

## Independent native snapshot, not a calibration acceptance

The molecular PVDF/SVP step25 checkpoint is64 atoms/C20H24F20, nonperiodic,
SHA256`ba4b0db9172b6f7a5b0f80d3f0a39088ce4ebae07577f9b643a25a934d92c1fa`.
The file hash is identical before and after the numeric read and matches the
native writer's checkpoint metadata. Nearest C/F separations range
1.3501652535198765–1.3696268068470914A (mean1.3605293280605644A).
This remains an UNCONVERGED finite molecular reference, Fmax0.03747009595109067
eV/A at the original1e-4 criterion, unlike the periodic model cell above.
It is context for the model diagnostic, not an accepted native minimum,
direct bulk comparison or permission to tune the force field. The all-four
SVP/TZVP reference admission remains incomplete. Electrical PH prints actual
iteration8 residual7.426e-15 and dielectric diagonal2.252595592/2.235385345/
2.447528138; its live handle has not supplied final Born tensors/native
receipt/source/replay completion. Neither native process was restarted.

## Controls and next work

Final Windows and Linux selected runs each pass38 (45.83s/45.44s): Cartesian
math/lifecycle/codec, phonon migration, energy owner and supercell owner.
Ten new controls cover full fitted field/flux/polarization chart pullbacks,
sampled atomic/all six log strains at two steps, all six independent physical
affine stress derivatives at two steps, exact constraints, bounded-memory
translation gauge, analytic atomic/full-shear elastic relaxation, separately
required forces/stresses, iteration limit/recovery/terminal failure and
source-permuted/doubled-period topology. Earlier Windows full topology/
copolymer/supercell selection passes74 in122.02s. The initial new-control
run had8 passes/one failure because a primitive Polymer was supplied for a
declared doubled independent period; the graph API now accepts the writer's
explicit PeriodicChain rather than inferring the graph from atom count.
No failed raw test payload and no loosened tolerances.

Next: diagnose/calibrate the model's geometry against qualified native
references before accepting periodic field/pre-strain seeds; complete
native geometry/electrical admissions. Full matched-geometry stability/BZ/
size and6x6 C/S/internal strain are still required after geometry admission.
All target identities/competing packing, actual field/pre-strain coordinates,
barriers/free energies/rates/cycling/loss and accurate viewer states remain
open. Neither optimizer traces nor a doubled Gamma sample establish dynamic
switching, bulk mechanical performance or all-target research completion.
