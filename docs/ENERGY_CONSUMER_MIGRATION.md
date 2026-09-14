# Complete-energy consumer migration

2026-09-13 PDT, isolated `physics-vector-repeat`, based on `e6150ed`.
C: native science dependencies stay unchanged; CPU-only/GPU reserved.
Previous turn made implementation progress. All five research stages remain
open; this change is not independently calibrated material evidence.

## One Hamiltonian, placement-only pullbacks

Canonical `energy` rows, analytic `energy_and_grad`, refinement's shape
Jacobians, mechanical shape relaxation, phonon forces, atom relaxation and
Gamma Hessians now consume `CrystalPacker.placed_energy_and_grad`.
The separate canonical complete-energy assembly and the separate phonon
complete-energy assembly are deleted, rather than retained behind a flag or
unreachable return. Canonical gradients now only differentiate placement of
atoms/lattice into a,b,gamma,phi1,phi2,dz, shared repeat coordinates and c.

Torsion comes from the actual periodic Cartesian quadruples and contributes
atom/repeat gradients once. The `e_torsion` keyword override is removed;
`chain_batch` supplies geometry only. Rigid-reference metadata/cache values
remain for introspection and rigid screens, not as independent Cartesian
energy facts. Refine/mechanics no longer add a second metadata torsion term
when differencing shape geometry against the complete gradient.
Phonon force diagnostics use placed valence atoms/signed repeat vectors,
not canonical local-frame reconstruction. Runtime notes now state that actual
torsion is included, without promoting Gamma/model frequencies to calibrated
material frequencies.

Canonical rows currently dispatch serially through the complete owner.
Pair arithmetic still uses the selected array backend and image working
memory is bounded. The former optimized batched direct screen must be
generalized through this owner before claiming its old GPU throughput;
retaining a parallel Hamiltonian for speed would not complete the migration.
The rigid table remains a screen, distinct from full deformable relaxation.

## Verification and numerical snapshot diagnosis

Four new controls prove canonical/phonon delegation with an observed owner,
show mechanics does not re-add angle metadata, compare a nonprimitive
torsional directional Hessian contribution at two steps, and compare the
entire 72-by-72 nonprimitive Gamma torsion contribution to an analytic matrix.
At the selected stationary trans torsions the isolated contribution is
`J_phi.T diag((V1+4V2+9V3)/2) J_phi`, assembled for both independent chains
with repeated atom IDs and reversed image translations. Complete-Hamiltonian
minus zero-torsion-Hamiltonian Hessians agree with this matrix at both steps.
Positive curvature of this isolated term is NOT full lattice stability.

Initial Windows owner controls:3 passed13.38s. After the added mechanics
metadata control,4 passed11.73s. Selected canonical gradient/owner controls:
10 passed50.00s. Selected pack/refine analytic controls:30 passed185.80s,
64 deselected. Phonon plus placed-energy controls:29 passed19.08s.
Initial Linux owner/placed/torsion/valence/phonon/geometry/Born selection:
70 passed46.30s. Final Linux same-source selection including the mechanics
metadata control:71 passed43.24s. Exact Windows mechanical unit references
and analytic stress control:5 passed8.98s,29 deselected. Updated exact Ewald
reference:1 passed0.81s,38 deselected. The broad Windows migration run is
exits with284 passed, one optional-native skip and one failure in353.18s:
the already-loaded former exact Ewald snapshot described below. Its updated
exact fixture passes in the separate current-source targeted run; this report
does not falsely label the original broad command green. All other controls
in that broad command pass. The final Linux selection includes the added
fourth owner/mechanics control; it is not the full inherited Linux suite.

One old Windows Ewald unit snapshot expected the former assembly's exact
bits: -8.645931161267756. Untouched C: reproduces that value. The new complete
owner gives -8.645931161267747, a difference of8.88e-15 kcal/mol from changed
image/summation order. Its new explicitly Windows reference fixture retains
an exact `==` assertion and passes; no tolerance/pass rule is loosened.
The corresponding mechanical Windows reference fixtures are likewise updated
as unit snapshots of this actual owner, retaining exact equality. Their
previous beta/alpha/gamma literals remain in git history; PE stays unchanged.
Largest mechanical snapshot change is about2.13e-13 kcal/mol. These scalar
unit snapshots are not independent native references or material calibration.
This does not claim inherited Linux bit-exact or PVDC RIS failures are fixed.
No failed log payload is retained beyond this diagnosis.

## Still required

The standalone supercell energy path has not migrated. It still has its own
assembly and value-only torsion, and its chain-count/reversal/order contract
must join the same topology/potential owner before using supercell/BZ/size
results as complete-energy evidence. Optimized batched screening is also open.

General independent-coordinate/full-six-strain minima, stable matched-state
Hessians/BZ/size convergence, full6x6 C/S/internal strain, independent matched-
Hamiltonian native calibration and all target chemistries remain open.
Actual field/barrier/rate/loss and the accurate all-target field/pre-strain
viewer remain required. Restricted canonical mechanics is not a full6x6
tensor; primitive two-backbone-site Gamma still cannot represent twisting.
Historical six-column data keep their actual `cb17be5` writer identity;
no historical result is recomputed or relabelled by this change.
Native sessions21123/8143 were re-polled live without restart/source edits.
