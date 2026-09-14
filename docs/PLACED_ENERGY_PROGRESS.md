# Complete placed-cell energy/force owner

2026-09-13 PDT, isolated D: `physics-vector-repeat`, based on `b7c2555`.
C: native science sources remain unchanged. Chemistry and all controls are
CPU-only; GPU reserved. This implementation is not native reference data,
material calibration or completion of a physics-research stage.

## Implemented boundary

`CrystalPacker.placed_energy_and_grad(P,H,flip)` owns a complete model energy
at independent Cartesian atoms P and lattice rows H. It returns kcal/mol per
cell, all Cartesian atom derivatives, and all nine lattice derivatives at
fixed stored Cartesian atoms. It does not reconstruct setting angles,
chain-local frames or an axial-only cell.

Included contributions use their existing defining owners:

| Contribution | Defining implementation |
|---|---|
| Short-range LJ/DSF | Packer's pair potential and analytic squared-distance derivative |
| Periodic geometry | PlacedCell nearest-fractional pair representatives and complete cutoff image box |
| Bond/angle | ChainValence full Cartesian/repeat derivative |
| Backbone torsion | ChainTorsion actual periodic geometry, not a metadata constant |
| Charge flux | FluxTopology complete Cartesian/repeat charge derivatives |
| Ewald/exclusions | Existing Ewald terms and periodic exclusion correction |
| Induced dipoles | Packer's stationary polarization solve and complete derivatives |
| Applied field | Permanent charge dipole plus the polarization owner's induced response |

The pair-table construction is extracted into `_pair_charge_tables`:
fixed-repeat, per-row/batched, new full-cell and existing phonon consumers all
use one implementation, including both cutoff shifts and the Ewald/DSF split.
No copied pair-potential or charge-flux implementation is introduced.
Placed charges and energy share one validated charge/derivative state owner.

For short-range image separations D = P_i-P_j+m_eff H, `m_eff` includes the
pair-specific nearest-fractional reindexing. The lattice derivative contracts
these effective integers, not just the small representative grid. Same-chain
bonded scales apply only to the actual axial topological images; reversal
reverses their image indices. Independent chains use signed H[2] for valence,
torsion, charge flux and charge exclusions. Their repeat derivatives contract
back into the physical lattice row with the same sign.

Every symmetric engineering strain has the affine derivative:

`dE/de_J = sum gP . (P S_J) + sum gH . (H S_J)`.

This derivative is not an elastic tensor, stress-free reference or intrinsic
pre-strain. Those require relaxation, stability, numerical/model/native gates
and declared volume/unit conversions.

## Verification and retained diagnosis

Fourteen new controls cover plain LJ/DSF, Ewald and fitted valence/charge-flux/
induced-dipole/field configurations; parallel and reversed chains; canonical
baseline energies; every Cartesian atom and every lattice component; all six
engineering strains at two steps; nonprimitive torsion forces; rotation of
atoms/cell/applied field; neutral translation and zero-field torque identities;
complete triclinic image coverage against a brute-force box; arbitrary stored
atom wrapping for the pair geometry; whole-chain lattice wrapping and its
fixed-coordinate derivative chain rule; invalid geometry and reversal refusal.

Initial failures were test setup, not weakened science gates: the fitted
context rebinds the module's packer class, so a class imported before that
context bypassed its required charge initialization and was correctly refused
by the existing zero-flux invariant. The test now constructs through the
module inside that context. A no-field covariance test also incorrectly
treated None as a vector; it now rotates an explicit zero vector.
No failed output payload is retained beyond this diagnosis.

Windows broad run: 282 passed, one optional-native test skipped in125.59s,
across placed energy, torsion, affine valence, phonon, geometry, Born,
polarization, Ewald, forcefield, packing and fitting. After the final
bounded-memory image-iterator refinement,24 placed-energy/geometry controls
pass again in14.85s. The iterator does not allocate a whole image grid.
Linux final same-source selection:67 passed in26.05s across placed energy,
torsion, affine valence, phonon, geometry and Born. Mounted-D: temporaries
use the previously diagnosed `--capture=sys`; no science criterion changes.
All test assertions and existing native/material acceptance criteria remain
unchanged. These are model/code controls, not calibrated material results.

## Remaining shared-consumer migration and physical requirements

The new owner is implemented and tested, but canonical packing/refinement and
the existing phonon assembly do NOT yet delegate their complete energies to
it. Their supplied torsion energy remains value-only; no accurate torsional
phonon or general all-atom mechanical result is promoted. The next change
must migrate those consumers through the complete owner, preserving the
distinction between rigid screening and exact deformable relaxation, rather
than retaining parallel potential assemblies or adding a compatibility flag.
The current CPU owner is chunked over images, not a claim of optimized batched
GPU coarse-screen throughput.

Then independent-coordinate/full-six-strain minima, stable Hessians/BZ/size
convergence, full C/S/internal-strain response, matched-Hamiltonian independent
native references, all target chemistries, actual field/barrier/rate/loss and
the accurate field/pre-strain viewer remain required. Historical six-column
artifacts retain their `cb17be5` writer identity; none is relabelled as this
solver. Native sessions21123 and8143 were re-polled live without restart or
writer-source edits. All five research stages remain open.
