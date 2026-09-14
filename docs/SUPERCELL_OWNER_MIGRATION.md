# Complete supercell energy-owner migration

2026-09-13 PDT, isolated `physics-vector-repeat`, based on `2c73693`.
The previous goal turn verified the migrated live D: processes; this turn
changes the model's defining ownership boundary. GPU remains reserved, all
chemistry and verification CPU-only, C: native source unchanged. No research
stage is closed by these implementation controls.

## One complete evaluation

`CrystalPacker.evaluate_chain_cell(P,H,layout)` now owns a single immutable
`CellEnergy`: complete terms, SOURCE-order atom derivatives and all nine
lattice derivatives. Existing `chain_energy_and_grad` pulls its tuple from
that result, preserving writable gradient outputs for existing consumers.
The complete assembly has not been copied into a decomposition helper.

`CellEnergyTerms` reports scaled pair, Cartesian torsion, valence, Ewald,
bare-charge bonded exclusion and permanent-dipole field energy, chain count
and monomers per chain. Ewald includes stationary induced-dipole energy and
its field response when enabled. Its `total` is the defining owner's actual
accumulation, not a caller's independently resummed Hamiltonian. The unhalved
pair bonded correction is a diagnostic already included in scaled pair;
it must NOT be added to total again.

`SupercellEnergy` delegates total, terms, pair/correction diagnostics,
gradients and image-convergence checks to that complete owner using
`ChainLayout.from_cell`. Its copied pair expression, tiled-table cache,
unscaled-minus-correction assembly, constant metadata torsion and stale
valence/flux/field refusals are removed. Reversals and noncontiguous source
atom order remain explicitly declared. Its chain property follows the
packer, which owns the model; a separate chain override is refused.
Its actual dipole delegates to the same complete charge/polarization owner.
The existing rigid-charge `cell_dipole` helper now also resolves the declared
local/source map, but does not pretend to include flux or induced moments.

The geometry owner supplies reciprocal-width image bounds and a lazy grid
for both production cutoff sums and expanded-shell checks. Explicit image
representatives must be finite, unique, representable integers. Effective
images, not raw representative-grid labels, still own bonded scales and
lattice derivatives. Real-space shell equality does not prove reciprocal
space convergence, a relaxed minimum, stable dynamics or physical validity.

## Verification

Nine new model controls cover independently distorted parallel/reversed
chains, triclinic strain, exact noncontiguous source-order force scatter,
dipoles/immutable results, complete observable delegation, actual per-chain
torsion and valence reconstruction, geometry/flux-dependent pair corrections,
image validation and topology ownership. Three sampled atomic components and
all nine lattice components are finite-differenced at two steps in both
harmonic/torsion/field and fitted valence/flux/Ewald/induced/field models.
All three field components at two steps independently satisfy
`dE/dfield = -EV_TO_KCAL*(permanent+induced dipole)` on actual nuclei.
The harmonic control is explicitly illustrative, not fitted material data.

Final-source Linux expanded selection: 78 pass in 108.57s, including all
existing supercell controls and 592-atom copolymer saved-structure readbacks,
layout, geometry, complete-owner and placed-energy controls. Final Windows
new-control run: 9 pass in 10.53s. Existing Windows Ewald/mechanical reference
and analytic-stress controls: 8 pass, 65 deselected, in 12.76s; exact stored
unit values are unchanged. The same expanded Windows selection passes all
78 in 157.01s. Final Linux new-control run, including the added field-dipole
derivative controls, passes all 9 in 15.91s. Neither expanded run deselects
cases; the reference-only Windows command explicitly deselects 65.

Development runs are not labelled green: the first selection had 33 pass/
2 failures (zero-term default harmonic test and differently accumulated
valence reconstruction). Initial broader Windows: 76 pass/2 failures in
155.31s; initial focused Linux: 40 pass/2 failures in 65.52s. The plain test
had also depended on import-order registration of a fitted preset. Controls
now use explicit nonzero harmonic parameters and identical sequential
reconstruction order instead of compensated `sum`. The corrected focused
Linux run passed 42 in 47.92s. No native/material criterion or physical
tolerance was changed, no failed raw payload retained, no historical native
or model record relabelled.

## Native observations and remaining research

Native sessions 21123 (molecular geometry) and 8143 (electrical recovery) were
re-polled live; neither restarted. Molecular PVDF/SVP latest accepted step25
has E=-75340.81518967928 eV and Fmax=0.03747009595109067 eV/A, topology intact,
converged=false at the original 1e-4 criterion. Other three basis/chemistry
references remain pending; the all-four geometry admission is incomplete.

Live `recover64.out` now prints actual electric-field iteration8, residual
7.426e-15 below the 1e-14 target, and dielectric diagonal
(2.252595592,2.235385345,2.447528138). This is an interim native observation:
the controller is running, Born tensors/final native completion and receipt/
source/replay verification remain pending. It is not a qualified material
dielectric or a completed electrical campaign. Pinned native writers and
their journals/lineage are not modified or committed by this model change.

The next model work is general independent-coordinate/full-six-strain
relaxation and matched-state stability, not the old restricted canonical
mechanical parameterization. Stable Hessians/BZ/size convergence, complete
6x6 C/S/internal strain, independent matched-Hamiltonian native calibration,
all target identities/conformations/packing, actual field/pre-strain states,
barriers/free energies/rates/cycling/loss and accurate viewer states remain
required. Primitive two-backbone-site Gamma cannot represent twisting.
Optimized batched screening is still open. No fabricated molecular motion,
many-MPa switching performance or full research completion is claimed.
