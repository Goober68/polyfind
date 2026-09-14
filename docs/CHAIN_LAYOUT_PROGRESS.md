# Declared independent-chain layout and generalized complete owner

2026-09-13 PDT, isolated `physics-vector-repeat`, based on `52ca399`.
Previous turn made implementation progress. C: native source stays unchanged;
all chemistry/controls CPU-only, GPU reserved. All five research stages remain
open. This is model infrastructure, not native/calibrated material evidence.

## Ownership boundary

`ChainLayout` owns copied immutable per-chain atom IDs in chemical-local order,
one boolean reversal per chain and every source atom's element label. It
requires every atom and every declared chain/local pair exactly once.
`from_cell` resolves the existing Cell's declared chain_of/local_of/reversed_of,
not a proximity detector or an assumption of contiguous storage. Canonical
placement owns its one/two-chain layout directly, including its scalar flip.

`CrystalPacker.chain_energy_and_grad(P,H,layout)` is the same complete model
assembly as the canonical placed energy owner, generalized to arbitrary
numbers of independent chains carrying this model's homogeneous complete
repeat topology. The old canonical entry point delegates to it; there is no
second potential implementation. Layout element/local order is checked against
the model's actual periodic topology before evaluation. Ordered working atoms
are gathered from the declared source IDs, and atom derivatives are scattered
back into the SOURCE order; the cell/lattice derivative holds those source
Cartesian coordinates fixed.

Pair tables and LJ coefficients cover the declared chain count. Valence,
torsion, charge flux and exclusions use each chain's own signed lattice repeat.
Polarization now derives the homogeneous polarizability/Thole arrays for the
declared cell size and uses every chain's explicit reversal, including reversed
chain zero and mixed patterns not expressible by a two-chain scalar flip.
A scalar canonical flip cannot silently infer reversals in a larger cell.
`chain_dipole` and the canonical placed dipole share the same ordered charge
state and generalized stationary polarization solve.

This support is for homogeneous complete periodic repeats; it does not claim
mixed different molecular blueprints, missing chemistry/stereochemistry inputs,
longitudinal replication inferred from partial data, or calibrated physics.

## Verification

Seven new controls cover eight independent chains, per-chain Cartesian
perturbations, arbitrary noncontiguous source permutations, exact source-order
force scatter/dipoles, repeated neutral energy/dipole per motif under2x1/1x2/
2x2 tiling, an eight-chain fitted valence/flux/Ewald/induced-dipole/field model,
sampled atom derivatives, all nine lattice components, all six engineering
strains at two steps, arbitrary mixed reversal patterns and metadata corruption.
They retain strict topology/order/immutability checks and existing tolerances.

Windows extended selection:119 passed70.97s, five expensive screen/predict/
antipolar cases deselected, across chain layout, shared owner, placed energy,
geometry, polarization, phonon, Born and Ewald. Linux focused same-source
selection:61 passed49.02s across layout, owner, placed energy, geometry,
phonon and Born; mounted-D: temporaries use the diagnosed `--capture=sys`.
Earlier new-layout-only controls pass7 in16.40s. Final same-source layout run,
including sampled first-chain/mixed-reversal forces, passes7 in17.83s.
No failed raw payload.
These are model/math/code controls, not material calibration or full stability.

## Remaining migration and research

The legacy SupercellEnergy adapter is NOT yet migrated. Its separate potential
assembly/value-only torsion, energy-term decomposition and cutoff/convergence
diagnostics must route through the declared-layout complete owner together;
changing only its total while leaving terms/diagnostics on another Hamiltonian
would violate ownership. Its current legacy results cannot be promoted as
full energy/force/field/supercell stability evidence. This turn establishes the
contract and complete-owner support required for that migration, not a falsely
completed adapter. Optimized batching is also still open.

General independent-coordinate/full-six-strain minima, stable matched-state
Hessians/BZ/size convergence, full6x6 C/S/internal response, independent matched-
Hamiltonian native/model calibration and all target chemistries remain open.
Actual field/barrier/rate/loss and the accurate field/pre-strain viewer remain
required. Primitive two-site Gamma cannot represent twisting. Earlier scalar
unit snapshots and six-column files retain their actual historical writer
identity; no native/material record is relabelled. Same native molecular and
electrical sessions21123/8143 were re-polled live without restart/source edits.
