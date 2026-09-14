# Periodic Cartesian torsion owner

2026-09-13 PDT, isolated D: `physics-vector-repeat`, based on `3092871`.
C: native science dependencies remain unchanged. All controls are CPU-only;
the GPU remains reserved. This is model infrastructure for the full five-stage
research, not independently calibrated material evidence.

## Geometry and derivative ownership

`torsion_geometry` owns one IUPAC angle convention, analytic four-point angle
gradient and three-term Fourier potential. The existing fitting gradient was
extracted here; fitting, finite-molecule energies, rigid packing values and
periodic Cartesian terms consume that owner rather than duplicate the potential.
The finite-molecule backend retains its coefficient dtype and arithmetic order.
Undefined zero-bond/collinear-arm geometries are refused, not clipped.

`chain_torsion` resolves periodic atom IDs and image indices from the existing
five-block chemical template. Each torsion representative starts in block zero;
its index matches the builder's dihedral index and the bond-type coefficient
phase. Angle metadata is not an energy input. `ChainTorsion` owns immutable
copied topology and coefficients, derives angles from complete image quadruples,
and returns energy, every Cartesian atom derivative and the full vector-repeat
derivative. Repeated atom IDs scatter-add all their image contributions.
The scalar axial derivative pulls back the same vector derivative.

For atom coordinates X, repeat r and engineering-strain generator S_J:

`dE/de_J = sum gX . (X S_J) + gr . (r S_J)`.

Changing an image representative leaves energy/forces/repeat derivative intact.
Wrapping independent stored atoms changes the fixed-coordinate repeat derivative;
the coordinate-chain-rule term restores the same physical affine derivative.

## Sampling diagnosis

The first stiffness test used two independent backbone sites, an all-trans
primitive repeat. Any two-site translation-periodic zigzag lies in the plane
spanned by its site separation and repeat vector. Arbitrary small Cartesian
perturbations and affine strains preserve its trans torsions. Correct image
forces cancel; zero primitive-Gamma torsion stiffness is not a gradient defect.

The nonzero-force/stiffness controls therefore use a declared four-backbone-site
repeat with independent site coordinates. The two-site case remains an explicit
zero-torsion regression. At the stationary four-site trans geometry, the tested
directional stiffness is about 24.223392 kcal/(mol A^2), agreeing at two steps
with the Fourier curvature contracted through the analytic angle gradient.
This is a selected model direction, not a full stable crystal Hessian or a
frequency. Primitive-Gamma support cannot replace chain-length/BZ convergence.
The inherited phonon documentation's blanket lower-bound claim is corrected
at its documentation/runtime-note owner, without changing calculated values.
Omitted Cartesian curvature need not be positive semidefinite away from a
torsion minimum; a nonstationary reference is not a qualified vibrational state.

## Verification

Eleven periodic-torsion tests cover builder phases for PVDF TT/TG+TG- and PE,
geometry rather than metadata, complete coordinate/repeat derivatives, two
finite-difference strain steps/all six engineering strains, force/torque
identities, frame/reversal/batch/axial semantics, stationary trans curvature,
primitive-repeat sampling, image/atom-wrapping gauges, shared fitting gradients,
immutability and invalid-input refusal. No physical criterion is weakened.

Windows fitting plus torsion: 45 passed, one optional-native test skipped,
4.15s. Linux torsion/affine-valence/phonon/geometry/Born/fitting selection:
82 passed, one optional-native test skipped, five reference/valence-design
cases deselected, 16.37s. The latter are not evidence that inherited Linux
bit-exact fixtures or the separate PVDC RIS reflection mismatch are fixed.
Linux uses the previously diagnosed `--capture=sys` on mounted D: temporaries.

Broad Windows final run: 268 passed, one optional-native test skipped in
102.80s, across torsion, affine valence, phonon, geometry, Born, polarization,
Ewald, forcefield, packing and fitting. After the documentation/runtime-note
correction, all 15 phonon controls pass again in 3.18s. These are model/code
controls, not independent native references or material calibration.

## Remaining implementation and physical scope

This term is not yet assembled into a shared full placed-cell Hamiltonian.
The existing canonical packer and phonon Cartesian assembly still treat the
supplied torsion energy as a value-only constant. Their passing controls do not
establish torsional phonon accuracy, stable all-atom minima or dynamics.

The next owner must unify independent placed coordinates, full lattice gradients,
short-range triclinic images, valence, charge-flux derivatives, Ewald/induced
dipoles and geometry-derived torsions. Canonical refinement and phonon consumers
must route through that owner instead of maintaining parallel assemblies.
Then full internal/strain relaxation, stability/BZ/size gates, six-by-six C/S,
matched-Hamiltonian native calibration, all target chemistries, actual barriers,
rates/loss and the accurate field/pre-strain viewer remain required.
Historical six-column files retain their actual `cb17be5` source identity;
they are not recomputed or relabelled by this extension.
