# Vector-repeat valence and phonon integration

2026-09-13 PDT continuation on isolated `physics-vector-repeat`, based on
`cb17be5`. C: live science sources remain unchanged; chemistry is CPU-only.
This extends the model's bond/angle geometry, not material calibration or
completion of any of the five physics-research stages.

## Implemented owner

`ChainValence` resolves periodic bond/angle topology once and now uses full
Cartesian repeat translations. `energy_and_repeat_grad` returns the complete
atom gradient (n,3) and repeat gradient (3,), including every image-index
contribution. Scalar `energy_and_grad` pulls that same derivative back to
axial c. Batched energies explicitly distinguish (M,) scalar repeats from
(M,3) vector repeats. `repeat_rows` refuses specifying both representations;
neither explicit fact silently overrides the other. Omitted c remains zero
for a nonperiodic geometry.

Public `bond_lengths` and `angle_values` own topology-aware metrics. Phonon
force diagnostics and built-bond-length preparation now call these instead
of inspecting internal vector helpers. Valence energies reuse the metrics;
analytic gradients differentiate the same bond/angle potential.

For row-vector atoms X and repeat r, each symmetric unit engineering strain
S_J has derivative:

`dE/de_J = sum_atoms gX . (X S_J) + gr . (r S_J)`.

Dropping the second term loses the periodic boundary's contribution, including
xz/yz shear. Repeat gradients are not stress in MPa: volume and the complete
cell Hamiltonian are absent from this single-chain valence calculation.

## Integration defect found and fixed

The preceding vector-polarization extension changed its private derivative
return from four pieces to `(gP, glat, gq)`, with all exclusion-image derivatives
owned by glat. The phonon consumer still unpacked four. Actual phonon test
failed `ValueError: expected 4, got 3`; its consumer is now updated, without
a scalar-gradient compatibility shim. The previously untested consumer is
included in the extended suite.

Changing valence's internal vectors also exposed the built-length helper's
scalar-array assumption (`IndexError`). Its callers now state their metric
intent through the valence owner, rather than constructing private image
arguments themselves. Both observed integration failures are corrected and
their existing tests pass. Failed runs retain this diagnosis only.

## Verification

Six new tests cover beta-like TT and alpha-like TG+TG- periodic PVDF repeats,
all atom components, all repeat components, both finite-difference strain
steps, all six engineering strains, neutrality-independent translation/torque
identities, arbitrary frame changes, reversal, batch semantics, axial pullback
and conflicting repeat specifications. They do not qualify either crystal
geometry or the fitted valence constants.

Final Windows Python3.12/NumPy2.4.6/SciPy1.18.1/pytest9.1.1 run:
223 passed in108.42s across affine valence, phonon, periodic geometry, Born,
polarization, Ewald, forcefield and packing tests. No assertions relaxed.
An earlier narrower corrected run passed88 in39.28s.

Final Linux Python3.12/NumPy2.5.3/SciPy1.18.1, CPU8-15/single-thread BLAS/OMP:
42 passed in13.63s across affine valence, phonon, periodic geometry and Born.
Use `--capture=sys` with D: Windows-mounted temporary storage as diagnosed in
the preceding report. This selected run does not claim the broader inherited
RIS/bit-exact fixture failures are fixed.

The earlier six-column model record remains immutable evidence of its actual
`cb17be5` writer, not a result computed by the new source. All14 stored defining
source hashes independently match `git show cb17be5:<path>`, with canonical
LF bytes. Its native-data/model-validity labels and raw files are unchanged.
It must not be relabelled as a current-source/native/material tensor.

## Next full-energy work

The general placed-cell energy/force owner must live with `CrystalPacker`.
Phonon and canonical refinement should consume it, not accumulate another
parallel potential assembly. Required contributions are:

| Contribution | Present support | Missing full-cell support |
|---|---|---|
| Bonds/angles | General repeat vectors and complete gradients | Assemble each independent chain and its full lattice derivative |
| Charge flux | General repeat vectors and complete gradients | Contract all coordinate/repeat charge derivatives into total energy/forces |
| Ewald/induced dipoles | General lattice and exclusion derivatives | Integrate through the shared complete energy owner |
| Short-range nonbonded | Existing canonical pair potential | Complete triclinic image selection and placed atom/lattice derivatives |
| Backbone torsion | Fourier value from supplied torsion metadata | Periodic topology-derived angles and full Cartesian/repeat derivatives |
| Relaxation/response | Restricted canonical variables | Independent coordinates, full six-strain relaxation, stable Hessians/C/S/internal response |

The existing supplied torsion energy only shifts the Cartesian kernel's value;
it cannot provide the missing torsional forces or stiffness. Adding a constant
or borrowing metadata angles would not repair that ownership. Torsion terms
must derive geometry from their own periodic topology and use consistent
coefficients/conventions across conformational and Cartesian paths.

No full six-strain mechanical tensor, relaxed-ion response, intrinsic many-MPa
pre-strain, field-induced trans/gauche conversion, barrier, rate, phonon accuracy,
loss or dynamics is established. Independent matched-Hamiltonian native
references, response/BZ/size/packing/model gates, all chemistries and the
accurate field/pre-strain viewer remain required. Live D: DFT/PH work continues
without source changes, restart or GPU use.
