# Complete affine clamped-ion model columns

2026-09-13 PDT, branch `physics-vector-repeat`, based on `a05a00c`.
This is a fitted Polyfind model calculation, not DFT calibration, an observed
material coefficient, relaxed-ion mechanics, dynamics, or loss prediction.
All existing producer material/response gates remain in force.

## Geometry ownership

The canonical packing/relaxation variables still cannot represent xz/yz
cell shear. Fixed-fractional-nuclei dipole evaluation now uses a separate
general Cartesian `PlacedCell`: atoms and lattice rows receive the same
affine map. `FluxTopology` owns bond-image indices and evaluates their full
Cartesian repeat translations. Reversing chain 2 reverses its repeat vector.
Both Ewald exclusion kernels retain vector-repeat derivatives; the induced
polarization solver includes exclusion derivatives in the full lattice
gradient rather than a separately asserted scalar c derivative.

Born and both clamped-ion examples reuse `CrystalPacker.placed_dipole`.
Complete response columns transform all strain/vector indices through
`strain_response`, including general rotations and engineering-shear factors.
Missing columns are refused, not skipped or filled with zero. Proper reduction
uses the same total charge-plus-induced polarization at the baseline.

## Calculated record

`beta_pvdf_clamped_ion_columns.json` schema 2 contains all six columns at
signed 0.25% and 0.5%, baseline geometry/cell/presets, charge-only/induced/total
polarization, source hashes, full frame weights and both proper/improper
definitions. `quantitatively_valid=false`; `relaxed_ion_full_tensor=false`.
No producer Berry values were reused or fitted.

Fine proper column vectors in producer axes (C/m2 per engineering strain):

| Strain | x | y | z | Two-amplitude vector change, C/m2 |
|---|---:|---:|---:|---:|
| xx | -6.41e-14 | -0.2960427853 | -9.88e-15 | 3.18e-7 |
| yy | -3.03e-15 | -0.1558552545 | -5.57e-15 | 3.27e-8 |
| zz | 7.84e-15 | -0.0108437092 | 3.42e-15 | 1.43e-6 |
| yz | 8.37e-15 | -8.01e-15 | -0.0149408730 | 1.51e-7 |
| xz | -3.34e-15 | 3.60e-15 | 1.26e-14 | 2.79e-14 |
| xy | -0.3271735819 | -7.14e-15 | -6.98e-15 | 1.47e-6 |

The largest relative change among resolved nonzero vectors is 0.0132% (zz).
xz is numerically tiny, not an omitted column or a fabricated exact zero;
no relative-amplitude acceptance is claimed at its unresolved denominator.
This does not clear the producer's failed xz amplitude gate or establish
physical accuracy of these model predictions.

Record SHA256:
`19a66432d7e7b7aa0286ea23666a045c46b0f6945bec5fb6aaee71824d25e767`.
XYZ SHA256:
`11869376645c0b093baadfaf871ca9a6cecdd3d2144642d8c9267999c81fe73a`.
All 14 defining source hashes match current canonical-LF bytes.

## Regression evidence

Windows Python 3.12 / NumPy 2.4.6 / SciPy 1.18.1 / pytest 9.1.1, CPU only:
84 geometry/Born/polarization/Ewald controls pass in 33.63s. The wider
forcefield/packing/mechanics run passes 149, fails 3 in 100.85s. Those three
bit-exact mechanics fixtures also fail in the untouched `a05a00c` checkout
with exactly the same computed values; assertions were not relaxed.

Linux Python 3.12 / NumPy 2.5.3 / SciPy 1.18.1, CPU 8-15, single-thread
BLAS/OMP: the wider run before the final additional frame test passes 225,
fails 10 in 160.12s. All ten failures reproduce with identical values in the
untouched checkout (10 failed, 8 passed selected controls, 11.41s). Nine are
bit-exact fixtures; one is the existing PVDC third-order reflection test's
50.0-degree mismatch. Its cause is not established by this baseline check;
it remains a separate RIS diagnosis, not physical validation. The final ten
new geometry/frame controls separately pass on Linux in 8.19s.

Linux pytest's default file capture fails during initialization on the D:
Windows-mounted temporary directory (`truncate` FileNotFoundError); no tests
run. `--capture=sys` removes that test-harness file operation, without
changing chemistry code, native jobs, geometry, or scientific criteria.

`baseline_control.json` retains complete old/new charge and induced dipoles
for 17 observations (zero and every signed shared-domain strain at both
amplitudes), both source maps, exactly one decoded baseline XYZ and all
component differences. Max difference is 1.5543122344752192e-15 e.A.
Control SHA256:
`ffdbc80a8b89ddad90921493a09ea3729c5b862c89cc15bb687757149ba318d7`.
Generated JSON/XYZ are LF byte records; scoped Git attributes preserve that
encoding across platforms. Source hashes use the existing canonical-LF owner.
The two independently relaxed old/new seeds differ slightly due numerical
optimizer paths; the matched-nuclei control, not that seed comparison,
establishes agreement of the shared evaluator domain.

## Reproduce

Force `POLYFIND_DEVICE=cpu`, `OMP_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1`.
Use the same runtime versions for byte-level artifact comparisons.

```text
python examples/clamped_ion_columns.py --out OUTPUT_DIRECTORY
```

```text
python examples/vector_repeat_baseline_control.py --legacy-repository LEGACY_CHECKOUT --record-dir OUTPUT_DIRECTORY --out CONTROL_JSON
```

## Remaining research

General affine cell energy/forces and independent chain coordinates must
extend the relaxation kernel before full six-strain internal response,
stable full C/S tensors or relaxed-ion piezoelectric coefficients exist.
Canonical mechanics is not silently made full by this dipole extension.
Stable matched-Hamiltonian references, native producer Born/epsilon,
cutoff/k-grid/BZ/response convergence, model calibration, all chemistries,
field/pre-strain motion, barriers/rates and low-loss dynamics remain required.
Neither live Sarco source tree nor either native process was changed.
