# PVDF C-F force balance at the rejected stationary cell

2026-09-13 PDT, isolated `physics-vector-repeat`, based on `58cf7ac`.
The preceding goal turn verified live native processes and the completed D:
migration. This continuation produces new force-balance and ablation evidence.
All five physics research stages remain open; GPU remains reserved.

## Finding

The complete model's rejected C-F length is a numerical stationary point whose
force agrees with the derivative of its energy. This audit does not establish
a positive full Hessian or a stable minimum. The radial balance is caused by
the fitted stretch term plus electrostatics, with a large contribution from
the subsequently fitted charge-flux response. Angle bending contributes no
resolved radial force on this terminal F.

The same zero-field, 24-atom/two-chain, four-backbone-site-per-chain protocol
in `CARTESIAN_RELAXATION_PROGRESS.md` again reaches accepted step 173 after
360 complete-owner evaluations. Energy is -28.39082734024751 kcal/mol/cell;
maximum atom force is 6.053288113469552e-5 kcal/(mol A), maximum selected
physical stress 0.01170084031469908 MPa. All six log strains are free.
The existing covalent-distance admission still fails on eight C-F edges.
No rejected geometry is retained or used for a field/pre-strain branch.

At source atoms C3/F5, the actual separation is 1.536630912505807 A.
Positive force below moves F outward along C-F. Forces are kcal/(mol A).

| Complete-model contribution | Radial force |
|---|---:|
| Stretch (analytic owner) | -17.9763583888 |
| Bend (analytic owner) | approximately 0 (less than 8e-16) |
| Ewald including stationary induced dipoles | +119.1200802054 |
| Bonded charge-charge exclusion correction | -99.3408221007 |
| Pair, force-shifted LJ in this Ewald run | -1.8029149699 |
| Torsion / permanent applied field | 0 / 0 |
| Complete analytic force | -0.0000152632 |

Thus electrostatics after its bonded correction contributes +19.7792581048;
LJ reduces the nonvalence remainder to +17.9763431. The correction is already
part of the Hamiltonian. The separately exposed unhalved `pair_correction`
diagnostic is not added again or counted as an independent energy term.

## Fixed-geometry factorial ablations

Only this F's radius is varied; all other nuclei and the cell stay at the
failed complete-model stationary geometry. Each ablation uses the same fitted
blueprint. Removing flux restores base charges and removes their geometry
derivatives; removing induction retains permanent electrostatics. These are
local force comparisons, not refits, relaxed ablation structures, native
references, switching paths or an assertion that one ablation is physical.

| C-F radius (A) | Current full model | Flux removed | Induction removed | Both removed |
|---|---:|---:|---:|---:|
| 1.3500000000 | +45.782298 | +29.651189 | +45.305601 | +29.473860 |
| 1.3605293281 | +43.163388 | +26.894914 | +42.585408 | +26.715477 |
| 1.4673820423, harmonic r0 | +16.842158 | -1.164701 | +15.074020 | -1.362542 |
| 1.5366309125, rejected stationary radius | -0.000015 | -19.478684 | -2.714633 | -19.685540 |

At the last radius, adding flux with induction present changes the force by
+19.4786689704. Adding induction with flux present changes it by +2.7146175503.
These conditional differences include interaction and must not be added as
independent contributions. Even with both additions removed, the original
periodic fitted model pushes outward by +26.715477 at 1.360529 A.

The fitted stretch parameters are unchanged: r0=1.4673820422843398 A and
k=259.5906378151254 kcal/(mol A^2). The historical `VALENCE_FIT.md` already
records this long r0 against its dataset's C-F mean of 1.358 A and attributes
the disagreement to reference geometries generated at a different level
from their PBE-D3 labels. This audit verifies the current model's forces; it
does not independently re-establish that historical explanation or certify
the original reference dataset. The 1.3605293281 A sample is the historical
unconverged finite native step-25 mean documented in the preceding report,
not a matched periodic equilibrium or a calibration target.

The assembled model also differs from the finite fitting Hamiltonian:
`FFParameters.applied` projects finite off-site fluorine charges onto nuclei
at the built chain geometry; periodic Ewald, the Born-fitted flux and the
induced response are then assembled. The Born coefficients were fitted to
response components, not jointly to this fully deformable model's energies
and forces. Their historical source admission remains a separate open gate.

## Implementation and verification

`ChainValence` owns pure stretch and bend evaluation through shared derivative
and scatter implementations. The complete path retains its original
bond-then-angle accumulation order. No parameters, potential expressions,
angle guards, native criteria or topology distance scales were changed.
The example has one `pvdf_model` constructor for the full run and ablations.
Term forces come from finite differences of the complete owner's named
energies, with coordinate-dependent charges re-evaluated at every displacement.
No second Hamiltonian is reconstructed in the diagnostic.

Reproduce with `POLYFIND_DEVICE=cpu`, `OMP_NUM_THREADS=1`,
`OPENBLAS_NUM_THREADS=1`, and this worktree's `src` on `PYTHONPATH`:

```text
python examples/cartesian_relaxation.py --maxiter 500 --bond-diagnosis
```

Actual Windows Python 3.12.10 / NumPy 2.4.6 / SciPy 1.18.1 run:

- Full radial force versus energy differences at h=1e-5 and 5e-6 A:
  maximum error 1.2594227883e-8 kcal/(mol A).
- All three ablations at both steps and four radii:
  maximum energy/force error 9.8381178759e-9 kcal/(mol A).
- Each named energy-term derivative at the two steps:
  maximum difference 1.4210854715e-8 kcal/(mol A).
- Five new component controls cover atom/all-repeat derivatives at two steps,
  component additivity, image-inclusive torque and transverse terminal-F bend.
- Windows: 30 selected valence/Cartesian/supercell-owner tests pass in 27.21 s.
  Linux on the same core source: the same 30 pass in 28.52 s.
- Six Windows exact energy-reference/analytic-stress controls pass unchanged
  in 10.18 s. Example codec control passes in 0.79 s. The final run after
  consolidating the radial displacement helper reproduces every diagnostic
  sample and ablation identically, including the step-173 stationary point.

Final actual output payload SHA256:
`2c13d7c3a9ef6c9b077a96b08d743daaa7baf3580b01465652953af893357511`.
Writer SHA256:
`4037bbc1f9834f47d543266ce4f95e8722fe5b1e29fa847255af7ee6b4ff15f1`.
Actual `pack.py` SHA256:
`24e55596a1803db6abc6e3af251f58e317531418ce3fe7c5a2c8df63f24f5786`.
These are execution file-byte identities; earlier report identities retain
their original writers. The compact diagnosis and reproduction protocol are
retained here; no failed raw geometry, trajectory or scratch is published.

## Native work and next scientific decision

Native handles 21123 and 8143 are confirmed live. Molecular PVDF/SVP has
accepted step 26, E=-75340.81554449623 eV and maximum free force
0.05519605938421371 eV/A, topology intact, still above the unchanged 1e-4
criterion. Its actual checkpoint SHA256 is
`33c06d087375728d524b23c52dd61330f9996624f5347338000e326fd7376838`.
The other three SVP/TZVP geometry cases are pending. Electrical PH still
has iteration 8 residual 7.426e-15 and its dielectric tensor, with final Born
tensors, completion and receipt/source/replay qualification pending.
C: native dependencies and running writers remain unchanged.

The next calibration must constrain the same complete Hamiltonian's energy,
full forces, zero-load geometry and electrical response together on admitted,
matched references. A response-only coefficient fit cannot establish those
other properties. Retain the current failed model as diagnosed; do not tune
r0 to the radius gate or disable flux to label a replacement calibrated.
Qualified native references, full stability/BZ/size/C/S/internal response,
all-target field/pre-strain coordinates, barriers/rates/loss and the accurate
viewer remain required. No research stage is closed by this local diagnosis.
