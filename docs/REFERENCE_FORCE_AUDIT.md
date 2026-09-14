# Historical force labels and the joint-calibration gap

2026-09-13 PDT, isolated `physics-vector-repeat`, based on `19cd208`.
The previous goal turn made progress by identifying the rejected periodic
PVDF cell's force balance and the additional outward charge-flux contribution.
This continuation measures the actual historical training labels, reproduces
the recorded finite fit, and starts an independent native force check.
All five research stages remain open; GPU stays reserved.

## CFE finite inputs and tighter-SCF zero checkpoint

Sarco980d2c71 publishes systems_boundary/cfe_ter_source_motif_v1 with verified
5/7/9-unit inputs,32/44/56 atoms. Source frame38 is the first conf-labelled
CFE-terpolymer frame, selected without energy/force ranking. All22 source
heavy atoms and stereocentres4/10/13 are preserved; actual fixed-core UFF
host/H initialization returns0 for the extensions. Shared preparation66fbaaeb
has34 passing tests, and existing CNEPO receipts replay unchanged. Read the
publication's PROTOCOL.md and VERIFICATION.md before using the coordinates.
They are prepared-only finite inputs, not field response or force references.
VDF/TrFE/CFE counts2/2/1,4/2/1,6/2/1 change composition as well as end distance.

The SCF1e-10 zero point passes full native replay at00:53:08PDT. Relative to
SCF1e-8, energy changes5.45697e-12 eV and maximum atomic force-vector change
is4.61987e-7 eV/A. Net force remains0.7196520756 versus0.7196520837 eV/A,
with identical actual ordinary grids. Tightening SCF alone does not resolve
the defect at this geometry. Displaced points remain required; session78188/
controller12087 continues with minus_001 child13855/start13342504. No refit,
force recentering or calibration admission. See Sarco
trainset_translation_y_scf_tight_probe_v1/MILESTONE_ZERO.md.

## Terminal translation failure and isolated SCF control, 2026-09-14

Sarco27efd41b archives all five native units, inputs/outputs and exact
terminal replay in trainset_translation_probe_v1. Session99954 exited.
Both energy-derived force errors exceed1e-4 eV/A:0.00152559685094 and
0.00043572297446, with inter-step difference0.00108987387648. Energy span
0.00120476604388 eV exceeds1e-6; maximum atomic force change0.020914606 eV/A
and net force0.725961418 eV/A exceed1e-4. All declared checks fail unchanged.
This proves sampled numerical origin sensitivity, not its particular cause.
Terminal result SHA256:
e9e62248f22aef1a466a166281ddd131c39e237a2f5580276b32e2d342e56b91.

The distinct SCF-only precision control is running under Sarco3fdf6c88:
trainset_translation_y_scf_tight_probe_v1/PROTOCOL.md. Session78188,
controller12087/start13272490, zero child12088/start13272892 starts00:41:33PDT;
kernel execution and native SCF iterations1-3 observed. Only EPS_SCF changes
1e-8 to1e-10, leaving the five translations, source, grid and other method
settings unchanged. CPU15/one thread/D:, GPU reserved;24 source files pinned.
The shared compiler/method/recipe owns precision, not per-call-site edits.
Thirty-five source/compiler/probe tests pass (one optional native parser test
skipped), plus15 original grid tests on their required CPU30; both historical
probe results replay with their original hashes and failed gates. Finer-XC
native-unit integration remains pending despite synthetic grid-role tests.
No corrected force labels, refit, field motion or research-stage admission.

## Source identity and first translated point, 2026-09-14 PDT

Sarco 607504e7 publishes
`materials/gpu_bundle/results/boundary_sensitivity/cfe_ter_chemistry_audit_v1/RESULT_REPORT.md`
and its exactly replayable compact result. All 15 cfe_ter frames (7 conf/8
scan) match the anchor-table graph in original atom order: C10H10ClF11,
32 atoms, VDF/TrFE/CFE/TrFE/VDF and S at source atoms 4, 10, 13. The explicit
backbone is 0,1,4,6,9,10,13,15,18,19. Polyfind's CFE homopolymer is a
different chemistry, not an interchangeable baseline. The shared graph
auditor b147245d has 22 passing new/existing controls. Result file SHA256:
91bbaa107504859a9337ef1e8615e3b8096d266708a091e6683c39423ce0759b.
This recovers constitution/stereochemistry, not the missing original
96-point scans, matched field/dipole histories, admitted native protocol,
relaxed packing or qualified energy/force labels. The original 18 A periodic
boundary metadata is retained, not silently reinterpreted as vacuum.

The live common-y experiment now has two whole-unit replay-verified points:
zero E=-8341.2939620562 eV, total Fy=-0.6039086187928724 eV/A;
minus_001 E=-8341.29456387585 eV, total Fy=-0.6028309705562616 eV/A.
The new zero's energy and all 32 force vectors exactly reproduce the preceding
radial zero at parsed precision. The -0.001 A common displacement changes
energy by -0.000601819650 eV: measured numerical origin sensitivity, already
larger than the 1e-6 eV span limit. This does not yet supply a central
derivative, identify the responsible numerical setting, or justify recentering
forces. Partial reduction payload SHA256:
beb11f6d9125108f6ebbe53434af2fa4c15fe7ab4470fb8cd903641e85037042.
Session 99954/controller 9192 remains live, +0.001 native child 10805 observed
using CPU at 00:11 PDT. Complete the unchanged five-point experiment.

## Common-translation control launched, 23:47 PDT

Sarco ab502461 generalizes the existing experiment through one immutable
ForceProbeCoordinate, which owns selected source atoms, displacement direction
and summed work-conjugate force. The radial path and common y-translation
reuse one source/execution/receipt/recovery/reduction implementation. Ten
controls pass in10.788s; the complete archived radial result replays with
unchanged state/result SHA256s and failed gate. Its original writer identity
is preserved, not relabelled as the new implementation.

The new experiment actually runs on D:/CPU15, one thread, GPU reserved:
session99954/controller9192/start12950992, first native child9193/start12951399,
case starts23:47:58PDT and SCF steps1-4 are observed. All32 atoms translate
equally along y at0 and+/-0.001/0.0005A. The zero input is byte-identical to
the preceding zero input; internal coordinates, cell and method stay fixed.
Read sarco trainset_translation_probe_v1/PROTOCOL.md and PHYSICS_RESEARCH.md
(launch relay0d0a7699). No translation point or admission has completed yet.
Its additional net-force/energy-span/force-change checks concern only this
sampled y path, not full origin/grid convergence or physical motion under field.

## Terminal native result, 23:35 PDT

Sarco 84cf6ec5 retains the full five-point native input/output/unit/terminal
result diagnosis and trainset_force_probe_v1/RESULTS.md. Exact replay passes;
session 44730 is exited, not running. The original scientific gate FAILS:
errors at h=0.001/0.0005 A are 0.0001429456973/0.00003271584893 eV/A,
and the difference between derivatives is 0.0001102298484 (limit 0.0001).
Both step values and the inter-step criterion remain part of that result.

Error reduction by 4.3693 on halving h, together with displaced native
force-integral checks, supports a substantial central-difference truncation
contribution. Simpson force integration differs from the energy-derived
interval force by 1.20186487e-5 and -3.41793485e-8 eV/A respectively.
Post-hoc Richardson force 1.1566269825683169 differs from analytic force by
-4.02743386e-6 eV/A; it neither changes the declared gate nor qualifies the
electronic/grid method. The raw native net force 0.7196521 eV/A remains a
separate concern. Next is a separately declared common-translation test
at the same method, followed by distinct numerical-convergence controls.
No labels/model parameters were changed, and no accurate field/pre-strain
geometry or completed material research stage follows from this diagnostic.

## Whole-dataset translation decomposition, 23:18 PDT

The audit now reports the orthogonal decomposition of each raw force array,
model force array and residual into common-translation and internal parts.
For N atoms, m = sum(F)/N and sum(|F|^2) = N|m|^2 + sum(|F-m|^2).
For any zero-net-force model the first term is an irreducible contribution
to unweighted squared force error. This is a diagnostic projection only:
no stored label, fitted parameter, force convention or acceptance gate changes.

All 566 raw frames have nonzero net force, with norms ranging from
0.1486006 to 2.7510042 eV/A, mean 0.9982068 and RMS 1.0935828 eV/A.
The recorded finite model is zero-net-force to numerical precision:
maximum net norm 4.624e-14 kcal/(mol A) over its 467 admitted frames.

| Partition | Frames | Original force RMS | Translation floor | Internal residual RMS | Translation share of squared residual |
|---|---:|---:|---:|---:|---:|
| Fit training | 312 | 4.6445787 | 0.4197092 | 4.6255763 | 0.816591% |
| Fit held out | 155 | 5.2312495 | 0.4413575 | 5.2125978 | 0.711820% |
| PVDF subset | 15 | 3.5979889 | 0.3663947 | 3.5792847 | 1.037000% |

RMS entries are force-component RMS in kcal/(mol A), weighted by actual
component count, not equally by frame. Original RMS squared equals the sum
of the two component RMS squares; these are not linearly subtractable errors.
Thus the net-force issue is widespread but does not explain most of the
fit mismatch. Merely subtracting mean force would neither qualify the native
labels nor resolve the rejected long C-F periodic equilibrium. CNEPO remains
outside the recorded fit, not silently scored with a different model.

Actual audit writer SHA256:
`b58508691791d95350431dd1f72b5b9d4499709d91c98f23ca4e2a6a39e5fab5`.
Source dataset SHA remains unchanged. Seven new synthetic controls verify
orthogonality, the zero-net-force lower bound, atom-count weighting,
non-mutation and malformed input rejection. Combined with existing finite-fit
tests: final 42 pass in 7.55s with the actual dataset path provided, including
the PVDF reference-topology control. The initial run passed 41 and skipped
that one path-dependent test. No source geometries or full output copy are
retained by the audit.

The third native point completed at 23:17:35 PDT, and all three completed
units passed the defining replay. At h=0.001 A, the central energy derivative
gives +1.1567739556994638 eV/A versus analytic +1.156631010002172.
Absolute error 0.00014294569729189632 eV/A exceeds the predeclared 0.0001
limit. This first-step comparison fails; it is not rounded into a pass.
The two h=0.0005 points remain running/pending, needed to distinguish
step-size sensitivity and complete the declared result. No protocol tuning
or label recalibration is performed mid-run. Partial reduction SHA256:
`1cf78168ec1c406c4b12b7741b3a2343544769c4fcd7b6a3ce0dab5925a299a5`.
Session 44730/controller 6211 is live; fourth native child 8071 is observed.

## First native checkpoint, 23:00 PDT

The zero point has now completed47 native SCF steps with clean exit and
passed the defining whole-input/output/receipt/unit replay. It gives outward
C-F1.1566310100eV/A versus source1.5998264914, with energy difference+0.0441423eV.
The source and new native forces therefore share the outward tendency, but
neither is admitted as accurate by this observation. The four displaced
points are still running/pending; their energy/force derivative gate is open.

Total atomic force is0.7196520837eV/A in the new point and1.0548348674 in the
stored label. The raw vectors are preserved. This introduces a separate
common-translation consistency concern; SCF convergence and a single-atom
energy/force check alone cannot establish origin-independent reference forces.
Native grid counts are243^3,144^3,81^3,48^3. Grid-origin error is a hypothesis
to test after the declared five points, not an established cause or a reason
to tune the current run. Full evidence is sarco's
`periodic_reference/trainset_force_probe_v1/MILESTONE_ZERO.md`.
First unit SHA2565c261e786b9717ef88310ce520ceb587824a112048ec56962984fc0cdab9fc80.
Session44730 remains live, now minus0.001A native child6218/start12665829.

## Measured source data

`examples/reference_force_audit.py` reads the existing 566-frame dataset;
it does not fit parameters, relax structures or retain copied geometries.
The defining `load_reference` filter selects467 frames over31 chemistries,
exactly the counts in `VALENCE_FIT.md`. There are19123 atoms,285 conformer
and281 scan frames over37 source systems. The excluded whole systems are
cnene,cnepo,dicnene,ene,epo,fa_tetra. In particular, the16 CNEPO source frames
were never part of this fit. The original fit's structural filter rejects
ring or unsaturated backbones; this is a coverage gap for the all-target goal.

Each terminal F force is projected along its original C->F bond. Positive
means outward. Values below use the reader's actual conversion23.0605
kcal/mol per eV; this factor is recorded explicitly by the audit.

| Source system | Frames | C-F probes | Mean length (A) | Mean label force, outward (kcal/mol/A) | Mean finite-fit force | In recorded fit |
|---|---:|---:|---:|---:|---:|---|
| PVDF | 15 | 150 | 1.3619047 | +29.5457033 | +26.9305824 | yes |
| VDCN | 16 | 160 | 1.3595151 | +28.5273157 | +28.6902463 | yes |
| CNEPO | 16 | 160 | 1.3587951 | +28.7741975 | not evaluated | no |
| FANOME | 17 | 170 | 1.3602662 | +29.3127950 | +28.6523912 | yes |

Every C-F label in these four source-system subsets points outward. For PVDF
the range is+19.5973548 to+39.6909457 kcal/(mol A), at lengths1.3475184 to
1.3745957 A. These are correlated samples from the stored scans/conformers,
not independent statistics or relaxed equilibria. A system label alone does
not establish a homopolymer composition or match a periodic candidate.

The original finite fit reproduces the old reported scores without refitting:
training energy RMS1.0224713 kcal/mol and force RMS4.6445787 kcal/(mol A);
held-out energy RMS1.3593459 and force RMS5.2312495, force correlation0.9343791.
The long fitted C-F equilibrium parameter is therefore consistent with a
large outward force already present in the interpreted training labels.
It is not solely an effect added by the periodic flux/induction model.
This observation does not establish whether the label origin, electronic
method, SCF state or geometry-generation mismatch caused that force.

## Boundary and protocol qualification

All566 extended-XYZ headers declare the same18 A cubic cell, `pbc="T T T"`,
and `Properties=species:S:1:pos:R:3:forces:R:3`. This layout matches the existing
reader's positional columns. `read_frames` discards cell/periodicity while
the recorded finite fit evaluates isolated molecular interactions. The bundle
README describes vacuum-box conformers and explicitly says periodic packing
frames are absent. These header facts do not prove the original CP2K boundary
conditions, converged SCF, basis/pseudopotentials, dispersion settings or
force-sign/unit conversion. Those native inputs/outputs and labeller code
were not located in the inspected local bundle/computational source tree.
No periodic packing qualification is inferred from a `pbc` flag.

The original import commit is sarco0d16b507, which records deliberate labeller
shutdown at566/740 frames. Its statement that the last frame was complete is
not a whole-corpus native energy/force/protocol validation. Reuse for joint
calibration requires source admission or new independently qualified labels,
and explicit boundary compatibility with the evaluating Hamiltonian.

## Exact-frame independent native check

The audit selects the largest absolute PVDF C-F radial residual of the recorded
finite fit: source frame13 (zero-based), C9/F10, length1.3711237546 A. Its label
force is36.8927988049 kcal/(mol A), or1.5998264913983513 eV/A under the reader's
conversion. The finite fit gives24.0311729334; central energy differences at
1e-5 and5e-6 A give24.0311729335 and24.0311729343. This verifies the model
derivative at that stored frame, not the native label derivative.

Sarco physics-native-provenance307c74e1 defines and launches a distinct native
five-point experiment at that exact32-atom C10H12F10 geometry and18 A cell:
PBE-D3(BJ), DZVP-MOLOPT-SR-GTH/GTH-PBE,500/60 Ry, SCF1e-8, full complex
unsymmetrized1x1x1 sampling, zero field. The original point and single-F
radial displacements+/-0.001 and+/-0.0005 A must agree in energy/force derivative
to1e-4 eV/A at both steps. This is an independent method check; the original
labelling protocol is unavailable, so it is not advertised as reproducing it.

Native source/protocol/receipt and lifecycle owners are reused. Publication:
`sarco/materials/gpu_bundle/periodic_reference/trainset_force_probe_v1/`.
The actual run starts22:51:34PDT, session44730, controller PID6211/start12612407,
first CP2K child6212/start12612868. Its first three SCF iterations are observed;
no native result or derivative gate has passed at this checkpoint. One CPU
thread on CPU15, D: storage, at least8 GiB available before each point.
Existing molecular, electrical and VDCN native runs remain separate and live.

## Reproduction and controls

Use the isolated worktree's `src` on `PYTHONPATH`, CPU backend and one BLAS thread:

```text
python examples/reference_force_audit.py D:/sarco-work/physics-native-provenance/materials/gpu_bundle/trainset_v1_566.xyz
```

Source dataset SHA256, verified unchanged before/after each audit:
`8787e90e37d04e8348fb92b04ad556c57c91ac72c5ab273a5213016df078baae`.
Audit writer SHA256:
`640b25b53afb9486594b3e4e3b8c70877c47470c2c79d9c8741b1a5fc152d452`.
Runtime: Windows Python3.12.10, NumPy2.4.6. Actual source/model derivative
error at the selected probe is below8.3e-10 kcal/(mol A). Three existing
finite-fit analytic/finite-difference controls pass in2.80s.

Four new native-preparation controls pass in5.948s: exact source/one-atom
displacement, compiled method/input binding, synthetic sign/two-step derivative
gates and missing-point refusal. They are preparation tests, not native physics
results. Initial test invocation lacked the existing tblite library path; the
standard environment fixes it without changing dependencies. A test using an
assumed higher-precision eV conversion failed by3.3e-6 eV/A; it was corrected
to the actual raw-source projection. The model conversion and native numerical
criteria were not altered. No ASE dependency was added to Windows.

Next: finish and replay all five native points, compare the new full forces
and energy derivatives with these labels, then decide which source data can
support the complete energy/force/response calibration. Qualified molecular
references and electrical tensors are still pending. Stable full mechanics,
all-target field/pre-strain coordinates, transition/rate/loss work and accurate
visualization retain their original scope.
