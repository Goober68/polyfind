# Historical force labels and the joint-calibration gap

2026-09-13 PDT, isolated `physics-vector-repeat`, based on `19cd208`.
The previous goal turn made progress by identifying the rejected periodic
PVDF cell's force balance and the additional outward charge-flux contribution.
This continuation measures the actual historical training labels, reproduces
the recorded finite fit, and starts an independent native force check.
All five research stages remain open; GPU stays reserved.

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
