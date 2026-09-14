# Resume: polyfind, state as of 2026-09-14

02:34 PDT: QS-only five-point control is COMPLETE and archived in Sarco
3af06c96; all six numerical limits still fail. Tighter SCF and tighter QS
alone both fail to remedy the defect; no recentering/refit. CFE is now
actually RUNNING from readye69274fd: session49815/process23124/start13938253,
CPU15/one thread. Its5mer starts02:33:14; native step2 projected force maximum
0.000469764 eV/A, not converged. 7/9 pending, GPU reserved. Read Sarco QS
RESULTS.md and CFE LAUNCH.md. The reference-data relay now lives at
D:/sarco-work/polyfind-reference-relay; merge27fb52d preserves our relay and
the concurrent upstream phonon-screen updateaad1362 without changing the
original C: source checkout. That new screen needs scientific review before
its model claims inform physical qualification.

Sarco e69274fd publishes the actual READY CFE source-motif5/7/9 reference
campaign, definition472ea450. All three are pending, no native calculation,
units or result. One shared input/reference path serves CNEPO and CFE;
explicit archive locations bind the original factory from D: without
rewriting scientific identities. The CFE core is released after preparation;
terminal clamps and collective azimuth remain. Same GFN2 accuracy0.001,
500-step cap and both fresh raw/projected1e-6eV/A force gates.
Read Sarco cfe_reference_candidates_v1/PREPARATION.md. READY SHA2adfd8d6...
Input admission96d54a10... differs from original factoryc65daac9....
22 shared/chemistry/basin tests and12 declaration/CLI tests pass; archive/
curvature regression25 pass/two opt-in skips. Exact original CNEPO input and
whole three-candidate reduction are preserved. Its historical native writer/
runtime relocation remains separate. CFE execution awaits the live final QS
point's CPU slot; GPU reserved. No physical model or field/strain admission.

Sarco c3c802dd publishes the exact neutral CNEPO5/7/9 dipole-axis analysis
(code c62a6693). Transverse magnitudes are 0.935212/0.698476/0.467867 e A;
signed axial components -0.433061/-0.531185/-0.680422 e A. Axial swiveling
cannot change the axial component and is prohibited by the archived held-
azimuth boundary. These are not field-induced changes, length convergence,
bulk polarization or a switching result. Actual archive replay passes with
new native calls forbidden; shared campaign and geometric controls pass.
See Sarco results/boundary_sensitivity/cnepo_dipole_axis_v1/REPORT.md.

QS-only zero/minus_001 both now replay; their forces change from original by
at most 5.94329e-7 eV/A, with net-force norms still 0.719652/0.725961 eV/A.
Three points remain; plus_001 child21418/controller19778 confirmed live.
See QS MILESTONE_TWO.md. No recentering/refit/physical admission; GPU reserved.

01:45 PDT: Sarcob13e3276 archives the completed five-point SCF-only control.
Full replay passes but all six numerical limits still fail; 100x tighter SCF
does not remedy this frame's defect. The unchanged QS-only ready control is
now genuinely running: session95548/controller19778/zero child19780, starts
01:44:28 PDT, CPU15/one thread, EPS_DEFAULT1e-16 with ORIGINAL EPS_SCF1e-8.
The explicit relocated archive also binds the real nine-source factory;
18 tests pass with zero force calls. Shared CFE orchestration and historical
native-receipt relocation remain pending. No refit or physical admission.

01:34 PDT: Sarcod3faa6d4 repairs the archive reader's original-identity/current-
storage boundary. All nine original source identities/coordinate hashes replay
from D: and match a separate unchanged-location read; seven tests pass41.063s.
Exact14-file/154820-byte Polyfind receipts are preserved because mixed line
endings cannot be reconstructed from Git-normalized history. Factory and
terminal-native relocation integration remain pending; CFE is not launched.
SCF-only control now has four independently replayed points, unchanged forces
within4.64e-7 eV/A. Its h=.001 derivative error remains.00152558426 eV/A.
Last plus_0005 child15174/controller12087 runs; QS stays ready, GPU reserved.

01:12 PDT: Sarco15d4a3f4 declares the independent QS-only translation control
(EPS_DEFAULT1e-16, ORIGINAL EPS_SCF1e-8). It is prepared/verified READY with
no native execution, pending completion of the live SCF control on CPU15.
SCF minus_001 now also reproduces original: max atomic-force change4.63481e-7
eV/A; two of five complete, plus_001 running. Native sources unchanged.
CFE native relaxation is held before execution by a reproduced archive-reader
relocation defect: historical C: identity paths vs D: storage are conflated.
Original receipts must remain intact; repair the archive location boundary,
not historical hashes. See Sarco CFE REFERENCE_LAUNCH_DIAGNOSIS.md. No refit.

Latest: Sarco980d2c71 publishes verified CFE-terpolymer source-motif5/7/9
inputs (32/44/56 atoms), using shared preparation66fbaaeb. All22 source heavy
atoms and three stereocentres are preserved; CNEPO archives replay unchanged.
These are prepared-only, changing VDF dilution with length, not equilibria.
The tighter-SCF zero point also completes: max force change4.61987e-7 eV/A
but net force remains0.7196521 eV/A. Read REFERENCE_FORCE_AUDIT.md; four
translated points remain, controller12087/current child13855, GPU reserved.

Latest00:41 PDT: original translation diagnostic is COMPLETE, archived in
Sarco27efd41b; every declared derivative/invariance gate fails. Read
REFERENCE_FORCE_AUDIT.md and Sarco trainset_translation_probe_v1/RESULTS.md.
The distinct tighter-SCF five-point control is now RUNNING, definition
3fdf6c88, session78188/controller12087/start13272490, native12088/start13272892.
Only EPS_SCF changes1e-8 to1e-10; same source, translations, method/grid and
gates. CPU15/one thread, D:, GPU reserved. Its24 defining files are pinned.
Both historical probe results replay unchanged after shared compiler/precision
integration. No new calibration result yet; all five stages remain open.

Current checkpoint (supersedes the dated launch/history below): Sarco
607504e7 recovers the actual cfe_ter source as VDF/TrFE/CFE/TrFE/VDF,
C10H10ClF11, with S at source atoms 4/10/13 in all 15 stored frames.
This is not Polyfind's CFE homopolymer. See REFERENCE_FORCE_AUDIT.md for
the source-linked report; identity is admitted, force/field labels are not.

The live translation diagnostic has two replay-verified points (zero and
-0.001 A); +0.001 A is running, child 10805 observed using CPU at 00:11 PDT.
The new zero reproduces the prior energy and all forces exactly. The first
rigid shift changes energy by -0.000601819650 eV despite unchanged internal
geometry/cell. This exceeds the declared 1e-6 eV span limit already; the
five-point derivative/force diagnosis remains incomplete. Do not stop or
modify the pinned experiment. Molecular PVDF/SVP reached step 28 with
Fmax 0.0454956382 eV/A, topology intact, still unconverged. Other jobs remain
live, GPU reserved, all five research stages open.

## Historical checkpoints (newest first)

The distinct common-y-translation experiment IS RUNNING: sarco definition
ab502461, launch relay 0d0a7699, publication trainset_translation_probe_v1.
Session 99954/controller 9192/start12950992, native zero child9193/start12951399,
case start23:47:58 PDT; SCF steps1-4 observed. All32 atoms move equally by
0 and +/-0.001/0.0005 A, same cell/internal geometry/electronic method, CPU15
one thread, D: and GPU reserved. New zero input matches the prior zero bytes.
Shared ForceProbeCoordinate owns displacement and total-force projection;
one existing execution/receipt/recovery path. Ten controls pass10.788s, and
the completed radial result replays unchanged, including its failed gate.
Translation source/protocol23-file pin is now live. No translation result
or material admission yet. Other native pins and all five stages stay open.

The radial native probe is now COMPLETE, not running: session 44730 exited,
five native units and terminal result pass exact replay. Sarco 84cf6ec5
archives the compact evidence and trainset_force_probe_v1/RESULTS.md.
The original force gate FAILS: h=0.001 error 1.42945697e-4 eV/A,
h=0.0005 error 3.27158489e-5, step difference 1.10229848e-4 (limit 1e-4).
Error ratio 4.369 and force-integral checks support a substantial finite-step
curvature effect. Post-hoc Richardson error -4.0274e-6 is not a replacement
gate. Raw net force 0.7196521 eV/A remains unqualified. Next: separately
declare/run common rigid translation at unchanged method, then distinguish
SCF/integration/XC-grid controls. No follow-up is launched at this checkpoint.
The radial experiment's source pin is released by terminal completion; other
live molecular/electrical/VDCN pins stay intact. All five stages remain open.

Latest diagnostic: all 566 historical frames have nonzero net force
(0.1486-2.7510 eV/A), but translation contributes only 0.71182% of held-out
squared force-fit error. The finite model's total force is numerically zero.
Audit v2 decomposes raw/model/residual forces without changing any labels;
read docs/REFERENCE_FORCE_AUDIT.md. Seven new controls plus finite-fit tests
pass 42 in 7.55s with the actual dataset, including reference topology.
Three native force points now pass receipt replay. The h=0.001 A derivative
misses the declared 1e-4 eV/A gate: error 1.42945697e-4. Do not round or
retune this into a pass. Two smaller-step points remain running/pending;
session 44730/controller 6211, fourth child 8071 live at 23:18 PDT.
All pinned native source/protocol owners remain unchanged. Molecular PVDF
accepted step 27, Fmax 0.0590822 eV/A, topology intact, still unconverged.
Full five-stage/all-target objective and GPU reservation remain unchanged.

Native force-probe checkpoint: zero point completes47 SCF steps, exact
source/receipt/unit replay passes. New radial C-F1.1566310eV/A versus
label1.5998265; energy/force derivative still awaits four displaced points.
Native net force0.7196521eV/A and label1.0548349 introduce an additional
common-translation/origin-consistency concern. Do not correct the raw force
sum or fit against it as calibrated. Read REFERENCE_FORCE_AUDIT.md's first
checkpoint and sarco trainset_force_probe_v1/MILESTONE_ZERO.md. Session44730
is live, minus0.001A child6218/start12665829; current pinned307c74e1 source
and protocol stay unchanged. Translation/grid investigation follows the
declared five points. Full goal and GPU reservation remain intact.

Latest reference audit: `docs/REFERENCE_FORCE_AUDIT.md` directly measures the
566-frame labels and reproduces the recorded467-frame finite fit. Every
sampled PVDF C-F label points outward (mean29.5457kcal/mol/A at1.361905A);
the base fit follows it. Six systems, including all16 CNEPO frames, are
excluded by the fit's structural filter. All headers declare18A periodic
vacuum boxes; original native labelling protocol remains unadmitted.
Distinct five-point native PBE-D3BJ force/energy check is actually running
on D:/CPU15, sarco307c74e1, session44730/controller6211/start12612407.
Read its exact protocol and wait on that handle; do not restart on quiet
output. No native comparison is complete at this checkpoint. All five stages
remain open, original native writers unchanged and GPU reserved.

Latest force diagnosis: `docs/CF_FORCE_BALANCE.md` records actual C-F radial
forces, complete-owner term derivatives and same-geometry flux/induction
ablations. The rejected 1.5366309125 A stationary radius balances inward
stretch17.97636 with outward electrostatics after exclusions19.77926 and
inward LJ1.80291 kcal/(mol A). Flux adds19.47867 conditional on induction;
even both response additions removed still push outward26.71548 at1.360529 A.
This directs joint energy/force/response calibration on admitted references,
not radius-gate retuning. Thirty focused controls pass on Windows and Linux;
six exact Windows reference/stress controls are unchanged. Native PVDF/SVP
step26 remains unconverged, PH Born completion pending, both handles live.
All five stages stay open; C: native source unchanged and GPU reserved.

Latest coupled relaxation: `docs/CARTESIAN_RELAXATION_PROGRESS.md` records
independent atoms/all six strains, exact positive-stretch pullbacks, explicit
translation gauge and a guarded convergence/checkpoint lifecycle. Actual
24-atom fitted PVDF model: numerical stationarity PASS but covalent-distance
admission FAIL (C/F1.5366A); no field branches or saved rejected geometry.
Diagnose calibration against qualified native references; do not loosen
the gate. Ten new controls/final38-pass selections are recorded there.
Native references remain unconverged/incomplete and all stages open.

Latest supercell migration: `docs/SUPERCELL_OWNER_MIGRATION.md` records the
complete owner's immutable terms/derivatives and shared image diagnostics.
Supercell total/terms/gradients/dipoles now delegate through explicit layout;
no copied pair expression, constant torsion or ignored field/flux remains.
Full general relaxation/stability/C/S and native/calibrated/all-target
field/barrier/loss/viewer requirements remain open. Older notes below are
historical checkpoints, not the current adapter state.

Latest layout continuation: `docs/CHAIN_LAYOUT_PROGRESS.md` records immutable
declared atom/local/reversal/element layout and the complete owner's arbitrary
homogeneous independent-chain count, force scatter and dipole/polarization
support. Legacy SupercellEnergy total/terms/diagnostics migration is still
pending; do not promote its old assembly as full supercell/field evidence.

Latest migration: `docs/ENERGY_CONSUMER_MIGRATION.md` records canonical
energy/refinement/mechanics/phonon delegation to the complete owner, removal
of duplicate assemblies and double torsion, entire nonprimitive Gamma torsion
matrix controls and exact Windows unit-snapshot updates. Standalone supercell
migration, optimized batching and full relaxed-cell/native/material research
requirements remain open. Historical extension notes below are checkpoints.

Latest full-energy continuation: `docs/PLACED_ENERGY_PROGRESS.md` records
the complete independent placed-cell model owner, all atom/nine lattice
derivatives, actual periodic torsion, general triclinic image geometry and
one pair-table construction. Canonical refinement/phonon consumer migration
is still pending; no new material result or calibrated dynamics is claimed.

Latest isolated continuation: `docs/CARTESIAN_TORSION_PROGRESS.md` records
the shared Fourier/dihedral owner, periodic Cartesian/repeat derivatives,
four-site trans stiffness and the two-site primitive-Gamma sampling limit.
This term is not yet assembled into a shared full-cell Hamiltonian; no
field-induced conformation, calibrated dynamics or mechanical tensor is claimed.

Further evening continuation: `docs/AFFINE_VALENCE_PROGRESS.md` records
full vector-repeat bond/angle gradients and metric ownership, corrected
phonon consumers,223 Windows/42 focused Linux passing controls and the
remaining shared full-cell energy/torsion/relaxation work. Historical
six-column artifacts still belong to `cb17be5`, not the newer source.

2026-09-13 evening extension: this is the isolated `physics-vector-repeat`
worktree on D: (`D:\sarco-work\polyfind-vector-repeat`), based on `a05a00c`.
Read `deliverables/clamped_ion_vector_repeat_v2/README.md` for general
Cartesian image/dipole ownership, all six affine model columns, matched
17-point legacy control and inherited platform/RIS diagnoses. Full relaxed
six-strain energy/forces/internal strain/C/S are not implemented or claimed.
The C: checkout stays on the collaboration branch with unchanged live source.
GPU is reserved; chemistry is CPU-only. Statements below describe the earlier
collaboration checkpoint, not this extension's runtime/platform evidence.

Read this first after a context clear. It is the map; the documents it points to
are the territory. Everything below is committed and pushed on branch
`claude/polymeric-stable-arrangements-uh06b1`; the working tree is clean, no
agents are running, no worktrees exist.

## What this is

`polyfind` searches for stable chain conformations and crystal packings of
semi-crystalline polymers, PVDF and its candidate relatives, and computes their
field and strain response. The novelty is the search: exact dynamic programming
over discrete torsional states instead of stochastic simulation, screw
decomposition to collapse packing to six variables, a tabulated chain-pair
interaction with FFT lattice sums, analytic gradients throughout. `DESIGN.md` is
the architecture; `README.md` is the entry point; `docs/BENCHMARK.md` is the
measured verdict against the goal.

The goal (`/goal`): a polymer solver that fits the field/strain measurement
needs and is 10-100x faster than existing methods. The `/goal` hook is active
and blocks stopping until met; see "Goal status" below for the honest position.

## The collaboration, and the rule that governs it

**The polyfind repository is the communication channel** with "Dipole", the
DFT/MLIP producer working in the sibling repo `E:\source\gh\sarco`
(`materials/gpu_bundle`). Both sides append dated sections to
`docs/REFERENCE_DATA_REQUEST.md`; that file is a live two-way exchange and its
last ~30 sections are the record of the past four days. Read it from the bottom
up before replying to anything. Never write into `sarco` except as an untracked
file; it has staged work in flight. Rebase-conflicts on that document are
routine (both sides append): resolve by keeping both sides, no markers.

Dipole now also edits Polyfind source directly (2270444 added
`chemical_graph.py` and made `Structure.bonds` a derived view). Review their
code diffs, run the suite, and record the verdict in the exchange.

Dipole works on a different machine with ROCm; this machine has no usable GPU
(`docs/PERFORMANCE_REVIEW.md` section 14). Large runs go to RunPod if ever
needed. The user (Niall) also sometimes drives Dipole's subagents directly.

## The governing requirement (Niall, 2026-09-13)

**The sarcomotor's biggest ask is low tan-delta and a high-frequency
response.** Any strain achieved under those conditions makes a workable
polymer; strain magnitude (d33 = -32) is secondary. Consequences: the
piezoelectric-shortfall chase is validation, not the gate; the
polar-versus-antipolar lattice margin (hard ferroelectric vs relaxor) is the
column that matters and it is the one the model cannot resolve; crystal
phonons (soft modes, lowest optical frequencies) are the intrinsic proxy the
solver can add. The tool computes no loss and no dynamics today. Recorded
in the exchange as the consumer note of 2026-09-13.

## Goal status, in one paragraph

**Speed: met**, well past target (`docs/BENCHMARK.md`): a full response tensor
in 0.3-7.8 s against the baseline's 87-143 s for a single finite-difference
derivative, and the baseline does no structure search at all. **Static side:
validated** against Dipole's periodic PBE-D3 to within a few percent on
polarization, transverse dielectric response and transverse Born charges (the
chain-axis components of both remain open in Dipole's convergence ladder).
**Piezoelectric coefficients: short**, d33 = -8.8 vs measured -32, d31 = +0.02
vs +20, and the reason has been narrowed by elimination to one of two things
neither side has: the electronic clamped-ion term, which no classical charge
model represents, or the 0.4-0.6 C/m^2 Berry-slope target itself, which came
from a sweep that failed its own gate. The narrowed request to Dipole is a
**clamped-ion piezoelectric tensor**, which discriminates. Until it arrives, the
tool's honest scope ends at the static side, and that boundary is sharper than
the goal asked for.

## What was eliminated, in order (do not re-run these)

1. Speed was not the constraint; the algorithm was rebuilt around exact DP and
   the pipeline fell from 289 s to 25 s warm. Nine defects found on the way,
   most by reformulating rather than testing (`docs/PERFORMANCE_REVIEW.md`).
2. The potential was unfitted, not merely cheap. A crystal-data fit failed
   (`DESIGN.md` 5.7); a DFT-data fit generalised (5.9); valence terms made the
   forces usable (5.10). The Berry-slope shortfall was never a fitting problem.
3. **Charges**: a Born-consistent charge-flux fit succeeded (4 params, 10
   observations, rms 0.08 e) and moved d33/d31 *away* from measurement, because
   on rigid pendants only the group Born sum acts under strain and it is 0.08 e.
   `docs/ELECTROMECHANICS.md` 5.9.
4. **Polarizability**: literature values, nothing fitted, reproduce the
   transverse dielectric tensor to 3% and close 23% of the d33 gap; cannot reach
   -32 at any physical value. `docs/ELECTROMECHANICS.md` 5.8.
5. **Internal strain / rigid pendants**: Dipole supplied a geometry-only
   Jacobian; the two models differ exactly where predicted (their C-F bonds
   stretch 0.24 A and turn 31 deg per unit strain, ours rigid to 1e-13) and the
   difference carries no polar dipole through either model's charges: +0.33 vs
   +0.34 C/m^2. `docs/INTERNAL_STRAIN.md`. Phrased as a difference between
   models, not real flexibility, at Dipole's request.
6. Electrostatic scale as a single fix: dead, optima 250x apart
   (`docs/ELECTROMECHANICS.md` 5.7).

## Screening verdicts that stand

`docs/SCREEN.md`, consolidated. The model **cannot decide polar vs antipolar**
for most chemistries: gaps sit inside its own 0.27 kcal/mol per monomer
resolution, and Dipole's three tiers (ours, MACE, vertical DFT) disagree on
the VDCN copolymer's sign, all within that resolution. Accessibility of the
polar phase is the one criterion that resolves; by it TrFE leads, which is
what industry uses (a positive control, not a discovery). Nitrile chemistries
(AN, VDCN, FANOME) cannot reach a polar all-trans phase cheaply and are marked
not rankable on polarity. The "nitrile chains misbehave" narrative was
**withdrawn**: it rested on a frozen-angle fit artefact (fixed:
`fit_ris(angles="relaxed")`, `docs/NITRILE_LANDSCAPE.md`) and on Dipole's
finite-chain gates, whose references were never at minima; on requalified
references the chemistries behave alike and the free-azimuth boundary is what
fails. Ewald is implemented and validated (Madelung to every digit); it changed
no polarity verdict, truncation was never the cause.

## Phonons (new, 2026-09-13)

`src/polyfind/phonon.py` gives Gamma-point lattice dynamics on the packer's own
potential: an all-atom energy assembled term for term (known-answer match to
1e-14 against `packer.energy`), analytic gradient, finite-difference Hessian,
acoustic sum rule to 7e-8 without imposing it. `docs/PHONONS.md` has the
results: at the all-atom minimum with stretch minima pinned to built bond
lengths, beta's softest modes are rigid-chain transverse librations at 34 to
49 cm^-1 (~1 THz), alpha 39, gamma 22. No torsional stiffness enters (the
packer's torsion term is a constant), so twist modes are lower bounds; the
constrained references are not all-atom stationary points (fitted stretch
r0 differ from built bond lengths by up to 0.12 A) and show imaginary
rigid-chain modes there. Reading for the requirement: the harmonic Gamma scale is
~1 THz, and that is all it says; it does not establish the absence of an
intrinsic sub-THz limit (Dipole's correction, accepted: no rates, linewidths,
field coupling or switching dynamics are in a Gamma curvature). Pinned-stretch
results are a changed Hamiltonian, labelled separately. Next validation: far-IR/Raman lattice modes of beta-PVDF, and
Dipole's bulk periodic curvature. `examples/crystal_phonons.py
--polymorphs beta,alpha,gamma` reproduces (gamma takes ~10 min).

## Dipole's open gates (theirs, not ours to chase)

- Berry polarization: **matrix running, three diagonal columns compared.**
  Their reducer is Vanderbilt-proper (`e = dP/de + d_jk P_i - d_ij P_k`); ours
  is dipole per reference volume, which lacks `- d_ij P_k` and so differs by
  P_y on the polar-strain column only. Clamped-ion e_y,jj, their frame, both
  P_y < 0, C/m^2, all `quantitatively_valid=false`: xx theirs -0.1675 / ours
  -0.2960; yy -0.1548 / -0.1559; zz +0.1125 / -0.0108. xy shear (x comp.) theirs -0.2905 / ours -0.3272; yz and xz shears
  are not expressible in our cell. Agreement on the polar column and the
  in-plane shear, a 0.12-0.13 gap on both transverse normal strains. Matrix
  is 25/25 complete on their side, still `quantitatively_valid=false` (xz
  amplitude gate fails on a ~1e-8 vector). Ours regenerated
  2026-09-13 by the scratch pattern around
  `examples/internal_strain_jacobian.py::clamped_ion` (the provider record
  file is gone from sarco, so the full example exits 2; call `clamped_ion`
  directly). Our reference P in their frame: charge-only -0.1187, induced
  -0.0248, total -0.1434. Reproducible record with geometry and source
  hashes: `deliverables/clamped_ion_columns/` from
  `examples/clamped_ion_columns.py` (Dipole asked for it; delivered 2026-09-13). **Do not fit to any Berry number until accepted.**
- **Open correctness item, ours:** `mechanics.piezoelectric` defines `e` as
  `(1/V0) dmu/de`, not Vanderbilt-proper; e_33 (polar strain) is off by P_y
  (~0.14 C/m^2, larger than our whole e_33 of 0.096) and the converse route
  agrees with it, so the field coupling carries the same term. Derive from
  the energy functional which definition the measured d33 is, then decide
  whether d33/d31 change. Secondary now that magnitude is not the gate, but
  it is an acceptance quantity and must be resolved before any d is quoted.
- Born/dielectric ladder: 4x8x16 passed its raw sum-rule gate (max 0.0065 e);
  the ladder as declared cannot pass because 4x8x8 failed; a new predeclared
  ladder may follow. Transverse eps (2.253, 2.235) is stable; chain-axis
  converging toward ~2.45.
- VDCN copolymer: polar fixed-cell stage completed under PBE-D3 (force 0.0042
  eV/A); antipolar and full-cell stages pending. Our starts (`deliverables/`)
  are what they are relaxing; the staggered variants there did not help and
  should not be spent compute on.
- Held-azimuth finite-chain ladder: size gates still 0/18 even with azimuth
  held. Finite chains do not converge; the periodic route is the only bulk one.
- Transverse stiffness anomaly (their 112 GPa vs our 17 and published ~20):
  still unisolated, queued behind other work; not on our critical path.

## Discipline that has kept this honest

Twelve-plus retractions across both sides, every one recorded in place rather
than overwritten. Rules that earned their keep: compute the like-for-like
quantity before touching a parameter; check a known-answer case first and
suspect the check before the claim when it fails; keep acceptance quantities
(d33, d31) out of every fit; state parameter count against observation count;
phrase results under the producer's declared bounds; never contract provisional
reference data as if validated. When an agent reports a plausible number
traceable to an invented constant, refuse it.

## Housekeeping

- Tests: `export PYTHONPATH="$PWD/src"; python -m pytest tests -q -p no:cacheprovider`, 599 pass, 5 skip (GPU). ~10 min. Dipole's runtime reports 3 last-digit
  frozen-literal failures in `test_mechanics.py` (PE/alpha/gamma energies at
  1e-15 relative) that do not reproduce here at any commit; machine noise,
  not a defect. A 1e-12 relative tolerance is the fix if ever needed.
- **Never set `POLYFIND_TABLE_CACHE=1`**: it is read as a directory name and once committed 85 MB of caches. It is gitignored now.
- Agents: launch in worktrees; they often pause on their own background jobs and report "waiting"; that is not stuck. Merge, run the suite, push, then remove the worktree.
- Commit attribution: `Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>` and the session URL, per the current system reminder.

## If continuing

Nothing on our side is unblocked. The xx column says the answer is "both":
a real electronic clamped-ion term of ~0.13 C/m^2 and an overstated target.
When the yy column and the accepted matrix land: write the scope boundary
into `DESIGN.md` and `docs/BENCHMARK.md` with the measured size of the
electronic term, and restate the piezoelectric shortfall against the accepted
Berry-derived total rather than the failed sweep's 0.4-0.6. Either way, reply in `docs/REFERENCE_DATA_REQUEST.md`, then re-run
`examples/electromechanics.py --preset flux` if any parameter changes.
