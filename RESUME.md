# Resume: polyfind, state as of 2026-09-13 (late)

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

## Dipole's open gates (theirs, not ours to chase)

- Berry polarization: **route repaired, matrix running, first column in.**
  All nine directional zero solves passed the 2e-7 gate; the 25-point
  clamped-ion matrix launched 2026-09-12 21:23 PDT. Signed xx pair complete:
  their preliminary proper clamped-ion e_y,xx = -0.1675 C/m^2, lab derivative
  +0.0387 (identity closes on their P_y -0.2062). Ours, regenerated
  2026-09-13 in their frame with the same sign of P: -0.2960. Difference 0.13
  C/m^2 (~5 pC/N): the electronic term we lack is real but about a quarter of
  the old 0.4-0.6 target, so both hypotheses are partly true, on one column,
  `quantitatively_valid=false`. The yy column (polar strain, the d33 one) is
  next. Recompute ours with the scratch script pattern in
  `examples/internal_strain_jacobian.py::clamped_ion` (the provider's
  `internal_strain_geometry.json` is no longer at its default path in sarco,
  so the full example exits 2; call `clamped_ion` directly). **Do not fit to
  any Berry number until they declare one accepted.**
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

- Tests: `export PYTHONPATH="$PWD/src"; python -m pytest tests -q -p no:cacheprovider`, 584 pass, 5 skip (GPU). ~9 min. Dipole's runtime reports 3 last-digit
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
