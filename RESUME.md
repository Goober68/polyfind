# Resume: polyfind, state as of 2026-09-13 (afternoon)

**Final status of the 14-hour window, 2026-11-07 19:06 UTC.** Closed today: the stage-7 problem solved by bounded
least squares (`pvdf-dft-valence-belt3`, 5.18-5.19); the sensitivity, screen and basin analyses (5.15, 5.17, 5.20); the
Gamma-route relaxed-ion tensor (5.16); the pendant-group finding and its fix (5.22: the rigid refinement cannot tilt the
CH2/CF2 groups as DFT's alpha does, the gap row carried 0.23 kcal/mol per monomer of that strain, `polyfind.pendants`
relaxes the groups as rigid tilts, the refit is `pvdf-dft-valence-belt4`, the sixth prediction set
`e1/e2_polyfind_belt4.json` carries the correction, the copolymer model reorders to TGTG' lowest); the soft C-C stretch
constant quantified on the Gamma route (beta unmoved, alpha's polar eps_0 3.47 -> 2.93 floored); stage 8 (shape term,
5.21) and stage 9 (tilts as refinement variables, 5.23) designed. Open: the lateral coupling (C_12 0.87, the one residual
past a tolerance at belt4), waiting on Dipole's rotated-setting cell for its anchor; Dipole's alpha and Gamma successors
(launched 05:10 UTC, no result by 19:06) -- when they land: `examples/dft_cell_geometry.py` on the alpha record for the
backbone angles and the pendant tilts (prediction CF2 wag 2.0-2.5, rock 1.1-1.4, twist 0.6-0.7 deg), the Gamma DFPT
against 2.23/4.61/2.73 and the Born charges, then re-evaluate belt4 and refit only if a residual leaves 1.5 tolerances
(task #16). Suite 724 passed, 5 skipped. Branch `claude/polymeric-stable-arrangements-uh06b1`, head 7f68362 before this status commit.

**Dipole, 2026-11-08 22:33 UTC.** A power outage on Dipole's host killed both DFT successors (alpha at BFGS 37/120
after 27 h, Gamma at representation 22/36); Dipole marked them lost and started fresh ones under the PI's standing
restart instruction, no checkpoint reuse: `pvdf_alpha_ref_fresh_20261108_v2` (ETA uncertain; the v1 pace says about
three days) and `beta_zero_gamma_fresh_20261108_v2` (40-46 h). Hourly lines on sarco main. Nothing from Dipole lands
before 2026-11-10; the stage-8 items that wait on it (task #16) stay parked.

**Correction, 2026-11-08 (read before the candidate paragraph below).** The candidate screen's first pass ran one
of an isotactic sequence's four TG conformers (TGTG'); for every chlorine and nitrile row the other handedness
(TG'TG) is 0.9-1.5 kcal/mol per monomer lower (a multi-start probe confirms it is the model's effect, not the
search's), so the defect rows' margins are +0.57 to +0.97, not near zero, and the defects read as TG-favouring
pins (CFE +7.9, FAN +7.4, CTFE +6.1, AN +3.9 kcal/mol per unit; CDFE does not pin) rather than margin tuners. The
near-degeneracy comes from the TrFE matrix (29-40 mol%), and a fourteen-monomer cell at the film's composition
(64/29/7, one CFE per fourteen) gives +0.11 where the ten-monomer 70/20/10 cell gives +0.97: the model reproduces
the relaxor at the real composition once the cell holds the matrix between the pins. `docs/CANDIDATE_MATERIALS.md`
sections 2-8 are rewritten on this (the 70/20/10 ranking is withdrawn; the kink-pair bracket of section 8 stands
since it does not depend on the margin). `examples/candidate_screen.py` carries the four conformers
(`TG_VARIANTS`, `conformation()` sized to the repeat, `derive()` takes the lowest TG). The 14-mer dynamics: the
71/29 matrix is degenerate (-0.05) and frozen at the drive (kink pair 5.5, channel 0.8 kHz, saturated at the bias);
one CFE per fourteen holds it off saturation (polar fraction 0.24) and cuts the kink to 1.8 (channel at the drive),
giving tan delta 0.07 for 0.125% AC strain at the operating point (`candidate_frequency_map_belt4_anchor14.md`):
the relaxor mechanism in one cell. The 14-mer matrix is degenerate within resolution from 29 to 36 mol% TrFE
(71/29 -0.05, 64/36 -0.06); CDFE at one per fourteen pins harder than CFE (+0.40 against +0.11; pin strengths are
cell properties: CFE +7.9 per unit in the 20 mol% ten-monomer cell, +2.3 in the 29 mol% fourteen-monomer one, CDFE
-0.3 and +6.3) with the fastest lattice (kink 1.1, channel 1.3 MHz) and no strain at the bias; CFE+CDFE per
fourteen adds the two pins (+0.64), so the pin-and-loosener hypothesis fails as posed. The kink pair's handedness
splits 5-9 kcal/mol in the pinned lattice and 0-0.5 in the relaxed one (section 8). The 14-mer brackets with the
pinned-from-own protocol (`--pinned-from-own`; the pinned-from-all-trans start sticks for long repeats): 64/29/7
1.8-3.5, matrix 5.5-7.6; CDFE-14 kink 1.1, CFE+CDFE-14 0.86 (loosest, hardest pinned), AN-14 2.0 with the cell
TT-ground (-0.15: AN favours the all-trans lattice at the film's dilution, FAN about +0.25 with tilts estimated),
so neither nitrile is the relaxor's pin. Deferred: the kink-position scan (the pinned protocol it needs exaggerates
handedness and sticks; the own basin cannot compare positions), the dilute kink (a mixed-conformation cell), the
migration barrier. CTFE at one per fourteen pins as CFE does (+0.135) with less stroke and an all-trans lattice
metastable against kinking (kink -0.6: dynamics undetermined). The pin density is the dial: margin -0.05 + 2.3/N for
one CFE per N in a 30 mol% matrix (one per twenty +0.055, kink 2.3, channel 157 kHz, 0.154% strain at tan delta
0.20; one per fourteen +0.11, kink 1.8, 364 kHz, 0.125% at 0.07), so under the objective the recommended controlled
sequence is a 30 mol% TrFE matrix with one CFE every fourteen to sixteen monomers, strictly periodic
(`docs/CANDIDATE_MATERIALS.md` section 5). All background computations of 2026-11-08 are complete.

**Where things stand, 2026-11-07 (read this first; the detail is in the dated paragraphs below and in `docs/CRYSTAL_STAGE.md`).**
Candidate pass, 2026-11-07 22:50 UTC (`docs/CANDIDATE_MATERIALS.md`, `examples/candidate_screen.py`,
`examples/candidate_loss_map.py`): ten-monomer sequences packed in TTTT and TGTG' under belt4 with pendant
tilts. TrFE sets the TT-TG margin in steps of about 0.3 kcal/mol per monomer per ten mol% (PVDF +0.74,
80/20 +0.18, 70/30 +0.15, 60/40 TT ground); one CFE per ten beside 20 mol% TrFE brings it to +0.03 (the
degeneracy the relaxor needs), CDFE leaves it, CTFE overshoots to TT ground, CFE without TrFE does
nothing; placement of the CFE moves it by 0.03 per monomer (the local-margin spread a controlled sequence
sets by choice). Loss per strain is not discriminated by this pass (tier 3, the barrier per chemistry,
is the gap). Ranking in the document's section 5; chlorine rows are extrapolated chemistry. Section 6 reads
Dipole's slate of other chemistries (VDCN, CNEPO, FANOMe, SCN, VDC, AN ...) against the same screen and against
Dipole's own registry gate and barrier table: the nitrile homopolymers' all-trans crystals are 6-9 kcal/mol per
monomer above their TGTG' forms (not conformational switches), a nitrile defect beside TrFE drives the lattice to a
saturated TT ferroelectric (the opposite of CFE), and Dipole's barrier table is single-bond torsional barriers (on
its loss term every slate candidate scores 2-2.6 times better than the baseline), which Dipole's own reading says do
not decide the loss: the field moves them 1-4%, switching is collective, and the collective number is section 8's. Nothing on the slate
displaces the VDF/TrFE/CFE family for this objective on current evidence (`examples/slate_screen.py`). Section 8
(`examples/kink_formation_screen.py`, `candidate_frequency_map.py`): the kink-pair formation energy in the packed
cell, the first piece of the dynamics the statics could not give, as a bracket between the relaxed kink crystal
(own packing basin) and the pinned all-trans lattice (`refine_crystal(cell_free=False)`): 4.3-9.3 kcal/mol in PVDF,
4.9-7.6 in the TrFE copolymers, 1.2-3.1 in the CFE terpolymer model, 0.4-0.6 with CDFE, so the copolymer-terpolymer
separation (3.7-4.4, a factor 500-1600 in rate) holds at both ends: the terpolymer's channel follows a 34-66 kHz
drive and the copolymer's is frozen there; loss per strain ranks CDFE then CFE at both ends. CTFE (-0.4 to +20.4)
and the nitrile rows (1.3-11.3, 0.4-12.2) are undetermined: free, their lattices relax to a kinked packing; pinned,
they have no room for a kink. Running: the pinned-then-released protocol (`kink_formation_belt4_pinned_released.json`,
the kink crystal in the basin nearest the all-trans cell). Open: the dilute kink (one kinked chain among all-trans
neighbours, a mixed-conformation cell the packer does not have) and the migration barrier.
The mission's two stages are closed and documented: the large-cell response path is 10x faster
(`docs/PERFORMANCE_REVIEW.md` 16), and the crystal terms of the valence potential are fitted to
Dipole's DFT anchors through seven stages. The current preset is `pvdf-dft-valence-belt3`
(fluorine belt sites, the PBE-D3 chain geometry, a fitted backbone-angle restraint; fitted by
bounded least squares, objective 9.11): against Dipole's records beta's cell is -0.1/-0.4/+0.7%,
the relaxed-ion stiffnesses 1.06/1.01/0.91, the Berry polarization 1.00, the polar Berry slope
0.92, the electronic permittivity 0.93-0.95, the alpha-beta gap and the backbone angles within
1.24 tolerances; its E1/E2 sets are `e1/e2_polyfind_belt3.json` (sections 5.18-5.19). The
pre-registered set on `pvdf-dft-valence` still stands as the primary one. Found 2026-11-07 (5.22):
the refinement's rigid CH2/CF2 groups cannot tilt as DFT's do (Dipole's alpha start geometry has
CF2 rock 2.3, wag 2.9, twist 1.3 deg; beta none by symmetry), and under belt3 alpha releases 0.23
kcal/mol per monomer through rigid group tilts alone, so the gap the fit matched carries that much
pendant strain; `polyfind.pendants.tilt_relaxation` evaluates it, `crystal_stage.py --pendant-tilts`
adds it to the gap row, and the refit from the belt3 point is the current preset
`pvdf-dft-valence-belt4` (objective 11.11; the gap back on the anchor with the relaxation inside
the fit, paid by the fluorine well -12%, the H-F pair +8% and the lateral coupling, C_12 now 0.87 of
DFT's and the largest residual at 1.29 tolerances; everything else within 1.1). The stage-8 shape
term starts from belt4. Prediction on record for Dipole's native alpha positions: CF2 wag 2.5, rock
1.1, twist 0.7 deg. Open on our side: the
lateral coupling (C_12 at 0.91, the lateral axis -0.4%) is held by the other anchors and the prior
together and no radial term adds a direction (5.15, 5.18, 5.20); a shape term (hydrogen belt
sites or a CF2 quadrupole) is the next form change and waits for Dipole's rotated-setting cell
(5.17, requested). Open on Dipole's side: the Gamma DFPT and alpha successors are running after a
WSL shutdown killed the first attempts; polyfind's like-for-like prediction for the DFPT tensor is
eps_0 = 2.23 polar, 4.49 lateral, 2.73 chain (5.16, 5.19), posted on the exchange; belt4 gives
2.23, 4.61, 2.73 (5.22), and it does not move when the valence fit's soft C-C stretch constant is
floored at a physical value; alpha's Gamma-route polar eps_0 does (3.47 to 2.93), so the floored
value is the one to quote for alpha (5.22). Suite: 724 passed, 5 skipped (2026-11-07 18:03 UTC). Fitting from
here: `examples/crystal_stage.py --method trf` (one evaluation per variable per iteration, 5.18).

2026-11-04 research note (cloud session, PI-nominated the same day to drive the electrostriction
research): the literature pass is `sarco/materials/computational/papers/2026-11-04_electrostriction_literature_pass.md`
(the giant strain is a field-induced TG/T3G-to-TT switch shared by crystals, OAF and tie molecules;
the relaxor dispersion is interfacial; loss, frequency and strain are one relaxation, with Ea 0.12 eV
and Tf the discriminating number); the plan is `docs/ELECTROSTRICTION_RESEARCH.md` (tiers: intrinsic
M/Q, switching margin and transition strain, switching barrier and tau(T), then the deformable
connectors), and the request to Dipole is `sarco/materials/ask-polyfind-dipole-electrostriction-verification-v1.md`
(E1 DFPT at strained cells, E2 two-phase PBE-D3 cells per chemistry, E3 one NEB, E4 the bench sweep
in Vogel-Fulcher form). New: `polyfind.electrostriction` (intrinsic electrostriction by the strain
and the polarization routes, agreeing to 5e-5 on force-shifted rigid beta; `TwoState` switch,
`transition_strain`: TGTG' to TTTT is +13.4% axial at the ideal angles), `examples/electrostriction.py`,
`tests/test_electrostriction.py`; `pack.polish`, `mechanics.relax` and `_state_at_strain` take
optional `gtol`/`ftol` (defaults unchanged). The PI's reminder that this seat also supports the
microreactor (fluoro-ribosome) programme is answered in section 10 of the plan: `polyfind.chainend`
gives the torsion-gate field per unit type (VDF end: 1.6 V/nm to write gauche, 0.2 to hold trans)
and beta's registry corrugation (0.48 kcal/mol per monomer, 41 pN per monomer peak) for D21/D24;
`examples/chain_end_gates.py`, `tests/test_chainend.py`. Level 1 of the loss model (PI: "go with level
1") is `polyfind.loss`: a bias-tilted Vogel-Fulcher switching channel giving tan delta, the strain per AC
field and the cycle efficiency at any (f, T, E_bias); on published terpolymer inputs the bias-hold point
at 34-66 kHz carries tan delta 0.06-0.12, not the assumed 0.0125, and crosses the 4 MW/m^3 heat budget
near 25 kHz (plan section 11; the bench sweep in VF form replaces every input). Level 2 ("Go") is
`polyfind.kmc`: the units on a lattice with coupling and random margins, rejection-free kinetic Monte Carlo
under a DC plus AC field, barrier and strength calibrated by measurement, small-signal permittivity at any
bias from the lattice's own fluctuations; first results (plan section 12; 20 x 20 and 32 x 32 agree): the bias axis has a knee near
15 V/um on the first parameters, the lattice's apparent activation (0.43 eV) is a fifth of the way from the bare
barrier to the VF slope. The three-state unit (the polar state's dipole points either way; margin `kT ln 2` for a
zero-field polar fraction of one half; the lattice calibrated on the centroid of its own loss peak) removes the
two-state lattice's over-saturation (strength at 55 V/um 45% of zero bias, measured 50%), puts the zero-bias strain
at 2f, and makes the stroke-at-budget surface rise from zero to a maximum near 20 V/um and decline to the bias-hold point
(by a third on the 32 x 32 lattice, six tenths on the 20 x 20) on the first parameters, the position set by the
polar fraction's rise with DC field, i.e. by the bench's quasi-static S(E_dc) curve (E4 item 7); the bias-hold
point is a question for the calibrated lattice. The coupling x spread grid
settles the temperature slope: 0.36-0.37 eV apparent at every disordered point against VF's 1.0, the lattice ordering
into a ferroelectric past J of about 0.03 eV, so the VF slope is an input (`at_temperature(T, tau_target=...)`
recalibrates the barrier at each T to the bench fit) and the lattice's deliverables are the bias and amplitude
dependence at each temperature; the spread is the knob for the suppression (section 12, last block). The PI runs the E4 bench items personally once an LCR meter is in hand (sarco ask, E4 status line).
Needed from the PI: push `kink-switching` and
`vector-repeat-merge`, which exist only on `E:\source\gh\polyfind`.

2026-11-05 update (the mission: a solver an order of magnitude faster than DFT at DFT-like
accuracy, worked in two stages; all on `claude/polymeric-stable-arrangements-uh06b1`).
**Speed** (`docs/PERFORMANCE_REVIEW.md` section 16): the 120-atom copolymer E1 job, which did
not finish in two hours, runs in 359 s on four cores (1,139 s serial), results identical;
neighbour-list Ewald real space and a BLAS reciprocal tensor (1e-13 agreement), one BLAS
thread per process (OpenBLAS's pool made a 360 x 360 solve 200-500x slower), independent
relaxations in forked workers (`mechanics.parallel_map`, `POLYFIND_WORKERS`), no Coulomb
work in the pair kernel under Ewald; serial profile 314 -> 69 s. Left: the dense
Lennard-Jones kernel, the serial reference build. **Accuracy** (`docs/CRYSTAL_STAGE.md`,
preset `pvdf-dft-valence-crystal`, opt-in): the eight nonbonded terms the conformer trainset
left free were fitted to Dipole's QE beta cell, the antipolar alpha cell, the two margins and
the Berry-phase P. Fixed-charge: P 0.115 -> 0.174 (DFT 0.18), alpha's long axis -6.2% ->
-2.7%, beta's polar axis only -4.5% -> -3.7% with the other axes drifting out by about one
tolerance; the isotropic Lennard-Jones form cannot open the polar H...F contact without
opening the lateral F...F one. With induced dipoles (the E1 setting) the polar axis stays at
-4.8% and P overshoots to 0.214, so the pre-registered set on `pvdf-dft-valence` remains the
better E1 prediction and both sets are kept (`docs/predictions/e1_polyfind_crystal.json`,
`e2_polyfind_crystal.json`). Next, in order: fit in the induced-dipole setting (plumb
`polarizable` through `fitting.predict`, or an antipolar branch on the mechanics path) with a
pair-specific H-F term, adding Dipole's six strained cells as anchors when they land.
**Afternoon, 2026-11-05.** Dipole's six +/-1% axial beta cells landed (exchange main 54897c21):
from their stresses the relaxed-ion stiffness is C_11 22.8, C_22 24.5, C_12 5.9, C_33 302.6 GPa
(record `beta_qe_axial_v1`); the pre-registered prediction has C_11 to 1%, C_22 24% low, C_12 47%
low. Stage 2 of the crystal fit (`docs/CRYSTAL_STAGE.md` section 5; the fit path now takes
`polarizable=`, `lj_cutoff=`, a pair-specific H-F Lennard-Jones term through `FFParameters.lj_pairs`,
and the stiffness as anchors) ran in the setting the predictions are made in: P and the margin
inside tolerance, beta's polar axis -5.0% -> -3.9% and no further, lateral +1.5%, C_12 still 57% of
DFT; preset `pvdf-dft-valence-crystal2`, prediction files `e1/e2_polyfind_crystal2.json`, PE's cell
degrades (not an anchor: a warning). The polar-axis scan says why: at the DFT spacing the
electrostatics pulls the axis in by 1.3-2.7 kcal/mol per A and the fits can only inflate the
repulsion. Charge penetration (`EwaldSpec(charge_damping=b)`, exclusion-consistent, tested) and
the PBE chain geometry (C-F 1.386 against the built 1.350) were each probed: -5.0% -> -4.3% and
-4.4% respectively, neither the lever alone; stage 3 (eleven variables, charge penetration free)
ended at the same floor as stage 2 with the damping off (section 5.6): three forms share the
floor, so the next lever is the fluorine's off-site charge reaching the lattice as a site, with the
PBE geometry and Dipole's DFPT tensors as anchors. **E2 correction**: the E2 script omitted the refined
backbone angles' strain energy, favouring the gauche conformers by about 0.45 kcal/mol per
monomer (PVDF TTTT margin +1.30 -> +0.85, copolymer +1.09 -> +0.46); fixed, files regenerated,
Dipole told.
**2026-11-06.** Stage 4 of the crystal fit broke the floor (`docs/CRYSTAL_STAGE.md` 5.7-5.9):
five probes showed that nothing isotropic moves beta's a:b ratio (0.530 against DFT's 0.558) and
that flattening the fluorine's repulsion in the plane perpendicular to C-F does; implemented as a
rigid pendant fragment with two massless, chargeless, dispersion-free "belt" sites 0.7 A either
side of each fluorine (pseudo-element X, `polymers.belt_fluorine`/`pvdf_belt`/`belt_variant`),
a repulsion-only X-H Lennard-Jones wall whose depth is the twelfth variable, inert within its own
chain so the conformer energetics are untouched. Result: objective 65.1 -> 16.3, every cell axis of
beta and alpha within 1% except the two chain axes (+1.0%, +1.4%: the backbone angle the lattice
settles at), P 0.179 (fitted to 0.18; Dipole's Berry 0.193 is the anchor from stage 5 on), margin
-0.760 (-0.726), C_11 25.1 (22.8), C_22 24.6 (24.5), C_12 4.9 (5.9) GPa. Preset
`pvdf-dft-valence-belt`, meaningful only with `--variant belt`. Dipole's Berry repeat landed
(0.19269 C/m^2 at the reference cell; the pre-registered P0 is -27% against it, stage 2 -2%, stage
4 -7%). Next: the belt preset's E1/E2 sets and their cross-check, stage 5 with the 0.193 anchor
from the stage-4 point, the PBE-D3 chain geometry on the belt chain, DFPT tensors as anchors.
**Later, 2026-11-06.** The belt preset's E1/E2 sets are written (`e1/e2_polyfind_belt.json`,
`docs/CRYSTAL_STAGE.md` 5.10; the copolymer E1 takes 1,356 s with the sites): beta against
Dipole cell -0.6/-0.1/+0.8%, C_11 1.10, C_22 1.00, C_12 0.83, Berry P 0.93, eps_inf 0.93-0.94 of
the DFPT electronic tensor; E2 PVDF TTTT margin +0.758 (DFT +0.726). Beta's lateral relaxed-ion
permittivity keeps rising across the stages (2.83 -> 3.25 -> 3.82; the chains' libration) and is
the set's sharpest number for Dipole's DFPT eps_0 to test. Cautions: the copolymer packing lands
in a different minimum (gamma 88 vs 92 deg, dipole tilt 38 vs 23 deg), so its rows compare
minima; PE's shear stiffness is 5x the pre-registered value (C_66 11.5 GPa), so a PE anchor comes
before the belt terms are trusted on a hydrocarbon. Stage 5 (the Berry anchor 0.1927, from the
stage-4 point; `crystal_stage_v5.json`) is running; `examples/crystal_stage.py --geometry pbe`
builds the belt chain at Dipole's PBE-D3 geometry for the stage after. Dipole's three-point polar
Berry pair landed (exchange main 6574b77d): d|P|/d eta_polar = -0.326 C/m^2, the relaxed-ion e_11
on one branch; polyfind's sets are 0.71 (pre-registered), 0.96 (stage 2), 0.87 (belt) of it and
0.98/0.98/0.93 relative to their own P (record `beta_qe_berry_axial_v1`, cross-check lines added;
candidate anchor for the next stage). Stage 6, that anchor on the belt form (`crystalfit.polar_slope`,
`STAGE6`, `--stage 6`; the slope comes from the stiffness block's strained states for nothing and
matches the E1 path to 1e-4), is implemented and queued to start from the stage-5 result
(`crystal_stage_v6.json`); write-up `docs/CRYSTAL_STAGE.md` 5.11.
**Exchange**: the Science 2022 main text makes FA the dehydrochlorinated CFE pair, i.e. the
same `-CH=CF-` unit as DB (the CF3 reading is withdrawn; E2 is six chemistries; record and
devstream amended, sarco a498326c). Dipole at 05:12 UTC: no accepted strained units yet, one
axial relaxation in SCF, alpha relaxing. The E1 report now carries the copolymer's
electrostriction column (intrinsic lattice Q about +3 m^4/C^2 along the polar axis). The first stage-5 run was cut by the session's two-hour background limit at evaluation 292 of 320 (best objective 12.14 from 18.15), and its result file was only written at the end, so the fitter now checkpoints its best-so-far to the output file on every improvement (`crystalfit.fit(checkpoint=...)`, `"complete": false` until the final write; `--x0` resumes from it). Stage 5 landed at objective 11.64 (from 18.15 at the stage-4 point): |P| 0.188 against the Berry 0.193, the alpha-beta gap on the anchor, every residual inside 1.2 tolerances, the chain axes (+0.7%, +1.0%) the one systematic left (`docs/CRYSTAL_STAGE.md` 5.11). Stage 6 (both Berry anchors) landed at 11.19 from 12.81: |P| 0.192 (z -0.10) and the slope -0.299 (z 0.90), the cell and the gap unchanged, the carbon diameter's excess gone (x1.016); not registered, stage 7 decides. Found and closed on the way: the E1/E2 scripts' refinement ran under `refine_crystal`'s default potential (truncated Coulomb, no induced dipoles), leaving c 0.2% short and the angles 0.3 deg closed against the fit path; they now pass the full potential (`refine_kw`) and reproduce the fit's cell to four decimals, and take `--geometry pbe` for the belt chain (5.10, last paragraph). The belt sets predate the fix and are left as they are. The PBE-D3 chain geometry at the same point makes the chain axes longer, not shorter: the lattice opens the backbone angle on either geometry, against the refinement's invented 105 kcal/mol/rad^2 restraint, and with the angles frozen at DFT's values the alpha-beta gap flips sign: the gap anchor was being matched through the alpha chain opening its angles (`docs/CRYSTAL_STAGE.md` 5.12). DFT's angles (beta 114.63 planar; alpha 118.2 at CH2, 116.9 at CF2 from the PBE-D3 start cell) are now anchors and the restraint's stiffness a variable, on the parameter set (`FFParameters.angle_stiffness`, bound for the refinement by `applied()`): stage 7, `--stage 7 --geometry pbe`, ran from the stage-6 result (5.13). It landed at 11.65 from 25.14 at that point: the restraint fitted at 178 kcal/mol/rad^2, alpha's angles on DFT's (118.05, 116.94), beta's 116.0 against 114.6 (its chain axis, +0.75%, is now that angle: the C-C bond is Dipole's), both Berry anchors and the gap held, every cell axis within 0.8%, C_12 (5.0 against 5.9) and beta's lateral axis (-0.6%) the residuals shared with stage 6 (5.13 has the like-for-like table). Registered as `pvdf-dft-valence-belt2` (`fitting.py`; the thirteenth variable is the restraint, bound for every refinement by `applied()`), which means something only on the belt chain at the PBE-D3 geometry (`--variant belt --geometry pbe`); its E1/E2 sets are `e1/e2_polyfind_belt2.json` (5.14): beta on Dipole's records cell -0.3/-0.6/+0.8%, C_11 1.03, C_22 0.99, C_12 0.85, Berry P 1.00, slope 0.92, eps_inf 0.95/0.93/0.94; E2 PVDF TTTT margin +0.726 on DFT's. Beta's lateral relaxed-ion permittivity is 4.39 (the libration: 2.83 -> 3.25 -> 3.82 -> 4.39 across the sets while P is pinned, i.e. the belt forms imply half the isotropic forms' libration stiffness), the number for Dipole's DFPT eps_0 to test. Cautions: the copolymer packing landed in a third arrangement (gamma 92.1, dipoles along b, 85 deg from a; read its rows in the dipole's frame; E2 puts its three conformations within 0.06 kcal/mol), and PE's cell is 9% smaller than under the belt set with eps_inf 3.04 and C_11, C_22 32, 34 GPa: Dipole's PE cell is the anchor needed before the C/H terms are trusted on a hydrocarbon. Suite after the preset: 710 passed, 5 skipped. **2026-11-07:** Dipole's native alpha relaxation and the beta Gamma DFPT job (SCF and D3 done, PH partial, no dynmat) were killed by a WSL shutdown and are terminal FAILED with zero accepted units; Dipole asked the PI for a recovery route and the PI chose fresh successors, relayed on the exchange devstream with the polyfind seat's priority (Gamma eps_0 lateral, alpha angles, PE cell). Stage 8 waits on those; the FA/DB stereogeometry choice for E2 is still open on Dipole's side. Under the PI's 14-hour standing instruction (05:06-19:06 UTC, hourly check-ins as a scheduled task): (1) the belt sites are out of the Ewald sums and the dipole solve (`CrystalPacker._es_idx`; `VDS98["X"] = 0`; belt copolymer energy call 0.64 -> 0.42 s, predictions moved < 4e-6; `docs/PERFORMANCE_REVIEW.md` 16.6); (2) the sensitivity analysis at the stage-7 point (`examples/crystal_stage_sensitivity.py`, `crystal_stage_v7_sensitivity.json`, 5.15): six of twelve directions constrained, no residual addressable alone, and the lateral coupling (C_12, beta b) reachable only by trading the alpha-beta gap to +2.3 tolerances and beta a to +1.2 (verified by a real evaluation, objective 22.0): the form needs a term that couples the lateral axes without moving the gap; (3) the relaxed-ion dielectric tensor by the Gamma route (`examples/gamma_dielectric.py`, 5.16: the fit-path reference, valence terms pinned at it, the backbone angles on the preset's restraint, belt sites as affine combinations of their CF2 group, all-atom Hessian + Born tensors): beta under belt2 eps_0 = 2.24 polar, 4.57 lateral, 2.74 chain against the rigid-chain route's 2.12, 4.39, 2.45; the libration at 46 cm^-1 carries 2.38 of the lateral 2.48; this is the number for Dipole's DFPT, posted on the exchange as a prediction; (4) basins (`examples/copolymer_basins.py`, 5.17): refining every minimum of the packing screen instead of the first finds, on every fitted potential, a second packing of beta with the chains turned 60 deg (dipoles 119 deg from a, gamma 88.5) degenerate with the DFT setting to 0.1 cal/mol per monomer, and the same pair for the copolymer (dipoles 24 and 85 deg from a, degenerate to 0.1 cal/mol), so the copolymer rows of the sets compared arrangements the screen happened to pick; the chain-orientation energy is unconstrained by every anchor in hand, and the DFT energy of beta in the rotated setting was requested on the exchange (sarco dd4b146c; geometries in `docs/predictions/basins_pvdf-dft-valence-belt2_beta_basin{0,1}.xyz`). The belt2 copolymer E1 cell re-run with the inert-site exclusion: 776 s against 1,356 s, numbers within the relaxation noise. `refined_reference(start=, start_params=)` and `e1_predictions.py --basin/--start-params` refine any screen minimum or a given arrangement; the belt2 copolymer in the pre-registered set's arrangement (`e1_polyfind_belt2_tilt24.json`) lines up with that set row by row (5.17), so the sets' copolymer rows compare once the arrangement is named. The crystal fitter has a bounded least-squares method (`--method trf`, residual vector = anchor z rows + ridge components, forward-difference Jacobian, one evaluation per variable per iteration, parallelisable); the stage-7 run from the stage-6 point reached 9.36 in 235 evaluations (31 min; passes Powell's 11.65 at evaluation 144) and its continuation converged at 9.11 in 53 more (`crystal_stage_v7_trf.json`, `crystal_stage_v7_trf2.json`; 5.18): a better minimum of the same objective, every residual within 1.24 tolerances, C_12 z -0.92 (Powell's point: -1.51), gap -0.21, carbon diameter at the valence value, restraint 228; Powell's stage-7 point was local and 5.15's trade-off holds only there. Registered as `pvdf-dft-valence-belt3`; its E1/E2 sets (both copolymer arrangements), Gamma-route eps_0 (2.23/4.49/2.73, within 2% of the posted belt2 prediction) and Jacobian are written (5.18, 5.19): against Dipole's records cell -0.1/-0.4/+0.7%, C_11 1.06, C_22 1.01, C_12 0.91, P 1.00, slope 0.92; E2 PVDF TTTT +0.736; the copolymer's three conformations within 0.06 kcal/mol on both belt presets. At the belt3 point the Jacobian says closing the lateral coupling buys 1.0 in the anchors' sum of squares and costs 0.8 in the other rows (P, beta a, the alpha angles) plus 3.4 of ridge, and a half-weight ridge probe does not move (7.06 -> 7.00, restraint 228 -> 238 only; `crystal_stage_v7_trf_ridge05.json`): the residual is held by the other anchors and the prior together, so a new term is needed, not a weaker prior (5.18); the requested rotated-setting cell is the datum to build it against. Screened against that Jacobian (`examples/stage8_screen.py`, 5.20): the belt distance and diameter and pair-specific F-F, H-H, C-C, C-F, C-H terms all lie in the span of the twelve variables (outside part 0-3% of their length), so no radial term adds a direction; a new term must have a shape (hydrogen belt sites, a CF2 quadrupole) and should wait for the orientation anchor. `--method trf` is the fitter from here. Suite after belt3: 715 passed, 5 skipped.

2026-11-04 update (cloud session, on `claude/polymeric-stable-arrangements-uh06b1`): the five
CPU performance patches of `docs/PERFORMANCE_REVIEW.md` section 15 are applied (commits
31a2456, bdbcb8f, e364e0b, 1700515; section 15.7 has the measurements on the merged tree):
vectorised FFT screen, single-build enumeration, fast chain geometry, cheaper table build,
closed-form chain Jacobian.  Funnel here 61.9 -> 24.9 s (2.5x with the refinement row under
cProfile, 2.1-2.2x without; section 15.7).  Table values move by < 1e-3 kcal/mol (clear any
on-disk `$POLYFIND_TABLE_CACHE`); line-group end points are unchanged to 2e-9; free-
parametrisation end points can stop elsewhere along their flat valleys (objective within
6e-4 kcal/mol per monomer, lattice energy and bend restraint moving up to 1e-2 each in
opposite directions; `refine_crystal(gradient="fd-geometry")` reproduces the previous
route).  The GPU port (G6) stays a patch.  IP: the solver's methods were assessed for
invention disclosure and found to be known art or implementation; the record is
`sarco/shared/ip/polyfind-methods-prior-art-2026-10-16.md`, and one consolidated defensive
row, row 88 of `sarco/shared/ip/idf-candidates.md`, was logged on the PI's instruction.

2026-10-09 update (on `claude/polymeric-stable-arrangements-uh06b1` unless named):
**the polar/antipolar margin is computed** (`docs/POLAR_MARGIN_DFT.md`): on Dipole's own QE
Hamiltonian (PBE-D3(BJ), SSSP precision, 90/360 Ry), fully relaxed, PVDF's gamma-free
antipolar cell is 1.11 meV (0.026 kcal/mol) per monomer BELOW beta, converged to 0.05 meV;
the sign is a +2.2 meV PBE vs -3.4 meV dispersion balance, so "degenerate" is the robust
reading. Our potential agrees to 0.02 kcal/mol. So for PVDF the static margin does not
separate hard from relaxor; barriers do (kink work below). Dipole's CP2K beta geometry is not
a minimum of their QE Hamiltonian (0.40 eV/A, -13.4 kbar, +8.7 meV/monomer); their Berry/Born
columns are evaluated there. **Far-IR** (`docs/PHONONS.md` last section, `polyfind.infrared`):
the "alpha" phonon row was the delta cell (P2_1cn); see the Phonons section. **Tests**: golden
values to 1e-12 relative, digests via stored arrays; 623 pass / 5 skip. **d figures** in
SCREEN.md/BENCHMARK.md rerun under the proper `e`. **Dipole's `physics-vector-repeat`** is
merged on the side branch `vector-repeat-merge` (8217385) only: the complete placed-cell owner
makes `CrystalPacker.energy` 6-11x and `energy_and_grad` 3-4x slower (energies equal to the
ulp), and Dipole's hold stands; asked in the exchange (consumer note 2026-10-09) whether they
accept a batched path held to the owner by a known-answer test. **Switching energetics**
(branch `kink-switching`, d6ce4f6/264115f off 8217385, `docs/SWITCHING.md` there; needs the
merged torsion): a beta kink is compact (180 deg within 1-2 monomers, lowest on the CF2-CH2
bond); a kink row costs 14.4-15.4 kcal/mol (0.63-0.67 eV) per kink, ~60x the resolution but
moved up to ~4 by the stretch choice; a lone kink pair untwists in a 2x2 host (no converged
lone-kink energy); the row hops one monomer over 8.22 kcal/mol (0.357 eV, first-order saddle,
L = 12, fixed cell only); the neutral 180-deg wall is -0.16 mJ/m^2, unresolvable; harmonic TST
with the model's Vineyard prefactor (4.78 THz) gives ~5 MHz per hop at 300 K, an estimate. On
this potential reversal goes by sheets, not single chains. The pinned host itself drifts (C-F
1.400 vs 1.350 built), a known limit carried by every number. That branch also adds an opt-in
`pair_list` to the complete owner, 7-8x faster at 192 atoms, default bit-for-bit unchanged.

2026-09-29 update (all on `claude/polymeric-stable-arrangements-uh06b1`, pushed): the `e`
definition item is closed -- `mechanics.piezoelectric` returns Vanderbilt's proper tensor and
its converse route holds the nominal field (`docs/ELECTROMECHANICS.md` 5.10; Born-fitted
d33 -9.86, d31 +0.42, not the cause of the shortfall). `phonon.relax_all_atom` owns its
acceptance and the screen phonons record is regenerated under it (tables unchanged).
`pack.polish` wraps `dz` inside its objective. `examples/film_d_budget.py`
(ELECTROMECHANICS 5.11): through this crystal's stiffness, d31 = +20 needs a chain-axis `e` of
6.3 C/m^2 (44x |P|), so d31 is a film-compliance property, not a crystal one; d33 = -32 needs
the polar column at -0.50 to -0.76 C/m^2 (2-3x ours), so the discriminating datum is the
provider's relaxed-ion polar column; the clamped-ion transverse differences move d by <1 pC/N.
New `polyfind.film` (ELECTROMECHANICS 5.12, `examples/film_response.py`): exact laminate of
crystal lamellae + passive amorphous layers along the draw axis. Film d31 +2.2..+6.3 (crystal
+0.42) from amorphous incompressibility; raising only the crystal polar column to -0.77
C/m^2 (3.2x) gives film d33 -32, d31 +12..+17.5, d32 +2.5..+3.1 (measured -32/+20/+1.5):
one crystal deficit explains both. Film inputs swept, never fitted; a measured draw modulus
and crystallinity for the same films would make it a prediction. BUT (5.13): the provider's own data
already estimate their crystal polar column at -0.34 (1.42x ours, via their larger Born charges;
kinematics agree, neither crystal relaxes internally under polar strain), which through the
laminate gives film d33 -14, d31 +6..+8: ~43% / 28-41% of measured, all signs right. The rest
(~2.3x) is outside the crystal (amorphous/interface/charge terms the laminate omits).
The d figures in `docs/SCREEN.md`'s response table and `docs/BENCHMARK.md` were rerun under
the proper `e` on 2026-10-08 (old values kept beside them); SCREEN.md's Step 3 correlations,
Rule 2 and candidates table still use the old `e` and are flagged in place.

2026-09-13 evening provider update: an isolated D: worktree now exists at
`D:\sarco-work\polyfind-vector-repeat`, branch `physics-vector-repeat`,
pushed `cb17be5`. It supplies vector-repeat placed-dipole geometry and a
complete six-column affine model record; see its
`deliverables/clamped_ion_vector_repeat_v2/README.md`. Full relaxed-ion
energy/forces/internal strain/C/S remain open. Do not merge source changes
into this C: checkout while the admitted Sarco native jobs depend on its
unchanged source. GPU is reserved; chemistry is CPU-only. The historical
platform/worktree statements below precede this update.

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
column that matters and it is the one the model cannot resolve (2026-10-09: for PVDF,
PBE-D3(BJ) puts it at -1.1 meV per monomer, degenerate, agreeing with the model, so the
discriminator moves to switching barriers); crystal
phonons (soft modes, lowest optical frequencies) are the intrinsic proxy the
solver can add. The tool computes no loss and no dynamics today. Recorded
in the exchange as the consumer note of 2026-09-13.

## Goal status, in one paragraph

**Speed: met**, well past target (`docs/BENCHMARK.md`): a full response tensor
in 0.3-7.8 s against the baseline's 87-143 s for a single finite-difference
derivative, and the baseline does no structure search at all. **Static side,
shipped preset `pvdf-dft-valence-flux-born` (rechecked 2026-10-08):** clamped-ion
dielectric within 5% of Dipole's DFPT (2.229/2.146/2.520 vs 2.235/2.253/2.601,
polar/long/chain), transverse Born charges within 0.19 e, C33 +6.4%, cell a/b/c
-3.6/+1.7/-1.2%; but |P| 0.143 is 19-24% BELOW the DFT 0.176-0.188 (the earlier
"within a few percent" was never true of this preset). Chain-axis components
remain open in Dipole's convergence ladder; every DFT reference is
`quantitatively_valid=false`.
**Piezoelectric coefficients: short**, crystal d33 = -9.86 vs measured -32, d31 = +0.42
vs +20 (Vanderbilt-proper since 2026-09-29; the earlier -8.8 / +0.02 were the
improper "film" bookkeeping, `docs/ELECTROMECHANICS.md` 5.10). Through the film laminate
(5.12) d31 becomes +2..+6; with the provider's estimated crystal polar column (5.13) the film
is d33 -14, d31 +6..+8, ~40% of measured with every sign right, and the remaining ~2.3x is
non-crystalline. The earlier account of the reason has been narrowed by elimination to one of two things
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
7. **Rigid chain rotation as a switching / loss proxy (2026-10-02)**: not usable in this cell.
   Clamped-cell in-phase rotation of beta's two chains: barrier 6.05 kcal/mol/monomer at
   40 deg, and the polar state stays a local minimum past 10 GV/m. Cell-relaxed continuation:
   1.5-2.3 kcal/mol/monomer but bound-limited and branch-tracking (not a true MEP), with
   1D-spinodal "coercive fields" of 14-17 GV/m that are just libration stiffness / dipole
   (the phonon column again). Real switching goes by kinks/domain walls along the chain,
   which a one-monomer, two-chain cell cannot hold. `docs/PHONONS.md`, last section.

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
results are a changed Hamiltonian, labelled separately. **First comparison with a measured
spectrum (2026-10-08, PHONONS.md last section, `polyfind.infrared`):** the "alpha" row above
is the packer's TGTG' reference, which relaxes to P2_1cn (delta-PVDF's polar arrangement),
not alpha's P2_1/c; the model's P2_1/c cell is 0.030 kcal/mol per monomer above it (inside
the resolution) with b 6.4% short. On the P2_1/c cell, the IR-active antiphase libration
(Rabolt and Johnson's 53 cm^-1 rotatory lattice mode, ~60 at 90 K) is 67.5 (68.1 with
Dipole's torsion), +8 high; nothing IR-active lies within 10 cm^-1 of the ~100/102 line;
the brightest mode is an internal Bu at 72.5 / 88.8 whose position is set by the torsion
term. No beta far-IR data below 100 cm^-1 were found. Still open: Dipole's bulk periodic
curvature. `examples/crystal_phonons.py
--polymorphs beta,alpha,gamma` reproduces (gamma takes ~10 min).

## The lattice-curvature column (new, 2026-09-14)

`examples/screen_phonons.py` and `docs/SCREEN_PHONONS.md`: Gamma-point phonons of every
screened chemistry's all-trans polar cell and exact antipolar cell, at gamma = 90 (the
screen's two cells, reproduced to the quoted digits) and with gamma free (the model's own
minima), each re-polished on the phonon packer and relaxed over every atom at fixed cell.
Record in `deliverables/screen_phonons/`. What it says: (1) the softest transverse
rigid-chain mode of the polar cell is 34.5 cm^-1 for beta-PVDF (construction check against
`docs/PHONONS.md`'s 34.0), 22 for CFE, 16 for TrFE-cand at its gamma-free cell, 46 for CDFE,
53-79 for the bulky-pendant candidates; VDCN's polar cell has 13 and 15 cm^-1 nitrile pendant
librations, the softest in the set and not rigid-chain. (2) The nitrile chemistries' rigid
all-trans repeats are 0.3-0.6 A from the all-atom minimum and their "antipolar" cells become
polar (0.24-0.27 e.A per monomer) once the pendants relax: the screen's rigid-chain polarity
column does not survive for AN, VDCN, VClCN-cand. (3) With gamma free, PVDF's +1.021 and
CFE's +1.501 polarity gaps collapse to -0.000 and +0.071 (DESIGN.md section 6 recorded this
for beta only); of the screen's three resolved-polar verdicts only CDFE's survives; after
all-atom relaxation PVDF, CFE and TrFE are degenerate to 0.02 kcal/mol per monomer. The
screen's polarity column should be re-measured with gamma free on both branches (not done;
`antipolar_cell_exact` screens at 90 only). (4) PVDF's gamma-free antipolar competitor
(a rectangular cell with chain 2 above chain 1 along the long axis) is a true minimum, stiffer
than beta (41/49/61 vs 34/40/49 cm^-1) and degenerate with it: beta is not near a soft-mode
instability on this potential, and polar-vs-antipolar for PVDF is a barrier question this
does not answer. Dipole reviewed the first record and found 12 of 36 relaxed cells a few 1e-6
above the requested force tolerance because `relax_all_atom` discarded termination; fixed at
the owner (`phonon.Relaxation`: unrounded force, termination evidence, Newton polish,
`converged`), consumers require acceptance; record regenerated 2026-09-29 (36/36 accepted, 13 via the Newton polish, every tabulated number unchanged to its printed digits). `pack.polish` now wraps `dz` (and the setting
angles) inside its objective, not only afterwards (2026-09-29): the packer's energy is only
periodic in `dz` for about two repeats outside `[0, c)`. The merge hold below is
about merging *their* branch into ours; Dipole's 2026-09-14 notes say their C: checkout is no
longer touched by incoming code (their relay is on D:), so our own source edits on this branch
are safe to push.

## Dipole's Polyfind extension branch (do not merge yet)

`origin/physics-vector-repeat` (e6150ed, 2026-09-13) is Dipole's isolated
extension of our source: a complete placed-cell energy with real Cartesian
torsion (`CrystalPacker.placed_energy_and_grad(P, H, flip)`), general
triclinic images, vector-repeat valence/flux/exclusion gradients, and a
six-column affine clamped-ion record identical to ours on the four shared
columns. Verified here: 640 pass / 5 skip / 0 fail. **Do not merge it into
the shared branch until Dipole releases the dependency**: their C: checkout
tracks our branch and their running native jobs depend on unchanged source.
Read it in a detached worktree. Torsion changes nothing at Gamma for beta
(primitive two-site repeat cannot twist); it raises alpha's torsional modes
by 4-11 cm^-1. Next lattice-dynamics step: beta with a doubled repeat.

## Dipole's open gates (theirs, not ours to chase)

- Berry polarization: **matrix running, three diagonal columns compared.**
  Their reducer is Vanderbilt-proper (`e = dP/de + d_jk P_i - d_ij P_k`); since
  2026-09-29 so is `mechanics.piezoelectric`'s `e` (the numbers below already
  applied the correction in `clamped_ion_columns.py`). Clamped-ion e_y,jj, their frame, both
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
- **Closed 2026-09-29, ours:** the `e` definition. Derived
  (`docs/ELECTROMECHANICS.md` 5.10): an electrode integrates the nominal
  polarization `F^-1 mu / V0` and a voltage fixes the nominal field `F^T E`,
  so the measured coefficient is Vanderbilt's. `mechanics.piezoelectric` now
  returns it as `e` (old ones kept as `e_reference_volume`, `e_improper`), and
  the converse route holds the nominal field, adding the field-stretch stress
  term both routes had omitted. Born-fitted preset: d33 -8.80 -> -9.86,
  d31 +0.02 -> +0.42, routes agree to 0.02 pC/N on both. Not the shortfall's
  cause (3% and 2% of the gaps); the proper dimensional term gives a
  *positive* d31 (-P S_13), so section 5.3's sign argument is withdrawn.
  Examples now label the Vanderbilt d as "film"; recorded deliverables that
  quote the old film numbers (screen, fit_born_flux runs) are stale until rerun.
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

- Tests: `export PYTHONPATH="$PWD/src"; python -m pytest tests -q -p no:cacheprovider`, 724 pass, 5 skip (GPU). ~5-7 min on an idle box. Dipole's runtime reports 3 last-digit
  frozen-literal failures in `test_mechanics.py` (PE/alpha/gamma energies at
  1e-15 relative) that do not reproduce here at any commit; machine noise,
  not a defect. Since 2026-10-08 every recorded golden value is held to 1e-12
  relative and the digest tests compare stored arrays (`tests/data/pendant_builds.npz`);
  `docs/PERFORMANCE_REVIEW.md` 15.6 has the per-test record.
- **Never set `POLYFIND_TABLE_CACHE=1`**: it is read as a directory name and once committed 85 MB of caches. It is gitignored now.
- Agents: launch in worktrees; they often pause on their own background jobs and report "waiting"; that is not stuck. Merge, run the suite, push, then remove the worktree.
- Commit attribution: `Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>` and the session URL, per the current system reminder.

## If continuing

On our side: the lattice-curvature column exists now (`docs/SCREEN_PHONONS.md`); its next
steps are the gamma-free re-measurement of the screen's polarity column (the alpha/delta
pair, 0.030 kcal/mol apart, is the same unresolved question), and Dipole's bulk periodic
curvature when it lands; the alpha far-IR comparison is done (above). On the piezoelectric side the xx column says the answer is "both":
a real electronic clamped-ion term of ~0.13 C/m^2 and an overstated target.
When the yy column and the accepted matrix land: write the scope boundary
into `DESIGN.md` and `docs/BENCHMARK.md` with the measured size of the
electronic term, and restate the piezoelectric shortfall against the accepted
Berry-derived total rather than the failed sweep's 0.4-0.6. Either way, reply in `docs/REFERENCE_DATA_REQUEST.md`, then re-run
`examples/electromechanics.py --preset flux` if any parameter changes.
