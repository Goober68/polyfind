# Valence terms, real forces, and where the wall is now

`docs/DFT_FIT.md` fitted the rigid potential — one Fourier torsion, 12–6 Lennard-Jones,
fixed atom-centred point charges — to 467 frames of PBE-D3 conformers and torsion scans.
It generalised (held-out energy error 3.32 → 1.76 kcal/mol over ten unseen chemistries)
and it fixed **one** of three acceptance tests. DESIGN.md 5.9 then said, plainly, that the
binding constraint had moved from the data to the functional form, and named two gaps:

* no bond or angle terms, so **93% of the reference force signal was unrepresentable** and
  a strained all-trans chain could not relieve a contact by opening an angle — the reason
  PVDC's planar zigzag carried 313 kcal/mol (74 now that its backbone angles are the
  measured 123/114 rather than an equal 114 — section 4) and FANOME's about a million;
* fixed atom-centred point charges, which the fit tried to escape by driving **fluorine's
  partial charge positive** once its bounds were released.

This document is what happened when both were closed. The code is `polyfind.fitting`
(`ValenceDesign`, `fit_valence`) and `polyfind.forcefield` (the `bond_terms`,
`angle_terms` and `charge_offsets` parameters of `SimpleFF`); the script is
`examples/fit_valence.py`; the result ships as the opt-in preset `pvdf-dft-valence`.
`SimpleFF()`'s defaults are bit-for-bit untouched and still the default.

**The headline, first, including the part that did not work.** With bond and angle terms
the forces become representable — held-out force error falls from 13.9 to 5.2 kcal/(mol Å)
against a reference whose own RMS is 14.0, and the correlation with the reference from
+0.10 to +0.93 — and the held-out energy error falls from 1.76 to 1.36 kcal/mol. Fluorine's
charge comes out **+0.114 e and off its bound**. The all-trans strain problem is largely
gone: with the angles free, PVDC falls from 47 to 12 kcal/mol and PVDC's own relaxed CH2
angle lands at 123.1°, the measured value.  (That first figure was 313 until PVDC's
backbone angles were themselves corrected to the measured 123/114 — see section 4.) But the forces **still do not help the
energies** (1.16 without them against 1.36 with), the off-site charge is within noise of
doing nothing, and of the three acceptance tests the score is unchanged at **one pass, two
failures**.

---

## 1. What was added

### Bond and angle terms

`SimpleFF` gains `bond_terms` and `angle_terms`: harmonic stretching `0.5 k (r - r0)^2` and
bending `0.5 k (theta - theta0)^2`, typed by the bond graph rather than by the polymer
(`bond_type_name`, `angle_type_name`), so one table serves a built `Structure`, a periodic
block and an inferred `Frame` alike. Both default to `None`, which is the rigid potential
every earlier number in this package was computed with.

Two typing decisions are measurements rather than conventions:

* a bond to a **one-coordinate N or O** is typed apart (`C=N`, not `C-N`): typed together,
  the reference set's "C–N" population spans 1.16 to 1.52 Å and no single harmonic
  describes it; split, the nitrile's 171 bonds sit at 1.140 ± 0.002 Å;
* an angle about a **two-coordinate centre** gets no term at all. Those are linear groups
  (a nitrile carbon) whose reference angles sit at 179–180°, where the derivative of
  `arccos` is singular and a harmonic in the angle is the wrong form. `polymers.py` already
  treats the nitrile as rigid; this is the same choice in a second place, not a new
  approximation. The largest angle that does carry a term is 128°.

Types, chosen by coverage of the reference set: stretch `C-C` (4561 bonds / 31 systems),
`C-H` (5081 / 31), `C-F` (5328 / 31), `C-Cl` (80 / 4), `C-O` (80 / 3), `C=N` (171 / 10) and
a wildcard for the remaining 235; bend `C-C-C`, `C-C-H`, `C-C-F`, `H-C-H`, `F-C-F`, `F-C-H`
— 95.7% of the 29,954 angles — and a wildcard. 28 new parameters, every equilibrium value
and every stiffness fitted.

**Inside a packing they change nothing, by construction.** `pack` builds rigid chains, so
every bond length and every bond angle is fixed and the valence energy is one additive
constant per chain: it cancels from every lattice-energy difference, every polarization and
every cell parameter. `FFParameters.applied` therefore does not forward them, and
`tests/test_forcefield.py` asserts the constancy over conformers rather than arguing it.

### Off-site charge for fluorine

`charge_offsets` slides an element's point charge along its own bond while its
Lennard-Jones centre stays on the nucleus. One parameter, `d_F`, and one sentence of
justification: **an atom-centred model forces one point to serve both fluorine's van der
Waals size and the centroid of the C–F dipole, and the fit had been resolving that conflict
by switching the dipole off.**

`pack` and `refine` carry one charge per atom and no off-site machinery, and this work was
not allowed to touch them, so a fitted offset reaches a crystal through
`project_offset_charges`: moving `q` a distance `d` along a bond of length `r` adds `q d` to
the dipole, and putting `q d / r` more charge on the atom and the same amount less on its
neighbour adds exactly the same thing. Same total charge, same total dipole, different
quadrupole — and at the 3 Å-and-longer separations a lattice sum is made of, the dipole is
the term that matters. The test asserts both invariants.

### Full Cartesian forces in the objective

`ValenceDesign` replaces the torsional projection with every force component of every atom:
32,277 numbers over the training set against 312 energies. Everything geometric is still
precomputed once — pair distances, bond and angle values, the analytic `dtheta/dx` and
`dphi/dx`, and, for the offset, the three constants that make the Coulomb site separation
`sqrt(c0 + c1 d + c2 d^2)` and the 3×3 projector that maps a site force back onto its two
atoms — so one evaluation is 43 ms and the 56-parameter fit is 260 s. `evaluate` reproduces
`SimpleFF.energy_frame` to 1e-13 and `SimpleFF.forces_frame` (plain central differences) to
1.5e-7 kcal/(mol Å); the tests assert it, because a design that drifted from the shipped
potential would make every number here meaningless.

## 2. The fit

56 parameters: the previous 27, unchanged and in the same order, plus 14 stretch, 14 bend
and the fluorine offset. Same split (21 training chemistries / 312 frames, 10 held out /
155), same weak ridge (0.3), same optimiser. Four starts, all four reaching cost 0.920769 —
the previous fit's four starts did not agree that closely. The valence block starts at
textbook GAFF-like values, deliberately *not* at averages of the reference geometries,
which would make the ridge a second and hidden fit to the same data.

## 3. Results

### Held out by chemistry, scored by one design

Energy RMS in kcal/mol after removing each system's own mean; force RMS in kcal/(mol Å);
`r` is the correlation between model and reference force components.

| model | train E | held E | train F | held F | r(F) |
|---|---|---|---|---|---|
| illustrative (rigid) | 3.156 | 3.316 | 13.93 | 14.37 | +0.06 |
| `pvdf-dft-fit` (rigid) | 1.498 | 1.757 | 13.70 | 14.06 | +0.10 |
| the rigid form refitted to these Cartesian forces | 1.580 | 1.995 | 13.48 | 13.86 | +0.17 |
| **this fit** | **1.022** | **1.359** | **4.64** | **5.23** | **+0.93** |
| this fit, force weight 0 | 0.859 | 1.157 | 41.98 | 44.29 | +0.77 |

The reference forces themselves have RMS 14.05 kcal/(mol Å) held out, which is what a model
predicting *zero* would score. Per held-out system: `vdcn` 0.50, `fan` 0.55, `cfe_ter` 0.61,
`tfe` 0.69, `fancl` 0.72, `pp` 0.91, `an` 0.93, `itfe` 0.96, `hfp` 1.04, `mvs` 3.32. PVDF
itself (in training) is 0.47, against 1.13 before. The one bad system is the sulfone again;
excluding it the held-out error is **0.80** kcal/mol over 137 frames, against 1.21 before.

### How much do the forces now contribute?

They contribute the forces, and nothing else. This is the clearest result here and the most
uncomfortable one.

| force weight | 0 | 0.003 | 0.01 | **0.03** | 0.1 | 0.3 | 1.0 |
|---|---|---|---|---|---|---|---|
| held-out energy (kcal/mol) | **1.157** | 1.353 | 1.360 | 1.359 | 1.390 | 1.538 | 1.940 |
| held-out force (kcal/(mol Å)) | 44.29 | 6.44 | 5.56 | 5.23 | 4.98 | 4.79 | 4.63 |
| correlation with reference forces | 0.768 | 0.910 | 0.927 | 0.934 | 0.940 | 0.944 | 0.947 |
| `q_C-F` (e) | 0.113 | 0.113 | 0.114 | 0.114 | 0.110 | 0.104 | 0.099 |
| `d_F` (Å) | 0.100 | 0.152 | 0.161 | 0.177 | 0.210 | 0.249 | 0.324 |

Three things to read off it.

**The valence terms are what made the forces usable, and that is now a measurement rather
than an inference.** Same objective, same data, same optimiser, only the functional form
differing: the rigid form fitted to full Cartesian forces still scores 13.86 against the
reference's own 14.05, i.e. no better than predicting nothing at all, with a correlation of
+0.17. With valence terms it scores 5.23 at +0.93. That is `DFT_FIT.md` section 3's "93% of
the force signal is bond stretching" cashed out.

**The forces still do not help the energies.** Held-out energy is *better* with the force
block switched off, 1.157 against 1.359, and worsens monotonically from 0.003 upward. The
diagnosis is the one the previous fit already reached and this fit now measures directly:
the reference geometries were relaxed at some other level of theory and labelled with
PBE-D3, so their forces record a systematic bond-length disagreement rather than a
conformational gradient. The fitted equilibrium lengths say so out loud — `r0(C-F)` comes
out 1.467 Å against a measured mean of 1.358, and `r0(C-C)` 1.458 against 1.535, i.e. the
forces want every C–F longer and every C–C shorter than the geometries are — and
reproducing that pulls the same parameters the relative energies want elsewhere. Section 2
of `DFT_FIT.md`'s "what would be next" asked for reference geometries relaxed at the
labelling level; this is the same request, now with a price attached: 0.2 kcal/mol of
held-out energy error.

**Weight 0.03 is shipped anyway**, for a reason the energy table does not show. At weight 0
the *conditioning* of the fit is bad: the Jacobian of the energy block alone, prior-scaled
over all 56 parameters, has a condition number of **19,800**, and adding the force block
drops it to **478** — comparable to the previous, 27-parameter fit's 504. Energies alone do
not determine 56 parameters; energies plus forces do. And a potential whose gradients are
worse than zero is not usable by anything that takes one, which `refine_crystal` does.

### Fitted values worth naming

| | fitted | note |
|---|---|---|
| `q_C-F` | **+0.114 e** | was 0.020, on its bound; now off it |
| `q_C-H` | −0.119 e | H positive, as electronegativity requires |
| `d_F` | +0.177 Å | the charge sits *beyond* the nucleus: a longer C–F dipole than the nuclei suggest |
| `scale14` | 0.923 | |
| V1,V2,V3 [CF2-CH2] | −0.459, −0.482, +1.724 | PVDF's only torsion type |
| `kb`, `r0` C–C | 100 ← bound, 1.458 Å | the stiffness ran to its floor |
| `kb`, `r0` C–H | 545, 1.100 Å | textbook-ish (the earlier diagnostic fitted 475, 1.104) |
| `kb`, `r0` C–F | 260, 1.467 Å | |
| `ka`, `theta0` C–C–C | 77, 108.8° | |
| `ka`, `theta0` H–C–H | 85, 110.9° | |
| `ka`, `theta0` F–C–F | 110, 109.9° | |

Six of 56 parameters sit on a bound: `x_H`, the three well depths `D_C`, `D_H`, `D_F` (all
the same four the previous fit had) and the C–C and wildcard stretch stiffnesses at their
floor of 100 kcal/(mol Å²). The last two are the data asking for a C–C bond four to six
times softer than any real one, which is worth stating plainly: the fit would rather the
C–C stretch did not contribute to relative conformer energies at all. Against that, **no
electrostatic parameter is on a bound any more**, which is the one that mattered.

### Does the fit still want a positive fluorine?

No, and it does not need the off-site charge to stop wanting one — which makes the
electrostatics change a partial success at best.

* At every force weight, and with `d_F` held at zero, `q_C-F` lands at +0.113 to +0.115 e.
  With `d_F` free the held-out error is 1.359 and `q_C-F` = +0.114; with `d_F` fixed at
  zero it is 1.347 and +0.115. The off-site site is, on this evidence, **a null result**:
  it costs nothing and buys nothing measurable.
* The decisive test is the one `DFT_FIT.md` used to expose the problem: take the bounds off
  and see what the fit chooses. Widening the Lennard-Jones box to 0.6–1.5 on the radii and
  0.05–4.0 on the well depths, and dropping the electronegativity signs entirely so every
  increment is free in ±0.6 e:

  | unbounded fit | `q_C-H` | `q_C-F` | `q_C-Cl` | `q_C-N` | held-out E |
  |---|---|---|---|---|---|
  | rigid form, energies only | **+0.264** | +0.177 | **−0.023** | **−0.019** | 1.432 |
  | rigid form, energies + forces | **+0.229** | +0.246 | +0.027 | +0.008 | 1.465 |
  | **valence + off-site** | **−0.098** | **+0.122** | **+0.102** | **+0.359** | **1.171** |

  Bold marks a sign that contradicts electronegativity. Unbounded, the rigid form gets
  three of six increments backwards — here it inverts *hydrogen* rather than fluorine, and
  which one flips depends on the start; that one flips does not. The new form, given the
  same freedom, chooses the electronegativity-correct sign for **every** increment without
  being told, and fits better while doing it.

So the honest statement is that **the valence terms, not the richer electrostatics, are
what stopped the fit fighting the charges.** With bonds and angles carrying the strain and
the forces, the point charges are no longer the only flexible thing in the model, and they
settle where chemistry says they should. The off-site fluorine site is along for the ride:
it is the right kind of change, it is measurably not degenerate with the charge, and it did
not need to be made. (That the unbounded fit also generalises better, 1.171 against 1.359,
says the shipped bounds now cost accuracy rather than buying physicality — an argument for
loosening them next time, made here rather than acted on.)

### Identifiability

Measured, not assumed (`valence_degeneracies`, and the tests):

* **charge scale against permittivity: still exactly degenerate.** Moving a charge site does
  not touch the argument — every energy still depends on the charges only through `q_i q_j`
  and on `eps_r` only as a divisor — so scaling every increment by 1.7 and the Coulomb term
  by its square moves no frame energy by more than 3e-13 kcal/mol. `eps_r` stays at 1. The
  new electrostatics does **not** lift this, and nothing but a polarization will.
* **charge against offset: correlated, not degenerate.** Only the product `q d` survives at
  long range, so this was the obvious new risk. Holding the C–F bond dipole fixed while
  trading one against the other moves frame energies by up to 2.79 kcal/mol at the optimum
  (9.44 at the starting point); their correlation in the fitted covariance is −0.54. The
  short range separates them, which is exactly the part of the model the offset was added
  for.
* **the worst-determined pair in the whole vector** is `kb(C=N)` against `r0(C=N)` at
  −0.999: 171 nitrile bonds that all sit within 0.002 Å of each other cannot separate a
  stiffness from an equilibrium length. Harmless — it is a nuisance type — but it is what
  the condition number of the energy-only block is mostly measuring.

### Acceptance tests

None of these is in the objective. Compared against `pvdf-dft-fit`, re-measured here rather
than quoted:

| | `pvdf-dft-fit` | this fit | target | |
|---|---|---|---|---|
| 1. E(alpha) − E(beta), kJ/mol per monomer | −4.81 | **−3.13** | −6.5 to −2.6 | **PASS** |
| 2. isolated-chain RIS, top candidate | TG+ | TG+ | anything but TG+ | **FAIL** |
| 3. alpha E(antipolar) − E(polar), kcal/mol per monomer | +0.092 | **+0.084** | ≤ 0 | **FAIL** |
| beta \|P\|, C/m² | 0.115 | 0.115 | 0.176–0.188 (DFT) | |

Re-measured after the batched geometry corrections of `REFERENCES.md` (PVDF's C–C bond
1.54 → 1.528 Å); the same rows read −4.54 / −2.96, +0.091 / +0.085 and 0.114 / 0.114
before them. The fit itself was run at the old geometry and has not been repeated.

**Test 1 still passes**, at −3.13 kJ/mol per monomer against the previous fit's −4.81. It is
inside the range four independent studies agree on, but nearer its edge, so this is a small
step in the wrong direction inside a passing test. The geometry correction helped it
slightly (−2.96 → −3.13, further from the −2.6 boundary).

**Test 2 still fails, and the geometry correction pushed it back the wrong way**: TG+ at
−4.37 kcal/mol per monomer against TG+G+G+ at −3.94, a gap of 0.43, where at the old
geometry it was 0.18 against TTG+G+. Measured on the same footing the illustrative
potential's gap is 1.58 and `pvdf-dft-fit`'s 0.49, so the fits still move this the right
way and still do not arrive — but the run of three successive narrowings recorded here is
not what a fresh measurement shows, and the 0.18 should not be quoted as the current
number.

**Test 3 still fails and barely moved**: +0.084 against +0.092. The electrostatics change
was aimed squarely at this test and did not shift it. What did change is that the reason
can no longer be "the fit refuses to have any C–F electrostatics": `q_C-F` is now +0.114
and beta's polarization is unchanged at 0.115 C/m², still well below the DFT range of
0.176–0.188. The remaining suspects are the ones the fit cannot see: the objective still
contains no polarization, and `applied` reaches the lattice with the dipole-equivalent
projection of the off-site model rather than the model itself.

**One finding came out of running test 2 rather than from the fit.** The new fit's adapted
gauche state sits at 60° where the previous two potentials adapt to 80°, and at 60° the
(G+, G−) pair is a hard steric clash under rigid geometry — about 1e6 kcal/mol. The
third-order RIS term is built by inclusion–exclusion at the state angles while the pair term
is measured at a relieved *basin minimum*, so subtracting one clash from another left a
spurious well that the ±50 clip turned into a flat −50 kcal/mol **bonus** for an impossible
conformation, and the ranking put `TG+G+G-` on top at −26.6 kcal/mol per monomer. That would
have been a "pass" for test 2 and it was an artifact. `_fit_third_order` now refuses to
apply the inclusion–exclusion across a steric overlap: a triple that overlaps scores `+cap`,
and a triple whose *subtracted* terms overlap scores 0. Where nothing overlaps — every
triple of the illustrative and `pvdf-dft-fit` PVDF fits — the arithmetic is unchanged, and
both of their rankings are bit-for-bit what they were. The general point is the one this
whole document is about: the frozen backbone angles that made the all-trans reference
unphysical also make the third-order RIS fit fragile, and the fix there was a guard, not a
cure.

## 4. The strain finding, re-measured

Lennard-Jones energy of a ten-bond all-trans oligomer, the same measurement
`polyfind.polymers` and `tests/test_pvdc.py` record, now with the backbone angles free to
relax under the fitted potential (`relax_backbone_angles`; one variable per backbone atom
of the repeat, the substituents following the changed backbone):

| polymer | LJ, rigid, unfitted | LJ, rigid, this fit | **LJ, angles relaxed** | E(relaxed) − E(rigid) | relaxed angles (°) |
|---|---|---|---|---|---|
| PVDF | 11 | 10 | **6** | −6 | 114.0, 107.3 |
| PVDC | 74 | 47 | **12** | −33 | **123.1**, 102.1 |
| CFE | 162 | 120 | **12** | −96 | 120.8, 103.2 |
| CDFE | 199 | 154 | **50** | −89 | 113.7, 124.2 |
| AN | 88 | 27 | **2** | −32 | 116.1, 102.9 |
| VDCN | 176 | 57 | **−1** | −65 | 118.6, 102.3 |
| FANOME | 1.0e6 | 1.4e5 | **4.8e3** | −1.4e5 | 135.0, 135.0 ← both at the bound |

Re-measured after the batched geometry corrections of `docs/REFERENCES.md`. The
PVDF row moved because its C–C bond went from 1.54 to 1.528 Å (9 / 8 / 5 / −5 /
113.7, 107.3 before), and the PVDC row because its backbone angles went from an
equal 114/114 to the measured 123/114 (313 / 231 / 12 / −199 / 123.1, 102.1
before). **The relaxed columns are unchanged for PVDC to three figures**, which is
the check that these are starting-point effects and not different minima: the
relaxation lands on the same 123.1 / 102.1° and the same 12 kcal/mol from either
start. What changed is how much of the drop the angle terms can claim. Most of
PVDC's celebrated 231 → 12 was the model starting from the wrong angle; the honest
figure for what a bend term buys is 47 → 12.

**This resolves the reference-state problem for five of the seven chemistries, and it is
worth saying plainly.** PVDF, PVDC, CFE, AN and VDCN all come into the same range as PVDF's
own strain once the angles can open — 12 kcal/mol or less on a ten-bond oligomer — so the
RIS convention of measuring energies from all-trans is measuring from a state those chains
could actually occupy. `docs/CHEMISTRY_EXTENSION.md`'s phase-3 follow-up concluded that
"variable backbone angles are no longer a refinement nicety; for anything past PVDF they
are a precondition for the model to describe a real molecule at all". That conclusion is
confirmed, and the precondition is now met for the chemistries that had the problem.

Two exceptions, stated rather than smoothed over. **CDFE** still carries 50 kcal/mol: its
relaxation opens the CF2 angle to 124° instead of relieving the Cl contact, so a
two-angle-per-repeat relaxation is not enough there. **FANOME** is still not a physical
structure — 4850 kcal/mol with both angles pinned at the 135° bound — which is the same
verdict as before, and it remains registered and deliberately not fitted.

**PVDC's angles against the measurement.** The published structure (Takahagi et al., 1988)
has C–CH2–C = 123° and C–CCl2–C = 114°. Relaxed under this potential PVDC gives
**123.1° and 102.1°**: the wide angle is right to a tenth of a degree, and the narrow one is
12° too closed. So the mechanism the audit proposed — that PVDC relieves its Cl···Cl
crowding by opening the CH2 angle — is reproduced quantitatively, while the compensating
angle is overdone. The relaxation is a two-parameter minimisation of one potential, not a
structure refinement, so the agreement on the first number is the meaningful part and the
second is a limit on how far it should be pushed.

## 5. What would be next

Unchanged from `DFT_FIT.md` in its first item, sharpened in the rest.

1. **A polarization or a molecular dipole in the training set.** Still the only observable
   that separates `eps_r` from the charge scale, still what test 3 turns on, and now the
   only remaining explanation for test 3 that has not been eliminated.
2. **Reference geometries relaxed at the labelling level.** The forces would then be
   conformational gradients rather than a record of a 0.1 Å bond-length disagreement, and
   the 32,000 force components would improve the energies instead of costing 0.2 kcal/mol
   to fit alongside them.
3. **Backbone angles as search variables, not only as a relaxation.** Section 4 shows the
   angles that matter move by 9 to 20 degrees from the frozen values. `refine_crystal`
   already relaxes them inside a packing; `fit_ris`, `enumerate_periodic` and the RIS model
   itself still assume they do not move, which is what made the third-order term fragile
   above.
4. **More than two angle variables per repeat for the bulky chemistries.** CDFE and FANOME
   are the cases where one angle per backbone atom of the repeat is demonstrably not
   enough.
