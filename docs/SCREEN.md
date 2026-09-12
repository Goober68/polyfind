# The first real screen: can this package rank electroactive polymers?

The campaign wants polymers that convert field energy into mechanical work: PVDF plus CFE,
CDFE, AN, VDCN, VDC, CNEPO and FANOME. This is the first time the whole funnel has been
pointed at that question rather than at PVDF alone.

Everything here is reproduced by `examples/screen_electroactive.py`, except the acceptance-test
table, which comes from `fitting.acceptance_tests`; the commands for both are at the end.

Potential: **`pvdf-dft-valence`** everywhere (RIS fit, packing, refinement, deformable mechanics)
with **`pvdf-dft-valence-flux`** added to the mechanics packer only, the same split
`examples/electromechanics.py` uses. It is the best available fitted preset: the only one with
valence terms, hence the only one on which `C_33` and the diagonal columns of `e` exist at all,
and the lowest held-out energy error of the three fits (1.36 kcal/mol against 1.76 and 3.32).
It is fitted to **PVDF** single-chain PBE-D3 energies and forces, so every number for a
chlorine or nitrile chemistry is a transfer with nothing behind it but the wildcard bond and
angle types of `fitting.VAL_BOND_TYPES`.

One column of this document — polar versus antipolar — has been wrong four times over, for four
different reasons, two of them found on this re-run. The history is kept, in one place, in "What
was wrong" below, because it is the reason that column is now reported the way it is.
**Everything before that section is the current state; every retraction is collected inside it.**

## Answer first

**The model cannot decide polarity for the phase that matters, and that is the screen's primary
output.** For the all-trans planar zigzag — the only conformation that can carry a
ferroelectric polarization, and the one the campaign is asking about — the polarity gap is
**inside the fitted potential's own error bar for five of the nine chemistries that have one**.
Of the four it resolves, three say polar (PVDF, CFE, CDFE) and one says antipolar (VDCN). So
the screen's honest answer on its headline column is "cannot tell" for the majority of the
sample, and the useful conclusion is not a table of numbers: **the potential has to improve
before polarity screening means anything.** What that needs is stated at the end.

**Yes, and decisively, for the structural half that is not polarity.** Whether a polar chain
conformation is *accessible* — how far the planar zigzag sits above the chain's own
conformational ground state, and above its own lattice ground state — separates the chemistries
by 0.4 to 21 kcal/mol per monomer against a 0.27 error bar. That column did need one correction
after all: VDCN's conformational row ("all-trans ground state, rank 1") was a frozen-angle
artefact of the RIS scan (`docs/NITRILE_LANDSCAPE.md`) and is regenerated below with the
backbone angles relaxed (all-trans 43rd of 67, +4.1), as are AN's, PVDC's, CFE's and CDFE's;
PVDF's stands bit-for-bit. The lattice half is an order of magnitude clear of the noise where
it matters, and it is what this screen is actually good for. It rejects PVDC (no polar
conformation closes at all), rejects the nitriles (their polar zigzag sits 4.1 and 16.3 kcal/mol
per monomer above their own lattice ground state), and picks out trifluoroethylene at the
lattice level, which is the known answer — though no longer at the conformational level, where
the relaxed scan puts its zigzag behind PVDF's.

**Yes, with a caveat, for the response ranking.** The ordering of the dipole-strain response
survives leave-one-out at Spearman `rho = +1.000` and the scale error is systematic to within a
factor of 1.7 — but the ranking is carried by the *fixed* charges, not by the fitted flux, so
what the leave-one-out licenses is coarser than it looks. See step 1.

**No for the response figures as a ranking of usefulness.** Work density is *positively*
correlated with how far the polar phase sits above its own ground state (+0.65), so ranking on it
alone picks the chemistries whose polar phase the polymer would never adopt. The other half of
that story as published — work density *anti*-correlated with the transverse chain dipole at
−0.69 — does not survive re-measurement: it is −0.16 on the corrected sample. The accessibility
column has to gate the response column; the dipole column does not enter.

**Of the three design rules this screen published, one stands.** Rule 2 (ranking on work
density alone inverts the answer) survives, because it rests on the columns that resolve. Rule 1
(a perpendicular dipole above 0.5 e.A per monomer forces antipolar packing) falls from four of
four to one of four, with the other three not resolved at all. Rule 3 (denser packing helps) is
demoted to an observation about a sample of eight. **The `trfe-cand` recommendation is withdrawn
as it was stated**; a weaker claim about a different column survives in its place.

## How to read any energy in this document: the error bar

Every energy here is a per-monomer lattice energy or a difference of two, and the question "is
this difference a prediction?" has one answer for all of them.

`pvdf-dft-valence`'s **held-out energy error is 1.359 kcal/mol RMS** over ten chemistries the
objective never saw (`docs/VALENCE_FIT.md`; re-measured from the shipped preset for this
document — 1.3593, max 5.75). That is a per-*frame* number and the frames are ten-bond
oligomers — eleven backbone carbons, five monomers — so the comparable quantity for a
per-monomer crystal energy is

    1.359 / 5 = 0.27 kcal/mol per monomer.

**A difference smaller than 0.27 kcal/mol per monomer is not a prediction of this model.** That
is `screen_electroactive.ERROR_BAR`, the script prints it beside every gap, and a gap inside it
is reported as `not resolved` rather than as a verdict. Three things to say before any table:

* It is generous to the model in one direction and harsh in the other. Generous, because
  dividing a total-energy error by the oligomer length assumes the error is extensive and
  uncorrelated, which is the most favourable assumption available. Harsh, because nothing
  guarantees that a *difference* between two cells of the same chain inherits the full error of
  an absolute energy — systematic error cancels in a difference. The honest reading is that 0.27
  is the right order of magnitude, and that the conservative reading, 1.36 per monomer, resolves
  **nothing at all** in this document.
* It is not the search error. The antipolar subspace search is converged to about 1e-3 kcal/mol
  per monomer, three orders of magnitude inside the error bar, and the search is not what limits
  anything here. Two defensible *protocols*, though, differ by more: the screen's α-PVDF gap is
  +0.027 and the acceptance test's is +0.084, the same quantity on the same chain, differing only
  in whether the polar branch's screen table was built inside or outside the
  `FFParameters.applied()` block. 0.06 is a fifth of the error bar and twice the effect of
  correct electrostatics, which is itself the point.
* It is not so large that it swallows every useful answer. The α/β polymorph gap, 0.747 kcal/mol
  per monomer, is comfortably outside it and lands inside the range four independent studies
  agree on. The error bar swallows the polarity column specifically.

## Step 1. Is the response ranking meaningful?

`docs/BENCHMARK.md` records the piezoelectric magnitudes as 5x and 9x short, from a charge-flux
coefficient fitted to twelve numbers at `R^2 = 0.39`. A screen survives that only if the error is
**systematic**, and the sibling project's `field_neighborhood_refined` is the one place that is
testable: four chemistries with an axial dipole response at two strains. The shipped coefficient
was fitted to all four, so comparing against them is circular. Refit with each chemistry held
out, predict the held-out one, compare orderings.

Reference `dmu/d(axial strain)`, e.A per unit strain (exploratory GFN2-xTB, finite pinned
two-chain pair in vacuum — **not** a bulk crystal):

| system | reference vector | \|ref\| | LOO prediction, shipped fit | pred/ref |
|---|---|---|---|---|
| AN | (+13.195, −18.645, −24.372) | 33.40 | 24.25 | 0.73 |
| VDCN | (+9.472, +10.804, +12.159) | 18.82 | 17.92 | 0.95 |
| CNEPO | (+3.485, −0.290, −4.282) | 5.53 | 3.26 | 0.59 |
| PVDF | (−0.106, +1.909, −4.143) | 4.56 | 2.49 | 0.55 |

**Verdict: the ranking survives.** `AN > VDCN > CNEPO > PVDF` is reproduced exactly under
leave-one-out, `rho = +1.000`, component-wise Pearson `r = +0.966`. It survives for two
alternative fit families too (`rho = +1.000` for angle-only-four-pairs and for angle+bond), even
though the held-out coefficients swing wildly — `k_angle(C-H)` runs over −1.62 .. −0.15 across
the four holdouts, a factor of eleven, and the ordering does not notice.

**The error is systematic, and here is its size.** Predicted/reference spans 0.55 .. 0.95: the
model under-predicts *every* chemistry, by between 5 % and 45 %. Take out that common
under-prediction and what is left is a factor of **1.74** between the best- and worst-predicted
chemistry. **Two chemistries whose true responses differ by less than about 1.7x can therefore
change places, and anything further apart will not.** The reference set spans 7.3x from PVDF to
AN, which is why its ordering is never in doubt and also why that ordering on its own is a weak
test. This spread, not the correlation, is what the screen can resolve.

**And the correlation is weaker evidence than it looks, for two separate reasons.**

*First, n = 4.* A random ordering of four systems gives `rho = 1` with probability 0.042 and
`rho >= 0.8` with probability 0.167. One perfect ordering is a one-in-twenty-four result, not a
demonstration.

*Second, and more important: the flux is not what does the work.* Run the same test with the
**fixed-increment** charges and no flux at all — a charge model with nothing whatever fitted to
this data — and it also gives `rho = +1.000`, with a *better* component-wise correlation
(`r = +0.989` against +0.966) and a *tighter* scale error (0.54 .. 0.78, a factor of 1.44). The
ranking is a property of the chemistries' geometry and of bond-charge increments fitted to
PBE-D3 energies; the flux changes magnitudes by up to 22 % and never changes the order.

That relocates the trust. The badly determined `R^2 = 0.39` coefficient is what makes β-PVDF's
`d_31` non-zero and positive at all, and that remains as poorly determined as
`docs/BENCHMARK.md` says. It is **not** what ranks chemistries. So the step-2 response columns
are usable as a coarse ranking and remain unusable as values.

**One limitation no amount of resampling fixes.** It is run on the *reference's own geometries*:
it answers "does a charge model calibrated without chemistry X reproduce X's dipole response at
a geometry someone else supplied?" It says nothing about whether polyfind's own predicted
structure for X is the right one, and the step-2 response columns depend on that too. So the
leave-one-out validates the charge model's transferability and **not** the pipeline end to end.
The structural half has its own independent validation — the PVDF polymorph ordering, cells and
densities below — and the two are separate pieces of evidence, not one.

## Step 2. The screen

### The conformational screen — the column that decides, and the one that still stands

Seven chemistries are expressible. **CNEPO is out of scope**: its epoxide bridges two backbone
carbons, which is a change of chain topology rather than a substituent, and
`polymers.BackboneAtom` carries pendants, not rings (`docs/CHEMISTRY_EXTENSION.md` phase 4). It
is reported as out of scope rather than approximated by something else.

Method per chemistry: `fit_ris` (step 10 deg, 6 monomers, third order, ~2900 conformers),
exhaustive `enumerate_periodic` to period 8, then packing of the top three by RIS energy **plus**
all-trans, TG+TG- and T3GT3G' whatever they rank — because screening only the top of the RIS
list lets a polymer look useless when its polar phase happens to rank 63rd, which for PVDF it
does.

**The scan's backbone angles are no longer frozen for the chemistries that need them free.**
`docs/NITRILE_LANDSCAPE.md` found that the frozen-angle one-bond scan has wells a chain that can
bend does not have — VDCN's at ±120 deg, 8.7 kcal/mol below planar trans — and that `fit_ris`'s
basin assignment folded that well into VDCN's trans state, so its published "all-trans ground
state at −11.48" was a chain built at 180 deg carrying the energy of a well at 120. `fit_ris`
now relaxes every backbone angle of the oligomer at each scan point (`angles="relaxed"`) for a
chemistry whose rigid profile has a minimum no relaxed minimum lies within 30 deg of, and the
answer is measured per chemistry rather than assigned (`forcefield.ANGLE_RELAXATION_DEFAULTS`,
`examples/angle_relaxation_defaults.py`): PVDF and PE keep the rigid scan, whose numbers every
table here depends on and which are bit-for-bit what they were; PVDC, CFE, CDFE, AN and VDCN
get the relaxed one, at 9-13 min per chemistry against 1 s. The rows below are the regenerated
ones; the previous values are kept in the last column so that the change is on record.

| polymer | scan | conformations | ground state | all-trans rank | ΔE(all-trans) | chain μ⊥ per monomer | previously (rigid) |
|---|---|---|---|---|---|---|---|
| VDCN | relaxed | 67 | TTG+G+ (−5.45) | 43 of 67 | **+4.10** | 0.645 e.A | TT (−11.48), **1 of 83, +0.00** |
| PVDC | relaxed | 65 | TTG+G+ (−4.64) | not enumerated | +3.13 | — | TTG+TTG+ (−8.80), +2.85 |
| PVDF | rigid | 84 | TG+TG+G+TG+G+ (−4.65) | 63 of 84 | +4.58 | 0.362 e.A | unchanged |
| AN | relaxed | 148 | TG−G−G− (−7.70) | 136 of 148 | +5.91 | 0.561 e.A | TG− (−13.05), 105 of 164, +6.14 |
| CFE | relaxed | 95 | G+G+ (−23.26) | 91 of 95 | +21.22 | 0.386 e.A | TG− (−30.90), 41 of 157, +11.16 |
| CDFE | relaxed | 147 | G+G+G+G−G−G+ (−19.07) | 146 of 147 | +17.71 | 0.325 e.A | TG+G+G+ (−30.42), 95 of 162, +16.67 |
| FANOME | — | not screened | — | — | — | — | — |

ΔE is kcal/mol per monomer above the conformational ground state; it is a *difference*, so the
known badness of the all-trans RIS reference cancels out of it even where the absolute energies
are meaningless (CFE's and CDFE's rigid −30 were relative to a state carrying 162 and 199 kcal/mol
of all-trans strain; relaxed, the reference is a chain that can bend and their absolute energies
are −23 and −19). μ⊥ is the **all-trans** repeat's transverse dipole per monomer, computed with
the fitted bond-charge increments rather than the polymer entry's illustrative charges; it is the
component that survives in a planar zigzag and cancels in a screw helix, so it is the one
β-PVDF's ferroelectricity is made of. PVDC has no entry because it has no all-trans repeat to
measure one on; its ΔE is that of the rebuilt all-trans under the model.

**What the regeneration did.** VDCN is the row that flips: with its trans state no longer
carrying the ±120 deg well (e1(T) goes from −8.7 to 0.0 on both bond types, the T basin minimum
sits at 180) its all-trans falls from first of 83 to 43rd of 67, 4.1 kcal/mol per monomer above a
TTG+G+ chain — the same verdict `docs/NITRILE_LANDSCAPE.md` reached with a per-*type* angle
relaxation (74th, +5.5), now with the per-*atom* one it asked for, which also puts the gauche well
at ±40 deg rather than that experiment's ±95. AN's verdict stands (a gauche-rich helix ground
state, all-trans 5.9 above) and its trans state is a basin edge at both levels, flagged by the
fit (`FitReport.edge_states`), so its ΔE is an upper bound on a well that is not there. PVDC's
ground state moves from TTG+TTG+ to TTG+G+ and its ΔE from +2.85 to +3.13 — neither changes the
geometric verdict below. CFE's and CDFE's all-trans move *up*, to 91st of 95 and 146th of 147:
their rigid trans states were ±140 and ±150 deg basin edges 20 kcal/mol below planar, and relaxed
they are −2.4 and −5.6.

**A caveat that applies to both scans, measured rather than assumed.** An RIS term is a basin
*minimum* — the pair term at its own best pair of angles, relaxed angles included — while a chain
is built at one adapted angle per state, so the model sits below the chain it describes. Against
directly relaxed periodic chains of their own top sequences the relaxed models are 3-7 kcal/mol
per monomer low for VDCN, AN and PVDC; the rigid PVDF model is 2.8 low on the same test against
rigid chains, so this is the fit's convention, not the relaxation. For CFE and CDFE the gap is
20-30 kcal/mol per monomer at *either* level (the rigid models are 14-42 low against rigid
chains): the chlorinated helices have non-additivity beyond three bonds that no third-order
model carries, and their rows — as the previous version of this document already said of their
absolute energies — cannot be read as a ranking. Their ΔE(all-trans) is quoted because it is what
the model says, not because the model is trusted there.

**This column's differences are still well clear of the error bar** — +3.13 to +21.22 kcal/mol
per monomer against 0.27 — but one of its rows had to be withdrawn and regenerated, and the
lesson is the one above: the rigid scan's numbers were a property of the frozen geometry for the
bulky-pendant chemistries. Whatever else this screen cannot decide, it can say whether a polar
chain conformation is accessible, and that is the criterion that separates the chemistries.

Two entries need their own sentence.

**PVDC has no polar conformation at all, and the reason is geometric rather than energetic.** Its
measured backbone angles are unequal — 123 deg at CH2, 114 deg at CCl2, both crystallographic —
so an ideal-angle repeat carries 123 − 114 = 9 deg of curl and every sequence is a circular arc
with zero rise. `enumerate_periodic` drops all-trans and TG+TG- both, on the rise-per-bond
filter. The screen rebuilds them anyway and searches for the torsion deflection that closes them:
**all-trans still does not close, TG+TG- does**, and the deflected glide is the only thing PVDC
packs. So PVDC has a crystal and it is not a polar one. **PVDC is not a candidate**, and the
model says so for the same reason the crystallography does — the real polymer is a TGTG' glide,
not a planar zigzag.

The closure is one equation in four unknowns, so its solutions are a manifold and which point of
it you land on depends on where you start. From the *fitted* model's own adapted state angles the
screen reaches 177.5 / 27.7 deg and, after refinement, c = 4.76 A against a measured 4.68. From
the ideal 180 / 60 it reaches 174.7 / 54.7 with c = 4.70, and `tests/test_pvdc.py` reaches
175.3 / 49.4 with c = 4.677 — that last being the published 175 / 49 and 4.68. The *fibre repeat*
is robust across all three to about 2 %; the torsion pair is not, and only the ideal-angle start
recovers it.

**FANOME could not be screened at all.** Its all-trans reference is an overlapping structure —
methyl hydrogens of methoxy groups on consecutive substituted carbons 0.80 A apart, ~1e6 kcal/mol
on a ten-bond oligomer, robust to the frozen C-O rotamer — so `fit_ris` would measure every RIS
energy from a state that is not a molecule. It is registered, its geometry is tested, and it is
deliberately not fitted (`docs/CHEMISTRY_EXTENSION.md` phase 3). Reporting a number for it would
be reporting noise.

### Packing: cells, densities and the polymorph ordering

Packing follows `fitting.predict`'s protocol, the one the α-below-β validation was measured with:
exhaustive table screen, exact polish, the antipolar branch searched explicitly, then refinement
of torsions, backbone angles and cell with the bond-angle strain added back into the reported
energy. The cells, densities and `ΔE lattice` below are the **truncated**-sum refinement, the
protocol the first version of this table used, so that only the polarity column changes its sum.
`ΔE lattice` is kcal/mol per monomer above each chemistry's own best packed phase.

| polymer | phase | a × b × c (A) | ρ | ΔE lattice | \|P\| |
|---|---|---|---|---|---|
| PVDF | TG+TG- (α) | 4.99 × 9.04 × 4.71 | 2.000 | +0.000 | 0.062 |
| PVDF | T3GT3G' (γ) | 5.28 × 8.89 × 9.34 | 1.939 | +0.470 | 0.069 |
| PVDF | TT (β) | 4.54 × 8.55 × 2.60 | 2.107 | +0.747 | 0.115 |
| PVDC | TG+TG- (deflected) | 6.23 × 11.23 × 4.76 | 1.932 | +0.000 | 0.023 |
| CFE | TG+TG- | 5.49 × 10.69 × 4.74 | 1.924 | +0.000 | 0.040 |
| CFE | T3GT3G' | 6.29 × 9.58 × 9.20 | 1.927 | +0.834 | 0.048 |
| CFE | TT | 5.49 × 8.25 × 2.69 | 2.190 | +12.981 | 0.115 |
| CDFE | TG+TG- | 5.53 × 10.95 × 5.03 | 2.149 | +0.000 | 0.007 |
| CDFE | T3GT3G' | 5.64 × 11.09 × 10.06 | 2.077 | +0.711 | 0.000 |
| CDFE | TT | 6.15 × 8.29 × 2.69 | 2.381 | +10.252 | 0.076 |
| AN | TG+TG- | 5.86 × 10.22 × 4.75 | 1.238 | +0.000 | 0.023 |
| AN | T3GT3G' | 6.19 × 10.39 × 9.31 | 1.178 | +1.141 | 0.048 |
| AN | TT | 6.63 × 6.62 × 2.66 | 1.510 | +5.106 | 0.025 |
| VDCN | TG+TG- | 6.17 × 11.55 × 4.72 | 1.543 | +0.000 | 0.059 |
| VDCN | T3GT3G' | 8.43 × 10.04 × 9.27 | 1.321 | +3.594 | 0.041 |
| VDCN | TT | 9.33 × 5.96 × 2.69 | 1.731 | +7.548 | 0.000 |

**The PVDF rows are the calibration and they behave.** The polymorph ordering comes out
α < γ < β, with `E(α) − E(β) = −0.747 kcal/mol =` **−3.13 kJ/mol per monomer** — inside the
−6.5 .. −2.6 kJ/mol range four independent studies across five functionals agree on, and
consistent with the −2.96 recorded for this preset in `fitting.py`. That is acceptance test 1
passing, through a different code path than the one that recorded it, and **it is a 0.747 against
a 0.27 error bar, so it is one of the few energy differences here that the model resolves.**
Cells against experiment (4.91 × 8.58 × 2.56, 4.96 × 9.64 × 4.62, 4.96 × 9.67 × 9.20): every one
of the nine edges within 8 %, five within 2 %. The worst are β's `a` (4.54 against 4.91, 7.5 %)
and γ's `b` (8.89 against 9.67, 8.1 %).

Three rows moved from the previous version of this table, all for the same reason: where the old
antipolar gap came out negative, the old script refined its "antipolar" cell, which was not one.
AN's TT was 11.27 × 4.20 × 2.66 at ρ 1.403 and is now 6.63 × 6.62 × 2.66 at 1.510; VDCN's TT was
6.06 × 10.57 × 2.69 at 1.504 and is now 9.33 × 5.96 × 2.69 at 1.731; VDCN's T3GT3G' moved
likewise. Every other cell and density is unchanged to the quoted digits.

### Polar or antipolar: the column the screen was built around, against its own error bar

The quantity is `E(best exactly-antipolar cell) − E(best cell overall)` per monomer, both branches
rigid, two chains, gamma = 90°, both under the same electrostatic sum, the antipolar branch
constructed from the chain's own moment by `antipolar_offsets` and searched by
`antipolar_cell_exact` inside the bounds the polar branch is screened in. Negative means the
crystal prefers antipolar and the net polarization is zero. It is judged on the **rigid**
comparison, which is the one that means something: the antipolar subspace is not stationary, so a
refinement started inside it slides back out and the refined cell's own flip flag is not the
answer. The zero dipole is checked at every answer rather than assumed — max |P| over every row
below is 7e-17 C/m².

Ewald, tinfoil, `pvdf-dft-valence`. **Every gap is printed beside the model's own resolution, and
the verdict is "not resolved" whenever it is inside it.**

| chemistry | phase | μ⊥ (e.A/mon) | gap | ± | verdict |
|---|---|---|---|---|---|
| pvdf | **TT (β)** | 0.362 | **+1.021** | 0.27 | **polar** |
| pvdf | TG+TG- (α) | 0.232 | +0.029 | 0.27 | *not resolved* |
| pvdf | T3GT3G' (γ) | 0.232 | +0.603 | 0.27 | polar |
| pvdc | TG+TG- (deflected) | 0.110 | +0.397 | 0.27 | polar |
| cfe | **TT** | 0.386 | **+1.501** | 0.27 | **polar** |
| cfe | TG+TG- | 0.136 | +0.247 | 0.27 | *not resolved* |
| cfe | T3GT3G' | 0.136 | +0.868 | 0.27 | polar |
| cdfe | **TT** | 0.325 | **+0.390** | 0.27 | **polar** |
| cdfe | TG+TG- | 0.235 | +0.269 | 0.27 | *not resolved* |
| cdfe | T3GT3G' | 0.235 | +0.251 | 0.27 | *not resolved* |
| an | **TT** | 0.561 | **+0.046** | 0.27 | ***not resolved*** |
| an | TG+TG- | 0.354 | +0.530 | 0.27 | polar |
| an | T3GT3G' | 0.354 | +1.280 | 0.27 | polar |
| vdcn | **TT** | 0.645 | **−0.645** | 0.27 | **antipolar** |
| vdcn | TG+TG- | 0.196 | +3.053 | 0.27 | polar |
| vdcn | T3GT3G' | 0.196 | +2.046 | 0.27 | polar |
| pvf-cand | **TT** | 0.302 | **+0.205** | 0.27 | ***not resolved*** |
| pvf-cand | TG+TG- | 0.181 | +0.526 | 0.27 | polar |
| pvf-cand | T3GT3G' | 0.181 | +0.795 | 0.27 | polar |
| trfe-cand | **TT** | 0.305 | **+0.051** | 0.27 | ***not resolved*** |
| trfe-cand | TG+TG- | 0.213 | +0.406 | 0.27 | polar |
| trfe-cand | T3GT3G' | 0.213 | +0.574 | 0.27 | polar |
| vfcn-cand | **TT** | 0.548 | **+0.265** | 0.27 | ***not resolved*** |
| vfcn-cand | TG+TG- | 0.214 | +1.083 | 0.27 | polar |
| vfcn-cand | T3GT3G' | 0.214 | +1.343 | 0.27 | polar |
| vclcn-cand | **TT** | 0.556 | **+0.134** | 0.27 | ***not resolved*** |
| vclcn-cand | TG+TG- | 0.193 | +1.791 | 0.27 | polar |
| vclcn-cand | T3GT3G' | 0.193 | +1.249 | 0.27 | polar |

**How to read the totals, because the honest number depends on which question is asked.**

*Across all 28 phases, 19 resolve.* That sounds better than it is. Eighteen of those nineteen say
"polar", and **fifteen of them are helical** phases whose transverse chain moment is 0.11 to 0.35
e.A per monomer and whose cells carry `|P|` of 0.007 to 0.069 C/m² — a polarity verdict on a cell
that is barely polar either way, and a screw helix cancels its transverse dipole over the
crystallographic repeat by construction. Whether α-PVDF or a TG+TG- CFE packs polar or antipolar is
not what the campaign is asking.

*For the all-trans planar zigzag — the only conformation that can carry a ferroelectric
polarization — the model resolves **four of nine** and cannot decide five.* Resolved: PVDF polar
by 1.02, CFE polar by 1.50, CDFE polar by 0.39, VDCN antipolar by 0.65. Not resolved: AN (+0.046),
trfe-cand (+0.051), vclcn-cand (+0.134), pvf-cand (+0.205), vfcn-cand (+0.265). **That is the
screen's primary output and it is a "cannot tell".** Two of the five sit below a fifth of the error
bar; vfcn-cand's +0.265 misses the threshold by 2.5 %, which is itself a reason not to read a
verdict off a number that close.

*And VDCN is the only antipolar verdict anywhere in the screen.* The previous version of this
document reported five.

### Ewald is not what decides any of this

Both branches were also measured with the truncated 8 Å damped-shifted-force sum on an identical
protocol, so that exactly one thing differs. Over the nine all-trans rows:

| chemistry | μ⊥ | truncated gap | Ewald gap | Ewald − truncated | verdict, both |
|---|---|---|---|---|---|
| pvdf | 0.362 | +1.018 | +1.021 | +0.003 | polar |
| cfe | 0.386 | +1.524 | +1.501 | −0.023 | polar |
| cdfe | 0.325 | +0.410 | +0.390 | −0.020 | polar |
| an | 0.561 | +0.065 | +0.046 | −0.019 | *not resolved* |
| vdcn | 0.645 | −0.558 | −0.645 | −0.087 | antipolar |
| pvf-cand | 0.302 | +0.175 | +0.205 | +0.030 | *not resolved* |
| trfe-cand | 0.305 | +0.034 | +0.051 | +0.017 | *not resolved* |
| vfcn-cand | 0.548 | +0.282 | +0.265 | −0.017 | *not resolved* |
| vclcn-cand | 0.556 | +0.103 | +0.134 | +0.031 | *not resolved* |

And on PVDF's three polymorphs, where the same quantity is measured under both sums through the
acceptance-test path as well: β +1.018 / +1.021, α +0.084 / +0.094, γ +0.584 / +0.603.

**Ewald moves every polarity gap by at most 0.09 kcal/mol per monomer and changes no verdict
anywhere.** A truncated dipole sum still cannot be defended as a matter of principle, and the
boundary convention still changes β's lattice energy by 1.17 kcal/mol per monomer (below). But on
these structures the truncated and the correct sum agree on the sign of every gap and nearly on
its size, and neither is what limits the answer. **The error bar is** — it is three to thirty
times the Ewald correction.

### The response figures — the weak half, flagged

Computed on the **polar** structure whether or not the polymer would adopt it, because the
question is what the polar phase would deliver. `d_33` and `d_31` are in the film convention
(3 = poling = the crystal's own polar axis, 1 = draw = the chain axis), generalised from
`fit_charge_flux.py`'s hard-coded x-axis to the actual polarization direction; for β-PVDF that
reduces to the same numbers, which is the check.

These columns come from the **truncated** sum and are the same measurement as before — nothing in
the antipolar correction touches them, since they never used the antipolar branch. Ewald is
deliberately not applied to them: `docs/BENCHMARK.md` records that under Ewald β's cell moves by
0.01 Å, its `|P|` by 0.0001 C/m² and `C_33` by 0.1 %, while γ's two independent routes to `d` stop
agreeing and the row would have to be reported as a non-measurement. Correct electrostatics buys
nothing here and costs one row.

| polymer | phase | \|P\| C/m² | C11 | C22 | C66 | C33 | d33 | d31 | blocking GPa | free strain | work kJ/m³ | routes | status |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| PVDF | TT | 0.150 | 28.4 | 20.1 | 5.1 | 330 | −1.88 | +2.33 | 0.072 | 2.3e−4 | 4.22 | 0.33 % | ok |
| CFE | TT | 0.168 | 32.9 | 27.6 | 6.1 | 636 | −1.20 | +1.36 | 0.021 | 1.8e−3 | 5.57 | 2.7 % | ok |
| CDFE | TT | 0.116 | 19.1 | 14.2 | 6.0 | 614 | −0.40 | +0.96 | 0.008 | 1.7e−3 | 6.64 | 1.1 % | ok |
| VDCN | TG+TG- | 0.074 | 36.9 | 25.2 | 16.3 | 191 | −3.94 | +1.71 | 0.020 | 3.9e−4 | 2.17 | 0.27 % | ok |
| AN | T3GT3G' | 0.058 | 9.8 | 19.5 | 7.1 | nan | — | — | 0.003 | 4.1e−4 | 0.43 | 0.37 % | **rigid path** |
| PVDC | TG+TG- | — | — | — | — | — | — | — | — | — | — | — | **failed** |

Elastic constants in GPa. "routes" is the disagreement between the two independent routes to `d` —
a dipole derivative and a zero-stress root find — the internal check that the flux entered the
energy gradient and the dipole derivative consistently. PVDF's row reproduces
`docs/BENCHMARK.md` exactly (`d_33 = −1.88`, `d_31 = +2.33`, `C_33 = 330` against a literature DFT
287), which is what says the generalised film convention and the packing protocol are the same
calculation the benchmark ran; CFE's, CDFE's and the candidates' rows reproduce the previous
table to the quoted digits.

**AN's and VDCN's rows changed phase, and that is a consequence of the polarity correction rather
than of the mechanics.** The script picks the response target as the packed phase with the largest
`|P|`. Once VDCN's all-trans correctly packs antipolar its `|P|` is zero, and once AN's all-trans
stops being refined from a mis-built antipolar cell its `|P|` falls below its γ-type phase's, so
both targets move off the zigzag and onto a helix. AN's lands on a sequence with no line group,
degrades to the rigid path where `eps_zz` is not a variable, and is reported as such. Neither
chemistry's response was ever load-bearing — 0.43 and 2.17 kJ/m³, the bottom of the table — but
the target rule is a heuristic and it should be read as one.

**PVDC's response could not be computed, in three documented steps.** No conformation closes at
ideal angles, so its chain has deflected torsions; a deflected chain breaks the glide pattern its
sequence is named for, so it has no line group and no `Shape`, so the deformable path is
unavailable; and the rigid fallback then fails its own strained-cell construction. Two of those
three are model limitations already on the roadmap (variable backbone angles, and a shape
parametrisation that does not need an exact symmetry pattern). It is reported as a failure rather
than forced.

**One methodological warning the screen produced by itself.** The work-density figure is dominated
by the softest reachable mode. The pvf candidate's C66 relaxes to 1.04 GPa, its shear `d` blows up
to 192 pC/N and its work density to 95.8 kJ/m³ — twenty times anything else — and the two routes
to `d` then disagree by 67 %, which is what caught it. It is reported as a non-measurement. Any
work density quoted from this package should be read together with its softest constant and its
route agreement.

### Wall time

The polarity column is now the cost. `antipolar_cell_exact` screens the antipolar subspace with a
grid whose resolution is set in angstroms, so the grid grows with the repeat, and its polish is
the exact kernel: 13 s for a one-monomer all-trans cell, 40–180 s for a two-monomer TG+TG-, and
480–580 s for a four-monomer T3GT3G'. Per chemistry, under Ewald, one process:

| polymer | fit | enumerate | pack + polarity | response | total |
|---|---|---|---|---|---|
| PVDF | 0.6 | 2.8 | 705.0 | 6.5 | **714.9 s** |
| PVDC | 0.4 | 3.6 | ~200 | 18 | ~225 s |
| AN | 1.3 | 10.2 | 907.6 | 96.5 | **1015.6 s** |
| VDCN | 1.7 | 9.7 | 1013.5 | 30.7 | **1055.5 s** |
| trfe-cand | 1.3 | 9.6 | 900.3 | 11.4 | **922.8 s** |

Six physical cores, `OMP_NUM_THREADS=1`, `POLYFIND_TABLE_PROCS=1`, one chemistry per process,
cold interaction-table cache throughout (`table_cache_dir=None`, deliberately — a disk table is
not keyed on the potential). The fit and the enumeration are 1 to 10 s; **the antipolar subspace
search is 80 % to 95 % of the cost**, and unlike the previous version of this table that cost buys
something: the four-point axial scan it replaces was aliased to the monomer period on every
four-monomer repeat here. Ten chemistries in about an hour of wall time on six cores.

The fit column above is the rigid scan's. With the angle-relaxed scan the bulky-pendant
chemistries now default to, the fit is 9-13 min per chemistry (measured with eight fits sharing
twelve cores: VDCN 750 s, AN 534 s, PVDC 568 s, CFE 657 s, CDFE 693 s, trfe-cand 592 s, vfcn-cand
718 s, vclcn-cand 655 s; 1.3-1.8 million calculator energies each, every one of the 2881 scan
conformers being a bounded minimisation over its 15 backbone angles), so for those chemistries the
fit is now comparable to the antipolar search rather than negligible beside it. PVDF, PE and
pvf-cand keep the 1 s rigid fit.

## The three acceptance tests, re-measured

These three have tracked this project since `DESIGN.md` 5.4 and none of them is in any fit's
objective. They stood at **one pass, two failures**. They are re-measured here under the 2x2 of
{truncated, Ewald} × {the old antipolar construction, the corrected one}, so that a change can be
attributed rather than asserted.

1. **α below β by 2.6 to 6.5 kJ/mol per monomer** — four independent studies across five
   exchange-correlation functionals agree on the sign and roughly on the magnitude.
2. **The isolated-chain RIS ranking does not put the TG+ 3/1 helix first** — PVDF does not form it.
3. **α's packed cell is antipolar** — the antipolar cell must not cost more than the polar one.

Reproduced from the shipped code, not from a one-off script —
`fitting.acceptance_tests(FITTED_VALENCE, coulomb=..., ewald=EwaldSpec(), antipolar=...)`, whose
defaults (`"dsf"`, `None`, `"legacy"`) are the fit's own and give the first row:

| sum | antipolar construction | E(α)−E(β) kJ/mol per monomer | test 1 | α antipolar gap | test 3 | β gap, as a control |
|---|---|---|---|---|---|---|
| truncated | old (`antipolar_cell`) | −3.13 | **PASS** | +0.084 | FAIL | +0.0000 (`\|P\|` = 0.115: not antipolar) |
| truncated | corrected | −3.13 | **PASS** | +0.084 | FAIL | +1.018 |
| Ewald | old | −3.11 | **PASS** | +0.095 | FAIL | +0.0000 (`\|P\|` = 0.115: the same defect) |
| Ewald | corrected | −3.11 | **PASS** | +0.094 | FAIL | +1.021 |

**They still stand at one of three, and nothing moved them.**

**Test 1 passes, unchanged, and is attributable to neither correction.** −3.13 kJ/mol per monomer
truncated, −3.11 under Ewald: Ewald shifts it by 0.02 kJ/mol (0.005 kcal/mol per monomer) and the
antipolar construction does not enter it at all, because α's polar cell is the one that gets
refined on every row. It is the one test whose margin is outside the error bar: 0.747 kcal/mol per
monomer against 0.27.

**Test 2 fails, unchanged, and structurally cannot be affected by either.** The isolated-chain RIS
ranking puts `TG+` first at −4.368 kcal/mol per monomer, ahead of `TG+G+G+` at −3.940 and
`TTG+G+` at −3.443 — a margin of 0.428, which is the figure `DESIGN.md`'s own acceptance
table records for this preset. No lattice sum
and no cell construction enter this test: it is a single chain. Attributable to neither, and it
could not have been otherwise.

**Test 3 fails, and the interesting part is by how little.** α's antipolar cell costs +0.094
kcal/mol per monomer more than its polar one under Ewald with the corrected construction. The
attribution is clean and it is nearly nothing in both directions:

* **the antipolar fix moves α by 0.000.** The old construction was *already correct for α*: the α
  helix's `m_x` is exactly zero, so "flip with equal setting angles" really is the antipolar
  subspace there (its cell's `|P|` is 1.3e-16). α is the one case the old helper got right, which
  is why `fitting.antipolar_cell` is left untouched for the fit's use.
* **Ewald moves α by +0.011**, from +0.084 to +0.095, in the *wrong* direction.
* the screen's own protocol gives +0.027 / +0.029 for the same quantity, differing from the
  acceptance test's only in whether the polar branch's screen table is built inside the
  `applied()` block. 0.06 of protocol, 0.01 of electrostatics, 0.00 of construction.

**So the honest statement about test 3 is not that the model gets α wrong; it is that the model
cannot decide α.** +0.027 to +0.095 against a 0.27 error bar is a third of the noise at most. The
same applies to any gap below about 0.3 per monomer, which is five of the nine all-trans rows
above. The test is written as a sign test and the model does not have the resolution to take it.

**What the corrected construction *did* change is the control column.** β's gap under the old
helper is `+0.0000` with `|P| = 0.115 C/m²` — it returns β's own polar minimum — and under the
corrected one it is +1.018. That is the defect, measured, on the structure it matters for.

## Step 3. What drives the answer, re-derived

Correlations of each structural feature against the predicted work density, over the **eight**
chemistries that produced a usable response (five incumbents and three candidates; PVDC failed and
pvf-cand is a non-measurement). **At this n these describe the sample; they do not establish
anything** — and this re-run is the demonstration of that, because the sample changed by three rows
and half the coefficients moved with it. The "as published" column is the previous version of this
table, on a sample of eight that overlaps this one in five rows.

| feature | Pearson | Spearman | as published | reading |
|---|---|---|---|---|
| packing density ρ | **+0.788** | **+0.857** | +0.81 / +0.91 | the strongest, and the most likely to be a confound |
| ΔE(all-trans above conformational ground) | **+0.793** | +0.571 | +0.63 / +0.43 (+0.648 / +0.429 on the rigid-scan column) | the response column rewards inaccessible phases — more clearly on the regenerated column |
| axial stiffness C33 (n = 7) | +0.588 | +0.500 | +0.35 / +0.19 | no rule |
| polarization \|P\| | +0.422 | +0.429 | +0.45 / +0.26 | no rule |
| transverse stiffness C11 | +0.415 | +0.381 | −0.18 / −0.29 | no rule, and it changed sign |
| gap to the nearest polymorph | +0.363 | +0.024 | +0.16 / −0.10 | no rule |
| transverse chain dipole μ⊥ | **−0.162** | **−0.119** | **−0.69 / −0.76** | **the anti-correlation is gone** |

Two of these moved enough to matter, and the sample is why: AN's response row is now its γ-type
phase and VDCN's its α-type, because the polarity correction moved which packing carries the
largest `|P|` (see the response table), and vclcn-cand's row is a measurement this time. So the
n = 8 sample is not the same n = 8, and **anything that moves when three of eight rows are
exchanged was never a rule.** That applies to C11's sign flip and, more importantly, to the μ⊥
anti-correlation, which was the headline of the published rule 1 at −0.69 / −0.76 and is −0.162 /
−0.119 here.

### Rule 1 — "a perpendicular dipole above 0.5 e.A per monomer forces antipolar packing". Falls.

As published it was exceptionless: every chemistry whose all-trans carries μ⊥ above 0.5 packs that
all-trans antipolar, four for four — VDCN 0.645, AN 0.561, vclcn 0.556, vfcn 0.548, with gaps of
−7.60, −1.86, −1.46 and −3.09. **All four of those gaps were comparisons between two polar cells,
and under two different potentials.** Re-measured:

| chemistry | μ⊥ | gap as published | gap, corrected (Ewald) | verdict |
|---|---|---|---|---|
| VDCN | 0.645 | −7.60 | **−0.645** | **antipolar** |
| AN | 0.561 | −1.86 | +0.046 | *not resolved* |
| vclcn-cand | 0.556 | −1.46 | +0.134 | *not resolved* |
| vfcn-cand | 0.548 | −3.09 | +0.265 | *not resolved* |

So the rule holds for **one chemistry of four** — and the part that matters more is that the three
that changed sides did not change to "polar". They changed to **"not resolved"**. A rule supported
by one case out of four, with the other three inside the model's own error bar, is not a rule.
**It must not be used as a reject filter**, and not because the electrostatics were truncated:
because the quantity it was fitted against was not measuring what it was named after, and because
the corrected quantity is below this model's resolution for three quarters of the cases it rested
on.

What is left of it is a sign, and only a sign. Pearson(μ⊥, all-trans polarity gap) over the nine
chemistries with an all-trans packing is **−0.511**, Spearman −0.450, against −0.49 / −0.51
recorded for the same quantity before — so that much is stable. And **none of the four
large-dipole chemistries comes out polar**: one is antipolar and three are unresolved. So the
tendency is not contradicted; it is *untested*, which is a different thing from supported, and a
filter that rejects a candidate has to be better than untested. The one-sided form the previous
version recommended — "μ⊥ above 0.5 is a reject filter, nothing below it is an accept filter" — is
withdrawn in both halves: the reject half rests on one case, and the accept half never had
support.

### Rule 2 — "ranking on work density alone inverts the answer". Stands.

It never used the antipolar branch. Its evidence is ΔE(all-trans above the conformational
ground state) against work density, **+0.648 Pearson / +0.429 Spearman**, against +0.63 / +0.43
as published — stable across a sample in which three of eight rows were exchanged. CDFE and CFE
still top the work-density table at 6.64 and 5.57 kJ/m³ while their polar phases sit 10.25 and
12.98 kcal/mol per monomer above their own lattice ground states, and trfe-cand is fourth of eight on work density
while being first of ten on accessibility — which is the inversion the rule is about. A screen that
ranks on work density alone picks the chemistries whose polar phase the polymer would never adopt,
and the accessibility column has to gate it. This is the one design rule of the original three that
survives, and it survives because
it rests on the two columns that are resolved — the conformational ΔE and the refined lattice
ranking — rather than on the polarity gap.

### Rule 3 — "denser packing helps". Demoted to an observation.

ρ is still the strongest single correlate, **+0.788 / +0.857** against +0.81 / +0.91 as
published, and it is the only correlation in the table that both survived the sample change and
exceeds 0.7. There is a plausible mechanism — work density is dipole per unit volume worked against
stiffness, and denser packing raises both — but in this sample density is nearly collinear with
halogen content, so it may be reading "has chlorine" rather than "is dense". Not to be used alone,
and at this n not to be stated as a rule.

### Newly supported: accessibility is the column to screen on, because it is the one that resolves

This is the rule the corrected data does support, and it is a rule about the method rather than
about the chemistry. Of the four structural columns this screen produces — conformational ΔE,
refined lattice ΔE, the polarity gap, and the response coefficients — the first two have
differences of 0.4 to 17 kcal/mol per monomer against a 0.27 error bar, the third has differences
of 0.03 to 3.05 of which the decision-relevant ones are mostly inside it, and the fourth is a
coarse ranking no finer than a factor of 1.7. **Screen on accessibility, gate with it, and treat
the polarity column as a flag for "compute this properly later" rather than as a verdict.** On
this sample that rule alone reproduces every conclusion the screen is willing to stand behind: it
rejects PVDC, rejects both nitriles, and ranks trifluoroethylene first.

### What the data does not support

No rule from the axial stiffness (+0.588 / +0.500, and it was +0.35 / +0.19 on the
previous sample — a correlation that nearly doubles when three rows are exchanged is describing the
rows), the transverse stiffness (+0.415, previously −0.18 — it changed *sign*), the polarization
magnitude (+0.422 / +0.429) or the energy gap to the nearest competing polymorph (+0.363 / +0.024).
Stating a chain-stiffness rule from these numbers would be inventing one. The polymorph-gap result
is the most disappointing, because a large gap to the nearest competitor is exactly what would make
a polar phase *stable*, and this sample cannot see it.

**And one thing that is not a rule but is worth stating as a warning.** Four of the seven
correlations above moved by more than 0.2 in Pearson when three of eight rows were exchanged, and
the transverse stiffness changed sign (two did in Spearman). At n = 8 a correlation coefficient is
a description of eight points, and this re-run measured how little that survives. The two that
did not move — density and accessibility — are the two this document is willing to quote, and even
those are quoted as directions rather than as coefficients.

## New candidates, re-assessed

**These are model suggestions, not predictions.** They are ordinary `Polymer` values constructed
inside `examples/screen_electroactive.py` and passed straight to the funnel — nothing in the
package changes and no default moves — with geometry and charges chosen exactly the way PVDF's,
PVDC's, CFE's and CDFE's were. The charges are illustrative and the potential is PVDF's.

They were proposed from three design rules, two of which have since fallen, so this is a
re-assessment rather than a recommendation carried forward.

| candidate | formula | ΔE(TT) above conformational ground | TT above lattice ground | TT polarity gap | verdict | ρ | \|P\| | work |
|---|---|---|---|---|---|---|---|---|
| **trfe-cand** | −(CHF−CF2)− | +5.41 (relaxed scan; was +3.96) | **+0.424** | +0.051 | *not resolved* | 2.310 | 0.108 | 2.49 |
| pvf-cand | −(CH2−CHF)− | +4.06 (rigid) | +0.984 | +0.205 | *not resolved* | 1.654 | 0.114 | not a measurement |
| vfcn-cand | −(CH2−C(F)(CN))− | +8.62 (relaxed; was +6.71) | +4.072 | +0.265 | *not resolved* | 1.824 | 0.166 | 2.28 |
| vclcn-cand | −(CH2−C(Cl)(CN))− | +4.13 (relaxed; was +9.89) | +16.342 | +0.134 | *not resolved* | 1.880 | 0.157 | 1.75 |
| *PVDF, for comparison* | −(CH2−CF2)− | +4.58 (rigid) | +0.747 | +1.021 | polar | 2.107 | 0.150 | 4.22 |

The ΔE(TT) column is the regenerated one: the candidates are not registered, so `fit_ris`'s
`angles="auto"` measured each on the fly, and three of the four have frozen-angle wells the
relaxed profile lacks (trfe-cand's at ±30 deg, +14 kcal/mol *above* trans and 25 deep behind a
barrier; vfcn-cand's at ±80; vclcn-cand's at ±110 and ±130, the latter 20 kcal/mol below trans)
while pvf-cand does not and keeps the rigid scan. The "TT above lattice ground" column is
packing, measured with the angles already relaxing, and was not re-run here; what this change
would alter there is which three RIS conformers are carried into packing, not the energy of a
packed cell.

Three of the four densities moved from the published table (pvf-cand 1.689 → 1.654, vfcn-cand
1.744 → 1.824, vclcn-cand 1.737 → 1.880) for the same reason AN's and VDCN's cells did: their
published all-trans gap was negative, so the old script refined a cell that was not antipolar.
trfe-cand's 2.310 is unchanged, because its gap was positive either way.

**Not one candidate's all-trans polarity is resolved.** That is the single most important line in
this section: the column every one of these four was recommended or rejected on is, for all four
of them, inside the model's error bar.

### The trfe-cand recommendation is withdrawn as it was stated

The previous version of this document recommended trifluoroethylene on two grounds, and they have
to be separated because one survives and one does not.

**Withdrawn: "the only chemistry here whose every packed phase prefers polar over antipolar by a
wide margin — +1.66, +1.87 and +1.96 for TT, TG+TG- and T3GT3G', the largest *minimum* of any
chemistry in this screen, so trfe is the one chemistry with no antipolar escape route."** All three
numbers came from `antipolar_cell`, which for a planar zigzag compares two polar cells, under a
packer carrying a different potential from the one the polar branch used; and the two helical
numbers were never re-measured when that was found. They have been re-measured now, with the
corrected construction and under Ewald: **+0.051 for TT, +0.406 for TG+TG-, +0.574 for T3GT3G'.**
The two helices are resolved polar, which is the direction the claim wanted — but the *minimum*
over phases is +0.051, and that minimum falls on the all-trans, the only phase that can carry the
polarization the application needs, and it is **not resolved**. The claim "no antipolar escape
route" rested on the minimum being large. It is not large; it is unmeasurable. Withdrawn.

It is also no longer a distinguishing property. On the corrected measurement the minimum over
phases is +0.029 for PVDF, +0.046 for AN, +0.051 for trfe-cand, +0.134 for vclcn-cand, +0.205 for
pvf-cand, +0.247 for CFE, +0.251 for CDFE and +0.265 for vfcn-cand — **eight of the ten chemistries
screened have their minimum inside the error bar**, spanning 0.24 in total. Only PVDC (+0.397,
polar) and VDCN (−0.645, antipolar) have a minimum the model resolves. "Largest minimum" ranks
eight numbers that are all the same number as far as this potential can tell.

**Stands: trfe-cand's all-trans sits closer to its own lattice ground state than anything else
screened.** That number never used the antipolar branch. Its TT is **+0.424 kcal/mol per monomer**
above its own best packed phase, against PVDF's +0.747, CFE's +12.98, CDFE's +10.25, AN's +5.11,
VDCN's +7.55, pvf's +0.984, vfcn's +4.07 — the smallest in the screen, and outside the error bar
by a factor of 1.6, so it is resolved. It is the one design criterion of the original three that
is both resolved by this model and satisfied by this candidate, and it is the criterion that
matters for a ferroelectric: the polar phase has to be the phase the polymer is in.

**And the conformational half of its accessibility moved against it.** With the angle-relaxed
scan trfe-cand's ΔE(TT) above its conformational ground state goes from +3.96 to +5.41 kcal/mol
per monomer, so on that column it is no longer ahead of PVDF (+4.58, rigid scan): the ordering
by conformational accessibility is now pvf-cand +4.06, VDCN +4.10, vclcn-cand +4.13, PVDF +4.58,
trfe-cand +5.41, AN +5.91, vfcn-cand +8.62, CDFE +17.71, CFE +21.22, where it was VDCN 0.00,
trfe-cand 3.96, pvf-cand 4.06, PVDF 4.58, AN 6.14, vfcn-cand 6.71, vclcn-cand 9.89, CFE 11.16,
CDFE 16.67. The lattice half, +0.424 against PVDF's +0.747, is what the claim above rests on and
it was not re-measured; read together, trfe-cand's lead is now a lattice-level statement only,
and a column that mixes a rigid scan (PVDF, pvf-cand) with a relaxed one (everything else) is
being compared at a resolution of a kcal/mol or two at best (the basin-minimum caveat in step 2).

So the honest statement about trfe-cand is **weaker than the one withdrawn and points the same
way**: the model says its polar zigzag is unusually accessible *at the lattice level*, and says
nothing trustworthy about whether that zigzag packs polar. That is still worth having — TrFE is the copolymer partner
already used to stabilise the polar phase of PVDF, and the accessibility column recovered it from
a potential fitted only to PVDF and from a rule extracted with no knowledge of it — but it is
evidence about the accessibility column, not about the polarity column, and it should be cited as
such.

### The nitrile candidates

vfcn-cand and vclcn-cand were proposed to test design rule 1 and were then reported as confirming
it. They do not confirm it and they do not refute it: their all-trans gaps are +0.265 and +0.134,
both inside the error bar, and the sign of both is *positive* — the opposite of the rule — which
is exactly why a sign inside the noise must not be read.

What does survive for them is the accessibility column, which is resolved and which rejects both:
vclcn-cand's polar zigzag sits +16.34 kcal/mol per monomer above its own lattice ground
state and vfcn-cand's +4.07, against PVDF's +0.747 and trfe's +0.424. Adding a nitrile buys dipole
and costs the crystal — which was rule 1's conclusion, reached through the column that is resolved
rather than the one that is not. (Their conformational ΔE(TT) moved with the relaxed scan, vfcn
up to +8.62 and vclcn *down* to +4.13 — its rigid trans state was a ±130 deg basin edge 20
kcal/mol below planar — but the rejection is the lattice number and that did not move.)

### And pvf-cand

Its all-trans polarity is not resolved (+0.205), its polar phase sits +0.984 above its own lattice
ground state — worse than PVDF and twice trfe's — and its response is a non-measurement, because
its `C66` relaxes to 1.04 GPa and the two routes to `d` disagree by 67 %. Nothing about it
survives as a recommendation.

## What was wrong, in order (this document's former addenda 1, 2 and 3)

One column — polar versus antipolar, and the "antipolar gap" that backs it — has been wrong four
times for four different reasons, and a design rule and a candidate recommendation went with it
each time. Two of the four were found on this re-run (rounds 3b and 3c). The episodes are
collected here so that the rest of the document can be read forwards.

**Round 1 — "fitting the potential broke β-PVDF's polarity". Withdrawn; the measurement was
wrong.** A hand check reported β packing antipolar under `pvdf-dft-valence-flux` by 0.125 kcal/mol
per monomer where the illustrative potential gave +1.73 the other way, and concluded that the fit
had broken the most application-relevant fact about the material. The check set
`phi2 = phi1 + 180` on a grid and then polished with **every cell variable free**. A polish does
not preserve that relation, so the "antipolar" start slid back into the polar basin and returned
something below the polar minimum it was being compared against. Constrained properly, β is polar
under both potentials and both summation methods. **The −0.125 is withdrawn.**

**Round 2 — "truncated electrostatics is the cause". Withdrawn; Ewald changes no verdict.** The
diagnosis attached to round 1 was that a polar/antipolar margin is decided by a dipole–dipole
lattice sum, that such sums are only conditionally convergent, and that this package's
damped-shifted-force sum at 8 Å therefore could not evaluate them at any cutoff. The first two
clauses are true and the conclusion did not follow. `polyfind.ewald` was written and validated
(rock-salt Madelung constant to every published digit, asserted without a tolerance; total
independent of the splitting parameter to ten digits; gradients against finite differences to
1e-9; a brute-force spherical sum to 100 Å recovering the vacuum surface term to four digits,
which measures the conditional convergence rather than asserting it) — and then **Ewald moved
every polarity gap by at most 0.09 kcal/mol per monomer and changed no verdict anywhere.**
Truncation was not the problem. This re-run re-confirms it on all nine all-trans rows and
adds PVDF's two helices, which had never been measured under both sums: α +0.027 / +0.029 and
γ +0.584 / +0.603 through the screen's own protocol, α +0.084 / +0.094 through the acceptance
test's. The other seventeen helical rows in the table above are Ewald only.

Two things Ewald did settle, and they are worth keeping. First, the **boundary convention matters
more than the summation method**: tinfoil (the bulk limit of a screened or short-circuited
crystal, the condition under which a ferroelectric's spontaneous polarization is defined,
measured, and computed in the Berry-phase references this package compares against) and vacuum (an
isolated unelectroded sample, paying its own depolarising field) differ by 1.17 kcal/mol per
monomer on β — enough to flip α right and β wrong. Tinfoil is the default, and no energy from that
module is quotable without its label. Second, a truncated dipole sum still cannot be *defended* in
principle; it just does not happen to be what was wrong here.

**Round 3 — the antipolar cell was not antipolar.** `fitting.antipolar_cell` builds the antipolar
branch as "chain 2 flipped, setting angles equal". Chain 1 contributes `Rz(phi1) m` to the cell
dipole and chain 2 contributes `Rz(phi2) M^flip m` with `M = diag(1, −1, −1)`, so the cell is
antipolar for every `(a, b, dz)` only when `Rz(phi2) M^flip m = −Rz(phi1) m`. Writing
`m = (m_t cos θ, m_t sin θ, m_z)`, the flip branch needs `dphi = 180 + 2θ`. Equal setting angles
means `dphi = 0`, which is that condition only at `θ = 90°`. **That holds for PVDF's α helix,
whose `m_x` is exactly zero, and fails for every planar zigzag**, whose moment lies along its own
`x`: there the flip is a *rotation* of the chain and reverses no dipole at all. On β-PVDF the
"antipolar" cell comes back carrying `|P| = 0.1155 C/m²` — the polarization of the polar minimum —
at an energy degenerate with it, a gap of +0.0000.

**So every all-trans antipolar gap this document used to print was a comparison between two polar
cells.** The correct offset, derived from each chain's own moment, is 180° for every all-trans
chain in the screen and 302° for CDFE's, and never the 0° the old construction used. The two
branches were not searched alike either: `antipolar_cell` polishes from two starts with no screen
and with `(a, b)` unbounded, while the polar branch is screened inside `default_bounds`.

**Round 3b — the corrected helper's own axial scan was a point count.** `antipolar_cell_exact`
sampled `dz` at four points across the repeat whatever the repeat was. On a one-monomer all-trans
cell that is 0.64 Å and harmless. On a **four-monomer γ-type repeat it is one sample per
monomer** — aliased exactly to the interchain registry the scan exists to search, so the grid can
only ever see one registry phase, and every γ-type row in this screen was screened that way. On the
30.8 Å repeat of the 11:1 VDF/VDCN copolymer it is 7.7 Å, and the helper returned −3.7152 kcal/mol
per monomer where a finer scan of the identical subspace returned −3.7514
(`examples/copolymer_readme.py`).

The resolution is now a length: `dz_step` (0.5 Å), `phi_step` (30°) and `ab_step`, with the `(a, b)`
grid still sized from the screen budget. The start selection is now two-pass as well — half the
starts spread over distinct cells as before, the rest over distinct *registries*, because the old
rule kept only one start per `(a, b)` region and so discarded every distinct axial registry.

**And then the measurement, which does not say what the premise said.** Two results, both
negative, and they are reported because the alternative is claiming credit for a fix that changed
no number.

*On every phase in this screen the dz resolution changes nothing.* Same subspace, same polish, the
old four-point scan against the new 0.5 Å one, truncated sum, PVDF:

| phase | repeat | old `n_dz` | gap | new `n_dz` | gap | `antipolar_cell` for comparison |
|---|---|---|---|---|---|---|
| TT (β) | 1 monomer, 2.56 Å | 4 | +1.0178 | 6 | +1.0178 | +0.0000, `\|P\|` = 0.1155 — **the defect** |
| TG+TG- (α) | 2 monomers, 4.52 Å | 4 | +0.0844 | 10 | +0.0844 | +0.0844, `\|P\|` = 0.0000 — correct |
| T3GT3G' (γ) | 4 monomers, 9.04 Å | 4 | +0.6649 | 19 | +0.6649 | +0.6649, `\|P\|` = 0.0000 — correct |

γ is the case the aliasing argument is about — four dz samples across four monomers is one sample
per monomer — and the exact polish recovers the same minimum from the aliased grid regardless. The
same table incidentally confirms the other half of round 3: `antipolar_cell` is **exactly right for
both helices** and wrong only for the zigzag, which is what its `m_x = 0` condition predicts.

*And on the copolymer, the 0.036 came from the start selection, not the dz step.* Re-run at the
polar cell's own `(a, b)`, the subspace minimum is −3.7514 at *either* dz resolution: the fine
grid's own screen minimum is −2.93 against the coarse grid's +1.01, a four-kcal/mol improvement in
the *screen*, and the polish arrives at the same cell from either start. What the old helper
actually lost on that case was one start per `(a, b)` region, so every distinct axial registry but
one was discarded before the polish ever saw it.

So the fix is justified on the aliasing argument, which is a provable blindness of the grid, and it
makes the search resolution-aware rather than repeat-length-dependent — but **it moved no number in
this document.** The claim that survives is the weaker one: any antipolar energy that helper
produced on a multi-monomer repeat before the change was a lower bound on the search rather than a
converged minimum, and on every case that has since been checked the bound was tight.

**Round 3c — and the script's antipolar branch was not even using the fitted potential.** Found on
this re-run, and it is the defect that made the published `antipolar gap` column the number it was.
`fitting.FFParameters.applied()` forces a fitted potential by **rebinding
`polyfind.pack.CrystalPacker`** for the duration of its block — so that a caller's `eps_r=1.0`
default cannot quietly override the fit. A name bound by `from polyfind.pack import CrystalPacker`
*before* that block therefore still points at the original class, and a packer built from it
carries the **illustrative** potential. `examples/screen_electroactive.py` imported the name at the
top of its packing section and built its antipolar packer from it, while its polar branch went
through `pack()`, which resolves the class through the module and so was fitted. **Every antipolar
gap that script printed was a difference between two different potentials**, on top of being a
difference between two polar cells.

It is measurable to the digit. Under the fitted potential the legacy `antipolar_cell` returns
β-PVDF's own polar minimum — gap `+0.0000`, `|P| = 0.1155 C/m²`. Build the same packer from the
unpatched class and it returns **−1.4835**, which is the **−1.48 the original table published**. So
the two defects compose exactly: the construction handed back a polar cell, and the potential
mismatch made that polar cell look 1.48 kcal/mol per monomer cheaper than the polar cell it was
being compared with. `fitting.predict` has always resolved the class through the module and carries
a comment saying why; the example script did not. It does now, and `antipolar_cell_exact`'s
docstring names the trap.

**What was withdrawn, in one table.**

| quantity, as published | what it actually was |
|---|---|
| β-PVDF antipolar by −1.48 | a polar cell under the illustrative potential against a polar cell under the fitted one; β is **polar by +1.02** |
| AN all-trans antipolar by −1.86 | the same two defects; **not resolved** (+0.046) |
| VDCN all-trans antipolar by −7.60 | the same two defects; antipolar by **−0.645**, the one surviving antipolar verdict |
| vfcn-cand −3.09, vclcn-cand −1.46, pvf-cand −0.26 | the same two defects; all three **not resolved** |
| the `arrangement` column, every row | re-measured above, with an error bar beside it |
| design rule 1, "μ⊥ > 0.5 ⇒ antipolar, four of four" | **one of four**, the other three not resolved |
| trfe-cand "+1.66 / +1.87 / +1.96, no antipolar escape route" | +0.051 / +0.406 / +0.574; the minimum is **not resolved** |
| "β packs antipolar under the fitted potential by 0.125" (round 1) | an unconstrained polish leaving the subspace |
| "truncated electrostatics is the cause" (round 2) | Ewald changes no verdict; ≤0.09 everywhere |

**And the one thing that is better than the record.** Every defect found was in the *comparison* —
the subspace, its search, and which potential each branch was scored with — and none was in the
potential or in the electrostatics. The potential's behaviour on this question was never
as bad as rounds 1 and 2 said: β polar, γ polar, α undecided by a third of the noise. What the
re-run establishes is not that the model is worse than thought but that **the question was being
asked below the model's resolution**, which is a different and more actionable problem.

## What to trust

**Stand behind:** the conformational enumeration and the ΔE column (exact given the RIS model, and
differences of 3.1 to 21.2 kcal/mol per monomer against a 0.27 error bar — with the RIS model
now fitted on an angle-relaxed scan for the bulky-pendant chemistries, after the VDCN row was
found to be a frozen-angle artefact and withdrawn; and with the caveat that the model sits 3-7
kcal/mol per monomer below its own chains for the nitriles and PVDC, 2.8 for PVDF, and 20-30 for
the chlorinated pair, whose rows are not a ranking at either level); the refined lattice
ranking where its gaps exceed 0.27; the cells and densities (every PVDF edge within 8 %, five of
nine within 2 %); the polymorph ordering α < γ < β with β 3.13 kJ/mol per monomer above α, inside the
2.6 to 6.5 kJ/mol range four independent studies across five functionals agree on; the axial
`C_33` (330 GPa for β-PVDF against a literature DFT 287, measured twice by independent routes —
stress and energy curvature — agreeing to 0.02 %); the Ewald implementation itself
(`DESIGN.md` 5.11); the nineteen polarity verdicts whose gaps exceed the error bar — of which
four are on an all-trans phase, and those four are the ones the campaign can use.

**Use only as a coarse ranking, no finer than a factor of 1.7:** `d_33`, `d_31`, `e`, the blocking
stress, the free strain and the work density. Step 1 licenses the *order*, not the values, and only
at the resolution the systematic-error spread allows — and it licenses that much only for the
charge model, not for the pipeline.

**Do not stand behind:** any polarity verdict whose gap is inside 0.27 kcal/mol per monomer, which
is five of the nine all-trans rows and nine of the twenty-eight rows overall — the table says
which; any piezoelectric value as a number; the transverse elastic constants beyond an order of
magnitude (`docs/ELECTROMECHANICS.md` section 8.1); anything for CFE or CDFE that depends on their
backbone angles, which are still PVDF's borrowed 114/114 and are the assumption PVDC's correction
already overturned once; and any comparison between chemistries closer together than a factor of
1.7 in the response columns.

**What would move this most.** The ordering has changed, because the binding constraint has.

1. **A potential that can resolve a polar/antipolar margin**, which means a held-out error well
   under 0.27 kcal/mol per monomer on *lattice* energies — not on oligomer conformers. Nothing else
   in this list matters until that exists, because the screen's headline column is currently
   "cannot tell" for the majority of the sample. Concretely it means first-principles crystal
   energies for the two arrangements of at least one chemistry, which nothing in this project has.
2. The bulk polarization-strain reference of `docs/REFERENCE_DATA_REQUEST.md`, which is the only
   thing that turns the response columns from an ordering into values. Ewald is now in place, so
   such a reference would no longer be confounded by the summation method — that was a prerequisite
   and it is met.
3. Relaxed backbone angles, without which CFE, CDFE and FANOME are scored against reference states
   they cannot occupy.
4. A shape parametrisation that does not require an exact line-group pattern, without which no
   deflected chain — PVDC's real structure among them — gets a response at all.

**What would not move it.** More search and more speed. The antipolar search is converged to about
1e-3 kcal/mol per monomer, three orders of magnitude inside the error bar, and the whole screen
runs in an hour on six cores. Neither is the limit on anything stated here. Two *protocol* choices
differ by 0.06 on α, which is larger than the Ewald correction and still a fifth of the error bar —
the same conclusion from the other side.

## Reproducing this

```
set PYTHONPATH=src
set OMP_NUM_THREADS=1
set POLYFIND_TABLE_PROCS=1          rem one table worker per screen process, not per core

rem the polarity column, under Ewald, one chemistry per process
python examples/screen_electroactive.py --screen     --coulomb ewald --only pvdf      --json pvdf.json
python examples/screen_electroactive.py --candidates --coulomb ewald --only trfe-cand --json trfe.json
rem ...and the same for pvdc, cfe, cdfe, an, vdcn and pvf-cand, vfcn-cand, vclcn-cand

rem the truncated twin of the same measurement, for the Ewald attribution
python examples/screen_electroactive.py --screen --coulomb dsf --only pvdf --json pvdf-dsf.json

rem everything in one process, truncated, as the first version of this document ran it
python examples/screen_electroactive.py --all --json out.json
```

The acceptance-test 2x2 comes from the package rather than from a script:

```python
from polyfind.ewald import EwaldSpec
from polyfind.fitting import FITTED_VALENCE, acceptance_tests, acceptance_table
for coulomb, anti in (("dsf", "legacy"), ("dsf", "exact"), ("ewald", "legacy"), ("ewald", "exact")):
    spec = None if coulomb == "dsf" else EwaldSpec()
    print(acceptance_table(acceptance_tests(FITTED_VALENCE, coulomb=coulomb, ewald=spec,
                                            antipolar=anti)))
```

`acceptance_tests()` with no arguments is unchanged and still gives the recorded score; the three
switches exist so that the attribution above is reproducible, and `beta_polar_gap` comes back as the
control that separates the two constructions.

Which chemistries get the angle-relaxed RIS scan, and why, is reproduced by

```
python examples/angle_relaxation_defaults.py --candidates
```

(3-5 s per chemistry: a rigid and a relaxed one-bond scan per bond type, and the rigid minima
with no relaxed minimum within 30 deg). The regenerated conformational rows come from the screen
itself — `fit_ris`'s default `angles="auto"` now reads `forcefield.ANGLE_RELAXATION_DEFAULTS` for
a registered polymer and measures a candidate on the fly — or, for a side-by-side of the two
scans on one chemistry, `fit_ris(..., angles="rigid")` against `angles="relaxed"`.

Do **not** set `POLYFIND_TABLE_CACHE=1`. That variable is read as a directory *name*, so the value
`1` makes a directory called `1`, and 85 MB of interaction tables were once committed because of
it. The screen passes `table_cache_dir=None` throughout and pays the build every time,
deliberately: a disk table is not keyed on the potential.

`--coulomb ewald` puts both branches of the polar/antipolar comparison on `polyfind.ewald` —
screened with the truncated twin, polished with Ewald, which is the division of labour
`pack(coulomb="ewald")` already draws because the reciprocal half of the sum is not pairwise and
cannot be tabulated. The refined cells, the densities, the lattice ranking and the response columns
stay on the truncated sum, which is the protocol the first version of this table was measured with,
so that one thing changes at a time. The two acceptance tests that depend on a crystal are measured
under both sums separately. `--error-bar` overrides the 0.27 threshold; `--anti-target`,
`--anti-polish` and `--anti-maxfev` size the antipolar search.

487 tests pass, 5 skipped (three added: that the axial resolution is a length, that the
subspace checks its own polarization, and that `predict`'s new switches default to the fit's own
choices bit for bit).
