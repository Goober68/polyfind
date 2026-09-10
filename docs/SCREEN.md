# The first real screen: can this package rank electroactive polymers?

The campaign wants polymers that convert field energy into mechanical work: PVDF plus CFE,
CDFE, AN, VDCN, VDC, CNEPO and FANOME. This is the first time the whole funnel has been
pointed at that question rather than at PVDF alone.

Everything here is reproduced by `python examples/screen_electroactive.py --all`, which is
where the numbers come from and which prints more of them than this file quotes. Potential:
**`pvdf-dft-valence`** everywhere (RIS fit, packing, refinement, deformable mechanics) with
**`pvdf-dft-valence-flux`** added to the mechanics packer only, the same split
`examples/electromechanics.py` uses. It is the best available fitted preset because it is the
only one with valence terms, hence the only one on which `C_33` and the diagonal columns of
`e` exist at all. It is fitted to **PVDF** single-chain PBE-D3 energies and forces, so every
number for a chlorine or nitrile chemistry is a transfer with nothing behind it but the
wildcard bond and angle types of `fitting.VAL_BOND_TYPES`.

## Answer first

**Yes for the structural half, and it is decisive.** The screen separates the six chemistries
it could screen cleanly on the only criterion that cannot be argued away — whether a polar
chain conformation is accessible and whether it packs polar — and the separation does not
depend on any of the weak numbers. Of the campaign's eight, one (CNEPO) is out of scope and
one (FANOME) could not be screened, both for reasons stated below.

**Yes, with a caveat, for the response ranking.** The ordering of the dipole-strain response
survives leave-one-out at Spearman `rho = +1.000`, and the scale error is systematic to within
a factor of 1.7. But the ranking is carried by the *fixed* charges, not by the fitted flux, so
what the leave-one-out actually licenses is coarser than it looks: see below.

**No for the response figures as a ranking of usefulness.** Measured across the eight
chemistries that produced one, work density is *anti*-correlated with the transverse chain dipole
(Spearman −0.76) and *positively* correlated with how far the polar phase sits above the
ground state (Pearson +0.63). Ranking on work density alone picks the chemistries whose polar
phase the polymer would never adopt. The two halves have to be read together.

---

## Step 1. Is the response ranking meaningful?

`docs/BENCHMARK.md` records the piezoelectric magnitudes as 5x and 9x short, from a
charge-flux coefficient fitted to twelve numbers at `R^2 = 0.39`. A screen survives that only
if the error is **systematic**, and the sibling project's `field_neighborhood_refined` is the
one place that is testable: four chemistries with an axial dipole response at two strains.
The shipped coefficient was fitted to all four, so comparing against them is circular. Refit
with each chemistry held out, predict the held-out one, and compare orderings.

Reference `dmu/d(axial strain)`, e.A per unit strain (exploratory GFN2-xTB, finite pinned
two-chain pair in vacuum — **not** a bulk crystal):

| system | reference vector | \|ref\| | LOO prediction, shipped fit | pred/ref |
|---|---|---|---|---|
| AN | (+13.195, −18.645, −24.372) | 33.40 | 24.25 | 0.73 |
| VDCN | (+9.472, +10.804, +12.159) | 18.82 | 17.92 | 0.95 |
| CNEPO | (+3.485, −0.290, −4.282) | 5.53 | 3.26 | 0.59 |
| PVDF | (−0.106, +1.909, −4.143) | 4.56 | 2.49 | 0.55 |

**Verdict: the ranking survives.** `AN > VDCN > CNEPO > PVDF` is reproduced exactly under
leave-one-out, `rho = +1.000`, component-wise Pearson `r = +0.966`. It survives for the two
alternative fit families too (`rho = +1.000` for angle-only-four-pairs and for angle+bond),
even though the held-out coefficients swing wildly — `k_angle(C-H)` runs over −1.62 .. −0.15
across the four holdouts, a factor of eleven, and the ordering does not notice.

**The error is systematic, and here is its size.** Predicted/reference spans 0.55 .. 0.95: the
model under-predicts *every* chemistry, by between 5 % and 45 %. Take out that common
under-prediction and what is left — the part that is not systematic — is a factor of **1.74**
between the best- and worst-predicted chemistry. **Two chemistries whose true responses differ
by less than about 1.7x can therefore change places, and anything further apart than that will
not.** The reference set spans 7.3x from PVDF to AN, which is why its ordering is never in
doubt, and it is also why that ordering on its own is a weak test. This spread, not the
correlation, is the number to quote when someone asks what the screen can resolve.

**And the correlation is weaker evidence than it looks, for two separate reasons.**

*First, n = 4.* A random ordering of four systems gives `rho = 1` with probability 0.042 and
`rho >= 0.8` with probability 0.167. One perfect ordering is a one-in-twenty-four result, not
a demonstration.

*Second, and more important: the flux is not what does the work.* Run the same test with the
**fixed-increment** charges and no flux at all — a charge model with nothing whatever fitted
to this data — and it also gives `rho = +1.000`, with a *better* component-wise correlation
(`r = +0.989` against the flux fit's +0.966) and a *tighter* scale error (0.54 .. 0.78, a
factor of 1.44). The ranking is a property of the chemistries' geometry and of bond-charge
increments fitted to PBE-D3 energies; the flux changes magnitudes by up to 22 % and never
changes the order.

That is a better result than the one that was asked for, and it relocates the trust. The
badly determined `R^2 = 0.39` coefficient is what makes beta-PVDF's `d_31` non-zero and
positive at all, and that remains as poorly determined as `docs/BENCHMARK.md` says. It is
**not** what ranks chemistries. So the step-2 response columns are usable as a coarse ranking
and remain unusable as values, exactly as before.

**One limitation of this test that no amount of resampling fixes.** It is run on the
*reference's own geometries*: the question it answers is "does a charge model calibrated
without chemistry X reproduce chemistry X's dipole response at a geometry someone else
supplied?" It says nothing about whether polyfind's own predicted structure for X is the right
one, and the step-2 response columns depend on that as well. So the leave-one-out validates the
charge model's transferability across chemistries and **not** the pipeline end to end. The
structural half of the screen has its own, independent validation — the PVDF polymorph
ordering, cells and densities below — and the two should be read as separate pieces of evidence
rather than as one.

## Step 2. The screen

Seven chemistries are expressible. **CNEPO is out of scope**: its epoxide bridges two backbone
carbons, which is a change of chain topology rather than a substituent, and
`polymers.BackboneAtom` carries pendants, not rings (`docs/CHEMISTRY_EXTENSION.md` phase 4).
It is reported as out of scope rather than approximated by something else.

Method per chemistry: `fit_ris` (step 10 deg, 6 monomers, third order, ~2900 evaluations),
exhaustive `enumerate_periodic` to period 8, then packing of the top three by RIS energy
**plus** all-trans, TG+TG- and T3GT3G' whatever they rank — because screening only the top of
the RIS list lets a polymer look useless when its polar phase happens to rank 63rd, which for
PVDF it does. Packing follows `fitting.predict`'s protocol, the one the alpha-below-beta
validation was measured with: exhaustive table screen, exact polish, the symmetric antipolar
branch searched explicitly, then refinement with the bond-angle strain added back.

### The conformational screen — the column that decides

| polymer | conformations | ground state | all-trans rank | ΔE(all-trans) | chain μ⊥ per monomer |
|---|---|---|---|---|---|
| VDCN | 83 | **TT** (−11.48) | **1 of 83** | **+0.00** | 0.645 e.A |
| PVDC | 81 | TTG+TTG+ (−8.80) | not enumerated | +2.85 | — |
| PVDF | 84 | TG+TG+G+TG+G+ (−4.65) | 63 of 84 | +4.58 | 0.362 e.A |
| AN | 164 | TG− (−13.05) | 105 of 164 | +6.14 | 0.561 e.A |
| CFE | 157 | TG− (−30.90) | 41 of 157 | +11.16 | 0.386 e.A |
| CDFE | 162 | TG+G+G+ (−30.42) | 95 of 162 | +16.67 | 0.325 e.A |
| FANOME | — | not screened | — | — | — |

ΔE is kcal/mol per monomer above the conformational ground state; it is a *difference*, so the
known badness of the all-trans RIS reference cancels out of it even where the absolute
energies are meaningless (CFE's and CDFE's −30 are relative to a state carrying 162 and 199
kcal/mol of all-trans strain). μ⊥ is the **all-trans** repeat's transverse dipole per monomer,
computed with the fitted bond-charge increments rather than the polymer entry's illustrative
charges; it is the component that survives in a planar zigzag and cancels in a screw helix,
so it is the one β-PVDF's ferroelectricity is made of. PVDC has no entry because it has no
all-trans repeat to measure one on.

Two entries need their own sentence.

**PVDC has no polar conformation at all, and the reason is geometric rather than energetic.**
Its measured backbone angles are unequal — 123 deg at CH2, 114 deg at CCl2, both
crystallographic — so an ideal-angle repeat carries 123 − 114 = 9 deg of curl and every
sequence is a circular arc with zero rise. `enumerate_periodic` drops all-trans and TG+TG-
both, on the rise-per-bond filter. The screen rebuilds them anyway and searches for the
torsion deflection that closes them: **all-trans still does not close, TG+TG- does**, and the
deflected glide is the only thing PVDC packs. So PVDC has a crystal and it is not a polar one.
**PVDC is not a candidate**, and the model says so for the same reason the crystallography
does — the real polymer is a TGTG' glide, not a planar zigzag.

The closure is one equation in four unknowns, so its solutions are a manifold and which point
of it you land on depends on where you start. Started from the *fitted* model's own adapted
state angles the screen reaches 177.5 / 27.7 deg and, after refinement, c = 4.76 A against a
measured 4.68. Started from the ideal 180 / 60 it reaches 174.7 / 54.7 with c = 4.70, and
`tests/test_pvdc.py` reaches 175.3 / 49.4 with c = 4.677 — that last one being the published
175 / 49 and 4.68. The *fibre repeat* is robust across all three to about 2 %; the torsion
pair is not, and only the ideal-angle start recovers it.

**FANOME could not be screened at all.** Its all-trans reference is an overlapping structure —
methyl hydrogens of methoxy groups on consecutive substituted carbons 0.80 A apart, ~1e6
kcal/mol on a ten-bond oligomer, robust to the frozen C-O rotamer — so `fit_ris` would measure
every RIS energy from a state that is not a molecule. It is registered, its geometry is
tested, and it is deliberately not fitted (`docs/CHEMISTRY_EXTENSION.md` phase 3). Reporting a
number for it would be reporting noise.

### Packing: cells, densities, and polar or antipolar

`arrangement` is judged on the **rigid** comparison, which is the one that means something:
`antipolar_cell` measures the symmetric antipolar subspace explicitly, while a refinement
started there slides back out of it. The gap is E(antipolar) − E(polar) per monomer, so
negative means the crystal wants to be antipolar and the net polarization is zero.

> **WITHDRAWN — the `arrangement` and `antipolar gap` columns of this table and of the
> candidate table below are not what they say.** `antipolar_cell`'s cell is not antipolar
> for a planar zigzag, so every all-trans row compares two polar cells. See addendum 2 at
> the end of this document for the corrected measurement.

| polymer | phase | a × b × c (A) | ρ | arrangement | antipolar gap | ΔE lattice |
|---|---|---|---|---|---|---|
| PVDF | TG+TG- (α) | 4.99 × 9.04 × 4.71 | 2.000 | polar | +0.09 | +0.000 |
| PVDF | T3GT3G' (γ) | 5.28 × 8.89 × 9.34 | 1.939 | polar | +1.20 | +0.470 |
| PVDF | TT (β) | 4.54 × 8.55 × 2.60 | 2.107 | antipolar | −1.48 | +0.747 |
| PVDC | TG+TG- (deflected) | 6.23 × 11.23 × 4.76 | 1.932 | polar | +3.46 | +0.000 |
| CFE | TG+TG- | 5.49 × 10.69 × 4.74 | 1.924 | polar | +2.03 | +0.000 |
| CFE | TT | 5.49 × 8.25 × 2.69 | 2.190 | polar | +0.67 | +12.981 |
| CDFE | TG+TG- | 5.53 × 10.95 × 5.03 | 2.149 | polar | +1.72 | +0.000 |
| CDFE | TT | 6.15 × 8.29 × 2.69 | 2.381 | polar | +1.50 | +10.252 |
| AN | TG+TG- | 6.28 × 9.74 × 4.73 | 1.220 | antipolar | −1.32 | +0.000 |
| AN | TT | 11.27 × 4.20 × 2.66 | 1.403 | antipolar | −1.86 | +6.974 |
| VDCN | T3GT3G' | 9.55 × 7.76 × 9.33 | 1.500 | antipolar | −2.84 | +0.000 |
| VDCN | TT | 6.06 × 10.57 × 2.69 | 1.504 | antipolar | **−7.60** | +7.400 |

**The PVDF rows are the calibration and they behave.** The polymorph ordering comes out
α < γ < β, with `E(α) − E(β) = −0.747 kcal/mol = ` **−3.13 kJ/mol per monomer** — inside the
−6.5 .. −2.6 kJ/mol range that four independent studies across five functionals agree on, and
consistent with the −2.96 recorded for this preset in `fitting.py`. That is acceptance test 1
passing, reproduced here through a different code path than the one that recorded it.
α's antipolar gap is +0.09, the same near-miss on acceptance test 3 that is
already on record. Cells against experiment (4.91 × 8.58 × 2.56, 4.96 × 9.64 × 4.62,
4.96 × 9.67 × 9.20): every one of the nine edges within 8 %, and five of the nine within 2 %.
The worst are β's `a` (4.54 against 4.91, 7.5 %) and γ's `b` (8.89 against 9.67, 8.1 %).

**VDCN is the sharpest result in the table.** Its all-trans is the conformational ground
state — the only chemistry here for which that is true — and it carries the largest transverse
chain dipole of anything screened, 0.645 e.A per monomer against PVDF's 0.362. And its
crystal wants to be antipolar by **7.60 kcal/mol per monomer**, five times PVDF's β and the
largest antipolar preference measured anywhere in this screen. The conformation is free; the
polarization is not.

### The response figures — the weak half, flagged

Computed on the **polar** structure whether or not the polymer would adopt it, because the
question is what the polar phase would deliver. `d_33` and `d_31` are in the film convention
(3 = poling = the crystal's own polar axis, 1 = draw = the chain axis), generalised from
`fit_charge_flux.py`'s hard-coded x-axis to the actual polarization direction; for beta-PVDF
that reduces to the same numbers, which is the check.

| polymer | phase | \|P\| C/m² | C11 | C22 | C66 | C33 | d33 | d31 | blocking GPa | free strain | work kJ/m³ | routes | status |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| PVDF | TT | 0.150 | 28.4 | 20.1 | 5.1 | 330 | −1.88 | +2.33 | 0.072 | 2.3e−4 | 4.22 | 0.33 % | ok |
| CFE | TT | 0.168 | 32.9 | 27.6 | 6.1 | 636 | −1.20 | +1.36 | 0.021 | 1.8e−3 | 5.57 | 2.7 % | ok |
| CDFE | TT | 0.116 | 19.1 | 14.2 | 6.0 | 614 | −0.40 | +0.96 | 0.008 | 1.7e−3 | 6.64 | 1.1 % | ok |
| VDCN | TT | 0.142 | 33.4 | 39.4 | 15.3 | 447 | −2.46 | +1.04 | 0.048 | 2.5e−4 | 2.17 | 1.3 % | ok |
| AN | TT | 0.023 | 41.0 | 26.0 | 10.6 | 517 | +0.85 | +0.19 | 0.002 | 2.1e−4 | 0.12 | 0.03 % | ok |
| PVDC | TG+TG- | — | — | — | — | — | — | — | — | — | — | — | **failed** |

Elastic constants in GPa. "routes" is the disagreement between the two independent routes to
`d` — a dipole derivative and a zero-stress root find — which is the internal check that the
flux entered the energy gradient and the dipole derivative consistently.

PVDF's row reproduces `docs/BENCHMARK.md` exactly (`d_33 = −1.88`, `d_31 = +2.33` proper,
`C_33 = 330` against a literature DFT 287), which is what says the generalised film convention
and the packing protocol are the same calculation the benchmark ran.

**PVDC's response could not be computed, in three documented steps.** No conformation closes
at ideal angles, so its chain has deflected torsions; a deflected chain breaks the glide
pattern its sequence is named for, so it has no line group and no `Shape`, so the deformable
path is unavailable; and the rigid fallback then fails its own strained-cell construction.
Two of those three are model limitations already on the roadmap (variable backbone angles, and
a shape parametrisation that does not need an exact symmetry pattern). It is reported as a
failure rather than forced.

**AN's row is not a contradiction, it is the point.** AN's TT reference re-packs to the
antipolar arrangement its own antipolar gap predicts (−1.86), so its polarization collapses
from the 0.122 C/m² of the polar branch to 0.023, and the work density with it. The
electromechanics is telling you what the packing already said.

### Wall time

| polymer | fit | enumerate | pack + refine | response | total |
|---|---|---|---|---|---|
| PVDF | 0.5 | 2.9 | 40.1 | 2.8 | **46.4 s** |
| PVDC | 0.4 | 3.7 | 20.5 | 17.9 | **42.5 s** |
| CFE | 0.6 | 5.8 | 46.5 | 4.4 | **57.3 s** |
| CDFE | 0.5 | 3.6 | 49.7 | 5.2 | **59.1 s** |
| AN | 0.5 | 4.3 | 71.6 | 5.7 | **82.1 s** |
| VDCN | 0.7 | 4.3 | 66.1 | 6.2 | **77.3 s** |
| FANOME | — | — | — | — | not screened |

Six chemistries, five to six conformations packed and refined each, one full response tensor
each: **six minutes for the whole incumbent set**, 637 s including the four new candidates and
the ranking audit. Seconds, not core-hours: six physical cores, `OMP_NUM_THREADS=1`, cold
interaction-table cache throughout (`table_cache_dir=None`, deliberately — a disk table is not
keyed on the potential). The fit and the enumeration are 1 to 6 s; **the packing is 80 % to
90 % of the cost**, and it is the part that grows with the crystallographic repeat: AN's and
VDCN's 16-bond gamma-type repeats are why they cost 80 s and PVDF 46 s.

## Step 3. What drives the answer

Correlations of each structural feature against the predicted work density, over the eight
chemistries that produced a usable response (six incumbents plus four candidates, minus PVDC
and pvf-cand which did not). **At n = 8 these describe the sample; they do not establish
anything.** They are reported because the *signs* are informative and two of them are the
opposite of what the campaign's intuition would be.

| feature | Pearson | Spearman | reading |
|---|---|---|---|
| packing density ρ | **+0.81** | **+0.91** | the strongest, and the most likely to be a confound |
| transverse chain dipole μ⊥ | **−0.69** | **−0.76** | **backwards from the intuition** |
| ΔE(all-trans above ground) | +0.63 | +0.43 | the response column rewards inaccessible phases |
| axial stiffness C33 | +0.35 | +0.19 | no rule |
| polarization \|P\| | +0.45 | +0.26 | no rule |
| transverse stiffness C11 | −0.18 | −0.29 | no rule |
| gap to the nearest polymorph | +0.16 | −0.10 | no rule |

**Rule 1, and it is the one to act on: a bigger pendant dipole does not buy more work, it buys
an antipolar crystal.** μ⊥ is anti-correlated with work density, and the mechanism is visible
directly rather than inferred. **Every chemistry in this screen whose all-trans carries μ⊥
above 0.5 e.A per monomer packs that all-trans antipolar, without exception** — VDCN 0.645,
AN 0.561, vclcn 0.557, vfcn 0.548, with antipolar gaps of −7.60, −1.86, −1.46 and −3.09. A
pendant big enough to carry a large dipole is big enough that the antipolar arrangement wins,
and the cell polarization goes to zero. The dipole you can build into a monomer is not the
dipole you get out of a crystal.

**The rule is one-sided and should be used that way.** A large μ⊥ is sufficient for an
antipolar all-trans; a small one is not sufficient for a polar one. Below 0.4 the all-trans
packing goes both ways — polar for CFE (0.386), CDFE (0.325) and trfe (0.305), antipolar for
PVDF (0.362) and pvf (0.303). So μ⊥ above 0.5 is usable as a **reject** filter and nothing
below it is usable as an accept filter; the antipolar gap has to be computed either way, and
computing it costs about a second on top of a packing that already ran.

**Rule 2: the response column ranks the wrong thing on its own.** ΔE(all-trans) correlates
*positively* with work density (+0.63): CDFE and CFE top the work-density table with 6.6 and
5.6 kJ/m³ while their polar phases sit 10.3 and 13.0 kcal/mol per monomer above their own
lattice ground states. Those are figures for a structure the polymer will not adopt. A screen
that ranks on work density alone inverts the answer; the accessibility column has to gate it.

**Rule 3, weaker: denser packing helps.** ρ is the strongest single correlate (+0.81 / +0.91)
and there is a plausible mechanism — work density is dipole per unit volume worked against
stiffness, and denser packing raises both. But in this sample density is nearly collinear with
halogen content, so it may be reading "has chlorine" rather than "is dense". Not to be used
alone.

**What the data does not support.** No rule from the axial stiffness (+0.35 / +0.19), the
transverse stiffness (−0.18), the polarization magnitude (+0.45 / +0.26) or the energy gap to
the nearest competing polymorph (+0.16 / −0.10). Stating a chain-stiffness rule from these
numbers would be inventing one. The gap result is the most disappointing, because a large gap
to the nearest competitor is exactly what would make a polar phase *stable*, and this sample
cannot see it.

**One methodological warning the screen produced by itself.** The work-density figure is
dominated by the softest reachable mode. The pvf candidate's C66 relaxes to 1.04 GPa, its
shear `d` blows up to 192 pC/N and its work density to 96 kJ/m³ — twenty times anything else —
and the two independent routes to `d` then disagree by 67 %, which is what caught it. It is
reported as a non-measurement. Any work density quoted from this package should be read
together with its softest constant and its route agreement.

## New candidates

Only run because step 1 came out positive. **These are model suggestions, not predictions.**
They are ordinary `Polymer` values constructed inside `examples/screen_electroactive.py` and
passed straight to the funnel — nothing in the package changes and no default moves — with
geometry and charges chosen exactly the way PVDF's, PVDC's, CFE's and CDFE's were. The charges
are illustrative and the potential is PVDF's.

The rules above suggest: keep the pendant small (rule 1), keep the polar conformation near the
ground state (rule 2), keep the packing dense (rule 3). That points at fluorine rather than
nitrile, which is a testable direction because it predicts the *known* answer.

| candidate | formula | ΔE(TT) above ground | TT above lattice ground | TT arrangement | antipolar gap | ρ | \|P\| | d31 | work |
|---|---|---|---|---|---|---|---|---|---|
| **trfe-cand** | −(CHF−CF2)− | +3.96 | **+0.42** | **polar** | **+1.66** | 2.310 | 0.108 | +1.69 | 2.49 |
| pvf-cand | −(CH2−CHF)− | +4.06 | +1.09 | antipolar | −0.26 | 1.689 | 0.114 | — | not a measurement |
| vfcn-cand | −(CH2−C(F)(CN))− | +6.71 | +5.05 | antipolar | −3.09 | 1.744 | 0.166 | +1.52 | 2.28 |
| vclcn-cand | −(CH2−C(Cl)(CN))− | +9.89 | +17.85 | antipolar | −1.46 | 1.737 | 0.157 | +1.08 | 1.75 |
| *PVDF, for comparison* | −(CH2−CF2)− | +4.58 | +0.75 | antipolar | −1.48 | 2.107 | 0.150 | +2.33 | 4.22 |

**None of them beats PVDF on the response figures, and one of them beats it on the criterion
that matters.** trfe-cand — trifluoroethylene, the CHF-CF2 repeat — is the only chemistry
anywhere in this screen whose all-trans phase is within half a kcal/mol per monomer of its own
lattice ground state (+0.42 against PVDF's +0.75) **and** whose every packed phase prefers
polar over antipolar by a wide margin: +1.66, +1.87 and +1.96 kcal/mol per monomer for its TT,
TG+TG- and T3GT3G'. That *minimum* of +1.66 is the largest of any chemistry in this screen —
CFE's is +0.67, CDFE's +1.03, and PVDF's is −1.48 because its β wants to be antipolar. So trfe
is the one chemistry here with no antipolar escape route. Its coefficients are
lower than PVDF's — `d_31 = +1.69` against +2.33, work 2.49 against 4.22 kJ/m³ — so the model
does not claim it is a better piezoelectric. It claims it is a better *ferroelectric*: the
polar phase is the one it would actually be in.

That is worth stating plainly because it is a positive control rather than a discovery. TrFE
is the copolymer partner already used to stabilise the polar phase of PVDF, and the screen
recovered exactly that, from a potential fitted only to PVDF and from rules extracted with no
knowledge of it. It is the strongest evidence in this document that the structural half of the
screen is doing real work.

The two nitrile candidates are rule 1 confirming itself: vfcn and vclcn carry the biggest
pendant dipoles among the candidates (0.548 and 0.557 e.A per monomer) and both pack their
all-trans antipolar, with vclcn's polar phase 17.8 kcal/mol per monomer above its own lattice
ground state — the worst accessibility in the whole screen. Adding a nitrile buys dipole and
loses the crystal. (vclcn's *helical* TG+TG- does pack polar, at +1.87; it is only the polar
zigzag that the nitrile pushes antipolar, which is exactly the phase the application needs.)

## What to trust

**Stand behind:** the conformational enumeration and the ΔE column (exact given the RIS
model); the polar/antipolar verdict and the antipolar gap (an explicit search of the symmetric
subspace, on the same protocol as the recorded PVDF validation); the cells and densities
(every PVDF edge within 8 %, five of nine within 2 %); the polymorph ordering (α < γ < β, β 3.13
kJ/mol above α, inside the 2.6 to 6.5 kJ/mol literature range); the axial `C_33` (330 GPa for
β-PVDF against a literature DFT 287, and measured twice by independent routes — stress and
energy curvature — agreeing to 0.02 %); the wall times.

**Use only as a coarse ranking, no finer than a factor of 1.7:** `d_33`, `d_31`, `e`, the
blocking stress, the free strain and the work density. Step 1 licenses the *order*, not the
values, and only at the resolution the systematic-error spread allows — and it licenses that
much only for the charge model, not for the pipeline (see the caveat at the end of step 1).

**Do not stand behind:** any piezoelectric value as a number; the transverse elastic constants
beyond an order of magnitude (`docs/ELECTROMECHANICS.md` section 8.1); anything for CFE or
CDFE that depends on their backbone angles, which are still PVDF's borrowed 114/114 and are
the assumption PVDC's correction already overturned once; and any comparison between
chemistries closer together than a factor of 1.7 in the response columns.

**What would move this most.** Not more search and not more speed. In order: the bulk
polarization-strain reference of `docs/REFERENCE_DATA_REQUEST.md`, which is the only thing
that turns the response columns from an ordering into values; relaxed backbone angles, without
which CFE, CDFE and FANOME are being scored against reference states they cannot occupy; and a
shape parametrisation that does not require an exact line-group pattern, without which no
deflected chain — PVDC's real structure among them — gets a response at all.


## Addendum: the polar/antipolar column is not trustworthy, and design rule 1 may be an artifact

Checked independently after this screen was produced, because the screen's
verdict that beta-PVDF packs antipolar contradicts the fact that beta *is* the
ferroelectric phase. Same test, both potentials, unconstrained best against the
best antipolar arrangement:

| potential | best \|P\| (C/m^2) | E(antipolar) - E(polar), kcal/mol per monomer | verdict |
|---|---|---|---|
| illustrative default | 0.1416 | **+1.73** | polar, correct |
| `pvdf-dft-valence-flux` | 0.1155 | **-0.125** | antipolar, **wrong** |

So fitting the potential flipped it. That fit improved held-out energies from
3.32 to 1.36 kcal/mol, corrected the alpha/beta ordering, and gave d31 its
experimentally correct sign - and it broke the most application-relevant
qualitative fact about the material.

**The right conclusion is not that the default potential is better.** The margin
is 0.125 kcal/mol per monomer against a 1.73 the other way, so this is a
near-degenerate balance being tipped, and the physics that decides it is a
long-range dipole-dipole lattice sum. Those decay as 1/r^3, which is
conditionally convergent: a truncated sum cannot evaluate them, whatever the
cutoff, and this package uses a damped-shifted-force sum at 8 A rather than
Ewald. DESIGN.md section 6 has listed that as a limitation since the beginning;
this is the first measurement of what it costs. The default getting beta right is
luck, not skill.

**Two consequences for how this screen should be read.**

The `arrangement` and `antipolar gap` columns are **not resolved by the model**
and should not be used. That is a change from what the screen's own summary says
it would stand behind.

More seriously, **design rule 1 is suspect for the same reason**. It says a large
pendant dipole forces antipolar packing, with no exceptions across eight
chemistries, and it is the rule the screen leans on hardest as a reject filter.
But chemistries with large dipoles are exactly where a truncated electrostatic
sum errs most, so the rule may be reporting the truncation rather than the
physics. It should be re-tested with Ewald before being used to reject anything.

Design rule 2 - that ranking by work density alone inverts the answer, because it
correlates positively with how inaccessible the polar phase is - does not depend
on the polar/antipolar determination and stands.

**What would settle it:** proper Ewald summation, which has been on the roadmap
as a known gap and now has a concrete symptom attached to it. It is also a
prerequisite for the bulk reference data requested in
`docs/REFERENCE_DATA_REQUEST.md` to be usable, since comparing a Berry-phase
polarization against a truncated-electrostatics model would confound the two.


## Addendum 2: Ewald, and why the column was wrong for a different reason

Ewald summation is now implemented and validated (`polyfind.ewald`, `DESIGN.md`
section 5.11: the rock-salt Madelung constant to every published digit, energies
independent of the splitting parameter to 1e-10, gradients against finite
differences to 1e-9). It was re-run over the same chemistries. **Three results,
and the last is the one that matters.**

### beta-PVDF is polar, under both potentials and under either sum

Addendum 1's table is superseded. Same quantity — the best exactly-antipolar cell
minus the best cell overall, rigid, two chains, gamma = 90° — with the antipolar
branch constructed from the chain's own moment rather than assumed:

| potential | sum | E(antipolar) − E(best), kcal/mol per monomer | verdict |
|---|---|---|---|
| illustrative default | truncated (8 Å DSF) | **+1.729** | polar, correct |
| illustrative default | Ewald, tinfoil | **+1.758** | polar, correct |
| illustrative default | Ewald, vacuum | −1.004 | antipolar |
| `pvdf-dft-valence-flux` | truncated (8 Å DSF) | **+1.018** | polar, correct |
| `pvdf-dft-valence-flux` | Ewald, tinfoil | **+1.021** | polar, correct |
| `pvdf-dft-valence-flux` | Ewald, vacuum | −0.769 | antipolar |

The fitted potential does **not** get beta wrong. Addendum 1's −0.125 came from
the old antipolar construction — which for a planar zigzag is not antipolar at all
(below) — together with a search of it that had not converged; the same
construction reproduces addendum 1's +1.73 for the default potential by
coincidence, because for that potential the converged value happens to agree.
Both rows of addendum 1's table are withdrawn. The one thing addendum 1 got right
is that the boundary condition matters: switch from the bulk (tinfoil) convention
to an isolated sample in vacuum and beta does flip to antipolar, because the
vacuum convention charges it 1.17 kcal/mol per monomer of depolarisation energy
for being polar. That is a statement about unelectroded slabs, not about the
crystal; `DESIGN.md` 5.11 says why tinfoil is the default.

### Ewald is not what was wrong

Three polymorphs whose experimental polarity is known, under Ewald (tinfoil), both
potentials: **beta polar +1.76 / +1.02 (correct), gamma polar +1.67 / +0.68
(correct), alpha polar +0.22 / +0.09 (wrong — it is antipolar)**. Alpha's margin
is an order of magnitude below the potential's own held-out error, so the model
does not decide alpha rather than getting it wrong with conviction, and correct
electrostatics did not move it (+0.20 / +0.08 truncated).

Across chemistries the same holds. Re-measured over the nine whose all-trans packs
at all, under `pvdf-dft-valence`, with the truncated sum and with Ewald (tinfoil),
on an identical protocol: **Ewald moves every antipolar gap by at most 0.09
kcal/mol per monomer and changes no verdict anywhere.**

| chemistry | μ⊥ (e.A/monomer) | truncated gap | Ewald gap | verdict, both |
|---|---|---|---|---|
| pvdf | 0.362 | +1.018 | +1.021 | polar |
| cfe | 0.386 | +1.524 | +1.501 | polar |
| cdfe | 0.325 | +0.410 | +0.390 | polar |
| an | 0.561 | +0.065 | +0.047 | polar |
| vdcn | 0.645 | −0.558 | −0.645 | **antipolar** |
| pvf-cand | 0.302 | +0.175 | +0.206 | polar |
| trfe-cand | 0.305 | +0.034 | +0.051 | polar |
| vfcn-cand | 0.548 | +0.282 | +0.264 | polar |
| vclcn-cand | 0.557 | +0.103 | +0.134 | polar |

So the suspicion recorded in addendum 1 — that the arrangement column was
reporting the truncation — **is not supported**. A truncated dipole sum still
cannot be defended as a matter of principle, and the boundary convention still
changes beta's lattice energy by 1.17 kcal/mol per monomer (above), but on these
structures the truncated and the correct sum agree on the sign of every gap and
nearly on its size.

### What was wrong was the antipolar cell itself

`fitting.antipolar_cell` builds the antipolar branch as "chain 2 flipped, setting
angles equal". That is the antipolar subspace only when the chain's transverse
dipole is perpendicular to the chain's own `x` axis. It is, for a TG+TG- helix,
whose `m_x` is exactly zero. **It is not, for any planar zigzag**: an all-trans
chain's moment lies along its own `x`, the flip `(x, −y, −z)` leaves that
component alone, and the "antipolar" cell comes back carrying the *full*
polarization of the polar minimum at an energy degenerate with it. Measured on
beta-PVDF: `|P| = 0.1416 C/m^2` and a gap of −0.0000.

**Every all-trans antipolar gap in this document was therefore a comparison
between two polar cells**, and the −1.48 for β-PVDF, −1.86 for AN and −7.60 for
VDCN are not antipolar gaps. That is not an inference from one case: the correct
offset, derived from each chain's own moment, is 180° for every all-trans chain in
the table above and 302° for CDFE's, and never the 0° the old construction uses.
(They are also not gaps between cells searched alike: `antipolar_cell` polishes
from two starts with `(a, b)` unbounded, while the polar branch is screened inside
`default_bounds`.) The corrected construction
is `fitting.antipolar_offsets` / `antipolar_cell_exact`, which derives the
subspace from the chain's own moment and searches it with a grid screen plus a
polish of the best distinct cells; the table above is measured with it, and the
zero dipole is checked at every answer rather than assumed (max |P| < 1e-15).

### So: does design rule 1 survive? No.

The rule was "every chemistry whose all-trans carries μ⊥ above 0.5 e.A per
monomer packs that all-trans antipolar, without exception" — four for four. Under
the corrected comparison it is **one for four** (VDCN), with the same numbers
under either sum. The Pearson correlation of μ⊥ against the antipolar gap is
−0.49 truncated, −0.51 Ewald: the *sign* of the tendency survives and the
exceptionless form does not.

**And the three that changed sides did not change to "polar" — they changed to
"not resolved".** AN +0.047, vclcn +0.134, vfcn +0.264 kcal/mol per monomer are
inside the noise of everything else in this model (the potential's own held-out
error is 1.36 kcal/mol, and the α/β polymorph gap it gets right is 0.75). Only
two chemistries here have a margin worth quoting: VDCN is antipolar by 0.6, and
PVDF and CFE are polar by 1.0 and 1.5.

**Rule 1 must not be used as a reject filter.** Not because the electrostatics
were truncated, but because the quantity it was fitted against was not measuring
what it was named after, and because the corrected quantity is below this model's
resolution for most of the sample. The `arrangement` and `antipolar gap` columns
of the tables above stay withdrawn.

Rule 2 is unaffected, as addendum 1 said. The trfe-cand result rests on the same
withdrawn column: its "+1.66 for TT, +1.87 for TG+TG-, +1.96 for T3GT3G', the
largest *minimum* of any chemistry here" was measured the old way, and re-measured
with the corrected construction its all-trans gap is +0.034 (truncated) / +0.051
(Ewald) — **not resolved**. The two helical phases were not re-measured, so the
claim that trfe has "no antipolar escape route" now rests on nothing that has been
checked, and is withdrawn with the column it came from. What survives for trfe is
the part that never used the antipolar branch: its all-trans sits +0.42 kcal/mol
per monomer above its own lattice ground state, the closest of anything screened.
