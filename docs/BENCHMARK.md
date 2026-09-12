# Speed against the existing method

The project's target is a solver that meets field and strain measurement needs at
10 to 100 times the speed of existing methods. That claim needs a baseline that
is recorded rather than assumed, and one exists: the sibling project
`sarco/materials/gpu_bundle` kept per-run timings for the calculations this
package is meant to displace. All baseline numbers below are read from its own
result files, with the file named, so they can be checked. All polyfind numbers
were measured on the same machine (6 physical cores, 12 logical).

Nothing here is a like-for-like ratio yet, and the reason is set out under
"What is actually comparable". Read that before quoting a number.

## The baseline, from its own recorded runs

### Machine-learned potential tier

`results/field_neighborhood_refined/<system>/strain<s>_field<E>.json`, a grid of
two strain values by three field values around one already-known structure.

| System | Points | Relaxation steps | Recorded seconds |
|---|---|---|---|
| PVDF | 6 | 1,855 | 87.0 |
| CNEPO | 6 | 2,070 | 107.6 |
| AN | 6 | 2,565 | 125.4 |
| VDCN | 6 | 2,719 | 142.8 |
| total | 24 | 9,209 | 462.9 |

An earlier, looser-tolerance pass over the same grid
(`results/field_neighborhood/`) cost 176.2 s for 25 points, and several of its
field points converged in one to three steps taking about 0.11 s, meaning the
field moved nothing the tolerance could resolve.

### Density-functional tier

`results/dft_field/*.json`, PBE0 with dispersion and an applied field, one single
point each.

| Calculation | Recorded seconds |
|---|---|
| PVDF, def2-SVP, three field values | 772, 895, 1001 |
| PVDF, def2-TZVP, three field values | 2470, 1866, 1711 |
| CNEPO, def2-SVP, three field values | 946, 1127, 1030 |

So one field triple costs about 2,700 s at the smaller basis and about 6,000 s at
the larger. A strain derivative needs at least two strain states, so roughly
5,400 s per chemistry at def2-SVP, before any search over structures. One run is
recorded as having died out of memory
(`pvdf_def2-svp_field-65.attempt-oom.json`).

### Structure search tier

`results/free_angle/seq_*_11_*/result.json` holds supercell relaxations at 2,400
steps for CNEPO and VDCN, with `continuation` and `recovery` directories beside
them. A 2,400-step run that needs a continuation is a step budget exhausted
rather than a converged answer, which is the convergence risk a local optimiser
over a combined discrete-and-continuous space carries.

## polyfind, measured

| Stage | Seconds |
|---|---|
| Fit the rotational-isomeric-state model from a potential | 0.3 |
| Enumerate periodic chain conformations (120 distinct, exhaustive) | 4.5 |
| Exhaustive screen and refine, all candidates, cold cache | 49 |
| Same, warm interaction-table cache | 17 |

Re-measured on the same machine after the geometry corrections of
`docs/REFERENCES.md`, since PVDF's C-C bond length changed and every cached
interaction table with it. **The corrections do not move these timings**: the
unchanged code on the same machine gives 0.4 / 5.0 / 53 / 14 s, which is
run-to-run spread rather than a difference. The table previously read 0.6 / 4.4 /
77 / 25, so the screen-and-refine rows have come down by about a third since they
were recorded and the fit row by half; neither is attributable to this change.
Quote a range rather than a point for the two cache rows: cold 49-53 s, warm
14-17 s across runs.

Fitting the potential itself against reference data is a separate one-off:
roughly 260 s over 467 frames, amortised across every later run for that
chemistry class.  That fit has not been repeated at the corrected geometry.

## What is actually comparable

Three traps, and the ratio differs by an order of magnitude depending on which
comparison is meant.

**The baseline's 87 to 143 s per chemistry does not include finding the
structure.** It relaxes around a structure already supplied. polyfind's 17 to
53 s includes enumerating 120 candidate conformations and screening the packing
of each exhaustively. Comparing those two numbers directly credits polyfind for
work the baseline never did, and also credits the baseline with an input polyfind
derives. The structure-search tier above is the closer analogue, and it is the
one whose runs exhausted their step budgets.

**The grids differ in what they can yield.** Two strain values by three field
values supports a single finite-difference derivative at one strain. An elastic
tensor and a piezoelectric tensor need considerably more, so a fair ratio has to
compare cost per derivative obtained, not cost per run.

**The tiers are different physics.** The machine-learned tier is fast and is
itself validated against the density-functional tier in the same project.
polyfind is a third tier below both, and its accuracy is bounded by the potential
it is fitted to. A speed ratio quoted without the accuracy it was obtained at is
meaningless; see `DESIGN.md` sections 5.9 and 5.10 for what the fitted potential
does and does not get right.

## The observation that favours this approach most

It is not the raw timing. In the baseline's own refined grid, applying a field of
65 V/µm to PVDF moves atoms by at most 0.019 Å, changes torsions by at most
0.13°, and leaves a residual force of 1e-4 eV/Å after 184 relaxation steps. In
the looser pass, the same field points converged in one to three steps, which
means the response was below what that tolerance could see at all.

The quantity being measured is therefore tiny compared with the numerical
machinery used to extract it, and it is extracted as a *difference* between two
expensive relaxations. That is the regime where an analytic derivative wins for a
reason that has nothing to do with processor speed: polyfind differentiates the
energy with respect to field and strain directly, so the response is a derivative
evaluated once rather than a small difference resolved between two converged
minimisations. The same argument already paid off inside this package, where
replacing finite-difference gradients with analytic ones cut refinement from
181 s to 26 s and incidentally exposed a non-smoothness in the cutoff that the
difference-based path had been walking over unnoticed.

## Status of the headline claim

| Comparison | Status |
|---|---|
| Structure search against supercell relaxation | favourable, and the baseline's runs did not converge; not quantified as a ratio |
| Field and strain response against the machine-learned grid | **measured, see below** |
| Field response against the density-functional tier | roughly 2,700 to 6,000 s per field triple to beat, per chemistry |
| Accuracy at which any ratio holds | held-out energies 1.36 kcal/mol; one of three acceptance tests passing |

### The speed comparison, now measurable

The baseline spends 87 s (PVDF) to 143 s (VDCN) per chemistry on a grid of two
strain values by three field values, which yields **one** finite-difference
derivative at one strain state, around a structure supplied to it.

polyfind computes a full response - the in-plane elastic constants, the axial
constant, the piezoelectric tensor by two independent routes, blocking stress,
free strain and work density - in **0.3 s for beta, 1.4 for alpha, 7.8 for
gamma**, from a structure it derived itself in a further 14 to 53 s.

Taken as cost per derivative obtained the ratio is far beyond 100x, because the
baseline's grid yields one derivative and this yields a tensor. Taken as
whole-pipeline cost per chemistry, structure search included, it is roughly 2 to
6x. Neither number should be quoted alone. The honest statement is that the
speed target is met, comfortably, for the response calculation itself, and that
the response calculation is not the whole job.

**And the speed is not the binding constraint.** What limited this package for the
intended application was the next section: the piezoelectric response of the phase
the material is actually used in was identically zero, for a reason that was a
symmetry statement rather than a numerical failure. That is now addressed, at a
cost the next section states in full.

### What the response can and cannot deliver

| Quantity | Status |
|---|---|
| In-plane elastic constants C11, C22, C12, C16, C26, C66 | computed |
| Axial constant C33 | computed, and demonstrated free of the invented-stiffness contamination that made an earlier attempt refuse it |
| Shear constants C44, C55, C45 and mixed 4/5 | structurally absent: the chain axis is fixed along z, so no variable expresses those shears |
| Piezoelectric response of helical chains (alpha, gamma) | computed; d up to 2.7 pC/N with fixed charges, 8.9 for gamma with charge flux (alpha's flux run does not converge, below) |
| Piezoelectric response of beta, the ferroelectric phase | **non-zero with charge flux**, with the measured signs and an order of magnitude short; exactly zero with fixed charges, and provably so |
| Polarization, blocking stress, free strain, work density | computed |

#### Why beta was exactly zero

With bond-charge-increment charges, a planar all-trans zigzag's dipole is *exactly*
independent of its backbone angle: the cell dipole is a sum over bonds of
`delta_ij (r_i - r_j)`, a backbone C-C bond carries no increment, the pendant bonds
have fixed lengths, and the zigzag's mirror pins each pendant pair's bisector
perpendicular to the chain axis whatever the backbone angle is. Measured, the
dipole held all ten printed digits (`mu_x = -0.7234964088` e.A) while the chain
repeat was driven from minus 1.17% to plus 1.14%. Axial strain changes the backbone
angle and nothing else, since bonds are rigid, so the dipole could not respond and
d33 and d31 were identically zero.

The dimensional term alone gave d33 = -4.43 pC/N against a measured -32. Its sign
is right, but that is arithmetic rather than physics: that term is negative on
every diagonal column for any stable crystal. The same term gave d31 = -0.14
against a measured +20, with the **wrong sign**.

#### What charge flux changes, and what it does not

`CrystalPacker(charge_flux=...)` lets the bond-charge increments depend on the local
geometry (`forcefield.FluxTopology`), which is the smallest change that breaks that
symmetry. The preset `pvdf-dft-valence-flux` carries two fitted coefficients;
`examples/fit_charge_flux.py` reproduces the fit and every number below.

**The symmetry does break, and by a measured amount.** Per unit `k_angle`, beta's
cell dipole gains +3.543 (C-H) and -6.469 (C-F) e.A per unit axial strain, against a
fixed-increment dipole that does not move at all. The two coefficients have opposite
signs because a zigzag's CH2 and CF2 bisectors point opposite ways -- the same
geometry that made the fixed model exactly zero.

**beta-PVDF, film convention** (axis 3 = poling = the crystal's polar x, axis 1 =
draw = the chain axis z), quoted with the poling axis along +P:

| | d33 | d32 | d31 |
|---|---|---|---|
| fixed charges, proper `d = e S` | 0 | 0 | 0 |
| fixed charges, plus the dimensional term | −4.43 | −5.90 | **−0.14** |
| charge flux, proper `d = e S` | −1.88 | −0.44 | **+2.33** |
| charge flux, plus the dimensional term | −6.29 | −7.09 | **+2.33** |
| measured (Nix and Ward 1986) | −32 | +1.5 | **+20** |

**d31's sign has turned positive**, which is the sharper test: d33's sign was
already right by arithmetic, while a positive d31 is something the dimensional term
cannot produce at all. The magnitudes are 5x (d33) and 9x (d31) too small. The two
independent routes to `d` -- a dipole derivative and a zero-stress root find --
agree to 0.33%, which is what says the flux entered the energy gradient and the
dipole derivative consistently rather than only one of them.

**And here is the honest reading of that +2.33.** Its sign is a property of the
fitting choice as much as of the data:

* The fit has **two parameters for twelve observations** and its `R^2` is **0.39**.
  The target is `dmu/d(axial strain)` for four chemistries from the sibling
  project's `field_neighborhood_refined`, an **exploratory GFN2-xTB** calculation on
  a *finite* two-chain pair in vacuum with the terminal backbone atoms pinned; its
  own `interpretation` field reads "finite-size, packing, stereochemistry and
  higher-level DFT validation outstanding. Strain is relative to fixed seed span,
  not a stress-free bulk lattice." It is not the PBE-D3 the rest of this potential
  is fitted to. It calibrates a mechanism and an order of magnitude, not a
  coefficient.
* Two things fit that data *better* and **neither is available to this model**: a
  bond-length flux channel (R² 0.86 on two parameters), inert here because
  `build_chain` places every atom at the polymer's own bond length -- measured, a
  bond-flux coefficient leaves beta's dipole holding every digit over the whole
  shape sweep -- and a plain charge-magnitude scale of x1.36 (R² 0.75 on one), which
  is not a flux at all but a statement that the increments are too small.
* **Fitted with a bond channel present the angle coefficient changes by a factor of
  six and both proper signs flip**: d33 = +0.03 and d31 = −0.06, dimensional total
  d31 = −0.33. What is shipped is the fit of *the model that is deployed* -- a
  rigid-bonded chain has only the angle channel -- so R² = 0.39 is honestly its own
  rather than borrowed from a channel that does nothing. The argument for that
  choice beyond self-consistency is that the angle-only family reproduces two
  experimental signs the fit never saw (d33 < 0 and d31 > 0) and the
  angle-plus-bond family reproduces neither.
* Leave one chemistry out and `k_angle(C-H)` moves over −1.62 .. −0.15, a factor of
  eleven. Every *converged* leave-one-out run keeps d31 positive (+0.50 to +2.16)
  and d33 negative, so the sign survives that resampling even though the magnitude
  spans a factor of four.
* Scaling the fitted coefficients by 0.25, 0.5, 1, 2 gives d31 = +0.52, +1.08,
  +2.33, +6.34. The magnitude is proportional to a badly determined number, and
  nothing about the value 2.33 means anything beyond its sign and its order.

**The flux is not free, and this is the part to read before using it.** The charges
enter the Coulomb sum, so the structure moves: beta's `c` goes 2.5469 -> 2.6052 A
(+2.3%), `a` 4.596 -> 4.492 A, `|P_x|` 0.116 -> 0.150 C/m², `C_11` 22.8 -> 28.4 and
`C_66` 3.4 -> 5.1 GPa, while `C_33` holds at 328 -> 330. Polyethylene stays
**exactly** zero (dipole below 1e-12 e.A, every `d` below 1e-9 pC/N) because its
repeat is centrosymmetric whatever the charges do, which is the null control the
flux had to pass and did. The model is **linear in the geometry with no
saturation**, so far from the reference it produces an unbounded electrostatic gain
from opening the backbone angle: at twice the fitted coefficient, and for
**alpha-PVDF at the fitted coefficient**, the relaxation runs to the line group's
+/-8 degree cap, the reference stops being axially stress-free and the two
piezoelectric routes disagree completely. Those rows are reported as non-measurements
rather than quoted.

#### What is still missing

* **Charge transfer along the backbone.** The increments are typed by element pair
  and a homonuclear pair has no orientation, so a C-C bond carries neither an
  increment nor a flux (`flux_topology` refuses one rather than orienting it by the
  arbitrary order of the bond list). That is the channel most likely to carry an
  *axial* dipole response, and it is outside this model.
* **Bond stretching**, hence the bond-flux channel, hence roughly half of a real
  chain modulus as well.
* **Electronic polarizability**: the dielectric constant is still 1 by construction.
* **The charge magnitude itself.** The single best one-parameter explanation of the
  reference residual is that the fitted increments are 36% too small, which is a
  statement about `docs/DFT_FIT.md`'s objective (no polarization data in it) rather
  than about flux.

So the honest statement is that beta's piezoelectric response has gone from
*structurally impossible* to *computable, with both measured signs and an order of
magnitude short*, and that the coefficient which gets it there is determined to
about a factor of four by data that does not describe a bulk crystal.

## Goal assessment

The target was a solver meeting field and strain measurement needs at 10 to 100
times existing methods. Split it, because the two halves land differently.

### Speed: met, on the axis it measures

| | Existing method | polyfind |
|---|---|---|
| Full response tensor for one chemistry | not obtainable from the 6-point grid | 0.3 s (beta) to 7.8 s (gamma) |
| One finite-difference derivative | 87 s (PVDF) to 143 s (VDCN) | included in the above |
| Structure search | 2,400-step budgets exhausted, continuations needed | 14 to 53 s, exhaustive |
| Field point at the density-functional tier | 772 to 2,470 s | not attempted; different accuracy tier |

Per derivative obtained the ratio is well past 100x. Per whole pipeline with the
structure search included, which the baseline does not do at all, it is 2 to 6x.
Both are true and neither should be quoted alone.

### Field and strain: qualitatively right for the first time, not quantitative

Charge flux makes the dipole respond to geometry, which fixed charges could not
do in principle. The out-of-sample test is the honest one, because the fit never
saw a measured piezoelectric constant:

| | d33 | d31 | interpretation |
|---|---|---|---|
| before | -4.43 | **-0.14** | d31 sign wrong |
| with flux | -6.29 | **+2.33** | both signs now match measurement |
| measured | -32 | +20 | magnitudes 5x and 9x short |

d31's sign is the meaningful result. d33 being negative proves nothing on its
own, since the dimensional term is negative for any stable crystal by
arithmetic; d31 turning positive requires a real mechanism, and it is stable
across the whole angle-only fit family including leave-one-chemistry-out.

The magnitudes should not be trusted, for reasons worth stating rather than
burying. The calibration data is exploratory semi-empirical work on a finite
pinned two-chain pair in vacuum, not a bulk crystal, and not the level the
energy model is fitted to. The fit's R-squared is 0.39. Two alternative
channels fit the *same* data better - a bond-stretch channel at 0.86 and a plain
charge rescale at 0.75 - and neither is usable here, one because bonds are rigid
and the other because it is not a mechanism. Fitting the angle channel alongside
a bond channel changes the coefficient sixfold and flips both signs back.
Leave-one-out spans the coefficient over a factor of eleven. The flux also
perturbs the structure it is meant to probe, moving beta's repeat by 2.3% and
C11 by 24%, and for two cases it drives the relaxation into the parametrisation's
angular cap, where the routes disagree completely; those are reported as
non-measurements rather than results.

**The magnitudes are not a mis-scaled electrostatics**, which was the obvious suspect
and has now been tested. Scanning the whole electrostatic strength over a factor of
4000 moves `d_33` and `d_31` the *opposite* way from the copolymer polar/antipolar
energy gap; the scale that would supply the magnitudes destroys the crystal before it
reaches them; and the gap has a non-electrostatic floor already 1.3x its reference.
`docs/ELECTROMECHANICS.md` section 5.7 has the scan. The conclusion is that these are
two deficiencies rather than one, and that the piezoelectric half is a missing
mechanism rather than a wrong coefficient.

### What Ewald costs

`CrystalPacker(coulomb="ewald")` (`DESIGN.md` 5.11) replaces the truncated Coulomb
sum with a proper lattice sum. Measured on the same machine, same structures,
tinfoil boundary, accuracy 1e-8 (which sets `alpha = 0.536 / A` at the packer's 8 Å
cutoff and `kmax = 4.61 / A`):

| | N atoms in the cell | k-vectors | real images | energy, one cell | with gradient | batched |
|---|---|---|---|---|---|---|
| beta | 12 | 164 | 315 | 1.04 → 4.35 ms (**x4.2**) | 3.14 → 9.45 ms (x3.0) | 0.85 → 3.33 ms/cfg (x3.9) |
| alpha | 24 | 348 | 245 | 3.05 → 9.76 ms (x3.2) | 5.47 → 24.7 ms (x4.5) | 1.82 → 10.0 ms/cfg (x5.5) |
| gamma | 48 | 718 | 125 | 16.2 → 28.8 ms (**x1.8**) | 35.0 → 72.4 ms (x2.1) | 15.4 → 31.0 ms/cfg (x2.0) |

So **two to five times the truncated kernel per configuration**, and the ratio
*falls* with cell size because the reciprocal sum is `O(N k)` against the real
sum's `O(N^2 I)`. Downstream:

| stage | truncated | Ewald | ratio |
|---|---|---|---|
| `refine_crystal`, beta | 0.1 s | 0.2 s | x1.6 |
| `refine_crystal`, alpha | 1.1 s | 1.9 s | x1.7 |
| `refine_crystal`, gamma | 3.9 s | 7.1 s | x1.8 |
| full response, beta | 1.33 s | 2.64 s | x2.0 |
| full response, alpha | 7.05 s | 14.0 s | x2.0 |
| full response, gamma | 79.3 s | 64.8 s | **x0.8** |

The response rows are the full deformable calculation of the section above
(`pvdf-dft-valence` + `pvdf-dft-valence-flux`, elastic tensor, both routes to `d`,
blocking stress, free strain), measured here with the same call both ways rather
than quoted against the 0.3 / 1.4 / 7.8 s recorded earlier on a differently loaded
machine. Gamma's **x0.8** is not Ewald being free: it is a relaxation that took a
different number of iterations, and it is the warning that these are wall times of
an optimisation, not of a kernel.

**What the response numbers do under Ewald**, which is the part that matters more
than the seconds:

| | C33 (GPa) | C11 (GPa) | \|P\| (C/m^2) | the two routes to `d` |
|---|---|---|---|---|
| beta | 330.4 → 330.6 | 28.4 → 26.0 | 0.1497 → 0.1543 | agree 0.33 % → **0.26 %** |
| alpha | 152.8 → 137.0 | 19.8 → 19.4 | 0.1022 → 0.1079 | 100 % → 100 %: a non-measurement either way, as already recorded |
| gamma | 99.3 → 102.3 | 13.4 → 15.9 | 0.1023 → 0.1075 | 0.51 % → **100 %**: a non-measurement *under Ewald only* |

`C33` — the one elastic constant this package stands behind — moves by 0.1 % for
beta and 3 % for gamma, and by 10 % for alpha. Beta's `d31` goes from +2.33 to
+2.43 pC/N, the same sign and the same order. The two routes to `d` agreeing
*better* for beta under Ewald (0.26 % against 0.33 %) is the sharpest available
check that the Ewald energy gradient and the Ewald dipole derivative are consistent
with each other.

**Gamma's route agreement is lost under Ewald and that is reported as a
non-measurement, not explained away.** It is not the gradient: finite differences
against `energy_and_grad` at gamma's own Ewald reference, with valence terms and
charge flux on, agree to 1e-9 on the cell, the coordinates and `c`. The failure is
in the converse route's zero-stress root find in a field, which is the same
parametrisation-cap failure already recorded for alpha's flux run; whether Ewald
merely moved gamma over that edge has not been established.

The refinement ratio is below the per-configuration ratio because a refinement
spends part of its time in the chain rebuild and the valence terms, which Ewald
does not touch. **The screen is unaffected**, by construction: the tabulated
chain-pair interaction needs a pairwise potential, Ewald's reciprocal half is not
one, and `pack(coulomb="ewald")` therefore screens with a truncated twin and
polishes with Ewald (`CrystalPacker._require_dsf` refuses any other arrangement).
So the 14-53 s structure-search rows above do not move at all.

**And what the extra cost buys is not what was expected.** Ewald shifts beta's
lattice energy by −1.12 kcal/mol per monomer, alpha's by −1.18 and gamma's by
−1.16, and moves the *structures* almost not at all: beta's refined cell goes from
4.59 x 8.59 x 2.6128 Å to 4.59 x 8.58 x 2.6141 and its `|P|` from 0.1408 to
0.1409 C/m^2. The α−β polymorph gap, which is an acceptance test, goes from +7.34
to +7.09 kJ/mol per monomer — a 3 % change, and **still outside** the −6.5 to −2.6
kJ/mol literature range, with the sign still backwards. (This row is the
illustrative potential, which fails that test either way. The sentence here
previously read "still inside the range it was already inside", which was wrong on
both counts; corrected 2026-09-11. On `pvdf-dft-valence`, the preset that passes,
the same switch moves the gap from −3.126 to −3.115 kJ/mol per monomer and it stays
inside.) The near-cancellation is why the truncated
sum got away with it for so long, and it is not a reason to keep using it: what
cancels is a nearly constant offset, not the differences a polar/antipolar
comparison is made of.

### The honest summary

For **screening** - ranking chemistries, getting the direction and rough scale of
a response, deciding which candidates deserve expensive validation - the solver
now does what was asked, at the speed asked. For **quantitative prediction** of
piezoelectric coefficients it does not, and the limiting factor is no longer the
search or the speed but the calibration data: a bulk, higher-level reference for
the dipole-strain response would be worth more than any further work inside this
package.

Still absent entirely: the shear constants involving the chain axis, which no
variable in this parametrisation can express; the bond-stretch channel;
polarizability; temperature; and domain switching.
