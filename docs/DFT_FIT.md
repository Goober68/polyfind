# Fitting `SimpleFF` to first-principles data

`DESIGN.md` section 5.7 fitted the built-in potential to crystal data and it did not
generalise.  The diagnosis there was about the *data*, not the optimiser: twelve
observations, several of them not independent, against five parameters, two of which were
structurally unidentifiable, and the quantities that would discriminate — torsion profiles
and the relative energies of conformers — were absent from the training set entirely.

This document records what happened when those quantities were supplied.  The code is
`polyfind.fitting` (the section after the crystal fit) and `polyfind.forcefield`
(`Frame` and friends); the script is `examples/fit_dft.py`; the result ships as the
opt-in preset `pvdf-dft-fit`, and `SimpleFF()`'s defaults are untouched.

**The headline, first.**  The fit generalises across chemistries — held-out energy error
falls from 3.32 to 1.76 kcal/mol on ten chemistries the objective never saw — and it fixes
**one of the three** known failures of the illustrative potential.  alpha-PVDF now comes
out below beta by 4.5 kJ/mol per monomer, inside the 2.6–6.5 range four independent
studies agree on and a sign change from the +7.4 the unfitted potential gives.  The
isolated-chain RIS ranking still prefers the 3/1 helix PVDF does not form, and alpha's
packed cell is still predicted polar.  Both of those are argued below to be limits of the
functional form rather than of the data.

---

## 1. The data

`$POLYFIND_TRAINSET` names an extended-XYZ file; nothing is vendored into this repository.
The set used here is 566 frames over 37 systems, every one of them a five-monomer
VDF-based oligomer (10 backbone carbons, capped `CH3-` and `-CF2H`), with `energy=`,
`system=` and `source=` on each comment line and per-atom `element x y z fx fy fz`.
Labels are PBE-D3.  Units in the file are eV and eV/A; `read_frames` converts both to
kcal/mol on the way in.

`source` is `scan` (a relaxed torsion scan) or `conf` (a conformer).  Both are used and
both are treated identically: what matters for the objective is that energies are only
comparable *within* a system.

### What is dropped, and why

`load_reference` removes a whole **system** — not a frame — when any of its frames has a
backbone that is unsaturated (a backbone C–C bond under 1.42 A) or passes through a ring.
A C=C torsion is a 60-plus kcal/mol barrier and a backbone inside a three-membered ring is
not a rotor at all; neither is expressible by one three-term Fourier series on a freely
rotating single bond, and including them would spend parameters on chemistry this
potential is not for.  Dropping frames rather than systems would be worse: it would leave
a system's survivors sampling only part of its own coordinate.

| dropped | shortest C–C (A) | energy spread (kcal/mol) | rule |
|---|---|---|---|
| `cnene` | 1.333 | 144.2 | backbone C=C |
| `ene` | 1.322 | 114.4 | backbone C=C |
| `dicnene` | 1.348 | 76.5 | backbone C=C |
| `fa_tetra` | 1.327 | 72.1 | backbone C=C |
| `cnepo` | 1.408 | 44.2 | backbone C=C |
| `epo` | 1.468 | 33.2 | backbone epoxide ring |

Everything kept spans 3.6 to 17.0 kcal/mol.  `cprop` is *kept* despite containing a ring:
its cyclopropyl is a pendant, and the test is on the backbone.  467 frames over 31 systems
remain.

## 2. Mapping a frame onto the potential

`SimpleFF` scores a `Structure` built by `build_chain`, whose topology comes from a
`Polymer`.  A reference frame names no polymer and its bond lengths and angles move.  So
`polyfind.forcefield.Frame` carries an arbitrary geometry with its topology *inferred*:

* `infer_bonds` — covalent-radius overlap, `r_ij < 1.25 (rcov_i + rcov_j)`.  Measured over
  the 467 kept frames, no distance falls anywhere near the cutoff: the longest bond sits at
  0.84 of its own (a C–F in `fan` at 1.389 A) and the closest non-bonded contact at 1.12
  (a C···I in `itfe` at 3.02 A), so a third of the scale is empty on either side of the
  decision.  (The margin is narrower on the *excluded* systems, 0.96 and 1.07, which is one
  more reason they are excluded.)
* `carbon_backbone` — the longest simple path through the carbons, by exhaustive
  depth-first search (the usual double-BFS trick is only correct on a tree, and a pendant
  cyclopropyl is not one).
* the backbone torsions are the consecutive quadruples along that path — the same set
  `Structure.dihedral_atoms` enumerates.

**Verification.**  `Frame.to_structure(PVDF)` permutes a frame's atoms into
`build_chain`'s own order and *checks* the result against a chain built from the polymer,
element by element and bond by bond, raising if either differs.  On the 15 pure-PVDF
frames it succeeds: 32 atoms, C10F10H12, 31 bonds, 10 backbone carbons alternating CH2 and
CF2 with a CH3 cap at one end and CF2H at the other — bond for bond what `build_chain`
builds for a five-monomer oligomer.  `SimpleFF.energy_frame(frame)` then equals
`SimpleFF.energy(frame.to_structure(PVDF))` to 1e-14 kcal/mol.  Both facts are asserted in
`tests/test_forcefield.py`; the second is what licenses fitting through the frame path at
all.

## 3. The objective

### Energies: relative, within each system, mean-centred

Absolute energies are not comparable across chemistries, so a per-system offset has to go.
The **mean** is subtracted, not the lowest frame.  The mean is the projection that removes
the offset optimally in least squares and treats every frame alike; referencing to the
lowest frame makes every residual of a system depend on that one frame's error, which
correlates them all and lets a single bad geometry shift a whole system's target.  The two
differ by a constant per system, which is exactly what centring removes, so nothing is
lost.

Residuals are divided by `sqrt(n_frames)`, so the reported number is a root-mean-square
error in kcal/mol and means what it says.

**Target accuracy.**  Matching DFT conformer energies to about 1 kcal/mol would be a good
classical potential.  Exact agreement is neither achievable nor the goal: this is a
rigid-geometry, fixed-charge potential with one three-term Fourier series per bond type
and no bonded terms at all.

### Forces: through the torsional projection, and worth almost nothing

Raw Cartesian forces cannot be fitted by a potential with no bond-stretch or angle-bend
terms.  Measured on this set, **93%** of the sum of squared reference force components is
accounted for by pure bond-stretch pair forces, and the fitted stretch constants are
textbook (C–H: k = 475 kcal/mol/A^2, r0 = 1.104 A, correlation between bond length and
force −0.84).  Every bond in the set is systematically *shorter* than the reference method
wants, which says plainly what the geometries are: relaxed at some other level of theory
and labelled with PBE-D3 single points.  Asking a bondless potential to reproduce that is
asking it to be wrong.

The projection that removes it is a **rigid rotation of one side of a bond**: it is the
only displacement that changes one torsion angle and leaves *every* bond length and *every*
bond angle untouched, so any bond or angle term contributes exactly zero to the projected
force, whatever its functional form.  `Frame.rotors()` builds those displacement fields
(normalised so `d phi_k / d theta_k = 1` exactly),
`Frame.reference_torques()` projects the reference forces onto them, and
`SimpleFF.frame_torques()` does the same for the model.  Both facts are tested:
a synthetic set of bond-direction force pairs projects to under 1e-10, and the rotor
changes torsion *k* by exactly the rotation angle and no other torsion at all.

What survives the projection is still not usable, and that is a measurement, not a guess:

| | value |
|---|---|
| reference torque, r.m.s. | 32.3 kcal/(mol rad) |
| illustrative potential's torque, r.m.s. | 9.7 |
| correlation between them | **+0.045** |
| is the *scanned* torsion's torque larger than the others'? | no (23.7 vs 23.0 for `adjfan`, and so on) |

A whole PVDF torsion barrier is about 4 kcal/mol, so a torsional gradient of 32
kcal/(mol rad) is not a torsion profile.  And along a properly relaxed scan the torque
about every *unconstrained* torsion must vanish, because the rigid rotor is a direction the
relaxation was free to move in; it does not vanish, which is the same conclusion the bond
lengths give.

So the forces enter at weight 0.01, chosen by held-out error, and the sweep says why:

| force weight | 0 | 0.003 | 0.01 | 0.03 | 0.1 | 0.3 |
|---|---|---|---|---|---|---|
| held-out energy r.m.s. (kcal/mol) | 1.763 | 1.760 | **1.757** | 1.773 | 1.883 | 2.042 |
| held-out torque r.m.s. | 32.45 | 32.45 | 32.44 | 32.43 | 32.45 | 32.49 |

The forces are worth a third of a percent, and the torque residual barely moves whatever
weight is put on it — the fitted model explains essentially none of the reference torque's
variance.  They are in the objective mostly so that the number is on the record.  The
richness of the force data is real (about 3,300 torsional projections against 467
energies) and it buys nothing, because at these geometries the signal is not the one this
potential can carry.

## 4. The split: by system, never by frame

Frames within one system are strongly correlated — a torsion scan is one molecule at a
dozen angles, sharing every error its geometry makes.  A random frame split would put nine
in train and three in test and report a small test error meaning nothing, which is how the
fit in section 5.7 was flattered.  `split_systems` therefore holds out whole chemistries:
every third system in alphabetical order, which is reproducible, has nothing to do with
how well any system fits, and spreads the held-out set over the substituent types.  PVDF
is forced into the training set — it is the chemistry every acceptance test is about.

* **train, 21 systems / 312 frames**: `adjfan btfe cf2h cprop ctfe_ter fan_wide fancf3
  fanome hfo ib isocy man mve nitro pmve pvdf scn tfp vdc vdf_trfe vsf`
* **held out, 10 systems / 155 frames**: `an cfe_ter fan fancl hfp itfe mvs pp tfe vdcn`

## 5. Parameters, and what the data can and cannot pin down

27 numbers (`REF_VARIABLES`):

* three Fourier coefficients for each of four **torsion types**, classified by the pendant
  elements on the two central backbone carbons: `CF2-CH2` (922 of the 3,269 backbone
  torsions, all 31 systems), `CF2-CHF` (1,402, 30 systems), `CF2-CF2` (101, 6 systems) and
  `other` (844, 28 systems — wherever a substituted carbon is central).  Every backbone
  bond of PVDF is `CF2-CH2`, so that one triple *is* the shipped preset's `torsion`;
* multipliers on the UFF Lennard-Jones minimum distance and well depth of C, H, F and Cl;
* six bond-charge increments (C–H, C–F, C–Cl, C–N, C–O, C–S);
* the 1–4 nonbonded scale.

The charge model deserves a note: `bci_charges` assigns `q_i = sum_j delta(e_i, e_j)` over
bonds, which is neutral by construction and *is* the charge model `polymers.py` already
uses, written as parameters instead of literals.  `{C-H: -0.10, C-F: +0.20}` reproduces
PVDF's table exactly (CH2 carbon −0.20, H +0.10, CF2 carbon +0.40, F −0.20), `{C-H: -0.06}`
reproduces polyethylene's and `{C-H: -0.10, C-Cl: +0.10}` PVDC's.  Because every
backbone–backbone bond is C–C and a same-element bond carries no increment, the charges of
one periodic block do not depend on the bonds crossing its boundary, which is what lets the
same numbers reach `pack` and `refine` unchanged.

### The identifiability traps of section 5.7

**Permittivity against charge scale: did not lift, and cannot.**  Section 5.7 hoped that
relaxed, off-ideal geometries might separate them.  They do not.  Every energy this
potential produces depends on the two only through `q_i q_j / eps_r`, so scaling every
charge by *s* and `eps_r` by *s²* leaves every energy — and therefore every force, which is
a derivative of an energy — identically unchanged.  Measured on the real frames:
multiplying every increment by 1.7 and dividing the Coulomb term by 1.7² moves no frame
energy by more than **1.2e-13 kcal/mol** and no torque by more than 4.8e-14.  This is a
property of the functional form, not of the geometries, and no amount of energy or force
data will break it; only a polarization can, and there is none here.  So `eps_r` is held at
1 and the increments carry the whole electrostatic scale.  `permittivity_degeneracy`
performs the measurement rather than asserting the conclusion.

**The three Fourier coefficients: did lift.**  At ideal torsion angles the third coefficient
vanishes identically and only `V1 + V2` enters, which is why the crystal fit could only
move one of them.  The reference scans sample the whole circle, so all three are separately
determined here, and the fit uses that: for `CF2-CH2` it lands on
`V1 = +1.740, V2 = -1.741, V3 = +2.564`.  Note what that means — `V1 + V2 = -0.001`, so at
ideal angles trans and gauche are *degenerate* and the whole gauche preference has moved
into the shape of the barrier, which is exactly the freedom the crystal fit did not have.

### Parameters that ran to their bounds — read this before trusting the numbers

Five did: `x_C`, `x_H`, `D_C` and `D_H` all at the low end (0.90 of UFF on the radii, half
of UFF on the well depths) and `q_C-F` at the low end (0.02 e).  That is a lot, and it is
worth being explicit that it is the bounds, not the data, holding four of them.  Unbounded,
the fit does something worse: it drives `x_C` and `x_H` to whatever floor it is given, cuts
the well depths of C, H and F to a fifth of UFF, and turns fluorine *positive* (`q_C-F`
goes to −0.07, reversing the CF2 dipole every polarization result in this package rests
on).  It can do all of that because relative conformer energies within one molecule
barely see the overall dispersion, so that term is nearly free to be traded against the
torsions — while a crystal, which is held together by exactly that dispersion, sees it
immediately.  The bounds are where this data stops being able to tell:

* Lennard-Jones, 10% on the minimum distance and a factor of two on the well depth: the
  same reasoning `CELL_TOL` uses in the other direction, that a rigid-geometry fixed-charge
  potential has no business claiming better than about ten percent on a van der Waals
  contact distance;
* charge increments, signed by electronegativity: carbon is more electronegative than
  hydrogen and less than F, Cl, N and O, so H comes out positive and those four negative.
  That is not a fact about any target here.  Sulfur and carbon are within 0.03 on the
  Pauling scale, so C–S is left free.

A parameter sitting on a bound is a parameter the data wanted to push somewhere
unphysical.  `q_C-F = 0.02` in particular means the fit would rather have almost no C–F
electrostatics at all, and that is the most likely reason acceptance test 3 still fails.

## 6. The fit, and what it costs

Bounded Levenberg–Marquardt (`scipy.optimize.least_squares`, trust-region reflective) from
four starting points — the illustrative potential plus three perturbations.  Unlike the
crystal fit, this objective is smooth and cheap: every geometry-dependent quantity
(distances, dihedral angles, exclusion masks, the torsional projections `dr_ij/dphi_k`) is
precomputed once by `ReferenceDesign`, so one evaluation is about 20 ms of flat array
arithmetic and the whole fit is 40 seconds.  `ReferenceDesign.evaluate` reproduces
`SimpleFF.energy_frame` and `SimpleFF.frame_torques` to 1e-8; the test asserts it, because
a design that drifted from the shipped potential would make every number below meaningless.

A weak ridge toward the illustrative values (weight 0.3, spread over the 27 parameters)
keeps parameters with almost no data behind them near where they started.  It does not cost
anything and it slightly helps, mostly by pulling parameters off their bounds:

| ridge | 0 | 0.05 | 0.3 | 1.0 | 3.0 |
|---|---|---|---|---|---|
| held-out energy r.m.s. | 1.774 | 1.766 | **1.757** | 1.742 | 1.736 |
| parameters sitting on a bound | 7 | 6 | 5 | 5 | 3 |

0.3 is what the shipped vector was fitted at.  Held-out error is *monotonically* better with
a stronger ridge, which is worth stating plainly: it means the last few percent of training
fit is overfitting, and that a still-tighter prior would be defensible.  The whole range is
two percent, below what this fit distinguishes, so nothing was chosen to make a number look
good; but nothing here argues for 0.3 over 1.0 either.

### Fitted values

| parameter | illustrative | fitted | | parameter | illustrative | fitted |
|---|---|---|---|---|---|---|
| V1 [CF2-CH2] | +1.300 | **+1.740** | | x_C | 1.000 | **0.900** ← bound |
| V2 [CF2-CH2] | −0.050 | **−1.741** | | x_H | 1.000 | **0.900** ← bound |
| V3 [CF2-CH2] | +2.500 | **+2.564** | | x_F | 1.000 | 1.059 |
| V1 [CF2-CHF] | +1.300 | +1.737 | | x_Cl | 1.000 | 0.993 |
| V2 [CF2-CHF] | −0.050 | −1.576 | | D_C | 1.000 | 0.500 ← bound |
| V3 [CF2-CHF] | +2.500 | +1.699 | | D_H | 1.000 | 0.500 ← bound |
| V1 [CF2-CF2] | +1.300 | −0.105 | | D_F | 1.000 | 0.576 |
| V2 [CF2-CF2] | −0.050 | −1.016 | | D_Cl | 1.000 | 1.030 |
| V3 [CF2-CF2] | +2.500 | +1.673 | | q C–H | −0.100 | −0.259 |
| V1 [other] | +1.300 | +1.072 | | q C–F | +0.200 | +0.020 ← bound |
| V2 [other] | −0.050 | −1.021 | | q C–Cl | +0.100 | +0.050 |
| V3 [other] | +2.500 | +2.302 | | q C–N | +0.050 | +0.111 |
| | | | | q C–O | +0.050 | +0.035 |
| | | | | q C–S | 0.000 | −0.101 |
| | | | | scale14 | 0.500 | **0.940** |

x multipliers are on the UFF minimum distance, D on the well depth; the shipped preset's
absolute values are C 3.466, H 2.597, F 3.563, Cl 3.918 A.  `D_C` and `D_H` are also at
their bound, which the note above covers.

The Jacobian of the (centred) energy block at the optimum has a condition number of 504
over all 27 parameters — no null direction, which is what "all three Fourier coefficients
are separately determined" means quantitatively.  Their pairwise correlations for the
shipped type are −0.65 (V1,V2), +0.26 (V1,V3) and −0.44 (V2,V3): correlated, as Fourier
coefficients on a bounded angular range always are, but not degenerate.

## 7. Results

### Held out by chemistry

Energy r.m.s. error, kcal/mol, after removing each system's own mean:

| | train (21 systems, 312 frames) | held out (10 systems, 155 frames) |
|---|---|---|
| illustrative | 3.156 | 3.316 |
| **fitted** | **1.498** | **1.757** |

Per held-out system: `tfe` 0.83, `itfe` 0.88, `an` 0.92, `fancl` 1.08, `cfe_ter` 1.08,
`vdcn` 1.12, `hfp` 1.43, `fan` 1.57, `pp` 1.67, `mvs` 3.93.  PVDF itself (in training) is
1.13.  The one bad system in each half is a sulfur one — `mvs` 3.93 held out and `vsf` 3.32
in training, both sulfone-containing, where the charge model has one increment for the C–S
bond and nothing for S=O.  Excluding `mvs` the held-out error is 1.21 kcal/mol over the
other 137 frames, which is the target this set out for; with it the number is 1.76, and it
is counted.

The torque tells a different story, and the cleanest way to put it is this: the held-out
reference torques have r.m.s. **32.56** kcal/(mol rad), so a model that simply predicted
*zero* would score 32.56.  The unfitted potential scores 33.83 — worse than predicting
nothing.  The fitted one scores 32.44, better than nothing by 0.4%, with a correlation
against the reference of +0.10 (up from +0.07).  Whatever is in those forces, this
functional form cannot represent it, for the reasons in section 3.

### Acceptance tests

None of these is in the objective; the fit never sees a crystal, a polymorph or an RIS
ranking.

| | illustrative | fitted | target | |
|---|---|---|---|---|
| 1. E(alpha) − E(beta), kJ/mol per monomer | +7.38 | **−4.54** | −6.5 to −2.6 | **PASS** |
| 2. isolated-chain RIS, top candidate | TG+ | TG+ | anything but the 3/1 helix | **FAIL** |
| 3. alpha E(antipolar) − E(polar), kcal/mol per monomer | +0.197 | +0.091 | ≤ 0 | **FAIL** |

**Test 1 passes**, and it is the one with the most agreement behind it: four independent
studies over five exchange-correlation functionals put beta 2.6 to 6.5 kJ/mol per monomer
above alpha (`REFERENCES.md` section 5).  The unfitted potential had the sign backwards and
the magnitude three to six times too large; the fit lands at 4.5 kJ/mol, near the middle of
the range.  This is DESIGN.md 5.4's headline failure, fixed by data the objective did see
(torsion profiles) applied to a quantity it did not.

**Test 2 fails**, but not by as much.  At the finer 10-degree scan step the illustrative
potential puts TG+ 1.64 kcal/mol per monomer clear of the next candidate; the fit narrows
that to 0.23 (TG+ −3.23, TTG+G+ −3.00).  The cause is visible in the RIS energies: the
gauche energy per bond goes from +0.19/+0.23 (CH2/CF2) to −1.59/−1.60, i.e. the fitted
potential makes gauche *strongly* favoured on the isolated chain, where the truth is near
degeneracy.  That is a 1.8 kcal/mol swing in a quantity whose real value is a few tenths,
made by a fit whose own per-frame error on PVDF is 1.13 kcal/mol.  The model has no
resolution left for it, and no amount of the same data would give it any.

**Test 3 fails**, halved.  The antipolar alpha cell still costs +0.091 kcal/mol per monomer
against +0.197 before.  The fitted charges are the obvious suspect: `q_C-F` runs to its
lower bound, beta's polarization falls from 0.140 to 0.114 C/m² (the DFT value is
0.176–0.188), and section 5.7's ablation already showed the electrostatics is what this
test turns on.  Intramolecular conformer energies simply do not constrain the
electrostatic scale — a polarization would, and there is none in this data.

### Which parameters did the work

Each group takes its own value from the full fit, everything else left illustrative:

| parameter set | train E | held-out E | E(a)−E(b) kJ | RIS top | E(anti)−E(polar) | beta \|P\| |
|---|---|---|---|---|---|---|
| illustrative | 3.156 | 3.316 | +7.38 | TG+ | +0.197 | 0.140 |
| torsions only | 2.614 | 3.115 | +2.61 | TG+ | +0.197 | 0.140 |
| Lennard-Jones only | 2.481 | 2.528 | +2.58 | **G+G+** | +0.059 | 0.144 |
| charges only | 2.752 | 3.182 | +4.64 | TG+ | +0.220 | 0.112 |
| 1–4 scale only | 3.060 | 3.360 | +7.38 | TG+ | +0.197 | 0.140 |
| **full fit** | **1.498** | **1.757** | **−4.54** | TG+ | +0.091 | 0.114 |

Two things come out of this, one comfortable and one not.

*Comfortable*: every group helps and none is close to the whole.  The full fit's held-out
error (1.757) is far below the best single group (2.528), which is the opposite of section
5.7's pattern, where two parameters transferred and three fitted noise.  This is a fit with
enough data behind it that its parameters are doing joint work.

*Not comfortable*: **different groups carry different acceptance tests, and nothing carries
both.**  The Lennard-Jones parameters alone dislodge the 3/1 helix from the top of the RIS
ranking (test 2) and halve the antipolar cost (test 3), while it takes the whole vector —
torsions included — to get the alpha/beta ordering right (test 1), and doing so puts TG+
back on top.  A three-term Fourier torsion plus a 12–6 Lennard-Jones plus fixed point
charges on rigid geometry appears not to be able to have both at once, and that is a
statement about the functional form, not about the data or the optimiser.

## 8. What would be next

Not more data of this kind.  Three specific things, in order of how much they would buy:

1. **A polarization or a dipole in the training set.**  It is the only observable that
   separates `eps_r` from the charge scale (section 5), and it is what test 3 turns on.
   Molecular dipoles of these same oligomers, from the same calculations, would cost
   nothing to produce and would pin the electrostatics the conformer energies cannot.
2. **Reference geometries relaxed at the labelling level.**  Then the forces would be
   torsional gradients rather than a record of a bond-length disagreement, the torque term
   would carry the roughly seven-to-one advantage in count that it structurally has, and
   the torsion profile would be pinned far harder than 3,300 useless projections pin it.
3. **Backbone angles as degrees of freedom.**  The frozen-angle approximation is what
   forces the same Lennard-Jones radius to serve an intramolecular 1–5 contact and an
   intermolecular packing distance, which is the tension section 7's ablation exposes.
   `docs/CHEMISTRY_EXTENSION.md` already proposes this for other reasons.

