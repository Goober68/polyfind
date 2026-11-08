# Candidate chains and defect placements for the loss / frequency / strain optimum

*First pass 2026-11-07, corrected 2026-11-08 for the handedness of the TG helix (section 2: the
first pass ran one of an isotactic sequence's four TG conformers, and for every chlorine and
nitrile row the other handedness is 0.9-1.5 kcal/mol per monomer lower, which turned the defect
rows' margins from near zero to +0.57-0.97 and the reading of the defects from margin tuners into
TG-favouring pins), on what the model can say today: the crystal terms of
`pvdf-dft-valence-belt4` (docs/CRYSTAL_STAGE.md 5.22), the pendant-group relaxation of the
same section, the loss model of docs/ELECTROSTRICTION_RESEARCH.md 11-12 and the 2025-2026
literature pass. Code: `examples/candidate_screen.py` (the crystal screen),
`examples/candidate_loss_map.py` (the mapping onto the loss model); results
`docs/predictions/candidate_screen_belt4.json`, `candidate_loss_map_belt4.json`. Every number here
is a model number at zero kelvin on a regular ten-monomer periodic sequence packed as two chains;
section 7 says what that leaves out and what decides each row, section 6 reads Dipole's slate of
other chemistries against the same screen, and section 8 adds the frequency: the kink-pair formation
energy in the lattice as the first piece of the dynamics.*

## 1. What "optimal tan delta, frequency and strain" means in the program's variables

The giant strain of the relaxor fluoropolymers is a field-driven conformational switch, TGTG' or
T3GT3G' to TTTT, inside and around the crystals (literature pass, finding 1); the material's
dielectric loss, its relaxation frequency and its strain are three readings of the same slow
polarization channel (finding 3). In the two-state picture the plan is built on (plan section 3):

* a unit of `n` monomers switches at the energy cost `n dE` and gains `n dmu` of dipole along
  the field; the **threshold field** `E* = dE / dmu` is intensive and decides whether a
  chemistry works at 50-100 MV/m at all;
* the **strain** is the switched fraction times the **transition strain** of the conformational
  change, a geometric number (axial +11% for TGTG' to TTTT in PVDF; transverse of the opposite
  sign and about half the size);
* the channel relaxes with `tau(T) = tau0 exp(Ea / k (T - Tf))`; below `f_relax = 1 / 2 pi tau`
  the loss per unit strain at the drive frequency is `omega tau`, so at fixed drive frequency the
  only lever on **loss per strain** is `tau`: the local barrier `Ea`, the cooperativity `Tf`, the
  attempt frequency;
* the operating point the program designs for is a DC bias that holds the polar fraction near its
  steepest point with a small AC drive on top (the bias-hold point: 55 V/um bias, 1.8 V/um AC,
  34-66 kHz, a 4 MW/m^3 heat budget; plan sections 11-12).

So a chemistry has three independent handles, and the pass below reads them separately:

| handle | sets | who moves it | what computes it today |
|---|---|---|---|
| margin `dE` and dipole gain `dmu` | `E*`, the slope of S(E) at the bias point, the polar fraction at bias | TrFE content; defect type, concentration and placement | the crystal screen (tier 2) |
| transition strain | the stroke per switched unit | the conformational pair, hence the chemistry's ground conformation | the crystal screen, exactly |
| `Ea`, `Tf`, `tau0` | `f_relax`, hence loss per strain at 34-66 kHz and the heat | defect pinning (Ea), the size of the cooperative unit (Tf), chain stiffness | the kink-pair formation energy in the lattice (section 8), the first piece, as a bracket between the relaxed kink crystal and the pinned lattice; the migration barrier and `Tf` are not computed (the bench's Vogel-Fulcher fit is the input) |

and a fourth that is a property of the sequence rather than of any unit: the **spread of local
margins** along the chain, which the driven lattice (plan section 12) found to be the knob that
suppresses the stroke at bias and broadens the loss peak. A controlled sequence is the one place
where that spread is a design variable rather than a statistical accident, which is what the
microreactor programme is for.

## 2. The screen

Ten-monomer model sequences, each packed as a two-chain crystal in TTTT (the polar, beta-like
conformation) and in the TGTG' helix (PVDF's non-polar ground) in its two handednesses and two
phases (TGTG', TG'TG, GTG'T, G'TGT: four conformers for an isotactic sequence, one for PVDF), under
`pvdf-dft-valence-belt4` with every CH2, CF2, CHF, CFCl and CHCl group relaxed as a rigid tilt
(5.22). From the pair: the margin `E(TTTT) - E(TG)` per monomer, the TG side being the lowest of
the four conformers that ran, the axial and areal transition strains, the dipole gained per
monomer (`P V / n` of the TTTT cell minus the TGTG' cell's) and `E* = dE / dmu`.

VDF and TrFE are the validated chemistry: their crystal terms were fitted to Dipole's PBE-D3
anchors and read back within 1.3 tolerances on everything anchored. CFE (CH2-CFCl), CTFE
(CF2-CFCl) and CDFE (CHCl-CF2) carry a chlorine whose nonbonded and charge parameters come from
the conformer valence fit alone; no crystal anchor has tested them, so those rows are
extrapolations. The E2 request to Dipole (docs/E1_E2_PREDICTIONS.md) is what validates them.

| sequence | chemistry | margin TT-TG (rigid) kcal/mol per monomer | TG tilt relaxation (lowest TG conformer) | axial, areal switch strain | P(TTTT) C/m^2 | mu(TTTT) e.A | E* (mu_TT) MV/m | E* (mu_TT - mu_TG) MV/m |
|---|---|---|---|---|---|---|---|---|
| pvdf | VDF/TrFE (validated) | +0.744 (+0.539) | -0.205 | +10.9%, -15.7% | 0.193 | 0.629 | 513 | 1208 |
| vdf7-trfe3 | VDF/TrFE (validated) | +0.145 (-0.172) | -0.357 | +10.6%, -17.2% | 0.167 | 0.560 | 112 | 324 |
| vdf8-trfe2 | VDF/TrFE (validated) | +0.177 (-0.127) | -0.343 | +10.5%, -16.4% | 0.174 | 0.578 | 133 | 448 |
| vdf6-trfe4 | VDF/TrFE (validated) | -0.329 (-0.735) | -0.441 | +10.3%, -16.7% | 0.162 | 0.546 | 0 | 0 |
| vdf5-trfe5 | VDF/TrFE (validated) | +0.752 (+0.526) | -0.253 | +10.9%, -19.6% | 0.160 | 0.537 | 607 | 1360 |
| vdf7-trfe2-cfe1 | chlorine (extrapolated) | +0.966 (+0.756) | -0.332 | +11.5%, -14.5% | 0.159 | 0.566 | 740 | 1163 |
| vdf7-trfe2-cfe1-adj | chlorine (extrapolated) | +0.961 (+0.758) | -0.326 | +11.5%, -14.6% | 0.159 | 0.566 | 736 | 1156 |
| vdf7-trfe2-ctfe1 | chlorine (extrapolated) | +0.790 (+0.439) | -0.455 | +11.1%, -12.3% | 0.136 | 0.506 | 676 | 1043 |
| vdf7-trfe2-cdfe1 | chlorine (extrapolated) | +0.146 (-0.267) | -0.447 | +11.3%, -8.8% | 0.130 | 0.532 | 119 | 240 |
| vdf9-cfe1 | chlorine (extrapolated) | +0.982 (+0.755) | -0.308 | +11.4%, -14.7% | 0.176 | 0.614 | 693 | 1425 |
| vdf8-cfe2 | chlorine (extrapolated) | +0.835 (+0.506) | -0.416 | +10.6%, -13.1% | 0.147 | 0.536 | 675 | 1832 |
| vdf7-trfe2-an1 | chlorine (extrapolated) | +0.568 (+0.334) | -0.288 | +11.6%, -16.9% | 0.146 | 0.558 | 441 | 802 |
| vdf7-trfe2-fan1 | chlorine (extrapolated) | +0.912 (+0.912) | +0.000 | +11.0%, -12.9% | 0.144 | 0.581 | 681 | 1155 |
| vdf9-an1 | chlorine (extrapolated) | +0.470 (+0.213) | -0.300 | +11.7%, -16.3% | 0.163 | 0.590 | 345 | 841 |
| vdf6-trfe3-cdfe1 | chlorine (extrapolated) | +0.805 (+0.480) | -0.362 | +10.8%, -8.4% | 0.127 | 0.523 | 668 | 1020 |
| vdf6-trfe3-cfe1 | chlorine (extrapolated) | +0.297 (+0.000) | -0.424 | +10.1%, -14.4% | 0.153 | 0.549 | 235 | 707 |
| vdf9-trfe4-cfe1-14 | chlorine (extrapolated) | +0.110 (-0.187) | -0.420 | +10.5%, -14.2% | 0.156 | 0.552 | 87 | 279 |
| vdf10-trfe4-14 | VDF/TrFE (validated) | -0.051 (-0.393) | -0.382 | +10.4%, -16.5% | 0.168 | 0.562 | 0 | 0 |
| vdf9-trfe5-14 | VDF/TrFE (validated) | -0.058 (-0.413) | -0.393 | +10.5%, -17.1% | 0.164 | 0.551 | 0 | 0 |
| vdf9-trfe4-cdfe1-14 | chlorine (extrapolated) | +0.402 (+0.077) | -0.381 | +10.8%, -9.6% | 0.134 | 0.529 | 329 | 574 |
| vdf8-trfe4-cfe1-cdfe1-14 | chlorine (extrapolated) | +0.641 (+0.303) | -0.458 | +11.3%, -14.0% | 0.129 | 0.526 | 529 | 1479 |

**Reading the screen.**

* **The helix handedness is the first lesson of this pass, and it corrected the pass.** The first
  run (2026-11-07) packed one TG conformer, TGTG'. For the VDF/TrFE copolymer rows that is enough:
  the 80/20 row's four conformers are degenerate to 0.0002 kcal/mol per monomer, and for a lone
  defect in a VDF run it hardly matters (VDF/CFE 90/10 and 80/20, VDF/AN 90/10: 0.01-0.05). For
  every chlorine and nitrile defect beside TrFE at one per ten it is not: TG'TG, the other handedness, is 0.9-1.5 kcal/mol per monomer below
  TGTG' (CFE 0.94, CFE next to TrFE 0.90, CTFE 0.89, AN 1.30, FAN 1.48), because the bulky
  substituent on a backbone carbon in a gauche segment fits one handedness of the helix and not
  the other (in the other, CDFE's chlorine comes within 1.8 A of the TrFE fluorine two monomers
  away and the chain cannot be built). The two handednesses are different conformers in this
  builder, not mirror images (RMSD 1.65 A against the mirror image), and the gap is not a search
  artefact: refining all three packing candidates of the CFE row's TGTG' chain leaves them 0.9-1.2
  above all three of its TG'TG chain. So the first pass's margins near zero for the defect rows
  (CFE +0.03, CFE adjacent +0.06, CTFE -0.10, AN -0.73) were the artefact of one handedness; the
  table above carries the lowest conformer and the corrected margins are +0.57 to +0.97. CDFE's
  other handedness cannot be built at the rigid geometry, so its +0.15 stands with that caveat.
* **TrFE content sets the margin, in steps of about 0.3 kcal/mol per monomer per ten mol%,** and
  this is handedness-free. PVDF's TGTG' is 0.74 below TTTT and its crystal threshold is 500 MV/m:
  the homopolymer cannot be field-switched, as is known (beta is made by drawing). At 20 and 30
  mol% TrFE the margin is +0.18 and +0.15 and the threshold 130 and 110 MV/m, the field window in
  which the giant-strain copolymers and terpolymers are driven; at 40 mol% TTTT is the ground
  conformation. The crossover sits between 30 and 40 mol% in this model.
* **The chlorine and nitrile defects are TG-favouring pins, not margin tuners.** One per ten
  beside 20 mol% TrFE takes the 80/20 copolymer's +0.18 to +0.97 (CFE), +0.79 (CTFE), +0.91 (FAN)
  and +0.57 (AN): the defect pins its neighbourhood in the non-polar helix by 4-8 kcal/mol per unit
  (section 4), and none of the ten-monomer defect cells switches at the operating field (E* 440-740
  MV/m). CDFE alone leaves the margin at +0.15, costing both conformations the same 17 kcal/mol per
  unit. In the real terpolymer the pins sit one per ten to fourteen monomers in a 30-35 mol% TrFE
  matrix, and the matrix between them is where the near-degeneracy lives; the ten-monomer cell
  with 20 mol% TrFE has no monomer more than three from a pin, so it is the wrong model of the
  relaxor's switching unit. The fourteen-monomer cells (below) put the pin at the film's dilution.
* **The pendant relaxation decides the sign of the small margins.** With rigid pendants the
  70/30 and 80/20 models had TTTT lowest (-0.17, -0.13); their TG cells release 0.34-0.44
  kcal/mol per monomer through their CF2 and CHF groups (CF2 wag about 3 deg, CHF wag and twist
  2-3 deg) and TTTT cells 0.03-0.04, which turns both margins around. The bend constants set that
  release to within about 20% (5.22), so the sign of a margin under 0.1 is not the model's to
  call: the 20-30 mol% window is near-degenerate, which is the physically right statement and the
  one Dipole's E2 cells decide.
* **The alternating 50/50 row is an artefact to keep in view.** A regular isotactic alternation
  (every CHF fluorine on the same side) packs its TGTG' 0.75 below TTTT where the real,
  statistical 50/50 copolymer is an all-trans ferroelectric; the model cell is the regular
  sequence, not the copolymer. It also says that sequence regularity moves the packing energy by
  half a kcal/mol per monomer at fixed composition: for a controlled-sequence material the pattern
  is a lever as large as the composition, and a model cell of a controlled sequence is then the
  right object, not a stand-in.
* **The switch strain is geometric and chemistry-insensitive:** +10 to +12% along the chain for
  every row, with a lateral contraction of 9-20% in the model.
  The lateral number is exaggerated: PBE-D3 packs beta 3% too dense against experiment and the
  model inherits that from its anchors, so the volumetric change of the switch (6-8% here against
  2.5% between the measured alpha and beta densities) is an upper bound.
* **The dipole gained per monomer** is 0.54-0.63 e.A for the fluorine chemistries and falls with
  chlorine (0.51-0.57); the TTTT polarization falls from 0.193 (PVDF) through 0.167 (70/30) to
  0.159 (CFE), 0.136 (CTFE) and 0.130 (CDFE) C/m^2. Every defect costs stroke per switched unit.

**The fourteen-monomer cells (`candidate_screen_belt4_14mer.json`).** The film the loss model is
anchored on is 63/30/7; one regular cell of that composition is fourteen monomers with one CFE,
`VVRVVRVCVVRVVR` (64/29/7), and its controls are the same cell with the CFE replaced by VDF
(71/29) and the matrix at 36 mol% TrFE (64/36). With TTTT, TG'TG and TGTG' (the two phases were
never the lowest above): 64/29/7 TTTT -7.138, TG'TG -7.248, TGTG' -7.246 (the two handednesses
degenerate to 0.003 here, where the ten-monomer CFE cell split them by 0.94: the split is a
property of the CFE's placement among the TrFE units, not of the chlorine alone), margin +0.11
kcal/mol per monomer, E* 280 MV/m on the TT-TG dipole difference. The matrix control 71/29 (the
same cell, CFE replaced by VDF): TTTT -7.666, TG'TG -7.614, margin -0.05 (TGTG' running): the 29
mol% matrix is degenerate within the model's resolution, and the pin tips the cell by +0.16 per
monomer, i.e. +2.3 kcal/mol per CFE unit (TT cost +7.4, TG cost +5.1 in this more open lattice,
against +6.9 and -1.0 in the ten-monomer 20 mol% cell): the pin is a gentler pin where the
matrix is near-degenerate. The 64/36 matrix control: TTTT -7.239, TGTG' -7.182, margin -0.06. So in
the fourteen-monomer cell the matrix is degenerate within the model's resolution from 29 to 36
mol% TrFE (the ten-monomer 70/30 and 80/20 cells gave +0.15 and +0.18, the 60/40 cell -0.33: the
cell length moves these numbers by 0.2, which is the sequence-regularity effect of the 50/50 row
again), and the pin sets where in that window the material sits. CDFE in place of the CFE
(`candidate_screen_belt4_14mer_b.json`): TTTT -6.250, TGTG' -6.652 (TG'TG not buildable, the F-Cl
clash), margin +0.40: at one per fourteen the CDFE unit tips the degenerate matrix by +0.45 per
monomer (TT cost +19.8, TG cost +13.5 per unit, pin +6.3) where CFE tips it by +0.16 (+2.3 per
unit). So the pin strengths of section 4's table are properties of the unit *in its cell*, not of
the unit: CFE pins at +7.9 per unit in the ten-monomer 20 mol% cell and +2.3 in the fourteen-monomer
29 mol% one, CDFE at -0.3 and +6.3. The fourteen-monomer cells at the film's composition are the
ones to read, and in them CDFE is the stronger tipper and leaves the cell at +0.40, off the steep
part of the curve (E* 570 MV/m), while CFE leaves it at +0.11. With one CFE and one CDFE per
fourteen (`DVRVVRVCVVRVVR`, 57/29/7/7): TTTT -5.928, TG'TG -6.569, margin +0.64. The two pins add
(-0.05 + 0.16 + 0.45 = +0.56 against +0.64 found), so at the film's dilution every chlorine unit
tips the degenerate matrix toward the helix and two of them tip it off the curve. So the pin diluted to one per fourteen in a 29 mol% matrix
leaves the cell near-degenerate where one per ten in a 20 mol% matrix pinned it at +0.97: at the
film's composition the model reproduces the relaxor's near-degeneracy, and the quantity that
carries it is the matrix between the pins.


## 3. The same candidates through the loss model

`examples/candidate_loss_map.py` puts each candidate's margin, dipole gain and axial switch strain
into the calibrated literature channel (the 63/30/7 terpolymer's Vogel-Fulcher relaxation,
`tau0 = 8.08e-13 s`, `Ea = 0.124 eV`, `Tf = 191 K`, dielectric strength 45, halved under 60 MV/m;
plan section 11) and evaluates the program's operating point. The relaxation is the same for every
row, so the table separates what the chemistry changes in this model (the polar fraction at bias,
the strain per field, the heat at a given strain) from what it does not yet (the loss per unit
strain, which is `omega tau` and is the bench's number until tier 3 computes a barrier per
chemistry).

| sequence | margin kcal/mol per monomer | E* MV/m | polar fraction at 55 V/um | S_ac at 1.8 V/um, 34 kHz | tan delta 34 / 66 kHz | heat MW/m^3 34 / 66 kHz | efficiency 34 / 66 kHz | stroke per lever 34 / 66 kHz |
|---|---|---|---|---|---|---|---|---|
| pvdf | +0.744 | 513 | 0.00 | 0.000% | 0.000 / 0.000 | 0.0 / 0.0 | 0.00 / 0.00 | 0 / 0 |
| vdf7-trfe3 | +0.145 | 112 | 0.11 | 0.068% | 0.051 / 0.098 | 5.3 / 18.8 | 0.13 / 0.07 | 23 / 45 |
| vdf8-trfe2 | +0.177 | 133 | 0.05 | 0.034% | 0.032 / 0.061 | 1.9 / 6.8 | 0.10 / 0.05 | 12 / 22 |
| vdf6-trfe4 | -0.329 | 0 | 0.88 | 0.071% | 0.054 / 0.104 | 5.8 / 20.8 | 0.13 / 0.07 | 24 / 47 |
| vdf5-trfe5 | +0.752 | 607 | 0.00 | 0.000% | 0.000 / 0.000 | 0.0 / 0.0 | 0.00 / 0.00 | 0 / 0 |
| vdf7-trfe2-cfe1 | +0.966 | 740 | 0.00 | 0.000% | 0.000 / 0.000 | 0.0 / 0.0 | 0.00 / 0.00 | 0 / 0 |
| vdf7-trfe2-cfe1-adj | +0.961 | 736 | 0.00 | 0.000% | 0.000 / 0.000 | 0.0 / 0.0 | 0.00 / 0.00 | 0 / 0 |
| vdf7-trfe2-ctfe1 | +0.790 | 676 | 0.00 | 0.000% | 0.000 / 0.000 | 0.0 / 0.0 | 0.00 / 0.00 | 0 / 0 |
| vdf7-trfe2-cdfe1 | +0.146 | 119 | 0.10 | 0.062% | 0.047 / 0.091 | 4.1 / 14.7 | 0.14 / 0.08 | 21 / 41 |
| vdf9-cfe1 | +0.982 | 693 | 0.00 | 0.000% | 0.000 / 0.000 | 0.0 / 0.0 | 0.00 / 0.00 | 0 / 0 |
| vdf8-cfe2 | +0.835 | 675 | 0.00 | 0.000% | 0.000 / 0.000 | 0.0 / 0.0 | 0.00 / 0.00 | 0 / 0 |
| vdf7-trfe2-an1 | +0.568 | 441 | 0.00 | 0.000% | 0.000 / 0.000 | 0.0 / 0.0 | 0.00 / 0.00 | 0 / 0 |
| vdf7-trfe2-fan1 | +0.912 | 681 | 0.00 | 0.000% | 0.000 / 0.000 | 0.0 / 0.0 | 0.00 / 0.00 | 0 / 0 |
| vdf9-an1 | +0.470 | 345 | 0.00 | 0.000% | 0.000 / 0.000 | 0.0 / 0.0 | 0.00 / 0.00 | 0 / 0 |
| vdf6-trfe3-cdfe1 | +0.805 | 668 | 0.00 | 0.000% | 0.000 / 0.000 | 0.0 / 0.0 | 0.00 / 0.00 | 0 / 0 |
| vdf6-trfe3-cfe1 | +0.297 | 235 | 0.00 | 0.001% | 0.001 / 0.001 | 0.0 / 0.0 | 0.02 / 0.01 | 0 / 1 |
| vdf9-trfe4-cfe1-14 | +0.110 | 87 | 0.24 | 0.125% | 0.074 / 0.143 | 13.7 / 47.8 | 0.16 / 0.10 | 42 / 82 |
| vdf10-trfe4-14 | -0.051 | 0 | 0.88 | 0.071% | 0.052 / 0.102 | 5.8 / 20.6 | 0.13 / 0.07 | 24 / 46 |
| vdf9-trfe5-14 | -0.058 | 0 | 0.88 | 0.072% | 0.053 / 0.103 | 5.8 / 20.8 | 0.13 / 0.08 | 24 / 47 |
| vdf9-trfe4-cdfe1-14 | +0.402 | 329 | 0.00 | 0.000% | 0.000 / 0.000 | 0.0 / 0.0 | 0.00 / 0.00 | 0 / 0 |
| vdf8-trfe4-cfe1-cdfe1-14 | +0.641 | 529 | 0.00 | 0.000% | 0.000 / 0.000 | 0.0 / 0.0 | 0.00 / 0.00 | 0 / 0 |

**Reading the mapping.** The operating point is 55 V/um bias, 1.8 V/um AC; the polar fraction
at bias is where the crystal margin puts the channel on its own switching curve, and the strain
per field is largest where that fraction is near one half: with `dmu` about 0.57 e.A the
margin that does it at 55 V/um is about +0.07 kcal/mol per monomer. After the handedness
correction no ten-monomer row sits there: the no-defect 70/30 and 80/20 and the CDFE row sit low
on the curve (polar fraction 0.05-0.11) and give 0.03-0.07% AC strain, the TT-ground 60/40 is
saturated (0.88) and gives the same, and every pinned defect row (CFE, CTFE, AN, FAN, margins
+0.57 to +0.97), PVDF, the alternating 50/50 and the VDF/CFE copolymers do not switch at all.
The rows that do sit on the curve are the fourteen-monomer cells of section 2: the 64/29/7 cell
(margin +0.11, polar fraction 0.24 at 55 V/um) carries the largest AC strain of the table (0.125%
at 34 kHz, tan delta 0.074, efficiency 0.16) and its 71/29 matrix control (margin -0.05) is
saturated (0.88) with half that strain, which is the pin doing its job: it holds the matrix off
saturation at the bias. Section 8 carries these two rows with their own kink energies. The loss per
unit strain is the same for every switchable row (efficiency 0.10-0.15 at 34 kHz, 0.05-0.09 at
66 kHz) because the relaxation is the same: a candidate with more strain at this bias also has
more loss and more heat, and all of them exceed the 4 MW/m^3 budget at 1.8 V/um at 34 kHz (the
level-1 and level-2 result of plan sections 11-12), so the AC amplitude, not the chemistry, sets
the heat in this picture. **The chemistry lever on loss per strain is `tau`, and this pass does
not compute it.** The literature channel's `f_relax` is 360 kHz at 300 K, a factor 5-10 above
the drive, so `tan delta` is about `omega tau` there (0.05-0.17 across the rows); halving it at
the same strain means doubling `f_relax`, i.e. a lower `Tf` (smaller cooperative units: more
defects, shorter trans runs) or a lower local barrier, and which defect does that is tier 3's
question (the E3 kink/flip barrier asked of Dipole, the bench's Vogel-Fulcher fit per film).

## 4. Defect type, concentration and placement: what the model says and what it cannot

**Type.** At one defect per ten monomers beside 20 mol% TrFE, every defect but CDFE pins its
neighbourhood in the non-polar helix. Against the 80/20 copolymer cell it replaces one VDF of
(TTTT -8.165, TG -8.341 kcal/mol per monomer), the cost of the defect in each conformation, per
defect unit, and the pin strength it leaves (TT cost minus TG cost: how much the unit prefers the
helix) are:

| defect, one per ten beside 20 mol% TrFE | TT cost per unit, kcal/mol | TG cost per unit | pin strength (TT - TG) | margin TT-TG per monomer |
|---|---|---|---|---|
| CFE (CH2-CFCl), two from either TrFE | +6.9 | -1.0 | +7.9 | +0.97 |
| CFE next to a TrFE | +7.9 | 0.0 | +7.9 | +0.96 |
| CTFE (CF2-CFCl) | +21.9 | +15.8 | +6.1 | +0.79 |
| CDFE (CHCl-CF2) | +17.2 | +17.6 | -0.3 (other handedness not buildable) | +0.15 |
| AN (CH2-CH(CN)) | -0.7 | -4.6 | +3.9 | +0.57 |
| FAN (CH2-CF(CN)) | +5.9 | -1.5 | +7.4 | +0.91 |

In this cell CFE, FAN and CTFE pin at 6-8 kcal/mol per unit, AN at 4, CDFE not at all: the chlorine
on the hydrogen-bearing carbon costs the all-trans and the helical lattice the same 17 kcal/mol, a
large unit that fits neither, where CFE's chlorine on the fluorine-bearing carbon costs TTTT 7 and
the helix nothing. These numbers are the unit's in a 20 mol% TrFE cell with a pin every ten; in the
fourteen-monomer 29 mol% cell the same CFE pins at +2.3 and CDFE at +6.3 (section 2), so a pin
strength is a property of the unit in its matrix and the fourteen-monomer cells are the ones to
design on. CTFE's two costs are both large (it is the bulkiest unit) and its pin is the
difference of two big numbers. The nitrile pins by packing *better* into the helix (TG cost -4.6)
rather than worse into the all-trans lattice (-0.7), the opposite route to the chlorines'. These
are the rows whose chemistry is extrapolated; the pin strength of CFE against CDFE is the
quantity most worth a DFT cell (two 120-atom cells each).

**Concentration.** The two levers are different in kind, as the first pass said, but the
division of labour is now the opposite of what it claimed: TrFE sets the margin of the matrix
(the near-degeneracy at 30-40 mol%), the pin sets how much of the chain is held in the helix and
how the polar order is broken into the relaxor's nanodomains, and the pin density decides whether
any matrix is left between pins: one CFE per ten in a 20 mol% matrix pins the whole cell (+0.97);
one per fourteen in a 29 mol% matrix leaves it at +0.11. CFE without TrFE does not do TrFE's job
(one or two per ten VDF leave the helix lowest by 0.84-0.98 with all four conformers run). The terpolymer recipe the literature converged on (about 30 mol% TrFE, 7-10 mol% CFE) is
what the fourteen-monomer cell reproduces; the ten-monomer 70/20/10 cell does not, and the margin
near zero it showed in the first pass was the handedness artefact, not the recipe.

**Placement.** Moving the single CFE from two monomers away from either TrFE to the position
next to one changes the sequence's margin by 0.005 kcal/mol per monomer: the TT cost of the unit
rises by 1.0 kcal/mol and its TG cost by the same, so the pin strength and the margin do not move.
The local-margin spread a statistical terpolymer carries from random placement is then not the
spread of the pin's position relative to TrFE (0.05 kcal/mol per unit) but the spread of the
distances between pins, i.e. of the matrix segment lengths, which the fourteen-against-ten
comparison puts at about 0.9 kcal/mol per monomer between a segment of three and one of six. That
is a large spread, the one the driven lattice of plan section 12 found to suppress the stroke at
bias and broaden the loss peak, and it is exactly what a controlled sequence sets by choice: a
strictly periodic pin spacing makes it zero and the switch as sharp as the chain's own
cooperativity allows, at the risk of ordering ferroelectrically when the coupling exceeds about
0.03 eV. The model cannot say which spread is optimal; the lattice model can, once the bench's
`S(E_dc)` curve fixes its coupling.

**Location in the larger sense.** The screen puts every defect inside a crystalline chain. The
literature's finding is that the giant strain is produced at and around the crystals, in the
oriented amorphous fraction and the tie molecules, where the margin is lower than the crystal's
and the barrier is the local one. The model has no row for a defect at a lamellar surface or in a
tie segment; that is tier 4, and for the controlled-sequence route it is the question of where
along a chain of a given length the defects go so that the chain folds with them at the surface.

## 5. Candidates, ranked for the sequence-controlled route

For a sequence-controlled material aimed at the bias-hold operating point, on this model and
with the confidence stated in section 7 (and after the handedness correction, which withdrew the
first pass's ranking):

1. **The matrix: VDF/TrFE at 29-36 mol% TrFE.** Degenerate within the model's resolution in the
   fourteen-monomer cells (71/29 -0.05, 64/36 -0.06; the ten-monomer 70/30 and 80/20 cells +0.15
   and +0.18), validated chemistry: this is
   where the near-degeneracy the relaxor needs comes from, as it does in every giant-strain
   terpolymer made so far. Without a pin it is a ferroelectric copolymer with hysteresis, which
   this equilibrium picture does not contain (the lattice model's coupling does).
2. **The pin: CFE at one per fourteen (7 mol%).** A TG pin (+7.9 kcal/mol per unit in the 20 mol%
   cell, +2.3 in the 29 mol% one) that holds the matrix off saturation at the bias (polar fraction
   0.24 against the pinless matrix's 0.88) and loosens the lattice for the kinks that switching
   needs (kink pair 5.5 kcal/mol in the 71/29 matrix, 1.8 with the pin: section 8), so the cell both
   switches at 55 V/um and follows a 34-66 kHz drive, with tan delta 0.07 for 0.125% AC strain in
   the mapping. Its placement relative to TrFE is not a lever for the
   margin (0.005 per monomer); the pin spacing is (section 4), and a controlled sequence fixes it.
   Chlorine chemistry: extrapolated until a DFT cell tests it.
3. **Alternative pins.** FAN (+7.4) and CTFE (+6.1) pin as CFE does; AN at half the strength
   (+3.9) would need twice the density for the same job, with the polarization cost that brings
   (section 6). CDFE does not pin in the ten-monomer 20 mol% cell (-0.3) but tips the fourteen-
   monomer 29 mol% cell by +0.45 per monomer, more than CFE (+0.16), leaving it at +0.40: at the
   film's composition it is a stronger pin than CFE and overshoots the steep part of the curve;
   its kink pair is the cheapest of the table (section 8). The first pass's CDFE rows were read
   wrongly on both points.
4. **A hypothesis tested and failed as posed:** that the pin (CFE) and a lattice loosener (CDFE,
   neutral in the ten-monomer cell with a free kink pair) do different jobs, so one of each per
   fourteen would hold the relaxor's disorder with the pin and buy its dynamics with the loosener.
   In the 29 mol% fourteen-monomer cell CDFE is not neutral (+0.45 per monomer, more than CFE's
   +0.16) and the two add: the CFE+CDFE cell sits at +0.64, off the curve, and CDFE alone at +0.40.
   What the idea needs is a unit that costs the all-trans and the helical lattice alike *at the
   film's dilution*, which neither chlorine is; the loosening both chlorines bring (kink pair 1.1-1.8
   against the matrix's 5.5) comes with the pinning, and at one per fourteen CFE's is the usable
   combination.
5. **Out:** PVDF and the alternating 50/50 (not switchable in the crystal), CFE without TrFE (the
   helix by 0.8-1.0), and the ten-monomer 70/20/10 cells as models of the switching unit.

Section 8 adds the frequency reading, which does not depend on the margin: a kink pair in the
all-trans lattice costs 4.3-9.3 kcal/mol in PVDF and 4.9-7.6 in the copolymers but 1.2-3.1 in the
CFE cell and 0.1-0.6 in the CDFE cell, so the chlorine loosens the lattice for the kinks that
switching needs, by a factor 200-1600 in rate. The operating-point loss and strain, which do
depend on the margin, now belong to the fourteen-monomer cell: its kink energy is the next
computation. What would change the ranking: the migration barrier and the cooperativity (not
computed); Dipole's E2 cells; and the bench's Vogel-Fulcher and bias-suppression numbers, which
replace the literature channel the mappings are built on.

## 6. The slate's other chemistries

Dipole's slate (`sarco/materials/computational/candidate_slate_v1.md`, `_v2.md`; ranking
`ranking_v3_correction_v1.md`, 2026-09-02; the PI's reframe `objective_function_v1.md`; the
registry gate `registry_gate_v1.md`; the barrier table `field_barrier_table_v1.json`) is a
single-chain conformer screen plus GFN1 packing, scored by a frozen figure of merit
`FOM = eps_s x (dipole/volume) x frustration / (B_mean x (1 + D_mean))`, a strain-per-loss proxy.
Its v3 ranking, relative to the FA tetrapolymer: VDCN 9.3, cyano-epoxide (CNEPO) 6.5, FANOMe 5.5,
thiocyanate 4.2-5.6, VDC 5.3, dicyano-alkene 5.3, VNC 3.5, adjFAN 3.0, AN 2.0, HFP 1.5, CDFE 1.2.
Three of Dipole's own later findings bear on it, and the PI's reframe (minimize tan delta subject
to floors on stroke rate, frequency response, crystallinity and radial stroke; "the FOM leader is
not automatically the recommendation") is the lens:

* **The registry gate.** Whether the polar, in-register VDF lattice is still the ground state with
  the defect in it, at 8.3 mol%: PVDF +0.56, CFE +0.24, CNEPO +0.31 (ground state), CDFE -0.06 and
  AN -0.11 (marginal), VDC -0.38, FANOMe -1.13 and VDCN -1.18 kcal/mol per monomer (the slipped
  lattice preferred by more than a kcal). The FOM's first and third place replace the polar
  crystal rather than frustrate it; Dipole's recommendation accordingly reordered to CNEPO, then
  AN and CDFE, with VDCN and FANOMe demoted.
* **The barrier table** is single-bond torsional barriers, not switching rates. Dipole's
  `field_barrier_v1.md` scores each pentad's defect bonds (approach, entry, exit) and a host
  reference bond for the trans-gauche rotation barrier B and the well difference D, at GFN2 and
  re-scored at DFT, and folds them into the FOM's loss term `L = B_mean (1 + D_mean)` over the
  defect bonds. At DFT: baseline terpolymer L 37.4, PVDF 27.3, CNEPO 16.7, AN 14.2, VDCN 14.5, so
  on this term every slate candidate scores 2-2.6 times better than the baseline; read bond by
  bond the picture is mixed (CNEPO's host-reference bond is 16.4 against PVDF's 7.4 while its
  flanking exit bond is 4.3; AN's host bond 10.2, its approach bond 7.6). Dipole's own conclusion
  from the same table is the one that matters here: the field moves these barriers by 1-4% at
  100 MV/m, so a single monomer cannot be switched by the field, switching is collective, and the
  single-chain barrier does not decide the loss. The quantity that does, the collective switching
  barrier of a trans run inside the lattice, has not been computed per chemistry by either seat
  (polyfind's one number is the kink-row hop in perfect beta, 0.357 eV); section 8 below starts on
  it with the kink-pair formation energy in the packed cell.
* **Orientation.** Mirroring an asymmetric defect's substituents in an otherwise unchanged
  lattice moves the packing energy by 0.3-0.5 kcal/mol per monomer (AN 0.53, CFE 0.35; SCN 1.7;
  the symmetric VDCN 0.006 as the control). That is the same statement section 2 makes from the
  alternating 50/50 row: a regular sequence's handedness pattern is a packing lever of the size of
  a composition change.

**What this seat adds: the crystal switch.** The FOM has no crystal switching margin in it; the
screen of section 2 is exactly that quantity. `examples/slate_screen.py` runs the slate's
polymers that polyfind can build through the same TTTT-against-TGTG' screen, under the plain
valence preset `pvdf-dft-valence` (no belt sites exist for these chemistries, so the fitted
crystal terms of belt4 do not apply; the PVDF control under the same preset shows the cost:
margin +1.00 against belt4's +0.74, P 0.140 against 0.193):

| polymer (valence preset, plain chains) | margin TT-TG, kcal/mol per monomer | switch strain axial, areal | P(TTTT) C/m^2 | threshold E*, MV/m |
|---|---|---|---|---|
| PVDF control | +1.00 | +10.6%, -14.6% | 0.140 | 990 |
| polyacrylonitrile (AN) | +5.76 | +12.4%, -27.7% | 0.043 | 15,800 |
| poly(vinylidene cyanide) (VDCN) | +9.31 | +14.3%, -19.2% | 0.113 | 7,400 |
| VDF/VDCN 11:1 (8.3 mol% VDCN) | +0.31 with tilts, -1.58 rigid: undetermined (its TGTG' cell releases 1.9 kcal/mol per monomer through pendant tilts, ten times a PVDF-family cell's, so the rigid geometry is not a reference) | +14.9%, -9.2% | 0.098 | 320, not meaningful |
| PVDC | not buildable: its 123 deg CH2 angle makes the all-trans chain non-periodic (its crystal chain is not planar) |

The nitrile homopolymers' all-trans crystals sit 6-9 kcal/mol per monomer above their TGTG' forms:
nitriles on alternate carbons cannot coexist in a planar zigzag (the strain finding of
docs/VALENCE_FIT.md section 4, with the backbone angles now free to relax), so as conformational
switches they are out by two orders of magnitude in threshold. That agrees with what is known
(polyacrylonitrile is no ferroelectric crystal; the VDCN copolymers are glassy dipolar
piezoelectrics, not switching crystals) and with Dipole's registry gate from the other side: a
unit that cannot sit in a polar lattice neither switches with it nor frustrates it, it replaces
it. As a dilute defect in a VDF/TrFE host the question is different and is the one Dipole's
recommendation asks (AN as "the makeable track"); the belt4 screen's rows for it:

| sequence (belt4 preset, as section 2) | margin TT-TG (rigid), kcal/mol per monomer | TGTG' tilt release | switch strain axial, areal | P(TTTT) C/m^2 | mu(TTTT) e.A | E* |
|---|---|---|---|---|---|---|
| VDF/TrFE/AN 70/20/10 | +0.57 (+0.33), TG'TG lowest; TGTG' alone gave -0.73 | -0.29 | +11.6%, -16.9% | 0.146 | 0.558 | 440 MV/m |
| VDF/TrFE/FAN 70/20/10 | +0.91 (+0.91; tilts not treated: the FAN carbon's pendants are a fluorine and a nitrile fragment), TG'TG lowest; TGTG' alone gave -0.57 | - | +11.0%, -12.9% | 0.144 | 0.581 | 680 MV/m |
| VDF/TrFE/CFE 70/20/10, for comparison | +0.97 (+0.76), TG'TG lowest; TGTG' alone gave +0.03 | -0.33 | +11.5%, -14.5% | 0.159 | 0.566 | 740 MV/m |

A nitrile defect beside 20 mol% TrFE is a TG-favouring pin of about half CFE's strength (section
4: +3.9 against +7.9 kcal/mol per unit), and it pins by the opposite route: where CFE's chlorine
costs the all-trans lattice 7 kcal/mol per unit and the helix nothing, the nitrile costs the
all-trans lattice nothing (-0.7) and packs 4.6 better into the helix (the TG'TG cell opens to 63-67
A^3 per monomer around it). With the lowest handedness its margin is +0.57; the first pass's TGTG'
alone gave -0.73 (TTTT the ground state) and the reading built on it, that a nitrile drives the
TrFE host to a saturated ferroelectric, is withdrawn. So a nitrile is a weaker relaxor pin than
CFE, needing about twice the density for the same hold, at a polarization beside TrFE of 0.146
C/m^2 against CFE's 0.159. In a VDF host without TrFE the VDF/AN 90/10 row stands with all four conformers run: one
acrylonitrile per ten VDF takes PVDF's +0.74 to +0.47 (threshold 510 to 840 MV/m on the TT-TG
dipole difference, 335 on the TTTT dipole; TTTT polarization 0.193 to 0.163 C/m^2), so a nitrile
is a weak trans-favouring unit in a VDF host and a weak helix-favouring pin beside TrFE: it
follows the lattice it is put in rather than setting it. Dipole's registry gate reads AN as
"marginal" (-0.11), a different quantity saying something compatible from the lattice's side: the
nitrile does not frustrate the polar crystal strongly. Its single-bond barriers do not settle its
speed either way (host bond 10.2 against CFE's 6.8, defect bonds' mean 9.6 against 12.2; see the
barrier paragraph above).

**What this seat cannot compute.** CNEPO, SCN, SF5, nitro, isocyanide, the cyano-alkenes and
cyclopropyl need sulfur, nitro, epoxide and sp2 backbone parameters polyfind does not have; for
those Dipole's single-chain and GFN1 numbers are the only ones, and the barrier column above is
the one that matters for this programme's objective.

## 8. Frequency and loss: the kink-pair formation energy

The tables above screen the statics. The objective is loss at the drive frequency, which is
`omega tau`, and `tau` is dynamics: the collective barrier to nucleate and move a kink pair in a
trans run inside the lattice, and the cooperativity behind `Tf`. Neither seat had computed any
of it per chemistry (section 6). `examples/kink_formation_screen.py` computes the first piece: the
energy of a kink pair (T G+ T G- T, one per chain, at bonds 0-2 of the twenty-bond repeat, in the
VDF run away from the TrFE and the defect) embedded in the all-trans two-chain crystal, with the
pendant tilts, against the all-trans cell, under the packing protocols below. Switching is nucleation-limited
once the single-bond barrier is ruled out as the rate (Dipole's own conclusion), so
`exp(-E_form / kT)` is the first-order rate factor between chemistries; the migration barrier is
not computed and is taken equal for every row. `examples/candidate_frequency_map.py` then gives
each candidate its own relaxation rate: the CFE 70/20/10 row is anchored on the measured
Vogel-Fulcher relaxation of the 63/30/7 film (364 kHz at 300 K) and every other row's rate is
scaled by its kink factor, with the screen's margin, dipole and switch strain as before, at the
program's operating point.

The kinked cell's packing is the protocol question, and the protocols bracket the number. The
**own basin** protocol packs the kinked chain by its own screen, refines it and relaxes the cell:
the fully relaxed kink crystal, every chain kinked in register and the cell following, the lowest
minimum the screen finds. The **same-basin start** (`--same-basin`) begins that refinement from the
all-trans cell's packing parameters, which is a start and not a constraint (the cell moves anyway;
kept for the record). The **pinned lattice** (`--pinned-lattice`, `refine_crystal(cell_free=False)`
then an internal-coordinate relaxation) holds `a, b, gamma` at the all-trans values and relaxes the
chains' placement, torsions and bond angles inside them: the lattice the surrounding crystal pins,
the upper bound for a dilute kink pair. **Pinned then released** (`--pinned-release`) lets the cell
go from the pinned solution by a local relaxation: the kink crystal in the basin nearest the
all-trans cell. Results `docs/predictions/kink_formation_belt4.json`, `_samebasin.json`,
`_pinned.json`, `_pinned_released.json`; kcal/mol per kink pair. The bracket's width is the lateral
relaxation the kink needs, which is a property of the chemistry's all-trans lattice in its own right.

| sequence | own basin (relaxed kink crystal) | same-basin start | pinned lattice (upper bound) | pinned then released | bracket, kcal/mol per kink pair |
|---|---|---|---|---|---|
| vdf10 | +4.26 | +4.59 | +9.31 | +6.80 | 4.3-9.3; the released cell (4.87 x 8.83 A) is the own-basin cell, 2.5 higher: the release moves the cell and the placement but not the torsions, so it is itself an upper bound on the nearest-basin kink crystal |
| vdf8-trfe2 | +4.92 | +6.89 | +7.58 | +5.16 | 4.9-7.6; released and own basin agree (5.2 against 4.9) |
| vdf7-trfe3 | +4.95 | +5.09 | +7.48 | +5.12 | 4.9-7.5; released and own basin agree (5.1 against 5.0) |
| vdf7-trfe2-cfe1 | +1.24 | +1.47 | +3.13 | +1.50 | 1.2-3.1; released and own basin agree (1.5 against 1.2) |
| vdf7-trfe2-cfe1-adj | +1.40 | +2.97 | +3.14 | +1.50 | 1.4-3.1; the placement does not change the kink cost (released 1.50 against 1.50) |
| vdf7-trfe2-ctfe1 | -0.42 | +5.68 | +20.42 | +13.63 | 13.6-20.4 in the all-trans cell's own basin; the own-basin -0.4 is a different crystal (5.10 x 9.50 A against the all-trans 5.69 x 8.13), so the all-trans CTFE lattice is metastable against a kinked packing and has no room for a kink itself |
| vdf7-trfe2-cdfe1 | +0.61 | +0.72 | +0.40 | +0.08 | 0.1-0.6, the narrowest: the kink pair is free within the protocols' resolution; the kink fits the all-trans CDFE lattice as it is (65.6 A^3 per monomer, gamma 96.7 deg) |
| vdf7-trfe2-an1 | +1.27 | +3.56 | +11.28 | +6.76 | 6.8-11.3 in the all-trans cell's basin (the own-basin 1.3 is another cell): the nitrile fills the interchain gap (section 6) and locks the lattice against kinks |
| vdf9-an1 | +0.38 | +1.35 | +12.17 | +7.11 | 7.1-12.2 in the all-trans cell's basin (the own-basin 0.4 is another cell, 5.15 x 9.58 A), as above |
| vdf9-trfe4-cfe1-14 (64/29/7, fourteen monomers) | +1.82 | - | +13.21 | +11.75 | 1.8 in its own basin; the pinned and released values are a stuck placement rather than a bound (the ten-monomer CFE cell's pinned kink is 3.1 in the same 57 A^3 lattice; releasing the cell recovers only 1.5), so the pinned lattice is being reached from the kink's own conformation instead (`--pinned-from-own`) |
| vdf10-trfe4-14 (71/29 matrix, fourteen monomers) | +5.47 | - | +7.27 | running | 5.5-7.3: the matrix without the pin has the copolymer's kink cost; the one CFE per fourteen cuts it to 1.8 |

The frequency map (`examples/candidate_frequency_map.py`) is run on both ends of the bracket, the
CFE 70/20/10 row anchored at the measured 364 kHz in each; the full tables are
`candidate_frequency_map_belt4.md` (own basin), `_samebasin.md`, `_pinned.md` and
`_pinned_released.md`. The statics columns carry the corrected margins (section 2): with them no
ten-monomer defect cell switches at 55 V/um, so their operating-point columns are empty and the
rate column is what this table says about them. The cell that carries the operating point is the
fourteen-monomer 64/29/7 cell, whose kink energy is the next computation; its row goes here when
it lands.

| sequence | f_relax proxy, lower / upper end | margin (corrected) | polar fraction at bias | S_ac at 34 kHz, lower / upper | tan delta at 34 kHz, lower / upper | tan delta per % strain, lower / upper | reading |
|---|---|---|---|---|---|---|---|
| vdf10 | 2.3 kHz / 12 Hz | +0.74 | 0.00 | 0.000% / 0.000% | 0.000 / 0.000 | - / - | not switchable at the bias (margin +0.74) |
| vdf8-trfe2 | 750 Hz / 210 Hz | +0.18 | 0.05 | 0.002% / 0.000% | 0.176 / 0.050 | 102.13 / 103.10 | frozen at the drive at both ends: the loss tail of a channel below 1 kHz, no strain |
| vdf7-trfe3 | 721 Hz / 249 Hz | +0.14 | 0.11 | 0.002% / 0.001% | 0.228 / 0.079 | 98.02 / 98.75 | the same |
| vdf7-trfe2-cfe1 | 364 kHz / 364 kHz | +0.97 | 0.00 | 0.000% / 0.000% | 0.000 / 0.000 | - / - | the anchor of the rate scale (364 kHz measured on the 63/30/7 film); with the corrected margin (+0.97) the ten-monomer cell itself does not switch at the bias |
| vdf7-trfe2-cfe1-adj | 276 kHz / 359 kHz | +0.96 | 0.00 | 0.000% / 0.000% | 0.000 / 0.000 | - / - | as above (margin +0.96) |
| vdf7-trfe2-ctfe1 | 5.8 MHz / frozen | +0.79 | 0.00 | 0.000% / 0.000% | 0.000 / 0.000 | - / - | pinned (margin +0.79) and, in the all-trans cell's basin, no room for a kink |
| vdf7-trfe2-cdfe1 | 1.0 MHz / 35.4 MHz | +0.15 | 0.10 | 0.062% / 0.062% | 0.016 / 0.000 | 0.27 / 0.01 | above the drive at both ends with the lowest loss per strain of the table; margin +0.15, polar fraction 0.10 |
| vdf7-trfe2-an1 | 343 kHz / frozen | +0.57 | 0.00 | 0.000% / 0.000% | 0.000 / 0.000 | - / - | pinned (margin +0.57); the kink is 6.8-11.3 in the all-trans cell's basin |
| vdf9-an1 | 1.5 MHz / frozen | +0.46 | 0.00 | 0.000% / 0.000% | 0.000 / 0.000 | - / - | not switchable at the bias (margin +0.46; other conformers running) |

**Anchored on the fourteen-monomer cell** (`candidate_frequency_map_belt4_anchor14.md`: the 64/29/7
cell takes the measured 364 kHz, every other row's rate is relative to its own-basin kink energy):

| sequence | kink pair (own basin) | f_relax proxy | margin | polar fraction at bias | S_ac at 34 kHz | tan delta 34 / 66 kHz | tan delta per % strain | efficiency 34 / 66 kHz |
|---|---|---|---|---|---|---|---|---|
| 64/29/7, one CFE per fourteen | +1.82 | 364 kHz (anchor) | +0.11 | 0.24 | 0.125% | 0.074 / 0.143 | 0.59 | 0.16 / 0.10 |
| 71/29, the matrix without the pin | +5.47 | 0.8 kHz | -0.05 | 0.88 | 0.003% | 0.262 / 0.136 | 101 | 0.00 / 0.00 |
| 64/29/7 with CDFE in place of CFE | +1.08 | 1.3 MHz | +0.40 | 0.00 | 0.000% | 0.000 / 0.000 | - | 0.02 / 0.01 |
| VDF/TrFE/CDFE 70/20/10, ten monomers | +0.61 | 2.8 MHz | +0.15 | 0.10 | 0.062% | 0.006 / 0.012 | 0.10 | 0.42 / 0.32 |
| VDF/TrFE 70/30, ten monomers | +4.95 | 2 kHz | +0.15 | 0.11 | 0.006% | 0.576 / 0.310 | 93 | 0.00 / 0.00 |

The matrix alone is a frozen ferroelectric at the drive: saturated at the bias (0.88) with its
channel at 0.8 kHz, so it carries loss and no strain at 34-66 kHz, the copolymer's known
behaviour. One CFE per fourteen does two things at once: it holds the matrix off saturation (polar
fraction 0.24) and it cuts the kink cost from 5.5 to 1.8 kcal/mol, a factor of 460 in rate, which
puts the channel at the drive. That is the relaxor terpolymer's mechanism in one cell, pin and
loosener in the same chlorine, and at the program's operating point the model gives tan delta 0.07
at 34 kHz for 0.125% AC strain, the largest strain of the table at the lowest loss per strain of
any row that switches. The ten-monomer CDFE row shows what a loosener without a pin does: a
cheaper kink still (0.6), a channel at 2.8 MHz far above the drive, a lower loss per strain (0.10)
but half the strain, because nothing holds its matrix on the steep part of the curve. CDFE at one
per fourteen in the 29 mol% matrix is the fastest lattice of the table (kink pair 1.1, channel 1.3
MHz) and switches nothing at the bias: at the film's dilution it pins harder than CFE (+0.40) and
overshoots the steep part, so the loosening it buys is not usable on its own. The cell with one
CFE and one CDFE per fourteen is pinned harder still (+0.64; section 5, item 4), its kink energy is
running.

**Reading.**

* **The relaxor mechanism, from inside the lattice.** A kink pair in a VDF run costs 4.3-9.3
  kcal/mol in PVDF and 4.9-7.6 in the 80/20 and 70/30 copolymers (TrFE stiffens the trans lattice
  against kinks), but 1.2-3.1 in the CFE terpolymer model and 0.4-0.6 with CDFE, although the kink
  sits in a VDF run away from the defect: one chlorine per ten monomers loosens the whole lattice
  (the model cell carries the defects in register) and cuts the nucleation cost by 3.7 kcal/mol at
  the lower ends of the brackets and 4.4 at the upper, a factor of 500-1600 in rate at 300 K. That
  separation does not depend on the protocol, and it is why the terpolymer switches at tens of
  kilohertz and the copolymer does not: at either end the copolymers' channel sits below 1 kHz
  against the terpolymer's 364 kHz, so at 34-66 kHz it is frozen, with the loss tail of a peak far
  below the drive (tan delta 0.05-0.23 for no strain). The kink energy is a property of the
  lattice around the kink and does not depend on the margin, so this reading survived the
  handedness correction that emptied the operating-point columns of the ten-monomer defect cells:
  a row needs both a near-zero margin (the matrix's, section 2) and a low kink energy (the
  chlorine's) to switch at the bias and follow the drive, and the two are set by different units.
* **Which end applies where** is an interpretation, not a computation: a relaxor operates in a
  lattice its own disorder has already expanded, so the terpolymer rows' operating number sits
  nearer the lower end; a ferroelectric copolymer's first kink nucleates in a pinned all-trans
  lattice, so its number sits nearer the upper end. Both readings widen the separation. The
  calculation that closes the bracket is a dilute kink: one kinked chain among all-trans
  neighbours that are free to displace, which needs a mixed-conformation cell the packer does not
  have (every chain in a cell is a copy of one periodic chain today).
* **The loss-per-strain column is the one to rank on**, and before the correction it put CDFE
  (0.01-0.27) ahead of CFE (0.58-0.76) at both ends with the no-defect copolymers out (loss without
  strain); after it the CFE cells have no strain at the bias (they are pinned, section 4) and CDFE
  keeps its 0.07-0.27, so the comparison now has to be made in the fourteen-monomer cell, where the
  CFE pin is diluted and the matrix switches.
  The CTFE row's rank is withdrawn: its kink energy runs from -0.4 (free, the lattice relaxes to a
  kinked packing 5.10 x 9.50 A from the all-trans 5.69 x 8.13) to +20.4 (pinned, no room for a
  kink), and no rate follows from that. The nitrile rows are undetermined for the same reason
  (1.3-11.3 beside TrFE, 0.4-12.2 in VDF/AN 90/10): the nitrile fills the interchain gap of the
  all-trans lattice (section 6) and locks it against kinks when it is pinned. The CFE placement
  does not change the kink cost (1.24 against 1.40 free, 3.13 against 3.14 pinned); the adjacent
  placement's higher margin (+0.06 against +0.03) gives it more strain at the bias and a lower loss
  per strain at the upper end.
* **Two cautions on the bracket itself.** The own-basin column is the lowest minimum the packing
  screen found, not a proven bound: for CDFE the pinned solution (0.40) and the released one (0.08)
  lie below it (0.61), so the screen missed kink-crystal minima there, and for the CTFE and nitrile
  rows it found a different crystal rather than the all-trans cell's basin (their released values,
  6.8-13.6, are the ones to read against the all-trans lattice). The released column agrees with
  the own basin for the copolymers and the CFE rows (within 0.3) and sits 2.5 above it for PVDF,
  whose release did not re-refine the torsions. And the same-basin start showed that a start is
  not a constraint: four of its rows moved
  to other cells anyway (PVDF to the a/b-swapped setting, CDFE to 5.10 x 10.17 A, AN beside TrFE
  to 5.19 x 9.49 A, CTFE to a cell 13% larger), and where the cell stayed the restarted refinement
  sat 0.1-2.0 kcal/mol per kink pair above the own-basin minimum.
* **The kink pair has a handedness too.** G+ T G- and G- T G+ are different conformers in a chiral
  sequence, as the helix's handednesses are, and the screen ran the first. The mirrored pair in
  the pinned lattice (`--mirror-kink`, `kink_formation_belt4_pinned_mirror.json`): PVDF 9.16
  against 9.31 (achiral, equal within the 0.15 resolution), CFE 70/20/10 12.41 against 3.13, CDFE
  5.20 against 0.40. So for a chiral sequence the two kink handednesses differ by 5-9 kcal/mol in
  the pinned lattice, the chain takes the cheaper one, and the table above happens to carry it for
  these rows; every chiral row's number is an upper bound on the cheaper handedness until both
  have run. In the relaxed kink crystal the difference is small: the mirrored CFE pair costs 1.70
  against 1.24 in its own basin, so the 5-9 kcal/mol split is the pinned lattice's (the mirrored
  kink's lateral offset collides where the lattice cannot give), not the kink's own. The mirrored CDFE pair in its own basin is -3.20: a kinked CDFE
  crystal 3.2 kcal/mol per kink pair *below* the all-trans one, so the ten-monomer CDFE all-trans
  lattice is metastable against kinking (as the CTFE one was, -0.4), the far end of the
  chlorine's loosening and one more reason the CDFE rows are read in the fourteen-monomer cell,
  not this one. In the fourteen-monomer 64/29/7 cell the two kink handednesses are equal within
  resolution in their own basin (1.79 against 1.82), so its 1.8 is the number.
* **What the proxy leaves out.** The migration barrier (the hop that moves the kink, 0.357 eV in
  perfect beta from an earlier polyfind computation) is assumed equal across rows; the
  cooperativity `Tf` is not computed and the anchoring on the measured film fixes it at the
  terpolymer's; the kinks of the two chains are in register; the chlorine and nitrile rows are
  extrapolated chemistry; and a kink pair in the model's trans run is a crystal-interior event
  where the film's switching also happens at interfaces. The direction of every statement above
  is the model's; the factor of 500-1600 between copolymer and terpolymer is the size of the effect,
  not a measurement. The frequency map's channel also carries the margin's effect on the
  relaxation (a biased two-level channel relaxes faster than a degenerate one), so a row's channel
  frequency at the operating point is not the proxy column's; the tan delta columns are the channel's.

## 7. Confidence, and what decides each row

* **Validated against DFT:** the VDF and TrFE crystal terms (cells, polarization, stiffness,
  polar Berry slope, backbone angles; 5.22). **Not validated:** every chlorine parameter; the
  TrFE unit's own margins (the 70/30 model's conformational ordering is a prediction awaiting
  Dipole's E2 cells, 5.22); the polar-fraction mapping of section 3 (calibrated to one published
  film).
* **The helix handedness, the lesson of 2026-11-08:** an isotactic sequence has four TG
  conformers (two handednesses, two phases) and the first pass ran one; for the chlorine and
  nitrile rows the other handedness is 0.9-1.5 kcal/mol per monomer lower and turned every defect
  row's margin from near zero to +0.57-0.97. The TG side of a margin is now the lowest conformer
  that ran, which for a real atactic chain, free to choose its handedness locally at each
  defect, is still an upper bound on the helix's energy (and the margin a lower bound) by the
  cost of the handedness changes. The protocols' resolution is 0.15 kcal/mol per kink pair
  (section 8's translation check) and about 0.1-0.2 kcal/mol per monomer between symmetry-related
  packings of the chlorine rows, so margins are not read below 0.1.
* **The model cell:** a regular, periodic, isotactic ten- or fourteen-monomer sequence packed as
  two chains at zero kelvin. Real copolymers are statistical and atactic, their crystals are lamellae a few
  nanometres thick with the switching at and near the interfaces, and the measured strain is
  the film's, not the crystal's (finding 1 of the literature pass). The screen's margins are
  the crystal's; the film's switch costs less at the defects and the interfaces, which is why
  the measured thresholds sit below the crystal's `E*`.
* **Known biases of the potential:** the lateral coupling is 13% soft (C_12, 5.22); the
  pendant-bend constants set the tilt relaxation to within about 20% (5.22); the TGTG' cells of
  the copolymers carry the pendant relaxation and the TTTT cells almost none, so the margins
  here are 0.2-0.4 kcal/mol per monomer lower than a rigid-pendant model would give (5.22).
* **What decides each row:** Dipole's E2 cells (the relaxed TGTG' and TTTT cells of the 70/30
  model first, then one CFE-bearing sequence) decide the margins; the bench's Vogel-Fulcher
  fit and bias-suppression ratio (E4) decide the mapping of section 3; the microreactor's D7
  and D8 (a single-addition unit; strictly periodic spacing) decide whether placement can be
  controlled at all.
