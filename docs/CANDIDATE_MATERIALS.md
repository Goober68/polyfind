# Candidate chains and defect placements for the loss / frequency / strain optimum

*First pass 2026-11-07, corrected 2026-11-08 for the handedness of the TG helix (section 2: the
first pass ran one of an isotactic sequence's four TG conformers, and for every chlorine and
nitrile row the other handedness is 0.9-1.5 kcal/mol per monomer lower, which turned the defect
rows' margins from near zero to +0.57-0.97 and the reading of the defects from margin tuners into
TG-favouring pins), and again 2026-11-10 for two properties of the packing search that the
mixed-conformation cell exposed (sections 2.1 and 8: the search holds a chiral chain at the
symmetry of its state pattern, so every chiral row's margin moves when both of its cells are
relaxed with every torsion free, 64/29/7 from +0.11 to +0.32; and it holds a kinked chain periodic
only by a penalty), on what the model can say today: the crystal terms of
`pvdf-dft-valence-belt4` (docs/CRYSTAL_STAGE.md 5.22), the pendant-group relaxation of the
same section, the loss model of docs/ELECTROSTRICTION_RESEARCH.md 11-12 and the 2025-2026
literature pass. Code: `examples/candidate_screen.py` (the crystal screen),
`examples/candidate_loss_map.py` (the mapping onto the loss model); results
`docs/predictions/candidate_screen_belt4.json`, `candidate_loss_map_belt4.json`. Every number here
is a model number at zero kelvin on a regular ten-monomer periodic sequence packed as two chains;
section 7 says what that leaves out and what decides each row, section 6 reads Dipole's slate of
other chemistries against the same screen, and section 8 adds the frequency: what a kink pair is
in the lattice (a sheet, never a lone defect) and the barrier of the step that moves it. **Sections
3-6 still carry the search's margins and kink energies; where section 2.1 has a polished margin
for a row, that is the number, and the readings of 3-6 that rest on a margin within 0.3 of zero
are suspended until the polished table is complete (section 5 says which).***

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
| vdf9-trfe4-an1-14 | chlorine (extrapolated) | -0.154 (-0.494) | -0.393 | +10.8%, -15.4% | 0.145 | 0.524 | 0 | 0 |
| vdf9-trfe4-fan1-14 | chlorine (extrapolated) | -0.091 (-0.091) | +0.000 | +10.6%, -12.2% | 0.141 | 0.537 | 0 | 0 |
| vdf9-trfe4-ctfe1-14 | chlorine (extrapolated) | +0.135 (-0.190) | -0.416 | +9.5%, -12.1% | 0.138 | 0.508 | 115 | 390 |
| vdf13-trfe6-cfe1-20 | chlorine (extrapolated) | +0.055 (-0.253) | -0.426 | +10.1%, -13.8% | 0.156 | 0.552 | 43 | 146 |

### 2.1 The margins with both cells relaxed alike (2026-11-10)

The table above takes each conformation's cell from the packing search, which refines a
recognised state pattern (all-trans, TGTG') under the line group of the pattern: one torsion value
per symmetry class. For PVDF that is the molecule's symmetry, and re-relaxing its cells with every
torsion free changes nothing. For a chiral (isotactic copolymer) sequence it is a constraint the
molecule does not have. `examples/polish_margins.py` re-relaxes both cells of a row in the mixed
cell (every torsion and backbone angle of both chains free, the cell free, each chain periodic
exactly, started half a degree off the search's torsions), then the pendant tilts as before:

| sequence | the table's margin | both cells polished | all-trans cell lower by | TG cell lower by | TG conformer | torsions off the pattern, TT / TG | axial switch strain |
|---|---|---|---|---|---|---|---|
| 64/29/7, one CFE per fourteen | +0.110 | **+0.317** | 0.206 | 0.533 | TG'TG | 10 / 21 deg | +12.1% |
| 71/29, the matrix without the pin | -0.051 | **+0.056** | 0.081 | 0.375 | TGTG' | 6 / 19 deg | +11.7% |
| 64/36 matrix | -0.058 | **+0.083** | 0.061 | 0.407 | TGTG' | 6 / 21 deg | +11.8% |
| 65/30/5, one CFE per twenty | +0.055 | **+0.258** | 0.175 | 0.517 | TGTG' | 10 / 21 deg | +12.1% |
| 64/29/7 with CDFE in place of CFE | +0.402 | **+0.685** | 0.150 | 0.637 | TGTG' | 14 / 25 deg | +12.1% |
| one CFE and one CDFE per fourteen | +0.641 | **+0.807** | 0.251 | 0.581 | TG'TG | 19 / 20 deg | +12.2% |
| 64/29/7 with CTFE in place of CFE | +0.135 | **+0.300** | 0.186 | 0.664 | TG'TG | 14 / 24 deg | +11.8% |
| 64/29/7 with AN in place of CFE | -0.154 | **+0.197** | 0.087 | 0.614 | TG'TG | 5 / 24 deg | +12.5% |
| 70/30, ten monomers (the E2 cell) | +0.145 | **+0.213** | 0.078 | 0.305 | TGTG' | 6 / 19 deg | +11.6% |
| 80/20, ten monomers | +0.177 | **+0.229** | 0.097 | 0.241 | G'TGT | 6 / 14 deg | +11.3% |
| 60/40, ten monomers | -0.329 | **-0.199** | 0.049 | 0.414 | TGTG' | 6 / 20 deg | +11.0% |
| 70/20/10 CFE, ten monomers | +0.966 | **+0.868** | 0.248 | 0.154 | TG'TG | 11 / 13 deg | +11.6% |
| 70/20/10 CDFE, ten monomers | +0.146 | **+0.392** | 0.220 | 0.719 | TGTG' | 16 / 23 deg | +12.8% |

kcal/mol per monomer, all-trans minus TG, pendant tilts included (`polished_margins_belt4.json`).
Every long cell is in; of the ten-monomer rows the five above, and the rest are not polished (the
FAN cells cannot be: their tilts are not treated).

* **Both cells fall, the helix by more.** The all-trans chain of a VDF/TrFE sequence is not
  planar: its torsions deflect by 2-6 degrees in a pattern that follows the TrFE units (ten at the
  CFE), which is the deflected-trans chain known for the TrFE-rich copolymers, and it gains
  0.05-0.25 kcal/mol per monomer. The TG chain gains 0.15-0.72: held to one trans and one gauche
  value it had to strain its backbone angles to stay commensurate (0.73 kcal/mol per monomer of
  bend in the 64/29/7 cell), and free to vary its torsions along the sequence (trans bonds
  160-180, gauche 47-76 degrees) it closes at half that bend. Nothing in either polished cell is at
  a bound, and the closest contacts are the all-trans cell's own (F-H 2.30 against 2.36 A).
* **Every margin but one moves toward the helix, by 0.05-0.35** (the ten-monomer 70/20/10 CFE cell,
  whose helix the handedness correction had already found, moves the other way by 0.10). The matrix stays degenerate within the
  model's resolution but on the other side of zero (+0.06 and +0.08 for 29 and 36 mol% TrFE), and
  the 64/29/7 cell is at +0.32, not +0.11: the CFE pin is +3.7 kcal/mol per unit at one per
  fourteen (it was +2.3), and the cell no longer sits on the steep part of the switching curve at
  55 V/um (E* near 300 MV/m on the all-trans dipole). The twenty-monomer cell with one CFE is at
  +0.26 (it was +0.055), the same pin per unit (3.5-4.0 against a 29-36 mol% matrix), so the 1/N
  law of the earlier pass holds with the new numbers: the margin of a 29-36 mol% matrix with one
  CFE per N monomers is about +0.06 + 3.7/N, and the steep part of the curve (+0.07, section 3) is
  the matrix with almost no pin.
* **Every unit tried at one per fourteen now tips the matrix toward the helix.** Against the 71/29
  matrix (+0.056): CFE +0.317 (3.7 kcal/mol per unit), CTFE +0.300 (3.4), AN +0.197 (2.0), CDFE
  +0.685 (8.8), CFE with CDFE +0.807 (the two pins add: 3.7 + 8.8 = 12.5 against 10.5 found). The
  nitrile's sign has changed: on the search's cells it left the all-trans lattice lowest (-0.154)
  and was read as saturating the matrix; polished, its helix gains 0.61 kcal/mol per monomer and its
  all-trans cell 0.09, and it is a weak pin of the chlorines' kind. Section 6's reading of AN as the
  one unit that favours the all-trans lattice at the film's dilution is withdrawn.
* **What this does to the pass.** The search's cells are constrained minima of the same potential
  and the polished cells are lower, so the polished margins supersede the table's wherever both
  exist. The readings that do not move: TrFE content sets the matrix's margin and brings it to
  zero near 30 mol%; chlorine and nitrile units pin the helix; the pin's strength is a property of
  the unit in its matrix. The reading that does: where the film's composition sits on the
  switching curve. The search's +0.11 put the 64/29/7 crystal at a polar fraction of 0.24 at the
  bias, which is about what the measured film does; the polished +0.32 puts the perfect crystal
  near zero there. Either the chlorine's parameters, which no crystal anchor has tested, overstate
  its pin, or the film switches where the crystal's margin does not apply (interfaces, the
  statistical sequence, the segments between distant pins). The first is decided by Dipole's E2
  cell with a CFE in it; this document cannot decide it, and section 5's recommendation of one
  CFE per fourteen to sixteen, which was read off the search's margins, is suspended with it.
* **The pre-registered E2 numbers** (`docs/E1_E2_PREDICTIONS.md`) were the search's cells. The
  polished 70/30 cells are registered there as an addendum before any DFT cell of that sequence
  exists, with their coordinates (`docs/predictions/polished_vdf7-trfe3_*.xyz`).

### 2.2 Reading the search's table (as written 2026-11-08)

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
tips the degenerate matrix toward the helix and two of them tip it off the curve. CTFE in place of
the CFE: TTTT -5.952, TGTG' -6.086, margin +0.135, a pin of CFE's strength at this dilution (+2.6
per unit against +2.3; it costs both lattices 21-24 kcal/mol per unit, the bulkiest unit) with 12%
less polarization (0.138 against 0.156 C/m^2) and 10% less switch strain (9.5% against 10.5%).

**The pin density is the dial (`candidate_screen_belt4_20mer.json`).** A twenty-monomer cell with
one CFE per twenty in a 30 mol% matrix (`VVRVVRVVRVCVVRVVRVVR`, 65/30/5): TTTT -7.175, TG'TG -7.230,
margin +0.055 (the two handednesses within 0.005). With the pinless matrix at -0.05, one per
twenty at +0.055 and one per fourteen at +0.11, the tip is linear in the pin density: about +2.3
kcal/mol per CFE unit spread over the monomers between pins, so the margin of a 30 mol% TrFE
matrix with one CFE per N monomers is -0.05 + 2.3/N kcal/mol per monomer in this model, and the
composition the loss model wants on the statics alone (margin near +0.07, section 3) is one CFE per
sixteen to twenty, 5-6 mol%; with the dynamics (section 8: the sparser pin leaves a stiffer lattice
and a slower channel) the optimum moves to one per fourteen to sixteen. So the pin diluted to one per fourteen in a 29 mol% matrix
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
| vdf9-trfe4-an1-14 | -0.154 | 0 | 0.87 | 0.076% | 0.055 / 0.107 | 5.9 / 21.0 | 0.14 / 0.08 | 26 / 50 |
| vdf9-trfe4-fan1-14 | -0.091 | 0 | 0.87 | 0.074% | 0.054 / 0.105 | 5.9 / 20.9 | 0.14 / 0.08 | 25 / 48 |
| vdf9-trfe4-ctfe1-14 | +0.135 | 115 | 0.12 | 0.060% | 0.052 / 0.101 | 4.8 / 17.4 | 0.11 / 0.06 | 20 / 39 |
| vdf13-trfe6-cfe1-20 | +0.055 | 43 | 0.60 | 0.157% | 0.086 / 0.167 | 20.6 / 71.1 | 0.17 / 0.10 | 53 / 102 |

**Reading the mapping.** The operating point is 55 V/um bias, 1.8 V/um AC; the polar fraction
at bias is where the crystal margin puts the channel on its own switching curve, and the strain
per field is largest where that fraction is near one half: with `dmu` about 0.57 e.A the
margin that does it at 55 V/um is about +0.07 kcal/mol per monomer. After the handedness
correction no ten-monomer row sits there: the no-defect 70/30 and 80/20 and the CDFE row sit low
on the curve (polar fraction 0.05-0.11) and give 0.03-0.07% AC strain, the TT-ground 60/40 is
saturated (0.88) and gives the same, and every pinned defect row (CFE, CTFE, AN, FAN, margins
+0.57 to +0.97), PVDF, the alternating 50/50 and the VDF/CFE copolymers do not switch at all.
The rows that do sit on the curve are the long cells of section 2: the 64/29/7 cell (margin +0.11,
polar fraction 0.24 at 55 V/um) carries 0.125% AC strain at 34 kHz (tan delta 0.074, efficiency
0.16), the 65/30/5 twenty-monomer cell (margin +0.055, polar fraction 0.60) the largest strain of
the table (0.157%, tan delta 0.086, efficiency 0.17), and the 71/29 matrix control (margin -0.05)
is saturated (0.88) with half that strain, which is the pin doing its job: it holds the matrix off
saturation at the bias, and the pin density sets where on the curve the cell sits. Section 8
carries these rows with their own kink energies. The loss per
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
fourteen-monomer 29 mol% cell the same CFE pins at +2.3, CTFE at +2.6, CDFE at +6.3, AN at -1.4 (it
favours the all-trans lattice there) and FAN at about +4 with rigid pendants (sections 2 and 6), so
a pin strength is a property of the unit in its matrix and the fourteen-monomer cells are the ones
to design on. CTFE's two costs are both large (it is the bulkiest unit) and its pin is the
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

*Suspended in part, 2026-11-10. Items 1-4 below were read off the search's margins and the
search's kink energies. With both cells of a row relaxed alike (2.1) the matrix sits at +0.05 to
+0.08 instead of -0.05, the 64/29/7 cell at +0.32 instead of +0.11 (a pin of +3.7 kcal/mol per CFE,
not +2.3) and the CDFE cell at +0.69 instead of +0.40; with the kink pair measured like for like
(8.2) it costs 6.9, 3.8 and 5.1 kcal/mol in the 71/29, 64/29/7 and PVDF lattices instead of 5.5,
1.8 and 4.3. What survives: item 1 (the matrix is where the near-degeneracy comes from), the
direction of item 2 (the chlorine pins the helix and cheapens the kink sheet, by 3.1 kcal/mol) and
item 5. What does not: the density. On the polished numbers one CFE per fourteen holds the perfect
crystal off the switching curve at 55 V/um altogether (polar fraction 0.00 in the mapping), and the
margin law +0.06 + 3.7/N puts the steep part at pins far sparser than any terpolymer that is known
to work, while the measured 63/30/7 film does switch there. The model and the film now disagree on
that one number, and the two explanations (the chlorine's untested parameters overstate its pin;
the film switches where the perfect crystal's margin does not apply) are decided by a DFT cell with
a CFE in it, not here. No pin density is recommended until that cell or the remaining polished
rows say which.*

For a sequence-controlled material aimed at the bias-hold operating point, on this model and
with the confidence stated in section 7 (and after the handedness correction, which withdrew the
first pass's ranking):

1. **The matrix: VDF/TrFE at 29-36 mol% TrFE.** Degenerate within the model's resolution in the
   fourteen-monomer cells (71/29 -0.05, 64/36 -0.06; the ten-monomer 70/30 and 80/20 cells +0.15
   and +0.18), validated chemistry: this is
   where the near-degeneracy the relaxor needs comes from, as it does in every giant-strain
   terpolymer made so far. Without a pin it is a ferroelectric copolymer with hysteresis, which
   this equilibrium picture does not contain (the lattice model's coupling does).
2. **The pin: CFE at one per fourteen to twenty (5-7 mol%).** A TG pin (+7.9 kcal/mol per unit in
   the 20 mol% cell, +2.3 in the 29 mol% one) that holds the matrix off saturation at the bias
   (polar fraction 0.24 at one per fourteen, 0.60 at one per twenty, against the pinless matrix's
   0.88) and loosens the lattice for the kinks that switching needs (kink pair 5.5 kcal/mol in the
   71/29 matrix, 1.8 with the pin: section 8), so the cell both switches at 55 V/um and follows a
   34-66 kHz drive. The pin density is the dial: the margin is -0.05 + 2.3/N for one CFE per N
   monomers, and the loss model's steep part (margin near +0.07) sits at one per sixteen to twenty
   on the statics; but the sparser pin leaves a stiffer lattice (kink pair 2.3 at one per twenty
   against 1.8 at one per fourteen) and a slower channel (157 against 364 kHz), so its loss at the
   drive is 0.20 against 0.07 for 0.154 against 0.125% strain (section 8). Under the objective (loss
   first, stroke as a floor) the controlled sequence this pass recommends is a 30 mol% TrFE matrix
   with one CFE every fourteen to sixteen monomers, strictly periodic: tan delta 0.07-0.1 for
   0.12-0.14% AC strain at the operating point in the mapping. Its placement relative to TrFE is not a lever for the
   margin (0.005 per monomer); the pin spacing is (section 4), and a controlled sequence fixes it.
   Chlorine chemistry: extrapolated until a DFT cell tests it.
3. **Alternative pins.** In the ten-monomer 20 mol% cell FAN (+7.4) and CTFE (+6.1) pin as CFE
   does and AN at half the strength (+3.9); at the film's dilution in the 29 mol% matrix AN favours
   the all-trans lattice instead (-1.4 per unit, the cell saturated at the bias) and FAN pins at
   about +4 (rigid pendants) and overshoots like CDFE, so neither nitrile is the relaxor's pin in
   a TrFE matrix (section 6). CTFE at one per fourteen pins as CFE does (+0.135 against +0.11) with
   12% less polarization and 10% less switch strain, and its kink pair is -0.6: the all-trans CTFE
   lattice is metastable against kinking at this dilution as it was at one per ten, so its polar
   phase is a kinked one and the rate proxy has no meaning for it (section 8). It ranks below CFE on
   stroke and its dynamics are undetermined. CDFE does not pin in the ten-monomer 20 mol% cell (-0.3) but tips the fourteen-
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
| 64/29/7 with AN in place of CFE, fourteen monomers | -0.15 (-0.49): TTTT the ground state, the matrix alone is -0.05 | -0.39 | +10.8%, -15.4% | 0.145 | 0.516 | 0 |
| 64/29/7 with FAN in place of CFE, fourteen monomers | -0.09 rigid (tilts not treated for the FAN cell); about +0.25 with the other units' tilts added by hand | - | +10.6%, -12.2% | 0.141 | 0.516 | 0 (rigid) |
| 64/29/7 with CFE, fourteen monomers, for comparison | +0.11 (-0.19) | -0.42 | +10.5%, -14.2% | 0.156 | 0.552 | 280 MV/m |

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

At the film's dilution, one per fourteen in the 29 mol% TrFE matrix, the nitrile does the opposite
of what a relaxor pin does: AN takes the degenerate matrix (-0.05) to -0.15, the all-trans lattice
lowest (per unit: TTTT 0.0, the helix +1.4), so it saturates the matrix at the bias (polar fraction
0.88, the pinless matrix's) instead of holding it on the steep part as CFE does (0.24). FAN at the
same dilution is -0.09 with rigid pendants, which the other units' tilt release would turn to about
+0.25 (the FAN cell's tilts cannot be treated: its carbon carries a fluorine and a nitrile
fragment), a pin of about CDFE's kind with the same overshoot; its chemistry is the least tested
of the table. So in this model the nitriles are not the relaxor's pin in a TrFE matrix. AN does loosen the
lattice it leaves ferroelectric (kink pair 2.0 against the matrix's 5.5, channel 285 kHz on the
fourteen-monomer anchor): a fast saturated ferroelectric, which under the bias-hold operating point
gives the pinless matrix's strain (0.076%) and loss (tan delta 0.055), not the relaxor's.

**What this seat cannot compute.** CNEPO, SCN, SF5, nitro, isocyanide, the cyano-alkenes and
cyclopropyl need sulfur, nitro, epoxide and sp2 backbone parameters polyfind does not have; for
those Dipole's single-chain and GFN1 numbers are the only ones, and the barrier column above is
the one that matters for this programme's objective.

## 8. Frequency and loss: what a kink pair is in the lattice, and what moves

*Rewritten 2026-11-10. The first form of this section (kept below as 8.7) took the formation
energy of a kink pair in the kinked two-chain cell as the rate factor between chemistries, and
bracketed a "dilute" kink between that cell and a pinned lattice. `polyfind.mixedcell` (a cell
whose chains differ, every chain periodic exactly; `examples/dilute_kink.py`,
`kink_migration.py`, `kink_site_scan.py`, `polish_margins.py`, tables by `kink_summary.py` in
`docs/predictions/kink_summary_belt4.md`) has now computed the dilute kink itself, the kinked
crystal against a reference relaxed the same way, and the barrier of the step that moves a kink.
Three of the earlier section's statements do not survive, and they are named where they fall.*

### 8.1 A lone kink pair is not a defect the crystal can hold

A kink pair T G+ T G- T makes its chain about 1.0 A shorter (0.68 A once it has relaxed inside a
kinked crystal). In an all-trans crystal the chain cannot give that length up, because every
chain is held in axial registry by its six neighbours. Sliding one chain rigidly along its axis in
the polished all-trans four-chain cell costs, per monomer of the sliding chain:

| lattice | barrier over one monomer period | at 1 A | curvature, kcal/mol/A^2 |
|---|---|---|---|
| PVDF | 6.80 | 5.70 | 14.2 |
| 71/29 (fourteen monomers) | 5.84 | 4.91 | 12.8 |
| 64/29/7 (fourteen monomers) | 4.56 | 3.77 | 10.0 |

Most of it is the two neighbours along the short axis (the stack the dipoles point along, 4.75 A
apart, CF2 meshed between CH2): sliding one whole sublattice, which keeps each chain in register
with those two, costs 3.58 in PVDF, and that number has an independent check. Relaxed with the
slip held and the cell free it is 0.51 kcal/mol per monomer of the cell (the cell opens b from
8.48 to 9.14 A), against 0.56 from the DFT seat's GFN1 registry gate for the same slip
(`sarco/materials/computational/registry_gate_v1.md`). The chain-axis shears C44 and C55, which
are this registry's elastic constants, are in no anchor the potential was fitted to
(`docs/BENCHMARK.md`: structurally absent from the fit path), so the single-chain numbers are
predictions; the request for the DFT values is filed (section 8.6).

**The same slide with the lattice responding (added 2026-11-10; `examples/chain_slide.py`,
`docs/predictions/chain_slide_belt4.json`).** The table above is a rigid lattice. One PVDF chain
held at the half-monomer slip in an `nx x ny` tiling of the polished cell, every other chain
relaxed at fixed cell, kcal/mol per monomer of the slid chain:

| chains in the cell | rigid lattice | neighbours turn about their axes | turn and move aside | the same with the lateral cell free |
|---|---|---|---|---|
| 4 | 6.80 | 6.80 | 6.07 | 2.86 |
| 16 | 6.80 | 6.03 | 3.70 | 3.09 |
| 36 | 6.80 | 6.02 | 3.36 | 3.12 |
| 64 | 6.80 | 6.02 | 3.26 | 3.13 |
| dilute limit | 6.80 | 6.02 | **3.15** | **3.15** |

* **The single-chain barrier is 3.15 kcal/mol per monomer, 46% of the rigid number.** At fixed
  cell the supercell's own strain bounds the dilute limit from above; with the lateral cell free
  (one slid chain per supercell at zero lateral stress) it is bounded from below. The two series
  are monotonic and extrapolate, linearly in one over the number of chains from 36 and 64, to 3.14
  and 3.15. Flexing chains add nothing (two-monomer repeat, sixteen chains: the chains stay planar
  and the number is the same to four figures) and the pendant tilts 0.01.
* **What gives way is the lattice's spacing, not the chains' orientation.** Turning alone is worth
  11%. In the 64-chain cell the two stack neighbours (along the short axis, 4.75 A away) move
  aside by 0.20 A without turning, the four diagonal neighbours (4.86 A) by 0.17 A and turn 5
  degrees, the next shell by under 0.09 A; the slid chain neither turns nor moves. A four-chain
  cell has no room for this (6.07), which matters for every other four-chain number of this
  section (below). The coupling of slip to spacing is what E5 (b) asks the DFT seat for, the cell
  opening by 7.7% under a sublattice slip, so that request tests this number too.
* **It is a saddle, and needs no path.** The all-trans chain has a mirror plane normal to its
  axis, so at the half-slip the axial force on every chain is zero whatever the lattice does
  laterally. The relaxed state is a minimum in every other coordinate (reached from perturbed
  starts as well as the symmetric one; in the four-chain cell the symmetric start stops on a
  stationary point 0.56 higher), and releasing the other chains' axial shifts leaves it where it
  is in every cell of sixteen chains or more (energy change under 1e-7; in the four-chain cell,
  where every chain is a neighbour of the slid one, it falls by 0.75, one more reason not to use
  that cell). A stationary point with one descending direction is the transition state of the
  slide. The relaxed energy rises to it smoothly, the same there and back
  (sixteen chains; slip in monomer periods, rigid then relaxed: 1/8 0.74, 0.66; 1/4 2.86, 2.04;
  3/8 5.52, 3.24; 1/2 6.80, 3.70).
* **The bottom of the well does not move.** The lattice's give is second order in the slip: the
  curvature at the minimum is 14.2 kcal/mol/A^2 per monomer rigid and 14.2 relaxed (13.5 by 0.2
  A), and at 1 A the relaxed energy is 3.0 (36 chains) against 5.7. For small displacements the
  registry is as stiff as the table says; at the top it is half as high.
* **A different softness, which is not a barrier.** If the neighbours' axial shifts are free too,
  the lattice held only far away, they follow a displaced chain (the first shell by 62% of its
  displacement in sixteen chains) and the curvature the chain sees falls to 5.8 (sixteen chains)
  and 4.9 (thirty-six): anti-plane shear of the lattice, whose compliance grows as the logarithm
  of the distance to whatever holds it. That is the elastic constant C44/C55, it vanishes at the
  half-slip (no axial forces there), and it leaves the saddle and the 3.15 where they are.
* **A long chain does not slide this way.** These are energies per monomer of a uniform slide: the
  amplitude of the registry. A stem of N monomers sliding uniformly pays N times as much, and
  beyond a few tens of monomers a localized defect is cheaper: its cost goes as the square root of
  this amplitude times the chain's stiffness (the Frenkel-Kontorova estimate used below), not as
  the amplitude times the length. A first sign is already in these cells: with a ten-monomer
  repeat in the four-chain cell the relaxation finds a state 0.04 lower in which the chains are no
  longer uniform along their length (torsions up to 13 degrees off trans) and which is not
  stationary along the slide.
* **The copolymer lattices give the same factor.** With one chain in four slid, the other chains'
  axial shifts held and the lateral cell free (the affordable protocol for fourteen-monomer
  chains, and not a saddle search; for PVDF it gives 2.86 against the dilute limit's 3.15), the
  highest point of the relaxed slide over one monomer period is 2.37
  kcal/mol per monomer in 71/29 (rigid 5.84) and 2.04 in 64/29/7 (rigid 4.56; the relaxed top is
  at 0.6 of a period, the rigid one at 0.5): 41% and 45% of the rigid lattice, PVDF 42%. The order
  of the three lattices is unchanged; the chlorine's loosening of the registry, 0.78 of the 71/29
  value in a rigid lattice, is 0.86 relaxed. After one whole period the slid chain sits in another
  sequence register with its neighbours, 0.95 (71/29) and 1.19 (64/29/7) above where it started.

`examples/dilute_kink.py` puts one kinked chain among all-trans neighbours: the polished
all-trans two-chain cell tiled 2 x 1 (four chains), the cell held at the all-trans values (the
boundary condition of a kink far from any other), every chain periodic exactly, relaxed in three
stages that each bound the next from above.

| sequence | chains x monomers | neighbours fixed | neighbours free to shift and turn as rigid chains | neighbours fully flexible |
|---|---|---|---|---|
| PVDF | 4 x 10 | +25.25 | +23.07 | reverts to all-trans |
| PVDF | 4 x 14 | +24.17 | +22.69 | reverts to all-trans |
| PVDF | 4 x 20 | +23.41 | +22.39 | not run |
| PVDF, sixteen chains (added 2026-11-10) | 16 x 10 | +25.61 | +22.2, a shelf and not a minimum (last bullet) | not run |
| PVDF, chain stiffness matched to DFT (see below) | 4 x 10 | +19.77 | +18.23 | reverts to all-trans |
| 64/29/7 | 4 x 14 | +21.16 | +20.13 | reverts to all-trans |
| 71/29 | 4 x 14 | reverts to all-trans | - | - |

kcal/mol per kink pair against the all-trans supercell (which stage `host_check` confirms is the
minimum with every chain free: 0.0000 below the tiled cell in all three lattices; the sixteen-chain
cell was not given that check, which costs a day there). The first column's numbers are minima.
The second column's are shelves: minima held by 0.02 kcal/mol in the four-chain cells and by
nothing in the sixteen-chain one (last bullet below).

* **It costs 20-25 kcal/mol where the neighbours hold still** (18-20 with the chain's stiffness
  corrected, below), 30-42 kT at 300 K: a population of 1e-13 per site at most. The strained kink's gauche bonds are flattened to 85-99 degrees and its
  backbone angles opened to 122 degrees.
* **The cost is a strain field, not the size of the cell.** Against the repeat length the PVDF
  numbers fall as 21.6 + 36/N (neighbours fixed) and 21.7 + 13/N (rigid bodies): the limit of an
  isolated kink pair on an infinite chain is 21.7, of which 5.3 is the kink pair itself (its energy
  in a sheet, 8.2) and 16 the stretch of the chain on either side of it, held by the registry over
  a decay length of about seven monomers (the Frenkel-Kontorova length from the model's own chain
  stiffness and registry curvature; the same estimate gives 12-28 for a 0.68-1.04 A deficit).
  *Added 2026-11-10:* all of these are four-chain cells, in which the lattice cannot give way
  laterally. That give halves the top of the registry (above) and leaves its bottom alone, and the
  kink's strain field lives at the bottom: in a sixteen-chain cell the ten-monomer number with
  the neighbours fixed is 25.6 (25.3 in four chains), and with the neighbours shifting and turning
  the energy is 22.25 when the kink has the shape of the four-chain minimum (gauche bonds at 99
  degrees), against that minimum's 23.07. Under one kcal/mol at equal shape, not a factor of two.
* **A known bias in that number, and its size.** The mixed cell's chains are built from torsions
  and backbone angles with rigid bonds, and with the preset's angle restraint (233 kcal/mol/rad^2,
  fitted to the backbone angles) their axial stiffness is 674 GPa, where the DFT seat's beta cell
  has C33 = 303 GPa and the potential's full valence form 330. The strain part of the lone kink
  scales as the square root of that stiffness. Run again with the restraint lowered to 100, which
  gives the chain 327 GPa (`--angle-stiffness 100 --angle-window 14`; the softer restraint shifts
  the all-trans backbone angle from 115.7 to 117.1 degrees, so this is a sensitivity run and not a
  better preset), the ten-monomer numbers are 19.8 and 18.2 where the preset gives 25.2 and 23.1:
  the lone kink with a chain of the right stiffness is about 18 at ten monomers and 16-17 at
  infinite repeat. The sheet energies and the step of 8.2-8.3 involve no stretch of this kind.
* **It is not even metastable once the neighbours respond.** With every chain's torsions and
  angles free the kink reverts to all-trans without a barrier, in PVDF and in the terpolymer
  lattice; in the 71/29 lattice it reverts with the neighbours fixed. *Checked 2026-11-10.* The
  parenthesis that stood here said the relaxation's damped descent before its quasi-Newton step
  made a reversion the basin's fate and not an optimiser's jump. That descent is short and proves
  nothing, so the reversions were run again by a damped descent alone, with steps of 0.05-0.1 (a
  degree of torsion, a tenth of an angstrom), which cannot cross a barrier
  (`dilute_kink.py --descend-to`; a return to all-trans that the descent did not itself reach is
  never accepted). They stand: PVDF with flexible neighbours, four chains, falls in every round,
  23.1, 15.6, 14.1 (gauche bonds at 124 degrees), 0.2, 0.03; 71/29 with the neighbours fixed falls
  from its strained start through 10.3 to 0.01. The terpolymer's and the fourteen-monomer PVDF's
  were not rerun.
* **How well the kinked state is held, measured.** The energy against the kink's own torsion: the
  two gauche bonds held at +tau and -tau, the rest of the stage relaxed, every point from the
  stage's stored state (`examples/kink_torsion_scan.py`; ten monomers, kcal/mol above all-trans,
  tau in degrees):
  - neighbours fixed, four chains: 93: 25.33; 99: 25.49; 105: 25.94; 112: 26.58; 116: 26.90; then
    the rest of the chain jumps to another branch (120: 19.12; 130: 14.09);
  - neighbours shifting and turning, four chains: 99: 23.08; 101: 23.09; 103: 23.09; 105: 22.93;
    108: 22.36; 112: 21.00; 120: 17.66; 130: 13.41;
  - neighbours shifting and turning, sixteen chains: 100: 22.17; 112: 20.21; 120: 17.03; 130: 12.98.

  With the neighbours fixed the kink is a real minimum (25.25, gauche bonds at 95 and -91
  degrees) behind a rim of at least 1.6 kcal/mol (26.90 at 116 degrees against the scan's own
  minimum of 25.33), 2.6 kT at 300 K. With the neighbours shifting
  and turning in the four-chain cell the minimum (23.07 at 99 degrees) has a lip of 0.02 kcal/mol
  at 102-103 degrees and nothing after it: a shelf. In sixteen chains there is no lip: the
  quasi-Newton relaxation returns the kink to all-trans, a damped descent with steps too short to
  jump crawls across the shelf without stopping (25.61, 22.25, 22.13, 22.07, 22.01 at 200
  evaluations apart, gauche bonds 94, 99, 102, 104, 106 degrees, when its budget ended), and the
  scan continues it down through the eclipsed position. So the second column of the table is the
  height of a shelf the strained kink crosses on its way back, and the lone kink is a state with
  a lifetime only where nothing around it can move. (An interim note here, and two commit
  messages, read the slow descent as a minimum and the quasi-Newton reversion as an optimiser's
  jump. The scan shows the relaxation was right and the descent merely slow on a flat stretch.)
* **So the "dilute kink" the earlier section asked for does not exist**, and its pinned-lattice
  column was not an upper bound on one: a kink pair on every chain of a cell held at the all-trans
  a and b (9.3 in PVDF) is still a sheet, with no chain shorter than its neighbour, and costs less
  than half of what one kink pair alone does. That statement is withdrawn.
* **Kink pairs exist in register**: on every chain at the same height, a sheet normal to the chain
  axis, which shortens the whole crystal by 0.68 A and strains nothing; or at a chain end, where
  the short side can slide. The edge of a sheet that covers part of a crystal is a dislocation line
  with its Burgers vector along the chain, and the lone kink is its smallest loop; the energies of
  partial sheets (two to eight chains kinked among sixteen) are the nucleation calculation and are
  not run.

### 8.2 The sheet, measured against a reference relaxed the same way

The kinked crystal of the earlier section is that sheet, one per repeat. Two things were wrong
with its numbers, both found by re-relaxing in the mixed cell (every torsion and backbone angle of
both chains free, the cell free, each chain periodic exactly):

1. **The closure was soft.** The packing search holds a chain periodic by a penalty, and a kinked
   chain keeps 0.33-0.40 degrees of residual rotation per repeat and gains 0.8 kcal/mol per kink
   pair from it. PVDF's +4.47 (rigid pendants) is +5.30 with the closure exact.
2. **The all-trans reference of a chiral sequence was not relaxed.** The search refines a
   recognised state pattern (all-trans, TGTG') under the line group of the pattern, one torsion
   value per symmetry class. That is the molecule's own symmetry for PVDF, whose polished cell is
   unchanged to 1e-4. A copolymer sequence has no such symmetry: its all-trans chain was held
   planar, and polished it is deflected by up to ten degrees per bond and lies 0.06-0.21 kcal/mol
   per monomer lower (section 2.1). The kinked chain, which has no recognised pattern, was refined
   with every torsion free. The comparison was not like for like, by 2.9 kcal/mol per chain in the
   64/29/7 cell.

| sequence | the search's cells | both cells polished, rigid pendants | with pendant tilts | kinked cell a, b, gamma, c | all-trans cell |
|---|---|---|---|---|---|
| PVDF | +4.47 | +5.30 | +5.06 | 4.958, 8.774, 89.18, 25.201 | 4.754, 8.481, 90.00, 25.877 |
| 64/29/7, kink pair in a VDF-VDF pair | +0.94 | +3.95 | +3.79 | 5.402, 8.426, 92.67, 35.496 | 5.302, 8.245, 91.70, 36.183 |
| 71/29, kink pair in a VDF-VDF pair | +5.68 | +7.09 | +6.85 | 5.150, 8.549, 92.16, 35.529 | 4.881, 8.459, 91.78, 36.187 |

kcal/mol per kink pair against the all-trans crystal.

* **The chlorine's loosening of the lattice for kinks survives, on higher numbers.** The earlier
  table had the kink pair at 4.3 in PVDF, 5.5 in the 71/29 matrix and 1.8 in the 64/29/7 cell.
  Like for like it is 5.1, 6.9 and 3.8 (5.3, 7.1 and 3.95 with rigid pendants): every number is higher, the copolymer
  matrix by most because its all-trans reference gained most from the polish relative to its kinked
  cell, and one CFE per fourteen still cuts the matrix's kink pair by 3.1 kcal/mol (it was 3.7) at
  this site, a VDF-VDF pair with a TrFE after it; the matrix has cheaper sites in its VDF run and
  between the cheapest sites of the two lattices the cut is 1.2 (8.4).
  The kinked cells are 2-6% wider along a than the all-trans cells (PVDF 4.75 to 4.96 A, 71/29 4.88
  to 5.15, 64/29/7 5.30 to 5.40): the terpolymer's all-trans lattice is already most of the way
  there, which is the loosening.
* **Kink pairs attract along the chain.** One kink pair per fourteen monomers costs the
  terpolymer chain +3.8, while the chain with a kink pair in every two monomers (the TG helix) lies
  4.4 below all-trans (section 2.1: +0.317 per monomer). The low-energy objects of this lattice are
  therefore not kinks but runs of helix and runs of trans with a wall between them, about 2
  kcal/mol per wall per chain, which is the two-state unit of the loss model with its coupling
  (plan sections 11-12) seen from inside the crystal.

### 8.3 The step that moves a kink pair

`examples/kink_migration.py` moves ONE chain's kink pair one monomer ahead of its neighbours' in
the kinked crystal, where every chain has a kink pair and none is shorter than another. The move
is two gauche migrations (counter-rotations of two second-neighbour bonds, the torsional move
that leaves both stems parallel) through the wide kink, each segment a climbing-image band on the
closure surface, every chain flexible, the cell fixed at the kinked crystal's:

    A   G+ T G- T T     in register with the neighbours' (the sheet)
    I   G+ T T T G-     the leading gauche bond moved two bonds on
    B   T T G+ T G-     the trailing one followed: one monomer ahead

| lattice | who moves | I | saddle A-I | saddle I-B | B | barrier ahead / back |
|---|---|---|---|---|---|---|
| PVDF | one sublattice (half the chains, each with its two stack neighbours) | +7.01 | 7.54 | 8.91 | +1.97 | 8.9 / 6.9 |
| PVDF | one chain in four (all six neighbours stay) | +9.85 | 10.27 | 12.46 | +7.88 | 12.5 / 4.6 |
| PVDF | one chain in sixteen (2026-11-12) | +9.96 | 10.73 (the band's top; see below) | 11.64 | +5.91 | 11.6 / 5.7 |
| 64/29/7, from a VDF-VDF pair onto VDF-TrFE | one sublattice | +6.33 | 6.43 | 11.48 | +8.98 | 11.5 / 2.5 |
| 71/29, from a VDF-VDF pair onto VDF-TrFE | one sublattice | +6.47 | 7.05 | 12.02 | +11.25 | 12.0 / 0.8 |
| PVDF | every chain together (the sheet in register), per chain | +5.06 | 7.38 | not settled (the band's top 7.9) | 0.00 | at least 7.4 each way, 14.8 per two-chain cell |

kcal/mol per moving chain against the sheet (`kink_migration_belt4.json`; tables with the bands'
tops beside the saddles in `kink_summary_belt4.md`).

**The saddles are stationary points, not band tops (2026-11-12).** Each band's highest image was
settled by minimising the energy in the hyperplane normal to the band's tangent and moving the
hyperplane until the force along the tangent vanished (`kink_migration.py --polish`,
`polyfind.mixedcell.saddle_polish`): under 0.01 kcal/mol per scaled unit along the path and across
it, with the curvature along the path negative. This replaced the bands' tops because a band's
force criterion does not settle a top that is soft across the path: the four-chain second top
stood at 12.66 when its band stopped at a force of 0.16, crept to 12.49 over 600 more iterations
at forces of 0.08-0.19, and is 12.46 settled; the sixteen-chain one was still falling
(11.84) when its band met the 0.1 criterion and is 11.64. The tops of the other
bands that had converged were within 0.03 of their saddles. The terpolymer's
second top, quoted before as 11.4 and 11.7 from two runs that stopped at their iteration limit, is
11.48. The in-register band had stopped at forces of 0.8-1.0: its first top was
0.86 too high and its second did not settle from that band.
The sixteen-chain first top is the one number of the table that is a band's (10.73 at a force of
0.09): its settling was stopped after 105 minutes in its first hyperplane (a second copy with a
progress log went from 10.73 to 10.69 in its first 70 evaluations and was still descending), and
since it lies a kcal/mol under the second saddle it does not enter the barrier.

* **The torsional step itself costs 7.5-8.9 kcal/mol in the PVDF lattice**: with a 3e12 /s attempt
  frequency, 1e6 steps per second at 300 K. The two saddles are sequential, not concerted: at the
  second the trailing gauche bond is at 110 degrees and the bond it is moving to has not started.
* **A sheet is stiff against one chain leaving it, and the cost is now measured in a cell where
  it is one chain** (2026-11-12: sixteen chains, one in sixteen ahead, every chain flexible). The
  completed step lies 5.9 kcal/mol up behind a barrier of 11.6, with 5.7 to come back; in the
  four-chain cell the same three numbers are 7.9, 12.5 and 4.6. So the cramped cell overstated the
  end state by 2.0 and the barrier by 0.8, and the 6.4 guessed on 2026-11-10 from rigid-chain
  neighbours is 5.91. With its stack moving too the step costs 2.0 behind 8.9 (row one), so one
  chain alone still pays three times what a chain of a moving stack does.
  - *Checks.* The sixteen-chain sheet, perturbed and relaxed with every chain and the cell free,
    returns to the tiled two-chain cell (2e-09 kcal/mol per chain, torsions within
    0.001 degrees). B was reached by two routes, from the rigid-chain-neighbour state and from
    the sheet with the moving chain's torsions reset: 5.911 and 5.911, the two structures within
    0.2 degrees in every torsion. In I and in B the moving chain's gauche bonds are where the step
    puts them and the other fifteen chains keep their kink pairs in the sheet.
  - *It is the completed step that needs the room, and it needs it along the stack.* B is
    +7.88, +7.51, +6.47 and +5.91 in cells of 2a x b, 2a x 2b, 4a x b and 4a x 2b (a the 4.96 A
    stack axis): doubling a is worth 1.4-1.6, doubling b 0.4-0.6. The wide kink is +9.85, +9.98,
    +10.02 and +9.96 in the same four cells, and its saddle 10.27 in four chains and at most 10.73 in
    sixteen: it does not gain from the room. Sixteen chains is therefore not the dilute limit
    of B; a wider cell along a was not run (an energy of the sixteen-chain cell is 4 s after
    the speed work of this date, and of a thirty-two-chain one an estimated five to six times
    that).
  - *How the lattice gives.* In B the moving chain has a monomer's length out of axial register
    with its neighbours (a kink pair shortens its chain by about 1 A, and this one sits a
    monomer further on than theirs), and it shares that out: the chain slides 0.32 A back
    along its axis against the rest, its six neighbours shift by up to 0.19 A along theirs,
    move aside by under 0.1 A and adjust their torsions by up to 11 degrees; beyond the first
    shell nothing moves by more than 0.06 A or 3 degrees. This is the single-chain slide of 8.1
    again (6.07 in four chains, 3.70 in sixteen, 3.15 dilute), which is why the two depend on
    the cell alike. In the wide kink the answer is one diagonal neighbour's torsions
    (16 degrees) and little else.
  - *The stack follows cheaply.* A second chain one monomer ahead beside the first along a
    brings the two to +8.35, 4.17 each, against 11.8 for two apart: the second costs
    2.4 where the first cost 5.9. A diagonal neighbour instead (5.07 A away, the other
    sublattice) brings them to +9.05, so it costs 3.1: the stack axis is the cheap direction,
    as it is for the registry (8.1). A whole sublattice ahead is 1.97 per chain. So the sheet
    advances by one chain stepping out, 1e4 /s per chain at 300 K with a 3e12 /s attempt
    frequency (2e3 /s on the four-chain barrier), returning at 2e8 /s unless a neighbour
    follows, which costs it less than half as much in energy. The barrier of that second step,
    and so the rate at which a stepped chain is joined rather than recalled, is not computed.
* **In the copolymer lattices the first saddle is lower and the completed step is not.** Into the
  wide kink the barrier is 6.4 in 64/29/7 and 7.0 in 71/29 against PVDF's 7.5, so the torsional
  event is no slower beside TrFE and chlorine. But the step computed here takes the kink pair from
  a VDF-VDF pair onto a TrFE unit, and it lands 9.0 (64/29/7) and 11.3 (71/29) kcal/mol up, far
  more than the 1.3 the whole sheet pays for the same move (8.4): half the chains with their kink
  pairs on a different unit from their neighbours' is a costly arrangement in a lattice of
  sequence-ordered chains, and the model cell, in which every chain carries the same sequence in
  register, makes it as costly as it can be. With that cost the barrier of the whole step is 11.5
  and 12.0 kcal/mol, so in this protocol the terpolymer's step is a factor of 2.5 faster than
  the matrix's and no more.
* **Moving every chain together is not the way round that cost.** In PVDF the concerted step of
  both chains of the cell has a first saddle of 7.4 kcal/mol per chain (14.8 for the cell; the band,
  stopped at a force of 1.0, had 8.2), and a second that the band puts at 7.9 and that did not
  settle from it: at least 14.8 for the cell against 8.9 for one sublattice followed by the other. The
  sheet moves one chain (or one stack) at a time. In
  the copolymer cells the concerted band does not settle on a minimum-energy path at all (before
  the climbing image is switched on its top is 12.7 per chain in 64/29/7, twice the sequential
  first saddle, and nine images cannot resolve what follows), so those rows are not reported. What
  would separate the chemistries cleanly is the sequential step in a cell whose chains do not carry
  the same sequence in register, which is a larger cell than the two-chain one and is not run.

### 8.4 In a copolymer the sheet sees the sequence

`examples/kink_site_scan.py` walks the sheet along a sequence one monomer at a time from the
polished kinked crystal, relaxing the cell at each position. Sheet energy in kcal/mol per kink
pair against the all-trans crystal, rigid pendants, the kink pair named by the two units it sits
on, in order along the chain (`kink_site_scan_belt4.json`; profiles in `kink_summary_belt4.md`):

| sequence, kink pair | profile along one repeat |
|---|---|
| 64/29/7 `VVRVVRVCVVRVVR`, G+ T G- | VV 3.95, VR 5.30, RV 4.56, VV 4.06, VR 5.25, RV 5.88, **VC 12.11, CV 11.41**, VV 4.75, VR 5.45, RV 4.47, VV 3.96, VR 5.33, RV 4.51 |
| 64/29/7, the mirrored pair G- T G+ | VV 4.18, VR 5.45, RV 3.73, VV 4.05, VR 3.57, RV 3.42, **VC 8.07, CV 8.28**, VV 4.14, VR 4.46, RV 3.83, VV 4.23, VR 4.36, RV 3.81 |
| 71/29 `VVRVVRVVVVRVVR`, G+ T G- | VV 7.09, VR 6.50, RV 7.20, VV 7.05, VR 6.47, RV 6.87, VV 4.81, VV 5.14, VV 7.09, VR 6.47, RV 7.20, VV 7.09, VR 6.50, RV 7.21 |
| 71/29, the mirrored pair | VV 5.79, VR 8.27, RV 5.42, VV 5.81, VR 8.28, RV 5.08, VV 4.66, VV 4.97, VV 5.74, VR 8.40, RV 5.41, VV 5.79, VR 8.28, RV 5.41 |
| 71/29 with the TrFE in pairs `VVVVVRRVVVVVRR`, G+ T G- | VV 5.05, VV 5.11, VV 5.40, VV 7.83, VR 10.41, RR 7.37, RV 7.16, and the same again (the two halves of the repeat agree to 0.001) |

Is a site's energy the site's, or the packing basin its relaxation fell into? In the 71/29 cell,
where neighbouring VDF-VDF positions come out 2 kcal/mol apart in cells 0.15 A apart in a, four
sites re-relaxed from the low sites' lattice returned to their scan energies and cells to three
decimals: the lattice follows the site.

* **In the matrix a kink pair costs 5 where VDF follows it and 7 where TrFE does.** A VDF-VDF pair
  whose next unit along the chain is VDF costs 4.7-5.4 (PVDF: 5.3); the same pair with a TrFE next
  costs 7.0-7.8, a pair with a TrFE in it 6.5-7.2, and one sitting on two TrFE units or leading into
  them 7.4-10.4. The kink pair is directional, and the mirrored one reads the same sequence
  differently (5.1-5.8 on the pairs the first finds at 7, 8.3-8.4 where the first finds 6.5, and
  the same 4.7-5.0 inside the VDF run), so a chain's cheapest kink away from a VDF run depends on
  its handedness as well as its place.
* **The chlorinated unit lowers every site of the matrix and is itself a wall.** In the 64/29/7
  lattice the TrFE-flanked VDF-VDF pairs that cost 7.1 in the matrix cost 3.95-4.06, and no matrix
  site is above 5.9 for either handedness: the wider lattice (8.2) takes the TrFE's penalty away.
  That is the loosening, 3.1 kcal/mol at the same site and 1.2 between the cheapest sites of the
  two lattices (3.42 against 4.66). With a gauche bond on either side of the CFE itself the sheet
  costs 11.4-12.1 for one handedness and 8.1-8.3 for the other, 4 to 8 above the matrix sites:
  neither kink pair sits at the CFE or passes it. The pin of the statics, which favours the helix
  around itself by 3.7 kcal/mol per unit (2.1), is to a single kink pair a wall, and a sheet lives
  in the segment between two of them. (The kink pair on the odd bonds is not scanned.)
* **What this says to a sequence designer.** At 71/29 the regular sequence gives the sheet two
  cheap sites in fourteen (its one four-VDF run) and a floor of 6.5-7.2 everywhere else, a
  corrugation of 0.7; grouping the TrFE in pairs gives it three cheap sites in seven (5.05-5.40)
  and a wall of 10.4 at each pair. Summed over sites, exp(-E/kT) at 300 K is 1.0e-3 per repeat and
  handedness for the regular sequence and 1.1e-3 for the paired one (one handedness scanned),
  against 11e-3 for the 64/29/7 cell: VDF runs are where kink sheets are cheap, TrFE units next to
  them are what they climb, and a chlorine does for every site what a VDF run does for its own. Whether cheap sites
  or a low wall matters more is the question of which rate model holds (8.5).

### 8.5 What the frequency map becomes

The map of the earlier section scaled each chemistry's relaxation rate by `exp(-E_form / kT)`
with the kinked-cell formation energies of the search, anchored on the measured 364 kHz. Two
things have changed under it: those energies carried the biases of 8.2, and the step of 8.3 now
exists as a second candidate for what limits the rate. Neither rate model is established, so the
map is given on both, for the three lattices that have a like-for-like sheet energy and a band,
with the polished margins (`candidate_frequency_map.py --kinks kink_formation_belt4_polished.json`
and `--steps kink_migration_belt4.json`; the 64/29/7 cell is the anchor at 364 kHz in both):

| lattice | kink sheet, kcal/mol per kink pair | f_relax if the sheet's energy sets the rate | step barrier (one sublattice) | f_relax if the step sets the rate | margin (polished) | polar fraction at 55 V/um |
|---|---|---|---|---|---|---|
| 64/29/7 | 3.79 | 364 kHz (anchor) | 11.5 | 364 kHz (anchor) | +0.317 | 0.00 |
| 71/29 | 6.85 | 2 kHz at this site; 33 kHz summed over the sequence's sites (8.4) | 12.0 | 147 kHz | +0.056 | 0.61 |
| PVDF | 5.06 | 43 kHz | 8.9 | 27 MHz | +0.744 | 0.00 |

(`candidate_frequency_map_belt4_polished.md` and `_steps.md`. At the operating point the 71/29 row
reads 0.012% AC strain at tan delta 1.0 on the first rate model and 0.18% at 0.21 on the second;
the other two rows carry no strain on either. The step barriers are on the saddles settled on
2026-11-12 (8.3): the terpolymer's is 11.48 where the band's top was read as 11.37, which moves the
71/29 row from 121 to 147 kHz and its tan delta from 0.26 to 0.21, and PVDF from 21 to 27 MHz.)

* **The two rate models bracket the matrix between slow and following.** If the population of
  kink sheets limits the channel (the earlier section's assumption, on corrected energies) the 71/29
  matrix sits at 2 kHz on the energy of the one site the sheet was first put on and at 33 kHz on
  the sum over all the sites of its sequence and both kink handednesses (8.4: 9% of the 64/29/7
  cell's), at the bottom of a 34-66 kHz drive band rather than far below it. If the step limits it, the matrix is within a factor
  of two and a half of the terpolymer, because the chlorine hardly lowers the step. The bench's Vogel-Fulcher
  fit on a copolymer film beside the terpolymer's decides between them; the model's own next piece
  is the wall of a helix run (8.2), whose creation and motion are the two quantities these proxies
  stand in for.
* **Only the pinless matrix is on the switching curve, and whether it delivers depends on the
  rate model.** The matrix sits on the steep part (polar fraction 0.61): frozen and lossy on the
  first rate model, and on the second the best row this document has had (0.18% at tan delta
  0.21, loss per strain 1.2), which is a ferroelectric copolymer and in a real film switches with
  hysteresis this equilibrium mapping does not contain; the 64/29/7 cell
  is fast and pinned off the curve (+0.32); PVDF does not switch. The earlier section's row for the
  64/29/7 cell (tan delta 0.07 for 0.125% strain) came from the search's margin of +0.11 and is
  withdrawn with it (section 5).
* **Absolute rates.** With a 3e12 /s attempt frequency the step barriers give 1e6 /s (PVDF), 1e4 /s
  (64/29/7) and 5e3 /s (71/29) per chain at 300 K, so the terpolymer's computed step is two
  orders slower than the measured channel (364 kHz is 2e6 /s). Its saddle is no longer the
  uncertain part (settled at 11.48, 8.3); what remains is that the model cell overstates the cost
  of a chain stepping onto a TrFE unit ahead of its neighbours, and that the film's channel need
  not be this step. These are the steps of half the chains at once, in a two-chain cell; one PVDF
  chain alone, in sixteen, has a barrier of 11.6 (1e4 /s), and the copolymers' single-chain steps
  are not computed. The earlier
  polyfind number for a kink-row hop in perfect beta, 8.2 kcal/mol (`ask ... verification-v1.md`,
  E3), is PVDF's 7.5-8.9 here by a different route.

### 8.6 What decides these numbers

Everything in 8.1-8.4 rests on one property no anchor of the potential has tested, the energy of
sliding a chain along its axis, and on the energy of a kink pair in a sheet. The request for both
is with the DFT seat (`sarco/materials/ask-polyfind-dipole-registry-kink-v1.md`, 2026-11-10; **E5
and the E2 addendum dispatched by the PI the same day**, E6 not; they wait for the alpha and Gamma
runs): **E5**, single points on the accepted beta cell with one sublattice and then one chain of
the doubled cell slid along c (polyfind registers +1.79 kcal/mol per monomer of the cell for the
half-monomer sublattice slip, +0.51 relaxed with the cell free, +1.70 for one chain in four); the
**E2 addendum**, the polished 70/30 cells (margin +0.213) and the polished 64/29/7 pair for
whenever the chlorine cell is declared (+0.317); and **E6**, the PVDF kink sheet (+5.30 per kink
pair). If the DFT registry is a factor of two softer than the model's, the lone kink falls to
about 13-16 kcal/mol and nothing above changes; if it is a factor of ten softer, kinks are free
and this section is wrong. The coordinates are in `docs/predictions/polished_*.xyz`
(`examples/export_mixed_cell.py`).

With the dispatch the E5 cells were built from Dipole's own record, atom for atom
(`docs/predictions/e5_structures/`, `examples/axial_registry.py`), and polyfind's numbers were
registered a second time in Dipole's cell, where its atoms sit within 0.010 A of the DFT's:
+1.53 for the half-monomer sublattice slip and +1.59 for one chain in four (6.3 per monomer of
the moved chain, 2.1 times the sublattice's), 14% and 7% under the values at polyfind's own
minimum because Dipole's b is 0.56% longer. Those are the like-for-like numbers for the single
points (`docs/predictions/e5_polyfind_belt4.json`; `docs/E1_E2_PREDICTIONS.md`, E5).

Two things the half-slipped lattice does in the model, found while registering (b), bear on 8.1.
At the fixed cell it is a saddle held by symmetry: released, the two sublattices turn 11.9
degrees about their axes in opposite senses and the cost falls from 1.79 to 1.44 kcal/mol per
monomer of the cell. With the cell free (the 0.51 of 8.1, which is the symmetric cell) the lowest
state has the chains deflected 7.4 degrees off trans and c 1.15% shorter, at 0.50. The
single-chain barriers of the table in 8.1 are rigid-lattice numbers, as E5's single points are;
with the lattice responding the half-slip of one chain costs 3.15, not 6.80 (8.1). The cheapest
DFT test of that is one relaxation of a cell the DFT seat already holds for E5 (c), with the
lateral cell free: **E5 d, dispatched by the PI on 2026-11-11**. polyfind's number for it is
registered for exactly that protocol (`docs/E1_E2_PREDICTIONS.md`, E5 d;
`docs/predictions/e5d_polyfind_belt4.json`): +0.708 kcal/mol per monomer of the cell against
1.668 rigid, a +5.1%, b +4.2%. If the DFT's relaxed number stays near its rigid one, the lattice
does not give way as the model says and the 3.15 of 8.1 goes back toward 6.8; if it comes out
near 0.45 of it, the mechanism stands.

### 8.7 The earlier form of this section (2026-11-08), for the record

*Superseded where 8.1-8.5 say so: the "upper bound for a dilute kink pair" (8.1), every kink
energy of a chiral row and the copolymer-to-terpolymer rate factor (8.2), and the reading of the
formation energy as the rate (8.5). The protocol comparisons between the search's own cells stand
as comparisons of those cells.*

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
| vdf9-trfe4-cfe1-14 (64/29/7, fourteen monomers) | +1.82 | - | +13.21 (stuck) | +11.75 (stuck) | 1.8-3.5: the pinned lattice reached from the kink's own conformation (`--pinned-from-own`, `kink_formation_belt4_14mer_pinned_from_own.json`) gives 3.52, the upper bound; the pinned and released values started from the all-trans placement are a stuck minimum (the ten-monomer CFE cell's pinned kink is 3.1 in the same 57 A^3 lattice) |
| vdf10-trfe4-14 (71/29 matrix, fourteen monomers) | +5.47 | - | +7.27 | +5.70 | 5.5-7.6 (pinned from own 7.56), the protocols agree: the matrix without the pin has the copolymer's kink cost; the one CFE per fourteen cuts it to 1.8-3.5 |

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
| 65/30/5, one CFE per twenty | +2.32 | 157 kHz | +0.055 | 0.60 | 0.154% | 0.199 / 0.383 | 1.29 | 0.09 / 0.05 |
| 71/29, the matrix without the pin | +5.47 | 0.8 kHz | -0.05 | 0.88 | 0.003% | 0.262 / 0.136 | 101 | 0.00 / 0.00 |
| 64/29/7 with CDFE in place of CFE | +1.08 | 1.3 MHz | +0.40 | 0.00 | 0.000% | 0.000 / 0.000 | - | 0.02 / 0.01 |
| one CFE and one CDFE per fourteen (57/29/7/7) | +0.86 | 1.8 MHz | +0.64 | 0.00 | 0.000% | 0.000 / 0.000 | - | 0.00 / 0.00 |
| 64/29/7 with AN in place of CFE | +1.97 | 285 kHz | -0.15 | 0.87 | 0.076% | 0.055 / 0.107 | 0.72 | 0.14 / 0.08 |
| 64/29/7 with CTFE in place of CFE | -0.57 | no meaning: the all-trans lattice is metastable against kinking | +0.135 | 0.12 | 0.060% | - | - | - |
| VDF/TrFE/CDFE 70/20/10, ten monomers | +0.61 | 2.8 MHz | +0.15 | 0.10 | 0.062% | 0.006 / 0.012 | 0.10 | 0.42 / 0.32 |
| VDF/TrFE 70/30, ten monomers | +4.95 | 2 kHz | +0.15 | 0.11 | 0.006% | 0.576 / 0.310 | 93 | 0.00 / 0.00 |

The matrix alone is a frozen ferroelectric at the drive: saturated at the bias (0.88) with its
channel at 0.8 kHz, so it carries loss and no strain at 34-66 kHz, the copolymer's known
behaviour. One CFE per fourteen does two things at once: it holds the matrix off saturation (polar
fraction 0.24) and it cuts the kink cost from 5.5-7.6 to 1.8-3.5 kcal/mol (own basin and pinned
lattice respectively), a factor of 460-900 in rate, which puts the channel at the drive. One CFE
per twenty does less of both: the margin falls to +0.055 and the strain rises to 0.154%, but the
kink costs 2.3 and the channel sits at 157 kHz, so at 34-66 kHz the loss is 0.20-0.38 and the loss
per strain 1.29, twice the one-per-fourteen cell's. The pin density therefore trades strain
against loss through the channel's position, and under the objective (loss first, stroke as a
floor) one per fourteen to sixteen is the model's optimum, not one per twenty. That is the relaxor terpolymer's mechanism in one cell, pin and
loosener in the same chlorine, and at the program's operating point the model gives tan delta 0.07
at 34 kHz for 0.125% AC strain, the largest strain of the table at the lowest loss per strain of
any row that switches. The ten-monomer CDFE row shows what a loosener without a pin does: a
cheaper kink still (0.6), a channel at 2.8 MHz far above the drive, a lower loss per strain (0.10)
but half the strain, because nothing holds its matrix on the steep part of the curve. CDFE at one
per fourteen in the 29 mol% matrix is the fastest lattice of the table (kink pair 1.1, channel 1.3
MHz) and switches nothing at the bias: at the film's dilution it pins harder than CFE (+0.40) and
overshoots the steep part, so the loosening it buys is not usable on its own. The cell with one
CFE and one CDFE per fourteen has the loosest lattice of all (kink pair 0.86) and is pinned harder
still (+0.64; section 5, item 4): in this model every chlorine buys loosening and pinning together,
and at the film's dilution one CFE per fourteen is the amount of chlorine that leaves the matrix on
the steep part of its curve.

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
