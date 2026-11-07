# Candidate chains and defect placements for the loss / frequency / strain optimum

*First pass, 2026-11-07, on what the model can say today: the crystal terms of
`pvdf-dft-valence-belt4` (docs/CRYSTAL_STAGE.md 5.22), the pendant-group relaxation of the
same section, the loss model of docs/ELECTROSTRICTION_RESEARCH.md 11-12 and the 2025-2026
literature pass. Code: `examples/candidate_screen.py` (the crystal screen),
`examples/candidate_loss_map.py` (the mapping onto the loss model); results
`docs/predictions/candidate_screen_belt4.json`, `candidate_loss_map_belt4.json`. Every number here
is a model number at zero kelvin on a regular ten-monomer periodic sequence packed as two chains;
section 6 says what that leaves out and what decides each row.*

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
| `Ea`, `Tf`, `tau0` | `f_relax`, hence loss per strain at 34-66 kHz and the heat | defect pinning (Ea), the size of the cooperative unit (Tf), chain stiffness | **not computed per chemistry yet** (tier 3); the bench's Vogel-Fulcher fit is the input |

and a fourth that is a property of the sequence rather than of any unit: the **spread of local
margins** along the chain, which the driven lattice (plan section 12) found to be the knob that
suppresses the stroke at bias and broadens the loss peak. A controlled sequence is the one place
where that spread is a design variable rather than a statistical accident, which is what the
microreactor programme is for.

## 2. The screen

Ten-monomer model sequences, each packed as a two-chain crystal in TTTT (the polar, beta-like
conformation) and TGTG' (PVDF's non-polar ground), under `pvdf-dft-valence-belt4` with every
CH2, CF2, CHF, CFCl and CHCl group relaxed as a rigid tilt (5.22). From the pair: the margin
`E(TTTT) - E(TGTG')` per monomer, the axial and areal transition strains, the dipole gained per
monomer (`P V / n` of the TTTT cell minus the TGTG' cell's) and `E* = dE / dmu`.

VDF and TrFE are the validated chemistry: their crystal terms were fitted to Dipole's PBE-D3
anchors and read back within 1.3 tolerances on everything anchored. CFE (CH2-CFCl), CTFE
(CF2-CFCl) and CDFE (CHCl-CF2) carry a chlorine whose nonbonded and charge parameters come from
the conformer valence fit alone; no crystal anchor has tested them, so those rows are
extrapolations. The E2 request to Dipole (docs/E1_E2_PREDICTIONS.md) is what validates them.

@SCREEN_TABLE@

@SCREEN_READING@

## 3. The same candidates through the loss model

`examples/candidate_loss_map.py` puts each candidate's margin, dipole gain and axial switch strain
into the calibrated literature channel (the 63/30/7 terpolymer's Vogel-Fulcher relaxation,
`tau0 = 8.08e-13 s`, `Ea = 0.124 eV`, `Tf = 191 K`, dielectric strength 45, halved under 60 MV/m;
plan section 11) and evaluates the program's operating point. The relaxation is the same for every
row, so the table separates what the chemistry changes in this model (the polar fraction at bias,
the strain per field, the heat at a given strain) from what it does not yet (the loss per unit
strain, which is `omega tau` and is the bench's number until tier 3 computes a barrier per
chemistry).

@LOSS_TABLE@

@LOSS_READING@

## 4. Defect type, concentration and placement: what the model says and what it cannot

@DEFECTS@

## 5. Candidates, ranked for the sequence-controlled route

@RANKING@

## 6. Confidence, and what decides each row

* **Validated against DFT:** the VDF and TrFE crystal terms (cells, polarization, stiffness,
  polar Berry slope, backbone angles; 5.22). **Not validated:** every chlorine parameter; the
  TrFE unit's own margins (the 70/30 model's conformational ordering is a prediction awaiting
  Dipole's E2 cells, 5.22); the polar-fraction mapping of section 3 (calibrated to one published
  film).
* **The model cell:** a regular, periodic, isotactic ten-monomer sequence packed as two chains at
  zero kelvin. Real copolymers are statistical and atactic, their crystals are lamellae a few
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
