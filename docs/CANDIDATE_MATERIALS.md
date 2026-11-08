# Candidate chains and defect placements for the loss / frequency / strain optimum

*First pass, 2026-11-07, on what the model can say today: the crystal terms of
`pvdf-dft-valence-belt4` (docs/CRYSTAL_STAGE.md 5.22), the pendant-group relaxation of the
same section, the loss model of docs/ELECTROSTRICTION_RESEARCH.md 11-12 and the 2025-2026
literature pass. Code: `examples/candidate_screen.py` (the crystal screen),
`examples/candidate_loss_map.py` (the mapping onto the loss model); results
`docs/predictions/candidate_screen_belt4.json`, `candidate_loss_map_belt4.json`. Every number here
is a model number at zero kelvin on a regular ten-monomer periodic sequence packed as two chains;
section 7 says what that leaves out and what decides each row, and section 6 reads Dipole's slate of
other chemistries against the same screen.*

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

| sequence | chemistry | margin TT-TG (rigid) kcal/mol per monomer | TGTG' tilt relaxation | axial, areal switch strain | P(TTTT) C/m^2 | mu(TTTT) e.A | E* (mu_TT) MV/m | E* (mu_TT - mu_TG) MV/m |
|---|---|---|---|---|---|---|---|---|
| pvdf | VDF/TrFE (validated) | +0.744 (+0.539) | -0.205 | +10.9%, -15.7% | 0.193 | 0.629 | 513 | 1208 |
| vdf7-trfe3 | VDF/TrFE (validated) | +0.145 (-0.172) | -0.357 | +10.6%, -17.2% | 0.167 | 0.560 | 112 | 324 |
| vdf8-trfe2 | VDF/TrFE (validated) | +0.177 (-0.127) | -0.343 | +10.5%, -16.4% | 0.174 | 0.578 | 133 | 448 |
| vdf6-trfe4 | VDF/TrFE (validated) | -0.329 (-0.735) | -0.441 | +10.3%, -16.7% | 0.162 | 0.546 | 0 | 0 |
| vdf5-trfe5 | VDF/TrFE (validated) | +0.752 (+0.526) | -0.253 | +10.9%, -19.6% | 0.160 | 0.537 | 607 | 1360 |
| vdf7-trfe2-cfe1 | chlorine (extrapolated) | +0.030 (-0.198) | -0.351 | +9.2%, -11.7% | 0.159 | 0.566 | 23 | 150 |
| vdf7-trfe2-cfe1-adj | chlorine (extrapolated) | +0.063 (-0.153) | -0.339 | +9.3%, -11.4% | 0.159 | 0.566 | 48 | 322 |
| vdf7-trfe2-ctfe1 | chlorine (extrapolated) | -0.098 (-0.383) | -0.390 | +8.9%, -10.1% | 0.136 | 0.506 | 0 | 0 |
| vdf7-trfe2-cdfe1 | chlorine (extrapolated) | +0.146 (-0.267) | -0.447 | +11.3%, -8.8% | 0.130 | 0.532 | 119 | 240 |
| vdf9-cfe1 | chlorine (extrapolated) | +0.976 (+0.747) | -0.310 | +11.3%, -14.8% | 0.176 | 0.614 | 689 | 1417 |
| vdf8-cfe2 | chlorine (extrapolated) | +0.835 (+0.506) | -0.416 | +10.6%, -13.1% | 0.147 | 0.536 | 675 | 1832 |
| vdf7-trfe2-an1 | chlorine (extrapolated) | -0.731 (-0.965) | -0.289 | +9.0%, -17.5% | 0.146 | 0.558 | 0 | 0 |
| vdf7-trfe2-fan1 | chlorine (extrapolated) | -0.572 (-0.572) | +0.000 | +9.1%, -13.1% | 0.144 | 0.581 | 0 | 0 |
| vdf9-an1 | chlorine (extrapolated) | +0.456 (+0.196) | -0.303 | +11.7%, -16.4% | 0.163 | 0.590 | 335 | 820 |

**Reading the screen.**

* **TrFE content sets the margin, in steps of about 0.3 kcal/mol per monomer per ten mol%.**
  PVDF's TGTG' is 0.74 below TTTT and its crystal threshold is 500 MV/m: the homopolymer cannot
  be field-switched, as is known (beta is made by drawing). At 20 and 30 mol% TrFE the margin is
  +0.18 and +0.15 and the threshold 130 and 110 MV/m, the field window in which the giant-strain
  copolymers and terpolymers are driven; at 40 mol% TTTT is the ground conformation. The crossover
  sits between 30 and 40 mol% in this model.
* **The pendant relaxation decides the sign of the small margins.** With rigid pendants the
  70/30 and 80/20 models had TTTT lowest (-0.17, -0.13); their TGTG' cells release 0.34-0.44
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
* **The switch strain is geometric and chemistry-insensitive:** +9 to +11% along the chain for
  every row (chlorine shaves it to 9%), with a lateral contraction of 10-20% in the model.
  The lateral number is exaggerated: PBE-D3 packs beta 3% too dense against experiment and the
  model inherits that from its anchors, so the volumetric change of the switch (6-8% here against
  2.5% between the measured alpha and beta densities) is an upper bound.
* **The dipole gained per monomer** is 0.54-0.63 e.A for the fluorine chemistries and falls with
  chlorine (0.51-0.57); the TTTT polarization falls from 0.193 (PVDF) through 0.167 (70/30) to
  0.159 (CFE), 0.136 (CTFE) and 0.130 (CDFE) C/m^2. Every defect costs stroke per switched unit.

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
| vdf7-trfe2-cfe1 | +0.030 | 23 | 0.77 | 0.110% | 0.073 / 0.141 | 13.4 / 47.4 | 0.13 / 0.08 | 37 / 72 |
| vdf7-trfe2-cfe1-adj | +0.063 | 48 | 0.56 | 0.153% | 0.088 / 0.170 | 22.3 / 77.5 | 0.15 / 0.09 | 52 / 100 |
| vdf7-trfe2-ctfe1 | -0.098 | 0 | 0.86 | 0.063% | 0.056 / 0.109 | 5.8 / 20.7 | 0.11 / 0.06 | 21 / 41 |
| vdf7-trfe2-cdfe1 | +0.146 | 119 | 0.10 | 0.062% | 0.047 / 0.091 | 4.1 / 14.7 | 0.14 / 0.08 | 21 / 41 |
| vdf9-cfe1 | +0.976 | 689 | 0.00 | 0.000% | 0.000 / 0.000 | 0.0 / 0.0 | 0.00 / 0.00 | 0 / 0 |
| vdf8-cfe2 | +0.835 | 675 | 0.00 | 0.000% | 0.000 / 0.000 | 0.0 / 0.0 | 0.00 / 0.00 | 0 / 0 |
| vdf7-trfe2-an1 | -0.731 | 0 | 0.88 | 0.062% | 0.053 / 0.102 | 5.7 / 20.4 | 0.10 / 0.06 | 21 / 40 |
| vdf7-trfe2-fan1 | -0.572 | 0 | 0.89 | 0.061% | 0.051 / 0.099 | 5.6 / 20.1 | 0.10 / 0.06 | 21 / 40 |
| vdf9-an1 | +0.456 | 335 | 0.00 | 0.000% | 0.000 / 0.000 | 0.0 / 0.0 | 0.00 / 0.00 | 0 / 0 |

**Reading the mapping.** The operating point is 55 V/um bias, 1.8 V/um AC; the polar fraction
at bias is where the crystal margin puts the channel on its own switching curve, and the strain
per field is largest where that fraction is near one half: with `dmu` about 0.57 e.A the
margin that does it at 55 V/um is about +0.07 kcal/mol per monomer. The CFE rows sit there
(+0.03, +0.06; polar fraction 0.77 and 0.56) and carry the largest AC strain (0.11-0.15%) and
stroke; the no-defect 70/30 and 80/20 and the CDFE row sit low on the curve (0.05-0.11) and give
half the strain; the TT-ground rows (60/40, CTFE) are saturated (0.86-0.88) and give the same
half; PVDF, the alternating 50/50 and the VDF/CFE copolymers do not switch at all. The loss per
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

**Type.** At one defect per ten monomers beside 20 mol% TrFE, the three chlorine units do
three different things to the margin: CFE (CH2-CFCl) takes the 80/20 copolymer's +0.18 to
+0.03 (it lowers the TGTG' side less than it lowers TTTT, i.e. it costs the all-trans packing
about 0.15 kcal/mol per monomer of the sequence, 1.5 kcal/mol per defect); CDFE (CHCl-CF2) leaves
it at +0.15 (the chlorine on the hydrogen-bearing carbon packs into the TTTT lattice as easily
as into TGTG'); CTFE (CF2-CFCl) overshoots to -0.10, TTTT becoming the ground state (its CF2-CFCl
pair prefers the trans packing). So in this model CFE is the defect that tunes a VDF-rich
copolymer onto the degeneracy, CDFE is inert for that purpose and only costs polarization, and
CTFE needs less TrFE beside it (a 75/15/10 rather than 70/20/10) to land in the same place.
These are the rows whose chemistry is extrapolated; the ordering CFE < CDFE in TT cost is the
one most worth a DFT cell.

**Concentration.** CFE without TrFE does not do TrFE's job: one or two CFE per ten VDF leave
TGTG' lowest by 0.8-1.0 kcal/mol per monomer, more than in PVDF. The two levers are different in
kind: TrFE lowers the all-trans conformation's energy chain by chain, CFE breaks the packing of
whichever conformation it sits in. The terpolymer recipe the literature converged on (about 30
mol% TrFE for the balance, 7-10 mol% CFE for the relaxor behaviour) is what the model reproduces
mechanistically. The CFE fraction that puts the margin on the operating point has only one value
here (10 mol% with 20 mol% TrFE); 5 and 15 mol% rows are the next thing to run.

**Placement.** Moving the single CFE from two monomers away from either TrFE to the position
next to one changes the sequence's margin by +0.03 kcal/mol per monomer, 0.3 kcal/mol per defect,
and the TTTT cell's energy by 0.10 per monomer. That is the size of the local-margin spread a
statistical terpolymer carries from its random placement (about 1 kT per 17-monomer unit, the
spread the driven lattice of plan section 12 found to suppress the stroke at bias and broaden the
loss peak), and it is the quantity a controlled sequence sets by choice: a strictly periodic
placement (D8) makes it zero and the switch as sharp as the chain's own cooperativity allows, at
the risk of ordering ferroelectrically when the coupling exceeds about 0.03 eV; a chosen pattern
of two placements sets a two-level spread. The model cannot say which spread is optimal; the
lattice model can, once the bench's `S(E_dc)` curve fixes its coupling.

**Location in the larger sense.** The screen puts every defect inside a crystalline chain. The
literature's finding is that the giant strain is produced at and around the crystals, in the
oriented amorphous fraction and the tie molecules, where the margin is lower than the crystal's
and the barrier is the local one. The model has no row for a defect at a lamellar surface or in a
tie segment; that is tier 4, and for the controlled-sequence route it is the question of where
along a chain of a given length the defects go so that the chain folds with them at the surface.

## 5. Candidates, ranked for the sequence-controlled route

For a sequence-controlled material aimed at the bias-hold operating point, on this model and
with the confidence stated in section 7:

1. **VDF/TrFE/CFE near 70/20/10, the CFE two monomers from either TrFE.** Margin +0.03, crystal
   threshold 20-50 MV/m, the largest strain per field at 55 V/um in the mapping; the statistical
   cousin of the 63/30/7 literature material. Chlorine chemistry: extrapolated until a DFT cell
   tests it.
2. **The same with the CFE next to a TrFE,** margin +0.06: the placement is the fine dial on the
   operating point and on the spread; both placements belong in the first controlled sequence
   so the dial is measured, not modelled.
3. **VDF/TrFE 70/30 to 75/25 without a defect.** Margin +0.15 to +0.18, threshold 110-130 MV/m:
   validated chemistry, but it needs twice the bias to reach the steep part of its curve and
   carries the ferroelectric copolymer's hysteresis, which this equilibrium picture does not
   contain (the lattice model's coupling does).
4. **VDF/TrFE/CTFE at 75/15/10.** CTFE is a stronger trans-favouring unit than CFE, so it wants
   less TrFE; its TTTT polarization (0.136) and switch strain (8.9%) are the lowest of the three
   defects, so it ranks below CFE unless its barrier or `Tf` is better, which nothing here says.
5. **CDFE** at 10 mol% does not move the margin and costs polarization: not favoured as the
   tuning defect.
6. **VDF/CFE without TrFE, PVDF, and the alternating 50/50:** not field-switchable in the crystal.

What would change the ranking: a barrier per chemistry (tier 3), which is the only lever on
loss per strain and which this pass does not have; Dipole's E2 cells for rows 1 and 3; and the
bench's Vogel-Fulcher and bias-suppression numbers, which replace the literature channel the
mapping of section 3 is built on.

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
| VDF/TrFE/AN 70/20/10 | -0.73 (-0.97): TTTT the ground state | -0.29 | +9.0%, -17.5% | 0.146 | 0.558 | 0 |
| VDF/TrFE/FAN 70/20/10 | -0.57 (-0.57; tilts not treated: the FAN carbon's pendants are a fluorine and a nitrile fragment) | - | +9.1%, -13.1% | 0.144 | 0.581 | 0 |
| VDF/TrFE/CFE 70/20/10, for comparison | +0.03 (-0.20) | -0.35 | +9.2%, -11.7% | 0.159 | 0.566 | 23 MV/m |

A nitrile defect beside 20 mol% TrFE makes TTTT the ground state by 0.6-0.7 kcal/mol per monomer:
the nitrile packs into the all-trans lattice's interchain gap (TTTT cell 58-62 A^3 per monomer,
gamma near 90 deg) and not into the TGTG' one (65 A^3, the cell sheared to 81-82 deg). That is the
opposite of CFE's action, which leaves the two conformations level. So a nitrile is not the
relaxor's tuning defect in a TrFE host, where it drives the material to a saturated ferroelectric
(polar fraction 0.88-0.89 at the 55 V/um bias in the mapping of section 3, half the CFE rows'
strain per field, tan delta 0.05); its place, if it has one, is in a VDF host without TrFE, where a
trans-favouring defect would bring PVDF's +0.74 toward the degeneracy from the other side. The
VDF/AN 90/10 row does that, but weakly: one acrylonitrile per ten VDF takes the margin from +0.74
to +0.46 kcal/mol per monomer (rigid +0.20; the TGTG' cell releases 0.30 through its tilts), the
threshold from 510 to 340 MV/m, and the TTTT polarization from 0.193 to 0.163 C/m^2. Per mole,
TrFE does the same job twice as well (10 mol% TrFE: +0.74 to +0.18) and keeps more polarization
(0.174), so a nitrile is the less efficient trans-favouring unit as well as the slower one; reaching
the degeneracy with AN alone would take 20-25 mol% on a linear extrapolation, at a polarization
cost of about 15%. Its polarization beside TrFE stays at 0.145 C/m^2, below CFE's 0.159. Dipole's registry gate reads AN as "marginal" (-0.11), a different quantity saying the
same thing from the lattice's side: the nitrile does not frustrate the polar crystal. Its
single-bond barriers do not settle its speed either way (host bond 10.2 against CFE's 6.8, defect
bonds' mean 9.6 against 12.2; see the barrier paragraph above).

**What this seat cannot compute.** CNEPO, SCN, SF5, nitro, isocyanide, the cyano-alkenes and
cyclopropyl need sulfur, nitro, epoxide and sp2 backbone parameters polyfind does not have; for
those Dipole's single-chain and GFN1 numbers are the only ones, and the barrier column above is
the one that matters for this programme's objective.

## 7. Confidence, and what decides each row

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
