# Electrostriction research plan: a better polymer, polyfind as the cheap surrogate, Dipole as the brute force

*2026-11-04. Written by the cloud polyfind session after the PI nominated it to continue the
research ("we're looking for a better polymer for electrostriction, understanding that tan
delta and frequency can trade off against strain, using/advancing polyfind to get us there as
needed; first make a pass on recent research; then drive Dipole to help run the brute-force
aspects and verify polyfind is accurate and can be used in lieu of weighty DFT"). The
literature pass is `sarco/materials/computational/papers/2026-11-04_electrostriction_literature_pass.md`;
the request to Dipole is `sarco/materials/ask-polyfind-dipole-electrostriction-verification-v1.md`;
the code this plan starts with is `polyfind.electrostriction`.*

## 1. The objective, as the program already states it

`sarco/materials/computational/objective_function_v1.md`: the plant is loss-limited. Minimise
tan delta subject to four floors: stroke rate `v = eps * lambda * f` (strain and drive
frequency trade along a hyperbola), `f_relax` above the drive frequency, enough crystallinity
for the recovery spring, and radial stroke `nu * eps`. Operating mode is bias-hold (50-65 V/um
DC, 1.5-3 V/um AC at tens of kHz, PI decision 2026-09-05). The frozen screening FOM is a
strain-per-loss proxy that encodes none of the floors, and two floor quantities, `f_relax` and
crystallinity, are outside the screen altogether.

"Better polymer" therefore means: more switchable polarization per monomer at the working
field (strain), with that polarization relaxing faster than the drive (frequency) and
dissipating less per cycle (loss), inside a lattice that still crystallises. The literature
pass says those three are one physical quantity seen three ways, and that quantity is
computable.

## 2. What the 2025-2026 literature changed (summary; the note has the sources)

1. **The giant strain is extrinsic.** The record strains (10.6% transverse in annealed
   P(VDF-TrFE-CFE) at 190 MV/m; 14.4% thickness with hydrogen-bonding polymer dots at 100
   MV/m; -3.5% at 50 MV/m with 1.9 mol% fluorinated alkyne) all come from a field-induced
   conformational transformation (TG / T3G to all-trans) that the crystal shares with the
   oriented amorphous fraction (OAF) and the taut-tie molecules between lamellae. The
   crystal-only ceiling for the alpha-to-beta transformation is 10.9% at full crystallinity;
   the measured strain exceeds the crystal share. Quasi-elastic neutron scattering puts the
   relaxor dispersion at the crystal-amorphous interfaces, not inside the crystals. The
   dipole *rotation* channel in the FE phase gives about 1.1% and Q33 < 3 m^4/C^2 (Bernholc's
   DFT in the Science 2022 paper), which is also what polyfind's crystal model computes.
2. **The strain lever is the margin of the polar conformation above the ground state, in the
   defect-bearing lattice.** Every successful route (alkyne or double-bond units, hydrogen-bond
   templating, dehydrochlorination to the FE critical endpoint, OAF-building anneals) brings
   the polar phase within reach of the field. The arithmetic: one VDF monomer dipole (about
   0.5 e.A) in 50 MV/m is 2.7 meV, 0.06 kcal/mol. A polar conformation more than a few
   hundredths of a kcal/mol per monomer above the ground state (per monomer of the
   cooperatively switching unit) does not switch at working fields.
3. **Loss, frequency and strain are the same relaxation.** The slow polarization channel that
   gives the giant strain is the lossy one, and it is gone above its relaxation frequency. For
   one Debye channel the loss per unit strain at drive frequency omega is omega*tau; at a
   fixed drive frequency the only lever on loss-per-strain is tau. Vogel-Fulcher fits give
   Ea = 0.12 eV for both the fluoropolymer terpolymer (Tf = 191 K, tau0 = 8e-13 s) and the
   fluorine-free liquid-crystal relaxor (Tf = 239 K); the 48 K in Tf is three orders of
   magnitude in f_relax at room temperature (about 0.4 MHz against a few hundred Hz). The
   cooperativity (Tf), i.e. how switching units interact, is the discriminating quantity.
4. **Coefficients are complex.** Measured M and Q carry a phase; the sign and magnitude depend
   on frequency; ab initio or model values compare only with proper coefficients. Every number
   this program quotes or computes should carry a frequency and, when measured, a phase.

## 3. The quantitative picture the plan is built on

**Two-state switch.** A unit of `n` monomers switches from a non-polar ground conformation to
the polar one at energy cost `n * dE` (kcal/mol per monomer, zero field) and gains `n * dmu`
of dipole along the field. The threshold field is `E* = dE / dmu` (independent of `n`), the
switched fraction is `x = 1 / (1 + exp(n (dE - dmu E) / kT))`, and the strain is `x` times the
transition strain of the conformational change (axial: the rise per monomer of the polar
chain over the ground chain minus one, +11% for TGTG' to TTTT in PVDF; transverse: the
change in cross-section per chain between the two packed cells). Three things follow. The
threshold is intensive and sets whether a chemistry works at 50-100 MV/m at all; it needs `dE`
to a few hundredths of a kcal/mol per monomer. The slope of S(E) at the bias point (the
bias-hold response) is `n`-dependent: a cooperative unit sharpens the switch. And the
transition strain is a geometric number polyfind already knows exactly.

**Relaxation.** The switch has a barrier `Eb` (kink nucleation and migration in a perfect
chain; a local TG/TT flip next to a defect or at the interface in a real relaxor) and an
attempt frequency `nu0`; interactions between units give the Vogel-Fulcher `Tf`. Then
`tau(T) = tau0 exp(Eb / k (T - Tf))`, `f_relax = 1 / 2 pi tau`, and for the switch channel
alone `tan delta(omega) = omega tau` below `f_relax`, with the strain amplitude falling as
`1 / (1 + omega^2 tau^2)` above it. The design target in these variables: `E* <= E_bias`,
`f_relax >> f_drive`, `dE` small but the ground state still crystalline.

**Intrinsic floor.** Above `f_relax` what remains is the lattice electrostriction `M^(0)`,
`Q^(0)`: the quadratic strain of the crystal with its chains rigid and its polarizability
electronic plus small rigid-body relaxation. This is the one part the crystal model computes
directly, and the first thing to verify against DFT because it is a clean, converged,
zero-kelvin second derivative on both sides.

## 4. Where polyfind stands (what is verified, what is not)

Verified against DFT or experiment, with the number:

* Lattice energetics on the fitted potential: alpha below beta by 4.5 kJ/mol per monomer
  (four studies: 2.6-6.5); held-out energy error 1.76 kcal/mol across ten chemistries
  (`docs/DFT_FIT.md`); polar/antipolar margin of all-trans PVDF agrees with periodic PBE-D3(BJ)
  to 0.02 kcal/mol per monomer on a quantity DFT puts at -0.026 (`docs/POLAR_MARGIN_DFT.md`);
  generic error bar 0.27 kcal/mol per monomer (`docs/SCREEN.md`).
* Cells: PE and alpha/beta/gamma-PVDF within about 7% of experiment.
* Dielectric: clamped-ion dielectric tensor within 3% of periodic DFPT with the polarizable
  model (`docs/ELECTROMECHANICS.md` 5.8).
* Piezoelectric: clamped-ion polar column agrees with the provider's Berry value to 0.001
  C/m^2; relaxed-ion polar column -0.24 against the provider's estimated -0.34 (their Born
  charges are 25% larger along the polar axis); film d33 through the exact laminate is 42-44%
  of the measurement, every sign right, and the rest is non-crystalline (5.13).
* Switching energetics (side branch `kink-switching`, on the PI's machine, not pushed): a
  kink row costs 14.4-15.4 kcal/mol per kink, a row hops one monomer over 8.22 kcal/mol
  (0.357 eV), harmonic TST gives about 5 MHz per hop at 300 K (`RESUME.md`, 2026-10-09).

Not verified or not computed: the isolated-chain RIS ranking still prefers the 3/1 helix
(DFT does too, in Bernholc's large cell, so this may be right for the wrong reason); the
packed alpha cell is predicted polar (delta); the C-F calibration failure (relaxed C-F 1.400
against 1.350 A) carries into every valence-path number; the transverse stiffness anomaly
(112 against 17 GPa) is open; **electrostriction was not computed at all** (section 9 of
`docs/ELECTROMECHANICS.md`: "the response reported is strictly linear in the field"); no
temperature, no amorphous or interfacial phase other than a passive laminate layer; the
chemistries the model expresses are VDF, VDC, CFE, CDFE, AN, VDCN, FANOME and PE
(`polyfind.polymers`), with N and O but not S in the nonbonded table, no sp2 backbone (no
double-bond unit), no alkyne.

Against this, the literature's levers are: `dE` to a few hundredths of a kcal/mol in a
defect-bearing lattice (the model is at 0.02 on one such margin and 0.27 generically), a
switching barrier (computed once, for a perfect crystal, three times the measured relaxor
Ea), and the OAF/interface (not modelled). So the research question for polyfind is not
"can it compute the giant strain" (the crystal part it can; the extrinsic part needs the
interface model of tier 4) but "can it compute the two numbers that decide the design, the
switching margin and the switching barrier, per chemistry, and be trusted".

## 5. The program, in tiers

Each tier names what polyfind computes, what Dipole computes to verify it, the acceptance
threshold, and the cost. Tiers 1-3 are crystal-level and sit on machinery that exists.

### Tier 1. Intrinsic electrostriction `M^(0)`, `Q^(0)` (this commit)

*polyfind:* `polyfind.electrostriction.electrostriction(ref, elastic)` returns the field
electrostriction tensor `M_iiJ` (nominal field along reference axis `i`, strain `J`) by two
routes that must agree because both are second derivatives of one fixed-voltage enthalpy:
the even part of the zero-stress strain at +E, 0, -E (the strain route), and half the mixed
derivative of the nominal polarization with respect to stress and field, from the state at
(+/- stress, +/- field) (the polarization route; the textbook `2M = d(eps0 eps)/d sigma` with
the *fixed-stress* permittivity, which differs from the fixed-strain one a DFPT run gives by
the piezoelectric term `e d / eps0`). `Q = M / (eps0 chi_free)^2`. The odd part of the same
strains is the converse piezoelectric coefficient and is checked against
`mechanics.piezoelectric`. Rigid or deformable path, fixed charges or polarizable, as the
reference is built. On the rigid, fixed-charge test crystal (illustrative potential,
force-shifted Lennard-Jones cutoff) the two routes agree to 5e-5 once the internal
relaxation is tightened: the polarization is linear in a leftover gradient of that
relaxation where the stress is quadratic, so the polarization route needs `gtol` five
orders tighter than the elastic constants do, and with the package's energy-shifted cutoff
the stress-targeted states land on either side of the force jump at the cutoff and the
route scatters by up to 3x. Both are now documented in the module.

*Dipole:* at strained, ion-relaxed cells (three strains along each of the polar, lateral and
chain axes) the relaxed-ion dielectric tensor by DFPT and the Berry polarization, for
beta-PVDF, PE (control), the delta/alpha cell and the P(VDF-TrFE) model cell, plus the
stiffness from the same strain set. The comparison is then made quantity by quantity on
both sides from the same inputs: `d eps_ii / d eta_J` (fixed strain), `e_iJ` from the
polarizations, `C_IJ`, and the `M` that follows from them with the same fixed-stress
conversion. Small cells, CPU, hours. Acceptance: polyfind `M` within a factor 1.5 of DFT on
the polar column, sign right on every column; PE's `M` small and of the lattice kind.

### Tier 2. Switching margin and transition strain per chemistry

*polyfind:* for each chemistry the packed, refined cells of the all-trans and the TG / T3G
conformations (`pack`, `refine_crystal`, existing), their energies per monomer, their
polarizations and their rise per monomer and cross-section per chain; from these `dE`, `dmu`,
the transition strain, `E*` and the two-state S(E) (`polyfind.electrostriction.TwoState`).
The screen (`docs/SCREEN.md`) already has the accessibility column for nine chemistries at
0.27 kcal/mol resolution; tier 2 asks for the same margins in the *defect-bearing* lattice
(one defect per `s` monomers, Dipole's pentad geometry) and at the resolution the switch
needs, which is what the polar/antipolar work reached.

*Dipole:* fully relaxed PBE-D3(BJ) cells of TTTT, T3GT3G' and TGTG' for PVDF (control: alpha
stable by a few tenths), the VDF/TrFE model cell, and the beta host with one CFE, CTFE, CDFE,
FA or DB unit; energy per monomer, Berry polarization, cell and `c` per monomer. Dipole's
GFN1 two-phase energies failed the absolute test (`field_response_v1.md` section 3: GFN1
over-binds the all-trans lattice by about 7 kcal/mol per monomer) and the MLIP route is
queued; PBE-D3 on 12-24-atom cells is hours on CPU. Acceptance: `dE` within 0.1 kcal/mol
per monomer of DFT on PVDF and the VDF/TrFE cell; ordering of the defect units preserved;
transition strain within 10%.

### Tier 3. Switching barrier, relaxation time, loss

*polyfind:* the kink-row hop (0.357 eV, perfect beta) is the starting number. Next: the same
climbing-image NEB for a local TG/TT flip beside one defect unit in a beta host and at the
end of a lamella (the chain leaving the crystal into a tether), on the fitted potential with
the complete placed-cell owner; a barrier distribution from the defect spacing; then
`tau(T)`, `f_relax`, `tan delta(omega)` from the Debye sum, and the VF form fitted to it so
it can be compared with the bench. This needs the `kink-switching` and `vector-repeat-merge`
branches, which exist only on the PI's machine.

*Dipole:* one DFT NEB on the smallest case (a 2x2-chain, 6-monomer beta supercell with one
CFE unit, local flip of one monomer pair) is days on CPU; the fine-tuned MLIP (566 PBE-D3
frames) is the intermediate and is Dipole's to run. Acceptance: barrier within 20% between
polyfind and the MLIP, MLIP within 20% of the one DFT case.

*Bench:* the small-signal eps'/tan delta(f, T) sweep already first on Dipole's bench list;
reported as VF parameters (`tau0`, `Ea`, `Tf`) so the computed barrier and cooperativity are
tested directly. This is also what decides the `f_relax` used by the operating-point
optimizer (assumed 30 kHz; the terpolymer's own VF fit says about 0.4 MHz at 300 K).

### Tier 4. The deformable connectors (OAF, tie molecules, interface)

The literature's extra factor of 1.5-2 in strain and the whole of the relaxor dispersion live
here. polyfind's film laminate (`polyfind.film`) has a passive amorphous layer; the model to
build is one where the tether layer switches too (a fraction of its monomers between TG and
TT under field, with the same `dE`, `dmu` and transition strain as the crystal but a
different barrier and no registry), so that the lamellar long period responds to field as
SAXS says it does (+10% at 150 MV/m) and the interface carries the loss. The
spontaneous-valency paper is the structural analogy: collective order from tether
deformability. No DFT verifies this tier; the SAXS and FTIR in-situ data of the ACS Nano paper
and the VF parameters do. Not started; it is the research, and it comes after tiers 1-3 have
the inputs it needs.

## 6. Candidates to run through tiers 1-3

In order of what the model can express today: PVDF (controls), CFE and CDFE units in a VDF
host (the baseline and its isomer), VDC, AN, VDCN, FANOME (the slate's leaders, already in
`polyfind.polymers`), the VDF/VDCN copolymer. To add: TrFE as a unit (one F on the CH2
carbon; two different substituents are expressible since phase 2), CTFE, the double-bond
unit (needs an sp2 backbone pair and a torsion profile across it), the alkyne (needs a
linear sp carbon pair) and the disulfonyl side chain (needs S in the nonbonded table and a
multi-atom pendant with its own torsions). The fluorine-free relaxor is worth a row because
it is a relaxor with no crystallinity (so it fails the recovery floor and the force column
as made: modulus about 1 MPa) but tests whether a dipole-pair side chain on a stiffer
backbone keeps the low `Tf`. *PI ruling, 2026-11-04: a test row only.* It is screened for
what the dipole-pair side chain does to the switching margin and the cooperativity on a
stiffer backbone; it is not a candidate as made (1 MPa modulus: blocking stress 0.04 MPa,
work density about 0.8 kJ/m^3, coupling about 0.002, no crystalline recovery spring, a
relaxation near a few hundred hertz at room temperature).

## 7. What Dipole is asked for (summary; the ask file has the specifications)

E1. Tier-1 references: DFPT dielectric tensors at nine strained, ion-relaxed cells for
beta-PVDF, PE, the alpha/delta cell and the P(VDF-TrFE) cell, with stiffness. CPU, hours.
E2. Tier-2 references: relaxed PBE-D3(BJ) cells, energies, polarizations and `c` per
monomer for the three conformations of PVDF, VDF/TrFE and the beta host with one CFE, CTFE,
CDFE, FA and DB unit. CPU, hours to a day.
E3. Tier-3: one MLIP NEB (local flip beside a CFE unit) and, if affordable, one DFT NEB on
the same path. Days.
E4. The bench sweep reported as VF parameters; the `f_relax` assumption re-examined.
Plus two housekeeping items: Study-IP WP0 is closed at the text level by the literature
note (section 1.4), and the libration column's reading as the relaxor response needs the
QENS caveat.

## 8. Immediate polyfind work

Done in this commit: `polyfind.electrostriction` (tier 1 routes, the two-state switch and
the transition strain), `tests/test_electrostriction.py`, `examples/electrostriction.py`;
the TrFE and CTFE units and a 70/30 VDF/TrFE model copolymer (`polyfind.polymers`:
`trfe`, `ctfe`, `vdf7-trfe3`, all chiral), `tests/test_trfe_ctfe.py`. First tier-1 numbers
on `pvdf-dft-valence` with induced dipoles and Ewald, rigid chain, force-shifted cutoff
(the consumer note of 2026-11-04 in `docs/REFERENCE_DATA_REQUEST.md` has the table):
beta-PVDF's polar column `Q_x,1` = +1.6 m^4/C^2 (the order of Bernholc's dipole-rotation
value, Q33 < 3), its lateral column +10.9 (the field turns the chains), PE's lattice value
-3.5 to -3.9 with an electronic `chi` of 1.9; the two routes agree to 2e-3 (beta) and 2e-5
(PE). Next, in order: tier-2 margins for the three conformations of PVDF, TrFE, CTFE, CFE,
CDFE and the VDF/TrFE cell at the polar/antipolar protocol's resolution; the cross-check
script that ingests Dipole's records (`examples/dft_crosscheck.py`, once the E1/E2 record
format is agreed); tier 3 once the two side branches are pushed.

## 9. Open questions for the PI

1. Push `kink-switching` and `vector-repeat-merge` from `E:\source\gh\polyfind` so the cloud
   session can continue the barrier work; they are not on the remote.
2. Priority of E1-E3 against D13/D19/D23: E1 and E2 are CPU-small and can run beside the GPU
   queue, as the earlier reference-data request argued; the ask says so but does not reorder
   anything.
3. ~~Whether the fluorine-free disulfonyl family should enter the slate as a row despite the
   force column.~~ Ruled 2026-11-04: a test row only (section 6).

## 10. The microreactor side: the same computations from the printer's end (2026-11-04)

The PI's reminder of 2026-11-04: this seat's work also supports the parallel microreactor
programme, the reversible and controllable synthesis at the ends of polymer chains
(`sarco/shared/ip/fluoro-ribosome-design-v1.md`; disclosure sections 10.16 and 10.23;
Dipole items D13 and D19-D25). The two programmes meet in the chain. The printer writes a
sequence and, through its torsion gate, a conformation; the muscle wants a sequence whose
polar conformation sits a few hundredths of a kcal/mol above its ground state and switches
fast. Tiers 2 and 3 of this plan, run per *printed sequence* rather than per random
copolymer, are the printer's product specification: which periodic spacing keeps relaxor
behaviour (Dipole's D8), which sequence holds a programmed torsion once the chain leaves
the tunnel (design 2.5(a): TrFE units lock trans, CFE units kink; the relaxor principle),
and what each such chain is worth as an actuator. The copolymer builder is the tool, the
70/30 VDF/TrFE cell added today is the first printable cell, and the CFE, CDFE and TrFE
units at a prescribed spacing in a VDF host are the next.

**What polyfind can say now, and what it cannot.** The bond-forming chemistry at the chain
end (depropagation D13, halogen-bond grip D19, block joining D20, oriented attack D22,
mid-chain fluorination D23) is electronic structure and stays with Dipole; the force field
supplies the conformational and lattice environment those steps happen in. Two of the
design's deciding numbers are of that kind, and `polyfind.chainend` computes them
(`examples/chain_end_gates.py`):

*Torsion gate (design 2.5, D21: at least 2 kT per bond at no more than 2 V/nm).* On the
fitted single-chain potential, an all-trans 8-mer, the bond entering the tunnel (15) and
the one inside it (14); the field couples to the transverse part of the dipole *change*,
0.5-0.7 e.A for a VDF end, not to the 0.44 e.A unit dipole:

| unit | bond | E_G - E_T (kcal/mol) | transverse dipole change (e.A) | field to write gauche (V/nm) | field to hold trans (V/nm) |
|---|---|---|---|---|---|
| VDF | 15 | +0.93 | 0.58 | 1.6 | 0.2 |
| VDF | 14 | -0.37 | 0.75 | 0.5 | 0.9 |
| TrFE | 15 | -1.42 | 0.55 | 0 | 2.0 |
| TrFE | 14 | -1.66 | 1.07 | 0 | 1.2 |
| CFE, CTFE | | not meaningful (frozen angles) | 0.5-0.7 | | |

Reading: for a VDF end the gate works within the limit in both directions, 1.6 V/nm to
write a gauche at the entering bond and 0.2 to hold trans. For a TrFE end the model puts
gauche 1.4 kcal/mol below trans at the end of an isolated chain, so holding trans costs
2.0 V/nm, at the limit; the potential is fitted to PVDF single chains, and the sign of this
margin is a transfer assumption the DFT provider's TrFE end models (D19) can check. The
chlorine-bearing ends carry the frozen-angle strain of `docs/CHEMISTRY_EXTENSION.md` phase 1
(all-trans CFE 25 kcal/mol above its gauche at this level), so their energies are not
numbers; the valence path with the angles relaxed is the next step for them.

*Registry ring and backtracking (design 2.4 and 10.2, D21 and D24).* Sliding one chain of
the refined beta-PVDF crystal along its axis through one repeat, cell held, fitted potential:
corrugation amplitude 0.48 kcal/mol per monomer (0.8 kT), the half-slipped lattice 0.48
above the registered one (the provider's GFN1 "registry" number for pure PVDF is +0.56,
`sarco/materials/computational/registry_gate_v1.md`, on a potential that over-binds), peak
slope 0.59 kcal/mol/A, i.e. 41 pN per monomer. Against the design: the ring's corrugation
of "no more than about 1 kT per notch (est.)" is confirmed in order for a commensurate
partner, and the D24 criterion of 20 pN to backtrack means a ring may hold less than half a
monomer's worth of crystalline contact, which is what the 2.3% mismatched armchair ring
does and the halogen-bond anchor is for. A BN ring or tube as the partner needs boron in
the nonbonded table (nitrogen is there); the machinery is the chain-pair one.

**Next on this side**, in order: the valence path for the chlorine-bearing ends; a
tube/ring partner for the corrugation; per-sequence tier-2 and tier-3 runs for the
printable cells (D8's question); and the held-end geometry with the anchor, which is
Dipole's.

## 11. Level 1 of the loss model: the derived small-signal tan delta (2026-11-04)

The PI's question of 2026-11-04: can tan delta be simulated or derived, since the target is
the optimum on the frequency-force-efficiency surface and that needs tan delta at each
performance level. The honest answer is that nothing in the record gave it (the one
earlier loss proxy, rigid chain rotation, was abandoned on 2026-10-02, and the provider's
own note says static margins do not bound loss without a dynamical model), but every input
of a derived model exists. Three levels were proposed: a closed-form relaxation model fed by
computed and bench-calibrated numbers (days); kinetic Monte Carlo on the sequence lattice
for the large-signal regime (weeks); molecular dynamics of the crystal-amorphous interface
for the dipole autocorrelation (months, GPU, Dipole's). The PI chose level 1, and it is
`polyfind.loss` (`examples/loss_level1.py`, `tests/test_loss.py`).

**The model.** Units of `n` monomers switch between a ground and a polar conformation
with dipole gain `dmu_u = n dmu` and cost `dE_u = n dE`; the polar fraction at nominal
field `E` is `x = 1 / (1 + exp((dE_u - dmu_u E) / kT))` and the channel's dielectric
strength is `N dmu_u^2 x(1-x) / (eps0 kT)`, so a degenerate unit gives the textbook
`N dmu_u^2 / (4 eps0 kT)` and the bias suppression `sech^2(dmu_u E / 2kT)`, and a unit
well above its ground state gives nothing at zero field. The relaxation time is
Vogel-Fulcher, `tau0 exp(Ea / k(T - Tf))`, tilted by the bias as `1 / cosh((dmu_u E -
dE_u) / 2kT)`. The permittivity is Debye or Cole-Cole on that time plus `eps_inf` and a
background loss; the strain of the channel follows the same spectral function with the
transition strain and `dx/dE` as its amplitude; the operating point is the per-cycle
bookkeeping of electrical loss (`pi eps0 eps'' E_ac^2`), mechanical loss, matched-load
work (`Y S_ac^2 / 4`), efficiency and heat (losses times frequency).

**Calibration, and which inputs are computed.** The bias suppression ratio gives the unit
dipole and hence `n` (with the computed `dmu` per monomer); the dielectric strength then
gives the switchable fraction; `tau0`, `Ea`, `Tf` come from the bench sweep's VF fit
(`loss.fit_vogel_fulcher` fits the peak temperature against frequency as the sweep
delivers it). On the computed side, `dmu` per monomer and the margin `dE` per sequence are
tier 2, the attempt time is the libration frequency (the 53-68 cm^-1 modes give about 1e-13
s, against the fitted 8e-13), the barrier is tier 3, and `Tf` has no computed value yet (the
polar/antipolar coupling between chains is 0.03 kcal/mol per monomer, far too small to be
the cooperativity; QENS puts it at the interface). A candidate is then the baseline channel
with its `dmu` scaled, its margin set, its barrier shifted and its cooperativity replaced
(`SwitchingChannel.scaled`).

**First numbers, on published inputs for the 63/30/7 terpolymer** (VF `tau0` 8.08e-13 s,
`Ea` 0.124 eV, `Tf` 191 K; dielectric strength 45; `eps` halving under 60 MV/m; polyfind's
`dmu` 0.45 e.A per monomer): the unit is `n` = 17 monomers (about 4 nm of chain, a
nanodomain), 49% of monomers switch, `f_relax` is 364 kHz unbiased and 491 kHz at 55 MV/m.

| f (kHz) | tan delta, Debye, 0 / 55 MV/m | Cole-Cole 0.3, 0 / 55 MV/m |
|---|---|---|
| 1 | 0.0025 / 0.0018 | 0.013 / 0.010 |
| 10 | 0.025 / 0.018 | 0.064 / 0.049 |
| 30 | 0.076 / 0.053 | 0.131 / 0.101 |
| 50 | 0.126 / 0.088 | 0.181 / 0.140 |
| 100 | 0.251 / 0.174 | 0.273 / 0.212 |

At the program's bias-hold point (55 V/um DC, 1.8 V/um AC, Y = 0.2 GPa, mechanical tan
delta 0.05) the channel gives an AC strain of 8e-4 and a blocking stress of 0.16 MPa at
every frequency below the relaxation, 32 J/m^3 of work per cycle, and:

| f (kHz) | tan delta | efficiency | heat (MW/m^3) |
|---|---|---|---|
| 1 | 0.002 | 0.57 | 0.02 |
| 10 | 0.018 | 0.33 | 0.7 |
| 34 | 0.060 | 0.16 | 5.9 |
| 50 | 0.088 | 0.12 | 12 |
| 66 | 0.115 | 0.09 | 21 |

Three readings, each to be tested by the sweep rather than believed. The optimizer's
`tan delta_0` of 0.0125 at 34-66 kHz is what this channel gives at 7 kHz; at 34-66 kHz it
gives 0.06-0.12, and the 4 MW/m^3 heat budget is crossed near 25 kHz, so on these inputs
the baseline's cooled bias-hold optimum sits near 20-25 kHz at tan delta 0.04-0.05 and an
efficiency near 0.2, not at 34-66 kHz. The assumed 30 kHz `f_relax` is a factor of ten
below this film's VF fit; the loss at the operating point is set by `omega tau`, so the
chemistry lever is `tau`, and a candidate with a three times shorter `tau` moves the whole
curve up by three times in frequency at the same loss. And the Arrhenius alternative with
polyfind's computed kink hop (0.357 eV, `Tf` = 0) lands in the same decade at 300 K (198
kHz against 364) but diverges from the VF form away from room temperature (24 against 0.2
kHz at 260 K), which is exactly what the sweep in temperature discriminates.

**What the bench must deliver for the model to be the film's and not the literature's**
(this sharpens E4 of the ask): the VF parameters from the peak temperatures at each
frequency; the dielectric strength and `eps_inf`; the small-signal `eps` under the
operating DC bias at 1-100 kHz (the suppression ratio); the DMA mechanical loss; and the
DC leakage, which this model does not contain. With those, `SwitchingChannel` is the
baseline, and the candidates are its scaled copies with tier-2 and tier-3 inputs.

## 12. Level 2 of the loss model: the driven lattice (2026-11-04)

The PI's "Go" after level 1. `polyfind.kmc` puts level 1's units on an `L x L` periodic
lattice (chain segments along `x`, neighbouring chains along `y`) with an Ising coupling
`J` and a margin per unit, a spread of margins being the random field that keeps an
interacting lattice from ordering into one domain (a relaxor rather than a ferroelectric).
Units flip singly at Arrhenius rates with the barrier top midway, under `E_dc + E_ac sin(wt)`
by rejection-free kinetic Monte Carlo. It is built from a level-1 channel and calibrated
twice by measurement: the barrier so that the lattice's own zero-field autocorrelation time
equals the channel's Vogel-Fulcher time (exact, since `exp(-Ea/kT)` multiplies every rate),
and the unit density so that its static susceptibility equals the measured dielectric
strength (the coupling's Curie-Weiss enhancement absorbed into the switchable fraction; the
dynamics do not depend on the density). Two estimators then read it. The small-signal
permittivity at any bias comes from the lattice's own fluctuations
(`IsingKMC.linear_response`, fluctuation-dissipation: static strength from the variance,
frequency dependence from the autocorrelation), which uses every fluctuation where a small
drive is noise-limited on a few hundred units; checked against the exact two-state
reference to 2% in strength, 12% in time and 10% in tan delta at and below the relaxation.
The large-signal response comes from driven runs (`IsingKMC.run`: trajectories, the
fundamental's `eps'`, `eps''`, tan delta, the loop area as the loss per cycle, the third
harmonic, the residual drift). `examples/loss_level2.py` runs the baseline terpolymer.

**What the lattice says on the first parameters** (20 x 20 units, `J` = 0.02 eV, margin
spread 0.08 eV, `n` = 17, `dmu_u` = 7.7 e.A, the channel's `tau0`; barrier by formula
0.261 eV, calibrated 0.274 eV; zero-bias strength 152 scaled to 45):

| bias (V/um) | polar fraction | chi_static | tau (s) | tan delta at 20 / 50 / 100 kHz | S_ac at 1.8 V/um | heat at 20 / 50 kHz (MW/m^3) |
|---|---|---|---|---|---|---|
| 0 | 0.56 | 36.5 | 4.0e-7 | 0.045 / 0.113 / 0.226 | 4.0e-3 | 13 / 45 |
| 10 | 0.77 | 44.6 | 1.1e-6 | 0.131 / 0.325 / 0.592 | 4.8e-3 | 26 / 91 |
| 20 | 0.90 | 11.3 | 2.7e-7 | 0.025 / 0.061 / 0.115 | 1.2e-3 | 1.6 / 6.5 |
| 30 | 0.93 | 4.1 | 9.1e-8 | 0.006 / 0.015 / 0.029 | 4.5e-4 | 0.21 / 0.84 |
| 55 | 0.97 | 1.4 | 4.0e-8 | 0.0013 / 0.0033 / 0.0066 | 1.6e-4 | 0.03 / 0.12 |

Three readings. First, the bias axis has a knee: below about 15 V/um the lattice is in its
switching range (polar fraction 0.5-0.8, strength 35-45, the slowest and lossiest dynamics
at 10 V/um where units sit at their thresholds); above it the units saturate, the strength,
the strain per AC volt and the loss all fall together, by 25x in strain and 400x in heat
between zero and 55 V/um. Second, the cycle efficiency at a given frequency is nearly
independent of bias (0.47-0.55 at 20 kHz, 0.34-0.46 at 50, 0.24-0.37 at 100 kHz, with the
assumed mechanical loss of 0.05 capping it near 0.6), because work and electrical loss
scale with the same susceptibility; what the bias moves is the stroke per volt and the
heat. Third, re-optimising the AC amplitude to the 4 MW/m^3 budget gives a broad optimum
in stroke rate across the bias axis (about 40-50 per second per unit lever at 20 V/um and
20-50 kHz, about 50 at 55 V/um and 100 kHz), so the bias-hold point is a choice of where on
the surface to sit, not a single optimum, and the lattice is the first object that gives
the surface with a derived loss rather than an assumed one.

**What is not yet calibrated, in order of consequence.** The lattice saturates ten times
harder than the film: at 55 V/um its strength is 4% of zero bias where the literature says
50% at 60 MV/m, so `n` (the unit size) must be recalibrated on the lattice itself rather
than by the level-1 sech^2 formula, and the knee will move up. Its temperature dependence
is flatter than the Vogel-Fulcher fit (1.4e-6 against 8.5e-6 s at 280 K, 3.6e-7 against
5.7e-8 at 320 K): the coupling is too weak, or the spread too wide, for the measured
cooperativity, and `J` against the sweep's slope in temperature is the calibration still to
do. The driven large-signal rows at this lattice size are noise-limited (the loss at the
program's amplitudes is below the statistical floor of 400 units; the loop areas come out
with either sign), so the small-signal column is the fluctuation-dissipation one and the
large-signal estimator needs larger lattices or many replicas. A 32 x 32 run was started
with the same parameters (zero-bias strength 108 before scaling); its result is the addendum below.

**Next**, in order: calibrate `n`, `J` and the spread jointly on the bench's three numbers
(the suppression under bias, the VF slope in temperature, the strength); larger lattices or
replica averaging for the large-signal loss; the margins drawn from a printed sequence
(tier 2's per-sequence `dE`), which is where the muscle's and the printer's questions become
one calculation; and level 3 for the interface, which this lattice does not contain.

**Addendum: the 32 x 32 run** (same parameters, 1024 units, 24 periods per driven row, 1.7
hours of CPU, `examples/loss_level2.py` without `--quick`). The larger lattice agrees with
the 20 x 20 one within a factor two in every column and moves no conclusion, so the shape of
the surface is not a finite-size effect. Calibration: barrier 0.281 eV (formula 0.261),
zero-field polar fraction 0.563, zero-bias strength 108 scaled to 45 (density factor
0.417); the fluctuation-dissipation strength at zero bias then reads 44.6 against the
measured 45.

| T (K) | lattice tau (s) | Vogel-Fulcher (s) | Arrhenius at 0.281 eV (s) |
|---|---|---|---|
| 270 | 4.1e-6 | 6.6e-5 | 1.4e-7 |
| 285 | 9.1e-7 | 3.6e-6 | 7.5e-8 |
| 300 | 4.9e-7 | 4.4e-7 | 4.2e-8 |
| 315 | 2.5e-7 | 8.9e-8 | 2.5e-8 |
| 330 | 1.5e-7 | 2.5e-8 | 1.6e-8 |

Over 270-330 K the lattice's apparent activation energy is 0.43 eV, the bare barrier's 0.28
and the VF fit's 1.0: the coupling and the random field add 0.15 eV of cooperativity, a
fifth of what the measured fit demands, so `J` must go up (towards the ordering point, which
the spread then has to hold off), and the calibration of `J` on the bench's slope is a
one-parameter fit with a clear target. At zero bias the lattice's spectrum is broader than
one Debye (its integral time 6.3e-7 s against the fitted 4.9e-7), which puts its tan delta
at 1.4x level 1's at every frequency (0.121 against 0.086 at 34 kHz, 0.223 against 0.166 at
66). At 55 V/um the strength is 2.2 against level 1's 24.7, the same saturation as on the
smaller lattice (a parameter, not a size).

Bias sweep in the program's band (fluctuation-dissipation, 1.8 V/um AC, `Y` 0.2 GPa,
mechanical tan delta 0.05; the two-state reading, see the correction below for the strain
columns):

| bias (V/um) | x | chi_s | tau (s) | tan delta 34 / 66 kHz | S_ac | efficiency 34 / 66 kHz | heat 34 / 66 kHz (MW/m^3) | S at the 4 MW/m^3 budget, 34 kHz |
|---|---|---|---|---|---|---|---|---|
| 0 | 0.56 | 42.6 | 4.3e-7 | 0.083 / 0.161 | 3.3e-3 | 0.44 / 0.35 | 23 / 65 | 1.4e-3 |
| 10 | 0.75 | 43.2 | 7.1e-7 | 0.138 / 0.266 | 3.3e-3 | 0.37 / 0.26 | 31 / 88 | 1.2e-3 |
| 20 | 0.86 | 19.0 | 4.0e-7 | 0.069 / 0.128 | 1.5e-3 | 0.34 / 0.24 | 7.1 / 21 | 1.1e-3 |
| 30 | 0.93 | 13.8 | 7.4e-7 | 0.119 / 0.209 | 1.0e-3 | 0.20 / 0.12 | 7.4 / 22 | 0.76e-3 |
| 45 | 0.97 | 3.5 | 9.1e-8 | 0.009 / 0.018 | 2.7e-4 | 0.31 / 0.21 | 0.29 / 0.93 | 1.0e-3 |
| 55 | 0.98 | 2.1 | 5.7e-8 | 0.004 / 0.008 | 1.65e-4 | 0.30 / 0.20 | 0.11 / 0.35 | 1.0e-3 |

The relaxation time is not monotonic in bias (slow near 10 and again near 30 V/um): with a
spread of margins, different sub-populations of units pass through their thresholds at
different biases, and each passage is a slow, lossy region. The driven rows at 55 V/um and
10 V/um AC (the amplitude the heat budget allows there) show the response still linear
(third harmonic 2-4.5% of the fundamental, drift below 1% per cycle, strain 9.3-9.8e-4 at
every frequency), so the fluctuation-dissipation reading extends to the budget amplitude at
the bias-hold point; their loss is at the noise floor even at 1024 units (tan delta
0.004-0.037 with scatter, the 1.8 V/um rows pure noise), which settles that the
large-signal loss wants replicas, not a larger lattice.

**A model limit found in reading the sweep, and its correction.** The lattice (and level
1's channel) is two-state: ground, or polar with its dipole along the drive axis. The real
unit has three states, ground, polar up and polar down. At zero bias the two polar states
are equally occupied, so the linear-in-field strain vanishes (the strain is at 2f, the
`Q P^2` of the literature) while the dielectric strength does not (both polar states
respond; in mean field `N dmu_u^2 x_polar / (eps0 kT)`, `1 / (1 - x)` times the two-state
value at the same polar fraction, which the strength calibration absorbs). In mean field the
three-state unit's linear strain per AC volt is the two-state's times `tanh(dmu_u E_dc /
kT)`, rising from zero to one over a bias scale `kT / dmu_u`, 34 V/um on these parameters,
while the loss is unchanged. So the two-state rows overstate the strain, the efficiency and
the stroke at budget at low bias, and the corrected surface is not broad but rising: stroke
at the 4 MW/m^3 budget and 34 kHz of 0 at zero bias, 12 per second per unit lever at 10
V/um, 20 at 20, 19 at 30, 30 at 45 and 31 at 55 (34 uncorrected), with the efficiency at
34 kHz 0, 0.07, 0.16, 0.12, 0.29 and 0.28 along the same biases: a plateau from about 45
V/um, and the program's bias-hold point of 50-65 V/um is on it. That mean-field estimate
compares units at the same polar fraction; once the three-state unit is calibrated on the
same bench numbers (zero-field polar fraction one half, margin `kT ln 2`, the same `n`), its
strain per field is the two-state's times `2 tanh(dmu_u E / 2kT)` instead, and the lattice's
own three-state reading, which supersedes both estimates, follows below. The knee and the
loss columns of the two-state tables stand.

**The three-state lattice** (same day, later; `IsingKMC(three_state=True)`, taken from a
`SwitchingChannel(states=3)`). Three things had to be settled before the lattice could carry
the third state. First, the unit's margin: a polar conformation with two dipole states has
twice the ground state's weight, so a unit with no margin is two thirds polar at zero field,
and on the lattice the conformational coupling took that to 92% polar, where the barrier
calibration ran into the rate floor. The calibrated channel now takes the zero-field polar
fraction (one half until the bench measures the all-trans content, E4 item 6) and sets the
margin that gives it, `kT ln 2` for three states; at one half the two- and three-state bias
suppressions coincide (`sech^2(dmu_u E / 2kT)`), so the same measurement calibrates the same
unit (`n` = 17, `dmu_u` 7.7 e.A), and the three-state channel's strain per field is the
two-state's times `2 tanh(dmu_u E / 2kT)`. Second, the calibration observable: the
polarization autocorrelation of the coupled three-state lattice is multi-exponential (units
at a conformational interface reverse fast, units inside a polar region slowly), so a fitted
time constant depends on the sampling step and the first calibration overshot threefold.
The lattice is now calibrated on the centroid in log frequency of its own loss peak
(`IsingKMC.loss_peak`), which is what the Vogel-Fulcher fit of a sweep is made on, scales
exactly with `exp(-Ea / kT)`, and is far less noisy than the maximum of a peak 1.8 decades
wide; the shift is iterated because the centroid of an off-centre peak is pulled toward the
grid's centre. Third, the strain estimator: the strain per field at a bias comes from the
cross-correlation of the polar fraction with the polarization (Kubo), beside the
autocorrelation that gives the permittivity, and the driven runs report the strain's second
harmonic. On the first lattice parameters (20 x 20, `J` 0.02 eV, spread 0.08 eV; barrier
0.289 eV, loss-peak width 1.8 decades, zero-field polar fraction 0.57, zero-bias strength
45.6 against the channel's 45 with no density scaling to speak of), 8000-sample
fluctuation-dissipation rows, 1.8 V/um AC, `Y` 0.2 GPa, mechanical tan delta 0.05:

| bias (V/um) | x | m | chi_s | tau_int (s) | tan delta 20 / 50 kHz | S_ac | efficiency 20 / 50 kHz | heat 20 / 50 kHz (MW/m^3) | stroke at budget 20 / 50 kHz | DC strain (switching population) |
|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 0.57 | 0.00 | 48 | 5.8e-7 | 0.067 / 0.165 | 1.1e-4 (noise) | 0.00 / 0.00 | 6.3 / 37 | 2 / 2 (noise) | 0 |
| 10 | 0.60 | 0.17 | 59 | 6.0e-7 | 0.071 / 0.176 | 9.8e-4 | 0.10 / 0.04 | 8.6 / 49 | 13 / 14 | 0.3% |
| 20 | 0.67 | 0.36 | 58 | 6.0e-7 | 0.071 / 0.176 | 1.9e-3 | 0.27 / 0.15 | 10 / 53 | 24 / 26 | 1.1% |
| 30 | 0.81 | 0.57 | 65 | 1.1e-6 | 0.123 / 0.291 | 2.5e-3 | 0.26 / 0.13 | 19 / 86 | 23 / 25 | 2.6% |
| 45 | 0.91 | 0.80 | 34 | 6.6e-7 | 0.074 / 0.181 | 8.6e-4 | 0.12 / 0.06 | 5.5 / 31 | 15 / 15 | 3.8% |
| 55 | 0.94 | 0.87 | 19 | 5.2e-7 | 0.054 / 0.132 | 3.6e-4 | 0.06 / 0.02 | 2.3 / 13 | 10 / 10 | 4.0% |
| 65 | 0.96 | 0.91 | 11 | 5.8e-7 | 0.053 / 0.129 | 2.3e-4 | 0.03 / 0.02 | 1.5 / 8.4 | 7.5 / 7.8 | 4.2% |

(`m` is the net polar fraction, the polarization over `N dmu_u`; the stroke at budget is
`S_ac f` with the AC amplitude re-scaled to the 4 MW/m^3 heat budget, per second per unit
of the clutch's lever; the DC strain is `s_switch (x - x_0)`, the quasi-static
electrostrictive strain of the switching population itself; the film's strain is that times
the population's share of the chain length, the channel's switchable fraction, 0.25 on this
calibration, since switched and unswitched segments lie in series along a chain.)

Readings. (a) The third state removes the saturation problem: the strength at 55 V/um is
45% of its zero-bias value (two states: 5%; measured: 50% at 60 MV/m), because the dipole
sign relaxes without a conformational change and the conformational coupling no longer
enhances the polarization fluctuations (zero-bias strength 48 before scaling against the
two-state lattice's 108-152). (b) The zero-bias strain is at twice the drive frequency, as
the real film's is: driven at 10 V/um AC the strain's second harmonic is 1.7e-3 and its
fundamental is at the noise floor. (c) The stroke-at-budget surface along the bias axis is
neither broad nor a plateau: it rises from zero to a maximum of about 25 per second per
unit lever at 20-30 V/um, where the net polarization is already 0.4-0.6 of full but the
polar fraction has not saturated, and falls to 10 at 55 V/um and 8 at 65, because the
coupled lattice switches cooperatively and is 94% polar by 55 V/um, leaving little to
switch; the efficiency follows the same shape (0.27 at 20 V/um and 20 kHz, 0.06 at 55).
Where that maximum sits is decided by the polar fraction's rise with the DC field, which
the lattice's DC strain column states and the bench's quasi-static strain curve `S(E_dc)`
measures directly. On these parameters the switching population's strain is 4% by 55 V/um
and flat beyond, which at a switchable fraction of 0.25 is about 1% of film strain, where
the neat terpolymer films reach 3.0-3.2% at 100 MV/m and are still rising: the lattice
saturates at too low a field, and its switching population or its strain per switch is too
small by about three times. Both are what `S(E_dc)` calibrates, the field scale through
`n` and the spread, the magnitude through the fraction and `s_switch` (E4 item 7, added
today), together with the suppression; the optimum will move up the bias axis with the field
scale, and the bias-hold point of 50-65 V/um is then a question for the calibrated lattice,
not a conclusion of this one. (d) The
loss at the bias-hold point is higher than level 1's: tan delta 0.09 / 0.20 / 0.30 at 20 /
50 / 100 kHz (level 1: 0.025 / 0.062 / 0.124), about 0.14-0.23 over 34-66 kHz, because
the biased lattice relaxes more slowly than the tilted well predicts (integral time 8.4e-7
s against 2.3e-7), and the 4 MW/m^3 budget at 1.8 V/um is crossed near 27 kHz (level 1 put
it near 25); the driven rows at 10 V/um AC, the amplitude the budget allows at 20 kHz, are
consistent with the small-signal column (tan delta 0.10 / 0.18 / 0.20, third harmonic 2-3%,
second harmonic of the strain 12-18%) and the 1.8 V/um rows are noise. (e) The temperature
dependence is as flat as the two-state lattice's: loss-peak times 1.2e-6, 3.6e-7 and
1.7e-7 s at 280, 300 and 320 K, an apparent activation of 0.39 eV against the
Vogel-Fulcher fit's 0.97 over the same window, and the first points of the coupling x
spread grid (`lattice_observables`, 270-330 K) put it at 0.36-0.37 eV for `J` of 0.02 and
0.035 eV alike, the lattice ordering into one conformation at `J` 0.035 with a spread of
0.08 eV (polar fraction 0.07, the strength then rising under bias as the field switches
the ordered lattice). The spread does move the suppression toward the measurement (0.34 at
a spread of 0.08 eV, 0.44 at 0.12, measured 0.50) without touching the slope. The
grid's verdict on whether this lattice can carry the Vogel-Fulcher slope follows.

**The coupling x spread grid, and what it settles** (three-state lattice, 20 x 20, `n` 17,
`lattice_observables` over 270-330 K; targets: apparent activation 1.0 eV from the VF
fit, suppression 0.50 at 60 MV/m):

| J (eV) | spread (eV) | barrier (eV) | zero-field polar fraction | strength | apparent activation (eV) | suppression at 60 MV/m | polar fraction at 60 MV/m | peak width (decades) |
|---|---|---|---|---|---|---|---|---|
| 0.020 | 0.08 | 0.289 | 0.56 | 44 | 0.36 | 0.34 | 0.95 | 1.6 |
| 0.020 | 0.12 | 0.284 | 0.54 | 41 | 0.37 | 0.44 | 0.80 | 1.6 |
| 0.035 | 0.08 | 0.261 | 0.07 (ordered, ground) | 7 | 0.37 | 1.2 (rises) | 1.00 | 2.6 |
| 0.035 | 0.12 | 0.234 | 0.60 | 46 | 0.37 | 0.22 | 0.97 | 2.4 |
| 0.050 | 0.08 | 0.271 | 0.00 (ordered, ground) | 0 | 0.27 | not switched at 60 MV/m | 0.00 | - |
| 0.050 | 0.12 | 0.209 | 0.07 (ordered, ground) | 6 | 0.37 | 2.2 (rises) | 1.00 | - |

The apparent activation energy is 0.36-0.37 eV at every disordered point and does not
rise with the coupling; past about `J` = 0.03 eV at these spreads the lattice orders into
one conformation, a ferroelectric with field-induced switching (its strength rises under
bias), not a stiffer relaxor. The 32 x 32 three-state run says the same over five
temperatures (loss-peak times 2.5e-6, 9.1e-7, 4.2e-7, 1.9e-7 and 1.2e-7 s at 270, 285,
300, 315 and 330 K: 0.39 eV apparent, peak widths 1.6-2.3 decades, the 300 K value within
4% of the calibration target). So this lattice, two-dimensional with a nearest-neighbour
conformational coupling, Gaussian margins and single-flip Arrhenius dynamics, does not
carry the Vogel-Fulcher slope: its temperature dependence is the barrier plus about 0.08 eV
of cooperativity whatever the coupling. The slope is an input from here on:
`IsingKMC.at_temperature(T, tau_target=tau_VF(T))` clones the lattice at `T` and
recalibrates its barrier to the bench fit's time there, and the lattice's deliverables are
the bias and amplitude dependence at each temperature, which is what the operating-point
question needs. What would carry the slope is a cooperative length growing toward a
freezing temperature, which is the interface and the chain connectivity of level 3, not
this lattice. The spread, by contrast, is a working knob: it moves the suppression toward
the measurement (0.34 to 0.44 for 0.08 to 0.12 eV at `J` 0.02) and broadens the peak; with
`n` it is calibrated on the suppression and the `S(E_dc)` curve when they arrive, `J`
staying below the ordering point. A size effect to carry along: at zero bias the 32 x 32
three-state lattice is twice as lossy as the 20 x 20 at the same calibration (tan delta
0.17 against 0.09 at 34 kHz; integral time 9.2e-7 against 4.3e-7 s), the slow tail of the
spectrum being the collective relaxation of polar regions that the smaller lattice
truncates; the biased rows of the 32 x 32 run (addendum below) do not share it.

**Addendum: the 32 x 32 three-state run** (same parameters; barrier 0.286 eV, loss-peak
time 3.9e-7 s against the 4.4e-7 target, zero-field polar fraction 0.55, zero-bias
strength 49 before scaling; it reached the end of the bias sweep in its two hours and was
stopped at the driven rows, for which the 20 x 20 rows above stand). Fluctuation-dissipation
rows, 1.8 V/um AC, `Y` 0.2 GPa, mechanical tan delta 0.05:

| bias (V/um) | x | m | chi_s | tan delta 34 / 66 kHz | S_ac | efficiency 34 / 66 kHz | heat 34 / 66 kHz (MW/m^3) | stroke at budget 34 / 66 kHz | DC strain (switching population) |
|---|---|---|---|---|---|---|---|---|---|
| 0 | 0.57 | 0.00 | 48 | 0.170 / 0.316 | 5.7e-5 (noise) | 0.00 / 0.00 | 26 / 84 | 1 / 1 (noise) | 0 |
| 10 | 0.59 | 0.17 | 50 | 0.154 / 0.279 | 1.1e-3 | 0.08 / 0.04 | 26 / 81 | 15 / 16 | 0.3% |
| 20 | 0.67 | 0.35 | 55 | 0.119 / 0.230 | 2.1e-3 | 0.23 / 0.14 | 26 / 84 | 28 / 30 | 1.1% |
| 30 | 0.78 | 0.55 | 52 | 0.146 / 0.271 | 2.0e-3 | 0.19 / 0.12 | 28 / 88 | 25 / 27 | 2.3% |
| 45 | 0.89 | 0.78 | 32 | 0.114 / 0.212 | 1.0e-3 | 0.12 / 0.08 | 13 / 44 | 20 / 21 | 3.6% |
| 55 | 0.94 | 0.87 | 20 | 0.083 / 0.159 | 7.2e-4 | 0.12 / 0.07 | 6.6 / 23 | 19 / 20 | 4.1% |
| 65 | 0.96 | 0.92 | 11 | 0.056 / 0.106 | 3.2e-4 | 0.06 / 0.03 | 2.6 / 9.3 | 13 / 14 | 4.4% |

Against the 20 x 20 lattice: the zero-bias loss is twice as high (the slow collective
tail), the loss at the bias-hold point is 30% lower (tan delta 0.061 / 0.147 against 0.087
/ 0.203 at 20 / 50 kHz; 2.3-2.5x level 1's), the strain per field there is twice as high
(3.6e-9 against 1.9e-9 per V/m; level 1 5.5e-9), and the 4 MW/m^3 budget at 1.8 V/um is
crossed at the same 26 kHz. The stroke-at-budget surface keeps its maximum at 20 V/um (28
per second per unit lever at 34 kHz) but its decline to the bias-hold point is gentle on
the larger lattice (19 at 55 V/um, against 10 on the 20 x 20), so the 50-65 V/um point
costs about a third of the maximum stroke on these parameters rather than six tenths, with
a third of the heat; the finite-size effect is in how fast the polar fraction saturates
(0.89 at 45 V/um here, 0.91 on the smaller lattice) and in the slow tail, and it is the
bench's `S(E_dc)` curve, not a larger lattice, that decides the position of the maximum.
