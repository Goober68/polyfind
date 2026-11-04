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
