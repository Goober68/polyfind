# Nitrile chains: a method gap in the RIS stage, not a fourth state, and not (on this potential) a chain without a minimum

**Verdict: the evidence supports the method reading, in a specific form — the frozen backbone
angles of the RIS scan, not a missing rotational state — with one reservation for the chemistry
reading that the screen has to carry.** Everything below is under one classical potential,
`pvdf-dft-valence`, fitted to PVDF data, whose every backbone torsion gets PVDF's CF2-CH2
Fourier triple and whose nitrile parameters were a minority of the reference set; the oligomer
path refuses charge flux by design, so the flux preset enters none of this.

What the four measurements say, in order of weight:

1. **There is no fourth torsional basin.** With the geometry relaxed, every bond adjacent to,
   one, two and three bonds from a nitrile carbon has the three canonical minima (a trans well
   deflected 10-20 deg off planar, and gauche wells at +/-50 deg) and nothing else. The only
   extra basin anywhere is a rigid-geometry well at +/-120 deg — in VDCN 8.7 kcal/mol *below*
   planar trans, behind a 13 kcal/mol barrier — and it disappears the moment the backbone
   angles are allowed to relax, relaxing to the +/-46 deg gauche minimum. That well is exactly
   what `fit_ris`'s frozen-angle scan folds into VDCN's "T" state: the fit places T at 180 deg
   and gives it the energy of the +/-120 deg conformer, which is why the screen reports an
   all-trans VDCN ground state 11.5 kcal/mol per monomer *below* the all-trans reference it is
   measured from.
2. **Relaxing the angles inside the RIS scan flips VDCN's conformational verdict** and leaves
   PVDF's and AN's qualitatively where they were: VDCN's all-trans goes from rank 1 of 83 to
   rank 74, 5.5 kcal/mol per monomer above a TG+ helix. The fix is on our side and it is one
   continuous freedom the package already has (`refine_crystal`'s per-atom backbone angles),
   applied one stage earlier.
3. **The single VDCN chain has a minimum as well-defined as PVDF's.** Perturbed from its own
   relaxed zigzag by up to 20 deg on every torsion, it returns eight times in eight to the same
   conformation within 0.001 kcal/mol per monomer, for every amplitude, where beta-PVDF already
   loses one run in eight at 20 deg; both scatter from 30 deg into gauche-containing chains
   1.3-2.4 kcal/mol per monomer below the zigzag. The strong chemistry reading — no minimum, a
   glassy ground state — is not what this potential shows. What looked like scatter from the
   RIS start was a 13 kcal/mol per monomer fall from the frozen-angle planar saddle, which is
   finding 1 again. The reservation that remains: VDCN's zigzag is not its own lowest periodic
   conformation (a chain with one gauche per four monomers is 0.6-2.4 below it, against
   PVDF's 1.5), the ordering among nitrile conformations moves by several kcal/mol with the
   angle treatment, and AN's all-trans neighbourhood has no well at all — so the nitrile rows
   of the conformational table cannot be read as a ranking before the potential's own error
   bar is even applied.
4. **The torsional correlation length does not discriminate** and is reported as ambiguous:
   the exact transfer-matrix length says VDCN is shortest (1.7 bonds at 300 K against PVDF's
   10.5), the ground-state-pattern run length says VDCN is longer (4.2 against 2.0), and both
   are properties of RIS models whose nitrile T state carries the artefact above.

**Recommendation for the screen.** Withdraw VDCN's row of the conformational table (its
"all-trans ground state, rank 1" is the rigid-scan artefact, not a result) and mark AN and
VDCN "not rankable by the rigid RIS stage as it stands", with FANOME staying "not screened"
for the reason it already carries. The freedom to add is the backbone bond angle, relaxed per
conformer inside `fit_ris`, with the state angles adapted from that relaxed scan; no fourth
state is warranted by anything here. Sarco's kinks 2-10 bonds from the nitrile are not
single-bond states of this potential at any distance we scanned, so if they are real they are
either cooperative or a property of the MLIP tier that this potential does not have; that is
the part of the chemistry reading this work cannot close.

Reproduce with `examples/nitrile_landscape.py` (`--part profiles | perturb | correlation |
risfit | all`, `--from-relaxed` for the perturbation control, `--json` to keep every number).
Wall time on one core: profiles 5 min, perturbation sweep 20 min (control 7 min), correlation
20 s, angle-relaxed RIS fits 11 min.

## The question

Four independent sightings (`docs/REFERENCE_DATA_REQUEST.md`, the 2026-09-11 exchange, and
`docs/GEOMETRY_GAPS.md`): sixteen of 192 dihedrals in Sarco's accepted VDCN/VDF reference more
than 30 deg from any rotational-isomeric state; five and fifteen such distortions in their two
basins from our all-trans starts, sitting 2-10 bonds from the nitrile rather than on it; the
polyacrylonitrile literature's loss of chain periodicity on relaxation regardless of tacticity;
and their finite-chain ladder, on which PVDF passes all six induced-dipole size gates and the
single-defect nitriles fail all six. Two readings: our three-state rigid-geometry search
cannot represent real minima (**method**: add a state or a continuous freedom), or nitrile
chains have no well-defined single-chain conformation (**chemistry**: no classical tier ranks
them, the honest output is "not rankable"). The consequences are opposite, so the question was
put to measurements rather than to argument.

## 1. Dihedral profiles: is there a basin at the kink angles?

One backbone dihedral of an eleven-bond oligomer is driven round the circle in 10 deg steps at
three levels. *Rigid*: every other bond trans and the polymer's frozen backbone angles — this is
`fit_ris`'s own 1D scan, reproduced to the digit. *Angles*: every backbone angle of the
oligomer relaxed (one per atom, bounded to 95-135 deg as `relax_backbone_angles` bounds them),
other bonds trans. *Full*: angles plus every other torsion relaxed, each grid point a local
minimisation from the angle-relaxed all-trans geometry so that the profile is one bond's with
the chain otherwise near its reference, not a hysteretic walk. Energies are kcal/mol for the
oligomer, relative to all-trans at the same level. Bonds scanned: PVDF's (control), both bond
types of VDCN and of AN (every bond of either homopolymer is adjacent to a nitrile carbon), and
in a VDF3-VDCN-VDF3 repeat the bonds 0, 1, 2 and 3 bonds from the nitrile-bearing carbon,
which is where Sarco's distortions sit.

| system, bond | rigid minima (deg: kcal/mol) | angles relaxed | fully relaxed |
|---|---|---|---|
| PVDF, CF2-CH2 (control) | 180: 0; +/-70: -1.94; +/-40: -1.29 | 180: 0; +/-50: -1.92 | 180: 0; +/-50: -3.32 |
| VDCN, both types | 180: 0; **+/-120: -8.7**; +/-30: -13.6 to -14.0 | 180: 0; +/-40: -3.9 to -4.5 | +/-160: -8.3; +/-50: -10.8 to -11.2 |
| AN, type 0 | -60: -11.9; +40: -6.6; +90: -10.0; **no minimum at 180** | -70: -5.2; +90: -3.6; no minimum at 180 | neighbours leave trans (drift 115-121 deg): not a one-bond profile |
| AN, type 1 | -90: -10.0; -40: -6.8; +60: -11.9; no minimum at 180 | -80: -3.8; +70: -5.1; no minimum at 180 | same (drift 113-131 deg) |
| copolymer, 0 bonds from CN (both) | 180: 0; **+/-120: +3.85**; +/-40: -1.91 | 180: 0; +/-50: -0.66 | +/-170: -0.20; +/-50: -2.59 |
| copolymer, 1 bond (both) | 180: 0; +/-90: -0.41; +/-40: +1.58 | 180: 0; +/-50: -1.14 | +/-170: -0.17; +/-50: -2.58 |
| copolymer, 2 bonds (both) | 180: 0; +/-70: -2.1; +/-40: -1.55 | 180: 0; +/-50: -2.2 | +/-170: -0.19; +/-50: -3.85 to -3.89 |
| copolymer, 3 bonds (both) | 180: 0; +/-70: -1.6; +/-40: -0.69 | 180: 0; +/-60: -1.38 | +/-170: -0.19; +/-50: -2.65 to -2.79 |

Bold entries are the only minima more than 30 deg from every canonical state (180, +/-60).
Both are rigid-geometry wells and both vanish at the next level. The basin test makes the same
point without a grid: relaxed with nothing held, VDCN's +/-120 deg rigid minimum goes to
+/-46 deg (the gauche well) and the copolymer's to +/-50 deg; PVDF's rigid minima go to 180 and
+/-51 deg. No relaxation of any rigid minimum ends in the kink region.

Three things the table says beyond the headline.

* **The frozen angles are the problem, and for VDCN they are a large one.** Relaxing the
  backbone angles of the all-trans oligomer takes VDCN from 255.7 to 176.2 kcal/mol (angles
  120 deg at CH2, 103 deg at C(CN)2), AN from 165.8 to 127.4 (117/103), PVDF from -80.8 to -87.5
  (115/107) — seven, three and half a kcal/mol per monomer. `docs/VALENCE_FIT.md` section 4
  measured the same thing; what was not measured there is that at rigid angles the *shape* of
  the nitrile profiles is wrong, not just their offset: the +/-120 deg well in VDCN is deeper
  than the trans well by 8.7 kcal/mol, and `fit_ris` cannot tell a deflected basin edge from a
  state. Its `argmin1` for VDCN's T state is -120 deg on one bond type and +120 on the other;
  the circular mean of those is 180, so the adapted T angle is 180 while `e1[T]` is -8.7. The
  screen's VDCN "TT ground state at -11.48 kcal/mol per monomer" is the sum of two such
  energies, and it describes a chain built at 180 deg with the energy of one at +/-120.
* **AN's trans state is not a minimum for a single bond at either level** — the profile
  descends monotonically from 180 to the gauche well on one side. The rigid fit again labels
  the basin edge (+/-120 deg) as T. At the fully relaxed level the neighbours leave trans
  entirely (115-131 deg of drift) at almost every grid point, so an all-trans AN
  neighbourhood is not a stable one in this potential; the fit's "T" is a state the chain does
  not have.
* **With an isolated VDCN unit the profiles are clean three-state at every distance**, and the
  gauche well is at +/-50 deg, not +/-60. The +/-120 rigid well is shallow (0.21 kcal/mol
  barrier) and 3.85 above trans there — it is the homopolymer's dense nitriles that make it
  deep. Nothing at 1, 2 or 3 bonds from the nitrile carbon differs from PVDF's own profile by
  more than the gauche well's depth (-2.6 to -3.9 against -3.3), which is the opposite of a
  nitrile-induced state.

## 2. Perturbation returns: is there a well-defined minimum?

A single periodic chain (one chain in a 30 A cell with an 8 A cutoff, so it sees only its own
axial images) starts from the RIS ground state at ideal angles, its torsions are perturbed by
Gaussian noise of width sigma, and `refine_crystal` relaxes torsions, backbone angles and cell
together with the fitted bend terms (`valence=`). Randomised torsions follow no line-group
pattern, so every run takes the refinement's documented `"free"` (penalty) path — the
parametrisation is the same for every endpoint. Eight seeds per sigma. Endpoints are compared
under the chain's own symmetries (shift by a repeat, reversal, mirror for an achiral repeat): a
cluster is a set of endpoints within 20 deg of each other on every torsion; "returned" means
within 20 deg of the reference minimum, the lowest endpoint at sigma = 2; a kink is a torsion
more than 30 deg from every canonical state, Sarco's criterion. VDCN's RIS ground state is TT;
PVDF's fitted one is a helix that is not commensurate within eight periods and cannot be built
as a repeat, so PVDF runs twice, all-trans (beta, the like-for-like control) and TG+TG- (alpha,
the lattice ground state). AN all-trans is included with the caveat that TT is not AN's RIS
ground state.

A run whose repeat transform is still a rotation of more than 0.5 deg at the end is a failed
closure of the penalty method (its energy is 80-90 kcal/mol per monomer and it is not a chain);
those are counted as "unclosed" and excluded. Energies are kcal/mol per monomer; the start is the
ideal-angle RIS chain before any relaxation.

| chain | sigma (deg) | closed / unclosed | distinct endpoints | returned | E min .. max (spread) | lowest vs reference | kinks per repeat, mean (max) | lowest endpoint |
|---|---:|---:|---:|---:|---|---:|---:|---|
| **VDCN TT**, start +25.02, reference +11.855 (deflected zigzag, +/-153 to +/-160) | | | | | | | | |
| | 2 | 8 / 0 | 2 | 88 % | +11.85 .. +12.14 (0.28) | +0.00 | 0.0 (0) | TTTTTTTT |
| | 5 | 8 / 0 | 3 | 38 % | +11.22 .. +12.14 (0.92) | -0.64 | 0.0 (0) | TG-TTTG+TT |
| | 10 | 8 / 0 | 3 | 50 % | +11.21 .. +12.14 (0.93) | -0.64 | 0.0 (0) | TTG-TTTG+T |
| | 20 | 8 / 0 | 3 | 50 % | +11.85 .. +15.53 (3.68) | -0.00 | 0.2 (2) | TTTTTTTT |
| | 30 | 6 / 2 | 4 | 33 % | +11.23 .. +15.03 (3.79) | -0.62 | 0.3 (2) | TTTG-TTTG+ |
| | 45 | 8 / 0 | 6 | 38 % | +11.03 .. +18.56 (7.54) | -0.83 | 1.1 (4) | G-TG+TTTTT |
| | 60 | 7 / 1 | 5 | 14 % | +10.68 .. +11.88 (1.20) | -1.18 | 1.0 (3) | G+G+KG-KTTK |
| **PVDF TT** (beta), start +5.62, reference +5.196 (deflected zigzag, +/-164 to +/-167) | | | | | | | | |
| | 2 | 8 / 0 | 2 | 88 % | +5.20 .. +5.25 (0.06) | +0.00 | 0.0 (0) | TTTTTTTT |
| | 5 | 8 / 0 | 2 | 38 % | +5.20 .. +5.25 (0.06) | +0.00 | 0.0 (0) | TTTTTTTT |
| | 10 | 8 / 0 | 2 | 75 % | +5.20 .. +5.35 (0.16) | +0.00 | 0.0 (0) | TTTTTTTT |
| | 20 | 8 / 0 | 1 | 75 % | +5.20 .. +5.37 (0.17) | +0.00 | 0.0 (0) | TTTTTTTT |
| | 30 | 8 / 0 | 5 | 50 % | +4.27 .. +6.87 (2.60) | -0.93 | 0.5 (2) | TTG-TG+TTT |
| | 45 | 8 / 0 | 7 | 12 % | +3.89 .. +7.76 (3.87) | -1.30 | 0.8 (2) | G+TTTG-TTT |
| | 60 | 8 / 0 | 6 | 0 % | +3.79 .. +9.55 (5.76) | -1.41 | 0.8 (3) | TG-TG-TG-TG- |
| **PVDF TG+TG-** (alpha), start +5.65, reference +3.969 | | | | | | | | |
| | 2 | 8 / 0 | 1 | 100 % | +3.97 .. +4.00 (0.03) | +0.00 | 0.0 (0) | TG+TG-TG+TG- |
| | 5 | 8 / 0 | 1 | 100 % | +3.97 .. +4.01 (0.04) | +0.00 | 0.0 (0) | TG+TG-TG+TG- |
| | 10 | 5 / 3 | 1 | 100 % | +3.97 .. +4.01 (0.04) | +0.00 | 0.0 (0) | TG+TG-TG+TG- |
| | 20 | 6 / 2 | 1 | 100 % | +3.97 .. +4.20 (0.23) | -0.00 | 0.0 (0) | TG+TG-TG+TG- |
| | 30 | 6 / 2 | 1 | 100 % | +3.97 .. +4.01 (0.04) | -0.00 | 0.0 (0) | TG+TG-TG+TG- |
| | 45 | 8 / 0 | 5 | 50 % | +3.74 .. +8.86 (5.12) | -0.23 | 0.8 (4) | TG+TG+TG+TT |
| | 60 | 7 / 1 | 5 | 0 % | +3.29 .. +6.20 (2.90) | -0.68 | 0.1 (1) | TG-TG-TG+TG+ |
| **AN TT**, start +15.27, reference +4.908 (G+KG-KTG-KG+: itself kinked) | | | | | | | | |
| | 2 | 8 / 0 | 8 | 12 % | +4.91 .. +6.47 (1.56) | +0.00 | 1.5 (3) | G+KG-KTG-KG+ |
| | 5 | 8 / 0 | 6 | 0 % | +5.06 .. +6.46 (1.40) | +0.15 | 1.2 (2) | KG-TKG-TG+G+ |
| | 10 | 8 / 0 | 5 | 0 % | +4.93 .. +8.34 (3.41) | +0.02 | 0.8 (2) | G+TG-TTG-TG+ |
| | 20 | 8 / 0 | 8 | 0 % | +4.92 .. +8.62 (3.70) | +0.01 | 1.9 (5) | G-TG+TTG+TG- |
| | 30 | 8 / 0 | 6 | 0 % | +5.47 .. +8.73 (3.27) | +0.56 | 2.0 (4) | KTG-TG+TG-T |
| | 45 | 8 / 0 | 5 | 12 % | +4.90 .. +8.98 (4.08) | -0.01 | 2.5 (3) | G-KG+TKG+KG- |
| | 60 | 8 / 0 | 8 | 0 % | +5.01 .. +8.61 (3.60) | +0.10 | 1.8 (4) | TTTG+TG+TG+ |

What it shows.

* **Alpha-PVDF is a well.** One endpoint, 100 % return, spread under 0.25, no kinks, up to 30 deg
  of noise on every torsion; only at 45-60 deg does it leave for other TG-type helices.
* **Beta-PVDF is a shallow well of a family**: two deflection patterns 0.06 apart, no kinks and
  at most two endpoints to 20 deg; from 30 deg it finds gauche-containing chains 0.9-1.4 below
  the zigzag, which is the known result that the planar zigzag is not PVDF's single-chain
  ground state, and its kinked endpoints are never its lowest.
* **VDCN's ideal-angle zigzag is left at 5 deg.** From 5 deg of noise three of eight runs
  land on a chain with one +/-87 deg gauche defect per four monomers, 0.64 below the deflected
  zigzag; from 20 deg kinked chains (torsions 30 deg or more from every state) appear, up to
  four per eight-bond repeat, and at 60 deg the *lowest* endpoint is a kinked one
  (G+G+KG-KTTK, 1.18 below the zigzag). The energy spread at 45 deg is 7.5 against
  beta-PVDF's 3.9. Read alone this would be the chemistry reading; the control below says it
  is the start, not the well.
* **AN all-trans is not a minimum at all**: 2 deg of noise gives eight distinct endpoints from
  eight seeds, all within 1.6 kcal/mol per monomer, the lowest of them kinked. This is the
  profile of section 1 seen from the other side — there is no trans well to return to — and it
  is also not AN's RIS ground state, so it says the all-trans neighbourhood is rugged, not that
  AN has no minimum.

The ideal-angle start is 13 kcal/mol per monomer above VDCN's relaxed zigzag and 0.4 above
PVDF's, so VDCN leaving at 5 deg where PVDF holds to 20 could be the height of the start rather
than the shallowness of the well. The control perturbs each chain from its own relaxed minimum:

| chain (perturbing its relaxed zigzag) | sigma (deg) | closed / unclosed | distinct endpoints | returned | E min .. max (spread) | lowest vs reference | kinks, mean (max) | lowest endpoint |
|---|---:|---:|---:|---:|---|---:|---:|---|
| **VDCN TT**, relaxed zigzag +11.855 (+/-153 to +/-161) | 2 | 8 / 0 | 1 | 100 % | +11.85 .. +11.86 (0.00) | +0.00 | 0.0 (0) | TTTTTTTT |
| | 5 | 8 / 0 | 1 | 100 % | +11.85 .. +11.85 (0.00) | -0.00 | 0.0 (0) | TTTTTTTT |
| | 10 | 8 / 0 | 1 | 100 % | +11.85 .. +11.86 (0.00) | +0.00 | 0.0 (0) | TTTTTTTT |
| | 20 | 8 / 0 | 1 | 100 % | +11.85 .. +11.85 (0.00) | +0.00 | 0.0 (0) | TTTTTTTT |
| | 30 | 8 / 0 | 4 | 50 % | +11.21 .. +12.17 (0.96) | -0.64 | 0.4 (3) | TTG-TTTG+T |
| | 45 | 8 / 0 | 6 | 12 % | +10.37 .. +12.50 (2.13) | -1.49 | 0.5 (4) | G+TTTG-TTT |
| | 60 | 7 / 1 | 6 | 0 % | +9.47 .. +19.77 (10.30) | -2.39 | 1.3 (3) | G-G-G-KG+G+G+K |
| **PVDF TT**, relaxed zigzag +5.196 (+/-164 to +/-167) | 2 | 8 / 0 | 1 | 100 % | +5.20 .. +5.20 (0.00) | +0.00 | 0.0 (0) | TTTTTTTT |
| | 5 | 8 / 0 | 1 | 100 % | +5.20 .. +5.20 (0.00) | -0.00 | 0.0 (0) | TTTTTTTT |
| | 10 | 8 / 0 | 1 | 100 % | +5.20 .. +5.20 (0.00) | -0.00 | 0.0 (0) | TTTTTTTT |
| | 20 | 8 / 0 | 2 | 75 % | +5.20 .. +5.71 (0.52) | +0.00 | 0.1 (1) | TTTTTTTT |
| | 30 | 8 / 0 | 8 | 12 % | +3.53 .. +6.48 (2.94) | -1.66 | 1.1 (3) | G+TG-TG-TG+T |
| | 45 | 8 / 0 | 6 | 25 % | +3.89 .. +9.02 (5.13) | -1.30 | 0.4 (1) | G+TTTG-TTT |
| | 60 | 8 / 0 | 6 | 0 % | +3.73 .. +5.69 (1.96) | -1.47 | 0.5 (2) | G+G+TG-G-TTT |

**This is the discriminating row of the whole document.** From its own relaxed minimum, VDCN's
zigzag returns eight times out of eight to the same conformation, to a thousandth of a kcal/mol,
for every perturbation up to 20 deg on every torsion — exactly as PVDF's does (PVDF already
loses one run in eight at 20 deg; VDCN loses none). Both leave from 30 deg, both find
gauche-containing chains 1.3-2.4 kcal/mol per monomer below their zigzag, both have kinked
local minima above it. The VDCN single chain has a minimum as well-defined as PVDF's; what made
it look frustrated from the ideal-angle start was a 13 kcal/mol per monomer fall from a
frozen-angle saddle, and that is section 1's finding again. The one residual asymmetry is the
spread at 60 deg (10.3 against 2.0), carried by a single high-energy kinked endpoint, and the
fact that VDCN's zigzag is not its own lowest periodic conformation by 2.4 kcal/mol per monomer
where PVDF's is not by 1.5 — a difference of degree.

## 3. Torsional correlation length: is it about periodicity?

The screen's own fitted RIS models (`fit_ris`, step 10, six monomers, third order) sampled
exactly with `polyfind.amorphous` (2000 chains of 400 bonds), and two readings of the
correlation length: exact, from the second eigenvalue of the repeat-unit transfer matrix
(`xi`), and the run length over which a sampled chain follows its own RIS ground-state pattern
at one phase (pooled over phases; "P(>= 8)" is the fraction of bonds in such a run of at least
eight bonds, a crystallographic repeat's worth).

| polymer (ground state) | T (K) | xi (bonds) | T fraction | trans run | pattern run (bonds) | P(>= 8) |
|---|---:|---:|---:|---:|---:|---:|
| PVDF (TG+TG+G+TG+G+) | 300 | 10.5 | 0.40 | 1.04 | 2.0 | 0.28 |
| | 400 | 5.8 | 0.41 | 1.08 | 2.0 | 0.21 |
| VDCN (TT) | 300 | 1.7 | 0.80 | 4.2 | 4.2 | 0.20 |
| | 400 | 1.8 | 0.80 | 3.9 | 3.9 | 0.17 |
| AN (TG-) | 300 | 20.9 | 0.50 | 1.05 | 20.9 | 0.46 |
| | 400 | 9.0 | 0.49 | 1.10 | 8.7 | 0.39 |

Ambiguous, and reported as such. By `xi` VDCN is the shortest by a factor of six against PVDF;
by the pattern run it is twice as long, because PVDF's fitted ground state is an eight-bond
helix the chain almost never follows for a full period (its T/G+/G- fractions are 0.40/0.30/0.30
at 300 K — close to random) while VDCN's is TT and the chain is 80 % trans. AN's one-handed TG-
helix persists for 21 bonds. None of the three metrics separates "cannot crystallise" from
"can": VDCN's 4-bond trans runs are shorter than a repeat but its `P(>= 8)` is within a third
of PVDF's. And every number here is a property of a model whose nitrile T state carries the
artefact of section 1 (VDCN's `e1[T]` is -8.7 rather than 0), so the VDCN row measures the
model rather than the chain.

## 4. The RIS fit with the angles relaxed: is the fix on our side?

The screen's `fit_ris` call, repeated with a calculator that relaxes the repeat's backbone
angles (one per backbone atom of the repeat, broadcast along the oligomer — the same freedom
`refine_crystal` has — bounded to 95-135 deg) for every conformer the scan scores. Everything
else is the screen's: step 10 deg, six monomers, third order, `enumerate_periodic` to period 8.

| polymer | fit | state angles | e1(T) / e1(G) per bond type | ground state (per monomer) | all-trans rank, above ground |
|---|---|---|---|---|---|
| PVDF | rigid | 180, +/-70 | 0 / -1.88, 0 / -1.94 | TG+TG+G+TG+G+ (-4.65) | 63 of 84, +4.58 |
| PVDF | angles relaxed | 180, +/-80 | 0 / -0.63, 0 / -0.71 | TG+ (-3.39) | 78, +3.39 |
| VDCN | rigid | 180, +/-30 | **-8.66** / -13.56, **-8.74** / -13.95 | **TT (-11.48)** | **1 of 83, +0.00** |
| VDCN | angles relaxed | 180, +/-95 | -0.22 / -1.02, -0.29 / -1.25 | **TG+ (-6.39)** | **74, +5.48** |
| AN | rigid | 180, +/-75 | -8.20 / -9.99, -11.90 | TG- (-13.05) | 105 of 164, +6.14 |
| AN | angles relaxed | 180, +/-80 | -1.46 / -3.11, -4.52 | TG- (-7.96) | 131, +6.73 |

PVDF and AN move quantitatively and keep their verdict: a TG-type helix ground state with
all-trans several kcal/mol per monomer above it. VDCN flips. Its rigid `e1(T)` of -8.7 is the
+/-120 deg well of section 1; relaxed, `e1(T)` is -0.2 and the all-trans chain is the 74th
conformation, 5.5 kcal/mol per monomer above a TG+ helix. The screen's one favourable nitrile
result — VDCN's polar zigzag being its own conformational ground state — was the frozen-angle
scan, not the chain. (The relaxed fit's gauche angle for VDCN, +/-95 deg, is itself a
per-*type* result: one angle per backbone atom of the repeat, broadcast along a twelve-monomer
oligomer, cannot accommodate one rotated bond locally the way the per-*atom* relaxation of
section 1 does, where the gauche well sits at +/-40 to 50 deg. The ranking's sensitivity to
that choice is the reservation in the verdict, and it is another reason not to read the
rigid numbers as a ranking.)

## What this changes, and what it does not

* **`docs/SCREEN.md`, conformational table, VDCN row**: "ground state TT, all-trans rank 1 of
  83, +0.00" is withdrawn as a result. It is the rigid scan's mislabelled +/-120 deg well.
  The lattice-level rejection of the nitriles (their polar zigzag 4.1 and 16.3 kcal/mol per
  monomer above their own lattice ground state) is not touched by this: it was measured on
  packed cells and the refinement there does relax the angles.
* **AN's row stands qualitatively** (all-trans +6.1 rigid, +6.7 relaxed above a TG- helix),
  but its "T" state is a basin edge at both levels and the fit's numbers should not be read
  as a ranking either.
* **No fourth state.** Adding one would be fitting a parameter to a frozen-angle artefact —
  the trap the project has fallen into before.
* **The continuous freedom to add is the backbone angle inside `fit_ris`**, and the state
  angles should be adapted from that relaxed scan. It costs 200-250 s per chemistry against 4 s
  rigid at the screen's settings, with the per-type angle broadcast; a per-atom relaxation
  would be the more faithful and the more expensive choice.
* **What stays open**: Sarco's distributed kinks are not single-bond states here at any
  distance from the nitrile, so this potential cannot say whether they are cooperative
  structures, packing accommodations, or MLIP-tier chemistry. Their finite-chain ladder
  failure for the nitriles is not addressed by anything above; it concerns the response, not
  the conformation.
