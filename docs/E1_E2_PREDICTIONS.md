# polyfind's pre-registered predictions for E1 and E2

*2026-11-04. Written and committed before any record of Dipole's E1 or E2 campaigns exists, so
that the comparison is blind on both sides. The PI dispatched E1 and E2 the same day ("go on
E1, E2"). The ask with the specifications and acceptance rules is
`sarco/materials/ask-polyfind-dipole-electrostriction-verification-v1.md`; this file is what
polyfind commits to, the cross-check script applies the rules, and `polyfind/docs/REFERENCE_DATA_REQUEST.md`
stays the exchange log.*

## What is being tested, in one paragraph

polyfind computes crystal-level quantities on a fitted polarizable potential in seconds: the
stiffness, the dielectric tensors, the piezoelectric constants, the electrostriction and the
margins between chain conformations. The plan (`docs/ELECTROSTRICTION_RESEARCH.md`) wants to
use those in place of DFT on the switching screen. E1 and E2 are the first chemistries and
quantities on which Dipole's periodic PBE-D3(BJ) and polyfind are compared quantity by
quantity with the thresholds fixed in advance: E1, sign on every compared component and
magnitude within a factor 1.5 on the polar-polar components; E2, margins within 0.1 kcal/mol
per monomer, the ordering of conformations preserved, the conformational strain within 10%.
A prediction that fails is a result, and says which quantities polyfind may stand in for and
which it may not.

## Definitions, so the two sides compute the same thing

- **Frame.** polyfind: x polar, y lateral, z chain; Voigt 1 = xx, 2 = yy, 6 = xy. Dipole's
  periodic reference: x = a (the lateral axis, 8.36 A), y = b (the polar axis, 4.73 A), z =
  chain. The cross-check swaps x and y; nothing else is rotated.
- **Potential.** `pvdf-dft-valence` (fitted to periodic PBE-D3(BJ) on PVDF), induced dipoles,
  Ewald, force-shifted Lennard-Jones, zero kelvin, the chain rigid along its axis. The chain
  axis strain (Voigt 3) is therefore unreachable here and E1's chain-axis column has no
  polyfind prediction; the ask still wants Dipole's.
- **Stiffness.** `C_IJ = d sigma_J / d eta_K`, relaxed-ion (internal coordinates relaxed at
  each strain), central differences at +/-0.2%, GPa. The reference cells are stress-free in
  the plane but not along the chain: with the chain rigid the axial stress at the reference
  is whatever the chain's own geometry gives (-3.8 GPa for beta, -3.3 for the TGTG' cell,
  +0.2 for PE), while Dipole's cells are fully relaxed. The in-plane quantities are compared
  as they are; the axial stress is reported so the difference is on record.
- **Dielectric tensors.** `eps_inf` is the clamped-ion (electronic) tensor, `eps_0` the
  relaxed-ion one at fixed cell. polyfind's ionic part is chain reorientation and the line
  group's shape parameters, not the optical phonons: expect polyfind's `eps_0` to fall below
  DFPT's where a stretching or bending mode carries dielectric weight (the polar axis of beta).
- **Strained dielectric tensor.** The *true* relaxed-ion `eps_ij(eta)` of the strained cell
  (dipole over the strained volume against the true field `F^-T Et`), at +/-1% normal strain
  along the polar and lateral axes, and its central-difference derivative `d eps_ii / d eta_J`.
  This is the quantity DFPT at a strained, ion-relaxed cell reports.
- **Piezoelectric constants.** `e_iJ` is Vanderbilt's proper constant `dPt_i / d eta_J`
  with `Pt = F^-1 mu / V0`. The true polarization at each strained cell is also given, so the
  improper `dP / d eta` of a Berry-phase slope can be compared with the right one of the two.
- **Electrostriction.** `M_iJ` (field along reference axis i, Voigt strain J) by two routes
  on the same enthalpy, the even part of the zero-stress strain at +/-50 MV/m and the mixed
  stress-field derivative of the nominal polarization; `Q = M / (eps0 chi_free)^2`.
- **E2 energies.** Fully relaxed cells (cell and internal coordinates) of the three
  conformations on the same potential; only differences per monomer carry meaning. The
  monomer is the vinyl unit (two backbone carbons). The 70/30 VDF/TrFE model chain is the
  isotactic one (the model cannot express syndiotactic or atactic placements), with the TrFE
  units at fixed positions in a 10-unit repeat.

## How to deliver, so the script ingests it

A JSON per system in Dipole's own frame, per the schema in `examples/crosscheck_e1_e2.py`:
`C_GPa.block` (3 x 3 normal block, x y z), `eps_inf.diag`, `eps_relaxed_ion.diag` at zero
strain, and `strained.{x|y|z}.{+0.010|-0.010}.{eps_diag, P_C_m2}`; for E2,
`conformations.{name}.{energy_per_monomer_kcal, P_C_m2, c_per_monomer_A}`. Each entry may
carry a `status` string, which the script quotes beside the comparison. Records in Dipole's
own format are fine too; the seat will transcribe them into this schema with provenance.

## Flags raised by the first cross-check, on what Dipole already has for beta-PVDF

The periodic reference already holds three beta-PVDF numbers that E1 asks for, from the
2026-09 exchange. Comparing them (the output is reproduced at the end of this file):

1. **Lateral stiffness.** Dipole's preliminary relaxed-ion stiffness along `a` (the lateral
   axis) is 112 GPa; polyfind's is 18.5. Along the polar axis the two agree within 10% (25.2
   against 22.6 GPa). A stiffness of 112 GPa along a direction held by van der Waals contact
   between chains is far above anything published for beta-PVDF (DFT and experiment put the
   two transverse moduli at 10-30 GPa), so the first thing to settle in E1 is whether the
   `a`-axis relaxation of the 13-cell manifold let the chains slide and rotate, or whether
   polyfind's Lennard-Jones contact is too soft. The reciprocity mismatch Dipole reports
   (1.5 GPa) does not cover a factor of six.
2. **Berry slope.** Dipole's polarization along the polar axis per unit lateral strain is
   0.80 C/m^2 on the 4x8x4 mesh and 0.56 on the dense one, and is not published as a
   piezoelectric constant because of that spread and a branch discontinuity. polyfind's
   true `dP/d eta` for the same component is 0.18 C/m^2 and its proper `e_12` is 0.04. If
   the converged Berry slope stays near 0.6-0.8, the lateral piezoelectric response of beta
   is three to four times what polyfind's rigid chain gives, which would point at the
   internal-strain (chain rotation under lateral strain) part of the response.
3. **Electronic dielectric tensor.** Agreement within 5% on all three axes (polar 2.20
   against 2.24, lateral 2.13 against 2.25, chain 2.51 against 2.60). polyfind's induced
   dipoles were fitted with beta's `eps_inf` in view, so this is a consistency check rather
   than an independent test; PE and the alpha cell are the independent ones.

## Dipole's first E1 record (2026-11-04 evening), and what it already tests

Dipole picked the ask up the same afternoon on the sarco branch `polyfind-electrostriction-reference`
and, rather than reuse the 2026-09 periodic reference, re-relaxed polyfind's own polar beta cell
(the polar-margin work's `pvdf_polar_pbe_d3bj.xyz`) with QE 7.6, PBE-D3(BJ), SSSP precision, 90/360
Ry, an 8x4x16 mesh **in polyfind's axis order** (x polar, y lateral, z chain), cell and ions free, to
forces of 1.1e-5 Ry/bohr and stresses of 0.05 kbar. Its records therefore need no axis swap; the
cross-check reads a `frame: polyfind` flag. Six +/-1% strained cells (ions free, cell fixed) are
declared from that accepted cell and three are running; the alpha (TGTG') cell is relaxing; DFPT with
the D3 Hessian and ionic-mode reduction, and branch-tracked Berry polarization, are the declared next
stages. Dipole marks the cell `quantitatively_valid: false` in its own scope (one local equilibrium,
phase identity and response not yet established).

The one quantity it already fixes is the stress-free cell, and the comparison is clean because both
sides start from the same geometry and the same Hamiltonian family:

| axis | polyfind (fitted potential) | Dipole QE 7.6 | polyfind / Dipole |
|---|---|---|---|
| a, polar (A) | 4.520 | 4.757 | -5.0% |
| b, lateral (A) | 8.521 | 8.528 | -0.1% |
| c, chain (A) | 2.600 | 2.581 | +0.7% |
| volume (A^3) | 100.1 | 104.7 | -4.4% |

polyfind's lateral axis matches QE to 0.1% and its chain axis is the rigid chain's own; the whole
4.4% volume deficit sits on the polar axis, where the fitted potential packs the chains 5% too
close along the dipoles. Two consequences for the quantities to come: polyfind's polar-axis
stiffness and permittivity are evaluated at a smaller spacing than Dipole's, which biases the
polar-polar components (`C_11`, `eps_xx`, `d eps_xx / d eta_1`) in a direction the strained cells will
show, and the polar-margin fit that set this potential was made on energies, not on the polar
spacing. Dipole's two DFT codes also differ between themselves: the 2026-09 CP2K cell had the lateral
axis at 8.358 A (2% below QE's 8.528) and the polar axis at 4.731 (0.5% below), so the lateral
stiffness flag of the previous section is read against the QE cell from here on.

## What the predictions themselves say, before any comparison

- **beta-PVDF.** Stretching the polar axis lowers the polar permittivity (`d eps_xx / d eta_1
  = -2.1`); opening the lateral spacing raises the lateral permittivity strongly (`d eps_yy
  / d eta_2 = +9.2`), because the chain libration softens, and lowers the other two. The
  polar electrostriction is `Q_x1 = +1.6 m^4/C^2`, the lateral `Q_y2 = +10.9`.
- **PE (the control).** No polarization and no piezoelectric response to numerical precision;
  every `d eps / d eta` is negative (a density effect, -2 to -4), and the lattice
  electrostriction is `Q = -3.5 to -3.9 m^4/C^2`. If Dipole finds a positive `d eps / d eta`
  or a `Q` far from these, the strain derivative of the polarizability is wrong on one side.
- **The TGTG' cell.** polyfind's two-chain TGTG' cell relaxes to a setting with a net
  polarization of about 0.07 C/m^2 along x, i.e. a delta-like (polar) arrangement rather
  than the antipolar alpha. Dipole should report the nominal polarization of each relaxed
  TGTG' cell so that like is compared with like; the alpha/delta margin itself is a separate,
  known question (`docs/POLAR_MARGIN_DFT.md` is the beta polar/antipolar analogue).
- **The 70/30 VDF/TrFE model cell (E1).** The isotactic TrFE placement makes the cell monoclinic
  (gamma 92.0 degrees) and gives the polarization a lateral component: `P0 = (0.100, 0.044, 0)`
  C/m^2, the polar component being beta's 0.140 diluted by the VDF content (0.70 x 0.140 = 0.098).
  The TrFE units soften the in-plane stiffness (`C_11` 16.5 against beta's 22.6 GPa, `C_66` 1.7
  against 4.3) and stiffen the chain libration: the lateral ionic contribution to `eps_0` is 0.36
  against beta's 0.70, and `d eps_yy / d eta_2` is +2.1 against beta's +9.2. The proper `e_11` is
  -1.65 times `P0` in both crystals (-0.164 against +0.233 C/m^2 with the opposite sign of `P0`),
  which is the rigid-dipole dilation term. Its electrostriction routes were not run (the cell has
  60 atoms per chain and the routes' tight relaxations dominated a two-hour job); `M` and `Q` for
  this cell are not predicted here.
- **E2, PVDF.** On this potential in the induced-dipole Ewald setting of these runs, TGTG' is the
  ground conformation, T3GT3G' lies 0.7 and TTTT 1.3 kcal/mol per monomer above it, and the
  TGTG' to TTTT conformational strain along the chain is +10.4%. Two cautions, found on
  2026-11-05 while re-reading the record. The TGTG' cell here relaxed to a polar (delta-like)
  packing with 0.07 C/m^2, not the antipolar alpha. And the setting matters: in the
  fixed-charge setting the potential was validated in (`docs/POLAR_MARGIN_DFT.md`), the
  antipolar TGTG' cell lies 0.747 kcal/mol per monomer below beta against Dipole's PBE-D3(BJ)
  -0.726, a match to 0.02; the extra 0.55 of the induced-dipole setting is the induced dipoles'
  own contribution to the two packings, which has never been checked against DFT. E2's PVDF row
  will therefore be compared in both settings, with the polar and antipolar TGTG' packings
  reported separately; the number that stands as the prediction of the validated setting is
  -0.747 for alpha (antipolar) against beta.
- **E2, the 70/30 VDF/TrFE model.** The same ordering, TGTG' below T3GT3G' below TTTT, with the
  all-trans penalty lowered from 1.30 to 1.09 kcal/mol per monomer and the helical penalty from
  0.73 to 0.56 by the TrFE units; the conformational strain is unchanged at +10.1%. The direction
  is the known one (TrFE stabilises the all-trans phase of the copolymer); the magnitude, 0.2
  kcal/mol per monomer at 30% TrFE, is what E2's VDF/TrFE row tests. The cells are large (20
  and 40 backbone bonds per chain) and their packing searches are the least converged of the set.
- **E2, the defect units.** The CFE, CTFE, CDFE and `-CH=CF-` cells in a beta host (the last
  is the unit both FA and DB name: a CFE that has lost HCl, per the Science 2022 main text read
  2026-11-05; the earlier CF3 reading of FA is withdrawn) are not predicted here: the potential's chlorine and the other defect parameters are unvalidated,
  and their validation is one of the things E2 is for. They will be added to this file, as
  predictions, before Dipole's corresponding records are read.


## E1 predictions

Potential: pvdf-dft-valence, induced dipoles, Ewald, rigid chain, force-shifted LJ, 0 K. Frame: x polar, y lateral, z chain; Voigt 1 xx, 2 yy, 6 xy; Dipole's frame has x and y swapped. Strain step +0.010 on each side.


### E1, cell and stiffness

| system | a (A, polar) | b (A, lateral) | gamma | c (A, chain) | V (A^3) | C_11 | C_12 | C_22 | C_66 (GPa) |
|---|---|---|---|---|---|---|---|---|---|
| beta-PVDF (TTTT) | 4.520 | 8.521 | 90.00 | 2.5999 | 100.12 | 22.57 | 3.09 | 18.55 | 4.34 |
| PE (control) | 4.061 | 7.088 | 90.00 | 2.5511 | 73.43 | 23.28 | 1.56 | 25.04 | 2.41 |
| alpha-PVDF (TGTG') | 4.977 | 9.038 | 90.00 | 4.7102 | 211.88 | 16.57 | 5.01 | 12.26 | 5.97 |
| VDF/TrFE 70/30 model (TTTT) | 4.793 | 8.852 | 91.96 | 25.9918 | 1102.19 | 16.51 | 1.20 | 16.96 | 1.74 |

### E1, dielectric tensors at the stress-free cell, proper piezoelectric constants, polarization

| system | eps_inf (x, y, z) | eps_0 relaxed-ion (x, y, z) | e_11, e_12, e_26 (C/m^2) | d_11, d_12 (pm/V) | P0 (C/m^2) |
|---|---|---|---|---|---|
| beta-PVDF (TTTT) | 2.1970, 2.1330, 2.5145 | 2.1970, 2.8311, 2.5145 | 0.2325, 0.0367, 0.0967 | 10.27, 0.27 | -0.1400, 6.38e-17, -4.57e-17 |
| PE (control) | 2.8782, 2.8613, 2.9821 | 2.8782, 2.8613, 2.9821 | 3.48e-14, 1.89e-14, 3.58e-14 | -1.35e-13, -2.00e-13 | -9.58e-17, -1.33e-16, -6.52e-18 |
| alpha-PVDF (TGTG') | 2.0749, 2.1084, 2.3804 | 2.0749, 2.2521, 2.3804 | -0.1370, -0.0079, -0.0776 | -9.21, 3.15 | 0.0726, -2.04e-17, -1.52e-17 |
| VDF/TrFE 70/30 model (TTTT) | 2.0430, 2.0078, 2.3794 | 2.0696, 2.3709, 2.3794 | -0.1637, -0.0436, -0.0188 | -10.23, 0.34 | 0.0999, 0.0435, -5.83e-08 |

### E1, the true relaxed-ion dielectric tensor against normal strain (central differences at +/-1%)

| system | strain axis | d eps_xx / d eta | d eps_yy / d eta | d eps_zz / d eta | dP_x / d eta (true, C/m^2) | eps diag at +1% | eps diag at -1% |
|---|---|---|---|---|---|---|---|
| beta-PVDF (TTTT) | xx (polar) | -2.055 | -0.029 | -1.565 | 0.2326 | 2.1769, 2.8314, 2.4991 | 2.2180, 2.8320, 2.5304 |
| beta-PVDF (TTTT) | yy (lateral) | -1.354 | 9.246 | -1.537 | 0.1767 | 2.1837, 2.9342, 2.4993 | 2.2108, 2.7492, 2.5301 |
| PE (control) | xx (polar) | -3.272 | -1.917 | -2.139 | -8.67e-16 | 2.8461, 2.8420, 2.9610 | 2.9116, 2.8804, 3.0038 |
| PE (control) | yy (lateral) | -2.337 | -4.095 | -2.059 | 1.17e-15 | 2.8553, 2.8209, 2.9617 | 2.9020, 2.9028, 3.0029 |
| alpha-PVDF (TGTG') | xx (polar) | -1.860 | 0.629 | -1.383 | -0.1370 | 2.0566, 2.2599, 2.3667 | 2.0938, 2.2473, 2.3943 |
| alpha-PVDF (TGTG') | yy (lateral) | -1.314 | -0.841 | -1.373 | -0.0806 | 2.0620, 2.2441, 2.3668 | 2.0882, 2.2609, 2.3942 |
| VDF/TrFE 70/30 model (TTTT) | xx (polar) | -1.585 | 0.239 | -1.408 | -0.1637 | 2.0541, 2.3737, 2.3655 | 2.0858, 2.3689, 2.3936 |
| VDF/TrFE 70/30 model (TTTT) | yy (lateral) | -0.727 | 2.061 | -1.396 | -0.1435 | 2.0626, 2.3915, 2.3656 | 2.0772, 2.3503, 2.3935 |

### E1, electrostriction (field along reference axis i, Voigt strain J), both routes

| system | field | M_strain (1e-18 m^2/V^2), J = 1, 2 | M_polarization, J = 1, 2 | chi_free | Q (m^4/C^2), J = 1, 2 |
|---|---|---|---|---|---|
| beta-PVDF (TTTT) | x (polar) | 0.00027, 0.00002 | 0.00027, 0.00002 | 1.468 | 1.59, 0.10 |
|  | y (lateral) | -0.00042, 0.00368 | -0.00043, 0.00368 | 2.075 | -1.26, 10.92 |
|  | z (chain) | -9.32e-06, -3.83e-06 | -9.32e-06, -3.83e-06 | 1.515 | -0.05, -0.02 |
| PE (control) | x (polar) | -0.00098, -0.00002 | -0.00098, -0.00002 | 1.878 | -3.54, -0.07 |
|  | y (lateral) | 0.00006, -0.00106 | 0.00006, -0.00106 | 1.861 | 0.22, -3.89 |
|  | z (chain) | -0.00003, -0.00001 | -0.00003, -0.00001 | 1.982 | -0.09, -0.04 |
| alpha-PVDF (TGTG') | x (polar) | -0.00015, 0.00007 | -0.00015, 0.00007 | 1.215 | -1.32, 0.58 |
|  | y (lateral) | 0.00145, -0.00134 | 0.00145, -0.00134 | 1.366 | 9.89, -9.14 |
|  | z (chain) | -1.74e-06, 3.40e-06 | -1.74e-06, 3.36e-06 | 1.380 | -0.01, 0.02 |
| VDF/TrFE 70/30 model (TTTT) | x (polar) | 0.00046, -2.56e-06 | 0.00046, -2.01e-06 | 1.382 | 3.06, -0.02 |
|  | y (lateral) | 0.00037, 0.00025 | 0.00037, 0.00025 | 1.404 | 2.41, 1.62 |
|  | z (chain) | -7.34e-06, -3.21e-06 | -7.33e-06, -3.21e-06 | 1.379 | -0.05, -0.02 |

## E2 predictions

Potential: pvdf-dft-valence, induced dipoles, Ewald, rigid chain, force-shifted LJ, 0 K. Energies per monomer on the same potential, so only differences carry meaning.


### E2, conformations

| system | conformation | E per monomer (kcal/mol) | margin vs ground | |P| (C/m^2) | c per monomer (A) | conformational strain vs ground | V per monomer (A^3) | cell a, b, gamma, c |
|---|---|---|---|---|---|---|---|---|
| PVDF | TTTT | -4.7436 | 1.304 | 0.1400 | 2.5999 | 10.4% | 50.06 | 4.520, 8.521, 90.0, 2.600 |
| PVDF | TGTG' | -6.0473 | 0.000 | 0.0726 | 2.3551 | 0.0% | 52.97 | 4.977, 9.038, 90.0, 4.710 |
| PVDF | T3GT3G' | -5.3225 | 0.725 | 0.0848 | 2.3356 | -0.8% | 53.94 | 5.287, 8.762, 85.6, 9.342 |
| VDF/TrFE 70/30 model | TTTT | -3.1461 | 1.087 | 0.1090 | 2.5992 | 10.1% | 55.11 | 4.793, 8.852, 92.0, 25.992 |
| VDF/TrFE 70/30 model | TGTG' | -4.2332 | 0.000 | 0.0642 | 2.3601 | 0.0% | 59.71 | 5.341, 9.476, 91.2, 23.601 |
| VDF/TrFE 70/30 model | T3GT3G' | -3.6733 | 0.560 | 0.0615 | 2.3418 | -0.8% | 59.55 | 5.565, 9.156, 93.6, 46.836 |

## First cross-check against what Dipole already has

Output of `examples/crosscheck_e1_e2.py` on `docs/predictions/e1_dipole_existing.json` (Dipole's statuses quoted inline):

```
E1: polyfind pvdf-dft-valence, induced dipoles, Ewald, rigid chain, force-shifted LJ, 0 K; Dipole frame: Dipole periodic reference: x = a (lateral), y = b (polar), z = chain

=== E1 beta-PVDF (TTTT)  <->  beta-PVDF, CP2K PBE-D3(BJ) zero-pressure cell (STATUS.md 2026-09-10) ===
  cell (polyfind a b c / Dipole a b c, A): 4.520 8.521 2.600 / 8.358 4.731 2.580  (polyfind's a is along x polar; Dipole's a is the lateral axis)
  C (GPa), Dipole status: preliminary relaxed-ion normal stiffness from the 13-cell strain manifold; reciprocity mismatch up to 1.5 GPa
  C_11 (polar,polar)                   polyfind      22.5655   Dipole      25.1850   ratio 0.90: within 1.5
  C_12 (polar,lateral)                 polyfind       3.0913   Dipole      10.0540   sign agrees
  C_21 (lateral,polar)                 polyfind       3.0913   Dipole       8.5470   sign agrees
  C_22 (lateral,lateral)               polyfind      18.5496   Dipole     111.9340   sign agrees
  C_33 (chain): polyfind rigid chain, none; Dipole 315.9
  eps_inf, Dipole status: QE DFPT electronic tensor at the accepted cell (BORN_RESPONSE_REPORT.md): complete evidence for its protocol, not an accepted quantitative reference
  eps_inf_x polar                      polyfind       2.1970   Dipole       2.2354   ratio 0.98: within 1.5
  eps_inf_y lateral                    polyfind       2.1330   Dipole       2.2526   sign agrees
  eps_inf_z chain                      polyfind       2.5145   Dipole       2.6005   sign agrees
  Berry slope dP_polar/d eta_lateral (Dipole, polarization along y per unit strain along x, centered 1%; 30% apart between k-grids and branch-discontinuous at the dense zero point: not published as a piezoelectric constant): grid_4x8x4 +0.803, grid_4x8x16 +0.563
    polyfind dP_polar/d eta_lateral (true) +0.177 C/m^2; proper e_12 +0.037

=== E1 beta-PVDF (TTTT)  <->  beta-PVDF, Dipole QE 7.6 PBE-D3(BJ) equilibrium reference pvdf_beta_ref_v1 (branch polyfind-electrostriction-reference, 2026-11-04) ===
  Dipole status: terminal complete 2026-11-04 23:06 UTC; the six +/-1% axial cells are complete (record beta_qe_axial_v1); no DFPT, polarization or Dipole-side stiffness yet
  cell (A), polyfind / Dipole, same axis order (x polar, y lateral, z chain):
  a (polar)                            polyfind       4.5196   Dipole       4.7568   -5.0%
  b (lateral)                          polyfind       8.5208   Dipole       8.5282   -0.1%
  c (chain)                            polyfind       2.5999   Dipole       2.5807   +0.7%
  volume (A^3)                         polyfind      100.125   Dipole      104.693   -4.4%
  Dipole total energy -137.903515 Eh per 2 monomers (differences only are comparable; E2)

=== E1 beta-PVDF (TTTT)  <->  beta-PVDF, Dipole QE 7.6 PBE-D3(BJ): six +/-1% axial strain cells from pvdf_beta_ref_v1 (exchange main 54897c21, 2026-11-05) ===
  Dipole status: six axial cells complete on exchange main at 54897c21 (2026-11-05); Dipole's lattice Hessian and Gamma response stages declared, not run
  cell (A), polyfind / Dipole, same axis order (x polar, y lateral, z chain):
  a (polar)                            polyfind       4.5196   Dipole       4.7568   -5.0%
  b (lateral)                          polyfind       8.5208   Dipole       8.5282   -0.1%
  c (chain)                            polyfind       2.5999   Dipole       2.5807   +0.7%
  volume (A^3)                         polyfind      100.125   Dipole      104.693   -4.4%
  C (GPa), Dipole status: derived by polyfind from Dipole's six strained-cell stresses and energies (relaxed ions, fixed strained cell); Dipole's own stiffness/DFPT stage is pending
  C_11 (polar,polar)                   polyfind      22.5655   Dipole      22.7690   ratio 0.99: within 1.5
  C_12 (polar,lateral)                 polyfind       3.0913   Dipole       5.8820   sign agrees
  C_21 (lateral,polar)                 polyfind       3.0913   Dipole       5.8820   sign agrees
  C_22 (lateral,lateral)               polyfind      18.5496   Dipole      24.5190   sign agrees
  C_33 (chain): polyfind rigid chain, none; Dipole 302.6
```


*Rendered by `examples/e1_e2_report.py` from `docs/predictions/e1_polyfind.json`, `docs/predictions/e2_polyfind.json` and `docs/predictions/e1_dipole_existing.json`; re-run it when a file changes.*
