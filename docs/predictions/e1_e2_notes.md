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
- **E2, PVDF.** On this potential TGTG' is the ground conformation, T3GT3G' lies 0.7 and
  TTTT 1.3 kcal/mol per monomer above it, and the TGTG' to TTTT conformational strain along
  the chain is +10.4%. Published PBE-D3 puts alpha below beta by a few tenths of a kcal/mol
  per monomer, so the sign is expected to hold and the magnitude is expected to fail the 0.1
  kcal/mol threshold: the fitted potential over-binds the helical conformations relative to
  all-trans by about a kcal/mol per monomer, which E2's PVDF row will measure and the
  refit will use.
- **E2, the 70/30 VDF/TrFE model.** The same ordering, TGTG' below T3GT3G' below TTTT, with the
  all-trans penalty lowered from 1.30 to 1.09 kcal/mol per monomer and the helical penalty from
  0.73 to 0.56 by the TrFE units; the conformational strain is unchanged at +10.1%. The direction
  is the known one (TrFE stabilises the all-trans phase of the copolymer); the magnitude, 0.2
  kcal/mol per monomer at 30% TrFE, is what E2's VDF/TrFE row tests. The cells are large (20
  and 40 backbone bonds per chain) and their packing searches are the least converged of the set.
- **E2, the defect units.** The CFE, CTFE, CDFE, FA and DB cells in a beta host are not
  predicted here: the potential's chlorine and the other defect parameters are unvalidated,
  and their validation is one of the things E2 is for. They will be added to this file, as
  predictions, before Dipole's corresponding records are read.
