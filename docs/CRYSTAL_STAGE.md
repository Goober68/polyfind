# The crystal stage of the fit: nonbonded terms against the DFT crystal anchors

*2026-11-05, kept current through 2026-11-07. Code: `polyfind.crystalfit` (bounded Powell and,
since 5.18, bounded least squares), driver `examples/crystal_stage.py`, results
`docs/predictions/crystal_stage_v*.json`. Sections 1 to 4 are the first stage and its preset
`pvdf-dft-valence-crystal`; section 5 is the stages after it in order, 5.1 to 5.22, and the
current preset is `pvdf-dft-valence-belt4` (5.22: the belt3 problem with the alpha-beta gap row
evaluated with the pendant groups relaxed, which the rigid refinement cannot do and DFT's alpha
does; belt3 is 5.18, 5.19). Every preset is opt-in; the E1/E2
predictions pre-registered on `pvdf-dft-valence` stand unchanged. `RESUME.md` opens with the
one-paragraph state.*

## 1. Why a second stage

`docs/VALENCE_FIT.md` fitted the torsions, bonds, angles, Lennard-Jones multipliers and
bond-charge increments of `pvdf-dft-valence` to 467 frames of PBE-D3 conformers. That data
cannot see a crystal. Relative conformer energies barely see the overall dispersion, so the
Lennard-Jones terms were bounded rather than determined (three of the four well-depth
multipliers sit at their 0.5 bound, section 5 of that document); and energies and forces
cannot separate a charge scale from a dielectric constant at all, which only a polarization
can. Measured against Dipole's periodic PBE-D3(BJ) cells on 2026-11-05 (the first E1 record,
`docs/E1_E2_PREDICTIONS.md`), the potential reproduces the alpha-beta margin to 0.02 kcal/mol
per monomer and the lateral and chain axes of beta to 0.1% and 0.7%, but packs beta's polar
axis 4.5% too close and alpha's long axis 6% too close, and carries a fixed-charge
polarization of 0.115 C/m^2 against the DFT 0.18.

Those three misses are all nonbonded: the polar axis is the distance between dipole-bearing
chains, the alpha long axis the same contact in the other packing, and the polarization the
charge the increments put on the C-F and C-H bonds. The mission (`RESUME.md`: a solver an
order of magnitude faster than DFT at DFT-like accuracy) needs them closed without undoing
what the conformer data did determine, so this stage moves only the terms the trainset left
free, holds them near the trainset values with a ridge, and leaves every torsion, bond and
angle term where the valence fit put it.

## 2. The anchors

All periodic PBE-D3(BJ), in polyfind's axis order (x polar, y lateral, z chain):

| anchor | value | source |
|---|---|---|
| beta cell | 4.757, 8.528, 2.581 A | Dipole's QE 7.6 `pvdf_beta_ref_v1` (2026-11-04): SSSP 90/360 Ry, 8x4x16 mesh, forces 1e-5 Ry/bohr, stress 0.05 kbar |
| alpha cell (antipolar TGTG' packing) | 4.945, 9.636, 4.660 A | the QE cell of `docs/POLAR_MARGIN_DFT.md` |
| E(alpha) - E(beta) | -0.726 kcal/mol per monomer | two-body D3, the shared anchor of 2026-10-09 |
| E(antipolar beta) - E(polar beta) | -0.026 kcal/mol per monomer | `docs/POLAR_MARGIN_DFT.md` (-1.15 meV) |
| beta polarization | 0.18 C/m^2 | published Berry-phase values, 0.178-0.181 |

Tolerances, which set the residual scale (`z = (polyfind - DFT) / tol`): 1% of each cell
axis, 0.01 C/m^2 on the polarization, 0.05 kcal/mol per monomer on the two margins. The
objective is the sum of squared `z` plus the ridge.

## 3. What moves, and how

Eight variables on top of `pvdf-dft-valence`: multipliers on the Lennard-Jones diameters of
C, H and F (bounds 0.90-1.10), multipliers on their well depths (0.5-2.0) and shifts of the
C-H and C-F bond-charge increments (+/-0.08 e). The ridge is Gaussian with widths 0.05 on
the diameters, 0.35 in the logarithm of the well depths and 0.05 e on the increments, weight
1, so a variable moved by one width costs one unit of objective, the same as a one-tolerance
miss on one anchor.

Every evaluation packs and refines polar beta (with its antipolar competitor, for the gap)
and the antipolar alpha cell exactly as the crystal fit of `docs/SCREEN.md` does
(`polyfind.fitting.predict`), in the fixed-charge Ewald setting the margins were validated
in; about 9 s on the 4-core box. Powell with bounds from the valence values; the fit is a
few hundred evaluations. The induced-dipole setting is evaluated afterwards, not fitted
(section 4), because the dielectric response the induced dipoles add is itself one of the
quantities E1 is testing.

## 4. Results

301 evaluations of Powell, 44 minutes on the 4-core box (while other jobs ran), objective
104.59 -> 41.54, of which 9.73 is the ridge. Fixed-charge Ewald setting throughout:

| quantity | DFT | `pvdf-dft-valence` | crystal stage | z before -> after |
|---|---|---|---|---|
| beta a (polar), A | 4.757 | 4.544 (-4.5%) | 4.581 (-3.7%) | -4.48 -> -3.70 |
| beta b (lateral), A | 8.528 | 8.545 (+0.2%) | 8.663 (+1.6%) | +0.19 -> +1.59 |
| beta c (chain), A | 2.581 | 2.601 (+0.8%) | 2.613 (+1.2%) | +0.78 -> +1.23 |
| alpha a, A | 4.945 | 4.991 (+0.9%) | 5.015 (+1.4%) | +0.92 -> +1.42 |
| alpha b, A | 9.636 | 9.036 (-6.2%) | 9.372 (-2.7%) | -6.23 -> -2.74 |
| alpha c, A | 4.660 | 4.713 (+1.1%) | 4.751 (+1.9%) | +1.13 -> +1.94 |
| beta polarization, C/m^2 | 0.18 | 0.115 | 0.174 | -6.52 -> -0.56 |
| E(alpha) - E(beta), kcal/mol per monomer | -0.726 | -0.744 | -0.748 | -0.37 -> -0.44 |
| E(antipolar) - E(polar) beta, rigid | -0.026 | 0.000 | 0.000 | +0.52 -> +0.52 |

The variables: diameters x_C 0.948, x_H 1.073, x_F 1.016; well depths D_C 0.580, D_H 1.239,
D_F 0.965; increments C-H -0.119 -> -0.195 e, C-F +0.114 -> +0.171 e. The ridge of 9.7 is a
combined move of about three widths, most of it in the two increments, the carbon well depth
and the hydrogen diameter.

**What closed.** The polarization: both bond charges grow by about half and `P` follows them,
0.115 -> 0.174 against 0.18. The alpha long axis: -6.2% -> -2.7%, from the larger hydrogen and
the smaller, shallower carbon. The two margins stay inside their tolerance, the alpha-beta one
by a wide margin: the torsion terms that carry it were not touched.

**What did not.** Beta's polar axis moves only from -4.5% to -3.7%, and the four axes that
were right drift out by 0.5-1.4% (lateral +1.6%, chain +1.2%, alpha a +1.4%, alpha c +1.9%).
The pattern is structural rather than a matter of more evaluations: the charges are now set
by `P`, and with them the electrostatic attraction between the head-to-tail dipoles of
neighbouring chains along the polar axis, which is what pulls that axis in; the only thing
that can push it back out is the H...F contact between a CH2 and the facing CF2 of the next
chain, while the lateral axis is set by F...F contacts between chains side by side. With one
isotropic Lennard-Jones term per element, the hydrogen and fluorine diameters that open the
polar contact also open the lateral one, and the fit stops where the two misses balance.
The antipolar-polar gap stays at zero because the rigid antipolar branch of `predict` is
degenerate with the polar one on any pairwise potential (``docs/POLAR_MARGIN_DFT.md``); it is
in the objective as a guard, not as a target this stage can reach.

**Verdict for the mission.** On the fixed-charge setting the stage brings the polarization and
both margins inside their tolerances and the six cell axes to within 1.2-3.7%, against
0.2-6.2% before: better on the three quantities that were far off, worse by about one
tolerance on the four that were right. That is not yet DFT-like accuracy on the cell. The next
lever is a form change, not a longer fit: a pair-specific H-F term (or an exponential-6
fluorine), which decouples the polar contact from the lateral one, fitted to the same anchors
plus the six strained-cell energies Dipole is computing (which pin the curvature, not just the
minimum).

### 4.1 The induced-dipole setting, the one E1 and E2 are predicted in

`examples/e1_predictions.py --preset pvdf-dft-valence-crystal`
(`docs/predictions/e1_polyfind_crystal.json`), the same construction as the pre-registered
set (induced dipoles, Ewald, rigid chain, force-shifted Lennard-Jones, 0 K):

| quantity | DFT | `pvdf-dft-valence` | crystal stage |
|---|---|---|---|
| beta a (polar), A | 4.757 | 4.520 (-5.0%) | 4.530 (-4.8%) |
| beta b (lateral), A | 8.528 | 8.521 (-0.1%) | 8.596 (+0.8%) |
| beta c (chain), A | 2.581 | 2.600 (+0.7%) | 2.610 (+1.1%) |
| beta polarization, C/m^2 | 0.18 | 0.140 | 0.214 |
| beta eps_inf (x, y, z) | 2.24, 2.25, 2.60 (Dipole's DFPT at the CP2K cell) | 2.197, 2.133, 2.515 | 2.171, 2.108, 2.495 |
| beta eps_0 relaxed-ion (x, y, z) | - | 2.197, 2.831, 2.515 | 2.171, 3.430, 2.495 |
| beta C_11, C_22, C_12, GPa | 22.8, 24.5, 5.9 (Dipole's six +/-1% QE cells, 2026-11-05) | 22.6, 18.6, 3.1 | 31.5, 23.5, 4.2 |
| beta C_66, GPa | - | 4.3 | 5.7 |
| beta e_11, C/m^2 | - | 0.233 | 0.356 |
| alpha (polar TGTG' packing) a, b, c, A | - | 4.977, 9.038, 4.710 | 4.984, 9.335, 4.744 |
| alpha polarization, C/m^2 | - | 0.073 | 0.104 |
| PE a, b, c, A | - | 4.061, 7.088, 2.551 | 4.202, 7.436, 2.544 |

(The sign of P and of e_11 differs between the two files because the packer chose the other
of the two symmetry-equivalent orientations of the CF2 dipoles along x; magnitudes compare.)

**Reading.** With induced dipoles the polar axis does not open: -5.0% -> -4.8%. The dipoles
the fixed-charge fit did not see add their own attraction along the polar axis and take
back most of the 0.8% the fit gained there, while the lateral and chain axes carry the
fit's drift (+0.8%, +1.1%). The polarization overshoots: the increments were set where the
induced dipoles were absent, and with them present P is 0.214, 19% above DFT, where the
valence preset was 22% below. The ionic part of the lateral dielectric constant nearly
doubles (eps_0 - eps_inf along y: 0.70 -> 1.32), as it must when it scales with the square of
the bond charges. Dipole's six +/-1% axial cells, which landed the same morning
(`docs/E1_E2_PREDICTIONS.md`, record `beta_qe_axial_v1`), give the first DFT stiffness to judge
by: C_11 22.8, C_22 24.5, C_12 5.9 GPa (stress and energy routes within 1.5%). The valence
preset has the polar stiffness C_11 to 1% and the lateral one 24% low; this preset has the
lateral one to 4% and the polar one 38% high; both have C_12 low (3.1 and 4.2). Neither
potential carries the right lateral-to-polar stiffness ratio, which is the same H...F against
F...F contact problem as the cell axes, seen through second derivatives. PE, which is not an anchor, moves 3.5-5% towards the measured cell (its
cross-section per chain 14.4 -> 15.6 A^2; experiment at room temperature 18.4).

In the setting E1 is judged in, the pre-registered set on `pvdf-dft-valence` therefore
remains the better prediction on every quantity but the magnitude of P; both sets stay in
the record and the cross-check runs on both. The lesson is the obvious one: fit in the
setting you predict in. The stage was fitted fixed-charge because the fit's packing and
refinement path (`fitting.predict`, through `pack` and `refine_crystal`) has no induced-dipole
option, while the mechanics path that has one (`refined_reference`) has no antipolar branch
for the alpha anchor; either needs plumbing before an induced-dipole fit, and that fit should
come with the form change of section 4, since the same isotropic terms would hit the same
ceiling on the polar axis.

**The copolymer cell and the E2 margins on this preset.** The 120-atom 70/30 VDF/TrFE
all-trans cell (`e1_polyfind_crystal.json`, 343 s with four workers): a 4.793 -> 4.765, b
8.852 -> 8.932, gamma 91.96 -> 92.68 deg, c per 10 monomers 25.99 -> 26.07 A; |P| 0.109 ->
0.169 C/m^2; eps_0 (x, y, z) 2.070, 2.371, 2.379 -> 2.085, 2.637, 2.374; C_11, C_22, C_66 16.5,
17.0, 1.7 -> 21.7, 20.6, 3.0 GPa; intrinsic Q along the polar axis 3.06 -> 3.92 m^4/C^2 and
M_x1 0.46 -> 0.74 x 1e-21 m^2/V^2 (chi_free 1.38 -> 1.55): the same stiffening, larger
polarization and larger ionic response as beta, with no DFT anchor yet to judge them against
(Dipole's copolymer cell is in E1). E2 (`e2_polyfind_crystal.json`, same script, same
construction; the numbers here are the corrected ones of section 5.4, with the backbone-angle
strain included): PVDF's ground conformation is still TGTG' and the margins open, TTTT +0.85
-> +0.98 and T3GT3G' +0.41 -> +0.52 kcal/mol per monomer against DFT's +0.726 for the alpha-beta
pair -- the fixed-charge fit reproduces that margin to 0.02 and the induced dipoles, on the
larger charges, add polar-lattice stabilisation that the fit never saw. The copolymer's
margins open the same way, TTTT +0.46 -> +0.53 and T3GT3G' +0.11 -> +0.18 kcal/mol per
monomer over TGTG' (conformational strain of the switch +9.9%, unchanged).


## 5. Stage 2: the setting the predictions are made in, a pair term, the stiffness

*2026-11-05, afternoon. `crystalfit.fit(stage=2)`, `examples/crystal_stage.py --stage 2`, result
`docs/predictions/crystal_stage_v2.json`.*

Section 4 ended with two lessons: fit in the setting you predict in, and change the form.
Stage 2 does both and adds the first DFT anchors on the response, Dipole's six strained cells.

**The setting.** Every evaluation now packs and refines beta (with its antipolar competitor)
and the antipolar alpha cell with induced dipoles, Ewald electrostatics and the
force-shifted Lennard-Jones term, through the same `fitting.predict` path as stage 1, which
took `polarizable=` and `lj_cutoff=` for the purpose (`pack` and `refine_crystal` with it).
At the start point this reproduces the pre-registered E1 construction to the digits printed
in `docs/E1_E2_PREDICTIONS.md`: beta 4.517, 8.521, 2.602 A, P 0.140 C/m^2, C_11 22.5, C_22
18.5, C_12 3.1 GPa.

**The form.** A tenth and ninth variable give the H-F pair its own Lennard-Jones diameter and
well depth, as multipliers on the geometric means of the (already multiplied) H and F
values, so (1, 1) is the mixing rule and changes nothing. This is the smallest change that
decouples the polar H...F contact from the lateral F...F one (section 4). It reaches the
lattice kernel's pair tables and the force field's own pair terms through
`FFParameters.lj_pairs` / `SimpleFF.lj_pairs` (bounds 0.85-1.15 on the diameter, 0.3-3 on
the depth; ridge widths 0.05 and 0.5 in the logarithm).

**The anchors.** The nine of stage 1 plus beta's relaxed-ion C_11, C_22 and C_12 (22.77,
24.52, 5.88 GPa, derived from Dipole's six +/-1% cells, section 4.1), tolerance 10% each,
computed at every evaluation from a `Reference` rebuilt from the refined beta prediction
exactly as the E1 script builds its own (`crystalfit.reference_from_prediction`), with the
six strained states relaxed in parallel. One evaluation is 11 s on the 4-core box, so the
fit costs about as much as stage 1 did.

### 5.1 Result

321 evaluations, 40 minutes on the 4-core box (7 s each once the tables are built), objective
112.99 -> 65.05 of which 4.87 is the ridge. Induced dipoles, Ewald, force-shifted Lennard-Jones
throughout, so the start column is the pre-registered E1 construction:

| quantity | DFT | start (`pvdf-dft-valence`) | stage 2 | z before -> after |
|---|---|---|---|---|
| beta a (polar), A | 4.757 | 4.517 (-5.0%) | 4.572 (-3.9%) | -5.04 -> -3.89 |
| beta b (lateral), A | 8.528 | 8.521 (-0.1%) | 8.656 (+1.5%) | -0.09 -> +1.50 |
| beta c (chain), A | 2.581 | 2.602 (+0.8%) | 2.610 (+1.1%) | +0.81 -> +1.12 |
| alpha a, A | 4.945 | 4.970 (+0.5%) | 5.031 (+1.7%) | +0.50 -> +1.73 |
| alpha b, A | 9.636 | 9.069 (-5.9%) | 9.339 (-3.1%) | -5.88 -> -3.09 |
| alpha c, A | 4.660 | 4.706 (+1.0%) | 4.742 (+1.8%) | +0.98 -> +1.77 |
| beta polarization, C/m^2 | 0.18 | 0.140 | 0.189 | -4.00 -> +0.93 |
| E(alpha) - E(beta), kcal/mol per monomer | -0.726 | -0.854 | -0.758 | -2.55 -> -0.64 |
| E(antipolar) - E(polar) beta, rigid | -0.026 | 0.000 | 0.000 | +0.52 -> +0.52 |
| beta C_11 (polar), GPa | 22.77 | 22.53 | 27.30 | -0.10 -> +1.99 |
| beta C_22 (lateral), GPa | 24.52 | 18.52 | 21.03 | -2.45 -> -1.42 |
| beta C_12, GPa | 5.88 | 3.10 | 3.36 | -4.72 -> -4.28 |

The variables: diameters x_C 1.001, x_H 1.019, x_F 1.011; well depths D_C 0.62, D_H 1.20, D_F
0.94; increments C-H -0.119 -> -0.173 e, C-F +0.114 -> +0.157 e; the H-F pair diameter x 1.037
and depth x 0.937 on their geometric means. The margin and the polarization are inside their
tolerances in the setting that matters, which stage 1 could not say. The polar axis moves as far
as stage 1 moved it (-3.9% here, -3.7% there) and no further, the lateral axis drifts out by the
same 1.5%, the polar stiffness overshoots once the charges are right, and C_12 stays at 57% of
DFT. The pair term, free to move the H-F diameter by 15%, moved it 3.7%: it is not where the
remaining error is.

### 5.2 Where the error is: a scan along the polar axis

The energy of the beta cell against `a` at the DFT `b` and `c`, internal coordinates relaxed,
split into its Lennard-Jones part (the same geometry on a charge-free twin) and the rest, which
is the electrostatics with the induced dipoles (kcal/mol per two-monomer cell, slopes at the DFT
`a = 4.757` A, positive meaning the term wants the axis shorter):

| preset | dE/da total | Lennard-Jones | electrostatic |
|---|---|---|---|
| `pvdf-dft-valence` | +1.99 | +0.65 | +1.34 |
| stage 2 | +1.93 | -0.74 | +2.67 |

On both potentials the electrostatics pulls the polar axis in by 1.3-2.7 kcal/mol per A at the
DFT spacing, and the fits respond the only way their variables allow, by inflating the
repulsion: the stage-2 Lennard-Jones term pushes outward at the DFT spacing, and still loses.
The electrostatic attraction is between the H of one chain's CH2 and the F of the facing CF2
across 2.5-3 A, where point charges are known to overbind: real charge distributions overlap
there and the interaction falls below 1/r (charge penetration). DFT has that screening; the
model has none, for the permanent charges (the induced dipoles already carry Thole damping).
Bigger diameters cannot supply it without opening every other contact, which is exactly the
lateral drift of both stages. The next change to the form is therefore a short-range damping
of the permanent-charge interaction, Gaussian-smeared charges in the Ewald real space
(`erf(beta r)/r` in place of `1/r`, one width to fit), not another Lennard-Jones variable.

### 5.3 The stage-2 preset as a prediction set

`examples/e1_predictions.py --preset pvdf-dft-valence-crystal2` (`docs/predictions/e1_polyfind_crystal2.json`,
the copolymer in 290 s with four workers; the E2 file beside it). Against the pre-registered set:

| quantity | DFT | `pvdf-dft-valence` | stage 2 |
|---|---|---|---|
| beta a, b, c, A | 4.757, 8.528, 2.581 | 4.520, 8.521, 2.600 | 4.576, 8.657, 2.606 |
| beta polarization, C/m^2 | 0.18 | 0.140 | 0.189 |
| beta eps_0 relaxed-ion (x, y, z) | - | 2.197, 2.831, 2.515 | 2.145, 3.253, 2.469 |
| beta C_11, C_22, C_12, C_66, GPa | 22.8, 24.5, 5.9, - | 22.6, 18.5, 3.1, 4.3 | 27.4, 21.1, 3.3, 4.7 |
| beta intrinsic Q_x1, m^4/C^2 | - | +1.59 | +3.45 |
| alpha (polar TGTG') a, b, c, A | - | 4.977, 9.038, 4.710 | 5.043, 9.276, 4.750 |
| alpha polarization, C/m^2 | - | 0.073 | 0.095 |
| copolymer a, b, c per 10 monomers, A | - | 4.793, 8.852, 25.99 | 4.794, 8.963, 26.04 |
| copolymer polarization, C/m^2 | - | 0.109 | 0.151 |
| copolymer C_11, C_22, C_66, GPa | - | 16.5, 17.0, 1.7 | 19.8, 18.7, 2.4 |
| copolymer intrinsic Q_x1, m^4/C^2 | - | +3.06 | +3.94 |
| PE (control) a, b, gamma, A, deg | - | 4.061, 7.088, 90.0 | 4.326, 6.865, 87.8 |
| PE C_11, C_22, C_66, GPa | - | 23.3, 25.0, 2.4 | 15.9, 14.6, 10.3 |

The fit path and the E1 path agree on the fitted beta cell to 0.004 A, as they should now that
they share one construction. PE, which is not an anchor, is the warning in this table: its cell
goes monoclinic and its stiffness halves, so the hydrogen and carbon terms the PVDF anchors
pulled on are not transferable to a chain without fluorine. A PE anchor (Dipole's PE control is
in the E1 request) belongs in the next objective.

**The chain geometry itself.** Dipole's reference gives the beta chain C-F 1.386 A, C-H 1.096,
C-C 1.533, F-C-F 105.1 deg, H-C-H 107.6 and a backbone angle of 114.6; polyfind builds its rigid
chain at the textbook C-F 1.350, C-H 1.090, C-C 1.528, 106.0, 108.0 and 114.0, and the lattice
then opens the backbone angle to 116.6 (c +0.7%). Rebuilding the chain at the PBE-D3 geometry
on the valence preset moves beta's polar axis from -5.0% to -4.4% and the lateral and chain
axes to +0.5% and +1.3%, with P and the stiffness unchanged: a real but minor part of the polar
deficit, worth adopting as the consistent choice against PBE-D3 anchors, not the cause.

### 5.4 E2, and a correction to the pre-registered margins

Comparing the E2 margins of the stage-2 preset with the fit's own alpha-beta margin exposed a
bookkeeping error in `examples/e2_predictions.py`: it reported the lattice energy of each
relaxed conformer cell and left out the strain energy of the refined backbone angles
(`RefineResult.angle_energy`, which the refinement keeps outside the lattice energy by
design). A DFT total energy carries that strain, so the margins compared with DFT must too.
The omission is 0.21 kcal/mol per monomer for beta (whose backbone angle the lattice opens to
116.6 deg) and 0.66 for the refined TGTG' cell, so it favoured the gauche conformers by about
0.45 kcal/mol per monomer. Regenerated with the strain included (both parts are recorded in
the files), the margins over TGTG' are:

| preset | PVDF TTTT | PVDF T3GT3G' | copolymer TTTT | copolymer T3GT3G' |
|---|---|---|---|---|
| pre-registered (strain omitted) | +1.304 | +0.725 | +1.087 | +0.560 |
| `pvdf-dft-valence`, corrected | +0.853 | +0.414 | +0.455 | +0.111 |
| stage 1 (`-crystal`) | +0.982 | +0.522 | +0.528 | +0.179 |
| stage 2 (`-crystal2`) | +0.756 | +0.323 | +0.346 | +0.022 |

The PVDF number is the one the fit path had been reporting all along (-0.854 for E(alpha) -
E(beta) in section 5.1's start column; DFT -0.726). The copolymer's correction is the large one,
because its refined TGTG' packing carries 0.84-1.10 kcal/mol per monomer of angle strain: the
all-trans phase of the 70/30 model is 0.35-0.46 kcal/mol per monomer above its ground state,
not 1.1, and on the stage-2 preset T3GT3G' is degenerate with TGTG'. That is the margin the
switching model (plan section 11) takes as its input, and it is now a third of what the
pre-registered file said. The pre-registered values are kept in
`docs/predictions/e1_e2_notes.md`; the correction goes to Dipole with these numbers.

### 5.5 Charge penetration: implemented, probed, in the fit

The permanent charges can now interact as ``1/r - (1 + b r/2) exp(-b r)/r`` in the lattice
sum (`EwaldSpec(charge_damping=b)`, b in 1/A; `FFParameters.charge_damping` carries it
through `applied()`; the bonded-exclusion correction removes the damped interaction, so
excluded pairs stay excluded; `tests/test_charge_damping.py` checks the off case bit for bit,
the on case against a direct sum, and the gradients against finite differences). Probed on
the valence preset at b = 1.5, 2.0, 3.0 /A: beta's polar axis -5.0% -> -4.3%, -4.7%, -5.0%, the
lateral axis +0.7%, +0.3%, 0.0%, the stiffness down (C_11 19.6, 21.1, 22.3 GPa), the margin
-0.78, -0.81, -0.85. Alone it is not the lever either: it weakens every contact, not the polar
one in particular. Stage 3 (`crystalfit.fit(stage=3)`, eleven variables, running as this is
written) lets the fit weigh it against the other ten; its result goes to section 5.6.

### 5.6 Stage 3: the fit does not want charge penetration

321 evaluations from the valence start with b free in 1-6 /A: objective 112.99 -> 64.86 against
stage 2's 65.05, and b ends at 3.95 /A, where the damping is 3e-4 of the pair interaction at 2.5 A,
i.e. point charges. Every residual is stage 2's to the second digit (beta a -3.9%, b +1.4%, P
0.189, margin -0.767, C_11 27.3, C_22 21.1, C_12 3.4). The result is saved as
`docs/predictions/crystal_stage_v3.json` and not registered as a preset: it is stage 2 with an
eleventh variable the data did not use.

**Where this leaves the accuracy stage.** Three forms have now been fitted to the same anchors in
the setting the predictions are made in: the isotropic Lennard-Jones terms and the bond charges
(stage 1, fixed-charge; its induced-dipole evaluation in 4.1), the same plus a pair-specific H-F
term (stage 2), the same plus charge penetration (stage 3). All three stop at the same place: the
polar axis 4% short, the lateral axis 1.5% long, C_12 at 57% of DFT, with the polarization, the
margin and the polar stiffness within tolerance. A floor that three different forms share is
set by what they have in common: atom-centred point charges (the fluorine's off-site charge of
the valence fit reaches the lattice as its dipole-equivalent projection onto the nucleus, so the
lattice sees no lone-pair anisotropy), isotropic repulsion per atom, and a rigid chain built at
the textbook geometry. The scan of 5.2 puts the excess attraction in the electrostatics of the
H...F contact along the polar axis, and the one lever not yet pulled that acts on that contact
and not on the F...F one is the shape of the fluorine's charge: forwarding the off-site site to
the lattice as a site, with the induced dipoles and the Ewald sum carrying it, instead of
projecting it. With it, the PBE-D3 chain geometry (5.3; 0.6% of the polar deficit), and
Dipole's DFPT dielectric and Born tensors as anchors when they land, that is the next stage. Until
then the pre-registered set on `pvdf-dft-valence` remains the primary prediction and the two
fitted presets its labelled alternatives.

### 5.7 Five probes on the morning of 2026-11-06: what moves the a:b ratio

Each probe keeps the valence preset and the pre-registered construction (induced dipoles,
Ewald, force-shifted Lennard-Jones) and changes one thing; the beta cell is re-relaxed. DFT:
a 4.757, b 8.528, a/b 0.558; the model: 4.520, 8.521, 0.530.

| change | a | b | a/b | reading |
|---|---|---|---|---|
| fluorine's off-site charge carried as a site (static, fixed-charge Ewald) | pull along a 0.948 -> 0.826 kcal/mol/A | - | - | 8% of the pull: not the lever (task 15's first item, dropped) |
| H-F repulsive wall x2, x4 (dispersion untouched) | -0.8%, +4.2% | +4.2%, +9.4% | 0.531, 0.531 | both axes move together; stiffness halves |
| all dispersion x0.7, x0.5 | -2.5%, -0.4% | +2.7%, +4.9% | 0.530, 0.529 | same curve; C_11 14, 10 GPa |
| carbon-pair dispersion x0.5 | -2.6% | +2.3% | 0.531 | same curve |
| chain at the PBE-D3 geometry (5.3) | -4.4% | +0.5% | 0.530 | same curve |
| repulsive site on the C-F axis beyond F, against other chains' H | +0.1% | +9.7% | 0.509 | wrong way: that approach is the diagonal contact |
| repulsive belt sites on F, perpendicular to C-F in the CF2 plane, d 0.7 A, weak / strong | +0.8% / +10.8% | +1.9% / +5.8% | **0.552 / 0.584** | the first change that moves the ratio, and in the right direction |

The contact geometry explains the table. At the model's own minimum the shortest interchain
contacts are F...H at 2.61 A between a chain and its neighbour along a (dx 2.26, dy 0.2, dz 1.3
A; 63 deg off the C-F bond: a sideways approach to the fluorine), then F...H at 2.64 A between
corner and centre chain (dy 2.3, dz 1.3; 46 deg off the bond: nearer head-on). In DFT's cell the
along-a contact is 2.76 A and the diagonal one 2.61. The model compresses the sideways contact
and keeps the head-on one, and every isotropic term, whatever its strength, treats the two
alike, which is why all of them slide the cell along one curve of constant a/b. A site on the
C-F axis stiffens the head-on approach and opens b; sites in the belt perpendicular to the bond
stiffen the sideways approach and open a. That is fluorine's polar flattening, the known
anisotropy of a bonded halogen (the density is drawn in along the bond axis and extended in
the belt), and it is what the point-atom Lennard-Jones form cannot express.

**The next form change is therefore anisotropic fluorine repulsion**, carried as rigid belt
sites on each fluorine (two massless repulsion-only sites in the CF2 plane, perpendicular to the
C-F bond, with a fitted distance and wall), through the polymer's pendant-fragment machinery so
that the chain builder, the packer, the Ewald sum (zero charge) and the induced dipoles (zero
polarizability) carry them unchanged. Task 15 is rewritten to that; the off-site charge as a
site is dropped on the first row of the table.

### 5.8 Stage 4: fluorine belt sites, implemented

The sites are carried by the polymer definition, not by the kernel. :func:`polyfind.polymers.belt_fluorine`
is a rigid pendant fragment: the fluorine at its bond length and two atoms of a pseudo-element
``X`` at ``(bond, +/-d, 0)`` in the pendant frame of :func:`polyfind.chain.pendant_frames`,
whose second axis is exactly the in-plane perpendicular to the C-F bond. :func:`pvdf_belt`,
:func:`trfe_belt` and :func:`vdf_trfe_7_3_belt` (and :func:`belt_variant` for the scripts) build
the polymers; the chain builder, the line groups, the packer, the Ewald sum and the induced
dipoles then carry four extra atoms per CF2 without knowing what they are. What ``X`` is made
to be: a zero well depth in the Lennard-Jones table until a preset gives it one (``UFF_LJ["X"]
= (2.6, 0.0)``), no dispersion term and no interaction at all within its own chain (the packer
zeroes its exclusion scales for every same-site image, and the force field's own pair terms skip
it, so the conformer energetics of the valence fit are untouched), zero charge, a covalent
radius that bond inference never reaches, no polarizability (chargeless and unpolarizable, the
packer leaves it out of the Ewald sums and the dipole solve altogether; until 2026-11-07 it carried
a 1e-6 A^3 placeholder and rode through both, 5.14 has the cost), zero mass. A chain with sites at zero depth
reproduces the plain chain's lattice energy, cell, polarization and stiffness exactly
(`tests/test_belt_sites.py`); the massless sites shift the chain's principal-axis frame by a
hair, which is why a fixed set of cell parameters at a rotated orientation is not bit-identical
and a relaxed cell is.

The wall is a Lennard-Jones repulsion of the ``X-H`` pair only (``X-C``, ``X-F``, ``X-X`` are
zeroed through ``lj_pairs``), its diameter fixed at 2.6 A and its depth the fitted variable, in
log10 (``crystalfit`` stage 4: the ten variables of stage 2 plus ``log10_D_X``, with the charge
penetration of stage 3 carried at its point-charge value; sites at d = 0.7 A). The probe that
sized the start point, on the valence preset with everything else fixed:

| depth D_X, kcal/mol | beta a | beta b | a/b | P | C_11, C_22, C_12 |
|---|---|---|---|---|---|
| 0 (plain chain) | 4.520 (-5.0%) | 8.521 (-0.1%) | 0.530 | 0.140 | 22.6, 18.5, 3.1 |
| 0.01 | 4.750 (-0.1%) | 8.636 (+1.3%) | 0.550 | 0.127 | 17.4, 16.5, 3.1 |
| 0.03 | 4.855 (+2.1%) | 8.708 (+2.1%) | 0.557 | 0.121 | 15.2, 15.3, 3.1 |
| 0.10 | 5.000 (+5.1%) | 8.827 (+3.5%) | 0.566 | 0.114 | 12.6, 13.5, 3.1 |
| DFT | 4.757 | 8.528 | 0.558 | 0.180 | 22.8, 24.5, 5.9 |

At the start point of the fit (D_X = 0.01, the other variables at the valence values) beta's
polar axis is -0.2% and alpha's long axis -0.2% from DFT, the first time either has been within
tolerance; what is left for the other ten variables is the polarization (0.127), the margin
(-1.06, the sites having changed the alpha-beta balance) and the stiffness (17, 16, 3.1 GPa),
which are the quantities the charges and well depths reach. One evaluation is 16 s.

### 5.9 Stage 4 result: the floor is broken

321 evaluations, 85 minutes on the 4-core box, objective 90.4 at the start point -> 16.29, of
which 8.00 is the ridge. Stage 2's 65.05 (4.87 ridge) is the comparison; same setting, same
anchors (the polarization anchor is still the published 0.18, the fit having started before the
Berry repeat of 2026-11-06 landed), two more variables, one of them held:

| quantity | DFT | stage 2 | stage 4 (belt) | z, stage 4 |
|---|---|---|---|---|
| beta a (polar), A | 4.757 | 4.572 (-3.9%) | 4.723 (-0.7%) | -0.71 |
| beta b (lateral), A | 8.528 | 8.656 (+1.5%) | 8.515 (-0.1%) | -0.15 |
| beta c (chain), A | 2.581 | 2.610 (+1.1%) | 2.606 (+1.0%) | +0.97 |
| alpha a, A | 4.945 | 5.031 (+1.7%) | 4.942 (-0.1%) | -0.06 |
| alpha b, A | 9.636 | 9.339 (-3.1%) | 9.595 (-0.4%) | -0.42 |
| alpha c, A | 4.660 | 4.742 (+1.8%) | 4.725 (+1.4%) | +1.39 |
| beta polarization, C/m^2 | 0.18 (Berry 0.193) | 0.189 | 0.179 | -0.09 (-1.4 against 0.193) |
| E(alpha) - E(beta), kcal/mol per monomer | -0.726 | -0.758 | -0.760 | -0.69 |
| E(antipolar) - E(polar) beta, rigid | -0.026 | 0.000 | 0.000 | +0.52 |
| beta C_11 (polar), GPa | 22.77 | 27.30 (+20%) | 25.08 (+10%) | +1.02 |
| beta C_22 (lateral), GPa | 24.52 | 21.03 (-14%) | 24.59 (+0.3%) | +0.03 |
| beta C_12, GPa | 5.88 | 3.36 (-43%) | 4.88 (-17%) | -1.71 |

Every cell axis of both polymorphs is inside the 1% tolerance for the first time, or within
1.4% of it on the two chain axes; the a:b ratio of beta is 0.5547 against DFT's 0.5578 (stage
2: 0.528). The lateral stiffness is exact, the polar one has come down from +20% to +10%, the
coupling C_12 has gone from 57% of DFT to 83%; the margin and the polarization stay inside their
tolerances. What did not move is what no variable here reaches: the antipolar-polar beta margin
of the rigid-chain construction (0.000 at every stage) and the two chain axes, which are the
backbone angle the lattice settles at (section 5.3: 116.6 deg against DFT's 114.6).

The variables, against stage 2's: diameters x_C 1.100 (1.001), x_H 1.027 (1.019), x_F 0.964
(1.011); well depths D_C 1.18 (0.62), D_H 0.72 (1.20), D_F 0.98 (0.94); increments C-H -0.180 e
(-0.173), C-F +0.146 e (+0.157); the H-F pair term at x 1.019 and depth x 1.0005 on its geometric
means, i.e. back at the mixing rule it had been bent away from (1.037, 0.937) when it was the
only term standing in for the anisotropy; the belt wall D_X = 10^-1.705 = 0.0197 kcal/mol, twice
the start point, at its fixed 2.6 A diameter; the charge penetration held at the point-charge
value. Stage 2 had inflated the hydrogen well (x1.20) and deflated carbon's (x0.62) to push the
polar H...F contact open with isotropic terms; stage 4 reverses both (hydrogen softer, carbon
larger and deeper) once the belt wall carries that contact. Half of the ridge is the carbon
diameter: +10% is two widths of its prior (3.96 of the 8.00); the C-H increment is 1.47, the
hydrogen well 0.89, the belt depth itself 0.09.

**Against the Berry polarization.** Read against Dipole's 0.1927 C/m^2 rather than the 0.18 it
was fitted to, the stage-4 P is 7% low. `crystalfit.DFT_ANCHORS["p_beta"]` is 0.1927 since the
repeat landed, so a re-run of stage 4 from this result is the next fit (stage 5 below); whether
the belt form can carry 0.193 without giving back the cell is the question it answers.

**What the two extra variables cost.** One evaluation went from 7 s (stage 2) to 16 s: four
massless sites per CF2 double the atom count the Ewald sum, the induced-dipole system and the
pair kernel see, and the strained-cell stiffness block (six relaxations) is the bulk of each
evaluation. The E1/E2 prediction runs on the belt polymers scale the same way (section 5.10).

### 5.10 The belt preset as a prediction set

`examples/e1_predictions.py --preset pvdf-dft-valence-belt --variant belt` and the E2 script
likewise (`docs/predictions/e1_polyfind_belt.json`, `e2_polyfind_belt.json`). The fit path and the
E1 path agree on the beta cell to 0.005 A and on the alpha-beta margin to 0.002 kcal/mol (the
polar and antipolar alpha cells are degenerate here to that level). Against the two earlier sets
and the DFT in hand:

| quantity | DFT | `pvdf-dft-valence` | stage 2 | belt |
|---|---|---|---|---|
| beta a, b, c, A | 4.757, 8.528, 2.581 | 4.520, 8.521, 2.600 | 4.576, 8.657, 2.606 | 4.728, 8.517, 2.602 |
| beta polarization, C/m^2 | 0.193 (Berry; published 0.18) | 0.140 | 0.189 | 0.179 |
| beta eps_inf (x, y, z) | 2.235, 2.253, 2.600 (DFPT, protocol evidence) | 2.197, 2.133, 2.515 | 2.145, 2.083, 2.469 | 2.108, 2.086, 2.446 |
| beta eps_0 relaxed-ion (x, y, z) | - | 2.197, 2.831, 2.515 | 2.145, 3.253, 2.469 | 2.108, 3.823, 2.446 |
| beta C_11, C_22, C_12, C_66, GPa | 22.8, 24.5, 5.9, - | 22.6, 18.5, 3.1, 4.3 | 27.4, 21.1, 3.3, 4.7 | 25.1, 24.6, 4.9, 5.3 |
| beta proper e_11, C/m^2; d_11, pm/V | - | 0.232; 10.3 | 0.312; 11.4 | 0.281; 11.2 |
| beta intrinsic Q_11, Q_22, m^4/C^2 | - | 1.59, 10.9 | 3.45, 13.4 | 3.83, 8.85 |
| beta d eps_yy / d eta_yy | - | 9.3 | 20.0 | 28.5 |
| alpha (polar TGTG') a, b, c, A | - | 4.977, 9.038, 4.710 | 5.043, 9.276, 4.750 | 4.957, 9.559, 4.732 |
| alpha polarization, C/m^2 | - | 0.073 | 0.095 | 0.094 |
| alpha C_11, C_22, C_12, C_66, GPa | - | 16.6, 12.3, 5.0, 6.0 | 19.2, 11.5, 5.2, 6.1 | 21.6, 14.3, 5.4, 6.7 |
| copolymer a, b, gamma, c per 10 monomers | - | 4.793, 8.852, 92.0, 25.99 | 4.794, 8.963, 92.1, 26.04 | 4.970, 8.401, 88.0, 25.99 |
| copolymer polarization, C/m^2 (tilt from a) | - | 0.109 (24 deg) | 0.151 (23 deg) | 0.153 (38 deg) |
| copolymer eps_0 relaxed-ion (x, y, z) | - | 2.070, 2.371, 2.379 | 2.069, 2.594, 2.361 | 2.411, 2.380, 2.402 |
| copolymer C_11, C_22, C_12, C_66, GPa | - | 16.5, 17.0, 1.2, 1.7 | 19.8, 18.7, 1.6, 2.4 | 24.6, 16.4, 6.1, 6.8 |
| copolymer intrinsic Q_11, m^4/C^2 | - | 3.06 | 3.94 | 1.63 |
| PE (control) a, b, gamma | - | 4.061, 7.088, 90.0 | 4.326, 6.865, 87.8 | 4.139, 7.231, 90.9 |
| PE C_11, C_22, C_12, C_66, GPa | - | 23.3, 25.0, 1.6, 2.4 | 15.9, 14.6, 9.0, 10.3 | 20.4, 17.3, 10.2, 11.5 |

E2, margins over TGTG' in kcal/mol per monomer with the backbone-angle strain included (5.4):

| | DFT | `pvdf-dft-valence` | stage 2 | belt |
|---|---|---|---|---|
| PVDF TTTT | +0.726 | +0.853 | +0.756 | +0.758 |
| PVDF T3GT3G' | - | +0.414 | +0.323 | +0.425 |
| copolymer TTTT | - | +0.455 | +0.346 | +0.172 |
| copolymer T3GT3G' | - | +0.111 | +0.022 | +0.185 |

**Against Dipole's records** (`examples/crosscheck_e1_e2.py --polyfind docs/predictions/e1_polyfind_belt.json`):
the beta cell -0.6%, -0.1%, +0.8% on the QE reference; C_11 at 1.10 of the six-cell value, C_22 at
1.00, C_12 at 0.83; the Berry polarization at 0.93; the electronic permittivity at 0.94, 0.93 and
0.94 of the DFPT tensor on all three axes (that tensor is protocol evidence on Dipole's side, not
an accepted reference, and polyfind's eps_inf is the induced-dipole model's, which every stage has
left 5-7% low). The stage-2 set was at 1.20, 0.86, 0.57, 1.05 (against 0.18) and 0.96, 0.92, 0.95. Dipole's
three-point polar Berry pair, landed the same morning, gives the relaxed-ion longitudinal
piezoelectric slope d|P|/d eta = -0.326 C/m^2 (`docs/predictions/e1_e2_notes.md`); the belt
set's -0.282 is 0.87 of it (pre-registered 0.71, stage 2 0.96), and relative to each set's own
P the slopes are -1.57 against -1.69 (0.93; the other two sets 0.98): the response's shape is
right and its size is P's. That slope, cheap on the fit path (P at the strained states the
stiffness block already relaxes), is the natural next anchor.

**What this set predicts that the others did not.** The lateral relaxed-ion permittivity of beta
rises with every stage, 2.83 -> 3.25 -> 3.82, and its strain derivative with it (9 -> 20 -> 29):
the ionic part is the chains' libration about their axes, which the belt wall leaves softer than
the isotropic forms did while the cell is now right. This is the belt set's sharpest testable
number once Dipole's DFPT eps_0 lands, and it is not an anchor. The polar electronic and
relaxed-ion permittivities are equal on every set (no ionic weight along the dipole in a rigid
chain), which DFPT's polar axis, where a bending mode carries weight, will not be.

**Two cautions.** The copolymer model's packing search did not land in the same arrangement:
under the belt preset the cell is monoclinic the other way (gamma 88 against 92 deg) and the
chain dipoles are tilted 38 deg from a rather than 23, so its rows compare minima, not only the
potential, and that search was already the least converged of the set (`docs/predictions/e1_e2_notes.md`).
Its E2 margin over TGTG' falls to +0.17 kcal/mol per monomer, the closest any set has come to
the known stabilisation of the all-trans phase by TrFE, but for the same reason not a clean
statement. PE is the second: its cell is orthorhombic again (gamma 90.9 against stage 2's 87.8),
but its shear stiffness C_66 is 11.5 GPa against the pre-registered 2.4 and its coupling C_12
10.2 against 1.6, with a non-zero C_16 of 7.2. The carbon diameter grew 10% and the hydrogen
well softened to carry PVDF's contacts, and PE, a chain of nothing but those two, is where that
shows. A PE anchor (Dipole's PE control is in the E1 request) is the next thing the objective
needs, before the belt terms are trusted on a hydrocarbon.

**Cost.** The belt copolymer E1 took 1,356 s on four cores against 290 s for the plain chain
(section 5.3): four sites per CF2 are 92 more atoms on a 120-atom cell, and the Ewald and
induced-dipole work grows with the square. Beta and alpha are 11 and 35 s. The E2 copolymer
T3GT3G' cell, the largest, took 1,011 s.

**The prediction path's refinement and the fit's (found 2026-11-06, closed).** `mechanics.refined_reference`
hands `refine_crystal` only what `refine_kw` says, and the E1/E2 scripts said nothing, so their
chain conformation and starting cell were refined under the refinement's own defaults (the
truncated Coulomb sum, no induced dipoles, the energy-shifted Lennard-Jones) and only then relaxed
in a, b and gamma under the full potential. At the stage-6 point that leaves a and b within 0.1%
of the fit path, c 0.2% short (2.5931 against 2.5978: `relax_reference` cannot move c) and the
backbone angles 0.3 deg closed (116.1 against 116.4). With the full potential passed to the
refinement the two paths agree to four decimals on every cell axis and the angles, and
`relax_reference` has nothing left to move; the scripts do that from here, and
`examples/e1_predictions.py --geometry` builds the belt chain at the geometry the preset was
fitted on. The belt sets above and the pre-registered sets were made before this; the differences
it would make (0.2% in c, 0.3 deg) are below the differences those sets report, so they are left
as they are and the next sets carry the consistent refinement. For the stage-7 preset this is not
optional: its restraint stiffness is fitted against the full potential's lateral pressure, and a
refinement under the truncated potential would relax the angles against a different pressure.

### 5.11 Stage 6: the polar Berry slope as an anchor (closed)

Dipole's three-point Berry pair (5.10; `docs/predictions/e1_e2_notes.md`) gives d|P|/d eta along
the polar axis, -0.3255 C/m^2, the relaxed-ion e_11 on one branch. On the fit path it costs
nothing: the strained states the stiffness block relaxes include the +/-0.2% polar pair, and
`crystalfit.polar_slope` takes the central difference of their true polarization's polar
component with the cell's orientation (the sign of P) removed, which is what makes it comparable
with a Berry-phase finite difference under either side's sign convention. `evaluate(slope=True)`
adds the residual at a tolerance of 0.03 C/m^2 (a tenth of the value, as for the stiffness);
`STAGE6` is stage 2's setting with it, on the twelve belt variables; `examples/crystal_stage.py
--stage 6`. At the stage-4 point the fit path gives -0.2814 against the E1 path's -0.2815 at
+/-1%, so the slope is linear to 1e-4 over that range and the two paths agree on it as they do on
the cell; z = 1.47. With the Berry P anchor (0.1927, z -1.36) the stage-4 point's objective is
20.30 on the stage-6 residual set, against its 16.29 on the stage-4 anchors.

**Stage 5: the Berry P anchor (0.1927 in place of the published 0.18), the twelve belt variables
from the stage-4 point.** The first run was cut by a two-hour session limit at evaluation 292 of
320 (best 12.14 from 18.15 at the stage-4 point on these anchors), and its result file was only
written at the end, which is why the fitter now checkpoints (section 7); the finishing run from
the logged best point converged in 73 evaluations to objective 11.64, of which 6.79 is the ridge
(`docs/predictions/crystal_stage_v5.json`):

| quantity | DFT | stage 4 (belt) | stage 5 (belt, Berry P) | z, stage 5 |
|---|---|---|---|---|
| beta a (polar), A | 4.757 | 4.723 (-0.7%) | 4.763 (+0.1%) | +0.12 |
| beta b (lateral), A | 8.528 | 8.515 (-0.1%) | 8.487 (-0.5%) | -0.48 |
| beta c (chain), A | 2.581 | 2.606 (+1.0%) | 2.599 (+0.7%) | +0.71 |
| alpha a, A | 4.945 | 4.942 (-0.1%) | 4.921 (-0.5%) | -0.48 |
| alpha b, A | 9.636 | 9.595 (-0.4%) | 9.718 (+0.8%) | +0.85 |
| alpha c, A | 4.660 | 4.725 (+1.4%) | 4.707 (+1.0%) | +1.01 |
| beta polarization, C/m^2 | 0.1927 (Berry) | 0.179 (-7%) | 0.188 (-2%) | -0.43 |
| E(alpha) - E(beta), kcal/mol per monomer | -0.726 | -0.760 | -0.726 | 0.00 |
| E(antipolar) - E(polar) beta, rigid | -0.026 | 0.000 | 0.000 | +0.52 |
| beta C_11 (polar), GPa | 22.77 | 25.08 (+10%) | 24.03 (+6%) | +0.55 |
| beta C_22 (lateral), GPa | 24.52 | 24.59 (+0.3%) | 24.74 (+0.9%) | +0.09 |
| beta C_12, GPa | 5.88 | 4.88 (-17%) | 5.19 (-12%) | -1.18 |
| ridge | | 8.00 | 6.79 | |
| objective on these anchors | | 18.15 | 11.64 | |

Every residual is inside 1.2 tolerances. The polarization came up to 2% of the Berry value and
the alpha-beta gap landed on the anchor without the cell giving anything back: the beta axes are
+0.1/-0.5/+0.7%, the alpha axes -0.5/+0.8/+1.0%. The move is in the charges and the hydrogen
well. The C-H and C-F increments grow (-0.068 e and +0.043 e against stage 4's -0.061 and
+0.032), which is where the extra polarization comes from; the hydrogen well returns to the
valence value (x0.98 from x0.72) and the carbon diameter's excess halves (+6% from +10%), the
two terms that carried half the stage-4 ridge; the fluorine well drops to x0.84 (x0.98) and the
carbon well to x1.12 (x1.18); the H-F pair term (x1.019, x0.997), the belt wall (10^-1.624 =
0.0238 kcal/mol) and the held penetration parameter are as before. The ridge falls to 6.79. The
one systematic left is the chain axes, +0.7% and +1.0%, the backbone angle again.

**On the stage-6 residual set** the stage-5 point gives d|P|/d eta = -0.2930 C/m^2 (z 1.08) against
the stage-4 point's -0.2814 (z 1.47): the slope moved with the polarization, and the slope per unit
polarization is -1.56 against Dipole's -1.69. Its objective there is 12.81.

**The PBE-D3 chain geometry at the same point** (`--geometry pbe`, one evaluation, 30 s, before
any refit): Dipole's bond lengths and angles make every chain axis longer rather than shorter (beta
c 2.615, +1.3%; alpha c 4.734, +1.6%), widen the alpha cell (b 9.787, +1.6%) and deepen the
alpha-beta gap (-0.790, z -1.28); the slope is unchanged (-0.2933) and the objective is 19.00. With
the C-C bond at 1.533 A, c = 2.615 is a backbone angle near 117 deg against the 114.6 the chain was
built at and DFT settles at (c = 2.581 is exactly the planar zigzag at 1.533 A and 114.6 deg). So the
chain-axis excess is not the textbook bond geometry's: the lattice opens the backbone angle on
either geometry, and a shorter C-C bond with longer C-F bonds opens it further. What opens it is
the next probe (5.12).

**Stage 6: both Berry anchors.** 301 evaluations, 100 minutes, from the stage-5 point: objective
12.81 to 11.19 on the stage-6 residual set, of which 5.00 is the ridge (6.79 before);
`docs/predictions/crystal_stage_v6.json`:

| quantity | DFT | stage 5 | stage 6 (belt, P and slope) | z, stage 6 |
|---|---|---|---|---|
| beta a (polar), A | 4.757 | 4.763 (+0.1%) | 4.750 (-0.1%) | -0.14 |
| beta b (lateral), A | 8.528 | 8.487 (-0.5%) | 8.479 (-0.6%) | -0.57 |
| beta c (chain), A | 2.581 | 2.599 (+0.7%) | 2.598 (+0.7%) | +0.65 |
| alpha a, A | 4.945 | 4.921 (-0.5%) | 4.913 (-0.6%) | -0.65 |
| alpha b, A | 9.636 | 9.718 (+0.8%) | 9.724 (+0.9%) | +0.92 |
| alpha c, A | 4.660 | 4.707 (+1.0%) | 4.700 (+0.9%) | +0.86 |
| beta polarization, C/m^2 | 0.1927 (Berry) | 0.188 (-2%) | 0.192 (-0.5%) | -0.10 |
| E(alpha) - E(beta), kcal/mol per monomer | -0.726 | -0.726 | -0.735 | -0.19 |
| E(antipolar) - E(polar) beta, rigid | -0.026 | 0.000 | 0.000 | +0.52 |
| beta C_11 (polar), GPa | 22.77 | 24.03 (+6%) | 24.16 (+6%) | +0.61 |
| beta C_22 (lateral), GPa | 24.52 | 24.74 (+0.9%) | 24.89 (+1.5%) | +0.15 |
| beta C_12, GPa | 5.88 | 5.19 (-12%) | 5.07 (-14%) | -1.38 |
| beta d\|P\|/d eta (polar), C/m^2 | -0.3255 | -0.2930 | -0.2986 | +0.90 |
| ridge | | 6.79 | 5.00 | |
| objective on the stage-6 set | | 12.81 | 11.19 | |

Both Berry anchors are met without the cell or the gap giving anything back: the polarization is
within half a percent of Dipole's and the slope within a tenth. The slope per unit polarization is
-1.56 against Dipole's -1.69, as at every earlier point, so the tenth that is left is the shape of
the response and the fit took the slope through the polarization. The ridge fell with the carbon
diameter: x_C 1.016 from 1.060 (1.100 at stage 4), the term that had been half of the stage-4
ridge; the polarization is now carried by the increments (C-F +0.053 e, C-H -0.057 e against
+0.043 and -0.068) with the hydrogen well above the valence value again (x1.12 from x0.98), the
fluorine well at x0.84, the carbon well at x1.14, the H-F pair at (1.023, 1.007) and the belt wall
at 10^-1.656 = 0.0221 kcal/mol. What is left is C_12 (-14%) and the chain axes (+0.7% and +0.9%),
and 5.12 says where the chain axes come from. The result is not registered as a preset: stage 7
(5.13) decides which of the two is.

### 5.12 The backbone angle: what opens it, and what DFT says it is

The refinement (`polyfind.refine.refine_crystal`) relaxes the two backbone angles of the repeat,
at the CH2 and the CF2 carbon, against a harmonic restraint of 105 kcal/mol/rad^2 about the angle
the chain was built at (114.0 textbook, 114.6 PBE-D3). The valence preset's fitted bend terms are
deliberately not forwarded to it (`FFParameters.applied` says why: for the rigid chains the packer
sees, the valence energy is a constant), so 105 is the one invented number left in the crystal
stage, and its strain is part of the energy per monomer the fit compares: 0.21 kcal/mol for beta
and 0.67 for alpha at the stage-5 point, i.e. 0.46 of the alpha-beta gap is restraint strain.

**The fitted bend terms cannot replace it.** An isolated TTTT chain relaxed under the valence
preset's own terms (`forcefield.relax_backbone_angles`, bend + intrachain nonbonded, the 1-4 scaling
0.92 they were fitted with) closes its angles to 113.8 and 107.3 deg (111.6 and 107.9 at the
stage-5 Lennard-Jones terms): a chain axis near 2.50 A, 3% short. The C-C-C term (k 77, theta_0
108.8) was fitted to oligomer energetics with the nonbonded terms doing the rest, and the rest
does not reproduce the backbone angle. The restraint is the better of the two.

**The chain's own nonbonded terms do not open it either.** Setting both angles together on the
isolated chain at the stage-5 terms, per monomer per degree near 116 deg: the Lennard-Jones term
falls by 0.10 kcal/mol, the Coulomb term rises by 0.08. They cancel. What opens the angle is
intermolecular.

**Frozen angles.** The same refinement with the angle freedom's bound set to zero
(`max_angle_change=0`; the line-group parametrisation bounds the angles, it does not take the
`refine_angles` flag), stage-5 parameters:

| | a | b | c | E_lattice + strain, kcal/mol per monomer | angles (CH2, CF2) | gauche torsion |
|---|---|---|---|---|---|---|
| beta, PBE-D3 geometry, refined | +0.6% | +0.2% | +1.3% | -10.037 + 0.187 = -9.850 | 117.0, 117.0 | |
| beta, PBE-D3 geometry, frozen | +1.5% | +0.5% | -0.0% | -9.643 | 114.6, 114.6 | |
| alpha, PBE-D3 geometry, refined | +0.2% | +1.6% | +1.6% | -11.222 + 0.583 = -10.639 | 119.5, 118.1 | 59.2 |
| alpha, PBE-D3 geometry, frozen | +3.6% | -0.4% | -0.2% | -9.552 | 114.6, 114.6 | 73.1 |
| beta, textbook, refined | +0.1% | -0.5% | +0.7% | -10.019 + 0.206 = -9.812 | 116.5, 116.5 | |
| beta, textbook, frozen | +1.1% | -0.2% | -0.7% | -9.581 | 114.0, 114.0 | |
| alpha, textbook, refined | -0.5% | +0.8% | +1.0% | -11.211 + 0.672 = -10.538 | 119.3, 117.7 | 59.5 |
| alpha, textbook, frozen | +3.2% | -0.8% | -0.9% | -9.286 | 114.0, 114.0 | 74.4 |

With the angles held, beta's chain axis is exact on the PBE-D3 geometry and the lateral pressure
appears where it was hiding: the polar axis goes to +1.5% (+1.1% textbook) for 0.21 kcal/mol per
monomer. Alpha pays five times that, 1.09 kcal/mol (1.25 textbook), its polar axis goes to +3.6%
and its gauche torsion to 73 deg; and E(alpha) - E(beta) turns from -0.79 (-0.73 textbook) to
+0.09 (+0.30): **the alpha-beta gap anchor is being matched through the alpha chain opening its
angles against the invented restraint**, not through the lattice energies of two chains at the
geometry DFT gives them.

**What DFT gives them.** Dipole's native beta cell (`pvdf_beta_ref_v1`) is a planar zigzag:
dihedral 180.00, both backbone angles 114.63 deg, C-C 1.5331 A, so that 2 r sin(theta/2) = 2.5807
is its chain axis exactly. The alpha reference's start geometry (the PBE-D3(BJ) relaxation of the
polar-margin deliverable, at the anchor cell 4.945/9.636/4.660; Dipole's native relaxation is
running) has backbone angles 118.23 at the CH2 carbon and 116.87 at the CF2 carbon, torsions
179.1 and 59.9, C-C 1.527 and 1.522. So alpha's opening is real, +3.6 and +2.2 deg from beta, and
polyfind has the direction, the ordering (CH2 wider than CF2) and the torsions (179.0, 59.5)
right; what it has wrong is the amount, both polymorphs 1 to 2.4 deg too open: beta 116.5-117.0
against 114.6, alpha 119.3-119.5 and 117.7-118.1 against 118.2 and 116.9. The restraint's
stiffness decides how much lateral pressure goes into the chain axes and how much of the gap the
alpha opening buys, and three DFT angles now pin it. That is stage 7.

### 5.13 Stage 7: the restraint as a variable, the backbone angles as anchors (closed; a lower minimum of the same problem is in 5.18)

The thirteenth variable is log10 of the restraint's multiplier on the refinement's default
(`refine.DEFAULT_ANGLE_STIFFNESS`, 105 kcal/mol/rad^2), bounded at a fifth and five times, with a
ridge width of a factor two. It lives on the parameter set (`FFParameters.angle_stiffness`) and
`applied()` binds it for the refinement (`refine.ANGLE_STIFFNESS`), so every refinement made under
a fitted preset uses the fitted value, the E1/E2 prediction path included; the valence preset and
the earlier crystal presets carry `None` and refine as before. The anchors are the three angles
above (beta's mean of two, alpha's CH2 and CF2) at a tolerance of 1.0 deg, about half of what the
1% chain-axis tolerance already allows (1% of c is 1.8 deg of backbone angle). Stage 7 is stage 6
with them, on the PBE-D3 chain geometry, which builds the chain at beta's DFT angle so that the
restraint's equilibrium is the DFT value and the bond lengths are Dipole's.

Two evaluations at the stage-5 point on that geometry show what the fit has to resolve:

| quantity | k = 105 | k = 210 | DFT |
|---|---|---|---|
| beta backbone angle | 117.02 (z 2.39) | 115.87 (1.24) | 114.63 |
| alpha angle at CH2 | 119.50 (1.27) | 117.92 (-0.31) | 118.23 |
| alpha angle at CF2 | 118.13 (1.26) | 116.86 (-0.01) | 116.87 |
| beta c, alpha c | +1.3%, +1.6% | +0.7%, +0.8% | |
| beta a, alpha a (polar) | +0.6%, +0.2% | +1.0%, +1.2% | |
| E(alpha) - E(beta) | -0.790 (-1.28) | -0.500 (4.52) | -0.726 |
| objective | 27.90 | 37.70 | |

Doubling the restraint puts alpha's angles on DFT's and halves the chain-axis excess, and the
gap opens to -0.50 while the polar axes grow 1%: once the chains cannot stretch, the lateral
potential at the stage-5 point is too repulsive and alpha is not bound enough against beta.
Whether the other twelve variables can hold the cell, the Berry anchors and the gap with the
angles at DFT's values was the question stage 7 had to answer, from the stage-6 result (5.11) on
the PBE-D3 geometry.

**Result.** Stage 7 closed at objective 11.65 from 25.14 at the stage-6 point (281 evaluations at 16
s; `crystal_stage_v7.json`), ridge 4.45, the restraint fitted at 178 kcal/mol/rad^2 (multiplier
1.70 on the default). The twelve other variables did hold everything: |P| 0.1928 against 0.1927,
the slope -0.300 against -0.326 (0.92; z 0.84), the gap -0.726 on the anchor (the k = 210
evaluation above had it at -0.50), every cell axis within 0.8%, C_11 23.5 and C_22 24.2 against
22.8 and 24.5. Alpha's angles are on DFT's (118.05 and 116.94 against 118.23 and 116.87); beta's
closed from 116.93 to 116.00 and is the largest angle residual left (z 1.37), and since the C-C
bond is now Dipole's its chain axis is that angle: 2 x 1.5331 x sin(116.00/2) = 2.6003 A, +0.75%.
The comparison with stage 6 is made under the same residual set, the stage-6 point evaluated on
its own (textbook) geometry with the three angle rows added (`evaluate(x6, **STAGE7)`, restraint at
the default 105):

| quantity | stage 6 (textbook, k = 105) | stage 7 (PBE-D3, k = 178) | DFT |
|---|---|---|---|
| beta a, b, c | 4.7503, 8.4790, 2.5978 (-0.14, -0.57, 0.65) | 4.7445, 8.4756, 2.6003 (-0.26, -0.61, 0.75) | 4.757, 8.528, 2.581 |
| alpha a, b, c | 4.9127, 9.7244, 4.6999 (-0.65, 0.92, 0.86) | 4.9264, 9.7073, 4.6815 (-0.38, 0.74, 0.46) | 4.945, 9.636, 4.660 |
| beta Berry P | 0.1917 (-0.10) | 0.1928 (0.01) | 0.1927 |
| E(alpha) - E(beta) | -0.735 (-0.19) | -0.726 (-0.01) | -0.726 |
| antipolar - polar (rigid) | -0.000 (0.52) | 0.000 (0.52) | -0.026 |
| C_11, C_22, C_12 | 24.16, 24.89, 5.07 (0.61, 0.15, -1.38) | 23.46, 24.19, 4.99 (0.30, -0.14, -1.51) | 22.77, 24.52, 5.88 |
| d|P|/d eta (polar) | -0.299 (0.90) | -0.300 (0.84) | -0.326 |
| beta backbone angle | 116.44 (1.81) | 116.00 (1.37) | 114.63 |
| alpha angle at CH2, CF2 | 119.17, 117.73 (0.94, 0.86) | 118.05, 116.94 (-0.18, 0.07) | 118.23, 116.87 |
| ridge, objective | 5.00, 16.08 | 4.45, 11.65 | |

Stage 7 is the better point on eleven of the sixteen rows, with the rigid antipolar margin a tie.
The four it loses are small, beta's three axes by at most a tenth of a tolerance and C_12 by 0.13,
and the a and c losses are one residual: a beta angle 1.4 deg open is a chain axis 0.75% long, and
the lateral pressure that opens it (5.12) is still there, resisted by 178 instead of 105. Beta's
lateral axis sits 0.6% short in both sets and C_12 0.9 GPa low in both; the belt wall that fixed
the a:b ratio (5.7) is a hard wall, and the coupling between the two lateral axes is what a hard
wall gives least of. Those two rows are the next thing the form, not the fit, has to change. In
parameters the restraint's rise was paid for by the carbon diameter (3.624 to 3.465 A, x 0.971 of
the valence value), hydrogen 2.651 to 2.572, fluorine 3.333 to 3.369 and the H-F pair 3.042 to
3.004; the increments barely moved (C-H -0.176 to -0.162, C-F +0.167 to +0.169). This point is
the preset `pvdf-dft-valence-belt2` (5.14).

### 5.14 The stage-7 preset as a prediction set (`pvdf-dft-valence-belt2`)

`examples/e1_predictions.py --preset pvdf-dft-valence-belt2 --variant belt --geometry pbe` and the
E2 script likewise (`docs/predictions/e1_polyfind_belt2.json`, `e2_polyfind_belt2.json`), with the
refinement under the full potential (5.10, last paragraph): the E1 path's beta cell is the fit
path's to four decimals (4.7445, 8.4756, 2.6003) and so is alpha's (4.9264, 9.7073, 4.6815), and
E2's PVDF TTTT margin is the fit's gap, +0.726. Against the belt set (5.10 has the two earlier
sets beside it) and the DFT in hand:

| quantity | DFT | belt (stage 4) | belt2 (stage 7) |
|---|---|---|---|
| beta a, b, c, A | 4.757, 8.528, 2.581 | 4.728, 8.517, 2.602 | 4.745, 8.476, 2.600 |
| beta polarization, C/m^2 | 0.193 (Berry) | 0.179 | 0.193 |
| beta eps_inf (x, y, z) | 2.235, 2.253, 2.600 (DFPT, protocol evidence) | 2.108, 2.086, 2.446 | 2.122, 2.095, 2.448 |
| beta eps_0 relaxed-ion (x, y, z) | - | 2.108, 3.823, 2.446 | 2.122, 4.386, 2.448 |
| beta C_11, C_22, C_12, C_66, GPa | 22.8, 24.5, 5.9, - | 25.1, 24.6, 4.9, 5.3 | 23.5, 24.2, 5.0, 5.5 |
| beta d\|P\|/d eta (polar), C/m^2 | -0.326 (Berry pair) | -0.282 | -0.300 |
| beta proper e_11, C/m^2; d_11, pm/V | - | 0.281; 11.2 | 0.300; 12.8 |
| beta intrinsic Q_11, Q_22, m^4/C^2 | - | 3.83, 8.85 | 5.58, 11.56 |
| beta d eps_yy / d eta_yy | - | 28.5 | 48.6 |
| alpha (polar TGTG') a, b, c, A | 4.945, 9.636, 4.660 (the antipolar anchor; degenerate here) | 4.957, 9.559, 4.732 | 4.926, 9.707, 4.682 |
| alpha polarization, C/m^2 | - | 0.094 | 0.099 |
| alpha C_11, C_22, C_12, C_66, GPa | - | 21.6, 14.3, 5.4, 6.7 | 21.1, 15.0, 5.2, 6.2 |
| copolymer a, b, gamma, c per 10 monomers | - | 4.970, 8.401, 88.0, 25.99 | 4.959, 8.397, 92.1, 25.89 |
| copolymer polarization, C/m^2 (tilt from a) | - | 0.153 (38 deg) | 0.166 (85 deg) |
| copolymer eps_0 relaxed-ion (x, y, z) | - | 2.411, 2.380, 2.402 | 2.769, 2.140, 2.417 |
| copolymer C_11, C_22, C_12, C_66, GPa | - | 24.6, 16.4, 6.1, 6.8 | 12.5, 21.1, 8.8, 9.7 |
| copolymer intrinsic Q_11, Q_22, m^4/C^2 | - | 1.63, -0.60 | -3.27, -0.78 |
| PE (control) a, b, gamma | - | 4.139, 7.231, 90.9 | 3.971, 6.887, 90.0 |
| PE eps_inf (x, y, z) | - | 2.771, 2.754, 2.904 | 3.038, 3.031, 3.094 |
| PE C_11, C_22, C_12, C_66, GPa | - | 20.4, 17.3, 10.2, 11.5 | 31.7, 33.8, 1.0, 4.9 |

E2, margins over TGTG' in kcal/mol per monomer with the backbone-angle strain included (5.4):

| | DFT | belt | belt2 |
|---|---|---|---|
| PVDF TTTT | +0.726 | +0.758 | +0.726 |
| PVDF T3GT3G' | - | +0.425 | +0.333 |
| copolymer TTTT | - | +0.172 | +0.059 |
| copolymer T3GT3G' | - | +0.185 | +0.030 |

**Against Dipole's records** (`examples/crosscheck_e1_e2.py --polyfind docs/predictions/e1_polyfind_belt2.json --e2 docs/predictions/e2_polyfind_belt2.json`):
the beta cell -0.3%, -0.6%, +0.8% on the QE reference (volume -0.1%); C_11 at 1.03 of the six-cell
value, C_22 at 0.99, C_12 at 0.85; the Berry polarization at 1.00; the polar Berry slope at 0.92,
and 0.92 relative to its own P as well (the belt set was 0.87 and 0.93); the electronic
permittivity at 0.95, 0.93 and 0.94 of the DFPT tensor, as every set has been, because no stage
has touched the polarizabilities. Every anchored quantity is where the fit put it, which is what
this check is for; the rows that are not anchors are the predictions.

**What this set predicts.** The ionic part of beta's lateral permittivity is the chains' libration
about their axes, and it has risen with every stage: eps_0 - eps_inf along y is 0.70, 1.17, 1.74,
2.29 (pre-registered, stage 2, belt, belt2). For a rigid-dipole libration that increment scales as
P^2 over the libration stiffness, and P is pinned by the Berry anchor, so the stiffness this form
implies is 0.60 of the isotropic forms' (belt) and 0.53 (belt2): the fluorine wall is hard where
the chains touch and leaves them free to rock between contacts. Its strain derivative says the
same, 49 per unit lateral strain against 29 and 20: the libration softens fast as the chains
separate. Dipole's DFPT eps_0 along the lateral axis is the number that tests this, and it is not
an anchor. The piezoelectric column follows P: proper e_11 0.300 C/m^2 (0.92 of the Berry slope,
which on one branch is this coefficient), d_11 12.8 pm/V, intrinsic Q_11 5.6 and Q_22 11.6
m^4/C^2 (3.8 and 8.9 on the belt set; the lateral Q is the libration's again).

**Two cautions, both sharper than in 5.10.** The copolymer model's packing search landed in a
third arrangement: the cell is back in the pre-registered setting (gamma 92.1 against the belt
set's 88.0) but the chain dipoles lie along b, 85 deg from a, so the frame's "polar" axis x is now
across the dipole and its rows are read in the dipole's frame (5.17 has the same cell refined in
the pre-registered set's arrangement, where the rows line up with that set's): the longitudinal proper
coefficient is e_22 = -0.253 C/m^2 (1.53 times P, the rigid-dipole dilation term of
`e1_e2_notes.md`), the stiffness along the dipole 21.1 GPa and across it 12.5, and the libration's
permittivity sits on x (2.77 - 2.04 = 0.73), the axis perpendicular to the dipole, as it must. Its
E2 margins put the three conformations within 0.06 kcal/mol per monomer, TTTT +0.059 over TGTG':
the TrFE stabilisation of the all-trans phase is nearly complete on this set, but the search's
minimum differs from set to set, so it is not a clean statement, and the copolymer rows compare
minima before they compare potentials. PE is the second and is worse than before: its cell volume
is 69.7 A^3 against 76.3 (belt) and 73.4 (pre-registered), its electronic permittivity 3.04 against
2.77, its in-plane stiffness 31.7 and 33.8 GPa against 20.4 and 17.3, with C_12 1.0 and C_66 4.9
(10.2 and 11.5 on the belt set). The restraint's rise was paid for by the carbon diameter (1.016
to 0.971 of the valence value), the hydrogen diameter (1.021 to 0.990) and the hydrogen well (1.116
to 1.012 of it), and PE, a chain of nothing but those two, is where it shows. Dipole's PE cell is
the anchor this objective needs next, before any belt set is trusted on a hydrocarbon.

**Cost.** Beta 28 s, PE 19 s, alpha 70 s, the copolymer E1 1,719 s and the E2 copolymer T3GT3G'
cell 967 s, all with the E1 and E2 runs and the test suite sharing the four cores (the belt set's
copolymer E1 took 1,356 s alone, 5.10). Since 2026-11-07 the packer leaves the belt sites out of the
Ewald sums and the dipole solve, being chargeless and unpolarizable: one energy-and-gradient call on
the belt copolymer cell fell from 0.64 to 0.42 s against 0.17 s for the plain chain, the Ewald share
back to the plain chain's and the remaining excess the Lennard-Jones pair arrays over 212 sites
(`docs/PERFORMANCE_REVIEW.md` 16.6). The sets above were made before this; the change moves no
prediction by more than 4e-6 relative (the stage-7 point re-evaluated), and the copolymer E1 cell
re-run with it takes 776 s on a quiet machine (601 s of build and tensors) against 1,356 s for the
belt set alone, its cell, polarization and permittivities unchanged and its stiffness within 2e-4
relative, which is the relaxation's own noise floor.

**What is left.** In beta, the lateral axis 0.6% short and C_12 0.9 GPa low are shared by every
belt set and did not move between stages 6 and 7: they are the coupling between the two lateral
axes, which a hard wall gives least of, and the form's next change. The backbone angle of beta
is 1.4 deg open against a restraint the fit itself set, so the chain-axis excess (+0.8%) is now
the lateral pressure's, not an invented constant's. Dipole's native alpha positions replace the
start geometry's angles as anchors when they land; the DFPT tensors (eps_0 lateral above all) and
the PE cell are the anchors that would let a stage 8 move the form rather than the fit.

### 5.15 Sensitivity at the stage-7 point: what the form can still move

The Jacobian of the sixteen residuals (tolerance units) with respect to the twelve free variables at
the stage-7 point, by central differences on the PBE-D3 belt geometry
(`examples/crystal_stage_sensitivity.py`, 25 evaluations, `docs/predictions/crystal_stage_v7_sensitivity.json`),
scaled by the ridge's prior widths so that an entry of 1 means one width of a variable moves the
residual by one tolerance:

| residual | x_C | x_H | x_F | D_C | D_H | D_F | q_CH | q_CF | x_HF | D_HF | log10_D_X | log10_k_angle | z at the point |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| beta a (polar) | -0.13 | 1.58 | 0.32 | -0.36 | 0.73 | -0.13 | 0.90 | -1.14 | 0.57 | 0.10 | 5.16 | 0.34 | -0.26 |
| beta b (lateral) | -0.04 | 1.62 | 1.40 | -0.28 | 0.72 | 0.38 | 0.66 | -1.46 | 2.05 | 1.27 | 2.37 | 0.10 | -0.61 |
| beta c (chain) | 0.03 | 0.00 | 0.18 | 0.01 | 0.00 | 0.09 | -0.12 | 0.14 | -0.02 | -0.01 | 0.00 | -0.50 | +0.75 |
| alpha a | 0.13 | 1.77 | 1.06 | -0.19 | 0.80 | 0.24 | 1.06 | -1.85 | 2.64 | 1.71 | 1.83 | 0.88 | -0.38 |
| alpha b | -0.46 | 0.87 | 0.43 | -0.47 | 0.36 | -0.09 | -0.04 | 0.11 | -0.42 | -0.56 | 4.39 | -0.79 | +0.74 |
| alpha c | 0.18 | 0.44 | 0.32 | 0.05 | 0.21 | 0.14 | 0.09 | -0.13 | 0.87 | 0.65 | -0.08 | -0.62 | +0.46 |
| beta Berry P | 0.04 | -0.89 | -0.49 | 0.18 | -0.40 | -0.08 | -2.69 | 4.21 | -0.69 | -0.35 | -2.15 | -0.01 | +0.01 |
| E(alpha) - E(beta) | 1.22 | 2.95 | -2.22 | 0.99 | 1.69 | -0.34 | 2.79 | -4.01 | 9.26 | 7.21 | -7.32 | 5.06 | -0.01 |
| beta C_11 | 0.53 | -0.48 | 0.22 | 0.90 | 0.08 | 0.76 | -1.74 | 2.31 | 0.49 | 1.07 | -3.35 | 0.05 | +0.30 |
| beta C_22 | 0.45 | -0.62 | -0.45 | 0.78 | 0.02 | 0.47 | -1.35 | 2.80 | -0.58 | 0.27 | -1.76 | 0.06 | -0.14 |
| beta C_12 | 0.25 | -1.06 | -0.57 | 0.33 | -0.35 | 0.00 | -1.82 | 3.10 | -2.68 | -1.61 | 0.99 | -0.04 | -1.51 |
| beta d\|P\|/d eta (polar) | -0.03 | 0.54 | 0.21 | -0.11 | 0.25 | 0.00 | 1.49 | -2.17 | 0.28 | 0.10 | 1.58 | 0.00 | +0.84 |
| beta backbone angle | 0.06 | 0.00 | 0.33 | 0.03 | 0.00 | 0.16 | -0.22 | 0.26 | -0.04 | -0.02 | -0.01 | -0.91 | +1.37 |
| alpha angle at CH2 | 0.14 | 0.36 | 0.82 | 0.05 | 0.19 | 0.39 | -0.09 | 0.04 | 0.71 | 0.55 | 0.02 | -1.34 | -0.18 |
| alpha angle at CF2 | -0.02 | 0.31 | 0.25 | -0.02 | 0.15 | 0.13 | -0.20 | 0.28 | 0.47 | 0.38 | 0.10 | -1.04 | +0.07 |

The rigid antipolar-polar margin row is zero in every column: under this form the two rigid
packings are degenerate whatever the parameters, and its z of 0.52 is a floor the fit cannot
touch. Three things the table says:

**Six directions of twelve.** The singular values of the scaled matrix are 17.3, 11.6, 4.7, 3.4,
1.35, 1.19, then 0.31, 0.29, 0.06, 0.04, 0.02, 0.00: sixteen residual rows constrain six directions
of the twelve-variable space to within a prior width, and the ridge decides the other six. That is
the parameter count against the observation count, stated: the anchors in hand determine half the
form, and every statement about a residual being "the form's" has to be read with the other half
set by the prior.

**No residual moves alone.** The minimum-norm step that moves one residual by one tolerance while
holding the fifteen others is tens to hundreds of widths for every row (14.7 for beta a, 68 for
C_12, 79 for beta b, 262 for P, 501 for the slope), almost all of it along the softest directions
(the hydrogen well depth): the anchors are coupled through the same few combinations of
variables, as a physical potential's should be. The lateral coupling in particular, C_12 and beta
b, is moved by the variables that set P (q_CF 3.10 and -1.46, q_CH -1.82), the gap (x_HF -2.68 and
2.05, D_HF -1.61 and 1.27, the wall depth 0.99 and 2.37) and the slope; the carbon diameter the
earlier stages leaned on gives C_12 0.25 per width and b nothing.

**The trade.** Asking the linear model for the best C_12 and beta b with every other residual held
within one tolerance and the variables inside their bounds puts both on their anchors, by x_C +2.2
widths (0.971 to 1.079), the restraint +0.7 widths (178 to 296 kcal/mol/rad^2), D_F -0.7 widths,
q_CH and x_HF -0.4, with three rows pinned at their limit: beta a, the gap and the beta backbone
angle. Evaluated for real at that point (`--verify`): C_12 z -0.01 and b +0.12 as predicted,
beta a +1.18 and the angle +1.05 close to the prediction, but the gap +2.31 against the model's
+0.98, the objective 22.0 (ridge 9.5) against 11.65. The gap is the nonlinear row, with the largest
sensitivities in the table (9.3 per width of the H-F diameter, 7.2 of its depth, 7.3 of the wall
depth), and it is the binding constraint: within this form the lateral coupling is bought with the
alpha-beta gap and beta's polar axis, there is no direction that improves C_12 and b without paying
in those, and the stage-7 optimum is where the tolerances and the ridge put that trade. So the
residual is the form's in the sense that matters: stage 8 has to add a direction, a term that
couples the two lateral axes without moving the gap, which the candidates of 5.14 (a softer or
longer-ranged wall, a C-F bond dipole or quadrupole, an H..F pair attraction) will be screened
against with this Jacobian before any is fitted.

*Correction (5.18).* This analysis is of the neighbourhood of Powell's point, which turned out to
be a local minimum: a least-squares run on the same objective found a point at 9.11 where C_12 is
at z -0.92 and the gap at -0.21, with the carbon diameter at the valence value. The qualitative
statement stands (C_12 and b move with the variables that set P, the gap and the slope; six
constrained directions of twelve), the quantitative trade-off does not: at the least-squares point
the same analysis finds the lateral coupling held by P, the polar axis and the ridge together, not
by the gap (5.18), and that Jacobian (`crystal_stage_v7_trf2_sensitivity.json`) is the one to
screen new terms with.

### 5.16 The relaxed-ion dielectric tensor by the Gamma route: the like-for-like for DFPT

Every eps_0 in the prediction sets is `mechanics.dielectric_tensor(relax=True)`: the chains
reorient under the field at fixed cell, which is the ionic response a rigid-chain parametrisation
can express. DFPT's relaxed-ion tensor lets every atom move. `examples/gamma_dielectric.py` makes
the comparison like for like: it builds the E1 reference on the fit path, gives the packer the
valence stretch and bend terms of `pvdf-dft-valence` pinned at that geometry (zero force at the
reference, the valence fit's stiffness for the stretches and the pendant bends) with the two
backbone angles carrying the preset's own restraint about the built angle instead (the term the
refinement balanced the lateral pressure against, so the reference is stationary over them too),
relaxes every atom at fixed cell, and forms eps_0 = eps_inf + (4 pi k_e / V) Z^T H^+ Z from the
all-atom Hessian, the Born tensors by the field route and the clamped-ion tensor, the translations
projected out. The belt sites are not coordinates: each is the fixed affine combination of its
carbon and that carbon's two fluorines that places it in the CF2 plane, which carries it exactly
through any rigid motion of the group (a first version that let each site follow its fluorine
alone misplaced it under the libration by the site's 0.7 A offset and gave 3.72 on the lateral
axis; the exact construction gives the number below).

| | pre-registered (`pvdf-dft-valence`, plain chain) | belt2 (`pvdf-dft-valence-belt2`, PBE-D3 belt chain) |
|---|---|---|
| eps_inf (x polar, y lateral, z chain) | 2.209, 2.135, 2.495 | 2.126, 2.095, 2.440 |
| eps_0, rigid-chain route (the E1 sets) | 2.196, 2.832, 2.515 | 2.122, 4.386, 2.448 |
| eps_0, Gamma route | 2.272, 2.870, 2.648 | 2.236, 4.573, 2.735 |
| ionic part, Gamma route | 0.063, 0.735, 0.153 | 0.110, 2.478, 0.295 |
| libration: frequency, its share of the lateral ionic part | 63.6 cm^-1, 0.687 | 46.5 cm^-1, 2.381 |
| other lateral carriers | 249 cm^-1 0.020, 891 cm^-1 0.022 | 248 cm^-1 0.043, 883 cm^-1 0.045 |
| polar carriers | 1063, 406, 335 cm^-1 (0.024, 0.017, 0.021) | 432, 331, 1176 cm^-1 (0.042, 0.034, 0.032) |
| chain-axis carriers | 200 cm^-1 0.091, 833 cm^-1 0.063 | 204 cm^-1 0.177, 823 cm^-1 0.118 |
| lowest six modes, cm^-1 | 31, 34, 40, 64, 77, 200 | 20, 41, 47, 50, 101, 204 |
| all-atom relaxation: max force before/after, displacement | 10.1 -> 6e-11 kcal/(mol A), 0.012 A | 6.5 -> 2e-12, 0.006 A |

The two routes agree on what the rigid-chain route was built to carry: the lateral axis is the
libration on both presets (0.69 of 0.73 and 2.38 of 2.48), and the full relaxation leaves it a
little softer than the rigid rotation (coupling to the pendant bends), 0.04 and 0.19 more. What
the rigid route cannot see is small and is now on the record: the C-F stretch and bend modes add
0.06 (plain) and 0.11 (belt2) along the polar axis, and two chain-axis modes 0.13 and 0.30 along
the chain. The rigid-chain route's field step was checked to be in the linear regime (eps_0
lateral 4.3867 to 4.3858 from 2.5e-4 to 4e-3 V/A). The result files also carry the Born tensors
per site at the relaxed geometry (`"born"`: Z[i][a][b] = d mu_a / d u_ib in e, polyfind's frame),
for the comparison with DFPT's when Dipole's Gamma successor delivers them. Under belt3 the mean
diagonals by element are H (+0.25, +0.23, +0.26), F (-0.27, -0.26, -0.27), C (+0.02, +0.03,
+0.01), the sum rule met to 4e-13; under the pre-registered potential H (+0.19, +0.16, +0.19), F
(-0.20, -0.18, -0.19). These presets carry no charge flux, so their Born tensors are the fixed
bond-charge increments dressed by the induced dipoles and nothing more; the flux-fitted presets of
`docs/ELECTROMECHANICS.md` exist because that falls short of DFT's Born tensors, and the same will
be true here. The dielectric and piezoelectric columns of this document do not rest on the Born
tensors' size (P and the slope are anchored directly), but a Born comparison with Dipole's DFPT is
the place where the fixed-charge form will be seen to end.

So the prediction for Dipole's DFPT tensor of beta at the reference cell, in polyfind's frame
(their a is our y, their b our x): **eps_0 = 2.24 (polar), 4.57 (lateral), 2.74 (chain)**, against
eps_inf 2.13, 2.10, 2.44, with the lateral ionic part a libration at 46 cm^-1. The pre-registered
set says 2.27, 2.87, 2.65 for the same cell under its own potential. The lateral number is the
form's statement about the chains' rocking stiffness between hard contacts (5.14), and DFPT's
relaxed-ion lateral permittivity is the test of it: if DFT finds the libration twice as stiff
(eps_0 lateral near 3.3), the belt wall is too hard where the chains touch and too open between
contacts; if it finds 4.5, the wall is right and the lateral-coupling residual is elsewhere. Beta's
lowest modes (20 to 50 cm^-1 on belt2 against 31 to 64 on the plain chain) are the second thing a
Gamma phonon calculation would test; they are predictions, not anchors.

*Alpha, for the record* (`gamma_dielectric_pvdf-dft-valence-belt3_alpha.json`, under belt3): the
rigid-chain route gives eps_0 2.02, 2.35, 2.31 and the Gamma route 3.32, 2.49, 2.82, the polar
axis gaining 1.30 from one mode at 50 cm^-1, a torsional rocking of the TGTG' chain that a rigid
conformation cannot have, and the chain axis 0.19 and 0.10 from modes at 39 and 75 cm^-1. Two
caveats before it is compared with anything (5.22 traces the relaxation to the valence fit's soft C-C
stretch constant: with the constants floored at a physical value alpha's polar eps_0 is 2.93 at belt4
instead of 3.47): the all-atom relaxation moved atoms by 0.13 A and
lowered the cell's energy by 1.5 kcal/mol (the line-group refinement keeps every monomer
equivalent, the atoms need not), so the tensor is that of the fully relaxed cell, not the E1
reference's; and the pinned valence bend terms of the pendants are the valence fit's stiffnesses,
which alpha's soft modes lean on more than beta's do. For alpha the rigid-chain route is not a
substitute for the full one, and the sets' alpha eps_0 rows should be read as electronic plus
libration only.

### 5.17 Basins: the packing search's other minima, refined

Every E1/E2 set stands on one refined start: `mechanics.refined_reference` takes the first of the
distinct minima the packing screen returns (`pack()`, the exhaustive table screen under the
truncated Coulomb kernel, the best few polished and deduplicated) and refines it under the full
potential. The copolymer model cell landed in a different arrangement under each preset (5.10,
5.14), which asked what the other minima are. `examples/copolymer_basins.py` refines each of the
screen's top candidates under the same full potential, relaxes the cell, and lists them by their
full-potential energy with their cells and dipole directions (`--xyz` writes each as extended XYZ;
`docs/predictions/basins_<preset>_<system>.json`).

**Beta has two.** On both presets the screen finds two minima and the full potential cannot tell
them apart:

| preset | basin | a, b, gamma, c | dipoles from a | E per monomer, kcal/mol |
|---|---|---|---|---|
| `pvdf-dft-valence` | 0 (the E1 cell, the DFT setting) | 4.517, 8.520, 90.0, 2.602 | 180 deg (along a) | -4.5297 |
| `pvdf-dft-valence` | 1 (chains rotated 60 deg) | 4.823, 8.007, 85.8, 2.600 | 118 deg | -4.5290 |
| `pvdf-dft-valence-belt2` | 0 (the E1 cell, the DFT setting) | 4.745, 8.476, 90.0, 2.600 | 180 deg | -9.3627 |
| `pvdf-dft-valence-belt2` | 1 (chains rotated 60 deg) | 4.857, 8.284, 88.5, 2.600 | 119 deg | -9.3626 |

(energies with the backbone-angle strain; the screen's truncated-kernel energies put basin 0
first by 0.10 and 0.016, which is why every set sits in the DFT setting). Basin 1 is a real
alternative packing of the same chains: both setting angles turned by about 60 deg, the dipoles
pointing 119 deg from a instead of along it, the cell a little longer along a and shorter along b
and sheared 1.5 deg. The form puts 0.1 cal/mol per monomer between the two settings (0.7 on the
pre-registered potential). That is not a statement about PVDF: it says the potential's energy as a
function of the chains' orientation about their axes has no preference at 60 deg, and that angular
energy is exactly what the libration (5.16), the lateral coupling (5.14, 5.15) and the polar
response are made of. DFT will have a preference, and it is one relaxed cell: **the energy of beta
in the rotated setting against the standard one is the cheapest anchor on the chain-orientation
energetics there is**, and `basins_pvdf-dft-valence-belt2_beta_basin1.xyz` is its start geometry
(polyfind's frame, no belt sites; the standard setting beside it as `_basin0.xyz`, which is
Dipole's cell to 0.3/0.6/0.8%). It is asked for on the exchange.

**The copolymer has the same two, and the sets sit in different ones.** Under the pre-registered
potential the screen lists three minima for the 70/30 model cell: the arrangement its E1 set stands
on (dipoles 24 deg from a), one with the chains turned (84 deg) and a flipped, antiparallel-chain
one. Under belt2 it lists two, the arrangement that set stands on (85 deg) and the flipped one; the
24-deg arrangement is not among them. Refining each arrangement under the other potential
(`--start-from`) closes the question:

| arrangement (dipoles from a) | `pvdf-dft-valence`: a, b, gamma; E per monomer | `pvdf-dft-valence-belt2`: a, b, gamma; E per monomer |
|---|---|---|
| 24 deg (the pre-registered and stage-2 sets) | 4.792, 8.851, 92.0; -2.9382 | 4.796, 8.674, 90.1; -7.4199 (transplanted) |
| 85 deg (the belt2 set) | 5.104, 8.310, 92.2; -2.9382 (transplanted; the screen's own 84-deg basin -2.9380) | 4.959, 8.397, 92.1; -7.4199 |
| flipped chains, 90 deg | 4.947, 8.637, 90.0; -2.8220 | 4.944, 8.780, 90.0; -7.0969 |

(kcal/mol per monomer with the backbone-angle strain; the belt set of 5.10 sits at gamma 88 and 38
deg, a third point on what is evidently the same flat valley, not re-run here; under belt3 the
same pair is -7.2715 and -7.2716 with the antiparallel one at -6.8840,
`basins_pvdf-dft-valence-belt3_vdftrfe.json`). On every potential the two polar arrangements are
degenerate to 0.0001 kcal/mol per monomer and the antiparallel one is 0.12 to 0.39 above them. The screen's truncated-kernel energies, which separate the two by 0.1
to 0.3, decide which one a set stands on, and that is an accident of the screen, not a property of
the potential: every copolymer row that differed between the sets (5.10, 5.14: the cell, the
dipole tilt, the stiffness, the Q column) compared two arrangements the potentials themselves do
not distinguish. Two consequences. For the sets, a copolymer row has to name its arrangement, and
the arrangement-by-arrangement comparison is the one that means something: the basins files hold
both for both potentials, and the next E1/E2 runs will refine both and report both. For the
physics, the finding is the same as beta's: on every fitted potential the energy of the chains'
orientation about their axes is flat to 0.1 cal/mol, for PVDF and for the copolymer alike, and
nothing in the anchors in hand constrains it. The rotated-setting cell requested above is the
measurement that does.

**The arrangement-matched row.** `examples/e1_predictions.py --start-params` refines a given cell
and placement instead of the screen's choice (through `mechanics.refined_reference(start_params=)`),
so the belt2 copolymer E1 cell was run in the pre-registered set's arrangement
(`e1_polyfind_belt2_tilt24.json`, 822 s). Beside the two sets as they stand:

| copolymer TTTT | pre-registered (24 deg) | belt2 as published (85 deg) | belt2 in the 24-deg arrangement |
|---|---|---|---|
| a, b, gamma, c | 4.793, 8.852, 92.0, 25.99 | 4.959, 8.397, 92.1, 25.89 | 4.797, 8.674, 90.1, 25.89 |
| P (x, y), \|P\| C/m^2 | (0.100, 0.044) 0.109 | (0.015, 0.166) 0.166 | (0.152, 0.068) 0.166 |
| eps_0 relaxed-ion (x, y, z) | 2.070, 2.371, 2.379 | 2.769, 2.140, 2.417 | 2.137, 2.773, 2.417 |
| C_11, C_22, C_12, C_66 GPa | 16.5, 17.0, 1.2, 1.7 | 12.5, 21.1, 8.8, 9.7 | 21.1, 22.3, 3.9, 4.8 |
| proper e_11, e_12 C/m^2; d_11 pm/V | -0.164, -0.044; -10.2 | -0.052, -0.097; -2.7 | -0.244, -0.085; -11.5 |
| intrinsic Q_11, Q_22 m^4/C^2 | 3.06, 1.62 | -3.27, -0.78 | 3.43, 1.49 |
| d eps_yy / d eta_yy | 2.1 | -1.9 | 3.1 |

In the same arrangement the two potentials tell one story: the libration sits on the lateral axis
(eps_0 y 2.77 against 2.37, the belt form's softer rocking again), the piezoelectric column
follows the polar component of P (e_11 -0.244 against -0.164, with P_x 0.152 against 0.100), the
intrinsic Q is the same within 15%, and the stiffness rises as beta's did. Every sign flip and
factor-of-three difference the published belt2 row showed against the pre-registered one (5.14)
was the 85-deg arrangement's frame, not the potential. So the copolymer rows of the sets are
usable once the arrangement is named, and the cautions of 5.10 and 5.14 reduce to one: the
copolymer's two polar arrangements are degenerate on every fitted potential, the screen picks one,
and the next sets will carry both. The second thing this row says is physical: in the 85-deg
arrangement the copolymer's E2 margins put TTTT within 0.06 kcal/mol of TGTG' (5.14), and with
the 24-deg arrangement as a second TTTT minimum at the same energy the all-trans phase of the
copolymer is a flat valley of arrangements rather than a single cell, which is what a relaxor-like
copolymer should look like on a potential and a hint that the TrFE content's effect on the
margins is already in the form.

### 5.18 Stage 7 by least squares: the method, and a better minimum (`pvdf-dft-valence-belt3`)

Every stage so far was fitted with bounded Powell, 240 to 320 sequential evaluations per run. The
objective is a sum of squares, so `crystalfit.fit(method="trf")` now fits it as one: the residual
vector is the z of every anchor row plus sqrt(RIDGE) times each ridge component
(`ridge_parts`, so its sum of squares is the objective), the Jacobian is a forward difference of
that vector with per-variable steps (`FD_STEPS`, the same steps as 5.15's analysis), and
`scipy.optimize.least_squares` (trust-region reflective, the variables' bounds) takes it from
there. One iteration costs one evaluation per free variable plus one or two steps, and those
twelve evaluations are independent, which Powell's line searches never are: on a machine with
the cores they run at the cost of one.

From the stage-6 point on the PBE-D3 belt geometry, where Powell's 281 evaluations (78 min) ended
at 11.65 (5.13), the least-squares run passes 11.65 at its 144th evaluation, stops on its
evaluation budget at 9.36 after 235 (31 min, `crystal_stage_v7_trf.json`), and a continuation from
there converges at **9.11** (ridge 4.09) in 53 more (7 min, `crystal_stage_v7_trf2.json`): 288
evaluations and 62 minutes for a lower minimum than 281 and 78 reached.

| quantity | Powell point (5.13, belt2) | least-squares point (belt3) | DFT |
|---|---|---|---|
| beta a, b, c | 4.7445, 8.4756, 2.6003 (-0.26, -0.61, 0.75) | 4.7525, 8.4902, 2.5984 (-0.09, -0.44, 0.68) | 4.757, 8.528, 2.581 |
| alpha a, b, c | 4.9264, 9.7073, 4.6815 (-0.38, 0.74, 0.46) | 4.9422, 9.6819, 4.6757 (-0.06, 0.48, 0.34) | 4.945, 9.636, 4.660 |
| beta Berry P | 0.1928 (0.01) | 0.1928 (0.01) | 0.1927 |
| E(alpha) - E(beta) | -0.726 (-0.01) | -0.736 (-0.21) | -0.726 |
| C_11, C_22, C_12 | 23.46, 24.19, 4.99 (0.30, -0.14, -1.51) | 24.17, 24.78, 5.34 (0.61, 0.11, -0.92) | 22.77, 24.52, 5.88 |
| d|P|/d eta (polar) | -0.300 (0.84) | -0.300 (0.84) | -0.326 |
| beta backbone angle | 116.00 (1.37) | 115.87 (1.24) | 114.63 |
| alpha angle at CH2, CF2 | 118.05, 116.94 (-0.18, 0.07) | 117.83, 116.60 (-0.40, -0.27) | 118.23, 116.87 |
| C, H, F diameters; F-H | 3.465, 2.572, 3.369; 3.004 | 3.615, 2.603, 3.418; 2.998 | |
| restraint, kcal/mol/rad^2 | 178 | 228 | |
| ridge, objective | 4.45, 11.65 | 4.09, 9.11 | |

The two points differ most in the carbon diameter, which the least-squares point leaves at the
valence value (3.615 A, x_C 1.013) where Powell had taken it down 4% (3.465), and in the
restraint (228 against 178); with those, C_12 comes from 4.99 to 5.34 GPa (z -1.51 to -0.92),
beta's lateral axis from -0.61% to -0.44%, both cells' polar axes onto the anchors, and the
alpha-beta gap moves 0.01 kcal/mol (z -0.21), not the 2.3 tolerances the linear model at Powell's
point charged for it (5.15). Every residual is within 1.24 tolerances, the largest being beta's
backbone angle (and so its chain axis), as before. Two corrections follow. The Powell point of
5.13 was a local minimum, and the statement in 5.15 that the lateral coupling could only be
bought with the gap holds for the neighbourhood of that point: at the least-squares point the
residual that the form could not reach is 0.6 of a tolerance smaller, and the Jacobian there is
the one to judge the form by (`crystal_stage_v7_trf2_sensitivity.json`, the same 25 evaluations
and analysis as 5.15). It says something different from 5.15's. The singular values are much as
before (six constrained directions of twelve), but the trade is now benign: the linear model
puts C_12 and beta b on their anchors with a step of half a prior width (x_C +0.54, x_F +0.63,
D_F -0.72, the restraint +0.52 widths) and only the beta backbone angle at its 1.0 limit, and the
real evaluation at that point confirms it almost row by row (C_12 z -0.00, b -0.05, the gap -0.19,
P +0.60, beta a +0.61, the angle +1.00): sum of squares of the anchors 5.3 against 5.0 at the
minimum, ridge 7.51 against 4.09, objective 12.83. So at the least-squares point closing the
lateral coupling buys 1.0 in the anchors' sum of squares (C_12 and b) and costs 0.8 in the other
rows (P and beta's polar axis to +0.6, the alpha angles to -0.8) plus 3.4 of ridge: the residual is
held by the other anchors and the prior together, each about as much as the other. A probe with
the ridge at half weight (`--ridge 0.5`, `crystal_stage_v7_trf_ridge05.json`, from the belt3 point)
says the same from the other side: the fit does not move, 7.06 to 7.00 in 27 evaluations with the
restraint 228 to 238 and every other variable where it was, C_12 still at -0.92. The prior alone is
not the lock, and weakening it would not close the coupling; a term that moves C_12 without P and
the polar axis would, and the chain-orientation anchor requested in 5.17 is the datum to build it
against. The point is the preset `pvdf-dft-valence-belt3`; its prediction sets, Gamma-route
dielectric tensor and cross-check are below (5.19). The belt2 preset and its sets stay as the
record of 5.13-5.17; `--method trf` is the fitter to use from here, and the Powell runs stand as
what they were.

### 5.19 The belt3 preset as a prediction set

`examples/e1_predictions.py --preset pvdf-dft-valence-belt3 --variant belt --geometry pbe` and the
E2 script likewise (`e1_polyfind_belt3.json`, `e2_polyfind_belt3.json`; the copolymer cell also in
the pre-registered set's arrangement, `e1_polyfind_belt3_tilt24.json`, by `--start-params`), the
refinement under the full potential, the belt sites out of the electrostatics. Beside belt2 and
the DFT in hand:

| quantity | DFT | belt2 (5.14) | belt3 |
|---|---|---|---|
| beta a, b, c, A | 4.757, 8.528, 2.581 | 4.745, 8.476, 2.600 | 4.753, 8.490, 2.598 |
| beta polarization, C/m^2 | 0.193 (Berry) | 0.193 | 0.193 |
| beta eps_inf (x, y, z) | 2.235, 2.253, 2.600 (DFPT, protocol evidence) | 2.122, 2.095, 2.448 | 2.119, 2.091, 2.443 |
| beta eps_0, rigid-chain route (x, y, z) | - | 2.122, 4.386, 2.448 | 2.119, 4.339, 2.443 |
| beta eps_0, Gamma route (5.16) | - | 2.236, 4.573, 2.735 | 2.232, 4.489, 2.731 |
| beta C_11, C_22, C_12, C_66, GPa | 22.8, 24.5, 5.9, - | 23.5, 24.2, 5.0, 5.5 | 24.2, 24.8, 5.3, 5.9 |
| beta d\|P\|/d eta (polar), C/m^2 | -0.326 (Berry pair) | -0.300 | -0.300 |
| beta proper e_11, C/m^2; d_11, pm/V | - | 0.300; 12.8 | 0.300; 12.5 |
| beta intrinsic Q_11, Q_22, m^4/C^2 | - | 5.58, 11.56 | 5.28, 11.00 |
| beta d eps_yy / d eta_yy | - | 48.6 | 45.7 |
| alpha (polar TGTG') a, b, c, A | 4.945, 9.636, 4.660 (the antipolar anchor) | 4.926, 9.707, 4.682 | 4.942, 9.682, 4.676 |
| alpha polarization, C/m^2; C_11, C_22, C_12, C_66 | - | 0.099; 21.1, 15.0, 5.2, 6.2 | 0.101; 22.3, 15.4, 5.3, 6.4 |
| copolymer, the screen's arrangement (dipoles 85 deg from a): a, b, gamma, c per 10 | - | 4.959, 8.397, 92.1, 25.89 | 4.986, 8.423, 91.8, 25.88 |
| copolymer (85 deg) polarization; eps_0; C_11, C_22, C_12, C_66 | - | 0.166; 2.769, 2.140, 2.417; 12.5, 21.1, 8.8, 9.7 | 0.165; 2.707, 2.125, 2.404; 13.3, 21.4, 8.9, 9.9 |
| copolymer, the pre-registered arrangement (24 deg): a, b, gamma, c per 10 | - | 4.797, 8.674, 90.1, 25.89 | 4.824, 8.700, 90.3, 25.88 |
| copolymer (24 deg) P (x, y); eps_0; C_11, C_22, C_12, C_66 | - | (0.152, 0.068); 2.137, 2.773, 2.417; 21.1, 22.3, 3.9, 4.8 | (0.150, 0.068); 2.122, 2.711, 2.404; 21.3, 23.5, 3.8, 4.9 |
| copolymer (24 deg) proper e_11; d_11; Q_11, Q_22 | - | -0.244; -11.5; 3.43, 1.49 | -0.244; -11.3; 3.44, 1.27 |
| PE (control) a, b, gamma; eps_inf (x, y, z) | - | 3.971, 6.887, 90.0; 3.038, 3.031, 3.094 | 4.048, 6.992, 90.0; 2.930, 2.926, 3.019 |
| PE C_11, C_22, C_12, C_66, GPa | - | 31.7, 33.8, 1.0, 4.9 | 29.6, 32.5, 1.7, 4.4 |

E2, margins in kcal/mol per monomer with the backbone-angle strain, relative to TGTG':

| | DFT | belt2 | belt3 |
|---|---|---|---|
| PVDF TTTT | +0.726 | +0.726 | +0.736 |
| PVDF T3GT3G' | - | +0.333 | +0.315 |
| copolymer TTTT | - | +0.059 | -0.036 |
| copolymer T3GT3G' | - | +0.030 | -0.057 |

**Against Dipole's records** (`examples/crosscheck_e1_e2.py --polyfind docs/predictions/e1_polyfind_belt3.json --e2 docs/predictions/e2_polyfind_belt3.json`):
beta's cell -0.1%, -0.4%, +0.7%; C_11 at 1.06 of the six-cell value, C_22 at 1.01, C_12 at 0.91
(belt2: 1.03, 0.99, 0.85); the Berry polarization at 1.00; the polar Berry slope at 0.92 and 0.92
relative to P; the electronic permittivity at 0.95, 0.93, 0.94. Every anchored quantity is where
the fit put it; the fit put them closer than before, which is 5.18.

**What changed, and what did not.** Almost nothing in the predictions: beta's libration is the
same mode at the same softness (eps_0 lateral 4.49 by the Gamma route against 4.57; the 47 cm^-1
mode carrying 2.30 of 2.40), its piezoelectric column the same (e_11 0.300, d_11 12.5), its
intrinsic Q 5% smaller; alpha is the same cell to 0.3%; PE is squeezed 4% less in volume with a
lower electronic permittivity (2.93 against 3.04) and the same stiffness, so the PE caution of
5.14 stands with a slightly smaller number. The copolymer under belt3 sits in the 85-deg
arrangement again by the screen's choice and reads like belt2's in it; in the pre-registered
arrangement it reads like the pre-registered set's (5.17), the two presets agreeing row by row in both (e_11 -0.244 in both, eps_0 lateral 2.71 against 2.77, Q_11 3.44 against 3.43). The one physical
change is E2's copolymer: under belt3 the T3GT3G' conformation is the lowest of the three model
cells, TTTT 0.021 above it and TGTG' 0.057, where belt2 had TGTG' lowest by 0.03 and 0.06; three
conformations within 0.06 kcal/mol per monomer on both presets, i.e. the copolymer's conformational
energy is flat at the level the fit resolves, which is the relaxor-like statement of 5.17 again and
the reason the copolymer rows of E2 should not be read as a ranking. The prediction posted for the
Gamma successor (5.16) stands within 2%.

### 5.20 What a new term would have to be: the radial candidates screened

5.14 named three candidates for the lateral coupling, and 5.18 said the residual is held by the
other anchors and the prior together, so a term helps only if it moves C_12 and beta's lateral axis
in a way the twelve fitted variables cannot already produce by moving. That is a question the
Jacobian answers before anything is fitted: perturb the candidate at the belt3 point, take the
change of the residual vector per step, and project it onto the column space of the Jacobian of
the twelve variables (`examples/stage8_screen.py`, `docs/predictions/stage8_screen_belt3.json`).
The part outside that span is the new direction; its value is how much of the residual the twelve
variables cannot reach it would take out.

| candidate (perturbation) | residual change per step | outside the span of the twelve | gain within one step, of 1.62 | what it moves most |
|---|---|---|---|---|
| belt distance d, 0.7 +- 0.1 A | 3.15 | 0.09 | 0.20 (the full 1.28 would need 12 steps, 1.2 A) | beta a +1.5, alpha b +1.6, gap -1.6, C_11 -1.0, C_12 +1.0 |
| belt site diameter, 2.6 +- 0.2 A | 4.42 | 0.00 | 0.002 | gap -2.8, beta a +2.0, alpha b +1.7, C_11 -1.3, b +0.9 |
| F-F pair diameter x 1.02 | 3.62 | 0.00 | 0.007 | gap -3.6, C_12 +0.3, alpha b +0.4 |
| F-F pair depth x 1.5 | 5.07 | 0.01 | 0.006 | gap -5.0, C_12 +0.6, alpha a -0.5 |
| H-H pair diameter x 1.02 | 0.03 | 0.00 | 0.001 | nothing |
| C-C pair diameter x 1.02 | 0.46 | 0.00 | 0.009 | gap +0.4 |
| C-F pair diameter x 1.02 | 0.56 | 0.01 | 0.016 | gap +0.4, alpha a +0.3, b +0.15 |
| C-H pair diameter x 1.02 | 0.20 | 0.00 | 0.002 | nothing |
| fluorine charge offset, 0.177 +- 0.03 A (an electrostatic shape the form already carries, fixed) | 0.51 | 0.00 | 0.000 | P +0.27, gap -0.23, C_12 +0.20 |
| polarizability scale x 1.05 (the induced-dipole response, fixed at the VDS98 values) | 0.55 | 0.15 | 0.22 (the full 0.65 would need the polarizabilities 23% smaller) | P +0.19, gap -0.39, slope -0.28 |

(the linear model's attainable minimum of the anchors' sum of squares with the twelve variables,
bounds ignored, is 1.62 of the 5.02 at the point; "gain within one step" is how much of that a
candidate moved by one step, the twelve variables re-optimised, would remove). Every radial contact
term, the belt geometry included, lies in the span of the twelve variables to within 0 to 3 percent
of its own length: the belt site's diameter is degenerate with its depth (a repulsive wall has one
effective radius), its distance from the fluorine with the fluorine and pair diameters, and a
pair-specific F-F, C-C, C-F, C-H or H-H term with the element diameters that already set every
contact. None of them is a new direction worth having: the one with any gain at all, the belt
distance, would need to move 1.2 A to deliver it. Of the two electrostatic quantities the form
carries fixed, the fluorine charge's offset from its nucleus is in the span as well (it is one more
way of setting P), while the polarizability scale is the only probed quantity with a direction of
its own, 0.15 of its 0.55 outside the span: it moves P, the gap and the polar slope in a combination
the twelve variables do not have, and would take 0.22 of the 1.62 within a 5% change. It is not a
lateral-coupling term, and the direction that gains points to smaller polarizabilities, which the
electronic permittivity, already 0.93 to 0.95 of DFPT's on every set, forbids; it stays fixed, and
is noted here as the one fixed quantity with a direction of its own. The anchors in hand (fifteen informative rows) fix six directions of the
twelve-variable space (5.15), so the problem is not a shortage of variables but that the lateral
coupling, P, the polar axis and the gap are moved by the same combinations of them: every radial
term pushes on the same combinations. What would not is a term with a shape rather than a radius:
the anisotropy the belt gave the fluorine, given to the hydrogen (sites on the CH2 group), or a
quadrupole on the CF2 group that changes how the chains' electrostatics couple the two lateral
axes without changing P. Both are a stage of work, and both should wait for the rotated-setting
cell (5.17): one more anchor row on the chain-orientation energy may change which direction is
missing, and it is the measurement of exactly the energy these terms would shape.

### 5.21 Stage 8, as designed and not yet built

What 5.20 leaves is a term with a shape. Two are concrete enough to name, with the mechanism each
would ride on and the test it has to pass before it is fitted.

**A hydrogen belt.** The anisotropy the belt sites gave the fluorine (5.7) given to the CH2 group:
massless, chargeless sites either side of each hydrogen in the HCH plane, a repulsion-only wall
against fluorine. It rides on the belt machinery as it stands (pendant sites, a pseudo-element
outside every bonded and electrostatic sum, a pair table with one non-zero row), but needs a
second pseudo-element, since the fluorine sites' wall is X-H and a hydrogen site's must be Y-F;
the packer's and force field's special cases for "X" become a set. Two variables (depth, as for
the fluorine belt; the sites' offset fixed at the probe's value). The doubt is written in 5.20:
the F...H contact is the same contact seen from the other side, and the screen may find this
direction in the span as well; it is the first thing to check, with `examples/stage8_screen.py`
on a tiny depth, before anything else is written.

**A CF2 quadrupole.** The electrostatics of the chain are bond-charge increments with the
fluorine's charge displaced along its bond (`charge_offsets`, 0.177 A). A quadrupole on the CF2
group, two further charge sites of opposite sign placed symmetrically about the C-F bond so the
group's dipole does not change, alters how neighbouring chains' charge distributions couple the
two lateral axes while leaving P, the Berry anchor and the slope where they are to first order.
It rides on the pendant-site and charge machinery, but the sites carry charge, so unlike the belt
they enter the Ewald sum, the induced-dipole solve and the dipole; the inert-site exclusion of
16.6 does not apply to them and the cost is that of two more charged atoms per CF2. One variable
(the quadrupole strength at a fixed geometry). It is the candidate the screen most likely passes,
because it is the only one of the three whose first-order effect on the anchors is not a contact.

**The order.** Neither is built before Dipole's rotated-setting cell lands (5.17): it adds the
first anchor row on the chain-orientation energy, which is the energy both terms shape, and with
it the Jacobian's picture of which direction is missing may change. When it does, the sequence is
the screen (two evaluations per candidate), then the stage-8 fit by least squares from the belt3
point with the surviving term as the fourteenth variable, then the sets and the Gamma-route tensor
as in 5.19.

One thing came before either (5.22): the refinement's rigid CH2 and CF2 groups, which cannot tilt
as DFT's do, and the 0.23 kcal/mol per monomer that costs the alpha-beta gap row. The stage-8 fit
starts from the corrected point, `pvdf-dft-valence-belt4`, whose lateral coupling is 1.3 tolerances
short, not from the belt3 one.

### 5.22 The pendant groups: what the rigid refinement cannot do, and what it costs the gap anchor

**The observation.** Dipole's PBE-D3 start geometry for the alpha successor
(`pvdf_alpha_ref_fresh_20261107_v1`, the TGTG' antipolar start whose native relaxation is pending;
read from its `state.json` by `examples/dft_cell_geometry.py`) has CH2 and CF2 groups that are not
where polyfind's chain builder puts them. The builder places both pendants of a backbone atom at
`cos(a/2) u +- sin(a/2) w`, `u` the outward bisector of the backbone angle and `w` the normal of
the backbone plane (`polyfind.chain.pendant_frames`); the refinement moves the cell, the setting,
the torsions and the backbone angles and never a pendant relative to its carbon. A real group can
differ from that placement by three small rigid rotations, measured against the local frame
`(u, v, w)` with `v` along the chain: *rock* (about `w`: the group's bisector tilts toward a
backbone neighbour), *wag* (about `v`: the bisector tilts out of the backbone plane) and *twist*
(about `u`: the X-X line turns toward the chain direction). The reader now reports all three, and
the plain-angle view of the same thing, the C-C-X asymmetries (rock: one X against the two
neighbours; wag: the two X against one neighbour; twist: both):

| cell | X-C-X (FF, HH) | CF2 rock, wag, twist | CH2 rock, wag, twist | largest C-C-X asymmetry |
|---|---|---|---|---|
| Dipole's alpha start geometry | 104.70, 107.17 | 2.29, 2.86, 1.29 | 1.01, 0.66, 0.88 | 4.35 deg |
| Dipole's beta reference (`pvdf_beta_ref_v1`) | 105.12, 107.64 | 0, 0, 0 | 0, 0, 0 | 0.00 |
| polyfind, belt3 fit path, alpha | 105.12, 107.64 | 0, 0, 0 | 0, 0, 0 | 0.00 |
| polyfind, belt3 fit path, beta | 105.12, 107.64 | 0, 0, 0 | 0, 0, 0 | 0.00 |

Beta's zeros are symmetry: the planar zigzag has a mirror through every group, so the symmetric
placement is exact and the refinement loses nothing by freezing it. Alpha's are not: in TGTG' the
two sides of every backbone atom are inequivalent (a trans neighbour one way, a gauche one the
other) and nothing forbids the group to turn. The DFT geometry turns it, by two to three degrees
on the CF2 groups, and the fit path cannot follow. (The X-C-X angles differ too, 104.70 against
beta's 105.12 at the CF2 and 107.17 against 107.64 at the CH2; the chain geometry is beta's, 5.11.)

**What the fitted potential wants.** `polyfind.pendants.tilt_relaxation` lets every group of a
refined reference turn rigidly about its carbon, three angles per group, with the carbons, the
cell, the bond lengths and the X-C-X angles fixed, the belt sites carried along as the affine
combinations of their CF2 atoms that place them (the Gamma route's construction, 5.16, now shared
in `polyfind.pendants`), and the valence bend terms pinned at the reference geometry so the C-C-X
bends resist with the valence fit's stiffness (78.6 kcal/mol/rad^2 for C-C-F, 60.4 for C-C-H, in
the kernel's 1/2 k convention) and the reference is a stationary point of the valence part. A
dozen L-BFGS iterations on the kernel's analytic forces, about a second per polymorph. Under
`pvdf-dft-valence-belt3`:

| polymorph | tilt-only (rigid groups) | pendant-only (H, F free) | all-atom (translation-free) |
|---|---|---|---|
| alpha, kcal/mol per monomer | -0.226 | -0.245 | -0.368 |
| alpha, CF2 rock, wag, twist (deg) | 1.17, 2.65, 0.73 | 1.19, 2.61, 0.75 | 1.37, 2.38, 0.91 |
| alpha, CH2 rock, wag, twist (deg) | 0.67, 2.48, 0.50 | 0.69, 2.48, 0.51 | 0.71, 2.36, 0.56 |
| beta, kcal/mol per monomer | 0.000 | -0.007 | -0.032 |
| E(alpha) - E(beta), from -0.736 rigid | -0.963 | -0.975 | -1.072 |

The tilt-only column is the one that means something: the pendant-only and all-atom relaxations
(the Gamma route's machinery, 5.16, with the carbons fixed or free) add bond stretching at the
valence fit's stiffnesses, and the all-atom one stretches every C-C bond from 1.533 to 1.553 A and
turns the gauche torsion from 60.7 to 56.0 deg against a C-C stretch constant that the valence fit
left at its bound (100 kcal/mol/A^2, five times softer than a real C-C bond): that is a weakness
of the pinned valence, not a property of the crystal, and it is also the origin of the 0.13-0.18 A
caveat on alpha's Gamma-route tensor in 5.16. The rigid tilts carry none of it. Their pattern is
DFT's: the CF2 groups wag most and rock second (2.65 and 1.17 against the start geometry's 2.86
and 2.29), the twist is smallest, and the signs follow the glide of TGTG' (the four CF2 groups
turn by the same magnitudes, two one way and two the other). The CH2 wag is the one component the
potential overdoes (2.48 against 0.66). Beta's groups do not move: the symmetric placement is a
stationary point, as it must be.

**What it costs the gap anchor.** The fit matched E(alpha) - E(beta) = -0.726 (DFT, relaxed cells
on both sides) with -0.736 from rigid-pendant refinements. With the tilts relaxed the same
parameters give -0.963, 4.7 tolerances off: the belt3 point carries 0.23 kcal/mol per monomer of
pendant strain in alpha that DFT's alpha does not have, and the other terms were fitted to make up
for it. The linear preview on the 5.18 Jacobian (`crystal_stage_v7_trf2_sensitivity.json`, the gap
row shifted by -4.53 and the least-squares problem re-solved with the ridge) puts the objective at
31.5 at the belt3 point and 10.7 after the re-solve, paid mainly by the fluorine well (D_F -8%),
the H-F pair depth (+7%), both charge increments (+0.019 and +0.014 e) and the restraint
(+0.04 dex), with the lateral coupling the row that gives (C_12 from -0.92 to -1.23). So the
correction moves the fit, and it moves the quantities 5.15 and 5.18 found held by the other anchors
and the prior together.

**What was done.** `examples/crystal_stage.py --pendant-tilts` evaluates the gap row with both
polymorphs' tilt relaxations added to their refined energies (`crystalfit.evaluate(...,
pendant_tilts=True)`; the relaxation energies and tilts are recorded in the result file under
`pendant_tilts`). An evaluation costs 15 s against 8 s without (two reference rebuilds and
relaxations). The refit by bounded least squares from the belt3 point is `crystal_stage_v8_tilts.json`.

**The refit.** 196 evaluations in 31 minutes (10 s each; the two reference rebuilds and tilt
relaxations add about 2 s to an evaluation): objective 11.11 (ridge 5.21), against 9.11 (4.09)
for the uncorrected problem at the belt3 point and 31.5 for the corrected one there. Registered as
`pvdf-dft-valence-belt4`, the current preset.

| quantity | DFT | belt3 (z) | belt4 (z) |
|---|---|---|---|
| beta a, b, c (A) | 4.757, 8.528, 2.581 | 4.753 (-0.09), 8.490 (-0.44), 2.598 (+0.68) | 4.755 (-0.04), 8.478 (-0.59), 2.596 (+0.58) |
| alpha a, b, c (A) | 4.945, 9.636, 4.660 | 4.942 (-0.06), 9.682 (+0.48), 4.676 (+0.34) | 4.957 (+0.24), 9.646 (+0.11), 4.680 (+0.43) |
| beta P (C/m^2) | 0.1927 | 0.1928 (+0.01) | 0.1927 (-0.00) |
| E(alpha) - E(beta), kcal/mol per monomer | -0.726 | -0.736 rigid (-0.21) | -0.744 with the tilts (-0.36) |
| beta C_11, C_22, C_12 (GPa) | 22.77, 24.52, 5.88 | 24.17 (+0.61), 24.78 (+0.11), 5.34 (-0.92) | 24.21 (+0.63), 25.07 (+0.22), 5.12 (-1.29) |
| beta d|P|/d eta (polar) | -0.3255 | -0.3004 (+0.84) | -0.2998 (+0.86) |
| backbone angles: beta, alpha CH2, alpha CF2 (deg) | 114.63, 118.23, 116.87 | 115.87 (+1.24), 117.83 (-0.40), 116.60 (-0.27) | 115.70 (+1.07), 117.56 (-0.67), 116.54 (-0.33) |
| alpha's tilt relaxation, kcal/mol per monomer | | -0.226 (not in the fit) | -0.205 (in the fit) |
| alpha CF2 rock, wag, twist; CH2 (deg) | 2.29, 2.86, 1.29; 1.01, 0.66, 0.88 (start geometry) | 1.17, 2.65, 0.73; 0.67, 2.48, 0.50 | 1.09, 2.47, 0.73; 0.72, 2.46, 0.53 |
| objective (ridge) | | 9.11 (4.09) | 11.11 (5.21) |

The gap is back on the anchor with the relaxation inside the fit, and what paid for it is what the
linear preview said would: the fluorine well 12% shallower (D_F 0.923 to 0.815 of the valence
value, which is where the ridge's extra 1.1 comes from), the H-F pair 8% deeper, the fluorine
diameter 1.7% smaller, the carbon and hydrogen wells 2-3% deeper, the charge increments by 0.002
and 0.004 e, the restraint 228 to 233 kcal/mol/rad^2. In the anchors the lateral coupling pays
(C_12 from 0.91 to 0.87 of DFT's, z -0.92 to -1.29, the largest residual now; the lateral axis
-0.4 to -0.6%) and the alpha CH2 angle a little; alpha's lateral axis (+0.5 to +0.1%), beta's polar
axis and the beta backbone angle (1.24 to 1.07) improve. So 5.15's reading sharpens: at the belt3
point the gap anchor was helping the lateral coupling through alpha's pendant strain, and with that
strain released the coupling is 1.3 tolerances short on its own. The stage-8 shape term (5.21)
starts from this point. The tilts at the new parameters are within 0.2 deg of the belt3 ones (CF2
wag 2.47, rock 1.09, twist 0.73), so the prediction for Dipole's native alpha positions is
unchanged by the refit. The E2 set with every conformation's groups relaxed
(`e2_polyfind_belt4.json`) puts TTTT +0.744 and T3GT3G' +0.402 above TGTG' (PVDF; T3GT3G'
releases 0.040), and reorders the VDF/TrFE 70/30 model's conformations: the copolymer's TGTG'
releases 0.357 kcal/mol per monomer (its CF2 groups wag 3.1 deg, its CHF groups wag 2.3 and twist
1.3), its T3GT3G' 0.122 and its TTTT 0.040 (the CHF groups break the zigzag's mirror), so TGTG' is
the ground conformation by 0.13 to 0.15 where the rigid belt3 set had T3GT3G' by 0.02 to 0.06; a
copolymer E2 record from Dipole would arbitrate, and the pre-registered set stands as the primary
one. The Gamma-route tensors at belt4 are in `gamma_dielectric_pvdf-dft-valence-belt4.json`: beta
eps_0 2.23 polar, 4.61 lateral, 2.73 chain (belt3 4.49 lateral), eps_inf 2.12, 2.09, 2.44, the
lateral libration at 46 cm^-1; alpha 3.47, 2.51, 2.81 with a 0.12 A all-atom relaxation.

**The soft C-C constant, quantified.** The all-atom routes (the Gamma route of 5.16, the pendant-only
and all-atom columns above) pin the valence fit's stretch constants, and the C-C one sits at that
fit's bound, 100 kcal/mol/A^2, five times softer than a real C-C bond. `gamma_dielectric.py
--k-bond-min 540` floors every stretch constant at a physical value
(`gamma_dielectric_pvdf-dft-valence-belt4_kbond540.json`): beta's tensor does not move (eps_0 2.23,
4.60, 2.73 against 2.23, 4.61, 2.73; its relaxation is 0.004 A either way), so the prediction posted
for Dipole's DFPT stands whatever the constant; alpha's polar component falls from 3.47 to 2.93
(its ionic part 1.45 to 0.90, the relaxation 0.12 to 0.08 A), the lateral and chain components by
0.03 and 0.11. The alpha caveat of 5.16 is therefore the soft constant itself: alpha's Gamma-route
polar eps_0 is 2.9 to 3.5 depending on it, and the floored value is the one to quote until the
stretch terms are refitted. The rigid tilts of 5.22 carry no bond stretching and are unaffected.

**The Jacobian at belt4** (`crystal_stage_v8_tilts_sensitivity.json`, taken with the same
correction, 25 evaluations) reads as 5.15 and 5.18 did: six of twelve directions are fixed by
the anchors within a prior width, the gap row is the stiffest (8.2 and 6.3 tolerances per width
of x_HF and D_HF, 7.5 of the belt wall), and the linear model closes C_12 and the lateral axis
together only by driving beta's polar axis, alpha's polar axis, C_11 and the beta backbone angle to
their full tolerance (step 1 to 2 widths in x_C, x_H, D_F). The lateral coupling is the form's
residual at belt4 as it was at belt3, and the stage-8 shape term is still what closes it. The 5.20
screen rerun against this Jacobian (`stage8_screen_belt4.json`, 20 evaluations) says the same of
the radial candidates: every pair diameter and depth, the belt diameter and the fluorine charge
offset lie in the span (outside part at most 0.006), the belt distance 0.056 with a gain of 0.06 per
step, and the polarizability scale is again the one fixed quantity with a direction of its own
(0.14) and again points to smaller polarizabilities, which eps_inf forbids.

**What it does not do.** The cell and the backbone do not re-relax with the tilted pendants, so
the corrected gap is complete to first order only; the alpha cell rows and the alpha angle rows are
still the rigid refinement's. Every alpha-type energy in the prediction sets carries the same
bias: E2's PVDF TTTT margin is the gap and takes the full 0.23, and the TGTG' and T3GT3G' rows of
every chemistry in `docs/E1_E2_PREDICTIONS.md` are rigid-pendant numbers. `examples/e2_predictions.py
--pendant-tilts` relaxes every conformation the same way; at the belt3 parameters PVDF's TGTG'
releases 0.226, T3GT3G' 0.052 (three of its four backbone atoms sit between trans bonds and barely
turn) and TTTT nothing, so the margins become TTTT +0.963 and T3GT3G' +0.489 above TGTG' against
the rigid set's +0.736 and +0.315 (the belt3 set stands as labelled; the refit's set carries the
correction). Beta's rows, the Berry anchors, the stiffnesses and the tensors are untouched, by
symmetry. The complete fix is the tilts as refinement variables, three per backbone atom of the
repeat, placed through the pendant frames' local coordinates (`polyfind.polymers.PendantAtom`
offsets, which the builder and `chainjac` already carry) with the valence C-C-X bends as their
restoring force; it is a stage-9 item, and the cheap consistent evaluation here is what the fit
can use now. Dipole's native alpha positions, when they land, test the pattern directly: the
prediction is CF2 wag 2.5, rock 1.1, twist 0.7 deg and CH2 wag 2.5 (belt4; belt3 within 0.2 deg),
against the start geometry's 2.9, 2.3, 1.3 and 0.7.

## 6. What it costs, and what it does not touch

(Stage 1; the same applies to stages 2 to 4, whose extra terms are likewise outside the
valence range; the belt sites have no interaction at all within their own chain, so the conformer
energetics are untouched by construction, `tests/test_belt_sites.py`.)

- **The conformer trainset.** The torsion, bond and angle terms are untouched, so the 467
  conformer frames of `docs/VALENCE_FIT.md` see this preset only through the nonbonded terms
  that act beyond the valence range (1-4 and further). The frames are not on this machine
  (`POLYFIND_TRAINSET` is unset), so the held-out scores and the three acceptance tests of
  that document have not been re-run on this preset; the ridge is the only statement about the
  cost so far (a combined three widths), and re-scoring is the first thing to do where the
  trainset is mounted.
- **Everything measured on `pvdf-dft-valence`** (the E1/E2 pre-registration, the loss-model
  calibrations, the level-2 lattice) is unchanged: this preset is opt-in by name, like every
  fitted preset in `polyfind.fitting`.
- **What it does not fit.** The strained-cell stiffness, the dielectric tensors and the
  piezoelectric constants are predictions of this preset, not inputs: Dipole's six axial cells
  and DFPT records are the test, through `examples/crosscheck_e1_e2.py --polyfind
  docs/predictions/e1_polyfind_crystal.json` (and the later sets' files likewise; from stage 2
  on the stiffness, and from stage 5 on the Berry polarization and its polar slope, are anchors
  rather than predictions, as each stage's section says).

## 7. Reproduce

```
python examples/crystal_stage.py --maxfev 300 --out docs/predictions/crystal_stage_v1.json
python examples/crystal_stage.py --stage 2 --maxfev 320 --out docs/predictions/crystal_stage_v2.json
python examples/crystal_stage.py --stage 3 --maxfev 320 --out docs/predictions/crystal_stage_v3.json
python examples/crystal_stage.py --stage 4 --maxfev 320 --out docs/predictions/crystal_stage_v4.json
python examples/crystal_stage.py --stage 4 --x0 docs/predictions/crystal_stage_v4.json --maxfev 320 --out docs/predictions/crystal_stage_v5.json   # the Berry P anchor
python examples/crystal_stage.py --stage 6 --x0 docs/predictions/crystal_stage_v5.json --maxfev 300 --out docs/predictions/crystal_stage_v6.json   # plus the polar Berry slope
python examples/crystal_stage.py --stage 7 --geometry pbe --x0 docs/predictions/crystal_stage_v6.json --maxfev 280 --out docs/predictions/crystal_stage_v7.json   # plus the restraint's stiffness, Dipole's backbone angles as anchors
python examples/e1_predictions.py --preset pvdf-dft-valence-belt --variant belt --out docs/predictions/e1_polyfind_belt.json
python examples/e2_predictions.py --preset pvdf-dft-valence-belt --variant belt --out docs/predictions/e2_polyfind_belt.json
python examples/crosscheck_e1_e2.py --polyfind docs/predictions/e1_polyfind_belt.json --e2 docs/predictions/e2_polyfind_belt.json
python examples/e1_predictions.py --preset pvdf-dft-valence-belt2 --variant belt --geometry pbe --out docs/predictions/e1_polyfind_belt2.json
python examples/e2_predictions.py --preset pvdf-dft-valence-belt2 --variant belt --geometry pbe --out docs/predictions/e2_polyfind_belt2.json
python examples/crosscheck_e1_e2.py --polyfind docs/predictions/e1_polyfind_belt2.json --e2 docs/predictions/e2_polyfind_belt2.json
python examples/crystal_stage_sensitivity.py --fit docs/predictions/crystal_stage_v7.json --out docs/predictions/crystal_stage_v7_sensitivity.json   # 5.15: the Jacobian at the stage-7 point (25 evaluations)
python examples/crystal_stage_sensitivity.py --analyze docs/predictions/crystal_stage_v7_sensitivity.json --verify   # its linear analysis and the real evaluation of the lateral-coupling trade
python examples/gamma_dielectric.py --preset pvdf-dft-valence --systems beta   # 5.16: the relaxed-ion tensor by the Gamma route, plain chain
python examples/gamma_dielectric.py --preset pvdf-dft-valence-belt2 --variant belt --geometry pbe --systems beta   # 5.16: the same under belt2 (docs/predictions/gamma_dielectric_<preset>.json)
python examples/crystal_stage.py --stage 7 --geometry pbe --method trf --x0 docs/predictions/crystal_stage_v6.json --maxfev 40 --out docs/predictions/crystal_stage_v7_trf.json   # 5.18: stage 7 by bounded least squares (235 evaluations, 9.36)
python examples/crystal_stage.py --stage 7 --geometry pbe --method trf --x0 docs/predictions/crystal_stage_v7_trf.json --maxfev 20 --out docs/predictions/crystal_stage_v7_trf2.json   # its continuation to convergence (53 evaluations, 9.11): the preset pvdf-dft-valence-belt3
python examples/crystal_stage_sensitivity.py --fit docs/predictions/crystal_stage_v7_trf2.json --out docs/predictions/crystal_stage_v7_trf2_sensitivity.json
python examples/crystal_stage.py --stage 7 --geometry pbe --method trf --ridge 0.5 --x0 docs/predictions/crystal_stage_v7_trf2.json --maxfev 30 --out docs/predictions/crystal_stage_v7_trf_ridge05.json   # 5.18: the half-weight ridge probe (the fit does not move)
python examples/stage8_screen.py --fit docs/predictions/crystal_stage_v7_trf2.json --sensitivity docs/predictions/crystal_stage_v7_trf2_sensitivity.json --out docs/predictions/stage8_screen_belt3.json   # 5.20: candidate terms against the Jacobian
python examples/e1_predictions.py --preset pvdf-dft-valence-belt3 --variant belt --geometry pbe --out docs/predictions/e1_polyfind_belt3.json
python examples/e2_predictions.py --preset pvdf-dft-valence-belt3 --variant belt --geometry pbe --out docs/predictions/e2_polyfind_belt3.json
python examples/gamma_dielectric.py --preset pvdf-dft-valence-belt3 --variant belt --geometry pbe --systems beta
python examples/dft_cell_geometry.py <sarco>/materials/electrostriction/results/qe_equilibrium_v1/pvdf_alpha_ref_fresh_20261107_v1/state.json   # 5.22: the pendant tilts of Dipole's alpha start geometry (rock, wag, twist against the builder's placement)
python examples/crystal_stage.py --stage 7 --geometry pbe --method trf --pendant-tilts --x0 docs/predictions/crystal_stage_v7_trf2.json --maxfev 160 --out docs/predictions/crystal_stage_v8_tilts.json   # 5.22: the gap row with both polymorphs' pendant groups relaxed as rigid tilts
python examples/e2_predictions.py --preset pvdf-dft-valence-belt4 --variant belt --geometry pbe --pendant-tilts --out docs/predictions/e2_polyfind_belt4.json   # 5.22: the E2 set with every conformation's pendant groups relaxed
python examples/gamma_dielectric.py --preset pvdf-dft-valence-belt4 --variant belt --geometry pbe --systems beta   # 5.22: the Gamma-route tensor at belt4
python examples/crystal_stage_sensitivity.py --fit docs/predictions/crystal_stage_v8_tilts.json --out docs/predictions/crystal_stage_v8_tilts_sensitivity.json   # 5.22: the Jacobian at belt4, with the pendant-tilt correction carried
python examples/crystal_stage_sensitivity.py --analyze docs/predictions/crystal_stage_v8_tilts_sensitivity.json   # its linear analysis: the lateral-coupling lock stands at belt4
python examples/stage8_screen.py --fit docs/predictions/crystal_stage_v8_tilts.json --sensitivity docs/predictions/crystal_stage_v8_tilts_sensitivity.json --out docs/predictions/stage8_screen_belt4.json   # 5.22: the 5.20 screen at belt4 (no radial candidate adds a direction)
python examples/gamma_dielectric.py --preset pvdf-dft-valence-belt4 --variant belt --geometry pbe --systems beta,alpha --out docs/predictions/gamma_dielectric_pvdf-dft-valence-belt4.json   # 5.22: both tensors at belt4
python examples/gamma_dielectric.py --preset pvdf-dft-valence-belt4 --variant belt --geometry pbe --systems beta,alpha --k-bond-min 540 --out docs/predictions/gamma_dielectric_pvdf-dft-valence-belt4_kbond540.json   # 5.22: the same with the stretch constants floored (beta unchanged, alpha's polar eps_0 3.47 -> 2.93)
```

The result file holds the variables, the residual table, the Lennard-Jones terms and
increments it implies, and the anchors and tolerances it was fitted to. It is also the run's
checkpoint: the best evaluation so far is written there each time it improves, with
`"complete": false` and the evaluation count, and the final write marks it complete. A run cut
short (the first stage-5 run was, by a two-hour session limit at evaluation 292 of 320) resumes
from its checkpoint with `--x0`; the Powell line searches start over, so the resumed run is a
restart from the best point, not a continuation.
