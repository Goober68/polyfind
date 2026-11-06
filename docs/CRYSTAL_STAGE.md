# The crystal stage of the fit: nonbonded terms against the DFT crystal anchors

*2026-11-05. Code: `polyfind.crystalfit`, driver `examples/crystal_stage.py`, result
`docs/predictions/crystal_stage_v1.json`. Preset: `pvdf-dft-valence-crystal` (opt-in; the
E1/E2 predictions pre-registered on `pvdf-dft-valence` stand unchanged).*

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
radius that bond inference never reaches, a polarizability of 1e-6 A^3 (a strict zero would
put an infinity on the dipole system's diagonal), zero mass. A chain with sites at zero depth
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

### 5.11 Stage 6: the polar Berry slope as an anchor (implemented; runs after stage 5)

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

Stage 6 (both Berry anchors, from the stage-5 point) is running; its result goes here.

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

### 5.13 Stage 7: the restraint as a variable, the backbone angles as anchors (implemented; runs after stage 6)

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
angles at DFT's values is the question stage 7 answers; it runs from the stage-6 result and its
result goes here.

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
  docs/predictions/e1_polyfind_crystal.json`.

## 7. Reproduce

```
python examples/crystal_stage.py --maxfev 300 --out docs/predictions/crystal_stage_v1.json
python examples/crystal_stage.py --stage 2 --maxfev 320 --out docs/predictions/crystal_stage_v2.json
python examples/crystal_stage.py --stage 3 --maxfev 320 --out docs/predictions/crystal_stage_v3.json
python examples/crystal_stage.py --stage 4 --maxfev 320 --out docs/predictions/crystal_stage_v4.json
python examples/crystal_stage.py --stage 4 --x0 docs/predictions/crystal_stage_v4.json --maxfev 320 --out docs/predictions/crystal_stage_v5.json   # the Berry P anchor
python examples/crystal_stage.py --stage 6 --x0 docs/predictions/crystal_stage_v5.json --maxfev 300 --out docs/predictions/crystal_stage_v6.json   # plus the polar Berry slope
python examples/crystal_stage.py --stage 7 --geometry pbe --x0 docs/predictions/crystal_stage_v6.json --maxfev 300 --out docs/predictions/crystal_stage_v7.json   # plus the restraint's stiffness, Dipole's backbone angles as anchors
python examples/e1_predictions.py --preset pvdf-dft-valence-belt --variant belt --out docs/predictions/e1_polyfind_belt.json
python examples/e2_predictions.py --preset pvdf-dft-valence-belt --variant belt --out docs/predictions/e2_polyfind_belt.json
python examples/crosscheck_e1_e2.py --polyfind docs/predictions/e1_polyfind_belt.json --e2 docs/predictions/e2_polyfind_belt.json
```

The result file holds the variables, the residual table, the Lennard-Jones terms and
increments it implies, and the anchors and tolerances it was fitted to. It is also the run's
checkpoint: the best evaluation so far is written there each time it improves, with
`"complete": false` and the evaluation count, and the final write marks it complete. A run cut
short (the first stage-5 run was, by a two-hour session limit at evaluation 292 of 320) resumes
from its checkpoint with `--x0`; the Powell line searches start over, so the resumed run is a
restart from the best point, not a continuation.
