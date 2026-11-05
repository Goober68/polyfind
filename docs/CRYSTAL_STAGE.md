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
construction): PVDF's ground conformation is still TGTG' and the margins open, TTTT +1.30 ->
+1.63 and T3GT3G' +0.73 -> +1.04 kcal/mol per monomer against DFT's +0.726 for the alpha-beta
pair -- the fixed-charge fit reproduces that margin to 0.02 and the induced dipoles, on the
larger charges, add 0.9 kcal/mol of polar-lattice stabilisation that the fit never saw. The
copolymer's margins open the same way, TTTT +1.09 -> +1.33 and T3GT3G' +0.56 -> +0.86 kcal/mol
per monomer over TGTG' (conformational strain of the switch +9.9%, unchanged), so on this
preset the all-trans phase is a quarter further from the ground state than the pre-registered
set says, in both polymers.


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

RESULTS2-TODO

## 6. What it costs, and what it does not touch

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
```

The result file holds the variables, the residual table, the Lennard-Jones terms and
increments it implies, and the anchors and tolerances it was fitted to.
