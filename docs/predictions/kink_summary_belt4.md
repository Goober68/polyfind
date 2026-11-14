**Lone kink pair in the all-trans crystal** (kcal/mol per kink pair; one chain of the supercell kinked, cell held at the all-trans values)

| sequence | chains x monomers | axial registry: barrier, at 1 A (per monomer, rigid slide) | neighbours fixed | neighbours rigid bodies | neighbours flexible |
|---|---|---|---|---|---|
| vdf10@0 | 4 x 10 | 6.80, 5.70 | +25.25 (gauche bonds 95, -91; 1 angles at the window) | +23.07 (gauche bonds 99, -98) | reverts to all-trans (+0.00) |
| vdf10@0, angle restraint 100 (chain stiffness matched to DFT) | 4 x 10 | 7.15, 5.93 | +20.32 (gauche bonds 88, -85; 5 angles at the window) | +18.55 (gauche bonds 89, -88; 3 angles at the window) | reverts to all-trans (+0.00) |
| vdf10@0, angle restraint 100 (chain stiffness matched to DFT), angle window 14 deg | 4 x 10 | - | +19.77 (gauche bonds 87, -82) | +18.23 (gauche bonds 89, -85) | - |
| vdf10@0 | 16 x 10 | - | +25.61 (gauche bonds 94, -90; 1 angles at the window) | +22.01 (gauche bonds 106, -105; a damped descent stopped at its budget, NOT a minimum: gradient 0.09) | - |
| vdf10-trfe4-14@0 | 4 x 14 | 5.84, 4.91 | reverts to all-trans (-0.00) | - | - |
| vdf14@0 | 4 x 14 | - | +24.17 (gauche bonds 92, -88) | +22.69 (gauche bonds 95, -92) | reverts to all-trans (+0.03) |
| vdf20@0 | 4 x 20 | - | +23.41 (gauche bonds 90, -87) | +22.39 (gauche bonds 92, -89) | - |
| vdf9-trfe4-cfe1-14@0 | 4 x 14 | 4.56, 3.77 | +21.16 (gauche bonds 85, -85) | +20.13 (gauche bonds 86, -86) | reverts to all-trans (+0.00) |

Energy against the kink's own torsion (the two gauche bonds held at +tau and -tau degrees, the rest of the stage relaxed from the stage's stored state; kcal/mol above all-trans; examples/kink_torsion_scan.py). A rise before 120 degrees is the rim that holds the kink; a break in the curve is another branch.

* `vdf10@0 2x1 kinked 0`, rigid host (4 chains; stored state +25.25): 91: 25.37; 93: 25.33; 95: 25.34; 97: 25.40; 99: 25.49; 101: 25.62; 103: 25.77; 105: 25.94; 108: 26.21; 112: 26.58; 116: 26.90; 120: 19.12; 130: 14.09
* `vdf10@0 2x1 kinked 0`, host rigid bodies (4 chains; stored state +23.07): 97: 23.10; 99: 23.08; 100: 23.08; 101: 23.09; 102: 23.09; 103: 23.09; 104: 23.04; 105: 22.93; 106: 22.78; 108: 22.36; 112: 21.00; 120: 17.66; 130: 13.41
* `vdf10@0 4x2 kinked 0`, host rigid bodies (16 chains; stored state +22.01, a descent's last point and not a minimum): 100: 22.17; 112: 20.21; 120: 17.03; 130: 12.98

**One chain slid half a monomer along its axis, the lattice responding** (kcal/mol per monomer of the slid chain against the unslid lattice; cell fixed unless stated; examples/chain_slide.py)

| repeat and cell | chains | rigid lattice | neighbours turn | turn and move aside | and flex | with pendant tilts | turn and move aside, lateral cell free | the two nearest shells at 'turn and move aside': largest move, turn |
|---|---|---|---|---|---|---|---|---|
| vdf1 2x1 | 4 | 6.798 | 6.798 | 6.074 | - | - | 2.861 (a +5.19%, b +4.20%) | 4.75 A: 0.142 A, 7.5 deg; 4.86 A: 0.128 A, 11.5 deg |
| vdf1 4x2 | 16 | 6.797 | 6.027 | 3.696 | - | - | 3.093 (a +1.35%, b +1.08%) | 4.75 A: 0.166 A, 0.0 deg; 4.86 A: 0.137 A, 5.0 deg |
| vdf1 4x3 | 24 | 6.797 | 6.017 | 3.537 | - | - | - | 4.75 A: 0.168 A, 0.0 deg; 4.86 A: 0.163 A, 5.1 deg |
| vdf1 6x2 | 24 | 6.798 | 6.026 | 3.503 | - | - | - | 4.75 A: 0.192 A, 0.0 deg; 4.86 A: 0.140 A, 5.1 deg |
| vdf1 6x3 | 36 | 6.798 | 6.016 | 3.359 | - | - | 3.122 (a +0.58%, b +0.48%) | 4.75 A: 0.193 A, 0.0 deg; 4.86 A: 0.162 A, 5.2 deg |
| vdf1 6x4 | 48 | 6.798 | 6.016 | 3.312 | - | - | - | 4.75 A: 0.188 A, 0.0 deg; 4.86 A: 0.171 A, 5.1 deg |
| vdf1 8x4 | 64 | 6.798 | 6.016 | 3.263 | - | - | 3.135 (a +0.32%, b +0.27%) | 4.75 A: 0.196 A, 0.0 deg; 4.86 A: 0.171 A, 5.1 deg |
| vdf2 2x1 | 4 | 6.798 | 6.798 | 6.074 | 6.074 | 5.934 | - | 4.75 A: 0.142 A, 7.5 deg; 4.86 A: 0.128 A, 11.5 deg |
| vdf2 4x2 | 16 | 6.797 | 6.027 | 3.696 | 3.696 | 3.684 | - | 4.75 A: 0.166 A, 0.0 deg; 4.86 A: 0.137 A, 5.0 deg |
| vdf10 2x1 | 4 | 6.798 | - | 6.074 | 6.030 (not stationary along the slide) | 5.896 | - | 4.75 A: 0.188 A, 7.5 deg; 4.86 A: 0.128 A, 11.5 deg |
| vdf1 2x1 c2.5807 | 4 | 6.672 | - | 5.968 | - | - | 2.833 (a +5.12%, b +4.17%) | 4.76 A: 0.142 A, 7.4 deg; 4.87 A: 0.127 A, 11.4 deg |

Dilute limit, linear in one over the number of chains from the two largest cells of the 2:1 family: fixed cell 3.140 (from 36 and 64 chains), lateral cell free 3.152 (from 36 and 64): **3.15**, 46% of the rigid lattice's 6.80. The fixed cell bounds it from above and the free cell from below at every size.

Along the slide (slip in monomer periods: rigid lattice, relaxed):

* `vdf1 4x2`, shift: 0.0625: 0.19, 0.18; 0.125: 0.74, 0.66; 0.188: 1.66, 1.33; 0.25: 2.86, 2.04; 0.312: 4.23, 2.70; 0.375: 5.52, 3.24; 0.438: 6.45, 3.58; 0.5: 6.80, 3.70. Highest relaxed 3.70 at 0.5 (rigid 6.80 at 0.5); the two sweeps differ by 0.000 at most.
* `vdf10-trfe4-14 2x1`, shift cell free: 0.1: 0.43, 0.37; 0.2: 1.66, 1.10; 0.3: 3.44, 1.75; 0.4: 5.11, 2.18; 0.5: 5.84, 2.37; 0.6: 5.24, 2.33; 0.7: 3.73, 2.07; 0.8: 2.17, 1.60; 0.9: 1.18, 1.09; 1: 1.01, 0.95. Highest relaxed 2.37 at 0.5 (rigid 5.84 at 0.5); the two sweeps differ by 0.000 at most.
* `vdf9-trfe4-cfe1-14 2x1`, shift cell free: 0.1: 0.33, 0.29; 0.2: 1.28, 0.88; 0.3: 2.64, 1.43; 0.4: 3.93, 1.81; 0.5: 4.56, 2.01; 0.6: 4.30, 2.04; 0.7: 3.47, 1.90; 0.8: 2.63, 1.62; 0.9: 2.12, 1.30; 1: 2.00, 1.19. Highest relaxed 2.04 at 0.6 (rigid 4.56 at 0.5); the two sweeps differ by 0.000 at most.

**Kink sheet: every chain kinked in register, own cell** (kcal/mol per kink pair against the all-trans crystal)

| sequence | the search's cells (line-group all-trans, soft closure) | both cells polished, rigid pendants | with pendant tilts | kinked cell a, b, gamma, c | all-trans cell |
|---|---|---|---|---|---|
| vdf10-trfe4-14-pairs@0 | +2.22 | +5.05 | +4.87 | 4.978, 8.688, 91.356, 35.516 | 4.848, 8.470, 91.791, 36.195 |
| vdf10-trfe4-14@0 | +5.68 | +7.09 | +6.85 | 5.150, 8.549, 92.162, 35.529 | 4.881, 8.459, 91.781, 36.187 |
| vdf10-trfe4-14@0~ | +3.70 | +5.79 | +5.70 | 4.984, 8.715, 90.352, 35.527 | 4.881, 8.459, 91.781, 36.187 |
| vdf10@0 | +4.47 | +5.30 | +5.06 | 4.958, 8.774, 89.178, 25.201 | 4.754, 8.481, 90.000, 25.877 |
| vdf9-trfe4-cfe1-14@0 | +0.94 | +3.95 | +3.79 | 5.402, 8.426, 92.669, 35.496 | 5.302, 8.245, 91.700, 36.183 |
| vdf9-trfe4-cfe1-14@0~ | +0.65 | +4.18 | +4.07 | 5.004, 9.075, 93.871, 35.529 | 5.302, 8.245, 91.700, 36.183 |

**One chain's kink pair one monomer ahead of its neighbours'** (kcal/mol per moving chain against the sheet; cell fixed). A saddle is given as the band's top and, in bold where it was done, the stationary point settled from it (minimum in the hyperplane normal to the band's tangent, zero force along it); the barrier and the rate are on the settled saddles where both exist.

| sequence | cell, moving chains | host | wide kink I | saddle A-I | saddle I-B | B (one monomer ahead) | barrier ahead / back | rate at 300 K, 3e12 /s attempt | bands' forces |
|---|---|---|---|---|---|---|---|---|---|
| vdf10-trfe4-14@0 ahead | 1x1, 1 of 2 | flexible | +6.47 | +7.02, **+7.05** | +12.02, **+12.02** | +11.25 | 12.02 / 0.77 | 5.2e+03 /s | 0.21, 0.10 |
| vdf10@0 ahead | 1x1, 1 of 2 | flexible | +7.01 | +7.55, **+7.54** | +8.94, **+8.91** | +1.97 | 8.91 / 6.95 | 9.6e+05 /s | 0.10, 0.10 |
| vdf10@0 ahead | 1x1, all 2 (the sheet in register) | flexible | +5.06 | +8.23, **+7.38** | +7.87, **+7.34** (polish not converged) | -0.00 | 8.23 / 8.23 | 3 /s | 1.04, 0.81 |
| vdf10@0 ahead | 2x1, 1 of 4 | flexible | +9.85 | +10.29, **+10.27** | +12.49, **+12.46** | +7.88 | 12.46 / 4.59 | 2.5e+03 /s | 0.03, 0.08 |
| vdf10@0 ahead | 4x2, 1 of 16 | flexible | +9.96 | +10.73, **+10.66** | +11.84, **+11.64** | +5.91 | 11.64 / 5.73 | 9.9e+03 /s | 0.09, 0.09 |
| vdf9-trfe4-cfe1-14@0 ahead | 1x1, 1 of 2 | flexible | +6.33 | +6.43, **+6.43** | +11.37, **+11.48** | +8.98 | 11.48 / 2.50 | 1.3e+04 /s | 0.10, 0.69 |

End points only (no band)

| sequence | cell, moving chains | host | wide kink I | B (one monomer ahead) |
|---|---|---|---|---|
| vdf10@0 ahead | 2x1, 1 of 4 | bodies | +17.81 | +11.51 |
| vdf10@0 ahead | 4x2, 1 of 16 | bodies | +17.36 | +9.32 |
| vdf10@0 ahead | 4x1, 1 of 8 | flexible | +10.02 | +6.47 |
| vdf10@0 ahead | 2x2, 1 of 8 | flexible | +9.98 | +7.51 |
| vdf10@0 ahead | 8x1, 1 of 16 | flexible | +10.05 | +6.19 |

**The single-chain step against the size of the cell** (PVDF, chain 0's kink pair one monomer ahead of every other chain's; kcal/mol against the sheet; cell fixed)

| cell (a x b tiling of the two-chain cell) | one chain in | neighbours rigid chains: I, B | neighbours flexible: I | saddle A-I | saddle I-B | barrier ahead / back | B |
|---|---|---|---|---|---|---|---|
| 1x1 | 2 (its a-neighbours move too: a stack) | - | +7.01 | +7.55, **+7.54** | +8.94, **+8.91** | 8.91 / 6.95 | +1.966 |
| 2x1 | 4 | +17.81, +11.51 | +9.85 | +10.29, **+10.27** | +12.49, **+12.46** | 12.46 / 4.59 | +7.876 |
| 4x1 | 8 | - | +10.02 | - | - | - | +6.470 |
| 2x2 | 8 | - | +9.98 | - | - | - | +7.512 |
| 8x1 | 16 | - | +10.05 | - | - | - | +6.191 |
| 4x2 | 16 | +17.36, +9.32 | +9.96 | +10.73, **+10.66** | +11.84, **+11.64** | 11.64 / 5.73 | +5.911 |
| 8x2 | 32 | - | +9.97 | **+10.66** (from the 4x2 cell's) | **+11.47** (from the 4x2 cell's) | 11.47 / 5.78 | +5.692 |
| 4x4 | 32 | - | - | - | - | - | +5.712 |
| 8x4 | 64 | - | +9.89 | **+10.63** (from the 8x2 cell's) | **+11.53** (from the 8x2 cell's) | 11.53 / 6.04 | +5.495 |
| 8x8 | 128 | - | - | - | - | - | +5.430 |

How the lattice answers in the 4x2 cell (16 chains, neighbours flexible): displacements against the kinked crystal, relative to the median of the chains whose kink pair stays

| state | chain or shell (distance from the moving chain) | chains | moves aside by, at most | turns, at most | shifts along the axis, at most | largest torsion change |
|---|---|---|---|---|---|---|
| wide kink I | the moving chain (gauche bonds at +95, -84 deg) | 1 | 0.17 A | -5.3 deg | -0.06 A | - |
| | 4.96 A | 6 | 0.07 A | 1.8 deg | -0.07 A | 16.2 deg |
| | 8.58 A | 5 | 0.03 A | 0.5 deg | +0.02 A | 3.1 deg |
| | 9.92 A | 2 | 0.01 A | 0.1 deg | +0.02 A | 0.9 deg |
| | 10.14 A | 1 | 0.00 A | 0.3 deg | -0.01 A | 4.8 deg |
| | 13.33 A | 1 | 0.01 A | 0.0 deg | +0.01 A | 0.6 deg |
| B, one monomer ahead | the moving chain (gauche bonds at +86, -64 deg) | 1 | 0.11 A | +3.9 deg | -0.32 A | - |
| | 4.96 A | 6 | 0.08 A | 0.9 deg | +0.19 A | 10.5 deg |
| | 8.58 A | 5 | 0.03 A | 0.4 deg | +0.06 A | 2.9 deg |
| | 9.92 A | 2 | 0.02 A | 0.2 deg | +0.04 A | 1.8 deg |
| | 10.14 A | 1 | 0.02 A | 0.4 deg | -0.01 A | 3.0 deg |
| | 13.33 A | 1 | 0.01 A | 0.1 deg | +0.02 A | 0.3 deg |

**Two chains one monomer ahead together** (end state only; kcal/mol against the sheet, neighbours flexible, cell fixed)

| cell | the two chains, and how far apart | B for the two | two separate single steps | gained by stepping side by side | per chain |
|---|---|---|---|---|---|
| 4x2 (16 chains) | 0,1, 5.07 A | +9.05 | +11.82 | +2.77 | +4.53 |
| 4x2 (16 chains) | 0,4, 4.96 A | +8.35 | +11.82 | +3.47 | +4.17 |

**A second chain stepping beside one that has stepped** (kcal/mol against the START, the sheet with the first chain one monomer on; every chain flexible, cell fixed; saddles as in the first table. The last columns set the neighbour's barrier to join against the first chain's own barrier to come back, at 300 K with the same attempt frequency for both)

| cell | the chain that steps; the one already ahead, and where it is | the start above the sheet | wide kink I | saddle A-I | saddle I-B | B (both ahead) | barrier to join / to part again | both ahead, above the sheet | the first chain's barrier to come back | it comes back, for each time this neighbour joins |
|---|---|---|---|---|---|---|---|---|---|---|
| 4x2 (16 chains) | 12; 0, 4.96 A away, at (+4.96, +0.00) A from it | +5.91 | none: no minimum there, the chain slides on to B | +11.67, **+10.96** | none (downhill) | +2.44 | 10.96 / 8.52 | +8.35 | 5.73 | 6.5e+03 times |
| 4x2 (16 chains) | 4; 0, 4.96 A away, at (-4.96, +0.00) A from it | +5.91 | +9.40 | +10.52, **+10.43** | +13.35, **+13.12** | +2.44 | 13.12 / 10.68 | +8.35 | 5.73 | 2.4e+05 times |

**The kink sheet along the sequence** (kcal/mol per kink pair against the all-trans crystal, rigid pendants; the kink pair's first unit is the monomer named, gauche bonds on the even bonds; `~` is the mirrored kink pair G- T G+)

* `VVVVVRRVVVVVRR` G+ T G- (14 of 14 positions): VV 5.05 VV 5.11 VV 5.40 VV 7.83 VR 10.41 RR 7.37 RV 7.16 VV 5.05 VV 5.11 VV 5.40 VV 7.83 VR 10.41 RR 7.37 RV 7.16
* `VVRVVRVVVVRVVR` G+ T G- (14 of 14 positions): VV 7.09 VR 6.50 RV 7.20 VV 7.05 VR 6.47 RV 6.87 VV 4.81 VV 5.14 VV 7.09 VR 6.47 RV 7.20 VV 7.09 VR 6.50 RV 7.21
* `VVRVVRVVVVRVVR` G- T G+ (14 of 14 positions): VV 5.79 VR 8.27 RV 5.42 VV 5.81 VR 8.28 RV 5.08 VV 4.66 VV 4.97 VV 5.74 VR 8.40 RV 5.41 VV 5.79 VR 8.28 RV 5.41
* `VVRVVRVCVVRVVR` G+ T G- (14 of 14 positions): VV 3.95 VR 5.30 RV 4.56 VV 4.06 VR 5.25 RV 5.88 VC 12.11 CV 11.41 VV 4.75 VR 5.45 RV 4.47 VV 3.96 VR 5.33 RV 4.51
* `VVRVVRVCVVRVVR` G- T G+ (14 of 14 positions): VV 4.18 VR 5.45 RV 3.73 VV 4.05 VR 3.57 RV 3.42 VC 8.07 CV 8.28 VV 4.14 VR 4.46 RV 3.83 VV 4.23 VR 4.36 RV 3.81

| sequence | scans (handednesses) | lowest site | highest site | sheets per repeat at 300 K, sum over sites of exp(-E/kT) | against 64/29/7 |
|---|---|---|---|---|---|
| vdf10-trfe4-14-pairs (VVVVVRRVVVVVRR) | 1 | 5.05 | 10.41 | 1.05e-03 per handedness | 0.095 |
| vdf10-trfe4-14 (VVRVVRVVVVRVVR) | 2 | 4.66 | 8.40 | 1.02e-03 per handedness | 0.092 |
| vdf9-trfe4-cfe1-14 (VVRVVRVCVVRVVR) | 2 | 3.42 | 12.11 | 1.12e-02 per handedness | 1.000 |

**Margins with both cells polished** (kcal/mol per monomer, all-trans minus TG; positive when TG is lower)

| sequence | the screen's margin | polished | all-trans cell lower by | TG cell lower by | TG conformer | torsions off the ideal pattern, TT / TG | axial transition strain |
|---|---|---|---|---|---|---|---|
| vdf9-trfe4-cfe1-14 (VVRVVRVCVVRVVR) | +0.110 | +0.317 | 0.206 | 0.533 | TG'TG | 10 / 21 deg | +12.1% |
| vdf10-trfe4-14 (VVRVVRVVVVRVVR) | -0.051 | +0.056 | 0.081 | 0.375 | TGTG' | 6 / 19 deg | +11.7% |
| vdf9-trfe5-14 (VVRVVRVRVVRVVR) | -0.058 | +0.083 | 0.061 | 0.407 | TGTG' | 6 / 21 deg | +11.8% |
| vdf9-trfe4-cdfe1-14 (VVRVVRVDVVRVVR) | +0.402 | +0.685 | 0.150 | 0.637 | TGTG' | 14 / 25 deg | +12.1% |
| vdf7-trfe3 (VVRVVRVVVR) | +0.145 | +0.213 | 0.078 | 0.305 | TGTG' | 6 / 19 deg | +11.6% |
| vdf8-trfe2 (VVRVVVVRVV) | +0.177 | +0.229 | 0.097 | 0.241 | G'TGT | 6 / 14 deg | +11.3% |
| vdf7-trfe2-cfe1 (VVRVVCVVRV) | +0.966 | +0.868 | 0.248 | 0.154 | TG'TG | 11 / 13 deg | +11.6% |
| vdf7-trfe2-cdfe1 (VVRVVDVVRV) | +0.146 | +0.392 | 0.220 | 0.719 | TGTG' | 16 / 23 deg | +12.8% |
| vdf13-trfe6-cfe1-20 (VVRVVRVVRVCVVRVVRVVR) | +0.055 | +0.258 | 0.175 | 0.517 | TGTG' | 10 / 21 deg | +12.1% |
| vdf6-trfe4 (VRVRVVRVRV) | -0.329 | -0.199 | 0.049 | 0.414 | TGTG' | 6 / 20 deg | +11.0% |
| vdf8-trfe4-cfe1-cdfe1-14 (DVRVVRVCVVRVVR) | +0.641 | +0.807 | 0.251 | 0.581 | TG'TG | 19 / 20 deg | +12.2% |
| vdf9-trfe4-ctfe1-14 (VVRVVRVKVVRVVR) | +0.135 | +0.300 | 0.186 | 0.664 | TG'TG | 14 / 24 deg | +11.8% |
| vdf9-trfe4-an1-14 (VVRVVRVAVVRVVR) | -0.154 | +0.197 | 0.087 | 0.614 | TG'TG | 5 / 24 deg | +12.5% |
