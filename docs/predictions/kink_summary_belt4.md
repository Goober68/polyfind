**Lone kink pair in the all-trans crystal** (kcal/mol per kink pair; one chain of the supercell kinked, cell held at the all-trans values)

| sequence | chains x monomers | axial registry: barrier, at 1 A (per monomer, rigid slide) | neighbours fixed | neighbours rigid bodies | neighbours flexible |
|---|---|---|---|---|---|
| vdf10@0 | 4 x 10 | 6.80, 5.70 | +25.25 (gauche bonds 95, -91; 1 angles at the window) | +23.07 (gauche bonds 99, -98) | reverts to all-trans (+0.00) |
| vdf10@0, angle restraint 100 (chain stiffness matched to DFT) | 4 x 10 | 7.15, 5.93 | +20.32 (gauche bonds 88, -85; 5 angles at the window) | +18.55 (gauche bonds 89, -88; 3 angles at the window) | reverts to all-trans (+0.00) |
| vdf10@0, angle restraint 100 (chain stiffness matched to DFT), angle window 14 deg | 4 x 10 | - | +19.77 (gauche bonds 87, -82) | +18.23 (gauche bonds 89, -85) | - |
| vdf10-trfe4-14@0 | 4 x 14 | 5.84, 4.91 | reverts to all-trans (-0.00) | - | - |
| vdf14@0 | 4 x 14 | - | +24.17 (gauche bonds 92, -88) | +22.69 (gauche bonds 95, -92) | reverts to all-trans (+0.03) |
| vdf20@0 | 4 x 20 | - | +23.41 (gauche bonds 90, -87) | +22.39 (gauche bonds 92, -89) | - |
| vdf9-trfe4-cfe1-14@0 | 4 x 14 | 4.56, 3.77 | +21.16 (gauche bonds 85, -85) | +20.13 (gauche bonds 86, -86) | reverts to all-trans (+0.00) |

**Kink sheet: every chain kinked in register, own cell** (kcal/mol per kink pair against the all-trans crystal)

| sequence | the search's cells (line-group all-trans, soft closure) | both cells polished, rigid pendants | with pendant tilts | kinked cell a, b, gamma, c | all-trans cell |
|---|---|---|---|---|---|
| vdf10-trfe4-14-pairs@0 | +2.22 | +5.05 | +4.87 | 4.978, 8.688, 91.356, 35.516 | 4.848, 8.470, 91.791, 36.195 |
| vdf10-trfe4-14@0 | +5.68 | +7.09 | +6.85 | 5.150, 8.549, 92.162, 35.529 | 4.881, 8.459, 91.781, 36.187 |
| vdf10-trfe4-14@0~ | +3.70 | +5.79 | +5.70 | 4.984, 8.715, 90.352, 35.527 | 4.881, 8.459, 91.781, 36.187 |
| vdf10@0 | +4.47 | +5.30 | +5.06 | 4.958, 8.774, 89.178, 25.201 | 4.754, 8.481, 90.000, 25.877 |
| vdf9-trfe4-cfe1-14@0 | +0.94 | +3.95 | +3.79 | 5.402, 8.426, 92.669, 35.496 | 5.302, 8.245, 91.700, 36.183 |
| vdf9-trfe4-cfe1-14@0~ | +0.65 | +4.18 | +4.07 | 5.004, 9.075, 93.871, 35.529 | 5.302, 8.245, 91.700, 36.183 |

**One chain's kink pair one monomer ahead of its neighbours'** (kcal/mol per moving chain against the sheet; cell fixed)

| sequence | cell, moving chains | host | wide kink I | saddle A-I | saddle I-B | B (one monomer ahead) | barrier ahead / back | rate at 300 K, 3e12 /s attempt | bands converged |
|---|---|---|---|---|---|---|---|---|---|
| vdf10-trfe4-14@0 ahead | 1x1, 1 of 2 | flexible | +6.47 | +7.02 | +12.02 | +11.25 | 12.02 / 0.77 | 5.2e+03 /s | A-I force 0.21, I-B force 0.10 |
| vdf10@0 ahead | 1x1, 1 of 2 | flexible | +7.01 | +7.55 | +8.94 | +1.97 | 8.94 / 6.97 | 9.2e+05 /s | yes |
| vdf10@0 ahead | 1x1, all 2 (the sheet in register) | flexible | +5.06 | +8.23 | +7.87 | -0.00 | 8.23 / 8.23 | 3 /s | A-I force 1.04, I-B force 0.81 |
| vdf10@0 ahead | 2x1, 1 of 4 | flexible | +9.85 | +10.35 | +12.66 | +7.88 | 12.66 / 4.78 | 1.8e+03 /s | A-I force 0.10, I-B force 0.16 |
| vdf9-trfe4-cfe1-14@0 ahead | 1x1, 1 of 2 | flexible | +6.33 | +6.43 | +11.37 | +8.98 | 11.37 / 2.39 | 1.6e+04 /s | A-I force 0.10, I-B force 0.69 |

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
