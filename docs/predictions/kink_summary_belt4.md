**Lone kink pair in the all-trans crystal** (kcal/mol per kink pair; one chain of the supercell kinked, cell held at the all-trans values)

| sequence | chains x monomers | axial registry: barrier, at 1 A (per monomer, rigid slide) | neighbours fixed | neighbours rigid bodies | neighbours flexible |
|---|---|---|---|---|---|
| vdf10@0 | 4 x 10 | 6.80, 5.70 | +25.25 (gauche bonds 95, -91) | +23.07 (gauche bonds 99, -98) | reverts to all-trans (+0.00) |
| vdf10-trfe4-14@0 | 4 x 14 | 5.84, 4.91 | reverts to all-trans (-0.00) | - | - |
| vdf14@0 | 4 x 14 | - | +24.17 (gauche bonds 92, -88) | +22.69 (gauche bonds 95, -92) | reverts to all-trans (+0.03) |
| vdf20@0 | 4 x 20 | - | +23.41 (gauche bonds 90, -87) | +22.39 (gauche bonds 92, -89) | - |
| vdf9-trfe4-cfe1-14@0 | 4 x 14 | 4.56, 3.77 | +21.16 (gauche bonds 85, -85) | +20.13 (gauche bonds 86, -86) | reverts to all-trans (+0.00) |

**Kink sheet: every chain kinked in register, own cell** (kcal/mol per kink pair against the all-trans crystal)

| sequence | the search's cells (line-group all-trans, soft closure) | both cells polished, rigid pendants | with pendant tilts | kinked cell a, b, gamma, c | all-trans cell |
|---|---|---|---|---|---|
| vdf10@0 | +4.47 | +5.30 | +5.06 | 4.958, 8.774, 89.178, 25.201 | 4.754, 8.481, 90.000, 25.877 |
| vdf9-trfe4-cfe1-14@0 | +0.94 | +3.95 | +3.79 | 5.402, 8.426, 92.669, 35.496 | 5.302, 8.245, 91.700, 36.183 |

**One chain's kink pair one monomer ahead of its neighbours'** (kcal/mol per moving chain against the sheet; cell fixed)

| sequence | cell, moving chains | host | wide kink I | saddle A-I | saddle I-B | B (one monomer ahead) | barrier ahead / back | rate at 300 K, 3e12 /s attempt | bands converged |
|---|---|---|---|---|---|---|---|---|---|

**Margins with both cells polished** (kcal/mol per monomer, all-trans minus TG; positive when TG is lower)

| sequence | the screen's margin | polished | all-trans cell lower by | TG cell lower by | TG conformer | torsions off the ideal pattern, TT / TG | axial transition strain |
|---|---|---|---|---|---|---|---|
| vdf9-trfe4-cfe1-14 (VVRVVRVCVVRVVR) | +0.110 | +0.317 | 0.206 | 0.533 | TG'TG | 10 / 21 deg | +12.1% |
| vdf10-trfe4-14 (VVRVVRVVVVRVVR) | -0.051 | +0.046 | 0.081 | 0.382 | TG'TG | 6 / 22 deg | +11.4% |
| vdf9-trfe5-14 (VVRVVRVRVVRVVR) | -0.058 | +0.083 | 0.061 | 0.407 | TGTG' | 6 / 21 deg | +11.8% |
