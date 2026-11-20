# A machine-learned potential as a second opinion on the margins

State as of 2026-11-20. The PI's direction, in order: "Go on 1" (a machine-learned interatomic
potential as a second opinion on polyfind's margins, tested before it is trusted), then, on
seeing pretrained ones tested: "You're using ots solutions? Thought we would train our own?"

Section 1 is the test. Section 2 is done: no pretrained potential passes it. Section 3 is our
own, and says for each part whether it is done, running or not started. Where it stands on a
real candidate is 3.8: for PVDF with reversed units, DFT meets polyfind's margin when the cell
is shaped as polyfind shapes it (+0.070 against +0.064) and puts it 0.33 lower when the cell
is free to shear along the chains, a freedom polyfind's cell lacks. In short: two short
fits on 76 PVDF frames were right to 0.1 to 0.2 kcal/mol per monomer where the pretrained ones
are off by 1.2 to 3.9 (3.3). Four long fits on a rented card, on 130 frames, are right to 0.01
on the packings that have frames near them and off by 0.07 to 0.21 on a packing they were not
shown; DFT single points on the cells they relax are right to 0.02 either way (3.6). So the
potential is used to relax a cell and DFT is asked once for its energy; the potential's own
energy is trusted only where it has been checked. Since 2026-11-19 (3.10) one fit on every
frame is in the repository, and two copolymer margins have been taken that way: a 25 mol% TrFE
cell, +0.678, which both of polyfind's presets meet to 0.04, and the same cell with a CFE
unit, +0.682, which one preset puts 0.30 lower and the other 0.38 higher. The error is at the
bonds of the CFCl carbon; with a torsion term of their own, fitted to DFT on lone chains, the
new preset gives +0.690 on that cell and passes on a CTFE cell too (3.11). A nitrile cell is the
other way round: DFT has all-trans below every helix (-0.386) where both presets have a helix
lowest, and the error is between chains as well as inside (3.12). The 25 rows whose units are
checked (27 once ethylene and DFE were) have their margins under the new preset in 3.13, with where each unit stands. Four
fourteen-monomer rows relaxed in the potential (3.14) agree with that preset on what a CFE does
and put every level a quarter of a kcal/mol lower. DFT on a six-monomer cell of the matrix's
motif says the potential has the level (-0.044 in DFT, -0.096 in the potential, +0.139 in the
preset): the preset's margins rest on two errors that cancel in the four-monomer cells and do
not in the six-monomer one (3.15). A TrFE term fitted to lone chains moves the error to the
four-monomer cells and is not adopted; no evaluator is within 0.1 on every cell, and DFT at
fourteen monomers needs more memory than this machine has (3.16). Kinks have been put to DFT
since (3.17 to 3.19), and the finding is about the structure before it is about any evaluator:
the kink sheets of polyfind's tables are upright, the kinks of all chains at one height, and a
sheet free to shape its cell inclines and ends 6.8 kcal/mol per kink pair lower in PVDF, below
the all-trans crystal (-1.96 in DFT, -1.66 in the potential, -1.17 and -0.78 in polyfind's two
presets). The frequencies and losses in `docs/CANDIDATE_MATERIALS.md` are made of the upright
ones.

## 1. The test

Polyfind ranks candidates on the margin, E(all-trans) - E(helix) per monomer. The switching
window at 55 V/um is 0.12 kcal/mol wide (a tenth of the cell polar at +0.15, half at +0.07), so
a second opinion has to resolve about **0.1 kcal/mol per monomer** to say anything about where
a candidate sits.

The project holds one such margin from DFT: PVDF, beta minus alpha, **+0.726 kcal/mol per
monomer** in PBE-D3(BJ) (`docs/POLAR_MARGIN_DFT.md`, "The alpha/beta control"; -31.46 meV as
alpha minus beta, of which D3 is +20.09). Polyfind's own value is +0.744. A potential that
misses this number by more than 0.1 cannot referee the candidates' margins, which are the same
quantity in cells whose pendants it has seen less of than PVDF's.

Published functionals put the same gap between 0.62 and 1.55 kcal/mol
(`docs/REFERENCES.md` section 5: 2.6 to 6.5 kJ/mol). That is the most a potential fitted to a
different functional could differ for an honest reason.

## 2. Pretrained potentials: none passes (done)

`examples/mlip_second_opinion.py`; numbers in `docs/predictions/mlip_second_opinion.json`;
`tests/test_mlip_second_opinion.py`. Two readings per potential: single points on the three
DFT-relaxed cells, and each cell relaxed in the potential itself (atoms and lattice, from the
DFT geometry), so that a potential is not charged for a geometry it would not choose.

Only potentials whose licence permits this project's use were run: the three MACE-MP
foundation models (MIT; fitted to PBE without dispersion, so D3(BJ) is added as the anchor has
it) and Rowan's Egret-1 family (MIT; fitted to wB97M-D3BJ on organic molecules, run as they
are).

| potential | margin on the DFT geometries | relaxed in the potential | error, in units of the 0.1 needed | antipolar - polar beta, meV | volume against DFT's, relaxed (alpha / beta) |
|---|---:|---:|---:|---:|---:|
| **DFT, PBE-D3(BJ)** | **+0.73** | **+0.73** | | **-1.1** | |
| polyfind | +0.74 | | 0.2 | -0.2 to -0.7 | |
| MACE-MPA-0 + D3 | +2.56 | +1.96 | 12 to 18 | -16 / +20 | +12% / +20% |
| MACE-MP-0b3 + D3 | -1.99 | -1.42 | 21 to 27 | +30 / -1 | +0% / -0% |
| MACE-MP-0a + D3 | +2.19 | +2.76 | 15 to 20 | -11 / -16 | -13% / -8% |
| Egret-1 | +3.96 | +3.72 | 30 to 32 | -185 / -192 | -7% / +4% |
| Egret-1T | +4.64 | +2.87 | 21 to 39 | -219 / -126 | +3% / +15% |

(kcal/mol per monomer; the antipolar column is single points / relaxed. Egret-1E returns
energies thousands of kcal/mol apart for these cells and is not usable on them at all.)

What the table says:

- Every one is off by 1.2 to 3.9 kcal/mol, twelve to thirty-nine times the resolution needed.
  They do not agree with each other on the sign: MACE-MP-0b3 puts all-trans below the helix.
- Relaxing in the potential does not rescue any of them. Three of the five move the crystal's
  volume by 8 to 20%, so their packing of the chains is not the DFT's.
- The errors are outside the 0.62 to 1.55 spread of published functionals, so they are the
  potentials', not a difference of functional.
- The D3 add-on alone gives +20.05 meV for alpha minus beta against the DFT margin's own D3
  part of +20.09. The long-range part is right; the error is in the network's part.
- The antipolar/polar pair, which DFT puts 1.1 meV apart, comes out 11 to 219 meV apart. That
  pair is one chain's rotation in the same lattice: what a pin does to a margin is a
  difference of that kind.

Not run, and why: the organic-chemistry models under an academic licence (MACE-OFF23) and the
other academic or gated foundation models (OMAT, MATPES, OMOL, UMA). Using them is the PI's
decision. MACE-OFF23 is fitted to the same family of data as Egret-1, so the Egret rows are
the best available guess at how it would do.

## 3. Our own

A potential fitted to this project's own DFT, in the anchor's recipe, on this project's own
chemistry. A pretrained network may still be the starting point (fitting on from one needs
fewer frames than starting from nothing), and both are to be tried on the same frames.

### 3.1 The DFT recipe here (done)

The anchor's recipe runs in this session: Quantum ESPRESSO 7.5 from conda-forge, the SSSP
1.3.0 PBE precision files (checksums equal to the anchor's record), PBE with two-body D3(BJ),
90/360 Ry. The anchor's own beta-PVDF input returns the anchor's energy to the last digit
(-275.80701219 Ry, D3 -0.05243146), in 52 s on four cores. `examples/mlip_frames.py check`
does this and no label is computed unless it passes.

### 3.2 Frames (PVDF: done, 76 frames, 2.5 hours of DFT; the other units: 3.4)

A frame is a cell with its DFT energy, the force on every atom and the stress. The project
held none: the stored relaxations keep their end points only, here and in the exchange.

`examples/mlip_frames.py make / label / collect`; the labelled frames are
`docs/predictions/mlip/pvdf_frames_pbe_d3bj.extxyz`. The PVDF set is 76 frames of 12 and 24
atoms from the four DFT-relaxed cells (helix alpha, polar beta, two antipolar betas):

| family | frames | what it carries |
|---|---:|---|
| `end` | 4 | the relaxed cells themselves |
| `strain` | 24 | 2% along each axis, both signs, and a lateral shear: the lattice's stiffness |
| `rattle` | 18 | a random strain up to 2%, every atom moved by a Gaussian of 0.05 A: bonds and angles |
| `md` | 30 | snapshots of dynamics at 500 K with a random strain: torsions, the chains' rotation and sliding |

The dynamics run in a pretrained potential (MACE-MP-0b3), used as a sampler only: it decides
which geometries are labelled, not what their labels are.

What this set does not hold, and so what a potential fitted to it cannot be asked: any unit
but VDF, and any geometry on the way between all-trans and the helix. A margin of a candidate
needs the first; a barrier needs the second.

### 3.3 Training and its tests (the two short fits; the long ones are in 3.6)

`examples/mlip_train.py prepare / train / evaluate / closure`; results in
`docs/predictions/mlip_own_pvdf_ft_0b3.json` and `..._0b3_b.json`; `tests/test_mlip_own.py`.

The network is fitted to the PBE part of each label; the D3(BJ) part is taken out before the
fit and added back analytically when the potential is used, as the pretrained MACE-MP models
are built. The fit starts from a pretrained network (MACE-MP-0b3, MIT licence) with every
weight free. Held out of the fit: the four relaxed cells, and one frame in nine of the rest
(64 frames in the fit, 12 out).

The tests, in order of what they decide:

1. beta minus alpha against +0.726 kcal/mol per monomer, to 0.1;
2. antipolar minus polar beta against -1.1 meV;
3. beta-PVDF's elastic constants against DFT's C_11 22.8, C_22 24.5, C_12 5.9, C_33 303 GPa
   (polyfind's rigid bonds cannot give C_33 at all);
4. the force error on frames left out of the fit;
5. closure: each cell relaxed in the potential is given back to DFT for a single point. The
   DFT force left there and the DFT energy above DFT's own minimum say whether the potential's
   minimum is the DFT's;
6. a packing no frame was made from: the helix cell with one chain turned half a turn, so the
   two chains' dipoles add (delta-PVDF's kind). Relaxed in the potential, then a DFT single
   point: its gap to beta by the potential against the same gap by DFT.

What the tests are worth. Tests 1 to 3 are on cells whose strained and shaken copies are in
the fit, so they show the potential can hold the resolution between frames it was shown
neighbours of: an interpolation. Tests 5 and 6 are on geometries the potential chose itself,
and 6 on a packing it never saw. None of them says a candidate's margin is right: that needs
frames with the candidate's units in them, and then a DFT single point on a candidate's own
cells as the test.

**The two fits.** The first: 60 epochs, one hour on this container's four cores, the last 24
with the energy weighted up. The second: 30 more epochs on from the first (37 minutes) with
energy weighted 300 times and stress 1000 times more, because in the first fit the two
together were under 1% of what was being minimised.

| test | DFT | first fit | second fit |
|---|---:|---:|---:|
| 1. beta minus alpha on the DFT geometries, kcal/mol | +0.726 | +0.671 | +0.617 |
| 1. the same, relaxed in the potential | +0.726 | +0.700 | +0.566 |
| 2. antipolar minus polar beta on the DFT geometries, meV | -1.1 | +2.7 | +4.2 |
| 3. C_11 / C_22 / C_12 / C_33, GPa | 22.8 / 24.5 / 5.9 / 303 | 39.4 / 25.7 / 9.8 / 280 | 31.1 / 19.0 / 9.1 / 271 |
| 3. volume of its relaxed cells against DFT's (alpha / beta / antipolar) | | -4.6% / -2.1% / -0.1% | -0.4% / +2.0% / +6.5% |
| 4. held out: energy, meV per atom; force, meV/A (of 1400 rms); stress, GPa | | 0.39; 37; 0.40 | 0.55; 39; 0.27 |
| 5. DFT force left on its relaxed cells, eV/A (on the turned-chain cell) | 0 | 0.01 to 0.02 (0.08) | 0.02 to 0.05 (0.06) |
| 5. DFT energy of its relaxed cells above DFT's minima, meV per monomer | 0 | 2.5 to 5.5 | 0.4 to 13.8 |
| 5. alpha minus beta on its own cells: the potential less DFT, kcal/mol | | -0.017 | +0.229 |
| 5. antipolar minus polar beta on its own cells: the potential less DFT | | +0.170 | -0.246 |
| 6. turned-chain helix minus beta on its own cells: the potential less DFT | | -0.011 | +0.105 |

What the table says:

- **An order of magnitude better than the pretrained potentials, and not yet inside 0.1.** The
  gaps between cells come out within 0.01 to 0.25 kcal/mol per monomer of DFT, where the
  pretrained potentials were off by 1.2 to 3.9. The first fit is inside 0.1 on the helix
  against all-trans, both for alpha and for the packing it never saw; the second is not; and
  neither holds the gap between two packings of the same all-trans chains.
- **The two fits differ from each other by as much as the resolution needed.** That is the
  honest size of the error at this amount of training: about 0.1 to 0.2, whichever fit is
  taken. The first fit's two readings of 0.01 to 0.02 are not to be quoted as its accuracy.
- **Its atoms are right and its lattice is not.** DFT finds almost no force on the atoms of
  the cells the potential relaxed (0.01 to 0.08 eV/A), but the cells are 0 to 7% off in volume
  and C_11 is 37 to 73% too stiff. The lattice across the chains is held by energies of a
  quarter of a meV per atom at the 2% strains in the frames, below either fit's energy error.
- **Weighting the energy and the stress up did not fix it.** The stress error fell by a third
  and the energies got worse. The fits are short: 60 epochs on 58 frames in batches of four is
  under a thousand steps of the optimiser, where fits of this kind are usually run for many
  thousands. On four cores a step takes four seconds.

A by-product. DFT puts the turned-chain helix cell (the two chains' dipoles adding, as in
delta-PVDF) level with alpha: 0.4 meV per monomer below the alpha cell's own DFT minimum on
the second fit's geometry, before any DFT relaxation of it, and 2.0 above on the first fit's.
The anchor's +0.726 is all-trans against alpha; against the lowest helix packing found it is
+0.735 at least. The two DFT single points are in
`docs/predictions/mlip/pvdf_probe_frames_pbe_d3bj.extxyz` with the six others.

### 3.4 The other units (done: 168 frames labelled, 13.5 hours of DFT here, no rented machine)

The labelled frames are `docs/predictions/mlip/units_frames_pbe_d3bj.extxyz`.

`examples/mlip_frames.py make --set units`. TrFE, CFE, CTFE, acrylonitrile, the reversed VDF
unit and TFE are put into the DFT-relaxed PVDF cells by replacing a unit's groups along their
own bonds; the cell is then relaxed and shaken in the pretrained sampler (atoms and lattice,
with D3), which again only chooses the geometries.

| cells | atoms | what they hold | frames |
|---|---:|---|---:|
| two monomers per chain, all-trans and helix, 12 compositions | 24 to 27 | each unit alone in PVDF (TrFE also at 50 mol%); pairs of units facing each other on neighbouring chains | 120 |
| four monomers per chain, all-trans and helix, 6 compositions | 48 to 51 | the leading cells' units together: a TrFE matrix; the nitrile with a reversed unit; CFE in TrFE; CTFE with a reversed unit; PVDF with reversed units; the nitrile with TFE | 48 |

Each cell gives its relaxed geometry, one shaken copy and two or three dynamics snapshots.

The sampler's cells are not to be trusted for their volume. It kept every chain's conformation
(all-trans stayed all-trans, the helix stayed a helix), but it put most all-trans cells holding
CTFE, the nitrile, the reversed unit or TFE 22 to 37% above PVDF's volume (64 to 72 A3 per
monomer, an amorphous density) while leaving the CFE and TrFE ones within 4%. So the order of
labelling is: each cell's relaxed frame first; then `mlip_frames.py aim` reads the DFT stress on
it and, with an assumed stiffness across the chains, moves the cell's other frames to half, all
and one and a half times the strain the stress points to; then those are labelled. The DFT
minimum is then inside what is labelled whichever way the sampler was wrong.

It was wrong both ways. On the 36 relaxed cells (labelled; 3.0 hours of DFT) the DFT stress
across the chains is -1.3 to -2.6 GPa on the all-trans cells the sampler left at PVDF's density
(TrFE and CFE: compressed; by the assumed stiffness they want 7 to 12% more along each lateral
axis) and -0.4 to +1.4 GPa on the all-trans cells it expanded (a little too large, mostly
along b). Along the chain every all-trans cell is under -0.2 to -2.2 GPa: DFT wants the chains
straighter than the sampler left them. The helix cells are within 1.9 GPa either way.

The steering did what it was for: in every one of the 36 cells the labelled frames have the
DFT stress on both sides of zero along both lateral axes, so each cell's DFT minimum is inside
what is labelled. The frames span 52.5 to 74.0 A3 per monomer; 66 hold chlorine and 46 nitrogen.

With the PVDF sets and the DFT single points of 3.6 this is 330 frames. `mlip_train.py prepare`
takes them all (282 into the fit, 48 out, among them every frame of one composition, the
nitrile with TFE, kept whole as a composition the fit has not seen). The D3 add-on's far tail
against the DFT program's is one number per composition, 1.01 to 1.31 meV A3 per atom over the
16 compositions and the same within each to 1%: the add-on is right for chlorine and nitrogen.

Cost, measured: a 24-atom frame takes 2.3 minutes on this container's four cores and a
48-atom frame 6.5, so the set is about ten hours of DFT here and costs nothing. (Building it
took 70 minutes of the sampler.) A rented
32-core machine would do it in about an hour and a half for a few dollars; that is the PI's
call and is not needed.

The cutoffs are the anchor's 90/360 Ry, set for H, C and F. The chlorine file is ultrasoft and
its table value is 100/800 Ry, so a chlorine pair and a nitrile pair (one CFE, and one
acrylonitrile, in four monomers; all-trans and helix) were rerun at 100/800 before the labels
are used (`mlip_frames.py cutoff`, `docs/predictions/mlip_cutoff_check.json`). Done, and the
anchor's cutoffs hold: the forces change by 0.9 to 1.6 meV/A rms (8 at the largest), the stress
by 0.01 GPa, and the gap between the two cells of the same atoms by 0.002 kcal/mol per monomer
for chlorine and 0.001 for the nitrile. (Those gaps are between sampler geometries and are not
margins.)

The test of the result is not in this set: DFT single points on a candidate's own all-trans
and helix cells (the 72-atom `vdf5-rvdf1-6` pair first: about twenty minutes each here),
against the potential on the same geometries.

### 3.5 What separated the short fits from the 0.1 (items 1 to 3 done: see 3.6)

1. **Training time, which is a GPU.** A fit of ten to thirty thousand steps is hours to days
   on this container's four cores and under an hour on one GPU. Runpod's Secure Cloud had an
   RTX A6000 in stock at $0.69 an hour when this was written; three to four hours for a
   handful of fits is about $3 of the $32 left. Nothing is rented without the PI's go.
2. **Several fits, not one.** The same frames fitted from different starting seeds give a
   spread, and the spread is the error bar a second opinion needs: a margin the fits agree on
   to 0.03 can be used, one they scatter on by 0.2 cannot.
3. **Frames near the minima.** Strains of 1% and 4% beside the 2%, two axes at once, and
   shakes of 0.01 A: the bottom of each well, where the margins and the lattice are decided.
   Done: 54 frames, 1.2 hours of DFT, `docs/predictions/mlip/pvdf_near_frames_pbe_d3bj.extxyz`
   (`mlip_frames.py make --set near`). No fit has used them yet.
4. The other units' frames (3.4), which the container can label while a GPU fits.

### 3.6 Four long fits on a rented card (done, 2026-11-16; PI: "Go on runpod")

One A40 on Runpod's Secure Cloud at $0.59 an hour, 3.7 hours, about $2.18 (the pod removed
itself; its record could not be read back from the session afterwards, see `RESUME.md`).
`examples/gpu_jobs/mlip_fit.sh`, queue `examples/gpu_jobs/queue_20261116_mlip.txt`; results
`docs/predictions/mlip_own_pvdf_s1.json`, `..._s2`, `mlip_own_pvdfb_s1.json`, `..._s2`.

The fits: on from MACE-MP-0b3 as before, on the 130 PVDF frames (the 76 and the 54 near the
minima; 112 in the fit, 18 out), 400 epochs each, about 10,000 steps of the optimiser where the
first fit had 900. Two weightings of the stress and of the second stage's energy (100 and 1000;
1000 and 10000; forces 100), two seeds of each.

| test | DFT | stress 100, seed 1 | seed 2 | stress 1000, seed 1 | seed 2 |
|---|---:|---:|---:|---:|---:|
| held out: energy, meV per atom | | 0.34 | 0.31 | 0.28 | 0.21 |
| held out: force, meV/A (of 1190 rms) | | 21 | 19 | 22 | 20 |
| held out: stress, GPa | | 0.10 | 0.11 | 0.07 | 0.07 |
| 1. beta minus alpha on the DFT geometries, kcal/mol | +0.726 | +0.732 | +0.722 | +0.732 | +0.721 |
| 2. antipolar minus polar beta, meV | -1.11 | -0.62 | -0.71 | -1.22 | -0.52 |
| 3. C_11, GPa | 22.8 | 21.4 | 23.7 | 21.8 | 22.3 |
| 3. C_22 | 24.5 | 24.8 | 26.1 | 25.6 | 26.1 |
| 3. C_12 | 5.9 | 8.2 | 8.5 | 8.5 | 8.3 |
| 3. C_33 | 303 | 284 | 285 | 298 | 299 |
| 5. DFT force left on its relaxed alpha, beta and antipolar cells, eV/A | 0 | 0.009 | 0.013 | 0.005 | 0.006 |
| 5. alpha minus beta on its own cells: the potential less DFT, kcal/mol | | -0.007 | -0.083 | -0.006 | +0.005 |
| 5. antipolar minus polar beta on its own cells: the potential less DFT | | +0.011 | +0.009 | -0.003 | +0.014 |
| 6. turned-chain helix minus beta on its own cells: the potential less DFT | | -0.210 | -0.103 | -0.146 | -0.068 |
| 6. the same gap by DFT on each fit's own cell, meV | | -31.19 | -30.59 | -31.56 | -31.15 |

What the table says:

- **Where frames are near, the potential is the DFT.** On the helix, polar and antipolar cells
  every gap is within 0.014 kcal/mol of DFT, with the antipolar pair's sign right in all four
  (the short fits had it wrong), except one fit that let its helix cell shrink by 1% and is
  off by 0.08 there. C_11, C_22 and C_33 are within 7%; C_12 is 40% high in all four.
- **On a packing it was not shown, it is not.** The helix cell with one chain turned (test 6)
  comes out 3 to 9 meV per monomer too low in every fit: 0.07 to 0.21 kcal/mol, outside the 0.1
  in three of four, and always the same way. The four fits differ among themselves by 0.14 on
  it, so their spread does flag the cell as uncertain, but their mean is still 0.13 off: the
  spread is a warning and not an error bar.
- **A DFT single point on the cell the potential relaxed is right either way.** The four fits
  relaxed the turned-chain cell to four slightly different geometries; DFT's energy on them is
  the same to 1 meV per monomer (0.02 kcal/mol), and DFT finds 0.03 to 0.04 eV/A of force left.
  The error is in the potential's energy, not in the geometry it finds.
- The stress weighting matters less than the length of the fit. The heavier one has the lower
  stress and energy errors on held-out frames and the better C_33.

The relaxations of test 5 stop at 0.01 eV/A, and by that measure a cell under 0.2 GPa has
arrived, so three of the four fits left DFT's cells where they were. Relaxed to 0.002 eV/A
instead (the two heavier-weighted fits; `mlip_own_pvdfb_s1_tight.json`, `..._s2_tight`), the
potential's own cells are within 0.8% of DFT's in volume, DFT's energy on them is within 0.5
meV per monomer of its own minima (on the helix cell 0.1 meV below: the anchor's own
relaxation had that much left), and the gaps on them are the potential's to 0.02 kcal/mol.
The turned-chain cell is as before: 0.07 and 0.15 off by the potential, the same to 0.3 meV by
DFT.

All 32 DFT single points on cells the fits relaxed are in
`docs/predictions/mlip/pvdf_probe_frames_pbe_d3bj.extxyz`. None is in a fit yet.

**How it is to be used, from this.** Relax a cell in the potential, then ask DFT once for the
energy of the result: one self-consistent calculation in place of the dozens to hundreds of a
DFT relaxation. The potential's own energy is taken only for cells of a kind it has been
checked on, and every DFT single point made this way goes back into the frames.

The by-product stands, sharper: on the four fits' turned-chain cells DFT puts that packing
-0.1 to +0.9 meV per monomer from the alpha cell's own minimum. The polar packing of the helix
(delta-PVDF's kind) and alpha are level in this recipe to within 1 meV.

What a card costs for this, measured: 7 seconds an epoch for one fit of 112 frames on the A40;
four fits side by side took four times as long each (the card's time is shared between
processes), so they should be run one after another or on separate cards.

### 3.7 Next (not started)

1. **The other units.** The 168 frames of 3.4 are labelled, so: one fit on everything (the
   330 frames), on a second pod, one fit at a time (`examples/gpu_jobs/queue_20261117_mlip.txt`).
   An epoch will be longer than the 7 s measured (two and a half times the frames, twice the
   atoms in the large cells), so the first job there is again a timing. Waits on the PI's go,
   and on the session reaching Runpod's control API again.
2. **The test that matters.** A candidate's own all-trans and helix cells (the 72-atom
   `vdf5-rvdf1-6` pair first) relaxed in that potential, then one DFT single point each: the
   margin by DFT, against polyfind's and against the potential's own.
3. Each round's DFT single points go back into the frames. The turned-chain cell is the first
   case: with its six single points in the fit, is the next unseen packing still off?

### 3.8 A candidate's own cells, first pass: a caution, not yet a margin (2026-11-17)

PVDF with one reversed unit in six (`vdf5-rvdf1-6`, 72 atoms, H, C and F only; polyfind's
margin +0.064 kcal/mol per monomer with the pendants tilted, -0.113 with them rigid). Its
cells are `docs/predictions/polished_vdf5-rvdf1-6_{TTTT,TGTGp}_tilted.xyz`, the ones Dipole
was asked for single points on (E7 a). Numbers in
`docs/predictions/mlip_candidate_vdf5-rvdf1-6.json`; the four DFT single points in
`docs/predictions/mlip/candidate_frames_pbe_d3bj.extxyz`. A 72-atom single point takes 15 to
18 minutes on this container's four cores.

| | all-trans minus helix by DFT, kcal/mol | DFT force left, largest (rms), eV/A | DFT stress, largest, GPa |
|---|---:|---:|---:|
| polyfind's cells as they are | +0.164 | 2.26 (0.36) and 1.62 (0.29) | 0.93 |
| the same cells relaxed in the PVDF-only potential (`pvdfb_s2`) | -0.276 | 0.85 (0.12) and 0.63 (0.11) | 0.23 |
| polyfind | +0.064 | | |

What it says, and what it does not:

- **A DFT single point on polyfind's cell is not a test of polyfind's margin at the 0.1
  level.** DFT finds up to 2.3 eV/A on the atoms of polyfind's cells, mostly on the carbons
  (the bonds and angles polyfind holds rigid, taken from PVDF, in a chain that now has CF2-CF2
  and CH2-CH2 neighbours). Relaxing lowers DFT's energy by 43 meV per monomer in the all-trans
  cell and 24 in the helix: 1.0 and 0.55 kcal/mol, unequal by 0.44. That is how the answers to
  E7's first step (single points on the files as given) have to be read.
- **The -0.276 is not the DFT margin either.** The relaxed cells still carry 0.6 to 0.9 eV/A in
  DFT, and they were relaxed by a potential that was never shown a reversed unit. On the 56
  frames of the units set made of H, C and F that potential misses the energy between frames
  of the same atoms by 0.57 kcal/mol rms and the forces by 9.5% (1.7% on PVDF); its own margin
  for this cell is -1.0. It moved the all-trans cell a long way (a sheared lattice, the two
  chains' defects no longer level with their neighbours'), and DFT agrees that this is lower,
  but not that it is a minimum.
- So the sign of this candidate's margin in DFT is open: +0.16 unrelaxed, -0.28 half relaxed,
  polyfind +0.06. Settling it is what the fit on all 330 frames is for (reversed units are in
  54 of them), followed by the same two steps: relax, then one DFT single point.

**Relaxed on in DFT itself (done, 2026-11-17).** `mlip_frames.py relax` takes both cells on
from where the potential left them, on DFT's own forces and stress, about twenty minutes a step
here; every step is a frame (`docs/predictions/mlip/candidate_relaxation_frames_pbe_d3bj.extxyz`,
17 of them, 5.3 hours of DFT). Seven steps for the all-trans cell and ten for the helix reached
the 0.05 eV/A asked for (0.035 and 0.047 left; stress under 0.25 GPa): **all-trans minus helix
-0.261 kcal/mol per monomer.** The first steps:

| DFT step | all-trans minus helix, kcal/mol | largest DFT force left, eV/A (all-trans, helix) |
|---:|---:|---:|
| 0 (as the potential left them) | -0.276 | 0.85, 0.63 |
| 1 | -0.281 | 0.18, 0.21 |
| 2 | -0.271 | 0.13, 0.14 |
| 3 | -0.269 | 0.20, 0.16 |
| 4 | -0.268 | 0.06, 0.09 |
| 5 | -0.267 | 0.08, 0.08 |

The two energies fell together, by 0.2 meV per monomer a step toward the end, and the gap never
moved by more than 0.02 from its first value. Both cells kept their conformation (every
backbone torsion of the all-trans cell within a degree of 180; the helix still TGTG'). In all,
DFT's energy is 47 meV per monomer below polyfind's own all-trans cell and 29 below its helix
cell. **Between the two minima reached from polyfind's own cells, DFT puts all-trans 0.26
kcal/mol per monomer below the helix, where polyfind puts it 0.06 above.** A difference of
0.33, against a switching window 0.12 wide.

What that is and is not:

- It is the gap between two local minima, each reached downhill from polyfind's cell. The
  all-trans cell moved a long way on the way (its lattice sheared: the chains' defects are no
  longer level with their neighbours'); the helix cell did not. A lower helix arrangement that
  a downhill relaxation does not find would bring the gap back toward polyfind's. **Eight were
  tried and none is lower by more than 1.2 meV per monomer** (0.03 kcal/mol): the helix cell
  with its neighbours along b moved one, two or three monomers along the chain (-0.9, -1.2,
  +0.3, +0.2 and -0.4 meV by DFT against the helix cell at the same stage), along a (+33 and
  +32), and polyfind's other helix conformer, the reversed unit at the other place in the helix
  (+8.5; polyfind has it +5.2). Each was relaxed in the potential and given one DFT single
  point (`docs/predictions/mlip/candidate_registry_frames_pbe_d3bj.extxyz`). If the lowest
  holds, the margin is -0.23 in place of -0.26. Not tried: both lattice vectors moved at once,
  other setting angles of the chains, and anything on the all-trans side beyond the one
  arrangement found (a lower one there only widens the gap).
- It is one composition, one reversed unit in six, in PBE-D3(BJ). Against PVDF in the same
  recipe (+0.726) the reversed units lower the margin by 0.99 in DFT and by 0.68 in polyfind
  (+0.744 to +0.064): 5.9 against 4.1 kcal/mol per reversed unit, if it scales with their
  number as polyfind's own two cells do. Whether it does at one in fourteen is not known.
- If it holds, every cell whose margin is placed by a reversed unit sits lower than polyfind
  has it: by about 0.13 per monomer at one in fourteen, more than the window is wide.

**Where the 0.33 comes from: a freedom polyfind's cell does not have (2026-11-18).** In
polyfind's cell the chain axis stands upright on the lattice plane, so a chain's reversed unit
is level with the reversed units of its own images along a and b: a column of them, 5 A apart.
The cell DFT prefers is sheared along the chains, the reversed units no longer level. To
separate that from polyfind's energies, the all-trans cell was relaxed once more in the same
potential with the shear forbidden (`mlip_train.py relax --mask "1 1 1 0 0 1"`) and given one
DFT single point, to set beside the cells relaxed freely at the same stage:

| all-trans cell | all-trans minus helix by DFT, kcal/mol |
|---|---:|
| held upright, as polyfind's is | **+0.070** |
| polyfind | +0.064 |
| free to shear | -0.276 |

- **Polyfind's energy is right at this registry.** With the cell shaped as polyfind shapes it,
  DFT gives polyfind's margin to 0.006 kcal/mol. This is the first of polyfind's candidate
  margins to be met by DFT, and on a chemistry no anchor was fitted to. (As first written this
  point read "polyfind's energies are right". They are not along the registry: 3.9.)
- **Polyfind's cell is not.** Free to shear, the all-trans cell is 15 meV per monomer lower in
  DFT (0.35 kcal/mol; 2.1 per reversed unit), and DFT's shear stress on the upright cell says
  the same (0.5 and 0.3 GPa). The helix cell gains nothing from the same freedom (the five
  shifts along b above: within 1.2 meV). So the whole difference is on the all-trans side.
- **Which cell is the material?** In a chain made by polymerisation the reversed units fall at
  random, and at one in six or one in fourteen a reversed unit's neighbours on the next chains
  are nearly always ordinary units. The upright cell puts every reversed unit beside reversed
  units; the sheared one puts it beside ordinary ones. On that count the sheared cell is the
  nearer to the material and the upright one is the artefact of a small periodic cell. That is
  an argument, not a calculation: a cell with the defects placed at random has not been run.
- **What it touches.** The column is there in every polyfind cell with one pin per chain
  period, whatever the pin and however long the period: its cost per pin does not thin out as
  the period grows, so at one reversed unit in fourteen it is still about 0.15 per monomer if
  the 2.1 per unit carries over. For the other pins (CFE, CTFE, the nitrile, TrFE) it has not
  been measured and need not have the same sign.
- **What would settle it in polyfind:** two more lattice freedoms (the lattice vectors a and b
  allowed a component along the chain), then the candidates polished again. The number to meet:
  0.35 kcal/mol per monomer for this cell's all-trans state, nothing for its helix. **Done in
  3.9: the freedoms are there, and polyfind does not meet the number.**

These comparisons are between cells relaxed in the potential only (0.6 to 0.9 eV/A left in
DFT); relaxing on in DFT moved the sheared margin by 0.015.

This reads on the margin alone, which places a cell on the strain axis; it says nothing yet of
loss, rate or blocked force.

### 3.9 The shear in polyfind's cell, and what polyfind charges for registry (2026-11-18; PI: "Yes")

**The freedom is there and exact.** `polyfind.mixedcell.MixedModel(shear=True)` has two more
variables: what the lattice vectors a and b carry along the chains, in A. A chain's image one a
away then sits that much higher. Nothing else in the cell changes: the Ewald sum and the dipole
solve already took any lattice, and the chain-pair sums now move each image along the chains as
well as across. It is held to a known answer: a cell sheared by half a repeat is an upright
cell twice as wide, which the old code evaluates, and the two agree to 1e-14 kcal/mol in every
term, both forms of the sums, four cells and four shears. At no shear the arithmetic is the
upright cell's bit for bit, so no stored number moves (`tests/test_mixedcell.py`).

**Polyfind does not take it.** The all-trans cell of the reversed-unit candidate, against its
own upright cell, kcal/mol per monomer (`examples/shear_registry.py candidate`,
`docs/predictions/shear_registry_belt4.json`):

| | shear along a, b (monomers) | DFT | polyfind |
|---|---|---:|---:|
| each free to shear in its own energy | DFT 0.30, -0.90; polyfind 0.026, -0.078 | -0.35 | **-0.03** |
| at the cell DFT relaxed to, shear held, the rest relaxed in polyfind | 0.30, -0.90 | -0.35 | **+0.43** |

- Polyfind leaves its upright cell in DFT's direction (one of a to minus three of b, exactly, in
  both) and stops a tenth of the way.
- DFT's relaxed cell goes into the mixed cell's variables to 0.015 A rms over its 72 atoms
  (0.03 at worst; the setting angles are polyfind's to 0.01 degree), so polyfind can be asked
  about that very structure. As fitted it is 1.33 above polyfind's upright cell; with
  everything but the shear relaxed, 0.43 above. DFT has it 0.35 below.
- A start at DFT's shear with polyfind's own chains runs back to polyfind's minimum from 1.58
  above it. Released from DFT's own structure, the shear settles in a second minimum of
  polyfind's, 0.20 above the upright cell (0.01, -0.93 monomers).

**Why: moving chains along each other costs three times too much in polyfind.** Plain PVDF, no
reversed unit: two chains of one monomer, sheared rigidly, polyfind and a PBE-D3(BJ) single
point on the same twelve atoms (`shear_registry.py frames`, `mlip_frames.py label`; 45 s of DFT
each). Two packings: polyfind's beta cell, and the candidate's all-trans packing with plain
chains in it. kcal/mol per monomer against the unsheared cell:

| shear, monomers | beta: polyfind | DFT | ratio | candidate's packing: polyfind | DFT | ratio |
|---|---:|---:|---:|---:|---:|---:|
| 1/8 along a | 0.223 | 0.085 | 2.6 | 0.096 | 0.040 | 2.4 |
| 1/4 along a | 0.868 | 0.296 | 2.9 | 0.371 | 0.148 | 2.5 |
| 3/8 along a | 1.734 | 0.533 | 3.3 | 0.759 | 0.286 | 2.7 |
| 1/2 along a | 2.343 | 0.692 | 3.4 | 1.132 | 0.405 | 2.8 |
| 1/8 along b | 0.049 | 0.018 | 2.7 | 0.047 | 0.019 | 2.5 |
| 1/4 along b | 0.195 | 0.069 | 2.8 | 0.185 | 0.076 | 2.5 |
| 3/8 along b | 0.430 | 0.147 | 2.9 | 0.411 | 0.160 | 2.6 |
| 1/2 along b | 0.732 | 0.240 | 3.1 | 0.707 | 0.259 | 2.7 |
| the candidate's own shear in DFT (0.30, -0.91) | 2.176 | 0.675 | 3.2 | 1.575 | 0.525 | 3.0 |

Eight more cells along that last line give the same (2.6 to 3.5). And the test registered for
this question, E5 (`docs/E1_E2_PREDICTIONS.md`: one sublattice of beta slid, and one chain of
four, in Dipole's accepted beta cell), which had been waiting on Dipole, run here in the same
DFT. The forces on Dipole's unslid cell are zero in it, so it is the same DFT. kcal/mol per
monomer of the cell:

| E5 cell | polyfind, registered like for like | DFT | ratio |
|---|---:|---:|---:|
| one sublattice slid c/8 | 0.172 | 0.066 | 2.6 |
| c/4 | 0.656 | 0.230 | 2.9 |
| 3c/8 | 1.251 | 0.396 | 3.2 |
| c/2 | 1.533 | 0.467 | 3.3 |
| one chain of four slid c/4 | 0.671 | 0.232 | 2.9 |
| c/2 | 1.585 | 0.469 | 3.4 |

- **The registered acceptance is not met.** Sign, order along the slide and the anisotropy
  hold (one chain of four costs 2.01 times one sublattice per moved monomer in DFT, 2.07 in
  polyfind). The magnitudes, asked to be within a factor 1.5, are 2.6 to 3.4 times DFT's.
- **Over all 32 cells: 2.4 to 3.5, median 2.9.** The ratio grows with the displacement in
  every series, so the wall the chains ride over is too steep and not only too high.
- **It is the repulsion and dispersion terms.** Split by term, polyfind's electrostatics
  (charges and induced dipoles) lower the sheared cell slightly in every one of the 26 cells,
  by at most an eighth of what the pair terms raise it. The pair terms alone are 2.7 to 3.7
  times DFT's whole cost.
- **Its handles are hydrogen's contacts, and they are the lattice's handles too.** One fitted
  variable stepped at a time (`shear_registry.py handles`), the cost of half a monomer of shear
  in the beta cell moves with four: the hydrogen diameter (+2%: +13.5% along a, +14.3% along
  b), the H-F pair's own diameter (+2%: +6.9%, +13.1%), the depth of the belt wall against H
  (+0.1 in its logarithm: +10.0%, +7.8%) and the fluorine diameter (+2%: +3.6%, +6.4%). The
  well depths, the carbon diameter and the charges move it by 1% or less. Each of the four
  also pushes the cell apart (the lateral forces on the cell go from -0.06 and -0.28 to
  between -0.3 and -0.7 kcal/mol/A). The fit was given the lattice and its stiffness and no
  registry, and set these variables for those. Taken as linear, the four steps move the two
  costs and the two lateral forces in nearly the same proportions (condition number 8000 for
  four targets on four variables): no moderate change of them lowers the cost and holds the
  lattice. The hydrogen diameter alone reaches the factor 2.9 at 17% off, and leaves forces
  of +5 and +3 kcal/mol/A closing the cell. So the likely need is a new term, a wall softer
  than r^-12 for hydrogen's contacts, not new values for the old ones. That is an estimate
  from eleven steps, not a fit.
- **It accounts for the candidate.** At DFT's shear, plain chains cost 1.57 in polyfind and
  0.52 in DFT: 1.05 apart. On the candidate itself the two are 0.77 apart (+0.43 and -0.35).
  What the reversed units gain by leaving their column is of one size in both, about one
  kcal/mol per monomer; polyfind charges three times too much to get there and stays.

**What stands and what does not.**

- *Stands:* energies compared at one registry. The beta-alpha gap (+0.744 against +0.726),
  the candidate's margin with both cells upright (+0.064 against +0.070), and the stiffness
  across the chains (C_11, C_22, C_12 of 24.4, 25.3, 5.2 GPa against DFT's 22.8, 24.5, 5.9).
- *Does not stand until redone:* anything in which a chain, or part of one, moves along its
  neighbours. That is the chain-slide barrier (3.15 kcal/mol per monomer of the slid chain),
  every kink step behind a relaxation frequency, and the two shear constants along the chains,
  which are the second derivative of exactly this energy (2.4 to 2.7 times too stiff at an
  eighth of a monomer). How much each relaxed barrier falls has not been computed; the
  direction is that polyfind overstates them, so its rates are too slow.
- *The margins of cells with pins:* a pin's registry with the pins of the next chains is worth
  about as much as the window is wide, DFT moves it where polyfind cannot, and the two-chain
  cell offers one registry to choose. One row is measured (this candidate: +0.06 to -0.26).

**Polyfind's own re-polish with the shear free (done; `examples/shear_polish.py`,
`docs/predictions/sheared_margins_belt4.json`).** All 168 polished cells of the 35 rows, each
relaxed again from its stored upright cell with the two shears free, every torsion and angle
free as before, to a gradient under 1e-4. Pendant tilts were not redone: these are rigid
margins, to be read as changes.

- **The all-trans cells stay upright.** Mean gain 0.004 kcal/mol per monomer; four of 35 gain
  more than 0.01; the largest is the reversed-unit candidate's 0.029, the cell DFT lowers by
  0.35.
- **The helix cells move.** Mean gain 0.038, the largest 0.294; 38 of 133 gain more than 0.05,
  half of those at about half a monomer of shear along a or b.
- **So polyfind's margins rise:** by 0.044 on average (median 0.018, from -0.028 to +0.265).
  Nine rows move by more than half the window and four by more than its 0.12; the lowest
  helix conformer changes in 16 rows.

| row | rigid margin, upright | shear free | change |
|---|---:|---:|---:|
| the nitrile with reversed VDF on 21 mol% TrFE (`vdf9-trfe3-an1-rvdf1-14`) | -0.103 | +0.163 | +0.265 |
| 60/40 (`vdf6-trfe4`) | -0.369 | -0.186 | +0.183 |
| `vdf6-trfe3-cfe1` | +0.266 | +0.407 | +0.141 |
| `vdf8-trfe4-an1-rvdf1-14` | -0.276 | -0.151 | +0.124 |
| the nitrile alone (`vdf9-trfe4-an1-14`) | +0.033 | +0.148 | +0.115 |
| the CFE terpolymer, 64/29/7 (`vdf9-trfe4-cfe1-14`) | +0.140 | +0.251 | +0.111 |
| PVDF with reversed units (`vdf5-rvdf1-6`) | -0.113 | -0.141 | -0.028 |

- **These are not corrections toward DFT.** On the one row DFT has, polyfind's change is
  -0.03 and DFT's is -0.33, and it comes from the side polyfind does not move. A potential
  that charges three times too much for registry finds the sheared cells it can afford, not
  the ones the material would take. What the table does show is the size of the thing: in
  polyfind's own energy the registry a cell is allowed moves a margin by as much as the
  window is wide, in a third of the rows.
- **The ranking by margin is therefore not settled by this table or by the one before it.**
  It needs the registry cost put right first. And it reads on the margin alone, one of the
  four parts.

**On the four parts of the objective.** This is a finding about the model, not a candidate's
evaluation, and it lands on them unequally. Strain (the margin's place in the window): off by
up to several tenths for rows whose cells shear in DFT, the sign per row unknown. Tan delta and
frequency: both come from barriers of the kind above, overstated by an amount not yet computed.
Blocked force per area: unaffected where it is taken through the crystal across the chains,
which is the constant DFT confirms.

### 3.10 One fit on every frame, and two copolymer margins by DFT (2026-11-19; PI: "Create, all 558")

**The fit.** One fit, `all_s1`, on all 558 labelled frames (the ten files of
`docs/predictions/mlip/` as they stood: PVDF, the other units, the candidate's cells and its
relaxation steps, the registry and shear sets, the chains pulled apart, the lone chains). 485 in
the fit, 73 out, among them all eight frames of one composition (the nitrile with TFE in a TrFE
matrix), which the fit never saw. On from MACE-MP-0b3 as before, 240 epochs, stress 1000 and
energy 10000. One A40 on Runpod's Secure Cloud at $0.59 an hour: 2.85 hours with the timing
run, about $1.68; the pod was removed and a read of the account confirmed it. Queue
`examples/gpu_jobs/queue_20261118_mlip.txt`, results `docs/predictions/mlip_own_all_s1.json`.
The model is in the repository this time, `docs/predictions/mlip/models/all_s1.model`
(14.5 MB): the DFT work below starts from it.

| test | DFT | `all_s1` |
|---|---:|---:|
| held out: energy between frames of the same atoms, meV per monomer rms | | 1.7 |
| the same, PVDF groups / chlorine / the nitrile / the composition never seen | | 0.6 to 1.3 / 3.4 / 4.5 / 2.9 |
| held out: force, meV/A (of 1229 rms) | | 28.8 |
| held out: stress, GPa | | 0.07 |
| beta minus alpha on the DFT geometries, kcal/mol | +0.726 | +0.724 |
| the same, each cell relaxed in the potential | | +0.725 |
| antipolar minus polar beta, meV | -1.11 | -0.31 |
| C_11 / C_22 / C_12 / C_33 of beta, GPa | 22.8 / 24.5 / 5.9 / 303 | 23.0 / 27.3 / 7.5 / 295 |

The anchor's cells are among the frames, so the two margin rows are the fit reproducing what it
was shown. The needed 0.1 kcal/mol is 4.3 meV per monomer: the nitrile pair is just outside it
(0.104), every other held-out group inside.

**What it was not shown: the margin of a TrFE cell.** The project held two margins from DFT,
PVDF and the reversed-unit candidate. Neither has a TrFE unit, and the leading rows are TrFE
matrices with a pin. `vdf3-trfe1-4` is the smallest cell that holds one (VVRV, 25 mol% TrFE,
48 atoms in two chains). Polyfind's cells for it, under `pvdf-dft-chain1` with the shear free
(`examples/rows_under_preset.py`), were written out as atoms
(`examples/row_cells_frames.py write`), relaxed in `all_s1` (`mlip_train.py relax`), and given
one DFT single point each (`mlip_frames.py label`); the two labelled cells of the same sequence
that the sampler had built for the frames went the same way.

| cell | DFT as polyfind has it, eV per monomer (force, eV/A) | relaxed in the potential (force) |
|---|---:|---:|
| all-trans | -2039.78865 (1.31) | -2039.81245 (0.07) |
| helix TG'TG, free to shear | -2039.80765 (1.27) | -2039.84186 (0.08) |
| the same helix held upright | -2039.79762 (1.38) | -2039.83433 (0.05) |
| the sampler's all-trans arrangement | -2039.77739 (0.56) | -2039.80482 (0.07) |
| the sampler's helix arrangement | -2039.81395 (0.33) | -2039.83920 (0.04) |

| all-trans minus helix, kcal/mol per monomer | DFT | chain1 | belt4 | `all_s1` alone |
|---|---:|---:|---:|---:|
| shear free | **+0.678** | +0.671 | +0.641 | +0.606 |
| helix held upright | +0.505 | +0.421 | +0.577 | +0.435 |
| the helix's gain from shearing | -0.174 | -0.250 | -0.064 | -0.171 |
| the sampler's arrangement of the same chains | +0.793 | | | +0.696 |

- **Both presets meet the first DFT margin of a TrFE cell**, by 0.007 (chain1) and 0.037
  (belt4), on a chemistry neither was fitted to. They split it differently: belt4 has the
  upright cells 0.07 too far apart and gives the helix a third of its gain from shearing;
  chain1 has them 0.08 too close and gives it half as much again.
- **Polyfind's arrangement of the chains is the lower one in DFT**, in both conformations (by
  7.6 meV per monomer all-trans, 2.7 helix, against the sampler's). The search is doing its job.
- **A DFT single point on polyfind's cell as it stands is again not a margin** (+0.438 with
  1.3 eV/A on the atoms, +0.678 relaxed), as on the candidate (3.8).
- **The potential finds the geometry; it is not the referee.** One relaxation in it leaves
  0.04 to 0.08 eV/A in DFT. Its own margins are 0.07 to 0.10 below DFT's on these three pairs:
  at the edge of the needed 0.1, on a composition it has 20 frames of. Its energy change along
  a relaxation is DFT's to 0.1 to 0.3 meV per monomer on three cells and to 1.0 and 2.9 on the
  other two.

**A chlorine cell: neither preset has it.** `vdf2-trfe1-cfe1-4` (VRVC, one CFE and one TrFE in
four, 48 atoms) is the cell above with a CFE in place of one VDF. The two presets disagree on
it by 0.68, so both presets' all-trans cell and lowest helix cell (the same placement, G'TGT,
in both) were relaxed in the potential; they arrive at the same two cells (0.06 meV apart), and
DFT was asked for those (0.06 eV/A left on each).

| all-trans minus helix, kcal/mol per monomer | DFT | chain1 | belt4 | `all_s1` alone |
|---|---:|---:|---:|---:|
| VRVC, shear free | **+0.682** | +0.381 | +1.057 | +0.768 |
| the CFE's effect: VRVC less VVRV | +0.004 | -0.289 | +0.416 | +0.162 |

- **chain1 is 0.30 low and belt4 0.38 high.** In this cell a CFE in place of a VDF leaves DFT's
  margin where it was. Belt4 makes it a strong pin for the helix; chain1 makes it a pin for
  all-trans. Neither preset's chlorine terms had a DFT energy of a chlorine cell among their
  targets; chain1's wall and charges were refitted on H, C and F only.
- **This is one cell, a dense one**: a CFE on every chain, one in four, where the leading rows
  have one in fourteen or twenty. In polyfind the unit's effect per unit is not the same at
  the two densities (belt4: +1.7 per unit here, +3.7 in `vdf9-trfe4-cfe1-14`; chain1: -1.2 and
  +0.5), so the +0.004 does not carry over to the terpolymer cell as a number. What carries
  over is that the presets' pins for chlorine are untested where they are right and wrong by
  0.3 to 0.4 where they have been tested.
- The potential alone is 0.086 high here, on a chlorine-to-chlorine contact between chains
  that none of its frames has.

**Where the chlorine error sits: inside the chain, in both presets.** One more DFT single
point per relaxed cell, on a chain of it alone in a 16 A box (`examples/chain_split.py`, as in
`docs/CRYSTAL_STAGE.md` 5.29), splits each margin into what is inside a chain and what is
between chains. Polyfind's own cells are split the same way (upright cells).

| all-trans minus helix | cell | inside a chain | between chains |
|---|---:|---:|---:|
| VVRV, helix upright: DFT | +0.505 | +2.935 | -2.430 |
| chain1 | +0.421 | +3.069 | -2.648 |
| belt4 | +0.577 | +3.664 | -3.087 |
| VRVC: DFT | +0.682 | +2.979 | -2.298 |
| chain1 | +0.408 | +2.864 | -2.457 |
| belt4 | +1.060 | +3.992 | -2.932 |
| the CFE's effect, VRVC less the upright VVRV: DFT | +0.177 | +0.044 | +0.133 |
| chain1 | -0.013 | -0.204 | +0.191 |
| belt4 | +0.483 | +0.328 | +0.155 |

(The partner here is the VVRV helix held upright, the cell polyfind's split is made on; the
VRVC cells do not shear in either preset or in the potential. Against the VVRV helix free to
shear, which is 0.174 lower, the unit's effect on the cell is the +0.004 of the table before.)

- **Between chains both presets have the CFE unit right**, to 0.06 (chain1) and 0.02 (belt4).
  The wall and the charges, where chlorine meets the next chain, are not the problem.
- **Inside the chain chain1 puts the unit 0.25 too low and belt4 0.28 too high**, per monomer
  of a four-monomer cell: about -1.0 and +1.1 kcal/mol per CFE unit, all-trans against helix.
  That is the bonds at the CFCl carbon: their torsion term is the shared triple that was
  fitted to CF2-CH2 bonds, and nothing chlorine-specific. It is the kind of error 5.29 found
  for the reversed unit, and the remedy there was a torsion triple of the unit's own bonds,
  fitted to DFT on lone chains (five minutes a single point).
- **If it is the unit's own bonds, the error goes with the unit and not with the density**: one
  CFE in fourteen would carry a fourteenth of it, -0.07 per monomer in chain1 and +0.08 in
  belt4. That would put the terpolymer cell's pin near +0.11 (chain1's +0.04 corrected) or
  +0.18 (belt4's +0.26 corrected), in place of two numbers a factor of seven apart. An
  inference from one cell, to be checked by the fit.
- **On the TrFE cell alone belt4 is right by cancellation again**: +0.73 inside the chain and
  -0.66 between chains, the same two errors as on the reversed-unit candidate. Chain1's are
  +0.13 and -0.22. So chain1 is the preset to build the chlorine term on.

**What this reads on.** The margin, which places a cell on the strain axis, and nothing else.
The loss, the rate and the blocked stress of these cells are not computed here. For the ranking
it says: rows that differ by a TrFE matrix alone can be compared in either preset; rows whose
margin is set by a chlorine unit (CFE, CTFE, CDFE: 15 of the 35) cannot, in either, until the
chlorine terms are fitted to DFT margins of cells like this one.

Frames: `docs/predictions/mlip/row_cells_vdf3-trfe1-4_frames_pbe_d3bj.extxyz` (10),
`row_cells_vdf2-trfe1-cfe1-4_frames_pbe_d3bj.extxyz` (2). Tables:
`docs/predictions/row_cells_dft_vdf3-trfe1-4.json`, `row_cells_dft_vdf2-trfe1-cfe1-4.json`
(`examples/row_cells_frames.py collect`); polyfind's side `margins_chain1_dft_rows.json`,
`margins_belt4_dft_rows.json`. The split: five lone-chain frames
`docs/predictions/mlip/row_cells_lone_chain_frames_pbe_d3bj.extxyz`, DFT's
`chain_split_dft_rows_pbe_d3bj.json`, polyfind's `chain_split_polyfind_dft_rows.json`. Tests: `tests/test_row_cells_dft.py`,
`tests/test_mlip_own.py::test_the_fit_on_every_frame_holds_the_anchor_and_a_composition_it_never_saw`.

### 3.11 A torsion term for the bonds at a CFCl carbon, and the CTFE cell (2026-11-19)

3.10 put chain1's chlorine error inside the chain, at the bonds of the CFCl carbon, which take
the triple fitted to CF2-CH2 bonds. This section gives those bonds a triple of their own, fitted
to DFT on lone chains, and tests it on two cells that were not in the fit.

**The lone chains.** The polymer with a CFE at every second unit, twelve atoms a repeat, alone
in a 16 A box, settled in DFT as the reversed unit's chains were (`docs/CRYSTAL_STAGE.md` 5.29;
`examples/chain_conformers.py`, chain `cfe2`): all-trans, and the helix. The chain reads the same
both ways, so the helix has two placements that differ: the CFCl-CH2 bond gauche one way or the
other, the chlorine on opposite sides of the fold. Forces left: 0.005 eV/A.

| against the all-trans chain, kcal/mol per monomer | DFT | chain1 | with the triple |
|---|---:|---:|---:|
| helix, CFCl-CH2 bond at -60 degrees | -2.475 | -1.783 | -2.403 |
| helix, the other placement (+50 degrees) | -3.189 | -2.835 | -3.293 |
| the second less the first | -0.714 | -1.052 | -0.890 |

- **chain1 has both helices too high**, by 0.69 and 0.35 per monomer: the gauche well of a CFCl
  bond is too shallow. That is the sign and the size 3.10 found in the cell (about 1 kcal/mol
  per CFE unit).
- **The triple** (`examples/unit_torsion_fit.py`, 25 s): (V1, V2, V3) = (+0.440, -2.180, +1.501)
  kcal/mol where the shared one is (+0.450, -0.522, +0.746). It lowers a gauche CFCl bond by
  about 1.25 kcal/mol against a trans one and leaves both rows within 0.11.
- **What the two rows hold and what they do not.** A torsion term is even in the angle, so both
  placements see the same one: the rows fix one combination of the three coefficients, the
  depth of the gauche well against the trans one. A ridge holds the rest to chain1's triple.
  **The barrier between the wells is not fitted.** It matters to how fast a kink passes a CFE
  unit, which is the rate axis, and it needs a DFT point on the barrier before it is trusted.
- What separates the two placements (0.71 in DFT) is not the torsion term but chlorine against
  fluorine across the fold. chain1 has it at 1.05; with the triple 0.89.
- A mistake, corrected: two of the four conformers I first built for this chain had two gauche
  bonds of opposite sense next to each other. That is not a conformer (the groups two bonds on
  meet); the chains broke up in the relaxation and the frames are discarded.

**The test: cells that were not in the fit.** The preset is `pvdf-dft-chain1-cfe`: chain1 with
that triple on the CH2-CFCl bonds and nothing else changed.

| all-trans minus helix, kcal/mol per monomer | DFT | chain1 | chain1-cfe | belt4 |
|---|---:|---:|---:|---:|
| VVRV (TrFE only) | +0.678 | +0.671 | +0.671 | +0.641 |
| VRVC (one CFE in four) | **+0.682** | +0.381 | **+0.690** | +1.057 |
| VRVK (one CTFE in four) | **+0.315** | +0.134 | **+0.400** | +0.734 |

- **The CFE cell passes**: 0.008 from DFT, where chain1 was 0.30 low, and in DFT's helix
  placement (G'TGT).
- **The CTFE cell passes, at the edge**: 0.085 high against a needed 0.1. CTFE is a new DFT
  margin here (`examples/unit_dft_check.py`, the whole check for one small row): +0.315, less
  than half the CFE cell's. Neither chain1 (0.18 low) nor belt4 (0.42 high) had it. The unit
  brings two bonds of its own: a CFCl-CH2 bond, which takes the new triple, and a CF2-CFCl bond,
  which has none and is trans in both cells, so this cell does not test it.
- A cell with no CFCl carbon is chain1's to the last digit.

The split into inside a chain and between chains, as in 3.10 (polyfind's cells upright, the same
helix placement as DFT's):

| all-trans minus helix | cell | inside a chain | between chains |
|---|---:|---:|---:|
| VRVC: DFT | +0.682 | +2.979 | -2.298 |
| chain1-cfe | +0.715 | +3.161 | -2.445 |
| VRVK: DFT | +0.315 | +1.091 | -0.776 |
| chain1 | +0.152 | +1.075 | -0.923 |
| chain1-cfe | +0.468 | +1.394 | -0.927 |
| belt4 | +0.861 | +2.346 | -1.484 |

- On the CFE cell the two parts are now off by +0.18 and -0.15, which are the two errors chain1
  has on the TrFE cell (+0.13 and -0.22). What the CFE unit adds to them is 0.05 inside the
  chain and 0.07 between chains.
- On the CTFE cell the inside of the chain is 0.30 high, more than the TrFE cell's 0.13, and
  the margin lands inside the 0.1 with help from the error between chains. That is a pass on
  the number that was asked for and not a clean one. CTFE rows carry that caution.
- The CTFE unit folds the chain and parts the chains far less than the CFE unit does in DFT
  (+1.09 and -0.78 against +2.98 and -2.30): its helix cell is denser than the others.

**What this reads on.** The margin, which places a cell on the strain axis, and nothing else.
The loss, the rate and the blocked stress are not computed here, and the rate of a chlorine row
will lean on the CFCl bond's barrier, which is not fitted. For the ranking: rows whose units are
VDF, TrFE, reversed VDF, CFE and CTFE (20 of the 35) now have a preset that meets DFT on a small
cell of each unit; the other 15 hold a unit that has not had its check.

Files: `docs/predictions/chain_conformers_cfe_pbe_d3bj.json` (the lone chains),
`unit_torsion_cfe_chain1.json` (the fit), `margins_chain1-cfe_dft_rows.json` (the cells under the
preset), `row_cells_dft_vdf2-trfe1-ctfe1-4.json` and `chain_split_dft_vdf2-trfe1-ctfe1-4.json`
(the CTFE cell in DFT), `mlip/row_cells_vdf2-trfe1-ctfe1-4_frames_pbe_d3bj.extxyz` (2 frames),
`chain_split_polyfind_dft_rows.json` (polyfind's splits, set `chain1-cfe`). Tests:
`tests/test_row_cells_dft.py`.

### 3.12 The nitrile: all-trans below every helix in DFT, both presets wrong, and no simple term fixes it (2026-11-19)

The same check for acrylonitrile (`examples/unit_dft_check.py --row vdf2-trfe1-an1-4`: VRVA, one
nitrile and one TrFE in four, fifty atoms): chain1's cells relaxed in our potential, then one DFT
single point on the all-trans cell and on each of the three helix placements the potential put
within 0.15 of its lowest (0.10 to 0.15 eV/A left).

| all-trans minus helix, kcal/mol per monomer | DFT | chain1 | belt4 | `all_s1` alone |
|---|---:|---:|---:|---:|
| helix TGTG' | **-0.386** | +0.258 | -0.224 | -0.466 |
| helix TG'TG | -0.404 | -0.776 | -0.818 | -0.352 |
| helix G'TGT | -0.678 | -0.360 | +0.354 | -0.383 |
| the margin (against the lowest helix) | **-0.386** | +0.258 | +0.354 | -0.352 |

- **In DFT the all-trans cell is the lowest, by 0.39 over the best helix placement.** With a
  VDF in the nitrile's place (VVRV, 3.10) the margin is +0.678: one nitrile in four moves it by
  -1.06 per monomer, -4.3 kcal/mol per nitrile unit. In this cell the nitrile pins all-trans.
- **Both presets have a helix cell lowest**, each on a different placement, and are off on single
  placements by up to 0.64 (chain1) and 1.03 (belt4). Seven of the 35 rows hold a nitrile.
- The potential, which was never shown this cell, has the sign and the margin to 0.03; on single
  placements it is off by 0.05 to 0.30.
- **This is one cell, a dense one**, as the CFE cell was: a nitrile on every chain, one in four.
  The rows have one in ten or fourteen. Whether the pin holds there is not known from this.

Split by lone chains cut from each cell (one more DFT point a cell), against polyfind's cells
held upright as in 3.10:

| all-trans minus helix | cell | inside a chain | between chains |
|---|---:|---:|---:|
| TGTG': DFT | -0.386 | +2.128 | -2.514 |
| chain1 | +0.060 | +1.984 | -1.924 |
| belt4 | -0.298 | +1.852 | -2.149 |
| TG'TG: DFT | -0.404 | +0.212 | -0.616 |
| chain1 | -0.790 | -0.240 | -0.551 |
| belt4 | -0.937 | -0.260 | -0.677 |
| G'TGT: DFT | -0.678 | +1.090 | -1.767 |
| chain1 | -0.355 | +0.748 | -1.103 |
| belt4 | +0.344 | +1.935 | -1.591 |

- **chain1 is wrong on both parts.** Between chains it is 0.59 and 0.66 short of DFT's binding of
  the all-trans cell on two placements (0.07 on the third). Inside the chain it is 0.14 to 0.45
  low. The CFE unit's error was inside the chain only.
- DFT's cells here are sheared along the chains (the TGTG' cell by 15 degrees) and polyfind's
  split holds its cells upright. With polyfind's own shear-free cells
  (`chain_split.py polyfind --states`) chain1's errors between chains are +0.81, +0.04 and
  +0.70, and inside -0.16, -0.41, -0.38. The fits below use the shear-free cells.
- **Why between chains.** Nitrogen's wall is at its UFF values in every preset (3.66 A, 0.069
  kcal/mol): the walls of H, C and F were refitted for chain1 and nitrogen's was not. The charge
  the C-N bond moves is 0.37 e, a bond dipole of 2.1 D, and it was never a fitting target.

The lone chains of the alternating VDF-AN polymer in DFT (`chain_conformers.py`, chain `an2`,
0.01 eV/A left), which are what a torsion term for the nitrile's bonds is fitted to:

| against the all-trans chain, kcal/mol per monomer | DFT | chain1 |
|---|---:|---:|
| helix, one placement (the CH-CH2 bond at -83 degrees) | -1.569 | -0.535 |
| helix, the other (+56 degrees) | +0.367 | +1.286 |

chain1 has both about 1.0 too high and their difference right (1.82 against 1.94).

**A first fit, and why it is not the fix.** Fitting the C-N increment and nitrogen's wall to
the three energies of the table (`examples/unit_nonbonded_fit.py`; 0.747 e, 3.70 A, 0.043)
brings them to -0.380, -0.358 and -0.716, within 0.05 of DFT. Split, the lowest placement is
then +0.50 too high inside the chain and 0.59 too low between chains: two errors that cancel,
which is how belt4 was right on TrFE and wrong elsewhere. On the lone chains that charge takes
the two placements 3.16 apart where DFT has 1.94. It is kept as
`docs/predictions/unit_nonbonded_an_margins_chain1.json` and not used.

**A second fit, in parts, with the energies held back as the test: it fails.** The charge and
nitrogen's wall were fitted to the between-chains parts alone (`--between`, polyfind's cells
shear-free), then a torsion triple for the nitrile's bonds to the lone chains on top of that.

| | charge and wall | with the torsion triple | DFT |
|---|---:|---:|---:|
| C-N increment; nitrogen's wall | 0.632 e; 3.445 A, 0.052 | the same | |
| between chains, three placements (fitted) | -2.815, -0.653, -1.595 | | -2.514, -0.616, -1.767 |
| lone chains, two placements | -1.182, +1.626 | -1.598, +0.379 (fitted) | -1.569, +0.367 |
| VRVA energies, three placements (**test**) | -0.185, -0.423, -0.597 | -0.018, +0.046, -0.237 | -0.386, -0.404, -0.678 |
| VRVA margin (**test**) | -0.185 | +0.046 | -0.386 |
| VVVA margin (**test**, the second cell) | -0.549 | +0.164 | -0.091 |

- The between-chains parts stay 0.30 and 0.17 off: as the charge moves the cells change basin,
  so one placement jumps by 0.19 for a change of 0.002 e, and a least-squares search has nothing
  smooth to follow.
- With the charge at 0.63 e the lone chains' two placements are 2.81 apart (DFT 1.94; chain1's
  0.37 e has 1.82). The torsion triple that then meets both rows has coefficients of 3 kcal/mol
  and does it through the 25 degrees between the two placements' angles. That is a fit of two
  numbers, not a description of the bond.
- **Both tests fail.** VRVA's margin is 0.20 high with the charge and wall alone and 0.43 high
  with the torsion; VVVA's is 0.46 low and 0.26 high. Neither is registered.
- A nitrogen-hydrogen contact term of its own (the kind chain1 needed for fluorine-hydrogen),
  tried at three settings (3.0 and 2.8 A, 0.04 and 0.08 kcal/mol; a scan, not stored), moves the
  between-chains parts by 0.15 at most where 0.7 to 0.8 is missing, and the G'TGT placement by
  0.02. Of the terms tried only the charge moves them.
- **The second nitrile cell in DFT** (VVVA, `vdf3-an1-4`, 0.19 eV/A left): all-trans minus helix
  **-0.091** = +0.878 inside a chain - 0.969 between chains. All-trans is the lower again, by
  less. chain1 has -0.292 on another helix placement and -0.751 on DFT's; belt4 +0.058. chain1's
  errors on this cell (-0.53 inside, -0.15 between on DFT's placement) are not its errors on
  VRVA (-0.16 to -0.41 inside, +0.04 to +0.81 between): there is no one correction to find
  with this form.
- **So polyfind does not have the nitrile, in any preset, and a bond charge with a wall and a
  torsion term is not enough to give it one.** What the form lacks is most likely the charge
  the nitrile draws from the carbon it sits on and that carbon's hydrogen, which one bond
  increment cannot move. That is a refit of the charge model on nitrile cells and lone chains,
  not a term added on.
- The potential has both nitrile cells' margins to 0.03 and 0.08 without having been shown
  either. For the seven nitrile rows it is the best evidence here (on the seven small cells of
  3.13 it is within 0.11 of DFT on six and 0.36 off on the TFE cell).
- **What it does to the ranking.** The leader of `docs/CANDIDATE_MATERIALS.md` 5.6 is a nitrile
  row (`vdf9-trfe3-an1-rvdf1-14`, belt4 +0.037). In belt4's fourteen-monomer cells the nitrile
  is a pin for the helix, +2.0 kcal/mol per unit (`vdf9-trfe4-an1-14` +0.197 against the
  matrix's +0.056); in DFT's four-monomer cells it is a pin for all-trans, -3.3 to -4.3 per
  unit. If DFT's number held at one in fourteen the leader would sit near -0.35. It does not
  hold: the row's cells relaxed in the potential at its own dilution give -0.055 (3.14).

**What this reads on.** The margin only: the strain axis. Nothing here is a loss, a rate or a
blocked stress.

Files: `docs/predictions/row_cells_dft_vdf2-trfe1-an1-4.json`,
`chain_split_dft_vdf2-trfe1-an1-4.json`, `mlip/row_cells_vdf2-trfe1-an1-4_frames_pbe_d3bj.extxyz`
(4 frames), `chain_conformers_an_pbe_d3bj.json`, `chain_split_polyfind_dft_rows.json` (sets
`chain1`, `belt4`, `an-fitted-to-margins`); the second cell `row_cells_dft_vdf3-an1-4.json`,
`chain_split_dft_vdf3-an1-4.json`, `mlip/row_cells_vdf3-an1-4_frames_pbe_d3bj.extxyz`; the fits
`unit_nonbonded_an_margins_chain1.json`, `unit_nonbonded_an_between_chain1.json`,
`unit_torsion_an_chain1.json` and their tests `margins_an_between_dft_rows.json`,
`margins_an_between_torsion_dft_rows.json`. Tests: `tests/test_row_cells_dft.py`.

### 3.13 The CDFE and TFE cells, where each unit stands, and 25 rows' margins under chain1-cfe (2026-11-19)

**CDFE.** The same check (`vdf2-trfe1-cdfe1-4`, VRVD; one of its four helix placements does not
build at its ideal torsions, and the polish now says so and goes on):

| all-trans minus helix | cell | inside a chain | between chains |
|---|---:|---:|---:|
| DFT (0.22 eV/A left) | **+0.968** | +2.394 | -1.426 |
| chain1, its cells shear-free | +0.975 | +3.091 | -2.117 |
| belt4 | +1.720 | | |

- **chain1 has the margin to 0.007 and the two parts 0.70 and 0.69 off, opposite ways.** The
  all-trans cell changes shape when it is relaxed in the potential (gamma from 94 to 105
  degrees), so DFT's split and polyfind's are of cells that differ, and the parts of a sum need
  not agree cell for cell. The lone chains settle it, because they have no cell:

| VDF-CDFE chain alone, against all-trans, kcal/mol per monomer | DFT | chain1 |
|---|---:|---:|
| helix, the CF2-CHCl bond at -62 degrees | -2.843 | -3.418 |
| helix, the other placement (+46 degrees) | +3.750 | +1.467 |

- chain1 has the low helix 0.58 too low and the high one 2.28 too low: the inside of the chain
  is wrong by itself, in the direction the cell's split shows. **So the margin of the cell is
  met by two errors that cancel, and CDFE is not counted as passed.** A torsion term for the
  CF2-CHCl bond would take about 0.3 off the cell's margin and leave it 0.3 low, with the error
  between chains still there. Like the nitrile, the unit needs both parts.

**Where the units stand** (one small cell each, all-trans minus helix against DFT):

| unit | rows of the 35 that hold it | DFT | chain1-cfe | reading |
|---|---:|---:|---:|---|
| TrFE (VVRV) | 30 (6 with nothing else) | +0.678 | +0.671 | passes |
| reversed VDF (the candidate of 3.8) | 7 | -0.261 | -0.199 | passes (`docs/CRYSTAL_STAGE.md` 5.29) |
| CFE (VRVC) | 9 | +0.682 | +0.690 | passes, with the term of 3.11 |
| CTFE (VRVK) | 3 | +0.315 | +0.400 | passes at the edge; its CF2-CFCl bond is not tested |
| CDFE (VRVD) | 4 | +0.968 | +0.975 | the margin by cancellation: not passed |
| acrylonitrile (VRVA, VVVA) | 7 | -0.386, -0.091 | +0.258, -0.292 | fails; no term of this form fixes it (3.12) |
| TFE (VRVF) | 2 | -0.076 | -0.060 | passes: the margin to 0.02, its parts 0.26 and 0.29 off, 0.1 more than on the TrFE cell |
| ethylene (VRVE) | 2 | -0.146 | -0.232 | passes at the edge (added 16:50 UTC): its parts +0.04 and -0.23 off; belt4 +0.500, the wrong sign; the potential -0.082 |
| VF (VRVY) | 1 | +0.290 | +0.479 | fails by 0.19 (added 20:40 UTC); belt4 +0.452, 0.16 high; the potential +0.274. DFT leaves 0.32 and 0.35 eV/A on the two cells, more than on any other pair, so its number is the least settled of the set |
| DFE (VRVX) | 1 | +0.961 | +0.955 | passes, and cleanly (added 21:10 UTC): its parts +0.04 and -0.07 off; belt4 +1.234, 0.27 high; **the potential +0.309, 0.65 low, its worst cell** (DFT leaves 0.29 eV/A on both cells) |

**TFE** (`vdf2-trfe1-tfe1-4`, VRVF; 0.07 eV/A left): DFT -0.076 = +0.895 inside a chain - 0.971
between chains, all-trans the lower by a little. chain1 -0.060 with the shear free; upright on the
same helix -0.115 = +1.150 - 1.264. Belt4 +0.224. chain1's two parts are off by +0.26 and -0.29,
which is its error on the TrFE cell (+0.13, -0.22) and about 0.1 more each way. It is counted as
passing, on the same standing as TrFE itself.

**The potential's own margins on the seven cells**, since it is what would be used where polyfind
has no preset:

| cell | DFT | `all_s1` alone | off |
|---|---:|---:|---:|
| TrFE (VVRV) | +0.678 | +0.606 | -0.07 |
| CFE (VRVC) | +0.682 | +0.768 | +0.09 |
| CTFE (VRVK) | +0.315 | +0.205 | -0.11 |
| nitrile (VRVA) | -0.386 | -0.466 | -0.08 |
| nitrile (VVVA) | -0.091 | -0.008 | +0.08 |
| CDFE (VRVD) | +0.968 | +0.935 | -0.03 |
| TFE (VRVF) | -0.076 | -0.435 | **-0.36** |

It has the sign of all seven and is within 0.11 on six. On the TFE cell it is 0.36 off: its own
margin is an indication and not a referee, as 3.6 found on a packing it had not been shown. (For
the nitrile rows of 3.12 that still makes it the better evidence: both presets have the wrong sign
there.)

**The 25 rows that hold none of the unchecked units and no nitrile**, each polished cell relaxed
again under `pvdf-dft-chain1-cfe` with the shear free and the pendants tilted
(`examples/rows_under_preset.py`; `docs/predictions/margins_chain1-cfe.json`). Margin is
all-trans minus the lowest helix placement, kcal/mol per monomer; the strain is the all-trans
cell's length along the chain over the helix cell's. Eight of the rows that have no CFCl carbon
are chain1's own run carried over (the two presets are the same potential for them).

| row | monomers | margin | cell held upright | strain along the chain | belt4 (upright) | its units |
|---|---:|---:|---:|---:|---:|---|
| `vdf5-rvdf1-6` | 6 | -0.199 | +0.042 | +9.0% | +0.064 | checked |
| `vdf9-trfe4-rvdf1-14` | 14 | -0.164 | -0.176 | +10.8% | -0.162 | checked |
| `vdf7-trfe2-cdfe1` | 10 | -0.095 | -0.130 | +11.6% | +0.392 | CDFE: not passed |
| `vdf8-trfe4-cfe1-rvdf1-14` | 14 | -0.070 | -0.060 | +11.7% | -0.044 | CFE: fitted term |
| `vdf12-rvdf2-14` | 14 | -0.030 | +0.159 | +9.3% | +0.162 | checked |
| `vdf8-trfe4-ctfe1-rvdf1-14` | 14 | +0.061 | +0.130 | +8.7% | +0.245 | CTFE: at the edge |
| `vdf6-trfe4` | 10 | +0.083 | -0.044 | +11.9% | -0.199 | checked |
| `vdf9-trfe4-cdfe1-14` | 14 | +0.089 | +0.068 | +11.3% | +0.685 | CDFE: not passed |
| `vdf9-trfe4-tfe1-14` | 14 | +0.147 | +0.062 | +12.2% | -0.026 | TFE: checked |
| `vdf10-trfe4-14` | 14 | +0.222 | +0.131 | +11.0% | +0.056 | checked |
| `vdf9-trfe5-14` | 14 | +0.232 | +0.149 | +11.6% | +0.083 | checked |
| `vdf9-trfe4-ctfe1-14` | 14 | +0.286 | +0.126 | +10.4% | +0.300 | CTFE: at the edge |
| `vdf13-trfe6-cfe1-20` | 20 | +0.307 | +0.199 | +11.0% | +0.258 | CFE: fitted term |
| `vdf7-trfe3` | 10 | +0.354 | +0.222 | +10.9% | +0.213 | checked |
| `vdf9-trfe4-cfe1-14` | 14 | +0.359 | +0.231 | +10.6% | +0.317 | CFE: fitted term |
| `vdf8-trfe2` | 10 | +0.364 | +0.219 | +10.8% | +0.229 | checked |
| `vdf6-trfe3-cdfe1` | 10 | +0.433 | +0.428 | +11.1% | +1.012 | CDFE: not passed |
| `vdf8-trfe4-cfe1-cdfe1-14` | 14 | +0.448 | +0.319 | +11.2% | +0.807 | CDFE: not passed |
| `vdf7-trfe2-ctfe1` | 10 | +0.478 | +0.449 | +10.9% | +0.900 | CTFE: at the edge |
| `vdf6-trfe3-cfe1` | 10 | +0.492 | +0.361 | +11.4% | +0.454 | CFE: fitted term |
| `vdf7-trfe2-cfe1` | 10 | +0.535 | +0.533 | +9.7% | +0.868 | CFE: fitted term |
| `vdf8-cfe2` | 10 | +0.634 | +0.634 | +9.5% | +0.812 | CFE: fitted term |
| `vdf7-trfe2-cfe1-adj` | 10 | +0.653 | +0.579 | +10.2% | +0.875 | CFE: fitted term |
| `vdf5-trfe5` | 10 | +0.740 | +0.618 | +10.5% | +0.736 | checked |
| `vdf9-cfe1` | 10 | +0.902 | +0.902 | +9.8% | +0.967 | CFE: fitted term |
| added 22:05 UTC, their units checked after the table was made: | | | | | | |
| `vdf9-trfe4-e1-14` | 14 | -0.046 | -0.046 | +10.2% | -0.045 | ethylene: at the edge |
| `vdf9-trfe4-dfe1-14` | 14 | +0.287 | +0.261 | +11.1% | +0.263 | DFE: checked |

- **This is the first of four numbers a row is ranked on, and the only one here.** The margin
  places a row on the strain axis: at 55 V/um a tenth of a cell is polar at +0.15 and half at
  +0.07 (section 1). The loss, the rate and the blocked stress have not been computed under this
  preset, so nothing in this table is a ranking.
- On that one axis: four rows sit in the window (+0.06 to +0.15), and two of the four lean on a
  unit that is not cleanly checked (CTFE at the edge, CDFE not passed). `vdf6-trfe4` (+0.083,
  TrFE alone) and `vdf9-trfe4-tfe1-14` (+0.147, a TFE in the matrix) are the two whose chemistry
  has passed.
- Five rows are below zero, all-trans the lower cell with no field: every row with a reversed
  VDF in it but one, and a CDFE row.
- The TrFE matrix is at +0.22 and the CFE terpolymer at +0.36 (+0.31 at one CFE in twenty):
  under this preset a CFE raises the margin by about 0.14 at one in fourteen. Belt4 had the
  matrix at +0.06 and the terpolymer at +0.32.
- None of these fourteen-monomer cells has a DFT energy of its own. What has been checked is
  each unit in a four-monomer cell.
- The CDFE rows have two or three helix placements in the polished table where the others have
  four, so their lowest helix is the lowest of fewer.
- The two rows added later: an ethylene in the matrix takes 0.27 off the margin (-3.8 kcal/mol
  per unit, a pin for all-trans) and a DFE adds 0.07 (+0.9 per unit). The levels carry what 3.15
  found: up to 0.2 high at this TrFE content. The eight rows that are not here are the seven
  with a nitrile, which no preset has, and the VF row, whose unit failed its check.

Files: `docs/predictions/row_cells_dft_vdf2-trfe1-cdfe1-4.json`,
`chain_split_dft_vdf2-trfe1-cdfe1-4.json`, `chain_conformers_cdfe_pbe_d3bj.json`,
`mlip/row_cells_vdf2-trfe1-cdfe1-4_frames_pbe_d3bj.extxyz` (2 frames),
`margins_chain1-cfe.json` (the 25 rows, and the two added later), `row_cells_dft_vdf2-trfe1-tfe1-4.json`,
`chain_split_dft_vdf2-trfe1-tfe1-4.json`. Tests: `tests/test_row_cells_dft.py`.

### 3.14 Four fourteen-monomer rows in the potential: the levels part from chain1's, the pins do not (2026-11-19)

DFT cannot do a cell of 168 atoms on this machine, and the potential can relax one in ten
minutes. So the four rows the ranking turns on were relaxed in it at their own dilution, each
row's five cells from polyfind's (`examples/mlip_train.py relax`, to 0.02 eV/A; the cells move
by 1.5% in volume at most), and its own energies collected
(`examples/row_cells_frames.py potential`, `docs/predictions/rows14_potential.json`). There is no
DFT on these cells. What the potential's own margin is worth is in 3.13: within 0.11 of DFT on six
small cells, 0.36 off on the seventh.

| all-trans minus the lowest helix, kcal/mol per monomer | `all_s1` | chain1-cfe | belt4 |
|---|---:|---:|---:|
| the TrFE matrix, `vdf10-trfe4-14` | -0.024 | +0.222 | +0.056 |
| the CFE terpolymer, `vdf9-trfe4-cfe1-14` | +0.078 | +0.359 | +0.317 |
| the nitrile alone, `vdf9-trfe4-an1-14` | -0.118 | no preset | +0.197 |
| the leader of the belt4 ranking, `vdf9-trfe3-an1-rvdf1-14` | -0.055 | no preset | +0.037 |
| what a CFE does at one in fourteen | +0.102 | +0.137 | +0.261 |
| what a nitrile does at one in fourteen | -0.094 | | +0.141 |

- **The potential and chain1-cfe agree on the CFE unit**: a pin for the helix of +0.10 and +0.14
  per monomer (+1.4 and +1.9 kcal/mol per unit). Belt4 has twice that.
- **The nitrile at one in fourteen is a pin for all-trans in the potential**, -0.09 per monomer,
  -1.3 kcal/mol per unit: the sign DFT has at one in four and a third of the size (-3.3 to
  -4.3). Belt4 has a pin for the helix of +2.0. So the per-unit number does not carry over
  between the two dilutions, and belt4 is wrong in sign at both.
- **The leader of the belt4 ranking comes out at -0.055**, all-trans the lower by a little: at
  the foot of the window or just under it, 0.09 below belt4's +0.037 and nowhere near the -0.35
  that DFT's per-unit number at one in four would have given. On this evidence it is neither
  confirmed nor thrown out.
- **The levels part.** The potential has the matrix 0.25 below chain1-cfe and the terpolymer
  0.28 below. On the 25 mol% TrFE cell of four monomers the two were 0.07 apart and DFT was
  with chain1 (3.10). One of them drifts between that cell and this one, and nothing here says
  which.
- An observation, not a test: the potential's levels are the ones the materials are known by.
  The 71/29 copolymer is a ferroelectric, polar with no field, and the potential has its cell
  at zero (-0.02); the CFE terpolymer is a relaxor, and the potential has it in the window
  (+0.08). chain1-cfe has both above the window (+0.22, +0.36). These are zero-kelvin cells of
  perfect crystals, so the match may be an accident.
- **What decides it**: DFT on a cell of the matrix's own motif that is still affordable here,
  VVRVVR, six monomers and seventy-two atoms (`vdf4-trfe2-6`), with polyfind and the potential
  on the same cell. Done: 3.15. DFT is with the potential.

**What this reads on.** The margin only, and from a potential that has been 0.36 off once. No
loss, rate or blocked stress.

Files: `docs/predictions/rows14_potential.json`. Test: `tests/test_row_cells_dft.py`.

### 3.15 The level, by DFT on a six-monomer cell of the matrix's motif: DFT is with the potential, and chain1 is 0.18 high (2026-11-19)

The cell is `vdf4-trfe2-6`, VVRVVR: 33 mol% TrFE, one unit in three, two chains of six monomers,
72 atoms. The check is the one every unit had (`examples/unit_dft_check.py`): polyfind's cells
under both presets, the five cells relaxed in the potential, DFT on the all-trans cell and the
helix cell as the potential leaves them (largest force left 0.09 eV/A), and each cell split into
its chains alone and the rest.

The four helix placements of this sequence are one cell: its constitution reads the same both
ways about the CHF carbon, and one repeat of the sequence is half a turn of the helix's pattern.
Polyfind's four cells have one energy to 1e-4 and the potential walks all four in the same 221
steps. So the cell has one margin and no choice of placement in it.

| all-trans minus helix, kcal/mol per monomer | cell | inside a chain | between chains |
|---|---:|---:|---:|
| DFT | **-0.044** | +1.603 | -1.647 |
| `all_s1` alone | -0.096 | | |
| belt4 (shear free -0.014; upright:) | -0.015 | +2.395 | -2.410 |
| chain1 (shear free **+0.139**; upright:) | +0.078 | +2.004 | -1.927 |

- **DFT is with the potential.** The potential is 0.05 low, belt4 0.03 high, chain1 0.18 high.
  The potential and chain1 are 0.235 apart on this cell and 0.246 apart on the
  fourteen-monomer matrix (3.14): the same gap, and here DFT says which side it is on.
- **Belt4 is right again by two errors that cancel**: 0.79 inside the chain and 0.76 between
  chains, as on the 25 mol% cell (0.7 and 0.7).
- **chain1 has not lost what it had; the cancellation it also rests on is incomplete here.**
  Its two parts against DFT's, cell by cell (polyfind's cell upright; the last column is the
  margin with the shear free, the one compared with DFT throughout):

| cell | inside a chain, off by | between chains, off by | margin, off by |
|---|---:|---:|---:|
| TrFE, VVRV | +0.13 | -0.22 | -0.01 |
| CFE, VRVC (with the term of 3.11) | +0.18 | -0.15 | +0.01 |
| CTFE, VRVK (with the term of 3.11) | +0.30 | -0.15 | +0.09 |
| TFE, VRVF | +0.26 | -0.29 | +0.02 |
| TrFE, VVRVVR | +0.40 | -0.28 | +0.18 |

  On every cell chain1 has the all-trans chain too high against the helix chain and the
  all-trans cell too strongly bound against the helix cell. Where the two are equal the margin
  passes; on the six-monomer cell the first has grown to 0.40 and the second has not. The
  passes of 3.10 to 3.13 are of the margin, and the margin is what a row is ranked on, but they
  are not passes of the two parts, and this cell is where that shows.
- **What is new in this cell is the TrFE on the other side of the fold.** In DFT's helix cell the
  two CF2-CHF bonds of a chain that are gauche sit at +65 and -56 degrees: one TrFE on each
  sense. The 25 mol% cell has one TrFE a chain, and its lowest helix, the only one DFT was asked
  for, has it on the favourable sense (in chain1 the other placement is 0.67 per monomer
  higher, 2.7 kcal/mol per unit, a number DFT has not been asked for). Inside the chain chain1
  is 0.5 kcal/mol off per four monomers there and 2.4 per six here. A torsion term is even in
  the angle, so it cannot tell the two senses apart; the difference between them is carried by
  the charges and walls inside the chain. Lone chains with a TrFE on each sense were put to DFT
  to test this directly, and they do not bear the suspicion out: 3.16.

**What it does to the 25 rows of 3.13.** Their levels under `pvdf-dft-chain1-cfe` are not good to
0.1: the preset is right at one TrFE in four and 0.18 high at one in three, and the
fourteen-monomer rows, at 29 mol%, are between the two. On the matrix the three evaluators,
each corrected by its own error on this cell, come to the same place:

| the fourteen-monomer rows, estimated | from `all_s1` (+0.05) | from chain1-cfe (-0.18) | from belt4 (-0.03) |
|---|---:|---:|---:|
| the TrFE matrix, `vdf10-trfe4-14` | +0.03 | +0.04 | +0.03 |
| the CFE terpolymer, `vdf9-trfe4-cfe1-14` | +0.13 | +0.18 | not used: belt4 has the CFE 0.38 high |

- These are estimates: one cell's offsets carried to a cell more than twice its size and to
  29 mol% from 33. No fourteen-monomer cell has a DFT energy. A fourth estimate, the preset
  with a TrFE term (3.16), comes out at +0.08 to +0.20 for the matrix: the agreement of these
  three is not the uncertainty.
- On them the copolymer sits just under the switching window (+0.07 to +0.15 at 55 V/um) and the
  terpolymer at its top or just above; on the fourth both are 0.1 higher. That is the order the materials are known in, and the
  CFE's pin is what moves one to the other: +0.10 to +0.14 per monomer at one in fourteen,
  where the potential and chain1-cfe agree (3.14).
- The differences a unit makes are the part of 3.13's table that stands. Its levels are to be
  read as up to 0.2 high at 29 mol% TrFE and above; 3.16 says why a term for the TrFE alone
  does not mend that.
- Along the chain the all-trans cell is 10.0% longer than the helix cell as the potential relaxed
  them (DFT leaves 0.12 GPa along the chain on the helix cell); chain1 has 11.6% and belt4
  11.5%.

**What this reads on.** The margin, which places a row on the strain axis, and the lattice strain
of one cell. The potential's margins of 3.14 have been put through the loss, rate and blocked
force tables with every other input as belt4 has it: that table and what it rests on are at the
head of `docs/CANDIDATE_MATERIALS.md`.

Files: `docs/predictions/row_cells_dft_vdf4-trfe2-6.json`,
`chain_split_dft_vdf4-trfe2-6.json`, `mlip/row_cells_vdf4-trfe2-6_frames_pbe_d3bj.extxyz` (2
frames), the row's entries in `margins_{chain1,belt4}_dft_rows.json` and
`chain_split_polyfind_dft_rows.json`. Test: `tests/test_row_cells_dft.py`.

### 3.16 Lone VDF-TrFE chains on each sense, a TrFE torsion term that is not adopted, and where the evaluators stand (2026-11-19)

**The lone chains.** The chain of alternating VDF and TrFE, alone in a box, all-trans and its
helix on each gauche sense of the CF2-CHF bond (`examples/chain_conformers.py`, chain `trfe2`;
each settled in DFT to 0.004 eV/A; `docs/predictions/chain_conformers_trfe_pbe_d3bj.json`):

| helix against all-trans, kcal/mol per monomer | DFT | chain1 | belt4 |
|---|---:|---:|---:|
| the CF2-CHF bond at -62 degrees (the favourable sense) | -2.802 | -3.071 | -3.741 |
| at +49 degrees (the other sense) | +0.631 | +0.551 | +0.439 |
| the difference between the senses | 3.43 | 3.62 | 4.18 |

- The two senses are 3.4 kcal/mol per monomer apart, 6.9 per TrFE unit: a TrFE all but forbids
  the helix on its wrong side.
- **chain1 has that difference to 0.19, and the suspicion of 3.15 is wrong.** Its error is on
  the favourable sense, the one the 25 mol% cell tested: the helix 0.27 too low there and 0.08
  too low on the other. Both helices are too low, so an even term can raise them.

**A torsion term for the CF2-CHF bond, fitted to those two rows** on `pvdf-dft-chain1-cfe`
(`examples/unit_torsion_fit.py --type CF2-CHF`; from (+0.450, -0.522, +0.746) to (+0.079,
+0.400, +0.326); the lone chains then at -2.850 and +0.639), and tried on the five cells with a
TrFE that DFT has, none of them in the fit:

| all-trans minus helix against DFT | DFT | `all_s1` | belt4 | chain1-cfe | chain1-cfe with the TrFE term |
|---|---:|---:|---:|---:|---:|
| TrFE, VVRV | +0.678 | -0.07 | -0.04 | -0.01 | -0.10 |
| TrFE, VVRVVR | -0.044 | -0.05 | +0.03 | +0.18 | +0.02 |
| CFE, VRVC | +0.682 | +0.09 | +0.38 | +0.01 | -0.09 |
| CTFE, VRVK | +0.315 | -0.11 | +0.42 | +0.09 | -0.05 |
| TFE, VRVF | -0.076 | -0.36 | +0.30 | +0.02 | -0.13 |
| the largest | | 0.36 | 0.42 | 0.18 | 0.13 |
| root mean square | | 0.18 | 0.29 | 0.09 | 0.09 |
| ethylene, VRVE (added after the five; not run with the term) | -0.146 | +0.06 | +0.65 | -0.09 | |
| VF, VRVY (added after the five; DFT's forces left 0.35 eV/A) | +0.290 | -0.02 | +0.16 | +0.19 | |
| DFE, VRVX (added after the five; DFT's forces left 0.29 eV/A) | +0.961 | **-0.65** | +0.27 | -0.01 | |

With the three cells added after the five, and the nitrile and CDFE cells of 3.12 and 3.13, every
unit of the 35 rows has now been put to DFT in a small cell. The potential misses the two units
it has fewest frames of, TFE by 0.36 and DFE by 0.65, and is within 0.11 on the other ten cells.
chain1-cfe is within 0.09 on seven cells, 0.18 to 0.19 off on two (the six-monomer TrFE cell,
VF), right by cancellation on CDFE and wrong on the two nitrile cells. Each is to be used where
it has passed, and the labelled cells of today are frames for the potential's next fit.

- **The term puts the six-monomer cell right and every four-monomer cell 0.05 to 0.13 low.**
  With the chain right, the error between chains (the all-trans cell 0.15 to 0.29 too strongly
  bound, 3.15) is no longer cancelled. The root mean square is the same with it and without.
  **It is not adopted**: no preset is registered, and the fit is kept as a record
  (`unit_torsion_trfe_chain1.json`, `margins_trfe_torsion_dft_rows.json`).
- What would have to change is between chains, in the walls and charges every unit shares, and
  that is a refit of the preset and not a term for one unit.
- **No evaluator is within 0.1 on all five.** The potential is within 0.11 on four and 0.36 off
  on the TFE cell. chain1-cfe is within 0.09 on four and 0.18 off on the six-monomer cell.
  Belt4 is within 0.04 on the two cells without a third unit and 0.3 to 0.4 off on the three
  with one. On the cells of 3.12 and 3.13 the potential is within 0.08 (two nitrile cells,
  CDFE) and polyfind is wrong in sign on the nitrile.
- So the level of the fourteen-monomer matrix stays an estimate from evaluators that disagree by
  more than the 0.1 that is needed: -0.024 in the potential, +0.222 in chain1-cfe, and +0.100
  in chain1-cfe with the TrFE term (the terpolymer +0.239 with it: a CFE still worth +0.14).
  The corrections of 3.15 bring the first two together at +0.03 to +0.04, but each is one
  cell's offset; the third, corrected by its own errors on the two TrFE cells (-0.10 at one in
  four, +0.02 at one in three), is +0.08 to +0.20. The matrix is somewhere from +0.03 to +0.2,
  and the terpolymer 0.10 to 0.14 above it. Only DFT on those cells narrows that.

**DFT at fourteen monomers does not fit this machine.** The leader's helix cell (170 atoms, 350
bands, 6 k-points) asks for 15.7 GB at the recipe's cutoffs and the machine has 15; the run was
stopped before it ran out. Two such cells for each of the leader, the matrix and the terpolymer
are a job for a rented machine with more memory.

Files: `docs/predictions/chain_conformers_trfe_pbe_d3bj.json`, `unit_torsion_trfe_chain1.json`,
`margins_trfe_torsion_dft_rows.json`, `margins_trfe_torsion_rows14.json` (the matrix and the
terpolymer under the term). Test: `tests/test_row_cells_dft.py`.

### 3.17 A kink's energy against DFT, for the first time (2026-11-19)

The frequency and the loss of every cell in `docs/CANDIDATE_MATERIALS.md` come from the energies
of kink pairs (T G+ T G- T), computed under belt4, and none had been put to DFT. This is the
first: the six-monomer cell of 3.15 (`vdf4-trfe2-6`, 72 atoms) with a kink pair on both chains,
in register, the lattice held at the all-trans cell's (`examples/dilute_kink.py`, stage
`flexible_host`; `examples/row_cells_frames.py kink` writes the cells as atoms; the potential
relaxes the atoms with the cell held; DFT is asked once for each relaxed cell). A kink pair
does not fit a six-monomer repeat: its gauche bonds flatten to about 150 and -140 degrees, and it
goes back to all-trans when the whole cell is freed. So this is a squeezed kink at a held
lattice, the same structure in every evaluator, and not the kink sheet of the fourteen-monomer
cells.

| kcal/mol per kink pair | volume per monomer, A3 | DFT | `all_s1` | polyfind |
|---|---:|---:|---:|---:|
| at belt4's all-trans lattice | 53.5 | **6.01** | 6.68 | belt4 6.60 |
| at the lattice the potential relaxes the all-trans cell to (DFT's stress on it under 0.12 GPa) | 57.5 | **4.15** | 4.75 | |
| at chain1's all-trans lattice | 57.8 | not asked | 5.01 | chain1 3.66 |

- **The lattice decides a third of it.** In DFT the same kink costs 6.0 in belt4's cell and
  4.15 in the cell DFT is at rest in. Belt4's all-trans cell is 7% too dense (0.31 kcal/mol
  per monomer above the potential's own), and a kink needs room.
- **Belt4's kink is 10% dear at its own lattice and 59% dear against DFT at DFT's.** The first
  is the evaluator, the second is what its tables use: 6.60 where DFT has 4.15.
- **chain1, whose cell is right to half a per cent, is 12% cheap**: 3.66 against 4.15.
- **The potential is 0.60 to 0.67 high at both lattices** (0.11 kcal/mol per monomer on the
  kinked cell, 0.00 on the all-trans one): it was never shown a kink. It has the lattice's part
  right, 1.92 against DFT's 1.85. The three labelled frames are in the repository for the next
  fit.
- **What this does not test.** The rates in the tables are ratios between cells: the kink
  sheet's energy in one cell against another's. One cell cannot say whether belt4's 59%
  applies to every cell alike (then every ratio is too wide by the same factor in the
  exponent) or differently. The four fourteen-monomer sheets are being relaxed in the potential
  for that.

- **Added 2026-11-20.** These are single points on the potential's geometry. A step toward
  DFT's own minimum (the potential tilted to DFT's forces, lattice held) expects DFT to move by
  0.04 to 0.05, so the numbers stand. What they are numbers of is an upright kink at a held
  lattice: 3.19 finds that a kink sheet free to shape its cell is a different structure with
  a different energy.

Files: `docs/predictions/kink_dft_vdf4-trfe2-6.json`,
`mlip/kink_vdf4-trfe2-6_frames_pbe_d3bj.extxyz` (4 frames). Test:
`tests/test_row_cells_dft.py`.

### 3.18 The kink sheets of four fourteen-monomer cells in the potential (2026-11-19)

**Withdrawn as a ranking, 2026-11-20 (3.19).** These eleven relaxations stopped at 0.02 eV/A,
where a cell does not answer a stress under 0.35 GPa. Ten of them stopped within 1.2 A of an
upright cell, which is not a minimum; the eleventh, the leader's 0.65, is the one whose cell
had sheared 3 to 4 A along the chains on the way to the inclined sheet of 3.19. So the 0.65
is not a property of the leader, the order of "lowest sheets" compares cells stopped at
different points of one slope, and the frequencies made from it are not results. The numbers
are left below as what was found.

The kink sheet is the structure the rates are made of: a kink pair on every chain, in register,
in a cell of its own, at one site of the sequence. Eleven sheets, two or three sites of each of
the four cells of 3.14 (the sites belt4 has cheapest, `examples/kink_site_scan.py`), were written
as atoms (`examples/row_cells_frames.py kink --scan`), relaxed in the potential with the cell
free, and set against the row's all-trans cell there (`sheets-potential`;
`docs/predictions/sheets14_potential.json`, `sheets14_potential_more_sites.json`,
`sheets14_potential_third_batch.json`). Ten kinks are still kinks after the relaxation (gauche
bonds near 80 and -70 degrees); one goes back to all-trans. There is no DFT on
these cells; on the one kink DFT has seen the potential is 0.6 high (3.17), which is taken to be
common. Each relaxation starts from belt4's cell for that site and stops at the nearest minimum,
so a value is an upper bound on that site's sheet.

| kink sheet, kcal/mol per kink pair | site (hand) | `all_s1` | belt4 | `all_s1` minus belt4 |
|---|---|---:|---:|---:|
| the 71/29 matrix | 12 (mirrored) | 4.54 | 4.66 | -0.12 |
| | 12 | **3.59** | 4.81 | -1.22 |
| the cell of the film's composition, 64/29/7 (CFE) | 10 (mirrored) | 4.59 | 3.42 | +1.17 |
| | 0 | **4.16** | 3.95 | +0.21 |
| | 8 (mirrored) | 6.60 | 3.57 | +3.03 |
| the nitrile alone | 8 (mirrored) | **2.95** | 2.15 | +0.80 |
| | 0 | 3.40 | 3.25 | +0.15 |
| | 10 (mirrored) | the kink does not hold: back to all-trans | 2.26 | |
| the leader (nitrile with reversed VDF, 21 mol% TrFE) | 8 (mirrored) | 2.49 | 1.39 | +1.09 |
| | 18 | 2.41 | 1.86 | +0.56 |
| | 18 (mirrored) | **0.65** | 1.94 | -1.28 |

- **Site by site the two evaluators do not agree**: the difference runs from -1.3 to +3.0, and
  one site belt4 has at 2.26 holds no kink in the potential. Belt4's profile of a sheet along
  the sequence, which its sums over sites are made of, is not confirmed. A first reading of
  this section, on the first four cells alone, had the leader 2.10 below the film's cell in
  the potential against belt4's 2.02 and called its speed confirmed to 0.1; that agreement was
  of one site each and does not survive the other seven.
- **What does survive is the order of the lowest sheets found**: the leader 0.65, the nitrile
  alone 2.95, the matrix 3.59, the film's cell 4.16.
- **The leader has one sheet at 0.65**, a kink pair for about kT, where its two other sites are
  at 2.4 to 2.5. If it is real the cell is soft in every way at once (its margin is at zero
  too, 3.14), and a population of sheets that is no longer dilute is outside what the rate
  model's exponential describes. It is one relaxation from one starting cell and has no DFT.
- **The CFE's cut of the kink sheet is not there.** Belt4 has a CFE making the matrix's lowest
  sheet 1.24 cheaper; in the potential the film's cell is 0.6 dearer than the matrix on the
  sites tried, three of one and two of the other. Belt4's chlorine is the unit it has 0.38
  wrong on the margin (3.10).
- With each cell's sites moved so that belt4's lowest sits at the lowest found here
  (`four_axes_at_margins.py --sheet-shift`), the rate model that sums over sites gives the
  matrix 697 kHz where belt4's sheets gave 33 and the nitrile alone 1.6 MHz for 1.8, against
  the film's 364 kHz. For the leader the formula gives 83 MHz from the 0.65 site, which is
  outside its range, and 4.4 MHz if its floor is the next site, 2.41; belt4's sheets gave 6.9
  (`docs/predictions/four_axes_potential_margins_sheets.json`; the table is at the head of
  `docs/CANDIDATE_MATERIALS.md`). All of this keeps belt4's profile above the floor, which the
  first bullet says is not to be trusted: the rates need every site of each cell in the
  potential, 112 relaxations of a quarter of an hour here, and DFT on the lowest of them.

Files: `docs/predictions/sheets14_potential.json`, `sheets14_potential_more_sites.json`,
`sheets14_potential_third_batch.json`, `four_axes_potential_margins_sheets.json`. Tests: `tests/test_four_axes_at_margins.py`.

### 3.19 A kink sheet in its own cell against DFT: the tables' upright sheet is not a minimum, and the inclined one is below all-trans (2026-11-20)

The first DFT on a kink sheet as the rate model means it: PVDF, ten monomers, two chains each
with one kink pair (120 atoms), the kinked crystal of belt4's packing search
(`kink_migration_belt4.json`, `vdf10@0`; `examples/row_cells_frames.py kink --sheet-of`) and
the all-trans cell, each relaxed in the potential with its cell free and given to DFT.

**The upright sheet.** Belt4's sheet is upright: the kinks of all the chains are at one height
and the cell has no shear along the chains.

| kcal/mol per kink pair | DFT | `all_s1` | DFT leaves on the kinked cell |
|---|---:|---:|---|
| the cell the potential stops at from belt4's (stop 0.02 eV/A, as every cell of 3.14 and 3.18) | 4.82 | 4.68 | 0.16 eV/A, 0.19 GPa |
| one step on, in the potential tilted to DFT's forces and stress there (0.57 A, 0.8 degrees of shear) | 4.14 (4.03 expected) | 4.43 | 0.06 eV/A, 0.18 GPa |

Belt4 has it at 5.30 in its own cell (5.06 with the pendants tilted). DFT is not at rest on
either cell, and the stress along the chains is as large after the step as before it: the
energy is on a slope, about a quarter of a kcal/mol per kink pair for each degree of shear.
The potential is on the same slope (0.29 GPa of its own stress on the cell it stopped at), and
its 4.68 is not its minimum: it is 0.25 lower on the stepped cell. The stop is what hid this.
The cell's share of a 0.02 eV/A stop is the stress times the volume over the number of atoms,
so no stress under 0.35 GPa moves a cell, whatever its size. That is worth a few hundredths of
a kcal/mol per monomer, which a margin survives; a kink pair's energy multiplies it by the ten
to fourteen monomers of the chain.

**The inclined sheet.** Relaxed to 0.003 eV/A the potential's cell goes on shearing and stops
(900 steps) at a = a0 - c0 and b = b0 - c0 of the all-trans lattice, c0 one monomer along the
chain: each chain's kink is one monomer along from the kink of its neighbour, in a and in b.
The kink is sharper (60, 172, -58 degrees, where the upright one has 81, 178, -68), the rest
of the chain is exactly trans, and the cell is denser (53.3 A3 a monomer; the upright sheet
55.1, all-trans 52.3). Polyfind's cell has this freedom (`MixedModel(shear=True)`, which
`examples/shear_polish.py` uses for the margins' cells) and its kink tools do not use it;
`examples/kink_sheet_shear.py` relaxes a stored sheet with it, from the upright cell and from
starts a monomer along a and b.

| kcal/mol per kink pair, against all-trans | DFT | `all_s1` | belt4 | chain1 |
|---|---:|---:|---:|---:|
| the upright sheet | 4.82 falling (above) | 4.68 falling | 5.30 (5.28 after a kick: a local minimum) | 5.29 |
| the inclined sheet | **-1.96** | -1.66 | -1.17 | -0.78 |

- **All four have the inclined sheet below the all-trans crystal**, 6.8 kcal/mol per kink pair
  under the upright one in DFT. DFT's number is a single point on the potential's cell, and
  here the potential's minimum is DFT's: 0.07 eV/A and 0.04 GPa are left, where the upright
  cell had 0.16 and 0.19. The all-trans cell is at its own minimum upright: slid by the
  same monomer it is the same crystal at the same energy, and no other slide is lower (eight
  starts in the primitive cell, in the potential; polyfind's shear relaxation agrees).
- **It is what PVDF's margin says it should be.** PVDF's helix is 0.73 kcal/mol per monomer
  below all-trans (the anchor), and a kink pair T G+ T G- is one turn of that helix, two
  monomers. Two in ten at the helix's energy is -1.45; the four numbers are -0.8 to -2.0. An
  inclined kink sheet is a slab of the other phase one turn thick, and its faces cost little.
- **In belt4 the upright sheet is a local minimum** (a kick of 0.05 A returns to 5.28), which
  is why its scans never left it. In DFT and in the potential it is not one.
- **Other staggers are dearer in belt4**: of nine starts a monomer or two along a and b,
  two end in this inclined sheet, one returns to the upright one and six end in other sheets
  at +0.4 to +4.4.

**What this does to the tables.** Every kink cell of `examples/dilute_kink.py`,
`kink_migration.py` and `kink_site_scan.py` is upright, so every kink-sheet energy in
`docs/CANDIDATE_MATERIALS.md`, and the frequency and the loss made from them, is the energy of
a structure 6 to 7 kcal/mol above the sheet's own minimum in PVDF. Whether the differences
between candidates survive is not known: an inclined sheet's energy follows the margin (about
minus two margins plus what its faces cost), and the upright one's did not. The first nine
sites of the four cells of 3.14 with the shear free, in belt4, already reorder them: the
leader's three sheets do not incline (1.36 to 1.99, as upright), the matrix's lowest goes from
4.66 to 1.14, and the film's-composition cell from 3.95 to -0.13 at one site (the scan of all
112 is running). The strain and the
blocked force do not use the sheets. The rate model itself assumes sheets that are rare, with
an energy of several kT; a sheet near or below zero is not an excitation of the polar crystal
but the start of the other phase, and the model has to be looked at again with the right
energies in hand.

**What 3.18 was.** Its eleven sheets were stopped at 0.02 eV/A: ten near upright, and the
leader's 0.65 part way down this slope (the only cell whose third vector had tilted, 3 to 4 A).
Its ranking is withdrawn there.

Files: `docs/predictions/kink_dft_vdf10.json`, `mlip/kink_vdf10_frames_pbe_d3bj.extxyz`
(4 frames: all-trans, the upright sheet, the step, the inclined sheet),
`kink_sheet_shear_belt4.json`. Tests: `tests/test_row_cells_dft.py`,
`tests/test_four_axes_at_margins.py`.

### 3.20 Every kink site of the four cells with the shear free: the film cell's sheets fall below all-trans and the leader's stay upright (belt4, 2026-11-20)

3.19 found on PVDF that the upright kink sheet of the tables is not the sheet's minimum. This is
the same freedom given to every kink site of the four fourteen-monomer cells of 3.14, in
polyfind under belt4 (`examples/kink_sheet_shear.py --key "ROW@0" --sites all`, and `"ROW@0~"`
for the other hand of the kink pair). A site is one position of the kink pair along the
fourteen monomers; the sheet has that kink on both chains of the cell. Each of the 108 sites at
which the upright scan kept its kink was relaxed from the upright cell (kicked off it) and from
four cells sheared by one monomer along a and along b (both signs of each), and the all-trans
cell was relaxed the same way. Energies are kcal/mol per kink pair against all-trans, pendants
as built. About six minutes a site on one core, eleven core-hours in all.

| cell | sites | cheapest sheet upright | cheapest with the shear free (units, shear along a, b in monomers) | sites below all-trans | sites that end sheared | mean fall |
|---|---:|---:|---|---:|---:|---:|
| the film's composition, VVRVVRVCVVRVVR | 28 | +3.42 | **-2.27** (V-CFE; 0.0, -1.8) | 9 | 28 | 4.3 |
| the nitrile alone, VVRVVRVAVVRVVR | 26 | +2.15 | **-2.08** (TrFE-V; 0.0, -1.8) | 6 | 24 | 3.2 |
| the leader, IVRVVRVAVVVVVR | 26 | +1.40 | **+0.65** (TrFE-V; 0.0, -1.7) | 0 | 6 | 0.8 |
| the matrix, VVRVVRVVVVRVVR | 28 | +4.66 | **+0.94** (V-TrFE; -0.9, -1.0) | 0 | 24 | 3.2 |

- **The order of the four reverses.** Upright, the leader had the cheapest sheet and the film
  cell the second dearest. With the shear free the film cell and the nitrile cell have sheets
  2 kcal/mol below their all-trans crystals, and the leader and the matrix have none below.
- **The sheets take three cells**: sheared by nearly two monomers along b (most sites of one
  hand in the film cell, the nitrile cell and four of the leader's), by one monomer along each
  of a and b (the other hand, and the matrix), or left upright (20 of the leader's 26).
- **The fall follows belt4's margin**, as 3.19 said it should. belt4's margins of the four are
  +0.32, +0.20, +0.04 and +0.06 kcal/mol per monomer (film cell, nitrile cell, leader, matrix:
  the helix lower in all four), and the cheapest sheets are -2.27, -2.08, +0.65 and +0.94. A
  sheet is a slab of the other phase; where the other phase is well below, so is the sheet.
  The helix cell is still the lowest of the three in every case (4.4, 2.8, 0.5 and 0.8 below
  all-trans per chain of fourteen), so no margin is displaced by a sheet.
- **The all-trans cell does not gain from the same shears**, so the fall is the kink's. From
  twelve starts sheared by whole monomers the upright all-trans cell is the lowest in all three
  cells tried: the others end 1.7 to 5.3 kcal/mol per cell higher in the film cell, 2.3 to 8.2
  in the matrix and 4.4 to 11.7 in the leader (`examples/shear_polish.py --starts`).
- **The margins do not move.** Every margin cell of the film cell and the leader (all-trans and
  four helix conformers) relaxed from the same twelve whole-monomer shears ends where the
  upright start had taken it: rigid margins +0.251 and +0.163, as
  `sheared_margins_belt4.json` has them. One helix conformer of the leader's that does not set
  its margin found a cell 1 kcal/mol lower.
- **Pendant tilts change little.** Relaxed on the sheet and on the all-trans cell
  (`polyfind.mixedcell.tilt_relaxation`), the film cell's -2.27 becomes -2.36, the leader's
  +0.65 becomes +0.50 and +1.36 becomes +1.05. The leader's all-trans cell sheared by two
  monomers is 7.4 above the upright one rigid and 7.6 with tilts.

**Why the leader's sheets stay upright.** A cell sheared by whole monomers has the same
backbone packing and another arrangement of the units between chains: upright, each unit sits
beside the same unit of the next chain; sheared, beside a different one. For the shear by two
monomers along b, which the film cell's and the nitrile cell's cheapest sheets take, the
all-trans cell pays 1.7 kcal/mol per cell in the film cell, 2.3 in the matrix, 2.5 in the
nitrile cell and **7.4 in the leader**. The leader's sheets do exist at that shear (nine of
one hand's fourteen sites reach it), and at the sites the two cells share they are 2.5 to 3.9
kcal/mol per kink pair above the nitrile cell's: the difference between what the two all-trans
cells pay, which is 2.5 per chain. So the leader's kinks gain what the others' gain and its
cell pays more for the arrangement. The leader is the cell with the reversed unit, whose CF2
and CH2 face the wrong neighbours unless the next chain's reversed unit is beside it.

**This is where belt4 is weakest, by 3.9.** belt4 charges about three times DFT's energy for
sliding plain PVDF chains along each other, and has the reversed-unit candidate's sheared
all-trans cell 0.43 kcal/mol per monomer above the upright one where DFT has it 0.35 below.
The cost that holds the leader's sheets upright is of that kind. Under chain1-cfe, whose wall
was fitted with those energies, the film cell's all-trans cell sheared by two monomers along b
is 0.9 kcal/mol per cell above the upright one (belt4 1.7) and the leader's 6.1 (belt4 8.7),
and most other starts slide back to the upright cell. So the preference is weaker there and
not gone. Neither preset has a DFT test on a fourteen-monomer cell.

**A periodic cell puts every unit beside its like.** In a statistical copolymer a reversed unit
has another reversed unit beside it about one time in fourteen. The upright cell is the
arrangement such a material cannot have throughout, and for the leader it is the arrangement
the all-trans cell needs. What the leader's all-trans phase costs when its reversed units are
not in register is a question about the model cell that the margins of every row with a
reversed unit carry, and it is open.

**On the four axes** (`examples/four_axes_at_margins.py --sites
docs/predictions/kink_site_scan_belt4_shear.json`; margins the potential's of 3.14, everything
else belt4's; `docs/predictions/four_axes_potential_margins_shear.json`). The rate model takes
the frequency as the film's measured 364 kHz times the ratio of the kink sites' populations to
the film cell's. With sheets below all-trans the sum of exp(-E/kT) is no longer a population
(a site cannot hold more than one sheet), so the table gives it two ways: the sum as before,
and each site's share of time with a sheet, 1/(1 + exp(E/kT)), which is the same number where
sheets are rare. 34 kHz drive, AC field at the heat budget:

| cell | margin (the potential's) | frequency, by the sum; by the share | tan delta | AC strain | blocked force per area through the crystal, pN/nm2 |
|---|---:|---|---|---|---|
| the film's composition (the anchor) | +0.078 | 364 kHz | 0.088 | 0.081% | 7.4 |
| the nitrile alone | -0.118 (taken as 0) | 306 kHz; 335 kHz | 0.066; 0.060 | 0.066%; 0.068% | 4.5; 4.7 |
| the leader | -0.055 (taken as 0) | 4 kHz; 30 kHz | 0.85; 0.65 | 0.007%; 0.020% | 0.4; 1.1 |
| the matrix | -0.024 (taken as 0) | 5 kHz; 34 kHz | 1.12; 0.53 | 0.008%; 0.021% | 0.9; 2.4 |
| the measured film's channel as calibrated | | 364 kHz | 0.060 | 0.066% | 6.5 |

On these sheets the leader is with the matrix: its relaxation is at or below the drive, so it
loses most of the field's work and follows little of it. The nitrile alone is level with the
film. None of the four is better than the film on the four axes together.

**The two cheapest sheets under two more evaluators (the same day).** belt4's cheapest sheet of
the film cell (site 12, at the chlorine) and of the leader (site 10), relaxed again elsewhere,
kcal/mol per kink pair against the same evaluator's all-trans cell:

| | belt4 | chain1-cfe | the potential, 0.004 eV/A |
|---|---:|---:|---:|
| the film cell's sheet at V-CFE | -2.27 | -1.65 | **-3.60** |
| the leader's sheet at TrFE-V | +0.65 | +1.66 | **+0.05** |
| the leader's above the film cell's | 2.92 | 3.31 | 3.65 |
| the nitrile cell's sheet at TrFE-V | -2.08 | not run | **+0.36** |
| the matrix's sheet at V-TrFE | +0.94 | not run | **+3.13** |

- **chain1-cfe** (`examples/kink_sheet_shear.py --preset pvdf-dft-chain1-cfe --starts ""
  --from-sheets ...`). Its own sheared starts slide back to the upright cell, because it charges
  a third of belt4's price for sliding, so each site is started from the sheets belt4 found.
  Nine sites of the film cell: eight are below all-trans (-0.36 to -1.65, and -4.50 at CFE-V of
  the other hand, where belt4 has -1.17). The one site of the leader's is the +1.66; its
  upright sheet loses the kink and two of its starts collapse (the preset has no nitrile term),
  and the leader was stopped there. `docs/predictions/kink_sheet_shear_chain1cfe_sites.json`.
- **The potential** (`examples/row_cells_frames.py kink --scan ... --sheared ...` writes the
  sheet, `mlip_train.py relax --fmax 0.004` relaxes it, `row_cells_frames.py sheets-potential
  --pairs` collects). 168 and 170 atoms, 560 and 370 steps from belt4's sheared cells; the
  sheets keep their kinks. The leader's all-trans cell was relaxed again to the same stop: it
  falls 1.15 kcal/mol per cell from where 3.14's stop of 0.02 eV/A left it, which is more than
  the leader's sheet (against the old reference the sheet would have read -0.53).
  `docs/predictions/sheets14_potential_sheared.json`, the cells in
  `mlip/sheets14_sheared_in_all_s1.extxyz`.

So the gap does not belong to belt4: three evaluators with different errors put the leader's
sheet 2.9 to 3.7 kcal/mol above the film cell's at these two sites. The sites are belt4's
choice; the potential's own cheapest site of either cell is not known (one site each took two
hours of this machine).

**The nitrile alone does not stay with the film.** belt4 had its cheapest sheet level with the
film cell's (-2.08); in the potential the same sheet, kink kept, is 0.36 above all-trans. That
is belt4's nitrile error again (3.12: it makes the nitrile a pin for the helix, and a sheet
follows the margin). So the row of the four-axis table above that puts the nitrile cell at 306
to 335 kHz is belt4's and is not supported: in the potential the film cell is the only one of
the four with a sheet below all-trans.

**Nor is the leader level with the matrix there.** The matrix's cheapest sheet is 3.13 in the
potential (belt4 0.94), so the potential's order at these four sites is the film cell (-3.60),
the leader (+0.05), the nitrile alone (+0.36), the matrix (+3.13): the leader three and a half
above the film cell and three below the matrix, where belt4 had the nitrile cell with the film
and the leader with the matrix. What all three evaluators share is the film cell first and the
leader about three behind it.

**A loose stop costs the potential's fourteen-monomer cells as much as their margins.** The
leader's all-trans cell moved 0.04 kcal/mol per monomer between 0.02 and 0.004 eV/A; the same
cell relaxed to 0.03 from belt4's denser cell stayed 0.18 per monomer above the tight one.
3.14's margins for these four cells (-0.12 to +0.08) came from helix cells stopped at 0.02 and
all-trans cells stopped at 0.02 or 0.005. Each cell's all-trans cell and lowest helix cell,
relaxed again to 0.005 or tighter (`docs/predictions/rows14_potential_tight.json`; each helix
is still a helix, the leader's with twelve of its fourteen gauche bonds):

| cell | margin in 3.14 | margin at the tight stop |
|---|---:|---:|
| the matrix | -0.024 | **+0.016** |
| the film's composition | +0.078 | **+0.129** |
| the nitrile alone | -0.118 | **-0.109** |
| the leader | -0.055 | **-0.070** |

The matrix changes sign and the film cell gains 0.05. Every cell that had been stopped at 0.02
fell 0.02 to 0.05 kcal/mol per monomer (the four helix cells and the two nitrile rows'
all-trans cells); the two all-trans cells already at 0.005 did not move. The film cell's row of the
four-axis table at +0.129 reads 364 kHz, tan delta 0.061, 0.071% and 6.2 pN/nm2
(`four_axes_potential_tight_margins_shear.json`), on the measured film's calibrated channel
(0.060, 0.066%, 6.5). The leader's and the nitrile cell's margins stay below zero and their
rows do not change; the matrix's row barely moves (5 to 34 kHz). The other three helix
placements of each cell have not been redone.

**In the potential the film cell's lowest state is the kinked one.** Per chain of fourteen
monomers against all-trans, all at the tight stop:

| cell | lowest helix cell | cheapest sheet tried (one kink pair a chain) |
|---|---:|---:|
| the film's composition | -1.81 | **-3.60** |
| the matrix | -0.23 | +3.13 |
| the leader | +0.98 | +0.05 |
| the nitrile alone | +1.53 | +0.36 |

One kink pair at the chlorine is worth twice the whole helix. Both presets have it the other
way (belt4: helix -4.4, sheet -2.3; chain1-cfe: helix -5.0, sheets to -4.5). If the potential
is right, the state the field has to convert in the film's composition is not the helix but
all-trans chains with a kink held at each chlorine, and the margin as this project defines it
(all-trans minus the lowest helix) does not name the film cell's lowest state. The potential
has been 0.3 kcal/mol off DFT on PVDF's sheet (3.19) and 0.05 per monomer off on a helix level
(3.15), which is not enough to reverse 1.8 per chain; it has never been tested on a chlorine
unit beside a kink. This is the first thing DFT should be asked at this size.

**What this is and is not.**

- It is belt4's energy surface searched properly for the first time, with two sites carried to
  two more evaluators. It is not yet a result on the materials: no DFT has seen a
  fourteen-monomer cell, and belt4's errors are large at these sites (the paragraph above on
  sliding, the chlorine's torsion of 3.11 at the film cell's cheapest site, and the nitrile's
  margin of 3.12).
- Four sheared starts a site. The film cell's sheets were found by them; a sheet of the
  leader's at a shear they do not reach is possible, though its all-trans cell's cost at every
  shear tried says there is little room.
- The rate model is one model. Its premise, that a kink sheet is a rare excitation whose
  population limits the rate, does not hold for a cell with sheets below all-trans; the share
  reading bounds the error and does not replace the premise.
- The leader's standing on frequency and loss in the tables before this (83 MHz, then 3.8 MHz)
  came from upright sheets and is withdrawn (3.19). This table is not its replacement until
  an evaluator checked against DFT has relaxed the same sheets.

Files: `docs/predictions/kink_sheet_shear_belt4_sites.json` (every site's starts, with the
relaxed states), `kink_site_scan_belt4_shear.json` (the scan with each site at its lowest
sheet), `four_axes_potential_margins_shear.json`, `kink_sheet_shear_chain1cfe_sites.json`,
`sheets14_potential_sheared.json`, `rows14_potential_tight.json`,
`four_axes_potential_tight_margins_shear.json`. Tests: `tests/test_four_axes_at_margins.py`.
