# A machine-learned potential as a second opinion on the margins

State as of 2026-11-16. The PI's direction, in order: "Go on 1" (a machine-learned interatomic
potential as a second opinion on polyfind's margins, tested before it is trusted), then, on
seeing pretrained ones tested: "You're using ots solutions? Thought we would train our own?"

Section 1 is the test. Section 2 is done: no pretrained potential passes it. Section 3 is our
own, and says for each part whether it is done, running or not started. Where it stands on a
real candidate is 3.8: not settled, and waiting on the fit of 3.7. In short: two short
fits on 76 PVDF frames were right to 0.1 to 0.2 kcal/mol per monomer where the pretrained ones
are off by 1.2 to 3.9 (3.3). Four long fits on a rented card, on 130 frames, are right to 0.01
on the packings that have frames near them and off by 0.07 to 0.21 on a packing they were not
shown; DFT single points on the cells they relax are right to 0.02 either way (3.6). So the
potential is used to relax a cell and DFT is asked once for its energy; the potential's own
energy is trusted only where it has been checked.

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

This reads on the margin alone, which places a cell on the strain axis; it says nothing yet of
loss, rate or blocked force.
