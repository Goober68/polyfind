# A machine-learned potential as a second opinion on the margins

State as of 2026-11-16. The PI's direction, in order: "Go on 1" (a machine-learned interatomic
potential as a second opinion on polyfind's margins, tested before it is trusted), then, on
seeing pretrained ones tested: "You're using ots solutions? Thought we would train our own?"

Section 1 is the test. Section 2 is done: no pretrained potential passes it. Section 3 is our
own, and says for each part whether it is done, running or not started. In one line: our own
potential, after two short fits on 76 PVDF frames, is right to 0.1 to 0.2 kcal/mol per monomer
where the pretrained ones are off by 1.2 to 3.9, and that is not yet inside the 0.1 needed
(3.3). What separates it from there is training time, which a rented GPU buys (3.5).

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

### 3.3 Training and its tests (two fits done; neither is inside the 0.1 with room to spare)

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

### 3.4 The other units (168 frames built; the 36 relaxed cells labelled, the other 132 being labelled)

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

It was wrong both ways. On the 36 relaxed cells (labelled; 3.2 hours of DFT) the DFT stress
across the chains is -1.3 to -2.6 GPa on the all-trans cells the sampler left at PVDF's density
(TrFE and CFE: compressed, they want 7 to 12% more along each lateral axis) and +0.1 to +1.4
GPa on the all-trans cells it expanded (a little too large, mostly along b). Along the chain
every all-trans cell is under -0.2 to -2.2 GPa: DFT wants the chains straighter than the
sampler left them. The helix cells are within 1.6 GPa either way.

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

### 3.5 What separates it from the 0.1 (recommended to the PI; not started)

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
