# The lattice-curvature column: Gamma-point phonons of every screened chemistry's polar and antipolar cell

`examples/screen_phonons.py`, 2026-09-14. Record: `deliverables/screen_phonons/screen_phonons.json`
(every number below, plus the modes, seeds and cell parameters behind it).

## Why this column exists

The governing requirement (`RESUME.md`, the consumer note of 2026-09-13 in
`docs/REFERENCE_DATA_REQUEST.md`) is low tan-delta and a high-frequency response. Neither side
computes loss. The two static quantities that bound it are the polar/antipolar lattice margin
and the lattice curvature. `docs/SCREEN.md` measures the first and finds it inside the fitted
potential's own 0.27 kcal/mol per monomer error bar for most chemistries. This document
measures the second, on the same all-trans polar zigzag and the same exactly-antipolar cell the
screen compares, with the phonon Hamiltonian of `docs/PHONONS.md`.

It is the static end of the loss question and no more: how stiff each lattice is against the
rigid-chain transverse motions that switching and relaxation would use. Everything
`docs/PHONONS.md` lists as absent from a Gamma-point curvature is absent here too: no torsional
stiffness about backbone bonds (twist modes are lower bounds), no LO-TO splitting, no
temperature, no rates, no linewidth, no field coupling. A Hessian at fixed cell also cannot see
a cell-shape instability, which is why the gamma-free cells are on the table beside the
gamma = 90 ones. And the potential is PVDF's (`pvdf-dft-valence-flux-born`, fitted to PVDF
single-chain PBE-D3 energies, forces and Born charges): on a chlorine or nitrile chemistry it is
the transfer assumption the screen already makes, and every number here inherits it.

## Construction

Four cells per chemistry, so that the screen's column and the model's own minima are both on
the table:

* **polar90 / anti90** are exactly the two cells `docs/SCREEN.md` compares: `pack()` at
  gamma = 90 (its lowest polar cell) and `fitting.antipolar_cell_exact` at gamma = 90, both on
  the screen's plain Ewald packer with the fitted charges. Their energies reproduce the
  screen's gap column to the quoted digits.
* **polar_free** is the lowest polar cell of `pack(gamma_free=True)`; **anti_free** is the
  antipolar subspace polished with gamma free, entered from seeds on either side of 90 deg
  (`antipolar_cell_exact` itself screens at gamma = 90 only; a gradient polish cannot leave
  90 deg because it is a symmetric point, hence the seeds).

Each cell is then re-polished on the phonon packer (valence, Born-fitted flux, induced
dipoles, Ewald): over the cell variables for the polar cells, and over the antipolar
subspace's own coordinates `(a, b[, gamma], phi1, dz)` with `phi2 - phi1` held for the
antipolar ones, which keeps the fixed-charge dipole exactly zero on any packer
(`fitting.antipolar_offsets`). The stretch minima are pinned to the built bond lengths
(`phonon.packer_with_built_bond_lengths`, the changed Hamiltonian `docs/PHONONS.md` labels
separately), the known-answer identity is checked on every cell, every atom is then let go at
fixed cell, and the Gamma-point Hessian is taken there with h = 1e-3 A and the acoustic sum
projected.

Two diagnostics travel with every cell. **"slid"** is the largest atom displacement when the
rigid cell's atoms were let go: a few hundredths of an Angstrom for a cell that was already a
minimum, tenths or more for a saddle, and the record keeps each chain's mean displacement
separately so a rigid slide is told from a pendant rearrangement. **The relaxed dipole per
monomer** (fixed charges) says whether a cell that was antipolar as a rigid packing stayed
antipolar once its atoms moved.

The rigid-chain cell here is not the deformable-path reference of `docs/PHONONS.md`, which also
relaxes torsions and backbone angles. Beta-PVDF's polar90 row is therefore a construction check
against that document, not a repeat of it: 34.5, 40.5, 49.2 cm^-1 here against 34.0, 40.7,
48.6 there, the same three rigid-chain modes in the same order.

## Results

### The margin, three ways (kcal/mol per monomer; verdicts against the screen's 0.27 error bar)

| chemistry | gap, gamma = 90 (the screen's) | verdict | gap, gamma free | verdict | anti_free gamma, density, closest contact | gap after all-atom relaxation, phonon packer |
|---|---:|---|---:|---|---|---:|
| pvdf | +1.021 | polar | -0.000 | not resolved | 61.7 deg, 2.116 g/cm^3, 2.64 A | -0.006 |
| pvdc | not measured: no all-trans repeat: no repeat, ideal or deflected (chain is not commensurate; cannot buil... | | | | | |
| cfe | +1.501 | polar | +0.071 | not resolved | 112.2 deg, 2.183 g/cm^3, 2.81 A | +0.014 |
| cdfe | +0.390 | polar | +0.391 | polar | 90.0 deg, 2.376 g/cm^3, 2.52 A | +0.444 |
| an | +0.046 | not resolved | +0.041 | not resolved | 90.8 deg, 1.540 g/cm^3, 2.38 A | +0.162 |
| vdcn | -0.645 | antipolar | -0.649 | antipolar | 89.5 deg, 1.781 g/cm^3, 2.34 A | -0.184 |
| fanome | not measured: FANOME's all-trans reference is an overlapping structure (methyl hydrogens of methoxy grou... | | | | | |
| pvf-cand | +0.205 | not resolved | +0.267 | not resolved | 85.8 deg, 1.698 g/cm^3, 2.50 A | +0.332 |
| trfe-cand | +0.051 | not resolved | +0.051 | not resolved | 90.8 deg, 2.357 g/cm^3, 2.57 A | -0.020 |
| vfcn-cand | +0.265 | not resolved | +0.046 | not resolved | 97.0 deg, 1.855 g/cm^3, 2.43 A | +0.177 |
| vclcn-cand | +0.134 | not resolved | +0.114 | not resolved | 89.0 deg, 1.964 g/cm^3, 2.60 A | -0.274 |

### The curvature, all atoms relaxed at fixed cell, stretch minima pinned (cm^-1)

RC = rigid-chain mode (both chains moving as rigid bodies), internal = a pendant or backbone motion; "slid" = largest atom displacement when the rigid cell's atoms were let go; the dipole is per monomer, fixed charges, after relaxation.

| chemistry | cell | lowest three optical | softest transverse | slid (A) | dipole/mon after (e.A) | imaginary |
|---|---|---|---|---:|---:|---:|
| pvdf | polar90 | 34.5, 40.5, 49.2 | 34.5 (RC) | 0.06 | 0.358 | 0 |
| pvdf | polar_free | 34.5, 40.5, 49.2 | 34.5 (RC) | 0.06 | 0.358 | 0 |
| pvdf | anti90 | 18.2, 22.5, 38.9 | 22.5 (RC) | 1.05 | 0.000 | 0 |
| pvdf | anti_free | 41.0, 49.4, 60.7 | 49.4 (RC) | 0.11 | 0.000 | 0 |
| cfe | polar90 | 22.0, 40.4, 66.0 | 22.0 (RC) | 0.05 | 0.386 | 0 |
| cfe | polar_free | 25.2, 42.6, 77.3 | 25.2 (RC) | 0.05 | 0.388 | 0 |
| cfe | anti90 | 23.0, 23.2, 42.4 | 23.0 (RC) | 0.22 | 0.000 | 0 |
| cfe | anti_free | 25.2, 47.4, 61.8 | 47.4 (RC) | 0.16 | 0.000 | 0 |
| cdfe | polar90 | 35.0, 46.0, 63.8 | 46.0 (RC) | 0.07 | 0.334 | 0 |
| cdfe | polar_free | 34.2, 46.3, 64.8 | 46.3 (RC) | 0.07 | 0.333 | 0 |
| cdfe | anti90 | 22.7, 39.6, 45.7 | 39.6 (RC) | 0.13 | 0.009 | 0 |
| cdfe | anti_free | 22.7, 39.6, 45.7 | 39.6 (RC) | 0.13 | 0.009 | 0 |
| an | polar90 | 63.9, 78.8, 83.6 | 78.8 (RC) | 0.34 | 0.287 | 0 |
| an | polar_free | 63.9, 78.8, 83.6 | 78.8 (RC) | 0.34 | 0.287 | 0 |
| an | anti90 | 54.8, 72.2, 88.7 | 54.8 (RC) | 0.57 | 0.254 | 0 |
| an | anti_free | 53.2, 71.3, 89.0 | 53.2 (RC) | 0.57 | 0.258 | 0 |
| vdcn | polar90 | 13.1, 14.6, 49.6 | 13.1 (internal) | 0.29 | 0.286 | 0 |
| vdcn | polar_free | 13.1, 14.6, 49.6 | 13.1 (internal) | 0.29 | 0.286 | 0 |
| vdcn | anti90 | 52.2, 66.7, 72.0 | 66.7 (RC) | 0.62 | 0.242 | 0 |
| vdcn | anti_free | 53.3, 66.6, 73.0 | 66.6 (RC) | 0.62 | 0.241 | 0 |
| pvf-cand | polar90 | 55.2, 63.3, 71.6 | 55.2 (RC) | 0.18 | 0.301 | 0 |
| pvf-cand | polar_free | 37.3, 53.7, 81.9 | 37.3 (RC) | 0.06 | 0.302 | 0 |
| pvf-cand | anti90 | 49.1, 55.2, 102.2 | 49.1 (RC) | 0.29 | 0.000 | 0 |
| pvf-cand | anti_free | 51.3, 52.0, 105.9 | 52.0 (RC) | 0.18 | 0.000 | 0 |
| trfe-cand | polar90 | 33.9, 44.4, 61.6 | 33.9 (RC) | 0.05 | 0.312 | 0 |
| trfe-cand | polar_free | 16.2, 42.0, 59.9 | 16.2 (RC) | 0.04 | 0.312 | 0 |
| trfe-cand | anti90 | 30.9, 46.8, 49.1 | 46.8 (RC) | 0.12 | 0.000 | 0 |
| trfe-cand | anti_free | 30.9, 46.6, 49.1 | 46.6 (RC) | 0.12 | 0.000 | 0 |
| vfcn-cand | polar90 | 36.9, 54.7, 65.9 | 54.7 (RC) | 0.58 | 0.243 | 0 |
| vfcn-cand | polar_free | 36.9, 54.7, 65.9 | 54.7 (RC) | 0.58 | 0.243 | 0 |
| vfcn-cand | anti90 | 23.6, 34.2, 65.1 | 34.2 (internal) | 0.30 | 0.000 | 0 |
| vfcn-cand | anti_free | 27.0, 44.6, 54.5 | 27.0 (internal) | 0.63 | 0.073 | 0 |
| vclcn-cand | polar90 | 40.2, 53.5, 57.4 | 53.5 (RC) | 0.53 | 0.218 | 0 |
| vclcn-cand | polar_free | 40.2, 53.5, 57.4 | 53.5 (RC) | 0.53 | 0.218 | 0 |
| vclcn-cand | anti90 | 29.0, 47.1, 61.8 | 47.1 (RC) | 0.78 | 0.267 | 0 |
| vclcn-cand | anti_free | 27.6, 46.2, 62.6 | 46.2 (RC) | 0.78 | 0.268 | 0 |

The record also carries, for every chiral chemistry (CFE, CDFE, AN and all four candidates), the package's standing warning that the chain built is the isotactic one and that syndiotactic and atactic chains are not expressible; for VFCN-cand and VClCN-cand, that reflection-with-reversal does not hold to 0.5 kcal/mol and the RIS fit was left unsymmetrised. PVDC has no all-trans repeat (ideal or deflected) and FANOME is not fitted (`docs/CHEMISTRY_EXTENSION.md`), exactly as in the screen.

## Reading, for the requirement

1. **The lattice-curvature column, polar all-trans cells** (softest transverse rigid-chain
   optical mode once every atom is relaxed at fixed cell; the lowest mode in brackets where it
   is something else): beta-PVDF 34.5; CFE 22.0 at gamma = 90, 25.2 at its gamma-free cell;
   CDFE 46.0 (35.0 is an axial slip); TrFE-cand 33.9 at gamma = 90 but **16.2** at its
   gamma-free cell (gamma 87.3, 0.6 kcal/mol per cell apart); PVF-cand 55.2 at gamma = 90,
   37.3 at its gamma-free cell (gamma 75.0, 0.37 kcal/mol per cell lower); VFCN-cand 54.7
   (36.9 axial); VClCN-cand 53.5 (40.2 axial); AN 78.8 (63.9 axial); VDCN **13.1 and 14.6**,
   and those two are not rigid-chain modes at all (rigid-chain share 0.28, carbon-dominated):
   they are nitrile pendant librations, the softest modes in the whole set. Read for the
   requirement: on this potential VDCN's polar lattice is the nearest to an instability and
   the instability is internal; CFE and TrFE have the softest rigid-chain transverse modes
   (16 to 25 cm^-1, 0.5 to 0.75 THz); PVDF and CDFE sit at 35 to 46; the bulky-pendant
   candidates are stiffest against rigid-chain motion but see point 2.
2. **Which cells are minima at all.** The polar cells of PVDF, CFE, CDFE, TrFE-cand and
   PVF-cand (gamma free) move by 0.04 to 0.07 A when let go: the rigid all-trans repeat is
   essentially the all-atom minimum and the frequencies belong to it. Every nitrile chemistry
   moves by 0.3 to 0.6 A, all of it pendant rearrangement (the chain means barely move), with
   two or three imaginary modes at the rigid cell that the rearrangement removes: the rigid
   repeat the screen packs is not near this potential's minimum for AN, VDCN, VFCN-cand or
   VClCN-cand, which is the same thing the campaign found for the nitrile finite-chain
   references (never at minima, `RESUME.md`). And once the pendants move, the *antipolar*
   cells of AN, VDCN and VClCN-cand carry 0.24 to 0.27 e.A per monomer of fixed-charge dipole,
   against 0.22 to 0.29 for their polar cells: **the rigid-chain polar/antipolar
   classification does not survive pendant relaxation for the nitrile chemistries.** Their
   "antipolar" column in `docs/SCREEN.md` is a statement about a rigid repeat that the
   potential itself will not keep.
3. **The margin, re-measured three ways.** At gamma = 90 the screen's column is reproduced to
   the quoted digits. With gamma free, PVDF's +1.021 becomes -0.000 and CFE's +1.501 becomes
   +0.071; TrFE-cand +0.051, VFCN-cand +0.046 and AN +0.041 are all inside the 0.27 error bar;
   only CDFE (+0.391, its antipolar minimum stays at gamma = 90) and VDCN (-0.649) keep a
   resolved sign, with PVF-cand (+0.267) on the bar. After every atom is relaxed on the phonon
   packer the same three fluorinated chemistries are degenerate again (PVDF -0.006, CFE
   +0.014, TrFE-cand -0.020), CDFE and PVF-cand stay polar (+0.444, +0.332) and the nitrile
   gaps are gaps between two polar states (point 2). For the requirement's discriminant, hard
   ferroelectric against relaxor, this potential therefore says: undecidable for PVDF, CFE
   and TrFE (the three chemistries industry uses), polar for CDFE and PVF-cand, and not a
   rigid-chain question for the nitriles.
4. **PVDF's degenerate antipolar competitor is a real, stiffer minimum.** The gamma = 90
   antipolar cell is a saddle (three imaginary rigid-chain modes at the rigid cell; the two
   chains slide 1.05 A in opposite directions when let go). The gamma-free one, a = 4.56,
   b = 9.64, gamma = 61.8, is a rectangular 4.56 x 8.50 A cell with chain 2 directly above
   chain 1 along the long axis instead of at the centre, antiparallel and half a repeat up.
   It moves 0.11 A when let go, has no imaginary modes, and its lowest modes (41.0, 49.4,
   60.7) are *stiffer* than beta's (34.5, 40.5, 49.2). So beta is not near a soft-mode
   instability on this potential; it has an antipolar competitor of the same energy at a
   different cell shape, reached by an in-plane chain slide, and whether the polymer can get
   there is a barrier question a Gamma-point curvature does not answer.
5. **An off-axis chain moment leaves an induced-dipole residual in the antipolar cell.**
   CDFE's moment lies 61 deg from its own x, so `phi2 = phi1 + 298.9` zeroes the fixed-charge
   dipole without being a crystal symmetry; the phonon packer's induced dipoles then leave
   0.0075 C/m^2 (against 0.113 in the polar cell) and the relaxed cell keeps 0.009 e.A per
   monomer. Recorded, not corrected: the screen's subspace is defined on the fixed charges.

None of this is a tan-delta. It is the static end of the loss question on a potential fitted
to PVDF, and the nitrile rows in particular are a statement about that transfer.

## What the gamma-free cells say about the screen's column

`DESIGN.md` section 6 already records that beta's polar and antipolar cells, distinguishable at
gamma = 90, become exactly degenerate when gamma is freed, and that the polarity verdict is
"a verdict about the parametrisation it was measured in". This screen sees the same thing from
inside the cell and extends it: the same is true of CFE (+1.501 at
gamma = 90, +0.071 with gamma free) and of VFCN-cand (+0.265 to +0.046), while CDFE's polar
verdict and VDCN's antipolar one survive the change. Of the screen's three resolved-polar
chemistries, two are gamma = 90 verdicts. Both `pack()` (`gamma_free=False` by default) and
`fitting.antipolar_cell_exact` (its grid rows carry `90.0` explicitly) measure at gamma = 90,
and `docs/SCREEN.md` states that protocol; re-measuring its polarity column with gamma free
on both branches is the screen-level change this implies, and it is not made here.

Two things found on the way, recorded for whoever touches the packer next. The packer's energy
is exactly periodic in the chain-2 axial offset `dz` for one or two repeats outside `[0, c)`
and drifts beyond that (1e-4 kcal/mol at three repeats, 1 kcal/mol at six, nonsense at nine;
measured on beta-PVDF under both sums). `pack.polish` leaves `dz` unbounded and wraps it
afterwards, so a polish that wanders far can be scored on a wrong energy; the tied polish here
wraps inside its objective (a VDCN seed once walked nine repeats out and reported an antipolar
cell 4.4 kcal/mol per monomer too low, at a cell the polarizable model could not evaluate).
That is a latent source fix, held with the other source changes under the merge hold. And
`antipolar_cell_exact`'s own polish is free in `(a, b)` and can return a cell just outside the
screen's bounds; a bounded polish started there projects onto the bound first, which for CDFE
raised the antipolar energy by 0.15 kcal/mol per monomer before it was caught.

## Reproducing this

```
set PYTHONPATH=src
python examples/screen_phonons.py --json screen_phonons.json          # every chemistry, ~1 h
python examples/screen_phonons.py --only pvdf                          # 2 min; the construction check
```
The relaxed-angle RIS fit (`docs/NITRILE_LANDSCAPE.md`) is most of the cost for CFE, CDFE, AN and
VDCN (four to seven minutes each); the packing, the two antipolar searches and the four
Hessians take about a minute per chemistry. `tests/test_screen_phonons.py` checks that the tied
polish never leaves the antipolar subspace, that a start outside the screen's bounds is not
projected uphill (which raised CDFE's antipolar energy by 0.15 kcal/mol per monomer before it
was caught), that the gamma-free search never returns above the gamma = 90 answer, and that the
summary reads the record.
