# Gamma-point phonons of the PVDF polymorphs on the fitted potential

*2026-09-13. Method in `src/polyfind/phonon.py`; report script `examples/crystal_phonons.py`;
tests `tests/test_phonon.py`. Motivated by the governing requirement (`RESUME.md`): low
tan-delta and high-frequency response. The intrinsic proxy a static lattice model can add is the
stiffness of the crystal's softest modes, and this is the first lattice dynamics in the package.*

## What is computed

The packer's energy is a function of seven cell parameters and one chain repeat; chain 2 is the
image of chain 1. `phonon.cell_energy_and_grad` assembles the same energy, term for term and in
the same units, for every atom of the cell independently (same-site column, valence bonds and
angles, chain-chain block and lateral images, Ewald with per-chain exclusion corrections, charge
flux re-solved per chain, induced dipoles, applied field) with an analytic all-atom gradient.
The Cartesian Hessian is a central difference of that gradient (h = 1e-3 A), mass-weighted with
`pack.MASS`, and diagonalised. Frequencies are in cm^-1; an imaginary mode is reported as a
negative wavenumber. The unit constant, sqrt(kcal/(mol A^2 amu)) = 108.5914 cm^-1, is derived
in the module from CODATA/SI values and rederived in a test.

**Known-answer check, mandatory.** At the packer's own placed coordinates `cell_energy` must
equal `packer.energy` to 1e-9 relative; it does to 3e-14 (beta), 1.5e-14 (alpha), 4.6e-14
(gamma), and for plain LJ/DSF, Ewald-only and one-chain packers. The folded all-atom gradient
matches `energy_and_grad`'s repeat gradient to 5e-13; the all-atom gradient matches central
differences to 1.4e-8; the gradient-difference Hessian matches energy second differences to
2e-6. The acoustic sum rule (three uniform translations) holds to 7e-8 of the largest Hessian
entry without being imposed; projecting it out moves optical modes by < 1e-10 cm^-1.

## What is not in it

* **No torsional stiffness about backbone bonds.** The packer evaluates its Fourier torsion
  term once from the chain's nominal dihedrals and adds it as a constant ("only shifts the
  value, never the gradient"); this module does the same so the known-answer identity holds.
  Chain-twist modes are therefore held only by 1-4 and longer nonbonded terms and by the
  neighbours, and are lower bounds. The rigid-chain translations and librations that carry the
  soft transverse modes below involve no torsion and are unaffected.
* **The reference states are not all-atom stationary points.** The deformable path relaxes cell,
  setting angles, line-group torsions and backbone angles, never bond lengths, and the fitted
  stretch minima are not the built bond lengths (C-F r0 1.467 vs built 1.350 A; C-C 1.458 vs
  1.528 A, `k` at its fit bound). The residual all-atom force at the beta reference is 41
  kcal/(mol A). A Hessian there shows imaginary rigid-chain transverse modes (beta: -134,
  -121 cm^-1). Pinning the stretch minima to the built lengths (stiffness kept; known-answer
  still 4e-15) leaves -69 and -20 cm^-1, so the mismatch is only part of it; the rest is the
  residual angle and nonbonded force at the constrained reference. Letting every atom relax at
  fixed cell (max displacement 0.057 A with pinned lengths) removes all imaginary modes and
  leaves exactly three zeros. The numbers quoted below are that state: all atoms relaxed at
  fixed cell, stretch minima pinned to the built bond lengths.
* Gamma point only; no LO-TO splitting (tinfoil Ewald, no non-analytic term); zero kelvin;
  no anharmonicity, so no linewidth and no loss.

## Results, `pvdf-dft-valence-flux-born` + induced dipoles + Ewald

Lowest optical modes, cm^-1, all atoms relaxed at fixed cell, stretch minima pinned. Character
from the mass-weighted eigenvector: T = transverse to the chain axis, A = axial, RC = rigid-chain
(both chains moving as rigid bodies).

| polymorph | atoms | imaginary | lowest six optical | character of the lowest |
|---|---:|---:|---|---|
| beta (TTTT) | 12 | 0 | 34.0, 40.7, 48.6, 90.8, 96.1, 205.2 | 34.0 T RC in-plane libration, 88% F; 40.7 A RC slip; 48.6 T RC along the polar axis |
| alpha (TGTG') | 24 | 0 | 39.0, 58.2, 59.2, 74.8, 85.0, 86.6 | 39.0 A; 58.2 T along x |
| gamma (T3GT3G') | 48 | 0 | 21.9, 26.9, 38.6, 42.7, 47.4, 52.0 | 21.9 T, F-dominated, not rigid-chain; 26.9 A |

*(2026-10-08: the alpha row is the packer's own TGTG' reference, which relaxes to space group
P2_1cn, the polar arrangement of delta-PVDF (form IV), not alpha's centrosymmetric P2_1/c. The
model's P2_1/c cell and its spectrum are in the last section; its lowest optical mode is 33.6
cm^-1, so item 3 of the reading below holds for the polar TGTG' cell only.)*

Beta, full spectrum (cm^-1): 0 0 0 34.0 40.7 48.6 90.8 96.1 205.2 206.1 235.7 241.6 311.7
322.6 341.4 343.0 443.9 453.1 611.0 633.6 634.3 664.1 821.8 840.6 906.3 959.1 963.9 964.0
1144.4 1144.5 1506.2 1507.9 2609.0 2609.7 2678.3 2679.5.

Sensitivity: h = 2.5e-4 or 4e-3 moves frequencies by < 0.02 cm^-1. Relaxing all atoms with the
fitted stretch minima instead (C-F stretched to 1.54 A, energy 16 kcal/mol lower) gives beta
41.9, 52.8, 73.8, 114.4, 116.5 cm^-1: the soft-mode ordering and scale survive, the values move
by 10 to 40 cm^-1. Gamma's reference-state Hessians show a 3e-5 relative asymmetry against 5e-7
elsewhere, consistent with a pair straddling the energy-shifted LJ cutoff; the sum rule is
unaffected.

## Reading, for the requirement

1. **The softest harmonic modes of the polar crystal are rigid-chain transverse librations and
   translations at 34 to 49 cm^-1, about 1 to 1.5 THz.** That is the observed harmonic scale
   and no more. It does not establish the absence of an intrinsic sub-THz response limit: a
   Gamma-point curvature gives no finite-temperature transition rate, no linewidth, no
   susceptibility or field coupling, no finite-q acoustic or homogeneous-strain response and
   no switching dynamics, and the absence of those mechanisms from a calculation is not
   evidence about them. (An earlier version of this paragraph classified any actuator-frequency
   limit as extrinsic; that was an overreach, pointed out by the producer, and is withdrawn.)
2. **Beta is not near a transverse instability in this potential once its atoms are at their
   minimum**, but the constrained reference sits on the unstable side of a soft rigid-chain
   direction (-69 cm^-1 with pinned lengths). That is a property of the packer's line-group
   constraint, which pins chain 2 to the cell centre; it says the in-plane chain-chain
   arrangement is soft, which is the same physics as the polar/antipolar margin the screen
   cannot resolve. A same-Hamiltonian DFT curvature of the bulk cell is the number that would
   settle how soft.
3. **Across polymorphs the softest mode is lowest for gamma (22 cm^-1) and highest for alpha
   (39 cm^-1)**, with beta between on its transverse modes. The ordering is what a mixed
   trans-gauche chain with a looser packing would give and is not yet checked against
   far-infrared or Raman data; that comparison is the next validation.
4. **What would change these numbers most**: a torsional stiffness term in the packer (raises
   twist modes, leaves the rigid-chain modes), a valence stretch fitted to bond lengths rather
   than to forces at fixed geometry, and Dipole's bulk periodic curvature for the like-for-like
   comparison; their finite-chain curvatures (0.002 to 0.009 eV/A^2 minimum internal curvature,
   AN/PVDF/VDCN 5/7/9) are chain-internal, not lattice, quantities and do not compare directly.

Nothing here is a tan-delta or a bandwidth. It is the static end of the loss question: how
stiff the lattice is against the motions that switching and relaxation would use. The
built-length-pinned stretch minima define a changed Hamiltonian; results under it are labelled
as such and are not to be mixed with fitted-stretch energies or modes.

## The relaxation's acceptance, 2026-09-14

`relax_all_atom` now returns a `Relaxation` rather than a bare geometry: the largest unrounded
residual force, L-BFGS-B's termination message and counts, up to three Newton polishes (a
finite-difference Hessian solved in the translation-free subspace, kept only when it lowers the
force) and `converged`, true only when the force is finite and within the requested `gtol`.
A Hessian consumer must require `converged`; `examples/crystal_phonons.py` and
`examples/screen_phonons.py` do. Dipole's review of the screen record found a third of its
cells a few 1e-6 above the tolerance because the old function discarded termination; the
criterion is unchanged, the acceptance is now the function's own.

## Across the screen, 2026-09-14

`docs/SCREEN_PHONONS.md` runs the same Hessian on the all-trans polar and exact antipolar
cells of every screened chemistry, at gamma = 90 and with gamma free, from
`examples/screen_phonons.py`. Beta's polar cell there (rigid-chain packing re-polished on this
Hamiltonian, no torsion or angle relaxation) gives 34.5, 40.5, 49.2 cm^-1 against the 34.0,
40.7, 48.6 above, the construction check. The gamma = 90 antipolar beta cell is a saddle whose
chains slide 1.05 A when let go; the gamma-free one is a true minimum, stiffer than beta and
degenerate with it.

## Rigid chain rotation as a switching path: why it is not a loss proxy (2026-10-02)

The next static quantity the loss requirement suggests after curvature is a switching barrier
and, from it, an intrinsic coercive field and hysteresis energy (`~4 P_s E_c` per cycle). It
was tried on beta-PVDF (`pvdf-dft-valence-flux-born`, induced dipoles, Ewald) by rotating the
two chains about their own axes from the polar state to the reversed one, and it does not give
a usable number. Recorded so it is not re-run.

* **Clamped cell, both chains in phase**: polar and reversed states degenerate as they must
  (`P` +0.1434 to −0.1434 C/m²); barrier **6.05 kcal/mol per monomer at 40°**, a steric peak.
  With a field against `P` the polar state is still a local minimum at 1 V/Å (10 GV/m): near
  the minimum a rigid rotation changes the dipole only at second order, so the restoring
  curvature wins at any sensible field.
* **Cell relaxed at each angle** (`a, b, gamma, dz` polished, continuation in angle): in phase
  1.66 kcal/mol per monomer at 105° with `b` on its 9.91 Å bound; with the bound widened to
  14 Å, 1.50 at 138° with `gamma` on its 120° bound and the path ending 1.15 kcal/mol above the
  reversed minimum, i.e. the continuation tracks a branch rather than a minimum-energy path.
  One chain only: 2.25 kcal/mol per monomer at 87°, ending in the antipolar cell, degenerate
  with polar (−0.008), consistent with section "Across the screen".
* **Coercive field along these paths** (field against `P`, one-dimensional spinodal
  `max U'(x) / (−dmu_P/dx)`): 14-17 GV/m. That is the rotational curvature at the minimum
  divided by the dipole -- the libration stiffness already measured above -- not a switching
  field; it is tens of times larger than reported intrinsic coercive fields of ferroelectric
  P(VDF-TrFE) films (recalled at ~0.5 GV/m, not verified here) and larger still than bulk.

**Why.** Switching in a real crystal proceeds by nucleated kinks and domain walls that rotate
segments of a chain, not whole infinite chains at once, and the per-monomer barrier of an
infinite rigid rotation is not the barrier of that process. A one-monomer-repeat, two-chain
cell cannot hold a kink. A switching or hysteresis proxy needs a supercell along the chain with
torsional freedom (the packer's torsion term is a constant today, section 2) and a proper
path method; a small-signal dielectric loss needs dynamics. Neither is a static-curvature
calculation, and at device frequencies both are expected to be dominated by the
non-crystalline phase, which `polyfind.film` treats as passive.

## Far-infrared activity of the alpha cells against the measured lines (2026-10-08)

*Method in `src/polyfind/infrared.py`; report script `examples/alpha_far_ir.py`; tests
`tests/test_infrared.py`; record `deliverables/alpha_far_ir/`. This is the first comparison of the
package's lattice dynamics with a measured spectrum. Nothing was fitted or tuned for it. Both
Hamiltonians are reported: the shared branch (torsion a constant of the nominal dihedrals, as in
every section above) and the merged `vector-repeat-merge` tree at 8217385 (Dipole's Cartesian
torsion in the Hessian). Everything else is section "Results": `pvdf-dft-valence-flux-born`,
induced dipoles, Ewald, stretch minima pinned to the built bond lengths, every atom relaxed at
fixed cell.*

### What is computed

A Gamma mode `k` with mass-weighted eigenvector `e_k` carries the cell dipole derivative
`d_k = sum_i Z*_i e_{k,i} / sqrt(m_i)`. Its infrared intensity is `|d_k|^2`, and the direction of
`d_k` is the polarization of the light that excites it. `Z*` is the model's own Born tensor at the
relaxed geometry, taken as the mixed second derivative of the same energy with respect to a
uniform applied field `F` and the displacement: `Z*_{i,ab} = -(1/EV_TO_KCAL) dg_{ib}/dF_a`, with
`g` the all-atom gradient of `phonon.cell_energy_and_grad`. The field enters that energy as
`-(sum q r) . F` with the fluxing charges and through the induced dipoles' permanent field. So
static charge, charge flux and electronic screening are all in `Z*`, and on either branch it is
the derivative of exactly the energy whose Hessian gives the modes.

Absolute intensities are `N_A / (12 eps0 c^2) |dmu/dQ|^2`, which is 974.9 km/mol per e^2/amu,
per mole of cells (four monomers). The comparison below uses only ratios and directions. The
space group is found from the relaxed coordinates: atoms matched modulo the lattice, tolerance
1e-3 A. Each mode's characters under it give its irreducible representation and the axes along
which symmetry allows its dipole.

**Checks, on all four calculations (two cells, two Hamiltonians):**

* The field route agrees to 7.4e-9 e with central differences of the cell dipole (each chain's
  charges re-derived, induced dipoles re-solved).
* Ten times the field step moves `Z*` by at most 4.5e-13 e: the energy is exactly quadratic in
  the field at fixed nuclei.
* The acoustic sum `|sum_i Z*_i|` is at most 4.7e-13 e, against a largest `|Z*|` of 1.90 e. The
  three translations carry at most 1.4e-27 of the strongest intensity.
* `Z*` respects the found operations to 1.5e-6 e or better, and every optical mode's character
  is +-1 to within 1.3e-7.
* Where the cell has an inversion, the strongest inversion-even mode carries at most 1.5e-10 of
  the strongest intensity. No mode has more than 2.1e-9 of it along a symmetry-forbidden axis.

The model's clamped-ion dielectric tensor, from the field curvature of the same energy, is
2.11/2.13/2.34 (x/y/z) for the antipolar cell and 2.13/2.16/2.36 for the polar cell. The tests
add a closed form (a diatomic's optical mode carries `q^2/mu`), `Z* = qI` for fixed charges,
agreement with `born.born_charges` at the placed geometry, and mutual exclusion in a C2h alpha
placement.

### The packer's alpha cell is not alpha's space group

**The packer's TGTG' reference is delta-like, not alpha.** This is the cell the "alpha" row
above and `examples/crystal_phonons.py` use: chain 2 flipped, setting angles 90 and 270 deg.
Relaxed, it has space group **P2_1cn**: a two-fold screw along a (x), a c-glide normal to b and
an n-glide normal to c, with atom-match errors of 5e-8 A or less. That is point group C2v, with
no inversion and a cell dipole along a of 1.033 e.A (|P| = 0.0795 C/m^2; 0.0783 with torsion).
P2_1cn is the space group of delta-PVDF (form IV, `docs/REFERENCES.md` 1.2), not alpha's P2_1/c.
`DESIGN.md` recorded long ago that the flip cancels the two chains' axial dipoles but not their
transverse ones (0.079 C/m^2), and here it shows in the phonons. The alpha row above is the
model's polar TGTG' cell. Its numbers stand as computed, but they are not alpha's.

**The model's alpha analogue is the antipolar cell.** It is the flip-and-equal-setting-angle
subspace of `fitting.antipolar_cell`, which has zero dipole for this helix by construction. It
was built in three steps:

1. The polar reference's chain was placed in that subspace.
2. It was polished over `(a, b, phi, dz)` from `antipolar_cell`'s two starts.
3. Cell and chain shape were relaxed together on the same deformable path
   (`mechanics.relax_reference_deformable`, with phi1 and phi2 free separately).

It stays in the subspace (phi1 = phi2 = 90.000 deg) and comes back **P2_1/c**: an inversion
(error 7e-8 A), a two-fold screw along b and a c-glide normal to b. That is alpha's group, in
the setting its cell is quoted in. The cell is a = 5.032, b = 9.019, c = 4.657 A, against the
measured 4.96, 9.64, 4.62 A (+1.4, -6.4, +0.8%).

With every atom relaxed at fixed cell it is a minimum on both Hamiltonians: no imaginary mode
and three zeros. It lies above the polar cell by:

* 0.072 kcal/mol per monomer at the cell level (fitted stretch);
* 0.030 kcal/mol per monomer (0.029 with torsion) after the all-atom relaxation with pinned
  stretch.

Both are inside the model's 0.27 kcal/mol per monomer resolution. The model therefore does not
choose between the alpha and delta packings. The sample's identity has to come from outside,
and this antipolar cell is the like-for-like one for alpha's lines.

### Modes below 200 cm^-1

**The P2_1/c (alpha) cell.** The intensity ("IR rel.") is relative to each calculation's
strongest mode below 200 cm^-1: 10.50 km/mol per cell without torsion and 12.86 with it. The
strongest mode of the whole spectrum is at 919.9 / 937.0 cm^-1 and carries 1024 / 1000 km/mol. "Along chain" is the share of `|d|^2` on the chain axis. "LO est." is
`infrared.lo_shift`, a first-order estimate with the wavevector along the mode's own dipole. The
character column comes from the `Mode` decomposition, with its rigid-chain share split here into
translations and libration about the chains' own axes. Modes are matched across the two
Hamiltonians by eigenvector overlap, which is 0.96 or more for every row.

| no torsion | torsion | irrep | IR rel. | along chain | LO est. | character |
|---:|---:|---|---|---|---|---|
| 33.6 | 33.1 | Ag | dark (Raman) | | | rigid-chain translation along the chain axis |
| 47.8 | 49.1 | Ag | dark (Raman) | | | rigid-chain translation along a |
| 57.1 | 56.4 | Bg | dark (Raman) | | | rigid-chain libration, co-rotating (0.99) |
| **67.5** | **68.1** | **Au** | **0.237 / 0.174** | 0.00 (along b) | 69.5 / 69.9 | **rigid-chain libration, counter-rotating (0.98)** |
| **72.5** | **88.8** | **Bu** | **1.000 / 1.000** | 0.25 / 0.22 (rest along a) | 80.1 / 96.4 | **internal** (rigid-chain 0.00), F 0.71 |
| 85.7 | 96.3 | Ag | dark (Raman) | | | internal, F 0.74 |
| 86.5 | 93.4 | Bg | dark (Raman) | | | rigid-chain translation along b |
| 162.1 | 167.9 | Bg | dark (Raman) | | | internal |
| 163.3 | 168.0 | Au | 0.009 / 0.018 | 0.00 | 163.3 / 168.0 | internal |
| 167.1 | 172.5 | Ag | dark (Raman) | | | internal |
| 177.3 | 175.3 | Au | 0.184 / 0.120 | 0.00 (along b) | 177.9 / 175.7 | internal, F 0.78 |
| 179.0 | 177.4 | Bg | dark (Raman) | | | internal |
| 189.2 | 193.9 | Bu | 0.347 / 0.269 | 0.93 / 0.89 | 190.2 / 194.8 | internal, F 0.78 |

Below 160 cm^-1 there are exactly two infrared-active modes. One is the counter-rotating
libration of the two chains (Au, polarized along b). The other is an internal, F-dominated mode
(Bu) with no rigid-chain content. The co-rotating libration is Bg, Raman-only. Torsion moves
the internal Bu mode by +16.2 cm^-1, the largest shift in the window, and the libration by
+0.6.

**The P2_1cn (polar, delta-like) cell, for the record.** These are the modes the earlier
"alpha" row describes. Eigenvector overlap across the Hamiltonians is 0.91 or more.

| no torsion | torsion | irrep | IR rel. | along chain | character |
|---:|---:|---|---|---|---|
| 39.0 | 39.2 | A1 | 0.001 / 0.000 | 0.00 | rigid-chain translation along the chain axis |
| 58.2 | 63.5 | B (z) | 0.026 / 0.015 | 1.00 | rigid-chain translation, mostly along a |
| 59.2 | 58.1 | A2 | dark | | rigid-chain libration, counter-rotating |
| 74.8 | 76.9 | B (y) | 0.295 / 0.244 | 0.00 | rigid-chain libration, co-rotating |
| 85.0 | 86.7 | A2 | dark | | rigid-chain translation along b |
| 86.6 | 99.3 | A1 | 1.000 / 1.000 | 0.00 (along a) | internal, F 0.75 |
| 94.4 | 105.4 | B (z) | 0.240 / 0.208 | 1.00 | internal |
| 173.1 | 171.9 | B (y) | 0.242 / 0.170 | 0.00 | internal |
| 186.0 | 189.5 | B (z) | 0.390 / 0.311 | 1.00 | internal |

The modes at 167.4/174.6, 170.2/174.9, 176.4/175.2 and 177.2/183.5 are dark or carry 0.02 or
less. In this cell the libration selection is inverted: the counter-rotating libration is A2 and
dark, and the co-rotating one is infrared-active.

### Measured lines

Access levels: "full text" means read in full here, or by the literature search where marked;
"abstract" means the abstract only; "second-hand" means quoted by a paper whose full text was
read.

| cm^-1 | T | technique | published assignment | source, and what was read |
|---:|---|---|---|---|
| 53 | room temperature | far-IR | "rotatory lattice mode", from the P2_1/c factor-group analysis, crystallinity and temperature dependence | Rabolt & Johnson 1973 (abstract); Chamorro-Posada 2016 section 3.2 (full text) |
| 60 | 90 K | far-IR | the 53 line, cooled from 300 K | Chamorro-Posada 2016 citing Rabolt & Johnson (second-hand). R&J's abstract reports the cooling shift without numbers and ties it to unit-cell contraction. |
| 78.93 | not stated | THz-TDS | "third resonance" | Chamorro-Posada 2016 citing Mori et al. 2015 (second-hand; Mori's full text not accessible) |
| 85 | not stated (spectra 13-373 K) | far-IR | intramolecular, single TGTG' chain | Rabolt & Johnson 1973 (abstract) |
| 100 | room temperature | far-IR, one of the "two main absorption lines" | intramolecular, single chain (R&J) | Chamorro-Posada 2016 citing R&J (second-hand). R&J's abstract prints **102**. |
| 175 | not stated | far-IR | intramolecular, single chain | Rabolt & Johnson 1973 (abstract) |

Further information, of more limited weight:

* **Chamorro-Posada's own computation.** He computes (PM6, 2x1x2 cell) a lattice libration with the two
  chains out of phase, which he likens to Kobayashi's R_c mode, and assigns it to the 60 cm^-1
  line. He
  says the R_c^0 libration "is not optically active" in the perfect crystal. He quotes
  Kobayashi, Tashiro & Tadokoro's force-constant value for the 100 line as 94 cm^-1 (all read
  here).
* **Polarization.** No per-line dichroism for form II below 200 cm^-1 was found. Luongo 1972
  (abstract) reports "strong perpendicular dichroism in a number of absorptions" between 1000
  and 50 cm^-1 in oriented film, naming no band. Meng et al. 2020 (full text, read by the
  literature search) see the ~1.5-1.6 THz loss-tangent feature unchanged when the film is
  rotated by 90 deg. The model's polarizations are therefore untested by anything read.
* **Search-result previews only, not used in any comparison below.**
  * Mori et al. 2015: peaks at 1.60, 2.36, 3.04 and 5.31 THz (53.4, 78.7, 101.4 and 177.1 cm^-1,
    converted here). The first and third are "clear", the second and fourth "small".
  * Kobayashi et al. 1975: librational lattice modes at 70, 53 and 84 cm^-1 for forms I, II and
    III.
  * Okada & Ando 2016: librational lattice modes at 53 and 65 cm^-1 for forms II and II_p.
* **Raman.** Nallasamy & Mohan 2005 (full text, read by the literature search; their FTIR stops
  at 200 cm^-1) give Raman 170 (Ag) and 177 (Bg), both very weak, assigned to CF2 torsion. The
  model's inversion-even modes nearest them are 167.1 Ag and 179.0 Bg (172.5 Ag and 177.4 Bg
  with torsion).

### The comparison

**On the P2_1/c cell, the like-for-like one:**

1. **The 53 cm^-1 line (60 at 90 K).** The model's counter-rotating libration lands at
   **67.5 cm^-1 (68.1 with torsion)**. It is infrared-active (Au), polarized perpendicular to the
   chain along b, at 0.24 (0.17) of the window's strongest intensity. Its character is the
   published one (R&J's rotatory lattice mode, Chamorro-Posada's out-of-phase libration), and the model's
   selection rule agrees with his statement that the in-phase libration is optically inactive:
   the model puts that libration at 57.1 Bg. Measured against the 90 K value the model is
   **+7.5 cm^-1 (+8.1), +12.5% (+13.5%)**. Against the room-temperature 53 it is +14.5 (+15.1).
2. **The 100 cm^-1 line.** No infrared-active model mode lies within 10 cm^-1 of it on either
   Hamiltonian. The model's strongest mode below 160 cm^-1 is the internal Bu mode at
   **72.5 cm^-1 (88.8 with torsion)**, mostly perpendicular to the chain (0.25 / 0.22 of `|d|^2`
   along it). It is **27.5 (11.2) cm^-1 below 100**, and 29.5 (13.2) below R&J's 102.
3. **The ~79 and 85 lines.** The same Bu mode is the only candidate: 72.5 is -6.4 from 78.93 and
   -12.5 from 85; 88.8 is +9.8 and +3.8. The measurement has two or three lines between 78 and
   102, and what was read does not settle whether 78.93 and 85 are one band. The model has one
   infrared-active mode there, so at most one of those lines is matched.
4. **The 175 line.** An internal Au mode at 177.3 (175.3), perpendicular to the chain, at 0.18
   (0.12): +2.3 (+0.3). The parallel-polarized Bu at 189.2 (193.9, 0.35 / 0.27) has no partner
   among the lines read.
5. **Relative intensities.** The measured "two main lines" are 53 and 100. In the model the
   brightest low-frequency mode is the internal Bu, and the libration carries a quarter (no
   torsion) to a sixth (torsion) of its intensity.

**On the polar P2_1cn cell**, the modes the earlier tables call alpha, the picture differs. With
torsion, its three infrared-active modes between 70 and 110 cm^-1 sit at 76.9 (co-rotating
libration), 99.3 (internal, strongest) and 105.4 (internal, parallel to the chain). Those are
-2.0 from 78.93 and -0.7 and +5.4 from 100. But the cell has nothing bright below 70 cm^-1:
58.2 / 63.5 cm^-1 carry 0.026 / 0.015, and its counter-rotating libration is dark. So it does
not reproduce the 53/60 line, and it has the wrong space group for the sample. Its closeness at
79 and 100 is not evidence for the model's alpha. What it does measure is how far the packing
alone moves the internal mode: 72.5 against 86.6 without torsion, and 88.8 against 99.3 with
it.

**Plainly.** On alpha's own cell, the model's two strongest low-frequency infrared modes land at
67.5 and 72.5 cm^-1 without torsion and at 68.1 and 88.8 with it. The libration is the one with
the 53-line's assignment, and it misses the 90 K anchor by +7.5 / +8.1 cm^-1. The brightest
mode misses the 100 line by -27.5 / -11.2 cm^-1. Nothing lands near 100. The 79-85 region gets
one mode where two or three lines are reported.

### What is and is not a fair comparison

* **Temperature.** The model is a 0 K harmonic calculation at its own static cell. The lines are
  room-temperature measurements, except the 60 cm^-1 value at 90 K, which was read second-hand.
  That 300 K to 90 K shift of the libration is +13%, the same size as the model's miss against
  the 90 K value. The 90 K value is the better anchor and is still not 0 K.
* **Cell.** The antipolar cell's b is 6.4% shorter than measured. R&J's abstract ties the
  libration's cooling blue-shift to cell contraction. The model's frequencies at the measured
  cell were not computed.
* **No LO-TO splitting.** The Hessian is the tinfoil, transverse one. The first-order LO
  estimate puts the bright Bu mode's LO partner 7.6-7.7 cm^-1 above its TO (80.1, or 96.4 with
  torsion) and the libration's about 2 cm^-1 above. Every other mode moves by at most 1.0 cm^-1. A
  semicrystalline film's absorption maximum can fall between TO and LO, depending on crystallite
  shape and surroundings, which nothing here models. So for that one mode the comparable model
  number spans its TO to about 8 cm^-1 higher. That would bring the torsion value to within 4-6
  cm^-1 of 100/102, but this is an estimate, not a result.
* **Pinned-stretch Hamiltonian.** As in every section above, the stretch minima are pinned to the
  built bond lengths, which makes it a changed Hamiltonian. On beta, relaxing with the fitted
  stretch minima instead moved the soft modes by 10 to 40 cm^-1 (section "Results", Sensitivity).
* **Torsion decides the internal mode.** The bright Bu mode moves 16 cm^-1 between the two
  Hamiltonians, so its position is as uncertain as the torsion term's curvature. That term was
  fitted to conformer energies, never to a curvature. The libration moves 0.6 cm^-1.
* **What the crystal model leaves out.** It is a perfect crystal at Gamma: no amorphous phase, no
  defects, and none of the disorder-activated "optically inactive" absorption that Chamorro-Posada
  invokes for 78.93. The intensities are harmonic and linear in the coordinate, from a classical
  charge model whose transverse Born charges are within 0.19 e of DFPT on beta (`RESUME.md`).
* **Assignments.** The species of the measured lines (Au or Bu) and their dichroism were not
  available in anything read. The assignments used here are the published ones: lattice
  libration for 53, intramolecular for 85, 102 and 175.

### Reading

The first experimental check of the lattice dynamics is partial.

**What the model gets right.** The rigid-chain libration is the mode class the sections above
call the softest. On alpha's cell the model gets its symmetry right: only the counter-rotating
libration is infrared-active in P2_1/c, and that is the published assignment of the 53 cm^-1
line. It places the mode within 8 cm^-1 (13%) of the 90 K value, on both Hamiltonians.

**What it does not get.** The intramolecular region between 79 and 102 cm^-1 is not reproduced.
The model has one bright mode there, 11 to 28 cm^-1 below the 100 line, and its position is set
mostly by the torsion term.

**What it implies for the earlier record.** The earlier "alpha" phonons are of the polar delta
arrangement. That record, and the cross-polymorph ordering drawn from it, should be read with
that label.

**What would sharpen this:**

* the full texts of Kobayashi et al. 1975, for species and dichroism, and of Mori et al. 2015,
  for line temperatures and widths;
* a measured far-IR spectrum of oriented alpha film, to test the polarizations;
* the model's spectrum at the measured cell, as a sensitivity, not a fit.

### Sources

**Read in full here:**

* P. Chamorro-Posada, arXiv:1604.03919v1 (2016), section 3.2, Table 2, Fig. 2 caption and
  references [34]-[36]. Published as J. Appl. Spectrosc. 85, 552 (2018); only the metadata of
  that version was seen.

**Abstract read here** (via OpenAlex):

* J. F. Rabolt and K. W. Johnson, "Low frequency vibrations in polyvinylidene fluoride (form
  II)", J. Chem. Phys. 59, 3710-3712 (1973), doi:10.1063/1.1680540.

**Read by the literature search:**

* Full text: N. Meng et al., J. Mater. Chem. C 8, 16436 (2020), doi:10.1039/D0TC04310A;
  P. Nallasamy and S. Mohan, Indian J. Pure Appl. Phys. 43, 821 (2005); A. Itoh, PhD thesis,
  arXiv:1405.5889 (tables 5.1-5.3, computed values only).
* Abstract only: Luongo, J. Polym. Sci. A-2 10, 1119 (1972); Latour & Rahmoune, Ferroelectrics
  159, 245 (1994).

**Not accessible** (paywall or bot block):

* M. Kobayashi, K. Tashiro, H. Tadokoro, Macromolecules 8, 158-171 (1975),
  doi:10.1021/ma60044a013;
* T. Mori et al., J. Mol. Struct. 1090, 93-97 (2015), doi:10.1016/j.molstruc.2014.12.004;
* R&J's full text;
* Okada & Ando, Polymer 86, 83 (2016).

**Not about PVDF.** The PCCP article offered as a source (doi:10.1039/C7CP06225G) is Krylov et
al., "Raman spectroscopy studies of the terahertz vibrational modes of a DUT-8 (Ni)
metal-organic framework", PCCP 19, 32099 (2017). Its Crossref metadata and reference list
contain no PVDF.

### Reproduce

```
D=deliverables/alpha_far_ir
python examples/alpha_far_ir.py --cell antipolar --out $D/antipolar_shared.json
python examples/alpha_far_ir.py --cell polar --out $D/polar_shared.json
PYTHONPATH=<vector-repeat-merge>/src python examples/alpha_far_ir.py --cell antipolar --out $D/antipolar_torsion.json
PYTHONPATH=<vector-repeat-merge>/src python examples/alpha_far_ir.py --cell polar --out $D/polar_torsion.json
python examples/alpha_far_ir.py --compare $D/antipolar_shared.json $D/antipolar_torsion.json $D/polar_shared.json $D/polar_torsion.json
```

Each computing run takes about 10 minutes. Against a tree that predates `polyfind.infrared`, the
example puts this checkout's copy on the package path, and the record names the file and its
SHA-256.
