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

1. **The softest intrinsic modes of the polar crystal are rigid-chain transverse librations and
   translations at 34 to 49 cm^-1, about 1 to 1.5 THz.** Nothing intrinsic to a well-ordered
   beta lattice limits the response below that. A frequency limit at kHz to MHz, if it exists,
   is extrinsic: the amorphous fraction, domain walls, defects and conduction, none of which
   this model contains.
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
stiff the lattice is against the motions that switching and relaxation would use.
