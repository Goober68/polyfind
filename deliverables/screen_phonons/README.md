# Lattice-curvature screen record

`screen_phonons.json` is the complete record behind `docs/SCREEN_PHONONS.md`: for every
screened chemistry, the screen's gamma = 90 polar and antipolar cells and their gap, the
gamma-free cells and seeds, each cell re-polished on the phonon packer (parameters, energy,
polarization, density, closest interchain contact), and for each cell the rigid-cell and
all-atom-relaxed Gamma-point spectra (lowest optical modes, the three lowest modes' axis,
transverse, rigid-chain and per-element shares, imaginary-mode counts, acoustic-sum residual,
chain mean displacements and relaxed dipole per monomer). Timings and every warning are kept.

Reproduce with `PYTHONPATH=src python examples/screen_phonons.py --json screen_phonons.json`
(about an hour; `--only pvdf` is two minutes). The potential is `pvdf-dft-valence-flux-born`
with induced dipoles and Ewald, stretch minima pinned to the built bond lengths; it is PVDF's
potential on every chemistry. Read the document's caveats before using any number here.
