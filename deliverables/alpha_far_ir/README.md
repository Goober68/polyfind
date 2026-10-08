# Far-infrared activity of the alpha-PVDF cells

The record behind the last section of `docs/PHONONS.md` (2026-10-08). There are four files, two
cells times two Hamiltonians, all from `examples/alpha_far_ir.py`:

| file | cell | Hamiltonian |
|---|---|---|
| `antipolar_shared.json` | flip + equal setting angles, relaxed to P2_1/c (C2h): the model's alpha | shared branch, torsion a constant of the nominal dihedrals |
| `antipolar_torsion.json` | the same | `vector-repeat-merge` at 8217385, Cartesian torsion in the Hessian |
| `polar_shared.json` | the packer's own TGTG' reference, P2_1cn (C2v), polar: delta's space group | shared branch |
| `polar_torsion.json` | the same | `vector-repeat-merge` at 8217385 |

## What each record holds

* **Provenance.** The `polyfind` path, git head and branch, the count of uncommitted paths, and
  the SHA-256 of `phonon.py`, `pack.py`, `born.py`, `forcefield.py` and the `infrared.py`
  actually used. The torsion tree predates the module, so the example loads this checkout's
  copy.
* **Geometry.** The cell parameters, the relaxed placed coordinates and the lattice.
* **Relaxation.** The relaxation's acceptance.
* **Spectrum.** Every Gamma frequency.
* **Born tensors.** The tensors (field route), their acoustic sum, the field-versus-dipole
  agreement and their symmetry error.
* **Electrostatics.** The clamped-ion dielectric tensor and the cell dipole.
* **Symmetry.** The symmetry operations found, the point group of the relaxed and unrelaxed
  cell, and the inversion-even and symmetry-forbidden intensity leaks.
* **Per mode, below 320 cm^-1.** Each mode stores:
  * the irrep, parity and symmetry-allowed dipole axes;
  * the dipole derivative (e/sqrt(amu)), km/mol per cell and the polarization shares;
  * the first-order LO estimate;
  * the `phonon.Mode` shares;
  * the mass-weighted eigenvector, which `--compare` uses to match modes across Hamiltonians
    and to split the rigid-chain share into translation and libration.

## Reproduce

Run from the repository root, about 10 minutes per computing run:

```
D=deliverables/alpha_far_ir
python examples/alpha_far_ir.py --cell antipolar --out $D/antipolar_shared.json
python examples/alpha_far_ir.py --cell polar --out $D/polar_shared.json
PYTHONPATH=<vector-repeat-merge>/src python examples/alpha_far_ir.py --cell antipolar --out $D/antipolar_torsion.json
PYTHONPATH=<vector-repeat-merge>/src python examples/alpha_far_ir.py --cell polar --out $D/polar_torsion.json
python examples/alpha_far_ir.py --compare $D/antipolar_shared.json $D/antipolar_torsion.json $D/polar_shared.json $D/polar_torsion.json
```

## Before using any number here

* The model is a 0 K harmonic, Gamma-only, TO-only calculation (no LO-TO) of a perfect crystal.
* The stretch minima are pinned to the built bond lengths.
* The measured lines are the ones listed in the example's `MEASURED`, with the access level of
  each source stated in `docs/PHONONS.md`.

Read that section's caveats before quoting anything from this record.
