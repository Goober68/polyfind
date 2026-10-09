# Polar/antipolar margin of all-trans PVDF, periodic PBE-D3(BJ) (Quantum ESPRESSO)

The record behind `docs/POLAR_MARGIN_DFT.md`. Read that document's caveats before using any
number here: this is a zero-kelvin static lattice energy on one Hamiltonian, not a barrier, not
a free energy, and by itself it says nothing about loss.

Hamiltonian, as in the provider's QE references (`sarco/materials/gpu_bundle/periodic_reference/beta_pvdf/born/scf.in`):
PBE, SSSP 1.3.0 PBE precision pseudopotentials (`H_ONCV_PBE-1.0.oncvpsp.upf`,
`C.pbe-n-kjpaw_psl.1.0.0.UPF`, `F.oncvpsp.upf`), 90/360 Ry, `vdw_corr='grimme-d3'`,
`dftd3_version=4` (BJ), `dftd3_threebody=.false.`, fixed occupations. Quantum ESPRESSO 7.5
(conda-forge), run in WSL on 2 x 3 MPI ranks.

**Result** (per monomer, relaxed cells, 8x4x16 mesh = the provider's 4x8x16):
E(antipolar) - E(polar) = **-1.113 meV (-0.0257 kcal/mol)** for polyfind's gamma-free
antipolar cell (`anti_free`), D3 part -3.361 meV, PBE part +2.248 meV. The screen's gamma = 90
antipolar cell (`anti90`) relaxes into the same packing, -1.151 meV. A symmetry-free
re-relaxation finds the antipolar basin's lowest point at -1.160 meV (4x2x8 mesh); beta comes
back to itself. Converged to better than 0.05 meV in k-mesh (4x2x8 .. 12x6x24; the 2x1x4
pre-relaxation mesh gets the sign wrong and is kept only as a warning), cutoff
(90/360 .. 130/520 Ry) and relaxation, and the provider's beta start reaches the same minimum
to 0.001 meV. Three-body D3 at fixed geometry: -0.10 meV. Polyfind's own all-atom-relaxed
margin for the same pair: -0.0056 kcal/mol. Relaxed beta: 8.528 x 4.757 x 2.581 A against
Hasegawa's 8.58 x 4.91 x 2.56. Control: E(alpha) - E(beta) = -31.5 meV = -3.04 kJ/mol
(reduced protocol), inside published dispersion-corrected DFT. The QE log files and checks
behind every number are summarised in `raw_stages.json`.

## Files

| path | contents |
|---|---|
| `polar_margin_dft.json` | the compact record: per cell the relaxed lattice, chain orientation and registry, the energy and D3 term of every final SCF (Ry per monomer), the k2 acceptance (forces, stress); the margins (Ry, meV, kcal/mol per monomer) at every mesh and cutoff with their D3 and PBE parts; the symmetry-free checks; lattices against experiment; polyfind's own margins; QE version and pseudopotential md5/SHA-256 |
| `raw_stages.json` | every stage's parsed observables (energy, D3, largest force component, stress tensor, irreducible k-points, symmetry operations, BFGS steps, acceptance), from `collect` |
| `starts/` | the DFT starting geometries (extended xyz, Cartesian A, lattice rows): polyfind's four `screen_phonons` cells for PVDF, both the rigid-chain and the all-atom-relaxed placement (`pvdf_<cell>_{rigid,relaxed}.xyz`), the provider's beta cell from `born/scf.in` with its SHA-256 (`pvdf_polar_dipole_start.xyz`), the alpha control start (`pvdf_alpha_start.xyz`); `polyfind_starts.json` holds polyfind's energies and margins for the same cells |
| `relaxed/` | the final PBE-D3(BJ) geometries, one per cell (`pvdf_<cell>_pbe_d3bj.xyz`; alpha's is its 2x1x2-mesh relaxation), and the symmetry-free re-relaxations (`pvdf_<cell>_nosym_pbe_d3bj.xyz`) |
| `kmesh_test_provider_geometry.json` | the pre-run k-mesh test on the provider's own beta geometry (4x8x16 against 2x4x8 and 3x6x12, their axis order), which also records the QE forces (0.0156 Ry/bohr) and stress (-13.4 kbar) at that geometry; inputs in `inputs/kmesh_test_provider_geometry/` |
| `inputs/<cell>/<stage>.in` | every pw.x input that was run, byte for byte (paths are this machine's WSL paths) |
| `runner/` | `run_pw.sh` (one niced pw.x run), `drive_cell.sh` (stage chaining and acceptance), `worker.sh` (job queue), `queue_as_run.txt` (the order things ran in) |

Not committed: pw.x outputs (kept on the WSL side under `~/polar_margin/runs`, 1-3 MB each),
wavefunctions and charge densities (deleted after each stage).

Cells: `polar` (beta), `anti_free` (polyfind's gamma-free antipolar competitor), `anti90`
(the screen's gamma = 90 antipolar cell), `polar_dipole` (beta from the provider's geometry),
`alpha` (TGTG' control). Stages: `relax0` (2x1x4), `relax1` (4x2x8) vc-relax; `scf_k0`,
`scf_k1`, `scf_k15`, `scf_k2`, `scf_k3` at meshes 2x1x4, 4x2x8, 6x3x12, 8x4x16, 12x6x24
(short in-plane x long in-plane x chain; 8x4x16 is the provider's 4x8x16 in their axis
order); `scf_k1_e110`, `scf_k1_e130` at 110/440 and 130/520 Ry; `scf_k1_3body` with the D3
three-body term; `relax_nosym`, a symmetry-free vc-relax from a perturbed final geometry.
A stage's `meets_criteria` is the relaxation acceptance (every force component < 1e-4 Ry/bohr,
every stress component < 0.5 kbar); it is expected to be false for the coarse-mesh,
higher-cutoff and three-body SCFs, which are evaluated at the two-body 90/360 geometry.

## Reproduce

```text
set PYTHONPATH=src
python examples/polar_margin_dft.py starts --alpha        # the starts (two minutes, plus a minute for alpha)
# in WSL, with QE 7.5 and the three SSSP files in ~/polar_margin/pseudo:
NP=3 NK=3 runner/drive_cell.sh polar RELAX FINALCHECK      # and likewise every line of runner/queue_as_run.txt
python examples/polar_margin_dft.py collect                # -> raw_stages.json, inputs/, relaxed/
python examples/polar_margin_dft.py report                 # -> polar_margin_dft.json
```
