# The polar/antipolar margin of all-trans PVDF from periodic PBE-D3(BJ)

*2026-10-08/09. Script `examples/polar_margin_dft.py`; record `deliverables/polar_margin_dft/`
(`polar_margin_dft.json`, every pw.x input, the relaxed geometries, the pseudopotential hashes).
Computed on this machine (Quantum ESPRESSO in WSL, 6 MPI ranks, no GPU) because the producer
who would normally supply it has been inactive for three weeks. It is the "same-Hamiltonian
PVDF polar/antipolar pair" the consumer notes of 2026-09-13 and 2026-09-14 in
`docs/REFERENCE_DATA_REQUEST.md` ask for, with the gamma-free antipolar cell as the competitor.*

## Answer first

On the provider's Hamiltonian (PBE-D3(BJ), two-body, SSSP 1.3.0 precision, 90/360 Ry), with
every cell fully relaxed (variable cell, all atoms) and the energies taken from fresh-basis
SCFs at the provider's 4x8x16 density:

| pair | E(antipolar) - E(polar), per monomer | D3 part | PBE part |
|---|---:|---:|---:|
| `anti_free` - beta | **-1.113 meV = -0.0257 kcal/mol** | -3.361 meV | +2.248 meV |
| `anti90` - beta (relaxed into the `anti_free` structure) | **-1.151 meV = -0.0265 kcal/mol** | -2.842 meV | +1.692 meV |

* **The antipolar packing is lower, by about 1.1 meV (0.026 kcal/mol) per monomer.** Polar
  beta and its antipolar competitor are degenerate on this Hamiltonian to a tenth of
  polyfind's 0.27 kcal/mol error bar and to a few hundredths of a kcal/mol absolutely. In the
  requirement's terms the static margin sits on the near-degenerate side, not the
  hard-ferroelectric side. That is a statement about two zero-kelvin lattice energies of
  perfect crystals on one Hamiltonian, not about loss; the caveats say what it does not cover.
* **Numerically converged** to better than 0.05 meV per monomer: k-mesh 4x2x8 to 6x3x12 to
  8x4x16 to 12x6x24 moves it by 0.004 meV; cutoff 90/360 to 110/440 to 130/520 Ry by at most
  0.04 meV; the two polar starts (polyfind's and the provider's) reach the same minimum to
  0.001 meV. The sign is not robust to the Hamiltonian, though: the margin is the difference
  of a +2.2 meV PBE part and a -3.4 meV dispersion part, so a 33% error in the D3 difference
  alone would flip it.
* **It is mostly dispersion.** At the relaxed geometries the D3 term favours the antipolar
  cell by 3.4 meV per monomer and the rest of PBE favours polar by 2.2 meV. (The split moves
  by 0.5 meV between the two antipolar cells while their totals agree to 0.04 meV, so only
  the total is a robust number; see "The D3 share".)
* **The gamma = 90 cell is not a separate minimum.** Started from polyfind's relaxed `anti90`
  (a = 4.49, b = 9.97 A, density 1.85), the DFT relaxation closes the long axis to 8.51 A and
  lands on the `anti_free` structure: same chain registry (chain 2 directly above chain 1 at
  half the long axis), same setting angles to 0.07 deg, cell within 0.2%, energy 0.04 meV
  per monomer lower through a 0.03 deg monoclinic shear its symmetry allows and
  `anti_free`'s does not.
* **Symmetry does not hold up the answer.** Relaxed again with no symmetry from a randomly
  perturbed geometry, `anti_free` stays antipolar in the same packing and ends 0.05 meV per
  monomer lower still. On the 4x2x8 mesh the lowest antipolar energy found is therefore
  -1.160 meV per monomer below beta, against -1.109 at the symmetric point: the margin is
  -1.11 to -1.16 meV (-0.026 to -0.027 kcal/mol) per monomer, depending only on which point
  of a very flat antipolar basin is quoted. Beta, put through the same symmetry-free
  relaxation, returns to itself (0.002 meV per monomer above its symmetric energy).
* **Polyfind agrees.** Its own all-atom-relaxed margin for the same pair is -0.0056 kcal/mol
  (-0.24 meV) per monomer, same sign, 0.02 kcal/mol from the DFT number; its gamma = 90 screen
  margin (+1.02 kcal/mol) is the artefact `docs/SCREEN_PHONONS.md` said it was.
* **Lattice.** Relaxed beta is a = 8.528, b = 4.757, c = 2.581 A against Hasegawa's 8.58,
  4.91, 2.56 (-0.6%, -3.1%, +0.8%): PBE-D3(BJ) shortens the polar axis by 3%.
* **The control comes out right.** On the same Hamiltonian, polyfind's antipolar TGTG'
  (alpha-like) cell sits 31.5 meV = 0.73 kcal/mol = 3.0 kJ/mol per monomer below beta. That
  is inside the 2.6-3.8 kJ/mol of published dispersion-corrected DFT and beside polyfind's own
  0.747 kcal/mol. This control ran on a reduced protocol: relaxed at a coarse mesh, energies
  at the beta cells' 4x2x8 density.

## Why this number

The governing requirement (`RESUME.md`) is low tan-delta and a high-frequency response. The
static quantity that bears most directly on it is the polar-versus-antipolar lattice energy
margin: polar packing well below antipolar is the hard-ferroelectric case, near-degenerate is
the relaxor case. Polyfind cannot resolve it for PVDF. Its error bar is 0.27 kcal/mol per
monomer (`docs/SCREEN.md`, "How to read any energy"), and on its own potential beta's polar
cell and its best antipolar cell are degenerate once gamma is free (screen packer, rigid
chains: -0.000 against the best cell, -0.017 against the polar one; phonon packer, every atom
relaxed at fixed cell: -0.006 kcal/mol per monomer;
`docs/SCREEN_PHONONS.md`). A first-principles margin for the same pair is the calibration
number the potential is missing.

## Hamiltonian: the provider's, checked against their files

The provider's Quantum ESPRESSO references for beta-PVDF
(`sarco/materials/gpu_bundle/periodic_reference/beta_pvdf/born/scf.in`, and the Berry inputs
generated by `berry_sweep.py`) use, and this calculation uses:

| setting | value |
|---|---|
| functional | `input_dft = 'PBE'` |
| dispersion | `vdw_corr = 'grimme-d3'`, `dftd3_version = 4` (Becke-Johnson damping), `dftd3_threebody = .false.` |
| pseudopotentials | SSSP 1.3.0 PBE precision: `H_ONCV_PBE-1.0.oncvpsp.upf`, `C.pbe-n-kjpaw_psl.1.0.0.UPF`, `F.oncvpsp.upf` (masses and species card copied from `born/scf.in`) |
| cutoffs | `ecutwfc = 90`, `ecutrho = 360` Ry |
| occupations | `fixed` (insulator) |

The pseudopotentials were downloaded from the Materials Cloud SSSP record
(doi:10.24435/materialscloud:f3-ym, `SSSP_1.3.0_PBE_precision.tar.gz`, archive md5
`fde94756886f32ada7bf597547557eb5` as listed by the record). Their md5 sums match the
`SSSP_1.3.0_PBE_precision.json` that ships with them, and their SHA-256 sums match,
digit for digit, the ones the provider recorded for their own copies in
`beta_pvdf/ELECTRICAL_DISPERSION_SCOPE.md`:

| file | md5 (= SSSP json) | SHA-256 (= provider's record) |
|---|---|---|
| `H_ONCV_PBE-1.0.oncvpsp.upf` | `1790becc920ee074925cf490c71280fe` | `69d63178...a5c09` |
| `C.pbe-n-kjpaw_psl.1.0.0.UPF` | `5d2aebdfa2cae82b50a7e79e9516da0f` | `9900d1ef...7f3d7` |
| `F.oncvpsp.upf` | `dbf8367619a749999a0442b33b80a7e7` | `ba6ca6af...7b8` |

Differences from the provider's runs, all stated:

* **QE version.** 7.5 (conda-forge build, OpenMPI 5, OpenBLAS); theirs is 7.6 built from
  source. Same D3 implementation family; not checked line by line.
* **Symmetry.** Theirs runs with `nosym`/`noinv` because Berry phases and DFPT need the full
  grid. Here QE finds and uses the symmetry of each cell (4 operations for every cell below,
  centrosymmetric for the antipolar ones); this changes k-point bookkeeping and constrains the
  relaxations to the symmetry the start already has. A symmetry-free relaxation from a
  perturbed geometry is the check on that constraint (below).
* **Three-body dispersion.** Their Born/DFPT SCF disables the Axilrod-Teller-Muto term (QE has
  no second derivative for it); their Berry SCFs enable it. The primary numbers here are
  two-body, as requested and as in `born/scf.in`; the three-body shift of the margin is
  reported separately at fixed geometry.
* **Bands.** Default `nbnd` (occupied bands only) instead of their 32; with fixed occupations
  this does not change the energy.

## Cells and starts

All four are two-chain, one-monomer-per-chain all-trans cells (12 atoms, 2 monomers).

| DFT cell | start | what it is |
|---|---|---|
| `polar` | polyfind `polar90` (= `polar_free` for PVDF), every atom relaxed at fixed cell on the phonon packer | beta-PVDF, chain 2 at (1/2, 1/2), both CF2 groups pointing the same way |
| `anti_free` | polyfind `anti_free`, relaxed likewise | the gamma-free antipolar competitor of `docs/SCREEN_PHONONS.md`: a = 4.56, b = 9.64, gamma = 61.8 deg, which is a rectangular 4.56 x 8.49 A cell with chain 2 directly above chain 1 along the long axis, antiparallel |
| `anti90` | polyfind `anti90`, relaxed likewise (the classical relaxation has already slid the chains 1.05 A) | the gamma = 90 exactly antipolar cell the screen compares, a saddle on the classical potential |
| `polar_dipole` | the provider's beta geometry, `born/scf.in` (their CP2K PBE-D3(BJ) zero-pressure cell) | a cross-check that a different start reaches the same DFT minimum |

The polyfind starts are regenerated by `examples/polar_margin_dft.py starts`, which calls
`examples/screen_phonons.py`'s own `one_chemistry` for PVDF with that script's defaults and
captures each cell's placement before and after `phonon.relax_all_atom` (accepted only when
its own criterion is met; all four were). It reproduces the record: all-atom-relaxed
anti_free minus polar -0.0056 kcal/mol per monomer, anti90 minus polar +0.597 at the
classical fixed cell. The in-plane basis is Lagrange-Gauss reduced before writing (Cartesian
positions unchanged), so the anti_free cell enters QE as the rectangular cell.

## Protocol

* **Relaxation.** `vc-relax`, `cell_dofree = 'all'`, BFGS for ions and cell, `press = 0`,
  `forc_conv_thr = 1e-4` Ry/bohr, `press_conv_thr = 0.5` kbar, `etot_conv_thr = 1e-6` Ry,
  `conv_thr = 1e-10` Ry. Two stages: a 2x1x4 pre-relaxation (mesh labels below are
  (short in-plane, long in-plane, chain axis)) and then 4x2x8 from its final geometry. Each
  vc-relax ends with QE's fresh-basis SCF at the final cell, and a stage is accepted only if
  that SCF has every force component below 1e-4 Ry/bohr and every stress component below
  0.5 kbar (otherwise it is repeated from its own final geometry; none needed it).
* **Final SCFs** at the relaxed geometry with identical settings and `conv_thr = 1e-11` Ry:
  meshes 2x1x4, 4x2x8, 6x3x12, 8x4x16 (the provider's 4x8x16 in their axis order) and 12x6x24;
  cutoffs 90/360, 110/440 and 130/520 Ry on the 4x2x8 mesh; and the three-body variant on the
  4x2x8 mesh. The 8x4x16 SCF is also the acceptance check of the relaxed geometry at the
  production mesh (if it failed the force/stress criteria the cell would be re-relaxed there;
  none needed it).
* **Same density for every cell.** The meshes are mapped onto each cell's own short and long
  in-plane axes. All cells hold two chains in nearly the same area, so the sampling density is
  matched without interpolation.
* Before the cells were run, the k-mesh was tested on the provider's own geometry (4x8x16,
  135 irreducible points, against 2x4x8 and 3x6x12): 3x6x12 reproduces 4x8x16 to 1.1e-7 Ry
  per cell, 2x4x8 to 1.9e-5 Ry, with stress within 0.01 kbar and forces within 3e-5 Ry/bohr.

## Results

### Relaxed cells

Every relaxation of the four beta-type cells and every 8x4x16 SCF met the criteria at its
final geometry (largest force component, largest stress component; alpha has its own
section):

| cell | a (short) | b (long) | c (chain) | gamma | V / monomer | QE symmetry | 8x4x16: force, stress | arrangement |
|---|---:|---:|---:|---:|---:|---|---|---|
| `polar` | 4.7568 | 8.5282 | 2.5808 | 90.00 | 52.346 | 4 ops, no inversion (+ C-centring) | 3.5e-5 Ry/bohr, 0.05 kbar | polar: both chains' CF2 the same way; chain 2 at (1/2, 1/2) |
| `polar_dipole` | 4.7594 | 8.5259 | 2.5808 | 90.00 | 52.361 | 4 ops, no inversion (+ C-centring) | 1.7e-5, 0.05 | polar, the same structure |
| `anti_free` | 4.7429 | 8.5190 | 2.5808 | 90.00 | 52.140 | 4 ops with inversion | 9.2e-5, 0.12 | antipolar: CF2 groups opposite; chain 2 at (0, 1/2) |
| `anti90` | 4.7523 | 8.5148 | 2.5807 | 89.97 | 52.213 | 4 ops with inversion | 8.5e-5, 0.15 | antipolar, the `anti_free` structure |

(A, deg, A^3.) The arrangement is read from the relaxed geometry, not assumed: per chain the
in-plane direction from the CF2 carbon to its two fluorines' midpoint (setting angle), and
chain 2's carbon-centroid offset from chain 1 in fractional coordinates
(`orientation` in the record). Polar cells: setting angles equal, offset (1/2, 1/2, 0).
Antipolar cells: setting angles 180 deg apart (to 0.07 deg in `anti90`), offset
(-0.004, 1/2, 0) and (-0.002, 1/2, 0). The polar cells stayed polar and the antipolar cells
stayed antipolar, which their symmetry guarantees once QE has found it (the antipolar cells
are centrosymmetric). Whether a symmetry-free relaxation keeps them there is the separate
check below.

**Against experiment** (beta, Hasegawa et al. 1972 via the provider's `hasegawa_planar.cif`):
a = 8.528 (8.58, -0.6%), b = 4.757 (4.91, -3.1%), c = 2.581 (2.56, +0.8%), volume -2.9%.
The polar axis is where PBE-D3(BJ) is short. For comparison, the provider's CP2K PBE-D3(BJ)
cell (`born/scf.in`) is 8.358 x 4.731 x 2.580, 2.0% shorter on the long axis than the QE
minimum of the same functional. The two codes' minima differ (Gaussian-basis CP2K against
converged plane waves); evaluated in QE, the provider's geometry carries forces up to
0.0156 Ry/bohr (0.40 eV/A, on fluorine) and a pressure of -13.4 kbar (stress -23.9 kbar along
the long axis, -15.4 along the polar axis), and its energy lies 8.7 meV per monomer above
the QE-relaxed beta (`kmesh_test_provider_geometry.json`). Their QE Berry and DFPT numbers are therefore
evaluated at a geometry that is not a minimum of the Hamiltonian they are evaluated with.
That is recorded here as a fact about the inputs; what it does to their response numbers was
not computed.

**Polyfind's starts** for comparison (phonon packer, all-atom relaxed at fixed cell):
polar 4.549 x 8.500 x 2.563, `anti_free` 4.559 x 8.488 x 2.563, `anti90` 4.487 x 9.970 x 2.563.
The DFT relaxation expands the short (polar) axis of both main cells by 4-5%. The chain registry
survives it: polyfind's `anti_free` has chain 2 at (0.06, 1/2), DFT's at (0.00, 1/2).

### The margin and its convergence

E(antipolar) - E(polar) per monomer, every entry at the same relaxed geometries
(two-body D3, 90/360 Ry, unless stated). Mesh labels are (short, long, chain axis) divisions.

| | `anti_free` - polar (meV) | (kcal/mol) | `anti90` - polar (meV) |
|---|---:|---:|---:|
| k-mesh 2x1x4 | +1.712 | +0.0395 | |
| k-mesh 4x2x8 | -1.109 | -0.0256 | -1.146 |
| k-mesh 6x3x12 | -1.113 | -0.0257 | |
| **k-mesh 8x4x16 (the provider's 4x8x16)** | **-1.113** | **-0.0257** | **-1.151** |
| k-mesh 12x6x24 | -1.113 | -0.0257 | |
| 4x2x8, 110/440 Ry | -1.148 | -0.0265 | -1.178 |
| 4x2x8, 130/520 Ry | -1.114 | -0.0257 | |
| 4x2x8, D3 three-body on | -1.213 | -0.0280 | -1.311 |

* **k-mesh.** From 4x2x8 on, the margin moves by 0.004 meV; 6x3x12, 8x4x16 and 12x6x24
  agree to 0.0002 meV (-1.1129, -1.1127, -1.1128). The total energies themselves agree to
  2e-8 Ry per cell between 8x4x16 and 12x6x24. The 2x1x4 row, the pre-relaxation mesh, has the wrong sign and is
  shown only as a warning. At this margin's scale the mesh matters, and a coarse-mesh
  answer would have said "polar by 1.7 meV".
* **Cutoff.** 90/360 to 110/440 to 130/520 Ry: -1.109, -1.148, -1.114. Not monotonic, a
  spread of 0.04 meV. The absolute energy per monomer moves by 2.1 meV over the same range,
  so the difference cancels to 2%.
* **Geometry tolerance.** Residual stresses of at most 0.15 kbar and forces of at most
  9.2e-5 Ry/bohr bound the energy error of an imperfect relaxation far below 0.01 meV
  (second order in the residual). The two polar starts, relaxed independently, agree in
  energy to 0.0007 meV per monomer at 8x4x16 while their cells differ by 0.003 A: the surface
  is that flat.
* **Converged to better than 0.05 meV (0.001 kcal/mol) per monomer** in k-mesh, cutoff and
  relaxation, against the requested 0.4 meV.
* **Three-body dispersion.** Switching on the Axilrod-Teller-Muto term at fixed geometry adds
  +19.11 meV per monomer to the polar cell's energy and +19.01 to `anti_free`'s (it is
  repulsive), which moves the margin by -0.10 meV (to -1.21); for `anti90` by -0.17 meV (to
  -1.31). It is not a relaxation with the three-body term on.

### The D3 share

QE prints the D3 energy separately; it depends on the nuclear positions only. At the relaxed
geometries, per monomer:

| pair | total | D3 part | everything else (PBE) |
|---|---:|---:|---:|
| `anti_free` - polar | -1.113 | -3.361 | +2.248 |
| `anti90` - polar | -1.151 | -2.842 | +1.692 |

PBE alone prefers the polar packing; dispersion prefers the antipolar one by more. The
split is not a robust number. The two antipolar cells are the same structure up to a 0.2%
strain and a 0.03 deg shear, and their totals agree to 0.04 meV, but their D3 parts differ by
0.52 meV. Along a flat direction of a relaxed cell the two parts trade off (that is what zero
stress means), so only their sum is pinned. What is robust is the order of magnitude: the
dispersion difference is three times the margin, so the margin's sign rests on D3(BJ) being
right to about 30% for the difference between two packings of the same chain. Nothing here
tests that.

### What the gamma = 90 antipolar cell became

On the classical potential it is a saddle. Let go at fixed cell, its chains slide 1.05 A in
opposite directions (`docs/SCREEN_PHONONS.md`), and that slid geometry was the DFT start. QE
finds four symmetry operations with inversion in it, so the relaxation is constrained to
keep it centrosymmetric (antipolar). With the cell free, the long axis closes from 9.97 to
8.51 A, the density rises from 1.85 to 2.04 g/cm^3, and the structure ends on `anti_free`:

| | `anti_free` | `anti90`, relaxed |
|---|---|---|
| cell (A, deg) | 4.7429 x 8.5190 x 2.5808, 90.00 | 4.7523 x 8.5148 x 2.5807, 89.97 |
| chain 2 offset (fractional) | (-0.0036, 0.5000, 0) | (-0.0022, 0.4998, 0) |
| setting angles (deg) | 0.00, 180.00 | -0.07, 179.93 |
| E, 8x4x16 (Ry per cell) | -275.80719396 | -275.80719954 |

It is the same packing. It is 0.038 meV per monomer lower through a small monoclinic shear
and rotation that its symmetry operations allow and `anti_free`'s exclude. That is a tenth
of the requested resolution. Whether `anti_free`'s orthorhombic point is a shallow saddle
along that direction is what the symmetry-free check below tests: it is, and the drop is
0.05 meV per monomer.

### Symmetry-free check

Each main cell's final geometry was perturbed (every atom displaced by a seeded random
0.03 A rms, the cell by a random symmetric 0.5% strain) and relaxed again with `nosym` on the
4x2x8 mesh under the same criteria (both passed: largest force 9.1e-5 Ry/bohr and stress
0.09 kbar for beta, 6.7e-5 and 0.06 for `anti_free`). The comparison is against the symmetric
relaxation's fresh-basis energy on the same mesh.

| cell | symmetric (Ry per cell) | symmetry-free (Ry per cell) | difference per monomer | where it went |
|---|---:|---:|---:|---|
| `polar` | -275.80701219 | -275.80701185 | +0.002 meV | back to beta: chain 2 at (0.499, 0.500, 0.003), setting angles equal to 0.07 deg, cell angles within 0.08 deg of 90 |
| `anti_free` | -275.80717525 | -275.80718272 | -0.051 meV | the same antipolar packing: chain 2 at (0.000, 0.500, 0.002), setting angles 180 deg apart to 0.06 deg, cell 4.7469 x 8.5166 x 2.5807 with angles within 0.03 deg of 90 |

Beta is a minimum, not a symmetric saddle: with no symmetry it comes back to itself, to the
relaxation tolerance. `anti_free`'s symmetric point is not quite the bottom of its basin. The
symmetry-free relaxation finds the same packing 0.05 meV per monomer lower, as `anti90` did
through its monoclinic shear (0.04 meV). The basin is that flat in those directions. Its
lowest point found puts the margin at **-1.160 meV (-0.0268 kcal/mol) per monomer** on the
4x2x8 mesh, against -1.109 at the symmetric point. No symmetry was needed to keep the polar
cell polar or the antipolar cell antipolar. (`anti90`'s own symmetry-free check was dropped
from the queue: it had relaxed into the `anti_free` structure, whose check covers it.)

### Polyfind's margins beside these

kcal/mol per monomer (meV in brackets), polyfind from `deliverables/polar_margin_dft/starts/polyfind_starts.json`
(regenerated here, reproducing `deliverables/screen_phonons/`):

| quantity | `anti_free` - polar | `anti90` - polar |
|---|---:|---:|
| polyfind screen packer, rigid chains (the screen's column) | -0.017 (-0.74) (gamma free; the record's "-0.000" is against the best cell, which is this antipolar one) | +1.021 (+44.3) (gamma = 90) |
| polyfind phonon packer, rigid chains, cell re-polished | -0.0076 (-0.33) | +1.166 (+50.6) |
| polyfind phonon packer, every atom relaxed, cell fixed | **-0.0056 (-0.24)** | +0.597 (+25.9) |
| **DFT, PBE-D3(BJ), cell and atoms relaxed** | **-0.0257 (-1.11)** | **-0.0265 (-1.15)** |
| polyfind's error bar | 0.27 (11.8) | |

For the pair that matters, polyfind and DFT agree in sign and to 0.02 kcal/mol per monomer.
The classical `anti90` numbers are at its fixed gamma = 90 cell (a = 4.49, b = 9.97 A). DFT
relaxed that cell, polyfind's protocol did not, so that column compares different cells. The
DFT result is what the classical gamma-free search (`antipolar_gamma_free`, every seed off
90 deg) already found: the gamma = 90 cell is not the antipolar minimum.

### The alpha/beta control

A check that this Hamiltonian, at these settings, gets a known polymorph gap right. It ran
on a **reduced protocol**, stated in full. The cell is polyfind's lowest antipolar TGTG'
two-chain cell (24 atoms; `pack(gamma_free=True)` under `pvdf-dft-valence` with Ewald,
`examples/polar_margin_dft.py starts --alpha`). It was relaxed with the same vc-relax settings
at the 2x1x2 mesh only (60 BFGS steps, mesh halved along the chain because the repeat is two
monomers). The 4x2x4 re-relaxation would have taken 4-8 h on the three free ranks and was
stopped by design. The energies are 4x2x4 SCFs (the beta cells' 4x2x8 density) on that
geometry, compared with beta's 4x2x8 energy. At its own mesh the relaxed geometry's
fresh-basis SCF has a largest force of 1.5e-4 Ry/bohr, above the 1e-4 the beta cells met, and
0.13 kbar stress. At 4x2x4 the same geometry carries 9.4e-4 Ry/bohr and 0.96 kbar, so it is
not a 4x2x4 minimum. An estimate of the error this leaves, not a measurement: on the beta-type
cells the coarse-mesh geometry sat 0.17-0.41 meV per monomer above the 4x2x8-relaxed one when
both were evaluated at 4x2x8, about 1% of the gap below.

QE finds 2 symmetry operations (no inversion), so this is not the refined P2_1/c alpha
structure but polyfind's antipolar TGTG' packing. It relaxes to a = 4.945, b = 9.636,
c = 4.660 A against alpha's 4.96, 9.64, 4.62 (-0.3%, -0.0%, +0.9%), with a 2.3 deg monoclinic
angle (beta = 87.7 deg) that the experimental cell does not have. The two chains' mean CF2
directions are antiparallel.

| | per monomer |
|---|---:|
| **E(alpha) - E(beta), PBE-D3(BJ) two-body** | **-31.46 meV = -0.726 kcal/mol = -3.04 kJ/mol** |
| D3 part / PBE part | +20.09 / -51.55 meV |
| with the three-body term, fixed geometry | -34.01 meV (-3.28 kJ/mol) |
| polyfind, `pvdf-dft-valence` screen (`docs/SCREEN.md`) | -0.747 kcal/mol (-3.13 kJ/mol) |
| published, beta above alpha (`docs/REFERENCES.md` section 5) | DFT-D2 3.5, vdW-DF 3.8, vdW-DF2 2.6, PBE0+dispersion (Itoh) 2.9, plain PBE 4.6-6.5 kJ/mol |

Alpha is lower than beta by 3.0 kJ/mol per monomer, inside the dispersion-corrected
literature range. Dispersion favours the denser beta here (+20 meV), opposite to its role in
the polar/antipolar margin. Without the D3 term the gap is -5.0 kJ/mol at these geometries,
in the plain-GGA range. That comparison is loose, since those studies relaxed with PBE alone.
The control behaves: the Hamiltonian that puts the two all-trans packings within 1.2 meV of
each other puts alpha 31 meV below beta, where the literature puts it.

## Caveats, in full

* **Zero kelvin, static, one Hamiltonian.** The number is a difference of two relaxed lattice
  energies. It contains no zero-point energy and no thermal free energy. A rough estimate,
  not computed: the classical potential puts `anti_free`'s three softest lattice modes 7 to
  12 cm^-1 above beta's (`docs/SCREEN_PHONONS.md`). Summed over those three modes the shift is
  27 cm^-1, worth about 1.7 meV per cell (0.8 meV per monomer) of zero-point energy alone in
  the polar cell's favour: the same order as the margin and of the opposite sign. The softer
  polar lattice also has the larger vibrational entropy, which favours it further with
  temperature. A vibrational free-energy difference could therefore decide the order
  at room temperature, and it was not computed.
* **Not a barrier.** Two relaxed minima say nothing about the path between them, how fast a
  polar domain could convert, or whether it ever would without a field. In particular the
  margin says nothing about loss by itself. Near-degeneracy is the static precondition for
  the relaxor-like case the requirement worries about, not a demonstration of it.
* **Only the cells that were asked for.** Two-chain, one-monomer-repeat, all-trans cells. No
  larger antipolar supercells, other stackings, kinks, defects or chain-end effects were
  searched. A lower antipolar packing would only make the margin more negative. A lower polar
  one is not ruled out either, though beta is the established polar structure.
* **The Hamiltonian's error is not in the error bar.** "Converged" above means converged in
  basis and sampling for PBE-D3(BJ). The margin is a third of its own dispersion component,
  and how far PBE-D3(BJ) is from the exact answer for a 1 meV packing difference was not
  measured. Published functionals spread the alpha/beta gap of the same polymer over a factor
  of 2.5 (`docs/REFERENCES.md` section 5). Read the sign as "this Hamiltonian's", and the
  magnitude, "within a few meV per monomer", as the robust statement.
* **Two-body D3** as requested and as in the provider's Born SCF. The three-body term is
  reported at fixed geometry only (-0.10 and -0.17 meV); the geometry was not relaxed with it.
* **Symmetry.** Relaxations keep the symmetry QE found in the start. The symmetry-free check
  is one perturbed relaxation per cell at 4x2x8, not a search. It found the antipolar
  cell's symmetric point 0.05 meV above a lower-symmetry variant of the same packing. A
  broader search could find more such variants, presumably of the same order, but none was run.
* **QE 7.5, not the provider's 7.6.** Same inputs, same pseudopotential bytes; the version
  difference was not tested.
* **The provider's geometry is not a QE minimum** (forces 0.40 eV/A, pressure -13.4 kbar in
  this Hamiltonian). This bears on their response numbers, not on this margin; it was not
  chased further.
* **Experiment.** Beta-PVDF is observed as a polar crystal in drawn and poled films. A static
  zero-kelvin preference of 1 meV per monomer for an antipolar packing, on one Hamiltonian
  and inside the vibrational and Hamiltonian uncertainties above, is not evidence against
  that. It is also not an explanation of it.

## Reproducing this

```
# Windows side, polyfind worktree (two minutes): the starts
set PYTHONPATH=src
python examples/polar_margin_dft.py starts [--alpha]

# WSL side: QE 7.5 from conda-forge, the SSSP 1.3.0 PBE precision files in ~/polar_margin/pseudo
micromamba create -n qe -c conda-forge qe=7.5 openmpi
NP=3 NK=3 deliverables/polar_margin_dft/runner/drive_cell.sh polar RELAX FINALCHECK scf_k0 scf_k1 ...
python examples/polar_margin_dft.py collect      # -> deliverables/polar_margin_dft/raw_stages.json, inputs/, relaxed/
python examples/polar_margin_dft.py report       # -> deliverables/polar_margin_dft/polar_margin_dft.json
```
Every pw.x input actually run is in `deliverables/polar_margin_dft/inputs/<cell>/<stage>.in`.
The runner scripts (`run_pw.sh`, `drive_cell.sh`, `worker.sh`, `alpha_reduced.sh`, and the
job queue in the order it ran) are in `deliverables/polar_margin_dft/runner/`. Paths in them
are this machine's. The pw.x outputs themselves are not committed; they remain on the WSL side
under `~/polar_margin/runs`, and every number above is parsed from them into
`raw_stages.json`.
