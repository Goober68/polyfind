# Electromechanical response: what `polyfind.mechanics` can and cannot compute

`pack.py` says which polymorph is stable and how polar it is. This document is about the
other half of the measurement need: how much mechanical work the crystal does in a field.
`src/polyfind/mechanics.py` computes the elastic stiffness, the piezoelectric coefficients
in both directions, and the three numbers an actuator is judged by — blocking stress, free
strain, work density. `examples/electromechanics.py` produces every table below:
`--preset illustrative|fitted|valence` for section 5, `--axial` for section 4.1,
`--stiffness-sweep` for section 4.3, `--cutoff` for section 6.

**Two paths, and the difference between them is one opt-in flag.**

* The **rigid path** is what this module did before, and it is still the default and the
  only thing available under the illustrative potential. The chain is a rigid body, three
  of the six Voigt strains are reachable, `C_33` is `nan`, and every diagonal column of the
  piezoelectric tensor is exactly zero.
* The **deformable path** forwards the fitted valence terms of `docs/VALENCE_FIT.md` into
  the lattice kernel (`CrystalPacker(valence=...)`, opt-in), which gives the chain a
  restoring force it never had. Four of six strains are then reachable, `C_33` is a
  constant of the potential rather than of an invented number, and the diagonal columns can
  be non-zero — though for beta-PVDF, for a reason that is itself a result, they are not.

Sections 1 to 3 set up both. Section 4 is the axial constant and the evidence that it has
stopped being an artifact. Section 5 has the tables. Sections 6 to 9 are convergence, cost,
literature and what is still missing.

## 1. Which strains this cell can express

A cell in `pack.py` is

    a_vec = (a, 0, 0)   b_vec = (b cos gamma, b sin gamma, 0)   c_vec = (0, 0, c)

`c_vec` is pinned to z, and z is also the chain axis. A homogeneous strain sends the rows
to `r_i (1 + eps)`, so the third row becomes `c (eps_xz, eps_yz, 1 + eps_zz)`.

| Voigt | component | rigid path | deformable path | why |
|---|---|---|---|---|
| 1, 2, 6 | `eps_xx`, `eps_yy`, `eps_xy` | **reachable** | **reachable** | carried by `(a, b, gamma)`. The rigid rotation about z that puts `a_vec` back on x is absorbed exactly by the setting angles `phi1, phi2`. |
| 4, 5 | `eps_yz`, `eps_xz` | **not expressible** | **not expressible** | these tilt `c_vec` off z. There is no cell variable for a third lattice vector that is not along the chain axis, so `C_44`, `C_55`, `C_45` and every mixed constant containing a 4 or a 5 are absent — not unconverged, absent. `strained_cell` raises rather than silently dropping them. |
| 3 | `eps_zz` | **not computable** | **reachable** | `c` is not a cell variable: it is the rise of the chain's repeat transform, so it moves only if the chain's *shape* does. Section 4. |

So the reachable stiffness is the 3×3 block `{11, 22, 12, 16, 26, 66}` rigid, and the 4×4
block that adds `{33, 13, 23, 36}` deformable. Those are the *true* constants and not a
plane-stress reduction — `C_11 = dsigma_1/deps_1` is defined at fixed `eps_2 … eps_6`, and
fixed `eps_zz` is what a rigid chain enforces for free and what an explicit constraint
enforces when it does not. On the rigid path the compliance `S` is the inverse of the 3×3
block only, so every `d` is an **axially clamped** coefficient; on the deformable path `S`
is the inverse of the 4×4 and `d` is clamped only in the two shears the cell cannot express
at all.

"Relaxed-ion" here means relaxed *chains*. The rigid-body internal coordinates — each
chain's setting angle and the axial offset `dz` — are relaxed at every strain state by
`pack.polish` on the analytic cell gradient. On the deformable path the chain's own
line-group shape parameters (`polyfind.linegroup`, the same parametrisation
`refine_crystal` uses) are relaxed with them, subject to an equality constraint that fixes
the repeat. That constraint is the whole trick: it is what lets the chain move and still
lets the calculation say which strain it is at.

## 2. Conventions

* Engineering Voigt strain `eps = (e_xx, e_yy, e_zz, 2 e_yz, 2 e_xz, 2 e_xy)`, tension positive.
* `sigma_J = (1/V0) dE/de_J` in GPa, tension positive, `V0` the reference cell volume.
* The electric enthalpy per unit reference volume is `h(eps, E) = u(eps) - E . P`, which is
  what the kernel minimises: it adds `-mu . E` to the cell energy. Then `sigma_J = dh/de_J`,
  `P_i = -dh/dE_i`, and

      e_iJ = (1/V0) dmu_i/de_J |_E  =  -dsigma_J/dE_i |_eps   (C/m^2)
      d_iJ = de_J/dE_i |_sigma=0    =  (e S)_iJ, S = C^-1     (pC/N = pm/V)

* `e` is the **proper** piezoelectric stress constant: the dipole per *reference* volume.
  The naive `dP/deps` with `P = mu/V` differs from it by `P_i` on the three diagonal
  columns. `Piezoelectric` reports both, and the difference is asserted to be `P_i` in
  `tests/test_mechanics.py` on both paths.
* **The axial stress is a Lagrange multiplier.** On the deformable path the relaxed energy
  at strain `eps_zz` is a minimum subject to `c(shape) = c0 (1 + eps_zz)`, and the
  multiplier of that constraint is `dE*/dc`, so `sigma_3 = (c0/V0) lambda`. It falls out of
  the stationarity condition `grad_shape E = lambda grad_shape c` at no extra cost, and the
  part of the shape gradient that is *not* along the constraint is reported as
  `Relaxed.residual` — a direct measurement of how converged each constrained relaxation is
  (it runs at 1e-16 for the states in section 5).
* Actuator at field `E`: free strain `eps_free = d^T E` at zero stress in the reachable
  components, blocking stress `sigma_block = C eps_free = e^T E`, the stress the crystal
  exerts on a rigid clamp (the external stress needed to *hold* `eps = 0` is its negative),
  work density `w = (1/4) sigma_block . eps_free` against a matched linear load. The field
  used for the tables is **0.01 V/A = 100 MV/m**, a realistic drive field for poled PVDF and
  about twice its coercive field; the figures are the linear response evaluated there, not a
  measurement at that field.
* Poling direction: the tables use the direction that maximises the work, the leading
  eigenvector of `d C d^T`.

One trap, recorded because it cost two wrong answers before it was found. A shear leaves
`a_vec` off the x axis, so the packer's canonical frame is rotated relative to the reference
by `theta = -e_xy/2`. The dipole has to be rotated back out of that frame before it is
differentiated, *and* an applied field has to be rotated into it. Each of those terms is
`P/2` — the same size as the coefficient being measured — so getting either wrong does not
degrade the answer, it replaces it. The direct/converse agreement in section 3 is what
caught both.

## 3. The two piezoelectric routes, and the check

`e` comes from the dipole of strained, internally relaxed cells. `d` comes from solving
`sigma(eps, E) = 0` — Newton on the analytic stress, the internal coordinates (and, on the
deformable path, the chain shape) relaxed at every iterate. The identity `d = e S` follows
from `d^2 h / deps dE` being symmetric, and the two sides are computed by disjoint
machinery, so agreement is evidence.

| | illustrative | `pvdf-dft-fit` | `pvdf-dft-valence` (deformable) |
|---|---|---|---|
| PVDF beta | 0.65% | 0.50% | 0.03% |
| PVDF alpha | 0.43% | 0.35% | 0.34% |
| PVDF gamma | 1.03% | 1.54% | 0.29% |
| PE | both routes zero, difference < 1e-11 pC/N | both zero, < 1e-11 pC/N | both zero, < 2e-12 pC/N |

(Largest `|d_from_e - d_direct|` as a fraction of the largest coefficient. The residual is
finite-difference error in `e` at `deps = 2e-3` and Newton tolerance in `d`. On the
deformable path the Newton stops at 20x the reference's own residual stress rather than at a
fixed 1e-8 GPa, because below the relaxation's noise floor it is chasing nothing: with the
fixed tolerance it ran to the iteration cap on every field axis, took 15 times as long, and
came back with a *worse* answer.)

## 4. Axial strain, and how it stopped being an artifact

### 4.1 Why the rigid path cannot have it

`c` *can* be moved. Bond lengths are rigid, but the backbone bond angles are not, and
shifting every backbone angle of the repeat together moves the repeat over a range of
several per cent in both directions:

| | reachable `eps_zz` (illustrative) | residual `sigma_zz` at fixed chain shape |
|---|---|---|
| PVDF beta | −6.5% … +2.4% | −4.51 GPa |
| PVDF alpha | −11.8% … +2.1% | −4.75 GPa |
| PVDF gamma | −10.3% … +3.5% | −3.59 GPa |
| PE | −5.8% … +3.6% | −2.27 GPa |

So the obstacle is not geometry. It is that **nothing in the rigid potential resists the
deformation.** Three separate facts said so, and they are still true of the rigid path.

**(a) The rigid kernel has no valence terms.** `pack.py` is Lennard-Jones + damped-shifted-force
Coulomb + a Fourier torsion term. `FFParameters.applied` forwards the fitted LJ, charges
and torsions into packing and does *not* forward `bond_terms` or `angle_terms`, on the
documented grounds that packing holds the chain rigid so they are an additive constant. The
only thing then resisting an opened angle is `refine_crystal`'s own harmonic restraint, at
`angle_stiffness = 105` kcal/mol/rad², described in its docstring as "UFF-like for sp3
carbon" and excluded from the reported lattice energy.

**(b) `C_33` was mostly that constant.** Sweeping it (illustrative potential, GPa):

| k (kcal/mol/rad²) | 0 | 52.5 | 105 | 210 | share of `C_33(105)` from k |
|---|---|---|---|---|---|
| PVDF beta | 83.8 | 242.4 | 401.1 | 718.4 | 79% |
| PVDF alpha | 203.4 | 266.1 | 328.8 | 454.3 | 38% |
| PVDF gamma | 792.6 | 851.9 | 911.1 | 1029.6 | 13% |
| PE | 60.5 | 220.5 | 380.5 | 700.4 | 84% |

**(c) The reference was axially stress-free *only because of* that restraint.** Along the
same path, `sigma_zz` (GPa):

| k | 0 | 52.5 | 105 | 210 |
|---|---|---|---|---|
| PVDF beta | −5.30 | −2.64 | **+0.03** | +5.35 |
| PE | −2.37 | −1.03 | **+0.30** | +2.97 |

The refinement's stationary point in `c` sat at `k = 105` (beta to 0.03 GPa, PE to 0.3 GPa)
and nowhere else: at `k = 0` the crystal wanted to contract along the chain at several GPa
and keep going, the lattice energy alone having no axial minimum at all. `axial_report`
still measures all of this and still returns `computable=False`, because for a packer
without valence terms it is still the right answer.

### 4.2 The fix: valence terms in the lattice kernel

`SimpleFF` has fitted harmonic stretch and bend terms — 7 bond types and 7 angle types,
fitted against DFT energies and Cartesian forces over 31 chemistries, with the provenance
and the held-out scores in `docs/VALENCE_FIT.md`. `CrystalPacker(valence=ff)` evaluates them
on the chain **as a periodic object**: `ChainValence` carries the image index of every atom
of every term, so a backbone bond joining the last atom of the repeat to the first atom of
the next one is a term whose energy depends on `c`. That is where the axial restoring force
comes from, and `energy_and_grad` returns its derivative in `g_coords` and `g_c`, which is
all the plumbing anything downstream needs.

Two things stay outside. **The tabulated chain-pair interaction stays rigid**: it is
tabulated for one set of chain coordinates and one energy-shifted pair potential, so
`_table_starts` refuses a packer carrying valence terms or a force-shifted cutoff, and
refuses a table whose chain is not the chain being packed. The table serves the *screen*;
deformation belongs to the direct kernel, which is what `refine_crystal` and this module
use. And **bond stretching is inert**: `build_chain` places every atom at the polymer's own
bond length, so `r - r0` is fixed by the chemistry and the stretch block contributes a
constant and a zero gradient. It is evaluated anyway, so that the reported energy is the
potential's, but the stretch channel — roughly half of a real chain modulus — is absent.

### 4.3 The evidence that the contamination is gone

The same sweep, run on the deformable path. `angle_stiffness` still moves the structure the
refinement hands over — the refined angles move by 0.85° for beta, 2.7° for PE, and the
refined `c` by up to 1.7% — and the relaxation against the *fitted* valence terms then has
to put it back:

| k (kcal/mol/rad²) | 0 | 52.5 | 105 | 210 | spread | refined `c` moved by | relaxed `c` |
|---|---|---|---|---|---|---|---|
| PVDF beta | 317.5990 | 317.5990 | 317.5990 | 317.5990 | 6e−6 GPa | 2.5632 → 2.5758 (+0.49%) | 2.56317 always |
| PVDF alpha | 250.1032 | 250.1032 | 250.1032 | 250.1032 | 2e−5 GPa | 4.6789 → 4.6624 (−0.35%) | 4.67892 always |
| PVDF gamma | 92.3652 | 92.3652 | 92.3651 | 92.3652 | 7e−5 GPa | 9.1174 → 9.1606 (+0.47%) | 9.07362 always |
| PE | 286.7916 | 286.7917 | 286.7916 | 286.7916 | 4e−5 GPa | 2.4934 → 2.5350 (+1.67%) | 2.49233 always |

`C_33` in GPa. The spread over the sweep is 1e−6 to 8e−5 *per cent* of the value — the
relaxation's own reproducibility, not a physical dependence — against a factor of 8.6 on the
rigid path. **The invented constant is no longer in the answer.** It cannot be: it enters
only through the starting structure, and the relaxation over `c` has a stationary point of
its own to find.

Three further checks that this number is what it says it is.

* **The reference is axially stress-free by measurement, not by construction.** The
  residual `sigma_zz` at the relaxed reference is below 1e-5 GPa for all four crystals,
  against several GPa on the rigid path (section 4.1). The multiplier route is what reports
  it, so this is also a check on the multiplier.
* **`C_33` is measured twice by different routes** — the derivative of the constraint
  multiplier, and the curvature of the relaxed energy itself — and they agree to 0.02-0.04%
  (beta 317.599 / 317.718, alpha 250.103 / 250.062, gamma 92.365 / 92.336, PE 286.792 /
  286.675).
* **The axial row and column of the 4×4 are also measured twice**: `C_J3` from the analytic
  cell gradient over the axial column, `C_3J` from the multiplier over the in-plane columns.
  The largest disagreement before symmetrising is 0.004 GPa for beta, 0.016 for alpha, 0.043
  for gamma and 0.199 for PE. PE's is 0.07% of its `C_33` but 3% of its `C_36 = 5.9`, which
  is the one entry of the set that should be read as two significant figures at most.

**What it is still not.** `C_33` here is an **upper bound**, on three counts, all of them
constraints that can only stiffen a path: bond lengths are rigid, so the stretch channel is
missing entirely; the line-group parametrisation admits only conformations of the repeat's
own symmetry; and beta's and PE's chains have exactly **one** free shape parameter after the
closure condition, so their axial deformation is a single-degree-of-freedom path. It is not
contaminated by `angle_stiffness`, by the refinement's starting point, or by the choice
between the two independent routes. It is limited by what the chain model can do.

## 5. Results

### 5.1 Elastic constants (GPa, relaxed-ion)

Illustrative potential, **rigid chain, `eps_zz` clamped**:

| | a, b, gamma, c (A, deg) | `C_11` | `C_22` | `C_12` | `C_66` | `C_16` | `C_26` |
|---|---|---|---|---|---|---|---|
| PVDF beta | 4.597, 8.587, 90.00, 2.630 | 44.42 | 37.76 | 2.25 | 3.10 | 0.00 | 0.00 |
| PVDF alpha | 5.096, 8.865, 90.00, 4.829 | 28.66 | 22.51 | 7.31 | 12.03 | 0.00 | 0.00 |
| PVDF gamma | 5.216, 8.892, 91.59, 9.566 | 29.56 | 31.58 | 9.15 | 9.25 | 2.19 | −0.48 |
| PE | 4.750, 7.080, 92.73, 2.575 | 27.46 | 30.46 | 16.14 | 15.83 | 4.57 | −4.29 |

`pvdf-dft-fit`, **rigid chain, `eps_zz` clamped**:

| | a, b, gamma, c (A, deg) | `C_11` | `C_22` | `C_12` | `C_66` | `C_16` | `C_26` |
|---|---|---|---|---|---|---|---|
| PVDF beta | 4.544, 8.621, 90.00, 2.618 | 21.60 | 15.47 | 1.84 | 2.90 | 0.00 | 0.00 |
| PVDF alpha | 5.062, 9.006, 90.00, 4.740 | 13.87 | 11.10 | 4.04 | 5.24 | 0.00 | 0.00 |
| PVDF gamma | 5.322, 8.869, 93.84, 9.410 | 11.56 | 12.72 | 3.44 | 4.30 | 0.47 | 0.37 |
| PE | 3.917, 6.788, 90.06, 2.534 | 36.22 | 32.43 | 16.69 | 24.04 | 14.76 | −11.48 |

`pvdf-dft-valence`, **deformable chain**, the full reachable 4×4:

| | a, b, gamma, c (A, deg) | `C_11` | `C_22` | `C_33` | `C_12` | `C_13` | `C_23` | `C_66` | `C_16` | `C_26` | `C_36` |
|---|---|---|---|---|---|---|---|---|---|---|---|
| PVDF beta | 4.603, 8.550, 90.00, 2.5632 | 22.46 | 17.59 | **317.6** | 2.17 | 12.57 | 2.28 | 3.36 | 0.00 | 0.00 | 0.00 |
| PVDF alpha | 5.072, 8.923, 90.00, 4.6789 | 14.87 | 11.48 | **250.1** | 5.17 | 0.63 | 7.73 | 5.45 | 0.00 | 0.00 | 0.00 |
| PVDF gamma | 5.447, 8.887, 95.96, 9.0736 | 12.46 | 12.36 | **92.4** | 3.61 | 10.77 | 4.45 | 4.14 | 0.33 | −0.43 | −4.38 |
| PE | 4.343, 6.766, 86.95, 2.4923 | 17.99 | 16.96 | **286.8** | 9.97 | 9.05 | 10.86 | 11.01 | −4.92 | 3.87 | 5.93 |

Beta and alpha come out orthorhombic in these axes (`C_16 = C_26 = 0` to rounding), a
symmetry the calculation was not told about and recovers. Gamma and PE relax to
`gamma != 90`, so they have genuine monoclinic `C_16`, `C_26`, `C_36`.

**What the deformation itself changes**, same potential and same relaxed reference, rigid
against deformable:

| | `C_11` rigid → deformable | `C_22` | `C_12` | largest \|e\| (C/m²) |
|---|---|---|---|---|
| PVDF beta | 22.464 → 22.464 | 17.594 → 17.594 | 2.166 → 2.166 | 0.0145 → 0.0145 |
| PVDF alpha | 16.327 → 14.867 (−9%) | 12.176 → 11.479 (−6%) | 4.223 → 5.171 | 0.0065 → 0.0992 (×15) |
| PVDF gamma | 13.873 → 12.455 (−10%) | 13.794 → 12.360 (−10%) | 4.429 → 3.613 | 0.0156 → 0.1900 (×12) |
| PE | 17.990 → 17.990 | 16.956 → 16.956 | 9.967 → 9.968 | 0 → 0 |

Beta and PE are *identical* to five figures, and that is not a coincidence or a
convergence accident — it is the mechanism. Their line group leaves exactly one free shape
parameter, and the constraint that holds `eps_zz = 0` fixes it, so under an in-plane strain
their chains have nowhere to relax to. Alpha and gamma have three and five parameters, so
freedom is left over, and letting them use it softens the in-plane block by about a tenth
and multiplies the piezoelectric response by an order of magnitude.

### 5.2 Piezoelectric coefficients

Only the non-zero entries. `e` in C/m², `d` in pC/N. Columns are the reachable Voigt set:
`eps_xx, eps_yy, eps_xy` on the rigid path, with `eps_zz` inserted third on the deformable
one. `d` is `e S`, with `d_direct` agreeing to the percentages in section 3.

| | illustrative `e` | illustrative `d` | `pvdf-dft-fit` `e` | `pvdf-dft-fit` `d` |
|---|---|---|---|---|
| PVDF beta | `e_y6 = +0.0564` | `d_y6 = +18.17` | `e_y6 = +0.0071` | `d_y6 = +2.46` |
| PVDF alpha | `e_y6 = −0.0204` | `d_y6 = −1.70` | `e_y6 = +0.0146` | `d_y6 = +2.79` |
| PVDF gamma | `e_x1 = +0.0277`, `e_x2 = −0.0275`, `e_x6 = +0.0170`, `e_y1 = −0.0155`, `e_y2 = +0.0154`, `e_y6 = −0.0095` | `d_x1 = +1.20`, `d_x2 = −1.19`, `d_x6 = +1.49`, `d_y1 = −0.67`, `d_y2 = +0.67`, `d_y6 = −0.83` | `e_x1 = +0.0081`, `e_x2 = −0.0080`, `e_x6 = +0.0053`, `e_y1 = −0.0043`, `e_y2 = +0.0043`, `e_y6 = −0.0028` | `d_x1 = +0.92`, `d_x2 = −0.91`, `d_x6 = +1.21`, `d_y1 = −0.49`, `d_y2 = +0.49`, `d_y6 = −0.65` |
| PE | all zero (< 1e-12 C/m²) | all zero (< 1e-11 pC/N) | all zero | all zero |

`pvdf-dft-valence`, deformable, full rows (`e` in C/m², rows `x, y, z`, columns
`eps_xx eps_yy eps_zz eps_xy`):

| | `e_x` | `e_y` | `e_z` |
|---|---|---|---|
| PVDF beta | 0, 0, 0, 0 | 0, 0, 0, **+0.01452** | 0, 0, 0, 0 |
| PVDF alpha | **+0.03844, −0.02522, +0.09915**, 0 | 0, 0, 0, **−0.00644** | 0, 0, 0, 0 |
| PVDF gamma | **+0.02920, −0.01081, −0.09725, +0.00502** | **+0.01896, +0.01730, −0.18995, −0.01789** | **−0.00693, +0.00338, +0.13565, +0.00714** |
| PE | 0, 0, 0, 0 | 0, 0, 0, 0 | 0, 0, 0, 0 |

and `d = e S` in pC/N, same layout:

| | `d_x` | `d_y` | `d_z` |
|---|---|---|---|
| PVDF beta | 0, 0, 0, 0 | 0, 0, 0, **+4.319** | 0, 0, 0, 0 |
| PVDF alpha | **+4.091, −4.391, +0.522**, 0 | 0, 0, 0, **−1.182** | 0, 0, 0, 0 |
| PVDF gamma | **+4.118, −1.568, −1.499, −0.866** | **+3.960, +1.034, −2.930, −7.637** | **−2.452, +0.436, +1.923, +4.004** |
| PE | 0, 0, 0, 0 | 0, 0, 0, 0 | 0, 0, 0, 0 |

Every "0" above is a machine zero: the largest of beta's diagonal `e` entries is 1e-11 C/m²
and PE's whole tensor is below 1e-14.

**Beta's diagonal columns are still exactly zero, and the reason has changed.** It is no
longer "the chain is rigid" — the chain deforms; `c` moves. It is that for a **planar
all-trans zigzag with bond-charge-increment charges, the cell dipole is exactly independent
of the backbone angle.** With BCI charges the dipole is a sum over bonds of `delta_ij
(r_i - r_j)`; a backbone C-C bond carries no increment; every C-H and C-F bond has a fixed
length; and the mirror symmetry of the zigzag pins the bisector of each pendant pair
perpendicular to the chain axis whatever the backbone angle is. So the one internal
coordinate beta and PE have moves `c` and leaves `mu` bit-for-bit where it was — measured
directly over a ±2° sweep of the shape parameter, `|Δmu| < 1e-10 e·A` while `c` moves by
2.3%, and asserted in
`tests/test_mechanics.py::test_a_planar_zigzag_cannot_change_its_dipole_by_bending`. Alpha
and gamma are helices, their conformations have no such mirror, and their diagonal columns
are the first non-zero ones this package has produced.

### 5.3 d33 and d31 for beta-PVDF: the answer, and its sign

Measured for a poled, uniaxially oriented beta-PVDF film, `d_33 = −32`, `d_31 = +20`,
`d_32 = +1.5` pC/N (section 8.4). In the film's convention axis 3 is the poling direction —
the crystal's polar axis, `x` in the packer's frame here — and axis 1 is the draw direction,
the chain axis `z`. So the film's `d_33` is this table's `d_x,xx`, its `d_31` is `d_x,zz`
and its `d_32` is `d_x,yy`.

**All three are exactly zero in the proper coefficient**, for the geometric reason above.
What the model *does* have on those columns is the **dimensional (thickness) term**: the
dipole per *current* volume changes under a dilation even when the dipole itself does not,
by `-P_i`, and contracting that with the compliance gives the charge an electrode sees when
a crystal of fixed dipoles changes shape. That is the term the Broadhurst-Davis model puts
at about two thirds of PVDF's activity (section 8.5), and it is worth quoting because it is
the only thing here with the right units to be compared:

| beta-PVDF, `pvdf-dft-valence` | film `d_33` | film `d_32` | film `d_31` |
|---|---|---|---|
| proper (this module's `d`) | 0 | 0 | 0 |
| dimensional term (`d_improper`, poling axis along +P) | **−4.46** | −5.96 | −0.14 |
| measured (Nix and Ward 1986) | −32 | +1.5 | +20 |

**The sign of `d_33` comes out right and the magnitude is seven times too small; `d_31`'s
sign comes out wrong.** Both are consequences of the same thing and neither is a
prediction. The dimensional term is negative on *every* diagonal column by construction —
it is `-P` times a positive compliance sum for any stable crystal — so getting `d_33 < 0`
is not evidence that the mechanism is right, it is arithmetic. The real `d_31 = +20` is
positive, which the dimensional term alone cannot produce at all, and which is exactly the
literature's argument against the dimensional model being the whole story. The intrinsic
part it is missing — a dipole that changes with strain — needs the chain's pendant geometry
or its bond lengths to relax, and neither is a degree of freedom in this model.

`d_improper` is reported by `Piezoelectric` and labelled, in the code and in the tables, as
not being a piezoelectric constant. `d = e S` is the coefficient; this is the term that
would be added to it by a model this package does not implement.

### 5.4 Blocking stress, free strain, work density (E = 0.01 V/A = 100 MV/m)

| potential | | poling direction | free strain | blocking stress (MPa) | work density (kJ/m³) |
|---|---|---|---|---|---|
| illustrative | PVDF beta | ±y | `e_xy = ±1.82e−3` | ±5.6 | 2.56 |
| | PVDF alpha | ±y | `e_xy = ∓1.70e−4` | ∓2.0 | 0.087 |
| | PVDF gamma | (0.873, −0.488, 0) | `e_xx = +1.37e−4`, `e_yy = −1.37e−4`, `e_xy = +1.71e−4` | +3.2, −3.1, +1.9 | 0.299 |
| | PE | any | 0 | 0 | 0 |
| `pvdf-dft-fit` | PVDF beta | ±y | `e_xy = ±2.46e−4` | ±0.7 | 0.044 |
| | PVDF alpha | ±y | `e_xy = ±2.79e−4` | ±1.5 | 0.102 |
| | PVDF gamma | (0.882, −0.471, 0) | `e_xx = +1.04e−4`, `e_yy = −1.03e−4`, `e_xy = +1.37e−4` | +0.9, −0.9, +0.6 | 0.068 |
| | PE | any | 0 | 0 | 0 |
| `pvdf-dft-valence` | PVDF beta | ±y | `e_xy = ±4.32e−4` | ±1.5 | 0.157 |
| | PVDF alpha | ±x | `e_xx = ±4.09e−4`, `e_yy = ∓4.39e−4`, `e_zz = ±5.2e−5` | ±3.8, ∓2.5, ±9.9 | 0.799 |
| | PVDF gamma | (−0.397, −0.777, +0.488) | `e_xx = −5.91e−4`, `e_yy = +3.2e−6`, `e_zz = +3.81e−4`, `e_xy = +8.23e−4` | −3.0, −0.8, +25.3, +1.5 | 3.16 |
| | PE | any | 0 | 0 | 0 |

The poling direction is the leading eigenvector of `d C d^T`, whose overall sign is
arbitrary; flipping the field flips the free strain and the blocking stress together and
leaves the work density alone, hence the ±. The full elastic triangle is twice the work
density. Gamma's 3.16 kJ/m³ is the largest figure in the set and it is the deformable path
that produces it: the same crystal gives 0.30 kJ/m³ rigid under the illustrative potential
and 0.068 under `pvdf-dft-fit`, because an axial channel it could not previously use is
where most of the work now comes from (`e_zz` and a blocking stress of 25 MPa along the
chain).

### 5.5 The polyethylene null

PE's atoms do carry charges (C −0.12, H +0.06 e). What vanishes is the cell dipole, and it
vanishes for *every* configuration, not only at the minimum: each CH2 group is neutral and
the all-trans repeat is centrosymmetric. The cancellation is exact to 1.7e−16 e·A over
random cells. So every piezoelectric number comes out at machine precision zero — `e` below
1e−12 C/m², `d` below 1e−11 pC/N by both routes on both paths (2e−12 deformable), work
density below 1e−24 kJ/m³ — while the elastic constants are
perfectly ordinary, `C_33 = 286.8` GPa included. It is a real test of the whole chain
(dipole, frame rotation, both derivative routes, and now the constrained relaxation as
well) rather than a multiplication by zero, and it is asserted on both paths in
`tests/test_mechanics.py`.

## 6. How converged these numbers are, and by what

Three knobs, and they do not matter equally.

**The strain step barely matters, and neither does the relaxation tolerance.** Between
`deps = 5e-4` and `4e-3`, `C_11`, `C_12`, `C_66` and `e` move by well under 1%, and the
pre-symmetrisation asymmetry of the stiffness block stays under 0.05 GPa for beta and 0.2 GPa
for PE. `pack.polish` stops at `gtol = 1e-5`, which looked like a loose tolerance to be taking
second derivatives through — the leftover gradient enters the stress difference divided by the
step — so it was checked rather than assumed: restarting L-BFGS-B from its own answer, once or
three more times, moves no constant of beta or PE by 0.002 GPa even at `deps = 5e-4`.

**The constrained relaxation is not the limit either.** Its own convergence measure — the
component of the shape gradient not along the constraint, divided by the whole — runs at
1e-16 for every state in section 5, the constraint itself is satisfied to 1e-15 relative,
and the two independent routes to `C_33` agree to 0.04%.

**The nonbonded cutoff is the accuracy limit, and it has two separate faults.** They need
separate fixes and it took measuring both to see that.

*Fault one: the force is discontinuous.* The Lennard-Jones term is energy-shifted but not
force-shifted, so its force jumps at `r = rc`. Straining the cell walks pairs across that
jump, and the finite difference sees a step the analytic derivative does not.
`CrystalPacker(lj_cutoff="force")` subtracts the tangent instead of the value —
`V(r) - V(rc) - (r - rc) V'(rc)` — so value and force both reach zero at `rc`. Since the
damped-shifted-force Coulomb term already carries a linear coefficient, the LJ force shift
rides in the same array and costs nothing at run time. Spread in `C_22` (largest minus
smallest over `deps = 5e-4, 2e-3, 4e-3`), GPa:

| | rc = 8, energy-shifted | rc = 8, force-shifted | rc = 12, energy-shifted | rc = 12, force-shifted |
|---|---|---|---|---|
| PVDF beta | **3.725** | 0.432 | 0.115 | 0.047 |
| PVDF alpha | 0.827 | 0.117 | 0.061 | 0.019 |
| PVDF gamma | 0.174 | 0.113 | 0.040 | 0.049 |
| PE | 0.092 | 0.032 | 0.345 | 0.022 |

Force-shifting improves it everywhere except gamma at 12 A, where both are already at the
0.04 GPa floor of the relaxation. PE at 12 A energy-shifted is the worst cell of the set
(0.345 GPa) and force-shifting takes it to 0.022, so this is not a beta-specific fix.

*Fault two: the tail is truncated,* and force-shifting does not help that at all — it is
slightly worse, because it removes a tail as well as a value. Beta, illustrative potential,
`deps = 2e-3`:

| cutoff (A) | shift | relaxed a, b | E/monomer | `C_11` | `C_22` | `C_12` | `C_66` |
|---|---|---|---|---|---|---|---|
| 8 (package default) | energy | 4.5973, 8.5867 | −8.153 | 44.42 | 37.77 | 2.25 | 3.10 |
| 12 | energy | 4.5882, 8.5532 | −9.181 | 45.72 | 42.09 | 2.56 | 3.36 |
| 16 | energy | 4.5850, 8.5465 | −9.407 | 46.55 | 42.73 | 2.58 | 3.38 |
| 20 | energy | 4.5840, 8.5445 | −9.486 | 46.82 | 42.95 | 2.60 | 3.41 |
| 8 | force | 4.6095, 8.6125 | −7.331 | 42.11 | 37.08 | 2.13 | 2.97 |
| 12 | force | 4.5917, 8.5604 | −8.906 | 44.90 | 41.25 | 2.56 | 3.35 |
| 16 | force | 4.5865, 8.5496 | −9.287 | 46.23 | 42.36 | 2.58 | 3.39 |
| 20 | force | 4.5847, 8.5460 | −9.424 | 46.63 | 42.76 | 2.59 | 3.40 |

Both forms converge to the same place. At the package default `C_22` is 12% low; at 12 A it
is 2% low. **So the force shift fixes the derivative and the cutoff length fixes the value,
and the recommendation is to use both**: `reference_from_chain(..., cutoff=12.0,
lj_cutoff="force")` gives a step spread of 0.05 GPa and a `C_22` within 4% of the 20 A
answer.

**What adopting either as a default would cost.** Both change the lattice energy, which is
why both are opt-in. Force-shifting at the default cutoff moves E/monomer by +0.82 kcal/mol
for beta (10%) and lengthening the cutoff to 12 A moves it by −1.03; the cells move by 0.3%
and 0.2% respectively. The relative quantity the package is actually judged on survives
both: alpha minus beta per monomer is +1.008 kcal/mol at the default, +1.032 at cutoff 12,
+0.965 force-shifted at 8, +1.019 force-shifted at 12 — a 4% spread on a number whose
acceptance window is much wider than that. What would have to be re-measured to make either
the default: the DESIGN 5.2 packing table and 5.3 funnel run (cell parameters and energies
per monomer for every polymorph), the fitted presets themselves (`docs/DFT_FIT.md`,
`docs/VALENCE_FIT.md` — their objectives contain lattice energies), the acceptance tests in
`polyfind.fitting`, and every table in this document. That is a large re-measurement for a
gain that is real but not urgent, so the default is unchanged and the knobs are documented
instead.

## 7. Cost

One full response calculation on one CPU core, wall clock, with pack and refine paid
separately:

| | beta | alpha | gamma | PE |
|---|---|---|---|---|
| response, rigid path | 0.32 s | 1.44 s | 7.8 s | 0.56 s |
| response, deformable path | 1.2 s | 6.7 s | 53 s | 1.3 s |
| kernel rows, deformable | 124 | 273 | 725 | 137 |
| pack + refine (once, beforehand) | 1.7-2.9 s | 2.7-5.0 s | 10-18 s | 1.7-3.2 s |

The rigid path is affordable because of the analytic gradient: `polish(...,
gradient="analytic")` spends one kernel row per function evaluation where the
finite-difference route spent `1 + 2 n_free`, and the stress itself is closed-form. The
deformable path costs a constrained SLSQP minimisation per strain state — nine for the
stiffness and `e`, plus two or three per field direction for the converse route — and each
objective evaluation is one kernel row plus one batched chain build. It scales with the
number of shape parameters (1 for beta and PE, 3 for alpha, 5 for gamma), which is why gamma
is 40 times beta. The shape gradient itself is free of kernel rows: the chain build is
differenced, not the energy, exactly as `refine_crystal` does it. Two hundred to seven
hundred kernel rows for a whole deformable response, against a hundred for a rigid one.

## 8. Literature, for order of magnitude only

Nothing above is a prediction of experiment, for reasons the package already records: the
illustrative potential is not quantitative (DESIGN.md section 5.4), and the fitted presets
are fitted to single-chain DFT energies and forces with no elastic or piezoelectric data in
the objective at all (docs/DFT_FIT.md, docs/VALENCE_FIT.md). The comparison below is a
sanity check on order of magnitude and nothing more.

Two conventions to keep straight before reading across. The published transverse numbers for
polymer crystals are usually **Young's moduli** `E_a`, `E_b`, not stiffness constants; `1/S_11`
is the comparable quantity here. On the rigid path that is an **axially clamped** modulus and
so an upper bound on the free `E_a`; on the deformable path the compliance includes `eps_zz`
and `1/S_11` is the free modulus. And the published *chain* modulus is `E_c = 1/S_33`, not
`C_33`; both are quoted below and they differ by 2-16%.

Engineering moduli, GPa. Illustrative | `pvdf-dft-fit` (both axially clamped) |
`pvdf-dft-valence` (free):

| | `E_a = 1/S_11` | `E_b = 1/S_22` | `E_c = 1/S_33` | `G_ab = 1/S_66` |
|---|---|---|---|---|
| PVDF beta | 44.3 \| 21.4 \| 21.7 | 37.7 \| 15.3 \| 17.4 | — \| — \| **310.5** | 3.10 \| 2.90 \| 3.36 |
| PVDF alpha | 26.3 \| 12.4 \| 12.5 | 20.7 \| 9.9 \| 9.5 | — \| — \| **244.2** | 12.03 \| 5.24 \| 5.45 |
| PVDF gamma | 26.3 \| 10.6 \| 10.2 | 28.6 \| 11.7 \| 11.2 | — \| — \| **77.7** | 9.04 \| 4.27 \| 3.85 |
| PE | 15.8 \| 6.3 \| 6.9 | 17.7 \| 6.2 \| 6.9 | — \| — \| **274.9** | 12.75 \| 4.52 \| 5.77 |

### 8.1 Beta-PVDF transverse stiffness

*Literature (calculated, DFT/PBE, VASP):* `C_11 = 19.83`, `C_22 = 24.50`, `C_33 = 287.33` GPa
from a stress-strain fit; `9.19`, `11.30`, `270.81` GPa from an energy-strain fit — the
factor-of-two spread between the two routes is the paper's own, not a transcription error.
Bulk modulus 11.18 GPa (PBE) / 25.25 GPa (PBE-D2).

> Pei, Y.; Zeng, X. C. "Elastic properties of poly(vinylidene fluoride) (PVDF) crystals: a
> density functional theory study." *J. Appl. Phys.* **109**, 093514 (2011).
> [10.1063/1.3574653](https://doi.org/10.1063/1.3574653)

*Against this package:* the illustrative potential's 44 / 38 GPa is 2-5x too stiff. Both
fitted presets land within about 1.5x of the stress-strain DFT numbers and bracket the
energy-strain ones — `pvdf-dft-fit` at 22 / 15 and `pvdf-dft-valence` at 22 / 18, which is a
better showing than either fit deserves, since no elastic data entered either objective.
Neither is a prediction; both are the right order. **Verdict: order of magnitude confirmed,
nothing finer.**

*Not found:* `C_12`, `C_13`, `C_23`, `C_44`, `C_55`, `C_66` for beta-PVDF. Pei and Zeng put the
full 6x6 matrices in AIP supplementary deposit E-JAPIAU-109-105107, which could not be
retrieved, so `C_12`, `C_13`, `C_23` and `C_66` have nothing to be compared against — and
`C_13 = 12.6` GPa is one of the new numbers, so this is a gap that matters more than it did.
The lattice-dynamics paper that would give both elastic *and* piezoelectric constants —
Tashiro, K.; Kobayashi, M.; Tadokoro, H.; Fukada, E., *Macromolecules* **13**, 691 (1980),
[10.1021/ma60075a040](https://doi.org/10.1021/ma60075a040) — is paywalled and was not read.

### 8.2 The chain-direction modulus, which is now a number

*Literature, beta-PVDF:* 287 GPa (DFT/PBE), 341 GPa (PBE-D2), against a **measured 177 GPa**
(Pei and Zeng, quoting K. Kaji, PhD thesis, Kyoto University, 1970 — read in Pei and Zeng, not
in Kaji).

*Literature, orthorhombic PE `E_c`:* 235 GPa (room temperature) and 254 GPa (−155 °C) by X-ray;
255 GPa (−165 °C) by X-ray; 260-358 GPa by Raman; 329 GPa by inelastic neutron scattering;
305-360 GPa across a dozen calculations, and 333 GPa from the DFT study below.

> Kurita, T.; Fukuda, Y.; Takahashi, M.; Sasanuma, Y. "Crystalline moduli of polymers,
> evaluated from density functional theory calculations under periodic boundary conditions."
> *ACS Omega* **3**, 4824 (2018).
> [10.1021/acsomega.8b00506](https://doi.org/10.1021/acsomega.8b00506) — its Table 3 is the
> compilation the figures above are taken from. Primary measurements therein: Nakamae, K.;
> Nishino, T.; Ohkubo, H., *J. Macromol. Sci. B* **30**, 1 (1991),
> [10.1080/00222349108245780](https://doi.org/10.1080/00222349108245780); Clements, J.;
> Jakeways, R.; Ward, I. M., *Polymer* **19**, 639 (1978),
> [10.1016/0032-3861(78)90116-7](https://doi.org/10.1016/0032-3861(78)90116-7).

*Against this package* (`pvdf-dft-valence`, deformable):

| | this package `E_c` | this package `C_33` | published |
|---|---|---|---|
| beta-PVDF | **310.5** | 317.6 | 287 (PBE), 341 (PBE-D2); measured 177 |
| PE | **274.9** | 286.8 | 235-255 X-ray, 329 neutron, 305-360 calculated, 333 DFT |
| alpha-PVDF | 244.2 | 250.1 | not found |
| gamma-PVDF | 77.7 | 92.4 | not found |

Beta sits between the two DFT values and 1.75x above the one measurement; PE sits inside the
measured range and below the calculated one. **Verdict: confirmed to the accuracy the
comparison can support**, which is not much — the calculated and measured chain moduli of PE
differ from each other by 50%, so anything in 230-360 GPa agrees with something.

**Read this before quoting either number.** Three things make it an upper bound and none of
them is small: bond lengths cannot stretch (in a real chain, stretching is roughly half the
axial compliance and its absence stiffens the answer), the line group admits only symmetric
conformations, and for beta and PE the chain has exactly one shape parameter, so the axial
deformation follows a single path rather than relaxing over a manifold. That the numbers
land near the DFT ones is therefore partly luck: a missing softening channel and a
transferable-but-not-fitted-here potential can cancel. What can be said without luck is
what section 4.3 measures — the number is no longer a report of an invented constant, it is
stationary, it is measured twice by different routes, and it does not move when the invented
constant is swept over a factor of two on either side.

### 8.3 Transverse moduli of polyethylene

*Literature (measured):* `E_a = 3.1`, `E_b = 3.8` GPa by X-ray at 20 °C (Sakurada, Ito and
Nakamae, via Kurita's Table 3); `E_a` and `E_b` both about 6 GPa by neutron scattering at 25 °C
(Holliday, 1971); `E_a = 9`, `E_b = 8` GPa by neutron at −196 °C (Twisleton, 1982).
*Calculated:* 10.9 / 7.8 GPa (Kurita 2018, B3LYP-D, 0 K); 6.9 / 8.6 GPa (Tashiro 1978);
5.9-13.7 GPa across older lattice-dynamics work.

*Against this package:* the illustrative potential gives 15.8 / 17.7 GPa, roughly 5x the X-ray
values and 1.5-2x the stiffest calculation. `pvdf-dft-fit` gives 6.3 / 6.2 and
`pvdf-dft-valence` 6.9 / 6.9, both sitting on top of the neutron measurement — and both
agreements should be **discounted entirely**, because the same presets squeeze PE's ab
cross-section to 26.6 A² and 29.3 A² against an experimental 36.7 (7.42 × 4.95; the a/b
labels are the packer's and need not match the crystallographic setting, so the area is the
comparison that means anything). A potential fitted to PVDF chains that gets PE's cell 20-28%
too dense has not earned a transverse modulus. **Verdict: unverified, and the apparent
agreement is not evidence.**

*Not found:* the Karasawa, Dasgupta and Goddard `C_ij` table (*J. Phys. Chem.* **95**, 2260
(1991), [10.1021/j100159a031](https://doi.org/10.1021/j100159a031)) is paywalled with no open
mirror, and no readable source gave PE's `C_66`. So the shear constant, the one place this
parametrisation is naturally strongest, has nothing to be checked against.

### 8.4 Piezoelectric coefficients

*Literature (measured, poled uniaxially oriented beta-PVDF film):* `d_31 = 20`, `d_32 = 1.5`,
`d_33 = −32`, `d_15 = −27`, `d_24 = −23` pC/N.

> Nix, E. L.; Ward, I. M. "The measurement of the shear piezoelectric coefficients of
> polyvinylidene fluoride." *Ferroelectrics* **67**, 137 (1986).
> [10.1080/00150198608245016](https://doi.org/10.1080/00150198608245016) — as tabulated in
> "Properties and applications of the beta phase poly(vinylidene fluoride)", *Polymers* **10**,
> 228 (2018), [10.3390/polym10030228](https://doi.org/10.3390/polym10030228). Corroborated at
> `d_31 = 20-28` pC/N by Harrison, J. S.; Ounaies, Z., "Piezoelectric polymers",
> NASA/CR-2001-211422 (2001).

*Against this package:* section 5.3 has the comparison in full. `d_33` and `d_31` are
identically zero in the proper coefficient, on every path and under every preset, and the
dimensional term that a Broadhurst-Davis model would add gives −4.5, −6.0 and −0.14 pC/N for
`d_33`, `d_32`, `d_31`. **Verdict: `d_33`'s sign reproduced but by arithmetic rather than by
mechanism (the dimensional term is negative on every diagonal column for any stable
crystal), its magnitude 7x too small, and `d_31`'s sign wrong.** `d_15` and `d_24` are among
the constants this cell cannot express at all (section 1).

The largest coefficient this module does produce for beta is `d_y6 = +18.2` pC/N
(illustrative), 2.5 (`pvdf-dft-fit`), 4.3 (`pvdf-dft-valence`) — the same order as a measured
`d_31`, and **not the same coefficient**. `d_y6` couples a transverse field to a shear.
Reading the agreement in magnitude as a match would be a category error.

*Not found:* `e_31` for PVDF in C/m² from any readable source, so the `e` tables have no
literature counterpart at all. Calculated single-crystal piezoelectric constants also could not
be verified: Nakhmanson, S. M.; Buongiorno Nardelli, M.; Bernholc, J., *Phys. Rev. Lett.* **92**,
115504 (2004), [10.1103/PhysRevLett.92.115504](https://doi.org/10.1103/PhysRevLett.92.115504)
is the right paper and was not opened; a search-engine summary quoting `e_33 = −0.332` for
beta-PVDF is **not** cited here because it could not be checked. Karasawa, N.; Goddard, W. A.
III, *Macromolecules* **25**, 7268 (1992) reports calculated elastic, dielectric *and*
piezoelectric constants for nine PVDF structures and is the single most useful paper for this
section; it is paywalled and was not read.

### 8.5 The mechanism, and why beta's diagonal columns are empty

*Literature:* the classic decomposition of PVDF's piezoelectric activity attributes about two
thirds to the **thickness (dimensional) effect** — the electrodes moving in the field of
constant crystal dipole moments — and the remaining third to a genuine change in the film's
dipole moment at constant thickness.

> Broadhurst, M. G.; Davis, G. T. "Physical basis for piezoelectricity in PVDF."
> *Ferroelectrics* **60**, 3 (1984).
> [10.1080/00150198408017504](https://doi.org/10.1080/00150198408017504); the underlying model
> is Broadhurst, M. G.; Davis, G. T.; McKinney, J. E.; Collins, R. E., *J. Appl. Phys.* **49**,
> 4992 (1978). Not an uncontested modern consensus: reviews since argue the dimensional model
> overestimates `d_33` for P(VDF-TrFE) and favour crystalline electrostriction or a
> crystalline-amorphous interfacial term instead.

*Against this package:* the dimensional term is present and is exactly the difference between
the proper and improper coefficients, quantified in section 5.3. The intrinsic third needs
the dipole to change with strain, and **the deformable path shows precisely which degrees of
freedom that requires.** Letting the backbone angles move is not enough for a planar zigzag:
the dipole is exactly invariant under it (section 5.2). It *is* enough for a helix, where
alpha and gamma get diagonal `e` entries of 0.1 to 0.19 C/m² and diagonal `d` of 0.5 to 2.9
pC/N. So the intrinsic channel is not absent from the model in principle — it is absent from
beta in particular, because beta's chain has no internal coordinate that changes its dipole.
Reaching it would need the pendant geometry or the bond lengths to relax, or polarizable
charges, none of which this model has.

## 9. What is not computed

* `C_44`, `C_55`, `C_45`, and `C_14`, `C_15`, `C_24`, `C_25`, `C_34`, `C_35`, `C_46`,
  `C_56`: the cell cannot express `eps_yz` or `eps_xz` — section 1. This is unchanged and
  will stay unchanged until the cell gains a general third lattice vector.
* `C_33`, `C_13`, `C_23`, `C_36` **under the illustrative potential and `pvdf-dft-fit`**:
  those potentials have no valence terms, so the deformable path refuses to run rather than
  reporting a number that would be a report of `refine_crystal`'s invented constant.
  `axial_report` measures what such a number would be made of and `electromechanical_response`
  reports `nan`.
* The bond-stretch contribution to `C_33`, roughly half of a real chain modulus:
  `polymer.bond_length` is a constant and `build_chain` places every atom at it, so the
  stretch terms are evaluated but contribute a constant and a zero gradient. Every axial
  constant here is an upper bound for that reason alone.
* Conformations outside the repeat's line group. The shape parameters are the line group's,
  so a symmetry-breaking relaxation under strain is not available — another way in which the
  constants are upper bounds, and the reason beta and PE have exactly one internal degree of
  freedom.
* The full compliance tensor, hence `d` as a truly free coefficient. `S` is the inverse of
  the reachable block (3×3 rigid, 4×4 deformable), so `d` is clamped in whatever is outside
  it.
* Electronic polarizability, depolarisation, and any field-induced change in the charges:
  the model is fixed point charges, so the dielectric constant is 1 by construction and the
  piezoelectric response has no electronic contribution.
* Any change in the *molecular* dipole with strain for a planar zigzag. Section 5.2 shows
  this is exact rather than approximate for bond-charge-increment charges, and it is what
  puts beta-PVDF's `d_33` and `d_31` out of reach.
* Temperature. Everything is a zero-kelvin second derivative; no phonons, no thermal
  expansion, no pyroelectric coefficient.
* Electrostriction. The response reported is strictly linear in the field; DESIGN.md section
  5.5's observation that a 0.5 V/A transverse field opens beta's `b` from 8.61 to 8.85 A is a
  quadratic effect and is not in any `d` here.
* Domain switching. A field opposing `P` finds the 180-degree-rotated cell at equal energy
  (DESIGN.md section 5.5); the linear response about one minimum says nothing about the path
  between them or about the coercive field.
