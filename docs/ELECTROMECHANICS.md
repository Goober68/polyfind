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
  `Relaxed.residual`, a measure of how converged each constrained relaxation is. It is a
  ratio, so it only means something where the shape actually moves: 1e-8 for alpha's axial
  column, 1e-3 for its in-plane ones, and meaningless (0/0) at the reference and at states
  whose symmetry forbids a conformational response. Section 6 says what was checked instead.
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
| PVDF beta | 0.65% | 0.52% | 0.03% |
| PVDF alpha | 0.37% | 0.51% | 0.72% |
| PVDF gamma | 0.17% | 0.59% | 0.52% |
| PE | both routes zero, difference < 5e-13 pC/N | both zero, < 7e-12 pC/N | both zero, < 2e-11 pC/N |

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
| PVDF beta | −6.59% … +2.30% | −4.97 GPa |
| PVDF alpha | −11.98% … +1.87% | −5.13 GPa |
| PVDF gamma | −10.51% … +3.28% | −3.92 GPa |
| PE | −5.76% … +3.57% | −2.27 GPa |

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
| PVDF beta | 89.4 | 251.1 | 412.9 | 736.4 | 78% |
| PVDF alpha | 218.1 | 282.2 | 346.4 | 474.8 | 37% |
| PVDF gamma | 996.7 | 1057.1 | 1117.5 | 1238.3 | 11% |
| PE | 60.5 | 220.5 | 380.5 | 700.4 | 84% |

**(c) The reference was axially stress-free *only because of* that restraint.** Along the
same path, `sigma_zz` (GPa):

| k | 0 | 52.5 | 105 | 210 |
|---|---|---|---|---|
| PVDF beta | −5.74 | −2.86 | **+0.02** | +5.78 |
| PE | −2.37 | −1.03 | **+0.30** | +2.97 |

The refinement's stationary point in `c` sat at `k = 105` (beta to 0.02 GPa, PE to 0.3 GPa)
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
| PVDF beta | 328.4079 | 328.4079 | 328.4079 | 328.4079 | 2e−5 GPa | 2.5469 → 2.5570 (+0.40%) | 2.54689 always |
| PVDF alpha | 260.4331 | 260.4327 | 260.4308 | 260.4307 | 2e−3 GPa | 4.6541 → 4.6352 (−0.41%) | 4.65408 always |
| PVDF gamma | 95.3393 | 95.3393 | 95.3392 | 95.3393 | 1e−4 GPa | 9.0516 → 9.0879 (+0.40%) | 9.01172 always |
| PE | 287.9797 | 287.9797 | 287.9797 | 287.9798 | 1e−4 GPa | 2.4928 → 2.5348 (+1.68%) | 2.49283 always |

`C_33` in GPa. The spread over the sweep is 6e−6 to 9e−4 *per cent* of the value — the
relaxation's own reproducibility, not a physical dependence — against a factor of 8.2 on the
rigid path. **The invented constant is no longer in the answer.** It cannot be: it enters
only through the starting structure, and the relaxation over `c` has a stationary point of
its own to find, which it reaches from every start in the sweep to five decimal places in
`c`.

Three further checks that this number is what it says it is.

* **The reference is axially stress-free by measurement, not by construction.** The
  residual `sigma_zz` at the relaxed reference is below 2e-5 GPa for all four crystals,
  against several GPa on the rigid path (section 4.1). The multiplier route is what reports
  it, so this is also a check on the multiplier.
* **`C_33` is measured twice by different routes** — the derivative of the constraint
  multiplier, and the curvature of the relaxed energy itself — and they agree to 0.004-0.07%
  (beta 328.408 / 328.395, alpha 260.433 / 260.284, gamma 95.339 / 95.322, PE 287.980 /
  288.173).
* **The axial row and column of the 4×4 are also measured twice**: `C_J3` from the analytic
  cell gradient over the axial column, `C_3J` from the multiplier over the in-plane columns.
  The largest disagreement before symmetrising is 0.004 GPa for beta, 0.11 for alpha, 0.054
  for gamma and 0.28 for PE — 0.1% or less of each `C_33`, but a tenth of PE's `C_66 = 1.11`,
  which is the one entry of the set that should be read as one significant figure.

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
| PVDF beta | 4.592, 8.593, 90.00, 2.6128 | 45.29 | 38.43 | 2.31 | 3.12 | 0.00 | 0.00 |
| PVDF alpha | 5.096, 8.840, 90.00, 4.8030 | 29.21 | 23.43 | 7.89 | 12.23 | 0.00 | 0.00 |
| PVDF gamma | 5.209, 8.887, 88.49, 9.5106 | 29.78 | 31.86 | 9.34 | 9.28 | −2.31 | 0.41 |
| PE | 4.750, 7.080, 87.27, 2.5753 | 27.46 | 30.46 | 16.14 | 15.83 | −4.57 | 4.28 |

`pvdf-dft-fit`, **rigid chain, `eps_zz` clamped**:

| | a, b, gamma, c (A, deg) | `C_11` | `C_22` | `C_12` | `C_66` | `C_16` | `C_26` |
|---|---|---|---|---|---|---|---|
| PVDF beta | 4.538, 8.627, 90.00, 2.6008 | 21.97 | 15.58 | 1.87 | 2.93 | 0.00 | 0.00 |
| PVDF alpha | 5.061, 8.997, 90.00, 4.7121 | 14.00 | 11.16 | 4.06 | 5.25 | 0.00 | 0.00 |
| PVDF gamma | 5.331, 8.868, 94.00, 9.3438 | 11.89 | 12.85 | 3.55 | 4.39 | 0.38 | 0.32 |
| PE | 3.917, 6.788, 89.94, 2.5339 | 36.22 | 32.43 | 16.69 | 24.04 | −14.76 | 11.48 |

`pvdf-dft-valence`, **deformable chain**, the full reachable 4×4:

| | a, b, gamma, c (A, deg) | `C_11` | `C_22` | `C_33` | `C_12` | `C_13` | `C_23` | `C_66` | `C_16` | `C_26` | `C_36` |
|---|---|---|---|---|---|---|---|---|---|---|---|
| PVDF beta | 4.596, 8.556, 90.00, 2.5469 | 22.79 | 17.92 | **328.4** | 2.21 | 12.72 | 2.30 | 3.41 | 0.00 | 0.00 | 0.00 |
| PVDF alpha | 5.069, 8.907, 90.00, 4.6541 | 15.10 | 11.33 | **260.4** | 5.27 | 0.49 | 8.10 | 5.53 | 0.00 | 0.00 | 0.00 |
| PVDF gamma | 5.453, 8.891, 84.00, 9.0117 | 12.69 | 12.49 | **95.3** | 3.70 | 10.50 | 4.47 | 4.22 | 0.28 | −0.43 | −4.18 |
| PE | 4.108, 7.130, 90.00, 2.4928 | 22.52 | 26.99 | **288.0** | 0.78 | 14.89 | 5.49 | 1.11 | 0.00 | 0.00 | 0.00 |

The overall sign of the polarization is the domain the packing happened to land in, and it
is arbitrary: flipping it flips every coefficient with an odd number of polar indices —
every `e` and `d` below, and none of the `C`. The signs are quoted as computed, and the one
place where a sign is a *result* rather than a convention (section 5.3) fixes the poling
axis along `+P` explicitly before quoting it.

Beta and alpha come out orthorhombic in these axes (`C_16 = C_26 = 0` to rounding), a
symmetry the calculation was not told about and recovers. Gamma and PE relax to
`gamma != 90`, so they have genuine monoclinic `C_16`, `C_26`, `C_36`.

**What the deformation itself changes**, same potential and same relaxed reference, rigid
against deformable:

| | `C_11` rigid → deformable | `C_22` | `C_12` | `C_66` | largest \|e\| (C/m²) |
|---|---|---|---|---|---|
| PVDF beta | 22.790 → 22.790 | 17.919 → 17.919 | 2.208 → 2.208 | 3.408 → 3.408 | 0.0149 → 0.0149 |
| PVDF alpha | 16.478 → 15.096 (−8%) | 12.098 → 11.332 (−6%) | 4.313 → 5.269 | 5.532 → 5.531 | 0.0077 → 0.1039 (×13) |
| PVDF gamma | 14.025 → 12.691 (−10%) | 13.982 → 12.493 (−11%) | 4.499 → 3.696 | 4.809 → 4.222 | 0.0154 → 0.1863 (×12) |
| PE | 22.514 → 22.515 | 26.991 → 26.991 | 0.778 → 0.778 | 1.110 → 1.111 | 0 → 0 |

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
| PVDF beta | `e_y6 = −0.0575` | `d_y6 = −18.41` | `e_y6 = −0.0067` | `d_y6 = −2.29` |
| PVDF alpha | `e_y6 = −0.0229` | `d_y6 = −1.87` | `e_y6 = −0.0157` | `d_y6 = −2.99` |
| PVDF gamma | `e_x1 = −0.0286`, `e_x2 = +0.0286`, `e_x6 = +0.0176`, `e_y1 = −0.0160`, `e_y2 = +0.0160`, `e_y6 = +0.0098` | `d_x1 = −1.23`, `d_x2 = +1.24`, `d_x6 = +1.53`, `d_y1 = −0.69`, `d_y2 = +0.70`, `d_y6 = +0.86` | `e_x1 = +0.0090`, `e_x2 = −0.0089`, `e_x6 = +0.0060`, `e_y1 = −0.0048`, `e_y2 = +0.0047`, `e_y6 = −0.0032` | `d_x1 = +1.01`, `d_x2 = −1.00`, `d_x6 = +1.36`, `d_y1 = −0.54`, `d_y2 = +0.54`, `d_y6 = −0.73` |
| PE | all zero (< 1e-12 C/m²) | all zero (< 1e-11 pC/N) | all zero | all zero |

`pvdf-dft-valence`, deformable, full rows (`e` in C/m², rows `x, y, z`, columns
`eps_xx eps_yy eps_zz eps_xy`):

| | `e_x` | `e_y` | `e_z` |
|---|---|---|---|
| PVDF beta | 0, 0, 0, 0 | 0, 0, 0, **−0.01490** | 0, 0, 0, 0 |
| PVDF alpha | **+0.03783, −0.02684, +0.10394**, 0 | 0, 0, 0, **−0.00770** | 0, 0, 0, 0 |
| PVDF gamma | **+0.02864, −0.01044, −0.09534, +0.00548** | **+0.01787, +0.01758, −0.18625, −0.01684** | **−0.00637, +0.00337, +0.13401, +0.00680** |
| PE | 0, 0, 0, 0 | 0, 0, 0, 0 | 0, 0, 0, 0 |

and `d = e S` in pC/N, same layout:

| | `d_x` | `d_y` | `d_z` |
|---|---|---|---|
| PVDF beta | 0, 0, 0, 0 | 0, 0, 0, **−4.372** | 0, 0, 0, 0 |
| PVDF alpha | **+4.114, −4.662, +0.536**, 0 | 0, 0, 0, **−1.392** | 0, 0, 0, 0 |
| PVDF gamma | **+3.833, −1.493, −1.372, −0.469** | **+3.448, +1.111, −2.681, −6.741** | **−2.176, +0.398, +1.783, +3.559** |
| PE | 0, 0, 0, 0 | 0, 0, 0, 0 | 0, 0, 0, 0 |

Every "0" above is a machine zero: the largest of beta's diagonal `e` entries is 8e-12 C/m²
and PE's whole tensor is below 2e-14.

**Beta's diagonal columns are still exactly zero, and the reason has changed.** It is no
longer "the chain is rigid" — the chain deforms; `c` moves. It is that for a **planar
all-trans zigzag with bond-charge-increment charges, the cell dipole is exactly independent
of the backbone angle.** With BCI charges the dipole is a sum over bonds of `delta_ij
(r_i - r_j)`; a backbone C-C bond carries no increment; every C-H and C-F bond has a fixed
length; and the mirror symmetry of the zigzag pins the bisector of each pendant pair
perpendicular to the chain axis whatever the backbone angle is. So the one internal
coordinate beta and PE have moves `c` and leaves `mu` where it was — measured directly over
a ±2° sweep of the shape parameter, over which `c` runs from −1.17% to +1.14% and beta's
`mu_x` stays at −0.7234964088 e·A in every digit printed (its transverse components sit at
1e-13, the numerical asymmetry of the chain build and not a response) — and asserted in
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
| dimensional term (`d_improper`, poling axis along +P) | **−4.43** | −5.90 | −0.14 |
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
| illustrative | PVDF beta | ±y | `e_xy = ∓1.84e−3` | ∓5.7 | 2.645 |
| | PVDF alpha | ±y | `e_xy = ∓1.87e−4` | ∓2.3 | 0.107 |
| | PVDF gamma | (0.872, 0.489, 0) | `e_xx = −1.41e−4`, `e_yy = +1.42e−4`, `e_xy = +1.76e−4` | −3.3, +3.3, +2.0 | 0.321 |
| | PE | any | 0 | 0 | 0 |
| `pvdf-dft-fit` | PVDF beta | ±y | `e_xy = ∓2.29e−4` | ∓0.7 | 0.038 |
| | PVDF alpha | ±y | `e_xy = ∓2.99e−4` | ∓1.6 | 0.117 |
| | PVDF gamma | (0.883, −0.470, 0) | `e_xx = +1.15e−4`, `e_yy = −1.14e−4`, `e_xy = +1.55e−4` | +1.0, −1.0, +0.7 | 0.084 |
| | PE | any | 0 | 0 | 0 |
| `pvdf-dft-valence` | PVDF beta | ±y | `e_xy = ∓4.37e−4` | ∓1.5 | 0.163 |
| | PVDF alpha | ±x | `e_xx = ±4.11e−4`, `e_yy = ∓4.66e−4`, `e_zz = ±5.4e−5` | ±3.8, ∓2.7, ±10.4 | 0.841 |
| | PVDF gamma | (−0.397, +0.775, −0.493) | `e_xx = +5.26e−4`, `e_yy = +7.3e−6`, `e_zz = −3.50e−4`, `e_xy = +7.16e−4` | +2.8, +0.8, −24.8, +1.4 | 2.80 |
| | PE | any | 0 | 0 | 0 |

The poling direction is the leading eigenvector of `d C d^T`, whose overall sign is
arbitrary; flipping the field flips the free strain and the blocking stress together and
leaves the work density alone, hence the ±. The full elastic triangle is twice the work
density. Gamma's 2.80 kJ/m³ is the largest figure in the set and it is the deformable path
that produces it: the same crystal gives 0.32 kJ/m³ rigid under the illustrative potential
and 0.084 under `pvdf-dft-fit`, because an axial channel it could not previously use is
where most of the work now comes from (`e_zz` and a blocking stress of 25 MPa along the
chain).

### 5.5 The polyethylene null

PE's atoms do carry charges (C −0.12, H +0.06 e). What vanishes is the cell dipole, and it
vanishes for *every* configuration, not only at the minimum: each CH2 group is neutral and
the all-trans repeat is centrosymmetric. The cancellation is exact to 1.7e−16 e·A over
random cells. So every piezoelectric number comes out at machine precision zero — `e` below
1e−12 C/m², `d` below 1e−11 pC/N by both routes on both paths (2e−11 deformable), work
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

**The constrained relaxation is not the limit either**, and it was checked in three ways
rather than trusted. The constraint itself is satisfied to 2e-14 relative or better on every
state. The component of the shape gradient not along the constraint runs at 1e-8 on alpha's
axial column and 1e-3 on its in-plane ones — that ratio is 0/0 at the reference and wherever
symmetry forbids a conformational response, so it is quoted where it means something and
ignored where it does not. And restarting a converged state from its own answer with three
times the iteration budget moves alpha's shape parameters by 1e-5 degrees, its stress by
1e-6 GPa and its dipole not at all.

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
| PVDF beta | **3.742** | 0.422 | 0.039 | 0.059 |
| PVDF alpha | 0.886 | 0.090 | 0.028 | 0.031 |
| PVDF gamma | 0.303 | 0.063 | 0.039 | 0.045 |
| PE | 0.092 | 0.030 | 0.054 | 0.022 |

The fix works where the fault is: at the default cutoff, where beta's `b = 8.59` and alpha's
lattice put pairs right on `rc`, force-shifting takes the spread down by a factor of nine
and ten. Everywhere else the spread is already 0.02-0.09 GPa, which is the relaxation's own
run-to-run scatter — repeating the sweep moves those entries by that much, so nothing should
be read into the differences between them.

*Fault two: the tail is truncated,* and force-shifting does not help that at all — it is
slightly worse, because it removes a tail as well as a value. Beta, illustrative potential,
`deps = 2e-3`:

| cutoff (A) | shift | relaxed a, b | E/monomer | `C_11` | `C_22` | `C_12` | `C_66` |
|---|---|---|---|---|---|---|---|
| 8 (package default) | energy | 4.5918, 8.5935 | −7.943 | 45.29 | 38.43 | 2.31 | 3.12 |
| 12 | energy | 4.5824, 8.5600 | −8.982 | 46.42 | 42.80 | 2.65 | 3.42 |
| 16 | energy | 4.5792, 8.5533 | −9.209 | 47.29 | 43.43 | 2.63 | 3.39 |
| 20 | energy | 4.5782, 8.5513 | −9.289 | 47.54 | 43.67 | 2.66 | 3.43 |
| 8 | force | 4.6039, 8.6192 | −7.115 | 42.50 | 37.56 | 2.16 | 2.97 |
| 12 | force | 4.5859, 8.5672 | −8.705 | 45.60 | 41.93 | 2.62 | 3.38 |
| 16 | force | 4.5808, 8.5564 | −9.089 | 46.95 | 43.07 | 2.65 | 3.41 |
| 20 | force | 4.5790, 8.5529 | −9.227 | 47.36 | 43.48 | 2.66 | 3.42 |

Both forms converge to the same place. At the package default `C_22` is 12% low; at 12 A it
is 2% low. **So the force shift fixes the derivative and the cutoff length fixes the value,
and the recommendation is to use both**: `reference_from_chain(..., cutoff=12.0,
lj_cutoff="force")` gives a step spread of 0.06 GPa and a `C_22` within 4% of the 20 A
answer.

**What adopting either as a default would cost.** Both change the lattice energy, which is
why both are opt-in. Force-shifting at the default cutoff moves E/monomer by +0.83 kcal/mol
for beta (10%) and lengthening the cutoff to 12 A moves it by −1.04; the cells move by 0.3%
and 0.2% respectively. The relative quantity the package is actually judged on survives
both: alpha minus beta per monomer is +0.929 kcal/mol at the default, +0.953 at cutoff 12,
+0.887 force-shifted at 8, +0.941 force-shifted at 12 — a 7% spread on a number whose
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
| response, rigid path | 0.8 s | 2.3 s | 18.8 s | 1.6 s |
| response, deformable path | 4.8 s | 18.5 s | 124 s | 2.3 s |
| kernel rows, deformable | 81 | 277 | 673 | 108 |
| pack + refine (once, beforehand) | 2.7-10 s | 4.9-13 s | 22-48 s | 3.6-6.8 s |

The rigid path is affordable because of the analytic gradient: `polish(...,
gradient="analytic")` spends one kernel row per function evaluation where the
finite-difference route spent `1 + 2 n_free`, and the stress itself is closed-form. The
deformable path costs a constrained SLSQP minimisation per strain state — nine for the
stiffness and `e`, plus two or three per field direction for the converse route — and each
objective evaluation is one kernel row plus one batched chain build. It scales with the
number of shape parameters (1 for beta and PE, 3 for alpha, 5 for gamma), which is why gamma
is 25 times beta. The shape gradient itself is free of kernel rows: the chain build is
differenced, not the energy, exactly as `refine_crystal` does it. Eighty to seven hundred
kernel rows for a whole deformable response, against a hundred for a rigid one. The timings
above were taken with other work on the same machine, so the ratios are steadier than the
absolute numbers.

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
| PVDF beta | 45.2 \| 21.8 \| 22.1 | 38.3 \| 15.4 \| 17.7 | — \| — \| **321.3** | 3.12 \| 2.93 \| 3.41 |
| PVDF alpha | 26.6 \| 12.5 \| 12.6 | 21.3 \| 10.0 \| 9.3 | — \| — \| **253.8** | 12.23 \| 5.25 \| 5.53 |
| PVDF gamma | 26.4 \| 10.9 \| 10.5 | 28.8 \| 11.8 \| 11.4 | — \| — \| **82.0** | 9.06 \| 4.37 \| 3.97 |
| PE | 15.8 \| 6.3 \| 21.7 | 17.8 \| 6.2 \| 26.9 | — \| — \| **277.2** | 12.76 \| 4.52 \| 1.11 |

(The rigid columns are `1/S` of the 3×3 in-plane block, so they are *axially clamped*; the
deformable column inverts the 4×4 and is the free modulus. PE's deformable column is not
comparable with its rigid ones at all: under `pvdf-dft-valence` PE packs into a different
cell from the one the other two potentials find — see section 8.3.)

### 8.1 Beta-PVDF transverse stiffness

*Literature (calculated, DFT/PBE, VASP):* `C_11 = 19.83`, `C_22 = 24.50`, `C_33 = 287.33` GPa
from a stress-strain fit; `9.19`, `11.30`, `270.81` GPa from an energy-strain fit — the
factor-of-two spread between the two routes is the paper's own, not a transcription error.
Bulk modulus 11.18 GPa (PBE) / 25.25 GPa (PBE-D2).

> Pei, Y.; Zeng, X. C. "Elastic properties of poly(vinylidene fluoride) (PVDF) crystals: a
> density functional theory study." *J. Appl. Phys.* **109**, 093514 (2011).
> [10.1063/1.3574653](https://doi.org/10.1063/1.3574653)

*Against this package:* the illustrative potential's 45 / 38 GPa is 2-5x too stiff. Both
fitted presets land within about 1.5x of the stress-strain DFT numbers and bracket the
energy-strain ones — `pvdf-dft-fit` at 22.0 / 15.6 and `pvdf-dft-valence` at 22.8 / 17.9,
which is a better showing than either fit deserves, since no elastic data entered either
objective. Neither is a prediction; both are the right order. **Verdict: order of magnitude
confirmed, nothing finer.**

*Not found:* `C_12`, `C_13`, `C_23`, `C_44`, `C_55`, `C_66` for beta-PVDF. Pei and Zeng put the
full 6x6 matrices in AIP supplementary deposit E-JAPIAU-109-105107, which could not be
retrieved, so `C_12`, `C_13`, `C_23` and `C_66` have nothing to be compared against — and
`C_13 = 12.7` GPa is one of the new numbers, so this is a gap that matters more than it did.
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
| beta-PVDF | **321.3** | 328.4 | 287 (PBE), 341 (PBE-D2); measured 177 |
| PE | **277.2** | 288.0 | 235-255 X-ray, 329 neutron, 305-360 calculated, 333 DFT |
| alpha-PVDF | 253.8 | 260.4 | not found |
| gamma-PVDF | 82.0 | 95.3 | not found |

Beta sits between the two DFT values (287 and 341) and 1.8x above the one measurement; PE
sits inside the measured range and below the calculated one. **Verdict: confirmed to the
accuracy the comparison can support**, which is not much — the calculated and measured chain
moduli of PE differ from each other by 50%, so anything in 230-360 GPa agrees with
something.

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

*Against this package:* the illustrative potential gives 15.8 / 17.8 GPa, roughly 5x the X-ray
values and 1.5-2x the stiffest calculation. `pvdf-dft-fit` gives 6.3 / 6.2, sitting on top of
the neutron measurement, and `pvdf-dft-valence` gives 21.7 / 26.9, three times the stiffest
calculation. Both should be **discounted entirely**, and the disagreement *between* them is
why: the two presets squeeze PE's ab cross-section to 26.6 A² and 29.3 A² against an
experimental 36.7 (7.42 × 4.95; the a/b labels are the packer's and need not match the
crystallographic setting, so the area is the comparison that means anything), and they do it
into two different cells — 3.92 × 6.79 at gamma 90 against 4.11 × 7.13 at gamma 90, with
`C_66` 24.0 against 1.11. A potential fitted to PVDF chains that gets PE's cell 20-28% too
dense has not earned a transverse modulus, and the axial one is the only number in this
column that is stable across the two packings (288.0 against 286.8 GPa, 0.4% apart, for two
cells whose transverse constants differ by a factor of three). **Verdict: unverified, and
the apparent agreement of one of them is not evidence.**

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
dimensional term that a Broadhurst-Davis model would add gives −4.43, −5.90 and −0.14 pC/N
for `d_33`, `d_32`, `d_31`. **Verdict: `d_33`'s sign reproduced but by arithmetic rather than
by mechanism (the dimensional term is negative on every diagonal column for any stable
crystal), its magnitude 7x too small, and `d_31`'s sign wrong.** `d_15` and `d_24` are among
the constants this cell cannot express at all (section 1).

The largest coefficient this module does produce for beta is `|d_y6| = 18.4` pC/N
(illustrative), 2.3 (`pvdf-dft-fit`), 4.4 (`pvdf-dft-valence`) — the same order as a measured
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
alpha and gamma get diagonal `e` entries of 0.10 to 0.19 C/m² and diagonal `d` of 0.5 to 2.7
pC/N. So the intrinsic channel is not absent from the model in principle — it is absent from
beta in particular, because beta's chain has no internal coordinate that changes its dipole.
Reaching it would need the pendant geometry or the bond lengths to relax, or polarizable
charges, none of which this model has.

*A caveat on the helices' diagonal columns.* Alpha's and gamma's largest diagonal `e` is the
`eps_zz` one, 0.10 and 0.19 C/m², and it comes from a conformation change under axial
strain. Their `d` from the two routes agree to 0.5-0.7%, and the entries are stable in sign
and to about 10% between the 1.54 A and 1.528 A backbone geometries (alpha `d_x,zz` 0.52 →
0.54, gamma `d_z,zz` 1.92 → 1.78). But alpha packs *polar* under this preset, which is one of
the three acceptance tests the fit fails (DESIGN.md section 5.10); a real alpha crystal is
antipolar and its diagonal response would cancel between the two chains. So alpha's numbers
describe the cell this potential produces, not the phase it is named after.

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
