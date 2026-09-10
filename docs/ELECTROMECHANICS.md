# Electromechanical response: what `polyfind.mechanics` can and cannot compute

`pack.py` says which polymorph is stable and how polar it is. This document is about the
other half of the measurement need: how much mechanical work the crystal does in a field.
`src/polyfind/mechanics.py` computes the elastic stiffness, the piezoelectric coefficients
in both directions, and the three numbers an actuator is judged by — blocking stress, free
strain, work density. `examples/electromechanics.py --axial` produces sections 4, 5 and 7
directly; the cutoff sweep in section 6 and the engineering moduli in section 8 are a few
lines on top of the same API (`elastic_constants` at several `cutoff`s, and `1/S_11` of
`Elastic.S`).

**The headline is a limitation, not a number.** Of the six Voigt strain components, three
are reachable, two cannot be expressed at all, and the sixth — the one along the chain —
is reachable geometrically but its energy is not the potential's. Everything else here is
downstream of that, so it comes first.

## 1. Which strains this cell can express

A cell in `pack.py` is

    a_vec = (a, 0, 0)   b_vec = (b cos gamma, b sin gamma, 0)   c_vec = (0, 0, c)

`c_vec` is pinned to z, and z is also the chain axis. A homogeneous strain sends the rows
to `r_i (1 + eps)`, so the third row becomes `c (eps_xz, eps_yz, 1 + eps_zz)`.

| Voigt | component | status | why |
|---|---|---|---|
| 1, 2, 6 | `eps_xx`, `eps_yy`, `eps_xy` | **reachable** | carried by `(a, b, gamma)`. The rigid rotation about z that puts `a_vec` back on x is absorbed exactly by the setting angles `phi1, phi2`. |
| 4, 5 | `eps_yz`, `eps_xz` | **not expressible** | these tilt `c_vec` off z. There is no cell variable for a third lattice vector that is not along the chain axis, so `C_44`, `C_55`, `C_45` and every mixed constant containing a 4 or a 5 are absent — not unconverged, absent. `strained_cell` raises rather than silently dropping them. |
| 3 | `eps_zz` | **geometrically reachable, energetically not computable** | `c` is not a cell variable at all: it is the rise of the chain's repeat transform. See section 4. |

So the reachable stiffness is the 3×3 block `{11, 22, 12, 16, 26, 66}`. Those are the
*true* constants and not a plane-stress reduction — `C_11 = dsigma_1/deps_1` is defined at
fixed `eps_2 … eps_6`, and fixed `eps_zz` is exactly what a rigid chain enforces. What
cannot be had from them is the full *compliance*: `S` is the inverse of the 3×3 block only,
so every `d` below is an **axially clamped** piezoelectric strain coefficient.

"Relaxed-ion" here means relaxed *chains*. The internal degrees of freedom of this
parametrisation are the rigid-body ones — each chain's setting angle and the axial offset
`dz` — and they are relaxed at every strain state by `pack.polish` on the analytic cell
gradient. The chain conformation is held rigid. Letting it move is not offered for the
in-plane block on purpose: it would change `c`, and the calculation would no longer be at
fixed `eps_zz`.

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
  columns, and for a rigid dipole array that difference is the *whole* of `dP/de_xx` —
  dilating a box of fixed dipoles changes `P` only through `V` and does nothing
  piezoelectric. `Piezoelectric` reports both and the difference is asserted to be `P_i` in
  `tests/test_mechanics.py`.
* Actuator at field `E`: free strain `eps_free = d^T E` at zero in-plane stress, blocking
  stress `sigma_block = C eps_free = e^T E`, the stress the crystal exerts on a rigid clamp
  (the external stress needed to *hold* `eps = 0` is its negative), work density
  `w = (1/4) sigma_block . eps_free` against a matched linear load. The field used for the
  tables is **0.01 V/A = 100 MV/m**, a realistic drive field for poled PVDF and about twice
  its coercive field; the figures are the linear response evaluated there, not a measurement
  at that field.
* Poling direction: the tables use the direction that maximises the work, the leading
  eigenvector of `d C d^T`. That is worth computing rather than assuming — see section 5.

One trap, recorded because it cost two wrong answers before it was found. A shear leaves
`a_vec` off the x axis, so the packer's canonical frame is rotated relative to the reference
by `theta = -e_xy/2`. The dipole has to be rotated back out of that frame before it is
differentiated, *and* an applied field has to be rotated into it. Each of those terms is
`P/2` — the same size as the coefficient being measured — so getting either wrong does not
degrade the answer, it replaces it. The direct/converse agreement in section 3 is what
caught both.

## 3. The two piezoelectric routes, and the check

`e` comes from the dipole of strained, internally relaxed cells. `d` comes from solving
`sigma(eps, E) = 0` — Newton on the analytic stress, the internal coordinates relaxed at
every iterate. The identity `d = e S` follows from `d^2 h / deps dE` being symmetric, and
the two sides are computed by disjoint machinery, so agreement is evidence.

| | illustrative | `pvdf-dft-fit` |
|---|---|---|
| PVDF beta | 0.65% | 0.50% |
| PVDF alpha | 0.43% | 0.35% |
| PVDF gamma | 1.03% | 1.54% |
| PE | both routes zero, difference < 1e-11 pC/N | both routes zero, difference < 1e-11 pC/N |

(Largest `|d_from_e - d_direct|` as a fraction of the largest coefficient. The residual is
finite-difference error in `e` at `deps = 2e-3` and Newton tolerance in `d`, not a
discrepancy in the physics.)

## 4. Axial strain: the answer is no, and here is the evidence

`c` *can* be moved. Bond lengths are rigid, but the backbone bond angles are not, and
shifting every backbone angle of the repeat together moves the repeat over a range of
several per cent in both directions:

| | reachable `eps_zz` (illustrative) | residual `sigma_zz` at fixed chain shape |
|---|---|---|
| PVDF beta | −6.5% … +2.4% | −4.51 GPa |
| PVDF alpha | −11.8% … +2.1% | −4.75 GPa |
| PVDF gamma | −10.3% … +3.5% | −3.59 GPa |
| PE | −5.8% … +3.6% | −2.27 GPa |

So the obstacle is not geometry. It is that **nothing in the potential resists the
deformation.** Three separate facts say so.

**(a) The kernel has no valence terms.** `pack.py` is Lennard-Jones + damped-shifted-force
Coulomb + a Fourier torsion term. `FFParameters.applied` forwards the fitted LJ, charges
and torsions into packing and *deliberately* does not forward `bond_terms` or `angle_terms`,
on the documented grounds that packing holds the chain rigid so they are an additive
constant. The only thing anywhere in the pipeline that resists opening an angle is
`refine_crystal`'s own harmonic restraint, at `angle_stiffness = 105` kcal/mol/rad²,
described in its docstring as "UFF-like for sp3 carbon" and excluded from the reported
lattice energy.

**(b) `C_33` is mostly that constant.** Sweeping it (illustrative potential, GPa):

| k (kcal/mol/rad²) | 0 | 52.5 | 105 | 210 | share of `C_33(105)` from k |
|---|---|---|---|---|---|
| PVDF beta | 83.8 | 242.4 | 401.1 | 718.4 | 79% |
| PVDF alpha | 203.4 | 266.1 | 328.8 | 454.3 | 38% |
| PVDF gamma | 792.6 | 851.9 | 911.1 | 1029.6 | 13% |
| PE | 60.5 | 220.5 | 380.5 | 700.4 | 84% |

`C_33` is exactly affine in an invented number. For the two all-trans cases — the ones
whose axial response is dominated by the angle, and the ones a chain modulus is usually
quoted for — four fifths of it is the constant.

**(c) The reference is axially stress-free *only because of* that restraint.** Along the
same path, `sigma_zz` (GPa):

| k | 0 | 52.5 | 105 | 210 |
|---|---|---|---|---|
| PVDF beta | −5.30 | −2.64 | **+0.03** | +5.35 |
| PE | −2.37 | −1.03 | **+0.02** | +2.97 |

The refinement's stationary point in `c` sits at `k = 105` and nowhere else. At `k = 0` the
crystal wants to contract along the chain at several GPa and keeps going; the lattice energy
alone has no axial minimum at all.

**(d) Bond stretching is absent entirely**, and in a real chain it is roughly half the
axial stiffness. Even a correctly fitted bend constant would leave that channel out, because
`polymer.bond_length` is a constant.

Two honest caveats on the table above. The path is a *uniform* shift of all backbone angles,
which for beta and PE is the direction the refinement actually relaxed along (hence the
`+0.03` and `+0.02`), but for alpha and gamma is only one direction of a larger shape
manifold — gamma's `sigma_zz` along it is −60 GPa, i.e. this is a constrained path and
gamma's `C_33` numbers are the stiffness of that constraint as much as of anything else. And
the torsions are held fixed along the path; letting them relax would soften every number.

**Verdict.** `electromechanical_response` reports `C_33 = nan`. `axial_report` returns
`computable=False` with the numbers above. A `C_33` of 400 GPa would have looked plausible
next to published values and would have meant nothing.

## 5. Results

### Elastic constants (GPa, relaxed-ion, rigid chain, `eps_zz` clamped)

Illustrative potential:

| | a, b, gamma, c (A, deg) | `C_11` | `C_22` | `C_12` | `C_66` | `C_16` | `C_26` |
|---|---|---|---|---|---|---|---|
| PVDF beta | 4.597, 8.587, 90.00, 2.630 | 44.42 | 37.76 | 2.25 | 3.10 | 0.00 | 0.00 |
| PVDF alpha | 5.096, 8.865, 90.00, 4.829 | 28.66 | 22.51 | 7.31 | 12.03 | 0.00 | 0.00 |
| PVDF gamma | 5.216, 8.892, 91.59, 9.566 | 29.56 | 31.58 | 9.15 | 9.25 | 2.19 | −0.48 |
| PE | 4.750, 7.080, 92.73, 2.575 | 27.46 | 30.46 | 16.14 | 15.83 | 4.57 | −4.29 |

`pvdf-dft-fit`:

| | a, b, gamma, c (A, deg) | `C_11` | `C_22` | `C_12` | `C_66` | `C_16` | `C_26` |
|---|---|---|---|---|---|---|---|
| PVDF beta | 4.544, 8.621, 90.00, 2.618 | 21.60 | 15.47 | 1.84 | 2.90 | 0.00 | 0.00 |
| PVDF alpha | 5.062, 9.006, 90.00, 4.740 | 13.87 | 11.10 | 4.04 | 5.24 | 0.00 | 0.00 |
| PVDF gamma | 5.322, 8.869, 93.84, 9.410 | 11.56 | 12.72 | 3.44 | 4.30 | 0.47 | 0.37 |
| PE | 3.917, 6.788, 90.06, 2.534 | 36.22 | 32.43 | 16.69 | 24.04 | 14.76 | −11.48 |

Beta and alpha come out orthorhombic in these axes (`C_16 = C_26 = 0` to rounding), which is
a symmetry the calculation was not told about and recovers. Gamma and PE relax to
`gamma != 90`, so they have genuine monoclinic `C_16`, `C_26`. PE under the fitted preset is
the outlier of the set: its ab cross-section collapses to 26.6 A² (3.92 × 6.79), against 33.6
(4.75 × 7.08) from the illustrative potential and 36.7 from experiment (7.42 × 4.95 — the a/b
labels are the packer's and need not match the crystallographic setting, so the area is the
comparison that means anything). That is 28% too dense, and the large `C_16`
follows the distorted cell. That preset is fitted to PVDF single-chain DFT; it is not a PE
potential and this is what that costs.

### Piezoelectric coefficients

Only the non-zero entries. `e` in C/m², `d` in pC/N, columns are `eps_xx`, `eps_yy`,
`eps_xy`; `d` is the axially clamped coefficient and the value shown is `e S`, with
`d_direct` agreeing to the percentages in section 3.

| | illustrative `e` | illustrative `d` | fitted `e` | fitted `d` |
|---|---|---|---|---|
| PVDF beta | `e_y6 = +0.0564` | `d_y6 = +18.17` | `e_y6 = +0.0071` | `d_y6 = +2.46` |
| PVDF alpha | `e_y6 = −0.0204` | `d_y6 = −1.70` | `e_y6 = +0.0146` | `d_y6 = +2.79` |
| PVDF gamma | `e_x1 = +0.0277`, `e_x2 = −0.0275`, `e_x6 = +0.0170`, `e_y1 = −0.0155`, `e_y2 = +0.0154`, `e_y6 = −0.0095` | `d_x1 = +1.20`, `d_x2 = −1.19`, `d_x6 = +1.49`, `d_y1 = −0.67`, `d_y2 = +0.67`, `d_y6 = −0.83` | `e_x1 = +0.0081`, `e_x2 = −0.0080`, `e_x6 = +0.0053`, `e_y1 = −0.0043`, `e_y2 = +0.0043`, `e_y6 = −0.0028` | `d_x1 = +0.92`, `d_x2 = −0.91`, `d_x6 = +1.21`, `d_y1 = −0.49`, `d_y2 = +0.49`, `d_y6 = −0.65` |
| PE | all zero (< 1e-12 C/m²) | all zero (< 1e-11 pC/N) | all zero | all zero |

**The shape of this table is itself a result, and it is a limitation.** Every `e_iJ` with
`J` diagonal is zero for beta, alpha and PE. That is not a bug. The charges are fixed point
charges on a rigid chain, so the cell dipole cannot change under a dilation at all, and the
proper constant removes the `1/V` term that the naive `dP/deps` would have shown. The only
piezoelectric channel the model has is **chain reorientation**: a field transverse to an
aligned dipole torques the chain, the chain rotates, and the lattice shears. Hence
`e_y6` for beta — an `E_y` produces `eps_xy`. Gamma is the one case with a full row, because
its chains are not aligned with a cell axis to begin with.

A direct corollary: for beta, the *polarization* direction is the worst poling direction
there is. A field along a dipole that is already aligned exerts no torque. The reported
actuator direction is therefore `y`, transverse to `P`, and a field along `P` gives a work
density smaller by 10^6 (asserted in the tests). In a real crystal the missing response
comes from the electronic polarizability, the bond-dipole change under strain, and the
chain-direction deformation — all three of which this potential does not have.

### Blocking stress, free strain, work density (E = 0.01 V/A = 100 MV/m)

| potential | | poling direction | free strain | blocking stress (MPa) | work density (kJ/m³) |
|---|---|---|---|---|---|
| illustrative | PVDF beta | ±y | `e_xy = ∓1.82e−3` | ∓5.6 | 2.56 |
| | PVDF alpha | ±y | `e_xy = ∓1.70e−4` | ∓2.0 | 0.087 |
| | PVDF gamma | (0.873, −0.488, 0) | `e_xx = +1.37e−4`, `e_yy = −1.37e−4`, `e_xy = +1.71e−4` | +3.2, −3.1, +1.9 | 0.299 |
| | PE | any | 0 | 0 | 0 |
| `pvdf-dft-fit` | PVDF beta | ±y | `e_xy = ±2.46e−4` | ±0.7 | 0.044 |
| | PVDF alpha | ±y | `e_xy = ±2.79e−4` | ±1.5 | 0.102 |
| | PVDF gamma | (0.882, −0.471, 0) | `e_xx = +1.04e−4`, `e_yy = −1.03e−4`, `e_xy = +1.37e−4` | +0.9, −0.9, +0.6 | 0.068 |
| | PE | any | 0 | 0 | 0 |

The poling direction is the leading eigenvector of `d C d^T`, whose overall sign is
arbitrary; flipping the field flips the free strain and the blocking stress together and
leaves the work density alone, hence the ±. The full elastic triangle is twice the work
density in each case. Sign convention as in section 2: the blocking stress is what the crystal pushes with, and the free strain is the
strain it reaches unloaded.

### The polyethylene null

PE's atoms do carry charges (C −0.12, H +0.06 e). What vanishes is the cell dipole, and it
vanishes for *every* configuration, not only at the minimum: each CH2 group is neutral and
the all-trans repeat is centrosymmetric. The cancellation is exact to 1.7e−16 e·A over
random cells. So every piezoelectric number comes out at machine precision zero — `e` below
1e−12 C/m², `d` below 1e−11 pC/N, work density below 1e−24 kJ/m³ — while the elastic constants are
perfectly ordinary. That makes it a real test of the whole chain (dipole, frame rotation,
both derivative routes) rather than a multiplication by zero, and it is asserted in
`tests/test_mechanics.py::test_polyethylene_has_exactly_no_piezoelectric_response`.

## 6. How converged these numbers are, and by what

Two knobs, and they do not matter equally.

**The strain step barely matters, and neither does the relaxation tolerance.** Between
`deps = 5e-4` and `4e-3`, `C_11`, `C_12`, `C_66` and `e` move by well under 1%, and the
pre-symmetrisation asymmetry of the stiffness block stays under 0.05 GPa for beta and 0.2 GPa
for PE. `pack.polish` stops at `gtol = 1e-5`, which looked like a loose tolerance to be taking
second derivatives through — the leftover gradient enters the stress difference divided by the
step — so it was checked rather than assumed: restarting L-BFGS-B from its own answer, once or
three more times, moves no constant of beta or PE by 0.002 GPa even at `deps = 5e-4`. The
suspicion was wrong, and the code does one pass.

**The nonbonded cutoff matters a lot, and it is the accuracy limit here.** The kernel's
Lennard-Jones term is energy-shifted but not force-shifted, so its force jumps at `r = rc`
— `pack.py` documents this in `_pair_energy_and_dv`. Straining the cell walks pairs across
that jump. Beta-PVDF's `b = 8.59 A` sits right on the package default `cutoff = 8.0`, and
the symptom is unmistakable: at 8 A, `C_22` still drifts 3.7 GPa between `deps = 5e-4` and
`4e-3`; at 12 A the same spread is 0.06 GPa. Beta, illustrative potential, `deps = 2e-3`:

| cutoff (A) | relaxed b | `C_11` | `C_22` | `C_12` | `C_66` | `e_y6` | `d_y6` |
|---|---|---|---|---|---|---|---|
| 8 (package default) | 8.587 | 44.42 | 37.76 | 2.25 | 3.10 | 0.0564 | 18.17 |
| 10 | 8.552 | 45.82 | 41.05 | 2.57 | 3.37 | 0.0586 | 17.38 |
| 12 | 8.553 | 45.72 | 42.09 | 2.56 | 3.36 | 0.0585 | 17.40 |
| 16 | 8.547 | 46.55 | 42.73 | 2.58 | 3.38 | 0.0587 | 17.39 |

So `C_22` at the default is ~11% low and `d_y6` ~4% high; everything is converged by 10-12 A.
The tables in section 5 are quoted at the package default, because that is the potential
every other number in this repository was computed with and quietly changing it here would
make them incomparable. Anyone wanting elastic constants for their own sake should pass
`cutoff=12.0` to `reference_from_chain`. `tests/test_mechanics.py::
test_step_dependence_of_C22_is_the_lennard_jones_cutoff` pins the diagnosis rather than the
symptom.

## 7. Cost

One full response calculation — reference relaxation, the six strained relaxations for the
3×3 stiffness and `e`, and the Newton solve for `d` along all three field axes — on one CPU
core:

| | beta | alpha | gamma | PE |
|---|---|---|---|---|
| response | 0.32 s | 1.44 s | 7.8 s | 0.56 s |
| kernel rows | 100 | 226 | 217 | 187 |
| pack + refine (once, beforehand) | 1.8 s | 2.9 s | 10.0 s | 1.7 s |

A hundred to two hundred kernel rows for a whole response. This is affordable only because of the
analytic gradient: `polish(..., gradient="analytic")` spends one kernel row per function
evaluation where the finite-difference route spent `1 + 2 n_free`, and the stress itself is
closed-form — the chain rule from the Voigt strain through `(a, b, gamma, phi1, phi2, dz, c)`
contracts the kernel's own `g_cell` and `g_c`, so the `3x3` stiffness costs six relaxations
rather than the thirty-odd energy evaluations a second difference of the energy would need.
`tests/test_mechanics.py::test_analytic_stress_matches_a_finite_difference_of_the_energy`
checks it against a finite difference of the energy to 1e−4.

## 8. Literature, for order of magnitude only

Nothing above is a prediction of experiment, for reasons the package already records: the
illustrative potential is not quantitative (DESIGN.md section 5.4), and the fitted presets
are fitted to single-chain DFT energies and forces with no elastic or piezoelectric data in
the objective at all (docs/DFT_FIT.md, docs/VALENCE_FIT.md). The comparison below is a
sanity check on order of magnitude and nothing more.

Two conventions to keep straight before reading across. The published transverse numbers for
polymer crystals are usually **Young's moduli** `E_a`, `E_b`, not stiffness constants; `1/S_11`
of the in-plane block is the comparable quantity here and is quoted below alongside `C_11`,
from which it differs by only a few per cent because `C_12` is small. And it is an **axially
clamped** modulus, so it is an upper bound on the free `E_a`.

Illustrative | `pvdf-dft-fit`, from section 5, GPa:

| | `C_11` | `E_a = 1/S_11` | `C_22` | `E_b` | `G_ab = 1/S_66` |
|---|---|---|---|---|---|
| PVDF beta | 44.4 \| 21.6 | 44.3 \| 21.4 | 37.8 \| 15.5 | 37.7 \| 15.3 | 3.10 \| 2.90 |
| PVDF alpha | 28.7 \| 13.9 | 26.3 \| 12.4 | 22.5 \| 11.1 | 20.7 \| 9.9 | 12.03 \| 5.24 |
| PVDF gamma | 29.6 \| 11.6 | 26.3 \| 10.6 | 31.6 \| 12.7 | 28.6 \| 11.7 | 9.04 \| 4.27 |
| PE | 27.5 \| 36.2 | 15.8 \| 6.3 | 30.5 \| 32.4 | 17.7 \| 6.2 | 12.75 \| 4.52 |

### 8.1 Beta-PVDF transverse stiffness

*Literature (calculated, DFT/PBE, VASP):* `C_11 = 19.83`, `C_22 = 24.50`, `C_33 = 287.33` GPa
from a stress-strain fit; `9.19`, `11.30`, `270.81` GPa from an energy-strain fit — the
factor-of-two spread between the two routes is the paper's own, not a transcription error.
Bulk modulus 11.18 GPa (PBE) / 25.25 GPa (PBE-D2).

> Pei, Y.; Zeng, X. C. "Elastic properties of poly(vinylidene fluoride) (PVDF) crystals: a
> density functional theory study." *J. Appl. Phys.* **109**, 093514 (2011).
> [10.1063/1.3574653](https://doi.org/10.1063/1.3574653)

*Against this package:* the illustrative potential's 44 / 38 GPa is 2-5x too stiff. The fitted
preset's 22 / 15 GPa lands within about 1.5x of the stress-strain DFT numbers and brackets the
energy-strain ones — a better showing than the fit deserves, since no elastic data entered its
objective. Neither is a prediction; both are the right order.

*Not found:* `C_12`, `C_13`, `C_23`, `C_44`, `C_55`, `C_66` for beta-PVDF. Pei and Zeng put the
full 6x6 matrices in AIP supplementary deposit E-JAPIAU-109-105107, which could not be
retrieved, so there is nothing to compare `C_12` or `C_66` against. The lattice-dynamics paper
that would give both elastic *and* piezoelectric constants — Tashiro, K.; Kobayashi, M.;
Tadokoro, H.; Fukada, E., *Macromolecules* **13**, 691 (1980),
[10.1021/ma60075a040](https://doi.org/10.1021/ma60075a040) — is paywalled and was not read.

### 8.2 Chain-direction modulus, and why section 4 matters

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

**This is the reason `C_33` is reported as `nan` and not as a number.** Section 4's angle path
gives 401 GPa for beta and 381 GPa for PE at `k = 105`. Those sit within about 1.4x of the
published chain moduli — close enough to pass an order-of-magnitude sanity check, close enough
to be quoted in a table, and they are 79% and 84% an invented bend constant with the
bond-stretch channel missing entirely. A plausible-looking wrong number is worse than a missing
one.

### 8.3 Transverse moduli of polyethylene

*Literature (measured):* `E_a = 3.1`, `E_b = 3.8` GPa by X-ray at 20 °C (Sakurada, Ito and
Nakamae, via Kurita's Table 3); `E_a` and `E_b` both about 6 GPa by neutron scattering at 25 °C
(Holliday, 1971); `E_a = 9`, `E_b = 8` GPa by neutron at −196 °C (Twisleton, 1982).
*Calculated:* 10.9 / 7.8 GPa (Kurita 2018, B3LYP-D, 0 K); 6.9 / 8.6 GPa (Tashiro 1978);
5.9-13.7 GPa across older lattice-dynamics work.

*Against this package:* the illustrative potential gives 15.8 / 17.7 GPa, roughly 5x the X-ray
values and 1.5-2x the stiffest calculation. The fitted preset gives 6.3 / 6.2 GPa, which sits
on top of the neutron measurement — and that agreement should be **discounted entirely**,
because the same preset squeezes PE's ab cross-section to 26.6 A² against an experimental
36.7 (section 5). A potential fitted to PVDF chains that gets PE's cell 28% too dense has
not earned a transverse modulus.

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

*Against this package:* the largest coefficient here is `d_y6 = 18.2` pC/N (illustrative) and
`2.5` pC/N (fitted) — the same order as a measured `d_31`, and **not the same coefficient**.
`d_y6` couples a transverse field to a shear; `d_31` couples the poling field to an extension
along the draw direction, and is identically zero in this model. The agreement in magnitude is
a coincidence of scale, and reading it as a match would be a category error.

*Not found:* `e_31` for PVDF in C/m² from any readable source, so the `e` table has no
literature counterpart at all. Calculated single-crystal piezoelectric constants also could not
be verified: Nakhmanson, S. M.; Buongiorno Nardelli, M.; Bernholc, J., *Phys. Rev. Lett.* **92**,
115504 (2004), [10.1103/PhysRevLett.92.115504](https://doi.org/10.1103/PhysRevLett.92.115504)
is the right paper and was not opened; a search-engine summary quoting `e_33 = −0.332` for
beta-PVDF is **not** cited here because it could not be checked. Karasawa, N.; Goddard, W. A.
III, *Macromolecules* **25**, 7268 (1992) reports calculated elastic, dielectric *and*
piezoelectric constants for nine PVDF structures and is the single most useful paper for this
section; it is paywalled and was not read.

### 8.5 The mechanism, and why the diagonal columns are empty

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

*Against this package:* this model has neither channel. The dipole per reference volume cannot
change under a dilation at all — fixed point charges on a rigid chain — which is why every
diagonal column of the proper `e` is zero, and the naive `dP/deps` that would have shown a
non-zero number there is nothing but the `1/V` derivative (section 2). The intrinsic third
needs the dipole to change with strain, which needs a polarizable or stretchable chain. What is
left is chain *reorientation*, a channel the classic decomposition does not list among its main
two. The actuator numbers in section 5 are still the physically meaningful quantity for an
actuator — free strain per applied field, obtained by direct relaxation and independent of the
proper/improper question — they are just the free strain of a mechanism that is not PVDF's
dominant one.


## 9. What is not computed

* `C_33`, and therefore the chain-direction modulus, Young's modulus along any direction
  with a `z` component, the bulk modulus, and Poisson ratios involving the chain axis —
  section 4.
* `C_44`, `C_55`, `C_45`, and `C_14`, `C_15`, `C_24`, `C_25`, `C_34`, `C_35`, `C_46`,
  `C_56`: the cell cannot express `eps_yz` or `eps_xz` — section 1.
* `C_13`, `C_23`, `C_36`: these need a strained-`c` state to be differentiated against, so
  they inherit section 4's problem. `dsigma_3/deps_J` is available from the same relaxations
  and is deliberately not reported as an elastic constant, because `sigma_3` at a fixed
  chain shape is not an axial stress — the only way a rigid chain's `c` can grow is for the
  repeats to slide apart across the periodic boundary, and the nonbonded sum excludes
  exactly the bonded pairs that cross it, so nothing resists.
* The full compliance tensor, hence `d` as an unclamped coefficient. `S` here is the inverse
  of the 3×3 in-plane block and `d` is the axially clamped strain coefficient.
* Any response that needs the chain conformation to relax at fixed strain — the torsions and
  bond angles are frozen through every strain state, so these are rigid-chain constants. A
  conformational relaxation would change `c` and break the `eps_zz` clamp; doing it properly
  needs `c` as a constrained variable, which needs section 4 solved first.
* Electronic polarizability, depolarisation, and any field-induced change in the charges:
  the model is fixed point charges, so the dielectric constant is 1 by construction and the
  piezoelectric response is purely orientational.
* Temperature. Everything is a zero-kelvin second derivative; no phonons, no thermal
  expansion, no pyroelectric coefficient.
* Electrostriction. The response reported is strictly linear in the field; DESIGN.md section
  5.5's observation that a 0.5 V/A transverse field opens beta's `b` from 8.61 to 8.85 A is a
  quadratic effect and is not in any `d` here.
* Domain switching. A field opposing `P` finds the 180-degree-rotated cell at equal energy
  (DESIGN.md section 5.5); the linear response about one minimum says nothing about the path
  between them or about the coercive field.
