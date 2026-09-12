# Internal strain: our Jacobian against the provider's, and what it does to the pendant hypothesis

`examples/internal_strain_jacobian.py` produces every number here; the reader is
`polyfind.born.load_internal_strain`. Nothing below is fitted to anything, and no physics or
parameter was changed.

## The verdict, phrased as the provider asked

**The two models' pendant kinematics differ exactly where the hypothesis said they would,
and the difference carries no polar dipole through the Born charges.** Under a strain along
the long lateral axis the provider's C–F bonds stretch by 0.24 Å and turn by 31° per unit
strain while ours stay rigid to 1e−13; but the internal-strain dipole that their displacement
pattern produces in our Born-consistent charge model is +0.33 C/m², the same +0.335 our own
rigid chain produces, and through their DFPT tensor it is +0.32. Read through Born charges of
either origin, the *total* atomic motion of their relaxed-ion geometries — affine plus
internal — gives between −0.004 and +0.06 C/m² of polar dipole per unit strain on every
normal axis, which is what a rigid chain gives identically. So the hypothesis that the
missing 0.4–0.6 C/m² is internal strain that rigid pendants cannot supply is **not supported
by this record**: it is a difference between the models' kinematics, not between their
piezoelectric responses. What can still differ is the electronic part of the clamped-ion
response, which neither dataset contains (the provider lists a same-Hamiltonian clamped-ion
response among its open gates), or the target itself, which rests on two Berry slopes that
both failed their numerical gates.

This establishes a difference between the two models' pendant kinematics. It does not
confirm real pendant flexibility, and cannot until the provider's reference gates clear:
position and force tolerance, basis, k-grid and geometry sensitivity, and the anomalous
transverse stiffness are all unisolated in the record, which is provisional and not an
accepted Jacobian.

## The data, and its caveats carried whole

Sarco `ccdf377`, `materials/gpu_bundle/periodic_reference/beta_pvdf/internal_strain_geometry.json`
and `INTERNAL_STRAIN_GEOMETRY_REPORT.md`, read only. Thirteen completed CP2K fixed-cell
relaxed-ion geometries (PBE-D3(BJ), `cp2k_2026.2_pbe_d3bj_relaxed_ion_fixed_cell_strain_v1`):
the zero cell and ±1 %, ±2 % normal strain along each axis, every source checked by the
provider against its geometry, output and force receipts (zero-cell maximum force 3.2e−4
eV/Å). The record carries `quantitatively_valid: false` and four open gates: atomic force
tolerance and numerical displacement sensitivity; CP2K basis, cutoff, k-grid and geometry
convergence; the transverse mechanical anomaly; and a same-Hamiltonian validated Born and
clamped-ion response, of which none is supplied. Neither failed Berry slope enters it. Their
zero geometry is the `beta_pvdf_equilibrium_zero-FINAL-1_45.xyz` the Born run used (SHA-256
`3b3a7714…`, checked against the receipt).

Their construction, which ours copies line for line: `u = r_relaxed − r_zero F`, evaluated as
the wrapped difference of fractional coordinates carried into the strained cell (a
displacement over 0.25 of a cell is refused as an ambiguous branch); the equal-atom mean
translation projected out and recorded; `J = (u(+h) − u(−h)) / 2h` in Å per unit engineering
strain. For each H or F, the bond vector from its nearest carbon at each state, its central
derivative *including the affine term*, the rate along the reference bond (length) and the
rate across it over the length (direction, degrees per unit strain). Their reason for the
bond vectors is the one that makes the test decisive: a rigid pendant cancels the affine
deformation of its bond, so nonaffine displacements alone cannot show flexibility. The
1 %-to-2 % change is `|J(2 %) − J(1 %)| / |J(1 %)|` in the Frobenius norm.

We fit nothing to it; we do not contract it with the DFPT Born tensor as though a validated
same-Hamiltonian result existed; the one contraction below is labelled for what it is.

## Frames and atom order

Their axes: Cartesian x, y, z are the row lattice vectors A, B, C; y is polar, z the chain, so
x is the long lateral axis; zero cell 8.358275 × 4.731416 × 2.580155 Å. Ours (the packer's):
x polar, y long, z chain. `polyfind.born` already canonicalises to the provider's `(a, b, c)` =
(long, polar, chain) for the Born comparison; here the same map is written as a proper rotation
by 90° about the chain axis, `(x, y, z)_provider = (−y, x, z)_packer`, with the polar sign fixed
as `born.canonical` fixes it — fluorine at +y of its carbon, as in their cell (the sign came
out +1; our P is −0.143 C/m² along −y, theirs likewise has the C→F direction as +y). Voigt:
our `eps_xx` (polar) is their `yy`, our `eps_yy` (long) is their `xx`, `zz` is `zz`.

Our reference is the `d_33`/`d_31` state of `docs/ELECTROMECHANICS.md` 5.9:
`pvdf-dft-valence-flux-born` with induced dipoles, cell and line-group shape relaxed to zero
stress on the deformable path (a = 4.5635 Å polar, b = 8.5029 long, γ = 90.00°, c = 2.5494;
setting angles 180°, chain offset dz = 0; shape parameter +0.17°). In their frame our lattice is
8.503 × 4.563 × 2.549 Å against their 8.358 × 4.731 × 2.580, V₀ 98.92 against 102.04 Å³. Both
cells have the same atom order per chain (C, H, H, C, F, F) and the same C-centred placement of
the second chain at (½, ½, 0); matching by fractional coordinates pairs every atom to within
0.064 Å (the H's and F's differ by the mirror, which is a symmetry of both cells).

Their atoms, our packer index, their position and ours in their frame at their fractional
coordinates (Å):

| their | element | ours | theirs (x, y, z) | ours |
|---:|---|---:|---|---|
| 1 | C(H2) | 0 | (0.000, −0.052, 0.000) | (0.000, −0.052, 0.000) |
| 2 | H | 1 | (+0.887, −0.698, 0.000) | (+0.867, −0.716, 0.000) |
| 3 | H | 2 | (−0.887, −0.698, 0.000) | (−0.867, −0.716, 0.000) |
| 4 | C(F2) | 3 | (0.000, +0.777, 1.290) | (0.000, +0.822, 1.290) |
| 5 | F | 5 | (+1.112, +1.627, 1.290) | (+1.060, +1.664, 1.290) |
| 6 | F | 4 | (−1.112, +1.627, 1.290) | (−1.060, +1.664, 1.290) |
| 7–12 | second chain | 6–11 | offset (4.179, 2.366, 0) | the same |

Ours is the deformable path: at every strain `deformable_state` relaxes the setting angles,
the chain offset and the line-group shape (one parameter for the beta zigzag, the backbone
angle) at fixed strain against the fitted valence terms; bond lengths, pendant angles and
pendant bond lengths are rigid, as they are everywhere in the shipped model. Every strained
state converged (`c` error ≤ 2e−16, no frame rotation), and the Born tensor at the reference
has acoustic sum 1.4e−8 e.

## 1. The Jacobian, atom by atom

Å per unit strain, their frame and order, both at the 1 % central difference. The "rigid"
column is what a cell of rigid chains pinned to their lattice sites gives under a transverse
strain, `−(r − r̄_chain)` along the strain axis and zero across it; `J − rigid` is the motion
within a chain, and it is exactly zero for us on both transverse axes. No such reference exists
for `zz`, which a rigid chain cannot take up at all. The second chain repeats the first to
0.002 Å per unit strain in their record and exactly in ours, so one chain is shown.

**Strain xx (long axis).** rms theirs 0.596 (their reported 0.596431), of which the rigid-chain
value on their geometry would be 0.821 and the within-chain remainder is 0.428; ours 0.804 =
rigid; rms difference 0.411.

| atom | theirs J | ours J | theirs J − rigid |
|---|---|---|---|
| C1 (CH2) | (0.000, +0.142, 0.000) | (0, 0, 0) | (0.000, +0.142, 0.000) |
| H2 | (−0.851, +0.142, +0.002) | (−0.882, 0, 0) | (+0.036, +0.142, +0.002) |
| H3 | (+0.852, +0.142, +0.002) | (+0.882, 0, 0) | (−0.035, +0.142, +0.002) |
| C4 (CF2) | (0.000, +0.157, 0.000) | (0, 0, 0) | (0.000, +0.157, 0.000) |
| F5 | (−0.463, −0.290, −0.001) | (−1.078, 0, 0) | (+0.649, −0.290, −0.001) |
| F6 | (+0.462, −0.292, −0.001) | (+1.078, 0, 0) | (−0.650, −0.292, −0.001) |

Where their atoms go: the fluorines follow 58 % of the affine x-stretch instead of none
(+0.65 Å per unit strain relative to rigid, each) and at the same time drop 0.29 Å per unit
strain towards the chain axis in y, while every other atom of the chain moves +0.14 to +0.16
in y; the CH2 and CF2 group centroids move ±0.142 in y, so the two groups approach each other
by 0.28 Å per unit strain along the polar axis. The carbons themselves stay 0.015 apart from
rigid, so the backbone barely changes: what moves is the F–C–F pendant, opening in x and
flattening in y. Ours: nothing within a chain.

**Strain yy (polar axis).** rms theirs 0.952 (rigid-chain value 0.980, remainder 0.034); ours
0.969 = rigid; rms difference 0.030.

| atom | theirs J | ours J | theirs J − rigid |
|---|---|---|---|
| C1 (CH2) | (0.000, +0.454, −0.001) | (0, +0.479, 0) | (0.000, −0.028, −0.001) |
| H2 | (+0.024, +1.098, −0.005) | (0, +1.119, 0) | (+0.024, −0.030, −0.005) |
| H3 | (−0.023, +1.098, −0.005) | (0, +1.119, 0) | (−0.023, −0.030, −0.005) |
| C4 (CF2) | (0.000, −0.319, −0.001) | (0, −0.364, 0) | (0.000, +0.028, −0.001) |
| F5 | (−0.015, −1.165, +0.007) | (0, −1.176, 0) | (−0.015, +0.031, +0.007) |
| F6 | (+0.015, −1.167, +0.007) | (0, −1.176, 0) | (+0.015, +0.030, +0.007) |

Along the polar axis their chain is, to 0.03 Å per unit strain, the rigid body ours is: the two
groups approach by 0.06 and the pendants tilt by 0.02. The 0.030 rms difference between the
models here is mostly the 3 % difference in the cell, not a difference in kinematics.

**Strain zz (chain axis).** rms theirs 0.114, ours 0.965, difference 0.857 — the largest
disagreement of the three, and not about pendants.

| atom | theirs J | ours J |
|---|---|---|
| C1 (CH2) | (0.000, +0.167, −0.005) | (0, +0.965, 0) |
| H2, H3 | (∓0.036, +0.079, −0.009) | (0, +0.965, 0) |
| C4 (CF2) | (0.000, −0.109, +0.001) | (0, −0.965, 0) |
| F5, F6 | (±0.001, −0.108, +0.011) | (0, −0.965, 0) |

Our chain can take up an axial strain only by opening the backbone angle, which moves the CH2
and CF2 groups apart along the polar axis by 1.93 Å per unit strain — the ∓1.99 of
`docs/ELECTROMECHANICS.md` 5.9 at the earlier structure. Their carbons move 0.28 apart, a
seventh of ours; on this record the rest of the strain is taken up by the C(H2)–C(F2) bond
itself, the channel our rigid bonds do not have. Their H's lag their carbon by 0.09 (the C–H
bond turns, below). This branch is also the one their own 26.6 % amplitude sensitivity sits
on, so the seventh is the least secure number here.

## 2. The bond vectors, affine term included

`d(r_pendant − r_carbon)/d eps`, Å per unit strain, total. A rigid pendant on a rigid chain
gives exactly zero here under every strain, because the nonaffine motion cancels the affine
term; a pendant that simply followed the lattice would show its reference vector's component
along the strain axis (the "affine" column) as its length rate.

**Strain xx**, per bond (chain 1; chain 2 identical to 0.001):

| bond | v₀ theirs (Å) | affine | theirs dv | length | direction | ours dv, length, direction |
|---|---|---:|---|---:|---:|---|
| C1–H2 | (+0.887, −0.646, 0) | +0.887 | (+0.035, −0.001, +0.002) | +0.029 | 1.07° | 0, 2e−14, 3e−13 |
| C1–H3 | (−0.887, −0.646, 0) | −0.887 | (−0.035, −0.001, +0.002) | +0.029 | 1.07° | 0, 9e−15, 3e−13 |
| C4–F5 | (+1.112, +0.850, 0) | +1.112 | (+0.649, −0.447, −0.001) | +0.244 | 30.7° | 0, 2e−14, 6e−13 |
| C4–F6 | (−1.112, +0.850, 0) | −1.112 | (−0.650, −0.449, −0.001) | +0.244 | 30.8° | 0, 2e−14, 6e−13 |

Summary over the three axes, mean over the four bonds of each type (theirs; ours at the
numerical floor in every entry):

| strain | C–F length rate (Å/strain, % of bond) | C–F direction (°/strain) | C–H length rate | C–H direction |
|---|---:|---:|---:|---:|
| xx (long) | +0.244 (17.4 %) | 30.7 | +0.029 (2.6 %) | 1.05 |
| yy (polar) | −0.010 (0.7 %) | 0.57 | +0.020 (1.8 %) | 0.66 |
| zz (chain) | +0.001 (0.1 %) | 0.44 | +0.023 (2.1 %) | 4.77 |
| ours, any | ≤ 4e−14 | ≤ 2e−9 | ≤ 3e−14 | ≤ 2e−9 |

So the answer to the specific question is: under the long-axis strain their C–F bonds stretch
and turn, by a great deal per unit strain (a 1 % strain lengthens the bond by 0.0024 Å and turns
it 0.3°), and ours do not move at all; under the polar strain theirs barely move either
(0.7 % and 0.6°), and under the chain strain the C–F bond is rigid in both models while their
C–H bond turns 4.8° per unit strain as the H's lag their carbon. The pendant that flexes is the
CF2 group under the long-axis strain, and it flexes by opening: both F's move outwards along x
and inwards along y, which is the F–C–F angle widening with the bond lengthening. That is also
the axis on which the provider reports an anomalous transverse stiffness; whether the two are
related is theirs to isolate.

## 3. The contraction — cross-Hamiltonian diagnostic, not a result

`(1/V₀) Σ_k Z_k · J_k` in C/m²: the dipole per reference volume that a displacement pattern
`J` produces through Born tensors `Z`, both in their frame and atom order, `V₀` the cell the
kinematics belong to. Their `J` is the provisional CP2K geometry; our `Z` is the
Born-consistent model's own tensor (`polyfind.born.born_charges`, acoustic sum 1.4e−8 e). It
is not a piezoelectric coefficient and not a same-Hamiltonian quantity; it answers one
question only: if our atoms moved as theirs do, what would our charges make of it. For scale,
our proper coefficient splits exactly as `e = e_clamped + (1/V₀) Σ Z · J` — `e_clamped` being
every atom displaced affinely with the cell, charges fluxed and dipoles re-solved, nothing
relaxed — and the closure of that identity on the same two states is 4e−5 C/m² or better on
every axis.

Polar (y) component; the x and z components are below 2e−3 everywhere:

| strain | our Z · their J (1 %) | (2 %) | our Z · our J | our e_clamped | our proper e, total | their Z(asr) · their J (1 %) |
|---|---:|---:|---:|---:|---:|---:|
| xx (long) | **+0.331** | +0.333 | +0.335 | −0.296 | +0.039 | +0.318 |
| yy (polar) | +0.401 | +0.404 | +0.395 | −0.299 | +0.096 | +0.495 |
| zz (chain) | −0.004 | +0.009 | +0.014 | −0.011 | +0.003 | −0.002 |

The last column contracts their own ASR-corrected DFPT tensor (QE, the Born run) with their
CP2K geometry — provisional on both sides and not same-Hamiltonian — and is shown only as a
bound on how much the choice of charges matters; the raw tensor gives the same to 1e−4.

By atom type, the polar component of the xx row: their `J` through our `Z` is C(H2) −0.008,
H +0.018, C(F2) +0.064, F +0.257; our rigid `J` through our `Z` is H +0.010, F +0.325; their
`J` through their `Z` is C(H2) −0.009, H −0.001, C(F2) +0.072, F +0.256. The F's move less along
x in their record than in ours (+0.65 against +1.08 nonaffine), which through `Z_yx` gives less
dipole, and at the same time move −0.29 in y and their carbon +0.16, which through `Z_yy` gives
it back: the pendant opening is polar-neutral to 0.004 C/m² in our charges and to 0.013 in
theirs.

**The Born-charge reading of the total motion.** Adding the affine term, `(1/V₀) Σ Z · (E_a r +
J)`, gives the ionic part of a proper coefficient in a Born-charge picture. For rigid chains
pinned to their sites it is identically zero, chain by chain, by the acoustic sum rule; the
entire +0.039 of our transverse coefficient is the induced dipoles' response to the strained
lattice, the part of `e_clamped` that Born charges do not describe.

| strain | our Z, our motion | our Z, their motion | their Z(asr), their motion |
|---|---:|---:|---:|
| xx | −0.335 + 0.335 = 0.000 | −0.335 + 0.331 = −0.004 | −0.259 + 0.318 = +0.059 |
| yy | −0.395 + 0.395 = 0.000 | −0.399 + 0.401 = +0.002 | −0.494 + 0.495 = +0.001 |
| zz | 0 + 0.014 = +0.014 | 0 − 0.004 = −0.004 | 0 − 0.002 = −0.002 |

**Reading.** The diagnostic was set up to say: near 0.4–0.6, the charges are consistent with the
motion being what is missing; away from it, the charges are implicated too. It lands at +0.33,
outside the band — but the number that matters is that it lands within 0.004 of what our own
rigid kinematics give, and within 0.06 of zero once the affine term is added, with either
tensor. Neither the charges nor the kinematics, through Born charges, supply the missing half
a C/m²: a rigid chain and the DFT relaxed-ion geometry agree on the ionic part of `e_y,xx` to
0.06 C/m². If a proper transverse coefficient of 0.4–0.6 is real, it is electronic — the part of
the clamped-ion response that is not the Born response to atomic motion — and no dataset on
either side has it. The provider's open gate for a same-Hamiltonian clamped-ion response is
exactly the datum, and the other possibility is that the target is not real: both Berry slopes
it was read from failed their gates.

## 4. One percent against two

`|J(2 %) − J(1 %)| / |J(1 %)|`, Frobenius, and the largest single component change in Å per
unit strain:

| strain | ours | theirs (their report) | theirs, largest component |
|---|---:|---:|---:|
| xx | 0.000 % (4e−14) | 6.412 % | 0.058 (on F, the x motion) |
| yy | 0.000 % (9e−15) | 0.867 % | 0.008 |
| zz | 0.113 % (1.1e−3) | 26.563 % | 0.035 |

Ours is linear to machine precision on both transverse axes, because a rigid chain's
nonaffine displacement is exactly linear in the strain, and to 0.11 % on the chain axis, where
the backbone angle reaches +1.754° / −1.714° at ±1 % and +3.551° / −3.392° at ±2 % (the
asymmetry is the whole of the nonlinearity). Their 26.6 % on the chain axis sits on a branch
whose displacements are small (rms 0.114 → 0.136 Å per unit strain between the two estimates)
against a zero-cell force tolerance of 3e−4 eV/Å, which is their first open gate; their 6.4 %
on the long axis is the pendant opening changing by 0.06 Å per unit strain between amplitudes
(within-chain remainder rms 0.428 → 0.390). We report ours as linear and theirs as their
diagnostic says, and draw nothing from the comparison beyond that.

## What this settles, and what it does not

* The kinematic difference between the models is real and located: under a long-axis strain
  the DFT relaxed-ion CF2 pendant opens (0.24 Å and 31° per unit strain in the C–F bond) and
  the chain's motion is 0.43 Å per unit strain rms away from a rigid body's; ours is a rigid
  body exactly. Under the polar
  strain both chains are rigid bodies to 0.03 Å per unit strain. Under the chain strain theirs
  stretches the backbone bond where ours can only open the angle.
* That difference does not supply the missing transverse response through Born charges of
  either origin: the ionic part of `e_y,xx` on their kinematics is between −0.004 and +0.06
  C/m². The internal-strain hypothesis, as stated on 2026-09-12 in `docs/BENCHMARK.md`, is
  not supported by the first data that could test it. What remains open is the electronic
  clamped-ion term or the validity of the 0.4–0.6 target.
* None of this confirms real pendant flexibility or a real coefficient. The record is
  provisional on the provider's own terms; every number above inherits its open gates; the
  contraction is cross-Hamiltonian and labelled so; nothing was fitted.

Reproduce with `python examples/internal_strain_jacobian.py --data <internal_strain_geometry.json>
--geometry <beta_pvdf_equilibrium_zero-FINAL-1_45.xyz> --json out.json` (about 30 s; the
`--born` default is the sibling checkout's `born_results.json`, found beside `--geometry` when
the checkout is elsewhere).
