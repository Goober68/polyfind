# Reference data request: bulk dipole response to strain

## Why this is the binding constraint

`polyfind` now reproduces both experimental piezoelectric **signs** for beta-PVDF
(d33 negative, d31 positive) where before d31 had the wrong sign, and it does so
out of sample, since the fit never saw a measured piezoelectric constant. The
**magnitudes** are 5x and 9x short and should not be trusted.

The limit is no longer the search, the speed, or the form of the potential. It is
the calibration data. The charge-flux coefficient that produces the response is
fitted against `sarco/materials/gpu_bundle/results/field_neighborhood_refined`,
which is exploratory GFN2-xTB on a **finite two-chain pair in vacuum with terminal
backbone atoms pinned**, at two strain values, with strain measured against a
fixed seed span rather than a stress-free lattice. Its own `interpretation` field
says as much.

Symptoms of that data being the limit, all measured (`docs/BENCHMARK.md`):

* the fit's R-squared is 0.39;
* a bond-stretch channel fits the same data better (0.86) and a plain charge
  rescale better still (0.75), and neither is a usable mechanism here;
* adding a bond channel alongside the angle channel changes the coefficient
  sixfold and **flips both signs back**;
* leave-one-chemistry-out spans the coefficient over a factor of eleven.

A bulk reference at the level the energy model is fitted to would settle all four.

## What is needed

**A periodic-crystal polarization response to strain.** Concretely, for
beta-PVDF:

1. **Periodic bulk cell**, not a finite cluster. The beta cell is about twelve
   atoms, so this is small.
2. **Polarization by Berry phase** (or Wannier centres). This is not optional: in
   a periodic system the dipole of a cell is not well defined, and a naive sum of
   charge times position depends on the choice of origin and cell. Take
   differences along a continuous strain path and watch for jumps by the
   polarization quantum `eR/V`; a jump misread as a response is the classic way
   this calculation goes wrong.
3. **Several strain points per component, not two.** Suggested
   `-0.02, -0.01, 0, +0.01, +0.02` for each of `eps_zz` (along the chain),
   `eps_xx` and `eps_yy`. Five points gives a linearity check and a curvature
   estimate; two gives a secant and no way to know it is one.
4. **Ions relaxed at each fixed strain**, with the cell held at the strained
   shape. The relaxed-ion response is the physical piezoelectric one. If the
   clamped-ion values come free, they are useful too, because the difference
   between them is exactly the internal-strain contribution that the charge-flux
   model is trying to represent.
5. **Report, per strain point**: the strain tensor, the Berry-phase polarization
   vector, the total energy, the relaxed atomic positions, and the residual
   forces and stress.
6. **Level of theory: PBE-D3**, to match the energy training set
   (`trainset_v1_566.xyz`). Matching matters more than absolute quality here: a
   flux fitted against one functional and an energy model fitted against another
   will not compose. If a hybrid is affordable, one point at PBE0-D3 as a spot
   check on the functional dependence would be valuable.

That is fifteen small periodic calculations, plus relaxations. On the hardware
described it should be hours, not days.

## Worth more than the above, if it is cheap to add

**Born effective charge tensors** `Z*_ij` per atom, from density-functional
perturbation theory at zero strain. These are the exact physical content of
charge flux: the difference between `Z*` and the static charge is precisely what
a fixed-charge model is missing, so they calibrate the mechanism directly instead
of through a fitted response. Most periodic codes produce them alongside the
piezoelectric tensor from the same run.

If the piezoelectric tensor itself comes out of that run, it is the direct answer
and the fit becomes a validation rather than a calibration.

## Controls that make the result trustworthy

* **Polyethylene, same protocol.** Its response must be exactly zero by symmetry.
  If it is not, something is wrong with the polarization branch tracking rather
  than with the polymer.
* **alpha-PVDF**, whose antipolar structure should give a much smaller response
  than beta, and **gamma**, which is polar but weaker. Getting the ordering right
  is a cheap sanity check on the whole protocol.
* Confirm the reference structure is **stress-free** before straining it. The
  existing data's strain is defined against a pinned span, which is part of why
  it is hard to use.

## Cross-check that costs nothing

Published work already computes polarization and Born effective charges for
beta-PVDF, including a comparison of the planar and alternately-deflected
structures (`docs/REFERENCES.md`, section 5.8's sources). Comparing any new
calculation against those before trusting it is worth the ten minutes.

## What we will do with it

Refit the charge-flux coefficients against `dPolarization/dStrain` for a bulk
crystal at the same level as the energy model, then recompute d33 and d31. The
test is already written and is out of sample: the fit does not see a measured
piezoelectric constant, so agreement with the experimental values, in sign and
then in magnitude, is a genuine prediction. If the magnitudes remain short after
this, the deficiency is in the model rather than the data, which is a different
and equally useful thing to learn.

## Cheaper fallbacks, in order

1. **Same GFN2 protocol but periodic and stress-free**, with five strain points.
   Keeps the level, removes the finite-size and pinning problems, and would show
   how much of the current scatter those cause.
2. **PBE-D3 at three strain points on beta only.** Enough to fit one coefficient
   against a bulk reference, which is the single most valuable number.
3. **Born effective charges at zero strain only.** No strain sweep, but it
   calibrates the mechanism directly and is one calculation.

Any of these is a large improvement on what the coefficient currently rests on.
