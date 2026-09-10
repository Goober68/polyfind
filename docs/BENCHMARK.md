# Speed against the existing method

The project's target is a solver that meets field and strain measurement needs at
10 to 100 times the speed of existing methods. That claim needs a baseline that
is recorded rather than assumed, and one exists: the sibling project
`sarco/materials/gpu_bundle` kept per-run timings for the calculations this
package is meant to displace. All baseline numbers below are read from its own
result files, with the file named, so they can be checked. All polyfind numbers
were measured on the same machine (6 physical cores, 12 logical).

Nothing here is a like-for-like ratio yet, and the reason is set out under
"What is actually comparable". Read that before quoting a number.

## The baseline, from its own recorded runs

### Machine-learned potential tier

`results/field_neighborhood_refined/<system>/strain<s>_field<E>.json`, a grid of
two strain values by three field values around one already-known structure.

| System | Points | Relaxation steps | Recorded seconds |
|---|---|---|---|
| PVDF | 6 | 1,855 | 87.0 |
| CNEPO | 6 | 2,070 | 107.6 |
| AN | 6 | 2,565 | 125.4 |
| VDCN | 6 | 2,719 | 142.8 |
| total | 24 | 9,209 | 462.9 |

An earlier, looser-tolerance pass over the same grid
(`results/field_neighborhood/`) cost 176.2 s for 25 points, and several of its
field points converged in one to three steps taking about 0.11 s, meaning the
field moved nothing the tolerance could resolve.

### Density-functional tier

`results/dft_field/*.json`, PBE0 with dispersion and an applied field, one single
point each.

| Calculation | Recorded seconds |
|---|---|
| PVDF, def2-SVP, three field values | 772, 895, 1001 |
| PVDF, def2-TZVP, three field values | 2470, 1866, 1711 |
| CNEPO, def2-SVP, three field values | 946, 1127, 1030 |

So one field triple costs about 2,700 s at the smaller basis and about 6,000 s at
the larger. A strain derivative needs at least two strain states, so roughly
5,400 s per chemistry at def2-SVP, before any search over structures. One run is
recorded as having died out of memory
(`pvdf_def2-svp_field-65.attempt-oom.json`).

### Structure search tier

`results/free_angle/seq_*_11_*/result.json` holds supercell relaxations at 2,400
steps for CNEPO and VDCN, with `continuation` and `recovery` directories beside
them. A 2,400-step run that needs a continuation is a step budget exhausted
rather than a converged answer, which is the convergence risk a local optimiser
over a combined discrete-and-continuous space carries.

## polyfind, measured

| Stage | Seconds |
|---|---|
| Fit the rotational-isomeric-state model from a potential | 0.6 |
| Enumerate periodic chain conformations (120 distinct, exhaustive) | 4.4 |
| Exhaustive screen and refine, all candidates, cold cache | 77 |
| Same, warm interaction-table cache | 25 |

Fitting the potential itself against reference data is a separate one-off:
roughly 260 s over 467 frames, amortised across every later run for that
chemistry class.

## What is actually comparable

Three traps, and the ratio differs by an order of magnitude depending on which
comparison is meant.

**The baseline's 87 to 143 s per chemistry does not include finding the
structure.** It relaxes around a structure already supplied. polyfind's 25 to
77 s includes enumerating 120 candidate conformations and screening the packing
of each exhaustively. Comparing those two numbers directly credits polyfind for
work the baseline never did, and also credits the baseline with an input polyfind
derives. The structure-search tier above is the closer analogue, and it is the
one whose runs exhausted their step budgets.

**The grids differ in what they can yield.** Two strain values by three field
values supports a single finite-difference derivative at one strain. An elastic
tensor and a piezoelectric tensor need considerably more, so a fair ratio has to
compare cost per derivative obtained, not cost per run.

**The tiers are different physics.** The machine-learned tier is fast and is
itself validated against the density-functional tier in the same project.
polyfind is a third tier below both, and its accuracy is bounded by the potential
it is fitted to. A speed ratio quoted without the accuracy it was obtained at is
meaningless; see `DESIGN.md` sections 5.9 and 5.10 for what the fitted potential
does and does not get right.

## The observation that favours this approach most

It is not the raw timing. In the baseline's own refined grid, applying a field of
65 V/µm to PVDF moves atoms by at most 0.019 Å, changes torsions by at most
0.13°, and leaves a residual force of 1e-4 eV/Å after 184 relaxation steps. In
the looser pass, the same field points converged in one to three steps, which
means the response was below what that tolerance could see at all.

The quantity being measured is therefore tiny compared with the numerical
machinery used to extract it, and it is extracted as a *difference* between two
expensive relaxations. That is the regime where an analytic derivative wins for a
reason that has nothing to do with processor speed: polyfind differentiates the
energy with respect to field and strain directly, so the response is a derivative
evaluated once rather than a small difference resolved between two converged
minimisations. The same argument already paid off inside this package, where
replacing finite-difference gradients with analytic ones cut refinement from
181 s to 26 s and incidentally exposed a non-smoothness in the cutoff that the
difference-based path had been walking over unnoticed.

## Status of the headline claim

| Comparison | Status |
|---|---|
| Structure search against supercell relaxation | favourable, and the baseline's runs did not converge; not yet quantified as a ratio |
| Field and strain response against the machine-learned grid | pending: needs polyfind's response calculation, in progress |
| Field response against the density-functional tier | roughly 2,700 to 6,000 s per field triple to beat, per chemistry |
| Accuracy at which any ratio holds | held-out energies 1.36 kcal/mol, one of three acceptance tests passing |

The 10 to 100 times target is not yet demonstrated and should not be quoted until
the middle row is filled in with a measurement at a stated accuracy.
