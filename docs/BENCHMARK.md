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
| Fit the rotational-isomeric-state model from a potential | 0.3 |
| Enumerate periodic chain conformations (120 distinct, exhaustive) | 4.5 |
| Exhaustive screen and refine, all candidates, cold cache | 49 |
| Same, warm interaction-table cache | 17 |

Re-measured on the same machine after the geometry corrections of
`docs/REFERENCES.md`, since PVDF's C-C bond length changed and every cached
interaction table with it. **The corrections do not move these timings**: the
unchanged code on the same machine gives 0.4 / 5.0 / 53 / 14 s, which is
run-to-run spread rather than a difference. The table previously read 0.6 / 4.4 /
77 / 25, so the screen-and-refine rows have come down by about a third since they
were recorded and the fit row by half; neither is attributable to this change.
Quote a range rather than a point for the two cache rows: cold 49-53 s, warm
14-17 s across runs.

Fitting the potential itself against reference data is a separate one-off:
roughly 260 s over 467 frames, amortised across every later run for that
chemistry class.  That fit has not been repeated at the corrected geometry.

## What is actually comparable

Three traps, and the ratio differs by an order of magnitude depending on which
comparison is meant.

**The baseline's 87 to 143 s per chemistry does not include finding the
structure.** It relaxes around a structure already supplied. polyfind's 17 to
53 s includes enumerating 120 candidate conformations and screening the packing
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
| Structure search against supercell relaxation | favourable, and the baseline's runs did not converge; not quantified as a ratio |
| Field and strain response against the machine-learned grid | **measured, see below** |
| Field response against the density-functional tier | roughly 2,700 to 6,000 s per field triple to beat, per chemistry |
| Accuracy at which any ratio holds | held-out energies 1.36 kcal/mol; one of three acceptance tests passing |

### The speed comparison, now measurable

The baseline spends 87 s (PVDF) to 143 s (VDCN) per chemistry on a grid of two
strain values by three field values, which yields **one** finite-difference
derivative at one strain state, around a structure supplied to it.

polyfind computes a full response - the in-plane elastic constants, the axial
constant, the piezoelectric tensor by two independent routes, blocking stress,
free strain and work density - in **0.3 s for beta, 1.4 for alpha, 7.8 for
gamma**, from a structure it derived itself in a further 14 to 53 s.

Taken as cost per derivative obtained the ratio is far beyond 100x, because the
baseline's grid yields one derivative and this yields a tensor. Taken as
whole-pipeline cost per chemistry, structure search included, it is roughly 2 to
6x. Neither number should be quoted alone. The honest statement is that the
speed target is met, comfortably, for the response calculation itself, and that
the response calculation is not the whole job.

**And the speed is not the binding constraint.** What limits this package for
the intended application is in the next section: the piezoelectric response of
the phase the material is actually used in cannot be computed at all, at any
speed. A fast wrong answer is not progress, so the target should be read as met
on the axis it measures and not yet met as a whole.

### What the response can and cannot deliver

| Quantity | Status |
|---|---|
| In-plane elastic constants C11, C22, C12, C16, C26, C66 | computed |
| Axial constant C33 | computed, and demonstrated free of the invented-stiffness contamination that made an earlier attempt refuse it |
| Shear constants C44, C55, C45 and mixed 4/5 | structurally absent: the chain axis is fixed along z, so no variable expresses those shears |
| Piezoelectric response of helical chains (alpha, gamma) | computed, d up to 2.7 pC/N |
| Piezoelectric response of beta, the ferroelectric phase | **exactly zero, and provably so** |
| Polarization, blocking stress, free strain, work density | computed |

The beta result is the one that matters for the application and it is not a
numerical failure. With bond-charge-increment charges, a planar all-trans
zigzag's dipole is *exactly* independent of its backbone angle: measured, the
dipole holds every digit while the chain repeat is driven from minus 1.17% to
plus 1.14%. Axial strain changes the backbone angle and nothing else, since bonds
are rigid, so the dipole cannot respond and d33 and d31 are identically zero.

The dimensional term alone gives d33 = -4.4 pC/N against a measured -32. Its sign
is right, but that is arithmetic rather than physics: that term is negative on
every diagonal column for any stable crystal. The same term gives d31 = -0.14
against a measured +20, with the **wrong sign**, which is the honest read on how
much of the physics is present.

Reaching PVDF's piezoelectricity therefore needs the dipole to respond to strain,
which requires bond stretching combined with charge flux, or charges that depend
on the backbone angle, or explicit polarizability. Fixed atom-centred charges on
a rigid-bonded chain cannot produce it in principle, not merely in practice.
