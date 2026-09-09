# Extending polyfind beyond PVDF: a design

Status: design only. None of this is implemented. It was produced as discussion
during the performance-optimisation session and is recorded here so it does not
disappear with the transcript. Claims about the current code have been checked
against the code and are marked as such; the compute estimate in section 4 is
an estimate and is marked as one.

Target chemistries: PVDF plus CFE, CDFE, AN, VDCN, VDC, CNEPO and FANOME, for
field-driven polarization-to-mechanical-work conversion, with a machine-learned
potential (MACE with D3) as the expensive validator.

## 1. What the current model cannot express

Verified against `src/polyfind/polymers.py` and `src/polyfind/pack.py`:

* `BackboneAtom` carries **one** `substituent` element, **one** `sub_bond`
  length and **one** `sub_charge`, all shared by both pendant atoms. The model
  is therefore restricted to a backbone atom bearing two identical, single-atom
  substituents.
* `UFF_LJ` defines parameters for carbon, hydrogen, fluorine and chlorine only.
  Nitrogen and oxygen are absent, so nitrile and methoxy groups cannot be
  scored at all.
* `CrystalPacker` places one or two chains per cell; `_place` has an explicit
  two-chain branch.
* Nothing anywhere couples the structure to an applied electric field.

VDF (CH2-CF2) and VDC (CH2-CCl2) fit the existing model. The rest do not: CFE
and FANOME need two different substituents on one carbon; AN needs a single
pendant rather than a pair; AN, VDCN and FANOME need multi-atom pendant groups;
and CNEPO's epoxide is a ring bridging two backbone carbons, which is a change
of topology rather than a substituent at all.

## 2. Phased extension

**Phase 1, symmetric drop-in (VDC, CH2-CCl2). Done.** No model change was
needed: chlorine was already parameterised, so PVDC is a `Polymer` entry with
bond lengths, angles and illustrative charges chosen the way PVDF's were. The
geometry builds exactly as specified, the whole funnel runs on it (fit,
enumerate, helix analysis, packing), and the flipped-chain fix holds for it.
`tests/test_pvdc.py` covers all of that.

It also produced the first real finding of this exercise, and it is a caution
about the phases that follow. With bond angles frozen, the all-trans PVDC chain
carries about 313 kcal/mol of Lennard-Jones strain, against 9 for PVDF: a
planar zigzag simply cannot relieve Cl...Cl contact when the angles cannot
open. Every rotation away from trans then lowers the energy, the fitted
first-order energies come out around -15 kcal/mol relative to a reference that
is not even metastable, and the RIS convention of measuring from all-trans
breaks down. The *direction* is correct, since PVDC is known not to adopt the
planar zigzag that PVDF's beta phase does, but the magnitude is an artifact of
rigid geometry rather than a property of the polymer.

The lesson generalises: the rigid-bond-angle assumption is already flagged as a
compromise for PVDF, where it is second order. For any substituent bulkier than
fluorine it becomes first order, and both the chlorine chemistries here and the
multi-atom pendant groups of phase 3 will need variable backbone angles from the
start, not merely at the final refinement stage. Performance-plan item 5 adds
exactly that capability; the chemistry work should build on it rather than
duplicate it. A second question it raises: for chemistries whose planar chain is
not metastable, the RIS reference state itself may need to be something other
than all-trans.

**Phase 2, asymmetric single-atom substituents (CFE CH2-C(F)(Cl), CDFE
CHCl-CF2).** Give `BackboneAtom` two independent substituent specifications
instead of one shared set. Default both to the same value so PE, PVDF and VDC
are untouched and every existing test still passes. The change is confined to
the dataclass and to wherever `chain.py` places the two pendants assuming they
are identical. Small to medium.

**Phase 3, multi-atom pendant groups (AN CH2-CH(CN), VDCN CH2-C(CN)2, FANOME
CH2-C(CN)(OCH3)).** A substituent becomes a small rigid fragment rather than an
atom. Nitrile is easy: it is linear, so it is two more atoms colinear with the
substituent bond and adds no rotational degree of freedom. Methoxy is the hard
one, because it has its own internal C-O rotation; the defensible first cut is
to freeze it at one representative rotamer, at the same level of simplification
as the rest of the rigid-geometry model, and to flag that explicitly rather
than adding an RIS dimension immediately. Needs nitrogen and oxygen nonbonded
parameters and illustrative charges for the nitrile and methoxy dipoles,
refitted through the existing `fit_ris` machinery exactly as PVDF's were.
Medium.

**Phase 4, ring backbone (CNEPO).** The epoxide bridges two backbone carbons.
The cheapest path is to reuse `RISModel`'s existing `clamp` mechanism, already
used by the sampler and the k-best search, to pin that bond's torsion to a
single fixed state, which removes the free-rotation assumption there without
any new dynamic-programming machinery. The ring's fixed local geometry, with
its bridging oxygen as a non-backbone atom, then becomes new construction logic
in `chain.py`. This is the most novel piece and the likeliest to need
iteration, but also the least urgent, since CNEPO reportedly already has
accepted relaxations and is not blocked on it.

## 3. Orthogonal to the phases

* **Copolymer composition.** A comparison model of one candidate unit per
  eleven VDF units (8.33 mol%) mostly falls out of the existing bond-type
  indexing, `t(i) = i mod B`, once `Polymer.backbone` can hold a full explicit
  24-bond sequence rather than a short repeating motif. Reuse VDF's fitted
  energies for the 22 background bonds and fit only the two candidate bonds and
  their junction terms. Real work, but configuration more than new mathematics.
* **More than two chains per cell.** A 2x2x1 eight-chain cell is what lets
  chains slip and register independently. `CrystalPacker` would have to
  generalise, and the search dimensionality grows quickly, so this should
  follow the tabulated chain-pair interaction of performance review section 1
  rather than duplicate it.
* **Field coupling.** The cheapest useful addition is a `-mu_cell . E` term in
  the lattice energy, built from the charges and coordinates already present,
  giving a fast way to scan field-induced conformational and packing shifts. It
  is a screening heuristic, not a substitute for a Berry-phase or DFT
  treatment.

Suggested order: phase 1, then phase 2, which unlocks VDC, CFE and CDFE with no
new physics, then phase 3, which unlocks AN, VDCN and FANOME, then phase 4.
Copolymer composition follows when a specific chemistry needs the like-for-like
comparison. Multi-chain packing and field coupling are their own later
initiative, each gated on the performance work it depends on.

## 4. The wider argument: calibrate cheap, search exhaustively, confirm expensive

The observation behind the phases is that the expensive campaign spends
machine-learned-potential compute on a search that does not need it.

Relaxing a full eight-chain, roughly 96-repeat supercell with MACE and D3 from
a handful of generic starting guesses per chemistry is expensive and offers no
convergence guarantee. It is a local optimiser searching a combined
discrete-torsion and continuous-packing space from generic starts, which is the
shape of failure behind starts that exhaust their step budgets. Every run
re-derives the whole energy landscape from scratch on a multi-thousand-atom
system, even though that landscape separates into a torsion problem that is a
Markov chain solvable exactly in O(N S^2), and a packing problem that is
pairwise between chains and therefore tabulatable and searchable by FFT. Both
are things polyfind already does for PVDF in seconds.

The alternative, per chemistry: fit an RIS and pair-interaction surrogate from
a few thousand single-point evaluations on a short oligomer, tens of atoms
rather than thousands; run the exact torsion search and the global packing
search classically; then hand the expensive potential two or three already-good
candidates to confirm and polish, which converge quickly because they start
near a true minimum rather than at a generic guess.

Rough accounting, and this is an estimate rather than a measurement: current
spend is of order 10^5 to 10^6 force evaluations on 10^3-atom systems per
chemistry, against 10^3 to 10^4 single points on 10 to 50-atom systems for
calibration plus 10^2 to 10^3 on the full system for confirmation. That is two
to three orders of magnitude, and it also stops missing minima the way a few
hand-chosen starts can.

The same restructuring answers several of the campaign's questions more cheaply
rather than merely faster. Competing structures are exactly what
`enumerate_periodic` plus an exhaustive packing search produce, classically and
for free once calibrated. A coarse barrier estimate can be read off the
tabulated energy surface, since a torsion flip or a slip coordinate is a table
lookup, leaving one expensive nudged-elastic-band run to confirm the winning
path instead of one per candidate pair per chemistry. Blocking stress and
recoverable strain fall out of a field term as a curvature calculation once
`-mu . E` is present. Questions about which electronic-structure method is
trustworthy at all stay out of scope, because no classical surrogate can
adjudicate them; the point is to shrink how often that method must run, not to
replace it.

**The load-bearing assumption**, stated plainly: that a short-oligomer fit
transfers to full-crystal energetics. polyfind already makes it for PVDF and it
holds qualitatively there, but it is untested for nitrile, ring and copolymer
chemistries. It is also directly checkable, by calibrating against relaxations
that already exist for PVDF and CFE before trusting the approach on the rest.
That check should come before the compute-saving claim is believed.
