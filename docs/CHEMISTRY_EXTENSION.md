# Extending polyfind beyond PVDF: a design

Status: phases 1, 2 and 3 are implemented; everything after them is design. The design was
produced as discussion during the performance-optimisation session and is
recorded here so it does not disappear with the transcript. Claims about the
current code have been checked against the code and are marked as such; the
compute estimate in section 4 is an estimate and is marked as one.

Target chemistries: PVDF plus CFE, CDFE, AN, VDCN, VDC, CNEPO and FANOME, for
field-driven polarization-to-mechanical-work conversion, with a machine-learned
potential (MACE with D3) as the expensive validator.

## 1. What the current model cannot express

Verified against `src/polyfind/polymers.py` and `src/polyfind/pack.py` as the code
stood before phase 2; the first bullet has since been lifted, see phase 2 below.

* `BackboneAtom` carried **one** `substituent` element, **one** `sub_bond`
  length and **one** `sub_charge`, all shared by both pendant atoms. The model
  was therefore restricted to a backbone atom bearing two identical, single-atom
  substituents. *(Lifted in phase 2: each of the three now takes one value or a
  pair. Still one atom per pendant, and still no tacticity.)*
* `UFF_LJ` defines parameters for carbon, hydrogen, fluorine and chlorine only.
  Nitrogen and oxygen are absent, so nitrile and methoxy groups cannot be
  scored at all. *(Lifted in phase 3: N and O added from the same UFF table.
  `pack.MASS` still has neither.)*
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
about the phases that follow. With bond angles frozen *at PVDF's 114 deg for both
backbone carbons*, the all-trans PVDC chain carries about 313 kcal/mol of
Lennard-Jones strain, against 9 for PVDF: a planar zigzag simply cannot relieve
Cl...Cl contact when the angles cannot open. Every rotation away from trans then
lowers the energy, the fitted first-order energies come out around -15 kcal/mol
relative to a reference that is not even metastable, and the RIS convention of
measuring from all-trans breaks down. The *direction* is correct, since PVDC is
known not to adopt the planar zigzag that PVDF's beta phase does, but the
magnitude is an artifact of rigid geometry rather than a property of the polymer.

That direction has since been checked against the crystallography and holds:
PVDC's crystal conformation is a *glide TGTG' form* with internal rotation
angles of 175 deg and 49 deg, in a monoclinic cell with the fibre axis at
4.68 A. The same structure also shows that the diagnosis above is right for a
concrete reason - the real chain relieves the Cl...Cl crowding by opening the
C-CH2-C backbone angle to 123 deg, against 114 deg at the CCl2 carbon, which is
precisely the degree of freedom this model freezes. See `docs/REFERENCES.md`
section 6.

**Those measured angles have now been adopted, and it settles the diagnosis.**
`polymers.py` carried 114 deg at both carbons, borrowed from PVDF's justification,
which does not transfer. With 123/114 the same ten-bond all-trans measurement
gives **74 kcal/mol against PVDF's 11** - a factor of 4.2 removed by nothing but
using the right angle - and the drop available by rotating every bond away from
trans falls from 79 kcal/mol to 4.0, so all-trans goes from grossly unphysical to
very nearly a local minimum. The residual 74 is what rigid geometry still costs;
the 239 that vanished was the wrong input.

Two consequences, one good and one a real loss of capability.

*The good one.* Ask the model which torsions make a TG+TG- repeat close under
these backbone angles, and it answers **175.3 deg and 49.4 deg with a chain repeat
of 4.677 A** - the published torsion pair to within half a degree and the measured
4.68 A fibre repeat to 0.06%. Under the old equal angles the same closure search
could only answer "trans exactly 180, gauche 49.6" and gave 4.491 A, 4.0% short,
and the published 175/49 pair did not close at all. So the wide CH2 angle is what
*creates* the trans deflection the crystallography reports: the 9 deg of curl per
repeat that unequal angles introduce is exactly what a 5 deg deflection of each
trans bond cancels. That is a structural prediction recovered from geometry alone.

*The loss.* That same 9 deg of curl means **no PVDC conformation closes at ideal
RIS angles any more**. An all-trans chain with unequal backbone angles is a
planar circular arc, not a straight stem, with zero rise per repeat; ideal
TG+TG- carries the same 9 deg. `periodic_chain` refuses every one of them, so the
phase-1 claim that "the whole funnel runs on PVDC" now holds only through
`periodic_chain_from_torsions` with the deflected torsions above. The RIS fit,
the enumeration and the helix analysis are unaffected; only ideal-angle packing
is lost, and it was packing a chain the polymer does not adopt.

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
CHCl-CF2). Done.** `BackboneAtom`'s three per-pendant fields (`substituent`,
`sub_bond`, `sub_charge`) each now accept *either* one value shared by both
pendants *or* an explicit `(first, second)` pair; the positional signature and
its old meaning are unchanged, so PE, PVDF and PVDC build byte-identical
coordinates and charges (asserted by digest in `tests/test_cfe_cdfe.py`).
`substituent_positions` takes the same one-or-two form for the bond length, which
is why the coordinate builders in `pack.py` and `linegroup.py` -- which forward
`spec.sub_bond` themselves rather than going through `build_chain` -- needed no
change. The two pendants keep an **equal split** of `sub_angle` about the
backbone bisector: any split preserves the specified inter-substituent angle, and
an unequal one (VSEPR would tilt the bisector towards the smaller pendant) would
need a parameter the monomer model does not carry. CFE and CDFE are registered
with illustrative geometry and charges chosen the way PVDF's and PVDC's were.

The finding, and it is the important part. A backbone atom with two different
pendants is a stereocentre, so the chain now has a tacticity, and the builder
silently picks one: pendants are placed on a fixed side of the local frame
`(previous, this, next)`, so the invariant `(s1 - C).[(prev - C) x (next - C)]`
has the same sign at every backbone atom, conformation-independently. **The chain
built is the isotactic one** (verified two ways in the tests: the sign invariant,
and the classical planar-zigzag same-side check). Syndiotactic and atactic chains
are not expressible; a `Polymer.backbone` holding an explicit multi-monomer
sequence, as already proposed for copolymers in section 3, is the natural way in.

*(Written during phase 2; the three places it lists as still assuming an achiral
repeat have since been fixed, see "Phase 2 outcome" below.)* The consequence is
that **`fit_ris`'s `symmetrize=True` default is no longer sound** for these
chemistries. Mirror symmetrisation is exact for an achiral
chain because reflecting a conformer maps `phi -> -phi` and returns the same
molecule; for a chiral repeat the reflection returns the *enantiomeric chain*, so
`E(phi)` and `E(-phi)` are energies of diastereomers and genuinely differ -- by
hundreds of kcal/mol with the illustrative potential, against exactly zero for
PVDF. Averaging them destroys precisely the asymmetry that makes an isotactic
chain select a one-handed helix. Fit these with `symmetrize=False`.
`Polymer.is_chiral` reports it, and `build_chain` warns once per polymer, because
`fit_ris` is where the unsound default lives and the warning has to reach a user
who never reads the polymer definition. Two further places still assume an
achiral repeat and are left as they are: `helix.canonical_sequence` (hence
`enumerate_periodic`) deduplicates a sequence against its G+/G- mirror, which for
a chiral chain discards a distinct conformer rather than a redundant one, and
`linegroup` counts a glide as a chain symmetry, which a chiral chain lacks.
Both belong with tacticity support rather than with this phase.

There is a correct replacement for `symmetrize`, and it is cheap. The symmetry a
chiral chain does possess is mirror *composed with chain reversal*, since
reversing the chain direction swaps `prev` and `next` and flips every
stereocentre's configuration back: `E(phi_1..phi_N) = E(-phi_N..-phi_1)`, verified
to machine precision here for CFE, CDFE and PVDF alike. In fitted-model terms,
measured for B = 2, that reads `e1[b, s] = e1[1 - b, m(s)]` and
`e2[b, s, s'] = e2[b, m(s'), m(s)]` -- a transpose and a bond-type swap on top of
the state mirror. Both hold to grid noise (< 0.1 kcal/mol at step = 30 deg) for
CFE and CDFE, while plain mirroring is violated by 12-45 kcal/mol; for PVDF both
hold, which is why the existing code is right for the achiral chemistries.
Teaching `symmetrize` that operation is a small change to `fit_ris`, left out of
this phase only because it touches a file phase 2 did not own.

**Phase 2 outcome: chirality, which the model had no way to express.** A
backbone carbon bearing two different substituents is a stereocentre, so CFE
and CDFE chains have a tacticity. The chain this code builds is unambiguously
the isotactic one: the first pendant always sits on the same side of the local
backbone frame, independent of conformation. Syndiotactic and atactic chains
are not expressible, and nothing warned about that before.

More seriously, reflecting a chiral chain gives its enantiomer rather than the
same molecule, so its G+ and G- conformers are genuinely inequivalent. The
fitting step averaged them by default, which for these two chemistries destroys
a real asymmetry of 24 kcal/mol (CFE) and 45 (CDFE) - and that asymmetry is
precisely what makes an isotactic chain choose a one-handed helix, so averaging
it away removes the physics the chemistry was added to study. `fit_ris` now
defaults to `symmetrize="auto"`, which mirrors only achiral polymers and warns
if mirroring is forced on a chiral one. `Polymer.is_chiral` exposes the
distinction and `build_chain` warns once per chiral polymer.

The two related places that also assumed achirality are now fixed as well, and
one of them corrected a claim made here. `helix.sequence_images` takes a
`chiral` flag that `enumerate_periodic` fills in from `Polymer.is_chiral`;
`linegroup.torsion_pattern` takes the same flag from `line_group`. The claim
that needed correcting was that chain reversal remains a symmetry of a chiral
chain on its own. It does not: reading an isotactic chain backwards swaps `prev`
and `next` at every stereocentre, so `E(phi_N..phi_1) == E(-phi_1..-phi_N)` --
the *enantiomer's* energy. Measured on 24-bond oligomers, reversal alone and
mirroring alone give bit-identical energies and both differ from the original by
up to 750 kcal/mol for CFE and 5100 for CDFE, against exactly zero for PE, PVDF
and PVDC. So the chiral group is shifts plus reflection-with-reversal and
nothing else, which is the same operation the fit uses. With it, PVDF's and
PVDC's candidate lists are unchanged to the name, CFE's grows from 69 to 90
candidates and CDFE's from 66 to 92, and the largest energy difference between a
candidate and its mirror partner is 31.1 kcal/mol per monomer for CFE and 38.8
for CDFE -- exactly zero for the achiral ones, as it must be.

A glide is an improper isometry and a chiral object has none, whatever its state
sequence reads like, so `torsion_pattern` no longer offers one to a chiral
repeat: `TG+TG-` and `TTTG+TTTG-` fall back to `free` (the penalty method),
while `TT` and `TG+` keep the screw they also have. Screws are proper isometries
and are unaffected.

The replacement symmetry is implemented and is the default for chiral
polymers. An earlier note here said its index mapping was unresolved, because a
first check left a large residual even for PVDF where it must vanish. That
check was wrong: it swapped the bond-type index on the pair term, which the
relation does not do. With the correct form -- a state mirror plus a bond-type
shift on the first-order term, a state mirror plus a transpose on the pair term
-- the residual is 0.02 to 0.09 kcal/mol at step = 20 deg for PE, PVDF, PVDC,
CFE and CDFE alike, against 32 to 47 for plain mirroring of the chiral pair.

`fit_ris` now averages over reflection-with-reversal whenever the polymer is
chiral, and validates the relation numerically before applying it: every
bond-type shift is tried, and the averaging is used only if the residual is
within `reversal_tol`, otherwise the fit is left unsymmetrised with a warning.
So the noise averaging is restored for chiral fits without the physics being
assumed.

The third-order terms and the `adapt_angles` mirror now follow the same
operation, each measured before it is used. Writing the reversal as `j -> K - j`
on the bonds, a term over `n` consecutive bonds maps to the term starting `n - 1`
bonds earlier, so each order carries a bond-type shift one lower than the order
below it, together with the state mirror and its own indices read backwards:
`e3[b, x, y, z] == e3[(c-2-b) % B, m(z), m(y), m(x)]`. At `step = 20 deg` that
residual is 0.00 kcal/mol for PVDF, CFE and CDFE, 0.02 for PVDC and 0.01 for PE,
against 40-100 (the cap) for every near miss -- mirroring without reversing the
triple, reversing without mirroring, or the right form at the wrong bond-type
shift -- on the chiral pair. The state angles obey `arg1[b, s] ==
-arg1[(c-b) % B, m(s)]` to 0.00 deg on a dense grid and 0.09 deg with the
adaptive scan, which is what makes the aggregated mirror-pair angles equal and
opposite for a chiral fit as well; the per-bond-type minima are 80 deg from
mirror-symmetric for CFE and CDFE, so it is the reversal that carries it, not
mirroring. Both checks are run at fit time against `reversal_tol` and
`angle_tol` and warn rather than guess: PVDC at `step = 45 deg` really does miss
the third-order check (residual 1.94 kcal/mol, from inclusion-exclusion between
two clashing conformers at +-45 deg), and the code declines to symmetrise there.

**Phase 3, multi-atom pendant groups (AN CH2-CH(CN), VDCN CH2-C(CN)2, FANOME
CH2-C(CN)(OCH3)). Done.** A substituent is now a `Pendant`: an ordered set of
`PendantAtom`s at fixed local coordinates in a per-pendant frame
(`chain.pendant_frames`) whose first axis is the pendant bond direction. A
single atom is the degenerate case -- one atom at `(bond, 0, 0)`, exactly where
`substituent_positions` always put it -- so `BackboneAtom`'s positional
signature is unchanged and PE, PVDF, PVDC, CFE and CDFE build byte-identical
coordinates and charges (asserted by digest in `tests/test_an_vdcn_fanome.py`).
A fragment may supply its own bond length and charges, in which case
`sub_bond`/`sub_charge` may be `None` and are filled in from it, which keeps
`spec.sub_bond` numeric for the callers in `pack`/`linegroup` that forward it.
`atoms_per_repeat` is now a sum over pendants rather than `3 * len(backbone)`.
Nitrile needs no new degrees of freedom (linear, both atoms on the pendant
axis); methoxy's two internal rotations are frozen -- methyl carbon anti to the
other pendant about the C-O bond, methyl itself staggered -- and the freezing is
stated in `methoxy`'s docstring and in the polymer comment. Nitrogen (3.660,
0.069) and oxygen (3.500, 0.060) come from the same UFF table as the existing
entries. `is_chiral`: AN and FANOME true, VDCN false -- two identical fragments
are placed as mirror images of one another, so a multi-atom pendant does not by
itself create a stereocentre.

The finding is the phase-1 one again, worse. All-trans Lennard-Jones strain, on
the same 10-bond oligomer that gives PVDF 11 kcal/mol and PVDC 74 (9 and 313
before the geometry corrections; AN, VDCN and FANOME are unaffected by them,
since neither PVDF's bond length nor PVDC's angles enter their definitions):
**AN 88, VDCN 176, FANOME 1.0e6**. AN and VDCN are strained but finite and comparable
with CFE (162) and CDFE (199); every rotation away from trans still lowers the
energy by tens of kcal/mol, so all-trans remains a poor RIS reference for them,
though a usable-if-suspect one in the same sense it is for the chlorine
chemistries. FANOME's all-trans chain is not a physical structure at all: methyl
hydrogens of methoxy groups on consecutive substituted carbons end up 0.80 A
apart. That is not an artifact of the frozen rotamer -- scanning the C-O
azimuth over 360 degrees gives between 7e3 and 2e8 kcal/mol, with no value
below 7e3 -- it is the frozen backbone angles again, and a pendant that reaches
2.4 A cannot get out of its 1-3 neighbour's way when they cannot open. So the
funnel is exercised end to end on VDCN and AN only; FANOME is registered and its
geometry tested, but fitting it against an all-trans reference would produce
numbers with no meaning. Variable backbone angles, and for FANOME a reference
state other than all-trans, are the prerequisites rather than the polish.

Two things elsewhere have not caught up, both out of the phase's file ownership
and both one-liners: `pack._batch_block_coords` and `linegroup._block_coords`
build coordinates themselves and index them as `3k + {0,1,2}`, so
`repeat_chains_from_torsions` and the continuous refinement stage still accept
single-atom pendants only; and `pack.MASS` has no N or O entry, so
`PeriodicChain.mass` and the densities derived from it raise for these three.
`build_chain`, `build_chain_batch`, the bond-angle override, `fit_ris`,
`enumerate_periodic`, `pack.periodic_chain` and `CrystalPacker.energy` are all
general.

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

### Phase 3 follow-up: the two gaps are closed

Phase 3 could not touch `pack.py` or `linegroup.py`, and left two things broken
for the new chemistries. Both are now fixed.

The batched block builders in both files laid atoms out at a fixed stride of
three per backbone atom, which silently assumed single-atom pendants. They now
compute the stride from each backbone atom's pendant sizes, matching
`build_chain`'s layout exactly, so `repeat_chains_from_torsions` and the
gradient path inside continuous refinement accept multi-atom pendants.
`pack.MASS` gained nitrogen and oxygen, without which any density involving a
nitrile or methoxy raised. Verified for AN and VDCN: chains build, masses and
densities come out, batched torsion perturbations reproduce the unperturbed
chain exactly, and the per-row gradient path returns finite energies.

The strain finding from phase 1 has now recurred twice more and is worth
stating as a general limit rather than a per-chemistry quirk. All-trans
Lennard-Jones strain, same measurement each time: PVDF 11, CFE 162, CDFE 199,
PVDC 74, AN 88, VDCN 176, FANOME about 1e6 kcal/mol. (PVDF and PVDC re-measured
after the geometry corrections, from 9 and 313; the other five are unchanged and
still carry the equal 114 deg backbone angles, which for CFE and CDFE in
particular is the same borrowed assumption PVDC's correction just overturned and
is therefore the obvious next thing to check.) Only PVDF has an all-trans chain
that is a sensible RIS reference, and PVDC is now the near miss rather than the
worst case. For everything bulkier the planar zigzag is strained enough that any
rotation lowers the energy, so the convention of measuring RIS energies from
all-trans is measuring from a state the polymer would never occupy.

FANOME is the extreme case and is instructive: its all-trans chain puts methyl
hydrogens of methoxy groups on consecutive substituted carbons 0.80 A apart.
That is not an artifact of freezing the methoxy rotamer, which was the obvious
suspect -- scanning the C-O azimuth through a full turn never drops the energy
below about 7e3 kcal/mol. It is the frozen backbone angles again, the same
limitation phase 1 identified, now severe enough to make the structure
unphysical rather than merely strained. FANOME is registered and its geometry
tested, but deliberately not fitted.

Two things follow for the roadmap. First, the RIS reference state needs to
become the relaxed chain rather than all-trans for these chemistries, or their
fitted energies are relative to nothing meaningful. Second, variable backbone
angles are no longer a refinement nicety; for anything past PVDF they are a
precondition for the model to describe a real molecule at all.
