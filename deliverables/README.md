# All-trans 11 VDF : 1 VDCN starting structures (8.33 mol%)

Produced by `polyfind` in answer to the request at the end of
`docs/NOTE_VDCN_CONVERGENCE.md`: *"an explicit periodic 11 VDF : 1 VDCN sequence with
topology-checked all-trans starts in independently constructed polar and antipolar
packings, screened for interchain close contacts"*, at the 592-atom eight-chain size.

**These are starting structures, not results.** Nothing here has been relaxed with a
machine-learned potential. They are meant to be fed straight into the fixed-cell then
full-cell MACE+D3 protocol.

**Read the "What we would not stand behind" section before using the energies.** The
structures are the deliverable; the energies are context.

## Reproduce

    PYTHONPATH=src python examples/copolymer_starts.py      # the aligned cells
    PYTHONPATH=src python examples/copolymer_staggered.py   # the staggered variants
    PYTHONPATH=src python examples/copolymer_readme.py      # writes this README from them

## Files


| file | contents |
|---|---|
| `vdf11_vdcn1_alltrans_polar_8chain_592atom.xyz` | polar, **592 atoms, eight chains** -- the requested size |
| `vdf11_vdcn1_alltrans_polar_2chain_148atom.xyz` | polar, 148 atoms, two chains -- the primitive cell the above tiles |
| `vdf11_vdcn1_alltrans_antipolar_8chain_592atom.xyz` | antipolar, **592 atoms, eight chains** -- the requested size |
| `vdf11_vdcn1_alltrans_antipolar_2chain_148atom.xyz` | antipolar, 148 atoms, two chains -- the primitive cell the above tiles |
| `vdf11_vdcn1_alltrans_polar_2chain.cif` | polar primitive as P1 CIF |
| `vdf11_vdcn1_alltrans_antipolar_2chain.cif` | antipolar primitive as P1 CIF |

Staggered variants, added 2026-09-11 (see *Staggered variants* below; these are the ones to run if you want the axial registry broken):

| file | contents |
|---|---|
| `vdf11_vdcn1_alltrans_polar_stagger_ladder_8chain_592atom.xyz` | polar, **592 atoms, eight chains**, `ladder` stagger: offsets [0, 3, 6, 9, 0, 3, 6, 9] monomers, 4 distinct axial heights, rho 1.6799 |
| `vdf11_vdcn1_alltrans_polar_stagger_spread_8chain_592atom.xyz` | polar, **592 atoms, eight chains**, `spread` stagger: offsets [0, 3, 6, 9, 1, 4, 7, 10] monomers, 8 distinct axial heights, rho 1.6789 |
| `vdf11_vdcn1_alltrans_antipolar_stagger_ladder_8chain_592atom.xyz` | antipolar, **592 atoms, eight chains**, `ladder` stagger: offsets [0, 3, 6, 9, 0, 3, 6, 9] monomers, 4 distinct axial heights, rho 1.5707 |
| `vdf11_vdcn1_alltrans_antipolar_stagger_spread_8chain_592atom.xyz` | antipolar, **592 atoms, eight chains**, `spread` stagger: offsets [0, 3, 6, 9, 1, 4, 7, 10] monomers, 8 distinct axial heights, rho 1.5697 |
| `vdf11_vdcn1_alltrans_polar_stagger_ladder_8chain.cif` | the same cell as a P1 CIF (fractional, wrapped) |
| `vdf11_vdcn1_alltrans_polar_stagger_spread_8chain.cif` | the same cell as a P1 CIF (fractional, wrapped) |
| `vdf11_vdcn1_alltrans_antipolar_stagger_ladder_8chain.cif` | the same cell as a P1 CIF (fractional, wrapped) |
| `vdf11_vdcn1_alltrans_antipolar_stagger_spread_8chain.cif` | the same cell as a P1 CIF (fractional, wrapped) |
| `staggered_summary.json` | every staggered number in this document, machine-readable |

And the machine-readable companions to both sets:

| file | contents |
|---|---|
| `summary.json` | every aligned number in this document, machine-readable |
| `ris_parameter_provenance.json` | per-term parameter source and match quality |

Extended XYZ: `Lattice="..."` on the comment line, `Properties=species:S:1:pos:R:3`, `pbc="T T T"`, plus the energy, density and minimum interchain distance as extra key-value pairs; the staggered files add `stagger_monomers="..."`, one integer per chain in the file's own chain order. Atom order is chain by chain, 74 atoms each, and within a chain it is backbone atom then its pendant atoms, repeating.

**Coordinates are not wrapped into the cell.** Each chain's crystallographic repeat is written whole, so a few pendant atoms sit just outside the `c` boundary. That is deliberate: it is the molecular choice of branch, and it is what makes the cell dipole quoted below a property of the cell rather than of where the origin was put. Any PBC-aware reader handles it; wrap if your tool insists. The CIFs are the exception: fractional coordinates are taken modulo one, because a CIF reader would otherwise place the overhanging atoms outside the box. Use the extended XYZ for anything that reads a dipole off the coordinates.

**What we could and could not validate about the format.** ASE is still not installed here, so these have not been round-tripped through the reader you will probably use. They have been round-tripped through an independent parser of our own, which is worth more than nothing and less than ASE: `tests/test_supercell.py` re-reads each staggered file as text and re-derives the energy, the density, the whole topology report and the stagger pattern from what comes back. Energies agree to 1e-8 kcal/mol per monomer, densities to six decimals, the topology verdict and safe window are identical, and the stagger recovered from the coordinates alone matches the declared `stagger_monomers`. Each header's extra key-value pairs are checked against the file's own geometry rather than against the generator's memory of it. Still spend the ten seconds to load one in your own stack before committing a long run.

## The sequence

`polyfind.polymers.VDF_VDCN_11_1` (registered as `"vdf11-vdcn1"`). `Polymer.backbone`
now holds an **explicit multi-monomer repeat** rather than one cycled monomer, so this is
one 24-bond periodic repeat of twelve monomers: eleven VDF and one VDCN, the VDCN at
monomer index 6 (backbone atoms 12 and 13). Every homopolymer is the degenerate
one-monomer case of the same field and is byte-for-byte unchanged.

One consequence is worth stating because every parameter claim below rests on it: in this
model **exactly one of the 24 backbone atoms differs from PVDF's**. VDCN's own CH2 entry is
field-for-field identical to PVDF's (same element, same two H pendants, 114/108 deg,
-0.20/+0.10 e), so the substitution is confined to the cyano carbon, atom 13, and the
junction it creates is two backbone bonds wide. The copolymer's backbone is therefore
*geometrically* PVDF's all-trans backbone, and its repeat is exactly twelve times the
homopolymer's c.

### Which bonds got which parameters

Assigned by `polyfind.ris.transfer_ris`, which matches each copolymer bond's backbone
environment against the comonomers' fitted homopolymer models -- the atoms the rotating
bonds join first, then the whole dihedral window -- and reports the match rather than
asserting it. Full table in `ris_parameter_provenance.json`.

| RIS term | bonds | source | fitted for this environment? |
|---|---|---|---|
| first order | 20 of 24 | PVDF, bond type `b % 2` | yes, exactly |
| first order | 11, 12 | **VDCN** types 1 and 0 | no: the rotating bond is the CH2-C(CN)2 / C(CN)2-CH2 bond that VDCN does have, but its far neighbour is CF2 rather than C(CN)2 (3 of 4 window atoms match) |
| first order | 10, 13 | **PVDF** types 0 and 1 | no: the rotating bond is a plain VDF CH2-CF2 bond, but the 1-4 partner across it is the cyano carbon (3 of 4 match) |
| second order | 19 of 24 | PVDF | yes, exactly |
| second order | 9, 10, 11, 12, 13 | PVDF / VDCN (see the JSON) | no |
| third order | 18 of 24 | PVDF | yes, exactly |
| third order | 8..13 | PVDF | no |

**The junction terms have no fitted values at all.** Bonds 10 and 13 are the two junction
bonds either side of the VDCN unit; bonds 11 and 12 are the two bonds the cyano carbon
itself sits on. Four of 24 first-order terms, five of 24 pair terms and six of 24 triple
terms are transferred rather than fitted. Fitting them needs dihedral scans of a
VDF-VDCN-VDF oligomer, which is one contained job and has not been done.

None of this affects the geometry in these files: an all-trans start is all-trans at
exactly 180 degrees whatever the torsional energies say. The provenance matters for
anything that later *ranks* copolymer conformations.

## The potential used

`polyfind.pack.CrystalPacker` at its defaults: UFF Lennard-Jones (Rappe et al. 1992,
including the N entry added with the phase-3 nitriles), damped-shifted-force Coulomb at an
8 A cutoff with alpha = 0.2 and eps_r = 1, plus the illustrative point charges that
`polymers.py` carries. The Ewald energies quoted alongside are the same geometry
re-evaluated with `coulomb="ewald"` and the default tinfoil boundary, because a truncated
sum cannot evaluate a dipole lattice sum and so cannot be trusted for a polar/antipolar
*energy* comparison.

The charges are illustrative, not fitted. The fitted PVDF potential
(`pvdf-dft-valence-flux`) has no nitrile parameters, so it could not be used here.

## How the two packings were built

The packer places **two** chains per cell, so the eight-chain 592-atom cell is a
**2x2x1 tiling** of the two-chain primitive. That reaches the size, the composition and the
formula exactly -- `C208H192F176N16`, 592 atoms, eight chains of 74 -- but the tiling is an
exact replication, so the four copies of the two-chain motif are **not independent**: a
chain cannot slip, rotate or register differently from its own image, because it is its own
image. If independent eight-chain registry is what you need, these are seeds for it rather
than an answer to it -- and *Staggered variants* below is that answer, built on top of these
files with eight independently registered chains.

One thing to be precise about, because we were not precise enough about it the first time:
the packer's `dz` slides chain 2 against chain 1, and at these cells it does so by
0.99 monomers (polar) and 1.49 (antipolar). So the eight chains sit at **two**
axial heights, not one. What a 2x2x1 tiling locks is the registry within each of the two
sublattices: the 4.80 A `b`-axis neighbours and the 10.5 A `a`-axis ones.

Both cells are orthorhombic with gamma fixed at 90 degrees. Your full-cell stage is what
releases that.

**Polar.** The PVDF homopolymer's all-trans two-chain cell was packed first (it is cheap and
exhaustive on a two-bond repeat) and used as a seed, because the copolymer's backbone is
geometrically identical to it. The chain is all but exactly twelve-fold periodic along z, so
`dz` and `dz + c/12` are two *different* VDCN-VDCN axial registries with the same backbone
packing, and a local optimiser cannot move between them: all twelve were taken as separate
starts, for each chain-direction flip. A 4000-cell uniform random screen over
`(a, b, phi1, phi2, dz)` was run alongside as insurance against a qualitatively different
packing the nitrile might prefer. The best distinct cells were polished with the exact
kernel and analytic gradients.

**Antipolar, constructed independently.** `polyfind.fitting.antipolar_offsets` derives the
exactly-antipolar subspace from the chain's own dipole moment rather than assuming one, and
`antipolar_cell_exact` searches it. A second, finer scan of the *same* subspace was run as
well, because `antipolar_cell_exact` samples `dz` at four points across the repeat and for a
twelve-monomer repeat that is a 7.7 A step -- coarser than the 2.56 A monomer period the
interchain registry actually varies on. The lower of the two was taken.

**That mattered, and it is worth reporting against our own helper.** `antipolar_cell_exact`
returned -3.7152 kcal/mol per monomer; the finer scan of the identical subspace returned
**-3.7514**, and the cell shipped here is the finer one. So the helper's four-point `dz` grid
did miss the basin on a repeat this long. Both cells are exactly antipolar -- the subspace is
the same and the constraint is exact either way -- and only the axial registry differs, so
this is a resolution limit in the screen rather than a correctness bug. It would bite any
copolymer, whose repeat is by construction many monomers long.

**The polarization was checked, not assumed.** A defect in exactly this area was found and
withdrawn in this repository (`docs/SCREEN.md`, addendum 3): the older
`antipolar_cell` helper used a flip with *equal* setting angles, which is antipolar only
when a chain's transverse moment is perpendicular to its own reference axis -- true of the
alpha helix, false of every planar zigzag, so on beta-PVDF it returned a cell with the full
polarization of the polar minimum. This chain is a planar zigzag and its moment is
`(+5.83, 0, 0)` e.A, lying along its own x, which is precisely the case that defect got wrong.
The dipole of the cell shipped here is reported below and is zero to floating point.


## The two cells

| | polar | antipolar |
|---|---|---|
| a (A) | 10.4598 | 11.3267 |
| b (A) | 4.8036 | 4.7451 |
| c (A), = chain repeat | 30.7557 | 30.7557 |
| gamma (deg) | 90.0 | 90.0 |
| phi1 (deg) | 90.00 | 90.00 |
| phi2 (deg) | 90.00 | 270.00 |
| dz (A) | 2.5304 | 3.8118 |
| chain 2 antiparallel | no | no |
| density (g/cm3) | 1.6816 | 1.5720 |
| E per monomer, truncated (kcal/mol) | -4.6082 | -3.7514 |
| E per monomer, Ewald (kcal/mol) | -6.9933 | -6.1225 |
| cell dipole (e.A) | (+2.018e-08, +11.66, -3.26e-13) | (+1.332e-15, -6.661e-16, -3.251e-13) |
| P (C/m2) | (+2.093e-10, +0.1209, -3.38e-15) | (+1.291e-17, -6.456e-18, -3.151e-15) |
| P magnitude (C/m2) | 1.209e-01 | 3.151e-15 |

E(antipolar) - E(polar) = **+0.8569** kcal/mol per monomer truncated, **+0.8708** with Ewald. Positive means the polar cell is the lower of the two under this potential.


**Which axis is which.** The polar cell's polarization lies along **b** (4.80 A), so that is the polar axis -- the short transverse one, as in beta-PVDF, whose polar axis this packer puts at 4.64 A. The long axis is the one the nitriles point along. The search does not know which axis is which and can return them either way round (`DESIGN.md` 5.2), so compare by length rather than by name.

**Which antipolar branch.** Both cells have `flip = 0`, i.e. the two chains run the *same* way. The antipolar cell is antipolar by setting angle, `phi2 = phi1 + 180`, not by an up-down chain pair. `antipolar_offsets` offers both branches for this chain (its moment has no axial component, so `flip = 0, dphi = 180` and `flip = 1, dphi = 180` are each exactly antipolar) and the search picked this one on energy.

The antipolar cell's dipole is 3.25e-13 e.A in its largest component, i.e. zero to floating point, and it was **measured rather than assumed** (see above). The polar cell's polarization magnitude is 0.1209 C/m2, which is what makes it polar.

## Density: read this before you use these

**These cells are looser than your accepted reference, and that is the safe direction for
your protocol, but you should know why.**

Reference points, all measured with the same default packer as above, all-trans, two chains
per cell:

| | a x b x c (A) | density | E per monomer |
|---|---|---|---|
| PVDF homopolymer | 4.64 x 8.62 x 2.563 | 2.073 | -7.029 |
| VDCN homopolymer | 9.50 x 6.23 x 2.58 | 1.696 | +24.932 |
| this copolymer, polar | 4.80 x 10.46 x 30.756 | 1.682 | -4.608 |
| your accepted kinked endpoint | 20.12 x 9.70 x 28.98, beta 110.4 | 1.96 | - |

(The two homopolymer rows reproduce with
`pack(periodic_chain(P, [0, 0], THREE_STATE), n_chains=2, screen="random", n_random=4000)`
for `P` in `PVDF`, `VDCN`. VDCN's `+24.9` is its own documented all-trans strain under rigid
bond angles, not a packing failure; see caveat 2.)

A mass-fraction mixing rule over those two homopolymer densities predicts **2.028** for
this composition. We get **1.682**, about **17% less dense**. The shape says
where it goes: the copolymer's short axis comes out **4.80 A**, 4% off PVDF's own
4.64 A, while its long axis comes out **10.46 A**, 21% above PVDF's 8.62 and
10% above the 9.50 A long axis of the *pure VDCN* cell. So the majority
monomer sets the short axis and the nitrile sets the long one, even though only one monomer
in twelve carries a nitrile.

> **Corrected 2026-09-11.** This paragraph used to continue: *"The reason is the tiling rather
> than the chemistry: the 2x2x1 replication puts every chain's nitrile at the same axial
> height, turning eight isolated bulky groups into a continuous plane of them that the whole
> structure has to clear."* **We tested that and it is wrong**, twice over. First, the tiled
> cells never did put every nitrile at the same height: the packer's own `dz` already offsets
> the second chain, and at the shipped cells it does so by 0.99 monomers (polar) and
> 1.49 (antipolar), so four chains sit at one height and four at another. What the tiling
> actually locks is the registry *within* each sublattice -- the 4.80 A `b`-axis neighbours and
> the 10.5 A `a`-axis ones. Second, breaking that lock does not recover the density: see
> *Staggered variants* below. The deficit is real and the sentence naming its cause was a
> guess; the sections below say what we now think it is.

This is not a search failure, and that was checked rather than assumed. Both checks are in
the generator and their numbers are in its log.

*Shape probe.* Fixing `(a, b)` at six shapes whose areas span 1.69-1.95 g/cm3 -- including
your accepted density -- and sampling the setting angles and the axial shift densely at each
leaves the dense shapes **repulsive**: the best of 120 samples is **+5.7** kcal/mol per
monomer at 1.954 g/cm3, +11.6 at 1.878 and +3.4 at 1.760. Only the shapes already close to
the reported cell sample attractively, and polishing the best three with every variable free
brings all three back to **the same cell** at -4.61 per monomer. So the polar cell above is
the minimum of this two-chain problem, not wherever the random screen happened to look.

*Registry sweep.* The twelve axial VDCN-VDCN registries at that cell span -4.61 to -3.11 per
monomer, and the one reported is the lowest: eleven of the twelve are within 0.11 of each
other and one is 1.5 higher. `dz` is the one cell variable a local polish cannot cross, so it
was swept rather than trusted -- but it is not a hidden lever here, because the nitriles point
into the wide inter-sheet gap and barely see each other.

What we would not conclude from this is that an all-trans copolymer cannot pack densely. The
tiling constraint alone could account for all of it, and an eight-chain search with
independent axial registries -- which `CrystalPacker` cannot do, it places two chains -- would
stagger the nitriles and might recover much of the 17%. It is, though, consistent with your
kink being how the chain makes room for the nitrile at high density, which would be a real
result about the copolymer rather than about either of our codes.

**That paragraph was a prediction, and it has since been tested.** The eight-chain search it
asks for is in *Staggered variants* below, and the answer is in that section's verdict table
rather than here. Read the two together: this section says where the deficit goes, that one
says whether breaking the registry gets it back.

For your protocol the practical consequence is the useful direction: a loose fixed-cell start
cannot do what your rejected seed did. The seed that failed was too *tight* -- chains close
enough to form interchain bonds -- and these are nowhere near that (see the contacts below).
Your full-cell stage is what should bring the volume down.


## Topology check

This is the point of the request, so it is done against a criterion like the consumer's
rather than only against our own geometry. `polyfind.topology` builds the *intended* bond
graph from the `Polymer` definition -- which monomers, which pendants, which intra-pendant
bonds, plus the one backbone bond that closes each chain onto its own next repeat -- and
compares it with the graph a distance-based detector would find:

    bonded  iff  d_ij <= scale * (r_i + r_j)

with `r` the Cordero covalent radii that `ase.data.covalent_radii` carries (C 0.76,
H 0.31, N 0.71, F 0.57 A). Every periodic image is examined, so nothing is hidden by a
minimum-image approximation.

A pass at one `scale` would not be worth much. What is reported is the **window of scales
over which the detected graph is exactly the intended graph**: below its lower edge an
intended bond goes undetected, above its upper edge a spurious bond appears. A wide window
means no reasonable choice of detector disagrees.


### Results

| | polar primitive | polar 8-chain | antipolar primitive | antipolar 8-chain |
|---|---|---|---|---|
| atoms | 148 | 592 | 148 | 592 |
| chains | 2 | 8 | 2 | 8 |
| formula | `C52H48F44N4` | `C208H192F176N16` | `C52H48F44N4` | `C208H192F176N16` |
| intended bonds | 148 | 592 | 148 | 592 |
| connected components at scale 1.2 | 2 x 74 | 8 x 74 | 2 x 74 | 8 x 74 |
| lost bonds (1.05-1.3) | 0 | 0 | 0 | 0 |
| new bonds (1.05-1.3) | 0 | 0 | 0 | 0 |
| new **interchain** bonds (1.05-1.3) | 0 | 0 | 0 | 0 |
| worst intended bond, d/(ri+rj) | 1.0187 | 1.0187 | 1.0187 | 1.0187 |
| closest non-bonded pair, d/(ri+rj) | 1.5844 | 1.5844 | 1.5844 | 1.5844 |
| **safe scale window** | 1.019 - 1.584 | 1.019 - 1.584 | 1.019 - 1.584 | 1.019 - 1.584 |
| **min interchain distance (A)** | **2.5900** | **2.5655** | **2.7904** | **2.5243** |
| min interchain d/(ri+rj) | 2.5392 | 2.5152 | 2.1800 | 2.4748 |
| closest intrachain non-bonded (A) | 1.7637 | 1.7637 | 1.7637 | 1.7637 |
| verdict | **PASS** | **PASS** | **PASS** | **PASS** |

**What that window means.** Any covalent scale between **1.019** and **1.584** recovers the intended graph of all four cells exactly: eight (or two) separate components of 74 atoms, no bond missing, no bond invented. The usual choices -- 1.05, 1.1, 1.15, 1.2, 1.3 times the sum of covalent radii -- all sit inside it. A bare **1.0** does not, and would lose every C-H bond; that is a property of the Cordero radii against a 1.09 A C-H bond and would do the same to any hydrocarbon, so if your detector uses 1.0 unscaled, check it against your own accepted structure before reading anything into it. The lower edge is set by the C-H bond at 1.09 A and the upper by a 1-3 backbone C...C at 2.41 A, both *intramolecular*, so the window is a property of the chain rather than of the packing; the packing has far more room than that. The closest interchain contact in either 592-atom cell sits at **2.47** times the sum of covalent radii, so a detector would have to be more than twice as generous as the most generous convention before it saw an interchain bond.


### Minimum interchain distance by element pair (A), 592-atom cells

| pair | polar | antipolar |
|---|---|---|
| H-N | 2.5655 | 2.5243 |
| F-N | 2.6268 | 2.7904 |
| C-H | 2.8167 | 2.7652 |
| F-H | 2.8323 | 2.7805 |
| C-N | 3.2028 | 3.4734 |
| C-F | 3.5754 | 3.5239 |
| C-C | 3.5886 | 3.5377 |
| N-N | 3.6083 | 4.2369 |
| F-F | 3.9008 | 3.7251 |
| H-H | 4.2171 | 4.0950 |

## Staggered variants: breaking the axial nitrile registry

You took up the offer at the end of our 2026-09-10 note, and this section is the answer.
A 2x2x1 tiling of a two-chain cell is an exact replication, so a chain cannot sit at a
different axial height from its own image -- it *is* its own image -- and we had argued that
the resulting nitrile registry was what cost the four files above about 17% density against a
mixing rule. The cells here slide each chain along its own axis by a whole number of monomers
so the comonomer units distribute along the repeat instead.

**The short version: in both polarities, breaking the registry does not recover the density (1.6816 -> 1.6809 g/cm3) and does not lower the energy (+0.0399 kcal/mol per monomer).** The verdict table below has the numbers and the honest
reading of them; the rest of this section is how we got there, because a result like this is
only worth anything if the search behind it is visible.

**Two corrections to what we sent you on 2026-09-10**, both found while building this:

1. The tiled cells never did put every nitrile at the same height. The packer's own `dz`
   offsets chain 2 against chain 1, and it does so by 0.99 monomers in the shipped polar
   cell and 1.49 in the antipolar one. Four chains sit at one height and four at another.
   What a tiling locks is the registry *within* each sublattice.
2. "It costs about 17% in density" was a guess at the cause, not a measurement. It is wrong.

**These are a third seed, not a prediction, and the caveats on the aligned files apply
unchanged.** The junction bonds either side of the VDCN unit still have no fitted torsional
parameters. The nitrile charges are still illustrative. And your two converged endpoints still
carry five and fifteen dihedrals outside any rotational-isomeric state, so every all-trans
start we can build -- aligned or staggered -- is **a different basin, not a better one**.

### Why a whole number of monomers

Because the copolymer's backbone is a PVDF backbone, a slide by `c / 12` is a *rigid* slide
of the skeleton and moves only the comonomer decoration. Measured rather than assumed:
a `c/12` slide maps **68 of 74** atoms onto an identical atom, worst displacement 7.5e-13 A. The 6 that move are exactly the comonomer site: `34` F -> C at 0.1297 A, `35` F -> C at 0.1297 A, `40` C -> F at 0.1297 A, `41` N -> F at 1.2817 A, `42` C -> F at 0.1297 A, `43` N -> F at 1.2817 A.

So an integer stagger changes which monomer of each chain carries the nitrile and changes
nothing else about the chain. A fractional slide would put the backbones of neighbouring
chains out of register as well, which is a different object; we release it at the end as a
half-monomer local relaxation and report what it buys.

### Which patterns, and why those

The eight chains of the 2x2x1 cell sit on a centred lattice, and the neighbour shells are not
what the cell edges suggest. Measured at the aligned polar cell, the close column pairs
are **4** pairs at 4.80 A and **16** pairs at 5.76 A. Two facts follow, and both pick the patterns. **Every** A-B pair is a close pair,
so no assignment of integer offsets can keep all close pairs more than 3 monomers apart --
that is a backtracking search over the constraint graph, not an estimate. And the A-B shell is
the one the packer's `dz` already moves, so the registry a tiling genuinely locks is the
within-sublattice one: `b` at 4.80 A and `a` at 10.5 A.

| pattern | offsets (monomers, chain order below) | min offset over close pairs | mean | distinct heights |
|---|---|---|---|---|
| `aligned` | 0 0 0 0 0 0 0 0 | 0 | 0.00 | 1 |
| `antiphase` | 0 6 0 6 0 6 0 6 | 0 | 4.80 | 2 |
| `sheet` | 0 0 0 0 6 6 6 6 | 0 | 2.40 | 2 |
| `ladder` | 0 3 6 9 0 3 6 9 | 3 **(optimal)** | 3.60 | 4 |
| `spread` | 0 3 6 9 1 4 7 10 | 2 | 3.60 | 8 |

Chain order is `build_cell`'s: `(0,0,0), (0,0,1), (0,1,0), (0,1,1), (1,0,0), (1,0,1), (1,1,0), (1,1,1)` as `(ix, iy, k)`, with `k = 1` the chain at
`(a + b) / 2`.

Rather than enumerate, the four non-trivial patterns **decompose the problem by neighbour
shell**, so that a null result can be attributed to a shell rather than to the whole idea:

* **`antiphase`** offsets the two sublattices by half a repeat and leaves both near shells
  aligned. It is *not* an independent structure -- `dz` already does exactly this -- so it is
  a check rather than a candidate, and it should come back identical to `aligned` after the
  registry sweep. It does, to every digit reported, which is how we know the eight-chain
  search is reproducing the two-chain one.
* **`sheet`** offsets the two sheets stacked along the long `a` axis by half a repeat and
  leaves the 4.80 A `b` neighbours aligned: the far shell alone.
* **`ladder`** ramps a quarter of the repeat per `b/2` layer, which puts the 4.80 A
  neighbours half a repeat apart -- the near shell, the one a tiling actually locks. It is
  affine in the transverse height, so it is a c-glide rather than an arbitrary decoration, and
  at 3 monomers it **reaches the bound above**: it is the maximally-separated arrangement
  for this lateral geometry. That was worth checking rather than assuming, because it is not
  obvious until the constraint graph is written down.
* **`spread`** adds one monomer per step along `a` on top of the ladder, so all eight chains
  sit at eight different heights -- the most axial slots occupied -- at the cost of dropping
  the nearest-neighbour minimum to 2 monomers. `ladder` and `spread` optimise different
  measures, which is why both are here and both are shipped.
* **`aligned`** is the 2x2x1 registry, run through the *identical* eight-chain search. It is
  the control that makes any difference a statement about the stagger rather than about the
  optimiser.

### The eight-chain energy

`CrystalPacker` places two chains, so it cannot relax a stagger at all.
`polyfind.supercell.SupercellEnergy` evaluates the same potential -- the same UFF
Lennard-Jones tables, the same damped-shifted-force Coulomb at 8 A, the same bonded
exclusions, the same charges, optionally the same Ewald sum -- on a cell of arbitrarily many
independent chains, as a direct sum over every periodic image within the cutoff. It is not
trusted on that description: the energy per monomer of a pair potential is invariant under
tiling, so it is checked against `CrystalPacker.energy` on the packer's own cells, and it
agrees to 2e-10 kcal/mol per monomer. The staggered cells and the aligned ones are
therefore scored by one function.

Each staggered cell was then relaxed over exactly the five variables the aligned cells were
-- `(a, b, phi1, phi2, dz)` at `gamma = 90` -- by the same sequence: a screen (the aligned
optimum at all twelve axial registries, nine fixed shapes spanning 1.68 to 2.10 g/cm3 with
the setting angles and `dz` sampled at each, and a uniform random screen), a polish of the
best five distinct cells, an axial-registry sweep, and a final tightening. The fixed dense
shapes matter: the question is whether a staggered cell *started dense* stays dense, and a
screen that only looked near the loose aligned optimum could not answer it.

The antipolar variants hold `phi2 = phi1 + 180` throughout the search rather than checking it
at the end, so they are exactly antipolar at every step. That is the same subspace
`antipolar_offsets` derives from the chain's own moment, and with four chains at each setting
angle the cell dipole cancels exactly.


### First, the cheapest possible measurement

Apply each pattern to the shipped cell and move nothing else. The number is the change in kcal/mol per monomer against the aligned registry at the *same* cell, so it isolates what the stagger does to the comonomer contacts from what relaxing the cell afterwards might win back.

| pattern | polar | antipolar |
|---|---|---|
| `aligned` | +0.0000 | +0.0000 |
| `antiphase` | +0.1054 | +0.0368 |
| `sheet` | +0.0402 | +0.0143 |
| `ladder` | +0.1282 | +0.0298 |
| `spread` | +0.1149 | +0.5315 |

**Every one of them is uphill.** So at the shipped shape the aligned registry is already the better arrangement for the nitriles, and the whole case for staggering rests on the cell being able to contract afterwards. The spread is small in absolute terms -- at most 0.531 kcal/mol per monomer, which is inside the 0.27 our own screen quotes as its resolution and about the size of the 0.135 your MACE+D3 run measured between the two polarities -- so read the signs, not the magnitudes.


### The staggered cells

| cell | a x b x c (A) | density | E/mon trunc. | E/mon Ewald | \|P\| (C/m2) | min offset |
|---|---|---|---|---|---|---|
| polar, **aligned (control)** | 20.920 x 9.607 x 30.756 | 1.6816 | -4.6082 | -6.9933 | 1.209e-01 | 0 |
| polar, antiphase | 20.920 x 9.607 x 30.756 | 1.6816 | -4.6082 | -6.9933 | 1.209e-01 | 0 |
| polar, sheet | 20.930 x 9.607 x 30.756 | 1.6809 | -4.5683 | -6.9405 | 1.209e-01 | 0 |
| polar, ladder | 20.938 x 9.609 x 30.756 | 1.6799 | -4.5236 | -6.9092 | 1.208e-01 | 3 |
| polar, spread | 20.951 x 9.608 x 30.756 | 1.6789 | -4.4944 | -6.8761 | 1.207e-01 | 2 |
| antipolar, **aligned (control)** | 22.653 x 9.490 x 30.756 | 1.5720 | -3.7514 | -6.1225 | 3.164e-15 | 0 |
| antipolar, antiphase | 22.653 x 9.490 x 30.756 | 1.5720 | -3.7514 | -6.1225 | 3.095e-15 | 0 |
| antipolar, sheet | 22.654 x 9.491 x 30.756 | 1.5719 | -3.7372 | -6.1014 | 3.271e-15 | 0 |
| antipolar, ladder | 22.665 x 9.493 x 30.756 | 1.5707 | -3.7234 | -6.0959 | 3.131e-15 | 3 |
| antipolar, spread | 22.678 x 9.494 x 30.756 | 1.5697 | -3.6957 | -6.0660 | 3.164e-15 | 2 |

Every cell is orthorhombic, 592 atoms, eight chains of 74, `C208H192F176N16`, c = 30.7557 A. `a x b` here is the **supercell**, twice the two-chain primitive's, so double the aligned `a` and `b` in *The two cells* above before comparing. The `aligned` rows are the shipped two-chain cells reproduced by the eight-chain search from scratch, and they come back to the reported figure, which is the control the rest of the table is read against.

**The antipolar cells' polarization was measured, not assumed**, for the reason the aligned report gives: the helper that builds antipolar cells has shipped two defects in this area. The largest magnitude across the antipolar rows above is 3.3e-15 C/m2. A stagger cannot change it -- translating a neutral chain moves its moment by `(sum q) d = 0`, which is also why polarity and stagger are independent choices here -- but it was checked rather than relied on, at every pattern.

One number to expect if you recompute that yourself: reading the polarization back off the extended-XYZ file gives about 2e-10 rather than 3e-15 C/m2, because the file carries eight decimals and the cancellation between the four chains up and the four down is exact only in the coordinates we computed it from. Both are zero against the polar cell's 0.12. The energies and densities do round-trip from the file: 1e-8 kcal/mol per monomer and exact to six decimals respectively.

**Releasing the per-chain slides.** Each chain was then allowed to slide off its integer monomer by up to half a monomer -- a degree of freedom the aligned cells could not have had, so reported separately rather than folded in:

| cell | E/mon with slides free | change | density | slides (monomers) |
|---|---|---|---|---|
| polar, antiphase | -4.6082 | +0.0000 | 1.6816 | +0.00 +0.00 +0.00 +0.00 +0.00 +0.00 +0.00 +0.00 |
| polar, sheet | -4.5683 | +0.0000 | 1.6809 | +0.00 +0.00 +0.00 +0.00 +0.00 +0.00 +0.00 +0.00 |
| polar, ladder | -4.5236 | +0.0000 | 1.6799 | +0.00 +0.00 +0.00 +0.00 +0.00 +0.00 +0.00 +0.00 |
| polar, spread | -4.4944 | -0.0001 | 1.6789 | +0.00 -0.01 -0.00 -0.01 +0.01 +0.00 +0.01 +0.00 |
| antipolar, antiphase | -3.7514 | +0.0000 | 1.5720 | +0.00 +0.00 +0.00 +0.00 +0.00 +0.00 +0.00 +0.00 |
| antipolar, sheet | -3.7372 | +0.0000 | 1.5719 | +0.00 +0.00 +0.00 +0.00 +0.00 +0.00 +0.00 +0.00 |
| antipolar, ladder | -3.7234 | -0.0000 | 1.5707 | +0.00 +0.00 +0.00 +0.00 +0.00 +0.00 +0.00 +0.00 |
| antipolar, spread | -3.6957 | -0.0001 | 1.5697 | +0.00 -0.01 -0.00 -0.01 +0.01 +0.00 +0.01 +0.00 |


### And does a staggered cell started dense stay dense?

The aligned deliverable answered its own density question with a fixed-shape probe: pin `(a, b)`, sample the setting angles and `dz` densely, and see whether the dense shapes are attractive at all. Here is the same probe run for every pattern, 40 samples per shape, polar branch, best kcal/mol per monomer found at each shape.

| shape a x b (A) | density | `aligned` | `antiphase` | `sheet` | `ladder` | `spread` |
|---|---|---|---|---|---|---|
| 10.46 x 4.80 | 1.6816 | -0.13 | +0.35 | +0.36 | -0.09 | -0.10 |
| 9.60 x 4.80 | 1.8336 | +7.12 | +7.03 | +7.04 | +7.13 | +7.19 |
| 9.00 x 4.80 | 1.9558 | +72.89 | +57.45 | +72.88 | +57.52 | +57.52 |
| 8.60 x 4.80 | 2.0468 | +73.35 | +73.26 | +73.31 | +73.33 | +73.32 |
| 8.20 x 4.80 | 2.1466 | +54.37 | +54.36 | +54.35 | +54.44 | +54.42 |
| 9.00 x 5.20 | 1.8054 | +3.87 | +3.89 | +3.85 | +3.87 | +3.88 |
| 8.60 x 5.60 | 1.7544 | +2.40 | +2.43 | +2.41 | +2.47 | +3.81 |
| 10.00 x 4.40 | 1.9202 | +12.69 | +12.47 | +12.53 | +12.49 | +12.46 |
| 7.80 x 5.00 | 2.1664 | +107.28 | +107.23 | +107.27 | +107.34 | +107.34 |

The pattern is the same one the aligned cells showed and the stagger does not change it: every shape at or above about 1.9 g/cm3 is **strongly repulsive whatever the stagger** -- tens of kcal/mol per monomer, not tenths -- and only the shapes already near the reported cell sample attractively. Staggering the comonomer units moves these numbers by a fraction of a kcal/mol against a wall tens of kcal/mol high. (5 of the 9 shapes are above 1.9.) Whatever is keeping this structure from the mixing rule's density, it is not the axial registry of the nitriles.


### Did breaking the registry recover the density?

Patterns excluded from this comparison because they relaxed back to the aligned cell *exactly*, `dz` having already reached that registry: `antiphase`. That is measured (same relaxed energy to 1e-6), not assumed, and it is the check that the eight-chain search reproduces the two-chain one.

| | polar | antipolar |
|---|---|---|
| genuine stagger patterns | `sheet`, `ladder`, `spread` | `sheet`, `ladder`, `spread` |
| aligned density | 1.6816 | 1.5720 |
| densest staggered | 1.6809 (`sheet`) | 1.5719 (`sheet`) |
| change in density | **-0.0007** | **-0.0001** |
| mixing-rule target | 2.028 | 2.028 |
| fraction of the deficit recovered | **-0.2%** | **-0.0%** |
| aligned E/mon | -4.6082 | -3.7514 |
| lowest staggered E/mon | -4.5683 (`sheet`) | -3.7372 (`sheet`) |
| change in E/mon | **+0.0399** | **+0.0141** |
| with the slides free | -4.5683 | -3.7372 |
| density recovered? | **no** | **no** |
| energy lowered? | **no** | **no** |

**It did not work, and that is worth saying plainly.** Breaking the axial
registry neither recovered the density nor lowered the energy, in either polarity. So the
aligned registry was **not** what was costing the 17%, and our own explanation of the deficit
was a guess that has now failed: the nitrile registry is not the binding constraint.

Two pieces of evidence in this section say why, and both were already visible in the aligned
report if we had read it properly. The fixed-geometry table shows every pattern is uphill at
the shipped shape, so the aligned registry was already the better arrangement for the
nitriles -- consistent with the aligned report's own observation that *"the nitriles point into
the wide inter-sheet gap and barely see each other"*. Something the nitriles barely see cannot
be what holds the cell open. And the fixed-shape probe shows the wall at high density is tens
of kcal/mol per monomer tall whatever the stagger, while the stagger moves the energy by
tenths. The registry was never the lever.

The leading remaining candidate is the **all-trans rigid-geometry constraint itself**: a chain
with three torsional states at 180 and +-60 degrees and fixed bond angles cannot make room for
a pendant nitrile the way a kinked chain can, and your converged endpoints do exactly that with
their five and fifteen out-of-state dihedrals. That is a hypothesis consistent with everything
we have, not a result -- testing it needs a relaxed-geometry search we have not run. What it
does mean is that the 17% is unlikely to be recoverable by any cell construction we can do on
top of an all-trans chain, so do not wait for a denser start from us.

**What to do with them anyway.** The spread between every cell in the table -- aligned and staggered, both polarities -- is at most 0.114 kcal/mol per monomer, against the 0.27 our own screen quotes as its resolution and the 0.135 your MACE+D3 run measured between the two polarities. Our ordering of these cells is therefore not information you should act on, and the reason to run the staggered ones is that they are a genuinely different starting registry whose relaxed endpoint we cannot predict -- the same argument that made the aligned pair worth running when your own seed failed on topology. They pass the same topology screen, they are looser rather than tighter than your rejected seed, and they break a symmetry the aligned cells could not. If they relax to the same basin as the aligned ones, that is a useful negative; if they do not, the 0.135 gap you measured was measured between two of several nearby basins rather than between the two phases.


### Topology of the staggered cells

This is the part of the handoff that has to be right, so it is **stricter** than the aligned report above: the same detected-graph-from-Cordero-radii test over all periodic images, at 7 covalent scale factors (1.05-1.35) rather than five, on all 10 relaxed cells including the controls.

| | polar aligned | polar antiphase | polar sheet | polar ladder | polar spread | antipolar aligned | antipolar antiphase | antipolar sheet | antipolar ladder | antipolar spread |
|---|---|---|---|---|---|---|---|---|---|---|
| atoms | 592 | 592 | 592 | 592 | 592 | 592 | 592 | 592 | 592 | 592 |
| components at 1.20 | 8 x 74 | 8 x 74 | 8 x 74 | 8 x 74 | 8 x 74 | 8 x 74 | 8 x 74 | 8 x 74 | 8 x 74 | 8 x 74 |
| lost bonds (1.05-1.35) | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| new bonds (1.05-1.35) | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| new **interchain** bonds (1.05-1.35) | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| worst intended bond d/(ri+rj) | 1.0187 | 1.0187 | 1.0187 | 1.0187 | 1.0187 | 1.0187 | 1.0187 | 1.0187 | 1.0187 | 1.0187 |
| closest non-bonded d/(ri+rj) | 1.5844 | 1.5844 | 1.5844 | 1.5844 | 1.5844 | 1.5844 | 1.5844 | 1.5844 | 1.5844 | 1.5844 |
| **safe scale window** | 1.019-1.584 | 1.019-1.584 | 1.019-1.584 | 1.019-1.584 | 1.019-1.584 | 1.019-1.584 | 1.019-1.584 | 1.019-1.584 | 1.019-1.584 | 1.019-1.584 |
| **min interchain d (A)** | **2.5655** | **2.5655** | **2.5599** | **2.5660** | **2.5643** | **2.5243** | **2.5243** | **2.5232** | **2.5254** | **2.5252** |
| min interchain d/(ri+rj) | 2.5152 | 2.5152 | 2.5097 | 2.5157 | 2.5140 | 2.4748 | 2.4748 | 2.4737 | 2.4759 | 2.4757 |
| verdict | **PASS** | **PASS** | **PASS** | **PASS** | **PASS** | **PASS** | **PASS** | **PASS** | **PASS** | **PASS** |

**Every covalent scale from 1.019 to 1.584 recovers the intended graph of all 10 cells exactly**: eight separate components of 74 atoms, zero bonds lost, zero invented, zero interchain. Both edges of that band are *intramolecular* -- the lower is the 1.09 A C-H bond against the Cordero radii, the upper a 1-3 backbone C...C -- so the verdict does not depend on the cutoff you pick, and staggering did not narrow it. The closest interchain contact anywhere in the set is **2.5232 A**, 2.47 times the sum of covalent radii, so a detector would have to be more than twice as generous as the most generous convention before it saw an interchain bond. A bare unscaled 1.0 is still outside the band and would lose every C-H bond, in these cells as in the aligned ones.


#### Minimum interchain distance by element pair (A)

| pair | polar aligned | polar antiphase | polar sheet | polar ladder | polar spread | antipolar aligned | antipolar antiphase | antipolar sheet | antipolar ladder | antipolar spread |
|---|---|---|---|---|---|---|---|---|---|---|
| H-N | 2.5655 | 2.5655 | 2.5599 | 2.5660 | 2.5643 | 2.5243 | 2.5243 | 2.5232 | 2.5254 | 2.5252 |
| F-N | 2.6268 | 2.6268 | 2.6250 | 2.6303 | 2.6318 | 2.7904 | 2.7904 | 2.7963 | 2.7918 | 2.8047 |
| C-H | 2.8167 | 2.8167 | 2.8152 | 2.8174 | 2.8167 | 2.7652 | 2.7652 | 2.7651 | 2.7666 | 2.7667 |
| F-H | 2.8323 | 2.8323 | 2.8313 | 2.8330 | 2.8325 | 2.7805 | 2.7805 | 2.7805 | 2.7819 | 2.7820 |
| C-N | 3.2028 | 3.2028 | 3.2018 | 3.2071 | 3.2096 | 3.4734 | 3.4734 | 3.4720 | 3.4745 | 3.4742 |
| C-F | 3.5754 | 3.5754 | 3.5718 | 3.5760 | 3.5748 | 3.5239 | 3.5239 | 3.5232 | 3.5252 | 3.5252 |
| C-C | 3.5886 | 3.5886 | 3.5847 | 3.5892 | 3.5879 | 3.5377 | 3.5377 | 3.5370 | 3.5391 | 3.5391 |
| N-N | 3.6083 | 3.6083 | 3.6174 | 3.6203 | 3.6259 | 4.2369 | 4.2369 | 4.2508 | 4.2356 | 4.2576 |
| F-F | 3.9008 | 3.9008 | 3.8986 | 3.9045 | 3.9061 | 3.7251 | 3.7251 | 3.7298 | 3.7272 | 3.7379 |
| H-H | 4.2171 | 4.2171 | 4.2147 | 4.2210 | 4.2227 | 4.0950 | 4.0950 | 4.0991 | 4.0972 | 4.1071 |


## What we would not stand behind

1. **The junction bonds have no fitted torsional parameters.** Four of 24 first-order
   terms, five of 24 pair terms and six of 24 triple terms are transferred from a
   homopolymer fitted in a different neighbourhood. They do not affect these all-trans
   coordinates at all, and they would affect any conformational ranking of this copolymer.

2. **The VDCN unit's own torsional preferences come from a homopolymer fit**, and that fit
   has a known reference-state problem: with rigid bond angles VDCN's all-trans chain
   carries about 176 kcal/mol of Lennard-Jones strain on a ten-bond oligomer (PVDF 11), so
   every rotation away from trans lowers the energy and the RIS convention of measuring from
   all-trans is measuring from a state the homopolymer would not occupy. With the angles
   free to relax that strain falls to about zero
   (`docs/VALENCE_FIT.md` section 4), which is why the copolymer is not obviously in the
   same trouble -- one nitrile in twelve monomers has no 1-3 nitrile neighbour to clash
   with. But we have not measured the copolymer's relaxed-angle strain, and we have not used
   the transferred model to rank anything here.

3. **Our charges are illustrative**, not fitted: the fitted PVDF potential
   (`pvdf-dft-valence-flux`) has no nitrile parameters, so it could not be used. Every
   energy, dipole and polarization above is therefore the
   rigid-ion value of an illustrative point-charge model, not a Berry-phase polarization,
   and carries no electronic contribution. Our own piezoelectric magnitudes are five to nine
   times short of measurement (`docs/BENCHMARK.md`).

4. **The antipolar/polar energy gap is a gap between two cells of an illustrative
   potential.** It is quoted with Ewald alongside the truncated sum for exactly that reason,
   and it is not evidence about which phase a real copolymer adopts.

5. **Eight chains in the four files at the top are four copies of two.** We used to write
   that this cost about 17% in density; it does not -- the staggered cells released exactly
   that constraint and came back no denser, and slightly higher in energy. The 17% deficit
   against the mixing rule is real, the tiling is not what causes it, and the best remaining
   candidate is the all-trans rigid-geometry constraint itself. See the correction in
   "Density" and the verdict in "Staggered variants".

5b. **`antipolar_cell_exact`'s own `dz` resolution was not good enough here** and a finer
   scan of the identical subspace found a lower cell. We took the lower one, but it means the
   antipolar branch of this search is only as converged as that finer scan, which is a
   0.5 A axial grid at one `(a, b)` rather than a global search of the subspace.

6. **Most importantly: this is a different basin from your accepted reference, not a better
   one.** Your accepted 8.33 mol% endpoint is not all-trans: of 192 mapped backbone
   dihedrals, 168 are T, eight are G+ and **sixteen are outside 30 degrees of any rotational
   isomeric state** -- one roughly +85-degree kink per symmetry-related chain. Our model
   cannot represent those sixteen at all. It has three states at 180, +60 and -60 degrees,
   and a rigid-geometry chain built from them; an 85-degree torsion is not in its vocabulary,
   and a start we build can never land in that basin. So these structures are **a second
   candidate for your electronic comparison**, constructed from the conformational side
   rather than found by relaxation, and they do not supersede your structure. If your kinked
   basin is lower under MACE+D3, that is a real finding about the copolymer and an equally
   real limit of the three-state rigid-geometry model, not a defect in these files.

7. **The backbone bond length is uniform at PVDF's 1.528 A.** VDCN's own entry carries the
   textbook 1.54 A, so the two bonds inside the VDCN unit are 0.8% short of what a
   VDCN-specific value would give them. The model carries one uniform backbone distance and
   the choice is declared rather than hidden.

8. **gamma is fixed at 90 degrees** and the chains are rigid. Your accepted cell is
   monoclinic (beta = 110.4 degrees). Releasing the cell is your full-cell stage's job.

9. **The staggered cells keep all eight chains at two setting angles and one `dz`.** They are
   relaxed over the same five variables the aligned cells were, which is what makes the
   comparison matched, but it means eight chains that are now axially independent are still
   orientationally a pair. A search with eight free setting angles is a bigger problem than
   this one and has not been run; if your full-cell stage rotates individual chains, that is a
   degree of freedom we did not explore rather than one we closed.

10. **The stagger is an integer number of monomers by construction.** The half-monomer
    release reported in that section is a local relaxation of the chosen pattern, not a search
    over fractional registries, and the patterns themselves were chosen from the neighbour-shell
    geometry rather than enumerated: with eight chains and twelve slots there are more
    arrangements than we tried, and `ladder` is optimal only for the measure stated (the minimum
    axial offset over close column pairs).

## Provenance

| | |
|---|---|
| repository | `polyfind`, branch `claude/polymeric-stable-arrangements-uh06b1` |
| generator, aligned | `examples/copolymer_starts.py` |
| generator, staggered | `examples/copolymer_staggered.py` |
| sequence | `polyfind.polymers.VDF_VDCN_11_1` |
| packer | `polyfind.pack.CrystalPacker` defaults (UFF LJ + DSF Coulomb, 8 A, illustrative charges) |
| eight-chain energy | `polyfind.supercell.SupercellEnergy`, checked against the packer |
| stagger construction | `polyfind.supercell.staggered_cell` / `stagger_pattern` |
| antipolar subspace | `polyfind.fitting.antipolar_offsets` / `antipolar_cell_exact` |
| topology check | `polyfind.topology.check_topology`, Cordero covalent radii |
| request | `docs/NOTE_VDCN_CONVERGENCE.md`, "Sarco response, 2026-09-10" |
| staggered request | `docs/REFERENCE_DATA_REQUEST.md`, "Consumer response, 2026-09-11" |
| scope | `docs/CHEMISTRY_EXTENSION.md` section 3, "Copolymer composition" |


Generated in 2626 s. Chain: 74 atoms, c = 30.7557 A, mass 782.448, 12 monomers per repeat, 8.33 mol% VDCN, uniform backbone bond 1.528 A. Staggered variants generated in 1752 s on top of that, 21141 eight-chain energy evaluations.
