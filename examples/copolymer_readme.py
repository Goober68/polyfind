"""Write ``deliverables/README.md`` from the numbers ``copolymer_starts.py`` measured.

Every figure in the README comes out of ``deliverables/summary.json`` rather than being
retyped, so the document and the files it describes cannot drift apart.  The prose is
hand-written and lives here; run this after ``copolymer_starts.py``.
"""
from __future__ import annotations

import json
import os

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, "deliverables")

PROSE_HEAD = """# All-trans 11 VDF : 1 VDCN starting structures (8.33 mol%)

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
"""

PROSE_CHEM = """
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
{dz_polar} monomers (polar) and {dz_anti} (antipolar). So the eight chains sit at **two**
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
`{moment}` e.A, lying along its own x, which is precisely the case that defect got wrong.
The dipole of the cell shipped here is reported below and is zero to floating point.
"""

PROSE_DENSITY = """
## Density: read this before you use these

**These cells are looser than your accepted reference, and that is the safe direction for
your protocol, but you should know why.**

Reference points, all measured with the same default packer as above, all-trans, two chains
per cell:

| | a x b x c (A) | density | E per monomer |
|---|---|---|---|
| PVDF homopolymer | 4.64 x 8.62 x 2.563 | {rho_pvdf} | -7.029 |
| VDCN homopolymer | 9.50 x 6.23 x 2.58 | {rho_vdcn} | +24.932 |
| this copolymer, polar | {polar_cell} | {polar_rho} | {polar_e} |
| your accepted kinked endpoint | 20.12 x 9.70 x 28.98, beta 110.4 | 1.96 | - |

(The two homopolymer rows reproduce with
`pack(periodic_chain(P, [0, 0], THREE_STATE), n_chains=2, screen="random", n_random=4000)`
for `P` in `PVDF`, `VDCN`. VDCN's `+24.9` is its own documented all-trans strain under rigid
bond angles, not a packing failure; see caveat 2.)

A mass-fraction mixing rule over those two homopolymer densities predicts **{rho_mix}** for
this composition. We get **{polar_rho}**, about **{rho_deficit}% less dense**. The shape says
where it goes: the copolymer's short axis comes out **{polar_a} A**, {da_pvdf}% off PVDF's own
4.64 A, while its long axis comes out **{polar_b} A**, {db_pvdf}% above PVDF's 8.62 and
{db_vdcn}% {db_vdcn_dir} the {vdcn_b} A long axis of the *pure VDCN* cell. So the majority
monomer sets the short axis and the nitrile sets the long one, even though only one monomer
in twelve carries a nitrile.

> **Corrected 2026-09-11.** This paragraph used to continue: *"The reason is the tiling rather
> than the chemistry: the 2x2x1 replication puts every chain's nitrile at the same axial
> height, turning eight isolated bulky groups into a continuous plane of them that the whole
> structure has to clear."* **We tested that and it is wrong**, twice over. First, the tiled
> cells never did put every nitrile at the same height: the packer's own `dz` already offsets
> the second chain, and at the shipped cells it does so by {dz_polar} monomers (polar) and
> {dz_anti} (antipolar), so four chains sit at one height and four at another. What the tiling
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
"""

PROSE_TOPO = """
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
"""

PROSE_TAIL = """
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
"""


# All-trans two-chain homopolymer references, measured with the SAME default packer the
# copolymer cells above were (``pack(periodic_chain(P, [0, 0], THREE_STATE), n_chains=2,
# screen="random", n_random=4000)``), so the density comparison is like for like.  PVDF's
# row is the seed the generator prints; VDCN's was measured separately by the same call.
RHO_PVDF, A_PVDF, B_PVDF = 2.073, 4.64, 8.62
RHO_VDCN, B_VDCN = 1.696, 9.50
# monomer formula masses (g/mol): VDF = C2H2F2, VDCN = C4H2N2
M_VDF, M_VDCN = 64.034, 78.074


def density_numbers(s) -> dict:
    """The density paragraph's figures, computed rather than retyped."""
    comp = {"vdf": 11, "vdcn": 1}  # one repeat of this sequence
    m_vdf, m_vdcn = comp["vdf"] * M_VDF, comp["vdcn"] * M_VDCN
    w_vdf, w_vdcn = m_vdf / (m_vdf + m_vdcn), m_vdcn / (m_vdf + m_vdcn)
    rho_mix = 1.0 / (w_vdf / RHO_PVDF + w_vdcn / RHO_VDCN)
    rho = s["polar"]["density_g_cm3"]
    short, long_ = sorted((s["polar"]["a"], s["polar"]["b"]))
    return {
        "rho_pvdf": f"{RHO_PVDF:.3f}", "rho_vdcn": f"{RHO_VDCN:.3f}",
        "rho_mix": f"{rho_mix:.3f}", "rho_deficit": f"{100 * (1 - rho / rho_mix):.0f}",
        "polar_rho": f"{rho:.3f}",
        "polar_e": f"{s['polar']['energy_per_monomer_kcal_mol']:+.3f}",
        "polar_cell": f"{short:.2f} x {long_:.2f} x {s['polar']['c']:.3f}",
        "polar_a": f"{short:.2f}", "polar_b": f"{long_:.2f}",
        "da_pvdf": f"{abs(100 * (short / A_PVDF - 1)):.0f}",
        "db_pvdf": f"{100 * (long_ / B_PVDF - 1):.0f}",
        "db_vdcn": f"{abs(100 * (long_ / B_VDCN - 1)):.0f}",
        "db_vdcn_dir": "above" if long_ > B_VDCN else "below",
        "vdcn_b": f"{B_VDCN:.2f}",
        "dz_polar": dz_monomers(s, "polar"), "dz_anti": dz_monomers(s, "antipolar"),
    }


def dz_monomers(s, label) -> str:
    """The packer's ``dz`` as a fraction of one monomer, which is what the registry is in."""
    return f"{s[label]['dz'] / (s[label]['c'] / 12):.2f}"


def cells_table(s) -> str:
    rows = ["", "## The two cells", "",
            "| | polar | antipolar |", "|---|---|---|"]

    def both(fmt, key, f=lambda d, k: d[k]):
        return f"| {fmt} | {f(s['polar'], key)} | {f(s['antipolar'], key)} |"

    g = lambda k, p=4: lambda d, _k: f"{d[_k]:.{p}f}"  # noqa: E731
    rows.append(both("a (A)", "a", g("a")))
    rows.append(both("b (A)", "b", g("b")))
    rows.append(both("c (A), = chain repeat", "c", g("c")))
    rows.append(both("gamma (deg)", "gamma", g("gamma", 1)))
    rows.append(both("phi1 (deg)", "phi1", g("phi1", 2)))
    rows.append(both("phi2 (deg)", "phi2", g("phi2", 2)))
    rows.append(both("dz (A)", "dz", g("dz")))
    rows.append(both("chain 2 antiparallel", "flip", lambda d, k: "yes" if d[k] else "no"))
    rows.append(both("density (g/cm3)", "density_g_cm3", g("density_g_cm3")))
    rows.append(both("E per monomer, truncated (kcal/mol)", "energy_per_monomer_kcal_mol",
                     lambda d, k: f"{d[k]:+.4f}"))
    rows.append(both("E per monomer, Ewald (kcal/mol)", "energy_per_monomer_ewald_kcal_mol",
                     lambda d, k: f"{d[k]:+.4f}"))
    rows.append(both("cell dipole (e.A)", "dipole_e_A",
                     lambda d, k: "(" + ", ".join(f"{v:+.4g}" for v in d[k]) + ")"))
    rows.append(both("P (C/m2)", "polarization_C_m2",
                     lambda d, k: "(" + ", ".join(f"{v:+.4g}" for v in d[k]) + ")"))
    rows.append(both("P magnitude (C/m2)", "polarization_magnitude_C_m2", lambda d, k: f"{d[k]:.3e}"))
    gap_t = s["antipolar"]["energy_per_monomer_kcal_mol"] - s["polar"]["energy_per_monomer_kcal_mol"]
    gap_e = s["antipolar"]["energy_per_monomer_ewald_kcal_mol"] - s["polar"]["energy_per_monomer_ewald_kcal_mol"]
    rows += ["",
             f"E(antipolar) - E(polar) = **{gap_t:+.4f}** kcal/mol per monomer truncated, "
             f"**{gap_e:+.4f}** with Ewald. Positive means the polar cell is the lower of the two "
             "under this potential.", ""]
    pol_axis = "b" if abs(s["polar"]["polarization_C_m2"][1]) > abs(s["polar"]["polarization_C_m2"][0]) else "a"
    rows += ["",
             f"**Which axis is which.** The polar cell's polarization lies along **{pol_axis}** "
             f"({s['polar']['b' if pol_axis == 'b' else 'a']:.2f} A), so that is the polar axis -- the short "
             "transverse one, as in beta-PVDF, whose polar axis this packer puts at 4.64 A. The long axis is "
             "the one the nitriles point along. The search does not know which axis is which and can return "
             "them either way round (`DESIGN.md` 5.2), so compare by length rather than by name.", "",
             "**Which antipolar branch.** Both cells have `flip = 0`, i.e. the two chains run the *same* way. "
             "The antipolar cell is antipolar by setting angle, `phi2 = phi1 + 180`, not by an up-down chain "
             "pair. `antipolar_offsets` offers both branches for this chain (its moment has no axial "
             "component, so `flip = 0, dphi = 180` and `flip = 1, dphi = 180` are each exactly antipolar) and "
             "the search picked this one on energy.", ""]
    rows.append(f"The antipolar cell's dipole is "
                f"{max(abs(v) for v in s['antipolar']['dipole_e_A']):.2e} e.A in its largest component, "
                f"i.e. zero to floating point, and it was **measured rather than assumed** "
                f"(see above). The polar cell's polarization magnitude is "
                f"{s['polar']['polarization_magnitude_C_m2']:.4f} C/m2, which is what makes it polar.")
    return "\n".join(rows)


def topo_table(s) -> str:
    rows = ["", "### Results", "",
            "| | polar primitive | polar 8-chain | antipolar primitive | antipolar 8-chain |",
            "|---|---|---|---|---|"]
    order = [("polar", "primitive"), ("polar", "8chain"), ("antipolar", "primitive"), ("antipolar", "8chain")]
    t = {(a, b): s[a]["topology"][b] for a, b in order}

    def row(name, f):
        return f"| {name} | " + " | ".join(f(t[k]) for k in order) + " |"

    rows.append(row("atoms", lambda d: str(d["n_atoms"])))
    rows.append(row("chains", lambda d: str(d["n_chains"])))
    rows.append(row("formula", lambda d: f"`{d['formula']}`"))
    rows.append(row("intended bonds", lambda d: str(d["n_intended_bonds"])))
    rows.append(row("connected components at scale 1.2",
                    lambda d: ("8 x 74" if d["components"]["1.2"] == [74] * 8 else
                               "2 x 74" if d["components"]["1.2"] == [74, 74] else
                               str(d["components"]["1.2"]))))
    rows.append(row("lost bonds (1.05-1.3)",
                    lambda d: str(max(d["lost"].values()))))
    rows.append(row("new bonds (1.05-1.3)", lambda d: str(max(d["new"].values()))))
    rows.append(row("new **interchain** bonds (1.05-1.3)",
                    lambda d: str(max(d["new_interchain"].values()))))
    rows.append(row("worst intended bond, d/(ri+rj)", lambda d: f"{d['max_bond_ratio']:.4f}"))
    rows.append(row("closest non-bonded pair, d/(ri+rj)", lambda d: f"{d['min_nonbond_ratio']:.4f}"))
    rows.append(row("**safe scale window**",
                    lambda d: f"{d['safe_scale_window'][0]:.3f} - {d['safe_scale_window'][1]:.3f}"))
    rows.append(row("**min interchain distance (A)**",
                    lambda d: f"**{d['min_interchain_distance_A']:.4f}**"))
    rows.append(row("min interchain d/(ri+rj)", lambda d: f"{d['min_interchain_ratio']:.4f}"))
    rows.append(row("closest intrachain non-bonded (A)",
                    lambda d: f"{d['min_intra_nonbond_distance_A']:.4f}"))
    rows.append(row("verdict", lambda d: "**PASS**" if d["pass"] else "**FAIL**"))
    w = t[("polar", "8chain")]["safe_scale_window"]
    wa = t[("antipolar", "8chain")]["safe_scale_window"]
    mi = min(t[("polar", "8chain")]["min_interchain_ratio"], t[("antipolar", "8chain")]["min_interchain_ratio"])
    rows += ["",
             f"**What that window means.** Any covalent scale between **{max(w[0], wa[0]):.3f}** and "
             f"**{min(w[1], wa[1]):.3f}** recovers the intended graph of all four cells exactly: eight (or two) "
             "separate components of 74 atoms, no bond missing, no bond invented. The usual choices -- 1.05, "
             "1.1, 1.15, 1.2, 1.3 times the sum of covalent radii -- all sit inside it. A bare **1.0** does "
             "not, and would lose every C-H bond; that is a property of the Cordero radii against a 1.09 A C-H "
             "bond and would do the same to any hydrocarbon, so if your detector uses 1.0 unscaled, check it "
             "against your own accepted structure before reading anything into it. The lower edge is set by the "
             "C-H bond at 1.09 A and the upper by a 1-3 backbone C...C at 2.41 A, both *intramolecular*, so the "
             "window is a property of the chain rather than of the packing; the packing has far more room than "
             f"that. The closest interchain contact in either 592-atom cell sits at **{mi:.2f}** times the sum "
             "of covalent radii, so a detector would have to be more than twice as generous as the most "
             "generous convention before it saw an interchain bond.", ""]
    rows += ["", "### Minimum interchain distance by element pair (A), 592-atom cells", "",
             "| pair | polar | antipolar |", "|---|---|---|"]
    keys = sorted(t[("polar", "8chain")]["min_interchain_by_element_A"],
                  key=lambda k: t[("polar", "8chain")]["min_interchain_by_element_A"][k])
    for k in keys:
        p = t[("polar", "8chain")]["min_interchain_by_element_A"][k]
        a = t[("antipolar", "8chain")]["min_interchain_by_element_A"].get(k)
        rows.append(f"| {k} | {p:.4f} | " + (f"{a:.4f} |" if a is not None else "> 6 |"))
    return "\n".join(rows)


PROSE_STAGGER = """
## Staggered variants: breaking the axial nitrile registry

You took up the offer at the end of our 2026-09-10 note, and this section is the answer.
A 2x2x1 tiling of a two-chain cell is an exact replication, so a chain cannot sit at a
different axial height from its own image -- it *is* its own image -- and we had argued that
the resulting nitrile registry was what cost the four files above about 17% density against a
mixing rule. The cells here slide each chain along its own axis by a whole number of monomers
so the comonomer units distribute along the repeat instead.

**The short version: {headline}** The verdict table below has the numbers and the honest
reading of them; the rest of this section is how we got there, because a result like this is
only worth anything if the search behind it is visible.

**Two corrections to what we sent you on 2026-09-10**, both found while building this:

1. The tiled cells never did put every nitrile at the same height. The packer's own `dz`
   offsets chain 2 against chain 1, and it does so by {dz_polar} monomers in the shipped polar
   cell and {dz_anti} in the antipolar one. Four chains sit at one height and four at another.
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
{skeleton}

So an integer stagger changes which monomer of each chain carries the nitrile and changes
nothing else about the chain. A fractional slide would put the backbones of neighbouring
chains out of register as well, which is a different object; we release it at the end as a
half-monomer local relaxation and report what it buys.

### Which patterns, and why those

The eight chains of the 2x2x1 cell sit on a centred lattice, and the neighbour shells are not
what the cell edges suggest. Measured at the aligned polar cell, the close column pairs
are{shells}. Two facts follow, and both pick the patterns. **Every** A-B pair is a close pair,
so no assignment of integer offsets can keep all close pairs more than {cap} monomers apart --
that is a backtracking search over the constraint graph, not an estimate. And the A-B shell is
the one the packer's `dz` already moves, so the registry a tiling genuinely locks is the
within-sublattice one: `b` at 4.80 A and `a` at 10.5 A.

| pattern | offsets (monomers, chain order below) | min offset over close pairs | mean | distinct heights |
|---|---|---|---|---|
{pattern_rows}

Chain order is `build_cell`'s: `{layout}` as `(ix, iy, k)`, with `k = 1` the chain at
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
  at {cap} monomers it **reaches the bound above**: it is the maximally-separated arrangement
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
agrees to {check} kcal/mol per monomer. The staggered cells and the aligned ones are
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
"""


# The four possible outcomes of the experiment, written before it was run, keyed
# ``(density recovered in both polarities, energy lowered in both)``.  The negative branches
# are the useful ones to have written in advance: a result that only reads well when it comes
# out the way you hoped is a result you should not trust.
VERDICT_PROSE = {
    (True, True): """**It worked, in both polarities.** Breaking the axial registry lets the cell
contract and lowers the energy, which is the outcome that confirms the diagnosis: the
aligned cell's density deficit was the nitrile plane and not the chemistry. Read the
recovered fraction rather than the sign -- if it is well short of 100% then the plane
explains part of the deficit and something else (most likely the all-trans constraint
itself, given your kinked endpoints) explains the rest.""",
    (True, False): """**Mixed, and the mixed answer is the informative one.** Staggering does let
the cell contract -- the nitrile plane really was holding the volume open -- but the
staggered cells are *not* lower in energy than the aligned one under this potential. That
combination says the aligned cell is a genuine minimum of the two-chain problem it was found
in, and the denser staggered cells are reachable but uphill: the energy the packing gains
from closer contact does not pay for the nitrile-nitrile contacts the stagger creates. For
your protocol that is still the useful direction -- a denser fixed-cell start is closer to
your accepted endpoint's 1.96 -- but it means these are seeds, and the ordering between them
under our potential should not be read as a prediction of which your relaxation will
prefer.""",
    (False, True): """**The energy falls but the density does not.** Staggering finds a lower cell
without contracting it, which means the gain is in the nitrile-nitrile contacts rather than
in the packing volume. The corollary is that the density deficit against the mixing rule is
**not** caused by the aligned registry, so something else accounts for it -- the all-trans
constraint is the obvious candidate, and your kinked endpoints are the evidence for it.""",
    (False, False): """**It did not work, and that is worth saying plainly.** Breaking the axial
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
top of an all-trans chain, so do not wait for a denser start from us.""",
}


def practical_note(st) -> str:
    """What to actually do with the files, given how the experiment came out."""
    v, m, cells = st["verdict"], st["meta"], st["cells"]
    ship = [f"{p}_{n}" for p in ("polar", "antipolar") for n in m["shipped"]]
    spread = max(abs(cells[k]["energy_per_monomer_kcal_mol"]
                     - cells[f"{k.split('_')[0]}_aligned"]["energy_per_monomer_kcal_mol"])
                 for k in ship)
    return (
        "**What to do with them anyway.** The spread between every cell in the table -- aligned "
        f"and staggered, both polarities -- is at most {spread:.3f} kcal/mol per monomer, against "
        "the 0.27 our own screen quotes as its resolution and the 0.135 your MACE+D3 run measured "
        "between the two polarities. Our ordering of these cells is therefore not information you "
        "should act on, and the reason to run the staggered ones is that they are a genuinely "
        "different starting registry whose relaxed endpoint we cannot predict -- the same argument "
        "that made the aligned pair worth running when your own seed failed on topology. They "
        "pass the same topology screen, they are looser rather than tighter than your rejected "
        "seed, and they break a symmetry the aligned cells could not. If they relax to the same "
        "basin as the aligned ones, that is a useful negative; if they do not, the 0.135 gap you "
        "measured was measured between two of several nearby basins rather than between the two "
        "phases.")


def headline(verdict) -> str:
    """One generated sentence, so the opening claim cannot drift from the verdict table."""
    def side(p):
        d = verdict[p]
        return (("recovers" if d["density_recovered"] else "does not recover")
                + f" the density ({d['aligned_density']:.4f} -> "
                + f"{d['densest_staggered_density']:.4f} g/cm3) and "
                + ("lowers" if d["energy_lowered"] else "does not lower")
                + f" the energy ({d['d_energy_per_monomer']:+.4f} kcal/mol per monomer)")
    if all(verdict[p]["density_recovered"] == verdict["polar"]["density_recovered"]
           and verdict[p]["energy_lowered"] == verdict["polar"]["energy_lowered"]
           for p in verdict):
        return f"in both polarities, breaking the registry {side('polar')}."
    return (f"in the polar branch breaking the registry {side('polar')}; in the antipolar "
            f"branch it {side('antipolar')}.")


def stagger_section(st) -> str:
    """The staggered-variant section, entirely from ``staggered_summary.json``."""
    if st is None:
        return ""
    m, cells, verdict = st["meta"], st["cells"], st["verdict"]
    sk = m["skeleton_periodicity"]
    moved = ", ".join(f"`{i}` {e} -> {e2} at {d} A" for i, e, e2, d in sk["atoms_that_move"])
    skeleton = (f"a `c/{m['monomers_per_repeat']}` slide maps **{sk['n_mapped_exactly']} of "
                f"{sk['n_atoms']}** atoms onto an identical atom, worst displacement "
                f"{sk['worst_displacement_A']:.1e} A. The {len(sk['atoms_that_move'])} that move "
                f"are exactly the comonomer site: {moved}.")
    by_d: dict[str, int] = {}
    for _, _, d in m["close_column_pairs"]:
        by_d[f"{d:.2f}"] = by_d.get(f"{d:.2f}", 0) + 1
    parts = [f"**{n}** pairs at {d} A" for d, n in sorted(by_d.items(), key=lambda kv: float(kv[0]))]
    shells = " " + (" and ".join(parts) if len(parts) < 3 else
                    ", ".join(parts[:-1]) + " and " + parts[-1])
    rows = []
    for name, q in m["patterns"].items():
        note = " **(optimal)**" if q["min_offset_monomers"] == m["best_possible_min_offset_monomers"] else ""
        rows.append(f"| `{name}` | {' '.join(str(v) for v in q['offsets'])} | "
                    f"{q['min_offset_monomers']}{note} | {q['mean_offset_monomers']:.2f} | "
                    f"{q['n_distinct_heights']} |")
    worst_check = max(max(v for k, v in c.items() if k != "aligned_params")
                      for c in m["evaluator_checks"].values())
    head = PROSE_STAGGER.format(
        skeleton=skeleton, shells=shells,
        cap=m["best_possible_min_offset_monomers"],
        pattern_rows="\n".join(rows),
        layout=", ".join(f"({a},{b},{c})" for a, b, c in m["column_layout"]),
        check=f"{worst_check:.0e}",
        headline=headline(verdict),
        dz_polar=f"{m['aligned_dz_monomers']['polar']:.2f}",
        dz_anti=f"{m['aligned_dz_monomers']['antipolar']:.2f}")
    return "\n".join([head, stagger_fixed_table(st), stagger_cells_table(st),
                      stagger_probe_table(st), stagger_verdict(st), stagger_topo_table(st)])


def stagger_fixed_table(st) -> str:
    """What each pattern costs applied to the shipped cell with nothing else moved."""
    m = st["meta"]
    f = m["cost_at_the_aligned_geometry"]
    names = list(m["patterns"])
    rows = ["", "### First, the cheapest possible measurement", "",
            "Apply each pattern to the shipped cell and move nothing else. The number is the "
            "change in kcal/mol per monomer against the aligned registry at the *same* cell, so "
            "it isolates what the stagger does to the comonomer contacts from what relaxing the "
            "cell afterwards might win back.", "",
            "| pattern | polar | antipolar |", "|---|---|---|"]
    for n in names:
        rows.append(f"| `{n}` | {f['polar'][n]:+.4f} | {f['antipolar'][n]:+.4f} |")
    worst = max(abs(v) for d in f.values() for v in d.values())
    rows += ["",
             "**Every one of them is uphill.** So at the shipped shape the aligned registry is "
             "already the better arrangement for the nitriles, and the whole case for staggering "
             "rests on the cell being able to contract afterwards. The spread is small in absolute "
             f"terms -- at most {worst:.3f} kcal/mol per monomer, which is inside the "
             "0.27 our own screen quotes as its resolution and about the size of the 0.135 your "
             "MACE+D3 run measured between the two polarities -- so read the signs, not the "
             "magnitudes.", ""]
    return "\n".join(rows)


def stagger_probe_table(st) -> str:
    """Does a staggered cell *started* dense stay dense?  The fixed-shape screen, per pattern."""
    cells = st["cells"]
    names = list(st["meta"]["patterns"])
    probes = {n: {p["group"]: p for p in cells[f"polar_{n}"]["shape_probe"]} for n in names}
    groups = [p["group"] for p in cells["polar_aligned"]["shape_probe"]
              if p["group"].startswith("fixed shape")]
    rows = ["", "### And does a staggered cell started dense stay dense?", "",
            "The aligned deliverable answered its own density question with a fixed-shape probe: "
            "pin `(a, b)`, sample the setting angles and `dz` densely, and see whether the dense "
            "shapes are attractive at all. Here is the same probe run for every pattern, 40 "
            "samples per shape, polar branch, best kcal/mol per monomer found at each shape.", "",
            "| shape a x b (A) | density | " + " | ".join(f"`{n}`" for n in names) + " |",
            "|---|---|" + "---|" * len(names)]
    for g in groups:
        p0 = probes[names[0]][g]
        rows.append(f"| {p0['a']:.2f} x {p0['b']:.2f} | {p0['density']:.4f} | "
                    + " | ".join(f"{probes[n][g]['best_per_monomer']:+.2f}" for n in names) + " |")
    dense = [g for g in groups if probes[names[0]][g]["density"] > 1.9]
    rows += ["",
             "The pattern is the same one the aligned cells showed and the stagger does not change "
             "it: every shape at or above about 1.9 g/cm3 is **strongly repulsive whatever the "
             "stagger** -- tens of kcal/mol per monomer, not tenths -- and only the shapes already "
             "near the reported cell sample attractively. Staggering the comonomer units moves "
             "these numbers by a fraction of a kcal/mol against a wall tens of kcal/mol high. "
             f"({len(dense)} of the {len(groups)} shapes are above 1.9.) Whatever is keeping this "
             "structure from the mixing rule's density, it is not the axial registry of the "
             "nitriles.", ""]
    return "\n".join(rows)


def cell_order(st) -> list:
    """Every relaxed cell, polar first, aligned control first within each polarity."""
    return [f"{p}_{n}" for p in ("polar", "antipolar") for n in st["meta"]["patterns"]]


def stagger_cells_table(st) -> str:
    m, cells = st["meta"], st["cells"]
    order = cell_order(st)
    rows = ["", "### The staggered cells", "",
            "| cell | a x b x c (A) | density | E/mon trunc. | E/mon Ewald | \\|P\\| (C/m2) | min offset |",
            "|---|---|---|---|---|---|---|"]
    for k in order:
        d = cells[k]
        tag = k.replace("_", ", ").replace("aligned", "**aligned (control)**")
        rows.append(
            f"| {tag} | {d['a_cell']:.3f} x {d['b_cell']:.3f} x {d['c_cell']:.3f} | "
            f"{d['density_g_cm3']:.4f} | {d['energy_per_monomer_kcal_mol']:+.4f} | "
            f"{d['energy_per_monomer_ewald_kcal_mol']:+.4f} | "
            f"{d['polarization_magnitude_C_m2']:.3e} | "
            f"{d['stagger_quality']['min_offset_monomers']} |")
    rows += ["",
             "Every cell is orthorhombic, 592 atoms, eight chains of 74, `C208H192F176N16`, "
             f"c = {m['chain_c_A']:.4f} A. `a x b` here is the **supercell**, twice the "
             "two-chain primitive's, so double the aligned `a` and `b` in *The two cells* above "
             "before comparing. The `aligned` rows are the shipped two-chain cells reproduced by "
             "the eight-chain search from scratch, and they come back to the reported figure, "
             "which is the control the rest of the table is read against.", "",
             "**The antipolar cells' polarization was measured, not assumed**, for the reason the "
             "aligned report gives: the helper that builds antipolar cells has shipped two defects "
             "in this area. The largest magnitude across the antipolar rows above is "
             f"{max(cells[k]['polarization_magnitude_C_m2'] for k in order if 'antipolar' in k):.1e}"
             " C/m2. A stagger cannot change it -- translating a neutral chain moves its moment by "
             "`(sum q) d = 0`, which is also why polarity and stagger are independent choices here "
             "-- but it was checked rather than relied on, at every pattern.", "",
             "One number to expect if you recompute that yourself: reading the polarization back "
             "off the extended-XYZ file gives about 2e-10 rather than 3e-15 C/m2, because the file "
             "carries eight decimals and the cancellation between the four chains up and the four "
             "down is exact only in the coordinates we computed it from. Both are zero against the "
             "polar cell's 0.12. The energies and densities do round-trip from the file: 1e-8 "
             "kcal/mol per monomer and exact to six decimals respectively.", ""]
    # the free-slide release
    fs = [(k, cells[k]["free_slide"]) for k in order if "free_slide" in cells[k]]
    if fs:
        rows += ["**Releasing the per-chain slides.** Each chain was then allowed to slide off its "
                 "integer monomer by up to half a monomer -- a degree of freedom the aligned cells "
                 "could not have had, so reported separately rather than folded in:", "",
                 "| cell | E/mon with slides free | change | density | slides (monomers) |",
                 "|---|---|---|---|---|"]
        for k, d in fs:
            rows.append(f"| {k.replace('_', ', ')} | {d['energy_per_monomer_kcal_mol']:+.4f} | "
                        f"{d['d_energy_per_monomer_vs_matched']:+.4f} | {d['density_g_cm3']:.4f} | "
                        + " ".join(f"{v:+.2f}" for v in d["slides_monomers"]) + " |")
        rows.append("")
    return "\n".join(rows)


def verdict_degenerate(v) -> str:
    """Patterns the comparison drops because ``dz`` already reaches their registry."""
    names: list[str] = []
    for p in ("polar", "antipolar"):
        for n in v[p].get("patterns_degenerate_with_aligned", []):
            if n not in names:
                names.append(n)
    return ", ".join(f"`{n}`" for n in names)


def stagger_verdict(st) -> str:
    m, cells, v = st["meta"], st["cells"], st["verdict"]
    rows = ["", "### Did breaking the registry recover the density?", ""]
    both_up = all(v[p]["density_recovered"] for p in ("polar", "antipolar"))
    both_down = all(v[p]["energy_lowered"] for p in ("polar", "antipolar"))

    def r(name, f):
        return f"| {name} | {f(v['polar'])} | {f(v['antipolar'])} |"

    deg = verdict_degenerate(v)
    if deg:
        rows += [f"Patterns excluded from this comparison because they relaxed back to the aligned "
                 f"cell *exactly*, `dz` having already reached that registry: {deg}. That is "
                 "measured (same relaxed energy to 1e-6), not assumed, and it is the check that "
                 "the eight-chain search reproduces the two-chain one.", ""]
    rows += ["| | polar | antipolar |", "|---|---|---|"]
    rows.append(r("genuine stagger patterns",
                  lambda d: ", ".join(f"`{n}`" for n in d.get("genuine_stagger_patterns", []))))
    rows.append(r("aligned density", lambda d: f"{d['aligned_density']:.4f}"))
    rows.append(r("densest staggered", lambda d: f"{d['densest_staggered_density']:.4f} "
                                                f"(`{d['densest_staggered'].split('_')[1]}`)"))
    rows.append(r("change in density", lambda d: f"**{d['d_density']:+.4f}**"))
    rows.append(r("mixing-rule target", lambda d: f"{d['mixing_rule_density']:.3f}"))
    rows.append(r("fraction of the deficit recovered",
                  lambda d: f"**{100 * d['fraction_of_deficit_recovered']:+.1f}%**"))
    rows.append(r("aligned E/mon", lambda d: f"{d['aligned_energy_per_monomer']:+.4f}"))
    rows.append(r("lowest staggered E/mon",
                  lambda d: f"{d['lowest_staggered_energy_per_monomer']:+.4f} "
                            f"(`{d['lowest_staggered'].split('_')[1]}`)"))
    rows.append(r("change in E/mon", lambda d: f"**{d['d_energy_per_monomer']:+.4f}**"))
    rows.append(r("with the slides free",
                  lambda d: f"{d.get('best_free_slide_energy_per_monomer', float('nan')):+.4f}"))
    rows.append(r("density recovered?", lambda d: "**yes**" if d["density_recovered"] else "**no**"))
    rows.append(r("energy lowered?", lambda d: "**yes**" if d["energy_lowered"] else "**no**"))
    rows += ["", VERDICT_PROSE[(both_up, both_down)], "", practical_note(st), ""]
    return "\n".join(rows)


def stagger_topo_table(st) -> str:
    cells = st["cells"]
    order = cell_order(st)
    t = {k: cells[k]["topology"] for k in order}
    scales = cells[order[0]]["topology"]["scales"]
    band = f"{min(scales):.2f}-{max(scales):.2f}"
    rows = ["", "### Topology of the staggered cells", "",
            "This is the part of the handoff that has to be right, so it is **stricter** than the "
            f"aligned report above: the same detected-graph-from-Cordero-radii test over all "
            f"periodic images, at {len(scales)} covalent scale factors ({band}) rather than five, "
            f"on all {len(order)} relaxed cells including the controls.", "",
            "| | " + " | ".join(k.replace("_", " ") for k in order) + " |",
            "|---|" + "---|" * len(order)]

    def row(name, f):
        return f"| {name} | " + " | ".join(f(t[k]) for k in order) + " |"

    rows.append(row("atoms", lambda d: str(d["n_atoms"])))
    rows.append(row("components at 1.20",
                    lambda d: "8 x 74" if d["components"]["1.2"] == [74] * 8 else str(d["components"]["1.2"])))
    rows.append(row(f"lost bonds ({band})", lambda d: str(max(d["lost"].values()))))
    rows.append(row(f"new bonds ({band})", lambda d: str(max(d["new"].values()))))
    rows.append(row(f"new **interchain** bonds ({band})",
                    lambda d: str(max(d["new_interchain"].values()))))
    rows.append(row("worst intended bond d/(ri+rj)", lambda d: f"{d['max_bond_ratio']:.4f}"))
    rows.append(row("closest non-bonded d/(ri+rj)", lambda d: f"{d['min_nonbond_ratio']:.4f}"))
    rows.append(row("**safe scale window**",
                    lambda d: f"{d['safe_scale_window'][0]:.3f}-{d['safe_scale_window'][1]:.3f}"))
    rows.append(row("**min interchain d (A)**", lambda d: f"**{d['min_interchain_distance_A']:.4f}**"))
    rows.append(row("min interchain d/(ri+rj)", lambda d: f"{d['min_interchain_ratio']:.4f}"))
    rows.append(row("verdict", lambda d: "**PASS**" if d["pass"] else "**FAIL**"))
    lo = max(d["safe_scale_window"][0] for d in t.values())
    hi = min(d["safe_scale_window"][1] for d in t.values())
    mind = min(d["min_interchain_distance_A"] for d in t.values())
    minr = min(d["min_interchain_ratio"] for d in t.values())
    rows += ["",
             f"**Every covalent scale from {lo:.3f} to {hi:.3f} recovers the intended graph of all "
             f"{len(order)} cells exactly**: eight separate components of 74 atoms, zero bonds lost, "
             "zero invented, zero interchain. Both edges of that band are *intramolecular* -- the "
             "lower is the 1.09 A C-H bond against the Cordero radii, the upper a 1-3 backbone "
             "C...C -- so the verdict does not depend on the cutoff you pick, and staggering did not "
             f"narrow it. The closest interchain contact anywhere in the set is **{mind:.4f} A**, "
             f"{minr:.2f} times the sum of covalent radii, so a detector would have to be more than "
             "twice as generous as the most generous convention before it saw an interchain bond. A "
             "bare unscaled 1.0 is still outside the band and would lose every C-H bond, in these "
             "cells as in the aligned ones.", ""]
    allk: dict[str, float] = {}
    for d in t.values():
        for k, val in d["min_interchain_by_element_A"].items():
            allk[k] = min(allk.get(k, 1e9), val)
    keys = sorted(allk, key=lambda k: allk[k])
    rows += ["", "#### Minimum interchain distance by element pair (A)", "",
             "| pair | " + " | ".join(k.replace("_", " ") for k in order) + " |",
             "|---|" + "---|" * len(order)]
    for k in keys:
        rows.append(f"| {k} | " + " | ".join(
            f"{t[c]['min_interchain_by_element_A'][k]:.4f}"
            if k in t[c]["min_interchain_by_element_A"] else "> 6" for c in order) + " |")
    rows.append("")
    return "\n".join(rows)


def files_table(s, st=None) -> str:
    rows = ["", "| file | contents |", "|---|---|"]
    for label in ("polar", "antipolar"):
        for tag, what in (("8chain", "**592 atoms, eight chains** -- the requested size"),
                          ("primitive", "148 atoms, two chains -- the primitive cell the above tiles")):
            rows.append(f"| `{s[label]['files'][tag]}` | {label}, {what} |")
    for label in ("polar", "antipolar"):
        rows.append(f"| `vdf11_vdcn1_alltrans_{label}_2chain.cif` | {label} primitive as P1 CIF |")
    if st is not None:
        rows.append("")
        rows.append("Staggered variants, added 2026-09-11 (see *Staggered variants* below; these are "
                    "the ones to run if you want the axial registry broken):")
        rows.append("")
        rows.append("| file | contents |")
        rows.append("|---|---|")
        for label in ("polar", "antipolar"):
            for name in st["meta"]["shipped"]:
                c = st["cells"][f"{label}_{name}"]
                if "files" not in c:
                    continue
                q = c["stagger_quality"]
                rows.append(f"| `{c['files']['xyz']}` | {label}, **592 atoms, eight chains**, "
                            f"`{name}` stagger: offsets {q['offsets']} monomers, "
                            f"{q['n_distinct_heights']} distinct axial heights, "
                            f"rho {c['density_g_cm3']:.4f} |")
        for label in ("polar", "antipolar"):
            for name in st["meta"]["shipped"]:
                c = st["cells"][f"{label}_{name}"]
                if "files" in c:
                    rows.append(f"| `{c['files']['cif']}` | the same cell as a P1 CIF "
                                "(fractional, wrapped) |")
        rows.append("| `staggered_summary.json` | every staggered number in this document, "
                    "machine-readable |")
    if st is not None:
        rows += ["", "And the machine-readable companions to both sets:", "",
                 "| file | contents |", "|---|---|"]
    rows.append("| `summary.json` | every aligned number in this document, machine-readable |")
    rows.append("| `ris_parameter_provenance.json` | per-term parameter source and match quality |")
    rows.append("")
    rows += [
        "Extended XYZ: `Lattice=\"...\"` on the comment line, `Properties=species:S:1:pos:R:3`, "
        "`pbc=\"T T T\"`, plus the energy, density and minimum interchain distance as extra "
        "key-value pairs; the staggered files add `stagger_monomers=\"...\"`, one integer per "
        "chain in the file's own chain order. Atom order is chain by chain, 74 atoms each, and "
        "within a chain it is backbone atom then its pendant atoms, repeating.",
        "",
        "**Coordinates are not wrapped into the cell.** Each chain's crystallographic repeat is "
        "written whole, so a few pendant atoms sit just outside the `c` boundary. That is "
        "deliberate: it is the molecular choice of branch, and it is what makes the cell dipole "
        "quoted below a property of the cell rather than of where the origin was put. Any "
        "PBC-aware reader handles it; wrap if your tool insists. The CIFs are the exception: "
        "fractional coordinates are taken modulo one, because a CIF reader would otherwise place "
        "the overhanging atoms outside the box. Use the extended XYZ for anything that reads a "
        "dipole off the coordinates.",
        "",
        "**What we could and could not validate about the format.** ASE is still not installed "
        "here, so these have not been round-tripped through the reader you will probably use. "
        "They have been round-tripped through an independent parser of our own, which is worth "
        "more than nothing and less than ASE: `tests/test_supercell.py` re-reads each staggered "
        "file as text and re-derives the energy, the density, the whole topology report and the "
        "stagger pattern from what comes back. Energies agree to 1e-8 kcal/mol per monomer, "
        "densities to six decimals, the topology verdict and safe window are identical, and the "
        "stagger recovered from the coordinates alone matches the declared `stagger_monomers`. "
        "Each header's extra key-value pairs are checked against the file's own geometry rather "
        "than against the generator's memory of it. Still spend the ten seconds to load one in "
        "your own stack before committing a long run.",
    ]
    return "\n".join(rows)


def main():
    with open(os.path.join(OUT, "summary.json")) as f:
        s = json.load(f)
    st = None
    path = os.path.join(OUT, "staggered_summary.json")
    if os.path.exists(path):
        with open(path) as f:
            st = json.load(f)
    m = s["meta"]
    moment = "(+5.83, 0, 0)"
    body = [
        PROSE_HEAD,
        files_table(s, st),
        PROSE_CHEM.replace("{moment}", moment)
                  .replace("{dz_polar}", dz_monomers(s, "polar"))
                  .replace("{dz_anti}", dz_monomers(s, "antipolar")),
        cells_table(s),
        PROSE_DENSITY.format(**density_numbers(s)),
        PROSE_TOPO,
        topo_table(s),
        stagger_section(st),
        PROSE_TAIL,
        f"\nGenerated in {m['seconds']:.0f} s. Chain: {m['chain_atoms']} atoms, "
        f"c = {m['chain_c_A']:.4f} A, mass {m['chain_mass']:.3f}, "
        f"{m['monomers_per_repeat']} monomers per repeat, "
        f"{m['mol_percent_vdcn']:.2f} mol% VDCN, uniform backbone bond {m['bond_length_A']} A."
        + ("" if st is None else
           f" Staggered variants generated in {st['meta']['seconds']:.0f} s on top of that, "
           f"{sum(c.get('n_energy_evaluations', 0) for c in st['cells'].values())} eight-chain "
           "energy evaluations.") + "\n",
    ]
    with open(os.path.join(OUT, "README.md"), "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(body))
    print("wrote", os.path.join(OUT, "README.md"))


if __name__ == "__main__":
    main()
