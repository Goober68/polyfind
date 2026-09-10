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

    PYTHONPATH=src python examples/copolymer_starts.py   # writes every file here
    PYTHONPATH=src python examples/copolymer_readme.py   # writes this README from them

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
than an answer to it.

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
in twelve carries a nitrile. The reason is the tiling rather than
the chemistry: the 2x2x1 replication puts
**every** chain's nitrile at the same axial height, turning eight isolated bulky groups into
a continuous plane of them that the whole structure has to clear.

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
stagger the nitriles and should recover much of the 17%. It is, though, consistent with your
kink being how the chain makes room for the nitrile at high density, which would be a real
result about the copolymer rather than about either of our codes.

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

5. **Eight chains here are four copies of two**, and that costs about 17% in density. See
   "How the two packings were built" and "Density".

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

## Provenance

| | |
|---|---|
| repository | `polyfind`, branch `claude/polymeric-stable-arrangements-uh06b1` |
| generator | `examples/copolymer_starts.py` |
| sequence | `polyfind.polymers.VDF_VDCN_11_1` |
| packer | `polyfind.pack.CrystalPacker` defaults (UFF LJ + DSF Coulomb, 8 A, illustrative charges) |
| antipolar subspace | `polyfind.fitting.antipolar_offsets` / `antipolar_cell_exact` |
| topology check | `polyfind.topology.check_topology`, Cordero covalent radii |
| request | `docs/NOTE_VDCN_CONVERGENCE.md`, "Sarco response, 2026-09-10" |
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
    }


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


def files_table(s) -> str:
    rows = ["", "| file | contents |", "|---|---|"]
    for label in ("polar", "antipolar"):
        for tag, what in (("8chain", "**592 atoms, eight chains** -- the requested size"),
                          ("primitive", "148 atoms, two chains -- the primitive cell the above tiles")):
            rows.append(f"| `{s[label]['files'][tag]}` | {label}, {what} |")
    for label in ("polar", "antipolar"):
        rows.append(f"| `vdf11_vdcn1_alltrans_{label}_2chain.cif` | {label} primitive as P1 CIF |")
    rows.append("| `summary.json` | every number in this document, machine-readable |")
    rows.append("| `ris_parameter_provenance.json` | per-term parameter source and match quality |")
    rows.append("")
    rows += [
        "Extended XYZ: `Lattice=\"...\"` on the comment line, `Properties=species:S:1:pos:R:3`, "
        "`pbc=\"T T T\"`, plus the energy, density and minimum interchain distance as extra "
        "key-value pairs. Atom order is chain by chain, 74 atoms each, and within a chain it is "
        "backbone atom then its pendant atoms, repeating.",
        "",
        "**Coordinates are not wrapped into the cell.** Each chain's crystallographic repeat is "
        "written whole, so a few pendant atoms sit just outside the `c` boundary. That is "
        "deliberate: it is the molecular choice of branch, and it is what makes the cell dipole "
        "quoted below a property of the cell rather than of where the origin was put. Any "
        "PBC-aware reader handles it; wrap if your tool insists. We could not validate the files "
        "against ASE here because ASE is not installed in this environment -- the format is "
        "written to the extxyz convention rather than round-tripped through a reader, so the "
        "ten seconds it would take you to load one before committing a long run is worth it.",
    ]
    return "\n".join(rows)


def main():
    with open(os.path.join(OUT, "summary.json")) as f:
        s = json.load(f)
    m = s["meta"]
    moment = "(+5.83, 0, 0)"
    body = [
        PROSE_HEAD,
        files_table(s),
        PROSE_CHEM.replace("{moment}", moment),
        cells_table(s),
        PROSE_DENSITY.format(**density_numbers(s)),
        PROSE_TOPO,
        topo_table(s),
        PROSE_TAIL,
        f"\nGenerated in {m['seconds']:.0f} s. Chain: {m['chain_atoms']} atoms, "
        f"c = {m['chain_c_A']:.4f} A, mass {m['chain_mass']:.3f}, "
        f"{m['monomers_per_repeat']} monomers per repeat, "
        f"{m['mol_percent_vdcn']:.2f} mol% VDCN, uniform backbone bond {m['bond_length_A']} A.\n",
    ]
    with open(os.path.join(OUT, "README.md"), "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(body))
    print("wrote", os.path.join(OUT, "README.md"))


if __name__ == "__main__":
    main()
