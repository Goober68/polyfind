"""Polymer definitions: rigid-geometry monomer templates and RIS state sets.

Geometry and charge parameters here are *illustrative* (UFF-like radii, textbook
bond lengths, modest point charges chosen to reproduce the CF2 dipole direction).
They are sufficient to exercise the pipeline and to reproduce known qualitative
behaviour; for production work the RIS energies should be re-fitted from a
machine-learned or ab initio potential via :func:`polyfind.forcefield.fit_ris`.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np


def pendant_pair(value, what: str) -> tuple:
    """One value shared by both pendants, or an explicit ``(first, second)`` pair.

    A bare scalar (or element symbol) means "both pendants are this", which is what
    every symmetric monomer -- CH2, CF2, CCl2 -- wants and what every existing call
    site passes positionally.  A 2-tuple/list gives the two pendants independently.
    Strings are scalars here, not sequences (``"Cl"`` is one element, not two).
    """
    if isinstance(value, (tuple, list)):
        if len(value) != 2:
            raise ValueError(f"{what} must be one value or a pair of two, got {value!r}")
        return (value[0], value[1])
    return (value, value)


# --------------------------------------------------------------- pendant fragments
@dataclass(frozen=True)
class PendantAtom:
    """One atom of a pendant group, at fixed local coordinates in the pendant frame.

    ``offset`` is ``(x, y, z)`` in angstrom in the orthonormal frame ``(e1, e2, e3)``
    that :func:`polyfind.chain.pendant_frames` builds at the backbone atom: ``e1`` points
    along the pendant bond, ``e2`` is perpendicular to it in the plane containing the
    backbone-angle bisector, and ``e3`` completes the frame with the handedness that
    makes the *second* pendant's frame the mirror image of the first's -- so two
    identical fragments give a locally mirror-symmetric backbone atom, exactly as two
    identical single atoms do.  Write a fragment's coordinates as if the frame were
    right-handed (which is what pendant 1 gets); a fragment with a mirror plane of its
    own, as nitrile and methoxy both have, is unaffected by landing on pendant 2.
    """

    element: str
    charge: float
    offset: tuple[float, float, float]


@dataclass(frozen=True)
class Pendant:
    """One pendant substituent: an ordered set of atoms placed rigidly on the pendant bond.

    A **single atom** is the degenerate case and is what every pre-phase-3 monomer has:
    one :class:`PendantAtom` at ``(bond, 0, 0)``, i.e. on the pendant bond axis at the
    bond length, which is exactly where :func:`polyfind.chain.substituent_positions` put
    it before fragments existed.  :func:`single_atom` builds that case, and
    :class:`BackboneAtom` builds it automatically from a plain element symbol, so the
    old positional spelling ``BackboneAtom("C", "F", 1.35, ...)`` is untouched.

    A **fragment** carries further atoms at fixed local coordinates.  Atom 0 is the one
    bonded to the backbone atom and must lie on the axis; the rest are rigid relative to
    the pendant bond direction, which is what "small rigid fragment" means here: a
    fragment adds atoms but no degrees of freedom.  Any internal torsion of the fragment
    is *frozen* at one representative rotamer, at the same level of approximation as the
    frozen backbone angles and bond lengths everywhere else in this model (see
    :func:`methoxy`, where the freezing is a real simplification rather than a symmetry).

    ``bonds`` lists the intra-fragment bonds as ``(i, j)`` local indices; the bond from
    the backbone atom to atom 0 is implicit.  They matter because
    :class:`polyfind.forcefield.SimpleFF` derives its nonbonded exclusions from the bond
    graph.
    """

    atoms: tuple[PendantAtom, ...]
    bonds: tuple[tuple[int, int], ...] = ()
    label: str = ""

    def __post_init__(self) -> None:
        if not self.atoms:
            raise ValueError("a pendant needs at least one atom")
        x, y, z = self.atoms[0].offset
        if x <= 0.0 or y != 0.0 or z != 0.0:
            raise ValueError(f"pendant atom 0 must sit on the pendant bond axis at (bond, 0, 0), got {self.atoms[0].offset!r}")
        n = len(self.atoms)
        for i, j in self.bonds:
            if not (0 <= i < n and 0 <= j < n) or i == j:
                raise ValueError(f"bad intra-fragment bond {(i, j)!r} for a {n}-atom pendant")

    def __len__(self) -> int:
        return len(self.atoms)

    @property
    def bond(self) -> float:
        """Backbone-to-pendant bond length (A): the distance to atom 0."""
        return float(self.atoms[0].offset[0])

    @property
    def charge(self) -> float:
        """Total charge of the pendant group (e).  For a single atom, its own charge."""
        return float(sum(a.charge for a in self.atoms))

    @property
    def elements(self) -> tuple[str, ...]:
        return tuple(a.element for a in self.atoms)


def single_atom(element: str, bond: float, charge: float) -> Pendant:
    """The one-atom pendant: the pre-phase-3 substituent, as a :class:`Pendant`."""
    return Pendant(atoms=(PendantAtom(element, charge, (float(bond), 0.0, 0.0)),), label=element)


def _place_atom(a, b, c, bond: float, angle_deg: float, dihedral_deg: float) -> np.ndarray:
    """NeRF placement for single points: |cd| = ``bond``, angle ``b-c-d``, dihedral ``a-b-c-d``.

    The same recursion as :func:`polyfind.chain.nerf`, restricted to one point and
    repeated here only so that :mod:`polyfind.polymers` need not import
    :mod:`polyfind.chain` (which imports this module).
    """
    a, b, c = (np.asarray(v, dtype=float) for v in (a, b, c))
    ang, dih = math.radians(angle_deg), math.radians(dihedral_deg)
    bc = (c - b) / np.linalg.norm(c - b)
    n = np.cross(b - a, bc)
    n = n / np.linalg.norm(n)
    m = np.cross(n, bc)
    return c + (-bond * math.cos(ang)) * bc + (bond * math.sin(ang) * math.cos(dih)) * m + (bond * math.sin(ang) * math.sin(dih)) * n


def nitrile(bond: float = 1.47, cn: float = 1.16, charges: tuple[float, float] = (+0.20, -0.35)) -> Pendant:
    """The -C#N pendant: a carbon and a nitrogen collinear with the pendant bond.

    Nitrile is the easy fragment and adds **no** degrees of freedom: the group is linear,
    so both atoms lie on the pendant bond axis and the fragment is rigid by construction
    rather than by approximation.  ``bond`` is the backbone C-C(N) distance and ``cn``
    the C#N triple bond.  Charges are illustrative (see the module docstring): the
    nitrogen carries the negative end of the C#N dipole.
    """
    return Pendant(
        atoms=(
            PendantAtom("C", charges[0], (bond, 0.0, 0.0)),
            PendantAtom("N", charges[1], (bond + cn, 0.0, 0.0)),
        ),
        bonds=((0, 1),),
        label="CN",
    )


def methoxy(
    bond: float = 1.41,
    oc: float = 1.43,
    coc: float = 111.5,
    ch: float = 1.09,
    och: float = 109.5,
    charges: tuple[float, float, float] = (-0.35, +0.15, +0.05),
) -> Pendant:
    """The -O-CH3 pendant, **frozen at one representative rotamer**.

    Methoxy is the hard fragment, and the simplification has to be stated rather than
    hidden.  Unlike nitrile it is *not* rigid: it has a real internal rotation about the
    backbone-C-O bond (which swings the methyl carbon around the pendant bond axis), and
    a second, nearly free three-fold rotation of the methyl about the O-C bond.  Neither
    is a degree of freedom this model can carry -- the whole rigid-geometry approximation
    is that only the backbone dihedrals vary -- so both are **frozen at one
    representative rotamer**, exactly as the backbone bond angles and bond lengths are
    frozen at one representative value:

    * the C-O torsion is frozen **staggered, with the methyl carbon anti to the backbone
      atom's other pendant**: the methyl sits at local ``-e2``, and because ``e2`` lies in
      the plane that also contains the other pendant's bond, the torsion
      (other pendant)-C-O-CH3 is then exactly 180 deg whichever of the two pendants the
      methoxy is.  Of the three staggered positions -- anti to the other pendant, anti to
      each backbone neighbour -- this is the one that is stateable without reference to a
      conformation, and on an all-trans chain it is also the least strained of the three
      by two orders of magnitude;
    * the methyl is frozen staggered with respect to the backbone carbon (dihedrals
      C(backbone)-O-C-H of 180, +60 and -60 deg).

    That is a choice, not a derivation.  A real treatment would give the C-O torsion its
    own rotational-isomeric-state dimension; this one trades that for a fragment that
    costs nothing in the existing machinery, and the cost is that the methoxy group's
    conformational entropy and its orientation-dependent sterics are simply absent.
    Charges are illustrative: ``(O, C, H)``, the three methyl hydrogens sharing the last.
    """
    o = np.array([bond, 0.0, 0.0])
    theta = math.radians(coc)
    c = o + oc * np.array([-math.cos(theta), -math.sin(theta), 0.0])  # staggered: methyl at -e2
    origin = np.zeros(3)  # the backbone atom, in pendant-frame coordinates
    q_o, q_c, q_h = charges
    atoms = [
        PendantAtom("O", q_o, tuple(float(v) for v in o)),
        PendantAtom("C", q_c, tuple(float(v) for v in c)),
    ]
    for psi in (180.0, 60.0, -60.0):  # staggered methyl, one H anti to the backbone atom
        h = _place_atom(origin, o, c, ch, och, psi)
        atoms.append(PendantAtom("H", q_h, tuple(float(v) for v in h)))
    return Pendant(atoms=tuple(atoms), bonds=((0, 1), (1, 2), (1, 3), (1, 4)), label="OCH3")


@dataclass(frozen=True)
class BackboneAtom:
    """One backbone atom of the repeat unit and its two pendant substituents.

    Each of ``substituent``, ``sub_bond`` and ``sub_charge`` is either a single value
    shared by both pendants (the symmetric case: CH2, CF2, CCl2 -- and the original
    meaning of the positional signature, kept unchanged) or an explicit
    ``(first, second)`` pair for a backbone atom bearing two *different* pendants
    (CFCl, CHCl).  Use the :attr:`substituents` / :attr:`sub_bonds` /
    :attr:`sub_charges` properties to read them as normalised 2-tuples; the raw
    fields are left as given so that consumers which only need the symmetric case
    (and the pair-aware :func:`polyfind.chain.substituent_positions`) can go on
    passing ``spec.sub_bond`` straight through.

    A pendant is either a **single atom**, named by its element symbol with its bond
    length in ``sub_bond`` and its charge in ``sub_charge`` -- exactly as before -- or a
    small rigid **fragment**, given as a :class:`Pendant` in place of the element symbol
    (``nitrile()``, ``methoxy()``).  A fragment carries its own bond length and charges,
    so the matching ``sub_bond`` / ``sub_charge`` entry may be ``None``; it is then filled
    in from the fragment, which is why ``spec.sub_bond`` stays a plain number (or pair of
    numbers) whatever the pendants are.  Giving a value explicitly is allowed but must
    agree with the fragment.  :attr:`pendants` is the general accessor and returns two
    :class:`Pendant` objects in every case, a single atom being the degenerate one.

    The three older accessors keep their meaning, generalised the obvious way:
    :attr:`substituents` gives the element of each pendant's *first* atom (the one bonded
    to the backbone), :attr:`sub_bonds` the backbone-pendant bond length, and
    :attr:`sub_charges` the *total* charge of each pendant group.  For single-atom
    pendants all three are literally what they were.

    ``sub_angle`` stays a single number: it is the angle *between* the two pendant bonds,
    which is one quantity however different the pendants are.  See
    :func:`polyfind.chain.substituent_positions` for how it is split between them and
    :func:`polyfind.chain.pendant_frames` for the frame the rest of a fragment hangs off.

    A backbone atom whose two pendants differ is a stereocentre
    (:attr:`is_stereocentre`); read :attr:`Polymer.is_chiral` before fitting an RIS
    model for such a polymer.
    """

    element: str
    substituent: str | Pendant | tuple  # pendant element(s) or fragment(s)
    sub_bond: float | None | tuple  # backbone-substituent bond length(s) (A); None from a fragment
    backbone_angle: float  # C(prev)-X-C(next) angle (deg)
    sub_angle: float  # substituent-X-substituent angle (deg)
    charge: float  # partial charge on the backbone atom (e)
    sub_charge: float | None | tuple  # total charge(s) of the substituent group(s) (e)

    def __post_init__(self) -> None:
        subs = pendant_pair(self.substituent, "substituent")
        bonds = pendant_pair(self.sub_bond, "sub_bond")
        charges = pendant_pair(self.sub_charge, "sub_charge")
        pendants = []
        for sub, bond, charge in zip(subs, bonds, charges):
            if isinstance(sub, Pendant):
                if bond is not None and not math.isclose(bond, sub.bond, rel_tol=0.0, abs_tol=1e-12):
                    raise ValueError(f"sub_bond {bond!r} contradicts pendant {sub.label!r} (bond {sub.bond})")
                if charge is not None and not math.isclose(charge, sub.charge, rel_tol=0.0, abs_tol=1e-12):
                    raise ValueError(f"sub_charge {charge!r} contradicts pendant {sub.label!r} (charge {sub.charge})")
                pendants.append(sub)
            elif bond is None or charge is None:
                raise ValueError(f"single-atom substituent {sub!r} needs its own sub_bond and sub_charge")
            else:
                pendants.append(single_atom(sub, bond, charge))
        object.__setattr__(self, "_pendants", tuple(pendants))
        # keep the raw fields numeric, so callers that forward ``spec.sub_bond`` straight
        # into ``substituent_positions`` (pack, linegroup) go on working unchanged
        b1, b2 = (p.bond for p in pendants)
        q1, q2 = (p.charge for p in pendants)
        object.__setattr__(self, "sub_bond", b1 if b1 == b2 else (b1, b2))
        object.__setattr__(self, "sub_charge", q1 if q1 == q2 else (q1, q2))

    @property
    def pendants(self) -> tuple[Pendant, Pendant]:
        """The two pendants, normalised; a single-atom pendant is the degenerate case."""
        return self._pendants  # type: ignore[attr-defined]

    @property
    def substituents(self) -> tuple[str, str]:
        """Element of each pendant's first atom (the one bonded to the backbone)."""
        return tuple(p.atoms[0].element for p in self.pendants)  # type: ignore[return-value]

    @property
    def sub_bonds(self) -> tuple[float, float]:
        return pendant_pair(self.sub_bond, "sub_bond")

    @property
    def sub_charges(self) -> tuple[float, float]:
        return pendant_pair(self.sub_charge, "sub_charge")

    @property
    def n_pendant_atoms(self) -> int:
        return sum(len(p) for p in self.pendants)

    @property
    def is_stereocentre(self) -> bool:
        """True when the two pendants differ, so the atom's configuration is a real
        degree of freedom (strictly a pseudo-asymmetric centre in an infinite chain,
        but it is what fixes the chain's tacticity).

        The comparison is on the whole :class:`Pendant`, so it covers fragments as well
        as single atoms: two nitriles are the same pendant (VDCN is achiral), a nitrile
        and a methoxy are not (FANOME is chiral).  Two *identical* fragments are placed
        as mirror images of one another (see :class:`PendantAtom`), so they leave the
        backbone atom locally mirror-symmetric just as two identical atoms do."""
        p1, p2 = self.pendants
        return p1 != p2


@dataclass(frozen=True)
class RISStates:
    """Rotational isomeric states for the backbone dihedrals."""

    names: tuple[str, ...]
    angles: tuple[float, ...]  # ideal dihedral angle per state (deg, IUPAC: trans = 180)
    mirror: tuple[int, ...]  # index of the mirror-image state (G+ <-> G-)

    @property
    def n(self) -> int:
        return len(self.names)

    def index(self, name: str) -> int:
        return self.names.index(name)


THREE_STATE = RISStates(names=("T", "G+", "G-"), angles=(180.0, 60.0, -60.0), mirror=(0, 2, 1))


@dataclass(frozen=True)
class Polymer:
    name: str
    backbone: tuple[BackboneAtom, ...]  # in chain order; one chemical repeat unit
    bond_length: float  # backbone-backbone bond length (A), assumed uniform
    states: RISStates = THREE_STATE
    formula: str = ""

    @property
    def bonds_per_repeat(self) -> int:
        return len(self.backbone)

    @property
    def atoms_per_repeat(self) -> int:
        """Backbone atoms plus every pendant atom of one repeat unit.

        Still exactly two pendants per backbone atom, but a pendant is a group of one or
        more atoms (phase 3), so this is no longer ``3 * len(backbone)`` -- it is that
        only when every pendant is a single atom, which is the case for PE, PVDF, PVDC,
        CFE and CDFE.
        """
        return sum(1 + a.n_pendant_atoms for a in self.backbone)

    @property
    def is_chiral(self) -> bool:
        """True when any backbone atom bears two *different* pendants.

        Such an atom is a stereocentre, so the chain has a tacticity.  The builder
        in :mod:`polyfind.chain` places pendant 1 on a fixed side of the local frame
        ``(previous, this, next)`` at every backbone atom, which gives the same
        configuration at every stereocentre relative to the chain direction: the
        chain it builds is the **isotactic** one.  Syndiotactic and atactic chains
        are not expressible in this model (see ``docs/CHEMISTRY_EXTENSION.md``).

        The consequence for RIS fitting, and the reason this flag exists: for an
        achiral chain, reflecting a conformer maps ``phi -> -phi`` and gives back
        the *same* molecule, so ``E(G+) == E(G-)`` exactly and
        :func:`polyfind.forcefield.fit_ris`'s ``symmetrize=True`` (which averages
        the fitted energies with their G+/G- mirror image, and forces mirror-pair
        state angles to equal magnitude) is exact.  For a chiral chain the mirror
        image is the *enantiomeric* chain, a different molecule, so ``E(G+)`` and
        ``E(G-)`` genuinely differ -- that difference is exactly what makes an
        isotactic chain choose a one-handed helix.  Averaging it away is unsound
        here, and ``symmetrize="auto"`` (the default) does not.

        What *does* survive, exactly, is mirror composed with chain reversal:
        ``E(phi_1..phi_N) == E(-phi_N..-phi_1)``, because reversing the chain
        direction swaps ``prev`` and ``next`` and so flips every stereocentre's
        configuration back.  In fitted-model terms (measured, B = 2), that reads
        ``e1[b, s] == e1[1 - b, m(s)]``, ``e2[b, s, s'] == e2[b, m(s'), m(s)]`` and
        ``e3[b, x, y, z] == e3[1 - b, m(z), m(y), m(x)]`` -- a state mirror with a
        bond-type shift that drops by one at each order, and the term's own indices
        read backwards.  All three hold to grid noise for CFE and CDFE while plain
        mirroring is violated by tens of kcal/mol, and
        :func:`polyfind.forcefield.fit_ris` averages a chiral fit over them after
        measuring each.  For an achiral chain reversal is a symmetry on its own,
        which is why plain mirroring is exact there and is what PE, PVDF and PVDC
        get.

        Note that plain *reversal* is no more a symmetry of a chiral chain than
        plain mirroring is: reading an isotactic chain backwards flips every
        stereocentre relative to the chain direction, so ``E(phi_N..phi_1)`` equals
        ``E(-phi_1..-phi_N)``, the enantiomer's energy, and both differ from
        ``E(phi_1..phi_N)`` (measured at up to 750 kcal/mol on a 24-bond CFE
        oligomer, against 0 for PVDF).  That is why
        :func:`polyfind.helix.sequence_images` -- and hence
        :func:`polyfind.enumerate.enumerate_periodic` -- uses only shifts and the
        composition for a chiral repeat, keeping a sequence and its mirror as two
        candidates, and why :func:`polyfind.linegroup.torsion_pattern` does not
        offer a chiral chain a glide (an improper isometry, which no chiral object
        has), leaving it the screws and the free/penalty fallback.
        """
        return any(a.is_stereocentre for a in self.backbone)


# Illustrative point charges; each repeat unit is neutral.
PVDF = Polymer(
    name="pvdf",
    formula="-(CH2-CF2)n-",
    # Queued improvement, deliberately not applied yet: the DFT Form I structure gives a
    # C-C distance of 1.528 A, which puts the computed chain repeat at 2.563 A against an
    # experimental 2.56, where the 1.54 A used here gives 2.583.  Changing it shifts every
    # PVDF energy, so it belongs with the next full re-measurement of the documented
    # tables rather than on its own; the potential is unfitted and contributes 3-6% cell
    # errors, which dwarf this 0.9%.
    bond_length=1.54,
    backbone=(
        # Equal backbone angles at both carbons, and exactly trans, is not a compromise:
        # it is what the accepted beta-PVDF structure has.  A DFT study of Form I gives a
        # C-C-C angle of 114.4 deg at *both* backbone carbons and an internal rotation
        # angle of exactly 180 deg (Nakhmanson-style periodic calculations agree), against
        # 112 deg and 2.534 A for polyethylene.  Hasegawa's alternately-deflected zigzag,
        # in which the CF2 carbons alternate out of the plane to relieve F...F crowding
        # (the F-F distance along the chain is 2.56 A against a van der Waals 2.70 A), is a
        # proposal rather than a refinement result: the deflection would double the chain
        # repeat to 5.12 A, and the intermediate layer line that doubling requires is not
        # present in the diffraction pattern.  Hasegawa's own fit used a statistically
        # disordered structure at c = 2.56 A.
        #
        # Every sentence of the paragraph above is sourced in docs/REFERENCES.md
        # sections 2.2 and 2.3.  Note the spread there: 114.4 deg is one DFT value,
        # and determinations range over 112.1-114.4 deg, so the 114.0 used here sits
        # near the top of the range rather than at a single agreed number.
        BackboneAtom("C", "H", 1.09, 114.0, 108.0, -0.20, +0.10),  # CH2
        BackboneAtom("C", "F", 1.35, 114.0, 106.0, +0.40, -0.20),  # CF2
    ),
)

PE = Polymer(
    name="pe",
    formula="-(CH2)n-",
    bond_length=1.54,
    backbone=(BackboneAtom("C", "H", 1.09, 112.0, 108.0, -0.12, +0.06),),
)

# Poly(vinylidene chloride).  The chlorine analogue of PVDF and the only other
# chemistry of interest that the present monomer model can express exactly: two
# identical single-atom substituents per backbone atom.  Chlorine is bulkier than
# fluorine (UFF x_i 3.95 vs 3.36 A) and its C-Cl bond much longer (1.77 vs 1.35 A),
# so steric crowding rather than electrostatics dominates its conformational
# preferences.  Charges are illustrative and each backbone atom is neutral, as for
# PVDF; the C-Cl dipole is set smaller than C-F, chlorine being the less
# electronegative.  Backbone angles are kept equal for the reason given above.
#
# That last sentence borrows PVDF's justification, and for PVDC it does not hold.
# The published structure (Takahagi, Chatani, Kusumoto & Tadokoro, Polym. J. 20,
# 883 (1988)) has *unequal* backbone angles, C-CH2-C = 123 deg and C-CCl2-C = 114
# deg, the wide one being exactly how the real chain relieves the Cl...Cl crowding
# that the strain note below measures.  Applying it would shift every PVDC energy,
# so it is queued rather than changed; see docs/REFERENCES.md, "Changes that need
# re-measurement".
#
# Caveat, measured rather than assumed: with these rigid angles the all-trans
# chain carries about 313 kcal/mol of Lennard-Jones strain (PVDF: 9), because a
# planar zigzag cannot relieve Cl...Cl contact when the angles cannot open.
# Every rotation away from trans lowers the energy, so all-trans is not a
# sensible RIS reference for this chemistry and its fitted energies should not
# be trusted.  That is qualitatively right -- PVDC does not adopt the planar
# zigzag that PVDF's beta phase does -- but the magnitude is an artifact of
# frozen bond angles.  Making the backbone angles refinement variables is the
# principled fix; see docs/CHEMISTRY_EXTENSION.md.
PVDC = Polymer(
    name="pvdc",
    formula="-(CH2-CCl2)n-",
    bond_length=1.54,
    backbone=(
        BackboneAtom("C", "H", 1.09, 114.0, 108.0, -0.20, +0.10),  # CH2
        BackboneAtom("C", "Cl", 1.77, 114.0, 110.0, +0.20, -0.10),  # CCl2
    ),
)

# --- Phase 2: asymmetric single-atom substituents (see docs/CHEMISTRY_EXTENSION.md).
#
# CFE and CDFE each carry a backbone atom with two *different* pendants, which the
# monomer model could not express before.  Geometry and charges are illustrative and
# chosen exactly the way PVDF's and PVDC's were: textbook C-F (1.35 A) and C-Cl
# (1.77 A) bond lengths, an F-C-Cl / H-C-Cl angle interpolated between PVDF's F-C-F
# (106 deg) and PVDC's Cl-C-Cl (110 deg), backbone angles kept equal for the reason
# given at PVDF, and modest point charges that leave every backbone atom's group
# (the atom plus its two pendants) exactly neutral, with the C-Cl dipole smaller
# than C-F because chlorine is the less electronegative.  They are good enough to
# exercise the pipeline; production numbers should be re-fitted with fit_ris.
#
# BOTH ARE CHIRAL (``Polymer.is_chiral``).  The pendant *order* below picks which
# enantiomer of the isotactic chain gets built: swapping the two entries of a pair
# gives the other one, whose energy landscape is this one's mirror image (the same
# minima, at mirrored conformations), so the choice is arbitrary but has to be made.
# What is not arbitrary is that G+ and G- are no longer equivalent for these chains:
# fit them with ``fit_ris(..., symmetrize=False)``.
CFE = Polymer(
    name="cfe",
    formula="-(CH2-CFCl)n-",
    bond_length=1.54,
    backbone=(
        BackboneAtom("C", "H", 1.09, 114.0, 108.0, -0.20, +0.10),  # CH2
        BackboneAtom("C", ("F", "Cl"), (1.35, 1.77), 114.0, 108.0, +0.30, (-0.20, -0.10)),  # CFCl
    ),
)

CDFE = Polymer(
    name="cdfe",
    formula="-(CHCl-CF2)n-",
    bond_length=1.54,
    backbone=(
        BackboneAtom("C", ("H", "Cl"), (1.09, 1.77), 114.0, 109.0, 0.00, (+0.10, -0.10)),  # CHCl
        BackboneAtom("C", "F", 1.35, 114.0, 106.0, +0.40, -0.20),  # CF2
    ),
)

# --- Phase 3: multi-atom pendant groups (see docs/CHEMISTRY_EXTENSION.md).
#
# AN, VDCN and FANOME need pendants that are small molecular fragments rather than single
# atoms.  Geometry is textbook and rigid: C-C#N 1.47 A with a 1.16 A triple bond, exactly
# linear, so the nitrile adds atoms but no degrees of freedom; C-O 1.41 A, O-CH3 1.43 A,
# C-O-C 111.5 deg for the methoxy, whose two internal rotations are FROZEN at one
# representative rotamer (the C-O torsion staggered, with the methyl carbon anti to the
# backbone atom's other pendant; the methyl itself staggered against the backbone carbon).
# That freezing is a real approximation, not a symmetry -- see :func:`methoxy` -- and it
# is at the same level as the frozen backbone angles the rest of the model already
# assumes.
#
# Charges are ILLUSTRATIVE, chosen exactly the way PVDF's and PVDC's were: modest point
# charges that leave every backbone atom's group (the atom plus both pendant groups)
# exactly neutral, with the negative end of each dipole on the more electronegative atom
# (nitrile N, ether O).  The nitrile group carries -0.15 e in total and the methoxy
# -0.05 e; the backbone atom takes whatever balances them.  Production numbers should be
# re-fitted with fit_ris against a real potential.
#
# CAVEAT, measured rather than assumed, and it is the phase-1 finding made worse.  With
# bond angles frozen a planar zigzag cannot relieve contact between substituents on 1-3
# backbone atoms, and a nitrile reaches 2.63 A from the backbone against chlorine's 1.77.
# All-trans Lennard-Jones strain on a 10-bond oligomer, against PVDF's 9 kcal/mol and
# PVDC's 313:  AN 88, VDCN 176, FANOME 1.0e6.  AN and VDCN are strained but finite and
# comparable with CFE (162) and CDFE (199); FANOME's all-trans chain is not a physical
# structure at all -- methyl hydrogens of methoxy groups on consecutive substituted
# carbons come 0.80 A apart, and that is robust to the frozen rotamer (every C-O azimuth
# gives between 7e3 and 2e8).  Rotating away from all-trans lowers the energy for all
# three, so all-trans is NOT a sensible RIS reference for any of them and their fitted
# energies should not be read as physics.  Variable backbone angles, and for FANOME a
# reference state other than all-trans, are the principled fix; see
# docs/CHEMISTRY_EXTENSION.md.
NITRILE = nitrile()
METHOXY = methoxy()

# Polyacrylonitrile: one hydrogen and one nitrile on the same backbone carbon, so that
# carbon is a stereocentre and the chain is CHIRAL (fit with the chiral path; see
# Polymer.is_chiral).  The pendant order (H first) picks the enantiomer, arbitrarily.
AN = Polymer(
    name="an",
    formula="-(CH2-CH(CN))n-",
    bond_length=1.54,
    backbone=(
        BackboneAtom("C", "H", 1.09, 114.0, 108.0, -0.20, +0.10),  # CH2
        BackboneAtom("C", ("H", NITRILE), (1.09, None), 114.0, 109.5, +0.05, (+0.10, None)),  # CH(CN)
    ),
)

# Poly(vinylidene cyanide): two identical nitriles on one carbon, so it is ACHIRAL --
# the phase-3 counterpart of PVDF and PVDC, and the case that shows a multi-atom pendant
# does not by itself create a stereocentre.
VDCN = Polymer(
    name="vdcn",
    formula="-(CH2-C(CN)2)n-",
    bond_length=1.54,
    backbone=(
        BackboneAtom("C", "H", 1.09, 114.0, 108.0, -0.20, +0.10),  # CH2
        BackboneAtom("C", NITRILE, None, 114.0, 110.0, +0.30, None),  # C(CN)2
    ),
)

# A nitrile and a methoxy on one carbon: two different pendants, so CHIRAL like AN, and
# the only chemistry here whose geometry rests on a frozen internal rotation.
FANOME = Polymer(
    name="fanome",
    formula="-(CH2-C(CN)(OCH3))n-",
    bond_length=1.54,
    backbone=(
        BackboneAtom("C", "H", 1.09, 114.0, 108.0, -0.20, +0.10),  # CH2
        BackboneAtom("C", (NITRILE, METHOXY), None, 114.0, 108.0, +0.20, None),  # C(CN)(OCH3)
    ),
)

POLYMERS: dict[str, Polymer] = {
    "pvdf": PVDF,
    "pe": PE,
    "pvdc": PVDC,
    "cfe": CFE,
    "cdfe": CDFE,
    "an": AN,
    "vdcn": VDCN,
    "fanome": FANOME,
}


def get_polymer(name: str) -> Polymer:
    try:
        return POLYMERS[name.lower()]
    except KeyError as e:
        raise KeyError(f"unknown polymer {name!r}; known: {sorted(POLYMERS)}") from e


# --- Nonbonded parameters (UFF: Rappe et al., JACS 1992, x_i = LJ minimum distance, D_i in kcal/mol)
UFF_LJ: dict[str, tuple[float, float]] = {
    "C": (3.851, 0.105),  # C_3
    "H": (2.886, 0.044),  # H_
    "N": (3.660, 0.069),  # N_3, added with the phase-3 nitrile groups
    "O": (3.500, 0.060),  # O_3, added with the phase-3 methoxy group
    "F": (3.364, 0.050),  # F_
    "Cl": (3.947, 0.227),  # Cl
}


def lj_params(elements) -> tuple[np.ndarray, np.ndarray]:
    """Per-atom (x_i, D_i) arrays for a list of element symbols."""
    x = np.array([UFF_LJ[e][0] for e in elements])
    d = np.array([UFF_LJ[e][1] for e in elements])
    return x, d
