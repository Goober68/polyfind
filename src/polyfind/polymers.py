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
class Monomer:
    """One chemical monomer unit of a repeat: its backbone atoms, in chain order.

    This is the unit of *composition*.  A homopolymer has one of them and its
    :attr:`Polymer.backbone` is that monomer's backbone cycled, which is what the model
    always assumed.  A copolymer has an explicit list of them, concatenated once into a
    multi-monomer :attr:`Polymer.backbone` that is then cycled as a whole -- so an
    11 VDF : 1 VDCN repeat is a 24-bond explicit repeat and every ``k % B`` / ``i % B``
    index in the package goes on meaning what it meant (see :func:`copolymer`).

    ``source`` names the polymer whose *fitted* parameters describe this unit, which is
    normally the homopolymer of the same monomer.  It carries no parameters itself; it is
    the provenance label that :func:`polyfind.ris.transfer_ris` reports against.
    """

    name: str
    backbone: tuple[BackboneAtom, ...]
    source: str = ""

    def __post_init__(self) -> None:
        if not self.backbone:
            raise ValueError(f"monomer {self.name!r} has no backbone atoms")
        object.__setattr__(self, "backbone", tuple(self.backbone))
        if not self.source:
            object.__setattr__(self, "source", self.name)

    @property
    def n_bonds(self) -> int:
        return len(self.backbone)


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
    """A polymer's repeat unit: the backbone atoms that are cycled to build a chain.

    ``backbone`` is the **explicit periodic repeat**, however many monomers long.  For a
    homopolymer it is one monomer (PVDF: CH2, CF2 -- two entries, B = 2) and that is the
    degenerate case of what it always was.  For a copolymer it is the whole composition
    sequence written out -- an 11 VDF : 1 VDCN repeat is 12 monomers, i.e. 24 entries and
    B = 24 -- and ``sequence`` records which monomer each stretch came from.  Nothing else
    in the package changes, because every consumer already indexes the repeat modulo ``B``:
    ``chain.build_chain`` takes a backbone atom's chemistry from ``backbone[k % B]`` and
    :class:`polyfind.ris.RISModel` takes a bond's energies from ``t(i) = i % B``.  Build one
    with :func:`copolymer` rather than by hand, so that ``backbone`` and ``sequence`` cannot
    disagree.

    ``sequence`` defaults to a single :class:`Monomer` covering the whole backbone, which is
    exactly the homopolymer case, so ``monomers_per_repeat == 1`` and every count derived
    from it is unchanged for PE, PVDF, PVDC, CFE, CDFE, AN, VDCN and FANOME.

    ``bond_length`` stays ONE number for the whole repeat: the model assumes a uniform
    backbone-backbone distance, and a copolymer whose comonomers have different fitted C-C
    distances has to pick one (:func:`copolymer` makes that explicit rather than silent).
    """

    name: str
    backbone: tuple[BackboneAtom, ...]  # in chain order; the explicit periodic repeat
    bond_length: float  # backbone-backbone bond length (A), assumed uniform
    states: RISStates = THREE_STATE
    formula: str = ""
    sequence: tuple[Monomer, ...] = ()  # the monomers of the repeat, in chain order

    def __post_init__(self) -> None:
        object.__setattr__(self, "backbone", tuple(self.backbone))
        if not self.sequence:
            object.__setattr__(self, "sequence", (Monomer(self.name, self.backbone),))
            return
        object.__setattr__(self, "sequence", tuple(self.sequence))
        flat = tuple(a for m in self.sequence for a in m.backbone)
        if flat != self.backbone:
            detail = (f"{len(flat)} backbone atoms against backbone's {len(self.backbone)}"
                      if len(flat) != len(self.backbone) else
                      "the same number of atoms but not the same atoms")
            raise ValueError(f"polymer {self.name!r}: the monomer sequence concatenates to {detail}")

    @property
    def bonds_per_repeat(self) -> int:
        return len(self.backbone)

    @property
    def monomers_per_repeat(self) -> int:
        """Monomers in one periodic repeat; 1 for every homopolymer."""
        return len(self.sequence)

    @property
    def is_copolymer(self) -> bool:
        return len({m.name for m in self.sequence}) > 1

    def monomer_count(self, n_bonds: int) -> int:
        """Monomers in a chain stretch of ``n_bonds`` backbone bonds.

        ``n_bonds // bonds_per_repeat`` for a homopolymer, which is what every caller
        used to compute inline; for a copolymer it multiplies by the monomers in the
        repeat, so an "energy per monomer" stays an energy per *monomer* rather than
        per twelve of them.
        """
        if n_bonds % self.bonds_per_repeat:
            raise ValueError(f"{n_bonds} bonds is not a whole number of {self.bonds_per_repeat}-bond repeats")
        return n_bonds // self.bonds_per_repeat * self.monomers_per_repeat

    @property
    def monomer_of_atom(self) -> tuple[int, ...]:
        """``k -> index into sequence`` for each backbone atom of the repeat."""
        return tuple(i for i, m in enumerate(self.sequence) for _ in m.backbone)

    def composition(self) -> dict[str, int]:
        """Monomer name -> count in one repeat."""
        out: dict[str, int] = {}
        for m in self.sequence:
            out[m.name] = out.get(m.name, 0) + 1
        return out

    def mol_percent(self) -> dict[str, float]:
        out = self.composition()
        n = float(sum(out.values()))
        return {k: 100.0 * v / n for k, v in out.items()}

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
    # APPLIED, with the full re-measurement of the DESIGN section 5 tables: the DFT
    # Form I structure gives a C-C distance of 1.528 A, not the textbook 1.54 A this
    # used to carry.  See docs/REFERENCES.md, "Changes that were re-measured", for the
    # before/after on every quantity it moved.
    bond_length=1.528,
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
# electronegative.
#
# The backbone angles are *unequal*, and that is measured rather than assumed.  An
# earlier version of this entry set both to 114 deg, borrowing the justification
# PVDF's comment gives for its equal angles; for PVDC that justification does not
# hold.  The published structure (Takahagi, Chatani, Kusumoto & Tadokoro, Polym. J.
# 20, 883 (1988)) has C-CH2-C = 123 deg and C-CCl2-C = 114 deg, the wide one being
# exactly how the real chain relieves the Cl...Cl crowding that the strain note
# below measures.  An independent check agrees: fitted against DFT data with no
# crystallographic input, relax_backbone_angles recovers 123.1 deg at the CH2 carbon
# (docs/VALENCE_FIT.md section 4).
#
# Caveat, measured rather than assumed: with rigid angles the all-trans chain is
# strained, because a planar zigzag cannot relieve Cl...Cl contact when the angles
# cannot open.  Opening the CH2 angle to its measured 123 deg takes three quarters of
# that away -- the ten-bond all-trans Lennard-Jones strain falls from 313 kcal/mol at
# the old 114/114 to 74 at 123/114, against PVDF's 11 -- which is the audit's
# diagnosis confirmed at first hand.  All-trans is still not a minimum: rotating every
# bond together still lowers the energy, but by 4.0 kcal/mol rather than the 79 the
# equal angles gave, so the reference state is now nearly metastable rather than
# grossly unphysical.  It is still not the *crystal* conformation, which is a glide
# TGTG' form, so PVDC's fitted first-order energies remain measured from a state the
# chain does not adopt.  Letting the angles relax under a potential that has a bend
# term is the principled fix; see docs/CHEMISTRY_EXTENSION.md and
# docs/VALENCE_FIT.md section 4.
#
# CONSEQUENCE FOR PACKING, and it is a real limitation.  Unequal backbone angles give
# a planar chain a net 123 - 114 = 9 deg of rotation per two-bond repeat, so at ideal
# torsions the all-trans PVDC chain is a circular arc rather than a linear stem, with
# zero rise; ideal TG+TG- carries the same curl.  No PVDC sequence is therefore
# commensurate at ideal RIS angles and ``periodic_chain`` refuses all of them.  What
# does close is the chain with the small torsion deflections the real polymer has, and
# it lands on the published structure: asking which torsions close a TG+TG- repeat
# under these angles gives 175.3 and 49.4 deg with c = 4.677 A, against Takahagi's
# 175 deg, 49 deg and a 4.68 A fibre repeat.  Build it with
# ``periodic_chain_from_torsions``; see tests/test_pvdc.py.
PVDC = Polymer(
    name="pvdc",
    formula="-(CH2-CCl2)n-",
    bond_length=1.54,
    backbone=(
        BackboneAtom("C", "H", 1.09, 123.0, 108.0, -0.20, +0.10),  # CH2
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
# All-trans Lennard-Jones strain on a 10-bond oligomer, against PVDF's 11 kcal/mol and
# PVDC's 74:  AN 88, VDCN 176, FANOME 1.0e6.  AN and VDCN are strained but finite and
# comparable with CFE (162) and CDFE (199); FANOME's all-trans chain is not a physical
# structure at all -- methyl hydrogens of methoxy groups on consecutive substituted
# carbons come 0.80 A apart, and that is robust to the frozen rotamer (every C-O azimuth
# gives between 7e3 and 2e8).  Rotating away from all-trans lowers the energy for all
# three, so all-trans is NOT a sensible RIS reference for any of them and their fitted
# energies should not be read as physics.  Variable backbone angles, and for FANOME a
# reference state other than all-trans, are the principled fix; see
# docs/CHEMISTRY_EXTENSION.md.
#
# Re-measured with a potential that has angle terms and the angles free to relax
# (docs/VALENCE_FIT.md section 4): AN falls from 88 to 2 kcal/mol and VDCN from 176 to
# about zero, so for those two all-trans becomes a state the chain could occupy and the
# reference-state objection is answered.  It is not answered for FANOME, which still
# carries 4.8e3 kcal/mol with both angles pinned at the 135 deg relaxation bound; the
# verdict on it is unchanged and it stays registered and deliberately unfitted.
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

# --- Copolymer composition (docs/CHEMISTRY_EXTENSION.md section 3, first bullet).
def monomer_of(polymer: Polymer, name: str | None = None) -> Monomer:
    """The repeat of a homopolymer, as a :class:`Monomer` for use in :func:`copolymer`.

    ``source`` is set to the homopolymer's own name, so the provenance of every parameter
    transferred onto a copolymer bond points back at the fit it came from.
    """
    if polymer.monomers_per_repeat != 1:
        raise ValueError(f"{polymer.name!r} is not a single-monomer repeat; take its sequence instead")
    return Monomer(name=name or polymer.name, backbone=polymer.backbone, source=polymer.name)


def copolymer(name: str, units, bond_length: float | None = None, formula: str = "",
              states: RISStates = THREE_STATE) -> Polymer:
    """A :class:`Polymer` whose repeat is an explicit sequence of monomers.

    ``units`` is the composition in chain order: :class:`Monomer` objects, or
    :class:`Polymer` homopolymers (converted by :func:`monomer_of`).  The monomers'
    backbones are concatenated into one explicit repeat, so ``bonds_per_repeat`` is the
    total and every ``% B`` index in the package continues to work unchanged -- that is the
    whole of the generalisation, and it is why an 11:1 sequence needs no new mathematics.

    ``bond_length`` must be given whenever the units disagree about it, because the model
    carries one uniform backbone-backbone distance and picking one silently would hide a
    real approximation.  ``states`` must be shared; the RIS state set is a property of the
    chain, not of a monomer.
    """
    mons = [u if isinstance(u, Monomer) else monomer_of(u) for u in units]
    if not mons:
        raise ValueError("a copolymer needs at least one monomer")
    if bond_length is None:
        lengths = {}
        for u in units:
            if isinstance(u, Polymer):
                lengths.setdefault(round(u.bond_length, 6), u.name)
        if len(lengths) != 1:
            raise ValueError(
                f"copolymer {name!r}: the units disagree about the backbone bond length ({lengths}); "
                "this model carries one uniform value, so pass bond_length= explicitly and say which"
            )
        bond_length = float(next(iter(lengths)))
    return Polymer(
        name=name,
        backbone=tuple(a for m in mons for a in m.backbone),
        bond_length=float(bond_length),
        states=states,
        formula=formula,
        sequence=tuple(mons),
    )


VDF_UNIT = monomer_of(PVDF, "vdf")
VDCN_UNIT = monomer_of(VDCN, "vdcn")

# The 11 VDF : 1 VDCN periodic copolymer asked for in docs/NOTE_VDCN_CONVERGENCE.md:
# twelve monomers, 8.33 mol% VDCN, one explicit 24-bond repeat.  The VDCN unit sits at
# monomer index 6, i.e. backbone atoms 12 (its CH2) and 13 (its C(CN)2), as far from the
# repeat boundary as the sequence allows; the choice is arbitrary for an infinite periodic
# chain and only makes the bond indices below easier to read.
#
# The ONLY backbone atom that differs chemically from PVDF's is atom 13.  VDCN's own CH2
# entry is field-for-field identical to PVDF's (same element, same H pendants, 114/108 deg,
# -0.20/+0.10 e), so the substitution is a single-atom substitution in this model and the
# junction it creates is two backbone bonds wide.
#
# BOND LENGTH, stated rather than hidden: PVDF carries a measured 1.528 A C-C (DFT Form I,
# docs/REFERENCES.md) and VDCN carries the textbook 1.54 A.  The repeat takes 1.528 A --
# PVDF's, which is the measured one and describes 22 of the 24 bonds -- so the two bonds
# inside the VDCN unit are 0.8% short of what VDCN's own entry would give them.
VDF_VDCN_11_1 = copolymer(
    name="vdf11-vdcn1",
    formula="-(CH2-CF2)6-(CH2-C(CN)2)-(CH2-CF2)5-",
    units=[VDF_UNIT] * 6 + [VDCN_UNIT] + [VDF_UNIT] * 5,
    bond_length=PVDF.bond_length,
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
    "vdf11-vdcn1": VDF_VDCN_11_1,
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
