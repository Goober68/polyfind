"""Polymer definitions: rigid-geometry monomer templates and RIS state sets.

Geometry and charge parameters here are *illustrative* (UFF-like radii, textbook
bond lengths, modest point charges chosen to reproduce the CF2 dipole direction).
They are sufficient to exercise the pipeline and to reproduce known qualitative
behaviour; for production work the RIS energies should be re-fitted from a
machine-learned or ab initio potential via :func:`polyfind.forcefield.fit_ris`.
"""
from __future__ import annotations

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

    ``sub_angle`` stays a single number: it is the angle *between* the two pendants,
    which is one quantity however different they are.  See
    :func:`polyfind.chain.substituent_positions` for how it is split between them.

    A backbone atom whose two pendants differ is a stereocentre
    (:attr:`is_stereocentre`); read :attr:`Polymer.is_chiral` before fitting an RIS
    model for such a polymer.
    """

    element: str
    substituent: str | tuple[str, str]  # pendant element(s) (H, F, Cl here)
    sub_bond: float | tuple[float, float]  # backbone-substituent bond length(s) (A)
    backbone_angle: float  # C(prev)-X-C(next) angle (deg)
    sub_angle: float  # substituent-X-substituent angle (deg)
    charge: float  # partial charge on the backbone atom (e)
    sub_charge: float | tuple[float, float]  # partial charge(s) on the substituent(s) (e)

    @property
    def substituents(self) -> tuple[str, str]:
        return pendant_pair(self.substituent, "substituent")

    @property
    def sub_bonds(self) -> tuple[float, float]:
        return pendant_pair(self.sub_bond, "sub_bond")

    @property
    def sub_charges(self) -> tuple[float, float]:
        return pendant_pair(self.sub_charge, "sub_charge")

    @property
    def is_stereocentre(self) -> bool:
        """True when the two pendants differ, so the atom's configuration is a real
        degree of freedom (strictly a pseudo-asymmetric centre in an infinite chain,
        but it is what fixes the chain's tacticity)."""
        e1, e2 = self.substituents
        b1, b2 = self.sub_bonds
        q1, q2 = self.sub_charges
        return e1 != e2 or b1 != b2 or q1 != q2


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
        # still exactly two pendants per backbone atom, identical or not
        return 3 * len(self.backbone)

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
    bond_length=1.54,
    backbone=(
        # Equal backbone angles: with rigid geometry and exact 180 deg dihedrals, unequal
        # angles (real PVDF: ~112 at CH2, ~116-118 at CF2) make the all-trans chain curve;
        # the real chain compensates by deflecting its dihedrals (~ +/-172 deg). That
        # deflection belongs to the continuous refinement stage, not the discrete RIS model.
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

POLYMERS: dict[str, Polymer] = {"pvdf": PVDF, "pe": PE, "pvdc": PVDC, "cfe": CFE, "cdfe": CDFE}


def get_polymer(name: str) -> Polymer:
    try:
        return POLYMERS[name.lower()]
    except KeyError as e:
        raise KeyError(f"unknown polymer {name!r}; known: {sorted(POLYMERS)}") from e


# --- Nonbonded parameters (UFF: Rappe et al., JACS 1992, x_i = LJ minimum distance, D_i in kcal/mol)
UFF_LJ: dict[str, tuple[float, float]] = {
    "C": (3.851, 0.105),
    "H": (2.886, 0.044),
    "F": (3.364, 0.050),
    "Cl": (3.947, 0.227),
}


def lj_params(elements) -> tuple[np.ndarray, np.ndarray]:
    """Per-atom (x_i, D_i) arrays for a list of element symbols."""
    x = np.array([UFF_LJ[e][0] for e in elements])
    d = np.array([UFF_LJ[e][1] for e in elements])
    return x, d
