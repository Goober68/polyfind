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


@dataclass(frozen=True)
class BackboneAtom:
    """One backbone atom of the repeat unit and its two pendant substituents."""

    element: str
    substituent: str  # element of both pendant atoms (H or F here)
    sub_bond: float  # backbone-substituent bond length (A)
    backbone_angle: float  # C(prev)-X-C(next) angle (deg)
    sub_angle: float  # substituent-X-substituent angle (deg)
    charge: float  # partial charge on the backbone atom (e)
    sub_charge: float  # partial charge on each substituent (e)


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
        return 3 * len(self.backbone)


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

POLYMERS: dict[str, Polymer] = {"pvdf": PVDF, "pe": PE}


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
