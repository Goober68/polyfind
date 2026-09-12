r"""Atomic polarizabilities and the induced-dipole model the lattice kernel can opt into.

The point charges of every potential in this package are fixed (or, with charge flux,
functions of the geometry) and never of the field: the dielectric constant of every crystal
here is exactly 1 by construction, which no material has.  ``docs/ELECTROMECHANICS.md``
section 5.7 measured what that costs -- the piezoelectric coefficients want electrostatics
several times stronger than the energies can tolerate -- and this module supplies the
missing response without touching the static charges: one isotropic point polarizability per
atom, induced dipoles solved self-consistently in the field of the charges, of each other and
of any applied field, interacting through the Ewald sum (:meth:`polyfind.ewald.Ewald.terms`
with ``dipoles``) and Thole-damped at short range (:class:`polyfind.ewald.Thole`).

**The energy.**  With ``alpha_i`` the polarizabilities, ``E0_i`` the permanent field (static
charges plus the applied field) and ``T`` the Ewald-summed, damped dipole interaction matrix,

    U(p) = sum_i |p_i|^2 / (2 alpha_i) - sum_i p_i . E0_i - (1/2) sum_ij p_i . T_ij . p_j

is minimised by ``(alpha^-1 - T) p = E0``, and at that minimum ``U = -(1/2) sum_i p_i . E0_i``:
*half* the dipole-field product, the other half having been spent polarising.  The factor
is the whole difference between a model that merely adds dipoles and one whose energy is
stationary in them, and :meth:`polyfind.pack.CrystalPacker.energy_and_grad` depends on that
stationarity -- by the Hellmann-Feynman argument the geometric gradient at the solution is
the partial derivative at fixed ``p``, which is why the solve is a direct linear solve and
not an iteration with a tolerance.

**The parameters are taken from the literature, not fitted here.**  :data:`VDS98` is Table 7
of van Duijnen and Swart, *Molecular and Atomic Polarizabilities: Thole's Model Revisited*,
J. Phys. Chem. A **102**, 2399-2407 (1998), doi:10.1021/jp980221f -- the exponential-damping
column, fitted to 52 experimental molecular polarizabilities with one damping constant
``a = 2.1304`` for every pair, converted here from atomic units (``1 bohr^3 = 0.148184711 A^3``).
It is the standard refit of Thole's 1981 model and the one that includes fluorine, which the
original did not.  The cubic form with ``a = 0.39`` is AMOEBA's choice (Ren and Ponder, J.
Phys. Chem. B **107**, 5933 (2003)) and is available for comparison; it was not refitted to
these polarizabilities and is not the default.  ``docs/REFERENCES.md`` entry 9 has the
provenance.

**Conventions.**  ``alpha`` in A^3 and dipoles in e.A, so with the kernel's Coulomb constant
``k = COULOMB / eps_r`` the field of a charge is ``k q / r^2`` and ``p = (alpha / k) E``; the
self energy is ``k |p|^2 / (2 alpha)``.  An applied field in V/A is ``EV_TO_KCAL`` times that
in the kernel's field unit.  The permanent field at a site excludes the chain's own 1-2 and
1-3 charges and scales its 1-4 ones by ``scale14``, exactly as the Coulomb energy does
(:func:`polyfind.ewald.charge_dipole_exclusion`); every pair of induced dipoles interacts.

Everything here is opt-in (``CrystalPacker(polarizable=Polarizable(), coulomb="ewald")``): the
default kernel has no induced dipoles and is bit-for-bit what it was.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .ewald import Thole

BOHR3_TO_A3 = 0.529177210903 ** 3  # 0.148184711 A^3

# van Duijnen & Swart 1998, Table 7, "exponential" column, atomic units, fitted to experiment.
_VDS98_AU = {"H": 2.7927, "C": 8.6959, "N": 6.5565, "O": 5.7494, "F": 2.6175,
             "S": 16.6984, "Cl": 16.1979, "Br": 23.5714, "I": 36.9880}
VDS98 = {e: round(v * BOHR3_TO_A3, 5) for e, v in _VDS98_AU.items()}  # A^3
VDS98_THOLE_A = 2.1304

# The same paper's "linear" column (Thole's original piecewise form), for the record only.
_VDS98_LINEAR_AU = {"H": 3.5020, "C": 10.1756, "N": 7.6048, "O": 6.3940, "F": 2.9413,
                    "S": 19.7422, "Cl": 16.1145, "Br": 22.6671, "I": 34.2434}
VDS98_LINEAR = {e: round(v * BOHR3_TO_A3, 5) for e, v in _VDS98_LINEAR_AU.items()}

TABLES = {"vds98": VDS98}


@dataclass(frozen=True)
class Polarizable:
    """What :class:`polyfind.pack.CrystalPacker` needs to carry induced dipoles.

    ``table``
        name of the polarizability table (``"vds98"``) or an explicit ``{element: A^3}``
        mapping; an explicit mapping may be partial and falls back to :data:`VDS98`.
    ``thole_a``, ``thole_form``
        the damping (:class:`polyfind.ewald.Thole`).  The defaults are the ones fitted
        *together with* the default table and should move with it.
    ``scale``
        multiplies every polarizability; a diagnostic knob for sensitivity, never a fit.
    ``exclude_bonded``
        whether the chain's own 1-2 / 1-3 charges are kept from polarising (and 1-4 scaled),
        as the Coulomb energy's exclusions are.  ``False`` lets every charge polarise every
        site and exists so that the effect of the choice can be measured.
    """

    table: str | dict = "vds98"
    thole_a: float = VDS98_THOLE_A
    thole_form: str = "exp"
    scale: float = 1.0
    exclude_bonded: bool = True

    def __post_init__(self) -> None:
        if isinstance(self.table, str) and self.table not in TABLES:
            raise ValueError(f"unknown polarizability table {self.table!r}; known: {sorted(TABLES)}")
        if not float(self.scale) > 0.0:
            raise ValueError(f"polarizability scale must be positive, got {self.scale!r}")
        Thole(self.thole_a, self.thole_form)  # validates

    def mapping(self) -> dict:
        if isinstance(self.table, str):
            return dict(TABLES[self.table])
        return {**VDS98, **{str(k): float(v) for k, v in self.table.items()}}

    def per_atom(self, elements) -> np.ndarray:
        """``alpha_i`` in A^3 for a list of element symbols, ``scale`` applied."""
        table = self.mapping()
        try:
            a = np.array([table[e] for e in elements], dtype=float)
        except KeyError as exc:
            raise KeyError(f"no polarizability for element {exc.args[0]!r}; pass it in "
                           "Polarizable(table={...})") from exc
        if np.any(a <= 0.0):
            raise ValueError("every polarizability must be positive")
        return a * float(self.scale)

    def thole(self) -> Thole:
        return Thole(float(self.thole_a), self.thole_form)

    def label(self) -> str:
        name = self.table if isinstance(self.table, str) else "custom"
        return (f"polarizable[{name}, {self.thole().label()}"
                + (f", x{self.scale:g}" if self.scale != 1.0 else "")
                + ("" if self.exclude_bonded else ", bonded field included") + "]")


# --- unit conversions for the dielectric constant --------------------------------------------
EPS0 = 8.8541878128e-12  # F/m
# (dmu/V) in e/A^2 per applied V/A  ->  dimensionless susceptibility  epsilon - 1
E_PER_A2_PER_V_PER_A_TO_CHI = 16.0218 / (EPS0 * 1e10)  # = 180.95 = 4 pi * 14.3996


def clausius_mossotti(alpha: float, volume: float) -> float:
    """Dielectric constant of a cubic lattice of one isotropic point polarizability per cell.

    ``(eps - 1) / (eps + 2) = 4 pi alpha / (3 V)`` exactly, for any Bravais lattice whose
    dipole lattice sum is the Lorentz value -- the cubic ones -- under tinfoil boundaries.
    The known answer for the self-consistent solve plus the Ewald dipole sum together.
    """
    x = 4.0 * np.pi * float(alpha) / (3.0 * float(volume))
    if x >= 1.0:
        raise ValueError("polarization catastrophe: 4 pi alpha / 3V >= 1")
    return (1.0 + 2.0 * x) / (1.0 - x)
