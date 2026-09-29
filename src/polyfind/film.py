"""From crystal to film: an exact laminate of crystal lamellae and amorphous layers.

What is measured on PVDF is a *film* coefficient, and ``docs/ELECTROMECHANICS.md`` section
5.11 shows that at least one of them cannot be a crystal property: ``d_31 = +20`` pC/N through
a 336 GPa chain axis would need a chain-column ``e`` of 6 C/m^2.  This module supplies the
missing layer.  In a uniaxially drawn semicrystalline film the chains run along the draw
direction and the crystal lamellae stack along it, separated by amorphous layers, so the
simplest geometry that respects that is a **laminate whose layer normal is the chain axis**
(``z`` in the packer's frame).

For a laminate the homogenisation is exact, not a mixing rule.  With the layer normal ``z`` and
a field *tangential* to the layers (the poling axis is in-plane), the quantities continuous
across every interface are the in-plane strains ``(eps_xx, eps_yy, eps_xy)``, the normal stress
``sigma_zz`` and the tangential field; the ones that differ from layer to layer are the in-plane
stresses and ``eps_zz``.  Each layer's strain-charge law ``eps = S sigma + d^T E`` is rewritten
with the continuous quantities as inputs,

    (sigma_p, eps_n) = M_k (eps_p, sigma_n, E),

those ``M_k`` are volume-averaged (the outputs are what averages), and the average is turned
back into a film ``S`` and ``d``.  The transverse shears ``sigma_yz``, ``sigma_xz`` are also
continuous but are outside what :mod:`polyfind.mechanics` can express; for an orthorhombic
crystal with the chain along ``z`` they do not couple to the four normal-and-in-plane-shear
components, so leaving them out changes nothing here (``Film.notes`` says so).

What the model leaves out, and says so: the amorphous phase is non-piezoelectric and
non-electrostrictive (no trapped charge, no interface polarization); crystallites are perfectly
aligned (polar axes along the poling direction, chains along the draw); the geometry is one
lamellar stack, not an aggregate; no field along the layer normal (that would need the phases'
permittivities, since ``D_z`` and not ``E_z`` is continuous).  Every one of those is a
reason a measured film can differ from this number, which is the point of stating them.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

# The four strain components mechanics.py expresses on the deformable path, in its order.
COMPONENTS = (0, 1, 2, 5)  # Voigt xx, yy, zz, xy
_P = [0, 1, 3]  # in-plane (continuous strain): xx, yy, xy -- positions in COMPONENTS
_N = [2]  # normal to the layers (continuous stress): zz
_TANGENTIAL = [0, 1]  # field components tangential to the layers: x, y
PC_PER_N_TO_PER_GV_PER_M = 1e-3  # d in pC/N = pm/V; times a field in GV/m gives strain


def isotropic_compliance(young_gpa: float, poisson: float) -> np.ndarray:
    """(4, 4) compliance (1/GPa) of an isotropic solid in the order xx, yy, zz, xy (engineering shear)."""
    if not young_gpa > 0.0:
        raise ValueError("Young's modulus must be positive")
    if not -1.0 < poisson < 0.5:
        raise ValueError("Poisson's ratio must lie in (-1, 0.5) for a stable isotropic solid")
    E, v = float(young_gpa), float(poisson)
    S = np.array([[1.0, -v, -v, 0.0], [-v, 1.0, -v, 0.0], [-v, -v, 1.0, 0.0], [0.0, 0.0, 0.0, 2.0 * (1.0 + v)]])
    return S / E


def _mixed(S: np.ndarray, d: np.ndarray) -> np.ndarray:
    """A layer's law with the continuous quantities as inputs.

    ``S`` (4, 4) in 1/GPa, ``d`` (2, 4) the tangential rows in strain per GV/m.  Returns the
    (4, 6) matrix taking ``(eps_p[3], sigma_n[1], E_t[2])`` to ``(sigma_p[3], eps_n[1])``.
    """
    Spp, Spn = S[np.ix_(_P, _P)], S[np.ix_(_P, _N)]
    Snp, Snn = S[np.ix_(_N, _P)], S[np.ix_(_N, _N)]
    dp, dn = d[:, _P].T, d[:, _N].T  # (3, 2), (1, 2): strain per unit field
    A = np.linalg.inv(Spp)
    top = np.hstack([A, -A @ Spn, -A @ dp])  # sigma_p
    bot = np.hstack([Snp @ A, Snn - Snp @ A @ Spn, dn - Snp @ A @ dp])  # eps_n
    return np.vstack([top, bot])


@dataclass
class Film:
    """The homogenised film: compliance, stiffness and tangential-field ``d``.

    ``S``, ``C`` (4, 4) in 1/GPa and GPa over ``COMPONENTS``; ``d`` (3, 4) in pC/N, rows the
    field direction ``x, y, z`` -- the ``z`` row (field along the layer normal) is ``nan``, see
    the module docstring.  ``crystallinity`` is the crystal volume fraction.  ``young_draw`` is
    the film's Young's modulus along the chain (draw) axis, a quantity a tensile test measures
    and so a check on the inputs.
    """

    S: np.ndarray
    C: np.ndarray
    d: np.ndarray
    crystallinity: float
    amorphous_young: float
    amorphous_poisson: float
    notes: list = field(default_factory=list)

    @property
    def young_draw(self) -> float:
        return 1.0 / float(self.S[2, 2])

    def film_coefficients(self, polar: int = 0) -> dict:
        """``d_33``, ``d_32``, ``d_31`` in the film convention: 3 = poling = ``polar`` (x or y),
        1 = draw = the chain axis ``z``, 2 = the remaining in-plane axis.  Signs as ``d`` has them;
        quote with the poling axis along +P by passing a ``d`` already oriented that way."""
        if polar not in _TANGENTIAL:
            raise ValueError("the poling axis must be tangential to the layers (x or y)")
        other = 1 - polar
        row = self.d[polar]
        return {"d33": float(row[polar]), "d32": float(row[other]), "d31": float(row[2])}


def laminate(crystal_S, crystal_d, crystallinity: float, amorphous_young: float,
             amorphous_poisson: float) -> Film:
    """Homogenise crystal lamellae and amorphous layers stacked along the chain axis ``z``.

    ``crystal_S`` (4, 4) is the crystal compliance over :data:`COMPONENTS` (1/GPa), e.g.
    ``Elastic.S`` from the deformable path; ``crystal_d`` (3, 4) its converse coefficients in
    pC/N (``Piezoelectric.d_from_e``), rows ``x, y, z``.  The amorphous layers are isotropic,
    elastic and electrically passive.  ``crystallinity`` is the crystal volume fraction.
    """
    Sc = np.asarray(crystal_S, dtype=float)
    dc = np.asarray(crystal_d, dtype=float)
    if Sc.shape != (4, 4) or dc.shape != (3, 4):
        raise ValueError("crystal_S must be (4, 4) and crystal_d (3, 4), over xx, yy, zz, xy")
    phi = float(crystallinity)
    if not 0.0 <= phi <= 1.0:
        raise ValueError("crystallinity is a volume fraction in [0, 1]")
    Sa = isotropic_compliance(amorphous_young, amorphous_poisson)
    da = np.zeros((2, 4))
    Mbar = (phi * _mixed(Sc, dc[_TANGENTIAL] * PC_PER_N_TO_PER_GV_PER_M)
            + (1.0 - phi) * _mixed(Sa, da))
    # back to strain form: sigma_p = A eps_p + B sigma_n + G E ; eps_n = D eps_p + F sigma_n + H E
    A, B, G = Mbar[:3, :3], Mbar[:3, 3:4], Mbar[:3, 4:]
    D, F, H = Mbar[3:, :3], Mbar[3:, 3:4], Mbar[3:, 4:]
    Ai = np.linalg.inv(A)
    # eps_p = Ai sigma_p - Ai B sigma_n - Ai G E ; eps_n = D Ai sigma_p + (F - D Ai B) sigma_n + (H - D Ai G) E
    S = np.zeros((4, 4))
    S[np.ix_(_P, _P)] = Ai
    S[np.ix_(_P, _N)] = -Ai @ B
    S[np.ix_(_N, _P)] = D @ Ai
    S[np.ix_(_N, _N)] = F - D @ Ai @ B
    dT = np.zeros((4, 2))  # strain per unit tangential field, per GV/m
    dT[_P] = -Ai @ G
    dT[_N] = H - D @ Ai @ G
    d = np.full((3, 4), np.nan)
    d[_TANGENTIAL] = dT.T / PC_PER_N_TO_PER_GV_PER_M
    S = 0.5 * (S + S.T)  # symmetric to rounding; the asymmetry is ~1e-17
    notes = [
        "laminate normal = chain (draw) axis; exact homogenisation for a tangential field",
        "amorphous phase isotropic, elastic, non-piezoelectric, non-electrostrictive",
        "crystallites perfectly aligned; one lamellar stack, not an aggregate",
        "eps_yz, eps_xz decouple for an orthorhombic crystal with the chain along z and are omitted",
        "no d for a field along the layer normal (needs the phases' permittivities)",
    ]
    return Film(S=S, C=np.linalg.inv(S), d=d, crystallinity=phi, amorphous_young=float(amorphous_young),
                amorphous_poisson=float(amorphous_poisson), notes=notes)
