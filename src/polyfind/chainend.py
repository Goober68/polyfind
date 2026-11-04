"""Chain-end and registry energetics for the pore-mouth machine.

The fluoro-ribosome design (``sarco/shared/ip/fluoro-ribosome-design-v1.md``; disclosure
sections 10.16 and 10.23) grows a PVDF-type chain one unit at a time at a pore mouth, with a
transverse electrode pair meant to write the torsion of the bond entering the tunnel
(design 2.5, Dipole item D21), a registry ring meant to notch the chain's axial position
(design 2.4, D21) and a backtracking editor that pulls the chain back one notch (design
10.2, D24).  The bond-forming chemistry at the end is electronic structure and stays with
the DFT provider (D13, D19, D20, D22, D23).  What the force field can say is the
conformational and lattice side of those steps, and this module says two things:

* :func:`torsion_gate`: at the end of an oligomer, the energy of a gauche against a trans
  entering bond, the change of the chain dipole transverse to the tunnel axis that comes
  with it, and the transverse field that biases the bond either way by a chosen margin.
  The field couples only to the transverse part of the dipole *change*, which is why the
  number is not the unit dipole times the field.
* :func:`registry_corrugation`: the lattice energy of a crystal as one chain slides along
  its axis through one repeat, per monomer: the corrugation a perfectly commensurate
  partner exerts, its half-slip value (the "registry" number of the provider's screen) and
  the peak force per monomer that a backtrack against it would have to overcome.

Both are what the fitted potential gives on frozen bond angles, with the caveats the
functions state; the plan's section 10 reads them against the design's deciding numbers.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .chain import build_chain
from .polymers import THREE_STATE

K_B_EV = 8.617333262e-5
KCAL_PER_EV = 23.060548
KCAL_MOL_A_TO_PN = 4184.0 / 6.02214076e23 / 1e-10 * 1e12  # 69.48 pN per kcal/mol/A


def chain_axis(chain, n_backbone: int | None = None) -> np.ndarray:
    """Unit vector along the chain: the principal axis of the first ``n_backbone`` backbone atoms."""
    bb = chain.coords[np.asarray(chain.backbone)[: n_backbone or len(chain.backbone)]]
    centred = bb - bb.mean(axis=0)
    _, _, vt = np.linalg.svd(centred, full_matrices=False)
    z = vt[0]
    if z @ (bb[-1] - bb[0]) < 0:
        z = -z
    return z / np.linalg.norm(z)


@dataclass
class TorsionGate:
    """One entering bond at a chain end, trans against gauche, in a transverse field.

    ``delta_e`` is ``E(gauche) - E(trans)`` in kcal/mol for the whole oligomer (positive:
    trans favoured at zero field), ``gauche`` the state angle that gave the lower gauche
    energy, ``dmu_perp`` the magnitude (e.A) of the dipole change transverse to the chain
    axis, ``field_to_gauche`` / ``field_to_trans`` the transverse field (V/A) at which the
    bond is biased the named way by ``margin_kcal``, zero when it already is; ``inf`` when
    there is no transverse dipole change to couple to.  ``field_axis`` is the direction of
    the dipole change in the chain frame (the gate's "+" polarity for gauche).
    """

    polymer: str
    bond: int
    gauche: float
    delta_e: float
    dmu_perp: float
    field_axis: np.ndarray
    margin_kcal: float
    field_to_gauche: float
    field_to_trans: float

    def table(self) -> str:
        return (f"{self.polymer:8s} bond {self.bond:3d} G{self.gauche:+.0f}: E_G - E_T = {self.delta_e:+8.3f} kcal/mol, "
                f"|dmu_perp| = {self.dmu_perp:.3f} e.A; field for a {self.margin_kcal:.2f} kcal/mol bias: "
                f"to gauche {self.field_to_gauche * 10:.2f} V/nm, to trans {self.field_to_trans * 10:.2f} V/nm")


def torsion_gate(polymer, ff, n_dihedrals: int = 16, bond: int | None = None, states=THREE_STATE,
                 margin_kt: float = 2.0, temperature: float = 298.15, axis_atoms: int | None = None) -> TorsionGate:
    """The trans/gauche gate at the end of an all-trans oligomer of ``polymer`` under ``ff``.

    ``bond`` is the dihedral index (default the last, the bond entering the tunnel);
    ``margin_kt`` the bias asked of the gate in units of ``kT``.  The chain axis is the
    principal axis of the backbone atoms before the gated bond (``axis_atoms`` overrides
    how many), so the transverse plane is the tunnel's.  Both gauche states are tried and
    the lower one kept; for a chiral unit they differ.  Bond angles are frozen, so for a
    chlorine-bearing unit the all-trans reference is strained and ``delta_e`` is not a
    physical number (``docs/CHEMISTRY_EXTENSION.md`` phase 1); the dipole change is still
    the geometry's.
    """
    import warnings

    n = int(n_dihedrals)
    k = n - 1 if bond is None else int(bond)
    margin_kcal = margin_kt * K_B_EV * temperature * KCAL_PER_EV
    trans = np.full(n, states.angles[0])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ch_t = build_chain(polymer, trans)
        e_t = float(ff.energy(ch_t))
        n_axis = axis_atoms if axis_atoms is not None else max(4, k + 1)
        z = chain_axis(ch_t, n_axis)
        mu_t = ch_t.charges @ ch_t.coords
        best = None
        for g in states.angles[1:]:
            d = trans.copy()
            d[k] = g
            ch_g = build_chain(polymer, d)
            e_g = float(ff.energy(ch_g))
            if best is None or e_g < best[0]:
                best = (e_g, float(g), ch_g.charges @ ch_g.coords)
    e_g, g, mu_g = best
    dmu = mu_g - mu_t
    dmu_perp = dmu - z * (dmu @ z)
    m = float(np.linalg.norm(dmu_perp))
    delta_e = e_g - e_t
    if m < 1e-9:
        to_g = to_t = float("inf")
        axis = np.zeros(3)
    else:
        axis = dmu_perp / m
        to_g = max(0.0, (delta_e + margin_kcal) / KCAL_PER_EV / m)
        to_t = max(0.0, (margin_kcal - delta_e) / KCAL_PER_EV / m)
    return TorsionGate(polymer=polymer.name, bond=k, gauche=g, delta_e=delta_e, dmu_perp=m, field_axis=axis,
                       margin_kcal=margin_kcal, field_to_gauche=to_g, field_to_trans=to_t)


@dataclass
class Corrugation:
    """Lattice energy per monomer against the axial slip of one chain, over one repeat."""

    dz: np.ndarray  # A, 0 .. c
    energy: np.ndarray  # kcal/mol per monomer, minimum subtracted
    c: float
    amplitude: float  # max - min
    half_slip: float  # E(c/2) - min: the provider's "registry" number
    peak_slope: float  # kcal/mol/A per monomer
    monomers_per_cell: int

    @property
    def peak_force_pn(self) -> float:
        return self.peak_slope * KCAL_MOL_A_TO_PN

    def table(self) -> str:
        return (f"corrugation over c = {self.c:.3f} A: amplitude {self.amplitude:.3f} kcal/mol per monomer, "
                f"half-slip {self.half_slip:.3f}, peak slope {self.peak_slope:.3f} kcal/mol/A = "
                f"{self.peak_force_pn:.0f} pN per monomer")


def registry_corrugation(ref, monomers_per_cell: int, n: int = 41) -> Corrugation:
    """Slide chain 2 of ``ref``'s cell along the chain axis through one repeat.

    ``ref`` is a :class:`polyfind.mechanics.Reference` (its ``params[5]`` is the offset
    ``dz``); every other cell parameter is held, so this is the corrugation at the
    crystal's own cell, not a re-relaxed slip.  Energies are per monomer
    (``monomers_per_cell`` is the packer's chain count times the monomers in its repeat).
    """
    p = np.asarray(ref.params, dtype=float).copy()
    dz = np.linspace(0.0, ref.c, int(n))
    rows = np.repeat(p[None], len(dz), axis=0)
    rows[:, 5] = p[5] + dz
    e = np.asarray(ref.packer.energy(rows, c=np.full(len(dz), ref.c)), dtype=float) / int(monomers_per_cell)
    e = e - e.min()
    half = int(np.argmin(np.abs(dz - ref.c / 2.0)))
    slope = np.gradient(e, dz)
    return Corrugation(dz=dz, energy=e, c=float(ref.c), amplitude=float(e.max()), half_slip=float(e[half]),
                       peak_slope=float(np.max(np.abs(slope))), monomers_per_cell=int(monomers_per_cell))
