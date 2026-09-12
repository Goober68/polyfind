"""Energy backends and RIS parameter fitting.

Two roles:

1. A pluggable :class:`Calculator` protocol.  :class:`SimpleFF` is a small,
   transparent intramolecular potential (Fourier torsion on the backbone
   dihedral + UFF Lennard-Jones + Coulomb for 1-4 and beyond) that exists so the
   whole pipeline runs and can be tested without heavy dependencies.  It can
   optionally carry harmonic bond and angle terms and an off-site charge position
   (``bond_terms``, ``angle_terms``, ``charge_offsets``), which is what makes the
   reference *forces* fittable and lets a strained chain relieve a contact by
   opening an angle; all three are off by default, so a default-constructed
   instance is the rigid, atom-centred potential it has always been.  See
   docs/VALENCE_FIT.md for what they bought and what they did not.  Its
   default parameters are *illustrative* -- they were chosen to be reasonable,
   not fitted -- but they are now all constructor arguments, and
   :mod:`polyfind.fitting` fits them to experimental crystal data and offers the
   result as a named preset (the defaults are untouched by that, so nothing
   changes for anyone who does not ask for one).  :class:`ASECalculator` wraps
   any ASE calculator (e.g. a MACE machine-learned potential) behind the same
   interface; per DESIGN.md 2.5 that is a validator at the end of the funnel and
   deliberately not a dependency, and fitting this potential is the productive
   direction instead.

2. :func:`fit_ris`: derive the RIS first- and second-order energies from a
   calculator by scanning one and two consecutive backbone dihedrals of a short
   oligomer.  This is the *only* place the (possibly expensive) potential is
   evaluated densely: B*(360/step) + B*(360/step)^2 single points, about 2.6k
   for step = 10 deg and B = 2.  The scan is embarrassingly parallel and goes
   through ``energy_batch`` so a GPU-backed calculator can evaluate all
   conformers at once.

3. :class:`Frame` and its helpers: an *arbitrary* all-atom geometry with its
   topology inferred from the coordinates, so that the same potential can be
   scored against reference (density-functional) data that was not produced by
   :func:`~polyfind.chain.build_chain`.  :mod:`polyfind.fitting` uses this to fit
   the potential to torsion scans and conformer energies; see that module for the
   objective and DESIGN.md 5.7 for why crystal data alone could not do the job.
"""
from __future__ import annotations

import os
import re
import warnings
from dataclasses import dataclass, field
from typing import Protocol, Sequence

import numpy as np

from . import backend as bk
from .chain import Structure, build_chain, build_chain_batch
from .polymers import Polymer, RISStates, THREE_STATE, UFF_LJ, lj_params
from .ris import RISModel

COULOMB = 332.0637  # kcal A / (mol e^2)

# A batch of conformers sharing one topology: (template, coords) with coords shape (M, n, 3),
# as produced by :func:`polyfind.chain.build_chain_batch`.
ConformerBatch = tuple[Structure, np.ndarray]


class Calculator(Protocol):
    """Energies in kcal/mol for one or many rigid-geometry oligomers.

    ``energy_batch`` accepts either a list of :class:`~polyfind.chain.Structure`
    (independent topologies) or a single :data:`ConformerBatch` ``(template, coords)``
    sharing one topology -- the form produced by :func:`~polyfind.chain.build_chain_batch`.
    The tuple form lets a GPU-resident or otherwise batched calculator (an MLIP, say)
    evaluate every row as one forward pass instead of looping; :class:`SimpleFF` does
    this via :meth:`SimpleFF.energy_coords`, and :class:`ASECalculator` loops but still
    accepts the tuple form so callers need not branch on it.
    """

    def energy(self, struct: Structure) -> float: ...

    def energy_batch(self, structs: Sequence[Structure] | ConformerBatch) -> np.ndarray: ...


def erfc_approx(x, xp=np):
    """Abramowitz-Stegun 7.1.26 (|err| < 1.5e-7), elementary ops only (runs on CuPy)."""
    ax = xp.abs(x)
    t = 1.0 / (1.0 + 0.3275911 * ax)
    poly = t * (0.254829592 + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))))
    y = poly * xp.exp(-ax * ax)
    return xp.where(x >= 0, y, 2.0 - y)


# ------------------------------------------------- arbitrary geometries (Frame)
#
# Everything in this package below :class:`~polyfind.chain.Structure` assumes rigid
# bond geometry: only the dihedrals vary, and the topology comes from the
# :class:`~polyfind.polymers.Polymer` that built the chain.  Reference data does not
# arrive that way -- a relaxed torsion scan moves bond lengths and angles too, and it
# names no polymer -- so a :class:`Frame` carries an arbitrary geometry with its
# topology *inferred from the coordinates*, and :class:`SimpleFF` scores it with the
# same three terms it uses everywhere else.

# Cordero et al., Dalton Trans. 2008, 2832 (single-bond covalent radii, A).
COVALENT_RADII: dict[str, float] = {
    "H": 0.31, "B": 0.84, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57, "Si": 1.11,
    "P": 1.07, "S": 1.05, "Cl": 1.02, "Br": 1.20, "I": 1.39,
}
BOND_TOLERANCE = 1.25  # a bond is r_ij < tolerance * (rcov_i + rcov_j)

# Environment variable naming the reference set :func:`read_frames` reads by default.
# The data is *not* vendored into this repository: it is large, it is not ours, and
# pinning it to one machine's path would be worse than asking for the path.
TRAINSET_ENV = "POLYFIND_TRAINSET"
EV_TO_KCAL = 23.0605


def infer_bonds(elements, coords, tolerance: float = BOND_TOLERANCE) -> list[tuple[int, int]]:
    """Bonds of a geometry, by the usual covalent-radius overlap rule.

    ``r_ij < tolerance * (rcov_i + rcov_j)``.  For the saturated fluoro-oligomers this
    package cares about the rule is unambiguous: the widest bond (C-C at 1.55 A) sits
    at 0.82 of its cutoff and the tightest non-bond (a geminal F...F at 2.2 A) at 1.9
    of its own, so nothing lands near the boundary.  A double bond is *not*
    distinguished from a single one -- callers that need to exclude unsaturated
    backbones test the bond length themselves (see :meth:`Frame.backbone_is_saturated`).
    """
    coords = np.asarray(coords, dtype=float)
    try:
        rad = np.array([COVALENT_RADII[e] for e in elements])
    except KeyError as e:  # pragma: no cover - guards a typo in a data file
        raise KeyError(f"no covalent radius for element {e.args[0]!r}") from e
    d = np.linalg.norm(coords[:, None, :] - coords[None, :, :], axis=-1)
    cut = tolerance * (rad[:, None] + rad[None, :])
    i, j = np.where(np.triu(d < cut, 1))
    return [(int(a), int(b)) for a, b in zip(i, j)]


def neighbour_lists(n_atoms: int, bonds) -> list[list[int]]:
    adj: list[list[int]] = [[] for _ in range(n_atoms)]
    for a, b in bonds:
        adj[a].append(b)
        adj[b].append(a)
    return adj


def carbon_backbone(elements, bonds) -> list[int]:
    """The longest simple path through the carbon atoms, in chain order.

    That is the backbone of a linear oligomer, and it is what
    :attr:`~polyfind.chain.Structure.backbone` holds for a built chain.  Exhaustive
    depth-first search rather than the usual double-breadth-first trick, because the
    latter is only correct on a tree and a pendant ring (cyclopropyl, say) is not one;
    with a dozen or two carbons the exhaustive search costs nothing.
    """
    carbons = [i for i, e in enumerate(elements) if e == "C"]
    if not carbons:
        return []
    if len(carbons) > 24:  # pragma: no cover - guard, not a supported case
        raise ValueError(f"exhaustive backbone search refuses {len(carbons)} carbons")
    among = set(carbons)
    adj = {c: [] for c in carbons}
    for a, b in bonds:
        if a in among and b in among:
            adj[a].append(b)
            adj[b].append(a)
    best: list[int] = []

    def walk(path, seen):
        nonlocal best
        if len(path) > len(best):
            best = list(path)
        for nxt in adj[path[-1]]:
            if nxt not in seen:
                seen.add(nxt)
                path.append(nxt)
                walk(path, seen)
                path.pop()
                seen.discard(nxt)

    for start in carbons:
        walk([start], {start})
    return best


def bci_charges(elements, bonds, increments) -> np.ndarray:
    """Point charges from bond-charge increments: ``q_i = sum_j delta(e_i, e_j)``.

    ``increments`` maps an *ordered* element pair to the charge the first element gains
    per bond to the second; the reverse pair gains the negative, so the molecule is
    neutral by construction and a bond between two atoms of the same element carries
    nothing.  Only pairs that appear need an entry; anything missing is zero.

    This is the charge model :mod:`polyfind.polymers` already uses, written as
    parameters instead of literals.  ``{("C", "H"): -0.10, ("C", "F"): +0.20}``
    reproduces PVDF's charges exactly (CH2 carbon -0.20, H +0.10, CF2 carbon +0.40,
    F -0.20), ``{("C", "H"): -0.06}`` reproduces polyethylene's, and
    ``{("C", "Cl"): +0.10}`` PVDC's.  Because every backbone-backbone bond is C-C and
    therefore carries no increment, the charges of a periodic block do not depend on
    the bonds that cross its boundary -- which is why the same function serves an
    oligomer and one repeat of a :class:`~polyfind.pack.PeriodicChain`.
    """
    q = np.zeros(len(elements))
    for a, b in bonds:
        ea, eb = elements[a], elements[b]
        d = increments.get((ea, eb))
        if d is None:
            d = -increments.get((eb, ea), 0.0)
        q[a] += d
        q[b] -= d
    return q


# ------------------------------------------------------------------ valence terms
#
# Harmonic bond stretching and angle bending, typed by the bond graph rather than by the
# polymer, so that one table serves a built :class:`~polyfind.chain.Structure`, a periodic
# block and an arbitrary :class:`Frame`.  They are **off by default**: a
# default-constructed :class:`SimpleFF` has no valence terms at all, which is the rigid
# potential every earlier number in this package was computed with.  What they buy is
# recorded in docs/VALENCE_FIT.md; the two things they are for are (i) making the
# reference *forces* representable, 93% of which are bond stretching, and (ii) letting a
# strained all-trans chain relieve itself by opening an angle instead of carrying hundreds
# of kcal/mol of Lennard-Jones contact.
#
# Inside a packing they change nothing, and that is by construction rather than by
# accident: :mod:`polyfind.pack` builds rigid chains, so every bond length and every bond
# angle is fixed and the valence energy is one additive constant per chain.  It therefore
# cancels from every lattice-energy *difference* and from every RIS energy, which are
# measured from all-trans.  ``tests/test_forcefield.py`` asserts exactly that.

# Elements that, when they carry a single bond in the graph, are terminal atoms of a
# multiple bond: a one-coordinate nitrogen is a nitrile and a one-coordinate oxygen a
# carbonyl or a sulfonyl.  They need their own stretch type because their length (C#N
# 1.14 A) has nothing to do with the single bond of the same element pair (C-N 1.47 A) --
# measured on the reference set, typing them together gives one "C-N" population spanning
# 1.16 to 1.52 A, which no single harmonic describes.
MULTIPLE_TERMINAL: tuple[str, ...] = ("N", "O")
WILDCARD = "*"  # the catch-all entry of a valence table


def bond_type_name(elements, adj, i: int, j: int) -> str:
    """Stretch type of the bond ``i-j``: ``"C-F"``, ``"C=N"`` (terminal multiple), ...

    The element pair in alphabetical order, joined by ``-`` for a single bond and ``=``
    for a bond to a one-coordinate :data:`MULTIPLE_TERMINAL` atom.  Purely a function of
    the graph and the elements, so it means the same thing for a built chain and for an
    inferred :class:`Frame`.
    """
    ei, ej = elements[i], elements[j]
    multiple = ((len(adj[i]) == 1 and ei in MULTIPLE_TERMINAL)
                or (len(adj[j]) == 1 and ej in MULTIPLE_TERMINAL))
    a, b = sorted((ei, ej))
    return f"{a}={b}" if multiple else f"{a}-{b}"


def angle_type_name(elements, adj, i: int, j: int, k: int) -> str | None:
    """Bend type of the angle ``i-j-k`` (``j`` central): ``"C-C-C"``, ``"F-C-H"``, ...

    The two outer elements in alphabetical order around the central one.  ``None`` for a
    **two-coordinate centre**, which is a linear group (a nitrile carbon, an isocyanate
    nitrogen): its reference angles sit at 179-180 degrees, where the derivative of
    ``arccos`` is singular and a harmonic in the angle is the wrong functional form.  Such
    groups are rigid in :mod:`polyfind.polymers` by construction too (``nitrile()`` adds
    atoms but no degrees of freedom), so leaving them without a bend term is the same
    modelling choice made in the same place twice, not a new approximation.
    """
    if len(adj[j]) < 3:
        return None
    a, b = sorted((elements[i], elements[k]))
    return f"{a}-{elements[j]}-{b}"


def valence_topology(elements, bonds) -> tuple[list, list]:
    """``(bonds, angles)`` with their type names: ``[(i, j, name)]`` and ``[(i, j, k, name)]``.

    Angles are every pair of neighbours of every atom with at least three of them; ``j``
    is the central atom.  Deterministic order (by central atom, then by neighbour index)
    so that two calls on the same topology give the same arrays.
    """
    adj = neighbour_lists(len(elements), bonds)
    bl = [(int(i), int(j), bond_type_name(elements, adj, int(i), int(j))) for i, j in bonds]
    al = []
    for j in range(len(elements)):
        nb = sorted(adj[j])
        for x in range(len(nb)):
            for y in range(x + 1, len(nb)):
                name = angle_type_name(elements, adj, nb[x], j, nb[y])
                if name is not None:
                    al.append((int(nb[x]), int(j), int(nb[y]), name))
    return bl, al


def _valence_lookup(table: dict, name: str, what: str) -> tuple[float, ...]:
    """``table[name]``, falling back to the :data:`WILDCARD` entry; raises if neither."""
    if name in table:
        return table[name]
    if WILDCARD in table:
        return table[WILDCARD]
    raise KeyError(f"no {what} parameters for type {name!r} and no {WILDCARD!r} entry; "
                   f"known: {sorted(table)}")


def offset_sites(elements, bonds, offsets) -> list[tuple[int, int, float]]:
    """``(atom, neighbour, distance)`` for every atom carrying an off-site charge.

    ``offsets`` maps an element symbol to a displacement along the bond *away from* its
    one neighbour (so positive moves the charge outward, beyond the nucleus).  Only
    one-coordinate atoms may carry one -- the construction is "slide the charge along the
    bond", and an atom with two bonds has no single bond to slide along -- and an element
    in the table that appears with more than one neighbour raises rather than being
    silently skipped.
    """
    if not offsets:
        return []
    adj = neighbour_lists(len(elements), bonds)
    out = []
    for i, e in enumerate(elements):
        d = offsets.get(e)
        if d is None:
            continue
        if len(adj[i]) != 1:
            raise ValueError(f"off-site charge asked for {e} atom {i}, which has "
                             f"{len(adj[i])} neighbours; only terminal atoms can carry one")
        out.append((int(i), int(adj[i][0]), float(d)))
    return out


def charge_site_coords(coords, sites) -> np.ndarray:
    """``coords`` with every :func:`offset_sites` atom moved along its bond (n, 3)."""
    coords = np.asarray(coords, dtype=float)
    if not sites:
        return coords
    out = coords.copy()
    for a, nb, d in sites:
        u = coords[a] - coords[nb]
        out[a] = coords[a] + d * u / np.linalg.norm(u)
    return out


def project_offset_charges(coords, charges, sites) -> np.ndarray:
    """Atom-centred charges with the same total dipole as the off-site model.

    Moving charge ``q`` a distance ``d`` along a bond of length ``r`` adds ``q d`` to the
    dipole along that bond; putting ``q d / r`` more charge on the atom and the same amount
    less on its neighbour adds exactly the same thing, with the total charge unchanged.  So
    this is the *exact* dipole-equivalent of the off-site model, and it differs from it only
    in the quadrupole and above.

    It exists because :mod:`polyfind.pack` and :mod:`polyfind.refine` carry one charge per
    atom and no off-site machinery, so a fitted off-site model reaches a crystal through
    its dipole-equivalent projection rather than exactly.  At the 3 A and longer separations
    a lattice sum is made of, the dipole is the term that matters; at the 2.5 A intramolecular
    contacts the fit is done on, it is not, which is the whole point of having the off-site
    site in the first place.  Both facts are measured in docs/VALENCE_FIT.md.
    """
    q = np.array(charges, dtype=float)
    coords = np.asarray(coords, dtype=float)
    for a, nb, d in sites:
        r = float(np.linalg.norm(coords[a] - coords[nb]))
        delta = q[a] * d / r
        q[a] += delta
        q[nb] -= delta
    return q


# ------------------------------------------------------------------ charge flux
#
# Fixed bond-charge increments on a rigid-bonded chain cannot produce a piezoelectric
# response for a planar all-trans zigzag, and the reason is a symmetry statement rather
# than a numerical limit: the cell dipole is a sum over bonds of ``delta_ij (r_i - r_j)``,
# a backbone C-C bond carries no increment, the pendant bonds have fixed lengths, and the
# zigzag's mirror pins each pendant pair's bisector perpendicular to the chain axis
# whatever the backbone angle is.  Axial strain moves only the backbone angle, so the
# dipole cannot move (docs/ELECTROMECHANICS.md 5.2, docs/BENCHMARK.md).
#
# **Charge flux** is the standard way out: let the increment itself depend on the local
# geometry, so that charge redistributes as the chain deforms.  The form here is
#
#     delta_ij = delta0(e_i, e_j) + ka(e_i, e_j) * G_ij + kb(e_i, e_j) * (r_ij - r0_ij)
#
# with ``i`` the atom that *gains* ``delta0`` (the first element of the ordered pair in
# :attr:`SimpleFF.charge_increments`, so C for C-H and C-F), and
#
#     G_ij = sum over k bonded to i, k != j, of (cos theta_kij - cos theta_tet)
#
# the **angle driver**: the sum of the angles this bond makes at its own atom, measured
# from the ideal sp3 value ``cos theta_tet = -1/3``.  Two properties earn it that form.
# It has a *natural zero* -- an ideal tetrahedral centre gives ``G = 0`` exactly, so the
# flux needs no fitted reference angle of its own and the base charges are the BCI ones at
# ideal geometry -- and it is a function of the bond graph and the coordinates alone, so
# the same expression serves a built chain, one repeat of a periodic chain and an
# arbitrary :class:`Frame`.  The **bond driver** ``r_ij - r0_ij`` is the *directional*
# channel: charge moves along the bond by the bond's own stretch, so a displacement of the
# terminal atom projected on the bond direction transfers charge and one across it does not.
# Its Born effective charge (:mod:`polyfind.born`) is ``-k_bond r u u^T`` on the terminal
# atom, ``u`` the unit bond vector -- large along the bond, zero across it, which is the
# anisotropy the periodic-DFPT tensor of beta-PVDF shows for fluorine (-1.3 e along the C-F
# bond, -0.4 e across it) and the angle driver, whose Born signature is the traceless
# ``u n^T`` with ``n`` perpendicular to the bond, cannot express.  Its zero depends on the
# path.  For a molecule (:meth:`SimpleFF.charges_at`) ``r0`` is the fitted stretch
# reference (:meth:`SimpleFF.bond_table`), because bonds there do stretch.  On the
# lattice path (:func:`polyfind.pack.chain_flux`) ``r0`` is **the bond's own length in the
# built chain**: :func:`polyfind.chain.build_chain` places every atom at the polymer's own
# bond length, the increments were fitted with the atoms there, and referencing the driver
# to a different fit's ``r0`` would shift every charge by an arbitrary constant
# (``k_bond (1.35 - 1.467)`` = -0.07 e per fluorine at ``k_bond = 0.6``).  With the built
# length as its zero the driver vanishes identically in a rigid-bond crystal -- no charge
# shift, no energy, no gradient, exactly as the stretch energy is inert
# (docs/ELECTROMECHANICS.md 4.2) -- and the channel is a pure response, seen only by a
# displacement that actually stretches a bond: the Born tensor.  That is what it is fitted
# to (``examples/fit_born_flux.py``).
#
# Neutrality is preserved for free: whatever ``delta_ij`` is, ``q_i += delta`` and
# ``q_j -= delta``, so the block's total charge is unchanged and the dipole stays a
# property of the block rather than of where the origin is.  Homonuclear pairs are
# **refused** rather than fluxed: the two atoms of a C-C bond would need an orientation to
# tell which one gains, and the only thing available to choose it with is the arbitrary
# index order of the bond list.

FLUX_COS_TET = -1.0 / 3.0  # cos of the ideal sp3 angle: the angle driver's natural zero


@dataclass
class FluxTopology:
    """Geometry-driven bond-charge increments of one block, resolved to index arrays.

    Built by :meth:`SimpleFF.flux_topology` for a molecule and by
    :func:`polyfind.pack.chain_flux` for one repeat of a periodic chain, which is why
    every atom carries an **image index**: a backbone bond joining the last atom of the
    repeat to the first atom of the next one has its partner at ``+1``, and the driver of
    such a bond depends on the repeat ``c`` as well as on the coordinates.  ``bsi == 0``
    marks the bonds whose gaining atom is in this repeat and ``bsj == 0`` those whose
    losing atom is; a bond wholly inside the repeat has both.

    ``off_*`` carry the off-site charge projection of :func:`project_offset_charges`,
    applied *after* the flux so that a fitted off-site model reaches a lattice as its
    exact dipole-equivalent, as it did before flux existed.
    """

    n_atoms: int
    bi: np.ndarray  # (nb,) the atom that gains the increment
    bj: np.ndarray  # (nb,) the other atom
    bsi: np.ndarray  # (nb,) image of bi, in units of c along z
    bsj: np.ndarray  # (nb,) image of bj
    d0: np.ndarray  # (nb,) base increment (e)
    ka: np.ndarray  # (nb,) angle-flux coefficient (e per unit of the dimensionless driver)
    kb: np.ndarray  # (nb,) bond-flux coefficient (e/A)
    r0: np.ndarray  # (nb,) reference bond length (A)
    ea: np.ndarray  # (ne,) which bond each angle entry drives
    ek: np.ndarray  # (ne,) the third atom of that angle
    eks: np.ndarray  # (ne,) its image
    off_a: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=int))
    off_n: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=int))
    off_d: np.ndarray = field(default_factory=lambda: np.empty(0))
    scale: float = 1.0  # SimpleFF.charge_scale, applied to the finished charges

    @property
    def n_bonds(self) -> int:
        return int(self.bi.size)

    @property
    def fluxes(self) -> bool:
        """Whether any coefficient is non-zero (a topology of all zeros is the BCI model)."""
        return bool(np.any(self.ka != 0.0) or np.any(self.kb != 0.0))

    # --- geometry -------------------------------------------------------------
    def _positions(self, X, cz):
        """Image-shifted positions of the three atom roles; ``X`` (M, n, 3), ``cz`` (M,)."""
        z = np.zeros(3)
        def at(idx, img):
            P = X[:, idx, :].copy()
            P[..., 2] += img[None, :] * cz[:, None]
            return P
        return at(self.bi, self.bsi), at(self.bj, self.bsj), at(self.ek, self.eks)

    def increments(self, coords, c) -> np.ndarray:
        """``delta_ij`` for every oriented bond, one row per row of ``coords`` (M, nb)."""
        X = np.asarray(coords, dtype=float)
        if X.ndim == 2:
            X = X[None]
        cz = np.broadcast_to(np.asarray(c, dtype=float).ravel(), (X.shape[0],))
        Pi, Pj, Pk = self._positions(X, cz)
        u = Pj - Pi
        r = np.linalg.norm(u, axis=-1)
        G = np.zeros((self.n_bonds, X.shape[0]))
        if self.ea.size:
            # The same expressions, in the same order, as :meth:`charges_and_grad` -- so the
            # two agree bit for bit and the kernel's value stays the number ``energy`` gives.
            uh = u / r[..., None]
            v = Pk - Pi[:, self.ea, :]
            vh = v / np.linalg.norm(v, axis=-1)[..., None]
            cosang = (uh[:, self.ea, :] * vh).sum(-1)
            np.add.at(G, self.ea, (cosang - FLUX_COS_TET).T)
        return self.d0[None, :] + self.ka[None, :] * G.T + self.kb[None, :] * (r - self.r0[None, :])

    def _project(self, q, X, dq=None, dc=None):
        """Apply the off-site dipole-equivalent projection, with its derivatives.

        ``q`` (n,), ``dq`` (n, n, 3), ``dc`` (n,) are modified in place and returned.  The
        projection moves ``q_a d / r`` from the terminal atom's neighbour onto the atom, so
        it depends on the geometry twice over -- through ``q_a``, which now fluxes, and
        through the bond length ``r`` -- and both derivatives are carried.
        """
        for a, nb, d in zip(self.off_a, self.off_n, self.off_d):
            w = X[a] - X[nb]
            r = float(np.linalg.norm(w))
            f = d / r
            if dq is not None:
                # delta = q_a f; d(delta) = f dq_a - q_a (f / r) d(r)
                dr = np.zeros_like(dq[a])
                dr[a] = w / r
                dr[nb] = -w / r
                dd = f * dq[a] - q[a] * f / r * dr
                dq[a] = dq[a] + dd
                dq[nb] = dq[nb] - dd
            if dc is not None:
                ddc = f * dc[a]  # the off-site bond never crosses the repeat boundary
                dc[a] = dc[a] + ddc
                dc[nb] = dc[nb] - ddc
            delta = q[a] * f
            q[a] += delta
            q[nb] -= delta
        return q, dq, dc

    def charges(self, coords, c=0.0) -> np.ndarray:
        """Fluxed point charges, one row per row of ``coords`` (M, n_atoms)."""
        X = np.asarray(coords, dtype=float)
        if X.ndim == 2:
            X = X[None]
        delta = self.increments(X, c)
        qT = np.zeros((self.n_atoms, X.shape[0]))
        mi, mj = self.bsi == 0, self.bsj == 0
        np.add.at(qT, self.bi[mi], delta[:, mi].T)
        np.add.at(qT, self.bj[mj], -delta[:, mj].T)
        q = np.ascontiguousarray(qT.T)
        if self.off_a.size:
            for m in range(X.shape[0]):
                self._project(q[m], X[m])
        return q * self.scale

    def charges_and_grad(self, coords, c=0.0):
        """``(q (n,), dq/dcoords (n, n, 3), dq/dc (n,))`` for ONE geometry.

        Closed form.  ``dq[a, b]`` is ``dq_a / dX_b``; ``c`` enters only through the image
        offsets, so its derivative is the z components of the same terms weighted by the
        image indices.  Verified against a central difference in ``tests/test_forcefield.py``
        and, through the lattice kernel, in ``tests/test_pack.py``.
        """
        X = np.asarray(coords, dtype=float).reshape(1, self.n_atoms, 3)
        cz = np.asarray(c, dtype=float).reshape(1)
        n, nb = self.n_atoms, self.n_bonds
        Pi, Pj, Pk = self._positions(X, cz)
        Pi, Pj, Pk = Pi[0], Pj[0], Pk[0]
        u = Pj - Pi
        r = np.linalg.norm(u, axis=-1)
        uh = u / r[:, None]
        G = np.zeros(nb)
        D = np.zeros((nb, n, 3))  # d(delta_b) / dX
        Dc = np.zeros(nb)  # d(delta_b) / dc
        rows = np.arange(nb)
        if self.kb.any():
            np.add.at(D, (rows, self.bj), self.kb[:, None] * uh)
            np.add.at(D, (rows, self.bi), -self.kb[:, None] * uh)
            Dc += self.kb * uh[:, 2] * (self.bsj - self.bsi)
        if self.ea.size:
            e = self.ea
            v = Pk - Pi[e]
            nv = np.linalg.norm(v, axis=-1)
            vh = v / nv[:, None]
            cosang = (uh[e] * vh).sum(-1)
            np.add.at(G, e, cosang - FLUX_COS_TET)
            ka = self.ka[e][:, None]
            gj = ka * (vh - cosang[:, None] * uh[e]) / r[e][:, None]
            gk = ka * (uh[e] - cosang[:, None] * vh) / nv[:, None]
            np.add.at(D, (e, self.bj[e]), gj)
            np.add.at(D, (e, self.ek), gk)
            np.add.at(D, (e, self.bi[e]), -(gj + gk))
            np.add.at(Dc, e, gj[:, 2] * self.bsj[e] + gk[:, 2] * self.eks
                      - (gj[:, 2] + gk[:, 2]) * self.bsi[e])
        delta = self.d0 + self.ka * G + self.kb * (r - self.r0)
        q = np.zeros(n)
        dq = np.zeros((n, n, 3))
        dqc = np.zeros(n)
        mi, mj = self.bsi == 0, self.bsj == 0
        np.add.at(q, self.bi[mi], delta[mi])
        np.add.at(q, self.bj[mj], -delta[mj])
        np.add.at(dq, self.bi[mi], D[mi])
        np.add.at(dq, self.bj[mj], -D[mj])
        np.add.at(dqc, self.bi[mi], Dc[mi])
        np.add.at(dqc, self.bj[mj], -Dc[mj])
        if self.off_a.size:
            self._project(q, X[0], dq, dqc)
        s = self.scale
        return q * s, dq * s, dqc * s


def flux_topology(elements, bonds, increments, flux, bond_r0=None, offsets=None,
                  images=None, scale: float = 1.0) -> FluxTopology:
    """Resolve a :class:`FluxTopology` for one block.

    ``increments`` is the bond-charge-increment mapping :func:`bci_charges` takes, ``flux``
    maps the same *ordered* element pairs to ``(k_angle, k_bond)``, and ``bond_r0`` is
    :meth:`SimpleFF.bond_table` -- ``{bond type: (k, r0)}``, of which only ``r0`` is read.
    ``images`` is ``{atom index: (position in the repeat, image)}`` for a periodic block;
    ``None`` means every atom is its own, at image 0, which is what a molecule wants.
    Atoms outside ``images`` are dropped, as are bonds with no endpoint in the repeat.
    """
    adj = neighbour_lists(len(elements), bonds)
    flux = {k: (float(a), float(b)) for k, (a, b) in (flux or {}).items()}
    for (ea, eb), (ka, kb) in flux.items():
        if ea == eb and (ka or kb):
            raise ValueError(
                f"charge flux asked for the homonuclear pair {ea}-{eb}: which of the two atoms "
                "gains the increment is then decided by the arbitrary order of the bond list, "
                "not by the chemistry, so it is refused rather than silently oriented")
    where = {i: (i, 0) for i in range(len(elements))} if images is None else dict(images)
    r0t = dict(bond_r0 or {})
    bi, bj, bsi, bsj, d0, ka, kb, r0 = [], [], [], [], [], [], [], []
    ea_, ek_, eks_ = [], [], []
    for a, b in bonds:
        if a not in where or b not in where:
            continue
        pa, pb = elements[a], elements[b]
        fwd = (pa, pb) in increments or (pa, pb) in flux
        rev = (pb, pa) in increments or (pb, pa) in flux
        if fwd and rev and pa != pb:
            raise ValueError(f"both ({pa}, {pb}) and ({pb}, {pa}) carry a charge increment or a "
                             "flux coefficient; one ordered pair per bond type")
        if fwd:
            i, j = a, b
        elif rev:
            i, j = b, a
        else:
            continue  # no increment and no flux: contributes nothing either way
        (pi, si), (pj, sj) = where[i], where[j]
        if si != 0 and sj != 0:
            continue  # the representative of this bond belongs to another repeat
        key = (elements[i], elements[j])
        fka, fkb = flux.get(key, (0.0, 0.0))
        name = bond_type_name(elements, adj, i, j)
        if fkb and name not in r0t and WILDCARD not in r0t:
            raise ValueError(f"bond charge flux on {key} needs a reference length for bond type "
                             f"{name!r}; give the potential bond_terms (or a {WILDCARD!r} entry)")
        nb_idx = len(bi)
        bi.append(pi), bj.append(pj), bsi.append(si), bsj.append(sj)
        d0.append(increments.get(key, 0.0)), ka.append(fka), kb.append(fkb)
        r0.append(r0t.get(name, r0t.get(WILDCARD, (0.0, 0.0)))[1] if fkb else 0.0)
        if fka:
            for k in adj[i]:
                if k == j or k not in where:
                    continue
                pk, sk = where[k]
                ea_.append(nb_idx), ek_.append(pk), eks_.append(sk)
    ints = lambda v: np.array(v, dtype=int)  # noqa: E731
    reals = lambda v: np.array(v, dtype=float)  # noqa: E731
    off = offset_sites(elements, bonds, offsets or {})
    off = [(where[a][0], where[nb][0], d) for a, nb, d in off
           if a in where and nb in where and where[a][1] == 0 and where[nb][1] == 0]
    return FluxTopology(
        n_atoms=len(elements) if images is None else 1 + max(p for p, _ in where.values()),
        bi=ints(bi), bj=ints(bj), bsi=reals(bsi), bsj=reals(bsj),
        d0=reals(d0), ka=reals(ka), kb=reals(kb), r0=reals(r0),
        ea=ints(ea_), ek=ints(ek_), eks=reals(eks_),
        off_a=ints([o[0] for o in off]), off_n=ints([o[1] for o in off]),
        off_d=reals([o[2] for o in off]), scale=float(scale))


def rotor_field(coords, far_mask: np.ndarray, b: int, c: int) -> np.ndarray:
    """Displacement field of a unit rotation about the bond ``b-c`` (n_atoms, 3).

    Atoms on ``c``'s side of the bond move with ``omega x (x - x_c)`` for ``omega`` the
    unit vector along ``b -> c``; the rest stand still.  This is the only
    displacement that changes a torsion angle without changing *any* bond length or
    *any* bond angle, which is exactly why it is the projection worth fitting forces
    through: a potential with no bonded terms -- this one -- cannot reproduce a
    Cartesian force component that a real molecule's bond and angle terms are holding
    up, but it can and must reproduce the torsional part.
    """
    coords = np.asarray(coords, dtype=float)
    axis = coords[c] - coords[b]
    axis = axis / np.linalg.norm(axis)
    u = np.cross(axis, coords - coords[c])
    return np.where(far_mask[:, None], u, 0.0)


def rotate_about_bond(coords, far_mask: np.ndarray, b: int, c: int, theta: float) -> np.ndarray:
    """Rotate the ``far_mask`` atoms by ``theta`` radians about the ``b -> c`` axis.

    An exact (Rodrigues) rotation, so bond lengths and bond angles are preserved to
    machine precision and a finite-difference derivative taken along it sees only the
    torsional part of the energy.
    """
    coords = np.asarray(coords, dtype=float)
    axis = coords[c] - coords[b]
    axis = axis / np.linalg.norm(axis)
    v = coords - coords[c]
    rot = (v * np.cos(theta) + np.cross(axis, v) * np.sin(theta)
           + np.outer(v @ axis, axis) * (1.0 - np.cos(theta)))
    return np.where(far_mask[:, None], coords[c] + rot, coords)


def _side_mask(n_atoms: int, adj, b: int, c: int) -> np.ndarray:
    """Atoms reachable from ``c`` without using the bond ``b-c`` (``c`` included)."""
    seen = {c}
    stack = [c]
    while stack:
        u = stack.pop()
        for v in adj[u]:
            if v in seen or (u == c and v == b):
                continue
            seen.add(v)
            stack.append(v)
    mask = np.zeros(n_atoms, dtype=bool)
    mask[list(seen)] = True
    return mask


def dihedral_angles(coords, torsions) -> np.ndarray:
    """Dihedral angles (deg) of an ``(n_tors, 4)`` index array; same convention as
    :meth:`SimpleFF._energy_from_coords` (IUPAC, trans = 180)."""
    coords = np.asarray(coords, dtype=float)
    t = np.asarray(torsions, dtype=int).reshape(-1, 4)
    a, b, c, d = (coords[t[:, k]] for k in range(4))
    b0, b1, b2 = b - a, c - b, d - c
    n1, n2 = np.cross(b0, b1), np.cross(b1, b2)
    x = (n1 * n2).sum(-1)
    y = np.linalg.norm(b1, axis=-1) * (b0 * n2).sum(-1)
    return np.degrees(np.arctan2(y, x))


@dataclass
class Frame:
    """One reference geometry: elements, coordinates, and whatever was computed for it.

    ``energy`` and ``forces`` are in kcal/mol and kcal/(mol A) -- converted on read, so
    nothing downstream has to remember what the file was in.  ``bonds``, ``backbone``
    and ``torsions`` are inferred from the coordinates by :func:`infer_bonds` and
    :func:`carbon_backbone`; ``torsions`` are the consecutive backbone quadruples, the
    same set :meth:`~polyfind.chain.Structure.dihedral_atoms` enumerates for a built
    chain.

    ``system`` names the chemistry and ``source`` how the frame was generated (a
    relaxed torsion ``scan`` or a ``conf`` ormer).  Both matter for the fit: energies
    are only comparable *within* a system, and a split that mixes frames of one system
    across train and test measures nothing (DESIGN.md 5.7).
    """

    system: str
    source: str
    elements: list[str]
    coords: np.ndarray
    energy: float | None = None
    forces: np.ndarray | None = None
    bonds: list[tuple[int, int]] = field(default_factory=list)
    backbone: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=int))
    torsions: np.ndarray = field(default_factory=lambda: np.empty((0, 4), dtype=int))
    _cache: dict = field(default_factory=dict, repr=False)

    @classmethod
    def from_geometry(cls, elements, coords, system: str = "", source: str = "",
                      energy: float | None = None, forces=None,
                      tolerance: float = BOND_TOLERANCE) -> "Frame":
        elements = list(elements)
        coords = np.asarray(coords, dtype=float)
        bonds = infer_bonds(elements, coords, tolerance)
        bb = carbon_backbone(elements, bonds)
        tors = np.array([bb[k:k + 4] for k in range(max(0, len(bb) - 3))], dtype=int).reshape(-1, 4)
        return cls(system=system, source=source, elements=elements, coords=coords,
                   energy=energy, forces=None if forces is None else np.asarray(forces, dtype=float),
                   bonds=bonds, backbone=np.array(bb, dtype=int), torsions=tors)

    @property
    def n_atoms(self) -> int:
        return len(self.elements)

    @property
    def n_dihedrals(self) -> int:
        return len(self.torsions)

    @property
    def adjacency(self) -> list[list[int]]:
        if "adj" not in self._cache:
            self._cache["adj"] = neighbour_lists(self.n_atoms, self.bonds)
        return self._cache["adj"]

    def with_coords(self, coords) -> "Frame":
        """The same molecule at different coordinates: topology kept, geometry replaced.

        The bond graph is *not* re-inferred, which is the point -- a displaced copy used
        for a finite-difference derivative must keep the topology of the frame it came
        from, or the derivative would step across a change of exclusions.  The cached
        graph distances and adjacency are shared; the geometry-dependent caches are not.
        """
        out = Frame(system=self.system, source=self.source, elements=self.elements,
                    coords=np.asarray(coords, dtype=float), energy=None, forces=None,
                    bonds=self.bonds, backbone=self.backbone, torsions=self.torsions)
        out._cache["adj"] = self.adjacency
        out._cache["dist"] = self.graph_distances()
        return out

    def graph_distances(self) -> np.ndarray:
        """Bond-count distances capped at 99, for the nonbonded exclusions."""
        if "dist" in self._cache:
            return self._cache["dist"]
        n, adj = self.n_atoms, self.adjacency
        dist = np.full((n, n), 99, dtype=int)
        for s in range(n):
            dist[s, s] = 0
            frontier = [s]
            for d in (1, 2, 3):
                nxt = []
                for u in frontier:
                    for v in adj[u]:
                        if dist[s, v] == 99:
                            dist[s, v] = d
                            nxt.append(v)
                frontier = nxt
        self._cache["dist"] = dist
        return dist

    def angles(self) -> np.ndarray:
        """The backbone dihedral angles (deg)."""
        return dihedral_angles(self.coords, self.torsions)

    def rotors(self) -> np.ndarray:
        """``(n_dihedrals, n_atoms, 3)`` rigid-rotation fields, one per backbone torsion.

        Normalised so that ``d phi_k / d theta_k = 1`` exactly (asserted in the tests):
        rotating about torsion ``k``'s central bond changes torsion ``k``'s angle by the
        rotation angle and leaves every other torsion, bond length and bond angle alone.
        """
        if "rotors" in self._cache:
            return self._cache["rotors"]
        adj = self.adjacency
        out = np.zeros((self.n_dihedrals, self.n_atoms, 3))
        for k, (a, b, c, d) in enumerate(np.asarray(self.torsions, dtype=int)):
            mask = _side_mask(self.n_atoms, adj, int(b), int(c))
            if mask[a]:  # the bond is in a ring: a rotation about it is not rigid
                raise ValueError(f"torsion {k} of system {self.system!r} sits in a ring")
            out[k] = rotor_field(self.coords, mask, int(b), int(c))
        self._cache["rotors"] = out
        return out

    def reference_torques(self) -> np.ndarray:
        """``-dE/dphi_k`` of the reference data, kcal/(mol rad), one per backbone torsion.

        The projection of the reference Cartesian forces onto :meth:`rotors`.  Because
        the rotation preserves every bond length and every bond angle, the bond and
        angle terms of whatever produced the forces contribute exactly zero here, which
        is what makes this comparable with a potential that has no such terms.
        """
        if self.forces is None:
            raise ValueError(f"frame {self.system!r} carries no forces")
        return np.einsum("kia,ia->k", self.rotors(), self.forces)

    def backbone_is_saturated(self, min_cc: float = 1.42) -> bool:
        """Every backbone C-C bond longer than ``min_cc`` (i.e. none of them multiple)."""
        c = self.coords
        return all(float(np.linalg.norm(c[i] - c[j])) >= min_cc
                   for i, j in zip(self.backbone, self.backbone[1:]))

    def backbone_is_acyclic(self) -> bool:
        """No backbone bond lies in a ring, so every backbone torsion is a real rotor."""
        adj = self.adjacency
        for i, j in zip(self.backbone, self.backbone[1:]):
            if _side_mask(self.n_atoms, adj, int(i), int(j))[int(i)]:
                return False
        return True

    def to_structure(self, polymer: Polymer) -> Structure:
        """This geometry as a :class:`~polyfind.chain.Structure` of ``polymer``.

        Only for frames whose topology really is ``polymer``'s: the atoms are reordered
        into :func:`~polyfind.chain.build_chain`'s order (each backbone atom followed by
        its pendants, then the two end caps) and the result is checked against a chain
        built from the polymer, element by element and bond by bond.  Raises if they
        disagree.  This is the bridge that proves the :class:`Frame` path and the
        :class:`Structure` path are scoring the same molecule.
        """
        order = _structure_order(self)
        elements = [self.elements[i] for i in order]
        ref = build_chain(polymer, [180.0] * (len(self.backbone) - 3))
        if elements != list(ref.elements):
            raise ValueError(f"{self.system!r} does not match {polymer.name}: elements "
                             f"{''.join(elements)} vs {''.join(ref.elements)}")
        pos = {old: new for new, old in enumerate(order)}
        bonds = sorted(tuple(sorted((pos[a], pos[b]))) for a, b in self.bonds)
        if bonds != sorted(tuple(sorted(b)) for b in ref.bonds):
            raise ValueError(f"{self.system!r} does not match {polymer.name}: bond graphs differ")
        return Structure(polymer=polymer, elements=elements, coords=self.coords[order],
                         charges=ref.charges.copy(), bonds=[tuple(b) for b in bonds],
                         backbone=np.array([pos[int(i)] for i in self.backbone], dtype=int),
                         n_dihedrals=len(self.torsions), dihedrals=self.angles(),
                         subs_of={pos[int(i)]: [pos[j] for j in self.adjacency[int(i)]
                                                if j not in set(int(x) for x in self.backbone)]
                                  for i in self.backbone})


def _structure_order(frame: Frame) -> list[int]:
    """Atom permutation putting ``frame`` into :func:`~polyfind.chain.build_chain` order.

    build_chain emits backbone atom ``k`` then its pendant atoms, for ``k`` along the
    chain, and finally the two cap hydrogens (first end, then last).  The pendants of
    one atom are ordered by element then by index so that the permutation is
    deterministic; for the CH2/CF2 atoms of a VDF oligomer the two pendants are
    identical elements, so any order gives the same element list and the same bond
    graph.  The backbone is taken in whichever direction puts the CH3 cap first, which
    is the direction :func:`~polyfind.chain.build_chain` builds.
    """
    bb = [int(i) for i in frame.backbone]
    adj = frame.adjacency
    bbset = set(bb)

    def pendants(k):
        return sorted((j for j in adj[k] if j not in bbset), key=lambda j: (frame.elements[j], j))

    # build_chain caps both ends with a single H, so the terminal backbone atoms have
    # one more pendant than the interior ones; the cap is the H that comes last.
    if len(pendants(bb[0])) < len(pendants(bb[-1])):
        bb = bb[::-1]
    order: list[int] = []
    caps: list[int] = []
    for pos, k in enumerate(bb):
        p = pendants(k)
        if pos in (0, len(bb) - 1):
            hydrogens = [j for j in p if frame.elements[j] == "H"]
            if not hydrogens:
                raise ValueError("terminal backbone atom carries no cap hydrogen")
            caps.append(hydrogens[-1])
            p = [j for j in p if j != hydrogens[-1]]
        order.append(k)
        order.extend(p)
    return order + caps


_COMMENT_KEY = re.compile(r'(\w+)=("[^"]*"|\S+)')


def read_frames(path: str | None = None, to_kcal: float = EV_TO_KCAL) -> list[Frame]:
    """Read an extended-XYZ trajectory into :class:`Frame` objects.

    ``path`` defaults to the ``POLYFIND_TRAINSET`` environment variable, so the
    reference set stays where it lives and no copy of it enters this repository.  Each
    comment line is parsed for ``energy=``, ``system=`` and ``source=``; each atom line
    is ``element x y z fx fy fz``.  ``to_kcal`` converts both the energy and the forces
    (default: from eV and eV/A, which is what the bundled set is in).
    """
    if path is None:
        path = os.environ.get(TRAINSET_ENV)
        if not path:
            raise FileNotFoundError(
                f"no reference set given and ${TRAINSET_ENV} is not set; point it at an "
                "extended-XYZ file with energy=, system= and source= on each comment line")
    with open(path) as fh:
        lines = fh.read().splitlines()
    frames: list[Frame] = []
    i = 0
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        n = int(lines[i].strip())
        info = {k: v.strip('"') for k, v in _COMMENT_KEY.findall(lines[i + 1])}
        rows = [ln.split() for ln in lines[i + 2:i + 2 + n]]
        elements = [r[0] for r in rows]
        arr = np.array([[float(x) for x in r[1:7]] for r in rows])
        energy = float(info["energy"]) * to_kcal if "energy" in info else None
        frames.append(Frame.from_geometry(
            elements, arr[:, :3], system=info.get("system", ""), source=info.get("source", ""),
            energy=energy, forces=arr[:, 3:6] * to_kcal if arr.shape[1] >= 6 else None))
        i += 2 + n
    return frames


@dataclass
class _Topology:
    pairs_i: np.ndarray
    pairs_j: np.ndarray
    lj_x: np.ndarray  # per pair
    lj_d: np.ndarray
    qq: np.ndarray  # per pair q_i q_j * scale
    torsions: np.ndarray  # (n_tors, 4) backbone quadruples
    tors_V: np.ndarray  # (n_tors, 3) Fourier coefficients of each of those dihedrals
    # Valence and off-site-charge terms.  Empty arrays mean "this potential has none",
    # which is the default and the only state every number recorded before them was
    # computed in; the kernel skips the blocks entirely rather than adding zeros.
    bond_idx: np.ndarray = field(default_factory=lambda: np.empty((0, 2), dtype=int))
    bond_k: np.ndarray = field(default_factory=lambda: np.empty(0))
    bond_r0: np.ndarray = field(default_factory=lambda: np.empty(0))
    angle_idx: np.ndarray = field(default_factory=lambda: np.empty((0, 3), dtype=int))
    angle_k: np.ndarray = field(default_factory=lambda: np.empty(0))
    angle_t0: np.ndarray = field(default_factory=lambda: np.empty(0))  # radians
    off_atom: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=int))
    off_nb: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=int))
    off_d: np.ndarray = field(default_factory=lambda: np.empty(0))


# The illustrative defaults, named so that "what the potential was before anyone fitted
# it" is a value and not a literal repeated in three modules.  Nothing here changes when
# a fitted preset is selected: a preset is a *different* SimpleFF, built on request.
# UFF Lennard-Jones entries for elements no ``Polymer`` in this package uses, so that a
# :class:`Frame` of an arbitrary chemistry can still be scored (Rappe et al., JACS 1992;
# same table and same convention as :data:`polyfind.polymers.UFF_LJ`, which wins wherever
# the two overlap -- today they do not overlap at all).
UFF_LJ_EXTRA: dict[str, tuple[float, float]] = {
    "S": (4.035, 0.274),   # S_3+2
    "Br": (3.732, 0.251),  # Br
    "I": (4.009, 0.339),   # I
    "P": (4.147, 0.305),   # P_3+3
    "Si": (4.295, 0.402),  # Si3
    "B": (4.083, 0.180),   # B_3
}

DEFAULT_TORSION: tuple[float, float, float] = (1.3, -0.05, 2.5)
DEFAULT_SCALE14 = 0.5
DEFAULT_EPS_R = 1.0

# Named parameter sets.  Each maps to :class:`SimpleFF` keyword arguments, so
# ``SimpleFF(**PRESETS[name])`` and :meth:`SimpleFF.from_preset` are the same thing.
# ``"illustrative"`` is the unfitted potential every default in this package still uses;
# the fitted entries are added by :mod:`polyfind.fitting`, which is where the data,
# the objective and the provenance of the numbers live.
PRESETS: dict[str, dict] = {
    "illustrative": {"torsion": DEFAULT_TORSION, "scale14": DEFAULT_SCALE14, "eps_r": DEFAULT_EPS_R},
}


@dataclass
class SimpleFF:
    """Small, transparent intramolecular potential (kcal/mol), with settable parameters.

    torsion = (V1, V2, V3) Fourier terms on each backbone C-C-C-C dihedral:
        V(phi) = 1/2 [V1 (1 + cos phi) + V2 (1 - cos 2 phi) + V3 (1 + cos 3 phi)]
    Nonbonded (LJ 12-6 with UFF radii/well depths, geometric combining; Coulomb with
    relative permittivity eps_r) for all pairs separated by >= 3 bonds, with 1-4 pairs
    scaled by ``scale14``.

    Every parameter is a constructor argument and every default is the illustrative
    value this potential has always had, so a default-constructed ``SimpleFF()`` is
    bit-for-bit the potential of every earlier run (``tests/test_forcefield.py``
    asserts exactly that against recorded energies):

    ``torsion``
        the shared Fourier triple.
    ``torsion_by_bond``
        one triple *per bond type* of the repeat, indexed by ``dihedral % B`` with
        ``B = polymer.bonds_per_repeat`` -- the same bond-type convention
        :func:`fit_ris` scans with.  ``None`` (the default) means "``torsion`` for
        every bond", which is what the potential did before this existed.
    ``eps_r``
        relative permittivity dividing every Coulomb term.
    ``charge_scale``
        multiplies the polymer's point charges, so Coulomb energies scale as
        ``charge_scale**2`` while dipoles and polarizations scale linearly.  It is
        exactly degenerate with ``eps_r`` in *energies* (only ``charge_scale**2 /
        eps_r`` enters) and only polarization data can tell the two apart; see
        :mod:`polyfind.fitting`.
    ``lj``
        per-element ``{symbol: (x_i, D_i)}`` overrides of the UFF table, merged over
        :data:`polyfind.polymers.UFF_LJ` for this instance only.  The module-level
        table is never touched, so two calculators with different parameters can be
        alive at once and nothing global changes underneath anyone.
    ``charge_increments``
        ``((element, element, delta), ...)`` bond-charge increments (see
        :func:`bci_charges`).  ``None`` (the default) means "use whatever charges the
        structure carries", which is the polymer's own table and what this potential
        has always done.  When set, charges are derived from the bond graph instead,
        so the *same* numbers describe a built chain, a periodic block and an
        arbitrary :class:`Frame`; ``charge_scale`` still multiplies the result.
    ``bond_terms``
        ``((type, k, r0), ...)`` harmonic bond stretching, ``0.5 k (r - r0)^2`` with
        ``k`` in kcal/(mol A^2) and ``r0`` in A, typed by :func:`bond_type_name` with
        :data:`WILDCARD` as the catch-all.  ``None`` (the default) means no stretch
        term, which is what a rigid-geometry potential has and what every number
        recorded in this package before them was computed with.
    ``angle_terms``
        ``((type, k, theta0), ...)`` harmonic angle bending, ``0.5 k (theta - theta0)^2``
        with ``k`` in kcal/(mol rad^2) and ``theta0`` in *degrees*, typed by
        :func:`angle_type_name`.  ``None`` (the default) means no bend term.
    ``charge_offsets``
        ``((element, distance), ...)`` off-site charge positions: the named element's
        point charge sits ``distance`` angstrom along its own bond, away from its
        neighbour, while its Lennard-Jones centre stays on the nucleus.  That
        separation is the point -- one atom-centred site has to serve both, and the
        C-F dipole is what suffers for it (docs/VALENCE_FIT.md).  ``None`` (the default)
        keeps every charge on its nucleus.  Only terminal atoms may carry one
        (:func:`offset_sites`).
    ``charge_flux``
        ``((element, element, k_angle, k_bond), ...)`` geometry dependence of the bond
        charge increments, on the same *ordered* pairs as ``charge_increments`` (see
        :class:`FluxTopology` for the form and for why it exists).  ``None`` (the
        default) is the fixed-increment model every earlier number was computed with.
        It is read by the **lattice** side -- ``CrystalPacker(charge_flux=ff)`` -- and by
        :meth:`charges_at`; the oligomer energy path resolves one fixed charge per
        topology and *refuses* a fluxing potential rather than silently using the base
        increments (:meth:`_topology`).
    """

    torsion: tuple[float, float, float] = DEFAULT_TORSION
    scale14: float = DEFAULT_SCALE14
    eps_r: float = DEFAULT_EPS_R
    torsion_by_bond: tuple[tuple[float, float, float], ...] | None = None
    charge_scale: float = 1.0
    lj: dict[str, tuple[float, float]] | None = None
    charge_increments: tuple[tuple[str, str, float], ...] | None = None
    bond_terms: tuple[tuple[str, float, float], ...] | None = None
    angle_terms: tuple[tuple[str, float, float], ...] | None = None
    charge_offsets: tuple[tuple[str, float], ...] | None = None
    charge_flux: tuple[tuple[str, str, float, float], ...] | None = None
    _cache: dict = field(default_factory=dict, repr=False)

    @classmethod
    def from_preset(cls, name: str) -> "SimpleFF":
        """The named parameter set from :data:`PRESETS` (e.g. ``"illustrative"``)."""
        try:
            kw = PRESETS[name]
        except KeyError as e:
            raise KeyError(f"unknown SimpleFF preset {name!r}; known: {sorted(PRESETS)}") from e
        return cls(**{k: v for k, v in kw.items() if not k.startswith("_")})

    # --- parameters -----------------------------------------------------------
    def lj_table(self) -> dict[str, tuple[float, float]]:
        """The effective ``{element: (x_i, D_i)}`` table: UFF with this instance's overrides.

        :data:`polyfind.polymers.UFF_LJ` carries only the elements the polymer
        definitions use; :data:`UFF_LJ_EXTRA` fills in the rest of the UFF main group so
        that a :class:`Frame` of some chemistry no ``Polymer`` describes can still be
        scored.  It is a fallback, never an override: an element in both takes
        ``UFF_LJ``'s value, and a caller's ``lj`` beats both.
        """
        return {**UFF_LJ_EXTRA, **UFF_LJ, **(self.lj or {})}

    def _lj_arrays(self, elements) -> tuple[np.ndarray, np.ndarray]:
        """Per-atom ``(x_i, D_i)`` with the overrides applied; ``UFF_LJ`` is left alone."""
        table = self.lj_table()
        try:
            pairs = [table[e] for e in elements]
        except KeyError as e:  # pragma: no cover - guards an unparameterised element
            raise KeyError(f"no Lennard-Jones parameters for element {e.args[0]!r}; "
                           "pass them in SimpleFF(lj=...)") from e
        arr = np.array(pairs, dtype=float)
        return arr[:, 0], arr[:, 1]

    def torsion_coefficients(self, polymer: Polymer, n_dihedrals: int) -> np.ndarray:
        """``(n_dihedrals, 3)`` Fourier coefficients, one row per backbone dihedral.

        Dihedral ``j`` belongs to bond type ``j % polymer.bonds_per_repeat`` -- the
        convention :func:`fit_ris` uses when it picks which bond to scan for each type.
        """
        if self.torsion_by_bond is None:
            return np.tile(np.asarray(self.torsion, dtype=float), (n_dihedrals, 1))
        V = np.asarray(self.torsion_by_bond, dtype=float)
        B = polymer.bonds_per_repeat
        if V.shape != (B, 3):
            raise ValueError(f"torsion_by_bond must be {B} triples for {polymer.name}, got shape {V.shape}")
        return V[[j % B for j in range(n_dihedrals)]]

    def increments(self) -> dict[tuple[str, str], float] | None:
        """:attr:`charge_increments` as the mapping :func:`bci_charges` wants."""
        if self.charge_increments is None:
            return None
        return {(a, b): float(d) for a, b, d in self.charge_increments}

    def charges_for(self, elements, bonds, fallback=None) -> np.ndarray:
        """Point charges for a topology: bond-charge increments if set, else ``fallback``.

        ``charge_scale`` is applied here, so this is the one place charges are decided
        and every caller -- :meth:`_topology`, :meth:`energy_frame`, and
        :meth:`polyfind.fitting.FFParameters.applied`'s packer -- gets the same answer.
        """
        inc = self.increments()
        q = np.asarray(fallback, dtype=float) if inc is None else bci_charges(elements, bonds, inc)
        return q * self.charge_scale if self.charge_scale != 1.0 else np.asarray(q, dtype=float)

    def bond_table(self) -> dict[str, tuple[float, float]]:
        """:attr:`bond_terms` as ``{type: (k, r0)}``; empty when there are none."""
        return {} if self.bond_terms is None else {t: (float(k), float(r0)) for t, k, r0 in self.bond_terms}

    def angle_table(self) -> dict[str, tuple[float, float]]:
        """:attr:`angle_terms` as ``{type: (k, theta0_deg)}``; empty when there are none."""
        return {} if self.angle_terms is None else {t: (float(k), float(t0)) for t, k, t0 in self.angle_terms}

    def offset_table(self) -> dict[str, float]:
        """:attr:`charge_offsets` as ``{element: distance}``; empty when there are none."""
        return {} if self.charge_offsets is None else {e: float(d) for e, d in self.charge_offsets}

    def has_valence(self) -> bool:
        return bool(self.bond_terms) or bool(self.angle_terms)

    # --- charge flux ----------------------------------------------------------
    def flux_table(self) -> dict[tuple[str, str], tuple[float, float]]:
        """:attr:`charge_flux` as ``{(element, element): (k_angle, k_bond)}``; empty when none."""
        if not self.charge_flux:
            return {}
        return {(a, b): (float(ka), float(kb)) for a, b, ka, kb in self.charge_flux}

    def has_flux(self) -> bool:
        """Whether any charge-flux coefficient is non-zero."""
        return any(ka or kb for ka, kb in self.flux_table().values())

    def flux_topology(self, elements, bonds, images=None, bond_r0=None) -> FluxTopology:
        """The :class:`FluxTopology` this potential resolves for one block.

        ``images`` is the periodic map :func:`polyfind.pack.chain_flux` supplies; ``None``
        treats the block as a molecule.  Needs :attr:`charge_increments`, because a flux is
        a perturbation of an increment and there is nothing to perturb without one.
        ``bond_r0`` (``{bond type: (k, r0)}``) is the stretch driver's zero; ``None`` takes
        the fitted stretch reference, :meth:`bond_table`, which is a molecule's choice --
        the lattice path passes the built chain's own lengths (see the note above
        :class:`FluxTopology`).
        """
        if self.charge_increments is None:
            raise ValueError("charge flux needs charge_increments: the flux is a geometry "
                             "dependence *of* the bond-charge increments, so there is nothing "
                             "for it to modify without them")
        return flux_topology(elements, bonds, self.increments(), self.flux_table(),
                             bond_r0=self.bond_table() if bond_r0 is None else bond_r0,
                             offsets=self.offset_table(), images=images, scale=self.charge_scale)

    def charges_at(self, elements, bonds, coords) -> np.ndarray:
        """Point charges of one geometry, flux included, ``charge_scale`` applied.

        The geometry-dependent counterpart of :meth:`charges_for`, and the *only* way to
        get a fluxed charge out of this class: everything that needs one has coordinates,
        and everything that does not is by definition using the fixed model.
        """
        top = self.flux_topology(elements, bonds)
        return top.charges(np.asarray(coords, dtype=float))[0]

    def _refuse_flux(self, what: str) -> None:
        if self.has_flux():
            raise ValueError(
                f"{what} resolves one fixed charge per topology and this potential has charge "
                "flux, whose charges depend on the coordinates.  Silently using the base "
                "increments would make the energy disagree with the dipole, so it is refused: "
                "use SimpleFF.charges_at(...) for the charges, or the lattice kernel "
                "(CrystalPacker(charge_flux=ff)), which carries the geometry dependence and "
                "its derivatives.")

    def valence_arrays(self, elements, bonds) -> tuple:
        """``(bond_idx, bond_k, bond_r0, angle_idx, angle_k, angle_t0_rad)`` for a topology.

        Empty arrays when this potential has no valence terms, which is what keeps the
        rigid path bit-for-bit what it was: the kernel then skips the blocks rather than
        adding zeros.
        """
        bt, at = self.bond_table(), self.angle_table()
        if not bt and not at:
            return (np.empty((0, 2), dtype=int), np.empty(0), np.empty(0),
                    np.empty((0, 3), dtype=int), np.empty(0), np.empty(0))
        bl, al = valence_topology(elements, bonds)
        if bt:
            b_idx = np.array([[i, j] for i, j, _ in bl], dtype=int).reshape(-1, 2)
            kr = np.array([_valence_lookup(bt, n, "bond") for _, _, n in bl], dtype=float).reshape(-1, 2)
        else:
            b_idx, kr = np.empty((0, 2), dtype=int), np.empty((0, 2))
        if at:
            a_idx = np.array([[i, j, k] for i, j, k, _ in al], dtype=int).reshape(-1, 3)
            kt = np.array([_valence_lookup(at, n, "angle") for _, _, _, n in al], dtype=float).reshape(-1, 2)
        else:
            a_idx, kt = np.empty((0, 3), dtype=int), np.empty((0, 2))
        return b_idx, kr[:, 0], kr[:, 1], a_idx, kt[:, 0], np.deg2rad(kt[:, 1])

    def offset_arrays(self, elements, bonds) -> tuple:
        """``(atoms, neighbours, distances)`` of the off-site charges; empty when none."""
        sites = offset_sites(elements, bonds, self.offset_table())
        if not sites:
            return np.empty(0, dtype=int), np.empty(0, dtype=int), np.empty(0)
        arr = np.array(sites, dtype=float)
        return arr[:, 0].astype(int), arr[:, 1].astype(int), arr[:, 2]

    def _param_key(self) -> tuple:
        """Everything a cached :class:`_Topology` depends on besides the structure."""
        lj = tuple(sorted((k, float(v[0]), float(v[1])) for k, v in (self.lj or {}).items()))
        tbb = None if self.torsion_by_bond is None else tuple(tuple(float(x) for x in t) for t in self.torsion_by_bond)
        inc = None if self.charge_increments is None else tuple(
            sorted((a, b, float(d)) for a, b, d in self.charge_increments))
        val = (tuple(sorted(self.bond_table().items())), tuple(sorted(self.angle_table().items())),
               tuple(sorted(self.offset_table().items())), tuple(sorted(self.flux_table().items())))
        return (tuple(float(x) for x in self.torsion), tbb, float(self.scale14),
                float(self.eps_r), float(self.charge_scale), lj, inc, val)

    # --- topology -------------------------------------------------------------
    def _pair_terms(self, elements, dist: np.ndarray, charges) -> tuple:
        """The nonbonded pair arrays shared by the :class:`Structure` and :class:`Frame`
        paths: indices, LJ minimum distance and well depth, and the Coulomb prefactor,
        for every pair at least three bonds apart with 1-4 pairs scaled."""
        n = len(elements)
        iu, ju = np.triu_indices(n, 1)
        keep = dist[iu, ju] >= 3
        iu, ju = iu[keep], ju[keep]
        scale = np.where(dist[iu, ju] == 3, self.scale14, 1.0)
        x, d = self._lj_arrays(elements)
        lj_x = np.sqrt(x[iu] * x[ju])
        lj_d = np.sqrt(d[iu] * d[ju]) * scale
        q = np.asarray(charges, dtype=float)
        qq = q[iu] * q[ju] * scale * COULOMB / self.eps_r
        return iu, ju, lj_x, lj_d, qq

    def _topology(self, struct: Structure) -> _Topology:
        self._refuse_flux("the oligomer energy path")
        key = (struct.polymer.name, struct.n_dihedrals, len(struct.elements), self._param_key())
        top = self._cache.get(key)
        if top is not None:
            return top
        n = struct.n_atoms
        adj = [[] for _ in range(n)]
        for a, b in struct.bonds:
            adj[a].append(b)
            adj[b].append(a)
        # graph distances up to 3 via BFS from each atom
        dist = np.full((n, n), 99, dtype=int)
        for s in range(n):
            dist[s, s] = 0
            frontier = [s]
            for d in (1, 2, 3):
                nxt = []
                for u in frontier:
                    for v in adj[u]:
                        if dist[s, v] == 99:
                            dist[s, v] = d
                            nxt.append(v)
                frontier = nxt
        q = self.charges_for(struct.elements, struct.bonds, struct.charges)
        iu, ju, lj_x, lj_d, qq = self._pair_terms(struct.elements, dist, q)
        tors = np.array([struct.dihedral_atoms(j) for j in range(struct.n_dihedrals)], dtype=int).reshape(-1, 4)
        tors_V = self.torsion_coefficients(struct.polymer, struct.n_dihedrals).reshape(-1, 3)
        top = _Topology(iu, ju, lj_x, lj_d, qq, tors, tors_V,
                        *self.valence_arrays(struct.elements, struct.bonds),
                        *self.offset_arrays(struct.elements, struct.bonds))
        self._cache[key] = top
        return top

    def frame_topology(self, frame: Frame, torsion=None) -> _Topology:
        """The resolved potential for one :class:`Frame`: pairs, exclusions and torsions.

        ``torsion`` is an ``(n_dihedrals, 3)`` array of Fourier coefficients, one row per
        backbone torsion; ``None`` broadcasts :attr:`torsion`.  Per-*type* coefficients
        are the caller's business -- :mod:`polyfind.fitting` types the torsions by the
        substituents on the two central carbons and passes the rows in -- because a
        frame names no polymer and so has no ``bonds_per_repeat`` to index by.
        """
        self._refuse_flux("the Frame energy path")
        if self.charge_increments is None:
            raise ValueError(
                "scoring a Frame needs charge_increments: a frame carries no polymer charge "
                "table to fall back on, so the charges have to come from the parameters "
                "(see bci_charges; importing polyfind.fitting registers a preset that has "
                "them, SimpleFF.from_preset('pvdf-dft-fit'))")
        q = self.charges_for(frame.elements, frame.bonds)
        iu, ju, lj_x, lj_d, qq = self._pair_terms(frame.elements, frame.graph_distances(), q)
        V = (np.tile(np.asarray(self.torsion, dtype=float), (frame.n_dihedrals, 1))
             if torsion is None else np.asarray(torsion, dtype=float).reshape(-1, 3))
        return _Topology(iu, ju, lj_x, lj_d, qq, np.asarray(frame.torsions, dtype=int).reshape(-1, 4), V,
                         *self.valence_arrays(frame.elements, frame.bonds),
                         *self.offset_arrays(frame.elements, frame.bonds))

    def energy_frame(self, frame: Frame, torsion=None) -> float:
        """Energy (kcal/mol) of an arbitrary geometry, topology inferred from coordinates.

        The same three terms as :meth:`energy`, evaluated through the same kernel; the
        only difference is where the topology came from.  Absolute values are not
        comparable with a reference method's -- this potential has no bonded terms and
        no zero of energy -- so only *differences within one chemistry* mean anything,
        which is exactly how :mod:`polyfind.fitting` uses them.
        """
        top = self.frame_topology(frame, torsion)
        return float(self._energy_from_coords(top, frame.coords[None], np)[0])

    def frame_torques(self, frame: Frame, torsion=None, step: float = 1e-5) -> np.ndarray:
        """``-dE/dphi_k`` (kcal/(mol rad)) about each backbone torsion of ``frame``.

        Central differences over a *rigid* rotation of one side of each central bond
        (:meth:`Frame.rotors`), which is the same displacement
        :meth:`Frame.reference_torques` projects the reference forces onto, so the two
        are directly comparable.  Finite differences rather than an analytic gradient
        because this is the reference implementation: it calls :meth:`energy_frame`
        and so cannot drift away from the energy it differentiates.
        ``polyfind.fitting`` carries an analytic version for its inner loop and the
        tests assert the two agree.
        """
        out = np.empty(frame.n_dihedrals)
        for k, (_, b, c, _) in enumerate(np.asarray(frame.torsions, dtype=int)):
            mask = _side_mask(frame.n_atoms, frame.adjacency, int(b), int(c))
            plus = frame.with_coords(rotate_about_bond(frame.coords, mask, int(b), int(c), step))
            minus = frame.with_coords(rotate_about_bond(frame.coords, mask, int(b), int(c), -step))
            out[k] = -(self.energy_frame(plus, torsion) - self.energy_frame(minus, torsion)) / (2 * step)
        return out

    # --- energies -------------------------------------------------------------
    def _energy_from_coords(self, top: _Topology, coords, xp):
        """coords: (M, n, 3) on the active backend -> (M,) energies."""
        pi, pj = xp.asarray(top.pairs_i), xp.asarray(top.pairs_j)
        ri = coords[:, pi]
        rj = coords[:, pj]
        r = xp.sqrt(((ri - rj) ** 2).sum(axis=-1))
        lj_x = xp.asarray(top.lj_x)
        lj_d = xp.asarray(top.lj_d)
        qq = xp.asarray(top.qq)
        s6 = (lj_x / r) ** 6
        e_lj = (lj_d * (s6 * s6 - 2.0 * s6)).sum(axis=1)
        if top.off_atom.size:
            # Off-site charges: the Coulomb term is evaluated between charge *sites*, the
            # Lennard-Jones one between nuclei.  Separating the two is the whole reason the
            # offset exists, so the two distances are genuinely different arrays here.
            cc = self._site_coords(top, coords, xp)
            rc = xp.sqrt(((cc[:, pi] - cc[:, pj]) ** 2).sum(axis=-1))
            e_c = (qq / rc).sum(axis=1)
        else:
            e_c = (qq / r).sum(axis=1)
        e_t = xp.zeros(coords.shape[0])
        if top.torsions.size:
            tors = xp.asarray(top.torsions)
            a, b, c, d = (coords[:, tors[:, k]] for k in range(4))
            b0, b1, b2 = b - a, c - b, d - c
            n1, n2 = xp.cross(b0, b1), xp.cross(b1, b2)
            x = (n1 * n2).sum(-1)
            y = xp.sqrt((b1 * b1).sum(-1)) * (b0 * n2).sum(-1)
            phi = xp.arctan2(y, x)
            # (n_tors,) coefficients, broadcast against phi's (M, n_tors): with one shared
            # triple every row holds the same number, so this is arithmetically identical
            # to the scalar form it replaces (the dtype cast keeps a float32 GPU batch in
            # float32, which a float64 coefficient array would silently promote).
            V = xp.asarray(top.tors_V, dtype=phi.dtype)
            V1, V2, V3 = V[:, 0], V[:, 1], V[:, 2]
            e_t = (0.5 * (V1 * (1 + xp.cos(phi)) + V2 * (1 - xp.cos(2 * phi)) + V3 * (1 + xp.cos(3 * phi)))).sum(axis=1)
        out = e_lj + e_c + e_t
        # The valence blocks are added only when there are any, so a potential without them
        # computes and returns exactly the expression it always did.
        if top.bond_idx.size:
            bi, bj = xp.asarray(top.bond_idx[:, 0]), xp.asarray(top.bond_idx[:, 1])
            rb = xp.sqrt(((coords[:, bi] - coords[:, bj]) ** 2).sum(axis=-1))
            db = rb - xp.asarray(top.bond_r0, dtype=rb.dtype)
            out = out + (0.5 * xp.asarray(top.bond_k, dtype=rb.dtype) * db * db).sum(axis=1)
        if top.angle_idx.size:
            ai, aj, ak = (xp.asarray(top.angle_idx[:, c]) for c in range(3))
            u = coords[:, ai] - coords[:, aj]
            v = coords[:, ak] - coords[:, aj]
            cos = ((u * v).sum(-1) / xp.sqrt((u * u).sum(-1) * (v * v).sum(-1)))
            theta = xp.arccos(xp.clip(cos, -1.0, 1.0))
            da = theta - xp.asarray(top.angle_t0, dtype=theta.dtype)
            out = out + (0.5 * xp.asarray(top.angle_k, dtype=theta.dtype) * da * da).sum(axis=1)
        return out

    def _site_coords(self, top: _Topology, coords, xp):
        """``coords`` with each off-site charge slid along its bond; batched, backend-agnostic."""
        a, nb = xp.asarray(top.off_atom), xp.asarray(top.off_nb)
        u = coords[:, a] - coords[:, nb]
        u = u / xp.sqrt((u * u).sum(axis=-1, keepdims=True))
        out = coords.copy()
        out[:, a] = coords[:, a] + xp.asarray(top.off_d, dtype=coords.dtype)[None, :, None] * u
        return out

    def energy(self, struct: Structure) -> float:
        top = self._topology(struct)
        return float(self._energy_from_coords(top, struct.coords[None], np)[0])

    def energy_coords(self, template: Structure, coords: np.ndarray) -> np.ndarray:
        """Energies for an (M, n, 3) coordinate array using ``template``'s cached topology.

        Same maths as :meth:`_energy_from_coords`/:meth:`energy_batch`, but skips building
        one :class:`Structure` per conformer -- the array typically comes straight from
        :func:`polyfind.chain.build_chain_batch`. Runs on GPU if available.
        """
        coords = np.asarray(coords)
        if coords.shape[0] == 0:
            return np.empty(0)
        xp = bk.get_backend()
        top = self._topology(template)
        coords_xp = xp.asarray(coords, dtype=bk.float_dtype())
        out = np.empty(coords.shape[0])
        chunk = 2048
        for k in range(0, coords.shape[0], chunk):
            out[k : k + chunk] = bk.to_numpy(self._energy_from_coords(top, coords_xp[k : k + chunk], xp))
        return out

    def energy_batch(self, structs: Sequence[Structure] | ConformerBatch) -> np.ndarray:
        """Vectorised over conformers (same topology); runs on GPU if available.

        Accepts either a list of :class:`Structure` or a ``(template, coords)`` tuple
        (see :class:`Calculator`); the latter skips restacking coordinates that are
        already a batched array.
        """
        if isinstance(structs, tuple):
            template, coords = structs
            return self.energy_coords(template, coords)
        structs = list(structs)
        if not structs:
            return np.empty(0)
        coords = np.stack([s.coords for s in structs])
        return self.energy_coords(structs[0], coords)

    def components(self, struct: Structure) -> dict:
        top = self._topology(struct)
        c = struct.coords
        r = np.linalg.norm(c[top.pairs_i] - c[top.pairs_j], axis=1)
        s6 = (top.lj_x / r) ** 6
        lj = float((top.lj_d * (s6 * s6 - 2 * s6)).sum())
        if top.off_atom.size:
            cc = charge_site_coords(c, list(zip(top.off_atom, top.off_nb, top.off_d)))
            coulomb = float((top.qq / np.linalg.norm(cc[top.pairs_i] - cc[top.pairs_j], axis=1)).sum())
        else:
            coulomb = float((top.qq / r).sum())
        bond = angle = 0.0
        if top.bond_idx.size:
            rb = np.linalg.norm(c[top.bond_idx[:, 0]] - c[top.bond_idx[:, 1]], axis=1)
            bond = float((0.5 * top.bond_k * (rb - top.bond_r0) ** 2).sum())
        if top.angle_idx.size:
            u = c[top.angle_idx[:, 0]] - c[top.angle_idx[:, 1]]
            v = c[top.angle_idx[:, 2]] - c[top.angle_idx[:, 1]]
            cos = (u * v).sum(1) / np.sqrt((u * u).sum(1) * (v * v).sum(1))
            angle = float((0.5 * top.angle_k * (np.arccos(np.clip(cos, -1, 1)) - top.angle_t0) ** 2).sum())
        return {
            "lj": lj,
            "coulomb": coulomb,
            "torsion": self.energy(struct) - lj - coulomb - bond - angle,
            "bond": bond,
            "angle": angle,
        }

    def forces_frame(self, frame: Frame, torsion=None, step: float = 1e-5) -> np.ndarray:
        """Cartesian forces ``-dE/dx`` (kcal/(mol A)) of a :class:`Frame`, shape (n, 3).

        Central differences on :meth:`energy_frame`, so this cannot drift away from the
        energy it differentiates: it is the reference implementation the analytic forces in
        :mod:`polyfind.fitting` are tested against, the same arrangement
        :meth:`frame_torques` already has.  With no valence terms these forces are the
        nonbonded-plus-torsion part only -- which is why fitting *Cartesian* forces needed
        the valence terms first, and is the measurement docs/DFT_FIT.md section 3 records.
        """
        out = np.zeros((frame.n_atoms, 3))
        for i in range(frame.n_atoms):
            for a in range(3):
                cp, cm = frame.coords.copy(), frame.coords.copy()
                cp[i, a] += step
                cm[i, a] -= step
                out[i, a] = -(self.energy_frame(frame.with_coords(cp), torsion)
                              - self.energy_frame(frame.with_coords(cm), torsion)) / (2 * step)
        return out


class ASECalculator:
    """Adapter: any ASE calculator (MACE, CHGNet, DFT codes...) as a polyfind Calculator.

    Energies are converted from eV to kcal/mol.  ``energy_batch`` loops by default;
    override for calculators with native batching.
    """

    EV_TO_KCAL = 23.0605

    def __init__(self, ase_calc):
        try:
            from ase import Atoms  # noqa: F401
        except ImportError as e:  # pragma: no cover
            raise ImportError("ASECalculator requires the 'ase' package") from e
        self.calc = ase_calc

    def _atoms_from(self, elements, coords):
        from ase import Atoms

        atoms = Atoms(symbols=elements, positions=coords)
        atoms.calc = self.calc
        return atoms

    def _atoms(self, struct: Structure):
        return self._atoms_from(struct.elements, struct.coords)

    def energy(self, struct: Structure) -> float:
        return float(self._atoms(struct).get_potential_energy()) * self.EV_TO_KCAL

    def energy_batch(self, structs: Sequence[Structure] | ConformerBatch) -> np.ndarray:
        """Loops (no native batching here); accepts the ``(template, coords)`` tuple form
        too, building one ``Atoms`` per row from the shared element list."""
        if isinstance(structs, tuple):
            template, coords = structs
            out = np.empty(len(coords))
            for i, c in enumerate(coords):
                atoms = self._atoms_from(template.elements, c)
                out[i] = float(atoms.get_potential_energy()) * self.EV_TO_KCAL
            return out
        return np.array([self.energy(s) for s in structs])


# ---------------------------------------------------------------------- fitting
def _basins(states: RISStates, grid: np.ndarray) -> list[np.ndarray]:
    """Grid indices belonging to each state's basin (nearest ideal angle on the circle)."""
    diff = np.abs((grid[:, None] - np.array(states.angles)[None, :] + 180.0) % 360.0 - 180.0)
    owner = diff.argmin(axis=1)
    return [np.where(owner == s)[0] for s in range(states.n)]


def _wrap180(x):
    """Wrap an angle (deg, array or scalar) to (-180, 180]."""
    return ((np.asarray(x, dtype=float) + 180.0) % 360.0) - 180.0


def _basin_bounds(states: RISStates, eps: float = 1e-6) -> list[tuple[float, float]]:
    """Per-state ``(lo, hi)`` offsets (deg) from the state's own ideal angle to the
    boundary with its nearest neighbour on each side of the circle, matching
    :func:`_basins`'s tie-break EXACTLY: ``argmin`` over ``|grid - ideal|`` hands an exact
    boundary point to the lowest-index state. So a boundary shared with a lower-indexed
    neighbour belongs to that neighbour, not this state -- this state's interval is
    nudged open on that side by ``eps``; a boundary shared with a higher-indexed
    neighbour belongs to this state, so that side stays closed. Without this, a
    continuous local search can sit exactly on a boundary the discrete (dense-grid)
    bucketing would have assigned to the other state, comparing two different basins.
    """
    angs = np.array(states.angles, dtype=float)
    bounds = []
    for s in range(states.n):
        hi, hi_sp = 180.0, None
        lo, lo_sp = -180.0, None
        for sp in range(states.n):
            if sp == s:
                continue
            d = float(_wrap180(angs[sp] - angs[s]))  # signed shortest distance a_s -> a_sp
            half = d / 2.0
            if d > 0 and half < hi:
                hi, hi_sp = half, sp
            elif d < 0 and half > lo:
                lo, lo_sp = half, sp
        if hi_sp is not None and hi_sp < s:
            hi -= eps
        if lo_sp is not None and lo_sp < s:
            lo += eps
        bounds.append((lo, hi))
    return bounds


def _clip_to_basin(x, ideal: float, lo: float, hi: float):
    """Clip angle(s) ``x`` to ``ideal``'s basin, expressed via the ``(lo, hi)`` offsets
    from :func:`_basin_bounds`; wraps back to a canonical angle."""
    off = np.clip(_wrap180(np.asarray(x, dtype=float) - ideal), lo, hi)
    return _wrap180(ideal + off)


def _quad_vertex(x0, x1, x2, y0, y1, y2, fallback):
    """Vertex of the parabola through three (possibly unequally spaced) points; falls
    back to ``fallback`` if the points are degenerate or the fit is not convex (no
    interior minimum, e.g. a flat or concave stencil)."""
    xs, ys = np.array([x0, x1, x2], dtype=float), np.array([y0, y1, y2], dtype=float)
    if len(set(np.round(xs, 9))) < 3:
        return fallback
    a, b, _ = np.polyfit(xs, ys, 2)
    if not np.isfinite(a) or a <= 1e-12:
        return fallback
    return float(-b / (2.0 * a))


def _refine_coord(polymer, calc, base, ref, j, other_j, other_val, states, bounds, s_idx, center, energy, h, score=None):
    """One coarse-to-fine refinement step along dihedral index ``j`` for a whole batch of
    1-D basins at once (one state each, flattened): evaluate ``center +/- h`` (clipped to
    each point's basin), fit a parabola per point and evaluate its vertex too (also
    clipped), and move to whichever of the four points evaluated has the lowest energy --
    so a bad fit (flat/concave stencil, or a vertex estimate that turns out worse) can
    never make things worse, only fail to help. ``other_j``/``other_val`` hold the other
    dihedral's index and (fixed, per-point) value when this is used for one direction of a
    2-D scan; ``None`` for a plain 1-D scan. Two batched ``energy_batch`` calls cover every
    basin regardless of how many there are. Returns ``(new_center, new_energy,
    n_evaluations)`` with the same shape as ``center``.
    """
    score = _rigid_scorer(polymer, calc) if score is None else score
    shape = center.shape
    c = center.reshape(-1)
    e = energy.reshape(-1)
    sidx = s_idx.reshape(-1)
    other = None if other_val is None else other_val.reshape(-1)
    n = c.shape[0]
    lo = np.array([_clip_to_basin(c[i] - h, states.angles[sidx[i]], *bounds[sidx[i]]) for i in range(n)])
    hi = np.array([_clip_to_basin(c[i] + h, states.angles[sidx[i]], *bounds[sidx[i]]) for i in range(n)])

    dihs = np.tile(base, (2 * n, 1))
    dihs[:n, j] = lo
    dihs[n:, j] = hi
    if other_j is not None:
        dihs[:n, other_j] = other
        dihs[n:, other_j] = other
    E = score(dihs) - ref
    e_lo, e_hi = E[:n], E[n:]

    vtx = np.array([_quad_vertex(lo[i], c[i], hi[i], e_lo[i], e[i], e_hi[i], c[i]) for i in range(n)])
    vtx = np.array([_clip_to_basin(vtx[i], states.angles[sidx[i]], *bounds[sidx[i]]) for i in range(n)])
    dihs_v = np.tile(base, (n, 1))
    dihs_v[:, j] = vtx
    if other_j is not None:
        dihs_v[:, other_j] = other
    e_v = score(dihs_v) - ref

    new_c, new_e = np.empty(n), np.empty(n)
    for i in range(n):
        candidates = [(e[i], c[i]), (e_lo[i], lo[i]), (e_hi[i], hi[i]), (e_v[i], vtx[i])]
        best_e, best_x = min(candidates, key=lambda t: t[0])
        new_c[i], new_e[i] = best_x, best_e
    return new_c.reshape(shape), new_e.reshape(shape), 3 * n


_STENCIL2D = [(dx, dy) for dx in (-1, 0, 1) for dy in (-1, 0, 1) if (dx, dy) != (0, 0)]  # 3x3 minus centre, incl. diagonals


def _refine_pair(polymer, calc, base, ref, j0, j1, states, bounds, s_idx, sp_idx, cx, cy, ce, h, score=None):
    """One coarse-to-fine pattern-search step for a whole batch of 2-D basins (state
    pairs) at once: evaluate the full 3x3 stencil (8 neighbours; the centre is already
    known) and move to whichever of the 9 points has the lowest energy.  The diagonal
    stencil directions let this track a tilted or curved valley (e.g. the pentane-effect
    G+G- basin, where the two dihedrals must move together to avoid a steric clash) in one
    step, unlike two independent 1-D scans; the plain greedy choice (no curve fit) is
    robust to the near-singular repulsive wall right next to that basin's minimum, which
    wrecked a quadratic/Newton fit (huge outlier values, degenerate points when a stencil
    arm clips onto the basin edge). One batched ``energy_batch`` call covers every basin
    regardless of how many there are. Returns ``(new_cx, new_cy, new_ce,
    n_evaluations)``.
    """
    score = _rigid_scorer(polymer, calc) if score is None else score
    shape = cx.shape
    x0, y0, e0 = cx.reshape(-1), cy.reshape(-1), ce.reshape(-1)
    si, pi = s_idx.reshape(-1), sp_idx.reshape(-1)
    n = x0.shape[0]
    n_off = len(_STENCIL2D)
    X = np.empty((n_off, n))
    Y = np.empty((n_off, n))
    for oi, (dx, dy) in enumerate(_STENCIL2D):
        X[oi] = [_clip_to_basin(x0[i] + dx * h, states.angles[si[i]], *bounds[si[i]]) for i in range(n)]
        Y[oi] = [_clip_to_basin(y0[i] + dy * h, states.angles[pi[i]], *bounds[pi[i]]) for i in range(n)]
    dihs = np.tile(base, (n_off * n, 1))
    for oi in range(n_off):
        dihs[oi * n : (oi + 1) * n, j0] = X[oi]
        dihs[oi * n : (oi + 1) * n, j1] = Y[oi]
    E = (score(dihs) - ref).reshape(n_off, n)

    out_x, out_y, out_e = x0.copy(), y0.copy(), e0.copy()
    for oi in range(n_off):
        better = E[oi] < out_e
        out_x[better] = X[oi][better]
        out_y[better] = Y[oi][better]
        out_e[better] = E[oi][better]
    n_eval = n_off * n
    return out_x.reshape(shape), out_y.reshape(shape), out_e.reshape(shape), n_eval


def _polish_pair(polymer, calc, base, ref, j0, j1, states, bounds, s_idx, sp_idx, cx, cy, ce, h, score=None):
    """A single quadratic-surface (Newton) polish on top of :func:`_refine_pair`'s
    pattern-search result: fit ``f0 + g.d + 1/2 d^T H d`` to the same 3x3 stencil, step to
    its stationary point if the fit is locally convex (else stay put), and keep whichever
    of the 10 points evaluated (9 stencil, already-known centre excepted, plus the Newton
    point) is lowest -- so this can only match or improve on the pattern-search centre,
    never lose it. Doing this once, after pattern search has already found the basin's
    true sub-region, is safe where iterating it from a coarse start was not (see
    :func:`_refine_pair`'s docstring): the neighbourhood is by now locally well-behaved, so
    the quadratic fit is no longer at the mercy of the repulsive wall or a boundary-clipped
    duplicate stencil point corrupting its curvature. Returns ``(cx, cy, ce,
    n_evaluations)``.
    """
    score = _rigid_scorer(polymer, calc) if score is None else score
    shape = cx.shape
    x0, y0, e0 = cx.reshape(-1), cy.reshape(-1), ce.reshape(-1)
    si, pi = s_idx.reshape(-1), sp_idx.reshape(-1)
    n = x0.shape[0]
    n_off = len(_STENCIL2D)
    X = np.empty((n_off, n))
    Y = np.empty((n_off, n))
    for oi, (dx, dy) in enumerate(_STENCIL2D):
        X[oi] = [_clip_to_basin(x0[i] + dx * h, states.angles[si[i]], *bounds[si[i]]) for i in range(n)]
        Y[oi] = [_clip_to_basin(y0[i] + dy * h, states.angles[pi[i]], *bounds[pi[i]]) for i in range(n)]
    dihs = np.tile(base, (n_off * n, 1))
    for oi in range(n_off):
        dihs[oi * n : (oi + 1) * n, j0] = X[oi]
        dihs[oi * n : (oi + 1) * n, j1] = Y[oi]
    E = (score(dihs) - ref).reshape(n_off, n)

    new_x, new_y = x0.copy(), y0.copy()
    for i in range(n):
        dx = np.concatenate([[0.0], X[:, i] - x0[i]])
        dy = np.concatenate([[0.0], Y[:, i] - y0[i]])
        de = np.concatenate([[e0[i]], E[:, i]])
        A = np.column_stack([np.ones_like(dx), dx, dy, dx**2, dx * dy, dy**2])
        try:
            coef, *_ = np.linalg.lstsq(A, de, rcond=None)
            _, gx, gy, hxx, hxy, hyy = coef
            H = np.array([[2 * hxx, hxy], [hxy, 2 * hyy]])
            if np.all(np.linalg.eigvalsh(H) > 1e-9):
                delta = np.linalg.solve(-H, [gx, gy])
                new_x[i] = _clip_to_basin(x0[i] + delta[0], states.angles[si[i]], *bounds[si[i]])
                new_y[i] = _clip_to_basin(y0[i] + delta[1], states.angles[pi[i]], *bounds[pi[i]])
        except np.linalg.LinAlgError:
            pass

    dihs_v = np.tile(base, (n, 1))
    dihs_v[:, j0] = new_x
    dihs_v[:, j1] = new_y
    e_v = score(dihs_v) - ref

    out_x, out_y, out_e = np.empty(n), np.empty(n), np.empty(n)
    for i in range(n):
        candidates = [(e0[i], x0[i], y0[i])]
        candidates += [(E[oi, i], X[oi, i], Y[oi, i]) for oi in range(n_off)]
        candidates.append((e_v[i], new_x[i], new_y[i]))
        best_e, best_x, best_y = min(candidates, key=lambda t: t[0])
        out_x[i], out_y[i], out_e[i] = best_x, best_y, best_e
    n_eval = n_off * n + n
    return out_x.reshape(shape), out_y.reshape(shape), out_e.reshape(shape), n_eval


def _dense_scan_bond(polymer, calc, states, base, j0, grid, basins, ref, score=None):
    """Bond-type ``b``'s share of the dense scan (batched conformer building): full
    ``len(grid)`` 1-D scan and ``len(grid)^2`` 2-D scan, basin minima on the grid exactly
    as before. ``score`` maps a ``(M, N)`` dihedral batch to absolute energies; ``None`` is
    the rigid-geometry scorer.  Returns ``(E1, E2, e1, arg1, e2, arg2, n_eval)``."""
    score = _rigid_scorer(polymer, calc) if score is None else score
    S = states.n
    j1 = j0 + 1
    dihs1 = np.tile(base, (len(grid), 1))
    dihs1[:, j0] = grid
    E1 = score(dihs1) - ref
    n_eval = len(grid)
    e1, arg1 = np.zeros(S), {}
    for s, idx in enumerate(basins):
        k = idx[np.argmin(E1[idx])]
        e1[s] = E1[k]
        arg1[states.names[s]] = float(grid[k])

    phis, psis = np.meshgrid(grid, grid, indexing="ij")
    dihs2 = np.tile(base, (len(grid) ** 2, 1))
    dihs2[:, j0] = phis.ravel()
    dihs2[:, j1] = psis.ravel()
    E2 = (score(dihs2) - ref).reshape(len(grid), len(grid))
    n_eval += len(grid) ** 2
    e2, arg2 = np.zeros((S, S)), {}
    for s, idx in enumerate(basins):
        for sp, idxp in enumerate(basins):
            sub = E2[np.ix_(idx, idxp)]
            k = np.unravel_index(np.argmin(sub), sub.shape)
            arg2[(states.names[s], states.names[sp])] = (float(grid[idx[k[0]]]), float(grid[idxp[k[1]]]))
            e2[s, sp] = sub[k]
    return E1, E2, e1, arg1, e2, arg2, n_eval


def _refine_schedule(coarse_step: float, step: float, rounds: int) -> list[float]:
    """Half-widths for the refinement rounds: the first bridges the coarse grid's spacing
    (a coarse-grid pick can be up to ``coarse_step/2`` away from the true basin minimum,
    e.g. when a steric-clash basin's minimum sits near the basin edge -- a stencil narrower
    than that would have to extrapolate rather than interpolate to find it), then each
    round halves down to ``step`` and stays there.
    """
    return [max(step, coarse_step / (2.0 * (r + 1))) for r in range(rounds)]


def _adaptive_scan_bond(polymer, calc, states, base, j0, coarse_grid, basins_coarse, bounds, ref, step, coarse_step, rounds=2, score=None):
    """Bond-type ``b``'s share of the coarse-to-fine scan: a coarse (``coarse_step``) grid
    to find each basin's approximate minimum, batched across every state (1-D) / state
    pair (2-D) at once, then ``rounds`` of refinement (see :func:`_refine_schedule` for the
    half-width schedule): :func:`_refine_coord` (quadratic-interpolation, safe by
    construction) for the 1-D basins; for the 2-D basins the pattern-search
    :func:`_refine_pair` (robust to the steep repulsive wall right next to a
    steric-clash basin's minimum) followed by one :func:`_polish_pair` quadratic step for
    sub-step precision. A couple of batched calculator calls per round cover every basin
    simultaneously. Returns the same shape as :func:`_dense_scan_bond`, with ``E1``/``E2``
    the coarse-grid energies (for diagnostics) and ``e1``/``e2``/``arg1``/``arg2`` the
    refined basin minima."""
    score = _rigid_scorer(polymer, calc) if score is None else score
    S = states.n
    j1 = j0 + 1
    n_eval = 0
    schedule = _refine_schedule(coarse_step, step, rounds)

    # ---- 1-D: coarse grid, then refine every state's basin minimum
    dihs1 = np.tile(base, (len(coarse_grid), 1))
    dihs1[:, j0] = coarse_grid
    E1 = score(dihs1) - ref
    n_eval += len(coarse_grid)
    center1, energy1, s_idx1 = np.empty(S), np.empty(S), np.arange(S)
    for s, idx in enumerate(basins_coarse):
        k = idx[np.argmin(E1[idx])]
        center1[s], energy1[s] = coarse_grid[k], E1[k]
    for h in schedule:
        center1, energy1, ne = _refine_coord(polymer, calc, base, ref, j0, None, None, states, bounds, s_idx1, center1, energy1, h, score=score)
        n_eval += ne
    arg1 = {states.names[s]: float(center1[s]) for s in range(S)}
    e1 = energy1

    # ---- 2-D: coarse grid, then refine every state-pair's basin minimum (phi then psi)
    phis, psis = np.meshgrid(coarse_grid, coarse_grid, indexing="ij")
    dihs2 = np.tile(base, (len(coarse_grid) ** 2, 1))
    dihs2[:, j0] = phis.ravel()
    dihs2[:, j1] = psis.ravel()
    E2 = (score(dihs2) - ref).reshape(len(coarse_grid), len(coarse_grid))
    n_eval += len(coarse_grid) ** 2
    cx, cy, ce = np.empty((S, S)), np.empty((S, S)), np.empty((S, S))
    s_idx_x, s_idx_y = np.empty((S, S), dtype=int), np.empty((S, S), dtype=int)
    for s, idx in enumerate(basins_coarse):
        for sp, idxp in enumerate(basins_coarse):
            sub = E2[np.ix_(idx, idxp)]
            k = np.unravel_index(np.argmin(sub), sub.shape)
            cx[s, sp], cy[s, sp], ce[s, sp] = coarse_grid[idx[k[0]]], coarse_grid[idxp[k[1]]], sub[k]
            s_idx_x[s, sp], s_idx_y[s, sp] = s, sp
    for h in schedule:
        cx, cy, ce, ne = _refine_pair(polymer, calc, base, ref, j0, j1, states, bounds, s_idx_x, s_idx_y, cx, cy, ce, h, score=score)
        n_eval += ne
    # one quadratic polish for sub-step precision, now that pattern search has found the
    # right sub-region (safe: see _polish_pair)
    cx, cy, ce, ne = _polish_pair(polymer, calc, base, ref, j0, j1, states, bounds, s_idx_x, s_idx_y, cx, cy, ce, schedule[-1] / 2.0, score=score)
    n_eval += ne
    arg2 = {(states.names[s], states.names[sp]): (float(cx[s, sp]), float(cy[s, sp])) for s in range(S) for sp in range(S)}
    e2 = ce
    return E1, E2, e1, arg1, e2, arg2, n_eval


# ------------------------------------------------------------ backbone-angle relaxation
#
# The RIS scan rotates one or two backbone bonds of an oligomer whose every other internal
# coordinate is frozen at the polymer's nominal value.  For PVDF that costs little: its
# nominal angles are near the potential's own optimum and the frozen one-bond profile has
# the three canonical minima at every level of relaxation.  For a chain with bulky
# pendants it is not harmless.  With the angles frozen VDCN's profile has *five* minima --
# 180, +/-120 behind a 13 kcal/mol barrier at -8.7, and +/-30 at -13.6 -- and the +/-120
# pair disappears the moment the backbone angles may relax (docs/NITRILE_LANDSCAPE.md,
# section 1).  ``_basins`` hands +/-120 to the trans label, so the rigid fit gave VDCN's T
# state the energy of a well that does not exist, and the screen's "all-trans VDCN ground
# state, rank 1" was that artefact.  The relaxed scan below relaxes every backbone angle of
# the oligomer at each scan point with the driven dihedral(s) held, so the basin minima it
# feeds to the fit are minima of a chain that can bend.

ANGLE_RELAXATION_BOUNDS = (95.0, 135.0)  # deg: the bounds relax_backbone_angles and refine_crystal use
FALLBACK_ANGLE_K = 120.0
"""Bend stiffness, kcal/(mol rad^2) in the ``0.5 k (theta - theta0)^2`` convention, restraining
the relaxed backbone angles to the polymer's nominal values when the calculator has no bend
term of its own.  It is the textbook C-C-C stiffness :data:`polyfind.fitting.VAL_ANGLE_K0`
starts the valence fit from (the fitted value is 77); without *some* restoring term nothing
resists opening an angle and the minimum runs to whichever bound relieves the most contact
(see :func:`relax_backbone_angles`)."""


def relax_angles_batch(polymer: Polymer, calc: Calculator, dihedrals, start=None,
                       lo: float = ANGLE_RELAXATION_BOUNDS[0], hi: float = ANGLE_RELAXATION_BOUNDS[1],
                       restraint_k: float = 0.0, h: float = 0.2, central: bool = True, gtol: float = 3e-3,
                       ftol: float = 1e-6, maxiter: int = 200, step_cap: float = 10.0,
                       stop_above: float | None = None, patience: int = 60) -> tuple:
    """Relax every backbone angle of a batch of oligomers at fixed torsions.

    ``dihedrals`` has shape ``(M, N)``: one conformer per row, and each conformer gets its
    own ``N + 3`` backbone angles -- one per backbone atom, the per-*atom* freedom of
    docs/NITRILE_LANDSCAPE.md section 1 rather than the per-*type* broadcast of
    :func:`relax_backbone_angles`, which cannot bend the chain locally around one rotated
    bond (with the per-type broadcast VDCN's relaxed gauche minimum lands at +/-95 deg; per
    atom it is at +/-40, where the fully relaxed chain has it).  ``start`` is the starting
    angles, ``(N + 3,)`` shared or ``(M, N + 3)`` per row; ``None`` starts from the polymer's
    nominal angles.  Every angle is bounded to ``[lo, hi]``.

    ``restraint_k`` adds ``0.5 k (theta - theta_nominal)^2`` per backbone angle, in
    kcal/(mol rad^2): the stand-in for a bend term when the calculator has none
    (:func:`fit_ris` passes :data:`FALLBACK_ANGLE_K` for a :class:`SimpleFF` without
    ``angle_terms``).  With the fitted valence terms present nothing extra is needed and the
    default of zero leaves the calculator's energy alone.

    Every conformer is minimised at once.  The objective is a sum of independent per-row
    terms, so one batched :func:`~polyfind.chain.build_chain_batch` call per
    finite-difference direction gives the gradient of every row -- ``2 (N + 3) + 1`` batched
    energy calls per iteration whatever ``M`` is (``N + 4`` with ``central=False``) -- while
    each row carries its own projected BFGS direction, Armijo backtracking and
    inverse-Hessian update, so rows converge independently and leave the batch as they do.
    Central differences are the default because the cheaper forward difference is not
    cheaper where it matters: on the stiff Lennard-Jones wall of a G+G- pentane-type
    contact its bias (``h/2`` times a curvature of kcal/mol/deg^2) misdirects BFGS into a
    crawl of hundreds of iterations, where the central stencil converges in under a
    hundred with a fifth of the energies.  ``gtol`` of 3e-3 kcal/mol/deg leaves a row within
    a few 1e-4 kcal/mol of its minimum, far inside anything the 10-degree grid or the
    potential's 0.27 kcal/mol error bar can resolve (1e-3 changes the fitted terms by
    2e-4).  Checked against a per-row ``scipy`` L-BFGS-B on the same objective: the same
    minima to 4e-5 kcal/mol over a 36-point VDCN scan, in less time.

    ``stop_above`` (with ``patience``) is the rule for hopeless rows: a conformer whose
    energy is still above ``stop_above`` after ``patience`` iterations is abandoned where
    it stands.  The rows that need it start ten or more orders of magnitude above the
    reference -- two pendants through each other -- and march down a wall for as many
    iterations as they are given without ever coming near a minimum; :func:`fit_ris` sets
    ``stop_above`` four caps above the relaxed all-trans reference, so an abandoned row is
    clipped to ``cap`` downstream exactly as it would have been.  ``None`` never abandons.

    Returns ``(energies, angles, status, n_energies)``: the relaxed energies ``(M,)``
    (restraint included), the relaxed angles ``(M, N + 3)``, a per-row status (0: the
    projected gradient is below ``gtol`` kcal/mol/deg; 1: stalled -- either no descent step
    is resolvable above the finite-difference noise or an accepted step lowered the energy
    by less than ``ftol`` kcal/mol, which is a minimum to that resolution or a row wedged
    against a wall taking microscopic steps; 2: out of iterations; 3: abandoned above
    ``stop_above``) and the number of calculator energies spent.
    """
    dihs = np.asarray(dihedrals, dtype=float)
    M, N = dihs.shape
    R = N + 3
    B = polymer.bonds_per_repeat
    nominal = np.array([polymer.backbone[k % B].backbone_angle for k in range(R)], dtype=float)
    if start is None:
        x = np.tile(nominal, (M, 1))
    else:
        s = np.asarray(start, dtype=float)
        x = np.tile(s, (M, 1)) if s.ndim == 1 else np.array(s, dtype=float, copy=True)
    x = np.clip(x, lo, hi)
    eye = np.eye(R)
    n_energies = 0

    def energies(dd, ang):
        nonlocal n_energies
        n_energies += len(dd)
        template, coords = build_chain_batch(polymer, dd, bond_angles=ang)
        E = calc.energy_batch((template, coords))
        if restraint_k:
            da = np.deg2rad(ang - nominal)
            E = E + 0.5 * restraint_k * (da * da).sum(axis=1)
        return E

    def value_and_gradient(idx):
        m = len(idx)
        xi = x[idx]
        if central:
            pert = np.empty((2 * R + 1, m, R))
            pert[0] = xi
            for r in range(R):
                pert[1 + 2 * r] = xi + h * eye[r]
                pert[2 + 2 * r] = xi - h * eye[r]
            E = energies(np.tile(dihs[idx], (2 * R + 1, 1)), pert.reshape(-1, R)).reshape(2 * R + 1, m)
            g = np.stack([(E[1 + 2 * r] - E[2 + 2 * r]) / (2.0 * h) for r in range(R)], axis=1)
            return E[0], g
        pert = np.empty((R + 1, m, R))
        pert[0] = xi
        for r in range(R):
            pert[1 + r] = xi + h * eye[r]
        E = energies(np.tile(dihs[idx], (R + 1, 1)), pert.reshape(-1, R)).reshape(R + 1, m)
        g = (E[1:] - E[0][None, :]).T / h
        return E[0], g

    def projected(xi, g):
        pg = g.copy()
        pg[(xi <= lo + 1e-9) & (g > 0)] = 0.0
        pg[(xi >= hi - 1e-9) & (g < 0)] = 0.0
        return pg

    active = np.arange(M)
    E_all, g_all = value_and_gradient(active)
    H0 = 5.0  # deg^2 per kcal/mol; rescaled per row after its first step (Nocedal & Wright 6.20)
    H = np.tile(eye * H0, (M, 1, 1))
    fresh = np.ones(M, dtype=bool)
    status = np.full(M, 2, dtype=int)
    for it in range(maxiter):
        pg = projected(x[active], g_all[active])
        done = np.abs(pg).max(axis=1) < gtol
        status[active[done]] = 0
        active = active[~done]
        if stop_above is not None and it >= patience and len(active):
            hopeless = E_all[active] > stop_above
            status[active[hopeless]] = 3
            active = active[~hopeless]
        if len(active) == 0:
            break
        xi, gi, Ei = x[active], g_all[active], E_all[active]
        pgi = projected(xi, gi)
        d = -np.einsum("mij,mj->mi", H[active], pgi)
        d[(xi <= lo + 1e-9) & (d < 0)] = 0.0  # never push against an active bound
        d[(xi >= hi - 1e-9) & (d > 0)] = 0.0
        slope = (d * gi).sum(axis=1)
        bad = slope >= 0  # not a descent direction: fall back to steepest descent
        if bad.any():
            H[active[bad]] = eye * H0
            fresh[active[bad]] = True
            d[bad] = -pgi[bad]
            slope[bad] = (d[bad] * gi[bad]).sum(axis=1)
        d *= np.minimum(1.0, step_cap / np.maximum(np.abs(d).max(axis=1), 1e-12))[:, None]
        t = np.ones(len(active))
        accepted = np.zeros(len(active), dtype=bool)
        x_new, E_new = xi.copy(), Ei.copy()
        trying = np.arange(len(active))
        for _ls in range(20):
            xt = np.clip(xi[trying] + t[trying, None] * d[trying], lo, hi)
            Et = energies(dihs[active[trying]], xt)
            ok = Et <= Ei[trying] + 1e-4 * ((xt - xi[trying]) * gi[trying]).sum(axis=1)
            x_new[trying[ok]], E_new[trying[ok]], accepted[trying[ok]] = xt[ok], Et[ok], True
            trying = trying[~ok]
            if len(trying) == 0:
                break
            t[trying] *= 0.5
        status[active[~accepted]] = 1
        # an accepted step that lowered the energy by less than ftol is a stall too: a row
        # wedged against a wall can go on taking microscopic Armijo steps for ever
        stagnant = accepted & (Ei - E_new < ftol)
        x[active[stagnant]], E_all[active[stagnant]] = x_new[stagnant], E_new[stagnant]
        status[active[stagnant]] = 1
        accepted &= ~stagnant
        acc = active[accepted]
        if len(acc):
            x[acc], E_all[acc] = x_new[accepted], E_new[accepted]
            E2, g2 = value_and_gradient(acc)
            s = x_new[accepted] - xi[accepted]
            y = g2 - gi[accepted]
            sy, yy = (s * y).sum(axis=1), (y * y).sum(axis=1)
            Hs = H[acc]
            first = fresh[acc] & (sy > 1e-10)
            Hs[first] = eye[None] * (sy[first] / yy[first])[:, None, None]
            fresh[acc] = False
            upd = sy > 1e-10
            if upd.any():
                rho = 1.0 / sy[upd]
                su, yu = s[upd], y[upd]
                A = eye[None] - rho[:, None, None] * su[:, :, None] * yu[:, None, :]
                Hs[upd] = A @ Hs[upd] @ np.transpose(A, (0, 2, 1)) + rho[:, None, None] * su[:, :, None] * su[:, None, :]
            H[acc], g_all[acc] = Hs, g2
        active = acc
    return E_all, x, status, n_energies


def _rigid_scorer(polymer: Polymer, calc: Calculator):
    """``dihedrals (M, N) -> energies (M,)`` at the polymer's frozen backbone angles: the two
    calls every scan of :func:`fit_ris` has always made, in the same order, so the rigid fit
    is bit-for-bit what it was."""

    def score(dihs):
        template, coords = build_chain_batch(polymer, dihs)
        return calc.energy_batch((template, coords))

    return score


class _RelaxedScorer:
    """``dihedrals (M, N) -> energies (M,)`` with every backbone angle relaxed per conformer
    (:func:`relax_angles_batch`), warm-started from the relaxed all-trans reference once
    :meth:`reference` has found it.  Counts what it spends."""

    def __init__(self, polymer, calc, restraint_k, lo, hi, options=None, cap: float = 50.0):
        self.polymer, self.calc = polymer, calc
        self.restraint_k, self.lo, self.hi = float(restraint_k), float(lo), float(hi)
        self.options = dict(options or {})
        self.cap = float(cap)
        self.start = None
        self.stop_above = None  # set by reference(): four caps above the relaxed all-trans energy
        self.n_conformers = 0
        self.n_energies = 0
        self.status = np.zeros(4, dtype=int)  # converged / stalled / out of iterations / abandoned

    def _relax(self, dihs, start, stop_above):
        E, ang, st, n = relax_angles_batch(self.polymer, self.calc, dihs, start=start, lo=self.lo, hi=self.hi,
                                           restraint_k=self.restraint_k, stop_above=stop_above, **self.options)
        self.n_conformers += len(dihs)
        self.n_energies += n
        self.status += np.bincount(st, minlength=4)[:4]
        self.last_angles = ang
        return E, ang

    def __call__(self, dihs, start=None):
        """Relaxed energies; ``start`` (``(M, N + 3)``) warm-starts each row, else the
        relaxed all-trans angles do.  The relaxed angles are left in ``last_angles``."""
        return self._relax(np.asarray(dihs, dtype=float), self.start if start is None else start, self.stop_above)[0]

    def reference(self, base) -> float:
        """Relax the all-trans reference (never abandoned); its angles warm-start every
        later conformer and its energy sets the abandonment threshold."""
        E, ang = self._relax(np.asarray(base, dtype=float)[None], None, None)
        self.start = ang[0]
        self.stop_above = float(E[0]) + 4.0 * self.cap
        return float(E[0])

    def status_counts(self) -> dict:
        return {"converged": int(self.status[0]), "stalled": int(self.status[1]),
                "out_of_iterations": int(self.status[2]), "abandoned": int(self.status[3])}


def _angle_restraint_for(calc: Calculator, angle_restraint: float | None) -> float:
    """The restraint stiffness :func:`fit_ris` relaxes under: the caller's number if given
    (``0`` for none); otherwise none for a :class:`SimpleFF` that carries bend terms, and
    :data:`FALLBACK_ANGLE_K` for one that does not.  A calculator of any other kind (an
    :class:`ASECalculator` wrapping an MLIP or a DFT code) is taken to be a complete
    potential with its own bending physics and gets no restraint."""
    if angle_restraint is not None:
        return float(angle_restraint)
    if isinstance(calc, SimpleFF):
        return 0.0 if calc.angle_table() else FALLBACK_ANGLE_K
    return 0.0


def _local_minima(E, min_depth: float = 0.0) -> list[tuple[int, float]]:
    """``(index, depth)`` of every local minimum of a cyclic profile at least ``min_depth``
    deep, the depth being the lower of the two barrier tops enclosing it.  A plateau of
    equal neighbours counts once (its first point)."""
    E = np.asarray(E, dtype=float)
    n = len(E)
    out = []
    for i in range(n):
        if not (E[i] < E[i - 1] and E[i] <= E[(i + 1) % n]):
            continue
        k = i
        while E[(k - 1) % n] >= E[k % n] and i - k < n:
            k -= 1
        left = E[k % n]
        k = i
        while E[(k + 1) % n] >= E[k % n] and k - i < n:
            k += 1
        right = E[k % n]
        depth = min(left, right) - E[i]
        if depth >= min_depth:
            out.append((i, float(depth)))
    return out


def _watershed_basins(states: RISStates, grid: np.ndarray, E) -> tuple[list[np.ndarray], list[bool]]:
    """Grid indices per state by descent on the one-bond profile ``E``, for the relaxed scan.

    Every local minimum of the profile is assigned to the state whose ideal angle is
    nearest (ties to the lowest-indexed state, as :func:`_basins` breaks them), and every
    grid point to the minimum it runs down to, so a state's basin is bounded by the
    profile's own barrier tops rather than by the midpoints between ideal angles.  The
    difference matters on a relaxed surface: VDCN's relaxed trans well is 1.2 kcal/mol deep
    with its barriers at +/-150, and beyond them the profile is already descending into
    the +/-40 gauche well, so the nearest-ideal T basin ``[120, 240]`` has its *lowest grid
    point* at its edge (+/-120, -0.44) rather than at its minimum (180, 0.0); read off that
    way the "trans" state lands at -120 deg carrying the slope's energy, which is the same
    basin-edge-for-state error the rigid scan makes with a different well.  By watershed
    each state's energy is the energy of a minimum of the surface, or nothing.

    A state with no minimum of its own (AN's trans, a monotonic descent at either level)
    keeps its nearest-ideal basin, whose lowest point is a basin edge; the second return
    value flags those states so the report can say so.  The rigid scan does not use this
    -- its numbers are pinned bit-for-bit -- and on PVDF the two assignments agree.
    """
    E = np.asarray(E, dtype=float)
    n = len(grid)
    minima = [i for i, _ in _local_minima(E)]
    ideal = np.array(states.angles, dtype=float)
    owner_of_min = {i: int(np.argmin(np.abs(_wrap180(grid[i] - ideal)))) for i in minima}
    listed = np.zeros(n, dtype=bool)
    listed[minima] = True
    end = np.empty(n, dtype=int)
    for i in range(n):
        k = i
        for _ in range(n):
            if listed[k]:
                break
            lo, hi = (k - 1) % n, (k + 1) % n
            # step to the lower neighbour; on a tie prefer a listed minimum, then the left
            e, _, nxt = min((E[lo], 0 if listed[lo] else 1, lo), (E[hi], 0 if listed[hi] else 1, hi))
            if e > E[k]:
                break
            k = nxt
        end[i] = k
    fallback = _basins(states, grid)
    basins, edge = [], []
    for s in range(states.n):
        mine = [i for i in minima if owner_of_min[i] == s]
        if mine:
            basins.append(np.where(np.isin(end, mine))[0])
            edge.append(False)
        else:
            basins.append(fallback[s])
            edge.append(True)
    return basins, edge


def _dense_scan_relaxed(polymer, states, base, N, grid, ref, score: "_RelaxedScorer"):
    """The dense 1-D and 2-D scans of every bond type on the relaxed surface, with
    watershed basins (:func:`_watershed_basins`) read off each bond type's own 1-D
    profile: the 2-D basin of a state pair on bonds ``(j0, j0 + 1)`` is the product of
    bond type ``b``'s basin for the first state and bond type ``(b + 1) % B``'s for the
    second, which is why every 1-D profile is taken before any 2-D scan.  Each 2-D row is
    warm-started from the relaxed angles of its 1-D row at the same first dihedral.
    Returns ``(scan1, scan2, e1, arg1, e2, arg2, n_eval, edge_states)`` in
    :func:`_dense_scan_bond`'s conventions plus, per bond type, the names of the states
    whose "minimum" is a basin edge."""
    B, S = polymer.bonds_per_repeat, states.n
    j0s = [(N // 2 - 1) - ((N // 2 - 1) - b) % B for b in range(B)]
    scan1, angles1, basins, edges = {}, {}, {}, {}
    n_eval = 0
    for b in range(B):
        dihs1 = np.tile(base, (len(grid), 1))
        dihs1[:, j0s[b]] = grid
        scan1[b] = score(dihs1) - ref
        angles1[b] = score.last_angles
        n_eval += len(grid)
        basins[b], edge = _watershed_basins(states, grid, scan1[b])
        edges[b] = [states.names[s] for s in range(S) if edge[s]]
    e1, e2 = np.zeros((B, S)), np.zeros((B, S, S))
    arg1, arg2, scan2 = {}, {}, {}
    for b in range(B):
        E1 = scan1[b]
        arg1[b] = {}
        for s, idx in enumerate(basins[b]):
            k = idx[np.argmin(E1[idx])]
            e1[b, s] = E1[k]
            arg1[b][states.names[s]] = float(grid[k])
        j0 = j0s[b]
        phis, psis = np.meshgrid(grid, grid, indexing="ij")
        dihs2 = np.tile(base, (len(grid) ** 2, 1))
        dihs2[:, j0] = phis.ravel()
        dihs2[:, j0 + 1] = psis.ravel()
        start = np.repeat(angles1[b], len(grid), axis=0)  # row (phi_i, psi_j) starts from the 1-D row phi_i
        E2 = (score(dihs2, start=start) - ref).reshape(len(grid), len(grid))
        n_eval += len(grid) ** 2
        bn = (b + 1) % B
        arg2[b] = {}
        for s, idx in enumerate(basins[b]):
            for sp, idxp in enumerate(basins[bn]):
                sub = E2[np.ix_(idx, idxp)]
                k = np.unravel_index(np.argmin(sub), sub.shape)
                arg2[b][(states.names[s], states.names[sp])] = (float(grid[idx[k[0]]]), float(grid[idxp[k[1]]]))
                e2[b, s, sp] = sub[k]
        scan2[b] = E2
    return scan1, scan2, e1, arg1, e2, arg2, n_eval, edges


def angle_relaxation_check(polymer: Polymer, calc: Calculator, step: float = 10.0, n_monomers: int = 5,
                           tol: float = 30.0, min_depth: float = 0.05, angle_restraint: float | None = None,
                           angle_bounds=ANGLE_RELAXATION_BOUNDS, relax_options: dict | None = None) -> dict:
    """Does the rigid one-bond profile have minima the angle-relaxed one lacks?

    The discriminator behind :data:`ANGLE_RELAXATION_DEFAULTS` and behind
    ``fit_ris(angles="auto")`` for a polymer not recorded there.  For each bond type it
    takes the same one-bond scan :func:`fit_ris` does (``step`` degrees round the circle,
    every other bond trans, the oligomer ``n_monomers`` long) twice -- once at the polymer's
    frozen backbone angles and once with every backbone angle relaxed per point
    (:func:`relax_angles_batch`, the same restraint and bounds :func:`fit_ris` would use) --
    finds the local minima of each (at least ``min_depth`` kcal/mol deep, so grid noise is
    not a well), and calls a rigid minimum an **orphan** when no relaxed minimum lies within
    ``tol`` degrees of it.  A chemistry needs the relaxed scan when it has an orphan on any
    bond type: that is a well the frozen geometry manufactures, and ``_basins`` would fold
    its energy into whichever state it lies nearest.

    ``tol`` of 30 deg is the kink criterion of docs/NITRILE_LANDSCAPE.md.  It is wide
    enough that PVDF's rigid double gauche well (+/-70 and +/-40, both within 30 of the
    relaxed +/-50) is not an orphan, and narrow enough that VDCN's +/-120 (80 deg from the
    relaxed +/-40, 60 from 180) is.

    Returns ``{"needs_relaxation": bool, "grid": ..., "bonds": {b: {...}}, "n_evaluations",
    "n_energies", "restraint_k"}``; each bond's entry carries both profiles (relative to
    their own all-trans), both minima lists as ``(angle, energy, depth)`` and the orphans.
    """
    B = polymer.bonds_per_repeat
    N = max(n_monomers * B, 2 * B + 6)
    base = np.full(N, 180.0)
    grid = np.arange(-180.0, 180.0, step)
    rigid = _rigid_scorer(polymer, calc)
    relaxed = _RelaxedScorer(polymer, calc, _angle_restraint_for(calc, angle_restraint),
                             angle_bounds[0], angle_bounds[1], relax_options)
    ref_rigid = calc.energy(build_chain(polymer, base))
    ref_relaxed = relaxed.reference(base)
    out = {"needs_relaxation": False, "grid": grid, "bonds": {}, "n_evaluations": 1 + 1,
           "restraint_k": relaxed.restraint_k, "relaxed_reference_angles": relaxed.start.copy()}
    for b in range(B):
        j0 = (N // 2 - 1) - ((N // 2 - 1) - b) % B
        dihs = np.tile(base, (len(grid), 1))
        dihs[:, j0] = grid
        E_rigid = rigid(dihs) - ref_rigid
        E_relaxed = relaxed(dihs) - ref_relaxed
        out["n_evaluations"] += 2 * len(grid)
        m_rigid = _local_minima(E_rigid, min_depth)
        m_relaxed = _local_minima(E_relaxed, min_depth)
        orphans = [(float(grid[i]), float(E_rigid[i]), d) for i, d in m_rigid
                   if all(abs(_wrap180(grid[i] - grid[k])) > tol for k, _ in m_relaxed)]
        out["bonds"][b] = {
            "rigid_profile": E_rigid, "relaxed_profile": E_relaxed,
            "rigid_minima": [(float(grid[i]), float(E_rigid[i]), d) for i, d in m_rigid],
            "relaxed_minima": [(float(grid[i]), float(E_relaxed[i]), d) for i, d in m_relaxed],
            "orphans": orphans,
        }
        if orphans:
            out["needs_relaxation"] = True
    out["n_energies"] = relaxed.n_energies + out["n_evaluations"]
    return out


# Which registered chemistries need the relaxed scan by default, measured rather than
# assigned: :func:`angle_relaxation_check` under the fitted potential (``pvdf-dft-valence``,
# step 10 deg, six monomers -- the screen's settings; ``examples/angle_relaxation_defaults.py``
# reproduces the table below and ``tests/test_forcefield.py`` re-measures the PVDF and VDCN
# rows).  A chemistry needs it when its rigid one-bond profile has a minimum no relaxed
# minimum lies within 30 deg of.  Bond type 0; minima are deg: kcal/mol relative to each
# level's own all-trans; bond type 1 gives the same verdict in every row.
#
#   polymer  rigid minima                                relaxed minima                  orphans     default
#   pvdf     180: 0; +/-70: -1.87; +/-40: -1.18           180: 0; +/-50: -1.79            none        rigid
#   pe       180: 0; +/-70: -0.39                         180: 0; +/-70: +0.01            none        rigid
#   pvdc     180: 0; +/-120: -2.44; +/-30: -11.29         180: 0; +/-30: -4.55            +/-120      relaxed
#   cfe      -140: -20.3; -40: -24.1; +30: -12.0; +80: -11.9   -150: -2.4; -40: -5.2; +40: -3.2   +80    relaxed
#   cdfe     -110: -14.2; -30: +68.5; +70: -25.0; +150: -19.2  -100: -4.5; +60: -11.2        -30, +150   relaxed
#   an       -60: -11.9; +40: -6.6; +90: -10.0            -70: -5.2; +90: -3.6            +40         relaxed
#   vdcn     180: 0; +/-120: -8.66; +/-30: -13.56         180: 0; +/-40: -3.86            +/-120      relaxed
#   fanome   (overlapping all-trans reference: every energy is -2e4 to -3e4 and means nothing; recorded
#             relaxed because the rule says so, and the screen does not fit it)             -80         relaxed
#
# Under the unfitted ``SimpleFF()`` with the :data:`FALLBACK_ANGLE_K` restraint the same
# rule gives the same verdict for every row except CDFE and FANOME (rigid there), so a test
# that fits CDFE with the illustrative potential relaxes angles it would not strictly need
# to; that costs time, not correctness.  ``vdf11-vdcn1`` is not recorded (nothing fits its
# 24-bond repeat directly; ``transfer_ris`` assembles it), so ``"auto"`` would measure it.
ANGLE_RELAXATION_DEFAULTS: dict[str, bool] = {
    "pvdf": False,
    "pe": False,
    "pvdc": True,
    "cfe": True,
    "cdfe": True,
    "an": True,
    "vdcn": True,
    "fanome": True,
}


@dataclass
class FitReport:
    grid: np.ndarray
    scan1: dict  # bond type -> (n_grid,) energies relative to all-trans
    scan2: dict  # bond type -> (n_grid, n_grid) energies relative to all-trans
    n_evaluations: int
    model: RISModel
    argmin1: dict  # bond type -> {state: angle at basin minimum}
    argmin2: dict  # bond type -> {(s, s'): (angle, angle)}
    angles: str = "rigid"  # "rigid" or "relaxed": how the backbone angles were treated
    relaxed_reference_angles: np.ndarray | None = None  # (N+3,) relaxed all-trans angles, relaxed scans only
    n_relaxation_energies: int = 0  # calculator energies spent inside the relaxations (0 when rigid)
    relaxation_status: dict | None = None  # counts of converged / stalled / out-of-iterations rows
    discriminator: dict | None = None  # angle_relaxation_check's result when "auto" had to measure
    third_order_guard: dict | None = None  # _fit_third_order's "capped" / "zeroed" masks and "raw" values
    edge_states: dict | None = None  # relaxed scans: bond type -> states whose "minimum" is a basin edge



def _type_shift(t: np.ndarray, k: int, B: int) -> np.ndarray:
    """``t[(k - b) % B]`` along the leading bond-type axis (the reversal's type map)."""
    return t[[(k - b) % B for b in range(B)]]


def _reversal_images(e1: np.ndarray, e2: np.ndarray, mirror, B: int):
    """Image of a fitted model under reflection composed with chain reversal.

    Reflection alone maps a chain to its enantiomer, so it is a symmetry only of an
    achiral repeat.  Reversing the chain direction swaps each stereocentre's
    neighbours and undoes that, so reflection-with-reversal is a symmetry of *any*
    linear repeat: E(phi_1..phi_N) == E(-phi_N..-phi_1).

    Writing the reversal as ``j -> K - j`` on the bonds, a term over ``n`` consecutive
    bonds starting at ``j`` maps to the term starting at ``K - n + 1 - j``, so each
    order gets its own bond-type shift, one lower than the order below it, together
    with the state mirror and a reversal of the term's own index order:

        e1[b, s]         == e1[(c - b) % B, m(s)]
        e2[b, s, s']     == e2[(c - 1 - b) % B, m(s'), m(s)]
        e3[b, s, s', s'']== e3[(c - 2 - b) % B, m(s''), m(s'), m(s)]

    ``c`` depends on where the repeat is phased, so every shift is tried and the one
    with the smallest residual is returned along with that residual and ``c`` itself;
    the caller applies the averaging only if the residual is at grid-noise level.
    For ``B = 2`` the shift is decided unambiguously (measured residuals 0.03-0.09
    kcal/mol at ``c = 1`` against 14-58 at ``c = 0``, for PVDF, PVDC, CFE and CDFE),
    and ``c = 1`` gives the swap-mirror on ``e1``, the plain transpose-mirror on
    ``e2`` and the swap-mirror-reverse on ``e3``.

    Returns ``(e1_img, e2_img, resid, c)``; :func:`_reversal_image3` builds the
    matching third-order image once ``e3`` has been fitted.
    """
    m = np.asarray(mirror)
    best = None
    for c in range(B):
        e1_img = _type_shift(e1, c, B)[:, m]
        e2_img = np.transpose(_type_shift(e2, c - 1, B)[:, m][:, :, m], (0, 2, 1))
        resid = max(float(np.abs(e1 - e1_img).max()), float(np.abs(e2 - e2_img).max()))
        if best is None or resid < best[2]:
            best = (e1_img, e2_img, resid, c)
    return best


def _reversal_image3(e3: np.ndarray, mirror, B: int, c: int):
    """Third-order image under reflection-with-reversal, and its residual.

    ``e3_img[b, x, y, z] = e3[(c - 2 - b) % B, m(z), m(y), m(x)]`` -- the state mirror,
    the triple read backwards (reversal reverses the order of a triple as well as
    mirroring its states) and a bond-type shift two below the first-order one.  The
    caller checks the residual before averaging over it.
    """
    m = np.asarray(mirror)
    img = np.transpose(_type_shift(e3, c - 2, B)[:, m][:, :, m][:, :, :, m], (0, 3, 2, 1))
    return img, float(np.abs(e3 - img).max())


def _reversal_angle_residual(argmin1: dict, states: RISStates, B: int, c: int) -> float:
    """How far the fitted basin minima are from ``arg1[b, s] == -arg1[(c - b) % B, m(s)]``.

    The 1-D scan of a bond of type ``b`` maps, under reflection-with-reversal, onto the
    scan of a bond of type ``(c - b) % B`` at the negated angle, so the basin minima of
    ``s`` and ``m(s)`` on the two types are negatives of each other.  Only the states
    :func:`fit_ris`'s ``adapt_angles`` averaging actually moves -- those with
    ``m(s) != s`` -- are checked: a self-mirror state (T) can sit in a symmetric double
    well whose two minima are degenerate, where the reported argmin is a tie-break
    rather than a measurement.  PVDC's trans basin is one: its planar zigzag is not even
    metastable, so the minimum runs out to both basin edges, at -120 and +120 deg, which
    differ by 4e-13 kcal/mol and are resolved opposite ways on the two bond types.
    Degrees.
    """
    worst = 0.0
    for b in range(B):
        for s, nm in enumerate(states.names):
            if states.mirror[s] == s:
                continue
            d = argmin1[b][nm] + argmin1[(c - b) % B][states.names[states.mirror[s]]]
            worst = max(worst, abs((d + 180.0) % 360.0 - 180.0))
    return worst


def fit_ris(
    polymer: Polymer,
    calc: Calculator,
    states: RISStates = THREE_STATE,
    step: float = 10.0,
    n_monomers: int = 5,
    name: str | None = None,
    symmetrize: bool | str = "auto",
    adapt_angles: bool = True,
    third_order: bool = False,
    cap: float = 50.0,
    reversal_tol: float = 0.5,
    angle_tol: float | None = None,
    scan: str = "dense",
    coarse_step: float = 30.0,
    angles: str = "auto",
    angle_restraint: float | None = None,
    angle_bounds=ANGLE_RELAXATION_BOUNDS,
    relax_options: dict | None = None,
) -> FitReport:
    """Derive an RIS model from dihedral scans of a short oligomer.

    First-order energies come from a 1D scan of one bond of each type (all other
    bonds trans); pair energies from a 2D scan of two consecutive bonds, minus the
    first-order terms.  Energies are basin minima, relative to all-trans.

    ``angles`` says what the backbone angles do while a bond is driven:

    * ``"rigid"``: frozen at the polymer's nominal values -- the scan as it always was,
      bit-for-bit (``scan="rigid"`` is accepted as shorthand for this with the dense grid).
    * ``"relaxed"``: every backbone angle of the oligomer relaxed per scan point with the
      driven dihedral(s) held (:func:`relax_angles_batch`; bounded to ``angle_bounds``),
      the all-trans reference relaxed the same way, and the basin assignment, the adapted
      state angles and the second- and third-order terms all read off the relaxed scan
      exactly as they read off the rigid one.  Under a :class:`SimpleFF` carrying bend
      terms (the ``pvdf-dft-valence`` preset) the relaxation is against those terms; under
      one without, a harmonic restraint of :data:`FALLBACK_ANGLE_K` about the nominal
      angles stands in (``angle_restraint`` overrides either way; ``0`` means none; any
      other calculator is trusted to have its own bending physics).  ``scan="relaxed"`` is
      shorthand for this with the dense grid.  Basins are assigned by descent on each bond
      type's relaxed profile (:func:`_watershed_basins`) rather than by nearest ideal angle,
      and the adaptive grid is refused.  It costs 9-13 minutes for a two-bond-type chemistry
      at the screen's settings (step 10 deg, six monomers, third order; measured with eight
      fits sharing twelve cores) against one second rigid -- every one of the ~2900 scan
      conformers is a bounded minimisation over its ``N + 3`` angles, 1.3-1.8 million
      calculator energies in all -- and ``relax_options`` (``h``, ``central``, ``gtol``,
      ``ftol``, ``maxiter``, ``step_cap``, ``patience``) is passed through.
    * ``"auto"`` (default): the recorded answer in :data:`ANGLE_RELAXATION_DEFAULTS` for a
      registered chemistry; for one not recorded there, :func:`angle_relaxation_check`
      is run first (a rigid and a relaxed one-bond scan per bond type, at this ``step``
      and ``n_monomers``) and its verdict used, kept in the report's ``discriminator``.

    Why this exists: with the angles frozen a chain with bulky pendants has one-bond
    minima that are not minima of a chain that can bend.  VDCN's rigid profile has wells
    at +/-120 deg, 8.7 kcal/mol *below* trans, that vanish the moment the angles relax;
    ``_basins`` hands them to the trans label, so the rigid fit gave VDCN a T state at 180
    deg carrying the energy of a well at 120, and every ranking built on it inherited that
    (docs/NITRILE_LANDSCAPE.md).  PVDF's profile has the same three minima at either level,
    which is why its rigid numbers, and every table built on them, stand.

    With
    ``symmetrize`` averages the model with its mirror image (G+ <-> G-).  That is
    exact for an achiral chain and removes grid/refinement noise, but it is *wrong*
    for a chiral one: reflecting a chiral chain gives its enantiomer, not the same
    molecule, so G+ and G- genuinely differ and averaging destroys a real asymmetry
    (measured at 24 kcal/mol for CFE and 45 for CDFE).  The default ``"auto"``
    therefore mirrors only when ``polymer.is_chiral`` is false; ``True`` and
    ``False`` force it either way, and ``True`` on a chiral polymer warns.  The
    symmetry that does survive reflection for a chiral chain is reflection composed
    with chain reversal, E(phi_1..phi_N) == E(-phi_N..-phi_1), and ``"auto"`` averages
    over that instead for a chiral polymer: a state mirror with a bond-type shift on
    the first-order term, a state mirror with a transpose on the pair term, and a state
    mirror with the triple read backwards on the third-order term (see
    :func:`_reversal_images`).  Each is checked numerically before use (every bond-type
    shift is tried, and the averaging is applied only if the residual is within
    ``reversal_tol``), so none can be applied where it does not hold.  Measured
    residuals at step = 20 deg are 0.03-0.09 kcal/mol on the first and second order and
    0.00-0.05 on the third for PE, PVDF, PVDC, CFE and CDFE alike, against 32-47 (pair)
    and 100 (triple, i.e. the ``cap``) for plain mirroring of the chiral pair.
    With ``adapt_angles``
    the state dihedral angles are moved to the 1D-scan basin minima (averaged over
    bond types and mirror pairs), as in classical RIS parametrisations.  Forcing a
    mirror pair to equal magnitude is valid under reflection-with-reversal too -- it
    follows from the same relation applied to the basin minima -- and is likewise
    validated first, against ``angle_tol`` degrees (default: the scan resolution
    ``step``, which is the quantum the basin minima are located to; measured 0.00 deg
    on a dense grid and <= 0.1 deg with ``scan="adaptive"``, against 80 deg for plain
    mirroring of CFE and CDFE).  With
    ``third_order`` triplet corrections are added (see :func:`_fit_third_order`),
    which is what distinguishes e.g. TG+TG+ (3/1 helix) from TG+TG- (alpha-PVDF).
    Energies above ``cap`` (steric overlap) are clipped to ``cap``.

    Every conformer of every scan is built in one batched :func:`~polyfind.chain.build_chain_batch`
    call and evaluated in one ``calc.energy_batch`` call per scan (or refinement round), so an
    expensive (MLIP/DFT) calculator sees few, large batches rather than many single points.

    ``scan`` selects the search strategy:

    * ``"dense"`` (default): the original ``step``-degree grid over the full 360 deg,
      basin minima read off the grid -- exact and reproducible, but most points are far
      from any basin.
    * ``"adaptive"``: a coarse (``coarse_step``, default 30 deg) grid locates each basin,
      then a bounded local refinement (not restricted to a grid) polishes each basin's
      minimum at ``step`` resolution, for several times fewer evaluations.

    Dense is the default because, once batched, it is already cheap with the built-in
    potential (0.29 s vs. adaptive's 0.22 s for a full PVDF fit) -- adaptive's ~2.8x fewer
    evaluations only pays for itself when the calculator is expensive (an MLIP or DFT).
    In exchange, adaptive is not confined to the grid at all: where a state pair's basin
    has no interior minimum (its energy just falls monotonically towards a neighbouring
    state), adaptive reports the true basin-edge infimum rather than the grid's best
    sample, which can materially change that pair's fitted energy. So treat ``"adaptive"``
    as an informed opt-in for expensive calculators, not a drop-in speedup.
    """
    B = polymer.bonds_per_repeat
    N = max(n_monomers * B, 2 * B + 6)
    base = np.full(N, 180.0)
    if scan in ("rigid", "relaxed"):
        scan, angles = "dense", scan
    if angles not in ("auto", "rigid", "relaxed"):
        raise ValueError(f"angles must be 'auto', 'rigid' or 'relaxed'; got {angles!r}")
    discriminator = None
    if angles == "auto":
        if polymer.name in ANGLE_RELAXATION_DEFAULTS:
            angles = "relaxed" if ANGLE_RELAXATION_DEFAULTS[polymer.name] else "rigid"
        else:
            discriminator = angle_relaxation_check(polymer, calc, step=step, n_monomers=n_monomers,
                                                  angle_restraint=angle_restraint, angle_bounds=angle_bounds,
                                                  relax_options=relax_options)
            angles = "relaxed" if discriminator["needs_relaxation"] else "rigid"
    relaxer = None
    if angles == "relaxed":
        relaxer = _RelaxedScorer(polymer, calc, _angle_restraint_for(calc, angle_restraint),
                                 angle_bounds[0], angle_bounds[1], relax_options, cap=cap)
        ref = relaxer.reference(base)
        score = relaxer
    else:
        ref = calc.energy(build_chain(polymer, base))
        score = _rigid_scorer(polymer, calc)
    n_eval = 1 + (0 if discriminator is None else discriminator["n_evaluations"])
    e1 = np.zeros((B, states.n))
    e2 = np.zeros((B, states.n, states.n))
    scan1, scan2, arg1, arg2 = {}, {}, {}, {}
    edge_states = None
    if scan == "dense" and angles == "relaxed":
        report_grid = np.arange(-180.0, 180.0, step)
        scan1, scan2, e1, arg1, e2, arg2, ne, edge_states = _dense_scan_relaxed(polymer, states, base, N, report_grid, ref, relaxer)
        n_eval += ne
    elif scan == "adaptive" and angles == "relaxed":
        raise ValueError("scan='adaptive' cannot be combined with angles='relaxed': the adaptive search "
                         "confines each state to its nearest-ideal-angle basin, whose edge on a relaxed "
                         "surface can lie on the slope into the neighbouring well (see _watershed_basins); "
                         "use the dense grid for a relaxed fit")
    elif scan == "dense":
        report_grid = np.arange(-180.0, 180.0, step)
        basins = _basins(states, report_grid)
        for b in range(B):
            # scanned bond j0 of type b, nearest the middle of the oligomer
            j0 = (N // 2 - 1) - ((N // 2 - 1) - b) % B
            E1, E2, e1_b, arg1_b, e2_b, arg2_b, ne = _dense_scan_bond(polymer, calc, states, base, j0, report_grid, basins, ref, score=score)
            scan1[b], scan2[b] = E1, E2
            e1[b], arg1[b] = e1_b, arg1_b
            e2[b], arg2[b] = e2_b, arg2_b
            n_eval += ne
    elif scan == "adaptive":
        # offset by half a step so grid points never land exactly on a basin boundary
        # (they would if coarse_step divides the state spacing evenly, e.g. 30 deg into
        # a 3-state model's 120 deg spacing); a boundary tie excludes the outer edge of
        # every basin from the coarse sample, which is exactly where a steric-clash
        # basin's true minimum can sit, right next to the repulsive wall.
        report_grid = np.arange(-180.0 + coarse_step / 2.0, 180.0, coarse_step)
        basins_coarse = _basins(states, report_grid)
        bounds = _basin_bounds(states)
        for b in range(B):
            j0 = (N // 2 - 1) - ((N // 2 - 1) - b) % B
            E1, E2, e1_b, arg1_b, e2_b, arg2_b, ne = _adaptive_scan_bond(
                polymer, calc, states, base, j0, report_grid, basins_coarse, bounds, ref, step, coarse_step, score=score
            )
            scan1[b], scan2[b] = E1, E2
            e1[b], arg1[b] = e1_b, arg1_b
            e2[b], arg2[b] = e2_b, arg2_b
            n_eval += ne
    else:
        raise ValueError(f"unknown scan mode {scan!r}; expected 'dense' or 'adaptive'")
    # remove first-order contributions from the pair scan (the 1D scans were done with the
    # partner trans, so e1 already includes the state's interaction with a trans partner)
    for b in range(B):
        bn = (b + 1) % B
        e2[b] = e2[b] - e1[b][:, None] - e1[bn][None, :]
    chiral = bool(getattr(polymer, "is_chiral", False))
    mode = symmetrize
    if mode == "auto":
        mode = "reversal" if chiral else "mirror"
    elif mode is True:
        mode = "mirror"
        if chiral:
            warnings.warn(
                f"symmetrize=True on chiral polymer {polymer.name!r}: reflection maps the chain "
                "to its enantiomer, so averaging G+ with G- destroys a real energy difference. "
                "Use symmetrize='auto' (the default), which averages over reflection-with-reversal "
                "instead, or False.",
                UserWarning,
                stacklevel=2,
            )
    elif mode is False:
        mode = "none"
    if mode not in ("mirror", "reversal", "none"):
        raise ValueError(f"symmetrize must be 'auto', 'mirror', 'reversal', True or False; got {symmetrize!r}")
    rev_shift = 0
    if mode == "mirror":
        m = np.array(states.mirror)
        e1 = 0.5 * (e1 + e1[:, m])
        e2 = 0.5 * (e2 + e2[:, m][:, :, m])
    elif mode == "reversal":
        e1_img, e2_img, resid, rev_shift = _reversal_images(e1, e2, states.mirror, B)
        if resid <= reversal_tol:
            e1 = 0.5 * (e1 + e1_img)
            e2 = 0.5 * (e2 + e2_img)
        else:
            warnings.warn(
                f"reflection-with-reversal does not hold for {polymer.name!r} to within "
                f"{reversal_tol} kcal/mol (residual {resid:.3f}); leaving the fit unsymmetrised.",
                UserWarning,
                stacklevel=2,
            )
            mode = "none"
    e1 = np.minimum(e1, cap)
    e2 = np.minimum(e2, cap)
    if adapt_angles:
        angs = []
        for s, nm in enumerate(states.names):
            vals = np.array([arg1[b][nm] for b in range(B)])
            # circular mean of the basin minima across bond types
            ang = np.degrees(np.arctan2(np.sin(np.deg2rad(vals)).mean(), np.cos(np.deg2rad(vals)).mean()))
            angs.append(ang)
        angs = np.array(angs)
        # Mirror pairs get equal magnitude and opposite sign.  Under plain mirroring that
        # is immediate.  Under reflection-with-reversal it follows from the per-bond-type
        # relation arg1[b, s] == -arg1[(c - b) % B, m(s)]: averaging it over b (a bijection
        # of the bond types) gives angle(m(s)) == -angle(s) for the aggregated angles too.
        # So the same averaging is right for a chiral fit -- but only if that relation
        # actually holds, which is measured here rather than assumed.
        pair_angles = mode == "mirror"
        if mode == "reversal":
            ang_tol = step if angle_tol is None else angle_tol
            ang_resid = _reversal_angle_residual(arg1, states, B, rev_shift)
            pair_angles = ang_resid <= ang_tol + 1e-9
            if not pair_angles:
                warnings.warn(
                    f"the state basin minima of {polymer.name!r} do not obey "
                    f"reflection-with-reversal to within {ang_tol} deg (residual "
                    f"{ang_resid:.2f} deg); leaving the adapted state angles unsymmetrised.",
                    UserWarning,
                    stacklevel=2,
                )
        if pair_angles:
            m = np.array(states.mirror)
            self_mirror = m == np.arange(states.n)
            angs = np.where(self_mirror, angs, np.sign(angs) * 0.5 * (np.abs(angs) + np.abs(angs[m])))
        # keep trans at exactly 180 (a basin minimum at -180 == 180)
        angs = np.where(np.abs(np.abs(angs) - 180.0) < 1e-6, 180.0, angs)
        states = RISStates(states.names, tuple(float(a) for a in angs), states.mirror)
    e3 = None
    guard = None
    if third_order:
        guard = {}
        e3 = np.clip(_fit_third_order(polymer, calc, states, N, base, ref=ref, cap=cap, score=score, diagnostics=guard), -cap, cap)
        n_eval += B * states.n ** 3 * 4
        if mode == "mirror":
            m = np.array(states.mirror)
            e3 = 0.5 * (e3 + e3[:, m][:, :, m][:, :, :, m])
        elif mode == "reversal":
            e3_img, resid3 = _reversal_image3(e3, states.mirror, B, rev_shift)
            if resid3 <= reversal_tol:
                e3 = 0.5 * (e3 + e3_img)
            else:
                warnings.warn(
                    f"reflection-with-reversal does not hold for the third-order terms of "
                    f"{polymer.name!r} to within {reversal_tol} kcal/mol (residual "
                    f"{resid3:.3f}); leaving them unsymmetrised.",
                    UserWarning,
                    stacklevel=2,
                )
    model = RISModel(states, B, e1, e2, e3, name=name or f"{polymer.name}-fit")
    report = FitReport(report_grid, scan1, scan2, n_eval, model, arg1, arg2, angles=angles, discriminator=discriminator,
                       third_order_guard=guard, edge_states=edge_states)
    if relaxer is not None:
        report.relaxed_reference_angles = relaxer.start.copy()
        report.n_relaxation_energies = relaxer.n_energies
        report.relaxation_status = relaxer.status_counts()
    return report


def _fit_third_order(polymer: Polymer, calc: Calculator, states: RISStates, N: int, base: np.ndarray,
                     ref: float = 0.0, cap: float = 50.0, score=None, diagnostics: dict | None = None) -> np.ndarray:
    """Triplet corrections by inclusion-exclusion at the state angles:

        e3[b, a, c, d] = E(a c d) - E(a c T) - E(T c d) + E(T c T)

    i.e. the part of the energy of three consecutive states not captured by the
    pairs (dominated by 1,6-type contacts).  4 * S^3 evaluations per bond type, built
    and evaluated as one batch per bond type.

    An energy above ``ref + cap`` is a steric overlap -- the same reading :func:`fit_ris`
    gives one for the first- and second-order terms -- and *differences* between overlaps
    carry no information, so the inclusion-exclusion is not applied across them:

    * if the triple itself overlaps, ``e3 = +cap``: the triple is impossible and the model
      should say so, whatever the sub-terms do;
    * if one of the *subtracted* terms overlaps but the triple does not, ``e3 = 0``: the
      correction is unresolvable and the pair terms are left to speak for themselves.

    This is not hypothetical bookkeeping.  With rigid geometry a clashing triple is worth
    about 1e6 kcal/mol, so subtracting one such number from another leaves a spurious well
    of 1e6 that the outer clip turns into a flat ``-cap`` -- a large **bonus** for a
    conformation that is in fact impossible.  It is exactly what a potential whose fitted
    gauche basin sits at 60 rather than 80 degrees produces for (G+, G+, G-), where the
    pair term is measured at a relieved *basin minimum* while the triple is evaluated at
    the fixed state angles, so the two disagree by whatever the clash is worth.  Where
    nothing overlaps -- which is every triple of the illustrative potential's PVDF fit --
    the arithmetic is unchanged.

    With ``score`` a relaxed-angle scorer (``fit_ris(angles="relaxed")``) the same rule is
    applied to the relaxed energies against the relaxed all-trans reference.  Relaxing the
    angles changes what counts as an overlap -- a contact the rigid chain could not relieve
    may now cost a few kcal/mol of bend strain instead of a wall -- so fewer triples trip
    the guard and the ones that still do are contacts no bend within the bounds relieves.
    """
    score = _rigid_scorer(polymer, calc) if score is None else score
    B, S = polymer.bonds_per_repeat, states.n
    t_idx = states.index("T") if "T" in states.names else 0
    e3 = np.zeros((B, S, S, S))
    # which triples the guard decided: the triple itself overlapped (capped) or only a
    # subtracted term did (zeroed), and the raw inclusion-exclusion value it would have given
    capped = np.zeros((B, S, S, S), dtype=bool)
    zeroed = np.zeros((B, S, S, S), dtype=bool)
    raw = np.zeros((B, S, S, S))
    for b in range(B):
        j0 = (N // 2 - 1) - ((N // 2 - 1) - b) % B
        rows, keys = [], []
        for a in range(S):
            for c in range(S):
                for d in range(S):
                    for (x, y, z) in ((a, c, d), (a, c, t_idx), (t_idx, c, d), (t_idx, c, t_idx)):
                        dd = base.copy()
                        dd[j0], dd[j0 + 1], dd[j0 + 2] = states.angles[x], states.angles[y], states.angles[z]
                        rows.append(dd)
                    keys.append((a, c, d))
        dihs = np.array(rows)
        E = score(dihs).reshape(-1, 4)
        over = E > ref + cap
        for (a, c, d), row, ov in zip(keys, E, over):
            raw[b, a, c, d] = row[0] - row[1] - row[2] + row[3]
            if ov[0]:
                e3[b, a, c, d] = cap
                capped[b, a, c, d] = True
            elif ov[1] or ov[2] or ov[3]:
                e3[b, a, c, d] = 0.0
                zeroed[b, a, c, d] = True
            else:
                e3[b, a, c, d] = row[0] - row[1] - row[2] + row[3]
    if diagnostics is not None:
        diagnostics.update({"capped": capped, "zeroed": zeroed, "raw": raw})
    return e3


def relax_backbone_angles(polymer: Polymer, calc: Calculator, dihedrals,
                          lo: float = 95.0, hi: float = 135.0, x0=None) -> dict:
    """Minimise ``calc`` over the repeat's backbone angles at fixed torsions.

    One variable per backbone atom of the repeat, broadcast along the chain the way
    :func:`~polyfind.chain.backbone_angle_lookup` does, so a five-monomer PVDF oligomer has
    two: the angle at every CH2 and the angle at every CF2.  The substituents follow the
    changed backbone automatically (:func:`~polyfind.chain.build_chain` places them from
    the actual neighbours), so this really is the chain relaxing rather than a rigid chain
    being re-labelled.

    Only meaningful for a calculator that has an angle term.  Without one nothing resists
    opening an angle and the minimum runs to whichever bound relieves the most contact,
    which is the same thing :func:`polyfind.refine.refine_crystal` says about the cell:
    the bend term plays for the angles the role the Fourier term plays for the torsions.

    Returns the starting and relaxed angles, energies and Lennard-Jones components.  The
    energy is the calculator's own, so with valence terms present it includes the strain
    the relaxation is *buying* the contact relief with, which is the honest accounting.
    """
    from scipy.optimize import minimize

    B = polymer.bonds_per_repeat
    a0 = (np.array([polymer.backbone[k].backbone_angle for k in range(B)], dtype=float)
          if x0 is None else np.asarray(x0, dtype=float))
    dih = np.asarray(dihedrals, dtype=float)

    def energy_at(a):
        return float(calc.energy(build_chain(polymer, dih, bond_angles=np.asarray(a, dtype=float))))

    res = minimize(energy_at, a0, method="Nelder-Mead",
                   bounds=[(lo, hi)] * B, options={"xatol": 1e-3, "fatol": 1e-6, "maxfev": 800})
    best = np.clip(np.asarray(res.x, dtype=float), lo, hi)
    if energy_at(best) > energy_at(a0):  # a failed search must not report a worse minimum
        best = a0
    start = build_chain(polymer, dih)
    end = build_chain(polymer, dih, bond_angles=best)
    return {"angles0": a0, "angles": best,
            "energy0": float(calc.energy(start)), "energy": float(calc.energy(end)),
            "lj0": float(calc.components(start)["lj"]), "lj": float(calc.components(end)["lj"]),
            "n_evaluations": int(res.nfev)}


def oligomer_energy_of_sequence(polymer: Polymer, calc: Calculator, seq, states: RISStates = THREE_STATE, n_periods: int = 6) -> float:
    """Energy per monomer of a finite oligomer built from a periodic sequence (sanity check)."""
    seq = list(seq)
    dih = np.array([states.angles[s] for s in seq * n_periods])
    return calc.energy(build_chain(polymer, dih)) / (len(dih) / polymer.bonds_per_repeat)
