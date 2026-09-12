"""Fit :class:`~polyfind.forcefield.SimpleFF` -- to crystal data, and to reference data.

Why this module exists is recorded in DESIGN.md section 2.5: what is wrong with
``SimpleFF`` is not that it is cheap, it is that it was never fitted, and no
machine-learned potential is or will be a dependency here.  The search around the
potential is now fast enough (a pack-and-refine of a reference chain is seconds) that
a parameter fit can sit on top of it, so that is what this module is.

**There are three fits here, and they are not equals.**  The first, described in the rest of
this docstring, fits five parameters to the experimental crystal data already in the
package; DESIGN.md 5.7 records that it did not generalise and why -- the discriminating
quantities, torsion profiles and conformer energies, were not in the training data.  The
second, in the section after it (:func:`fit_reference`, :class:`ReferenceDesign`,
:func:`acceptance_tests`), fits 27 parameters to first-principles torsion scans and
conformers of VDF oligomers, split by chemistry rather than by frame, and it does
generalise: held-out energy error over ten unseen chemistries falls from 3.32 to 1.76
kcal/mol, and alpha-PVDF finally comes out below beta by an amount the literature agrees
with.  It ships as the preset ``pvdf-dft-fit``; ``examples/fit_dft.py`` reproduces it and
docs/DFT_FIT.md says what it is and is not worth.

The third (:class:`ValenceDesign`, :func:`fit_valence`, the section after that one) adds what
DESIGN.md 5.9 said the *functional form* was missing: harmonic bond and angle terms, full
Cartesian forces in the objective instead of a torsional projection, and an off-site charge
site for fluorine.  56 parameters, same split.  Held-out energy error falls to 1.36 kcal/mol
and held-out force error from 13.9 to 5.2 kcal/(mol A) against a reference whose own RMS is
14.0; fluorine's charge comes out physical and off its bound; the all-trans strain that made
the RIS reference state meaningless for five chemistries is relieved by letting the angles
relax.  Two of the three acceptance tests still fail.  It ships as ``pvdf-dft-valence``;
``examples/fit_valence.py`` reproduces it and docs/VALENCE_FIT.md is blunt about which
numbers did not move.  **No fit here changes any default.**

**What is fitted** (five numbers, :data:`VARIABLES`)::

    eps_r         relative permittivity dividing every Coulomb term
    charge_scale  multiplier on the polymer's point charges
    x_H, x_F      multipliers on the UFF LJ minimum distances of H and F
    V1            the one-fold Fourier torsion coefficient (V2, V3 held at their
                  illustrative values)

**What they are fitted to** (:data:`FITTED_CASES`), all of it data already in the
package plus the polarizations of DESIGN.md 5.5:

* PE orthorhombic, beta-PVDF and alpha-PVDF: the cell edges, the repeat ``c`` and the
  density from :data:`polyfind.pipeline.EXPERIMENTAL_CELLS`;
* beta-PVDF's spontaneous polarization, about 0.13 C/m^2 (but see the caveat at
  ``P_BETA_EXPERIMENT``: that figure is a model estimate, not a measurement);
* that alpha-PVDF's ground state is the *antipolar* cell (|P| = 0), expressed as the
  energy of the antipolar cell relative to the polar one -- the discriminator for the
  electrostatics, since alpha is the one entry of DESIGN.md 5.5 that the illustrative
  potential gets wrong;
* the near-degeneracy of alpha and beta with alpha slightly favoured.

**What is held out** (:data:`HELD_OUT_CASES`): gamma-PVDF, entirely -- its cell, its
density and its polarization are never seen by the objective and exist to answer the
only question that matters here, whether a five-parameter fit to about a dozen numbers
generalises at all.  :func:`fit` reports gamma's error before and after.

**What is deliberately *not* a target: PVDC.**  Its experimental cell is not in this
repository, and the model of it that is here cannot support one: with frozen backbone
angles the all-trans PVDC chain carries ~313 kcal/mol of Lennard-Jones strain
(``polyfind.polymers``, ``tests/test_pvdc.py``), the real polymer does not adopt the
planar zigzag at all, and its crystal is not a two-chain cell with chain 2 at
(1/2, 1/2).  Fitting to a number the model cannot express would flatter the fit and
teach us nothing, so PVDC is left out and said to be left out.

**Identifiability, stated up front because it bounds what any fit here can mean.**
``charge_scale`` and ``eps_r`` enter every *energy* only as ``charge_scale**2 /
eps_r``: they are exactly degenerate in the energy, and only the polarization targets
(which scale linearly with ``charge_scale`` and not at all with ``eps_r``) separate
them.  ``V3`` cannot be fitted from crystal data with ideal-angle chains at all -- its
term vanishes at both 180 and +/-60 degrees -- and ``V2`` only ever appears added to
``V1``, so one torsion coefficient is fitted and the others are left alone.
"""
from __future__ import annotations

import time
from contextlib import contextmanager
from dataclasses import dataclass, field, replace
from typing import Iterable, Sequence

import numpy as np
from scipy.optimize import least_squares, minimize

from .forcefield import (COULOMB, DEFAULT_SCALE14, DEFAULT_TORSION, PRESETS, TRAINSET_ENV,
                         UFF_LJ_EXTRA, Frame, SimpleFF, offset_sites, project_offset_charges,
                         read_frames, valence_topology)
from .pack import CrystalPacker, PackResult, pack, periodic_chain
from .pipeline import EXPERIMENTAL_CELLS
from .polymers import PE, PVDF, UFF_LJ, Polymer, THREE_STATE
from .refine import refine_crystal

T, GP, GM = 0, 1, 2

# Spontaneous polarizations (C/m^2).  Alpha is the antipolar phase, whose polarization
# is zero by symmetry.
#
# The beta value is a fit target and its provenance is now known to be shakier than the
# name suggests: 0.13 is a *rigid-dipole model* estimate, not a measurement (docs/
# REFERENCES.md, entry on the beta polarization).  Measured remanent polarizations for
# the polymer span roughly 0.05 to 0.14 depending on orientation and preparation, while
# perfect-crystal DFT gives 0.176 to 0.188 -- and a perfect crystal is what this package
# actually computes, so the DFT range is arguably the right comparison.  The value is
# left alone because changing it moves fitted parameters, and the fit it feeds did not
# generalise anyway; it is queued for the next re-measurement.
P_BETA_EXPERIMENT = 0.13  # see the caveat above before trusting this as "experimental"


# --------------------------------------------------------------------- parameters
@dataclass(frozen=True)
class FFParameters:
    """One complete parameter set for the potential, in both places it is used.

    ``SimpleFF`` (the oligomer/RIS side) and ``CrystalPacker``/``refine_crystal`` (the
    lattice side) evaluate the *same* Lennard-Jones, Coulomb and torsion terms, so one
    object carries them for both; :meth:`simple_ff` builds the calculator and
    :meth:`applied` makes the lattice side use the same numbers.
    """

    torsion: tuple[float, float, float] = DEFAULT_TORSION
    eps_r: float = 1.0
    charge_scale: float = 1.0
    lj: tuple[tuple[str, float, float], ...] = ()  # (element, x_i, D_i) overrides; hashable
    scale14: float = 0.5
    charge_increments: tuple[tuple[str, str, float], ...] = ()  # bond-charge increments
    # Valence terms and off-site charges (docs/VALENCE_FIT.md).  Empty means "none", which
    # is the rigid, atom-centred potential everything above this fit was measured with.
    bond_terms: tuple[tuple[str, float, float], ...] = ()  # (type, k, r0)
    angle_terms: tuple[tuple[str, float, float], ...] = ()  # (type, k, theta0 deg)
    charge_offsets: tuple[tuple[str, float], ...] = ()  # (element, offset A)
    # Charge flux (element, element, k_angle, k_bond); see FITTED_VALENCE_FLUX and
    # forcefield.FluxTopology.  Like the valence terms it is *not* forwarded by
    # :meth:`applied` -- it is opt-in per packer, ``CrystalPacker(charge_flux=...)``.
    charge_flux: tuple[tuple[str, str, float, float], ...] = ()

    @property
    def lj_dict(self) -> dict[str, tuple[float, float]]:
        return {e: (x, d) for e, x, d in self.lj}

    def simple_ff(self) -> SimpleFF:
        return SimpleFF(torsion=tuple(self.torsion), scale14=self.scale14, eps_r=self.eps_r,
                        charge_scale=self.charge_scale, lj=self.lj_dict or None,
                        charge_increments=self.charge_increments or None,
                        bond_terms=self.bond_terms or None, angle_terms=self.angle_terms or None,
                        charge_offsets=self.charge_offsets or None,
                        charge_flux=self.charge_flux or None)

    def preset_kwargs(self) -> dict:
        """:class:`~polyfind.forcefield.SimpleFF` keyword arguments, for :data:`PRESETS`."""
        kw = {"torsion": tuple(float(v) for v in self.torsion), "scale14": float(self.scale14),
              "eps_r": float(self.eps_r), "charge_scale": float(self.charge_scale)}
        if self.lj:
            kw["lj"] = {e: (float(x), float(d)) for e, x, d in self.lj}
        if self.charge_increments:
            kw["charge_increments"] = tuple((a, b, float(d)) for a, b, d in self.charge_increments)
        for name in ("bond_terms", "angle_terms"):
            v = getattr(self, name)
            if v:
                kw[name] = tuple((t, float(a), float(b)) for t, a, b in v)
        if self.charge_offsets:
            kw["charge_offsets"] = tuple((e, float(d)) for e, d in self.charge_offsets)
        if self.charge_flux:
            kw["charge_flux"] = tuple((a, b, float(ka), float(kb)) for a, b, ka, kb in self.charge_flux)
        return kw

    def describe(self) -> str:
        lj = ", ".join(f"{e} x_i={x:.3f}" for e, x, _ in self.lj) or "UFF unchanged"
        q = ("; ".join(f"d{a}{b}={d:+.3f}" for a, b, d in self.charge_increments)
             or "polymer charges")
        val = (f" [{len(self.bond_terms)} stretch, {len(self.angle_terms)} bend types]"
               if self.has_valence() else " [rigid: no valence terms]")
        off = "".join(f" [{e} charge {d:+.3f} A off its nucleus]" for e, d in self.charge_offsets)
        flx = ("" if not self.charge_flux else " [charge flux: " + "; ".join(
            f"{a}-{b} ka={ka:+.4f} kb={kb:+.4f}" for a, b, ka, kb in self.charge_flux) + "]")
        return (f"torsion=({self.torsion[0]:+.3f}, {self.torsion[1]:+.3f}, {self.torsion[2]:+.3f}) "
                f"eps_r={self.eps_r:.3f} charge_scale={self.charge_scale:.3f} [{lj}] [{q}]{val}{off}{flx}")

    def has_valence(self) -> bool:
        return bool(self.bond_terms) or bool(self.angle_terms)

    @contextmanager
    def applied(self):
        """Make :mod:`polyfind.pack` and :mod:`polyfind.refine` use these parameters.

        Those two modules read the potential from three places that are not arguments:
        the module-level ``UFF_LJ`` table (via ``lj_params``), the charges carried on the
        ``PeriodicChain``, and ``CrystalPacker``'s default torsion triple.  Rather than
        mutate the global table -- which would change the potential for every other
        object alive, including calculators built earlier -- this rebinds two names for
        the duration of the block and puts them back afterwards: ``pack.lj_params`` and
        the ``CrystalPacker`` class that ``pack`` and ``refine`` construct.

        The replacement class forces *this* parameter set on every packer built inside
        the block, so a caller's ``eps_r=1.0`` default cannot quietly override the fit.

        One thing must stay outside: :func:`polyfind.lattice_table.pair_table` keys its
        (on-disk) cache on the chain and on ``(cutoff, alpha, eps_r)``, not on the LJ
        table, so a table built in here could be cached under a key that does not
        describe it.  Screen tables are therefore built by :func:`screen_table` before
        the block is entered and passed in explicitly; the table only chooses starts and
        the polish that follows is exact, so nothing depends on the table matching the
        potential exactly.

        Charges follow the same route.  When :attr:`charge_increments` is set the block's
        charges are recomputed from *its own* bond graph, inferred from the coordinates
        the chain carries, so a fitted charge model reaches the lattice unchanged.  That
        inference is safe here for a reason worth stating: every backbone-backbone bond is
        C-C and a same-element bond carries no increment, so the bonds that cross the
        periodic boundary contribute nothing and a single block gives the same charges as
        an infinite chain would.

        :attr:`bond_terms` and :attr:`angle_terms` are deliberately *not* forwarded, and that
        is not an omission.  ``pack`` builds rigid chains, so inside a packing every bond
        length and every bond angle is fixed and the valence energy is one additive constant
        per chain: it cancels from every lattice-energy difference, from every polarization and
        from every cell parameter, and forwarding it would change nothing but the absolute
        number.  The one place bond angles do move is :func:`polyfind.refine.refine_crystal`,
        which has carried its own harmonic bend term since before this fit and is left exactly
        as it was, so that nothing here can be confused with a change to the refinement stage.

        A caller who *does* want them in the lattice energy -- because the chain is going to be
        deformed, and then they are not a constant -- asks for them by name and by object:
        ``CrystalPacker(chain, valence=SimpleFF.from_preset("pvdf-dft-valence"))``, and
        ``refine_crystal(..., valence=...)``.  That is opt-in, it changes every energy it
        touches, and it is refused by the tabulated screen, which is rigid.  See
        :mod:`polyfind.pack` and ``docs/ELECTROMECHANICS.md`` section 4.2.
        """
        from . import pack as pack_mod
        from . import refine as refine_mod
        from .forcefield import bci_charges, infer_bonds

        ff = self.simple_ff()
        q_scale, torsion, eps_r = self.charge_scale, tuple(self.torsion), self.eps_r
        increments = ff.increments()
        offsets = ff.offset_table()

        class _FittedPacker(pack_mod.CrystalPacker):
            def __init__(self, chain, **kw):
                if increments is not None:
                    bonds = infer_bonds(chain.elements, chain.coords)
                    q = bci_charges(chain.elements, bonds, increments)
                    if offsets:
                        # pack carries one charge per atom and no off-site machinery, so an
                        # off-site model reaches the lattice as its exact dipole-equivalent
                        # projection.  See project_offset_charges: same total charge, same
                        # total dipole, different quadrupole -- and at the 3 A and longer
                        # separations a lattice sum is made of, the dipole is the term that
                        # matters.  docs/VALENCE_FIT.md measures what the difference costs.
                        q = project_offset_charges(chain.coords, q,
                                                   offset_sites(chain.elements, bonds, offsets))
                    chain = replace(chain, charges=q * q_scale)
                elif q_scale != 1.0:
                    chain = replace(chain, charges=np.asarray(chain.charges, dtype=float) * q_scale)
                kw["eps_r"] = eps_r
                kw["torsion"] = torsion
                super().__init__(chain, **kw)

        old_lj, old_packer_pack, old_packer_refine = pack_mod.lj_params, pack_mod.CrystalPacker, refine_mod.CrystalPacker
        pack_mod.lj_params = ff._lj_arrays
        pack_mod.CrystalPacker = _FittedPacker
        refine_mod.CrystalPacker = _FittedPacker
        try:
            yield ff
        finally:
            pack_mod.lj_params = old_lj
            pack_mod.CrystalPacker = old_packer_pack
            refine_mod.CrystalPacker = old_packer_refine


ILLUSTRATIVE = FFParameters()  # exactly today's defaults

# The fitted vector.  Five parameters for ~11 independent observations; see the module
# docstring for why it is not more (V2/V3 and one of eps_r/charge_scale are not
# identifiable from this data) and section "generalisation" of the fit report for what
# even five costs.
VARIABLES = ("eps_r", "charge_scale", "x_H", "x_F", "V1")
BOUNDS = np.array([(1.0, 6.0), (0.40, 1.30), (0.90, 1.15), (0.90, 1.15), (-1.5, 3.0)])
X0 = np.array([1.0, 1.0, 1.0, 1.0, DEFAULT_TORSION[0]])  # = the illustrative potential

UFF_X = {"H": 2.886, "F": 3.364}  # the values the multipliers multiply (polymers.UFF_LJ)
UFF_D = {"H": 0.044, "F": 0.050}


def parameters_from_vector(x) -> FFParameters:
    """Map the fit vector (:data:`VARIABLES`) onto a :class:`FFParameters`."""
    eps_r, q, sH, sF, V1 = (float(v) for v in x)
    return FFParameters(
        torsion=(V1, DEFAULT_TORSION[1], DEFAULT_TORSION[2]),
        eps_r=eps_r,
        charge_scale=q,
        lj=(("H", UFF_X["H"] * sH, UFF_D["H"]), ("F", UFF_X["F"] * sF, UFF_D["F"])),
    )


def vector_from_parameters(p: FFParameters) -> np.ndarray:
    lj = p.lj_dict
    return np.array([p.eps_r, p.charge_scale,
                     lj.get("H", (UFF_X["H"], 0))[0] / UFF_X["H"],
                     lj.get("F", (UFF_X["F"], 0))[0] / UFF_X["F"],
                     p.torsion[0]])


# ------------------------------------------------------------------------- cases
@dataclass(frozen=True)
class CrystalCase:
    """One packed reference structure and the experimental numbers it is compared with."""

    key: str
    polymer: Polymer
    label: str  # key into EXPERIMENTAL_CELLS[polymer.name]
    sequence: tuple[int, ...]
    polarization: float | None = None  # experimental |P| (C/m^2), None if not a target
    antipolar: bool = False  # the phase is antipolar: its ground-state cell must have P = 0

    @property
    def cell(self) -> tuple[float, float, float, float]:
        return EXPERIMENTAL_CELLS[self.polymer.name][self.label]


FITTED_CASES: tuple[CrystalCase, ...] = (
    CrystalCase("pe", PE, "orthorhombic (all-trans)", (T,)),
    CrystalCase("beta", PVDF, "beta (TTTT)", (T, T), polarization=P_BETA_EXPERIMENT),
    CrystalCase("alpha", PVDF, "alpha/delta (TGTG')", (T, GP, T, GM), polarization=0.0, antipolar=True),
)

HELD_OUT_CASES: tuple[CrystalCase, ...] = (
    CrystalCase("gamma", PVDF, "gamma/epsilon (T3GT3G')", (T, T, T, GP, T, T, T, GM)),
)

ALL_CASES = FITTED_CASES + HELD_OUT_CASES


# ------------------------------------------------------------------ predictions
_CHAINS: dict[str, object] = {}
_TABLES: dict[str, object] = {}


def reference_chain(case: CrystalCase):
    """The rigid ideal-angle periodic chain of ``case`` (cached; potential-independent)."""
    if case.key not in _CHAINS:
        _CHAINS[case.key] = periodic_chain(case.polymer, list(case.sequence), THREE_STATE)
    return _CHAINS[case.key]


def screen_table(case: CrystalCase, verbose: bool = False):
    """The screen table for ``case``, built *outside* any :meth:`FFParameters.applied` block.

    It selects the starts the exact-kernel polish then refines, so it is a search
    accelerator and not part of the potential: building it once with the illustrative
    parameters and reusing it for every trial parameter set keeps the inner loop at a
    few seconds and cannot change an energy, only which basins are visited.  That is
    also its one honest weakness, recorded here: a parameter set whose global minimum
    lies in a basin this screen never proposes will be scored at its best *proposed*
    minimum.  The alpha antipolar branch, the one place where that is known to matter,
    is therefore searched explicitly by :func:`antipolar_cell` rather than left to the
    screen.
    """
    from .lattice_table import pair_table
    from .pack import SCREEN_TABLE

    if case.key not in _TABLES:
        _TABLES[case.key] = pair_table(reference_chain(case), cutoff=8.0, alpha=0.2, eps_r=1.0,
                                       verbose=verbose, **SCREEN_TABLE)
    return _TABLES[case.key]


@dataclass
class Prediction:
    """What the potential says about one reference structure."""

    key: str
    a: float
    b: float
    c: float
    density: float
    energy_per_monomer: float  # lattice + torsion + bond-angle strain, per monomer
    polarization: float
    polar_gap: float | None = None  # E(antipolar cell) - E(polar cell), rigid, per monomer
    rigid: PackResult | None = None
    refined: PackResult | None = None
    torsions: np.ndarray | None = None
    seconds: float = 0.0

    @property
    def axes(self) -> tuple[float, float]:
        """(a, b) sorted: the search does not know which axis is which (DESIGN.md 5.2)."""
        return tuple(sorted((self.a, self.b)))


def chain_moment(packer: CrystalPacker) -> np.ndarray:
    """``mu`` of ONE chain in its own frame (e.A): ``sum_i q_i r_i`` over the repeat.

    Independent of the cell, of ``dz`` and of which image the chain is written at, because
    the repeat is neutral.  It is what decides the antipolar subspace
    (:func:`antipolar_offsets`).
    """
    q = np.asarray(packer._q_cell[: packer.n], dtype=float)
    return np.einsum("n,nc->c", q, np.asarray(packer.chain.coords, dtype=float))


def antipolar_offsets(packer: CrystalPacker, tol: float = 1e-9) -> list[tuple[int, float]]:
    r"""``[(flip, dphi)]`` for which ``phi2 = phi1 + dphi`` makes the cell dipole exactly zero.

    Chain 1 contributes ``Rz(phi1) m`` to the cell dipole and chain 2 contributes
    ``Rz(phi2) M^flip m``, where ``m`` is :func:`chain_moment` and ``M = diag(1, -1, -1)``
    is the mirror a flip applies.  So the cell is antipolar for *every* ``(a, b, dz)`` iff
    ``Rz(phi2) M^flip m = -Rz(phi1) m``.  Writing ``m = (m_t cos th, m_t sin th, m_z)``:

    * ``flip = 1`` sends ``th -> -th`` and ``m_z -> -m_z``, so ``dphi = 180 + 2 th``
      satisfies it for **any** chain;
    * ``flip = 0`` leaves ``m_z`` alone, so ``dphi = 180`` satisfies it only when
      ``m_z = 0`` -- true of any planar zigzag, whose moment is perpendicular to its axis.

    **This is a correction.**  :func:`antipolar_cell` uses ``dphi = 0``, which is the above
    only when ``th = 90 deg``.  That holds for PVDF's alpha helix, whose ``m_x`` is zero,
    and fails for its beta zigzag, whose moment lies along its own ``x``: there a flip is a
    *rotation* of the chain (``DESIGN.md`` 5.6 notes the same symmetry for the energy) and
    reverses no dipole at all, so ``antipolar_cell`` returns a cell with the full
    polarization of the polar minimum and an energy degenerate with it.  Every all-trans
    "antipolar gap" obtained through ``antipolar_cell`` is therefore a comparison between
    two polar cells; see ``docs/SCREEN.md``.
    """
    m = chain_moment(packer)
    th = np.degrees(np.arctan2(m[1], m[0]))
    out = [(1, float((180.0 + 2.0 * th) % 360.0))]
    if abs(m[2]) < tol:
        out.append((0, 180.0))
    return out


def antipolar_cell_exact(packer: CrystalPacker, bounds: dict | None = None,
                         screen_packer: CrystalPacker | None = None, target: int = 6000,
                         n_polish: int = 8, maxfev: int = 900, dz_step: float = 0.5,
                         phi_step: float = 30.0, ab_step: float | None = None,
                         pol_tol: float = 1e-9) -> tuple[np.ndarray, float, float]:
    """Lowest cell in the exactly-antipolar subspace of :func:`antipolar_offsets`.

    Searched the way ``pack()`` searches the unconstrained space, which
    :func:`antipolar_cell` does not: a grid over the whole subspace scored with a cheap
    truncated ``screen_packer`` (the packer itself if none is given), then an exact polish
    of the best *distinct* cells.  ``antipolar_cell`` polishes from two starts with no
    screen and with ``(a, b)`` unbounded, so its branch and the polar branch were neither
    bounded alike nor screened alike, and on a near-degenerate pair of phases that decides
    the answer.

    **The axial scan is a length, not a point count, and that is a correction.**  ``dz`` used
    to be sampled at four points across the repeat whatever the repeat was.  On a one-monomer
    all-trans cell that is a 0.64 A step and harmless.  On a **four-monomer gamma-type repeat
    it is one sample per monomer** -- aliased exactly to the interchain registry the scan
    exists to search, so the grid can only ever see one registry phase.  On the 30.8 A repeat
    of the 11:1 VDF/VDCN copolymer it is 7.7 A, three times the 2.56 A monomer period, and
    this function returned -3.7152 kcal/mol per monomer where a finer scan of the identical
    subspace returned -3.7514 (``examples/copolymer_readme.py``).

    ``dz_step`` (A), ``phi_step`` (deg) and ``ab_step`` (A) now set the resolution in absolute
    units, so the grid grows with the repeat instead of thinning out over it.  Two things were
    wrong and both are fixed, and the measurement says which did what: re-run at the polar
    cell's own ``(a, b)``, the copolymer's subspace minimum is -3.7514 at *either* dz
    resolution -- the exact polish recovers the basin from a start the coarse grid scored at
    +1.01 -- while the fine grid's own screen minimum is -2.93 against the coarse grid's
    +1.01.  So on that case the 0.036 this function missed came from its **start selection**,
    which kept only one start per ``(a, b)`` region and so discarded every distinct axial
    registry, and not from the dz step alone.  The selection is now two-pass (below) and the
    dz step no longer aliases.

    **And the change moved no number on any phase of ``docs/SCREEN.md``'s screen**, which is
    recorded here rather than quietly: beta-PVDF's gap is +1.0178 at four dz points and +1.0178
    at six, alpha's +0.0844 at four and at ten, gamma's +0.6649 at four and at nineteen -- and
    gamma is exactly the aliased case, four samples across four monomers.  The exact polish
    recovers the basin from the aliased grid every time it has been checked.  So the fix rests
    on the aliasing being a provable blindness of the *grid* rather than on any number it
    changed, and the claim that survives is the weaker one: **any antipolar energy this
    function produced on a multi-monomer repeat before the change is a lower bound on the
    search rather than a converged minimum, and on every case checked since, the bound was
    tight.**

    The grid is scored with the exact kernel (or with ``screen_packer``'s, which an Ewald
    caller should make its truncated twin, the division of labour ``pack(coulomb="ewald")``
    already draws), so **the screen's cost now grows linearly with the repeat** -- about
    20 ms per cell on a 24-atom chain and 280 on a 74-atom one, 5000 cells on a two-monomer
    repeat and 8000 on a four-monomer one.  ``target`` is the budget the ``(a, b)`` grid is
    sized from, and ``dz_step``, ``phi_step`` and ``ab_step`` (A, deg, A; the last defaults
    to the budget, never coarser than six points) are what a caller trades against it on a
    long repeat, together with a narrowed ``bounds``.  The tabulated screen of
    :mod:`polyfind.lattice_table` is deliberately *not* used here: its ``dz`` axis is itself
    eight points across the repeat, which is the same blindness on a long repeat that this
    function just stopped having.

    ``bounds`` sizes the **screen**; the polish that follows is free in ``(a, b)`` and may
    return a cell slightly outside them -- measured on beta-PVDF, 0.6 % past the upper edge of
    ``b``.  That is deliberate, because the object wanted is the true minimum of the subspace
    rather than the best cell inside a search heuristic, but it does leave the antipolar branch
    searched over a marginally wider domain than ``pack()``'s bounded polish gives the polar
    one.  That asymmetry can only matter if the polar optimum sits *on* a bound, which is
    checkable and for beta-PVDF is not the case (its ``(4.58, 8.55)`` is interior in both axes
    against bounds of ``(4.33, 9.89)``).  Check it before quoting a gap from a chemistry whose
    polar cell is at an edge.

    Returns ``(params, energy_per_monomer, max|P|)``.  The subspace is antipolar
    analytically for every ``(a, b, dz)``, so a non-zero dipole means the construction is
    broken rather than that the caller chose badly: it is checked here against ``pol_tol``
    and raised, *and* returned so a caller can assert on it as well.

    **A caller trap worth naming, because it cost this repository a table.**  ``packer`` must
    be built *inside* the :meth:`FFParameters.applied` block **and resolved through the
    module** (``polyfind.pack.CrystalPacker``, not a name imported earlier), because
    ``applied()`` forces the fitted potential by rebinding that class.  A packer built from an
    earlier-imported name carries the illustrative potential, and then this function's energy
    and the polar branch's ``pack()`` energy are differences of two *different* potentials.
    ``examples/screen_electroactive.py`` did that until ``docs/SCREEN.md``'s re-run found it.
    """
    from .pack import default_bounds

    score = (screen_packer or packer).energy
    c = float(packer.chain.c)
    bd = bounds or default_bounds(packer.chain)
    # widened at the low edge: antiparallel dipoles can pack denser than the polar screen went
    lo_a, lo_b = 0.85 * bd["a"][0], 0.85 * bd["b"][0]
    branches = antipolar_offsets(packer)
    # resolution in angstroms and degrees, not in points: a four-point dz scan is fine on a
    # 2.6 A repeat and blind on a 30.8 A one
    n_dz = max(4, int(np.ceil(c / float(dz_step))))
    n_phi = max(4, int(np.ceil(360.0 / float(phi_step))))
    if ab_step is None:
        n_ab = max(6, int(np.sqrt(target / (len(branches) * n_dz * n_phi))))
        n_a = n_b = n_ab
    else:
        n_a = max(2, int(np.ceil((bd["a"][1] - lo_a) / float(ab_step))) + 1)
        n_b = max(2, int(np.ceil((bd["b"][1] - lo_b) / float(ab_step))) + 1)
    avs = np.linspace(lo_a, bd["a"][1], n_a)
    bvs = np.linspace(lo_b, bd["b"][1], n_b)
    phis = np.linspace(0.0, 360.0, n_phi, endpoint=False)
    dzs = np.linspace(0.0, c, n_dz, endpoint=False)
    dz_tol = max(1.0, 2.0 * float(dz_step))

    def near_ab(g, t) -> bool:
        return abs(g[0] - t[0]) + abs(g[1] - t[1]) < 1.0

    def near_registry(g, t) -> bool:
        return (abs((g[3] - t[3] + 180.0) % 360.0 - 180.0) < 20.0
                and abs((g[5] - t[5] + 0.5 * c) % c - 0.5 * c) < dz_tol)

    rows = []
    for flip, dphi in branches:
        grid = np.array([[a, b, 90.0, p, p + dphi, z, float(flip)]
                         for a in avs for b in bvs for p in phis for z in dzs])
        e = score(grid)
        order = np.argsort(e)
        taken = []
        # Two passes, because the two things a start can be distinct in are not the same
        # thing.  Pass 1 spreads half the starts over different *cells*, as before, since
        # neighbouring (a, b) usually differ only in which (phi, dz) the same basin was
        # entered at.  Pass 2 then fills the rest with distinct *registries* -- a different
        # (phi, dz) at an (a, b) already taken -- which pass 1 throws away and which on a
        # long repeat is the basin that was being missed.
        passes = ((max(1, n_polish // 2), lambda g: any(near_ab(g, t) for t in taken)),
                  (n_polish, lambda g: any(near_ab(g, t) and near_registry(g, t) for t in taken)))
        for limit, rejected in passes:
            for i in order:
                if len(taken) >= limit:
                    break
                g = grid[i]
                if rejected(g):
                    continue
                taken.append(g)
                rows.append((float(e[i]), g, flip, dphi))
    rows.sort(key=lambda r: r[0])
    best, best_e = None, np.inf
    for _, g, flip, dphi in rows[:n_polish]:
        def obj(v, _f=flip, _d=dphi):
            a, b, phi, dz = v
            if a < 2.0 or b < 2.0:
                return 1e6
            return float(packer.energy(np.array([a, b, 90.0, phi, phi + _d, dz, float(_f)])[None])[0])

        r = minimize(obj, [g[0], g[1], g[3], g[5]], method="Nelder-Mead",
                     options={"maxfev": maxfev, "xatol": 1e-3, "fatol": 1e-5})
        if r.fun < best_e:
            best = np.array([r.x[0], r.x[1], 90.0, r.x[2], r.x[2] + dphi, r.x[3], float(flip)])
            best_e = float(r.fun)
    if best is None:
        raise RuntimeError("no cell in the antipolar subspace could be polished")
    pol = float(np.abs(packer.polarization(best[None])).max())
    if not pol < pol_tol:
        raise RuntimeError(
            f"the antipolar subspace returned a cell with max|P| = {pol:.3e} C/m^2 > {pol_tol:g}: "
            f"the construction is wrong, not the cell (branches {branches}, "
            f"chain moment {np.round(chain_moment(packer), 6).tolist()} e.A)")
    return best, best_e / (packer.n_chains * packer.chain.n_monomers), pol


def antipolar_cell(packer: CrystalPacker, start: PackResult, maxfev: int = 250) -> tuple[np.ndarray, float]:
    """Best cell with chain 2 flipped *and* the two setting angles equal.

    **Only antipolar when the chain's transverse moment is perpendicular to its own x
    axis**, which is true of PVDF's alpha helix and false of any planar zigzag; use
    :func:`antipolar_cell_exact`, which derives the subspace from the chain's moment
    instead of assuming it, for anything new.  This function is kept exactly as it was
    because the fit, its acceptance tests and the numbers recorded in ``DESIGN.md`` 5.7,
    5.9 and 5.10 were all measured with it.

    Flipping chain 2 reverses its axial dipole and equal setting angles cancel the
    transverse part *if* that part sits along the chain's own y -- and then the cell dipole
    is exactly zero for any (a, b, dz), which the tests assert for the alpha helix.  It is a
    symmetric subspace of the seven cell parameters, so an unconstrained polish would slide
    out of it towards the polar minimum; tying phi2 to phi1 keeps the search inside it.

    Returns ``(params, energy_per_monomer)``.
    """
    c = packer.chain.c
    best, best_e = None, np.inf
    for phi0, dz0 in ((float(start.phi1), 0.0), (float(start.phi1) + 90.0, c / 2)):
        def obj(v, _p=None):
            a, b, phi, dz = v
            if a < 2.0 or b < 2.0:
                return 1e6
            return float(packer.energy(np.array([a, b, 90.0, phi, phi, dz, 1.0])[None])[0])

        r = minimize(obj, [start.a, start.b, phi0, dz0], method="Nelder-Mead",
                     options={"maxfev": maxfev, "xatol": 1e-3, "fatol": 1e-4})
        if r.fun < best_e:
            best, best_e = r.x, float(r.fun)
    a, b, phi, dz = best
    params = np.array([a, b, 90.0, phi, phi, dz, 1.0])
    return params, best_e / (packer.n_chains * packer.chain.n_monomers)


def predict(case: CrystalCase, params: FFParameters = ILLUSTRATIVE, n_refine: int = 4,
            refine: bool = True, verbose: bool = False, coulomb: str = "dsf", ewald=None,
            antipolar: str = "legacy", always_gap: bool = False) -> Prediction:
    """Pack and refine ``case`` under ``params``: the inner loop of the fit.

    Table screen (prebuilt, see :func:`screen_table`) -> exact-kernel polish -> continuous
    refinement of the torsions, the backbone angles and the cell.  For an antipolar phase
    the symmetric branch of :func:`antipolar_cell` is searched as well and the lower of
    the two rigid branches is the one that gets refined, so the fit cannot be fooled by a
    screen that only ever proposes polar starts.

    **The defaults are the fit's and they do not move.**  ``coulomb="dsf"``, ``ewald=None``,
    ``antipolar="legacy"`` reproduce every number this module has ever recorded, bit for bit.
    The switches exist so that ``docs/SCREEN.md``'s acceptance-test attribution is
    reproducible from the shipped code rather than from a one-off script:

    ``coulomb="ewald"`` (with an optional ``ewald=EwaldSpec(...)``) puts the polish, the
        antipolar branch *and* the refinement on :mod:`polyfind.ewald`; the table screen stays
        truncated either way, which is the division of labour :func:`polyfind.pack.pack`
        already draws.
    ``antipolar="exact"`` measures the antipolar subspace with :func:`antipolar_cell_exact`,
        which derives it from the chain's own moment, instead of :func:`antipolar_cell`, which
        assumes a perpendicular moment and is therefore only correct for a helix.
    ``always_gap=True`` computes the gap for a case not marked ``antipolar`` as well, which is
        what makes beta and gamma usable as controls on the construction.
    """
    from . import pack as pack_mod

    chain = reference_chain(case)
    table = screen_table(case)  # outside the patch, deliberately
    t0 = time.time()
    with params.applied():
        results = pack(chain, table=table, n_refine=n_refine, verbose=verbose,
                       coulomb=coulomb, ewald=ewald)
        rigid = results[0]
        gap = None
        if case.antipolar or always_gap:
            # resolved through the module so that the patched class is the one built
            packer = pack_mod.CrystalPacker(chain, n_chains=2, coulomb=coulomb, ewald=ewald)
            if antipolar == "legacy":
                anti_params, anti_e = antipolar_cell(packer, rigid)
            elif antipolar == "exact":
                screen_pk = None if coulomb == "dsf" else pack_mod.CrystalPacker(chain, n_chains=2)
                anti_params, anti_e, _ = antipolar_cell_exact(packer, screen_packer=screen_pk)
            else:
                raise ValueError(f"unknown antipolar {antipolar!r} (expected 'legacy' or 'exact')")
            gap = anti_e - rigid.energy_per_monomer
            if gap < 0.0:  # the antipolar cell is the ground state: refine that one
                rigid = packer.result(anti_params)
        if not refine:
            e = rigid.energy_per_monomer
            out = Prediction(case.key, rigid.a, rigid.b, rigid.c, rigid.density, e,
                             rigid.polarization_magnitude, gap, rigid, None, rigid.dihedrals)
            out.seconds = time.time() - t0
            return out
        ref = refine_crystal(case.polymer, rigid, eps_r=params.eps_r,
                             coulomb=coulomb, ewald=ewald)
        r = ref.result
    out = Prediction(case.key, r.a, r.b, r.c, r.density,
                     r.energy_per_monomer + ref.angle_energy, r.polarization_magnitude,
                     gap, rigid, r, ref.torsions)
    out.seconds = time.time() - t0
    return out


def predict_all(cases: Iterable[CrystalCase], params: FFParameters = ILLUSTRATIVE, **kw) -> dict[str, Prediction]:
    return {case.key: predict(case, params, **kw) for case in cases}


# -------------------------------------------------------------------- objective
@dataclass(frozen=True)
class Residual:
    name: str
    target: float
    value: float
    tol: float
    weight: float

    @property
    def z(self) -> float:
        return (self.value - self.target) / self.tol

    @property
    def contribution(self) -> float:
        return self.weight * self.z ** 2


# Weights and tolerances, and why they are what they are.
#
# Every residual is (predicted - target) / tol, so tol sets the scale on which a term is
# "one unit wrong" and weight says how many such units the term is worth.  The scales
# come from what the measurement means, not from the spread of the numbers:
#
# * CELL_TOL = 2% of the experimental edge.  The unfitted potential is already within
#   about 7% on the edges and 5% on the densities (DESIGN.md 5.2), and a rigid-geometry,
#   fixed-charge potential has no business claiming better than a couple of percent on a
#   van der Waals contact distance, so 2% is roughly "as good as this model can be".
# * DENSITY_WEIGHT = 0.25, because the density is *not* an independent observation: with
#   gamma fixed at 90 degrees it is exactly the cell mass over a*b*c, all three of which
#   are already targets.  Dropping it entirely would be defensible; keeping it at a
#   quarter weight lets a volume error speak once without letting it speak four times.
# * POLARIZATION: tol 0.02 C/m^2 (about 15% of beta's 0.13) at weight 1.  Beta's value is
#   the one experimental number here that the illustrative potential already reproduces,
#   so this term mostly acts as a constraint that the fit not break it -- which matters,
#   because it is the only thing that pins charge_scale separately from eps_r.
# * POLAR_GAP: weight 2 on a one-sided (hinge) residual, tol 0.1 kcal/mol per monomer.
#   Alpha is the antipolar phase; DESIGN.md 5.5 records that the illustrative potential
#   makes its cell polar (0.078 C/m^2) instead, and 5.4 attributes that to unscreened
#   Coulomb over-rewarding dipole alignment.  This is the discriminator the task calls
#   for, so it carries twice the weight of a cell edge, and it is one-sided: once the
#   antipolar cell wins, winning by more is not better.
# * ENERGY_GAP: weight 2, tol 0.5 kcal/mol per monomer, target -0.05 (alpha slightly
#   below beta).  This is the headline failure of DESIGN.md 5.4 (beta 3.9 kcal/mol per
#   monomer below alpha where experiment has them near-degenerate with alpha favoured).
#   The tolerance is deliberately loose: "near-degenerate" is a statement about tenths of
#   a kcal/mol, and pretending to fit the sign of a 0.05 difference with a fixed-charge
#   potential would be a fiction.
CELL_TOL = 0.02  # relative
DENSITY_WEIGHT = 0.25
POLARIZATION_TOL = 0.02
POLARIZATION_RATIO_TOL = 0.15  # |P_alpha| should be well under a fifth of |P_beta|
POLAR_GAP_TOL = 0.1
POLAR_GAP_WEIGHT = 2.0
ENERGY_GAP_TOL = 0.5
ENERGY_GAP_WEIGHT = 2.0
ENERGY_GAP_TARGET = -0.05  # kcal/mol per monomer, E(alpha) - E(beta)
BIG = 1e4  # objective value for a parameter set that fails to pack or refine


def case_residuals(case: CrystalCase, pred: Prediction) -> list[Residual]:
    """Cell, density and polarization residuals for one structure."""
    ea, eb, ec, erho = case.cell
    ea, eb = sorted((ea, eb))
    pa, pb = pred.axes
    out = [
        Residual(f"{case.key} a", ea, pa, CELL_TOL * ea, 1.0),
        Residual(f"{case.key} b", eb, pb, CELL_TOL * eb, 1.0),
        Residual(f"{case.key} c", ec, pred.c, CELL_TOL * ec, 1.0),
        Residual(f"{case.key} rho", erho, pred.density, CELL_TOL * erho, DENSITY_WEIGHT),
    ]
    if case.polarization is not None and not case.antipolar:
        out.append(Residual(f"{case.key} |P|", case.polarization, pred.polarization, POLARIZATION_TOL, 1.0))
    if case.antipolar and pred.polar_gap is not None:
        # one-sided: E(antipolar) - E(polar) <= 0 is the target, and being further below
        # it is neither better nor worse.
        z = max(0.0, pred.polar_gap)
        out.append(Residual(f"{case.key} antipolar-polar", 0.0, z, POLAR_GAP_TOL, POLAR_GAP_WEIGHT))
    return out


def all_residuals(preds: dict[str, Prediction], cases: Iterable[CrystalCase] = FITTED_CASES) -> list[Residual]:
    """Every residual of the objective, including the alpha/beta energy ordering."""
    cases = list(cases)
    out: list[Residual] = []
    for case in cases:
        if case.key in preds:
            out += case_residuals(case, preds[case.key])
    if "alpha" in preds and "beta" in preds:
        gap = preds["alpha"].energy_per_monomer - preds["beta"].energy_per_monomer
        out.append(Residual("alpha-beta energy", ENERGY_GAP_TARGET, gap, ENERGY_GAP_TOL, ENERGY_GAP_WEIGHT))
        # Alpha's polarization as a *fraction of beta's*, not as an absolute.  Alpha
        # should be near zero and beta near 0.13, and the physical statement is that one
        # phase is antipolar and the other ferroelectric.  Scoring |P_alpha| directly
        # would let the fit score well by shrinking every charge, which is not that
        # statement and would be the wrong mechanism; the ratio is invariant under a
        # uniform charge rescale, so it can only be improved by finding a different cell.
        pb = preds["beta"].polarization
        ratio = preds["alpha"].polarization / pb if pb > 1e-9 else 0.0
        out.append(Residual("alpha |P| / beta |P|", 0.0, ratio, POLARIZATION_RATIO_TOL, 1.0))
    return out


def objective(x, cases: Iterable[CrystalCase] = FITTED_CASES, refine: bool = True,
              verbose: bool = False) -> float:
    """Weighted sum of squared normalised residuals for the parameter vector ``x``."""
    params = parameters_from_vector(np.clip(x, BOUNDS[:, 0], BOUNDS[:, 1]))
    try:
        preds = predict_all(cases, params, refine=refine)
    except Exception as exc:  # a parameter set that cannot be packed is simply bad
        if verbose:
            print(f"  objective: {type(exc).__name__}: {exc}")
        return BIG
    res = all_residuals(preds, cases)
    value = float(sum(r.contribution for r in res))
    if not np.isfinite(value):
        return BIG
    return value


# --------------------------------------------------------------------------- fit
@dataclass
class FitResult:
    parameters: FFParameters
    x: np.ndarray
    objective_before: float
    objective_after: float
    n_evaluations: int
    seconds: float
    before: dict[str, Prediction] = field(default_factory=dict)
    after: dict[str, Prediction] = field(default_factory=dict)
    history: list[tuple[np.ndarray, float]] = field(default_factory=list)

    def table(self) -> str:
        return comparison_table(self.before, self.after)


def fit(n_screen: int = 16, maxfev: int = 70, seed: int = 0, refine: bool = True,
        cases: Iterable[CrystalCase] = FITTED_CASES, verbose: bool = True) -> FitResult:
    """Fit :data:`VARIABLES` to :data:`FITTED_CASES`, holding :data:`HELD_OUT_CASES` out.

    Optimiser: a Latin-hypercube screen of ``n_screen`` parameter sets followed by
    Nelder-Mead (adaptive) from the best of them.  Five parameters, an objective that
    costs seconds and is *not* smooth -- the packing minimum can jump between basins, and
    the refinement is itself an iterative local search, so the objective has small
    discontinuities and no reliable gradient -- rule out anything derivative-based, and
    the budget (order 10^2 evaluations) rules out a population method.  Nelder-Mead needs
    no derivatives, tolerates the noise, and with five parameters its simplex is small
    enough to move; the screen in front of it is what protects against its one real
    weakness, converging into whichever basin it started in.
    """
    cases = list(cases)
    rng = np.random.default_rng(seed)
    t0 = time.time()
    history: list[tuple[np.ndarray, float]] = []

    def f(x) -> float:
        x = np.clip(np.asarray(x, dtype=float), BOUNDS[:, 0], BOUNDS[:, 1])
        v = objective(x, cases, refine=refine)
        history.append((x.copy(), v))
        if verbose:
            print(f"  [{len(history):3d}] f={v:9.3f}  " +
                  "  ".join(f"{k}={val:.3f}" for k, val in zip(VARIABLES, x)), flush=True)
        return v

    if verbose:
        print(f"screening {n_screen} parameter sets ({len(cases)} structures each)")
    lo, hi = BOUNDS[:, 0], BOUNDS[:, 1]
    n = len(VARIABLES)
    # Latin hypercube: one sample per stratum per axis, axes shuffled independently.
    cube = (rng.permuted(np.tile(np.arange(n_screen), (n, 1)), axis=1).T + rng.random((n_screen, n))) / n_screen
    screen = lo + cube * (hi - lo)
    screen[0] = X0  # the illustrative potential is always one of the screened points
    values = [f(x) for x in screen]
    x_start = screen[int(np.argmin(values))]

    if verbose:
        print(f"Nelder-Mead from f={min(values):.3f} at " +
              "  ".join(f"{k}={v:.3f}" for k, v in zip(VARIABLES, x_start)))
    res = minimize(f, x_start, method="Nelder-Mead",
                   options={"maxfev": maxfev, "adaptive": True, "xatol": 1e-3, "fatol": 1e-3})
    x_best = np.clip(res.x, lo, hi)
    if min(v for _, v in history) < res.fun:  # Nelder-Mead can end above a point it visited
        x_best, _ = min(history, key=lambda hv: hv[1])
    best = parameters_from_vector(x_best)
    seconds = time.time() - t0

    all_cases = list(cases) + [c for c in HELD_OUT_CASES if c.key not in {c2.key for c2 in cases}]
    before = predict_all(all_cases, ILLUSTRATIVE, refine=refine)
    after = predict_all(all_cases, best, refine=refine)
    return FitResult(
        parameters=best, x=x_best,
        objective_before=float(sum(r.contribution for r in all_residuals(before, cases))),
        objective_after=float(sum(r.contribution for r in all_residuals(after, cases))),
        n_evaluations=len(history), seconds=seconds, before=before, after=after, history=history,
    )


# ------------------------------------------------------------------- reporting
def comparison_table(before: dict[str, Prediction], after: dict[str, Prediction],
                     cases: Iterable[CrystalCase] = ALL_CASES) -> str:
    """Every target's value and error before and after the fit, held-out ones marked."""
    fitted_keys = {c.key for c in FITTED_CASES}
    rows = ["", f"{'target':<26s} {'experiment':>10s} {'before':>10s} {'err':>8s} {'after':>10s} {'err':>8s}  fitted?"]
    rows.append("-" * 84)

    def line(name, target, b, a, fmt="{:10.3f}", err_fmt="{:+7.1f}%", fitted=True, rel=True):
        if rel and target not in (0.0, None):
            eb, ea = 100 * (b - target) / target, 100 * (a - target) / target
            e_b, e_a = err_fmt.format(eb), err_fmt.format(ea)
        else:
            e_b, e_a = f"{b - target:+7.3f}", f"{a - target:+7.3f}"
        rows.append(f"{name:<26s} {fmt.format(target):>10s} {fmt.format(b):>10s} {e_b:>8s} "
                    f"{fmt.format(a):>10s} {e_a:>8s}  {'yes' if fitted else 'HELD OUT'}")

    for case in cases:
        if case.key not in before or case.key not in after:
            continue
        pb, pa = before[case.key], after[case.key]
        ea, eb, ec, erho = case.cell
        ea, eb = sorted((ea, eb))
        fitted = case.key in fitted_keys
        line(f"{case.key} a (A)", ea, pb.axes[0], pa.axes[0], fitted=fitted)
        line(f"{case.key} b (A)", eb, pb.axes[1], pa.axes[1], fitted=fitted)
        line(f"{case.key} c (A)", ec, pb.c, pa.c, fitted=fitted)
        line(f"{case.key} rho (g/cm3)", erho, pb.density, pa.density, fitted=fitted)
        if case.polarization is not None:
            line(f"{case.key} |P| (C/m2)", case.polarization, pb.polarization, pa.polarization,
                 fitted=fitted, rel=False)
        elif case.key == "gamma":
            rows.append(f"{'gamma |P| (C/m2)':<26s} {'polar':>10s} {pb.polarization:10.3f} {'':>8s} "
                        f"{pa.polarization:10.3f} {'':>8s}  HELD OUT")
        if case.antipolar and pb.polar_gap is not None and pa.polar_gap is not None:
            line(f"{case.key} E(anti)-E(polar)", 0.0, pb.polar_gap, pa.polar_gap, rel=False, fitted=fitted)
    if {"alpha", "beta"} <= set(before) and {"alpha", "beta"} <= set(after):
        gb = before["alpha"].energy_per_monomer - before["beta"].energy_per_monomer
        ga = after["alpha"].energy_per_monomer - after["beta"].energy_per_monomer
        line("E(alpha)-E(beta) /mon", ENERGY_GAP_TARGET, gb, ga, rel=False)
    return "\n".join(rows)


def residual_table(preds: dict[str, Prediction], cases: Iterable[CrystalCase] = FITTED_CASES) -> str:
    res = all_residuals(preds, cases)
    rows = [f"{'residual':<26s} {'target':>9s} {'value':>9s} {'z':>7s} {'w':>5s} {'w z^2':>8s}"]
    for r in res:
        rows.append(f"{r.name:<26s} {r.target:9.3f} {r.value:9.3f} {r.z:7.2f} {r.weight:5.2f} {r.contribution:8.3f}")
    rows.append(f"{'TOTAL':<26s} {'':>9s} {'':>9s} {'':>7s} {'':>5s} {sum(r.contribution for r in res):8.3f}")
    return "\n".join(rows)


def chain_check(params: FFParameters = ILLUSTRATIVE, step: float = 20.0, n_monomers: int = 5,
                max_period: int = 4) -> dict:
    """A generalisation probe the objective never sees: the *isolated-chain* RIS ranking.

    Nothing in the crystal objective constrains the shape of the torsion profile away
    from 180 and +/-60 degrees, so a fit is free to buy the alpha/beta ordering with a
    torsion coefficient that ruins the conformational statistics.  DESIGN.md 5.4 lists
    the illustrative potential's other failure -- its RIS ranking prefers the TG+ 3/1
    helix, which PVDF does not form -- so this refits the RIS model with the fitted
    potential and reports the gauche energies and the ranking, as a check for damage.
    """
    from .enumerate import enumerate_periodic
    from .forcefield import fit_ris

    # third order is what separates TG+TG+ (the 3/1 helix) from TG+TG- (alpha) at all:
    # with pair terms only they are exactly degenerate (DESIGN.md 2.1).
    rep = fit_ris(PVDF, params.simple_ff(), step=step, n_monomers=n_monomers, third_order=True)
    m = rep.model
    gp = THREE_STATE.index("G+")
    cands = enumerate_periodic(PVDF, m, max_period=max_period, k_per_period=20)
    return {
        "gauche_CH2": float(m.first_order[0, gp]),
        "gauche_CF2": float(m.first_order[1, gp]),
        "ranking": [(c.name, float(c.energy_per_monomer), c.known_as) for c in cands[:3]],
        "n_evaluations": rep.n_evaluations,
    }


def register_preset(name: str, params: FFParameters) -> None:
    """Add ``params`` to :data:`polyfind.forcefield.PRESETS` under ``name``."""
    PRESETS[name] = params.preset_kwargs()


# The result of the fit that ``examples/fit_potential.py`` reproduces (86 objective
# evaluations, 1,516 s on 4 CPU cores; weighted objective 100.80 -> 56.47).  It is a
# *named* alternative, never a default: ``SimpleFF()`` is still the illustrative
# potential everywhere in the package, and nothing selects this unless a caller asks for
# it by name.
#
# READ THIS BEFORE USING IT.  The fit improves the fitted targets and does not, on the
# evidence, produce a better potential:
#
# * Almost all of the improvement is in the three terms that have a parameter of their
#   own -- PE's a (-6.8% -> +0.0%, which is x_H), alpha's b (-8.0% -> -3.3%, x_F) and the
#   alpha/beta energy gap (+1.81 -> +0.36 kcal/mol per monomer, which is V1 and nothing
#   else, since at ideal angles the gauche penalty is exactly 0.75 (V1 + V2) per gauche
#   bond).  One knob per target is fitting, not learning.
# * Terms with no knob of their own got *worse*: beta's |P| moved from 0.140 to 0.160
#   C/m^2 against an experimental 0.13, and beta's b and c and alpha's a and c each moved
#   1-3% further out.
# * The discriminator did not move.  Alpha's antipolar cell still costs more than its
#   polar one (+0.197 -> +0.166 kcal/mol per monomer), so the phase is still predicted
#   polar, and no point in the five-parameter box reaches zero.  The fit in fact chose
#   *stronger* net electrostatics (charge_scale^2 / eps_r = 1.42), which is the opposite
#   of the screening DESIGN.md 5.4 blames -- so that diagnosis is not supported by this
#   data through this parameterisation.
# * Held out, gamma-PVDF's cell edges improved (RMS 5.8% -> 3.5%) while its density got
#   much worse (-1.2% -> -9.6%), because the old density was two large edge errors of
#   opposite sign cancelling.  A held-out result that improves on one measure and
#   degrades on another is not evidence of a better potential.
# * :data:`ABLATIONS` separates the two halves and settles which is which.  The two LJ
#   radii alone carry essentially all of the held-out improvement (gamma RMS 3.73%
#   against the full fit's 3.48%) and leave beta's polarization alone (0.128 C/m^2);
#   the other three parameters alone lower the *fitted* objective (100.80 -> 88.96)
#   while making the held-out structure worse than the unfitted potential (gamma RMS
#   5.96% against 5.80%) and beta's polarization worse (0.175).  Two transferable
#   parameters and three fitted ones, on twelve independent observations.
FITTED_PVDF = FFParameters(
    torsion=(0.6441793459233638, -0.05, 2.5),
    eps_r=1.078786016517853,
    charge_scale=1.2391244521272469,
    lj=(("H", 3.0022163314038894, 0.044), ("F", 3.6358797109024574, 0.050)),
)
register_preset("pvdf-crystal-fit", FITTED_PVDF)

# The two halves of the fitted vector, held-out gamma being the thing to watch.  Run by
# ``examples/fit_potential.py --ablate``; the measured numbers are in FITTED_PVDF's note.
ABLATIONS: dict[str, np.ndarray] = {
    "LJ radii only": np.array([1.0, 1.0, vector_from_parameters(FITTED_PVDF)[2],
                               vector_from_parameters(FITTED_PVDF)[3], DEFAULT_TORSION[0]]),
    "everything but the radii": np.array([FITTED_PVDF.eps_r, FITTED_PVDF.charge_scale, 1.0, 1.0,
                                          FITTED_PVDF.torsion[0]]),
}


# =====================================================================================
# Fitting to first-principles reference data (docs/DFT_FIT.md)
# =====================================================================================
#
# Everything above this line fits the potential to *crystal* data, and DESIGN.md 5.7
# records why that did not work: twelve observations, several of them not independent,
# against five parameters, of which two were structurally unidentifiable; and the
# quantities that would discriminate -- torsion profiles and the relative energies of
# conformers -- were not in the training set at all.
#
# This section is the missing ingredient.  It fits the same potential to relaxed torsion
# scans and conformers of VDF-based oligomers labelled with density-functional theory
# (PBE-D3), read from an extended-XYZ file named by ``$POLYFIND_TRAINSET``.  The data is
# not vendored: it is not ours, it is large, and pinning it to one machine's path would
# be worse than asking for the path.
#
# **What is fitted** (:data:`REF_VARIABLES`, 27 numbers):
#
# * three Fourier coefficients for each of four *torsion types*, classified by the
#   substituents on the two central backbone carbons (:data:`REF_TORSION_TYPES`).  Only
#   ``CF2-CH2`` -- every backbone bond of PVDF -- reaches the shipped preset;
# * a multiplier on the UFF Lennard-Jones minimum distance and well depth of C, H, F and
#   Cl (:data:`REF_LJ_ELEMENTS`);
# * six bond-charge increments (:data:`REF_INCREMENTS`), which is the charge model
#   :mod:`polyfind.polymers` already uses written as parameters instead of literals --
#   see :func:`polyfind.forcefield.bci_charges`;
# * the 1-4 nonbonded scale.
#
# **What is NOT fitted, and why.**  ``eps_r`` is held at 1.  It is not an oversight and
# it is not a choice that data could settle: every energy this potential produces depends
# on the charges and the permittivity only through ``q_i q_j / eps_r``, so multiplying
# every charge by ``s`` and ``eps_r`` by ``s^2`` leaves every energy, and therefore every
# force, *identically* unchanged.  Energies and forces cannot separate them however many
# of them there are; only a polarization can, and there is no polarization in this data.
# :func:`permittivity_degeneracy` demonstrates it numerically rather than asserting it.
#
# **Objective**: relative energies within each system plus torsional forces; see
# :class:`ReferenceDesign` for the exact definition and :func:`fit_reference` for the
# force weight and the split.


REFERENCE_ENV = TRAINSET_ENV

# The four torsion types.  A backbone dihedral is typed by the two *central* carbons'
# pendant elements: "FF" is a CF2 carbon, "HH" a CH2, "FH" a CHF.  Those three pairings are
# the VDF/TrFE family and between them cover 74% of the 3,269 backbone torsions in the
# reference set -- CF2-CH2 922 of them over all 31 systems, CF2-CHF 1,402 over 30, CF2-CF2
# 101 over 6.  Everything with a substituted carbon in the middle (nitrile, methyl,
# chlorine, ether...) shares the fourth, 844 torsions over 28 systems, which exists so that
# those systems can contribute to the *shared* parameters without each buying three of
# their own.
REF_TORSION_TYPES: tuple[str, ...] = ("CF2-CH2", "CF2-CHF", "CF2-CF2", "other")
_CLASS_PAIR_TO_TYPE = {("FF", "HH"): 0, ("FF", "FH"): 1, ("FF", "FF"): 2}

# Elements whose Lennard-Jones parameters are fitted, as multipliers on the UFF values.
# The rest (N, O, S, Br, I) stay at UFF: they appear in a handful of systems each and
# fitting them would be fitting those systems rather than the chemistry.
REF_LJ_ELEMENTS: tuple[str, ...] = ("C", "H", "F", "Cl")

# Bond-charge increments that are fitted.  Pairs not listed carry zero, so a nitro group
# or a sulfone is described with charges only on its C-N / C-S bond; those systems are
# nuisance data here, not targets.
REF_INCREMENTS: tuple[tuple[str, str], ...] = (
    ("C", "H"), ("C", "F"), ("C", "Cl"), ("C", "N"), ("C", "O"), ("C", "S"))

REF_VARIABLES: tuple[str, ...] = tuple(
    [f"{v}[{t}]" for t in REF_TORSION_TYPES for v in ("V1", "V2", "V3")]
    + [f"x_{e}" for e in REF_LJ_ELEMENTS]
    + [f"D_{e}" for e in REF_LJ_ELEMENTS]
    + [f"q_{a}-{b}" for a, b in REF_INCREMENTS]
    + ["scale14"])

_NT = len(REF_TORSION_TYPES)
_NL = len(REF_LJ_ELEMENTS)
_NI = len(REF_INCREMENTS)

# The starting point is the illustrative potential: the same Fourier triple on every
# torsion type, UFF radii and well depths untouched, and exactly the charges
# :mod:`polyfind.polymers` gives PVDF and PVDC written as increments.
REF_X0 = np.array([*DEFAULT_TORSION] * _NT + [1.0] * _NL + [1.0] * _NL
                  + [-0.10, 0.20, 0.10, 0.05, 0.05, 0.0] + [DEFAULT_SCALE14])

# Bounds, and why they are where they are.  Two of the three groups are physical
# constraints rather than convenience, and it is worth being explicit about which:
#
# * torsion coefficients: +/- 8 kcal/mol, which no C-C single-bond torsion approaches.
#   Effectively unbounded; the data decides.
# * Lennard-Jones: the minimum distance may move by 10% of its UFF value and the well
#   depth by a factor of two.  This is the same reasoning CELL_TOL uses in the other
#   direction: a rigid-geometry fixed-charge potential has no business claiming better
#   than about ten percent on a van der Waals contact distance, and an *unbounded* fit
#   here does something worse than claim it -- measured, it drives x_C and x_H to the
#   floor and the well depths of C, H and F to a fifth of UFF, i.e. it switches the
#   dispersion off.  It can do that because relative conformer energies within one
#   molecule barely see the overall attraction, so the term is nearly free to be traded
#   against the torsions; a crystal, which is held together by exactly that attraction,
#   sees it immediately.  The bound is where this data stops being able to tell.
# * charge increments: signed by electronegativity, which is not a fact about any target
#   here.  Carbon is more electronegative than hydrogen and less than fluorine, chlorine,
#   nitrogen and oxygen, so H is positive and F, Cl, N, O negative -- an unbounded fit
#   turns fluorine *positive* (measured: q_C-F goes from +0.20 to -0.07), which reverses
#   the CF2 dipole that every polarization result in this package rests on, for a gain in
#   held-out energy error of a few hundredths of a kcal/mol.  Sulfur and carbon are within
#   0.03 of each other on the Pauling scale, so C-S is left free to take either sign.
REF_BOUNDS = np.array(
    [(-8.0, 8.0)] * (3 * _NT)
    + [(0.90, 1.15)] * _NL + [(0.50, 2.00)] * _NL
    + [(-0.60, -0.02), (0.02, 0.60), (0.02, 0.60), (0.02, 0.60), (0.02, 0.60), (-0.30, 0.30)]
    + [(0.0, 1.0)])
# Prior widths: how far a parameter may drift before the (weak) ridge notices.  They are
# scales, not beliefs -- wide enough that the data decides, narrow enough that a
# parameter with almost no data behind it stays near the value it started at.
REF_SIGMA = np.array([4.0] * (3 * _NT) + [0.10] * _NL + [0.60] * _NL + [0.35] * _NI + [0.40])
REF_RIDGE = 0.30  # weight of the whole ridge term, shared over the parameters

# The force weight, in the units :meth:`ReferenceDesign.residuals` defines: how many
# kcal^2/mol^2 of energy error one kcal^2/(mol rad)^2 of torque error is worth.  0.01 is
# what the held-out energy error picks, and it picks it very weakly -- the measured sweep
# is 1.763, 1.758, 1.774, 1.86, 2.05, 2.31 kcal/mol at weights 0, 0.01, 0.03, 0.1, 0.3, 1.
# The forces are worth about a percent, and are in the objective mostly to say so.  Why so
# little, when there are seven torque numbers per frame against one energy, is measured
# rather than guessed and is recorded in docs/DFT_FIT.md: at these geometries the
# reference forces are dominated by a bond-length mismatch this potential structurally
# cannot carry, and what survives the torsional projection is several times larger than
# any single-bond torsion barrier and essentially uncorrelated with what the potential
# says (r = +0.05).
REF_FORCE_WEIGHT = 0.01

_UFF_ALL: dict[str, tuple[float, float]] = {**UFF_LJ_EXTRA, **UFF_LJ}


def carbon_class(frame: Frame, i: int, backbone: set) -> str:
    """The pendant elements of backbone carbon ``i``, sorted and joined: ``"HH"``,
    ``"FF"``, ``"FH"``, ``"ClCl"``, ``"CH"``..."""
    return "".join(sorted(frame.elements[j] for j in frame.adjacency[i] if j not in backbone))


def torsion_types(frame: Frame) -> np.ndarray:
    """Index into :data:`REF_TORSION_TYPES` for each backbone torsion of ``frame``."""
    bb = set(int(x) for x in frame.backbone)
    out = np.empty(frame.n_dihedrals, dtype=int)
    for k, (_, b, c, _) in enumerate(np.asarray(frame.torsions, dtype=int)):
        pair = tuple(sorted((carbon_class(frame, int(b), bb), carbon_class(frame, int(c), bb))))
        out[k] = _CLASS_PAIR_TO_TYPE.get(pair, _NT - 1)
    return out


def load_reference(path: str | None = None) -> list[Frame]:
    """The reference set, with the systems this potential cannot represent removed.

    A whole *system* is dropped -- not a frame -- if any of its frames has a backbone
    that is unsaturated (a C-C bond under 1.42 A) or passes through a ring.  Both are
    outside what a single three-term Fourier torsion on a freely rotating single bond can
    express: a C=C torsion is a 60-plus kcal/mol barrier and a backbone inside a
    three-membered ring is not a rotor at all.  Dropping frames rather than systems would
    leave a system's remaining frames sampling only part of its own coordinate, which is
    worse than not having it.  The measured rejects are ``ene``, ``cnene``, ``dicnene``
    and ``fa_tetra`` (backbone C=C, energy spreads of 72 to 144 kcal/mol against 4 to 17
    for everything kept) and ``epo`` and ``cnepo`` (backbone epoxide).
    """
    frames = read_frames(path)
    bad = {f.system for f in frames if not (f.backbone_is_saturated() and f.backbone_is_acyclic())}
    return [f for f in frames if f.system not in bad]


def split_systems(frames: Iterable[Frame], anchor: str = "pvdf", stride: int = 3,
                  offset: int = 1) -> tuple[list[str], list[str]]:
    """Split the systems -- never the frames -- into train and test.

    Frames within one system are strongly correlated: a torsion scan is one molecule at
    twelve angles, and its energies share every error the geometry makes.  A random frame
    split would put nine of them in train and three in test and report a small test error
    that means nothing, which is exactly how DESIGN.md 5.7's fit was flattered.  So whole
    chemistries are held out: every ``stride``-th system in alphabetical order, which is
    reproducible, has nothing to do with how well any system fits, and spreads the test
    set over the substituent types rather than clustering it.  ``anchor`` (PVDF) is forced
    into the training set: it is the chemistry every acceptance test is about.
    """
    systems = sorted({f.system for f in frames})
    test = [s for i, s in enumerate(systems) if i % stride == offset and s != anchor]
    return [s for s in systems if s not in test], test


class ReferenceDesign:
    """Every geometry-dependent quantity of a set of frames, precomputed once.

    The frames never move, so distances, dihedral angles, exclusion masks and the
    torsional projections ``dr_ij / d phi_k`` are constants; only the parameters change.
    Precomputing them turns one objective evaluation into a handful of array operations
    over flat arrays -- about 20 ms for 300 frames, against roughly a minute if each
    evaluation rebuilt topologies and differentiated numerically -- which is what makes a
    27-parameter least-squares fit with numerical derivatives affordable at all.

    :meth:`evaluate` reproduces :meth:`polyfind.forcefield.SimpleFF.energy_frame` and
    :meth:`~polyfind.forcefield.SimpleFF.frame_torques` exactly; ``tests/test_fitting.py``
    asserts that against the reference implementation rather than trusting it.
    """

    def __init__(self, frames: Sequence[Frame]):
        frames = list(frames)
        if not frames:
            raise ValueError("no frames")
        self.frames = frames
        self.systems = sorted({f.system for f in frames})
        sidx = {s: i for i, s in enumerate(self.systems)}
        el_index = {e: i for i, e in enumerate(sorted(_UFF_ALL))}
        self.elements = sorted(_UFF_ALL)
        self.el_index = el_index

        p_i, p_j, p_s14, p_r, p_frame = [], [], [], [], []
        d_frame, d_phi, d_type, d_tau = [], [], [], []
        s_pair, s_dih, s_val = [], [], []
        a_el, a_M = [], []
        energies, f_sys = [], []
        n_atoms = n_pairs = n_dih = 0
        for fi, f in enumerate(frames):
            base = n_atoms
            a_el.extend(el_index[e] for e in f.elements)
            a_M.append(_increment_matrix(f))
            dist = f.graph_distances()
            iu, ju = np.triu_indices(f.n_atoms, 1)
            keep = dist[iu, ju] >= 3
            iu, ju = iu[keep], ju[keep]
            r = np.linalg.norm(f.coords[iu] - f.coords[ju], axis=1)
            p_i.append(iu + base)
            p_j.append(ju + base)
            p_s14.append(dist[iu, ju] == 3)
            p_r.append(r)
            p_frame.append(np.full(len(iu), fi))
            rhat = (f.coords[iu] - f.coords[ju]) / r[:, None]
            rot = f.rotors()
            tau = f.reference_torques() if f.forces is not None else np.zeros(f.n_dihedrals)
            phi = np.radians(f.angles())
            types = torsion_types(f)
            for k in range(f.n_dihedrals):
                d_frame.append(fi)
                d_phi.append(phi[k])
                d_tau.append(tau[k])
                d_type.append(int(types[k]))
                S = (rhat * (rot[k][iu] - rot[k][ju])).sum(1)
                nz = np.flatnonzero(np.abs(S) > 1e-12)
                s_pair.append(nz + n_pairs)
                s_dih.append(np.full(len(nz), n_dih))
                s_val.append(S[nz])
                n_dih += 1
            n_pairs += len(iu)
            n_atoms += f.n_atoms
            energies.append(f.energy)
            f_sys.append(sidx[f.system])

        self.pair_i = np.concatenate(p_i)
        self.pair_j = np.concatenate(p_j)
        self.pair_s14 = np.concatenate(p_s14)
        self.inv_r = 1.0 / np.concatenate(p_r)
        self.pair_frame = np.concatenate(p_frame).astype(int)
        self.dih_frame = np.array(d_frame, dtype=int)
        self.dih_phi = np.array(d_phi)
        self.dih_type = np.array(d_type, dtype=int)
        self.torque_ref = np.array(d_tau)
        self.s_pair = np.concatenate(s_pair) if s_pair else np.empty(0, dtype=int)
        self.s_dih = np.concatenate(s_dih) if s_dih else np.empty(0, dtype=int)
        self.s_val = np.concatenate(s_val) if s_val else np.empty(0)
        self.atom_el = np.array(a_el, dtype=int)
        self.increment_matrix = np.vstack(a_M)
        self.energy_ref = np.array(energies)
        self.frame_system = np.array(f_sys, dtype=int)
        self.n_frames = len(frames)
        self.n_torsions = n_dih
        self.system_count = np.bincount(self.frame_system, minlength=len(self.systems))
        self.energy_ref_centred = self.centre(self.energy_ref)

    def centre(self, x: np.ndarray) -> np.ndarray:
        """Subtract each system's own mean.

        The *mean* and not the lowest frame.  Absolute energies are not comparable across
        chemistries, so some per-system offset has to go; the mean is the projection that
        removes it optimally in least squares, and it treats every frame alike.  Referencing
        to the lowest frame instead would make every residual of a system depend on that one
        frame's error, correlating them all and letting a single bad geometry shift a whole
        system's target.  Nothing is lost: the two differ by a constant per system, which is
        precisely what the centring removes.
        """
        m = np.bincount(self.frame_system, weights=x, minlength=len(self.systems)) / self.system_count
        return x - m[self.frame_system]

    def evaluate(self, p: dict) -> tuple[np.ndarray, np.ndarray]:
        """``(energies, torques)`` for the unpacked parameter dict of
        :func:`reference_unpack`: one energy per frame (kcal/mol) and one torque per
        backbone torsion (kcal/(mol rad))."""
        x = p["lj_x"][self.atom_el]
        d = p["lj_d"][self.atom_el]
        q = self.increment_matrix @ p["increments"]
        scale = np.where(self.pair_s14, p["scale14"], 1.0)
        lj_x = np.sqrt(x[self.pair_i] * x[self.pair_j])
        lj_d = np.sqrt(d[self.pair_i] * d[self.pair_j]) * scale
        qq = q[self.pair_i] * q[self.pair_j] * scale * COULOMB
        s6 = (lj_x * self.inv_r) ** 6
        e_pair = lj_d * (s6 * s6 - 2.0 * s6) + qq * self.inv_r
        # d(pair energy)/dr, for the torsional projection
        de_dr = (12.0 * lj_d * self.inv_r) * (s6 - s6 * s6) - qq * self.inv_r ** 2
        V = p["torsion"][self.dih_type]
        phi = self.dih_phi
        e_tors = 0.5 * (V[:, 0] * (1 + np.cos(phi)) + V[:, 1] * (1 - np.cos(2 * phi))
                        + V[:, 2] * (1 + np.cos(3 * phi)))
        dt_dphi = 0.5 * (-V[:, 0] * np.sin(phi) + 2 * V[:, 1] * np.sin(2 * phi)
                         - 3 * V[:, 2] * np.sin(3 * phi))
        energy = (np.bincount(self.pair_frame, weights=e_pair, minlength=self.n_frames)
                  + np.bincount(self.dih_frame, weights=e_tors, minlength=self.n_frames))
        torque = -(np.bincount(self.s_dih, weights=de_dr[self.s_pair] * self.s_val,
                               minlength=self.n_torsions) + dt_dphi)
        return energy, torque

    def residuals(self, x, force_weight: float, ridge: float = REF_RIDGE) -> np.ndarray:
        """The least-squares residual vector: relative energies, torques, ridge.

        Energies are centred per system and divided by ``sqrt(n_frames)``, so the energy
        block's sum of squares is the mean squared error in kcal^2/mol^2 and its square
        root is a number in kcal/mol that means what it says.  The torque block is scaled
        the same way and multiplied by ``sqrt(force_weight)``, so ``force_weight`` is
        literally "how many kcal^2/mol^2 of energy error one kcal^2/(mol rad)^2 of torque
        error is worth".
        """
        p = reference_unpack(x)
        energy, torque = self.evaluate(p)
        blocks = [(self.centre(energy) - self.energy_ref_centred) / np.sqrt(self.n_frames)]
        if force_weight > 0:
            blocks.append(np.sqrt(force_weight) * (torque - self.torque_ref) / np.sqrt(self.n_torsions))
        if ridge > 0:
            blocks.append(np.sqrt(ridge / len(x)) * (np.asarray(x, dtype=float) - REF_X0) / REF_SIGMA)
        return np.concatenate(blocks)

    def errors(self, x) -> dict:
        """Root-mean-square energy error (kcal/mol) and torque error (kcal/(mol rad)),
        overall and per system."""
        p = reference_unpack(x)
        energy, torque = self.evaluate(p)
        de = self.centre(energy) - self.energy_ref_centred
        dt = torque - self.torque_ref
        per = {}
        for i, s in enumerate(self.systems):
            m = self.frame_system == i
            per[s] = float(np.sqrt((de[m] ** 2).mean()))
        return {"energy_rms": float(np.sqrt((de ** 2).mean())),
                "torque_rms": float(np.sqrt((dt ** 2).mean())),
                "energy_max": float(np.abs(de).max()),
                "per_system": per}


def _increment_matrix(frame: Frame) -> np.ndarray:
    """``(n_atoms, len(REF_INCREMENTS))``: charges are ``M @ increments``, exactly."""
    M = np.zeros((frame.n_atoms, _NI))
    lookup = {pair: t for t, pair in enumerate(REF_INCREMENTS)}
    for i, j in frame.bonds:
        ei, ej = frame.elements[i], frame.elements[j]
        if (ei, ej) in lookup:
            t = lookup[(ei, ej)]
            M[i, t] += 1.0
            M[j, t] -= 1.0
        elif (ej, ei) in lookup:
            t = lookup[(ej, ei)]
            M[j, t] += 1.0
            M[i, t] -= 1.0
    return M


def reference_unpack(x) -> dict:
    """The fit vector as the arrays :meth:`ReferenceDesign.evaluate` wants."""
    x = np.asarray(x, dtype=float)
    torsion = x[:3 * _NT].reshape(_NT, 3)
    off = 3 * _NT
    names = sorted(_UFF_ALL)
    lj_x = np.array([_UFF_ALL[e][0] for e in names])
    lj_d = np.array([_UFF_ALL[e][1] for e in names])
    for k, e in enumerate(REF_LJ_ELEMENTS):
        i = names.index(e)
        lj_x[i] = _UFF_ALL[e][0] * x[off + k]
        lj_d[i] = _UFF_ALL[e][1] * x[off + _NL + k]
    off += 2 * _NL
    return {"torsion": torsion, "lj_x": lj_x, "lj_d": lj_d,
            "increments": x[off:off + _NI], "scale14": float(x[-1])}


def reference_ff_parameters(x, eps_r: float = 1.0) -> FFParameters:
    """The fit vector as an :class:`FFParameters` for PVDF and its immediate relatives.

    Every backbone bond of PVDF is a ``CF2-CH2`` torsion, so that one type's triple *is*
    the potential's ``torsion``; the other three types describe bonds PVDF does not have
    and are carried only by :data:`REF_VARIABLES`.  The Lennard-Jones multipliers and the
    charge increments transfer as they stand.
    """
    p = reference_unpack(x)
    names = sorted(_UFF_ALL)
    lj = tuple((e, float(p["lj_x"][names.index(e)]), float(p["lj_d"][names.index(e)]))
               for e in REF_LJ_ELEMENTS)
    inc = tuple((a, b, float(v)) for (a, b), v in zip(REF_INCREMENTS, p["increments"]))
    return FFParameters(torsion=tuple(float(v) for v in p["torsion"][0]), eps_r=eps_r,
                        charge_scale=1.0, lj=lj, scale14=float(p["scale14"]),
                        charge_increments=inc)


def permittivity_degeneracy(design: "ReferenceDesign", x=None, s: float = 1.7) -> dict:
    """Measure, rather than assert, that ``eps_r`` and the charge scale are degenerate.

    Multiply every charge by ``s`` and the permittivity by ``s**2`` and nothing an energy
    or a force can see changes.  ``eps_r`` is not a variable of :meth:`residuals` -- it
    is fixed at 1 and the increments carry the whole electrostatic scale -- so the check
    is done directly on the Coulomb term: it reports the largest change in any frame
    energy and any torque, which should be at rounding level.
    """
    x = REF_X0 if x is None else np.asarray(x, dtype=float)
    p = reference_unpack(x)
    e0, t0 = design.evaluate(p)
    q = p.copy()
    q["increments"] = p["increments"] * s
    e1, t1 = design.evaluate(q)
    # undo the s**2 the Coulomb term picked up, which is exactly what eps_r = s**2 does
    coul0 = e0 - _no_charge(design, p)[0]
    coul1 = e1 - _no_charge(design, q)[0]
    tc0 = t0 - _no_charge(design, p)[1]
    tc1 = t1 - _no_charge(design, q)[1]
    return {"scale": s,
            "max_energy_change": float(np.abs(coul1 / s ** 2 - coul0).max()),
            "max_torque_change": float(np.abs(tc1 / s ** 2 - tc0).max())}


def _no_charge(design: "ReferenceDesign", p: dict):
    q = dict(p)
    q["increments"] = np.zeros_like(p["increments"])
    return design.evaluate(q)


@dataclass
class ReferenceFit:
    """One fit of :data:`REF_VARIABLES` to reference data, with its held-out score."""

    x: np.ndarray
    force_weight: float
    train_systems: list[str]
    test_systems: list[str]
    train: dict
    test: dict
    n_evaluations: int
    seconds: float

    @property
    def parameters(self) -> FFParameters:
        return reference_ff_parameters(self.x)

    def table(self) -> str:
        rows = [f"{'system':<12s} {'frames':>6s} {'E rms (kcal/mol)':>18s}  set"]
        rows.append("-" * 48)
        for s in self.train_systems:
            if s in self.train["per_system"]:
                rows.append(f"{s:<12s} {'':>6s} {self.train['per_system'][s]:18.3f}  train")
        for s in self.test_systems:
            if s in self.test["per_system"]:
                rows.append(f"{s:<12s} {'':>6s} {self.test['per_system'][s]:18.3f}  HELD OUT")
        rows.append("-" * 48)
        rows.append(f"{'train':<12s} {'':>6s} {self.train['energy_rms']:18.3f}")
        rows.append(f"{'held out':<12s} {'':>6s} {self.test['energy_rms']:18.3f}")
        return "\n".join(rows)


def fit_reference(frames: Sequence[Frame] | None = None, force_weight: float = REF_FORCE_WEIGHT,
                  path: str | None = None, n_starts: int = 3, seed: int = 0,
                  free: Sequence[int] | None = None, x0=None, ridge: float = REF_RIDGE,
                  verbose: bool = True) -> ReferenceFit:
    """Fit :data:`REF_VARIABLES` to the reference set, holding whole systems out.

    Optimiser: bounded Levenberg-Marquardt (``scipy.optimize.least_squares`` with a
    trust-region reflective step) from ``n_starts + 1`` starting points -- the
    illustrative potential plus perturbations of it.  Unlike the crystal fit of
    DESIGN.md 5.7, this objective *is* smooth and cheap (no packing, no basin hopping,
    about 20 ms per evaluation), so a derivative-based least-squares method is the right
    tool and a few hundred evaluations are nothing.

    ``free`` restricts the fit to a subset of the parameter indices, the rest held at
    ``x0``; :data:`REF_ABLATIONS` uses it.
    """
    if frames is None:
        frames = load_reference(path)
    frames = list(frames)
    train_sys, test_sys = split_systems(frames)
    design = ReferenceDesign([f for f in frames if f.system in train_sys])
    held = ReferenceDesign([f for f in frames if f.system in test_sys])
    x0 = REF_X0 if x0 is None else np.asarray(x0, dtype=float)
    idx = np.arange(len(REF_X0)) if free is None else np.asarray(free, dtype=int)
    lo, hi = REF_BOUNDS[idx, 0], REF_BOUNDS[idx, 1]

    def expand(v):
        full = x0.copy()
        full[idx] = v
        return full

    calls = {"n": 0}

    def f(v):
        calls["n"] += 1
        return design.residuals(expand(v), force_weight, ridge)

    rng = np.random.default_rng(seed)
    starts = [x0[idx]] + [np.clip(x0[idx] + rng.normal(0, REF_SIGMA[idx]), lo, hi)
                          for _ in range(n_starts)]
    t0 = time.time()
    best, best_cost = None, np.inf
    for k, s in enumerate(starts):
        res = least_squares(f, s, bounds=(lo, hi), x_scale=REF_SIGMA[idx], max_nfev=6000)
        if verbose:
            print(f"  start {k}: cost {res.cost:.6f} in {res.nfev} evaluations", flush=True)
        if res.cost < best_cost:
            best, best_cost = res.x.copy(), float(res.cost)
    x = expand(best)
    return ReferenceFit(x=x, force_weight=force_weight, train_systems=train_sys,
                        test_systems=test_sys, train=design.errors(x), test=held.errors(x),
                        n_evaluations=calls["n"], seconds=time.time() - t0)


def reference_table(x, design: "ReferenceDesign", held: "ReferenceDesign") -> str:
    """Per-system energy error, train and held out, in kcal/mol."""
    tr, te = design.errors(x), held.errors(x)
    rows = ["", f"{'system':<12s} {'frames':>7s} {'E rms':>9s} {'E max':>9s}  set", "-" * 50]
    for d, label, err in ((design, "train", tr), (held, "HELD OUT", te)):
        for i, s in enumerate(d.systems):
            n = int((d.frame_system == i).sum())
            rows.append(f"{s:<12s} {n:7d} {err['per_system'][s]:9.3f} {'':>9s}  {label}")
    rows.append("-" * 50)
    rows.append(f"{'train':<12s} {design.n_frames:7d} {tr['energy_rms']:9.3f} {tr['energy_max']:9.3f}")
    rows.append(f"{'held out':<12s} {held.n_frames:7d} {te['energy_rms']:9.3f} {te['energy_max']:9.3f}")
    return "\n".join(rows)


# --------------------------------------------------------------- acceptance tests
#
# These matter more than the objective value.  Each is a documented failure of the
# illustrative potential that a *better* potential should fix, and none of them is in
# the fit's objective -- the fit never sees a crystal, a polymorph or an RIS ranking.
#
# 1. alpha below beta by 2.6-6.5 kJ/mol per monomer.  Four independent studies across
#    five exchange-correlation functionals agree on the sign and roughly on the
#    magnitude (docs/REFERENCES.md section 5); the illustrative potential has the sign
#    backwards and the magnitude three to six times too large.
# 2. The isolated-chain RIS ranking does not put the TG+ 3/1 helix first.  PVDF does not
#    form it (DESIGN.md 5.4).
# 3. alpha's packed cell is antipolar: the antipolar cell must not cost more than the
#    polar one.
ALPHA_BETA_RANGE_KJ = (-6.5, -2.6)  # E(alpha) - E(beta), kJ/mol per monomer
KJ_PER_KCAL = 4.184


def acceptance_tests(params: FFParameters = ILLUSTRATIVE, step: float = 20.0,
                     verbose: bool = False, coulomb: str = "dsf", ewald=None,
                     antipolar: str = "legacy") -> dict:
    """Run the three acceptance tests and report each as a pass or a fail with numbers.

    ``coulomb``, ``ewald`` and ``antipolar`` go straight to :func:`predict` and default to the
    fit's own choices, so calling this with no arguments reproduces every recorded score.
    ``docs/SCREEN.md`` reports the 2x2 over ``coulomb`` and ``antipolar``, which is how test 3's
    failure is attributed: neither switch moves it, and its margin is a third of the potential's
    own resolution.  Beta's gap comes back too (``beta_polar_gap``), as the control that shows
    what the two constructions differ on -- ``+0.0000`` with a polarization of 0.115 C/m^2 for
    the legacy one, which is beta's own polar minimum, against ``+1.018`` for the exact one.
    """
    preds = predict_all([c for c in ALL_CASES if c.key in ("alpha", "beta")], params,
                        verbose=verbose, coulomb=coulomb, ewald=ewald, antipolar=antipolar,
                        always_gap=True)
    gap = preds["alpha"].energy_per_monomer - preds["beta"].energy_per_monomer
    gap_kj = gap * KJ_PER_KCAL
    chk = chain_check(params, step=step)
    top = chk["ranking"][0][0]
    polar_gap = preds["alpha"].polar_gap
    return {
        "alpha_beta_kj": gap_kj,
        "alpha_beta_pass": bool(ALPHA_BETA_RANGE_KJ[0] <= gap_kj <= ALPHA_BETA_RANGE_KJ[1]),
        "ris_top": top,
        "ris_ranking": chk["ranking"],
        "ris_pass": top != "TG+",
        "alpha_polar_gap": polar_gap,
        "alpha_antipolar_pass": bool(polar_gap is not None and polar_gap <= 0.0),
        "alpha_polarization": preds["alpha"].polarization,
        "beta_polarization": preds["beta"].polarization,
        # beta is the control, not a test: it is the ferroelectric phase, so its gap must come
        # out positive, and the two antipolar constructions differ on it by 1.02 kcal/mol per
        # monomer (docs/SCREEN.md).  Reported rather than asserted, since no acceptance test
        # was ever written for it and adding one now would be changing the score.
        "beta_polar_gap": preds["beta"].polar_gap,
        "coulomb": coulomb if coulomb == "dsf" else f"ewald[{getattr(ewald, 'boundary', 'tinfoil')}]",
        "antipolar_construction": antipolar,
        "predictions": preds,
    }


def acceptance_table(result: dict) -> str:
    ok = lambda b: "PASS" if b else "FAIL"
    rank = ", ".join(f"{n} ({e:+.2f})" for n, e, _ in result["ris_ranking"])
    return "\n".join([
        f"  1. E(alpha) - E(beta)      {result['alpha_beta_kj']:+8.2f} kJ/mol per monomer "
        f"(want {ALPHA_BETA_RANGE_KJ[0]:+.1f} to {ALPHA_BETA_RANGE_KJ[1]:+.1f})   {ok(result['alpha_beta_pass'])}",
        f"  2. isolated-chain RIS      top = {result['ris_top']:<8s} "
        f"(want anything but TG+)                     {ok(result['ris_pass'])}",
        f"     ranking: {rank}",
        f"  3. alpha E(anti)-E(polar)  {result['alpha_polar_gap']:+8.3f} kcal/mol per monomer "
        f"(want <= 0)                {ok(result['alpha_antipolar_pass'])}",
        f"     |P| alpha {result['alpha_polarization']:.4f}, beta {result['beta_polarization']:.4f} C/m^2"
        + ("" if result.get("beta_polar_gap") is None else
           f";  beta E(anti)-E(polar) {result['beta_polar_gap']:+.3f} (control, wants > 0)"),
        f"     sum {result.get('coulomb', 'dsf')}, antipolar construction "
        f"{result.get('antipolar_construction', 'legacy')!r}",
    ])


# ------------------------------------------------------------------ the result
#
# The fit ``examples/fit_dft.py`` reproduces: 27 parameters, force weight 0.01, four
# starts, 21 training chemistries (312 frames) and 10 held out (155), about 40 s.
# Numbers, and what they are worth, in docs/DFT_FIT.md.  Like ``pvdf-crystal-fit`` this is
# a *named* alternative and never a default: ``SimpleFF()`` is still the illustrative
# potential everywhere, and nothing selects this unless a caller asks for it by name.
#
# READ THIS BEFORE USING IT.
#
# * It **generalises**, which the crystal fit did not.  Held-out energy error over ten
#   chemistries the objective never saw falls from 3.32 to 1.76 kcal/mol (train 3.16 ->
#   1.50), and no group of parameters gets near that alone -- torsions only 3.12,
#   Lennard-Jones only 2.53, charges only 3.18, the 1-4 scale only 3.36, which is the
#   opposite of DESIGN.md 5.7's finding that three of five parameters fit noise.
# * It fixes **one of three** acceptance tests, none of which is in the objective.
#   alpha-PVDF comes out 4.54 kJ/mol per monomer below beta, inside the 2.6-6.5 range four
#   independent studies agree on and a change of sign from the +7.38 the illustrative
#   potential gives.  That is DESIGN.md 5.4's headline failure, fixed.
# * It does **not** fix the other two.  The isolated-chain RIS ranking still puts the TG+
#   3/1 helix first, though its margin over the next candidate falls from 1.64 to 0.23
#   kcal/mol per monomer; and alpha's antipolar cell still costs more than its polar one
#   (+0.091 kcal/mol per monomer against +0.197 unfitted), so the phase is still predicted
#   polar.  Both move the right way and neither arrives.
# * Three parameters sit on their bounds -- ``x_C``, ``x_H`` and ``q_C-F`` -- which is the
#   data asking for something unphysical and being stopped.  ``q_C-F`` at 0.02 e means the
#   fit would rather have almost no C-F electrostatics at all; beta's polarization drops
#   from 0.140 to 0.114 C/m^2 as a result, away from the DFT value of 0.176-0.188.  That is
#   the most likely reason the antipolar test still fails, and it is a limit of what
#   *intramolecular* data can say about charges: conformer energies within one molecule
#   barely constrain the electrostatic scale, and a crystal's polarization is precisely
#   what would.  See ``REF_ABLATIONS`` and docs/DFT_FIT.md section 7.
REF_FITTED_X = np.array([
    1.7402782671587158,
    -1.7413263287763914,
    2.5637472343117422,
    1.7365311864572492,
    -1.575841782083297,
    1.6994672735105334,
    -0.10502580323091212,
    -1.0164094798299586,
    1.6729222712411376,
    1.0719973718259224,
    -1.0206766363085034,
    2.302350263500768,
    0.9000000000000001,
    0.900000000000494,
    1.0592165135041003,
    0.9927143193267619,
    0.5000000000000001,
    0.5000000000013645,
    0.5755166899867618,
    1.029994372637237,
    -0.25945214320699067,
    0.020000000000000004,
    0.04969246744992755,
    0.11119545724160368,
    0.03478498537117752,
    -0.10112507882084698,
    0.9396060965909317,
])

FITTED_DFT = reference_ff_parameters(REF_FITTED_X)
register_preset("pvdf-dft-fit", FITTED_DFT)

# Which group of parameters did the work?  Measured, with each group's own value taken
# from the full fit and everything else left illustrative:
#
#   parameter set        train E  held E   E(a)-E(b)   RIS top   E(anti)-E(polar)  beta |P|
#   illustrative           3.156   3.316    +7.38 kJ     TG+           +0.197        0.140
#   torsions only          2.614   3.115    +2.61        TG+           +0.197        0.140
#   Lennard-Jones only     2.481   2.528    +2.58       G+G+           +0.059        0.144
#   charges only           2.752   3.182    +4.64        TG+           +0.220        0.112
#   1-4 scale only         3.060   3.360    +7.38        TG+           +0.197        0.140
#   full fit               1.498   1.757    -4.54        TG+           +0.091        0.114
#
# Two things to take from it.  First, every group helps and none of them is close to the
# whole: the held-out error of the full fit (1.757) is far below the best single group
# (2.528), so this is not the crystal fit's pattern of one knob per target.  Second, and
# less comfortably, **different groups carry different acceptance tests and no setting
# carries both**: the Lennard-Jones parameters alone are what dislodge the 3/1 helix from
# the top of the RIS ranking and what halve the antipolar cost, while it takes the whole
# vector -- torsions included -- to get the alpha/beta ordering right, and doing so puts
# TG+ back on top.  A potential of this form appears not to be able to have both at once.
REF_ABLATIONS: dict[str, np.ndarray] = {}


def reference_ablations() -> dict[str, np.ndarray]:
    """Parameter vectors with one group of :data:`REF_VARIABLES` fitted and the rest at
    their illustrative values, for asking which group transfers."""
    groups = {
        "torsions only": np.arange(3 * _NT),
        "Lennard-Jones only": np.arange(3 * _NT, 3 * _NT + 2 * _NL),
        "charges only": np.arange(3 * _NT + 2 * _NL, 3 * _NT + 2 * _NL + _NI),
        "1-4 scale only": np.array([len(REF_X0) - 1]),
    }
    out = {}
    for name, idx in groups.items():
        v = REF_X0.copy()
        v[idx] = REF_FITTED_X[idx]
        out[name] = v
    return out


# =====================================================================================
# Valence terms, full Cartesian forces and off-site charges (docs/VALENCE_FIT.md)
# =====================================================================================
#
# DESIGN.md 5.9 ends by saying where the binding constraint had moved to: not the data any
# more but the functional form, which lacked (i) bond and angle terms, without which 93% of
# the reference force signal is unrepresentable and an all-trans chain cannot relieve a
# contact by opening an angle, and (ii) electrostatics richer than fixed atom-centred point
# charges, with the fit wanting a *positive* fluorine as the evidence.  This section is
# both, fitted together:
#
# * harmonic bond stretching and angle bending, typed by the bond graph
#   (:data:`VAL_BOND_TYPES`, :data:`VAL_ANGLE_TYPES`), with every equilibrium value and
#   every stiffness a fitted parameter rather than an invented one;
# * an **off-site charge for fluorine**: the charge slides along the C-F bond while the
#   Lennard-Jones centre stays on the nucleus, which is one parameter (``d_F``) and the one
#   thing atom-centred charges structurally cannot do -- the same point had to carry
#   fluorine's van der Waals size and the C-F dipole's centroid, and the fit resolved that
#   conflict by switching the dipole off;
# * the objective's force block is now the **full Cartesian forces**, every component of
#   every atom, not the torsional projection the rigid form had to fall back on.
#
# The first 27 entries of the vector are exactly :data:`REF_VARIABLES` in the same order,
# so "the previous fit" is this fit with the new blocks switched off, and
# :func:`valence_ablations` can say what each block bought.

# Bond stretch types (:func:`polyfind.forcefield.bond_type_name`), chosen by coverage of
# the reference set: C-C 4561 bonds / 31 systems, C-H 5081 / 31, C-F 5328 / 31, C-Cl 80 / 4,
# C-O 80 / 3, C=N (the nitrile, typed apart from single C-N because its length is 1.14 A
# against 1.47 and no one harmonic covers both) 171 / 10.  The wildcard takes the 235
# remaining bonds -- C-S, S=O, N=O, C-I, C-Br, the isocyanate's C-N -- which belong to
# two or three systems of nuisance data each and do not deserve a type of their own.
VAL_BOND_TYPES: tuple[str, ...] = ("C-C", "C-H", "C-F", "C-Cl", "C-O", "C=N", "*")

# Angle bend types (:func:`polyfind.forcefield.angle_type_name`).  The six carbon-centred
# ones cover 95.7% of the 29,954 angles in the set; the wildcard takes the rest.  Angles
# about a two-coordinate centre (a nitrile carbon) get no term at all -- see
# ``angle_type_name`` for why a harmonic in the angle is the wrong form there.
VAL_ANGLE_TYPES: tuple[str, ...] = ("C-C-C", "C-C-H", "C-C-F", "H-C-H", "F-C-F", "F-C-H", "*")

# The element whose charge is allowed off its nucleus.  Fluorine, because it is the element
# the diagnosis is about and the one this package's polarization results turn on; and one
# element means one parameter, which is the smallest change that can answer the question.
VAL_OFFSET_ELEMENT = "F"

_NB, _NA = len(VAL_BOND_TYPES), len(VAL_ANGLE_TYPES)

VAL_VARIABLES: tuple[str, ...] = tuple(
    list(REF_VARIABLES)
    + [f"kb_{t}" for t in VAL_BOND_TYPES] + [f"r0_{t}" for t in VAL_BOND_TYPES]
    + [f"ka_{t}" for t in VAL_ANGLE_TYPES] + [f"t0_{t}" for t in VAL_ANGLE_TYPES]
    + [f"d_{VAL_OFFSET_ELEMENT}"])

# Starting values for the valence block: textbook (GAFF/AMBER-like) stiffnesses and
# equilibrium lengths and angles, *not* averages of the reference geometries -- starting at
# the data's own means would make the weak ridge toward X0 a second, hidden fit to it.  In
# the convention used here the energy is 0.5 k (r - r0)^2, so k is twice the "K" of a force
# field that writes K (r - r0)^2.
VAL_BOND_K0 = (620.0, 680.0, 740.0, 440.0, 640.0, 1400.0, 500.0)
VAL_BOND_R0 = (1.526, 1.090, 1.380, 1.766, 1.410, 1.150, 1.600)
VAL_ANGLE_K0 = (120.0, 100.0, 120.0, 70.0, 160.0, 100.0, 110.0)
VAL_ANGLE_T0 = (112.7, 109.5, 109.0, 107.8, 107.0, 107.0, 109.5)
VAL_OFFSET0 = 0.0  # start with the charge on the nucleus: the old model is the starting point

VAL_X0 = np.concatenate([REF_X0, VAL_BOND_K0, VAL_BOND_R0, VAL_ANGLE_K0, VAL_ANGLE_T0,
                         [VAL_OFFSET0]])

# Bounds.  The first 27 are REF_BOUNDS unchanged, so a comparison with the previous fit is a
# comparison of forms and not of boxes.  For the new block:
#
# * stretch stiffness 100-3000 and bend stiffness 10-500 kcal/(mol A^2 | rad^2): loose
#   enough to be uninformative for the types with thousands of bonds behind them, present
#   only so the wildcard type cannot run away on 200 bonds of nuisance data;
# * equilibrium lengths 0.9-2.6 A and angles 90-135 deg: the range of a real single bond and
#   of a real bond angle.  These are the parameters the fit is being *asked* for, so the
#   bounds are deliberately far from where any of them is expected to land;
# * the fluorine charge offset +/- 0.6 A, up to about 40% of a C-F bond in either direction.
#   Zero is the atom-centred model, so the sign the fit chooses is the physical statement:
#   negative pulls the charge back towards carbon (a shorter dipole than the nuclei
#   suggest), positive pushes it beyond the nucleus (a longer one).
VAL_BOUNDS = np.vstack([REF_BOUNDS,
                        np.array([(100.0, 3000.0)] * _NB + [(0.90, 2.60)] * _NB
                                 + [(10.0, 500.0)] * _NA + [(90.0, 135.0)] * _NA
                                 + [(-0.60, 0.60)])])
VAL_SIGMA = np.concatenate([REF_SIGMA, [300.0] * _NB, [0.08] * _NB,
                            [80.0] * _NA, [4.0] * _NA, [0.15]])

# How many kcal^2/mol^2 of energy error one kcal^2/(mol A)^2 of *Cartesian force* error is
# worth.  Chosen by held-out energy error, measured in docs/VALENCE_FIT.md; unlike the rigid
# fit's torque weight, this one is not a formality.
VAL_FORCE_WEIGHT = 0.03


def valence_unpack(x) -> dict:
    """The fit vector as the arrays :meth:`ValenceDesign.evaluate` wants.

    The first 27 entries go through :func:`reference_unpack` unchanged -- same parameters,
    same order, same meaning -- and the rest are the valence block and the charge offset.
    """
    x = np.asarray(x, dtype=float)
    p = reference_unpack(x[:len(REF_X0)])
    off = len(REF_X0)
    p["bond_k"] = x[off:off + _NB]
    p["bond_r0"] = x[off + _NB:off + 2 * _NB]
    off += 2 * _NB
    p["angle_k"] = x[off:off + _NA]
    p["angle_t0"] = np.deg2rad(x[off + _NA:off + 2 * _NA])
    p["offset"] = float(x[-1])
    return p


def valence_ff_parameters(x, eps_r: float = 1.0) -> FFParameters:
    """The fit vector as an :class:`FFParameters`: valence terms and charge offset included."""
    base = reference_ff_parameters(np.asarray(x, dtype=float)[:len(REF_X0)], eps_r=eps_r)
    p = valence_unpack(x)
    bonds = tuple((t, float(k), float(r0))
                  for t, k, r0 in zip(VAL_BOND_TYPES, p["bond_k"], p["bond_r0"]))
    angles = tuple((t, float(k), float(np.degrees(t0)))
                   for t, k, t0 in zip(VAL_ANGLE_TYPES, p["angle_k"], p["angle_t0"]))
    offs = () if p["offset"] == 0.0 else ((VAL_OFFSET_ELEMENT, float(p["offset"])),)
    return replace(base, bond_terms=bonds, angle_terms=angles, charge_offsets=offs)


def _unit_rows(v: np.ndarray) -> np.ndarray:
    return v / np.linalg.norm(v, axis=-1, keepdims=True)


def _angle_gradient(coords, idx: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """``(theta, dtheta/dx)`` for angle triples ``idx`` (n, 3): radians and (n, 3, 3).

    Analytic: with ``u = x_i - x_j``, ``v = x_k - x_j`` and ``s = sin theta``,
    ``dtheta/dx_i = (cos theta * uhat - vhat) / (|u| s)``, the mirror expression for
    ``x_k``, and minus their sum for the central atom -- which is what makes the gradient of
    a rigid translation vanish.  ``tests/test_fitting.py`` checks it against the
    finite-difference forces of :meth:`polyfind.forcefield.SimpleFF.forces_frame`.
    """
    u = coords[idx[:, 0]] - coords[idx[:, 1]]
    v = coords[idx[:, 2]] - coords[idx[:, 1]]
    nu, nv = np.linalg.norm(u, axis=1), np.linalg.norm(v, axis=1)
    uh, vh = u / nu[:, None], v / nv[:, None]
    cos = np.clip((uh * vh).sum(1), -1.0, 1.0)
    theta = np.arccos(cos)
    s = np.maximum(np.sqrt(1.0 - cos ** 2), 1e-9)
    gi = (cos[:, None] * uh - vh) / (nu * s)[:, None]
    gk = (cos[:, None] * vh - uh) / (nv * s)[:, None]
    return theta, np.stack([gi, -(gi + gk), gk], axis=1)


def _dihedral_gradient(coords, idx: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """``(phi, dphi/dx)`` for torsion quadruples ``idx`` (n, 4): radians and (n, 4, 3).

    The standard result for the dihedral of ``i-j-k-l``, in the same sign convention as
    :func:`polyfind.forcefield.dihedral_angles` (IUPAC, trans = 180).  Checked against
    finite differences in the tests, because a sign error here would quietly corrupt every
    fitted torsion coefficient.
    """
    b1 = coords[idx[:, 1]] - coords[idx[:, 0]]
    b2 = coords[idx[:, 2]] - coords[idx[:, 1]]
    b3 = coords[idx[:, 3]] - coords[idx[:, 2]]
    n1, n2 = np.cross(b1, b2), np.cross(b2, b3)
    nb2 = np.linalg.norm(b2, axis=1)
    s1, s2 = (n1 * n1).sum(1), (n2 * n2).sum(1)
    phi = np.arctan2((nb2[:, None] * b1 * n2).sum(1), (n1 * n2).sum(1))
    gi = -nb2[:, None] * n1 / s1[:, None]
    gl = nb2[:, None] * n2 / s2[:, None]
    f = ((b1 * b2).sum(1) / nb2 ** 2)[:, None]
    g = ((b3 * b2).sum(1) / nb2 ** 2)[:, None]
    gj = -(1.0 + f) * gi + g * gl
    gk = f * gi - (1.0 + g) * gl
    return phi, np.stack([gi, gj, gk, gl], axis=1)


class ValenceDesign:
    """Every geometry-dependent quantity of a set of frames, for the valence fit.

    The same idea as :class:`ReferenceDesign` -- the frames never move, so distances,
    angles, dihedrals and all their Cartesian derivatives are constants and only the
    parameters change -- extended with what the new terms need: bond and angle lists with
    their analytic gradients, the full Cartesian reference forces, and the geometry of the
    off-site charge.

    The charge offset deserves a note, because it is the one parameter that moves an
    interaction *site*.  Writing the site separation as ``A + d U``, with ``A`` the nuclear
    separation and ``U`` the difference of the two unit offset directions, the Coulomb
    distance is ``sqrt(c0 + c1 d + c2 d^2)`` for three constants per pair -- so a parameter
    that changes the geometry still costs no geometry rebuild.  The Cartesian force needs
    the chain rule through the site's own dependence on its two atoms, which is a constant
    3x3 projector per off-site atom; both are precomputed here.

    :meth:`evaluate` reproduces :meth:`polyfind.forcefield.SimpleFF.energy_frame` and
    :meth:`~polyfind.forcefield.SimpleFF.forces_frame` to rounding, and
    ``tests/test_fitting.py`` asserts it rather than trusting it.
    """

    def __init__(self, frames: Sequence[Frame]):
        frames = list(frames)
        if not frames:
            raise ValueError("no frames")
        self.frames = frames
        self.systems = sorted({f.system for f in frames})
        sidx = {s: i for i, s in enumerate(self.systems)}
        self.elements = sorted(_UFF_ALL)
        el_index = {e: i for i, e in enumerate(self.elements)}
        b_type = {t: i for i, t in enumerate(VAL_BOND_TYPES)}
        a_type = {t: i for i, t in enumerate(VAL_ANGLE_TYPES)}

        p_i, p_j, p_s14, p_r, p_frame, p_A, p_U = [], [], [], [], [], [], []
        b_i, b_j, b_t, b_r, b_u = [], [], [], [], []
        a_idx, a_t, a_theta, a_grad = [], [], [], []
        d_frame, d_type, d_phi, d_idx, d_grad, d_tau = [], [], [], [], [], []
        r_dih, r_atom, r_vec = [], [], []
        o_atom, o_nb, o_proj = [], [], []
        a_el, a_M, f_ref, energies, f_sys, atom_frame = [], [], [], [], [], []
        n_atoms = n_dih = 0
        for fi, f in enumerate(frames):
            base = n_atoms
            c = f.coords
            a_el.extend(el_index[e] for e in f.elements)
            atom_frame.extend([fi] * f.n_atoms)
            a_M.append(_increment_matrix(f))
            f_ref.append(np.zeros((f.n_atoms, 3)) if f.forces is None else f.forces)
            # --- nonbonded pairs, and the off-site geometry they need
            dist = f.graph_distances()
            iu, ju = np.triu_indices(f.n_atoms, 1)
            keep = dist[iu, ju] >= 3
            iu, ju = iu[keep], ju[keep]
            p_i.append(iu + base)
            p_j.append(ju + base)
            p_s14.append(dist[iu, ju] == 3)
            p_r.append(np.linalg.norm(c[iu] - c[ju], axis=1))
            p_frame.append(np.full(len(iu), fi))
            u_dir = np.zeros((f.n_atoms, 3))
            for at, nb, _ in offset_sites(f.elements, f.bonds, {VAL_OFFSET_ELEMENT: 1.0}):
                rv = c[at] - c[nb]
                n = float(np.linalg.norm(rv))
                uh = rv / n
                u_dir[at] = uh
                o_atom.append(at + base)
                o_nb.append(nb + base)
                o_proj.append((np.eye(3) - np.outer(uh, uh)) / n)
            p_A.append(c[iu] - c[ju])
            p_U.append(u_dir[iu] - u_dir[ju])
            # --- valence lists
            bl, al = valence_topology(f.elements, f.bonds)
            if bl:
                ii = np.array([b[0] for b in bl])
                jj = np.array([b[1] for b in bl])
                b_i.append(ii + base)
                b_j.append(jj + base)
                b_t.append(np.array([b_type.get(b[2], b_type["*"]) for b in bl]))
                b_r.append(np.linalg.norm(c[ii] - c[jj], axis=1))
                b_u.append(_unit_rows(c[ii] - c[jj]))
            if al:
                idx = np.array([[a[0], a[1], a[2]] for a in al])
                theta, grad = _angle_gradient(c, idx)
                a_idx.append(idx + base)
                a_t.append(np.array([a_type.get(a[3], a_type["*"]) for a in al]))
                a_theta.append(theta)
                a_grad.append(grad)
            # --- torsions: angle, type, Cartesian gradient, and the rotor projection
            tors = np.asarray(f.torsions, dtype=int).reshape(-1, 4)
            if len(tors):
                phi, grad = _dihedral_gradient(c, tors)
                d_frame.append(np.full(len(tors), fi))
                d_type.append(torsion_types(f))
                d_phi.append(phi)
                d_idx.append(tors + base)
                d_grad.append(grad)
                d_tau.append(f.reference_torques() if f.forces is not None else np.zeros(len(tors)))
                rot = f.rotors()
                for k in range(len(tors)):
                    nz = np.flatnonzero(np.abs(rot[k]).sum(1) > 1e-12)
                    r_dih.append(np.full(len(nz), n_dih + k))
                    r_atom.append(nz + base)
                    r_vec.append(rot[k][nz])
                n_dih += len(tors)
            n_atoms += f.n_atoms
            energies.append(f.energy)
            f_sys.append(sidx[f.system])

        cat = lambda xs: (np.concatenate(xs) if xs else np.empty(0))  # noqa: E731
        self.n_frames, self.n_atoms, self.n_torsions = len(frames), n_atoms, n_dih
        self.pair_i, self.pair_j = np.concatenate(p_i), np.concatenate(p_j)
        self.pair_s14 = np.concatenate(p_s14)
        self.pair_r = np.concatenate(p_r)
        self.inv_r = 1.0 / self.pair_r
        A, U = np.concatenate(p_A), np.concatenate(p_U)
        self.pair_rhat = _unit_rows(A)
        self.pair_frame = np.concatenate(p_frame).astype(int)
        self.site_A, self.site_U = A, U
        self.site_c0 = (A * A).sum(1)
        self.site_c1 = 2.0 * (A * U).sum(1)
        self.site_c2 = (U * U).sum(1)
        self.bond_i = cat(b_i).astype(int)
        self.bond_j = cat(b_j).astype(int)
        self.bond_type = cat(b_t).astype(int)
        self.bond_r = cat(b_r)
        self.bond_u = np.concatenate(b_u) if b_u else np.empty((0, 3))
        self.angle_idx = (np.concatenate(a_idx) if a_idx else np.empty((0, 3), dtype=int)).astype(int)
        self.angle_type = cat(a_t).astype(int)
        self.angle_theta = cat(a_theta)
        self.angle_grad = np.concatenate(a_grad) if a_grad else np.empty((0, 3, 3))
        self.dih_frame = cat(d_frame).astype(int)
        self.dih_type = cat(d_type).astype(int)
        self.dih_phi = cat(d_phi)
        self.dih_idx = (np.concatenate(d_idx) if d_idx else np.empty((0, 4), dtype=int)).astype(int)
        self.dih_grad = np.concatenate(d_grad) if d_grad else np.empty((0, 4, 3))
        self.torque_ref = cat(d_tau)
        self.rot_dih = cat(r_dih).astype(int)
        self.rot_atom = cat(r_atom).astype(int)
        self.rot_vec = np.concatenate(r_vec) if r_vec else np.empty((0, 3))
        self.off_atom = np.array(o_atom, dtype=int)
        self.off_nb = np.array(o_nb, dtype=int)
        self.off_proj = np.array(o_proj).reshape(-1, 3, 3)
        self.atom_el = np.array(a_el, dtype=int)
        self.atom_frame = np.array(atom_frame, dtype=int)
        self.increment_matrix = np.vstack(a_M)
        self.force_ref = np.vstack(f_ref)
        self.energy_ref = np.array(energies)
        self.frame_system = np.array(f_sys, dtype=int)
        self.system_count = np.bincount(self.frame_system, minlength=len(self.systems))
        self.energy_ref_centred = self.centre(self.energy_ref)

    # --- the objective's pieces -------------------------------------------------
    def centre(self, x: np.ndarray) -> np.ndarray:
        """Subtract each system's own mean (see :meth:`ReferenceDesign.centre`)."""
        m = np.bincount(self.frame_system, weights=x, minlength=len(self.systems)) / self.system_count
        return x - m[self.frame_system]

    def _scatter(self, idx_list, val_list) -> np.ndarray:
        """Sum per-interaction 3-vectors onto atoms: one bincount per Cartesian component."""
        idx = np.concatenate(idx_list)
        val = np.concatenate(val_list)
        return np.stack([np.bincount(idx, weights=val[:, c], minlength=self.n_atoms)
                         for c in range(3)], axis=1)

    def evaluate(self, p: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """``(energies, forces, torques)``: per frame (kcal/mol), per atom (kcal/(mol A))
        and per backbone torsion (kcal/(mol rad), the rotor projection of those forces)."""
        x = p["lj_x"][self.atom_el]
        d = p["lj_d"][self.atom_el]
        q = self.increment_matrix @ p["increments"]
        scale = np.where(self.pair_s14, p["scale14"], 1.0)
        lj_x = np.sqrt(x[self.pair_i] * x[self.pair_j])
        lj_d = np.sqrt(d[self.pair_i] * d[self.pair_j]) * scale
        qq = q[self.pair_i] * q[self.pair_j] * scale * COULOMB
        s6 = (lj_x * self.inv_r) ** 6
        e_lj = lj_d * (s6 * s6 - 2.0 * s6)
        dlj_dr = (12.0 * lj_d * self.inv_r) * (s6 - s6 * s6)  # d(LJ)/dr
        f_lj = (-dlj_dr)[:, None] * self.pair_rhat
        # Coulomb between charge sites: the separation is A + d U, so the distance and its
        # direction both move with the one offset parameter and no geometry is rebuilt.
        dq = p["offset"]
        if dq == 0.0:
            r_c, sep = self.pair_r, self.site_A
        else:
            r_c = np.sqrt(self.site_c0 + dq * self.site_c1 + dq * dq * self.site_c2)
            sep = self.site_A + dq * self.site_U
        e_c = qq / r_c
        f_c = (qq / r_c ** 3)[:, None] * sep  # force on site i, i.e. -dE/dS_i
        energy = (np.bincount(self.pair_frame, weights=e_lj, minlength=self.n_frames)
                  + np.bincount(self.pair_frame, weights=e_c, minlength=self.n_frames))
        # Coulomb forces act on sites; map them onto atoms through the site's own Jacobian
        forces = self._scatter([self.pair_i, self.pair_j], [f_c, -f_c])
        idx_list, val_list = [self.pair_i, self.pair_j], [f_lj, -f_lj]
        if dq != 0.0 and self.off_atom.size:
            corr = dq * np.einsum("nab,nb->na", self.off_proj, forces[self.off_atom])
            idx_list += [self.off_atom, self.off_nb]
            val_list += [corr, -corr]
        forces = forces + self._scatter(idx_list, val_list)
        # --- valence
        if self.bond_i.size:
            k, r0 = p["bond_k"][self.bond_type], p["bond_r0"][self.bond_type]
            dr = self.bond_r - r0
            energy += np.bincount(self.atom_frame[self.bond_i], weights=0.5 * k * dr * dr,
                                  minlength=self.n_frames)
            fb = (-k * dr)[:, None] * self.bond_u
            forces += self._scatter([self.bond_i, self.bond_j], [fb, -fb])
        if self.angle_idx.size:
            k, t0 = p["angle_k"][self.angle_type], p["angle_t0"][self.angle_type]
            dt = self.angle_theta - t0
            energy += np.bincount(self.atom_frame[self.angle_idx[:, 1]], weights=0.5 * k * dt * dt,
                                  minlength=self.n_frames)
            g = (-k * dt)[:, None, None] * self.angle_grad  # (n, 3 atoms, 3 components)
            forces += self._scatter([self.angle_idx[:, c] for c in range(3)],
                                    [g[:, c] for c in range(3)])
        # --- torsion
        V = p["torsion"][self.dih_type]
        phi = self.dih_phi
        e_t = 0.5 * (V[:, 0] * (1 + np.cos(phi)) + V[:, 1] * (1 - np.cos(2 * phi))
                     + V[:, 2] * (1 + np.cos(3 * phi)))
        dt_dphi = 0.5 * (-V[:, 0] * np.sin(phi) + 2 * V[:, 1] * np.sin(2 * phi)
                         - 3 * V[:, 2] * np.sin(3 * phi))
        energy += np.bincount(self.dih_frame, weights=e_t, minlength=self.n_frames)
        gt = (-dt_dphi)[:, None, None] * self.dih_grad
        forces += self._scatter([self.dih_idx[:, c] for c in range(4)], [gt[:, c] for c in range(4)])
        torque = np.bincount(self.rot_dih,
                             weights=(forces[self.rot_atom] * self.rot_vec).sum(1),
                             minlength=self.n_torsions)
        return energy, forces, torque

    def residuals(self, x, force_weight: float = VAL_FORCE_WEIGHT,
                  ridge: float = REF_RIDGE) -> np.ndarray:
        """Relative energies, full Cartesian forces and the ridge, in the usual scaled units.

        Each block is divided by the square root of its own length, so the energy block's sum
        of squares is a mean squared error in kcal^2/mol^2 and ``force_weight`` reads as "how
        many of those one kcal^2/(mol A)^2 of force error is worth".
        """
        p = valence_unpack(x)
        energy, forces, _ = self.evaluate(p)
        blocks = [(self.centre(energy) - self.energy_ref_centred) / np.sqrt(self.n_frames)]
        if force_weight > 0:
            blocks.append(np.sqrt(force_weight) * (forces - self.force_ref).ravel()
                          / np.sqrt(3 * self.n_atoms))
        if ridge > 0:
            blocks.append(np.sqrt(ridge / len(x)) * (np.asarray(x, dtype=float) - VAL_X0) / VAL_SIGMA)
        return np.concatenate(blocks)

    def errors(self, x) -> dict:
        """Energy, force and torque root-mean-square errors, overall and per system."""
        p = valence_unpack(x)
        energy, forces, torque = self.evaluate(p)
        de = self.centre(energy) - self.energy_ref_centred
        df = forces - self.force_ref
        dt = torque - self.torque_ref
        per, per_f = {}, {}
        for i, s in enumerate(self.systems):
            m = self.frame_system == i
            per[s] = float(np.sqrt((de[m] ** 2).mean()))
            sel = m[self.atom_frame]
            per_f[s] = float(np.sqrt((df[sel] ** 2).mean()))
        return {"energy_rms": float(np.sqrt((de ** 2).mean())),
                "energy_max": float(np.abs(de).max()),
                "force_rms": float(np.sqrt((df ** 2).mean())),
                "torque_rms": float(np.sqrt((dt ** 2).mean())),
                "force_ref_rms": float(np.sqrt((self.force_ref ** 2).mean())),
                "force_corr": float(np.corrcoef(forces.ravel(), self.force_ref.ravel())[0, 1]),
                "per_system": per, "per_system_force": per_f}


@dataclass
class ValenceFit:
    """One fit of :data:`VAL_VARIABLES` to reference data, with its held-out score."""

    x: np.ndarray
    force_weight: float
    train_systems: list[str]
    test_systems: list[str]
    train: dict
    test: dict
    n_evaluations: int
    seconds: float

    @property
    def parameters(self) -> FFParameters:
        return valence_ff_parameters(self.x)

    def table(self) -> str:
        rows = [f"{'system':<12s} {'E rms':>9s} {'F rms':>9s}  set", "-" * 44]
        for s in self.train_systems:
            if s in self.train["per_system"]:
                rows.append(f"{s:<12s} {self.train['per_system'][s]:9.3f} "
                            f"{self.train['per_system_force'][s]:9.3f}  train")
        for s in self.test_systems:
            if s in self.test["per_system"]:
                rows.append(f"{s:<12s} {self.test['per_system'][s]:9.3f} "
                            f"{self.test['per_system_force'][s]:9.3f}  HELD OUT")
        rows.append("-" * 44)
        rows.append(f"{'train':<12s} {self.train['energy_rms']:9.3f} {self.train['force_rms']:9.3f}")
        rows.append(f"{'held out':<12s} {self.test['energy_rms']:9.3f} {self.test['force_rms']:9.3f}")
        return "\n".join(rows)


def fit_valence(frames: Sequence[Frame] | None = None, force_weight: float = VAL_FORCE_WEIGHT,
                path: str | None = None, n_starts: int = 3, seed: int = 0,
                free: Sequence[int] | None = None, x0=None, ridge: float = REF_RIDGE,
                designs: tuple["ValenceDesign", "ValenceDesign"] | None = None,
                verbose: bool = True) -> ValenceFit:
    """Fit :data:`VAL_VARIABLES` to the reference set, holding whole chemistries out.

    The same optimiser, the same split and the same ridge as :func:`fit_reference` -- what
    changes is the model (valence terms, off-site charge) and the force block (full
    Cartesian forces instead of the torsional projection).  ``free`` restricts the fit to a
    subset of the parameter indices, the rest held at ``x0``, which is how
    :func:`valence_ablations` asks what each block bought; ``designs`` reuses a
    (train, held-out) pair that has already been built, which is worth doing because
    building one is a second of work and a sweep builds none of it twice.
    """
    if designs is None:
        if frames is None:
            frames = load_reference(path)
        frames = list(frames)
        train_sys, test_sys = split_systems(frames)
        design = ValenceDesign([f for f in frames if f.system in train_sys])
        held = ValenceDesign([f for f in frames if f.system in test_sys])
    else:
        design, held = designs
        train_sys, test_sys = list(design.systems), list(held.systems)
    x0 = VAL_X0 if x0 is None else np.asarray(x0, dtype=float)
    idx = np.arange(len(VAL_X0)) if free is None else np.asarray(free, dtype=int)
    lo, hi = VAL_BOUNDS[idx, 0], VAL_BOUNDS[idx, 1]
    calls = {"n": 0}

    def expand(v):
        full = x0.copy()
        full[idx] = v
        return full

    def f(v):
        calls["n"] += 1
        return design.residuals(expand(v), force_weight, ridge)

    rng = np.random.default_rng(seed)
    starts = [x0[idx]] + [np.clip(x0[idx] + rng.normal(0, VAL_SIGMA[idx]), lo, hi)
                          for _ in range(n_starts)]
    t0 = time.time()
    best, best_cost = None, np.inf
    for k, s in enumerate(starts):
        res = least_squares(f, s, bounds=(lo, hi), x_scale=VAL_SIGMA[idx], max_nfev=6000)
        if verbose:
            print(f"  start {k}: cost {res.cost:.6f} in {res.nfev} evaluations", flush=True)
        if res.cost < best_cost:
            best, best_cost = res.x.copy(), float(res.cost)
    x = expand(best)
    return ValenceFit(x=x, force_weight=force_weight, train_systems=train_sys,
                      test_systems=test_sys, train=design.errors(x), test=held.errors(x),
                      n_evaluations=calls["n"], seconds=time.time() - t0)


def rigid_vector(x27=None) -> np.ndarray:
    """A :data:`VAL_VARIABLES` vector that *is* the old rigid form: no valence, no offset.

    Zero stiffnesses and a zero charge offset reduce the new potential to exactly the one
    :data:`REF_VARIABLES` describes, which is what lets :class:`ValenceDesign` score the
    previous fit and this one on the same footing -- same frames, same centring, same
    Cartesian force block -- instead of comparing two numbers from two objectives.
    """
    x = VAL_X0.copy()
    if x27 is not None:
        x[:len(REF_X0)] = np.asarray(x27, dtype=float)
    off = len(REF_X0)
    x[off:off + _NB] = 0.0
    x[off + 2 * _NB:off + 2 * _NB + _NA] = 0.0
    x[-1] = 0.0
    return x


def valence_degeneracies(design: "ValenceDesign", x=None, s: float = 1.7) -> dict:
    """Measure the degeneracies of the *new* electrostatics rather than assuming them.

    Three questions, each answered by moving parameters and looking at what changes:

    * **charge scale against permittivity** -- exactly degenerate in the old form.  Moving
      the charge site does not touch that argument (every energy still depends on the
      charges only through ``q_i q_j``, and on ``eps_r`` only as a divisor), so the
      prediction is that it is *still* exact.  Reported as the largest change in any frame
      energy after multiplying every increment by ``s`` and dividing the Coulomb term by
      ``s^2``.
    * **charge against offset** -- new, and the one worth worrying about: at long range only
      the product ``q d`` (the dipole) enters, so a fit could trade one against the other.
      Reported as the largest energy change when ``q_C-F`` is scaled by ``s`` and ``d_F``
      divided by it, holding the C-F bond dipole fixed.  A number far above rounding means
      the short-range part separates them; a number near rounding means it does not.
    * **the near-linear angle guard** -- reported as the largest bend angle carrying a term,
      so that a reader can see the harmonic is never being asked to work near 180 degrees.
    """
    x = VAL_X0 if x is None else np.asarray(x, dtype=float)
    p = valence_unpack(x)
    e0, f0, _ = design.evaluate(p)
    no_q = dict(p)
    no_q["increments"] = np.zeros_like(p["increments"])
    e_bare = design.evaluate(no_q)[0]
    q = dict(p)
    q["increments"] = p["increments"] * s
    e1 = design.evaluate(q)[0]
    e1_bare = design.evaluate({**q, "increments": np.zeros_like(p["increments"])})[0]
    eps = float(np.abs((e1 - e1_bare) / s ** 2 - (e0 - e_bare)).max())
    # q_C-F against d_F at fixed bond dipole.  q d is what survives at long range; if the
    # short range separates them, holding q d fixed still moves the energy.
    i_qf = REF_VARIABLES.index("q_C-F")
    y = x.copy()
    y[i_qf] = x[i_qf] * s
    d0 = x[-1] if x[-1] != 0.0 else 0.10
    y[-1] = d0 / s
    z = x.copy()
    z[-1] = d0
    dip = float(np.abs(design.centre(design.evaluate(valence_unpack(y))[0])
                       - design.centre(design.evaluate(valence_unpack(z))[0])).max())
    return {"scale": s, "charge_permittivity_max_energy_change": eps,
            "charge_offset_max_energy_change": dip,
            "max_bent_angle_deg": float(np.degrees(design.angle_theta.max())),
            "n_angles": int(design.angle_theta.size), "n_bonds": int(design.bond_r.size)}


# ------------------------------------------------------------------ the result
#
# The fit ``examples/fit_valence.py`` reproduces: 56 parameters, force weight 0.03, four
# starts (all four land on the same cost to six figures), 21 training chemistries and 10
# held out, about 260 s.  Numbers and provenance in docs/VALENCE_FIT.md.  Like every other
# fit here it is a *named* alternative and never a default: ``SimpleFF()`` is still the
# illustrative rigid potential and nothing selects this unless a caller asks for it.
#
# READ THIS BEFORE USING IT.  Scored by the same design over the same frames:
#
#   model                                   train E  held E  train F  held F   r(F)
#   illustrative (rigid)                      3.156   3.316    13.93   14.37  +0.06
#   pvdf-dft-fit (rigid)                      1.498   1.757    13.70   14.06  +0.10
#   the rigid form refitted to these forces   1.580   1.995    13.48   13.86  +0.17
#   this fit                                  1.022   1.359     4.64    5.23  +0.93
#   this fit, force weight 0                  0.859   1.157    41.98   44.29  +0.77
#
# (energies kcal/mol, forces kcal/(mol A); the reference forces themselves have RMS 14.05
# held out, which is what predicting *zero* would score.)
#
# * **The valence terms are what made the forces usable.**  Same objective, same optimiser,
#   only the functional form differs: the rigid form fitted to full Cartesian forces still
#   scores 13.86 against the reference's own 14.05, i.e. no better than predicting nothing,
#   while this one scores 5.23 with a correlation of +0.93.  That is the 93%-of-the-signal
#   claim of docs/DFT_FIT.md section 3 cashed out.
# * **The forces still do not help the energies.**  Held-out energy is *better* without
#   them (1.157 against 1.359), and the sweep is monotone from 0.003 upwards.  The reason is
#   the one the previous fit already diagnosed and this fit now measures: the reference
#   geometries were relaxed at some other level and labelled with PBE-D3, so the forces
#   record a systematic bond-length disagreement -- the fitted C-F equilibrium comes out
#   0.11 A *longer* than the geometries and C-C 0.08 A shorter -- and reproducing that pulls
#   the same parameters the relative energies want elsewhere.  Weight 0.03 is shipped anyway
#   because a potential whose gradients are wrong is not usable by anything that takes one,
#   and ``refine_crystal`` does; the energy it costs, 0.2 kcal/mol, is stated rather than
#   hidden.
# * **Fluorine's charge is physical again**: ``q_C-F`` comes out +0.114 e, off its bound,
#   where the previous fit sat on the 0.02 floor and an unbounded version of it turned
#   fluorine positive.  Six of 56 parameters sit on a bound (``x_H``, the three well depths,
#   and two stretch stiffnesses at their floor), against five of 27 before -- but none of
#   the six is an electrostatic parameter now, which is the one that mattered.
# * **The off-site charge is nearly free and buys nearly nothing.**  With ``d_F`` held at
#   zero the held-out energy is 1.347 against 1.359, and ``q_C-F`` is +0.115 rather than
#   +0.114.  So it is the valence terms, not the off-site site, that made the electrostatics
#   behave; the offset earns its keep only as the thing that shows the charge is no longer
#   fighting the repulsion, and it is reported as a null result.
# * **Acceptance**: test 1 passes (-2.96 kJ/mol per monomer, inside -6.5 to -2.6, against
#   the previous fit's -4.54).  Test 2 still fails: the isolated-chain ranking still puts
#   TG+ first, though its margin over the next candidate falls from 1.64 (illustrative) to
#   0.23 (previous fit) to 0.18 here.  Test 3 still fails, and barely moved: +0.085 kcal/mol
#   per monomer against +0.091.  Two of three still fail and the electrostatics change did
#   not fix the one it was aimed at.
VAL_FITTED_X = np.array([
    -0.45895730339776053, -0.4822495465560268, 1.7238990579737528,      # V1-3 [CF2-CH2]
    -0.26594448635740686, -0.5101634076419305, 0.81360756199377,        # V1-3 [CF2-CHF]
    -1.7452730586791623, 0.10165433247810902, 1.7683612752865407,       # V1-3 [CF2-CF2]
    -0.1458670286033234, -0.6126486958804321, 1.5119133250038446,       # V1-3 [other]
    0.9265637493571174, 0.9000000000000001, 1.054557911972751, 0.9751545700904893,   # x
    0.5000000000000001, 0.5000000000000001, 0.5000000000000001, 1.062609075380776,   # D
    -0.11925595971322667, 0.11364575244497159, 0.12446382203875808,     # q C-H, C-F, C-Cl
    0.36992912294510616, 0.1038657564132131, -0.14782354168138115,      # q C-N, C-O, C-S
    0.9229004452808719,                                                 # scale14
    100.00000000000001, 544.5777681680139, 259.5906378151254, 415.5067435646117,
    520.699033323524, 1385.4077330667221, 100.0000000000219,            # kb
    1.4575381757261119, 1.0998399522884275, 1.4673820422843398, 1.7861811855100675,
    1.4230352561726245, 1.1933496481543355, 1.7567790900302314,         # r0
    77.10439889761436, 60.36021994798029, 78.62988041171613, 84.67686115253565,
    110.28664551905766, 129.53166130565626, 43.91212924577346,          # ka
    108.81342550717888, 111.02523291095487, 109.11136328641616, 110.92307973586375,
    109.90112341671663, 109.28749731412486, 107.3248298127191,          # theta0
    0.17745728103468422,                                                # d_F
])

FITTED_VALENCE = valence_ff_parameters(VAL_FITTED_X)
register_preset("pvdf-dft-valence", FITTED_VALENCE)


# --------------------------------------------------------------------- charge flux
#
# ``pvdf-dft-valence`` with the geometry dependence of its bond-charge increments turned
# on (:class:`polyfind.forcefield.FluxTopology`).  Everything else -- torsions,
# Lennard-Jones, the increments themselves, the valence terms, the off-site charge -- is
# untouched, and the flux is opt-in per packer (``CrystalPacker(charge_flux=...)``) exactly
# as the valence terms are, so ``pvdf-dft-valence`` and every number measured with it are
# unchanged.
#
# WHY IT EXISTS.  With fixed increments a planar all-trans zigzag's dipole is *exactly*
# independent of its backbone angle, and the backbone angle is the only internal coordinate
# an axial strain moves when bond lengths are rigid.  So beta-PVDF -- the phase the material
# is used in -- has ``d_33 = d_31 = 0`` identically, and no amount of relaxation changes it
# (docs/ELECTROMECHANICS.md 5.2).  A charge that moves with the geometry is the smallest
# change that breaks that symmetry, and measured on beta it does: per unit ``k_angle`` the
# cell dipole gains ``+3.543`` (C-H) and ``-6.469`` (C-F) e.A per unit axial strain, against
# a fixed-increment dipole that holds all ten printed digits over the same sweep.
#
# HOW IT WAS FITTED, AND WHAT THAT IS WORTH.  ``examples/fit_charge_flux.py`` reproduces it.
# The target is ``dmu/d(axial strain)`` for four chemistries (PVDF, VDCN, AN, CNEPO) from
# ``sarco/materials/gpu_bundle/results/field_neighborhood_refined``: twelve numbers, three
# vector components each.  **Two parameters for twelve observations, 6:1** -- deliberately,
# after DESIGN.md 5.7, where five parameters on twelve dependent observations fitted noise.
#
#   * The reference is an **exploratory GFN2-xTB** calculation on a *finite* two-chain pair
#     in vacuum with the terminal backbone atoms pinned, and its own ``interpretation``
#     field says "finite-size, packing, stereochemistry and higher-level DFT validation
#     outstanding.  Strain is relative to fixed seed span, not a stress-free bulk lattice."
#     It calibrates a mechanism and an order of magnitude, not a bulk coefficient.  It is
#     also not the PBE-D3 the rest of this potential is fitted to.
#   * **The fit is poor: R^2 = 0.39** against the fixed-increment residual (rms 3.43 -> 2.68
#     e.A per unit strain).  Two things fit that residual *better* and neither is available
#     to this model: a bond-length flux channel (R^2 = 0.86 on two parameters), which is
#     inert here because ``build_chain`` places every atom at the polymer's own bond length,
#     and a plain charge-magnitude scale of x1.36 (R^2 = 0.75 on one), which is not a flux
#     at all but a statement that the increments are too small.
#   * **The angle coefficient is not separable from the bond one.**  Fitted with a bond
#     channel present it comes out ``k_angle(C-H) = -0.213`` instead of ``-1.214``, and the
#     sign of beta's ``d_31`` follows the choice: positive with this fit, negative with that
#     one.  What is shipped is the fit of *the model that is deployed* -- a rigid-bonded
#     chain has only the angle channel, so only the angle channel is fitted -- and R^2 =
#     0.39 is then honestly its own rather than borrowed from a channel that does nothing.
#   * Leave-one-chemistry-out moves ``k_angle(C-H)`` over ``-1.62 .. -0.15``, a factor of
#     eleven, though the *sign* of every converged beta ``d_31`` stays positive.
#   * The backbone C-C bond carries neither an increment nor a flux, because the increments
#     are typed by element pair and a homonuclear pair has no orientation
#     (:func:`polyfind.forcefield.flux_topology` refuses one).  Charge transfer *along* the
#     backbone is therefore outside this model, and it is the channel most likely to carry
#     an axial response.
#
# WHAT IT CHANGES IN THE CRYSTAL, which is not only the dipole.  The charges enter the
# Coulomb sum, so beta's relaxed structure moves: ``c`` 2.5469 -> 2.6052 A (+2.3%), ``a``
# 4.596 -> 4.492, ``|P_x|`` 0.116 -> 0.150 C/m^2, while ``C_33`` holds at 328 -> 330 GPa.
# The relaxation is a genuine stationary point (shape gradient 2e-9, residual axial stress
# 3e-8 GPa); at coefficients about a third larger the line group's +/-8 degree cap binds and
# the reference stops being stress-free, which ``examples/fit_charge_flux.py`` reports as a
# residual axial stress rather than hiding.
FLUX_ANGLE_CH = -1.2139  # e per unit of the dimensionless angle driver
FLUX_ANGLE_CF = -0.0975
FITTED_VALENCE_FLUX = replace(FITTED_VALENCE, charge_flux=(
    ("C", "H", FLUX_ANGLE_CH, 0.0),
    ("C", "F", FLUX_ANGLE_CF, 0.0),
))
register_preset("pvdf-dft-valence-flux", FITTED_VALENCE_FLUX)


# --------------------------------------------------------- charge flux fitted to Born charges
#
# The same two bond types, both channels, fitted to the periodic PBE-D3 DFPT Born effective
# charges of beta-PVDF (``examples/fit_born_flux.py``; docs/ELECTROMECHANICS.md 5.9).  The
# stretch channel's zero on the lattice path is the built chain's own bond length
# (:func:`polyfind.pack.built_bond_lengths`), so at the reference geometry these charges are
# ``pvdf-dft-valence``'s to the last digit and ``k_bond`` is a pure response.
#
# WHAT WAS FITTED.  The four transverse components (aa, bb; ab, ba for the pendants) of the
# four atom types, twelve numbers of which the acoustic sum rule -- satisfied identically by
# any neutral flux -- ties two: **four parameters for ten independent observations, 2.5:1**.
# The chain-axis components were kept out (the provider's acoustic-sum correction moves each
# of them by 0.052 e on a raw residual of -0.62 e) and are reported as a check.  The fit is
# linear (at fixed geometry the dipole is linear in every coefficient), rms 0.083 e on
# components of order 1 e, from 0.63 e with no flux and 1.62 e with the GFN2-fitted
# ``pvdf-dft-valence-flux``; raw and corrected reference agree to 1e-4 in every coefficient,
# transposing the off-diagonal convention moves them by 4e-3, and leave-one-type-out keeps
# ``k_bond(C-F)`` within 0.58..0.61.  Two iterations: fitted at the ``-flux`` structure
# (0.570, rms 0.076), then refitted at its own relaxed structure (these values), which moved
# no coefficient by more than the leave-one-out spread.
#
# WHAT IT SAYS.  Fluorine's dynamical charge is carried by charge flowing *along* the C-F
# bond as it stretches (``k_bond(C-F)`` = 0.57 e/A, Born charge -1.3 e along the bond
# against -0.4 e across it); the C-H bond carries almost none (hydrogen's is isotropic to
# 0.005 e in the reference).  The GFN2-fitted angle flux had put a dynamical charge of -4.4 e
# on the CH2 carbon where the reference has -0.12, and that channel -- ``k_angle(C-H)`` =
# -1.21, now +0.10 -- is what beta's proper ``d_33`` and ``d_31`` were made of.  Read
# section 5.9 before quoting a piezoelectric coefficient from this preset: the Born tensor is
# the *only* observable the stretch channel enters on a rigid-bond chain, so it cannot move
# ``d`` at all, and the refitted angle channel moves it the wrong way.
BORN_FLUX = (
    ("C", "H", 0.1119, 0.0422),  # (k_angle, k_bond)
    ("C", "F", 0.0963, 0.5849),
)
FITTED_VALENCE_FLUX_BORN = replace(FITTED_VALENCE, charge_flux=BORN_FLUX)
register_preset("pvdf-dft-valence-flux-born", FITTED_VALENCE_FLUX_BORN)


def valence_ablations() -> dict[str, np.ndarray]:
    """``free`` index sets for :func:`fit_valence`: which block of the new form did the work.

    Everything outside the named block is held at :data:`VAL_X0` -- the illustrative
    nonbonded potential plus *textbook* valence terms -- so these say what fitting a block
    buys, not what having it buys.  For the latter, :func:`rigid_vector` gives the form with
    no valence terms at all, and the note above ``VAL_FITTED_X`` records what that scores.

    Measured held-out energy / force RMS, one start each, force weight 0.03:

        nonbonded + torsion only          1.432 / 20.27   (valence present but not fitted)
        + valence, no off-site charge     1.347 /  5.24
        + off-site charge, no valence     1.435 / 20.27
        valence only, nonbonded left as is 2.729 /  7.08
        everything                        1.359 /  5.23

    The pattern is not the previous fit's.  There, different parameter groups carried
    different acceptance tests and none carried two; here one block -- the valence terms --
    carries the entire force result and most of the energy one, and the off-site charge is
    within noise of doing nothing at all.
    """
    off = len(REF_X0)
    groups = {
        "nonbonded + torsion only": np.arange(len(REF_X0)),
        "valence, no off-site charge": np.arange(len(VAL_X0) - 1),
        "off-site charge, no valence": np.concatenate([np.arange(len(REF_X0)), [len(VAL_X0) - 1]]),
        "valence only": np.arange(off, len(VAL_X0) - 1),
    }
    return {name: idx for name, idx in groups.items()}
