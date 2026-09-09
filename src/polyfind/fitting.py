"""Fit :class:`~polyfind.forcefield.SimpleFF` to the crystal data already in this repo.

Why this module exists is recorded in DESIGN.md section 2.5: what is wrong with
``SimpleFF`` is not that it is cheap, it is that it was never fitted, and no
machine-learned potential is or will be a dependency here.  The search around the
potential is now fast enough (a pack-and-refine of a reference chain is seconds) that
a parameter fit can sit on top of it, so that is what this module is.

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
* beta-PVDF's spontaneous polarization, about 0.13 C/m^2;
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
from typing import Iterable

import numpy as np
from scipy.optimize import minimize

from .forcefield import DEFAULT_TORSION, PRESETS, SimpleFF
from .pack import CrystalPacker, PackResult, pack, periodic_chain
from .pipeline import EXPERIMENTAL_CELLS
from .polymers import PE, PVDF, Polymer, THREE_STATE
from .refine import refine_crystal

T, GP, GM = 0, 1, 2

# Spontaneous polarizations (C/m^2).  Experimental value for beta-PVDF; alpha is the
# antipolar phase, whose polarization is zero by symmetry.  DESIGN.md section 5.5.
P_BETA_EXPERIMENT = 0.13


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

    @property
    def lj_dict(self) -> dict[str, tuple[float, float]]:
        return {e: (x, d) for e, x, d in self.lj}

    def simple_ff(self) -> SimpleFF:
        return SimpleFF(torsion=tuple(self.torsion), scale14=self.scale14, eps_r=self.eps_r,
                        charge_scale=self.charge_scale, lj=self.lj_dict or None)

    def preset_kwargs(self) -> dict:
        """:class:`~polyfind.forcefield.SimpleFF` keyword arguments, for :data:`PRESETS`."""
        kw = {"torsion": tuple(float(v) for v in self.torsion), "scale14": float(self.scale14),
              "eps_r": float(self.eps_r), "charge_scale": float(self.charge_scale)}
        if self.lj:
            kw["lj"] = {e: (float(x), float(d)) for e, x, d in self.lj}
        return kw

    def describe(self) -> str:
        lj = ", ".join(f"{e} x_i={x:.3f}" for e, x, _ in self.lj) or "UFF unchanged"
        return (f"torsion=({self.torsion[0]:+.3f}, {self.torsion[1]:+.3f}, {self.torsion[2]:+.3f}) "
                f"eps_r={self.eps_r:.3f} charge_scale={self.charge_scale:.3f} [{lj}]")

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
        """
        from . import pack as pack_mod
        from . import refine as refine_mod

        ff = self.simple_ff()
        q_scale, torsion, eps_r = self.charge_scale, tuple(self.torsion), self.eps_r

        class _FittedPacker(pack_mod.CrystalPacker):
            def __init__(self, chain, **kw):
                if q_scale != 1.0:
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


def antipolar_cell(packer: CrystalPacker, start: PackResult, maxfev: int = 250) -> tuple[np.ndarray, float]:
    """Best cell with chain 2 flipped *and* the two setting angles equal.

    That configuration is the antipolar one: flipping chain 2 reverses its axial dipole
    and equal setting angles cancel the transverse part, so the cell dipole is exactly
    zero for any (a, b, dz) -- asserted in the tests, not assumed.  It is a symmetric
    subspace of the seven cell parameters, so an unconstrained polish would slide out of
    it towards the polar minimum; tying phi2 to phi1 keeps the search inside it and
    measures what the antipolar phase costs under the current potential.

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
            refine: bool = True, verbose: bool = False) -> Prediction:
    """Pack and refine ``case`` under ``params``: the inner loop of the fit.

    Table screen (prebuilt, see :func:`screen_table`) -> exact-kernel polish -> continuous
    refinement of the torsions, the backbone angles and the cell.  For an antipolar phase
    the symmetric branch of :func:`antipolar_cell` is searched as well and the lower of
    the two rigid branches is the one that gets refined, so the fit cannot be fooled by a
    screen that only ever proposes polar starts.
    """
    from . import pack as pack_mod

    chain = reference_chain(case)
    table = screen_table(case)  # outside the patch, deliberately
    t0 = time.time()
    with params.applied():
        results = pack(chain, table=table, n_refine=n_refine, verbose=verbose)
        rigid = results[0]
        gap = None
        if case.antipolar:
            # resolved through the module so that the patched class is the one built
            packer = pack_mod.CrystalPacker(chain, n_chains=2)
            anti_params, anti_e = antipolar_cell(packer, rigid)
            gap = anti_e - rigid.energy_per_monomer
            if gap < 0.0:  # the antipolar cell is the ground state: refine that one
                rigid = packer.result(anti_params)
        if not refine:
            e = rigid.energy_per_monomer
            out = Prediction(case.key, rigid.a, rigid.b, rigid.c, rigid.density, e,
                             rigid.polarization_magnitude, gap, rigid, None, rigid.dihedrals)
            out.seconds = time.time() - t0
            return out
        ref = refine_crystal(case.polymer, rigid, eps_r=params.eps_r)
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
