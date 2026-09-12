"""Screen the electroactive-polymer campaign's chemistries, and first test whether the ranking means anything.

    python examples/screen_electroactive.py --ranking       # step 1: is the response ranking meaningful?
    python examples/screen_electroactive.py --screen        # step 2: the funnel on every expressible chemistry
    python examples/screen_electroactive.py --candidates    # step 3: new substituent combinations
    python examples/screen_electroactive.py --all --json out.json
    python examples/screen_electroactive.py --screen --coulomb ewald --only pvdf   # polarity under Ewald

``docs/SCREEN.md`` is written from this script's output and quotes it.

**The polarity column is reported against the model's own error bar, and mostly comes out
"cannot tell".**  ``--coulomb ewald`` puts *both* branches of the polar/antipolar comparison
on the Ewald sum of :mod:`polyfind.ewald` (screened with the truncated twin and polished with
Ewald, the division of labour ``pack(coulomb="ewald")`` draws); the refined cells, densities
and lattice ranking stay on the truncated sum, which is the protocol the earlier table was
measured with, so that one thing changes at a time.  Every gap is printed beside
:data:`ERROR_BAR`, the fitted potential's own held-out energy error per monomer, and a gap
inside it is reported as ``not resolved`` rather than as a verdict -- see :func:`verdict`.

**Step 1 is a gate, not a formality.**  ``docs/BENCHMARK.md`` says the piezoelectric
magnitudes are 5x and 9x short, from a charge-flux coefficient fitted to twelve numbers
with ``R^2 = 0.39``.  A screen survives that *only if the error is systematic*, and the
sibling project's ``results/field_neighborhood_refined`` is the one place that is testable:
four chemistries with an axial dipole response at two strains.  The shipped coefficient was
fitted to all four, so comparing against them directly is circular.  ``--ranking`` refits
with each chemistry held out, predicts the held-out one, and reports the Spearman
correlation of the predicted ordering against the reference ordering -- together with the
control that matters more than the fit itself: what the *fixed-increment* charges alone,
with no flux and nothing fitted to this data at all, predict for the same four.

**Step 2 screens what the package is good at.**  Exact conformational enumeration, exhaustive
packing and polymorph ordering are the trustworthy half (``docs/BENCHMARK.md``); the
response numbers are the weak half and are flagged as such in every table.  The single most
decision-relevant column is not a coefficient: it is whether a **polar chain conformation is
accessible at all**, and how far above the conformational ground state it sits.  A polymer
whose all-trans phase is 5 kcal/mol per monomer up is not a ferroelectric whatever its
coefficients say.

**What is out of scope, and why.**  CNEPO's epoxide bridges two backbone carbons, which is a
change of chain *topology* rather than a substituent; ``polyfind.polymers.BackboneAtom``
carries pendants, not rings, so CNEPO cannot be built (``docs/CHEMISTRY_EXTENSION.md``
phase 4).  It is reported as out of scope rather than approximated by something else.

**Potential.**  ``pvdf-dft-valence`` everywhere (the RIS fit, the packing, the refinement and
the deformable mechanics), with ``pvdf-dft-valence-flux`` added to the *mechanics* packer
only -- the same split ``examples/electromechanics.py`` uses and for the same reason.  It is
the best available fitted preset: it is the only one with valence terms, hence the only one
on which ``C_33`` and the diagonal columns of ``e`` are computable at all.  It is fitted to
**PVDF** single-chain PBE-D3 energies and forces, so applying it to a chlorine or nitrile
chemistry is a transfer assumption with nothing behind it but the wildcard bond and angle
types (``fitting.VAL_BOND_TYPES`` ends in ``"*"``).  Say so wherever one of these numbers is
quoted.
"""
from __future__ import annotations

import argparse
import itertools
import json
import os
import sys
import time
import warnings

import numpy as np

REF_DATA = os.environ.get(
    "POLYFIND_FIELD_NEIGHBORHOOD",
    "E:/source/gh/sarco/materials/gpu_bundle/results/field_neighborhood_refined",
)


# ==================================================================== step 1: the ranking
SYSTEMS = ("pvdf", "vdcn", "an", "cnepo")
TAGS = ("strain+0.000_field+0", "strain+0.020_field+0")


def _flux_module():
    """``examples/fit_charge_flux.py``, imported as a module (it is the fit this audits)."""
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import fit_charge_flux

    return fit_charge_flux


def spearman(a, b) -> float:
    a, b = np.asarray(a, float), np.asarray(b, float)
    ra = np.argsort(np.argsort(a)).astype(float)
    rb = np.argsort(np.argsort(b)).astype(float)
    ra -= ra.mean()
    rb -= rb.mean()
    return float(ra @ rb / np.sqrt((ra @ ra) * (rb @ rb)))


def ranking_audit(root: str = REF_DATA) -> dict:
    """Leave-one-chemistry-out: does the predicted ordering of ``dmu/d(strain)`` survive?

    The fit targets the *residual* of the fixed-increment model, so the prediction for a
    held-out system is ``base + A @ k`` with ``k`` fitted on the other three.  Three fit
    families are audited, and so is the no-flux control.
    """
    from polyfind.forcefield import SimpleFF
    import polyfind.fitting  # noqa: F401 - registers the presets

    fcf = _flux_module()
    base = SimpleFF.from_preset("pvdf-dft-valence")
    data = fcf.load(root)
    P = fcf.PAIRS
    t = fcf.target(data)
    b0, _ = fcf.columns(data, base, [])
    r = t - b0

    def mags(v):
        return np.array([float(np.linalg.norm(v[3 * n:3 * n + 3])) for n in range(4)])

    ref_mag, fix_mag = mags(t), mags(b0)
    out = {"reference": {s: list(np.round(t[3 * n:3 * n + 3], 4)) for n, s in enumerate(SYSTEMS)},
           "reference_magnitude": dict(zip(SYSTEMS, np.round(ref_mag, 3))),
           "fixed_increment_magnitude": dict(zip(SYSTEMS, np.round(fix_mag, 3))),
           "variants": {}}

    print("=" * 100)
    print("STEP 1  Is the response ranking meaningful?  Leave-one-chemistry-out on the flux fit.")
    print("=" * 100)
    print("Reference dmu/d(axial strain), e.A per unit strain, from field_neighborhood_refined")
    print("(exploratory GFN2-xTB, finite pinned two-chain pair in vacuum -- not a bulk crystal).")
    print(f"  {'system':8s} {'reference vector':30s} {'|ref|':>7s}  {'fixed-increment vector':30s} {'|fix|':>7s}")
    for n, s in enumerate(SYSTEMS):
        sl = slice(3 * n, 3 * n + 3)
        print(f"  {s:8s} {str(np.round(t[sl], 3)):30s} {ref_mag[n]:7.2f}  "
              f"{str(np.round(b0[sl], 3)):30s} {fix_mag[n]:7.2f}")
    ref_order = [SYSTEMS[i] for i in np.argsort(-ref_mag)]
    print(f"\n  reference ordering: {' > '.join(ref_order)}")
    out["reference_order"] = ref_order

    variants = {
        "angle only (C-H, C-F)  [SHIPPED]": [("angle", P[0]), ("angle", P[1])],
        "angle only, four pairs": [("angle", p) for p in P],
        "angle + bond (C-H, C-F)": [("angle", p) for p in P[:2]] + [("bond", p) for p in P[:2]],
    }
    ctl = {"order": [SYSTEMS[i] for i in np.argsort(-fix_mag)],
           "spearman": round(spearman(ref_mag, fix_mag), 4),
           "pearson_components": round(float(np.corrcoef(t, b0)[0, 1]), 4),
           "ratio_to_reference": dict(zip(SYSTEMS, np.round(fix_mag / ref_mag, 3)))}
    print(f"  CONTROL, fixed increments, no flux, nothing fitted to this data:")
    print(f"    ordering {' > '.join(ctl['order'])}   Spearman rho = {ctl['spearman']:+.3f}   "
          f"component-wise Pearson r = {ctl['pearson_components']:+.3f}")
    print(f"    predicted/reference magnitude: "
          + "  ".join(f"{s}={v:.2f}" for s, v in ctl["ratio_to_reference"].items()))
    out["control_fixed_increments"] = ctl

    for name, keys in variants.items():
        _, A = fcf.columns(data, base, keys)
        k_all, _, r2_all = fcf.fit(A, r)
        insample = mags(b0 + A @ k_all)
        pred = np.zeros(12)
        coeffs = {}
        for n, s in enumerate(SYSTEMS):
            m = np.ones(12, dtype=bool)
            m[3 * n:3 * n + 3] = False
            k, _, _ = fcf.fit(A[m], r[m])
            coeffs[s] = [round(float(v), 4) for v in k]
            pred[3 * n:3 * n + 3] = b0[3 * n:3 * n + 3] + A[3 * n:3 * n + 3] @ k
        loo = mags(pred)
        ratio = loo / ref_mag
        rec = {
            "n_params": len(keys),
            "r2_all_four": round(float(r2_all), 4),
            "coefficients_all_four": [round(float(v), 4) for v in k_all],
            "in_sample_order": [SYSTEMS[i] for i in np.argsort(-insample)],
            "in_sample_spearman": round(spearman(ref_mag, insample), 4),
            "loo_order": [SYSTEMS[i] for i in np.argsort(-loo)],
            "loo_spearman": round(spearman(ref_mag, loo), 4),
            "loo_pearson_components": round(float(np.corrcoef(t, pred)[0, 1]), 4),
            "loo_magnitude": dict(zip(SYSTEMS, np.round(loo, 3))),
            "loo_ratio_to_reference": dict(zip(SYSTEMS, np.round(ratio, 3))),
            "loo_ratio_spread": round(float(ratio.max() / ratio.min()), 3),
            "loo_coefficients": coeffs,
        }
        out["variants"][name] = rec
        print(f"\n  --- {name}: {len(keys)} parameters, R^2 = {r2_all:+.3f} on all four")
        print(f"      {'system':8s} {'|ref|':>7s} {'|LOO pred|':>11s} {'pred/ref':>9s}   held-out coefficients")
        for n, s in enumerate(SYSTEMS):
            print(f"      {s:8s} {ref_mag[n]:7.2f} {loo[n]:11.2f} {ratio[n]:9.2f}   "
                  + "  ".join(f"{v:+.4f}" for v in coeffs[s]))
        print(f"      LOO ordering {' > '.join(rec['loo_order'])}   "
              f"Spearman rho = {rec['loo_spearman']:+.3f}   "
              f"(in sample {rec['in_sample_spearman']:+.3f})")
        print(f"      component-wise Pearson r = {rec['loo_pearson_components']:+.3f};  "
              f"scale error spans a factor of {rec['loo_ratio_spread']:.2f} across the four")

    perms = np.array([spearman([3, 2, 1, 0], p) for p in itertools.permutations(range(4))])
    out["null"] = {"p_rho_eq_1": round(float((perms > 0.999).mean()), 4),
                   "p_rho_ge_0.8": round(float((perms >= 0.8).mean()), 4)}
    print(f"\n  Null model: with n = 4 a random ordering gives rho = 1 with probability "
          f"{out['null']['p_rho_eq_1']:.3f} and rho >= 0.8 with probability "
          f"{out['null']['p_rho_ge_0.8']:.3f}.  One perfect ordering over four systems is")
    print("  weak evidence on its own; the useful statistics are the scale-error spread and")
    print("  the fixed-increment control.")
    return out


# ================================================================= step 2/3: the screen
def _presets():
    from polyfind.fitting import FITTED_VALENCE
    from polyfind.forcefield import SimpleFF

    return FITTED_VALENCE, SimpleFF.from_preset("pvdf-dft-valence"), SimpleFF.from_preset("pvdf-dft-valence-flux")


def all_trans_seq(polymer) -> tuple:
    return tuple([0] * polymer.bonds_per_repeat)


def close_torsions(polymer, seq, states, max_dev: float = 20.0, w: float = 2e-3):
    """Torsions near ``seq``'s ideal angles that make a **one-period** repeat close, or ``None``.

    Some chemistries have no conformation that closes at ideal RIS angles at all.  PVDC is
    the measured case: its backbone angles are unequal (123 deg at CH2, 114 at CCl2, both
    crystallographic), so an ideal repeat carries 123 - 114 = 9 deg of curl, every sequence
    is a circular arc rather than a stem, and ``pack.periodic_chain`` refuses all of them.
    The polymer would then drop out of the screen for a reason that is an artifact of
    *ideal torsions* rather than a property of the chain
    (``docs/CHEMISTRY_EXTENSION.md`` phase 1).  What the real chain does instead is deflect
    its torsions a few degrees, so this searches for that: least squares on the repeat's own
    ``rotation_error`` with a weak pull towards the ideal angles, bounded to ``max_dev``.

    The closure condition is one equation in ``len(seq)`` unknowns, so its solutions form a
    manifold and which point of it comes back depends on where the search starts -- which is
    ``states.angles``, and those are the *fitted* model's adapted state angles, not the ideal
    180 / 60.  For PVDC's TG+TG-: from the ideal angles it lands on 174.7 / 54.7 deg with
    c = 4.705 A; from the fitted model's own angles, on 177.5 / 27.7 with c = 4.76 after
    refinement; and ``tests/test_pvdc.py``, from the ideal angles with a different
    regularisation, on 175.3 / 49.4 with c = 4.677 -- which is the published 175 deg, 49 deg
    and 4.68 A (Takahagi et al. 1988).  So the **fibre repeat** is robust across the manifold
    to about 2 % and the **torsion pair** is not; only quote the latter from the ideal-angle
    start.  ``refine_crystal`` relaxes the torsions properly afterwards either way.
    """
    from scipy.optimize import least_squares

    from polyfind.pack import periodic_chain_from_torsions

    x0 = np.array([states.angles[s] for s in seq], dtype=float)
    name = "".join(states.names[s] for s in seq)

    def resid(x):
        try:
            ch = periodic_chain_from_torsions(polymer, name, x)
        except (ValueError, RuntimeError):
            return np.concatenate([[1e3], w * (x - x0)])
        return np.concatenate([[ch.rotation_error], w * (x - x0)])

    sol = least_squares(resid, x0, bounds=(x0 - max_dev, x0 + max_dev),
                        xtol=1e-13, ftol=1e-13, gtol=1e-13)
    try:
        ch = periodic_chain_from_torsions(polymer, name, sol.x)
    except (ValueError, RuntimeError):
        return None
    return (sol.x, ch) if ch.rotation_error < 0.02 else None


MAX_ATOMS_PER_CHAIN = 48  # PipelineConfig's own rule: packing cost goes as atoms^2


def build_repeat(polymer, seq, states, helix=None, allow_deflection: bool = True):
    """``(chain, how)`` for one conformation: the ideal-angle repeat, or a closed deflected one.

    A conformation whose crystallographic repeat needs more than
    :data:`MAX_ATOMS_PER_CHAIN` atoms is refused rather than packed -- the same cut
    ``polyfind.pipeline.PipelineConfig.max_atoms_per_chain`` makes, and for the same reason.
    """
    from polyfind.pack import periodic_chain

    if helix is not None and helix.c is not None and helix.periods_per_repeat:
        n_atoms = polymer.atoms_per_repeat * len(seq) * helix.periods_per_repeat // polymer.bonds_per_repeat
        if n_atoms > MAX_ATOMS_PER_CHAIN:
            return None, (f"crystallographic repeat is {helix.periods_per_repeat} sequence periods "
                          f"({n_atoms} atoms, c = {helix.c:.1f} A) > {MAX_ATOMS_PER_CHAIN}: too large to pack")
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return periodic_chain(polymer, seq, states, helix=helix), "ideal angles"
    except (ValueError, RuntimeError) as e:
        first = str(e)
    if not allow_deflection or len(seq) > 4:
        # The closure search below forces a *one-period* repeat, so it is a repair for the
        # rigid-angle artifact on a low-index conformation that should be a straight stem
        # -- PVDC's all-trans and TG+TG-, whose 9 deg of curl per repeat comes entirely
        # from the unequal backbone angles.  Applied to a long helix it would not deflect
        # that conformation, it would replace it with a different one, so it is not
        # attempted past a four-bond sequence.
        return None, f"no ideal-angle repeat ({first})"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        got = close_torsions(polymer, seq, states)
    if got is None:
        return None, f"no repeat, ideal or deflected ({first})"
    x, ch = got
    return ch, f"closed by deflecting torsions to {np.round(x, 1).tolist()} deg"


def chain_dipole(chain, packer_cls_kw=None) -> tuple:
    """``(transverse, axial)`` repeat dipole per monomer, e.A, with the fitted charges.

    Transverse is the component perpendicular to the chain axis.  It is the one that
    matters: a helix with a screw symmetry cancels its transverse dipole over the
    crystallographic repeat and keeps only the axial one, while a planar zigzag adds it
    coherently -- which is the whole of beta-PVDF's ferroelectricity.
    """
    from polyfind.pack import CrystalPacker

    pk = CrystalPacker(chain, n_chains=2, **(packer_cls_kw or {}))
    q = pk.chain.charges
    mu = q @ pk.chain.coords
    n = pk.chain.n_monomers
    return float(np.linalg.norm(mu[:2]) / n), float(mu[2] / n)


def refined_reference_of(polymer, chain, label, valence, flux, angle_stiffness=None):
    """``mechanics.refined_reference`` for an already-built chain (which may be a deflected one)."""
    from polyfind import mechanics as M
    from polyfind.linegroup import repeat_chains
    from polyfind.pack import pack
    from polyfind.refine import refine_crystal

    res = pack(chain, n_chains=2, table_cache_dir=None)[0]
    kw = {"valence": valence}
    if angle_stiffness is not None:
        kw["angle_stiffness"] = angle_stiffness
    rr = refine_crystal(polymer, res, **kw)
    angles0 = np.array([polymer.backbone[k].backbone_angle for k in range(polymer.bonds_per_repeat)])
    frame = repeat_chains(polymer, res.chain, np.asarray(res.dihedrals, float)[None], angles0)[0]
    built = repeat_chains(polymer, res.chain, rr.torsions[None], rr.angles, align_to=frame)[0]
    ref = M.reference_from_chain(built, rr.result.params, n_chains=2, label=label,
                                 valence=valence, charge_flux=flux)
    return ref, rr, res


def film_coefficients(resp):
    """``(d33, d31, e_polar_zz)`` in the film convention, along the crystal's *own* polar axis.

    The PVDF literature quotes a poled film's ``d_33`` (field and strain both along the
    poling direction) and ``d_31`` (field along the poling direction, strain along the draw
    direction, which is the chain axis).  ``examples/fit_charge_flux.py`` hard-codes the
    poling axis as the packer's ``x``, which is right for beta-PVDF and is *not* something
    to assume for a chemistry whose packing has not been looked at: the polar axis is
    wherever the packing puts it, anywhere in the plane perpendicular to the chains.

    So the poling direction here is the reference polarization itself, ``n = P/|P|``, which
    lies in the ``xy`` plane; the film's 3-axis is ``n`` and its 1-axis is ``z``.  Writing
    ``n = (cos t, sin t, 0)``, a strain along ``n`` is
    ``cos^2 t eps_1 + sin^2 t eps_2 + cos t sin t eps_6`` in engineering components, and the
    field projection onto the ``d`` rows is ``n`` itself.  For beta-PVDF, where ``P`` is
    along ``-x``, this reduces to ``-d[x, xx]`` and ``-d[x, zz]`` -- the same numbers
    ``fit_charge_flux`` prints, which is the check that the generalisation is one.

    Returns ``(None, None, None)`` for a crystal with no polarization to pole along.
    """
    el, pz = resp.elastic, resp.piezo
    P = np.asarray(resp.polarization, dtype=float)
    if float(np.linalg.norm(P[:2])) < 1e-6:
        return None, None, None
    n = np.array([P[0], P[1], 0.0]) / float(np.linalg.norm(P[:2]))
    cols = list(el.reachable)
    if 2 not in cols:
        return None, None, None
    ct, st = float(n[0]), float(n[1])
    w = np.zeros(len(cols))          # the strain-along-n combination, in the reachable order
    w[cols.index(0)] = ct * ct
    w[cols.index(1)] = st * st
    if 5 in cols:
        w[cols.index(5)] = ct * st
    d = pz.d_from_e
    d33 = float(n @ d @ w)
    d31 = float(n @ d[:, cols.index(2)])
    e_pz = float(n @ pz.e[:, cols.index(2)])
    return d33, d31, e_pz


def response_of(polymer, chain, label, valence, flux):
    """Pack, refine and take the electromechanical response, deformable where that is possible.

    The deformable path needs a :class:`~polyfind.mechanics.Shape`, which needs a line
    group, and a line group needs the repeat's torsions to follow the sequence's own
    symmetry pattern exactly.  A chain closed by :func:`close_torsions` does not: its
    deflections break the glide it was named for.  Such a repeat therefore falls back to the
    **rigid** path, where ``eps_zz`` is not a variable at all -- ``C_33`` is ``nan``, every
    diagonal column of ``e`` is exactly zero, and the only piezoelectric channel left is
    chain reorientation.  That is a much weaker result and it is labelled as one rather than
    printed alongside the deformable rows as if it were the same measurement.
    """
    from polyfind import mechanics as M
    from polyfind.linegroup import LineGroupError

    ref, rr, res = refined_reference_of(polymer, chain, label, valence, flux)
    try:
        shape = M.shape_of(polymer, ref, angles=rr.angles)
    except LineGroupError as e:
        resp = M.electromechanical_response(ref, shape=None)
        return resp, rr, res, f"rigid chain: no line group for this repeat ({e})"
    resp = M.electromechanical_response(ref, shape=shape)
    return resp, rr, res, None


# The model's own resolution on a per-monomer lattice energy difference, which every
# polarity verdict in this script is measured against.
#
# ``pvdf-dft-valence``'s held-out energy error is 1.359 kcal/mol RMS over ten chemistries
# the objective never saw (``docs/VALENCE_FIT.md``; re-measured from the shipped preset, it
# is 1.3593).  That is a *per-frame* number and the frames are ten-bond oligomers -- eleven
# backbone carbons, five monomers -- so the comparable quantity for a per-monomer crystal
# energy is 1.359 / 5 = 0.27 kcal/mol per monomer.  A polarity gap smaller than that is not
# a prediction: the model does not resolve it.  Both figures are reported, because the
# per-frame one is the conservative reading and under it nothing here resolves at all.
ERROR_BAR_PER_FRAME = 1.3593      # kcal/mol, held-out energy RMS of pvdf-dft-valence
MONOMERS_PER_FRAME = 5            # ten-bond oligomers, eleven backbone carbons
ERROR_BAR = ERROR_BAR_PER_FRAME / MONOMERS_PER_FRAME  # 0.272 kcal/mol per monomer


def verdict(gap: float, error_bar: float = ERROR_BAR) -> str:
    """``'antipolar'``, ``'polar'`` or ``'not resolved'`` for one polarity gap.

    The gap is ``E(best antipolar cell) - E(best cell overall)`` per monomer, so a negative
    gap means the crystal prefers antipolar.  Inside +/- ``error_bar`` the model is not
    deciding anything and saying which side of zero it fell on would be reporting noise.
    """
    if abs(float(gap)) <= error_bar:
        return "not resolved"
    return "antipolar" if gap < 0.0 else "polar"


def _packers(chain, coulomb: str, boundary: str):
    """``(packer, screen_packer, ewald_spec)`` for the polar/antipolar comparison.

    **Resolved through the module, and that is load-bearing.**
    :meth:`polyfind.fitting.FFParameters.applied` forces the fitted potential by *rebinding
    ``polyfind.pack.CrystalPacker``* for the duration of its block, so a name imported with
    ``from polyfind.pack import CrystalPacker`` before the block -- or at module scope --
    still points at the original class and builds a packer carrying the **illustrative**
    potential.  This script used to do exactly that: its polar branch went through
    ``pack()``, which resolves the class through the module and so was fitted, while its
    antipolar branch was built from the module-level name and so was not.  **Every antipolar
    gap it printed was a difference between two different potentials** on top of being a
    difference between two polar cells (``docs/SCREEN.md``).  ``fitting.predict`` has the
    same note and does it the same way.

    An Ewald packer screens with its own truncated twin, which is the division of labour
    ``pack(coulomb="ewald")`` already draws and for the same reason: the reciprocal half of
    the sum is not pairwise, so it cannot be tabulated or cheaply gridded, and the screen
    only chooses starts while the polish is exact.
    """
    import polyfind.pack as pack_mod
    from polyfind.ewald import EwaldSpec

    if coulomb == "dsf":
        return pack_mod.CrystalPacker(chain, n_chains=2), None, None
    spec = EwaldSpec(boundary=boundary)
    return (pack_mod.CrystalPacker(chain, n_chains=2, coulomb="ewald", ewald=spec),
            pack_mod.CrystalPacker(chain, n_chains=2), spec)


def screen_one(polymer, max_period: int = 8, k_per_period: int = 40, top_pack: int = 3,
               fit_step: float = 10.0, fit_monomers: int = 6, mech: int = 1,
               skip_fit_reason: str | None = None, coulomb: str = "dsf",
               boundary: str = "tinfoil", anti_polish: int = 8, anti_maxfev: int = 900,
               anti_target: int = 6000, error_bar: float = ERROR_BAR) -> dict:
    """The whole funnel for one chemistry.  Returns a record; prints as it goes."""
    from polyfind.enumerate import KNOWN_CHAINS, enumerate_periodic
    from polyfind.forcefield import fit_ris
    from polyfind.helix import canonical_sequence
    from polyfind.polymers import THREE_STATE

    FITTED_VALENCE, VAL, FLUX = _presets()
    rec = {"polymer": polymer.name, "formula": polymer.formula, "chiral": bool(polymer.is_chiral),
           "atoms_per_repeat": polymer.atoms_per_repeat, "timings": {}, "notes": []}
    print("\n" + "=" * 100)
    print(f"{polymer.name.upper()}  {polymer.formula}   "
          f"({'chiral, isotactic' if polymer.is_chiral else 'achiral'}, "
          f"{polymer.atoms_per_repeat} atoms per repeat)")
    print("=" * 100)
    t_chem = time.time()

    if skip_fit_reason:
        rec["screened"] = False
        rec["reason"] = skip_fit_reason
        print(f"  NOT SCREENED: {skip_fit_reason}")
        rec["timings"]["total"] = round(time.time() - t_chem, 2)
        return rec

    # --- 1. the RIS fit ------------------------------------------------------------------
    t0 = time.time()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        fit = fit_ris(polymer, VAL, step=fit_step, n_monomers=fit_monomers, third_order=True)
    for w in caught:
        msg = str(w.message).split("\n")[0]
        if msg not in rec["notes"]:
            rec["notes"].append(msg)
    model = fit.model
    rec["timings"]["fit"] = round(time.time() - t0, 2)
    rec["fit_evaluations"] = int(fit.n_evaluations)
    # which scan the chemistry got: rigid (the frozen-angle scan) or relaxed (every backbone
    # angle relaxed per conformer, the default for chemistries whose rigid profile has wells
    # the relaxed one lacks -- docs/NITRILE_LANDSCAPE.md, forcefield.ANGLE_RELAXATION_DEFAULTS)
    rec["fit_angles"] = fit.angles
    rec["fit_relaxation_energies"] = int(fit.n_relaxation_energies)
    print(f"  [1] RIS fit ({fit.angles} backbone angles): {fit.n_evaluations} conformers"
          + (f", {fit.n_relaxation_energies} energies in relaxation" if fit.angles == "relaxed" else "")
          + f", {rec['timings']['fit']:.1f} s"
          + (f"   (symmetrised over reflection+reversal: chiral)" if polymer.is_chiral else ""))

    # --- 2. exhaustive conformational enumeration -----------------------------------------
    t0 = time.time()
    cands = enumerate_periodic(polymer, model, max_period=max_period, k_per_period=k_per_period)
    rec["timings"]["enumerate"] = round(time.time() - t0, 2)
    rec["n_conformations"] = len(cands)
    e0 = cands[0].energy_per_monomer
    print(f"  [2] {len(cands)} distinct periodic conformations to period {max_period} "
          f"({rec['timings']['enumerate']:.1f} s).  Ground state {cands[0].name} at "
          f"{e0:+.3f} kcal/mol per monomer.")
    rec["conformations"] = [
        {"rank": c.rank, "name": c.name, "period": c.period,
         "e_per_monomer": round(c.energy_per_monomer, 4),
         "de_above_ground": round(c.energy_per_monomer - e0, 4),
         "c": None if c.helix.c is None else round(float(c.helix.c), 3),
         "label": c.helix.label, "known_as": c.known_as}
        for c in cands[:12]]
    for c in cands[:6]:
        print("      " + c.row())

    # --- 3. which conformations to carry forward ------------------------------------------
    # The top few by RIS energy, plus the three conformations that are *structurally*
    # decisive for this application whatever they rank: all-trans (the polar zigzag that
    # beta-PVDF's ferroelectricity is), TG+TG- (its antipolar competitor, alpha) and
    # T3GT3G' (gamma, polar but weaker).  Screening only the top of the RIS list would let a
    # polymer look useless because its polar phase happened to rank fifth.
    wanted, seen = [], set()

    def add(cand, why):
        if cand is None or cand.seq in seen:
            return
        seen.add(cand.seq)
        wanted.append((cand, why))

    def by_name(s):
        """The enumerated candidate for ``s``, or a synthetic one when it was filtered out.

        ``enumerate_periodic`` drops any sequence whose rise per bond is below 0.3 A,
        because such a chain does not propagate as a stem.  For PVDC that filter removes
        *all-trans and TG+TG- both*: with the measured 123 / 114 deg backbone angles an
        ideal PVDC repeat is a circular arc with zero rise.  Dropping them is right for the
        RIS ranking and wrong for this screen, because the deflected TG+TG- glide is the
        structure PVDC actually has.  So a structurally decisive sequence that the
        enumeration filtered out is rebuilt here, with its own RIS energy, and handed to
        ``build_repeat``, which will find the deflection that closes it or say it cannot.
        """
        from polyfind.enumerate import Candidate
        from polyfind.helix import helix_parameters

        try:
            key = canonical_sequence(model.parse(s), model.states, polymer.bonds_per_repeat,
                                     chiral=bool(polymer.is_chiral))
        except (ValueError, KeyError):
            return None
        hit = next((c for c in cands if c.seq == key), None)
        if hit is not None:
            return hit
        h = helix_parameters(polymer, np.array(key), model.states, max_m=8)
        e = float(model.energy(list(key), periodic=True))
        return Candidate(seq=key, name="".join(model.states.names[x] for x in key),
                         period=len(key), energy_per_period=e,
                         energy_per_monomer=e / (len(key) / polymer.bonds_per_repeat),
                         helix=h, rank=-1)

    # the all-trans conformation, wherever it sits: the single most decision-relevant row
    at_cand = by_name("T" * polymer.bonds_per_repeat)
    if at_cand.rank < 0:
        rec["all_trans"] = {"in_list": False,
                            "de_above_ground": round(at_cand.energy_per_monomer - e0, 4)}
        print(f"      all-trans ({at_cand.name}): NOT among the {len(cands)} enumerated "
              f"conformations -- it does not propagate as a stem (rise per bond below 0.3 A). "
              f"Its RIS energy is {at_cand.energy_per_monomer - e0:+.3f} kcal/mol per monomer "
              f"above the ground state.")
    else:
        rec["all_trans"] = {"in_list": True, "rank": at_cand.rank,
                            "de_above_ground": round(at_cand.energy_per_monomer - e0, 4),
                            "c": None if at_cand.helix.c is None else round(float(at_cand.helix.c), 3)}
        print(f"      all-trans ({at_cand.name}): rank {at_cand.rank} of {len(cands)}, "
              f"{at_cand.energy_per_monomer - e0:+.3f} kcal/mol per monomer above the ground state")
    add(at_cand, "all-trans: the polar zigzag")
    for c in cands[:top_pack]:
        add(c, f"RIS rank {c.rank}")
    if polymer.bonds_per_repeat == 2:
        add(by_name("TG+TG-"), "TG+TG-: the alpha-type competitor")
        add(by_name("TTTG+TTTG-"), "T3GT3G': the gamma-type competitor")
    for lbl, s in KNOWN_CHAINS.get(polymer.name, {}).items():
        add(by_name(s), f"known: {lbl}")

    # --- 4. packing -----------------------------------------------------------------------
    # The protocol is :func:`polyfind.fitting.predict`'s, because that is the one the
    # alpha-below-beta validation in ``docs/BENCHMARK.md`` was measured with: exhaustive
    # table screen, exact polish, the *symmetric antipolar branch* searched explicitly (a
    # screen that only ever proposes polar starts would otherwise miss it), then continuous
    # refinement of torsions, backbone angles and cell, with the bond-angle strain added
    # back into the reported energy.
    from polyfind.fitting import antipolar_cell_exact
    from polyfind.pack import pack as _pack
    from polyfind.refine import refine_crystal

    t0 = time.time()
    packed = []
    rec["coulomb"] = coulomb if coulomb == "dsf" else f"ewald[{boundary}]"
    rec["error_bar_per_monomer"] = round(error_bar, 4)
    print(f"  [3] packing {len(wanted)} conformations (exhaustive screen, antipolar branch, "
          f"refinement; 2 chains per cell; polarity under "
          f"{'the truncated 8 A DSF sum' if coulomb == 'dsf' else f'Ewald, {boundary}'}):")
    with FITTED_VALENCE.applied():
        for cand, why in wanted:
            chain, how = build_repeat(polymer, cand.seq, model.states, helix=cand.helix)
            if chain is None:
                packed.append({"name": cand.name, "why": why, "packed": False, "reason": how})
                print(f"      {cand.name:<16s} NOT PACKED: {how}   [{why}]")
                continue
            try:
                packer, screen_pk, spec = _packers(chain, coulomb, boundary)
                # Both branches of the comparison use the same sum, and an Ewald run screens
                # both with the truncated twin and polishes both with Ewald -- otherwise the
                # gap would be a difference of two different summations rather than of two
                # cells.
                rigid = _pack(chain, n_chains=2, table_cache_dir=None,
                              coulomb=coulomb, ewald=spec)[0]
                # ``antipolar_cell_exact``, not ``antipolar_cell``: the latter's
                # "flip and equal setting angles" is the antipolar subspace only when the
                # chain's transverse moment is perpendicular to its own x, which is false
                # for every planar zigzag, so it compared two *polar* cells for every
                # all-trans row this script printed before docs/SCREEN.md's correction.
                # Its axial scan is now a length rather than four points, which matters here:
                # a four-point dz scan of a four-monomer gamma-type repeat lands one sample
                # per monomer and is aliased to the registry it is supposed to be searching.
                t_anti = time.time()
                anti_params, anti_e, anti_pol = antipolar_cell_exact(
                    packer, screen_packer=screen_pk, target=anti_target,
                    n_polish=anti_polish, maxfev=anti_maxfev)
                assert anti_pol < 1e-9, f"antipolar branch is not antipolar: |P| = {anti_pol:.2e}"
                t_anti = time.time() - t_anti
                polar_gap = anti_e - rigid.energy_per_monomer
                if polar_gap < 0.0:
                    rigid = packer.result(anti_params)
                rr = refine_crystal(polymer, rigid)
                res = rr.result
                e_mon = res.energy_per_monomer + rr.angle_energy
            except (AssertionError, ValueError, RuntimeError, KeyError, np.linalg.LinAlgError) as e:
                packed.append({"name": cand.name, "why": why, "packed": False, "reason": str(e)})
                print(f"      {cand.name:<16s} NOT PACKED: {e}")
                continue
            perp, axial = chain_dipole(chain)
            row = {"name": cand.name, "why": why, "packed": True, "build": how,
                   "rank": cand.rank, "de_above_ground": round(cand.energy_per_monomer - e0, 4),
                   "a": round(res.a, 3), "b": round(res.b, 3), "c": round(res.c, 4),
                   "gamma": round(res.gamma, 2), "density": round(res.density, 4),
                   # Judged on the *rigid* comparison, which is the one that means something:
                   # ``antipolar_cell_exact`` measures the symmetric antipolar subspace explicitly,
                   # while a refinement started there slides back out of it (the subspace is
                   # not stationary), so the refined cell's own flip flag is not the answer.
                   # And it is judged against the model's own resolution: inside +/- ERROR_BAR
                   # the sign of the gap is not a prediction.
                   "arrangement": verdict(polar_gap, error_bar),
                   "arrangement_sign": "antipolar" if polar_gap < 0.0 else "polar",
                   "antipolar_gap": round(float(polar_gap), 4),
                   "antipolar_e_per_monomer": round(float(anti_e), 4),
                   "polar_e_per_monomer": round(float(rigid.energy_per_monomer), 4),
                   "antipolar_max_P": float(f"{anti_pol:.3e}"),
                   "antipolar_seconds": round(t_anti, 1),
                   "flip": int(res.flip),
                   "e_per_monomer": round(e_mon, 4),
                   "angle_strain": round(float(rr.angle_energy), 4),
                   "P": round(res.polarization_magnitude, 4),
                   "chain_dipole_perp_per_monomer": round(perp, 4),
                   "chain_dipole_axial_per_monomer": round(axial, 4)}
            packed.append(row)
            print(f"      {cand.name:<16s} {res.row()}   [{why}; {how}]")
            print(f"      {'':<16s} E/mon with angle strain {e_mon:+8.3f}; antipolar minus polar "
                  f"{polar_gap:+.3f} +/- {error_bar:.2f} kcal/mol per monomer -> "
                  f"{row['arrangement'].upper()} ({t_anti:.0f} s, max|P| {anti_pol:.1e}); "
                  f"chain dipole per monomer: transverse {perp:.4f}, axial {axial:+.4f} e.A")
    rec["timings"]["pack"] = round(time.time() - t0, 2)
    rec["packings"] = packed
    ok = [p for p in packed if p["packed"]]
    if ok:
        best = min(ok, key=lambda p: p["e_per_monomer"])
        rec["lattice_ranking"] = [
            {"name": p["name"], "de_per_monomer": round(p["e_per_monomer"] - best["e_per_monomer"], 4),
             "arrangement": p["arrangement"]}
            for p in sorted(ok, key=lambda p: p["e_per_monomer"])]
        print(f"      lattice energy ranking (kcal/mol per monomer, refined, angle strain "
              f"included, relative to the best):")
        for p in rec["lattice_ranking"]:
            print(f"        {p['name']:<16s} {p['de_per_monomer']:+7.3f}   {p['arrangement']}")

    # --- 5. the electromechanical response ------------------------------------------------
    t0 = time.time()
    rec["responses"] = []
    # The polar phase first -- it is the one the application needs, whether or not it is
    # the lattice ground state -- then the ground state if that is a different structure.
    # ``response_of`` re-packs without the antipolar constraint, so the polar phase is the
    # structure that is measured even where the antipolar branch is lower.
    targets = []
    polar = [p for p in ok if p["P"] > 1e-2]
    if polar:
        targets.append(max(polar, key=lambda p: p["P"]))
    if ok:
        b = min(ok, key=lambda p: p["e_per_monomer"])
        if not any(t["name"] == b["name"] for t in targets):
            targets.append(b)
    # ``mech`` defaults to 1 -- the polar phase and nothing else -- because the second
    # target is where the cost blows up: a helical phase on the flux path is exactly the
    # case ``docs/BENCHMARK.md`` records as not converging (alpha-PVDF's relaxation runs to
    # the line group's +/-8 deg cap and the two routes to d then disagree completely), and
    # it costs many minutes to arrive at a row that has to be reported as a
    # non-measurement.  ``--mech 2`` asks for it anyway.
    targets = targets[:mech]
    print(f"  [4] electromechanical response on {len(targets)} packing(s) "
          f"(deformable path: valence terms in the kernel, charge flux in the mechanics packer):")
    sys.stdout.flush()
    for tgt in targets:
        cand = next(c for c, _ in wanted if c.name == tgt["name"])
        with FITTED_VALENCE.applied():
            chain, _ = build_repeat(polymer, cand.seq, model.states, helix=cand.helix)
            try:
                resp, rr, res, degraded = response_of(polymer, chain, tgt["name"], VAL, FLUX)
            except (ValueError, RuntimeError, np.linalg.LinAlgError) as e:
                rec["responses"].append({"name": tgt["name"], "ok": False, "reason": str(e)})
                print(f"      {tgt['name']}: FAILED ({e})")
                continue
        el, pz, ac = resp.elastic, resp.piezo, resp.actuator
        d33, d31, e_pz = film_coefficients(resp)
        row = {
            "name": tgt["name"], "ok": True, "path": "rigid" if degraded else "deformable",
            "degraded": degraded,
            "P_magnitude": round(float(np.linalg.norm(resp.polarization)), 4),
            "P_vector": [round(float(v), 4) for v in resp.polarization],
            "cell": [round(float(resp.reference.params[0]), 3), round(float(resp.reference.params[1]), 3),
                     round(float(resp.reference.params[2]), 2), round(float(resp.reference.c), 4)],
            "volume": round(resp.volume, 2),
            "C11": round(float(el.C[0, 0]), 2), "C22": round(float(el.C[1, 1]), 2),
            "C12": round(float(el.C[0, 1]), 2), "C66": round(float(el.C[5, 5]), 2),
            "C33": round(float(el.C[2, 2]), 1),
            "C33_energy_route": round(float(el.c33_from_energy), 1),
            "residual_axial_stress": float(f"{el.residual_stress[2]:.3e}"),
            "d33_film": None if d33 is None else round(d33, 3),
            "d31_film": None if d31 is None else round(d31, 3),
            "e_polar_zz": None if e_pz is None else round(e_pz, 5),
            "reachable_voigt": [int(k) + 1 for k in el.reachable],
            "e_matrix": [[round(float(v), 5) for v in r] for r in pz.e],
            "d_from_e": [[round(float(v), 4) for v in r] for r in pz.d_from_e],
            "d_direct": [[round(float(v), 4) for v in r] for r in pz.d_direct],
            "routes_agree_pct": round(float(pz.relative_difference * 100), 3),
            "actuator_direction": [round(float(v), 3) for v in ac.direction],
            "free_strain_max": float(f"{np.abs(ac.free_strain).max():.3e}"),
            "blocking_stress_max": round(float(np.abs(ac.blocking_stress).max()), 5),
            "work_density": round(float(ac.work_density), 4),
            "seconds": round(resp.seconds, 1),
        }
        rec["responses"].append(row)
        print(f"      --- {tgt['name']}"
              + (f"   [DEGRADED to the {row['path']} path: {degraded}]" if degraded else ""))
        print("      " + resp.summary().replace("\n", "\n      "))
        # Two checks that decide whether the row is a measurement at all, both of them the
        # ones ``docs/BENCHMARK.md`` applies to alpha-PVDF's flux run: the reference must be
        # axially stress-free (otherwise the line group's +/-8 deg cap has bound and the
        # relaxation stopped at a wall, not a minimum), and the dipole-derivative route and
        # the zero-stress root find must agree (they are independent, so a disagreement
        # means the flux did not reach one of them consistently).
        bad = []
        if degraded is None and abs(el.residual_stress[2]) > 1e-3:
            bad.append(f"residual axial stress {el.residual_stress[2]:+.3e} GPa: the reference is "
                       f"not stress-free (the +/-8 deg angle cap has bound)")
        if pz.relative_difference > 0.05:
            bad.append(f"the two routes to d disagree by {pz.relative_difference * 100:.0f}% "
                       f"(max |dd| = {pz.max_abs_difference:.3g} pC/N)")
        if bad:
            row["ok"] = False
            row["reason"] = "; ".join(bad)
            print(f"      NOT A MEASUREMENT: {row['reason']}")
    rec["timings"]["response"] = round(time.time() - t0, 2)
    rec["timings"]["total"] = round(time.time() - t_chem, 2)
    rec["screened"] = True
    print(f"  wall time for {polymer.name}: {rec['timings']['total']:.1f} s "
          f"(fit {rec['timings']['fit']:.1f}, enumerate {rec['timings']['enumerate']:.1f}, "
          f"pack {rec['timings']['pack']:.1f}, response {rec['timings']['response']:.1f})")
    return rec


# ------------------------------------------------------------------- the candidate set
def incumbents():
    """The campaign's chemistries, with the ones that cannot be screened marked and why."""
    from polyfind.polymers import get_polymer

    out = []
    for n in ("pvdf", "pvdc", "cfe", "cdfe", "an", "vdcn"):
        out.append((get_polymer(n), None))
    out.append((get_polymer("fanome"),
                "FANOME's all-trans reference is an overlapping structure (methyl hydrogens of "
                "methoxy groups on consecutive substituted carbons 0.80 A apart, ~1e6 kcal/mol on a "
                "ten-bond oligomer, robust to the frozen C-O rotamer), so fit_ris would measure "
                "every RIS energy from a state that is not a molecule.  It is registered and its "
                "geometry is tested, and deliberately not fitted "
                "(docs/CHEMISTRY_EXTENSION.md phase 3)."))
    return out


def new_candidates():
    """Substituent combinations expressible in the current model, built here rather than registered.

    These are **model suggestions, not predictions**.  They are ordinary
    :class:`polyfind.polymers.Polymer` values constructed in this script and passed straight
    to the funnel, so nothing in the package changes and no default moves: geometry and
    charges are chosen exactly the way PVDF's, PVDC's, CFE's and CDFE's were (textbook bond
    lengths, an inter-pendant angle interpolated between the parents', backbone angles at
    114 deg, modest point charges leaving every backbone group neutral).  The charges are
    *illustrative*, so every dipole below inherits that.
    """
    from polyfind.polymers import BackboneAtom, Polymer, nitrile

    CN = nitrile()
    CH2 = BackboneAtom("C", "H", 1.09, 114.0, 108.0, -0.20, +0.10)
    return [
        # Poly(vinyl fluoride): the smallest polar repeat this model can express.  One
        # fluorine per two backbone carbons, so half PVDF's dipole per monomer in a repeat
        # of nearly the same length -- the test of whether dipole *density* or dipole per
        # monomer is what the work density follows.  Chiral (a stereocentre at the CHF).
        Polymer(name="pvf-cand", formula="-(CH2-CHF)n-", bond_length=1.528,
                backbone=(CH2, BackboneAtom("C", ("H", "F"), (1.09, 1.35), 114.0, 107.0,
                                            +0.10, (+0.10, -0.20)))),
        # Trifluoroethylene, the VDF copolymer partner that is known to stabilise the polar
        # all-trans phase.  Included as a *positive control* on the design rules as much as a
        # candidate: if the rules are worth anything they should like it.
        Polymer(name="trfe-cand", formula="-(CHF-CF2)n-", bond_length=1.528,
                backbone=(BackboneAtom("C", ("H", "F"), (1.09, 1.35), 114.0, 107.0,
                                       +0.10, (+0.10, -0.20)),
                          BackboneAtom("C", "F", 1.35, 114.0, 106.0, +0.40, -0.20))),
        # A fluorine and a nitrile on the same carbon: VDCN's dipole strength with one
        # nitrile's worth of the steric bulk that costs VDCN its all-trans phase.
        Polymer(name="vfcn-cand", formula="-(CH2-C(F)(CN))n-", bond_length=1.54,
                backbone=(CH2, BackboneAtom("C", ("F", CN), (1.35, None), 114.0, 108.0,
                                            +0.35, (-0.20, None)))),
        # A chlorine and a nitrile: the same idea with a bulkier, less electronegative
        # halogen, to separate "big dipole" from "small pendant" in the rules.
        Polymer(name="vclcn-cand", formula="-(CH2-C(Cl)(CN))n-", bond_length=1.54,
                backbone=(CH2, BackboneAtom("C", ("Cl", CN), (1.77, None), 114.0, 109.0,
                                            +0.25, (-0.10, None)))),
    ]


# ------------------------------------------------------------------------- design rules
def design_rules(records: list) -> dict:
    """Correlate the structural columns against the predicted work density.

    Deliberately blunt about ``n``: with a handful of chemistries a correlation coefficient
    is a description of the sample, not evidence.  Anything reported here is a *rule the
    data is consistent with*, and the report says which ones the data cannot support.
    """
    rows = []
    for r in records:
        if not r.get("screened"):
            continue
        good = [x for x in r.get("responses", []) if x.get("ok")]
        if not good:
            continue
        best = max(good, key=lambda x: x["work_density"])
        pk = next((p for p in r["packings"] if p["packed"] and p["name"] == best["name"]), None)
        if pk is None:
            continue
        gap = None
        lat = r.get("lattice_ranking") or []
        if len(lat) > 1:
            here = next((i for i, p in enumerate(lat) if p["name"] == best["name"]), None)
            if here is not None:
                others = [abs(p["de_per_monomer"] - lat[here]["de_per_monomer"])
                          for i, p in enumerate(lat) if i != here]
                gap = round(min(others), 4) if others else None
        rows.append({
            "polymer": r["polymer"], "phase": best["name"],
            "work_density": best["work_density"],
            # |P| is already a dipole per unit volume (C/m^2), which is the feature the
            # question asks about; the per-monomer transverse chain dipole is the same
            # quantity before the packing decides how densely the chains sit.
            "P": best["P_magnitude"],
            "chain_dipole_perp": pk["chain_dipole_perp_per_monomer"],
            "density": pk["density"], "C33": best["C33"], "C11": best["C11"],
            "d31": best["d31_film"], "d33": best["d33_film"],
            "de_all_trans_above_ground": (r.get("all_trans") or {}).get("de_above_ground"),
            "polymorph_gap": gap,
            "arrangement": pk["arrangement"],
        })
    if not rows:
        return {"rows": [], "correlations": {}}
    keys = ("P", "chain_dipole_perp", "density", "C33", "C11",
            "de_all_trans_above_ground", "polymorph_gap")
    w = np.array([r["work_density"] for r in rows], float)
    corr = {}
    for k in keys:
        v = np.array([np.nan if r[k] is None else r[k] for r in rows], float)
        m = np.isfinite(v) & np.isfinite(w)
        if m.sum() >= 3 and np.ptp(v[m]) > 0 and np.ptp(w[m]) > 0:
            corr[k] = {"n": int(m.sum()),
                       "pearson": round(float(np.corrcoef(v[m], w[m])[0, 1]), 3),
                       "spearman": round(spearman(v[m], w[m]), 3)}
        else:
            corr[k] = {"n": int(m.sum()), "pearson": None, "spearman": None}
    return {"rows": rows, "correlations": corr}


def print_rules(rules: dict) -> None:
    print("\n" + "=" * 100)
    print("STEP 3  What drives the answer: structural features against the predicted work density")
    print("=" * 100)
    if not rules["rows"]:
        print("  no screened chemistry produced a usable response; no rule can be extracted")
        return
    hdr = (f"  {'polymer':12s} {'phase':<12s} {'w kJ/m3':>8s} {'|P| C/m2':>9s} {'mu_perp':>8s} "
           f"{'rho':>6s} {'C33':>7s} {'C11':>7s} {'d31':>7s} {'dE(TT)':>7s} {'gap':>7s}  arrangement")
    print(hdr)
    # Any of these can legitimately be None -- ``film_coefficients`` returns None for a cell
    # with no polarization to pole along, which is exactly what happens when a polar phase
    # re-packs antipolar, and ``C_33`` is nan on the rigid fallback path.  Formatting them was
    # a latent crash (AN's TT hit it on this re-run once its polarization collapsed).
    def num(v, w, p):
        return (" " * (w - 3) + "n/a") if v is None or not np.isfinite(v) else f"{v:{w}.{p}f}"

    for r in sorted(rules["rows"], key=lambda x: -x["work_density"]):
        print(f"  {r['polymer']:12s} {r['phase']:<12s} {num(r['work_density'], 8, 3)} "
              f"{num(r['P'], 9, 4)} {num(r['chain_dipole_perp'], 8, 4)} {num(r['density'], 6, 3)} "
              f"{num(r['C33'], 7, 1)} {num(r['C11'], 7, 2)} {num(r['d31'], 7, 2)} "
              f"{num(r['de_all_trans_above_ground'], 7, 2)} {num(r['polymorph_gap'], 7, 2)}  "
              f"{r['arrangement']}")
    print(f"\n  correlation of each feature with the work density (n = {len(rules['rows'])}; "
          f"at this n these describe the sample, they do not establish a rule):")
    for k, v in rules["correlations"].items():
        if v["pearson"] is None:
            print(f"    {k:28s} n={v['n']}  not computable")
        else:
            print(f"    {k:28s} n={v['n']}  Pearson {v['pearson']:+.3f}  Spearman {v['spearman']:+.3f}")


def print_polarity(records: list, error_bar: float = ERROR_BAR) -> dict:
    """The polarity column against the model's own error bar, one row per phase.

    This is the table the screen's primary output turns out to be, so it is printed on its
    own rather than buried in the packing rows: a gap, the resolution of the model that
    produced it, and whether the verdict is a prediction or a "cannot tell".
    """
    rows = []
    for r in records:
        for p in r.get("packings", []):
            if p.get("packed"):
                rows.append((r["polymer"], p))
    print("\n" + "=" * 100)
    print("POLARITY against the model's own error bar")
    print("=" * 100)
    print(f"  E(best antipolar cell) - E(best cell), rigid, 2 chains, per monomer.  Negative = the")
    print(f"  crystal prefers antipolar.  Error bar {error_bar:.2f} kcal/mol per monomer: "
          f"{ERROR_BAR_PER_FRAME:.3f} held-out")
    print(f"  energy RMS over ten unseen chemistries, on ten-bond ({MONOMERS_PER_FRAME}-monomer) oligomer frames.")
    if not rows:
        print("  nothing packed")
        return {"rows": [], "n_resolved": 0, "n_total": 0}
    print(f"\n  {'chemistry':12s} {'phase':<14s} {'mu_perp':>8s} {'gap':>8s} {'+/-':>6s} "
          f"{'|P| anti':>9s}  verdict")
    out = []
    for name, p in rows:
        out.append({"polymer": name, "phase": p["name"], "gap": p["antipolar_gap"],
                    "error_bar": round(error_bar, 4), "verdict": p["arrangement"],
                    "sign": p["arrangement_sign"],
                    "mu_perp": p["chain_dipole_perp_per_monomer"]})
        print(f"  {name:12s} {p['name']:<14s} {p['chain_dipole_perp_per_monomer']:8.3f} "
              f"{p['antipolar_gap']:+8.3f} {error_bar:6.2f} {p['antipolar_max_P']:9.1e}  "
              f"{p['arrangement']}"
              + ("" if p["arrangement"] != "not resolved" else f"  (sign alone says {p['arrangement_sign']})"))
    n_res = sum(1 for r in out if r["verdict"] != "not resolved")
    print(f"\n  {n_res} of {len(out)} phases resolved; {len(out) - n_res} inside the error bar.")
    by_chem = {}
    for r in out:
        by_chem.setdefault(r["polymer"], []).append(r)
    n_chem_res = sum(1 for v in by_chem.values() if any(x["verdict"] != "not resolved" for x in v))
    print(f"  {n_chem_res} of {len(by_chem)} chemistries have at least one resolved phase.")
    return {"rows": out, "n_resolved": n_res, "n_total": len(out),
            "n_chemistries_resolved": n_chem_res, "n_chemistries": len(by_chem),
            "error_bar_per_monomer": round(error_bar, 4),
            "error_bar_per_frame": ERROR_BAR_PER_FRAME}


# --------------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--ranking", action="store_true", help="step 1: the leave-one-out ranking audit")
    ap.add_argument("--screen", action="store_true", help="step 2: the funnel on the campaign's chemistries")
    ap.add_argument("--candidates", action="store_true", help="step 3: new substituent combinations")
    ap.add_argument("--all", action="store_true", help="all three")
    ap.add_argument("--only", default=None, help="comma-separated polymer names to screen")
    ap.add_argument("--mech", type=int, default=1,
                    help="how many packings per chemistry get a full response (1 = the polar phase)")
    ap.add_argument("--max-period", type=int, default=8)
    ap.add_argument("--k-per-period", type=int, default=40)
    ap.add_argument("--coulomb", default="dsf", choices=("dsf", "ewald"),
                    help="electrostatic sum for the polar/antipolar comparison (both branches)")
    ap.add_argument("--boundary", default="tinfoil", choices=("tinfoil", "vacuum"),
                    help="Ewald boundary convention; tinfoil is the bulk limit of a screened crystal")
    ap.add_argument("--error-bar", type=float, default=ERROR_BAR,
                    help="kcal/mol per monomer below which a polarity gap is 'not resolved'")
    ap.add_argument("--anti-polish", type=int, default=8, help="antipolar starts polished")
    ap.add_argument("--anti-maxfev", type=int, default=900, help="function evaluations per antipolar polish")
    ap.add_argument("--anti-target", type=int, default=6000, help="antipolar screen budget in cells")
    ap.add_argument("--data", default=REF_DATA, help="path to results/field_neighborhood_refined")
    ap.add_argument("--json", default=None, help="write the whole record here")
    args = ap.parse_args()
    if args.all:
        args.ranking = args.screen = args.candidates = True
    if not (args.ranking or args.screen or args.candidates):
        ap.error("nothing to do: pass --ranking, --screen, --candidates or --all")

    out = {"generated": time.strftime("%Y-%m-%d %H:%M:%S"),
           "preset": "pvdf-dft-valence (+ pvdf-dft-valence-flux in the mechanics packer)",
           "coulomb": args.coulomb if args.coulomb == "dsf" else f"ewald[{args.boundary}]",
           "error_bar_per_monomer": round(args.error_bar, 4),
           "error_bar_per_frame": ERROR_BAR_PER_FRAME}
    kw = dict(max_period=args.max_period, k_per_period=args.k_per_period, mech=args.mech,
              coulomb=args.coulomb, boundary=args.boundary, error_bar=args.error_bar,
              anti_polish=args.anti_polish, anti_maxfev=args.anti_maxfev,
              anti_target=args.anti_target)
    t_all = time.time()

    if args.ranking:
        out["ranking"] = ranking_audit(args.data)

    records = []
    if args.screen:
        print("\n" + "=" * 100)
        print("STEP 2  The funnel on every chemistry the model can express")
        print("=" * 100)
        print("CNEPO is out of scope: its epoxide bridges two backbone carbons, which is a change of")
        print("chain topology rather than a substituent, and BackboneAtom carries pendants, not rings")
        print("(docs/CHEMISTRY_EXTENSION.md phase 4).  It is not approximated by anything else.")
        out["out_of_scope"] = {"cnepo": "ring backbone (epoxide bridging two backbone carbons); "
                                        "not expressible by BackboneAtom, see CHEMISTRY_EXTENSION.md phase 4"}
        only = set(args.only.split(",")) if args.only else None
        for polymer, skip in incumbents():
            if only and polymer.name not in only:
                continue
            records.append(screen_one(polymer, skip_fit_reason=skip, **kw))
        out["incumbents"] = records

    cand_records = []
    if args.candidates:
        print("\n" + "=" * 100)
        print("STEP 3b  New substituent combinations, run through the identical funnel")
        print("=" * 100)
        print("MODEL SUGGESTIONS, NOT PREDICTIONS.  Illustrative charges, PVDF-fitted potential,")
        print("rigid backbone angles, isotactic-only chirality, two chains per cell.")
        only = set(args.only.split(",")) if args.only else None
        for polymer in new_candidates():
            if only and polymer.name not in only:
                continue
            cand_records.append(screen_one(polymer, **kw))
        out["candidates"] = cand_records

    if records or cand_records:
        out["polarity"] = print_polarity(records + cand_records, args.error_bar)
        rules = design_rules(records + cand_records)
        print_rules(rules)
        out["design_rules"] = rules

    out["wall_seconds"] = round(time.time() - t_all, 1)
    print(f"\ntotal wall time {out['wall_seconds']:.1f} s")
    if args.json:
        with open(args.json, "w") as fh:
            json.dump(out, fh, indent=2, default=float)
        print(f"record written to {args.json}")


if __name__ == "__main__":
    main()
