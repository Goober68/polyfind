"""Scan the overall electrostatic strength against three discrepancies at once.

Reproduces ``docs/ELECTROMECHANICS.md`` section 5.7.  The hypothesis under test was that a
single over-strong electrostatics explains all three of: the VDF/VDCN copolymer's
polar/antipolar energy gap coming out 3-6x larger than an MLIP reference, and beta-PVDF's
``d_33`` and ``d_31`` coming out 5x and 9x smaller than measurement.  It does not; the
numbers are in the document and the verdict is that these are two deficiencies, not one.

Every energy depends on the charges and the permittivity only through ``q^2 / eps_r``
(``docs/DFT_FIT.md`` section 5), so one scale ``s`` on that ratio is the whole electrostatic
strength *for energies*.  A polarization is linear in ``q`` and blind to ``eps_r``, so the
same ``s`` has two inequivalent realisations and only a piezoelectric coefficient can tell
them apart:

* ``--ray eps``: ``eps_r = 1/s`` at fixed charges -- energies scale, ``P`` does not;
* ``--ray q``: ``q -> sqrt(s) q`` at ``eps_r = 1`` -- energies scale identically and
  ``P`` scales as ``sqrt(s)``.

Usage (``PYTHONPATH=src``), each mode independent::

    python examples/electrostatic_scale.py beta --ray q            # d_33, d_31, C, cell
    python examples/electrostatic_scale.py beta --ray q --cap 25   # ... with the shape cap widened
    python examples/electrostatic_scale.py copolymer               # the polar/antipolar gap
    python examples/electrostatic_scale.py polymorph               # the three acceptance tests
    python examples/electrostatic_scale.py screening               # an explicit screening length
    python examples/electrostatic_scale.py coulomb-content         # the bound on the held-out error

``beta`` is about 6 s a point, ``polymorph`` 8 s, and ``copolymer`` about 200 s a point: it
relaxes both branches of a 148-atom two-chain cell under a full Ewald sum at every scale,
because the hypothesis is about the landscape's stiffness and a clamped cell would not test
it.  Nothing here is a default and no preset is registered: every parameter set is built as
a local :class:`~polyfind.fitting.FFParameters` and applied inside a context manager.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from dataclasses import replace

import numpy as np
from scipy.optimize import minimize

import polyfind.pack as pack_mod
from polyfind import mechanics as M
from polyfind.chain import build_chain
from polyfind.ewald import EwaldSpec
from polyfind.fitting import FITTED_VALENCE, acceptance_tests
from polyfind.forcefield import SimpleFF
from polyfind.lattice_table import clear_pair_table_cache
from polyfind.pack import periodic_chain, polish
from polyfind.polymers import AN, PVDC, PVDF, THREE_STATE, VDCN, VDF_VDCN_11_1

T = 0
BETA = (PVDF, [T, T], "PVDF beta (TTTT)")
SCALES = (0.0625, 0.125, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 4.0)
# The two delivered copolymer cells (``deliverables/summary.json``), which every scale
# starts from so that the scan compares relaxations of one basin rather than two searches.
CP = VDF_VDCN_11_1
C_REPEAT = 30.755727067694547
POLAR = np.array([10.45975535996956, 4.803557094896295, 90.0,
                  89.99999958340592, 90.0000002182697, 2.5304053003427294, 0.0])
ANTI = np.array([11.326746983379785, 4.745132768836892, 90.0,
                 90.00024128693593, 270.00024128693593, 3.8118291735798397, 0.0])
LO = np.array([4.0, 4.0, 90.0, 0.0, 0.0, 0.0])
HI = np.array([14.0, 14.0, 90.0, 360.0, 360.0, C_REPEAT])
FREE = [0, 1, 3, 4, 5]


def knobs(ray: str, s: float) -> tuple[float, float]:
    """``(eps_r, charge_scale)`` realising energy scale ``s`` along ``ray``."""
    if ray == "eps":
        return 1.0 / s, 1.0
    if ray == "q":
        return 1.0, float(np.sqrt(s))
    raise ValueError(f"unknown ray {ray!r} (expected 'eps' or 'q')")


def scaled(ray: str, s: float):
    """``(FFParameters, valence SimpleFF, flux SimpleFF)`` at electrostatic scale ``s``.

    The flux potential carries its own ``charge_scale`` because a packer given
    ``charge_flux`` recomputes the charges from *it* rather than from the chain, so scaling
    only the :class:`FFParameters` would leave the fluxing charges unscaled.
    """
    eps_r, q = knobs(ray, s)
    params = replace(FITTED_VALENCE, eps_r=eps_r, charge_scale=q)
    val = replace(SimpleFF.from_preset("pvdf-dft-valence"), eps_r=eps_r, charge_scale=q, _cache={})
    flx = replace(SimpleFF.from_preset("pvdf-dft-valence-flux"), eps_r=eps_r, charge_scale=q, _cache={})
    return params, val, flx


# ------------------------------------------------------------------ beta: d_33 and d_31
def beta_point(ray: str, s: float, cap: float | None = None, alpha: float = 0.2) -> dict:
    t0 = time.time()
    params, val, flx = scaled(ray, s)
    clear_pair_table_cache()
    with params.applied():
        polymer, seq, label = BETA
        ref, rr = M.refined_reference(polymer, seq, label=label, valence=val,
                                      refine_kw={"valence": val}, charge_flux=flx,
                                      pack_kw={"table_cache_dir": None, "alpha": alpha},
                                      alpha=alpha)
        shape = M.shape_of(polymer, ref, angles=rr.angles)
        if cap is not None:
            shape = replace(shape, max_angle_change=cap)
        # how far the one shape variable moves, and whether it is on its bound: past the cap
        # the reference is not a stationary point and neither d nor C_33 is a measurement
        rel = M.relax_deformable(ref, shape, ref.params, shape.x0, free=M._CELL_FREE,
                                 c_target=None, maxiter=300)
        resp = M.electromechanical_response(ref, polymer=polymer, shape=shape)
        P, el, pz = resp.polarization, resp.elastic, resp.piezo
        # film convention (docs/ELECTROMECHANICS.md section 5.3): axis 3 = poling = packer x,
        # axis 1 = draw = chain axis z; reachable columns are Voigt (1, 2, 3, 6)
        sgn = 1.0 if P[0] >= 0 else -1.0
        return {"ray": ray, "s": s, "cap": shape.max_angle_change, "alpha": alpha,
                "eps_r": params.eps_r, "charge_scale": params.charge_scale,
                "shape_x": float(rel.x[0]), "shape_grad": float(rel.shape_gradient),
                "on_cap": bool(abs(abs(rel.x[0]) - shape.max_angle_change) < 1e-6),
                "a": float(resp.reference.params[0]), "b": float(resp.reference.params[1]),
                "c": float(resp.reference.c), "absP": float(np.linalg.norm(P)),
                "C11": float(el.C[0, 0]), "C22": float(el.C[1, 1]), "C33": float(el.C[2, 2]),
                "d33_film": float(sgn * pz.d_improper[0, 0]),
                "d32_film": float(sgn * pz.d_improper[0, 1]),
                "d31_film": float(sgn * pz.d_improper[0, 2]),
                "route_rel_diff": float(pz.relative_difference),
                "res_zz": float(el.residual_stress[2]), "seconds": time.time() - t0}


def beta_row(r: dict) -> str:
    flag = "  NOT A MEASUREMENT" if (r["on_cap"] or r["route_rel_diff"] > 0.05) else ""
    return (f"{r['s']:8.4f} {r['eps_r']:8.3f} {r['charge_scale']:6.3f} {r['a']:7.3f} "
            f"{r['c']:7.4f} {r['absP']:7.4f} {r['C11']:7.2f} {r['C33']:8.2f} "
            f"{r['d33_film']:+9.3f} {r['d31_film']:+8.3f} {r['shape_x']:7.3f} "
            f"{r['route_rel_diff'] * 100:6.2f}%{flag}")


BETA_HEAD = (f"{'s':>8} {'eps_r':>8} {'q':>6} {'a':>7} {'c':>7} {'|P|':>7} {'C11':>7} "
             f"{'C33':>8} {'d33':>9} {'d31':>8} {'shape':>7} {'routes':>7}")


# --------------------------------------------------------- copolymer: polar vs antipolar
def relax_polar(pk, ch, p0, n_keep=2):
    """Twelve-point VDCN registry sweep, then an L-BFGS polish of the best few."""
    n = CP.monomers_per_repeat
    grid = np.repeat(np.asarray(p0, float)[None], n, axis=0)
    grid[:, 5] = (p0[5] + np.arange(n) * ch.c / n) % ch.c
    e = pk.energy(grid)
    best, best_e = None, np.inf
    for i in np.argsort(e)[:n_keep]:
        p = polish(pk, grid[i], LO, HI, FREE, maxfev=400, method="lbfgs", gradient="analytic")
        ee = float(pk.energy(p[None])[0])
        if ee < best_e:
            best, best_e = p, ee
    return best, best_e


def relax_anti(pk, ch, p0, n_keep=2):
    """The same sweep inside the antipolar subspace: ``phi2`` tied to ``phi1``, so ``P``
    stays exactly zero and the branch cannot drift into the polar one."""
    n = CP.monomers_per_repeat
    dphi, flip = float(p0[4] - p0[3]), float(p0[6])
    grid = np.repeat(np.asarray(p0, float)[None], n, axis=0)
    grid[:, 5] = (p0[5] + np.arange(n) * ch.c / n) % ch.c
    e = pk.energy(grid)

    def obj(v):
        a, b, phi, dz = v
        if a < 2.0 or b < 2.0:
            return 1e6
        return float(pk.energy(np.array([a, b, 90.0, phi, phi + dphi, dz, flip])[None])[0])

    best, best_e = None, np.inf
    for i in np.argsort(e)[:n_keep]:
        g = grid[i]
        r = minimize(obj, [g[0], g[1], g[3], g[5]], method="Nelder-Mead",
                     options={"maxfev": 400, "xatol": 1e-4, "fatol": 1e-6})
        p = np.array([r.x[0], r.x[1], 90.0, r.x[2], r.x[2] + dphi, r.x[3], flip])
        ee = float(pk.energy(p[None])[0])
        if ee < best_e:
            best, best_e = p, ee
    return best, best_e


def copolymer_point(ray: str, s: float, coulomb: str = "ewald", alpha: float = 0.2,
                    relaxed: bool = True) -> dict:
    """The polar minus antipolar energy per monomer, both branches relaxed at this scale.

    Ewald is the only sum in which this difference means anything (``polyfind.pack``), with
    the exception ``screening`` exercises: under a short enough screening length the dipole
    lattice sum converges absolutely and a damped truncated sum is legitimate.
    """
    t0 = time.time()
    eps_r, q = knobs(ray, s)
    params = replace(FITTED_VALENCE, eps_r=eps_r, charge_scale=q)
    with params.applied():
        ch = periodic_chain(CP, [0] * CP.bonds_per_repeat, THREE_STATE)
        kw = dict(n_chains=2, alpha=alpha)
        if coulomb == "ewald":
            kw.update(coulomb="ewald", ewald=EwaldSpec())
        pk = pack_mod.CrystalPacker(ch, **kw)
        nm = 2 * ch.n_monomers
        if relaxed:
            pp, ep = relax_polar(pk, ch, POLAR)
            pa, ea = relax_anti(pk, ch, ANTI)
        else:
            pp, pa = POLAR, ANTI
            ep, ea = float(pk.energy(pp[None])[0]), float(pk.energy(pa[None])[0])
        return {"ray": ray, "s": s, "coulomb": coulomb, "alpha": alpha, "relaxed": relaxed,
                "E_polar": ep / nm, "E_anti": ea / nm, "gap": (ep - ea) / nm,
                "polar_a": float(pp[0]), "polar_b": float(pp[1]),
                "anti_a": float(pa[0]), "anti_b": float(pa[1]),
                "rho_polar": float(pk.density(pp[0], pp[1], pp[2])),
                "rho_anti": float(pk.density(pa[0], pa[1], pa[2])),
                "absP_polar": float(np.linalg.norm(pk.polarization(pp))),
                "absP_anti": float(np.linalg.norm(pk.polarization(pa))),
                "seconds": time.time() - t0}


# --------------------------------------------------------------- the cost: what is right
def polymorph_point(ray: str, s: float) -> dict:
    """The three acceptance tests at this scale; identical along both rays, by construction."""
    params, _, _ = scaled(ray, s)
    clear_pair_table_cache()
    t0 = time.time()
    r = acceptance_tests(params)
    beta = r["predictions"]["beta"]
    return {"ray": ray, "s": s, "alpha_beta_kj": float(r["alpha_beta_kj"]),
            "alpha_beta_pass": bool(r["alpha_beta_pass"]), "ris_top": r["ris_top"],
            "alpha_polar_gap": float(r["alpha_polar_gap"]),
            "beta_absP": float(r["beta_polarization"]),
            "beta_a": float(beta.a), "beta_b": float(beta.b), "beta_c": float(beta.c),
            "seconds": time.time() - t0}


def coulomb_content(polymer, ff, n=240, n_dih=7, seed=0, window=20.0) -> dict:
    """RMS of the mean-centred Coulomb part of relative conformer energies, kcal/mol.

    The held-out energy error cannot be re-measured without ``$POLYFIND_TRAINSET``, which is
    not vendored.  This bounds what the scan does to it: scaling the electrostatics by ``s``
    shifts every frame energy by ``(s - 1) E_coul``, and the fit's residual is on energies
    with each system's own mean removed, so this number times ``|s - 1|`` is the RMS shift
    imposed on exactly the quantity the 1.36 kcal/mol measures.  ``window`` keeps only
    conformers within 20 kcal/mol of the minimum: the reference set's kept systems span 3.6
    to 17 kcal/mol, so a rigid-geometry clash is not a frame it contains.
    """
    rng = np.random.default_rng(seed)
    dih = np.array((180.0, 70.0, -70.0))[rng.choice(3, size=(n, n_dih))]
    tot, coul = [], []
    for d in dih:
        c = ff.components(build_chain(polymer, d))
        e = sum(c.values())
        if np.isfinite(e) and e < 1e4:
            tot.append(e)
            coul.append(c["coulomb"])
    tot, coul = np.array(tot), np.array(coul)
    keep = tot - tot.min() <= window
    tot, coul = tot[keep] - tot[keep].mean(), coul[keep] - coul[keep].mean()
    return {"polymer": polymer.name, "n_kept": int(len(tot)),
            "rms_total_centred": float(np.sqrt((tot ** 2).mean())),
            "rms_coulomb_centred": float(np.sqrt((coul ** 2).mean()))}


# ----------------------------------------------------------------------------- the modes
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["beta", "copolymer", "polymorph", "screening",
                                     "coulomb-content"])
    ap.add_argument("--ray", default="q", choices=["eps", "q"])
    ap.add_argument("--scales", type=float, nargs="*", default=list(SCALES))
    ap.add_argument("--cap", type=float, default=None,
                    help="Shape.max_angle_change in degrees; 8 is the shipped value and "
                         "every scale past about 2 sits on it")
    ap.add_argument("--screening", type=float, nargs="*", default=[0.1, 0.2, 0.4, 0.6, 1.0],
                    help="damped-shifted-force alpha, i.e. 1/(screening length in A)")
    ap.add_argument("--frozen", action="store_true",
                    help="do NOT relax the copolymer cells; they are then the delivered ones")
    ap.add_argument("--json", default=None, help="also write every row to this file")
    args = ap.parse_args()

    rows = []
    if args.mode == "beta":
        print(f"# beta-PVDF, pvdf-dft-valence(-flux), ray={args.ray}, cap="
              f"{'8 (shipped)' if args.cap is None else args.cap}")
        print(BETA_HEAD)
        for s in args.scales:
            r = beta_point(args.ray, s, cap=args.cap)
            rows.append(r)
            print(beta_row(r), flush=True)
    elif args.mode == "copolymer":
        print(f"# 11VDF:1VDCN polar minus antipolar, Ewald, pvdf-dft-valence, "
              f"{'cells as delivered' if args.frozen else 'both branches relaxed'}")
        print(f"{'s':>8} {'E_polar':>9} {'E_anti':>9} {'gap':>8} {'rho_p':>6} {'rho_a':>6} "
              f"{'a_p':>7} {'a_a':>7} {'|P|_p':>7} {'sec':>5}")
        for s in args.scales:
            r = copolymer_point(args.ray, s, relaxed=not args.frozen)
            rows.append(r)
            print(f"{r['s']:8.4f} {r['E_polar']:9.4f} {r['E_anti']:9.4f} {r['gap']:8.4f} "
                  f"{r['rho_polar']:6.3f} {r['rho_anti']:6.3f} {r['polar_a']:7.3f} "
                  f"{r['anti_a']:7.3f} {r['absP_polar']:7.4f} {r['seconds']:5.0f}", flush=True)
        if len(rows) > 2:
            g1, g0 = np.polyfit([r["s"] for r in rows], [r["gap"] for r in rows], 1)
            print(f"  linear in s: gap = {g0:+.4f} {g1:+.4f} s -- the intercept is the gap "
                  f"with the electrostatics switched off", flush=True)
    elif args.mode == "polymorph":
        print("# the three acceptance tests against the same scale (identical on both rays)")
        print(f"{'s':>8} {'a-b kJ':>8} {'pass':>5} {'RIS top':>9} {'anti':>7} {'|P|beta':>8} "
              f"{'beta a':>7} {'beta b':>7} {'beta c':>7}")
        for s in args.scales:
            r = polymorph_point(args.ray, s)
            rows.append(r)
            print(f"{r['s']:8.4f} {r['alpha_beta_kj']:8.2f} {str(r['alpha_beta_pass']):>5} "
                  f"{r['ris_top']:>9} {r['alpha_polar_gap']:+7.3f} {r['beta_absP']:8.4f} "
                  f"{r['beta_a']:7.3f} {r['beta_b']:7.3f} {r['beta_c']:7.4f}", flush=True)
    elif args.mode == "screening":
        print("# an explicit screening length: qq erfc(alpha r)/r leaves r -> 0 alone and "
              "kills the tail")
        print(f"{'1/alpha':>8} {'gap (dsf)':>10} {'d33':>9} {'d31':>8} {'C11':>7} {'|P|':>7}")
        for a in args.screening:
            b = beta_point("q", 1.0, alpha=a)
            c = copolymer_point("q", 1.0, coulomb="dsf", alpha=a)
            rows.append({"alpha": a, **{f"beta_{k}": v for k, v in b.items()},
                         **{f"copo_{k}": v for k, v in c.items()}})
            print(f"{1 / a:8.2f} {c['gap']:10.4f} {b['d33_film']:+9.3f} {b['d31_film']:+8.3f} "
                  f"{b['C11']:7.2f} {b['absP']:7.4f}", flush=True)
    else:
        ff = SimpleFF.from_preset("pvdf-dft-valence")
        print("# RMS of the mean-centred Coulomb part of relative conformer energies "
              "(kcal/mol)")
        rows = [coulomb_content(p, ff) for p in (PVDF, PVDC, AN, VDCN)]
        for r in rows:
            print(f"  {r['polymer']:6s} n={r['n_kept']:4d}  total {r['rms_total_centred']:7.3f}"
                  f"  coulomb {r['rms_coulomb_centred']:7.3f}")
        print()
        print(f"{'s':>8} " + " ".join(f"{r['polymer']:>8}" for r in rows)
              + "   RMS shift imposed on the held-out energy residual")
        for s in (0.0625, 0.25, 0.5, 2.0, 4.0):
            print(f"{s:8.4f} " + " ".join(
                f"{abs(s - 1) * r['rms_coulomb_centred']:8.2f}" for r in rows))

    if args.json:
        with open(args.json, "w") as f:
            json.dump(rows, f, indent=1)


if __name__ == "__main__":
    sys.exit(main())
