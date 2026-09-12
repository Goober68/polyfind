"""Nitrile chains: is the kinked geometry a missing RIS state or a frustrated chain?

Four independent sightings say nitrile-bearing backbones (VDCN, AN) leave the three-state
rotational-isomeric description where PVDF does not: Sarco's kinked reference, their two
converged basins with distortions 2-10 bonds from the nitrile, the polyacrylonitrile
literature's loss of chain periodicity, and their finite-chain ladder on which PVDF passes all
six size gates and the nitriles fail all six (docs/REFERENCE_DATA_REQUEST.md, 2026-09-11).
Two readings with opposite consequences: a **method** gap (a fourth state, or continuous
torsions, would represent these structures) or a **chemistry** fact (no well-defined
single-chain minimum, so no classical tier ranks them).  This script runs the three
measurements that discriminate, all under ``pvdf-dft-valence``:

    python examples/nitrile_landscape.py --part profiles     # 1. relaxed dihedral profiles
    python examples/nitrile_landscape.py --part perturb      # 2. perturb-and-relax returns
    python examples/nitrile_landscape.py --part correlation  # 3. torsional correlation length
    python examples/nitrile_landscape.py --part risfit       # 4. the RIS fit with angles relaxed
    python examples/nitrile_landscape.py --part all --json out.json
    python examples/nitrile_landscape.py --part perturb --quick   # 3 amplitudes x 3 seeds

1. **Profiles.**  One backbone dihedral is driven round the circle on an eleven-bond oligomer
   at three levels: rigid (every other bond trans, the polymer's frozen backbone angles --
   exactly ``fit_ris``'s own 1D scan), angle-relaxed (every backbone angle of the oligomer
   free, others trans) and fully relaxed (angles plus every other torsion free, each point a
   local minimisation from the angle-relaxed all-trans geometry, so the profile is that of one
   bond with the chain otherwise near its reference rather than a hysteretic walk).  Local
   minima are listed; a minimum more than 30 deg from 180 and from +/-60 is a basin the
   three-state model has no state for.  Then every rigid-profile minimum is relaxed with
   nothing held, to see whether it survives as a minimum or was frozen-angle strain.  Bonds
   scanned: PVDF's (control), both bond types of VDCN and AN (every bond of either homopolymer
   is adjacent to a nitrile carbon), and in a VDF3-VDCN-VDF3 repeat the bonds 0, 1, 2 and 3
   bonds from the nitrile-bearing carbon -- the distances at which Sarco finds the kinks.
2. **Perturbation returns.**  A single periodic chain (one chain in a 30 A cell, 8 A cutoff, so
   it sees only its own axial images) starts from the RIS ground state and is relaxed by
   :func:`polyfind.refine.refine_crystal` after its torsions are perturbed by Gaussian noise
   of increasing width.  The randomised torsions follow no line-group pattern, so the
   refinement takes its documented ``"free"`` (penalty) path for every run; the bend terms are
   the fitted ones (``valence=``).  Endpoints are compared under the chain's own symmetries
   (shift by a repeat, reversal, and for an achiral repeat the mirror): how many distinct
   conformations, how wide their energy spread, how often the reference minimum is reached,
   and how many torsions end more than 30 deg from every RIS state (Sarco's kink criterion).
   VDCN's RIS ground state is all-trans (the screen's TT).  PVDF's fitted RIS ground state is
   a helix that is not commensurate within eight periods and cannot be built as a repeat, so
   PVDF is run twice: all-trans (beta; the like-for-like control for VDCN's TT) and TG+TG-
   (alpha, the known lattice ground state).  AN all-trans is included as a second nitrile,
   with the caveat that TT is not AN's RIS ground state.
3. **Correlation length.**  The screen's own fitted RIS models (``fit_ris``, step 10, six
   monomers, third order) are sampled exactly with :mod:`polyfind.amorphous` at 200-500 K, and
   the torsional correlation length is read two ways: exactly, from the second eigenvalue of
   the repeat-unit transfer matrix; and as the run length over which a sampled chain follows
   its own RIS ground-state pattern.  A crystal needs that run to be longer than a repeat.
4. **The RIS fit with the angles relaxed.**  Part 1 shows what the frozen backbone angles do
   to a nitrile scan, so the screen's own ``fit_ris`` is repeated with a calculator that
   relaxes the repeat's backbone angles for every conformer it scores (the freedom
   ``refine_crystal`` already has, applied inside the scan) and the fitted state angles,
   first-order energies, ground state and all-trans rank are compared with the rigid fit's.
   If the fix on our side gives a self-consistent model, that is the method reading shown
   rather than argued.

Caveat that applies to every number here: one classical potential, fitted to PVDF data, whose
nitrile torsion gets PVDF's CF2-CH2 Fourier triple (``SimpleFF.torsion`` is one triple) and
whose C-N increment and C=N stretch were fitted to a reference set in which the nitriles were
a minority.  The oligomer energy path refuses a fluxing potential by design
(``SimpleFF._refuse_flux``), so parts 1 and 2 use ``pvdf-dft-valence`` -- the preset the
screen's RIS fit, packing and refinement used; the flux enters the screen only in the
mechanics packer.  docs/NITRILE_LANDSCAPE.md records the results and the verdict.
"""
from __future__ import annotations

import argparse
import json
import time
import warnings

import numpy as np
from scipy.optimize import minimize

import polyfind.fitting  # noqa: F401  (registers the fitted presets)
from polyfind.amorphous import ensemble_stats, sample_ensemble
from polyfind.chain import build_chain_batch
from polyfind.enumerate import enumerate_periodic
from polyfind.fitting import FITTED_VALENCE
from polyfind.forcefield import SimpleFF, fit_ris
from polyfind.pack import periodic_chain
from polyfind.polymers import PVDF, THREE_STATE, VDCN, copolymer, get_polymer, monomer_of
from polyfind import refine as refine_mod
from polyfind.refine import refine_crystal

warnings.simplefilter("ignore")
VAL = SimpleFF.from_preset("pvdf-dft-valence")
IDEAL = np.array([180.0, 60.0, -60.0])  # the standard RIS state angles Sarco's criterion uses
KINK = 30.0  # deg from every state
CLOSURE_TOL = 0.5  # deg of residual repeat rotation above which a "free" refinement did not close


def wrap(x):
    return (np.asarray(x, dtype=float) + 180.0) % 360.0 - 180.0


def dist_to_states(phi, states=IDEAL):
    phi = np.asarray(phi, dtype=float)
    return np.abs(wrap(phi[..., None] - np.asarray(states)[None, :])).min(axis=-1)


# --------------------------------------------------------------------------- part 1
def oligomer_energies(polymer, dih, ang):
    """Energies of a batch of oligomers with per-row torsions ``(M, N)`` and per-atom
    backbone angles ``(M, N+3)``."""
    tpl, coords = build_chain_batch(polymer, dih, bond_angles=ang)
    return VAL.energy_batch((tpl, coords))


def relax_oligomer(polymer, N, fixed, tors0, ang0, free_torsions=True, lo=95.0, hi=135.0, maxiter=400):
    """Minimise the oligomer energy over the free torsions and every backbone angle.

    ``fixed`` is ``{dihedral index: value}``; the other torsions start from ``tors0`` and are
    free when ``free_torsions``; every backbone atom of the oligomer carries its own angle,
    bounded to ``[lo, hi]`` as :func:`polyfind.forcefield.relax_backbone_angles` bounds them.
    L-BFGS-B on a batched central-difference gradient (one ``build_chain_batch`` per step).
    """
    tors0 = np.asarray(tors0, dtype=float)
    ang0 = np.asarray(ang0, dtype=float)
    free = [j for j in range(N) if j not in fixed] if free_torsions else []
    nt = len(free)
    x0 = np.concatenate([tors0[free], ang0])
    h = np.concatenate([np.full(nt, 0.5), np.full(N + 3, 0.2)])

    def unpack(X):
        X = np.atleast_2d(X)
        dih = np.repeat(tors0[None], X.shape[0], axis=0)
        for j, v in fixed.items():
            dih[:, j] = v
        if nt:
            dih[:, free] = X[:, :nt]
        return dih, X[:, nt:]

    def fg(x):
        rows = [x]
        for j in range(len(x)):
            for s in (1.0, -1.0):
                y = x.copy()
                y[j] += s * h[j]
                rows.append(y)
        dih, ang = unpack(np.array(rows))
        E = oligomer_energies(polymer, dih, ang)
        return float(E[0]), (E[1::2] - E[2::2]) / (2 * h)

    bounds = [(None, None)] * nt + [(lo, hi)] * (N + 3)
    res = minimize(fg, x0, jac=True, method="L-BFGS-B", bounds=bounds,
                   options={"maxiter": maxiter, "ftol": 1e-10, "gtol": 1e-4})
    dih, ang = unpack(res.x[None])
    return float(res.fun), dih[0], ang[0], int(res.nit)


def local_minima(grid, E):
    """Indices of the local minima of a periodic profile and their escape barriers."""
    n = len(E)
    out = []
    for i in range(n):
        if E[i] < E[(i - 1) % n] and E[i] <= E[(i + 1) % n]:
            barriers = []
            for step in (1, -1):
                k, top = i, E[i]
                while True:
                    k = (k + step) % n
                    top = max(top, E[k])
                    if E[k] > E[(k + step) % n] and k != i:  # passed a maximum
                        break
                    if k == i:
                        break
                barriers.append(top - E[i])
            out.append((i, float(grid[i]), float(E[i]), float(min(barriers))))
    return out


def profile(polymer, N, j0, grid):
    B = polymer.bonds_per_repeat
    ang_poly = np.array([polymer.backbone[k % B].backbone_angle for k in range(N + 3)])
    trans = np.full(N, 180.0)
    e_rigid_ref = float(oligomer_energies(polymer, trans[None], ang_poly[None])[0])
    e_ref, _, ang_ref, _ = relax_oligomer(polymer, N, {}, trans, ang_poly, free_torsions=False)
    rigid, angles, full, drift = [], [], [], []
    rows = np.repeat(trans[None], len(grid), axis=0)
    rows[:, j0] = grid
    rigid = oligomer_energies(polymer, rows, np.repeat(ang_poly[None], len(grid), axis=0)) - e_rigid_ref
    for phi in grid:
        e_a, _, _, _ = relax_oligomer(polymer, N, {j0: float(phi)}, trans, ang_ref, free_torsions=False)
        angles.append(e_a - e_ref)
        e_f, d_f, _, _ = relax_oligomer(polymer, N, {j0: float(phi)}, trans, ang_ref, free_torsions=True)
        full.append(e_f - e_ref)
        dev = np.abs(wrap(d_f - 180.0))
        dev[j0] = 0.0
        drift.append(float(dev.max()))
    out = {"grid": grid.tolist(), "rigid": np.asarray(rigid).tolist(), "angles": angles, "full": full,
           "neighbour_drift": drift, "e_rigid_all_trans": e_rigid_ref, "e_relaxed_all_trans": e_ref,
           "relaxed_angles_all_trans": ang_ref.tolist()}
    for level in ("rigid", "angles", "full"):
        out[level + "_minima"] = local_minima(grid, np.asarray(out[level]))
    # the basin test: relax every rigid minimum with nothing held
    basin = []
    for _, phi, e_r, _ in out["rigid_minima"]:
        t0 = trans.copy()
        t0[j0] = phi
        e_f, d_f, _, nit = relax_oligomer(polymer, N, {}, t0, ang_ref, free_torsions=True)
        dev = np.abs(wrap(d_f - 180.0))
        dev[j0] = 0.0
        basin.append({"rigid_phi": phi, "rigid_dE": e_r, "relaxed_phi": float(wrap(d_f[j0])),
                      "relaxed_dE": e_f - e_ref, "others_drift": float(dev.max()), "iterations": nit})
    out["basin_test"] = basin
    return out


def profile_systems():
    vdf, vdcn = monomer_of(PVDF, "vdf"), monomer_of(VDCN, "vdcn")
    co = copolymer("vdf3-vdcn-vdf3", [vdf] * 3 + [vdcn] + [vdf] * 3, bond_length=PVDF.bond_length,
                   formula="-(CH2-CF2)3-(CH2-C(CN)2)-(CH2-CF2)3-")
    N = 11  # 14 backbone atoms: exactly one copolymer repeat, nitrile carbon at atom 7
    systems = [("pvdf", PVDF, 5, "control: PVDF, CF2-CH2 bond")]
    for name in ("vdcn", "an"):
        p = get_polymer(name)
        for j0 in (4, 5):
            systems.append((name, p, j0, f"{name}: bond type {j0 % 2} (adjacent to the nitrile carbon)"))
    for j0 in range(2, 10):
        d_bond = min(abs(j0 + 1 - 7), abs(j0 + 2 - 7))  # central bond (j0+1, j0+2) to atom 7
        systems.append(("vdf3-vdcn-vdf3", co, j0, f"copolymer: bond {d_bond} from the nitrile carbon "
                        f"(four-atom window {max(0, d_bond - 1)} away)"))
    return N, systems


def run_profiles(step=10.0, verbose=True):
    grid = np.arange(-180.0, 180.0, step)
    N, systems = profile_systems()
    out = []
    for name, p, j0, label in systems:
        t0 = time.time()
        prof = profile(p, N, j0, grid)
        prof.update({"polymer": name, "dihedral": j0, "label": label, "seconds": round(time.time() - t0, 1)})
        out.append(prof)
        if verbose:
            print(f"\n{label}   [dihedral {j0} of an {N}-bond oligomer, {prof['seconds']} s]")
            print(f"  all-trans: rigid {prof['e_rigid_all_trans']:8.2f}, angle-relaxed {prof['e_relaxed_all_trans']:8.2f} kcal/mol"
                  f"  (relaxed angles {', '.join(f'{a:.0f}' for a in prof['relaxed_angles_all_trans'][1:-1])})")
            for level in ("rigid", "angles", "full"):
                mins = prof[level + "_minima"]
                txt = "; ".join(f"{phi:+5.0f} deg: {e:+6.2f} (barrier {b:4.2f})" for _, phi, e, b in mins)
                kinks = [phi for _, phi, _, _ in mins if dist_to_states(phi) > KINK]
                print(f"  {level:7s} minima: {txt}" + (f"   <- {len(kinks)} in the kink region" if kinks else ""))
            print(f"  full profile, neighbour drift from trans: max {max(prof['neighbour_drift']):.0f} deg")
            for b in prof["basin_test"]:
                print(f"  basin test: rigid minimum at {b['rigid_phi']:+5.0f} ({b['rigid_dE']:+6.2f}) relaxes to "
                      f"{b['relaxed_phi']:+6.1f} deg at {b['relaxed_dE']:+6.2f}; others moved {b['others_drift']:.0f} deg")
    return out


# --------------------------------------------------------------------------- part 2
def torsion_images(tors, B, chiral):
    """Images of a repeat's torsion vector under shift by a repeat unit, reversal, and (achiral) mirror.

    Reversal maps dihedral ``j`` (about the bond between atoms ``j+1`` and ``j+2``) to
    ``(P - 3 - j) mod P`` and leaves its value unchanged; the mirror negates every torsion.  A
    chiral repeat has neither alone, only their composition (:func:`polyfind.helix.sequence_images`).
    """
    t = np.asarray(tors, dtype=float)
    P = len(t)
    rev = t[[(P - 3 - j) % P for j in range(P)]]
    base = [t, -rev] if chiral else [t, rev, -t, -rev]
    return [np.roll(b, -r) for b in base for r in range(0, P, B)]


def torsion_distance(a, b, B, chiral):
    return min(float(np.abs(wrap(img - np.asarray(b))).max()) for img in torsion_images(a, B, chiral))


def ris_string(tors, states=IDEAL, names=("T", "G+", "G-")):
    t = np.asarray(tors, dtype=float)
    d = np.abs(wrap(t[:, None] - np.asarray(states)[None, :]))
    near = d.argmin(axis=1)
    return "".join("K" if d[i, near[i]] > KINK else names[near[i]] for i in range(len(t)))


def perturbation_returns(polymer, seq_name, n_bonds=8, sigmas=(2, 5, 10, 20, 30, 45, 60), seeds=8,
                         maxiter=400, from_relaxed=False, verbose=True):
    """Perturb-and-relax sweep of one periodic chain.

    ``from_relaxed=False`` perturbs the RIS chain at its ideal angles, as the question was put.
    ``from_relaxed=True`` first relaxes that chain (four seeds at 2 deg, lowest taken) and
    perturbs the relaxed torsions instead: the ideal-angle start sits 13 kcal/mol per monomer
    above VDCN's relaxed zigzag and 0.4 above PVDF's, so a chain that scatters from the ideal
    start might only be falling a long way, and this variant is the control for that.
    """
    B = polymer.bonds_per_repeat
    seq = THREE_STATE.names
    idx = [seq.index(s) for s in _parse(seq_name)]
    seq_idx = (idx * (n_bonds // len(idx)))[:n_bonds]
    chiral = bool(polymer.is_chiral)
    runs = []
    with FITTED_VALENCE.applied():
        ch = periodic_chain(polymer, seq_idx, THREE_STATE)
        pk = refine_mod.CrystalPacker(ch, n_chains=1, cutoff=8.0, valence=VAL)
        start = pk.result(np.array([30.0, 30.0, 90.0, 0.0, 0.0, 0.0, 0.0]))
        e_start = start.energy_per_monomer
        base = ch.dihedrals.copy()

        def one(sigma, s, base_tors):
            rng = np.random.default_rng(100000 * int(sigma) + 17 * s + 1)
            tors = base_tors + rng.normal(0.0, sigma, len(base_tors))
            t0 = time.time()
            rr = refine_crystal(polymer, start, torsions=tors, valence=VAL, parametrisation="free",
                                max_torsion_change=180.0, max_angle_change=15.0,
                                maxiter=maxiter, maxfev=4 * maxiter)
            return {"sigma": sigma, "seed": s, "e": float(rr.result.energy_per_monomer),
                    "tors": [float(x) for x in wrap(rr.torsions)], "angles": [float(x) for x in rr.angles],
                    "rotation_error": float(rr.rotation_error), "evaluations": int(rr.n_evaluations),
                    "mode": rr.parametrisation, "seconds": round(time.time() - t0, 1),
                    "start_tors": [float(x) for x in wrap(tors)]}

        if from_relaxed:
            pre = [one(2, s, base) for s in range(4)]
            pre = [r for r in pre if r["rotation_error"] <= CLOSURE_TOL] or pre
            best = min(pre, key=lambda r: r["e"])
            base = np.asarray(best["tors"], dtype=float)
        if verbose:
            print(f"\n{polymer.name} {seq_name} x {n_bonds // len(idx)}: RIS start (ideal angles, unrelaxed) "
                  f"E/mon = {e_start:+.3f} kcal/mol, c = {start.c:.3f} A"
                  + (f"; perturbing the RELAXED chain at {best['e']:+.3f} ({ris_string(base)}, "
                     + " ".join(f"{t:+.0f}" for t in base) + ")" if from_relaxed else ""))
        for sigma in sigmas:
            for s in range(seeds):
                try:
                    runs.append(one(sigma, s, base))
                except Exception as e:  # noqa: BLE001  (record the failure rather than lose the sweep)
                    runs.append({"sigma": sigma, "seed": s, "failed": str(e)})
    # A run whose repeat transform is still a rotation is a failed closure of the penalty
    # method, not a minimum of the chain: it is counted and excluded, never averaged in.
    ok = [r for r in runs if "e" in r and r["rotation_error"] <= CLOSURE_TOL]
    ref_pool = [r for r in ok if r["sigma"] == min(sigmas)] or ok
    ref = min(ref_pool, key=lambda r: r["e"])
    summary = []
    for sigma in sigmas:
        rows = [r for r in ok if r["sigma"] == sigma]
        n_unclosed = sum(1 for r in runs if r["sigma"] == sigma and r not in rows)
        if not rows:
            continue
        # greedy clustering under the chain's symmetries, 20 deg tolerance
        clusters = []
        for r in rows:
            for c in clusters:
                if torsion_distance(r["tors"], c["tors"], B, chiral) <= 20.0:
                    c["n"] += 1
                    c["e_min"] = min(c["e_min"], r["e"])
                    break
            else:
                clusters.append({"tors": r["tors"], "n": 1, "e_min": r["e"], "ris": ris_string(r["tors"])})
        es = np.array([r["e"] for r in rows])
        returned = [torsion_distance(r["tors"], ref["tors"], B, chiral) <= 20.0 for r in rows]
        kinks = [int((dist_to_states(r["tors"]) > KINK).sum()) for r in rows]
        strings = {min(ris_string(img) for img in torsion_images(r["tors"], B, chiral)) for r in rows}
        row = {"sigma": sigma, "n": len(rows), "n_unclosed": n_unclosed, "n_distinct": len(clusters), "n_ris_strings": len(strings),
               "return_fraction": float(np.mean(returned)), "e_min": float(es.min()), "e_max": float(es.max()),
               "e_spread": float(es.max() - es.min()), "e_below_ref": float(es.min() - ref["e"]),
               "kinks_mean": float(np.mean(kinks)), "kinks_max": int(max(kinks)),
               "rotation_error_max": float(max(r["rotation_error"] for r in rows)),
               "clusters": sorted(clusters, key=lambda c: c["e_min"])}
        summary.append(row)
        if verbose:
            print(f"  sigma={sigma:3.0f}: {row['n']:2d} closed runs ({n_unclosed} not closed) -> {row['n_distinct']:2d} distinct endpoints "
                  f"({row['n_ris_strings']} RIS strings), returned {100 * row['return_fraction']:3.0f}%, "
                  f"E/mon min {row['e_min']:+.3f} max {row['e_max']:+.3f} (spread {row['e_spread']:.3f}, "
                  f"lowest {row['e_below_ref']:+.3f} vs reference), kinks per repeat mean {row['kinks_mean']:.1f} "
                  f"max {row['kinks_max']}, max rot. residual {row['rotation_error_max']:.2e} deg")
            for c in row["clusters"][:4]:
                print(f"      x{c['n']:<2d} {c['ris']:<10s} E/mon {c['e_min']:+.3f}  tors "
                      + " ".join(f"{t:+5.0f}" for t in c["tors"]))
    return {"polymer": polymer.name, "sequence": seq_name, "n_bonds": n_bonds, "chiral": chiral,
            "from_relaxed": from_relaxed, "base_tors": [float(x) for x in base],
            "e_start": float(e_start), "reference": {"e": ref["e"], "tors": ref["tors"], "ris": ris_string(ref["tors"]),
                                                     "sigma": ref["sigma"], "seed": ref["seed"]},
            "summary": summary, "runs": runs}


def _parse(name):
    out, i = [], 0
    for_names = sorted(THREE_STATE.names, key=len, reverse=True)
    while i < len(name):
        for nm in for_names:
            if name.startswith(nm, i):
                out.append(nm)
                i += len(nm)
                break
        else:
            raise ValueError(name)
    return out


def run_perturb(quick=False, from_relaxed=False, verbose=True):
    sigmas = (5, 20, 45) if quick else (2, 5, 10, 20, 30, 45, 60)
    seeds = 3 if quick else 8
    systems = ([("vdcn", "TT"), ("pvdf", "TT")] if from_relaxed
               else [("vdcn", "TT"), ("pvdf", "TT"), ("pvdf", "TG+TG-"), ("an", "TT")])
    return [perturbation_returns(get_polymer(n), s, sigmas=sigmas, seeds=seeds, from_relaxed=from_relaxed,
                                 verbose=verbose) for n, s in systems]


# --------------------------------------------------------------------------- part 3
def transfer_correlation_length(model, T):
    """Correlation length (bonds) from the two leading eigenvalues of the repeat-unit transfer matrix."""
    m = model._aug() if model.third_order is not None else model
    lw1, lw2 = m._log_weights(T)
    U = np.eye(m.S)
    for b in range(m.B):
        U = U @ np.exp(lw1[b][:, None] + lw2[b])
    lam = np.sort(np.abs(np.linalg.eigvals(U)))[::-1]
    ratio = lam[1] / lam[0]
    xi = np.inf if ratio >= 1.0 - 1e-15 else -1.0 / np.log(ratio)
    return float(xi * model.B), float(ratio)


def pattern_runs(states, pattern, B):
    """Runs (bonds) over which a sampled chain follows ``pattern`` at one fixed phase.

    Pooled over the phases (shifts by a repeat unit); also the per-bond longest run, for the
    fraction of bonds sitting in a run of at least ``k`` bonds.
    """
    states = np.asarray(states)
    P, (n_chains, n) = len(pattern), states.shape
    runs, longest = [], np.zeros(states.shape, dtype=int)
    for sh in range(0, P, B):
        pat = np.array([pattern[(i + sh) % P] for i in range(n)])
        match = states == pat[None, :]
        for c in range(n_chains):
            x = np.concatenate([[False], match[c], [False]]).astype(int)
            d = np.diff(x)
            starts, ends = np.where(d == 1)[0], np.where(d == -1)[0]
            for a, b in zip(starts, ends):
                runs.append(b - a)
                longest[c, a:b] = np.maximum(longest[c, a:b], b - a)
    return np.array(runs), longest


def run_correlation(temperatures=(200.0, 300.0, 400.0, 500.0), n_chains=2000, n_bonds=400, verbose=True):
    out = []
    for name in ("pvdf", "vdcn", "an"):
        p = get_polymer(name)
        # the document's correlation lengths are those of the screen's *rigid* models (the
        # nitriles now default to the angle-relaxed scan, which section 4 motivated)
        fit = fit_ris(p, VAL, step=10.0, n_monomers=6, third_order=True, angles="rigid")
        m = fit.model
        gs = enumerate_periodic(p, m, max_period=8, k_per_period=40)[0]
        rec = {"polymer": name, "ground_state": gs.name, "period": gs.period,
               "state_angles": [float(a) for a in m.states.angles], "rows": []}
        if verbose:
            print(f"\n{name}: fitted RIS ground state {gs.name} (period {gs.period} bonds), "
                  f"state angles {tuple(round(a) for a in m.states.angles)}")
        for T in temperatures:
            xi, ratio = transfer_correlation_length(m, T)
            states, _ = sample_ensemble(p, m, n_bonds=n_bonds, n_chains=n_chains, T=T, build=False)
            st = ensemble_stats(p, m, states, None, T)
            runs, longest = pattern_runs(states, list(gs.seq), p.bonds_per_repeat)
            row = {"T": T, "xi_transfer_bonds": xi, "lambda_ratio": ratio,
                   "state_fractions": {k: round(v, 4) for k, v in st.state_fractions.items()},
                   "trans_run_mean": st.trans_run_mean, "p_trans_run_ge_8": st.trans_run_p_ge[8],
                   "pattern_run_mean": float(runs.mean()) if runs.size else 0.0,
                   "p_bond_in_pattern_run_ge_8": float((longest >= 8).mean()),
                   "p_bond_in_pattern_run_ge_16": float((longest >= 16).mean())}
            rec["rows"].append(row)
            if verbose:
                print(f"  T={T:4.0f} K: xi = {xi:7.2f} bonds (|l2/l1| = {ratio:.3f}); "
                      f"T fraction {row['state_fractions'].get('T', 0):.3f}; trans run {st.trans_run_mean:5.2f}; "
                      f"ground-pattern run {row['pattern_run_mean']:6.2f} bonds, "
                      f"P(bond in pattern run >= 8) {row['p_bond_in_pattern_run_ge_8']:.3f}, >= 16: {row['p_bond_in_pattern_run_ge_16']:.3f}")
        out.append(rec)
    return out


# --------------------------------------------------------------------------- part 4
class AngleRelaxedCalculator:
    """``VAL`` with the repeat's backbone angles relaxed for every conformer it scores.

    The freedom is the one :func:`polyfind.refine.refine_crystal` and
    :func:`polyfind.forcefield.relax_backbone_angles` already give a chain -- one angle per
    backbone atom of the repeat, broadcast along the oligomer, bounded to [95, 135] deg --
    applied inside the RIS scan instead of after it.  Each conformer's torsions are read back
    off the batch's coordinates, so the wrapper drops into :func:`fit_ris` unchanged.
    """

    def __init__(self, polymer, lo=95.0, hi=135.0, maxiter=60):
        self.polymer, self.lo, self.hi, self.maxiter = polymer, lo, hi, maxiter
        self.B = polymer.bonds_per_repeat
        self.a0 = np.array([polymer.backbone[k].backbone_angle for k in range(self.B)], dtype=float)
        self.n_relaxations = 0

    def _dihedrals(self, template, coords):
        from polyfind.chain import dihedral
        bb = template.backbone
        return np.array([[dihedral(c, bb[j], bb[j + 1], bb[j + 2], bb[j + 3]) for j in range(template.n_dihedrals)]
                         for c in coords])

    def energy_of(self, dih):
        """Angle-relaxed energy of one torsion vector."""
        dih = np.asarray(dih, dtype=float)
        h = 0.2

        def fg(a):
            rows = [a] + [a + s * h * np.eye(self.B)[j] for j in range(self.B) for s in (1.0, -1.0)]
            tpl, coords = build_chain_batch(self.polymer, np.repeat(dih[None], len(rows), axis=0),
                                            bond_angles=np.array(rows))
            E = VAL.energy_batch((tpl, coords))
            return float(E[0]), (E[1::2] - E[2::2]) / (2 * h)

        res = minimize(fg, self.a0, jac=True, method="L-BFGS-B", bounds=[(self.lo, self.hi)] * self.B,
                       options={"maxiter": self.maxiter, "ftol": 1e-9, "gtol": 1e-3})
        self.n_relaxations += 1
        return float(res.fun)

    def energy(self, struct):
        return self.energy_of(struct.dihedrals)

    def energy_batch(self, structs):
        if isinstance(structs, tuple):
            template, coords = structs
            dihs = self._dihedrals(template, np.asarray(coords))
        else:
            dihs = np.array([s.dihedrals for s in structs])
        return np.array([self.energy_of(d) for d in dihs])


def run_risfit(verbose=True):
    """The screen's RIS fit, rigid against angle-relaxed, for PVDF and the two nitriles.

    Three rows per polymer.  The two the document records: ``rigid`` (the frozen-angle scan)
    and ``angle-relaxed`` (this script's per-*type* wrapper, one angle per backbone atom of
    the repeat broadcast along the oligomer).  And the row that came out of them:
    ``fit_ris(angles="relaxed")``, the package's own per-*atom* relaxation with watershed
    basins, which is what VDCN and AN now get by default.  ``angles="rigid"`` is passed on
    the first two explicitly, since the default has moved.
    """
    out = []
    for name in ("pvdf", "vdcn", "an"):
        p = get_polymer(name)
        rec = {"polymer": name}
        for label, calc, angles in (("rigid", VAL, "rigid"), ("angle-relaxed", AngleRelaxedCalculator(p), "rigid"),
                                    ("fit_ris(angles='relaxed')", VAL, "relaxed")):
            t0 = time.time()
            fit = fit_ris(p, calc, step=10.0, n_monomers=6, third_order=True, angles=angles)
            m = fit.model
            cands = enumerate_periodic(p, m, max_period=8, k_per_period=40)
            e0 = cands[0].energy_per_monomer
            tt = "T" * p.bonds_per_repeat
            at = next((c for c in cands if c.name == tt), None)
            scans = {b: {"min_phi": float(fit.grid[np.argmin(fit.scan1[b])]), "min_E": float(np.min(fit.scan1[b])),
                         "E_at_180": float(fit.scan1[b][0]), "argmin": fit.argmin1[b],
                         "profile": [float(x) for x in fit.scan1[b]]} for b in fit.scan1}
            rec[label] = {"state_angles": [float(a) for a in m.states.angles],
                          "e1": np.round(m.first_order, 3).tolist(), "seconds": round(time.time() - t0, 1),
                          "ground_state": cands[0].name, "e_ground": e0,
                          "top": [(c.name, round(c.energy_per_monomer - e0, 3)) for c in cands[:6]],
                          "all_trans_rank": None if at is None else at.rank,
                          "all_trans_above_ground": None if at is None else at.energy_per_monomer - e0,
                          "scan1": scans}
            if verbose:
                r = rec[label]
                print(f"\n{name} RIS fit, {label} ({r['seconds']} s): state angles "
                      f"{tuple(round(a) for a in r['state_angles'])}; e1 {r['e1']}")
                for b, s in scans.items():
                    print(f"  1D scan bond type {b}: minimum {s['min_E']:+.2f} at {s['min_phi']:+.0f} deg; "
                          f"basin minima {s['argmin']}")
                print(f"  ground state {r['ground_state']} ({e0:+.3f}); all-trans rank {r['all_trans_rank']} "
                      f"at {r['all_trans_above_ground'] if at is None else round(r['all_trans_above_ground'], 3)} above")
                print("  top: " + ", ".join(f"{n} {d:+.2f}" for n, d in r["top"]))
        out.append(rec)
    return out


# --------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--part", choices=("profiles", "perturb", "correlation", "risfit", "all"), default="all")
    ap.add_argument("--quick", action="store_true", help="perturbation sweep with 3 amplitudes x 3 seeds")
    ap.add_argument("--from-relaxed", action="store_true",
                    help="perturbation sweep from each chain's own relaxed minimum (VDCN and PVDF all-trans only)")
    ap.add_argument("--step", type=float, default=10.0, help="profile grid step (deg)")
    ap.add_argument("--json", default=None, help="write every result to this file")
    args = ap.parse_args()
    t0 = time.time()
    results = {"preset": "pvdf-dft-valence", "torsion_triple": [float(v) for v in VAL.torsion]}
    if args.part in ("profiles", "all"):
        print("=" * 100 + "\nPART 1  relaxed dihedral profiles (kcal/mol for the oligomer; rigid level relative to rigid\n"
              "        all-trans, the other two relative to angle-relaxed all-trans)\n" + "=" * 100)
        results["profiles"] = run_profiles(step=args.step)
    if args.part in ("perturb", "all"):
        print("\n" + "=" * 100 + "\nPART 2  perturb the RIS ground state, relax the periodic chain (free-torsion refinement)\n" + "=" * 100)
        results["perturb"] = run_perturb(quick=args.quick, from_relaxed=args.from_relaxed)
    if args.part in ("correlation", "all"):
        print("\n" + "=" * 100 + "\nPART 3  torsional correlation length of the fitted RIS models\n" + "=" * 100)
        results["correlation"] = run_correlation()
    if args.part in ("risfit", "all"):
        print("\n" + "=" * 100 + "\nPART 4  the RIS fit with the backbone angles relaxed inside the scan\n" + "=" * 100)
        results["risfit"] = run_risfit()
    print(f"\ntotal {time.time() - t0:.0f} s")
    if args.json:
        with open(args.json, "w") as f:
            json.dump(results, f, indent=1)
        print(f"wrote {args.json}")


if __name__ == "__main__":
    main()
