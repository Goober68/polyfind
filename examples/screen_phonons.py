"""Gamma-point lattice curvature of the polar and antipolar all-trans cells of every screened chemistry.

    python examples/screen_phonons.py                          # every screenable chemistry and candidate
    python examples/screen_phonons.py --only pvdf,cfe          # a subset
    python examples/screen_phonons.py --json out.json          # keep the whole record

The governing requirement is low loss and a high-frequency response.  Neither side computes a
tan-delta; the two static quantities that bound it are the polar/antipolar lattice margin and the
lattice curvature.  ``docs/SCREEN.md`` measures the first and finds it inside the potential's
own error bar for most chemistries.  This script measures the second, on the same all-trans
polar zigzag and the same exactly-antipolar cell the screen compares, with the phonon
Hamiltonian of ``examples/crystal_phonons.py`` (``pvdf-dft-valence-flux-born``, induced
dipoles, Ewald, stretch minima pinned to the built bond lengths).  For each chemistry it reports:

* the softest optical modes of the polar cell once every atom is relaxed at fixed cell -- the
  stiffness of the transverse rigid-chain motions that switching and relaxation would use;
* whether each cell is a stationary point of the all-atom potential at all, or a saddle: the
  screen judges the *symmetric antipolar subspace* at ``gamma = 90``, which ``DESIGN.md``
  section 6 records is not where the model's antipolar minimum is (freeing ``gamma`` makes
  beta's polar and antipolar cells degenerate).  The Hessian sees that from inside the cell:
  the atoms slide when let go.

Four cells per chemistry, so that the screen's column and the model's own minima are both on
the table.  ``polar90`` and ``anti90`` are exactly the two cells ``docs/SCREEN.md`` compares
(``pack()`` and ``fitting.antipolar_cell_exact`` at ``gamma = 90``, plain Ewald packer).
``polar_free`` is the lowest polar cell of ``pack(gamma_free=True)``; ``anti_free`` is the
antipolar subspace polished with ``gamma`` free from several seeds
(:func:`polish_antipolar`).  Each cell is re-polished on the phonon packer -- over the cell
variables for the polar cells, over the subspace's own coordinates ``(a, b[, gamma], phi1, dz)``
with ``phi2 - phi1`` held for the antipolar ones -- and the atoms are then let go at fixed
cell.  The rigid-chain cell is not the deformable-path reference of ``docs/PHONONS.md`` (which
relaxes torsions and backbone angles as well), so beta-PVDF's ``polar90`` numbers here are a
construction check against that document, not a repeat of it.

Everything ``docs/PHONONS.md`` says is not in a Gamma-point curvature is not in this either: no
torsional stiffness (twist modes are lower bounds), no LO-TO splitting, no temperature,
no rates, no loss; and a Hessian at fixed cell cannot see a cell-shape instability, which is
why the gamma-free cells are on the table beside the ``gamma = 90`` ones.  The potential is
PVDF's: applying it to a chlorine or nitrile chemistry is the same transfer assumption the
screen makes, stated wherever a number is quoted.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
import traceback
import warnings

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import screen_electroactive as SE  # noqa: E402

from polyfind.ewald import EwaldSpec  # noqa: E402
from polyfind.fitting import antipolar_cell_exact  # noqa: E402
from polyfind.forcefield import SimpleFF, fit_ris  # noqa: E402
from polyfind.lattice_table import clear_pair_table_cache  # noqa: E402
from polyfind.phonon import (check_known_answer, force_terms, packer_with_built_bond_lengths,  # noqa: E402
                             phonons_gamma, placed_coordinates, relax_all_atom)
from polyfind.polarizability import Polarizable  # noqa: E402

PRESET_FLUX = "pvdf-dft-valence-flux-born"
POLAR_MIN_P = 1e-2  # C/m^2: below this a cell counts as non-polar, as in the screen
GAMMA_SEEDS = (65.0, 75.0, 85.0, 95.0, 105.0, 115.0)
CELL_KEYS = ["a", "b", "gamma", "phi1", "phi2", "dz"]


def _mode_row(m) -> dict:
    return {"cm1": round(float(m.freq_cm1), 2), "thz": round(float(m.freq_thz), 4),
            "transverse": round(m.transverse, 2), "axial": round(m.axial, 2),
            "rigid_chain": round(float(m.rigid_chain), 2),
            "axis": [round(float(v), 2) for v in m.axis],
            "by_element": {k: round(float(v), 2) for k, v in m.by_element.items()}}


def _softest(ph, want_transverse: bool, want_rigid: bool = False):
    """The lowest optical :class:`~polyfind.phonon.Mode` with the asked character, or ``None``."""
    for m in ph._modes_optical:
        if want_transverse and m.transverse < 0.5:
            continue
        if want_rigid and m.rigid_chain < 0.5:
            continue
        return m
    return None


def _lohi(chain, gamma_free: bool):
    from polyfind.pack import default_bounds

    b = default_bounds(chain, gamma_free=gamma_free)
    return np.array([b[k][0] for k in CELL_KEYS]), np.array([b[k][1] for k in CELL_KEYS])


def polish_antipolar(packer, params, lo, hi, gamma_free: bool = False, maxiter: int = 400) -> tuple:
    """``(params, E)``: L-BFGS over ``(a, b[, gamma], phi1, dz)`` with ``phi2 - phi1`` held.

    :func:`polyfind.fitting.antipolar_offsets` shows the cell dipole vanishes for every
    ``(a, b, gamma, dz)`` and every common rotation once ``phi2 = phi1 + dphi`` with the right
    flip, so these are the subspace's own coordinates: every point of the polish is exactly
    antipolar, on whatever packer it runs.  ``antipolar_cell_exact`` itself screens at
    ``gamma = 90`` only.
    """
    from scipy.optimize import minimize

    p0 = np.asarray(params, dtype=float).copy()
    dphi = p0[4] - p0[3]
    idx = [0, 1, 2, 3, 5] if gamma_free else [0, 1, 3, 5]
    c = float(packer.chain.c)

    def fg(x):
        p = p0.copy()
        p[idx] = x
        p[4] = p[3] + dphi
        # The packer's image bookkeeping assumes the chain-2 offset inside its own repeat; its
        # energy is periodic in dz there and nonsense far outside (a polish once walked VDCN's
        # dz nine repeats out and reported an energy 4 kcal/mol per monomer too low).  The
        # objective is evaluated on the wrapped offset; the gradient is unchanged by the wrap.
        p[5] = p[5] % c
        p[3], p[4] = p[3] % 360.0, p[4] % 360.0
        E, g, _, _ = packer.energy_and_grad(p)
        gi = np.asarray(g, dtype=float).copy()
        gi[3] = g[3] + g[4]
        return float(E), gi[idx]

    # ``antipolar_cell_exact``'s polish is free in (a, b) and may return a cell slightly outside
    # the screen's bounds; L-BFGS-B would project such a start onto the bound before its first
    # step, which is a move uphill, and a start *on* a bound with the gradient pointing out
    # cannot move at all.  The (a, b) bounds are widened by a full Angstrom around the start.
    lo = np.minimum(np.asarray(lo, dtype=float), p0[:6] - 1.0)
    hi = np.maximum(np.asarray(hi, dtype=float), p0[:6] + 1.0)
    bounds = [(lo[i], hi[i]) if i in (0, 1, 2) else (None, None) for i in idx]
    res = minimize(fg, p0[idx], jac=True, method="L-BFGS-B", bounds=bounds,
                   options={"ftol": 1e-10, "gtol": 1e-6, "maxiter": maxiter})
    p = p0.copy()
    p[idx] = res.x
    p[4] = p[3] + dphi
    p[3], p[4], p[5] = p[3] % 360.0, p[4] % 360.0, p[5] % c
    if p[5] > c - 1e-9:
        p[5] = 0.0
    return p, float(res.fun)


def antipolar_gamma_free(packer, anti90, lo, hi) -> tuple:
    """The lowest exactly-antipolar cell with ``gamma`` free, from the ``gamma = 90`` answer and seeds.

    A gradient polish from ``gamma = 90`` cannot leave it (it is a symmetric point), so the
    subspace is entered from :data:`GAMMA_SEEDS` on either side, at the ``gamma = 90`` axial
    registry and at the one half a repeat away.  Returns ``(params, E_per_monomer, seeds)``.
    """
    c = float(packer.chain.c)
    n_mon = packer.chain.n_monomers * packer.n_chains
    best_p, best_e = np.asarray(anti90, dtype=float).copy(), float(packer.energy(np.asarray(anti90, float)[None])[0])
    seeds = []
    for g0 in (90.0,) + GAMMA_SEEDS:
        for dz0 in (anti90[5], (anti90[5] + c / 2.0) % c):
            p0 = np.asarray(anti90, dtype=float).copy()
            p0[2], p0[5] = g0, dz0
            p, E = polish_antipolar(packer, p0, lo, hi, gamma_free=True)
            seeds.append({"gamma0": g0, "dz0": round(float(dz0), 3), "gamma": round(float(p[2]), 2),
                          "dz": round(float(p[5]), 3), "e_per_monomer": round(E / n_mon, 4)})
            if E < best_e - 1e-9:
                best_p, best_e = p, E
    return best_p, best_e / n_mon, seeds


def closest_contact(packer, params) -> float:
    """The shortest interatomic distance between different chains (lateral and axial images included), A.

    A sanity number for a packed cell: a polarizable point-dipole model can run into a
    polarization catastrophe when two polarizable atoms sit closer than the damping reaches,
    and a suspiciously dense cell is easier to judge with its closest contact beside it.
    """
    Pn, latn = placed_coordinates(packer, params)
    n, nc = packer.n, packer.n_chains
    best = np.inf
    for i in range(nc):
        Xi = Pn[i * n:(i + 1) * n]
        for j in range(nc):
            Xj = Pn[j * n:(j + 1) * n]
            for n1 in (-1, 0, 1):
                for n2 in (-1, 0, 1):
                    for n3 in (-1, 0, 1):
                        if i == j and (n1, n2) == (0, 0):
                            continue  # the chain's own axial images are bonded neighbours, not contacts
                        shift = n1 * latn[0] + n2 * latn[1] + n3 * latn[2]
                        d = np.linalg.norm(Xi[:, None, :] - (Xj + shift)[None, :, :], axis=-1)
                        best = min(best, float(d.min()))
    return best


def phonon_stage(full, plain, params, label: str, h: float, n_report: int = 12) -> dict:
    """Curvature of one cell on the pinned-stretch phonon packer: at the rigid-chain cell, then relaxed."""
    t0 = time.time()
    pk = packer_with_built_bond_lengths(full)
    e_cell, e_pack, rel = check_known_answer(pk, params)
    P0, latn = placed_coordinates(pk, params)
    N = pk.N
    q = np.asarray(plain._q_cell, dtype=float)[:N]  # the fixed charges, for a dipole diagnostic only
    n_mon = pk.chain.n_monomers * pk.n_chains
    mu0 = q @ P0
    ft = force_terms(pk, params, P0, latn)
    ph0 = phonons_gamma(pk, params, h=h, Pn=P0, latn=latn, asr="project", n_report=n_report)
    Pn, E, gmax = relax_all_atom(pk, params, Pn=P0, latn=latn)
    mu1 = q @ Pn
    # where the atoms went, chain by chain: the mean displacement of each chain (a rigid slide or
    # shift) and the largest internal displacement once that mean is removed
    n = pk.n
    slides, internal = [], []
    for s in range(pk.n_chains):
        d = (Pn - P0)[s * n:(s + 1) * n]
        m = d.mean(axis=0)
        slides.append([round(float(v), 4) for v in m])
        internal.append(round(float(np.abs(d - m).max()), 4))
    ph1 = phonons_gamma(pk, params, h=h, Pn=Pn, latn=latn, asr="project", n_report=n_report)
    soft_t = _softest(ph1, True)
    soft_rc = _softest(ph1, True, True)
    rec = {"label": label, "atoms": int(N), "known_answer_rel": float(f"{rel:.2e}"),
           "params": [round(float(v), 4) for v in params],
           "rigid": {"energy": round(float(e_pack), 6), "max_force": float(f"{ft['total']:.3e}"),
                     "bond_strain": round(float(ft["bond_strain"]), 4), "angle_strain_deg": round(float(ft["angle_strain"]), 2),
                     "n_imaginary": int(ph0.n_imaginary), "n_zero": int(ph0.n_zero),
                     "lowest_optical": [round(float(v), 2) for v in ph0.optical[:6]],
                     "modes": [_mode_row(m) for m in ph0._modes_optical[:3]],
                     "dipole_per_monomer": round(float(np.linalg.norm(mu0)) / n_mon, 4)},
           "relaxed": {"energy": round(float(E), 6), "lowered_by": round(float(e_pack - E), 6),
                       "max_force": float(f"{gmax:.1e}"), "max_displacement": round(float(np.abs(Pn - P0).max()), 4),
                       "chain_slides": slides, "max_internal_displacement": internal,
                       "n_imaginary": int(ph1.n_imaginary), "n_zero": int(ph1.n_zero),
                       "asr_max_relative": float(f"{ph1.asr_max / ph1.hessian_scale:.1e}"),
                       "lowest_optical": [round(float(v), 2) for v in ph1.optical[:6]],
                       "modes": [_mode_row(m) for m in ph1._modes_optical[:3]],
                       "softest_transverse": None if soft_t is None else _mode_row(soft_t),
                       "softest_rigid_chain_transverse": None if soft_rc is None else _mode_row(soft_rc),
                       "dipole_per_monomer": round(float(np.linalg.norm(mu1)) / n_mon, 4)},
           "seconds": round(time.time() - t0, 1)}
    print(f"    [{label}] {N} atoms; known-answer {rel:.1e}; rigid cell: max force {ft['total']:.2e}, "
          f"{ph0.n_imaginary} imaginary, lowest optical {' '.join(f'{v:.1f}' for v in ph0.optical[:4])} cm^-1; "
          f"|mu|/mon {rec['rigid']['dipole_per_monomer']:.4f} e.A")
    print(f"    [{label}] relaxed at fixed cell: lowered {e_pack - E:.4f} kcal/mol, max displacement "
          f"{np.abs(Pn - P0).max():.4f} A, {ph1.n_imaginary} imaginary; |mu|/mon {rec['relaxed']['dipole_per_monomer']:.4f} e.A; "
          f"lowest optical {' '.join(f'{v:.1f}' for v in ph1.optical[:6])} cm^-1  ({rec['seconds']} s)")
    print(f"    [{label}] chain mean displacements (x, y, z) A: " + "; ".join(str(v) for v in slides)
          + f"; largest internal displacement after removing them: {max(internal):.4f} A")
    for m in ph1._modes_optical[:3]:
        print("      " + m.describe())
    return rec


def one_chemistry(polymer, h: float, anti_target: int, anti_polish: int, anti_maxfev: int,
                  fit_step: float = 10.0, fit_monomers: int = 6) -> dict:
    from polyfind import pack as pack_mod
    from polyfind.pack import polish

    FV, VAL, _ = SE._presets()
    FLX = SimpleFF.from_preset(PRESET_FLUX)
    rec = {"polymer": polymer.name, "formula": polymer.formula, "atoms_per_repeat": polymer.atoms_per_repeat,
           "preset": PRESET_FLUX, "timings": {}, "notes": []}
    print("\n" + "=" * 100)
    print(f"{polymer.name.upper()}  {polymer.formula}")
    print("=" * 100)
    t_all = time.time()

    # --- the all-trans chain, built as the screen builds it --------------------------------
    t0 = time.time()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        fit = fit_ris(polymer, VAL, step=fit_step, n_monomers=fit_monomers, third_order=True)
    for w in caught:
        msg = str(w.message).split("\n")[0]
        if msg not in rec["notes"]:
            rec["notes"].append(msg)
    rec["timings"]["fit"] = round(time.time() - t0, 1)
    seq = SE.all_trans_seq(polymer)
    chain, how = SE.build_repeat(polymer, seq, fit.model.states, helix=None)
    rec["build"] = how
    if chain is None:
        rec["screened"] = False
        rec["reason"] = f"no all-trans repeat: {how}"
        print(f"  NOT MEASURED: {rec['reason']}")
        return rec
    print(f"  all-trans repeat: {how}; c = {chain.c:.4f} A, {chain.n_atoms} atoms per repeat "
          f"(RIS fit {rec['timings']['fit']} s, {fit.angles} angles)")
    n_mon = chain.n_monomers * 2

    spec = EwaldSpec()
    clear_pair_table_cache()
    with FV.applied():
        plain = pack_mod.CrystalPacker(chain, n_chains=2, coulomb="ewald", ewald=spec)
        screen_pk = pack_mod.CrystalPacker(chain, n_chains=2)
        lo90, hi90 = _lohi(chain, False)
        lof, hif = _lohi(chain, True)

        # --- (1) the two cells the screen compares: gamma = 90, plain Ewald packer ----------
        t0 = time.time()
        results = pack_mod.pack(chain, n_chains=2, table_cache_dir=None, coulomb="ewald", ewald=spec)
        best90 = results[0]
        polar90 = next((r for r in results if r.polarization_magnitude > POLAR_MIN_P), None)
        rec["timings"]["pack"] = round(time.time() - t0, 1)
        t0 = time.time()
        anti90, anti90_e, anti90_pol = antipolar_cell_exact(plain, screen_packer=screen_pk, target=anti_target,
                                                            n_polish=anti_polish, maxfev=anti_maxfev)
        rec["timings"]["antipolar"] = round(time.time() - t0, 1)
        gap90 = float(anti90_e - best90.energy_per_monomer)
        print(f"  gamma = 90 (the screen's column): best {best90.row()}")
        if polar90 is not None and polar90 is not best90:
            print(f"                                 polar {polar90.row()}")
        print(f"                                 antipolar E/mon {anti90_e:+.4f} (max|P| {anti90_pol:.1e}); antipolar minus "
              f"best {gap90:+.3f} +/- {SE.ERROR_BAR:.2f} kcal/mol per monomer -> {SE.verdict(gap90).upper()}  "
              f"(pack {rec['timings']['pack']} s, antipolar {rec['timings']['antipolar']} s)")

        # --- (2) the same two branches with gamma free: the model's own minima --------------
        t0 = time.time()
        results_f = pack_mod.pack(chain, n_chains=2, table_cache_dir=None, coulomb="ewald", ewald=spec, gamma_free=True)
        bestf = results_f[0]
        polarf = next((r for r in results_f if r.polarization_magnitude > POLAR_MIN_P), None)
        antif, antif_e, seeds = antipolar_gamma_free(plain, anti90, lof, hif)
        ra_f = plain.result(antif)
        rec["timings"]["gamma_free"] = round(time.time() - t0, 1)
        gapf = float(antif_e - bestf.energy_per_monomer)
        print(f"  gamma free: best {bestf.row()}")
        if polarf is not None and polarf is not bestf:
            print(f"              polar {polarf.row()}")
        print(f"              antipolar E/mon {antif_e:+.4f} at gamma {antif[2]:.1f} dz {antif[5]:.3f} (max|P| "
              f"{ra_f.polarization_magnitude:.1e}; from gamma = 90's {anti90_e:+.4f}); antipolar minus best "
              f"{gapf:+.3f} -> {SE.verdict(gapf).upper()}  ({rec['timings']['gamma_free']} s)")
        rec["screen"] = {
            "gamma90": {"best": best90.row(), "best_e_per_monomer": round(float(best90.energy_per_monomer), 4),
                        "best_P": round(best90.polarization_magnitude, 4),
                        "polar": None if polar90 is None else polar90.row(),
                        "polar_e_per_monomer": None if polar90 is None else round(float(polar90.energy_per_monomer), 4),
                        "antipolar_e_per_monomer": round(float(anti90_e), 4), "antipolar_max_P": float(f"{anti90_pol:.2e}"),
                        "antipolar_params": [round(float(v), 4) for v in anti90],
                        "antipolar_gap": round(gap90, 4), "arrangement": SE.verdict(gap90)},
            "gamma_free": {"best": bestf.row(), "best_e_per_monomer": round(float(bestf.energy_per_monomer), 4),
                           "best_P": round(bestf.polarization_magnitude, 4),
                           "polar": None if polarf is None else polarf.row(),
                           "polar_e_per_monomer": None if polarf is None else round(float(polarf.energy_per_monomer), 4),
                           "antipolar_e_per_monomer": round(float(antif_e), 4),
                           "antipolar_max_P": float(f"{ra_f.polarization_magnitude:.2e}"),
                           "antipolar_params": [round(float(v), 4) for v in antif], "antipolar_seeds": seeds,
                           "antipolar_gap": round(gapf, 4), "arrangement": SE.verdict(gapf),
                           "antipolar_density": round(float(ra_f.density), 4),
                           "antipolar_closest_contact": round(closest_contact(plain, antif), 4)},
            "error_bar": SE.ERROR_BAR}
        print(f"              gamma-free antipolar cell: density {ra_f.density:.3f} g/cm^3, closest interchain contact "
              f"{rec['screen']['gamma_free']['antipolar_closest_contact']:.3f} A")

        # --- (3) each cell re-polished on the phonon packer ----------------------------------
        full = pack_mod.CrystalPacker(chain, n_chains=2, valence=VAL, charge_flux=FLX,
                                      polarizable=Polarizable(), coulomb="ewald", ewald=spec)
        cells = []

        def sane(label, rp, p) -> bool:
            """False, with a note, when the phonon packer cannot evaluate the cell (polarization catastrophe)."""
            e, P = float(rp.energy_per_monomer), float(rp.polarization_magnitude)
            if np.isfinite(e) and abs(e) < 1e4 and np.isfinite(P) and P < 10.0:
                return True
            rec["notes"].append(f"{label}: the phonon packer's induced-dipole solve diverges at this cell "
                                f"(E/mon {e:.3g}, |P| {P:.3g}; closest interchain contact {closest_contact(plain, p):.3f} A): "
                                "not measured")
            print("    " + rec["notes"][-1])
            return False

        def polar_cell(label, r, lo, hi):
            if r is None:
                rec["notes"].append(f"{label}: pack() found no polar cell")
                return
            free = [i for i in range(6) if hi[i] > lo[i]]
            p = polish(full, r.params, lo, hi, free, maxiter=300)
            rp = full.result(p)
            e0 = float(full.energy_per_monomer(r.params[None])[0])
            print(f"  {label} on the phonon packer: E/mon {e0:+.4f} -> {rp.energy_per_monomer:+.4f}, a {r.a:.3f}->{rp.a:.3f} "
                  f"b {r.b:.3f}->{rp.b:.3f} gamma {r.gamma:.1f}->{rp.gamma:.1f}, |P| {rp.polarization_magnitude:.3f}")
            if not sane(label, rp, p):
                return
            if rp.polarization_magnitude > POLAR_MIN_P:
                cells.append((label, p, rp))
            else:
                rec["notes"].append(f"{label}: lost its polarization under the phonon packer's cell polish")
                print(f"    {label} is no longer polar after the polish: not measured")

        def anti_cell(label, p0, lo, hi, gamma_free):
            # The subspace is the screen's: phi2 - phi1 fixed by the *fixed* charges' moment
            # (antipolar_offsets on the plain packer).  Holding it keeps the fixed-charge dipole
            # exactly zero on any packer; the phonon packer's flux charges tilt the moment
            # slightly and its induced dipoles need not cancel unless the arrangement is a
            # crystal symmetry (it is for a moment along a symmetry axis, PVDF; it is not for
            # an off-axis moment, CDFE), so its residual |P| is recorded, not asserted on.
            p0 = np.asarray(p0, dtype=float)
            if not sane(label, full.result(p0), p0):
                return
            p, E = polish_antipolar(full, p0, lo, hi, gamma_free=gamma_free)
            rp = full.result(p)
            e0 = float(full.energy_per_monomer(p0[None])[0])
            e_chk = float(full.energy(p[None])[0])
            print(f"  {label} on the phonon packer (phi2 - phi1 held): E/mon {e0:+.4f} -> {rp.energy_per_monomer:+.4f}, "
                  f"a {p0[0]:.3f}->{rp.a:.3f} b {p0[1]:.3f}->{rp.b:.3f} gamma {p0[2]:.1f}->{rp.gamma:.1f} "
                  f"dz {p0[5]:.3f}->{rp.dz:.3f}; |P| fixed charges {plain.result(p).polarization_magnitude:.1e}, "
                  f"phonon packer (flux + induced) {rp.polarization_magnitude:.1e}")
            if plain.result(p).polarization_magnitude > 1e-9:
                raise RuntimeError(f"{label}: tied polish left the fixed-charge antipolar subspace: "
                                   f"|P| = {plain.result(p).polarization_magnitude:.2e}")
            if abs(e_chk - E) > 1e-6 * max(1.0, abs(E)) or rp.energy_per_monomer > e0 + 1e-6:
                rec["notes"].append(f"{label}: tied polish did not lower the energy consistently (objective {E:.6f}, "
                                    f"energy {e_chk:.6f}, per monomer {e0:+.4f} -> {rp.energy_per_monomer:+.4f}); "
                                    "the rigid start is used instead")
                print(f"    {label}: polish inconsistent, keeping the start: " + rec["notes"][-1])
                p, rp = p0.copy(), full.result(p0)
            cells.append((label, p, rp))

        polar_cell("polar90", polar90, lo90, hi90)
        anti_cell("anti90", anti90, lo90, hi90, False)
        polar_cell("polar_free", polarf, lof, hif)
        anti_cell("anti_free", antif, lof, hif, True)
        rec["phonon_packer"] = {c[0]: {"row": c[2].row(), "e_per_monomer": round(float(c[2].energy_per_monomer), 4),
                                       "P": float(f"{c[2].polarization_magnitude:.3e}"), "density": round(float(c[2].density), 4),
                                       "closest_contact": round(closest_contact(plain, c[1]), 4),
                                       "params": [round(float(v), 4) for v in c[1]]} for c in cells}
        for c in cells:
            print(f"  {c[0]}: density {c[2].density:.3f} g/cm^3, closest interchain contact "
                  f"{rec['phonon_packer'][c[0]]['closest_contact']:.3f} A")

        # --- (4) curvature -------------------------------------------------------------------
        rec["cells"] = {}
        for label, p, _ in cells:
            try:
                rec["cells"][label] = phonon_stage(full, plain, p, label, h)
            except Exception as e:  # noqa: BLE001 - one cell failing must not lose the chemistry
                rec["cells"][label] = {"label": label, "error": f"{type(e).__name__}: {e}"}
                print(f"    [{label}] FAILED: {type(e).__name__}: {e}")
                traceback.print_exc()
    rec["n_monomers_per_cell"] = int(n_mon)
    rec["screened"] = True
    rec["timings"]["total"] = round(time.time() - t_all, 1)
    return rec


def relaxed_gap(r: dict, polar: str, anti: str):
    """Antipolar minus polar, kcal/mol per monomer, on the phonon packer with every atom relaxed at fixed cell."""
    cp, ca = r.get("cells", {}).get(polar, {}), r.get("cells", {}).get(anti, {})
    if "relaxed" not in cp or "relaxed" not in ca:
        return None
    return (ca["relaxed"]["energy"] - cp["relaxed"]["energy"]) / r["n_monomers_per_cell"]


def _low3(c: dict) -> str:
    if "relaxed" not in c:
        return c.get("error", "-")[:22]
    return " ".join(f"{v:.1f}" for v in c["relaxed"]["lowest_optical"][:3])


def _slid(c: dict) -> str:
    return "-" if "relaxed" not in c else f"{c['relaxed']['max_displacement']:.2f}"


def summary_table(records: list) -> str:
    hdr = (f"{'chemistry':<11s} {'gap90':>7s} {'gapfree':>8s} {'verdict(free)':<13s} {'gapRelax':>8s}  "
           f"{'polar90 lowest 3':<22s} {'slid':>5s}  {'polar_free lowest 3':<22s} {'slid':>5s}  "
           f"{'anti90 slid':>11s}  {'anti_free lowest 3':<22s} {'slid':>5s}")
    lines = [hdr]
    for r in records:
        if not r.get("screened"):
            lines.append(f"{r['polymer']:<11s} not measured: {r.get('reason', '')[:110]}")
            continue
        s90, sf = r["screen"]["gamma90"], r["screen"]["gamma_free"]
        c = r["cells"]
        g = relaxed_gap(r, "polar_free", "anti_free")
        lines.append(f"{r['polymer']:<11s} {s90['antipolar_gap']:+7.3f} {sf['antipolar_gap']:+8.3f} {sf['arrangement']:<13s} "
                     f"{'-' if g is None else f'{g:+.3f}':>8s}  "
                     f"{_low3(c.get('polar90', {})):<22s} {_slid(c.get('polar90', {})):>5s}  "
                     f"{_low3(c.get('polar_free', {})):<22s} {_slid(c.get('polar_free', {})):>5s}  "
                     f"{_slid(c.get('anti90', {})):>11s}  "
                     f"{_low3(c.get('anti_free', {})):<22s} {_slid(c.get('anti_free', {})):>5s}")
    return "\n".join(lines)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--only", default=None, help="comma-separated polymer names")
    ap.add_argument("--h", type=float, default=1e-3, help="finite-difference step, A")
    ap.add_argument("--anti-target", type=int, default=6000)
    ap.add_argument("--anti-polish", type=int, default=8)
    ap.add_argument("--anti-maxfev", type=int, default=900)
    ap.add_argument("--no-candidates", action="store_true", help="skip the screen's unregistered candidates")
    ap.add_argument("--json", default=None)
    args = ap.parse_args(argv)
    only = set(args.only.split(",")) if args.only else None
    todo = []
    for polymer, skip in SE.incumbents():
        if only and polymer.name not in only:
            continue
        todo.append((polymer, skip))
    if not args.no_candidates:
        for polymer in SE.new_candidates():
            if only and polymer.name not in only:
                continue
            todo.append((polymer, None))
    print(f"phonon screen of {len(todo)} chemistries; preset {PRESET_FLUX} + induced dipoles + Ewald; h = {args.h:g} A")
    records = []
    for polymer, skip in todo:
        if skip:
            records.append({"polymer": polymer.name, "formula": polymer.formula, "screened": False, "reason": skip})
            print(f"\n{polymer.name.upper()}: NOT MEASURED: {skip[:120]}...")
            continue
        try:
            records.append(one_chemistry(polymer, args.h, args.anti_target, args.anti_polish, args.anti_maxfev))
        except Exception as e:  # noqa: BLE001
            records.append({"polymer": polymer.name, "formula": polymer.formula, "screened": False,
                            "reason": f"{type(e).__name__}: {e}"})
            print(f"  FAILED: {type(e).__name__}: {e}")
            traceback.print_exc()
        if args.json:
            with open(args.json, "w") as f:
                json.dump({"preset": PRESET_FLUX, "h": args.h, "records": records}, f, indent=1)
    print("\n" + "=" * 100)
    print("SUMMARY  gap90 = antipolar minus best at gamma = 90 (docs/SCREEN.md's column), gapfree = the same with")
    print("         gamma free, both on the screen packer, kcal/mol per monomer; gapRelax = anti_free minus")
    print("         polar_free on the phonon packer after every atom relaxed at fixed cell.  Curvature columns:")
    print("         lowest three optical modes (cm^-1) of that relaxed state, stretch minima pinned; 'slid' = the")
    print("         largest atom displacement (A) when the rigid cell's atoms were let go: a few hundredths for a")
    print("         cell that was already a minimum, an Angstrom for a saddle.  The potential is PVDF's throughout.")
    print("=" * 100)
    print(summary_table(records))
    return 0


if __name__ == "__main__":
    sys.exit(main())
