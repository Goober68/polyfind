"""Staggered eight-chain 11 VDF : 1 VDCN starts: break the axial nitrile registry.

The cells shipped by ``copolymer_starts.py`` are a 2x2x1 tiling of a twelve-monomer
two-chain repeat, so every one of the eight chains carries its VDCN unit at the *same*
axial height.  A mass-fraction mixing rule over the two homopolymers predicts 2.028 g/cm3
for this composition and the tiled cell reaches 1.682, about 17% less dense, and that was
shown to be a construction artifact rather than a search failure (six fixed-shape starts at
about 1.96 all polished back to the same cell).  A real 8.33 mol% solid would stagger the
units between chains.

This script builds the staggered variants, relaxes them under the same potential, checks
their topology at least as strictly as before, and writes them into ``deliverables/``.
Three things make it more than a coordinate shuffle:

* the chain is measured, not assumed, to be monomer-periodic in its skeleton, so a
  whole-monomer slide moves only the comonomer decoration (:func:`skeleton_periodicity`);
* the eight chains are relaxed as eight *independent* columns
  (:class:`polyfind.supercell.SupercellEnergy`), because a tiling cannot relax a stagger --
  and that evaluator is checked against ``CrystalPacker`` on the packer's own cell, so the
  staggered cells and the aligned ones are scored by one function;
* the ``aligned`` pattern is run through the identical search as a control, so "the
  staggered cell is lower" is a statement about the stagger and not about the optimiser.

Run with ``PYTHONPATH=src python examples/copolymer_staggered.py`` after
``copolymer_starts.py`` (it reads the aligned cells from ``deliverables/summary.json``).
"""
from __future__ import annotations

import itertools
import json
import os
import sys
import time

import numpy as np
from scipy.optimize import minimize

from polyfind.ewald import EwaldSpec
from polyfind.pack import CrystalPacker, periodic_chain
from polyfind.polymers import THREE_STATE, VDF_VDCN_11_1 as CP
from polyfind.supercell import (
    STAGGER_PATTERNS, SupercellEnergy, cell_density, cell_dipole, cell_polarization,
    cell_to_cif, column_layout, stagger_pattern, staggered_cell,
)
from polyfind.topology import build_cell, check_topology

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, "deliverables")
NX = NY = 2
PATTERNS = ("aligned", "antiphase", "sheet", "ladder", "spread")
SHIP = ("ladder", "spread")
# The fixed shapes the screen pins (a, b) at.  The first is filled in from the aligned cell;
# the rest run from its density up past the consumer's accepted 1.96 g/cm3, because the whole
# question is whether a *staggered* cell started dense stays dense.
SHAPES = ((9.6, 4.8), (9.0, 4.8), (8.6, 4.8), (8.2, 4.8), (9.0, 5.2), (8.6, 5.6),
          (10.0, 4.4), (7.8, 5.0))
# A band, not a point.  The five the aligned deliverable reported, plus 1.25 and the two
# coarser ones, so the new cells are judged at least as strictly as the old ones.
SCALES = (1.05, 1.10, 1.15, 1.20, 1.25, 1.30, 1.35)
MIXING_RULE_DENSITY = 2.028  # mass-fraction rule over the two homopolymer cells
CLOSE_COLUMN_A = 7.0  # column-axis separation below which two chains are "close neighbours"


def log(*a):
    print(*a, flush=True)


# ------------------------------------------------------------ 1. is a slide a slide?
def skeleton_periodicity(ch):
    """Which atoms a one-monomer axial slide maps onto an identical atom, and how exactly.

    The stagger is an integer number of monomers because the copolymer's backbone *is* a
    PVDF backbone: sliding by ``c / 12`` should map the whole carbon/hydrogen/fluorine
    skeleton onto itself and move only the six atoms that distinguish the comonomer site.
    Should, not does -- if the cyano carbon carried a different backbone angle the skeleton
    would not be monomer-periodic at all and a slide would be a distortion.  So it is
    measured.
    """
    X, c, n = ch.coords, ch.c, ch.n_atoms
    els = list(ch.elements)
    shifted = X + np.array([0.0, 0.0, c / ch.n_monomers])
    matched, moved, worst = 0, [], 0.0
    for i in range(n):
        best = (np.inf, -1)
        for k in (-1, 0, 1):
            d = np.linalg.norm(X - (shifted[i] + np.array([0.0, 0.0, k * c])), axis=1)
            j = int(np.argmin(d))
            if d[j] < best[0]:
                best = (float(d[j]), j)
        if best[0] < 1e-8 and els[best[1]] == els[i]:
            matched += 1
            worst = max(worst, best[0])
        else:
            moved.append((i, els[i], els[best[1]], round(best[0], 4)))
    return {"n_atoms": n, "n_mapped_exactly": matched,
            "worst_displacement_A": worst,
            "atoms_that_move": moved}


# --------------------------------------------------- 2. which stagger, and how separated
def column_offsets(a: float, b: float, gamma: float, nx=NX, ny=NY, k_chains=2) -> np.ndarray:
    """Lateral position of each chain's axis, in the order ``build_cell`` emits."""
    g = np.deg2rad(gamma)
    av, bv = np.array([a, 0.0]), np.array([b * np.cos(g), b * np.sin(g)])
    return np.array([ix * av + iy * bv + k * 0.5 * (av + bv)
                     for ix, iy, k in column_layout(nx, ny, k_chains)])


def close_column_pairs(a, b, gamma, cut=CLOSE_COLUMN_A, nx=NX, ny=NY):
    """``(t, u, separation)`` for every chain pair whose axes come within ``cut``.

    Distances are minimised over the *supercell's* lateral images, so a chain's neighbour
    across the cell boundary counts, which is the whole point: the aligned cell's problem is
    a plane of nitriles and a plane does not stop at the cell edge.
    """
    o = column_offsets(a, b, gamma, nx, ny)
    g = np.deg2rad(gamma)
    A = np.array([nx * a, 0.0])
    B = np.array([ny * b * np.cos(g), ny * b * np.sin(g)])
    out = []
    for t, u in itertools.combinations_with_replacement(range(len(o)), 2):
        d = min(np.linalg.norm(o[u] - o[t] + i * A + j * B)
                for i in (-1, 0, 1) for j in (-1, 0, 1)
                if not (t == u and i == 0 and j == 0))
        if d < cut:
            out.append((t, u, float(d)))
    return sorted(out, key=lambda r: r[2])


def stagger_quality(offsets, pairs, monomers=12):
    """How far apart a pattern keeps the comonomer units of neighbouring chains.

    ``min_offset`` is the smallest circular monomer offset over every close column pair --
    the quantity a stagger is for.  ``n_heights`` is how many of the ``monomers`` axial slots
    the cell's chains occupy, which is the complementary measure: a pattern can keep every
    near pair well separated and still put several chains at one height.
    """
    m = np.asarray(offsets, int)
    circ = lambda d: min(d % monomers, (-d) % monomers)  # noqa: E731
    gaps = [circ(int(m[u]) - int(m[t])) for t, u, _ in pairs]
    return {"min_offset_monomers": int(min(gaps)) if gaps else monomers,
            "mean_offset_monomers": float(np.mean(gaps)) if gaps else float(monomers),
            "n_distinct_heights": int(len(set(m.tolist()))),
            "offsets": m.tolist()}


def best_possible_min_offset(pairs, n_chains, monomers=12):
    """The largest ``min_offset`` *any* assignment of integer offsets can reach.

    Backtracking over the close-pair constraint graph, largest target first.  Reported so
    that "maximally separated" is a measured claim about this lateral geometry rather than an
    assertion: for the centred cell here every A-B pair is a close pair, which caps it.
    """
    circ = lambda d: min(d % monomers, (-d) % monomers)  # noqa: E731
    adj = [[] for _ in range(n_chains)]
    for t, u, _ in pairs:
        if t != u:
            adj[t].append(u)
            adj[u].append(t)

    def feasible(d):
        m = [None] * n_chains

        def go(i):
            if i == n_chains:
                return True
            for v in range(monomers):
                if all(m[j] is None or circ(v - m[j]) >= d for j in adj[i]):
                    m[i] = v
                    if go(i + 1):
                        return True
                    m[i] = None
            return False

        return go(0)

    for d in range(monomers // 2, 0, -1):
        if feasible(d):
            return d
    return 0


# --------------------------------------------------------------- 3. the eight-chain relax
class Objective:
    """The eight-chain lattice energy as a function of the free cell variables.

    ``antipolar`` ties ``phi2 = phi1 + 180``, which is the exactly-antipolar branch
    ``polyfind.fitting.antipolar_offsets`` derives for this chain (its moment has no axial
    component, so both ``flip`` branches are exact and the shipped cell took ``flip = 0``).
    Holding the constraint rather than penalising it is what keeps the cell antipolar to
    floating point at every step of the search, not only at the end.
    """

    ZMAX = 0.5  # monomers: how far a per-chain slide may wander off the chosen pattern

    def __init__(self, se, ch, pk, stagger, antipolar: bool, bounds, free_z: bool = False):
        self.se, self.ch, self.pk = se, ch, pk
        self.stagger = np.asarray(stagger, int)
        self.antipolar = antipolar
        self.lo, self.hi = bounds
        self.n_cell = len(self.stagger)
        self.n_base = 4 if antipolar else 5
        self.free_z = free_z
        self.i_dz = self.n_base - 1
        self.n_calls = 0

    def params(self, x):
        if self.antipolar:
            a, b, phi1, dz = x[:4]
            phi2 = phi1 + 180.0
        else:
            a, b, phi1, phi2, dz = x[:5]
        return np.array([a, b, 90.0, phi1, phi2, dz, 0.0], float)

    def offsets(self, x):
        """The per-chain axial offset in monomers: the integer pattern plus any free slide.

        Chain 0 is held, because a rigid translation of the whole cell along ``c`` is a
        symmetry of the energy and leaving it free would only make the optimiser wander.
        """
        m = self.stagger.astype(float)
        if self.free_z:
            m = m + np.concatenate([[0.0], np.asarray(x[self.n_base :], float)])
        return m

    def x0(self, a, b, phi1, phi2, dz, z=None):
        base = [a, b, phi1, dz] if self.antipolar else [a, b, phi1, phi2, dz]
        if not self.free_z:
            return np.array(base, float)
        return np.concatenate([base, np.zeros(self.n_cell - 1) if z is None else np.asarray(z, float)])

    def cell(self, x):
        return staggered_cell(self.ch, self.params(x), self.offsets(x), nx=NX, ny=NY, packer=self.pk)

    def __call__(self, x):
        a, b = x[0], x[1]
        if not (self.lo[0] <= a <= self.hi[0] and self.lo[1] <= b <= self.hi[1]):
            return 1e6
        if self.free_z and np.any(np.abs(x[self.n_base :]) > self.ZMAX):
            return 1e6
        self.n_calls += 1
        return self.se.energy(self.cell(x))

    def per_monomer(self, x):
        return self(x) / (self.n_cell * self.ch.n_monomers)


def screen(obj, ch, seed_params, rng, n_random=360, n_shape=40):
    """Candidate starts: the aligned optimum's registries, dense fixed shapes, and random.

    The dense shapes are the point of the exercise.  The aligned cell's density deficit was
    attributed to the nitrile plane, so the test is whether a staggered cell *started dense*
    stays dense; a screen that only looks near the loose aligned optimum could not answer
    that however many samples it took.
    """
    c = ch.c
    a0, b0, phi1_0, phi2_0, dz0 = seed_params
    rows, groups = [], []
    # (i) the aligned optimum, at all twelve axial registries of the chain pair
    start = 0
    for j in range(ch.n_monomers):
        rows.append(obj.x0(a0, b0, phi1_0, phi2_0, (dz0 + j * c / ch.n_monomers) % c))
    groups.append(("aligned optimum, 12 registries", a0, b0, start, len(rows)))
    # (ii) fixed shapes spanning the aligned density up past the consumer's 1.96, setting
    # angles and dz sampled densely at each
    for a, b in ((a0, b0),) + SHAPES:
        start = len(rows)
        for _ in range(n_shape):
            rows.append(obj.x0(a, b, rng.uniform(0, 360), rng.uniform(0, 360), rng.uniform(0, c)))
        groups.append((f"fixed shape {a:.2f} x {b:.2f}", a, b, start, len(rows)))
    # (iii) a plain random screen, as insurance against a qualitatively different packing
    start = len(rows)
    for _ in range(n_random):
        rows.append(obj.x0(rng.uniform(6.0, 14.0), rng.uniform(4.0, 8.0),
                           rng.uniform(0, 360), rng.uniform(0, 360), rng.uniform(0, c)))
    groups.append(("uniform random", 0.0, 0.0, start, len(rows)))
    return np.array(rows), groups


def distinct(rows, energies, n, i_dz, sep_ab=0.5, sep_dz=1.0):
    """The ``n`` lowest starts that are not near-duplicates of one already taken."""
    taken = []
    for i in np.argsort(energies):
        x = rows[i]
        if any(abs(x[0] - t[0]) + abs(x[1] - t[1]) < sep_ab and abs(x[i_dz] - t[i_dz]) < sep_dz
               for t in taken):
            continue
        taken.append(x)
        if len(taken) >= n:
            break
    return taken


def polish(obj, x0, maxfev=260):
    r = minimize(obj, np.asarray(x0, float), method="Nelder-Mead",
                 options={"maxfev": maxfev, "xatol": 1e-3, "fatol": 1e-5})
    return np.asarray(r.x, float), float(r.fun)


def registry_sweep(obj, ch, x, e, n_keep=3):
    """Re-polish the twelve axial registries of a polished cell.

    ``dz`` is the one variable a local polish cannot cross: the chain is monomer-periodic in
    its skeleton, so ``dz`` and ``dz + c/12`` are two different comonomer registries with the
    same backbone packing, separated by a barrier the optimiser will not climb.  The aligned
    deliverable swept them; so does this, for every pattern.
    """
    c = ch.c
    i = obj.i_dz
    grid = np.repeat(np.asarray(x, float)[None], ch.n_monomers, axis=0)
    grid[:, i] = (x[i] + np.arange(ch.n_monomers) * c / ch.n_monomers) % c
    ee = np.array([obj(g) for g in grid])
    per = obj.n_cell * ch.n_monomers
    log(f"      registry sweep: twelve dz offsets span {ee.min() / per:+.4f} to {ee.max() / per:+.4f}")
    out, out_e = np.asarray(x, float), e
    for k in np.argsort(ee)[:n_keep]:
        p, pe = polish(obj, grid[k])
        if pe < out_e:
            out, out_e = p, pe
    return out, out_e


def relax(se, ch, pk, stagger, antipolar, seed_params, bounds, rng, label, free_z=True):
    """Screen, polish, sweep the axial registry, then release the per-chain slides.

    The first four stages give the staggered cell exactly the five free variables the aligned
    cells were relaxed over -- ``(a, b, phi1, phi2, dz)`` at ``gamma = 90`` -- so the
    comparison is matched.  The last stage then lets each chain slide off its integer
    monomer by up to half a monomer, which is a degree of freedom the aligned cells could not
    have had, and it is reported separately: without it, "staggering did not help" could be a
    statement about the parameterisation rather than about the stagger.
    """
    obj = Objective(se, ch, pk, stagger, antipolar, bounds)
    per = obj.n_cell * ch.n_monomers
    rows, groups = screen(obj, ch, seed_params, rng)
    t0 = time.time()
    e = np.array([obj(x) for x in rows])
    log(f"    {label}: screened {len(rows)} cells in {time.time() - t0:.0f} s, "
        f"best {e.min() / per:+.4f} kcal/mol per monomer")
    probe = []
    for name, a, b, lo, hi in groups:
        rho = pk.density(a, b, 90.0) if a else float("nan")
        probe.append({"group": name, "a": a, "b": b, "density": rho,
                      "best_per_monomer": float(e[lo:hi].min() / per), "n": hi - lo})
        log(f"      {name:34s} rho {rho:6.4f}  best {e[lo:hi].min() / per:+9.4f}")
    obj.probe = probe
    best, best_e = None, np.inf
    for x0 in distinct(rows, e, 5, obj.i_dz):
        p, pe = polish(obj, x0)
        if pe < best_e:
            best, best_e = p, pe
    best, best_e = registry_sweep(obj, ch, best, best_e)
    best, best_e = polish(obj, best, maxfev=400)  # final tightening
    log(f"    {label}: {best_e / per:+.4f} kcal/mol per monomer (matched five variables) "
        f"after {obj.n_calls} energy evaluations")
    # Skipped for the aligned control on purpose: a free per-chain slide *is* a stagger, so
    # releasing it there would stop the control being a control.
    if not free_z or not np.any(np.asarray(stagger)):
        return obj, best, best_e, None
    fobj = Objective(se, ch, pk, stagger, antipolar, bounds, free_z=True)
    fx, fe = polish(fobj, fobj.x0(*_unpack(obj, best)), maxfev=1100)
    log(f"    {label}: {fe / per:+.4f} with the per-chain slides free "
        f"({(fe - best_e) / per:+.4f}); slides "
        f"{np.round(fx[fobj.n_base:], 3)} monomers")
    return obj, best, best_e, (fobj, fx, fe)


def _unpack(obj, x):
    p = obj.params(x)
    return (p[0], p[1], p[3], p[4], p[5])


# ------------------------------------------------------------------------- 4. reporting
def describe(se, se_ew, obj, ch, x, stagger, label):
    cell = obj.cell(x)
    p = obj.params(x)
    mu = cell_dipole(cell, ch.charges)
    pol = cell_polarization(cell, ch.charges)
    per = cell.n_chains * ch.n_monomers
    return {
        "label": label,
        "stagger_monomers": np.asarray(stagger, int).tolist(),
        "a_primitive": float(p[0]), "b_primitive": float(p[1]),
        "a_cell": float(NX * p[0]), "b_cell": float(NY * p[1]), "c_cell": float(ch.c),
        "gamma": 90.0, "phi1": float(p[3] % 360), "phi2": float(p[4] % 360),
        "dz": float(p[5] % ch.c), "flip": 0,
        "n_chains": cell.n_chains, "n_atoms": cell.n_atoms,
        "energy_per_cell_kcal_mol": float(se.energy(cell)),
        "energy_per_monomer_kcal_mol": float(se.energy(cell) / per),
        "energy_per_monomer_ewald_kcal_mol": float(se_ew.energy(cell) / per),
        "density_g_cm3": float(cell_density(cell, ch.mass)),
        "dipole_e_A": [float(v) for v in mu],
        "polarization_C_m2": [float(v) for v in pol],
        "polarization_magnitude_C_m2": float(np.linalg.norm(pol)),
        "image_sum_converged_to": float(se.energy_is_converged(cell)),
    }


def topology_block(cell, tag):
    t0 = time.time()
    rep = check_topology(cell, CP, scales=SCALES)
    log(f"    topology {tag} ({time.time() - t0:.0f} s)")
    log("      " + rep.describe().replace("\n", "\n      "))
    return rep, {
        "n_atoms": rep.n_atoms, "n_chains": rep.n_chains, "formula": rep.formula,
        "n_intended_bonds": rep.n_intended_bonds,
        "scales": list(SCALES),
        "max_bond_ratio": rep.max_bond_ratio,
        "min_nonbond_ratio": rep.min_nonbond_ratio,
        "safe_scale_window": list(rep.safe_scale_window),
        "min_interchain_distance_A": rep.min_interchain_distance,
        "min_interchain_ratio": rep.min_interchain_ratio,
        "min_interchain_pair": list(rep.min_interchain_pair),
        "min_intra_nonbond_distance_A": rep.min_intra_nonbond_distance,
        "min_interchain_by_element_A": {f"{a}-{b}": v for (a, b), v in rep.min_interchain_by_element.items()},
        "components": {str(k): v for k, v in rep.components.items()},
        "lost": {str(k): v for k, v in rep.missing.items()},
        "new": {str(k): v for k, v in rep.extra.items()},
        "new_interchain": {str(k): v for k, v in rep.extra_interchain.items()},
        "pass": rep.ok,
    }


def main():
    t_start = time.time()
    with open(os.path.join(OUT, "summary.json")) as f:
        aligned = json.load(f)

    log("== 1. the chain, and whether a monomer slide is a rigid slide ==")
    ch = periodic_chain(CP, [0] * CP.bonds_per_repeat, THREE_STATE)
    per_cell = NX * NY * 2 * ch.n_monomers
    sk = skeleton_periodicity(ch)
    log(f"  chain: {ch.n_atoms} atoms, c = {ch.c:.4f} A, {ch.n_monomers} monomers, "
        f"mass {ch.mass:.3f}")
    log(f"  a c/{ch.n_monomers} slide maps {sk['n_mapped_exactly']} of {sk['n_atoms']} atoms onto an "
        f"identical atom, worst displacement {sk['worst_displacement_A']:.2e} A")
    log(f"  the {len(sk['atoms_that_move'])} that move are the comonomer site: "
        + ", ".join(f"{i}{e}->{e2} {d} A" for i, e, e2, d in sk["atoms_that_move"]))

    log("== 2. the evaluator, checked against CrystalPacker ==")
    pk = CrystalPacker(ch, n_chains=2)
    pk_ew = CrystalPacker(ch, n_chains=2, coulomb="ewald", ewald=EwaldSpec())
    se, se_ew = SupercellEnergy(pk), SupercellEnergy(pk_ew)
    checks = {}
    for lab in ("polar", "antipolar"):
        d = aligned[lab]
        p = np.array([d["a"], d["b"], d["gamma"], d["phi1"], d["phi2"], d["dz"], d["flip"]], float)
        checks[lab] = {
            "dsf_2chain": se.check_against_packer(p, 1, 1),
            "dsf_8chain": se.check_against_packer(p, NX, NY),
            "ewald_8chain": se_ew.check_against_packer(p, NX, NY),
            "aligned_params": p.tolist(),
        }
        log(f"  {lab}: |supercell - packer| per monomer = "
            f"{checks[lab]['dsf_2chain']:.2e} (2-chain), {checks[lab]['dsf_8chain']:.2e} (8-chain), "
            f"{checks[lab]['ewald_8chain']:.2e} (8-chain, Ewald)")

    log("== 3. the stagger patterns, and how separated they keep the nitriles ==")
    a0 = aligned["polar"]["a"]
    b0 = aligned["polar"]["b"]
    pairs = close_column_pairs(a0, b0, 90.0)
    log(f"  close column pairs (axes within {CLOSE_COLUMN_A} A at the aligned polar cell): "
        + ", ".join(f"{t}-{u} {d:.2f} A" for t, u, d in pairs))
    cap = best_possible_min_offset(pairs, NX * NY * 2, ch.n_monomers)
    log(f"  no assignment of integer offsets can keep every close pair more than "
        f"{cap} monomers apart")
    quality = {}
    for name in PATTERNS:
        m = stagger_pattern(name, NX, NY, 2, ch.n_monomers)
        quality[name] = stagger_quality(m, pairs, ch.n_monomers)
        q = quality[name]
        log(f"  {name:10s} {m}  min offset {q['min_offset_monomers']} monomers "
            f"(mean {q['mean_offset_monomers']:.2f}), {q['n_distinct_heights']} distinct heights"
            + ("  <- optimal" if q["min_offset_monomers"] == cap else ""))

    log("== 3b. what the stagger costs at the aligned cell's own geometry ==")
    # The cleanest separation there is: apply each pattern to the shipped cell without moving
    # anything else.  A positive number means the aligned registry was already the better one
    # for the nitriles at that shape, so any gain from staggering has to come from the cell
    # being allowed to contract afterwards.
    fixed = {}
    for lab in ("polar", "antipolar"):
        p = np.array(checks[lab]["aligned_params"], float)
        base = se.energy(staggered_cell(ch, p, np.zeros(NX * NY * 2, int), nx=NX, ny=NY, packer=pk))
        fixed[lab] = {}
        for name in PATTERNS:
            m = stagger_pattern(name, NX, NY, 2, ch.n_monomers)
            e = se.energy(staggered_cell(ch, p, m, nx=NX, ny=NY, packer=pk))
            fixed[lab][name] = (e - base) / per_cell
            log(f"  {lab} {name:10s} {fixed[lab][name]:+.4f} kcal/mol per monomer "
                f"against the aligned registry at the same cell")
    log(f"  the shipped cells' own dz is already "
        + ", ".join(f"{lab} {aligned[lab]['dz'] / (ch.c / ch.n_monomers):.3f}"
                    for lab in ("polar", "antipolar"))
        + " monomers, so the two sublattices were never at the same height")

    log("== 4. relax each pattern, both polarities ==")
    bounds = (np.array([4.0, 4.0]), np.array([14.0, 14.0]))
    results: dict[str, dict] = {}
    for lab in ("polar", "antipolar"):
        d = aligned[lab]
        seed = (d["a"], d["b"], d["phi1"], d["phi2"], d["dz"])
        anti = lab == "antipolar"
        for name in PATTERNS:
            rng = np.random.default_rng(20260911)
            m = stagger_pattern(name, NX, NY, 2, ch.n_monomers)
            key = f"{lab}_{name}"
            t0 = time.time()
            obj, x, e, free = relax(se, ch, pk, m, anti, seed, bounds, rng, key)
            r = describe(se, se_ew, obj, ch, x, m, key)
            r["stagger_quality"] = quality[name]
            r["seconds"] = time.time() - t0
            r["n_energy_evaluations"] = obj.n_calls
            r["shape_probe"] = obj.probe
            r["cost_at_the_aligned_geometry"] = fixed[lab][name]
            if free is not None:
                fobj, fx, fe = free
                fr = describe(se, se_ew, fobj, ch, fx, m, key + " + free slides")
                fr["slides_monomers"] = [0.0] + [float(v) for v in fx[fobj.n_base:]]
                fr["d_energy_per_monomer_vs_matched"] = (
                    fr["energy_per_monomer_kcal_mol"] - r["energy_per_monomer_kcal_mol"])
                r["free_slide"] = fr
            results[key] = r
            log(f"    {key}: a x b x c = {r['a_cell']:.4f} x {r['b_cell']:.4f} x {r['c_cell']:.4f}, "
                f"rho {r['density_g_cm3']:.4f}, E/mon {r['energy_per_monomer_kcal_mol']:+.4f} "
                f"(Ewald {r['energy_per_monomer_ewald_kcal_mol']:+.4f}), "
                f"|P| {r['polarization_magnitude_C_m2']:.3e} C/m2")
            results[key]["_x"] = np.asarray(x, float).tolist()
            results[key]["_obj"] = obj

    log("== 5. the aligned control ==")
    for lab in ("polar", "antipolar"):
        got = results[f"{lab}_aligned"]
        want = aligned[lab]
        log(f"  {lab}: the eight-chain search with no stagger returns "
            f"E/mon {got['energy_per_monomer_kcal_mol']:+.4f} rho {got['density_g_cm3']:.4f} "
            f"against the shipped two-chain cell's {want['energy_per_monomer_kcal_mol']:+.4f} / "
            f"{want['density_g_cm3']:.4f}")
        got["matches_shipped_aligned_cell"] = {
            "d_energy_per_monomer": got["energy_per_monomer_kcal_mol"] - want["energy_per_monomer_kcal_mol"],
            "d_density": got["density_g_cm3"] - want["density_g_cm3"],
        }

    log("== 6. topology, and export ==")
    for lab in ("polar", "antipolar"):
        for name in PATTERNS:
            key = f"{lab}_{name}"
            r = results[key]
            obj = r.pop("_obj")
            cell = obj.cell(np.array(r["_x"]))
            rep, block = topology_block(cell, key)
            r["topology"] = block
            if name not in SHIP:
                continue
            fname = f"vdf11_vdcn1_alltrans_{lab}_stagger_{name}_8chain_592atom.xyz"
            extra = (f'polyfind_label="{lab} staggered-{name} all-trans 11VDF:1VDCN" '
                     f'stagger_monomers="{" ".join(str(v) for v in r["stagger_monomers"])}" '
                     f'energy_per_monomer_kcal_mol={r["energy_per_monomer_kcal_mol"]:.6f} '
                     f'density_g_cm3={r["density_g_cm3"]:.6f} '
                     f'min_interchain_distance_A={rep.min_interchain_distance:.4f}')
            with open(os.path.join(OUT, fname), "w") as f:
                f.write(cell.to_extxyz(extra))
            cname = f"vdf11_vdcn1_alltrans_{lab}_stagger_{name}_8chain.cif"
            with open(os.path.join(OUT, cname), "w") as f:
                f.write(cell_to_cif(cell, title=f"vdf11_vdcn1_{lab}_stagger_{name}"))
            r["files"] = {"xyz": fname, "cif": cname}
            log(f"    wrote {fname} and {cname}")

    log("== 7. did breaking the registry recover the density? ==")
    verdict = {}
    for lab in ("polar", "antipolar"):
        base = results[f"{lab}_aligned"]
        # A pattern the packer's own ``dz`` can already reach is not a stagger, and including
        # it would let the comparison answer itself: ``antiphase`` relaxes back to the aligned
        # cell exactly, so "the best staggered cell equals the aligned cell" would be a
        # statement about that degeneracy rather than about staggering.  Degeneracy is decided
        # by measurement -- same relaxed energy to 1e-6 -- not by which name we gave a pattern.
        e0 = base["energy_per_monomer_kcal_mol"]
        genuine = [n for n in PATTERNS if n != "aligned"
                   and abs(results[f"{lab}_{n}"]["energy_per_monomer_kcal_mol"] - e0) > 1e-6]
        degenerate = [n for n in PATTERNS if n != "aligned" and n not in genuine]
        if degenerate:
            log(f"  {lab}: {', '.join(degenerate)} relaxed back to the aligned cell exactly "
                f"(dz already reaches that registry), so excluded from the comparison")
        best = min((results[f"{lab}_{n}"] for n in genuine),
                   key=lambda r: r["energy_per_monomer_kcal_mol"])
        dense = max((results[f"{lab}_{n}"] for n in genuine),
                    key=lambda r: r["density_g_cm3"])
        verdict[lab] = {
            "genuine_stagger_patterns": genuine,
            "patterns_degenerate_with_aligned": degenerate,
            "aligned_density": base["density_g_cm3"],
            "aligned_energy_per_monomer": base["energy_per_monomer_kcal_mol"],
            "lowest_staggered": best["label"],
            "lowest_staggered_energy_per_monomer": best["energy_per_monomer_kcal_mol"],
            "lowest_staggered_density": best["density_g_cm3"],
            "densest_staggered": dense["label"],
            "densest_staggered_density": dense["density_g_cm3"],
            "d_energy_per_monomer": best["energy_per_monomer_kcal_mol"] - base["energy_per_monomer_kcal_mol"],
            "d_density": dense["density_g_cm3"] - base["density_g_cm3"],
            "mixing_rule_density": MIXING_RULE_DENSITY,
            "fraction_of_deficit_recovered": (
                (dense["density_g_cm3"] - base["density_g_cm3"])
                / (MIXING_RULE_DENSITY - base["density_g_cm3"])),
            "energy_lowered": best["energy_per_monomer_kcal_mol"] < base["energy_per_monomer_kcal_mol"],
            "density_recovered": dense["density_g_cm3"] > base["density_g_cm3"],
        }
        fs = [results[f"{lab}_{n}"]["free_slide"] for n in genuine
              if "free_slide" in results[f"{lab}_{n}"]]
        if fs:
            bf = min(fs, key=lambda r: r["energy_per_monomer_kcal_mol"])
            verdict[lab]["best_free_slide"] = bf["label"]
            verdict[lab]["best_free_slide_energy_per_monomer"] = bf["energy_per_monomer_kcal_mol"]
            verdict[lab]["best_free_slide_density"] = bf["density_g_cm3"]
            verdict[lab]["free_slide_beats_aligned"] = (
                bf["energy_per_monomer_kcal_mol"] < base["energy_per_monomer_kcal_mol"])
        v = verdict[lab]
        log(f"  {lab}: lowest staggered {v['lowest_staggered']} at "
            f"{v['lowest_staggered_energy_per_monomer']:+.4f} vs aligned "
            f"{v['aligned_energy_per_monomer']:+.4f} ({v['d_energy_per_monomer']:+.4f}); "
            f"densest {v['densest_staggered']} at {v['densest_staggered_density']:.4f} vs "
            f"{v['aligned_density']:.4f} ({v['d_density']:+.4f}, "
            f"{100 * v['fraction_of_deficit_recovered']:+.1f}% of the deficit to the mixing rule)")

    out = {
        "meta": {
            "polymer": CP.name, "formula": CP.formula,
            "monomers_per_repeat": CP.monomers_per_repeat,
            "mol_percent_vdcn": CP.mol_percent()["vdcn"],
            "chain_atoms": ch.n_atoms, "chain_c_A": ch.c, "chain_mass": ch.mass,
            "nx": NX, "ny": NY, "n_chains": NX * NY * 2, "monomers_per_cell": per_cell,
            "scales": list(SCALES),
            "mixing_rule_density": MIXING_RULE_DENSITY,
            "column_layout": column_layout(NX, NY, 2),
            "close_column_pairs": [[t, u, d] for t, u, d in pairs],
            "best_possible_min_offset_monomers": cap,
            "skeleton_periodicity": sk,
            "evaluator_checks": checks,
            "patterns": {k: quality[k] for k in PATTERNS},
            "cost_at_the_aligned_geometry": fixed,
            "aligned_dz_monomers": {lab: aligned[lab]["dz"] / (ch.c / ch.n_monomers)
                                    for lab in ("polar", "antipolar")},
            "shipped": list(SHIP),
            "seconds": time.time() - t_start,
        },
        "cells": {k: {kk: vv for kk, vv in v.items() if not kk.startswith("_")}
                  for k, v in results.items()},
        "verdict": verdict,
    }
    with open(os.path.join(OUT, "staggered_summary.json"), "w") as f:
        json.dump(out, f, indent=1)
    log(f"done in {time.time() - t_start:.0f} s")


if __name__ == "__main__":
    sys.exit(main())
