"""Build the topology-checked all-trans 11 VDF : 1 VDCN starting structures.

Produces, under ``deliverables/``, a polar and an independently constructed antipolar
packing of the 8.33 mol% copolymer, each as a two-chain primitive cell and as the
592-atom eight-chain cell asked for in ``docs/NOTE_VDCN_CONVERGENCE.md``, with the full
topology and close-contact report and a README naming the provenance.

Run with ``PYTHONPATH=src python examples/copolymer_starts.py``.
"""
from __future__ import annotations

import json
import os
import sys
import time
import warnings

import numpy as np

from polyfind.ewald import EwaldSpec
from polyfind.fitting import antipolar_cell_exact, antipolar_offsets, chain_moment
from polyfind.forcefield import SimpleFF, fit_ris
from polyfind.pack import CrystalPacker, default_bounds, pack, periodic_chain, polish, to_cif
from polyfind.polymers import PVDF, THREE_STATE, VDCN, VDF_VDCN_11_1
from polyfind.ris import transfer_ris, transfer_summary
from polyfind.topology import build_cell, check_topology, repeat_bond_graph

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, "deliverables")
CP = VDF_VDCN_11_1
KEYS = ["a", "b", "gamma", "phi1", "phi2", "dz"]


def log(*a):
    print(*a, flush=True)


# --------------------------------------------------------------------- 1. the chain
def build_reference_chain():
    ch = periodic_chain(CP, [0] * CP.bonds_per_repeat, THREE_STATE)
    starts, bonds, n_rep = repeat_bond_graph(CP)
    assert np.array_equal(np.asarray(ch.backbone), starts), (ch.backbone, starts)
    assert n_rep == ch.n_atoms == 74, (n_rep, ch.n_atoms)
    # every intended bond, including the one that closes the repeat onto its own image, is a
    # real bond length in the built coordinates -- checked here because the topology report
    # later leans on the same edge list
    bb = set(int(k) for k in starts)
    worst_bb = 0.0
    for i, j, dk in bonds:
        d = float(np.linalg.norm(ch.coords[i] - (ch.coords[j] + np.array([0.0, 0.0, dk * ch.c]))))
        if i in bb and j in bb:
            worst_bb = max(worst_bb, abs(d - CP.bond_length))
        else:
            assert 0.9 < d < 1.9, (i, j, dk, d)
    assert worst_bb < 1e-6, worst_bb
    return ch


# -------------------------------------------------------- 2. RIS parameter provenance
def ris_provenance():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        vdf = fit_ris(PVDF, SimpleFF(), step=10.0, n_monomers=5, third_order=True, adapt_angles=False)
        # angles="rigid": VDCN now defaults to the angle-relaxed scan (docs/NITRILE_LANDSCAPE.md;
        # its rigid T state is a +/-120 deg basin edge), which is the better model but a
        # different one from the one every result this script recorded was built on, and a
        # hundred times slower.  Pinned so the script goes on reproducing those results.
        cn = fit_ris(VDCN, SimpleFF(), step=10.0, n_monomers=5, third_order=True, adapt_angles=False, angles="rigid")
        model, report = transfer_ris(CP, [(PVDF, vdf.model), (VDCN, cn.model)], name=CP.name)
    return model, report


# ------------------------------------------------------------------ 3. polar packing
def pvdf_seed():
    """The beta-like all-trans two-chain cell of the PVDF homopolymer, as a seed.

    The copolymer's backbone is geometrically identical to PVDF's all-trans backbone (same
    1.528 A bond, same 114 deg angles, the only chemical difference being one carbon's
    pendants), so the homopolymer's packing is the right neighbourhood to start from, and
    it is found exhaustively and cheaply on a two-bond repeat.
    """
    ch = periodic_chain(PVDF, [0, 0], THREE_STATE)
    res = pack(ch, n_chains=2, screen="random", n_random=4000, n_refine=6)
    return res[0], ch


def polar_cell(pk, ch, seed, lo, hi, free, n_random=2000, seconds_budget=None):
    c = ch.c
    c_mon = c / CP.monomers_per_repeat
    starts = []
    # seeded: the homopolymer cell, with every distinct axial registry of the VDCN units.
    # The chain is all but exactly 12-fold periodic along z, so dz and dz + c/12 are two
    # different VDCN-VDCN registries and the same backbone packing: twelve real starts that
    # a local polish cannot move between.
    for flip in (0.0, 1.0):
        for j in range(CP.monomers_per_repeat):
            starts.append([seed.a, seed.b, 90.0, seed.phi1, seed.phi2, (seed.dz + j * c_mon) % c, flip])
    n_seed = len(starts)
    rng = np.random.default_rng(20260910)
    for flip in (0.0, 1.0):
        u = rng.random((n_random, 6))
        g = lo + u * (hi - lo)
        starts.extend(np.concatenate([g, np.full((n_random, 1), flip)], axis=1))
    grid = np.array(starts, dtype=float)
    t0 = time.time()
    e = pk.energy(grid)
    log(f"    polar screen: {len(grid)} cells ({n_seed} seeded) in {time.time() - t0:.0f} s; "
        f"best seeded {e[:n_seed].min() / (2 * ch.n_monomers):.4f}, best random "
        f"{e[n_seed:].min() / (2 * ch.n_monomers):.4f} kcal/mol per monomer")
    order = np.argsort(e)
    taken = []
    for i in order:
        g = grid[i]
        if any(abs(g[0] - t[0]) + abs(g[1] - t[1]) < 0.6 and abs(g[6] - t[6]) < 0.5 for t in taken):
            continue
        taken.append(g)
        if len(taken) >= 10:
            break
    best, best_e = None, np.inf
    for g in taken:
        p = polish(pk, g, lo, hi, free, maxfev=600, method="lbfgs", gradient="analytic")
        ee = float(pk.energy(p[None])[0])
        if ee < best_e:
            best, best_e = p, ee
    best, best_e = registry_sweep(pk, ch, best, best_e, lo, hi, free)
    return shape_probe(pk, ch, best, best_e, lo, hi, free)


def shape_probe(pk, ch, best, best_e, lo, hi, free, n_sample=60, n_polish=3):
    """Convergence check: is the cell found also reached from denser, differently shaped starts?

    A random screen over five continuous variables is sparse, and the quantity a consumer
    will look at first is the density.  So this fixes ``(a, b)`` at several shapes whose cell
    areas would give 1.75-1.95 g/cm3 -- including the consumer's own accepted density --
    samples the setting angles and the axial shift densely at each, and polishes the best few
    with every variable free.  If they all come back to the same cell, the cell is the
    minimum of this two-chain problem rather than wherever the screen happened to look; if
    the dense shapes are all repulsive, the loose cell is what the protruding nitrile costs
    and not a search failure.  Both outcomes are logged.
    """
    c = ch.c
    rng = np.random.default_rng(7)
    shapes = [(4.6, 9.4), (5.0, 9.0), (5.4, 8.4), (6.0, 8.0), (5.2, 9.6), (4.8, 10.0)]
    rows = []
    for a, b in shapes:
        line = []
        for flip in (0.0, 1.0):
            g = np.column_stack([np.full(n_sample, a), np.full(n_sample, b), np.full(n_sample, 90.0),
                                 rng.uniform(0, 360, n_sample), rng.uniform(0, 360, n_sample),
                                 rng.uniform(0, c, n_sample), np.full(n_sample, flip)])
            e = pk.energy(g) / (2 * ch.n_monomers)
            i = int(np.argmin(e))
            rows.append((float(e[i]), g[i]))
            line.append(f"flip{flip:.0f} {e[i]:+8.3f}")
        log(f"    shape probe a={a:5.2f} b={b:5.2f} (rho {pk.density(a, b, 90.0):.3f}): " + "  ".join(line))
    rows.sort(key=lambda r: r[0])
    out, out_e = best, best_e
    for e0, g in rows[:n_polish]:
        p = polish(pk, g, lo, hi, free, maxfev=600, method="lbfgs", gradient="analytic")
        ee = float(pk.energy(p[None])[0])
        log(f"    shape probe polish {e0:+8.3f} -> {ee / (2 * ch.n_monomers):+8.4f} per monomer  "
            f"a={p[0]:.4f} b={p[1]:.4f} rho={pk.density(p[0], p[1], p[2]):.4f}")
        if ee < out_e:
            out, out_e = p, ee
    log("    shape probe " + ("found nothing lower: the cell above is the minimum of this "
                              "two-chain problem" if out_e >= best_e - 1e-9 else
                              "IMPROVED on the screened cell"))
    return out, out_e


def registry_sweep(pk, ch, best, best_e, lo, hi, free, dphi=None, flip=None, n_keep=4):
    """Re-polish the twelve axial VDCN registries of a polished cell.

    ``dz`` is the one cell variable a local polish cannot cross.  The chain is all but
    exactly twelve-fold periodic along z, so ``dz`` and ``dz + c/12`` are two different
    VDCN-VDCN registries with the same backbone packing, separated by a barrier the
    optimiser will not climb.  Sweeping them at the cell the search actually found -- rather
    than at a guessed one -- is the cheap way to be sure the registry is the best one for
    that packing.  With ``dphi`` given the sweep stays inside an antipolar subspace by
    tying ``phi2`` to ``phi1``.
    """
    c = ch.c
    grid = np.repeat(np.asarray(best, dtype=float)[None], CP.monomers_per_repeat, axis=0)
    grid[:, 5] = (best[5] + np.arange(CP.monomers_per_repeat) * c / CP.monomers_per_repeat) % c
    e = pk.energy(grid)
    log(f"    registry sweep: twelve dz offsets span "
        f"{e.min() / (2 * ch.n_monomers):+.4f} to {e.max() / (2 * ch.n_monomers):+.4f} kcal/mol per monomer")
    out, out_e = best, best_e
    for i in np.argsort(e)[:n_keep]:
        if dphi is None:
            p = polish(pk, grid[i], lo, hi, free, maxfev=600, method="lbfgs", gradient="analytic")
            ee = float(pk.energy(p[None])[0])
        else:
            from scipy.optimize import minimize

            def obj(v):
                a, b, phi, dz = v
                if a < 2.0 or b < 2.0:
                    return 1e6
                return float(pk.energy(np.array([a, b, 90.0, phi, phi + dphi, dz, float(flip)])[None])[0])

            r = minimize(obj, [grid[i][0], grid[i][1], grid[i][3], grid[i][5]], method="Nelder-Mead",
                         options={"maxfev": 300, "xatol": 1e-3, "fatol": 1e-5})
            p = np.array([r.x[0], r.x[1], 90.0, r.x[2], r.x[2] + dphi, r.x[3], float(flip)])
            ee = float(r.fun)
        if ee < out_e:
            out, out_e = p, ee
    return out, out_e


# -------------------------------------------------------------- 4. antipolar packing
def antipolar_fine(pk, ch, a0, b0, lo, hi, n_phi=24, n_dz=60):
    """A finer screen of the *same* exactly-antipolar subspace antipolar_cell_exact uses.

    When this was written ``antipolar_cell_exact`` screened ``dz`` on four points across the
    whole repeat, which for a twelve-monomer repeat is a 7.7 A step -- coarser than the
    2.56 A monomer period the interchain registry actually varies on -- and it returned
    -3.7152 where this scan returned -3.7514.  The helper's axial resolution is now a length
    (``dz_step``, 0.5 A) rather than a point count, so that particular blindness is gone; its
    remaining cost on a 74-atom, 30.8 A repeat is the ``(a, b)`` axis, at 284 ms per
    exact-kernel cell.  This stays as the cheap independent check it always was: a 0.5 A
    ``dz`` scan and a 15 deg ``phi`` scan at the polar cell's own ``(a, b)``, inside the
    identical subspace, handed to the same constrained polish.  The two are compared and the
    lower is taken.
    """
    from scipy.optimize import minimize

    c = ch.c
    branches = antipolar_offsets(pk)
    rows = []
    for flip, dphi in branches:
        grid = np.array([[a0, b0, 90.0, p, p + dphi, z, float(flip)]
                         for p in np.linspace(0.0, 360.0, n_phi, endpoint=False)
                         for z in np.linspace(0.0, c, n_dz, endpoint=False)])
        e = pk.energy(grid)
        taken = []
        for i in np.argsort(e):
            g = grid[i]
            if any(abs((g[3] - t[3] + 180) % 360 - 180) < 20.0 and abs(g[5] - t[5]) < 1.0 for t in taken):
                continue
            taken.append(g)
            rows.append((float(e[i]), g, flip, dphi))
            if len(taken) >= 4:
                break
    rows.sort(key=lambda r: r[0])
    best, best_e = None, np.inf
    for _, g, flip, dphi in rows[:6]:
        def obj(v, _f=flip, _d=dphi):
            a, b, phi, dz = v
            if a < 2.0 or b < 2.0:
                return 1e6
            return float(pk.energy(np.array([a, b, 90.0, phi, phi + _d, dz, float(_f)])[None])[0])

        r = minimize(obj, [g[0], g[1], g[3], g[5]], method="Nelder-Mead",
                     options={"maxfev": 300, "xatol": 1e-3, "fatol": 1e-5})
        if r.fun < best_e:
            best = np.array([r.x[0], r.x[1], 90.0, r.x[2], r.x[2] + dphi, r.x[3], float(flip)])
            best_e = float(r.fun)
    return best, best_e


# ------------------------------------------------------------------------- 5. report
def describe_cell(pk, ch, params, label):
    res = pk.result(params)
    pol = pk.polarization(np.asarray(params)[None])[0]
    return {
        "label": label,
        "a": res.a, "b": res.b, "c": res.c, "gamma": res.gamma,
        "phi1": res.phi1, "phi2": res.phi2, "dz": res.dz, "flip": res.flip,
        "n_chains_in_primitive": res.n_chains,
        "energy_per_cell_kcal_mol": res.energy_per_cell,
        "energy_per_monomer_kcal_mol": res.energy_per_monomer,
        "density_g_cm3": res.density,
        "dipole_e_A": [float(v) for v in res.dipole],
        "polarization_C_m2": [float(v) for v in pol],
        "polarization_magnitude_C_m2": float(np.linalg.norm(pol)),
    }


def main():
    os.makedirs(OUT, exist_ok=True)
    t_start = time.time()
    log("== 1. chain ==")
    ch = build_reference_chain()
    log(f"  {CP.name}: B={CP.bonds_per_repeat} bonds, {CP.monomers_per_repeat} monomers per repeat, "
        f"{CP.composition()}, {CP.mol_percent()['vdcn']:.2f} mol% VDCN")
    log(f"  chain: {ch.n_atoms} atoms, c = {ch.c:.4f} A, mass {ch.mass:.3f}, radius {ch.radius:.3f} A, "
        f"helix '{ch.helix.label}', charge sum {float(ch.charges.sum()):+.2e} e")

    log("== 2. RIS parameter provenance ==")
    model, report = ris_provenance()
    log(transfer_summary(report))
    with open(os.path.join(OUT, "ris_parameter_provenance.json"), "w") as f:
        json.dump({"polymer": CP.name,
                   "terms": [vars(t) | {"exact": t.exact} for t in report]}, f, indent=1)

    pk = CrystalPacker(ch, n_chains=2)
    m = chain_moment(pk)
    log(f"  chain moment (e.A) = ({m[0]:+.4f}, {m[1]:+.4f}, {m[2]:+.4f}); "
        f"antipolar branches (flip, dphi) = {antipolar_offsets(pk)}")
    bd = default_bounds(ch)
    bd["a"], bd["b"] = (4.0, 14.0), (4.0, 14.0)  # widened: default_bounds scales off the
    # chain radius, which one protruding nitrile in twelve monomers inflates, and it would
    # otherwise exclude the beta-PVDF-like 4.6 A polar axis entirely
    lo = np.array([bd[k][0] for k in KEYS])
    hi = np.array([bd[k][1] for k in KEYS])
    free = [0, 1, 3, 4, 5]

    log("== 3. polar packing ==")
    seed, pvdf_ch = pvdf_seed()
    log(f"  PVDF all-trans seed: {seed.row()}")
    p_polar, e_polar = polar_cell(pk, ch, seed, lo, hi, free)
    log(f"  polar: {np.round(p_polar, 4)}  E/cell {e_polar:.4f}")

    log("== 4. antipolar packing (independently constructed) ==")
    t0 = time.time()
    p_a1, e_a1, pol_a1 = antipolar_cell_exact(pk, bounds=bd, target=1200, n_polish=4, maxfev=300)
    log(f"  antipolar_cell_exact: E/mon {e_a1:.4f}  max|P| {pol_a1:.2e}  ({time.time() - t0:.0f} s)")
    t0 = time.time()
    p_a2, e_a2 = antipolar_fine(pk, ch, p_polar[0], p_polar[1], lo, hi)
    log(f"  fine subspace scan:   E/mon {e_a2 / (2 * ch.n_monomers):.4f}  ({time.time() - t0:.0f} s)")
    p_anti = p_a1 if e_a1 * (2 * ch.n_monomers) <= e_a2 else p_a2
    which = "antipolar_cell_exact" if p_anti is p_a1 else "fine subspace scan"
    log(f"  taking the {which}: {np.round(p_anti, 4)}")
    # the same axial-registry sweep the polar branch gets, kept inside the subspace by tying
    # phi2 to phi1 with the branch's own dphi
    e_anti = float(pk.energy(np.asarray(p_anti)[None])[0])
    dphi = float(p_anti[4] - p_anti[3])
    p_anti, e_anti = registry_sweep(pk, ch, p_anti, e_anti, lo, hi, free,
                                    dphi=dphi, flip=float(p_anti[6]))
    log(f"  after the registry sweep: {np.round(p_anti, 4)}  "
        f"E/mon {e_anti / (2 * ch.n_monomers):+.4f}")

    log("== 5. cells, energies, polarizations ==")
    pk_ew = CrystalPacker(ch, n_chains=2, coulomb="ewald", ewald=EwaldSpec())
    summary = {}
    for label, params in (("polar", p_polar), ("antipolar", p_anti)):
        d = describe_cell(pk, ch, params, label)
        d["energy_per_monomer_ewald_kcal_mol"] = float(
            pk_ew.energy(np.asarray(params)[None])[0] / (2 * ch.n_monomers))
        summary[label] = d
        log(f"  {label}: a={d['a']:.4f} b={d['b']:.4f} c={d['c']:.4f} gamma={d['gamma']:.1f} "
            f"phi=({d['phi1']:.2f},{d['phi2']:.2f}) dz={d['dz']:.4f} flip={d['flip']}")
        log(f"    rho={d['density_g_cm3']:.4f} g/cm3  E/mon={d['energy_per_monomer_kcal_mol']:+.4f} "
            f"(Ewald {d['energy_per_monomer_ewald_kcal_mol']:+.4f})  "
            f"|P|={d['polarization_magnitude_C_m2']:.6f} C/m2  P={np.round(d['polarization_C_m2'], 6)}")

    log("== 6. topology and close contacts ==")
    for label, params in (("polar", p_polar), ("antipolar", p_anti)):
        for tag, (nx, ny) in (("primitive", (1, 1)), ("8chain", (2, 2))):
            cell = build_cell(ch, params, nx=nx, ny=ny, packer=pk)
            t0 = time.time()
            rep = check_topology(cell, CP)
            log(f"  {label} / {tag} ({time.time() - t0:.0f} s)")
            log("    " + rep.describe().replace("\n", "\n    "))
            summary[label].setdefault("topology", {})[tag] = {
                "n_atoms": rep.n_atoms, "n_chains": rep.n_chains, "formula": rep.formula,
                "n_intended_bonds": rep.n_intended_bonds,
                "max_bond_ratio": rep.max_bond_ratio,
                "min_nonbond_ratio": rep.min_nonbond_ratio,
                "safe_scale_window": list(rep.safe_scale_window),
                "min_interchain_distance_A": rep.min_interchain_distance,
                "min_interchain_ratio": rep.min_interchain_ratio,
                "min_intra_nonbond_distance_A": rep.min_intra_nonbond_distance,
                "min_interchain_by_element_A": {f"{a}-{b}": v for (a, b), v in rep.min_interchain_by_element.items()},
                "components": {str(k): v for k, v in rep.components.items()},
                "lost": {str(k): v for k, v in rep.missing.items()},
                "new": {str(k): v for k, v in rep.extra.items()},
                "new_interchain": {str(k): v for k, v in rep.extra_interchain.items()},
                "pass": rep.ok,
            }
            name = f"vdf11_vdcn1_alltrans_{label}_{'8chain_592atom' if tag == '8chain' else '2chain_148atom'}.xyz"
            extra = (f'polyfind_label="{label} all-trans 11VDF:1VDCN" '
                     f'energy_per_monomer_kcal_mol={summary[label]["energy_per_monomer_kcal_mol"]:.6f} '
                     f'density_g_cm3={summary[label]["density_g_cm3"]:.6f} '
                     f'min_interchain_distance_A={rep.min_interchain_distance:.4f}')
            with open(os.path.join(OUT, name), "w") as f:
                f.write(cell.to_extxyz(extra))
            summary[label].setdefault("files", {})[tag] = name
            log(f"    wrote {name}")
        with open(os.path.join(OUT, f"vdf11_vdcn1_alltrans_{label}_2chain.cif"), "w") as f:
            f.write(to_cif(pk.result(params), title=f"vdf11_vdcn1_{label}"))

    summary["meta"] = {
        "polymer": CP.name, "formula": CP.formula,
        "bonds_per_repeat": CP.bonds_per_repeat,
        "monomers_per_repeat": CP.monomers_per_repeat,
        "mol_percent_vdcn": CP.mol_percent()["vdcn"],
        "chain_atoms": ch.n_atoms, "chain_c_A": ch.c, "chain_mass": ch.mass,
        "bond_length_A": CP.bond_length,
        "pvdf_seed": {"a": seed.a, "b": seed.b, "c": seed.c, "phi1": seed.phi1, "phi2": seed.phi2,
                      "dz": seed.dz, "flip": seed.flip, "energy_per_monomer": seed.energy_per_monomer},
        "seconds": time.time() - t_start,
    }
    with open(os.path.join(OUT, "summary.json"), "w") as f:
        json.dump(summary, f, indent=1)
    log(f"done in {time.time() - t_start:.0f} s")


if __name__ == "__main__":
    sys.exit(main())
