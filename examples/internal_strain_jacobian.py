"""Internal strain of beta-PVDF: our relaxed-ion Jacobian by the provider's construction, against theirs.

    python examples/internal_strain_jacobian.py                    # the comparison, tables to stdout
    python examples/internal_strain_jacobian.py --json out.json    # and every number to a file

The hypothesis under test (``docs/ELECTROMECHANICS.md`` 5.9, ``docs/BENCHMARK.md`` 2026-09-12): with
the Born charges now matching periodic DFPT, the piezoelectric shortfall must sit in *which atoms
move* under a macroscopic strain -- internal strain -- and a chain whose pendant groups ride
rigidly on their carbons supplies none of it.  The provider's geometry-only Jacobian (sarco
``ccdf377``, ``internal_strain_geometry.json``: thirteen CP2K fixed-cell relaxed-ion geometries
at +/-1% and +/-2% normal strain) is the first data that can test that, and this script does
what its report asks: build *our* Jacobian by the same construction, in their frame and atom
order, then compare the total pendant-relative bond-vector derivatives, affine term included --
the comparison that can tell a rigid pendant from a flexible one, because a rigid pendant
cancels the affine deformation of its bond and the nonaffine entries alone cannot show that.

**The construction, theirs and ours alike.**  ``u = r_relaxed - r_zero F`` with ``F`` the
homogeneous deformation, evaluated as the wrapped difference of fractional coordinates carried
into the strained cell; the equal-atom mean is projected out; ``J = (u(+h) - u(-h)) / 2h``.  For
each H or F, the bond vector from its nearest carbon at every state, its central derivative, the
rate along the reference bond (length) and the rate across it over the length (direction, degrees
per unit strain).  Ours is the deformable path of :mod:`polyfind.mechanics`
(:func:`~polyfind.mechanics.deformable_state`: setting angles, chain offset and the line-group
shape relaxed at fixed strain against the fitted valence terms) at the
``pvdf-dft-valence-flux-born`` structure with induced dipoles on -- the reference state of the
``d_33`` / ``d_31`` tables in ``docs/ELECTROMECHANICS.md`` 5.9.

**Frames.**  The packer's axes are ``x`` = polar, ``y`` = long lateral, ``z`` = chain; the
provider's are ``x`` = long (A, 8.358 A), ``y`` = polar (B, 4.731 A), ``z`` = chain (C, 2.580 A),
the ``(a, b, c)`` that :mod:`polyfind.born` already uses for the Born comparison.  The map is the
rotation by 90 degrees about the chain axis, ``(x, y, z)_provider = (-y, x, z)_packer`` up to the
polar sign, which is fixed as :func:`polyfind.born.canonical` fixes it: the fluorines sit at
``+y`` of their carbon, as in the provider's cell.  Atoms are then matched by fractional
coordinates.  Our Voigt ``eps_xx`` (polar) is their ``yy``, our ``eps_yy`` (long) their ``xx``,
``zz`` is ``zz``.

**What is and is not concluded.**  The provider's record is provisional -- force and position
tolerance, basis, k-grid and geometry sensitivity and the anomalous transverse stiffness are
unisolated -- so nothing is fitted to it, and it is not contracted with the DFPT Born tensor as
though a validated same-Hamiltonian result existed.  The one contraction made is labelled:
*their* kinematics with *our* Born-consistent charges, a cross-Hamiltonian diagnostic of whether
our charges would carry the missing response if our atoms moved as the reference's do.  Every
conclusion is a difference between the two models' pendant kinematics, not a confirmation of
real pendant flexibility, until the provider's gates clear.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np

from polyfind import backend as bk
from polyfind import mechanics as M
from polyfind.born import _carbon_of, born_charges, chain_frame_coords, load_internal_strain
from polyfind.fitting import FITTED_VALENCE
from polyfind.forcefield import SimpleFF
from polyfind.pack import E_PER_A2_TO_C_PER_M2
from polyfind.polarizability import Polarizable

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fit_born_flux import BETA, DEFAULT_DATA, beta_reference  # noqa: E402  (the 5.9 reference state)

AXES = ("x", "y", "z")  # the provider's: long, polar, chain
PACKER_VOIGT = {"x": 1, "y": 0, "z": 2}  # provider normal strain -> our Voigt index of the same strain
AMPLITUDES = (0.01, 0.02)
CONTINUITY = 0.25  # the provider's cap on a wrapped fractional displacement
DEFAULT_RECORD = os.path.join(os.path.dirname(DEFAULT_DATA), "internal_strain_geometry.json")
TARGET = (0.4, 0.6)  # C/m^2: the Berry slope less the dimensional |P|, docs/ELECTROMECHANICS.md 5.9


# --- geometry, the provider's way -----------------------------------------------------------
def frac(r: np.ndarray, cell: np.ndarray) -> np.ndarray:
    """Fractional coordinates for lattice vectors in the rows of ``cell``."""
    return r @ np.linalg.inv(cell)


def internal_displacement(zero: tuple, strained: tuple) -> tuple:
    """``(u - mean, mean, continuous fractional coordinates)``: ``u = r - r_zero F`` image-matched."""
    r0, c0 = zero
    r1, c1 = strained
    f0 = frac(r0, c0)
    d = frac(r1, c1) - f0
    d -= np.round(d)
    if np.abs(d).max() > CONTINUITY:
        raise RuntimeError(f"a fractional displacement of {np.abs(d).max():.3f} is not a continuous branch")
    u = d @ c1
    t = u.mean(axis=0)
    return u - t, t, f0 + d


def parents_of(elements, r0: np.ndarray, c0: np.ndarray) -> list:
    """``[(pendant, carbon)]``: every H or F with its nearest carbon under the minimum image."""
    f0 = frac(r0, c0)
    carbons = [i for i, e in enumerate(elements) if e == "C"]
    out = []
    for i, e in enumerate(elements):
        if e == "C":
            continue
        best, dbest = -1, np.inf
        for j in carbons:
            d = f0[i] - f0[j]
            d -= np.round(d)
            dist = float(np.linalg.norm(d @ c0))
            if dist < dbest:
                best, dbest = j, dist
        if dbest > 1.8:
            raise RuntimeError(f"atom {i} ({e}) has no carbon within 1.8 A")
        out.append((i, best))
    return out


def pair_response(zero: tuple, minus: tuple, plus: tuple, h: float, parents: list) -> dict:
    """One axis at one amplitude, exactly as the provider's ``pair_response``."""
    um, tm, fm = internal_displacement(zero, minus)
    up, tp, fp = internal_displacement(zero, plus)
    J = (up - um) / (2.0 * h)
    r0, c0 = zero
    f0 = frac(r0, c0)
    pend = []
    for atom, parent in parents:
        image = np.round(f0[atom] - f0[parent])
        v0 = (f0[atom] - f0[parent] - image) @ c0
        vm = (fm[atom] - fm[parent] - image) @ minus[1]
        vp = (fp[atom] - fp[parent] - image) @ plus[1]
        length = float(np.linalg.norm(v0))
        direction = v0 / length
        dv = (vp - vm) / (2.0 * h)
        dl = float(direction @ dv)
        rot = float(np.linalg.norm(dv - dl * direction)) / length
        pend.append({"atom": atom, "parent": parent, "dv": dv, "dlength": dl, "rotation_deg": float(np.degrees(rot)),
                     "length": length, "v0": v0})
    return {"J": J, "rms": float(np.sqrt(np.mean(np.sum(J ** 2, axis=1)))),
            "translation": {"negative": tm, "positive": tp}, "pendants": pend,
            "max_midpoint": float(np.linalg.norm(0.5 * (up + um), axis=1).max())}


def rigid_prediction(elements, r0: np.ndarray, axis: int, n_chain: int) -> np.ndarray:
    """``J`` of a cell whose chains are rigid bodies pinned to their own lattice sites.

    Nothing inside a chain moves, so every atom's nonaffine displacement is minus the affine
    one relative to its chain's mean: ``-(r_k - r_chain)`` along the strain axis, zero across
    it.  The remainder ``J - J_rigid`` is the motion *within* a chain, which is what the pendant
    hypothesis is about.
    """
    J = np.zeros((len(elements), 3))
    for s in range(len(elements) // n_chain):
        sl = slice(s * n_chain, (s + 1) * n_chain)
        J[sl, axis] = -(r0[sl, axis] - r0[sl, axis].mean())
    return J


# --- our cell in the provider's frame ---------------------------------------------------------
def placed(ref, params, chain) -> tuple:
    """``(positions (N, 3), lattice rows (3, 3))`` of ``chain`` at cell ``params``, packer frame."""
    P, lat = ref.packer._place(np.asarray(params, dtype=float)[None], coords=np.asarray(chain.coords, dtype=float)[None],
                               c=np.array([float(chain.c)]))
    return np.asarray(bk.to_numpy(P), dtype=float)[0], np.asarray(bk.to_numpy(lat), dtype=float)[0]


def frame_map(elements, r0: np.ndarray, n_chain: int, cz: float) -> np.ndarray:
    """``L`` (3, 3), ``provider = L @ packer``: the 90-degree turn about the chain axis that takes the
    packer's (polar, long, chain) to the provider's (long, polar, chain), with the polar sign such
    that the heavy pendants sit at ``+y`` of their carbon, as :func:`polyfind.born.canonical` has it."""
    carbon = {i: _carbon_of(elements, r0, i, n_chain, cz) for i in range(len(elements)) if elements[i] != "C"}
    heavy = [i for i in carbon if elements[i] != "H"]
    s = float(np.sign(np.mean([r0[i, 0] - r0[carbon[i], 0] for i in heavy])) or 1.0) if heavy else 1.0
    return np.array([[0.0, -s, 0.0], [s, 0.0, 0.0], [0.0, 0.0, 1.0]])


def to_provider(r: np.ndarray, cell: np.ndarray, L: np.ndarray) -> tuple:
    """Positions and lattice rows turned by ``L``, the rows re-ordered and signed so that row ``i``
    lies along ``+axis i`` (a lattice vector and its negative span the same lattice)."""
    rp = r @ L.T
    rows = cell @ L.T
    out = np.zeros((3, 3))
    for i in range(3):
        k = int(np.argmax(np.abs(rows[:, i])))
        out[i] = rows[k] * np.sign(rows[k, i])
    return rp, out


def match_atoms(el_ours, f_ours: np.ndarray, el_theirs, f_theirs: np.ndarray, cell_theirs: np.ndarray) -> tuple:
    """``(perm, max mismatch A)``: ``perm[i]`` is our atom for the provider's atom ``i``, by fractional
    coordinates after the translation that puts our first atom on theirs."""
    t = f_theirs[0] - f_ours[0]
    perm, used, worst = [], set(), 0.0
    for i in range(len(el_theirs)):
        best, dbest = -1, np.inf
        for j in range(len(el_ours)):
            if j in used or el_ours[j] != el_theirs[i]:
                continue
            d = f_ours[j] + t - f_theirs[i]
            d -= np.round(d)
            dist = float(np.linalg.norm(d @ cell_theirs))
            if dist < dbest:
                best, dbest = j, dist
        if best < 0:
            raise RuntimeError(f"no unmatched {el_theirs[i]} for the provider's atom {i + 1}")
        perm.append(best)
        used.add(best)
        worst = max(worst, dbest)
    return perm, worst


def our_states(ref, shape, amplitudes=AMPLITUDES) -> dict:
    """``{None: zero, (axis, +/-h): state}`` with each state ``(positions, rows, StrainState | None)``
    in the packer frame, the chain relaxed at fixed strain on the deformable path."""
    out = {None: placed(ref, ref.params, ref.packer.chain) + (None,)}
    for ax in AXES:
        for h in amplitudes:
            for sgn in (-1.0, 1.0):
                eps = np.zeros(6)
                eps[PACKER_VOIGT[ax]] = sgn * h
                st = M.deformable_state(ref, shape, eps)
                if st.relaxed.c_error > 1e-8 or abs(st.cell.theta) > 1e-9:
                    raise RuntimeError(f"strain {ax} {sgn * h:+.3f}: c error {st.relaxed.c_error:.1e}, theta {st.cell.theta:.1e}")
                out[(ax, sgn * h)] = placed(ref, st.params, st.relaxed.chain) + (st,)
    return out


# --- the contraction ------------------------------------------------------------------------
def contraction(Z: np.ndarray, J: np.ndarray, V0: float) -> np.ndarray:
    """``(1/V0) sum_k Z_k . J_k`` in C/m^2: the dipole per reference volume that displacements ``J``
    (A per unit strain) produce through Born tensors ``Z`` (N, 3, 3), both in the same frame and
    atom order.  A uniform translation drops out through the acoustic sum rule."""
    return np.einsum("kib,kb->i", Z, J) / V0 * E_PER_A2_TO_C_PER_M2


def by_type(Z: np.ndarray, J: np.ndarray, labels, V0: float) -> dict:
    """The polar (``y``) component of :func:`contraction`, split by atom type."""
    out = {}
    for lab in dict.fromkeys(labels):
        idx = [i for i, l in enumerate(labels) if l == lab]
        out[lab] = float(contraction(Z[idx], J[idx], V0)[1])
    return out


def clamped_ion(ref, h: float) -> dict:
    """Our **clamped-ion** proper coefficient, packer frame: ``{Voigt K: (3,) C/m^2}``.

    Every atom is displaced affinely with the cell and nothing is relaxed; the fluxed charges
    and the induced dipoles are re-solved at the deformed geometry and the dipole per reference
    volume is differenced.  With ``J`` the relaxed-ion internal Jacobian the proper coefficient
    splits as ``e = e_clamped + (1/V0) sum_k Z_k . J_k`` to second order in ``h``, which is the
    check that :func:`contraction` is the internal-strain term and nothing else.  (A rigid chain
    cannot be deformed affinely, so this state exists only as a dipole evaluation, never as a
    relaxation.)
    """
    packer = ref.packer
    P0, lat0 = placed(ref, ref.params, packer.chain)
    flux = getattr(packer, "_flux", None)
    n, nc = packer.n, packer.n_chains
    out = {}
    for K in (0, 1, 2):
        mus = []
        for sgn in (1.0, -1.0):
            eps = np.zeros(6)
            eps[K] = sgn * h
            sc = M.strained_cell(ref.params, ref.c, eps)
            F = np.eye(3) + M.strain_tensor(eps)
            Pn, latn = P0 @ F, lat0 @ F
            q = np.array(packer._q_cell, dtype=float)
            if flux is not None:
                for s in range(nc):
                    sl = slice(s * n, (s + 1) * n)
                    q[sl] = flux.charges(chain_frame_coords(sc.params, latn, Pn[sl], s), float(latn[2, 2]))[0]
            mu = q @ Pn
            if packer.polarizable is not None:
                mu = mu + packer._polarize(Pn, q, latn, float(sc.params[6]))[1].sum(axis=0)
            mus.append(mu)
        out[K] = (mus[0] - mus[1]) / (2.0 * h) / ref.volume * E_PER_A2_TO_C_PER_M2
    return out


def their_born(path: str, record) -> dict | None:
    """Per-atom DFPT Born tensors in the record's atom order, or ``None`` if the two do not share
    the zero geometry.  Provisional on the provider's own terms; used only in a labelled diagnostic."""
    if not os.path.exists(path):
        return None
    with open(path) as fh:
        d = json.load(fh)
    zero = record.receipts.get("zero", {})
    if d.get("geometry") != zero.get("geometry") or d.get("geometry_sha256") != zero.get("geometry_sha256"):
        return None
    atoms = d["atoms"]
    if [a["element"] for a in atoms] != record.elements:
        return None
    return {"corrected": np.array([a["born_effective_charge_asr_corrected_e"] for a in atoms], dtype=float),
            "raw": np.array([a["born_effective_charge_e"] for a in atoms], dtype=float),
            "protocol": d.get("protocol", "")}


# --- report ---------------------------------------------------------------------------------
def fmt3(v, w: int = 7, p: int = 3) -> str:
    return "(" + ", ".join(f"{float(x):+{w}.{p}f}" for x in v) + ")"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--data", default=DEFAULT_RECORD, help="the provider's internal_strain_geometry.json")
    ap.add_argument("--geometry", default=None, help="its zero-strain extended xyz, if not beside --data")
    ap.add_argument("--born", default=DEFAULT_DATA, help="the provider's born_results.json (labelled diagnostic only)")
    ap.add_argument("--preset", default="pvdf-dft-valence-flux-born")
    ap.add_argument("--json", default=None)
    args = ap.parse_args(argv)
    if not os.path.exists(args.data):
        print(f"{args.data} not found: it is sarco ccdf377's periodic_reference/beta_pvdf/internal_strain_geometry.json; "
              "pass --data (and --geometry for the zero xyz beside born_results.json)", file=sys.stderr)
        return 2
    rec = load_internal_strain(args.data, geometry=args.geometry)
    if rec.positions is None:  # the receipt's zero geometry beside born_results.json, then
        name = os.path.basename(rec.receipts.get("zero", {}).get("geometry", ""))
        cand = os.path.join(os.path.dirname(args.born), name)
        if name and os.path.exists(cand):
            rec = load_internal_strain(args.data, geometry=cand)
    born_path = args.born
    if not os.path.exists(born_path) and args.geometry:
        cand = os.path.join(os.path.dirname(args.geometry), "born_results.json")
        born_path = cand if os.path.exists(cand) else born_path
    if rec.positions is None or rec.cell is None:
        print("the record's zero geometry was not found; --geometry must name it", file=sys.stderr)
        return 2
    N = rec.n_atoms
    V0_theirs = float(abs(np.linalg.det(rec.cell)))
    print(f"# provider record: {os.path.basename(args.data)}, quantitatively_valid={rec.valid}, axes '{rec.convention}'")
    print(f"#   zero cell diag {np.diag(rec.cell)} A, V0 {V0_theirs:.3f} A^3, zero max force {rec.receipts['zero']['max_force_norm_eV_A']:.2e} eV/A")
    print("#   open gates: " + "; ".join(rec.gates))
    print(f"#   their diagnostics at 1%: atom rms x/y/z = " + "/".join(f"{rec.rms[(ax, 0.01)]:.6f}" for ax in AXES)
          + " A per unit strain; 1%->2% change " + ", ".join(f"{100 * rec.amplitude_change[ax]:.3f}%" for ax in AXES))
    out = {"record": {"path": args.data, "valid": rec.valid, "gates": rec.gates, "V0": V0_theirs,
                      "cell_diag": np.diag(rec.cell).tolist(), "receipts": rec.receipts}}

    # --- our reference state, the 5.9 one --------------------------------------------------
    t0 = time.time()
    packer, params, shape, rel = beta_reference(SimpleFF.from_preset(args.preset), Polarizable())
    ref = M.Reference(packer=packer, params=params, c=float(packer.chain.c), label=BETA[2])
    with FITTED_VALENCE.applied():
        ref, shape = M.relax_reference_deformable(ref, shape)
        states = our_states(ref, shape)
        born = born_charges(ref.packer, ref.params)
    n_chain = ref.packer.n
    el_ours = list(ref.packer.elements) * ref.packer.n_chains
    r0, c0, _ = states[None]
    L = frame_map(el_ours, r0, n_chain, float(ref.c))
    P_packer = ref.polarization()
    P_prov = L @ P_packer
    print(f"\n# our reference: {args.preset} + induced dipoles, deformable path, {time.time() - t0:.0f} s; "
          f"a={ref.params[0]:.4f} (polar) b={ref.params[1]:.4f} (long) gamma={ref.params[2]:.2f} c={ref.c:.4f}; "
          f"phi1={ref.params[3]:.1f} phi2={ref.params[4]:.1f} dz={ref.params[5]:.1e}; shape {shape.labels()} at {rel.x}")
    print(f"#   frame: provider = L . packer, L rows {L.tolist()} (det {np.linalg.det(L):+.0f}); "
          f"P packer {fmt3(P_packer, 7, 4)} -> provider {fmt3(P_prov, 7, 4)} C/m^2")
    V0_ours = ref.volume

    # the provider frame, and the atom matching
    rp0, cp0 = to_provider(r0, c0, L)
    perm, mismatch = match_atoms(el_ours, frac(rp0, cp0), rec.elements, frac(rec.positions, rec.cell), rec.cell)
    print(f"#   our lattice rows in their frame: diag {np.diag(cp0)} A (theirs {np.diag(rec.cell)}), V0 {V0_ours:.3f} A^3 "
          f"(theirs {V0_theirs:.3f}); atoms matched by fractional coordinates, worst {mismatch:.3f} A")
    print("#   atom map (their index: element, our packer index, their position / ours in their frame at their fractional):")
    f_ours = frac(rp0, cp0)
    f_theirs = frac(rec.positions, rec.cell)
    tshift = f_theirs[0] - f_ours[0]
    for i in range(N):
        j = perm[i]
        fo = f_ours[j] + tshift
        fo -= np.round(fo - f_theirs[i])
        print(f"      {i + 1:2d}: {rec.elements[i]:2s} <- {j:2d}  theirs {fmt3(rec.positions[i])}  ours {fmt3(fo @ rec.cell)}")
    out["frame"] = {"L": L.tolist(), "perm": perm, "mismatch_A": mismatch, "V0_ours": V0_ours,
                    "cell_diag_ours": np.diag(cp0).tolist(), "P_provider_frame": P_prov.tolist()}

    # --- our Jacobian, their construction, their frame and order ------------------------------
    prov = {k: to_provider(v[0], v[1], L) for k, v in states.items()}
    parents_ours = parents_of(el_ours, *prov[None])
    ours = {}
    for ax in AXES:
        for h in AMPLITUDES:
            r = pair_response(prov[None], prov[(ax, -h)], prov[(ax, h)], h, parents_ours)
            # into their atom order
            inv = {j: i for i, j in enumerate(perm)}
            r["J"] = r["J"][perm]
            r["rms"] = float(np.sqrt(np.mean(np.sum(r["J"] ** 2, axis=1))))
            for p in r["pendants"]:
                p["atom"], p["parent"] = inv[p["atom"]], inv[p["parent"]]
            r["pendants"].sort(key=lambda p: p["atom"])
            # our proper piezoelectric row from the same two states, for context
            sp, sm = states[(ax, h)][2], states[(ax, -h)][2]
            r["e_total"] = L @ ((sp.m - sm.m) / (2.0 * h))
            ours[(ax, h)] = r
    rp_theirs = rec.positions

    # --- 1. the Jacobians, atom by atom -------------------------------------------------------
    print("\n# === internal-strain Jacobian du/d eps, A per unit strain, provider frame (x long, y polar, z chain), their atom order ===")
    print("# 'rigid' is what a cell of rigid chains pinned to their lattice sites gives under a transverse strain")
    print("# (-(r - r_chain) along the strain axis); J - rigid is the motion within a chain.  No such reference exists for")
    print("# zz, which a rigid chain cannot take up at all.  Ours at 1% central difference; theirs at 1%.")
    out["jacobian"] = {}
    rms = lambda X: float(np.sqrt(np.mean(np.sum(X ** 2, axis=1))))  # noqa: E731
    for ax in AXES:
        a = AXES.index(ax)
        Jt, Jo = rec.J[(ax, 0.01)], ours[(ax, 0.01)]["J"]
        if ax == "z":
            Rt, Ro = np.zeros_like(Jt), np.zeros_like(Jo)
            Ft, Fo = Jt, Jo
            print(f"\n# strain zz:  rms theirs {rms(Jt):.4f}  ours {rms(Jo):.4f}  |theirs - ours| rms {rms(Jt - Jo):.4f}"
                  "  (all of J is within-chain motion)")
            print(f"  {'atom':6s} {'theirs J':>27s} {'ours J':>27s}")
            for i in range(N):
                print(f"  {i + 1:2d} {rec.elements[i]:2s}  {fmt3(Jt[i])}  {fmt3(Jo[i])}")
        else:
            Rt = rigid_prediction(rec.elements, rp_theirs, a, n_chain)
            Ro = rigid_prediction(rec.elements, rp0[perm], a, n_chain)
            Ft, Fo = Jt - Rt, Jo - Ro
            print(f"\n# strain {ax}{ax}:  rms theirs {rms(Jt):.4f} (rigid-chain value {rms(Rt):.4f}, within-chain remainder {rms(Ft):.4f})"
                  f"  ours {rms(Jo):.4f} (rigid {rms(Ro):.4f}, remainder {rms(Fo):.4f})  |theirs - ours| rms {rms(Jt - Jo):.4f}")
            print(f"  {'atom':6s} {'theirs J':>27s} {'ours J':>27s} {'theirs J - rigid':>27s} {'ours J - rigid':>27s}")
            for i in range(N):
                print(f"  {i + 1:2d} {rec.elements[i]:2s}  {fmt3(Jt[i])}  {fmt3(Jo[i])}  {fmt3(Ft[i])}  {fmt3(Fo[i])}")
        # where their atoms go: group centroids of the within-chain remainder
        groups = []
        for s in range(N // n_chain):
            base = s * n_chain
            for name, idx in (("CH2", [base, base + 1, base + 2]), ("CF2", [base + 3, base + 4, base + 5])):
                groups.append((s + 1, name, Ft[idx].mean(axis=0), Fo[idx].mean(axis=0)))
        print("  within-chain remainder, group centroids (theirs | ours): " + "; ".join(
            f"chain {s} {nm} {fmt3(gt)} | {fmt3(go)}" for s, nm, gt, go in groups))
        out["jacobian"][ax] = {"theirs": Jt.tolist(), "ours": Jo.tolist(), "rigid_theirs": Rt.tolist(), "rigid_ours": Ro.tolist(),
                               "rms": {"theirs": rms(Jt), "ours": rms(Jo), "rigid_theirs": rms(Rt), "rigid_ours": rms(Ro),
                                       "remainder_theirs": rms(Ft), "remainder_ours": rms(Fo), "difference": rms(Jt - Jo)},
                               "ours_2pc": ours[(ax, 0.02)]["J"].tolist(), "theirs_2pc": rec.J[(ax, 0.02)].tolist(),
                               "group_remainder": [{"chain": s, "group": nm, "theirs": gt.tolist(), "ours": go.tolist()}
                                                   for s, nm, gt, go in groups]}

    # --- 2. the bond vectors, affine term included --------------------------------------------
    print("\n# === total pendant-relative bond-vector derivative d(r_pendant - r_carbon)/d eps, A per unit strain, affine included ===")
    print("# length: rate along the reference bond (A per unit strain); direction: rotation rate (deg per unit strain).")
    print("# A rigid pendant on a rigid chain gives exactly zero here whatever the strain; affine motion alone would give")
    print("# v0 along the strain axis (the 'affine' column) for its length.")
    out["bonds"] = {}
    for ax in AXES:
        a = AXES.index(ax)
        pt = {p["atom"]: p for p in rec.pendants[(ax, 0.01)]}
        po = {p["atom"]: p for p in ours[(ax, 0.01)]["pendants"]}
        print(f"\n# strain {ax}{ax}:")
        print(f"  {'bond':10s} {'v0 (theirs)':>27s} {'affine':>8s} {'theirs dv':>27s} {'length':>8s} {'dir deg':>8s} {'ours dv':>27s} {'length':>8s} {'dir deg':>8s}")
        rows = []
        for i in sorted(pt):
            p, q = pt[i], po[i]
            # their reference bond vector from the zero geometry (same construction as theirs)
            f0 = frac(rp_theirs, rec.cell)
            image = np.round(f0[i] - f0[p["parent"]])
            v0t = (f0[i] - f0[p["parent"]] - image) @ rec.cell
            affine = float(v0t[a])
            lab = f"{rec.elements[p['parent']]}{p['parent'] + 1}-{p['element']}{i + 1}"
            print(f"  {lab:10s} {fmt3(v0t)} {affine:+8.3f} {fmt3(p['dv'])} {p['dlength']:+8.4f} {p['rotation_deg']:+8.3f} "
                  f"{fmt3(q['dv'])} {q['dlength']:+8.1e} {q['rotation_deg']:+8.1e}")
            rows.append({"bond": lab, "v0_theirs": v0t.tolist(), "affine": affine, "theirs": {"dv": p["dv"].tolist(),
                         "dlength": p["dlength"], "rotation_deg": p["rotation_deg"]},
                         "ours": {"dv": q["dv"].tolist(), "dlength": q["dlength"], "rotation_deg": q["rotation_deg"],
                                  "v0": q["v0"].tolist(), "length": q["length"]}})
        for el in ("F", "H"):
            sel = [i for i in pt if pt[i]["element"] == el]
            print(f"  C-{el}: theirs mean |length rate| {np.mean([abs(pt[i]['dlength']) for i in sel]):.4f} A/strain "
                  f"({100 * np.mean([abs(pt[i]['dlength']) for i in sel]) / np.mean([np.linalg.norm(rows[k]['v0_theirs']) for k, i2 in enumerate(sorted(pt)) if i2 in sel]):.1f}% of the bond per unit strain), "
                  f"mean rotation {np.mean([pt[i]['rotation_deg'] for i in sel]):.3f} deg/strain; "
                  f"ours {np.max([abs(po[i]['dlength']) for i in sel]):.1e} A/strain, {np.max([po[i]['rotation_deg'] for i in sel]):.1e} deg/strain")
        out["bonds"][ax] = rows

    # --- 3. the contraction, labelled ----------------------------------------------------------
    Zp = np.einsum("ab,nbc,dc->nad", L, born.Z, L)[perm]  # our Born tensors, provider frame, their order
    labels = [born.labels[j] for j in perm]
    with FITTED_VALENCE.applied():
        clamped = clamped_ion(ref, 0.01)
    tb = their_born(born_path, rec)
    print("\n# === CROSS-HAMILTONIAN DIAGNOSTIC, NOT A RESULT: their displacements per unit strain through our Born-consistent")
    print("# charges, (1/V0) sum_k Z_k . J_k, C/m^2 -- the internal-strain dipole *their* kinematics would produce in *our*")
    print("# charge model.  Provisional CP2K geometry, our charges; nothing fitted; not a piezoelectric coefficient.")
    print(f"# Our Born tensor at this structure: |asr| {np.abs(born.asr).max():.1e} e.  V0 = theirs ({V0_theirs:.2f} A^3) for their J,")
    print(f"# ours ({V0_ours:.2f}) for ours.  'clamped' is our clamped-ion proper coefficient (every atom affine, nothing relaxed),")
    print("# and total = clamped + internal is the closure check on the same two states.")
    print(f"# Polar (y) sign: our P_y = {P_prov[1]:+.4f}; a positive entry increases P_y, so 'along P' = same sign as P_y.")
    print(f"# Target for the transverse response (docs/ELECTROMECHANICS.md 5.9): |e_y,xx| ~ {TARGET[0]}-{TARGET[1]} C/m^2 against our 0.04.")
    out["contraction"] = {}
    print(f"  {'strain':7s} {'ours Z . their J (1%)':>27s} {'(2%)':>27s} {'ours Z . our J (1%)':>27s} "
          f"{'our clamped-ion e':>27s} {'our proper e (total)':>27s} {'closure':>9s}")
    for ax in AXES:
        ct1 = contraction(Zp, rec.J[(ax, 0.01)], V0_theirs)
        ct2 = contraction(Zp, rec.J[(ax, 0.02)], V0_theirs)
        co = contraction(Zp, ours[(ax, 0.01)]["J"], V0_ours)
        et = ours[(ax, 0.01)]["e_total"]
        ec = L @ clamped[PACKER_VOIGT[ax]]
        print(f"  {ax + ax:7s} {fmt3(ct1, 8, 4)} {fmt3(ct2, 8, 4)} {fmt3(co, 8, 4)} {fmt3(ec, 8, 4)} {fmt3(et, 8, 4)} "
              f"{float(np.abs(et - ec - co).max()):9.1e}")
        entry = {"ours_Z_their_J_1pc": ct1.tolist(), "ours_Z_their_J_2pc": ct2.tolist(), "ours_Z_our_J_1pc": co.tolist(),
                 "our_clamped_e_1pc": ec.tolist(), "our_proper_e_1pc": et.tolist(),
                 "by_type_y": {"ours_Z_their_J": by_type(Zp, rec.J[(ax, 0.01)], labels, V0_theirs),
                               "ours_Z_our_J": by_type(Zp, ours[(ax, 0.01)]["J"], labels, V0_ours)}}
        if tb is not None:
            entry["their_Zasr_their_J_1pc"] = contraction(tb["corrected"], rec.J[(ax, 0.01)], V0_theirs).tolist()
            entry["their_Zraw_their_J_1pc"] = contraction(tb["raw"], rec.J[(ax, 0.01)], V0_theirs).tolist()
            entry["by_type_y"]["their_Zasr_their_J"] = by_type(tb["corrected"], rec.J[(ax, 0.01)], labels, V0_theirs)
        out["contraction"][ax] = entry
    if tb is not None:
        print(f"# the same J through the provider's own DFPT tensor ({tb['protocol']}) -- provisional on both sides and not")
        print("# same-Hamiltonian (QE Born, CP2K geometry), so a bound on how much the choice of charges matters, not a result:")
        for ax in AXES:
            e = out["contraction"][ax]
            print(f"  {ax + ax:7s} their Z(asr) . their J (1%) {fmt3(e['their_Zasr_their_J_1pc'], 8, 4)}   "
                  f"Z(raw) {fmt3(e['their_Zraw_their_J_1pc'], 8, 4)}")
    print("# polar (y) component by atom type, C/m^2:")
    for ax in AXES:
        bt = out["contraction"][ax]["by_type_y"]
        print(f"  {ax + ax:7s} " + "  ".join(
            f"{k}: " + ", ".join(f"{lab} {bt[k][lab]:+.3f}" for lab in bt[k]) for k in bt))
    # the Born-charge reading of the TOTAL atomic motion: affine + internal.  For rigid chains pinned to
    # their sites it vanishes identically (the acoustic sum rule holds chain by chain), so whatever is
    # left of a proper coefficient is the part of the clamped-ion response Born charges do not describe:
    # in our model the induced dipoles' response to the strained lattice, in DFT the electronic term.
    print("# Born-charge reading of the TOTAL atomic motion, polar (y) component, C/m^2: affine + internal = total.  Zero by the")
    print("# acoustic sum rule for rigid chains; what remains of any proper coefficient is the electronic part of the clamped-ion")
    print("# response that Born charges do not describe (in our model the induced dipoles' lattice response; in DFT, unmeasured).")
    r_ours = rp0[perm]
    for ax in AXES:
        a = AXES.index(ax)
        aff_o = np.zeros_like(r_ours)
        aff_o[:, a] = r_ours[:, a]
        aff_t = np.zeros_like(rp_theirs)
        aff_t[:, a] = rp_theirs[:, a]
        bo = contraction(Zp, aff_o, V0_ours)[1]
        bt_o = contraction(Zp, aff_t, V0_theirs)[1]
        io, it = out["contraction"][ax]["ours_Z_our_J_1pc"][1], out["contraction"][ax]["ours_Z_their_J_1pc"][1]
        line = (f"  {ax + ax:7s} ours Z, our motion {bo:+.4f} {io:+.4f} = {bo + io:+.4f};  ours Z, their motion "
                f"{bt_o:+.4f} {it:+.4f} = {bt_o + it:+.4f}")
        out["contraction"][ax]["born_affine_y"] = {"ours_Z_our_geometry": float(bo), "ours_Z_their_geometry": float(bt_o)}
        out["contraction"][ax]["born_total_y"] = {"ours_Z_our_motion": float(bo + io), "ours_Z_their_motion": float(bt_o + it)}
        if tb is not None:
            bt_t = contraction(tb["corrected"], aff_t, V0_theirs)[1]
            itt = out["contraction"][ax]["their_Zasr_their_J_1pc"][1]
            line += f";  their Z(asr), their motion {bt_t:+.4f} {itt:+.4f} = {bt_t + itt:+.4f}"
            out["contraction"][ax]["born_affine_y"]["their_Zasr_their_geometry"] = float(bt_t)
            out["contraction"][ax]["born_total_y"]["their_Zasr_their_motion"] = float(bt_t + itt)
        print(line)
    ey, eo = out["contraction"]["x"]["ours_Z_their_J_1pc"][1], out["contraction"]["x"]["ours_Z_our_J_1pc"][1]
    print(f"# polar response to the long-axis strain: their J through our charges gives an internal-strain term of {ey:+.4f} C/m^2 "
          f"({'inside' if TARGET[0] <= abs(ey) <= TARGET[1] else 'outside'} the {TARGET[0]}-{TARGET[1]} band; "
          f"{'along' if ey * P_prov[1] > 0 else 'against'} P), our own rigid kinematics {eo:+.4f}: the two displacement "
          f"patterns differ by {rms(rec.J[('x', 0.01)] - ours[('x', 0.01)]['J']):.3f} A per unit strain rms and by "
          f"{ey - eo:+.4f} C/m^2 of polar dipole.  Our clamped-ion e_y,xx is {out['contraction']['x']['our_clamped_e_1pc'][1]:+.4f}; "
          f"the total {out['contraction']['x']['our_proper_e_1pc'][1]:+.4f}.")

    # --- 4. 1% versus 2% --------------------------------------------------------------------------
    print("\n# === 1% versus 2%: |J(2%) - J(1%)| / |J(1%)| (Frobenius) and the largest component change, A per unit strain ===")
    out["amplitude"] = {}
    for ax in AXES:
        J1, J2 = ours[(ax, 0.01)]["J"], ours[(ax, 0.02)]["J"]
        n1 = float(np.linalg.norm(J1))
        chg = float(np.linalg.norm(J2 - J1) / n1) if n1 > 1e-10 else float("nan")
        mx = float(np.abs(J2 - J1).max())
        t1, t2 = rec.J[(ax, 0.01)], rec.J[(ax, 0.02)]
        tchg = rec.amplitude_change[ax]
        tmx = float(np.abs(t2 - t1).max())
        # where their change sits: the within-chain remainder at the two amplitudes
        a = AXES.index(ax)
        Rt = np.zeros_like(t1) if ax == "z" else rigid_prediction(rec.elements, rp_theirs, a, n_chain)
        rem1, rem2 = t1 - Rt, t2 - Rt
        what = "their J rms" if ax == "z" else "their within-chain remainder rms"
        print(f"  {ax}{ax}: ours {100 * chg:.3f}% (max {mx:.1e});  theirs {100 * tchg:.3f}% (max {tmx:.3f}); "
              f"{what} {rms(rem1):.4f} -> {rms(rem2):.4f}; "
              f"our shape parameter at +/-1%: {states[(ax, 0.01)][2].relaxed.x} / {states[(ax, -0.01)][2].relaxed.x}, "
              f"+/-2%: {states[(ax, 0.02)][2].relaxed.x} / {states[(ax, -0.02)][2].relaxed.x} deg")
        out["amplitude"][ax] = {"ours": chg, "ours_max": mx, "theirs": tchg, "theirs_max": tmx,
                                "shape_x": {f"{s:+.2f}": states[(ax, s)][2].relaxed.x.tolist() for s in (0.01, -0.01, 0.02, -0.02)}}
        for h in AMPLITUDES:
            for sgn in (1.0, -1.0):
                st = states[(ax, sgn * h)][2]
                out["amplitude"][ax][f"state_{sgn * h:+.2f}"] = {"iterations": st.relaxed.iterations, "residual": st.relaxed.residual,
                                                                 "c_error": st.relaxed.c_error, "energy": st.energy}

    if args.json:
        with open(args.json, "w") as fh:
            json.dump(out, fh, indent=1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
