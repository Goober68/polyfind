"""Far-infrared activity of alpha-PVDF's Gamma-point modes, on the model's own Born charges, against the measured lines.

    D=deliverables/alpha_far_ir
    python examples/alpha_far_ir.py --cell antipolar --out $D/antipolar_shared.json      # ~10 min each
    python examples/alpha_far_ir.py --cell polar --out $D/polar_shared.json
    PYTHONPATH=<tree with torsion>/src python examples/alpha_far_ir.py --cell antipolar --out $D/antipolar_torsion.json
    PYTHONPATH=<tree with torsion>/src python examples/alpha_far_ir.py --cell polar --out $D/polar_torsion.json
    python examples/alpha_far_ir.py --compare $D/antipolar_shared.json $D/antipolar_torsion.json \
        $D/polar_shared.json $D/polar_torsion.json

Each computing line takes one Hamiltonian -- whichever ``polyfind`` is on ``PYTHONPATH`` -- and
one cell, and writes a record; ``--compare`` reads records, pairs the two Hamiltonians of each
cell side by side (modes matched by eigenvector overlap) and sets every IR-active mode against
the measured lines.  The construction is ``examples/crystal_phonons.py``'s stage (3) for alpha
(TGTG'): ``pvdf-dft-valence-flux-born``, induced dipoles, Ewald, cell and shape relaxed on the
deformable path, stretch minima pinned to the built bond lengths, every atom relaxed at fixed
cell, the Hessian a central difference of the analytic all-atom gradient, translations
projected.  At that geometry the Born tensors are the mixed field-displacement derivative of
the same energy (:func:`polyfind.infrared.born_charges_by_field`), checked against central
differences of the cell dipole and against the acoustic sum rule; each mode's dipole
derivative, intensity and polarization follow (:func:`polyfind.infrared.infrared`), with the
cell's space group found from the relaxed coordinates and every mode's characters under it.

**Two cells, because the packer's alpha is not alpha's space group.**  ``--cell polar`` is the
packer's own TGTG' reference, the one ``docs/PHONONS.md`` tabulates: chain 2 is flipped, but the
two chains' transverse dipoles align (``DESIGN.md``'s correction) and the cell is ``P2_1cn``,
point group ``C2v``, polar -- the space group of delta-PVDF (form IV).  ``--cell antipolar`` puts
the same relaxed chain in the flip-and-equal-setting-angle subspace of
:func:`polyfind.fitting.antipolar_cell`, polishes it there and relaxes cell and shape on the same
deformable path; it comes back ``P2_1/c`` (``C2h``, 2_1 along the long axis b), alpha's own
group, and is the like-for-like cell for alpha's measured lines.

The tree is identified by what its packer offers: a ``CrystalPacker.placed_energy_and_grad``
is the placed-cell owner with real Cartesian torsion in the Hessian; without it the torsion
term is the constant of the nominal dihedrals (``docs/PHONONS.md``).  Nothing is fitted or
tuned here; the measured lines are carried with their sources and only compared.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# Measured far-infrared / THz lines of alpha-PVDF (form II) below 200 cm^-1, as read in the
# sources; docs/PHONONS.md (section of 2026-10-08) carries the full citations and what could and
# could not be read.  RJ73 = Rabolt & Johnson, J. Chem. Phys. 59, 3710 (1973), abstract only;
# CP16 = Chamorro-Posada, arXiv:1604.03919 section 3.2, full text, quoting RJ73 and Mori et al.,
# J. Mol. Struct. 1090, 93 (2015), whose full text could not be read.  "n.s." = the abstract
# reports spectra from about 13 to 373 K without saying which temperature a line was read at.
# Never filled from the model.
MEASURED: list = [
    {"cm1": 53.0, "T": "RT", "what": "IR; rotatory lattice mode (RJ73 abs.; CP16)"},
    {"cm1": 60.0, "T": "90 K", "what": "the 53 line cooled from 300 K (CP16, citing RJ73)"},
    {"cm1": 78.93, "T": "n.s.", "what": "THz-TDS 'third resonance' (CP16, citing Mori 2015)"},
    {"cm1": 85.0, "T": "n.s.", "what": "IR; intramolecular, single TGTG' chain (RJ73 abs.)"},
    {"cm1": 100.0, "T": "RT", "what": "IR main line (CP16, citing RJ73); '102' in RJ73's abstract"},
    {"cm1": 175.0, "T": "n.s.", "what": "IR; intramolecular, single TGTG' chain (RJ73 abs.)"},
]


def _sha(path: str) -> str:
    with open(path, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def _provenance() -> dict:
    import polyfind

    src = os.path.dirname(os.path.abspath(polyfind.__file__))
    root = os.path.dirname(os.path.dirname(src))

    def git(*a):
        try:
            return subprocess.run(["git", "-C", root, *a], capture_output=True, text=True, timeout=30).stdout.strip()
        except Exception:  # noqa: BLE001 - provenance is best effort
            return ""

    files = {f: _sha(os.path.join(src, f)) for f in ("phonon.py", "pack.py", "born.py", "infrared.py", "forcefield.py")
             if os.path.exists(os.path.join(src, f))}
    return {"polyfind": src, "git_head": git("rev-parse", "HEAD"), "git_branch": git("rev-parse", "--abbrev-ref", "HEAD"),
            "uncommitted_paths": len([ln for ln in git("status", "--porcelain").splitlines() if ln.strip()]),
            "source_sha256": files}


def dipole_function(packer, params, latn):
    """The cell dipole ``sum_j q_j r_j + sum_j p_j`` (e.A) at placed coordinates, charges re-derived per chain.

    The independent route for the check on :func:`polyfind.infrared.born_charges_by_field`: the
    packer's own ``placed_dipole`` where the tree has one, else the chain-frame flux charges and
    the induced-dipole solve exactly as :func:`polyfind.phonon.cell_energy_and_grad` builds them.
    """
    if hasattr(packer, "placed_dipole"):
        return lambda P: np.sum(packer.placed_dipole(P, latn, float(params[6])), axis=0)
    from polyfind.born import chain_frame_coords

    n, nc, cz = packer.n, packer.n_chains, float(latn[2, 2])
    flux = getattr(packer, "_flux", None)

    def mu(P):
        if flux is not None:
            q = np.concatenate([np.asarray(flux.charges(chain_frame_coords(params, latn, P[s * n:(s + 1) * n], s), cz)[0],
                                           dtype=float) for s in range(nc)])
        else:
            q = np.asarray(packer._q_cell, dtype=float)
        m = q @ P
        if packer.polarizable is not None:
            m = m + packer._polarize(P, q, latn, float(params[6]))[1].sum(axis=0)
        return m

    return mu


def _infrared_module():
    """``polyfind.infrared``; for a tree on ``PYTHONPATH`` that predates the module (the torsion
    branch before this one is merged into it), the copy beside this example is put on the
    package's path and imported into that tree's ``polyfind`` -- every other module still comes
    from the tree under test, and the record says which file was used."""
    import importlib

    import polyfind

    try:
        return importlib.import_module("polyfind.infrared")
    except ImportError:
        here = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src", "polyfind")
        if not os.path.exists(os.path.join(here, "infrared.py")):
            raise
        polyfind.__path__.append(here)
        return importlib.import_module("polyfind.infrared")


def reference_and_shape(polymorph, preset: str):
    """``examples/crystal_phonons.py:reference``, step for step, also returning the :class:`~polyfind.mechanics.Reference`
    and :class:`~polyfind.mechanics.Shape` so that the antipolar cell can be relaxed on the same path."""
    import crystal_phonons as CP
    from polyfind import mechanics as M
    from polyfind.fitting import FITTED_VALENCE
    from polyfind.forcefield import SimpleFF
    from polyfind.lattice_table import clear_pair_table_cache
    from polyfind.polarizability import Polarizable

    val = SimpleFF.from_preset("pvdf-dft-valence")
    flx = SimpleFF.from_preset(preset)
    clear_pair_table_cache()
    polymer, seq, label = polymorph
    with FITTED_VALENCE.applied():
        ref, rr = M.refined_reference(polymer, seq, label=label, valence=val, refine_kw={"valence": val},
                                      charge_flux=flx, pack_kw={"table_cache_dir": None},
                                      polarizable=Polarizable(), **CP.EW)
        shape = M.shape_of(polymer, ref, angles=rr.angles)
        rel = M.relax_deformable(ref, shape, ref.params, shape.x0, free=M._CELL_FREE, c_target=None, maxiter=300)
        ref.packer.update_chain(rel.chain)
        ref = M.Reference(packer=ref.packer, params=rel.params, c=float(ref.packer.chain.c), label=label)
        ref, shape = M.relax_reference_deformable(ref, shape)
    return ref.packer, np.asarray(ref.params, dtype=float), ref, shape


def antipolar_reference(ref, shape, maxfev: int):
    """The model's antipolar TGTG' cell: chain 2 flipped *and* the two setting angles equal.

    For the alpha helix that subspace has exactly zero cell dipole (:func:`polyfind.fitting.antipolar_cell`,
    whose docstring says why); its centrosymmetric member is the model's analogue of alpha's
    ``P2_1/c``.  The polar reference's chain is placed there and the four free numbers
    ``(a, b, phi, dz)`` polished from ``antipolar_cell``'s two starts, then cell *and* chain shape are
    relaxed together on the deformable path exactly as the polar reference was
    (:func:`polyfind.mechanics.relax_reference_deformable`).  That relaxation frees ``phi1`` and
    ``phi2`` separately; the symmetric point is stationary by symmetry and the run reports whether
    it stayed there (the space group below is the test).
    """
    from types import SimpleNamespace

    from polyfind import mechanics as M
    from polyfind.fitting import antipolar_cell

    packer, params = ref.packer, np.asarray(ref.params, dtype=float)
    p_ap, _ = antipolar_cell(packer, SimpleNamespace(a=params[0], b=params[1], phi1=params[3]), maxfev=maxfev)
    ref_ap = M.Reference(packer=packer, params=p_ap, c=float(packer.chain.c), label=ref.label + " antipolar")
    ref_ap, _ = M.relax_reference_deformable(ref_ap, shape)
    return ref_ap.packer, np.asarray(ref_ap.params, dtype=float)


def compute(args) -> dict:
    import crystal_phonons as CP
    from polyfind import phonon as PH
    from polyfind.pack import CrystalPacker

    IR = _infrared_module()
    torsion = hasattr(CrystalPacker, "placed_energy_and_grad")
    label = ("placed-cell owner, Cartesian torsion in the Hessian" if torsion
             else "shared branch, torsion a constant of the nominal dihedrals")
    prov = _provenance()
    prov["infrared_module"] = os.path.abspath(IR.__file__)
    prov["infrared_sha256"] = _sha(IR.__file__)
    print(f"# Hamiltonian: {label}\n# polyfind {prov['polyfind']} @ {prov['git_head'][:10]} ({prov['git_branch']}, "
          f"{prov['uncommitted_paths']} uncommitted paths); infrared module {prov['infrared_module']}")
    t0 = time.time()
    packer, params, ref, shape = reference_and_shape(CP.ALPHA, args.preset)
    print(f"reference built in {time.time() - t0:.0f} s: a={params[0]:.4f} b={params[1]:.4f} gamma={params[2]:.2f} "
          f"c={packer.chain.c:.4f} A; phi1={params[3]:.2f} phi2={params[4]:.2f} dz={params[5]:.4f} flip={int(params[6])}")
    e_polar = float(packer.energy(params[None])[0]) / (packer.n_chains * packer.chain.n_monomers)
    print(f"  packer energy {e_polar:.5f} kcal/mol per monomer (fitted stretch r0, cell level)")
    e_antipolar = None
    if args.cell == "antipolar":
        t0 = time.time()
        packer, params = antipolar_reference(ref, shape, args.maxfev)
        e_antipolar = float(packer.energy(params[None])[0]) / (packer.n_chains * packer.chain.n_monomers)
        print(f"antipolar cell (flip, phi1 = phi2) relaxed over cell and shape in {time.time() - t0:.0f} s: a={params[0]:.4f} "
              f"b={params[1]:.4f} gamma={params[2]:.2f} c={packer.chain.c:.4f} A; phi1={params[3]:.3f} phi2={params[4]:.3f} "
              f"dz={params[5]:.4f} flip={int(params[6])}; packer energy {e_antipolar:.5f} kcal/mol per monomer, "
              f"{e_antipolar - e_polar:+.5f} against the polar reference")
    pk = PH.packer_with_built_bond_lengths(packer)
    e_cell, e_pack, rel_ka = PH.check_known_answer(pk, params)
    print(f"known answer (stretch minima pinned): cell_energy {e_cell:.10f} vs packer.energy {e_pack:.10f}, rel {rel_ka:.1e}")
    P0, latn = PH.placed_coordinates(pk, params)
    rel = PH.relax_all_atom(pk, params, Pn=P0, latn=latn)
    print("relaxation " + rel.summary())
    if not rel.converged:
        raise RuntimeError("the all-atom relaxation was not accepted; a Hessian here would not be at a stationary point")
    Pn = rel.Pn
    print(f"max displacement from the reference {np.abs(Pn - P0).max():.4f} A")
    t0 = time.time()
    ph = PH.phonons_gamma(pk, params, h=args.h, Pn=Pn, latn=latn, asr="project")
    print(f"Hessian in {time.time() - t0:.0f} s; asymmetry {ph.hessian.asymmetry:.1e}; ASR residual "
          f"{ph.asr_max / ph.hessian_scale:.1e} of the largest entry; imaginary {ph.n_imaginary}, zero {ph.n_zero}")

    # --- Born tensors, two routes, and the checks ---------------------------------------------
    Z = IR.born_charges_by_field(pk, params, Pn, latn, step=args.field_step)
    Z2 = IR.born_charges_by_field(pk, params, Pn, latn, step=10 * args.field_step)
    Zd = IR.born_charges_by_dipole(dipole_function(pk, params, latn), Pn, h=1e-4)
    Zd2 = IR.born_charges_by_dipole(dipole_function(pk, params, latn), Pn, h=5e-5)
    asr = Z.sum(axis=0)
    print(f"Born tensors: field route vs dipole route max |dZ| {np.abs(Z - Zd).max():.1e} e; field step x10 moves Z by "
          f"{np.abs(Z - Z2).max():.1e}; dipole step /2 by {np.abs(Zd - Zd2).max():.1e}; largest |Z| {np.abs(Z).max():.3f} e")
    print("acoustic sum rule sum_i Z*_i (e):\n" + np.array2string(asr, precision=2, floatmode="maxprec", suppress_small=False))
    eps = IR.electronic_dielectric(pk, params, Pn, latn)
    print("clamped-ion dielectric tensor (x polar, y long lateral, z chain):\n" + np.array2string(eps, precision=4))
    V = abs(float(np.linalg.det(latn)))
    mu = np.asarray(dipole_function(pk, params, latn)(Pn), dtype=float)
    print(f"cell dipole at the relaxed geometry (charges + induced) {np.array2string(mu, precision=5)} e.A; "
          f"|P| = {np.linalg.norm(mu) / V * 16.0217663:.4f} C/m^2 (V = {V:.3f} A^3)")

    # --- symmetry -------------------------------------------------------------------------------
    elements = list(pk.elements) * pk.n_chains
    ops = IR.space_group(Pn, latn, elements, tol=args.sym_tol)
    ops0 = IR.space_group(P0, latn, elements, tol=args.sym_tol)
    print(f"space group of the relaxed cell (atom-match tolerance {args.sym_tol:g} A): point group "
          f"{IR.point_group_name(ops)}; of the unrelaxed reference: {IR.point_group_name(ops0)}")
    for op in ops:
        print("  " + op.describe())
    ir = IR.infrared(ph, Z, ops)
    acoustic = list(ph.acoustic)
    ac_rel = max(ir.modes[i].intensity for i in acoustic) / max(m.intensity for m in ir.modes)
    print(f"Born tensors' own symmetry error max |Z[g i] - R Z[i] R^T| {ir.born_symmetry:.1e} e; largest intensity on an "
          f"inversion-even optical mode {ir.g_leak(acoustic):.1e} of the strongest; on the acoustic modes {ac_rel:.1e}; "
          f"along a symmetry-forbidden axis {ir.forbidden_leak(acoustic):.1e}")
    chars = ir.chi
    if chars is not None:
        optical = [k for k in range(chars.shape[1]) if k not in acoustic]  # the projected translations are degenerate
        dev = float(np.abs(np.abs(chars[:, optical]) - 1.0).max())
        print(f"characters: max ||chi| - 1| over every optical mode and operation {dev:.1e}")

    strongest = max(ir.modes, key=lambda m: m.intensity if m.index not in acoustic else -1)
    print(f"\nstrongest mode of the whole spectrum: {strongest.freq_cm1:.1f} cm^-1 {strongest.irrep}, "
          f"{strongest.km_per_mol:.1f} km/mol per cell")
    print(f"\nmodes up to {args.fmax:g} cm^-1 (relative to the strongest optical mode in that window):")
    print(ir.table(args.fmax, acoustic=acoustic))

    modes = []
    for m in ir.modes:
        if m.index in acoustic or m.freq_cm1 > args.fstore:
            continue
        md = ph.mode(m.index)
        lo = IR.lo_shift(m.freq_cm1, m.dipole, eps, V)
        modes.append({"index": m.index, "freq_cm1": round(m.freq_cm1, 3), "irrep": m.irrep, "allowed": m.allowed,
                      "parity": None if m.parity is None else round(m.parity, 6),
                      "km_per_mol": m.km_per_mol, "dipole_e_per_sqrt_amu": [float(v) for v in m.dipole],
                      "dipole_shares_xyz": [round(float(v), 4) for v in m.shares],
                      "lo_estimate_cm1": round(lo, 2),
                      "eigvec_shares_xyz": [round(float(v), 4) for v in md.axis], "rigid_chain": round(md.rigid_chain, 4),
                      "by_type": {k: round(v, 4) for k, v in md.by_type.items()},
                      "eigenvector": [round(float(v), 8) for v in ph.modes[:, m.index]]})
    print("\nfirst-order LO estimate for the IR-active modes (wavevector along the mode's own dipole; "
          "not in the Hessian, see polyfind.infrared.lo_shift):")
    for r in modes:
        if r["freq_cm1"] <= args.fmax and r["allowed"] != "":
            print(f"  TO {r['freq_cm1']:7.2f} -> LO {r['lo_estimate_cm1']:7.2f} cm^-1  ({r['irrep']}, {r['km_per_mol']:.3f} km/mol)")
    return {"hamiltonian": label, "torsion_in_hessian": torsion, "cell": args.cell, "provenance": prov, "preset": args.preset,
            "packer_energy_polar_per_monomer": e_polar, "packer_energy_antipolar_per_monomer": e_antipolar,
            "relaxed_energy_per_monomer_pinned": float(rel.energy) / (pk.n_chains * pk.chain.n_monomers),
            "h": args.h, "field_step": args.field_step,
            "params": [float(v) for v in params], "c": float(pk.chain.c), "volume_A3": V,
            "known_answer_rel": rel_ka, "relaxation": rel.summary(), "max_displacement_A": float(np.abs(Pn - P0).max()),
            "hessian_asymmetry": ph.hessian.asymmetry, "n_imaginary": ph.n_imaginary,
            "freq_cm1": [round(float(v), 3) for v in ph.freq_cm1],
            "born_field_vs_dipole_max_e": float(np.abs(Z - Zd).max()), "born_asr_e": asr.tolist(),
            "born_asr_max_e": float(np.abs(asr).max()), "born_symmetry_error_e": ir.born_symmetry,
            "eps_inf": eps.tolist(), "point_group": ir.point_group, "point_group_unrelaxed": IR.point_group_name(ops0),
            "symmetry_ops": [op.describe() for op in ops], "g_leak": ir.g_leak(acoustic), "acoustic_leak": ac_rel,
            "forbidden_leak": ir.forbidden_leak(acoustic), "cell_dipole_eA": mu.tolist(),
            "strongest": {"freq_cm1": strongest.freq_cm1, "irrep": strongest.irrep, "km_per_mol": strongest.km_per_mol},
            "elements": elements, "positions_A": [[round(float(x), 8) for x in r] for r in Pn],
            "lattice_A": latn.tolist(), "born_charges_e": np.round(Z, 6).tolist(), "modes": modes}


# ------------------------------------------------------------------------- the comparison
def _load(path: str) -> dict:
    with open(path) as fh:
        return json.load(fh)


def _match(a: dict, b: dict) -> dict:
    """``index in a -> (index in b, overlap |e_a . e_b|^2)``, a one-to-one assignment on eigenvector overlap."""
    from scipy.optimize import linear_sum_assignment

    Ea = np.array([m["eigenvector"] for m in a["modes"]]).T
    Eb = np.array([m["eigenvector"] for m in b["modes"]]).T
    O = (Ea.T @ Eb) ** 2
    r, c = linear_sum_assignment(-O)
    return {a["modes"][i]["index"]: (b["modes"][j]["index"], float(O[i, j])) for i, j in zip(r, c)}


def _pol(m: dict) -> str:
    x, y, z = m["dipole_shares_xyz"]
    return f"x{x:.2f} y{y:.2f} z{z:.2f}"


def _tag(r: dict) -> str:
    return f"{r['cell']}/{'torsion' if r['torsion_in_hessian'] else 'shared'}"


def _rigid(rec: dict, m: dict) -> tuple:
    """``(translation share, libration share, "co"/"counter")`` of a stored mode's mass-weighted eigenvector.

    Per chain, the three uniform translations and the rotation about the chain's own axis through
    its mass centre (orthogonal to them in the mass-weighted metric), as :class:`polyfind.phonon.Mode`
    builds its rigid-chain basis -- split here so that a libration can be told from a translation.
    ``co`` when the two chains turn the same way about the cell's ``+z``, ``counter`` otherwise
    (``-`` when the libration share is below 0.05)."""
    from polyfind.pack import MASS

    E = np.asarray(m["eigenvector"], dtype=float).reshape(-1, 3)
    P = np.asarray(rec["positions_A"], dtype=float)
    w = np.sqrt(np.array([MASS[e] for e in rec["elements"]]))
    n = len(w) // 2
    trans = rot = 0.0
    coef = []
    for s in range(2):
        sl = slice(s * n, (s + 1) * n)
        com = (w[sl, None] ** 2 * P[sl]).sum(axis=0) / (w[sl] ** 2).sum()
        for d in range(3):
            v = np.zeros((n, 3))
            v[:, d] = w[sl]
            trans += float((E[sl] * v).sum() / np.linalg.norm(v)) ** 2
        r = P[sl] - com
        v = np.stack([-r[:, 1], r[:, 0], np.zeros(n)], axis=1) * w[sl, None]
        c = float((E[sl] * v).sum() / np.linalg.norm(v))
        rot += c * c
        coef.append(c)
    sense = "-" if rot < 0.05 else ("co" if coef[0] * coef[1] > 0 else "counter")
    return trans, rot, sense


def _side_by_side(a: dict, b: dict, fmax: float) -> None:
    mt = _match(a, b)
    bi = {m["index"]: m for m in b["modes"]}
    win = [m for m in a["modes"] if m["freq_cm1"] <= fmax]
    ref_a = max(m["km_per_mol"] for m in win)
    ref_b = max(m["km_per_mol"] for m in b["modes"] if m["freq_cm1"] <= fmax)
    print(f"\n{a['cell']} cell: modes up to {fmax:g} cm^-1, A = {_tag(a)}, B = {_tag(b)}, matched by eigenvector overlap "
          f"(ovl); intensity relative to each one's strongest mode in the window ({ref_a:.3f} / {ref_b:.3f} km/mol per "
          "cell); pol = shares of |d|^2 along x/y/z (z = chain); character (A) = eigenvector shares x/y/z, rigid-chain "
          "(RC, of which translation T and libration L about the chains' own axes, co/counter-rotating), F share")
    print(f"{'A cm^-1':>8s} {'irr':>3s} {'rel':>5s} {'pol (A)':>15s} | {'B cm^-1':>8s} {'irr':>3s} {'rel':>5s} "
          f"{'pol (B)':>15s} | {'ovl':>4s} {'shift':>6s} | character (A)")
    for m in win:
        j, ov = mt[m["index"]]
        n = bi[j]
        e = m["eigvec_shares_xyz"]
        tr, lib, sense = _rigid(a, m)
        print(f"{m['freq_cm1']:8.2f} {m['irrep']:>3s} {m['km_per_mol'] / ref_a:5.3f} {_pol(m):>15s} | "
              f"{n['freq_cm1']:8.2f} {n['irrep']:>3s} {n['km_per_mol'] / ref_b:5.3f} {_pol(n):>15s} | {ov:4.2f} "
              f"{n['freq_cm1'] - m['freq_cm1']:+6.1f} | {e[0]:.2f}/{e[1]:.2f}/{e[2]:.2f} RC {m['rigid_chain']:.2f} "
              f"(T {tr:.2f} L {lib:.2f} {sense}) F {m['by_type'].get('F', 0.0):.2f}")


def compare(paths, fmax: float, active: float) -> None:
    """Every record's header; side-by-side tables of the two Hamiltonians for each cell; the measured lines."""
    recs = [_load(p) for p in paths]
    for r in recs:
        pv = r["provenance"]
        mu = r.get("cell_dipole_eA")
        print(f"{_tag(r):18s} {r['hamiltonian']} ({pv['git_branch']} @ {pv['git_head'][:10]}, {pv['uncommitted_paths']} "
              f"uncommitted paths): point group {r['point_group']}; imaginary {r['n_imaginary']}; Born ASR "
              f"{r['born_asr_max_e']:.1e} e, field vs dipole {r['born_field_vs_dipole_max_e']:.1e} e, symmetry "
              f"{r['born_symmetry_error_e']:.1e} e; forbidden-axis leak {r.get('forbidden_leak', float('nan')):.1e}; "
              f"eps_inf " + "/".join(f"{r['eps_inf'][i][i]:.3f}" for i in range(3))
              + ("" if mu is None else f"; cell dipole ({mu[0]:+.4f}, {mu[1]:+.4f}, {mu[2]:+.4f}) e.A")
              + f"; E {r['relaxed_energy_per_monomer_pinned']:.5f} kcal/mol per monomer (pinned, all-atom)")
    by_cell: dict = {}
    for r in recs:
        by_cell.setdefault(r["cell"], []).append(r)
    for cell, rs in by_cell.items():
        rs.sort(key=lambda r: r["torsion_in_hessian"])
        if len(rs) == 2:
            _side_by_side(rs[0], rs[1], fmax)
    print(f"\nthe IR-active modes (relative intensity >= {active:g} of the strongest up to {fmax:g} cm^-1), strongest first, "
          "each against the nearest measured line (model: 0 K harmonic TO; see docs/PHONONS.md for every caveat)")
    for r in recs:
        win = [m for m in r["modes"] if m["freq_cm1"] <= fmax]
        ref = max(m["km_per_mol"] for m in win)
        act = sorted((m for m in win if m["km_per_mol"] / ref >= active), key=lambda m: -m["km_per_mol"])
        print(f"  {_tag(r)}:")
        for m in act:
            near = min(MEASURED, key=lambda ln: abs(ln["cm1"] - m["freq_cm1"])) if MEASURED else None
            pol = "par" if m["dipole_shares_xyz"][2] > 0.5 else "perp"
            out = f"    {m['freq_cm1']:7.2f} {m['irrep']:>3s} rel {m['km_per_mol'] / ref:5.3f} {pol:>4s} chain ({_pol(m)})"
            if near is not None:
                out += f"  nearest measured {near['cm1']:g} cm^-1 ({near['T']}): model - measured {m['freq_cm1'] - near['cm1']:+6.1f}"
            print(out)
    if MEASURED:
        print("\nmeasured lines (sources in docs/PHONONS.md):")
        for ln in MEASURED:
            print(f"  {ln['cm1']:7.1f} cm^-1  {ln['T']:>6s}  {ln['what']}")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", help="write this Hamiltonian's record (JSON) here")
    ap.add_argument("--compare", nargs="+", metavar="RECORD", help="print records side by side (pairs by cell)")
    ap.add_argument("--preset", default="pvdf-dft-valence-flux-born")
    ap.add_argument("--cell", choices=("polar", "antipolar"), default="polar",
                    help="the packer's own alpha reference (polar; docs/PHONONS.md's) or the antipolar flip-and-equal-angle cell")
    ap.add_argument("--maxfev", type=int, default=250, help="Nelder-Mead evaluations per start for the antipolar polish")
    ap.add_argument("--h", type=float, default=1e-3, help="Hessian finite-difference step, A")
    ap.add_argument("--field-step", type=float, default=1e-3, help="applied-field step for the Born tensors, V/A")
    ap.add_argument("--sym-tol", type=float, default=1e-3, help="atom-match tolerance for the space group, A")
    ap.add_argument("--fmax", type=float, default=200.0, help="report modes up to this wavenumber, cm^-1")
    ap.add_argument("--fstore", type=float, default=320.0, help="store eigenvectors of modes up to this wavenumber")
    ap.add_argument("--active", type=float, default=0.05, help="relative intensity counted as IR-active in the comparison")
    args = ap.parse_args(argv)
    if args.compare:
        compare(args.compare, fmax=args.fmax, active=args.active)
        return 0
    rec = compute(args)
    if args.out:
        os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
        with open(args.out, "w") as fh:
            json.dump(rec, fh, indent=1)
        print(f"\nrecord written to {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
