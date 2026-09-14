"""Our clamped-ion piezoelectric columns for beta-PVDF, as a reproducible record for the producer.

    python examples/clamped_ion_columns.py                 # prints the record and writes deliverables/clamped_ion_columns/
    python examples/clamped_ion_columns.py --out DIR       # elsewhere

The producer's Berry-phase clamped-ion matrix (``docs/REFERENCE_DATA_REQUEST.md``, 2026-09-13) is
compared column by column against this model's clamped-ion response.  The producer asked that the
exact compared baseline geometry, cell and source hashes be retained with a reproducible
four-column record rather than as numbers in the exchange; this script is that record.

**The quantity.**  Every atom is displaced affinely with the cell (fixed fractional nuclei),
nothing is relaxed, the charge-flux charges and the induced dipoles are re-solved at the deformed
geometry, and the cell dipole per *reference* volume is central-differenced,
``(1/V0) dmu_i/de_J``.  That is the dipole-per-reference-volume derivative of
:func:`examples.internal_strain_jacobian.clamped_ion`.  The producer's reducer follows
Vanderbilt's proper piezoelectric tensor, ``e_ijk = dP_i/de_jk + delta_jk P_i - delta_ij P_k``
(J. Phys. Chem. Solids 61, 147 (2000), eqs. 15/24), which equals the dipole-per-reference-volume
derivative minus ``delta_ij P_k``; for an engineering shear ``gamma_jk`` (strain tensor
off-diagonal ``gamma/2``) the correction is ``-(delta_ij P_k + delta_ik P_j) / 2``.  Both
definitions are written out.  The ``P`` subtracted is the total polarization of the same dipole,
charge-flux charges plus induced dipoles; the charge-only and induced parts are recorded
separately as the producer asked.

**Frame.**  The packer's axes are ``x`` = polar, ``y`` = long lateral, ``z`` = chain; the
producer's are ``x`` = long lateral, ``y`` = polar, ``z`` = chain.  The map is the rotation by
90 degrees about the chain axis used by the Born and internal-strain comparisons
(:func:`examples.internal_strain_jacobian.frame_map`), and the polar sign is fixed the same way
(fluorines at ``+y`` of their carbon, so ``P_y < 0``, as in the producer's cell).  Voigt columns
are given in the producer's frame; the two shears that involve the chain axis (``yz``, ``xz``)
cannot be applied by this cell parametrisation (the chain axis is fixed along ``z``) and are
recorded as not expressible.

**What this is not.**  Not a fit, not a piezoelectric coefficient of the material, not a
same-Hamiltonian quantity, and not a comparison in itself: the producer's numbers carry their
``quantitatively_valid=false`` label and are not reproduced here.  The record carries the SHA-256
of every source file that defines the calculation, the preset names, the NumPy version, and the
SHA-256 of the placed geometry it was evaluated at.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))

from polyfind import mechanics as M  # noqa: E402

_spec = importlib.util.spec_from_file_location("internal_strain_jacobian", os.path.join(ROOT, "examples", "internal_strain_jacobian.py"))
ISJ = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ISJ)

PRESET = "pvdf-dft-valence-flux-born"
OUT = os.path.join(ROOT, "deliverables", "clamped_ion_columns")
SOURCES = ("examples/clamped_ion_columns.py", "examples/internal_strain_jacobian.py", "examples/fit_born_flux.py",
           "src/polyfind/mechanics.py", "src/polyfind/pack.py", "src/polyfind/born.py", "src/polyfind/forcefield.py",
           "src/polyfind/polarizability.py", "src/polyfind/fitting.py", "src/polyfind/ewald.py", "src/polyfind/polymers.py",
           "src/polyfind/chain.py")
PROVIDER_VOIGT = ("xx", "yy", "zz", "yz", "xz", "xy")
IDX = [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)]


def sha256_text(path: str) -> str:
    data = open(path, "rb").read().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return hashlib.sha256(data).hexdigest()


def cell_dipole(packer, params, Pn, latn) -> tuple:
    """``(mu_charges, mu_induced)`` in e.A for placed coordinates, charges re-solved per chain."""
    flux = getattr(packer, "_flux", None)
    n, nc = packer.n, packer.n_chains
    q = np.array(packer._q_cell, dtype=float)
    if flux is not None:
        for s in range(nc):
            sl = slice(s * n, (s + 1) * n)
            q[sl] = flux.charges(ISJ.chain_frame_coords(params, latn, Pn[sl], s), float(latn[2, 2]))[0]
    mu_q = q @ Pn
    mu_ind = np.zeros(3)
    if packer.polarizable is not None:
        mu_ind = packer._polarize(Pn, q, latn, float(params[6]))[1].sum(axis=0)
    return mu_q, mu_ind


def clamped_columns(ref, h: float) -> dict:
    """Packer-frame ``{K: (3,) C/m^2}`` for the expressible Voigt columns, dipole per reference volume."""
    packer = ref.packer
    P0, lat0 = ISJ.placed(ref, ref.params, packer.chain)
    out = {}
    for K in range(6):
        try:
            eps = np.zeros(6)
            eps[K] = h
            M.strained_cell(ref.params, ref.c, eps)
        except ValueError:
            continue
        mus = []
        for sgn in (1.0, -1.0):
            eps = np.zeros(6)
            eps[K] = sgn * h
            sc = M.strained_cell(ref.params, ref.c, eps)
            F = np.eye(3) + M.strain_tensor(eps)
            mu_q, mu_ind = cell_dipole(packer, sc.params, P0 @ F, lat0 @ F)
            mus.append(mu_q + mu_ind)
        out[K] = (mus[0] - mus[1]) / (2.0 * h) / ref.volume * ISJ.E_PER_A2_TO_C_PER_M2
    return out


def to_provider_columns(cols_pk: dict, L: np.ndarray, P_prov: np.ndarray) -> dict:
    """Map packer Voigt columns to the producer's frame and apply the Vanderbilt correction."""
    rec = {}
    for K, ev_pk in cols_pk.items():
        e = np.zeros(6)
        e[K] = 1.0
        eps_pv = L @ M.strain_tensor(e) @ L.T
        vals = [eps_pv[i, j] for i, j in IDX]
        J = int(np.argmax(np.abs(vals)))
        sgn = float(np.sign(vals[J]))
        ev = L @ ev_pk * sgn
        i_, j_ = IDX[J]
        corr = np.zeros(3)
        if i_ == j_:
            corr[i_] = P_prov[i_]
        else:
            corr[i_] += 0.5 * P_prov[j_]
            corr[j_] += 0.5 * P_prov[i_]
        rec[PROVIDER_VOIGT[J]] = {"packer_voigt": K, "sign_from_frame_map": sgn,
                                   "dipole_per_reference_volume_C_m2": ev.tolist(),
                                   "vanderbilt_proper_C_m2": (ev - corr).tolist(),
                                   "improper_lab_dP_de_C_m2": (ev - (P_prov if i_ == j_ else 0.0)).tolist() if i_ == j_ else None}
    for name in PROVIDER_VOIGT:
        rec.setdefault(name, "not expressible: this cell parametrisation fixes the chain axis along z and has no variable that tilts it")
    return rec


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--preset", default=PRESET)
    args = ap.parse_args(argv)

    packer, params, shape, rel = ISJ.beta_reference(ISJ.SimpleFF.from_preset(args.preset), ISJ.Polarizable())
    ref = M.Reference(packer=packer, params=params, c=float(packer.chain.c), label=ISJ.BETA[2])
    with ISJ.FITTED_VALENCE.applied():
        ref, shape = M.relax_reference_deformable(ref, shape)
        states = ISJ.our_states(ref, shape)
        cols = {h: clamped_columns(ref, h) for h in (0.0025, 0.01)}
        P0, lat0 = ISJ.placed(ref, ref.params, packer.chain)
        mu_q, mu_ind = cell_dipole(packer, ref.params, P0, lat0)
    n_chain = packer.n
    el = list(packer.elements) * packer.n_chains
    r0, c0, _ = states[None]
    L = ISJ.frame_map(el, r0, n_chain, float(ref.c))
    conv = ISJ.E_PER_A2_TO_C_PER_M2 / ref.volume
    Pq, Pind = L @ (mu_q * conv), L @ (mu_ind * conv)
    P = Pq + Pind
    rp, cp = ISJ.to_provider(P0, lat0, L)

    xyz_lines = [str(len(el)), f'Lattice="{" ".join(f"{v:.10f}" for v in np.asarray(cp).reshape(-1))}" Properties=species:S:1:pos:R:3 '
                 f"frame=producer(x=long,y=polar,z=chain) units=A"]
    xyz_lines += [f"{e:2s} {x:16.10f} {y:16.10f} {z:16.10f}" for e, (x, y, z) in zip(el, np.asarray(rp))]
    xyz = "\n".join(xyz_lines) + "\n"
    geometry_sha = hashlib.sha256(xyz.encode("ascii")).hexdigest()

    record = {
        "schema_version": 1,
        "purpose": "Polyfind clamped-ion piezoelectric columns of beta-PVDF for like-for-like comparison with the producer's Berry-phase clamped-ion matrix",
        "generator": "examples/clamped_ion_columns.py",
        "status": "model quantity on the fitted potential; not a fit, not a material coefficient, not same-Hamiltonian; the producer's values are not reproduced here",
        "calculation_owners_sha256": {s: sha256_text(os.path.join(ROOT, s)) for s in SOURCES},
        "hash_canonicalization": "UTF-8/ASCII text with CRLF and CR normalized to LF",
        "numpy_version": np.__version__,
        "presets": {"charge_flux": args.preset, "valence": "pvdf-dft-valence", "polarizable": "Polarizable() literature values", "coulomb": "ewald"},
        "reference_state": {
            "construction": "examples/fit_born_flux.py:beta_reference then mechanics.relax_reference_deformable, as examples/internal_strain_jacobian.py",
            "packer_params_a_b_gamma_phi1_phi2_dz_flip": np.asarray(ref.params, dtype=float).tolist(),
            "c_A": float(ref.c), "volume_A3": float(ref.volume),
            "lattice_rows_producer_frame_A": np.asarray(cp).tolist(),
            "frame_map_L_rows_producer_from_packer": L.tolist(),
            "geometry_xyz_sha256": geometry_sha, "geometry_file": "beta_pvdf_reference_producer_frame.xyz",
        },
        "polarization_producer_frame_C_m2": {"charge_only_with_flux": Pq.tolist(), "induced_dipoles": Pind.tolist(), "total": P.tolist(),
                                              "note": "the total is what the Vanderbilt correction subtracts; the clamped-ion dipoles include the induced part"},
        "definitions": {
            "dipole_per_reference_volume": "(1/V0) dmu_i/de_J, every atom affine with the cell, charges fluxed and induced dipoles re-solved, central difference at +/-h",
            "vanderbilt_proper": "dipole_per_reference_volume - delta_ij P_k (normal strains); - (delta_ij P_k + delta_ik P_j)/2 per unit engineering shear",
            "shear_convention": "engineering gamma; strain tensor off-diagonal = gamma/2 (mechanics.strain_tensor)",
            "voigt_order_producer_frame": list(PROVIDER_VOIGT),
        },
        "columns": {f"h={h}": to_provider_columns(cols[h], L, P) for h in cols},
    }
    os.makedirs(args.out, exist_ok=True)
    with open(os.path.join(args.out, "beta_pvdf_reference_producer_frame.xyz"), "w", newline="\n") as fh:
        fh.write(xyz)
    with open(os.path.join(args.out, "beta_pvdf_clamped_ion_columns.json"), "w", newline="\n") as fh:
        json.dump(record, fh, indent=1)
        fh.write("\n")

    print(f"reference: {record['reference_state']['packer_params_a_b_gamma_phi1_phi2_dz_flip']} c={ref.c:.4f} V0={ref.volume:.3f}; geometry sha256 {geometry_sha[:16]}...")
    print(f"P producer frame: charge-only {Pq.round(5)} induced {Pind.round(5)} total {P.round(5)} C/m^2")
    for h in cols:
        print(f"h = {h}:")
        for name in PROVIDER_VOIGT:
            c = record["columns"][f"h={h}"][name]
            if isinstance(c, str):
                print(f"  {name}: {c}")
            else:
                print(f"  {name}: (1/V0)dmu/de {np.round(c['dipole_per_reference_volume_C_m2'], 5)}  Vanderbilt proper {np.round(c['vanderbilt_proper_C_m2'], 5)}")
    print(f"written to {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
