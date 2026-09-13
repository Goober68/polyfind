"""Gamma-point phonons of PVDF polymorphs on the fitted potential: soft modes and the lowest optical modes.

    python examples/crystal_phonons.py                            # beta-PVDF, the 5.9 reference state
    python examples/crystal_phonons.py --polymorphs beta,alpha    # and alpha (TGTG'), same construction
    python examples/crystal_phonons.py --polymorphs gamma         # gamma (T3GT3G'), 48 atoms: minutes
    python examples/crystal_phonons.py --h 5e-4                   # finite-difference step (A)

The governing requirement is low loss and a high-frequency response, and the intrinsic proxy a
lattice model can add is the *stiffness of the transverse modes*: a polar crystal near an
instability has a soft transverse optical mode, and its frequency is what this script reports.
:mod:`polyfind.phonon` assembles the packer's own energy for every atom of the cell
independently, checks it against :meth:`~polyfind.pack.CrystalPacker.energy` (the known-answer
identity, printed first), takes the Cartesian Hessian by central differences of the analytic
all-atom gradient, mass-weights it and diagonalises.  The reference state is built exactly as
``examples/internal_strain_jacobian.py`` builds it: ``pvdf-dft-valence-flux-born``, induced
dipoles, Ewald, cell and shape relaxed on the deformable path.

Two Hessians are reported for each polymorph.  The first is at the reference state itself,
which is stationary over the packer's *own* variables (cell, setting angles, line-group shape)
but not necessarily over every atom independently -- the residual all-atom force is printed,
and a Hessian there can show imaginary modes that are the unrelaxed pendants, not a lattice
instability.  The second is after the atoms are let go at fixed cell
(:func:`polyfind.phonon.relax_all_atom`), which is the geometry a vibrational frequency belongs
to.  Both raw and acoustic-sum-projected results are printed; the projection changes nothing
but the three acoustic modes when the residual is small, and the residual is printed.

Read with the module docstring's caveat: the packer's Fourier torsion term is a constant of the
chain's nominal dihedrals, so the Hessian carries no torsional stiffness about the backbone
bonds.  The chain-twist modes are held by the nonbonded terms alone and are lower bounds.
"""
from __future__ import annotations

import argparse
import sys
import time

import numpy as np

from polyfind import mechanics as M
from polyfind.ewald import EwaldSpec
from polyfind.fitting import FITTED_VALENCE
from polyfind.forcefield import SimpleFF
from polyfind.lattice_table import clear_pair_table_cache
from polyfind.phonon import (SQRT_KCAL_A2_AMU_TO_CM1, SQRT_KCAL_A2_AMU_TO_THZ, check_known_answer, force_terms,
                             packer_with_built_bond_lengths, phonons_gamma, placed_coordinates, relax_all_atom)
from polyfind.polarizability import Polarizable
from polyfind.polymers import PVDF

T, GP, GM = 0, 1, 2
BETA = (PVDF, [T, T], "PVDF beta (TTTT)")
ALPHA = (PVDF, [T, GP, T, GM], "PVDF alpha (TGTG')")
GAMMA = (PVDF, [T, T, T, GP, T, T, T, GM], "PVDF gamma (T3GT3G')")
PRESET = "pvdf-dft-valence-flux-born"
EW = dict(coulomb="ewald", ewald=EwaldSpec())


def reference(polymorph, preset: str = PRESET):
    """``(packer, params)``: the polymorph packed, refined and relaxed over cell and shape, as
    ``examples/fit_born_flux.py:beta_reference`` and ``examples/internal_strain_jacobian.py`` do it."""
    val = SimpleFF.from_preset("pvdf-dft-valence")
    flx = SimpleFF.from_preset(preset)
    clear_pair_table_cache()
    polymer, seq, label = polymorph
    with FITTED_VALENCE.applied():
        ref, rr = M.refined_reference(polymer, seq, label=label, valence=val, refine_kw={"valence": val},
                                      charge_flux=flx, pack_kw={"table_cache_dir": None},
                                      polarizable=Polarizable(), **EW)
        shape = M.shape_of(polymer, ref, angles=rr.angles)
        rel = M.relax_deformable(ref, shape, ref.params, shape.x0, free=M._CELL_FREE, c_target=None, maxiter=300)
        ref.packer.update_chain(rel.chain)
        ref = M.Reference(packer=ref.packer, params=rel.params, c=float(ref.packer.chain.c), label=label)
        ref, shape = M.relax_reference_deformable(ref, shape)
    return ref.packer, np.asarray(ref.params, dtype=float)


def stage(packer, params, Pn, latn, h: float, n_lowest: int) -> None:
    t0 = time.time()
    ft = force_terms(packer, params, Pn, latn)
    print(f"residual all-atom force, max |dE/dP| kcal/(mol A): total {ft['total']:.3e}; valence bonds {ft['bonds']:.3e}, "
          f"valence angles {ft['angles']:.3e}, everything else {ft['other']:.3e}; "
          f"max |r - r0| {ft['bond_strain']:.4f} A, max |theta - theta0| {ft['angle_strain']:.2f} deg")
    raw = phonons_gamma(packer, params, h=h, Pn=Pn, latn=latn, asr="none", n_report=n_lowest)
    print(f"[raw]  ({time.time() - t0:.1f} s)")
    print(raw.report(n_lowest))
    proj = phonons_gamma(packer, params, h=h, Pn=Pn, latn=latn, asr="project", n_report=n_lowest)
    opt_raw, opt_proj = raw.optical, proj.optical
    print(f"[asr=project]  optical modes move by at most {np.abs(opt_raw - opt_proj).max():.3e} cm^-1; "
          f"acoustic now {' '.join(f'{proj.freq_cm1[i]:.2e}' for i in proj.acoustic)} cm^-1")
    print(f"lowest three optical (projected): " + ", ".join(f"{v:.2f}" for v in opt_proj[:3]) + " cm^-1")
    for m in proj.lowest_optical[:3]:
        print("  " + m.describe())


def report(packer, params, label: str, h: float, n_lowest: int = 6) -> None:
    print(f"\n=== {label}: a={params[0]:.4f} b={params[1]:.4f} gamma={params[2]:.2f} c={packer.chain.c:.4f} A; "
          f"phi1={params[3]:.2f} phi2={params[4]:.2f} dz={params[5]:.3f} flip={int(params[6])}; "
          f"{packer.N} atoms, {3 * packer.N} modes; {packer._ewald.describe()}")
    e_cell, e_pack, rel = check_known_answer(packer, params)
    print(f"known-answer check: cell_energy {e_cell:.10f} vs packer.energy {e_pack:.10f} kcal/mol, "
          f"relative difference {rel:.1e} (must be < 1e-9)")
    print(f"unit conversion: sqrt(kcal/(mol A^2 amu)) = {SQRT_KCAL_A2_AMU_TO_CM1:.4f} cm^-1 = {SQRT_KCAL_A2_AMU_TO_THZ:.5f} THz")
    P0, latn = placed_coordinates(packer, params)
    print("\n--- (1) at the reference state: stationary over the packer's variables (cell, setting angles, "
          "torsions, backbone angles), not over every atom")
    stage(packer, params, P0, latn, h, n_lowest)
    print("\n--- (2) every atom relaxed at fixed cell against the potential as fitted")
    Pn, E, gmax = relax_all_atom(packer, params, Pn=P0, latn=latn)
    print(f"E {E:.6f} kcal/mol (reference {e_pack:.6f}, lowered by {e_pack - E:.6f}), max |dE/dP| {gmax:.1e}, "
          f"max displacement {np.abs(Pn - P0).max():.4f} A -- the fitted stretch r0 values are not bond potentials "
          "(see polyfind.phonon.packer_with_built_bond_lengths), so this geometry has its C-F bonds stretched")
    stage(packer, params, Pn, latn, h, n_lowest)
    print("\n--- (3) every atom relaxed at fixed cell with the stretch r0 pinned to the built bond lengths "
          "(stiffness kept): the geometry stays at the reference's bond lengths")
    pk3 = packer_with_built_bond_lengths(packer)
    e3 = check_known_answer(pk3, params)
    print(f"known-answer check on the modified packer: relative difference {e3[2]:.1e}")
    Pn, E, gmax = relax_all_atom(pk3, params, Pn=P0, latn=latn)
    print(f"E {E:.6f} kcal/mol (its reference {e3[1]:.6f}, lowered by {e3[1] - E:.6f}), max |dE/dP| {gmax:.1e}, "
          f"max displacement {np.abs(Pn - P0).max():.4f} A")
    stage(pk3, params, Pn, latn, h, n_lowest)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--polymorphs", default="beta", help="comma-separated subset of beta,alpha,gamma")
    ap.add_argument("--h", type=float, default=1e-3, help="finite-difference step, A")
    ap.add_argument("--preset", default=PRESET)
    args = ap.parse_args(argv)
    known = {"beta": BETA, "alpha": ALPHA, "gamma": GAMMA}
    polymorphs = [known[k.strip()] for k in args.polymorphs.split(",") if k.strip()]
    for poly in polymorphs:
        t0 = time.time()
        packer, params = reference(poly, args.preset)
        print(f"# {poly[2]}: reference built in {time.time() - t0:.0f} s ({args.preset} + induced dipoles, deformable path)")
        report(packer, params, poly[2], args.h)
    return 0


if __name__ == "__main__":
    sys.exit(main())
