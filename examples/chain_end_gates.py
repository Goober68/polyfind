"""What the force field says to the pore-mouth machine: the torsion gate and the registry ring.

    python examples/chain_end_gates.py                # fitted single-chain potential, four units
    python examples/chain_end_gates.py --illustrative # the built-in potential beside it

Prints, for VDF, TrFE, CFE and CTFE ends of an 8-mer, the trans/gauche energy of the entering
bond, the transverse dipole change and the field that biases the bond either way by 2 kT
(design 2.5 / D21: the deciding number is 2 kT per bond at no more than 2 V/nm); then the
registry corrugation of beta-PVDF as one chain slides a repeat (design 2.4 / D21 and the
backtracking force of D24).  ``docs/ELECTROSTRICTION_RESEARCH.md`` section 10 reads them.
"""
import argparse
import contextlib

import polyfind.fitting  # noqa: F401  (registers the fitted SimpleFF presets)
from polyfind import chainend as CE
from polyfind import mechanics as M
from polyfind.forcefield import SimpleFF
from polyfind.polymers import CFE, CTFE, PVDF, TRFE

T = 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--illustrative", action="store_true", help="also run the built-in potential")
    args = ap.parse_args()
    potentials = [("pvdf-dft-fit", SimpleFF.from_preset("pvdf-dft-fit"))]
    if args.illustrative:
        potentials.append(("illustrative", SimpleFF()))
    print("torsion gate at the end of an all-trans 8-mer (bond 15 enters the tunnel; bond 14 is one inside):")
    for label, ff in potentials:
        print(f"  [{label}]")
        for polymer in (PVDF, TRFE, CFE, CTFE):
            for bond in (15, 14):
                print("   " + CE.torsion_gate(polymer, ff, bond=bond).table())
    print("  frozen bond angles: the chlorine-bearing ends' E_G - E_T are strain artefacts; the dipole changes are not.")

    print("\nregistry corrugation of beta-PVDF (chain 2 slid along c, cell held):")
    from polyfind.fitting import FITTED_DFT
    for label, ctx in (("illustrative", contextlib.nullcontext()), ("pvdf-dft-fit", FITTED_DFT.applied())):
        with ctx:
            ref, _ = M.refined_reference(PVDF, [T, T], label="beta", lj_cutoff="force", pack_kw={"table_cache_dir": None})
            ref = M.relax_reference(ref)
            cor = CE.registry_corrugation(ref, monomers_per_cell=2)
            print(f"  [{label}] " + cor.table())


if __name__ == "__main__":
    main()
