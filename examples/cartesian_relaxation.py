"""CPU-only exploratory complete-model stationary cells, never dynamics.

Produces a self-validating compact JSON unit on stdout. Failed/unconverged
work keeps protocol/source identity and diagnosis, not rejected geometry.
Accepted geometry remains a model stationary point, NOT stable/calibrated
material evidence. A seed unit supplies actual coordinates and declared maps.
"""
import argparse
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import platform

import numpy as np
import scipy

from polyfind import pack as pack_mod
from polyfind.cartesian_mechanics import CartesianChart,CartesianRelaxer,CartesianTolerance,RelaxationPhase
from polyfind.ewald import EwaldSpec
from polyfind.fitting import FITTED_VALENCE
from polyfind.forcefield import SimpleFF
from polyfind.polarizability import Polarizable
from polyfind.polymers import PVDF,THREE_STATE
from polyfind.topology import Cell,build_cell,check_topology


class PilotRecord:
    """Own canonical JSON bytes, version and complete integrity coverage."""
    SCHEMA = "complete_cartesian_stationary_model_pilot_v1"

    @staticmethod
    def canonical(value):
        return json.dumps(value,sort_keys=True,separators=(",",":"),allow_nan=False).encode()

    @classmethod
    def encode(cls,payload):
        if payload.get("schema") != cls.SCHEMA:
            raise ValueError("model pilot schema differs")
        return cls.canonical(dict(payload=payload,payload_sha256=hashlib.sha256(cls.canonical(payload)).hexdigest()))

    @classmethod
    def decode(cls,raw):
        envelope = json.loads(raw)
        if set(envelope) != {"payload","payload_sha256"}:
            raise ValueError("whole model pilot envelope required")
        payload = envelope["payload"]
        if cls.encode(payload) != raw:
            raise ValueError("whole model pilot integrity/schema differs")
        return payload


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--field",nargs=3,type=float,default=(0.,0.,0.),metavar=("EX","EY","EZ"))
    parser.add_argument("--prestrain",type=float,default=0.,help="actual fractional stretch along the seed c row")
    parser.add_argument("--fixed-cell",action="store_true")
    parser.add_argument("--maxiter",type=int,default=500)
    parser.add_argument("--seed",type=Path)
    args = parser.parse_args()
    if not np.isfinite([*args.field,args.prestrain]).all() or args.prestrain <= -1.:
        raise ValueError("finite field and positive chain stretch required")
    with FITTED_VALENCE.applied():
        chain = pack_mod.periodic_chain(PVDF,(0,0,0,0),THREE_STATE)
        pk = pack_mod.CrystalPacker(chain,coulomb="ewald",ewald=EwaldSpec(),
                                    valence=SimpleFF.from_preset("pvdf-dft-valence"),
                                    charge_flux=SimpleFF.from_preset("pvdf-dft-valence-flux-born"),
                                    polarizable=Polarizable(),lj_cutoff="force",field=args.field)
    if pk.xp is not np:
        raise RuntimeError("GPU reserved: this pilot requires the CPU backend")
    package = Path(pack_mod.__file__).resolve().parent
    sources = {path.name:hashlib.sha256(path.read_bytes()).hexdigest() for path in sorted(package.glob("*.py"))}
    cell = build_cell(chain,[4.6,8.6,90.,0.,0.,.4,0.],packer=pk)
    initial_kind = "two-chain doubled-chemical-period trans seed with independent 0.003 A perturbations, RNG7"
    if args.seed is not None:
        seed = PilotRecord.decode(args.seed.read_bytes())
        if seed["sources_sha256"] != sources or seed["admission"] != "stationary_model_geometry_only":
            raise ValueError("accepted same-source complete model seed required")
        state = seed["geometry"]
        cell = Cell(state["elements"],np.asarray(state["coords_A"]),np.asarray(state["lattice_rows_A"]),
                    np.asarray(state["chain_of"]),np.asarray(state["local_of"]),state["n_per_chain"],
                    np.asarray(state["reversed_of"],dtype=bool))
        initial_kind = "accepted actual seed unit "+hashlib.sha256(args.seed.read_bytes()).hexdigest()
    else:
        cell = replace(cell,coords=cell.coords+np.random.default_rng(7).normal(size=cell.coords.shape)*.003)
    axis = cell.lattice[2]/np.linalg.norm(cell.lattice[2])
    stretch = np.eye(3)+args.prestrain*np.outer(axis,axis)
    cell = replace(cell,coords=cell.coords@stretch,lattice=cell.lattice@stretch)
    initial = dict(elements=cell.elements,coords_A=cell.coords.tolist(),lattice_rows_A=cell.lattice.tolist(),
                   chain_of=cell.chain_of.tolist(),local_of=cell.local_of.tolist(),
                   n_per_chain=cell.n_per_chain,reversed_of=cell.reversed_of.tolist())
    free = np.full(6,not args.fixed_cell,dtype=bool)
    chart = CartesianChart(cell,free)
    tolerance = CartesianTolerance()
    payload = dict(schema=PilotRecord.SCHEMA,quantitatively_valid=False,
                   scope="stationary complete model only; not stability/BZ/size/native calibration/material/dynamics/barrier/loss",
                   sources_sha256=sources,writer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                   runtime=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,device="cpu"),
                   model=dict(polymer="pvdf",chemical_periods=2,independent_chains=2,
                              valence="pvdf-dft-valence",flux="pvdf-dft-valence-flux-born",polarizable=True,
                              coulomb="ewald",ewald_accuracy=pk._ewald.spec.accuracy,lj_cutoff="force"),
                   conditions=dict(field_V_A=list(args.field),initial_chain_axis=axis.tolist(),
                                   applied_seed_chain_stretch=args.prestrain,free_log_strain=free.tolist(),
                                   macro_polar_rotation="clamped relative to seed",voigt_order=["xx","yy","zz","yz","xz","xy"],
                                   atomic_force_tolerance_kcal_mol_A=tolerance.atomic_force_kcal_mol_A,
                                   stress_tolerance_MPa=tolerance.stress_MPa,maxiter=args.maxiter),
                   initial_kind=initial_kind,initial_geometry_sha256=hashlib.sha256(PilotRecord.canonical(initial)).hexdigest())
    driver = CartesianRelaxer(pk,chart,tolerance)
    try:
        result = driver.run(maxiter=args.maxiter)
        checkpoint = result.checkpoint
        report = checkpoint.stationarity
        payload.update(phase=checkpoint.phase.value,admission="not_admitted",accepted_steps=checkpoint.accepted_steps,
                       energy_kcal_mol=checkpoint.energy_kcal_mol,evaluations=result.evaluations,
                       optimizer_status=result.optimizer_status,message=result.message,
                       max_atomic_force_kcal_mol_A=report.max_atomic_force_kcal_mol_A,
                       stress_MPa=report.stress_MPa.tolist(),log_conjugate_stress_MPa=report.log_conjugate_stress_MPa.tolist(),
                       max_selected_physical_stress_MPa=report.max_selected_physical_stress_MPa,
                       max_free_log_stress_MPa=report.max_free_log_stress_MPa)
        if checkpoint.phase == RelaxationPhase.CONVERGED:
            final = chart.cell(checkpoint.parameters)
            topology = check_topology(final,chain)
            payload["topology_ok"] = topology.ok
            payload["safe_covalent_scale_window"] = list(topology.safe_scale_window)
            if not topology.ok:
                payload["topology_diagnosis"] = dict(worst_intended_bond=topology.max_bond_pair,
                    max_bond_covalent_radius_ratio=topology.max_bond_ratio,
                    min_nonbond_covalent_radius_ratio=topology.min_nonbond_ratio,
                    missing_by_scale=topology.missing,extra_by_scale=topology.extra,
                    components_by_scale=topology.components,
                    remaining_decision="diagnose model geometry/calibration; do not change covalent-distance gate")
            if topology.ok:
                permanent,induced = pk.chain_dipole(final.coords,final.lattice,chart.layout)
                payload.update(admission="stationary_model_geometry_only",geometry=dict(
                    elements=final.elements,coords_A=final.coords.tolist(),lattice_rows_A=final.lattice.tolist(),
                    chain_of=final.chain_of.tolist(),local_of=final.local_of.tolist(),n_per_chain=final.n_per_chain,
                    reversed_of=final.reversed_of.tolist()),dipole_eA=dict(permanent=permanent.tolist(),induced=induced.tolist()))
    except Exception as error:
        payload.update(phase=driver.phase.value,admission="not_admitted",diagnosis=f"{type(error).__name__}: {error}")
    print(PilotRecord.encode(payload).decode(),end="")


if __name__ == "__main__":
    main()
