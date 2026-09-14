"""Read-only C-F label audit and reproduction of the recorded finite-molecule fit.

This is an input-data diagnostic, not native source/SCF/protocol admission.
No parameters are fitted and no molecular geometries are retained in output.
"""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import platform

import numpy as np

from polyfind import fitting as F
from polyfind.forcefield import EV_TO_KCAL,read_frames


def statistics(values):
    values = np.asarray(values,dtype=float)
    if values.size == 0:
        return None
    if not np.isfinite(values).all():
        raise ValueError("nonfinite reference statistic")
    return dict(count=int(values.size),minimum=float(values.min()),maximum=float(values.max()),
                mean=float(values.mean()),rms=float(np.sqrt(np.mean(values*values))),
                positive_fraction=float(np.mean(values > 0.)))


def terminal_cf(frame):
    """Atom-labelled radial probes; no assumption that labels are equilibrium forces."""
    for fluorine,element in enumerate(frame.elements):
        neighbours = frame.adjacency[fluorine]
        if element != "F" or len(neighbours) != 1:
            continue
        carbon = neighbours[0]
        if frame.elements[carbon] != "C":
            continue
        vector = frame.coords[fluorine]-frame.coords[carbon]
        length = float(np.linalg.norm(vector))
        if length == 0.:
            raise ValueError("zero C-F separation")
        yield carbon,fluorine,length,vector/length


def audit(path):
    path = Path(path).resolve()
    source_sha = hashlib.sha256(path.read_bytes()).hexdigest()
    raw = read_frames(str(path))
    fitted = F.load_reference(str(path))
    admitted_systems = {frame.system for frame in fitted}
    # The defining fit admission drops whole systems, preserving source order.
    source_indices = [i for i,frame in enumerate(raw) if frame.system in admitted_systems]
    if len(source_indices) != len(fitted):
        raise ValueError("fit admission no longer preserves whole source systems")
    design = F.ValenceDesign(fitted)
    params = F.valence_unpack(F.VAL_FITTED_X)
    _,predicted,_ = design.evaluate(params)
    predictions = {source:predicted[design.atom_frame == i] for i,source in enumerate(source_indices)}
    rows,worst = {},None
    for index,frame in enumerate(raw):
        if frame.forces is None or frame.energy is None or not np.isfinite(frame.forces).all():
            raise ValueError("complete finite energy/force labels required")
        row = rows.setdefault(frame.system,dict(frames=0,atoms=0,lengths=[],reference_radial=[],
                                               fitted_radial=[],radial_error=[],force_norm=[]))
        row["frames"] += 1
        row["atoms"] += frame.n_atoms
        row["force_norm"].extend(np.linalg.norm(frame.forces,axis=1))
        for carbon,fluorine,length,direction in terminal_cf(frame):
            measured = float(frame.forces[fluorine]@direction)
            row["lengths"].append(length)
            row["reference_radial"].append(measured)
            if index not in predictions:
                continue
            predicted_force = float(predictions[index][fluorine]@direction)
            error = predicted_force-measured
            row["fitted_radial"].append(predicted_force)
            row["radial_error"].append(error)
            if frame.system == "pvdf" and (worst is None or abs(error) > worst["absolute_error"]):
                worst = dict(absolute_error=abs(error),source_frame=index,carbon_atom=carbon,
                             fluorine_atom=fluorine,length_A=length,reference_force=measured,
                             fitted_force=predicted_force,direction=direction)
    summary = {}
    for name,row in sorted(rows.items()):
        summary[name] = dict(frames=row["frames"],atoms=row["atoms"],
                             included_in_recorded_fit=name in admitted_systems,
                             cf_length_A=statistics(row["lengths"]),
                             reference_outward_force_kcal_mol_A=statistics(row["reference_radial"]),
                             fitted_outward_force_kcal_mol_A=statistics(row["fitted_radial"]),
                             radial_error_kcal_mol_A=statistics(row["radial_error"]),
                             reference_atom_force_norm_kcal_mol_A=statistics(row["force_norm"]))
    if worst is not None:
        frame = raw[worst["source_frame"]]
        ff = F.FITTED_VALENCE.simple_ff()
        torsion = params["torsion"][F.torsion_types(frame)]
        controls = []
        for h in (1e-5,5e-6):
            energies = []
            for step in (-h,h):
                coords = frame.coords.copy()
                coords[worst["fluorine_atom"]] += step*worst["direction"]
                energies.append(ff.energy_frame(frame.with_coords(coords),torsion=torsion))
            force = -(energies[1]-energies[0])/(2*h)
            if abs(force-worst["fitted_force"]) > 2e-5:
                raise ValueError("recorded finite fit energy/force derivative differs")
            controls.append(dict(step_A=h,force_kcal_mol_A=force))
        worst.pop("direction")
        worst["finite_difference_controls"] = controls
    train,test = F.split_systems(fitted)
    if hashlib.sha256(path.read_bytes()).hexdigest() != source_sha:
        raise ValueError("reference bytes changed during audit")
    return dict(schema="finite_reference_cf_force_audit_v1",quantitatively_valid=False,
        scope="source-file label and finite fitted-model diagnostic only; no native SCF/force-sign/unit/boundary/source admission or bulk calibration",
        force_convention="positive terminal-F force points from C toward F; labels interpreted as eV/A by the existing reader, converted to kcal/(mol A)",
        input_eV_to_kcal_mol=EV_TO_KCAL,
        sampling="statistics over stored frames and bonds; correlated scan samples, not independent uncertainty estimates",
        source_path=str(path),source_sha256=source_sha,
        writer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        runtime=dict(python=platform.python_version(),numpy=np.__version__),
        counts=dict(frames=len(raw),atoms=sum(f.n_atoms for f in raw),systems=len(rows),
                    fit_frames=len(fitted),fit_systems=len(admitted_systems),
                    sources=dict(sorted(Counter(f.source for f in raw).items()))),
        excluded_systems=sorted(set(rows)-admitted_systems),systems=summary,
        recorded_fit_train=F.ValenceDesign([f for f in fitted if f.system in train]).errors(F.VAL_FITTED_X),
        recorded_fit_held=F.ValenceDesign([f for f in fitted if f.system in test]).errors(F.VAL_FITTED_X),
        worst_pvdf_radial_error=worst)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("path",type=Path)
    print(json.dumps(audit(parser.parse_args().path),sort_keys=True,separators=(",",":"),allow_nan=False))
