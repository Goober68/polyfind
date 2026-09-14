"""Compare old and new dipole evaluators on exactly the same decoded nuclei.

The legacy checkout is read-only. Both workers use their own existing
clamped-ion dipole implementation; this control does not copy either model.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np


def worker(repository, record_dir):
    repository, record_dir = Path(repository).resolve(), Path(record_dir).resolve()
    sys.path.insert(0, str(repository / "src"))
    sys.path.insert(0, str(repository / "examples"))
    import clamped_ion_columns as model

    record = json.loads((record_dir / "beta_pvdf_clamped_ion_columns.json").read_text())
    xyz_bytes = (record_dir / record["reference_state"]["geometry_file"]).read_bytes()
    if hashlib.sha256(xyz_bytes).hexdigest() != record["reference_state"]["geometry_xyz_sha256"]:
        raise ValueError("reference geometry bytes changed")
    lines = xyz_bytes.decode("ascii").splitlines()
    count = int(lines[0])
    if len(lines) != count + 2:
        raise ValueError("reference geometry atom count changed")
    Xp = np.array([[float(v) for v in line.split()[1:]] for line in lines[2:]])
    state = record["reference_state"]
    L = np.asarray(state["frame_map_L_rows_producer_from_packer"], dtype=float)
    np.testing.assert_allclose(L @ L.T, np.eye(3), atol=1e-14)
    # Row coordinates in provider axes are X @ L.T.
    X = Xp @ L
    H = np.asarray(state["lattice_rows_producer_frame_A"], dtype=float) @ L
    params = np.asarray(state["packer_params_a_b_gamma_phi1_phi2_dz_flip"], dtype=float)
    c = float(state["c_A"])
    packer, _, _, _ = model.ISJ.beta_reference(model.ISJ.SimpleFF.from_preset(record["presets"]["charge_flux"]),
                                              model.ISJ.Polarizable())
    if X.shape != (packer.N, 3):
        raise ValueError("reference atoms do not match the model")
    points = {}
    # Precisely the common domain; no legacy chain-axis shear is synthesized.
    for label, K, amplitude in [("zero", None, 0.)] + [
        (f"{K}:{h}:{sign}", K, sign * h)
        for K in (0, 1, 2, 5) for h in (.0025, .005) for sign in (-1, 1)
    ]:
        eps = np.zeros(6)
        if K is not None:
            eps[K] = amplitude
        F = np.eye(3) + model.M.strain_tensor(eps)
        sc = model.M.strained_cell(params, c, eps)
        points[label] = np.asarray(model.cell_dipole(packer, sc.params, X @ F, H @ F)).tolist()
    owners = model.SOURCES
    print(json.dumps({"points_charge_and_induced_eA": points,
                      "input_geometry_sha256": hashlib.sha256(xyz_bytes).hexdigest(),
                      "source_sha256": {p: model.sha256_text(repository / p) for p in owners}}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--legacy-repository", required=True)
    parser.add_argument("--record-dir", required=True)
    parser.add_argument("--worker-repository")
    parser.add_argument("--out")
    args = parser.parse_args()
    if args.worker_repository:
        worker(args.worker_repository, args.record_dir)
        return
    current = Path(__file__).resolve().parents[1]
    env = dict(os.environ, POLYFIND_DEVICE="cpu", OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    records = []
    for repository in (args.legacy_repository, current):
        out = subprocess.check_output([sys.executable, str(Path(__file__).resolve()),
                                       "--legacy-repository", str(args.legacy_repository),
                                       "--record-dir", args.record_dir,
                                       "--worker-repository", str(repository)], env=env, text=True)
        records.append(json.loads(out))
    old, new = records
    if old["input_geometry_sha256"] != new["input_geometry_sha256"] or old["points_charge_and_induced_eA"].keys() != new["points_charge_and_induced_eA"].keys():
        raise ValueError("baseline control inputs differ")
    differences = {key: (np.array(new["points_charge_and_induced_eA"][key]) - old["points_charge_and_induced_eA"][key]).tolist()
                   for key in old["points_charge_and_induced_eA"]}
    sys.path.insert(0, str(current / "src"))
    sys.path.insert(0, str(current / "examples"))
    from clamped_ion_columns import sha256_text

    report = {"scope": "same 10-decimal XYZ decoded nuclei; 17 complete dipole observations in the shared four-column domain; not material/DFT qualification",
              "control_source_sha256": sha256_text(Path(__file__)),
              "source_hash_canonicalization": "UTF-8/ASCII text with CRLF and CR normalized to LF",
              "numpy_version": np.__version__, "legacy": old, "vector_repeat": new,
              "differences_eA": differences,
              "max_abs_component_difference_eA": max(float(np.max(np.abs(v))) for v in differences.values()),
              "quantitatively_valid": False}
    if args.out:
        output = Path(args.out)
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(report, indent=1) + "\n", encoding="utf-8", newline="\n")
        print(json.dumps({"output": str(output.resolve()),
                          "points": len(differences),
                          "max_abs_component_difference_eA": report["max_abs_component_difference_eA"]}))
    else:
        print(json.dumps(report, indent=1))


if __name__ == "__main__":
    main()
