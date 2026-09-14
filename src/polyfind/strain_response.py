"""Cartesian frame changes and proper reduction of complete affine responses.

Columns are derivatives per engineering Voigt strain (xx,yy,zz,yz,xz,xy).
This is tensor arithmetic, not a material-validity or numerical pass policy.
"""
import numpy as np

from .mechanics import strain_tensor


def _real(value, shape):
    array = np.asarray(value)
    if array.shape != shape or array.dtype.kind not in "iuf" or not np.isfinite(array).all():
        raise ValueError(f"expected finite real array with shape {shape}")
    return np.array(array, dtype=float, copy=True)


def frame_weights(frame) -> np.ndarray:
    """(6,6): new-frame unit strains expressed in old engineering components."""
    L = _real(frame, (3, 3))
    if not np.allclose(L @ L.T, np.eye(3), rtol=0., atol=1e-12):
        raise ValueError("frame must be orthogonal; it is not repaired")
    result = []
    for component in np.eye(6):
        E = L.T @ strain_tensor(component) @ L
        result.append([E[0, 0], E[1, 1], E[2, 2], 2*E[1, 2], 2*E[0, 2], 2*E[0, 1]])
    return np.asarray(result)


def rotate_columns(columns, frame) -> np.ndarray:
    """Complete (3,6) response -> arbitrary orthogonal Cartesian frame."""
    response = _real(columns, (3, 6))
    weights = frame_weights(frame)
    return np.asarray(frame, dtype=float) @ response @ weights.T


def proper_columns(dipole_per_reference_volume, polarization) -> np.ndarray:
    """Vanderbilt proper response from (1/V0)dmu/de and SAME total P.

    Subtract eps_J @ P for each symmetric unit engineering strain. No column
    is omitted or zero-filled, including chain-axis shears.
    """
    response = _real(dipole_per_reference_volume, (3, 6))
    P = _real(polarization, (3,))
    correction = np.stack([strain_tensor(component) @ P for component in np.eye(6)], axis=1)
    return response - correction
