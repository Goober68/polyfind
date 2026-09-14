"""Periodic image geometry, independent of the crystal's Cartesian frame.

Legacy scalar repeats mean ``c * z_hat``. Explicit vectors carry the physical
translation of one chain repeat; they need not be lattice-normal or positive.
"""
from __future__ import annotations

import numpy as np
from dataclasses import dataclass


@dataclass(frozen=True, eq=False)
class PlacedCell:
    """Cartesian atoms and lattice rows; no canonical-axis restriction."""

    coords: np.ndarray
    lattice: np.ndarray

    def __post_init__(self):
        for name in ("coords", "lattice"):
            value = np.asarray(getattr(self, name))
            if value.dtype.kind not in "iuf" or not np.isfinite(value).all():
                raise ValueError("placed geometry must be finite and real")
            value = np.array(value, dtype=float, copy=True)
            value.setflags(write=False)
            object.__setattr__(self, name, value)
        if self.coords.ndim != 2 or self.coords.shape[1] != 3 or not len(self.coords):
            raise ValueError("placed coordinates must have shape (N,3), N > 0")
        if self.lattice.shape != (3, 3):
            raise ValueError("lattice rows must have shape (3,3)")
        volume = np.linalg.det(self.lattice)
        if not np.isfinite(volume) or volume <= 0.0:
            raise ValueError("lattice rows must define a nonsingular right-handed cell")

    def deformed(self, deformation) -> PlacedCell:
        """Apply the same row-vector affine map to every atom and lattice row."""
        F = np.asarray(deformation)
        if F.shape != (3, 3) or F.dtype.kind not in "iuf" or not np.isfinite(F).all():
            raise ValueError("deformation must be a finite real (3,3) matrix")
        return PlacedCell(self.coords @ F, self.lattice @ F)


def repeat_vector(c) -> np.ndarray:
    """One scalar axial repeat or one Cartesian translation -> (3,) in A."""
    value = np.asarray(c)
    if value.dtype.kind not in "iuf" or not np.isfinite(value).all():
        raise ValueError("repeat must be finite and real")
    if value.ndim == 0:
        return np.array([0.0, 0.0, float(value)])
    if value.shape != (3,):
        raise ValueError("repeat must be a scalar or a Cartesian (3,) vector")
    return np.array(value, dtype=float, copy=True)


def repeat_gradient(c, gradient):
    """Pull a Cartesian repeat derivative back to the caller's repeat coordinates."""
    repeat_vector(c)
    g = np.asarray(gradient, dtype=float).reshape(3)
    return float(g[2]) if np.asarray(c).ndim == 0 else g.copy()


def repeat_rows(c, rows: int, *, repeat=None) -> np.ndarray:
    """Batched translations, with unambiguous legacy scalar-row semantics.

    ``c`` is a scalar or (M,) axial repeats; ``repeat`` explicitly supplies
    (3,) or (M,3) translations instead. A three-row scalar batch is never
    mistaken for a Cartesian vector.
    """
    value = np.asarray(c if repeat is None else repeat)
    if value.dtype.kind not in "iuf" or not np.isfinite(value).all():
        raise ValueError("repeat must be finite and real")
    if repeat is not None:
        if value.shape not in ((3,), (rows, 3)):
            raise ValueError("explicit repeats must have shape (3,) or (M,3)")
        return np.broadcast_to(value, (rows, 3)).astype(float, copy=True)
    if value.ndim > 1 or value.size not in (1, rows):
        raise ValueError("axial repeats must be scalar or (M,)")
    out = np.zeros((rows, 3))
    out[:, 2] = np.broadcast_to(value.reshape(-1), (rows,))
    return out
