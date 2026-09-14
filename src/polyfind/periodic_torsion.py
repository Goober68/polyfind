"""Topology-owned Cartesian Fourier torsions of one translated chain repeat."""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np

from .periodic_geometry import repeat_rows
from .torsion_geometry import dihedral_radians, dihedral_value_gradient, fourier_terms, fourier_derivative


@dataclass(frozen=True, eq=False)
class ChainTorsion:
    n_atoms: int
    atoms: np.ndarray
    images: np.ndarray
    coefficients: np.ndarray

    def __post_init__(self):
        if isinstance(self.n_atoms, bool) or not isinstance(self.n_atoms, (int, np.integer)) or self.n_atoms <= 0:
            raise ValueError("torsion topology needs a positive atom count")
        atoms, images, V = np.asarray(self.atoms), np.asarray(self.images), np.asarray(self.coefficients)
        if atoms.ndim != 2 or atoms.shape[1] != 4 or not len(atoms) or atoms.dtype.kind not in "iu":
            raise ValueError("torsion atoms must be nonempty integer (K,4)")
        if np.any(atoms < 0) or np.any(atoms >= self.n_atoms):
            raise ValueError("torsion atom index outside the repeat")
        if images.shape != atoms.shape or images.dtype.kind not in "iu":
            raise ValueError("torsion image indices must be integer (K,4)")
        if int(atoms.max()) > np.iinfo(np.int64).max or int(images.max()) > np.iinfo(np.int64).max or int(images.min()) < np.iinfo(np.int64).min:
            raise ValueError("torsion indices must be representable as int64")
        if V.shape not in ((3,), (len(atoms), 3)) or V.dtype.kind not in "iuf" or not np.isfinite(V).all():
            raise ValueError("torsion coefficients must be finite real (3,) or (K,3)")
        for name, value in (("atoms", atoms), ("images", images), ("coefficients", np.broadcast_to(V, (len(atoms), 3)))):
            copy = np.array(value, dtype=float if name == "coefficients" else np.int64, copy=True)
            copy.setflags(write=False)
            object.__setattr__(self, name, copy)

    def _points(self, coords, c, repeat):
        X = np.asarray(coords)
        if X.dtype.kind not in "iuf" or not np.isfinite(X).all():
            raise ValueError("repeat coordinates must be finite and real")
        if X.shape == (self.n_atoms, 3):
            X = X[None]
        if X.ndim != 3 or X.shape[1:] != (self.n_atoms, 3) or not len(X):
            raise ValueError("repeat coordinates must have shape (N,3) or (M,N,3)")
        X = np.asarray(X, dtype=float)
        repeats = repeat_rows(c, len(X), repeat=repeat)
        points = X[:, self.atoms] + self.images[None, ..., None]*repeats[:, None, None, :]
        return points

    def angles(self, coords, c=None, *, repeat=None):
        """(M,K) angles in radians, derived from complete periodic quadruples."""
        return dihedral_radians(self._points(coords, c, repeat))

    def energy(self, coords, c=None, *, repeat=None):
        """(M,) Fourier energy, kcal/mol per repeat, not a metadata constant."""
        return fourier_terms(self.angles(coords, c, repeat=repeat), self.coefficients).sum(axis=1)

    def energy_and_repeat_grad(self, coords, repeat):
        """(E, gcoords (N,3), grepeat (3,)) of ONE complete repeat geometry."""
        points = self._points(coords, None, repeat)
        if len(points) != 1:
            raise ValueError("repeat gradient requires one geometry")
        phi, dphi = dihedral_value_gradient(points[0])
        forces = fourier_derivative(phi, self.coefficients)[:, None, None]*dphi
        gradient = np.zeros((self.n_atoms, 3))
        np.add.at(gradient, self.atoms.ravel(), forces.reshape(-1, 3))
        repeat_gradient = (forces*self.images[..., None]).sum(axis=(0, 1))
        return float(fourier_terms(phi, self.coefficients).sum()), gradient, repeat_gradient

    def energy_and_grad(self, coords, c):
        """Scalar-axial repeat pullback of the complete Cartesian derivative."""
        E, gX, gr = self.energy_and_repeat_grad(coords, repeat_rows(c, 1)[0])
        return E, gX, float(gr[2])
