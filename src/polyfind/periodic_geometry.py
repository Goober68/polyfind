"""Periodic image geometry, independent of the crystal's Cartesian frame.

Legacy scalar repeats mean ``c * z_hat``. Explicit vectors carry the physical
translation of one chain repeat; they need not be lattice-normal or positive.
"""
from __future__ import annotations

import numpy as np
from dataclasses import dataclass
from itertools import islice


@dataclass(frozen=True, eq=False)
class ChainLayout:
    """Declared complete repeat atom IDs per chain, in chemical-local order."""

    atoms: np.ndarray
    reversed_of: np.ndarray
    elements: tuple[str, ...]

    def __post_init__(self):
        atoms = np.asarray(self.atoms)
        reversal = np.asarray(self.reversed_of)
        if atoms.ndim != 2 or not all(atoms.shape) or atoms.dtype.kind not in "iu":
            raise ValueError("chain atoms must be nonempty integer (chains,local_atoms)")
        if not np.array_equal(np.sort(atoms.ravel()),np.arange(atoms.size)):
            raise ValueError("chain layout must own every cell atom exactly once")
        if reversal.shape != (len(atoms),) or reversal.dtype.kind != "b":
            raise ValueError("one explicit boolean reversal is required per chain")
        elements = tuple(self.elements)
        if len(elements) != atoms.size or any(not isinstance(e,str) or not e for e in elements):
            raise ValueError("one element label is required per cell atom")
        for name,value in (("atoms",atoms.astype(np.int64)),("reversed_of",reversal)):
            value = np.array(value,copy=True)
            value.setflags(write=False)
            object.__setattr__(self,name,value)
        object.__setattr__(self,"elements",elements)

    @classmethod
    def canonical(cls,elements,n_chains,flip):
        if isinstance(n_chains,bool) or not isinstance(n_chains,(int,np.integer)) or n_chains not in (1,2):
            raise ValueError("canonical placement has one or two chains")
        if not np.isscalar(flip) or not np.isfinite(flip) or flip not in (0.,1.):
            raise ValueError("chain reversal must be exactly 0 or 1")
        labels = tuple(elements)*int(n_chains)
        atoms = np.arange(len(labels)).reshape(int(n_chains),len(elements))
        reversal = np.zeros(int(n_chains),dtype=bool)
        if n_chains == 2:
            reversal[1] = flip == 1.
        return cls(atoms,reversal,labels)

    @classmethod
    def from_cell(cls,cell):
        chains,local = np.asarray(cell.chain_of),np.asarray(cell.local_of)
        n = cell.n_per_chain
        N = len(cell.elements)
        if isinstance(n,bool) or not isinstance(n,(int,np.integer)) or n <= 0 or N == 0 or N % n:
            raise ValueError("cell must contain complete positive-size chain repeats")
        nc = N//n
        if chains.shape != (N,) or local.shape != (N,) or chains.dtype.kind not in "iu" or local.dtype.kind not in "iu":
            raise ValueError("cell chain/local maps must be integer (N,)")
        if np.any(chains < 0) or np.any(chains >= nc) or np.any(local < 0) or np.any(local >= n):
            raise ValueError("cell chain/local indices are outside the declared layout")
        keys = chains.astype(np.int64)*n+local.astype(np.int64)
        if not np.array_equal(np.sort(keys),np.arange(N)):
            raise ValueError("each chain/local pair must occur exactly once")
        atoms = np.empty(N,dtype=np.int64)
        atoms[keys] = np.arange(N)
        return cls(atoms.reshape(nc,n),cell.reversed_of,tuple(cell.elements))


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

    def pair_image_chunks(self, cutoff, chunk_elems=60000):
        """Yield full pair separations and effective lattice integers (I,N,N,3).

        A nearest fractional representative bounds each component by 1/2.
        Reciprocal-row widths bound every image that can enter the cutoff,
        including triclinic cells and atoms outside the stored primary cell.
        Effective integers, not the representative grid alone, own dD/dH.
        """
        cutoff = np.asarray(cutoff)
        if cutoff.shape != () or cutoff.dtype.kind not in "iuf" or not np.isfinite(cutoff) or cutoff <= 0:
            raise ValueError("cutoff must be a positive finite scalar")
        if isinstance(chunk_elems, bool) or not isinstance(chunk_elems, (int, np.integer)) or chunk_elems <= 0:
            raise ValueError("chunk_elems must be a positive integer")
        inverse = np.linalg.inv(self.lattice)
        fractional = self.coords @ inverse
        differences = fractional[:, None]-fractional[None, :]
        nearest = np.round(differences)
        base = (differences-nearest) @ self.lattice
        bounds = np.ceil(cutoff*np.linalg.norm(inverse, axis=0)+.5)
        if not np.isfinite(bounds).all() or np.any(bounds >= np.iinfo(np.int64).max):
            raise ValueError("cutoff image bounds are not finite representable integers")
        bounds = bounds.astype(np.int64)
        axes = [range(-int(b),int(b)+1) for b in bounds]
        grid = ((i,j,k) for i in axes[0] for j in axes[1] for k in axes[2])
        step = max(1, int(chunk_elems)//len(self.coords)**2)
        while True:
            integers = np.asarray(list(islice(grid,step)))
            if not len(integers):
                break
            yield base[None]+(integers@self.lattice)[:,None,None,:], integers[:,None,None,:]-nearest[None]


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

    ``c`` is a scalar or (M,) axial repeats (None means zero); ``repeat`` exclusively supplies
    (3,) or (M,3) translations instead. A three-row scalar batch is never
    mistaken for a Cartesian vector.
    """
    if repeat is not None and c is not None:
        raise ValueError("specify axial c or Cartesian repeat, not both")
    value = np.asarray((0.0 if c is None else c) if repeat is None else repeat)
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
