"""One IUPAC dihedral/Fourier implementation for fitting and molecular kernels."""
from __future__ import annotations

import numpy as np


def _geometry(points, xp):
    points = xp.asarray(points)
    if points.shape[-2:] != (4, 3) or points.dtype.kind not in "iuf" or not bool(xp.isfinite(points).all()):
        raise ValueError("dihedral points must be finite real (...,4,3)")
    if points.dtype.kind in "iu":
        points = points.astype(float)
    a, b, c, d = (points[..., k, :] for k in range(4))
    b0, b1, b2 = b-a, c-b, d-c
    n1, n2 = xp.cross(b0, b1), xp.cross(b1, b2)
    norm = xp.sqrt((b1*b1).sum(-1))
    s1, s2 = (n1*n1).sum(-1), (n2*n2).sum(-1)
    if bool(xp.any(norm == 0)) or bool(xp.any(s1 == 0)) or bool(xp.any(s2 == 0)):
        raise ValueError("dihedral is undefined at a zero bond or collinear arm")
    if not bool(xp.isfinite(norm).all()) or not bool(xp.isfinite(s1).all()) or not bool(xp.isfinite(s2).all()):
        raise ValueError("dihedral geometry is nonfinite")
    phi = xp.arctan2(norm*(b0*n2).sum(-1), (n1*n2).sum(-1))
    return phi, b0, b1, b2, n1, n2, norm, s1, s2


def dihedral_radians(points, xp=np):
    """(...,4,3) ordered points -> (...) radians, trans = +/-pi."""
    return _geometry(points, xp)[0]


def dihedral_value_gradient(points):
    """(...,4,3) -> angles in radians and all four point gradients."""
    points = np.asarray(points)
    if points.shape[-2:] != (4, 3) or points.dtype.kind not in "iuf" or not np.isfinite(points).all():
        raise ValueError("dihedral points must be finite real (...,4,3)")
    phi, b0, b1, b2, n1, n2, norm, s1, s2 = _geometry(points, np)
    gi = -norm[..., None]*n1/s1[..., None]
    gl = norm[..., None]*n2/s2[..., None]
    f = ((b0*b1).sum(-1)/norm**2)[..., None]
    g = ((b2*b1).sum(-1)/norm**2)[..., None]
    gj = -(1+f)*gi + g*gl
    gk = f*gi - (1+g)*gl
    gradient = np.stack([gi, gj, gk, gl], axis=-2)
    if not np.isfinite(gradient).all():
        raise ValueError("dihedral gradient is nonfinite")
    return phi, gradient


def indexed_dihedral_gradient(coords, indices):
    """(N,3) coordinates and (K,4) indices -> (phi, dphi/dpoints)."""
    return dihedral_value_gradient(np.asarray(coords)[np.asarray(indices)])


def fourier_terms(phi, coefficients, xp=np):
    """Per-angle kcal/mol; coefficient triples broadcast against radian angles."""
    V = xp.asarray(coefficients, dtype=phi.dtype)
    V1, V2, V3 = V[..., 0], V[..., 1], V[..., 2]
    return .5*(V1*(1+xp.cos(phi)) + V2*(1-xp.cos(2*phi)) + V3*(1+xp.cos(3*phi)))


def fourier_derivative(phi, coefficients):
    """dE/dphi in kcal/(mol rad), same Fourier convention as fourier_terms."""
    V = np.asarray(coefficients, dtype=np.asarray(phi).dtype)
    return .5*(-V[..., 0]*np.sin(phi) + 2*V[..., 1]*np.sin(2*phi) - 3*V[..., 2]*np.sin(3*phi))
