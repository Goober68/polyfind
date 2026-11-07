"""Pendant-group tilts: the relaxation the fit path's rigid CH2 and CF2 groups are denied.

The chain builder (:func:`polyfind.chain.pendant_frames`) puts the two pendants of a backbone
atom at ``cos(a/2) u +- sin(a/2) w`` -- ``u`` the outward bisector of the backbone angle, ``w``
the normal of the backbone plane -- and the crystal refinement moves cell, setting, torsions and
backbone angles but never a pendant relative to its carbon.  A DFT relaxation does: in a chain
whose two sides of a backbone atom are inequivalent (alpha's TGTG', not beta's planar zigzag) the
group turns a little as a rigid body.  Against the local frame ``(u, v, w)`` -- ``v`` along the
chain, from one backbone neighbour toward the other -- the turn has three components:

* **rock**, about ``w``: the group's bisector tilts toward a backbone neighbour;
* **wag**, about ``v``: the bisector tilts out of the backbone plane;
* **twist**, about ``u``: the X-X line turns toward the chain direction.

Dipole's PBE-D3 alpha start geometry has CF2 rock 2.3, wag 2.9, twist 1.3 deg (CH2 1.0, 0.7,
0.9) and beta none (``examples/dft_cell_geometry.py``); the fit path's cells have none anywhere.

:func:`tilt_relaxation` gives the energy the fitted potential releases when every group of a
refined reference may turn rigidly about its carbon, with the carbons, the cell, the bond lengths
and the X-C-X angles fixed, and with the valence bend terms pinned at the reference geometry so
the C-C-X bends resist with the valence fit's stiffness (``pvdf-dft-valence``: C-C-F 78.6,
C-C-H 60.4 kcal/mol/rad^2 in the kernel's 1/2 k convention).  It is exactly what a refinement
carrying the tilts as variables would recover, to first order in the coupled re-relaxation of
cell and backbone, and it is cheap: a dozen L-BFGS iterations on 3 angles per group with the
kernel's analytic forces.  docs/CRYSTAL_STAGE.md 5.22 has the numbers; under
``pvdf-dft-valence-belt3`` alpha releases 0.226 kcal/mol per monomer and beta nothing, which
moves the alpha-beta gap the fit matched with rigid pendants by 4.5 tolerances.

:func:`pinned_valence_packer` and :func:`belt_site_constraint` are shared with the Gamma-route
dielectric of ``examples/gamma_dielectric.py``, which first built them.
"""
from __future__ import annotations

import copy
from dataclasses import dataclass, replace

import numpy as np
from scipy.optimize import minimize

from . import phonon as PH
from .forcefield import SimpleFF
from .pack import chain_valence
from .refine import resolve_angle_stiffness

VALENCE_PRESET = "pvdf-dft-valence"
# a bond, by element pair (carbon first), is a distance under this (A); the reader of
# examples/dft_cell_geometry.py uses the same thresholds
BOND_MAX = {("C", "C"): 1.75, ("C", "H"): 1.25, ("C", "F"): 1.55, ("C", "Cl"): 1.95}


def pinned_valence_packer(ref, polymer, k_backbone: float | None = None, k_bond_min: float | None = None):
    """A copy of the reference's packer carrying valence terms stationary at its geometry.

    Bond lengths and angles of the valence preset's tables, with every equilibrium value pinned at
    the reference's own (the backbone C-C-C angles at the polymer's built angle, with
    ``k_backbone`` -- the refinement's restraint, so the all-atom energy agrees with the fit path's
    -- else the preset's stiffness), so the reference is a stationary point of the valence part.
    ``k_bond_min`` floors every stretch constant (kcal/mol/A^2, the kernel's 1/2 k convention): the
    valence fit left the C-C constant at its bound, 100, five times softer than a real C-C bond,
    and a floor shows what that softness does to an all-atom result (5.22).
    Returns ``(packer, n_backbone_angles, mean_backbone_angle_deg)``."""
    from . import fitting  # noqa: F401  (registers the valence preset)
    chain = ref.packer.chain
    val = SimpleFF.from_preset(VALENCE_PRESET)
    cv = chain_valence(chain, val.bond_table(), val.angle_table())
    X = np.asarray(chain.coords, dtype=float)
    cz = float(chain.c)
    r = np.linalg.norm(cv._bond_vectors(X[None], np.array([cz]))[0], axis=-1)
    ez = np.array([0.0, 0.0, 1.0])
    ri = X[cv.ang_i] + cv.ang_si[:, None] * cz * ez
    rk = X[cv.ang_k] + cv.ang_sk[:, None] * cz * ez
    rj = X[cv.ang_j]
    u, v = ri - rj, rk - rj
    theta = np.arccos(np.clip((u * v).sum(1) / np.linalg.norm(u, axis=1) / np.linalg.norm(v, axis=1), -1.0, 1.0))
    el = np.asarray(list(chain.elements))
    backbone = (el[cv.ang_i] == "C") & (el[cv.ang_j] == "C") & (el[cv.ang_k] == "C")
    built = np.array([polymer.backbone[m].backbone_angle for m in range(polymer.bonds_per_repeat)], dtype=float)
    if np.ptp(built) > 1e-9:
        raise ValueError(f"the polymer's backbone angles differ along the repeat: {built}")
    if k_backbone is None:
        k_backbone = resolve_angle_stiffness(None, None)
    t0 = np.where(backbone, np.deg2rad(float(built[0])), theta)
    kk = np.where(backbone, float(k_backbone), cv.ang_kk)
    bk = np.asarray(cv.bond_k, dtype=float)
    if k_bond_min is not None:
        bk = np.maximum(bk, float(k_bond_min))
    pk = copy.copy(ref.packer)
    pk._valence = replace(cv, bond_r0=r.copy(), bond_k=bk, ang_t0=t0, ang_kk=kk)
    pk.update_chain(chain)
    return pk, int(backbone.sum()), float(np.rad2deg(theta[backbone].mean())) if backbone.any() else None


def belt_site_constraint(pk, Pn):
    """``(J (3N, 3Na), atoms)``: the real atoms are free coordinates and each belt site is the fixed
    affine combination of its carbon and that carbon's two pendant atoms that places it (it lies in
    the pendant plane: the two fluorines of a CF2, the fluorine and hydrogen of a CHF), so a rigid
    motion of the group carries it exactly and a deformation to first order."""
    el = np.asarray(list(pk.elements) * pk.n_chains)
    atoms = np.flatnonzero(el != "X")
    col = {int(i): a for a, i in enumerate(atoms)}
    N, Na = pk.N, len(atoms)
    J = np.zeros((3 * N, 3 * Na))
    for i in atoms:
        J[3 * i:3 * i + 3, 3 * col[int(i)]:3 * col[int(i)] + 3] = np.eye(3)
    cz = float(pk.chain.c)
    shifts = [np.array([0.0, 0.0, k * cz]) for k in (-1, 0, 1)]

    def nearest(x, cands):
        best = None
        for i in cands:
            for sh in shifts:
                d = float(np.linalg.norm(Pn[x] - Pn[i] - sh))
                if best is None or d < best[0]:
                    best = (d, i, sh)
        return best

    for x in np.flatnonzero(el == "X"):
        s = x // pk.n
        chain_atoms = [i for i in atoms if i // pk.n == s]
        _, c, shc = nearest(x, [i for i in chain_atoms if el[i] == "C"])
        pc = Pn[c] + shc
        fl = sorted(((float(np.linalg.norm(pc - Pn[i] - sh)), int(i), sh) for i in chain_atoms if el[i] != "C" for sh in shifts),
                    key=lambda t: (t[0], t[1]))[:2]
        if fl[1][0] > 1.6:
            raise ValueError(f"belt site {x}: its carbon has no two pendant atoms within 1.6 A ({fl})")
        (_, f1, sh1), (_, f2, sh2) = fl
        A = np.stack([Pn[f1] + sh1 - pc, Pn[f2] + sh2 - pc], axis=1)  # (3, 2)
        ab, *_ = np.linalg.lstsq(A, Pn[x] - pc, rcond=None)
        if np.linalg.norm(A @ ab - (Pn[x] - pc)) > 1e-6:
            raise ValueError(f"belt site {x} is off the CF2 plane")
        for d in range(3):
            J[3 * x + d, 3 * col[int(c)] + d] = 1.0 - ab.sum()
            J[3 * x + d, 3 * col[f1] + d] = ab[0]
            J[3 * x + d, 3 * col[f2] + d] = ab[1]
    return J, atoms


def rotation(axis, angle_rad: float) -> np.ndarray:
    """Rodrigues: the proper rotation by ``angle_rad`` about ``axis``."""
    k = np.asarray(axis, dtype=float)
    k = k / np.linalg.norm(k)
    K = np.array([[0.0, -k[2], k[1]], [k[2], 0.0, -k[0]], [-k[1], k[0], 0.0]])
    return np.eye(3) + np.sin(angle_rad) * K + (1.0 - np.cos(angle_rad)) * (K @ K)


def _bonds(lattice, P, sym):
    """Bonded pairs ``(i, j, image)`` by distance over the 27 neighbouring images, each once."""
    lattice = np.asarray(lattice, dtype=float)
    images = np.array([[a, b, c] for a in (-1, 0, 1) for b in (-1, 0, 1) for c in (-1, 0, 1)])
    out = []
    for i in range(len(sym)):
        for j in range(i, len(sym)):
            key = tuple(sorted((sym[i], sym[j]), key=lambda e: (e != "C", e)))
            if key not in BOND_MAX:
                continue
            for img in images:
                if i == j and (not img.any() or tuple(img) < tuple(-img)):
                    continue
                if np.linalg.norm(P[j] + img @ lattice - P[i]) < BOND_MAX[key]:
                    out.append((i, j, img.copy()))
    return out


@dataclass
class PendantGroup:
    carbon: int          # index into the real atoms
    x: list              # the two pendant atoms' indices into the real atoms
    substituents: str    # "FF", "HH", ...
    frame: np.ndarray    # (3, 3) rows u (bisector), v (chain), w (normal)
    d: np.ndarray        # (2, 3) the pendant vectors from the carbon at the reference


def pendant_groups(lattice, P, sym) -> list[PendantGroup]:
    """Every backbone atom with two backbone neighbours and two single-atom pendants, with its frame."""
    nb = {i: [] for i in range(len(sym))}
    for i, j, img in _bonds(lattice, P, sym):
        nb[i].append((j, img))
        nb[j].append((i, -img))
    lattice = np.asarray(lattice, dtype=float)
    groups = []
    for c in range(len(sym)):
        if sym[c] != "C":
            continue
        cc = [(j, img) for j, img in nb[c] if sym[j] == "C"]
        xx = [(j, img) for j, img in nb[c] if sym[j] != "C"]
        if len(cc) != 2 or len(xx) != 2:
            continue
        n1 = P[cc[0][0]] + cc[0][1] @ lattice - P[c]
        n2 = P[cc[1][0]] + cc[1][1] @ lattice - P[c]
        n1, n2 = n1 / np.linalg.norm(n1), n2 / np.linalg.norm(n2)
        u = -(n1 + n2)
        u /= np.linalg.norm(u)
        v = n2 - n1
        v /= np.linalg.norm(v)
        w = np.cross(u, v)
        d = np.array([P[j] + img @ lattice - P[c] for j, img in xx])
        groups.append(PendantGroup(c, [j for j, _ in xx], "".join(sorted(sym[j] for j, _ in xx)), np.array([u, v, w]), d))
    return groups


def group_displacements(groups, tilts_rad, n_atoms: int) -> np.ndarray:
    """Displacements (3 n_atoms,) of the pendant atoms when each group turns rigidly by its
    ``(rock, wag, twist)`` about its carbon: rock about ``w``, wag about ``v``, twist about ``u``."""
    D = np.zeros((n_atoms, 3))
    for g, (rk, wg, tw) in zip(groups, np.asarray(tilts_rad, dtype=float).reshape(len(groups), 3)):
        u, v, w = g.frame
        R = rotation(u, tw) @ rotation(v, wg) @ rotation(w, rk)
        for j, d in zip(g.x, g.d):
            D[j] = R @ d - d
    return D.ravel()


@dataclass
class TiltRelaxation:
    dE_per_cell: float            # kcal/mol per cell, <= 0
    dE_per_monomer: float
    n_monomers: int
    tilts_deg: np.ndarray         # (n_groups, 3): rock, wag, twist, signed in each group's own frame
    groups: list
    P0: np.ndarray                # placed coordinates before (all sites, belt sites included)
    P1: np.ndarray                # and after
    lattice: np.ndarray
    elements: np.ndarray
    n_iterations: int
    message: str

    @property
    def tilts_by_substituents(self) -> dict:
        """Mean magnitude of rock, wag and twist (deg) over the groups of each substituent pair."""
        out = {}
        for g, t in zip(self.groups, self.tilts_deg):
            out.setdefault(g.substituents, []).append(np.abs(t))
        return {s: {"rock": float(np.mean([t[0] for t in ts])), "wag": float(np.mean([t[1] for t in ts])),
                    "twist": float(np.mean([t[2] for t in ts])), "n": len(ts)} for s, ts in sorted(out.items())}


def tilt_relaxation(ref, polymer, k_backbone: float | None = None, gtol: float = 1e-9) -> TiltRelaxation:
    """Let every CH2/CF2 group of the refined reference turn rigidly about its carbon.

    Carbons, cell, bond lengths and X-C-X angles stay; the belt sites follow their CF2 groups
    exactly (:func:`belt_site_constraint`); the valence bends are pinned at the reference
    (:func:`pinned_valence_packer`), so the C-C-X bends resist with the valence stiffness and the
    reference is a stationary point of the valence part.  Call under the parameter set's
    ``applied()`` so the packer is the fitted one."""
    pk, _, _ = pinned_valence_packer(ref, polymer, k_backbone)
    params = np.asarray(ref.params, dtype=float)
    P0, latn = PH.placed_coordinates(pk, params)
    P0, latn = np.asarray(P0, dtype=float), np.asarray(latn, dtype=float)
    J, atoms = belt_site_constraint(pk, P0)
    el = np.asarray(list(pk.elements) * pk.n_chains)
    sym = [str(e) for e in el[atoms]]
    groups = pendant_groups(latn, P0[atoms], sym)
    if not groups:
        raise ValueError("no backbone atom with two backbone neighbours and two single-atom pendants")
    n_mon = pk.n_chains * pk.chain.n_monomers
    E0, _ = PH.cell_energy_and_grad(pk, params, P0, latn)
    E0 = float(E0)
    h = 1e-6

    def eg(t):
        dX = group_displacements(groups, t, len(atoms))
        P = P0 + (J @ dX).reshape(-1, 3)
        E, g = PH.cell_energy_and_grad(pk, params, P, latn)
        gX = J.T @ np.asarray(g, dtype=float).ravel()
        grad = np.empty_like(t)
        for i in range(len(t)):
            tp, tm = t.copy(), t.copy()
            tp[i] += h
            tm[i] -= h
            grad[i] = gX @ (group_displacements(groups, tp, len(atoms)) - group_displacements(groups, tm, len(atoms))) / (2 * h)
        return float(E), grad

    res = minimize(eg, np.zeros(3 * len(groups)), jac=True, method="L-BFGS-B",
                   options={"maxiter": 2000, "gtol": gtol, "ftol": 1e-15})
    E1 = float(res.fun)
    if E1 > E0 + 1e-9:  # the start is admissible; never report a rise
        E1, x = E0, np.zeros_like(res.x)
    else:
        x = res.x
    P1 = P0 + (J @ group_displacements(groups, x, len(atoms))).reshape(-1, 3)
    return TiltRelaxation(E1 - E0, (E1 - E0) / n_mon, int(n_mon), np.degrees(x.reshape(len(groups), 3)), groups,
                          P0, P1, latn, el, int(res.nit), str(res.message))
