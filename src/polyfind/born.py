r"""Born effective charges of a packed cell, from the model's own dipole.

``Z*_{i,ab} = d mu_a / d u_{i,b}``: the change of the cell dipole ``mu = sum_j q_j r_j + sum_j
p_j`` when ONE atom ``i`` of the cell -- together with every one of its lattice images, which is
what a Gamma-point displacement is -- moves by ``u`` along ``b``, at zero macroscopic field (the
Ewald sum's tinfoil boundary), in units of ``e``.  It is the quantity a DFPT code prints as the
Born effective charge, and it is the like-for-like comparison for a charge model whose charges
follow the geometry (:class:`polyfind.forcefield.FluxTopology`) and whose atoms carry induced
dipoles (:mod:`polyfind.polarizability`): both respond to a displacement, and both are in here.

Three contributions add:

* the static charge, ``q_i I`` (an off-site charge, :func:`polyfind.forcefield.project_offset_charges`,
  makes it ``q_i (r + d) / r`` perpendicular to its bond and ``q_i`` along it);
* the charge flux, ``sum_j r_j (dq_j / du_i)``: a bond-stretch flux ``k_bond`` puts
  ``-k_bond r u u^T`` on the terminal atom (``u`` the unit bond vector: charge flowing *along*
  the bond, the response that is large along the bond and zero across it), while the angle
  flux puts ``u n^T`` terms with ``n`` perpendicular to the bond -- traceless, and antisymmetric
  in the off-diagonals;
* the induced dipoles, ``sum_j dp_j / du_i``: the electronic screening of the moving charge by
  every polarizable neighbour, resolved through the self-consistent solve of
  :meth:`polyfind.pack.CrystalPacker._polarize`.

The acoustic sum rule ``sum_i Z*_i = 0`` holds identically here -- a rigid translation moves no
charge and induces no dipole, and the cell is neutral -- so it is a check on the arithmetic and
never on the model.  The derivative is a central difference of the exact dipole; the dipole is
linear algebra in the displaced coordinates, so the error is ``O(h^2)`` and
``examples/fit_born_flux.py`` measures it by halving ``h``.

Frame: the packer's, ``x`` = polar, ``y`` = the long lateral axis, ``z`` = the chain axis.
:func:`canonical` maps to the DFPT provider's ``(a = long, b = polar, c = chain)`` and reflects
every pendant atom onto its representative (the one on the ``+a`` side of its carbon, with the
polar axis pointing from the carbon towards its heaviest pendant), so that symmetry-equivalent
atoms can be averaged and compared with one row of the reference table.  Only the sign of the
``ab``/``ba`` off-diagonals depends on that convention; the diagonals do not.
"""
from __future__ import annotations

import json
import os
from dataclasses import dataclass

import numpy as np

from . import backend as bk

REFERENCE_AXES = ("a", "b", "c")  # long lateral, polar, chain -- the provider's frame
_SWAP_XY = np.array([[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])  # packer -> provider axes
COMPONENTS = (("aa", 0, 0), ("bb", 1, 1), ("cc", 2, 2), ("ab", 0, 1), ("ba", 1, 0))
TRANSVERSE = ("aa", "bb", "ab", "ba")


def _rotz(X, deg: float) -> np.ndarray:
    c, s = np.cos(np.deg2rad(deg)), np.sin(np.deg2rad(deg))
    return np.stack([c * X[..., 0] - s * X[..., 1], s * X[..., 0] + c * X[..., 1], X[..., 2]], axis=-1)


def chain_frame_coords(params, latn, P_chain, s: int) -> np.ndarray:
    """Undo :meth:`polyfind.pack.CrystalPacker._place` for chain ``s``.

    Placed coordinates of one chain -> that chain's own frame (the one its
    :class:`~polyfind.forcefield.FluxTopology` was resolved in, with the repeat along ``+z``).
    Chain 2 is the chain rotated by ``phi2`` and offset by ``(a + b) / 2 + dz``, after the
    optional ``(x, -y, -z)`` flip, which is its own inverse.
    """
    a, b, gam, phi1, phi2, dz, flip = (float(v) for v in params)
    if s == 0:
        return _rotz(np.asarray(P_chain, dtype=float), -phi1)
    off = 0.5 * (latn[0] + latn[1]) + np.array([0.0, 0.0, dz])
    X = _rotz(np.asarray(P_chain, dtype=float) - off[None, :], -phi2)
    if flip > 0.5:
        X = X * np.array([1.0, -1.0, -1.0])
    return X


def cell_dipole(packer, params, Pn, latn, moved_chain: int | None = None) -> np.ndarray:
    """``mu`` (e.A) of the cell at placed coordinates ``Pn`` (N, 3), induced dipoles included.

    With charge flux, ``moved_chain`` says which chain's charges must be re-derived from its
    (displaced) geometry; every other chain keeps the packer's own charges, which are exactly
    what the flux gives at the undisplaced geometry (asserted by :func:`born_charges`).
    """
    n = packer.n
    cz = float(latn[2, 2])
    q = np.array(packer._q_cell, dtype=float)
    flux = getattr(packer, "_flux", None)
    if flux is not None and moved_chain is not None:
        sl = slice(moved_chain * n, (moved_chain + 1) * n)
        q[sl] = flux.charges(chain_frame_coords(params, latn, Pn[sl], moved_chain), cz)[0]
    mu = q @ Pn
    if packer.polarizable is not None:
        mu = mu + packer._polarize(Pn, q, latn, float(params[6]))[1].sum(axis=0)
    return mu


def atom_labels(chain) -> list:
    """One label per atom of the repeat: the element, and for carbon its pendants -- ``C(F2)``."""
    from .pack import _repeat_images

    tpl, where = _repeat_images(chain)
    adj: dict = {}
    for a, b in tpl.bonds:
        adj.setdefault(int(a), []).append(int(b))
        adj.setdefault(int(b), []).append(int(a))
    labels = [None] * chain.n_atoms
    for idx, (p, img) in where.items():
        if img != 0:
            continue
        e = tpl.elements[idx]
        if e == "C":
            pend = sorted(tpl.elements[k] for k in adj.get(idx, []) if tpl.elements[k] != "C")
            labels[p] = "C(" + "".join(f"{x}{pend.count(x) if pend.count(x) > 1 else ''}" for x in sorted(set(pend))) + ")" if pend else "C"
        else:
            labels[p] = e
    if any(v is None for v in labels):
        raise RuntimeError("could not label every atom of the repeat")
    return labels


@dataclass
class BornCharges:
    """``Z[i]`` (3, 3) per placed atom, ``Z[i, a, b] = d mu_a / d u_{i b}`` in e, packer frame."""

    Z: np.ndarray  # (N, 3, 3)
    elements: list
    labels: list
    positions: np.ndarray  # (N, 3) placed coordinates
    n_chain: int  # atoms per chain
    h: float

    @property
    def asr(self) -> np.ndarray:
        """``sum_i Z_i``: identically zero for this model up to the finite difference."""
        return self.Z.sum(axis=0)

    def canonical(self) -> "CanonicalBorn":
        return canonical(self)


def born_charges(packer, params, h: float = 1e-4) -> BornCharges:
    """The Born effective charges of every atom of the cell at ``params`` (7-vector).

    Central differences of :func:`cell_dipole` in the placed coordinates, one atom of one
    chain at a time -- not the repeat coordinate, which would move the atom's copy in every
    chain at once.
    """
    packer._require_neutral("the Born effective charges")
    params = np.asarray(params, dtype=float).reshape(7)
    P, lat = packer._place(params[None])
    Pn = np.asarray(bk.to_numpy(P), dtype=float)[0]
    latn = np.asarray(bk.to_numpy(lat), dtype=float)[0]
    N, n = packer.N, packer.n
    flux = getattr(packer, "_flux", None)
    if flux is not None:
        # the undisplaced chain, read back through the inverse placement, must reproduce the
        # packer's own charges -- otherwise the frame inversion is wrong and so is everything
        for s in range(packer.n_chains):
            sl = slice(s * n, (s + 1) * n)
            q = flux.charges(chain_frame_coords(params, latn, Pn[sl], s), float(latn[2, 2]))[0]
            err = float(np.abs(q - packer._q_cell[sl]).max())
            if err > 1e-9:
                raise RuntimeError(f"chain {s}: charges re-derived from the placed geometry differ "
                                   f"from the packer's by {err:.2e} e")
    Z = np.zeros((N, 3, 3))
    for i in range(N):
        s = i // n
        for b in range(3):
            Pp, Pm = Pn.copy(), Pn.copy()
            Pp[i, b] += h
            Pm[i, b] -= h
            Z[i, :, b] = (cell_dipole(packer, params, Pp, latn, s)
                          - cell_dipole(packer, params, Pm, latn, s)) / (2.0 * h)
    return BornCharges(Z=Z, elements=list(packer.elements) * packer.n_chains,
                       labels=atom_labels(packer.chain) * packer.n_chains, positions=Pn,
                       n_chain=n, h=float(h))


# --- the provider's axes, and the representative of each symmetry-equivalent set -----------
def _carbon_of(elements, positions, i: int, n_chain: int, cz: float) -> int:
    """The carbon a pendant atom ``i`` is bonded to (nearest carbon of its own chain, z-images included)."""
    s = i // n_chain
    best, dbest = -1, np.inf
    for j in range(s * n_chain, (s + 1) * n_chain):
        if elements[j] != "C":
            continue
        for k in (-1, 0, 1):
            d = float(np.linalg.norm(positions[i] - positions[j] - np.array([0.0, 0.0, k * cz])))
            if d < dbest:
                best, dbest = j, d
    if best < 0 or dbest > 2.0:
        raise RuntimeError(f"atom {i} ({elements[i]}) has no carbon within 2 A")
    return best


@dataclass
class CanonicalBorn:
    """Per-type Born tensors in the provider's ``(a, b, c)`` axes, symmetry-equivalent atoms averaged.

    ``mean[label]`` is the (3, 3) tensor of the representative atom of that type; ``spread``
    the largest deviation of any equivalent atom from it after the reflections, which is the
    symmetry check (it must be at the finite-difference floor for a symmetric cell);
    ``count`` how many atoms of the cell carry the label.
    """

    mean: dict
    spread: dict
    count: dict
    per_atom: np.ndarray  # (N, 3, 3), every atom after its own reflection

    def component(self, label: str, name: str) -> float:
        for nm, a, b in COMPONENTS:
            if nm == name:
                return float(self.mean[label][a, b])
        raise KeyError(name)


def _to_reference_axes(Z_packer: np.ndarray, elements, positions, n_chain: int, cz: float) -> np.ndarray:
    """Every atom's tensor in the provider's axes, reflected onto its representative (N, 3, 3)."""
    N = len(elements)
    Zd = np.einsum("ab,nbc,cd->nad", _SWAP_XY, Z_packer, _SWAP_XY)
    pos = positions @ _SWAP_XY  # (a, b, c)
    carbon = {i: _carbon_of(elements, positions, i, n_chain, cz) for i in range(N) if elements[i] != "C"}
    # the polar sign: + from a carbon towards its heaviest (non-hydrogen) pendant, as in the
    # provider's cell, where the fluorines sit at +b of their carbon; hydrogen-only chains
    # (polyethylene) have no polar axis and keep the packer's sign
    heavy = [i for i in carbon if elements[i] != "H"]
    s_pol = 1.0
    if heavy:
        s_pol = float(np.sign(np.mean([pos[i, 1] - pos[carbon[i], 1] for i in heavy])) or 1.0)
    out = np.empty_like(Zd)
    for i in range(N):
        s_long = 1.0
        if i in carbon:
            off = pos[i, 0] - pos[carbon[i], 0]
            s_long = 1.0 if abs(off) < 1e-6 else float(np.sign(off))
        R = np.diag([s_long, s_pol, 1.0])
        out[i] = R @ Zd[i] @ R
    return out


def canonical_from(born: BornCharges, cz: float) -> CanonicalBorn:
    """:class:`CanonicalBorn` of ``born``; ``cz`` is the repeat along the chain axis (A)."""
    per = _to_reference_axes(born.Z, born.elements, born.positions, born.n_chain, cz)
    mean, spread, count = {}, {}, {}
    for lab in dict.fromkeys(born.labels):
        idx = [i for i, l in enumerate(born.labels) if l == lab]
        m = per[idx].mean(axis=0)
        mean[lab] = m
        spread[lab] = float(np.abs(per[idx] - m[None]).max())
        count[lab] = len(idx)
    return CanonicalBorn(mean=mean, spread=spread, count=count, per_atom=per)


def canonical(born: BornCharges, cz: float | None = None) -> CanonicalBorn:
    """:class:`CanonicalBorn` of ``born``.  ``cz`` (the chain repeat, A) defaults to the z-extent
    of the placed cell plus one bond, which is enough to pick each pendant's own carbon."""
    if cz is None:
        z = born.positions[:, 2]
        cz = float(z.max() - z.min()) + 1.6
    return canonical_from(born, cz)


# --- the provider's reference ---------------------------------------------------------------
@dataclass
class BornReference:
    """The periodic-DFPT Born tensors, per type, in the provider's own axes and representative.

    ``corrected`` carries the acoustic-sum-corrected tensors and ``raw`` the ones before the
    correction; the difference is the whole of what the correction did, and for this data set
    it lives almost entirely in ``cc``.  ``asr_raw`` is the raw acoustic sum.
    """

    corrected: dict
    raw: dict
    count: dict
    asr_raw: np.ndarray
    asr_corrected: np.ndarray
    eps_inf: np.ndarray
    protocol: str
    path: str

    def component(self, label: str, name: str, corrected: bool = True) -> float:
        table = self.corrected if corrected else self.raw
        for nm, a, b in COMPONENTS:
            if nm == name:
                return float(table[label][a, b])
        raise KeyError(name)


def _read_xyz(path: str) -> tuple:
    with open(path) as fh:
        lines = fh.read().splitlines()
    n = int(lines[0].split()[0])
    el, X = [], []
    for ln in lines[2:2 + n]:
        parts = ln.split()
        el.append(parts[0])
        X.append([float(v) for v in parts[1:4]])
    return el, np.array(X, dtype=float)


def load_reference(path: str) -> BornReference:
    """Read ``born_results.json`` (and the geometry it names, from the same directory).

    The tensors are stored as printed by ``ph.x``'s ``(d Force / dE)`` block, whose row is the
    field direction and column the displacement: ``[alpha][beta] = dP_alpha / du_beta``, the
    same index order :func:`born_charges` uses.  Each pendant atom is reflected onto the
    representative on the ``+a`` side of its carbon before averaging, exactly as
    :func:`canonical` does for the model.
    """
    with open(path) as fh:
        d = json.load(fh)
    el, X = _read_xyz(os.path.join(os.path.dirname(path), d["geometry"]))
    atoms = d["atoms"]
    if len(atoms) != len(el):
        raise ValueError("born_results.json and its geometry disagree on the atom count")
    raw = np.array([a["born_effective_charge_e"] for a in atoms], dtype=float)
    cor = np.array([a["born_effective_charge_asr_corrected_e"] for a in atoms], dtype=float)
    # label and reflect as the model is: nearest carbon, polar sign from the heavy pendants
    N = len(el)
    carbon = {}
    for i in range(N):
        if el[i] == "C":
            continue
        dist = [np.linalg.norm(X[i] - X[j]) if el[j] == "C" else np.inf for j in range(N)]
        carbon[i] = int(np.argmin(dist))
    labels = []
    for i in range(N):
        if el[i] != "C":
            labels.append(el[i])
            continue
        pend = sorted(el[j] for j, c in carbon.items() if c == i)
        labels.append("C(" + "".join(f"{x}{pend.count(x) if pend.count(x) > 1 else ''}" for x in sorted(set(pend))) + ")" if pend else "C")
    heavy = [i for i in carbon if el[i] != "H"]
    s_pol = float(np.sign(np.mean([X[i, 1] - X[carbon[i], 1] for i in heavy])) or 1.0) if heavy else 1.0
    tables = {"raw": {}, "corrected": {}}
    count = {}
    for key, T in (("raw", raw), ("corrected", cor)):
        per = np.empty_like(T)
        for i in range(N):
            s_long = 1.0
            if i in carbon:
                off = X[i, 0] - X[carbon[i], 0]
                s_long = 1.0 if abs(off) < 1e-6 else float(np.sign(off))
            R = np.diag([s_long, s_pol, 1.0])
            per[i] = R @ T[i] @ R
        for lab in dict.fromkeys(labels):
            idx = [i for i, l in enumerate(labels) if l == lab]
            tables[key][lab] = per[idx].mean(axis=0)
            count[lab] = len(idx)
    return BornReference(corrected=tables["corrected"], raw=tables["raw"], count=count,
                         asr_raw=np.array(d["acoustic_sum_tensor_e"], dtype=float),
                         asr_corrected=np.array(d["acoustic_sum_tensor_asr_corrected_e"], dtype=float),
                         eps_inf=np.array(d["dielectric_tensor_electronic"], dtype=float),
                         protocol=str(d.get("protocol", "")), path=path)


def comparison_table(model: CanonicalBorn, ref: BornReference | None = None, labels=None,
                     static=None) -> str:
    """Rows ``type component model [static] [ref corrected, ref raw]`` in e."""
    labels = list(model.mean) if labels is None else list(labels)
    head = f"{'type':8s} {'comp':4s} {'model':>9s}"
    if static is not None:
        head += f" {'static':>8s}"
    if ref is not None:
        head += f" {'DFPT asr':>9s} {'DFPT raw':>9s} {'diff':>8s}"
    lines = [head]
    for lab in labels:
        for nm, a, b in COMPONENTS:
            if nm in ("ab", "ba") and abs(model.mean[lab][a, b]) < 5e-4 and (
                    ref is None or lab not in ref.corrected or abs(ref.corrected[lab][a, b]) < 5e-4):
                continue
            row = f"{lab:8s} {nm:4s} {model.mean[lab][a, b]:+9.4f}"
            if static is not None:
                row += f" {static.get(lab, float('nan')):+8.4f}" if a == b else f" {'':8s}"
            if ref is not None and lab in ref.corrected:
                rc, rr = ref.corrected[lab][a, b], ref.raw[lab][a, b]
                row += f" {rc:+9.4f} {rr:+9.4f} {model.mean[lab][a, b] - rc:+8.4f}"
            lines.append(row)
        lines.append(f"{lab:8s} spread over {model.count[lab]} equivalent atoms: {model.spread[lab]:.1e} e")
    return "\n".join(lines)
