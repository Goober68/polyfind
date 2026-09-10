"""Does a packed periodic cell have the bond graph it was meant to have?

A packing search scores a cell against a nonbonded potential.  It does *not* check that
the cell it returns is still the molecule it started from: two chains can be placed close
enough that a downstream consumer -- anything that derives bonds from interatomic
distances, which is what every machine-learned-potential driver and every structure
validator does -- sees covalent bonds between them, or fails to see one that should be
there.  That is not a scoring error; an energy can be perfectly finite while the topology
is wrong, and the wrongness only shows up later as a relaxation that changes connectivity.

So this module asks the question a consumer's detector asks, and asks it the way a
consumer asks it: a pair is bonded when

    d_ij <= scale * (r_i + r_j)

with ``r`` a covalent radius.  The useful output is not a yes/no at one ``scale`` -- a
deliverable that passes at 1.2 and fails at 1.3 is not a deliverable -- but the **window of
scales over which the detected graph is exactly the intended one**:

* :attr:`TopologyReport.max_bond_ratio`, the worst intended bond's ``d / (r_i + r_j)``:
  below this an intended bond goes undetected (a "lost bond");
* :attr:`TopologyReport.min_nonbond_ratio`, the closest pair that is *not* an intended
  bond: above this a spurious bond appears (a "new bond"), and when that pair is on two
  different chains it is a new *interchain* bond, the failure mode that merges separate
  chains into one connected component.

Anything strictly between the two reproduces the intended graph exactly, so a wide window
means no reasonable choice of detector disagrees.  The intended graph itself comes from the
:class:`~polyfind.polymers.Polymer` definition (:func:`repeat_bond_graph`) rather than from
the code that produced the coordinates, so the check is independent of the builder it is
checking.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass

import numpy as np

from .polymers import Polymer

# Covalent radii in A (Cordero et al., Dalton Trans. 2008): the same table ASE's
# ``ase.data.covalent_radii`` carries, and therefore the same numbers a consumer's
# distance-based bond detector will be using.
COVALENT_RADIUS: dict[str, float] = {
    "H": 0.31,
    "C": 0.76,
    "N": 0.71,
    "O": 0.66,
    "F": 0.57,
    "Cl": 1.02,
}


def repeat_bond_graph(polymer: Polymer) -> tuple[np.ndarray, list[tuple[int, int, int]], int]:
    """``(starts, bonds, n_atoms)`` for one crystallographic repeat, in block-local indices.

    ``starts[k]`` is the local index of backbone atom ``k`` of the repeat; its pendant atoms
    follow it in pendant order, which is exactly the layout
    :func:`polyfind.chain.build_chain` produces and :func:`polyfind.pack.periodic_chain`
    slices out.  ``bonds`` is ``(i, j, dk)`` with ``dk`` the number of ``c`` translations
    between the two atoms: ``0`` for every bond inside the repeat and ``1`` for the single
    backbone bond that closes the chain onto its own next repeat.

    Derived from the polymer's chemistry alone -- which monomers, which pendants, which
    intra-pendant bonds -- so it is the *intended* graph rather than a re-reading of the
    built coordinates.
    """
    B = polymer.bonds_per_repeat
    widths = [1 + polymer.backbone[k].n_pendant_atoms for k in range(B)]
    starts = np.concatenate([[0], np.cumsum(widths)]).astype(int)
    bonds: list[tuple[int, int, int]] = []
    for k in range(B):
        base = int(starts[k])
        if k + 1 < B:
            bonds.append((base, int(starts[k + 1]), 0))
        else:
            bonds.append((base, int(starts[0]), 1))  # closes onto the next repeat, +c away
        off = base + 1
        for pendant in polymer.backbone[k].pendants:
            bonds.append((base, off, 0))
            bonds.extend((off + i, off + j, 0) for i, j in pendant.bonds)
            off += len(pendant)
    return starts[:-1], bonds, int(starts[-1])


@dataclass
class Cell:
    """A periodic all-atom cell, with each atom's chain and within-chain index."""

    elements: list[str]
    coords: np.ndarray  # (n, 3)
    lattice: np.ndarray  # (3, 3), rows are the lattice vectors
    chain_of: np.ndarray  # (n,) which chain each atom belongs to
    local_of: np.ndarray  # (n,) index of the atom within its chain's repeat
    n_per_chain: int
    # per chain: True when the chain was placed antiparallel, i.e. with its z reversed.  The
    # bond that closes a chain onto its own next repeat then runs towards -c rather than +c,
    # so the intended graph's image labels flip for it (:func:`intended_edges`).  This is
    # construction data, not a measurement: it comes from the packer's ``flip``.
    reversed_of: np.ndarray = None  # type: ignore[assignment]

    def __post_init__(self) -> None:
        if self.reversed_of is None:
            self.reversed_of = np.zeros(int(self.chain_of.max()) + 1, dtype=bool)

    @property
    def n_atoms(self) -> int:
        return len(self.elements)

    @property
    def n_chains(self) -> int:
        return int(self.chain_of.max()) + 1

    @property
    def volume(self) -> float:
        return float(abs(np.linalg.det(self.lattice)))

    def formula(self) -> str:
        c = Counter(self.elements)
        return "".join(f"{e}{c[e]}" for e in ("C", "H", "F", "N", "O", "Cl") if c[e])

    def to_extxyz(self, comment_extra: str = "") -> str:
        """Extended XYZ with the lattice on the comment line (ASE/extxyz convention)."""
        lat = " ".join(f"{v:.8f}" for v in self.lattice.reshape(-1))
        head = f'Lattice="{lat}" Properties=species:S:1:pos:R:3 pbc="T T T"'
        if comment_extra:
            head += " " + comment_extra
        lines = [str(self.n_atoms), head]
        for e, (x, y, z) in zip(self.elements, self.coords):
            lines.append(f"{e:<2s} {x:16.8f} {y:16.8f} {z:16.8f}")
        return "\n".join(lines) + "\n"


def build_cell(chain, params, nx: int = 1, ny: int = 1, packer=None) -> Cell:
    """The packed cell of ``params`` (:attr:`~polyfind.pack.PackResult.params`), tiled.

    ``nx``/``ny`` tile the two-chain cell in the ab plane, which is how an eight-chain cell
    is reached from a packer that places two: ``nx = ny = 2`` gives four times the packer's
    chains in a cell of ``(nx a, ny b, c)``.  **The tiling is an exact replication**, so the
    four copies of the two-chain motif are not independent -- a chain cannot slip, rotate or
    register differently from its own image, because it *is* its own image.  Tiling reaches
    the size and the composition; it does not add packing degrees of freedom.
    """
    from .pack import CrystalPacker

    pk = packer or CrystalPacker(chain, n_chains=2)
    params = np.asarray(params, dtype=float)
    P, lat = pk._place(params[None])
    base = np.asarray(P[0], dtype=float)
    lat = np.asarray(lat[0], dtype=float)
    n = chain.n_atoms
    k_chains = base.shape[0] // n
    flipped = float(params[6]) > 0.5
    elements, coords, chain_of, local_of, rev = [], [], [], [], []
    cid = 0
    for ix in range(nx):
        for iy in range(ny):
            shift = ix * lat[0] + iy * lat[1]
            for k in range(k_chains):
                coords.append(base[k * n : (k + 1) * n] + shift)
                elements.extend(chain.elements)
                chain_of.extend([cid] * n)
                local_of.extend(range(n))
                rev.append(k == 1 and flipped)  # ``_place`` flips chain 2 only
                cid += 1
    super_lat = np.array([nx * lat[0], ny * lat[1], lat[2]])
    return Cell(elements=elements, coords=np.concatenate(coords), lattice=super_lat,
                chain_of=np.array(chain_of), local_of=np.array(local_of), n_per_chain=n,
                reversed_of=np.array(rev, dtype=bool))


@dataclass
class TopologyReport:
    n_atoms: int
    n_chains: int
    n_per_chain: int
    formula: str
    n_intended_bonds: int
    max_bond_ratio: float
    max_bond_pair: tuple
    min_nonbond_ratio: float
    min_nonbond_pair: tuple
    min_interchain_ratio: float
    min_interchain_distance: float
    min_interchain_pair: tuple
    min_intra_nonbond_distance: float
    min_interchain_by_element: dict
    components: dict  # scale -> sorted component sizes
    missing: dict  # scale -> intended bonds the detector would not see
    extra: dict  # scale -> detected bonds that are not intended
    extra_interchain: dict  # scale -> of those, how many join two different chains

    @property
    def safe_scale_window(self) -> tuple[float, float]:
        """Every covalent scale strictly inside this range reproduces the intended graph."""
        return (self.max_bond_ratio, self.min_nonbond_ratio)

    @property
    def ok(self) -> bool:
        """The detected graph is the intended one over a non-empty range of scales, every
        chain is its own connected component at every scale tried, and no interchain bond
        appears at any of them."""
        if not (self.min_nonbond_ratio > self.max_bond_ratio):
            return False
        want = sorted([self.n_per_chain] * self.n_chains)
        return all(self.components[s] == want and not self.missing[s] and not self.extra[s]
                   for s in self.components)

    def describe(self) -> str:
        lo, hi = self.safe_scale_window
        lines = [
            f"{self.formula}: {self.n_atoms} atoms, {self.n_chains} chains of {self.n_per_chain}, "
            f"{self.n_intended_bonds} intended bonds",
            f"  worst intended bond       d/(ri+rj) = {self.max_bond_ratio:.4f}   {self.max_bond_pair}",
            f"  closest non-bonded pair   d/(ri+rj) = {self.min_nonbond_ratio:.4f}   {self.min_nonbond_pair}",
            f"  closest interchain pair   d/(ri+rj) = {self.min_interchain_ratio:.4f}   "
            f"d = {self.min_interchain_distance:.4f} A   {self.min_interchain_pair}",
            f"  closest intrachain non-bonded pair  d = {self.min_intra_nonbond_distance:.4f} A",
            f"  intended graph reproduced for every covalent scale in ({lo:.4f}, {hi:.4f})",
            "  minimum interchain distance by element pair (A):",
        ]
        for pair in sorted(self.min_interchain_by_element, key=lambda p: self.min_interchain_by_element[p]):
            lines.append(f"    {pair[0]:>2s}-{pair[1]:<2s} {self.min_interchain_by_element[pair]:7.4f}")
        lines.append("  detected graph per covalent scale:")
        for s in sorted(self.components):
            comp = self.components[s]
            tag = "OK" if (comp == sorted([self.n_per_chain] * self.n_chains)
                           and not self.missing[s] and not self.extra[s]) else "FAIL"
            lines.append(f"    scale {s:4.2f}: components {comp}  lost {self.missing[s]}  "
                         f"new {self.extra[s]} (interchain {self.extra_interchain[s]})  {tag}")
        lines.append(f"  verdict: {'PASS' if self.ok else 'FAIL'}")
        return "\n".join(lines)


def _image_range(lattice: np.ndarray, cut: float) -> list[int]:
    """How many images along each axis are needed to cover ``cut``."""
    out = []
    for i in range(3):
        j, k = (i + 1) % 3, (i + 2) % 3
        nrm = np.cross(lattice[j], lattice[k])
        h = abs(lattice[i] @ nrm) / np.linalg.norm(nrm)  # interplanar spacing
        out.append(int(np.ceil(cut / h)))
    return out


def _key(i: int, j: int, img: tuple) -> tuple:
    """Canonical edge key: lower atom index first, image measured from it.

    One physical pair therefore gets one key whichever image loop reaches it.  A pair of an
    atom with its own image (``i == j``) has no index order to break the tie, so the image
    is normalised to the lexicographically larger of ``+img`` and ``-img``.
    """
    if i < j:
        return (i, j, img)
    neg = (-img[0], -img[1], -img[2])
    if i > j:
        return (j, i, neg)
    return (i, j, max(img, neg))


def intended_edges(cell: Cell, polymer: Polymer) -> set:
    """The bonds the cell is *supposed* to have, as canonical ``(i, j, image)`` keys."""
    _, bonds, n_rep = repeat_bond_graph(polymer)
    if n_rep != cell.n_per_chain:
        raise ValueError(f"{polymer.name!r} has {n_rep} atoms per repeat but the cell has "
                         f"{cell.n_per_chain} per chain")
    out = set()
    for t in range(cell.n_chains):
        off = t * cell.n_per_chain
        sgn = -1 if cell.reversed_of[t] else 1
        for i, j, dk in bonds:
            out.add(_key(off + i, off + j, (0, 0, sgn * dk)))
    return out


def check_topology(cell: Cell, polymer: Polymer, scales=(1.05, 1.1, 1.15, 1.2, 1.3),
                   cut: float = 6.0) -> TopologyReport:
    """Compare ``cell``'s distance-derived bond graph against the intended one.

    Every pair within ``cut`` of every periodic image is examined, so nothing is hidden by a
    minimum-image approximation.  ``scales`` are the covalent-radius multipliers at which
    the detected graph is compared, edge by edge, with the intended one.
    """
    want = intended_edges(cell, polymer)
    X, n = cell.coords, cell.n_atoms
    elements = np.array(cell.elements)
    r = np.array([COVALENT_RADIUS[e] for e in cell.elements])
    rsum = r[:, None] + r[None, :]
    lat = cell.lattice
    na, nb, nc = _image_range(lat, cut)
    same_chain = cell.chain_of[:, None] == cell.chain_of[None, :]
    els = sorted(set(cell.elements))
    code = np.array([els.index(e) for e in cell.elements])
    el_masks = {(a, b): (code == a)[:, None] & (code == b)[None, :]
                for a in range(len(els)) for b in range(a, len(els))}
    el_min = {k: np.inf for k in el_masks}

    detected = {s: set() for s in scales}
    best_non = (np.inf, None)
    best_inter = (np.inf, None, np.inf)
    best_intra_non = np.inf
    for ia in range(-na, na + 1):
        for ib in range(-nb, nb + 1):
            for ic in range(-nc, nc + 1):
                img = (ia, ib, ic)
                shift = ia * lat[0] + ib * lat[1] + ic * lat[2]
                d = np.linalg.norm(X[:, None, :] - (X[None, :, :] + shift), axis=2)
                if img == (0, 0, 0):
                    np.fill_diagonal(d, np.inf)
                # blank out the intended bonds so what is left is everything that must NOT
                # look like a bond; done by direct indexing rather than by testing pairs
                dn = d.copy()
                for i, j, e_img in ((i, j, e) for (i, j, e) in want if e == img or e == tuple(-v for v in img)):
                    if e_img == img:
                        dn[i, j] = np.inf
                    if tuple(-v for v in e_img) == img:
                        dn[j, i] = np.inf
                near = dn <= cut
                if near.any():
                    ratio = dn / rsum
                    k = int(np.argmin(np.where(near, ratio, np.inf)))
                    i, j = divmod(k, n)
                    if ratio[i, j] < best_non[0]:
                        best_non = (float(ratio[i, j]),
                                    (cell.elements[i], cell.elements[j], int(i), int(j), img, round(float(dn[i, j]), 4)))
                    inter = near & ~same_chain
                    if inter.any():
                        k = int(np.argmin(np.where(inter, dn, np.inf)))
                        i, j = divmod(k, n)
                        if dn[i, j] < best_inter[2]:
                            best_inter = (float(ratio[i, j]),
                                          (cell.elements[i], cell.elements[j], int(i), int(j), img),
                                          float(dn[i, j]))
                        dw = np.where(inter, dn, np.inf)
                        for key, m in el_masks.items():
                            mm = m | m.T if key[0] != key[1] else m
                            v = dw[mm]
                            if v.size:
                                el_min[key] = min(el_min[key], float(v.min()))
                    intra = near & same_chain
                    if intra.any():
                        best_intra_non = min(best_intra_non, float(np.where(intra, dn, np.inf).min()))
                for s in scales:
                    hit = d <= s * rsum
                    for i, j in zip(*np.nonzero(hit)):
                        detected[s].add(_key(int(i), int(j), img))

    # worst intended bond, straight off the edge list
    bond_ratio, bond_pair = -np.inf, None
    for i, j, img in want:
        shift = img[0] * lat[0] + img[1] * lat[1] + img[2] * lat[2]
        dd = float(np.linalg.norm(X[i] - (X[j] + shift)))
        rr = dd / (r[i] + r[j])
        if rr > bond_ratio:
            bond_ratio, bond_pair = rr, (cell.elements[i], cell.elements[j], int(i), int(j), img, round(dd, 4))

    by_el = {(els[a], els[b]): v for (a, b), v in el_min.items() if np.isfinite(v)}

    components, missing, extra, extra_inter = {}, {}, {}, {}
    for s in scales:
        edges = detected[s]
        components[s] = _component_sizes(n, edges)
        missing[s] = len(want - edges)
        ex = edges - want
        extra[s] = len(ex)
        extra_inter[s] = sum(1 for i, j, _ in ex if cell.chain_of[i] != cell.chain_of[j])
    return TopologyReport(
        n_atoms=n, n_chains=cell.n_chains, n_per_chain=cell.n_per_chain, formula=cell.formula(),
        n_intended_bonds=len(want),
        max_bond_ratio=float(bond_ratio), max_bond_pair=bond_pair,
        min_nonbond_ratio=best_non[0], min_nonbond_pair=best_non[1],
        min_interchain_ratio=best_inter[0], min_interchain_pair=best_inter[1],
        min_interchain_distance=best_inter[2],
        min_intra_nonbond_distance=best_intra_non,
        min_interchain_by_element=by_el,
        components=components, missing=missing, extra=extra, extra_interchain=extra_inter,
    )


def _component_sizes(n: int, edges) -> list[int]:
    """Sizes of the connected components of the periodic bond graph (union-find)."""
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for i, j, _ in edges:
        ri, rj = find(i), find(j)
        if ri != rj:
            parent[ri] = rj
    sizes: dict[int, int] = {}
    for x in range(n):
        rx = find(x)
        sizes[rx] = sizes.get(rx, 0) + 1
    return sorted(sizes.values())
