"""Every chain of a periodic cell turned about its own axis by an angle, the cell kept: a start for a relaxation in
another setting of the chains (docs/MLIP_SECOND_OPINION.md 3.23).

    python -I examples/turn_chains.py FRAMES.extxyz FRAME DEGREES OUT.xyz NEWNAME

The chains are the connected components of the bonds; each is unwrapped along its bonds and turned about the line through
its centre along the third cell vector.  Prints the closest contact between chains in the start."""
import sys
import numpy as np
from ase.io import read, write
from ase.neighborlist import neighbor_list
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

src, frame, angle, out, name = sys.argv[1], sys.argv[2], float(sys.argv[3]), sys.argv[4], sys.argv[5]
a = next(x for x in read(src, index=":") if str(x.info["frame"]) == frame)
n = len(a)
i, j, S = neighbor_list("ijS", a, {("C", "C"): 1.75, ("C", "H"): 1.3, ("C", "F"): 1.6, ("C", "Cl"): 2.0, ("C", "N"): 1.4})
ncomp, lab = connected_components(coo_matrix((np.ones(len(i)), (i, j)), shape=(n, n)), directed=False)
assert ncomp == 2, ncomp
cell = a.cell.array
c = cell[2] / np.linalg.norm(cell[2])
pos = a.positions.copy()
new = pos.copy()
nb = {k: [] for k in range(n)}
for ii, jj, s in zip(i, j, S):
    nb[ii].append((jj, s))
th = np.radians(angle)
for comp in range(ncomp):
    atoms = np.flatnonzero(lab == comp)
    un = {atoms[0]: pos[atoms[0]]}          # unwrap along the bonds: a spanning tree, each atom placed once
    stack = [atoms[0]]
    while stack:
        k = stack.pop()
        for m, s in nb[k]:
            if m not in un:
                un[m] = un[k] + (pos[m] + s @ cell - pos[k])
                stack.append(m)
    r = np.array([un[k] for k in atoms])
    p0 = r.mean(axis=0)
    d = r - p0
    par = np.outer(d @ c, c)
    perp = d - par
    turned = perp * np.cos(th) + np.cross(c, perp) * np.sin(th)
    new[atoms] = p0 + par + turned
b = a.copy()
b.positions = new
b.wrap()
b.info = {k: v for k, v in a.info.items() if k in ("family", "n_monomers", "kmesh", "conformation")}
b.info["frame"] = name
b.info["note"] = f"{frame} with every chain turned {angle:+.0f} degrees about its own axis, the cell kept: a start, not relaxed"
for k in list(b.arrays):
    if k not in ("numbers", "positions"):
        del b.arrays[k]
b.calc = None
write(out, b, format="extxyz")
i2, j2, d2 = neighbor_list("ijd", b, 2.2)
close = min(dd for ii, jj, dd in zip(i2, j2, d2) if lab[ii] != lab[jj]) if any(lab[ii] != lab[jj] for ii, jj in zip(i2, j2)) else None
print(name, "written; the closest contact between chains", close)
