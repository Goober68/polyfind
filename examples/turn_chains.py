"""Every chain of a periodic cell turned about its own axis by an angle, the cell kept: a start for a relaxation in
another setting of the chains (docs/MLIP_SECOND_OPINION.md 3.23; examples/registry_cells.py has the turn and builds whole
sets of such starts).

    python -I examples/turn_chains.py FRAMES.extxyz FRAME DEGREES OUT.xyz NEWNAME

The chains are the connected components of the bonds; each is unwrapped along its bonds and turned about the line through
its centre along the third cell vector.  Prints the closest contact between chains in the start."""
import importlib.util
import sys
from pathlib import Path

from ase.io import read, write

spec = importlib.util.spec_from_file_location("registry_cells", Path(__file__).resolve().parent / "registry_cells.py")
registry = importlib.util.module_from_spec(spec)
spec.loader.exec_module(registry)

src, frame, angle, out, name = sys.argv[1], sys.argv[2], float(sys.argv[3]), sys.argv[4], sys.argv[5]
a = next(x for x in read(src, index=":") if str(x.info["frame"]) == frame)
b = registry.turned(a, angle)
b.info = {k: v for k, v in a.info.items() if k in ("family", "n_monomers", "kmesh", "conformation")}
b.info["frame"] = name
b.info["note"] = f"{frame} with every chain turned {angle:+.0f} degrees about its own axis, the cell kept: a start, not relaxed"
for k in list(b.arrays):
    if k not in ("numbers", "positions"):
        del b.arrays[k]
b.calc = None
write(out, b, format="extxyz")
print(name, "written; the closest contact between chains", registry.closest_between_chains(b, reach=2.2))
