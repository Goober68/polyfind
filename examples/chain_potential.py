"""The electrostatic potential around one chain: polyfind's charges and induced dipoles against DFT's.

The rigid cells of examples/interchain_frames.py say that polyfind's electrostatic energy between
chains is right for plain PVDF and more than twice too large for the chain in which every second
unit is reversed (-CH2-CF2-CF2-CH2-), and a fit of energies cannot say which part of the charge
model is at fault: every way of weakening it fits about as well, and most end on a bound.  The
potential a chain makes around itself is the direct measurement.  This script takes one chain out
of a packing, puts it alone in a wide box, has DFT give the electrostatic potential on a grid
(``pp.x``, ``plot_num = 11``), and compares it, at the distances where a neighbouring chain's atoms
sit, with the potential of polyfind's sources: its point charges and the dipoles they induce
along the chain, summed over the same periodic box.

The charges are split as examples/interchain_wall.py splits them: the CH2 groups' charges, the CF2
groups' charges, and a third pattern polyfind's model does not have, charge moved along each
backbone bond between a CH2 carbon and a CF2 carbon.  Induction is linear, so the potential is
linear in the three multipliers and a least-squares fit gives what DFT's potential asks of each.

    python examples/chain_potential.py frames --out /opt/dft/chainpot
    python examples/chain_potential.py dft --frames /opt/dft/chainpot --pw .../pw.x --pseudo .../pseudo
    python examples/chain_potential.py collect --frames /opt/dft/chainpot
    python examples/chain_potential.py fit --samples docs/predictions/chain_potential_pbe_samples.npz
"""
import argparse
import importlib.util
import json
import os
import subprocess
import time
from pathlib import Path

import numpy as np
from scipy.special import erfc

from polyfind.fitting import PRESET_PARAMETERS


def _load(name):
    spec = importlib.util.spec_from_file_location(name, Path(__file__).with_name(f"{name}.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


IW = _load("interchain_wall")
IF, SR = IW.IF, IW.SR
PRESET = SR.PRESET
PRED = Path("docs/predictions")
COULOMB = 332.0637133  # kcal A / (mol e^2)
RY_TO_KCAL = 313.7547
# The chains: name -> (packing of examples/interchain_frames.py it is taken from, what it is).
CHAINS = {"ht": ("beta", "plain PVDF, all-trans"), "hh": ("hh", "every second unit reversed, all-trans"),
          "tgtg": ("alpha", "plain PVDF, TGTG'")}
PATTERNS = ("s_CH", "s_CF", "d_CC")
BOX = 16.0   # A, the box's lateral edge: the chain's images are this far apart
APART = 30.0  # A, the least lateral edge of the cell polyfind's dipoles are solved in


def isolated(pack: str, table: dict, check: bool = True):
    """One chain of a packing with the others far away: ``(elements, positions, patterns (3, n),
    induced dipoles (3, n, 3) in e A, repeat along the chain in A)`` of its sites, the dipoles
    those each charge pattern induces along the chain itself."""
    model, x0, ch, _ = IF.model_of(pack, table[pack].get("cell"), check)
    x = np.array(x0, dtype=float)
    x[0] *= np.ceil(APART / x[0])
    x[1] *= np.ceil(APART / x[1])
    latn, _, _, Pn, _ = model.build(x)
    if abs(latn[2, 0]) > 1e-9 or abs(latn[2, 1]) > 1e-9:
        raise RuntimeError("the chain axis is not along z")
    n = model.n
    el = [str(e) for e in model.pk.elements]
    pat = IW.charge_patterns(el, np.asarray(model.pk._q_cell, dtype=float)[:n], model.nc)
    fl = model.flips.astype(float)
    dip = np.array([model.pk._polarize(Pn, p, latn, fl)[1][:n] for p in pat])
    return el, Pn[:n].copy(), pat[:, :n], dip, float(latn[2, 2])


def pattern_names(n: int) -> tuple:
    """The names of a chain's charge patterns: three, or four for a chain built with bond sites."""
    return PATTERNS if n == 3 else PATTERNS + ("M_H",)


def ewald_potential(points, sites, q, cell, alpha: float = 0.4, tol: float = 1e-9):
    """The potential (e/A) at ``points`` of point charges ``q`` at ``sites`` repeated on the
    orthorhombic ``cell`` (three edges), the cell neutral, no term at G = 0 (the convention of a
    plane-wave code's Hartree potential)."""
    cell = np.asarray(cell, dtype=float)
    if abs(q.sum()) > 1e-9:
        raise ValueError("the cell is not neutral")
    rc = np.sqrt(-np.log(tol)) / alpha
    out = np.zeros(len(points))
    reach = np.ceil(rc / cell).astype(int)
    shifts = np.array([[i, j, k] for i in range(-reach[0], reach[0] + 1) for j in range(-reach[1], reach[1] + 1)
                       for k in range(-reach[2], reach[2] + 1)], dtype=float) * cell
    for s, qs in zip(sites, q):
        d = points[:, None, :] - (s + shifts)[None, :, :]
        r = np.sqrt((d * d).sum(axis=-1))
        out += qs * np.where(r < rc, erfc(alpha * r) / np.maximum(r, 1e-12), 0.0).sum(axis=1)
    gmax = 2.0 * alpha * np.sqrt(-np.log(tol))
    m = np.ceil(gmax * cell / (2.0 * np.pi)).astype(int)
    G = np.array([[i, j, k] for i in range(-m[0], m[0] + 1) for j in range(-m[1], m[1] + 1) for k in range(-m[2], m[2] + 1)
                  if (i, j, k) != (0, 0, 0)], dtype=float) * (2.0 * np.pi / cell)
    g2 = (G * G).sum(axis=1)
    keep = g2 < gmax * gmax
    G, g2 = G[keep], g2[keep]
    S = (q[None, :] * np.exp(-1j * (G @ sites.T))).sum(axis=1)
    coef = 4.0 * np.pi / np.prod(cell) * np.exp(-g2 / (4.0 * alpha * alpha)) / g2 * S
    for a in range(0, len(points), 4000):
        out[a:a + 4000] += (np.exp(1j * (points[a:a + 4000] @ G.T)) @ coef).real
    return out


def source_potential(points, sites, q, dipoles, cell, h: float = 0.01, **kw):
    """The potential (kcal/mol per e) of charges and point dipoles (e A), the dipoles as pairs of
    charges ``h`` A either side of their sites."""
    pos, ch = [np.asarray(sites, dtype=float)], [np.asarray(q, dtype=float)]
    norm = np.linalg.norm(dipoles, axis=1)
    on = norm > 1e-14
    unit = dipoles[on] / norm[on, None]
    pos += [sites[on] + h * unit, sites[on] - h * unit]
    ch += [norm[on] / (2.0 * h), -norm[on] / (2.0 * h)]
    return COULOMB * ewald_potential(points, np.concatenate(pos), np.concatenate(ch), cell, **kw)


def frames(args):
    from ase import Atoms
    from ase.io import write

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    table = json.loads(Path(args.cells).read_text())
    params = PRESET_PARAMETERS[PRESET]
    if args.sites:
        # the same chains with a charge site on every C-H bond, charged so that the lattice sums carry it:
        # its potential, the dipoles it induces included, is then polyfind's own and not an extra source
        from dataclasses import replace

        params = replace(params, bond_sites=(("H", -1.0),))
        SR.C.SITES_FOR = params
    with params.applied():
        for name in [c for c in args.chains.split(",") if c]:
            el, P, pat, dip, c = isolated(CHAINS[name][0], table, check=not args.sites)
            atom = np.array([e in SR.ATOMS for e in el])
            centre = P[atom, :2].mean(axis=0)
            P = P - np.array([centre[0] - BOX / 2.0, centre[1] - BOX / 2.0, 0.0])
            a = Atoms(symbols=[e for e, k in zip(el, atom) if k], positions=P[atom], cell=[BOX, BOX, c], pbc=True)
            a.info.update(frame=f"chain_{name}", kmesh=f"1 1 {max(2, int(round(20.0 / c)))}", source=f"polyfind {PRESET}, one chain: {CHAINS[name][1]}")
            start = out / f"chain_{name}.start.xyz"
            if args.sites and start.exists():  # the atoms DFT was run on: the sited chain's must be the same
                from ase.io import read
                if np.abs(read(start).positions - a.positions).max() > 1e-8:
                    raise RuntimeError(f"{name}: the chain with bond sites does not have its atoms where the chain without them does")
            else:
                write(start, a, format="extxyz")
            np.savez(out / f"chain_{name}.model{'_sites' if args.sites else ''}.npz", elements=np.array(el), sites=P, patterns=pat, dipoles=dip,
                     cell=np.array([BOX, BOX, c]))
            mu = (pat[:2].sum(axis=0)[:, None] * P).sum(axis=0), dip[:2].sum(axis=(0, 1))
            print(f"{name}: {int(atom.sum())} atoms, repeat {c:.4f} A, k {a.info['kmesh']}; dipole per repeat of the charges {np.round(mu[0], 3).tolist()}, "
                  f"of the induced dipoles {np.round(mu[1], 3).tolist()} e A")
    print(f"written {out}")


def dft(args):
    from ase.io import read

    MF = _load("mlip_frames")
    src = Path(args.frames)
    for start in sorted(src.glob("chain_*.start.xyz")):
        name = start.name[:-len(".start.xyz")]
        cube = src / f"{name}.esp.cube"
        if cube.exists():
            continue
        a = read(start)
        work = src / "work" / name
        text = MF.pw_input(a, args.pseudo, work / "scratch", args.conv_thr).replace("  tstress = .true.", "  tstress = .false.")
        t0 = time.time()
        out = work / f"{name}.out"
        finished = lambda: out.exists() and all(k in out.read_text(errors="replace") for k in ("JOB DONE", "convergence has been achieved"))  # noqa: E731
        if not finished():
            MF.run_pw(args.pw, args.np, text, work, name)
        if not finished():
            print(f"{name}: the SCF did not finish ({out})", flush=True)
            continue
        (work / "pp.in").write_text("\n".join([
            "&INPUTPP", "  prefix = 'frame'", f"  outdir = '{work / 'scratch'}'", f"  filplot = '{work / 'esp.dat'}'", "  plot_num = 11", "/",
            "&PLOT", "  iflag = 3", "  output_format = 6", f"  fileout = '{cube}'", "/", ""]))
        env = {**os.environ, "OMP_NUM_THREADS": "1", "PATH": f"{Path(args.pw).parent}:{os.environ.get('PATH', '/usr/bin:/bin')}"}
        with open(work / "pp.out", "w") as fh:
            subprocess.run([str(Path(args.pw).with_name("pp.x")), "-in", "pp.in"], cwd=work, stdout=fh, stderr=subprocess.STDOUT, check=False, env=env)
        print(f"{name}: {'potential written' if cube.exists() else 'pp.x FAILED'}, {time.time() - t0:.0f} s", flush=True)


def samples(cube: Path, atoms_at, cell, r_min: float, r_max: float, step: float):
    """``(points, DFT potential in kcal/mol per e)``: grid points between ``r_min`` and ``r_max``
    from the nearest nucleus (its images along the chain counted), about ``step`` apart."""
    from ase.io.cube import read_cube_data

    data, _ = read_cube_data(str(cube))
    n = np.array(data.shape)
    stride = np.maximum(1, np.round(step / (cell / n)).astype(int))
    idx = [np.arange(0, n[k], stride[k]) for k in range(3)]
    grid = np.stack(np.meshgrid(*idx, indexing="ij"), axis=-1).reshape(-1, 3)
    pts = grid / n * cell
    near = np.full(len(pts), np.inf)
    for s in atoms_at:
        d = pts - s
        d -= cell * np.round(d / cell)
        near = np.minimum(near, np.sqrt((d * d).sum(axis=1)))
    keep = (near >= r_min) & (near <= r_max)
    # the code's potential is the energy of an electron in Ry: the potential a positive charge feels is its negative
    return pts[keep], -RY_TO_KCAL * data[grid[keep, 0], grid[keep, 1], grid[keep, 2]], near[keep]


EXTRA = IW.EXTRA
START = {"s_CH": 1.0, "s_CF": 1.0, "d_CC": 0.0, **{k: 0.0 for k in EXTRA}}


def sampled_blocks(src: Path, names, r_min: float, r_max: float, spacing: float, sites: bool = False) -> dict:
    """``{chain: sites, patterns, dipoles, cell, points, DFT potential, distance to the nearest
    nucleus}`` from the working directory's potential files."""
    out = {}
    for name in names:
        if not (src / f"chain_{name}.esp.cube").exists():
            continue
        m = np.load(src / f"chain_{name}.model{'_sites' if sites else ''}.npz")
        el, where, cell = [str(e) for e in m["elements"]], m["sites"], m["cell"]
        atom = np.array([e in SR.ATOMS for e in el])
        pts, v, near = samples(src / f"chain_{name}.esp.cube", where[atom], cell, r_min, r_max, spacing)
        out[name] = {"elements": el, "sites": where, "patterns": m["patterns"], "dipoles": m["dipoles"], "cell": cell, "points": pts, "v": v, "near": near}
    return out


def collect(args):
    """The sampled potentials in one compressed file the fit can be redone from (the grids
    themselves are 30 MB a chain)."""
    data = sampled_blocks(Path(args.frames), [c for c in args.chains.split(",") if c], args.r_min, args.r_max, args.spacing, args.sites)
    flat = {"chains": np.array(list(data)), "r_min": args.r_min, "r_max": args.r_max, "spacing": args.spacing}
    for name, b in data.items():
        flat.update({f"{name}_elements": np.array(b["elements"]), f"{name}_sites": b["sites"], f"{name}_patterns": b["patterns"], f"{name}_dipoles": b["dipoles"],
                     f"{name}_cell": b["cell"], f"{name}_points": b["points"].astype(np.float32), f"{name}_v": b["v"].astype(np.float32),
                     f"{name}_near": b["near"].astype(np.float32)})
    np.savez_compressed(args.out, **flat)
    print(f"written {args.out}: " + ", ".join(f"{n} {len(b['v'])} points" for n, b in data.items()))


def stored_blocks(path: Path, names, r_min: float) -> dict:
    """:func:`sampled_blocks` from :func:`collect`'s file, the points nearer than ``r_min`` dropped."""
    z = np.load(path)
    out = {}
    for name in [str(c) for c in z["chains"] if str(c) in names]:
        keep = z[f"{name}_near"] >= r_min - 1e-6
        out[name] = {"elements": [str(e) for e in z[f"{name}_elements"]], "sites": z[f"{name}_sites"], "patterns": z[f"{name}_patterns"],
                     "dipoles": z[f"{name}_dipoles"], "cell": z[f"{name}_cell"], "points": z[f"{name}_points"][keep].astype(float),
                     "v": z[f"{name}_v"][keep].astype(float), "near": z[f"{name}_near"][keep].astype(float)}
    return out


def potential_blocks(data: dict, extras: bool = True) -> dict:
    """Per chain: DFT's potential at the sampled points and every source's potential there at unit
    strength, each with its own mean over the points removed (the zero of a periodic potential is
    a convention).  Without ``extras`` only the chain's own charge patterns are summed."""
    blocks = {}
    for name, d in data.items():
        pts, sites, pat, dip, cell = d["points"], d["sites"], d["patterns"], d["dipoles"], d["cell"]
        own = pattern_names(len(pat))
        basis = {own[k]: source_potential(pts, sites, pat[k], dip[k], cell) for k in range(len(pat))}
        bare = sum(source_potential(pts, sites, pat[k], 0.0 * dip[k], cell) for k in range(2))
        for key, (pos, q) in (IW.extra_sources(d["elements"], sites).items() if extras else ()):
            if key not in basis:  # a chain built with bond sites has the site's potential as polyfind's own, induction in it
                basis[key] = COULOMB * ewald_potential(pts, pos, q, cell)
        blocks[name] = {"v": d["v"] - d["v"].mean(), "basis": {k: b - b.mean() for k, b in basis.items()}, "bare": bare - bare.mean(), "near": d["near"],
                        "points": len(pts)}
    return blocks


def solve(blocks: dict, group, free) -> dict:
    """Least squares of the ``free`` multipliers on the chains of ``group``, the rest held at
    :data:`START`; each chain weighs the same whatever its number of points, and the errors are
    absolute (kcal/mol per e), not relative to the chain's own potential."""
    held = {k: v for k, v in START.items() if k not in free}
    model = lambda b, w: sum(w[k] * b["basis"][k] for k in w if w[k] != 0.0)  # noqa: E731
    A = np.concatenate([np.stack([blocks[n]["basis"][k] for k in free], axis=1) / np.sqrt(blocks[n]["points"]) for n in group])
    y = np.concatenate([(blocks[n]["v"] - model(blocks[n], held)) / np.sqrt(blocks[n]["points"]) for n in group])
    return {**held, **dict(zip(free, (float(c) for c in np.linalg.lstsq(A, y, rcond=None)[0])))}


def rms_error(b: dict, w: dict) -> float:
    return float(np.sqrt(((sum(w[k] * b["basis"][k] for k in w if w[k] != 0.0) - b["v"]) ** 2).mean()))


def fit(args):
    names = [c for c in args.chains.split(",") if c]
    data = (stored_blocks(Path(args.samples), names, args.r_min) if args.samples
            else sampled_blocks(Path(args.frames), names, args.r_min, args.r_max, args.spacing, args.sites))
    blocks = potential_blocks(data)

    def model(b, w):
        return sum(w[k] * b["basis"][k] for k in w if w[k] != 0.0)

    def score(b, w):
        mdl = model(b, w)
        return float(np.sqrt(((mdl - b["v"]) ** 2).mean())), float((mdl * b["v"]).sum() / (b["v"] ** 2).sum())

    out = {"what": "the electrostatic potential around one chain alone in a box, polyfind's charges and induced dipoles against PBE (pp.x plot_num 11)",
           "potential": "kcal/mol per e, each chain's own mean over the sampled points removed", "preset": PRESET,
           "sampled": f"grid points {args.r_min} to {args.r_max} A from the nearest nucleus, about {args.spacing} A apart", "chains": {}, "fits": {}}
    for name, b in blocks.items():
        rms = float(np.sqrt((b["v"] ** 2).mean()))
        e, slope = score(b, START)
        eb = float(np.sqrt(((b["bare"] - b["v"]) ** 2).mean()))
        shells = {}
        for lo, hi in ((args.r_min, 2.6), (2.6, 3.2), (3.2, 4.0), (4.0, args.r_max)):
            k = (b["near"] >= lo) & (b["near"] < hi)
            if k.sum() > 10:
                shells[f"{lo:.1f}-{hi:.1f}"] = {"points": int(k.sum()), "rms_dft": float(np.sqrt((b["v"][k] ** 2).mean())),
                                                "rms_error": float(np.sqrt(((model(b, START)[k] - b["v"][k]) ** 2).mean()))}
        out["chains"][name] = {"what": CHAINS[name][1], "points": b["points"], "rms_dft": rms, "rms_error": e, "model_over_dft": slope,
                               "charges_alone_rms_error": eb, "by_distance": shells}
        print(f"{name} ({CHAINS[name][1]}): {b['points']} points, DFT rms {rms:.3f} kcal/mol per e; polyfind as it is: rms error {e:.3f}, model over DFT {slope:.2f} "
              f"(its charges without the induced dipoles: {eb:.3f})")
        print("     by distance from the nearest nucleus: " + "; ".join(f"{k} A: DFT rms {s['rms_dft']:.3f}, error {s['rms_error']:.3f}" for k, s in shells.items()))
    everyone = list(blocks)
    trials = [(f"{n} alone, the three charge patterns", [n], list(PATTERNS)) for n in blocks]
    trials += [("together, the three charge patterns", everyone, list(PATTERNS))]
    trials += [(f"together, charges held, {k}", everyone, [k]) for k in EXTRA]
    trials += [(f"together, the three charge patterns and {k}", everyone, list(PATTERNS) + [k]) for k in EXTRA]
    trials += [("together, the CH2 and CF2 groups' charges and M_H", everyone, ["s_CH", "s_CF", "M_H"]),
               ("together, the three charge patterns, M_H and Q_F", everyone, list(PATTERNS) + ["M_H", "Q_F"]),
               ("together, everything", everyone, list(PATTERNS) + list(EXTRA))]
    for label, group, free in trials:
        w = solve(blocks, group, free)
        scores = {n: score(b, w) for n, b in blocks.items()}
        out["fits"][label] = {"fitted_to": group, "free": free, "multipliers": w, "rms_error": {n: s[0] for n, s in scores.items()},
                              "model_over_dft": {n: s[1] for n, s in scores.items()}}
        print(f"{label}: " + ", ".join(f"{k} {w[k]:+.3f}" for k in free) + ";  rms error " + ", ".join(f"{n} {s[0]:.3f}" for n, s in scores.items()))
    if args.out:
        Path(args.out).write_text(json.dumps(out, indent=1))
        print(f"written {args.out}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="step", required=True)
    s = sub.add_parser("frames")
    s.add_argument("--out", required=True, help="working directory for the DFT, outside the repository")
    s.add_argument("--cells", default="/opt/dft/interchain/polyfind.json", help="examples/interchain_frames.py's table, for the packings' cells")
    s.add_argument("--chains", default="ht,hh,tgtg")
    s.add_argument("--sites", action="store_true", help="the chains with a charge site on every C-H bond (beside the plain ones, same atoms)")
    s = sub.add_parser("dft")
    s.add_argument("--frames", required=True)
    s.add_argument("--pw", required=True)
    s.add_argument("--pseudo", required=True)
    s.add_argument("--np", type=int, default=4)
    s.add_argument("--conv-thr", type=float, default=1e-9)
    s = sub.add_parser("collect")
    s.add_argument("--frames", required=True)
    s.add_argument("--chains", default="ht,hh,tgtg")
    s.add_argument("--r-min", type=float, default=2.2)
    s.add_argument("--r-max", type=float, default=5.0)
    s.add_argument("--spacing", type=float, default=0.3)
    s.add_argument("--sites", action="store_true", help="the chains built with bond sites (frames --sites)")
    s.add_argument("--out", default=str(PRED / "chain_potential_pbe_samples.npz"))
    s = sub.add_parser("fit")
    s.add_argument("--frames", default="", help="the working directory with the potential files")
    s.add_argument("--samples", default="", help="collect's file, instead of the working directory")
    s.add_argument("--sites", action="store_true", help="with --frames: the chains built with bond sites")
    s.add_argument("--chains", default="ht,hh,tgtg")
    s.add_argument("--r-min", type=float, default=2.2, help="A from the nearest nucleus: nearer than a neighbouring chain's atoms come")
    s.add_argument("--r-max", type=float, default=5.0)
    s.add_argument("--spacing", type=float, default=0.3, help="A between sampled grid points")
    s.add_argument("--out", default="")
    args = ap.parse_args()
    {"frames": frames, "dft": dft, "collect": collect, "fit": fit}[args.step](args)


if __name__ == "__main__":
    main()
