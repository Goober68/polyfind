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
    python examples/chain_potential.py fit --frames /opt/dft/chainpot
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


def isolated(pack: str, table: dict):
    """One chain of a packing with the others far away: ``(elements, positions, patterns (3, n),
    induced dipoles (3, n, 3) in e A, repeat along the chain in A)`` of its sites, the dipoles
    those each charge pattern induces along the chain itself."""
    model, x0, ch, _ = IF.model_of(pack, table[pack].get("cell"))
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
    with PRESET_PARAMETERS[PRESET].applied():
        for name in [c for c in args.chains.split(",") if c]:
            el, P, pat, dip, c = isolated(CHAINS[name][0], table)
            atom = np.array([e in SR.ATOMS for e in el])
            centre = P[atom, :2].mean(axis=0)
            P = P - np.array([centre[0] - BOX / 2.0, centre[1] - BOX / 2.0, 0.0])
            a = Atoms(symbols=[e for e, k in zip(el, atom) if k], positions=P[atom], cell=[BOX, BOX, c], pbc=True)
            a.info.update(frame=f"chain_{name}", kmesh=f"1 1 {max(2, int(round(20.0 / c)))}", source=f"polyfind {PRESET}, one chain: {CHAINS[name][1]}")
            write(out / f"chain_{name}.start.xyz", a, format="extxyz")
            np.savez(out / f"chain_{name}.model.npz", elements=np.array(el), sites=P, patterns=pat, dipoles=dip, cell=np.array([BOX, BOX, c]))
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
        out, _ = MF.run_pw(args.pw, args.np, text, work, name)
        if MF.read_pw(out) is None:
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


def fit(args):
    src = Path(args.frames)
    names = [c for c in args.chains.split(",") if c and (src / f"chain_{c}.esp.cube").exists()]
    blocks = {}
    for name in names:
        m = np.load(src / f"chain_{name}.model.npz")
        el, sites, pat, dip, cell = [str(e) for e in m["elements"]], m["sites"], m["patterns"], m["dipoles"], m["cell"]
        atom = np.array([e in SR.ATOMS for e in el])
        pts, v, near = samples(src / f"chain_{name}.esp.cube", sites[atom], cell, args.r_min, args.r_max, args.step)
        basis = np.stack([source_potential(pts, sites, pat[k], dip[k], cell) for k in range(3)], axis=1)
        bare = np.stack([source_potential(pts, sites, pat[k], 0.0 * dip[k], cell) for k in range(3)], axis=1)
        blocks[name] = {"v": v - v.mean(), "basis": basis - basis.mean(axis=0), "bare": bare - bare.mean(axis=0), "near": near, "points": len(pts)}

    def score(b, w, key="basis"):
        model = b[key] @ w
        return float(np.sqrt(((model - b["v"]) ** 2).mean())), float((model * b["v"]).sum() / (b["v"] ** 2).sum())

    def solve(bs, free):
        A = np.concatenate([b["basis"][:, free] for b in bs])
        rest = [k for k in range(3) if k not in free]
        y = np.concatenate([b["v"] - b["basis"][:, rest] @ np.array([1.0, 1.0, 0.0])[rest] for b in bs])
        w = np.array([1.0, 1.0, 0.0])
        w[free] = np.linalg.lstsq(A, y, rcond=None)[0]
        return w

    out = {"what": "the electrostatic potential around one chain alone in a box, polyfind's charges and induced dipoles against PBE (pp.x plot_num 11)",
           "potential": "kcal/mol per e, each chain's own mean over the sampled points removed", "preset": PRESET,
           "sampled": f"grid points {args.r_min} to {args.r_max} A from the nearest nucleus, about {args.step} A apart", "chains": {}, "fits": {}}
    w0 = np.array([1.0, 1.0, 0.0])
    for name, b in blocks.items():
        rms = float(np.sqrt((b["v"] ** 2).mean()))
        e, slope = score(b, w0)
        eb, slope_b = score(b, w0, "bare")
        shells = {}
        for lo, hi in ((args.r_min, 2.6), (2.6, 3.2), (3.2, 4.0), (4.0, args.r_max)):
            k = (b["near"] >= lo) & (b["near"] < hi)
            if k.sum() > 10:
                mdl = b["basis"][k] @ w0
                shells[f"{lo:.1f}-{hi:.1f}"] = {"points": int(k.sum()), "rms_dft": float(np.sqrt((b["v"][k] ** 2).mean())),
                                                "model_over_dft": float((mdl * b["v"][k]).sum() / (b["v"][k] ** 2).sum())}
        out["chains"][name] = {"what": CHAINS[name][1], "points": b["points"], "rms_dft": rms, "rms_error": e, "model_over_dft": slope,
                               "charges_alone_rms_error": eb, "charges_alone_over_dft": slope_b, "by_distance": shells}
        print(f"{name} ({CHAINS[name][1]}): {b['points']} points, DFT rms {rms:.3f} kcal/mol per e; polyfind as it is: rms error {e:.3f}, model over DFT {slope:.2f} "
              f"(its charges without the induced dipoles: {eb:.3f}, {slope_b:.2f})")
        print("     by distance from the nearest nucleus: " + "; ".join(f"{k} A: DFT rms {s['rms_dft']:.3f}, model over DFT {s['model_over_dft']:.2f}" for k, s in shells.items()))
    trials = [("each chain alone, the two groups' charges", None, [0, 1]), ("each chain alone, all three", None, [0, 1, 2])]
    trials += [("all chains together, the two groups' charges", list(blocks), [0, 1]), ("all chains together, the charge along C-C alone", list(blocks), [2]),
               ("all chains together, all three", list(blocks), [0, 1, 2])]
    if "ht" in blocks:
        trials += [("fitted to plain all-trans PVDF alone, all three", ["ht"], [0, 1, 2])]
    for label, which, free in trials:
        for group in ([[n] for n in blocks] if which is None else [which]):
            w = solve([blocks[n] for n in group], free)
            scores = {n: score(b, w) for n, b in blocks.items()}
            key = f"{label}: {','.join(group)}"
            out["fits"][key] = {"fitted_to": group, "free": [PATTERNS[k] for k in free], **{PATTERNS[k]: float(w[k]) for k in range(3)},
                                "rms_error": {n: s[0] for n, s in scores.items()}, "model_over_dft": {n: s[1] for n, s in scores.items()}}
            print(f"{key}: " + ", ".join(f"{PATTERNS[k]} {w[k]:+.3f}" for k in range(3)) + "; rms error "
                  + ", ".join(f"{n} {s[0]:.3f} (x{s[1]:.2f})" for n, s in scores.items()))
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
    s = sub.add_parser("dft")
    s.add_argument("--frames", required=True)
    s.add_argument("--pw", required=True)
    s.add_argument("--pseudo", required=True)
    s.add_argument("--np", type=int, default=4)
    s.add_argument("--conv-thr", type=float, default=1e-9)
    s = sub.add_parser("fit")
    s.add_argument("--frames", required=True)
    s.add_argument("--chains", default="ht,hh,tgtg")
    s.add_argument("--r-min", type=float, default=2.2, help="A from the nearest nucleus: nearer than a neighbouring chain's atoms come")
    s.add_argument("--r-max", type=float, default=5.0)
    s.add_argument("--step", type=float, default=0.3)
    s.add_argument("--out", default="")
    args = ap.parse_args()
    {"frames": frames, "dft": dft, "fit": fit}[args.step](args)


if __name__ == "__main__":
    main()
