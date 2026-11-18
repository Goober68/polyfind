"""A lone chain pulled along its axis: its energy against its repeat, DFT and polyfind.

The helix of plain PVDF held at a longer repeat opens its gauche bonds until it is the all-trans
chain; held at a shorter one it closes them.  The energy along that line is what the inside of a
chain is, with nothing between chains in it, and it is the line a chain follows when the
material switches.  DFT's line is made a point at a time: the last labelled point is relaxed at
the next repeat in our own potential tilted to DFT's forces there (examples/mlip_train.py relax
--dft-gradient --repeat), and DFT is asked once for the result.

    python examples/chain_stretch.py scan --start LONE.dft.xyz --repeats 4.70,4.75 --tag s --model MODEL --out DIR --pw PW --pseudo PSEUDO
    python examples/chain_stretch.py collect --frames DIR --branch helix --names "START,up_*,down_*" --out docs/predictions/chain_stretch_pbe_d3bj.json
    python examples/chain_stretch.py polyfind --dft docs/predictions/chain_stretch_pbe_d3bj.json --fit name=FILE.json --out FILE.json
"""
import argparse
import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
EV = 23.060548


def _load(name):
    spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def backbone(atoms):
    """``(indices, unwrapped positions)`` of the carbons along the chain, one repeat and three more."""
    from ase.neighborlist import neighbor_list

    i, j, S = neighbor_list("ijS", atoms, {("C", "C"): 1.85})
    cell = atoms.cell.array
    nb = {}
    for p, q, s in zip(i, j, S):
        nb.setdefault(int(p), []).append((int(q), s))
    cur, off, order, pos = min(nb), np.zeros(3, dtype=int), [], []
    for _ in range(len(nb) + 3):
        order.append(cur)
        pos.append(atoms.positions[cur] + off @ cell)
        q, s = max(nb[cur], key=lambda t: (atoms.positions[t[0]] + (off + t[1]) @ cell)[2])  # up the chain
        cur, off = q, off + s
    return order, np.array(pos)


def dihedrals(atoms):
    """Backbone dihedrals in degrees, one per backbone bond of the repeat, in the order of :func:`backbone`."""
    _, P = backbone(atoms)
    out = []
    for k in range(len(P) - 3):
        b0, b1, b2 = P[k] - P[k + 1], P[k + 2] - P[k + 1], P[k + 3] - P[k + 2]
        u = b1 / np.linalg.norm(b1)
        v, w = b0 - (b0 @ u) * u, b2 - (b2 @ u) * u
        out.append(float(np.degrees(np.arctan2(np.cross(u, v) @ w, v @ w))))
    return out


def scan(args):
    """From the labelled start, one repeat after another: relax in the tilted potential, label, go on from there."""
    from ase.io import read

    args.out.mkdir(parents=True, exist_ok=True)
    last = args.start
    for c in [float(v) for v in args.repeats.split(",") if v]:
        frame = f"{args.tag}_{int(round(1000 * c)):05d}"
        done = args.out / f"{frame}.dft.xyz"
        if not done.exists():
            if not (args.out / f"{frame}.start.xyz").exists():
                subprocess.run([sys.executable, "-I", "-W", "ignore", str(HERE / "mlip_train.py"), "relax", "--model", str(args.model), "--name", args.tag,
                                "--dft-gradient", "--mask", "0 0 0 0 0 0", "--repeat", f"{c:.6f}", "--frame", frame, "--threads", str(args.threads),
                                "--fmax", str(args.fmax), "--cells", str(last), "--out", str(args.out)], check=True)
            subprocess.run([sys.executable, "-W", "ignore", str(HERE / "mlip_frames.py"), "label", "--frames", str(args.out), "--pw", args.pw, "--pseudo", args.pseudo,
                            "--np", str(args.np), "--names", frame] + (["--pools", str(args.pools)] if args.pools else []), check=True)
        if not done.exists():
            raise SystemExit(f"{frame}: no label; the scan stops here")
        a = read(done)
        n = a.info["n_monomers"]
        print(f"{frame}: repeat {a.cell[2, 2]:.4f} A, DFT {a.info['REF_energy'] / n:.5f} eV per monomer (predicted {a.info['dft_predicted_eV_per_monomer']:.5f}), "
              f"largest force {np.abs(a.arrays['REF_forces']).max():.3f} eV/A, dihedrals " + " ".join(f"{d:+.0f}" for d in dihedrals(a)[:2 * n]), flush=True)
        last = done


# A branch of a line, and where polyfind starts it: the chain of a packing (examples/interchain_frames.py),
# or a chain built at a conformer's nominal torsions (examples/chain_conformers.py: the letters of
# its repeat and the angles its bonds of each type take in turn).  A file holds one polymer's
# branches, on one scale.
BRANCHES = {"helix": {"pack": "alpha"}, "trans": {"pack": "beta"},
            "hh_tt": {"letters": "VI", "pattern": {}},
            "hh_special": {"letters": "VI", "pattern": {"CF2-CF2": (60.0,), "other": (-60.0,)}},
            "hh_plain": {"letters": "VI", "pattern": {"CF2-CH2": (60.0, -60.0)}}}


def collect(args):
    """The labelled points of one branch of the line, by repeat, merged into the output file: energy
    per monomer against the lowest point of all branches, the tension along the chain (from DFT's
    stress), what is left of the forces, the backbone dihedrals."""
    import fnmatch

    from ase.io import read

    out = json.loads(args.out.read_text()) if args.out.exists() else {
        "what": "a lone chain in a 16 A box held at a repeat, PBE-D3(BJ), relaxed at each (examples/chain_stretch.py): energy per monomer "
                "against the lowest point of all branches (kcal/mol), the tension dE/dc per monomer (kcal/mol/A, from DFT's stress), the forces left "
                "on the atoms (eV/A). Each branch is one conformer followed along its repeat; a file holds one polymer.", "branches": {}}
    points = []
    for f in sorted(args.frames.glob("*.dft.xyz")):
        if not any(fnmatch.fnmatch(f.name[:-8], pat) for pat in args.names.split(",")):
            continue
        a = read(f)
        n = int(a.info["n_monomers"])
        frc = a.arrays["REF_forces"]
        szz = float(np.asarray(a.info["REF_stress"]).reshape(3, 3)[2, 2])
        points.append({"frame": a.info["frame"], "repeat": float(a.cell[2, 2]), "energy_eV_per_monomer": float(a.info["REF_energy"]) / n,
                       "predicted_minus_dft": (float(a.info["dft_predicted_eV_per_monomer"]) - float(a.info["REF_energy"]) / n) * EV,
                       "largest_force": float(np.abs(frc).max()), "tension": szz * a.get_volume() / float(a.cell[2, 2]) / n * EV,
                       "dihedrals": dihedrals(a)[:2 * n]})
    points.sort(key=lambda q: q["repeat"])
    out["branches"][args.branch] = {"monomers": n, "start": BRANCHES[args.branch], "points": points}
    low = min(q["energy_eV_per_monomer"] for br in out["branches"].values() for q in br["points"])
    out["lowest_eV_per_monomer"] = low
    for name, br in out["branches"].items():
        for q in br["points"]:
            q["energy"] = (q["energy_eV_per_monomer"] - low) * EV
            if name == args.branch:
                print(f"{q['frame']:28s} repeat {q['repeat']:.4f} A: {q['energy']:+.3f} kcal/mol per monomer, tension {q['tension']:+.3f}, "
                      f"largest force {q['largest_force']:.3f} eV/A, predicted {q['predicted_minus_dft']:+.3f}, dihedrals " + " ".join(f"{d:+.0f}" for d in q["dihedrals"]))
    args.out.write_text(json.dumps(out, indent=1))
    print(f"written {args.out}")


def walk(one, y0, repeats, tilts: bool = True) -> list:
    """A one-chain model held at each repeat (A, the model's own repeat unit) in turn, outward from
    the one nearest ``y0``'s: its torsions and backbone angles free, then its pendant groups tilted.
    ``[{repeat, energy, torsions, tilt, ...}]`` in the order given, kcal/mol per chain repeat."""
    from polyfind import mixedcell as MC

    free = one.index("internal")
    repeats = [float(c) for c in repeats]
    first = int(np.argmin([abs(c - y0[3]) for c in repeats]))
    order = sorted(range(len(repeats)), key=lambda k: repeats[k])
    at = order.index(first)
    out = {}
    for leg in (order[at:], order[:at][::-1]):
        y = np.array(y0, dtype=float)
        for k in leg:
            y[3] = repeats[k]
            one.reset_multipliers()
            y, info = one.relax(y, free=free, closed=[0], maxiter=800, patience=60, patience_tol=1e-6)
            if info["stalled"] or info["gradient_max"] > 1e-3:
                y, info = one.relax_on_surface(y, free=free, closed=[0], gtol=1e-4, maxiter=300, patience=60)
            e = one.physical_energy(y)
            tilt = float(MC.tilt_relaxation(one, y)[0]) if tilts else 0.0
            parts = one.terms(y)
            out[k] = {"repeat": repeats[k], "energy": e + tilt, "tilt": tilt, "torsions": [float(v) for v in one.state(y).tors[0]],
                      "angles": [float(v) for v in one.state(y).angs[0]], "c_mismatch": float(np.abs(parts["c_mismatch"]).max()),
                      "torsion_energy": float(parts["torsion"]), "bend_energy": float(parts["bend"])}
    return [out[k] for k in range(len(repeats))]


def polyfind_line(SP, model, x, ch, repeats, tilts: bool = True, chain: int = 0) -> list:
    """:func:`walk` for one chain of a cell, alone on the 16 A lattice under the parameter set in force."""
    from polyfind import mixedcell as MC

    CS = _load("chain_split")
    pk1 = SP.pack_mod.CrystalPacker(ch, n_chains=1, **SP.DK.KW)
    st = model.state(x)
    one = MC.MixedModel(model.polymer, model.name, pk1, model.flips[chain:chain + 1], [ch], angle_stiffness=model.angle_stiffness)
    y0 = np.concatenate([[CS.BOX, CS.BOX, 90.0, st.lattice[3]], [0.0, 0.0, st.place[chain, 2], 0.0], st.tors[chain], st.angs[chain]])
    return walk(one, y0, repeats, tilts)


def polyfind_line_from_torsions(SP, polymer, torsions, repeats, tilts: bool = True):
    """:func:`walk` for a chain built at ``torsions`` (the model's bond order); ``(line, monomers in the model's repeat)``."""
    from polyfind import mixedcell as MC
    from polyfind.pack import periodic_chain_from_torsions

    CS = _load("chain_split")
    tors = np.asarray(torsions, dtype=float)
    ch = periodic_chain_from_torsions(polymer, "lone", tors)
    pk1 = SP.pack_mod.CrystalPacker(ch, n_chains=1, **SP.DK.KW)
    one = MC.MixedModel(polymer, "lone", pk1, np.array([False]), [ch], angle_stiffness=SP.DK.angle_stiffness())
    y0 = np.concatenate([[CS.BOX, CS.BOX, 90.0, ch.c], np.zeros(4), tors, one.angles0])
    return walk(one, y0, repeats, tilts), int(ch.n_monomers)


def polyfind_branches(SP, IF, table, branches: dict, tilts: bool = True) -> dict:
    """``{branch: [energy per monomer at each of its repeats]}`` under the parameter set in force,
    every branch on one scale (the chains are one polymer), with the gauche angle beside it."""
    out = {}
    for name, br in branches.items():
        spec = BRANCHES[name]
        if "pack" in spec:
            model, x0, ch, _ = IF.model_of(spec["pack"], table[spec["pack"]].get("cell"), check=False)
            n = ch.n_monomers
            line = polyfind_line(SP, model, x0, ch, [q["repeat"] * n / br["monomers"] for q in br["points"]], tilts=tilts)
        else:
            CC = _load("chain_conformers")
            polymer = SP.C.make_polymer(name, spec["letters"])
            tors = CC.torsions_of(CC.chain_types(polymer, 2 * len(spec["letters"])), spec["pattern"])
            line, n = polyfind_line_from_torsions(SP, polymer, tors, [q["repeat"] * len(spec["letters"]) / br["monomers"] for q in br["points"]], tilts=tilts)
        out[name] = [{"repeat": q["repeat"], "energy": r["energy"] / n, "tilt": r["tilt"] / n, "c_mismatch": r["c_mismatch"],
                      "gauche": float(np.abs(np.mod(np.array(r["torsions"]) + 180.0, 360.0) - 180.0).min()), "angles": r["angles"],
                      "torsions": [float(v) for v in (np.mod(np.array(r["torsions"]) + 180.0, 360.0) - 180.0)]}
                     for q, r in zip(br["points"], line)]
    return out


def polyfind(args):
    from polyfind import crystalfit as C
    from polyfind import fitting as F

    CA = _load("candidate_acceptance")
    SP = CA.SP
    IF = _load("interchain_frames")
    table = json.loads((HERE.parent / "docs" / "predictions" / "interchain_belt4.json").read_text())["packings"]
    dft = json.loads(args.dft.read_text())["branches"]
    flat = [(name, q) for name, br in dft.items() for q in br["points"]]
    ref = int(np.argmin([q["energy"] for _, q in flat]))
    out = json.loads(args.out.read_text()) if args.out and args.out.exists() else {"what": "polyfind's lone chains of plain PVDF held at DFT's repeats (examples/chain_stretch.py): "
                                                                                           "energy per monomer against the point where DFT is lowest, kcal/mol", "sets": {}}
    sets = {}
    for item in args.fit:
        name, path = item.split("=", 1)
        sets[name] = F.PRESET_PARAMETERS[path] if path in F.PRESET_PARAMETERS else C.parameters(np.array(json.loads(Path(path).read_text())["x"], dtype=float))
    for name, prm in sets.items():
        previous = SP.C.SITES_FOR
        SP.C.SITES_FOR = prm
        try:
            with prm.applied():
                lines = polyfind_branches(SP, IF, table, dft)
        finally:
            SP.C.SITES_FOR = previous
        mine = [r for br in dft for r in lines[br]]
        zero = mine[ref]["energy"]
        for r in mine:
            r["energy"] -= zero
        out["sets"][name] = lines
        err = np.array([r["energy"] - q["energy"] for r, (_, q) in zip(mine, flat)])
        print(f"{name}: rms error along the line {np.sqrt((err ** 2).mean()):.3f} kcal/mol per monomer")
        for r, (br, q) in zip(mine, flat):
            print(f"   {br:5s} repeat {q['repeat']:.3f}: polyfind {r['energy']:+.3f} (gauche {r['gauche']:.0f}), DFT {q['energy']:+.3f} (gauche {np.abs(q['dihedrals']).min():.0f})")
    if args.out:
        args.out.write_text(json.dumps(out, indent=1))
        print(f"written {args.out}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="step", required=True)
    s = sub.add_parser("scan")
    s.add_argument("--start", type=Path, required=True, help="a DFT-labelled lone chain (examples/chain_split.py, examples/mlip_frames.py label)")
    s.add_argument("--repeats", required=True, help="A, comma list, in the order they are visited")
    s.add_argument("--tag", required=True, help="the frames' prefix")
    s.add_argument("--model", type=Path, required=True)
    s.add_argument("--out", type=Path, required=True)
    s.add_argument("--pw", required=True)
    s.add_argument("--pseudo", required=True)
    s.add_argument("--np", type=int, default=4)
    s.add_argument("--pools", type=int, default=0)
    s.add_argument("--threads", type=int, default=1)
    s.add_argument("--fmax", type=float, default=0.005)
    s = sub.add_parser("collect")
    s.add_argument("--frames", type=Path, required=True)
    s.add_argument("--names", required=True, help="the branch's frames, comma list, shell patterns allowed (up_*)")
    s.add_argument("--branch", required=True, choices=sorted(BRANCHES))
    s.add_argument("--out", type=Path, required=True)
    s = sub.add_parser("polyfind")
    s.add_argument("--dft", type=Path, required=True, help="the collected DFT line")
    s.add_argument("--fit", action="append", default=[], help="name=crystal-stage result file, or name=preset")
    s.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()
    {"scan": scan, "collect": collect, "polyfind": polyfind}[args.step](args)


if __name__ == "__main__":
    main()
