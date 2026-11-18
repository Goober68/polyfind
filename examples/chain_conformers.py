"""Lone chains in named conformations: DFT's relaxed energy and torsions for each, and polyfind held at those torsions.

What a torsion term of a bond type is fitted to.  A reversed unit puts one CF2-CF2 and one
CH2-CH2 bond into a chain, and in a candidate's helix cells the two are always gauche together
or trans together, so the cells give one number for both.  The polymer in which every second
unit is reversed can turn them apart: with two monomers a repeat, both at once or neither; with
four, the CF2-CF2 bonds alone (T3G T3G') or the CH2-CH2 bonds alone.

DFT's chain is settled in rounds (with ``--fixed-repeat`` at the repeat it was built with, and
then walked along its repeat by examples/chain_stretch.py scan): relaxed in our own potential, labelled, then relaxed in the
potential tilted to DFT's forces there and labelled again (examples/mlip_train.py relax
--dft-gradient; no atom may move more than a set distance in a round, since a tilt has no end
along a soft turn).  Polyfind's chain is held at DFT's torsions with its backbone angles, its
repeat and its pendant groups free: nothing can hop to another conformer, and at DFT's minimum
an error in where polyfind's own minimum sits costs it at second order only.

    python examples/chain_conformers.py frames --out DIR
    python examples/chain_conformers.py settle --frames DIR --names hh2_tt,hh2_special --model MODEL --pw PW --pseudo PSEUDO [--rounds 2]
    python examples/chain_conformers.py collect --frames DIR --out docs/predictions/chain_conformers_pbe_d3bj.json
    python examples/chain_conformers.py polyfind --dft docs/predictions/chain_conformers_pbe_d3bj.json --fit name=FILE.json
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
BOX = 16.0
# chain: (the letters of its repeat, {conformer: {bond type: the angles its bonds of that type take in turn}}); a bond not named is trans
CHAINS = {
    "hh2": ("VI", {"tt": {}, "special": {"CF2-CF2": (60.0,), "other": (-60.0,)}, "plain": {"CF2-CH2": (60.0, -60.0)}}),
    "hh4": ("VIVI", {"ff": {"CF2-CF2": (60.0, -60.0)}, "hh": {"other": (60.0, -60.0)}}),
}
REFERENCE = "hh2_tt"   # the zero of the energies
TYPE_OF = {"CF2-CH2": "CF2-CH2", "CH2-CF2": "CF2-CH2", "CF2-CF2": "CF2-CF2", "CH2-CH2": "other"}


def _load(name):
    spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def torsions_of(types, pattern: dict) -> np.ndarray:
    """A conformer's torsions for a chain whose bonds have these types (names), from its pattern."""
    seen, out = {}, []
    for t in types:
        angles = pattern.get(t)
        if angles is None:
            out.append(180.0)
        else:
            out.append(angles[seen.get(t, 0) % len(angles)])
            seen[t] = seen.get(t, 0) + 1
    return np.array(out)


def model_order(dihedrals, bonds, types) -> np.ndarray:
    """DFT's dihedrals (walking up the chain from some atom, either way round) put on a model's
    bonds: the turn of the list, forward or reversed, whose bond types are the model's.  Chains
    that differ by where the walk began or by its direction are the same chain."""
    mine = [TYPE_OF[b] for b in bonds]
    d = np.asarray(dihedrals, dtype=float)
    n = len(types)
    if len(d) != n:
        raise ValueError(f"{len(d)} dihedrals for a chain of {n} bonds")
    for seq, val in ((mine, d), (mine[::-1], d[::-1])):
        for s in range(n):
            if all(seq[(k + s) % n] == types[k] for k in range(n)):
                return np.array([val[(k + s) % n] for k in range(n)])
    raise ValueError(f"no turn of the bonds {bonds} is {list(types)}")


def chain_types(polymer, n):
    from polyfind import fitting as F

    return [F.REF_TORSION_TYPES[k] for k in F.polymer_torsion_types(polymer, n)]


def frames(args):
    """Each conformer as polyfind builds it at its nominal torsions, alone in a 16 A box."""
    from ase import Atoms
    from ase.io import write

    from polyfind.fitting import PRESET_PARAMETERS
    from polyfind.pack import periodic_chain_from_torsions

    C = _load("candidate_screen")
    args.out.mkdir(parents=True, exist_ok=True)
    with PRESET_PARAMETERS["pvdf-dft-valence-belt4"].applied():
        for chain, (letters, confs) in CHAINS.items():
            polymer = C.make_polymer(chain, letters)
            types = chain_types(polymer, 2 * len(letters))
            for conf, pattern in confs.items():
                tors = torsions_of(types, pattern)
                ch = periodic_chain_from_torsions(polymer, conf, tors)
                keep = [k for k, e in enumerate(ch.elements) if e in ("C", "H", "F")]
                a = Atoms(symbols=[ch.elements[k] for k in keep], positions=ch.coords[keep] + [BOX / 2, BOX / 2, 0.0], cell=[BOX, BOX, ch.c], pbc=True)
                name = f"{chain}_{conf}"
                a.info.update(frame=name, family="chain", n_monomers=len(letters), kmesh=f"1 1 {max(2, int(round(20.0 / ch.c)))}",
                              source=f"polyfind's {letters} chain, {conf}: " + " ".join(f"{t}:{v:+.0f}" for t, v in zip(types, tors)))
                write(args.out / f"{name}.start.xyz", a, format="extxyz")
                print(f"{name}: {len(a)} atoms, repeat {ch.c:.4f} A, closes to {ch.rotation_error:.2f} deg; " + " ".join(f"{t}:{v:+.0f}" for t, v in zip(types, tors)))


def settle(args):
    """Round 0: relax the start in the potential and label it.  Each later round: relax the last
    labelled frame in the potential tilted to DFT's gradient there, the repeat free, and label."""
    from ase.io import read

    CS = _load("chain_split")
    for name in [n for n in args.names.split(",") if n]:
        last = args.frames / f"{name}.start.xyz"
        for r in range(args.rounds + 1):
            frame = f"{name}_r{r}"
            done = args.frames / f"{frame}.dft.xyz"
            if not done.exists():
                if not (args.frames / f"{frame}.start.xyz").exists():
                    subprocess.run([sys.executable, "-I", "-W", "ignore", str(HERE / "mlip_train.py"), "relax", "--model", str(args.model), "--name", "lone",
                                    "--mask", "0 0 0 0 0 0" if args.fixed_repeat else "0 0 1 0 0 0",
                                    "--frame", frame, "--threads", str(args.threads), "--fmax", str(args.fmax), "--cells", str(last), "--out", str(args.frames)]
                                   + (["--dft-gradient", "--max-move", str(args.max_move)] if r else []), check=True)
                subprocess.run([sys.executable, "-W", "ignore", str(HERE / "mlip_frames.py"), "label", "--frames", str(args.frames), "--pw", args.pw, "--pseudo", args.pseudo,
                                "--np", str(args.np), "--names", frame] + (["--pools", str(args.pools)] if args.pools else []), check=True)
            if not done.exists():
                raise SystemExit(f"{frame}: no label")
            a = read(done)
            n = a.info["n_monomers"]
            t = CS.backbone_torsions(a)
            print(f"{frame}: repeat {a.cell[2, 2]:.4f} A, DFT {a.info['REF_energy'] / n:.5f} eV per monomer"
                  + (f" (predicted {a.info['dft_predicted_eV_per_monomer']:.5f})" if "dft_predicted_eV_per_monomer" in a.info else "")
                  + f", largest force {np.abs(a.arrays['REF_forces']).max():.3f} eV/A; " + " ".join(f"{b}:{d:+.0f}" for b, d in zip(t["bonds"], t["dihedrals"])), flush=True)
            last = done


def collect(args):
    """Each conformer's last labelled round: energy per monomer against :data:`REFERENCE`, its torsions, what is left of the forces."""
    from ase.io import read

    CS = _load("chain_split")
    rows = {}
    for chain, (letters, confs) in CHAINS.items():
        for conf in confs:
            name = f"{chain}_{conf}"
            steps = sorted(args.frames.glob(f"{name}_r*.dft.xyz"))
            if not steps:
                continue
            a = read(steps[-1])
            n = int(a.info["n_monomers"])
            rows[name] = {"chain": chain, "letters": letters, "conformer": conf, "rounds": len(steps) - 1, "repeat": float(a.cell[2, 2]),
                          "energy_eV_per_monomer": float(a.info["REF_energy"]) / n, "largest_force": float(np.abs(a.arrays["REF_forces"]).max()),
                          "predicted_minus_dft": (float(a.info["dft_predicted_eV_per_monomer"]) - float(a.info["REF_energy"]) / n) * EV if "dft_predicted_eV_per_monomer" in a.info else None,
                          **CS.backbone_torsions(a)}
    zero = rows[REFERENCE]["energy_eV_per_monomer"]
    for name, r in rows.items():
        r["energy"] = (r["energy_eV_per_monomer"] - zero) * EV
        print(f"{name:14s} {r['energy']:+.3f} kcal/mol per monomer, repeat {r['repeat']:.4f} A, largest force {r['largest_force']:.3f} eV/A, {r['rounds']} rounds"
              + ("" if r["predicted_minus_dft"] is None else f", predicted {r['predicted_minus_dft']:+.3f}") + "; " + " ".join(f"{b}:{d:+.0f}" for b, d in zip(r["bonds"], r["dihedrals"])))
    out = {"what": "lone chains of the polymer with every second unit reversed, in a 16 A box, each conformer settled in DFT (PBE-D3(BJ); examples/chain_conformers.py): "
                   f"energy per monomer against {REFERENCE} (kcal/mol), backbone dihedrals (degrees) with the two groups each turns between", "conformers": rows}
    args.out.write_text(json.dumps(out, indent=1))
    print(f"written {args.out}")


def polyfind_at_torsions(SP, polymer, torsions, tilts: bool = True) -> dict:
    """A chain alone on the 16 A lattice under the parameter set in force, its torsions held at
    ``torsions`` (the model's bond order), its backbone angles and its repeat free, then its
    pendant groups tilted.  ``{energy (kcal/mol per chain repeat), repeat, angles, n_monomers}``."""
    from polyfind import mixedcell as MC
    from polyfind.pack import periodic_chain_from_torsions

    tors = np.asarray(torsions, dtype=float)
    ch = periodic_chain_from_torsions(polymer, "lone", tors)
    pk1 = SP.pack_mod.CrystalPacker(ch, n_chains=1, **SP.DK.KW)
    one = MC.MixedModel(polymer, "lone", pk1, np.array([False]), [ch], angle_stiffness=SP.DK.angle_stiffness())
    y = np.concatenate([[BOX, BOX, 90.0, ch.c], np.zeros(4), tors, one.angles0])
    free = np.concatenate([one.index("c"), one.index("angs")])
    y, info = one.relax(y, free=free, closed=[0], maxiter=800, patience=60, patience_tol=1e-6)
    if info["stalled"] or info["gradient_max"] > 1e-3:
        y, info = one.relax_on_surface(y, free=free, closed=[0], gtol=1e-4, maxiter=300, patience=60)
    e = one.physical_energy(y)
    tilt = float(MC.tilt_relaxation(one, y)[0]) if tilts else 0.0
    st = one.state(y)
    return {"energy": e + tilt, "tilt": tilt, "repeat": float(st.lattice[3]), "angles": [float(v) for v in st.angs[0]], "n_monomers": int(ch.n_monomers),
            "closure": float(np.abs(one.terms(y)["rotation_error"]).max())}


def polyfind_conformers(SP, C, rows: dict, tilts: bool = True) -> dict:
    """``{name: energy per monomer}`` of DFT's conformers under the parameter set in force (an
    ``applied()`` block), each chain held at DFT's torsions."""
    out, polymers = {}, {}
    for name, r in rows.items():
        letters = r["letters"]
        if letters not in polymers:
            polymers[letters] = C.make_polymer(r["chain"], letters)
        polymer = polymers[letters]
        tors = model_order(r["dihedrals"], r["bonds"], chain_types(polymer, 2 * len(letters)))
        q = polyfind_at_torsions(SP, polymer, tors, tilts=tilts)
        out[name] = {"energy": q["energy"] / q["n_monomers"], "repeat": q["repeat"], "closure": q["closure"], "tilt": q["tilt"] / q["n_monomers"]}
    return out


def polyfind(args):
    from polyfind import crystalfit as CF
    from polyfind import fitting as F

    SP = _load("shear_polish")
    rows = json.loads(args.dft.read_text())["conformers"]
    for item in args.fit:
        name, path = item.split("=", 1)
        prm = F.PRESET_PARAMETERS[path] if path in F.PRESET_PARAMETERS else CF.parameters(np.array(json.loads(Path(path).read_text())["x"], dtype=float))
        previous = SP.C.SITES_FOR
        SP.C.SITES_FOR = prm
        try:
            with prm.applied():
                mine = polyfind_conformers(SP, SP.C, rows)
        finally:
            SP.C.SITES_FOR = previous
        zero = mine[REFERENCE]["energy"]
        print(name)
        for k, r in rows.items():
            print(f"   {k:14s} polyfind {mine[k]['energy'] - zero:+.3f}, DFT {r['energy']:+.3f}; repeat {mine[k]['repeat']:.3f} against {r['repeat']:.3f} A (closes to {mine[k]['closure']:.1e} deg)")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="step", required=True)
    s = sub.add_parser("frames")
    s.add_argument("--out", type=Path, required=True)
    s = sub.add_parser("settle")
    s.add_argument("--frames", type=Path, required=True)
    s.add_argument("--names", required=True)
    s.add_argument("--model", type=Path, required=True)
    s.add_argument("--pw", required=True)
    s.add_argument("--pseudo", required=True)
    s.add_argument("--rounds", type=int, default=2, help="rounds in the tilted potential after the first relaxation")
    s.add_argument("--max-move", type=float, default=0.3)
    s.add_argument("--fixed-repeat", action="store_true", help="hold the chain at the start's repeat: a folded conformer that is not a minimum of our "
                                                               "potential opens out to the all-trans chain if its repeat is free")
    s.add_argument("--np", type=int, default=4)
    s.add_argument("--pools", type=int, default=0)
    s.add_argument("--threads", type=int, default=1)
    s.add_argument("--fmax", type=float, default=0.005)
    s = sub.add_parser("collect")
    s.add_argument("--frames", type=Path, required=True)
    s.add_argument("--out", type=Path, required=True)
    s = sub.add_parser("polyfind")
    s.add_argument("--dft", type=Path, required=True)
    s.add_argument("--fit", action="append", default=[], help="name=crystal-stage result file, or name=preset")
    args = ap.parse_args()
    {"frames": frames, "settle": settle, "collect": collect, "polyfind": polyfind}[args.step](args)


if __name__ == "__main__":
    main()
