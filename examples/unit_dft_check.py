"""One small row against DFT, start to finish: polyfind's margin under two presets, DFT's, and where they part.

A unit's parameters are tested by the margin of the smallest cell that holds the unit
(docs/MLIP_SECOND_OPINION.md 3.10): four monomers a chain, forty-eight atoms or so.  This runs
the steps for one such row, each skipped when its result is there, so the same command goes on
after a stop:

1. polyfind: the packing search, the polish (belt4), then every cell relaxed under each preset
   with the shear free (examples/candidate_screen.py, polish_margins.py, rows_under_preset.py);
2. the first preset's all-trans cell and helix cells written as atoms and relaxed in our
   potential (examples/row_cells_frames.py write, mlip_train.py relax);
3. DFT single points on the relaxed all-trans cell and on the helix cells the potential puts
   within ``--window`` of its lowest (examples/mlip_frames.py label);
4. a lone chain cut from the all-trans cell and from DFT's lowest helix cell, one DFT single
   point each (examples/chain_split.py frames), which splits the margin into what is inside a
   chain and what is between chains;
5. the tables (examples/row_cells_frames.py collect, chain_split.py dft and polyfind).

    python examples/unit_dft_check.py --row vdf2-trfe1-ctfe1-4 --work /opt/dft/unit4 --model docs/predictions/mlip/models/all_s1.model \\
        --pw /opt/mm/qe/bin/pw.x --pseudo /opt/pseudo

The row must be in examples/candidate_screen.py's SEQUENCES.  DFT wants the machine to itself:
the steps run one after another.
"""
import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PRED = HERE.parent / "docs" / "predictions"
CONFS = "TTTT,TGTG',TG'TG,GTG'T,G'TGT"
KCAL = 23.0605
PRESETS = {"chain1": "pvdf-dft-chain1", "belt4": "pvdf-dft-valence-belt4"}


def run(cmd, log, isolated=False):
    """One step as its own process (``-I`` for the ones that read frame files), its output appended to ``log``."""
    line = [sys.executable] + (["-I"] if isolated else []) + ["-W", "ignore"] + [str(c) for c in cmd]
    with open(log, "a") as fh:
        fh.write("\n$ " + " ".join(line[1:]) + "\n")
        fh.flush()
        done = subprocess.run(line, stdout=fh, stderr=subprocess.STDOUT)
    if done.returncode:
        raise SystemExit(f"step failed ({done.returncode}), see {log}: {' '.join(str(c) for c in cmd[:3])}")


def safe(conf: str) -> str:
    return conf.replace("'", "p")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--row", required=True)
    ap.add_argument("--work", type=Path, required=True, help="working directory for the DFT, outside the repository; the row gets its own folder in it")
    ap.add_argument("--model", type=Path, required=True)
    ap.add_argument("--pw", required=True)
    ap.add_argument("--pseudo", required=True)
    ap.add_argument("--kmesh", default="4 2 2")
    ap.add_argument("--window", type=float, default=0.15, help="kcal/mol per monomer: helix cells this close to the potential's lowest get a DFT single point too")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--np", type=int, default=4)
    ap.add_argument("--check", type=Path, default=None, help="a labelled directory whose anchor check is copied (saves recomputing it)")
    ap.add_argument("--polyfind-only", action="store_true", help="stop after step 1")
    args = ap.parse_args()
    from ase.io import read

    row, work = args.row, args.work / args.row
    work.mkdir(parents=True, exist_ok=True)
    if args.check and not (work / "check").exists():
        shutil.copytree(args.check, work / "check")
    log = work / "steps.log"
    tables = ["--polished", PRED / "polished_margins_belt4_dft_rows.json", "--references", PRED / "dilute_kink_belt4_dft_rows.json"]
    margins = {name: PRED / f"margins_{name}_dft_rows.json" for name in PRESETS}

    def has_margin(path):
        return path.exists() and "margin" in json.loads(path.read_text()).get("rows", {}).get(row, {})

    # 1. polyfind
    screen = PRED / "candidate_screen_belt4_trfe_small.json"
    if not (screen.exists() and json.loads(screen.read_text())["sequences"].get(row, {}).get("derived")):
        run([HERE / "candidate_screen.py", "--sequences", row, "--conformations", CONFS, "--out", screen], log)
    if not has_margin(tables[1]):
        run([HERE / "polish_margins.py", "--rows", row, "--tg", CONFS.split(",", 1)[1], "--references", tables[3], "--out", tables[1]], log)
    for name, preset in PRESETS.items():
        if not has_margin(margins[name]):
            run([HERE / "rows_under_preset.py", "--preset", preset, "--rows", row, *tables, "--out", margins[name]], log)
    mine = {name: json.loads(margins[name].read_text())["rows"][row] for name in PRESETS}
    print(f"{row}: polyfind " + ", ".join(f"{n} {e['margin']:+.3f} (upright {e['margin_upright']:+.3f}, helix {e['helix_sheared']})" for n, e in mine.items()), flush=True)
    if args.polyfind_only:
        return

    # 2. the first preset's cells as atoms, relaxed in the potential
    first = next(iter(PRESETS))
    confs = list(mine[first]["cells"])
    name_of = {c: f"row_{row}_{safe(c)}_{first}" for c in confs}
    if not all((work / f"{n}.start.xyz").exists() for n in name_of.values()):
        run([HERE / "row_cells_frames.py", "write", "--margins", margins[first], "--row", row, "--kmesh", args.kmesh, "--out", work], log)
    todo = [work / f"{n}.start.xyz" for n in name_of.values() if not (work / f"{n}_in_all_s1.start.xyz").exists()]
    if todo:
        run([HERE / "mlip_train.py", "relax", "--model", args.model, "--name", "all_s1", "--cells", *todo, "--out", work, "--threads", args.threads], log, isolated=True)
    own = {c: read(work / f"{n}_in_all_s1.start.xyz").info["potential_eV_per_monomer"] for c, n in name_of.items()}
    helix = {c: v for c, v in own.items() if c != "TTTT"}
    low = min(helix.values())
    chosen = [c for c, v in sorted(helix.items(), key=lambda kv: kv[1]) if KCAL * (v - low) <= args.window]
    # mirror placements relax to the same energy: one of each pair is enough
    kept = []
    for c in chosen:
        if not any(abs(KCAL * (helix[c] - helix[k])) < 2e-3 for k in kept):
            kept.append(c)
    print(f"{row}: in the potential, kcal/mol per monomer above its lowest helix cell: " + ", ".join(f"{c} {KCAL * (v - low):+.3f}" for c, v in own.items())
          + f"; DFT on TTTT and {kept}", flush=True)

    # 3. DFT on the relaxed cells
    relaxed = {c: f"{name_of[c]}_in_all_s1" for c in ["TTTT"] + kept}
    missing = [n for n in relaxed.values() if not (work / f"{n}.dft.xyz").exists()]
    if missing:
        run([HERE / "mlip_frames.py", "label", "--frames", work, "--pw", args.pw, "--pseudo", args.pseudo, "--np", args.np, "--names", ",".join(missing)], log, isolated=True)
    dft = {c: read(work / f"{n}.dft.xyz") for c, n in relaxed.items()}
    energy = {c: a.info["REF_energy"] / a.info["n_monomers"] for c, a in dft.items()}
    best = min(kept, key=lambda c: energy[c])
    print(f"{row}: DFT all-trans minus helix " + ", ".join(f"{c} {KCAL * (energy['TTTT'] - energy[c]):+.3f}" for c in kept)
          + f" kcal/mol per monomer; largest force left {max(abs(a.arrays['REF_forces']).max() for a in dft.values()):.2f} eV/A", flush=True)

    # 4. lone chains of the all-trans cell and of DFT's lowest helix cell
    split = work / "split"
    split.mkdir(exist_ok=True)
    if not (split / "check").exists() and (work / "check").exists():
        shutil.copytree(work / "check", split / "check")
    pair = [relaxed["TTTT"], relaxed[best]]
    if not all((split / f"{n}_chain0.start.xyz").exists() for n in pair):
        run([HERE / "chain_split.py", "frames", "--cells", *[work / f"{n}.dft.xyz" for n in pair], "--out", split], log, isolated=True)
    if not all((split / f"{n}_chain0.dft.xyz").exists() for n in pair):
        run([HERE / "mlip_frames.py", "label", "--frames", split, "--pw", args.pw, "--pseudo", args.pseudo, "--np", args.np, "--pools", 1], log, isolated=True)

    # 5. the tables
    pairs = [f"{'free' if c == best else safe(c)}={name_of['TTTT']},{name_of[c]}" for c in kept]
    run([HERE / "row_cells_frames.py", "collect", "--frames", work, "--row", row, "--pairs", *pairs, "--margins", *[f"{n}={p}" for n, p in margins.items()],
         "--out", PRED / f"row_cells_dft_{row}.json", "--extxyz", PRED / "mlip" / f"row_cells_{row}_frames_pbe_d3bj.extxyz"], log)
    run([HERE / "chain_split.py", "dft", "--frames", split, "--pairs", f"{pair[0]}={pair[1]}", "--out", PRED / f"chain_split_dft_{row}.json"], log, isolated=True)
    run([HERE / "chain_split.py", "polyfind", "--row", row, "--confs", f"{best},TTTT", *[f"--fit={n}={p}" for n, p in PRESETS.items()], *tables,
         "--out", PRED / "chain_split_polyfind_dft_rows.json"], log)
    d = json.loads((PRED / f"chain_split_dft_{row}.json").read_text())["pairs"]
    q = next(iter(d.values()))
    print(f"{row}: DFT {q['cell']:+.3f} = inside a chain {q['inside_a_chain']:+.3f} + between chains {q['between_chains']:+.3f} (helix {best})", flush=True)
    p = json.loads((PRED / "chain_split_polyfind_dft_rows.json").read_text())["sets"]
    for n in PRESETS:
        r = p[n][row]["pairs"][f"TTTT minus {best}"]
        print(f"{row}: {n} upright, the same helix: {r['cell']:+.3f} = {r['inside_a_chain']:+.3f} {r['between_chains']:+.3f}; off by {r['inside_a_chain'] - q['inside_a_chain']:+.3f} inside, "
              f"{r['between_chains'] - q['between_chains']:+.3f} between", flush=True)


if __name__ == "__main__":
    main()
