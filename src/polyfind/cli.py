"""Command-line interface: ``polyfind <command> [options]``."""
from __future__ import annotations

import argparse
import sys

import numpy as np


def _add_common(p):
    p.add_argument("--polymer", default="pvdf", help="pvdf or pe")


def cmd_fit(args):
    from .forcefield import SimpleFF, fit_ris
    from .polymers import get_polymer

    rep = fit_ris(get_polymer(args.polymer), SimpleFF(eps_r=args.eps_r), step=args.step, n_monomers=args.monomers, third_order=args.third_order)
    print(rep.model.describe())
    print(f"{rep.n_evaluations} energy evaluations")
    if args.out:
        rep.model.save(args.out)
        print(f"saved {args.out}")


def _model(args):
    from .ris import RISModel
    from .forcefield import SimpleFF, fit_ris
    from .polymers import get_polymer

    if args.model:
        return RISModel.load(args.model)
    return fit_ris(get_polymer(args.polymer), SimpleFF(), step=10.0, n_monomers=6, third_order=True).model


def cmd_enumerate(args):
    from .enumerate import enumerate_periodic, table
    from .polymers import get_polymer

    model = _model(args)
    cands = enumerate_periodic(get_polymer(args.polymer), model, max_period=args.max_period, k_per_period=args.k)
    print(table(cands, top=args.top))


def cmd_pack(args):
    from .pack import pack, periodic_chain, to_cif
    from .polymers import get_polymer
    from .refine import refine_crystal

    model = _model(args)
    poly = get_polymer(args.polymer)
    seq = model.parse(args.sequence)
    chain = periodic_chain(poly, seq, model.states)
    print(f"chain {chain.name}: c = {chain.c:.3f} A, {chain.n_atoms} atoms per repeat")
    res = pack(chain, n_random=args.n_random, n_refine=args.n_refine, rng=np.random.default_rng(args.seed), verbose=True)
    for r in res:
        print("  " + r.row())
    best = res[0]
    if args.refine:
        rr = refine_crystal(poly, best, maxfev=args.maxfev)
        print(rr.summary())
        best = rr.result
    if args.cif:
        open(args.cif, "w").write(to_cif(best, title=chain.name))
        print(f"wrote {args.cif}")


def cmd_sample(args):
    from .amorphous import ensemble_stats, sample_ensemble
    from .polymers import get_polymer

    model = _model(args)
    poly = get_polymer(args.polymer)
    states, coords = sample_ensemble(poly, model, args.n_bonds, args.n_chains, args.temperature, rng=np.random.default_rng(args.seed))
    print(ensemble_stats(poly, model, states, coords, args.temperature).summary())


def cmd_pipeline(args):
    from .pipeline import PipelineConfig, run_pipeline

    cfg = PipelineConfig(polymer=args.polymer, max_period=args.max_period, top_k_pack=args.top, n_random=args.n_random, refine=not args.no_refine, temperature=args.temperature, seed=args.seed, workers=args.workers)
    if args.model:
        from .ris import RISModel

        cfg.ris_model = RISModel.load(args.model)
    res = run_pipeline(cfg, verbose=True)
    print()
    print(res.report())


def main(argv=None):
    ap = argparse.ArgumentParser(prog="polyfind", description="stable chain conformations and crystal packings of semi-crystalline polymers")
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("fit", help="fit RIS energies from dihedral scans of the built-in potential")
    _add_common(p)
    p.add_argument("--step", type=float, default=10.0)
    p.add_argument("--monomers", type=int, default=6)
    p.add_argument("--eps-r", type=float, default=1.0)
    p.add_argument("--third-order", action="store_true")
    p.add_argument("--out", default=None)
    p.set_defaults(func=cmd_fit)

    p = sub.add_parser("enumerate", help="rank periodic chain conformations")
    _add_common(p)
    p.add_argument("--model", default=None, help="RIS model json (else fit on the fly)")
    p.add_argument("--max-period", type=int, default=8)
    p.add_argument("-k", type=int, default=60)
    p.add_argument("--top", type=int, default=15)
    p.set_defaults(func=cmd_enumerate)

    p = sub.add_parser("pack", help="pack one periodic chain conformation into a 2-chain cell")
    _add_common(p)
    p.add_argument("sequence", help="e.g. TT, TG+TG-, TTTG+TTTG-")
    p.add_argument("--model", default=None)
    p.add_argument("--n-random", type=int, default=3000)
    p.add_argument("--n-refine", type=int, default=4)
    p.add_argument("--refine", action="store_true", help="refine torsions + cell after packing")
    p.add_argument("--maxfev", type=int, default=2000)
    p.add_argument("--cif", default=None)
    p.add_argument("--seed", type=int, default=0)
    p.set_defaults(func=cmd_pack)

    p = sub.add_parser("sample", help="exact Boltzmann sampling of an amorphous ensemble")
    _add_common(p)
    p.add_argument("--model", default=None)
    p.add_argument("--n-bonds", type=int, default=200)
    p.add_argument("--n-chains", type=int, default=2000)
    p.add_argument("--temperature", type=float, default=450.0)
    p.add_argument("--seed", type=int, default=0)
    p.set_defaults(func=cmd_sample)

    p = sub.add_parser("pipeline", help="run the whole funnel and print a report")
    _add_common(p)
    p.add_argument("--model", default=None)
    p.add_argument("--max-period", type=int, default=8)
    p.add_argument("--top", type=int, default=3)
    p.add_argument("--n-random", type=int, default=3000)
    p.add_argument("--no-refine", action="store_true")
    p.add_argument("--temperature", type=float, default=450.0)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--workers", type=int, default=None, help="processes for packing+refinement; default min(candidates, cpu_count-1); 1 = serial")
    p.set_defaults(func=cmd_pipeline)

    args = ap.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main(sys.argv[1:])
