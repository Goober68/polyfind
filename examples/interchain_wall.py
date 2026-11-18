"""Which wall between chains reproduces DFT's energies of rigidly strained and sheared cells: a trial outside the kernels.

The cells of examples/interchain_frames.py and examples/shear_registry.py differ, within a packing,
only in where whole rigid chains are: the energy differences DFT gives for them are energy between
chains.  Polyfind's energy between chains is its electrostatics (the Ewald sum of the fixed charges
with the induced dipoles), which this trial leaves as it is, plus a sum over pairs of sites on
different chains of a wall-and-dispersion term, which is what is on trial.  With the pair
distances of every cell extracted once, any form of that term is a few array operations per cell,
and it can be fitted to DFT directly before a line of polyfind's kernels is changed.

``pairs --frames DIR [--registry DIR]``
    For every cell: the distances between sites of different chains inside the cutoff, by element
    pair, with polyfind's electrostatic energy and (where labelled) DFT's energy.  Checks that the
    present form on those distances is polyfind's own repulsion and dispersion to rounding.

``fit --form 12-6 | n-6 | exp-6 [--free ...]``
    Least squares of a form's parameters on DFT's differences minus polyfind's electrostatic ones.

The forms.  ``12-6``: polyfind's, ``D [(x/r)^12 - 2 (x/r)^6]`` per pair with geometric mixing, the
H-F pair's own term and the belt sites' wall against hydrogen.  ``n-6``: the same with the wall's
exponent a variable, ``D [6/(n-6) (x/r)^n - n/(n-6) (x/r)^6]`` (the minimum still at ``x`` with
depth ``D``).  ``exp-6``: ``D [6/(z-6) exp(z (1 - r/x)) - z/(z-6) (x/r)^6]``.  All force-shifted at
the cutoff as the packer's is.
"""
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np

from polyfind import fitting as F
from polyfind.fitting import PRESET_PARAMETERS


def _load(name):
    spec = importlib.util.spec_from_file_location(name, Path(__file__).with_name(f"{name}.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


IF = _load("interchain_frames")
SR = IF.SR
PRESET = SR.PRESET
PRED = Path("docs/predictions")
EV_TO_KCAL = SR.EV_TO_KCAL
ELEMENTS = ("C", "H", "F", "X")
PAIRS = tuple((a, b) for i, a in enumerate(ELEMENTS) for b in ELEMENTS[i:])


def interchain_distances(Pn, latn, n: int, rc: float) -> dict:
    """``{pair index (i, j) of the repeat's sites: distances}`` between sites of different chains
    inside ``rc``, every unordered pair of sites once.  ``Pn`` is the placed cell, chain after
    chain, ``n`` sites each; a chain's images across the chains count as different chains, its
    images along ``c`` (the third lattice vector) alone do not."""
    N = Pn.shape[0]
    chain = np.arange(N) // n
    site = np.arange(N) % n
    inv = np.linalg.inv(latn)
    widths = abs(np.linalg.det(latn)) / np.linalg.norm(np.cross(np.roll(latn, -1, axis=0), np.roll(latn, -2, axis=0)), axis=1)
    reach = np.ceil(rc / widths + 1.0).astype(int)
    frac = Pn @ inv
    out: dict = {}
    for i in range(-reach[0], reach[0] + 1):
        for j in range(-reach[1], reach[1] + 1):
            for k in range(-reach[2], reach[2] + 1):
                shift = np.array([i, j, k], dtype=float) @ latn
                d = Pn[:, None, :] - Pn[None, :, :] - shift
                r = np.sqrt((d * d).sum(axis=-1))
                other = (chain[:, None] != chain[None, :]) | (i != 0) | (j != 0)
                ia, ib = np.nonzero(other & (r < rc))
                for a, b, rr in zip(site[ia], site[ib], r[ia, ib]):
                    out.setdefault((int(a), int(b)), []).append(float(rr))
    # every unordered pair of sites was met twice (a against b at +shift, b against a at -shift)
    return {k: np.array(v) for k, v in out.items()}


def by_element(dist: dict, elements) -> dict:
    """The distances of :func:`interchain_distances` gathered by element pair (sorted as :data:`PAIRS`)."""
    order = {e: i for i, e in enumerate(ELEMENTS)}
    out = {p: [] for p in PAIRS}
    for (a, b), r in dist.items():
        ea, eb = elements[a], elements[b]
        key = (ea, eb) if order[ea] <= order[eb] else (eb, ea)
        out[key].append(r)
    return {p: (np.concatenate(v) if v else np.zeros(0)) for p, v in out.items()}


def wall(r, x, D, form: str, steep: float, rc: float, dispersion: bool = True):
    """The pair term of ``form`` at distances ``r``, minimum ``D`` deep at ``x``, force-shifted at
    ``rc``; ``steep`` is the wall's exponent ``n`` (``n-6``) or the exponential's ``z`` (``exp-6``).
    Without ``dispersion`` it is the wall alone (the belt sites')."""
    def raw(rr):
        s6 = (x / rr) ** 6
        if form == "12-6":
            rep, drep = D * s6 * s6, -12.0 * D * s6 * s6 / rr
            c6 = 2.0
        elif form == "n-6":
            n = steep
            rep, drep = D * 6.0 / (n - 6.0) * (x / rr) ** n, -n * D * 6.0 / (n - 6.0) * (x / rr) ** n / rr
            c6 = n / (n - 6.0)
        elif form == "exp-6":
            z = steep
            rep = D * 6.0 / (z - 6.0) * np.exp(z * (1.0 - rr / x))
            drep = -z / x * rep
            c6 = z / (z - 6.0)
        else:
            raise ValueError(f"unknown form {form!r}")
        if not dispersion:
            return rep, drep
        return rep - c6 * D * s6, drep + 6.0 * c6 * D * s6 / rr

    v, _ = raw(np.asarray(r, dtype=float))
    vc, dvc = raw(np.array(rc))
    return v - vc - (np.asarray(r, dtype=float) - rc) * dvc


def lj_tables(params: F.FFParameters) -> dict:
    """``{element pair: (x, D)}`` of a parameter set as the packer mixes it: geometric means, the
    pair overrides, and the belt site seeing hydrogen alone."""
    lj = {e: (float(x), float(d)) for e, x, d in params.lj}
    out = {}
    for a, b in PAIRS:
        out[(a, b)] = (np.sqrt(lj[a][0] * lj[b][0]), np.sqrt(lj[a][1] * lj[b][1]))
    for a, b, x, d in params.lj_pairs:
        key = (a, b) if (a, b) in out else (b, a)
        if key in out:
            out[key] = (float(x), float(d))
    return out


def pair_energy(cell: dict, tables: dict, form: str = "12-6", steep: float = 12.0, rc: float = 8.0, belt_form: str | None = None,
                belt_steep: float | None = None) -> float:
    """The energy between chains of one cell, kcal/mol per cell: half the sum over ordered pairs."""
    e = 0.0
    for p in PAIRS:
        r = cell["r"][p]
        x, D = tables[p]
        if r.size == 0 or D == 0.0:
            continue
        if "X" in p:  # the belt wall: repulsion only
            e += wall(r, x, D, belt_form or form, belt_steep or steep, rc, dispersion=False).sum()
        else:
            e += wall(r, x, D, form, steep, rc).sum()
    return 0.5 * float(e)


def pairs(args):
    from ase.io import read

    src = Path(args.frames)
    reg = Path(args.registry) if args.registry else None
    cells = []
    with PRESET_PARAMETERS[PRESET].applied():
        table = json.loads((src / "polyfind.json").read_text())
        alpha = table["alpha"]["cell"]
        tabs = lj_tables(PRESET_PARAMETERS[PRESET])
        worst = 0.0
        for pack in ("beta", "candidate", "alpha"):
            model, x0, ch, n_mon = IF.model_of(pack, alpha)
            rc = float(model.pk.rc)
            el = [str(e) for e in model.pk.elements]
            monomer = float(model.state(x0).lattice[3]) / (n_mon // 2)
            todo = [(row["frame"], (row["scale"][0], row["scale"][1], row["shear_monomers"][0], row["shear_monomers"][1]), src, row) for row in table[pack]["frames"]]
            if reg is not None and pack != "alpha":
                rt = json.loads((reg / "polyfind.json").read_text())[pack]["shears"]
                todo += [(row["frame"], (1.0, 1.0, row["shear_monomers"][0], row["shear_monomers"][1]), reg, row) for row in rt]
            first = None
            for name, v, where, row in todo:
                x = IF.moved(model, x0, ch, v, monomer)
                latn, _, _, Pn, _ = model.build(x)
                r = by_element(interchain_distances(Pn, latn, model.n, rc), el)
                e, es = model.physical_energy(x), SR.electrostatics(model, x)
                f = where / f"{name}.dft.xyz"
                dft = float(read(f).info["REF_energy"]) * EV_TO_KCAL if f.exists() else None
                cell = {"frame": name, "packing": pack, "variant": list(v), "monomers": n_mon, "rc": rc, "r": r, "polyfind": e, "electrostatics": es, "dft": dft}
                mine = pair_energy(cell, tabs, rc=rc)
                first = (e - es, mine) if first is None else first
                worst = max(worst, abs(((e - es) - first[0]) - (mine - first[1])))
                cells.append(cell)
            print(f"{pack}: {len(todo)} cells, {n_mon} monomers; site pairs inside the cutoff in the unstrained cell: "
                  + ", ".join(f"{a}-{b} {cells[-len(todo)]['r'][(a, b)].size}" for a, b in PAIRS if cells[-len(todo)]['r'][(a, b)].size))
    print(f"the present form on these distances against polyfind's own repulsion and dispersion, as differences within a packing: "
          f"largest disagreement {worst:.1e} kcal/mol per cell")
    out = Path(args.out)
    np.save(out, np.array(cells, dtype=object), allow_pickle=True)
    print(f"written {out}: {len(cells)} cells, {sum(c['dft'] is not None for c in cells)} with DFT")


VARIABLES = ("x_C", "x_H", "x_F", "D_C", "D_H", "D_F", "x_HF", "D_HF", "log10_D_X", "steep", "steep_H")
BOUNDS = {"x_C": (0.8, 1.4), "x_H": (0.6, 1.3), "x_F": (0.8, 1.3), "D_C": (0.2, 5.0), "D_H": (0.2, 5.0), "D_F": (0.2, 5.0),
          "x_HF": (0.7, 1.3), "D_HF": (0.2, 5.0), "log10_D_X": (-6.0, 0.5), "steep": (7.0, 16.0), "steep_H": (7.0, 16.0)}


def start_vector(form: str) -> dict:
    """The present preset in the trial's variables: multipliers of one on its diameters and depths
    (the H-F pair's on the pair's own values), the belt wall's depth, and the form's steepness
    (12 for a power, 13.8 for the exponential, the value at which exp-6 has the 12-6 curvature at
    its minimum).  ``steep_H`` is the steepness of the pairs that involve hydrogen or a belt site,
    when it is freed separately from the rest."""
    belt = dict((e, d) for e, _, d in PRESET_PARAMETERS[PRESET].lj)["X"]
    z0 = 13.772 if form == "exp-6" else 12.0
    return {"x_C": 1.0, "x_H": 1.0, "x_F": 1.0, "D_C": 1.0, "D_H": 1.0, "D_F": 1.0, "x_HF": 1.0, "D_HF": 1.0,
            "log10_D_X": float(np.log10(belt)), "steep": z0, "steep_H": z0}


def trial_tables(v: dict) -> dict:
    """``{element pair: (x, D)}`` of the trial variables ``v``, mixed as the packer mixes."""
    base = PRESET_PARAMETERS[PRESET]
    lj = {e: (float(x), float(d)) for e, x, d in base.lj}
    hf = next((x, d) for a, b, x, d in base.lj_pairs if {a, b} == {"F", "H"})
    hf_rule = (np.sqrt(lj["H"][0] * lj["F"][0]), np.sqrt(lj["H"][1] * lj["F"][1]))
    el = {"C": (lj["C"][0] * v["x_C"], lj["C"][1] * v["D_C"]), "H": (lj["H"][0] * v["x_H"], lj["H"][1] * v["D_H"]),
          "F": (lj["F"][0] * v["x_F"], lj["F"][1] * v["D_F"]), "X": (lj["X"][0], 10.0 ** v["log10_D_X"])}
    out = {(a, b): (np.sqrt(el[a][0] * el[b][0]), np.sqrt(el[a][1] * el[b][1])) for a, b in PAIRS}
    # the H-F pair keeps the ratio the preset gave it to the mixing rule, times the trial's own multipliers
    mixed = out[("H", "F")]
    out[("H", "F")] = (mixed[0] * hf[0] / hf_rule[0] * v["x_HF"], mixed[1] * hf[1] / hf_rule[1] * v["D_HF"])
    for p in (("C", "X"), ("F", "X"), ("X", "X")):
        out[p] = (out[p][0], 0.0)
    return out


def trial_energy(cell: dict, v: dict, form: str, split: bool) -> float:
    """The energy between chains of a cell under the trial variables, kcal/mol per cell; with
    ``split`` the pairs with a hydrogen or a belt site take ``steep_H`` and the rest ``steep``."""
    tabs = trial_tables(v)
    e = 0.0
    for p in PAIRS:
        r = cell["r"][p]
        x, D = tabs[p]
        if r.size == 0 or D == 0.0:
            continue
        steep = v["steep_H"] if (split and ("H" in p or "X" in p)) else v["steep"]
        e += wall(r, x, D, form, steep, cell["rc"], dispersion="X" not in p).sum()
    return 0.5 * float(e)


def groups(cells) -> list:
    """``[(label, reference cell, cells)]``: each packing's strained cells against its unstrained
    one, and its sheared cells (examples/shear_registry.py's) against theirs, labelled cells only."""
    out = []
    for pack in ("beta", "candidate", "alpha"):
        for kind, prefix in (("strain", "inter_"), ("shear", "shear_")):
            mine = [c for c in cells if c["packing"] == pack and c["frame"].startswith(prefix) and c["dft"] is not None]
            if len(mine) > 1:
                out.append((f"{pack} {kind}", mine[0], mine[1:]))
    return out


def residuals(v: dict, cells, form: str, split: bool, sigma=(0.02, 0.10)):
    """``[(group, frame, model, DFT, z)]`` over the labelled cells: differences against each
    group's reference in kcal/mol per monomer; ``z`` with a tolerance of ``sigma[0]`` and
    ``sigma[1]`` of DFT's value in quadrature."""
    out = []
    for label, ref, mine in groups(cells):
        e_ref = trial_energy(ref, v, form, split) + ref["electrostatics"]
        for c in mine:
            model = (trial_energy(c, v, form, split) + c["electrostatics"] - e_ref) / c["monomers"]
            dft = (c["dft"] - ref["dft"]) / c["monomers"]
            out.append((label, c["frame"], model, dft, (model - dft) / np.hypot(sigma[0], sigma[1] * dft)))
    return out


def summary(rows) -> dict:
    out = {}
    for label in dict.fromkeys(r[0] for r in rows):
        d = np.array([r[2] - r[3] for r in rows if r[0] == label])
        ref = np.array([r[3] for r in rows if r[0] == label])
        mod = np.array([r[2] for r in rows if r[0] == label])
        out[label] = {"cells": int(d.size), "rms_error": float(np.sqrt((d * d).mean())), "rms_dft": float(np.sqrt((ref * ref).mean())),
                      "slope_model_over_dft": float((mod * ref).sum() / (ref * ref).sum())}
    d = np.array([r[2] - r[3] for r in rows])
    out["all"] = {"cells": int(d.size), "rms_error": float(np.sqrt((d * d).mean())), "chi2": float(sum(r[4] ** 2 for r in rows))}
    return out


def fit(args):
    from scipy.optimize import least_squares

    cells = list(np.load(args.pairs, allow_pickle=True))
    form = args.form
    v0 = start_vector(form)
    free = [n for n in args.free.split(",") if n]
    for n in free:
        if n not in VARIABLES:
            raise SystemExit(f"unknown variable {n!r}: {VARIABLES}")
    split = "steep_H" in free
    for item in [s for s in args.set.split(",") if s]:
        k, val = item.split("=")
        v0[k] = float(val)

    def vec(y):
        v = dict(v0)
        v.update(dict(zip(free, y)))
        return v

    def fun(y):
        return np.array([r[4] for r in residuals(vec(y), cells, form, split)])

    start = summary(residuals(v0, cells, form, split))
    print(f"form {form}; {start['all']['cells']} labelled cells; before: rms error {start['all']['rms_error']:.4f} kcal/mol per monomer, chi2 {start['all']['chi2']:.1f}")
    if free:
        lo, hi = np.array([BOUNDS[n][0] for n in free]), np.array([BOUNDS[n][1] for n in free])
        y0 = np.clip(np.array([v0[n] for n in free]), lo + 1e-9, hi - 1e-9)
        res = least_squares(fun, y0, bounds=(lo, hi), x_scale=np.maximum(np.abs(y0), 0.1) * 0.05, xtol=1e-10, ftol=1e-10)
        v = vec(res.x)
        at_bound = [n for n, y, a, b in zip(free, res.x, lo, hi) if min(y - a, b - y) < 1e-6 * (b - a)]
    else:
        v, at_bound = dict(v0), []
    rows = residuals(v, cells, form, split)
    end = summary(rows)
    print(f"after: rms error {end['all']['rms_error']:.4f}, chi2 {end['all']['chi2']:.1f}; " + ", ".join(f"{n} {v[n]:.4f}" for n in free) + (f"  AT A BOUND: {at_bound}" if at_bound else ""))
    for label, sm in end.items():
        if label != "all":
            print(f"   {label:18s} {sm['cells']:3d} cells: rms error {sm['rms_error']:.4f} of rms {sm['rms_dft']:.4f}; model over DFT {sm['slope_model_over_dft']:.2f}"
                  f"   (before: {start[label]['rms_error']:.4f}, {start[label]['slope_model_over_dft']:.2f})")
    if args.rows:
        for r in rows:
            print(f"      {r[0]:18s} {r[1]:22s} model {r[2]:+.4f} DFT {r[3]:+.4f} z {r[4]:+.2f}")
    if args.out:
        out_path = Path(args.out)
        disk = json.loads(out_path.read_text()) if out_path.exists() else {"what": "trial walls between rigid chains fitted to DFT outside polyfind's kernels (examples/interchain_wall.py)",
                                                                           "energy": "kcal/mol per monomer", "preset": PRESET, "trials": {}}
        disk["trials"][args.name or f"{form} free {args.free or 'nothing'}"] = {
            "form": form, "free": free, "set": args.set, "variables": v, "at_bound": at_bound, "before": start, "after": end,
            "rows": [{"group": r[0], "frame": r[1], "model": r[2], "dft": r[3]} for r in rows]}
        out_path.write_text(json.dumps(disk, indent=1))
        print(f"written {out_path}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="step", required=True)
    s = sub.add_parser("pairs")
    s.add_argument("--frames", required=True, help="the working directory of examples/interchain_frames.py")
    s.add_argument("--registry", default="", help="the working directory of examples/shear_registry.py, for its sheared cells too")
    s.add_argument("--out", default="/opt/dft/interchain/pairs.npy")
    s = sub.add_parser("fit")
    s.add_argument("--pairs", default="/opt/dft/interchain/pairs.npy")
    s.add_argument("--form", default="12-6", choices=("12-6", "n-6", "exp-6"))
    s.add_argument("--free", default="", help=f"comma list of {VARIABLES}; none scores the start")
    s.add_argument("--set", default="", help="start values, name=value comma list")
    s.add_argument("--rows", action="store_true", help="print every cell")
    s.add_argument("--name", default="")
    s.add_argument("--out", default="")
    args = ap.parse_args()
    {"pairs": pairs, "fit": fit}[args.step](args)


if __name__ == "__main__":
    main()
