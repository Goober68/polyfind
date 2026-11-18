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

``fit --form 12-6 | n-6 | exp-6 | b14-7 [--free ...]``
    Least squares of a form's parameters on DFT's differences minus polyfind's electrostatic ones.

The forms.  ``12-6``: polyfind's, ``D [(x/r)^12 - 2 (x/r)^6]`` per pair with geometric mixing, the
H-F pair's own term and the belt sites' wall against hydrogen.  ``n-6``: the same with the wall's
exponent a variable, ``D [6/(n-6) (x/r)^n - n/(n-6) (x/r)^6]`` (the minimum still at ``x`` with
depth ``D``).  ``exp-6``: ``D [6/(z-6) exp(z (1 - r/x)) - z/(z-6) (x/r)^6]``.  ``b14-7``: Halgren's
buffered 14-7, ``D ((1+d)/(r/x+d))^7 ((1+g)/((r/x)^7+g) - 2)``, which stays finite at ``r = 0``
where exp-6 turns over.  All force-shifted at the cutoff as the packer's is.
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
COULOMB = 332.0637133  # kcal A / (mol e^2)
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


def charge_patterns(elements, q, n_chains: int) -> np.ndarray:
    """``(3, sites of the cell)``: the cell's charges split into the part on the CH2 groups, the
    part on the CF2 groups, and a third pattern polyfind's charge model does not have: one
    electron moved along every backbone bond between a CH2 carbon and a CF2 carbon, from the
    first to the second (the CF2 carbon +1, the CH2 carbon -1 per such bond).  ``elements`` and
    ``q`` are one chain's, backbone carbons in order, each followed by what it carries; the
    chain closes on its own image.

    A chain built with charge sites on its C-H bonds (``M``, polyfind.polymers.bond_site_atom)
    gives a fourth row: -1 on every site and +1 on its hydrogen.  What the sites carry in ``q`` is
    then taken out of the first row, which is the atoms' own charges as it is without sites."""
    el = [str(e) for e in elements]
    q = np.asarray(q, dtype=float)
    carbons = [i for i, e in enumerate(el) if e == "C"]
    ends = carbons[1:] + [len(el)]
    kind = ["F" if "F" in el[c + 1:e] else "H" for c, e in zip(carbons, ends)]
    sited = "M" in el
    out = np.zeros((4 if sited else 3, len(el)))
    for c, e, k in zip(carbons, ends, kind):
        out[0 if k == "H" else 1, c:e] = q[c:e]
    for k in range(len(carbons)):
        a, b = k, (k + 1) % len(carbons)
        if kind[a] != kind[b]:
            out[2, carbons[a if kind[a] == "F" else b]] += 1.0
            out[2, carbons[b if kind[a] == "F" else a]] -= 1.0
    if sited:
        for i in [k for k, e in enumerate(el) if e == "M"]:
            out[3, i], out[3, i - 1] = -1.0, 1.0
            out[0, i - 1] += out[0, i]   # the hydrogen as it is without its site
            out[0, i] = 0.0
    return np.tile(out, (1, n_chains))


def site_charge(elements, q) -> float:
    """What each hydrogen has given its bond's site (e, positive), zero for a chain without sites."""
    el = [str(e) for e in elements]
    return float(-np.mean([q[i] for i, e in enumerate(el) if e == "M"])) if "M" in el else 0.0


def electrostatic_form(model, x, patterns) -> np.ndarray:
    """``K`` with the cell's electrostatic energy (fixed charges and the dipoles they induce)
    equal to ``v K v`` for charges ``v @ patterns``: the energy is a homogeneous quadratic in the
    charges, so six evaluations give it exactly."""
    latn, _, _, Pn, _ = model.build(x)
    fl = model.flips.astype(float)
    e = lambda q: float(model.pk._polarize(Pn, q, latn, fl)[0])  # noqa: E731
    m = len(patterns)
    K = np.zeros((m, m))
    for i in range(m):
        K[i, i] = e(patterns[i])
    for i in range(m):
        for j in range(i + 1, m):
            K[i, j] = K[j, i] = 0.5 * (e(patterns[i] + patterns[j]) - K[i, i] - K[j, j])
    return K


EXTRA = ("X_belt", "Q_F", "mu_F", "Q_H", "mu_H", "M_H")
BOND_SITE = 0.65  # where the C-H bond's own charge site sits, as a fraction of the way from the carbon to the hydrogen


def extra_sources(el, sites, h: float = 0.25) -> dict:
    """Sources polyfind's lattice does not have, each ``(positions, charges)`` of unit strength:

    ``X_belt``  on every fluorine +2 and on each of its two belt sites -1: a quadrupole across the
                C-F bond in the F-C-F plane, on sites the chain already carries;
    ``Q_F``, ``Q_H``  a quadrupole along the bond on every fluorine, every hydrogen (+1 at ``h`` A
                either side, -2 on the nucleus: ``2 h^2`` e A^2);
    ``mu_F``, ``mu_H``  a dipole along the bond on every fluorine, every hydrogen, pointing away
                from the carbon (1 e A);
    ``M_H``     a charge site on every C-H bond, :data:`BOND_SITE` of the way from the carbon to the
                hydrogen, carrying -1 with +1 on the hydrogen: the one of these a code that sums
                over point sites can carry as it is.

    ``sites`` are whole chains, not wrapped into the cell: each atom's carbon is the nearest one."""
    el = np.asarray([str(e) for e in el])
    sites = np.asarray(sites, dtype=float)
    C = np.nonzero(el == "C")[0]
    X = np.nonzero(el == "X")[0]
    out = {k: ([], []) for k in EXTRA}

    def add(key, pos, q):
        out[key][0].extend(pos)
        out[key][1].extend(q)

    for i in np.nonzero((el == "F") | (el == "H"))[0]:
        c = C[np.argmin(np.linalg.norm(sites[C] - sites[i], axis=1))]
        u = (sites[i] - sites[c]) / np.linalg.norm(sites[i] - sites[c])
        e = str(el[i])
        add(f"Q_{e}", [sites[i] + h * u, sites[i], sites[i] - h * u], [1.0, -2.0, 1.0])
        add(f"mu_{e}", [sites[i] + h * u, sites[i] - h * u], [0.5 / h, -0.5 / h])
        if e == "H":
            add("M_H", [sites[c] + BOND_SITE * (sites[i] - sites[c]), sites[i]], [-1.0, 1.0])
        if e == "F" and X.size >= 2:
            near = X[np.argsort(np.linalg.norm(sites[X] - sites[i], axis=1))[:2]]
            add("X_belt", [sites[i], sites[near[0]], sites[near[1]]], [2.0, -1.0, -1.0])
    return {k: (np.array(v[0]).reshape(-1, 3), np.array(v[1])) for k, v in out.items()}


def ewald_green(pos, lattice, alpha: float = 0.45, tol: float = 1e-10) -> np.ndarray:
    """``G`` (1/A) with the electrostatic energy of unit point charges at ``pos`` repeated on
    ``lattice`` (rows) equal to ``(1/2) q G q`` for a neutral set ``q``: every pair's interaction
    over all images, a charge's with its own images on the diagonal, no term at the zero wave
    vector.  A plain Ewald sum for a cell of a few hundred sources."""
    pos = np.asarray(pos, dtype=float)
    lat = np.asarray(lattice, dtype=float)
    vol = abs(np.linalg.det(lat))
    widths = vol / np.linalg.norm(np.cross(np.roll(lat, -1, axis=0), np.roll(lat, -2, axis=0)), axis=1)
    cut = np.sqrt(-np.log(tol))
    rc = cut / alpha
    frac = pos @ np.linalg.inv(lat)
    span = np.ceil(frac.max(axis=0) - frac.min(axis=0)).astype(int)  # whole chains are not wrapped: pairs can be cells apart
    reach = np.ceil(rc / widths).astype(int) + span
    from scipy.special import erfc

    d0 = pos[:, None, :] - pos[None, :, :]
    G = np.zeros((len(pos), len(pos)))
    same = np.eye(len(pos), dtype=bool)
    for i in range(-reach[0], reach[0] + 1):
        for j in range(-reach[1], reach[1] + 1):
            for k in range(-reach[2], reach[2] + 1):
                d = d0 + np.array([i, j, k], dtype=float) @ lat
                r = np.sqrt((d * d).sum(axis=-1))
                ok = (r < rc) & ~(same & (i == 0) & (j == 0) & (k == 0)) & (r > 1e-12)
                G += np.where(ok, erfc(alpha * r) / np.where(ok, r, 1.0), 0.0)
    rec = 2.0 * np.pi * np.linalg.inv(lat).T
    gmax = 2.0 * alpha * cut
    m = np.ceil(gmax / (2.0 * np.pi) * np.linalg.norm(lat, axis=1)).astype(int)
    idx = np.array([[i, j, k] for i in range(-m[0], m[0] + 1) for j in range(-m[1], m[1] + 1) for k in range(-m[2], m[2] + 1)
                    if (i, j, k) != (0, 0, 0)], dtype=float)
    g = idx @ rec
    g2 = (g * g).sum(axis=1)
    keep = g2 < gmax * gmax
    g, g2 = g[keep], g2[keep]
    ph = np.exp(1j * (pos @ g.T))
    G += ((ph * (4.0 * np.pi / vol * np.exp(-g2 / (4.0 * alpha * alpha)) / g2)) @ ph.conj().T).real
    G[same] -= 2.0 * alpha / np.sqrt(np.pi)
    return G


def first_order_form(el, Pn, latn, patterns) -> np.ndarray:
    """The electrostatic energy (kcal/mol per cell) of fixed point sources, nothing induced, as a
    quadratic form over the three charge patterns of :func:`charge_patterns` followed by the
    sources of :func:`extra_sources` at unit strength.  A dipole or quadrupole is a cluster of
    charges, whose own energy is a constant of the chain and leaves every difference between
    rigid cells."""
    extra = extra_sources(el, Pn)
    pos = [Pn] + [extra[k][0] for k in EXTRA]
    start = np.cumsum([0] + [len(p) for p in pos])
    allpos = np.concatenate(pos)
    Q = np.zeros((3 + len(EXTRA), len(allpos)))
    Q[:3, :len(Pn)] = patterns[:3]
    for k, name in enumerate(EXTRA):
        Q[3 + k, start[k + 1]:start[k + 2]] = extra[name][1]
    return 0.5 * COULOMB * Q @ ewald_green(allpos, latn) @ Q.T


def by_element(dist: dict, elements) -> dict:
    """The distances of :func:`interchain_distances` gathered by element pair (sorted as :data:`PAIRS`)."""
    order = {e: i for i, e in enumerate(ELEMENTS)}
    out = {p: [] for p in PAIRS}
    for (a, b), r in dist.items():
        ea, eb = elements[a], elements[b]
        if ea not in order or eb not in order:  # a bond's charge site: no pair term with anything
            continue
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
        elif form == "b14-7":  # Halgren's buffered 14-7, finite at r = 0; ``steep`` is its delta (0.07), gamma 0.12
            dl, gm = steep, 0.12
            def b147(q):
                rho = q / x
                return D * ((1.0 + dl) / (rho + dl)) ** 7 * ((1.0 + gm) / (rho ** 7 + gm) - (2.0 if dispersion else 0.0))
            h = 1e-6
            return b147(rr), (b147(rr + h) - b147(rr - h)) / (2.0 * h)
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
    params = PRESET_PARAMETERS[PRESET]
    if args.charges:
        # the same rigid cells under other charges: "dCH=..,dCF=..,M=.." are the two bond increments and the
        # charge of a site on every C-H bond (the chains are then built with those sites)
        from dataclasses import replace

        new = dict(item.split("=") for item in args.charges.split(","))
        incr = tuple((a, b, float(new.get(f"d{a}{b}", d))) for a, b, d in params.charge_increments)
        params = replace(params, charge_increments=incr, bond_sites=((("H", float(new["M"])),) if "M" in new else ()),
                         backbone_pull=((("F", float(new["pull_F"])),) if "pull_F" in new else ()))
        SR.C.SITES_FOR = params
        print(params.describe())
    with params.applied():
        table = json.loads((src / "polyfind.json").read_text())
        tabs = lj_tables(params)
        worst = worst_es = worst_1 = 0.0
        for pack in table:
            model, x0, ch, n_mon = IF.model_of(pack, table[pack].get("cell"), check=not args.charges)
            rc = float(model.pk.rc)
            el = [str(e) for e in model.pk.elements]
            monomer = float(model.state(x0).lattice[3]) / (n_mon // 2)
            todo = [(row["frame"], (row["scale"][0], row["scale"][1], row["shear_monomers"][0], row["shear_monomers"][1]), src, row) for row in table[pack]["frames"]]
            if reg is not None and pack in ("beta", "candidate"):
                rt = json.loads((reg / "polyfind.json").read_text())[pack]["shears"]
                todo += [(row["frame"], (1.0, 1.0, row["shear_monomers"][0], row["shear_monomers"][1]), reg, row) for row in rt]
            first = None
            pat = charge_patterns(el, np.asarray(model.pk._q_cell, dtype=float)[:model.n], model.nc)
            q_site = site_charge(el, np.asarray(model.pk._q_cell, dtype=float)[:model.n])
            for name, v, where, row in todo:
                x = IF.moved(model, x0, ch, v, monomer)
                latn, _, _, Pn, _ = model.build(x)
                r = by_element(interchain_distances(Pn, latn, model.n, rc), el)
                e, es = model.physical_energy(x), SR.electrostatics(model, x)
                K = electrostatic_form(model, x, pat)
                own = np.array([1.0, 1.0, 0.0, q_site][:len(pat)])
                worst_es = max(worst_es, abs(own @ K @ own - es))
                K1 = first_order_form(list(model.pk.elements) * model.nc, Pn, latn, pat)
                idx = model.pk._es_idx
                q_all = np.asarray(model.pk._q_cell, dtype=float)
                fixed = float(model.pk._ewald.terms(Pn if idx is None else Pn[idx], q_all if idx is None else q_all[idx], latn).total)
                if not q_site:  # the plain sum has the atoms' charges alone
                    first_fixed = (fixed, K1[:2, :2].sum()) if first is None else first_fixed
                    worst_1 = max(worst_1, abs((fixed - first_fixed[0]) - (K1[:2, :2].sum() - first_fixed[1])))
                f = where / f"{name}.dft.xyz"
                dft = float(read(f).info["REF_energy"]) * EV_TO_KCAL if f.exists() else None
                cell = {"frame": name, "packing": pack, "variant": list(v), "monomers": n_mon, "rc": rc, "r": r, "polyfind": e, "electrostatics": es, "dft": dft, "K": K, "K1": K1}
                mine = pair_energy(cell, tabs, rc=rc)
                first = (e - es, mine) if first is None else first
                worst = max(worst, abs(((e - es) - first[0]) - (mine - first[1])))
                cells.append(cell)
            print(f"{pack}: {len(todo)} cells, {n_mon} monomers; site pairs inside the cutoff in the unstrained cell: "
                  + ", ".join(f"{a}-{b} {cells[-len(todo)]['r'][(a, b)].size}" for a, b in PAIRS if cells[-len(todo)]['r'][(a, b)].size))
    print(f"the present form on these distances against polyfind's own repulsion and dispersion, as differences within a packing: "
          f"largest disagreement {worst:.1e} kcal/mol per cell")
    print(f"the quadratic form of the charges at polyfind's own charges against its electrostatic energy: largest disagreement {worst_es:.1e} kcal/mol per cell")
    print(f"a plain Ewald sum of polyfind's fixed charges against polyfind's own, as differences within a packing: largest disagreement {worst_1:.1e} kcal/mol per cell")
    out = Path(args.out)
    np.save(out, np.array(cells, dtype=object), allow_pickle=True)
    print(f"written {out}: {len(cells)} cells, {sum(c['dft'] is not None for c in cells)} with DFT")


VARIABLES = ("x_C", "x_H", "x_F", "D_C", "D_H", "D_F", "x_HF", "D_HF", "log10_D_X", "steep", "steep_H", "x_FF", "D_FF", "x_HH", "D_HH", "es", "s_CH", "s_CF", "d_CC") + EXTRA
BOUNDS = {"x_C": (0.8, 1.4), "x_H": (0.6, 1.3), "x_F": (0.8, 1.3), "D_C": (0.2, 5.0), "D_H": (0.2, 5.0), "D_F": (0.2, 5.0),
          "x_HF": (0.7, 1.3), "D_HF": (0.2, 5.0), "log10_D_X": (-6.0, 0.5), "steep": (7.0, 16.0), "steep_H": (7.0, 16.0),
          "x_FF": (0.7, 1.3), "D_FF": (0.2, 5.0), "x_HH": (0.7, 1.3), "D_HH": (0.2, 5.0), "es": (0.0, 3.0),
          "s_CH": (0.0, 2.0), "s_CF": (0.0, 2.0), "d_CC": (-0.4, 0.4), **{k: (-3.0, 3.0) for k in EXTRA}}


def start_vector(form: str) -> dict:
    """The present preset in the trial's variables: multipliers of one on its diameters and depths
    (the H-F pair's on the pair's own values), the belt wall's depth, and the form's steepness
    (12 for a power, 13.8 for the exponential, the value at which exp-6 has the 12-6 curvature at
    its minimum).  ``steep_H`` is the steepness of the pairs that involve hydrogen or a belt site,
    when it is freed separately from the rest.  ``es`` multiplies polyfind's electrostatic energy:
    a diagnostic, not a parameter of polyfind (its charges and polarisabilities are not on trial
    here), which says whether what a wall cannot fit has the electrostatics' shape.  ``s_CH`` and
    ``s_CF`` multiply the charges of the CH2 and the CF2 groups, and ``d_CC`` is a charge moved
    along each backbone bond from a CH2 carbon to a CF2 carbon (:func:`charge_patterns`), in e."""
    belt = dict((e, d) for e, _, d in PRESET_PARAMETERS[PRESET].lj)["X"]
    z0 = {"exp-6": 13.772, "b14-7": 0.07}.get(form, 12.0)
    return {"x_C": 1.0, "x_H": 1.0, "x_F": 1.0, "D_C": 1.0, "D_H": 1.0, "D_F": 1.0, "x_HF": 1.0, "D_HF": 1.0,
            "log10_D_X": float(np.log10(belt)), "steep": z0, "steep_H": z0, "x_FF": 1.0, "D_FF": 1.0, "x_HH": 1.0, "D_HH": 1.0, "es": 1.0,
            "s_CH": 1.0, "s_CF": 1.0, "d_CC": 0.0, **{k: 0.0 for k in EXTRA}}


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
    # like pairs may leave the element's own values too (multipliers on them; one changes nothing)
    out[("F", "F")] = (out[("F", "F")][0] * v.get("x_FF", 1.0), out[("F", "F")][1] * v.get("D_FF", 1.0))
    out[("H", "H")] = (out[("H", "H")][0] * v.get("x_HH", 1.0), out[("H", "H")][1] * v.get("D_HH", 1.0))
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
    for pack in dict.fromkeys(c["packing"] for c in cells):
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
        w = np.array([v.get("s_CH", 1.0), v.get("s_CF", 1.0), v.get("d_CC", 0.0)])
        u = np.array([v.get(k, 0.0) for k in EXTRA])

        def es(c):
            # polyfind's own energy of the charge patterns (induction in it), and the sources it does not have
            # to first order: with the charges and with each other, nothing induced by them
            if c["K"].shape[0] == 4:  # cells built with the bond sites: their charge is polyfind's own fourth pattern
                w4 = np.append(w, v.get("M_H", 0.0))
                ux = np.where(np.array(EXTRA) == "M_H", 0.0, u)
                more = 2.0 * float(w @ c["K1"][:3, 3:] @ ux) + float(ux @ c["K1"][3:, 3:] @ ux) if np.any(ux) else 0.0
                return v.get("es", 1.0) * (float(w4 @ c["K"] @ w4) + more)
            more = 2.0 * float(w @ c["K1"][:3, 3:] @ u) + float(u @ c["K1"][3:, 3:] @ u) if np.any(u) else 0.0
            return v.get("es", 1.0) * (float(w @ c["K"] @ w) + more)

        e_ref = trial_energy(ref, v, form, split) + es(ref)
        for c in mine:
            model = (trial_energy(c, v, form, split) + es(c) - e_ref) / c["monomers"]
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
    if args.packings:
        cells = [c for c in cells if c["packing"] in args.packings.split(",")]
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
        box = dict(BOUNDS, **({"steep": (0.0, 0.6), "steep_H": (0.0, 0.6)} if form == "b14-7" else {}))
        lo, hi = np.array([box[n][0] for n in free]), np.array([box[n][1] for n in free])
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
    predicted = {}
    for pack in [p for p in args.predict.split(",") if p]:
        # the fitted variables on a packing's cells whether DFT has them or not: against its first cell, kcal/mol per monomer
        mine = [c for c in cells if c["packing"] == pack and c["frame"].startswith("inter_")]
        both = residuals(v, [dict(c, dft=0.0) for c in mine], form, split)
        predicted[pack] = {r[1]: r[2] for r in both}
        print(f"   predicted, {pack}: " + ", ".join(f"{c['variant'][2]:+.3f},{c['variant'][3]:+.3f}: {predicted[pack][c['frame']]:+.4f}"
                                                   + (f" (DFT {(c['dft'] - mine[0]['dft']) / c['monomers']:+.4f})" if c["dft"] is not None and mine[0]["dft"] is not None else "")
                                                   for c in mine[1:]))
    if args.out:
        out_path = Path(args.out)
        disk = json.loads(out_path.read_text()) if out_path.exists() else {"what": "trial walls between rigid chains fitted to DFT outside polyfind's kernels (examples/interchain_wall.py)",
                                                                           "energy": "kcal/mol per monomer", "preset": PRESET, "trials": {}}
        disk["trials"][args.name or f"{form} free {args.free or 'nothing'}"] = {
            "form": form, "free": free, "set": args.set, "variables": v, "at_bound": at_bound, "before": start, "after": end, "predicted": predicted,
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
    s.add_argument("--charges", default="", help="other charges for the same cells: dCH=..,dCF=..,M=..,pull_F=.. (a site on every C-H bond; what each F draws along the backbone)")
    s = sub.add_parser("fit")
    s.add_argument("--pairs", default="/opt/dft/interchain/pairs.npy")
    s.add_argument("--form", default="12-6", choices=("12-6", "n-6", "exp-6", "b14-7"))
    s.add_argument("--free", default="", help=f"comma list of {VARIABLES}; none scores the start")
    s.add_argument("--set", default="", help="start values, name=value comma list")
    s.add_argument("--packings", default="", help="comma list: fit and score these packings' cells alone")
    s.add_argument("--rows", action="store_true", help="print every cell")
    s.add_argument("--predict", default="", help="comma list of packings: the fitted variables on all their cells, labelled or not")
    s.add_argument("--name", default="")
    s.add_argument("--out", default="")
    args = ap.parse_args()
    {"pairs": pairs, "fit": fit}[args.step](args)


if __name__ == "__main__":
    main()
