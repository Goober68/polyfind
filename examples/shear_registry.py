"""What it costs to move neighbouring chains along each other: polyfind against DFT on the same atoms.

Three steps, each writing or reading files; the DFT between them is ``examples/mlip_frames.py label``.

``frames --out DIR``
    Plain PVDF, two chains of one monomer, in two packings: polyfind's beta cell (the stored
    ``vdf10`` all-trans cell) and the packing of the reversed-unit candidate's all-trans cell
    (``vdf5-rvdf1-6``: other setting angles, another cell angle) with plain chains in it.  Each is
    sheared along its chains rigidly (``polyfind.mixedcell.MixedModel(shear=True)``): along ``a``,
    along ``b``, and along the line the candidate's cell took in DFT, to that cell's own shear.
    Writes one start frame per shear, polyfind's energy of each (split into its electrostatics,
    fixed charges and induced dipoles, and its repulsion and dispersion), and the seven registered E5
    cells (docs/predictions/e5_structures: one sublattice, or one chain of four, slid rigidly in
    Dipole's beta cell) as start frames too.

``candidate``
    The DFT-relaxed all-trans cell of the reversed-unit candidate put into the mixed cell's
    variables (least squares on the atoms, which the frame keeps in the exported order), and
    polyfind asked about it: its energy there, with the shear held and everything else relaxed,
    and with the shear free.

``collect --frames DIR --out JSON``
    The labelled frames against polyfind's numbers: a table per packing, the E5 cells against the
    registered like-for-like predictions, and the labelled frames gathered into one file.

    python examples/shear_registry.py frames --out /opt/dft/shear
    python examples/mlip_frames.py label --frames /opt/dft/shear --pw .../pw.x --pseudo .../pseudo
    python examples/shear_registry.py candidate
    python examples/shear_registry.py collect --frames /opt/dft/shear
"""
import argparse
import importlib.util
import json
import time
from pathlib import Path

import numpy as np

from polyfind import mixedcell as MC
from polyfind import pack as pack_mod
from polyfind.fitting import PRESET_PARAMETERS
from polyfind.linegroup import repeat_chains

DS = importlib.util.spec_from_file_location("dk", Path(__file__).with_name("dilute_kink.py"))
DK = importlib.util.module_from_spec(DS)
DS.loader.exec_module(DK)
C = DK.C
PRED = Path("docs/predictions")
EV_TO_KCAL = 23.060548
PRESET = "pvdf-dft-valence-belt4"
CANDIDATE = "vdf5-rvdf1-6"
DFT_SHEAR = (0.303, -0.909)  # the candidate's DFT-relaxed all-trans cell, in monomers along a and b
SHEARS = ([(0.0, 0.0)] + [(k / 8, 0.0) for k in (1, 2, 3, 4)] + [(0.0, k / 8) for k in (1, 2, 3, 4)]
          + [(0.1, -0.3), (0.2, -0.6), DFT_SHEAR, (0.4, -1.2), (0.5, -1.5)])
PACKINGS = {"beta": ("vdf10", "polyfind's beta cell"), "candidate": (CANDIDATE, "the packing of the reversed-unit candidate's all-trans cell, plain chains in it")}
E5 = {"e5_beta_a_slide_c8": ("a", "c/8"), "e5_beta_a_slide_c4": ("a", "c/4"), "e5_beta_a_slide_3c8": ("a", "3c/8"), "e5_beta_a_slide_c2": ("a", "c/2"),
      "e5_beta_2a_slide_c4": ("c", "c/4"), "e5_beta_2a_slide_c2": ("c", "c/2")}
ATOMS = ("C", "H", "F", "Cl", "N", "O")


def plain_model(pack: str):
    """``(model, x, chain, description)``: plain PVDF, two chains of one monomer, in the packing of
    the stored all-trans cell of ``PACKINGS[pack]``: its lateral cell, setting angles and heights
    (as fractions of a monomer), a one-monomer all-trans chain with that cell's mean backbone angles."""
    ref = json.load(open(PRED / "dilute_kink_belt4.json"))["all_trans_references"][PACKINGS[pack][0]]
    polymer = C.make_polymer("vdf1", "V")
    P = polymer.bonds_per_repeat
    n_ref = len(ref["tors"]) // P
    st_ref = MC.MixedState.from_vector(ref["state"], 2, len(ref["tors"]))
    cm = float(st_ref.lattice[3]) / n_ref
    angs = np.asarray(st_ref.angs[0]).reshape(-1, P).mean(axis=0)
    ch = repeat_chains(polymer, "TT", np.full((1, P), 180.0), angs)[0]
    pk = pack_mod.CrystalPacker(ch, n_chains=2, **DK.KW)
    place = st_ref.place.copy()
    place[:, 3] = (np.mod(place[:, 3] / cm + 0.5, 1.0) - 0.5) * ch.c
    st = MC.MixedState(np.array([st_ref.lattice[0], st_ref.lattice[1], st_ref.lattice[2], ch.c]), place,
                       np.full((2, P), 180.0), np.tile(angs, (2, 1)), np.zeros(2))
    model = MC.MixedModel(polymer, "TT", pk, np.asarray(ref["flips"]) > 0.5, [ch, ch], angle_stiffness=DK.angle_stiffness(), shear=True)
    return model, st.vector(), ch


def electrostatics(model, x) -> float:
    """The part of the cell's energy that moves with the chains' registry and is electrostatic:
    the Ewald sum of the fixed charges with the induced dipoles' energy (kcal/mol per cell).  The
    bonded-exclusion corrections belong to each chain alone and do not change in a rigid shear."""
    latn, _, _, Pn, _ = model.build(x)
    return float(model.pk._polarize(Pn, np.asarray(model.pk._q_cell, dtype=float), latn, model.flips.astype(float))[0])


def real_atoms(model):
    el = np.array([str(e) for e in list(model.pk.elements) * model.nc])
    return el, np.array([e in ATOMS for e in el])


def frames(args):
    from ase import Atoms
    from ase.io import read, write

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    table = {}
    with PRESET_PARAMETERS[PRESET].applied():
        for pack in PACKINGS:
            model, x0, ch = plain_model(pack)
            e0 = model.physical_energy(x0)
            es0 = electrostatics(model, x0)
            el, keep = real_atoms(model)
            st = model.state(x0)
            table[pack] = {"what": PACKINGS[pack][1], "cell_A_deg": st.lattice.tolist(), "setting_angles_deg": st.place[:, 2].tolist(),
                           "second_chain_flipped": bool(model.flips[1]), "polyfind_kcal_per_monomer": e0 / 2, "shears": []}
            for i, (fa, fb) in enumerate(SHEARS):
                x = x0.copy()
                x[model.index("shear")] = [fa * ch.c, fb * ch.c]
                e = model.physical_energy(x)
                latn, _, _, Pn, _ = model.build(x)
                a = Atoms(symbols=list(el[keep]), positions=Pn[keep], cell=latn, pbc=True)
                a.wrap()
                name = f"shear_{pack}_{i:02d}"
                a.info.update(family="shear", frame=name, n_monomers=2, kmesh="4 2 8", source=f"polyfind {PRESET}, plain PVDF, {PACKINGS[pack][1]}, sheared rigidly",
                              shear_monomers_a=fa, shear_monomers_b=fb)
                write(out / f"{name}.start.xyz", a, format="extxyz")
                es = electrostatics(model, x)
                table[pack]["shears"].append({"frame": name, "shear_monomers": [fa, fb], "polyfind": (e - e0) / 2, "polyfind_electrostatics": (es - es0) / 2,
                                              "polyfind_repulsion_dispersion": ((e - es) - (e0 - es0)) / 2})
            print(f"{pack}: cell {np.round(st.lattice, 4).tolist()}, {len(SHEARS)} frames; polyfind at half a monomer along a {table[pack]['shears'][4]['polyfind']:+.3f}, "
                  f"along b {table[pack]['shears'][8]['polyfind']:+.3f} kcal/mol per monomer")
    for f in sorted((PRED / "e5_structures").glob("*.xyz")):
        a = read(f)
        name = "e5_" + f.stem
        a.info = {"family": "e5", "frame": name, "n_monomers": len(a) // 6, "kmesh": "4 2 8" if len(a) == 12 else "2 2 8", "source": str(f)}
        write(out / f"{name}.start.xyz", a, format="extxyz")
    (out / "polyfind.json").write_text(json.dumps(table, indent=1))
    print(f"written {out}: {2 * len(SHEARS)} sheared cells and the seven E5 cells; label them with examples/mlip_frames.py label --frames {out}")


def upright_frame(a):
    """``(lattice rows, positions)`` of an ASE frame turned so that its third lattice vector is
    along z and the first one's part across it along x: the mixed cell's frame."""
    cell, pos = np.array(a.cell), a.get_positions()
    ez = cell[2] / np.linalg.norm(cell[2])
    ax = cell[0] - (cell[0] @ ez) * ez
    ex = ax / np.linalg.norm(ax)
    R = np.array([ex, np.cross(ez, ex), ez])
    return cell @ R.T, pos @ R.T


def candidate(args):
    from ase.io import read
    from scipy.optimize import least_squares

    out_path = Path(args.out)
    a = [f for f in read(PRED / "mlip" / "candidate_relaxation_frames_pbe_d3bj.extxyz", ":") if f.info["frame"].startswith("rvdf6_TT")][-1]
    latn, pos = upright_frame(a)
    A, B = float(np.hypot(*latn[0, :2])), float(np.hypot(*latn[1, :2]))
    gam, c, sa, sb = float(np.degrees(np.arctan2(latn[1, 1], latn[1, 0]))), float(latn[2, 2]), float(latn[0, 2]), float(latn[1, 2])
    inv = np.linalg.inv(latn)
    t0 = time.time()
    with PRESET_PARAMETERS[PRESET].applied():
        stored = json.load(open(PRED / "dilute_kink_belt4.json"))["all_trans_references"][CANDIDATE]
        d = {k: (np.asarray(v, dtype=float) if isinstance(v, list) else v) for k, v in stored.items()}
        letters, _ = DK.sequence(CANDIDATE)
        polymer = C.make_polymer(CANDIDATE, letters)
        ch = DK.reference_chain(polymer, d)
        pk = pack_mod.CrystalPacker(ch, n_chains=2, **DK.KW)
        st = MC.MixedState.from_vector(d["state"], 2, len(d["tors"]))
        st.shear = np.zeros(2)
        model = MC.MixedModel(polymer, "ref", pk, np.asarray(d["flips"]) > 0.5, [ch, ch], angle_stiffness=DK.angle_stiffness(), shear=True)
        x0 = st.vector()
        per = 2 * ch.n_monomers
        cm = float(st.lattice[3]) / ch.n_monomers
        e0, g0 = model.value_and_grad(x0)
        e0 = model.physical_energy(x0)
        el, keep = real_atoms(model)
        if list(el[keep]) != list(a.get_chemical_symbols()):
            raise RuntimeError("the DFT frame's atoms are not in the exported cell's order")

        def miss(v, base, free):
            x = base.copy()
            x[free] = v
            df = (model.build(x)[3][keep] - pos) @ inv
            return ((df - np.round(df)) @ latn).ravel()

        xs = x0.copy()
        xs[:4] = [A, B, gam, c]
        xs[model.index("shear")] = [sa, sb]
        place = np.concatenate([model.index("phi"), model.index("xy"), model.index("z")])
        best = None
        for dphi in (0.0, 180.0):  # the relaxed cell may have turned and moved: rigid chains first, from a spread of starts
            for dz in np.arange(0.0, c, c / 12):
                xt = xs.copy()
                xt[model.index("phi")] += dphi
                xt[model.index("z", 1)] += dz
                r = least_squares(miss, xt[place], args=(xt, place), method="lm")
                if best is None or r.cost < best[0]:
                    best = (r.cost, xt, r.x)
        xt = best[1].copy()
        xt[place] = best[2]
        every = np.concatenate([place, model.index("internal")])
        r = least_squares(miss, xt[every], args=(xt, every), method="lm")
        xf = xt.copy()
        xf[every] = r.x
        dev = np.sqrt((r.fun ** 2).reshape(-1, 3).sum(axis=1))
        xp = model.project(xf, model.index("internal"), [0, 1])
        e_fit = model.physical_energy(xp)
        free = np.setdiff1d(np.arange(model.n_var), np.concatenate([model.index("xy", 1), model.index("z", 1)]))
        rigid = np.concatenate([model.index("lattice"), model.index("shear"), model.index("phi"), model.index("xy", 0), model.index("z", 0)])
        x_held, i_held = model.relax(xp, free=np.setdiff1d(free, model.index("shear")), closed=[0, 1], maxiter=800)
        x_free, i_free = model.relax(x_held, free=free, closed=[0, 1], maxiter=800)
        # from the upright cell: the shear polyfind's own energy takes, and where a start at DFT's shear ends
        x_loc, i_loc = model.relax(x0, free=free, closed=[0, 1], maxiter=800)
        xd = x0.copy()
        xd[model.index("shear")] = [sa, sb]
        e_rigid = model.physical_energy(xd)
        xd, _ = model.relax(xd, free=rigid, closed=[], maxiter=300)
        x_back, i_back = model.relax(xd, free=free, closed=[0, 1], maxiter=800)
        dft = json.load(open(PRED / f"mlip_candidate_{CANDIDATE}.json"))["all_trans_cell_held_upright"]
        gain = lambda x: (model.physical_energy(x) - e0) / per  # noqa: E731
        shear = lambda x: [float(v) for v in model.state(x).shear]  # noqa: E731
        out = {
            "what": f"the all-trans cell of {CANDIDATE} with polyfind's cell free to shear along its chains, and polyfind's energy at the cell DFT relaxed to",
            "energy": "kcal/mol per monomer against polyfind's polished upright cell, rigid pendants", "preset": PRESET,
            "monomer_A": cm, "upright_shear_gradient_kcal_per_A_per_cell": [float(v) for v in g0[model.index("shear")]],
            "polyfind_shear_free_from_upright": {"gain": gain(x_loc), "shear_A": shear(x_loc), "shear_monomers": [v / cm for v in shear(x_loc)],
                                                 "cell": model.state(x_loc).lattice.tolist(), "evaluations": int(i_loc["evaluations"]), "gradient_max": float(i_loc["gradient_max"])},
            "dft_relaxed_cell": {"frame": a.info["frame"], "cell_A_deg": [A, B, gam, c], "shear_A": [sa, sb], "shear_monomers": [sa / cm, sb / cm],
                                 "fit_rms_A": float(np.sqrt((dev ** 2).mean())), "fit_worst_A": float(dev.max()),
                                 "setting_angles_deg": model.state(xf).place[:, 2].tolist(), "polyfind_setting_angles_deg": st.place[:, 2].tolist()},
            "polyfind_at_dft_cell": {"as_fitted": (e_fit - e0) / per, "shear_held_rest_relaxed": gain(x_held), "shear_free": gain(x_free),
                                     "shear_free_ends_at_A": shear(x_free), "shear_free_ends_at_monomers": [v / cm for v in shear(x_free)]},
            "polyfind_from_dft_shear_with_its_own_chains": {"rigid_at_start": (e_rigid - e0) / per, "relaxed": gain(x_back), "ends_at_A": shear(x_back)},
            "dft": {"sheared_below_upright": -float(dft["above_the_sheared_all_trans_cell_kcal"]), "source": f"mlip_candidate_{CANDIDATE}.json, all_trans_cell_held_upright"},
            "seconds": time.time() - t0,
        }
    disk = json.load(open(out_path)) if out_path.exists() else {}
    disk["candidate"] = out
    out_path.write_text(json.dumps(disk, indent=1))
    p, q = out["polyfind_shear_free_from_upright"], out["polyfind_at_dft_cell"]
    print(f"{CANDIDATE} all-trans: polyfind with the shear free {p['gain']:+.4f} at {np.round(p['shear_monomers'], 3).tolist()} monomers; "
          f"DFT's cell is at {np.round(out['dft_relaxed_cell']['shear_monomers'], 3).tolist()} (fit {out['dft_relaxed_cell']['fit_rms_A']:.3f} A rms), "
          f"where polyfind is {q['shear_held_rest_relaxed']:+.4f} with the shear held and DFT {out['dft']['sheared_below_upright']:+.4f}")


def collect(args):
    from ase.io import read, write

    src, out_path = Path(args.frames), Path(args.out)
    table = json.loads((src / "polyfind.json").read_text())
    kept = []
    ratios = []
    for pack, block in table.items():
        zero = None
        for row in block["shears"]:
            f = src / f"{row['frame']}.dft.xyz"
            if not f.exists():
                raise SystemExit(f"{f} is not labelled yet")
            a = read(f)
            kept.append(a)
            e = float(a.info["REF_energy"]) / a.info["n_monomers"]
            zero = e if zero is None else zero
            row["dft"] = (e - zero) * EV_TO_KCAL
            row["dft_largest_force_eV_per_A"] = float(np.abs(a.arrays["REF_forces"]).max())
            if row["dft"] > 1e-3:
                row["ratio"] = row["polyfind"] / row["dft"]
                ratios.append(row["ratio"])
        block["ratio_range"] = [min(r["ratio"] for r in block["shears"] if "ratio" in r), max(r["ratio"] for r in block["shears"] if "ratio" in r)]
        print(f"{pack}: polyfind over DFT {block['ratio_range'][0]:.2f} to {block['ratio_range'][1]:.2f}")
        for row in block["shears"]:
            print(f"   shear ({row['shear_monomers'][0]:+.3f}, {row['shear_monomers'][1]:+.3f}) monomers: polyfind {row['polyfind']:+.4f}, DFT {row['dft']:+.4f}"
                  + (f", ratio {row['ratio']:.2f}" if "ratio" in row else ""))
    registered = json.load(open(PRED / "e5_polyfind_belt4.json"))["dipole_cell"]
    e5 = {}
    f_ref = src / "e5_beta_2a_reference.dft.xyz"
    if f_ref.exists():
        ref = read(f_ref)
        kept.append(ref)
        e_ref = float(ref.info["REF_energy"]) / ref.info["n_monomers"]
        for name, (part, slide) in E5.items():
            a = read(src / f"{name}.dft.xyz")
            kept.append(a)
            dft = (float(a.info["REF_energy"]) / a.info["n_monomers"] - e_ref) * EV_TO_KCAL
            pred = float(registered[part]["slides"][slide]["per_cell_monomer_kcal"])
            e5[name] = {"part": part, "slide": slide, "atoms": len(a), "polyfind_registered": pred, "dft": dft, "ratio": pred / dft}
            ratios.append(pred / dft)
            print(f"E5 {part} {slide:5s}: polyfind (registered, like for like) {pred:+.4f}, DFT {dft:+.4f} kcal/mol per monomer of the cell, ratio {pred / dft:.2f}")
    disk = json.load(open(out_path)) if out_path.exists() else {}
    disk.update(
        what="the energy of moving neighbouring chains along each other, polyfind against PBE-D3(BJ) single points on the same atoms",
        energy="kcal/mol per monomer of the cell, against the cell at no shear (no slide); rigid chains, nothing relaxed",
        preset=PRESET, dft="Quantum ESPRESSO 7.5, PBE-D3(BJ), SSSP 1.3.0 precision, 90/360 Ry, 4x2x8 (2x2x8 on the doubled cell): the anchor's setup",
        packings=table, e5_like_for_like=e5,
        e5_note="the registered test's zero is the unslid cell; here it is the doubled reference cell at the same k-point density, which is the same crystal",
        ratio_range=[min(ratios), max(ratios)], ratio_median=float(np.median(ratios)), cells=len(ratios))
    out_path.write_text(json.dumps(disk, indent=1))
    write(args.extxyz, kept, format="extxyz")
    print(f"polyfind over DFT on {len(ratios)} cells: {min(ratios):.2f} to {max(ratios):.2f}, median {np.median(ratios):.2f}; written {out_path} and {args.extxyz} ({len(kept)} frames)")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="step", required=True)
    s = sub.add_parser("frames")
    s.add_argument("--out", required=True, help="working directory for the DFT, outside the repository")
    s = sub.add_parser("candidate")
    s.add_argument("--out", default=str(PRED / "shear_registry_belt4.json"))
    s = sub.add_parser("collect")
    s.add_argument("--frames", required=True)
    s.add_argument("--out", default=str(PRED / "shear_registry_belt4.json"))
    s.add_argument("--extxyz", default=str(PRED / "mlip" / "shear_registry_frames_pbe_d3bj.extxyz"))
    args = ap.parse_args()
    {"frames": frames, "candidate": candidate, "collect": collect}[args.step](args)


if __name__ == "__main__":
    main()
