"""E5: the axial registry of beta-PVDF, on the structures the DFT runs.

A conformational kink pair shortens its chain by about 1 A, and in polyfind's mixed-conformation
cell what confines kink pairs to sheets is the energy of sliding a chain along its axis against
its neighbours (docs/CANDIDATE_MATERIALS.md section 8).  No anchor of the potential tests that
energy.  E5 of sarco/materials/ask-polyfind-dipole-registry-kink-v1.md asks Dipole for it at
PBE-D3(BJ), as single points on its accepted beta cell (``pvdf_beta_ref_v1``):

    a   the chain at (1/2, 1/2) of the two-chain cell slid rigidly along c by c/8, c/4, 3c/8, c/2;
    b   the c/2 cell of (a) relaxed with the slip held (ions and cell);
    c   ONE chain of the cell doubled along its short axis (four chains) slid by c/4 and c/2.

This script does two things.

``--write-structures DIR`` builds those cells from Dipole's record (``--dipole-result``), atom
for atom: the accepted cell with the moved chain's z shifted, nothing else touched, as extended
XYZ and as one JSON in the layout of Dipole's case files (``cell_A``, ``positions_A``,
``symbols``).

Without it, polyfind's numbers are computed and written to ``--out``, twice:

``own_cell``     the registered protocol: polyfind's own relaxed beta cell (the polished all-trans
                 reference of examples/dilute_kink.py), slid the same way;
``dipole_cell``  like for like: the cell set to Dipole's accepted one, polyfind's chains relaxed
                 in it (the atoms then sit on Dipole's to the printed deviation), and the same
                 rigid slides.  This removes the difference of the two minima from the comparison.

(b) has one polyfind number, own minimum to own constrained minimum, since the DFT's (b) is a
relaxation and ends at the DFT's.  Also stored, not requested: the c/2 slip with the chains
relaxed at the fixed cell, which lies between (a) and (b).

The c/2 cell is symmetric (every atom on a mirror plane normal to c, the two sublattices related
by a body-centring translation), so a relaxation started on it stays on it, here as in a DFT
code that keeps the symmetry it finds.  That is the registered protocol.  Each relaxation is
therefore run twice, as started and from a start with the symmetry broken (setting angles,
torsions and the sublattice offset perturbed, ``--seed``), and both are stored: the second is
what a relaxation without symmetry would reach.

Energies are kcal/mol with rigid pendants unless stated: per monomer of the cell (the ask's
unit), per monomer of the moved chain, and per cell in kcal/mol and mEh, the last being what
two SCF energies give directly.
"""
import argparse
import hashlib
import importlib.util
import json
import time
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment

from polyfind.fitting import PRESET_PARAMETERS

DS = importlib.util.spec_from_file_location("dk", Path(__file__).with_name("dilute_kink.py"))
DK = importlib.util.module_from_spec(DS)
DS.loader.exec_module(DK)
GS = importlib.util.spec_from_file_location("dg", Path(__file__).with_name("dft_cell_geometry.py"))
DG = importlib.util.module_from_spec(GS)
GS.loader.exec_module(DG)
MC = DK.MC

KCAL_PER_EH = 627.5094740631
ROW = "vdf10"
A_FRACTIONS = (0.125, 0.25, 0.375, 0.5)
C_FRACTIONS = (0.25, 0.5)
CURVE = tuple(k / 16 for k in range(1, 9))
ATOMS = ("C", "H", "F")


# ----------------------------------------------------------------------------- Dipole's cell
def dipole_reference(path: str, case: str) -> dict:
    """Dipole's accepted cell: lattice rows, positions, symbols, energy, and the record's hash."""
    raw = Path(path).read_bytes()
    d = json.loads(raw)
    rec = d["payload"]["records"][case]["reference"]
    return {"case": case, "cell_A": np.asarray(rec["cell_A"], dtype=float), "positions_A": np.asarray(rec["positions_A"], dtype=float),
            "symbols": [str(s) for s in rec["symbols"]], "energy_Eh": float(rec["energy_Eh"]), "result_sha256": hashlib.sha256(raw).hexdigest()}


def molecules(cell, pos, sym) -> list:
    """The atoms of each bonded chain (periodic bonds), in the order of each chain's first atom."""
    n = len(sym)
    parent = list(range(n))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i, j, _, _ in DG.bonds(cell, pos, sym):
        parent[find(i)] = find(j)
    groups = {}
    for i in range(n):
        groups.setdefault(find(i), []).append(i)
    return sorted(groups.values(), key=lambda g: g[0])


def slid_structures(ref: dict) -> dict:
    """``name -> {cell_A, positions_A, symbols, ...}`` of every E5 cell, from Dipole's reference."""
    cell, pos, sym = ref["cell_A"], ref["positions_A"], ref["symbols"]
    c = float(cell[2, 2])
    mol = molecules(cell, pos, sym)
    if [len(m) for m in mol] != [6, 6] or mol[0] != list(range(6)) or mol[1] != list(range(6, 12)):
        raise RuntimeError(f"expected two chains of six atoms in file order, found {mol}")
    frac = (pos[mol[1]].mean(axis=0) - pos[mol[0]].mean(axis=0)) @ np.linalg.inv(cell)
    if not np.allclose(frac[:2] % 1.0, [0.5, 0.5], atol=1e-3):
        raise RuntimeError(f"the second chain is not at (1/2, 1/2) of the first: {frac}")
    out = {}
    tag = {0.125: "c8", 0.25: "c4", 0.375: "3c8", 0.5: "c2"}
    for f in A_FRACTIONS:
        p = pos.copy()
        p[mol[1], 2] += f * c
        out[f"beta_a_slide_{tag[f]}"] = {"part": "a", "slide_fraction_of_c": f, "moved_atoms": mol[1], "cell_A": cell.copy(), "positions_A": p, "symbols": list(sym),
                                         "mesh": [8, 4, 16], "run": "single-point SCF" + ("; also the start of (b)" if f == 0.5 else ""), "monomers": 2, "moved_monomers": 1}
    cell2 = cell.copy()
    cell2[0] *= 2.0
    pos2 = np.concatenate([pos, pos + cell[0]])
    sym2 = list(sym) * 2
    out["beta_2a_reference"] = {"part": "c", "slide_fraction_of_c": 0.0, "moved_atoms": [], "cell_A": cell2.copy(), "positions_A": pos2.copy(), "symbols": sym2, "mesh": [4, 4, 16],
                                "run": "single-point SCF, the zero of (c); twice the accepted cell's energy to the SCF threshold if the meshes are equivalent", "monomers": 4, "moved_monomers": 0}
    for f in C_FRACTIONS:
        p = pos2.copy()
        p[mol[0], 2] += f * c
        out[f"beta_2a_slide_{tag[f]}"] = {"part": "c", "slide_fraction_of_c": f, "moved_atoms": mol[0], "cell_A": cell2.copy(), "positions_A": p, "symbols": sym2,
                                          "mesh": [4, 4, 16], "run": "single-point SCF", "monomers": 4, "moved_monomers": 1}
    return out


def write_structures(ref: dict, out_dir: Path, source: str) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    cells = slid_structures(ref)
    index = {}
    for name, s in cells.items():
        lat = " ".join(f"{v:.12f}" for v in s["cell_A"].ravel())
        moved = ",".join(str(i) for i in s["moved_atoms"]) or "none"
        lines = [str(len(s["symbols"])),
                 f'Lattice="{lat}" Properties=species:S:1:pos:R:3 pbc="T T T" source="{source}; E5 part {s["part"]}: atoms {moved} (0-based) moved along z by '
                 f'{s["slide_fraction_of_c"]} c; nothing else changed" mesh="{" ".join(str(m) for m in s["mesh"])}"']
        lines += [f"{e:2s} {p[0]:18.12f} {p[1]:18.12f} {p[2]:18.12f}" for e, p in zip(s["symbols"], s["positions_A"])]
        path = out_dir / f"{name}.xyz"
        path.write_text("\n".join(lines) + "\n")
        index[name] = {"file": path.name, "sha256": hashlib.sha256(path.read_bytes()).hexdigest(), "part": s["part"], "slide_fraction_of_c": s["slide_fraction_of_c"],
                       "moved_atoms": s["moved_atoms"], "atoms": len(s["symbols"]), "monomers": s["monomers"], "moved_monomers": s["moved_monomers"], "mesh": s["mesh"], "run": s["run"],
                       "initial": {"cell_A": s["cell_A"].tolist(), "positions_A": s["positions_A"].tolist(), "symbols": s["symbols"]}}
        print(f"  {path}: {len(s['symbols'])} atoms, part {s['part']}, slide {s['slide_fraction_of_c']} c, mesh {s['mesh']}")
    doc = {"schema": "polyfind_e5_structures_v1", "source": source, "reference_case": ref["case"], "reference_result_sha256": ref["result_sha256"],
           "reference_energy_Eh": ref["energy_Eh"], "frame": "x polar (short axis), y lateral (long axis), z chain; Dipole's own coordinates, unwrapped",
           "method": "as pvdf_beta_ref_v1: QE 7.6, PBE-D3(BJ) two-body, SSSP 1.3.0 precision, 90/360 Ry, SCF 1e-11 Ry",
           "reference": {"cell_A": ref["cell_A"].tolist(), "positions_A": ref["positions_A"].tolist(), "symbols": ref["symbols"], "mesh": [8, 4, 16],
                         "run": "none: the accepted cell, the zero of (a) and (b)"},
           "structures": index}
    (out_dir / "e5_structures.json").write_text(json.dumps(doc, indent=1) + "\n")
    print(f"  {out_dir / 'e5_structures.json'}")
    return doc


# ----------------------------------------------------------------------------- polyfind's cell
def real_atoms(model, x):
    """``(lattice rows, positions, symbols, chain index)`` of the real atoms of a mixed cell."""
    latn, _, _, Pn, _ = model.build(x)
    el = np.array([str(e) for e in list(model.pk.elements) * model.nc])
    chain = np.repeat(np.arange(model.nc), model.n)
    keep = np.array([e in ATOMS for e in el])
    return latn, Pn[keep], el[keep], chain[keep]


def axial_phase(model, x, c_mono: float) -> np.ndarray:
    """The height (A, modulo one monomer period) of each chain's CH2 carbons."""
    latn, P, el, chain = real_atoms(model, x)
    out = np.empty(model.nc)
    for s in range(model.nc):
        idx = np.where(chain == s)[0]
        ch2 = [i for k, i in enumerate(idx[:-1]) if el[i] == "C" and el[idx[k + 1]] == "H"]
        z = P[ch2, 2]
        out[s] = (np.angle(np.exp(2j * np.pi * z / c_mono).mean()) / (2 * np.pi) * c_mono) % c_mono
    return out


def slip(model, x, c_mono: float, moved: int, against: int) -> float:
    """The axial offset of chain ``moved`` from chain ``against``, in (-c/2, c/2]."""
    ph = axial_phase(model, x, c_mono)
    u = (ph[moved] - ph[against]) % c_mono
    return float(u - c_mono if u > c_mono / 2 + 1e-9 else u)


def deviation_from(model, x, ref: dict, n_mon: int) -> dict:
    """How far polyfind's atoms (cell set to Dipole's) sit from Dipole's: the largest and the
    r.m.s. distance over a one-to-one assignment of like atoms, the rigid shift of the whole
    cell and the two-fold turn about the chain axis being free."""
    latn, P, el, _ = real_atoms(model, x)
    cell = ref["cell_A"].copy()
    cell[2] *= n_mon
    Q = np.concatenate([ref["positions_A"] + k * ref["cell_A"][2] for k in range(n_mon)])
    qs = np.array(ref["symbols"] * n_mon)
    inv = np.linalg.inv(cell)
    best = None
    for turn in (np.diag([1.0, 1.0, 1.0]), np.diag([-1.0, -1.0, 1.0])):
        Pt = P @ turn.T
        i0 = int(np.where(el == "C")[0][0])
        for j0 in np.where(qs == "C")[0]:
            shifted = Pt + (Q[j0] - Pt[i0])
            total, worst = 0.0, 0.0
            for e in ATOMS:
                a, b = np.where(el == e)[0], np.where(qs == e)[0]
                if len(a) != len(b):
                    raise RuntimeError(f"atom counts differ for {e}: {len(a)} against {len(b)}")
                d = shifted[a][:, None, :] - Q[b][None, :, :]
                f = d @ inv
                d = (f - np.round(f)) @ cell
                dist = np.linalg.norm(d, axis=2)
                r, c = linear_sum_assignment(dist)
                total += float((dist[r, c] ** 2).sum())
                worst = max(worst, float(dist[r, c].max()))
            if best is None or worst < best["largest_A"]:
                best = {"largest_A": worst, "rms_A": float(np.sqrt(total / len(el)))}
    return best


def energies(dE_cell: float, n_cell_mon: int, n_moved_mon: int) -> dict:
    return {"dE_cell_kcal": float(dE_cell), "dE_cell_mEh": float(dE_cell / KCAL_PER_EH * 1e3), "per_cell_monomer_kcal": float(dE_cell / n_cell_mon),
            "per_moved_monomer_kcal": float(dE_cell / n_moved_mon)}


def rigid_slides(model, x, chain: int, c_mono: float, fractions, n_mon: int, reduce_to: int) -> dict:
    """Energies of chain ``chain`` slid rigidly by each fraction of the monomer period.  The model
    cell holds ``n_mon`` monomers per chain; ``reduce_to`` is the number of chains of the cell the
    DFT runs (one monomer per chain), so energies are per that cell."""
    e0 = model.physical_energy(x)
    iz = model.index("z", chain)[0]
    out = {}
    for f in fractions:
        y = x.copy()
        y[iz] += f * c_mono
        dE = (model.physical_energy(y) - e0) / n_mon  # per DFT cell: one monomer of every chain
        out[f"{f:.4f}"] = energies(dE, reduce_to, 1)
    return out


def relax_fixed_cell(model, x, maxiter: int, held_z=()):
    """Chains relaxed in the cell as it stands (closure exact); ``held_z`` chains keep their axial shift."""
    nc = model.nc
    held = [model.index("lattice"), model.index("c"), model.index("xy", nc - 1), model.index("z", nc - 1)] + [model.index("z", s) for s in held_z]
    free = np.setdiff1d(np.arange(model.n_var), np.concatenate(held))
    everyone = list(range(nc))
    x = model.project(x, model.index("internal"), everyone)
    xr, info = model.relax(x, free=free, closed=everyone, maxiter=maxiter)
    for _ in range(3):
        if "success" in info["message"].lower():
            break
        xr, more = model.relax(xr, free=free, closed=everyone, maxiter=maxiter)
        info = {**more, "iterations": info["iterations"] + more["iterations"]}
    return xr, info


def perturbed(model, x, seed: int = 0):
    """``x`` with its symmetry broken: setting angles within 3 degrees, torsions within 1 degree and
    the first chain's lateral position within a hundredth of the cell, closure restored."""
    rng = np.random.default_rng(seed)
    y = np.asarray(x, dtype=float).copy()
    y[model.index("phi")] += rng.uniform(-3.0, 3.0, model.nc)
    y[model.index("tors")] += rng.uniform(-1.0, 1.0, model.nc * model.P)
    y[model.index("xy", 0)] += rng.uniform(-0.01, 0.01, 2)
    return model.project(y, model.index("internal"), list(range(model.nc)))


def shape(model, x) -> dict:
    """What a broken symmetry shows in: the chains' turn about their axes and the torsions off trans."""
    st = model.state(x)
    turn = DK.wrap(st.place[:, 2] - 180.0 * np.round(st.place[:, 2] / 180.0))
    return {"setting_angle_turn_deg": [float(v) for v in turn], "largest_torsion_off_trans_deg": float(np.abs(DK.wrap(st.tors - 180.0)).max()),
            "lateral_offset_fraction": [float(v) for v in st.place[0, :2] - np.round(st.place[0, :2] * 2) / 2]}


def relax_held_slip(model, x, moved: int, against: int, target: float, c_of, maxiter: int, log):
    """Cell and chains relaxed with the slip of chain ``moved`` held at ``target`` monomer periods
    (``c_of(x)`` is the monomer period at ``x``): the axial shifts are held during each
    relaxation and corrected to the relaxed repeat between relaxations."""
    nc = model.nc
    held = np.concatenate([model.index("xy", against), model.index("z")])
    free = np.setdiff1d(np.arange(model.n_var), held)
    everyone = list(range(nc))
    info, total = {}, 0
    for sweep in range(8):
        x, info = model.relax(x, free=free, closed=everyone, maxiter=maxiter)
        total += info["iterations"]
        cm = c_of(x)
        u = slip(model, x, cm, moved, against)
        want = target * cm
        err = ((want - u + cm / 2) % cm) - cm / 2
        log(f"      sweep {sweep}: E {model.physical_energy(x):.5f}, cell {np.round(model.state(x).lattice, 4).tolist()}, slip {u:.6f} A of {want:.6f}; {info['message']}")
        if abs(err) < 2e-6 and "success" in info["message"].lower():
            break
        x = x.copy()
        x[model.index("z", moved)[0]] += err
    return x, {**info, "iterations": total, "slip_A": u, "slip_error_A": float(err)}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--dipole-result", default="../sarco/materials/electrostriction/results/qe_equilibrium_v1/pvdf_beta_ref_v1/result.json")
    ap.add_argument("--case", default="pvdf_beta_ref_v1")
    ap.add_argument("--write-structures", default="", help="directory: write the E5 cells built from Dipole's record and stop")
    ap.add_argument("--source", default="Dipole pvdf_beta_ref_v1 (sarco materials/electrostriction/results/qe_equilibrium_v1/pvdf_beta_ref_v1/result.json)")
    ap.add_argument("--preset", default="pvdf-dft-valence-belt4")
    ap.add_argument("--maxiter", type=int, default=600)
    ap.add_argument("--seed", type=int, default=0, help="seed of the symmetry-breaking perturbation")
    ap.add_argument("--no-broken", action="store_true", help="skip the relaxations from a symmetry-broken start")
    ap.add_argument("--references", default="docs/predictions/dilute_kink_belt4.json")
    ap.add_argument("--structures", default="docs/predictions/e5_structures/e5_structures.json", help="the index written by --write-structures, for the hashes")
    ap.add_argument("--out", default="docs/predictions/e5_polyfind_belt4.json")
    args = ap.parse_args()
    log = lambda s: print(s, flush=True)  # noqa: E731
    ref = dipole_reference(args.dipole_result, args.case)
    if args.write_structures:
        write_structures(ref, Path(args.write_structures), args.source)
        return
    t_start = time.time()
    letters, _ = DK.sequence(ROW)
    cell_d = ref["cell_A"]
    if not np.allclose(cell_d, np.diag(np.diag(cell_d)), atol=1e-6):
        raise RuntimeError("Dipole's cell is expected orthorhombic with c along z")
    with PRESET_PARAMETERS[args.preset].applied():
        polymer = DK.C.make_polymer(ROW, letters)
        stored = json.load(open(args.references))["all_trans_references"][ROW]
        d = DK.all_trans_reference(polymer, ROW, stored)
        model, x0, ch = DK.make_model(polymer, d, 1, 1)
        n_mon = ch.n_monomers
        c_of = lambda x: float(model.state(x).lattice[3]) / n_mon  # noqa: E731
        results = {"schema": "polyfind_e5_prediction_v1",
                   "potential": f"{args.preset}: induced dipoles, Ewald, force-shifted LJ, fluorine belt sites, PBE-D3 chain geometry, rigid bonds, 0 K; rigid pendants unless stated",
                   "frame": "x polar (short axis a), y lateral (long axis b), z chain",
                   "energy": "kcal/mol. per_cell_monomer: the cell's energy change over the monomers in the cell (the ask's unit); per_moved_monomer: over the monomers of the moved chain; "
                             "dE_cell: per DFT cell (12 atoms in a and b, 24 in c), also in mEh",
                   "model_cell": f"two chains of {n_mon} VDF monomers ({ROW}), the polished all-trans reference of examples/dilute_kink.py; energies reduced to one monomer per chain",
                   "dipole_reference": {"case": ref["case"], "result_sha256": ref["result_sha256"], "energy_Eh": ref["energy_Eh"],
                                        "cell_A": [float(cell_d[0, 0]), float(cell_d[1, 1]), 90.0, float(cell_d[2, 2])]}}
        if Path(args.structures).exists():
            idx = json.load(open(args.structures))
            results["structures"] = {k: {kk: v[kk] for kk in ("file", "sha256", "part", "slide_fraction_of_c", "atoms", "mesh")} for k, v in idx["structures"].items()}

        def block(label, model, x, describe):
            """Parts a and c (rigid slides) at the two-chain state ``x``, and the c/2 slip relaxed at the fixed cell."""
            cm = c_of(x)
            st = model.state(x)
            u0 = slip(model, x, cm, 0, 1)
            log(f"{label}: cell a {st.lattice[0]:.5f}, b {st.lattice[1]:.5f}, gamma {st.lattice[2]:.4f}, c per monomer {cm:.5f} A; "
                f"E {model.physical_energy(x) / (2 * n_mon):.5f} kcal/mol per monomer; the sublattices' axial offset {u0:+.2e} A")
            if abs(u0) > 1e-3:
                raise RuntimeError("the reference's two chains are not in register")
            out = {"description": describe, "cell": [float(st.lattice[0]), float(st.lattice[1]), float(st.lattice[2]), cm],
                   "energy_per_monomer_kcal": float(model.physical_energy(x) / (2 * n_mon))}
            t0 = time.time()
            curve = rigid_slides(model, x, 0, cm, CURVE, n_mon, 2)
            out["a"] = {"cell": "two chains, 12 atoms", "moved": "one chain (one sublattice), rigid", "zero": "the unslid cell",
                        "slides": {k: curve[f"{f:.4f}"] for f, k in zip(A_FRACTIONS, ("c/8", "c/4", "3c/8", "c/2"))}, "curve_fraction_of_c": curve}
            log("   a (one sublattice slid, rigid): " + ", ".join(f"{k} {v['per_cell_monomer_kcal']:+.3f}" for k, v in out["a"]["slides"].items())
                + f" kcal/mol per monomer of the cell ({v_last(out['a']['slides'])['dE_cell_mEh']:+.3f} mEh per cell at c/2); {time.time() - t0:.0f} s")
            # c: one chain of the cell doubled along its short axis
            t0 = time.time()
            d4 = dict(d, state=np.asarray(x, dtype=float), energy=float(model.physical_energy(x)))
            m4, x4, _ = DK.make_model(polymer, d4, 2, 1)
            c4 = rigid_slides(m4, x4, 0, cm, C_FRACTIONS, n_mon, 4)
            out["c"] = {"cell": "four chains (doubled along a), 24 atoms", "moved": "one chain of four, rigid", "zero": "the unslid doubled cell",
                        "slides": {k: c4[f"{f:.4f}"] for f, k in zip(C_FRACTIONS, ("c/4", "c/2"))}}
            log("   c (one chain of four slid, rigid): " + ", ".join(f"{k} {v['per_cell_monomer_kcal']:+.3f} ({v['per_moved_monomer_kcal']:+.2f} per monomer of the moved chain)"
                                                                   for k, v in out["c"]["slides"].items()) + f"; {time.time() - t0:.0f} s")
            out["anisotropy_c_over_a_per_moved_monomer_at_half"] = float(out["c"]["slides"]["c/2"]["per_moved_monomer_kcal"] / out["a"]["slides"]["c/2"]["per_moved_monomer_kcal"])
            # not requested: the c/2 slip with the chains relaxed in the fixed cell
            t0 = time.time()
            y = x.copy()
            y[model.index("z", 0)[0]] += 0.5 * cm
            yr, info = relax_fixed_cell(model, y, args.maxiter, held_z=(0,))
            dE = (model.physical_energy(yr) - model.physical_energy(x)) / n_mon
            out["half_slip_ions_relaxed_fixed_cell"] = {**energies(dE, 2, 1), "slip_A": abs(slip(model, yr, cm, 0, 1)), "requested": False, "symmetry": "kept (as started)",
                                                        "message": info["message"], "iterations": info["iterations"], **shape(model, yr)}
            log(f"   c/2 with the chains relaxed at the fixed cell (not requested): {dE / 2:+.3f} per monomer of the cell; slip {abs(slip(model, yr, cm, 0, 1)):.5f} A; "
                f"{info['iterations']} it, {info['message']}; {time.time() - t0:.0f} s")
            if not args.no_broken:
                t0 = time.time()
                yq, info = relax_fixed_cell(model, perturbed(model, yr, args.seed), args.maxiter, held_z=(0,))
                dEq = (model.physical_energy(yq) - model.physical_energy(x)) / n_mon
                sh = shape(model, yq)
                out["half_slip_ions_relaxed_fixed_cell"]["symmetry_broken"] = {**energies(dEq, 2, 1), "below_symmetric_per_cell_monomer_kcal": float((dE - dEq) / 2),
                                                                             "slip_A": abs(slip(model, yq, cm, 0, 1)), "message": info["message"], "iterations": info["iterations"], **sh}
                log(f"      from a symmetry-broken start: {dEq / 2:+.3f} ({(dEq - dE) / 2:+.3f}); chains turned {np.round(sh['setting_angle_turn_deg'], 2).tolist()} deg, "
                    f"torsions up to {sh['largest_torsion_off_trans_deg']:.2f} deg off trans; {info['iterations']} it, {info['message']}; {time.time() - t0:.0f} s")
            return out

        def v_last(slides):
            return slides["c/2"]

        # --- polyfind's own minimum
        own = block("own cell", model, x0, "polyfind's own relaxed beta cell (the registered protocol)")
        # b: the c/2 slip relaxed with the cell free, the slip held
        t0 = time.time()
        cm0 = c_of(x0)
        y = x0.copy()
        y[model.index("z", 0)[0]] += 0.5 * cm0
        log("   b (c/2, slip held, ions and cell relaxed):")
        yb, info = relax_held_slip(model, y, 0, 1, 0.5, c_of, args.maxiter, log)
        sb, s0 = model.state(yb), model.state(x0)
        dE = (model.physical_energy(yb) - model.physical_energy(x0)) / n_mon
        tilt_b, tilts, _ = MC.tilt_relaxation(model, yb)
        tilt_0 = float(d.get("tilt_per_cell", 0.0) or 0.0)
        dE_t = dE + (tilt_b - tilt_0) / n_mon
        def cell_change(st):
            return {"a": float(st.lattice[0] / s0.lattice[0] - 1), "b": float(st.lattice[1] / s0.lattice[1] - 1), "gamma_deg": float(st.lattice[2] - s0.lattice[2]),
                    "c": float(st.lattice[3] / s0.lattice[3] - 1)}

        own["b"] = {"cell": "two chains, 12 atoms", "held": "the slip (every atom's coordinate along c, in units of c)", "free": "the cell (a, b, gamma, c) and everything else",
                    "zero": "the relaxed unslid cell", "symmetry": "kept (as started): the registered protocol", **energies(dE, 2, 1), "with_pendant_tilts": energies(dE_t, 2, 1),
                    "relaxed_cell": [float(sb.lattice[0]), float(sb.lattice[1]), float(sb.lattice[2]), float(sb.lattice[3]) / n_mon], "cell_change": cell_change(sb),
                    "largest_tilt_deg": float(np.abs(tilts).max()), "slip_A": abs(info["slip_A"]), "message": info["message"], "iterations": info["iterations"], **shape(model, yb)}
        log(f"   b: {dE / 2:+.3f} per monomer of the cell rigid pendants, {dE_t / 2:+.3f} with pendant tilts (largest {np.abs(tilts).max():.2f} deg); a {own['b']['cell_change']['a']:+.2%}, "
            f"b {own['b']['cell_change']['b']:+.2%}, c {own['b']['cell_change']['c']:+.2%}, gamma {own['b']['cell_change']['gamma_deg']:+.3f} deg; {time.time() - t0:.0f} s")
        if not args.no_broken:
            t0 = time.time()
            log("   b from a symmetry-broken start:")
            yq, info = relax_held_slip(model, perturbed(model, yb, args.seed), 0, 1, 0.5, c_of, args.maxiter, log)
            sq = model.state(yq)
            dEq = (model.physical_energy(yq) - model.physical_energy(x0)) / n_mon
            tilt_q, tilts_q, _ = MC.tilt_relaxation(model, yq)
            sh = shape(model, yq)
            own["b"]["symmetry_broken"] = {**energies(dEq, 2, 1), "with_pendant_tilts": energies(dEq + (tilt_q - tilt_0) / n_mon, 2, 1), "below_symmetric_per_cell_monomer_kcal": float((dE - dEq) / 2),
                                           "relaxed_cell": [float(sq.lattice[0]), float(sq.lattice[1]), float(sq.lattice[2]), float(sq.lattice[3]) / n_mon], "cell_change": cell_change(sq),
                                           "largest_tilt_deg": float(np.abs(tilts_q).max()), "slip_A": abs(info["slip_A"]), "message": info["message"], "iterations": info["iterations"], **sh}
            log(f"   b from a symmetry-broken start: {dEq / 2:+.3f} ({(dEq - dE) / 2:+.3f}), {(dEq + (tilt_q - tilt_0) / n_mon) / 2:+.3f} with pendant tilts; chains turned "
                f"{np.round(sh['setting_angle_turn_deg'], 2).tolist()} deg, torsions up to {sh['largest_torsion_off_trans_deg']:.2f} deg off trans; a {cell_change(sq)['a']:+.2%}, "
                f"b {cell_change(sq)['b']:+.2%}, c {cell_change(sq)['c']:+.2%}; {info['iterations']} it, {info['message']}; {time.time() - t0:.0f} s")
        results["own_cell"] = own

        # --- like for like: Dipole's cell, polyfind's chains relaxed in it
        t0 = time.time()
        xd = x0.copy()
        xd[0], xd[1], xd[2], xd[3] = float(cell_d[0, 0]), float(cell_d[1, 1]), 90.0, float(cell_d[2, 2]) * n_mon
        xd, info = relax_fixed_cell(model, xd, args.maxiter)
        dev = deviation_from(model, xd, ref, n_mon)
        log(f"Dipole's cell: polyfind's chains relaxed in it ({info['iterations']} it, {info['message']}, {time.time() - t0:.0f} s); atoms sit {dev['largest_A']:.4f} A at most "
            f"({dev['rms_A']:.4f} r.m.s.) from Dipole's; energy above polyfind's own minimum {(model.physical_energy(xd) - model.physical_energy(x0)) / (2 * n_mon):+.4f} kcal/mol per monomer")
        dip = block("Dipole's cell", model, xd, "the cell set to Dipole's accepted one, polyfind's chains relaxed in it (like for like with the DFT single points)")
        dip["deviation_from_dipole_positions"] = dev
        dip["above_own_minimum_per_monomer_kcal"] = float((model.physical_energy(xd) - model.physical_energy(x0)) / (2 * n_mon))
        results["dipole_cell"] = dip
    results["acceptance"] = {"a, c": "sign, the ordering along the slide, and c per moved monomer above a per moved monomer; magnitudes within a factor 1.5",
                             "b": "sign; magnitude within a factor 1.5; b opening and a closing"}
    results["seconds"] = time.time() - t_start
    Path(args.out).write_text(json.dumps(results, indent=1) + "\n")
    log(f"written {args.out} ({time.time() - t_start:.0f} s)")


if __name__ == "__main__":
    main()
