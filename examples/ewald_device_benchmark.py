"""A large cell's evaluation on a GPU: is it right there, what does it buy, and is the card busy.

    python examples/ewald_device_benchmark.py                 # 16, 32 and 64 chains
    python examples/ewald_device_benchmark.py --chains 16     # a few minutes; start here
    python examples/ewald_device_benchmark.py --chains 16 32 64 128 --host-max-chains 16

One energy-and-gradient evaluation of a cell of many chains (``polyfind.mixedcell``) is, on a
CPU, four fifths the Ewald sums of the induced-dipole model (``polyfind.ewald.Ewald.terms``,
twice: the field and the dipole tensor for the solve, then the gradient at the solved dipoles)
and a tenth the solve between them (measured 2026-11-13, sixteen and thirty-two chains).  Those
three steps run on a device when the Ewald object is given a device array module, and the model
then evaluates everything else in its batched form: the pair sums of all chain pairs on the
device, and what is left for the host -- building the chains, the bonded lists, carrying the
gradient back to the torsions -- as a fixed number of array operations.  The first A100 run
(docs/predictions/ewald_device_a100_20261113) had only the three steps on the card and showed
why that is not enough: they ran 46 to 65 times faster and the evaluation 5 to 10, because 70%
of it was still the host's loops and the card sat idle meanwhile.

For each cell size this script

1. builds the cell (a tiling of the kinked crystal of ``docs/predictions/kink_migration_belt4.json``,
   every variable moved a little so that no two chains are copies of each other) and takes the
   arguments its evaluation hands to the Ewald sum;
2. checks the device against the host on exactly those arguments -- every energy term, the
   field, the tensor, the solved dipoles, every gradient -- and on the whole evaluation; and
   makes the two checks that need no host run, for cells the host is too slow for: the
   untouched tiling must have the stored two-chain cell's energy per chain, and the gradient
   must agree with the device's own energy differences;
3. times the three steps and the whole evaluation on both, the clock read only after the
   device has finished, the first device call (which compiles kernels, once per machine)
   reported apart; the host is timed at one BLAS thread, which is what the package sets, at
   eight, and at every core it is allotted, and the device is compared with the fastest;
4. reports the device memory the evaluation needed, the evaluation part by part (which parts
   are the host's), and -- the direct answer to "is the GPU being used" -- the card's own
   utilisation figure sampled while evaluations run back to back, for the batched form and, with
   ``--loop-form``, for the loop form it replaces;
5. times two things the port does not use yet, to decide whether it should: the symmetric
   rank-k update against the general matrix product it replaces, and a Cholesky solve against
   the LU solve.

It prints a table and writes everything to a JSON file (``--out``).  With more than one
device, ``--all-devices`` runs the same evaluation on all of them at once, one process each,
which is how a band of images or a set of finite differences would use them.

Nothing here needs a GPU to *import*: on a machine without one the host columns are measured
and the device columns say why they are empty.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import os
import platform
import subprocess
import sys
import time
import traceback
from pathlib import Path

import numpy as np
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parent.parent
TILINGS = {2: (1, 1), 4: (2, 1), 8: (2, 2), 16: (4, 2), 32: (8, 2), 64: (8, 4), 128: (16, 4), 256: (16, 8)}
ROW, PRESET = "vdf10", "pvdf-dft-valence-belt4"


def log(msg: str = "") -> None:
    print(msg, flush=True)


# --------------------------------------------------------------------------- the machine
def allotted_cores() -> dict:
    """How many CPUs this process may really use, and what says so.

    A container can see every CPU of its host and still be held to a few of them by a quota: a
    Runpod pod on a 256-thread machine reports 256 usable CPUs and is allotted 32, and a BLAS
    given 256 threads there ran a sixteen-chain evaluation slower than on one (measured
    2026-11-13).  So the count is the affinity mask cut down by the control group's quota and by
    the pod's own statement of its allotment."""
    out = {}
    try:
        out["affinity"] = len(os.sched_getaffinity(0))
    except AttributeError:
        out["affinity"] = os.cpu_count() or 1
    try:  # control groups, version 2 then version 1
        quota, period = open("/sys/fs/cgroup/cpu.max").read().split()[:2]
        if quota != "max":
            out["cgroup_quota"] = float(quota) / float(period)
    except Exception:
        try:
            quota = float(open("/sys/fs/cgroup/cpu/cpu.cfs_quota_us").read())
            period = float(open("/sys/fs/cgroup/cpu/cpu.cfs_period_us").read())
            if quota > 0:
                out["cgroup_quota"] = quota / period
        except Exception:
            pass
    env = os.environ.get("RUNPOD_CPU_COUNT", "")
    if env.isdigit() and int(env) > 0:
        out["RUNPOD_CPU_COUNT"] = int(env)
    limits = [out["affinity"]]
    if "cgroup_quota" in out:
        limits.append(max(1, int(out["cgroup_quota"])))
    if "RUNPOD_CPU_COUNT" in out:
        limits.append(out["RUNPOD_CPU_COUNT"])
    out["cores"] = max(1, min(limits))
    return out


def machine(cp) -> dict:
    info = {"python": platform.python_version(), "numpy": np.__version__, "platform": platform.platform()}
    try:
        import scipy

        info["scipy"] = scipy.__version__
    except Exception:
        pass
    try:
        info["cpu_count"] = os.cpu_count()
        info["cpus_usable"] = len(os.sched_getaffinity(0))
    except Exception:
        pass
    info["cpus_allotted"] = allotted_cores()
    try:
        with open("/proc/cpuinfo") as fh:
            names = {line.split(":", 1)[1].strip() for line in fh if line.startswith("model name")}
        info["cpu"] = "; ".join(sorted(names))
    except Exception:
        pass
    try:
        with open("/proc/meminfo") as fh:
            info["memory_gb"] = round(int(fh.readline().split()[1]) / 1e6, 1)
    except Exception:
        pass
    try:
        from threadpoolctl import threadpool_info

        info["blas"] = [{k: lib.get(k) for k in ("internal_api", "num_threads", "version", "threading_layer")}
                        for lib in threadpool_info() if lib.get("user_api") == "blas"]
    except Exception:
        pass
    try:
        out = subprocess.run(["nvidia-smi", "--query-gpu=index,name,memory.total,driver_version", "--format=csv,noheader"],
                             capture_output=True, text=True, timeout=30)
        info["nvidia_smi"] = [line.strip() for line in out.stdout.splitlines() if line.strip()]
    except Exception as exc:
        info["nvidia_smi"] = f"not available ({type(exc).__name__})"
    if cp is not None:
        info["device_module"] = f"{cp.__name__} {getattr(cp, '__version__', '')}".strip()
        try:
            dev = cp.cuda.Device()
            free, total = dev.mem_info
            info["device"] = {"id": int(dev.id), "name": cp.cuda.runtime.getDeviceProperties(dev.id)["name"].decode(),
                              "memory_gb": round(total / 1e9, 1), "free_gb": round(free / 1e9, 1),
                              "compute_capability": str(dev.compute_capability), "visible_devices": int(cp.cuda.runtime.getDeviceCount())}
        except Exception as exc:
            info["device"] = f"not readable ({type(exc).__name__}: {exc})"
    return info


def device_module(want: str):
    """``(cupy, None)`` when a device can be used, ``(None, why not)`` otherwise."""
    if want == "none":
        return None, "not asked for (--device none)"
    if want == "double":  # the tests' stand-in: NumPy underneath, a device's strictness on top
        spec = importlib.util.spec_from_file_location("device_double", ROOT / "tests" / "device_double.py")
        dd = importlib.util.module_from_spec(spec)
        sys.modules["device_double"] = dd
        spec.loader.exec_module(dd)
        return dd, None
    try:
        import cupy as cp
    except Exception as exc:
        return None, f"CuPy cannot be imported ({type(exc).__name__}: {exc})"
    try:
        n = cp.cuda.runtime.getDeviceCount()
    except Exception as exc:
        return None, f"no CUDA device ({type(exc).__name__}: {exc})"
    if n < 1:
        return None, "no CUDA device"
    if want != "auto":
        cp.cuda.Device(int(want)).use()
    return cp, None


# --------------------------------------------------------------------------- the cell
def load_drivers():
    spec = importlib.util.spec_from_file_location("dk", ROOT / "examples" / "dilute_kink.py")
    dk = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(dk)
    return dk


def build(dk, nc: int, seed: int = 0):
    """``(model, x, x_tiled, e_two)``: the tiling with ``nc`` chains of the kinked crystal; ``x`` is
    the tiling with every variable moved by a hundredth of its scale (a tenth of a degree in a
    torsion, a hundredth of an angstrom), ``x_tiled`` the tiling itself, and ``e_two`` the stored
    energy of the two-chain cell it is a tiling of (kcal/mol per two chains)."""
    from polyfind import mixedcell as MC
    from polyfind import pack as pack_mod

    nx, ny = TILINGS[nc]
    stored = json.load(open(ROOT / "docs" / "predictions" / "kink_migration_belt4.json"))
    d = stored["kinked_references"][f"{ROW}@0"]
    letters, _ = dk.sequence(ROW)
    polymer = dk.C.make_polymer(ROW, letters)
    ch = dk.reference_chain(polymer, d)
    pk = pack_mod.CrystalPacker(ch, n_chains=nc, **dk.KW)
    st, flips = MC.tile_state(MC.MixedState.from_vector(d["state"], 2, len(d["tors"])), np.asarray(d["flips"]) > 0.5, nx, ny)
    model = MC.MixedModel(polymer, "ref", pk, flips, [ch] * nc, angle_stiffness=dk.angle_stiffness())
    x_tiled = st.vector()
    x = x_tiled + 0.01 * np.random.default_rng(seed).standard_normal(x_tiled.size) * model.scale(x_tiled)
    return model, x, x_tiled, float(d["energy"])


def evaluate(model, x, xp=np, repeat: int = 1):
    """``(caught, energy, gradient, seconds)`` of the fastest of ``repeat`` whole evaluations.

    The dipole step is caught on its way through: ``caught`` holds the arguments the evaluation
    hands to the Ewald sum, what the step returned, and the seconds it took inside this very
    evaluation.  Before each evaluation the Ewald object's kept phases are dropped: a relaxation
    never evaluates one geometry twice, and a repeat that found the last repeat's phases waiting
    would be timing something no run does."""
    from polyfind import backend as bk

    pk = model.pk
    best = None
    for _ in range(max(1, repeat)):
        caught = {}
        inner = pk._polarize

        def spy(Pn, q_cell, latn, flip, grad=False, charge_grad=False, **kw):
            caught.update(Pn=np.array(Pn), q=np.array(q_cell), lat=np.array(latn), flip=np.array(flip))
            bk.synchronize(xp)
            t0 = time.perf_counter()
            out = inner(Pn, q_cell, latn, flip, grad=grad, charge_grad=charge_grad, **kw)
            bk.synchronize(xp)
            caught["seconds"] = time.perf_counter() - t0
            caught["out"] = out
            return out

        pk._ewald.release()
        pk._polarize = spy
        try:
            bk.synchronize(xp)
            t0 = time.perf_counter()
            e, g = model.value_and_grad(x)
            bk.synchronize(xp)
            seconds = time.perf_counter() - t0
        finally:
            del pk._polarize  # back to the class's own
        if best is None or seconds < best[3]:
            best = (caught, float(e), np.asarray(g), seconds)
    return best


# --------------------------------------------------------------------------- timing
def clock(fn, xp, repeat: int):
    """``(best seconds, all seconds, last result)``; the device is waited for before each read."""
    from polyfind import backend as bk

    times, out = [], None
    for _ in range(max(1, repeat)):
        bk.synchronize(xp)
        t0 = time.perf_counter()
        out = fn()
        bk.synchronize(xp)
        times.append(time.perf_counter() - t0)
    return min(times), times, out


def steps(pk, caught, xp, repeat: int):
    """The three steps of the dipole evaluation on the array module ``xp``, timed one by one on the
    captured arguments, and their results on the host for comparison."""
    from polyfind import backend as bk
    from polyfind.forcefield import COULOMB

    ew = pk._ewald
    idx = pk._es_idx
    Pa, qa = (caught["Pn"], caught["q"]) if idx is None else (caught["Pn"][idx], caught["q"][idx])
    na = Pa.shape[0]
    scale = pk._thole_scale_es if xp is np else xp.asarray(pk._thole_scale_es)
    pref = COULOMB / pk.eps_r
    E0 = np.asarray(caught["out"][2])
    E0a = E0 if idx is None else E0[idx]
    keep = {}

    def tensor_step():
        ew.release()  # a new geometry in a real run: the phases are built here
        keep["t0"] = ew.terms(Pa, qa, caught["lat"], thole=pk._thole, thole_scale=scale, field=True, tensor=True, device_tensor=True)
        return keep["t0"]

    def solve_step():
        A = -keep["t0"].tensor.transpose(0, 2, 1, 3).reshape(3 * na, 3 * na)
        diag = xp.arange(3 * na)
        A[diag, diag] += xp.asarray(np.repeat(pref / pk._alpha_es, 3))
        keep["A"] = A
        rhs = xp.asarray(E0a.ravel())
        sol = bk.cholesky_solve(xp, A, rhs) if xp is not np else None  # as CrystalPacker._polarize solves it
        keep["cholesky"] = sol is not None
        keep["p"] = bk.to_numpy(xp.linalg.solve(A, rhs) if sol is None else sol).reshape(na, 3)
        return keep["p"]

    def gradient_step():
        return ew.terms(Pa, qa, caught["lat"], dipoles=keep["p"], thole=pk._thole, thole_scale=scale, grad=True)

    out = {"sites": int(na)}
    out["tensor_s"], out["tensor_all_s"], t0 = clock(tensor_step, xp, repeat)
    out["solve_s"], out["solve_all_s"], p = clock(solve_step, xp, repeat)
    out["gradient_s"], out["gradient_all_s"], t1 = clock(gradient_step, xp, repeat)
    out["k_vectors"], out["images"] = int(t0.n_kvec), int(t0.n_image)
    out["solve_method"] = "cholesky" if keep["cholesky"] else "lu"
    results = {"energy_terms": {n: float(getattr(t0, n)) for n in ("real", "recip", "self_", "background", "surface")},
               "field": np.asarray(t0.field), "tensor": t0.tensor, "dipoles": p,
               "gradient_terms": {n: float(getattr(t1, n)) for n in ("real", "recip", "self_", "background", "surface")},
               "grad_coords": np.asarray(t1.grad_coords), "grad_lattice": np.asarray(t1.grad_lattice), "A": keep["A"]}
    return out, results


def differences(host: dict, dev: dict, xp) -> dict:
    """Largest difference of every result, as a fraction of the result's own largest magnitude."""
    from polyfind import backend as bk

    out = {}
    for group in ("energy_terms", "gradient_terms"):
        for name, a in host[group].items():
            b = dev[group][name]
            out[f"{group}.{name}"] = abs(a - b) / max(1.0, abs(a))
    for name in ("field", "dipoles", "grad_coords", "grad_lattice"):
        a, b = host[name], dev[name]
        out[name] = float(np.abs(a - b).max() / max(1e-300, np.abs(a).max()))
    # the tensor is compared where it is: the host's goes over once
    a = xp.asarray(host["tensor"])
    out["tensor"] = float(bk.to_numpy(abs(a - dev["tensor"]).max()) / float(np.abs(host["tensor"]).max()))
    out["largest"] = max(out.values())
    return out


def directional_check(model, x, g, seed: int = 1, h: float = 2e-3) -> dict:
    """The gradient against the energy's own fourth-order central difference along one random
    direction of every variable (scaled as the minimiser sees them).

    The agreement to expect is about one part in 1e5, not rounding: the Ewald real-space sum
    is cut at ``erfc(alpha rc) ~ 1e-8``, so the energy jumps by about 1e-7 kcal/mol whenever a
    pair crosses the cutoff, and a difference over a short step sees those jumps (measured on the
    host: 3e-6 at four chains and 2e-5 at sixteen with this step, 1e-4 at a tenth of it)."""
    d = np.random.default_rng(seed).standard_normal(x.size) * model.scale(x)
    d /= np.linalg.norm(d / model.scale(x))
    f = lambda t: float(model.value_and_grad(x + t * d)[0])  # noqa: E731
    num = (-f(2 * h) + 8.0 * f(h) - 8.0 * f(-h) + f(-2 * h)) / (12.0 * h)
    ana = float(g @ d)
    return {"analytic": ana, "difference": float(num), "relative_error": abs(num - ana) / max(1e-12, abs(ana))}


HOST_PARTS = ("build", "valence", "exclusions", "assembly")  # of the batched form; the rest is on the array module


def parts_of(model, x, xp) -> dict:
    """One batched evaluation with the clock read at the end of each part (the device waited for
    each time, which an ordinary evaluation does not do)."""
    model.pk._ewald.release()
    model.timings = {}
    try:
        model.value_and_grad(x)
        parts = dict(model.timings)
    finally:
        model.timings = None
    parts["host_parts"] = sum(parts.get(k, 0.0) for k in HOST_PARTS)
    parts["all"] = sum(v for k, v in parts.items() if k != "host_parts")
    return parts


def busy_check(model, x, xp, seconds: float) -> dict:
    """Evaluations back to back for ``seconds``, each at a new geometry as a relaxation's are, with
    the card's own utilisation figure (``nvidia-smi``: the share of the last sample period in
    which a kernel was running) read five times a second meanwhile."""
    from polyfind import backend as bk

    proc = None
    if hasattr(xp, "cuda"):
        try:
            proc = subprocess.Popen(["nvidia-smi", "--query-gpu=index,utilization.gpu", "--format=csv,noheader,nounits", "-lms", "200"],
                                    stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
        except OSError:
            proc = None
    rng = np.random.default_rng(5)
    sc = model.scale(x)
    model.value_and_grad(x)  # warm
    bk.synchronize(xp)
    time.sleep(1.2)  # let the figure forget whatever ran before
    t0 = time.perf_counter()
    n = 0
    while time.perf_counter() - t0 < seconds:
        model.value_and_grad(x + 1e-4 * rng.standard_normal(x.size) * sc)
        n += 1
    bk.synchronize(xp)
    wall = time.perf_counter() - t0
    out = {"seconds": wall, "evaluations": n, "per_evaluation_s": wall / n}
    if proc is not None:
        proc.terminate()
        per_card: dict = {}
        for line in proc.communicate()[0].splitlines():
            try:
                card, value = (v.strip() for v in line.split(","))
                per_card.setdefault(int(card), []).append(float(value))
            except ValueError:
                continue
        # the samples from one period after the start (the figure lags by its period) to the end
        skip = 6 + 5
        busiest = max(per_card.values(), key=lambda v: sum(v[skip:]) / max(1, len(v[skip:])), default=[])
        kept = busiest[skip:]
        if kept:
            out.update(utilization_mean=sum(kept) / len(kept), utilization_min=min(kept), utilization_max=max(kept),
                       utilization_samples=len(kept))
    return out


def busy_child(args) -> None:
    """One of several processes on a card: build the cell, warm up, wait for the common start, then
    evaluate back to back until the common end and write how many evaluations that was."""
    from polyfind import backend as bk
    from polyfind.fitting import PRESET_PARAMETERS

    bk.set_backend("cpu")
    cp, why = device_module(args.device)
    if cp is None:
        raise SystemExit(f"no device: {why}")
    dk = load_drivers()
    start, nc = args.busy_child, args.chains[0]
    with PRESET_PARAMETERS[PRESET].applied():
        model, x, _, _ = build(dk, nc)
        ew = model.pk._ewald
        ew.xp, ew.device_min_sites = cp, 0
        rng = np.random.default_rng(os.getpid())
        sc = model.scale(x)
        for _ in range(2):
            model.value_and_grad(x + 1e-4 * rng.standard_normal(x.size) * sc)
        bk.synchronize(cp)
        ready = time.time()
        while time.time() < start:
            time.sleep(0.01)
        t0 = time.perf_counter()
        n = 0
        while time.time() < start + args.busy_seconds:
            model.value_and_grad(x + 1e-4 * rng.standard_normal(x.size) * sc)
            n += 1
        bk.synchronize(cp)
        wall = time.perf_counter() - t0
    Path(args.out).write_text(json.dumps({"evaluations": n, "seconds": wall, "ready_before_start_s": start - ready}))


def several_processes(args, nc: int, count: int) -> dict:
    """``count`` processes on one card at once, each its own cell of ``nc`` chains and its own
    evaluations: what the card delivers when one process leaves it partly idle.  The processes
    are started fresh (a process that has used the device cannot fork ones that use it), given
    time to build and warm up, and timed over one common window."""
    lead = 14.0 + 2.0 * count + (20.0 if nc >= 64 else 0.0)
    start = time.time() + lead
    procs = []
    for i in range(count):
        out = Path(args.out).with_suffix(f".proc{nc}_{count}_{i}.json")
        cmd = [sys.executable, __file__, "--chains", str(nc), "--device", args.device, "--busy-child", repr(start),
               "--busy-seconds", str(args.busy_seconds), "--out", str(out)]
        procs.append((out, subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)))
    meter = None
    try:
        meter = subprocess.Popen(["nvidia-smi", "--query-gpu=index,utilization.gpu", "--format=csv,noheader,nounits", "-lms", "200"],
                                 stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
    except OSError:
        pass
    t_meter = time.time()
    runs = []
    for out, p in procs:
        err = p.communicate()[1]
        runs.append(json.loads(out.read_text()) if p.returncode == 0 and out.exists() else {"failed": (err or "")[-400:]})
        out.unlink(missing_ok=True)
    res = {"processes": count, "chains": nc, "runs": runs}
    if meter is not None:
        meter.terminate()
        per_card: dict = {}
        for line in meter.communicate()[0].splitlines():
            try:
                card, value = (v.strip() for v in line.split(","))
                per_card.setdefault(int(card), []).append(float(value))
            except ValueError:
                continue
        first = int((start - t_meter + 1.5) / 0.2)  # the samples of the common window, less the meter's lag
        last = int((start - t_meter + args.busy_seconds) / 0.2)
        kept = max(per_card.values(), key=lambda v: sum(v[first:last]), default=[])[first:last]
        if kept:
            res.update(utilization_mean=sum(kept) / len(kept), utilization_min=min(kept), utilization_max=max(kept))
    good = [r for r in runs if "evaluations" in r]
    if good:
        res["evaluations_per_second"] = sum(r["evaluations"] / r["seconds"] for r in good)
        res["per_process_s"] = [r["seconds"] / max(1, r["evaluations"]) for r in good]
    return res


def tiling_check(model, x_tiled, e_two: float) -> dict:
    """The untouched tiling against the two-chain cell it tiles: the same energy per chain, whatever
    the size of the cell.  The stored number is the host's, so this holds a device to the host at
    sizes the host is not run at."""
    e = float(model.physical_energy(x_tiled)) * 2.0 / model.nc
    return {"per_two_chains": e, "stored_two_chain_cell": e_two, "difference": abs(e - e_two)}


def extras(cp, results_dev, n_sites: int, n_k: int) -> dict:
    """Two candidates the port does not use, timed on this cell's sizes with the device waited for."""
    from polyfind import backend as bk

    out = {}
    try:  # the symmetric product of the reciprocal tensor: general product against rank-k update
        rng = cp.random.default_rng(0)
        M = rng.standard_normal((n_sites, n_k), dtype=cp.float64)
        t_gemm, _, G = clock(lambda: M @ M.T, cp, 3)
        from cupy import cublas

        # out := a.T a for the (K, N) column-ordered view of M, no copy; one triangle is filled and
        # the other left at zero, so the whole matrix is the sum with its transpose less the diagonal
        t_syrk, _, S = clock(lambda: cublas.syrk("T", M.T), cp, 3)
        S = S + S.T - cp.diag(cp.diag(S))
        out["symmetric_product"] = {"shape": [n_sites, n_k], "gemm_s": t_gemm, "syrk_s": t_syrk,
                                    "gflops_gemm": 2.0 * n_sites * n_sites * n_k / t_gemm / 1e9,
                                    "largest_difference": float(bk.to_numpy(abs(S - G).max()) / float(bk.to_numpy(abs(G).max())))}
        M = G = S = None  # free them before the solve's matrices are copied
    except Exception as exc:
        out["symmetric_product"] = f"failed: {type(exc).__name__}: {exc}"
    try:  # the solve: LU as the port does it, against Cholesky (the matrix is symmetric positive definite)
        A = results_dev["A"]
        b = cp.random.default_rng(1).standard_normal(A.shape[0], dtype=cp.float64)
        t_lu, _, x_lu = clock(lambda: cp.linalg.solve(A, b), cp, 3)
        from cupyx import lapack

        asym = float(bk.to_numpy(abs(A - A.T).max()) / float(bk.to_numpy(abs(A).max())))
        t_ch, _, x_ch = clock(lambda: lapack.posv(A, b), cp, 3)  # (it copies both, as the LU solve does)
        out["solve"] = {"unknowns": int(A.shape[0]), "lu_s": t_lu, "cholesky_s": t_ch, "asymmetry_of_matrix": asym,
                        "largest_difference": float(bk.to_numpy(abs(x_lu - x_ch).max()) / float(bk.to_numpy(abs(x_lu).max())))}
    except Exception as exc:
        out["solve"] = f"failed: {type(exc).__name__}: {exc}"
    return out


# --------------------------------------------------------------------------- one cell size
def one_size(dk, nc: int, cp, args) -> dict:
    row: dict = {"chains": nc, "tiling": list(TILINGS[nc])}
    t0 = time.perf_counter()
    model, x, x_tiled, e_two = build(dk, nc)
    pk, ew = model.pk, model.pk._ewald
    row["variables"] = int(x.size)
    row["build_s"] = time.perf_counter() - t0
    log(f"\n== {nc} chains ({TILINGS[nc][0]}x{TILINGS[nc][1]}), {x.size} variables; built in {row['build_s']:.1f} s")
    do_host = nc <= args.host_max_chains
    host_res = caught = None
    if do_host:
        # The package holds the BLAS to one thread (polyfind/__init__.py: its parallelism is
        # processes, and many threads on small solves were found to be slower).  For a cell this
        # size threads help, so the host is timed at several counts and a device is compared with
        # the fastest of them.
        ew.xp = np
        row["host_runs"] = []
        for threads in args.host_threads:
            with threadpool_limits(limits=threads, user_api="blas"):
                c, e_host, g_host, whole = evaluate(model, x)
                rep = args.repeat if whole < 20.0 else 1
                if rep > 1:  # the fastest of the repeats, its dipole step with it
                    again = evaluate(model, x, np, rep - 1)
                    if again[3] < whole:
                        c, whole = again[0], again[3]
                run = {"blas_threads": threads, "whole_s": whole, "dipole_step_s": c["seconds"], "energy": e_host}
                log(f"   host, {threads} BLAS thread{'s' if threads > 1 else ''}: whole evaluation {whole:.2f} s, of which the dipole step "
                    f"{c['seconds']:.2f} s; energy {e_host:.6f}")
                st, res = steps(pk, c, np, rep)
                run.update(st)
                log(f"      tensor and field {st['tensor_s']:.2f} s, solve {st['solve_s']:.2f} s, gradient {st['gradient_s']:.2f} s "
                    f"({st['sites']} sites, {st['k_vectors']} k vectors)")
                if caught is None:  # the first run's arguments and results are the reference
                    caught, host_res = c, res
                    host_res.pop("A")
                    row["sites"], row["k_vectors"] = st["sites"], st["k_vectors"]
                    if whole < 20.0:
                        run["gradient_check"] = directional_check(model, x, g_host)
                        run["tiling_check"] = tiling_check(model, x_tiled, e_two)
                        log(f"      gradient against its own energy difference: relative error {run['gradient_check']['relative_error']:.1e}; "
                            f"tiling against the two-chain cell: {run['tiling_check']['difference']:.1e} kcal/mol")
                del res
                row["host_runs"].append(run)
                ew.release()
        row["host"] = min(row["host_runs"], key=lambda r: r["whole_s"])  # the faster one: what a device has to beat
        # the batched form on the host, at the thread count the loops were fastest with: the same
        # numbers to rounding, and what the host itself gains from the reorganisation
        with threadpool_limits(limits=row["host"]["blas_threads"], user_api="blas"):
            model.xp = np
            try:
                c, e_b, g_b, whole_b = evaluate(model, x, np, 1 if row["host"]["whole_s"] > 20.0 else args.repeat)
                row["host_batched"] = {"blas_threads": row["host"]["blas_threads"], "whole_s": whole_b, "dipole_step_s": c["seconds"],
                                       "energy": e_b, "parts_s": parts_of(model, x, np),
                                       "difference_from_loops": {"energy": abs(e_b - e_host) / max(1.0, abs(e_host)),
                                                                 "gradient": float(np.abs(g_b - g_host).max() / np.abs(g_host).max())}}
            finally:
                model.xp = None
            hb = row["host_batched"]
            log(f"   host, batched form, {hb['blas_threads']} BLAS thread(s): whole evaluation {whole_b:.2f} s, of which the dipole step "
                f"{c['seconds']:.2f} s; against the loops: energy {hb['difference_from_loops']['energy']:.1e}, gradient "
                f"{hb['difference_from_loops']['gradient']:.1e}")
            ew.release()
    if cp is None:
        return row
    pool = cp.get_default_memory_pool() if hasattr(cp, "get_default_memory_pool") else None
    free_pool = pool.free_all_blocks if pool is not None else (lambda: None)
    pool_gb = (lambda: pool.total_bytes() / 1e9) if pool is not None else (lambda: float("nan"))
    try:
        ew.xp, ew.device_min_sites = cp, 0
        free_pool()
        first = evaluate(model, x, cp)[3]  # compiles the kernels this machine has not compiled before
        caught_dev, e_dev, g_dev, whole = evaluate(model, x, cp, args.repeat)
        polar = caught_dev["seconds"]
        if caught is None:  # no host run: the arguments are the device evaluation's own
            caught = caught_dev
        dev = {"first_call_s": first, "whole_s": whole, "dipole_step_s": polar, "energy": float(e_dev),
               "pool_gb_after_evaluation": pool_gb()}
        log(f"   device: whole evaluation {whole:.3f} s (first call {first:.1f} s), of which the dipole step {polar:.3f} s; energy {float(e_dev):.6f}")
        st, dev_res = steps(pk, caught, cp, args.repeat)
        dev.update(st)
        dev["pool_gb_peak"] = pool_gb()
        if hasattr(cp, "cuda"):  # everything the process holds on the card: the pool, the libraries' workspaces, the context
            free, total = cp.cuda.Device().mem_info
            dev["device_gb_in_use"], dev["device_gb_total"] = (total - free) / 1e9, total / 1e9
        row["sites"], row["k_vectors"] = st["sites"], st["k_vectors"]
        log(f"   device: tensor and field {st['tensor_s']:.3f} s, solve {st['solve_s']:.3f} s, gradient {st['gradient_s']:.3f} s; "
            f"memory pool {dev['pool_gb_peak']:.2f} GB")
        dev["parts_s"] = parts_of(model, x, cp)
        ps = dev["parts_s"]
        log("   device: the evaluation part by part (ms): " + ", ".join(f"{k} {1e3 * v:.1f}" for k, v in ps.items() if k not in ("host_parts", "all"))
            + f"; the host's parts ({', '.join(HOST_PARTS)}) {1e3 * ps['host_parts']:.1f} of {1e3 * ps['all']:.1f}")
        if args.busy_seconds > 0:
            dev["busy"] = busy_check(model, x, cp, args.busy_seconds)
            b = dev["busy"]
            log(f"   device: {b['evaluations']} evaluations back to back in {b['seconds']:.1f} s ({1e3 * b['per_evaluation_s']:.0f} ms each)"
                + (f"; card utilisation {b['utilization_mean']:.0f}% (lowest sample {b['utilization_min']:.0f}%, highest {b['utilization_max']:.0f}%, "
                   f"{b['utilization_samples']} samples)" if "utilization_mean" in b else "; no utilisation figure on this machine"))
            if args.loop_form:  # the form this replaces, on the same card: sums on the device, loops on the host
                model.xp = "loops"
                try:
                    whole_l = evaluate(model, x, cp, args.repeat)
                    dev["loop_form"] = {"whole_s": whole_l[3], "dipole_step_s": whole_l[0]["seconds"],
                                        "busy": busy_check(model, x, cp, args.busy_seconds)}
                finally:
                    model.xp = None
                lf = dev["loop_form"]
                log(f"   device, loop form (the sums on the card, the rest in host loops, as before): whole evaluation {lf['whole_s']:.3f} s"
                    + (f"; card utilisation {lf['busy']['utilization_mean']:.0f}%" if "utilization_mean" in lf["busy"] else ""))
        dev["gradient_check"] = directional_check(model, x, g_dev)
        dev["tiling_check"] = tiling_check(model, x_tiled, e_two)
        log(f"   device: gradient against its own energy difference: relative error {dev['gradient_check']['relative_error']:.1e}; "
            f"tiling against the two-chain cell: {dev['tiling_check']['difference']:.1e} kcal/mol")
        if host_res is not None:
            diff = differences(host_res, dev_res, cp)
            diff["whole.energy"] = abs(e_dev - row["host"]["energy"]) / max(1.0, abs(row["host"]["energy"]))
            diff["whole.gradient"] = float(np.abs(g_dev - g_host).max() / np.abs(g_host).max())
            diff["largest"] = max(v for k, v in diff.items() if k != "largest")
            dev["difference_from_host"] = diff
            dev["speedup_whole"] = row["host"]["whole_s"] / whole  # against the faster host run
            dev["speedup_dipole_step"] = row["host"]["dipole_step_s"] / polar
            dev["speedup_against_blas_threads"] = row["host"]["blas_threads"]
            log(f"   device against host: largest relative difference {diff['largest']:.1e} "
                f"(energy {diff['whole.energy']:.1e}, gradient {diff['whole.gradient']:.1e}, tensor {diff['tensor']:.1e})")
            log(f"   speedup against the host on {row['host']['blas_threads']} thread(s): whole evaluation {dev['speedup_whole']:.1f}x, "
                f"dipole step {dev['speedup_dipole_step']:.1f}x; outside the dipole step: {whole - polar:.3f} s of {whole:.3f} s")
        if not args.no_extras:
            dev["extras"] = extras(cp, dev_res, st["sites"], st["k_vectors"])
            log(f"   extras: {json.dumps(dev['extras'], default=float)}")
        row["device"] = dev
        del dev_res
    except Exception as exc:  # out of memory at this size, or anything else: say so and go on
        row["device"] = {"failed": f"{type(exc).__name__}: {exc}", "traceback": traceback.format_exc()}
        log(f"   device: FAILED at {nc} chains: {type(exc).__name__}: {exc}\n{traceback.format_exc()}")
    finally:
        ew.release()
        ew.xp, ew.device_min_sites = np, 256
        pk._thole_scale_dev = None
        free_pool()
    return row


def table(rows) -> str:
    lines = ["(host: the faster of its runs, with its BLAS threads; speedup: against that run)",
             "chains  sites  k vectors | host: whole  dipole step  threads | device: whole  dipole step  first call  memory | speedup  largest difference",
             "-" * 133]
    for r in rows:
        h, d = r.get("host"), r.get("device")
        left = f"{r['chains']:6d} {r.get('sites', 0):6d} {r.get('k_vectors', 0):10d} | "
        left += f"{h['whole_s']:10.2f} s {h['dipole_step_s']:10.2f} s {h['blas_threads']:8d} | " if h else f"{'-':>12} {'-':>12} {'-':>8} | "
        if d is None:
            left += "no device"
        elif "failed" in d:
            left += d["failed"][:70]
        else:
            left += f"{d['whole_s']:12.3f} s {d['dipole_step_s']:10.3f} s {d['first_call_s']:9.1f} s {d['pool_gb_peak']:6.2f} GB | "
            left += f"{d['speedup_whole']:6.1f}x  {d['difference_from_host']['largest']:.1e}" if "speedup_whole" in d else \
                f"{'-':>7}  (no host run; tiling {d['tiling_check']['difference']:.1e} kcal/mol, gradient {d['gradient_check']['relative_error']:.1e})"
        lines.append(left)
    return "\n".join(lines)


def all_devices(args, n: int) -> list:
    """The same run on every device at once, one process each."""
    procs = []
    for i in range(n):
        out = Path(args.out).with_suffix(f".device{i}.json")
        cmd = [sys.executable, __file__, "--chains", *map(str, args.chains), "--host-max-chains", "0", "--repeat", str(args.repeat),
               "--device", "0", "--no-extras", "--busy-seconds", "0", "--out", str(out)]
        with open(out.with_suffix(".log"), "w") as fh:  # the child keeps its own copy of the handle
            procs.append((i, out, subprocess.Popen(cmd, env={**os.environ, "CUDA_VISIBLE_DEVICES": str(i)},
                                                   stdout=fh, stderr=subprocess.STDOUT)))
    rows = []
    for i, out, p in procs:
        p.wait()
        rows.append({"device": i, "returncode": p.returncode, "rows": json.load(open(out))["rows"] if out.exists() else None})
    return rows


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--chains", type=int, nargs="+", default=[16, 32, 64], choices=sorted(TILINGS))
    ap.add_argument("--host-max-chains", type=int, default=64, help="largest cell the host is timed and compared on (a 64-chain "
                    "host evaluation takes minutes and about 12 GB)")
    ap.add_argument("--repeat", type=int, default=3)
    ap.add_argument("--host-threads", type=int, nargs="+", default=None, help="BLAS thread counts the host is timed at "
                    "(default: 1, the package's own setting; 8; and every core this process is allotted, which in a "
                    "container is the quota and not what the machine has; 0 stands for that last count)")
    ap.add_argument("--device", default="auto", help="'auto', a device number, 'none' for the host columns alone, or 'double' "
                    "for the tests' NumPy-backed stand-in (no speed in it: it checks the script and the port's bookkeeping "
                    "on a machine without a GPU)")
    ap.add_argument("--no-extras", action="store_true", help="skip the two side measurements")
    ap.add_argument("--busy-seconds", type=float, default=8.0, help="how long evaluations run back to back while the card's "
                    "utilisation is sampled (0 to skip)")
    ap.add_argument("--processes", type=int, nargs="+", default=None, metavar="N",
                    help="after the run, for each size: N processes at once on the one card, each evaluating its own cell back to "
                    "back for --busy-seconds, for every N given (1 2 4, say): the card's throughput when it is shared")
    ap.add_argument("--processes-max-chains", type=int, default=32, help="the largest cell --processes is run at (each process "
                    "holds its own cell on the card: 21 GB at sixty-four chains)")
    ap.add_argument("--busy-child", type=float, default=None, help=argparse.SUPPRESS)  # one process of --processes: the common start time
    ap.add_argument("--loop-form", action="store_true", help="also time, and sample the card under, the loop form with the sums "
                    "on the device (the organisation before the batched form): host work with the card nearly idle, which "
                    "is the point of the comparison and the reason it is not run unless asked for")
    ap.add_argument("--all-devices", action="store_true", help="after the run, repeat it on every visible device at once")
    ap.add_argument("--out", default="ewald_device_benchmark.json")
    args = ap.parse_args()
    if args.busy_child is not None:
        busy_child(args)
        return

    from polyfind import backend as bk
    from polyfind.fitting import PRESET_PARAMETERS

    # A mixed cell is evaluated in float64.  Left to itself the package moves its *screening*
    # kernels to a GPU in float32 wherever it finds one; that backend is for the screen, not for
    # this, and the cell refuses a packer built for it.  The device this script is about is the
    # one the Ewald object is given below, which the model's batched form follows.
    bk.set_backend("cpu")
    cores = allotted_cores()["cores"]
    if args.host_threads is None:
        args.host_threads = sorted({1, min(8, cores), cores})
    args.host_threads = sorted({cores if t == 0 else t for t in args.host_threads})  # 0: every allotted core
    cp, why = device_module(args.device)
    result = {"machine": machine(cp), "preset": PRESET, "row": ROW, "batched_kernels": "host, float64",
              "host_blas_threads_timed": args.host_threads, "rows": []}
    log(json.dumps(result["machine"], indent=1))
    if cp is None:
        log(f"\nNo device: {why}.  Host columns only.")
        result["no_device"] = why
    dk = load_drivers()
    with PRESET_PARAMETERS[PRESET].applied():
        for nc in sorted(args.chains):
            result["rows"].append(one_size(dk, nc, cp, args))
            Path(args.out).write_text(json.dumps(result, indent=1, default=float))
    log("\n" + table(result["rows"]))
    if args.processes and cp is not None and hasattr(cp, "cuda"):
        result["several_processes"] = []
        for nc in sorted(c for c in args.chains if c <= args.processes_max_chains):
            for count in args.processes:
                res = several_processes(args, nc, count)
                result["several_processes"].append(res)
                Path(args.out).write_text(json.dumps(result, indent=1, default=float))
                failed = sum("failed" in r for r in res["runs"])
                log(f"   {nc} chains, {count} process{'es' if count > 1 else ''} on the card: "
                    + (f"{res['evaluations_per_second']:.1f} evaluations per second in all, "
                       f"{1e3 * sum(res['per_process_s']) / len(res['per_process_s']):.0f} ms each" if "evaluations_per_second" in res else "none ran")
                    + (f"; card utilisation {res['utilization_mean']:.0f}%" if "utilization_mean" in res else "")
                    + (f"; {failed} FAILED: {next(r['failed'] for r in res['runs'] if 'failed' in r)[-200:]}" if failed else ""))
    if args.all_devices and cp is not None and hasattr(cp, "cuda"):
        n = int(cp.cuda.runtime.getDeviceCount())
        if n > 1:
            log(f"\n== all {n} devices at once, one process each")
            t0 = time.perf_counter()
            result["all_devices"] = {"devices": n, "runs": all_devices(args, n), "wall_s": time.perf_counter() - t0}
            for run in result["all_devices"]["runs"]:
                for r in run["rows"] or []:
                    d = r.get("device", {})
                    log(f"   device {run['device']}: {r['chains']} chains, whole evaluation {d.get('whole_s', float('nan')):.3f} s")
        else:
            log("\n--all-devices: one device visible, nothing to repeat")
    Path(args.out).write_text(json.dumps(result, indent=1, default=float))
    log(f"\nwritten to {args.out}")


if __name__ == "__main__":
    main()
