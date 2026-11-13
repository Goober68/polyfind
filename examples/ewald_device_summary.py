"""Tables of an Ewald device benchmark, from the JSON files examples/ewald_device_benchmark.py wrote.

    python examples/ewald_device_summary.py docs/predictions/ewald_device_a100_20261113/*.json

Prints Markdown: the machine, one row per cell size (the host at each BLAS thread count it was
timed at, the device, what stays on the host, the agreement, the card's memory), the three ported
steps apart, and the side measurements.  Every figure is read from the files; nothing is typed in.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path


def load(paths) -> tuple:
    machine, rows, seen = None, [], set()
    for p in paths:
        d = json.loads(Path(p).read_text())
        if "rows" not in d:
            continue
        machine = machine or d.get("machine")
        for r in d["rows"]:
            if r["chains"] not in seen:
                seen.add(r["chains"])
                rows.append(r)
    return machine, sorted(rows, key=lambda r: r["chains"])


def seconds(x) -> str:
    if x is None:
        return "-"
    return f"{x:.1f} s" if x >= 10 else (f"{x:.2f} s" if x >= 0.1 else f"{1e3 * x:.1f} ms")


def host_run(r, threads=None):
    runs = r.get("host_runs") or []
    if threads is None:
        return min(runs, key=lambda h: h["whole_s"]) if runs else None
    return next((h for h in runs if h["blas_threads"] == threads), None)


def main(paths) -> None:
    machine, rows = load(paths)
    if not rows:
        raise SystemExit("no benchmark rows in " + ", ".join(map(str, paths)))
    dev = (machine or {}).get("device") if isinstance((machine or {}).get("device"), dict) else None
    allotted = (machine or {}).get("cpus_allotted", {})
    print(f"Device: {dev['name'] if dev else '-'}, {dev['memory_gb'] if dev else '-'} GB; {machine.get('device_module', '-')}; "
          f"host {machine.get('cpu', '-')}, {allotted.get('cores', '-')} cores allotted "
          f"(the container sees {allotted.get('affinity', machine.get('cpus_usable', '-'))}); NumPy {machine.get('numpy')}.\n")

    threads = sorted({h["blas_threads"] for r in rows for h in r.get("host_runs", [])})
    head = ["chains", "sites", "k vectors"] + [f"host, {t} thread{'s' if t > 1 else ''}" for t in threads] + \
           ["device", "of it on the host", "device against fastest host", "largest difference", "card memory"]
    print("One whole evaluation (energy and gradient):\n")
    print("| " + " | ".join(head) + " |")
    print("|" + "---|" * len(head))
    for r in rows:
        d = r.get("device") or {}
        cells = [str(r["chains"]), f"{r.get('sites', '-'):,}" if r.get("sites") else "-", f"{r.get('k_vectors', '-'):,}" if r.get("k_vectors") else "-"]
        for t in threads:
            h = host_run(r, t)
            cells.append(seconds(h["whole_s"]) if h else "-")
        if "failed" in d:
            cells += ["failed: " + d["failed"][:60], "-", "-", "-", "-"]
        elif d:
            best = host_run(r)
            left = d["whole_s"] - d["dipole_step_s"]
            diff = d.get("difference_from_host", {}).get("largest")
            cells += [seconds(d["whole_s"]), f"{seconds(left)} ({100 * left / d['whole_s']:.0f}%)",
                      f"{best['whole_s'] / d['whole_s']:.1f}x" if best else "-",
                      f"{diff:.1e}" if diff is not None else "no host run",
                      f"{d['pool_gb_peak']:.1f} GB" if d.get("pool_gb_peak") == d.get("pool_gb_peak") else "-"]
        else:
            cells += ["-"] * 5
        print("| " + " | ".join(cells) + " |")

    print("\nThe three ported steps, timed apart (fastest host run against the device):\n")
    print("| chains | tensor and field | solve | gradient | the three together | device against host |")
    print("|---|---|---|---|---|---|")
    for r in rows:
        d, h = r.get("device") or {}, host_run(r)
        if "tensor_s" not in d:
            continue
        three_d = d["tensor_s"] + d["solve_s"] + d["gradient_s"]
        if h:
            three_h = h["tensor_s"] + h["solve_s"] + h["gradient_s"]
            print(f"| {r['chains']} | {seconds(h['tensor_s'])} / {seconds(d['tensor_s'])} | {seconds(h['solve_s'])} / {seconds(d['solve_s'])} | "
                  f"{seconds(h['gradient_s'])} / {seconds(d['gradient_s'])} | {seconds(three_h)} / {seconds(three_d)} | {three_h / three_d:.0f}x |")
        else:
            print(f"| {r['chains']} | - / {seconds(d['tensor_s'])} | - / {seconds(d['solve_s'])} | - / {seconds(d['gradient_s'])} | - / {seconds(three_d)} | - |")

    print("\nChecks that need no host run, and the first call:\n")
    print("| chains | tiling against the stored two-chain cell | gradient against the device's own energy difference | first device evaluation | device memory in use |")
    print("|---|---|---|---|---|")
    for r in rows:
        d = r.get("device") or {}
        if "whole_s" not in d:
            continue
        print(f"| {r['chains']} | {d['tiling_check']['difference']:.1e} kcal/mol | {d['gradient_check']['relative_error']:.1e} | "
              f"{seconds(d['first_call_s'])} | {d.get('device_gb_in_use', float('nan')):.1f} of {d.get('device_gb_total', float('nan')):.1f} GB |")

    print("\nSide measurements (not used by the port):\n")
    print("| chains | general product / rank-k update | double-precision rate of the product | LU solve / Cholesky solve | unknowns |")
    print("|---|---|---|---|---|")
    for r in rows:
        x = (r.get("device") or {}).get("extras") or {}
        sp, so = x.get("symmetric_product"), x.get("solve")
        if not isinstance(sp, dict) and not isinstance(so, dict):
            continue
        a = f"{seconds(sp['gemm_s'])} / {seconds(sp['syrk_s'])}" if isinstance(sp, dict) else str(sp)[:40]
        b = f"{sp['gflops_gemm'] / 1e3:.1f} TFLOP/s" if isinstance(sp, dict) else "-"
        c = f"{seconds(so['lu_s'])} / {seconds(so['cholesky_s'])}" if isinstance(so, dict) else str(so)[:40]
        print(f"| {r['chains']} | {a} | {b} | {c} | {so['unknowns'] if isinstance(so, dict) else '-':,} |" if isinstance(so, dict)
              else f"| {r['chains']} | {a} | {b} | {c} | - |")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    main(sys.argv[1:])
