#!/usr/bin/env bash
# The device tests, then the benchmark on the card alone: every size's evaluation part by part and
# the card's utilisation, and several processes sharing the card at the sizes one process leaves it
# partly idle.  A job of examples/runpod_jobs_pod.sh; arguments: the chain counts (default 16 32 64),
# then, after a literal "large", one more size run by itself (it may not fit the card).
set -uo pipefail
SIZES=()
LARGE=""
while [ $# -gt 0 ]; do
  if [ "$1" = "large" ]; then LARGE="${2:-}"; break; fi
  SIZES+=("$1"); shift
done
[ ${#SIZES[@]} -gt 0 ] || SIZES=(16 32 64)
DEVICE="${BENCH_DEVICE:-auto}"
# (the tests choose their own array modules and start from the package's default, the host:
# the switch the pod sets for its jobs is taken away for them)
env -u POLYFIND_EWALD_DEVICE "$PY" -m pytest tests/test_ewald_device.py tests/test_mixed_batched.py tests/test_gpu.py -q \
  -k "ewald or double or device or dipole or batched or mixed or gram" -p no:cacheprovider || exit 1
"$PY" examples/ewald_device_benchmark.py --device "$DEVICE" --chains "${SIZES[@]}" --host-max-chains 0 \
  --processes ${BENCH_PROCESSES:-1 2 3 4} --processes-max-chains "${BENCH_PROCESSES_MAX:-32}" \
  --out "$JOB_OUT/ewald_device_benchmark.json" || exit 1
if [ -n "$LARGE" ]; then
  "$PY" examples/ewald_device_benchmark.py --device "$DEVICE" --chains "$LARGE" --host-max-chains 0 \
    --out "$JOB_OUT/ewald_device_benchmark_$LARGE.json" || echo "== the $LARGE-chain run did not finish (exit $?)"
fi
"$PY" examples/ewald_device_summary.py "$JOB_OUT"/ewald_device_benchmark*.json > "$JOB_OUT/summary.md" || true
