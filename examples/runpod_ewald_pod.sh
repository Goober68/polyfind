#!/usr/bin/env bash
# examples/runpod_ewald_bench.sh on a RunPod pod that nobody is logged in to.
#
# Given to the pod as its command:
#
#   bash -c "curl -sL https://raw.githubusercontent.com/Goober68/polyfind/claude/polymeric-stable-arrangements-uh06b1/examples/runpod_ewald_pod.sh | bash"
#
# with port 8000 exposed as an HTTP port.  While the benchmark runs, and for a while after, its
# log and its JSON files can be read at
#
#   https://<pod id>-8000.proxy.runpod.net/<BENCH_TOKEN>/
#
# (that address is public to whoever has it, so the token is the lock: make it long and random;
# the root of the server lists nothing).  The pod removes itself: BENCH_LINGER_S after the
# benchmark has ended, and after BENCH_LIMIT_S whatever has happened, so that a pod whose
# watcher has gone away does not go on being billed.
#
#   BENCH_TOKEN      required: the path the results are served under
#   BENCH_LIMIT_S    the pod's longest life in seconds (default 10800, three hours)
#   BENCH_LINGER_S   how long the results stay up after the run (default 2700)
#   BENCH_KEEP=1     do not remove the pod
#   BENCH_PORT       default 8000
#   POLYFIND_BRANCH and the EWALD_BENCH_* variables go on to the benchmark script
set -u

TOKEN="${BENCH_TOKEN:?BENCH_TOKEN is not set}"
PORT="${BENCH_PORT:-8000}"
ROOT="${BENCH_ROOT:-/srv/bench}"
BRANCH="${POLYFIND_BRANCH:-claude/polymeric-stable-arrangements-uh06b1}"
STATUS="$ROOT/$TOKEN/STATUS"
PY="$(command -v python3 || command -v python)"

mkdir -p "$ROOT/$TOKEN"
: > "$ROOT/index.html"   # so that the root shows an empty page and not the token
(cd "$ROOT" && exec "$PY" -m http.server "$PORT" --bind 0.0.0.0 >/dev/null 2>&1) &
SERVER=$!

remove_pod() {
  [ "${BENCH_KEEP:-0}" = "1" ] && return 0
  [ -n "${RUNPOD_POD_ID:-}" ] || return 0
  # With the pod's own key.  The command line tool has spelled this three ways over its
  # versions, and an image carries whichever it was built with; then the API itself; and if
  # the pod cannot be removed, stopping it at least ends the billing of the GPU.
  runpodctl config --apiKey "${RUNPOD_API_KEY:-}" >/dev/null 2>&1
  for attempt in "runpodctl pod remove" "runpodctl pod delete" "runpodctl remove pod"; do
    if $attempt "$RUNPOD_POD_ID" >/dev/null 2>&1; then echo "   removed with: $attempt" >> "$STATUS"; return 0; fi
  done
  code=$(curl -s -o /dev/null -w "%{http_code}" -X DELETE "https://rest.runpod.io/v1/pods/$RUNPOD_POD_ID" \
           -H "Authorization: Bearer ${RUNPOD_API_KEY:-}")
  echo "   command line removal failed; API delete answered $code" >> "$STATUS"
  case "$code" in 2*) return 0 ;; esac
  for attempt in "runpodctl pod stop" "runpodctl stop pod"; do
    if $attempt "$RUNPOD_POD_ID" >/dev/null 2>&1; then echo "   could not remove; stopped with: $attempt" >> "$STATUS"; return 0; fi
  done
  echo "   THE POD COULD NOT REMOVE OR STOP ITSELF: remove it by hand" >> "$STATUS"
  return 1
}

echo "started $(date -u +%FT%TZ) pod ${RUNPOD_POD_ID:-none}" > "$STATUS"
( sleep "${BENCH_LIMIT_S:-10800}"; echo "time limit reached $(date -u +%FT%TZ): removing the pod" >> "$STATUS"; remove_pod ) &

cd "${BENCH_WORKDIR:-/root}" || exit 1
export EWALD_BENCH_OUT="$ROOT/$TOKEN/run"
if [ -n "${BENCH_SCRIPT:-}" ]; then
  cp "$BENCH_SCRIPT" runpod_ewald_bench.sh
else
  curl -sL "https://raw.githubusercontent.com/Goober68/polyfind/$BRANCH/examples/runpod_ewald_bench.sh" -o runpod_ewald_bench.sh
fi
bash runpod_ewald_bench.sh
CODE=$?
echo "finished $(date -u +%FT%TZ) exit $CODE" >> "$STATUS"

sleep "${BENCH_LINGER_S:-2700}"
echo "removing the pod $(date -u +%FT%TZ)" >> "$STATUS"
remove_pod
if [ -n "${RUNPOD_POD_ID:-}" ]; then
  sleep infinity   # a command that ends is started again: if the removal did not take, do not run it all twice
fi
kill "$SERVER" 2>/dev/null
exit "$CODE"
