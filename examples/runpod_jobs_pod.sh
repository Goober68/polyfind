#!/usr/bin/env bash
# A RunPod pod that works through a queue of jobs kept in the repository, with nobody logged in.
#
# Given to the pod as its command:
#
#   bash -c "curl -sL https://raw.githubusercontent.com/Goober68/polyfind/<commit>/examples/runpod_jobs_pod.sh | bash"
#
# with port 8000 exposed as an HTTP port and BENCH_TOKEN and BENCH_QUEUE set.
#
# The queue is a text file of the branch (BENCH_QUEUE: its path in the repository), one job a line:
#
#   <script under examples/gpu_jobs/> [arguments]
#
# The pod installs the package once (examples/runpod_ewald_bench.sh, set-up only), then follows
# the branch: every ten seconds it fetches it and runs the next line of the queue it has not run,
# with the code of the commit that is the branch's head at that moment.  So a job is started --
# and one that failed is corrected and started again, as a new line -- by a push.  Nothing is sent
# to the pod and nothing on it listens for commands; the queue is append-only (the pod counts
# lines).  A line `stop` ends the run.  Empty lines and lines starting with # are not jobs.
#
# Job N writes into <served>/jobs/NN_<script>/ : its output as log.txt and whatever it puts in
# $JOB_OUT.  The served directory is at
#
#   https://<pod id>-8000.proxy.runpod.net/<BENCH_TOKEN>/
#
# (public to whoever has the address, so the token is the lock; the root lists nothing), and
# STATUS there has one line per event.  The pod removes itself: BENCH_LINGER_S after `stop`,
# after BENCH_IDLE_S with nothing new in the queue, and after BENCH_LIMIT_S whatever has happened.
#
#   BENCH_TOKEN      required: the path the results are served under
#   BENCH_QUEUE      required: the queue file, e.g. examples/gpu_jobs/queue_20261113.txt
#   BENCH_LIMIT_S    the pod's longest life in seconds (default 7200)
#   BENCH_IDLE_S     how long it waits for a new line before it removes itself (default 900)
#   BENCH_LINGER_S   how long the results stay up after `stop` (default 300)
#   BENCH_KEEP=1     do not remove the pod
#   BENCH_PORT       default 8000
#   POLYFIND_BRANCH  the branch followed (default the working branch)
#
# A job sees: PY (the Python of the installation), JOB_OUT (its directory), POLYFIND_DEVICE=cpu and
# POLYFIND_EWALD_DEVICE=cuda (the large cells' evaluation on the card), and runs in the repository.
set -u

TOKEN="${BENCH_TOKEN:?BENCH_TOKEN is not set}"
QUEUE="${BENCH_QUEUE:?BENCH_QUEUE is not set}"
PORT="${BENCH_PORT:-8000}"
ROOT="${BENCH_ROOT:-/srv/bench}"
BRANCH="${POLYFIND_BRANCH:-claude/polymeric-stable-arrangements-uh06b1}"
SERVED="$ROOT/$TOKEN"
STATUS="$SERVED/STATUS"
PY0="$(command -v python3 || command -v python)"

mkdir -p "$SERVED/jobs"
: > "$ROOT/index.html"   # so that the root shows an empty page and not the token
(cd "$ROOT" && exec "$PY0" -m http.server "$PORT" --bind 0.0.0.0 >/dev/null 2>&1) &
SERVER=$!
say() { echo "$(date -u +%FT%TZ) $*" >> "$STATUS"; }

remove_pod() {
  [ "${BENCH_KEEP:-0}" = "1" ] && return 0
  [ -n "${RUNPOD_POD_ID:-}" ] || return 0
  runpodctl config --apiKey "${RUNPOD_API_KEY:-}" >/dev/null 2>&1
  for attempt in "runpodctl pod remove" "runpodctl pod delete" "runpodctl remove pod"; do
    if $attempt "$RUNPOD_POD_ID" >/dev/null 2>&1; then say "removed with: $attempt"; return 0; fi
  done
  code=$(curl -s -o /dev/null -w "%{http_code}" -X DELETE "https://rest.runpod.io/v1/pods/$RUNPOD_POD_ID" \
           -H "Authorization: Bearer ${RUNPOD_API_KEY:-}")
  say "command line removal failed; API delete answered $code"
  case "$code" in 2*) return 0 ;; esac
  for attempt in "runpodctl pod stop" "runpodctl stop pod"; do
    if $attempt "$RUNPOD_POD_ID" >/dev/null 2>&1; then say "could not remove; stopped with: $attempt"; return 0; fi
  done
  say "THE POD COULD NOT REMOVE OR STOP ITSELF: remove it by hand"
  return 1
}

say "started pod ${RUNPOD_POD_ID:-none}, queue $QUEUE of $BRANCH"
( sleep "${BENCH_LIMIT_S:-7200}"; say "time limit reached: removing the pod"; remove_pod ) &

cd "${BENCH_WORKDIR:-/root}" || exit 1
if [ -n "${BENCH_SCRIPT_DIR:-}" ]; then   # a local copy of the tree, for trying this script without a pod
  cp "$BENCH_SCRIPT_DIR/examples/runpod_ewald_bench.sh" runpod_ewald_bench.sh
else
  curl -sL "https://raw.githubusercontent.com/Goober68/polyfind/$BRANCH/examples/runpod_ewald_bench.sh" -o runpod_ewald_bench.sh
fi
if EWALD_BENCH_SETUP_ONLY=1 EWALD_BENCH_OUT="$SERVED/setup" bash runpod_ewald_bench.sh; then
  say "set up"
else
  say "SET-UP FAILED (see setup/log.txt); finished exit 1"
  sleep "${BENCH_LINGER_S:-300}"; remove_pod
  [ -n "${RUNPOD_POD_ID:-}" ] && sleep infinity
  exit 1
fi
[ -d polyfind ] && cd polyfind
PY="$(cat .ewald_bench_python)"
export PY POLYFIND_DEVICE=cpu POLYFIND_EWALD_DEVICE="${BENCH_EWALD_DEVICE:-cuda}" POLYFIND_BRANCH="$BRANCH"

follow() {  # bring the tree to the branch's head; a fetch that fails is tried again next time round
  [ -n "${BENCH_NO_FETCH:-}" ] && return 0
  git fetch -q --depth 1 origin "$BRANCH" 2>/dev/null && git reset -q --hard FETCH_HEAD
}

done_lines=0
last=$(date +%s)
while true; do
  follow
  # the queue's job lines, in order
  mapfile -t LINES < <(grep -vE '^[[:space:]]*(#|$)' "$QUEUE" 2>/dev/null)
  if [ "${#LINES[@]}" -gt "$done_lines" ]; then
    line="${LINES[$done_lines]}"
    done_lines=$((done_lines + 1))
    last=$(date +%s)
    if [ "$line" = "stop" ]; then say "stop: finished exit 0"; break; fi
    read -r -a WORDS <<< "$line"
    script="examples/gpu_jobs/${WORDS[0]}"
    name="$(printf '%02d' "$done_lines")_$(basename "${WORDS[0]}" .sh)"
    export JOB_OUT="$SERVED/jobs/$name"
    mkdir -p "$JOB_OUT"
    if [ ! -f "$script" ]; then say "job $done_lines: NO SUCH SCRIPT $script"; continue; fi
    say "job $done_lines started: $line (commit $(git rev-parse --short HEAD 2>/dev/null))"
    bash "$script" "${WORDS[@]:1}" > "$JOB_OUT/log.txt" 2>&1
    say "job $done_lines ended: exit $? ($name)"
    last=$(date +%s)
  else
    if [ $(( $(date +%s) - last )) -ge "${BENCH_IDLE_S:-900}" ]; then say "nothing new in the queue for ${BENCH_IDLE_S:-900} s: finished exit 0"; break; fi
    sleep 10
  fi
done

sleep "${BENCH_LINGER_S:-300}"
say "removing the pod"
remove_pod
if [ -n "${RUNPOD_POD_ID:-}" ]; then
  sleep infinity   # a command that ends is started again: if the removal did not take, do not run it all twice
fi
kill "$SERVER" 2>/dev/null
exit 0
