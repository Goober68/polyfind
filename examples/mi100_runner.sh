#!/usr/bin/env bash
# The PI's own machine with an AMD card (an Instinct MI100) working through a queue of jobs kept in
# this repository, as the rented pods did (examples/runpod_jobs_pod.sh), with nobody logged in.
#
# One command on that machine, either
#
#   curl -sL https://raw.githubusercontent.com/Goober68/polyfind/claude/polymeric-stable-arrangements-uh06b1/examples/mi100_runner.sh | MI100_FOLLOW=1 bash
#
# or, from a checkout of this branch there,
#
#   MI100_FOLLOW=1 bash examples/mi100_runner.sh
#
# (under tmux, screen or nohup if the terminal will be closed).  It needs Linux on x86_64 with curl,
# tar and git, about 12 GB of disk, and no root: everything goes under MI100_ROOT.
#
# WHAT IT DOES TO THE MACHINE.  It makes its own clone of the repository and its own programs under
# MI100_ROOT (examples/mi100_setup.sh) and runs the lines of the queue file in order, one job a
# line:
#
#   <script under examples/gpu_jobs/> [arguments]
#
# With MI100_FOLLOW=1 it then keeps following the branch: every 30 seconds it fetches it and runs
# the next line it has not run, with the code of the branch's head at that moment.  A job is started
# -- and one that failed is corrected and started again, as a new line -- by a push to the branch.
# THAT IS REMOTE CONTROL OF THIS MACHINE BY WHOEVER CAN PUSH TO THE BRANCH: a pushed line runs as
# the user who started this script.  Without MI100_FOLLOW it runs the queue as it stands when
# started and ends.  A line `stop` ends it either way; so does Ctrl-C, or killing it.  Empty lines
# and lines starting with # are not jobs; the queue is append-only (lines are counted).
#
# WHAT COMES BACK.  Job N writes into MI100_ROOT/results/jobs/NN_<script>/ : its output as log.txt
# and whatever it puts in $JOB_OUT (small files: a job keeps its bulk under MI100_ROOT/triage).
# After the set-up and after every job the results folder is committed and pushed to the branch
# MI100_RESULTS_BRANCH of the same repository, with whatever right to push git has on this machine.
# If the push is refused the results stay in the folder and STATUS says so; nothing else is pushed
# anywhere, and no file over 20 MB is committed.  The repository is public, so the machine's name
# is not written and the home folder's path is replaced by ~ in the logs before they are committed.
#
#   MI100_ROOT             everything it makes (default ~/polyfind_mi100)
#   MI100_FOLLOW=1         keep following the branch (see above)
#   MI100_IDLE_S           with MI100_FOLLOW: end after this long with nothing new (default 86400)
#   MI100_QUEUE            the queue file (default examples/gpu_jobs/queue_mi100.txt)
#   MI100_REPO             the repository's address (default: the checkout's own origin, or
#                          https://github.com/Goober68/polyfind.git when piped from curl)
#   MI100_RESULTS_BRANCH   default mi100-results
#   MI100_NO_PUSH=1        keep the results on the machine
#   MI100_NP               processes for DFT (default the physical cores)
#   POLYFIND_BRANCH        the branch followed (default the working branch)
#
# A job sees what MI100_ROOT/env.sh exports (PY, TRIAGE_PW, TRIAGE_PSEUDO, TRIAGE_PSEUDO_NC,
# TRIAGE_MODEL, TRIAGE_NP, TRIAGE_FRAMES) and JOB_OUT, and runs in the clone.
set -u

ROOT="${MI100_ROOT:-$HOME/polyfind_mi100}"
BRANCH="${POLYFIND_BRANCH:-claude/polymeric-stable-arrangements-uh06b1}"
QUEUE="${MI100_QUEUE:-examples/gpu_jobs/queue_mi100.txt}"
RESULTS_BRANCH="${MI100_RESULTS_BRANCH:-mi100-results}"
REPO="${MI100_REPO:-}"
if [ -z "$REPO" ] && [ -n "${BASH_SOURCE[0]:-}" ] && [ -f "${BASH_SOURCE[0]}" ]; then
  REPO="$(git -C "$(dirname "${BASH_SOURCE[0]}")" config --get remote.origin.url 2>/dev/null || true)"
fi
REPO="${REPO:-https://github.com/Goober68/polyfind.git}"
CODE="$ROOT/polyfind"
RESULTS="$ROOT/results"
STATUS="$RESULTS/STATUS"

mkdir -p "$ROOT"
exec 9> "$ROOT/runner.lock"
if command -v flock >/dev/null && ! flock -n 9; then
  echo "another runner already has $ROOT (its lock is held): not starting a second one"
  exit 1
fi
for tool in curl tar git; do command -v "$tool" >/dev/null || { echo "no $tool on this machine"; exit 1; }; done

# the code: a clone of its own, so that a checkout the PI works in is never reset
if [ ! -d "$CODE/.git" ]; then
  git clone -q --depth 1 --branch "$BRANCH" "$REPO" "$CODE" || { echo "could not clone $REPO"; exit 1; }
fi
follow() {  # bring the clone to the branch's head; a fetch that fails is tried again next time round
  git -C "$CODE" fetch -q --depth 1 origin "$BRANCH" 2>/dev/null && git -C "$CODE" reset -q --hard FETCH_HEAD
}
follow

# the results: a second clone, on a branch of its own that holds nothing but results
if [ ! -d "$RESULTS/.git" ]; then
  if git ls-remote --exit-code --heads "$REPO" "$RESULTS_BRANCH" >/dev/null 2>&1; then
    git clone -q --depth 1 --branch "$RESULTS_BRANCH" "$REPO" "$RESULTS"
  else
    git init -q "$RESULTS"
    git -C "$RESULTS" remote add origin "$REPO"
    git -C "$RESULTS" checkout -q --orphan "$RESULTS_BRANCH"
  fi
fi
# the results' commits carry the runner's name, not whoever's git is set up on this machine
git -C "$RESULTS" config user.name "mi100 runner"
git -C "$RESULTS" config user.email "mi100-runner@localhost"
mkdir -p "$RESULTS/jobs" "$RESULTS/setup"
say() { echo "$(date -u +%FT%TZ) $*" | tee -a "$STATUS"; }

publish() {  # commit the results folder and push it; a refusal is said once a state, not hidden
  ( cd "$RESULTS" || exit 0
    # the repository is public: the home folder's name (a user's name) is taken out of what is committed
    find . -path ./.git -prune -o -type f \( -name '*.txt' -o -name STATUS \) -size -20M -print0 | xargs -0 -r sed -i "s|$HOME|~|g"
    find . -path ./.git -prune -o -type f -size +20M -print | sed 's|^\./||' > .git/info/exclude
    git add -A >/dev/null 2>&1
    git diff --cached --quiet && exit 0
    git commit -q -m "$1" >/dev/null 2>&1 )
  [ -n "${MI100_NO_PUSH:-}" ] && return 0
  if git -C "$RESULTS" push -q origin "HEAD:$RESULTS_BRANCH" 2>"$ROOT/push.err"; then
    PUSH_FAILED=""
  elif [ -z "${PUSH_FAILED:-}" ]; then
    PUSH_FAILED=1
    say "THE RESULTS COULD NOT BE PUSHED ($(tail -n 1 "$ROOT/push.err" | cut -c1-160)): they are in $RESULTS"
  fi
}

say "started, queue $QUEUE of $BRANCH at $(git -C "$CODE" rev-parse --short HEAD), follow ${MI100_FOLLOW:-0}"
if MI100_ROOT="$ROOT" bash "$CODE/examples/mi100_setup.sh" > "$RESULTS/setup/log.txt" 2>&1; then
  say "set up"
else
  say "SET-UP FAILED (setup/log.txt: $(tail -n 1 "$RESULTS/setup/log.txt" | cut -c1-200))"
  publish "set-up failed"
  exit 1
fi
publish "set up"
# shellcheck disable=SC1091
. "$ROOT/env.sh"
cd "$CODE" || exit 1

done_lines=0
last=$(date +%s)
while true; do
  [ -n "${MI100_FOLLOW:-}" ] && follow
  mapfile -t LINES < <(grep -vE '^[[:space:]]*(#|$)' "$QUEUE" 2>/dev/null)
  if [ "${#LINES[@]}" -gt "$done_lines" ]; then
    line="${LINES[$done_lines]}"
    done_lines=$((done_lines + 1))
    if [ "$line" = "stop" ]; then say "stop: finished"; break; fi
    read -r -a WORDS <<< "$line"
    script="examples/gpu_jobs/${WORDS[0]}"
    name="$(printf '%02d' "$done_lines")_$(basename "${WORDS[0]}" .sh)"
    export JOB_OUT="$RESULTS/jobs/$name"
    mkdir -p "$JOB_OUT"
    if [ ! -f "$script" ]; then say "job $done_lines: NO SUCH SCRIPT $script"; publish "job $done_lines: no such script"; continue; fi
    say "job $done_lines started: $line (commit $(git rev-parse --short HEAD 2>/dev/null))"
    publish "job $done_lines started"
    bash "$script" "${WORDS[@]:1}" > "$JOB_OUT/log.txt" 2>&1
    say "job $done_lines ended: exit $? ($name)"
    publish "job $done_lines ended: $line"
    last=$(date +%s)
  elif [ -z "${MI100_FOLLOW:-}" ]; then
    say "the queue is done: finished (start with MI100_FOLLOW=1 to have it wait for more)"
    break
  elif [ $(( $(date +%s) - last )) -ge "${MI100_IDLE_S:-86400}" ]; then
    say "nothing new in the queue for ${MI100_IDLE_S:-86400} s: finished"
    break
  else
    sleep 30
  fi
done
publish "finished"
exit 0
