#!/usr/bin/env bash
# The potential fitted on from a model of ours on the card (examples/gpu_jobs/mlip_fit.sh, whose
# arguments these are), while the processors, which a fit hardly uses, build the DFT program for
# the card (examples/gpu_jobs/qe_gpu_build.sh; 39 minutes when it was first done).  A job of
# examples/runpod_jobs_pod.sh:
#
#   refit.sh NAME EPOCHS SEED [frames FILE...] [prepare ARGUMENT...] [train ARGUMENT...]
#
# Into $JOB_OUT: what mlip_fit.sh leaves (NAME_sSEED.model, its log and results) and the build's
# log (qe_gpu_build.log, with the build's own files).  The fit's failure is the job's; the
# build's is not (the cells are then labelled on the processors alone).  REFIT_NO_BUILD=1 fits
# without building.
set -uo pipefail
HERE="$(dirname "$0")"
if [ "${REFIT_NO_BUILD:-0}" != "1" ] && [ ! -x /opt/qe_gpu/bin/pw.x ]; then
  bash "$HERE/qe_gpu_build.sh" "${REFIT_BUILD_JOBS:-5}" > "$JOB_OUT/qe_gpu_build.log" 2>&1 &
  BUILD=$!
fi
bash "$HERE/mlip_fit.sh" "$@"
code=$?
if [ -n "${BUILD:-}" ]; then
  wait "$BUILD" && echo "the DFT program for the card is built" || echo "THE BUILD FOR THE CARD FAILED (qe_gpu_build.log): DFT is the processors' alone"
fi
exit $code
