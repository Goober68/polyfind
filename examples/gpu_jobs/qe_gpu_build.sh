#!/usr/bin/env bash
# Quantum ESPRESSO built for the pod's card: pw.x with OpenACC and CUDA, by NVIDIA's compiler.
# A job of examples/runpod_jobs_pod.sh; tried on a pod only (it needs a card and 25 GB of disk).
#
#   qe_gpu_build.sh [JOBS]       JOBS: make's parallel jobs (default 3: the processors are DFT's)
#
# The processors' pw.x (conda-forge, examples/gpu_jobs/dft_cells.sh) takes an hour for a
# 168-atom cell and leaves the card idle.  This fetches the NVIDIA HPC SDK (NVHPC_PACKAGE,
# default nvhpc-25-7, from NVIDIA's apt repository: 5.6 GB) and Quantum ESPRESSO 7.5's source
# (the version the labels were computed with), configures as the source's README_GPU.md says and
# makes pw.  Into /opt/qe_gpu/bin go two small scripts, pw.x and mpirun, that set the compiler's
# libraries and call the real ones, so that examples/mlip_frames.py finds an mpirun beside pw.x
# as it does in the conda installation.  Nothing is believed until qe_gpu_check.sh has run.
#
# Into $JOB_OUT: versions.txt, configure.log, make.log (their ends are in log.txt), built.txt
# when pw.x exists.  The job fails when it does not.
set -euo pipefail
JOBS="${1:-3}"
PREFIX="${QE_GPU_PREFIX:-/opt/qe_gpu}"
PACKAGE="${NVHPC_PACKAGE:-nvhpc-25-7}"
QE_TAG="${QE_GPU_TAG:-qe-7.5}"
mkdir -p "$PREFIX/bin"
t0=$(date +%s)
nvidia-smi --query-gpu=name,compute_cap,driver_version,memory.total --format=csv,noheader | tee "$JOB_OUT/versions.txt"
CC_CAP=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader | head -n 1 | tr -d '. ')

export DEBIAN_FRONTEND=noninteractive
if ! ls /opt/nvidia/hpc_sdk/Linux_x86_64/2*/compilers/bin/nvfortran >/dev/null 2>&1; then
  curl -fsSL https://developer.download.nvidia.com/hpc-sdk/ubuntu/DEB-GPG-KEY-NVIDIA-HPC-SDK | gpg --dearmor -o /usr/share/keyrings/nvidia-hpcsdk-archive-keyring.gpg
  echo 'deb [signed-by=/usr/share/keyrings/nvidia-hpcsdk-archive-keyring.gpg] https://developer.download.nvidia.com/hpc-sdk/ubuntu/amd64 /' > /etc/apt/sources.list.d/nvhpc.list
  apt-get update -y -q > "$JOB_OUT/apt.log" 2>&1
  apt-get install -y -q --no-install-recommends "$PACKAGE" wget >> "$JOB_OUT/apt.log" 2>&1 || { tail -n 30 "$JOB_OUT/apt.log"; echo "the compiler did not install"; exit 1; }
  apt-get clean
fi
echo "compiler installed in $(( $(date +%s) - t0 )) s"
NV=$(ls -d /opt/nvidia/hpc_sdk/Linux_x86_64/2*/ | sort -V | tail -n 1); NV="${NV%/}"
CUDA_DIR=$(ls -d "$NV"/cuda/[0-9]* | sort -V | tail -n 1)
CUDA_RT=$(basename "$CUDA_DIR")
MPI=""
for d in "$NV/comm_libs/mpi" "$NV"/comm_libs/openmpi4/openmpi-* "$NV"/comm_libs/openmpi/openmpi-* "$NV"/comm_libs/*/openmpi4/openmpi-*; do
  if [ -x "$d/bin/mpif90" ] && [ -x "$d/bin/mpirun" ]; then MPI="$d"; break; fi
done
[ -n "$MPI" ] || { echo "no MPI with mpif90 and mpirun under $NV/comm_libs:"; ls "$NV/comm_libs"; exit 1; }
{ echo "sdk $NV"; echo "cuda toolkit $CUDA_DIR (runtime $CUDA_RT), compute capability $CC_CAP"; echo "mpi $MPI"; "$NV/compilers/bin/nvfortran" --version | head -n 2; } | tee -a "$JOB_OUT/versions.txt"

cat > "$PREFIX/env.sh" <<ENV
export PATH="$MPI/bin:$NV/compilers/bin:\$PATH"
export LD_LIBRARY_PATH="$NV/compilers/lib:$NV/math_libs/lib64:$CUDA_DIR/lib64:$MPI/lib:\${LD_LIBRARY_PATH:-}"
ENV
source "$PREFIX/env.sh"

SRC="$PREFIX/q-e-$QE_TAG"
if [ ! -d "$SRC" ]; then
  curl -fsSL "https://gitlab.com/QEF/q-e/-/archive/$QE_TAG/q-e-$QE_TAG.tar.gz" | tar -xz -C "$PREFIX"
fi
cd "$SRC"
if ! ./configure --with-cuda="$CUDA_DIR" --with-cuda-runtime="$CUDA_RT" --with-cuda-cc="$CC_CAP" --enable-openmp --with-scalapack=no \
     MPIF90=mpif90 F90=nvfortran CC=nvc > "$JOB_OUT/configure.log" 2>&1; then
  tail -n 40 "$JOB_OUT/configure.log"; echo "configure failed"; exit 1
fi
tail -n 25 "$JOB_OUT/configure.log"
grep -E "^(DFLAGS|MPIF90|F90|CUDA_|LD_LIBS|BLAS_LIBS|FFT_LIBS|FFLAGS) " make.inc | tee -a "$JOB_OUT/versions.txt" || true
t1=$(date +%s)
if ! make -j"$JOBS" pw > "$JOB_OUT/make.log" 2>&1; then
  grep -nE "Error|error:|NVFORTRAN-[SF]-|undefined reference" "$JOB_OUT/make.log" | tail -n 30 || true
  tail -n 30 "$JOB_OUT/make.log"; echo "make pw failed after $(( $(date +%s) - t1 )) s"; exit 1
fi
[ -x "$SRC/bin/pw.x" ] || { tail -n 30 "$JOB_OUT/make.log"; echo "make ended and there is no pw.x"; exit 1; }
echo "pw.x made in $(( $(date +%s) - t1 )) s"

cat > "$PREFIX/bin/pw.x" <<WRAP
#!/usr/bin/env bash
source "$PREFIX/env.sh"
exec "$SRC/bin/pw.x" "\$@"
WRAP
cat > "$PREFIX/bin/mpirun" <<WRAP
#!/usr/bin/env bash
source "$PREFIX/env.sh"
exec "$MPI/bin/mpirun" "\$@"
WRAP
chmod +x "$PREFIX/bin/pw.x" "$PREFIX/bin/mpirun"
ldd "$SRC/bin/pw.x" | grep -iE "not found|cuda|cufft|cusolver|cublas|acc" | head -n 20 || true
echo "$PREFIX/bin/pw.x built in $(( $(date +%s) - t0 )) s in all" | tee "$JOB_OUT/built.txt"
