#!/bin/bash
#SBATCH --job-name=senga2_512_8n
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=128
#SBATCH --time=02:30:00
#SBATCH --partition=standard
#SBATCH --qos=standard
#SBATCH --account=e647
#SBATCH --mem=0
#SBATCH --output=results/slurm-%j.out
#SBATCH --error=results/slurm-%j.err

# 512^3 / 10-step / 8-node (1024 rank) tiling + ops_execute sweep.
#
# BEFORE submitting, rebuild the tiled binary at 512^3:
#   source /work/e647/e647/ireguly5/source_archer2_gcc
#   cd /work/e647/e647/ireguly5/SENGA2_OPS/src
#   make senga2_mpi_tiled SENGA_NX=512 -j8
#
# Then from the app root (this script's cwd on Archer2):
#   mkdir -p results
#   sbatch scripts/submit_archer2_8node_512.sh
#
# After the job, copy back ONE tarball:
#   scp archer2: .../results/archer2_8node_512_results.tgz .
#   tar xzf archer2_8node_512_results.tgz
#   python3 scripts/collate_senga_logs.py results

set -u

source /work/e647/e647/ireguly5/source_archer2_gcc

# Change if this clone lives somewhere else.
ROOT="${SLURM_SUBMIT_DIR:-/work/e647/e647/ireguly5/SENGA2_OPS}"
cd "$ROOT"

BIN="$ROOT/senga2_mpi_tiled"
if [[ ! -x "$BIN" && -x "$ROOT/src/senga2_mpi_tiled" ]]; then
  BIN="$ROOT/src/senga2_mpi_tiled"
fi
if [[ ! -x "$BIN" ]]; then
  echo "ERROR: senga2_mpi_tiled not found. Build with: make senga2_mpi_tiled SENGA_NX=512"
  exit 1
fi

mkdir -p results
export OMP_NUM_THREADS=1
# HDF5 dumps of 512^3 are tens of GB; leave off for the sweep.
unset SENGA_DUMP_LAST || true

# Install the 512^3 / 10-step control file for the duration of the job.
if [[ ! -f input/cont_512_10.dat ]]; then
  echo "ERROR: input/cont_512_10.dat missing"
  exit 1
fi
cp -f input/cont.dat "results/cont.dat.bak.$$"
cp -f input/cont_512_10.dat input/cont.dat

restore_cont() {
  if [[ -f "results/cont.dat.bak.$$" ]]; then
    mv -f "results/cont.dat.bak.$$" input/cont.dat
  fi
}
trap restore_cont EXIT

SRUN="srun --hint=nomultithread --distribution=block:block --unbuffered"

# name | SENGA_OPS_EXEC (empty = compiled default, ispec off) | extra OPS args
# Tiling axis first (default exec), then exec-site axis at CACHE_SIZE=2.
CONFIGS=(
  "base|| -OPS_DIAGS=2"
  "largetile|| -OPS_DIAGS=2 OPS_TILING_MAXDEPTH=6 OPS_TILESIZE_X=1000 OPS_TILESIZE_Y=1000 OPS_TILESIZE_Z=1000"
  "cache1|| -OPS_DIAGS=2 OPS_TILING_MAXDEPTH=6 OPS_CACHE_SIZE=1"
  "cache2|| -OPS_DIAGS=2 OPS_TILING_MAXDEPTH=6 OPS_CACHE_SIZE=2"
  "cache4|| -OPS_DIAGS=2 OPS_TILING_MAXDEPTH=6 OPS_CACHE_SIZE=4"
  "cache2_ispec|+ispec| -OPS_DIAGS=2 OPS_TILING_MAXDEPTH=6 OPS_CACHE_SIZE=2"
  "cache2_all|all| -OPS_DIAGS=2 OPS_TILING_MAXDEPTH=6 OPS_CACHE_SIZE=2"
)

echo "sweep start $(date -Iseconds) host=$(hostname) bin=$BIN" | tee results/sweep.log
echo "SLURM_JOB_ID=${SLURM_JOB_ID:-none} nodes=${SLURM_NNODES:-?} ntasks=${SLURM_NTASKS:-?}" | tee -a results/sweep.log

run_one() {
  local name="$1"
  local exec_spec="$2"
  local extra="$3"
  local log="results/${name}.out"
  echo "===== $(date -Iseconds) START $name SENGA_OPS_EXEC='${exec_spec}' ${extra} =====" | tee -a results/sweep.log
  if [[ -n "$exec_spec" ]]; then
    export SENGA_OPS_EXEC="$exec_spec"
  else
    unset SENGA_OPS_EXEC || true
  fi
  # shellcheck disable=SC2086
  time $SRUN "$BIN" $extra > "$log" 2>&1
  local rc=$?
  echo "===== $(date -Iseconds) END $name rc=$rc =====" | tee -a results/sweep.log
  if grep -q 'nxglbl=\s*256\b' "$log"; then
    echo "ERROR: binary compiled at 256^3 (nxglbl=256). Rebuild with SENGA_NX=512 and resubmit." | tee -a results/sweep.log
    return 99
  fi
  return $rc
}

fail=0
for spec in "${CONFIGS[@]}"; do
  IFS='|' read -r name exec extra <<< "$spec"
  run_one "$name" "$exec" "$extra"
  rc=$?
  if [[ $rc -ne 0 ]]; then
    fail=1
    if [[ $rc -eq 99 ]]; then
      echo "aborting remaining configs" | tee -a results/sweep.log
      break
    fi
    echo "RUN FAILED: $name (continuing)" | tee -a results/sweep.log
  fi
done

if [[ -f scripts/collate_senga_logs.py ]]; then
  python3 scripts/collate_senga_logs.py results | tee -a results/sweep.log
fi

cp -f scripts/collate_senga_logs.py results/ 2>/dev/null || true
tar -czf results/archer2_8node_512_results.tgz \
  -C results --exclude='archer2_8node_512_results.tgz' .
echo "sweep done $(date -Iseconds) tarball=results/archer2_8node_512_results.tgz fail=$fail" | tee -a results/sweep.log
exit $fail
