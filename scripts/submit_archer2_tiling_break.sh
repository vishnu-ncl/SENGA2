#!/bin/bash
#SBATCH --job-name=senga2_tb_512
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=128
#SBATCH --time=03:00:00
#SBATCH --partition=standard
#SBATCH --qos=standard
#SBATCH --account=e647
#SBATCH --mem=0
#SBATCH --output=results/slurm-%j.out
#SBATCH --error=results/slurm-%j.err

# 512^3 / 10-step / 8-node (1024 rank) tiling_break performance + correctness.
#
# Default (ispec off) uses ops_tiling_break: one halo for the fused species
# chain, cache-blocked segments at the old +ispec points.
#
# BEFORE submitting, rebuild OPS (feature/tiling_comm_separate) and the app:
#   source /work/e647/e647/ireguly5/source_archer2_gcc
#   cd .../ops/fortran && make f_mpi
#   cd .../SENGA2_OPS/src && make senga2_mpi_tiled SENGA_NX=512 -j8
#
# Then from the app root:
#   mkdir -p results
#   sbatch scripts/submit_archer2_tiling_break.sh
#
# After the job:
#   scp archer2:.../results/archer2_tiling_break_results.tgz .
#   python3 scripts/collate_senga_logs.py results
#   cat results/correctness.txt
#
# 512^3 last-step dumps are tens of GB. They stay under results/dumps/ and
# are excluded from the tarball. Override CONT / STEPS if needed.

set -u

source /work/e647/e647/ireguly5/source_archer2_gcc

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

CONT_IN="${CONT_IN:-input/cont_512_10.dat}"
STEPS="${STEPS:-10}"
if [[ ! -f "$CONT_IN" ]]; then
  echo "ERROR: $CONT_IN missing"
  exit 1
fi

H5DIFF="${H5DIFF:-}"
if [[ -z "$H5DIFF" ]]; then
  if command -v h5diff >/dev/null 2>&1; then
    H5DIFF="$(command -v h5diff)"
  elif [[ -n "${HDF5_DIR:-}" && -x "${HDF5_DIR}/bin/h5diff" ]]; then
    H5DIFF="${HDF5_DIR}/bin/h5diff"
  fi
fi

mkdir -p results results/dumps
export OMP_NUM_THREADS=1

cp -f input/cont.dat "results/cont.dat.bak.$$"
cp -f "$CONT_IN" input/cont.dat

restore_cont() {
  if [[ -f "results/cont.dat.bak.$$" ]]; then
    mv -f "results/cont.dat.bak.$$" input/cont.dat
  fi
}
trap restore_cont EXIT

SRUN="srun --hint=nomultithread --distribution=block:block --unbuffered"
DUMP="output/timestep$(printf '%08d' "$STEPS").h5"
CHECK=results/correctness.txt
: > "$CHECK"

echo "tiling_break job start $(date -Iseconds) host=$(hostname) bin=$BIN cont=$CONT_IN steps=$STEPS" | tee results/sweep.log
echo "SLURM_JOB_ID=${SLURM_JOB_ID:-none} nodes=${SLURM_NNODES:-?} ntasks=${SLURM_NTASKS:-?}" | tee -a results/sweep.log
echo "h5diff=${H5DIFF:-MISSING}" | tee -a results/sweep.log

run_one() {
  local name="$1"
  local exec_spec="$2"
  local extra="$3"
  local dump_mode="$4"
  local log="results/${name}.out"
  echo "===== $(date -Iseconds) START $name SENGA_OPS_EXEC='${exec_spec}' ${extra} dump=${dump_mode} =====" | tee -a results/sweep.log
  if [[ -n "$exec_spec" ]]; then
    export SENGA_OPS_EXEC="$exec_spec"
  else
    unset SENGA_OPS_EXEC || true
  fi
  if [[ "$dump_mode" == "dump" ]]; then
    export SENGA_DUMP_LAST=1
    rm -f "$DUMP"
  else
    unset SENGA_DUMP_LAST || true
  fi
  # shellcheck disable=SC2086
  time $SRUN "$BIN" $extra > "$log" 2>&1
  local rc=$?
  echo "===== $(date -Iseconds) END $name rc=$rc =====" | tee -a results/sweep.log
  if grep -q 'nxglbl=\s*256\b' "$log"; then
    echo "ERROR: binary compiled at 256^3 (nxglbl=256). Rebuild with SENGA_NX=512." | tee -a results/sweep.log "$CHECK"
    return 99
  fi
  if [[ "$dump_mode" == "dump" ]]; then
    if [[ -f "$DUMP" ]]; then
      cp -f "$DUMP" "results/dumps/${name}.h5"
    else
      echo "ERROR: no dump $DUMP after $name" | tee -a results/sweep.log "$CHECK"
      return 1
    fi
  fi
  return $rc
}

note() { echo "$*" | tee -a results/sweep.log "$CHECK"; }

# Performance: no HDF5. Default exec is tiling_break (ispec off).
PERF=(
  "base|| -OPS_DIAGS=2"
  "largetile|| -OPS_DIAGS=2 OPS_TILING_MAXDEPTH=6 OPS_TILESIZE_X=1000 OPS_TILESIZE_Y=1000 OPS_TILESIZE_Z=1000"
  "cache1|| -OPS_DIAGS=2 OPS_TILING_MAXDEPTH=6 OPS_CACHE_SIZE=1"
  "cache2|| -OPS_DIAGS=2 OPS_TILING_MAXDEPTH=6 OPS_CACHE_SIZE=2"
  "cache4|| -OPS_DIAGS=2 OPS_TILING_MAXDEPTH=6 OPS_CACHE_SIZE=4"
  "cache2_ispec|+ispec| -OPS_DIAGS=2 OPS_TILING_MAXDEPTH=6 OPS_CACHE_SIZE=2"
  "cache2_all|all| -OPS_DIAGS=2 OPS_TILING_MAXDEPTH=6 OPS_CACHE_SIZE=2"
)

# Correctness: last-step dump, DIAGS=3 so we see Tiling breaks / sanity lines.
CORRECT=(
  "correct_break|| -OPS_DIAGS=3 OPS_TILING_MAXDEPTH=6 OPS_CACHE_SIZE=2"
  "correct_ispec|+ispec| -OPS_DIAGS=3 OPS_TILING_MAXDEPTH=6 OPS_CACHE_SIZE=2"
)

fail=0
abort=0
for spec in "${PERF[@]}"; do
  IFS='|' read -r name exec extra <<< "$spec"
  run_one "$name" "$exec" "$extra" "nodump"
  rc=$?
  if [[ $rc -ne 0 ]]; then
    fail=1
    if [[ $rc -eq 99 ]]; then
      abort=1
      echo "aborting remaining configs" | tee -a results/sweep.log
      break
    fi
    echo "RUN FAILED: $name (continuing)" | tee -a results/sweep.log
  fi
done

if [[ $abort -eq 0 ]]; then
  for spec in "${CORRECT[@]}"; do
    IFS='|' read -r name exec extra <<< "$spec"
    run_one "$name" "$exec" "$extra" "dump"
    rc=$?
    if [[ $rc -ne 0 ]]; then
      fail=1
      if [[ $rc -eq 99 ]]; then
        abort=1
        echo "aborting remaining configs" | tee -a results/sweep.log
        break
      fi
      echo "RUN FAILED: $name (continuing)" | tee -a results/sweep.log
    fi
  done
fi

note ""
note "## Correctness checks"
note ""

sanity=$(grep -h "terminal read dependency sanity check fail" results/*.out 2>/dev/null | wc -l | tr -d ' ')
if [[ "$sanity" -eq 0 ]]; then
  note "PASS: no terminal-read sanity check fails"
else
  note "FAIL: $sanity terminal-read sanity check lines (rebuild OPS with the edge-tile fix)"
  fail=1
fi

if [[ -f results/correct_break.out ]]; then
  if grep -q 'Tiling breaks: one halo for 411 loops, 10 cache-blocked segments' results/correct_break.out; then
    note "PASS: default printed one halo / 10 cache-blocked segments"
  else
    note "FAIL: default missing 'Tiling breaks: one halo for 411 loops, 10 cache-blocked segments'"
    fail=1
  fi
  if grep -q 'unblocked loops: 0/40' results/correct_break.out; then
    note "PASS: default executed 40-loop cache-blocked species segments"
  else
    note "FAIL: default missing 40-loop cache-blocked species segments"
    fail=1
  fi
else
  note "FAIL: results/correct_break.out missing"
  fail=1
fi

if [[ -f results/correct_ispec.out ]]; then
  if grep -q 'Tiling breaks:' results/correct_ispec.out; then
    note "FAIL: +ispec should not print Tiling breaks"
    fail=1
  else
    note "PASS: +ispec has no Tiling breaks (per-species ops_execute)"
  fi
fi

def_msgs="$(sed -n 's/.*Average number of MPI messages per process: *//p' results/correct_break.out 2>/dev/null | tail -1)"
isp_msgs="$(sed -n 's/.*Average number of MPI messages per process: *//p' results/correct_ispec.out 2>/dev/null | tail -1)"
def_mb="$(sed -n 's/.*Average volume of MPI communications per process: *//p' results/correct_break.out 2>/dev/null | awk '{print $1}' | tail -1)"
isp_mb="$(sed -n 's/.*Average volume of MPI communications per process: *//p' results/correct_ispec.out 2>/dev/null | awk '{print $1}' | tail -1)"

if [[ -n "${def_msgs:-}" && -n "${isp_msgs:-}" ]]; then
  if python3 - "$def_msgs" "$isp_msgs" "${def_mb:-0}" "${isp_mb:-0}" >> "$CHECK" 2>> "$CHECK" <<'PY'
import sys
def_msgs, isp_msgs, def_mb, isp_mb = map(float, sys.argv[1:])
if def_msgs <= 0 or isp_msgs <= 0:
    raise SystemExit("FAIL: missing MPI message counts")
ratio_m = isp_msgs / def_msgs
ratio_v = isp_mb / def_mb if def_mb > 0 else 0.0
print(f"msgs default={def_msgs:.0f} +ispec={isp_msgs:.0f} ratio={ratio_m:.3f}")
print(f"MB    default={def_mb:.2f} +ispec={isp_mb:.2f} ratio={ratio_v:.3f}")
if not (1.8 <= ratio_m <= 2.2):
    raise SystemExit(f"FAIL: expected +ispec/~2x messages, got {ratio_m:.3f}")
if def_mb > 0 and not (1.8 <= ratio_v <= 2.2):
    raise SystemExit(f"FAIL: expected +ispec/~2x volume, got {ratio_v:.3f}")
print("PASS: +ispec is ~2x MPI messages and volume vs tiling_break")
PY
  then
    tail -5 "$CHECK" | tee -a results/sweep.log >/dev/null
  else
    note "FAIL: MPI volume comparison"
    fail=1
  fi
else
  note "FAIL: could not parse MPI message counts from correctness logs"
  fail=1
fi

if [[ -z "${H5DIFF:-}" ]]; then
  note "FAIL: h5diff not found (set H5DIFF=...)"
  fail=1
elif [[ -f results/dumps/correct_break.h5 && -f results/dumps/correct_ispec.h5 ]]; then
  if "$H5DIFF" -p 1e-10 results/dumps/correct_break.h5 results/dumps/correct_ispec.h5 \
      > results/h5diff.out 2>&1; then
    note "PASS: h5diff -p 1e-10 correct_break vs correct_ispec"
  else
    note "FAIL: h5diff -p 1e-10 correct_break vs correct_ispec (see results/h5diff.out)"
    fail=1
  fi
else
  note "FAIL: correctness dumps missing under results/dumps/"
  fail=1
fi

if [[ $fail -eq 0 ]]; then
  note "CORRECTNESS PASSED"
else
  note "CORRECTNESS FAILED"
fi

if [[ -f scripts/collate_senga_logs.py ]]; then
  python3 scripts/collate_senga_logs.py results | tee -a results/sweep.log
fi

cp -f scripts/collate_senga_logs.py results/ 2>/dev/null || true
tar -czf results/archer2_tiling_break_results.tgz \
  -C results --exclude='archer2_tiling_break_results.tgz' --exclude='dumps' .
echo "job done $(date -Iseconds) tarball=results/archer2_tiling_break_results.tgz fail=$fail" | tee -a results/sweep.log
if [[ $fail -ne 0 ]]; then
  exit 1
fi
exit 0
