#!/bin/bash
# Check that default SENGA (ispec off) decouples halo from cache blocking:
#   * one halo for the fused 411-loop species chain (comms of old -ispec)
#   * 10 cache-blocked segments at the old +ispec points (tiles of old +ispec)
#   * last-step dump matches SENGA_OPS_EXEC=+ispec
#
# Needs a tiled binary built against OPS with ops_tiling_break, e.g.
#   make senga2_mpi_tiled SENGA_NX=256
#
# Usage (from the SENGA2 root):
#   scripts/test_tiling_break.sh
#
# Optional env: NP (default 32), BIN, H5DIFF, STEPS (default 3)

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

BIN="${BIN:-}"
if [[ -z "$BIN" ]]; then
  if [[ -x "$ROOT/senga2_mpi_tiled" ]]; then
    BIN="$ROOT/senga2_mpi_tiled"
  elif [[ -x "$ROOT/src/senga2_mpi_tiled" ]]; then
    BIN="$ROOT/src/senga2_mpi_tiled"
  else
    echo "ERROR: senga2_mpi_tiled not found. Build with: make senga2_mpi_tiled SENGA_NX=256"
    exit 1
  fi
fi

STEPS="${STEPS:-3}"
NP="${NP:-32}"
CONT_IN="input/cont_256_${STEPS}.dat"
if [[ ! -f "$CONT_IN" ]]; then
  echo "ERROR: $CONT_IN missing"
  exit 1
fi

H5DIFF="${H5DIFF:-}"
if [[ -z "$H5DIFF" ]]; then
  if command -v h5diff >/dev/null 2>&1; then
    H5DIFF="$(command -v h5diff)"
  elif [[ -x /home/shared/software/op-dsl/hdf5/bin/h5diff ]]; then
    H5DIFF=/home/shared/software/op-dsl/hdf5/bin/h5diff
  fi
fi

LOGDIR="$ROOT/tiling_split_logs/tiling_break_test"
mkdir -p "$LOGDIR"
DUMP="output/timestep$(printf '%08d' "$STEPS").h5"

COMMON="-OPS_DIAGS=3 OPS_TILING_MAXDEPTH=6 OPS_CACHE_SIZE=2"
export OMP_NUM_THREADS=1
export SENGA_DUMP_LAST=1

cp -f input/cont.dat "$LOGDIR/cont.dat.bak"
cp -f "$CONT_IN" input/cont.dat
restore_cont() {
  if [[ -f "$LOGDIR/cont.dat.bak" ]]; then
    mv -f "$LOGDIR/cont.dat.bak" input/cont.dat
  fi
}
trap restore_cont EXIT

run_one() {
  local name="$1"
  local exec_spec="$2"
  echo "===== $(date -Iseconds) START $name SENGA_OPS_EXEC='$exec_spec' ====="
  rm -f "$DUMP"
  if [[ -n "$exec_spec" ]]; then
    export SENGA_OPS_EXEC="$exec_spec"
  else
    unset SENGA_OPS_EXEC || true
  fi
  mpirun -np "$NP" --use-hwthread-cpus -bind-to core "$BIN" $COMMON \
    > "$LOGDIR/${name}.out" 2>&1
  if [[ -f "$DUMP" ]]; then
    cp -f "$DUMP" "$LOGDIR/${name}.h5"
  else
    echo "ERROR: no dump $DUMP after $name"
    return 1
  fi
  echo "===== $(date -Iseconds) END $name ====="
}

fail() { echo "FAIL: $*" >&2; exit 1; }
pass() { echo "PASS: $*"; }

run_one default ""
run_one ispec "+ispec"

DEF="$LOGDIR/default.out"
ISP="$LOGDIR/ispec.out"

grep -q 'total_ttime' "$DEF" || fail "default run did not finish"
grep -q 'total_ttime' "$ISP" || fail "+ispec run did not finish"
grep -q 'ispec= F' "$DEF" || fail "default should have ispec off"
grep -q 'ispec= T' "$ISP" || fail "+ispec should have ispec on"

if ! grep -q 'Tiling breaks: one halo for 411 loops, 10 cache-blocked segments' "$DEF"; then
  fail "default missing 'Tiling breaks: one halo for 411 loops, 10 cache-blocked segments'"
fi
if grep -q 'Tiling breaks:' "$ISP"; then
  fail "+ispec should not print Tiling breaks (it uses ops_execute per species)"
fi

# Species segments are 40-loop cache-blocked plans (not the fused 411-loop skew).
if ! grep -q 'unblocked loops: 0/40' "$DEF"; then
  fail "default should execute 40-loop cache-blocked species segments"
fi
if ! grep -q 'unblocked loops: 173/411' "$DEF"; then
  echo "NOTE: fused 411-loop plan is still built for the halo list (expected)"
fi

def_msgs="$(sed -n 's/.*Average number of MPI messages per process: *//p' "$DEF" | tail -1)"
isp_msgs="$(sed -n 's/.*Average number of MPI messages per process: *//p' "$ISP" | tail -1)"
def_mb="$(sed -n 's/.*Average volume of MPI communications per process: *//p' "$DEF" | awk '{print $1}' | tail -1)"
isp_mb="$(sed -n 's/.*Average volume of MPI communications per process: *//p' "$ISP" | awk '{print $1}' | tail -1)"

python3 - "$def_msgs" "$isp_msgs" "$def_mb" "$isp_mb" <<'PY'
import sys
def_msgs, isp_msgs, def_mb, isp_mb = map(float, sys.argv[1:])
if def_msgs <= 0 or isp_msgs <= 0:
    raise SystemExit("missing MPI message counts")
ratio_m = isp_msgs / def_msgs
ratio_v = isp_mb / def_mb if def_mb > 0 else 0.0
print(f"msgs default={def_msgs:.0f} +ispec={isp_msgs:.0f} ratio={ratio_m:.3f}")
print(f"MB    default={def_mb:.2f} +ispec={isp_mb:.2f} ratio={ratio_v:.3f}")
# Coupled +ispec re-exchanges each species: about 2x messages and volume.
if not (1.8 <= ratio_m <= 2.2):
    raise SystemExit(f"expected +ispec/~2x messages, got {ratio_m:.3f}")
if not (1.8 <= ratio_v <= 2.2):
    raise SystemExit(f"expected +ispec/~2x volume, got {ratio_v:.3f}")
PY

if [[ -z "$H5DIFF" ]]; then
  fail "h5diff not found (set H5DIFF=...)"
fi
if ! "$H5DIFF" -p 1e-10 "$LOGDIR/default.h5" "$LOGDIR/ispec.h5"; then
  fail "h5diff -p 1e-10 default vs +ispec"
fi
pass "h5diff -p 1e-10 default (tiling_break) matches +ispec"

echo "TEST PASSED"
