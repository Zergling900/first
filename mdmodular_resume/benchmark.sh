#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")" && pwd)"
THREADS="${THREADS:-8}"
WORK_DIR="${WORK_DIR:-/tmp/mdbench_$(date +%Y%m%d_%H%M%S)_$$}"
BIN_PATH="${BIN_PATH:-$WORK_DIR/run.out}"
CXX="${CXX:-g++}"
CXXFLAGS="${CXXFLAGS:--std=c++17 -O2 -fopenmp}"

if [[ "${1:-}" == "--help" ]]; then
    cat <<'USAGE'
Usage:
  ./benchmark.sh
Environment variables:
  THREADS   OpenMP threads (default: 8)
  WORK_DIR  temporary benchmark workspace
  BIN_PATH  output binary path (default: $WORK_DIR/run.out)
  CXX       compiler (default: g++)
  CXXFLAGS  compiler flags (default: -std=c++17 -O2 -fopenmp)
USAGE
    exit 0
fi

mkdir -p "$WORK_DIR"
rm -rf "$WORK_DIR/mdmodular_resume"
cp -a "$ROOT_DIR" "$WORK_DIR/mdmodular_resume"

# Build
$CXX $CXXFLAGS "$ROOT_DIR"/*.cpp -o "$BIN_PATH"
cp "$BIN_PATH" "$WORK_DIR/mdmodular_resume/run.out"

# Short fixed benchmark profile (same profile used in optimization checks).
perl -pi -e 's/^endtime\s*=.*/endtime = 350.0 #fs benchmark/; s/^dT\s*=.*/dT = 200 #K benchmark/; s/^T_end\s*=.*/T_end = 350 #K benchmark/; s/^T_time\s*=.*/T_time = 100 #fs benchmark/; s/^output_Data_time\s*=.*/output_Data_time = 1000.0 #fs benchmark/; s/^output_Et_time\s*=.*/output_Et_time = 1000 #fs benchmark/' "$WORK_DIR/mdmodular_resume/parameter.p1"
perl -pi -e 's/^extra_steps_max\s*=.*/extra_steps_max = 0/; s/^plateau_blocks\s*=.*/plateau_blocks = 0/' "$WORK_DIR/mdmodular_resume/parameter.p3"

pushd "$WORK_DIR/mdmodular_resume" >/dev/null
OMP_NUM_THREADS="$THREADS" /usr/bin/time -f 'WALL=%e' ./run.out > bench.log 2>&1
popd >/dev/null

STEP_S="$(rg -o 'step/s = [0-9.]+ steps/s' "$WORK_DIR/mdmodular_resume/bench.log" | awk '{print $3}' | tail -n 1)"
TIME_ALL_MS="$(rg -o 'time_all = [0-9.]+ ms' "$WORK_DIR/mdmodular_resume/bench.log" | awk '{print $3}' | tail -n 1)"
WALL_S="$(rg -o 'WALL=[0-9.]+' "$WORK_DIR/mdmodular_resume/bench.log" | awk -F= '{print $2}' | tail -n 1)"

echo "benchmark_workdir=$WORK_DIR/mdmodular_resume"
echo "threads=$THREADS"
echo "time_all_ms=${TIME_ALL_MS:-NA}"
echo "step_per_s=${STEP_S:-NA}"
echo "wall_s=${WALL_S:-NA}"
