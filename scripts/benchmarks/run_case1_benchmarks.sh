#!/usr/bin/env bash
# =============================================================================
# MIIC benchmarks on case study 1 (authors' server job eaP9Br0llvzxN07zZZEC), one after another (never in parallel,
# so timings are comparable): stock miic against the patched "deterministic" variant, to check that the variant gives
# the server's network and to compare running times with the server's (935 s of MIIC on 2 threads, whole_run.log).
# Runs already done are skipped, so the script can be stopped and started again.
# Usage (from the root of the repository, in WSL/Linux, after scripts/case1/download_case1.sh and
# build_miic_variant.sh):
#   bash scripts/benchmarks/run_case1_benchmarks.sh 2>&1 | tee -a results/benchmarks/queue.log
# =============================================================================

set -uo pipefail

runs=(
  # variant       dataset  threads
  "baseline       case1    1"   # reference: stock miic, one thread
  "deterministic  case1    1"   # patches, same thread count
  "deterministic  case1    2"   # the server's thread count
  "deterministic  case1    8"
  "baseline       case1    2"   # the server's set-up with stock miic (may crash with >1 thread)
  "deterministic  case1    4"   # scaling with threads (12 = all logical CPUs of the i7-1365U: 2 P-cores x 2 + 8 E-cores)
  "deterministic  case1    12"
)

for run in "${runs[@]}"; do
  read -r variant dataset threads <<< "$run"
  bash scripts/benchmarks/bench_miic.sh "$variant" "$dataset" "$threads"
done
echo "$(date '+%F %T') all case 1 runs done"
