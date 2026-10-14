#!/usr/bin/env bash
# =============================================================================
# The MI benchmarks of notebook 10 (change 3), one after another, each in its own R process timed with GNU time.
# Run when no other benchmark is running. Runs already done are skipped.
#
# Usage (from the root of the repository, in WSL/Linux):
#   bash scripts/benchmarks/run_mi_benchmarks.sh 2>&1 | tee -a results/benchmarks/queue.log
# Writes results/benchmarks/mi/<method>_<cell type>_b<block size>_r<repeat>.{rds,time.txt}
# =============================================================================

set -uo pipefail
out=results/benchmarks/mi
mkdir -p "$out"

runs=(
  # method   cell type  block  repeats
  "original  DCs        100    3"   # reference
  "fast      DCs        100    3"   # same blocks: effect of sparse matrix + undetected genes skipped
  "fast      DCs        500    1"   # block size sweep: larger blocks turned out slower
  "fast      DCs        1000   1"
  "fast      DCs        50     1"   # so smaller blocks are tried too
  "fast      DCs        200    1"
  "original  all        100    1"   # the 3 cell types, as in notebook 05
  "fast      all        100    1"
)

for run in "${runs[@]}"; do
  read -r method cell_type block repeats <<< "$run"
  for r in $(seq 1 "$repeats"); do
    name="$out/${method}_${cell_type}_b${block}_r${r}"
    if [[ -f "$name.rds" ]]; then echo "Already done: $name"; continue; fi
    echo "$(date '+%F %T') start $name"
    /usr/bin/time -v -o "$name.time.txt" Rscript scripts/benchmarks/bench_mi.R "$method" "$cell_type" "$block" "$name.rds" \
      > "$name.log" 2>&1 || echo "FAILED: $name (see $name.log)"
  done
done
echo "$(date '+%F %T') MI benchmarks done"
