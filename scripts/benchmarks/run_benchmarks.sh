#!/usr/bin/env bash
# =============================================================================
# The MIIC benchmarks of notebook 10, one after another (never in parallel, so timings are comparable).
# Runs already done are skipped, so the script can be stopped and started again.
#
# Variants needed (build_miic_variant.sh): baseline, no_rcout, early_conf, combined, baseline_idle and
# combined_again (copies of baseline and combined, so that the runs get their own folders), deterministic.
# Usage (from the root of the repository, in WSL/Linux, after building the variants):
#   bash scripts/benchmarks/run_benchmarks.sh 2>&1 | tee -a results/benchmarks/queue.log
# =============================================================================

set -uo pipefail

runs=(
  # variant     dataset  threads
  "baseline     rebuilt  1"   # reference
  "baseline     rebuilt  4"   # does the stock version crash with 4 threads?
  "no_rcout     rebuilt  1"   # change 1, outputs must be identical
  "no_rcout     rebuilt  4"   # change 1, multithreading
  "no_rcout     rebuilt  8"
  "early_conf   rebuilt  1"   # change 2 alone
  "combined     rebuilt  8"   # changes 1 + 2
  "baseline     authors  1"   # final check: reference on the authors' input
  "combined     authors  8"   # final check: must give the authors' network
  # Added after the first runs (notebook 10):
  "baseline_idle    rebuilt  1"   # reference again: the first one shared the machine with a notebook render
  "no_rcout         authors  8"   # orientations changed on the authors' input: is change 2 involved?
  "combined         dedup    8"   # ... are the duplicated columns the cause?
  "combined_again   authors  8"   # ... is the 8-thread run reproducible?
  "deterministic    authors  8"   # change 3: exact ties broken by position
  "deterministic    rebuilt  8"
)

for run in "${runs[@]}"; do
  read -r variant dataset threads <<< "$run"
  bash scripts/benchmarks/bench_miic.sh "$variant" "$dataset" "$threads"
done
echo "$(date '+%F %T') all runs done"
