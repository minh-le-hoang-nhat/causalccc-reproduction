#!/bin/bash
# Render 04_causalccc.qmd for the fetal and the postnatal run, one after the other.
# Each run: gene selection (about 1 hour) + MIIC on 8 threads (about 20 minutes); about 2.5 hours for both.
# Results that already exist (selected_genes.rds, miic_result.rds) are reloaded, so an interrupted run can be restarted.
#
# Usage, from the root of the repository in WSL/Linux (keep the terminal open and the PC awake until the end):
#   bash notebooks/thymus_case3_xcl1_age/run_both.sh                # both runs
#   bash notebooks/thymus_case3_xcl1_age/run_both.sh postnatal      # one run only
#   bash notebooks/thymus_case3_xcl1_age/run_both.sh --dry-run      # check both runs in about a minute (no gene selection, no MIIC)
set -uo pipefail

cd "$(dirname "$0")"
dry_run=false
runs=()
for arg in "$@"; do
  case "$arg" in
    --dry-run)        dry_run=true ;;
    fetal|postnatal)  runs+=("$arg") ;;
    *) echo "Unknown argument: $arg (expected fetal, postnatal or --dry-run)"; exit 1 ;;
  esac
done
[ ${#runs[@]} -eq 0 ] && runs=(fetal postnatal)

log_dir="../../results/thymus_case3_xcl1_age/04_causalccc/logs"
mkdir -p "$log_dir"
status=0
for run in "${runs[@]}"; do
  suffix=$([ "$dry_run" = true ] && echo "_dryrun" || echo "")
  log="$log_dir/${run}${suffix}.log"
  echo "$(date '+%F %T')  start  $run$suffix  (log: $log)"
  start=$(date +%s)
  if quarto render 04_causalccc.qmd -P run:"$run" -P dry_run:"$dry_run" --output "04_causalccc_${run}${suffix}.html" > "$log" 2>&1; then
    echo "$(date '+%F %T')  done   $run$suffix  in $(( ($(date +%s) - start) / 60 )) min -> 04_causalccc_${run}${suffix}.html"
  else
    echo "$(date '+%F %T')  FAILED $run$suffix  after $(( ($(date +%s) - start) / 60 )) min, see $log"
    status=1
  fi
done
exit $status
