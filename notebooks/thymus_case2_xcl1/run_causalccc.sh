#!/bin/bash
# Render 04_causalccc.qmd: CausalCCC run CD8αα(I)/NK -> DC1/aDC1 with sequencing depth (n_counts) as a variable.
# The gene selection (selected_genes.rds) is reloaded, so only MIIC runs: about 30 minutes on 8 threads.
# The first run, without depth, is kept in results/thymus_case2_xcl1/04_causalccc/miic_input_without_depth/.
#
# Usage, from the root of the repository in WSL/Linux (keep the terminal open and the PC awake until the end):
#   bash notebooks/thymus_case2_xcl1/run_causalccc.sh
set -uo pipefail

cd "$(dirname "$0")"
log_dir="../../results/thymus_case2_xcl1/04_causalccc/logs"
mkdir -p "$log_dir"
log="$log_dir/with_depth.log"
echo "$(date '+%F %T')  start  case 2 with depth  (log: $log)"
start=$(date +%s)
if quarto render 04_causalccc.qmd > "$log" 2>&1; then
  echo "$(date '+%F %T')  done   in $(( ($(date +%s) - start) / 60 )) min -> 04_causalccc.html"
else
  echo "$(date '+%F %T')  FAILED after $(( ($(date +%s) - start) / 60 )) min, see $log"
  exit 1
fi
