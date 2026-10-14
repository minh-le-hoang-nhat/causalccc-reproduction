#!/usr/bin/env bash
# =============================================================================
# One MIIC benchmark run: a miic variant (built with build_miic_variant.sh) on a humanSkin input (case study 2),
# with the settings of scripts/run_miic.R, timed with GNU time.
#
# Usage (from the root of the repository):
#   bash scripts/benchmarks/bench_miic.sh <variant> <dataset> <n_threads>
#   dataset: rebuilt (notebook 06 input, 217 variables), authors (authors' server input, 238 variables) or dedup
#   (authors' input without its 14 duplicated columns, notebook 07, 224 variables), or case1 (authors' server input of
#   case study 1, 144 variables)
# Writes results/benchmarks/miic/<dataset>/<variant>_t<n_threads>/:
#   miic_result.rds, run.log (MIIC messages), time.txt (GNU time: wall clock, peak memory), status.txt
# A run whose status.txt says "ok" is skipped when the script is called again.
# =============================================================================

set -uo pipefail

variant="$1"; dataset="$2"; n_threads="$3"
lib="$HOME/R/miic_variants/$variant"
out="results/benchmarks/miic/$dataset/${variant}_t${n_threads}"

case "$dataset" in
  rebuilt)
    dir="results/case2/06_rebuilt_network"
    files=("$dir/input.tsv" "$dir/state_order.tsv" "$dir/interact_edges.tsv") ;;
  authors)
    dir="data/case2/causalccc_server_w0788tM1wzXGHLzNaKKa"
    files=("$dir/input.tsv" "$dir/clean_state_order.tsv" "$dir/clean_interact_edges.tsv") ;;
  dedup)   # the authors' input without its 14 duplicated columns (notebook 07)
    dir="results/case2/07_differences/authors_dedup"
    files=("$dir/input.tsv" "$dir/state_order.tsv" "$dir/interact_edges.tsv") ;;
  case1)   # authors' server input of case study 1 (144 variables, no duplicated columns)
    dir="data/case1/causalccc_server_eaP9Br0llvzxN07zZZEC"
    files=("$dir/input.tsv" "$dir/clean_state_order.tsv" "$dir/clean_interact_edges.tsv") ;;
  *) echo "Unknown dataset: $dataset" >&2; exit 1 ;;
esac

if [[ ! -d "$lib" ]]; then echo "Variant not built: $lib" >&2; exit 1; fi
if [[ -f "$out/status.txt" && "$(cat "$out/status.txt")" == "ok" ]]; then
  echo "Already done: $out"
  exit 0
fi
mkdir -p "$out"

# Which miic is loaded, to be sure the variant is the one used
R_LIBS="$lib" Rscript -e 'cat(find.package("miic"), "\n")' > "$out/miic_path.txt" 2>&1
if ! grep -q "$lib" "$out/miic_path.txt"; then
  echo "The variant is not the miic loaded first: $(cat "$out/miic_path.txt")" >&2
  exit 1
fi

echo "$(date '+%F %T') start $out"
R_LIBS="$lib" /usr/bin/time -v -o "$out/time.txt" \
  Rscript scripts/run_miic.R "${files[@]}" "$out" "$n_threads" > "$out/run.log" 2>&1
status=$?
if [[ $status -eq 0 && -f "$out/miic_result.rds" ]]; then echo ok > "$out/status.txt"; else echo "failed ($status)" > "$out/status.txt"; fi
echo "$(date '+%F %T') end $out: $(cat "$out/status.txt")"
