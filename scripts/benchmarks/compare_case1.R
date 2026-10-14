# Case study 1: compare the local MIIC runs of run_case1_benchmarks.sh with the authors' server job
# (eaP9Br0llvzxN07zZZEC) and with the stock one-thread run, and tabulate the running times.
# Usage (from the root of the repository, in WSL/Linux, after run_case1_benchmarks.sh):
#   Rscript scripts/benchmarks/compare_case1.R
# Writes results/benchmarks/case1_server_vs_local.tsv and prints it as a Markdown table.

suppressPackageStartupMessages(library(tidyverse))
source("scripts/benchmarks/compare_outputs.R")

job_dir   <- "data/case1/causalccc_server_eaP9Br0llvzxN07zZZEC"
bench_dir <- "results/benchmarks/miic/case1"

# Wall-clock time (seconds) and peak memory (MB) from a GNU time -v report (as in notebook 10)
read_gnu_time <- function(file) {
  lines <- readLines(file)
  wall  <- str_match(lines, "Elapsed \\(wall clock\\) time \\(h:mm:ss or m:ss\\): (.*)$")[, 2] |> na.omit()
  parts <- as.numeric(str_split_1(wall, ":"))
  rss   <- str_match(lines, "Maximum resident set size \\(kbytes\\): ([0-9]+)")[, 2] |> na.omit()
  tibble(wall_s = sum(parts * 60^rev(seq_along(parts) - 1)), peak_mb = as.numeric(rss) / 1024)
}

# The server's network as written by the run itself (21 January 2025); the April 2025 edit of
# edgesList.miic.summary.txt only changes the values of the 18 L-R links (data/README.md)
server <- list(all.edges.summary = read_tsv(file.path(job_dir, "edgesList.miic.summary.txt_save"),
                                            show_col_types = FALSE))
# Server time: miic prints "TIME : <s>" at the end of its search (2 threads, as in script.sh)
server_log <- readLines(file.path(job_dir, "whole_run.log"))
server_s   <- as.numeric(str_match(server_log, "R : TIME : ([0-9.]+)")[, 2] |> na.omit())

runs <- tibble(run = list.dirs(bench_dir, full.names = FALSE, recursive = FALSE)) |>
  mutate(dir    = file.path(bench_dir, run),
         status = map_chr(file.path(dir, "status.txt"), \(f) if (file.exists(f)) readLines(f)[1] else "not run")) |>
  separate_wider_regex(run, c(variant = ".*", "_t", threads = "[0-9]+"), cols_remove = FALSE) |>
  mutate(threads = as.integer(threads)) |>
  arrange(variant, threads)

ok_runs <- filter(runs, status == "ok")
result  <- set_names(map(file.path(ok_runs$dir, "miic_result.rds"), readRDS), ok_runs$run)
reference <- if ("baseline_t1" %in% names(result)) result$baseline_t1 else NULL

local_rows <- ok_runs |>
  mutate(time = map(file.path(dir, "time.txt"), read_gnu_time),
         miic_s = map_dbl(run, \(r) result[[r]]$time[["total"]]),
         vs_server = map(run, \(r) compare_miic(server, result[[r]])),
         vs_stock_t1 = map(run, \(r) if (is.null(reference)) tibble(verdict = NA_character_)
                                    else compare_miic(reference, result[[r]]))) |>
  unnest_wider(time) |>
  mutate(links = map_int(vs_server, "links_new"),
         verdict_vs_server = map_chr(vs_server, "verdict"),
         orientations_changed_vs_server = map_int(vs_server, "orientations_changed"),
         causal_changed_vs_server = map_int(vs_server, "causal_changed"),
         verdict_vs_stock_t1 = map_chr(vs_stock_t1, "verdict"))

server_links <- miic_links(server)
table <- bind_rows(
  tibble(setup = "CausalCCC web server (Curie, 2025-01-21)", miic = "server", threads = 2L, miic_s = server_s,
         links = sum(server_links$type == "P"), verdict_vs_server = "(reference)"),
  local_rows |>
    transmute(setup = "local WSL (this laptop)",
              miic = if_else(variant == "baseline", "stock (causalccc branch)", "patched (deterministic)"),
              threads, miic_s, wall_s, peak_mb, links, verdict_vs_server, orientations_changed_vs_server,
              causal_changed_vs_server, verdict_vs_stock_t1)) |>
  mutate(speedup_vs_server = server_s / miic_s,
         speedup_vs_stock_t1 = if ("baseline_t1" %in% ok_runs$run)
           local_rows$miic_s[local_rows$run == "baseline_t1"] / miic_s else NA_real_)

failed <- filter(runs, status != "ok")
if (nrow(failed) > 0) message("Runs not ok: ", paste(failed$run, failed$status, sep = " = ", collapse = "; "))

write_tsv(table, "results/benchmarks/case1_server_vs_local.tsv")
fmt <- table |>
  mutate(across(c(miic_s, wall_s), \(x) if_else(is.na(x), "", sprintf("%.0f s (%.1f min)", x, x / 60))),
         across(c(speedup_vs_server, speedup_vs_stock_t1), \(x) if_else(is.na(x), "", sprintf("x%.2f", x))),
         peak_mb = if_else(is.na(peak_mb), "", sprintf("%.0f", peak_mb)))
print(knitr::kable(fmt, format = "pipe"))
