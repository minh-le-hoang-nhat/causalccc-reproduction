# Run MIIC with the settings of the CausalCCC web server on three input files
# (mosaic table, state order, L-R pairs), e.g. on a computing cluster.
#
# Usage, from the repository root:
#   Rscript scripts/run_miic.R <input.tsv> <state_order.tsv> <interact_edges.tsv> <output_dir> [n_threads]
# Example:
#   Rscript scripts/run_miic.R results/case2/06_rebuilt_network/input.tsv \
#     results/case2/06_rebuilt_network/state_order.tsv results/case2/06_rebuilt_network/interact_edges.tsv \
#     results/case2/06_rebuilt_network
#
# Keep n_threads = 1: a run of the case study 2 input with 4 threads crashed (C stack errors, notebook 06).
# Writes <output_dir>/miic_result.rds, which the notebooks reload instead of running MIIC again.
# The settings are those of the server runs (notebooks 02 and 04), with the same seed as the notebooks.

suppressPackageStartupMessages({
  library(readr)
  library(miic)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop("Usage: Rscript scripts/run_miic.R <input> <state_order> <interact_edges> <output_dir> [n_threads]")

n_threads <- if (length(args) >= 5) as.integer(args[5]) else 1
dir.create(args[4], recursive = TRUE, showWarnings = FALSE)

mosaic         <- read_tsv(args[1], guess_max = Inf, show_col_types = FALSE)
state_order    <- read_tsv(args[2], show_col_types = FALSE)
interact_edges <- read_tsv(args[3], col_names = c("ligand", "receptor"), show_col_types = FALSE)

set.seed(2026)
start <- Sys.time()
miic_result <- miic(
  input_data      = as.data.frame(mosaic),
  state_order     = as.data.frame(state_order),
  interact_edges  = as.data.frame(interact_edges),
  cplx            = "nml",
  orientation     = TRUE,
  latent          = "orientation",
  ori_proba_ratio = 0.01,
  n_shuffles      = 100,
  conf_threshold  = 0.01,
  test_mar        = TRUE,
  consistent      = "no",
  n_eff           = -1,
  mode            = "S",
  n_threads       = n_threads
)
attr(miic_result, "run_time") <- Sys.time() - start

saveRDS(miic_result, file.path(args[4], "miic_result.rds"))
message("Done in ", format(attr(miic_result, "run_time")), ": ", file.path(args[4], "miic_result.rds"))
