# Run MIIC with the settings of the CausalCCC web server on three input files
# (mosaic table, state order, L-R pairs), e.g. on a computing cluster.
#
# Usage, from the repository root:
#   Rscript scripts/run_miic.R <input.tsv> <state_order.tsv> <interact_edges.tsv> <output_dir> [n_threads] [settings.rds]
# Example:
#   Rscript scripts/run_miic.R results/case2/06_rebuilt_network/input.tsv \
#     results/case2/06_rebuilt_network/state_order.tsv results/case2/06_rebuilt_network/interact_edges.tsv \
#     results/case2/06_rebuilt_network
#
# With the installed miic (branch causalccc), keep n_threads = 1: with several threads it crashes (C stack errors,
# notebooks 06 and 10). With the patched "deterministic" variant (docs/INSTALL.md, section 4), 8 threads give the
# one-thread links, orientations and genuine causations about 3 times faster (notebook 10):
#   R_LIBS=~/R/miic_variants/deterministic Rscript scripts/run_miic.R <input> <state_order> <edges> <output_dir> 8
# settings.rds (optional): a named list of miic() arguments that replace the settings below, e.g. those of another
# server job (tutorial, section 8b); n_threads always comes from the command line. Without it, the run is unchanged.
# Writes <output_dir>/miic_result.rds, which the notebooks reload instead of running MIIC again; the settings used are
# stored in its attribute "settings".
# The settings are those of the server runs (notebooks 02 and 04), with the same seed as the notebooks.

suppressPackageStartupMessages({
  library(readr)
  library(miic)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript scripts/run_miic.R <input> <state_order> <interact_edges> <output_dir> [n_threads] [settings.rds]")
}

n_threads <- if (length(args) >= 5) as.integer(args[5]) else 1
dir.create(args[4], recursive = TRUE, showWarnings = FALSE)

settings <- list(
  cplx            = "nml",
  orientation     = TRUE,
  latent          = "orientation",
  ori_proba_ratio = 0.01,
  n_shuffles      = 100,
  conf_threshold  = 0.01,
  test_mar        = TRUE,
  consistent      = "no",
  n_eff           = -1,
  mode            = "S"
)
if (length(args) >= 6) {
  custom <- readRDS(args[6])
  unknown <- setdiff(names(custom), names(formals(miic)))
  if (length(unknown) > 0) stop("Not miic() arguments in ", args[6], ": ", paste(unknown, collapse = ", "))
  settings <- modifyList(settings, custom[names(custom) != "n_threads"])
}
settings$n_threads <- n_threads

mosaic         <- read_tsv(args[1], guess_max = Inf, show_col_types = FALSE)
state_order    <- read_tsv(args[2], show_col_types = FALSE)
interact_edges <- read_tsv(args[3], col_names = c("ligand", "receptor"), show_col_types = FALSE)

set.seed(2026)
start <- Sys.time()
miic_result <- do.call(miic, c(list(input_data     = as.data.frame(mosaic),
                                     state_order    = as.data.frame(state_order),
                                     interact_edges = as.data.frame(interact_edges)),
                               settings))
attr(miic_result, "run_time") <- Sys.time() - start
attr(miic_result, "settings") <- settings

saveRDS(miic_result, file.path(args[4], "miic_result.rds"))
message("Done in ", format(attr(miic_result, "run_time")), ": ", file.path(args[4], "miic_result.rds"))
