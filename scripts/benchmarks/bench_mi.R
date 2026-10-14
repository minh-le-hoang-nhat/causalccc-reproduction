# One MI benchmark run (change 3 of notebook 10): the MI table of one cell type of humanSkin (case study 2), computed
# either by causalCCC.MIselection() ("original", with the Seurat 5 fix of notebook 05) or by compute_mi_table()
# ("fast", scripts/helpers.R). Same cells, genes and features of interest as the MI tables of notebook 05.
#
# Usage, from the repository root (wrapped in /usr/bin/time by bench_mi.sh):
#   Rscript scripts/benchmarks/bench_mi.R <original|fast> <cell type|all> <block_size> <output.rds>
# The output keeps the table (variables, foi, value, cell_type) and, as attribute "seconds", the time of the MI
# computation alone (without loading the data).

suppressPackageStartupMessages({
  library(tidyverse)
  library(Matrix)
  library(Seurat)
  library(miic)
  library(causalCCC)
})
source("scripts/helpers.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) stop("Usage: Rscript scripts/benchmarks/bench_mi.R <original|fast> <cell type|all> <block_size> <output.rds>")
method     <- args[1]
cell_types <- if (args[2] == "all") c("Fibroblasts", "T_cell", "DCs") else args[2]
block_size <- as.integer(args[3])

rebuilt  <- readRDS("results/case2/05_rebuild_input/rebuild_inputs.rds")
so       <- CreateSeuratObject(counts = rebuilt$counts, meta.data = rebuilt$meta)
metadata <- c("labels", "condition")
# Genes of interest of each cell type: the features of the MI tables of notebook 05, without the metadata
reference <- readRDS("results/case2/05_rebuild_input/mi_tables.rds")
goi_of    <- \(t) setdiff(unique(as.character(reference$foi[reference$cell_type == t])), metadata)

MIselection <- causalCCC.MIselection   # Seurat 5 fix of notebook 05
body(MIselection) <- parse(text = gsub('slot = "counts"', 'layer = "counts"', deparse(body(MIselection)), fixed = TRUE))[[1]]

counts <- GetAssayData(so, assay = "RNA", layer = "counts")
start  <- Sys.time()
tables <- map(cell_types, \(t) {
  mi <- if (method == "original") {
    invisible(capture.output(suppressMessages(
      out <- MIselection(so, interact_ident = "cell_type", oneinteract = t, goi = goi_of(t),
                         metadata_list = metadata, return_full = TRUE))))
    out
  } else {
    compute_mi_table(counts, so@meta.data, so$cell_type == t, goi = goi_of(t), metadata = metadata,
                     block_size = block_size)
  }
  as_tibble(mi) |> mutate(foi = as.character(foi), cell_type = t)
}) |> list_rbind()
attr(tables, "seconds") <- as.numeric(difftime(Sys.time(), start, units = "secs"))

dir.create(dirname(args[4]), recursive = TRUE, showWarnings = FALSE)
saveRDS(tables, args[4])
message(method, " ", args[2], " block ", block_size, ": ", round(attr(tables, "seconds"), 1), " s")
