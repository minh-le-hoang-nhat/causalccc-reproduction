# Case study 2 – input for the "duplicates" experiment (notebook 09).
#
# Rebuilds the mosaic table with the selection settings found in notebook 07 (one run per sender -> receiver direction,
# 25 genes per run, no cap of 8 on the ligands/receptors, both metadata as features of interest), then ADDS the 14
# duplicated columns of the authors' table (a copy of the original column, placed right after it).
# The cells are in the authors' row order; the columns are in OUR order (the authors' column order is unknown).
#
# Usage, from the repository root:
#   Rscript scripts/case2/build_duplicates_input.R
# Writes results/case2/09_duplicates_experiment/{input.tsv, state_order.tsv, interact_edges.tsv}.
# Then run MIIC on them with scripts/run_miic.R (about 3 hours, one thread).

suppressPackageStartupMessages({
  library(tidyverse)
  library(Matrix)   # methods for the sparse count matrix
})

root <- getwd()
if (!file.exists(file.path(root, "scripts", "helpers.R"))) stop("Run this script from the root of the repository.")
source(file.path(root, "scripts", "helpers.R"))

job_dir <- file.path(root, "data", "case2", "causalccc_server_w0788tM1wzXGHLzNaKKa")
nb05    <- file.path(root, "results", "case2", "05_rebuild_input")
out_dir <- file.path(root, "results", "case2", "09_duplicates_experiment")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cell_types <- c("Fibroblasts", "T_cell", "DCs")
strip      <- function(v) str_remove_all(v, "_(Fibroblasts|T_cell|DCs)")
metadata   <- c("labels", "condition")

rebuilt     <- readRDS(file.path(nb05, "rebuild_inputs.rds"))
mi_tables   <- readRDS(file.path(nb05, "mi_tables.rds"))
counts      <- rebuilt$counts
meta        <- rebuilt$meta
match_cells <- arrange(rebuilt$match_cells, row)   # the authors' row order
lr_roles    <- rebuilt$lr_roles

mosaic_authors <- read_tsv(file.path(job_dir, "input.tsv"), guess_max = Inf, show_col_types = FALSE)
state_authors  <- read_tsv(file.path(job_dir, "clean_state_order.tsv"), show_col_types = FALSE)
lr_pairs       <- read_tsv(file.path(job_dir, "clean_interact_edges.tsv"),
                           col_names = c("ligand", "receptor"), show_col_types = FALSE)

row_type_a <- case_when(!is.na(mosaic_authors$labels_Fibroblasts) ~ "Fibroblasts",
                        !is.na(mosaic_authors$labels_T_cell)      ~ "T_cell",
                        !is.na(mosaic_authors$labels_DCs)         ~ "DCs")
col_type_a <- column_cell_type(mosaic_authors, row_type_a)

# ---- 1. Gene selection with the final settings ---------------------------------------------------------------------
mi_by_type <- split(mi_tables, mi_tables$cell_type)
directions <- lr_pairs |>
  transmute(sender = col_type_a[ligand], receiver = col_type_a[receptor],
            ligand = strip(ligand), receptor = strip(receptor))

one_side <- function(t, goi) select_genes(mi_by_type[[t]], foi = c(unique(goi), metadata), n_genes = 25)
runs <- directions |>
  group_by(sender, receiver) |>
  summarise(sender_genes   = list(one_side(first(sender), ligand)),
            receiver_genes = list(one_side(first(receiver), receptor)),
            .groups = "drop")
selection <- map(set_names(cell_types), \(t) {
  unique(c(unlist(runs$sender_genes[runs$sender == t]), unlist(runs$receiver_genes[runs$receiver == t])))
})

genes <- map(set_names(cell_types), \(t) {
  unique(c(lr_roles$gene[lr_roles$cell_type == t], setdiff(selection[[t]], metadata)))
})

# Check: the same genes as the authors' table (without its duplicates)?
authors_genes <- map(set_names(cell_types), \(t) {
  unique(strip(names(mosaic_authors)[col_type_a == t & !str_detect(names(mosaic_authors), "^(labels|condition)_")]))
})
cat("\n== Genes: ours (final settings) vs the authors' ==\n")
print(tibble(cell_type = cell_types,
             ours = lengths(genes), authors = lengths(authors_genes),
             missing = map_int(cell_types, \(t) length(setdiff(authors_genes[[t]], genes[[t]]))),
             extra   = map_int(cell_types, \(t) length(setdiff(genes[[t]], authors_genes[[t]])))))
stopifnot(all(map_lgl(cell_types, \(t) setequal(genes[[t]], authors_genes[[t]]))))

# ---- 2. Mosaic table (one column per gene), cells in the authors' order ------------------------------------------------
blocks <- map(cell_types, \(t) {
  in_type <- match_cells$cell_type == t
  cells   <- match_cells$cell[in_type]
  values  <- as.matrix(counts[genes[[t]], cells, drop = FALSE]) |> t() |> as_tibble()
  names(values) <- paste0(genes[[t]], "_", t)
  values[[paste0("labels_", t)]]    <- meta$labels[cells]
  values[[paste0("condition_", t)]] <- meta$condition[cells]
  full <- values[rep(NA_integer_, nrow(match_cells)), ]
  full[in_type, ] <- values
  full
})
mosaic <- bind_cols(blocks)
row_type <- match_cells$cell_type
col_type <- column_cell_type(mosaic, row_type)

# ---- 3. The duplicated columns of the authors' table, added back -----------------------------------------------------
dup_second <- duplicated_columns(mosaic_authors, row_type_a, col_type_a)   # the later copy of each pair
dup_keys <- paste0(strip(dup_second), "_", col_type_a[dup_second])
stopifnot(length(dup_keys) == 14, all(dup_keys %in% names(mosaic)))
cat("\n== Duplicated genes added ==\n"); print(dup_keys)

mosaic_with_dups <- mosaic
for (k in dup_keys) mosaic_with_dups[[paste0(k, "_dup")]] <- mosaic[[k]]
column_order <- unlist(map(names(mosaic), \(v) c(v, if (v %in% dup_keys) paste0(v, "_dup"))))
mosaic_with_dups <- mosaic_with_dups[, column_order]
cat("\nDimensions:", dim(mosaic_with_dups), "(the authors' table:", dim(mosaic_authors), ")\n")

# ---- 4. Check the values against the authors' columns --------------------------------------------------------------
key_authors <- tibble(variable_authors = names(mosaic_authors), cell_type = col_type_a) |>
  mutate(key = if_else(str_detect(variable_authors, "^(labels|condition)_"),
                       variable_authors, paste0(strip(variable_authors), "_", cell_type)))
cat("\n== Values: every column of the authors' table against the column of ours with the same gene ==\n")
key_authors |>
  mutate(identical = map2_lgl(variable_authors, key, \(a, k) identical(as.character(mosaic_authors[[a]]),
                                                                       as.character(mosaic[[k]])))) |>
  count(identical) |> print()

# ---- 5. State order and L-R pairs (same rules as notebook 06) -------------------------------------------------------
role_of <- lr_roles |> transmute(variable = paste0(gene, "_", cell_type), role)
group_colours <- state_authors |> distinct(group, group_color)

state_order <- tibble(var_names = names(mosaic), cell_type = col_type) |>
  left_join(role_of, by = c("var_names" = "variable")) |>
  mutate(is_metadata = str_detect(var_names, "^(labels|condition)_"),
         group = case_when(is_metadata    ~ "metadata",
                           role == "ligand"   ~ "ligands",
                           role == "receptor" ~ "receptors",
                           .default = paste(cell_type, "genes")),
         var_type = if_else(map_lgl(mosaic, is.character) | map_int(mosaic, n_distinct) <= 5, 0, 1),
         is_contextual = 0) |>
  left_join(select(state_authors, var_names, levels_increasing_order), by = "var_names") |>
  left_join(group_colours, by = "group") |>
  select(var_names, var_type, levels_increasing_order, is_contextual, group, group_color)

# a copy of a variable has the state of the original
state_dups <- state_order |> filter(var_names %in% dup_keys) |> mutate(var_names = paste0(var_names, "_dup"))
state_order <- bind_rows(state_order, state_dups) |> arrange(match(var_names, column_order))
stopifnot(identical(state_order$var_names, names(mosaic_with_dups)))

interact_edges <- lr_pairs |>
  mutate(ligand   = paste0(strip(ligand), "_", col_type_a[ligand]),
         receptor = paste0(strip(receptor), "_", col_type_a[receptor])) |>
  select(ligand, receptor)
stopifnot(all(c(interact_edges$ligand, interact_edges$receptor) %in% names(mosaic_with_dups)))

write_tsv(mosaic_with_dups, file.path(out_dir, "input.tsv"), na = "NA")
write_tsv(state_order, file.path(out_dir, "state_order.tsv"), na = "NA")
write_tsv(interact_edges, file.path(out_dir, "interact_edges.tsv"), col_names = FALSE)
write_lines(dup_keys, file.path(out_dir, "duplicated_genes.txt"))
cat("\nWritten to", out_dir, "\n")
