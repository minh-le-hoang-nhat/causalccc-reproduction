# Test 1 of the benchmarks (notebook 10): are the outputs of a modified run identical to the reference run?
# Functions, sourced by notebook 10, and a command line:
#   Rscript scripts/benchmarks/compare_outputs.R <reference miic_result.rds> <new miic_result.rds>
# Requires dplyr (loaded by the command line; the notebook loads the tidyverse).

# One row per link of a MIIC edge summary, genes sorted alphabetically and the orientation recoded accordingly
# (2 = first -> second, -2 = second -> first, 1 = not oriented, 6 = bidirected). Each L-R pair is listed twice in this
# miic version, once as an empty placeholder computed on 0 cells (notebook 02, section 4): the placeholder is dropped.
miic_links <- function(miic_result) {
  miic_result$all.edges.summary |>
    as_tibble() |>
    mutate(first = pmin(x, y), second = pmax(x, y),
           link = paste(first, second, sep = " | "),
           orientation = if_else(x == first | infOrt %in% c(1, 6), infOrt, -infOrt)) |>
    arrange(link, Nxy_ai %in% 0) |>
    distinct(link, .keep_all = TRUE) |>
    select(link, type, orientation, is_causal, info_shifted, confidence)
}

# Largest absolute difference between two numeric vectors (NA on both sides counts as equal)
max_diff <- function(a, b) {
  both_na <- is.na(a) & is.na(b)
  if (any(is.na(a) != is.na(b))) return(Inf)
  if (all(both_na)) return(0)
  max(abs(a[!both_na] - b[!both_na]))
}

# Compares two MIIC results. Returns a one-row tibble with the verdict:
#   "identical"                    same retained links, orientations, genuine causations and values;
#   "same decisions, values differ" same links, orientations and genuine causations, some values differ;
#   "different network"            otherwise (see the Jaccard index and the counts of changes).
compare_miic <- function(reference, new) {
  ref <- miic_links(reference)
  cmp <- miic_links(new)
  kept_ref <- ref$link[ref$type == "P"]
  kept_cmp <- cmp$link[cmp$type == "P"]
  shared   <- inner_join(filter(ref, type == "P"), filter(cmp, type == "P"), by = "link", suffix = c("_ref", "_new"))
  all_rows <- inner_join(ref, cmp, by = "link", suffix = c("_ref", "_new"))

  out <- tibble(
    links_reference        = length(kept_ref),
    links_new              = length(kept_cmp),
    links_jaccard          = length(intersect(kept_ref, kept_cmp)) / length(union(kept_ref, kept_cmp)),
    orientations_changed   = sum(shared$orientation_ref != shared$orientation_new),
    causal_changed         = sum(shared$is_causal_ref != shared$is_causal_new),
    max_diff_info_shifted  = max_diff(all_rows$info_shifted_ref, all_rows$info_shifted_new),
    max_diff_confidence    = max_diff(all_rows$confidence_ref, all_rows$confidence_new),
    links_confidence_diff  = sum(!mapply(identical, all_rows$confidence_ref, all_rows$confidence_new)))

  same_decisions <- setequal(kept_ref, kept_cmp) && out$orientations_changed == 0 && out$causal_changed == 0
  same_values    <- out$max_diff_info_shifted == 0 && out$links_confidence_diff == 0
  out$verdict <- if (same_decisions && same_values) "identical" else
    if (same_decisions) "same decisions, values differ" else "different network"
  relocate(out, verdict)
}

# Compares two MI tables (columns variables, foi, value, cell_type), and the genes select_genes() (scripts/helpers.R)
# picks from each, for every cell type with all its features of interest and n_genes = 25 (the exact settings of
# notebook 07) or the default number. Returns a one-row tibble; verdict "identical" when every value is equal.
compare_mi <- function(reference, new) {
  both <- full_join(reference, new, by = c("cell_type", "foi", "variables"), suffix = c("_ref", "_new"))
  same_genes <- map_lgl(unique(reference$cell_type), \(t) {
    foi <- unique(reference$foi[reference$cell_type == t])
    goi <- setdiff(foi, c("labels", "condition"))
    all(map_lgl(list(25, default_n_genes(goi)), \(n) identical(
      select_genes(filter(reference, cell_type == t), foi, n),
      select_genes(filter(new, cell_type == t), foi, n))))
  })
  tibble(rows_reference        = nrow(reference),
         rows_new              = nrow(new),
         rows_missing          = sum(is.na(both$value_new)) + sum(is.na(both$value_ref)),
         max_diff_value        = max_diff(both$value_ref, both$value_new),
         values_differing      = sum(both$value_ref != both$value_new, na.rm = TRUE),
         selected_genes_same   = all(same_genes)) |>
    mutate(verdict = if_else(rows_missing == 0 & max_diff_value == 0 & selected_genes_same, "identical",
                             if_else(selected_genes_same, "same selection, values differ", "different selection"))) |>
    relocate(verdict)
}

if (sys.nframe() == 0) {
  suppressPackageStartupMessages(library(dplyr))
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 2) stop("Usage: Rscript scripts/benchmarks/compare_outputs.R <reference.rds> <new.rds>")
  print(as.data.frame(compare_miic(readRDS(args[1]), readRDS(args[2]))), row.names = FALSE)
}
