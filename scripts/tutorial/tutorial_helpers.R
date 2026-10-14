# Functions of the CausalCCC tutorial (notebooks/tutorial/causalccc_tutorial.qmd).
# Usage, at the top of the notebook, after scripts/helpers.R:
#   source(file.path(root, "scripts", "tutorial", "tutorial_helpers.R"))
# Requires the tidyverse, Seurat, miic and causalCCC to be loaded (miic before causalCCC).

# ---- Fix for SeuratObject >= 5.3 --------------------------------------------------------------------------------
# causalCCC.MIselection() and causalCCC.mosaic() read the counts with GetAssayData(slot = "counts"); the `slot`
# argument is defunct since SeuratObject 5.3. This function makes copies of the two functions with layer = "counts",
# and a copy of causalCCC.wrapper() (with a second fix, see below) that calls these copies (an R function looks up the functions it calls in its own
# environment, so a patched copy of the wrapper needs an environment that holds the patched copies).
# The copies are assigned in the global environment, for this R session only; the installed package is unchanged.
patch_causalccc_all <- function() {
  ns        <- asNamespace("causalCCC")
  patch_env <- new.env(parent = ns)
  for (fn_name in c("causalCCC.MIselection", "causalCCC.mosaic")) {
    fn <- get(fn_name, envir = ns)
    body(fn) <- str2lang(gsub('slot = "counts"', 'layer = "counts"', deparse1(body(fn), collapse = "\n"), fixed = TRUE))
    environment(fn) <- ns
    assign(fn_name, fn, envir = patch_env)
  }
  # Second fix, in the wrapper: with two or more metadata, MIselection() can return a metadata column among the
  # selected "genes" (for example labels, which shares a lot of MI with condition); causalCCC.mosaic() then looks for
  # it in the count matrix and stops ("column not found: [labels]"). The metadata are removed from the selected genes.
  wrapper <- get("causalCCC.wrapper", envir = ns)
  code <- deparse1(body(wrapper), collapse = "\n")
  for (side in c("senders", "receivers")) {
    old <- sprintf("MI_%s_genes, goi_%s", side, side)
    if (!grepl(old, code, fixed = TRUE)) stop("causalCCC.wrapper() code has changed: patch not applied")
    code <- gsub(old, sprintf("setdiff(MI_%s_genes, names(metadata_%s)), goi_%s", side, side, side), code, fixed = TRUE)
  }
  body(wrapper) <- str2lang(code)
  environment(wrapper) <- patch_env
  assign("causalCCC.wrapper", wrapper, envir = patch_env)
  for (fn_name in ls(patch_env)) assign(fn_name, get(fn_name, envir = patch_env), envir = globalenv())
  invisible(ls(patch_env))
}

# ---- Raw counts from log-normalised values ----------------------------------------------------------------------
# norm: genes x cells sparse matrix (dgCMatrix) of log(1 + scale * count / library size). In each cell the smallest
# non-zero value is a count of 1, which gives the library size, then every count (notebook case2/05, section 1).
recover_counts <- function(norm, scale = 1e4) {
  cell_of_value <- rep(seq_len(ncol(norm)), diff(norm@p))
  library_size  <- as.numeric(scale / expm1(tapply(norm@x, cell_of_value, min)))
  counts   <- norm
  counts@x <- expm1(norm@x) * library_size[cell_of_value] / scale
  attr(counts, "max_rounding_error") <- max(abs(counts@x - round(counts@x)))
  counts@x <- round(counts@x)
  counts
}

# ---- Gene selection: one run per sender -> receiver direction ---------------------------------------------------
# What causalCCC.wrapper() does for one direction, applied to precomputed MI tables (columns variables, foi, value,
# cell_type) and repeated for every direction of `directions` (columns sender, receiver, ligand, receptor).
#   cap: number of ligands / receptors used as features of interest (the wrapper uses the first 8)
#   n_genes: function of the genes of interest giving the number of genes (default: the formula of MIselection)
#   drop: regular expression of genes removed from the candidates before the selection (NULL: none)
# Returns, for each cell type, the selected genes pooled over the directions in which it takes part.
select_by_direction <- function(mi_tables, directions, meta_sender = NULL, meta_receiver = NULL, cap = 8,
                                n_genes = default_n_genes, drop = NULL) {
  mi_by_type <- split(mi_tables, mi_tables$cell_type)
  if (!is.null(drop)) mi_by_type <- map(mi_by_type, \(mi) filter(mi, !str_detect(variables, drop)))
  one_side <- function(t, goi, meta) {
    goi <- head(unique(goi), cap)
    select_genes(mi_by_type[[t]], foi = c(goi, meta), n_genes = n_genes(goi))
  }
  runs <- directions |>
    group_by(sender, receiver) |>
    summarise(sender_genes   = list(one_side(first(sender), ligand, meta_sender)),
              receiver_genes = list(one_side(first(receiver), receptor, meta_receiver)),
              .groups = "drop")
  cell_types <- unique(c(directions$sender, directions$receiver))
  map(set_names(cell_types), \(t) {
    setdiff(unique(c(unlist(runs$sender_genes[runs$sender == t]), unlist(runs$receiver_genes[runs$receiver == t]))),
            c("labels", "condition"))
  })
}

# Selected genes compared with a reference selection, per cell type. `exclude` (named list per cell type): genes not
# counted (the ligands and receptors, which are always kept).
score_selection <- function(selected, reference, exclude = list()) {
  imap(reference, \(ref, t) {
    ours <- setdiff(selected[[t]], exclude[[t]])
    ref  <- setdiff(ref, exclude[[t]])
    tibble(cell_type = t, selected = length(ours), shared = length(intersect(ours, ref)),
           missing = length(setdiff(ref, ours)), extra = length(setdiff(ours, ref)),
           jaccard = length(intersect(ours, ref)) / length(union(ours, ref)))
  }) |> list_rbind()
}

# ---- Mosaic table and state order for any number of cell types --------------------------------------------------
# counts: genes x cells; meta: one row per cell with columns cell_type and the metadata in `metadata`;
# genes: named list (cell type -> genes). Every variable is named <gene>_<cell type> (or <metadata>_<cell type>), so
# each gene appears exactly once per cell type, whatever the number of directions it was selected in.
# Rows: the cells of each type, in a random order (as causalCCC.mosaic() does; set the seed before calling).
build_mosaic <- function(counts, meta, genes, metadata = c("labels", "condition")) {
  cell_types <- names(genes)
  cells_of   <- map(set_names(cell_types), \(t) sample(which(meta$cell_type == t)))
  n_rows     <- sum(lengths(cells_of))
  offset     <- set_names(cumsum(c(0, lengths(cells_of)))[seq_along(cells_of)], cell_types)   # rows before each type
  blocks <- imap(cells_of, \(cells, t) {
    values <- as_tibble(t(as.matrix(counts[genes[[t]], cells, drop = FALSE])))
    names(values) <- paste0(genes[[t]], "_", t)
    for (m in metadata) values[[paste0(m, "_", t)]] <- as.character(meta[[m]][cells])
    full <- values[rep(NA_integer_, n_rows), ]
    full[offset[[t]] + seq_along(cells), ] <- values
    full
  })
  mosaic <- bind_cols(blocks)
  attr(mosaic, "row_type") <- rep(cell_types, lengths(cells_of))
  mosaic
}

# State order with the rules of causalCCC.state_order(): metadata discrete; a gene discrete if it has 5 or fewer
# distinct values *counting NA*; groups ligands / receptors / "<cell type> genes" / metadata.
# roles: tibble(variable, role) for the ligands and receptors; levels: named list (metadata -> ordered levels).
build_state_order <- function(mosaic, roles, levels = list()) {
  tibble(var_names = names(mosaic)) |>
    left_join(roles, by = c("var_names" = "variable")) |>
    mutate(metadata  = str_extract(var_names, paste0("^(", paste(names(levels), collapse = "|"), ")(?=_)")),
           cell_type = str_remove(var_names, "^.*?_"),
           group = case_when(!is.na(metadata)  ~ "metadata",
                             role == "ligand"   ~ "ligands",
                             role == "receptor" ~ "receptors",
                             .default = paste(cell_type, "genes")),
           var_type = if_else(map_lgl(mosaic, is.character) | map_int(mosaic, n_distinct) <= 5, 0, 1),
           levels_increasing_order = map_chr(metadata, \(m) if (is.na(m)) NA_character_ else
                                               paste(levels[[m]], collapse = ",")),
           is_contextual = 0) |>
    select(var_names, var_type, levels_increasing_order, is_contextual, group)
}

# ---- MIIC runs --------------------------------------------------------------------------------------------------
# The settings of the authors' server jobs (notebooks case1/02 and case2/04); `...` overrides any of them.
authors_miic_settings <- list(cplx = "nml", orientation = TRUE, latent = "orientation", ori_proba_ratio = 0.01,
                              n_shuffles = 100, conf_threshold = 0.01, test_mar = TRUE, consistent = "no",
                              n_eff = -1, mode = "S", n_threads = 1)

# Runs miic() unless `file` exists, then returns the saved result. interact_edges must have the columns `ligand` and
# `receptor` (miic reads interact_edges$ligand and interact_edges$receptor).
run_miic_cached <- function(file, mosaic, state_order, interact_edges, seed = 2026, ...) {
  if (file.exists(file)) return(readRDS(file))
  settings <- modifyList(authors_miic_settings, list(...))
  set.seed(seed)
  start <- Sys.time()
  # capture.output() hides the progress bars printed by the C++ core
  invisible(capture.output(suppressMessages(
    result <- do.call(miic, c(list(input_data = as.data.frame(mosaic), state_order = as.data.frame(state_order),
                                   interact_edges = as.data.frame(interact_edges)), settings)))))
  attr(result, "run_time") <- Sys.time() - start
  attr(result, "settings") <- settings
  # miic() replaces an invalid value (for example cplx = "bic") by its default with only a warning: check that the
  # settings it used are those requested
  used <- c(cplx = result$params$cplx, latent = result$params$latent, consistent = result$params$consistent)
  asked <- unlist(settings[names(used)])
  if (!identical(unname(used), unname(asked)))
    stop("miic() did not use the requested settings: asked ", paste(names(asked), asked, collapse = ", "),
         "; used ", paste(names(used), used, collapse = ", "))
  saveRDS(result, file)
  result
}

# Retained links of a MIIC result, one row per link, without the placeholder rows written for each L-R pair
# (Nxy_ai = 0, notebook case1/02, section 4). Orientation recoded for genes sorted alphabetically:
# 2 = first -> second, -2 = second -> first, 1 = not oriented, 6 = bidirected.
retained_links <- function(miic_result, interact_edges) {
  lr <- link_name(interact_edges$ligand, interact_edges$receptor)
  as_tibble(miic_result$all.edges.summary) |>
    filter(type == "P", !(link_name(x, y) %in% lr & Nxy_ai %in% 0)) |>
    mutate(link = link_name(x, y), is_lr = link %in% lr,
           orientation = if_else(x == pmin(x, y) | infOrt %in% c(1, 6), infOrt, -infOrt)) |>
    distinct(link, .keep_all = TRUE)
}

# One row per run: size of the network (L-R links excluded) and agreement with the reference run.
compare_runs <- function(runs, reference, interact_edges) {
  ref <- retained_links(runs[[reference]], interact_edges) |> filter(!is_lr)
  imap(runs, \(run, name) {
    net    <- retained_links(run, interact_edges) |> filter(!is_lr)
    common <- inner_join(net, ref, by = "link", suffix = c("", "_ref"))
    minutes <- round(as.numeric(attr(run, "run_time"), units = "mins"), 1)   # before tibble(): its column `run` would hide the result
    tibble(run = name,
           minutes          = minutes,
           links            = nrow(net),
           not_oriented     = sum(net$infOrt == 1),
           oriented         = sum(abs(net$infOrt) == 2),
           bidirected       = sum(net$infOrt == 6),
           genuine_causes   = sum(net$is_causal == "Y"),
           jaccard_links    = nrow(common) / length(union(net$link, ref$link)),
           same_orientation = mean(common$orientation == common$orientation_ref),
           same_causal      = mean(common$is_causal == common$is_causal_ref))
  }) |> list_rbind()
}
