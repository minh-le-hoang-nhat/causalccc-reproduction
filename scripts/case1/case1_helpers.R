# Functions shared by the case study 1 notebooks 03 and 04 (CausalCCC tutorial run with two sets of parameters).
# Usage, at the top of a notebook, after scripts/helpers.R:  source(file.path(root, "scripts", "case1", "case1_helpers.R"))
# Requires the tidyverse, Seurat, miic and causalCCC to be loaded (miic before causalCCC).

senders   <- "CD4 T"
receivers <- "CD8 T"

# ---- Fix for SeuratObject >= 5.3 --------------------------------------------------------------------------------
# causalCCC.MIselection() and causalCCC.mosaic() call GetAssayData(slot = "counts"). The `slot` argument is defunct in
# SeuratObject >= 5.3, so copies of the two functions with `layer = "counts"` are created in the global environment,
# for this R session only; the installed package is not modified.
patch_causalccc <- function() {
  for (fn_name in c("causalCCC.MIselection", "causalCCC.mosaic")) {
    fn <- get(fn_name, envir = asNamespace("causalCCC"))
    body(fn) <- str2lang(gsub('slot = "counts"', 'layer = "counts"', deparse1(body(fn), collapse = "\n"), fixed = TRUE))
    assign(fn_name, fn, envir = globalenv())
  }
}

# ---- Our cells compared with the authors' input -----------------------------------------------------------------
# For each variable of the authors' mosaic table: is there a gene of our object whose counts, in the cells of the same
# type, are exactly the authors' values? The authors' table has no cell names, so values are compared as sorted sets.
# Returns one row per variable: the authors' gene name, the matching gene of our object (NA if none).
match_authors_genes <- function(t_cells, authors_dir) {
  mosaic        <- read_tsv(file.path(authors_dir, "input.tsv"), show_col_types = FALSE)
  counts        <- GetAssayData(t_cells, assay = "RNA", layer = "counts")
  sender_cols   <- str_subset(names(mosaic), "_senders$")
  is_sender_row <- rowSums(!is.na(mosaic[, sender_cols])) > 0

  same_values <- function(gene, cells, values) identical(sort(as.numeric(counts[gene, cells])), sort(values))

  tibble(variable = names(mosaic)) |>
    mutate(cell_type = if_else(variable %in% sender_cols, senders, receivers),
           gene      = str_remove(variable, "_(senders|receivers)$")) |>
    mutate(object_gene = map2_chr(variable, cell_type, \(variable, cell_type) {
      cells  <- which(t_cells$celltype == cell_type)
      values <- mosaic[[variable]][if (cell_type == senders) is_sender_row else !is_sender_row]
      gene   <- str_remove(variable, "_(senders|receivers)$")
      if (gene %in% rownames(counts) && same_values(gene, cells, values)) return(gene)
      # Otherwise: the genes with the same total count, then the one with exactly the same values
      candidates <- rownames(counts)[Matrix::rowSums(counts[, cells]) == sum(values)]
      found <- keep(candidates, \(g) same_values(g, cells, values))
      if (length(found) == 1) found else NA_character_
    }))
}

# ---- Building the CausalCCC input files (tutorial steps) --------------------------------------------------------
# settings: list(remove_ribo, n_genes, sender_genes_in_receivers), see notebooks 03 and 04.
# The gene selection (causalCCC.MIselection, about 10 minutes per cell type) is saved in out_dir/selected_genes.rds
# and reloaded afterwards; delete that file to run it again. Writes the four files uploaded to the web server in
# out_dir/causalccc_files/, and the MI tables in out_dir/MI_tables/. Returns the mosaic table, state order and L-R pairs.
build_causalccc_input <- function(t_cells, interact_edges, settings, out_dir) {
  files_dir <- file.path(out_dir, "causalccc_files")
  mi_dir    <- file.path(out_dir, "MI_tables")
  dir.create(files_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(mi_dir, recursive = TRUE, showWarnings = FALSE)

  ligands   <- unique(interact_edges$ligands)
  receptors <- unique(interact_edges$receptors)

  # Gene selection: the genes sharing the most mutual information (MI) with the ligands (sender cells) or with the
  # receptors (receiver cells). MI is computed for each gene pair, so removing genes does not change the MI of the others.
  selected_file <- file.path(out_dir, "selected_genes.rds")
  if (file.exists(selected_file)) {
    selected <- readRDS(selected_file)
  } else {
    mi_input <- if (settings$remove_ribo) t_cells[!grepl("^Rp[sl]", rownames(t_cells)), ] else t_cells
    select_for <- function(cell_type, goi, colour) {
      causalCCC.MIselection(data_input = mi_input, assay_name = "RNA", interact_ident = "celltype",
                            oneinteract = cell_type, goi = goi,
                            n_genes = settings$n_genes %||% default_n_genes(goi),
                            metadata_list = NULL, save = TRUE, output_dir = mi_dir,
                            color_heatmap = colour, plot = TRUE)
    }
    selected <- list(senders   = select_for(senders, ligands, "darkgreen"),
                     receivers = select_for(receivers, receptors, "darkorange"))
    saveRDS(selected, selected_file)
  }

  genes_senders   <- unique(c(ligands, selected$senders))
  genes_receivers <- unique(c(receptors, selected$receivers))
  if (settings$sender_genes_in_receivers) genes_receivers <- unique(c(genes_receivers, genes_senders))

  # Mosaic table: one row per cell; sender genes are NA in receiver cells and vice versa.
  # Genes present in both lists get a _senders / _receivers suffix.
  mosaic <- causalCCC.mosaic(data_input = t_cells, assay_name = "RNA", interact_ident = "celltype",
                             senders_name = senders, receivers_name = receivers,
                             genes_senders = genes_senders, genes_receivers = genes_receivers,
                             metadata_senders = NULL, metadata_receivers = NULL)

  # State order: type (continuous / discrete) and group (ligands, receptors, senders / receivers genes) of each variable
  state_order <- causalCCC.state_order(mosaic_data_table = mosaic,
                                       genes_senders = genes_senders, genes_receivers = genes_receivers,
                                       ligands = ligands, receptors = receptors,
                                       metadata_senders = NULL, metadata_receivers = NULL)

  # L-R pairs with the same suffixed names as the mosaic table for genes present in both cell types
  edges_out <- interact_edges |>
    mutate(ligands   = if_else(ligands %in% genes_receivers, paste0(ligands, "_senders"), ligands),
           receptors = if_else(receptors %in% genes_senders, paste0(receptors, "_receivers"), receptors))

  writeLines(causalCCC.layout(state_order, network_height = 8), file.path(files_dir, "causalCCC_layout.json"))
  write.table(mosaic, file.path(files_dir, "causalCCC_df.csv"), quote = FALSE, sep = ",", row.names = FALSE)
  write.table(state_order, file.path(files_dir, "causalCCC_st.tsv"), quote = FALSE, sep = "\t", row.names = FALSE)
  write.table(edges_out, file.path(files_dir, "causalCCC_interactEdges.tsv"), quote = FALSE, sep = "\t",
              row.names = FALSE)

  list(mosaic = mosaic, state_order = state_order, interact_edges = edges_out)
}

# ---- Running MIIC with the settings of the authors' server job --------------------------------------------------
# Writes the three files in the format of the server (tab-separated mosaic, state order, L-R pairs without header) to
# out_dir/miic_input/, then runs scripts/run_miic.R (settings of the authors' job, notebook 02) unless
# out_dir/miic_input/miic_result.rds already exists. Returns the MIIC result.
# One thread: with 2 threads (the authors' server setting) MIIC crashed on our input ("C stack usage is too close to
# the limit"), as with 4 threads in case study 2. The number of threads does not change the result.
run_miic_authors_settings <- function(built, out_dir, root, n_threads = 1) {
  miic_dir <- file.path(out_dir, "miic_input")
  dir.create(miic_dir, recursive = TRUE, showWarnings = FALSE)
  result_file <- file.path(miic_dir, "miic_result.rds")
  if (!file.exists(result_file)) {
    files <- file.path(miic_dir, c("input.tsv", "state_order.tsv", "interact_edges.tsv"))
    write_tsv(built$mosaic, files[1], na = "NA")
    write_tsv(built$state_order, files[2], na = "NA")
    write_tsv(built$interact_edges, files[3], col_names = FALSE)
    status <- system2("Rscript", c(shQuote(file.path(root, "scripts", "run_miic.R")), shQuote(files),
                                   shQuote(miic_dir), n_threads))
    if (status != 0) stop("MIIC failed, see the messages above.")
  }
  readRDS(result_file)
}

# ---- Reading networks -------------------------------------------------------------------------------------------
# aliases: named vector, our gene name = the authors' name (genes renamed between annotation versions)
harmonise    <- function(genes, aliases) unname(if_else(genes %in% names(aliases), aliases[genes], genes))
strip_suffix <- function(x) str_remove(x, "_(senders|receivers)$")

# A node is a gene on one side: "B2m(S)" in sender cells, "B2m(R)" in receiver cells
read_nodes <- function(state_order, aliases = character()) {
  state_order |>
    transmute(var   = var_names,
              gene  = harmonise(strip_suffix(var_names), aliases),
              side  = if_else(group %in% c("senders genes", "ligands"), "S", "R"),
              group,
              node  = paste0(gene, "(", side, ")"))
}

# Retained links (type "P") of a MIIC edge list, with their orientation
# (infOrt: 1 = undirected, 2 = x -> y, -2 = y -> x, 6 = bidirected).
# Our miic version lists each L-R pair twice: an empty placeholder computed on 0 cells, and the link added at the end
# (notebook 02, section 4). For pairs listed twice, the placeholder (Nxy_ai = 0) is dropped.
read_network <- function(edges, nodes) {
  id <- setNames(nodes$node, nodes$var)
  edges |>
    filter(type == "P") |>
    arrange(Nxy_ai %in% 0) |>
    transmute(x = unname(id[x]), y = unname(id[y]),
              key = paste(pmin(x, y), pmax(x, y), sep = " -- "),
              orientation = case_when(infOrt ==  2 ~ paste(x, "->", y),
                                      infOrt == -2 ~ paste(y, "->", x),
                                      infOrt ==  6 ~ "bidirected",
                                      .default     = "undirected"),
              is_causal) |>
    distinct(key, .keep_all = TRUE)
}

# ---- Comparing with the authors ---------------------------------------------------------------------------------
summarise_nodes <- function(nodes, authors_nodes) {
  tibble(sender_genes     = sum(nodes$side == "S"),
         receiver_genes   = sum(nodes$side == "R"),
         ribosomal_genes  = sum(str_detect(nodes$gene, "^Rp[sl]")),
         authors_genes_found = length(intersect(nodes$node, authors_nodes$node)))
}

compare_networks <- function(net, authors_net) {
  both <- intersect(net$key, authors_net$key)
  tibble(links                  = nrow(net),
         authors_links          = nrow(authors_net),
         authors_links_found    = length(both),
         with_same_orientation  = sum(net$orientation[match(both, net$key)] ==
                                      authors_net$orientation[match(both, authors_net$key)]))
}

# Links described in the Results section of the paper (case study 1, Dupuis et al. 2025, doi:10.1093/nar/gkaf404):
# first-degree neighbours of the H2 ligands, Stat1 -> B2m, Stat1 and IFN genes, and downstream of Cd8a / Klrd1
paper_links <- tribble(
  ~a,         ~b,
  "H2-K1(S)", "B2m(S)",
  "H2-K1(S)", "Psmb8(S)",
  "H2-K1(S)", "Ifi27l2b(S)",
  "H2-K1(S)", "Slfn1(S)",
  "H2-K1(S)", "Shisa5(S)",
  "H2-D1(S)", "B2m(S)",
  "H2-D1(S)", "Irf7(S)",
  "H2-D1(S)", "Mir703(S)",
  "Stat1(S)", "B2m(S)",
  "Stat1(S)", "Ifi27l2b(S)",
  "Stat1(S)", "Isg20(S)",
  "Stat1(S)", "Irf7(S)",
  "H2-K1(S)", "Cd8a(R)",
  "H2-K1(S)", "Klrd1(R)",
  "Cd8a(R)",  "Cd8b1(R)",
  "Cd8b1(R)", "Ly6c2(R)",
  "Ly6c2(R)", "Ly6c1(R)") |>
  mutate(key = paste(pmin(a, b), pmax(a, b), sep = " -- "))

# Orientation of each paper link in each network ("absent" if the link is not retained)
check_paper_links <- function(networks) {
  out <- select(paper_links, a, b)
  for (run in names(networks)) {
    net <- networks[[run]]
    out[[run]] <- if_else(paper_links$key %in% net$key, net$orientation[match(paper_links$key, net$key)], "absent")
  }
  out
}

# Genes directly linked to `node`
neighbours <- function(net, node) sort(c(net$y[net$x == node], net$x[net$y == node]))
