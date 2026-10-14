# Helper functions shared by the notebooks (from notebook 03 on).
# Notebooks 01 and 02 define the same functions inline; they are gathered here unchanged
# or generalised (to more than two cell types) so that the case study 2 notebooks can reuse them.
# Usage, at the top of a notebook:  source(file.path(root, "scripts", "helpers.R"))
# Requires the tidyverse to be loaded.

# The repository root is the folder that contains both README.md and notebooks/
find_repo_root <- function(dir = getwd()) {
  while (!(file.exists(file.path(dir, "README.md")) && dir.exists(file.path(dir, "notebooks")))) {
    parent <- dirname(dir)
    if (parent == dir) stop("Repository root not found. Open the repository folder in Positron.")
    dir <- parent
  }
  dir
}

# A name for a link that does not depend on which gene is written as x and which as y
link_name <- function(a, b) paste(pmin(a, b), pmax(a, b), sep = " | ")

# Genes directly linked to `node` in an edge list with columns x, y, infOrt, is_causal
neighbors_of <- function(node, net) {
  net |>
    filter(x == node | y == node) |>
    transmute(node, neighbor = if_else(x == node, y, x), infOrt, is_causal,
              # direction seen from `node`: "out" = node -> neighbor, "in" = neighbor -> node
              direction = case_when(infOrt %in% c(1, 6)                 ~ if_else(infOrt == 1, "none", "bidirected"),
                                    (infOrt == 2) == (x == node)         ~ "out",
                                    .default                             = "in"))
}

# Number of links between `from` and every gene it can reach, ignoring directions
# (breadth-first search: neighbours at step 1, their neighbours at step 2, ...)
graph_distance <- function(net, from) {
  distance <- setNames(0, from)
  frontier <- from
  step     <- 0
  repeat {
    touching  <- net |> filter(x %in% frontier | y %in% frontier)
    new_genes <- setdiff(unique(c(touching$x, touching$y)), names(distance))
    if (length(new_genes) == 0) break
    step     <- step + 1
    distance <- c(distance, setNames(rep(step, length(new_genes)), new_genes))
    frontier <- new_genes
  }
  distance
}

# Cell type of each column of a mosaic table: the cell type whose rows hold its values.
# `row_type` gives the cell type of each row.
column_cell_type <- function(mosaic, row_type) {
  map_chr(mosaic, \(column) {
    types <- unique(row_type[!is.na(column)])
    if (length(types) == 1) types else paste(sort(types), collapse = "+")
  })
}

# Retained links in a common format: genes sorted alphabetically, orientation recoded accordingly
# (2 = first -> second, -2 = second -> first, 1 = not oriented, 6 = bidirected).
# `cell_type` is a named vector giving the cell type of each variable.
standardise <- function(edge_list, cell_type) {
  edge_list |>
    filter(type == "P") |>
    mutate(first       = pmin(x, y),
           second      = pmax(x, y),
           link        = link_name(x, y),
           orientation = if_else(x == first | infOrt %in% c(1, 6), infOrt, -infOrt),
           type_first  = cell_type[first],
           type_second = cell_type[second],
           category    = if_else(type_first == type_second,
                                 paste("within", type_first),
                                 "ligand-receptor")) |>
    select(link, first, second, category, orientation, is_causal, info_shifted, confidence)
}

# End of causalCCC.MIselection(), reproduced line by line (notebook 05, section 4), applied to a precomputed MI table
# (columns variables, foi, value): the 3 genes with the highest MI with each feature of interest, plus the n_genes
# genes with the highest MI with any feature of interest, keeping MI above bit_threshold.
select_genes <- function(mi, foi, n_genes, bit_threshold = 0.001) {
  mi <- mi |> filter(foi %in% !!foi)
  # top 3 genes for each feature of interest, in the order of the features
  top3 <- map(foi, \(f) {
    mi |> filter(foi == f, value > bit_threshold) |> arrange(desc(value)) |> head(3) |> pull(variables)
  }) |> unlist()
  # genes ranked by their highest MI with any feature of interest
  top_n <- mi |>
    group_by(variables) |>
    slice_max(value, n = 1, with_ties = FALSE) |>
    ungroup() |>
    filter(value > bit_threshold) |>
    arrange(desc(value)) |>
    head(n_genes) |>
    pull(variables)
  unique(c(top3, top_n))
}

# Default number of genes of causalCCC.MIselection()
default_n_genes <- function(goi) max(15, min(length(goi) * 4, 100))

# The MI table of causalCCC.MIselection(return_full = TRUE) for one cell type, computed faster (notebook 10).
# Each block of genes goes through the same miic() call as in MIselection (causalCCC commit 0aa8586, l. 556-630), so
# the values are identical, but:
# - the counts stay sparse and only one block of genes at a time is made dense (MIselection makes the whole matrix
#   dense, then copies it again for every feature of interest);
# - genes never detected in these cells are not computed: MIselection gives them 0;
# - blocks of block_size genes instead of 100, so fewer miic() calls (the MI of a gene with a feature does not depend on
#   the other genes of its block, notebook 07);
# - with cache_dir (one folder per cell type), the table of each feature is saved, and reused by later calls.
# counts: genes x cells sparse matrix of raw counts; meta: one row per cell, same order; cells: the cells of the type
# (logical or indices); goi: genes of interest; metadata: columns of meta used both as features of interest and as
# variables, like metadata_list in MIselection. Requires miic to be loaded.
# Returns a tibble (variables, foi, value) with the rows of MIselection(return_full = TRUE) in the same order.
compute_mi_table <- function(counts, meta, cells, goi, metadata = NULL, block_size = 1000, n_threads = 8,
                             cache_dir = NULL) {
  counts    <- counts[, cells, drop = FALSE]
  meta      <- meta[cells, , drop = FALSE]
  metadata  <- metadata[!map_lgl(metadata, \(m) all(is.na(meta[[m]])))]   # as MIselection, skip all-NA metadata
  variables <- c(rownames(counts), metadata)
  detected  <- c(rownames(counts)[Matrix::rowSums(counts != 0) > 0], metadata)
  if (!is.null(cache_dir)) dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  one_feature <- function(f) {
    cache_file <- if (!is.null(cache_dir)) file.path(cache_dir, paste0(make.names(f), ".rds"))
    if (!is.null(cache_file) && file.exists(cache_file)) return(readRDS(cache_file))

    values <- if (f %in% metadata) meta[[f]] else as.numeric(counts[f, ])
    todo   <- setdiff(detected, f)
    mi     <- setNames(rep(0, length(variables)), variables)
    for (block in split(todo, ceiling(seq_along(todo) / block_size))) {
      genes   <- setdiff(block, metadata)
      df_loop <- as.data.frame(as.matrix(Matrix::t(counts[genes, , drop = FALSE])))
      for (m in intersect(block, metadata)) df_loop[[m]] <- meta[[m]]
      df_loop <- df_loop[, block, drop = FALSE]
      df_loop$onefoi <- values
      # as in MIselection: discrete if not numeric or 5 distinct values or fewer; genes are consequences only, so
      # miic computes the gene-feature pairs only
      df_st <- data.frame(var_names      = names(df_loop),
                          var_type       = map_dbl(df_loop, \(v) if (!is.numeric(v) || n_distinct(v) <= 5) 0 else 1),
                          is_contextual  = 0,
                          is_consequence = c(rep(1, length(block)), 0))
      invisible(capture.output(capture.output(suppressMessages(suppressWarnings(
        res <- miic(input_data = df_loop, state_order = df_st, orientation = FALSE, latent = "no",
                    n_threads = n_threads))), type = "message")))
      res <- res$all.edges.summary
      mi[if_else(res$x == "onefoi", res$y, res$x)] <- round(res$info_shifted, 6)
    }
    mi[is.na(mi)] <- 0
    table <- tibble(variables = variables, foi = f, value = unname(mi))
    if (!is.null(cache_file)) saveRDS(table, cache_file)
    table
  }

  map(unique(c(goi, metadata)), one_feature) |> list_rbind()
}

# ---- Drawing networks -------------------------------------------------------------------------------------------
# Three functions: network_data() chooses the nodes and links to draw; plot_network_interactive() draws them with
# visNetwork (in the HTML notebooks); plot_network_static() draws them with ggplot2 (PNG files for slides).
# Conventions in both drawings: colour = cell type; diamonds = centres; dashed = given L-R link; arrowhead =
# orientation found by MIIC (infOrt); green = genuine causation; thicker = stronger link (info_shifted).

cell_type_colours <- c(Fibroblasts = "#b5d99c", T_cell = "#f2c57c", DCs = "#9cc3e6")

# Label of a variable in the drawings: gene without cell-type suffixes, "subtype" for the labels metadata,
# and a star for the variables listed in `copies` (e.g. the duplicated columns of the authors' table)
gene_label <- function(v, copies = character()) {
  label <- str_remove_all(v, "_(Fibroblasts|T_cell|DCs)")
  label <- if_else(label == "labels", "subtype", label)
  if_else(v %in% copies, paste0(label, "*"), label)
}

# Columns of a mosaic table identical to an earlier column of the same cell type (notebook 03, section 3)
duplicated_columns <- function(mosaic, row_type, col_type) {
  unlist(map(unique(col_type), \(t) {
    block <- mosaic[row_type == t, col_type == t]
    names(block)[duplicated(as.list(block))]
  }))
}
cell_type_borders <- c(Fibroblasts = "#4f7a28", T_cell = "#a8742a", DCs = "#2f6491")

# edges: MIIC edge list (x, y, type, infOrt, is_causal, info_shifted); lr_links: link_name() of the L-R pairs;
# cell_type_of: named vector giving the cell type of each variable; label_of: function turning a variable name into
# its label. With centres = NULL: the whole network. Otherwise: the variables within max_distance intracellular links
# of a centre, the retained links among them, and the L-R links between centres (like Figure 3C of the paper).
network_data <- function(edges, lr_links, cell_type_of, label_of = identity, centres = NULL, max_distance = 2) {
  retained <- edges |>
    filter(type == "P") |>
    mutate(is_lr = link_name(x, y) %in% lr_links)

  if (is.null(centres)) {
    keep  <- unique(c(retained$x, retained$y))
    links <- retained
  } else {
    intracellular <- retained |> filter(!is_lr)
    keep <- unique(unlist(map(centres, \(centre) {
      d <- graph_distance(intracellular, centre)
      names(d)[d <= max_distance]
    })))
    links <- retained |>
      filter(x %in% keep, y %in% keep) |>
      filter(!is_lr | (x %in% centres & y %in% centres))
  }

  strongest <- max(links$info_shifted[!links$is_lr])
  list(
    nodes = tibble(id = keep, label = label_of(keep), cell_type = unname(cell_type_of[keep]),
                   is_centre = keep %in% centres),
    links = links |>
      transmute(from = x, to = y, infOrt, is_causal, is_lr, info_shifted,
                # 0 to 1; L-R links drawn thin, since their strength is set artificially (notebook 04, section 6)
                strength = if_else(is_lr, 0.15, pmin(info_shifted / strongest, 1)))
  )
}

plot_network_interactive <- function(net, title = NULL, seed = 1, height = "600px") {
  nodes <- net$nodes |>
    transmute(id, label, group = cell_type, title = id,
              shape = if_else(is_centre, "diamond", "dot"),
              size  = if_else(is_centre, 25, 12))
  edges <- net$links |>
    transmute(from, to,
              arrows = case_when(infOrt == 2 ~ "to", infOrt == -2 ~ "from", infOrt == 6 ~ "to;from", .default = ""),
              dashes = is_lr,
              color  = if_else(is_causal == "Y", "#2e7d32", "#9e9e9e"),
              width  = 1 + 5 * strength,
              title  = paste0(from, " — ", to, "<br>strength (info_shifted): ", round(info_shifted),
                              "<br>orientation (infOrt): ", infOrt, "<br>genuine causation: ", is_causal))

  plot <- visNetwork::visNetwork(nodes, edges, main = title, width = "100%", height = height)
  for (t in names(cell_type_colours)) {
    plot <- visNetwork::visGroups(plot, groupname = t,
                                  color = list(background = cell_type_colours[[t]], border = cell_type_borders[[t]]))
  }
  plot <- plot |>
    visNetwork::visLegend(position = "right", main = "cell type") |>
    visNetwork::visOptions(highlightNearest = list(enabled = TRUE, degree = 1))
  # Large networks: a precomputed layout (igraph) instead of the slow physics simulation of the browser
  if (nrow(nodes) > 80) {
    visNetwork::visIgraphLayout(plot, layout = "layout_with_fr", randomSeed = seed)
  } else {
    visNetwork::visLayout(plot, randomSeed = seed)
  }
}

plot_network_static <- function(net, title = NULL, seed = 1, label_size = 3) {
  graph <- igraph::graph_from_data_frame(select(net$links, from, to), directed = FALSE,
                                         vertices = select(net$nodes, id))
  set.seed(seed)
  xy <- igraph::layout_with_fr(graph)
  xy <- (xy - min(xy)) / (max(xy) - min(xy))   # scale to 0-1
  pos <- tibble(id = igraph::V(graph)$name, px = xy[, 1], py = xy[, 2])

  # Segments shortened at both ends so that arrowheads stop at the edge of the nodes;
  # every oriented link is drawn from its tail to its head
  gap <- 0.012
  segments <- net$links |>
    left_join(rename(pos, from = id, x0 = px, y0 = py), by = "from") |>
    left_join(rename(pos, to = id, x1 = px, y1 = py), by = "to") |>
    mutate(swap = infOrt == -2,
           xa = if_else(swap, x1, x0), ya = if_else(swap, y1, y0),
           xb = if_else(swap, x0, x1), yb = if_else(swap, y0, y1),
           len = sqrt((xb - xa)^2 + (yb - ya)^2), ux = (xb - xa) / len, uy = (yb - ya) / len,
           xa = xa + gap * ux, ya = ya + gap * uy, xb = xb - gap * ux, yb = yb - gap * uy,
           kind = case_when(abs(infOrt) == 2 ~ "oriented", infOrt == 6 ~ "bidirected", .default = "none"),
           colour = if_else(is_causal == "Y", "genuine causation", "other"),
           linetype = if_else(is_lr, "L-R link (given)", "learnt link"))

  nodes <- net$nodes |> left_join(pos, by = "id")
  head_one  <- grid::arrow(length = grid::unit(1.6, "mm"), type = "closed")
  head_both <- grid::arrow(length = grid::unit(1.6, "mm"), type = "closed", ends = "both")
  segment_aes <- aes(x = xa, y = ya, xend = xb, yend = yb, colour = colour, linetype = linetype, linewidth = strength)

  ggplot() +
    geom_segment(data = filter(segments, kind == "none"), segment_aes, alpha = 0.7) +
    geom_segment(data = filter(segments, kind == "oriented"), segment_aes, arrow = head_one) +
    geom_segment(data = filter(segments, kind == "bidirected"), segment_aes, arrow = head_both) +
    geom_point(data = nodes, aes(x = px, y = py, fill = cell_type, shape = is_centre, size = is_centre),
               colour = "grey20", stroke = 0.4) +
    ggrepel::geom_text_repel(data = nodes, aes(x = px, y = py, label = label), size = label_size,
                             max.overlaps = Inf, segment.colour = "grey70", segment.size = 0.2, seed = seed) +
    scale_fill_manual(values = cell_type_colours, name = "cell type") +
    scale_shape_manual(values = c("FALSE" = 21, "TRUE" = 23), guide = "none") +
    scale_size_manual(values = c("FALSE" = 2.5, "TRUE" = 5), guide = "none") +
    scale_colour_manual(values = c("genuine causation" = "#2e7d32", "other" = "grey60"), name = NULL) +
    scale_linetype_manual(values = c("L-R link (given)" = "dashed", "learnt link" = "solid"), name = NULL) +
    scale_linewidth(range = c(0.2, 1.2), guide = "none") +
    guides(fill = guide_legend(override.aes = list(shape = 21, size = 4))) +
    coord_equal(clip = "off") +
    labs(title = title) +
    theme_void() +
    theme(plot.background = element_rect(fill = "white", colour = NA), legend.position = "bottom",
          plot.title = element_text(face = "bold"))
}

# Shortcut kept for notebook 06: interactive drawing of the neighbourhood of `centres`
plot_neighbourhood <- function(edges, lr_links, centres, cell_type_of, label_of = identity,
                               max_distance = 2, title = NULL, seed = 1) {
  network_data(edges, lr_links, cell_type_of, label_of, centres, max_distance) |>
    plot_network_interactive(title = title, seed = seed)
}
