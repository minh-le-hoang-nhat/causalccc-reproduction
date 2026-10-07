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

# Interactive drawing of the neighbourhood of some nodes, like Figure 3C of the paper (requires visNetwork).
# edges: MIIC edge list (x, y, type, infOrt, is_causal, info_shifted); lr_links: link_name() of the L-R pairs;
# centres: the nodes around which to draw; cell_type_of: named vector giving the cell type of each variable;
# label_of: function turning a variable name into the label shown.
# Nodes: all variables within max_distance intracellular links of a centre (diamonds = centres), coloured by cell type.
# Links: all retained links among these nodes, plus the L-R links between centres (dashed). Arrowheads follow MIIC's
# orientation (infOrt); genuine causations in green; width grows with the strength of the link (info_shifted).
plot_neighbourhood <- function(edges, lr_links, centres, cell_type_of, label_of = identity,
                               max_distance = 2, title = NULL, seed = 1) {
  retained      <- edges |> filter(type == "P")
  intracellular <- retained |> filter(!link_name(x, y) %in% lr_links)
  keep <- unique(unlist(map(centres, \(centre) {
    d <- graph_distance(intracellular, centre)
    names(d)[d <= max_distance]
  })))

  links <- retained |>
    filter(x %in% keep, y %in% keep) |>
    mutate(is_lr = link_name(x, y) %in% lr_links) |>
    filter(!is_lr | (x %in% centres & y %in% centres))
  strongest <- max(links$info_shifted[!links$is_lr])

  nodes <- tibble(id = keep) |>
    mutate(label = label_of(id),
           group = unname(cell_type_of[id]),
           shape = if_else(id %in% centres, "diamond", "dot"),
           size  = if_else(id %in% centres, 25, 12),
           title = id)

  vis_edges <- links |>
    transmute(from = x, to = y,
              arrows = case_when(infOrt == 2 ~ "to", infOrt == -2 ~ "from", infOrt == 6 ~ "to;from", .default = ""),
              dashes = is_lr,
              color  = if_else(is_causal == "Y", "#2e7d32", "#9e9e9e"),
              width  = 1 + 5 * pmin(info_shifted / strongest, 1),
              title  = paste0(x, " — ", y, "<br>strength (info_shifted): ", round(info_shifted),
                              "<br>orientation (infOrt): ", infOrt, "<br>genuine causation: ", is_causal))

  visNetwork::visNetwork(nodes, vis_edges, main = title, width = "100%", height = "600px") |>
    visNetwork::visGroups(groupname = "Fibroblasts", color = list(background = "#b5d99c", border = "#4f7a28")) |>
    visNetwork::visGroups(groupname = "T_cell", color = list(background = "#f2c57c", border = "#a8742a")) |>
    visNetwork::visGroups(groupname = "DCs", color = list(background = "#9cc3e6", border = "#2f6491")) |>
    visNetwork::visLegend(position = "right", main = "cell type") |>
    visNetwork::visOptions(highlightNearest = list(enabled = TRUE, degree = 1)) |>
    visNetwork::visLayout(randomSeed = seed)
}
