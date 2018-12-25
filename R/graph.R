#' Functions depending on igraph
#'
#' @details
#' `make_igraph` converts raw population tbl into graph.
#' @param population tbl
#' @rdname graph
#' @export
make_igraph = function(population) {
  population %>%
    dplyr::transmute(
      from = .data$ancestor,
      to = .data$id
    ) %>%
    dplyr::filter(.data$from > 0L) %>%
    dplyr::arrange(.data$to) %>%
    dplyr::mutate_all(as.character) %>%
    igraph::graph_from_data_frame()
}

#' @details
#' `subtree` extracts subgraph among terminal nodes.
#' @param graph igraph
#' @param nodes integer cell IDs
#' @rdname graph
#' @export
subtree = function(graph, nodes = integer(0L)) {
  vs = igraph::V(graph)
  idx = as_idx(nodes, as_ids(vs))
  idx = igraph::ego(graph, order = 1073741824L, nodes = idx, mode = "in") %>%
    purrr::flatten_dbl() %>%
    unique()
  igraph::induced_subgraph(graph, idx)
}

#' @details
#' `internal_nodes` selects major common ancestors above threshold.
#' @param sensitivity minimum allele frequency
#' @rdname graph
#' @export
internal_nodes = function(graph, nodes, sensitivity) {
  n = length(nodes)
  counts = paths_to_source(graph, nodes) %>%
    purrr::flatten_int() %>%
    table()
  as.integer(names(counts)[(counts / n) > sensitivity])
}

as_ids = function(vs) {
  ns = names(vs)
  if (is.null(ns)) as.vector(vs) else as.integer(ns)
}

as_idx = function(nodes, ids) {
  idx = match(nodes, ids)
  is_na = is.na(idx)
  if (any(is_na)) {
    warning("Node not found: ", paste(nodes[is_na], collapse = ", "))
    idx = idx[!is_na]
  }
  idx
}

paths_to_sink = function(graph, nodes) {
  vs = igraph::V(graph)
  ids = as_ids(vs)
  idx = as_idx(nodes, ids)
  igraph::ego(graph, order = 1073741824L, nodes = idx, mode = "out") %>%
    lapply(function(x) ids[x])
}

paths_to_source = function(graph, nodes = integer(0L)) {
  vs = igraph::V(graph)
  ids = as_ids(vs)
  idx = as_idx(nodes, ids)
  igraph::ego(graph, order = 1073741824L, nodes = idx, mode = "in") %>%
    lapply(function(x) ids[x])
}

distances_from_origin = function(graph, nodes = integer(0L)) {
  vs = igraph::V(graph)
  idx = as_idx(nodes, as_ids(vs))
  igraph::distances(graph, 1, idx, mode = "out", weights = NA, algorithm = "unweighted") %>%
    as.integer()
}

# NOTE: (vcount + 1) / 2 cannot be used if death rate > 0
count_sink = function(graph, nodes = integer(0L)) {
  edges = igraph::as_edgelist(graph, names = FALSE)
  if (length(nodes) > 0L) {
    vs = igraph::V(graph)
    idx = as_idx(nodes, as_ids(vs))
    egos = igraph::ego(graph, order = 1073741824L, nodes = idx, mode = "out")
    sink = edges[!edges[,2L] %in% edges[,1L], 2L]
    vapply(egos, function(x) {sum(x %in% sink)}, integer(1L), USE.NAMES = FALSE)
  } else {
    NROW(edges) - sum(edges[,2L] %in% edges[,1L])
  }
}
