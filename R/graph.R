#' Functions depending on igraph
#'
#' @details
#' `make_igraph` converts raw population tbl into graph.
#' @param population tbl
#' @rdname graph
#' @export
make_igraph = function(population) {
  as_symbolic_edgelist(population) %>%
    graph_from_symbolic_edgelist()
}

as_symbolic_edgelist = function(population) {
  is_not_origin = (population$ancestor > 0L)
  cbind(
    population[is_not_origin, "ancestor"],
    population[is_not_origin, "id"]
  )
}

# simpler version of igraph::graph_from_edgelist
graph_from_symbolic_edgelist = function(el, directed = TRUE) {
  edges = as.character(t(el))
  labels = unique(edges)
  ids = seq_along(labels)
  names(ids) = labels
  g = igraph::make_graph(ids[edges], directed = directed)
  igraph::V(g)$name = labels
  g
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
  origin = vs[igraph::degree(graph, mode = "in") == 0L]
  idx = as_idx(nodes, as_ids(vs))
  igraph::distances(graph, origin, idx, mode = "out", weights = NA, algorithm = "unweighted") %>%
    as.integer()
}

sink_nodes = function(graph) {
  vs = igraph::V(graph)
  deg = igraph::degree(graph, vs, mode = "out", loops = FALSE)
  as_ids(vs)[deg == 0L]
}

sinks = function(graph, nodes) {
  .sink = sink_nodes(graph)
  paths = paths_to_sink(graph, nodes)
  lapply(paths, function(x) x[x %in% .sink])
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
    nrow(edges) - sum(edges[,2L] %in% edges[,1L])
  }
}
