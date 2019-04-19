#' Functions depending on igraph
#'
#' @details
#' `make_igraph` converts raw population tbl into graph.
#' @param population tbl
#' @rdname graph
#' @export
make_igraph = function(population) {
  as_symbolic_edgelist(population) %>%
    igraphlite::graph_from_data_frame()
}

as_symbolic_edgelist = function(population) {
  is_not_origin = (population[["ancestor"]] > 0L)
  cbind(
    population[is_not_origin, "ancestor"],
    population[is_not_origin, "id"]
  )
}

#' @details
#' `subtree` extracts subgraph among terminal nodes.
#' @param graph igraph
#' @param nodes integer cell IDs
#' @rdname graph
#' @export
subtree = function(graph, nodes = integer(0L)) {
  vids = igraphlite::as_vids(graph, nodes)
  vids = upstream_vertices(graph, vids, to_mrca = TRUE)
  igraphlite::induced_subgraph(graph, vids)
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

upstream_vertices = function(graph, vids, to_mrca = TRUE) {
  vlist = igraphlite::neighborhood(graph, vids, order = 1073741824L, mode = 2L)
  vids = unique(unlist(vlist, use.names = FALSE))
  if (to_mrca) {
    vids = setdiff(vids, Reduce(intersect, vlist)[-1L])
  }
  vids
}

paths_to_sink = function(graph, nodes) {
  vids = igraphlite::as_vids(graph, nodes)
  vnames = igraphlite::Vnames(graph)
  igraphlite::neighborhood(graph, vids, order = 1073741824L, mode = 1L) %>%
    lapply(function(x) vnames[x])
}

paths_to_source = function(graph, nodes = integer(0L)) {
  vids = igraphlite::as_vids(graph, nodes)
  vnames = igraphlite::Vnames(graph)
  igraphlite::neighborhood(graph, vids, order = 1073741824L, mode = 2L) %>%
    lapply(function(x) vnames[x])
}

distances_from_origin = function(graph, nodes = integer(0L)) {
  vids = igraphlite::as_vids(graph, nodes)
  origin = graph$source
  igraphlite::shortest_paths(graph, origin, vids, mode = 1L, algorithm = "unweighted") %>%
    as.integer()
}

sink_nodes = function(graph) {
  igraphlite::Vnames(graph)[graph$is_sink]
}

sinks = function(graph, nodes) {
  .sink = sink_nodes(graph)
  paths = paths_to_sink(graph, nodes)
  lapply(paths, function(x) x[x %in% .sink])
}

# NOTE: (vcount + 1) / 2 cannot be used if death rate > 0
count_sink = function(graph, nodes = integer(0L)) {
  if (length(nodes) > 0L) {
    vids = igraphlite::as_vids(graph, nodes)
    egos = igraphlite::neighborhood(graph, vids, order = 1073741824L, mode = 1L)
    sink = graph$sink
    vapply(egos, function(x) sum(x %in% sink), integer(1L), USE.NAMES = FALSE)
  } else {
    sum(graph$is_sink)
  }
}
