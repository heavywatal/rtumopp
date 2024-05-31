#' Functions depending on igraph
#'
#' @details
#' `make_igraph` converts raw population tbl into graph.
#' @param population tbl
#' @rdname graph
#' @export
make_igraph = function(population) {
  el = as_symbolic_edgelist(population)
  igraphlite::graph_from_data_frame(el)
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
#' @param trim whether to remove common ancestors older than MRCA.
#' @rdname graph
#' @export
subtree = function(graph, nodes = integer(0L), trim = FALSE) {
  vids = igraphlite::as_vids(graph, nodes)
  vids = upstream_vertices(graph, vids, trim = trim)
  igraphlite::induced_subgraph(graph, vids)
}

#' @details
#' `trim_root()` removes ancestral nodes older than MRCA.
#' @rdname graph
#' @export
trim_root = function(graph) {
  vn = igraphlite::Vnames(graph)
  is_branching = igraphlite::degree(graph, mode = 1L) > 1L
  mrca = min(vn[is_branching])
  older = vn < mrca
  if (any(older)) {
    graph = igraphlite::induced_subgraph(graph, igraphlite::V(graph)[!older])
  }
  graph
}

#' @details
#' `internal_nodes` selects major common ancestors above threshold.
#' @param sensitivity minimum allele frequency
#' @rdname graph
#' @export
internal_nodes = function(graph, nodes, sensitivity) {
  n = length(nodes)
  counts = paths_to_source(graph, nodes) |>
    purrr::list_c() |>
    table()
  as.integer(names(counts)[(counts / n) > sensitivity])
}

upstream_vertices = function(graph, vids, trim = FALSE) {
  vlist = neighborhood_in(graph, vids)
  vids = unique(unlist(vlist, use.names = FALSE))
  if (trim) {
    ca = Reduce(intersect, vlist)
    vids = setdiff(vids, ca[-1L]) # keep MRCA in graph
  }
  vids
}

neighborhood_out = function(graph, vids) {
  igraphlite::neighborhood(graph, vids, order = 1073741824L, mode = 1L)
}

neighborhood_in = function(graph, vids = numeric(0)) {
  igraphlite::neighborhood(graph, vids, order = 1073741824L, mode = 2L)
}

paths_to_sink = function(graph, nodes) {
  vids = igraphlite::as_vids(graph, nodes)
  vnames = igraphlite::Vnames(graph)
  res = neighborhood_out(graph, vids)
  lapply(res, function(x) vnames[x])
}

paths_to_source = function(graph, nodes = integer(0L)) {
  vids = igraphlite::as_vids(graph, nodes)
  vnames = igraphlite::Vnames(graph)
  res = neighborhood_in(graph, vids)
  lapply(res, function(x) vnames[x])
}

shortest_dist_from_source = function(graph, vids = numeric(0)) {
  res = igraphlite::distances(graph, vids, igraphlite::Vsource(graph), mode = 2L, algorithm = "unweighted")
  as.integer(res)
}

distances_from_origin = function(graph, nodes = integer(0L)) {
  vids = igraphlite::as_vids(graph, nodes)
  shortest_dist_from_source(graph, vids)
}

sink_nodes = function(graph) {
  igraphlite::Vnames(graph)[igraphlite::is_sink(graph)]
}

sinks = function(graph, nodes) {
  .sink = sink_nodes(graph)
  paths = paths_to_sink(graph, nodes)
  lapply(paths, function(x) x[x %in% .sink])
}

mean_branch_length = function(graph, from = igraphlite::Vsink(graph), to = from) {
  # TODO: Avoid creating huge matrix
  m = igraphlite::distances(graph, from, to, mode = 3L, algorithm = "unweighted")
  # Exclude self-comparison only if "within"
  n = length(m)
  if (identical(from, to)) {
    n = n - length(from)
  }
  sum(m) / n
}

# NOTE: (vcount + 1) / 2 cannot be used if death rate > 0
count_sink = function(graph, nodes = integer(0L)) {
  if (length(nodes) > 0L) {
    vids = igraphlite::as_vids(graph, nodes)
    egos = igraphlite::neighborhood(graph, vids, order = 1073741824L, mode = 1L)
    .sink = igraphlite::Vsink(graph)
    vapply(egos, function(x) sum(x %in% .sink), integer(1L), USE.NAMES = FALSE)
  } else {
    sum(igraphlite::is_sink(graph))
  }
}
