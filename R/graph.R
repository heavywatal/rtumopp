#' Functions depending on igraph
#'
#' @description
#' `make_igraph` converts raw population tbl into graph
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

#' `subtree` extracts subgraph among terminal nodes
#' @param graph igraph
#' @param nodes igraph vertices
#' @rdname graph
#' @export
subtree = function(graph, nodes=character(0L)) {
  paths_to_source(graph, nodes) %>%
    purrr::flatten_chr() %>%
    unique() %>%
    {igraph::induced_subgraph(graph, .)}
}

paths_to_sink = function(graph, nodes) {
  igraph::ego(graph, order = 1073741824L, nodes = nodes, mode = "out") %>%
    purrr::map(names)
}

paths_to_source = function(graph, nodes=character(0L)) {
  igraph::ego(graph, order = 1073741824L, nodes = nodes, mode = "in") %>%
    purrr::map(names)
}

distances_from_origin = function(graph, nodes=character(0L)) {
  igraph::distances(graph, "1", nodes, mode = "out", weights = NA, algorithm = "unweighted") %>%
    as.integer()
}

#' `layout_genealogy` returns coordinates of nodes and edges for plotting
#' @rdname graph
#' @export
layout_genealogy = function(population) {
  extra_cols = population %>%
    dplyr::transmute(name = as.character(.data$id), extant = .data$death == 0, .data$clade)
  make_igraph(population) %>%
    wtl::igraph_layout(igraph::as_tree(flip.y = FALSE)) %>%
    dplyr::rename(pos = "x", age = "y", posend = "xend", ageend = "yend") %>%
    dplyr::left_join(extra_cols, by = c(to = "name"))
}
