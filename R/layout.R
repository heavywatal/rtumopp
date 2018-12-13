#' Layout cell genealogy
#'
#' @details
#' `layout_genealogy` returns coordinates of nodes and edges for plotting.
#' @param population tbl
#' @rdname layout
#' @export
layout_genealogy = function(population) {
  extra_cols = population %>%
    dplyr::transmute(name = as.character(.data$id), extant = .data$death == 0, .data$clade)
  make_igraph(population) %>%
    igraph_layout(igraph::as_tree(flip.y = FALSE)) %>%
    dplyr::rename(pos = "x", age = "y", posend = "xend", ageend = "yend") %>%
    dplyr::left_join(extra_cols, by = c(to = "name"))
}

igraph_layout = function(graph, layout = igraph::nicely(), ...) {
  nodes = igraph_layout_nodes(graph, layout, ...)
  to_nodes = dplyr::rename(nodes, xend = "x", yend = "y")
  igraph::as_data_frame(graph, "edges") %>%
    tibble::as_tibble() %>%
    dplyr::left_join(nodes, by = c(from = "name")) %>%
    dplyr::left_join(to_nodes, by = c(to = "name"))
}

igraph_layout_nodes = function(graph, layout = igraph::nicely(), ...) {
  igraph::layout_(graph, layout, ...) %>%
    magrittr::set_colnames(c("x", "y")) %>%
    tibble::as_tibble() %>%
    dplyr::mutate(name = igraph::V(graph)$name)
}
