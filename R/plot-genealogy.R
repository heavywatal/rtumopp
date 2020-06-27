#' Visualize cell genealogy
#'
#' @details
#' `augment_genealogy` calculates coordinates of nodes and edges for plotting.
#' @param graph igraph
#' @inheritParams igraphlite::augment.Rcpp_IGraph
#' @rdname plot-genealogy
#' @export
augment_genealogy = function(graph, layout = igraphlite::layout_reingold_tilford) {
  igraphlite::augment(graph, layout = layout) %>%
    dplyr::rename(pos = "x", age = "y", parent_pos = "xend", parent_age = "yend") %>%
    dplyr::mutate(extant = .data$to %in% graph$sink)
}

#' @details
#' `gggenealogy()` creates a basic ggplot object.
#' @param data tbl from `augment_genealogy()`
#' @param mapping `ggplot2::aes()`
#' @param ... passed to `ggplot2::geom_segment()`
#' @rdname plot-genealogy
#' @export
gggenealogy = function(data, mapping = ggplot2::aes(), ...) {
  ggplot2::ggplot(data) +
    utils::modifyList(ggplot2::aes_(~age, ~pos), mapping) +
    ggplot2::geom_segment(ggplot2::aes_(xend = ~parent_age, yend = ~parent_pos), ...)
}
