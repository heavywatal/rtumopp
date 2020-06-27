#' Layout cell genealogy
#'
#' @details
#' `layout_genealogy` returns coordinates of nodes and edges for plotting.
#' @param population tbl
#' @rdname layout
#' @export
layout_genealogy = function(population) {
  extra_cols = population %>%
    dplyr::transmute(name = .data$id, extant = .data$death == 0, .data$clade)
  make_igraph(population) %>%
    igraphlite::augment(igraphlite::layout_reingold_tilford) %>%
    dplyr::rename(pos = "x", age = "y", posend = "xend", ageend = "yend") %>%
    dplyr::left_join(extra_cols, by = c(to = "name"))
}
