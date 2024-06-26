#' Visualize cell genealogy
#'
#' @details
#' `augment_genealogy` calculates coordinates of nodes and edges for plotting.
#' @param graph igraph
#' @param lengths lengths of edges. [edge_lengths()] is useful for this.
#' The edge attribute "weight" is used if `TRUE`.
#' @rdname plot-genealogy
#' @export
augment_genealogy = function(graph, lengths = numeric(0)) {
  vnames = igraphlite::Vnames(graph)
  tips = vnames[igraphlite::is_sink(graph)]
  root = vnames[igraphlite::is_source(graph)]
  res = igraphlite::augment(graph, layout = igraphlite::layout_reingold_tilford) |>
    dplyr::rename(
      parent = "from", node = "to",
      d = "y", x_parent = "xend", d_parent = "yend"
    ) |>
    dplyr::mutate(node_type = dplyr::case_when(
      .data$node %in% tips ~ "tip",
      .data$node %in% root ~ "root",
      TRUE ~ "internal"
    ))
  if (length(lengths) > 0) {
    udist = distances_upstream(graph, weights = lengths)
    named_dist = stats::setNames(udist, vnames)
    res$d = named_dist[as.character(res$node)]
    res$d_parent = named_dist[as.character(res$parent)]
  }
  res
}

#' @details
#' `gggenealogy()` creates a basic ggplot object.
#' @param data tbl from [augment_genealogy()]
#' @param mapping [ggplot2::aes()]
#' @param ... passed to `ggplot2::geom_segment()`
#' @rdname plot-genealogy
#' @export
gggenealogy = function(data, mapping = ggplot2::aes(), ...) {
  ggplot2::ggplot(data) +
    utils::modifyList(ggplot2::aes(.data[["d"]], .data[["x"]]), mapping) +
    ggplot2::geom_segment(ggplot2::aes(xend = .data[["d_parent"]], yend = .data[["x_parent"]]), ...)
}

remove_tandem_ancestors = function(data) {
  data |>
    dplyr::group_by(.data$x, .data$x_parent) |>
    dplyr::summarize(
      idx = which.max(!!as.name("d")),
      parent = (!!as.name("parent"))[which.min(!!as.name("d"))],
      node = (!!as.name("node"))[!!as.name("idx")],
      d = max(!!as.name("d")),
      d_parent = min(!!as.name("d_parent")),
      node_type = (!!as.name("node_type"))[!!as.name("idx")]
    ) |>
    dplyr::ungroup() |>
    dplyr::select(!"idx") |>
    dplyr::relocate("x", "x_parent", .after = "d_parent") |>
    dplyr::arrange(.data$node)
}
