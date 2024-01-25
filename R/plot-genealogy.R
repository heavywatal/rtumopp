#' Visualize cell genealogy
#'
#' @details
#' `augment_genealogy` calculates coordinates of nodes and edges for plotting.
#' @param graph igraph
#' @param mu mutation rate per cell division.
#' @param accel assumption of accelerated mutation.
#' @rdname plot-genealogy
#' @export
augment_genealogy = function(graph, mu = 0, accel = 0) {
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
  if (mu > 0 || accel > 0) {
    named_dist = genetic_distance(graph, igraphlite::V(graph), mu = mu, accel = accel)
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

genetic_distance = function(graph, vids = igraphlite::V(graph), mu = 0, accel = 0) {
  age = distances_from_origin(graph)
  segments = (1 + accel)**age
  if (mu > 0) {
    segments = stats::rpois(graph$vcount, mu * segments)
  }
  segments[igraphlite::is_source(graph)] = 0
  ancestors = neighborhood_in(graph, vids)
  names(ancestors) = igraphlite::Vnames(graph)[vids]
  purrr::map_dbl(ancestors, \(v) sum(segments[v]))
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
