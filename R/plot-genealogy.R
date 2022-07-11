#' Visualize cell genealogy
#'
#' @details
#' `augment_genealogy` calculates coordinates of nodes and edges for plotting.
#' @param graph igraph
#' @rdname plot-genealogy
#' @export
augment_genealogy = function(graph) {
  vnames = igraphlite::Vnames(graph)
  tips = vnames[graph$is_sink]
  root = vnames[graph$is_source]
  igraphlite::augment(graph, layout = igraphlite::layout_reingold_tilford) |>
    dplyr::rename(
      parent = "from", node = "to",
      d = "y", x_parent = "xend", d_parent = "yend"
    ) |>
    dplyr::mutate(node_type = dplyr::case_when(
      .data$node %in% tips ~ "tip",
      .data$node %in% root ~ "root",
      TRUE ~ "internal"
    ))
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
    utils::modifyList(ggplot2::aes_(~d, ~x), mapping) +
    ggplot2::geom_segment(ggplot2::aes_(xend = ~d_parent, yend = ~x_parent), ...)
}

# @param distance named vector from [genetic_distance()]
mutate_distance = function(data, distance) {
  data$d = distance[as.character(data$node)]
  data$d_parent = distance[as.character(data$parent)]
  data
}

genetic_distance = function(graph, vids = graph$V, mu = 0, accel = 0) {
  age = distances_from_origin(graph)
  segments = (1 + accel)**age
  if (mu > 0) {
    segments = stats::rpois(graph$vcount, mu * segments)
  }
  segments[graph$is_source] = 0
  ancestors = neighborhood_in(graph, vids)
  names(ancestors) = igraphlite::Vnames(graph)[vids]
  purrr::map_dbl(ancestors, function(x) {
    sum(segments[x])
  })
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
    dplyr::select(-.data$idx) |>
    dplyr::relocate(.data$x, .data$x_parent, .after = .data$d_parent) |>
    dplyr::arrange(.data$node)
}
