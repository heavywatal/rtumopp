#' Functions to deal with driver mutations
#'
#' @details
#' `propagate_drivers` calculates the effect of drivers in descendant cells
#' @param drivers data.frame
#' @param graph igraph
#' @rdname drivers
#' @export
propagate_drivers = function(drivers, graph) {
  drivers |>
    dplyr::mutate(id = paths_to_sink(graph, .data$id)) |>
    dplyr::mutate(coef = 1 + .data$coef) |>
    tidyr::unnest("id") |>
    dplyr::group_by(.data$id, .data$type) |>
    dplyr::summarize(value = prod(!!as.name("coef"))) |>
    tidyr::pivot_wider(names_from = "type", values_from = "value")
}
