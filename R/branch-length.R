#' Calculate genetic distance with igraph
#'
#' @description
#' `mean_branch_length` calculates mean branch length within/between sub-graphs.
#' @param graph igraph
#' @param from,to igraph vertices
#' @rdname branch-length
#' @export
mean_branch_length = function(graph, from = igraph::V(graph), to = from) {
  .d = igraph::distances(graph, from, to, mode = "all", weights = NA, algorithm = "unweighted")
  .n = length(from) * length(to) - sum(from %in% to)
  sum(.d) / .n
}

#' @description
#' `within_between_samples` summarizes branch lengths.
#' @param regions output of `sample_uniform_regions()`
#' @rdname branch-length
#' @export
within_between_samples = function(graph, regions) {
  regions %>%
    tibble::rowid_to_column() %>%
    dplyr::mutate(
      id = lapply(.data$id, as_idx, vs = igraph::V(graph)),
      within = purrr::map_dbl(.data$id, ~mean_branch_length(graph, .x))
    ) %>%
    purrr::transpose() %>%
    purrr::cross2(., ., .filter = ~.x$rowid >= .y$rowid) %>%
    purrr::map_dfr(~{
      row_i = .x[[1L]]
      row_j = .x[[2L]]
      tibble::tibble(
        region_i = row_i$rowid,
        region_j = row_j$rowid,
        within_i = row_i$within,
        within_j = row_j$within,
        euclidean = dist_euclidean(row_i, row_j),
        between = mean_branch_length(graph, row_i$id, row_j$id)
      )
    }) %>%
    dplyr::mutate(
      within = 0.5 * (.data$within_i + .data$within_j),
      fst = fst_HBK(.data$within, .data$between)
    )
}

# Translate cell ID to igraph vertex index
as_idx = function(cell_id, vs) {
  match(cell_id, as.integer(as_ids(vs)))
}

# Kst by Hudson, Boos, and Kaplan (1992).
fst_HBK = function(within, between, n=2) {
  (between - within) / (between + within / (n - 1))
}
