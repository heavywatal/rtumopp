#' Calculate genetic distance with igraph
#'
#' @details
#' `mean_branch_length` calculates mean branch length within/between sub-graphs.
#' @param graph igraph
#' @param from,to igraph vertices (not cell ID)
#' @rdname branch-length
#' @export
mean_branch_length = function(graph, from = graph$sink, to = from) {
  # TODO: Avoid creating huge matrix
  m = igraphlite::shortest_paths(graph, from, to, mode = 3L, algorithm = "unweighted")
  if (length(from)) {
    nzero = sum(from %in% to)
  } else {
    nzero = nrow(m)
  }
  sum(m) / (nrow(m) * ncol(m) - nzero)
}

#' @details
#' `within_between_samples` summarizes branch lengths.
#' @param regions output of `sample_uniform_regions()`
#' @rdname branch-length
#' @export
within_between_samples = function(graph, regions) {
  rows = regions %>%
    tibble::rowid_to_column() %>%
    dplyr::mutate(
      id = lapply(.data$id, igraphlite::as_vids, graph = graph),
      within = purrr::map_dbl(.data$id, ~ mean_branch_length(graph, .x))
    ) %>%
    purrr::transpose()
  purrr::cross2(rows, rows, .filter = ~ .x$rowid >= .y$rowid) %>%
    purrr::map_dfr(~ {
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

# Kst by Hudson, Boos, and Kaplan (1992).
fst_HBK = function(within, between, n = 2) {
  (between - within) / (between + within / (n - 1))
}

Kst = function(graph, regions) {
  vids = lapply(regions$id, igraphlite::as_vids, graph = graph)
  t_within = purrr::map_dbl(vids, ~ mean_branch_length(graph, .x))
  t_total = mean_branch_length(graph, unique(unlist(vids)))
  1 - mean(t_within) / t_total
}
