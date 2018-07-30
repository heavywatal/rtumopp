#' Sample cells from a population
#'
#' @description
#' `combn_ids` makes combinations of cell IDs for various number of samples.
#' @param x list of ID vectors
#' @param m number of regions to sample
#' @return tibble
#' @rdname sample-combn
#' @export
combn_ids = function(x, m = seq_along(x)) {
  tibble::tibble(
    nsam = m,
    id = purrr::map(m, combn_int_list, x = x)
  ) %>%
    tidyr::unnest()
}

# Make union of integers from list
union_int = function(x) {
  unique(purrr::flatten_int(x))
}

# Shortcut with different defaults
combn_int_list = function(x, m, FUN = union_int, simplify = FALSE, ...) {
  utils::combn(x, m, FUN = FUN, simplify = simplify, ...)
}

#' @description
#' `summarize_capture_rate` calculates expected allele capture rate
#' on various combinations of samples.
#' @param combinations nested tibble from `combn_ids`
#' @param population tibble
#' @inheritParams filter_common_ancestors
#' @rdname sample-combn
#' @export
summarize_capture_rate = function(combinations, population, threshold = 0.01) {
  ids = filter_common_ancestors(population, threshold = threshold)$id
  len_ids = length(ids)
  dplyr::transmute(
    combinations,
    threshold,
    .data$nsam,
    capture_rate = purrr::map_dbl(.data$id, ~sum(.x %in% ids)) / len_ids
  )
}
