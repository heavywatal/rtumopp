#' Utilities for variant allele frequency
#'
#' @description
#' `make_vaf` is a shortcut to make neutral VAF pattern.
#' @inheritParams subtree
#' @inheritParams mutate_clades
#' @rdname vaf
#' @export
make_vaf = function(graph, samples, mu, threshold = 0.05) {
  mutated = subtree(graph, purrr::flatten_chr(samples)) %>%
    mutate_clades(mu = mu) %>%
    purrr::map(as.integer)
  tally_vaf(samples, mutated) %>%
    filter_detectable(threshold) %>%
    sort_vaf() %>%
    tidy_vaf()
}

#' @description
#' `tally_vaf` evaluates overlap of sampled and mutated cells.
#' @param samples list of integer IDs
#' @param sites list of integer IDs
#' @rdname vaf
#' @export
tally_vaf = function(samples, sites) {
  purrr::map_dfc(samples, function(region) {
    purrr::map_int(sites, function(holders) {
      sum(region %in% holders)
    }) / length(region)
  })
}

#' @description
#' `tidy_vaf` transforms vaf table.
#' @param tbl output of `tally_vaf`
#' @rdname vaf
#' @export
tidy_vaf = function(tbl) {
  col_names = seq_len(ncol(tbl))
  tbl %>%
    rlang::set_names(col_names) %>%
    tibble::rowid_to_column(var = "site") %>%
    tidyr::gather("sample", "frequency", -"site") %>%
    dplyr::filter(.data$frequency > 0) %>%
    dplyr::mutate(sample = as.integer(.data$sample))
}

#' @description
#' `filter_detectable` removes sites where freq < threshold.
#' @param threshold minimum detectable frequency
#' @rdname vaf
#' @export
filter_detectable = function(tbl, threshold) {
  tbl[tbl < threshold] = 0
  dplyr::filter(tbl, rowSums(tbl) > 0)
}

#' @description
#' `sort_vaf` reorders rows and columns of VAF table.
#' @param method passed to `stats::hclust`
#' @rdname vaf
#' @export
sort_vaf = function(tbl, method = c("average", "ward.D2", "complete", "single")) {
  method = match.arg(method)
  if (nrow(tbl) < 2L) return(tbl)
  rsums = rowSums(tbl > 0)
  w_nonzero = 10
  w_shared = ifelse(rsums < 2L, 0, 1000)
  tbl_weighted = dplyr::mutate_all(tbl, ~{
    ifelse(.x > 0, w_nonzero, 0) + w_shared + .x
  })
  d_rows = stats::dist(tbl_weighted, method = "euclidean")
  order_rows = stats::hclust(d_rows, method = method)$order
  tbl_shared = dplyr::filter(tbl_weighted, rsums > 1L)
  if (nrow(tbl_shared) > 1L) {
    d_cols = stats::dist(t(tbl_shared), method = "euclidean")
    order_cols = stats::hclust(d_cols, method = method)$order
    tbl[order_rows, order_cols]
  } else {
    tbl[order_rows, ]
  }
}
