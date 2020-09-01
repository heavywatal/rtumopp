#' Utilities for variant allele frequency
#'
#' @details
#' `make_vaf` is a shortcut to make neutral VAF pattern.
#' @inheritParams subtree
#' @inheritParams mutate_clades
#' @rdname vaf
#' @export
make_vaf = function(graph, samples, mu, threshold = 0.05) {
  mutated = subtree(graph, purrr::flatten_int(samples)) %>%
    mutate_clades(mu = mu)
  tally_vaf(samples, mutated$carriers) %>%
    filter_detectable(threshold) %>%
    sort_vaf() %>%
    tidy_vaf()
}

#' @details
#' `tally_vaf` evaluates overlap of sampled and mutated cells.
#' @param samples list of integer IDs
#' @param sites list of integer IDs
#' @rdname vaf
#' @export
tally_vaf = function(samples, sites) {
  names(samples) = seq_along(samples)
  lapply(samples, function(region) {
    vapply(sites, function(holders) {
      sum(region %in% holders)
    }, integer(1), USE.NAMES = FALSE) / length(region)
  }) %>% tibble::as_tibble()
}

#' @details
#' `tidy_vaf` transforms vaf table.
#' @param tbl output of `tally_vaf`
#' @rdname vaf
#' @export
tidy_vaf = function(tbl) {
  tbl %>%
    tibble::rowid_to_column(var = "site") %>%
    tidyr::gather("sample", "frequency", -"site") %>%
    dplyr::filter(.data$frequency > 0) %>%
    dplyr::mutate(sample = as.integer(.data$sample))
}

#' @details
#' `filter_detectable` removes sites where freq < threshold.
#' @param threshold minimum detectable frequency
#' @rdname vaf
#' @export
filter_detectable = function(tbl, threshold) {
  tbl[tbl < threshold] = 0
  dplyr::filter(tbl, Reduce(`+`, tbl) > 0)
}

#' @details
#' `sort_vaf` reorders rows and columns of VAF table.
#' @param method passed to `stats::hclust`
#' @rdname vaf
#' @export
sort_vaf = function(tbl, method = c("average", "ward.D2", "complete", "single")) {
  method = match.arg(method)
  if (nrow(tbl) < 2L) {
    return(tbl)
  }
  rsums = Reduce(`+`, tbl > 0)
  w_nonzero = 10
  w_shared = ifelse(rsums < 2L, 0, 1000)
  covr_issue377 = function(.x) {
    ifelse(.x > 0, w_nonzero, 0) + w_shared + .x
  }
  tbl_weighted = dplyr::mutate(tbl, dplyr::across(dplyr::everything(), covr_issue377))
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

#' @details
#' `dist_vaf` summarizes diversity within/between regions.
#' @inheritParams sample_uniform_regions
#' @rdname vaf
#' @export
#' @seealso [dist_genealogy()]
dist_vaf = function(tbl, ncell = Inf) {
  n_regions = ncol(tbl)
  seq_n = seq_len(n_regions)
  mtrx = matrix(0, n_regions, n_regions, dimnames = list(seq_n, seq_n))
  for (i in seq_n) {
    for (j in seq.int(i, n_regions)) {
      p_i = tbl[, i, drop = TRUE]
      p_j = tbl[, j, drop = TRUE]
      mtrx[i, j] = mtrx[j, i] = mean(p_hetero(p_i, p_j))
    }
  }
  diag(mtrx) = diag(mtrx) / (1 - 1 / ncell)
  mtrx
}

p_hetero = function(p1, p2 = p1) {
  p1 * (1 - p2) + (1 - p1) * p2
}
