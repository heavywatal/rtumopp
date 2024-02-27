#' Utilities for variant allele frequency
#'
#' @details
#' `make_vaf` generates neutral VAF pattern from cell genealogy.
#' @param samples list of integer IDs
#' @inheritParams subtree
#' @inheritParams make_sample
#' @rdname vaf
#' @export
make_vaf = function(graph, samples, mu) {
  df_mut = mutate_clades(graph, mu = mu)
  names(samples) = seq_along(samples)
  purrr::map(samples, function(sampled) {
    freqs = vapply(df_mut$carriers, function(carriers) {
      sum(carriers %in% sampled)
    }, integer(1), USE.NAMES = FALSE)
    rep(freqs / length(sampled), times = df_mut$number)
  }) |> tibble::as_tibble()
}

#' @details
#' `make_longer_vaf` is a shortcut to make VAF in longer format.
#' @rdname vaf
#' @export
make_longer_vaf = function(graph, samples, mu, threshold = 0.05) {
  if (sum(igraphlite::is_sink(graph)) > sum(lengths(samples))) {
    graph = subtree(graph, unlist(samples))
  }
  make_vaf(graph, samples, mu) |>
    filter_detectable(threshold) |>
    sort_vaf() |>
    longer_vaf()
}

#' @details
#' `longer_vaf` transforms vaf table.
#' @param vaf output of [make_vaf()]
#' @rdname vaf
#' @export
longer_vaf = function(vaf) {
  vaf |>
    tibble::rowid_to_column(var = "site") |>
    tidyr::pivot_longer(!"site", names_to = "sample", values_to = "frequency")
}

#' @details
#' `filter_detectable` removes sites where freq < threshold.
#' @param threshold minimum detectable frequency
#' @rdname vaf
#' @export
filter_detectable = function(vaf, threshold) {
  vaf[vaf < threshold] = 0
  dplyr::filter(vaf, Reduce(`+`, vaf) > 0)
}

#' @details
#' `sort_vaf` reorders rows and columns of VAF table.
#' @param method passed to `stats::hclust`
#' @rdname vaf
#' @export
sort_vaf = function(vaf, method = c("average", "ward.D2", "complete", "single")) {
  method = match.arg(method)
  if (nrow(vaf) < 2L) {
    return(vaf)
  }
  rsums = Reduce(`+`, vaf > 0)
  w_nonzero = 10
  w_shared = ifelse(rsums < 2L, 0, 1000)
  tbl_weighted = dplyr::mutate(vaf, dplyr::across(dplyr::everything(), \(.x) {
    ifelse(.x > 0, w_nonzero, 0) + w_shared + .x
  }))
  d_rows = stats::dist(tbl_weighted, method = "euclidean")
  order_rows = stats::hclust(d_rows, method = method)$order
  tbl_shared = dplyr::filter(tbl_weighted, rsums > 1L)
  if (nrow(tbl_shared) > 1L) {
    d_cols = stats::dist(t(tbl_shared), method = "euclidean")
    order_cols = stats::hclust(d_cols, method = method)$order
    vaf[order_rows, order_cols]
  } else {
    vaf[order_rows, ]
  }
}
