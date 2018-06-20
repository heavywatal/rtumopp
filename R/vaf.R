#' Utilities for variant allele frequency
#'
#' @description
#' `tally_vaf` evaluates overlap of sampled and mutated cells
#' @param samples list of integer IDs
#' @param sites list of integer IDs
#' @rdname vaf
#' @export
tally_vaf = function(samples, sites) {
  purrr::map_dfc(samples, function(region) {
    purrr::map_int(sites, function(holders) {
      sum(region %in% holders)
    })
  })
}

#' `tidy_vaf` transforms vaf table
#' @param tbl output of `tally_vaf`
#' @rdname vaf
#' @export
tidy_vaf = function(tbl) {
  tbl %>%
    rlang::set_names(seq_len(ncol(.))) %>%
    tibble::rowid_to_column(var = "site") %>%
    tidyr::gather("sample", "frequency", -"site") %>%
    dplyr::filter(.data$frequency > 0) %>%
    dplyr::mutate(sample = as.integer(.data$sample))
}
