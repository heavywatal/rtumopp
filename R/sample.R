#' Sample cells from a population
#'
#' @description
#' `sample_uniform_regions` sample uniformly distributed regions
#' @param tbl tibble of extant cells with id, x, y, z columns
#' @param nsam number of regions to sample
#' @param ncell number of cells per specimen
#' @rdname sample
#' @export
sample_uniform_regions = function(tbl, nsam=2L, ncell=10L) {
  tbl %>%
    dplyr::select(.data$id, .data$x, .data$y, .data$z) %>%
    sample_regions(kmeans_centers(., nsam), ncell = ncell)
}

#' `sample_random_regions` samples multiple regions at random
#' @rdname sample
#' @export
sample_random_regions = function(tbl, nsam=2L, ncell=10L) {
  centers = tbl %>%
    dplyr::select(.data$x, .data$y, .data$z) %>%
    dplyr::sample_n(nsam)
  sample_regions(tbl, centers, ncell = ncell)
}

#' `sample_bulk` samples a bulk of cells near the specified center
#' @param center named (x, y, z) vector, list, or tibble
#' @rdname sample
#' @export
sample_bulk = function(tbl, center=c(x=0, y=0, z=0), ncell=10L) {
  d = dist_euclidean(tbl, center)
  tbl$id[utils::head(order(d), ncell)]
}

kmeans_centers = function(tbl, centers, iter.max = 32L) {
  tbl %>%
    dplyr::select(.data$x, .data$y, .data$z) %>%
    stats::kmeans(centers = centers, iter.max = iter.max) %>%
    {tibble::as_tibble(.$centers)}
}

sample_regions = function(tbl, centers, ncell = 10L) {
  centers %>%
    dplyr::mutate(id = purrr::pmap(., function(x, y, z) {
      sample_bulk(tbl, center = c(x = x, y = y, z = z), ncell = ncell)
    }))
}

tidy_regions = function(regions) {
  tibble::tibble(id = regions$id) %>%
    tibble::rowid_to_column(var = "region") %>%
    tidyr::unnest()
}
