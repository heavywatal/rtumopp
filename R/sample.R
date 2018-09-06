#' Sample cells from a population
#'
#' @description
#' `sample_uniform_regions` sample uniformly distributed regions.
#' @param tbl tibble of extant cells with id, x, y, z columns
#' @param nsam number of regions to sample
#' @param ncell number of cells per specimen
#' @param jitter amount of random variations on x and y axes
#' @rdname sample
#' @export
sample_uniform_regions = function(tbl, nsam = 2L, ncell = 10L, jitter = 0) {
  centers = kmeans_centers(tbl, nsam)
  if (jitter > 0) {
    centers = mutate_jitter(centers, .data$x, .data$y, amount = jitter)
  }
  sample_regions(tbl, centers, ncell = ncell)
}

#' @description
#' `sample_random_regions` samples multiple regions at random.
#' @rdname sample
#' @export
sample_random_regions = function(tbl, nsam = 2L, ncell = 10L) {
  centers = tbl %>%
    dplyr::select(.data$x, .data$y, .data$z) %>%
    dplyr::sample_n(nsam)
  sample_regions(tbl, centers, ncell = ncell)
}

#' @description
#' `sample_bulk` samples a bulk of cells near the specified center.
#' @param center named (x, y, z) vector, list, or tibble
#' @rdname sample
#' @export
sample_bulk = function(tbl, center = c(x = 0, y = 0, z = 0), ncell = 10L) {
  d = dist_euclidean(tbl, center)
  tbl$id[utils::head(order(d), ncell)]
}

kmeans_centers = function(tbl, centers, iter.max = 32L) {
  tbl %>%
    dplyr::select(.data$x, .data$y, .data$z) %>%
    stats::kmeans(centers = centers, iter.max = iter.max) %>%
    purrr::pluck("centers") %>%
    tibble::as_tibble()
}

sample_regions = function(tbl, centers, ncell = 10L) {
  dplyr::mutate(centers, id = purrr::pmap(centers, function(x, y, z) {
    sample_bulk(tbl, center = c(x = x, y = y, z = z), ncell = ncell)
  }))
}

tidy_regions = function(regions) {
  tibble::tibble(id = regions$id) %>%
    tibble::rowid_to_column(var = "region") %>%
    tidyr::unnest()
}

mutate_jitter = function(.data, ..., amount) {
  #TODO receive values from ... instead of amount
  dplyr::mutate_at(.data, dplyr::vars(...), function(x) {
    x + stats::runif(length(x), -amount, amount)
  })
}

#' @description
#' `evaluate_mrs` is a shortcut to evaluate multi-region sampling.
#' @param population tbl
#' @inheritParams filter_common_ancestors
#' @inheritParams internal_nodes
#' @rdname sample
#' @export
evaluate_mrs = function(population, nsam, ncell, threshold = 0.05, sensitivity = 0.05, jitter = 0) {
  graph = make_igraph(population)
  regions = filter_extant(population) %>%
    sample_uniform_regions(nsam = nsam, ncell = ncell, jitter = jitter)
  sampled = purrr::flatten_int(regions$id)
  detectable = internal_nodes(graph, sampled, sensitivity = sensitivity)
  major_ca = population %>%
    add_node_property(graph) %>%
    filter_common_ancestors(threshold = threshold) %>%
    dplyr::pull("id")
  sum(detectable %in% major_ca) / length(major_ca)
}
