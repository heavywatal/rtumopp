#' Sample cells from a population
#'
#' @details
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
    centers = dplyr::mutate(centers, dplyr::across(c("x", "y"), function(x) {
      x + stats::runif(length(x), -jitter, jitter)
    }))
  }
  sample_regions(tbl, centers, ncell = ncell)
}

#' @details
#' `sample_random_regions` samples multiple regions at random.
#' @rdname sample
#' @export
sample_random_regions = function(tbl, nsam = 2L, ncell = 10L) {
  if (ncell > 1L) {
    centers = dplyr::slice_sample(tbl[c("x", "y", "z")], n = nsam)
    sample_regions(tbl, centers, ncell = ncell)
  } else {
    tbl[c("x", "y", "z", "id")] |>
      dplyr::slice_sample(n = nsam) |>
      dplyr::mutate(id = as.list(.data$id))
  }
}

#' @details
#' `sample_bulk` samples a bulk of cells near the specified center.
#' @param center named (x, y, z) vector, list, or tibble
#' @rdname sample
#' @export
sample_bulk = function(tbl, center = c(x = 0, y = 0, z = 0), ncell = 10L) {
  d = dist_euclidean(tbl, center)
  tbl$id[d < nth_element(d, ncell)]
}

kmeans_centers = function(tbl, centers, iter.max = 32L) {
  tbl = tbl[c("x", "y", "z")]
  result = stats::kmeans(tbl, centers = centers, iter.max = iter.max)
  tibble::as_tibble(result[["centers"]])
}

sample_regions = function(tbl, centers, ncell = 10L) {
  covr_issue377 = function(x, y, z) {
    sample_bulk(tbl, center = c(x = x, y = y, z = z), ncell = ncell)
  }
  centers$id = purrr::pmap(centers, covr_issue377)
  centers
}

tidy_regions = function(regions) {
  tibble::tibble(id = regions$id) |>
    tibble::rowid_to_column(var = "region") |>
    tidyr::unnest("id")
}

#' @details
#' `evaluate_mrs` is a shortcut to evaluate multi-region sampling.
#' @param population tbl
#' @inheritParams filter_common_ancestors
#' @inheritParams internal_nodes
#' @rdname sample
#' @export
evaluate_mrs = function(population, nsam, ncell, threshold = 0.05, sensitivity = 0.05, jitter = 0) {
  graph = make_igraph(population)
  regions = filter_extant(population) |>
    sample_uniform_regions(nsam = nsam, ncell = ncell, jitter = jitter)
  sampled = purrr::flatten_int(regions$id)
  subgraph = subtree(graph, sampled)
  detectable = internal_nodes(subgraph, sampled, sensitivity = sensitivity)
  major_ca = population |>
    add_node_property(graph) |>
    filter_common_ancestors(threshold = threshold) |>
    dplyr::pull("id")
  sum(detectable %in% major_ca) / length(major_ca)
}

#' @details
#' `distances_mrs` is a shortcut for sampling and calculation
#' @rdname sample
#' @export
distances_mrs = function(population, nsam, ncell, jitter = 0) {
  graph = make_igraph(population)
  regions = filter_extant(population) |>
    sample_uniform_regions(nsam = nsam, ncell = ncell, jitter = jitter)
  sampled = purrr::flatten_int(regions$id)
  subgraph = subtree(graph, sampled)
  pairwise_distances(subgraph, regions)
}
