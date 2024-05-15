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
  tbl = tbl[c("x", "y", "z", "id")]
  km = kmeans_xyz(tbl, nsam)
  centers = km[["centers"]]
  if (jitter > 0) {
    xy = centers[, c("x", "y")]
    centers[, c("x", "y")] = xy + stats::runif(length(xy), -jitter, jitter)
  }
  id = tbl |>
    dplyr::mutate(cluster = km[["cluster"]]) |>
    tidyr::nest(.by = "cluster") |>
    dplyr::arrange("cluster") |>
    purrr::pmap(\(cluster, data) {
      sample_bulk(data, centers[cluster, ], ncell = ncell)
    })
  as.data.frame(centers) |>
    tibble::new_tibble() |>
    dplyr::mutate(id = id)
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
  stopifnot(ncell < nrow(tbl))
  d = dist_euclidean(tbl, center)
  tbl$id[d < nth_element(d, ncell)]
}

kmeans_xyz = function(tbl, centers, iter.max = 2L) {
  tbl = tbl[c("x", "y", "z")]
  suppressWarnings(stats::kmeans(tbl, centers, iter.max = iter.max))
}

sample_regions = function(tbl, centers, ncell = 10L) {
  centers$id = purrr::pmap(centers, \(x, y, z) {
    sample_bulk(tbl, center = c(x = x, y = y, z = z), ncell = ncell)
  })
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
  sampled = purrr::list_c(regions$id)
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
  sampled = purrr::list_c(regions$id)
  subgraph = subtree(graph, sampled)
  pairwise_distances(subgraph, regions)
}
