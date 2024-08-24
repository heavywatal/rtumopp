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
  centers = kmeans_xyz(tbl, nsam)
  if (jitter > 0) {
    xy = centers[, c("x", "y")]
    centers[, c("x", "y")] = xy + stats::runif(length(xy), -jitter, jitter)
  }
  centers = as.data.frame(centers) |> tibble::new_tibble()
  centers$id = purrr::pmap(centers, \(...) {
    sample_bulk(tbl, c(...), ncell = ncell)
  })
  centers
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
  tbl$id[d <= nth_element(d, ncell)][seq_len(ncell)]
}

kmeans_xyz = function(tbl, centers, iter.max = 2L) {
  tbl = tbl[c("x", "y", "z")]
  if (requireNamespace("ClusterR", quietly = TRUE)) {
    kmeans_arma(tbl, centers, iter.max = iter.max)
  } else {
    kmeans_base(tbl, centers, iter.max = iter.max)[["centers"]]
  }
}

kmeans_arma = function(tbl, centers, iter.max = 2L) {
  # ClusterR calls set.seed() and the default seed is 1.
  res = ClusterR::KMeans_arma(tbl, centers, n_iter = iter.max, seed = runif.int(1L))
  class(res) = NULL
  colnames(res) = colnames(tbl)
  res
}

kmeans_base = function(tbl, centers, iter.max = 2L) {
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

#' @rdname sample
#' @export
sim_sample_single_cell = function(graph, nrep = 1L, ncell = 100L, mu = 0, accel = 0) {
  on.exit(stats::runif(1L)) # proxy of parallel::nextRNGStream(.Random.seed)
  vsink = igraphlite::Vsink(graph)  # slower than you think
  parallel::mclapply(seq_len(nrep), \(i) {
    sampling_single_cell(graph, vsink, ncell, mu, accel)
  }) |>
    purrr::list_rbind(names_to = "repl") |>
    mutate_chisq(mu = mu, accel = accel)
}

sampling_single_cell = function(graph, vsink, ncell, mu, accel = 0) {
  vsam = sample(vsink, ncell, replace = FALSE)
  graph = subgraph_upstream(graph, vsam, trim = FALSE)
  vsam = igraphlite::Vsink(graph)
  lens = edge_lengths(graph, mu = mu, accel = accel)
  dista = distances_upstream(graph, vsam, weights = lens, trim = TRUE)
  summary_row(dista)
}

#' @param nrep number of replications.
#' @inheritParams edge_lengths
#' @rdname sample
#' @export
sim_sample_biopsy = function(graph, population, nrep = 1L, nsam = 5L, ncell = 100L, mu = 0, accel = 0) {
  on.exit(stats::runif(1L)) # proxy of parallel::nextRNGStream(.Random.seed)
  extant = population |> tumopp::filter_extant()
  parallel::mclapply(seq_len(nrep), \(i) {
    sampling_biopsy(graph, extant, nsam = nsam, ncell = ncell, mu = mu, accel = accel)
  }) |>
    purrr::list_rbind(names_to = "repl") |>
    mutate_chisq(mu = mu, accel = accel)
}

sampling_biopsy = function(graph, extant, nsam, ncell, mu, accel = 0) {
  regions = sample_uniform_regions(extant, nsam = nsam, ncell = ncell)
  subgraph = subtree(graph, unlist(regions$id, use.names = FALSE), trim = TRUE)
  if (accel > 0) warning("accel > 0 on a trimmed tree")
  vaf = make_vaf(subgraph, regions$id, mu = mu, accel = accel)
  dista = colSums(vaf)
  summary_row(dista)
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
