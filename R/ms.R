#' Sprinkle mutations on genealogy
#'
#' `mutate_clades` samples nodes to mutate and returns their descendants.
#' @param graph igraph
#' @param mu mutation rate per cell division (ignored if segsites is given)
#' @param segsites number of segregating sites
#' @rdname ms
#' @export
mutate_clades = function(graph, mu = NULL, segsites = NULL) {
  vs = igraph::V(graph)
  indegree = igraph::degree(graph, vs, mode = "in", loops = FALSE)
  nodes = as_ids(vs)[indegree > 0L]
  # TODO: remove low-freq variants?
  number = if (is.null(mu)) {
    if (is.null(segsites)) stop("specify either mu or segsites")
    stats::rmultinom(1L, segsites, rep(1, length(nodes))) %>% as.vector()
  } else if (mu > 0) {
    if (!is.null(segsites)) stop("segsites is ignored if mu is given")
    stats::rpois(length(nodes), mu)
  } else {
    rep(1L, length(nodes))
  }
  tibble::tibble(origin = nodes, number = number) %>%
    dplyr::filter(.data$number > 0) %>%
    dplyr::mutate(carriers = paths_to_sink(graph, .data$origin))
}

#' @details
#' `make_sample` creates a genotype matrix from a given genealogy tree.
#' @param nsam number of cells to sample
#' @rdname ms
#' @export
make_sample = function(graph, nsam = 0L, mu = NULL, segsites = NULL) {
  vs = igraph::V(graph)
  outdegree = igraph::degree(graph, vs, mode = "out", loops = FALSE)
  nodes = as_ids(vs)[outdegree == 0L]
  if (nsam > 0L) {
    nodes = sample(nodes, nsam, replace = FALSE)
    graph = subtree(graph, nodes)
  }
  mutations = mutate_clades(graph, mu = mu, segsites = segsites)
  origins = purrr::map(mutations$carriers, ~ as.integer(nodes %in% .x))
  m = t(do.call(rbind, origins) * mutations$number)
  rownames(m) = nodes
  colnames(m) = mutations$origin
  m
}
