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
  if (is.null(segsites)) {
    if (is.null(mu)) stop("specify either mu or segsites")
    if (mu > 0) {
      segsites = stats::rpois(1L, length(nodes) * mu)
    }
  } else if (!is.null(mu)) warning("mu is ignored if segsites is given")
  mutants = if (is.null(segsites)) {
    nodes # if (mu <= 0)
  } else {
    sample(nodes, segsites, replace = TRUE)
  }
  # TODO: remove internal nodes?
  paths_to_sink(graph, mutants)
}

#' @description
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
  }
  subgraph = subtree(graph, nodes)
  segsites = mutate_clades(subgraph, mu = mu, segsites = segsites)
  cols = purrr::map(segsites, ~ as.integer(nodes %in% .x))
  m = do.call(cbind, cols)
  rownames(m) = nodes
  m
}
