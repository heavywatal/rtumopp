#' Sprinkle mutations on genealogy
#'
#' `mutate_clades` samples nodes to mutate and returns their descendants
#' @param graph igraph
#' @param mu mutation rate per cell division (ignored if segsites is given)
#' @param segsites number of segregating sites
#' @rdname ms
#' @export
mutate_clades = function(graph, mu=NULL, segsites=NULL) {
  nodes = names(igraph::V(graph))
  nodes = nodes[igraph::degree(graph, mode = "out") > 0L]
  # remove singletons
  # TODO: remove low-freq variants?
  if (is.null(segsites)) {
    if (is.null(mu)) stop("specify either mu or segsites")
    if (mu > 0) {
      segsites = stats::rpois(1L, length(nodes) * mu)
    }
  } else if (!is.null(mu)) warning("mu is ignored if segsites is given")
  mutants = if (is.null(segsites)) {
    nodes  # if (mu <= 0)
  } else {
    sample(nodes, segsites, replace = TRUE)
  }
  # TODO: remove internal nodes?
  paths_to_sink(graph, mutants)
}
