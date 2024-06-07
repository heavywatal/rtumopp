#' Sprinkle mutations on genealogy
#'
#' @details
#' `make_sample` creates a genotype matrix from a given genealogy tree.
#' @param graph igraph
#' @param nsam number of cells to sample
#' @param mu mutation rate per cell division (ignored if segsites is given)
#' @param segsites number of segregating sites
#' @rdname ms
#' @export
make_sample = function(graph, nsam = 0L, mu = NULL, segsites = NULL) {
  nodes = igraphlite::Vnames(graph)[igraphlite::is_sink(graph)]
  if (nsam > 0L) {
    nodes = sample(nodes, nsam, replace = FALSE)
    graph = subtree(graph, nodes)
  }
  mutations = mutate_clades(graph, mu = mu, segsites = segsites)
  origins = purrr::map(mutations$carriers, \(x) as.integer(nodes %in% x))
  m = t(do.call(rbind, origins) * mutations$number)
  rownames(m) = nodes
  colnames(m) = mutations$origin
  m
}

mutate_clades = function(graph, mu = NULL, segsites = NULL) {
  nodes = igraphlite::Vnames(graph)[!igraphlite::is_source(graph)]
  number = if (is.null(mu)) {
    if (is.null(segsites)) stop("specify either mu or segsites")
    stats::rmultinom(1L, segsites, rep_len(1L, length(nodes))) |> as.vector()
  } else if (mu > 0) {
    if (!is.null(segsites)) stop("segsites is ignored if mu is given")
    stats::rpois(length(nodes), mu)
  } else {
    rep_len(1L, length(nodes))
  }
  tibble::tibble(origin = nodes, number = number) |>
    dplyr::filter(.data$number > 0) |>
    dplyr::mutate(carriers = paths_to_sink(graph, .data$origin))
}

#' @details
#' `edge_lengths()` calculates lengths of edges on a given genealogy tree.
#' @param accel assumption undocumented yet
#' @rdname ms
#' @export
edge_lengths = function(graph, mu = 0, accel = 0) {
  # distances(): "unweighted" is slower than "dijkstra"
  ecount = igraphlite::ecount(graph)
  lens = rep_len(1L, ecount)
  if (accel > 0) {
    lens = (1 + accel)**edge_ages(graph)
  }
  if (mu > 0) {
    lens = stats::rpois(ecount, mu * lens)
  }
  lens
}

edge_ages = function(graph) {
  age = igraphlite::Eattr(graph, "age")
  if (is.null(age)) {
    age = distances_upstream(graph, igraphlite::igraph_to(graph), trim = FALSE)
    igraphlite::Eattr(graph, "age") = age
  }
  age
}

mutate_edges = function(graph, mu = 0, accel = 0) {
  igraphlite::Eattr(graph, "weight") = edge_lengths(graph, mu = mu, accel = accel)
  graph
}
