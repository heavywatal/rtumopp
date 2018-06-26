#' Calculate MATH (mutant-allele tumor heterogeneity) score from variant allele frequencies
#' @inheritParams stats::mad
#' @return numeric
#' @rdname frequency
#' @export
math_score = function(x, constant=1.4826, na.rm=FALSE) {
  med = stats::median(x, na.rm = na.rm)
  mad = stats::mad(x, center = med, constant = constant, na.rm = na.rm)
  mad / med
}

#' Extract cells whose expected allele frequencies are above threshold
#' @param population tibble with descendants and id column
#' @param threshold lowerbound of detectable allele frequency
#' @return integer IDs
#' @rdname frequency
#' @export
detectable_mutants_all = function(population, threshold) {
  population$id[population$allelefreq > threshold]
}

#' @param graph igraph
#' @param nodes igraph vertices
#' @rdname frequency
#' @export
detectable_mutants = function(graph, nodes, threshold) {
  n = length(nodes)
  counts = paths_to_source(graph, nodes) %>%
    purrr::flatten_chr() %>%
    table()
  counts[(counts / n) > threshold] %>% names() %>% as.integer()
}
