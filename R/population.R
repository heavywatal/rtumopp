#' Functions to modify population data.frame
#'
#' `filter_extant` collects the extant cells at the end of the simulation.
#' @param population tibble
#' @rdname population
#' @export
filter_extant = function(population) {
  dplyr::filter(population, .data$death == 0)
}

#' @description
#' `count_extant` counts the number of extant cells
#' @rdname population
#' @export
count_extant = function(population) {
  sum(population$death == 0)
}

#' @description
#' `filter_common_ancestors` collects major common ancestors.
#' @param threshold minimum frequency of detectable alleles
#' @rdname population
#' @export
filter_common_ancestors = function(population, threshold = 0.05) {
  dplyr::filter(population, .data$allelefreq >= threshold)
}

# Add graph-related columns
add_node_property = function(population, graph, num_clades = 4L) {
  extant = filter_extant(population)
  clade_data = sort_clades(population, graph, n = num_clades)
  freq_data = count_extant_descendants(graph, extant$id) %>%
    dplyr::transmute(.data$id, allelefreq = .data$n / nrow(extant))
  population %>%
    dplyr::mutate(age = distances_from_origin(graph, .data$id)) %>%
    dplyr::left_join(freq_data, by = "id") %>%
    dplyr::left_join(clade_data, by = "id")
}

sort_clades = function(population, graph, n = 4L) {
  founders = list_clade_founders(population, n)
  tibble::tibble(
    clade = factor(founders),
    id = paths_to_sink(graph, founders)
  ) %>% tidyr::unnest()
}

count_extant_descendants = function(graph, nodes) {
  paths = paths_to_source(graph, nodes)
  tibble::tibble(id = purrr::flatten_int(paths)) %>%
    dplyr::count(.data$id)
}
