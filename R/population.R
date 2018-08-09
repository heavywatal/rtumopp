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
#' `filter_common_ancestors` collects major common ancestors.
#' @param threshold minimum frequency of detectable alleles
#' @rdname population
#' @export
filter_common_ancestors = function(population, threshold = 0.05) {
  dplyr::filter(population, .data$allelefreq >= threshold)
}

# Add graph-related columns
add_node_property = function(population, graph, num_clades = 4L) {
  .order = length(population$death == 0)
  clade_data = sort_clades(population, graph, n = num_clades)
  population %>%
    dplyr::mutate(
      age = distances_from_origin(graph, .data$id),
      allelefreq = count_sink(graph, .data$id) / .order
    ) %>%
    dplyr::left_join(clade_data, by = "id")
}

sort_clades = function(population, graph, n = 4L) {
  founders = list_clade_founders(population, n)
  tibble::tibble(
    clade = factor(founders),
    id = paths_to_sink(graph, founders)
  ) %>% tidyr::unnest()
}
