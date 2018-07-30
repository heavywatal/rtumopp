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

# Add various columns to the raw population data
modify_population = function(population, graph, coord, dimensions, ..., num_clades = 4L) {
  extant = filter_extant(population)
  strelem = get_se(coord, dimensions)
  col_surface = detect_surface(extant, strelem) %>%
    dplyr::select(.data$id, .data$surface)
  if (coord == "hex") {
    population = trans_coord_hex(population)
  }
  max_phi = c(hex = 12L, moore = 27L, neumann = 6L)[coord]
  population %>%
    add_node_property(graph, num_clades) %>%
    dplyr::mutate(r = dist_euclidean(.), phi = .data$phi / max_phi) %>%
    dplyr::left_join(col_surface, by = "id")
}

# Add graph-related columns
add_node_property = function(population, graph, num_clades) {
  .nodes = as.character(population$id)
  .size = count_sink(graph)
  founders = list_clade_founders(population, num_clades = num_clades)
  clade_data = tibble::tibble(
    clade = factor(founders),
    id = paths_to_sink(graph, as.character(founders)) %>% purrr::map(as.integer)
  ) %>% tidyr::unnest()
  population %>%
    dplyr::mutate(
      age = distances_from_origin(graph, .nodes),
      allelefreq = count_sink(graph, .nodes) / .size
    ) %>%
    dplyr::left_join(clade_data, by = "id")
}
