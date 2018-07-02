#' Functions to modify population data.frame
#'
#' `modify_population` add various columns to the raw population data
#' @param population tibble
#' @param coord string
#' @param dimensions integer
#' @param ... ignored
#' @param num_clades integer
#' @rdname population
modify_population = function(population, coord, dimensions, ..., num_clades = 4L) {
  extant = filter_extant(population)
  strelem = get_se(coord, dimensions)
  col_surface = detect_surface(extant, strelem) %>%
    dplyr::select(.data$id, .data$surface)
  if (coord == "hex") {
    population = trans_coord_hex(population)
  }
  max_phi = c(hex = 12L, moore = 27L, neumann = 6L)[coord]
  population %>%
    set_graph_property(num_clades = num_clades) %>%
    dplyr::mutate(r = dist_euclidean(.), phi = .data$phi / max_phi) %>%
    dplyr::left_join(col_surface, by = "id")
}

#' `filter_extant` collects the extant cells at the end of the simulation
#' @rdname population
#' @export
filter_extant = function(population) {
  dplyr::filter(population, .data$death == 0)
}

#' `filter_common_ancestors` collects major common ancestors
#' @param threshold minimum frequency of detectable alleles
#' @rdname population
#' @export
filter_common_ancestors = function(population, threshold = 0.05) {
  dplyr::filter(population, .data$allelefreq >= threshold)
}

# Add age and clade column
set_graph_property = function(population, num_clades) {
  .graph = make_igraph(population)
  .nodes = as.character(population$id)
  .size = count_sink(.graph)
  .out = dplyr::mutate(
    population,
    age = distances_from_origin(.graph, .nodes),
    allelefreq = count_sink(.graph, .nodes) / .size
  )
  founders = list_clade_founders(.out, num_clades = num_clades)
  clade_data = founders %>%
    as.character() %>%
    rlang::set_names() %>%
    purrr::map_dfr(~{
      tibble::tibble(id = igraph::subcomponent(.graph, .x, mode = "out")$name)
    }, .id = "clade") %>%
    dplyr::mutate(
      id = as.integer(.data$id),
      clade = factor(.data$clade, levels = as.character(founders))
    )
  dplyr::left_join(.out, clade_data, by = "id")
}
