#' Modify population table
#' @param population tibble
#' @param coord string
#' @param dimensions integer
#' @param ... ignored
#' @param num_clades integer
#' @rdname population
modify_population = function(population, coord, dimensions, ..., num_clades=4L) {
  extant = filter_extant(population)
  strelem = get_se(coord, dimensions)
  col_surface = detect_surface(extant, strelem) %>%
    dplyr::select(.data$id, .data$surface)
  if (coord == "hex") {
    population = trans_coord_hex(population)
  }
  max_phi = c(hex = 12L, moore = 27L, neumann = 6L)[coord]
  population %>%
    set_graph_property() %>%
    dplyr::mutate(r = dist_euclidean(.), phi = .data$phi / max_phi) %>%
    dplyr::left_join(col_surface, by = "id")
}

#' Filter extant cells
#' @rdname population
#' @export
filter_extant = function(population) {
  dplyr::filter(population, .data$death == 0)
}

# Add age and clade column
set_graph_property = function(population) {
  .graph = make_igraph(population)
  .nodes = as.character(population$id)
  .out = dplyr::mutate(
    population,
    age = distances_from_origin(.graph, .nodes),
    allelefreq = allele_freqs(.graph, .nodes)
  )
  founders = list_clade_founders(.out, 4L)
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

list_clade_founders = function(population, num_clades) {
  origin = sum(population$age == 0L)
  stopifnot(num_clades >= origin)
  num_divisions = num_clades - origin
  roots = utils::head(population$id, num_divisions)
  seq_len(num_divisions + num_clades) %>% setdiff(roots)
}
