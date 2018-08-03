#' Read TSV files in given directories
#'
#' `read_confs` reads config files.
#' @param indirs a string vector
#' @rdname read
#' @export
read_confs = function(indirs = getwd()) {
  file.path(indirs, "program_options.conf") %>%
    stats::setNames(indirs) %>%
    purrr::map_dfr(wtl::read_boost_ini, .id = "directory")
}

#' @description
#' `read_populations` reads populations.
#' @rdname read
#' @export
read_populations = function(indirs = getwd()) {
  file.path(indirs, "population.tsv.gz") %>%
    purrr::map(readr::read_tsv, col_types = readr::cols(
      beta = readr::col_double(),
      delta = readr::col_double(),
      rho = readr::col_double()
    ))
}

#' @description
#' `read_results` reads confs and populations as a nested tibble,
#' @rdname read
#' @export
read_results = function(indirs = getwd()) {
  autohex = getOption("tumopp.autohex", TRUE)
  read_confs(indirs) %>% dplyr::mutate(
    population = read_populations(indirs) %>%
      purrr::map_if((.data$coord == "hex") & autohex, trans_coord_hex),
    graph = purrr::map(.data$population, make_igraph)
  )
}

#' @description
#' `read_snapshots` calls `read_results` and reads snapshots.
#' @rdname read
#' @export
read_snapshots = function(indirs = getwd()) {
  read_results(indirs) %>%
    dplyr::mutate(snapshots = purrr::map2(.data$directory, .data$population, ~{
      meta_info = dplyr::select(.y, .data$id, .data$x, .data$y, .data$z, .data$clade)
      readr::read_tsv(file.path(.x, "snapshots.tsv.gz")) %>%
        dplyr::select(-.data$x, -.data$y, -.data$z) %>%
        dplyr::left_join(meta_info, by = "id")
    }))
}
