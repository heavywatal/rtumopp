#' Read TSV files in given directories
#'
#' `read_confs` reads config files.
#' @param indirs a string vector
#' @rdname read
#' @export
read_confs = function(indirs = getwd()) {
  stats::setNames(indirs, indirs) %>%
    parallel::mclapply(.read_conf) %>%
    dplyr::bind_rows(.id = "directory")
}

.read_conf = function(indir) {
  tsv = file.path(indir, "program_options.tsv.gz")
  if (file.exists(tsv)) {
    readr::read_tsv(tsv)
  } else {
    wtl::read_boost_ini(file.path(indir, "program_options.conf"))
  }
}

#' @description
#' `read_populations` reads populations.
#' @rdname read
#' @export
read_populations = function(indirs = getwd()) {
  file.path(indirs, "population.tsv.gz") %>%
    parallel::mclapply(read_tumopp)
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
    graph = parallel::mclapply(.data$population, make_igraph)
  )
}

#' @description
#' `read_snapshots` calls `read_results` and reads snapshots.
#' @rdname read
#' @export
read_snapshots = function(indirs = getwd()) {
  read_results(indirs) %>%
    dplyr::mutate(snapshots = purrr::map2(indirs, .data$population, ~{
      meta_info = dplyr::select(.y, .data$id, .data$clade)
      read_tumopp(file.path(.x, "snapshots.tsv.gz")) %>%
        dplyr::left_join(meta_info, by = "id")
    }))
}

read_tumopp = function(file) {
  readr::read_tsv(file, col_types = readr::cols(
    x = readr::col_integer(),
    y = readr::col_integer(),
    z = readr::col_integer(),
    id = readr::col_integer(),
    ancestor = readr::col_integer(),
    birth = readr::col_double(),
    death = readr::col_double(),
    omega = readr::col_integer()
  ), progress = FALSE)
}
