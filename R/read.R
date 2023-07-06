#' Read TSV files in given directories
#'
#' `read_confs` reads config files.
#' @param indirs a string vector
#' @param mc.cores The number of cores to use for concurrent execution.
#' @rdname read
#' @export
read_confs = function(indirs = getwd(), mc.cores = getOption("mc.cores", 1L)) {
  stats::setNames(indirs, indirs) |>
    parallel::mclapply(.read_conf, mc.cores = mc.cores) |>
    dplyr::bind_rows(.id = "directory")
}

.read_conf = function(indir) {
  json = file.path(indir, "config.json")
  if (file.exists(json)) {
    return(read_json(json))
  }
  tsv = file.path(indir, "program_options.tsv.gz")
  if (file.exists(tsv)) {
    return(readr::read_tsv(tsv))
  }
  read_boost_ini(file.path(indir, "program_options.conf"))
}

read_json = function(file) {
  jsonlite::read_json(file) |> tibble::as_tibble()
}

from_json = function(conf) {
  jsonlite::fromJSON(conf) |> tibble::as_tibble()
}

read_boost_ini = function(file) {
  readr::read_delim(file, "=", col_names = c("key", "val"), comment = "#", trim_ws = TRUE) |>
    dplyr::summarize(dplyr::across(dplyr::everything(), paste0, collapse = "\t")) |>
    paste0(collapse = "\n") |>
    readr::read_tsv()
}

#' @details
#' `read_results` reads confs and populations as a nested tibble.
#' @param graph add graph column if TRUE
#' @rdname read
#' @export
read_results = function(indirs = getwd(), mc.cores = getOption("mc.cores", 1L), graph = TRUE) {
  autohex = getOption("tumopp.autohex", TRUE)
  pop = read_populations(indirs, mc.cores = mc.cores)
  results = read_confs(indirs) |> dplyr::mutate(
    population = purrr::map_if(pop, (.data$coord == "hex") & autohex, trans_coord_hex)
  )
  if (graph) {
    results = dplyr::mutate(results, graph = parallel::mclapply(results$population, make_igraph, mc.cores = mc.cores))
  }
  results
}

#' @details
#' `read_populations` reads populations.
#' @rdname read
#' @export
read_populations = function(indirs = getwd(), mc.cores = getOption("mc.cores", 1L)) {
  file.path(indirs, "population.tsv.gz") |>
    parallel::mclapply(read_tumopp, mc.cores = mc.cores)
}

#' @details
#' `read_snapshots` calls [read_results()] and reads snapshots.
#' @rdname read
#' @export
read_snapshots = function(indirs = getwd(), mc.cores = getOption("mc.cores", 1L)) {
  read_results(indirs, mc.cores = mc.cores) |>
    dplyr::mutate(snapshots = purrr::map2(indirs, .data$population, \(x, y) {
      meta_info = dplyr::select(y, "id", "clade")
      read_tumopp(file.path(x, "snapshots.tsv.gz")) |>
        dplyr::left_join(meta_info, by = "id")
    }))
}

#' @details
#' `read_tumopp` is an alias of [readr::read_tsv()] with some options
#' to read result files: `population.tsv.gz` and `snapshots.tsv.gz`.
#' @param file passed as the first argument of [readr::read_tsv()]
#' @rdname read
#' @export
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
