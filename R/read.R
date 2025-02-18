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
    paste(collapse = "\n") |>
    readr::read_tsv()
}

#' @details
#' `read_results` reads configs and populations as a nested tibble.
#' @param graph add graph column if TRUE
#' @rdname read
#' @export
read_results = function(
    indirs = getwd(),
    graph = getOption("tumopp.graph", TRUE),
    mc.cores = getOption("mc.cores", 1L)) {
  parallel::mclapply(indirs, .read_result, mc.cores = mc.cores, graph = graph) |> purrr::list_rbind()
}

.read_result = function(indir, graph) {
  res = .read_conf(indir)
  population = read_tumopp(file.path(indir, "population.tsv.gz"))
  transforming = (res$coord == "hex") && getOption("tumopp.autohex", TRUE)
  if (transforming) {
    population = trans_coord_hex(population)
  }
  res[["population"]] = list(population)
  if (graph) {
    res[["graph"]] = list(make_igraph(population))
  }
  snapshots = file.path(indir, "snapshots.tsv.gz")
  if (file.exists(snapshots)) {
    snapshots = read_tumopp(snapshots)
    if (transforming) {
      snapshots = trans_coord_hex(snapshots)
    }
    res[["snapshots"]] = list(snapshots)
  }
  drivers = file.path(indir, "drivers.tsv.gz")
  if (file.exists(drivers)) {
    res[["drivers"]] = list(readr::read_tsv(drivers))
  }
  benchmark = file.path(indir, "benchmark.tsv.gz")
  if (file.exists(benchmark)) {
    res[["benchmark"]] = list(readr::read_tsv(benchmark))
  }
  res
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
