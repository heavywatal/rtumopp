#' Run tumopp
#'
#' `tumopp()` returns full results with config columns in a data.frame.
#' See `tumopp("-h")` or <https://heavywatal.github.io/tumopp/group__params.html>
#' for the list of command options.
#'
#' A `population` data.frame includes ancestral cells.
#' Extant cells can be extracted by filtering with `death == 0` or [filter_extant()] function.
#' The sampling time, i.e., the end of a simulation, is typically the maximum value of birth time.
#'
#' The default unit of time (`birth` and `death` columns) is the average cell cycle of newborn cells
#' (given the parameter `-b`/`--beta0` is set to 1).
#' For example, step-wise tumor growth and integer values in the birth column will be observed
#' if `-k`/`--shape` parameter is set to a very large value like `10**6`.
#' If you are considering some cell line whose average cell cycle is 4 days for example,
#' then the unit of those columns can be interpreted as 4 days,
#' or you can set `--beta0=0.25` to change the unit to a day.
#'
#' The `omega` column denotes the number of cell divisions allowed for each cell.
#' Negative values denote unlimited proliferation potential.
#' @param args command line arguments as a string vector or list of strings.
#' @param ... not used.
#' @export
tumopp = function(args, ...) UseMethod("tumopp")

#' @param graph add graph column if TRUE
#' @param cache A parent directory to cache results.
#' `TRUE` is equivalent to `getOption("tumopp.cache", "~/.cache/tumopp")`.
#' [tempdir()] is used in other cases including the default (`NULL`),
#' which is discarded at the end of an R session.
#' Set `FALSE` to force `tumopp` to run and overwrite previous results if any.
#' @rdname tumopp
#' @export
tumopp.default = function(
    args = character(0L), ...,
    graph = getOption("tumopp.graph", TRUE),
    cache = NULL) {
  if (length(args) == 1L) {
    args = stringr::str_split_1(args, "\\s+")
  }
  if ("--version" %in% args) {
    return(system2(tumopp_path(), "--version", stdout = TRUE))
  }
  if (any(c("-h", "--help") %in% args)) {
    msg = system2(tumopp_path(), args, stdout = TRUE)
    return(invisible(message(paste(msg, collapse = "\n"))))
  }
  args = append_seed(args)
  cache_dir = sanitize_cache_dir(cache)
  outdir = file.path(cache_dir, cache_name(args))
  if (isFALSE(cache) || !dir.exists(outdir)) {
    ret = system2(tumopp_path(), c(args, "-o", outdir))
    stopifnot(ret == 0L)
  }
  .read_result(outdir, graph = graph)
}

sanitize_cache_dir = function(cache) {
  if (isTRUE(cache)) {
    cache = getOption("tumopp.cache", "~/.cache/tumopp")
    dir.create(cache, showWarnings = FALSE, recursive = TRUE, mode = "0755")
  } else if (length(cache) != 1L || !is.character(cache) || is.na(cache) || !nzchar(cache)) {
    cache = tempdir()
  }
  cache
}

cache_name = function(args) {
  x = stringr::str_flatten(args) |>
    stringr::str_remove_all("[ =.]+") |>
    stringr::str_replace_all("-+", "-")
  paste0(".tumopp", x)
}

#' @param mc.cores The number of cores to use for concurrent execution.
#' @rdname tumopp
#' @export
tumopp.list = function(args, ..., graph = TRUE, mc.cores = getOption("mc.cores", 1L)) {
  args = purrr::map(args, append_seed)
  out = parallel::mclapply(args, tumopp, ..., graph = FALSE, mc.cores = mc.cores) |>
    dplyr::bind_rows(.id = "args")
  if (graph) {
    out$graph = lapply(out[["population"]], make_igraph)
  }
  out
}

#' @rdname tumopp
#' @export
tumopp.data.frame = function(args, ..., graph = TRUE, mc.cores = getOption("mc.cores", 1L)) {
  tumopp(vectorize_args(args), ..., graph = graph, mc.cores = mc.cores)
}

#' @details
#' `make_args()` returns argument combinations in a tibble.
#' @param alt named list of altered arguments.
#' @param const named list of constant arguments.
#' @param times,each passed to [rep()]
#' @rdname tumopp
#' @export
make_args = function(alt, const = NULL, times = 1L, each = 1L) {
  .grid = rlang::exec(tidyr::crossing, !!!alt) |> filter_valid_LP()
  idx = rep(seq_len(nrow(.grid)), times = times, each = each)
  dplyr::slice(.grid, idx) |>
    dplyr::mutate(!!!const)
}

vectorize_args = function(.tbl) {
  # suppress scientific notation like "-N1e+06"
  withr::with_options(list(scipen = 255L), {
    purrr::pmap(.tbl, function(...) {
      .params = list(...) |>
        purrr::discard(isFALSE) |>
        unlist()
      .names = names(.params)
      .template = ifelse(nchar(.names) > 1L, "--%s=%s", "-%s%s")
      sprintf(.template, .names, .params) |>
        stringr::str_remove("=?TRUE$")
    })
  })
}

# prior is a named list of generating functions
generate_args = function(prior, const = NULL, n = 1L) {
  generate_valid(prior, n = n) |>
    vectorize_args() |>
    purrr::map(\(x) c(const, x))
}

generate_valid = function(prior, n = 1L) {
  output = NULL
  while (n > 0L) {
    generated = purrr::map(prior, rlang::exec, n) |>
      tibble::as_tibble()
    generated = suppressMessages(filter_valid_LP(generated))
    output = dplyr::bind_rows(output, generated)
    n = n - nrow(output)
  }
  output
}

filter_valid_LP = function(x) {
  if (all(c("L", "P") %in% names(x))) {
    nrow_orig = nrow(x)
    x = dplyr::inner_join(x, valid_LP_combinations(), by = c("L", "P"))
    num_excluded = nrow_orig - nrow(x)
    if (num_excluded > 0L) {
      message("Note: ", num_excluded, " L-P combinations were invalid and excluded")
    }
  }
  x
}

valid_LP_combinations = function() {
  rbind(
    tidyr::crossing(L = c("const", "linear", "step"), P = c("random", "mindrag")),
    tidyr::crossing(L = "const", P = c("roulette", "minstraight", "stroll"))
  )
}

append_seed = function(args, seed = NULL) {
  if (!any(stringr::str_starts(args, "--seed"))) {
    args = c(args, paste0("--seed=", seed %||% runif.int(1L)))
  }
  args
}

runif.int = function(n, min = -.Machine$integer.max, max = .Machine$integer.max) {
  .offset = min - 1
  as.integer(sample.int(max - .offset, n, replace = TRUE) + .offset)
}
