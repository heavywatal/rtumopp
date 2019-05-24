#' Run tumopp
#'
#' `tumopp()` returns full results with config columns in a data.frame.
#' @param args command line arguments as a string vector or list of strings.
#' @param ... not used.
#' @export
tumopp = function(args, ...) UseMethod("tumopp")

#' @param graph add graph column if TRUE
#' @rdname tumopp
#' @export
tumopp.default = function(args = character(0L), ..., graph = TRUE) {
  if (length(args) == 1L) {
    args = stringr::str_split(args, "\\s+") %>% purrr::flatten_chr()
  }
  result = cpp_tumopp(args)
  if (length(result) == 0L) return(invisible(NULL))
  .out = from_json(result["config"])
  .pop = read_tumopp(result["population"])
  transforming = ((.out$coord == "hex") && getOption("tumopp.autohex", TRUE))
  if (transforming) {
    .pop = trans_coord_hex(.pop)
  }
  .out$population = list(.pop)
  if (graph) {
    .out$graph = list(make_igraph(.pop))
  }
  .snapshots = result["snapshots"]
  if (nzchar(.snapshots)) {
    .snapshots = read_tumopp(.snapshots)
    if (transforming) {
      .snapshots = trans_coord_hex(.snapshots)
    }
    .out$snapshots = list(.snapshots)
  }
  .drivers = result["drivers"]
  if (nzchar(.drivers)) {
    .out$drivers = list(readr::read_tsv(.drivers))
  }
  .benchmark = result["benchmark"]
  if (nzchar(.benchmark)) {
    .out$benchmark = list(readr::read_tsv(.benchmark))
  }
  .out
}

#' @param mc.cores The number of cores to use for concurrent execution.
#' @rdname tumopp
#' @export
tumopp.list = function(args, ..., graph = TRUE, mc.cores = getOption("mc.cores", 1L)) {
  out = parallel::mclapply(args, tumopp, ..., graph = FALSE, mc.cores = mc.cores) %>%
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
#' @param const named vector of constant arguments.
#' @param times,each passed to `rep()`
#' @rdname tumopp
#' @export
make_args = function(alt, const = NULL, times = 1L, each = 1L) {
  now = format(Sys.time(), "%Y%m%d_%H%M%S")
  purrr::invoke(expand.grid, alt, stringsAsFactors = FALSE) %>%
    filter_valid_LP() %>%
    dplyr::slice(rep(seq_len(nrow(.)), times = times, each = each)) %>%
    dplyr::mutate(o = paste(now, dplyr::row_number(), sep = "_")) %>%
    dplyr::mutate(!!!const)
}

vectorize_args = function(.tbl) {
  purrr::pmap(.tbl, function(...) {
    .params = c(...)
    .names = names(.params)
    .template = ifelse(nchar(.names) > 1L, "--%s=%s", "-%s%s")
    sprintf(.template, .names, .params)
  })
}

# prior is a named list of generating functions
generate_args = function(prior, const = NULL, n = 1L) {
  generate_valid(prior, n = n) %>%
    vectorize_args() %>%
    purrr::map(~ c(const, .x))
}

generate_valid = function(prior, n = 1L) {
  output = NULL
  k = 0L
  while (k < n) {
    generated = purrr::map_dfc(prior, function(f, x) {
      f(x)
    }, x = n - k)
    generated = suppressMessages(filter_valid_LP(generated))
    output = dplyr::bind_rows(output, generated)
    k = nrow(output)
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
    tidyr::crossing(L = c("const"), P = c("roulette", "minstraight", "stroll"))
  )
}
