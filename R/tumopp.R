#' Run tumopp
#'
#' `tumopp()` returns full results with config columns in a data.frame.
#' @param args command line arguments as a string vector or list of strings.
#' @param ... not used.
#' @export
tumopp = function(args, ...) UseMethod("tumopp")

#' @rdname tumopp
#' @export
tumopp.default = function(args = character(0L), ...) {
  if (length(args) == 1L) {
    args = stringr::str_split(args, "\\s+") %>% purrr::flatten_chr()
  }
  message(paste(args, collapse = " "))
  result = cpp_tumopp(args)
  if (length(result) == 0L) return(invisible(NULL))
  .out = from_json(result["config"])
  .pop = read_tumopp(result["population"])
  .snapshots = read_tumopp(result["snapshots"])
  if ((.out$coord == "hex") && getOption("tumopp.autohex", TRUE)) {
    .pop = trans_coord_hex(.pop)
    .snapshots = trans_coord_hex(.snapshots)
  }
  .out = .out %>% dplyr::mutate(
    population = list(.pop),
    graph = list(make_igraph(.pop))
  )
  if (nrow(.snapshots) > 0L) {
    .out = .out %>% dplyr::mutate(snapshots = list(.snapshots))
  }
  .drivers = readr::read_tsv(result["drivers"])
  if (nrow(.drivers) > 0L) {
    .out = .out %>% dplyr::mutate(drivers = list(.drivers))
  }
  .out
}

#' @param mc.cores The number of cores to use for concurrent execution.
#' @rdname tumopp
#' @export
tumopp.list = function(args, ..., mc.cores = getOption("mc.cores", 1L)) {
  parallel::mclapply(args, tumopp, ..., mc.cores = mc.cores) %>%
    dplyr::bind_rows(.id = "args")
}

#' @description
#' `make_args()` returns argument combinations in a list.
#' @param alt named list of altered arguments.
#' @param const unnamed vector of constant arguments.
#' @param nreps number of repeats.
#' @rdname tumopp
#' @export
make_args = function(alt, const = NULL, nreps = 1L) {
  purrr::invoke(expand.grid, alt, stringsAsFactors = FALSE) %>%
    filter_valid_LP() %>%
    vectorize_args() %>%
    append_o() %>%
    purrr::map(~c(const, .x))
}

vectorize_args = function(.tbl) {
  purrr::pmap(.tbl, function(...) {
    .params = c(...)
    .names = names(.params)
    .template = ifelse(nchar(.names) > 1L, "--%s=%s", "-%s%s")
    sprintf(.template, .names, .params)
  })
}

append_o = function(args, fmt = "%Y%m%d_%H%M%S_") {
  prefix = format(Sys.time(), fmt)
  o = paste0("-o", prefix, seq_along(args))
  purrr::map2(args, o, c)
}

# prior is a named list of generating functions
generate_args = function(prior, const = NULL, n = 1L) {
  generate_valid(prior, n = n) %>%
    vectorize_args() %>%
    purrr::map(~c(const, .x))
}

generate_valid = function(prior, n = 1L) {
  output = NULL
  k = 0L
  while (k < n) {
    generated = purrr::map_dfc(prior, function(f, x) {f(x)}, x = n - k)
    generated = suppressMessages(filter_valid_LP(generated))
    output = dplyr::bind_rows(output, generated)
    k = NROW(output)
  }
  output
}

filter_valid_LP = function(x) {
  if (all(c("L", "P") %in% names(x))) {
    nrow_orig = nrow(x)
    x = dplyr::inner_join(x, valid_LP_combinations(), by = c("L", "P"))
    num_excluded = nrow_orig - nrow(x)
    if (num_excluded > 0L) {
      message(num_excluded, " L-P combinations were invalid and excluded")
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
