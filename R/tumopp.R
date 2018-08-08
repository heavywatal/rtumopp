#' Run tumopp
#'
#' `tumopp()` returns full results with config columns in a data.frame.
#' @param args command line arguments as a string vector or list of strings.
#' @param ... not used.
#' @export
tumopp = function(args, ...) UseMethod("tumopp")

#' @param nsam number of samples for ms-like output.
#' @rdname tumopp
#' @export
tumopp.default = function(args = character(0L), nsam = 0L, ...) {
  if (length(args) == 1L) {
    args = stringr::str_split(args, "\\s+") %>% purrr::flatten_chr()
  }
  message(paste(args, collapse = " "))
  nrep = as.integer(nsam > 0L)
  result = cpp_tumopp(c(nsam, nrep, args))
  if (length(result) == 0L) return(invisible(NULL))
  .out = wtl::read_boost_ini(result["config"])
  .pop = read_tumopp(result["specimens"])
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
  .dist = readr::read_tsv(result["distances"])
  if (nrow(.dist) > 0L) {
    .out = .out %>% dplyr::mutate(distances = list(.dist))
  }
  if (nsam > 0L) {
    .out$ms = wtl::parse_ms(strsplit(result["ms"], "\n")[[1L]])
  }
  .out
}

#' @param mc.cores The number of cores to use for concurrent execution.
#' @rdname tumopp
#' @export
tumopp.list = function(args, nsam = 0L, ..., mc.cores = getOption("mc.cores", 1L)) {
  parallel::mclapply(args, tumopp, nsam = nsam, mc.cores = mc.cores) %>%
    dplyr::bind_rows(.id = "args")
}

#' @description
#' `mslike()` returns only binary genotypes in ms-like format.
#' @rdname tumopp
#' @export
mslike = function(nsam = 20L, args = character(0L)) {
  strsplit(cpp_tumopp_ms(nsam, args), "\n")[[1L]]
}

#' @description
#' `make_args()` returns argument combinations in a list.
#' @param alt named list of altered arguments.
#' @param const unnamed vector of constant arguments.
#' @param nreps number of repeats.
#' @rdname tumopp
#' @export
make_args = function(alt, const = NULL, nreps = 1L) {
  altered = purrr::invoke(expand.grid, alt, stringsAsFactors = FALSE) %>%
    filter_valid_LP()
  prefix = format(Sys.time(), "%Y%m%d_%H%M_")
  paste0(prefix, seq_len(nreps)) %>%
    purrr::map_dfr(~dplyr::mutate(altered, o = .x)) %>%
    purrr::pmap(function(...) {
      .params = c(...)
      .names = names(.params)
      .template = ifelse(nchar(.names) > 1L, "--%s=%s", "-%s%s")
      c(const, sprintf(.template, .names, .params))
    }) %>%
    stats::setNames(purrr::map_chr(., paste, collapse = " "))
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
