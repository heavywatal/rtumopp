#' Run tumopp
#'
#' `tumopp()` returns full results with config columns in a data.frame
#' @param args command line arguments as a string vector or list of strings
#' @param npair number of samples to measure genetic and physical distance
#' @param nsam number of samples for ms-like output
#' @export
tumopp = function(args = character(0L), npair = 0L, nsam = 0L) {
  if (is.list(args)) {
    purrr::map_dfr(args, tumopp, .id = "args")
  } else {
    if (length(args) == 1L) {
      args = stringr::str_split(args, "\\s+") %>% purrr::flatten_chr()
    }
    message(paste(args, collapse = " "))
    nrep = as.integer(nsam > 0L)
    result = cpp_tumopp(c(nsam, nrep, args), npair = npair)
    if (length(result) == 0L) return(invisible(NULL))
    .out = wtl::read_boost_ini(result["config"]) %>%
      dplyr::mutate(population = list(readr::read_tsv(result["specimens"]))) %>%
      dplyr::mutate(graph = purrr::map(population, make_igraph)) %>%
      dplyr::mutate(population = purrr::pmap(., modify_population)) %>%
      dplyr::mutate(drivers = list(readr::read_tsv(result["drivers"])))
    if (npair > 0L) {
      .dist = readr::read_tsv(result["distances"])
      .out = .out %>% dplyr::mutate(distances = list(.dist))
    }
    if (nsam > 0L) {
      .out$ms = wtl::parse_ms(strsplit(result["ms"], "\n")[[1L]])
    }
    .out
  }
}

#' `mslike()` returns only binary genotypes in ms-like format
#' @rdname tumopp
#' @export
mslike = function(nsam = 20L, args = character(0L)) {
  strsplit(cpp_tumopp_ms(nsam, args), "\n")[[1L]]
}

#' `make_args()` returns argument combinations in a list
#' @param alt named list of altered arguments
#' @param const unnamed vector of constant arguments
#' @param nreps number of repeats
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
