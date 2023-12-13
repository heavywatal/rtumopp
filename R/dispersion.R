#' Utilities for dispersion

#' `summary_row()` and `mutate_chisq()` calculate various dispersion stats.
#' `stats_dispersion()` is a shortcut to call them both at once.
#' @param x A numeric vector.
#' @rdname dispersion
#' @export
stats_dispersion = function(x) {
  summary_row(x) |> mutate_chisq()
}

#' @rdname dispersion
#' @export
summary_row = function(x) {
  tibble::tibble_row(
    nsam = length(x),
    mean = mean(x),
    var = stats::var(x)
  )
}

#' @param .data A data.frame.
#' @param ... passed to [dplyr::mutate()].
#' @rdname dispersion
#' @export
mutate_chisq = function(.data, ...) {
  dplyr::mutate(.data,
    R = !!as.name("var") / !!as.name("mean"),
    chisq = !!as.name("R") * (!!as.name("nsam") - 1),
    ...
  )
}

#' @details
#' `rpois_dispersion()` summarizes the results of [stats::rpois()].
#' @param n,lambda passed to [stats::rpois()].
#' @param nrep number of replications.
#' @rdname dispersion
#' @export
rpois_dispersion = function(n, lambda, nrep = 1L) {
  seq_len(nrep) |>
    purrr::map(\(.i) summary_row(stats::rpois(n, lambda))) |>
    purrr::list_rbind(names_to = "repl") |>
    mutate_chisq()
}

#' @details
#' `tibble_dchisq()` generates a tibble.
#' @rdname dispersion
#' @export
tibble_dchisq = function(x, n) {
  tibble::tibble(chisq = x, density = stats::dchisq(x, df = n - 1))
}
