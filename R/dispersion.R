#' Utilities for dispersion

#' `stats_dispersion()` calculates various dispersion stats.
#' @param x A numeric vector
#' @rdname dispersion
#' @export
stats_dispersion = function(x) {
  tibble::tibble_row(
    nsam = length(x),
    mean = mean(x),
    devsq = sum((x - .data$mean) ** 2),
    chisq = .data$devsq / .data$mean,
    var = .data$devsq / (.data$nsam - 1),
    R = .data$var / .data$mean
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
    purrr::map(\(.i) stats_dispersion(stats::rpois(n, lambda))) |>
    purrr::list_rbind(names_to = "repl")
}

#' @details
#' `tibble_dchisq()` generates a tibble.
#' @rdname dispersion
#' @export
tibble_dchisq = function(x, n) {
  tibble::tibble(chisq = x, density = stats::dchisq(x, df = n - 1))
}
