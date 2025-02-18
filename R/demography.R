#' Functions for event history
#'
#' `extract_history` extracts birth and death events from raw population data.
#' @param population data.frame
#' @rdname demography
#' @export
extract_history = function(population) {
  population |>
    dplyr::select("id", "ancestor", "birth", "death") |>
    tidyr::pivot_longer(c("birth", "death"), names_to = "event", values_to = "time") |>
    dplyr::filter(!(.data$event == "death" & .data$time == 0)) |> # alive
    dplyr::mutate(event = factor(.data$event, levels = c("death", "birth"))) |>
    dplyr::arrange(.data$time, .data$event) |>
    dplyr::mutate(size = cumsum(ifelse(.data$event == "birth", 1L, -1L)))
}

#' @details
#' `summarize_demography` simplifies history.
#' @param history data.frame from [extract_history()]
#' @rdname demography
#' @export
summarize_demography = function(history) {
  history |>
    dplyr::group_by(.data$time) |>
    dplyr::summarize(size = max(!!as.name("size")))
}

list_clade_founders = function(population, n) {
  stopifnot(n < 16L)
  .history = extract_history(population) |> utils::head(10000L)
  demography = summarize_demography(.history)
  indices = which(demography$size == n)
  the_time = demography$time[utils::tail(indices, 1L)]
  .history = dplyr::filter(.history, .data$time <= the_time)
  id_born = dplyr::filter(.history, .data$event == "birth")$id
  id_dead = dplyr::filter(.history, .data$event == "death")$id
  setdiff(id_born, id_dead)
}
