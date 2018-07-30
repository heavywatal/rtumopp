#' Functions for event history
#'
#' `extract_history` extracts birth and death events from raw population data
#' @param population data.frame
#' @rdname demography
#' @export
extract_history = function(population) {
  population %>%
    dplyr::select(.data$id, .data$ancestor, .data$birth, .data$death) %>%
    tidyr::gather("event", "time", c("birth", "death")) %>%
    dplyr::filter(!(.data$event == "death" & .data$time == 0)) %>% # alive
    dplyr::mutate(event = factor(.data$event, levels = c("death", "birth"))) %>%
    dplyr::arrange(.data$time, .data$event) %>%
    dplyr::mutate(size = cumsum(ifelse(.data$event == "birth", 1L, -1L)))
}

#' @description
#' `summarise_demography` simplifies history
#' @param history data.frame from `extract_history`
#' @rdname demography
#' @export
summarise_demography = function(history) {
  history %>%
    dplyr::group_by(.data$time) %>%
    dplyr::summarise(size = utils::tail(.data$size, 1L))
}

list_clade_founders = function(population, num_clades) {
  history = extract_history(population)
  demography = summarise_demography(history)
  indices = which(demography$size == num_clades)
  the_time = demography$time[utils::tail(indices, 1L)]
  history = dplyr::filter(history, .data$time <= the_time)
  id_born = dplyr::filter(history, .data$event == "birth")$id
  id_dead = dplyr::filter(history, .data$event == "death")$id
  setdiff(id_born, id_dead)
}
