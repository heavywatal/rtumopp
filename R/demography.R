#' Extract demography from raw population data
#'
#' @param population data.frame
#' @rdname demography
#' @export
extract_demography = function(population) {
  population %>%
    dplyr::select(.data$birth, .data$death) %>%
    tidyr::gather("event", "time", c("birth", "death")) %>%
    dplyr::filter(!(.data$time == 0 & .data$event == "death")) %>% # alive
    dplyr::mutate(event = factor(.data$event, levels = c("death", "birth"))) %>%
    dplyr::arrange(.data$time, .data$event) %>%
    dplyr::mutate(dn = ifelse(.data$event == "birth", 1, -1)) %>%
    dplyr::group_by(.data$time) %>%
    dplyr::summarise(dn = sum(.data$dn)) %>%
    dplyr::mutate(size = cumsum(.data$dn))
}
