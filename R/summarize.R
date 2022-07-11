#' Summarizing functions
#'
#' `genetic_stats` calculates summary statistics of genealogy.
#' @param extant a data.frame of extant cells
#' @rdname summarize
#' @export
genetic_stats = function(extant) {
  dplyr::summarize(
    extant,
    mean_age = mean(.data$age), sd_age = stats::sd(.data$age),
    max_age = max(.data$age), min_age = min(.data$age),
    evenness = evenness(.data$clade)
  )
}

#' @details
#' `morphological_stats` calculates summary statistics of morphology.
#' @param coord switch to normalize phi
#' @rdname summarize
#' @export
morphological_stats = function(extant, coord = "") {
  max_phi = switch(coord,
    hex = 12,
    moore = 27,
    neumann = 6,
    1
  )
  extant_surface = dplyr::filter(extant, .data$surface)
  extant_surface |>
    dplyr::mutate(
      r = dist_euclidean(extant_surface),
      phi = .data$phi / max_phi
    ) |>
    dplyr::summarise(
      phi_mean = mean(.data$phi), phi_sd = stats::sd(.data$phi),
      r_mean = mean(.data$r), r_sd = stats::sd(.data$r)
    ) |>
    dplyr::mutate(surface = sum(extant$surface) / nrow(extant))
}

# H / H_max
evenness = function(species) {
  freqs = table(species)
  shannon_index(freqs) / log(length(freqs))
}

shannon_index = function(freqs, base = exp(1)) {
  freqs = freqs / sum(freqs)
  -sum(freqs * log(freqs, base))
}
