#' Evaluate heterogeneity
#'
#' `math_score` calculates MATH (mutant-allele tumor heterogeneity) score from variant allele frequencies.
#' @inheritParams stats::mad
#' @return numeric
#' @rdname heterogeneity
#' @export
math_score = function(x, constant = 1.4826, na.rm = FALSE) {
  med = stats::median(x, na.rm = na.rm)
  mad = stats::mad(x, center = med, constant = constant, na.rm = na.rm)
  mad / med
}
