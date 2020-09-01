#' Calculate F_ST from matrix
#'
#' @details
#' `fst_between` estimates F_ST.
#' - theta: Weir and Cockerham (1984); Weir and Hill (2002).
#' - N_ST: Lynch and Crease (1990).
#' - <F_ST>: Hudson, Slatkin, and Maddison (1992) Genetics.
#' @param m matrix output from [dist_genealogy()] or [dist_vaf()]
#' @param sample_sizes cell numbers of each region
#' @rdname fst
#' @export
fst_between = function(m, sample_sizes = NULL) {
  if (is.null(sample_sizes)) {
    avg_t_within = mean(diag(m, names = FALSE))
    avg_t_between = mean(lower_tri(m))
  } else {
    w = num_pairs(sample_sizes)
    avg_t_within = weighted_mean(diag(m, names = FALSE), diag(w, names = FALSE))
    avg_t_between = weighted_mean(lower_tri(m), lower_tri(w))
  }
  1 - avg_t_within / avg_t_between
}

#' @details
#' `fst_total` estimates F_ST.
#' - G_ST, gamma_ST: Nei (1973, 1977, 1982).
#' - K_ST: Hudson, Boos, and Kaplan (1992) MBE.
#' @rdname fst
#' @export
fst_total = function(m, sample_sizes = NULL) {
  if (is.null(sample_sizes)) {
    avg_t_within = mean(diag(m, names = FALSE))
    avg_t_total = mean(m)
  } else {
    w = num_pairs(sample_sizes)
    avg_t_within = weighted_mean(diag(m, names = FALSE), diag(w, names = FALSE))
    avg_t_total = weighted_mean(m, w)
  }
  1 - avg_t_within / avg_t_total
}

#' @rdname fst
#' @export
pairwise_fst = function(m, sample_sizes = NULL) {
  ltri = lower.tri(m)
  col_tri = col(m)[ltri]
  row_tri = row(m)[ltri]
  diag_m = diag(m, names = FALSE)
  if (is.null(sample_sizes)) {
    within = 0.5 * (diag_m[col_tri] + diag_m[row_tri])
  } else {
    diag_w = diag(sample_sizes, names = FALSE)
    wm = diag_w * diag_m
    within = (wm[col_tri] + wm[row_tri]) / (diag_w[col_tri] + diag_w[row_tri])
  }
  between = m[ltri]
  m[ltri] = 1 - within / between
  stats::as.dist(m)
}

num_pairs = function(sample_sizes) {
  m = sample_sizes %*% t(sample_sizes)
  diag(m) = diag(m, names = FALSE) * (sample_sizes - 1) / sample_sizes
  m
}

lower_tri = function(x) x[lower.tri(x)]

weighted_mean = function(x, w, na.rm = TRUE) {
  sum(x * w, na.rm = na.rm) / sum(w, na.rm = na.rm)
}
