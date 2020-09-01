#' Calculate F_ST from matrix
#'
#' @details
#' `fst` calculates Hudson's F_ST estimator.
#' - theta: Weir and Cockerham (1984); Weir and Hill (2002).
#' - N_ST: Lynch and Crease (1990).
#' - <F_ST>: Hudson, Slatkin, and Maddison (1992) Genetics.
#' @param m matrix output from [dist_genealogy()] or [dist_vaf()]
#' @param sample_sizes cell numbers of each region
#' @seealso [dist_genealogy()]
#' @rdname fst
#' @export
fst = function(m, sample_sizes = NULL) {
  if (is.null(sample_sizes)) {
    h_within = mean(diag(m, names = FALSE))
    h_between = mean(lower_tri(m))
  } else {
    w = num_pairs(sample_sizes)
    h_within = weighted_mean(diag(m, names = FALSE), diag(w, names = FALSE))
    h_between = weighted_mean(lower_tri(m), lower_tri(w))
  }
  1 - h_within / h_between
}

#' @details
#' `gst` calculates Nei's F_ST estimator.
#' - G_ST, gamma_ST: Nei (1973, 1977, 1982).
#' - K_ST: Hudson, Boos, and Kaplan (1992) MBE.
#' @rdname fst
#' @export
gst = function(m, sample_sizes = NULL) {
  if (is.null(sample_sizes)) {
    h_within = mean(diag(m, names = FALSE))
    h_total = mean(m)
  } else {
    w = num_pairs(sample_sizes)
    h_within = weighted_mean(diag(m, names = FALSE), diag(w, names = FALSE))
    h_total = weighted_mean(m, w)
  }
  1 - h_within / h_total
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
