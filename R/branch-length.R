#' Calculate genetic distance with igraph
#'
#' @details
#' `mean_branch_length` calculates mean branch length within/between sub-graphs.
#' @param graph igraph
#' @param from,to igraph vertices (not cell ID)
#' @rdname branch-length
#' @export
mean_branch_length = function(graph, from = graph$sink, to = from) {
  # TODO: Avoid creating huge matrix
  m = igraphlite::shortest_paths(graph, from, to, mode = 3L, algorithm = "unweighted")
  # Exclude self-comparison only if "within"
  n = length(m) - if (identical(from, to)) {length(from)} else {0}
  sum(m) / n
}

#' @details
#' `dist_genealogy` summarizes branch lengths within/between regions.
#' @param cell_ids id column of sampled regions
#' @rdname branch-length
#' @export
dist_genealogy = function(graph, cell_ids) {
  vids = lapply(cell_ids, igraphlite::as_vids, graph = graph)
  n_regions = length(vids)
  .dn = list(seq_len(n_regions), seq_len(n_regions))
  mtrx = matrix(0, n_regions, n_regions, dimnames = .dn)
  for (i in seq_len(n_regions)) {
    for (j in seq.int(i, n_regions)) {
      mtrx[i, j] = mean_branch_length(graph, vids[[i]], vids[[j]])
      mtrx[j, i] = mtrx[i, j]
    }
  }
  mtrx
}

#' @details
#' `fst_between` estimates F_ST.
#' theta: Weir and Cockerham (1984); Weir and Hill (2002)
#' N_ST: Lynch and Crease (1990)
#' <F_ST>: Hudson, Slatkin, and Maddison (1992) Genetics
#' @param m output from `dist_genealogy()`
#' @param sample_sizes cell numbers of each region
#' @rdname branch-length
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
#' G_ST, gamma_ST: Nei (1973, 1977, 1982)
#' K_ST: Hudson, Boos, and Kaplan (1992) MBE
#' @rdname branch-length
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

#' @details
#' `pairwise_distances` summarizes physical and genetic distances.
#' @param regions output of `sample_uniform_regions()`
#' @rdname branch-length
#' @export
pairwise_distances = function(graph, regions) {
  m = dist_genealogy(graph, regions[["id"]])
  .dim = dim(m)
  ltri = lower.tri(m)
  tibble::tibble(
    i = .col(.dim)[ltri],
    j = .row(.dim)[ltri],
    euclidean = c(stats::dist(regions[c("x", "y", "z")])),
    fst = c(pairwise_fst(m))
  )
}

pairwise_fst = function(x, weight = NULL) {
  ltri = lower.tri(x)
  col_tri = col(x)[ltri]
  row_tri = row(x)[ltri]
  diag_x = diag(x, names = FALSE)
  if (is.null(weight)) {
    within = 0.5 * (diag_x[col_tri] + diag_x[row_tri])
  } else {
    diag_w = diag(weight, names = FALSE)
    wx = diag_w * diag_x
    within = (wx[col_tri] + wx[row_tri]) / (diag_w[col_tri] + diag_w[row_tri])
  }
  between = x[ltri]
  x[ltri] = 1 - within / between
  stats::as.dist(x)
}

num_pairs = function(sample_sizes) {
  m = sample_sizes %*% t(sample_sizes)
  diag(m) = diag(m, names = FALSE) * (sample_sizes - 1) / sample_sizes
  m
}

lower_tri = function(x) x[lower.tri(x)]

offdiag = function(x) x[lower.tri(x) | upper.tri(x)]

weighted_mean = function(x, w, na.rm = TRUE) {
  sum(x * w, na.rm = na.rm) / sum(w, na.rm = na.rm)
}
