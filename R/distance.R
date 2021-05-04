#' Calculate genetic distance from genealogy or VAF

#' @details
#' `dist_genealogy` summarizes branch lengths within/between regions.
#' @param graph igraph
#' @param cell_ids id column of sampled regions
#' @seealso [fst()]
#' @rdname distance
#' @export
dist_genealogy = function(graph, cell_ids) {
  vids = lapply(cell_ids, igraphlite::as_vids, graph = graph)
  n_regions = length(vids)
  seq_n = seq_len(n_regions)
  mtrx = matrix(0, n_regions, n_regions, dimnames = list(seq_n, seq_n))
  for (i in seq_n) {
    for (j in seq.int(i, n_regions)) {
      mtrx[i, j] = mean_branch_length(graph, vids[[i]], vids[[j]])
      mtrx[j, i] = mtrx[i, j]
    }
  }
  mtrx
}

#' @details
#' `dist_vaf` summarizes diversity within/between regions.
#' @param vaf data.frame or matrix. e.g., output from [make_vaf()].
#' @inheritParams sample_uniform_regions
#' @rdname distance
#' @export
dist_vaf = function(vaf, ncell = Inf) {
  n_regions = ncol(vaf)
  seq_n = seq_len(n_regions)
  mtrx = matrix(0, n_regions, n_regions, dimnames = list(seq_n, seq_n))
  for (i in seq_n) {
    for (j in seq.int(i, n_regions)) {
      p_i = vaf[, i, drop = TRUE]
      p_j = vaf[, j, drop = TRUE]
      mtrx[i, j] = mtrx[j, i] = mean(p_hetero(p_i, p_j), na.rm = TRUE)
    }
  }
  diag(mtrx) = diag(mtrx) / (1 - 1 / ncell)
  mtrx
}

p_hetero = function(p1, p2 = p1) {
  p1 * (1 - p2) + (1 - p1) * p2
}

#' @details
#' `pairwise_distances` summarizes physical and genetic distances.
#' @param regions output of [sample_uniform_regions()]
#' @rdname distance
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
