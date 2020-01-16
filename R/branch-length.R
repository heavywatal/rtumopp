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
  if (length(from)) {
    nzero = sum(from %in% to)
  } else {
    nzero = nrow(m)
  }
  sum(m) / (nrow(m) * ncol(m) - nzero)
}

#' @details
#' `within_between_samples` summarizes branch lengths.
#' @param regions output of `sample_uniform_regions()`
#' @rdname branch-length
#' @export
within_between_samples = function(graph, regions) {
  rows = regions %>%
    tibble::rowid_to_column() %>%
    dplyr::mutate(
      id = lapply(.data$id, igraphlite::as_vids, graph = graph),
      within = purrr::map_dbl(.data$id, ~ mean_branch_length(graph, .x))
    ) %>%
    purrr::transpose()
  purrr::cross2(rows, rows, .filter = ~ .x$rowid >= .y$rowid) %>%
    purrr::map_dfr(~ {
      row_i = .x[[1L]]
      row_j = .x[[2L]]
      tibble::tibble(
        region_i = row_i$rowid,
        region_j = row_j$rowid,
        within_i = row_i$within,
        within_j = row_j$within,
        euclidean = dist_euclidean(row_i, row_j),
        between = mean_branch_length(graph, row_i$id, row_j$id)
      )
    }) %>%
    dplyr::mutate(
      within = 0.5 * (.data$within_i + .data$within_j),
      # HSM = fst_HSM(.data$within, .data$between),
      fst = fst_HBK(.data$within, .data$between)
    )
}

# Fst by Hudson, Slatkin, and Maddison (1992).
fst_HSM = function(within, between) {
  1 - within / between
}

# Kst by Hudson, Boos, and Kaplan (1992).
fst_HBK = function(within, between, n = 2) {
  (between - within) / (between + within / (n - 1))
}

Kst = function(graph, regions) {
  vids = lapply(regions$id, igraphlite::as_vids, graph = graph)
  t_within = purrr::map_dbl(vids, ~ mean_branch_length(graph, .x))
  t_total = mean_branch_length(graph, unique(unlist(vids)))
  1 - mean(t_within) / t_total
}

fst_HBK_roa = function(vaf) {
  avg_H_S = vaf %>%
    dplyr::mutate_all(Hexp) %>%
    dplyr::summarise_all(mean) %>%
    unlist()
  combinations(vaf) %>%
    dplyr::mutate(
      H_S = (avg_H_S[.data$region_i] + avg_H_S[.data$region_j]) / 2,
      H_T = purrr::map2_dbl(.data$region_i, .data$region_j, function(region_i, region_j) {
        p_T = (vaf[[region_i]] + vaf[[region_j]]) / 2
        mean(Hexp(p_T))
      }),
      fst = 1 - .data$H_S / .data$H_T
    )
}

fst_HBK_aor = function(vaf) {
  H_S = dplyr::mutate_all(vaf, Hexp)
  df = combinations(vaf)
  df$fst = purrr::pmap_dbl(df, function(region_i, region_j) {
    p_T = (vaf[[region_i]] + vaf[[region_j]]) / 2
    H_T = Hexp(p_T)
    H_Sij = (H_S[[region_i]] + H_S[[region_j]]) / 2
    mean(1 - H_Sij / H_T, na.rm = TRUE)
    # roa
    # 1 - mean(H_Sij, na.rm = TRUE) / mean(H_T, na.rm = TRUE)
  })
  df
}

# recommended
fst_HSM_roa = function(vaf, n = Inf) {
  df = combinations(vaf)
  df$fst = purrr::pmap_dbl(df, function(region_i, region_j) {
    p_i = vaf[[region_i]]
    p_j = vaf[[region_j]]
    N = bhatia_numerator(p_i, p_j, n)
    D = bhatia_denominator(p_i, p_j)
    mean(N, na.rm = TRUE) / mean(D, na.rm = TRUE)
  })
  df
}

# for testing
fst_HSM_aor = function(vaf, n = Inf) {
  df = combinations(vaf)
  df$fst = purrr::pmap_dbl(df, function(region_i, region_j) {
    p_i = vaf[[region_i]]
    p_j = vaf[[region_j]]
    fst = bhatia_numerator(p_i, p_j, n) / bhatia_denominator(p_i, p_j)
    mean(fst, na.rm = TRUE)
  })
  df
}

Hexp = function(p) 2 * p * (1 - p)

combinations = function(vaf) {
  tidyr::crossing(region_i = colnames(vaf), region_j = colnames(vaf)) %>%
    dplyr::filter(.data$region_i < .data$region_j)
}

bhatia_numerator = function(p1, p2, n1 = Inf, n2 = n1) {
  (p1 - p2) ** 2 - p1 * (1 - p1) / (n1 - 1) - p2 * (1 - p2) / (n2 - 1)
}

bhatia_denominator = function(p1, p2) {
  p1 * (1 - p2) + p2 * (1 - p1)
}

fst_HSM_bhatia = function(p1, p2, n1 = Inf, n2 = n1) {
  bhatia_numerator(p1, p2, n1, n2) / bhatia_denominator(p1, p2)
}
