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
    dplyr::mutate(within = 0.5 * (.data$within_i + .data$within_j))
}

#' @details
#' `pairwise_branch_length` summarizes lengths within/between regions.
#' @param regions output of `sample_uniform_regions()`
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

if (FALSE) {
# #######1#########2#########3#########4#########5#########6#########7#########

load_all()

result = tumopp::tumopp("-N20000 -D3 -Chex -k10 -Lconst")
population = result$population[[1L]]
extant = population %>% tumopp::filter_extant()
graph = tumopp::make_igraph(population)
regions = tumopp::sample_uniform_regions(extant, nsam = 5L, ncell = 200L)
subgraph = tumopp::subtree(graph, purrr::flatten_int(regions$id))

eu = dist(regions[,-4])
m = dist_genealogy(subgraph, regions[["id"]])
w = num_pairs(lengths(regions[["id"]]))

bench::mark(
  within_between_samples(subgraph, regions),
  dist(regions[,-4]),
  dist_genealogy(subgraph, regions[["id"]]),
  check = FALSE
)

dist_genealogy(subgraph, regions[["id"]])
mean(pop_spec_fst(subgraph, regions[["id"]], "between"))
mean(pop_spec_fst(subgraph, regions[["id"]], "approx"))
mean(pop_spec_fst(subgraph, regions[["id"]], "total"))
mean(pop_spec_fst(subgraph, regions[["id"]], "unbiased"))

bench::mark(
  fst_between(subgraph, regions[["id"]]),
  fst_tfor(subgraph, regions[["id"]]),
  fst_total(subgraph, regions[["id"]]),
  fst_unbiased(subgraph, regions[["id"]]),
  check = FALSE
)

plot(c(dist(regions[,-4])), m[lower.tri(m)])

# #######1#########2#########3#########4#########5#########6#########7#########
}

num_pairs = function(sample_sizes) {
  m = sample_sizes %*% t(sample_sizes)
  diag(m) = diag(m) * (sample_sizes - 1) / sample_sizes
  m
}

offdiag = function(x) x[lower.tri(x) | upper.tri(x)]

pop_spec_fst = function(graph, cell_ids, reference = c("between", "total", "approx", "unbiased")) {
  reference = match.arg(reference)
  sample_sizes = lengths(cell_ids)
  m = dist_genealogy(graph, cell_ids)
  w = num_pairs(sample_sizes)
  avg_t_ref = if (reference %in% c("total", "unbiased")) {
    weighted_mean(m, w)
  } else {
    weighted_mean(offdiag(m), offdiag(w))
  }
  if (reference == "approx") {
    (avg_t_ref - diag(m)) / (avg_t_ref + diag(m) / (length(cell_ids) - 1))
  } else if (reference == "unbiased") {
    (avg_t_ref - diag(m)) / (avg_t_ref - diag(m) / length(cell_ids))
  } else {
    1 - diag(m) / avg_t_ref
  }
}

fst_unbiased = function(graph, cell_ids) {
  sample_sizes = lengths(cell_ids)
  m = dist_genealogy(graph, cell_ids)
  w = num_pairs(sample_sizes)
  avg_t_ref = weighted_mean(m, w)
  (avg_t_ref - diag(m)) / (avg_t_ref - diag(m) / length(cell_ids))
}

fst_between = function(graph, cell_ids, weight = FALSE) {
  sample_sizes = lengths(cell_ids)
  m = dist_genealogy(graph, cell_ids)
  w = num_pairs(sample_sizes)
  avg_t_between = weighted_mean(offdiag(m), offdiag(w))
  avg_t_within = if (weight) {
    weighted_mean(diag(m), sample_sizes)
  } else {
    mean(diag(m))
  }
  gamma_Slatkin(avg_t_within, avg_t_between, length(cell_ids))
}

fst_tfor = function(graph, cell_ids, weight = FALSE) {
  sample_sizes = lengths(cell_ids)
  m = dist_genealogy(graph, cell_ids)
  w = num_pairs(sample_sizes)
  avg_t_total = weighted_mean(m, w)
  avg_t_within = if (weight) {
    weighted_mean(diag(m), sample_sizes)
  } else {
    mean(diag(m))
  }
  1 - avg_t_within / avg_t_total
}

# Eq. 8 in Slatkin 1991 Genome Res.
fst_total = function(graph, cell_ids, weight = FALSE) {
  vids = lapply(cell_ids, igraphlite::as_vids, graph = graph)
  avg_t_total = mean_branch_length(graph, unique(unlist(vids)))
  t_within = purrr::map_dbl(vids, mean_branch_length, graph = graph)
  avg_t_within = if (weight) {
    # K_ST by Hudson, Boos, and Kaplan (1992) MBE
    weighted_mean(t_within, lengths(cell_ids))
  } else {
    mean(t_within)
  }
  1 - avg_t_within / avg_t_total
}

weighted_mean = function(x, w, na.rm = TRUE) {
  sum(x * w, na.rm = na.rm) / sum(w, na.rm = na.rm)
}

pairwise_fst = function(x, weight = NULL) {
  ltri = lower.tri(x)
  col_tri = col(x)[ltri]
  row_tri = row(x)[ltri]
  diag_x = diag(x, names = FALSE)
  if (missing(weight)) {
    within = 0.5 * (diag_x[col_tri] + diag_x[row_tri])
  } else {
    w = diag(weight, names = FALSE)
    wx = w * diag_x
    within = (wx[col_tri] + wx[row_tri]) / (w[col_tri] + w[row_tri])
  }
  between = x[ltri]
  x[ltri] = gamma_Slatkin(within, between, 2)
  as.dist(x)
}

# F_ST by Hudson, Slatkin, and Maddison (1992).
# gamma_ST by Nei (1982).
gamma_Slatkin = function(within, between, n = Inf) {
  (between - within) / (between + within / (n - 1))
}

hsm_Slatkin = function(within, between) {
  1 - within / between
}

is_within = function(sample_sizes) {
  num_inds = sum(sample_sizes)
  v = rep.int(seq_along(sample_sizes), sample_sizes)
  row_pop = matrix(v, num_inds, num_inds)
  row_pop == t(row_pop)
}
