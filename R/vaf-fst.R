gamma_np = function(vaf) {
  pi_within = dplyr::mutate_all(vaf, Hexp)
  pi_total = Hexp(rowMeans(vaf))
  gamma_i = 1 - colMeans(pi_within) / mean(pi_total)
  mean(gamma_i)
}

hsm_np = function(vaf) {
  pi_within = dplyr::mutate_all(vaf, Hexp)
  for (i in seq_along(vaf)) {
    vaf_o = vaf[-i]
    pi_o = dplyr::mutate_all(vaf_o, Hexp)
    avg_pi_within = 0.5 * (pi_within[[i]] + pi_o)
    pi_between = vaf[[i]] * (1 - vaf_o) + (1 - vaf[[i]]) * vaf_o
    print(str(avg_pi_within))
    print(str(pi_between))
    mean(1 - avg_pi_within / pi_between)
  }
  gamma_i = 1 - colMeans(pi_within) / mean(pi_total)
  mean(gamma_i)
}

pairwise_fst_HBK = function(vaf) {
  H_S = dplyr::mutate_all(vaf, Hexp)
  df = combinations(vaf)
  out = purrr::pmap_dfr(df, function(region_i, region_j) {
    p_T = (vaf[[region_i]] + vaf[[region_j]]) / 2
    H_T = Hexp(p_T)
    H_Sij = (H_S[[region_i]] + H_S[[region_j]]) / 2
    data.frame(
      fst = 1 - mean(H_Sij, na.rm = TRUE) / mean(H_T, na.rm = TRUE),
      fst_aor = mean(1 - H_Sij / H_T, na.rm = TRUE)
    )
  })
  dplyr::bind_cols(df, out)
}

pairwise_fst_HSM = function(vaf, n = Inf) {
  df = combinations(vaf)
  out = purrr::pmap_dfr(df, function(region_i, region_j) {
    p_i = vaf[[region_i]]
    p_j = vaf[[region_j]]
    N = bhatia_numerator(p_i, p_j, n)
    D = bhatia_denominator(p_i, p_j)
    data.frame(
      fst = mean(N, na.rm = TRUE) / mean(D, na.rm = TRUE),
      fst_aor = mean(N / D, na.rm = TRUE)
    )
  })
  dplyr::bind_cols(df, out)
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
