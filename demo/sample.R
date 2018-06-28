library(tidyverse)
library(wtl)
library(tumopp)

wtl::refresh("rtumopp")

(.result = tumopp(str_split("-N40000 -D2 -Chex -k24 -Lconst", " ")[[1]]))
(.population = .result$population[[1]])
(.extant = .population %>% filter_extant())
(.graph = make_igraph(.population))

(.regions = sample_uniform_regions(.extant, 8L, 100L))

.extant %>%
  dplyr::left_join(tumopp:::tidy_regions(.regions), by = "id") %>%
  plot_lattice2d(size = 0.3) +
  geom_point(data = function(x) {dplyr::filter(x, !is.na(region))}, aes(x, y), size = 0.3, alpha = 0.4) +
  scale_colour_brewer(palette = "Spectral", guide = FALSE) +
  theme(axis.title = element_blank())

# #######1#########2#########3#########4#########5#########6#########7#########

.tidy = make_vaf(.graph, .regions$id, -1) %>% print()

.subgraph = tumopp::subtree(.graph, purrr::flatten_chr(.regions$id))
.mutated = mutate_clades(.subgraph, mu = 1)
.mutated = mutate_clades(.subgraph, mu = -1)
.mutated = mutate_clades(.subgraph, segsites=1000L)

.vaf = tally_vaf(.regions$id, .mutated %>% purrr::map(as.integer)) %>% print()
.tidy = .vaf %>%
  filter_detectable(0.05) %>%
  sort_vaf() %>%
  # dplyr::filter(rowSums(. > 0) > 1L) %>%
  tidy_vaf() %>%
  print()

ggplot(.tidy, aes(sample, site)) +
  geom_tile(aes(fill = frequency)) +
  scale_fill_distiller(palette = "Spectral", limit = c(0, 1), guide = FALSE) +
  coord_cartesian(expand = FALSE)

# #######1#########2#########3#########4#########5#########6#########7#########

wtl::refresh("rtumopp")

.threshold = 0.01
.detectable = .regions$id %>% purrr::map(~shared_ancestors(.graph, as.character(.x), .threshold) %>% as.integer())
.combn_biopsy = .detectable %>% combn_ids() %>% print()
.tbl = summarize_capture_rate(.combn_biopsy, .population, .threshold) %>% print()
.tbl %>% plot_capture_rate()

.combn_capture_rate = function(population, regions, graph, thr_range = c(0.01, 0.03, 0.05), ...) {
  purrr::map_dfr(thr_range, function(thr) {
    regions$id %>%
      purrr::map(~shared_ancestors(graph, as.character(.x), thr) %>% as.integer()) %>%
      combn_ids() %>%
      summarize_capture_rate(population, thr)
  })
}
.tbl_capture = .combn_capture_rate(.population, .regions, .graph) %>% print()
.tbl_capture %>% plot_capture_rate() + facet_wrap(~threshold)

.plot_allelefreq_biopsy = function(population, sample_ancestors, threshold) {
  .df = population %>%
    dplyr::filter(.data$allelefreq >= threshold) %>%
    dplyr::transmute(
      .data$id,
      .data$allelefreq,
      captured = .data$id %in% sample_ancestors
    )
  ggplot(.df, aes(allelefreq, group=captured))+
    geom_histogram(aes(fill=captured), binwidth=0.02)+
    theme_bw()
}
.plot_allelefreq_biopsy(.population, .combn_biopsy$id[[3]], .threshold)

# #######1#########2#########3#########4

Rprof()
.distances = within_between_samples(.subgraph, .regions) %>% print()
Rprof(NULL)
summaryRprof()

.xmax = max(.distances$euclidean)
.ymax = max(max(.distances$fst), 0.6)
.distances %>%
  ggplot(aes(euclidean, fst)) +
  geom_point() +
  stat_smooth(method = lm, formula = y ~ x + 0) +
  coord_cartesian(xlim = c(0, .xmax), ylim = c(0, .ymax)) +
  theme_wtl()

# #######1#########2#########3#########4#########5#########6#########7#########

.tr_P = c(random='Push-1', roulette='Push-2', mindrag='Push-3', minstraight='Push-4')
.tr_L = c(const='Constant-rate', step='Step-function', linear='Linear-function')

setwd("~/working/tumopp/spatial")
result_dirs = fs::dir_ls(type = "directory")
raw_results = tumopp::read_results(result_dirs)

# saveRDS(raw_results, "raw_results.rds")

results = raw_results %>%
  select_if(~ n_distinct(.x) > 1L) %>%
  dplyr::select(-outdir, -seed, -drivers) %>%
  print()

.population = results$population[[1]]

df_sampled = results %>%
  dplyr::mutate(
    local = factor(.tr_L[local], levels=.tr_L),
    path = factor(.tr_P[path], levels=.tr_P),
    extant = purrr::map(population, filter_extant),
    regions = purrr::map(extant, sample_uniform_regions, nsam=6L, ncell=100L),
    graph = purrr::map(population, make_igraph),
    subgraph = purrr::map2(graph, regions, ~tumopp::subtree(.x, purrr::flatten_chr(.y$id)))
  ) %>% print()

# saveRDS(df_sampled, "df_sampled.rds")

df_distance = df_sampled %>%
  dplyr::mutate(distance = purrr::map2(subgraph, regions, within_between_samples)) %>%
  print()

df_distance$distance[[1L]] %>%
  ggplot(aes(euclidean, fst)) +
  geom_point() +
  geom_smooth(method = lm, formula = y ~ 0 + x, se = FALSE)

.p_fst = df_distance %>%
  dplyr::select(local, path, shape, distance) %>%
  dplyr::group_by(local, path, shape) %>%
  dplyr::mutate(replicate = seq_along(local)) %>%
  tidyr::unnest() %>%
  ggplot(aes(euclidean, fst)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = lm, formula = y ~ 0 + x, se = FALSE, colour = "#000000", alpha = 0.5) +
  facet_grid(local + path ~ shape + replicate) +
  coord_cartesian(ylim = c(0, 1))
.p_fst
ggsave('fst-6.png', .p_fst, width = 12, height = 12)


df_vaf = df_sampled %>%
  dplyr::mutate(tidy_vaf = purrr::map2(subgraph, regions, ~make_vaf(.x, .y$id, mu = -1))) %>%
  print()

.p_vaf = df_vaf %>%
  dplyr::select(local, path, shape, tidy_vaf) %>%
  dplyr::group_by(local, path, shape) %>%
  dplyr::mutate(replicate = seq_along(local)) %>%
  tidyr::unnest() %>%
  ggplot(aes(sample, site)) +
  geom_tile(aes(fill = frequency)) +
  scale_fill_distiller(palette = "Spectral", limit = c(0, 1), guide = FALSE) +
  coord_cartesian(expand = FALSE) +
  facet_grid(local + path ~ shape + replicate, scales = "free_y")
.p_vaf
ggsave('vaf-6.png', .p_vaf, width = 12, height = 12)


df_extant = df_sampled %>%
  dplyr::mutate(extant = purrr::map2(extant, regions, ~{
    .x %>%
      dplyr::left_join(tumopp:::tidy_regions(.y), by = "id") %>%
      dplyr::mutate(clade = as.factor(as.integer(clade)))
  })) %>%
  dplyr::select(local, path, shape, extant) %>%
  dplyr::group_by(local, path, shape) %>%
  dplyr::mutate(replicate = seq_along(local)) %>%
  tidyr::unnest() %>%
  print()

.p_sampled = df_extant %>%
  sample_n(40000) %>%
  plot_lattice2d(size = 0.3) +
  geom_point(data = function(x) {dplyr::filter(x, !is.na(region))}, aes(x, y), size = 0.3, alpha = 0.4) +
  scale_colour_brewer(palette = "Spectral", guide = FALSE) +
  facet_grid(local + path ~ shape + replicate) +
  wtl::erase(axis.title, axis.text, axis.ticks) +
  theme(panel.spacing = unit(0, 'mm'))
.p_sampled
ggsave('samples-6.png', .p_sampled, width = 12, height = 12)

df_capture = df_sampled %>%
  dplyr::mutate(
    capture_tbl = purrr::pmap(., .combn_capture_rate, thr_range = c(0.05))
  ) %>% print()

# saveRDS(df_capture, "df_capture.rds")

df_capture_tidy = df_capture %>%
  dplyr::select(local, path, shape, capture_tbl) %>%
  dplyr::group_by(local, path, shape) %>%
  dplyr::mutate(replicate = seq_along(local)) %>%
  tidyr::unnest() %>%
  print()

df_capture_tidy %>%
  plot_capture_rate() +
  facet_grid(local + path ~ shape + replicate + threshold) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  theme(panel.spacing = unit(0, 'mm'))

.xrange = data.frame(nsam = range(df_capture_tidy$nsam) + c(-0.9, 0.9), capture_rate = 0.5)
.p_capture = df_capture_tidy %>%
  ggplot2::ggplot(ggplot2::aes_(~nsam, ~capture_rate)) +
    ggplot2::geom_ribbon(data=.xrange, ymin = 0.9, ymax = 1.0, fill = "dodgerblue", alpha = 0.5) +
    ggplot2::geom_ribbon(data=.xrange, ymin = 0.8, ymax = 0.9, fill = "orange", alpha = 0.5) +
    ggplot2::stat_summary(fun.y = mean, geom = "bar", fill = "#666666") +
    # ggplot2::stat_summary(fun.data = wtl::mean_sd, geom = "errorbar", width = 0.2, size = 0.3, colour = "#333333") +
    ggplot2::coord_cartesian(ylim = c(0, 1), xlim = .xrange$nsam, expand = FALSE) +
  facet_grid(local + path ~ shape + replicate) +
  scale_y_continuous(breaks = c(0, 0.5, 1,0)) +
  theme(panel.grid = element_blank())
.p_capture
ggsave("capture_rate-6.png", .p_capture, width=11, height=12)
