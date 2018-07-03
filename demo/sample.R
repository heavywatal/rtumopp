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

.distances = within_between_samples(.subgraph, .regions) %>% print()

.xmax = max(.distances$euclidean)
.ymax = max(max(.distances$fst), 0.6)
.distances %>%
  ggplot(aes(euclidean, fst)) +
  geom_point() +
  stat_smooth(method = lm, formula = y ~ x + 0) +
  coord_cartesian(xlim = c(0, .xmax), ylim = c(0, .ymax)) +
  theme_wtl()
