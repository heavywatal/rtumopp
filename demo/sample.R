library(tidyverse)
library(tumopp)

# load_all()

(result = tumopp("-N40000 -D2 -Chex -k24 -Lconst"))
(population = result$population[[1L]])
(graph = result$graph[[1L]])
(extant = population |> filter_extant())

(regions = sample_uniform_regions(extant, 8L, 100L))
subgraph = tumopp::subtree(graph, unlist(regions$id))

extant |>
  dplyr::left_join(tumopp:::tidy_regions(regions), by = "id") |>
  plot_lattice2d(size = 0.3) +
  geom_point(data = function(x) {
    dplyr::filter(x, !is.na(region))
  }, size = 0.3, alpha = 0.4) +
  theme(axis.title = element_blank())

# #######1#########2#########3#########4#########5#########6#########7#########

mutated = mutate_clades(subgraph, mu = 1)
mutated = mutate_clades(subgraph, mu = -1)
mutated = mutate_clades(subgraph, segsites = 1000L)

.vaf = make_vaf(subgraph, regions$id, mu = -1) |> print()
.tidy = .vaf |>
  filter_detectable(0.05) |>
  sort_vaf() |>
  longer_vaf() |>
  print()

.tidy = make_longer_vaf(graph, regions$id, -1) |> print()

ggplot(.tidy) +
  aes(sample, site) +
  geom_tile(aes(fill = frequency)) +
  scale_fill_distiller(palette = "Spectral", limit = c(0, 1), guide = FALSE) +
  coord_cartesian(expand = FALSE)


testdf = tibble::tibble(mu = rep(c(1, 4, 16, 64), each = 200)) |>
  dplyr::mutate(fst = wtl::mcmap_dbl(mu, ~ tumopp::make_vaf(subgraph, regions$id, mu = .x) |>
    tumopp::dist_vaf(ncell) |>
    tumopp::fst())) |>
  print()

# load_all()
xi = tumopp::make_vaf(subgraph, regions$id, mu = -1) |>
  tumopp::dist_vaf(ncell) |>
  tumopp::fst()
ggplot(testdf) +
  aes(fst) +
  geom_histogram(bins = 30) +
  geom_vline(xintercept = xi) +
  facet_wrap(~mu)


# #######1#########2#########3#########4#########5#########6#########7#########

distances = pairwise_distances(subgraph, regions) |> print()

.xmax = max(distances$euclidean)
.ymax = max(max(distances$fst), 0.6)
distances |>
  ggplot(aes(euclidean, fst)) +
  geom_point() +
  stat_smooth(method = lm, formula = y ~ x + 0) +
  coord_cartesian(xlim = c(0, .xmax), ylim = c(0, .ymax))

m = dist_genealogy(subgraph, regions[["id"]])
fst(m)
gst(m)

w = num_pairs(lengths(regions[["id"]]))
