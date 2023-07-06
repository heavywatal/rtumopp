# load_all()

(result = tumopp::tumopp("-N40000 -D2 -Chex -k24 -Lconst"))
(population = result$population[[1L]])
(graph = result$graph[[1L]])
(extant = population |> tumopp::filter_extant())

(regions = tumopp::sample_uniform_regions(extant, 8L, 100L))
subgraph = tumopp::subtree(graph, unlist(regions$id))

extant |>
  dplyr::left_join(tumopp:::tidy_regions(regions), by = "id") |>
  tumopp::plot_lattice2d(size = 0.3) +
  ggplot2::geom_point(data = function(x) {
    dplyr::filter(x, !is.na(region))
  }, size = 0.3, alpha = 0.4) +
  ggplot2::theme(axis.title = ggplot2::element_blank())

# #######1#########2#########3#########4#########5#########6#########7#########

mutated = tumopp:::mutate_clades(subgraph, mu = 1)
mutated = tumopp:::mutate_clades(subgraph, mu = -1)
mutated = tumopp:::mutate_clades(subgraph, segsites = 1000L)

.vaf = tumopp::make_vaf(subgraph, regions$id, mu = -1) |> print()
.tidy = .vaf |>
  tumopp::filter_detectable(0.05) |>
  tumopp::sort_vaf() |>
  tumopp::longer_vaf() |>
  print()

.tidy = tumopp::make_longer_vaf(graph, regions$id, -1) |> print()

ggplot2::ggplot(.tidy) +
  ggplot2::aes(sample, site) +
  ggplot2::geom_tile(ggplot2::aes(fill = frequency)) +
  ggplot2::scale_fill_distiller(palette = "Spectral", limit = c(0, 1), guide = FALSE) +
  ggplot2::coord_cartesian(expand = FALSE)

ncell = 100L
testdf = tibble::tibble(mu = rep(c(1, 4, 16, 64), each = 200)) |>
  dplyr::mutate(fst = wtl::mcmap_dbl(mu, \(x) tumopp::make_vaf(subgraph, regions$id, mu = x) |>
    tumopp::dist_vaf(ncell) |>
    tumopp::fst())) |>
  print()

# load_all()
xi = tumopp::make_vaf(subgraph, regions$id, mu = -1) |>
  tumopp::dist_vaf(ncell) |>
  tumopp::fst()
ggplot2::ggplot(testdf) +
  ggplot2::aes(fst) +
  ggplot2::geom_histogram(bins = 30) +
  ggplot2::geom_vline(xintercept = xi) +
  ggplot2::facet_wrap(ggplot2::vars(mu))


# #######1#########2#########3#########4#########5#########6#########7#########

distances = tumopp::pairwise_distances(subgraph, regions) |> print()

.xmax = max(distances$euclidean)
.ymax = max(max(distances$fst), 0.6)
distances |>
  ggplot2::ggplot() +
  ggplot2::aes(euclidean, fst) +
  ggplot2::geom_point() +
  ggplot2::stat_smooth(method = lm, formula = y ~ x + 0) +
  ggplot2::coord_cartesian(xlim = c(0, .xmax), ylim = c(0, .ymax))

m = tumopp::dist_genealogy(subgraph, regions[["id"]])
tumopp::fst(m)
tumopp::gst(m)

w = tumopp:::num_pairs(lengths(regions[["id"]]))
