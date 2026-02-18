test_that("make_vaf works", {
  set.seed(42L)
  result = tumopp("-N256 -k1e6")
  population = result$population[[1L]]
  extant = filter_extant(population)
  graph = make_igraph(population)
  nsam = 4L
  ncell = 4L
  regions = sample_uniform_regions(extant, nsam = nsam, ncell = ncell) |>
    expect_named(c("x", "y", "z", "id"))
  vaf_tidy = make_longer_vaf(graph, regions$id, mu = 1) |>
    expect_named(c("site", "sample", "frequency"))
  subgraph = subtree(graph, unlist(regions$id)) |>
    expect_silent()
  df_distances = pairwise_distances(subgraph, regions) |>
    expect_silent() |>
    expect_named(c("i", "j", "euclidean", "fst")) |>
    expect_shape(dim = c(as.integer(choose(nsam, 2L)), 4L)) |>
    expect_identical(pairwise_distances(graph, regions))
  expect_type(
    internal_nodes(subgraph, unlist(regions$id), sensitivity = 0.10),
    "integer"
  )
})
