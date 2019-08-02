test_that("make_vaf works", {
  result = tumopp("--seed=42 -N256 -k1e6")
  population = result$population[[1L]]
  extant = filter_extant(population)
  graph = make_igraph(population)
  nsam = 4L
  ncell = 4L
  .colnames = c("x", "y", "z", "id")
  expect_named({
    regions = sample_uniform_regions(extant, nsam = nsam, ncell = ncell)
  }, .colnames)
  .colnames = c("site", "sample", "frequency")
  expect_named({
    vaf = make_vaf(graph, regions$id, mu = 1)
  }, .colnames)
  sampled = purrr::flatten_int(regions$id)
  expect_silent({
    subgraph = subtree(graph, sampled)
  })
  expect_silent({
    df_distances = within_between_samples(subgraph, regions)
  })
  expect_named(df_distances, c("region_i", "region_j", "within_i", "within_j", "euclidean", "between", "within", "fst"))
  expect_equal(nrow(df_distances), choose(nsam, 2L))
  expect_equal(df_distances, within_between_samples(graph, regions))
  expect_type(internal_nodes(subgraph, sampled, sensitivity = 0.10), "integer")
})
