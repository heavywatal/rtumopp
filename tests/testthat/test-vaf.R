context("test-vaf")

test_that("make_vaf works", {
  result = tumopp("-N1024")
  population = result$population[[1L]]
  extant = filter_extant(population)
  graph = make_igraph(population)
  .colnames = c("x", "y", "z", "id")
  expect_named({regions = sample_uniform_regions(extant, nsam = 4L, ncell = 4L)}, .colnames)
  .colnames = c("site", "sample", "frequency")
  expect_named({vaf = make_vaf(graph, regions$id, mu = 1)}, .colnames)
})
