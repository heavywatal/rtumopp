test_that("fst works", {
  set.seed(42L)
  result = tumopp("-N256 -k1e6")
  population = result$population[[1L]]
  extant = filter_extant(population)
  graph = make_igraph(population)
  nsam = 4L
  ncell = 4L
  regions = sample_uniform_regions(extant, nsam = nsam, ncell = ncell)
  subgraph = subtree(graph, unlist(regions$id))
  vaf = make_vaf(subgraph, regions$id, mu = -1)
  expect_silent({
    d_vaf = dist_vaf(vaf, ncell)
  })
  expect_identical(dim(d_vaf), c(4L, 4L))
  expect_silent({
    d_tree = dist_genealogy(subgraph, regions$id)
  })
  expect_identical(dim(d_tree), c(4L, 4L))
  expect_equal(fst(d_vaf), fst(d_tree), tolerance = 1e-9)
  expect_equal(gst(d_vaf), gst(d_tree), tolerance = 1e-9)
})
