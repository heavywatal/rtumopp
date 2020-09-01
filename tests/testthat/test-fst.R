test_that("make_vaf works", {
  result = tumopp("--seed=42 -N256 -k1e6")
  population = result$population[[1L]]
  extant = filter_extant(population)
  graph = make_igraph(population)
  nsam = 4L
  ncell = 4L
  regions = sample_uniform_regions(extant, nsam = nsam, ncell = ncell)
  sampled = purrr::flatten_int(regions$id)
  subgraph = subtree(graph, sampled)
  mutated = mutate_clades(subgraph, regions$id, mu = -1)
  vaf = tally_vaf(regions$id, mutated$carriers)
  expect_silent({
    d_vaf = dist_vaf(vaf, ncell)
  })
  expect_equal(dim(d_vaf), c(4L, 4L))
  expect_silent({
    d_tree = dist_genealogy(subgraph, regions$id)
  })
  expect_equal(dim(d_tree), c(4L, 4L))
  expect_equal(fst_between(d_vaf), fst_between(d_tree))
})
