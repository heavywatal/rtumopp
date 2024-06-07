test_that("edge_lengths works", {
  result = tumopp("--seed=42 -N32 -k1e6")
  population = result$population[[1L]]
  graph = make_igraph(population)
  exp_len = nrow(population) - 1L
  expect_vector(edge_lengths(graph), integer(), exp_len)
  expect_vector(edge_lengths(graph, mu = 1), integer(), exp_len)
  expect_vector(edge_lengths(graph, accel = 1), numeric(), exp_len)
})
