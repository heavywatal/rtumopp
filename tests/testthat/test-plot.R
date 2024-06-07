test_that("augment_genealogy() works", {
  result = tumopp("--seed=42 -N16 -k1e6")
  graph = result$graph[[1L]]

  genealogy = augment_genealogy(graph) |>
    expect_s3_class("data.frame")
  expect_s3_class(gggenealogy(genealogy), "gg")

  lens = edge_lengths(graph, mu = 3)
  augment_genealogy(graph, lengths = lens) |>
    gggenealogy() |>
    expect_s3_class("gg")
})
