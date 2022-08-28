test_that("write-read roundtrip works", {
  result = tumopp(c("-N8", "-o", tempfile()))
  expect_message(write_results(result), "outdir")
  expect_silent({
    roundtrip = read_results(result$outdir, graph = TRUE)
  })
  roundtrip$directory = NULL
  expect_identical(dplyr::select(result, where(is.atomic)), dplyr::select(roundtrip, where(is.atomic)))
  expect_identical(result$population[[1]], roundtrip$population[[1]])
})
