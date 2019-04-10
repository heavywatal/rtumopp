context("test-io")

test_that("write-read roundtrip works", {
  result = tumopp(c("-N8", "-o", tempfile()))
  expect_message(write_results(result), "outdir")
  expect_silent({
    roundtrip = read_results(result$outdir, graph = TRUE)
  })
  roundtrip$directory = NULL
  expect_equal(dplyr::select_if(result, is.atomic), dplyr::select_if(roundtrip, is.atomic))
  expect_equal(result$population[[1]], roundtrip$population[[1]])
})
