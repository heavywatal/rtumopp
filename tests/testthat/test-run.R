context("test-run")

test_that("tumopp runs", {
  expect_s3_class(tumopp('-N8'), "data.frame")
})
