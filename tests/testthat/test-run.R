test_that("tumopp runs", {
  expect_s3_class(tumopp("-N8"), "data.frame")

  .const = c(D = 3, C = "hex", k = 1, N = 8)
  .alt = list(
    L = c("const", "step"),
    P = c("random", "mindrag", "roulette")
  )
  expect_message({
    argslist = make_args(alt = .alt, const = .const, each = 2L)
  }, "invalid and excluded")
  expect_equal(nrow(argslist), 2L * prod(lengths(.alt)) - 2L)
  expect_s3_class(argslist, "data.frame")

  expect_s3_class(tumopp(argslist[1:2], mc.cores = 1L), "data.frame")
})
